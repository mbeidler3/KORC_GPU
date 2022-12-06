program main

   USE local, ONLY: i4, r8
   USE eval_field_mod
   USE mpi

   implicit none

   !!KORC variables
   INTEGER,PARAMETER 	:: rp = KIND(0.d0)
   INTEGER :: it,pp,t_steps
   INTEGER :: nRE
   REAL(rp),PARAMETER :: C_PI = 4.0_rp*ATAN(1.0_rp) !< Definition of @f$\pi@f$
   REAL(rp),PARAMETER :: C_E = 1.602176E-19_rp !< Absolute value of electron charge in Coulombs (C).
   REAL(rp),PARAMETER :: C_ME = 9.109382E-31_rp !< Electron mass in kg
   REAL(rp),PARAMETER :: C_C = 299792458.0_rp !< Light speed in m/s
   REAL(rp),ALLOCATABLE,DIMENSION(:) :: X_X,X_Y,X_Z,xguess,yguess
   REAL(rp),DIMENSION(3) :: X_0
   REAL(rp),ALLOCATABLE,DIMENSION(:) :: V_X,V_Y,V_Z
   REAL(rp),ALLOCATABLE,DIMENSION(:) :: B_X,B_Y,B_Z
   REAL(rp),ALLOCATABLE,DIMENSION(:) :: E_X,E_Y,E_Z
   REAL(rp),ALLOCATABLE,DIMENSION(:) :: rnd1,gam
   REAL(rp)	:: dt,simulation_time
   REAL(rp) :: Eo,gam0,v0,eta0,chi0
   REAL(rp) :: BMAG,b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z,v1,v2,v3,BMAG1
   REAL(rp) :: v_norm,B_norm,t_norm,x_norm
   CHARACTER(100) :: path_to_inputs,path_to_outputs,nimdumpin,field_model,tmp_str
   REAL(rp)    :: U_X,U_Y,U_Z
   REAL(rp)    :: U_hs_X,U_hs_Y,U_hs_Z
   REAL(rp)    :: tau_X,tau_Y,tau_Z
   REAL(rp)    :: t_X,t_Y,t_Z
   REAL(rp)    :: up_X,up_Y,up_Z
   REAL(rp)    :: cross_X,cross_Y,cross_Z
   REAL(rp)    :: sigma,us,gp,g0,s
   INTEGER,PARAMETER 	:: default_unit_open = 101
   INTEGER,PARAMETER 	:: output_write = 202,data_write = 102
   INTEGER  :: argn,read_stat
   INTEGER  :: c1,c2,cr
   REAL  :: rate
   INTEGER  :: mpierr,mpi_rank,nmpi
   LOGICAL :: mpiinit=.FALSE.,mpi_process_initialized=.FALSE.,all_mpis_initialized=.FALSE.
   !!NIMROD variables
   CHARACTER(64) :: dumpname
   REAL(r8), DIMENSION(1:3) :: RZP
   INTEGER(i4) :: ifail
   REAL(r8), DIMENSION(2) :: xyguess
   INTEGER(i4) , PARAMETER :: nqty=3
   REAL(r8), DIMENSION(1:nqty) :: fval
   !!KORC namelist
   NAMELIST /input_parameters/ nRE,simulation_time,field_model,nimdumpin,X_0,Eo,eta0,chi0

   !! setting up mpi

   call MPI_INIT(mpierr)
   if (mpierr .NE. MPI_SUCCESS) then
      write(6,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
      write(6,'(/," ERROR: Initializing MPI. Aborting... ")')
      write(6,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
      call MPI_ABORT(MPI_COMM_WORLD, -10, mpierr)
   end if

   call MPI_INITIALIZED(mpi_process_initialized,mpierr)

   write(6,*) 'mpi_process_initialized',mpi_process_initialized

   call MPI_REDUCE(mpi_process_initialized,all_mpis_initialized,1, &
   MPI_LOGICAL,MPI_LAND,0,MPI_COMM_WORLD,mpierr)

   call MPI_BCAST(all_mpis_initialized,1, &
   MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)

   call MPI_COMM_SIZE(MPI_COMM_WORLD, nmpi, mpierr)
   if (mpierr .NE. MPI_SUCCESS) then
      write(6,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
      write(6,'(/," ERROR: Obtaining size of communicator. Aborting... ")')
      write(6,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
      call MPI_ABORT(MPI_COMM_WORLD, -10, mpierr)
   end if

   call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, mpierr)
   if (mpierr .NE. MPI_SUCCESS) then
      write(6,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
      write(6,'(/," ERROR: Obtaining MPI rank. Aborting... ")')
      write(6,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
      call MPI_ABORT(MPI_COMM_WORLD, -10, mpierr)
   end if

   write(6,*) 'all_mpis_initialized',all_mpis_initialized

   write(6,*) 'mpi_num',nmpi
   write(6,*) 'mpi_rank',mpi_rank

   if (all_mpis_initialized) then
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

      !! find input/output file
      argn = command_argument_count()

      call get_command_argument(1,path_to_inputs)
      call get_command_argument(2,path_to_outputs)
      
      !! open log and data output files
      write(tmp_str,'(I3.3)') mpi_rank
      OPEN(UNIT=data_write, &
      FILE=TRIM(path_to_outputs)//"data_"//TRIM(tmp_str)//".korc", &
      STATUS='UNKNOWN',FORM='FORMATTED',POSITION='REWIND')
   
      if (mpi_rank.EQ.0) then
         OPEN(UNIT=output_write, &
         FILE=TRIM(path_to_outputs)//"output.korc", &
         STATUS='UNKNOWN',FORM='FORMATTED',POSITION='REWIND')
      endif

      if (mpi_rank.EQ.0) then
         write(output_write,'(/,"* * * * * COMMUNICATIONS  * * * * *")')
         write(output_write,'(/,"  MPI communications initialized!  ")')
         write(output_write,'(/,"  Number of MPI processes: ",I5)') nmpi
         write(output_write,'(/,"* * * * * * * * * * * * * * * * * *")')
      end if
   else
      if (mpi_rank.EQ.0) then
         write(6,'(/,"* * * * * * * COMMUNICATIONS * * * * * * *")')
         write(6,'(/," ERROR: MPI not initialized. Aborting... ")')
         write(6,'(/,"* * * * * * * * * ** * * * * * * * * * * *")')
         call MPI_ABORT(MPI_COMM_WORLD, -10, mpierr)
         STOP
      end if
   endif

   !! initialize system_clock
   CALL system_clock(count_rate=cr)
   rate = REAL(cr)

   CALL system_clock(c1)

   !! set defaults for inputs, open input file, read from input file

   nRE=1
   simulation_time=0._rp
   field_model='ANALYTIC' ! NIMROD
   nimdumpin = 'dump.00000.h5'
   X_0=(/1.6,0.0,0.0/)
   Eo=10E6_rp
   eta0=10._rp
   chi0=0._rp

   OPEN(UNIT=default_unit_open,FILE=TRIM(path_to_inputs), &
   STATUS='OLD',FORM='formatted',POSITION='REWIND')
   read(default_unit_open,nml=input_parameters,IOSTAT=read_stat)
   close(default_unit_open)

   if (field_model=='NIMROD') then
      CALL nimfield_init(nimdumpin,'b')
   endif

   !! Allocate solution fields

   ALLOCATE(X_X(nRE))
   ALLOCATE(X_Y(nRE))
   ALLOCATE(X_Z(nRE))
   ALLOCATE(V_X(nRE))
   ALLOCATE(V_Y(nRE))
   ALLOCATE(V_Z(nRE))
   ALLOCATE(B_X(nRE))
   ALLOCATE(B_Y(nRE))
   ALLOCATE(B_Z(nRE))
   ALLOCATE(E_X(nRE))
   ALLOCATE(E_Y(nRE))
   ALLOCATE(E_Z(nRE))
   ALLOCATE(rnd1(nRE))
   ALLOCATE(gam(nRE))
   if (field_model=='NIMROD') then
      ALLOCATE(xguess(nRE))
      ALLOCATE(yguess(nRE))
   endif

   !! Initialize locations
   X_X=X_0(1)
   X_Y=X_0(2)
   X_Z=X_0(3)

   write(output_write,*) '* * * * * * * * * Fields * * * * * * * * *'

   !! Initialize fields
   if (field_model=='ANALYTICAL') then
      B_X=0._rp
      B_Y=0._rp
      B_Z=1._rp

   else if (field_model=='NIMROD') then

      xguess = 1.0
      yguess = 1.0
      ifail = 0

      !$acc parallel loop
      do pp=1,nRE

         RZP(1)=SQRT(X_X(pp)**2+X_Y(pp)**2)
         RZP(2)=X_Z(pp)
         RZP(3)=ATAN2(X_X(pp),X_Y(pp))

         xyguess(1)=xguess(pp)
         xyguess(2)=yguess(pp)

         CALL eval_field(RZP, 'b', fval, nqty, xyguess, ifail)

         xguess=xyguess(1)
         yguess=xyguess(2)
         B_X(pp)=fval(1)*SIN(RZP(3))+fval(3)*COS(RZP(3))
         B_Y(pp)=fval(1)*COS(RZP(3))-fval(3)*SIN(RZP(3))
         B_Z(pp)=fval(2)
         !WRITE(6,*) "(X,Y,Z)",  X_X(pp),X_Y(pp),X_Z(pp)
         !WRITE(6,*) "(R,Z,PHI)",  RZP
         !WRITE(6,*) "(B_R,B_Z,B_PHI)=",  fval
         !WRITE(6,*) "(B_X,B_Y,B_Z)=",  B_X(pp),B_Y(pp),B_Z(pp)
         !WRITE(6,*) "xyguess=",  xyguess
      end do
      !$acc end parallel loop
   endif

   E_X=0._rp
   E_Y=0._rp
   E_Z=0._rp

   if (mpi_rank.EQ.0) then
      write(output_write,*) 'B',B_X(1),B_Y(1),B_Z(1)
      write(output_write,*) 'E',E_X(1),E_Y(1),E_Z(1)
   endif

   !! Set kinetic energy and pitch, random gyrophase, and then velocity
   eta0=eta0*C_PI/180._rp
   gam0=1._rp+(Eo*C_E/(C_ME*C_C**2))
   v0=C_C*sqrt(1._rp-1/gam0**2)

   !call RANDOM_NUMBER(rnd1)
   !chi0=2*C_PI*rnd1(1)
   chi0=chi0*C_PI/180._rp

   !write(6,*) eta0,gam0,v0,chi0

   BMAG=sqrt(B_X(1)**2+B_Y(1)**2+B_Z(1)**2)
   b1x=B_X(1)/BMAG
   b1y=B_Y(1)/BMAG
   b1z=B_Z(1)/BMAG

   b2x=(b1y-b1z)
   b2y=(b1z-b1x)
   b2z=(b1x-b1y)
   BMAG1=sqrt(b2x**2+b2y**2+b2z**2)
   b2x=b2x/BMAG1
   b2y=b2y/BMAG1
   b2z=b2z/BMAG1

   b3x=(b1y*b2z-b1z*b2y)
   b3y=(b1z*b2x-b1x*b2z)
   b3z=(b1x*b2y-b1y*b2x)

   v1=v0*cos(eta0)
   v2=v0*sin(eta0)*cos(chi0)
   v3=v0*sin(eta0)*sin(chi0)
   
   V_X=v1*b1x+v2*b2x+v3*b3x
   V_Y=v1*b1y+v2*b2y+v3*b3y
   V_Z=v1*b1z+v2*b2z+v3*b3z
   
   !write(6,*) 'b1',b1x,b1y,b1z
   !write(6,*) 'b2',b2x,b2y,b2z
   !write(6,*) 'b3',b3x,b3y,b3z

   !write(6,*) 'v1',v1,v2,v3
   !write(6,*) 'VV',V_X,V_Y,V_Z

   gam=gam0

   v_norm=C_C
   B_norm=B_Y(1)
   t_norm=C_ME/(C_E*B_norm)
   x_norm=V_norm*t_norm

   X_X=X_X/x_norm
   X_Y=X_Y/x_norm
   X_Z=X_Z/x_norm

   V_X=V_X/v_norm
   V_Y=V_Y/v_norm
   V_Z=V_Z/v_norm

   B_X=B_X/b_norm
   B_Y=B_Y/b_norm
   B_Z=B_Z/b_norm

   E_X=E_X/(b_norm*v_norm)
   E_Y=E_Y/(b_norm*v_norm)
   E_Z=E_Z/(b_norm*v_norm)

   if (mpi_rank.EQ.0) then
      write(output_write,*) '* * * * * * * * * Initial Conditions * * * * * * * * *'

      write(output_write,*) 'Number of electrons: ',nRE
      write(output_write,*) 'gam0,eta0,chi0: ',gam0,eta0,chi0
      write(output_write,*) 'X0: ',X_X(1)*x_norm,X_Y(1)*x_norm,X_Z(1)*x_norm
      write(output_write,*) 'V0: ',V_X(1)*v_norm,V_Y(1)*v_norm,V_Z(1)*v_norm
      write(data_write,*) 'X0: ',X_X(1)*x_norm,X_Y(1)*x_norm,X_Z(1)*x_norm
      write(data_write,*) 'V0: ',V_X(1)*v_norm,V_Y(1)*v_norm,V_Z(1)*v_norm
   endif

   !! Set timestep to resolve relativistic gyrofrequency, simulation time and
   !! number of time steps
   dt = 0.01_rp*(2.0_rp*C_PI/(C_E*BMAG/( gam0*C_ME )))/t_norm

   simulation_time=simulation_time/t_norm
   t_steps=ceiling(simulation_time/dt)
   it=simulation_time/float(t_steps)

   if (mpi_rank.EQ.0) then
      write(output_write,*) '* * * * * * * * * Timings * * * * * * * * *'

      write(output_write,*) 'simulation time:',simulation_time*t_norm
      write(output_write,*) 'dt:',dt*t_norm
      write(output_write,*) 't_steps:',t_steps
   endif

   !! Particle push

   !! Allocating work arrays

   call system_clock(c2)

   if (mpi_rank.EQ.0) then
      write(output_write,*) 'Setup time:',(c2-c1)/rate

      write(output_write,*) '* * * * * * * * * Begin Orbits * * * * * * * * *'
   endif

   !$acc parallel loop
   do pp=1,nRE

      !! Initial half step
      X_X(pp) = X_X(pp)+dt/2*V_X(pp)
      X_Y(pp) = X_Y(pp)+dt/2*V_Y(pp)
      X_Z(pp) = X_Z(pp)+dt/2*V_Z(pp)

      write(data_write,*) 'X1/2',X_X*x_norm,X_Y*x_norm,X_Z*x_norm

      !! Main iteration loop
      do it=1,t_steps

         g0=gam(pp)

         U_X = g0*V_X(pp)
         U_Y = g0*V_Y(pp)
         U_Z = g0*V_Z(pp)

         !! Interpolate NIMROD magnetic field
         if (field_model=='NIMROD') then
            RZP(1)=SQRT(X_X(pp)**2+X_Y(pp)**2)
            RZP(2)=X_Z(pp)
            RZP(3)=ATAN2(X_X(pp),X_Y(pp))

            xyguess(1)=xguess(pp)
            xyguess(2)=yguess(pp)

            CALL eval_field(RZP, 'b', fval, nqty, xyguess, ifail)
            
            xguess=xyguess(1)
            yguess=xyguess(2)
            B_X(pp)=fval(1)*SIN(RZP(3))+fval(3)*COS(RZP(3))
            B_Y(pp)=fval(1)*COS(RZP(3))-fval(3)*SIN(RZP(3))
            B_Z(pp)=fval(2)
         endif

         ! LEAP-FROG SCHEME FOR LORENTZ FORCE !

         cross_X=V_Y(pp)*B_Z(pp)-V_Z(pp)*B_Y(pp)
         cross_Y=V_Z(pp)*B_X(pp)-V_X(pp)*B_Z(pp)
         cross_Z=V_X(pp)*B_Y(pp)-V_Y(pp)*B_X(pp)

         !write(6,*) 'vcrossB',cross_X,cross_Y,cross_Z

         U_hs_X = U_X + 0.5*dt*(E_X(pp) +cross_X)
         U_hs_Y = U_Y + 0.5*dt*(E_Y(pp) +cross_Y)
         U_hs_Z = U_Z+ 0.5*dt*(E_Z(pp) +cross_Z)

         !write(6,*) 'half step',0.5*a*cross_X(pp),0.5*a*cross_Y(pp),0.5*a*cross_Z(pp)

         tau_X = 0.5*dt*B_X(pp)
         tau_Y = 0.5*dt*B_Y(pp)
         tau_Z = 0.5*dt*B_Z(pp)

         up_X = U_hs_X + 0.5*dt*E_X(pp)
         up_Y = U_hs_Y + 0.5*dt*E_Y(pp)
         up_Z = U_hs_Z + 0.5*dt*E_Z(pp)

         gp = SQRT( 1.0 + up_X*up_X+up_Y*up_Y+up_Z*up_Z )

         sigma = gp*gp - (tau_X*tau_X+tau_Y*tau_Y+tau_Z*tau_Z)

         us = up_X*tau_X+up_Y*tau_Y+up_Z*tau_Z
         ! variable 'u^*' in Vay, J.-L. PoP (2008)

         gam(pp) = SQRT( 0.5*(sigma + SQRT(sigma*sigma + &
         4.0*(tau_X*tau_X+tau_Y*tau_Y+tau_Z*tau_Z + us*us))) )

         t_X = tau_X/gam(pp)
         t_Y = tau_Y/gam(pp)
         t_Z = tau_Z/gam(pp)

         s = 1.0/(1.0 + t_X*t_X+t_Y*t_Y+t_Z*t_Z)
         ! variable 's' in Vay, J.-L. PoP (2008)

         cross_X=up_Y*t_Z-up_Z*t_Y
         cross_Y=up_Z*t_X-up_X*t_Z
         cross_Z=up_X*t_Y-up_Y*t_X

         U_X = s*(up_X + (up_X*t_X+up_Y*t_Y+up_Z*t_Z)*t_X + cross_X)
         U_Y = s*(up_Y + (up_X*t_X+up_Y*t_Y+up_Z*t_Z)*t_Y + cross_Y)
         U_Z = s*(up_Z + (up_X*t_X+up_Y*t_Y+up_Z*t_Z)*t_Z + cross_Z)
         ! LEAP-FROG SCHEME FOR LORENTZ FORCE !

         V_X(pp) = U_X/gam(pp)
         V_Y(pp) = U_Y/gam(pp)
         V_Z(pp) = U_Z/gam(pp)

         X_X(pp) = X_X(pp) + dt*V_X(pp)
         X_Y(pp) = X_Y(pp) + dt*V_Y(pp)
         X_Z(pp) = X_Z(pp) + dt*V_Z(pp)

         write(data_write,*) 'step',it
         write(data_write,*) 'V',V_X*v_norm,V_Y*v_norm,V_Z*v_norm
         write(data_write,*) 'X',X_X*x_norm,X_Y*x_norm,X_Z*x_norm

      end do
   end do
   !$acc end parallel loop

   write(data_write,*) 'V: ',V_X(1)*v_norm,V_Y(1)*v_norm,V_Z(1)*v_norm
   write(data_write,*) 'X: ',X_X(1)*x_norm,X_Y(1)*x_norm,X_Z(1)*x_norm

   call system_clock(c1)

   if (mpi_rank.EQ.0) then
      write(output_write,*) 'Pusher time:',(c1-c2)/rate

      ! * * * FINALIZING SIMULATION * * *
      write(output_write,'("KORC ran suppessfully!")')
   endif
   close(output_write)
   close(data_write)

   call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
   call MPI_FINALIZE(mpierr)

end program main
