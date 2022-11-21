program main

   implicit none

   INTEGER,PARAMETER 	:: rp = KIND(0.d0)
   INTEGER :: it,pp,t_steps
   INTEGER :: nRE
   REAL(rp),PARAMETER :: C_PI = 4.0_rp*ATAN(1.0_rp) !< Definition of @f$\pi@f$
   REAL(rp),PARAMETER :: C_E = 1.602176E-19_rp !< Absolute value of electron charge in Coulombs (C).
   REAL(rp),PARAMETER :: C_ME = 9.109382E-31_rp !< Electron mass in kg
   REAL(rp),PARAMETER :: C_C = 299792458.0_rp !< Light speed in m/s
   REAL(rp),ALLOCATABLE,DIMENSION(:) :: X_X,X_Y,X_Z
   REAL(rp),ALLOCATABLE,DIMENSION(:) :: V_X,V_Y,V_Z
   REAL(rp),ALLOCATABLE,DIMENSION(:) :: B_X,B_Y,B_Z
   REAL(rp),ALLOCATABLE,DIMENSION(:) :: E_X,E_Y,E_Z
   REAL(rp),ALLOCATABLE,DIMENSION(:) :: rnd1,gam
   REAL(rp)	:: dt,simulation_time
   REAL(rp)  :: Eo,gam0,v0,eta0,chi0
   REAL(rp)  :: v_norm,B_norm,t_norm,x_norm
   CHARACTER(100) :: path_to_inputs,path_to_outputs
   REAL(rp)     :: U_X,U_Y,U_Z
   REAL(rp)    :: U_hs_X,U_hs_Y,U_hs_Z
   REAL(rp)    :: tau_X,tau_Y,tau_Z
   REAL(rp)     :: t_X,t_Y,t_Z
   REAL(rp)    :: up_X,up_Y,up_Z
   REAL(rp)     :: cross_X,cross_Y,cross_Z
   REAL(rp)     :: sigma,us,gp,g0,s
   INTEGER,PARAMETER 	:: default_unit_open = 101
   INTEGER,PARAMETER 	:: output_write = 202,data_write = 102
   INTEGER  :: argn,read_stat
   INTEGER  :: c1,c2,cr
   REAL  :: rate
   NAMELIST /input_parameters/ nRE,simulation_time

   !! initialize system_clock
   CALL system_clock(count_rate=cr)
   rate = REAL(cr)

   CALL system_clock(c1)

   !! find input/output file
   argn = command_argument_count()

   call get_command_argument(1,path_to_inputs)
   call get_command_argument(2,path_to_outputs)

   !! open log and data output files
   OPEN(UNIT=data_write, &
   FILE=TRIM(path_to_outputs)//"data.korc", &
   STATUS='UNKNOWN',FORM='FORMATTED',POSITION='REWIND')

   OPEN(UNIT=output_write, &
   FILE=TRIM(path_to_outputs)//"output.korc", &
   STATUS='UNKNOWN',FORM='FORMATTED',POSITION='REWIND')

   !! set defaults for inputs, open input file, read from input file

   nRE=1
   simulation_time=0._rp

   OPEN(UNIT=default_unit_open,FILE=TRIM(path_to_inputs), &
   STATUS='OLD',FORM='formatted',POSITION='REWIND')
   read(default_unit_open,nml=input_parameters,IOSTAT=read_stat)
   close(default_unit_open)

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

   !! Initialize fields
   B_X=0._rp
   B_Y=0._rp
   B_Z=1._rp

   E_X=0._rp
   E_Y=0._rp
   E_Z=0._rp

   write(output_write,*) '* * * * * * * * * Fields * * * * * * * * *'

   write(output_write,*) 'B',B_X(1),B_Y(1),B_Z(1)
   write(output_write,*) 'E',E_X(1),E_Y(1),E_Z(1)

   !! Initialize location
   X_X=0._rp
   X_Y=0._rp
   X_Z=0._rp

   !! Set kinetic energy and pitch, random gyrophase, and then velocity
   Eo=10E6_rp
   eta0=90._rp*C_PI/180._rp
   gam0=1._rp+(Eo*C_E/(C_ME*C_C**2))
   v0=C_C*sqrt(1._rp-1/gam0**2)

   !call RANDOM_NUMBER(rnd1)
   !chi0=2*C_PI*rnd1(1)
   chi0=0

   V_X=v0*sin(eta0)*cos(chi0)
   V_Y=v0*sin(eta0)*sin(chi0)
   V_Z=v0*cos(eta0)
   gam=gam0

   v_norm=C_C
   B_norm=B_Z(1)
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

   write(output_write,*) '* * * * * * * * * Initial Conditions * * * * * * * * *'

   write(output_write,*) 'Number of electrons: ',nRE
   write(output_write,*) 'gam0,eta0,chi0: ',gam0,eta0,chi0
   write(output_write,*) 'X0: ',X_X(1)*x_norm,X_Y(1)*x_norm,X_Z(1)*x_norm
   write(output_write,*) 'V0: ',V_X(1)*v_norm,V_Y(1)*v_norm,V_Z(1)*v_norm
   write(data_write,*) 'X0: ',X_X(1)*x_norm,X_Y(1)*x_norm,X_Z(1)*x_norm
   write(data_write,*) 'V0: ',V_X(1)*v_norm,V_Y(1)*v_norm,V_Z(1)*v_norm

   !! Set timestep to resolve relativistic gyrofrequency, simulation time and
   !! number of time steps
   dt = 0.01_rp*(2.0_rp*C_PI/(C_E*B_Z(1)/( gam0*C_ME )))/t_norm

   simulation_time=simulation_time/t_norm
   t_steps=ceiling(simulation_time/dt)
   it=simulation_time/float(t_steps)

   write(output_write,*) '* * * * * * * * * Timings * * * * * * * * *'

   write(output_write,*) 'simulation time:',simulation_time*t_norm
   write(output_write,*) 'dt:',dt*t_norm
   write(output_write,*) 't_steps:',t_steps

   !! Particle push

   !! Allocating work arrays

   call system_clock(c2)

   write(output_write,*) 'Setup time:',(c2-c1)/rate

   write(output_write,*) '* * * * * * * * * Begin Orbits * * * * * * * * *'

   !$acc parallel loop
   do pp=1,nRE

      !! Initial half step
      X_X(pp) = X_X(pp)+dt/2*V_X(pp)
      X_Y(pp) = X_Y(pp)+dt/2*V_Y(pp)
      X_Z(pp) = X_Z(pp)+dt/2*V_Z(pp)

      !! Main iteration loop
      do it=1,t_steps

         g0=gam(pp)

         U_X = g0*V_X(pp)
         U_Y = g0*V_Y(pp)
         U_Z = g0*V_Z(pp)

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

      end do
   end do
   !$acc end parallel loop

   write(output_write,*) '* * * * * * * * * Initial Conditions * * * * * * * * *'
   write(data_write,*) 'V: ',V_X(1)*v_norm,V_Y(1)*v_norm,V_Z(1)*v_norm
   write(data_write,*) 'X: ',X_X(1)*x_norm,X_Y(1)*x_norm,X_Z(1)*x_norm

   call system_clock(c1)

   write(output_write,*) 'Pusher time:',(c1-c2)/rate

   ! * * * FINALIZING SIMULATION * * *
   write(output_write,'("KORC ran suppessfully!")')
   close(output_write)
   close(data_write)

end program main
