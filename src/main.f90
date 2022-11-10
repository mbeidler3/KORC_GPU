program main

  implicit none

  INTEGER,PARAMETER 	:: rp = KIND(0.d0)
  INTEGER :: it,pp,cc,t_steps
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
  INTEGER :: pchunk
  REAL(rp),ALLOCATABLE,DIMENSION(:)     :: U_X,U_Y,U_Z
  REAL(rp),ALLOCATABLE,DIMENSION(:)     :: U_hs_X,U_hs_Y,U_hs_Z
  REAL(rp),ALLOCATABLE,DIMENSION(:)     :: tau_X,tau_Y,tau_Z
  REAL(rp),ALLOCATABLE,DIMENSION(:)     :: t_X,t_Y,t_Z
  REAL(rp),ALLOCATABLE,DIMENSION(:)     :: up_X,up_Y,up_Z
  REAL(rp),ALLOCATABLE,DIMENSION(:)     :: cross_X,cross_Y,cross_Z
  REAL(rp),ALLOCATABLE,DIMENSION(:)     :: sigma,us,gp,g0,s,Bmag
  INTEGER,PARAMETER 	:: default_unit_open = 101
  INTEGER,PARAMETER 	:: output_write = 202
  INTEGER  :: argn,read_stat

  NAMELIST /input_parameters/ nRE,pchunk,simulation_time

  !! find input/output file
  argn = command_argument_count()

  call get_command_argument(1,path_to_inputs)
  call get_command_argument(2,path_to_outputs)

  !! open output file
  OPEN(UNIT=output_write, &
       FILE=TRIM(path_to_outputs)//"output.korc", &
       STATUS='UNKNOWN',FORM='FORMATTED',POSITION='REWIND')

  !! set defaults for inputs, open input file, read from input file

  nRE=1
  pchunk=1
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

  write(output_write,*) 'B',B_X,B_Y,B_Z
  write(output_write,*) 'E',E_X,E_Y,E_Z

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

  write(output_write,*) 'gam0,eta0,chi0',gam0,eta0,chi0
  write(output_write,*) 'X0',X_X*x_norm,X_Y*x_norm,X_Z*x_norm
  write(output_write,*) 'V0',V_X*v_norm,V_Y*v_norm,V_Z*v_norm

  !! Set timestep to resolve relativistic gyrofrequency, simulation time and
  !! number of time steps
  dt = 0.01_rp*(2.0_rp*C_PI/(C_E*B_Z(1)/( gam0*C_ME )))/t_norm

  simulation_time=simulation_time/t_norm
  t_steps=ceiling(simulation_time/dt)
  dt=simulation_time/float(t_steps)

  write(output_write,*) '* * * * * * * * * Timings * * * * * * * * *'

  write(output_write,*) 'simulation time:',simulation_time*t_norm
  write(output_write,*) 'dt:',dt*t_norm
  write(output_write,*) 't_steps:',t_steps

  write(output_write,*) '* * * * * * * * * Begin Orbits * * * * * * * * *'

  !! Particle push

  !! Allocating work arrays

  ALLOCATE(U_X(pchunk))
  ALLOCATE(U_Y(pchunk))
  ALLOCATE(U_Z(pchunk))
  ALLOCATE(U_hs_X(pchunk))
  ALLOCATE(U_hs_Y(pchunk))
  ALLOCATE(U_hs_Z(pchunk))
  ALLOCATE(tau_X(pchunk))
  ALLOCATE(tau_Y(pchunk))
  ALLOCATE(tau_Z(pchunk))
  ALLOCATE(t_X(pchunk))
  ALLOCATE(t_Y(pchunk))
  ALLOCATE(t_Z(pchunk))
  ALLOCATE(up_X(pchunk))
  ALLOCATE(up_Y(pchunk))
  ALLOCATE(up_Z(pchunk))
  ALLOCATE(cross_X(pchunk))
  ALLOCATE(cross_Y(pchunk))
  ALLOCATE(cross_Z(pchunk))
  ALLOCATE(sigma(pchunk))
  ALLOCATE(us(pchunk))
  ALLOCATE(gp(pchunk))
  ALLOCATE(g0(pchunk))
  ALLOCATE(s(pchunk))
  ALLOCATE(Bmag(pchunk))

  !$OMP PARALLEL DO &
  !$OMP& FIRSTPRIVATE(dt,t_steps) &
  !$OMP& PRIVATE(pp,cc,it,g0,U_X,U_Y,U_Z,cross_X,cross_Y,cross_Z, &
  !$OMP& U_hs_X,U_hs_Y,U_hs_Z,tau_X,tau_Y,tau_Z,up_X,up_Y,up_Z, &
  !$OMP& gp,sigma,us,t_X,t_Y,t_Z,s) &
  !$OMP& SHARED(X_X,X_Y,X_Z,V_X,V_Y,V_Z,B_X,B_Y,B_Z,E_X,E_Y,E_Z,gam)
  do pp=1,nRE,pchunk

     !! Initial half step
     !$OMP SIMD
     do cc=1,pchunk
        X_X(pp-1+cc) = X_X(pp-1+cc)+dt/2*V_X(pp-1+cc)
        X_Y(pp-1+cc) = X_Y(pp-1+cc)+dt/2*V_Y(pp-1+cc)
        X_Z(pp-1+cc) = X_Z(pp-1+cc)+dt/2*V_Z(pp-1+cc)
     end do
     !$OMP END SIMD

     write(output_write,*) 'X1/2',X_X*x_norm,X_Y*x_norm,X_Z*x_norm

     !! Main iteration loop
     do it=1,t_steps

        !$OMP SIMD
        do cc=1,pchunk

           g0(cc)=gam(pp-1+cc)

           U_X(cc) = g0(cc)*V_X(pp-1+cc)
           U_Y(cc) = g0(cc)*V_Y(pp-1+cc)
           U_Z(cc) = g0(cc)*V_Z(pp-1+cc)

           ! Magnitude of magnetic field
           Bmag(cc) = SQRT(B_X(pp-1+cc)*B_X(pp-1+cc)+B_Y(pp-1+cc)*B_Y(pp-1+cc)+ &
               B_Z(pp-1+cc)*B_Z(pp-1+cc))

           ! LEAP-FROG SCHEME FOR LORENTZ FORCE !

           cross_X(cc)=V_Y(pp-1+cc)*B_Z(pp-1+cc)-V_Z(pp-1+cc)*B_Y(pp-1+cc)
           cross_Y(cc)=V_Z(pp-1+cc)*B_X(pp-1+cc)-V_X(pp-1+cc)*B_Z(pp-1+cc)
           cross_Z(cc)=V_X(pp-1+cc)*B_Y(pp-1+cc)-V_Y(pp-1+cc)*B_X(pp-1+cc)

           !write(6,*) 'vcrossB',cross_X,cross_Y,cross_Z

           U_hs_X(cc) = U_X(cc) + 0.5*dt*(E_X(pp-1+cc) +cross_X(cc))
           U_hs_Y(cc) = U_Y(cc) + 0.5*dt*(E_Y(pp-1+cc) +cross_Y(cc))
           U_hs_Z(cc) = U_Z(cc) + 0.5*dt*(E_Z(pp-1+cc) +cross_Z(cc))

           !write(6,*) 'half step',0.5*a*cross_X(cc),0.5*a*cross_Y(cc),0.5*a*cross_Z(cc)

           tau_X(cc) = 0.5*dt*B_X(pp-1+cc)
           tau_Y(cc) = 0.5*dt*B_Y(pp-1+cc)
           tau_Z(cc) = 0.5*dt*B_Z(pp-1+cc)

           up_X(cc) = U_hs_X(cc) + 0.5*dt*E_X(pp-1+cc)
           up_Y(cc) = U_hs_Y(cc) + 0.5*dt*E_Y(pp-1+cc)
           up_Z(cc) = U_hs_Z(cc) + 0.5*dt*E_Z(pp-1+cc)

           gp(cc) = SQRT( 1.0 + up_X(cc)*up_X(cc)+up_Y(cc)*up_Y(cc)+ &
                up_Z(cc)*up_Z(cc) )

           sigma(cc) = gp(cc)*gp(cc) - (tau_X(cc)*tau_X(cc)+ &
                tau_Y(cc)*tau_Y(cc)+tau_Z(cc)*tau_Z(cc))

           us(cc) = up_X(cc)*tau_X(cc)+up_Y(cc)*tau_Y(cc)+ &
                up_Z(cc)*tau_Z(cc)
           ! variable 'u^*' in Vay, J.-L. PoP (2008)

           gam(pp-1+cc) = SQRT( 0.5*(sigma(cc) + SQRT(sigma(cc)*sigma(cc) + &
                4.0*(tau_X(cc)*tau_X(cc)+tau_Y(cc)*tau_Y(cc)+ &
                tau_Z(cc)*tau_Z(cc) + us(cc)*us(cc)))) )

           t_X(cc) = tau_X(cc)/gam(pp-1+cc)
           t_Y(cc) = tau_Y(cc)/gam(pp-1+cc)
           t_Z(cc) = tau_Z(cc)/gam(pp-1+cc)

           s(cc) = 1.0/(1.0 + t_X(cc)*t_X(cc)+t_Y(cc)*t_Y(cc)+ &
                t_Z(cc)*t_Z(cc))
           ! variable 's' in Vay, J.-L. PoP (2008)

           cross_X(cc)=up_Y(cc)*t_Z(cc)-up_Z(cc)*t_Y(cc)
           cross_Y(cc)=up_Z(cc)*t_X(cc)-up_X(cc)*t_Z(cc)
           cross_Z(cc)=up_X(cc)*t_Y(cc)-up_Y(cc)*t_X(cc)

           U_X(cc) = s(cc)*(up_X(cc) + (up_X(cc)*t_X(cc)+ &
                up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_X(cc) + cross_X(cc))
           U_Y(cc) = s(cc)*(up_Y(cc) + (up_X(cc)*t_X(cc)+ &
                up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_Y(cc) + cross_Y(cc))
           U_Z(cc) = s(cc)*(up_Z(cc) + (up_X(cc)*t_X(cc)+ &
                up_Y(cc)*t_Y(cc)+up_Z(cc)*t_Z(cc))*t_Z(cc) + cross_Z(cc))
           ! LEAP-FROG SCHEME FOR LORENTZ FORCE !

           V_X(pp-1+cc) = U_X(cc)/gam(pp-1+cc)
           V_Y(pp-1+cc) = U_Y(cc)/gam(pp-1+cc)
           V_Z(pp-1+cc) = U_Z(cc)/gam(pp-1+cc)

           X_X(pp-1+cc) = X_X(pp-1+cc) + dt*V_X(pp-1+cc)
           X_Y(pp-1+cc) = X_Y(pp-1+cc) + dt*V_Y(pp-1+cc)
           X_Z(pp-1+cc) = X_Z(pp-1+cc) + dt*V_Z(pp-1+cc)
        end do
        !$OMP END SIMD

        write(output_write,*) 'step',it
        write(output_write,*) 'V',V_X*v_norm,V_Y*v_norm,V_Z*v_norm
        write(output_write,*) 'X',X_X*x_norm,X_Y*x_norm,X_Z*x_norm

     end do
  end do
  !$OMP END PARALLEL DO

  ! * * * FINALIZING SIMULATION * * *
  write(output_write,'("KORC ran successfully!")')
  close(output_write)

end program main
