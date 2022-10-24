program main

  implicit none

  INTEGER :: it,pp,cc,t_steps
  INTEGER,PARAMETER :: nRE=1
  REAL, PARAMETER :: C_PI = 4.0*ATAN(1.0) !< Definition of @f$\pi@f$
  REAL, PARAMETER :: C_E = 1.602176E-19 !< Absolute value of electron charge in Coulombs (C).
  REAL, PARAMETER :: C_ME = 9.109382E-31 !< Electron mass in kg
  REAL, PARAMETER :: C_C = 299792458.0 !< Light speed in m/s
  REAL,DIMENSION(nRE) :: X_X,X_Y,X_Z
  REAL,DIMENSION(nRE) :: V_X,V_Y,V_Z
  REAL,DIMENSION(nRE) :: Y_R,Y_PHI,Y_Z
  REAL,DIMENSION(nRE) :: B_X,B_Y,B_Z
  REAL,DIMENSION(nRE) :: E_X,E_Y,E_Z
  REAL,DIMENSION(nRE) 			:: rnd1,gam
  REAL			:: dt,a,simulation_time
  REAL  :: Eo,gam0,v0,eta0
  CHARACTER(100) :: path_to_outputs
  INTEGER,PARAMETER :: pchunk=1
  REAL,DIMENSION(pchunk)     :: U_L_X,U_L_Y,U_L_Z
  REAL,DIMENSION(pchunk)     :: U_X,U_Y,U_Z
  REAL,DIMENSION(pchunk)     :: U_RC_X,U_RC_Y,U_RC_Z
  REAL,DIMENSION(pchunk)     :: U_os_X,U_os_Y,U_os_Z
  REAL,DIMENSION(pchunk)     :: U_hs_X,U_hs_Y,U_hs_Z
  REAL,DIMENSION(pchunk)     :: tau_X,tau_Y,tau_Z
  REAL,DIMENSION(pchunk)     :: t_X,t_Y,t_Z
  REAL,DIMENSION(pchunk)     :: up_X,up_Y,up_Z
  REAL,DIMENSION(pchunk)     :: cross_X,cross_Y,cross_Z
  REAL,DIMENSION(pchunk)     :: sigma,us,gp,g0,s,Bmag
  INTEGER,PARAMETER 	:: output_write = 202

  !! open output file
  call get_command_argument(1,path_to_outputs)

  OPEN(UNIT=output_write, &
       FILE=TRIM(path_to_outputs)//"output.korc", &
       STATUS='UNKNOWN',FORM='FORMATTED',POSITION='REWIND')

  !! Initialize fields
  B_X=0.
  B_Y=0.
  B_Z=1.

  E_X=0.
  E_Y=0.
  E_Z=0.

  !! Initialize location
  X_X=0.
  X_Y=0.
  X_Z=0.

  !! Set kinetic energy and pitch, random gyrophase, and then velocity
  Eo=10E6
  eta0=10*C_PI/180
  gam0=1+(Eo*C_E)
  v0=C_C*sqrt(1-1/gam0**2)

  call RANDOM_NUMBER(rnd1)

  V_X=sin(eta0)*cos(2*C_PI*rnd1)
  V_Y=sin(eta0)*sin(2*C_PI*rnd1)
  V_Z=v0*cos(eta0)/C_C
  gam=gam0

  !! Set timestep to resolve relativistic gyrofrequency, simulation time and
  !! number of time steps
  dt = 0.01*(2.0*C_PI/(C_E*B_Z(1)/( gam0*C_ME )))

  simulation_time=1E-6
  t_steps=ceiling(simulation_time/dt)
  dt=simulation_time/float(t_steps)

  a=dt*C_E/C_ME

  !! Particle push
  !$OMP PARALLEL DO &
  !$OMP& FIRSTPRIVATE(dt,a,t_steps) &
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

           U_hs_X(cc) = U_X(cc) + 0.5*a*(E_X(pp-1+cc) +cross_X(cc))
           U_hs_Y(cc) = U_Y(cc) + 0.5*a*(E_Y(pp-1+cc) +cross_Y(cc))
           U_hs_Z(cc) = U_Z(cc) + 0.5*a*(E_Z(pp-1+cc) +cross_Z(cc))

           !write(6,*) 'half step',0.5*a*cross_X(cc),0.5*a*cross_Y(cc),0.5*a*cross_Z(cc)

           tau_X(cc) = 0.5*a*B_X(pp-1+cc)
           tau_Y(cc) = 0.5*a*B_Y(pp-1+cc)
           tau_Z(cc) = 0.5*a*B_Z(pp-1+cc)

           up_X(cc) = U_hs_X(cc) + 0.5*a*E_X(pp-1+cc)
           up_Y(cc) = U_hs_Y(cc) + 0.5*a*E_Y(pp-1+cc)
           up_Z(cc) = U_hs_Z(cc) + 0.5*a*E_Z(pp-1+cc)

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

     end do
  end do
  !$OMP END PARALLEL DO

  ! * * * FINALIZING SIMULATION * * *
  write(output_write,'("KORC ran successfully!")')
  close(output_write)

end program main
