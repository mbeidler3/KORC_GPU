module params_gpu

IMPLICIT NONE

INTEGER,PARAMETER :: rp = KIND(0.d0)
INTEGER,PARAMETER :: is = KIND(INT(1,1))
INTEGER :: it,jt,pp,t_steps
INTEGER :: nRE
REAL(rp),PARAMETER :: C_PI = 4.0_rp*ATAN(1.0_rp) !< Definition of @f$\pi@f$
REAL(rp),PARAMETER :: C_E = 1.602176E-19_rp !< Absolute value of electron charge in Coulombs (C).
REAL(rp),PARAMETER :: C_ME = 9.109382E-31_rp !< Electron mass in kg
REAL(rp),PARAMETER :: C_C = 299792458.0_rp !< Light speed in m/s
REAL(rp),ALLOCATABLE,DIMENSION(:) :: rnd1,gam
REAL(rp)	:: dt,simulation_time
REAL(rp)  :: Eo,gam0,v0,eta0,chi0
CHARACTER(100) :: path_to_inputs,path_to_outputs,field_type
REAL(rp)     :: X_X_loop,X_Y_loop,X_Z_loop
REAL(rp)     :: V_X_loop,V_Y_loop,V_Z_loop
REAL(rp)     :: B_X_loop,B_Y_loop,B_Z_loop
REAL(rp)     :: E_X_loop,E_Y_loop,E_Z_loop
REAL(rp)     :: U_X,U_Y,U_Z
REAL(rp)    :: U_hs_X,U_hs_Y,U_hs_Z
REAL(rp)    :: tau_X,tau_Y,tau_Z
REAL(rp)     :: t_X,t_Y,t_Z
REAL(rp)    :: up_X,up_Y,up_Z
REAL(rp)     :: cross_X,cross_Y,cross_Z
REAL(rp)     :: sigma,us,gp,gam_loop,s
INTEGER,PARAMETER 	:: default_unit_open = 101
INTEGER,PARAMETER 	:: output_write = 202,data_write = 102
INTEGER  :: argn,read_stat
INTEGER  :: c1,c2,cr
REAL(rp)  :: rate
REAL(rp)  :: v_norm,B_norm,t_norm,x_norm

end module params_gpu