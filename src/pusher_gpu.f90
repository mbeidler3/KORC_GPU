module pusher_gpu

use params_gpu
use interp_gpu 
  
IMPLICIT NONE

CONTAINS

subroutine FO_push(nRE,dt,t_steps,field_type,x_norm,v_norm,X_X,X_Y,X_Z,V_X,V_Y,V_Z, &
   gam,B_X,B_Y,B_Z,E_X,E_Y,E_Z)
  REAL(rp),DIMENSION(nRE),INTENT(INOUT) :: X_X,X_Y,X_Z
  REAL(rp),DIMENSION(nRE),INTENT(INOUT) :: V_X,V_Y,V_Z
  REAL(rp),DIMENSION(nRE),INTENT(INOUT) :: gam
  REAL(rp),DIMENSION(nRE),INTENT(IN) :: B_X,B_Y,B_Z
  REAL(rp),DIMENSION(nRE),INTENT(IN) :: E_X,E_Y,E_Z
  INTEGER,INTENT(IN) :: nRE,t_steps
  REAL(rp),INTENT(IN) :: dt,x_norm,v_norm
  CHARACTER(100),INTENT(IN) :: field_type

#ifdef ACC 
  !$acc routine (intper_fields) seq
#endif

#ifdef ACC  
  !$acc  parallel loop &
  !$acc& private(X_X_loop,X_Y_loop,X_Z_loop,V_X_loop, &
  !$acc& V_Y_loop,V_Z_loop,gam_loop)
#endif ACC

#ifdef OMP
  !$omp  parallel do &
  !$omp& default(none) &
  !$omp& firstprivate(dt,nRE,t_steps,field_type) &
  !$omp& shared(X_X,X_Y,X_Z,V_X,V_Y,V_Z,gam,B_X,B_Y,B_Z,E_X,E_Y,E_Z) &
  !$omp& private(X_X_loop,X_Y_loop,X_Z_loop,V_X_loop,V_Y_loop,V_Z_loop,gam_loop, &
  !$omp& B_X_loop,B_Y_loop,B_Z_loop,E_X_loop,E_Y_loop,E_Z_loop,U_X,U_Y,U_Z, &
  !$omp& cross_X,cross_Y,cross_z,U_hs_X,U_hs_Y,U_hs_Z,tau_X,tau_Y,tau_Z, &
  !$omp& up_X,up_Y,up_Z,gp,sigma,us,t_X,t_Y,t_Z,s) 
#endif OMP
  do pp=1,nRE

    X_X_loop=X_X(pp)
    X_Y_loop=X_Y(pp)
    X_Z_loop=X_Z(pp)
 
    V_X_loop=V_X(pp)
    V_Y_loop=V_Y(pp)
    V_Z_loop=V_Z(pp)
 
    B_X_loop=B_X(pp)
    B_Y_loop=B_Y(pp)
    B_Z_loop=B_Z(pp)
 
    E_X_loop=E_X(pp)
    E_Y_loop=E_Y(pp)
    E_Z_loop=E_Z(pp)
 
    gam_loop=gam(pp)
 
    !! Initial half step
    X_X_loop = X_X_loop+dt/2*V_X_loop
    X_Y_loop = X_Y_loop+dt/2*V_Y_loop
    X_Z_loop = X_Z_loop+dt/2*V_Z_loop
 
 !   !! Main iteration loop
#ifdef ACC  
    !$acc loop seq
#endif ACC
    do it=1,t_steps
 
       U_X = gam_loop*V_X_loop
       U_Y = gam_loop*V_Y_loop
       U_Z = gam_loop*V_Z_loop
 
#ifdef PSPLINE
       if (field_type.eq.'PSPLINE') then
            call interp_fields(X_X_loop,X_Y_loop,B_X_loop,B_Y_loop,B_Z_loop, &
              E_X_loop,E_Y_loop,E_Z_loop)
       endif
#endif PSPLINE

       ! LEAP-FROG SCHEME FOR LORENTZ FORCE !
 
       cross_X=V_Y_loop*B_Z_loop-V_Z_loop*B_Y_loop
       cross_Y=V_Z_loop*B_X_loop-V_X_loop*B_Z_loop
       cross_Z=V_X_loop*B_Y_loop-V_Y_loop*B_X_loop
 
       !write(6,*) 'vcrossB',cross_X,cross_Y,cross_Z
 
       U_hs_X = U_X + 0.5*dt*(E_X_loop +cross_X)
       U_hs_Y = U_Y + 0.5*dt*(E_Y_loop +cross_Y)
       U_hs_Z = U_Z + 0.5*dt*(E_Z_loop +cross_Z)
 
       !write(6,*) 'half step',0.5*a*cross_X(pp),0.5*a*cross_Y(pp),0.5*a*cross_Z(pp)
 
       tau_X = 0.5*dt*B_X_loop
       tau_Y = 0.5*dt*B_Y_loop
       tau_Z = 0.5*dt*B_Z_loop
 
       up_X = U_hs_X + 0.5*dt*E_X_loop
       up_Y = U_hs_Y + 0.5*dt*E_Y_loop
       up_Z = U_hs_Z + 0.5*dt*E_Z_loop
 
       gp = SQRT( 1.0 + up_X*up_X+up_Y*up_Y+up_Z*up_Z )
 
       sigma = gp*gp - (tau_X*tau_X+tau_Y*tau_Y+tau_Z*tau_Z)
 
       us = up_X*tau_X+up_Y*tau_Y+up_Z*tau_Z
       ! variable 'u^*' in Vay, J.-L. PoP (2008)
 
       gam_loop = SQRT( 0.5*(sigma + SQRT(sigma*sigma + &
       4.0*(tau_X*tau_X+tau_Y*tau_Y+tau_Z*tau_Z + us*us))) )
 
       t_X = tau_X/gam_loop
       t_Y = tau_Y/gam_loop
       t_Z = tau_Z/gam_loop
 
       s = 1.0/(1.0 + t_X*t_X+t_Y*t_Y+t_Z*t_Z)
       ! variable 's' in Vay, J.-L. PoP (2008)
 
       cross_X=up_Y*t_Z-up_Z*t_Y
       cross_Y=up_Z*t_X-up_X*t_Z
       cross_Z=up_X*t_Y-up_Y*t_X
 
       U_X = s*(up_X + (up_X*t_X+up_Y*t_Y+up_Z*t_Z)*t_X + cross_X)
       U_Y = s*(up_Y + (up_X*t_X+up_Y*t_Y+up_Z*t_Z)*t_Y + cross_Y)
       U_Z = s*(up_Z + (up_X*t_X+up_Y*t_Y+up_Z*t_Z)*t_Z + cross_Z)
       ! LEAP-FROG SCHEME FOR LORENTZ FORCE !
 
       V_X_loop = U_X/gam_loop
       V_Y_loop = U_Y/gam_loop
       V_Z_loop = U_Z/gam_loop
 
       X_X_loop = X_X_loop + dt*V_X_loop
       X_Y_loop = X_Y_loop + dt*V_Y_loop
       X_Z_loop = X_Z_loop + dt*V_Z_loop
 
       !write(data_write,*) 'V: ',V_X_loop*v_norm,V_Y_loop*v_norm,V_Z_loop*v_norm
       !write(data_write,*) 'X: ',X_X_loop*x_norm,X_Y_loop*x_norm,X_Z_loop*x_norm
 
    end do
!#ifdef ACC  
!    !$acc end loop seq
!#endif ACC
 
    X_X(pp)=X_X_loop
    X_Y(pp)=X_Y_loop
    X_Z(pp)=X_Z_loop
 
    V_X(pp)=V_X_loop
    V_Y(pp)=V_Y_loop
    V_Z(pp)=V_Z_loop
 
    gam(pp)=gam_loop
    
  end do
#ifdef ACC  
  !$acc end parallel loop
#endif ACC

#ifdef OMP
  !$omp end parallel do
#endif

end subroutine FO_push

end module pusher_gpu