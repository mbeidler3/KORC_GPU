!-----------------------------------------------------------------------
!     $Id: get_fields.f90 7623 2022-01-18 03:55:24Z tbechtel $
!     get arbitrary lagr_quad field at a given location, either
!     including the equilibrium field or not, and including gradient
!     or not
!-----------------------------------------------------------------------
      MODULE get_field_mod
      USE input
      USE global
      USE fields
      USE map_mod
      IMPLICIT NONE

      INTERFACE get_field
        MODULE PROCEDURE                                                &
     &    get_field_woeq,get_field_weq,get_field_eq,get_field_wderiv,   &
     &    get_field_eq_wderiv,get_field_woeq_wderiv
      END INTERFACE

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram get_field_woeq
!     compute the field strength at the specified (x,y,phi) coord
!     we presume that x and y are in rb_cel(1)
!     output does not include equilibrium quantities
!-----------------------------------------------------------------------
      SUBROUTINE get_field_woeq(coord,laq,outfld)
      USE physdat
      TYPE(lagr_quad_type) :: laq
      REAL(r8), DIMENSION(:), INTENT(IN) :: coord
      REAL(r8), DIMENSION(:), INTENT(OUT) :: outfld

      REAL(r8) :: x,y,gcosv,gsinv,angle
      COMPLEX(r8) :: ff(laq%nqty,nmodes_total)
      INTEGER(i4) :: imode
!-----------------------------------------------------------------------
!     Check for consistent inputs
!-----------------------------------------------------------------------
      IF (SIZE(coord)<3)                                                &
     &  CALL nim_stop('get_field: coord must contain x, y, and phi')
      IF (SIZE(outfld)/=laq%nqty)                                       &
     &  CALL nim_stop('get_field: output does not match intended size')
      IF (laq%nfour/=nmodes_total)                                      &
     &  CALL nim_stop('get_field: processor must have all modes')
!-----------------------------------------------------------------------
!     Find the field within the block
!-----------------------------------------------------------------------
      x=coord(1)
      y=coord(2)
      angle=coord(3)
      IF (geom=='lin') angle=twopi*angle/per_length
      CALL lagr_quad_eval(laq,x,y,0_i4,ff)
!-----------------------------------------------------------------------
!     Sum the perturbed fields over number of Fourier modes:
!       ff(1,1)-------r-comp, 1st mode
!       ff(2,1)-------z-comp, 1st mode
!       ff(3,1)-------phi-comp, 1st mode
!       ff(1,2)-------r-comp, 2nd mode
!       ff(2,2)-------z-comp, 2nd mode
!       ff(3,2)-------phi-comp, 2nd mode
!       ff(1,3)-------r-comp, 3rd mode
!       ...
!
!       F(R,Z,phi)=sum_over_modes[ff_n(R,Z)*exp( i*n*phi)
!                                 +conjg(ff_n(R,Z)*exp(-i*n*phi)]
!-----------------------------------------------------------------------
      outfld=0._r8
      ! Sum Fourier components.
      DO imode=1,nmodes_total
        IF (keff_total(imode)==0) THEN
          ! Ensure n=0 terms are purely real
          outfld(:)=outfld(:)+REAL(ff(:,imode),r8)
        ELSE
          ! n/=0 terms
          gcosv=COS(angle*keff_total(imode))
          gsinv=SIN(angle*keff_total(imode))
          outfld(:)=outfld(:)+2._r8*(REAL(  ff(:,imode),r8)*gcosv       &
     &                               -AIMAG(ff(:,imode)   )*gsinv)
        END IF
      END DO

      RETURN
      END SUBROUTINE get_field_woeq
!-----------------------------------------------------------------------
!     subprogram get_field_weq
!     compute the field strength at the specified (x,y,phi) coord with
!       equilibrium field included
!     we presume that x and y are in rb_cel(1)
!     the toroidal component of the eq field is modified by bigr_exp
!-----------------------------------------------------------------------
      SUBROUTINE get_field_weq(coord,laqeq,laq,outfld,bigr,bigr_exp)
      USE physdat
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laqeq
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      REAL(r8), INTENT(IN) :: bigr_exp
      REAL(r8), DIMENSION(:), INTENT(IN) :: coord
      REAL(r8), INTENT(OUT) :: bigr
      REAL(r8), DIMENSION(:), INTENT(OUT) :: outfld

      REAL(r8) :: x,y,gcosv,gsinv,angle,rzf(2)
      COMPLEX(r8) :: ff(laq%nqty,nmodes_total)
      INTEGER(i4) :: imode
!-----------------------------------------------------------------------
!     Check for consistent inputs
!-----------------------------------------------------------------------
      IF (SIZE(coord)<3)                                                &
     &  CALL nim_stop('get_field: coord must contain x, y, and phi')
      IF (SIZE(outfld)/=laq%nqty)                                       &
     &  CALL nim_stop('get_field: output does not match intended size')
      IF (laqeq%nqty/=laq%nqty)                                         &
     &  CALL nim_stop('get_field: eq and perturbed fields do not match')
      IF (laq%nfour/=nmodes_total)                                      &
     &  CALL nim_stop('get_field: processor must have all modes')
!-----------------------------------------------------------------------
!     Get coordinates
!-----------------------------------------------------------------------
      x=coord(1)
      y=coord(2)
      angle=coord(3)
      IF (geom=='lin') THEN
        angle=twopi*angle/per_length
        bigr=1.
      ELSE
        CALL lagr_quad_eval(rb_cel(1)%rz,x,y,0_i4,rzf)
        bigr=rzf(1)
      END IF
!-----------------------------------------------------------------------
!     Find the equilibrium field within the block
!-----------------------------------------------------------------------
      CALL lagr_quad_eval(laqeq,x,y,0_i4,outfld)
      IF(SIZE(outfld)==3) outfld(3) = outfld(3)*(bigr**bigr_exp)
!-----------------------------------------------------------------------
!     Sum the perturbed fields over number of Fourier modes:
!       ff(1,1)-------r-comp, 1st mode
!       ff(2,1)-------z-comp, 1st mode
!       ff(3,1)-------phi-comp, 1st mode
!       ff(1,2)-------r-comp, 2nd mode
!       ff(2,2)-------z-comp, 2nd mode
!       ff(3,2)-------phi-comp, 2nd mode
!       ff(1,3)-------r-comp, 3rd mode
!       ...
!
!       F(R,Z,phi)=sum_over_modes[ff_n(R,Z)*exp( i*n*phi)
!                                 +conjg(ff_n(R,Z)*exp(-i*n*phi)]
!-----------------------------------------------------------------------
      CALL lagr_quad_eval(laq,x,y,0_i4,ff)
      ! Sum Fourier components.
      DO imode=1,nmodes_total
        IF (keff_total(imode)==0) THEN
          ! Ensure n=0 terms are purely real
          outfld(:)=outfld(:)+REAL(ff(:,imode),r8)
        ELSE
          ! n/=0 terms
          gcosv=COS(angle*keff_total(imode))
          gsinv=SIN(angle*keff_total(imode))
          outfld(:)=outfld(:)+2._r8*(REAL(  ff(:,imode),r8)*gcosv       &
     &                                 -AIMAG(ff(:,imode)   )*gsinv)
        END IF
      END DO

      RETURN
      END SUBROUTINE get_field_weq
!-----------------------------------------------------------------------
!     subprogram get_field_wderiv
!     compute the field strength and gradient at the specified
!       (x,y,phi) coord with equilibrium field included
!     we presume that x and y are in rb_cel(1)
!     the toroidal component of the eq field is modified by bigr_exp
!     syntax for the output is:
!       outfld(:,0)=fields
!       outfld(:,1:3)=gradient
!     note: gradient has the 1/R factor in d/dphi
!-----------------------------------------------------------------------
      SUBROUTINE get_field_wderiv(coord,laqeq,laq,outfld,bigr,bigr_exp)
      USE physdat
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laqeq
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      REAL(r8), INTENT(IN) :: bigr_exp
      REAL(r8), DIMENSION(:), INTENT(IN) :: coord
      REAL(r8), INTENT(OUT) :: bigr
      REAL(r8), DIMENSION(:,0:), INTENT(OUT) :: outfld

      REAL(r8) :: x,y,gcosv,gsinv,angle
      REAL(r8) :: rzf(2),drzdx(2),drzdy(2)
      REAL(r8) :: dxdr,dxdz,dydr,dydz
      COMPLEX(r8) :: ff(laq%nqty,nmodes_total)
      COMPLEX(r8) :: dffx(laq%nqty,nmodes_total)
      COMPLEX(r8) :: dffy(laq%nqty,nmodes_total)
      COMPLEX(r8) :: derivs(laq%nqty,nmodes_total,3)
      REAL(r8) :: dfqx(laq%nqty)
      REAL(r8) :: dfqy(laq%nqty)
      INTEGER(i4) :: imode
!-----------------------------------------------------------------------
!     Check for consistent inputs
!-----------------------------------------------------------------------
      IF (SIZE(coord)<3)                                                &
     &  CALL nim_stop('get_field: coord must contain x, y, and phi')
      IF (SIZE(outfld,1)/=laq%nqty)                                     &
     &  CALL nim_stop('get_field: output does not match intended size')
      IF (SIZE(outfld,2)/=4)                                            &
     &  CALL nim_stop('get_field: output does not match intended size')
      IF (laqeq%nqty/=laq%nqty)                                         &
     &  CALL nim_stop('get_field: eq and perturbed fields do not match')
      IF (laq%nfour/=nmodes_total)                                      &
     &  CALL nim_stop('get_field: processor must have all modes')
!-----------------------------------------------------------------------
!     Get coordinates and metric
!-----------------------------------------------------------------------
      x=coord(1)
      y=coord(2)
      angle=coord(3)
      CALL lagr_quad_eval(rb_cel(1)%rz,x,y,1_i4,rzf,drzdx,drzdy)
      IF (geom=='tor') THEN
        bigr=rzf(1)
      ELSE
        bigr=1._r8
        angle=twopi*angle/per_length
      END IF
      CALL calcderivs(drzdx,drzdy,dxdr,dxdz,dydr,dydz)
!-----------------------------------------------------------------------
!     Find the equilibrium field within the block
!-----------------------------------------------------------------------
      outfld=0._r8
      derivs=0._r8
      CALL lagr_quad_eval(laqeq,x,y,1_i4,outfld(:,0),dfqx,dfqy)
      derivs(:,1,1)=(dfqx*dxdr+dfqy*dydr)
      derivs(:,1,2)=(dfqx*dxdz+dfqy*dydz)
      IF(SIZE(outfld,1)==3) THEN
        outfld(3,0) = outfld(3,0)*(bigr**bigr_exp)
        derivs(3,1,1)=(derivs(3,1,1)                                    &
     &               +bigr_exp*bigr**(-bigr_exp-1)*outfld(3,0))         &
     &                *(bigr**bigr_exp)
        derivs(3,1,2)=derivs(3,1,2)*(bigr**bigr_exp)
      ENDIF

!-----------------------------------------------------------------------
!     Sum the perturbed fields over number of Fourier modes:
!	      ff(1,1)-------r-comp, 1st mode
!	      ff(2,1)-------z-comp, 1st mode
!	      ff(3,1)-------phi-comp, 1st mode
!	      ff(1,2)-------r-comp, 2nd mode
!	      ff(2,2)-------z-comp, 2nd mode
!	      ff(3,2)-------phi-comp, 2nd mode
!	      ff(1,3)-------r-comp, 3rd mode
!       ...
!
!       F(R,Z,phi)=sum_over_modes[ff_n(R,Z)*exp( i*n*phi)
!                                 +conjg(ff_n(R,Z)*exp(-i*n*phi)]
!
!     Derivatives:
! outfield(1,1)=∂_R(F_R), outfield(1,2)=∂_Z(F_R), outfield(1,3)=∂_p(F_R)
! outfield(2,1)=∂_R(F_Z), outfield(2,2)=∂_Z(F_Z), outfield(2,3)=∂_p(F_Z)
! outfield(3,1)=∂_R(F_p), outfield(3,2)=∂_Z(F_p), outfield(3,3)=∂_p(F_p)
!-----------------------------------------------------------------------
      CALL lagr_quad_eval(laq,x,y,1_i4,ff,dffx,dffy)

      ! Derivatives: d/dR,d/dZ,d/dphi in last index
      derivs(:,:,1)=derivs(:,:,1)+dffx*dxdr+dffy*dydr
      derivs(:,:,2)=derivs(:,:,2)+dffx*dxdz+dffy*dydz
      DO imode=1,nmodes_total
        derivs(:,imode,3)=(0,1)*keff_total(imode)*ff(:,imode)/bigr
      END DO

      ! Sum Fourier components.
      DO imode=1,nmodes_total
        IF (keff_total(imode)==0) THEN
          ! Ensure n=0 terms are purely real
          outfld(:,0)=outfld(:,0)+REAL(ff(:,imode),r8)
          outfld(:,1:3)=outfld(:,1:3)+REAL(derivs(:,imode,:),r8)
        ELSE
          ! n/=0 terms
          gcosv=COS(angle*keff_total(imode))
          gsinv=SIN(angle*keff_total(imode))
          outfld(:,0)=outfld(:,0)+2._r8*(REAL(  ff(:,imode),r8)*gcosv &
     &                                   -AIMAG(ff(:,imode)   )*gsinv)
          outfld(:,1:3)=outfld(:,1:3)                                 &
     &                     +2._r8*(  REAL(derivs(:,imode,:),r8)*gcosv &
     &                             -AIMAG(derivs(:,imode,:)   )*gsinv)
        END IF
      END DO

      RETURN

      CONTAINS
        SUBROUTINE calcderivs(drzdx,drzdy,dxdr,dxdz,dydr,dydz)
        REAL(r8), DIMENSION(2), INTENT(IN) :: drzdx,drzdy
        REAL(r8), INTENT(OUT) :: dxdr,dxdz,dydr,dydz
        REAL(r8) :: jac
        ! find the jacobian, then perform a matrix inversion.
        jac = drzdx(1)*drzdy(2)-drzdx(2)*drzdy(1)
        dxdr= drzdy(2)/jac; dydr=-drzdx(2)/jac
        dxdz=-drzdy(1)/jac; dydz= drzdx(1)/jac
        RETURN
        END SUBROUTINE calcderivs

      END SUBROUTINE get_field_wderiv
!-----------------------------------------------------------------------
!     subprogram get_field_n0
!     get the n=0 (perturbed) field strength at the specified
!       (x,y) coord with equilibrium field included
!     we presume that x and y are in rb_cel(1)
!     the toroidal component of the eq field is modified by bigr_exp
!-----------------------------------------------------------------------
      SUBROUTINE get_field_n0(coord,laqeq,laq,outfld,bigr,bigr_exp)
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laqeq
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      REAL(r8), INTENT(IN) :: bigr_exp
      REAL(r8), DIMENSION(:), INTENT(IN) :: coord
      REAL(r8), INTENT(OUT) :: bigr
      REAL(r8), DIMENSION(:), INTENT(OUT) :: outfld

      REAL(r8) :: x,y,rzf(2)
      COMPLEX(r8) :: ff(laq%nqty,laq%nfour)
      INTEGER(i4) :: imode
!-----------------------------------------------------------------------
!     Check for consistent inputs
!-----------------------------------------------------------------------
      IF (SIZE(coord)<2)                                                &
     &  CALL nim_stop('get_field: coord must contain x and  y')
      IF (SIZE(outfld)/=laq%nqty)                                       &
     &  CALL nim_stop('get_field: output does not match intended size')
      IF (laqeq%nqty/=laq%nqty)                                         &
     &  CALL nim_stop('get_field: eq and perturbed fields do not match')
      IF (.NOT.nonlinear.AND.lin_nmax>=lin_nmodes)                      &
     &  CALL nim_stop('get_field: n=0 mode does not exist')
!-----------------------------------------------------------------------
!     Get coordinates
!-----------------------------------------------------------------------
      x=coord(1)
      y=coord(2)
      IF (geom=='tor') THEN
        CALL lagr_quad_eval(rb_cel(1)%rz,x,y,0_i4,rzf)
        bigr=rzf(1)
      ELSE
        bigr=1.
      END IF
!-----------------------------------------------------------------------
!     Find the equilibrium field within the block
!-----------------------------------------------------------------------
      CALL lagr_quad_eval(laqeq,x,y,0_i4,outfld)
      IF(SIZE(outfld)==3) outfld(3) = outfld(3)*(bigr**bigr_exp)
!-----------------------------------------------------------------------
!     Add the perturbed n=0 field:
!	      ff(1,1)-------r-comp, 1st mode
!	      ff(2,1)-------z-comp, 1st mode
!	      ff(3,1)-------phi-comp, 1st mode
!-----------------------------------------------------------------------
      CALL lagr_quad_eval(laq,x,y,0_i4,ff)
      DO imode=1,laq%nfour
        ! Ensure n=0 term is purely real
        IF (keff_total(imode)==0) THEN
          outfld=outfld+REAL(ff(:,imode),r8)
          RETURN
        END IF
      END DO

      RETURN
      END SUBROUTINE get_field_n0
!-----------------------------------------------------------------------
!     subprogram get_field_eq
!     get the equilibrium field strength at the specified (x,y) coord
!     we presume that x and y are in rb_cel(1)
!     the toroidal component is modified by bigr_exp
!-----------------------------------------------------------------------
      SUBROUTINE get_field_eq(coord,laqeq,outfld,bigr,bigr_exp)
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laqeq
      REAL(r8), INTENT(IN) :: bigr_exp
      REAL(r8), DIMENSION(:), INTENT(IN) :: coord
      REAL(r8), INTENT(OUT) :: bigr
      REAL(r8), DIMENSION(:), INTENT(OUT) :: outfld

      REAL(r8) :: x,y,rzf(2)
!-----------------------------------------------------------------------
!     Check for consistent inputs
!-----------------------------------------------------------------------
      IF (SIZE(coord)<2)                                                &
     &   CALL nim_stop('get_field: coord must contain x and  y')
      IF (SIZE(outfld)/=laqeq%nqty)                                     &
     &   CALL nim_stop('get_field: output does not match intended size')
!-----------------------------------------------------------------------
!     Get coordinates
!-----------------------------------------------------------------------
      x=coord(1)
      y=coord(2)
      IF (geom=='tor') THEN
        CALL lagr_quad_eval(rb_cel(1)%rz,x,y,0_i4,rzf)
        bigr=rzf(1)
      ELSE
        bigr=1.
      END IF
!-----------------------------------------------------------------------
!     Find the equilibrium field within the block
!-----------------------------------------------------------------------
      CALL lagr_quad_eval(laqeq,x,y,0_i4,outfld)
      IF(SIZE(outfld)==3) outfld(3) = outfld(3)*(bigr**bigr_exp)

      RETURN
      END SUBROUTINE get_field_eq
!-----------------------------------------------------------------------
!     subprogram get_field_eq_wderiv
!     compute the equilibrium field strength and derivatives at the
!       specified (x,y) coord
!     we presume that x and y are in rb_cel(1)
!     the toroidal component of the eq field is modified by bigr_exp
!     syntax for the output is:
!       outfld(:,0)=eq fields
!       outfld(:,1:3)=eq derivatives
!-----------------------------------------------------------------------
      SUBROUTINE get_field_eq_wderiv(coord,laqeq,outfld,bigr,bigr_exp)
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laqeq
      REAL(r8), INTENT(IN) :: bigr_exp
      REAL(r8), DIMENSION(:), INTENT(IN) :: coord
      REAL(r8), INTENT(OUT) :: bigr
      REAL(r8), DIMENSION(:,0:), INTENT(OUT) :: outfld

      REAL(r8) :: x,y
      REAL(r8) :: rzf(2),drzdx(2),drzdy(2)
      REAL(r8) :: dxdr,dxdz,dydr,dydz
      REAL(r8) :: dfqx(laqeq%nqty)
      REAL(r8) :: dfqy(laqeq%nqty)
!-----------------------------------------------------------------------
!     Check for consistent inputs
!-----------------------------------------------------------------------
      IF (SIZE(coord)<2)                                                &
     &   CALL nim_stop('get_field: coord must contain x and  y')
      IF (SIZE(outfld,1)/=laqeq%nqty)                                   &
     &   CALL nim_stop('get_field: output does not match intended size')
      IF (SIZE(outfld,2)/=4)                                            &
     &   CALL nim_stop('get_field: output does not match intended size')
!-----------------------------------------------------------------------
!     Get coordinates
!-----------------------------------------------------------------------
      x=coord(1)
      y=coord(2)
      CALL lagr_quad_eval(rb_cel(1)%rz,x,y,1_i4,rzf,drzdx,drzdy)
      IF (geom=='tor') THEN
        bigr=rzf(1)
      ELSE
        bigr=1.
      END IF
      CALL calcderivs(drzdx,drzdy,dxdr,dxdz,dydr,dydz)
!-----------------------------------------------------------------------
!     Find the equilibrium field within the block
!-----------------------------------------------------------------------
      CALL lagr_quad_eval(laqeq,x,y,1_i4,outfld(:,0),dfqx,dfqy)
      ! n=0 derivatives: d/dR,d/dZ,d/dphi in last index
      outfld(:,1)=(dfqx*dxdr+dfqy*dydr)
      outfld(:,2)=(dfqx*dxdz+dfqy*dydz)
      outfld(:,3)=0.
      IF(SIZE(outfld,1)==3) THEN
        outfld(3,0) = outfld(3,0)*(bigr**bigr_exp)
        outfld(3,1)=(outfld(3,1)                                        &
     &               +bigr_exp*bigr**(-bigr_exp-1)*outfld(3,0))         &
     &                *(bigr**bigr_exp)
        outfld(3,2)=outfld(3,2)*(bigr**bigr_exp)
      ENDIF

      RETURN

      CONTAINS
        SUBROUTINE calcderivs(drzdx,drzdy,dxdr,dxdz,dydr,dydz)
        REAL(r8), DIMENSION(2), INTENT(IN) :: drzdx,drzdy
        REAL(r8), INTENT(OUT) :: dxdr,dxdz,dydr,dydz
        REAL(r8) :: jac
        ! find the jacobian, then perform a matrix inversion.
        jac = drzdx(1)*drzdy(2)-drzdx(2)*drzdy(1)
        dxdr= drzdy(2)/jac; dydr=-drzdx(2)/jac
        dxdz=-drzdy(1)/jac; dydz= drzdx(1)/jac
        RETURN
        END SUBROUTINE calcderivs

      END SUBROUTINE get_field_eq_wderiv
!-----------------------------------------------------------------------
!     subprogram get_field_woeq_wderiv
!     compute the field strength and gradient at the specified
!       (x,y,phi) coord without the equilibrium field included
!     we presume that x and y are in rb_cel(1)
!     the toroidal component of the eq field is modified by bigr_exp
!     syntax for the output is:
!       outfld(:,0)=fields
!       outfld(:,1:3)=gradient
!     note: gradient has the 1/R factor in d/dphi
!-----------------------------------------------------------------------
      SUBROUTINE get_field_woeq_wderiv(coord,laq,outfld,bigr,bigr_exp)
      USE physdat
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      REAL(r8), INTENT(IN) :: bigr_exp
      REAL(r8), DIMENSION(:), INTENT(IN) :: coord
      REAL(r8), INTENT(OUT) :: bigr
      REAL(r8), DIMENSION(:,0:), INTENT(OUT) :: outfld

      REAL(r8) :: x,y,gcosv,gsinv,angle
      REAL(r8) :: rzf(2),drzdx(2),drzdy(2)
      REAL(r8) :: dxdr,dxdz,dydr,dydz
      COMPLEX(r8) :: ff(laq%nqty,nmodes_total)
      COMPLEX(r8) :: dffx(laq%nqty,nmodes_total)
      COMPLEX(r8) :: dffy(laq%nqty,nmodes_total)
      COMPLEX(r8) :: derivs(laq%nqty,nmodes_total,3)
      INTEGER(i4) :: imode
!-----------------------------------------------------------------------
!     Check for consistent inputs
!-----------------------------------------------------------------------
      IF (SIZE(coord)<3)                                                &
     &  CALL nim_stop('get_field: coord must contain x, y, and phi')
      IF (SIZE(outfld,1)/=laq%nqty)                                     &
     &  CALL nim_stop('get_field: output does not match intended size')
      IF (SIZE(outfld,2)/=4)                                            &
     &  CALL nim_stop('get_field: output does not match intended size')
      IF (laq%nfour/=nmodes_total)                                      &
     &  CALL nim_stop('get_field: processor must have all modes')
!-----------------------------------------------------------------------
!     Get coordinates and metric
!-----------------------------------------------------------------------
      x=coord(1)
      y=coord(2)
      angle=coord(3)
      CALL lagr_quad_eval(rb_cel(1)%rz,x,y,1_i4,rzf,drzdx,drzdy)
      IF (geom=='tor') THEN
        bigr=rzf(1)
      ELSE
        bigr=1._r8
        angle=twopi*angle/per_length
      END IF
      CALL calcderivs(drzdx,drzdy,dxdr,dxdz,dydr,dydz)
!-----------------------------------------------------------------------
!     Sum the perturbed fields over number of Fourier modes:
!	      ff(1,1)-------r-comp, 1st mode
!	      ff(2,1)-------z-comp, 1st mode
!	      ff(3,1)-------phi-comp, 1st mode
!	      ff(1,2)-------r-comp, 2nd mode
!	      ff(2,2)-------z-comp, 2nd mode
!	      ff(3,2)-------phi-comp, 2nd mode
!	      ff(1,3)-------r-comp, 3rd mode
!       ...
!
!       F(R,Z,phi)=sum_over_modes[ff_n(R,Z)*exp( i*n*phi)
!                                 +conjg(ff_n(R,Z)*exp(-i*n*phi)]
!
!     Derivatives:
! outfield(1,1)=∂_R(F_R), outfield(1,2)=∂_Z(F_R), outfield(1,3)=∂_p(F_R)
! outfield(2,1)=∂_R(F_Z), outfield(2,2)=∂_Z(F_Z), outfield(2,3)=∂_p(F_Z)
! outfield(3,1)=∂_R(F_p), outfield(3,2)=∂_Z(F_p), outfield(3,3)=∂_p(F_p)
!-----------------------------------------------------------------------
      outfld=0._r8
      derivs=0._r8
      CALL lagr_quad_eval(laq,x,y,1_i4,ff,dffx,dffy)

      ! Derivatives: d/dR,d/dZ,d/dphi in last index
      derivs(:,:,1)=derivs(:,:,1)+dffx*dxdr+dffy*dydr
      derivs(:,:,2)=derivs(:,:,2)+dffx*dxdz+dffy*dydz
      DO imode=1,nmodes_total
        derivs(:,imode,3)=(0,1)*keff_total(imode)*ff(:,imode)/bigr
      END DO

      ! Sum Fourier components.
      DO imode=1,nmodes_total
        IF (keff_total(imode)==0) THEN
          ! Ensure n=0 terms are purely real
          outfld(:,0)=outfld(:,0)+REAL(ff(:,imode),r8)
          outfld(:,1:3)=outfld(:,1:3)+REAL(derivs(:,imode,:),r8)
        ELSE
          ! n/=0 terms
          gcosv=COS(angle*keff_total(imode))
          gsinv=SIN(angle*keff_total(imode))
          outfld(:,0)=outfld(:,0)+2._r8*(REAL(  ff(:,imode),r8)*gcosv   &
     &                                   -AIMAG(ff(:,imode)   )*gsinv)
          outfld(:,1:3)=outfld(:,1:3)                                   &
     &                     +2._r8*(  REAL(derivs(:,imode,:),r8)*gcosv   &
     &                             -AIMAG(derivs(:,imode,:)   )*gsinv)
        END IF
      END DO

      RETURN

      CONTAINS
        SUBROUTINE calcderivs(drzdx,drzdy,dxdr,dxdz,dydr,dydz)
        REAL(r8), DIMENSION(2), INTENT(IN) :: drzdx,drzdy
        REAL(r8), INTENT(OUT) :: dxdr,dxdz,dydr,dydz
        REAL(r8) :: jac
        ! find the jacobian, then perform a matrix inversion.
        jac = drzdx(1)*drzdy(2)-drzdx(2)*drzdy(1)
        dxdr= drzdy(2)/jac; dydr=-drzdx(2)/jac
        dxdz=-drzdy(1)/jac; dydz= drzdx(1)/jac
        RETURN
        END SUBROUTINE calcderivs

      END SUBROUTINE get_field_woeq_wderiv
!-----------------------------------------------------------------------
!     subprogram get_field_dm
!     compute the field strength, 1st and 2nd derivatives at the
!       specified (x,y,phi) coord with equilibrium field included
!     we presume that x and y are in rb_cel(1)
!     the toroidal component of the eq field is modified by bigr_exp
!     syntax for the output is:
!       outfld(:,0)=fields
!       outfld(:,1:3)=1st derivatives
!       outfld(:,4:9)=2dn derivatives
!
!     see comments in get_comp_field_dm for details on how the 2nd
!     order derivatives are computed
!-----------------------------------------------------------------------
      SUBROUTINE get_field_dm(coord,laqeq,laq,outfld,dmode,n0only)
      USE physdat
      USE math_tran
      TYPE(lagr_quad_2D_type), INTENT(INOUT), OPTIONAL :: laqeq
      TYPE(lagr_quad_type), INTENT(INOUT), OPTIONAL :: laq
      REAL(r8), DIMENSION(:), INTENT(IN), OPTIONAL :: coord
      REAL(r8), DIMENSION(:,0:), INTENT(OUT), OPTIONAL :: outfld
      INTEGER(i4), INTENT(IN), OPTIONAL :: dmode
      LOGICAL, INTENT(IN), OPTIONAL :: n0only

      REAL(r8) :: x,y,gcosv,gsinv,angle,bigr
      REAL(r8) :: rzf(2),drzdx(2),drzdy(2)
      REAL(r8) :: d2rzdxx(2),d2rzdyy(2),d2rzdxy(2)
      REAL(r8) :: dxdr,dxdz,dydr,dydz
      INTEGER(i4) :: imode,dm
      LOGICAL :: n0

      REAL(r8), DIMENSION(:), ALLOCATABLE :: dfqx,dfqy,d2fqxx,d2fqxy,   &
     &                                       d2fqyy
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rhs,sol_real,sol_imag
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: mat
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: ff,dffx,dffy,d2ffxx,  &
     &                                            d2ffxy,d2ffyy
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: derivs,dderivs
!-----------------------------------------------------------------------
!     Check for consistent inputs
!-----------------------------------------------------------------------
      IF (.NOT. PRESENT(coord))                                         &
     &  CALL nim_stop('get_field: coord must be passed in.')
      IF (.NOT. (PRESENT(laqeq) .OR. PRESENT(laq)))                     &
     &  CALL nim_stop('get_field: equilibrium or perturbed field is'//  &
     &                ' not specified.')
      IF (.NOT. PRESENT(outfld))                                        &
     &  CALL nim_stop('get_field: outfld is not specified.')
      IF (SIZE(coord)<3)                                                &
     &  CALL nim_stop('get_field: coord must contain x, y, and phi')
      IF (PRESENT(laq) .AND. SIZE(outfld,1)/=laq%nqty)                  &
     &  CALL nim_stop('get_field: output does not match intended size')
      IF (PRESENT(laqeq) .AND. PRESENT(laq) .AND. laqeq%nqty/=laq%nqty) &
     &  CALL nim_stop('get_field: eq and perturbed fields do not match')
      IF (PRESENT(laq) .AND. laq%nfour/=nmodes_total)                                      &
     &  CALL nim_stop('get_field: processor must have all modes')
      dm=0
      IF (PRESENT(dmode)) dm=dmode
      n0=.FALSE.
      IF (PRESENT(n0only)) n0=n0only
!-----------------------------------------------------------------------
!     Allocate local arrays
!-----------------------------------------------------------------------
      IF (PRESENT(laqeq)) THEN
        ALLOCATE(dfqx(laqeq%nqty),dfqy(laqeq%nqty))
        ALLOCATE(d2fqxx(laqeq%nqty),d2fqxy(laqeq%nqty),                 &
     &           d2fqyy(laqeq%nqty))
        ALLOCATE(derivs(laqeq%nqty,nmodes_total,3))
        ALLOCATE(dderivs(laqeq%nqty,nmodes_total,6))
        ALLOCATE(mat(laqeq%nqty,nmodes_total,3,3),                      &
     &           rhs(laqeq%nqty,nmodes_total,3),                        &
     &           sol_real(laqeq%nqty,nmodes_total,3),                   &
     &           sol_imag(laqeq%nqty,nmodes_total,3))
      ENDIF
      IF (PRESENT(laq)) THEN
        ALLOCATE(ff(laq%nqty,nmodes_total),dffx(laq%nqty,nmodes_total), &
     &       dffy(laq%nqty,nmodes_total),d2ffxx(laq%nqty,nmodes_total), &
     &     d2ffxy(laq%nqty,nmodes_total),d2ffyy(laq%nqty,nmodes_total))
        IF (.NOT. ALLOCATED(derivs)) THEN
          ALLOCATE(derivs(laq%nqty,nmodes_total,3))
          ALLOCATE(dderivs(laq%nqty,nmodes_total,6))
          ALLOCATE(mat(laq%nqty,nmodes_total,3,3),                      &
     &             rhs(laq%nqty,nmodes_total,3),                        &
     &             sol_real(laq%nqty,nmodes_total,3),                   &
     &             sol_imag(laq%nqty,nmodes_total,3))
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     Get coordinates and metric
!-----------------------------------------------------------------------
      x=coord(1)
      y=coord(2)
      angle=coord(3)
      CALL lagr_quad_eval(rb_cel(1)%rz,x,y,2_i4,rzf,drzdx,drzdy,        &
     &                    d2rzdxx,d2rzdxy,d2rzdyy)
      IF (geom=='tor') THEN
        bigr=rzf(1)
      ELSE
        bigr=1._r8
        angle=twopi*angle/per_length
      END IF
      IF (dm>0) THEN
        CALL calcderivs (drzdx,drzdy,dxdr,dxdz,dydr,dydz)
        CALL calc2derivs(drzdx,drzdy,mat)
        CALL math_solve_3x3(mat,sol_real,rhs,'factor')
      ENDIF
!-----------------------------------------------------------------------
!     Find the equilibrium field within the block
!-----------------------------------------------------------------------
      outfld =0._r8
      derivs =0._r8
      dderivs=0._r8
      IF (PRESENT(laqeq)) THEN
        CALL lagr_quad_eval(laqeq,x,y,2_i4,outfld(:,0),dfqx,dfqy,       &
     &                      d2fqxx,d2fqxy,d2fqyy)
        IF (dm>0) THEN
          derivs(:,1,1)=(dfqx*dxdr+dfqy*dydr)
          derivs(:,1,2)=(dfqx*dxdz+dfqy*dydz)
        ENDIF
        IF (dm>1) THEN
          ! Solving for equilibrium part
          rhs(:,1,1)=d2fqxx-derivs(:,1,1)*d2rzdxx(1)-                   &
     &                      derivs(:,1,2)*d2rzdxx(2)
          rhs(:,1,2)=d2fqyy-derivs(:,1,1)*d2rzdyy(1)-                   &
     &                      derivs(:,1,2)*d2rzdyy(2)
          rhs(:,1,3)=d2fqxy-derivs(:,1,1)*d2rzdxy(1)-                   &
     &                      derivs(:,1,2)*d2rzdxy(2)
          CALL math_solve_3x3(mat,sol_real,rhs,'solve ')

          dderivs(:,1,1)=sol_real(:,1,1)
          dderivs(:,1,2)=sol_real(:,1,3)
          dderivs(:,1,4)=sol_real(:,1,2)
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     Sum the perturbed fields over number of Fourier modes:
!	      ff(1,1)-------r-comp, 1st mode
!	      ff(2,1)-------z-comp, 1st mode
!	      ff(3,1)-------phi-comp, 1st mode
!	      ff(1,2)-------r-comp, 2nd mode
!	      ff(2,2)-------z-comp, 2nd mode
!	      ff(3,2)-------phi-comp, 2nd mode
!	      ff(1,3)-------r-comp, 3rd mode
!       ...
!
!       F(R,Z,phi)=sum_over_modes[ff_n(R,Z)*exp( i*n*phi)
!                                 +conjg(ff_n(R,Z)*exp(-i*n*phi)]
!
!     Derivatives:
! outfield(1,1)=∂_R(F_R), outfield(1,2)=∂_Z(F_R), outfield(1,3)=∂_p(F_R)
! outfield(2,1)=∂_R(F_Z), outfield(2,2)=∂_Z(F_Z), outfield(2,3)=∂_p(F_Z)
! outfield(3,1)=∂_R(F_p), outfield(3,2)=∂_Z(F_p), outfield(3,3)=∂_p(F_p)
!-----------------------------------------------------------------------
      IF (PRESENT(laq)) THEN
        CALL lagr_quad_eval(laq,x,y,2_i4,ff,dffx,dffy,                  &
     &                      d2ffxx,d2ffxy,d2ffyy)

        IF (dm > 0) THEN
          ! Derivatives: d/dR,d/dZ,d/dphi in last index
          derivs(:,:,1)=derivs(:,:,1)+dffx*dxdr+dffy*dydr
          derivs(:,:,2)=derivs(:,:,2)+dffx*dxdz+dffy*dydz
          DO imode=1,nmodes_total
            derivs(:,imode,3)=(0,1)*keff_total(imode)*ff(:,imode)/bigr
          END DO
        ENDIF
        IF (dm > 1) THEN
          ! Second derivatives: d^2/dR^2, d^2/dRdZ, d^2/dRdphi
          !                     d^2/dZ^2, d^2/dZdphi
          !                     d^2/dphi^2  in last index
          ! Solving for real part
          rhs=0._r8
          rhs(:,:,1)=REAL(d2ffxx-derivs(:,:,1)*d2rzdxx(1)-              &
     &                           derivs(:,:,2)*d2rzdxx(2))
          rhs(:,:,2)=REAL(d2ffyy-derivs(:,:,1)*d2rzdyy(1)-              &
     &                           derivs(:,:,2)*d2rzdyy(2))
          rhs(:,:,3)=REAL(d2ffxy-derivs(:,:,1)*d2rzdxy(1)-              &
     &                           derivs(:,:,2)*d2rzdxy(2))
          CALL math_solve_3x3(mat,sol_real,rhs,'solve ')
          ! Solving for imaginary part
          rhs=0._r8
          rhs(:,:,1)=AIMAG(d2ffxx-derivs(:,:,1)*d2rzdxx(1)-             &
     &                            derivs(:,:,2)*d2rzdxx(2))
          rhs(:,:,2)=AIMAG(d2ffyy-derivs(:,:,1)*d2rzdyy(1)-             &
     &                            derivs(:,:,2)*d2rzdyy(2))
          rhs(:,:,3)=AIMAG(d2ffxy-derivs(:,:,1)*d2rzdxy(1)-             &
     &                            derivs(:,:,2)*d2rzdxy(2))
          CALL math_solve_3x3(mat,sol_imag,rhs,'solve ')

          dderivs(:,:,1)=dderivs(:,:,1)+sol_real(:,:,1)+                &
     &                            (0,1)*sol_imag(:,:,1)
          dderivs(:,:,2)=dderivs(:,:,2)+sol_real(:,:,3)+                &
     &                            (0,1)*sol_imag(:,:,3)
          dderivs(:,:,4)=dderivs(:,:,4)+sol_real(:,:,2)+                &
     &                            (0,1)*sol_imag(:,:,2)
          DO imode=1,nmodes_total
            dderivs(:,imode,3)=dderivs(:,imode,3)+                      &
     &                         (0,1)*keff_total(imode)*derivs(:,imode,1)/bigr
            dderivs(:,imode,5)=dderivs(:,imode,5)+                      &
     &                         (0,1)*keff_total(imode)*derivs(:,imode,2)/bigr
            dderivs(:,imode,6)=dderivs(:,imode,6)-                      &
     &                         ff(:,imode)*(keff_total(imode)/bigr)**2
          ENDDO
        ENDIF
      ENDIF

      ! Sum Fourier components.
      DO imode=1,nmodes_total
        IF (keff_total(imode)==0) THEN
          ! Ensure n=0 terms are purely real
          IF (PRESENT(laq))                                             &
     &      outfld(:,0)=outfld(:,0)+REAL(ff(:,imode),r8)
          IF (dm>0)                                                     &
     &      outfld(:,1:3)=outfld(:,1:3)+REAL( derivs(:,imode,:),r8)
          IF (dm>1)                                                     &
     &      outfld(:,4:9)=outfld(:,4:9)+REAL(dderivs(:,imode,:),r8)
        ELSE
          IF (PRESENT(laq).AND..NOT.n0) THEN
            ! n/=0 terms
            gcosv=COS(angle*keff_total(imode))
            gsinv=SIN(angle*keff_total(imode))
            outfld(:,0)=outfld(:,0)+2._r8*( REAL(ff(:,imode),r8)*gcosv  &
     &                                    -AIMAG(ff(:,imode)   )*gsinv)
            IF (dm>0) outfld(:,1:3)=outfld(:,1:3)                       &
     &                     +2._r8*(  REAL(derivs(:,imode,:),r8)*gcosv   &
     &                             -AIMAG(derivs(:,imode,:)   )*gsinv)
            IF (dm>1) outfld(:,4:9)=outfld(:,4:9)                       &
     &                     +2._r8*(  REAL(dderivs(:,imode,:),r8)*gcosv  &
     &                             -AIMAG(dderivs(:,imode,:)   )*gsinv)
          ENDIF
        ENDIF
      ENDDO

      IF (PRESENT(laqeq)) DEALLOCATE(dfqx,dfqy,d2fqxx,d2fqxy,d2fqyy)
      IF (PRESENT(laq)) DEALLOCATE(ff,dffx,dffy,d2ffxx,d2ffxy,d2ffyy)
      IF (dm>0) DEALLOCATE(derivs,dderivs,mat,rhs,sol_real,sol_imag)

      RETURN

      CONTAINS
        SUBROUTINE calcderivs(drzdx,drzdy,dxdr,dxdz,dydr,dydz)
        REAL(r8), DIMENSION(2), INTENT(IN) :: drzdx,drzdy
        REAL(r8), INTENT(OUT) :: dxdr,dxdz,dydr,dydz
        REAL(r8) :: jac
        ! find the jacobian, then perform a matrix inversion.
        jac = drzdx(1)*drzdy(2)-drzdx(2)*drzdy(1)
        dxdr= drzdy(2)/jac; dydr=-drzdx(2)/jac
        dxdz=-drzdy(1)/jac; dydz= drzdx(1)/jac
        RETURN
        END SUBROUTINE calcderivs

        SUBROUTINE calc2derivs(drzdx,drzdy,mat)
        REAL(r8), DIMENSION(2), INTENT(IN) :: drzdx,drzdy
        REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mat
        ! find the jacobian matrix for 2nd derivatives
        mat(:,:,1,1)=drzdx(1)**2
        mat(:,:,1,2)=drzdx(2)**2
        mat(:,:,1,3)=2._r8*drzdx(1)*drzdx(2)
        mat(:,:,2,1)=drzdy(1)**2
        mat(:,:,2,2)=drzdy(2)**2
        mat(:,:,2,3)=2._r8*drzdy(1)*drzdy(2)
        mat(:,:,3,1)=drzdx(1)*drzdy(1)
        mat(:,:,3,2)=drzdx(2)*drzdy(2)
        mat(:,:,3,3)=drzdx(1)*drzdy(2)+drzdx(2)*drzdy(1)
        RETURN
        END SUBROUTINE calc2derivs

      END SUBROUTINE get_field_dm
!-----------------------------------------------------------------------
!     subprogram get_comp_field_dm
!     return the fourier coefficient of the field and derivatives
!     Assume that x and y are in rb_cel(1)
!
!     Assume that the factors of R have been stripped out of
!     beq and jeq as done in eval_field: nimfield_init
!
!     output is array with all coefficient
!     outfld(:,n,0) imode=n fields, imode=0 stores eq (set to 0 in no eq)
!     outfld(:,n,1:3) imode=n 1st order derivatives
!     outfld(:,n,4:9) imode=n 2nd order derivatives
!-----------------------------------------------------------------------
! Fields (suprresing middle index = Fourier mode)
! outfield(i,0)=F_i
! 1st Derivatives
! outfield(i,1)=∂_R(F_i), outfield(i,2)=∂_Z(F_i), outfield(i,3)=∂_p(F_i)
! 2nd derivatives
! outfield(i,4)=∂_RR(F_i), outfield(i,5)=∂_RZ(F_i), outfield(i,6)=∂_Rp(F_i)
! outfield(i,7)=∂_ZZ(F_i), outfield(i,8)=∂_Zp(F_i), outfield(i,9)=∂_pp(F_i)
!
! derivs:
! derivs(i,1) = ∂_R(F_i), derivs(i,2) = ∂_Z(F_i), derivs(i,3) = ∂_p(F_i)
!
! dderivs:
! dderivs(:,1) = ∂_RR(F), dderivs(:,2) = ∂_RZ(F), dderivs(:,3) = ∂_Rp(F)
! dderivs(:,4) = ∂_ZZ(F), dderivs(:,5) = ∂_Zp(F), dderivs(:,6) = ∂_pp(F)
!
! to 2nd order derivatives are calculated from the equation
! r = C^{-1} [s-Dt]
!
! r^T = [∂_RR(F),∂_ZZ(F),∂_RZ(F)]
! s^T = [∂_xx(F),∂_yy(F),∂_xy(F)]
! t^T = [∂_R(F),∂_Z(F)]
!
!     | (∂_x(R))^2   (∂_x(Z))^2   2 ∂_x(R)∂_x(Z)            |
! C = | (∂_y(R))^2   (∂_y(Z))^2   2 ∂_y(R)∂_y(Z)            |
!     | ∂_x(R)∂_y(R) ∂_x(Z)∂_y(Z) ∂_x(R)∂_y(Z)+∂_y(R)∂_x(Z) |
!
!
!     | ∂_xx(R)  ∂_xx(Z) |
! D = | ∂_yy(R)  ∂_yy(Z) |
!     | ∂_xy(R)  ∂_xy(Z) |
!-----------------------------------------------------------------------
      SUBROUTINE get_comp_field_dm(coord,laq,outfld,dmode,laqeq)
      USE physdat
      USE math_tran
      REAL(r8), DIMENSION(:), INTENT(IN) :: coord
      TYPE(lagr_quad_type), INTENT(INOUT):: laq
      COMPLEX(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: outfld
      INTEGER(i4), INTENT(IN), OPTIONAL :: dmode
      TYPE(lagr_quad_2D_type), INTENT(INOUT), OPTIONAL :: laqeq

      REAL(r8) :: x,y,bigr
      REAL(r8) :: rzf(2),drzdx(2),drzdy(2)
      REAL(r8) :: d2rzdxx(2),d2rzdyy(2),d2rzdxy(2)
      REAL(r8) :: dxdr,dxdz,dydr,dydz
      INTEGER(i4) :: imode,dm

      REAL(r8), DIMENSION(laq%nqty) :: dfqx,dfqy,d2fqxx,d2fqxy,   &
     &                                 d2fqyy

      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rhs,sol_real,sol_imag
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: mat
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: ff,dffx,dffy,d2ffxx,  &
     &                                            d2ffxy,d2ffyy
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: derivs,dderivs

      REAL(r8), DIMENSION(laq%nqty) :: eqfld
      REAL(r8), DIMENSION(laq%nqty,3) :: eqderivs
!-----------------------------------------------------------------------
!     Check for consistent inputs
!-----------------------------------------------------------------------
      IF (SIZE(coord)<2)                                                &
     &  CALL nim_stop('get_field: coord must contain x, y')
      IF (SIZE(outfld,1)/=laq%nqty)                                     &
     &  CALL nim_stop('get_field: output nqty not match laq nqty')
      IF (SIZE(outfld,2)/=laq%nfour+1)                                  &
     &  CALL nim_stop('get_field: output nfour not match laq nfour')
      IF(PRESENT(laqeq))THEN
        IF(laq%nqty/=laqeq%nqty) THEN
          CALL nim_stop('get_field: laq nqty not match laqeq nqty')
        ENDIF
      ENDIF
      IF (laq%nfour/=nmodes_total)                                      &
     &  CALL nim_stop('get_field: processor must have all modes')
      dm=0
      IF (PRESENT(dmode)) dm=dmode
!-----------------------------------------------------------------------
!     Allocate local arrays
!-----------------------------------------------------------------------
      ALLOCATE(ff(laq%nqty,nmodes_total),dffx(laq%nqty,nmodes_total),   &
     &       dffy(laq%nqty,nmodes_total),d2ffxx(laq%nqty,nmodes_total), &
     &     d2ffxy(laq%nqty,nmodes_total),d2ffyy(laq%nqty,nmodes_total))
      ALLOCATE( derivs(laq%nqty,nmodes_total,3))
      ALLOCATE(dderivs(laq%nqty,nmodes_total,6))
      ALLOCATE(mat(laq%nqty,nmodes_total,3,3),                          &
     &         rhs(laq%nqty,nmodes_total,3),                            &
     &    sol_real(laq%nqty,nmodes_total,3),                            &
     &    sol_imag(laq%nqty,nmodes_total,3))
!-----------------------------------------------------------------------
!     Get coordinates and metric
!-----------------------------------------------------------------------
      x=coord(1)
      y=coord(2)
      CALL lagr_quad_eval(rb_cel(1)%rz,x,y,2_i4,rzf,drzdx,drzdy,        &
     &                    d2rzdxx,d2rzdxy,d2rzdyy)
      IF (geom=='tor') THEN
        bigr=rzf(1)
      ELSE
        bigr=1._r8
      END IF
      IF (dm>0) THEN
        CALL calcderivs (drzdx,drzdy,dxdr,dxdz,dydr,dydz)
        IF (dm>1) THEN
          CALL calc2derivs(drzdx,drzdy,mat)
          CALL math_solve_3x3(mat,sol_real,rhs,'factor')
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     Find the equilibrium field and it's derivatives
!-----------------------------------------------------------------------
      outfld =(0._r8,0._r8)
      IF (PRESENT(laqeq)) THEN
        CALL lagr_quad_eval(laqeq,x,y,2_i4,eqfld,dfqx,dfqy,       &
     &                      d2fqxx,d2fqxy,d2fqyy)
        outfld(:,0,0)=(1.0_r8,0.0_r8)*eqfld
        IF (dm>0) THEN
          eqderivs(:,1)=(dfqx*dxdr+dfqy*dydr)
          eqderivs(:,2)=(dfqx*dxdz+dfqy*dydz)
          outfld(:,0,1:2)=(1.0_r8,0.0_r8)*eqderivs(:,1:2)
        ENDIF
        IF (dm>1) THEN
          rhs(:,1,1)=d2fqxx-eqderivs(:,1)*d2rzdxx(1)-                 &
     &                      eqderivs(:,2)*d2rzdxx(2)
          rhs(:,1,2)=d2fqyy-eqderivs(:,1)*d2rzdyy(1)-                 &
     &                      eqderivs(:,2)*d2rzdyy(2)
          rhs(:,1,3)=d2fqxy-eqderivs(:,1)*d2rzdxy(1)-                 &
     &                      eqderivs(:,2)*d2rzdxy(2)
          CALL math_solve_3x3(mat,sol_real,rhs,'solve ')
          outfld(:,0,4)=(1._r8,0._r8)*sol_real(:,1,1)
          outfld(:,0,5)=(1._r8,0._r8)*sol_real(:,1,3)
          outfld(:,0,6)=(0._r8,0._r8)
          outfld(:,0,7)=(1._r8,0._r8)*sol_real(:,1,2)
          outfld(:,0,8)=(0._r8,0._r8)
          outfld(:,0,9)=(0._r8,0._r8)
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     Calculate the perturbed fields and their derivatives
!-----------------------------------------------------------------------
      derivs =(0._r8,0._r8)
      dderivs=(0._r8,0._r8)
      CALL lagr_quad_eval(laq,x,y,2_i4,ff,dffx,dffy,                    &
     &                    d2ffxx,d2ffxy,d2ffyy)
!-----------------------------------------------------------------------
!     Set the perturbed fields
!-----------------------------------------------------------------------
      DO imode=1,nmodes_total
        outfld(:,imode,0)=ff(:,imode)
      ENDDO
!-----------------------------------------------------------------------
!     Calculate the first order derivatives
!-----------------------------------------------------------------------
      IF (dm > 0) THEN
        derivs(:,:,1)=dffx*dxdr+dffy*dydr
        derivs(:,:,2)=dffx*dxdz+dffy*dydz
        DO imode=1,nmodes_total
          derivs(:,imode,3)=(0,1)*keff_total(imode)*ff(:,imode)/bigr
          outfld(:,imode,1)=derivs(:,imode,1)
          outfld(:,imode,2)=derivs(:,imode,2)
          outfld(:,imode,3)=derivs(:,imode,3)
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     Calculate the second order derivatives
!-----------------------------------------------------------------------
      IF (dm > 1) THEN
        rhs=0._r8
        rhs(:,:,1)=REAL(d2ffxx-derivs(:,:,1)*d2rzdxx(1)-              &
     &                           derivs(:,:,2)*d2rzdxx(2))
        rhs(:,:,2)=REAL(d2ffyy-derivs(:,:,1)*d2rzdyy(1)-              &
     &                           derivs(:,:,2)*d2rzdyy(2))
        rhs(:,:,3)=REAL(d2ffxy-derivs(:,:,1)*d2rzdxy(1)-              &
     &                           derivs(:,:,2)*d2rzdxy(2))
        CALL math_solve_3x3(mat,sol_real,rhs,'solve ')
        rhs=0._r8
        rhs(:,:,1)=AIMAG(d2ffxx-derivs(:,:,1)*d2rzdxx(1)-             &
     &                            derivs(:,:,2)*d2rzdxx(2))
        rhs(:,:,2)=AIMAG(d2ffyy-derivs(:,:,1)*d2rzdyy(1)-             &
     &                            derivs(:,:,2)*d2rzdyy(2))
        rhs(:,:,3)=AIMAG(d2ffxy-derivs(:,:,1)*d2rzdxy(1)-             &
     &                            derivs(:,:,2)*d2rzdxy(2))
        CALL math_solve_3x3(mat,sol_imag,rhs,'solve ')

        dderivs(:,:,1)=dderivs(:,:,1)+sol_real(:,:,1)+                &
     &                            (0,1)*sol_imag(:,:,1)
        dderivs(:,:,2)=dderivs(:,:,2)+sol_real(:,:,3)+                &
     &                            (0,1)*sol_imag(:,:,3)
        dderivs(:,:,4)=dderivs(:,:,4)+sol_real(:,:,2)+                &
     &                            (0,1)*sol_imag(:,:,2)
        DO imode=1,nmodes_total
          dderivs(:,imode,3)=dderivs(:,imode,3)+                      &
     &                         (0,1)*keff_total(imode)*derivs(:,imode,1)/bigr
          dderivs(:,imode,5)=dderivs(:,imode,5)+                      &
     &                         (0,1)*keff_total(imode)*derivs(:,imode,2)/bigr
          dderivs(:,imode,6)=dderivs(:,imode,6)-                      &
     &                         ff(:,imode)*(keff_total(imode)/bigr)**2
          outfld(:,imode,4)=dderivs(:,imode,1)
          outfld(:,imode,5)=dderivs(:,imode,2)
          outfld(:,imode,6)=dderivs(:,imode,3)
          outfld(:,imode,7)=dderivs(:,imode,4)
          outfld(:,imode,8)=dderivs(:,imode,5)
          outfld(:,imode,9)=dderivs(:,imode,6)
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     Clean up and exit
!-----------------------------------------------------------------------
      DEALLOCATE(ff,dffx,dffy,d2ffxx,d2ffxy,d2ffyy)
      DEALLOCATE(derivs,dderivs,mat,rhs,sol_real,sol_imag)

      RETURN
      CONTAINS
        SUBROUTINE calcderivs(drzdx,drzdy,dxdr,dxdz,dydr,dydz)
        REAL(r8), DIMENSION(2), INTENT(IN) :: drzdx,drzdy
        REAL(r8), INTENT(OUT) :: dxdr,dxdz,dydr,dydz
        REAL(r8) :: jac
        ! find the jacobian, then perform a matrix inversion.
        ! This could be replaced with a call to math_inv_2x2
        jac = drzdx(1)*drzdy(2)-drzdx(2)*drzdy(1)
        dxdr= drzdy(2)/jac; dydr=-drzdx(2)/jac
        dxdz=-drzdy(1)/jac; dydz= drzdx(1)/jac
        RETURN
        END SUBROUTINE calcderivs

        SUBROUTINE calc2derivs(drzdx,drzdy,mat)
        REAL(r8), DIMENSION(2), INTENT(IN) :: drzdx,drzdy
        REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: mat
        ! find the jacobian matrix for 2nd derivatives
        ! sets up the C matrix: see get_comp_field_dm documentation
        mat(:,:,1,1)=drzdx(1)**2
        mat(:,:,1,2)=drzdx(2)**2
        mat(:,:,1,3)=2._r8*drzdx(1)*drzdx(2)
        mat(:,:,2,1)=drzdy(1)**2
        mat(:,:,2,2)=drzdy(2)**2
        mat(:,:,2,3)=2._r8*drzdy(1)*drzdy(2)
        mat(:,:,3,1)=drzdx(1)*drzdy(1)
        mat(:,:,3,2)=drzdx(2)*drzdy(2)
        mat(:,:,3,3)=drzdx(1)*drzdy(2)+drzdx(2)*drzdy(1)
        RETURN
        END SUBROUTINE calc2derivs

      END SUBROUTINE get_comp_field_dm
!-----------------------------------------------------------------------
      END MODULE get_field_mod
