module interp_gpu

use params_gpu
use pspline_gpu

IMPLICIT NONE
  
#ifdef PSPLINE

TYPE :: KORC_2D_FIELDS_INTERPOLANT
    !! @note Derived type containing 2-D PSPLINE interpolants for
    !! cylindrical components of vector fields \(\mathbf{F}(R,Z) =
    !! F_R\hat{e}_R + F_\phi\hat{e}_phi+ F_Z\hat{e}_Z\).
    !! Real precision of 8 bytes. @endnote
    TYPE(EZspline2)    :: X
    !! Interpolant of \(F_X(R,Z)\).
    TYPE(EZspline2)    :: Y
    !! Interpolant of \(F_Y(R,Z)\).
    TYPE(EZspline2)    :: Z
    !! Interpolant of \(F_Z(R,Z)\).
    INTEGER               :: NX
    !! Size of mesh containing the field data along the \(R\)-axis.
    INTEGER               :: NY
    !! Size of mesh containing the field data along the \(Z\)-axis.
    INTEGER, DIMENSION(2) :: BCSR = (/ 0, 0 /)
    !! Not-a-knot boundary condition for the interpolants at both
    !! ends of the \(R\) direction.
    INTEGER, DIMENSION(2) :: BCSZ = (/ 0, 0 /)
    !! Not-a-knot boundary condition for the interpolants at both
    !! ends of the \(Z\) direction.
END TYPE KORC_2D_FIELDS_INTERPOLANT 

#endif

CONTAINS

#ifdef PSPLINE

subroutine initialize_interpolants(XF,YF,BF_X,BF_Y,BF_Z,EF_X,EF_Y,EF_Z,bfield_2d,efield_2d)
    REAL(rp),DIMENSION(20), INTENT(IN) :: XF,YF
    REAL(rp),DIMENSION(20,20), INTENT(IN) :: BF_X,BF_Y,BF_Z,EF_X,EF_Y,EF_Z
    TYPE(KORC_2D_FIELDS_INTERPOLANT),INTENT(INOUT) :: bfield_2d,efield_2d
    INTEGER                                        :: ezerr
  
    write(output_write,'("* * * * INITIALIZING FIELDS INTERPOLANT * * * *")')

    bfield_2d%NX = size(XF)
    bfield_2d%NY = size(YF)

    efield_2d%NX = size(XF)
    efield_2d%NY = size(YF)

    call EZspline_init2(bfield_2d%X,bfield_2d%NX,bfield_2d%NY, &
        bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
    call EZspline_error(ezerr)

    call EZspline_init2(bfield_2d%Y,bfield_2d%NX,bfield_2d%NY, &
        bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
    call EZspline_error(ezerr)

    call EZspline_init2(bfield_2d%Z,bfield_2d%NX,bfield_2d%NY, &
    bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
    call EZspline_error(ezerr)

    call EZspline_init2(efield_2d%X,efield_2d%NX,efield_2d%NY, &
    efield_2d%BCSR,efield_2d%BCSZ,ezerr)
    call EZspline_error(ezerr)

    call EZspline_init2(efield_2d%Y,efield_2d%NX,efield_2d%NY, &
    efield_2d%BCSR,efield_2d%BCSZ,ezerr)
    call EZspline_error(ezerr)

    call EZspline_init2(efield_2d%Z,efield_2d%NX,efield_2d%NY, &
    efield_2d%BCSR,efield_2d%BCSZ,ezerr)
    call EZspline_error(ezerr)

    bfield_2d%X%x1 = XF
    bfield_2d%X%x2 = YF

    bfield_2d%Y%x1 = XF
    bfield_2d%Y%x2 = YF

    bfield_2d%Z%x1 = XF
    bfield_2d%Z%x2 = YF

    efield_2d%X%x1 = XF
    efield_2d%X%x2 = YF

    efield_2d%Y%x1 = XF
    efield_2d%Y%x2 = YF

    efield_2d%Z%x1 = XF
    efield_2d%Z%x2 = YF

    call EZspline_setup2(bfield_2d%X,BF_X,ezerr,.TRUE.)
    call EZspline_error(ezerr)

    call EZspline_setup2(bfield_2d%Y,BF_Y,ezerr,.TRUE.)
    call EZspline_error(ezerr)

    call EZspline_setup2(bfield_2d%Z,BF_Z,ezerr,.TRUE.)
    call EZspline_error(ezerr)

    call EZspline_setup2(efield_2d%X,EF_X,ezerr,.TRUE.)
    call EZspline_error(ezerr)

    call EZspline_setup2(efield_2d%Y,EF_Y,ezerr,.TRUE.)
    call EZspline_error(ezerr)

    call EZspline_setup2(efield_2d%Z,EF_Z,ezerr,.TRUE.)
    call EZspline_error(ezerr)

end subroutine initialize_interpolants
  
subroutine interp_fields(XX,YY,BX,BY,BZ,EX,EY,EZ,bfield_2d,efield_2d)
    !$acc routine seq
    REAL(rp),INTENT(IN)   :: XX,YY
    REAL(rp),INTENT(OUT)   :: BX,BY,BZ,EX,EY,EZ
    TYPE(KORC_2D_FIELDS_INTERPOLANT),INTENT(IN) :: bfield_2d,efield_2d
    INTEGER                                        :: ezerr

    !$acc routine (EZspline_interp2_FOvars_cloud) seq
    !$acc routine (EZspline_error) seq
  
    call EZspline_interp2_FOvars_cloud(bfield_2d%X,bfield_2d%Y,bfield_2d%Z, &
        efield_2d%X,efield_2d%Y,efield_2d%Z,XX,YY,BX,BY,BZ,EX,EY,EZ,ezerr)
    call EZspline_error(ezerr)
  
end subroutine
  
subroutine finalize_interpolants(bfield_2d,efield_2d)  
   TYPE(KORC_2D_FIELDS_INTERPOLANT),INTENT(IN) :: bfield_2d,efield_2d
   INTEGER                                        :: ezerr

    call Ezspline_free2(bfield_2d%X, ezerr)
    call Ezspline_free2(bfield_2d%Y, ezerr)
    call Ezspline_free2(bfield_2d%Z, ezerr)

    call Ezspline_free2(efield_2d%X, ezerr)
    call Ezspline_free2(efield_2d%Y, ezerr)
    call Ezspline_free2(efield_2d%Z, ezerr)
  
end subroutine finalize_interpolants

#endif PSPLINE
  
end module interp_gpu