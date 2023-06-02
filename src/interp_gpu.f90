module interp_gpu

use params_gpu

#ifdef PSPLINE
use EZspline_obj	! psplines module
use EZspline		! psplines module
#endif

IMPLICIT NONE
  
#ifdef PSPLINE
TYPE, PRIVATE :: KORC_2D_FIELDS_INTERPOLANT
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
TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: bfield_2d
!! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
!! the magnetic field.
TYPE(KORC_2D_FIELDS_INTERPOLANT), PRIVATE      :: efield_2d
!! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
!! the electric field.
INTEGER                                        :: ezerr
!! Error status during PSPLINE interpolations.

PUBLIC :: interp_fields,initialize_interpolants,finalize_interpolants
    
CONTAINS

subroutine initialize_interpolants(XF,YF,BF_X,BF_Y,BF_Z,EF_X,EF_Y,EF_Z)
    REAL(rp),DIMENSION(20), INTENT(IN) :: XF,YF
    REAL(rp),DIMENSION(20,20), INTENT(IN) :: BF_X,BF_Y,BF_Z,EF_X,EF_Y,EF_Z
  
    write(output_write,'("* * * * INITIALIZING FIELDS INTERPOLANT * * * *")')

    bfield_2d%NX = size(XF)
    bfield_2d%NY = size(YF)

    efield_2d%NX = size(XF)
    efield_2d%NY = size(YF)

    call EZspline_init(bfield_2d%X,bfield_2d%NX,bfield_2d%NY, &
        bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
    call EZspline_error(ezerr)

    call EZspline_init(bfield_2d%Y,bfield_2d%NX,bfield_2d%NY, &
        bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
    call EZspline_error(ezerr)

    call EZspline_init(bfield_2d%Z,bfield_2d%NX,bfield_2d%NY, &
    bfield_2d%BCSR,bfield_2d%BCSZ,ezerr)
    call EZspline_error(ezerr)

    call EZspline_init(efield_2d%X,efield_2d%NX,efield_2d%NY, &
    efield_2d%BCSR,efield_2d%BCSZ,ezerr)
    call EZspline_error(ezerr)

    call EZspline_init(efield_2d%Y,efield_2d%NX,efield_2d%NY, &
    efield_2d%BCSR,efield_2d%BCSZ,ezerr)
    call EZspline_error(ezerr)

    call EZspline_init(efield_2d%Z,efield_2d%NX,efield_2d%NY, &
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
  
    call EZspline_setup(bfield_2d%X,BF_X,ezerr,.TRUE.)
    call EZspline_error(ezerr)

    call EZspline_setup(bfield_2d%Y,BF_Y,ezerr,.TRUE.)
    call EZspline_error(ezerr)

    call EZspline_setup(bfield_2d%Z,BF_Z,ezerr,.TRUE.)
    call EZspline_error(ezerr)

    call EZspline_setup(efield_2d%X,EF_X,ezerr,.TRUE.)
    call EZspline_error(ezerr)

    call EZspline_setup(efield_2d%Y,EF_Y,ezerr,.TRUE.)
    call EZspline_error(ezerr)

    call EZspline_setup(efield_2d%Z,EF_Z,ezerr,.TRUE.)
    call EZspline_error(ezerr)
  
    write(output_write,'("* * * * * * INTERPOLANT INITIALIZED * * * * * *",/)')

end subroutine initialize_interpolants 
  
subroutine interp_fields(XX,YY,BX,BY,BZ,EX,EY,EZ)
    REAL(rp),INTENT(IN)   :: XX,YY
    REAL(rp),DIMENSION(1),INTENT(OUT)   :: BX,BY,BZ
    REAL(rp),DIMENSION(1),INTENT(OUT)   :: EX,EY,EZ
  
    call EZspline_interp(bfield_2d%X,bfield_2d%Y,bfield_2d%Z, &
        efield_2d%X,efield_2d%Y,efield_2d%Z,1,(/XX/),(/YY/),BX,BY,BZ,EX,EY,EZ,ezerr)
    call EZspline_error(ezerr)
  
end subroutine
  
subroutine finalize_interpolants  

    call Ezspline_free(bfield_2d%X, ezerr)
    call Ezspline_free(bfield_2d%Y, ezerr)
    call Ezspline_free(bfield_2d%Z, ezerr)

    call Ezspline_free(efield_2d%X, ezerr)
    call Ezspline_free(efield_2d%Y, ezerr)
    call Ezspline_free(efield_2d%Z, ezerr)
  
end subroutine finalize_interpolants
#endif
  
end module interp_gpu