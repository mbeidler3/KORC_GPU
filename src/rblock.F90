!-----------------------------------------------------------------------
!     $Id: rblock.F90 6435 2019-02-21 17:53:45Z charlson $
!     subprograms for handling finite element computations on logically
!     rectangular grid-blocks.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.
!-----------------------------------------------------------------------
!     1.  rblock_set.
!     2.  rblock_make_real_matrix.
!     3.  rblock_make_comp_matrix.
!     4.  rblock_get_real_rhs.
!     5.  rblock_get_comp_rhs.
!     6.  rblock_get_comp_rhs_q.
!     7.  rblock_basis_set.
!     9.  rblock_bicube_set.
!     10. rblock_real_qp_update.
!     11. rblock_comp_qp_update.
!     12. rblock_real_qp_alloc.
!     13. rblock_comp_qp_alloc.
!     14. rblock_real_qp_dealloc.
!     15. rblock_comp_qp_dealloc.
!     17. rblock_real_qpe_update.
!     18. rblock_comp_qpe_update.
!     19. rblock_real_qpe_alloc.
!     20. rblock_comp_qpe_alloc.
!     21. rblock_real_qpe_dealloc.
!     22. rblock_comp_qpe_dealloc.
!-----------------------------------------------------------------------
!     declarations.
!-----------------------------------------------------------------------
      MODULE rblock
      USE local
      USE rblock_type_mod
      IMPLICIT NONE

      INTERFACE rblock_make_matrix
        MODULE PROCEDURE rblock_make_real_matrix,rblock_make_comp_matrix
      END INTERFACE

      INTERFACE rblock_get_rhs
        MODULE PROCEDURE rblock_get_real_rhs,rblock_get_comp_rhs
      END INTERFACE

      INTERFACE rblock_qp_update
        MODULE PROCEDURE rblock_real_qp_update,rblock_comp_qp_update
      END INTERFACE

      INTERFACE rblock_qp_alloc
        MODULE PROCEDURE rblock_real_qp_alloc,rblock_comp_qp_alloc
      END INTERFACE

      INTERFACE rblock_qp_dealloc
        MODULE PROCEDURE rblock_real_qp_dealloc,rblock_comp_qp_dealloc
      END INTERFACE

      INTERFACE rblock_qp_fft_save
        MODULE PROCEDURE qp_fft_save,qp_fft_save_zflr,  &
          &         qp_fft_noeq_save,qp_fft_noeq_save_zflr 
      END INTERFACE

      INTERFACE rblock_qpe_update
        MODULE PROCEDURE rblock_real_qpe_update,rblock_comp_qpe_update
      END INTERFACE

      INTERFACE rblock_qpe_alloc
        MODULE PROCEDURE rblock_real_qpe_alloc,rblock_comp_qpe_alloc
      END INTERFACE

      INTERFACE rblock_qpe_dealloc
        MODULE PROCEDURE rblock_real_qpe_dealloc,rblock_comp_qpe_dealloc
      END INTERFACE

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 1. rblock_set.
!     set the locations and weights for quadratures.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_set(ngr,poly_degree,int_formula,surf_int_f,rb)

      INTEGER(i4), INTENT(IN) :: ngr,poly_degree
      CHARACTER(8), INTENT(IN) :: int_formula,surf_int_f
      TYPE(rblock_type), INTENT(INOUT) :: rb

      REAL(r8), DIMENSION(ngr+poly_degree-1) :: xg1d,wg1d
      INTEGER(i4) :: ix,iy,ig
!-----------------------------------------------------------------------
!     set number of quadrature points and weights according to input.
!     the number of points is now adjusted automatically with
!     poly_degree.
!-----------------------------------------------------------------------
      IF (int_formula=='gaussian'.OR.rb%degenerate.OR.rb%r0block) THEN
        CALL gauleg(0._r8,1._r8,xg1d,wg1d,ngr+poly_degree-1_i4)
      ELSE
        CALL lobleg(0._r8,1._r8,xg1d,wg1d,ngr+poly_degree-1_i4)
      ENDIF
      rb%ng=(ngr+poly_degree-1)**2
      rb%nge=ngr+poly_degree-1
      ALLOCATE(rb%xg(rb%ng),rb%yg(rb%ng),rb%wg(rb%ng),rb%xge(rb%nge),   &
     &         rb%wge(rb%nge))
      ig=0
      DO iy=1,rb%nge
        DO ix=1,rb%nge
          ig=ig+1
          rb%xg(ig)=xg1d(ix)
          rb%yg(ig)=xg1d(iy)
          rb%wg(ig)=wg1d(ix)*wg1d(iy)
        ENDDO
      ENDDO
      IF (surf_int_f=='lobatto') THEN
        CALL lobleg(0._r8,1._r8,xg1d,wg1d,ngr+poly_degree-1_i4)
      ELSE
        CALL gauleg(0._r8,1._r8,xg1d,wg1d,ngr+poly_degree-1_i4)
      ENDIF
      DO ix=1,rb%nge
        rb%xge(ix)=xg1d(ix)
        ! TODO:: think about this.
        rb%wge(ix)=wg1d(ix) !*wg1d(ngr+poly_degree-1)
      ENDDO
!-----------------------------------------------------------------------
!     terminate rblock_set.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_set
!-----------------------------------------------------------------------
!     subprogram 2. rblock_make_real_matrix.
!     computes a linear response matrix for a supplied integrand
!     subprogram that fits the interface block at the beginning
!     of this subprogram.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_make_real_matrix(rb,mat_str,integrand,nqty)
      USE tblock_type_mod
      USE matrix_type_mod
      USE pardata, ONLY: global2local
      USE global, ONLY: mpsq_block

      TYPE(rblock_type), INTENT(INOUT) :: rb
      TYPE(rbl_mat_type), INTENT(INOUT) :: mat_str
      INTEGER(i4), INTENT(IN) :: nqty

      INTEGER(i4) :: ityp,jtyp,idof,jdof,idofst,jdofst,                 &
     &               jxoff,jyoff,ipol,ix,iy,ix0,iy0,mx,my,mpseudo
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER :: mat

      TYPE (tblock_type) :: tdum

      REAL(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: integr
!-----------------------------------------------------------------------
!     interface block for the integrand computation.
!-----------------------------------------------------------------------
#include "integrand_real.finc"
!-----------------------------------------------------------------------
!     preliminary computations.
!-----------------------------------------------------------------------
      IF(ALLOCATED(mpsq_block)) THEN
        mpseudo=mpsq_block(global2local(rb%id))
      ELSE
        mpseudo=1_i4
      ENDIF
      mx=rb%mx
      my=rb%my
!-----------------------------------------------------------------------
!     initialize matrix arrays.
!-----------------------------------------------------------------------
      DO jtyp=1,mat_str%nbtype
        DO ityp=1,mat_str%nbtype
          mat_str%mat(ityp,jtyp)%arr=0._r8
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     flag the tblock as a dummy.
!-----------------------------------------------------------------------
      tdum%tgeom%mvert=-1
!-----------------------------------------------------------------------
!     quadrature-point looping is now done inside the integr routine,
!     and integr is summed for each basis function, element by
!     element.
!
!     factors of Jacobian and quadrature weight are already in the basis
!     function arrays.
!-----------------------------------------------------------------------
      ALLOCATE(integr(nqty,nqty,rb%mx*rb%my,mat_str%nbasis_el,       &
     &                   mat_str%nbasis_el))
      CALL integrand(integr,rb%bigr,rb,tdum,mpseudo)
!-----------------------------------------------------------------------
!     assemble and accumulate the contributions from each element.
!     the degree-of-freedom arrays are now used to replace
!     complicated looping with pre-computed information.
!-----------------------------------------------------------------------
      idofst=1
      DO ityp=1,mat_str%nbtype
        IF (mat_str%nq_type(ityp)==0) CYCLE
        jdofst=1
        DO jtyp=1,mat_str%nbtype
          IF (mat_str%nq_type(jtyp)==0) CYCLE
          mat=>mat_str%mat(jtyp,ityp)%arr
          DO ipol=1,mx*my
            iy0=(ipol-1)/mx+1
            ix0=ipol-mx*(iy0-1)
            DO idof=idofst,mat_str%den_type(ityp)
              iy=iy0+mat_str%dof_iy(idof)
              ix=ix0+mat_str%dof_ix(idof)
              DO jdof=jdofst,mat_str%den_type(jtyp)
                jxoff=mat_str%dof_ix(jdof)-mat_str%dof_ix(idof)
                jyoff=mat_str%dof_iy(jdof)-mat_str%dof_iy(idof)
                mat(mat_str%dof_iq(jdof),jxoff,jyoff,                   &
     &              mat_str%dof_iq(idof),ix,iy)=                        &
     &            mat(mat_str%dof_iq(jdof),jxoff,jyoff,                 &
     &                mat_str%dof_iq(idof),ix,iy)+                      &
     &            integr(mat_str%dof_iv(jdof),mat_str%dof_iv(idof),     &
     &                 ipol,mat_str%dof_ib(jdof),mat_str%dof_ib(idof))
              ENDDO
            ENDDO
          ENDDO
          jdofst=mat_str%den_type(jtyp)+1
        ENDDO
        idofst=mat_str%den_type(ityp)+1
      ENDDO

      DEALLOCATE(integr)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_make_real_matrix
!-----------------------------------------------------------------------
!     subprogram 3. rblock_make_comp_matrix.
!     computes a linear response matrix for a supplied integrand
!     subprogram that fits the interface block at the beginning
!     of this subprogram.
!
!     complex version.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_make_comp_matrix(rb,mat_str,integrand,nqty)
      USE tblock_type_mod
      USE matrix_type_mod
      USE pardata, ONLY: global2local
      USE global, ONLY: mpsq_block

      TYPE(rblock_type), INTENT(INOUT), TARGET :: rb
      TYPE(rbl_comp_mat_type), INTENT(INOUT) :: mat_str
      INTEGER(i4), INTENT(IN) :: nqty

      INTEGER(i4) :: ityp,jtyp,idof,jdof,idofst,jdofst,                 &
     &               jxoff,jyoff,ipol,ix,iy,ix0,iy0,mx,my,mpseudo
      COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER :: mat

      TYPE (tblock_type) :: tdum

      COMPLEX(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: integr
!-----------------------------------------------------------------------
!     interface block for the integrand computation.
!-----------------------------------------------------------------------
#include "integrand_comp.finc"
!-----------------------------------------------------------------------
!     preliminary computations.
!-----------------------------------------------------------------------
      IF(ALLOCATED(mpsq_block)) THEN
        mpseudo=mpsq_block(global2local(rb%id))
      ELSE
        mpseudo=1_i4
      ENDIF
      mx=rb%mx
      my=rb%my
!-----------------------------------------------------------------------
!     initialize matrix arrays.
!-----------------------------------------------------------------------
      DO jtyp=1,mat_str%nbtype
        DO ityp=1,mat_str%nbtype
          mat_str%mat(ityp,jtyp)%arr=0._r8
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     flag the tblock as a dummy.
!-----------------------------------------------------------------------
      tdum%tgeom%mvert=-1
!-----------------------------------------------------------------------
!     quadrature-point looping is now done inside the integr routine,
!     and integr is summed for each basis function, element by
!     element.
!
!     factors of Jacobian and quadrature weight are already in the basis
!     function arrays.
!-----------------------------------------------------------------------
      ALLOCATE(integr(nqty,nqty,rb%mx*rb%my,mat_str%nbasis_el,          &
     &                   mat_str%nbasis_el))
      CALL integrand(integr,rb%bigr,rb,tdum,mpseudo)
!-----------------------------------------------------------------------
!     assemble and accumulate the contributions from each element.
!     the degree-of-freedom arrays are now used to replace
!     complicated looping with pre-computed information.
!-----------------------------------------------------------------------
      idofst=1
      DO ityp=1,mat_str%nbtype
        IF (mat_str%nq_type(ityp)==0) CYCLE
        jdofst=1
        DO jtyp=1,mat_str%nbtype
          IF (mat_str%nq_type(jtyp)==0) CYCLE
          mat=>mat_str%mat(jtyp,ityp)%arr
          DO ipol=1,mx*my
            iy0=(ipol-1)/mx+1
            ix0=ipol-mx*(iy0-1)
            DO idof=idofst,mat_str%den_type(ityp)
              iy=iy0+mat_str%dof_iy(idof)
              ix=ix0+mat_str%dof_ix(idof)
              DO jdof=jdofst,mat_str%den_type(jtyp)
                jxoff=mat_str%dof_ix(jdof)-mat_str%dof_ix(idof)
                jyoff=mat_str%dof_iy(jdof)-mat_str%dof_iy(idof)
                mat(mat_str%dof_iq(jdof),jxoff,jyoff,                   &
     &              mat_str%dof_iq(idof),ix,iy)=                        &
     &            mat(mat_str%dof_iq(jdof),jxoff,jyoff,                 &
     &                mat_str%dof_iq(idof),ix,iy)+                      &
     &            integr(mat_str%dof_iv(jdof),mat_str%dof_iv(idof),     &
     &                 ipol,mat_str%dof_ib(jdof),mat_str%dof_ib(idof))
              ENDDO
            ENDDO
          ENDDO
          jdofst=mat_str%den_type(jtyp)+1
        ENDDO
        idofst=mat_str%den_type(ityp)+1
      ENDDO

      DEALLOCATE(integr)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_make_comp_matrix
!-----------------------------------------------------------------------
!     subprogram 4. rblock_get_real_rhs.
!     performs finite-element integrations for a rhs of an equation
!     producing real data, where the integrand is computed with a
!     supplied subroutine.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_get_real_rhs(rb,rhs,integrand,nq)
      USE tblock_type_mod
      USE vector_type_mod
      USE pardata, ONLY: global2local
      USE global, ONLY: mpsq_block

      INTEGER(i4), INTENT(IN) :: nq
      TYPE(rblock_type), INTENT(INOUT), TARGET :: rb
      TYPE(vector_type), INTENT(INOUT) :: rhs

      INTEGER(i4) :: ix,iy,ipol,mx,my,mxm1,mym1,iq,ig,iv,nv,ib
      INTEGER(i4) :: start_horz,start_vert,start_int,                   &
     &               n_grid,n_horz,n_vert,n_int,mpseudo

      TYPE (tblock_type) :: tdum

      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: integr
      REAL(r8), DIMENSION(:,:,:), POINTER :: rhsg
      REAL(r8), DIMENSION(:,:,:,:), POINTER :: rhsh,rhsv,rhsi
!-----------------------------------------------------------------------
!     interface block for the integrand computation.
!-----------------------------------------------------------------------
#include "integrand_real_rhs.finc"
!-----------------------------------------------------------------------
!     preliminary computations.
!-----------------------------------------------------------------------
      IF(ALLOCATED(mpsq_block)) THEN
        mpseudo=mpsq_block(global2local(rb%id))
      ELSE
        mpseudo=1_i4
      ENDIF
      mx=rb%mx
      my=rb%my
      mxm1=mx-1
      mym1=my-1
!-----------------------------------------------------------------------
!     examine the vector to determine what bases are used.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN
        n_grid=4
        rhsg=>rhs%arr
        rhsg(1:nq,:,:)=0._r8
      ELSE
        n_grid=0
      ENDIF
      iv=n_grid+1
      IF (ASSOCIATED(rhs%arrh)) THEN
        start_horz=iv
        n_horz=SIZE(rhs%arrh,2)
        rhsh=>rhs%arrh
        rhsh(1:nq,:,:,:)=0._r8
      ELSE
        n_horz=0
      ENDIF
      iv=iv+2*n_horz
      IF (ASSOCIATED(rhs%arrv)) THEN
        start_vert=iv
        n_vert=SIZE(rhs%arrv,2)
        rhsv=>rhs%arrv
        rhsv(1:nq,:,:,:)=0._r8
      ELSE
        n_vert=0
      ENDIF
      iv=iv+2*n_vert
      IF (ASSOCIATED(rhs%arri)) THEN
        start_int=iv
        n_int=SIZE(rhs%arri,2)
        rhsi=>rhs%arri
        rhsi(1:nq,:,:,:)=0._r8
      ELSE
        n_int=0
      ENDIF
      iv=iv+n_int-1
!-----------------------------------------------------------------------
!     flag the tblock as a dummy and allocate int.
!-----------------------------------------------------------------------
      tdum%tgeom%mvert=-1
      ALLOCATE(integr(nq,mx*my,iv))
!-----------------------------------------------------------------------
!       evaluate the integrand all quadrature points.
!-----------------------------------------------------------------------
        CALL integrand(integr,rb%bigr,rb,tdum,mpseudo)
!-----------------------------------------------------------------------
!       assemble and accumulate the contributions from each element
!       into the correct arrays.
!       grid vertex-centered bases first.
!
!       factors of Jacobian and quadrature weight are already in the
!       test function arrays.
!-----------------------------------------------------------------------
        IF (n_grid==4) THEN
          ipol=1
          DO iy=0,mym1
            DO ix=0,mxm1
              rhsg(1:nq,ix,iy)=rhsg(1:nq,ix,iy)                         &
     &          +integr(:,ipol,1)
              rhsg(1:nq,ix+1,iy)=rhsg(1:nq,ix+1,iy)                     &
     &          +integr(:,ipol,2)
              rhsg(1:nq,ix,iy+1)=rhsg(1:nq,ix,iy+1)                     &
     &          +integr(:,ipol,3)
              rhsg(1:nq,ix+1,iy+1)=rhsg(1:nq,ix+1,iy+1)                 &
     &          +integr(:,ipol,4)
              ipol=ipol+1
            ENDDO
          ENDDO
        ENDIF
!-----------------------------------------------------------------------
!       horizontal side-centered bases.
!-----------------------------------------------------------------------
        IF (n_horz>0) THEN
          iv=start_horz
          DO ib=1,n_horz
            ipol=1
            DO iy=0,mym1
              DO ix=1,mx
                rhsh(1:nq,ib,ix,iy)=rhsh(1:nq,ib,ix,iy)                 &
     &            +integr(:,ipol,iv)
                rhsh(1:nq,ib,ix,iy+1)=rhsh(1:nq,ib,ix,iy+1)             &
     &            +integr(:,ipol,iv+1)
                ipol=ipol+1
              ENDDO
            ENDDO
            iv=iv+2
          ENDDO
        ENDIF
!-----------------------------------------------------------------------
!       vertical side-centered bases.
!-----------------------------------------------------------------------
        IF (n_vert>0) THEN
          iv=start_vert
          DO ib=1,n_vert
            ipol=1
            DO iy=1,my
              DO ix=0,mxm1
                rhsv(1:nq,ib,ix,iy)=rhsv(1:nq,ib,ix,iy)                 &
     &            +integr(:,ipol,iv)
                rhsv(1:nq,ib,ix+1,iy)=rhsv(1:nq,ib,ix+1,iy)             &
     &            +integr(:,ipol,iv+1)
                ipol=ipol+1
              ENDDO
            ENDDO
            iv=iv+2
          ENDDO
        ENDIF
!-----------------------------------------------------------------------
!       element interior-centered bases.
!-----------------------------------------------------------------------
        IF (n_int>0) THEN
          iv=start_int
          DO ib=1,n_int
            ipol=1
            DO iy=1,my
              DO ix=1,mx
                rhsi(1:nq,ib,ix,iy)=rhsi(1:nq,ib,ix,iy)                 &
     &            +integr(:,ipol,iv)
                ipol=ipol+1
              ENDDO
            ENDDO
            iv=iv+1
          ENDDO
        ENDIF

      DEALLOCATE(integr)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_get_real_rhs
!-----------------------------------------------------------------------
!     subprogram 5. rblock_get_comp_rhs.
!     performs finite-element integrations for a rhs of an equation
!     producing complex data, where the integrand is computed with a
!     supplied subroutine.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_get_comp_rhs(rb,rhs,integrand,nq,nfour)
      USE tblock_type_mod
      USE vector_type_mod
      USE pardata, ONLY: global2local
      USE global, ONLY: mpsq_block

      INTEGER(i4), INTENT(IN) :: nq,nfour
      TYPE(rblock_type), INTENT(INOUT), TARGET :: rb
      TYPE(cvector_type), INTENT(INOUT) :: rhs

      INTEGER(i4) :: ix,iy,mx,my,mxm1,mym1,iq,jf,ig,iv,nv,ib,ipol
      INTEGER(i4) :: start_horz,start_vert,start_int,                   &
     &               n_grid,n_horz,n_vert,n_int,mpseudo

      TYPE (tblock_type) :: tdum

      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: integr
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: rhsg
      COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: rhsh,rhsv,rhsi
!-----------------------------------------------------------------------
!     interface block for the integrand computation.
!-----------------------------------------------------------------------
#include "integrand_comp_rhs.finc"
!-----------------------------------------------------------------------
!     preliminary computations.
!-----------------------------------------------------------------------
      IF(ALLOCATED(mpsq_block)) THEN
        mpseudo=mpsq_block(global2local(rb%id))
      ELSE
        mpseudo=1_i4
      ENDIF
      mx=rb%mx
      my=rb%my
      mxm1=mx-1
      mym1=my-1
!-----------------------------------------------------------------------
!     examine the vector to determine what bases are used.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN
        n_grid=4
        rhsg=>rhs%arr
        rhsg(1:nq,:,:,1:nfour)=0._r8
      ELSE
        n_grid=0
      ENDIF
      iv=n_grid+1
      IF (ASSOCIATED(rhs%arrh)) THEN
        start_horz=iv
        n_horz=SIZE(rhs%arrh,2)
        rhsh=>rhs%arrh
        rhsh(1:nq,:,:,:,1:nfour)=0._r8
      ELSE
        n_horz=0
      ENDIF
      iv=iv+2*n_horz
      IF (ASSOCIATED(rhs%arrv)) THEN
        start_vert=iv
        n_vert=SIZE(rhs%arrv,2)
        rhsv=>rhs%arrv
        rhsv(1:nq,:,:,:,1:nfour)=0._r8
      ELSE
        n_vert=0
      ENDIF
      iv=iv+2*n_vert
      IF (ASSOCIATED(rhs%arri)) THEN
        start_int=iv
        n_int=SIZE(rhs%arri,2)
        rhsi=>rhs%arri
        rhsi(1:nq,:,:,:,1:nfour)=0._r8
      ELSE
        n_int=0
      ENDIF
      iv=iv+n_int-1
!-----------------------------------------------------------------------
!     flag the tblock as a dummy and allocate int.
!-----------------------------------------------------------------------
      tdum%tgeom%mvert=-1
      ALLOCATE(integr(nq,mx*my,iv,nfour))
!-----------------------------------------------------------------------
!       evaluate the integrand all quadrature points.
!-----------------------------------------------------------------------
        CALL integrand(integr,rb%bigr,rb,tdum,mpseudo)
!-----------------------------------------------------------------------
!       assemble and accumulate the contributions from each element
!       into the correct arrays.
!       grid vertex-centered bases first.
!
!       factors of Jacobian and quadrature weight are already in the
!       test function arrays.
!-----------------------------------------------------------------------
        IF (n_grid==4) THEN
          DO jf=1,nfour
            ipol=1
            DO iy=0,mym1
              DO ix=0,mxm1
                rhsg(1:nq,ix,iy,jf)=rhsg(1:nq,ix,iy,jf)                 &
     &            +integr(:,ipol,1,jf)
                rhsg(1:nq,ix+1,iy,jf)=rhsg(1:nq,ix+1,iy,jf)             &
     &            +integr(:,ipol,2,jf)
                rhsg(1:nq,ix,iy+1,jf)=rhsg(1:nq,ix,iy+1,jf)             &
     &            +integr(:,ipol,3,jf)
                rhsg(1:nq,ix+1,iy+1,jf)=rhsg(1:nq,ix+1,iy+1,jf)         &
     &            +integr(:,ipol,4,jf)
                ipol=ipol+1
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!-----------------------------------------------------------------------
!       horizontal side-centered bases.
!-----------------------------------------------------------------------
        IF (n_horz>0) THEN
          DO jf=1,nfour
            iv=start_horz
            DO ib=1,n_horz
              ipol=1
              DO iy=0,mym1
                DO ix=1,mx
                  rhsh(1:nq,ib,ix,iy,jf)=rhsh(1:nq,ib,ix,iy,jf)         &
     &              +integr(:,ipol,iv,jf)
                  rhsh(1:nq,ib,ix,iy+1,jf)=rhsh(1:nq,ib,ix,iy+1,jf)     &
     &              +integr(:,ipol,iv+1,jf)
                  ipol=ipol+1
                ENDDO
              ENDDO
              iv=iv+2
            ENDDO
          ENDDO
        ENDIF
!-----------------------------------------------------------------------
!       vertical side-centered bases.
!-----------------------------------------------------------------------
        IF (n_vert>0) THEN
          DO jf=1,nfour
            iv=start_vert
            DO ib=1,n_vert
              ipol=1
              DO iy=1,my
                DO ix=0,mxm1
                  rhsv(1:nq,ib,ix,iy,jf)=rhsv(1:nq,ib,ix,iy,jf)         &
     &              +integr(:,ipol,iv,jf)
                  rhsv(1:nq,ib,ix+1,iy,jf)=rhsv(1:nq,ib,ix+1,iy,jf)     &
     &              +integr(:,ipol,iv+1,jf)
                  ipol=ipol+1
                ENDDO
              ENDDO
              iv=iv+2
            ENDDO
          ENDDO
        ENDIF
!-----------------------------------------------------------------------
!       element interior-centered bases.
!-----------------------------------------------------------------------
        IF (n_int>0) THEN
          DO jf=1,nfour
            iv=start_int
            DO ib=1,n_int
              ipol=1
              DO iy=1,my
                DO ix=1,mx
                  rhsi(1:nq,ib,ix,iy,jf)=rhsi(1:nq,ib,ix,iy,jf)         &
     &              +integr(:,ipol,iv,jf)
                  ipol=ipol+1
                ENDDO
              ENDDO
              iv=iv+1
            ENDDO
          ENDDO
        ENDIF

      DEALLOCATE(integr)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_get_comp_rhs
!-----------------------------------------------------------------------
!     subprogram 6. rblock_get_comp_rhs_q.
!     performs finite-element integrations for a rhs of an equation
!     producing complex data, where the integr is computed with a
!     supplied subroutine.  this is the same as rblock_get_comp_rhs,
!     except that the rhs array is assumed to have the quantity and
!     Fourier component indices dimensioned correctly for this equation
!     for efficiency.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_get_comp_rhs_q(rb,rhs,integrand)
      USE tblock_type_mod
      USE vector_type_mod
      USE pardata, ONLY: global2local
      USE global, ONLY: mpsq_block

      TYPE(rblock_type), INTENT(INOUT) :: rb
      TYPE(cvector_type), INTENT(INOUT) :: rhs

      INTEGER(i4) :: ix,iy,mx,my,mxm1,mym1,iq,jf,ig,iv,nv,ib,ipol
      INTEGER(i4) :: start_horz,start_vert,start_int,start_disc,mpseudo,&
     &               n_grid,n_horz,n_vert,n_int,nq,nfour,nqd,n_disc

      TYPE (tblock_type) :: tdum

      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: integr
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: rhsg
      COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: rhsh,rhsv,rhsi,rhsd
!-----------------------------------------------------------------------
!     interface block for the integrand computation.
!-----------------------------------------------------------------------
#include "integrand_comp_rhs.finc"
!-----------------------------------------------------------------------
!     preliminary computations.
!-----------------------------------------------------------------------
      IF(ALLOCATED(mpsq_block)) THEN
        mpseudo=mpsq_block(global2local(rb%id))
      ELSE
        mpseudo=1_i4
      ENDIF
      mx=rb%mx
      my=rb%my
      mxm1=mx-1
      mym1=my-1
!-----------------------------------------------------------------------
!     examine the vector to determine what bases are used.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN
        n_grid=4
        rhsg=>rhs%arr
        nq=SIZE(rhsg,1)
        nfour=SIZE(rhsg,4)
        rhsg=0._r8
      ELSE
        n_grid=0
      ENDIF
      iv=n_grid+1
      IF (ASSOCIATED(rhs%arrh)) THEN
        start_horz=iv
        n_horz=SIZE(rhs%arrh,2)
        rhsh=>rhs%arrh
        nq=SIZE(rhsh,1)
        nfour=SIZE(rhsh,5)
        rhsh=0._r8
      ELSE
        n_horz=0
      ENDIF
      iv=iv+2*n_horz
      IF (ASSOCIATED(rhs%arrv)) THEN
        start_vert=iv
        n_vert=SIZE(rhs%arrv,2)
        rhsv=>rhs%arrv
        nq=SIZE(rhsv,1)
        nfour=SIZE(rhsv,5)
        rhsv=0._r8
      ELSE
        n_vert=0
      ENDIF
      iv=iv+2*n_vert
      IF (ASSOCIATED(rhs%arri)) THEN
        start_int=iv
        n_int=SIZE(rhs%arri,2)
        rhsi=>rhs%arri
        nq=SIZE(rhsi,1)
        nfour=SIZE(rhsi,5)
        rhsi=0._r8
      ELSE
        n_int=0
      ENDIF
      iv=iv+n_int
!-----------------------------------------------------------------------
!     a separate set of element-centered computations (for discontinuous
!     auxiliary fields) may be used.  if so, it is assumed that the
!     number of components at each node is no larger than nq.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arrtmp)) THEN
        start_disc=iv
        n_disc=SIZE(rhs%arrtmp,2)
        rhsd=>rhs%arrtmp
        nqd=SIZE(rhsd,1)
        rhsd=0._r8
      ELSE
        n_disc=0
        nqd=0
      ENDIF
      iv=iv+n_disc-1
!-----------------------------------------------------------------------
!     flag the tblock as a dummy and allocate int.
!-----------------------------------------------------------------------
      tdum%tgeom%mvert=-1
      ALLOCATE(integr(nq,mx*my,iv,nfour))
      integr=0
!-----------------------------------------------------------------------
!       evaluate the integr all quadrature points.
!-----------------------------------------------------------------------
      CALL integrand(integr,rb%bigr,rb,tdum,mpseudo)
!-----------------------------------------------------------------------
!       assemble and accumulate the contributions from each element
!       into the correct arrays.
!       grid vertex-centered bases first.
!
!       factors of Jacobian and quadrature weight are already in the
!       test function arrays.
!-----------------------------------------------------------------------
      IF (n_grid==4) THEN
        DO jf=1,nfour
          ipol=1
          DO iy=0,mym1
            DO ix=0,mxm1
              rhsg(:,ix,iy,jf)=rhsg(:,ix,iy,jf)                         &
     &          +integr(:,ipol,1,jf)
              rhsg(:,ix+1,iy,jf)=rhsg(:,ix+1,iy,jf)                     &
     &          +integr(:,ipol,2,jf)
              rhsg(:,ix,iy+1,jf)=rhsg(:,ix,iy+1,jf)                     &
     &          +integr(:,ipol,3,jf)
              rhsg(:,ix+1,iy+1,jf)=rhsg(:,ix+1,iy+1,jf)                 &
     &          +integr(:,ipol,4,jf)
              ipol=ipol+1
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!       horizontal side-centered bases.
!-----------------------------------------------------------------------
      IF (n_horz>0) THEN
        DO jf=1,nfour
          iv=start_horz
          DO ib=1,n_horz
            ipol=1
            DO iy=0,mym1
              DO ix=1,mx
                rhsh(:,ib,ix,iy,jf)=rhsh(:,ib,ix,iy,jf)                 &
     &            +integr(:,ipol,iv,jf)
                rhsh(:,ib,ix,iy+1,jf)=rhsh(:,ib,ix,iy+1,jf)             &
     &            +integr(:,ipol,iv+1,jf)
                ipol=ipol+1
              ENDDO
            ENDDO
            iv=iv+2
          ENDDO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!       vertical side-centered bases.
!-----------------------------------------------------------------------
      IF (n_vert>0) THEN
        DO jf=1,nfour
          iv=start_vert
          DO ib=1,n_vert
            ipol=1
            DO iy=1,my
              DO ix=0,mxm1
                rhsv(:,ib,ix,iy,jf)=rhsv(:,ib,ix,iy,jf)                 &
     &            +integr(:,ipol,iv,jf)
                rhsv(:,ib,ix+1,iy,jf)=rhsv(:,ib,ix+1,iy,jf)             &
     &            +integr(:,ipol,iv+1,jf)
                ipol=ipol+1
              ENDDO
            ENDDO
            iv=iv+2
          ENDDO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!       element interior-centered bases.
!-----------------------------------------------------------------------
      IF (n_int>0) THEN
        DO jf=1,nfour
          iv=start_int
          DO ib=1,n_int
            ipol=1
            DO iy=1,my
              DO ix=1,mx
                rhsi(:,ib,ix,iy,jf)=rhsi(:,ib,ix,iy,jf)                 &
     &            +integr(:,ipol,iv,jf)
                ipol=ipol+1
              ENDDO
            ENDDO
            iv=iv+1
          ENDDO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!       auxiliary computations for discontinuous auxiliary fields.
!-----------------------------------------------------------------------
        IF (n_disc>0) THEN
          DO jf=1,nfour
            iv=start_disc
            DO ib=1,n_disc
              ipol=1
              DO iy=1,my
                DO ix=1,mx
                  rhsd(:,ib,ix,iy,jf)=rhsd(:,ib,ix,iy,jf)               &
     &              +integr(1:nqd,ipol,iv,jf)
                  ipol=ipol+1
                ENDDO
              ENDDO
              iv=iv+1
            ENDDO
          ENDDO
        ENDIF

      DEALLOCATE(integr)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_get_comp_rhs_q
!-----------------------------------------------------------------------
!     subprogram 7. rblock_basis_set.
!     evaluate basis values and derivatives at quadrature points.
!
!     the input array, pd_arr, is a list of basis polynomial degree
!     values needed for continuous (nodal) expansions.
!
!     the input arrays, pdm_arr and pdmmin_arr, list of basis polynomial
!     degree values needed for discontinuous modal expansions.  these
!     expansions may be incomplete, so the minimum degree is provided
!     in the pdmmin_arr list.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_basis_set(rb,pd_arr,pdm_arr,pdmmin_arr,         &
     &                            pdmmax_arr,met_spl,geom)

      USE math_tran
      USE lagr_quad_mod

      TYPE(rblock_type), INTENT(INOUT) :: rb
      INTEGER(i4), DIMENSION(:), INTENT(IN) :: pd_arr,pdm_arr,          &
     &             pdmmin_arr,pdmmax_arr
      CHARACTER(*), INTENT(IN) :: met_spl,geom

      REAL(r8), DIMENSION(:), ALLOCATABLE :: alpha,alphax,alphay
      REAL(r8), DIMENSION(2,rb%mx,rb%my) :: rz,drzdx,drzdy
      REAL(r8), DIMENSION(2,rb%ng,rb%mx*rb%my) :: drzdxq,drzdyq
      REAL(r8), DIMENSION(1,1,1) :: dc
      REAL(r8) :: dx,dy
      INTEGER(i4), DIMENSION(SIZE(pd_arr)) :: use_pd
      INTEGER(i4), DIMENSION(SIZE(pdm_arr)) :: use_pdm
      INTEGER(i4) :: iv,nv,ig,mx,my,ib,ix,iy,ipol,num_bases,iset,       &
     &               num_basesm
!-----------------------------------------------------------------------
!     determine the number of unique basis sets are needed according
!     to the pd_arr array, which lists the polynomial degree of
!     different fields.
!-----------------------------------------------------------------------
      use_pd=1_i4
      DO ib=2,SIZE(pd_arr)
        IF (MINVAL(ABS(pd_arr(1:ib-1)-pd_arr(ib)))==0) use_pd(ib)=0
      ENDDO
      num_bases=SUM(use_pd)
      ALLOCATE(rb%base_pd(num_bases))
!-----------------------------------------------------------------------
!     do the same for modal discontinuous bases.  negative values
!     for the polynomial degree indicated a field that is not used.
!     a unique basis is determined by the pd value and by the minimum
!     and maximum degrees for the limited representation (all three
!     values have to match for two expansions to have the same basis).
!-----------------------------------------------------------------------
      use_pdm=1_i4
      IF (pdm_arr(1)<0) use_pdm(1)=0
      DO ib=2,SIZE(pdm_arr)
        IF (MINVAL(ABS(pdm_arr(1:ib-1)-pdm_arr(ib)))==0.AND.            &
     &      MINVAL(ABS(pdmmin_arr(1:ib-1)-pdmmin_arr(ib)))==0.AND.      &
     &      MINVAL(ABS(pdmmax_arr(1:ib-1)-pdmmax_arr(ib)))==0)          &
     &    use_pdm(ib)=0
        IF (pdm_arr(ib)<0) use_pdm(ib)=0
      ENDDO
      num_basesm=SUM(use_pdm)
      ALLOCATE(rb%base_modal(num_basesm))
!-----------------------------------------------------------------------
!     allocate arrays used for saving the basis function values, basis
!     function gradient values, and grid derivatives at each of the
!     Gaussian quadrature points.
!
!     the quadrature points are now indexed first, and (ix,iy) indices
!     are combined.
!-----------------------------------------------------------------------
      mx=rb%mx
      my=rb%my
      iset=0
      DO ib=1,SIZE(pd_arr)
        IF (use_pd(ib)==0) CYCLE
        iset=iset+1
        rb%base_pd(iset)%poly_deg_basis=pd_arr(ib)
        nv=(pd_arr(ib)+1)**2
        ALLOCATE(rb%base_pd(iset)%alpha(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_pd(iset)%dalpdr(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_pd(iset)%dalpdz(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_pd(iset)%dalpdrc(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_pd(iset)%alpham(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_pd(iset)%dalpmdr(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_pd(iset)%dalpmdz(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_pd(iset)%dalpmdrc(rb%ng,mx*my,nv))
      ENDDO
      iset=0
      DO ib=1,SIZE(pdm_arr)
        IF (use_pdm(ib)==0) CYCLE
        iset=iset+1
        rb%base_modal(iset)%poly_deg_basis=pdm_arr(ib)
        rb%base_modal(iset)%poly_degmin_basis=pdmmin_arr(ib)
        rb%base_modal(iset)%poly_degmax_basis=pdmmax_arr(ib)
        nv=(pdmmax_arr(ib)-pdmmin_arr(ib)+1)*                           &
     &     (2*pdm_arr(ib)+pdmmin_arr(ib)-pdmmax_arr(ib)+1)
        ALLOCATE(rb%base_modal(iset)%alpha(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_modal(iset)%dalpdr(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_modal(iset)%dalpdz(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_modal(iset)%dalpdrc(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_modal(iset)%alpham(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_modal(iset)%dalpmdr(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_modal(iset)%dalpmdz(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_modal(iset)%dalpmdrc(rb%ng,mx*my,nv))
      ENDDO
      ALLOCATE(rb%dxdr(rb%ng,mx*my))
      ALLOCATE(rb%dxdz(rb%ng,mx*my))
      ALLOCATE(rb%dydr(rb%ng,mx*my))
      ALLOCATE(rb%dydz(rb%ng,mx*my))
      ALLOCATE(rb%bigr(rb%ng,mx*my))
      ALLOCATE(rb%wjac(rb%ng,mx*my))
      ALLOCATE(rb%jac2d(rb%ng,mx*my))
!-----------------------------------------------------------------------
!     evaluate and save the basis function weights at the quadrature
!     points.  evaluate the grid derivatives to construct the
!     gradients of the basis functions.  finally, save the quadrature
!     weight times the coordinate-mapping Jacobian for efficiency during
!     finite element computations.
!
!     for toroidal geometry dalpdrc holds d(alpha)/dr + alpha/r, and
!     for linear geometry it's just d(alpha)/dx.
!-----------------------------------------------------------------------
      DO ig=1,rb%ng
        dx=rb%xg(ig)
        dy=rb%yg(ig)
        SELECT CASE(met_spl)
        CASE('pcnst')
          rz=0.25*(rb%rz%fs(:,0:mx-1,0:my-1)                            &
     &            +rb%rz%fs(:,1:mx  ,0:my-1)                            &
     &            +rb%rz%fs(:,0:mx-1,1:my  )                            &
     &            +rb%rz%fs(:,1:mx  ,1:my  ))
          drzdx=0.5*(rb%rz%fs(:,1:mx  ,1:my  )                          &
     &              -rb%rz%fs(:,0:mx-1,1:my  )                          &
     &              +rb%rz%fs(:,1:mx  ,0:my-1)                          &
     &              -rb%rz%fs(:,0:mx-1,0:my-1))
          drzdy=0.5*(rb%rz%fs(:,1:mx  ,1:my  )                          &
     &              -rb%rz%fs(:,1:mx  ,0:my-1)                          &
     &              +rb%rz%fs(:,0:mx-1,1:my  )                          &
     &              -rb%rz%fs(:,0:mx-1,0:my-1))
        CASE('liner','linear','bilinear')
          rz=(1-dx)*(1-dy)*rb%rz%fs(:,0:mx-1,0:my-1)                    &
     &      +   dx *(1-dy)*rb%rz%fs(:,1:mx  ,0:my-1)                    &
     &      +(1-dx)*   dy *rb%rz%fs(:,0:mx-1,1:my  )                    &
     &      +   dx *   dy *rb%rz%fs(:,1:mx  ,1:my  )
          drzdx=   dy *(rb%rz%fs(:,1:mx  ,1:my  )                       &
     &                 -rb%rz%fs(:,0:mx-1,1:my  ))                      &
     &         +(1-dy)*(rb%rz%fs(:,1:mx  ,0:my-1)                       &
     &                 -rb%rz%fs(:,0:mx-1,0:my-1))
          drzdy=   dx *(rb%rz%fs(:,1:mx  ,1:my  )                       &
     &                 -rb%rz%fs(:,1:mx  ,0:my-1))                      &
     &         +(1-dx)*(rb%rz%fs(:,0:mx-1,1:my  )                       &
     &                 -rb%rz%fs(:,0:mx-1,0:my-1))
        CASE('iso','lagrz')
          CALL lagr_quad_all_eval(rb%rz,dx,dy,rz,drzdx,drzdy,1_i4)
        CASE DEFAULT
          CALL nim_stop('Rblock_basis_set: '//TRIM(met_spl)//           &
     &                  ' is not a valid option for met_spl.')
        END SELECT

        ipol=1
        DO iy=1,rb%my
          DO ix=1,rb%mx
            drzdxq(:,ig,ipol)=drzdx(:,ix,iy)
            drzdyq(:,ig,ipol)=drzdy(:,ix,iy)
            rb%bigr(ig,ipol)=rz(1,ix,iy)
            ipol=ipol+1
          ENDDO
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     invert the Jacobian matrix and save partial derivatives.
!-----------------------------------------------------------------------
      CALL math_grid('all',drzdxq,drzdyq,rb%jac2d(:,:),                 &
     &               rb%dxdr(:,:),rb%dxdz(:,:),                         &
     &               rb%dydr(:,:),rb%dydz(:,:))
      rb%wjac(:,:)=rb%jac2d(:,:)
      IF (geom=='tor') THEN
        rb%wjac=rb%wjac*rb%bigr
      ELSE
        rb%bigr=1._r8
      ENDIF
!-----------------------------------------------------------------------
!     multiply the jacobian by the quadrature weight.
!-----------------------------------------------------------------------
      DO ig=1,rb%ng
        rb%wjac(ig,:)=rb%wg(ig)*rb%wjac(ig,:)
      ENDDO
!-----------------------------------------------------------------------
!     load the alpha arrays for all quadrature points with the index
!     reordering.  this is needed for each unique basis set.
!
!     basis/test functions used in the matrix integrand routines are
!     multipied by the square root of (Jacobian*quadrature-weight).
!
!     test functions used in rhs integrand routines are multiplied by
!     Jacobian*quadrature-weight (no square root).
!-----------------------------------------------------------------------
      DO ig=1,rb%ng
        dx=rb%xg(ig)
        dy=rb%yg(ig)
        DO iset=1,num_bases
          nv=(rb%base_pd(iset)%poly_deg_basis+1)**2
          ALLOCATE(alpha(nv),alphax(nv),alphay(nv))
          CALL lagr_quad_bases(dx,dy,alpha,alphax,alphay,1_i4)

          DO iv=1,nv
            rb%base_pd(iset)%alpha (ig,:,iv)=alpha(iv)
            rb%base_pd(iset)%dalpdr(ig,:,iv)=rb%dxdr(ig,:)*alphax(iv)   &
     &                                      +rb%dydr(ig,:)*alphay(iv)
            rb%base_pd(iset)%dalpdz(ig,:,iv)=rb%dxdz(ig,:)*alphax(iv)   &
     &                                      +rb%dydz(ig,:)*alphay(iv)
            IF (geom=='tor') THEN
              rb%base_pd(iset)%dalpdrc(ig,:,iv)=                        &
     &          rb%base_pd(iset)%dalpdr(ig,:,iv)+alpha(iv)/rb%bigr(ig,:)
            ELSE
              rb%base_pd(iset)%dalpdrc(ig,:,iv)=                        &
     &          rb%base_pd(iset)%dalpdr(ig,:,iv)
            ENDIF

            rb%base_pd(iset)%alpham(ig,:,iv)=                           &
     &        rb%base_pd(iset)%alpha(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_pd(iset)%dalpmdr(ig,:,iv)=                          &
     &        rb%base_pd(iset)%dalpdr(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_pd(iset)%dalpmdz(ig,:,iv)=                          &
     &        rb%base_pd(iset)%dalpdz(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_pd(iset)%dalpmdrc(ig,:,iv)=                         &
     &        rb%base_pd(iset)%dalpdrc(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_pd(iset)%alpha(ig,:,iv)=                            &
     &        rb%base_pd(iset)%alpha(ig,:,iv)*rb%wjac(ig,:)
            rb%base_pd(iset)%dalpdr(ig,:,iv)=                           &
     &        rb%base_pd(iset)%dalpdr(ig,:,iv)*rb%wjac(ig,:)
            rb%base_pd(iset)%dalpdz(ig,:,iv)=                           &
     &        rb%base_pd(iset)%dalpdz(ig,:,iv)*rb%wjac(ig,:)
            rb%base_pd(iset)%dalpdrc(ig,:,iv)=                          &
     &        rb%base_pd(iset)%dalpdrc(ig,:,iv)*rb%wjac(ig,:)
          ENDDO
          DEALLOCATE(alpha,alphax,alphay)
        ENDDO

        DO iset=1,num_basesm
          nv=(rb%base_modal(iset)%poly_degmax_basis                     &
     &       -rb%base_modal(iset)%poly_degmin_basis+1)*                 &
     &       (2*rb%base_modal(iset)%poly_deg_basis+1                    &
     &       -rb%base_modal(iset)%poly_degmax_basis                     &
     &       +rb%base_modal(iset)%poly_degmin_basis)
          ALLOCATE(alpha(nv),alphax(nv),alphay(nv))
          CALL modal_disc_bases(dx,dy,alpha,alphax,alphay,1_i4,         &
     &                          rb%base_modal(iset)%poly_deg_basis,     &
     &                          rb%base_modal(iset)%poly_degmin_basis,  &
     &                          rb%base_modal(iset)%poly_degmax_basis)

          DO iv=1,nv
            rb%base_modal(iset)%alpha (ig,:,iv)=alpha(iv)
            rb%base_modal(iset)%dalpdr(ig,:,iv)=rb%dxdr(ig,:)*alphax(iv)&
     &                                         +rb%dydr(ig,:)*alphay(iv)
            rb%base_modal(iset)%dalpdz(ig,:,iv)=rb%dxdz(ig,:)*alphax(iv)&
     &                                         +rb%dydz(ig,:)*alphay(iv)
            IF (geom=='tor') THEN
              rb%base_modal(iset)%dalpdrc(ig,:,iv)=                     &
     &          rb%base_modal(iset)%dalpdr(ig,:,iv)+                    &
     &          alpha(iv)/rb%bigr(ig,:)
            ELSE
              rb%base_modal(iset)%dalpdrc(ig,:,iv)=                     &
     &          rb%base_modal(iset)%dalpdr(ig,:,iv)
            ENDIF

            rb%base_modal(iset)%alpham(ig,:,iv)=                        &
     &        rb%base_modal(iset)%alpha(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_modal(iset)%dalpmdr(ig,:,iv)=                       &
     &        rb%base_modal(iset)%dalpdr(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_modal(iset)%dalpmdz(ig,:,iv)=                       &
     &        rb%base_modal(iset)%dalpdz(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_modal(iset)%dalpmdrc(ig,:,iv)=                      &
     &        rb%base_modal(iset)%dalpdrc(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_modal(iset)%alpha(ig,:,iv)=                         &
     &        rb%base_modal(iset)%alpha(ig,:,iv)*rb%wjac(ig,:)
            rb%base_modal(iset)%dalpdr(ig,:,iv)=                        &
     &        rb%base_modal(iset)%dalpdr(ig,:,iv)*rb%wjac(ig,:)
            rb%base_modal(iset)%dalpdz(ig,:,iv)=                        &
     &        rb%base_modal(iset)%dalpdz(ig,:,iv)*rb%wjac(ig,:)
            rb%base_modal(iset)%dalpdrc(ig,:,iv)=                       &
     &        rb%base_modal(iset)%dalpdrc(ig,:,iv)*rb%wjac(ig,:)
          ENDDO
          DEALLOCATE(alpha,alphax,alphay)
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_basis_set
!-----------------------------------------------------------------------
!     subprogram 9. rblock_bicube_set.
!     evaluate bicubic spline data at the gaussian quadrature points for
!     this block.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_bicube_set(bc,qbc,rb,d_order)
      USE bicube

      TYPE(bicube_type), INTENT(INOUT) :: bc
      TYPE(rb_real_qp_type), INTENT(OUT) :: qbc
      TYPE(rblock_type), INTENT(INOUT) :: rb
      CHARACTER(*), INTENT(IN) :: d_order

      INTEGER(i4) :: ig,bcmode,mxb,myb,iq,ix,iy,ipol
      REAL(r8), DIMENSION(bc%nqty,rb%mx,rb%my) :: f,fx,fy
      REAL(r8), DIMENSION(1,1,1) :: db
      REAL(r8) :: dx,dy
!-----------------------------------------------------------------------
!     allocate quadrature storage space as needed.
!-----------------------------------------------------------------------
      mxb=rb%mx
      myb=rb%my
      ALLOCATE(qbc%qpf(bc%nqty,rb%ng,mxb*myb))
      IF (d_order/='values') THEN
        ALLOCATE(qbc%qpfr(bc%nqty,rb%ng,mxb*myb))
        ALLOCATE(qbc%qpfz(bc%nqty,rb%ng,mxb*myb))
      ELSE
        ALLOCATE(qbc%qpfr(1,1,1))
        ALLOCATE(qbc%qpfz(1,1,1))
      ENDIF
!-----------------------------------------------------------------------
!     set evaluation mode.
!-----------------------------------------------------------------------
      SELECT CASE(d_order)
      CASE('values')
        bcmode=0
      CASE('1st derivs')
        bcmode=1
      CASE DEFAULT
        CALL nim_stop                                                   &
     &    ('Unrecognized derivative order in rblock_bicube_set.')
      END SELECT
!-----------------------------------------------------------------------
!     loop over quadrature points, evaluate, and save.
!-----------------------------------------------------------------------
      DO ig=1,rb%ng
        dx=rb%xg(ig)
        dy=rb%yg(ig)
        CALL bicube_all_eval(bc,dx,dy,f,fx,fy,db,db,db,bcmode)
        ipol=1
        DO iy=1,myb
          DO ix=1,mxb
            qbc%qpf(:,ig,ipol)=f(:,ix,iy)
            ipol=ipol+1
          ENDDO
        ENDDO
        IF (bcmode>0) THEN
          ipol=1
          DO iy=1,myb
            DO ix=1,mxb
              qbc%qpfr(:,ig,ipol)=rb%dxdr(ig,ipol)*fx(:,ix,iy)          &
     &                           +rb%dydr(ig,ipol)*fy(:,ix,iy)
              qbc%qpfz(:,ig,ipol)=rb%dxdz(ig,ipol)*fx(:,ix,iy)          &
     &                           +rb%dydz(ig,ipol)*fy(:,ix,iy)
              ipol=ipol+1
            ENDDO
          ENDDO
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
!     deallocate bicube evaluation matrix.
!-----------------------------------------------------------------------
      DEALLOCATE(bc%cmats)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_bicube_set
!-----------------------------------------------------------------------
!     subprogram 10. rblock_real_qp_update.
!     evaluate 2D lagrange_quad data at the gaussian quadrature points
!     for this block.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_real_qp_update(laq,qlaq,rb)

      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq
      TYPE(rb_real_qp_type), INTENT(INOUT) :: qlaq
      TYPE(rblock_type), INTENT(IN) :: rb

      INTEGER(i4) :: ig,mxb,myb,iq,ix,iy,ipol
      REAL(r8), DIMENSION(laq%nqty,rb%mx,rb%my) :: f,fx,fy
      REAL(r8) :: dx,dy
!-----------------------------------------------------------------------
!     loop over quadrature points, evaluate, and save.
!-----------------------------------------------------------------------
      mxb=rb%mx
      myb=rb%my
      DO ig=1,rb%ng
        dx=rb%xg(ig)
        dy=rb%yg(ig)
        CALL lagr_quad_all_eval(laq,dx,dy,f,fx,fy,1_i4)
        ipol=1
        DO iy=1,myb
          DO ix=1,mxb
            qlaq%qpf (:,ig,ipol)=f(:,ix,iy)
            qlaq%qpfr(:,ig,ipol)=rb%dxdr(ig,ipol)*fx(:,ix,iy)           &
     &                          +rb%dydr(ig,ipol)*fy(:,ix,iy)
            qlaq%qpfz(:,ig,ipol)=rb%dxdz(ig,ipol)*fx(:,ix,iy)           &
     &                          +rb%dydz(ig,ipol)*fy(:,ix,iy)
            ipol=ipol+1
          ENDDO
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_real_qp_update
!-----------------------------------------------------------------------
!     subprogram 11. rblock_comp_qp_update.
!     evaluate 3D lagrange_quad data at the gaussian quadrature points
!     for this block.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_comp_qp_update(laq,qlaq,rb)

      TYPE(lagr_quad_type), INTENT(IN) :: laq
      TYPE(rb_comp_qp_type), INTENT(INOUT) :: qlaq
      TYPE(rblock_type), INTENT(IN) :: rb

      INTEGER(i4) :: ig,mxb,myb,iq,im,ix,iy,ipol
      COMPLEX(r8), DIMENSION(laq%nqty,rb%mx,rb%my,laq%nfour) :: f,fx,fy
      REAL(r8) :: dx,dy
!-----------------------------------------------------------------------
!     loop over quadrature points, evaluate, and save.
!-----------------------------------------------------------------------
      mxb=rb%mx
      myb=rb%my
      DO ig=1,rb%ng
        dx=rb%xg(ig)
        dy=rb%yg(ig)
        CALL lagr_quad_all_eval(laq,dx,dy,f,fx,fy,1_i4)
        DO im=1,laq%nfour
          ipol=1
          DO iy=1,myb
            DO ix=1,mxb
              qlaq%qpf (:,ig,ipol,im)=f(:,ix,iy,im)
              qlaq%qpfr(:,ig,ipol,im)=rb%dxdr(ig,ipol)*fx(:,ix,iy,im)   &
     &                               +rb%dydr(ig,ipol)*fy(:,ix,iy,im)
              qlaq%qpfz(:,ig,ipol,im)=rb%dxdz(ig,ipol)*fx(:,ix,iy,im)   &
     &                               +rb%dydz(ig,ipol)*fy(:,ix,iy,im)
              ipol=ipol+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_comp_qp_update
!-----------------------------------------------------------------------
!     subprogram 12. rblock_real_qp_alloc.
!     allocate space for 2D data at the gaussian quadrature points for
!     this block.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_real_qp_alloc(qlaq,rb,nqty)

      TYPE(rb_real_qp_type), INTENT(OUT) :: qlaq
      TYPE(rblock_type), INTENT(IN) :: rb
      INTEGER(i4), INTENT(IN) :: nqty

      ALLOCATE(qlaq%qpf (nqty,rb%ng,rb%mx*rb%my)); qlaq%qpf=0._r8
      ALLOCATE(qlaq%qpfr(nqty,rb%ng,rb%mx*rb%my)); qlaq%qpfr=0._r8
      ALLOCATE(qlaq%qpfz(nqty,rb%ng,rb%mx*rb%my)); qlaq%qpfz=0._r8
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_real_qp_alloc
!-----------------------------------------------------------------------
!     subprogram 13. rblock_comp_qp_alloc.
!     allocate space for 3D data at the gaussian quadrature points for
!     this block.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_comp_qp_alloc(qlaq,rb,nqty,nfour)

      TYPE(rb_comp_qp_type), INTENT(OUT) :: qlaq
      TYPE(rblock_type), INTENT(IN) :: rb
      INTEGER(i4), INTENT(IN) :: nqty,nfour

      ALLOCATE(qlaq%qpf(nqty,rb%ng,rb%mx*rb%my,nfour)); qlaq%qpf=0._r8
      ALLOCATE(qlaq%qpfr(nqty,rb%ng,rb%mx*rb%my,nfour)); qlaq%qpfr=0._r8
      ALLOCATE(qlaq%qpfz(nqty,rb%ng,rb%mx*rb%my,nfour)); qlaq%qpfz=0._r8
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_comp_qp_alloc
!-----------------------------------------------------------------------
!     subprogram 14. rblock_real_qp_dealloc.
!     deallocate space for 2D data at the gaussian quadrature points for
!     this block.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_real_qp_dealloc(qlaq)

      TYPE(rb_real_qp_type), INTENT(INOUT) :: qlaq

      DEALLOCATE(qlaq%qpf)
      IF (ALLOCATED(qlaq%qpfr)) DEALLOCATE(qlaq%qpfr)
      IF (ALLOCATED(qlaq%qpfz)) DEALLOCATE(qlaq%qpfz)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_real_qp_dealloc
!-----------------------------------------------------------------------
!     subprogram 15. rblock_comp_qp_dealloc.
!     deallocate space for 3D data at the gaussian quadrature points for
!     this block.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_comp_qp_dealloc(qlaq)

      TYPE(rb_comp_qp_type), INTENT(INOUT) :: qlaq

      DEALLOCATE(qlaq%qpf)
      IF (ALLOCATED(qlaq%qpfr)) DEALLOCATE(qlaq%qpfr)
      IF (ALLOCATED(qlaq%qpfz)) DEALLOCATE(qlaq%qpfz)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_comp_qp_dealloc
!-----------------------------------------------------------------------
!     subprogram 7. qp_fft_save.                                        
!     save real data as a function of toroidal angle at quadrature      
!     points to reduce the number of ffts called during 3D matrix       
!     iterations.                                                       
!-----------------------------------------------------------------------
      SUBROUTINE qp_fft_save(fcmp,fphi,lx,ly,mps,nq,ng,feq) 
      USE local 
      USE global 
      USE input 
      USE fft_mod 
      IMPLICIT NONE 
                                                                        
      INTEGER(i4), INTENT(IN) :: lx,ly,mps,nq,ng 
      REAL(r8), DIMENSION(nq,mps,nphi), INTENT(INOUT) :: fphi 
      REAL(r8), DIMENSION(nq,ng,lx*ly), INTENT(IN) :: feq 
      COMPLEX(r8), DIMENSION(nq,ng,lx*ly,nmodes), INTENT(IN) :: fcmp 
                                                                        
      INTEGER(i4) :: ig,im,iq,ix,iy,nqm,ipol,iphi 
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: fcmp2 
!-----------------------------------------------------------------------
!     reordering the quadrature-point index eliminates the data shuffle 
!     that was here, and mps includes a combination of poloidal indices 
!     and quadrature-point indices.                                     
!-----------------------------------------------------------------------
      ALLOCATE(fcmp2(nq,ng,lx*ly,nmodes)) 
      DO im=1,nmodes 
        IF (keff(im)==0) THEN 
          fcmp2(:,:,:,im)=fcmp(:,:,:,im)+feq 
        ELSE 
          fcmp2(:,:,:,im)=fcmp(:,:,:,im) 
        ENDIF 
      ENDDO 
!-----------------------------------------------------------------------
!     perform ffts of the packed data.                                  
!-----------------------------------------------------------------------
      CALL fft_nim('inverse',ng*lx*ly,mps,nphi,nq,fcmp2,fphi,dealiase) 
      DEALLOCATE(fcmp2) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE qp_fft_save 
!-----------------------------------------------------------------------
!     subprogram 7. qp_fft_save_zflr                                     
!     save real data as a function of toroidal angle at quadrature      
!     points to reduce the number of ffts called during 3D matrix       
!     iterations.                                                       
!-----------------------------------------------------------------------
      SUBROUTINE qp_fft_save_zflr(fcmp,fphi,lx,ly,mps,nq,ng,feq,zflr) 
      USE local 
      USE global 
      USE input 
      USE fft_mod 
      IMPLICIT NONE 
                                                                        
      INTEGER(i4), INTENT(IN) :: lx,ly,mps,nq,ng 
      REAL(r8), DIMENSION(nq,mps,nphi), INTENT(INOUT) :: fphi 
      REAL(r8), DIMENSION(nq,ng,lx*ly), INTENT(IN) :: feq 
      COMPLEX(r8), DIMENSION(nq,ng,lx*ly,nmodes), INTENT(IN) :: fcmp 
      LOGICAL, INTENT(IN) :: zflr
                                                                        
      INTEGER(i4) :: ig,im,iq,ix,iy,nqm,ipol,iphi 
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: fcmp2 
!-----------------------------------------------------------------------
!     reordering the quadrature-point index eliminates the data shuffle 
!     that was here, and mps includes a combination of poloidal indices 
!     and quadrature-point indices.                                     
!-----------------------------------------------------------------------
      ALLOCATE(fcmp2(nq,ng,lx*ly,nmodes)) 
      DO im=1,nmodes 
        IF (keff(im)==0) THEN 
          fcmp2(:,:,:,im)=fcmp(:,:,:,im)+feq 
        ELSE 
          fcmp2(:,:,:,im)=fcmp(:,:,:,im) 
        ENDIF 
      ENDDO 
!-----------------------------------------------------------------------
!     perform ffts of the packed data.                                  
!-----------------------------------------------------------------------
      CALL fft_nim('inverse',ng*lx*ly,mps,nphi,nq,fcmp2,fphi,dealiase) 
      fphi=MAX(fphi,smallnum)
      DEALLOCATE(fcmp2) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE qp_fft_save_zflr
!-----------------------------------------------------------------------
!     subprogram 8. qp_fft_noeq_save.                                  
!     save real data as a function of toroidal angle at quadrature      
!     points to reduce the number of ffts called during 3D matrix       
!     iterations.  this version does not have an 'equilibrium'          
!     component to be added.                                            
!-----------------------------------------------------------------------
      SUBROUTINE qp_fft_noeq_save(fcmp,fphi,lx,ly,mps,nq,ng) 
      USE local 
      USE global 
      USE input 
      USE fft_mod 
      IMPLICIT NONE 
                                                                        
      INTEGER(i4), INTENT(IN) :: lx,ly,mps,nq,ng 
      REAL(r8), DIMENSION(nq,mps,nphi), INTENT(INOUT) :: fphi 
      COMPLEX(r8), DIMENSION(nq,ng,lx*ly,nmodes), INTENT(INOUT) :: fcmp 
!-----------------------------------------------------------------------
!     reordering the quadrature-point index eliminates the data shuffle 
!     that was here, and mps includes a combination of poloidal indices 
!     and quadrature-point indices.  this routine now calls the fft     
!     routine directly.                                                 
!-----------------------------------------------------------------------
      CALL fft_nim('inverse',ng*lx*ly,mps,nphi,nq,fcmp,fphi,dealiase) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE qp_fft_noeq_save
!-----------------------------------------------------------------------
!     subprogram 8. qp_fft_noeq_save_zflr.                                  
!     save real data as a function of toroidal angle at quadrature      
!     points to reduce the number of ffts called during 3D matrix       
!     iterations.  this version does not have an 'equilibrium'          
!     component to be added.                                            
!-----------------------------------------------------------------------
      SUBROUTINE qp_fft_noeq_save_zflr(fcmp,fphi,lx,ly,mps,nq,ng,zflr) 
      USE local 
      USE global 
      USE input 
      USE fft_mod 
      IMPLICIT NONE 
                                                                        
      INTEGER(i4), INTENT(IN) :: lx,ly,mps,nq,ng 
      REAL(r8), DIMENSION(nq,mps,nphi), INTENT(INOUT) :: fphi 
      COMPLEX(r8), DIMENSION(nq,ng,lx*ly,nmodes), INTENT(INOUT) :: fcmp 
      LOGICAL, INTENT(IN) :: zflr
!-----------------------------------------------------------------------
!     reordering the quadrature-point index eliminates the data shuffle 
!     that was here, and mps includes a combination of poloidal indices 
!     and quadrature-point indices.  this routine now calls the fft     
!     routine directly.                                                 
!-----------------------------------------------------------------------
      CALL fft_nim('inverse',ng*lx*ly,mps,nphi,nq,fcmp,fphi,dealiase) 
      fphi=MAX(fphi,smallnum)
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE qp_fft_noeq_save_zflr
!-----------------------------------------------------------------------
!     subprogram 16. rblock_real_qpe_update.
!     evaluate 2D lagrange_quad data at the gaussian quadrature points
!     for this block.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_real_qpe_update(lae,qlae,rb)

      TYPE(lagr_edge_1D_type), INTENT(IN) :: lae
      TYPE(rb_real_qpe_type), INTENT(INOUT) :: qlae
      TYPE(rblock_type), INTENT(IN) :: rb

      INTEGER(i4) :: ig,iq
      REAL(r8), DIMENSION(2,rb%mxe) :: rz, drzds
      REAL(r8), DIMENSION(lae%nqty,rb%mxe) :: fx
      REAL(r8) :: dx
!-----------------------------------------------------------------------
!     loop over quadrature points, evaluate, and save.
!-----------------------------------------------------------------------
      DO ig=1,rb%nge
        dx=rb%xge(ig)
        CALL lagr_edge_all_eval(lae,dx,qlae%qpf(:,:,ig),fx,1_i4)
        CALL lagr_edge_all_eval(rb%rzedge,dx,rz,drzds,1_i4)
        FORALL(iq=1:lae%nqty)
         WHERE(drzds(1,:).NE.0) qlae%qpfr(iq,:,ig)=fx(iq,:)/drzds(1,:)
         WHERE(drzds(2,:).NE.0) qlae%qpfz(iq,:,ig)=fx(iq,:)/drzds(2,:)
        ENDFORALL
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_real_qpe_update
!-----------------------------------------------------------------------
!     subprogram 16. rblock_comp_qpe_update.
!     evaluate 3D lagrange_quad data at the gaussian quadrature points
!     for this block.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_comp_qpe_update(lae,qlae,rb)

      TYPE(lagr_edge_type), INTENT(IN) :: lae
      TYPE(rb_comp_qpe_type), INTENT(INOUT) :: qlae
      TYPE(rblock_type), INTENT(IN) :: rb

      INTEGER(i4) :: ig,iq,im
      COMPLEX(r8), DIMENSION(lae%nqty,rb%mxe,lae%nfour) :: fx,fy
      REAL(r8), DIMENSION(2,rb%mxe) :: rz, drzds
      REAL(r8) :: dx
!-----------------------------------------------------------------------
!     loop over quadrature points, evaluate, and save.
!-----------------------------------------------------------------------
      DO ig=1,rb%nge
        dx=rb%xge(ig)
        CALL lagr_edge_all_eval(lae,dx,qlae%qpf(:,:,:,ig),fx,1_i4)
        CALL lagr_edge_all_eval(rb%rzedge,dx,rz,drzds,1_i4)
        FORALL(im=1:lae%nfour,iq=1:lae%nqty)
          WHERE (drzds(1,:).NE.0)                                       &
     &       qlae%qpfr(iq,:,im,ig)=fx(iq,:,im)/drzds(1,:)
          WHERE (drzds(2,:).NE.0)                                       &
     &       qlae%qpfz(iq,:,im,ig)=fx(iq,:,im)/drzds(2,:)
        ENDFORALL
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_comp_qpe_update
!-----------------------------------------------------------------------
!     subprogram 12. rblock_real_qpe_alloc.
!     allocate space for real edge data at the gaussian quadrature point
!     for this block.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_real_qpe_alloc(qlae,rb,nqty)

      TYPE(rb_real_qpe_type), INTENT(OUT) :: qlae
      TYPE(rblock_type), INTENT(IN) :: rb
      INTEGER(i4), INTENT(IN) :: nqty

      ALLOCATE(qlae%qpf(nqty,rb%mxe,rb%nge));  qlae%qpf=0._r8
      ALLOCATE(qlae%qpfr(nqty,rb%mxe,rb%nge)); qlae%qpfr=0._r8
      ALLOCATE(qlae%qpfz(nqty,rb%mxe,rb%nge)); qlae%qpfz=0._r8
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_real_qpe_alloc
!-----------------------------------------------------------------------
!     subprogram 18. rblock_comp_qpe_alloc.
!     allocate space for 3D data at the gaussian quadrature points for
!     this block.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_comp_qpe_alloc(qlae,rb,nqty,nfour)

      TYPE(rb_comp_qpe_type), INTENT(OUT) :: qlae
      TYPE(rblock_type), INTENT(IN) :: rb
      INTEGER(i4), INTENT(IN) :: nqty,nfour

      ALLOCATE(qlae%qpf( nqty,rb%mxe,nfour,rb%nge)); qlae%qpf=0._r8
      ALLOCATE(qlae%qpfr(nqty,rb%mxe,nfour,rb%nge)); qlae%qpfr=0._r8
      ALLOCATE(qlae%qpfz(nqty,rb%mxe,nfour,rb%nge)); qlae%qpfz=0._r8
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_comp_qpe_alloc
!-----------------------------------------------------------------------
!     subprogram 19. rblock_real_qpe_dealloc.
!     deallocate space for 2D data at the gaussian quadrature points for
!     this block.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_real_qpe_dealloc(qlae)

      TYPE(rb_real_qpe_type), INTENT(INOUT) :: qlae

      DEALLOCATE(qlae%qpf)
      IF (ALLOCATED(qlae%qpfr)) DEALLOCATE(qlae%qpfr)
      IF (ALLOCATED(qlae%qpfz)) DEALLOCATE(qlae%qpfz)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_real_qpe_dealloc
!-----------------------------------------------------------------------
!     subprogram 20. rblock_comp_qpe_dealloc.
!     deallocate space for 3D data at the gaussian quadrature points for
!     this block.
!-----------------------------------------------------------------------
      SUBROUTINE rblock_comp_qpe_dealloc(qlae)

      TYPE(rb_comp_qpe_type), INTENT(INOUT) :: qlae

      DEALLOCATE(qlae%qpf)
      IF (ALLOCATED(qlae%qpfr)) DEALLOCATE(qlae%qpfr)
      IF (ALLOCATED(qlae%qpfz)) DEALLOCATE(qlae%qpfz)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_comp_qpe_dealloc
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      END MODULE rblock
