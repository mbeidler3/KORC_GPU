!-----------------------------------------------------------------------
!     $Id: vector_type_mod.F90 6130 2018-07-20 22:03:49Z jking $
!     contains a module that defines structures for the linear algebra
!     vectors.
!
!     this has been broken into separate modules for different data
!     types to facilitate compilation on some machines.
!-----------------------------------------------------------------------
#include "config.f"
!-----------------------------------------------------------------------
!     code organization.
!-----------------------------------------------------------------------
!     module vector_defn_type_mod
!     module rvector_type_mod
!     1. vector_rtype_alloc.
!     2. vector_rtype_dealloc.
!     3. vector_assign_rsc.
!     4. vector_assign_csc.
!     5. vector_assign_int.
!     6. vector_assign_vec.
!     7. vector_assignq_vec
!     8. vector_assignp_vec
!     9. vector_assignq_cvec2
!     10. vector_assign_cvec.
!     11. vector_assign_vec3.
!     12. vector_ptassign_bc.
!     13. vector_ptassign_laq2.
!     13.1. vector_ptassign_modq2.
!     14. vector_ptassign_tl2.
!     15. vector_add_vec.
!     16. vector_mult_rsc.
!     17. vector_mult_int.
!     module cvector_type_mod
!     18. vector_ctype_alloc.
!     19. vector_ctype_dealloc.
!     20. cvector_assign_rsc.
!     20.1 cvector_assign_rsc_nq.
!     21. cvector_assign_csc.
!     22. cvector_assign_int.
!     23. cvector_assign_cvec.
!     24. cvector_assignq_cvec.
!     25. cvector_assignp_cvec.
!     26. cvector_assign_vec.
!     27. cvector_assign_cvec2.
!     28. cvector_assignp_cvec2.
!     29. cvector_ptassign_laq.
!     29.1. cvector_ptassign_modq.
!     30. cvector_ptassign_tl.
!     31. cvector_add_cvec.
!     32. cvector_add_iq.
!     32.1 cvector_add_nq.
!     33. cvector_addc_cvec.
!     34. cvector_mult_rsc.
!     35. cvector_mult_csc.
!     36. cvector_mult_int.
!     37. cvector_real_comp.
!     module cvector_2D_type_mod
!     38. vector_2D_ctype_alloc.
!     39. vector_2D_ctype_dealloc.
!     40. cvector_2D_assign_rsc.
!     41. cvector_2D_assign_csc.
!     42. cvector_2D_assign_int.
!     43. cvector_2D_assign_cvec2.
!     44. cvector_2D_assignq_cvec2
!     45. cvector_2D_assignp_cvec2
!     46. cvector_2D_assign_vec
!     47. cvector_2D_assignq_vec
!     48. cvector_2D_assignp_vec.
!     49. cvector_2D_assign_cvec.
!     50. cvector_2D_assignp_cvec.
!     51. cvector_2D_add_cvec2.
!     52. cvector_2D_addc_cvec2.
!     53. cvector_2D_mult_rsc.
!     54. cvector_2D_mult_csc.
!     55. cvector_2D_mult_int.
!     55.1. cvector_2D_pack_cvec.
!     55.2. cvector_2D_pack_cvec2.
!     55.3. cvector_2D_unpack_cvec.
!     55.4. cvector_2D_unpack_add_cvec.
!     55.5. cvector_2D_unpack_cvec2.
!     module rvector_3D_type_mod
!     56. vector_3D_rtype_alloc.
!     57. vector_3D_rtype_dealloc.
!     58. vector_3D_assign_rsc.
!     59. vector_3D_assign_csc.
!     60. vector_3D_assign_int.
!     61. vector_3D_assign_vec.
!     62. vector_3D_assign_vec3.
!     63. vector_3D_add_vec.
!     64. vector_3D_add_vec3.
!     65. vector_3D_mult_rsc.
!     66. vector_3D_mult_int.
!     module vector_type_mod
!-----------------------------------------------------------------------
!     module for type definitions.
!-----------------------------------------------------------------------
      MODULE vector_defn_type_mod
      USE local
      IMPLICIT NONE
!-----------------------------------------------------------------------
!     the vector_type is set-up for 2D arrays of real vector quantities
!     with element side and interior centerings, as well as grid
!     vertices.
!
!     the arrtmp arrays are extra space for temporary fields
!     and are not affected by assignment and algebraic operations.
!-----------------------------------------------------------------------
      TYPE :: vector_type
        REAL(r8), DIMENSION(:,:,:), POINTER :: arr
        REAL(r8), DIMENSION(:,:,:,:), POINTER :: arrh,arrv,arri,arrtmp
        INTEGER(i4), POINTER :: mem_id
      END TYPE vector_type
!-----------------------------------------------------------------------
!     the cvector_type is set-up for 3D arrays of complex vector
!     quantities with element side and interior centerings, as well
!     as grid vertices.
!-----------------------------------------------------------------------
      TYPE :: cvector_type
        COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: arr
        COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: arrh,arrv,arri,   &
     &               arrtmp
        INTEGER(i4), POINTER :: mem_id
      END TYPE cvector_type
!-----------------------------------------------------------------------
!     the cvector_2D_type is 2D complex vector array used for interim
!     computations.
!-----------------------------------------------------------------------
      TYPE :: cvector_2D_type
        COMPLEX(r8), DIMENSION(:,:,:), POINTER :: arr
        COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: arrh,arrv,arri,     &
     &               arrtmp
        INTEGER(i4), POINTER :: mem_id
      END TYPE cvector_2D_type
!-----------------------------------------------------------------------
!     the vector_3D_type is set-up for 3D arrays of real vector
!     quantities with the last index running over indices for the
!     periodic coordinate.
!-----------------------------------------------------------------------
      TYPE :: vector_3D_type
        REAL(r8), DIMENSION(:,:,:,:), POINTER :: arr
        REAL(r8), DIMENSION(:,:,:,:,:), POINTER :: arrh,arrv,arri,      &
     &               arrtmp
        INTEGER(i4), POINTER :: mem_id
      END TYPE vector_3D_type
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      END MODULE vector_defn_type_mod

!-----------------------------------------------------------------------
!     module for 2D-in-space arrays of real vectors quantities.
!-----------------------------------------------------------------------
      MODULE rvector_type_mod
      USE vector_defn_type_mod
      IMPLICIT NONE

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 1. vector_rtype_alloc.
!     allocates space for a real vector_type structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_rtype_alloc(rvt,poly_degree,mx,my,nqty,nbt,nqt)
#ifdef OBJ_MEM_PROF
      USE memlog, only: memlogger
#endif

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,poly_degree
      INTEGER(i4), INTENT(IN), OPTIONAL :: nbt,nqt
      TYPE(vector_type), INTENT(OUT) :: rvt

      INTEGER(i4) :: sz
!-----------------------------------------------------------------------
!     allocate space according to the basis functions needed.
!     if poly_degree is non-positive, storage for discontinuous-field
!     coefficients is allocated.
!
!-PRE triangles will need something here, too.
!-----------------------------------------------------------------------
      SELECT CASE(poly_degree)
      CASE(:0)  !  piecewise continuous fields
        ALLOCATE(rvt%arri(nqty,(poly_degree-1)**2,mx,my))
        NULLIFY(rvt%arr,rvt%arrh,rvt%arrv)
      CASE(1)  !  linear elements
        ALLOCATE(rvt%arr(nqty,0:mx,0:my))
        NULLIFY(rvt%arri,rvt%arrh,rvt%arrv)
      CASE(2:) !  higher-order elements
        ALLOCATE(rvt%arr(nqty,0:mx,0:my))
        ALLOCATE(rvt%arrh(nqty,poly_degree-1,1:mx,0:my))
        ALLOCATE(rvt%arrv(nqty,poly_degree-1,0:mx,1:my))
        ALLOCATE(rvt%arri(nqty,(poly_degree-1)**2,1:mx,1:my))
      END SELECT
!-----------------------------------------------------------------------
!     if the optional input for temporary arrays is present, allocate
!     the arrtmp arrays; otherwise, nullify them.  note that nbt is
!     the number of temporary bases, and nqt is the number of quantities
!     at each basis.
!-----------------------------------------------------------------------
      IF (PRESENT(nbt)) THEN
        ALLOCATE(rvt%arrtmp(nqt,nbt,1:mx,1:my))
      ELSE
        NULLIFY(rvt%arrtmp)
      ENDIF
!-----------------------------------------------------------------------
!     register this object.
!-----------------------------------------------------------------------
      NULLIFY(rvt%mem_id)
#ifdef OBJ_MEM_PROF
      SELECT CASE(poly_degree)
      CASE(0)
        sz=SIZEOF(rvt%arri)
      CASE(1)
        sz=SIZEOF(rvt%arr)
      CASE(2:)
        sz= SIZEOF(rvt%arri)+SIZEOF(rvt%arr)                            &
     &     +SIZEOF(rvt%arrh)+SIZEOF(rvt%arrv)
      END SELECT
      CALL memlogger%update(rvt%mem_id,'rvector_2D','unknown',sz)
#endif
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_rtype_alloc
!-----------------------------------------------------------------------
!     subprogram 2. vector_rtype_dealloc.
!     deallocates space for a real vector_type structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_rtype_dealloc(rvt)
#ifdef OBJ_MEM_PROF
      USE memlog, only: memlogger
#endif

      TYPE(vector_type), INTENT(INOUT) :: rvt

      IF (ASSOCIATED(rvt%arr)) THEN
        DEALLOCATE(rvt%arr)
        NULLIFY(rvt%arr)
      ENDIF
      IF (ASSOCIATED(rvt%arrh)) THEN
        DEALLOCATE(rvt%arrh)
        NULLIFY(rvt%arrh)
      ENDIF
      IF (ASSOCIATED(rvt%arrv)) THEN
        DEALLOCATE(rvt%arrv)
        NULLIFY(rvt%arrv)
      ENDIF
      IF (ASSOCIATED(rvt%arri)) THEN
        DEALLOCATE(rvt%arri)
        NULLIFY(rvt%arri)
      ENDIF
      IF (ASSOCIATED(rvt%arrtmp)) THEN
        DEALLOCATE(rvt%arrtmp)
        NULLIFY(rvt%arrtmp)
      ENDIF
!-----------------------------------------------------------------------
!     unregister this object.
!-----------------------------------------------------------------------
#ifdef OBJ_MEM_PROF
      CALL memlogger%update(rvt%mem_id,' ',' ',0,resize=.TRUE.)
#endif
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_rtype_dealloc
!-----------------------------------------------------------------------
!     subprogram 3. vector_assign_rsc.
!     assign a real scalar value to a vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_assign_rsc(vec,rscalar)

      TYPE(vector_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rscalar

!-----------------------------------------------------------------------
!     if the grid vertex-centered array is allocated, treat as a
!     standard element.  If, not the structure represents piecewise
!     constant.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=rscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rscalar
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_assign_rsc: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assign_rsc
!-----------------------------------------------------------------------
!     subprogram 4. vector_assign_csc.
!     assign a complex scalar value to a vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_assign_csc(vec,cscalar)

      TYPE(vector_type), INTENT(INOUT) :: vec
      COMPLEX(r8), INTENT(IN) :: cscalar

!-----------------------------------------------------------------------
!     if the grid vertex-centered array is allocated, treat as a
!     standard element.  If, not the structure represents piecewise
!     constant.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=cscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=cscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=cscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=cscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=cscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=cscalar
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_assign_csc: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assign_csc
!-----------------------------------------------------------------------
!     subprogram 5. vector_assign_int.
!     assign a integer value to a vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_assign_int(vec,int)

      TYPE(vector_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

!-----------------------------------------------------------------------
!     if the grid vertex-centered array is allocated, treat as a
!     standard element.  If, not the structure represents piecewise
!     constant.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int
        IF (ASSOCIATED(vec%arri)) vec%arri=int
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_assign_int: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assign_int
!-----------------------------------------------------------------------
!     subprogram 6. vector_assign_vec.
!     set one vector structure equal to another.
!-----------------------------------------------------------------------
      SUBROUTINE vector_assign_vec(vec1,vec2)

      TYPE(vector_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(:,:,:)=vec2%arr(:,:,:)
        IF (ASSOCIATED(vec1%arrh)) vec1%arrh(:,:,:,:)=vec2%arrh(:,:,:,:)
        IF (ASSOCIATED(vec1%arrv)) vec1%arrv(:,:,:,:)=vec2%arrv(:,:,:,:)
        IF (ASSOCIATED(vec1%arri)) vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp(:,:,:,:)=vec2%arrtmp(:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_assign_vec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assign_vec
!-----------------------------------------------------------------------
!     subprogram 7. vector_assignq_vec
!     set one vector structure equal to another. limit the quantity.
!-----------------------------------------------------------------------
      SUBROUTINE vector_assignq_vec(vec1,vec2,nqty)

      TYPE(vector_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: nqty

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(1:nqty,:,:)=vec2%arr(1:nqty,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &     vec1%arrh(1:nqty,:,:,:)=vec2%arrh(1:nqty,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &     vec1%arrv(1:nqty,:,:,:)=vec2%arrv(1:nqty,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &     vec1%arri(1:nqty,:,:,:)=vec2%arri(1:nqty,:,:,:)
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arrtmp))          &
     &     vec1%arri(1:nqty,:,:,:)=vec2%arri(1:nqty,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(1:nqty,:,:,:)=vec2%arri(1:nqty,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_assignq_vec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assignq_vec
!-----------------------------------------------------------------------
!     subprogram 8. vector_assignp_vec
!     set part of a complex 2D vector structure equal to part of a
!     complex vector structure. Assign
!       vec1(nq1:sz-1+nq1)=vec2(nq2:sz-1+nq2).
!     No checks on compatibility are performed.
!-----------------------------------------------------------------------
      SUBROUTINE vector_assignp_vec(vec1,vec2,nq1,nq2,sz)

      TYPE(vector_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: nq1,nq2,sz

      INTEGER(i4) :: ne1,ne2

      ne1=nq1+sz-1
      ne2=nq2+sz-1
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(nq1:ne1,:,:)=vec2%arr(nq2:ne2,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &    vec1%arrh(nq1:ne1,:,:,:)=vec2%arrh(nq2:ne2,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &    vec1%arrv(nq1:ne1,:,:,:)=vec2%arrv(nq2:ne2,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &    vec1%arri(nq1:ne1,:,:,:)=vec2%arri(nq2:ne2,:,:,:)
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arrtmp))          &
     &    vec1%arrtmp(nq1:ne1,:,:,:)=vec2%arrtmp(nq2:ne2,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(nq1:ne1,:,:,:)=vec2%arri(nq2:ne2,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_assign_part_vec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assignp_vec
!-----------------------------------------------------------------------
!     subprogram 9. vector_assignq_cvec2
!     set one 2D complex vector structure equal to real vector.
!-----------------------------------------------------------------------
      SUBROUTINE vector_assignq_cvec2(vec1,vec2,nqty)

      TYPE(vector_type),     INTENT(INOUT) :: vec1
      TYPE(cvector_2D_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: nqty

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(1:nqty,:,:)=vec2%arr(1:nqty,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &     vec1%arrh(1:nqty,:,:,:)=vec2%arrh(1:nqty,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &     vec1%arrv(1:nqty,:,:,:)=vec2%arrv(1:nqty,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &     vec1%arri(1:nqty,:,:,:)=vec2%arri(1:nqty,:,:,:)
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arrtmp))          &
     &     vec1%arrtmp(1:nqty,:,:,:)=vec2%arrtmp(1:nqty,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(1:nqty,:,:,:)=vec2%arri(1:nqty,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_assignq_cvec2: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assignq_cvec2
!-----------------------------------------------------------------------
!     subprogram 10. vector_assign_cvec.
!     set one real vector structure equal to part of a complex vector
!     structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_assign_cvec(vec1,vec2,r_i,fcomp,nqty)

      TYPE(vector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2
      CHARACTER(*), INTENT(IN) :: r_i
      INTEGER(i4), INTENT(IN) :: fcomp
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqty

!-----------------------------------------------------------------------
!     if the number of vector components is specified with nqty, limit
!     transfer.
!-----------------------------------------------------------------------
      IF (PRESENT(nqty)) THEN
        IF (ASSOCIATED(vec1%arr)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arr(1:nqty,:,:)=vec2%arr(1:nqty,:,:,fcomp)
            IF (ASSOCIATED(vec1%arrh))                                  &
     &        vec1%arrh(1:nqty,:,:,:)=vec2%arrh(1:nqty,:,:,:,fcomp)
            IF (ASSOCIATED(vec1%arrv))                                  &
     &        vec1%arrv(1:nqty,:,:,:)=vec2%arrv(1:nqty,:,:,:,fcomp)
            IF (ASSOCIATED(vec1%arri))                                  &
     &        vec1%arri(1:nqty,:,:,:)=vec2%arri(1:nqty,:,:,:,fcomp)
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arrtmp))      &
     &        vec1%arrtmp(1:nqty,:,:,:)=vec2%arrtmp(1:nqty,:,:,:,fcomp)
          CASE ('imag','IMAG')
            vec1%arr(1:nqty,:,:)=AIMAG(vec2%arr(1:nqty,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrh)) vec1%arrh(1:nqty,:,:,:)=         &
     &        AIMAG(vec2%arrh(1:nqty,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrv)) vec1%arrv(1:nqty,:,:,:)=         &
     &        AIMAG(vec2%arrv(1:nqty,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arri)) vec1%arri(1:nqty,:,:,:)=         &
     &        AIMAG(vec2%arri(1:nqty,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arrtmp))      &
     &        vec1%arri(1:nqty,:,:,:)=                                  &
     &        AIMAG(vec2%arri(1:nqty,:,:,:,fcomp))
          CASE DEFAULT
            CALL nim_stop                                               &
     &        ('Vector_assign_cvec: '//r_i//' flag not recognized.')
          END SELECT
!-----------------------------------------------------------------------
!       cell-centered data only.
!-----------------------------------------------------------------------
        ELSE IF (ASSOCIATED(vec1%arri)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arri(1:nqty,:,:,:)=vec2%arri(1:nqty,:,:,:,fcomp)
          CASE ('imag','IMAG')
            vec1%arri(1:nqty,:,:,:)=AIMAG(vec2%arri(1:nqty,:,:,:,fcomp))
          CASE DEFAULT
            CALL nim_stop                                               &
     &        ('Vector_assign_cvec: '//r_i//' flag not recognized.')
          END SELECT
        ELSE
          CALL nim_stop                                                 &
     &      ('Vector_assign_cvec: vector arrays not associated.')
        ENDIF
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!     the r12mi3 and i12r3 flags are used in several of the
!     nimrod management routines.  r12mi3 means transfer the real
!     1 & 2 vector components and minus the third imaginary  comp.
!     i12r3 means transfer the imaginary 1 & 2 and the real third.
!-----------------------------------------------------------------------
      ELSE
        IF (ASSOCIATED(vec1%arr)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arr=vec2%arr(:,:,:,fcomp)
            IF (ASSOCIATED(vec1%arrh))                                  &
     &        vec1%arrh=vec2%arrh(:,:,:,:,fcomp)
            IF (ASSOCIATED(vec1%arrv))                                  &
     &        vec1%arrv=vec2%arrv(:,:,:,:,fcomp)
            IF (ASSOCIATED(vec1%arri))                                  &
     &        vec1%arri=vec2%arri(:,:,:,:,fcomp)
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arrtmp))      &
     &        vec1%arrtmp=vec2%arrtmp(:,:,:,:,fcomp)
          CASE ('imag','IMAG')
            vec1%arr=AIMAG(vec2%arr(:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrh))                                  &
     &        vec1%arrh=AIMAG(vec2%arrh(:,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrv))                                  &
     &        vec1%arrv=AIMAG(vec2%arrv(:,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arri))                                  &
     &        vec1%arri=AIMAG(vec2%arri(:,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arrtmp))      &
     &        vec1%arrtmp=AIMAG(vec2%arrtmp(:,:,:,:,fcomp))
          CASE ('r12mi3')
            vec1%arr(1:2,:,:)=vec2%arr(1:2,:,:,fcomp)
            vec1%arr(3,:,:)=-AIMAG(vec2%arr(3,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrh)) THEN
              vec1%arrh(1:2,:,:,:)=vec2%arrh(1:2,:,:,:,fcomp)
              vec1%arrh(3,:,:,:)=-AIMAG(vec2%arrh(3,:,:,:,fcomp))
            ENDIF
            IF (ASSOCIATED(vec1%arrv)) THEN
              vec1%arrv(1:2,:,:,:)=vec2%arrv(1:2,:,:,:,fcomp)
              vec1%arrv(3,:,:,:)=-AIMAG(vec2%arrv(3,:,:,:,fcomp))
            ENDIF
            IF (ASSOCIATED(vec1%arri)) THEN
              vec1%arri(1:2,:,:,:)=vec2%arri(1:2,:,:,:,fcomp)
              vec1%arri(3,:,:,:)=-AIMAG(vec2%arri(3,:,:,:,fcomp))
            ENDIF
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arrtmp)) THEN
              vec1%arrtmp(1:2,:,:,:)=vec2%arrtmp(1:2,:,:,:,fcomp)
              vec1%arrtmp(3,:,:,:)=-AIMAG(vec2%arrtmp(3,:,:,:,fcomp))
            ENDIF
          CASE ('i12r3')
            vec1%arr(1:2,:,:)=AIMAG(vec2%arr(1:2,:,:,fcomp))
            vec1%arr(3,:,:)=vec2%arr(3,:,:,fcomp)
            IF (ASSOCIATED(vec1%arrh)) THEN
              vec1%arrh(1:2,:,:,:)=AIMAG(vec2%arrh(1:2,:,:,:,fcomp))
              vec1%arrh(3,:,:,:)=vec2%arrh(3,:,:,:,fcomp)
            ENDIF
            IF (ASSOCIATED(vec1%arrv)) THEN
              vec1%arrv(1:2,:,:,:)=AIMAG(vec2%arrv(1:2,:,:,:,fcomp))
              vec1%arrv(3,:,:,:)=vec2%arrv(3,:,:,:,fcomp)
            ENDIF
            IF (ASSOCIATED(vec1%arri)) THEN
              vec1%arri(1:2,:,:,:)=AIMAG(vec2%arri(1:2,:,:,:,fcomp))
              vec1%arri(3,:,:,:)=vec2%arri(3,:,:,:,fcomp)
            ENDIF
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arrtmp)) THEN
              vec1%arrtmp(1:2,:,:,:)=AIMAG(vec2%arrtmp(1:2,:,:,:,fcomp))
              vec1%arrtmp(3,:,:,:)=vec2%arrtmp(3,:,:,:,fcomp)
            ENDIF
          CASE DEFAULT
            CALL nim_stop                                               &
     &        ('Vector_assign_cvec: '//r_i//' flag not recognized.')
          END SELECT
!-----------------------------------------------------------------------
!       cell-centered data only.
!-----------------------------------------------------------------------
        ELSE IF (ASSOCIATED(vec1%arri)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arri=vec2%arri(:,:,:,:,fcomp)
          CASE ('imag','IMAG')
            vec1%arri=AIMAG(vec2%arri(:,:,:,:,fcomp))
          CASE ('r12mi3')
            vec1%arri(1:2,:,:,:)=vec2%arri(1:2,:,:,:,fcomp)
            vec1%arri(3,:,:,:)=-AIMAG(vec2%arri(3,:,:,:,fcomp))
          CASE ('i12r3')
            vec1%arri(1:2,:,:,:)=AIMAG(vec2%arri(1:2,:,:,:,fcomp))
            vec1%arri(3,:,:,:)=vec2%arri(3,:,:,:,fcomp)
          CASE DEFAULT
            CALL nim_stop                                               &
     &        ('Vector_assign_cvec: '//r_i//' flag not recognized.')
          END SELECT
        ELSE
          CALL nim_stop                                                 &
     &      ('Vector_assign_cvec: vector arrays not associated.')
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assign_cvec
!-----------------------------------------------------------------------
!     subprogram 11. vector_assign_vec3.
!     set a vector structure equal to a 3D vector structure at the
!     specified index of the periodic coordinate.
!-----------------------------------------------------------------------
      SUBROUTINE vector_assign_vec3(vec1,vec2,ip)

      TYPE(vector_type), INTENT(INOUT) :: vec1
      TYPE(vector_3D_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: ip
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!--------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(:,:,:)=vec2%arr(:,:,:,ip)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &    vec1%arrh(:,:,:,:)=vec2%arrh(:,:,:,:,ip)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &    vec1%arrv(:,:,:,:)=vec2%arrv(:,:,:,:,ip)
        IF (ASSOCIATED(vec1%arri))                                      &
     &    vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:,ip)
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arrtmp))          &
     &    vec1%arrtmp(:,:,:,:)=vec2%arrtmp(:,:,:,:,ip)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:,ip)
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_assign_vec3: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assign_vec3
!-----------------------------------------------------------------------
!     subprogram 12. vector_ptassign_bc.
!     make a pointer assignment of real vector data to bicube data.
!-----------------------------------------------------------------------
      SUBROUTINE vector_ptassign_bc(vec,bc)
      USE bicube

      TYPE(vector_type), INTENT(OUT) :: vec
      TYPE(bicube_type), INTENT(IN) :: bc

      IF (ASSOCIATED(bc%fs)) THEN
        vec%arr=>bc%fs
      ELSE
        NULLIFY(vec%arr)
      ENDIF
      NULLIFY(vec%arrh,vec%arrv,vec%arri,vec%arrtmp)

      RETURN
      END SUBROUTINE vector_ptassign_bc
!-----------------------------------------------------------------------
!     subprogram 13. vector_ptassign_laq2.
!     make a pointer assignment of real vector data to 2D lagrange
!     quadrilateral data.
!-----------------------------------------------------------------------
      SUBROUTINE vector_ptassign_laq2(vec,laq)
      USE lagr_quad_mod

      TYPE(vector_type), INTENT(OUT) :: vec
      TYPE(lagr_quad_2D_type), INTENT(IN), TARGET :: laq

      IF (ALLOCATED(laq%fs)) THEN
        vec%arr=>laq%fs
      ELSE
        NULLIFY(vec%arr)
      ENDIF

      IF (ALLOCATED(laq%fsh)) THEN
        vec%arrh=>laq%fsh
      ELSE
        NULLIFY(vec%arrh)
      ENDIF

      IF (ALLOCATED(laq%fsv)) THEN
        vec%arrv=>laq%fsv
      ELSE
        NULLIFY(vec%arrv)
      ENDIF

      IF (ALLOCATED(laq%fsi)) THEN
        vec%arri=>laq%fsi
      ELSE
        NULLIFY(vec%arri)
      ENDIF

      NULLIFY(vec%arrtmp)

      RETURN
      END SUBROUTINE vector_ptassign_laq2
!-----------------------------------------------------------------------
!     subprogram 13.1. vector_ptassign_modq2.
!     make a pointer assignment of real vector data to 2D modal
!     quadrilateral data.
!-----------------------------------------------------------------------
      SUBROUTINE vector_ptassign_modq2(vec,modq)
      USE modal_type_mod

      TYPE(vector_type), INTENT(OUT) :: vec
      TYPE(modal_quad_2D_type), INTENT(IN) :: modq

      IF (ASSOCIATED(modq%fsi)) THEN
        vec%arri=>modq%fsi
      ELSE
        NULLIFY(vec%arri)
      ENDIF

      NULLIFY(vec%arr,vec%arrh,vec%arrv,vec%arrtmp)

      RETURN
      END SUBROUTINE vector_ptassign_modq2
!-----------------------------------------------------------------------
!     subprogram 14. vector_ptassign_tl2.
!     make a pointer assignment of real vector data to 2D tri_linear
!     data.
!-----------------------------------------------------------------------
      SUBROUTINE vector_ptassign_tl2(vec,tl2)
      USE tri_linear

      TYPE(vector_type), INTENT(OUT) :: vec
      TYPE(tri_linear_2D_type), INTENT(IN) :: tl2

      IF (ASSOCIATED(tl2%fs)) THEN
        vec%arr=>tl2%fs
      ELSE
        NULLIFY(vec%arr)
      ENDIF
      NULLIFY(vec%arrh,vec%arrv,vec%arri)

      RETURN
      END SUBROUTINE vector_ptassign_tl2
!-----------------------------------------------------------------------
!     subprogram 15. vector_add_vec.
!     add one vector structure to another.
!-----------------------------------------------------------------------
      SUBROUTINE vector_add_vec(vec1,vec2,v1fac,v2fac)

      TYPE(vector_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      REAL(r8), INTENT(IN), OPTIONAL :: v1fac,v2fac

      REAL(r8) :: v1f,v2f
!-----------------------------------------------------------------------
!     set coefficients to input if used.
!-----------------------------------------------------------------------
      IF (PRESENT(v1fac)) THEN
        v1f=v1fac
      ELSE
        v1f=1
      ENDIF
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1
      ENDIF
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr=v1f*vec1%arr+v2f*vec2%arr
        IF (ASSOCIATED(vec1%arrh)) vec1%arrh=v1f*vec1%arrh+v2f*vec2%arrh
        IF (ASSOCIATED(vec1%arrv)) vec1%arrv=v1f*vec1%arrv+v2f*vec2%arrv
        IF (ASSOCIATED(vec1%arri)) vec1%arri=v1f*vec1%arri+v2f*vec2%arri
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp=v1f*vec1%arrtmp+v2f*vec2%arrtmp
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri=v1f*vec1%arri+v2f*vec2%arri
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_add_vec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_add_vec
!-----------------------------------------------------------------------
!     subprogram 15.1 vector_mult_vec.
!     multiply one vector structure to another.
!-----------------------------------------------------------------------
      SUBROUTINE vector_mult_vec(vec1,vec2,iq1in,iq2in)

      TYPE(vector_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN), OPTIONAL :: iq1in,iq2in

      INTEGER(i4) :: iq1,iq2
!-----------------------------------------------------------------------
!     set coefficients to input if used.
!-----------------------------------------------------------------------
      IF (PRESENT(iq1in)) THEN
        iq1=iq1in
      ELSE
        iq1=1
      ENDIF
      IF (PRESENT(iq2in)) THEN
        iq2=iq2in
      ELSE
        iq2=1
      ENDIF
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(iq1,:,:)=vec1%arr(iq1,:,:)*vec2%arr(iq2,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &    vec1%arrh(iq1,:,:,:)=vec1%arrh(iq1,:,:,:)*vec2%arrh(iq2,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &    vec1%arrv(iq1,:,:,:)=vec1%arrv(iq1,:,:,:)*vec2%arrv(iq2,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &    vec1%arri(iq1,:,:,:)=vec1%arri(iq1,:,:,:)*vec2%arri(iq2,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp(iq1,:,:,:)=                                       &
     &    vec1%arrtmp(iq1,:,:,:)*vec2%arrtmp(iq2,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(iq1,:,:,:)=vec1%arri(iq1,:,:,:)*vec2%arri(iq2,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_mult_vec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_mult_vec
!-----------------------------------------------------------------------
!     subprogram 15.2 vector_div_vec.
!     divide one vector structure to another.
!-----------------------------------------------------------------------
      SUBROUTINE vector_div_vec(vec1,vec2,iq1in,iq2in)

      TYPE(vector_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN), OPTIONAL :: iq1in,iq2in

      INTEGER(i4) :: iq1,iq2
!-----------------------------------------------------------------------
!     set coefficients to input if used.
!-----------------------------------------------------------------------
      IF (PRESENT(iq1in)) THEN
        iq1=iq1in
      ELSE
        iq1=1
      ENDIF
      IF (PRESENT(iq2in)) THEN
        iq2=iq2in
      ELSE
        iq2=1
      ENDIF
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(iq1,:,:)=vec1%arr(iq1,:,:)/vec2%arr(iq2,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &    vec1%arrh(iq1,:,:,:)=vec1%arrh(iq1,:,:,:)/vec2%arrh(iq2,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &    vec1%arrv(iq1,:,:,:)=vec1%arrv(iq1,:,:,:)/vec2%arrv(iq2,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &    vec1%arri(iq1,:,:,:)=vec1%arri(iq1,:,:,:)/vec2%arri(iq2,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp(iq1,:,:,:)=                                       &
     &    vec1%arrtmp(iq1,:,:,:)/vec2%arrtmp(iq2,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(iq1,:,:,:)=vec1%arri(iq1,:,:,:)/vec2%arri(iq2,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_div_vec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_div_vec
!-----------------------------------------------------------------------
!     subprogram 16. vector_mult_rsc.
!     multiply a vector structure by a real scalar.
!-----------------------------------------------------------------------
      SUBROUTINE vector_mult_rsc(vec,rsc)

      TYPE(vector_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rsc

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rsc*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rsc*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rsc*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=rsc*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rsc*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rsc*vec%arri
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_mult_rsc: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_mult_rsc
!-----------------------------------------------------------------------
!     subprogram 17. vector_mult_int.
!     multiply a vector structure by a real scalar.
!-----------------------------------------------------------------------
      SUBROUTINE vector_mult_int(vec,int)

      TYPE(vector_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=int*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int*vec%arri
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_mult_int: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_mult_int
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      END MODULE rvector_type_mod


!-----------------------------------------------------------------------
!     module for 3D-in-space arrays of vectors quantities.
!-----------------------------------------------------------------------
      MODULE cvector_type_mod
      USE vector_defn_type_mod
      IMPLICIT NONE

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 18. vector_ctype_alloc.
!     allocates space for a complex vector_type structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_ctype_alloc(cvt,poly_degree,mx,my,nqty,nfour,   &
     &                              nbt,nqt)
#ifdef OBJ_MEM_PROF
      USE memlog, only: memlogger
#endif

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,nfour,poly_degree
      INTEGER(i4), INTENT(IN), OPTIONAL :: nbt,nqt
      TYPE(cvector_type), INTENT(OUT) :: cvt

      INTEGER(i4) :: sz
!-----------------------------------------------------------------------
!     allocate space according to the basis functions needed.
!-----------------------------------------------------------------------
      SELECT CASE(poly_degree)
      CASE(:0)  !  piecewise continuous fields
        ALLOCATE(cvt%arri(nqty,(poly_degree-1)**2,mx,my,nfour))
        NULLIFY(cvt%arr,cvt%arrh,cvt%arrv)
      CASE(1) !  linear elements
        ALLOCATE(cvt%arr(nqty,0:mx,0:my,nfour))
        NULLIFY(cvt%arri,cvt%arrh,cvt%arrv)
      CASE(2:) !  higher-order elements
        ALLOCATE(cvt%arr(nqty,0:mx,0:my,nfour))
        ALLOCATE(cvt%arrh(nqty,poly_degree-1,1:mx,0:my,nfour))
        ALLOCATE(cvt%arrv(nqty,poly_degree-1,0:mx,1:my,nfour))
        ALLOCATE(cvt%arri(nqty,(poly_degree-1)**2,1:mx,1:my,nfour))
      END SELECT
!-----------------------------------------------------------------------
!     if the optional input for temporary arrays is present, allocate
!     the arrtmp arrays; otherwise, nullify them.  note that nbt is
!     the number of temporary bases, and nqt is the number of quantities
!     at each basis.
!-----------------------------------------------------------------------
      IF (PRESENT(nbt)) THEN
        ALLOCATE(cvt%arrtmp(nqt,nbt,1:mx,1:my,nfour))
      ELSE
        NULLIFY(cvt%arrtmp)
      ENDIF
!-----------------------------------------------------------------------
!     register this object.
!-----------------------------------------------------------------------
      NULLIFY(cvt%mem_id)
#ifdef OBJ_MEM_PROF
      SELECT CASE(poly_degree)
      CASE(0)
        sz=SIZEOF(cvt%arri)
      CASE(1)
        sz=SIZEOF(cvt%arr)
      CASE(2:)
        sz= SIZEOF(cvt%arri)+SIZEOF(cvt%arr)                            &
     &     +SIZEOF(cvt%arrh)+SIZEOF(cvt%arrv)
      END SELECT
      CALL memlogger%update(cvt%mem_id,'cvector_3D','unknown',sz)
#endif
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_ctype_alloc
!-----------------------------------------------------------------------
!     subprogram 19. vector_ctype_dealloc.
!     deallocates space for a complex vector_type structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_ctype_dealloc(cvt)
#ifdef OBJ_MEM_PROF
      USE memlog, only: memlogger
#endif

      TYPE(cvector_type), INTENT(INOUT) :: cvt

      IF (ASSOCIATED(cvt%arr)) THEN
        DEALLOCATE(cvt%arr)
        NULLIFY(cvt%arr)
      ENDIF
      IF (ASSOCIATED(cvt%arrh)) THEN
        DEALLOCATE(cvt%arrh)
        NULLIFY(cvt%arrh)
      ENDIF
      IF (ASSOCIATED(cvt%arrv)) THEN
        DEALLOCATE(cvt%arrv)
        NULLIFY(cvt%arrv)
      ENDIF
      IF (ASSOCIATED(cvt%arri)) THEN
        DEALLOCATE(cvt%arri)
        NULLIFY(cvt%arri)
      ENDIF
      IF (ASSOCIATED(cvt%arrtmp)) THEN
        DEALLOCATE(cvt%arrtmp)
        NULLIFY(cvt%arrtmp)
      ENDIF
!-----------------------------------------------------------------------
!     unregister this object.
!-----------------------------------------------------------------------
#ifdef OBJ_MEM_PROF
      CALL memlogger%update(cvt%mem_id,' ',' ',0,resize=.TRUE.)
#endif
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_ctype_dealloc
!-----------------------------------------------------------------------
!     subprogram 20. cvector_assign_rsc.
!     assign a real scalar value to a complex vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_assign_rsc(vec,rscalar)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rscalar

!-----------------------------------------------------------------------
!     if the grid vertex-centered array is allocated, treat as a
!     standard element.  If, not the structure represents piecewise
!     constant.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=rscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rscalar
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_assign_rsc: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assign_rsc
!-----------------------------------------------------------------------
!     subprogram 20.1 cvector_assign_rsc_nq.
!     assign a real scalar value to a complex vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_assign_rsc_nq(vec,rscalar,i1,n)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rscalar
      INTEGER(i4), INTENT(IN) :: i1,n

      INTEGER(i4) :: i11
      i11 = i1 + n - 1
!-----------------------------------------------------------------------
!     if the grid vertex-centered array is allocated, treat as a
!     standard element.  If, not the structure represents piecewise
!     constant.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr(i1:i11,:,:,:)=rscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh(i1:i11,:,:,:,:)=rscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv(i1:i11,:,:,:,:)=rscalar
        IF (ASSOCIATED(vec%arri)) vec%arri(i1:i11,:,:,:,:)=rscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp(i1:i11,:,:,:,:)=rscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri(i1:i11,:,:,:,:)=rscalar
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_assign_rsc_nq: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assign_rsc_nq
!-----------------------------------------------------------------------
!     subprogram 21. cvector_assign_csc.
!     assign a complex scalar value to a complex vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_assign_csc(vec,cscalar)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      COMPLEX(r8), INTENT(IN) :: cscalar

!-----------------------------------------------------------------------
!     if the grid vertex-centered array is allocated, treat as a
!     standard element.  If, not the structure represents piecewise
!     constant.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=cscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=cscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=cscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=cscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=cscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=cscalar
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_assign_csc: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assign_csc
!-----------------------------------------------------------------------
!     subprogram 22. cvector_assign_int.
!     assign a integer value to a complex vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_assign_int(vec,int)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

!-----------------------------------------------------------------------
!     if the grid vertex-centered array is allocated, treat as a
!     standard element.  If, not the structure represents piecewise
!     constant.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int
        IF (ASSOCIATED(vec%arri)) vec%arri=int
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_assign_int: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assign_int
!-----------------------------------------------------------------------
!     subprogram 23. cvector_assign_cvec.
!     set one complex vector structure equal to another.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_assign_cvec(vec1,vec2)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(:,:,:,:)=vec2%arr(:,:,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &    vec1%arrh(:,:,:,:,:)=vec2%arrh(:,:,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &    vec1%arrv(:,:,:,:,:)=vec2%arrv(:,:,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &    vec1%arri(:,:,:,:,:)=vec2%arri(:,:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp(:,:,:,:,:)=vec2%arrtmp(:,:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(:,:,:,:,:)=vec2%arri(:,:,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_assign_cvec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assign_cvec
!-----------------------------------------------------------------------
!     subprogram 24. cvector_assignq_cvec.
!     set a quantity range of one complex vector structure equal to
!     that of another.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_assignq_cvec(vec1,vec2,nqty)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: nqty

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(1:nqty,:,:,:)=vec2%arr(1:nqty,:,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &    vec1%arrh(1:nqty,:,:,:,:)=vec2%arrh(1:nqty,:,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &    vec1%arrv(1:nqty,:,:,:,:)=vec2%arrv(1:nqty,:,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &    vec1%arri(1:nqty,:,:,:,:)=vec2%arri(1:nqty,:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp(1:nqty,:,:,:,:)=vec2%arrtmp(1:nqty,:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(1:nqty,:,:,:,:)=vec2%arri(1:nqty,:,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_assignq_cvec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assignq_cvec
!----------------------------------------------------------------------
!     subprogram 25. cvector_assignp_cvec.
!     set part of the first complex vector structure equal to part of
!     second.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_assignp_cvec(vec1,vec2,nq1,nq2,sz)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: nq1,nq2,sz

      INTEGER(i4) :: ne1,ne2

      ne1=nq1+sz-1
      ne2=nq2+sz-1
!-----------------------------------------------------------------------
!     copy only one of the vector components.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(nq1:ne1,:,:,:)=vec2%arr(nq2:ne2,:,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &    vec1%arrh(nq1:ne1,:,:,:,:)=vec2%arrh(nq2:ne2,:,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &    vec1%arrv(nq1:ne1,:,:,:,:)=vec2%arrv(nq2:ne2,:,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &    vec1%arri(nq1:ne1,:,:,:,:)=vec2%arri(nq2:ne2,:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp(nq1:ne1,:,:,:,:)=vec2%arrtmp(nq2:ne2,:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(nq1:ne1,:,:,:,:)=vec2%arri(nq2:ne2,:,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_assignp_cvec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assignp_cvec
!-----------------------------------------------------------------------
!     subprogram 26. cvector_assign_vec.
!     set one component of a complex vector structure equal to a real
!     vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_assign_vec(vec1,vec2,r_i,fcomp,nqty)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      CHARACTER(*), INTENT(IN) :: r_i
      INTEGER(i4), INTENT(IN) :: fcomp
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqty

!-----------------------------------------------------------------------
!     if the number of vector components is specified with nqty, limit
!     transfer.
!-----------------------------------------------------------------------
      IF (PRESENT(nqty)) THEN
        IF (ASSOCIATED(vec1%arr)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arr(1:nqty,:,:,fcomp)=vec2%arr(1:nqty,:,:)             &
     &        +(0,1)*AIMAG(vec1%arr(1:nqty,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrh))                                  &
     &        vec1%arrh(1:nqty,:,:,:,fcomp)=vec2%arrh(1:nqty,:,:,:)     &
     &          +(0,1)*AIMAG(vec1%arrh(1:nqty,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrv))                                  &
     &        vec1%arrv(1:nqty,:,:,:,fcomp)=vec2%arrv(1:nqty,:,:,:)     &
     &          +(0,1)*AIMAG(vec1%arrv(1:nqty,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arri))                                  &
     &        vec1%arri(1:nqty,:,:,:,fcomp)=vec2%arri(1:nqty,:,:,:)     &
     &          +(0,1)*AIMAG(vec1%arri(1:nqty,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))    &
     &        vec1%arrtmp(1:nqty,:,:,:,fcomp)=vec2%arrtmp(1:nqty,:,:,:) &
     &          +(0,1)*AIMAG(vec1%arrtmp(1:nqty,:,:,:,fcomp))
          CASE ('imag','IMAG')
            vec1%arr(1:nqty,:,:,fcomp)=(0,1)*vec2%arr(1:nqty,:,:)       &
     &        +REAL(vec1%arr(1:nqty,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arrh))                                  &
     &        vec1%arrh(1:nqty,:,:,:,fcomp)=                            &
     &          (0,1)*vec2%arrh(1:nqty,:,:,:)                           &
     &          +REAL(vec1%arrh(1:nqty,:,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arrv))                                  &
     &        vec1%arrv(1:nqty,:,:,:,fcomp)=                            &
     &          (0,1)*vec2%arrv(1:nqty,:,:,:)                           &
     &          +REAL(vec1%arrv(1:nqty,:,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arri))                                  &
     &        vec1%arri(1:nqty,:,:,:,fcomp)=                            &
     &          (0,1)*vec2%arri(1:nqty,:,:,:)                           &
     &          +REAL(vec1%arri(1:nqty,:,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))    &
     &        vec1%arrtmp(1:nqty,:,:,:,fcomp)=                          &
     &          (0,1)*vec2%arrtmp(1:nqty,:,:,:)                         &
     &          +REAL(vec1%arrtmp(1:nqty,:,:,:,fcomp),r8)
          CASE DEFAULT
            CALL nim_stop                                               &
     &        ('Cvector_assign_vec: '//r_i//' flag not recognized.')
          END SELECT
!-----------------------------------------------------------------------
!       cell-centered data only.
!-----------------------------------------------------------------------
        ELSE IF (ASSOCIATED(vec1%arri)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arri(1:nqty,:,:,:,fcomp)=vec2%arri(1:nqty,:,:,:)       &
     &        +(0,1)*AIMAG(vec1%arri(1:nqty,:,:,:,fcomp))
          CASE ('imag','IMAG')
            vec1%arri(1:nqty,:,:,:,fcomp)=(0,1)*vec2%arri(1:nqty,:,:,:) &
     &        +REAL(vec1%arri(1:nqty,:,:,:,fcomp),r8)
          CASE DEFAULT
            CALL nim_stop                                               &
     &        ('Cvector_assign_vec: '//r_i//' flag not recognized.')
          END SELECT
        ELSE
          CALL nim_stop                                                 &
     &      ('Cvector_assign_vec: vector arrays not associated.')
        ENDIF
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!     the r12mi3 and i12r3 flags are used in several of the
!     nimrod management routines.  r12mi3 means transfer the real
!     1 & 2 vector components and minus the third imaginary  comp.
!     i12r3 means transfer the imaginary 1 & 2 and the real third.
!-----------------------------------------------------------------------
      ELSE
        IF (ASSOCIATED(vec1%arr)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arr(:,:,:,fcomp)=vec2%arr                              &
     &               +(0,1)*AIMAG(vec1%arr(:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrh))                                  &
     &        vec1%arrh(:,:,:,:,fcomp)=vec2%arrh                        &
     &                    +(0,1)*AIMAG(vec1%arrh(:,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrv))                                  &
     &        vec1%arrv(:,:,:,:,fcomp)=vec2%arrv                        &
     &                    +(0,1)*AIMAG(vec1%arrv(:,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arri))                                  &
     &        vec1%arri(:,:,:,:,fcomp)=vec2%arri                        &
     &                    +(0,1)*AIMAG(vec1%arri(:,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))    &
     &        vec1%arrtmp(:,:,:,:,fcomp)=vec2%arrtmp                    &
     &                    +(0,1)*AIMAG(vec1%arrtmp(:,:,:,:,fcomp))
          CASE ('imag','IMAG')
            vec1%arr(:,:,:,fcomp)=(0,1)*vec2%arr                        &
     &                            +REAL(vec1%arr(:,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arrh))                                  &
     &        vec1%arrh(:,:,:,:,fcomp)=(0,1)*vec2%arrh                  &
     &                                +REAL(vec1%arrh(:,:,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arrv))                                  &
     &        vec1%arrv(:,:,:,:,fcomp)=(0,1)*vec2%arrv                  &
     &                                +REAL(vec1%arrv(:,:,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arri))                                  &
     &        vec1%arri(:,:,:,:,fcomp)=(0,1)*vec2%arri                  &
     &                                +REAL(vec1%arri(:,:,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))    &
     &        vec1%arrtmp(:,:,:,:,fcomp)=(0,1)*vec2%arrtmp              &
     &                              +REAL(vec1%arrtmp(:,:,:,:,fcomp),r8)
          CASE ('r12mi3')
            vec1%arr(1:2,:,:,fcomp)=vec2%arr(1:2,:,:)                   &
     &               +(0,1)*AIMAG(vec1%arr(1:2,:,:,fcomp))
            vec1%arr(3,:,:,fcomp)=-(0,1)*vec2%arr(3,:,:)                &
     &               +REAL(vec1%arr(3,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arrh)) THEN
              vec1%arrh(1:2,:,:,:,fcomp)=vec2%arrh(1:2,:,:,:)           &
     &                 +(0,1)*AIMAG(vec1%arrh(1:2,:,:,:,fcomp))
              vec1%arrh(3,:,:,:,fcomp)=-(0,1)*vec2%arrh(3,:,:,:)        &
     &                 +REAL(vec1%arrh(3,:,:,:,fcomp),r8)
            ENDIF
            IF (ASSOCIATED(vec1%arrv)) THEN
              vec1%arrv(1:2,:,:,:,fcomp)=vec2%arrv(1:2,:,:,:)           &
     &                 +(0,1)*AIMAG(vec1%arrv(1:2,:,:,:,fcomp))
              vec1%arrv(3,:,:,:,fcomp)=-(0,1)*vec2%arrv(3,:,:,:)        &
     &                 +REAL(vec1%arrv(3,:,:,:,fcomp),r8)
            ENDIF
            IF (ASSOCIATED(vec1%arri)) THEN
              vec1%arri(1:2,:,:,:,fcomp)=vec2%arri(1:2,:,:,:)           &
     &                 +(0,1)*AIMAG(vec1%arri(1:2,:,:,:,fcomp))
              vec1%arri(3,:,:,:,fcomp)=-(0,1)*vec2%arri(3,:,:,:)        &
     &                 +REAL(vec1%arri(3,:,:,:,fcomp),r8)
            ENDIF
            IF(ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp)) THEN
              vec1%arrtmp(1:2,:,:,:,fcomp)=vec2%arrtmp(1:2,:,:,:)       &
     &                 +(0,1)*AIMAG(vec1%arrtmp(1:2,:,:,:,fcomp))
              vec1%arrtmp(3,:,:,:,fcomp)=-(0,1)*vec2%arrtmp(3,:,:,:)    &
     &                 +REAL(vec1%arrtmp(3,:,:,:,fcomp),r8)
            ENDIF
          CASE ('i12r3')
            vec1%arr(1:2,:,:,fcomp)=(0,1)*vec2%arr(1:2,:,:)             &
     &               +REAL(vec1%arr(1:2,:,:,fcomp),r8)
            vec1%arr(3,:,:,fcomp)=vec2%arr(3,:,:)                       &
     &               +(0,1)*AIMAG(vec1%arr(3,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrh)) THEN
              vec1%arrh(1:2,:,:,:,fcomp)=(0,1)*vec2%arrh(1:2,:,:,:)     &
     &                  +REAL(vec1%arrh(1:2,:,:,:,fcomp),r8)
              vec1%arrh(3,:,:,:,fcomp)=vec2%arrh(3,:,:,:)               &
     &                  +(0,1)*AIMAG(vec1%arrh(3,:,:,:,fcomp))
            ENDIF
            IF (ASSOCIATED(vec1%arrv)) THEN
              vec1%arrv(1:2,:,:,:,fcomp)=(0,1)*vec2%arrv(1:2,:,:,:)     &
     &                  +REAL(vec1%arrv(1:2,:,:,:,fcomp),r8)
              vec1%arrv(3,:,:,:,fcomp)=vec2%arrv(3,:,:,:)               &
     &                  +(0,1)*AIMAG(vec1%arrv(3,:,:,:,fcomp))
            ENDIF
            IF (ASSOCIATED(vec1%arri)) THEN
              vec1%arri(1:2,:,:,:,fcomp)=(0,1)*vec2%arri(1:2,:,:,:)     &
     &                  +REAL(vec1%arri(1:2,:,:,:,fcomp),r8)
              vec1%arri(3,:,:,:,fcomp)=vec2%arri(3,:,:,:)               &
     &                  +(0,1)*AIMAG(vec1%arri(3,:,:,:,fcomp))
            ENDIF
            IF(ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp)) THEN
              vec1%arrtmp(1:2,:,:,:,fcomp)=(0,1)*vec2%arrtmp(1:2,:,:,:) &
     &                  +REAL(vec1%arrtmp(1:2,:,:,:,fcomp),r8)
              vec1%arrtmp(3,:,:,:,fcomp)=vec2%arrtmp(3,:,:,:)           &
     &                  +(0,1)*AIMAG(vec1%arrtmp(3,:,:,:,fcomp))
            ENDIF
          CASE DEFAULT
            CALL nim_stop                                               &
     &        ('Cvector_assign_vec: '//r_i//' flag not recognized.')
          END SELECT
!-----------------------------------------------------------------------
!       cell-centered data only.
!-----------------------------------------------------------------------
        ELSE IF (ASSOCIATED(vec1%arri)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arri(:,:,:,:,fcomp)=vec2%arri                          &
     &                  +(0,1)*AIMAG(vec1%arri(:,:,:,:,fcomp))
          CASE ('imag','IMAG')
            vec1%arri(:,:,:,:,fcomp)=(0,1)*vec2%arri                    &
     &                               +REAL(vec1%arri(:,:,:,:,fcomp),r8)
          CASE ('r12mi3')
            vec1%arri(1:2,:,:,:,fcomp)=vec2%arri(1:2,:,:,:)             &
     &               +(0,1)*AIMAG(vec1%arri(1:2,:,:,:,fcomp))
            vec1%arri(3,:,:,:,fcomp)=-(0,1)*vec2%arri(3,:,:,:)          &
     &               +REAL(vec1%arri(3,:,:,:,fcomp),r8)
          CASE ('i12r3')
            vec1%arri(1:2,:,:,:,fcomp)=(0,1)*vec2%arri(1:2,:,:,:)       &
     &                +REAL(vec1%arri(1:2,:,:,:,fcomp),r8)
            vec1%arri(3,:,:,:,fcomp)=vec2%arri(3,:,:,:)                 &
     &                +(0,1)*AIMAG(vec1%arri(3,:,:,:,fcomp))
          CASE DEFAULT
            CALL nim_stop                                               &
     &        ('Cvector_assign_vec: '//r_i//' flag not recognized.')
          END SELECT
        ELSE
          CALL nim_stop                                                 &
     &      ('Cvector_assign_vec: vector arrays not associated.')
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assign_vec
!-----------------------------------------------------------------------
!     subprogram 27. cvector_assign_cvec2.
!     set one component of a complex vector structure equal to a complex
!     2D vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_assign_cvec2(vec1,vec2,fcomp,nqty)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_2D_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: fcomp
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqty

!-----------------------------------------------------------------------
!     if the number of vector components is specified with nqty, limit
!     transfer.
!-----------------------------------------------------------------------
      IF (PRESENT(nqty)) THEN
        IF (ASSOCIATED(vec1%arr)) THEN
          vec1%arr(1:nqty,:,:,fcomp)=vec2%arr(1:nqty,:,:)
          IF (ASSOCIATED(vec1%arrh))                                    &
     &      vec1%arrh(1:nqty,:,:,:,fcomp)=vec2%arrh(1:nqty,:,:,:)
          IF (ASSOCIATED(vec1%arrv))                                    &
     &      vec1%arrv(1:nqty,:,:,:,fcomp)=vec2%arrv(1:nqty,:,:,:)
          IF (ASSOCIATED(vec1%arri))                                    &
     &      vec1%arri(1:nqty,:,:,:,fcomp)=vec2%arri(1:nqty,:,:,:)
          IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))      &
     &      vec1%arrtmp(1:nqty,:,:,:,fcomp)=vec2%arrtmp(1:nqty,:,:,:)
        ELSE IF (ASSOCIATED(vec1%arri)) THEN
          vec1%arri(1:nqty,:,:,:,fcomp)=vec2%arri(1:nqty,:,:,:)
        ELSE
          CALL nim_stop                                                 &
     &      ('Cvector_assign_cvec2: vector arrays not associated.')
        ENDIF
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      ELSE
        IF (ASSOCIATED(vec1%arr)) THEN
          vec1%arr(:,:,:,fcomp)=vec2%arr
          IF (ASSOCIATED(vec1%arrh)) vec1%arrh(:,:,:,:,fcomp)=vec2%arrh
          IF (ASSOCIATED(vec1%arrv)) vec1%arrv(:,:,:,:,fcomp)=vec2%arrv
          IF (ASSOCIATED(vec1%arri)) vec1%arri(:,:,:,:,fcomp)=vec2%arri
          IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))      &
     &      vec1%arrtmp(:,:,:,:,fcomp)=vec2%arrtmp 
        ELSE IF (ASSOCIATED(vec1%arri)) THEN
          vec1%arri(:,:,:,:,fcomp)=vec2%arri
        ELSE
          CALL nim_stop                                                 &
     &      ('Cvector_assign_cvec2: vector arrays not associated.')
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assign_cvec2
!-----------------------------------------------------------------------
!     subprogram 28. cvector_assignp_cvec2.
!     set part of a complex vector structure equal to part of a
!     complex 2D vector structure
!-----------------------------------------------------------------------
      SUBROUTINE cvector_assignp_cvec2(vec1,vec2,fcomp,nq1,nq2,sz)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_2D_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: fcomp,nq1,nq2,sz

      INTEGER(i4) :: ne1,ne2

      ne1=nq1+sz-1
      ne2=nq2+sz-1
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(nq1:ne1,:,:,fcomp)=vec2%arr(nq2:ne2,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &    vec1%arrh(nq1:ne1,:,:,:,fcomp)=vec2%arrh(nq2:ne2,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &    vec1%arrv(nq1:ne1,:,:,:,fcomp)=vec2%arrv(nq2:ne2,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &    vec1%arri(nq1:ne1,:,:,:,fcomp)=vec2%arri(nq2:ne2,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp(nq1:ne1,:,:,:,fcomp)=vec2%arrtmp(nq2:ne2,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(nq1:ne1,:,:,:,fcomp)=vec2%arri(nq2:ne2,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_assignp_cvec2: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assignp_cvec2
!-----------------------------------------------------------------------
!     subprogram 29. cvector_ptassign_laq.
!     make a pointer assignment of complex vector data to lagrange
!     quadrilateral data.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_ptassign_laq(vec,laq)
      USE lagr_quad_mod

      TYPE(cvector_type), INTENT(OUT) :: vec
      TYPE(lagr_quad_type), INTENT(IN), TARGET :: laq

      IF (ALLOCATED(laq%fs)) THEN
         vec%arr=>laq%fs
      ELSE
        NULLIFY(vec%arr)
      ENDIF

      IF (ALLOCATED(laq%fsh)) THEN
        vec%arrh=>laq%fsh
      ELSE
        NULLIFY(vec%arrh)
      ENDIF

      IF (ALLOCATED(laq%fsv)) THEN
        vec%arrv=>laq%fsv
      ELSE
        NULLIFY(vec%arrv)
      ENDIF

      IF (ALLOCATED(laq%fsi)) THEN
        vec%arri=>laq%fsi
      ELSE
        NULLIFY(vec%arri)
      ENDIF

      NULLIFY(vec%arrtmp)

      RETURN
      END SUBROUTINE cvector_ptassign_laq
!-----------------------------------------------------------------------
!     subprogram 29.1. cvector_ptassign_modq.
!     make a pointer assignment of complex vector data to modal
!     quadrilateral data.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_ptassign_modq(vec,modq)
      USE modal_type_mod

      TYPE(cvector_type), INTENT(OUT) :: vec
      TYPE(modal_quad_type), INTENT(IN) :: modq

      IF (ASSOCIATED(modq%fsi)) THEN
        vec%arri=>modq%fsi
      ELSE
        NULLIFY(vec%arri)
      ENDIF

      NULLIFY(vec%arr,vec%arrh,vec%arrv,vec%arrtmp)

      RETURN
      END SUBROUTINE cvector_ptassign_modq
!-----------------------------------------------------------------------
!     subprogram 30. cvector_ptassign_tl.
!     make a pointer assignment of complex vector data to 3D tri_linear
!     data.
!-PRE this will eventually need the additional arrays.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_ptassign_tl(vec,tl)
      USE tri_linear

      TYPE(cvector_type), INTENT(OUT) :: vec
      TYPE(tri_linear_type), INTENT(IN) :: tl

      IF (ASSOCIATED(tl%fs)) THEN
        vec%arr=>tl%fs
      ELSE
        NULLIFY(vec%arr)
      ENDIF
      NULLIFY(vec%arrh,vec%arrv,vec%arri,vec%arrtmp)

      RETURN
      END SUBROUTINE cvector_ptassign_tl
!-----------------------------------------------------------------------
!     subprogram 31. cvector_add_cvec.
!     add one complex vector structure to another.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_add_cvec(vec1,vec2,v1fac,v2fac)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2
      REAL(r8), INTENT(IN), OPTIONAL :: v1fac,v2fac

      REAL(r8) :: v1f,v2f
!-----------------------------------------------------------------------
!     set coefficients to input if used.
!-----------------------------------------------------------------------
      IF (PRESENT(v1fac)) THEN
        v1f=v1fac
      ELSE
        v1f=1
      ENDIF
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1
      ENDIF
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr=v1f*vec1%arr+v2f*vec2%arr
        IF (ASSOCIATED(vec1%arrh)) vec1%arrh=v1f*vec1%arrh+v2f*vec2%arrh
        IF (ASSOCIATED(vec1%arrv)) vec1%arrv=v1f*vec1%arrv+v2f*vec2%arrv
        IF (ASSOCIATED(vec1%arri)) vec1%arri=v1f*vec1%arri+v2f*vec2%arri
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp=v1f*vec1%arrtmp+v2f*vec2%arrtmp
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri=v1f*vec1%arri+v2f*vec2%arri
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_add_cvec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_add_cvec
!-----------------------------------------------------------------------
!     subprogram 32. cvector_add_iq.
!     add i1 component of a complex vector structure to i2
!     component of another complex vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_add_iq(vec1,vec2,i1,i2,v1fac,v2fac)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: i1,i2
      REAL(r8), INTENT(IN), OPTIONAL :: v1fac,v2fac

      REAL(r8) :: v1f,v2f
!-----------------------------------------------------------------------
!     set coefficients to input if used.
!-----------------------------------------------------------------------
      IF (PRESENT(v1fac)) THEN
        v1f=v1fac
      ELSE
        v1f=1
      ENDIF
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1
      ENDIF

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(i1,:,:,:)=v1f*vec1%arr(i1,:,:,:)+v2f*vec2%arr(i2,:,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &    vec1%arrh(i1,:,:,:,:)=v1f*vec1%arrh(i1,:,:,:,:)+              &
     &                          v2f*vec2%arrh(i2,:,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &    vec1%arrv(i1,:,:,:,:)=v1f*vec1%arrv(i1,:,:,:,:)+              &
     &                          v2f*vec2%arrv(i2,:,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &    vec1%arri(i1,:,:,:,:)=v1f*vec1%arri(i1,:,:,:,:)+              &
     &                          v2f*vec2%arri(i2,:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp(i1,:,:,:,:)=v1f*vec1%arrtmp(i1,:,:,:,:)+          &
     &                            v2f*vec2%arrtmp(i2,:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(i1,:,:,:,:)=v1f*vec1%arri(i1,:,:,:,:)+                &
     &                        v2f*vec2%arri(i2,:,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('cvector_add_iq: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_add_iq
!-----------------------------------------------------------------------
!     subprogram 32.1 cvector_add_nq.
!     add n components of two complex vector structures together
!     starting with indicies i1 and i2 respectively
!-----------------------------------------------------------------------
      SUBROUTINE cvector_add_nq(vec1,vec2,i1,i2,n)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: i1,i2,n

      INTEGER(i4) :: i11,i22
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      i11 = i1 + n - 1
      i22 = i2 + n - 1
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(i1:i11,:,:,:)=vec1%arr(i1:i11,:,:,:)+                  &
     &                         vec2%arr(i2:i22,:,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &    vec1%arrh(i1:i11,:,:,:,:)=vec1%arrh(i1:i11,:,:,:,:)+          &
     &                              vec2%arrh(i2:i22,:,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &    vec1%arrv(i1:i11,:,:,:,:)=vec1%arrv(i1:i11,:,:,:,:)+          &
     &                              vec2%arrv(i2:i22,:,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &    vec1%arri(i1:i11,:,:,:,:)=vec1%arri(i1:i11,:,:,:,:)+          &
     &                              vec2%arri(i2:i22,:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp(i1:i11,:,:,:,:)=vec1%arrtmp(i1:i11,:,:,:,:)+      &
     &                                vec2%arrtmp(i2:i22,:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(i1:i11,:,:,:,:)=vec1%arri(i1:i11,:,:,:,:)+            &
     &                            vec2%arri(i2:i22,:,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('cvector_add_nq: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_add_nq
!-----------------------------------------------------------------------
!     subprogram 33. cvector_addc_cvec.
!     add one complex vector structure to another a complex scalar
!     with complex coefficients.  here, only v2fac is optional, so that
!     there is no confusion with cvector_add_cvec in the module
!     procedure definition.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_addc_cvec(vec1,vec2,v1f,v2fac)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2
      COMPLEX(r8), INTENT(IN) :: v1f
      COMPLEX(r8), INTENT(IN), OPTIONAL :: v2fac

      COMPLEX(r8) :: v2f
!-----------------------------------------------------------------------
!     set coefficient to input if used.
!-----------------------------------------------------------------------
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1._r8
      ENDIF
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr=v1f*vec1%arr+v2f*vec2%arr
        IF (ASSOCIATED(vec1%arrh)) vec1%arrh=v1f*vec1%arrh+v2f*vec2%arrh
        IF (ASSOCIATED(vec1%arrv)) vec1%arrv=v1f*vec1%arrv+v2f*vec2%arrv
        IF (ASSOCIATED(vec1%arri)) vec1%arri=v1f*vec1%arri+v2f*vec2%arri
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp=v1f*vec1%arrtmp+v2f*vec2%arrtmp
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri=v1f*vec1%arri+v2f*vec2%arri
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_addc_cvec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_addc_cvec
!-----------------------------------------------------------------------
!     subprogram 34. cvector_mult_rsc.
!     multiply a complex vector structure by a real scalar.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_mult_rsc(vec,rsc)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rsc

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rsc*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rsc*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rsc*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=rsc*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rsc*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rsc*vec%arri
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_mult_rsc: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_mult_rsc
!-----------------------------------------------------------------------
!     subprogram 35. cvector_mult_csc.
!     multiply a complex vector structure by a complex scalar.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_mult_csc(vec,csc)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      COMPLEX(r8), INTENT(IN) :: csc

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=csc*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=csc*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=csc*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=csc*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=csc*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=csc*vec%arri
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_mult_csc: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_mult_csc
!-----------------------------------------------------------------------
!     subprogram 36. cvector_mult_int.
!     multiply a complex vector structure by a real scalar.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_mult_int(vec,int)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=int*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int*vec%arri
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_mult_int: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_mult_int
!-----------------------------------------------------------------------
!     subprogram 37. cvector_real_comp.
!     set one component of a complex vector structure its real part.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_real_comp(vec1,fcomp)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      INTEGER(i4), INTENT(IN) :: fcomp

!-----------------------------------------------------------------------
!     use r8 in the intrinsic calls to maintain accuracy.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr))                                         &
     &  vec1%arr(:,:,:,fcomp)=REAL(vec1%arr(:,:,:,fcomp),r8)
      IF (ASSOCIATED(vec1%arrh))                                        &
     &  vec1%arrh(:,:,:,:,fcomp)=REAL(vec1%arrh(:,:,:,:,fcomp),r8)
      IF (ASSOCIATED(vec1%arrv))                                        &
     &  vec1%arrv(:,:,:,:,fcomp)=REAL(vec1%arrv(:,:,:,:,fcomp),r8)
      IF (ASSOCIATED(vec1%arri))                                        &
     &  vec1%arri(:,:,:,:,fcomp)=REAL(vec1%arri(:,:,:,:,fcomp),r8)
      IF (ASSOCIATED(vec1%arrtmp))                                      &
     &  vec1%arrtmp(:,:,:,:,fcomp)=REAL(vec1%arrtmp(:,:,:,:,fcomp),r8)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_real_comp
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      END MODULE cvector_type_mod


!-----------------------------------------------------------------------
!     module for 2D-in-space arrays of complex vectors quantities.
!-----------------------------------------------------------------------
      MODULE cvector_2D_type_mod
      USE local
      USE vector_defn_type_mod
      IMPLICIT NONE

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 38. vector_2D_ctype_alloc.
!     allocates space for a 2D complex vector_type structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_2D_ctype_alloc(cvt,poly_degree,mx,my,nqty,      &
     &                                 nbt,nqt)
#ifdef OBJ_MEM_PROF
      USE memlog, only: memlogger
#endif

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,poly_degree
      INTEGER(i4), INTENT(IN), OPTIONAL :: nbt,nqt
      TYPE(cvector_2D_type), INTENT(OUT) :: cvt

      INTEGER(i4) :: sz
!-----------------------------------------------------------------------
!     allocate space according to the basis functions needed.
!-----------------------------------------------------------------------
      SELECT CASE(poly_degree)
      CASE(:0)  !  piecewise continuous fields
        ALLOCATE(cvt%arri(nqty,(poly_degree-1)**2,mx,my))
        NULLIFY(cvt%arr,cvt%arrh,cvt%arrv)
      CASE(1) !  linear elements
        ALLOCATE(cvt%arr(nqty,0:mx,0:my))
        NULLIFY(cvt%arri,cvt%arrh,cvt%arrv)
      CASE(2:) !  higher-order elements
        ALLOCATE(cvt%arr(nqty,0:mx,0:my))
        ALLOCATE(cvt%arrh(nqty,poly_degree-1,1:mx,0:my))
        ALLOCATE(cvt%arrv(nqty,poly_degree-1,0:mx,1:my))
        ALLOCATE(cvt%arri(nqty,(poly_degree-1)**2,1:mx,1:my))
      END SELECT
!-----------------------------------------------------------------------
!     if the optional input for temporary arrays is present, allocate
!     the arrtmp arrays; otherwise, nullify them.  note that nbt is
!     the number of temporary bases, and nqt is the number of quantities
!     at each basis.
!-----------------------------------------------------------------------
      IF (PRESENT(nbt)) THEN
        ALLOCATE(cvt%arrtmp(nqt,nbt,1:mx,1:my))
      ELSE
        NULLIFY(cvt%arrtmp)
      ENDIF
!-----------------------------------------------------------------------
!     register this object.
!-----------------------------------------------------------------------
      NULLIFY(cvt%mem_id)
#ifdef OBJ_MEM_PROF
      SELECT CASE(poly_degree)
      CASE(0)
        sz=SIZEOF(cvt%arri)
      CASE(1)
        sz=SIZEOF(cvt%arr)
      CASE(2:)
        sz= SIZEOF(cvt%arri)+SIZEOF(cvt%arr)                            &
     &     +SIZEOF(cvt%arrh)+SIZEOF(cvt%arrv)
      END SELECT
      CALL memlogger%update(cvt%mem_id,'cvector_2D','unknown',sz)
#endif
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_2D_ctype_alloc
!-----------------------------------------------------------------------
!     subprogram 39. vector_2D_ctype_dealloc.
!     deallocates space for a 2D complex vector_type structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_2D_ctype_dealloc(cvt)
#ifdef OBJ_MEM_PROF
      USE memlog, only: memlogger
#endif

      TYPE(cvector_2D_type), INTENT(INOUT) :: cvt

      IF (ASSOCIATED(cvt%arr)) THEN
        DEALLOCATE(cvt%arr)
        NULLIFY(cvt%arr)
      ENDIF
      IF (ASSOCIATED(cvt%arrh)) THEN
        DEALLOCATE(cvt%arrh)
        NULLIFY(cvt%arrh)
      ENDIF
      IF (ASSOCIATED(cvt%arrv)) THEN
        DEALLOCATE(cvt%arrv)
        NULLIFY(cvt%arrv)
      ENDIF
      IF (ASSOCIATED(cvt%arri)) THEN
        DEALLOCATE(cvt%arri)
        NULLIFY(cvt%arri)
      ENDIF
      IF (ASSOCIATED(cvt%arrtmp)) THEN
        DEALLOCATE(cvt%arrtmp)
        NULLIFY(cvt%arrtmp)
      ENDIF
!-----------------------------------------------------------------------
!     unregister this object.
!-----------------------------------------------------------------------
#ifdef OBJ_MEM_PROF
      CALL memlogger%update(cvt%mem_id,' ',' ',0,resize=.TRUE.)
#endif
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_2D_ctype_dealloc
!-----------------------------------------------------------------------
!     subprogram 40. cvector_2D_assign_rsc.
!     assign a real scalar value to a 2D complex vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assign_rsc(vec,rscalar)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rscalar

!-----------------------------------------------------------------------
!     if the grid vertex-centered array is allocated, treat as a
!     standard element.  If, not the structure represents piecewise
!     constant.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=rscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rscalar
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_2D_assign_rsc: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assign_rsc
!-----------------------------------------------------------------------
!     subprogram 41. cvector_2D_assign_csc.
!     assign a 2D complex scalar value to a 2D complex vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assign_csc(vec,cscalar)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec
      COMPLEX(r8), INTENT(IN) :: cscalar

!-----------------------------------------------------------------------
!     if the grid vertex-centered array is allocated, treat as a
!     standard element.  If, not the structure represents piecewise
!     constant.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=cscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=cscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=cscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=cscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=cscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=cscalar
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_2D_assign_csc: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assign_csc
!-----------------------------------------------------------------------
!     subprogram 42. cvector_2D_assign_int.
!     assign a integer value to a 2D complex vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assign_int(vec,int)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

!-----------------------------------------------------------------------
!     if the grid vertex-centered array is allocated, treat as a
!     standard element.  If, not the structure represents piecewise
!     constant.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int
        IF (ASSOCIATED(vec%arri)) vec%arri=int
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_2D_assign_int: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assign_int
!-----------------------------------------------------------------------
!     subprogram 43. cvector_2D_assign_cvec2.
!     set one 2D complex vector structure equal to another.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assign_cvec2(vec1,vec2)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(cvector_2D_type), INTENT(IN) :: vec2

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(:,:,:)=vec2%arr(:,:,:)
        IF (ASSOCIATED(vec1%arrh)) vec1%arrh(:,:,:,:)=vec2%arrh(:,:,:,:)
        IF (ASSOCIATED(vec1%arrv)) vec1%arrv(:,:,:,:)=vec2%arrv(:,:,:,:)
        IF (ASSOCIATED(vec1%arri)) vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp(:,:,:,:)=vec2%arrtmp(:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_2D_assign_cvec2: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assign_cvec2
!-----------------------------------------------------------------------
!     subprogram 44. cvector_2D_assignq_cvec2
!     set one 2D complex vector structure equal to another.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assignq_cvec2(vec1,vec2,nqty)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(cvector_2D_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: nqty

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(1:nqty,:,:)=vec2%arr(1:nqty,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &     vec1%arrh(1:nqty,:,:,:)=vec2%arrh(1:nqty,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &     vec1%arrv(1:nqty,:,:,:)=vec2%arrv(1:nqty,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &     vec1%arri(1:nqty,:,:,:)=vec2%arri(1:nqty,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &     vec1%arrtmp(1:nqty,:,:,:)=vec2%arrtmp(1:nqty,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(1:nqty,:,:,:)=vec2%arri(1:nqty,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_2D_assignq_cvec2: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assignq_cvec2
!-----------------------------------------------------------------------
!     subprogram 45. cvector_2D_assignp_cvec2
!     set part of a 2D complex vector structure equal to another.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assignp_cvec2(vec1,vec2,nq1,nq2,sz)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(cvector_2D_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: nq1,nq2,sz

      INTEGER(i4) :: ne1,ne2

      ne1=nq1+sz-1_i4
      ne2=nq2+sz-1_i4
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(nq1:ne1,:,:)=vec2%arr(nq2:ne2,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &     vec1%arrh(nq1:ne1,:,:,:)=vec2%arrh(nq2:ne2,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &     vec1%arrv(nq1:ne1,:,:,:)=vec2%arrv(nq2:ne2,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &     vec1%arri(nq1:ne1,:,:,:)=vec2%arri(nq2:ne2,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &     vec1%arrtmp(nq1:ne1,:,:,:)=vec2%arrtmp(nq2:ne2,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(nq1:ne1,:,:,:)=vec2%arri(nq2:ne2,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_2D_assignp_cvec2: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assignp_cvec2
!-----------------------------------------------------------------------
!     subprogram 46. cvector_2D_assign_vec
!     set one 2D complex vector structure equal to another.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assign_vec(vec1,vec2)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(:,:,:)=vec2%arr(:,:,:)
        IF (ASSOCIATED(vec1%arrh)) vec1%arrh(:,:,:,:)=vec2%arrh(:,:,:,:)
        IF (ASSOCIATED(vec1%arrv)) vec1%arrv(:,:,:,:)=vec2%arrv(:,:,:,:)
        IF (ASSOCIATED(vec1%arri)) vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_2D_assign_vec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assign_vec
!-----------------------------------------------------------------------
!     subprogram 47. cvector_2D_assignq_vec
!     set nqty of one 2D complex vector structure equal to another.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assignq_vec(vec1,vec2,nqty)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: nqty

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(1:nqty,:,:)=vec2%arr(1:nqty,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &     vec1%arrh(1:nqty,:,:,:)=vec2%arrh(1:nqty,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &     vec1%arrv(1:nqty,:,:,:)=vec2%arrv(1:nqty,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &     vec1%arri(1:nqty,:,:,:)=vec2%arri(1:nqty,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &     vec1%arrtmp(1:nqty,:,:,:)=vec2%arrtmp(1:nqty,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(1:nqty,:,:,:)=vec2%arri(1:nqty,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_2D_assignq_vec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assignq_vec
!-----------------------------------------------------------------------
!     subprogram 48. cvector_2D_assignp_vec.
!     set one complex 2D vector structure equal to part of a
!     vector structure. Assign vec1(nq1:sz-1+nq1)=vec2(nq2:sz-1+nq2).
!     No checks on compatibility are performed.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assignp_vec(vec1,vec2,nq1,nq2,sz)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: nq1,nq2,sz

      INTEGER(i4) :: ne1,ne2

      ne1=nq1+sz-1
      ne2=nq2+sz-1
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(nq1:ne1,:,:)=vec2%arr(nq2:ne2,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &    vec1%arrh(nq1:ne1,:,:,:)=vec2%arrh(nq2:ne2,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &    vec1%arrv(nq1:ne1,:,:,:)=vec2%arrv(nq2:ne2,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &    vec1%arri(nq1:ne1,:,:,:)=vec2%arri(nq2:ne2,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp(nq1:ne1,:,:,:)=vec2%arrtmp(nq2:ne2,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(nq1:ne1,:,:,:)=vec2%arri(nq2:ne2,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_2D_assign_part_cvec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assignp_vec
!-----------------------------------------------------------------------
!     subprogram 49. cvector_2D_assign_cvec.
!     set one complex 2D vector structure equal to part of a complex
!     vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assign_cvec(vec1,vec2,fcomp,nqty)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: fcomp
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqty

!-----------------------------------------------------------------------
!     if the number of vector components is specified with nqty, limit
!     transfer.
!-----------------------------------------------------------------------
      IF (PRESENT(nqty)) THEN
        IF (ASSOCIATED(vec1%arr)) THEN
          vec1%arr(1:nqty,:,:)=vec2%arr(1:nqty,:,:,fcomp)
          IF (ASSOCIATED(vec1%arrh))                                    &
     &      vec1%arrh(1:nqty,:,:,:)=vec2%arrh(1:nqty,:,:,:,fcomp)
          IF (ASSOCIATED(vec1%arrv))                                    &
     &      vec1%arrv(1:nqty,:,:,:)=vec2%arrv(1:nqty,:,:,:,fcomp)
          IF (ASSOCIATED(vec1%arri))                                    &
     &      vec1%arri(1:nqty,:,:,:)=vec2%arri(1:nqty,:,:,:,fcomp)
          IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))      &
     &      vec1%arrtmp(1:nqty,:,:,:)=vec2%arrtmp(1:nqty,:,:,:,fcomp)
        ELSE IF (ASSOCIATED(vec1%arri)) THEN
          vec1%arri(1:nqty,:,:,:)=vec2%arri(1:nqty,:,:,:,fcomp)
        ELSE
          CALL nim_stop                                                 &
     &      ('Cvector_2D_assign_cvec: vector arrays not associated.')
        ENDIF
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      ELSE
        IF (ASSOCIATED(vec1%arr)) THEN
          vec1%arr=vec2%arr(:,:,:,fcomp)
          IF (ASSOCIATED(vec1%arrh))                                    &
     &      vec1%arrh=vec2%arrh(:,:,:,:,fcomp)
          IF (ASSOCIATED(vec1%arrv))                                    &
     &      vec1%arrv=vec2%arrv(:,:,:,:,fcomp)
          IF (ASSOCIATED(vec1%arri))                                    &
     &      vec1%arri=vec2%arri(:,:,:,:,fcomp)
          IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))      &
     &      vec1%arrtmp=vec2%arrtmp(:,:,:,:,fcomp)
        ELSE IF (ASSOCIATED(vec1%arri)) THEN
          vec1%arri=vec2%arri(:,:,:,:,fcomp)
        ELSE
          CALL nim_stop                                                 &
     &      ('Cvector_2D_assign_cvec: vector arrays not associated.')
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assign_cvec
!-----------------------------------------------------------------------
!     subprogram 50. cvector_2D_assignp_cvec.
!     set one complex 2D vector structure equal to part of a complex
!     vector structure. Assign vec1(nq1:sz-1+nq1)=vec2(nq2:sz-1+nq2).
!     No checks on compatibility are performed.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assignp_cvec(vec1,vec2,fcomp,nq1,nq2,sz)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: fcomp,nq1,nq2,sz

      INTEGER(i4) :: ne1,ne2

      ne1=nq1+sz-1
      ne2=nq2+sz-1
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(nq1:ne1,:,:)=vec2%arr(nq2:ne2,:,:,fcomp)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &    vec1%arrh(nq1:ne1,:,:,:)=vec2%arrh(nq2:ne2,:,:,:,fcomp)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &    vec1%arrv(nq1:ne1,:,:,:)=vec2%arrv(nq2:ne2,:,:,:,fcomp)
        IF (ASSOCIATED(vec1%arri))                                      &
     &    vec1%arri(nq1:ne1,:,:,:)=vec2%arri(nq2:ne2,:,:,:,fcomp)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp(nq1:ne1,:,:,:)=vec2%arrtmp(nq2:ne2,:,:,:,fcomp)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(nq1:ne1,:,:,:)=vec2%arri(nq2:ne2,:,:,:,fcomp)
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_2D_assign_part_cvec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assignp_cvec
!-----------------------------------------------------------------------
!     subprogram 51. cvector_2D_add_cvec2.
!     add one 2D complex vector structure to another.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_add_cvec2(vec1,vec2,v1fac,v2fac)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(cvector_2D_type), INTENT(IN) :: vec2
      REAL(r8), INTENT(IN), OPTIONAL :: v1fac,v2fac

      REAL(r8) :: v1f,v2f
!-----------------------------------------------------------------------
!     set coefficients to input if used.
!-----------------------------------------------------------------------
      IF (PRESENT(v1fac)) THEN
        v1f=v1fac
      ELSE
        v1f=1
      ENDIF
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1
      ENDIF
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr=v1f*vec1%arr+v2f*vec2%arr
        IF (ASSOCIATED(vec1%arrh)) vec1%arrh=v1f*vec1%arrh+v2f*vec2%arrh
        IF (ASSOCIATED(vec1%arrv)) vec1%arrv=v1f*vec1%arrv+v2f*vec2%arrv
        IF (ASSOCIATED(vec1%arri)) vec1%arri=v1f*vec1%arri+v2f*vec2%arri
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp=v1f*vec1%arrtmp+v2f*vec2%arrtmp
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri=v1f*vec1%arri+v2f*vec2%arri
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_2D_add_cvec2: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_add_cvec2
!-----------------------------------------------------------------------
!     subprogram 52. cvector_2D_addc_cvec2.
!     add one 2D complex vector structure to another with a complex
!     coefficients.  here, only v2fac is optional, so that
!     there is no confusion with cvector_add_cvec in the module
!     procedure definition.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_addc_cvec2(vec1,vec2,v1f,v2fac)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(cvector_2D_type), INTENT(IN) :: vec2
      COMPLEX(r8), INTENT(IN) :: v1f
      COMPLEX(r8), INTENT(IN), OPTIONAL :: v2fac

      COMPLEX(r8) :: v2f
!-----------------------------------------------------------------------
!     set coefficient to input if used.
!-----------------------------------------------------------------------
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1._r8
      ENDIF
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr=v1f*vec1%arr+v2f*vec2%arr
        IF (ASSOCIATED(vec1%arrh)) vec1%arrh=v1f*vec1%arrh+v2f*vec2%arrh
        IF (ASSOCIATED(vec1%arrv)) vec1%arrv=v1f*vec1%arrv+v2f*vec2%arrv
        IF (ASSOCIATED(vec1%arri)) vec1%arri=v1f*vec1%arri+v2f*vec2%arri
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp=v1f*vec1%arrtmp+v2f*vec2%arrtmp
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri=v1f*vec1%arri+v2f*vec2%arri
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_2D_addc_cvec2: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_addc_cvec2
!-----------------------------------------------------------------------
!     subprogram 53. cvector_2D_mult_rsc.
!     multiply a complex 2D vector structure by a real scalar.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_mult_rsc(vec,rsc,nq1,nq2)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rsc
      INTEGER(i4), INTENT(IN), OPTIONAL :: nq1,nq2

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (PRESENT(nq2)) THEN
        IF (ASSOCIATED(vec%arr)) THEN
          vec%arr(nq1:nq2,:,:)=rsc*vec%arr(nq1:nq2,:,:)
          IF (ASSOCIATED(vec%arrh))                                     &
     &      vec%arrh(nq1:nq2,:,:,:)=rsc*vec%arrh(nq1:nq2,:,:,:)
          IF (ASSOCIATED(vec%arrv))                                     &
     &      vec%arrv(nq1:nq2,:,:,:)=rsc*vec%arrv(nq1:nq2,:,:,:)
          IF (ASSOCIATED(vec%arri))                                     &
     &      vec%arri(nq1:nq2,:,:,:)=rsc*vec%arri(nq1:nq2,:,:,:)
          IF (ASSOCIATED(vec%arrtmp))                                   &
     &      vec%arrtmp(nq1:nq2,:,:,:)=rsc*vec%arrtmp(nq1:nq2,:,:,:)
        ELSE IF (ASSOCIATED(vec%arri)) THEN
          vec%arri(nq1:nq2,:,:,:)=rsc*vec%arri(nq1:nq2,:,:,:)
        ELSE
          CALL nim_stop                                                 &
     &      ('Cvector_2D_mult_rsc: vector arrays not associated.')
        ENDIF
      ELSE
        IF (ASSOCIATED(vec%arr)) THEN
          vec%arr=rsc*vec%arr
          IF (ASSOCIATED(vec%arrh)) vec%arrh=rsc*vec%arrh
          IF (ASSOCIATED(vec%arrv)) vec%arrv=rsc*vec%arrv
          IF (ASSOCIATED(vec%arri)) vec%arri=rsc*vec%arri
        ELSE IF (ASSOCIATED(vec%arri)) THEN
          vec%arri=rsc*vec%arri
        ELSE
          CALL nim_stop                                                 &
     &      ('Cvector_2D_mult_rsc: vector arrays not associated.')
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_mult_rsc
!-----------------------------------------------------------------------
!     subprogram 54. cvector_2D_mult_csc.
!     multiply a complex 2D vector structure by a complex scalar.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_mult_csc(vec,csc,nq1,nq2)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec
      COMPLEX(r8), INTENT(IN) :: csc
      INTEGER(i4), INTENT(IN), OPTIONAL :: nq1,nq2

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (PRESENT(nq2)) THEN
        IF (ASSOCIATED(vec%arr)) THEN
          vec%arr(nq1:nq2,:,:)=csc*vec%arr(nq1:nq2,:,:)
          IF (ASSOCIATED(vec%arrh))                                     &
     &      vec%arrh(nq1:nq2,:,:,:)=csc*vec%arrh(nq1:nq2,:,:,:)
          IF (ASSOCIATED(vec%arrv))                                     &
     &      vec%arrv(nq1:nq2,:,:,:)=csc*vec%arrv(nq1:nq2,:,:,:)
          IF (ASSOCIATED(vec%arri))                                     &
     &      vec%arri(nq1:nq2,:,:,:)=csc*vec%arri(nq1:nq2,:,:,:)
          IF (ASSOCIATED(vec%arrtmp))                                     &
     &      vec%arrtmp(nq1:nq2,:,:,:)=csc*vec%arrtmp(nq1:nq2,:,:,:)
        ELSE IF (ASSOCIATED(vec%arri)) THEN
          vec%arri(nq1:nq2,:,:,:)=csc*vec%arri(nq1:nq2,:,:,:)
        ELSE
          CALL nim_stop                                                 &
     &      ('Cvector_2D_mult_csc: vector arrays not associated.')
        ENDIF
      ELSE
        IF (ASSOCIATED(vec%arr)) THEN
          vec%arr=csc*vec%arr
          IF (ASSOCIATED(vec%arrh)) vec%arrh=csc*vec%arrh
          IF (ASSOCIATED(vec%arrv)) vec%arrv=csc*vec%arrv
          IF (ASSOCIATED(vec%arri)) vec%arri=csc*vec%arri
          IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=csc*vec%arrtmp
        ELSE IF (ASSOCIATED(vec%arri)) THEN
          vec%arri=csc*vec%arri
        ELSE
          CALL nim_stop                                                 &
     &      ('Cvector_2D_mult_csc: vector arrays not associated.')
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_mult_csc
!-----------------------------------------------------------------------
!     subprogram 55. cvector_2D_mult_int.
!     multiply a complex 2D vector structure by a real scalar.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_mult_int(vec,int)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=int*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int*vec%arri
      ELSE
        CALL nim_stop                                                   &
     &    ('Cvector_2D_mult_int: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_mult_int
!-----------------------------------------------------------------------
!     subprogram 55.1. cvector_2D_pack_cvec.
!     packs the interior and discontinuous coefficients of a cvector
!     into the interior coefficients of a complex 2D vector.
!
!     the arri array of the 2D vector should be allocated to hold all
!     of the quantity and basis indices as a flat array in its first
!     array index.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_pack_cvec(cvt,cv2p,nqi,nbi,nqd,nbd,fcomp)

      INTEGER(i4), INTENT(IN) :: nqi,nbi,nqd,nbd,fcomp
      TYPE(cvector_type), INTENT(IN) :: cvt
      TYPE(cvector_2D_type), INTENT(OUT) :: cv2p

      INTEGER(i4) :: ix,iy,iv,ib,iq
!-----------------------------------------------------------------------
!     at each element, loop over the interior and discontinuous bases
!     separately.
!-----------------------------------------------------------------------
      DO iy=1,SIZE(cvt%arri,4)
        DO ix=1,SIZE(cvt%arri,3)
          iv=1
          DO ib=1,nbi
            DO iq=1,nqi
              cv2p%arri(iv,1,ix,iy)=cvt%arri(iq,ib,ix,iy,fcomp)
              iv=iv+1
            ENDDO
          ENDDO
          DO ib=1,nbd
            DO iq=1,nqd
              cv2p%arri(iv,1,ix,iy)=cvt%arrtmp(iq,ib,ix,iy,fcomp)
              iv=iv+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_pack_cvec
!-----------------------------------------------------------------------
!     subprogram 55.2. cvector_2D_pack_cvec2.
!     packs the interior and discontinuous coefficients of a 2D complex
!     vector into the interior coefficients of another complex 2D
!     vector.
!
!     the arri array of the 2D vector should be allocated to hold all
!     of the quantity and basis indices as a flat array in its first
!     array index.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_pack_cvec2(cv2in,cv2p,nqi,nbi,nqd,nbd)

      INTEGER(i4), INTENT(IN) :: nqi,nbi,nqd,nbd
      TYPE(cvector_2D_type), INTENT(IN) :: cv2in
      TYPE(cvector_2D_type), INTENT(OUT) :: cv2p

      INTEGER(i4) :: ix,iy,iv,ib,iq
!-----------------------------------------------------------------------
!     at each element, loop over the interior and discontinuous bases
!     separately.
!-----------------------------------------------------------------------
      DO iy=1,SIZE(cv2in%arri,4)
        DO ix=1,SIZE(cv2in%arri,3)
          iv=1
          DO ib=1,nbi
            DO iq=1,nqi
              cv2p%arri(iv,1,ix,iy)=cv2in%arri(iq,ib,ix,iy)
              iv=iv+1
            ENDDO
          ENDDO
          DO ib=1,nbd
            DO iq=1,nqd
              cv2p%arri(iv,1,ix,iy)=cv2in%arrtmp(iq,ib,ix,iy)
              iv=iv+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_pack_cvec2
!-----------------------------------------------------------------------
!     subprogram 55.3. cvector_2D_unpack_cvec.
!     unpacks the interior coefficients of a complex 2D vector into
!     another complex 2D vector and a cvector, respectively.
!
!     the arri array of the packed 2D vector is allocated to hold all
!     of the quantity and basis indices as a flat array in its first
!     array index.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_unpack_cvec(cv2p,cv2i,cvt,nqi,nbi,nqd,nbd,  &
     &                                  fcomp)

      INTEGER(i4), INTENT(IN) :: nqi,nbi,nqd,nbd,fcomp
      TYPE(cvector_2D_type), INTENT(IN) :: cv2p
      TYPE(cvector_2D_type), INTENT(OUT) :: cv2i
      TYPE(cvector_type), INTENT(OUT) :: cvt

      INTEGER(i4) :: ix,iy,iv,ib,iq
!-----------------------------------------------------------------------
!     at each element, loop over the interior and discontinuous bases
!     separately.
!-----------------------------------------------------------------------
      DO iy=1,SIZE(cv2p%arri,4)
        DO ix=1,SIZE(cv2p%arri,3)
          iv=1
          DO ib=1,nbi
            DO iq=1,nqi
              cv2i%arri(iq,ib,ix,iy)=cv2p%arri(iv,1,ix,iy)
              iv=iv+1
            ENDDO
          ENDDO
          DO ib=1,nbd
            DO iq=1,nqd
              cvt%arrtmp(iq,ib,ix,iy,fcomp)=cv2p%arri(iv,1,ix,iy)
              iv=iv+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_unpack_cvec
!-----------------------------------------------------------------------
!     subprogram 55.4. cvector_2D_unpack_add_cvec.
!     unpacks the interior coefficients of a complex 2D vector into
!     another complex 2D vector and a cvector, respectively.
!     here the result is added to the interior of the cvector.
!
!     the arri array of the packed 2D vector is allocated to hold all
!     of the quantity and basis indices as a flat array in its first
!     array index.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_unpack_add_cvec(cv2p,cv2i,cvt,nqi,nbi,nqd,  &
     &                                      nbd,fcomp)

      INTEGER(i4), INTENT(IN) :: nqi,nbi,nqd,nbd,fcomp
      TYPE(cvector_2D_type), INTENT(IN) :: cv2p
      TYPE(cvector_2D_type), INTENT(OUT) :: cv2i
      TYPE(cvector_type), INTENT(INOUT) :: cvt

      INTEGER(i4) :: ix,iy,iv,ib,iq
!-----------------------------------------------------------------------
!     at each element, loop over the interior and discontinuous bases
!     separately.
!-----------------------------------------------------------------------
      DO iy=1,SIZE(cv2p%arri,4)
        DO ix=1,SIZE(cv2p%arri,3)
          iv=1
          DO ib=1,nbi
            DO iq=1,nqi
              cv2i%arri(iq,ib,ix,iy)=cv2p%arri(iv,1,ix,iy)
              iv=iv+1
            ENDDO
          ENDDO
          DO ib=1,nbd
            DO iq=1,nqd
              cvt%arri(iq,ib,ix,iy,fcomp)=                              &
     &          cvt%arri(iq,ib,ix,iy,fcomp)+cv2p%arri(iv,1,ix,iy)
              iv=iv+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_unpack_add_cvec
!-----------------------------------------------------------------------
!     subprogram 55.5. cvector_2D_unpack_cvec2.
!     unpacks the interior coefficients of a complex 2D vector into
!     another complex 2D vector, skipping the indices from a
!     discontinuous field if unpd is false.
!
!     the arri array of the packed 2D vector is allocated to hold all
!     of the quantity and basis indices as a flat array in its first
!     array index.
!-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_unpack_cvec2(cv2p,cv2i,nqi,nbi,nqd,nbd,unpd)

      INTEGER(i4), INTENT(IN) :: nqi,nbi,nqd,nbd
      LOGICAL, INTENT(IN) :: unpd
      TYPE(cvector_2D_type), INTENT(IN) :: cv2p
      TYPE(cvector_2D_type), INTENT(OUT) :: cv2i

      INTEGER(i4) :: ix,iy,iv,ib,iq
!-----------------------------------------------------------------------
!     at each element, loop over the interior bases.
!-----------------------------------------------------------------------
      DO iy=1,SIZE(cv2p%arri,4)
        DO ix=1,SIZE(cv2p%arri,3)
          iv=1
          DO ib=1,nbi
            DO iq=1,nqi
              cv2i%arri(iq,ib,ix,iy)=cv2p%arri(iv,1,ix,iy)
              iv=iv+1
            ENDDO
          ENDDO
          IF (unpd) THEN
            DO ib=1,nbd
              DO iq=1,nqd
                cv2i%arrtmp(iq,ib,ix,iy)=cv2p%arri(iv,1,ix,iy)
                iv=iv+1
              ENDDO
            ENDDO
          ELSE
            iv=iv+nqd*nbd
          ENDIF
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_unpack_cvec2
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      END MODULE cvector_2D_type_mod


!-----------------------------------------------------------------------
!     module for 3D-in-space arrays of real vectors quantities.
!-----------------------------------------------------------------------
      MODULE rvector_3D_type_mod
      USE local
      USE vector_defn_type_mod
      IMPLICIT NONE

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 56. vector_3D_rtype_alloc.
!     allocates space for a real 3D vector_type structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_3D_rtype_alloc(rvt,poly_degree,mx,my,nqty,nph,  &
     &                                 nbt,nqt)
#ifdef OBJ_MEM_PROF
      USE memlog, only: memlogger
#endif

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,poly_degree,nph
      INTEGER(i4), INTENT(IN), OPTIONAL :: nbt,nqt
      TYPE(vector_3D_type), INTENT(OUT) :: rvt

      INTEGER(i4) :: sz
!-----------------------------------------------------------------------
!     allocate space according to the basis functions needed.
!     if poly_degree is non-positive, storage for discontinuous-field
!     coefficients is allocated.
!-PRE triangles will need something here, too.
!-----------------------------------------------------------------------
      SELECT CASE(poly_degree)
      CASE(:0)  !  piecewise continuous fields
        ALLOCATE(rvt%arri(nqty,(poly_degree-1)**2,mx,my,nph))
        NULLIFY(rvt%arr,rvt%arrh,rvt%arrv)
      CASE(1) !  linear elements
        ALLOCATE(rvt%arr(nqty,0:mx,0:my,nph))
        NULLIFY(rvt%arri,rvt%arrh,rvt%arrv)
      CASE(2:) !  higher-order elements
        ALLOCATE(rvt%arr(nqty,0:mx,0:my,nph))
        ALLOCATE(rvt%arrh(nqty,poly_degree-1,1:mx,0:my,nph))
        ALLOCATE(rvt%arrv(nqty,poly_degree-1,0:mx,1:my,nph))
        ALLOCATE(rvt%arri(nqty,(poly_degree-1)**2,1:mx,1:my,nph))
      END SELECT
!-----------------------------------------------------------------------
!     if the optional input for temporary arrays is present, allocate
!     the arrtmp arrays; otherwise, nullify them.  note that nbt is
!     the number of temporary bases, and nqt is the number of quantities
!     at each basis.
!-----------------------------------------------------------------------
      IF (PRESENT(nbt)) THEN
        ALLOCATE(rvt%arrtmp(nqt,nbt,1:mx,1:my,nph))
      ELSE
        NULLIFY(rvt%arrtmp)
      ENDIF
!-----------------------------------------------------------------------
!     register this object.
!-----------------------------------------------------------------------
      NULLIFY(rvt%mem_id)
#ifdef OBJ_MEM_PROF
      SELECT CASE(poly_degree)
      CASE(0)
        sz=SIZEOF(rvt%arri)
      CASE(1)
        sz=SIZEOF(rvt%arr)
      CASE(2:)
        sz= SIZEOF(rvt%arri)+SIZEOF(rvt%arr)                            &
     &     +SIZEOF(rvt%arrh)+SIZEOF(rvt%arrv)
      END SELECT
      CALL memlogger%update(rvt%mem_id,'rvector_3D','unknown',sz)
#endif
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_rtype_alloc
!-----------------------------------------------------------------------
!     subprogram 57. vector_3D_rtype_dealloc.
!     deallocates space for a real 3D vector_type structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_3D_rtype_dealloc(rvt)
#ifdef OBJ_MEM_PROF
      USE memlog, only: memlogger
#endif

      TYPE(vector_3D_type), INTENT(INOUT) :: rvt

      IF (ASSOCIATED(rvt%arr)) THEN
        DEALLOCATE(rvt%arr)
        NULLIFY(rvt%arr)
      ENDIF
      IF (ASSOCIATED(rvt%arrh)) THEN
        DEALLOCATE(rvt%arrh)
        NULLIFY(rvt%arrh)
      ENDIF
      IF (ASSOCIATED(rvt%arrv)) THEN
        DEALLOCATE(rvt%arrv)
        NULLIFY(rvt%arrv)
      ENDIF
      IF (ASSOCIATED(rvt%arri)) THEN
        DEALLOCATE(rvt%arri)
        NULLIFY(rvt%arri)
      ENDIF
      IF (ASSOCIATED(rvt%arrtmp)) THEN
        DEALLOCATE(rvt%arrtmp)
        NULLIFY(rvt%arrtmp)
      ENDIF
!-----------------------------------------------------------------------
!     unregister this object.
!-----------------------------------------------------------------------
#ifdef OBJ_MEM_PROF
      CALL memlogger%update(rvt%mem_id,' ',' ',0,resize=.TRUE.)
#endif
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_rtype_dealloc
!-----------------------------------------------------------------------
!     subprogram 58. vector_3D_assign_rsc.
!     assign a real scalar value to a 3D vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_3D_assign_rsc(vec,rscalar)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rscalar

!-----------------------------------------------------------------------
!     if the grid vertex-centered array is allocated, treat as a
!     standard element.  If, not the structure represents piecewise
!     constant.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=rscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rscalar
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_3D_assign_rsc: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_assign_rsc
!-----------------------------------------------------------------------
!     subprogram 59. vector_3D_assign_csc.
!     assign the real part of a complex scalar value to a 3D vector
!     structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_3D_assign_csc(vec,cscalar)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec
      COMPLEX(r8), INTENT(IN) :: cscalar

!-----------------------------------------------------------------------
!     if the grid vertex-centered array is allocated, treat as a
!     standard element.  If, not the structure represents piecewise
!     constant.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=cscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=cscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=cscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=cscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=cscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=cscalar
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_3D_assign_csc: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_assign_csc
!-----------------------------------------------------------------------
!     subprogram 60. vector_3D_assign_int.
!     assign an integer value to a 3D vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_3D_assign_int(vec,int)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

!-----------------------------------------------------------------------
!     if the grid vertex-centered array is allocated, treat as a
!     standard element.  If, not the structure represents piecewise
!     constant.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int
        IF (ASSOCIATED(vec%arri)) vec%arri=int
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_3D_assign_int: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_assign_int
!-----------------------------------------------------------------------
!     subprogram 61. vector_3D_assign_vec.
!     set a 3D vector structure equal to a vector structure at the
!     specified index of the periodic coordinate.
!-----------------------------------------------------------------------
      SUBROUTINE vector_3D_assign_vec(vec1,vec2,ip)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: ip
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(:,:,:,ip)=vec2%arr(:,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &    vec1%arrh(:,:,:,:,ip)=vec2%arrh(:,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &    vec1%arrv(:,:,:,:,ip)=vec2%arrv(:,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &    vec1%arri(:,:,:,:,ip)=vec2%arri(:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp(:,:,:,:,ip)=vec2%arrtmp(:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(:,:,:,:,ip)=vec2%arri(:,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_3D_assign_vec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_assign_vec
!-----------------------------------------------------------------------
!     subprogram 62. vector_3D_assign_vec3.
!     set a 3D vector structure equal to another.
!-----------------------------------------------------------------------
      SUBROUTINE vector_3D_assign_vec3(vec1,vec2)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec1
      TYPE(vector_3D_type), INTENT(IN) :: vec2
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(:,:,:,:)=vec2%arr(:,:,:,:)
        IF (ASSOCIATED(vec1%arrh))                                      &
     &    vec1%arrh(:,:,:,:,:)=vec2%arrh(:,:,:,:,:)
        IF (ASSOCIATED(vec1%arrv))                                      &
     &    vec1%arrv(:,:,:,:,:)=vec2%arrv(:,:,:,:,:)
        IF (ASSOCIATED(vec1%arri))                                      &
     &    vec1%arri(:,:,:,:,:)=vec2%arri(:,:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp(:,:,:,:,:)=vec2%arrtmp(:,:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(:,:,:,:,:)=vec2%arri(:,:,:,:,:)
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_3D_assign_vec3: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_assign_vec3
!-----------------------------------------------------------------------
!     subprogram 63. vector_3D_add_vec.
!     add a 2D vector structure to every index of the periodic
!     coordinate of a 3D vector structure.
!-----------------------------------------------------------------------
      SUBROUTINE vector_3D_add_vec(vec1,vec2,v1fac,v2fac)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      REAL(r8), INTENT(IN), OPTIONAL :: v1fac,v2fac

      INTEGER(i4) :: ip,np
      REAL(r8) :: v1f,v2f
!-----------------------------------------------------------------------
!     set coefficients to input if used.
!-----------------------------------------------------------------------
      IF (PRESENT(v1fac)) THEN
        v1f=v1fac
      ELSE
        v1f=1
      ENDIF
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1
      ENDIF
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        np=SIZE(vec1%arr,4)
        DO ip=1,np
          vec1%arr(:,:,:,ip)=v1f*vec1%arr(:,:,:,ip)+v2f*vec2%arr
        ENDDO
        IF (ASSOCIATED(vec1%arrh)) THEN
          DO ip=1,np
            vec1%arrh(:,:,:,:,ip)=v1f*vec1%arrh(:,:,:,:,ip)+            &
     &                            v2f*vec2%arrh
          ENDDO
        ENDIF
        IF (ASSOCIATED(vec1%arrv)) THEN
          DO ip=1,np
            vec1%arrv(:,:,:,:,ip)=v1f*vec1%arrv(:,:,:,:,ip)+            &
     &                            v2f*vec2%arrv
          ENDDO
        ENDIF
        IF (ASSOCIATED(vec1%arri)) THEN
          DO ip=1,np
            vec1%arri(:,:,:,:,ip)=v1f*vec1%arri(:,:,:,:,ip)+            &
     &                            v2f*vec2%arri
          ENDDO
        ENDIF
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp)) THEN
          DO ip=1,np
            vec1%arrtmp(:,:,:,:,ip)=v1f*vec1%arrtmp(:,:,:,:,ip)+        &
     &                            v2f*vec2%arrtmp
          ENDDO
        ENDIF
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        np=SIZE(vec1%arri,5)
        DO ip=1,np
          vec1%arri(:,:,:,:,ip)=v1f*vec1%arri(:,:,:,:,ip)+v2f*vec2%arri
        ENDDO
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_3D_add_vec: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_add_vec
!-----------------------------------------------------------------------
!     subprogram 64. vector_3D_add_vec3.
!     add one 3D vector structure to another.
!-----------------------------------------------------------------------
      SUBROUTINE vector_3D_add_vec3(vec1,vec2,v1fac,v2fac)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec1
      TYPE(vector_3D_type), INTENT(IN) :: vec2
      REAL(r8), INTENT(IN), OPTIONAL :: v1fac,v2fac

      REAL(r8) :: v1f,v2f
!-----------------------------------------------------------------------
!     set coefficients to input if used.
!-----------------------------------------------------------------------
      IF (PRESENT(v1fac)) THEN
        v1f=v1fac
      ELSE
        v1f=1
      ENDIF
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1
      ENDIF
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr=v1f*vec1%arr+v2f*vec2%arr
        IF (ASSOCIATED(vec1%arrh)) vec1%arrh=v1f*vec1%arrh+v2f*vec2%arrh
        IF (ASSOCIATED(vec1%arrv)) vec1%arrv=v1f*vec1%arrv+v2f*vec2%arrv
        IF (ASSOCIATED(vec1%arri)) vec1%arri=v1f*vec1%arri+v2f*vec2%arri
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))        &
     &    vec1%arrtmp=v1f*vec1%arrtmp+v2f*vec2%arrtmp
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri=v1f*vec1%arri+v2f*vec2%arri
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_3D_add_vec3: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_add_vec3
!-----------------------------------------------------------------------
!     subprogram 65. vector_3D_mult_rsc.
!     multiply a 3D vector structure by a real scalar.
!-----------------------------------------------------------------------
      SUBROUTINE vector_3D_mult_rsc(vec,rsc)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rsc

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rsc*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rsc*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rsc*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=rsc*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rsc*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rsc*vec%arri
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_3D_mult_rsc: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_mult_rsc
!-----------------------------------------------------------------------
!     subprogram 66. vector_3D_mult_int.
!     multiply a 3D vector structure by an integer.
!-----------------------------------------------------------------------
      SUBROUTINE vector_3D_mult_int(vec,int)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=int*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int*vec%arri
      ELSE
        CALL nim_stop                                                   &
     &    ('Vector_3D_mult_int: vector arrays not associated.')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_mult_int
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      END MODULE rvector_3D_type_mod


!-----------------------------------------------------------------------
!     the generic module contains overloaded assignment and operator
!     and interfaces.
!-----------------------------------------------------------------------
      MODULE vector_type_mod
      USE rvector_type_mod
      USE cvector_type_mod
      USE cvector_2D_type_mod
      USE rvector_3D_type_mod
      IMPLICIT NONE

      INTERFACE ASSIGNMENT (=)
        MODULE PROCEDURE vector_assign_rsc,vector_assign_csc,           &
     &    vector_assign_int,vector_assign_vec,cvector_assign_rsc,       &
     &    cvector_assign_csc,cvector_assign_int,cvector_assign_cvec,    &
     &    cvector_2D_assign_rsc,cvector_2D_assign_csc,                  &
     &    cvector_2D_assign_int,cvector_2D_assign_cvec2,                &
     &    vector_3D_assign_rsc,vector_3D_assign_csc,                    &
     &    vector_3D_assign_int,vector_3D_assign_vec3
      END INTERFACE

      INTERFACE vector_add
        MODULE PROCEDURE vector_add_vec,cvector_add_cvec,               &
     &    cvector_2D_add_cvec2,cvector_addc_cvec,                       &
     &    cvector_2D_addc_cvec2,vector_3D_add_vec,vector_3D_add_vec3
      END INTERFACE

      INTERFACE vector_mult
        MODULE PROCEDURE vector_mult_rsc,vector_mult_int,               &
     &    cvector_mult_rsc,cvector_mult_int,                            &
     &    cvector_2D_mult_rsc,cvector_2D_mult_int,                      &
     &    cvector_mult_csc,cvector_2D_mult_csc,vector_3D_mult_rsc,      &
     &    vector_3D_mult_int
      END INTERFACE

      INTERFACE vector_type_alloc
        MODULE PROCEDURE vector_rtype_alloc,vector_ctype_alloc,         &
     &    vector_2D_ctype_alloc,vector_3D_rtype_alloc
      END INTERFACE

      INTERFACE vector_type_dealloc
        MODULE PROCEDURE vector_rtype_dealloc,vector_ctype_dealloc,     &
     &    vector_2D_ctype_dealloc,vector_3D_rtype_dealloc
      END INTERFACE
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      END MODULE vector_type_mod
