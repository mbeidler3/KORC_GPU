!-----------------------------------------------------------------------
!     $Id: matrix_type_mod.F90 7623 2022-01-18 03:55:24Z tbechtel $
!     contains a module that defines structures for saving matrices
!     and their factors.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.
!-----------------------------------------------------------------------
!     0. matrix_type_mod.
!     1. matrix_rbl_real_alloc
!     2. matrix_rbl_real_dealloc
!     3. matrix_rbl_comp_alloc
!     4. matrix_rbl_comp_dealloc
!     5. matrix_rbl_make_real
!     6. matrix_rbl_dof_init
!     7. matrix_tbl_real_alloc
!     8. matrix_tbl_real_dealloc
!     9. matrix_tbl_comp_alloc
!     10. matrix_tbl_comp_dealloc
!     11. matrix_tbl_make_real
!     12. factor_ptr_copy_rmat_info
!     13. factor_dealloc_rmat_info
!     14. factor_ptr_copy_cmat_info
!     15. factor_dealloc_cmat_info
!     16. matrix_factor_real_log
!     17. matrix_factor_comp_log
!     18. sparsity_pattern_log.
!-----------------------------------------------------------------------
!     declarations.
!-----------------------------------------------------------------------
#include "config.f"
#ifdef HAVE_PASTIX
#include "pastix_fortran.h"
#endif
#ifdef HAVE_PARDISO
#ifdef HAVE_MPI
#include "mkl_cluster_sparse_solver.f90"
#else
#include "mkl_pardiso.f90"
#endif
#endif
      MODULE matrix_type_mod
      USE local
#ifdef HAVE_SUPERLU_DIST
      USE superlu_dist_mod, ONLY: slu_matrix_storage
#endif
#ifdef HAVE_PARDISO
#ifdef HAVE_MPI
      USE mkl_cluster_sparse_solver, ONLY:                              &
     &                          mkl_cluster_sparse_solver_handle
#else
      USE mkl_pardiso, ONLY: mkl_pardiso_handle
#endif
#endif
      IMPLICIT NONE
#ifdef HAVE_MUMPS
      INCLUDE "zmumps_struc.h"
      INCLUDE "dmumps_struc.h"
#endif
!-----------------------------------------------------------------------
!     types used for defining matrix structures for a single Fourier
!     component.
!-----------------------------------------------------------------------
      TYPE :: global_matrix_type
        TYPE(rbl_mat_type), DIMENSION(:), POINTER :: rbl_mat
        TYPE(tbl_mat_type), DIMENSION(:), POINTER :: tbl_mat
        INTEGER(i4) :: fcomp,foff
        INTEGER(i4) :: nqty,nqdis
        CHARACTER(1), DIMENSION(:), POINTER :: vcomp
        LOGICAL :: eliminated
        LOGICAL :: symmetric
        REAL(r8) :: diag_scale
        CHARACTER(16) :: essential_cond
      END TYPE global_matrix_type
      TYPE :: complex_matrix_type
        TYPE(rbl_comp_mat_type), DIMENSION(:), POINTER :: rbl_mat
        TYPE(tbl_comp_mat_type), DIMENSION(:), POINTER :: tbl_mat
        INTEGER(i4) :: fcomp,foff
        INTEGER(i4) :: nqty,nqdis
        CHARACTER(1), DIMENSION(:), POINTER :: vcomp
        LOGICAL :: eliminated
        LOGICAL :: hermitian
        REAL(r8) :: diag_scale
        CHARACTER(16) :: essential_cond
      END TYPE complex_matrix_type
!-----------------------------------------------------------------------
!     types used for holding matrix elements in rblocks.
!     nbasis_el is the number of nonzero basis functions per element.
!     nbasis_cont and nbasis_disc are the numbers of continuous and
!     discontinuous bases, respectively.
!
!     the mat arrays are dimensioned nbtypeXnbtype to cover
!     the matrix coupling among the different basis types.
!
!     ndof_el is the number of degrees of freedom per element, and
!     the pointer arrays dof_ib, dof_iv, dof_ix, dof_iy, and dof_iq hold
!     information that helps transfer integrals from element-based
!     storage to structured rblock storage.  for each degree of freedom,
!
!       dof_ib = element basis-function index
!       dof_iv = element (physical-field) vector-component index
!       dof_ix = block ix offset relative to element (=-1 or 0)
!       dof_iy = block iy offset relative to element (=-1 or 0)
!       dof_iq = block quantity index
!
!     also, for each basis type, den_type holds the ending DOF index.
!-----------------------------------------------------------------------
      TYPE :: rbl_mat_type
        INTEGER(i4) :: nbasis_el
        INTEGER(i4) :: nbasis_cont
        INTEGER(i4) :: nbasis_disc
        INTEGER(i4) :: nbtype
        INTEGER(i4) :: ndof_el
        INTEGER(i4) :: mx,my
        INTEGER(i4), DIMENSION(:), POINTER :: ix0,iy0
        INTEGER(i4), DIMENSION(:), POINTER :: nb_type,nq_type
        INTEGER(i4), DIMENSION(:), POINTER :: dof_ib,dof_iv,dof_ix,     &
     &               dof_iy,dof_iq
        INTEGER(i4), DIMENSION(:), POINTER :: den_type
        TYPE(arr_6d_type), DIMENSION(:,:), POINTER :: mat
        INTEGER(i4), POINTER :: mem_id
      END TYPE rbl_mat_type
      TYPE :: rbl_comp_mat_type
        INTEGER(i4) :: nbasis_el
        INTEGER(i4) :: nbasis_cont
        INTEGER(i4) :: nbasis_disc
        INTEGER(i4) :: nbtype
        INTEGER(i4) :: ndof_el
        INTEGER(i4) :: mx,my
        INTEGER(i4), DIMENSION(:), POINTER :: ix0,iy0
        INTEGER(i4), DIMENSION(:), POINTER :: nb_type,nq_type
        INTEGER(i4), DIMENSION(:), POINTER :: dof_ib,dof_iv,dof_ix,     &
     &               dof_iy,dof_iq
        INTEGER(i4), DIMENSION(:), POINTER :: den_type
        TYPE(comp_arr_6d_type), DIMENSION(:,:), POINTER :: mat
        INTEGER(i4), POINTER :: mem_id
      END TYPE rbl_comp_mat_type
      TYPE :: rbl_mat_info
        INTEGER(i4) :: nbasis_el
        INTEGER(i4) :: nbtype
        INTEGER(i4) :: mx,my
        INTEGER(i4), DIMENSION(:), POINTER :: ix0,iy0
        INTEGER(i4), DIMENSION(:), POINTER :: nb_type,nq_type
      END TYPE rbl_mat_info
!-----------------------------------------------------------------------
!     types used for holding rblock matrix elements.
!-----------------------------------------------------------------------
      TYPE :: arr_6d_type
        REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER :: arr
      END TYPE arr_6d_type
      TYPE :: comp_arr_6d_type
        COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER :: arr
      END TYPE comp_arr_6d_type
!-----------------------------------------------------------------------
!     types used for holding matrix elements in rblocks.
!-----------------------------------------------------------------------
      TYPE :: tbl_mat_type
        INTEGER(i4), DIMENSION(:), POINTER :: ix0,iy0
        TYPE(matrix_element_type3), DIMENSION(:), POINTER :: lmat
      END TYPE tbl_mat_type
      TYPE :: tbl_comp_mat_type
        INTEGER(i4), DIMENSION(:), POINTER :: ix0,iy0
        TYPE(comp_matrix_element_type3), DIMENSION(:), POINTER :: lmat
      END TYPE tbl_comp_mat_type
!-----------------------------------------------------------------------
!     types used for holding tblock matrix elements.  it was formerly
!     part of tblock_type_mod.
!-----------------------------------------------------------------------
      TYPE :: matrix_element_type3
        REAL(r8), DIMENSION(:,:,:), POINTER :: element
        INTEGER(i4), DIMENSION(:), POINTER :: from_vert
      END TYPE matrix_element_type3
      TYPE :: comp_matrix_element_type3
        COMPLEX(r8), DIMENSION(:,:,:), POINTER :: element
        INTEGER(i4), DIMENSION(:), POINTER :: from_vert
      END TYPE comp_matrix_element_type3
!-----------------------------------------------------------------------
!     type for holding the distributed col/row pointer for sent data.
!-----------------------------------------------------------------------
      TYPE :: send_start_acc
        INTEGER(i4), ALLOCATABLE :: data(:,:)
      END TYPE send_start_acc
!-----------------------------------------------------------------------
!     type for holding the distributed row pointer for sent data.
!-----------------------------------------------------------------------
      TYPE :: send_start_row
        INTEGER(i4), ALLOCATABLE :: irow(:)
      END TYPE send_start_row
!-----------------------------------------------------------------------
!     type for holding the sparsity pattern for direct solves.
!-----------------------------------------------------------------------
      TYPE :: sparsity_pattern
        !> running total of unique nodes by block (1:nbl_total)
        !> where nrow=irowst_block(nbl_total) and the starting row for
        !> any block is given by irowst_block(id-1)
        INTEGER(i4), ALLOCATABLE :: irowst_block(:)
        !> global (local if sparsity distributed) row/col index
        INTEGER(iSolve), POINTER :: j_acc(:)
        !> global (local if sparsity distributed) matrix col/row ptr
        INTEGER(iSolve), POINTER :: start_acc(:)
        !> distributed (CSC/CSR) matrix: 0-index col/row ptr
        INTEGER(iSolve), POINTER :: start_loc(:)
        !> save location for j_acc if overwritten by the external solver
        INTEGER(iSolve), ALLOCATABLE :: j_acc_save(:)
        !> save location if overwritten by the external solver
        INTEGER(iSolve), ALLOCATABLE :: start_acc_save(:)
        INTEGER(iSolve) :: nrow   !< number of rows (global mat)
        INTEGER(iSolve) :: nnz    !< number of non-zeros (global mat)
        INTEGER(iSolve) :: nbw    !< maximum width between irow and icol
        INTEGER(iSolve) :: mloc   !< number of local rows
        INTEGER(iSolve) :: fstrow !< first local row
        INTEGER(iSolve) :: lstrow !< last local row
        INTEGER(iSolve) :: nnzloc !< local number of non-zeros
        INTEGER(i4) :: indst      !< local starting index of val/j_acc
        INTEGER(i4) :: indend     !< local ending index of val/j_acc
        INTEGER(i4) :: nsnd       !< node number of sends
        INTEGER(i4) :: nrcv       !< node number of receives
        !> list of nodes that send to this node (size nrcv)
        INTEGER(i4), ALLOCATABLE :: recvlist(:)
        !> list of size of receives (size nrcv)
        INTEGER(i4), ALLOCATABLE :: recvtot(:)
        !> location in acc of received data (size max(recvtot),nrcv)
        INTEGER(i4), ALLOCATABLE :: recvjentry(:,:)
        !> list of nodes that this node is sends to (size nsnd)
        INTEGER(i4), ALLOCATABLE :: sendlist(:)
        !> list of size of sends (size nsnd)
        INTEGER(i4), ALLOCATABLE :: sendtot(:)
        !> location in jentry() where send starts (size nsnd)
        INTEGER(i4), ALLOCATABLE :: sendst(:)
        !> row index where send starts (size nsnd)
        INTEGER(i4), ALLOCATABLE :: sendrowst(:)
        !> col/row ptr for sends
        TYPE(send_start_acc), ALLOCATABLE :: sendstacc(:)
        !> location in acc of sent data (col ind if dist sparsity)
        INTEGER(i4), ALLOCATABLE :: jentry(:)
        !> asynchronous receive request handle
        INTEGER(i4), ALLOCATABLE :: recv_req(:)
        !> asynchronous send request handle
        INTEGER(i4), ALLOCATABLE :: send_req(:)
        !> for global matrix gather communication
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: algdispl,algcount
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: algdsplr,algcntr
        !> true if matrix has been factored previously
        LOGICAL :: acc_lustored
        !> true if the external solver distributed interface is used
        LOGICAL :: matrix_distributed
        !> true if the sparsity pattern is distributed
        LOGICAL :: sparsity_distributed
        !> true if a direct solver is being used
        LOGICAL :: direct_solver
        !> SuperLU handle for bridge routine storage
        INTEGER(i8) :: acc_handle
#ifdef HAVE_SUPERLU_DIST
        !> SuperLU matrix data
        TYPE(slu_matrix_storage) :: slu_matrix_data
#endif
        !> memory ids for external solvers
        INTEGER(i4), POINTER :: extsolver_mem_id1
        INTEGER(i4), POINTER :: extsolver_mem_id2
        INTEGER(i4), POINTER :: extsolver_mem_id3
        !> integer options
        INTEGER(i4) :: iopts(solve_nopts)
        !> double precision options
        REAL(r8) :: dopts(solve_nopts)
        !> factor to reduce threads for external solver
        REAL(r8) :: thread_reduct_fac
#ifdef HAVE_PASTIX
        !> data structure for PaStiX solver
        pastix_data_ptr_t :: pastix_data
        !> integer pastix input/output parameters
        pastix_int_t :: iparm(IPARM_SIZE)
        !> floating point pastix input/output parameters
        double precision :: dparm(DPARM_SIZE)
        pastix_int_t, POINTER, DIMENSION(:) :: perm
        pastix_int_t, POINTER, DIMENSION(:) :: invp
        pastix_int_t, POINTER, DIMENSION(:) :: loc2glb
#endif /* HAVE_PASTIX */
#ifdef HAVE_MUMPS
        !> data structure for MUMPS solver
        TYPE(zmumps_struc) :: zmumps_par
        TYPE(dmumps_struc) :: dmumps_par
        ! Mumps distribution solution communication arrays
        !> running total of unique FE nodes by mpi node (0:nprocs_layer)
        !> this is a condensed version of irowst_block
        INTEGER(i4), ALLOCATABLE :: irowst_node(:)
        !> # of sends and receives for mumps dist soln (1:nprocs_layer)
        INTEGER(i4), ALLOCATABLE :: mumps_nsend(:),mumps_nrecv(:)
        !> inode,indices in sol_loc of send soln (2,1:nsend)
        TYPE(send_start_row), ALLOCATABLE :: mumps_sisol(:)
        !> local indices of recv soln (1:nprocs)
        TYPE(send_start_row), ALLOCATABLE :: mumps_risol(:)
        !> map of sol_loc to local indices of loc soln
        INTEGER(i4), ALLOCATABLE :: mumps_lisol(:,:)
        !> local solution values array
        COMPLEX(r8), POINTER :: mumps_sol_loc(:)
        !> local solution index array
        INTEGER(i4), POINTER :: mumps_isol_loc(:)
        !> number of local values in sol_loc
        INTEGER(i4) :: nlsol
        !> max value of mumps_nsend
        INTEGER(i4) :: maxsend
        !> max value of mumps_nrecv
        INTEGER(i4) :: maxrecv
        !> number of non-zero sized receives
        INTEGER(i4) :: ncrecv
        !> boolean for reset comm arrays
        LOGICAL :: reset_comm_arrays
#endif
#ifdef HAVE_HYPRE
        !> pointers to the hypre storage of the matrix, RHS, soln and
        !> solver.
        INTEGER(i8) :: AA,bb,xx,solver,precond
        !> max number of entries per row
        INTEGER(i4) :: rownnzmax
        !> global row indices of the 2x2 real local matrix
        INTEGER(iSolve), ALLOCATABLE :: row2(:)
        !> nqty of the matrix
        INTEGER(i4) :: nqty
        !> iteration count, number of times hypre solve called in count
        INTEGER(i4) :: its2d,navg
#endif /* HAVE_HYPRE */
#ifdef HAVE_PARDISO
        !> pardiso work space
#ifdef HAVE_MPI
        TYPE(mkl_cluster_sparse_solver_handle) :: pt(64)
#else
        TYPE(mkl_pardiso_handle) :: pt(64)
#endif /* HAVE_MPI */
        !> pardiso input/output array
        INTEGER(i4) :: iparm(64)
#endif /* HAVE_PARDISO */
      END TYPE sparsity_pattern
!-----------------------------------------------------------------------
!     types used to define arrays for matrix factors and working
!     pointers.
!-----------------------------------------------------------------------
      TYPE :: matrix_factor_type
        TYPE(bl_fac_type),      DIMENSION(:), POINTER :: bl_fac
        TYPE(bl_sparsity_type), DIMENSION(:), POINTER :: bl_spp
        TYPE(rbl_mat_info),     DIMENSION(:), POINTER :: mat_info
        TYPE(sparsity_pattern) :: spp
        INTEGER(i4) :: nqty,ilu_fill
        !> factors for precon=='lapack'
        REAL(r8), DIMENSION(:,:), ALLOCATABLE :: a11
        !> memory id for memory logging.
        INTEGER(i4), POINTER :: mem_id
        CHARACTER(1) :: mat_name
      END TYPE matrix_factor_type
      TYPE :: complex_factor_type
        TYPE(bl_comp_fac_type), DIMENSION(:), POINTER :: bl_fac
        TYPE(bl_sparsity_type), DIMENSION(:), POINTER :: bl_spp
        TYPE(rbl_mat_info),     DIMENSION(:), POINTER :: mat_info
        TYPE(sparsity_pattern) :: spp
        INTEGER(i4) :: nqty,ilu_fill
        !> factors for precon=='lapack'
        COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: a11
        !> memory id for memory logging.
        INTEGER(i4), POINTER :: mem_id
        CHARACTER(1) :: mat_name
      END TYPE complex_factor_type
!-----------------------------------------------------------------------
!     block by block matrix factor arrays.  no one preconditioner uses
!     all of these arrays, and only those needed are allocated.
!-----------------------------------------------------------------------
      TYPE :: bl_fac_type
        REAL(r8), DIMENSION(:,:,:), POINTER :: lt_dg0,ut_dg0
        REAL(r8), DIMENSION(:,:,:,:), POINTER :: lt_dg1,ut_dg1,         &
     &            xelim,yelim,xlc,xcl,ylc,ycl
        REAL(r8), DIMENSION(:,:,:,:,:), POINTER :: lt_cen,lt_off,lt_per,&
     &            ut_cen,ut_off,ut_per,adix,adiy
        TYPE(arr_4d_type), DIMENSION(:), POINTER :: pmat
        INTEGER(i4) :: nv,mx,my
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ix0,iy0
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: nb_type,nq_type
        LOGICAL :: perblock,degenerate
        !> memory id for memory logging.
        INTEGER(i4), POINTER :: mem_id
      END TYPE bl_fac_type
      TYPE :: bl_comp_fac_type
        COMPLEX(r8), DIMENSION(:,:,:), POINTER :: lt_dg0,ut_dg0
        COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: lt_dg1,ut_dg1,      &
     &               xelim,yelim,xlc,xcl,ylc,ycl
        COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: lt_cen,lt_off,    &
     &            lt_per,ut_cen,ut_off,ut_per,adix,adiy
        TYPE(comp_arr_4d_type), DIMENSION(:), POINTER :: pmat
        INTEGER(i4) :: nv,mx,my
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ix0,iy0
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: nb_type,nq_type
        LOGICAL :: perblock,degenerate
        !> memory id for memory logging.
        INTEGER(i4), POINTER :: mem_id
      END TYPE bl_comp_fac_type
!-----------------------------------------------------------------------
!     types used for holding inverses of diagonal point-block.
!-----------------------------------------------------------------------
      TYPE :: arr_4d_type
        REAL(r8), DIMENSION(:,:,:,:), POINTER :: arr
      END TYPE arr_4d_type
      TYPE :: comp_arr_4d_type
        COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: arr
      END TYPE comp_arr_4d_type
!-----------------------------------------------------------------------
!     type for block specific information for direct solves.
!-----------------------------------------------------------------------
      TYPE :: bl_sparsity_type
        TYPE(row_type), DIMENSION(:), POINTER :: row_ind
        LOGICAL :: perblock,degenerate
        !> memory id for memory logging.
        INTEGER(i4), POINTER :: mem_id
      END TYPE bl_sparsity_type
!-----------------------------------------------------------------------
!     type for assigning row numbers and organizing direct solves.
!-----------------------------------------------------------------------
      TYPE :: row_type
        INTEGER(i4), DIMENSION(:,:,:), POINTER :: rarr
      END TYPE row_type

!-----------------------------------------------------------------------
!     interfaces for the allocation and deallocation routines.
!-----------------------------------------------------------------------
      INTERFACE matrix_rbl_alloc
        MODULE PROCEDURE matrix_rbl_real_alloc,matrix_rbl_comp_alloc
      END INTERFACE

      INTERFACE matrix_rbl_dealloc
        MODULE PROCEDURE matrix_rbl_real_dealloc,matrix_rbl_comp_dealloc
      END INTERFACE

      INTERFACE matrix_tbl_alloc
        MODULE PROCEDURE matrix_tbl_real_alloc,matrix_tbl_comp_alloc
      END INTERFACE

      INTERFACE matrix_tbl_dealloc
        MODULE PROCEDURE matrix_tbl_real_dealloc,matrix_tbl_comp_dealloc
      END INTERFACE

      INTERFACE factor_ptr_copy_mat_info
        MODULE PROCEDURE factor_ptr_copy_rmat_info,                     &
     &                   factor_ptr_copy_cmat_info
      END INTERFACE

      INTERFACE factor_dealloc_mat_info
        MODULE PROCEDURE factor_dealloc_rmat_info,                      &
     &                   factor_dealloc_cmat_info
      END INTERFACE

#ifdef OBJ_MEM_PROF
      INTERFACE matrix_factor_log
        MODULE PROCEDURE matrix_factor_real_log,matrix_factor_comp_log
      END INTERFACE
#endif

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 1. matrix_rbl_real_alloc
!     allocate arrays needed for a real rblock matrix.
!-----------------------------------------------------------------------
      SUBROUTINE matrix_rbl_real_alloc(rmat,mx,my,nqc,pdc,nqdin,nbdin)
#ifdef OBJ_MEM_PROF
      USE memlog, only: memlogger
#endif
      TYPE(rbl_mat_type), INTENT(OUT) :: rmat
      INTEGER(i4), INTENT(IN) :: mx,my,nqc,pdc
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqdin,nbdin

      INTEGER(i4) :: id,jd,ixst,iyst,jxst,jyst,x0off,x1off,y0off,y1off, &
     &               sz,nqd,nbd

      IF (PRESENT(nqdin)) THEN
        nqd=nqdin
      ELSE 
        nqd=0_i4
      ENDIF
      IF (PRESENT(nbdin)) THEN
        nbd=nbdin
      ELSE 
        nbd=0_i4
      ENDIF
!-----------------------------------------------------------------------
!     the mat array allows multiple basis function types (grid vertex,
!     horizontal element side, vertical element side, and interior-
!     centered).
!
!     the structures accommodate equations with nqc continuous fields
!     of polynomial degree pdc and nqd discontinuous fields with
!     nbd basis functions.  the storage for 'element interiors'
!     is used for both the interior coefficients of continuous fields
!     and for all coefficients of discontinuous fields.
!
!     if poly_degree=1, there is only one basis type.  otherwise, there
!     are 4.
!-----------------------------------------------------------------------
      rmat%nbasis_el=MIN(nqc,1_i4)*(pdc+1)**2+MIN(nqd,1_i4)*nbd
      rmat%nbasis_cont=MIN(nqc,1_i4)*(pdc+1)**2
      rmat%nbasis_disc=MIN(nqd,1_i4)*nbd
      rmat%ndof_el=nqc*(pdc+1)**2+nqd*nbd
      IF (pdc==1.AND.nqd==0) THEN
        rmat%nbtype=1
      ELSE
        rmat%nbtype=4
      ENDIF
      rmat%mx=mx
      rmat%my=my
      ALLOCATE(rmat%mat(rmat%nbtype,rmat%nbtype))
      ALLOCATE(rmat%ix0(rmat%nbtype))
      ALLOCATE(rmat%iy0(rmat%nbtype))
      ALLOCATE(rmat%nb_type(rmat%nbtype))
      ALLOCATE(rmat%nq_type(rmat%nbtype))
      ALLOCATE(rmat%den_type(rmat%nbtype))
      ALLOCATE(rmat%dof_ib(rmat%ndof_el))
      ALLOCATE(rmat%dof_iv(rmat%ndof_el))
      ALLOCATE(rmat%dof_ix(rmat%ndof_el))
      ALLOCATE(rmat%dof_iy(rmat%ndof_el))
      ALLOCATE(rmat%dof_iq(rmat%ndof_el))
!-----------------------------------------------------------------------
!     logical indices for each of the basis types.
!       index 1 is grid vertex-centered.
!       index 2 is horizontal side-centered.
!       index 3 is vertical side-centered.
!       index 4 is interior-centered and all discontinuous bases.
!
!     this version has the 6D matrix array indices defined
!     (col_comp,col_x_off,col_y_off,row_comp,row_x_index,row_y_index),
!     where comp is vector component and basis for types with multiple
!     bases (vector component varying faster) and off is the offset
!     from the row index.
!-----------------------------------------------------------------------
      DO id=1,rmat%nbtype
        SELECT CASE(id)
        CASE(1)
          rmat%ix0(id)=0
          rmat%iy0(id)=0
          rmat%nb_type(id)=MIN(nqc,1_i4)
          rmat%nq_type(id)=nqc
          rmat%den_type(id)=4*nqc
        CASE(2)
          rmat%ix0(id)=1
          rmat%iy0(id)=0
          rmat%nb_type(id)=MIN(nqc,1_i4)*(pdc-1)
          rmat%nq_type(id)=nqc*(pdc-1)
          rmat%den_type(id)=rmat%den_type(id-1)+2*nqc*(pdc-1)
        CASE(3)
          rmat%ix0(id)=0
          rmat%iy0(id)=1
          rmat%nb_type(id)=MIN(nqc,1_i4)*(pdc-1)
          rmat%nq_type(id)=nqc*(pdc-1)
          rmat%den_type(id)=rmat%den_type(id-1)+2*nqc*(pdc-1)
        CASE(4)
          rmat%ix0(id)=1
          rmat%iy0(id)=1
          rmat%nb_type(id)=MIN(nqc,1_i4)*(pdc-1)**2+                    &
     &                     MIN(nqd,1_i4)*nbd
          rmat%nq_type(id)=nqc*(pdc-1)**2+nqd*nbd
          rmat%den_type(id)=rmat%den_type(id-1)+rmat%nq_type(id)
        END SELECT
      ENDDO
      DO id=1,rmat%nbtype
        ixst=rmat%ix0(id)
        iyst=rmat%iy0(id)
        DO jd=1,rmat%nbtype
          jxst=rmat%ix0(jd)
          jyst=rmat%iy0(jd)
          x0off=jxst-1
          x1off=1-ixst
          y0off=jyst-1
          y1off=1-iyst
          ALLOCATE(rmat%mat(jd,id)%                                     &
     &      arr(rmat%nq_type(jd),x0off:x1off,y0off:y1off,               &
     &          rmat%nq_type(id),ixst:mx,iyst:my))
            rmat%mat(jd,id)%arr=0
        ENDDO
      ENDDO
      CALL matrix_rbl_dof_init(rmat%dof_ib,rmat%dof_iv,rmat%dof_ix,     &
     &                         rmat%dof_iy,rmat%dof_iq,nqc,pdc,nqd,nbd)
!-----------------------------------------------------------------------
!     register this object.
!-----------------------------------------------------------------------
      NULLIFY(rmat%mem_id)
#ifdef OBJ_MEM_PROF
      sz= SIZEOF(rmat%mat)+SIZEOF(rmat%ix0)+SIZEOF(rmat%iy0)            &
     &   +SIZEOF(rmat%nb_type)+SIZEOF(rmat%nq_type)+SIZEOF(rmat)        &
     &   +SIZEOF(rmat%dof_ib)+SIZEOF(rmat%dof_iv)+SIZEOF(rmat%dof_ix)   &
     &   +SIZEOF(rmat%dof_iy)+SIZEOF(rmat%dof_iq)+SIZEOF(rmat%den_type)
      DO id=1,rmat%nbtype
        DO jd=1,rmat%nbtype
          sz=sz+SIZEOF(rmat%mat(jd,id)%arr)
        ENDDO
      ENDDO
      CALL memlogger%update(rmat%mem_id,'rmat','unknown',sz)
#endif
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_rbl_real_alloc
!-----------------------------------------------------------------------
!     subprogram 2. matrix_rbl_real_dealloc
!     deallocate arrays needed for a real rblock matrix.
!-----------------------------------------------------------------------
      SUBROUTINE matrix_rbl_real_dealloc(rmat)
#ifdef OBJ_MEM_PROF
      USE memlog, only: memlogger
#endif
      TYPE(rbl_mat_type), INTENT(INOUT) :: rmat

      INTEGER(i4) :: id,jd

!-----------------------------------------------------------------------
!     loop over different basis combinations and deallocte.
!-----------------------------------------------------------------------
      DO jd=1,SIZE(rmat%mat,2)
        DO id=1,SIZE(rmat%mat,1)
          DEALLOCATE(rmat%mat(id,jd)%arr)
        ENDDO
      ENDDO
      DEALLOCATE(rmat%mat,rmat%ix0,rmat%iy0)
      DEALLOCATE(rmat%nb_type,rmat%nq_type)
      DEALLOCATE(rmat%den_type,rmat%dof_ib,rmat%dof_iv,                 &
     &           rmat%dof_ix,rmat%dof_iy,rmat%dof_iq)
!-----------------------------------------------------------------------
!     unregister this object.
!-----------------------------------------------------------------------
#ifdef OBJ_MEM_PROF
      CALL memlogger%update(rmat%mem_id,'rmat',' ',0,resize=.TRUE.)
#endif
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_rbl_real_dealloc
!-----------------------------------------------------------------------
!     subprogram 3. matrix_rbl_comp_alloc
!     allocate arrays needed for a comp rblock matrix.
!-----------------------------------------------------------------------
      SUBROUTINE matrix_rbl_comp_alloc(cmat,mx,my,nqc,pdc,nqdin,nbdin)
#ifdef OBJ_MEM_PROF
      USE memlog, only: memlogger
#endif
      TYPE(rbl_comp_mat_type), INTENT(OUT) :: cmat
      INTEGER(i4), INTENT(IN) :: mx,my,nqc,pdc
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqdin,nbdin

      INTEGER(i4) :: id,jd,ixst,iyst,jxst,jyst,x0off,x1off,y0off,y1off, &
     &               sz,nqd,nbd

      IF (PRESENT(nqdin)) THEN
        nqd=nqdin
      ELSE 
        nqd=0_i4
      ENDIF
      IF (PRESENT(nbdin)) THEN
        nbd=nbdin
      ELSE 
        nbd=0_i4
      ENDIF
!-----------------------------------------------------------------------
!     the mat array allows multiple basis function types (grid vertex,
!     horizontal element side, vertical element side, and interior-
!     centered).
!
!     the structures accommodate equations with nqc continuous fields
!     of polynomial degree pdc and nqd discontinuous fields with
!     nbd basis functions.  the storage for 'element interiors'
!     is used for both the interior coefficients of continuous fields
!     and for all coefficients of discontinuous fields.
!
!     if poly_degree=1, there is only one basis type.  otherwise, there
!     are 4.
!-----------------------------------------------------------------------
      cmat%nbasis_el=MIN(nqc,1_i4)*(pdc+1)**2+MIN(nqd,1_i4)*nbd
      cmat%nbasis_cont=MIN(nqc,1_i4)*(pdc+1)**2
      cmat%nbasis_disc=MIN(nqd,1_i4)*nbd
      cmat%ndof_el=nqc*(pdc+1)**2+nqd*nbd
      cmat%mx=mx
      cmat%my=my
      IF (pdc==1.AND.nqd==0) THEN
        cmat%nbtype=1
      ELSE
        cmat%nbtype=4
      ENDIF
      ALLOCATE(cmat%mat(cmat%nbtype,cmat%nbtype))
      ALLOCATE(cmat%ix0(cmat%nbtype))
      ALLOCATE(cmat%iy0(cmat%nbtype))
      ALLOCATE(cmat%nb_type(cmat%nbtype))
      ALLOCATE(cmat%nq_type(cmat%nbtype))
      ALLOCATE(cmat%den_type(cmat%nbtype))
      ALLOCATE(cmat%dof_ib(cmat%ndof_el))
      ALLOCATE(cmat%dof_iv(cmat%ndof_el))
      ALLOCATE(cmat%dof_ix(cmat%ndof_el))
      ALLOCATE(cmat%dof_iy(cmat%ndof_el))
      ALLOCATE(cmat%dof_iq(cmat%ndof_el))
!-----------------------------------------------------------------------
!     logical indices for each of the basis types.
!       index 1 is grid vertex-centered.
!       index 2 is horizontal side-centered.
!       index 3 is vertical side-centered.
!       index 4 is interior-centered and all discontinuous bases.
!
!     this version has the 6D matrix array indices defined
!     (col_comp,col_x_off,col_y_off,row_comp,row_x_index,row_y_index),
!     where comp is vector component and basis for types with multiple
!     bases (vector component varying faster) and off is the offset
!     from the row index.
!-----------------------------------------------------------------------
      DO id=1,cmat%nbtype
        SELECT CASE(id)
        CASE(1)
          cmat%ix0(id)=0
          cmat%iy0(id)=0
          cmat%nb_type(id)=MIN(nqc,1_i4)
          cmat%nq_type(id)=nqc
          cmat%den_type(id)=4*nqc
        CASE(2)
          cmat%ix0(id)=1
          cmat%iy0(id)=0
          cmat%nb_type(id)=MIN(nqc,1_i4)*(pdc-1)
          cmat%nq_type(id)=nqc*(pdc-1)
          cmat%den_type(id)=cmat%den_type(id-1)+2*nqc*(pdc-1)
        CASE(3)
          cmat%ix0(id)=0
          cmat%iy0(id)=1
          cmat%nb_type(id)=MIN(nqc,1_i4)*(pdc-1)
          cmat%nq_type(id)=nqc*(pdc-1)
          cmat%den_type(id)=cmat%den_type(id-1)+2*nqc*(pdc-1)
        CASE(4)
          cmat%ix0(id)=1
          cmat%iy0(id)=1
          cmat%nb_type(id)=MIN(nqc,1_i4)*(pdc-1)**2+                    &
     &                     MIN(nqd,1_i4)*nbd
          cmat%nq_type(id)=nqc*(pdc-1)**2+nqd*nbd
          cmat%den_type(id)=cmat%den_type(id-1)+cmat%nq_type(id)
        END SELECT
      ENDDO
      DO id=1,cmat%nbtype
        ixst=cmat%ix0(id)
        iyst=cmat%iy0(id)
        DO jd=1,cmat%nbtype
          jxst=cmat%ix0(jd)
          jyst=cmat%iy0(jd)
          x0off=jxst-1
          x1off=1-ixst
          y0off=jyst-1
          y1off=1-iyst
          ALLOCATE(cmat%mat(jd,id)%                                     &
     &      arr(cmat%nq_type(jd),x0off:x1off,y0off:y1off,               &
     &          cmat%nq_type(id),ixst:mx,iyst:my))
            cmat%mat(jd,id)%arr=0
        ENDDO
      ENDDO
      CALL matrix_rbl_dof_init(cmat%dof_ib,cmat%dof_iv,cmat%dof_ix,     &
     &                         cmat%dof_iy,cmat%dof_iq,nqc,pdc,nqd,nbd)
!-----------------------------------------------------------------------
!     register this object.
!-----------------------------------------------------------------------
      NULLIFY(cmat%mem_id)
#ifdef OBJ_MEM_PROF
      sz= SIZEOF(cmat%mat)+SIZEOF(cmat%ix0)+SIZEOF(cmat%iy0)            &
     &   +SIZEOF(cmat%nb_type)+SIZEOF(cmat%nq_type)+SIZEOF(cmat)        &
     &   +SIZEOF(cmat%dof_ib)+SIZEOF(cmat%dof_iv)+SIZEOF(cmat%dof_ix)   &
     &   +SIZEOF(cmat%dof_iy)+SIZEOF(cmat%dof_iq)+SIZEOF(cmat%den_type)
      DO id=1,cmat%nbtype
        DO jd=1,cmat%nbtype
          sz=sz+SIZEOF(cmat%mat(jd,id)%arr)
        ENDDO
      ENDDO
      CALL memlogger%update(cmat%mem_id,'cmat','unknown',sz)
#endif
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_rbl_comp_alloc
!-----------------------------------------------------------------------
!     subprogram 4. matrix_rbl_comp_dealloc
!     deallocate arrays needed for a comp rblock matrix.
!-----------------------------------------------------------------------
      SUBROUTINE matrix_rbl_comp_dealloc(cmat)
#ifdef OBJ_MEM_PROF
      USE memlog, only: memlogger
#endif
      TYPE(rbl_comp_mat_type), INTENT(INOUT) :: cmat

      INTEGER(i4) :: id,jd

!-----------------------------------------------------------------------
!     loop over different basis combinations and deallocte.
!-----------------------------------------------------------------------
      DO jd=1,SIZE(cmat%mat,2)
        DO id=1,SIZE(cmat%mat,1)
          DEALLOCATE(cmat%mat(id,jd)%arr)
        ENDDO
      ENDDO
      DEALLOCATE(cmat%mat,cmat%ix0,cmat%iy0)
      DEALLOCATE(cmat%nb_type,cmat%nq_type)
      DEALLOCATE(cmat%den_type,cmat%dof_ib,cmat%dof_iv,                 &
     &           cmat%dof_ix,cmat%dof_iy,cmat%dof_iq)
!-----------------------------------------------------------------------
!     unregister this object.
!-----------------------------------------------------------------------
#ifdef OBJ_MEM_PROF
      CALL memlogger%update(cmat%mem_id,'cmat',' ',0,resize=.TRUE.)
#endif
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_rbl_comp_dealloc
!-----------------------------------------------------------------------
!     subprogram 5. matrix_rbl_make_real
!     set matrix elements of a comp rblock matrix to their real parts.
!-----------------------------------------------------------------------
      SUBROUTINE matrix_rbl_make_real(cmat)

      TYPE(rbl_comp_mat_type), INTENT(INOUT) :: cmat

      INTEGER(i4) :: id,jd

!-----------------------------------------------------------------------
!     loop over different basis combinations.
!-----------------------------------------------------------------------
      DO jd=1,SIZE(cmat%mat,2)
        DO id=1,SIZE(cmat%mat,1)
          cmat%mat(id,jd)%arr=REAL(cmat%mat(id,jd)%arr,r8)
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_rbl_make_real
!-----------------------------------------------------------------------
!     subprogram 6. matrix_rbl_dof_init
!     initialize arrays that help translate element-based data.
!
!     see matrix_rbl_real_alloc for definitions of nqc, pdc, nqd,
!     and nbd.
!-----------------------------------------------------------------------
      SUBROUTINE matrix_rbl_dof_init(dib,div,dix,diy,diq,               &
     &                               nqc,pdc,nqd,nbd)

      INTEGER(i4), DIMENSION(:), INTENT(OUT) :: dib,div,dix,diy,diq
      INTEGER(i4), INTENT(IN) :: nqc,pdc,nqd,nbd

      INTEGER(i4) :: ix,iy,iq,ii,id
!-----------------------------------------------------------------------
!     fill the dof arrays that are used to generate indices
!     when transferring information into rmat arrays.
!-----------------------------------------------------------------------
      ii=1
      DO iy=0,1        !     4 continuous grid-centered bases
        DO ix=0,1
          DO iq=1,nqc
            dib(ii)=ix+2*iy+1
            div(ii)=iq
            dix(ii)=ix-1
            diy(ii)=iy-1
            diq(ii)=iq
            ii=ii+1
          ENDDO
        ENDDO
      ENDDO
      DO id=0,pdc-2    !     continuous horizontal side bases
        DO iy=0,1
          DO iq=1,nqc
            dib(ii)=iy+2*id+5
            div(ii)=iq
            dix(ii)=0
            diy(ii)=iy-1
            diq(ii)=iq+nqc*id
            ii=ii+1
          ENDDO
        ENDDO
      ENDDO
      DO id=0,pdc-2    !     continuous vertical side bases
        DO ix=0,1
          DO iq=1,nqc
            dib(ii)=ix+2*id+5+2*(pdc-1)
            div(ii)=iq
            dix(ii)=ix-1
            diy(ii)=0
            diq(ii)=iq+nqc*id
            ii=ii+1
          ENDDO
        ENDDO
      ENDDO
      DO id=0,(pdc-1)**2-1  !  interior bases for continuous expansions
        DO iq=1,nqc
          dib(ii)=id+1+4*pdc
          div(ii)=iq
          dix(ii)=0
          diy(ii)=0
          diq(ii)=iq+nqc*id
          ii=ii+1
        ENDDO
      ENDDO
      DO id=0,nbd-1    !  bases for discontinuous expansions
        DO iq=1,nqd
          dib(ii)=id+1+MIN(nqc,1_i4)*(pdc+1)**2  ! continue interior #s
          div(ii)=iq
          dix(ii)=0
          diy(ii)=0
          diq(ii)=iq+nqd*id+nqc*(pdc-1)**2
          ii=ii+1
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      END SUBROUTINE matrix_rbl_dof_init
!-----------------------------------------------------------------------
!     subprogram 7. matrix_tbl_real_alloc
!     allocate arrays needed for a real tblock matrix.
!-----------------------------------------------------------------------
      SUBROUTINE matrix_tbl_real_alloc(tmat,tg,nq)
      USE tri_linear

      TYPE(tbl_mat_type), INTENT(OUT) :: tmat
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: nq

      INTEGER(i4) :: iv,nn
!-----------------------------------------------------------------------
!-PRE allocate matrix elements for linear bases only.
!-----------------------------------------------------------------------
      ALLOCATE(tmat%lmat(0:tg%mvert))
      DO iv=0,tg%mvert
        nn=SIZE(tg%neighbor(iv)%vertex)-1
        ALLOCATE(tmat%lmat(iv)%element(nq,nq,0:nn))
        ALLOCATE(tmat%lmat(iv)%from_vert(0:nn))
        tmat%lmat(iv)%from_vert=tg%neighbor(iv)%vertex
      ENDDO
      ALLOCATE(tmat%ix0(1))
      ALLOCATE(tmat%iy0(1))
      tmat%ix0=0
      tmat%iy0=0
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_tbl_real_alloc
!-----------------------------------------------------------------------
!     subprogram 8. matrix_tbl_real_dealloc
!     deallocate arrays needed for a real tblock matrix.
!-----------------------------------------------------------------------
      SUBROUTINE matrix_tbl_real_dealloc(tmat)

      TYPE(tbl_mat_type), INTENT(INOUT) :: tmat

      INTEGER(i4) :: iv
!-----------------------------------------------------------------------
!-PRE deallocate linear bases only.
!-----------------------------------------------------------------------
      DO iv=0,SIZE(tmat%lmat)-1
        DEALLOCATE(tmat%lmat(iv)%element)
        DEALLOCATE(tmat%lmat(iv)%from_vert)
      ENDDO
      DEALLOCATE(tmat%lmat,tmat%ix0,tmat%iy0)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_tbl_real_dealloc
!-----------------------------------------------------------------------
!     subprogram 9. matrix_tbl_comp_alloc
!     allocate arrays needed for a comp tblock matrix.
!-----------------------------------------------------------------------
      SUBROUTINE matrix_tbl_comp_alloc(tmat,tg,nq)
      USE tri_linear

      TYPE(tbl_comp_mat_type), INTENT(OUT) :: tmat
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: nq

      INTEGER(i4) :: iv,nn
!-----------------------------------------------------------------------
!-PRE allocate matrix elements for linear bases only.
!-----------------------------------------------------------------------
      ALLOCATE(tmat%lmat(0:tg%mvert))
      DO iv=0,tg%mvert
        nn=SIZE(tg%neighbor(iv)%vertex)-1
        ALLOCATE(tmat%lmat(iv)%element(nq,nq,0:nn))
        ALLOCATE(tmat%lmat(iv)%from_vert(0:nn))
        tmat%lmat(iv)%from_vert=tg%neighbor(iv)%vertex
      ENDDO
      ALLOCATE(tmat%ix0(1))
      ALLOCATE(tmat%iy0(1))
      tmat%ix0=0
      tmat%iy0=0
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_tbl_comp_alloc
!-----------------------------------------------------------------------
!     subprogram 10. matrix_tbl_comp_dealloc
!     deallocate arrays needed for a comp tblock matrix.
!-----------------------------------------------------------------------
      SUBROUTINE matrix_tbl_comp_dealloc(tmat)

      TYPE(tbl_comp_mat_type), INTENT(INOUT) :: tmat

      INTEGER(i4) :: iv
!-----------------------------------------------------------------------
!-PRE deallocate linear bases only.
!-----------------------------------------------------------------------
      DO iv=0,SIZE(tmat%lmat)-1
        DEALLOCATE(tmat%lmat(iv)%element)
        DEALLOCATE(tmat%lmat(iv)%from_vert)
      ENDDO
      DEALLOCATE(tmat%lmat,tmat%ix0,tmat%iy0)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_tbl_comp_dealloc
!-----------------------------------------------------------------------
!     subprogram 11. matrix_tbl_make_real
!     set matrix elements of a comp tblock matrix to their real parts.
!-----------------------------------------------------------------------
      SUBROUTINE matrix_tbl_make_real(tmat)

      TYPE(tbl_comp_mat_type), INTENT(INOUT) :: tmat

      INTEGER(i4) :: iv
!-----------------------------------------------------------------------
!-PRE using linear bases only.
!-----------------------------------------------------------------------
      DO iv=0,SIZE(tmat%lmat)-1
        tmat%lmat(iv)%element=REAL(tmat%lmat(iv)%element,r8)
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_tbl_make_real
!-----------------------------------------------------------------------
!     subprogram 12. factor_ptr_copy_rmat_info
!     log the memory of the matrix factor structure.
!-----------------------------------------------------------------------
      SUBROUTINE factor_ptr_copy_rmat_info(fac,rmat)
      TYPE(matrix_factor_type), INTENT(INOUT) :: fac
      TYPE(rbl_mat_type), INTENT(IN) :: rmat(:)
      INTEGER(i4) :: ibl,nbl

      nbl=SIZE(rmat)
      ALLOCATE(fac%mat_info(nbl))
      DO ibl=1,nbl
        fac%mat_info(ibl)%nbasis_el=rmat(ibl)%nbasis_el
        fac%mat_info(ibl)%nbtype=rmat(ibl)%nbtype
        fac%mat_info(ibl)%mx=rmat(ibl)%mx
        fac%mat_info(ibl)%my=rmat(ibl)%my
        fac%mat_info(ibl)%ix0=>rmat(ibl)%ix0
        fac%mat_info(ibl)%iy0=>rmat(ibl)%iy0
        fac%mat_info(ibl)%nb_type=>rmat(ibl)%nb_type
        fac%mat_info(ibl)%nq_type=>rmat(ibl)%nq_type
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE factor_ptr_copy_rmat_info
!-----------------------------------------------------------------------
!     subprogram 13. factor_dealloc_rmat_info
!     log the memory of the matrix factor structure.
!-----------------------------------------------------------------------
      SUBROUTINE factor_dealloc_rmat_info(fac)
      TYPE(matrix_factor_type), INTENT(INOUT) :: fac
      INTEGER(i4) :: ibl

      DEALLOCATE(fac%mat_info)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE factor_dealloc_rmat_info
!-----------------------------------------------------------------------
!     subprogram 14. factor_ptr_copy_cmat_info
!     log the memory of the matrix factor structure.
!-----------------------------------------------------------------------
      SUBROUTINE factor_ptr_copy_cmat_info(fac,cmat)
      TYPE(complex_factor_type), INTENT(INOUT) :: fac
      TYPE(rbl_comp_mat_type), INTENT(IN) :: cmat(:)
      INTEGER(i4) :: ibl,nbl

      nbl=SIZE(cmat)
      ALLOCATE(fac%mat_info(nbl))
      DO ibl=1,nbl
        fac%mat_info(ibl)%nbasis_el=cmat(ibl)%nbasis_el
        fac%mat_info(ibl)%nbtype=cmat(ibl)%nbtype
        fac%mat_info(ibl)%mx=cmat(ibl)%mx
        fac%mat_info(ibl)%my=cmat(ibl)%my
        fac%mat_info(ibl)%ix0=>cmat(ibl)%ix0
        fac%mat_info(ibl)%iy0=>cmat(ibl)%iy0
        fac%mat_info(ibl)%nb_type=>cmat(ibl)%nb_type
        fac%mat_info(ibl)%nq_type=>cmat(ibl)%nq_type
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE factor_ptr_copy_cmat_info
!-----------------------------------------------------------------------
!     subprogram 15. factor_dealloc_cmat_info
!     log the memory of the matrix factor structure.
!-----------------------------------------------------------------------
      SUBROUTINE factor_dealloc_cmat_info(fac)
      TYPE(complex_factor_type), INTENT(INOUT) :: fac
      INTEGER(i4) :: ibl

      DEALLOCATE(fac%mat_info)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE factor_dealloc_cmat_info

#ifdef OBJ_MEM_PROF
!-----------------------------------------------------------------------
!     subprogram 16. matrix_factor_real_log
!     log the memory of the matrix factor structure.
!-----------------------------------------------------------------------
      SUBROUTINE matrix_factor_real_log(fac)
      USE memlog, ONLY: memlogger
      TYPE(matrix_factor_type), INTENT(INOUT) :: fac
      INTEGER(i4) :: sz,ii,ij
      TYPE(bl_fac_type), POINTER :: blptr
      TYPE(bl_sparsity_type), POINTER :: bsptr

      sz=SIZEOF(fac)
      IF(ALLOCATED(fac%a11)) sz=sz+SIZEOF(fac%a11)
      sz=sz+sparsity_pattern_log(fac%spp)
      CALL memlogger%update(fac%mem_id,'rfac','unknown',sz)
      IF(ASSOCIATED(fac%bl_fac)) THEN
        DO ii=1,SIZE(fac%bl_fac)
          sz=SIZEOF(fac%bl_fac(ii))
          blptr=>fac%bl_fac(ii)
          IF(ASSOCIATED(blptr%lt_dg0)) sz=sz+SIZEOF(blptr%lt_dg0)
          IF(ASSOCIATED(blptr%ut_dg0)) sz=sz+SIZEOF(blptr%ut_dg0)
          IF(ASSOCIATED(blptr%lt_dg1)) sz=sz+SIZEOF(blptr%lt_dg1)
          IF(ASSOCIATED(blptr%ut_dg1)) sz=sz+SIZEOF(blptr%ut_dg1)
          IF(ASSOCIATED(blptr%xelim)) sz=sz+SIZEOF(blptr%xelim)
          IF(ASSOCIATED(blptr%yelim)) sz=sz+SIZEOF(blptr%yelim)
          IF(ASSOCIATED(blptr%xlc)) sz=sz+SIZEOF(blptr%xlc)
          IF(ASSOCIATED(blptr%xcl)) sz=sz+SIZEOF(blptr%xcl)
          IF(ASSOCIATED(blptr%ylc)) sz=sz+SIZEOF(blptr%ylc)
          IF(ASSOCIATED(blptr%ycl)) sz=sz+SIZEOF(blptr%ycl)
          IF(ASSOCIATED(blptr%lt_cen)) sz=sz+SIZEOF(blptr%lt_cen)
          IF(ASSOCIATED(blptr%lt_off)) sz=sz+SIZEOF(blptr%lt_off)
          IF(ASSOCIATED(blptr%lt_per)) sz=sz+SIZEOF(blptr%lt_per)
          IF(ASSOCIATED(blptr%ut_cen)) sz=sz+SIZEOF(blptr%ut_cen)
          IF(ASSOCIATED(blptr%ut_off)) sz=sz+SIZEOF(blptr%ut_off)
          IF(ASSOCIATED(blptr%ut_per)) sz=sz+SIZEOF(blptr%ut_per)
          IF(ASSOCIATED(blptr%adix)) sz=sz+SIZEOF(blptr%adix)
          IF(ASSOCIATED(blptr%adiy)) sz=sz+SIZEOF(blptr%adiy)
          IF(ALLOCATED(blptr%ix0)) sz=sz+SIZEOF(blptr%ix0)
          IF(ALLOCATED(blptr%iy0)) sz=sz+SIZEOF(blptr%iy0)
          IF(ALLOCATED(blptr%nb_type)) sz=sz+SIZEOF(blptr%nb_type)
          IF(ALLOCATED(blptr%nq_type)) sz=sz+SIZEOF(blptr%nq_type)
          IF(ASSOCIATED(blptr%pmat)) THEN
            DO ij=LBOUND(blptr%pmat,DIM=1),UBOUND(blptr%pmat,DIM=1)
              IF(ASSOCIATED(blptr%pmat(ij)%arr))                        &
     &           sz=sz+SIZEOF(blptr%pmat(ij)%arr)
            ENDDO
          ENDIF
          CALL memlogger%update(blptr%mem_id,'rfac_bl','unknown',sz)
        ENDDO
      ENDIF
      IF(ASSOCIATED(fac%bl_spp)) THEN
        DO ii=1,SIZE(fac%bl_spp)
          sz=SIZEOF(fac%bl_spp(ii))
          bsptr=>fac%bl_spp(ii)
          IF(ASSOCIATED(bsptr%row_ind)) THEN
            DO ij=LBOUND(bsptr%row_ind,DIM=1),                          &
                  UBOUND(bsptr%row_ind,DIM=1)
              IF(ASSOCIATED(bsptr%row_ind(ij)%rarr))                    &
     &           sz=sz+SIZEOF(bsptr%row_ind(ij)%rarr)
            ENDDO
          ENDIF
          CALL memlogger%update(bsptr%mem_id,'rfac_bl_sparse',          &
     &                                       'unknown',sz)
        ENDDO
      ENDIF
      END SUBROUTINE matrix_factor_real_log
!-----------------------------------------------------------------------
!     subprogram 17. matrix_factor_comp_log
!     log the memory of the matrix factor structure.
!-----------------------------------------------------------------------
      SUBROUTINE matrix_factor_comp_log(fac)
      USE memlog, ONLY: memlogger
      USE pardata, ONLY: node
      TYPE(complex_factor_type), INTENT(INOUT) :: fac
      INTEGER(i4) :: sz,ii,ij
      TYPE(bl_comp_fac_type), POINTER :: blptr
      TYPE(bl_sparsity_type), POINTER :: bsptr

      sz=SIZEOF(fac)
      IF(ALLOCATED(fac%a11)) sz=sz+SIZEOF(fac%a11)
      sz=sz+sparsity_pattern_log(fac%spp)
      CALL memlogger%update(fac%mem_id,'cfac','unknown',sz)
      IF(ASSOCIATED(fac%bl_fac)) THEN
        DO ii=1,SIZE(fac%bl_fac)
          sz=SIZEOF(fac%bl_fac(ii))
          blptr=>fac%bl_fac(ii)
          IF(ASSOCIATED(blptr%lt_dg0)) sz=sz+SIZEOF(blptr%lt_dg0)
          IF(ASSOCIATED(blptr%ut_dg0)) sz=sz+SIZEOF(blptr%ut_dg0)
          IF(ASSOCIATED(blptr%lt_dg1)) sz=sz+SIZEOF(blptr%lt_dg1)
          IF(ASSOCIATED(blptr%ut_dg1)) sz=sz+SIZEOF(blptr%ut_dg1)
          IF(ASSOCIATED(blptr%xelim)) sz=sz+SIZEOF(blptr%xelim)
          IF(ASSOCIATED(blptr%yelim)) sz=sz+SIZEOF(blptr%yelim)
          IF(ASSOCIATED(blptr%xlc)) sz=sz+SIZEOF(blptr%xlc)
          IF(ASSOCIATED(blptr%xcl)) sz=sz+SIZEOF(blptr%xcl)
          IF(ASSOCIATED(blptr%ylc)) sz=sz+SIZEOF(blptr%ylc)
          IF(ASSOCIATED(blptr%ycl)) sz=sz+SIZEOF(blptr%ycl)
          IF(ASSOCIATED(blptr%lt_cen)) sz=sz+SIZEOF(blptr%lt_cen)
          IF(ASSOCIATED(blptr%lt_off)) sz=sz+SIZEOF(blptr%lt_off)
          IF(ASSOCIATED(blptr%lt_per)) sz=sz+SIZEOF(blptr%lt_per)
          IF(ASSOCIATED(blptr%ut_cen)) sz=sz+SIZEOF(blptr%ut_cen)
          IF(ASSOCIATED(blptr%ut_off)) sz=sz+SIZEOF(blptr%ut_off)
          IF(ASSOCIATED(blptr%ut_per)) sz=sz+SIZEOF(blptr%ut_per)
          IF(ASSOCIATED(blptr%adix)) sz=sz+SIZEOF(blptr%adix)
          IF(ASSOCIATED(blptr%adiy)) sz=sz+SIZEOF(blptr%adiy)
          IF(ALLOCATED(blptr%ix0)) sz=sz+SIZEOF(blptr%ix0)
          IF(ALLOCATED(blptr%iy0)) sz=sz+SIZEOF(blptr%iy0)
          IF(ALLOCATED(blptr%nb_type)) sz=sz+SIZEOF(blptr%nb_type)
          IF(ALLOCATED(blptr%nq_type)) sz=sz+SIZEOF(blptr%nq_type)
          IF(ASSOCIATED(blptr%pmat)) THEN
            DO ij=LBOUND(blptr%pmat,DIM=1),UBOUND(blptr%pmat,DIM=1)
              IF(ASSOCIATED(blptr%pmat(ij)%arr))                        &
     &           sz=sz+SIZEOF(blptr%pmat(ij)%arr)
            ENDDO
          ENDIF
          CALL memlogger%update(blptr%mem_id,'cfac_bl','unknown',sz)
        ENDDO
      ENDIF
      IF(ASSOCIATED(fac%bl_spp)) THEN
        DO ii=1,SIZE(fac%bl_spp)
          sz=SIZEOF(fac%bl_spp(ii))
          bsptr=>fac%bl_spp(ii)
          IF(ASSOCIATED(bsptr%row_ind)) THEN
            DO ij=LBOUND(bsptr%row_ind,DIM=1),                          &
                  UBOUND(bsptr%row_ind,DIM=1)
              IF(ASSOCIATED(bsptr%row_ind(ij)%rarr))                    &
     &           sz=sz+SIZEOF(bsptr%row_ind(ij)%rarr)
            ENDDO
          ENDIF
          CALL memlogger%update(bsptr%mem_id,'cfac_bl_sparse',          &
     &                                       'unknown',sz)
        ENDDO
      ENDIF
      END SUBROUTINE matrix_factor_comp_log
!-----------------------------------------------------------------------
!     subprogram 18. sparsity_pattern_log.
!     log the memory of the sparsity pattern structure.
!-----------------------------------------------------------------------
      REAL(r8) FUNCTION sparsity_pattern_log(spp) RESULT(sz)
      USE memlog, ONLY: memlogger
      TYPE(sparsity_pattern), INTENT(IN) :: spp

      sz=0
      sz=sz+SIZEOF(spp)
      IF(ALLOCATED(spp%irowst_block)) sz=sz+SIZEOF(spp%irowst_block)
      IF(ASSOCIATED(spp%j_acc)) sz=sz+SIZEOF(spp%j_acc)
      IF(ASSOCIATED(spp%start_acc)) sz=sz+SIZEOF(spp%start_acc)
      IF(ASSOCIATED(spp%start_loc)) sz=sz+SIZEOF(spp%start_loc)
      IF(ALLOCATED(spp%j_acc_save)) sz=sz+SIZEOF(spp%j_acc_save)
      IF(ALLOCATED(spp%start_acc_save)) sz=sz+SIZEOF(spp%start_acc_save)
      IF(ALLOCATED(spp%recvlist)) sz=sz+SIZEOF(spp%recvlist)
      IF(ALLOCATED(spp%recvtot)) sz=sz+SIZEOF(spp%recvtot)
      IF(ALLOCATED(spp%sendlist)) sz=sz+SIZEOF(spp%sendlist)
      IF(ALLOCATED(spp%sendtot)) sz=sz+SIZEOF(spp%sendtot)
      IF(ALLOCATED(spp%sendst)) sz=sz+SIZEOF(spp%sendst)
      IF(ALLOCATED(spp%sendst)) sz=sz+SIZEOF(spp%sendrowst)
      IF(ALLOCATED(spp%recv_req)) sz=sz+SIZEOF(spp%recv_req)
      IF(ALLOCATED(spp%send_req)) sz=sz+SIZEOF(spp%send_req)
      IF(ALLOCATED(spp%recvjentry)) sz=sz+SIZEOF(spp%recvjentry)
      IF(ALLOCATED(spp%jentry)) sz=sz+SIZEOF(spp%jentry)
      IF(ALLOCATED(spp%algcount)) sz=sz+SIZEOF(spp%algcount)
      IF(ALLOCATED(spp%algdispl)) sz=sz+SIZEOF(spp%algdispl)
      IF(ALLOCATED(spp%algcntr)) sz=sz+SIZEOF(spp%algcntr)
      IF(ALLOCATED(spp%algdsplr)) sz=sz+SIZEOF(spp%algdsplr)
#ifdef HAVE_PASTIX
      IF(ASSOCIATED(spp%perm)) sz=sz+SIZEOF(spp%perm)
      IF(ASSOCIATED(spp%invp)) sz=sz+SIZEOF(spp%invp)
      IF(ASSOCIATED(spp%loc2glb)) sz=sz+SIZEOF(spp%loc2glb)
#endif /* HAVE_PASTIX */
      END FUNCTION sparsity_pattern_log
#endif /* OBJ_MEM_PROF */
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      END MODULE matrix_type_mod
