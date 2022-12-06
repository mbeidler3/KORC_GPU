!#include "config.f"
!-----------------------------------------------------------------------
!     code organization.
!-----------------------------------------------------------------------
!     1. nimfield_init
!     2. eval_field
!-----------------------------------------------------------------------
MODULE eval_field_mod
   IMPLICIT NONE
   LOGICAL :: have_b, have_j, have_v, have_n, have_t, have_p, have_f,      &
  &           have_e, have_neut, have_diff, have_fsa_beq2

CONTAINS
!-----------------------------------------------------------------------
!     subprogram 1. nimfield_init
!     A simple initialization routine for getting NIMROD data.
!-----------------------------------------------------------------------
   SUBROUTINE nimfield_init(dumpname, fieldlist)
      USE local, ONLY: i4, r8
      USE fields, ONLY: rb, nrbl, nrbl_total
      USE global, ONLY: istep, t, keff, keff_total, nmodes, nmodes_total,   &
     &                  smallnum
      USE io, ONLY: ofname, nim_wr, out_unit
      USE dump, ONLY: dump_read
      USE pardata, ONLY: nprocs, node
      USE physdat, ONLY: physdat_set, zeff, kboltz
      USE input, ONLY: dealiase, dump_file, dump_ja, nlayers, nonlinear,    &
     &                 h5io_block_proc_stride, lin_nmodes, mm_mdmap, nphi,  &
     &                 poly_degree, transfer_eq, transfer_eq_fields, geom,  &
     &                 set_phys_constants, chrg_input, zeff_input,          &
     &                 mi_input, me_input, gam_input, kblz_input, &
     &                 mu0_input, c_input, dump_fsa_beq2

      USE map_mod, ONLY: spline_guess, map_type_init
      USE lagr_quad_mod, ONLY: lagr_quad_alloc, lagr_quad_eval,             &
     &                         lagr_quad_basis_assign_loc
#ifdef HAVE_FC_HDF5
      USE mpi_nim, ONLY: comm_nimrod
      USE io, ONLY: h5in
#endif
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: dumpname, fieldlist
      INTEGER(i4) :: ierror, ns_groups, color, ii, ibl, ix, iy, ix0, iy0, ibasis
      INTEGER(i4) :: mxb, myb
      REAL(r8), DIMENSION(1) :: te, ti
      LOGICAL :: file_stat
      LOGICAL, SAVE :: open_file = .FALSE.
      CHARACTER(64) :: msg, fld
      ns_groups = 1_i4

!-----------------------------------------------------------------------
!     set output file name
!-----------------------------------------------------------------------
      ofname = 'nimfield.out'
!-----------------------------------------------------------------------
!     read the namelist file and get dump file from command-line
!-----------------------------------------------------------------------
      INQUIRE (FILE='nimrod.in', EXIST=file_stat)
      IF (.NOT. file_stat) THEN
         WRITE (nim_wr, *) 'The input file, nimrod.in, does not exist.'
      END IF

      OPEN (UNIT=out_unit, FILE=TRIM(ofname), STATUS='UNKNOWN',          &
     &     POSITION='APPEND')
      CALL read_namelist('nimrod.in', .FALSE., 'nimfl')
      IF (set_phys_constants) THEN
         CALL physdat_set(chrg_input, zeff_input, mi_input, me_input,       &
      &    gam_input, kblz_input, mu0_input, c_input) !,misp_input,zisp_input)
      ELSE
         CALL physdat_set()
      END IF
      CLOSE (out_unit)
!-----------------------------------------------------------------------
!       Adjust input
!-----------------------------------------------------------------------
      nlayers = 1  ! no layer parallelization
      h5io_block_proc_stride = nprocs/nlayers
#ifdef HAVE_FC_HDF5
      h5in%comm = comm_nimrod
#endif
!-----------------------------------------------------------------------
!     read restart dump and complete wavenumber initialization.
!-----------------------------------------------------------------------

      IF (nonlinear) THEN
         IF (dealiase < 3) THEN
            nmodes_total = nphi/2 + 1
         ELSE
            nmodes_total = nphi/dealiase + 1
         END IF
      ELSE
         nphi = 0
         nmodes_total = lin_nmodes
      END IF
      IF (nprocs == 1) THEN
         nlayers = 1
         nrbl_total = nrbl
      END IF
      INQUIRE (FILE=TRIM(dumpname), EXIST=file_stat)
      IF (.NOT. file_stat) THEN
         msg = 'File does not exist: '//TRIM(dumpname)
         CALL nim_stop(msg)
      END IF

      dump_file = dumpname
      CALL dump_read(nmodes, nmodes_total, keff, keff_total, t, istep)
      open_file = .TRUE.

!-----------------------------------------------------------------------
!     Initialize map_mod
!-----------------------------------------------------------------------
      spline_guess = .FALSE.
      IF (mm_mdmap > 0) spline_guess = .TRUE.
      DO ii = 1, 64
         fld(ii:ii) = 'n'
      END DO
      have_b = .FALSE.
      have_j = .FALSE.
      have_v = .FALSE.
      have_p = .FALSE.
      have_n = .FALSE.
      have_t = .FALSE.
      have_f = .FALSE.
      have_e = .FALSE.
      have_neut = .FALSE.
      have_diff = .FALSE.
      have_fsa_beq2 = .FALSE.

      IF (SCAN(fieldlist, "Ab", .TRUE.) /= 0) THEN
         fld(1:1) = 'y'
         have_b = .TRUE.
         IF (geom == 'tor') THEN
            DO ibl = 1, nrbl
               rb(ibl)%be_eq%fs(3, :, :) =                                    &
        &        rb(ibl)%be_eq%fs(3, :, :)/rb(ibl)%rz%fs(1, :, :)
               IF (ALLOCATED(rb(ibl)%be_eq%fsh)) THEN
                  rb(ibl)%be_eq%fsh(3, :, :, :) =                               &
         &          rb(ibl)%be_eq%fsh(3, :, :, :)/rb(ibl)%rz%fsh(1, :, :, :)
               END IF
               IF (ALLOCATED(rb(ibl)%be_eq%fsv)) THEN
                  rb(ibl)%be_eq%fsv(3, :, :, :) =                               &
         &          rb(ibl)%be_eq%fsv(3, :, :, :)/rb(ibl)%rz%fsv(1, :, :, :)
               END IF
               IF (ALLOCATED(rb(ibl)%be_eq%fsi)) THEN
                  rb(ibl)%be_eq%fsi(3, :, :, :) =                               &
         &          rb(ibl)%be_eq%fsi(3, :, :, :)/rb(ibl)%rz%fsi(1, :, :, :)
               END IF
            END DO
         END IF
      END IF
      IF (SCAN(fieldlist, "Aj", .TRUE.) /= 0) THEN
         fld(2:2) = 'y'
         have_j = .TRUE.
         IF (.NOT. dump_ja) THEN
            CALL nim_stop('Evaluation of j requires dump_ja=T')
         END IF
         IF (geom == 'tor') THEN
            DO ibl = 1, nrbl
               rb(ibl)%ja_eq%fs(3, :, :) =                                    &
        &        rb(ibl)%ja_eq%fs(3, :, :)*rb(ibl)%rz%fs(1, :, :)
               IF (ALLOCATED(rb(ibl)%ja_eq%fsh)) THEN
                  rb(ibl)%ja_eq%fsh(3, :, :, :) =                               &
         &          rb(ibl)%ja_eq%fsh(3, :, :, :)*rb(ibl)%rz%fsh(1, :, :, :)
               END IF
               IF (ALLOCATED(rb(ibl)%ja_eq%fsv)) THEN
                  rb(ibl)%ja_eq%fsv(3, :, :, :) =                               &
         &          rb(ibl)%ja_eq%fsv(3, :, :, :)*rb(ibl)%rz%fsv(1, :, :, :)
               END IF
               IF (ALLOCATED(rb(ibl)%ja_eq%fsi)) THEN
                  rb(ibl)%ja_eq%fsi(3, :, :, :) =                               &
         &          rb(ibl)%ja_eq%fsi(3, :, :, :)*rb(ibl)%rz%fsi(1, :, :, :)
               END IF
            END DO
         END IF
      END IF
      IF (SCAN(fieldlist, "Av", .TRUE.) /= 0) THEN
         fld(3:3) = 'y'
         have_v = .TRUE.
      END IF
      IF (SCAN(fieldlist, "Ap", .TRUE.) /= 0) THEN
         fld(4:4) = 'y'
         have_p = .TRUE.
      END IF
      IF (SCAN(fieldlist, "An", .TRUE.) /= 0) THEN
         fld(5:5) = 'y'
         have_n = .TRUE.
      END IF
      IF (SCAN(fieldlist, "At", .TRUE.) /= 0) THEN
         fld(6:6) = 'y'
         have_t = .TRUE.
         DO ibl = 1, nrbl
            mxb = rb(ibl)%mx
            myb = rb(ibl)%my
            CALL lagr_quad_alloc(rb(ibl)%tele_eq, mxb, myb, 1_i4,            &
       &                         poly_degree, 'teleeq', (/'teleeq'/))
            CALL lagr_quad_alloc(rb(ibl)%tion_eq, mxb, myb, 1_i4,            &
       &                         poly_degree, 'tioneq', (/'tioneq'/))
            DO ibasis = 1, SIZE(rb(ibl)%tele_eq%dx)
               ix0 = rb(ibl)%tele_eq%ix0(ibasis)
               iy0 = rb(ibl)%tele_eq%iy0(ibasis)
               DO iy = iy0, myb
                  DO ix = ix0, mxb
                     CALL lagr_quad_eval(rb(ibl)%pres_eq,                    &
          &                           ix - ix0 + rb(ibl)%tele_eq%dx(ibasis),     &
          &                           iy - iy0 + rb(ibl)%tele_eq%dy(ibasis), 0_i4)
                     CALL lagr_quad_eval(rb(ibl)%prese_eq,                   &
          &                           ix - ix0 + rb(ibl)%tele_eq%dx(ibasis),     &
          &                           iy - iy0 + rb(ibl)%tele_eq%dy(ibasis), 0_i4)
                     CALL lagr_quad_eval(rb(ibl)%nd_eq,                      &
          &                           ix - ix0 + rb(ibl)%tele_eq%dx(ibasis),     &
          &                           iy - iy0 + rb(ibl)%tele_eq%dy(ibasis), 0_i4)
                     IF (transfer_eq .AND. (INDEX(transfer_eq_fields, 'n')      &
          &                         + INDEX(transfer_eq_fields, 'N')) /= 0) THEN
                        te = 0._r8 !smallnum
                        ti = 0._r8 !smallnum
                     ELSE
                        te = MAX(smallnum,                                      &
           &                rb(ibl)%prese_eq%f(1)/(kboltz*rb(ibl)%nd_eq%f(1)))
                        ti = MAX(smallnum,                                      &
           &                 (rb(ibl)%pres_eq%f(1) - rb(ibl)%prese_eq%f(1))*    &
           &                 zeff/(kboltz*rb(ibl)%nd_eq%f(1)))
                     END IF
                     CALL lagr_quad_basis_assign_loc                         &
          &            (rb(ibl)%tele_eq, te, ibasis, ix, iy)
                     CALL lagr_quad_basis_assign_loc                         &
          &            (rb(ibl)%tion_eq, ti, ibasis, ix, iy)
                  END DO
               END DO
            END DO
         END DO
      END IF
      IF (SCAN(fieldlist, "Af", .TRUE.) /= 0) THEN
         fld(7:7) = 'y'
         have_f = .TRUE.
      END IF
      IF (SCAN(fieldlist, "Ae", .TRUE.) /= 0) THEN
         fld(8:8) = 'y'
         have_e = .TRUE.
      END IF
      IF (SCAN(fieldlist, "Ad", .TRUE.) /= 0) THEN
         fld(10:10) = 'y'
         have_diff = .TRUE.
      END IF

      IF (have_b .AND. dump_fsa_beq2) THEN
         fld(11:11) = 'y'
         have_fsa_beq2 = .TRUE.
      END IF

      CALL map_type_init(b=fld(1:1), j=fld(2:2), v=fld(3:3), p=fld(4:4),   &
     &                   pe=fld(4:4), n=fld(5:5), ti=fld(6:6), te=fld(6:6),&
     &                   f=fld(7:7), e=fld(8:8), diff_shape=fld(10:10),   &
     &                   fsa_beq2=fld(11:11))

      RETURN
   END SUBROUTINE nimfield_init
!-----------------------------------------------------------------------
!     subprogram 2. eval_field
!     A simple routine to evaluate a field at a NIMROD RZPhi location
!-----------------------------------------------------------------------
   SUBROUTINE eval_field(RZP, fname, fval, nqty, xyguess, ifail)
      USE local, ONLY: i4, r8
      USE map_mod, ONLY: rb_cel, rz_to_xy
      USE get_field_mod, ONLY: get_field
      IMPLICIT NONE

      REAL(r8), DIMENSION(1:3), INTENT(IN) :: RZP
      CHARACTER(*), INTENT(IN) :: fname
      REAL(r8), DIMENSION(1:nqty), INTENT(OUT) :: fval
      INTEGER(i4), INTENT(IN) :: nqty
      REAL(r8), DIMENSION(2), INTENT(INOUT) :: xyguess
      INTEGER(i4), INTENT(INOUT) :: ifail
      REAL(r8), DIMENSION(1:3) :: xyp
      REAL(r8), DIMENSION(1:nqty) :: pval
      REAL(r8) :: bigr, err, xc, yc
      LOGICAL :: fail
!-----------------------------------------------------------------------
!     Switch to logical coordinate, do lagrange quad eval
!-----------------------------------------------------------------------
      ifail = 0
      fail = .TRUE.
      xc = xyguess(1)
      yc = xyguess(2)
      CALL rz_to_xy(RZP(1), RZP(2), xc, yc, fail, err)
      xyguess(1) = xc; xyguess(2) = yc
      xyp(1) = xc; xyp(2) = yc; xyp(3) = RZP(3);

      IF (fname == 'b') THEN
         IF (.NOT. have_b) CALL nim_stop('b not initialized for eval')
         CALL get_field(xyp, rb_cel(1)%be_eq, rb_cel(1)%be,                &
      &                 fval, bigr, 0._r8)
      ELSEIF (fname == 'j') THEN
         IF (.NOT. have_j) CALL nim_stop('j not initialized for eval')
         CALL get_field(xyp, rb_cel(1)%ja_eq, rb_cel(1)%ja,                &
      &                 fval, bigr, 0._r8)
      ELSEIF (fname == 'v') THEN
         IF (.NOT. have_v) CALL nim_stop('v not initialized for eval')
         CALL get_field(xyp, rb_cel(1)%ve_eq, rb_cel(1)%ve,                &
      &                 fval, bigr, 0._r8)
      ELSEIF (fname == 'n') THEN
         IF (.NOT. have_n) CALL nim_stop('n not initialized for eval')
         CALL get_field(xyp, rb_cel(1)%nd_eq, rb_cel(1)%nd,                &
      &                 fval, bigr, 0._r8)
      ELSEIF (fname == 'p') THEN
         IF (.NOT. have_p) CALL nim_stop('p not initialized for eval')
         CALL get_field(xyp, rb_cel(1)%pres_eq, rb_cel(1)%pres,            &
      &                 fval, bigr, 0._r8)
      ELSEIF (fname == 'pe') THEN
         IF (.NOT. have_p) CALL nim_stop('p not initialized for eval')
         CALL get_field(xyp, rb_cel(1)%prese_eq, rb_cel(1)%prese,          &
      &                 fval, bigr, 0._r8)
      ELSEIF (fname == 'pi') THEN
         IF (.NOT. have_p) CALL nim_stop('p not initialized for eval')
         CALL get_field(xyp, rb_cel(1)%pres_eq, rb_cel(1)%pres,            &
      &                 pval, bigr, 0._r8)
         CALL get_field(xyp, rb_cel(1)%prese_eq, rb_cel(1)%prese,          &
      &                 fval, bigr, 0._r8)
         fval = pval - fval
      ELSEIF (fname == 'ti') THEN
         IF (.NOT. have_t) CALL nim_stop('t not initialized for eval')
         CALL get_field(xyp, rb_cel(1)%tion_eq, rb_cel(1)%tion,            &
      &                 fval, bigr, 0._r8)
      ELSEIF (fname == 'te') THEN
         IF (.NOT. have_t) CALL nim_stop('t not initialized for eval')
         CALL get_field(xyp, rb_cel(1)%tele_eq, rb_cel(1)%tele,            &
      &                 fval, bigr, 0._r8)
      ELSEIF (fname == 'e') THEN
         IF (.NOT. have_e) CALL nim_stop('e not initialized for eval')
         CALL get_field(xyp, rb_cel(1)%eef_eq, rb_cel(1)%eef,              &
      &                 fval, bigr, 0._r8)
      ELSEIF (fname == 'd') THEN
         IF (.NOT. have_diff)                                             &
    &      CALL nim_stop('diff_shape not initialized for eval')
         CALL get_field(xyp, rb_cel(1)%diff_shape, fval, bigr, 0._r8)
      ELSEIF (fname == 'fsa_beq2') THEN
         IF (.NOT. have_fsa_beq2)                                         &
    &      CALL nim_stop('fsa_beq2 not initialized for eval')
         CALL get_field(xyp, rb_cel(1)%fsa_beq2, fval, bigr, 0._r8)
      ELSE
         CALL nim_stop('field '//TRIM(fname)//' not found for eval')
      END IF

      RETURN
   END SUBROUTINE eval_field
!-----------------------------------------------------------------------
!     end module
!-----------------------------------------------------------------------
END MODULE eval_field_mod

!-----------------------------------------------------------------------
!   Main program
!-----------------------------------------------------------------------
!PROGRAM main
!   USE local, ONLY: i4, r8
!   USE eval_field_mod
!   IMPLICIT NONE

!   CHARACTER(64) :: dumpname, fname
!   REAL(r8), DIMENSION(1:3)  :: RZP
!   INTEGER(i4) :: ifail
!   REAL(r8), DIMENSION(2) :: xyguess
!   INTEGER(i4) , PARAMETER :: nqty=3
!   REAL(r8), DIMENSION(1:nqty) :: fval

!   dumpname = 'dump.00000.h5'
!   fname = 'b'
!   RZP(1) = 1.7
!   RZP(2) = 0.0
!   RZP(3) = 0.0
!   xyguess(1) = 1.0
!   xyguess(2) = 1.0
!   ifail = 0

!   CALL nimfield_init(dumpname, fname)
!   CALL eval_field(RZP, fname, fval, nqty, xyguess, ifail)
!   WRITE(6,*) "B field [Teslas]=",  fval

!END PROGRAM main
