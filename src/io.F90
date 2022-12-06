!-----------------------------------------------------------------------
!     $Id: io.F90 6534 2019-04-03 21:45:39Z jking $
!     module containing fortran unit numbers for io.
!-----------------------------------------------------------------------
#include "config.f"
      MODULE io
#ifdef HAVE_FC_HDF5
      USE hdf5_api
#endif /* HAVE_FC_HDF5 */
      IMPLICIT NONE

      CHARACTER(64) :: ofname="none"         !< output file name

      INTEGER, PARAMETER :: in_unit = 1      !< input file unit
      INTEGER, PARAMETER :: out_unit = 2     !< output file unit
      INTEGER, PARAMETER :: xdr_unit = 3     !< unit for xdraw output
      INTEGER, PARAMETER :: nim_rd = 5       !< stdin unit for reads
      INTEGER, PARAMETER :: nim_wr = 6       !< stdout unit for writes
      INTEGER, PARAMETER :: dx1_unit = 8     !< nimrod.dx unit
      INTEGER, PARAMETER :: dx2_unit = 9     !< nimrod.bin unit
      INTEGER, PARAMETER :: txt_unit = 10    !< unit for .txt or .dat output
      INTEGER, PARAMETER :: hst_unit = 11    !< unit for nimhist.bin
      INTEGER, PARAMETER :: it_unit = 12     !< unit for iter.out
      INTEGER, PARAMETER :: en_unit = 13     !< unit for energy.bin
      INTEGER, PARAMETER :: dis_unit = 14    !< unit for discharge.bin
      INTEGER, PARAMETER :: eq_unit = 16     !< equilibrium file
      INTEGER, PARAMETER :: grd_unit = 17    !< unit for grid check
      INTEGER, PARAMETER :: tec2d = 18       !< tecplot output unit
      INTEGER, PARAMETER :: tec1d=19
      INTEGER, PARAMETER :: iota_unit=20     !< iotabar file unit
      INTEGER, PARAMETER :: dump_unit = 25   !< dump file unit
      INTEGER, PARAMETER :: rstrt_unit = 26  !< restart file unit
      INTEGER, PARAMETER :: pie_unit = 27    !< pie file unit
      INTEGER, PARAMETER :: rim_unit = 28    !< rim file unit
      INTEGER, PARAMETER :: xy_unit = 31     !< xy binary output
      INTEGER, PARAMETER :: xt_unit = 32     !< xt binary output
      INTEGER, PARAMETER :: yt_unit = 33     !< yt binary output
      INTEGER, PARAMETER :: grid_unit = 34   !< grid binary output
      INTEGER, PARAMETER :: con_unit = 35    !< contour plot output
      INTEGER, PARAMETER :: ascii_unit = 41  !< ascii diagnostic file
      INTEGER, PARAMETER :: binary_unit = 42 !< binary diagnostic file
      INTEGER, PARAMETER :: eq1_unit=43
      INTEGER, PARAMETER :: eq2_unit=44
      INTEGER, PARAMETER :: pack_ascii_unit=45
      INTEGER, PARAMETER :: pack_bin_unit=46
      INTEGER, PARAMETER :: pack_in_unit=47
      INTEGER, PARAMETER :: pack_out_unit=48
      INTEGER, PARAMETER :: pack_spline_unit=49
      INTEGER, PARAMETER :: dcut_unit = 51
      INTEGER, PARAMETER :: memlog_unit = 52 !< memory logger output
      INTEGER, PARAMETER :: prb_in_unit = 54
      INTEGER, PARAMETER :: prb_out_unit = 55
      INTEGER, PARAMETER :: temp_unit = 91   !< temporary files
      INTEGER, PARAMETER :: xycel_unit = 92  !< xy binary output
      INTEGER, PARAMETER :: xtcel_unit = 93  !< xt binary output
      INTEGER, PARAMETER :: aux_eq_unit = 94 !< aux_eq_write file

#ifdef HAVE_FC_HDF5
      TYPE(hdf5ErrorType) :: h5err
      TYPE(hdf5InOpts), SAVE :: h5in
      INTEGER(HID_T) :: fileid,rootgid

      CONTAINS

!-----------------------------------------------------------------------
!     subprogram 1. fch5init.
!     initialize fcio h5 information for writes and reads.
!-----------------------------------------------------------------------
      SUBROUTINE fch5init
      USE mpi_nim
      IMPLICIT NONE
      LOGICAL, SAVE :: h5init=.FALSE.

      IF (.NOT.h5init) THEN
        CALL vshdf5_fcinit
        CALL vshdf5_inith5vars(h5in, h5err)
        h5in%comm=comm_stride_layer
        h5in%info=mpi_info_null
        ! h5in%data_xfer_mode=H5FD_MPIO_COLLECTIVE_F does not work
        h5in%data_xfer_mode=H5FD_MPIO_INDEPENDENT_F
        h5in%verbose=.FALSE.
        h5in%debug=.FALSE.
        h5init=.TRUE.
      ENDIF
      END SUBROUTINE fch5init

!-----------------------------------------------------------------------
!     subprogram 2. fch5resetvars.
!     reset h5 information for writes and reads.
!-----------------------------------------------------------------------
      SUBROUTINE fch5resetvars
      IMPLICIT NONE

      h5in%mesh =  " "
      h5in%vsAxisLabels = " "
      h5in%units =  " "
      h5in%vsCentering =  " "
      h5in%vsMD =  " "
      h5in%vsIndexOrder = " "
      h5in%vsLabels = " "

      END SUBROUTINE fch5resetvars

#endif /* HAVE_FC_HDF5 */

      END MODULE io
