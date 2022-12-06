!     $Id: time.f90 5640 2018-02-08 00:53:05Z charlson $
! routines for timing statistics
! storage for accumulating incremental times

      MODULE time
      USE local
      IMPLICIT NONE

! cummulative times and counters

      INTEGER(i4) :: seamcount = 0
      INTEGER(i4) :: segcount = 0
      INTEGER(i4) :: seamcount_loop,segcount_loop

      REAL(r8) :: time_io = 0.0
      REAL(r8) :: time_seam = 0.0
      REAL(r8) :: time_seg = 0.0
      REAL(r8) :: time_iter = 0.0
      REAL(r8) :: time_fac = 0.0
      REAL(r8) :: time_slu = 0.0
      REAL(r8) :: time_mumps = 0.0
      REAL(r8) :: time_hypre = 0.0
      REAL(r8) :: time_pardiso = 0.0
      REAL(r8) :: time_fft = 0.0
      REAL(r8) :: time_fftcm = 0.0
      REAL(r8) :: time_mat = 0.0
      REAL(r8) :: time_rhs = 0.0
      REAL(r8) :: time_line = 0.0
      REAL(r8) :: time_stcon = 0.0
      REAL(r8) :: time_cel = 0.0
      REAL(r8) :: time_part = 0.0
      REAL(r8) :: time_sep = 0.0
      REAL(r8) :: time_sol = 0.0
      REAL(r8) :: time_imp_rad = 0.0
      REAL(r8) :: time_imp_lsd = 0.0
      REAL(r8) :: time_imp_SPI = 0.0
      REAL(r8) :: time_io_loop,time_seam_loop,time_seg_loop
      REAL(r8) :: time_iter_loop,time_fft_loop,time_mat_loop
      REAL(r8) :: time_fac_loop,time_line_loop,time_stcon_loop
      REAL(r8) :: time_rhs_loop,time_slu_loop,time_fftcm_loop
      REAL(r8) :: time_mumps_loop,time_sep_loop,time_sol_loop
      REAL(r8) :: time_hypre_loop,time_part_loop,time_cel_loop
      REAL(r8) :: time_pardiso_loop
      REAL(r8) :: time_imp_rad_loop,time_imp_lsd_loop
      REAL(r8) :: time_imp_SPI_loop

      REAL(r8) :: time_total_start,time_total_end
      REAL(r8) :: time_loop_start,time_loop_end
      REAL(r8) :: timestart,timeend

      CONTAINS

! initialize cummulative values

      SUBROUTINE timer_init
      IMPLICIT NONE

      seamcount = 0
      segcount = 0
      time_io = 0.0
      time_seam = 0.0
      time_seg = 0.0
      time_iter = 0.0
      time_fac = 0.0
      time_slu = 0.0
      time_mumps = 0.0
      time_hypre = 0.0
      time_pardiso = 0.0
      time_fft = 0.0
      time_fftcm = 0.0
      time_mat = 0.0
      time_rhs = 0.0
      time_line = 0.0
      time_stcon = 0.0
      time_cel = 0.0
      time_part = 0.0
      time_sep = 0.0
      time_sol = 0.0
      time_imp_rad = 0.0
      time_imp_lsd = 0.0
      time_imp_SPI = 0.0

      RETURN
      END SUBROUTINE timer_init

! store io & seam times/counters from within main timestepping loop

      SUBROUTINE timer_close
      IMPLICIT NONE

      seamcount_loop = seamcount
      segcount_loop = segcount
      time_io_loop = time_io
      time_seam_loop = time_seam
      time_seg_loop = time_seg
      time_iter_loop = time_iter
      time_fac_loop = time_fac
      time_slu_loop = time_slu
      time_mumps_loop = time_mumps
      time_hypre_loop = time_hypre
      time_pardiso_loop = time_pardiso
      time_fft_loop = time_fft
      time_fftcm_loop = time_fftcm
      time_mat_loop = time_mat
      time_rhs_loop = time_rhs
      time_line_loop = time_line
      time_stcon_loop = time_stcon
      time_cel_loop = time_cel
      time_part_loop = time_part
      time_sep_loop = time_sep
      time_sol_loop = time_sol
      time_imp_rad_loop = time_imp_rad
      time_imp_lsd_loop = time_imp_lsd
      time_imp_SPI_loop = time_imp_SPI

      RETURN
      END SUBROUTINE timer_close

! print-out times, averaged across procs

      SUBROUTINE timer_outs(thetime,thestep,txtflg,binflg)
      USE mpi_nim
      USE pardata
      USE io
      IMPLICIT NONE
      REAL, INTENT(IN) :: thetime
      INTEGER, INTENT(IN) :: thestep
      LOGICAL, INTENT(IN) :: txtflg,binflg
      REAL(r8) :: tmp,time_setup,therealstep
      REAL(r8) :: time_total,time_loop
      INTEGER(i4) :: ierror,uwrite,iwr
      CHARACTER(6), SAVE :: hist_pos='REWIND'

      CALL timer_close

! average cummulative timers across all procs

      IF (nprocs>1) THEN
        CALL mpi_allreduce(time_io_loop,tmp,1,mpi_nim_real,mpi_sum,     &
     &       comm_nimrod,ierror)
        time_io_loop = tmp/nprocs
        CALL mpi_allreduce(time_seam_loop,tmp,1,mpi_nim_real,mpi_sum,   &
     &       comm_nimrod,ierror)
        time_seam_loop = tmp/nprocs
        CALL mpi_allreduce(time_seg_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_seg_loop = tmp/nprocs
        CALL mpi_allreduce(time_iter_loop,tmp,1,mpi_nim_real,mpi_sum,   &
     &       comm_nimrod,ierror)
        time_iter_loop = tmp/nprocs
        CALL mpi_allreduce(time_fac_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_fac_loop = tmp/nprocs
        CALL mpi_allreduce(time_slu_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_slu_loop = tmp/nprocs
        CALL mpi_allreduce(time_mumps_loop,tmp,1,mpi_nim_real,mpi_sum,  &
     &       comm_nimrod,ierror)
        time_mumps_loop = tmp/nprocs
        CALL mpi_allreduce(time_hypre_loop,tmp,1,mpi_nim_real,mpi_sum,  &
     &       comm_nimrod,ierror)
        time_hypre_loop = tmp/nprocs
        CALL mpi_allreduce(time_pardiso_loop,tmp,1,mpi_nim_real,mpi_sum,&
     &       comm_nimrod,ierror)
        time_pardiso_loop = tmp/nprocs
        CALL mpi_allreduce(time_fft_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_fft_loop = tmp/nprocs
        CALL mpi_allreduce(time_fftcm_loop,tmp,1,mpi_nim_real,mpi_sum,  &
     &       comm_nimrod,ierror)
        time_fftcm_loop = tmp/nprocs
        CALL mpi_allreduce(time_mat_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_mat_loop = tmp/nprocs
        CALL mpi_allreduce(time_rhs_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_rhs_loop = tmp/nprocs
        CALL mpi_allreduce(time_line_loop,tmp,1,mpi_nim_real,mpi_sum,   &
     &       comm_nimrod,ierror)
        time_line_loop = tmp/nprocs
        CALL mpi_allreduce(time_stcon_loop,tmp,1,mpi_nim_real,mpi_sum,  &
     &       comm_nimrod,ierror)
        time_stcon_loop = tmp/nprocs
        CALL mpi_allreduce(time_cel_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_cel_loop = tmp/nprocs
        CALL mpi_allreduce(time_part_loop,tmp,1,mpi_nim_real,mpi_sum,   &
     &       comm_nimrod,ierror)
        time_part_loop = tmp/nprocs
        CALL mpi_allreduce(time_sep_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_sep_loop = tmp/nprocs
        CALL mpi_allreduce(time_sol_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_sol_loop = tmp/nprocs
        CALL mpi_allreduce(time_imp_rad_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_imp_rad_loop = tmp/nprocs
        CALL mpi_allreduce(time_imp_lsd_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_imp_lsd_loop = tmp/nprocs
        CALL mpi_allreduce(time_imp_SPI_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_imp_SPI_loop = tmp/nprocs
      ENDIF

! print-out

   15 format(a,es12.5)

      if (node == 0) then

        IF (txtflg) THEN
        OPEN(UNIT=out_unit,FILE=TRIM(ofname),STATUS='UNKNOWN',          &
     &       POSITION='APPEND')

        do iwr=1,2

          if (iwr == 1) then
            uwrite=nim_wr
          else
            uwrite=out_unit
          endif

          write (uwrite,*)
          write (uwrite,*)                                              &
     &       'Timing statistics:'
          if (time_seam_loop>0)                                         &
     &      write (uwrite,15)'    Seam  time      = ',time_seam_loop
          if (time_seg_loop>0)                                          &
     &      write (uwrite,15)'    Seg   time      = ',time_seg_loop
          if (time_io_loop>0)                                           &
     &      write (uwrite,15)'    I/O   time      = ',time_io_loop
          if (time_iter_loop>0)                                         &
     &      write (uwrite,15)'    Iteration time  = ',time_iter_loop
          if (time_fac_loop>0)                                          &
     &      write (uwrite,15)'    Factoring time  = ',time_fac_loop
          if (time_slu_loop>0)                                          &
     &      write (uwrite,15)'    SuperLU time    = ',time_slu_loop
          if (time_mumps_loop>0)                                        &
     &      write (uwrite,15)'    Mumps time      = ',time_mumps_loop
          if (time_hypre_loop>0)                                        &
     &      write (uwrite,15)'    Hypre time      = ',time_hypre_loop
          if (time_pardiso_loop>0)                                      &
     &      write (uwrite,15)'    Pardiso time    = ',time_pardiso_loop
          if (time_line_loop>0)                                         &
     &      write (uwrite,15)'    Line_comm time  = ',time_line_loop
          if (time_fft_loop>0)                                          &
     &      write (uwrite,15)'    FFT   time      = ',time_fft_loop
          if (time_fftcm_loop>0)                                        &
     &      write (uwrite,15)'    FFT_comm time   = ',time_fftcm_loop
          if (time_mat_loop>0)                                          &
     &      write (uwrite,15)'    FE_matrix time  = ',time_mat_loop
          if (time_rhs_loop>0)                                          &
     &      write (uwrite,15)'    FE_vector time  = ',time_rhs_loop
          if (time_stcon_loop>0)                                        &
     &      write (uwrite,15)'    Static_con time = ',time_stcon_loop
          if (time_cel_loop>0)                                          &
     &     write (uwrite,15)'     CEL_comp time   = ',time_cel_loop
          if (time_part_loop>0)                                         &
     &      write (uwrite,15)'    Part_push time  = ',time_part_loop
          if (time_sep_loop>0)                                          &
     &      write (uwrite,15)'    Contour time    = ',time_sep_loop
          if (time_sol_loop>0)                                          &
     &      write (uwrite,15)'    Sol_find time   = ',time_sol_loop
          if (time_imp_rad_loop>0)                                      &
     &      write (uwrite,15)'    Imp_rad time    = ',time_imp_rad_loop
          if (time_imp_lsd_loop>0)                                      &
     &      write (uwrite,15)'    Imp_lsode time  = ',time_imp_lsd_loop
          if (time_imp_SPI_loop>0)                                      &
     &      write (uwrite,15)'    Imp_SPI time    = ',time_imp_SPI_loop

        enddo

        CLOSE(UNIT=out_unit)
        endif

!        if (binflg) then
!          CALL open_bin(xdr_unit,'timestat.bin','UNKNOWN',hist_pos,32_i4)
!          hist_pos='APPEND'
!          therealstep=REAL(thestep)
!          WRITE(xdr_unit) REAL( (/ therealstep , thetime ,              &
!     &       time_io_loop , time_seam_loop , time_seg_loop ,            &
!     &       time_iter_loop , time_fac_loop , time_slu_loop ,           &
!     &       time_mumps_loop , time_hypre_loop , time_pardiso_loop ,    &
!     &       time_fft_loop , time_fftcm_loop , time_mat_loop ,          &
!     &       time_rhs_loop , time_line_loop , time_stcon_loop ,         &
!     &       time_cel_loop , time_part_loop , time_sep_loop ,           &
!     &       time_sol_loop ,                                            &
!     &       time_imp_rad_loop , time_imp_lsd_loop , time_imp_SPI_loop  &
!            & /) , 4 )
!          WRITE(xdr_unit)
!          CALL close_bin(xdr_unit,'timestat.bin' )
!        endif 
      endif


      RETURN
      END SUBROUTINE timer_outs

! print-out stats, averaged across procs

      SUBROUTINE timer_stats
      USE mpi_nim
      USE pardata
      USE io
      IMPLICIT NONE
      REAL(r8) :: tmp,time_setup
      REAL(r8) :: time_total,time_loop
      INTEGER(i4) :: ierror,uwrite,iwr

! average cummulative timers across all procs

      IF (nprocs>1) THEN
        CALL mpi_allreduce(time_io_loop,tmp,1,mpi_nim_real,mpi_sum,     &
     &       comm_nimrod,ierror)
        time_io_loop = tmp/nprocs
        CALL mpi_allreduce(time_seam_loop,tmp,1,mpi_nim_real,mpi_sum,   &
     &       comm_nimrod,ierror)
        time_seam_loop = tmp/nprocs
        CALL mpi_allreduce(time_seg_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_seg_loop = tmp/nprocs
        CALL mpi_allreduce(time_iter_loop,tmp,1,mpi_nim_real,mpi_sum,   &
     &       comm_nimrod,ierror)
        time_iter_loop = tmp/nprocs
        CALL mpi_allreduce(time_fac_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_fac_loop = tmp/nprocs
        CALL mpi_allreduce(time_slu_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_slu_loop = tmp/nprocs
        CALL mpi_allreduce(time_mumps_loop,tmp,1,mpi_nim_real,mpi_sum,  &
     &       comm_nimrod,ierror)
        time_mumps_loop = tmp/nprocs
        CALL mpi_allreduce(time_hypre_loop,tmp,1,mpi_nim_real,mpi_sum,  &
     &       comm_nimrod,ierror)
        time_hypre_loop = tmp/nprocs
        CALL mpi_allreduce(time_pardiso_loop,tmp,1,mpi_nim_real,mpi_sum,&
     &       comm_nimrod,ierror)
        time_pardiso_loop = tmp/nprocs
        CALL mpi_allreduce(time_fft_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_fft_loop = tmp/nprocs
        CALL mpi_allreduce(time_fftcm_loop,tmp,1,mpi_nim_real,mpi_sum,  &
     &       comm_nimrod,ierror)
        time_fftcm_loop = tmp/nprocs
        CALL mpi_allreduce(time_mat_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_mat_loop = tmp/nprocs
        CALL mpi_allreduce(time_rhs_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_rhs_loop = tmp/nprocs
        CALL mpi_allreduce(time_line_loop,tmp,1,mpi_nim_real,mpi_sum,   &
     &       comm_nimrod,ierror)
        time_line_loop = tmp/nprocs
        CALL mpi_allreduce(time_stcon_loop,tmp,1,mpi_nim_real,mpi_sum,  &
     &       comm_nimrod,ierror)
        time_stcon_loop = tmp/nprocs
        CALL mpi_allreduce(time_cel_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_cel_loop = tmp/nprocs
        CALL mpi_allreduce(time_part_loop,tmp,1,mpi_nim_real,mpi_sum,   &
     &       comm_nimrod,ierror)
        time_part_loop = tmp/nprocs
        CALL mpi_allreduce(time_sep_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_sep_loop = tmp/nprocs
        CALL mpi_allreduce(time_sol_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_sol_loop = tmp/nprocs
        CALL mpi_allreduce(time_imp_rad_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_imp_rad_loop = tmp/nprocs
        CALL mpi_allreduce(time_imp_lsd_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_imp_lsd_loop = tmp/nprocs
        CALL mpi_allreduce(time_imp_SPI_loop,tmp,1,mpi_nim_real,mpi_sum,    &
     &       comm_nimrod,ierror)
        time_imp_SPI_loop = tmp/nprocs
      ENDIF

! compute derivative times

      time_total = MAX(time_total_end - time_total_start,               &
     &             SQRT(TINY(time_total)))
      time_loop = MAX(time_loop_end - time_loop_start,                  &
     &             SQRT(TINY(time_total)))
      time_setup = time_total - time_loop

! print-out

   10 format(a,2es12.5)

      if (node == 0) then

        OPEN(UNIT=out_unit,FILE=TRIM(ofname),STATUS='UNKNOWN',          &
     &       POSITION='APPEND')

        do iwr=1,2

          if (iwr == 1) then
            uwrite=nim_wr
          else
            uwrite=out_unit
          endif

          write (uwrite,*)
          write (uwrite,*)                                              &
     &       'Timing statistics:'
          write (uwrite,*)                                              &
     &       '                        CPU secs    % of loop'
          if (time_seam_loop>0)                                         &
     &      write (uwrite,10)'    Seam  time      = ',time_seam_loop,   &
     &         time_seam_loop/time_loop*100
          if (time_seg_loop>0)                                          &
     &      write (uwrite,10)'    Seg   time      = ',time_seg_loop,    &
     &         time_seg_loop/time_loop*100
          if (time_io_loop>0)                                           &
     &      write (uwrite,10)'    I/O   time      = ',time_io_loop,     &
     &         time_io_loop/time_loop*100
          if (time_iter_loop>0)                                         &
     &      write (uwrite,10)'    Iteration time  = ',time_iter_loop,   &
     &         time_iter_loop/time_loop*100
          if (time_fac_loop>0)                                          &
     &      write (uwrite,10)'    Factoring time  = ',time_fac_loop,    &
     &         time_fac_loop/time_loop*100
          if (time_slu_loop>0)                                          &
     &      write (uwrite,10)'    SuperLU time    = ',time_slu_loop,    &
     &         time_slu_loop/time_loop*100
          if (time_mumps_loop>0)                                        &
     &      write (uwrite,10)'    Mumps time      = ',time_mumps_loop,  &
     &         time_mumps_loop/time_loop*100
          if (time_hypre_loop>0)                                        &
     &      write (uwrite,10)'    Hypre time      = ',time_hypre_loop,  &
     &         time_hypre_loop/time_loop*100
          if (time_pardiso_loop>0)                                      &
     &      write (uwrite,10)'    Pardiso time    = ',                  &
     &         time_pardiso_loop,time_pardiso_loop/time_loop*100
          if (time_line_loop>0)                                         &
     &      write (uwrite,10)'    Line_comm time  = ',time_line_loop,   &
     &         time_line_loop/time_loop*100
          if (time_fft_loop>0)                                          &
     &      write (uwrite,10)'    FFT   time      = ',time_fft_loop,    &
     &         time_fft_loop/time_loop*100
          if (time_fftcm_loop>0)                                        &
     &      write (uwrite,10)'    FFT_comm time   = ',time_fftcm_loop,  &
     &         time_fftcm_loop/time_loop*100
          if (time_mat_loop>0)                                          &
     &      write (uwrite,10)'    FE_matrix time  = ',time_mat_loop,    &
     &         time_mat_loop/time_loop*100
          if (time_rhs_loop>0)                                          &
     &      write (uwrite,10)'    FE_vector time  = ',time_rhs_loop,    &
     &         time_rhs_loop/time_loop*100
          if (time_stcon_loop>0)                                        &
     &      write (uwrite,10)'    Static_con time = ',time_stcon_loop,  &
     &         time_stcon_loop/time_loop*100
          if (time_cel_loop>0)                                          &
     &     write (uwrite,10)'    CEL_comp time    = ',time_cel_loop,    &
     &         time_cel_loop/time_loop*100
          if (time_part_loop>0)                                         &
     &      write (uwrite,10)'    Part_push time  = ',time_part_loop,   &
     &         time_part_loop/time_loop*100
          if (time_sep_loop>0)                                          &
     &      write (uwrite,10)'    Contour time    = ',time_sep_loop,    &
     &         time_sep_loop/time_loop*100
          if (time_sol_loop>0)                                          &
     &      write (uwrite,10)'    Sol_find time   = ',time_sol_loop,    &
     &         time_sol_loop/time_loop*100
          if (time_imp_rad_loop>0)                                      &
     &      write (uwrite,10)'    Imp_rad time    = ',time_imp_rad_loop,&
     &         time_imp_rad_loop/time_loop*100
          if (time_imp_lsd_loop>0)                                      &
     &      write (uwrite,10)'    Imp_lsode time  = ',time_imp_lsd_loop,&
     &         time_imp_lsd_loop/time_loop*100
          if (time_imp_SPI_loop>0)                                      &
     &      write (uwrite,10)'    Imp_SPI time    = ',time_imp_SPI_loop,&
     &         time_imp_SPI_loop/time_loop*100
          write (uwrite,*)'    # of Seamings in loop =',seamcount_loop
          write (uwrite,*)'    # of SegSeams in loop =',segcount_loop
          write (uwrite,*)                                              &
     &       '                  CPU secs   % of total'
          write (uwrite,10)'    Loop  time      = ',time_loop,          &
     &       time_loop/time_total*100
          write (uwrite,10)'    Setup time      = ',time_setup,         &
     &       time_setup/time_total*100
          write (uwrite,10)'    Total time      = ',time_total,100._r8

        enddo

        CLOSE(UNIT=out_unit)

      endif

      RETURN
      END SUBROUTINE timer_stats


      END MODULE time
