!-----------------------------------------------------------------------
!     $Id: local.F90 4535 2015-08-20 04:27:38Z jking $
!     module containing defintions of real and integer kinds for linux
!     computers.
!-----------------------------------------------------------------------
#include "config.f"
      MODULE local
      IMPLICIT NONE

      INTEGER, PARAMETER ::                                             &
     &     i4=SELECTED_INT_KIND(9),                                     &
     &     i8=SELECTED_INT_KIND(18),                                    &
     &     r4=SELECTED_REAL_KIND(6,37),                                 &
     &     r8=SELECTED_REAL_KIND(13,307)

#ifdef USE_SOLVE64BIT
      INTEGER, PARAMETER :: iSolve=SELECTED_INT_KIND(18)
#else
      INTEGER, PARAMETER :: iSolve=SELECTED_INT_KIND(9)
#endif

#ifdef __Aix
      LOGICAL, PARAMETER :: rewind_namel=.true.,single_pr=.false.
#else
      LOGICAL, PARAMETER :: rewind_namel=.false.,single_pr=.false.
#endif

      INTEGER(i4), PARAMETER :: solve_nopts = 15

      END MODULE local

!-----------------------------------------------------------------------
!     System-dependent timer - returns elapsed CPU seconds
!     the 8-byte integer variables and count_tmp test are used to deal
!     with system-clock resettings and the limited maximum count.
!-----------------------------------------------------------------------
      SUBROUTINE timer(time)
      USE local

      IMPLICIT NONE

#ifdef HAVE_MPI
      INCLUDE "mpif.h"
#endif
      REAL(r8), INTENT(INOUT) :: time

      INTEGER(i4), SAVE :: count_rate,count_max,count
      INTEGER(i8), SAVE :: last_count,count_tmp
      LOGICAL, SAVE :: first_call=.true.,warned=.false.

#ifdef HAVE_MPI
!     mpi timer:
      time = MPI_Wtime()
#else
!     f90 intrinsic timer:
      IF (first_call) THEN
        CALL SYSTEM_CLOCK(count=count,count_rate=count_rate,            &
     &                    count_max=count_max)
        count_tmp=count
        first_call=.false.
      ELSE
        CALL SYSTEM_CLOCK(count=count)
        count_tmp=count
        DO WHILE (count_tmp<last_count-count_max/2)
          count_tmp=count_tmp+count_max
        ENDDO
      ENDIF
      time=REAL(count_tmp)/count_rate
      last_count=count_tmp
#endif

      RETURN
      END SUBROUTINE timer
#ifdef __Opteron
!-----------------------------------------------------------------------
!     binary open--pathf90 needs to call assign to generate standard dat
!     formats:  plot files are written as 32 bit, but dump files are
!     written as 64 bit.
!-----------------------------------------------------------------------
      SUBROUTINE open_bin(funit,fname,fstat,fpos,fbit)
      USE local

      CHARACTER(*), INTENT(IN) :: fname,fstat,fpos
      INTEGER, INTENT(IN) :: funit,fbit

      CHARACTER(256) :: assign_instr,msg
      CHARACTER(229) :: path
      CHARACTER(128) :: directory
      INTEGER(i4) :: ierror

!-----------------------------------------------------------------------
!     assign standard format to the full path name of the binary file.
!     do not call nim_stop if there is an error--this leads to an
!     infinite loop.
!-----------------------------------------------------------------------
      path=ADJUSTL(fname)
      IF (fbit==32.OR.fbit==64) THEN
!        IF (path(1:1)/='/') THEN
!          CALL system('pwd > temporary_path_name')
!          OPEN(UNIT=temp_unit,FILE='temporary_path_name')
!          READ(temp_unit,'(a)') directory
!          CLOSE(temp_unit)
!          CALL system('rm temporary_path_name')
!          path=TRIM(ADJUSTL(directory))//'/'//path
!        ENDIF
!        WRITE(assign_instr,'(2a)')'assign -F f77.mips f:',TRIM(path)
        WRITE(assign_instr,'(2a)')                                      &
     &    'assign -F f77.mips -N be f:',TRIM(path)
!        WRITE(assign_instr,'(a,i2,2a)')
!     $    'assign -F f77.mips -N ieee_',fbit,' f:',TRIM(path)
        CALL assign(TRIM(assign_instr),ierror)
        IF (ierror/=0) THEN
          WRITE(msg,'(3a)') 'Open_bin unable to assign file ',          &
     &         TRIM(ADJUSTL(fname)),'.'
          OPEN(UNIT=out_unit,FILE=TRIM(ofname),STATUS='UNKNOWN',        &
     &         POSITION='APPEND')
          WRITE(out_unit,'(2a)') 'OPEN_BIN => ', TRIM(msg)
          WRITE(nim_wr,'(2a)') 'OPEN_BIN => ', TRIM(msg)
          CLOSE(UNIT=out_unit)
          STOP
        ENDIF
      ENDIF

      OPEN(UNIT=funit,FILE=TRIM(path),STATUS=fstat,POSITION=fpos,       &
     &     FORM='UNFORMATTED')

      RETURN
      END SUBROUTINE open_bin
!-----------------------------------------------------------------------
!     binary close--pathf90 needs to call assign.
!-----------------------------------------------------------------------
      SUBROUTINE close_bin(funit,fname)
      USE local

      CHARACTER(*), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: funit

      CHARACTER(256) :: assign_instr,msg
      CHARACTER(244) :: path
      CHARACTER(128) :: directory
      INTEGER(i4) :: ierror

      CLOSE(UNIT=funit)

      path=ADJUSTL(fname)
!      IF (path(1:1)/='/') THEN
!        CALL system('pwd > temporary_path_name')
!        OPEN(UNIT=temp_unit,FILE='temporary_path_name')
!        READ(temp_unit,'(a)') directory
!        CLOSE(temp_unit)
!        CALL system('rm temporary_path_name')
!        path=TRIM(ADJUSTL(directory))//'/'//path
!      ENDIF
      WRITE(assign_instr,'(2a)') 'assign -R f:',TRIM(path)
      CALL assign(TRIM(assign_instr),ierror)

      IF (ierror/=0) THEN
        OPEN(UNIT=out_unit,FILE=TRIM(ofname),STATUS='UNKNOWN',          &
     &       POSITION='APPEND')
        WRITE(out_unit,'(3a)') 'Close_bin unable to reassign file ',    &
     &       TRIM(ADJUSTL(fname)),'!!!'
        CLOSE(UNIT=out_unit)
        WRITE(nim_wr,'(3a)') 'Close_bin unable to reassign file ',      &
     &       TRIM(ADJUSTL(fname)),'!!!'
      ENDIF

      RETURN
      END SUBROUTINE close_bin
#else
!-----------------------------------------------------------------------
!     binary open--c90 needs to call assign, Linux doesn't.
!-----------------------------------------------------------------------
      SUBROUTINE open_bin(funit,fname,fstat,fpos,fbit)
      USE local

      CHARACTER(*), INTENT(IN) :: fname,fstat,fpos
      INTEGER, INTENT(IN) :: funit,fbit

      OPEN(UNIT=funit,FILE=fname,STATUS=fstat,POSITION=fpos,            &
     &     FORM='UNFORMATTED')

      RETURN
      END SUBROUTINE open_bin
!-----------------------------------------------------------------------
!     binary close--c90 needs to call assign, Linux doesn't.
!-----------------------------------------------------------------------
      SUBROUTINE close_bin(funit,fname)
      USE local

      CHARACTER(*), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: funit

      CLOSE(UNIT=funit)

      RETURN
      END SUBROUTINE close_bin
#endif
!-----------------------------------------------------------------------
!     issue a shell command by calling a local c routine.
!-----------------------------------------------------------------------
!      SUBROUTINE system_call(command)
!      USE local
!
!      CHARACTER(*), INTENT(IN) :: command
!
!#ifdef __pathscale
!      CALL system(command)
!#else
!!#ifdef __xlf
!!      CALL nim_stop('system calls not supported with XL compiler')
!!#else
!      CHARACTER(256) :: char_to_c
!      char_to_c=TRIM(command)
!      CALL local_system(char_to_c)
!!#endif
!#endif
!
!      RETURN
!      END SUBROUTINE system_call
!-----------------------------------------------------------------------
!     Get the command line arguments
!-----------------------------------------------------------------------
      SUBROUTINE get_arg(i,arg)
      USE local

      INTEGER(i4), INTENT(IN) :: i
      CHARACTER(*), INTENT(OUT) :: arg

      ! Fortran2003 standard
      CALL get_command_argument(i, arg)

      ! Common extension
      !call getarg(i,arg)

      RETURN
      END SUBROUTINE get_arg
!-----------------------------------------------------------------------
!     Get the number of command line arguments
!-----------------------------------------------------------------------
      SUBROUTINE get_arg_count(numargs)
      USE local

      INTEGER(i4), INTENT(OUT) :: numargs

       ! Fortran2003 standard
       numargs = command_argument_count()

       ! Common extension
       !numargs = iargc()

      RETURN
      END SUBROUTINE get_arg_count
!-----------------------------------------------------------------------
!     Change system dir
!     This works for every compiler but gfortran and xlf
!-----------------------------------------------------------------------
      SUBROUTINE sys_chdir(path, ierror)
      USE local
      IMPLICIT none
      CHARACTER(len=*),INTENT(IN)  :: path ! chdir to this dir
      INTEGER(i4),INTENT(OUT) :: ierror    ! return code
      INTEGER(i4)             :: lenpath   ! length of path
#ifndef __gfortran
      INTEGER, EXTERNAL    :: chdir
#endif
!-----------------------------------------------------------------------
      lenpath=len_trim(path)
#ifndef __xlf
#ifdef __gfortran
      CALL chdir(path(1:lenpath))
#else
      ierror=chdir(path(1:lenpath))
#endif
#endif
      END SUBROUTINE sys_chdir
