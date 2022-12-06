!-----------------------------------------------------------------------
!  $Id: memlog.F90 4517 2015-07-29 23:08:48Z jking $
!> @file memlog.F90
!! @brief routines for logging the use of system heap memory.
!-----------------------------------------------------------------------
#include "config.f"
!-----------------------------------------------------------------------
!> @brief this module includes the singleton object used to register
!!  heap memory allocation and deallocation by NIMROD objects.  
!! @details this library is meant to be a simple interface to logging
!!  heap allocation by NIMROD objects. there are two essential 
!!  routines: report and update. Object allocation and deallocation
!!  functions should call the update routine. The report statement
!!  can be placed anywhere in the code to report the memory usage.
!!  Example report:
!!  > USE memlog, ONLY: memlogger
!!  > CALL memlogger%report(1,istep,'time_step_loop')
!!  Example update (allocate):
!!  > USE memlog, only: memlogger
!!  > CALL memlogger%update(rvt%mem_id,'rvec','unknown',SIZEOF(arri))
!!  NOTE: mem_id is a member of the type, it is an integer pointer
!!  used to track the object internally in memlog.
!!  Example update (deallocate):
!!  > USE memlog, only: memlogger
!!  > CALL memlogger%update(rvt%mem_id,' ',' ',0,resize=.TRUE.)
!-----------------------------------------------------------------------
MODULE memlog
  USE local, ONLY: i4,r8
  IMPLICIT NONE
  PRIVATE

  INTEGER(i4), PARAMETER :: num_obj_max=1e5_i4
  INTEGER(i4), PARAMETER :: strlen = 16
  CHARACTER(16), SAVE :: rootdir = 'none'

!-----------------------------------------------------------------------
!> @brief the empty_list is a linked list of the unallocated entries
!!  of mem_table.
!-----------------------------------------------------------------------
  TYPE empty_list
    PRIVATE
    INTEGER(i4) :: mem_id
    TYPE(empty_list), POINTER :: next
    TYPE(empty_list), POINTER :: last
  END TYPE empty_list

!-----------------------------------------------------------------------
!> @brief the mem_table is a table of the memory allocated to 
!!  registered objects. 
!! @details mem_table contains tables of the size, type, description 
!!  and ids of register objects, as well as a linked list (first_empty)
!!  that indicates empties in the tables where objects have been
!!  allocated and then subsequently deallocated.
!-----------------------------------------------------------------------
  TYPE mem_table
    PRIVATE
    REAL(r8), DIMENSION(:), ALLOCATABLE :: obj_size
    INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: obj_list ! linked list
    CHARACTER(len=strlen), DIMENSION(:,:), ALLOCATABLE :: obj_desc 
    INTEGER(i4) :: current_size=0_i4
    INTEGER(i4) :: current_fill=0_i4
    TYPE(empty_list), POINTER :: first_empty
    INTEGER(i4), POINTER :: my_mem_id=>NULL()
  CONTAINS
    PROCEDURE, PASS :: add_entry
    PROCEDURE, PASS :: modify_entry
    PROCEDURE, PASS :: delete_entry
    PROCEDURE, PASS :: sort_list
    ! Final keyword is not yet supported with gfortran. The free 
    ! routine of mem_logger can be eliminated when it is.
    !FINAL :: free_list
    PROCEDURE, PASS :: free_list
  END TYPE mem_table

!-----------------------------------------------------------------------
!> @brief the mem_logger type stores, provides an update interface 
!!  to, and reports data from the mem_table table.
!! @details the mem_logger contains the user facing routines update
!!  and report. the former increments and decrements the memory 
!!  counters, while the latter is used to report the held statistics.
!-----------------------------------------------------------------------
  TYPE mem_logger
    PRIVATE
    TYPE(mem_table) :: memtable
    INTEGER(i4) :: node=-1_i4
    INTEGER(i4) :: nprocs=1_i4
    ! Pointers to the total and resident virtual memory pulled from a 
    ! system call.
    INTEGER(i4), POINTER :: vsiz_mem_id=>NULL()
    INTEGER(i4), POINTER :: vrss_mem_id=>NULL()
  CONTAINS
    PROCEDURE, PASS :: update
    PROCEDURE, PASS :: report
    !FINAL :: finalize - can be eliminated when FINAL is supported.
    PROCEDURE, PASS :: finalize
  END TYPE mem_logger

!-----------------------------------------------------------------------
!> @brief this is the single instance of the memlogger.
!-----------------------------------------------------------------------
  TYPE(mem_logger), SAVE, PUBLIC :: memlogger

CONTAINS

!-----------------------------------------------------------------------
!> @brief add an entry to the linked list, and return the new id.
!! @param[inout] mem_id 
!! @param[in] obj_type
!! @param[in] obj_name
!! @payyram[in] mem_size
!-----------------------------------------------------------------------
  SUBROUTINE add_entry(this,mem_id,obj_type,obj_name,mem_size)
    CLASS(mem_table), INTENT(INOUT), TARGET :: this
    INTEGER(i4), POINTER, INTENT(INOUT) :: mem_id
    CHARACTER(len=*), INTENT(IN) :: obj_type
    CHARACTER(len=*), INTENT(IN) :: obj_name
    REAL(r8), INTENT(IN) :: mem_size
    
    INTEGER(i4) :: slen,ii
    LOGICAL :: update_logmem

    !-------------------------------------------------------------------
    ! resize arrays if needed. minimize the times the code encounters
    ! the critical section.
    !-------------------------------------------------------------------
    IF (this%current_size==0_i4) THEN
      !$omp critical
      IF (this%current_size==0_i4) THEN
        this%current_size=num_obj_max
        ALLOCATE(this%obj_size(num_obj_max))
        this%obj_size(:)=0_i4
        ALLOCATE(this%obj_list(2_i4,num_obj_max))
        DO ii=1,this%current_size
          this%obj_list(1:2,ii)=ii
        ENDDO
        ALLOCATE(this%obj_desc(2_i4,num_obj_max))
        this%obj_desc=''
        ALLOCATE(this%first_empty)
        this%first_empty%mem_id=1_i4
        NULLIFY(this%first_empty%next)
        NULLIFY(this%first_empty%last)
        update_logmem=.TRUE. ! exit omp critical section and register.
      ELSE
        update_logmem=.FALSE.
      ENDIF
      !$omp end critical
    ELSE
      update_logmem=.FALSE.
    ENDIF

    IF (update_logmem) THEN
      CALL memlogger%update(this%my_mem_id,'mem_table','singleton',     &
                            INT(SIZEOF(this%obj_size)                   &
                            +SIZEOF(this%obj_list)                      &
                            +SIZEOF(this%obj_desc),i4))
    ENDIF
    IF (this%current_fill>=num_obj_max) THEN
      CALL nim_stop('memlogger: memtable full, increase num_obj_max')
    ENDIF

    !-------------------------------------------------------------------
    ! assign this object a mem_id, increment the empty list, and 
    ! store object values.
    !-------------------------------------------------------------------
    !$omp critical
    mem_id=>this%obj_list(2,this%obj_list(1,this%first_empty%mem_id))
    this%current_fill=this%current_fill+1_i4
    IF (.NOT. ASSOCIATED(this%first_empty%next)) THEN
      this%first_empty%mem_id=this%first_empty%mem_id+1_i4
    ELSE
      this%first_empty=>this%first_empty%next
      DEALLOCATE(this%first_empty%last)
      NULLIFY(this%first_empty%last)
    ENDIF
    !$omp end critical

    slen=MIN(strlen,LEN(obj_type))
    this%obj_desc(1,mem_id)=TRIM(obj_type(1:slen))
    slen=MIN(strlen,LEN(obj_name))
    this%obj_desc(2,mem_id)=TRIM(obj_name(1:slen))
    this%obj_size(mem_id)=mem_size
#ifdef DEBUG
    IF (mem_size<0) THEN
      WRITE(nim_wr,'(a)') 'mem_logger warning: registering a negative'  & 
                 //' memory object!'
    ENDIF
#endif

    RETURN
  END SUBROUTINE add_entry

!-----------------------------------------------------------------------
!> @brief modify an entry of the linked list.
!! @param[in] mem_id 
!! @param[in] mem_size
!-----------------------------------------------------------------------
  REAL(r8) FUNCTION modify_entry(this,mem_id,mem_size,resize)        &
              RESULT(newmem)
    CLASS(mem_table), INTENT(INOUT) :: this
    INTEGER(i4), POINTER, INTENT(IN) :: mem_id
    REAL(r8), INTENT(IN) :: mem_size
    LOGICAL :: resize

    !$omp critical
    IF (resize) THEN
      this%obj_size(mem_id)=mem_size
    ELSE
      this%obj_size(mem_id)=this%obj_size(mem_id)+mem_size
    ENDIF
    newmem=this%obj_size(mem_id)
    !$omp end critical

    RETURN
  END FUNCTION modify_entry

!-----------------------------------------------------------------------
!> @brief delete an entry and update the linked list,
!! @param[inout] mem_id 
!-----------------------------------------------------------------------
  SUBROUTINE delete_entry(this,mem_id)
    CLASS(mem_table), INTENT(INOUT) :: this
    INTEGER(i4), POINTER, INTENT(INOUT) :: mem_id
    LOGICAL :: first
 
    !$omp critical
    this%obj_size(mem_id)=0
    this%obj_desc(:,mem_id)=''
    this%current_fill=this%current_fill-1_i4

    ALLOCATE(this%first_empty%last)
    this%first_empty%last%next=>this%first_empty
    this%first_empty=>this%first_empty%last
    this%first_empty%mem_id=mem_id
    NULLIFY(this%first_empty%last)
    !$omp end critical

    RETURN
  END SUBROUTINE delete_entry

!-----------------------------------------------------------------------
!> @brief remove unused entries and sort the mem_list
!-----------------------------------------------------------------------
  SUBROUTINE sort_list(this)
    CLASS(mem_table), INTENT(INOUT) :: this
    
    INTEGER(i4) :: shift,ii,tmp,icnt

    !-------------------------------------------------------------------
    ! fill in arrays.
    !-------------------------------------------------------------------
    shift=0_i4
    DO ii=1,this%current_size
      IF (this%obj_size(ii)==0_i4) THEN
        shift=shift+1_i4
        IF (ii-shift==this%current_fill) EXIT
        this%obj_list(2,ii)=this%current_fill+shift
        CYCLE
      ENDIF
      IF (shift/=0_i4) THEN
        this%obj_size(ii-shift)=this%obj_size(ii)
        this%obj_size(ii)=0
        tmp=this%obj_list(1,ii-shift)
        this%obj_list(1,ii-shift)=this%obj_list(1,ii)
        this%obj_list(1,ii)=tmp
        this%obj_list(2,ii)=ii-shift
        this%obj_desc(:,ii-shift)=this%obj_desc(:,ii)
        this%obj_desc(:,ii)=''
      ENDIF
    ENDDO

    !-------------------------------------------------------------------
    ! reset the mem_id linked list.
    !-------------------------------------------------------------------
    DO WHILE (ASSOCIATED(this%first_empty%next))
      this%first_empty=>this%first_empty%next
      DEALLOCATE(this%first_empty%last)
    ENDDO
    NULLIFY(this%first_empty%last)
    this%first_empty%mem_id=this%current_fill+1

    !-------------------------------------------------------------------
    ! sort the packed list using a quick sort.
    !-------------------------------------------------------------------
    CALL qsort(this%obj_size(:this%current_fill),                       &
               this%obj_list(1,:this%current_fill),this%obj_list(2,:),  &
               this%obj_desc(:,:this%current_fill),1_i4)

    RETURN ! do only first index for now.     
    !-------------------------------------------------------------------
    ! now sort based on the second index.
    !-------------------------------------------------------------------
    icnt=1
    shift=1
    DO 
      IF (this%obj_desc(1,icnt)==this%obj_desc(1,shift)) THEN
        icnt=icnt+1
      ELSE
        CALL qsort(this%obj_size(shift:icnt),                           &
                   this%obj_list(1,shift:icnt),this%obj_list(2,:),      &
                   this%obj_desc(:,shift:icnt),2_i4)
        shift=icnt
      ENDIF
      IF (this%obj_size(icnt)==0) EXIT
    ENDDO
     
    RETURN
    CONTAINS

    !-------------------------------------------------------------------
    ! recursive quick-sort algorithm.
    !-------------------------------------------------------------------
    PURE RECURSIVE SUBROUTINE qsort(arrsz,arrl1,arrl2,arrstr,ip)
    REAL(r8), DIMENSION(:), INTENT(INOUT) :: arrsz
    INTEGER(i4), DIMENSION(:), INTENT(INOUT) :: arrl1
    INTEGER(i4), DIMENSION(:), INTENT(INOUT) :: arrl2
    CHARACTER(len=strlen), DIMENSION(:,:), INTENT(INOUT) :: arrstr
    INTEGER(i4), INTENT(IN) :: ip
    INTEGER(i4) :: iq

    IF (SIZE(arrsz)>1) THEN
      CALL partition(arrsz,arrl1,arrl2,arrstr,iq,ip)
      CALL qsort(arrsz(:iq-1),arrl1(:iq-1),arrl2,arrstr(:,:iq-1),ip)
      CALL qsort(arrsz(iq:),  arrl1(iq:),  arrl2,arrstr(:,iq:),  ip)
    ENDIF
    RETURN
    END SUBROUTINE qsort

    !-------------------------------------------------------------------
    ! partition into two sorted sub-arrays greater and less than the 
    ! pivot.
    !-------------------------------------------------------------------
    PURE SUBROUTINE partition(arrsz,arrl1,arrl2,arrstr,iq,ip)
    REAL(r8), DIMENSION(:), INTENT(INOUT) :: arrsz
    INTEGER(i4), DIMENSION(:), INTENT(INOUT) :: arrl1
    INTEGER(i4), DIMENSION(:), INTENT(INOUT) :: arrl2
    CHARACTER(len=strlen), DIMENSION(:,:), INTENT(INOUT) :: arrstr
    INTEGER(i4), INTENT(OUT) :: iq    ! split index
    INTEGER(i4), INTENT(IN) :: ip     ! pivot index of arr2
    
    INTEGER(i4) :: ii,jj,tmp,ipivot
    REAL(r8) :: rtmp
    INTEGER(i4), DIMENSION(2) :: temp1
    CHARACTER(strlen), DIMENSION(2) :: stmp
    CHARACTER(strlen) :: pivot

    ipivot=SIZE(arrstr,DIM=2)/2
    pivot=TRIM(arrstr(ip,ipivot))
    ii=1
    jj=SIZE(arrstr,DIM=2)
    DO WHILE (ii <= jj)
      DO WHILE (arrstr(ip,jj) > pivot)
        jj=jj-1
      ENDDO
      DO WHILE (arrstr(ip,ii) < pivot) 
        ii=ii+1
      ENDDO
      IF (ii < jj) THEN ! exchange ii and jj
        rtmp=arrsz(ii);    arrsz(ii)=arrsz(jj);       arrsz(jj)=rtmp
        tmp=arrl1(ii);     arrl1(ii)=arrl1(jj);       arrl1(jj)=tmp
        stmp=arrstr(:,ii); arrstr(:,ii)=arrstr(:,jj); arrstr(:,jj)=stmp
        tmp=arrl2(arrl1(ii)) 
        arrl2(arrl1(ii))=arrl2(arrl1(jj))
        arrl2(arrl1(jj))=tmp
        jj=jj-1
        ii=ii+1
      ELSEIF (ii==jj) THEN
        jj=jj-1
        ii=ii+1
      ENDIF
    ENDDO
    iq=ii
    END SUBROUTINE partition

  END SUBROUTINE sort_list
 
!-----------------------------------------------------------------------
!> @brief destructor to free the empty_list.
!-----------------------------------------------------------------------
  SUBROUTINE free_list(this)
    CLASS(mem_table), INTENT(INOUT) :: this 
 
    DO WHILE(ASSOCIATED(this%first_empty%next))
      this%first_empty=>this%first_empty%next
      DEALLOCATE(this%first_empty%last)
    ENDDO
    DEALLOCATE(this%first_empty)
    NULLIFY(this%first_empty)
    IF (ALLOCATED(this%obj_size)) DEALLOCATE(this%obj_size)
    IF (ALLOCATED(this%obj_list)) DEALLOCATE(this%obj_list)
    IF (ALLOCATED(this%obj_desc)) DEALLOCATE(this%obj_desc)
    this%current_fill=0

  END SUBROUTINE free_list

!-----------------------------------------------------------------------
!> @brief the update routine registers and tracks memory allocation.
!! @details the update routine registers an object by name and type, 
!!  increments or decrements the memory in that object. the mem_id
!!  should initially be set to zero, which indicates to the logger
!!  that the object is not registered. after registering the object
!!  it will return a unique identifier to be stored in the object and
!!  used during subsequent calls. mem_size is the size of the memory
!!  being allocated (positive sign) or deallocated (negative sign) in
!!  bytes. to avoid overflows for large objects, set inMB=true and 
!!  pass mem_size in megabytes.
!! @param[inout] mem_id 
!! @param[in] obj_type
!! @param[in] obj_name
!! @param[in] mem_size
!! @param[in] resize
!! @param[in] inMB
!-----------------------------------------------------------------------
  SUBROUTINE update(this,mem_id,obj_type,obj_name,mem_size,resize,inMb)
    CLASS(mem_logger), INTENT(INOUT) :: this
    INTEGER(i4), POINTER, INTENT(INOUT) :: mem_id
    CHARACTER(len=*), INTENT(IN) :: obj_type
    CHARACTER(len=*), INTENT(IN) :: obj_name
    INTEGER(i4), INTENT(IN) :: mem_size
    LOGICAL, OPTIONAL :: resize
    LOGICAL, OPTIONAL :: inMB

    REAL(r8) :: total_mem,size_MB
    LOGICAL :: rsz

    rsz=.FALSE.
    IF (PRESENT(resize)) rsz=resize
 
    IF (PRESENT(inMB)) THEN
      IF (inMB) THEN
        size_MB=mem_size
      ELSE
        size_MB=mem_size/1e6_r8
      ENDIF
    ELSE
      size_MB=mem_size/1e6_r8
    ENDIF

    IF (.NOT. ASSOCIATED(mem_id)) THEN
      CALL this%memtable%add_entry(mem_id,obj_type,obj_name,size_MB)
    ELSE
      total_mem=this%memtable%modify_entry(mem_id,size_MB,rsz)
      IF (total_mem < 1e-6) THEN
#ifdef DEBUG
        IF (total_mem < 0) THEN
          WRITE(*,*) 'mem_logger warning: object memory less than zero'
        ENDIF
#endif
        CALL this%memtable%delete_entry(mem_id)
        NULLIFY(mem_id)
      ENDIF
    ENDIF

    RETURN
  END SUBROUTINE update 

!-----------------------------------------------------------------------
!> @brief report memory usage statistics in MegaBytes.
!> @details verbosity level indicates can be set to
!!   (0) aggregate statistics over all MPI processes, and
!!   (1) full statistics on each registered object by MPI process 
!!       are output to file at each report call (step and locstr
!!       must be passed).
!! @param[in] verbosity
!! @param[in] step
!! @param[in] locstr
!-----------------------------------------------------------------------
  SUBROUTINE report(this,verbosity,step,locstr)
    USE iso_c_binding
    USE io, ONLY: nim_wr,memlog_unit
#ifdef HAVE_OPENMP
      USE omp_lib
#endif 
#ifdef HAVE_MPI
    INCLUDE "mpif.h"
#endif
    CLASS(mem_logger), INTENT(INOUT) :: this
    INTEGER(i4), OPTIONAL :: verbosity
    INTEGER(i4), OPTIONAL :: step
    CHARACTER(LEN=*), OPTIONAL :: locstr
   
    INTEGER(i4) :: ii,jj,shift,icnt,verbose,curlen,node,nprocs,ierr,it
    INTEGER(i4) :: iexit,allexit
    REAL(r8) :: pagesize,VmRSS,VmSize
    REAL(r8), DIMENSION(2) :: arr ! num obj, memory (Mb)
    REAL(r8), ALLOCATABLE :: aggarr(:,:), sndarr(:)
    CHARACTER(80) :: strtmp,repstr,dirstr
    CHARACTER(strlen) :: obj_type
    CHARACTER(strlen), ALLOCATABLE :: chararr(:)
    
    INTERFACE
      FUNCTION mkdir(path,mode) BIND(c,NAME="mkdir")
        USE iso_c_binding
        INTEGER(c_int) :: mkdir
        CHARACTER(KIND=c_char,LEN=1) :: path(*)
        INTEGER(c_int16_t), VALUE :: mode
      END FUNCTION mkdir
    END INTERFACE

!-----------------------------------------------------------------------
!   we're going to query /proc/self/statm which contains the memory
!   used by pid on unix based systems. this is not very portable,
!   so undefine USE_STATM if it is not compiling on your system.
!-----------------------------------------------------------------------
#define USE_STATM
#ifdef USE_STATM
    INTERFACE
      FUNCTION getpagesize() BIND(c,NAME="getpagesize")
        USE iso_c_binding
        INTEGER(c_int) :: getpagesize
      END FUNCTION getpagesize
    END INTERFACE
#endif

#ifdef HAVE_OPENMP
      it=OMP_GET_THREAD_NUM()
#else
      it=0
#endif
    !-------------------------------------------------------------------
    ! initializations.
    !-------------------------------------------------------------------
    IF (.NOT. ALLOCATED(this%memtable%obj_size)) RETURN
    IF (it /= 0) RETURN ! master thread only
    verbose=0_i4
    IF (PRESENT(verbosity)) verbose=verbosity
#ifdef HAVE_MPI
    node=this%node
    nprocs=this%nprocs
    IF (node<0) THEN 
      CALL mpi_comm_rank(mpi_comm_world,node,ierr)
      CALL mpi_comm_size(mpi_comm_world,nprocs,ierr)
      this%node=node
      this%nprocs=nprocs
    ENDIF
    IF (node==0) THEN
      ALLOCATE(aggarr(nprocs,2),sndarr(2*nprocs))
      ALLOCATE(chararr(nprocs))
    ENDIF
#else
    node=0
#endif
#ifdef USE_STATM
    pagesize=REAL(getpagesize(),r8)
    OPEN(UNIT=memlog_unit,STATUS='OLD',FILE='/proc/self/statm',         &
         FORM='FORMATTED')
    READ(memlog_unit,*) VmSize, VmRss
    CLOSE(UNIT=memlog_unit)
    !-------------------------------------------------------------------
    ! Resident virtual memory (RAM space used, or virtual memory 
    ! which has caused a page fault). Conversion to real avoids
    ! overflows.
    !-------------------------------------------------------------------
    VmRSS=VmRSS*pagesize/1e6_r8
    IF (VmRSS<=0.OR.VmRSS>1e15) VmRSS=1 ! Ensure sanity
    CALL this%update(this%vrss_mem_id,'VmRSS','VmRSS',                  &
                     INT(VmRSS,i4),resize=.TRUE.,inMB=.TRUE.)
    !-------------------------------------------------------------------
    ! Total virtual memory used.
    !-------------------------------------------------------------------
    VmSize=VmSize*pagesize/1e6_r8
    IF (VmSize<=0.OR.VmSize>1e15) VmSize=1 ! Ensure sanity
    CALL this%update(this%vsiz_mem_id,'VmSize','VmSize',                &
                     INT(VmSize,i4),resize=.TRUE.,inMB=.TRUE.)    
#endif
    !-------------------------------------------------------------------
    ! create an output directory for verbose output by node.
    ! if verbose info, each node opens anad writes a file with its
    ! state information. these files are placed in the directory
    ! obj_mem_#/step_locstring/
    !-------------------------------------------------------------------
    IF (verbose>=1) THEN 
      IF (.NOT. PRESENT(step) .OR. .NOT. PRESENT(locstr))               &
        CALL nim_stop('mem_logger::report : step and locstr are '//     &
                      'required for verbose output.')
      IF (rootdir(1:4)=='none') THEN
        IF (node==0) THEN
          icnt=0
          DO 
            WRITE(strtmp,'(i4)') icnt
            rootdir='mem_log_'//TRIM(ADJUSTL(strtmp))
            ii=mkdir(TRIM(rootdir)//c_null_char,INT(o'755',c_int16_t))
            IF (ii/=0_i4.AND.icnt<10000) THEN
              icnt=icnt+1
            ELSE 
              EXIT
            ENDIF
          ENDDO 
        ENDIF
#ifdef HAVE_MPI
        CALL mpi_bcast(rootdir,16,mpi_character,0,mpi_comm_world,ierr)
#endif
      ENDIF
      WRITE(strtmp,'(i16)') step
      dirstr=TRIM(rootdir)//'/'//TRIM(ADJUSTL(strtmp))//'_'//          &
             TRIM(locstr)//'/'
      IF (node==0) THEN
        ii=mkdir(TRIM(dirstr)//c_null_char,INT(o'755',c_int16_t))
        IF (ii/=0_i4) CALL nim_stop('mem_logger I/O error')
      ENDIF
#ifdef HAVE_MPI
      CALL mpi_barrier(mpi_comm_world,ierr)
#endif
      WRITE(strtmp,'(i16)') node
      OPEN(UNIT=memlog_unit,STATUS='UNKNOWN',POSITION='APPEND',        &
          FILE=TRIM(dirstr)//'node'//TRIM(ADJUSTL(strtmp))//'.out')
    ENDIF
 
    !-------------------------------------------------------------------
    ! sort the arrays for aggregation.
    !-------------------------------------------------------------------
    CALL this%memtable%sort_list()

    !-------------------------------------------------------------------
    ! report memory usage by object type and instance. after the sort
    ! the arrays are densely packed.
    !-------------------------------------------------------------------
    shift=1
    ii=1
    arr=0
    IF (verbose>=1_i4) WRITE(memlog_unit,*) '!-- '//                    &
                       TRIM(this%memtable%obj_desc(1,shift))//' --!'
#ifdef HAVE_MPI
    IF (node==0) THEN
      repstr='Memory statistics from MPI processes'
      IF (PRESENT(locstr)) repstr=TRIM(repstr)//' at '//TRIM(locstr)
      WRITE(nim_wr,'(a)') ''
      WRITE(nim_wr,'(a)') repstr
      repstr=''
      repstr(strlen+2:)=                                                &
                  '# objs                 memory (Mb)                  '
      WRITE(nim_wr,'(a)') repstr
      repstr=''
      repstr(1:4)='type'
      repstr(strlen+2:)=                                                &
                  'total  / high / low  | total  / high   / low '
      WRITE(nim_wr,'(a)') repstr
    ENDIF
    iexit=1
    !-------------------------------------------------------------------
    ! use an mpi_gather / broadcast to choose the first object
    !-------------------------------------------------------------------
    obj_type=this%memtable%obj_desc(1,1)
    CALL mpi_gather(obj_type,strlen,mpi_character,chararr,strlen,       &
                    mpi_character,0_i4,mpi_comm_world,ierr)
    IF (node==0) THEN
      DO jj=2,nprocs       
        IF (obj_type > chararr(jj) .AND. chararr(jj)/=' ')              &
              obj_type=chararr(jj)
      ENDDO
    ENDIF
    CALL mpi_bcast(obj_type,strlen,mpi_character,0,mpi_comm_world,ierr)
#else
    repstr='Memory statistics'
    IF (PRESENT(locstr)) repstr=TRIM(repstr)//' at '//TRIM(locstr)
    WRITE(nim_wr,'(a)') ''
    WRITE(nim_wr,'(a)') repstr
    repstr=''
    repstr(1:4)='type'
    repstr(strlen+2:)='# objs   memory Mb'
    WRITE(nim_wr,'(a)') repstr
    obj_type=this%memtable%obj_desc(1,1)
#endif
    obj_loop: DO
      !-----------------------------------------------------------------
      ! increment arrays if the same object type.
      !-----------------------------------------------------------------
      IF ( ii<=this%memtable%current_fill                               &
                    .AND. this%memtable%obj_desc(1,ii)==obj_type ) THEN
        arr(1)=arr(1)+1
        arr(2)=arr(2)+this%memtable%obj_size(ii)
        IF (verbose>=1_i4) THEN
          WRITE(strtmp,'(es10.2)') this%memtable%obj_size(ii)
          WRITE(memlog_unit,*) TRIM(this%memtable%obj_desc(2,ii))//     &
                               ' '//TRIM(ADJUSTL(strtmp))//' Mb'
        ENDIF
        ii=ii+1
      !-----------------------------------------------------------------
      ! the object type differs. send data the restart the counters.
      !-----------------------------------------------------------------
      ELSE
        !---------------------------------------------------------------
        ! write out statistics.
        !---------------------------------------------------------------
#ifdef HAVE_MPI
        !---------------------------------------------------------------
        ! send data to node 0, for write operation.
        ! Then handle the case of a process with 
        ! no objects of a given type. Could be node=0 or node>0.
        !---------------------------------------------------------------
        CALL mpi_gather(arr,2_i4,mpi_double_precision,sndarr,2_i4,      &
                        mpi_double_precision,0_i4,mpi_comm_world,ierr)
        IF (node==0) THEN
          aggarr=RESHAPE(sndarr,(/nprocs,2/),ORDER=(/2,1/))
          repstr=obj_type//' '
          WRITE(strtmp,'(i6)') INT(SUM(aggarr(:,1)),i4)
          repstr=repstr(1:strlen)//' '//ADJUSTL(strtmp)
          WRITE(strtmp,'(i6)') INT(MAXVAL(aggarr(:,1)),i4)
          repstr=repstr(1:strlen+8)//' '//ADJUSTL(strtmp)
          WRITE(strtmp,'(i6)') INT(MINVAL(aggarr(:,1)),i4)
          repstr=repstr(1:strlen+15)//' '//ADJUSTL(strtmp)
          WRITE(strtmp,'(es8.2)') SUM(aggarr(:,2))
          repstr=repstr(1:strlen+22)//' '//ADJUSTL(strtmp)
          WRITE(strtmp,'(es8.2)') MAXVAL(aggarr(:,2))
          repstr=repstr(1:strlen+31)//' '//ADJUSTL(strtmp)
          WRITE(strtmp,'(es8.2)') MINVAL(aggarr(:,2))
          repstr=repstr(1:strlen+40)//' '//ADJUSTL(strtmp)
          WRITE(nim_wr,'(a)') repstr
        ENDIF
        !---------------------------------------------------------------
        ! use an mpi_allreduce to determine exit condition.
        !---------------------------------------------------------------
        IF (ii-1 >= this%memtable%current_fill) iexit=0
        CALL mpi_allreduce(iexit,allexit,1,mpi_integer,mpi_max,         &
                           mpi_comm_world,ierr)
        IF (allexit==0) EXIT
        !---------------------------------------------------------------
        ! use an mpi_gather / broadcast to choose next object.
        !---------------------------------------------------------------
        obj_type=this%memtable%obj_desc(1,ii)
        CALL mpi_gather(obj_type,strlen,mpi_character,chararr,strlen,   &
                        mpi_character,0_i4,mpi_comm_world,ierr)
        IF (node==0) THEN
          DO jj=2,nprocs       
            IF (obj_type > chararr(jj) .AND. chararr(jj)>' ')          &
              obj_type=chararr(jj)
          ENDDO
        ENDIF
        CALL mpi_bcast(obj_type,strlen,mpi_character,0,                 &
                       mpi_comm_world,ierr)
#else 
        !---------------------------------------------------------------
        ! serial write.
        !---------------------------------------------------------------
        repstr=this%memtable%obj_desc(1,shift)//' '
        strtmp=''; WRITE(strtmp,'(i6)') INT(arr(1))
        repstr(1:strlen+9)=repstr(1:strlen)//' '//ADJUSTL(strtmp)
        strtmp=''; WRITE(strtmp,'(es8.2)') REAL(arr(2))/REAL(1.e6,r8)
        repstr=repstr(1:strlen+9)//' '//ADJUSTL(strtmp)
        WRITE(nim_wr,'(a)') repstr
        IF (ii-1>=this%memtable%current_fill) EXIT
        obj_type=this%memtable%obj_desc(1,ii)
#endif
        arr=0
        shift=ii
        IF (verbose>=1_i4) WRITE(memlog_unit,*) '!-- '//                &
                           TRIM(this%memtable%obj_desc(1,shift))//' --!'
      ENDIF

    ENDDO obj_loop

    IF (verbose>=1) CLOSE(UNIT=memlog_unit)
    
    RETURN
  END SUBROUTINE report

!-----------------------------------------------------------------------
!> @brief call destructors. can be removed with FINAL implementation
!!  in gfortran.
!-----------------------------------------------------------------------
  SUBROUTINE finalize(this)
    CLASS(mem_logger), INTENT(INOUT) :: this 
 
    CALL this%memtable%free_list

  END SUBROUTINE finalize

END MODULE memlog

!-----------------------------------------------------------------------
!> @brief this module is a wrapper of the update function as a 
!!  subroutine to be used in calls from C for SuperLU. 
!! @details the structure of the module avoids passing pointers 
!!  between C and Fortran at the cost of a small pointer array.
!-----------------------------------------------------------------------
MODULE memlog_c
  USE iso_c_binding
  USE local, ONLY: i4
  IMPLICIT NONE
  PUBLIC memlog_update
  PRIVATE

  INTEGER(i4), PARAMETER :: strlen=16 ! Should be the same as memlog

  TYPE cptr_arr
    INTEGER(i4), POINTER :: cptr
  END TYPE cptr_arr

  TYPE(cptr_arr), ALLOCATABLE :: cptrarr(:)
  INTEGER(i4) :: cptrfill

CONTAINS
!-----------------------------------------------------------------------
!> @brief function is a wrapper of the memlog_update function for C.
!! @details the proper declaration in C is
!!  > extern void* memlog_update(int *mem_id, char *obj_type,
!!  >                char *obj_name, int *mem_size, int *c_resize);
!-----------------------------------------------------------------------
  SUBROUTINE memlog_update(mem_id,obj_type,obj_name,mem_size,c_resize)  &
    BIND(C,name='memlog_update') 
    USE memlog, ONLY: memlogger
    INTEGER(c_int), INTENT(INOUT) :: mem_id
    CHARACTER(KIND=c_char,LEN=1), INTENT(IN) :: obj_type(*)
    CHARACTER(KIND=c_char,LEN=1), INTENT(IN) :: obj_name(*)
    INTEGER(c_int), INTENT(IN) :: mem_size
    INTEGER(c_int), INTENT(IN) :: c_resize

    INTEGER(i4), POINTER :: fmem_id  
    INTEGER(i4) :: ii
    CHARACTER(LEN=strlen) :: ftype,fname
    LOGICAL :: fresize
    TYPE(cptr_arr), ALLOCATABLE :: tmparr(:)

    TYPE(c_ptr) :: nullptr

    !-------------------------------------------------------------------
    ! handle allocation and resizing
    !-------------------------------------------------------------------
    IF (.NOT. ALLOCATED(cptrarr)) THEN 
      ALLOCATE(cptrarr(10))
      DO ii=1,SIZE(cptrarr)
        cptrarr(ii)%cptr=>NULL()
      ENDDO
      cptrfill=0
    ENDIF
    IF (cptrfill==SIZE(cptrarr)) THEN
      ALLOCATE(tmparr(SIZE(cptrarr)))
      tmparr=cptrarr
      DEALLOCATE(cptrarr)
      ALLOCATE(cptrarr(2*SIZE(cptrarr)))
      cptrarr(1:SIZE(tmparr))=tmparr
      DO ii=SIZE(tmparr)+1,SIZE(cptrarr)
        cptrarr(ii)%cptr=>NULL()
      ENDDO
      DEALLOCATE(tmparr)
    ENDIF

    IF (mem_id==0) THEN
      cptrfill=cptrfill+1
      mem_id=cptrfill
    ENDIF 

    !-------------------------------------------------------------------
    ! convert c to fortran 
    !-------------------------------------------------------------------
    ii=1; ftype=''
    DO WHILE(obj_type(ii)/=c_null_char .AND. ii<strlen)
      ftype=TRIM(ftype)//obj_type(ii)
      ii=ii+1
    ENDDO
    ii=1; fname=''
    DO WHILE(obj_name(ii)/=c_null_char .AND. ii<strlen)
      fname=TRIM(fname)//obj_name(ii)
      ii=ii+1
    ENDDO
    IF (c_resize==0) THEN
      fresize=.FALSE.
    ELSE
      fresize=.TRUE.
    ENDIF

    !-------------------------------------------------------------------
    ! call the update routine
    !-------------------------------------------------------------------
    CALL memlogger%update(cptrarr(mem_id)%cptr,ftype,fname,mem_size,    &
                          resize=fresize)

  END SUBROUTINE memlog_update

END MODULE memlog_c

