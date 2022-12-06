!----------------------------------------------------------------------
! $Id: mpi_serial.f90 4517 2015-07-29 23:08:48Z jking $
! this file is to be used on serial machines.  it contains two parts.
! the first is a module of parameters and the second is a set of
! dummy routines to replace mpi routines.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! module for interface to F77 MPI include file
! defines all MPI variables via parameter statements
! use this module to define machine-specific MPI datatypes
! so that NIMROD source will not have to change for a new
! machine or if system.f SELECTED_KINDS are changed
!-----------------------------------------------------------------------

      MODULE mpi_nim

      USE local
      implicit none

      INTEGER(i4), PARAMETER :: mpi_nim_int = 1_i4
      INTEGER(i4), PARAMETER :: mpi_nim_real = 1_i4
      INTEGER(i4), PARAMETER :: mpi_nim_double = 2_i4
      INTEGER(i4), PARAMETER :: mpi_nim_comp = 1_i4
      INTEGER(i4), PARAMETER :: mpi_nim_logical = 1_i4
      INTEGER(i4), PARAMETER :: mpi_nim_char = 1_i4

      INTEGER(i4), PARAMETER :: mpi_comm_world = 1_i4
      INTEGER(i4), PARAMETER :: mpi_max = 1_i4
      INTEGER(i4), PARAMETER :: mpi_min = 1_i4
      INTEGER(i4), PARAMETER :: mpi_sum = 1_i4
      INTEGER(i4), PARAMETER :: mpi_land = 1_i4
      INTEGER(i4), PARAMETER :: mpi_lor = 1_i4
      INTEGER(i4), PARAMETER :: mpi_source = 1_i4
      INTEGER(i4), PARAMETER :: mpi_any_source = 1_i4
      INTEGER(i4), PARAMETER :: mpi_any_tag = 1_i4
      INTEGER(i4), PARAMETER :: mpi_tag = 1_i4
      INTEGER(i4), PARAMETER :: mpi_status_size = 1_i4
      INTEGER(i4), PARAMETER :: mpi_packed = 1_i4
      INTEGER(i4), PARAMETER :: mpi_double_precision = 1_i4
      INTEGER(i4), PARAMETER :: mpi_double_complex = 1_i4
      INTEGER(i4), PARAMETER :: mpi_integer = 1_i4
      INTEGER(i4), PARAMETER :: mpi_request_null = 1_i4
      INTEGER(i4), PARAMETER :: mpi_thread_single = 0_i4
      INTEGER(i4), PARAMETER :: mpi_thread_funneled = 1_i4
      INTEGER(i4), PARAMETER :: mpi_thread_serialized = 2_i4
      INTEGER(i4), PARAMETER :: mpi_thread_multiple = 3_i4
      INTEGER(i4), PARAMETER :: mpi_info_null = 0_i4
      REAL(r8) :: mpi_wtime

! MPI communicators
! comm_nimrod = communicator for nimrod calculation.
!              Corresponds to mpi_comm_world in general case.
! comm_layer = within a layer, across blocks                            
!              all the procs owning blocks in the same layer
! comm_node0 = only node0, all other processes together
! comm_stride_layer = within a layer, with h5io_block_proc_stride

      integer(i4) :: comm_nimrod
      integer(i4) :: comm_layer
      integer(i4) :: comm_node0
      integer(i4) :: comm_stride_layer

      END MODULE mpi_nim

!------------------------------------------------------------------------
! stubs for MPI calls - use on serial machine to replace real MPI
! library so NIMROD will still compile
! except for mpi_comm_rank, mpi_comm_size and mpi_allreduce
! these are all no-operation rouines
!-----------------------------------------------------------------------
      subroutine mpi_testany(nrecv,recv_request,irecv,flag,status,      &
     &                       ierror)
      implicit none
      real :: recv_request
      real :: status
      integer nrecv,irecv,ierror
      dimension recv_request(nrecv),status(1)
      logical flag

      return
      end

!-----------------------------------------------------------------------
      subroutine comm_create_(nsend,procsend,comm,nrecv,plan)
      USE local
      implicit none
      integer nsend,nrecv,procsend,comm
      real(r8) plan

      nrecv=nsend
      plan=REAL(nsend)

      return
      end

!-----------------------------------------------------------------------
      subroutine comm_destroy_(plan)
      USE local
      implicit none
      real(r8) plan

      return
      end

!-----------------------------------------------------------------------
      subroutine comm_do_(plan,sendbuf,n,recvbuf)
      USE local
      implicit none
      integer n,sendbuf(*),recvbuf(*)
      real(r8) plan

      integer i

      do 10 i=1,nint(plan)
        recvbuf(i)=sendbuf(i)
   10 continue

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_init_nim(ierror)
      implicit none
      integer ierror,nprocs

      nprocs=1

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_finalize_nim(ierror)
      implicit none
      integer ierror

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_init_thread(tlvlin,tlvlout,ierror)
      implicit none
      integer tlvlin,tlvlout,ierror

      return
      end

!-----------------------------------------------------------------------
! return processor id = 0

      subroutine mpi_comm_rank_nim(comm,rank,ierror)
      implicit none
      integer comm,rank,ierror

      rank = 0

      return
      end

!-----------------------------------------------------------------------
! return # of processors = 1

      subroutine mpi_comm_nim(comm,size,ierror)
      implicit none
      integer comm,size,ierror

      size = 1

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_bcast_nim(buffer,count,datatype,root,comm,ierror)
      implicit none
      integer buffer,count,datatype,root,comm,ierror

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_barrier_nim(comm,ierror)
      implicit none
      integer comm,ierror

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_allreduce(sendbuf,recvbuf,count,                   &
     &     datatype,op,comm,ierror)
      implicit none
      integer sendbuf,recvbuf,count,datatype,op,comm,ierror
      dimension sendbuf(count),recvbuf(count)

      integer i

      do 10 i=1,count
        recvbuf(i)=sendbuf(i)
   10 continue

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_reduce_nim(sendbuf,recvbuf,count,                      &
     &     datatype,op,root,comm,ierror)
      implicit none
      integer sendbuf,recvbuf,count,datatype,op,root,comm,ierror
      dimension sendbuf(count),recvbuf(count)

      integer i

      do 10 i=1,count
        recvbuf(i)=sendbuf(i)
   10 continue

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_reduce_scatter(sendbuf,recvbuf,recvcounts,         &
     &     datatype,op,comm,ierror)
      implicit none
      integer sendbuf,recvbuf,recvcounts,datatype,op,comm,ierror

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_send(buf,count,datatype,dest,tag,comm,ierror)
      implicit none
      integer buf,count,datatype,dest,tag,comm,ierror

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_recv(buf,count,datatype,source,tag,comm,           &
     &     status,ierror)
      implicit none
      integer buf,count,datatype,source,tag,comm,status,ierror

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_irecv(buf,count,datatype,source,tag,comm,          &
     &     request,ierror)
      implicit none
      integer buf,count,datatype,source,tag,comm,request,ierror

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_isend(buf,count,datatype,source,tag,comm,          &
     &     request,ierror)
      implicit none
      integer buf,count,datatype,source,tag,comm,request,ierror

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_sendrecv(sendbuf,scounts,stypes,dest,sendtag,      &
     &        recvbuf,rcounts,rtypes,source,recvtag,comm,status,ierror)
      implicit none
      integer sendbuf,scounts,stypes,dest,sendtag,recvbuf,rcounts,      &
     &        rtypes,source,recvtag,comm,status,ierror
      dimension sendbuf(scounts),recvbuf(rcounts)
      integer i

      do 10 i=1,min(rcounts,scounts)
        recvbuf(i)=sendbuf(i)
   10 continue

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_test(request,flag,status,ierror)
      implicit none
      logical flag
      integer request,status,ierror

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_comm_split(datatype1,ilayer,ii,datatype2,ierror)
      implicit none
      integer datatype1,datatype2,ilayer,ierror,ii

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_wait(recv_request,status,ierror)
      implicit none
      real :: recv_reqest
      real :: recv_request
      real :: status
      real :: statuses
      integer nrecv,ierror
      dimension recv_reqest(1),statuses(1)

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_waitall(nrecv,recv_request,statuses,ierror)
      implicit none
      integer nrecv,ierror,recv_request,statuses
      dimension recv_request(nrecv),statuses(1,1)

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_waitany(nrecv,recv_request,irecv,status,ierror)
      implicit none
      integer nrecv,irecv,ierror,recv_request,status
      dimension recv_request(nrecv),status(1)

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_allgather(sendbuf,counts,datatypes,                &
     &        recvbuf,countr,displs,datatyper,comm,ierror)
      implicit none
      integer sendbuf,recvbuf,counts,countr,datatypes,datatyper,        &
     &        displs,comm,ierror
      dimension sendbuf(counts),recvbuf(counts)

      integer i

      do 10 i=1,counts
        recvbuf(i)=sendbuf(i)
   10 continue

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_gather(sendbuf,counts,datatypes,                   &
     &        recvbuf,countr,displs,datatyper,comm,ierror)
      implicit none
      integer sendbuf,recvbuf,counts,countr,datatypes,datatyper,        &
     &        displs,comm,ierror
      dimension sendbuf(counts),recvbuf(counts)

      integer i

      do 10 i=1,counts
        recvbuf(i)=sendbuf(i)
   10 continue

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_allgatherv(sendbuf,counts,datatypes,               &
     &        recvbuf,countr,displs,datatyper,comm,ierror)
      implicit none
      integer sendbuf,recvbuf,counts,countr,datatypes,datatyper,        &
     &        displs,comm,ierror
      dimension sendbuf(counts),recvbuf(counts)

      integer i

      do 10 i=1,counts
        recvbuf(i)=sendbuf(i)
   10 continue

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_alltoallv(sendbuf,scounts,sdispls,datatypes,       &
     &        recvbuf,rcounts,rdispls,datatyper,comm,ierror)
      implicit none
      integer sendbuf,scounts,sdispls,recvbuf,rcounts,rdispls,          &
     &        datatypes,datatyper,comm,ierror
      dimension sendbuf(scounts),recvbuf(rcounts)

      integer i

      do 10 i=1,min(rcounts,scounts)
        recvbuf(i)=sendbuf(i)
   10 continue

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_abort_nim(comm,code,ierror)
      implicit none
      integer comm,code,ierror

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_get_address(datum,addrs,ierror)
      integer datum,addrs,ierror

      addrs=0

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_address(datum,addrs,ierror)
      implicit none
      integer datum,addrs,ierror

      addrs=0

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_pack(inbuf,incount,dtype,outbuf,outcount,pos,      &
     &        comm,ierror)
      implicit none
      integer inbuf,incount,dtype,outbuf,outcount,pos,comm,ierror

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_unpack(inbuf,insize,pos,outbuf,outcount,dtype,     &
     &        comm,ierror)
      implicit none
      integer inbuf,insize,dtype,outbuf,outcount,pos,comm,ierror

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_pack_size(incount,dtype,comm,insize,ierror)
      implicit none
      integer incount,dtype,comm,insize,ierror

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_type_struct(incount,blength,displ,atype,           &
     &        newtype,ierror)
      implicit none
      integer incount,blength,displ,atype,newtype,ierror

      newtype=incount*blength

      return
      end

!-----------------------------------------------------------------------
      subroutine mpi_type_commit(dtype,ierror)
      implicit none
      integer dtype,ierror

      return
      end 
