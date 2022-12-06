!-----------------------------------------------------------------------
!     $Id: pardata.F90 7620 2022-01-17 21:09:29Z jking $
!     data structures for performing block decomposition and
!       seam communication in parallel
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     data structures.
!-----------------------------------------------------------------------
#include "config.f"
      MODULE pardata
      USE local
#ifdef HAVE_FFTW3
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      INCLUDE 'fftw3.f03'

      TYPE :: fftw_arr
        REAL(C_DOUBLE), POINTER :: rarr(:,:)
        COMPLEX(C_DOUBLE_COMPLEX), POINTER :: carr(:,:)
      END TYPE fftw_arr
#else
      IMPLICIT NONE
#endif

! Processor info for mpi_comm_world

      integer(i4) :: global_nprocs=1
      integer(i4) :: global_node=0

! Processor info for comm_nimrod
! nprocs = total # of processors (assigned to blocks and layers)
! node = processor id of me (0 to nprocs-1)

      integer(i4) :: nprocs=1
      integer(i4) :: node=0

! Processor info for comm_closure

      integer(i4) :: closure_nprocs=1
      integer(i4) :: closure_node=0

      integer(i4) :: ns_nprocs=1 ! for integrating over s
      integer(i4) :: ns_group=0  ! # of processor groups with separate s
      integer(i4) :: ns_layer=1  ! for integrating over s
      integer(i4) :: ns_lo(3)       ! lower speed index of ns_group
      integer(i4) :: ns_hi(3)       ! upper speed index of ns_group
      integer(i4) :: ns_num(3)      ! number of speed points per group
      integer(i4) :: xi_dof(3)      ! pitch-angle degrees of freedom

! Processor info for comm_edge

      integer(i4) :: edge_nprocs=1
      integer(i4) :: edge_node=0

! layer2proc(0:nlayers-1) = proc ids of all layers for the
!   local grid blocks.

      integer(i4), dimension(:), allocatable :: layer2proc

! nprocs_layer = # of procs assigned to each layer
! node_layer = I am this proc within the layer (0 to nprocs_layer-1)
! ilayer = which layer this proc belong to (0 to nlayers-1)
! mode_lo = my 1st mode corresponds to this global mode (1 to nmodes_total)
! mode_hi = my last mode corresponds to this global mode (1 to nmodes_total)
! mode2layer = maps each global mode (1 to nmodes_total) to a unique layer

      integer(i4) :: nprocs_layer
      integer(i4) :: node_layer
      integer(i4) :: ilayer
      integer(i4) :: mode_lo,mode_hi
      integer(i4), dimension(:), allocatable :: mode2layer

! MPI communicators
! comm_mode =  withih a block, across modes (layers)
!              all the procs owning modes for the same block(s)
! comm_xbl = for kprad closures
! comm_edge = within a layer, across edge blocks blocks
! comm_edge_io = all nodes with edge_node=0 for I/O
! comm_closure = group sharing same rblock but different speed

      integer(i4) :: comm_mode
      integer(i4) :: comm_xbl
      integer(i4) :: comm_edge
      integer(i4) :: comm_edge_io
      integer(i4) :: comm_closure

! block2proc(1:nbl_total)
!   block2proc(i) = proc id of owner of global block i (within my layer)
! global2local(1:nbl_total)
!   global2local(i) = local index (1:nbl) of global block i
!                     on proc who owns it (within my layer)
! loc2glob(1:nbl)
!   loc2glob(i) = global index (1:nbl_total) of local block i
! block_sizes(6,1:nbl_total)
!   block_sizes(1,i) = nvert (# of seam verts) in global block i (r or tbl)
!   block_sizes(2,i) = npt (# of seam pts) in global block i (r or tbl)
!   block_sizes(3,i) = mx (# of x vertices) in global rblock i (0 for tb)
!   block_sizes(4,i) = my (# of y vertices) in global rblock i (0 for tb)
!   block_sizes(5,i) = degenerate flag in global rblock i (0 for tb)
!   block_sizes(6,i) = number of edge vertices in global rblock i (0 for tb)

      integer(i4), dimension(:), allocatable :: block2proc
      integer(i4), dimension(:), allocatable :: global2local
      integer(i4), dimension(:), allocatable :: loc2glob
      integer(i4), dimension(:,:), allocatable :: block_sizes

! data stucture for message sending of SEAM DATA
! nsend = # of messages I will send
! allocate send(1:nsend) of send_type,
!   one for each proc a message will be sent to
! proc = processor id to send message to
! count = number of datums to send, a "datum" = nqty values
! block(count) = which of my local blocks (1:nbl) the datum comes from
! vertex(count) = which seam vertex (1:nvert) the datum comes from
! image(count) = which vertex image (1:nimage) the datum corresponds to,
! data(nqty*count) = all values packed into one message,
!                    is a 1-d vector so that data can be packed contiguously
!                    regardless of what nqty is on a particular seaming call
! cdata(nqty*count) = a complex version of data

      type :: send_type
        integer(i4) :: proc
        integer(i4) :: count
        integer(i4), dimension(:), pointer :: block
        integer(i4), dimension(:), pointer :: vertex
        integer(i4), dimension(:), pointer :: image
        real(r8), dimension(:), pointer :: data
        complex(r8), dimension(:), pointer :: cdata
      end type send_type

      integer(i4) :: nsend
      type(send_type), dimension(:), pointer :: send

! data stucture for message receiving of SEAM DATA
! nrecv = # of messages I will receive
! allocate recv(1:nrecv) of recv_type,
!   one for each proc a message will come from
! proc = processor id of who will send me the messaeg
! count = number of datums I will receive, a "datum" = nqty values
! block(count) = which of my local blocks (1:nbl) the datum gets summed to
! vertex(count) = which seam vertex in my block the datum corresponds to
! order(count) = unique ordering of this image in all (except degenerate
!		 points) representations so sums produce identical round-off
! data(nqty*count) = all values packed into one message,
!                    received data will be packed contiguously
!                    regardless of what nqty is on a particular seaming call
! cdata(nqty*count) = a complex version of data
! recv_request = array of requests for posting asynchronous receives

      type :: recv_type
        integer(i4) :: proc
        integer(i4) :: count
        integer(i4), dimension(:), pointer :: block
        integer(i4), dimension(:), pointer :: vertex
        integer(i4), dimension(:), pointer :: order
        real(r8), dimension(:), pointer :: data
        complex(r8), dimension(:), pointer :: cdata
      end type recv_type

      integer(i4) :: nrecv
      type(recv_type), dimension(:), pointer :: recv
      integer(i4), dimension(:), allocatable :: recv_request

! data stucture for performing SEAMING of values between blocks I own
! nself = # of seam points whose images are also owned by me
! not needed when nprocs = 1, since can just loop over all seams
! allocate self(1:nself) of self type,
!   one for each point, note that a pair of
!   self-referencing pts will be stored twice
! block_out = which of my local blocks (1:nbl) the datum gets summed to
! vertex_out = which seam vertex (1:nvert) the datum gets summed to
! block_in = which of my local blocks (1:nbl) the datum comes from
! vertex_in = which seam vertex (1:nvert) the datum comes from
! order = unique ordering of this image in all (except degenerate
!	  points) representations so sums produce identical round-off

      type :: self_type
        integer(i4) block_out,vertex_out
        integer(i4) block_in,vertex_in
        integer(i4) order
      end type self_type

      integer(i4) :: nself
      type(self_type), dimension(:), pointer :: self

! data stucture for message sending of SEGMENT DATA
! nsendseg = # of messages I will send
! allocate sendseg(1:nsendseg) of sendseg_type,
!   one for each proc a message will be sent to
! proc = processor id to send message to
! count = number of datums to send, a "datum" = nqty x nqty values
! block(count) = which of my local blocks (1:nbl) the datum comes from
! segment(count) = which segment (1:nvert) the datum comes from
! data(nqty*count) = all values packed into one message,
!                    is a 1-d vector so that data can be packed contiguously
! cdata(nqty*count) = a complex version of data
! mat_data(nmat*nqty**2*count) = holds nmat off-diagonal matrix
!                                elements without symmetry assumptions
! cmat_data(nmat*nqty**2*count) = complex version of mat_data.

      type :: sendseg_type
        integer(i4) :: proc
        integer(i4) :: count
        integer(i4), dimension(:), pointer :: block
        integer(i4), dimension(:), pointer :: segment
        real(r8), dimension(:), pointer :: data
        complex(r8), dimension(:), pointer :: cdata
        real(r8), dimension(:), pointer :: mat_data
        complex(r8), dimension(:), pointer :: cmat_data
      end type sendseg_type

      integer(i4) :: nsendseg
      type(sendseg_type), dimension(:), pointer :: sendseg

! data stucture for message receiving of SEGMENT DATA
! nrecvseg = # of messages I will receive
! allocate recvseg(1:nrecvseg) of recvseg_type,
!   one for each proc a message will come from
! proc = processor id of who will send me the messaeg
! count = number of datums I will receive, a "datum" = nqty x nqty values
! block(count) = which of my local blocks (1:nbl) the datum gets summed to
! segment(count) = which segment in my block the datum corresponds to
! data(nqty*count) = all values packed into one message,
!                    received data will be packed contiguously
! cdata(nqty*count) = a complex version of data
! mat_data(nmat*nqty**2*count) = holds off-diagonal matrix elements
!                                without symmetry assumptions
! cmat_data(nmat*nqty**2*count) = complex version of mat_data.
! recvseg_request = array of requests for posting asynchronous receives

      type :: recvseg_type
        integer(i4) :: proc
        integer(i4) :: count
        integer(i4), dimension(:), pointer :: block
        integer(i4), dimension(:), pointer :: segment
        real(r8), dimension(:), pointer :: data
        complex(r8), dimension(:), pointer :: cdata
        real(r8), dimension(:), pointer :: mat_data
        complex(r8), dimension(:), pointer :: cmat_data
      end type recvseg_type

      integer(i4) :: nrecvseg
      type(recvseg_type), dimension(:), pointer :: recvseg
      integer(i4), dimension(:), allocatable :: recvseg_request

! data stucture for seaming SEGMENTS of values between blocks I own
! nselfseg = # of segments whose images are also owned by me
! not needed when nprocs = 1, since can just loop over all segments
! allocate selfseg(1:nselfseg) of selfseg type,
!   one for each segment, note that a pair of
!   self-referencing segments will be stored twice
! block_out = which of my local blocks (1:nbl) the datum gets summed to
! segment_out = which segment (1:nvert) the datum gets summed to
! block_in = which of my local blocks (1:nbl) the datum comes from
! segment_in = which segment (1:nvert) the datum comes from

      type :: selfseg_type
        integer(i4) block_out,segment_out
        integer(i4) block_in,segment_in
      end type selfseg_type

      integer(i4) :: nselfseg
      type(selfseg_type), dimension(:), pointer :: selfseg

! the line structure contains data required to initialize communication
! from blocks to global lines for preconditioning:

! data stored by line:
! nlinks = # of blocks contributing to this line.
! perpst & perpen = start and end indices in direction perpendicular to
!   the line.
! parast & paraen = start and end indices in the parallel direction for
!   computations with the line as a whole.  parast=1 for periodic lines
!   and paraen=0 for nonperiodic lines.
! ncomm = # of off-processor communications in the line.
! glb_bl_recv = global block number for the block with the same index
!   as this line.
! bl2line_parast & bl2line_paraen = parallel start and end indices
!   within a block that are sent to a line.
! line2bl_parast & line2bl_paraen = parallel start and end indices
!   within a block that are received from a line.
! par_dir = normal rblock index that represents the parallel direction

! data stored by link or segment of a line.
! node = node from which this line receives this segment.
! glb_bl = global block index holding this segment.
! loc_bl = local block index holding this segment.
! mpara & mperp = cell dimensions of the block containing this segment
!   in the directions parallel and perpendicular to this line.
! bl2line_segst & bl2line_segen = line parallel indices filled by this
!   segment.
! line2bl_segst & line2bl_segen = line parallel indices sent to this
!   segment.
! bl2line_count & line2bl_count = # of block vertices collected
!   from this segment and # shipped back, respectively.
! recv_req & send_req = request arrays used for nonblocking mpi calls
!   (both dimensioned ncomm for this line).
! req_index = converts from segment index to request array index.
! bl_perpst & bl_perpen = start and end perpendicular indices for
!   a segment of the block with this line's index.

      TYPE :: line_type
        INTEGER(i4) :: nlinks,perpst,perpen,parast,paraen,ncomm,        &
     &                 glb_bl_recv,bl2line_parast,bl2line_paraen,       &
     &                 line2bl_parast,line2bl_paraen,par_dir
        INTEGER(i4), DIMENSION(:), POINTER :: node,glb_bl,loc_bl,       &
     &               mpara,mperp,bl2line_segst,bl2line_segen,           &
     &               bl2line_count,line2bl_segst,line2bl_segen,         &
     &               line2bl_count,recv_req,send_req,req_index,         &
     &               bl_perpst,bl_perpen
        LOGICAL :: periodic
      END TYPE line_type

      TYPE(line_type), DIMENSION(:), ALLOCATABLE :: linex,liney

! threaded ffts used these shared arrays to permit MPI communication
! by the master thread.

      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: f_coefth,compth
      INTEGER(i4) :: fft_sync1,fft_sync2,fft_sync3,fft_sync4,fft_syncA
#ifdef HAVE_FFTW3
      TYPE(C_PTR), ALLOCATABLE :: plan_r2c(:),plan_c2r(:)
      TYPE(fftw_arr), ALLOCATABLE :: arrs(:)
      LOGICAL, ALLOCATABLE :: thread_alloc(:)
#endif

!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------

      END MODULE pardata
