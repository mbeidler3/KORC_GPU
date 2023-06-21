module pspline_gpu
  
  use params_gpu

  IMPLICIT NONE
    
#ifdef PSPLINE
  integer, parameter :: fp = rp
  real(fp), parameter :: ezspline_twopi = 6.2831853071795865_fp
  
  type :: EZspline2
      ! 2-d Spline/Akima Hermite/Piecewise Linear interpolation
      ! Grid
      real(fp), dimension(:), allocatable :: x1, x2
      ! The boundary condition values (for slope and 2nd derivative).
      ! Can be optionally set by the user. Not used for periodic and
      ! not a knot boundary conditions.
      real(fp), dimension(:), allocatable :: bcval1min, bcval1max
      real(fp), dimension(:), allocatable :: bcval2min, bcval2max
      ! Select between spline (0) and Akima spline (1); default=0 (spline)
      integer :: isHermite  ! set after EZspline_init call...
      ! set =0 for Spline, Akima or Hybrid; =1 for piecewise linear: this is set
      ! by EZspline_init, EZhybrid_init, or EZlinear_init; DO NOT SET DIRECTLY:
      integer :: isLinear
      ! set =0 by init routines other than EZhybrid_init which sets it =1:
      integer :: isHybrid
      ! the following is set by EZhybrid_init; other EZ*_init routines clear:
      integer :: hspline(2)  ! interpolation code along each dimension
      !        -1: zonal step fcn; =0: pc linear; =1: Akima Hermite; =2: Spline
      !
      ! Grid sizes (set during EZ*_init call).
      integer :: n1, n2
      ! Grid zone lookup method
      integer :: klookup1,klookup2
      ! Type of boundary conditions (set during EZspline_init call) on left
      ! and right hand side. Possible values are:
      !
      ! -1 periodic
      ! 0 not a knot
      ! +1 1st derivative imposed
      ! +2 2nd derivative imposed
      !
      ! For instance, ibctype1 =(/1, 0/) for 1st derivative set on left-hand
      ! and not a knot boundary conditions on right-hand side. The values of
      ! the derivatives are set via  bcval1min. (See above)
      integer ibctype1(2), ibctype2(2)
      ! Grid lengths. DO NOT SET.
      real(fp) :: x1min, x1max, x2min, x2max
      ! Compact cubic coefficient arrays. DO NOT SET.
      real(fp), dimension(:,:,:), allocatable :: fspl
      ! Control/Other. DO NOT SET
      integer :: isReady
      integer :: ilin1, ilin2
      real(fp), dimension(:,:), allocatable :: x1pkg, x2pkg
      integer :: nguard
  end type EZspline2
  
  PUBLIC :: EZspline_init2,EZspline_error,EZspline_setup2, &
    EZspline_interp2_FOvars_cloud,Ezspline_free2

#endif
  
  CONTAINS
  
#ifdef PSPLINE
  
  subroutine EZspline_init2(spline_o, n1, n2, BCS1, BCS2, ier)
      implicit none
      type(EZspline2) spline_o
      integer, intent(in) :: n1, n2
      integer, intent(in) :: BCS1(2), BCS2(2)
      ! ier:
      ! 0=ok
      ! 1=allocation error
      ! 2=wrong BCS1 code
      ! 3=wrong BCS2 code
      ! 99=something strange happened in EZspline_init
      integer, intent(out) :: ier
      integer i, iok
    
      call EZspline_preInit2(spline_o)
      ier = 0
    
      if(EZspline_allocated2(spline_o)) then
         ier = 100  ! allocated already
         return
      else
         call EZspline_preInit2(spline_o)
      end if
    
      spline_o%n1 = n1
      spline_o%n2 = n2
    
      spline_o%klookup1 = 3    ! default lookup algorithm selection
      spline_o%klookup2 = 3    ! default lookup algorithm selection
    
      spline_o%isHermite = 0 ! default is spline interpolation
      spline_o%isLinear  = 0 ! spline, not piecewise linear interpolation
      spline_o%isHybrid  = 0
      spline_o%hspline = 0
    
      iok = 0
      allocate(spline_o%x1(n1), stat=iok)
      if(iok /= 0) ier = 1
      allocate(spline_o%x2(n2), stat=iok)
      if(iok /= 0) ier = 1
      allocate(spline_o%fspl(4,n1,n2), stat=iok)
      if(iok /= 0) ier = 1
      allocate(spline_o%bcval1min(n2), stat=iok)
      if(iok /= 0) ier = 1
      allocate(spline_o%bcval1max(n2), stat=iok)
      if(iok /= 0) ier = 1
      allocate(spline_o%bcval2min(n1), stat=iok)
      if(iok /= 0) ier = 1
      allocate(spline_o%bcval2max(n1), stat=iok)
      if(iok /= 0) ier = 1
      allocate(spline_o%x1pkg(n1, 4), stat=iok)
      if(iok /= 0) ier = 1
      allocate(spline_o%x2pkg(n2, 4), stat=iok)
      if(iok /= 0) ier = 1
    
      do i = 1, 2
    
         spline_o%bcval1min(1:n2) = 0.0_fp
         spline_o%bcval1max(1:n2) = 0.0_fp
         select case(BCS1(i))
         case (-1)
            spline_o%ibctype1(i) = -1
         case (0)
            spline_o%ibctype1(i) = 0
         case (1)
            spline_o%ibctype1(i) = 1
         case (2)
            spline_o%ibctype1(i) = 2
         case default
            ier = 2
            spline_o%ibctype1(i) = 0
         end select
    
         spline_o%bcval2min(1:n1) = 0.0_fp
         spline_o%bcval2max(1:n1) = 0.0_fp
         select case(BCS2(i))
         case (-1)
            spline_o%ibctype2(i) = -1
         case (0)
            spline_o%ibctype2(i) = 0
         case (1)
            spline_o%ibctype2(i) = 1
         case (2)
            spline_o%ibctype2(i) = 2
         case default
            ier = 3
            spline_o%ibctype2(i) = 0
         end select
    
      end do
      if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) then
         spline_o%ibctype1(1)=-1
         spline_o%ibctype1(2)=-1
      end if
      if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) then
         spline_o%ibctype2(1)=-1
         spline_o%ibctype2(2)=-1
      end if
    
      !
      ! default grid
      spline_o%x1min = 0.0_fp
      spline_o%x1max = 1.0_fp
      if(BCS1(2)==-1) spline_o%x1max = ezspline_twopi
    
      spline_o%x1 = spline_o%x1min + (spline_o%x1max - spline_o%x1min)* &
           (/ (real(i-1,fp)/real(spline_o%n1-1, fp), i=1,spline_o%n1 ) /)
      spline_o%x2min = 0.0_fp
      spline_o%x2max = 1.0_fp
      if(BCS2(2)==-1) spline_o%x2max = ezspline_twopi
      spline_o%x2 = spline_o%x2min + (spline_o%x2max - spline_o%x2min)* &
           (/ (real(i-1,fp)/real(spline_o%n2-1, fp), i=1,spline_o%n2 ) /)
    
      spline_o%isReady = 0
    
      return
    end subroutine EZspline_init2
  
    subroutine EZspline_preInit2(spline_o)
     type(EZspline2) :: spline_o
     spline_o%nguard=123456789
     spline_o%isReady=0
     spline_o%ibctype1=0 ; spline_o%ibctype2=0
   end subroutine EZspline_preInit2
  
    subroutine EZspline_setup2(spline_o, f, ier, exact_dim)
      implicit none
      type(EZspline2) spline_o
      real(fp), dimension(:,:), intent(in) :: f
      ! ier:
      ! 0=ok
      ! 98=some error occurred in EZspline_setup
      integer, intent(out) :: ier
    
      logical, intent(in), OPTIONAL :: exact_dim  ! set T to require exact
      !  dimensioning match between f and spline_o%fspl; default is F
    
      logical :: iexact
      integer ifail
    
      integer :: ipx, ipy
      integer iper, imsg, itol, inum, ii, jj, in0, in1, in2
      real(fp) ztol, df1, df2
    
      !-------------------------
      iexact=.FALSE.
      if(present(exact_dim)) iexact = exact_dim
    
      if( .not.EZspline_allocated2(spline_o) ) then
         ier = 98
         return
      end if
    
      in0 = size(spline_o%fspl,1)
      in1 = size(spline_o%fspl,2)
      in2 = size(spline_o%fspl,3)
    
      ier = 57
      if(size(f,1).lt.in1) return
      if(size(f,2).lt.in2) return
    
      if(iexact) then
         ier = 58
         if(size(f,1).gt.in1) return
         if(size(f,2).gt.in2) return
      end if
      ier = 0
    
      !
      ! recompute min/max in case user changed the grid manually
      spline_o%x1max = maxval(spline_o%x1)
      spline_o%x2max = maxval(spline_o%x2)
      spline_o%x1min = minval(spline_o%x1)
      spline_o%x2min = minval(spline_o%x2)
    
      !
      ! set xpkg, useful for cloud and array interpolation
      imsg=0
      itol=0        ! range tolerance option
      ztol=5.e-7_fp ! range tolerance, if itol is set
      iper=0
      if(spline_o%ibctype1(1)==-1 .OR. spline_o%ibctype1(2)==-1) iper=1
      call genxpkg(spline_o%n1,spline_o%x1(1),spline_o%x1pkg(1,1),&
           iper,imsg,itol,ztol,spline_o%klookup1,ifail)
      if(ifail/=0) ier=27
      iper=0
      if(spline_o%ibctype2(1)==-1 .OR. spline_o%ibctype2(2)==-1) iper=1
      call genxpkg(spline_o%n2,spline_o%x2(1),spline_o%x2pkg(1,1),&
           iper,imsg,itol,ztol,spline_o%klookup2,ifail)
      if(ifail/=0) ier=27
    
      spline_o%isReady = 0
    
      spline_o%fspl(1,1:in1,1:in2) = f(1:in1,1:in2)
    
      ! this fixes a VMS f90 compiler optimizer problem:
      if(ztol.eq.-1.2345d30) &
           write(6,*) 'spline_o%fspl(1,1,1) = ', spline_o%fspl(1,1,1)
    
      call mkbicub(             &
      spline_o%x1(1), spline_o%n1, &
      spline_o%x2(1), spline_o%n2, &
      spline_o%fspl(1,1,1), spline_o%n1, &
      spline_o%ibctype1(1), spline_o%bcval1min(1), &
      spline_o%ibctype1(2), spline_o%bcval1max(1), &
      spline_o%ibctype2(1), spline_o%bcval2min(1), &
      spline_o%ibctype2(2), spline_o%bcval2max(1), &
      spline_o%ilin1, spline_o%ilin2, ifail)
 
      if(ifail /= 0) then
        ier = 98
      else
        spline_o%isReady = 1
      end if
    
      return
    end subroutine EZspline_setup2  
  
subroutine genxpkg(nx,x,xpkg,iper,imsg,itol,ztol,ialg,ier)
  !
  !  from an x axis assemble a "package":
  !
  !  MOD DMC Feb 2010: handle very small grids: nx=2, nx=3...
  !  interchanged meaning of xpkg(1,4) and xpkg(3,4);
  !  xpkg(3,4) only set if nx.ge.3;
  !  xpkg(4,4) only set if nx.ge.4.
  !
  !  there are corresponding changes in xlookup: simplified code lookup
  !  code for the cases nx=2 and nx=3.
  !
  !     xpkg(j,1) = x(j), j = 1 to nx    ... nx.ge.2
  !
  !     if(abs(ialg).ne.3) then...
  !       for j=1:nx-1
  !       xpkg(j,2) = h(j) = x(j+1)-x(j), j=1 to nx-1
  !     else
  !       for j=1:nx-1
  !       xpkg(j,2) = index location, with linear offset, of
  !                   xpkg(1,1)+<h>*(j-1) in the original x(1:nx)
  !            (piecewise linear indexing function)
  !     end if
  !     xpkg(nx,2) = <h> = (x(nx)-x(1))/(nx-1)
  !
  !     xpkg(j,3) = 1/h(j)
  !     xpkg(nx,3) = 1/<h>
  !
  !     xpkg(1,4) = +/- tolerance epsilon for out-of-range warnings/errors
  !                 +if message is to be written on out of range errors
  !                 -if no message to be written.  In either case, a
  !                        warning flag is set.
  !
  !     xpkg(2,4) = 1.0 if x is *periodic*, else 0.0
  !
  !  only set if nx.ge.3:
  !     xpkg(3,4) = 0.0 if x is *evenly spaced* (to ztol or 5.e-7 rel) else:
  !                   = 1.0 => use (1/h) newton's method like search algorithm
  !                   = 2.0 => use binary search algorithm
  !                   = 3.0 => use piecewise-linear indexing function...
  !
  !  only set if nx.ge.4:
  !     xpkg(4,4) = 0.0 -- do not use past search result as start point
  !                        for next search;
  !               = 1.0 -- do use past search result as start point for
  !                        next search.
  !
  !  tolerance epsilon means:
  !     if xtest is within epsilon*max(abs(x(1)),abs(x(nx))) outside the
  !     range [x(1),x(nx)] treat xtest as being at the endpoint x(1) or
  !     x(nx) whichever is closer.
  !
  ! input arguments:
  !
  !============
  implicit none
  integer ialgu,iabs,ix,ixp,itest,i
  !============
  real(fp) :: ztolr,ztola,zh,xtest
  !============
  integer nx                        ! length of x, .ge.4
  real(fp) :: x(nx)                        ! x axis vector, strict ascending
  !
  integer iper                      ! =1 if x is periodic
  integer imsg                      ! =1 for range error messages
  integer itol                      ! =1 to specify tolerance, else default
  !
  ! default tolerance is 5.0e-7
  !
  real(fp) :: ztol                         ! range tolerance, if itol=1
  !
  ! lookup algorithm selection for uneven grids:
  !
  integer ialg                      ! = +/- 1:  <1/h(j)> "Newton" method
  !                                       ! = +/- 2:  binary search
  !                                       ! = +/- 3:  indexing function
  !
  !       to use past search result as init. cond. for next search, give
  !       ialg .lt. 0; otherwise give ialg .gt. 0.
  !
  ! output arguments:
  !
  real(fp) :: xpkg(nx,4)                   ! xpkg, as described above
  integer ier                       ! completion code, 0=OK
  !
  !------------------------------------------------
  !
  if(nx.lt.2) then
     write(6,*) ' %genxpkg:  nx.ge.2 required!'
     ier=1
     return
  else
     ier=0
  end if
  !
  ialgu=ialg
  if(ialgu.eq.0) ialgu=3
  if(iabs(ialgu).gt.3) ialgu=3
  !
  !  get tolerance for end point range check & even spacing check
  !
  if(itol.eq.1) then
     ztolr=abs(ztol)
  else
     ztolr=5.0E-7_fp
  end if
  !
  ztola=max(abs(x(1)),abs(x(nx)))*ztolr
  !
  !  assume even spacing for now...
  !
  if(nx.ge.3) then
     xpkg(3,4)=0.0_fp
  end if
  !
  !  mark if x axis is a periodic coordinate
  !
  if(iper.eq.1) then
     xpkg(2,4)=1.0_fp
  else
     xpkg(2,4)=0.0_fp
  end if
  !
  !  store tolerance parameter
  !
  xpkg(1,4)=ztola
  !
  !  mark if messages are to be written if range lookup errors occur
  !
  if(imsg.eq.1) then
     continue                       ! xpkg(1,4) left .gt. 0.0
  else
     xpkg(1,4)=-xpkg(1,4)
  end if
  !
  !  OK check linearity and spacing
  !
  ier=0
  !
  xpkg(nx,2)=(x(nx)-x(1))/(nx-1)    ! average spacing
  !
  do ix=1,nx
     xpkg(ix,1)=x(ix)
     if((ier.eq.0).and.(ix.lt.nx)) then
        if(x(ix+1).le.x(ix)) then
           ier=1
           write(6,*) ' %genxpkg:  x axis not strict ascending!'
        else
           zh=x(ix+1)-x(ix)
           !
           !  checking even spacing now...
           !
           if(nx.ge.3) then
              if(abs(zh-xpkg(nx,2)).gt.ztola) xpkg(3,4)=1.0_fp
           end if
           !
           xpkg(ix,2)=zh
           xpkg(ix,3)=1.0_fp/zh
           !
        end if
     end if
  end do
  !
  if(ier.ne.0) return
  !
  !  OK, store inverse average spacing
  !
  xpkg(nx,3)=1.0_fp/xpkg(nx,2)
  !
  !  if even spacing is detected, redefine x axis slightly, for
  !  improved regularity of behaviour
  !
  if(nx.ge.3) then
     !
     if(xpkg(3,4).eq.0.0_fp) then
        do ix=1,nx-2
           ixp=ix+1
           xpkg(ixp,1)=xpkg(1,1)+ix*xpkg(nx,2)
        end do
        if(nx.gt.3) then
           if(ialgu.lt.0) then
              xpkg(4,4)=1.0_fp ! check init guess
           else
              xpkg(4,4)=0.0_fp
           end if
        end if
     end if
     !
     !  if uneven spacing is detected, must use an algorithm...
     !
     if(xpkg(3,4).ne.0.0_fp) then
        xpkg(3,4)=abs(ialgu)
        if(nx.gt.3) then
           if(ialgu.lt.0) then
              xpkg(4,4)=1.0_fp ! check init guess
           else
              xpkg(4,4)=0.0_fp
           end if
        end if
        !
        if(abs(ialgu).eq.3) then
           !
           !  construct a piecewise linear indexing function
           !
           xpkg(1,2)=1.0_fp
           xtest=xpkg(1,1)
           itest=1
           do i=2,nx-1
              xtest=xtest+xpkg(nx,2) ! x1 + (i-1)*<h>
10            continue
              if((xpkg(itest,1).le.xtest).and. &
                   (xtest.le.xpkg(itest+1,1))) then
                 xpkg(i,2)=itest+(xtest-xpkg(itest,1))/ &
                      (xpkg(itest+1,1)-xpkg(itest,1))
              else
                 itest=itest+1
                 go to 10
              end if
           end do
           !
           !  (implicitly, xpkg(nx,2)=nx*1.0; but we leave <h> in xpkg(nx,2)
           !   instead...)
           !
        end if
     end if
     !
  end if  ! nx.ge.3
  !
  return
end subroutine genxpkg

subroutine mkbicub(x,nx,y,ny,f,nf2, &
  ibcxmin,bcxmin,ibcxmax,bcxmax, &
  ibcymin,bcymin,ibcymax,bcymax, &
  ilinx,iliny,ier)
!
!  setup bicubic spline, dynamic allocation of workspace
!  fortran-90 fixed form source
!
!  ** see mkbicubw.f90 **
!  arguments are identical, except the workspace argument is
!  omitted (created here)
!
!  --NOTE-- dmc 22 Feb 2004 -- rewrite for direct calculation of
!  coefficients, to avoid large transient use of memory.  mkbicubw is
!  no longer called (but reference to mkbicubw comments is still correct).
!
!  abbreviated description of arguments, mkbicubw has details:
!
!
!  input:
implicit none
integer nx                        ! length of x vector
integer ny                        ! length of y vector
real(fp) :: x(nx)                        ! x vector, strict ascending
real(fp) :: y(ny)                        ! y vector, strict ascending
!
integer nf2                       ! 2nd dimension of f, nf2.ge.nx
!  input/output:
real(fp) :: f(4,nf2,ny)                  ! data & spline coefficients
!
!  input bdy condition data:
integer ibcxmin                   ! bc flag for x=xmin
real(fp) :: bcxmin(ny)                   ! bc data vs. y at x=xmin
integer ibcxmax                   ! bc flag for x=xmax
real(fp) :: bcxmax(ny)                   ! bc data vs. y at x=xmax
!
integer ibcymin                   ! bc flag for y=ymin
real(fp) :: bcymin(nx)                   ! bc data vs. x at y=ymin
integer ibcymax                   ! bc flag for y=ymax
real(fp) :: bcymax(nx)                   ! bc data vs. x at y=ymax
!
!  output linear grid flags and completion code (ier=0 is normal):
!
integer ilinx                     ! =1: x grid is "nearly" equally spaced
integer iliny                     ! =1: y grid is "nearly" equally spaced
!
integer ier                       ! =0 on exit if there is no error.
!
!-----------------------
integer ierx,iery
!
real(fp), dimension(:,:), allocatable :: fwk
real(fp) :: zbcmin,zbcmax
integer ix,iy,ibcmin,ibcmax
!
real(fp), dimension(:,:,:), allocatable :: fcorr
integer iflg2
real(fp) :: zdiff(2),hy
!
!-----------------------
!
!  see if 2nd pass is needed due to inhomogeneous d/dy bdy cond.
!
iflg2=0
if(ibcymin.ne.-1) then
  if((ibcymin.eq.1).or.(ibcymin.eq.2)) then
     do ix=1,nx
        if (bcymin(ix).ne.0.0_fp) iflg2=1
     end do
  end if
  if((ibcymax.eq.1).or.(ibcymax.eq.2)) then
     do ix=1,nx
        if (bcymax(ix).ne.0.0_fp) iflg2=1
     end do
  end if
end if
!
!  check boundary condition specifications
!
ier=0
!
call ibc_ck(ibcxmin,'bcspline','xmin',-1,7,ier)
if(ibcxmin.ge.0) call ibc_ck(ibcxmax,'bcspline','xmax',0,7,ier)
call ibc_ck(ibcymin,'bcspline','ymin',-1,7,ier)
if(ibcymin.ge.0) call ibc_ck(ibcymax,'bcspline','ymax',0,7,ier)
!
!  check ilinx & x vector
!
call splinck(x,nx,ilinx,1.0E-3_fp,ierx)
if(ierx.ne.0) ier=2
!
if(ier.eq.2) then
  write(6,'('' ?bcspline:  x axis not strict ascending'')')
end if
!
!  check iliny & y vector
!
call splinck(y,ny,iliny,1.0E-3_fp,iery)
if(iery.ne.0) ier=3
!
if(ier.eq.3) then
  write(6,'('' ?bcspline:  y axis not strict ascending'')')
end if
!
if(ier.ne.0) return
!
!------------------------------------
allocate(fwk(2,max(nx,ny)))
!
!  evaluate fxx (spline in x direction)
!
zbcmin=0
zbcmax=0
do iy=1,ny
  fwk(1,1:nx) = f(1,1:nx,iy)
  if((ibcxmin.eq.1).or.(ibcxmin.eq.2)) zbcmin=bcxmin(iy)
  if((ibcxmax.eq.1).or.(ibcxmax.eq.2)) zbcmax=bcxmax(iy)
  call mkspline(x,nx,fwk, &
       ibcxmin,zbcmin,ibcxmax,zbcmax,ilinx,ier)
  if(ier.ne.0) return
  f(2,1:nx,iy)=fwk(2,1:nx)
end do
!
!  evaluate fyy (spline in y direction)
!  use homogeneous boundary condition; correction done later if necessary
!
zbcmin=0
zbcmax=0
ibcmin=ibcymin
ibcmax=ibcymax
do ix=1,nx
  fwk(1,1:ny) = f(1,ix,1:ny)
  if(iflg2.eq.1) then
     if((ibcymin.eq.1).or.(ibcymin.eq.2)) ibcmin=0
     if((ibcymax.eq.1).or.(ibcymax.eq.2)) ibcmax=0
  end if
  call mkspline(y,ny,fwk, &
       ibcmin,zbcmin,ibcmax,zbcmax,iliny,ier)
  if(ier.ne.0) return
  f(3,ix,1:ny)=fwk(2,1:ny)
end do
!
!  evaluate fxxyy (spline fxx in y direction; BC simplified; avg
!  d2(d2f/dx2)/dy2 and d2(df2/dy2)/dx2
!
zbcmin=0
zbcmax=0
ibcmin=ibcymin
ibcmax=ibcymax
do ix=1,nx
  fwk(1,1:ny) = f(2,ix,1:ny)
  if(iflg2.eq.1) then
     if((ibcymin.eq.1).or.(ibcymin.eq.2)) ibcmin=0
     if((ibcymax.eq.1).or.(ibcymax.eq.2)) ibcmax=0
  end if
  call mkspline(y,ny,fwk, &
       ibcmin,zbcmin,ibcmax,zbcmax,iliny,ier)
  if(ier.ne.0) return
  f(4,ix,1:ny)= fwk(2,1:ny)
end do
!
if(iflg2.eq.1) then
  allocate(fcorr(2,nx,ny))
  !
  !  correct for inhomogeneous y boundary condition
  !
  do ix=1,nx
     !  the desired inhomogenous BC is the difference btw the
     !  requested derivative (1st or 2nd) and the current value

     zdiff(1)=0.0_fp
     if(ibcymin.eq.1) then
        hy=y(2)-y(1)
        zdiff(1)=(f(1,ix,2)-f(1,ix,1))/hy + &
             hy*(-2*f(3,ix,1)-f(3,ix,2))/6
        zdiff(1)=bcymin(ix)-zdiff(1)
     else if(ibcymin.eq.2) then
        zdiff(1)=bcymin(ix)-f(3,ix,1)
     end if

     zdiff(2)=0.0_fp
     if(ibcymax.eq.1) then
        hy=y(ny)-y(ny-1)
        zdiff(2)=(f(1,ix,ny)-f(1,ix,ny-1))/hy + &
             hy*(2*f(3,ix,ny)+f(3,ix,ny-1))/6
        zdiff(2)=bcymax(ix)-zdiff(2)
     else if(ibcymax.eq.2) then
        zdiff(2)=bcymax(ix)-f(3,ix,ny)
     end if
     !
     fwk(1,1:ny)=0.0_fp  ! values are zero; only BC is not
     call mkspline(y,ny,fwk,ibcymin,zdiff(1),ibcymax,zdiff(2), &
          iliny,ier)
     if(ier.ne.0) return
     fcorr(1,ix,1:ny)=fwk(2,1:ny)  ! fyy-correction
  end do
  !
  zbcmin=0
  zbcmax=0
  do iy=1,ny
     fwk(1,1:nx)=fcorr(1,1:nx,iy)
     call mkspline(x,nx,fwk,ibcxmin,zbcmin,ibcxmax,zbcmax, &
          ilinx,ier)
     if(ier.ne.0) return
     fcorr(2,1:nx,iy)=fwk(2,1:nx)  ! fxxyy-correction
  end do
  !
  f(3:4,1:nx,1:ny)=f(3:4,1:nx,1:ny)+fcorr(1:2,1:nx,1:ny)
  !
  deallocate(fcorr)
end if
!
!  correction spline -- f=fxx=zero; fyy & fxxyy are affected
!
deallocate(fwk)
!------------------------------------
!
!  thats all
!
return
end subroutine mkbicub

subroutine mkspline(x,nx,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax, &
  ilinx,ier)
!
!  make a 2-coefficient 1d spline
!
!  only 2 coefficients, the data and its 2nd derivative, are needed to
!  fully specify a spline.  See e.g. Numerical Recipies in Fortran-77
!  (2nd edition) chapter 3, section on cubic splines.
!
!  input:
!============
implicit none
integer i,inwk
!============
real(fp) :: bcxmax
!============
integer nx                        ! no. of data points
real(fp) :: x(nx)                        ! x axis data, strict ascending order
!
!  input/output:
real(fp) :: fspl(2,nx)                   ! f(1,*):  data in; f(2,*):  coeffs out
!
!     f(1,j) = f(x(j))  on input (unchanged on output)
!     f(2,j) = f''(x(j)) (of interpolating spline) (on output).
!
!  ...boundary conditions...
!
!  input:
!
integer ibcxmin                   ! b.c. flag @ x=xmin=x(1)
real(fp) :: bcxmin                       ! b.c. data @xmin
!
integer ibcxmax                   ! b.c. flag @ x=xmax=x(nx)
!
!  ibcxmin=-1 -- periodic boundary condition
!                (bcxmin,ibcxmax,bcxmax are ignored)
!
!                the output spline s satisfies
!                s'(x(1))=s'(x(nx)) ..and.. s''(x(1))=s''(x(nx))
!
!  if non-periodic boundary conditions are used, then the xmin and xmax
!  boundary conditions can be specified independently:
!
!  ibcxmin (ibcxmax) = 0 -- this specifies a "not a knot" boundary
!                condition, see "cubsplb.f90".  This is a common way
!                for inferring a "good" spline boundary condition
!                automatically from data in the vicinity of the
!                boundary.  ... bcxmin (bcxmax) are ignored.
!
!  ibcxmin (ibcxmax) = 1 -- boundary condition is to have s'(x(1))
!                ( s'(x(nx)) ) match the passed value bcxmin (bcxmax).
!
!  ibcxmin (ibcxmax) = 2 -- boundary condition is to have s''(x(1))
!                ( s''(x(nx)) ) match the passed value bcxmin (bcxmax).
!
!  ibcxmin (ibcxmax) = 3 -- boundary condition is to have s'(x(1))=0.0
!                ( s'(x(nx))=0.0 )
!
!  ibcxmin (ibcxmax) = 4 -- boundary condition is to have s''(x(1))=0.0
!                ( s''(x(nx))=0.0 )
!
!  ibcxmin (ibcxmax) = 5 -- boundary condition is to have s'(x(1))
!                ( s'(x(nx)) ) match the 1st divided difference
!                e.g. at x(1):  d(1)/h(1), where
!                           d(j)=f(1,j+1)-f(1,j)
!                           h(j)=x(j+1)-x(j)
!
!  ibcxmin (ibcxmax) = 6 -- BC is to have s''(x(1)) ( s''(x(nx)) )
!                match the 2nd divided difference
!                e.g. at x(1):
!                     e(1) = [d(2)/h(2) - d(1)/h(1)]/(0.5*(h(1)+h(2)))
!
!  ibcxmin (ibcxmax) = 7 -- BC is to have s'''(x(1)) ( s'''(x(nx)) )
!                match the 3rd divided difference
!                e.g. at x(1): [e(2)-e(1)]/(0.33333*(h(1)+h(2)+h(3)))
!
!  output:
!
integer ilinx                     ! =1: hint, x axis is ~evenly spaced
!
!  let dx[avg] = (x(nx)-x(1))/(nx-1)
!  let dx[j] = x(j+1)-x(j), for all j satisfying 1.le.j.lt.nx
!
!  if for all such j, abs(dx[j]-dx[avg]).le.(1.0e-3*dx[avg]) then
!  ilinx=1 is returned, indicating the data is (at least nearly)
!  evenly spaced.  Even spacing is useful, for speed of zone lookup,
!  when evaluating a spline.
!
!  if the even spacing condition is not satisfied, ilinx=2 is returned.
!
integer ier                       ! exit code, 0=OK
!
!  an error code is returned if the x axis is not strict ascending,
!  or if nx.lt.4, or if an invalid boundary condition specification was
!  input.
!
!------------------------------------
!
!  this routine calls traditional 4 coefficient spline software, and
!  translates the result to 2 coefficient form.
!
!  this could be done more efficiently but we decided out of conservatism
!  to use the traditional software.
!
!------------------------------------
!  workspaces -- f90 dynamically allocated
!
real(fp), dimension(:,:), allocatable :: fspl4 ! traditional 4-spline
real(fp), dimension(:), allocatable :: wk ! cspline workspace
!
!------------------------------------
!
allocate(fspl4(4,nx),wk(nx))
!
!  make the traditional call
!
do i=1,nx
  fspl4(1,i)=fspl(1,i)
  fspl(2,i)=0.0_fp  ! for now
end do
!
inwk=nx
!
!  boundary conditions imposed by cspline...
!
call cspline(x,nx,fspl4,ibcxmin,bcxmin,ibcxmax,bcxmax, &
    wk,inwk,ilinx,ier)
!
if(ier.eq.0) then
  !
  !  copy the output -- careful of end point.
  !
  do i=1,nx-1
     fspl(2,i)=2.0_fp*fspl4(3,i)
  end do
  fspl(2,nx)=2.0_fp*fspl4(3,nx-1) + &
       (x(nx)-x(nx-1))*6.0_fp*fspl4(4,nx-1)
end if
!
deallocate(fspl4,wk)
!
return
end subroutine mkspline

!  cspline -- dmc 15 Feb 1999
!
!  a standard interface to the 1d spline setup routine
!    modified dmc 3 Mar 2000 -- to use Wayne Houlberg's v_spline code.
!    new BC options added.
!
subroutine cspline(x,nx,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax, &
  wk,iwk,ilinx,ier)
!
!============
implicit none
integer iwk,nx,ierx,inum,i
!============
real(fp) :: half,sixth
!============
real(fp) :: x(nx)                        ! x axis (in)
real(fp) :: fspl(4,nx)                   ! spline data (in/out)
integer ibcxmin                   ! x(1) BC flag (in, see comments)
real(fp) :: bcxmin                       ! x(1) BC data (in, see comments)
integer ibcxmax                   ! x(nx) BC flag (in, see comments)
real(fp) :: bcxmax                       ! x(nx) BC data (in, see comments)
real(fp) :: wk(iwk)                      ! workspace of size at least nx
integer ilinx                     ! even spacing flag (out)
integer ier                       ! output, =0 means OK
!
!  ** note wk(...) array is not used unless ibcxmin=-1 (periodic spline
!  evaluation)
!
!  this routine computes spline coefficients for a 1d spline --
!  evaluation of the spline can be done by cspeval.f90 subroutines
!  or directly by inline code.
!
!  the input x axis x(1...nx) must be strictly ascending, i.e.
!  x(i+1).gt.x(i) is required for i=1 to nx-1.  This is checked and
!  ier=1 is set and the routine exits if the test is not satisfied.
!
!  on output, ilinx=1 is set if, to a reasonably close tolerance,
!  all grid spacings x(i+1)-x(i) are equal.  This allows a speedier
!  grid lookup algorithm on evaluation of the spline.  If on output
!  ilinx=2, this means the spline x axis is not evenly spaced.
!
!  the input data for the spline are given in f[j] = fspl(1,j).  The
!  output data are the spline coefficients fspl(2,j),fspl(3,j), and
!  fspl(4,j), j=1 to nx.  The result is a spline s(x) satisfying the
!  boundary conditions and with the properties
!
!     s(x(j)) = fspl(1,j)
!     s'(x) is continuous even at the grid points x(j)
!     s''(x) is continuous even at the grid points x(j)
!
!  the formula for evaluation of s(x) is:
!
!     let dx = x-x(i), where x(i).le.x.le.x(i+1).  Then,
!     s(x)=fspl(1,i) + dx*(fspl(2,i) +dx*(fspl(3,i) + dx*fspl(4,i)))
!
!  ==>boundary conditions.  Complete specification of a 1d spline
!  requires specification of boundary conditions at x(1) and x(nx).
!
!  this routine provides 4 options:
!
! -1 ***** PERIODIC BC
!  ibcxmin=-1  --  periodic boundary condition.  This means the
!    boundary conditions s'(x(1))=s'(x(nx)) and s''(x(1))=s''(x(nx))
!    are imposed.  Note that s(x(1))=s(x(nx)) (i.e. fspl(1,1)=fspl(1,nx))
!    is not required -- that is determined by the fspl array input data.
!    The periodic boundary condition is to be preferred for periodic
!    data.  When splining periodic data f(x) with period P, the relation
!    x(nx)=x(1)+n*P, n = the number of periods (usually 1), should hold.
!    (ibcxmax, bcxmin, bcxmax are ignored).
!
!  if a periodic boundary condition is set, this covers both boundaries.
!  for the other types of boundary conditions, the type of condition
!  chosen for the x(1) boundary need not be the same as the type chosen
!  for the x(nx) boundary.
!
!  0 ***** NOT A KNOT BC
!  ibcxmin=0 | ibcxmax=0 -- this specifies a "not a knot" boundary
!    condition -- see cubsplb.f90.  This is a common way for inferring
!    a "good" spline boundary condition automatically from data in the
!    vicinity of the boundary.  (bcxmin | bcxmax are ignored).
!
!  1 ***** BC:  SPECIFIED SLOPE
!  ibcxmin=1 | ibcxmax=1 -- boundary condition is to have s'(x(1)) |
!    s'(x(nx)) match the passed value (bcxmin | bcxmax).
!
!  2 ***** BC:  SPECIFIED 2nd DERIVATIVE
!  ibcxmin=2 | ibcxmax=2 -- boundary condition is to have s''(x(1)) |
!    s''(x(nx)) match the passed value (bcxmin | bcxmax).
!
!  3 ***** BC:  SPECIFIED SLOPE = 0.0
!  ibcxmin=3 | ibcxmax=3 -- boundary condition is to have s'(x(1)) |
!    s'(x(nx)) equal to ZERO.
!
!  4 ***** BC:  SPECIFIED 2nd DERIVATIVE = 0.0
!  ibcxmin=4 | ibcxmax=4 -- boundary condition is to have s''(x(1)) |
!    s''(x(nx)) equal to ZERO.
!
!  5 ***** BC:  1st DIVIDED DIFFERENCE
!  ibcxmin=5 | ibcxmax=5 -- boundary condition is to have s'(x(1)) |
!    s'(x(nx)) equal to the slope from the 1st|last 2 points
!
!  6 ***** BC:  2nd DIVIDED DIFFERENCE
!  ibcxmin=6 | ibcxmax=6 -- boundary condition is to have s''(x(1)) |
!    s''(x(nx)) equal to the 2nd derivative from the 1st|last 3 points
!
!  7 ***** BC:  3rd DIVIDED DIFFERENCE
!  ibcxmin=7 | ibcxmax=7 -- boundary condition is to have s'''(x(1)) |
!    s'''(x(nx)) equal to the 3rd derivative from the 1st|last 4 points
!
!---------------------------------------------------------------------
half = 0.5_fp
sixth = 0.166666666666666667_fp
!
!  error checks
!
ier = 0
if(nx.lt.2) then
  write(6,'('' ?cspline:  at least 2 x points required.'')')
  ier=1
end if
call ibc_ck(ibcxmin,'cspline','xmin',-1,7,ier)
if(ibcxmin.ge.0) call ibc_ck(ibcxmax,'cspline','xmax',0,7,ier)
!
!  x axis check
!
call splinck(x,nx,ilinx,1.0E-3_fp,ierx)
if(ierx.ne.0) ier=2
!
if(ier.eq.2) then
  write(6,'('' ?cspline:  x axis not strict ascending'')')
end if
!
if(ibcxmin.eq.-1) then
  inum=nx
  if(iwk.lt.inum) then
     write(6,1009) inum,iwk,nx
1009    format( &
          ' ?cspline:  workspace too small.  need:  ',i6,' got:  ',i6/ &
          '  (need = nx, nx=',i6)
     ier=3
  end if
end if
!
if(ier.ne.0) return
!
!  OK -- evaluate spline
!
if(ibcxmin.eq.1) then
  fspl(2,1)=bcxmin
else if(ibcxmin.eq.2) then
  fspl(3,1)=bcxmin
end if
!
if(ibcxmax.eq.1) then
  fspl(2,nx)=bcxmax
else if(ibcxmax.eq.2) then
  fspl(3,nx)=bcxmax
end if
!
call v_spline(ibcxmin,ibcxmax,nx,x,fspl,wk)
!
do i=1,nx
  fspl(3,i)=half*fspl(3,i)
  fspl(4,i)=sixth*fspl(4,i)
end do
!
return
end subroutine cspline

subroutine v_spline(k_bc1,k_bcn,n,x,f,wk)
  !***********************************************************************
  !V_SPLINE evaluates the coefficients for a 1d cubic interpolating spline
  !References:
  !  Forsythe, Malcolm, Moler, Computer Methods for Mathematical
  !    Computations, Prentice-Hall, 1977, p.76
  !  Engeln-Muellges, Uhlig, Numerical Algorithms with Fortran, Springer,
  !    1996, p.251
  !  W.A.Houlberg, D.McCune 3/2000
  !Input:
  !  k_bc1-option for BC at x1 = x(1)
  !       =-1 periodic, ignore k_bcn
  !       =0 not-a-knot
  !       =1 s'(x1) = input value of f(2,1)
  !       =2 s''(x1) = input value of f(3,1)
  !       =3 s'(x1) = 0.0
  !       =4 s''(x1) = 0.0
  !       =5 match first derivative to first 2 points
  !       =6 match second derivative to first 3 points
  !       =7 match third derivative to first 4 points
  !       =else use not-a-knot
  !  k_bcn-option for boundary condition at xn = x(n)
  !       =0 not-a-knot
  !       =1 s'(xn) = input value of f(2,n)
  !       =2 s''(xn) = input value of f(3,n)
  !       =3 s'(xn) = 0.0
  !       =4 s''(xn) = 0.0
  !       =5 match first derivative to last 2 points
  !       =6 match second derivative to lasst 3 points
  !       =7 match third derivative to last 4 points
  !       =else use knot-a-knot
  !  n-number of data points or knots-(n.ge.2)
  !  x(n)-abscissas of the knots in strictly increasing order
  !  f(1,i)-ordinates of the knots
  !  f(2,1)-input value of s'(x1) for k_bc1=1
  !  f(2,n)-input value of s'(xn) for k_bcn=1
  !  f(3,1)-input value of s''(x1) for k_bc1=2
  !  f(3,n)-input value of s''(xn) for k_bcn=2
  !  wk(n)-scratch work area for periodic BC
  !Output:
  !  f(2,i)=s'(x(i))
  !  f(3,i)=s''(x(i))
  !  f(4,i)=s'''(x(i))
  !Comments:
  !  s(x)=f(1,i)+f(2,i)*(x-x(i))+f(3,i)*(x-x(i))**2/2!
  !       +f(4,i)*(x-x(i))**3/3! for x(i).le.x.le.x(i+1)
  !  W_SPLINE can be used to evaluate the spline and its derivatives
  !  The cubic spline is twice differentiable (C2)
  !
  !  modifications -- dmc 18 Feb 2010:
  !    Deal with n.lt.4 -- the general tridiagonal spline method
  !    does not have the right formulation for n.eq.3 "not a knot" or periodic
  !    boundary conditions, nor for n.eq.2 with any boundary conditions.
  !
  !    Apply boundary conditions even for n=2, when the "spline" is really
  !    just a single cubic polynomial.
  !    In this case, several boundary condition (BC) options are mapped to
  !      BC option 5, 1st divided difference.  If 5 is used on both sides
  !      of an n=2 "spline" you get a linear piece which is what the old
  !      code always gave, regardless of BC option settings.  The following
  !      BC controls are mapped to 5:
  !        periodic (-1)
  !        not a knot (0) (for n=2 no grid point exists for knot location).
  !        option (5) is preserved
  !        options 6 and 7 -- mapped to (5); higher divided differences
  !          need n>2; in fact 7 needs n>3; for n=3 option 6 is substituted.
  !
  !    The boundary condition screening is done at the start of the code;
  !    passed controls k_bc1 and k_bcn are mapped to i_bc1 and i_bcn.
  !
  !    ALSO: for n=3, "not a knot" from both left and right needs special
  !      interpretation, since the 2 boundary conditions overlap.  The
  !      chosen interpretation is a parabolic fit to the 3 given data points.
  !      and so f''' = 0 and f'' = constant.  If "not a knot" is used on
  !      one side only, the solution is a single cubic piece and special
  !      code is also needed.
  !    ALSO: for n=3, "periodic" boundary condition needs special code; this
  !      is added.
  !
  !  bugfixes -- dmc 24 Feb 2004:
  !    (a) fixed logic for not-a-knot:
  !          !    Set f(3,1) for not-a-knot
  !                    IF(k_bc1.le.0.or.k_bc1.gt.7) THEN ...
  !        instead of
  !          !    Set f(3,1) for not-a-knot
  !                    IF(k_bc1.le.0.or.k_bc1.gt.5) THEN ...
  !        and similarly for logic after cmt
  !          !    Set f(3,n) for not-a-knot
  !        as required since k_bc*=6 and k_bc*=7 are NOT not-a-knot BCs.
  !
  !    (b) the BCs to fix 2nd derivative at end points did not work if that
  !        2nd derivative were non-zero.  The reason is that in those cases
  !        the off-diagonal matrix elements nearest the corners are not
  !        symmetric; i.e. elem(1,2).ne.elem(2,1) and
  !        elem(n-1,n).ne.elem(n,n-1) where I use "elem" to refer to
  !        the tridiagonal matrix elements.  The correct values for the
  !        elements is:   elem(1,2)=0, elem(2,1)=x(2)-x(1)
  !                       elem(n,n-1)=0, elem(n-1,n)=x(n)-x(n-1)
  !        the old code in effect had these as all zeroes.  Since this
  !        meant the wrong set of linear equations was solved, the
  !        resulting spline had a discontinuity in its 1st derivative
  !        at x(2) and x(n-1).  Fixed by introducing elem21 and elemnn1
  !        to represent the non-symmetric lower-diagonal values.  Since
  !        elem21 & elemnn1 are both on the lower diagonals, logic to
  !        use them occurs in the non-periodic forward elimination loop
  !        only.  DMC 24 Feb 2004.
  !***********************************************************************
  !Declaration of input variables
  implicit none
  integer :: k_bc1, k_bcn, n
  real(fp) ::  x(*), wk(*), f(4,*)
  !Declaration in local variables
  integer :: i, ib, imax, imin
  real(fp) :: a1, an, b1, bn, q, t, hn
  real(fp) ::       elem21,                  elemnn1

  integer :: i_bc1,i_bcn  ! screened BC controls

  integer :: iord1,iord2  ! used for n=2 only
  real(fp) :: h,f0,fh         ! used for n=2,3 only
  real(fp) :: h1,h2,h3,dels   ! used for n=3 special cases
  real(fp) :: f1,f2,f3,aa,bb  ! used for n=3

  integer :: i3knots      ! for n=3, number of not-a-knot BCs
  integer :: i3perio      ! for n=3, periodic BC

  !------------------------------------------------------------
  !  screen the BC options (DMC Feb. 2010...)

  i_bc1=k_bc1
  i_bcn=k_bcn

  if((i_bc1.lt.-1).or.(i_bc1.gt.7)) i_bc1=0  ! outside [-1:7] -> not-a-knot
  if((i_bcn.lt.0).or.(i_bcn.gt.7)) i_bcn=0   ! outside [0:7] -> not-a-knot

  if(i_bc1.eq.-1) i_bcn=-1  ! periodic BC

  i3knots=0
  i3perio=0
  if(n.eq.3) then
     i_bc1=min(6,i_bc1)
     i_bcn=min(6,i_bcn)
     if(i_bc1.eq.0) i3knots = i3knots + 1
     if(i_bcn.eq.0) i3knots = i3knots + 1
     if(i_bc1.eq.-1) i3perio = 1
  end if

  if(n.eq.2) then
     if(i_bc1.eq.-1) then
        i_bc1=5
        i_bcn=5
     end if
     if((i_bc1.eq.0).or.(i_bc1.gt.5)) i_bc1=5
     if((i_bcn.eq.0).or.(i_bcn.gt.5)) i_bcn=5

     if((i_bc1.eq.1).or.(i_bc1.eq.3).or.(i_bc1.eq.5)) then
        iord1=1  ! 1st derivative match on LHS
     else
        iord1=2  ! 2nd derivative match on LHS
     end if

     if((i_bcn.eq.1).or.(i_bcn.eq.3).or.(i_bcn.eq.5)) then
        iord2=1  ! 1st derivative match on RHS
     else
        iord2=2  ! 2nd derivative match on RHS
     end if
  end if

  !Set default range
  imin=1
  imax=n
  !Set first and second BC values
  a1=0.0_fp
  b1=0.0_fp
  an=0.0_fp
  bn=0.0_fp
  IF(i_bc1.eq.1) THEN
     a1=f(2,1)
  ELSEIF(i_bc1.eq.2) THEN
     b1=f(3,1)
  ELSEIF(i_bc1.eq.5) THEN
     a1=(f(1,2)-f(1,1))/(x(2)-x(1))
  ELSEIF(i_bc1.eq.6) THEN
     b1=2.0_fp*((f(1,3)-f(1,2))/(x(3)-x(2)) &
          -(f(1,2)-f(1,1))/(x(2)-x(1)))/(x(3)-x(1))
  end if
  IF(i_bcn.eq.1) THEN
     an=f(2,n)
  ELSEIF(i_bcn.eq.2) THEN
     bn=f(3,n)
  ELSEIF(i_bcn.eq.5) THEN
     an=(f(1,n)-f(1,n-1))/(x(n)-x(n-1))
  ELSEIF(i_bcn.eq.6) THEN
     bn=2.0_fp*((f(1,n)-f(1,n-1))/(x(n)-x(n-1)) &
          -(f(1,n-1)-f(1,n-2))/(x(n-1)-x(n-2)))/(x(n)-x(n-2))
  end if
  !Clear f(2:4,n)
  f(2,n)=0.0_fp
  f(3,n)=0.0_fp
  f(4,n)=0.0_fp
  IF(n.eq.2) THEN
     if((i_bc1.eq.5).and.(i_bcn.eq.5)) then
        !Coefficients for n=2 (this was the original code)
        f(2,1)=(f(1,2)-f(1,1))/(x(2)-x(1))
        f(3,1)=0.0_fp
        f(4,1)=0.0_fp
        f(2,2)=f(2,1)
        f(3,2)=0.0_fp
        f(4,2)=0.0_fp
     else if((iord1.eq.1).and.(iord2.eq.1)) then
        ! LHS: match a1 for 1st deriv; RHS: match an for 1st deriv.
        f(2,1)=a1
        f(2,2)=an
        h = (x(2)-x(1))
        f0 = f(1,1)
        fh = f(1,2)

        ! setting xx = x-x(1),
        ! f = c1*xx**3 + c2*xx**2 + a1*xx + f0
        !   -->  c1*h**3   + c2*h**2 = fh - f0 - a1*h
        !   and  3*c1*h**2 + 2*c2*h  = an - a1
        ! f' = 3*c1*xx*2 + 2*c2*xx + a1
        ! f'' = 6*c1*xx + 2*c2
        ! f''' = 6*c1

        ! solve 2x2 system for c1 -> f(4,1)/6 and c2 -> f(3,1)/2

        f(3,1) = (3*(fh-f0)/(h*h) - (2*a1 + an)/h)*2       ! 2*c2
        f(4,1) = (-2*(fh-f0)/(h*h*h) + (a1 + an)/(h*h))*6  ! 6*c1

        f(4,2) = f(4,1)
        f(3,2) = f(4,1)*h + f(3,1)

     else if((iord1.eq.1).and.(iord2.eq.2)) then
        ! LHS: match a1 for 1st deriv; RHS: match bn for 2nd deriv.
        f(2,1)=a1
        f(3,2)=bn
        h = (x(2)-x(1))
        f0 = f(1,1)
        fh = f(1,2)

        ! setting xx = x-x(1),
        ! f = c1*xx**3 + c2*xx**2 + a1*xx + f0
        !   -->  c1*h**3   + c2*h**2 = fh - f0 - a1*h
        !   and  6*c1*h    + 2*c2    = bn
        ! f' = 3*c1*xx*2 + 2*c2*xx + a1
        ! f'' = 6*c1*xx + 2*c2
        ! f''' = 6*c1

        ! solve 2x2 system for c1 -> f(4,1)/6 and c2 -> f(3,1)/2

        f(3,1) = (-bn/4 + 3*(fh-f0)/(2*h*h) - 3*a1/(2*h))*2       ! 2*c2
        f(4,1) = (bn/(4*h) - (fh-f0)/(2*h*h*h) + a1/(2*h*h))*6    ! 6*c1

        f(4,2) = f(4,1)
        f(2,2) = f(4,1)*h*h/2 + f(3,1)*h + a1
     else if((iord1.eq.2).and.(iord2.eq.1)) then
        ! LHS: match b1 for 2nd deriv; RHS: match an for 1st deriv.
        f(3,1)=b1
        f(2,2)=an
        h = (x(2)-x(1))
        f0 = f(1,1)
        fh = f(1,2)

        ! setting xx = x-x(1),
        ! f = c1*xx**3 + (b1/2)*xx**2 + c3*xx + f0
        !   -->  c1*h**3   + c3*h = fh - f0 - b1*h**2/2
        !   and  3*c1*h**2 + c3   = an - b1*h
        ! f' = 3*c1*xx*2 + b1*xx + c3
        ! f'' = 6*c1*xx + b1
        ! f''' = 6*c1

        ! solve 2x2 system for c1 -> f(4,1)/6 and c2 -> f(3,1)/2

        f(2,1) = 3*(fh-f0)/(2*h) - b1*h/4 - an/2                  ! c3
        f(4,1) = (an/(2*h*h) - (fh-f0)/(2*h*h*h) - b1/(4*h))*6    ! 6*c1

        f(4,2) = f(4,1)
        f(3,2) = f(4,1)*h + f(3,1)
     else if((iord1.eq.2).and.(iord2.eq.2)) then
        ! LHS: match b1 for 2nd deriv; RHS: match bn for 2nd deriv.
        f(3,1)=b1
        f(3,2)=bn
        h = (x(2)-x(1))
        f0 = f(1,1)
        fh = f(1,2)

        ! setting xx = x-x(1),
        ! f = c1*xx**3 + (b1/2)*xx**2 + c3*xx + f0
        !   -->  c1*h**3   + c3*h = fh - f0 - b1*h**2/2
        !   and  6*c1*h           = bn - b1
        ! f' = 3*c1*xx*2 + b1*xx + c3
        ! f'' = 6*c1*xx + b1
        ! f''' = 6*c1

        ! solve 2x2 system for c1 -> f(4,1)/6 and c2 -> f(3,1)/2

        f(2,1) = (fh-f0)/h -b1*h/3 -bn*h/6                 ! c3
        f(4,1) = (bn-b1)/h                                 ! 6*c1

        f(4,2) = f(4,1)
        f(2,2) = f(4,1)*h*h/2 + b1*h + f(2,1)
     end if

  ELSE IF(i3perio.eq.1) then
     !Special case: nx=3 periodic spline
     h1=x(2)-x(1)
     h2=x(3)-x(2)
     h=h1+h2

     dels=(f(1,3)-f(1,2))/h2 - (f(1,2)-f(1,1))/h1

     f(2,1)= (f(1,2)-f(1,1))/h1 + (h1*dels)/h
     f(3,1)= -6*dels/h
     f(4,1)= 12*dels/(h1*h)

     f(2,2)= (f(1,3)-f(1,2))/h2 - (h2*dels)/h
     f(3,2)= 6*dels/h
     f(4,2)= -12*dels/(h2*h)

     f(2,3)=f(2,1)
     f(3,3)=f(3,1)
     f(4,3)=f(4,1)


  ELSE IF(i3knots.eq.2) then
     !Special case: nx=3, not-a-knot on both sides
     h1=x(2)-x(1)
     h2=x(3)-x(2)
     h=h1+h2
     ! this is just a quadratic fit through 3 pts
     f1=f(1,1)-f(1,2)
     f2=f(1,3)-f(1,2)

     !  quadratic around origin at (x(2),f(1,2))
     !          aa*h1**2 - bb*h1 = f1
     !          aa*h2**2 + bb*h2 = f2

     aa = (f2*h1 + f1*h2)/(h1*h2*h)
     bb = (f2*h1*h1 - f1*h2*h2)/(h1*h2*h)

     f(4,1:3)=0.0_fp  ! f''' = 0 (quadratic polynomial)
     f(3,1:3)=2*aa  ! f'' = const

     f(2,1)=bb-2*aa*h1
     f(2,2)=bb
     f(2,3)=bb+2*aa*h2

  ELSE IF(i3knots.eq.1) then
     !Special cases: nx=3, not-a-knot on single side
     if((i_bc1.eq.1).or.(i_bc1.eq.3).or.(i_bc1.eq.5)) then
        ! f' LHS condition; not-a-knot RHS
        !  a1 = f' LHS BC
        h2=x(2)-x(1)
        h3=x(3)-x(1)

        f2=f(1,2)-f(1,1)
        f3=f(1,3)-f(1,1)

        !  find cubic aa*xx**3 + bb*xx**2 + a1*xx
        !    satisfying aa*h2**3 + bb*h2**2 + a1*h2 = f2
        !           and aa*h3**3 + bb*h3**2 + a1*h3 = f3

        aa=a1/(h2*h3) + f3/(h3*h3*(h3-h2)) - f2/(h2*h2*(h3-h2))
        bb=-a1*(h3*h3-h2*h2)/(h2*h3*(h3-h2)) &
             + f2*h3/(h2*h2*(h3-h2)) - f3*h2/(h3*h3*(h3-h2))

        f(2,1)=a1
        f(3,1)=2*bb
        f(4,1)=6*aa

        f(2,2)=3*aa*h2*h2 + 2*bb*h2 + a1
        f(3,2)=6*aa*h2 + 2*bb
        f(4,2)=6*aa

        f(2,3)=3*aa*h3*h3 + 2*bb*h3 + a1
        f(3,3)=6*aa*h3 + 2*bb
        f(4,3)=6*aa

     else if((i_bc1.eq.2).or.(i_bc1.eq.4).or.(i_bc1.eq.6)) then
        ! f'' LHS condition; not-a-knot RHS
        !  b1 = f'' LHS BC
        h2=x(2)-x(1)
        h3=x(3)-x(1)

        f2=f(1,2)-f(1,1)
        f3=f(1,3)-f(1,1)

        !  find cubic aa*xx**3 + (b1/2)*xx**2 + bb*xx
        !    satisfying aa*h2**3 + bb*h2 = f2 -(b1/2)*h2**2
        !           and aa*h3**3 + bb*h3 = f3 -(b1/2)*h3**2

        aa= -(b1/2)*(h3-h2)/(h3*h3-h2*h2) &
             -f2/(h2*(h3*h3-h2*h2)) + f3/(h3*(h3*h3-h2*h2))
        bb= -(b1/2)*h2*h3*(h3-h2)/(h3*h3-h2*h2) &
             +f2*h3*h3/(h2*(h3*h3-h2*h2)) &
             -f3*h2*h2/(h3*(h3*h3-h2*h2))

        f(2,1)=bb
        f(3,1)=b1
        f(4,1)=6*aa

        f(2,2)=3*aa*h2*h2 + b1*h2 + bb
        f(3,2)=6*aa*h2 + b1
        f(4,2)=6*aa

        f(2,3)=3*aa*h3*h3 + b1*h3 + bb
        f(3,3)=6*aa*h3 + b1
        f(4,3)=6*aa

     else if((i_bcn.eq.1).or.(i_bcn.eq.3).or.(i_bcn.eq.5)) then
        ! f' RHS condition; not-a-knot LHS
        !  an = f' RHS BC
        h2=x(2)-x(3)
        h3=x(1)-x(3)

        f2=f(1,2)-f(1,3)
        f3=f(1,1)-f(1,3)

        !  find cubic aa*xx**3 + bb*xx**2 + an*xx
        !    satisfying aa*h2**3 + bb*h2**2 + an*h2 = f2
        !           and aa*h3**3 + bb*h3**2 + an*h3 = f3

        aa=an/(h2*h3) + f3/(h3*h3*(h3-h2)) - f2/(h2*h2*(h3-h2))
        bb=-an*(h3*h3-h2*h2)/(h2*h3*(h3-h2)) &
             + f2*h3/(h2*h2*(h3-h2)) - f3*h2/(h3*h3*(h3-h2))

        f(2,3)=an
        f(3,3)=2*bb
        f(4,3)=6*aa

        f(2,2)=3*aa*h2*h2 + 2*bb*h2 + an
        f(3,2)=6*aa*h2 + 2*bb
        f(4,2)=6*aa

        f(2,1)=3*aa*h3*h3 + 2*bb*h3 + an
        f(3,1)=6*aa*h3 + 2*bb
        f(4,1)=6*aa

     else if((i_bcn.eq.2).or.(i_bcn.eq.4).or.(i_bcn.eq.6)) then
        ! f'' RHS condition; not-a-knot LHS
        !  bn = f'' RHS BC
        h2=x(2)-x(3)
        h3=x(1)-x(3)

        f2=f(1,2)-f(1,3)
        f3=f(1,1)-f(1,3)

        !  find cubic aa*xx**3 + (bn/2)*xx**2 + bb*xx
        !    satisfying aa*h2**3 + bb*h2 = f2 -(bn/2)*h2**2
        !           and aa*h3**3 + bb*h3 = f3 -(bn/2)*h3**2

        aa= -(bn/2)*(h3-h2)/(h3*h3-h2*h2) &
             -f2/(h2*(h3*h3-h2*h2)) + f3/(h3*(h3*h3-h2*h2))
        bb= -(bn/2)*h2*h3*(h3-h2)/(h3*h3-h2*h2) &
             +f2*h3*h3/(h2*(h3*h3-h2*h2)) &
             -f3*h2*h2/(h3*(h3*h3-h2*h2))

        f(2,3)=bb
        f(3,3)=bn
        f(4,3)=6*aa

        f(2,2)=3*aa*h2*h2 + bn*h2 + bb
        f(3,2)=6*aa*h2 + bn
        f(4,2)=6*aa

        f(2,1)=3*aa*h3*h3 + bn*h3 + bb
        f(3,1)=6*aa*h3 + bn
        f(4,1)=6*aa

     end if
  ELSE IF(n.gt.2) THEN
     !Set up tridiagonal system for A*y=B where y(i) are the second
     !  derivatives at the knots
     !  f(2,i) are the diagonal elements of A
     !  f(4,i) are the off-diagonal elements of A
     !  f(3,i) are the B elements/3, and will become c/3 upon solution
     f(4,1)=x(2)-x(1)
     f(3,2)=(f(1,2)-f(1,1))/f(4,1)
     DO i=2,n-1
        f(4,i)=x(i+1)-x(i)
        f(2,i)=2.0_fp*(f(4,i-1)+f(4,i))
        f(3,i+1)=(f(1,i+1)-f(1,i))/f(4,i)
        f(3,i)=f(3,i+1)-f(3,i)
     end do
     !
     !  (dmc): save now:
     !
     elem21=f(4,1)
     elemnn1=f(4,n-1)
     !
     !  BC's
     !    Left
     IF(i_bc1.eq.-1) THEN
        f(2,1)=2.0_fp*(f(4,1)+f(4,n-1))
        f(3,1)=(f(1,2)-f(1,1))/f(4,1)-(f(1,n)-f(1,n-1))/f(4,n-1)
        wk(1)=f(4,n-1)
        DO i=2,n-3
           wk(i)=0.0_fp
        end do
        wk(n-2)=f(4,n-2)
        wk(n-1)=f(4,n-1)
     ELSEIF(i_bc1.eq.1.or.i_bc1.eq.3.or.i_bc1.eq.5) THEN
        f(2,1)=2.0_fp*f(4,1)
        f(3,1)=(f(1,2)-f(1,1))/f(4,1)-a1
     ELSEIF(i_bc1.eq.2.or.i_bc1.eq.4.or.i_bc1.eq.6) THEN
        f(2,1)=2.0_fp*f(4,1)
        f(3,1)=f(4,1)*b1/3.0_fp
        f(4,1)=0.0_fp  ! upper diagonal only (dmc: cf elem21)
     ELSEIF(i_bc1.eq.7) THEN
        f(2,1)=-f(4,1)
        f(3,1)=f(3,3)/(x(4)-x(2))-f(3,2)/(x(3)-x(1))
        f(3,1)=f(3,1)*f(4,1)**2/(x(4)-x(1))
     ELSE                             ! not a knot:
        imin=2
        f(2,2)=f(4,1)+2.0_fp*f(4,2)
        f(3,2)=f(3,2)*f(4,2)/(f(4,1)+f(4,2))
     end if
     !    Right
     IF(i_bcn.eq.1.or.i_bcn.eq.3.or.i_bcn.eq.5) THEN
        f(2,n)=2.0_fp*f(4,n-1)
        f(3,n)=-(f(1,n)-f(1,n-1))/f(4,n-1)+an
     ELSEIF(i_bcn.eq.2.or.i_bcn.eq.4.or.i_bcn.eq.6) THEN
        f(2,n)=2.0_fp*f(4,n-1)
        f(3,n)=f(4,n-1)*bn/3.0_fp
        !xxx          f(4,n-1)=0.0  ! dmc: preserve f(4,n-1) for back subst.
        elemnn1=0.0_fp  !  lower diaganol only (dmc)
     ELSEIF(i_bcn.eq.7) THEN
        f(2,n)=-f(4,n-1)
        f(3,n)=f(3,n-1)/(x(n)-x(n-2))-f(3,n-2)/(x(n-1)-x(n-3))
        f(3,n)=-f(3,n)*f(4,n-1)**2/(x(n)-x(n-3))
     ELSEIF(i_bc1.ne.-1) THEN         ! not a knot:
        imax=n-1
        f(2,n-1)=2.0_fp*f(4,n-2)+f(4,n-1)
        f(3,n-1)=f(3,n-1)*f(4,n-2)/(f(4,n-1)+f(4,n-2))
     end if
     IF(i_bc1.eq.-1) THEN
        !Solve system of equations for second derivatives at the knots
        !  Periodic BC
        !    Forward elimination
        DO i=2,n-2
           t=f(4,i-1)/f(2,i-1)
           f(2,i)=f(2,i)-t*f(4,i-1)
           f(3,i)=f(3,i)-t*f(3,i-1)
           wk(i)=wk(i)-t*wk(i-1)
           q=wk(n-1)/f(2,i-1)
           wk(n-1)=-q*f(4,i-1)
           f(2,n-1)=f(2,n-1)-q*wk(i-1)
           f(3,n-1)=f(3,n-1)-q*f(3,i-1)
        end do
        !    Correct the n-1 element
        wk(n-1)=wk(n-1)+f(4,n-2)
        !    Complete the forward elimination
        !    wk(n-1) and wk(n-2) are the off-diag elements of the lower corner
        t=wk(n-1)/f(2,n-2)
        f(2,n-1)=f(2,n-1)-t*wk(n-2)
        f(3,n-1)=f(3,n-1)-t*f(3,n-2)
        !    Back substitution
        f(3,n-1)=f(3,n-1)/f(2,n-1)
        f(3,n-2)=(f(3,n-2)-wk(n-2)*f(3,n-1))/f(2,n-2)
        DO ib=3,n-1
           i=n-ib
           f(3,i)=(f(3,i)-f(4,i)*f(3,i+1)-wk(i)*f(3,n-1))/f(2,i)
        end do
        f(3,n)=f(3,1)
     ELSE
        !  Non-periodic BC
        !    Forward elimination
        !    For Not-A-Knot BC the off-diagonal end elements are not equal
        DO i=imin+1,imax
           IF((i.eq.n-1).and.(imax.eq.n-1)) THEN
              t=(f(4,i-1)-f(4,i))/f(2,i-1)
           ELSE
              if(i.eq.2) then
                 t=elem21/f(2,i-1)
              else if(i.eq.n) then
                 t=elemnn1/f(2,i-1)
              else
                 t=f(4,i-1)/f(2,i-1)
              end if
           end if
           IF((i.eq.imin+1).and.(imin.eq.2)) THEN
              f(2,i)=f(2,i)-t*(f(4,i-1)-f(4,i-2))
           ELSE
              f(2,i)=f(2,i)-t*f(4,i-1)
           end if
           f(3,i)=f(3,i)-t*f(3,i-1)
        end do
        !    Back substitution
        f(3,imax)=f(3,imax)/f(2,imax)
        DO ib=1,imax-imin
           i=imax-ib
           IF((i.eq.2).and.(imin.eq.2)) THEN
              f(3,i)=(f(3,i)-(f(4,i)-f(4,i-1))*f(3,i+1))/f(2,i)
           ELSE
              f(3,i)=(f(3,i)-f(4,i)*f(3,i+1))/f(2,i)
           end if
        end do
        !    Reset d array to step size
        f(4,1)=x(2)-x(1)
        f(4,n-1)=x(n)-x(n-1)
        !    Set f(3,1) for not-a-knot
        IF(i_bc1.le.0.or.i_bc1.gt.7) THEN
           f(3,1)=(f(3,2)*(f(4,1)+f(4,2))-f(3,3)*f(4,1))/f(4,2)
        end if
        !    Set f(3,n) for not-a-knot
        IF(i_bcn.le.0.or.i_bcn.gt.7) THEN
           f(3,n)=f(3,n-1)+(f(3,n-1)-f(3,n-2))*f(4,n-1)/f(4,n-2)
        end if
     end if
     !f(3,i) is now the sigma(i) of the text and f(4,i) is the step size
     !Compute polynomial coefficients
     DO i=1,n-1
        f(2,i)=(f(1,i+1)-f(1,i))/f(4,i)-f(4,i)*(f(3,i+1)+2.0_fp*f(3,i))
        f(4,i)=(f(3,i+1)-f(3,i))/f(4,i)
        f(3,i)=6.0_fp*f(3,i)
        f(4,i)=6.0_fp*f(4,i)
     end do
     IF(i_bc1.eq.-1) THEN
        f(2,n)=f(2,1)
        f(3,n)=f(3,1)
        f(4,n)=f(4,1)
     ELSE
        hn=x(n)-x(n-1)
        f(2,n)=f(2,n-1)+hn*(f(3,n-1)+0.5_fp*hn*f(4,n-1))
        f(3,n)=f(3,n-1)+hn*f(4,n-1)
        f(4,n)=f(4,n-1)
        IF(i_bcn.eq.1.or.i_bcn.eq.3.or.i_bcn.eq.5) THEN
           f(2,n)=an
        ELSE IF(i_bcn.eq.2.or.i_bcn.eq.4.or.i_bcn.eq.6) THEN
           f(3,n)=bn
        end if
     end if
  end if
  RETURN
END subroutine v_spline

subroutine splinck(x,inx,ilinx,ztol,ier)
  !
  !  check if a grid is strictly ascending and if it is evenly spaced
  !  to w/in ztol
  !
  !============
  implicit none
  integer inx,ix
  !============
  real(fp) :: dxavg,zeps,zdiffx,zdiff
  !============
  real(fp) :: x(inx)                       ! input -- grid to check
  !
  integer ilinx                     ! output -- =1 if evenly spaced =2 O.W.
  !
  real(fp) :: ztol                         ! input -- spacing check tolerance
  !
  integer ier                       ! output -- =0 if OK
  !
  !  ier=1 is returned if x(1...inx) is NOT STRICTLY ASCENDING...
  !
  !-------------------------------
  !
  ier=0
  ilinx=1
  if(inx.le.1) return
  !
  dxavg=(x(inx)-x(1))/(inx-1)
  zeps=abs(ztol*dxavg)
  !
  do ix=2,inx
     zdiffx=(x(ix)-x(ix-1))
     if(zdiffx.le.0.0_fp) ier=2
     zdiff=zdiffx-dxavg
     if(abs(zdiff).gt.zeps) then
        ilinx=2
     end if
  end do
10 continue
  !
  return
end subroutine splinck

subroutine ibc_ck(ibc,slbl,xlbl,imin,imax,ier)
   implicit none
   ! Check that spline routine ibc flag is in range
   integer, intent(in)  :: ibc          ! flag value
   character(len=*), intent(in) :: slbl ! subroutine name
   character(len=*), intent(in) :: xlbl ! axis label
   integer, intent(in)  :: imin         ! min allowed value
   integer, intent(in)  :: imax         ! max allowed value
   integer, intent(out) :: ier          ! output -- set =1 if error detected

   if((ibc.lt.imin).or.(ibc.gt.imax)) then
      ier=1
 !      write(6,1001) slbl,xlbl,ibc,imin,imax
 !1001  format(' ?',a,' -- ibc',a,' = ',i9,' out of range ',i2,' to ',i2)
   end if

   return
end subroutine ibc_ck

subroutine EZspline_interp2_FOvars_cloud(spline_oBR, spline_oBPHI, &
      spline_oBZ, spline_oER, spline_oEPHI, spline_oEZ, p1, p2, fBR, &
      fBPHI, fBZ, fER, fEPHI, fEZ, ier)
   !$acc routine seq
   type(EZspline2) spline_oBR,spline_oBPHI,spline_oBZ
   type(EZspline2) spline_oER,spline_oEPHI,spline_oEZ
   real(fp), intent(in) :: p1, p2
   real(fp), intent(out):: fBR, fBPHI, fBZ
   real(fp), intent(out):: fER, fEPHI, fEZ
   integer, intent(out) :: ier
   integer :: ifail
   integer:: iwarn = 0

   !$acc routine (EZspline_allocated2) seq
   !$acc routine (evbicub_FOvars) seq

   ier = 0
   ifail = 0
  
   if( .not.EZspline_allocated2(spline_oBR) .or. spline_oBR%isReady /= 1) then
      ier = 94
      return
   endif

   call evbicub_FOvars(p1, p2,  &
      spline_oBR%x1(1), spline_oBR%n1, &
      spline_oBR%x2(1), spline_oBR%n2, &
      spline_oBR%ilin1, spline_oBR%ilin2, &
      spline_oBR%fspl(1,1,1), spline_oBPHI%fspl(1,1,1), &
      spline_oBZ%fspl(1,1,1), spline_oER%fspl(1,1,1), &
      spline_oEPHI%fspl(1,1,1), spline_oEZ%fspl(1,1,1), &
      spline_oBR%n1, &
      fBR, fBPHI, fBZ, fER, fEPHI, fEZ, ifail)

   if(ifail /= 0) ier = 97
  
end subroutine EZspline_interp2_FOvars_cloud
  
logical function EZspline_allocated2(spline_o)
   !$acc routine seq
   type(EZspline2) spline_o
   EZspline_allocated2 = allocated(spline_o%fspl) &
         .and. allocated(spline_o%x1) .and. allocated(spline_o%x1pkg) &
         .and. allocated(spline_o%x2) .and. allocated(spline_o%x2pkg) &
         .and. (spline_o%nguard == 123456789) ! check that ezspline_init has been called
end function EZspline_allocated2
  
subroutine evbicub_FOvars(xget,yget,x,nx,y,ny,ilinx,iliny, &
   fBR,fBPHI,fBZ,fER,fEPHI,fEZ,inf2,fvalBR,fvalBPHI,fvalBZ, &
   fvalER,fvalEPHI,fvalEZ,ier)
   !$acc routine seq
   !
   !  evaluate a 2d cubic Spline interpolant on a rectilinear
   !  grid -- this is C2 in both directions.
   !
   !  this subroutine calls two subroutines:
   !     herm2xy  -- find cell containing (xget,yget)
   !     fvbicub  -- evaluate interpolant function and (optionally) derivatives
   !
   !  input arguments:
   !  ================
   !
   !============
   implicit none
   integer inf2
   !============
   integer,intent(in) :: nx,ny                     ! grid sizes
   real(fp) :: xget,yget                    ! target of this interpolation
   real(fp) :: x(nx)                        ! ordered x grid
   real(fp) :: y(ny)                        ! ordered y grid
   integer ilinx                     ! ilinx=1 => assume x evenly spaced
   integer iliny                     ! iliny=1 => assume y evenly spaced
   !
   real(fp) :: fBR(0:3,inf2,ny)               ! function data
   real(fp) :: fBPHI(0:3,inf2,ny)
   real(fp) :: fBZ(0:3,inf2,ny)
   real(fp) :: fER(0:3,inf2,ny)
   real(fp) :: fEPHI(0:3,inf2,ny)
   real(fp) :: fEZ(0:3,inf2,ny)
   !
   !       f 2nd dimension inf2 must be .ge. nx
   !       contents of f:
   !
   !  f(0,i,j) = f @ x(i),y(j)
   !  f(1,i,j) = d2f/dx2 @ x(i),y(j)
   !  f(2,i,j) = d2f/dy2 @ x(i),y(j)
   !  f(3,i,j) = d4f/dx2dy2 @ x(i),y(j)
   !
   !      (these are spline coefficients selected for continuous 2-
   !      diffentiability, see mkbicub[w].f90)
   !
   !
   !  ict(1)=1 -- return f  (0, don't)
   !  ict(2)=1 -- return df/dx  (0, don't)
   !  ict(3)=1 -- return df/dy  (0, don't)
   !  ict(4)=1 -- return d2f/dx2  (0, don't)
   !  ict(5)=1 -- return d2f/dy2  (0, don't)
   !  ict(6)=1 -- return d2f/dxdy (0, don't)
   !                   the number of non zero values ict(1:6)
   !                   determines the number of outputs...
   !
   !  new dmc December 2005 -- access to higher derivatives (even if not
   !  continuous-- but can only go up to 3rd derivatives on any one coordinate.
   !     if ict(1)=3 -- want 3rd derivatives
   !          ict(2)=1 for d3f/dx3
   !          ict(3)=1 for d3f/dx2dy
   !          ict(4)=1 for d3f/dxdy2
   !          ict(5)=1 for d3f/dy3
   !               number of non-zero values ict(2:5) gives no. of outputs
   !     if ict(1)=4 -- want 4th derivatives
   !          ict(2)=1 for d4f/dx3dy
   !          ict(3)=1 for d4f/dx2dy2
   !          ict(4)=1 for d4f/dxdy3
   !               number of non-zero values ict(2:4) gives no. of outputs
   !     if ict(1)=5 -- want 5th derivatives
   !          ict(2)=1 for d5f/dx3dy2
   !          ict(3)=1 for d5f/dx2dy3
   !               number of non-zero values ict(2:3) gives no. of outputs
   !     if ict(1)=6 -- want 6th derivatives
   !          d6f/dx3dy3 -- one value is returned.
   !
   ! output arguments:
   ! =================
   !
   real(fp) :: fvalBR                      ! output data
   real(fp) :: fvalBPHI
   real(fp) :: fvalBZ
   real(fp) :: fvalER
   real(fp) :: fvalEPHI
   real(fp) :: fvalEZ

   integer ier                       ! error code =0 ==> no error
   !
   !  fval(1) receives the first output (depends on ict(...) spec)
   !  fval(2) receives the second output (depends on ict(...) spec)
   !  fval(3) receives the third output (depends on ict(...) spec)
   !  fval(4) receives the fourth output (depends on ict(...) spec)
   !  fval(5) receives the fourth output (depends on ict(...) spec)
   !  fval(6) receives the fourth output (depends on ict(...) spec)
   !
   !  examples:
   !    on input ict = [1,1,1,0,0,1]
   !   on output fval= [f,df/dx,df/dy,d2f/dxdy], elements 5 & 6 not referenced.
   !
   !    on input ict = [1,0,0,0,0,0]
   !   on output fval= [f] ... elements 2 -- 6 never referenced.
   !
   !    on input ict = [0,0,0,1,1,0]
   !   on output fval= [d2f/dx2,d2f/dy2] ... elements 3 -- 6 never referenced.
   !
   !    on input ict = [0,0,1,0,0,0]
   !   on output fval= [df/dy] ... elements 2 -- 6 never referenced.
   !
   !  ier -- completion code:  0 means OK
   !-------------------
   !  local:
   !
   integer i,j                      ! cell indices
   !
   !  normalized displacement from (x(i),y(j)) corner of cell.
   !    xparam=0 @x(i)  xparam=1 @x(i+1)
   !    yparam=0 @y(j)  yparam=1 @y(j+1)
   !
   real(fp) :: xparam,yparam
   !
   !  cell dimensions and
   !  inverse cell dimensions hxi = 1/(x(i+1)-x(i)), hyi = 1/(y(j+1)-y(j))
   !
   real(fp) :: hx,hy
   real(fp) :: hxi,hyi
   !
   !  0 .le. xparam .le. 1
   !  0 .le. yparam .le. 1
   !
   !  ** the interface is very similar to herm2ev.f90; can use herm2xy **
   !---------------------------------------------------------------------
   !$acc routine (herm2xy) seq
   !$acc routine (fvbicub) seq
   !
   i=0
   j=0
   call herm2xy(xget,yget,x,nx,y,ny,ilinx,iliny, &
      i,j,xparam,yparam,hx,hxi,hy,hyi,ier)
   if(ier.ne.0) return
   !
   call fvbicub(fvalBR,i,j,xparam,yparam,hx,hxi,hy,hyi,fBR,inf2,ny)
   call fvbicub(fvalBPHI,i,j,xparam,yparam,hx,hxi,hy,hyi,fBPHI,inf2,ny)
   call fvbicub(fvalBZ,i,j,xparam,yparam,hx,hxi,hy,hyi,fBZ,inf2,ny)
   call fvbicub(fvalER,i,j,xparam,yparam,hx,hxi,hy,hyi,fER,inf2,ny)
   call fvbicub(fvalEPHI,i,j,xparam,yparam,hx,hxi,hy,hyi,fEPHI,inf2,ny)
   call fvbicub(fvalEZ,i,j,xparam,yparam,hx,hxi,hy,hyi,fEZ,inf2,ny)
   !
   return
end subroutine evbicub_FOvars

subroutine herm2xy(xget,yget,x,nx,y,ny,ilinx,iliny, &
   i,j,xparam,yparam,hx,hxi,hy,hyi,ier)
   !---------------------------------------------------------------------
   !  herm2xy -- look up x-y zone
   !
   !  this is the "first part" of herm2ev, see comments, above.
   !
   !
   !  input of herm2xy
   !  ================
   !
   !============
   !$acc routine seq
   implicit none
   integer nxm,nym,ii,jj
   !============
   real(fp) :: zxget,zyget,zxtol,zytol
   !============
   integer nx,ny                     ! array dimensions
   !
   real(fp) :: xget,yget                    ! target point
   real(fp) :: x(nx),y(ny)                  ! indep. coords., strict ascending
   !
   integer ilinx                     ! =1:  x evenly spaced
   integer iliny                     ! =1:  y evenly spaced
   !
   !  output of herm2xy
   !  =================
   integer i,j                       ! index to cell containing target pt
   !          on exit:  1.le.i.le.nx-1   1.le.j.le.ny-1
   !
   !  normalized position w/in (i,j) cell (btw 0 and 1):
   !
   real(fp) :: xparam                       ! (xget-x(i))/(x(i+1)-x(i))
   real(fp) :: yparam                       ! (yget-y(j))/(y(j+1)-y(j))
   !
   !  grid spacing
   !
   real(fp) :: hx                           ! hx = x(i+1)-x(i)
   real(fp) :: hy                           ! hy = y(j+1)-y(j)
   !
   !  inverse grid spacing:
   !
   real(fp) :: hxi                          ! 1/hx = 1/(x(i+1)-x(i))
   real(fp) :: hyi                          ! 1/hy = 1/(y(j+1)-y(j))
   !
   integer ier                       ! return ier.ne.0 on error
   !
   !------------------------------------
   !$acc routine (zonfind) seq
   !
   ier=0
   !
   !  range check
   !
   zxget=xget
   zyget=yget
   if((xget.lt.x(1)).or.(xget.gt.x(nx))) then
      zxtol=4.0E-7_fp*max(abs(x(1)),abs(x(nx)))
      if((xget.lt.x(1)-zxtol).or.(xget.gt.x(nx)+zxtol)) then
         ier=1
 !        write(6,1001) xget,x(1),x(nx)
 !  1001    format(' ?herm2ev:  xget=',1pe11.4,' out of range ', &
 !           1pe11.4,' to ',1pe11.4)
      else
         if((xget.lt.x(1)-0.5_fp*zxtol).or. &
            (xget.gt.x(nx)+0.5_fp*zxtol)) &
 !           write(6,1011) xget,x(1),x(nx)
 !  1011    format(' %herm2ev:  xget=',1pe15.8,' beyond range ', &
 !           1pe15.8,' to ',1pe15.8,' (fixup applied)')
         if(xget.lt.x(1)) then
            zxget=x(1)
         else
            zxget=x(nx)
         end if
      end if
   end if
   if((yget.lt.y(1)).or.(yget.gt.y(ny))) then
      zytol=4.0E-7_fp*max(abs(y(1)),abs(y(ny)))
      if((yget.lt.y(1)-zytol).or.(yget.gt.y(ny)+zytol)) then
         ier=1
 !        write(6,1002) yget,y(1),y(ny)
 !  1002    format(' ?herm2ev:  yget=',1pe11.4,' out of range ', &
 !           1pe11.4,' to ',1pe11.4)
      else
         if((yget.lt.y(1)-0.5_fp*zytol).or. &
            (yget.gt.y(ny)+0.5_fp*zytol)) &
 !           write(6,1012) yget,y(1),y(ny)
 !  1012    format(' %herm2ev:  yget=',1pe15.8,' beyond range ', &
 !           1pe15.8,' to ',1pe15.8,' (fixup applied)')
         if(yget.lt.y(1)) then
            zyget=y(1)
         else
            zyget=y(ny)
         end if
      end if
   end if
   if(ier.ne.0) return
   !
   !  now find interval in which target point lies..
   !
   nxm=nx-1
   nym=ny-1
   !
   if(ilinx.eq.1) then
      ii=1+nxm*(zxget-x(1))/(x(nx)-x(1))
      i=min(nxm, ii)
      if(zxget.lt.x(i)) then
         i=i-1
      else if(zxget.gt.x(i+1)) then
         i=i+1
      end if
   else
      if((1.le.i).and.(i.lt.nxm)) then
         if((x(i).le.zxget).and.(zxget.le.x(i+1))) then
            continue  ! already have the zone
         else
            call zonfind(x,nx,zxget,i)
         end if
      else
         i=nx/2
         call zonfind(x,nx,zxget,i)
      end if
   end if
   !
   if(iliny.eq.1) then
      jj=1+nym*(zyget-y(1))/(y(ny)-y(1))
      j=min(nym, jj)
      if(zyget.lt.y(j)) then
         j=j-1
      else if(zyget.gt.y(j+1)) then
         j=j+1
      end if
   else
      if((1.le.j).and.(j.lt.nym)) then
         if((y(j).le.zyget).and.(zyget.le.y(j+1))) then
            continue  ! already have the zone
         else
            call zonfind(y,ny,zyget,j)
         end if
      else
         j=ny/2
         call zonfind(y,ny,zyget,j)
      end if
   end if
   !
   hx=(x(i+1)-x(i))
   hy=(y(j+1)-y(j))
   !
   hxi=1.0_fp/hx
   hyi=1.0_fp/hy
   !
   xparam=(zxget-x(i))*hxi
   yparam=(zyget-y(j))*hyi
   !
   return
end subroutine herm2xy

subroutine zonfind(x,nx,zxget,i)
   !$acc routine seq
   implicit none
   integer nx,nxm,i1,i2,ij,ii
   real(fp) :: dx
   real(fp) :: x(nx),zxget
   integer i
   !
   !  find index i such that x(i).le.zxget.le.x(i+1)
   !
   !  x(1...nx) is strict increasing and x(1).le.zxget.le.x(nx)
   !  (this is assumed to already have been checked -- no check here!)
   !
   nxm=nx-1
   if((i.lt.1).or.(i.gt.nxm)) then
      i1=1
      i2=nx-1
      go to 10
   end if
 
   if(x(i).gt.zxget) then
      ! look down
      dx=x(i+1)-x(i)
      if((x(i)-zxget).gt.4*dx) then
         i1=1
         i2=i-1
         go to 10
      else
         i2=i-1
         do ij=i2,1,-1
            if((x(ij).le.zxget).and.(zxget.le.x(ij+1))) then
               i=ij
               return
            end if
         end do
         i=1
         return
      end if
   else if(x(i+1).lt.zxget) then
      ! look up
      dx=x(i+1)-x(i)
      if((zxget-x(i+1)).gt.4*dx) then
         i1=i+1
         i2=nxm
         go to 10
      else
         i2=i+1
         do ij=i2,nxm
            if((x(ij).le.zxget).and.(zxget.le.x(ij+1))) then
               i=ij
               return
            end if
         end do
         ij=nxm
         return
      end if
   else
      ! already there...
      return
   end if
 
   !---------------------------
   !  binary search
   !
 10 continue
 
   if(i1.eq.i2) then
      ! found by proc. of elimination
      i=i1
      return
   end if
 
   ii=(i1+i2)/2
 
   if(zxget.lt.x(ii)) then
      i2=ii-1
   else if(zxget.gt.x(ii+1)) then
      i1=ii+1
   else
      ! found
      i=ii
      return
   end if
 
   go to 10
 
   return
end subroutine zonfind
  
subroutine fvbicub(fval,i,j,xparam,yparam,hx,hxi,hy,hyi, &
   fin,inf2,ny)
   !$acc routine seq
   !
   !============
   implicit none
   integer ny,inf2
   !============
   real(fp) :: z36th,xp,xpi,xp2,xpi2,cx,cxi,hx2,yp,ypi,yp2,ypi2,cy
   real(fp) :: cyi,hy2,cxd,cxdi,cyd,cydi
   !============
   !
   integer i,j         ! target cells (i,j)
   real(fp) :: xparam,yparam
   ! normalized displacements from (i,j) corners
   !
   real(fp) :: hx,hy           ! grid spacing, and
   real(fp) :: hxi,hyi         ! inverse grid spacing 1/(x(i+1)-x(i))
   ! & 1/(y(j+1)-y(j))
   !
   real(fp) :: fin(0:3,inf2,ny)             ! interpolant data (cf "evbicub")
   !
   real(fp) :: fval                ! output returned
   !
   !  for detailed description of fin, ict and fval see subroutine
   !  evbicub comments.  Note ict is not vectorized; the same output
   !  is expected to be returned for all input vector data points.
   !
   !  note that the index inputs ii,jj and parameter inputs
   !     xparam,yparam,hx,hxi,hy,hyi are vectorized, and the
   !     output array fval has a vector ** 1st dimension ** whose
   !     size must be given as a separate argument
   !
   !  to use this routine in scalar mode, pass in ivec=ivecd=1
   !
   !---------------
   !  Spline evaluation consists of a "mixing" of the interpolant
   !  data using the linear functionals xparam, xpi = 1-xparam,
   !  yparam, ypi = 1-yparam, and the cubic functionals
   !  xparam**3-xparam, xpi**3-xpi, yparam**3-yparam, ypi**3-ypi ...
   !  and their derivatives as needed.
   !
   real(fp) :: sum
   !
   real(fp), parameter :: sixth = 0.166666666666666667_fp
   !
   !---------------
   !
   z36th=sixth*sixth
      !
   !  in x direction...
   !
   xp=xparam
   xpi=1.0_fp-xp
   xp2=xp*xp
   xpi2=xpi*xpi
   !
   cx=xp*(xp2-1.0_fp)
   cxi=xpi*(xpi2-1.0_fp)
   hx2=hx*hx
   !
   !   ...and in y direction
   !
   yp=yparam
   ypi=1.0_fp-yp
   yp2=yp*yp
   ypi2=ypi*ypi
   !
   cy=yp*(yp2-1.0_fp)
   cyi=ypi*(ypi2-1.0_fp)
   hy2=hy*hy
   !
   sum=xpi*(ypi*fin(0,i,j)  +yp*fin(0,i,j+1))+ &
         xp*(ypi*fin(0,i+1,j)+yp*fin(0,i+1,j+1))
   !
   sum=sum+sixth*hx2*( &
         cxi*(ypi*fin(1,i,j)  +yp*fin(1,i,j+1))+ &
         cx*(ypi*fin(1,i+1,j)+yp*fin(1,i+1,j+1)))
   !
   sum=sum+sixth*hy2*( &
         xpi*(cyi*fin(2,i,j)  +cy*fin(2,i,j+1))+ &
         xp*(cyi*fin(2,i+1,j)+cy*fin(2,i+1,j+1)))
   !
   sum=sum+z36th*hx2*hy2*( &
         cxi*(cyi*fin(3,i,j)  +cy*fin(3,i,j+1))+ &
         cx*(cyi*fin(3,i+1,j)+cy*fin(3,i+1,j+1)))
   !
   fval=sum
   !
   return
end subroutine fvbicub
  
subroutine EZspline_free2(spline_o, ier)
   implicit none
   type(EZspline2) spline_o
   ! ier:
   ! 101= warning, spline object was never allocated
   integer, intent(out) :: ier
   integer ifail
   
   ier = 0
   if(.not.EZspline_allocated2(spline_o)) ier=101
   
   deallocate(spline_o%x1, stat=ifail)
   deallocate(spline_o%x2, stat=ifail)
   deallocate(spline_o%fspl, stat=ifail)
   deallocate(spline_o%bcval1min, stat=ifail)
   deallocate(spline_o%bcval1max, stat=ifail)
   deallocate(spline_o%bcval2min, stat=ifail)
   deallocate(spline_o%bcval2max, stat=ifail)
   deallocate(spline_o%x1pkg, stat=ifail)
   deallocate(spline_o%x2pkg, stat=ifail)
   
   call EZspline_preInit2(spline_o)
   
   spline_o%n1 = 0
   spline_o%n2 = 0
   
   return
end subroutine EZspline_free2
  
subroutine EZspline_error(ier)
   !$acc routine seq
   !
   ! Error handling routine. Maps error ier code to a meaningful message.
   ! Note: does not abort nor stop if ier/=0.
   !
   implicit none
   integer, intent(in) :: ier

   if(ier == 0) return
   write(6,*) '**EZspline** ERROR/WARNING #', ier,' occurred'

   select case(ier)
   case(1)
      write(6,*) '**EZspline** allocation error'
   case(2)
      write(6,*) '**EZspline** wrong BCS1 code'
   case(3)
      write(6,*) '**EZspline** wrong BCS2 code'
   case(4)
      write(6,*) '**EZspline** wrong BCS3 code'
   case(5)
      write(6,*) '**EZspline** Que??'
   case(6)
      write(6,*) '**EZspline** out of interval p1 < min(x1)'
   case(7)
      write(6,*) '**EZspline** out of interval p1 > max(x1)'
   case(8)
      write(6,*) '**EZspline** out of interval p2 < min(x2)'
   case(9)
      write(6,*) '**EZspline** out of interval p2 > max(x2)'
   case(10)
      write(6,*) '**EZspline** out of interval p3 < min(x3)'
   case(11)
      write(6,*) '**EZspline** out of interval p3 > max(x3)'
   case(12)
      write(6,*) '**EZspline** negative derivative order'
   case(13)
      write(6,*) '**EZspline** derivative order too high'
   case(14)
      write(6,*) '**EZspline** x1 grid is not strictly increasing'
   case(15)
      write(6,*) '**EZspline** x2 grid is not strictly increasing'
   case(16)
      write(6,*) '**EZspline** x3 grid is not strictly increasing'
   case(17)
      write(6,*) '**EZspline** could not save spline object in file '
   case(18)
      write(6,*) '**EZspline** memory allocation failure in coefficient setup'

   case(20)
      write(6,*) '**EZspline** attempt to load spline object with wrong rank.'
   case(21)
      write(6,*) '**EZspline** could not load spline object from file '
   case(22)
      write(6,*) '**EZspline** loaded spline object from file but failed at coefficient set-up'
   case(23)
      write(6,*) '**EZspline** failed to free spline object'
   case(24)
      write(6,*) '**EZspline** 2nd order derivative not supported for Akima-Hermite (isHermite=1)'
   case(25)
      write(6,*) '**EZspline** not supported for Akima-Hermite (isHermite=1)'
   case(26)
      write(6,*) '**EZspline** memory allocation error in EZspline_interp'
   case(27)
      write(6,*) '**EZspline** an error ocurred in genxpkg'
   case(28)
      write(6,*) '**EZspline** memory allocation failure in ezspline_interp'
   case(29)
      write(6,*) '**EZspline** memory deallocation failure in ezspline_interp'
   case(30)
      write(6,*) '**EZspline** memory allocation error in EZspline_gradient'
   case(31)
      write(6,*) '**EZspline** memory deallocation error in EZspline_gradient'
   case(32)
      write(6,*) '**EZspline** memory allocation error in EZspline_derivative'
   case(33)
      write(6,*) '**EZspline** memory deallocation error in EZspline_derivative'
   case(34)
      write(6,*) '**EZspline** could not open netCDF file in EZspline_2netcdf'
   case(35)
      write(6,*) '**EZspline** could not write into netCDF file in EZspline_2netcdf'
   case(36)
      write(6,*) '**EZspline** could not read from netCDF file in EZspline_2netcdf'
   case(37)
      write(6,*) '**EZspline** could not close netCDF file in EZspline_2netcdf'
   case(38)
      write(6,*) '**EZspline** could not define variable (cdfDefVar) in EZspline_2netcdf'
   case(39)
      write(6,*) '**EZspline** could not open netCDF file in EZspline_save'
   case(40)
      write(6,*) '**EZspline** could not write into netCDF file in EZspline_save'
   case(41)
      write(6,*) '**EZspline** could not close netCDF file in EZspline_save'
   case(42)
      write(6,*) '**EZspline** could not define variable (cdfDefVar) in EZspline_save'
   case(43)
      write(6,*) '**EZspline** could not open netCDF file in EZspline_load'
   case(44)
      write(6,*) '**EZspline** could not read from netCDF file in EZspline_load'
   case(45)
      write(6,*) '**EZspline** could not close netCDF file in EZspline_load'
   case(46)
      write(6,*) '**EZspline** 2nd order derivative not supported for Piecewise Linear Interpolation (isLinear=1)'
   case(47)
      write(6,*) '**EZspline** not supported for Piecewise Linear Interpolation (isLinear=1)'
   case(50)
      write(6,*) '**EZspline** ezspline_save (optional) spline name is blank.'
   case(51)
      write(6,*) '**EZspline** ezspline_save (optional) spline name too long (max 20 characters).'
   case(52)
      write(6,*) '**EZspline** ezspline_save (optional) spline name contains'
      write(6,*) '             imbedded blanks or other illegal characters.'
   case(53)
      write(6,*) '**EZspline** attempt to write named spline object to NetCDF'
      write(6,*) '             file with change of dimensionality or data type.'
   case(54)
      write(6,*) '**EZspline** hybrid interpolation specification not in range -1:2'
      write(6,*) '             error in EZhybrid_init.'
   case(55)
      write(6,*) '**EZspline** hybrid interpolation cannot mix Hermite and Spline interpolation.'
      write(6,*) '             hspline(i)=1 and hspline(j)=2 in EZhybrid_init.'
   case(56)
      write(6,*) '**EZspline** non-default boundary condition unavailable: zonal or piecewise linear dimension.'
      write(6,*) '             in EZhybrid_init.'
   case(57)
      write(6,*) '**EZspline** dimension of "f" smaller than corresponding "fspl"'
      write(6,*) '             dimension in "spline_o".'
   case(58)
      write(6,*) '**EZspline** dimension of "f" larger than corresponding "fspl"'
      write(6,*) '             dimension in "spline_o".'
   case(90)
      write(6,*) '**EZspline** an error occurred after attempting to evaluate the'
      write(6,*) '             Hermite polynomials'
   case(91)
      write(6,*) '**EZspline** an error occurred after attempting to set up the'
      write(6,*) '             Hermite polynomial coefficients'
   case(92)
      write(6,*) '**EZspline** warning in EZspline_load. Looks like saved object '
      write(6,*) '             was not properly set-up (isReady=0).'
   case(93)
      write(6,*) '**EZspline** warning in EZspline_save. Looks like saved object '
      write(6,*) '             was not properly set-up (isReady=0).'
   case(94)
      write(6,*) '**EZspline** an error occurred in EZspline_interp. Did you forget'
      write(6,*) '             to set up the cubic spline coefficients by calling'
      write(6,*) '             call EZspline_setup(spline_o, f, ier)'
      write(6,*) '             ?'
   case(95)
      write(6,*) '**EZspline** some error occurred in EZspline_gradient'
   case(96)
      write(6,*) '**EZspline** some error occurred in EZspline_derivative'
   case(97)
      write(6,*) '**EZspline** some error occurred in EZspline_interp apparently'
      write(6,*) '             related to a PSPLINE routine. Check if argument is '
      write(6,*) '             outside interpolation domain by calling'
      write(6,*) '             call EZspline_isInDomain(spline_o, [[k1, k2, k3,] .OR. k,] p1, p2, p3, ier ,ier)'
      write(6,*) '             call EZspline_error(ier)'
   case(98)
      write(6,*) '**EZspline** error occurred in EZspline_setup'
      write(6,*) '  if no other explanation-- ezspline_init call never made.'
   case(99)
      write(6,*) '**EZspline** some error occurred in EZspline_init, EZhybrid_init,  or EZlinear_init'
   case(100)
      write(6,*) '**EZSPLINE** EZspline_init, EZhybrid_init,  or EZlinear_init -- object already allocated.'
   case(101)
      write(6,*) '**EZSPLINE** object was never allocated.'
   case default
      write(6,*) '**EZspline** '
   end select

   return
end subroutine EZspline_error
  
#endif PSPLINE
    
end module pspline_gpu