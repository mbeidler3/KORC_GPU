!-----------------------------------------------------------------------
!     $Id: fft2d.F90 4620 2015-12-03 19:18:40Z jking $
!-----------------------------------------------------------------------
!     code organization.                                                
!-----------------------------------------------------------------------
!     0. fft2d_mod.                                                       
!     1. fft2d.                                                         
!     2. cpftv.       
!     3. nim_fftw
!-----------------------------------------------------------------------
#include "config.f"
      MODULE fft2d_mod
#ifdef HAVE_OPENMP
      USE omp_lib
#endif 
      IMPLICIT NONE
      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 1. fft2d                                               
!                                                                       
!     slices data into sections for cray vector machines.  an option    
!     for non-vector machines may be desirable.                         
!                                                                       
!     ###### CRAY-2 de-aliased version. ######                          
!                                                                       
!     Performs the FFT of a two-dimensional array over                  
!     the second dimension, for all points in the first                 
!     dimension.  The operations over this first dimension are          
!     vectorized.  This routine uses Oscar Buneman's program            
!     VCFT.  Load the binary VCFT on the Cray-2.                        
!     (This binary is available from FILEM 1505 .BUTILITY)              
!     WARNING: The CAL VCFT on the Cray-2 requires that MY be           
!     between 4 and 21.                                                 
!                                                                       
!     compr and compi have been combined into one array with both       
!     real and imaginary parts to reduce data rearranging in nimrod.    
!     (CRS, 8/14/98)                                                    
!                                                                       
!     The work array is placed here, so that the loop-splitting         
!     can be easily varied by changing the incf paramter.  The          
!     indices of the 'imaginary' part of the data have been changed,    
!     and the work array is now complex for optimization.               
!     (CRS, 12/10/98)                                                   
!                                                                       
!     This version of fft2d uses a complex comp array directly.         
!     (CRS, 1/13/99)                                                    
!                                                                       
!     This program can be compiled with CIVIC, CFT2, or CFT77.          
!                                                                       
!     The definition of the arguments is as follows:                    
!                                                                       
!       IFLAG:   An integer flag which determines the direction of      
!                the FFT.  Fourier analysis is performed when IFLAG=1,  
!                and Fourier synthesis (i.e., back to real space) is    
!                performed when IFLAG.ne.1.                             
!                                                                       
!       NX:      The number of points in the non-FFT direction          
!                (over which the vectorization is performed).           
!                                                                       
!       MY:      The power-of-two which gives the number of points      
!                in the FFT direction.  Thus the number of points in    
!                the FFT direction is NY=2**MY.  Note that MY needs     
!                to be between 2 and 10, inclusive.                     
!                                                                       
!       RE:      Real array dimensioned NX by NY which contains the     
!                data in real space. This array is used as input        
!                when IFLAG=1, and is output when IFLAG.ne.1.           
!                                                                       
!       COMP:    Real and imaginary part of complex array dimensioned   
!                NX by 2*(NY/3+1) which contains the non-aliased        
!                Fourier modes.  This array is used as input when       
!                IFLAG.ne.1, and is output when IFLAG=1.                
!                                                                       
!       DEALIASE: When this optional input is present and set to F, the 
!                 routine uses all possible Fourier coefficients.  See  
!                 below.                                                
!                                                                       
!     Note that RE is not overwritten in the call with IFLAG=1,         
!     and similarly, COMP is not overwritten in the call with           
!     IFLAG.ne.1.                                                       
!                                                                       
!     The operation of this routine can be summarized by:               
!                                                                       
!                           NY                                          
!       COMP(I,M) = 1/NY * SUM RE(I,J)*EXP[-2*pi*i*(J-1)*(M-1)/NY]      
!                          J=1                                          
!                                                                       
!                         for I=1,2,...,NX and M=1,2,...,NY/2+1         
!                         when IFLAG=1, and                             
!                                                                       
!                  NY                                                   
!       RE(I,J) = SUM COMP(I,M)*EXP[2*pi*i*(J-1)*(M-1)/NY]              
!                 M=1                                                   
!                                                                       
!                         for I=1,2,...,NX and J=1,2,...,NY             
!                         when IFLAG.ne.1, where the elements COMP(I,M) 
!                         for M=NY/2+2,...,NY are not stored in array   
!                         COMP, but obey COMP(I,M)=conj(COMP(I,NY-M+2)).
!                                                                       
!                         Standard operation for this routine de-aliases
!                         the COMP array such that only the             
!                         modes with M=1,2,...,NY/3+1 are nonzero.      
!                         However, if the optional dealiase flag is     
!                         present and set to false, all (1<=M<=NY/2+1)  
!                         are used.                                     
!-----------------------------------------------------------------------
      SUBROUTINE fft2d(re,comp,nx,my,iflag,dealiase) 
      USE local 
      IMPLICIT NONE 
                                                                        
      INTEGER(i4), INTENT(IN) :: nx,my,iflag 
      REAL(r8), DIMENSION(nx,*), INTENT(INOUT) :: re 
      COMPLEX(r8), DIMENSION(nx,*), INTENT(INOUT) :: comp 
      LOGICAL, INTENT(IN), OPTIONAL :: dealiase 
                                                                        
      INTEGER(i4) :: ny,nyd3,n,n2,n2p,i,j,m,k,mm,i1,i2 
      INTEGER(i4), PARAMETER :: incf=24,istride=2*incf-2 
      COMPLEX(r8), DIMENSION(incf,2**my+1) :: work 
      COMPLEX(r8) :: ctmp 
      REAL(r8) :: xnrm,xnrm2 
!-----------------------------------------------------------------------
!                                                                       
      ny=2**my 
      nyd3=ny/3+1 
      if (present(dealiase)) then 
        if (.not.dealiase) nyd3=ny/2+1 !  800 loop is skipped entirely 
      endif 
      xnrm=1._r8/ny 
      xnrm2=xnrm*0.5 
 
      if (iflag==1) then 
 
!     FORWARD TRANSFORM.                                                
 
      do 400 i=1,nx,istride 
      n=MIN(istride,nx-i+1_i4) 
      n2=n/2 
      n2p=n-n2 
 
!     Pack two x-lines at a time.                                       
 
      do 100 j=1,ny 
      do 100 k=1,n2 
      i1=i-2+2*k 
      i2=i1+1 
  100 work(k,j)=re(i1,j)+(0,1)*re(i2,j) 
 
      if (n2<n2p) then 
        do 110 j=1,ny 
  110   work(n2p,j)=re(i-1+n,j) 
      endif 
 
      call cpftv(work,ny,incf,n2p,1_i4,-1_i4) 
 
!     Unpack the two x-lines (normalization moved here).                
 
!dir$ ivdep                                                             
      do 200 k=1,n2 
      i1=i-2+2*k 
      i2=i1+1 
      ctmp=xnrm*work(k,1) 
      comp(i1,1)=REAL(ctmp,r8) 
  200 comp(i2,1)=AIMAG(ctmp) 
 
      do 300 m=2,nyd3 
      mm=ny+2-m 
!dir$ ivdep                                                             
      do 300 k=1,n2 
      i1=i-2+2*k 
      i2=i1+1 
      ctmp=CONJG(work(k,mm)) 
      comp(i1,m)=xnrm2*(ctmp+work(k,m)) 
  300 comp(i2,m)=xnrm2*(ctmp-work(k,m))*(0,1) 
 
      if (n2<n2p) then 
        i1=i-1+n 
        comp(i1,1)=xnrm*REAL(work(n2p,1),r8) 
        do 310 m=2,nyd3 
        mm=ny+2-m 
  310   comp(i1,m)=xnrm2*(CONJG(work(n2p,mm))+work(n2p,m)) 
      endif 
 
  400 continue 
 
      else 
 
!     INVERSE TRANSFORM.                                                
 
      do 950 i=1,nx,istride 
      n=MIN(istride,nx-i+1_i4) 
      n2=n/2 
      n2p=n-n2 
 
!     Pack two x-lines at atime.                                        
 
!     m=0                                                               
 
      do 500 k=1,n2 
      i1=i-2+2*k 
      i2=i1+1 
  500 work(k,1)=REAL(comp(i1,1),r8)+(0,1)*REAL(comp(i2,1),r8) 
 
!     Positive and negative m.                                          
 
      do 600 m=2,nyd3 
      mm=ny+2-m 
!dir$ ivdep                                                             
      do 600 k=1,n2 
      i1=i-2+2*k 
      i2=i1+1 
      work(k,m )=comp(i1,m)+(0,1)*comp(i2,m) 
  600 work(k,mm)=CONJG(comp(i1,m))+(0,1)*CONJG(comp(i2,m)) 
 
      if (n2<n2p) then 
        i1=i-1+n 
        work(n2p,1)=REAL(comp(i1,1),r8) 
        do 610 m=2,nyd3 
        mm=ny+2-m 
        work(n2p,m )=comp(i1,m) 
  610   work(n2p,mm)=CONJG(comp(i1,m)) 
      endif 
 
!     De-aliased m.                                                     
 
      do 800 m=nyd3+1,ny-nyd3+1 
      do 800 k=1,n2p 
  800 work(k,m)=0. 
 
      call cpftv(work,ny,incf,n2p,1_i4,1_i4) 
 
      do 900 j=1,ny 
!dir$ ivdep                                                             
      do 900 k=1,n2 
      i1=i-2+2*k 
      i2=i1+1 
      re(i1,j)=work(k,j) 
  900 re(i2,j)=AIMAG(work(k,j)) 
 
      if (n2<n2p) then 
        do 910 j=1,ny 
        i1=i-1+n 
  910   re(i1,j)=work(n2p,j) 
      endif 
 
  950 continue 
 
      end if 
 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE fft2d 
!-----------------------------------------------------------------------
!     subprogram 2. cpftv                                               
!     perform transforms.                                               
!                                                                       
!     This subroutine performs an FFT in one direction of a             
!     two-dimensional array, for each row in the non-FFT direction.     
!     The loops are vectorized in this perpendicular direction.         
!     This routine is a modified version of Langdon's CPFT              
!     FFT routine (which is based on Singleton's FFT algorithm).        
!     This routine is unnormalized, in the sense that if it is          
!     called twice in succession (with the same array), first with      
!     ISIGN=1 and then with ISIGN=-1, the array is multiplied by N.     
!                                                                       
!     The meaning of the arguments is as follows:                       
!     CA    = complex array containing the the complex                  
!             data to be Fourier-transformed,                           
!     N     = length of the FFT (must be a power of 2),                 
!     INCF  = increment between successive elements of R (and I)        
!             in the FFT direction when it is regarded                  
!             as a real one-dimensional array,                          
!     NV    = number of rows over which to perform the FFT              
!             (i.e., in the vector direction),                          
!     INCV  = increment between successive elements of R (and I)        
!             in the perpendicular direction when it is regarded        
!             as a real one-dimensional array,                          
!     ISIGN = sign to be used in the exponential of the FFT             
!             (use +1 or -1).                                           
!                                                                       
!#### Modified by Z. Mikic, SAIC, La Jolla, August 18, 1985.            
!                                                                       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE cpftv (ca,n,incf,nv,incv,isign) 
      USE local 
      implicit none
                                                                        
      COMPLEX(r8), DIMENSION(:,:), INTENT(INOUT) :: ca 
      INTEGER(i4), INTENT(IN) :: n,incf,nv,incv,isign 
                                                                        
      INTEGER(i4) :: ji,ij,inc,is,ninc,n1,n2,nvl,span,rc,it,iv,k0,k1 
      INTEGER(i4), PARAMETER :: log2nx=15 
      REAL(r8), DIMENSION(log2nx), SAVE :: sines=0 
      REAL(r8) :: i0,i1,qt,qq,cosx,sinx,c,s,t,sgn 
      COMPLEX(r8) :: c0,c1,ct 
!-----------------------------------------------------------------------
!                                                                       
      if(sines(1)==1.) go to 1 
      sines(1)=1. 
      qt=1. 
      qt=atan(qt) 
      do 2 is=2,log2nx 
      qq=sin(qt) 
      sines(is)=qq 
    2 qt=qt*.5 
    1 continue 
 
      if(n==1) return 
      nvl=1+(nv-1)*incv 
!*ZM* For 2-d arrays R and I, set INC=1.                                
      inc=1 
      sgn=isign 
      ninc=n*inc 
      span=ninc 
      it=n/2 
      do 1000 is=1,log2nx 
!... (2000=recur)                                                       
      if(it==1) go to 2000 
 1000 it=it/2 
 
!  if truncated rather than rounded arithmetic is used,                 
!  singleton's magnitude correcton should be applied to cos and sin.    
 1500 t=sinx+(s*cosx-c*sinx) 
      cosx=cosx-(c*cosx+s*sinx) 
      sinx=t 
!... (3000=repl)                                                        
 3000 k1=k0+span 
!*ZM* Vector loops introduced all end with 01.                          
!dir$ ivdep                                                             
      do 101 iv=1,nvl,incv 
      c0=ca(iv,1+k0) 
      c1=ca(iv,1+k1) 
      ca(iv,1+k0)=c0+c1 
      c0=c0-c1 
      ca(iv,1+k1)=(cosx+(0,1)*sinx)*c0 
  101 continue 
      k0=k1+span 
      if(k0<ninc) go to 3000 
      k1=k0-ninc 
      cosx=-cosx 
      k0=span-k1 
      if(k1<k0) go to 3000 
      k0=k0+inc 
      k1=span-k0 
      if(k0<k1) go to 1500 
 2000 continue 
      span=span/2 
      k0=0 
!... (4000=zero)                                                        
 4000 k1=k0+span 
!dir$ ivdep                                                             
      do 201 iv=1,nvl,incv 
      c0=ca(iv,1+k0) 
      c1=ca(iv,1+k1) 
      ca(iv,1+k0)=c0+c1 
      ca(iv,1+k1)=c0-c1 
  201 continue 
      k0=k1+span 
      if(k0<ninc) go to 4000 
      if(span==inc) go to 5000 
      k0=span/2 
 4500 k1=k0+span 
!dir$ ivdep                                                             
      do 301 iv=1,nvl,incv 
      c0=ca(iv,1+k0) 
      c1=ca(iv,1+k1) 
      ca(iv,1+k0)=c0+c1 
      ca(iv,1+k1)=sgn*(0,1)*(c0-c1) 
  301 continue 
      k0=k1+span 
      if(k0<ninc) go to 4500 
      k1=inc+inc 
      if(span==k1) go to 2000 
      c=2.*sines(is)**2 
      is=is-1 
      sinx=sign(sines(is),sgn) 
      s=sinx 
      cosx=1.-c 
      k0=inc 
      go to 3000 
 
 5000 n1=ninc-inc 
      n2=ninc/2 
      ij=0 
      ji=0 
      rc=0 
      if(n2==inc) return 
      go to 5020 
!... (5010=even)                                                        
 5010 ij=n1-ij 
      ji=n1-ji 
!dir$ ivdep                                                             
      do 401 iv=1,nvl,incv 
      ct=ca(iv,1+ij) 
      ca(iv,1+ij)=ca(iv,1+ji) 
      ca(iv,1+ji)=ct 
  401 continue 
      if(ij>n2) go to 5010 
!... (5020=odd)                                                         
 5020 ij=ij+inc 
      ji=ji+n2 
!dir$ ivdep                                                             
      do 501 iv=1,nvl,incv 
      ct=ca(iv,1+ij) 
      ca(iv,1+ij)=ca(iv,1+ji) 
      ca(iv,1+ji)=ct 
  501 continue 
      it=n2 
!... (6000=incrv)                                                       
 6000 it=it/2 
      rc=rc-it 
      if(rc>=0) go to 6000 
      rc=rc+2*it 
      ji=rc 
      ij=ij+inc 
      if(ij<=ji) go to 5010 
      if(ij<n2) go to 5020 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE cpftv 
#ifdef HAVE_FFTW3
!-----------------------------------------------------------------------
!     subprogram 3. nim_fftw
!                                                                       
!     The definition of the arguments is as follows:                    
!                                                                       
!       IFLAG:   An integer flag which determines the direction of      
!                the FFT.  Fourier analysis is performed when IFLAG=1,  
!                and Fourier synthesis (i.e., back to real space) is    
!                performed when IFLAG.ne.1.                             
!                                                                       
!       NX:      The number of points in the non-FFT direction          
!                (over which the vectorization is performed).           
!                                                                       
!       NY:      The number of points in the FFT direction.
!                                                                       
!       RE:      Real array dimensioned NX by NY which contains the     
!                data in real space. This array is used as input        
!                when IFLAG=1, and is output when IFLAG.ne.1.           
!                                                                       
!       COMP:    Real and imaginary part of complex array dimensioned   
!                NX by 2*(NY/3+1) which contains the non-aliased        
!                Fourier modes.  This array is used as input when       
!                IFLAG.ne.1, and is output when IFLAG=1.                
!                                                                       
!       DEALIASE: When this optional input is present and set to F, the 
!                 routine uses all possible Fourier coefficients.  See  
!                 below.                                                
!                                                                       
!     Note that RE is not overwritten in the call with IFLAG=1,         
!     and similarly, COMP is not overwritten in the call with           
!     IFLAG.ne.1.                                                       
!                                                                       
!     The operation of this routine can be summarized by:               
!                                                                       
!                           NY                                          
!       COMP(I,M) = 1/NY * SUM RE(I,J)*EXP[-2*pi*i*(J-1)*(M-1)/NY]      
!                          J=1                                          
!                                                                       
!                         for I=1,2,...,NX and M=1,2,...,NY/2+1         
!                         when IFLAG=1, and                             
!                                                                       
!                  NY                                                   
!       RE(I,J) = SUM COMP(I,M)*EXP[2*pi*i*(J-1)*(M-1)/NY]              
!                 M=1                                                   
!                                                                       
!                         for I=1,2,...,NX and J=1,2,...,NY             
!                         when IFLAG.ne.1, where the elements COMP(I,M) 
!                         for M=NY/2+2,...,NY are not stored in array   
!                         COMP, but obey COMP(I,M)=conj(COMP(I,NY-M+2)).
!                                                                       
!                         Standard operation for this routine de-aliases
!                         the COMP array such that only the             
!                         modes with M=1,2,...,NY/3+1 are nonzero.      
!                         However, if the optional dealiase flag is     
!                         present and set to false, all (1<=M<=NY/2+1)  
!                         are used.                                     
!-----------------------------------------------------------------------
      SUBROUTINE nim_fftw(re,comp,nx,ny,iflag,indealiase) 
      USE pardata
      IMPLICIT NONE 
                                                                        
      INTEGER(i4), INTENT(IN) :: nx,ny,iflag 
      REAL(r8), DIMENSION(nx,ny), INTENT(INOUT) :: re 
      COMPLEX(r8), DIMENSION(nx,*), INTENT(INOUT) :: comp 
      INTEGER(i4), INTENT(IN), OPTIONAL :: indealiase 

#ifdef __xlf
      INTEGER(i4), PARAMETER :: vecstride=4_i4 ! vectorization stride
#else
      INTEGER(i4), PARAMETER :: vecstride=16_i4 ! vectorization stride
#endif
      INTEGER(i4), SAVE :: fftw_ny
      LOGICAL, SAVE :: first_call=.TRUE.

      INTEGER(i4) :: it,nt,mt,nmodes,nm_tot,rdim(1),cdim(1)
      INTEGER(i4) :: iv,ii,ichunk,iend,vs
      TYPE(C_PTR) :: ptr
      INTEGER(i4) :: dealiase

!-----------------------------------------------------------------------
!     determine OpenMP threading information.
!-----------------------------------------------------------------------
#ifdef HAVE_OPENMP
      it=OMP_GET_THREAD_NUM()+1
      nt=OMP_GET_NUM_THREADS()
      mt=OMP_GET_MAX_THREADS()
#else
      it=1
      nt=1
      mt=1
#endif
!-----------------------------------------------------------------------
!     determine mode information from ny.
!-----------------------------------------------------------------------
      dealiase=3 ! default
      IF (PRESENT(indealiase)) dealiase=indealiase
      IF (dealiase<3) THEN
        nmodes=ny/2+1
      ELSE
        nmodes=ny/dealiase+1
      ENDIF
      nm_tot=ny/2+1
!-----------------------------------------------------------------------
!     allocate if needed
!-----------------------------------------------------------------------
      IF (first_call) THEN
        IF (it==1) THEN 
          ALLOCATE(thread_alloc(mt))
          ALLOCATE(plan_r2c(mt),plan_c2r(mt),arrs(mt))
          thread_alloc=.FALSE.
          first_call=.FALSE.
        ENDIF
        !$omp critical
        fft_syncA=fft_syncA+1
        !$omp end critical
        DO WHILE (fft_syncA/=nt)
          !$omp flush(fft_syncA)
        ENDDO
        fftw_ny=ny
      ENDIF
!-----------------------------------------------------------------------
!     the plan size is never recomputed, check that it has not changed.
!-----------------------------------------------------------------------
      IF (fftw_ny/=ny) THEN
        CALL nim_stop("FFTW plan size changed!")
      ENDIF
!-----------------------------------------------------------------------
!     setup plans if needed - use fftw_alloc for arrays alignment
!-----------------------------------------------------------------------
      IF (.NOT.thread_alloc(it)) THEN
        thread_alloc(it)=.TRUE.
#ifdef __xlf
        ALLOCATE(arrs(it)%rarr(ny,vecstride))
        ALLOCATE(arrs(it)%carr(nm_tot,vecstride))
#else
        ptr =  fftw_alloc_real(int(vecstride*ny, C_SIZE_T))
        CALL C_F_POINTER(ptr, arrs(it)%rarr, [ny,vecstride])
        ptr = fftw_alloc_complex(int(vecstride*nm_tot, C_SIZE_T))
        CALL C_F_POINTER(ptr, arrs(it)%carr, [nm_tot,vecstride])
#endif
        rdim(1)=ny
        cdim(1)=nm_tot
        !$omp critical
        plan_r2c(it) = fftw_plan_many_dft_r2c(                          &
     &    1_i4, rdim, vecstride, arrs(it)%rarr, rdim, 1_i4, ny,         &
     &    arrs(it)%carr, cdim, 1_i4, nm_tot, FFTW_MEASURE)
        plan_c2r(it) = fftw_plan_many_dft_c2r(                          &
     &    1_i4, rdim, vecstride, arrs(it)%carr, cdim, 1_i4, nm_tot,     &
     &    arrs(it)%rarr, rdim, 1_i4, ny, FFTW_MEASURE)
        !$omp end critical
      ENDIF
!-----------------------------------------------------------------------
!     excute the plan in chunks of vecstride*nt.
!-----------------------------------------------------------------------
      ichunk = 1
      DO WHILE (ichunk <= nx)
        iend = ichunk + vecstride - 1
        vs=vecstride
        IF (iend > nx) THEN
          iend=nx
          vs=iend-ichunk+1
          IF (iflag==1) THEN
            arrs(it)%rarr=0._r8
          ELSE
            arrs(it)%carr=0._r8
          ENDIF
        ENDIF
!-----------------------------------------------------------------------
!       load data into aligned arrays.
!-----------------------------------------------------------------------
        IF (iflag==1) THEN
          arrs(it)%rarr(:,1:vs)=TRANSPOSE(re(ichunk:iend,1:ny))
        ELSE
          IF (dealiase>=3) THEN
            arrs(it)%carr(1:nmodes,1:vs)=                               &
     &        TRANSPOSE(comp(ichunk:iend,1:nmodes))
            arrs(it)%carr(nmodes+1:nm_tot,:)=0_i4
          ELSE
            arrs(it)%carr(:,1:vs)=TRANSPOSE(comp(ichunk:iend,1:nm_tot))
          ENDIF
        ENDIF
!-----------------------------------------------------------------------
!       execute the ffts.
!-----------------------------------------------------------------------
        IF (iflag==1) THEN
          CALL fftw_execute_dft_r2c(plan_r2c(it), arrs(it)%rarr,        &
     &                              arrs(it)%carr)
        ELSE
          CALL fftw_execute_dft_c2r(plan_c2r(it), arrs(it)%carr,        &
     &                              arrs(it)%rarr)
        ENDIF
!-----------------------------------------------------------------------
!       unload data from aligned arrays.
!-----------------------------------------------------------------------
        IF (iflag==1) THEN
          comp(ichunk:iend,1:nmodes)=                                   &
     &      TRANSPOSE(arrs(it)%carr(1:nmodes,1:vs))
        ELSE
          re(ichunk:iend,:)=TRANSPOSE(arrs(it)%rarr(:,1:vs))
        ENDIF
        ichunk = ichunk + vecstride
      ENDDO
!-----------------------------------------------------------------------
!     apply the normalization.
!-----------------------------------------------------------------------
      IF (iflag==1) THEN
        comp(:,1:nmodes)=comp(:,1:nmodes)/ny
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE nim_fftw
#endif /* HAVE_FFTW3 */
!-----------------------------------------------------------------------
!     close module.                                                     
!-----------------------------------------------------------------------
      END MODULE fft2d_mod

