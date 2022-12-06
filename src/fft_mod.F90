!-----------------------------------------------------------------------
!     $Id: fft_mod.F90 4620 2015-12-03 19:18:40Z jking $
!     module containing fast Fourier transform routines for use with    
!     nimrod.  this is a hybrid version having the newer wrapper        
!     retrofitted with Zoran Mikic's FFT routines.                      
!     (CRS and SJP--last revised 8/14/98 to make the ffts and           
!     configuration-space operations scale with nlayers.  note that     
!     the passed parameters, nf and ng are different than the           
!     previous nx and ny.)                                              
!                                                                       
!     modified 12/10/98 to use a complex work array in cpftv.  CRS      
!                                                                       
!     modified 1/13/99 to handle complex Fourier coefficient arrays, as 
!     directly as possible.  CRS                                        
!                                                                       
!     modified 9/16/99 for separating vector components and Fourier     
!     components into separate array indices and having the toroidal    
!     index last in the returned real-space arrays.                     
!                                                                       
!     the option to perform ffts in nimrod's data structures without    
!     dealiasing is being added.  CRS, 5/30/08                          
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.                                                
!-----------------------------------------------------------------------
!     0. fft_mod.                                                       
!     1. fft_nim.                                                       
!-----------------------------------------------------------------------
#include "config.f"
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      MODULE fft_mod 
      USE local 
#ifdef HAVE_FFTW3
      USE fft2d_mod, ONLY: fft2d,nim_fftw
#else
      USE fft2d_mod, ONLY: fft2d
#endif 
#ifdef HAVE_OPENMP
      USE omp_lib
#endif 
      IMPLICIT NONE 
                                                                        
      CONTAINS 
!-----------------------------------------------------------------------
!     subprogram 1. fft_nim.                                            
!     performs array packing, and calls the actual fft routine.         
!     this version uses a complex array for the passed Fourier          
!     components.  the data is three-dimensional, and only the last is  
!     transformed.                                                      
!                                                                       
!     the arguments for this subprogram are                             
!                                                                       
!	direction:  either 'forward' or 'inverse'.  the former sums           
!		functions of configuration space with exp{-inx}, and                 
!		divides by the number of cells for normalization.  the               
!		latter sums Fourier coefficients with exp{+inx}.                     
!                                                                       
!	nf:  dimension of the non-transformed direction of the array          
!		of Fourier coefficients.                                             
!                                                                       
!	nr:  dimension of the non-transformed direction of the array          
!		of configuration-space data.  If nr<nf, different                    
!		segments of the non-transformed dimension are allocated              
!		to the different processor layers, and the config-space              
!		data on a layer does not cover a full grid-block.                    
!                                                                       
!	nphi: the number of cells of the transformable dimension.
!                                                                       
!	nq:  the number of quantities (vector components) represented         
!		in the data field.                                                   
!                                                                       
!	f_coef: complex 3D array (1:nq,1:nf,(N+1)), where N=nphi/3,
!               (see below if dealiase /= 3)        
!		representing the n>=0 Fourier coefficients of nq real                
!		functions of configuration space.                                    
!                                                                       
!	re: 3D array (1:nq,1:nr,1:nphi) for the data as                    
!		a function of configuration space.                                   
!                                                                       
!       dealiase: if this optional argument is present and set to
!		zero, then the maximum Fourier coefficient is                       
!		N=nphi/2 in the fft routine and in the last array                 
!		dimension of f_coef.  this implies a different layer                 
!		decomposition for the complex coefficients and for the               
!		real data. otherwise it is used to dealiase the 
!               such that N=nphi/dealiase.
!
!     the last two arrays are now assumed-size dummy arrays, so that the
!     non-transformed directions are packed together for better         
!     efficiency.                                                       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE fft_nim(direction,nf,nr,nphi,nq,f_coef,re,dealiase)
      USE time 
      USE mpi_nim 
      USE pardata
      IMPLICIT NONE 

      CHARACTER(*), INTENT(IN) :: direction 
      INTEGER(i4), INTENT(IN) :: nf,nr,nphi,nq 
      COMPLEX(r8), DIMENSION(nq,nf,*), INTENT(INOUT) :: f_coef 
      REAL(r8), DIMENSION(nq,nr,*), INTENT(INOUT) :: re 
      INTEGER(i4), INTENT(IN), OPTIONAL :: dealiase
                                                                        
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: f_ctmp 
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: comp 
      REAL(r8) :: timestart_fft,timeend_fft,timestart_comm,timeend_comm
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: recvcounts,displs,      &
     &             scounts,sdispls,rcounts,rdispls                      
      INTEGER(i4) :: isign,lphi,iq,im,im0,idp,idm,ierror,rep,nfl,       &
     &               nm_tot,il,nl,nm,rem,ip,ifst,ifen,npo,jfc,          &
     &               nmodes,ifc,it,nt,soff,stoff,smoff,sfc,nmil,imst
      LOGICAL :: dealflag
!-----------------------------------------------------------------------
!     start timer and allocate arrays.                                  
!-----------------------------------------------------------------------
      CALL timer(timestart_fft)
#ifdef HAVE_OPENMP
      it=OMP_GET_THREAD_NUM()
      nt=OMP_GET_NUM_THREADS()
#else
      it=0
      nt=1
#endif
      nl=nprocs/nprocs_layer                   ! nlayers
      lphi=NINT(LOG(REAL(nphi,r8))/LOG(2._r8)) ! log2(nphi)
      nm_tot=nphi/3+1                          ! nmodes total 
      dealflag=.TRUE.
      IF (PRESENT(dealiase)) THEN 
        IF (dealiase<3) THEN 
          nm_tot=nphi/2+1 
          dealflag=.FALSE.
        ELSE
#ifdef HAVE_FFTW3
          nm_tot=nphi/dealiase+1
#else
          nm_tot=nphi/3+1
#endif
        ENDIF
      ENDIF
      nm=nm_tot/nl                          ! nmodes for this layer
      nfl=nf/nl                             ! number of RZ pts per layer
      rem=MODULO(nm_tot,nl)                 ! catch uneven layer decomp 
      rep=MODULO(nf,nl)                     ! catch uneven RZ pt decomp
      nmodes=mode_hi-mode_lo+1              ! nmodes for this layer
!-----------------------------------------------------------------------
!     check direction.                                                  
!-----------------------------------------------------------------------
      IF (direction/='forward'.AND.direction/='inverse') CALL nim_stop  &
     &  ('fft_nim does not recognize direction '//direction//'.')
!-----------------------------------------------------------------------
!     pack arrays.                                                      
!-----------------------------------------------------------------------
      dir_if: IF (direction=='forward') THEN 
!-----------------------------------------------------------------------
!       if nphi<=2, the only component is n=0.  there should only be    
!       one processor-layer in these cases.                             
!-----------------------------------------------------------------------
        IF (nphi<=2) THEN 
          f_coef(:,:,1)=re(:,:,1) 
          IF (it==0) THEN
            CALL timer(timeend_fft) 
            time_fft = time_fft + timeend_fft-timestart_fft 
          ENDIF
          RETURN 
        ENDIF 
        isign=1 
        ALLOCATE(comp(nq,nr,nm_tot)) 
      ELSE dir_if 
!-----------------------------------------------------------------------
!       if nphi<=2, the only component is n=0.                          
!-----------------------------------------------------------------------
        IF (nphi<=2) THEN 
          DO ip=1,nphi 
            re(:,:,ip)=f_coef(:,:,1) 
          ENDDO 
          IF (it==0) THEN
            CALL timer(timeend_fft) 
            time_fft = time_fft + timeend_fft-timestart_fft 
          ENDIF
          RETURN 
        ENDIF 
        isign=-1
        ALLOCATE(comp(nq,nr,nm_tot))
!-----------------------------------------------------------------------
!       there are 3 possibilities:  1) Fourier coefficients are located 
!       on separate layers but all layers need all of the configuration-
!       space data, 2) coefficients on separate layers and config-space 
!       data is broken up into separate portions across the poloidal    
!       plane (within a block), 3) there is only one layer.             
!                                                                       
!       1) if different Fourier coefficients are located on different   
!       processor layers, collect all before packing.                   
!       recvcounts = # of values from each proc across layers           
!       displs = offsets into comp array for each proc's section        
!-----------------------------------------------------------------------
        layer_if: IF (nl>1.AND.nr==nf.AND.nt==1) THEN
          ALLOCATE(recvcounts(0:nl-1),displs(0:nl-1)) 
          displs(0)=0 
          DO il=0,nl-1 
            IF (il<rem) THEN 
              recvcounts(il)=nq*nf*(nm+1) 
            ELSE 
              recvcounts(il)=nq*nf*nm 
            ENDIF 
            IF (il>0) displs(il)=displs(il-1)+recvcounts(il-1) 
          ENDDO 
          CALL timer(timestart_comm)
          CALL mpi_allgatherv(f_coef,recvcounts(ilayer),mpi_nim_comp,   &
     &         comp,recvcounts,displs,mpi_nim_comp,                     &
     &         comm_mode,ierror)
          CALL timer(timeend_comm)
          time_fftcm = time_fftcm + timeend_comm-timestart_comm
          DEALLOCATE(recvcounts,displs) 
!-----------------------------------------------------------------------
!       1a) threaded communication.
!       Similar to (1) but uses f_coefth and compth from pardata,
!       which are shared between threads to reduce the mpi_allgatherv
!       to be between only the master threads. The master thread
!       allocates and communicates, and all threads load/unload their
!       contributions individually.  Four barriers are required to
!       synchronize this operation between threads.
!       PRE- this block of code (unused)
!-----------------------------------------------------------------------
        ELSE IF (nl>1.AND.nr==nf) THEN layer_if
          CALL nim_stop("Houston, we have a (FFT) problem.")
          ! master
          IF (it==0) THEN
            ALLOCATE(recvcounts(0:nl-1),displs(0:nl-1),                 &
     &                f_coefth(nq*nf*nt*nmodes),compth(nq*nr*nt*nm_tot))
            displs(0)=0
            DO il=0,nl-1
              IF (il<rem) THEN
                recvcounts(il)=nq*nf*(nm+1)*nt
              ELSE
                recvcounts(il)=nq*nf*nm*nt
              ENDIF
              IF (il>0) displs(il)=displs(il-1)+recvcounts(il-1)
            ENDDO
          ENDIF
          ! end master
          !f_coefth(:,it*nf+1:(it+1)*nf,:)=f_coef(:,:,1:nmodes)
          ! master
          IF (it==0) THEN
            CALL timer(timestart_comm)
            CALL mpi_allgatherv(f_coefth,recvcounts(ilayer),mpi_nim_comp, &
     &           compth,recvcounts,displs,mpi_nim_comp,comm_mode,ierror)
            CALL timer(timeend_comm)
            time_fftcm = time_fftcm + timeend_comm-timestart_comm
          ENDIF
          ! end master
          !comp(:,:,:)=compth(:,it*nr+1:(it+1)*nr,:)
          ! master
          IF (it==0) THEN
            DEALLOCATE(recvcounts,displs,f_coefth,compth)
          ENDIF
          ! end master
!-----------------------------------------------------------------------
!       2) transpose type communication.  First make a 1D
!       array holding contiguous information for communication and      
!       the send & recv displacement and count arrays.
!-----------------------------------------------------------------------
        ELSE IF (nf>nr.AND.nt==1) THEN layer_if
          ALLOCATE(f_ctmp(nq*nf*nmodes)) 
          ALLOCATE(sdispls(0:nl-1),scounts(0:nl-1),                     &
     &             rdispls(0:nl-1),rcounts(0:nl-1))                     
          sdispls(0)=0 
          rdispls(0)=0 
          jfc=1 
          DO il=0,nl-1 
            ifst=il*nfl+1+MIN(rep,il) 
            ifen=(il+1)*nfl+MIN(rep,il+1_i4) 
            npo=ifen-ifst+1 
            scounts(il)=npo*nq*nmodes 
            IF (il<rem) THEN 
              rcounts(il)=(nm+1)*nq*nr 
            ELSE 
              rcounts(il)=nm*nq*nr 
            ENDIF 
            IF (il>0) THEN 
              sdispls(il)=sdispls(il-1)+scounts(il-1) 
              rdispls(il)=rdispls(il-1)+rcounts(il-1) 
            ENDIF 
            DO im=1,nmodes 
              DO ifc=ifst,ifen 
                f_ctmp(jfc:jfc+nq-1)=f_coef(1:nq,ifc,im) 
                jfc=jfc+nq 
              ENDDO 
            ENDDO 
          ENDDO 
          CALL timer(timestart_comm)
          CALL mpi_alltoallv(f_ctmp,scounts,sdispls,mpi_nim_comp,       &
     &         comp,rcounts,rdispls,mpi_nim_comp,                       &
     &         comm_mode,ierror)
          CALL timer(timeend_comm)
          time_fftcm = time_fftcm + timeend_comm-timestart_comm
          DEALLOCATE(f_ctmp,scounts,sdispls,rcounts,rdispls) 
!-----------------------------------------------------------------------
!       2a) threaded communication.
!       Similar to (2) but uses f_coefth and compth from pardata,
!       which are shared between threads to reduce the mpi_alltoallv
!       to be between only the master threads. The master thread
!       allocates and communicates, and all threads load/unload their
!       contributions individually.  Four barriers are required to
!       synchronize this operation between threads.
!-----------------------------------------------------------------------
        ELSE IF (nf>nr) THEN layer_if
          IF (it==0) THEN
            ALLOCATE(f_coefth(nq*nf*nt*nmodes),compth(nq*nr*nt*nm_tot))
            ALLOCATE(sdispls(0:nl-1),scounts(0:nl-1),                   &
     &               rdispls(0:nl-1),rcounts(0:nl-1))
            sdispls(0)=0
            rdispls(0)=0
            DO il=0,nl-1
              npo=nfl
              IF (il<rep) npo=npo+1
              scounts(il)=npo*nq*nt*nmodes
              IF (il<rem) THEN
                rcounts(il)=(nm+1)*nq*nt*nr
              ELSE
                rcounts(il)=nm*nq*nt*nr
              ENDIF
              IF (il>0) THEN
                sdispls(il)=sdispls(il-1)+scounts(il-1)
                rdispls(il)=rdispls(il-1)+rcounts(il-1)
              ENDIF
            ENDDO
          ENDIF
          !$omp critical
          fft_sync1=fft_sync1+1
          !$omp end critical
          DO WHILE (fft_sync1/=nt)
            !$omp flush(fft_sync1)
          ENDDO
          soff=1
          DO il=0,nl-1
            ifst=il*nfl+1+MIN(rep,il)
            ifen=(il+1)*nfl+MIN(rep,il+1_i4)
            npo=ifen-ifst+1
            stoff=soff+nq*npo*nmodes*it
            DO im=1,nmodes
              smoff=stoff+nq*npo*(im-1)
              DO ifc=ifst,ifen
                sfc=smoff+(ifc-ifst)*nq
                f_coefth(sfc:sfc+nq-1)=f_coef(1:nq,ifc,im)
              ENDDO
            ENDDO
            soff=soff+nq*npo*nmodes*nt
          ENDDO
          !$omp critical
          fft_sync2=fft_sync2+1
          !$omp end critical
          DO WHILE (fft_sync2/=nt)
            !$omp flush(fft_sync2)
          ENDDO
          IF (it==0) THEN
            CALL timer(timestart_comm)
            CALL mpi_alltoallv(f_coefth,scounts,sdispls,mpi_nim_comp,   &
     &           compth,rcounts,rdispls,mpi_nim_comp,comm_mode,ierror)
            CALL timer(timeend_comm)
            time_fftcm = time_fftcm + timeend_comm-timestart_comm
            fft_sync1=0_i4; fft_sync4=0_i4
          ENDIF
          !$omp critical
          fft_sync3=fft_sync3+1
          !$omp end critical
          DO WHILE (fft_sync3/=nt)
            !$omp flush(fft_sync3)
          ENDDO
          soff=1
          imst=1
          DO il=0,nl-1
            nmil=nm 
            IF (il<rem) nmil=nmil+1
            sfc=soff+nq*nr*nmil*it
            DO im=imst,nmil+imst-1
              DO ifc=1,nr
                comp(1:nq,ifc,im)=compth(sfc:sfc+nq-1)
                sfc=sfc+nq
              ENDDO
            ENDDO
            imst=imst+nmil
            soff=soff+nq*nr*nmil*nt
          ENDDO
          !$omp critical
          fft_sync4=fft_sync4+1
          !$omp end critical
          DO WHILE (fft_sync4/=nt)
            !$omp flush(fft_sync4)
          ENDDO
          IF (it==0) THEN
            DEALLOCATE(f_coefth,compth,scounts,sdispls,rcounts,rdispls)
            fft_sync2=0_i4; fft_sync3=0_i4
          ENDIF
        ENDIF layer_if
      ENDIF dir_if 
!-----------------------------------------------------------------------
!     perform the transform.                                            
!-----------------------------------------------------------------------
#ifdef HAVE_FFTW3
      IF (nl==1) THEN 
        CALL nim_fftw(re,f_coef,nr*nq,nphi,isign,dealiase) 
      ELSE 
        CALL nim_fftw(re,comp,  nr*nq,nphi,isign,dealiase) 
      ENDIF 
#else /* HAVE_FFTW3 */
      IF (nl==1) THEN 
        CALL fft2d(re,f_coef,nr*nq,lphi,isign,dealflag) 
      ELSE 
        CALL fft2d(re,comp  ,nr*nq,lphi,isign,dealflag) 
      ENDIF 
#endif /* HAVE_FFTW3 */
!-----------------------------------------------------------------------
!     put Fourier component data into nimrod form if this is a          
!     forward transform.                                                
!-----------------------------------------------------------------------
      IF (direction=='forward') THEN 
!-----------------------------------------------------------------------
!       for multi-layer runs, collect all portions of the block in a 1D 
!       array, then sort.  the first part of the if block rearranges    
!       from layered configuration-space data (splitting the poloidal   
!       directions) to layered Fourier coefficients (splitting          
!       components).  the second part goes from unsplit config-space    
!       data to split coefficients.                                     
!-----------------------------------------------------------------------
        layer_if2: IF (nf>nr.AND.nt==1) THEN
          ALLOCATE(f_ctmp(nq*nf*nmodes)) 
          ALLOCATE(sdispls(0:nl-1),scounts(0:nl-1),                     &
     &             rdispls(0:nl-1),rcounts(0:nl-1))                     
          rdispls(0)=0 
          sdispls(0)=0 
          DO il=0,nl-1 
            ifst=il*nfl+1+MIN(rep,il) 
            ifen=(il+1)*nfl+MIN(rep,il+1_i4) 
            npo=ifen-ifst+1 
            rcounts(il)=npo*nq*nmodes 
            IF (il<rem) THEN 
              scounts(il)=(nm+1)*nq*nr 
            ELSE 
              scounts(il)=nm*nq*nr 
            ENDIF 
            IF (il>0) THEN 
              rdispls(il)=rdispls(il-1)+rcounts(il-1) 
              sdispls(il)=sdispls(il-1)+scounts(il-1) 
            ENDIF 
          ENDDO 
!-----------------------------------------------------------------------
!         gather this layer's Fourier components from all block         
!         portions in comp.                                             
!-----------------------------------------------------------------------
          CALL timer(timestart_comm)
          CALL mpi_alltoallv(comp,scounts,sdispls,mpi_nim_comp,         &
     &         f_ctmp,rcounts,rdispls,mpi_nim_comp,                     &
     &         comm_mode,ierror)
          CALL timer(timeend_comm)
          time_fftcm = time_fftcm + timeend_comm-timestart_comm
          DEALLOCATE(scounts,sdispls,rcounts,rdispls) 
!-----------------------------------------------------------------------
!         select the data for this Fourier component layer.             
!-----------------------------------------------------------------------
          jfc=1 
          DO il=0,nl-1 
            ifst=il*nfl+1+MIN(rep,il) 
            ifen=(il+1)*nfl+MIN(rep,il+1_i4) 
            npo =ifen-ifst+1 
            DO im=1,nmodes 
              DO ifc=ifst,ifen 
                f_coef(1:nq,ifc,im)=f_ctmp(jfc:jfc+nq-1) 
                jfc=jfc+nq 
              ENDDO 
            ENDDO 
          ENDDO 
          DEALLOCATE(f_ctmp) 
          IF (ilayer==0) f_coef(:,:,1)=REAL(f_coef(:,:,1),r8) 
!-----------------------------------------------------------------------
!       Case for threaded communication by the master node.
!-----------------------------------------------------------------------
        ELSE IF (nf>nr) THEN
          IF (it==0) THEN
            ALLOCATE(f_coefth(nq*nf*nt*nmodes),compth(nq*nr*nt*nm_tot))
            ALLOCATE(sdispls(0:nl-1),scounts(0:nl-1),                   &
     &               rdispls(0:nl-1),rcounts(0:nl-1))
            sdispls(0)=0
            rdispls(0)=0
            DO il=0,nl-1
              npo=nfl
              IF (il<rep) npo=npo+1
              rcounts(il)=npo*nq*nt*nmodes
              IF (il<rem) THEN
                scounts(il)=(nm+1)*nq*nt*nr
              ELSE
                scounts(il)=nm*nq*nt*nr
              ENDIF
              IF (il>0) THEN
                sdispls(il)=sdispls(il-1)+scounts(il-1)
                rdispls(il)=rdispls(il-1)+rcounts(il-1)
              ENDIF
            ENDDO
          ENDIF
          !$omp critical
          fft_sync1=fft_sync1+1
          !$omp end critical
          DO WHILE (fft_sync1/=nt)
            !$omp flush(fft_sync1)
          ENDDO
          soff=1
          imst=1
          DO il=0,nl-1
            nmil=nm
            IF (il<rem) nmil=nmil+1
            sfc=soff+nq*nr*nmil*it
            DO im=imst,nmil+imst-1
              DO ifc=1,nr
                compth(sfc:sfc+nq-1)=comp(1:nq,ifc,im)
                sfc=sfc+nq
              ENDDO
            ENDDO
            imst=imst+nmil
            soff=soff+nq*nr*nmil*nt
          ENDDO
          !$omp critical
          fft_sync2=fft_sync2+1
          !$omp end critical
          DO WHILE (fft_sync2/=nt)
            !$omp flush(fft_sync2)
          ENDDO
          IF (it==0) THEN
            CALL timer(timestart_comm)
            CALL mpi_alltoallv(compth,scounts,sdispls,mpi_nim_comp,     &
     &           f_coefth,rcounts,rdispls,mpi_nim_comp,comm_mode,ierror)
            CALL timer(timeend_comm)
            time_fftcm = time_fftcm + timeend_comm-timestart_comm
            fft_sync1=0_i4; fft_sync4=0_i4
          ENDIF
          !$omp critical
          fft_sync3=fft_sync3+1
          !$omp end critical
          DO WHILE (fft_sync3/=nt)
            !$omp flush(fft_sync3)
          ENDDO
          soff=1
          DO il=0,nl-1
            ifst=il*nfl+1+MIN(rep,il)
            ifen=(il+1)*nfl+MIN(rep,il+1_i4)
            npo=ifen-ifst+1
            stoff=soff+nq*npo*nmodes*it
            DO im=1,nmodes
              smoff=stoff+nq*npo*(im-1)
              DO ifc=ifst,ifen
                sfc=smoff+(ifc-ifst)*nq
                f_coef(1:nq,ifc,im)=f_coefth(sfc:sfc+nq-1)
              ENDDO
            ENDDO
            soff=soff+nq*npo*nmodes*nt
          ENDDO
          !$omp critical
          fft_sync4=fft_sync4+1
          !$omp end critical
          DO WHILE (fft_sync4/=nt)
            !$omp flush(fft_sync4)
          ENDDO
          IF (it==0) THEN
            DEALLOCATE(f_coefth,compth,scounts,sdispls,rcounts,rdispls)
            fft_sync2=0_i4; fft_sync3=0_i4
          ENDIF
          IF (ilayer==0) f_coef(:,:,1)=REAL(f_coef(:,:,1),r8)
!-----------------------------------------------------------------------
!       complete config-space data to component on layers.              
!-----------------------------------------------------------------------
        ELSE IF (nl>1) THEN layer_if2 
          f_coef(1:nq,1:nf,1:nmodes)=comp(:,:,mode_lo:mode_hi) 
          IF (ilayer==0) f_coef(:,:,1)=REAL(f_coef(:,:,1),r8) 
        ENDIF layer_if2 
      ENDIF 
      DEALLOCATE(comp) 
!-----------------------------------------------------------------------
!     complete timing.                                                  
!-----------------------------------------------------------------------
      IF (it==0) THEN
        CALL timer(timeend_fft) 
        time_fft = time_fft + timeend_fft-timestart_fft 
      ENDIF
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE fft_nim
!-----------------------------------------------------------------------
!     close module.                                                     
!-----------------------------------------------------------------------
      END MODULE fft_mod 
