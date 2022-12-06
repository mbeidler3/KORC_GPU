!-----------------------------------------------------------------------
!     Routines useful for exchanging data from a local lagr_quad_type to
!     a global lagr_quad_type on the rblocks                            
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.                                                
!-----------------------------------------------------------------------
!     0.  exchange_mod.                                                 
!     1.  exchange_init                                                                                          
!     3.  exch_lagr_2D                                                  
!     4.  exch_lagr_3D                                                                                                 
!-----------------------------------------------------------------------
!     subprogram 0. exchange_mod                                        
!-----------------------------------------------------------------------
      MODULE exchange_mod 
      USE local 
      USE global 
      USE input 
      USE rblock_type_mod 
      USE lagr_quad_mod 
      USE fields 
      USE mpi_nim 
      USE pardata 
      USE send_rblock_mod 
      IMPLICIT NONE 
      
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: x0,y0,mxarr,myarr 
      REAL(r8), PRIVATE, PARAMETER :: z99 = -99. 
      COMPLEX(r8), PRIVATE, PARAMETER :: zc99=(-99.,-99.) 
 
      INTERFACE exch_lagr 
        MODULE PROCEDURE exch_lagr_2D,exch_lagr_3D 
      END INTERFACE 
 
      CONTAINS 
!-----------------------------------------------------------------------
!     subprogram 1. exchange_init                                       
!     Initialize the arrays used for data exchange                      
!-----------------------------------------------------------------------
      SUBROUTINE exchange_init 
      IMPLICIT NONE 
      INTEGER(i4) :: nb,ibl,i,ierr 
 
! Put mx and my's from each block into a global array                   
      nb=nrbl_total 
      IF (.NOT.ALLOCATED(x0))ALLOCATE(x0(nb),y0(nb),mxarr(nb),myarr(nb)) 
      DO ibl = 1,nrbl 
        x0(ibl)=rb(ibl)%mx 
        y0(ibl)=rb(ibl)%my 
      ENDDO 
      IF (nprocs>1) THEN 
         CALL mpi_allgather(x0(1:nrbl),nrbl,mpi_nim_int,mxarr,nrbl,     &
     &                                 mpi_nim_int,comm_layer,ierr)     
         CALL mpi_allgather(y0(1:nrbl),nrbl,mpi_nim_int,myarr,nrbl,     &
     &                                 mpi_nim_int,comm_layer,ierr)     
      ELSE 
         mxarr=x0 
         myarr=y0 
      ENDIF 
 
! For each block ibl, calculate its starting location                   
      DO ibl = 1,nxbl*nybl 
       i = mod(ibl-1,nybl) 
       IF(ibl>nybl) THEN 
          x0(ibl) = SUM(mxarr(1:ibl-nybl:nybl)) 
       ELSE 
          x0(ibl) = 0 
       ENDIF 
       y0(ibl) = SUM(myarr(1:i)) 
      ENDDO 
 
      RETURN 
      END SUBROUTINE exchange_init 
!-----------------------------------------------------------------------
!     subprogram 3. exch_lagr_2D(laq1,laq2,ibl_id)                      
!     Place the data from a block into a 1D array, communicate it if    
!     necessary, and then unpack it into a global array                 
!-----------------------------------------------------------------------
      SUBROUTINE exch_lagr_2D(laq1,laq2,ibl_id) 
 
      IMPLICIT NONE 
      INTEGER(i4) :: ibl_id 
      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq1 
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq2 
      INTEGER(i4) :: mxl,myl,ibl 
      INTEGER(i4) :: status(mpi_status_size) 
      INTEGER(i4) :: i,xn,yn,qty,n_side,n_int,length,ierr 
      INTEGER(i4) :: myzfs,myzfsh,myzfsv,myzfsi,mylength 
      INTEGER(i4), DIMENSION(nprocs) :: zfs,zfsh,zfsv,zfsi,iblar 
      REAL(r8), ALLOCATABLE :: data1(:),data2(:,:) 
      LOGICAL, ALLOCATABLE :: lfs(:,:,:) 
      LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: lfsh,lfsv,lfsi 
!-----------------------------------------------------------------------
!     initialize global lagrange quad variables                         
!-----------------------------------------------------------------------
      laq2%title = laq1%title;      laq2%name = laq1%name 
      qty    = laq2%nqty
      n_side = laq2%n_side
      n_int  = laq2%n_int
                                                                        
      IF (nprocs > 1) THEN 
        CALL mpi_allgather(ibl_id,1,mpi_nim_int,iblar,                  &
     &                   1,mpi_nim_int,comm_nimrod,ierr)                
      ELSE 
        iblar=ibl_id 
      ENDIF 
      mxl=MAXVAL(mxarr);  myl=MAXVAL(myarr) 
      ALLOCATE(lfs(1:qty,0:mxl,0:myl)) 
                                                ! Filter for PACK/UNPACK
      lfs=.true. 
!-----------------------------------------------------------------------
!     Exchange the data from block_decomposed data to single block data 
!-----------------------------------------------------------------------
      ! Poly_degree=1                                                   
      IF (n_side==0) THEN 
        myzfs=SIZE(laq1%fs) 
        IF (nprocs > 1) THEN 
          CALL mpi_allgather(myzfs,1,mpi_nim_int,zfs,                   &
     &                       1,mpi_nim_int,comm_nimrod,ierr)            
        ELSE 
          zfs=myzfs 
        ENDIF 
        length=MAXVAL(zfs) 
        ALLOCATE(data1(length)); data1=0 
        data1(1:myzfs)=PACK(laq1%fs(:,:,:),.true.) 
        ALLOCATE(data2(length,nprocs)) 
        IF (nprocs > 1) THEN 
          CALL mpi_allgather(data1,length,mpi_nim_real,data2,           &
     &                       length,mpi_nim_real,comm_nimrod,ierr)      
        ELSE 
          data2(:,1)=data1 
        ENDIF 
        DEALLOCATE(data1) 
        parallel_block_loop1: DO i=1,nprocs 
          ibl=iblar(i) 
          xn = x0(ibl) 
          yn = y0(ibl) 
          laq2%fs(1:qty,xn:xn+mxarr(ibl),yn:yn+myarr(ibl))=             &
     &     UNPACK(data2(1:zfs(i),i),lfs(1:qty,0:mxarr(ibl),             &
     &            0:myarr(ibl)),z99)                                    
        ENDDO parallel_block_loop1 
        DEALLOCATE(data2,lfs) 
                                                                        
      ! Poly_degree>1                                                   
      ELSE 
 
        ALLOCATE(lfsh(1:qty,1:n_side,1:mxl,0:myl)) 
        ALLOCATE(lfsv(1:qty,1:n_side,0:mxl,1:myl)) 
        ALLOCATE(lfsi(1:qty,1:n_int,1:mxl,1:myl)) 
        lfsh=.true.;  lfsv=.true.;  lfsi=.true. 
        myzfs=SIZE(laq1%fs);    myzfsh=SIZE(laq1%fsh); 
        myzfsv=SIZE(laq1%fsv);  myzfsi=SIZE(laq1%fsi) 
        mylength=myzfs+myzfsh+myzfsv+myzfsi 
        IF (nprocs > 1) THEN 
          CALL mpi_allgather(myzfs,1,mpi_nim_int,zfs,                   &
     &                     1,mpi_nim_int,comm_nimrod,ierr)              
          CALL mpi_allgather(myzfsh,1,mpi_nim_int,zfsh,                 &
     &                     1,mpi_nim_int,comm_nimrod,ierr)              
          CALL mpi_allgather(myzfsv,1,mpi_nim_int,zfsv,                 &
     &                     1,mpi_nim_int,comm_nimrod,ierr)              
          CALL mpi_allgather(myzfsi,1,mpi_nim_int,zfsi,                 &
     &                     1,mpi_nim_int,comm_nimrod,ierr)              
          CALL mpi_allreduce(mylength,length,1,mpi_nim_int,             &
     &                     MPI_MAX,comm_nimrod,ierr)                    
        ELSE 
          zfs=myzfs;  zfsh=myzfsh;  zfsv=myzfsv; zfsi=myzfsi 
          length=mylength 
        ENDIF 
        ALLOCATE(data1(length)); data1=0 
        data1(1:myzfs)=PACK(laq1%fs(:,:,:),.true.) 
        data1(myzfs+1:myzfs+myzfsh)=PACK(laq1%fsh(:,:,:,:),.true.) 
        data1(myzfs+myzfsh+1:myzfs+myzfsh+myzfsv)=                      &
     &        PACK(laq1%fsv(:,:,:,:),.true.)                            
        data1(myzfs+myzfsh+myzfsv+1:myzfs+myzfsh+myzfsv+myzfsi)=        &
     &        PACK(laq1%fsi(:,:,:,:),.true.)                            
        ALLOCATE(data2(length,nprocs)) 
        IF (nprocs > 1) THEN 
          CALL mpi_allgather(data1,length,mpi_nim_real,data2,           &
     &                     length,mpi_nim_real,comm_nimrod,ierr)        
        ELSE 
          data2(:,1)=data1 
        ENDIF 
        DEALLOCATE(data1) 
        parallel_block_loop2: DO i=1,nprocs 
          ibl=iblar(i) 
          xn = x0(ibl) 
          yn = y0(ibl) 
          laq2%fs(1:qty,xn:xn+mxarr(ibl),yn:yn+myarr(ibl))=             &
     &     UNPACK(data2(1:zfs(i),i),                                    &
     &                        lfs(1:qty,0:mxarr(ibl),0:myarr(ibl)),z99) 
          laq2%fsh(1:qty,1:n_side,xn+1:xn+mxarr(ibl),yn:yn+myarr(ibl))= &
     &     UNPACK(data2(zfs(i)+1:zfs(i)+zfsh(i),i),                     &
     &              lfsh(1:qty,1:n_side,1:mxarr(ibl),0:myarr(ibl)),z99) 
          laq2%fsv(1:qty,1:n_side,xn:xn+mxarr(ibl),yn+1:yn+myarr(ibl))= &
     &     UNPACK(data2(zfs(i)+zfsh(i)+1:zfs(i)+zfsh(i)+zfsv(i),i),     &
     &              lfsv(1:qty,1:n_side,0:mxarr(ibl),1:myarr(ibl)),z99) 
          laq2%fsi(1:qty,1:n_int,xn+1:xn+mxarr(ibl),yn+1:yn+myarr(ibl))=&
     &     UNPACK(data2(zfs(i)+zfsh(i)+zfsv(i)+1:                       &
     &                                zfs(i)+zfsh(i)+zfsv(i)+zfsi(i),i),&
     &               lfsi(1:qty,1:n_int,1:mxarr(ibl),1:myarr(ibl)),z99) 
        ENDDO parallel_block_loop2 
        DEALLOCATE(data2,lfs,lfsh,lfsv,lfsi) 
      ENDIF 
 
      RETURN 
      END SUBROUTINE exch_lagr_2D 
!-----------------------------------------------------------------------
!     subprogram 4. exch_lagr_3D(laq1,laq2,ibl_id)                      
!     Place the data from a block into a 1D array, communicate it if    
!     necessary, and then unpack it into a global array                 
!-----------------------------------------------------------------------
      SUBROUTINE exch_lagr_3D(laq1,laq2,ibl_id) 
      IMPLICIT NONE 
      INTEGER(i4) :: ibl_id 
      TYPE(lagr_quad_type), INTENT(IN) :: laq1 
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq2 
      INTEGER(i4) :: mxl,myl,ibl,iblar(nprocs) 
      INTEGER(i4) :: status(mpi_status_size) 
      INTEGER(i4) :: i,xn,yn,qty,n_side,n_int,nfour,length,ierr 
      INTEGER(i4) :: myzfs,myzfsh,myzfsv,myzfsi,mylength,modelo,modehi 
      INTEGER(i4), DIMENSION(nprocs) :: zfs,zfsh,zfsv,zfsi,nfar,ml,mh 
      COMPLEX(r8), ALLOCATABLE :: data1(:),data2(:,:) 
      LOGICAL, ALLOCATABLE :: lfs(:,:,:,:) 
      LOGICAL, DIMENSION(:,:,:,:,:), ALLOCATABLE :: lfsh,lfsv,lfsi 
!-----------------------------------------------------------------------
!     initialize global lagrange quad variables                         
!-----------------------------------------------------------------------
      IF (laq1%nqty /= laq2%nqty)                                       &
     &  CALL nim_stop('Incompatible array sizes')
      laq2%title = laq1%title;    laq2%name = laq1%name 
      qty    = laq2%nqty
      n_side = laq2%n_side
      n_int  = laq2%n_int
      nfour=nmodes_total
      laq2%nfour=nfour 
 
      modelo=ilayer*nmodes+1 
      modehi=modelo+nmodes-1 
      IF (nprocs > 1) THEN 
        CALL mpi_allgather(ibl_id,1,mpi_nim_int,iblar,                  &
     &                   1,mpi_nim_int,comm_nimrod,ierr)                
        CALL mpi_allgather(laq1%nfour,1,mpi_nim_int,nfar,               &
     &                   1,mpi_nim_int,comm_nimrod,ierr)                
        CALL mpi_allgather(modelo,1,mpi_nim_int,ml,                     &
     &                   1,mpi_nim_int,comm_nimrod,ierr)                
        CALL mpi_allgather(modehi,1,mpi_nim_int,mh,                     &
     &                   1,mpi_nim_int,comm_nimrod,ierr)                
      ELSE 
        iblar=ibl_id 
        nfar=laq1%nfour 
        ml=modelo 
        mh=modehi 
      ENDIF 
      mxl=MAXVAL(mxarr);  myl=MAXVAL(myarr) 
      ALLOCATE(lfs(1:qty,0:mxl,0:myl,1:nfour)) 
                                                ! Filter for PACK/UNPACK
      lfs=.true. 
!-----------------------------------------------------------------------
!     Exchange the data from block_decomposed data to single block data 
!-----------------------------------------------------------------------
      ! poly_degree=1                                                   
      IF (n_side==0) THEN 
        myzfs=SIZE(laq1%fs) 
        IF (nprocs > 1) THEN 
          CALL mpi_allgather(myzfs,1,mpi_nim_int,zfs,                   &
     &                     1,mpi_nim_int,comm_nimrod,ierr)              
        ELSE 
          zfs=myzfs 
        ENDIF 
        length=MAXVAL(zfs) 
        ALLOCATE(data1(length)); data1=0 
        data1(1:myzfs)=PACK(laq1%fs(:,:,:,:),.true.) 
        ALLOCATE(data2(length,nprocs)) 
        IF (nprocs > 1) THEN 
          CALL mpi_allgather(data1,length,mpi_nim_comp,data2,           &
     &                    length,mpi_nim_comp,comm_nimrod,ierr)         
        ELSE 
          data2(:,1)=data1 
        ENDIF 
        DEALLOCATE(data1) 
        parallel_block_loop1: DO i=1,nprocs 
          ibl=iblar(i) 
          xn = x0(ibl) 
          yn = y0(ibl) 
          laq2%fs(1:qty,xn:xn+mxarr(ibl),yn:yn+myarr(ibl),ml(i):mh(i))= &
     &     UNPACK(data2(1:zfs(i),i),lfs(1:qty,0:mxarr(ibl),0:myarr(ibl),&
     &            1:nfar(i)),zc99)                                      
        ENDDO parallel_block_loop1 
        DEALLOCATE(data2,lfs) 
                                                                        
      ! poly_degree=1                                                   
      ELSE 
        ALLOCATE(lfsh(1:qty,1:n_side,1:mxl,0:myl,1:nfour)) 
        ALLOCATE(lfsv(1:qty,1:n_side,0:mxl,1:myl,1:nfour)) 
        ALLOCATE(lfsi(1:qty,1:n_int,1:mxl,1:myl,1:nfour)) 
        lfsh=.true. 
        lfsv=.true. 
        lfsi=.true. 
        myzfs=SIZE(laq1%fs) 
        myzfsh=SIZE(laq1%fsh) 
        myzfsv=SIZE(laq1%fsv) 
        myzfsi=SIZE(laq1%fsi) 
        mylength=myzfs+myzfsh+myzfsv+myzfsi 
        IF (nprocs > 1) THEN 
          CALL mpi_allgather(myzfs,1,mpi_nim_int,zfs,                   &
     &                     1,mpi_nim_int,comm_nimrod,ierr)              
          CALL mpi_allgather(myzfsh,1,mpi_nim_int,zfsh,                 &
     &                     1,mpi_nim_int,comm_nimrod,ierr)              
          CALL mpi_allgather(myzfsv,1,mpi_nim_int,zfsv,                 &
     &                     1,mpi_nim_int,comm_nimrod,ierr)              
          CALL mpi_allgather(myzfsi,1,mpi_nim_int,zfsi,                 &
     &                     1,mpi_nim_int,comm_nimrod,ierr)              
          CALL mpi_allreduce(mylength,length,1,mpi_nim_int,             &
     &                     MPI_MAX,comm_nimrod,ierr)                    
        ELSE 
          zfs=myzfs 
          zfsh=myzfsh 
          zfsv=myzfsv 
          zfsi=myzfsi 
          length=mylength 
        ENDIF 
        ALLOCATE(data1(length)); data1=0 
        data1(1:myzfs)=PACK(laq1%fs(:,:,:,:),.true.) 
        data1(myzfs+1:myzfs+myzfsh)=PACK(laq1%fsh(:,:,:,:,:),.true.) 
        data1(myzfs+myzfsh+1:myzfs+myzfsh+myzfsv)=                      &
     &        PACK(laq1%fsv(:,:,:,:,:),.true.)                          
        data1(myzfs+myzfsh+myzfsv+1:myzfs+myzfsh+myzfsv+myzfsi)=        &
     &        PACK(laq1%fsi(:,:,:,:,:),.true.)                          
        ALLOCATE(data2(length,nprocs)) 
        IF (nprocs > 1) THEN 
          CALL mpi_allgather(data1,length,mpi_nim_comp,data2,           &
     &                     length,mpi_nim_comp,comm_nimrod,ierr)        
        ELSE 
          data2(:,1)=data1 
        ENDIF 
        DEALLOCATE(data1) 
        parallel_block_loop2: DO i=1,nprocs 
          ibl=iblar(i) 
          xn = x0(ibl);  yn = y0(ibl) 
          laq2%fs(1:qty,xn:xn+mxarr(ibl),yn:yn+myarr(ibl),ml(i):mh(i))= &
     &       UNPACK(data2(1:zfs(i),i),                                  &
     &       lfs(1:qty,0:mxarr(ibl),0:myarr(ibl),1:nfar(i)),zc99)       
          laq2%fsh(1:qty,1:n_side,xn+1:xn+mxarr(ibl),yn:yn+myarr(ibl)   &
     &             ,ml(i):mh(i))=                                       &
     &       UNPACK(data2(zfs(i)+1:zfs(i)+zfsh(i),i),                   &
     &              lfsh(1:qty,1:n_side,1:mxarr(ibl),0:myarr(ibl),      &
     &                   1:nfar(i)),zc99)                               
          laq2%fsv(1:qty,1:n_side,xn:xn+mxarr(ibl),yn+1:yn+myarr(ibl)   &
     &             ,ml(i):mh(i))=                                       &
     &       UNPACK(data2(zfs(i)+zfsh(i)+1:zfs(i)+zfsh(i)+zfsv(i),i),   &
     &              lfsv(1:qty,1:n_side,0:mxarr(ibl),1:myarr(ibl),      &
     &                   1:nfar(i)),zc99)                               
          laq2%fsi(1:qty,1:n_int,xn+1:xn+mxarr(ibl),yn+1:yn+myarr(ibl)  &
     &             ,ml(i):mh(i))=                                       &
     &       UNPACK(data2(zfs(i)+zfsh(i)+zfsv(i)+1:                     &
     &                    zfs(i)+zfsh(i)+zfsv(i)+zfsi(i),i),            &
     &              lfsi(1:qty,1:n_int,1:mxarr(ibl),1:myarr(ibl),       &
     &                   1:nfar(i)),zc99)                               
        ENDDO parallel_block_loop2 
        DEALLOCATE(data2,lfs,lfsh,lfsv,lfsi) 
      ENDIF 
 
      RETURN 
      END SUBROUTINE exch_lagr_3D 
!-----------------------------------------------------------------------
!     subprogram 4. allreduce_lagr_2D(laq2)                      
!     allreduce (MPI_SUM) global 2D rblock data.
!-----------------------------------------------------------------------
      SUBROUTINE allreduce_lagr_2D(laq2) 
      IMPLICIT NONE 
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq2 
      INTEGER(i4) :: status(mpi_status_size) 
      INTEGER(i4) :: qty,n_side,n_int,length,ierr 
      INTEGER(i4) :: nfs,nfsh,nfsv,nfsi
      REAL(r8), ALLOCATABLE :: data1(:),data2(:) 
      LOGICAL, ALLOCATABLE :: lfs(:,:,:) 
      LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: lfsh,lfsv,lfsi 
!-----------------------------------------------------------------------
!     abbreviate local (to this subroutine) variables
!-----------------------------------------------------------------------
      qty    = laq2%nqty
      n_side = laq2%n_side
      n_int  = laq2%n_int
      ! poly_degree=1                                                   
      IF (n_side==0) THEN 
        length=SIZE(laq2%fs)
        ALLOCATE(data1(length)); data1=0
        data1=PACK(laq2%fs(:,:,:),.true.)
        ALLOCATE(data2(length)); data2=0
        CALL mpi_allreduce(data1,data2,length,mpi_nim_real,MPI_SUM,     &
     &                                                 comm_nimrod,ierr)
        DEALLOCATE(data1)
        ALLOCATE(lfs(1:qty,0:mx,0:my)) 
        lfs=.true.  ! Filter for PACK/UNPACK
        laq2%fs(1:qty,0:mx,0:my)=UNPACK(data2,lfs(1:qty,0:mx,0:my),z99)
        DEALLOCATE(data2,lfs)
      ! poly_degree=1                                                   
      ELSE 
        nfs=SIZE(laq2%fs) 
        nfsh=SIZE(laq2%fsh) ; nfsh=nfsh+nfs
        nfsv=SIZE(laq2%fsv) ; nfsv=nfsv+nfsh
        nfsi=SIZE(laq2%fsi) ; nfsi=nfsi+nfsv
        length = nfsi
        ALLOCATE(data1(length)); data1=0 
        data1(1:nfs)      =PACK(laq2%fs(:,:,:),.true.) 
        data1(nfs+1:nfsh) =PACK(laq2%fsh(:,:,:,:),.true.) 
        data1(nfsh+1:nfsv)=PACK(laq2%fsv(:,:,:,:),.true.)
        data1(nfsv+1:nfsi)=PACK(laq2%fsi(:,:,:,:),.true.)
        ALLOCATE(data2(length)); data2=0
        CALL mpi_allreduce(data1,data2,length,mpi_nim_real,MPI_SUM,     &
     &                                                 comm_nimrod,ierr)
        DEALLOCATE(data1)
        ALLOCATE(lfs(1:qty,0:mx,0:my)) 
        ALLOCATE(lfsh(1:qty,1:n_side,1:mx,0:my)) 
        ALLOCATE(lfsv(1:qty,1:n_side,0:mx,1:my)) 
        ALLOCATE(lfsi(1:qty,1:n_int,1:mx,1:my)) 
        lfs=.true. ; lfsh=.true. ; lfsv=.true. ; lfsi=.true. 
        laq2%fs(1:qty,0:mx,0:my)=                                       &
     &     UNPACK(data2(1:nfs),lfs(1:qty,0:mx,0:my),z99)
        laq2%fsh(1:qty,1:n_side,1:mx,0:my)=                             &
     &     UNPACK(data2(nfs+1:nfsh),                                    &
     &                         lfsh(1:qty,1:n_side,1:mx,0:my),z99)
        laq2%fsv(1:qty,1:n_side,0:mx,1:my)=                             &
     &     UNPACK(data2(nfsh+1:nfsv),                                   &
     &                         lfsv(1:qty,1:n_side,0:mx,1:my),z99)
        laq2%fsi(1:qty,1:n_int,1:mx,1:my)=                              &
     &     UNPACK(data2(nfsv+1:nfsi),                                   &
     &                         lfsi(1:qty,1:n_int,1:mx,1:my),z99)
        DEALLOCATE(data2,lfs,lfsh,lfsv,lfsi)
      ENDIF
 
      RETURN 
      END SUBROUTINE allreduce_lagr_2D 
!-----------------------------------------------------------------------
!     subprogram 4. allreduce_lagr_3D(laq2)                      
!     allreduce (MPI_SUM) global 3D rblock data.
!-----------------------------------------------------------------------
      SUBROUTINE allreduce_lagr_3D(laq2) 
      IMPLICIT NONE 
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq2 
      INTEGER(i4) :: status(mpi_status_size) 
      INTEGER(i4) :: i,j,qty,n_side,n_int,nfour,length,ierr 
      INTEGER(i4) :: nfs,nfsh,nfsv,nfsi
      INTEGER(i4), DIMENSION(nprocs) :: ml,mh 
      COMPLEX(r8), ALLOCATABLE :: data1(:),data2(:) 
      LOGICAL, ALLOCATABLE :: lfs(:,:,:,:) 
      LOGICAL, DIMENSION(:,:,:,:,:), ALLOCATABLE :: lfsh,lfsv,lfsi 
!-----------------------------------------------------------------------
!     abbreviate local (to this subroutine) variables
!-----------------------------------------------------------------------
      qty    = laq2%nqty
      n_side = laq2%n_side
      n_int  = laq2%n_int
      nfour=nmodes_total
      ! poly_degree=1                                                   
      IF (n_side==0) THEN 
        length=SIZE(laq2%fs)
        ALLOCATE(data1(length)); data1=0
        data1=PACK(laq2%fs(:,:,:,:),.true.)
        ALLOCATE(data2(length)); data2=0
        CALL mpi_allreduce(data1,data2,length,mpi_nim_comp,MPI_SUM,     &
     &                                                 comm_nimrod,ierr)
        DEALLOCATE(data1)
        ALLOCATE(lfs(1:qty,0:mx,0:my,1:nfour)) 
        lfs=.true.  ! Filter for PACK/UNPACK
        laq2%fs(1:qty,0:mx,0:my,1:nfour)=                               &
     &     UNPACK(data2,lfs(1:qty,0:mx,0:my,1:nfour),zc99)
        DEALLOCATE(data2,lfs)
      ! poly_degree=1                                                   
      ELSE 
        nfs=SIZE(laq2%fs) 
        nfsh=SIZE(laq2%fsh) ; nfsh=nfsh+nfs
        nfsv=SIZE(laq2%fsv) ; nfsv=nfsv+nfsh
        nfsi=SIZE(laq2%fsi) ; nfsi=nfsi+nfsv
        length = nfsi
        ALLOCATE(data1(length)); data1=0 
        data1(1:nfs)      =PACK(laq2%fs(:,:,:,:),.true.) 
        data1(nfs+1:nfsh) =PACK(laq2%fsh(:,:,:,:,:),.true.) 
        data1(nfsh+1:nfsv)=PACK(laq2%fsv(:,:,:,:,:),.true.)
        data1(nfsv+1:nfsi)=PACK(laq2%fsi(:,:,:,:,:),.true.)
        ALLOCATE(data2(length)); data2=0
        CALL mpi_allreduce(data1,data2,length,mpi_nim_comp,MPI_SUM,     &
     &                                        comm_nimrod,ierr)
        DEALLOCATE(data1)
        ALLOCATE(lfs(1:qty,0:mx,0:my,1:nfour)) 
        ALLOCATE(lfsh(1:qty,1:n_side,1:mx,0:my,1:nfour)) 
        ALLOCATE(lfsv(1:qty,1:n_side,0:mx,1:my,1:nfour)) 
        ALLOCATE(lfsi(1:qty,1:n_int,1:mx,1:my,1:nfour)) 
        lfs=.true. ; lfsh=.true. ; lfsv=.true. ; lfsi=.true. 
        laq2%fs(1:qty,0:mx,0:my,1:nfour)=                               &
     &     UNPACK(data2(1:nfs),lfs(1:qty,0:mx,0:my,1:nfour),zc99)
        laq2%fsh(1:qty,1:n_side,1:mx,0:my,1:nfour)=                     &
     &     UNPACK(data2(nfs+1:nfsh),                                    &
     &               lfsh(1:qty,1:n_side,1:mx,0:my,1:nfour),zc99)
        laq2%fsv(1:qty,1:n_side,0:mx,1:my,1:nfour)=                     &
     &     UNPACK(data2(nfsh+1:nfsv),                                   &
     &               lfsv(1:qty,1:n_side,0:mx,1:my,1:nfour),zc99)
        laq2%fsi(1:qty,1:n_int,1:mx,1:my,1:nfour)=                      &
     &     UNPACK(data2(nfsv+1:nfsi),                                   &
     &                lfsi(1:qty,1:n_int,1:mx,1:my,1:nfour),zc99)
        DEALLOCATE(data2,lfs,lfsh,lfsv,lfsi)
      ENDIF 
 
      RETURN 
      END SUBROUTINE allreduce_lagr_3D 
!-----------------------------------------------------------------------
!     subprogram 4. reduce_lagr_3D(laq1,laq2,ibl_id)                      
!     reduce (MPI_SUM) global rblock data in laq2 to root node then 
!     communicate local rblocks back.
!-----------------------------------------------------------------------
      SUBROUTINE reduce_lagr_3D(laq1,laq2,ibl_id) 
      IMPLICIT NONE 
      INTEGER(i4) :: ibl_id 
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq1 
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq2 
      INTEGER(i4) :: ibl,iblar(nprocs) 
      INTEGER(i4) :: status(mpi_status_size) 
      INTEGER(i4) :: i,j,xn,yn,qty,n_side,n_int,nfour,length,ierr 
      INTEGER(i4) :: nfs,nfsh,nfsv,nfsi,modelo,modehi 
      INTEGER(i4), DIMENSION(nprocs) :: ml,mh 
      COMPLEX(r8), ALLOCATABLE :: data1(:),data2(:) 
      COMPLEX(r8), ALLOCATABLE :: fst(:,:,:,:) 
      COMPLEX(r8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: fsth,fstv,fsti
      LOGICAL, ALLOCATABLE :: lfs(:,:,:,:) 
      LOGICAL, DIMENSION(:,:,:,:,:), ALLOCATABLE :: lfsh,lfsv,lfsi 
!-----------------------------------------------------------------------
!     initialize global lagrange quad variables                         
!-----------------------------------------------------------------------
      laq2%title = laq1%title;    laq2%name = laq1%name 
      qty    = laq1%nqty;         laq2%nqty   = qty 
      n_side = laq1%n_side;       laq2%n_side = n_side 
      n_int  = laq1%n_int;        laq2%n_int  = n_int 
      nfour=nmodes_total
      laq2%nfour=nfour 
 
      modelo=ilayer*nmodes+1 
      modehi=modelo+nmodes-1 
      CALL mpi_allgather(ibl_id,1,mpi_nim_int,iblar,                    &
     &                 1,mpi_nim_int,comm_nimrod,ierr)                  
      CALL mpi_allgather(modelo,1,mpi_nim_int,ml,                       &
     &                 1,mpi_nim_int,comm_nimrod,ierr)                  
      CALL mpi_allgather(modehi,1,mpi_nim_int,mh,                       &
     &                 1,mpi_nim_int,comm_nimrod,ierr)                
!-----------------------------------------------------------------------
!     Exchange the data from block_decomposed data to single block data 
!-----------------------------------------------------------------------
      ! poly_degree=1                                                   
      IF (n_side==0) THEN 
        length=SIZE(laq2%fs)
        ALLOCATE(data1(length)); data1=0
        data1=PACK(laq2%fs(:,:,:,:),.true.)
        IF (node==0) THEN
          ALLOCATE(data2(length))
          data2=0
        ENDIF
        CALL mpi_reduce(data1,data2,length,mpi_nim_comp,MPI_SUM,0,      &
     &                                        comm_nimrod,ierr)
        DEALLOCATE(data1)
        IF (node==0) THEN
          ALLOCATE(lfs(1:qty,0:mx,0:my,1:nfour)) 
          lfs=.true.  ! Filter for PACK/UNPACK
          laq2%fs(1:qty,0:mx,0:my,1:nfour)=                             &
     &     UNPACK(data2,lfs(1:qty,0:mx,0:my,1:nfour),zc99)
          DEALLOCATE(data2,lfs)
          DO i=1,nprocs
            ibl=iblar(i) 
            xn = x0(ibl) 
            yn = y0(ibl) 
            ALLOCATE(fst(1:qty,0:mxarr(ibl),0:myarr(ibl),1:nmodes))
            fst=laq2%fs(1:qty,xn:xn+mxarr(ibl),yn:yn+myarr(ibl),        &
     &                                              ml(i):mh(i))
            IF (i==1) THEN
               laq1%fs(1:qty,0:mxarr(1),0:myarr(1),1:nmodes)=fst
            ELSE
               CALL mpi_send(fst,SIZE(fst),mpi_nim_comp,i-1_i4,0,       &
     &                                             comm_nimrod,ierr)
            ENDIF
            DEALLOCATE(fst)
          ENDDO
        ELSE
          CALL mpi_recv(laq1%fs,SIZE(laq1%fs),mpi_nim_comp,0,0,         &
     &         comm_nimrod,status,ierr)
        ENDIF
      ! poly_degree=1                                                   
      ELSE 
        nfs=SIZE(laq2%fs) 
        nfsh=SIZE(laq2%fsh) ; nfsh=nfsh+nfs
        nfsv=SIZE(laq2%fsv) ; nfsv=nfsv+nfsh
        nfsi=SIZE(laq2%fsi) ; nfsi=nfsi+nfsv
        length = nfsi
        ALLOCATE(data1(length)); data1=0 
        data1(1:nfs)      =PACK(laq2%fs(:,:,:,:),.true.) 
        data1(nfs+1:nfsh) =PACK(laq2%fsh(:,:,:,:,:),.true.) 
        data1(nfsh+1:nfsv)=PACK(laq2%fsv(:,:,:,:,:),.true.)
        data1(nfsv+1:nfsi)=PACK(laq2%fsi(:,:,:,:,:),.true.)
        ALLOCATE(data2(length)); data2=0
        CALL mpi_reduce(data1,data2,length,mpi_nim_comp,MPI_SUM,0,      &
     &                                        comm_nimrod,ierr)
        DEALLOCATE(data1)
        IF (node==0) THEN
          ALLOCATE(lfs(1:qty,0:mx,0:my,1:nfour)) 
          ALLOCATE(lfsh(1:qty,1:n_side,1:mx,0:my,1:nfour)) 
          ALLOCATE(lfsv(1:qty,1:n_side,0:mx,1:my,1:nfour)) 
          ALLOCATE(lfsi(1:qty,1:n_int,1:mx,1:my,1:nfour)) 
          lfs=.true. ; lfsh=.true. ; lfsv=.true. ; lfsi=.true. 
          laq2%fs(1:qty,0:mx,0:my,1:nfour)=                             &
     &     UNPACK(data2(1:nfs),lfs(1:qty,0:mx,0:my,1:nfour),zc99)
          laq2%fsh(1:qty,1:n_side,1:mx,0:my,1:nfour)=                   &
     &     UNPACK(data2(nfs+1:nfsh),                                    &
     &               lfsh(1:qty,1:n_side,1:mx,0:my,1:nfour),zc99)
          laq2%fsv(1:qty,1:n_side,0:mx,1:my,1:nfour)=                   &
     &     UNPACK(data2(nfsh+1:nfsv),                                   &
     &               lfsv(1:qty,1:n_side,0:mx,1:my,1:nfour),zc99)
          laq2%fsi(1:qty,1:n_int,1:mx,1:my,1:nfour)=                    &
     &     UNPACK(data2(nfsv+1:nfsi),                                   &
     &                lfsi(1:qty,1:n_int,1:mx,1:my,1:nfour),zc99)
          DEALLOCATE(data2,lfs,lfsh,lfsv,lfsi)
          DO i=1,nprocs
            ibl=iblar(i) 
            xn = x0(ibl) 
            yn = y0(ibl) 
! need to make fsth,fstv and fsti conformal yet!
            ALLOCATE(fst(1:qty,0:mxarr(ibl),0:myarr(ibl),1:nmodes))
            ALLOCATE(fsth(1:qty,1:n_side,1:mxarr(ibl),0:myarr(ibl),1:nmodes))
            ALLOCATE(fstv(1:qty,1:n_side,0:mxarr(ibl),1:myarr(ibl),1:nmodes))
            ALLOCATE(fsti(1:qty,1:n_int,1:mxarr(ibl),1:myarr(ibl),1:nmodes))
            fst=laq2%fs(1:qty,xn:xn+mxarr(ibl),                         &
     &                                   yn:yn+myarr(ibl),ml(i):mh(i))
            fsth=laq2%fsh(1:qty,1:n_side,xn+1:xn+mxarr(ibl),            &
     &                                   yn:yn+myarr(ibl),ml(i):mh(i))
            fstv=laq2%fsv(1:qty,1:n_side,xn:xn+mxarr(ibl),              &
     &                                   yn+1:yn+myarr(ibl),ml(i):mh(i))
            fsti=laq2%fsi(1:qty,1:n_int,xn+1:xn+mxarr(ibl),              &
     &                                  yn+1:yn+myarr(ibl),ml(i):mh(i))
            IF (i==1) THEN
              laq1%fs(1:qty,0:mxarr(1),0:myarr(1),1:nmodes)=fst
              laq1%fsh(1:qty,1:n_side,1:mxarr(1),0:myarr(1),1:nmodes)=fsth
              laq1%fsv(1:qty,1:n_side,0:mxarr(1),1:myarr(1),1:nmodes)=fstv
              laq1%fsi(1:qty,1:n_int,1:mxarr(1),1:myarr(1),1:nmodes)=fsti
            ELSE
              CALL mpi_send(fst,SIZE(fst),                              &
     &                      mpi_nim_comp,i-1_i4,0,comm_nimrod,ierr)
              CALL mpi_send(fsth,SIZE(fsth),                            &
     &                      mpi_nim_comp,i-1_i4,0,comm_nimrod,ierr)
              CALL mpi_send(fstv,SIZE(fstv),                            &
     &                      mpi_nim_comp,i-1_i4,0,comm_nimrod,ierr)
              CALL mpi_send(fsti,SIZE(fsti),                            &
     &                      mpi_nim_comp,i-1_i4,0,comm_nimrod,ierr)
            ENDIF
            DEALLOCATE(fst,fsth,fstv,fsti)
          ENDDO
        ELSE
          CALL mpi_recv(laq1%fs,SIZE(laq1%fs),mpi_nim_comp,0,0,         &
     &         comm_nimrod,status,ierr)
          CALL mpi_recv(laq1%fsh,SIZE(laq1%fsh),mpi_nim_comp,0,0,       &
     &         comm_nimrod,status,ierr)
          CALL mpi_recv(laq1%fsv,SIZE(laq1%fsv),mpi_nim_comp,0,0,       &
     &         comm_nimrod,status,ierr)
          CALL mpi_recv(laq1%fsi,SIZE(laq1%fsi),mpi_nim_comp,0,0,       &
     &         comm_nimrod,status,ierr)
        ENDIF
      ENDIF 
 
      RETURN 
      END SUBROUTINE reduce_lagr_3D 
!-----------------------------------------------------------------------
      END MODULE exchange_mod 
