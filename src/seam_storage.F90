!-----------------------------------------------------------------------
!     $Id: seam_storage.F90 7005 2020-01-09 21:19:16Z held $
!     module containing the seam structures used for communications     
!     among grid blocks.                                                
!-----------------------------------------------------------------------                                  
!     2. dump_read_seam.                                                                  
!     4. h5_dump_seam.                                    
!     5. h5_dump_sdum.                                    
!     6. dump_read_seam.   
!-----------------------------------------------------------------------
#include "config.f"
      MODULE seam_storage_mod 
      USE edge_type_mod 
      IMPLICIT NONE 
                                                                        
      TYPE(edge_type), POINTER :: seam0 
      TYPE(edge_type), DIMENSION(:), POINTER :: seam 
      INTEGER(i4), DIMENSION(:), POINTER :: exblock_list 
      INTEGER(i4) :: max_imags 
                                                                       
      CONTAINS 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     subprogram 2. dump_read_seam.                                     
!     reads seam data to dump file.                                     
!-----------------------------------------------------------------------
      SUBROUTINE dump_read_seam(seam) 
      USE io 
                                                                        
      TYPE (edge_type), INTENT(INOUT) :: seam 
                                                                        
      INTEGER(i4) :: i,nv,np,ip,iq,ivert
      REAL(r8) :: iread 
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ia1read 
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: ia2read 
!-----------------------------------------------------------------------
!     read block descriptors.                                           
!-----------------------------------------------------------------------
      seam%name='seam' 
      READ(rstrt_unit) iread 
      seam%id=NINT(iread) 
      READ(rstrt_unit) iread 
      seam%nvert=NINT(iread) 
!-----------------------------------------------------------------------
!     read seam data.                                                   
!-----------------------------------------------------------------------
      ALLOCATE(seam%vertex(seam%nvert)) 
      ALLOCATE(ia1read(2)) 
      seam%npt=0
      DO i=1,seam%nvert 
         READ(rstrt_unit) iread 
         np=NINT(iread) 
         seam%npt=seam%npt+np
         ALLOCATE(seam%vertex(i)%ptr(2,np)) 
         ALLOCATE(ia2read(2,np)) 
         READ(rstrt_unit) ia2read 
         seam%vertex(i)%ptr=NINT(ia2read) 
         DEALLOCATE(ia2read) 
         IF (seam%id>0) THEN 
           READ(rstrt_unit) ia1read 
           seam%vertex(i)%intxy=NINT(ia1read) 
         ENDIF 
      ENDDO 
      DEALLOCATE(ia1read) 
!-----------------------------------------------------------------------
!     read external corner data.                                        
!-----------------------------------------------------------------------
      IF(seam%id==0.AND.seam%nvert>0)THEN 
         ALLOCATE(seam%excorner(seam%nvert)) 
         DO i=1,seam%nvert 
            READ(rstrt_unit) iread 
            seam%excorner(i)=(NINT(iread)>0) 
         ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     nullify pointers
!-----------------------------------------------------------------------
      DO ivert=1,seam%nvert
        NULLIFY(seam%vertex(ivert)%ptr2)
        NULLIFY(seam%vertex(ivert)%order)
        NULLIFY(seam%vertex(ivert)%seam_in)
        NULLIFY(seam%vertex(ivert)%seam_out)
        NULLIFY(seam%vertex(ivert)%seam_hold)
        NULLIFY(seam%vertex(ivert)%seam_save)
        NULLIFY(seam%vertex(ivert)%seam_cin)
        NULLIFY(seam%vertex(ivert)%seam_cout)
        NULLIFY(seam%vertex(ivert)%seam_chold)
        NULLIFY(seam%vertex(ivert)%seam_csave)
      ENDDO
      NULLIFY(seam%segment)
      NULLIFY(seam%expoint)
      NULLIFY(seam%r0point)
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE dump_read_seam 
#ifdef HAVE_FC_HDF5
!-----------------------------------------------------------------------
!     subprogram 4. h5_dump_seam.                                    
!     writes seam data to dump file.                                    
!-----------------------------------------------------------------------
      SUBROUTINE h5_dump_seam(seam,gid)
      USE io 
                                                                        
      TYPE (edge_type), INTENT(IN) :: seam 
      INTEGER(HID_T), INTENT(IN) :: gid
                                                                        
      INTEGER(i4) :: iv,nv,np,ipt
      INTEGER(i4), ALLOCATABLE :: exc(:),bptr(:,:),bxy(:,:)
      CHARACTER(64) :: group_name,seam_name
      INTEGER(HID_T) :: sid
!-----------------------------------------------------------------------
!     setup seam group.
!-----------------------------------------------------------------------
      WRITE(seam_name,fmt='(i4.4)') seam%id
      CALL make_group(gid,TRIM(seam_name),sid,h5in,h5err)
!-----------------------------------------------------------------------
!     write block descriptors.                                          
!-----------------------------------------------------------------------
      CALL write_attribute(sid,"id",seam%id,h5in,h5err)
      CALL write_attribute(sid,"nvert",seam%nvert,h5in,h5err)
!-----------------------------------------------------------------------
!     write seam data, combine to write larger arrays.
!-----------------------------------------------------------------------
      ALLOCATE(exc(seam%nvert))
      DO iv=1,seam%nvert 
        np=SIZE(seam%vertex(iv)%ptr,2)
        exc(iv)=np
      ENDDO
      CALL dump_h5(sid,"np",exc,h5in,h5err)
      ALLOCATE(bptr(2,seam%npt))
      IF (seam%id>0) ALLOCATE(bxy(2,seam%nvert))
      ipt=1_i4
      DO iv=1,seam%nvert 
        np=SIZE(seam%vertex(iv)%ptr,2)
        bptr(:,ipt:ipt+np-1)=seam%vertex(iv)%ptr
        ipt=ipt+np
        IF (seam%id>0) bxy(:,iv)=seam%vertex(iv)%intxy
      ENDDO 
      CALL dump_h5(sid,"vertex",bptr,h5in,h5err)
      DEALLOCATE(bptr)
      IF (seam%id>0) THEN
        CALL dump_h5(sid,"intxy",bxy,h5in,h5err)
        DEALLOCATE(bxy)
      ENDIF
!-----------------------------------------------------------------------
!     write external corner data.                                       
!-----------------------------------------------------------------------
      IF(seam%id==0.AND.seam%nvert>0) THEN 
        DO iv=1,seam%nvert 
          IF (seam%excorner(iv)) THEN 
            exc(iv) = 1_i4
          ELSE 
            exc(iv) = 0_i4
          ENDIF 
        ENDDO 
        CALL dump_h5(sid,"excorner",exc,h5in,h5err)
      ENDIF 
      DEALLOCATE(exc)
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      CALL close_group(TRIM(seam_name),sid,h5err)
      RETURN 
      END SUBROUTINE h5_dump_seam
!-----------------------------------------------------------------------
!     subprogram 5. h5_dump_sdum.                                    
!     dummy write seam data to dump file.                                    
!-----------------------------------------------------------------------
      SUBROUTINE h5_dump_sdum(id,gid)
      USE io 
      USE pardata, ONLY: block_sizes
                                                                        
      INTEGER(i4), INTENT(IN) :: id
      INTEGER(HID_T), INTENT(IN) :: gid
                                                                        
      INTEGER(i4) :: iv
      INTEGER(i4) :: dims1(1), dims2(2)
      CHARACTER(64) :: group_name,seam_name
      INTEGER(HID_T) :: sid
!-----------------------------------------------------------------------
!     setup seam group.
!-----------------------------------------------------------------------
      WRITE(seam_name,fmt='(i4.4)') id
      CALL make_group(gid,TRIM(seam_name),sid,h5in,h5err)
!-----------------------------------------------------------------------
!     write block descriptors.                                          
!-----------------------------------------------------------------------
      CALL write_attribute(sid,"id",id,h5in,h5err)
      CALL write_attribute(sid,"nvert",block_sizes(1,id),h5in,h5err)
!-----------------------------------------------------------------------
!     write seam data, combine to write larger arrays.
!-----------------------------------------------------------------------
      dims1(1)=block_sizes(1,id)
      CALL dump_h5_int_dum(sid,"np",dims1,h5in,h5err)
      dims2(1)=2; dims2(2)=block_sizes(2,id)
      CALL dump_h5_int_dum(sid,"vertex",dims2,h5in,h5err)
      IF (id>0) THEN
        dims2(2)=block_sizes(1,id)
        CALL dump_h5_int_dum(sid,"intxy",dims2,h5in,h5err)
      ENDIF
!-----------------------------------------------------------------------
!     write external corner data.                                       
!-----------------------------------------------------------------------
      IF(id==0.AND.block_sizes(1,id)>0) THEN 
        CALL dump_h5_int_dum(sid,"excorner",dims1,h5in,h5err)
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      CALL close_group(TRIM(seam_name),sid,h5err)
      RETURN 
      END SUBROUTINE h5_dump_sdum 
!-----------------------------------------------------------------------
!     subprogram 6. dump_read_seam.                                     
!     reads seam data to dump file.                                     
!-----------------------------------------------------------------------
      SUBROUTINE h5_read_seam(seam,id,gid) 
      USE io 
                                                                        
      TYPE (edge_type), INTENT(INOUT) :: seam 
      INTEGER(i4), INTENT(IN) :: id
      INTEGER(HID_T), INTENT(IN) :: gid
       
      INTEGER(i4) :: iv,nv,np,ipt,ivert
      INTEGER(i4), ALLOCATABLE :: exc(:),bptr(:,:),bxy(:,:)
      CHARACTER(64) :: group_name,seam_name
      INTEGER(HID_T) :: sid
!-----------------------------------------------------------------------
!     open seam group.
!-----------------------------------------------------------------------
      WRITE(seam_name,fmt='(i4.4)') id
      CALL open_group(gid,TRIM(seam_name),sid,h5err)
!-----------------------------------------------------------------------
!     read block descriptors.                                           
!-----------------------------------------------------------------------
      seam%name='seam' 
      CALL read_attribute(sid,"id",seam%id,h5in,h5err)
      CALL read_attribute(sid,"nvert",seam%nvert,h5in,h5err)
!-----------------------------------------------------------------------
!     read seam data.                                                   
!-----------------------------------------------------------------------
      ALLOCATE(exc(seam%nvert),seam%vertex(seam%nvert)) 
      CALL read_h5(sid,"np",exc,h5in,h5err)
      seam%npt=SUM(exc)
      ALLOCATE(bptr(2,seam%npt))
      CALL read_h5(sid,"vertex",bptr,h5in,h5err)
      IF (seam%id>0) THEN
        ALLOCATE(bxy(2,seam%nvert))
        CALL read_h5(sid,"intxy",bxy,h5in,h5err)
      ENDIF
      ipt=1_i4
      DO iv=1,seam%nvert 
        np=exc(iv)
        ALLOCATE(seam%vertex(iv)%ptr(2,np))
        seam%vertex(iv)%ptr=bptr(:,ipt:ipt+np-1)
        ipt=ipt+np
        IF (seam%id>0) seam%vertex(iv)%intxy=bxy(:,iv)
      ENDDO 
!-----------------------------------------------------------------------
!     read external corner data.                                        
!-----------------------------------------------------------------------
      IF(seam%id==0.AND.seam%nvert>0)THEN 
        ALLOCATE(seam%excorner(seam%nvert))
        CALL read_h5(sid,"excorner",exc,h5in,h5err)
        DO iv=1,seam%nvert 
          seam%excorner(iv)=(exc(iv)>0) 
        ENDDO 
      ENDIF 
      DEALLOCATE(exc)
!-----------------------------------------------------------------------
!     nullify pointers
!-----------------------------------------------------------------------
      DO ivert=1,seam%nvert
        NULLIFY(seam%vertex(ivert)%ptr2)
        NULLIFY(seam%vertex(ivert)%order)
        NULLIFY(seam%vertex(ivert)%seam_in)
        NULLIFY(seam%vertex(ivert)%seam_out)
        NULLIFY(seam%vertex(ivert)%seam_hold)
        NULLIFY(seam%vertex(ivert)%seam_save)
        NULLIFY(seam%vertex(ivert)%seam_cin)
        NULLIFY(seam%vertex(ivert)%seam_cout)
        NULLIFY(seam%vertex(ivert)%seam_chold)
        NULLIFY(seam%vertex(ivert)%seam_csave)
      ENDDO
      NULLIFY(seam%segment)
      NULLIFY(seam%expoint)
      NULLIFY(seam%r0point)
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      CALL close_group(TRIM(seam_name),sid,h5err)
      RETURN 
      END SUBROUTINE h5_read_seam 
#endif /* HAVE_FC_HDF5 */

!-----------------------------------------------------------------------
!     close module                                                      
!-----------------------------------------------------------------------
      END MODULE seam_storage_mod 
