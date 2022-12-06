!-----------------------------------------------------------------------
!     $Id: rt_swap.f90 5452 2017-08-30 15:55:45Z held $
!     routines to swap EQ fields on dump file reads/writes
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.                                                
!-----------------------------------------------------------------------
!     1. eq_swap
!     2. swap_back
!-----------------------------------------------------------------------
      MODULE rt_swap_mod

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 1. eq_swap                                            
!     when requested, move the equilibrium fields to the n=0 part of the
!     solution vector to use them as initial conditions.                
!-----------------------------------------------------------------------
      SUBROUTINE eq_swap(nmodes,keff) 
      USE local 
      USE fields 
      USE lagr_quad_mod 
      USE input 
      USE physdat 
      IMPLICIT NONE 
                                                                        
      INTEGER(i4), INTENT(IN) :: nmodes 
      REAL(r8), DIMENSION(nmodes), INTENT(IN) :: keff 
                                                                        
      INTEGER(i4) :: ibl,imode,ix,mxb,myb,mb,mv,mc 
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: btmp,jtmp 
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE ::                      &
     &                             bvtmp,jvtmp ,bhtmp,jhtmp ,bitmp,jitmp
      LOGICAL :: n0_found=.false. 
!-----------------------------------------------------------------------
!     Allocate storage space for the variables on all processors.       
!-----------------------------------------------------------------------
      IF (.NOT. ALLOCATED(rb(1)%be_sv%fs)) THEN 
        DO ibl=1,nrbl 
         mxb=rb(ibl)%mx 
         myb=rb(ibl)%my 
                                                                        
         IF (.NOT.ALLOCATED(rb(ibl)%ja%fs)) THEN 
          CALL lagr_quad_alloc(rb(ibl)%ja,mxb,myb,3_i4,nmodes,          &
     &                         poly_degree,'ja',(/'  ja  '/))           
         ENDIF 
         CALL lagr_quad_alloc(rb(ibl)%be_sv,mxb,myb,3_i4,poly_degree) 
         CALL lagr_quad_alloc(rb(ibl)%ja_sv,mxb,myb,3_i4,poly_degree) 
         CALL lagr_quad_alloc(rb(ibl)%ve_sv,mxb,myb,3_i4,poly_degree) 
         CALL lagr_quad_alloc(rb(ibl)%pres_sv,mxb,myb,1_i4,poly_degree) 
         CALL lagr_quad_alloc(rb(ibl)%prese_sv,mxb,myb,1_i4,poly_degree) 
         CALL lagr_quad_alloc(rb(ibl)%tele_sv,mxb,myb,1_i4,poly_degree) 
         CALL lagr_quad_alloc(rb(ibl)%tion_sv,mxb,myb,1_i4,poly_degree) 
        ENDDO 
        DO ibl=nrbl+1,nbl 
          mv=tb(ibl)%mvert 
          mc=tb(ibl)%mcell 
          CALL tri_linear_alloc(tb(ibl)%be_sv,mv,3_i4) 
          CALL tri_linear_alloc(tb(ibl)%ja_sv,mv,3_i4) 
          CALL tri_linear_alloc(tb(ibl)%ve_sv,mv,3_i4) 
          CALL tri_linear_alloc(tb(ibl)%pres_sv,mv,1_i4) 
          CALL tri_linear_alloc(tb(ibl)%prese_sv,mv,1_i4) 
          CALL tri_linear_alloc(tb(ibl)%tele_sv,mv,1_i4) 
          CALL tri_linear_alloc(tb(ibl)%tion_sv,mv,1_i4) 
        ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     Calculate temperatures right away since they aren't in eq vars    
!-----------------------------------------------------------------------
      DO ibl=1,nrbl 
       rb(ibl)%tele_sv%fs=rb(ibl)%prese_eq%fs/rb(ibl)%nd_eq%fs/kboltz 
       rb(ibl)%tion_sv%fs=(rb(ibl)%pres_eq%fs-rb(ibl)%prese_eq%fs)/     &
     &                       (rb(ibl)%nd_eq%fs/zeff*kboltz)             
       IF( poly_degree > 1) THEN 
        rb(ibl)%tele_sv%fsh=                                            &
     &                  rb(ibl)%prese_eq%fsh/(rb(ibl)%nd_eq%fsh*kboltz) 
        rb(ibl)%tele_sv%fsv=                                            &
     &                  rb(ibl)%prese_eq%fsv/(rb(ibl)%nd_eq%fsv*kboltz) 
        rb(ibl)%tele_sv%fsi=                                            &
     &                  rb(ibl)%prese_eq%fsi/(rb(ibl)%nd_eq%fsi*kboltz) 
                                                                        
        rb(ibl)%tion_sv%fsh=                                            &
     &                  (rb(ibl)%pres_eq%fsh-rb(ibl)%prese_eq%fsh)/     &
     &                  (rb(ibl)%nd_eq%fsh/zeff*kboltz)                 
        rb(ibl)%tion_sv%fsv=                                            &
     &                  (rb(ibl)%pres_eq%fsv-rb(ibl)%prese_eq%fsv)/     &
     &                  (rb(ibl)%nd_eq%fsv/zeff*kboltz)                 
        rb(ibl)%tion_sv%fsi=                                            &
     &                  (rb(ibl)%pres_eq%fsi-rb(ibl)%prese_eq%fsi)/     &
     &                  (rb(ibl)%nd_eq%fsi/zeff*kboltz)                 
       ENDIF 
      ENDDO 
      DO ibl=nrbl+1,nbl 
        tb(ibl)%tele_sv%fs=tb(ibl)%prese_eq%fs/(tb(ibl)%nd_eq%fs*kboltz) 
        tb(ibl)%tion_sv%fs=(tb(ibl)%pres_eq%fs-tb(ibl)%prese_eq%fs)/    &
     &                       (tb(ibl)%nd_eq%fs*kboltz)                  
      ENDDO 
!-----------------------------------------------------------------------
!     loop over modes to find the n=0                                   
!-----------------------------------------------------------------------
      DO imode=1,nmodes 
        IF (keff(imode)==0) THEN 
          n0_found=.true. 
          EXIT 
        ENDIF 
      ENDDO 
!-----------------------------------------------------------------------
!     loop over blocks and transfer data to n=0.                        
!-----------------------------------------------------------------------
      IF (n0_found) THEN 
       DO ibl=1,nrbl 
        mxb=SIZE(rb(1)%be_eq%fs,2) 
        myb=SIZE(rb(1)%be_eq%fs,3) 
        ALLOCATE(btmp(3,mxb,myb),jtmp(3,mxb,myb)) 
        btmp=rb(ibl)%be_eq%fs 
        jtmp=rb(ibl)%ja_eq%fs 
        IF (geom=='tor') THEN 
          WHERE(rb(ibl)%rz%fs(1,:,:)==0) 
            btmp(3,:,:)=0. 
          ELSEWHERE 
            btmp(3,:,:)=btmp(3,:,:)/rb(ibl)%rz%fs(1,:,:) 
          ENDWHERE 
          jtmp(3,:,:)=jtmp(3,:,:)*rb(ibl)%rz%fs(1,:,:) 
        ENDIF 
        rb(ibl)%be%fs(:,:,:,imode)=rb(ibl)%be%fs(:,:,:,imode)+btmp 
        rb(ibl)%ja%fs(:,:,:,imode)=rb(ibl)%ja%fs(:,:,:,imode)+jtmp 
        rb(ibl)%ve%fs(:,:,:,imode)=rb(ibl)%ve%fs(:,:,:,imode)+          &
     &                             rb(ibl)%ve_eq%fs                     
        rb(ibl)%pres%fs(:,:,:,imode)=rb(ibl)%pres%fs(:,:,:,imode)+      &
     &                               rb(ibl)%pres_eq%fs                 
        rb(ibl)%prese%fs(:,:,:,imode)=rb(ibl)%prese%fs(:,:,:,imode)+    &
     &                                rb(ibl)%prese_eq%fs               
        rb(ibl)%tele%fs(:,:,:,imode)=rb(ibl)%tele%fs(:,:,:,imode)+      &
     &                             rb(ibl)%tele_sv%fs                   
        rb(ibl)%tion%fs(:,:,:,imode)=rb(ibl)%tion%fs(:,:,:,imode)+      &
     &                             rb(ibl)%tion_sv%fs                   
        IF (poly_degree>1) THEN 
          mb=SIZE(rb(1)%be_eq%fsv,2); ALLOCATE(bvtmp(3,mb,mxb,myb-1)) 
          mb=SIZE(rb(1)%be_eq%fsh,2); ALLOCATE(bhtmp(3,mb,mxb-1,myb)) 
          mb=SIZE(rb(1)%be_eq%fsi,2); ALLOCATE(bitmp(3,mb,mxb-1,myb-1)) 
          mb=SIZE(rb(1)%ja_eq%fsv,2); ALLOCATE(jvtmp(3,mb,mxb,myb-1)) 
          mb=SIZE(rb(1)%ja_eq%fsh,2); ALLOCATE(jhtmp(3,mb,mxb-1,myb)) 
          mb=SIZE(rb(1)%ja_eq%fsi,2); ALLOCATE(jitmp(3,mb,mxb-1,myb-1)) 
          bvtmp=rb(ibl)%be_eq%fsv;    jvtmp=rb(ibl)%ja_eq%fsv 
          bhtmp=rb(ibl)%be_eq%fsh;    jhtmp=rb(ibl)%ja_eq%fsh 
          bitmp=rb(ibl)%be_eq%fsi;    jitmp=rb(ibl)%ja_eq%fsi 
          IF (geom=='tor') THEN 
            WHERE(rb(ibl)%rz%fsv(1,:,:,:)==0) 
              bvtmp(3,:,:,:)=0. 
            ELSEWHERE 
              bvtmp(3,:,:,:)=bvtmp(3,:,:,:)/rb(ibl)%rz%fsv(1,:,:,:) 
            ENDWHERE 
            bhtmp(3,:,:,:)=bhtmp(3,:,:,:)/rb(ibl)%rz%fsh(1,:,:,:) 
            bitmp(3,:,:,:)=bitmp(3,:,:,:)/rb(ibl)%rz%fsi(1,:,:,:) 
            jvtmp(3,:,:,:)=jvtmp(3,:,:,:)*rb(ibl)%rz%fsv(1,:,:,:) 
            jhtmp(3,:,:,:)=jhtmp(3,:,:,:)*rb(ibl)%rz%fsh(1,:,:,:) 
            jitmp(3,:,:,:)=jitmp(3,:,:,:)*rb(ibl)%rz%fsi(1,:,:,:) 
          ENDIF 
          rb(ibl)%be%fsv(:,:,:,:,imode)=rb(ibl)%be%fsv(:,:,:,:,imode)+  &
     &                                  bvtmp                           
          rb(ibl)%be%fsh(:,:,:,:,imode)=rb(ibl)%be%fsh(:,:,:,:,imode)+  &
     &                                  bhtmp                           
          rb(ibl)%be%fsi(:,:,:,:,imode)=rb(ibl)%be%fsi(:,:,:,:,imode)+  &
     &                                  bitmp                           
          rb(ibl)%ja%fsv(:,:,:,:,imode)=rb(ibl)%ja%fsv(:,:,:,:,imode)+  &
     &                                  jvtmp                           
          rb(ibl)%ja%fsh(:,:,:,:,imode)=rb(ibl)%ja%fsh(:,:,:,:,imode)+  &
     &                                  jhtmp                           
          rb(ibl)%ja%fsi(:,:,:,:,imode)=rb(ibl)%ja%fsi(:,:,:,:,imode)+  &
     &                                  jitmp                           
          rb(ibl)%ve%fsh(:,:,:,:,imode)=rb(ibl)%ve%fsh(:,:,:,:,imode)+  &
     &                                  rb(ibl)%ve_eq%fsh               
          rb(ibl)%ve%fsv(:,:,:,:,imode)=rb(ibl)%ve%fsv(:,:,:,:,imode)+  &
     &                                  rb(ibl)%ve_eq%fsv               
          rb(ibl)%ve%fsi(:,:,:,:,imode)=rb(ibl)%ve%fsi(:,:,:,:,imode)+  &
     &                                  rb(ibl)%ve_eq%fsi               
          rb(ibl)%pres%fsh(:,:,:,:,imode)=rb(ibl)%pres_eq%fsh+          &
     &                                   rb(ibl)%pres%fsh(:,:,:,:,imode)
          rb(ibl)%pres%fsv(:,:,:,:,imode)=rb(ibl)%pres_eq%fsv+          &
     &                                   rb(ibl)%pres%fsv(:,:,:,:,imode)
          rb(ibl)%pres%fsi(:,:,:,:,imode)=rb(ibl)%pres_eq%fsi+          &
     &                                   rb(ibl)%pres%fsi(:,:,:,:,imode)
          rb(ibl)%prese%fsh(:,:,:,:,imode)=rb(ibl)%prese_eq%fsh+        &
     &                                  rb(ibl)%prese%fsh(:,:,:,:,imode)
          rb(ibl)%prese%fsv(:,:,:,:,imode)=rb(ibl)%prese_eq%fsv+        &
     &                                  rb(ibl)%prese%fsv(:,:,:,:,imode)
          rb(ibl)%prese%fsi(:,:,:,:,imode)=rb(ibl)%prese_eq%fsi+        &
     &                                  rb(ibl)%prese%fsi(:,:,:,:,imode)
                                                                        
          rb(ibl)%tele%fsh(:,:,:,:,imode)=rb(ibl)%tele_sv%fsh+          &
     &            rb(ibl)%tele%fsh(:,:,:,:,imode)                       
          rb(ibl)%tele%fsv(:,:,:,:,imode)=rb(ibl)%tele_sv%fsv+          &
     &            rb(ibl)%tele%fsv(:,:,:,:,imode)                       
          rb(ibl)%tele%fsi(:,:,:,:,imode)=rb(ibl)%tele_sv%fsi+          &
     &            rb(ibl)%tele%fsi(:,:,:,:,imode)                       
                                                                        
          rb(ibl)%tion%fsh(:,:,:,:,imode)=rb(ibl)%tion_sv%fsh+          &
     &                  rb(ibl)%tion%fsh(:,:,:,:,imode)                 
          rb(ibl)%tion%fsv(:,:,:,:,imode)=rb(ibl)%tion_sv%fsv+          &
     &                  rb(ibl)%tion%fsv(:,:,:,:,imode)                 
          rb(ibl)%tion%fsi(:,:,:,:,imode)=rb(ibl)%tion_sv%fsi+          &
     &                  rb(ibl)%tion%fsi(:,:,:,:,imode)                 
          DEALLOCATE(bvtmp,bitmp,bhtmp,jvtmp,jitmp,jhtmp) 
        ENDIF 
        DEALLOCATE(btmp,jtmp) 
       ENDDO 
!-----------------------------------------------------------------------
!-PRE tblocks, linear only for now.                                     
!-----------------------------------------------------------------------
       DO ibl=nrbl+1,nbl 
        tb(ibl)%be%fs(1:2,:,:,imode)=tb(ibl)%be_eq%fs(1:2,:,:)          &
     &                                   +tb(ibl)%be%fs(1:2,:,:,imode)  
        IF (geom=='tor') THEN 
          DO ix=0,tb(ibl)%mvert 
            IF (tb(ibl)%tgeom%xs(ix)==0) THEN 
              tb(ibl)%be%fs(3,ix,0,imode)=0 
            ELSE 
              tb(ibl)%be%fs(3,ix,0,imode)=tb(ibl)%be_eq%fs(3,ix,0)/     &
     &                                    tb(ibl)%tgeom%xs(ix)          
            ENDIF 
          ENDDO 
        ELSE 
          tb(ibl)%be%fs(3,:,:,imode)=tb(ibl)%be_eq%fs(3,:,:)            &
     &                                   +tb(ibl)%be%fs(3,:,:,imode)    
        ENDIF 
        tb(ibl)%ve%fs(:,:,:,imode)=tb(ibl)%ve_eq%fs                     &
     &                                   +tb(ibl)%ve%fs(:,:,:,imode)    
        tb(ibl)%pres%fs(:,:,:,imode)=tb(ibl)%pres_eq%fs                 &
     &                                   +tb(ibl)%pres%fs(:,:,:,imode)  
        tb(ibl)%prese%fs(:,:,:,imode)=tb(ibl)%prese_eq%fs               &
     &                                   +tb(ibl)%prese%fs(:,:,:,imode) 
        tb(ibl)%tele%fs(:,:,:,imode)=tb(ibl)%tele_sv%fs                 &
     &                                   +tb(ibl)%tele%fs(:,:,:,imode)  
        tb(ibl)%tion%fs(:,:,:,imode)=tb(ibl)%tion_sv%fs(:,:,:)          &
     &                                   +tb(ibl)%tion%fs(:,:,:,imode)  
       ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     Now save the equilibrium arrays and zero them out                 
!     This is done for every processor                                  
!-----------------------------------------------------------------------
      DO ibl=1,nrbl 
        ! Don't do nd_eq.  Screws up continuity options                 
        ! We've already calculate temperatures here                     
        rb(ibl)%be_sv=rb(ibl)%be_eq; rb(ibl)%be_eq=0 
        rb(ibl)%ja_sv=rb(ibl)%ja_eq; rb(ibl)%ja_eq=0 
        rb(ibl)%ve_sv=rb(ibl)%ve_eq; rb(ibl)%ve_eq=0 
        rb(ibl)%pres_sv=rb(ibl)%pres_eq; rb(ibl)%pres_eq=0 
        rb(ibl)%prese_sv=rb(ibl)%prese_eq; rb(ibl)%prese_eq=0 
      ENDDO 
      DO ibl=nrbl+1,nbl 
        tb(ibl)%be_sv=tb(ibl)%be_eq; tb(ibl)%be_eq=0 
        tb(ibl)%ja_sv=tb(ibl)%ja_eq; tb(ibl)%ja_eq=0 
        tb(ibl)%ve_sv=tb(ibl)%ve_eq; tb(ibl)%ve_eq=0 
        tb(ibl)%pres_sv=tb(ibl)%pres_eq; tb(ibl)%pres_eq=0 
        tb(ibl)%prese_sv=tb(ibl)%prese_eq; tb(ibl)%prese_eq=0 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE eq_swap 
!-----------------------------------------------------------------------
!     subprogram 2.  swap_back                                          
!     If we've done an eq_swap, then before writing to the dump files,  
!     go ahead and swap back                                            
!-----------------------------------------------------------------------
      SUBROUTINE swap_back(nmodes,keff) 
      USE local 
      USE fields 
      USE lagr_quad_mod 
      USE input 
      IMPLICIT NONE 
                                                                        
      INTEGER(i4), INTENT(IN) :: nmodes 
      REAL(r8), DIMENSION(nmodes), INTENT(IN) :: keff 
                                                                        
      INTEGER(i4) :: ibl,im,ix,mxb,myb,mb 
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: btmp,jtmp 
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE ::                      &
     &                             bvtmp,jvtmp ,bhtmp,jhtmp ,bitmp,jitmp
      LOGICAL :: n0_found=.false. 
!-----------------------------------------------------------------------
!     loop over modes to find the n=0                                   
!-----------------------------------------------------------------------
      DO im=1,nmodes 
        IF (keff(im)==0) THEN 
          n0_found=.true. 
          EXIT 
        ENDIF 
      ENDDO 
!-----------------------------------------------------------------------
!     loop over blocks and subtract out the equilibrium fields which    
!     were saved                                                        
!-----------------------------------------------------------------------
      DO ibl=1,nrbl 
        mxb=SIZE(rb(1)%be_sv%fs,2) 
        myb=SIZE(rb(1)%be_sv%fs,3) 
        ALLOCATE(btmp(3,mxb,myb),jtmp(3,mxb,myb)) 
        btmp=rb(ibl)%be_sv%fs 
        jtmp=rb(ibl)%ja_sv%fs 
        IF (geom=='tor') THEN 
          WHERE(rb(ibl)%rz%fs(1,:,:)==0) 
            btmp(3,:,:)=0. 
          ELSEWHERE 
            btmp(3,:,:)=btmp(3,:,:)/rb(ibl)%rz%fs(1,:,:) 
          ENDWHERE 
          jtmp(3,:,:)=jtmp(3,:,:)*rb(ibl)%rz%fs(1,:,:) 
        ENDIF 
        rb(ibl)%be%fs(:,:,:,im)=rb(ibl)%be%fs(:,:,:,im)-btmp 
        rb(ibl)%ja%fs(:,:,:,im)=rb(ibl)%ja%fs(:,:,:,im)-jtmp 
        rb(ibl)%ve%fs(:,:,:,im)=rb(ibl)%ve%fs(:,:,:,im)-                &
     &                             rb(ibl)%ve_sv%fs                     
        rb(ibl)%pres%fs(:,:,:,im)=rb(ibl)%pres%fs(:,:,:,im)-            &
     &                               rb(ibl)%pres_sv%fs                 
        rb(ibl)%prese%fs(:,:,:,im)=rb(ibl)%prese%fs(:,:,:,im)-          &
     &                                rb(ibl)%prese_sv%fs               
        rb(ibl)%tele%fs(:,:,:,im)=rb(ibl)%tele%fs(:,:,:,im)-            &
     &                             rb(ibl)%tele_sv%fs                   
        rb(ibl)%tion%fs(:,:,:,im)=rb(ibl)%tion%fs(:,:,:,im)-            &
     &                             rb(ibl)%tion_sv%fs                   
        IF (poly_degree>1) THEN 
          mb=SIZE(rb(1)%be_sv%fsv,2); ALLOCATE(bvtmp(3,mb,mxb,myb-1)) 
          mb=SIZE(rb(1)%be_sv%fsh,2); ALLOCATE(bhtmp(3,mb,mxb-1,myb)) 
          mb=SIZE(rb(1)%be_sv%fsi,2); ALLOCATE(bitmp(3,mb,mxb-1,myb-1)) 
          mb=SIZE(rb(1)%ja_sv%fsv,2); ALLOCATE(jvtmp(3,mb,mxb,myb-1)) 
          mb=SIZE(rb(1)%ja_sv%fsh,2); ALLOCATE(jhtmp(3,mb,mxb-1,myb)) 
          mb=SIZE(rb(1)%ja_sv%fsi,2); ALLOCATE(jitmp(3,mb,mxb-1,myb-1)) 
          bvtmp=rb(ibl)%be_sv%fsv;    jvtmp=rb(ibl)%ja_sv%fsv 
          bhtmp=rb(ibl)%be_sv%fsh;    jhtmp=rb(ibl)%ja_sv%fsh 
          bitmp=rb(ibl)%be_sv%fsi;    jitmp=rb(ibl)%ja_sv%fsi 
          IF (geom=='tor') THEN 
            WHERE(rb(ibl)%rz%fsv(1,:,:,:)==0) 
              bvtmp(3,:,:,:)=0. 
            ELSEWHERE 
              bvtmp(3,:,:,:)=bvtmp(3,:,:,:)/rb(ibl)%rz%fsv(1,:,:,:) 
            ENDWHERE 
            bhtmp(3,:,:,:)=bhtmp(3,:,:,:)/rb(ibl)%rz%fsh(1,:,:,:) 
            bitmp(3,:,:,:)=bitmp(3,:,:,:)/rb(ibl)%rz%fsi(1,:,:,:) 
            jvtmp(3,:,:,:)=jvtmp(3,:,:,:)*rb(ibl)%rz%fsv(1,:,:,:) 
            jhtmp(3,:,:,:)=jhtmp(3,:,:,:)*rb(ibl)%rz%fsh(1,:,:,:) 
            jitmp(3,:,:,:)=jitmp(3,:,:,:)*rb(ibl)%rz%fsi(1,:,:,:) 
          ENDIF 
          rb(ibl)%be%fsv(:,:,:,:,im)=rb(ibl)%be%fsv(:,:,:,:,im)-        &
     &                                  bvtmp                           
          rb(ibl)%be%fsh(:,:,:,:,im)=rb(ibl)%be%fsh(:,:,:,:,im)-        &
     &                                  bhtmp                           
          rb(ibl)%be%fsi(:,:,:,:,im)=rb(ibl)%be%fsi(:,:,:,:,im)-        &
     &                                  bitmp                           
          rb(ibl)%ja%fsv(:,:,:,:,im)=rb(ibl)%ja%fsv(:,:,:,:,im)-        &
     &                                  jvtmp                           
          rb(ibl)%ja%fsh(:,:,:,:,im)=rb(ibl)%ja%fsh(:,:,:,:,im)-        &
     &                                  jhtmp                           
          rb(ibl)%ja%fsi(:,:,:,:,im)=rb(ibl)%ja%fsi(:,:,:,:,im)-        &
     &                                  jitmp                           
          rb(ibl)%ve%fsh(:,:,:,:,im)=rb(ibl)%ve%fsh(:,:,:,:,im)-        &
     &                                  rb(ibl)%ve_sv%fsh               
          rb(ibl)%ve%fsv(:,:,:,:,im)=rb(ibl)%ve%fsv(:,:,:,:,im)-        &
     &                                  rb(ibl)%ve_sv%fsv               
          rb(ibl)%ve%fsi(:,:,:,:,im)=rb(ibl)%ve%fsi(:,:,:,:,im)-        &
     &                                  rb(ibl)%ve_sv%fsi               
                                                                        
         rb(ibl)%pres%fsh(:,:,:,:,im)=rb(ibl)%pres%fsh(:,:,:,:,im)-     &
     &                                rb(ibl)%pres_sv%fsh               
         rb(ibl)%pres%fsv(:,:,:,:,im)=rb(ibl)%pres%fsv(:,:,:,:,im)-     &
     &                                rb(ibl)%pres_sv%fsv               
         rb(ibl)%pres%fsi(:,:,:,:,im)=rb(ibl)%pres%fsi(:,:,:,:,im)-     &
     &                                rb(ibl)%pres_sv%fsi               
                                                                        
         rb(ibl)%prese%fsh(:,:,:,:,im)=rb(ibl)%prese%fsh(:,:,:,:,im)-   &
     &                                 rb(ibl)%prese_sv%fsh             
         rb(ibl)%prese%fsv(:,:,:,:,im)=rb(ibl)%prese%fsv(:,:,:,:,im)-   &
     &                                 rb(ibl)%prese_sv%fsv             
         rb(ibl)%prese%fsi(:,:,:,:,im)=rb(ibl)%prese%fsi(:,:,:,:,im)-   &
     &                                 rb(ibl)%prese_sv%fsi             
                                                                        
          rb(ibl)%tele%fsh(:,:,:,:,im)=rb(ibl)%tele%fsh(:,:,:,:,im)-    &
     &                                 rb(ibl)%tele_sv%fsh              
          rb(ibl)%tele%fsv(:,:,:,:,im)=rb(ibl)%tele%fsv(:,:,:,:,im)-    &
     &                                 rb(ibl)%tele_sv%fsv              
          rb(ibl)%tele%fsi(:,:,:,:,im)=rb(ibl)%tele%fsi(:,:,:,:,im)-    &
     &                                 rb(ibl)%tele_sv%fsi              
                                                                        
          rb(ibl)%tion%fsh(:,:,:,:,im)=rb(ibl)%tion%fsh(:,:,:,:,im)-    &
     &                                 rb(ibl)%tion_sv%fsh              
          rb(ibl)%tion%fsv(:,:,:,:,im)=rb(ibl)%tion%fsv(:,:,:,:,im)-    &
     &                                 rb(ibl)%tion_sv%fsv              
          rb(ibl)%tion%fsi(:,:,:,:,im)=rb(ibl)%tion%fsi(:,:,:,:,im)-    &
     &                                 rb(ibl)%tion_sv%fsi              
          DEALLOCATE(bvtmp,bitmp,bhtmp,jvtmp,jitmp,jhtmp) 
        ENDIF 
        DEALLOCATE(btmp,jtmp) 
      ENDDO 
!-----------------------------------------------------------------------
!-PRE tblocks, linear only for now.                                     
!-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl 
        tb(ibl)%be%fs(1:2,:,:,im)=tb(ibl)%be%fs(1:2,:,:,im)-            &
     &                            tb(ibl)%be_sv%fs(1:2,:,:)             
        IF (geom=='tor') THEN 
          DO ix=0,tb(ibl)%mvert 
            IF (tb(ibl)%tgeom%xs(ix)==0) THEN 
              tb(ibl)%be%fs(3,ix,0,im)=0 
            ELSE 
              tb(ibl)%be%fs(3,ix,0,im)=tb(ibl)%be%fs(3,ix,0,im)-        &
     &                tb(ibl)%be_sv%fs(3,ix,0)/tb(ibl)%tgeom%xs(ix)     
            ENDIF 
          ENDDO 
        ELSE 
          tb(ibl)%be%fs(3,:,:,im)=tb(ibl)%be%fs(3,:,:,im)-              &
     &                            tb(ibl)%be_sv%fs(3,:,:)               
        ENDIF 
        tb(ibl)%ve%fs(:,:,:,im)=tb(ibl)%ve%fs(:,:,:,im)-                &
     &                          tb(ibl)%ve_sv%fs                        
        tb(ibl)%ja%fs(:,:,:,im)=tb(ibl)%ja%fs(:,:,:,im)-                &
     &                          tb(ibl)%ja_sv%fs                        
        tb(ibl)%pres%fs(:,:,:,im)=tb(ibl)%pres%fs(:,:,:,im)-            &
     &                            tb(ibl)%pres_sv%fs                    
        tb(ibl)%prese%fs(:,:,:,im)=tb(ibl)%prese%fs(:,:,:,im)-          &
     &                             tb(ibl)%prese_sv%fs                  
        tb(ibl)%tele%fs(:,:,:,im)=tb(ibl)%tele%fs(:,:,:,im)-            &
     &                            tb(ibl)%tele_sv%fs                    
        tb(ibl)%tion%fs(:,:,:,im)=tb(ibl)%tion%fs(:,:,:,im)-            &
     &                            tb(ibl)%tion_sv%fs(:,:,:)             
      ENDDO 
!-----------------------------------------------------------------------
!     Now put the equilibrium arrays to the way they were               
!     This is done for every processor                                  
!-----------------------------------------------------------------------
      DO ibl=1,nrbl 
        ! Don't do nd_eq.  Screws up continuity options                 
        ! We've already calculate temperatures here                     
        rb(ibl)%be_eq=rb(ibl)%be_sv; rb(ibl)%be_sv=0 
        rb(ibl)%ja_eq=rb(ibl)%ja_sv; rb(ibl)%ja_sv=0 
        rb(ibl)%ve_eq=rb(ibl)%ve_sv; rb(ibl)%ve_sv=0 
        rb(ibl)%pres_eq=rb(ibl)%pres_sv; rb(ibl)%pres_sv=0 
        rb(ibl)%prese_eq=rb(ibl)%prese_sv; rb(ibl)%prese_sv=0 
        rb(ibl)%tele_eq=rb(ibl)%tele_sv; rb(ibl)%tele_sv=0 
        rb(ibl)%tion_eq=rb(ibl)%tion_sv; rb(ibl)%tion_sv=0 
      ENDDO 
      DO ibl=nrbl+1,nbl 
        tb(ibl)%be_eq=tb(ibl)%be_sv; tb(ibl)%be_sv=0 
        tb(ibl)%ja_eq=tb(ibl)%ja_sv; tb(ibl)%ja_sv=0 
        tb(ibl)%ve_eq=tb(ibl)%ve_sv; tb(ibl)%ve_sv=0 
        tb(ibl)%pres_eq=tb(ibl)%pres_sv; tb(ibl)%pres_sv=0 
        tb(ibl)%prese_eq=tb(ibl)%prese_sv; tb(ibl)%prese_sv=0 
        tb(ibl)%tele_eq=tb(ibl)%tele_sv; tb(ibl)%tele_sv=0 
        tb(ibl)%tion_eq=tb(ibl)%tion_sv; tb(ibl)%tion_sv=0 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE swap_back 
      END MODULE rt_swap_mod
