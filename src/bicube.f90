!-----------------------------------------------------------------------
!     $Id: bicube.f90 4517 2015-07-29 23:08:48Z jking $
!     fits functions to bicubic splines.                                
!     Reference: H. Spaeth, "Spline Algorithms for Curves and Surfaces,"
!     Translated from the German by W. D. Hoskins and H. W. Sager.      
!     Utilitas Mathematica Publishing Inc., Winnepeg, 1974.             
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.                                                
!-----------------------------------------------------------------------
!      0. bicube_type definition.                                       
!      1. bicube_alloc.                                                                                           
!      3. bicube_fit.                                                   
!      4. bicube_eval.                                                  
!      5. bicube_getco.                                                 
!      6. bicube_all_eval.                                              
!      7. bicube_all_getco.                                             
!      8. bicube_write_xy.                                              
!      9. bicube_write_yx.                                              
!     10. bicube_write_arrays.                                                                                            
!     12. bicube_assign_rsc.                                            
!     13. bicube_assign_bc.                                             
!     14. bicube_assign_int.                                            
!     15. dump_read_bicube.                                             
!     16. dump_write_bicube.                                            
!-----------------------------------------------------------------------
!     subprogram 0. bicube_type definition.                             
!     defines bicube_type.                                              
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      MODULE bicube 
      USE local 
      USE spline 
      IMPLICIT NONE 
                                                                        
      TYPE :: bicube_type 
      INTEGER(i4) :: mx,my,nqty,ix,iy 
      REAL(r8), DIMENSION(:), POINTER :: xs,ys 
      REAL(r8), DIMENSION(:,:,:), POINTER :: fs,fsx,fsy,fsxy 
      REAL(r8), DIMENSION(:), POINTER :: f,fx,fy,fxx,fxy,fyy 
      REAL(r8), DIMENSION(:,:,:,:,:), POINTER :: cmats 
      CHARACTER(6), DIMENSION(:), POINTER :: title 
      CHARACTER(6) :: name 
      END TYPE bicube_type 
                                                                        
!-----------------------------------------------------------------------
!     overload asignment and operators.                                 
!-----------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=) 
        MODULE PROCEDURE bicube_assign_rsc,bicube_assign_bc,            &
     &    bicube_assign_int                                             
      END INTERFACE 
                                                                        
      CONTAINS 
!-----------------------------------------------------------------------
!     subprogram 1. bicube_alloc.                                       
!     allocates space for bicube_type.                                  
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE bicube_alloc(bcs,mx,my,nqty) 
                                                                        
      INTEGER(i4), INTENT(IN) :: mx,my,nqty 
      TYPE(bicube_type), INTENT(OUT) :: bcs 
!-----------------------------------------------------------------------
!     allocate space.                                                   
!-----------------------------------------------------------------------
      bcs%mx=mx 
      bcs%my=my 
      bcs%ix=0 
      bcs%iy=0 
      bcs%nqty=nqty 
      ALLOCATE(bcs%xs(0:mx)) 
      ALLOCATE(bcs%ys(0:my)) 
      ALLOCATE(bcs%fs(nqty,0:mx,0:my)) 
      ALLOCATE(bcs%fsx(nqty,0:mx,0:my)) 
      ALLOCATE(bcs%fsy(nqty,0:mx,0:my)) 
      ALLOCATE(bcs%fsxy(nqty,0:mx,0:my)) 
      ALLOCATE(bcs%title(nqty)) 
      ALLOCATE(bcs%f(nqty)) 
      ALLOCATE(bcs%fx(nqty)) 
      ALLOCATE(bcs%fy(nqty)) 
      ALLOCATE(bcs%fxx(nqty)) 
      ALLOCATE(bcs%fxy(nqty)) 
      ALLOCATE(bcs%fyy(nqty)) 
      NULLIFY(bcs%cmats) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE bicube_alloc 
!-----------------------------------------------------------------------
!     subprogram 3. bicube_fit.                                         
!     fits functions to bicubic splines.                                
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE bicube_fit(bcs,endmode1,endmode2) 
                                                                        
      TYPE(bicube_type), INTENT(INOUT) :: bcs 
      CHARACTER(*), INTENT(IN) :: endmode1,endmode2 
                                                                        
      INTEGER(i4) :: iqty 
      TYPE(spline_type) :: spl 
!-----------------------------------------------------------------------
!     evaluate y derivatives.                                           
!-----------------------------------------------------------------------
      CALL spline_alloc(spl,bcs%my,bcs%mx+1_i4) 
      spl%xs=bcs%ys 
      DO iqty=1,bcs%nqty 
         spl%fs=TRANSPOSE(bcs%fs(iqty,:,:)) 
         CALL spline_fit(spl,endmode2) 
         bcs%fsy(iqty,:,:)=TRANSPOSE(spl%fs1) 
      ENDDO 
      CALL spline_dealloc(spl) 
!-----------------------------------------------------------------------
!     evaluate x derivatives.                                           
!-----------------------------------------------------------------------
      CALL spline_alloc(spl,bcs%mx,bcs%my+1_i4) 
      spl%xs=bcs%xs 
      DO iqty=1,bcs%nqty 
         spl%fs=bcs%fs(iqty,:,:) 
         CALL spline_fit(spl,endmode1) 
         bcs%fsx(iqty,:,:)=spl%fs1 
      ENDDO 
!-----------------------------------------------------------------------
!     evaluate mixed derivatives.                                       
!-----------------------------------------------------------------------
      DO iqty=1,bcs%nqty 
         spl%fs=bcs%fsy(iqty,:,:) 
         CALL spline_fit(spl,endmode1) 
         bcs%fsxy(iqty,:,:)=spl%fs1 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      CALL spline_dealloc(spl) 
      RETURN 
      END SUBROUTINE bicube_fit 
!-----------------------------------------------------------------------
!     subprogram 4. bicube_eval.                                        
!     evaluates bicubic spline function.                                
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE bicube_eval(bcs,x,y,mode) 
                                                                        
      TYPE(bicube_type), INTENT(INOUT) :: bcs 
      REAL(r8), INTENT(IN) :: x,y 
      INTEGER(i4), INTENT(IN) :: mode 
                                                                        
      INTEGER(i4) :: i 
      REAL(r8) :: dx,dy 
      REAL(r8), DIMENSION (4,4,bcs%nqty) :: c 
!-----------------------------------------------------------------------
!     find x interval.                                                  
!-----------------------------------------------------------------------
      bcs%ix=max(bcs%ix,0_i4) 
      bcs%ix=min(bcs%ix,bcs%mx-1_i4) 
      DO 
         IF(x >= bcs%xs(bcs%ix) .OR. bcs%ix <= 0)EXIT 
         bcs%ix=bcs%ix-1 
      ENDDO 
      DO 
         IF(x < bcs%xs(bcs%ix+1) .OR. bcs%ix >= bcs%mx-1)EXIT 
         bcs%ix=bcs%ix+1 
      ENDDO 
!-----------------------------------------------------------------------
!     find y interval.                                                  
!-----------------------------------------------------------------------
      bcs%iy=max(bcs%iy,0_i4) 
      bcs%iy=min(bcs%iy,bcs%my-1_i4) 
      DO 
         IF(y >= bcs%ys(bcs%iy) .OR. bcs%iy <= 0)EXIT 
         bcs%iy=bcs%iy-1 
      ENDDO 
      DO 
         IF(y < bcs%ys(bcs%iy+1) .OR. bcs%iy >= bcs%my-1)EXIT 
         bcs%iy=bcs%iy+1 
      ENDDO 
!-----------------------------------------------------------------------
!     find offsets and compute local coefficients.                      
!-----------------------------------------------------------------------
      dx=x-bcs%xs(bcs%ix) 
      dy=y-bcs%ys(bcs%iy) 
      IF(ASSOCIATED(bcs%cmats))THEN 
         c=bcs%cmats(:,:,:,bcs%ix+1,bcs%iy+1) 
      ELSE 
         c=bicube_getco(bcs) 
      ENDIF 
!-----------------------------------------------------------------------
!     evaluate f.                                                       
!-----------------------------------------------------------------------
      bcs%f=0 
      DO i=4,1,-1 
         bcs%f=bcs%f*dx                                                 &
     &        +((c(i,4,:)*dy                                            &
     &        +c(i,3,:))*dy                                             &
     &        +c(i,2,:))*dy                                             &
     &        +c(i,1,:)                                                 
      ENDDO 
      IF(mode == 0)RETURN 
!-----------------------------------------------------------------------
!     evaluate first derivatives of f                                   
!-----------------------------------------------------------------------
      bcs%fx=0 
      bcs%fy=0 
      DO i=4,1,-1 
         bcs%fy=bcs%fy*dx                                               &
     &        +(c(i,4,:)*3*dy                                           &
     &        +c(i,3,:)*2)*dy                                           &
     &        +c(i,2,:)                                                 
         bcs%fx=bcs%fx*dy                                               &
     &        +(c(4,i,:)*3*dx                                           &
     &        +c(3,i,:)*2)*dx                                           &
     &        +c(2,i,:)                                                 
      ENDDO 
      IF(mode == 1)RETURN 
!-----------------------------------------------------------------------
!     evaluate second derivatives of f                                  
!-----------------------------------------------------------------------
      bcs%fxx=0 
      bcs%fyy=0 
      bcs%fxy=0 
      DO i=4,1,-1 
         bcs%fyy=bcs%fyy*dx                                             &
     &        +(c(i,4,:)*3*dy                                           &
     &        +c(i,3,:))*2                                              
         bcs%fxx=bcs%fxx*dy                                             &
     &        +(c(4,i,:)*3*dx                                           &
     &        +c(3,i,:))*2                                              
      ENDDO 
      DO i=4,2,-1 
         bcs%fxy=bcs%fxy*dx                                             &
     &        +((c(i,4,:)*3*dy                                          &
     &        +c(i,3,:)*2)*dy                                           &
     &        +c(i,2,:))*(i-1)                                          
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE bicube_eval 
!-----------------------------------------------------------------------
!     subprogram 5. bicube_getco.                                       
!     computes coefficient matrices.                                    
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      FUNCTION bicube_getco(bcs) RESULT(cmat) 
                                                                        
      TYPE(bicube_type), INTENT(IN) :: bcs 
      REAL(r8), DIMENSION(4,4,bcs%nqty) :: cmat 
                                                                        
      REAL(r8) :: hxfac,hxfac2,hxfac3 
      REAL(r8) :: hyfac,hyfac2,hyfac3 
      REAL(r8), DIMENSION(3:4,4) :: gxmat,gymat 
      REAL(r8), DIMENSION(4,4,bcs%nqty) :: temp 
!-----------------------------------------------------------------------
!     compute gxmat.                                                    
!-----------------------------------------------------------------------
      hxfac=1/(bcs%xs(bcs%ix+1)-bcs%xs(bcs%ix)) 
      hxfac2=hxfac*hxfac 
      hxfac3=hxfac2*hxfac 
      gxmat(3,1)=-3*hxfac2 
      gxmat(3,2)=-2*hxfac 
      gxmat(3,3)=3*hxfac2 
      gxmat(3,4)=-hxfac 
      gxmat(4,1)=2*hxfac3 
      gxmat(4,2)=hxfac2 
      gxmat(4,3)=-2*hxfac3 
      gxmat(4,4)=hxfac2 
!-----------------------------------------------------------------------
!     compute gymat.                                                    
!-----------------------------------------------------------------------
      hyfac=1/(bcs%ys(bcs%iy+1)-bcs%ys(bcs%iy)) 
      hyfac2=hyfac*hyfac 
      hyfac3=hyfac2*hyfac 
      gymat(3,1)=-3*hyfac2 
      gymat(3,2)=-2*hyfac 
      gymat(3,3)=3*hyfac2 
      gymat(3,4)=-hyfac 
      gymat(4,1)=2*hyfac3 
      gymat(4,2)=hyfac2 
      gymat(4,3)=-2*hyfac3 
      gymat(4,4)=hyfac2 
!-----------------------------------------------------------------------
!     compute smat.                                                     
!-----------------------------------------------------------------------
      cmat(1,1,:)=bcs%fs(:,bcs%ix,bcs%iy) 
      cmat(1,2,:)=bcs%fsy(:,bcs%ix,bcs%iy) 
      cmat(1,3,:)=bcs%fs(:,bcs%ix,bcs%iy+1) 
      cmat(1,4,:)=bcs%fsy(:,bcs%ix,bcs%iy+1) 
      cmat(2,1,:)=bcs%fsx(:,bcs%ix,bcs%iy) 
      cmat(2,2,:)=bcs%fsxy(:,bcs%ix,bcs%iy) 
      cmat(2,3,:)=bcs%fsx(:,bcs%ix,bcs%iy+1) 
      cmat(2,4,:)=bcs%fsxy(:,bcs%ix,bcs%iy+1) 
      cmat(3,1,:)=bcs%fs(:,bcs%ix+1,bcs%iy) 
      cmat(3,2,:)=bcs%fsy(:,bcs%ix+1,bcs%iy) 
      cmat(3,3,:)=bcs%fs(:,bcs%ix+1,bcs%iy+1) 
      cmat(3,4,:)=bcs%fsy(:,bcs%ix+1,bcs%iy+1) 
      cmat(4,1,:)=bcs%fsx(:,bcs%ix+1,bcs%iy) 
      cmat(4,2,:)=bcs%fsxy(:,bcs%ix+1,bcs%iy) 
      cmat(4,3,:)=bcs%fsx(:,bcs%ix+1,bcs%iy+1) 
      cmat(4,4,:)=bcs%fsxy(:,bcs%ix+1,bcs%iy+1) 
!-----------------------------------------------------------------------
!     multiply by gymat^T.                                              
!-----------------------------------------------------------------------
      temp(:,1:2,:)=cmat(:,1:2,:) 
      temp(:,3,:)                                                       &
     &     =cmat(:,1,:)*gymat(3,1)                                      &
     &     +cmat(:,2,:)*gymat(3,2)                                      &
     &     +cmat(:,3,:)*gymat(3,3)                                      &
     &     +cmat(:,4,:)*gymat(3,4)                                      
      temp(:,4,:)                                                       &
     &     =cmat(:,1,:)*gymat(4,1)                                      &
     &     +cmat(:,2,:)*gymat(4,2)                                      &
     &     +cmat(:,3,:)*gymat(4,3)                                      &
     &     +cmat(:,4,:)*gymat(4,4)                                      
!-----------------------------------------------------------------------
!     multiply by gxmat.                                                
!-----------------------------------------------------------------------
      cmat(1:2,:,:)=temp(1:2,:,:) 
      cmat(3,:,:)                                                       &
     &     =gxmat(3,1)*temp(1,:,:)                                      &
     &     +gxmat(3,2)*temp(2,:,:)                                      &
     &     +gxmat(3,3)*temp(3,:,:)                                      &
     &     +gxmat(3,4)*temp(4,:,:)                                      
      cmat(4,:,:)                                                       &
     &     =gxmat(4,1)*temp(1,:,:)                                      &
     &     +gxmat(4,2)*temp(2,:,:)                                      &
     &     +gxmat(4,3)*temp(3,:,:)                                      &
     &     +gxmat(4,4)*temp(4,:,:)                                      
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END FUNCTION bicube_getco 
!-----------------------------------------------------------------------
!     subprogram 6. bicube_all_eval.                                    
!     evaluates bicubic splines in all intervals for equal spacing.     
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE bicube_all_eval(bcs,dx,dy,f,fx,fy,fxx,fyy,fxy,mode) 
                                                                        
      TYPE(bicube_type), INTENT(INOUT) :: bcs 
      REAL(r8), INTENT(IN) :: dx,dy 
      REAL(r8), INTENT(OUT), DIMENSION(:,:,:) ::                        &
     &     f,fx,fy,fxx,fyy,fxy                                          
      INTEGER(i4), INTENT(IN) :: mode 
                                                                        
      INTEGER(i4) :: i,ix,iy 
      REAL(r8), DIMENSION(bcs%mx) :: dxv 
      REAL(r8), DIMENSION(bcs%my) :: dyv 
!-----------------------------------------------------------------------
!     compute local displacements and coefficients.                     
!-----------------------------------------------------------------------
      dxv=(bcs%xs(1:bcs%mx)-bcs%xs(0:bcs%mx-1))*dx 
      dyv=(bcs%ys(1:bcs%my)-bcs%ys(0:bcs%my-1))*dy 
      CALL bicube_all_getco(bcs) 
!-----------------------------------------------------------------------
!     evaluate f.                                                       
!-----------------------------------------------------------------------
      f=0 
      DO i=4,1,-1 
         IF(i /= 4)THEN 
            DO ix=1,bcs%mx 
               f(:,ix,:)=f(:,ix,:)*dxv(ix) 
            ENDDO 
         ENDIF 
         DO iy=1,bcs%my 
            f(:,:,iy)=f(:,:,iy)                                         &
     &           +((bcs%cmats(i,4,:,:,iy)*dyv(iy)                       &
     &            +bcs%cmats(i,3,:,:,iy))*dyv(iy)                       &
     &            +bcs%cmats(i,2,:,:,iy))*dy                            &
     &            +bcs%cmats(i,1,:,:,iy)                                
         ENDDO 
      ENDDO 
      IF(mode == 0)RETURN 
!-----------------------------------------------------------------------
!     evaluate fx.                                                      
!-----------------------------------------------------------------------
      fx=0 
      DO i=4,1,-1 
         IF(i /= 4)THEN 
            DO iy=1,bcs%my 
               fx(:,:,iy)=fx(:,:,iy)*dyv(iy) 
            ENDDO 
         ENDIF 
         DO ix=1,bcs%mx 
            fx(:,ix,:)=fx(:,ix,:)                                       &
     &           +(bcs%cmats(4,i,:,ix,:)*3*dxv(ix)                      &
     &           +bcs%cmats(3,i,:,ix,:)*2)*dxv(ix)                      &
     &           +bcs%cmats(2,i,:,ix,:)                                 
         ENDDO 
      ENDDO 
!-----------------------------------------------------------------------
!     evaluate fy.                                                      
!-----------------------------------------------------------------------
      fy=0 
      DO i=4,1,-1 
         IF(i /= 4)THEN 
            DO ix=1,bcs%mx 
               fy(:,ix,:)=fy(:,ix,:)*dxv(ix) 
            ENDDO 
         ENDIF 
         DO iy=1,bcs%my 
            fy(:,:,iy)=fy(:,:,iy)                                       &
     &           +(bcs%cmats(i,4,:,:,iy)*3*dyv(iy)                      &
     &           +bcs%cmats(i,3,:,:,iy)*2)*dyv(iy)                      &
     &           +bcs%cmats(i,2,:,:,iy)                                 
         ENDDO 
      ENDDO 
      IF(mode == 1)RETURN 
!-----------------------------------------------------------------------
!     evaluate fxx.                                                     
!-----------------------------------------------------------------------
      fxx=0 
      DO i=4,1,-1 
         IF(i /= 4)THEN 
            DO iy=1,bcs%my 
               fxx(:,:,iy)=fxx(:,:,iy)*dyv(iy) 
            ENDDO 
         ENDIF 
         DO ix=1,bcs%mx 
            fxx(:,ix,:)=fxx(:,ix,:)                                     &
     &           +(bcs%cmats(4,i,:,ix,:)*3*dxv(ix)                      &
     &           +bcs%cmats(3,i,:,ix,:))*2                              
         ENDDO 
      ENDDO 
!-----------------------------------------------------------------------
!     evaluate fyy and fxy                                              
!-----------------------------------------------------------------------
      fyy=0 
      DO i=4,1,-1 
         IF(i /= 4)THEN 
            DO ix=1,bcs%mx 
               fyy(:,ix,:)=fyy(:,ix,:)*dxv(ix) 
               fxy(:,ix,:)=fxy(:,ix,:)*dxv(ix) 
            ENDDO 
         ENDIF 
         DO iy=1,bcs%my 
            fyy(:,:,iy)=fyy(:,:,iy)                                     &
     &           +(bcs%cmats(i,4,:,:,iy)*3*dyv(iy)                      &
     &           +bcs%cmats(i,3,:,:,iy))*2                              
            fxy(:,:,iy)=fxy(:,:,iy)                                     &
     &           +((bcs%cmats(i,4,:,:,iy)*3*dyv(iy)                     &
     &           +bcs%cmats(i,3,:,:,iy)*2)*dyv(iy)                      &
     &           +bcs%cmats(i,2,:,:,iy))*(i-1)                          
         ENDDO 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE bicube_all_eval 
!-----------------------------------------------------------------------
!     subprogram 7. bicube_all_getco.                                   
!     computes coefficient matrices.                                    
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE bicube_all_getco(bcs) 
                                                                        
      TYPE(bicube_type), INTENT(INOUT) :: bcs 
                                                                        
      INTEGER(i4) :: ix,iy 
      REAL(r8), DIMENSION(bcs%mx) :: hxfac,hxfac2,hxfac3 
      REAL(r8), DIMENSION(bcs%my) :: hyfac,hyfac2,hyfac3 
      REAL(r8), DIMENSION(3:4,4,bcs%mx) :: gxmat 
      REAL(r8), DIMENSION(3:4,4,bcs%my) :: gymat 
      REAL(r8), DIMENSION(4,4,bcs%nqty,bcs%mx,bcs%my) :: temp 
!-----------------------------------------------------------------------
!     allocate space.                                                   
!-----------------------------------------------------------------------
      IF(ASSOCIATED(bcs%cmats))THEN 
         RETURN 
      ELSE 
         ALLOCATE(bcs%cmats(4,4,bcs%nqty,bcs%mx,bcs%my)) 
      ENDIF 
!-----------------------------------------------------------------------
!     compute gxmat.                                                    
!-----------------------------------------------------------------------
      hxfac=1/(bcs%xs(1:bcs%mx)-bcs%xs(0:bcs%mx-1)) 
      hxfac2=hxfac*hxfac 
      hxfac3=hxfac2*hxfac 
      gxmat(3,1,:)=-3*hxfac2 
      gxmat(3,2,:)=-2*hxfac 
      gxmat(3,3,:)=3*hxfac2 
      gxmat(3,4,:)=-hxfac 
      gxmat(4,1,:)=2*hxfac3 
      gxmat(4,2,:)=hxfac2 
      gxmat(4,3,:)=-2*hxfac3 
      gxmat(4,4,:)=hxfac2 
!-----------------------------------------------------------------------
!     compute gymat.                                                    
!-----------------------------------------------------------------------
      hyfac=1/(bcs%ys(1:bcs%my)-bcs%ys(0:bcs%my-1)) 
      hyfac2=hyfac*hyfac 
      hyfac3=hyfac2*hyfac 
      gymat(3,1,:)=-3*hyfac2 
      gymat(3,2,:)=-2*hyfac 
      gymat(3,3,:)=3*hyfac2 
      gymat(3,4,:)=-hyfac 
      gymat(4,1,:)=2*hyfac3 
      gymat(4,2,:)=hyfac2 
      gymat(4,3,:)=-2*hyfac3 
      gymat(4,4,:)=hyfac2 
!-----------------------------------------------------------------------
!     compute smat.                                                     
!-----------------------------------------------------------------------
      bcs%cmats(1,1,:,:,:)=bcs%fs(:,0:bcs%mx-1,0:bcs%my-1) 
      bcs%cmats(1,2,:,:,:)=bcs%fsy(:,0:bcs%mx-1,0:bcs%my-1) 
      bcs%cmats(1,3,:,:,:)=bcs%fs(:,0:bcs%mx-1,1:bcs%my) 
      bcs%cmats(1,4,:,:,:)=bcs%fsy(:,0:bcs%mx-1,1:bcs%my) 
      bcs%cmats(2,1,:,:,:)=bcs%fsx(:,0:bcs%mx-1,0:bcs%my-1) 
      bcs%cmats(2,2,:,:,:)=bcs%fsxy(:,0:bcs%mx-1,0:bcs%my-1) 
      bcs%cmats(2,3,:,:,:)=bcs%fsx(:,0:bcs%mx-1,1:bcs%my) 
      bcs%cmats(2,4,:,:,:)=bcs%fsxy(:,0:bcs%mx-1,1:bcs%my) 
      bcs%cmats(3,1,:,:,:)=bcs%fs(:,1:bcs%mx,0:bcs%my-1) 
      bcs%cmats(3,2,:,:,:)=bcs%fsy(:,1:bcs%mx,0:bcs%my-1) 
      bcs%cmats(3,3,:,:,:)=bcs%fs(:,1:bcs%mx,1:bcs%my) 
      bcs%cmats(3,4,:,:,:)=bcs%fsy(:,1:bcs%mx,1:bcs%my) 
      bcs%cmats(4,1,:,:,:)=bcs%fsx(:,1:bcs%mx,0:bcs%my-1) 
      bcs%cmats(4,2,:,:,:)=bcs%fsxy(:,1:bcs%mx,0:bcs%my-1) 
      bcs%cmats(4,3,:,:,:)=bcs%fsx(:,1:bcs%mx,1:bcs%my) 
      bcs%cmats(4,4,:,:,:)=bcs%fsxy(:,1:bcs%mx,1:bcs%my) 
!-----------------------------------------------------------------------
!     multiply by gymat^T.                                              
!-----------------------------------------------------------------------
      temp(:,1:2,:,:,:)=bcs%cmats(:,1:2,:,:,:) 
      DO iy=1,bcs%my 
         temp(:,3,:,:,iy)                                               &
     &        =bcs%cmats(:,1,:,:,iy)*gymat(3,1,iy)                      &
     &        +bcs%cmats(:,2,:,:,iy)*gymat(3,2,iy)                      &
     &        +bcs%cmats(:,3,:,:,iy)*gymat(3,3,iy)                      &
     &        +bcs%cmats(:,4,:,:,iy)*gymat(3,4,iy)                      
         temp(:,4,:,:,iy)                                               &
     &        =bcs%cmats(:,1,:,:,iy)*gymat(4,1,iy)                      &
     &        +bcs%cmats(:,2,:,:,iy)*gymat(4,2,iy)                      &
     &        +bcs%cmats(:,3,:,:,iy)*gymat(4,3,iy)                      &
     &        +bcs%cmats(:,4,:,:,iy)*gymat(4,4,iy)                      
      ENDDO 
!-----------------------------------------------------------------------
!     multiply by gxmat.                                                
!-----------------------------------------------------------------------
      bcs%cmats(1:2,:,:,:,:)=temp(1:2,:,:,:,:) 
      DO ix=1,bcs%mx 
         bcs%cmats(3,:,:,ix,:)                                          &
     &        =gxmat(3,1,ix)*temp(1,:,:,ix,:)                           &
     &        +gxmat(3,2,ix)*temp(2,:,:,ix,:)                           &
     &        +gxmat(3,3,ix)*temp(3,:,:,ix,:)                           &
     &        +gxmat(3,4,ix)*temp(4,:,:,ix,:)                           
         bcs%cmats(4,:,:,ix,:)                                          &
     &        =gxmat(4,1,ix)*temp(1,:,:,ix,:)                           &
     &        +gxmat(4,2,ix)*temp(2,:,:,ix,:)                           &
     &        +gxmat(4,3,ix)*temp(3,:,:,ix,:)                           &
     &        +gxmat(4,4,ix)*temp(4,:,:,ix,:)                           
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE bicube_all_getco 
!-----------------------------------------------------------------------
!     subprogram 8. bicube_write_xy.                                    
!     produces ascii and binarx output for bicubic spline fits.         
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE bicube_write_xy(bcs,out,bin,iua,iub) 
                                                                        
      TYPE(bicube_type), INTENT(INOUT) :: bcs 
      LOGICAL, INTENT(IN) :: out,bin 
      INTEGER, INTENT(IN) :: iua,iub 
                                                                        
      INTEGER(i4) :: ix,iy,jx,jy,iqty 
      REAL(r8) :: x,y,dx,dy 
                                                                        
      LOGICAL, PARAMETER :: interp=.TRUE. 
      CHARACTER(80) :: format2,format1 
!-----------------------------------------------------------------------
!     write formats.                                                    
!-----------------------------------------------------------------------
   10 FORMAT(1x,'iy = ',i3,', y = ',1p,e11.3) 
   20 FORMAT(1x,'iy = ',i3,', jy = ',i1,', y = ',1p,e11.3) 
   30 FORMAT('(/4x,"ix",6x,"x",4x,',i3.3,'(4x,"f(",i3.3,")",1x)/)') 
   40 FORMAT('(i6,1p,e11.3,',i3.3,'e11.3)') 
!-----------------------------------------------------------------------
!     abort.                                                            
!-----------------------------------------------------------------------
      IF(.NOT. (out .OR. bin))RETURN 
!-----------------------------------------------------------------------
!     write input data.                                                 
!-----------------------------------------------------------------------
      IF(out)THEN 
         WRITE(iua,'(1x,a/)')"input data" 
         WRITE(format1,30)bcs%nqty 
         WRITE(format2,40)bcs%nqty 
      ENDIF 
      DO iy=0,bcs%my 
         y=bcs%ys(iy) 
         IF(out)then 
            WRITE(iua,10)iy,bcs%ys(iy) 
            WRITE(iua,format1)(iqty,iqty=1,bcs%nqty) 
         ENDIF 
         DO ix=0,bcs%mx 
            x=bcs%xs(ix) 
            bcs%f=bcs%fs(ix,iy,:) 
            IF(out)WRITE(iua,format2)ix,x,bcs%f 
            IF(bin)WRITE(iub)REAL(x,4),REAL(bcs%f,4) 
         ENDDO 
         IF(out)WRITE(iua,format1)(iqty,iqty=1,bcs%nqty) 
         IF(bin)WRITE(iub) 
      ENDDO 
!-----------------------------------------------------------------------
!     begin loops over y for interpolated data.                         
!-----------------------------------------------------------------------
      IF(interp)THEN 
         IF(out)WRITE(iua,'(1x,a/)')"interpolated data" 
         DO iy=0,bcs%my-1 
            dy=(bcs%ys(iy+1)-bcs%ys(iy))/4 
            DO jy=0,4 
               y=bcs%ys(iy)+dy*jy 
               IF(out)then 
                  WRITE(iua,20)iy,jy,y 
                  WRITE(iua,format1)(iqty,iqty=1,bcs%nqty) 
               ENDIF 
!-----------------------------------------------------------------------
!     begin loops over x for interpolated data.                         
!-----------------------------------------------------------------------
               DO ix=0,bcs%mx-1 
                  dx=(bcs%xs(ix+1)-bcs%xs(ix))/4 
                  DO jx=0,4 
                     x=bcs%xs(ix)+dx*jx 
                     CALL bicube_eval(bcs,y,x,0_i4) 
                     IF(out)WRITE(iua,format2)ix,x,bcs%f 
                     IF(bin)WRITE(iub)REAL(x,4),REAL(bcs%f,4) 
                  ENDDO 
                  IF(out)WRITE(iua,'()') 
               ENDDO 
!-----------------------------------------------------------------------
!     complete loops over y.                                            
!-----------------------------------------------------------------------
               IF(out)WRITE(iua,format1)(iqty,iqty=1,bcs%nqty) 
               IF(bin)WRITE(iub) 
            ENDDO 
         ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate routine.                                                
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE bicube_write_xy 
!-----------------------------------------------------------------------
!     subprogram 9. bicube_write_yx.                                    
!     produces ascii and binary output for bicubic spline fits.         
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE bicube_write_yx(bcs,out,bin,iua,iub) 
                                                                        
      TYPE(bicube_type), INTENT(INOUT) :: bcs 
      LOGICAL, INTENT(IN) :: out,bin 
      INTEGER, INTENT(IN) :: iua,iub 
                                                                        
      INTEGER(i4) :: ix,iy,jx,jy,iqty 
      REAL(r8) :: x,y,dx,dy 
                                                                        
      LOGICAL, PARAMETER :: interp=.TRUE. 
      CHARACTER(80) :: format1,format2 
!-----------------------------------------------------------------------
!     write formats.                                                    
!-----------------------------------------------------------------------
   10 FORMAT(1x,'ix = ',i3,', x = ',1p,e11.3) 
   20 FORMAT(1x,'ix = ',i3,', jx = ',i1,', x = ',1p,e11.3) 
   30 FORMAT('(/4x,"iy",6x,"y",4x,',i3.3,'(4x,"f(",i3.3,")",1x)/)') 
   40 FORMAT('(i6,1p,e11.3,',i3.3,'e11.3)') 
!-----------------------------------------------------------------------
!     abort.                                                            
!-----------------------------------------------------------------------
      IF(.NOT. (out .OR. bin))RETURN 
!-----------------------------------------------------------------------
!     write input data.                                                 
!-----------------------------------------------------------------------
      IF(out)THEN 
         WRITE(iua,'(1x,a/)')"input data" 
         WRITE(format1,30)bcs%nqty 
         WRITE(format2,40)bcs%nqty 
      ENDIF 
      DO ix=0,bcs%mx 
         x=bcs%xs(ix) 
         IF(out)then 
            WRITE(iua,10)ix,bcs%xs(ix) 
            WRITE(iua,format1)(iqty,iqty=1,bcs%nqty) 
         ENDIF 
         DO iy=0,bcs%my 
            y=bcs%ys(iy) 
            bcs%f=bcs%fs(ix,iy,:) 
            IF(out)WRITE(iua,format2)iy,y,bcs%f 
            IF(bin)WRITE(iub)REAL(y,4),REAL(bcs%f,4) 
         ENDDO 
         IF(out)WRITE(iua,format1)(iqty,iqty=1,bcs%nqty) 
         IF(bin)WRITE(iub) 
      ENDDO 
!-----------------------------------------------------------------------
!     begin loops over x for interpolated data.                         
!-----------------------------------------------------------------------
      IF(interp)THEN 
         IF(out)WRITE(iua,'(1x,a/)')"interpolated data" 
         DO ix=0,bcs%mx-1 
            dx=(bcs%xs(ix+1)-bcs%xs(ix))/4 
            DO jx=0,4 
               x=bcs%xs(ix)+dx*jx 
               IF(out)then 
                  WRITE(iua,20)ix,jx,x 
                  WRITE(iua,format1)(iqty,iqty=1,bcs%nqty) 
               ENDIF 
!-----------------------------------------------------------------------
!     begin loops over y for interpolated data.                         
!-----------------------------------------------------------------------
               DO iy=0,bcs%my-1 
                  dy=(bcs%ys(iy+1)-bcs%ys(iy))/4 
                  DO jy=0,4 
                     y=bcs%ys(iy)+dy*jy 
                     CALL bicube_eval(bcs,x,y,0_i4) 
                     IF(out)WRITE(iua,format2)iy,y,bcs%f 
                     IF(bin)WRITE(iub)REAL(y,4),REAL(bcs%f,4) 
                  ENDDO 
                  IF(out)WRITE(iua,'()') 
               ENDDO 
!-----------------------------------------------------------------------
!     complete loops over x.                                            
!-----------------------------------------------------------------------
               IF(out)WRITE(iua,format1)(iqty,iqty=1,bcs%nqty) 
               IF(bin)WRITE(iub) 
            ENDDO 
         ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate routine.                                                
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE bicube_write_yx 
!-----------------------------------------------------------------------
!     subprogram 10. bicube_write_arrays.                               
!     produces ascii and binary output for bicubic spline fits.         
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE bicube_write_arrays(bcs,out,iua,iqty) 
                                                                        
      TYPE(bicube_type), INTENT(INOUT) :: bcs 
      LOGICAL, INTENT(IN) :: out 
      INTEGER, INTENT(IN) :: iua,iqty 
                                                                        
      CHARACTER(80) :: format1,format2 
      INTEGER(i4) :: ix,iy 
!-----------------------------------------------------------------------
!     formats.                                                          
!-----------------------------------------------------------------------
   10 FORMAT('(/2x,"ix/iy",',i3.3,'(3x,i3.3,5x)/)') 
   20 FORMAT('(i5,1p,',i3.3,'e11.3)') 
!-----------------------------------------------------------------------
!     abort.                                                            
!-----------------------------------------------------------------------
      IF(.NOT. out)RETURN 
!-----------------------------------------------------------------------
!     write fs.                                                         
!-----------------------------------------------------------------------
      WRITE(format1,10)bcs%my+1 
      WRITE(format2,20)bcs%my+1 
      WRITE(iua,"(a)")"fs" 
      WRITE(iua,format1)(iy,iy=0,bcs%my) 
      WRITE(iua,format2)(ix,(bcs%fs(iqty,ix,iy),iy=0,bcs%my),           &
     &     ix=0,bcs%mx)                                                 
      WRITE(iua,format1)(iy,iy=0,bcs%my) 
      WRITE(iua,"(a/)")"fs" 
!-----------------------------------------------------------------------
!     write fsx.                                                        
!-----------------------------------------------------------------------
      WRITE(iua,"(a)")"fsx" 
      WRITE(iua,format1)(iy,iy=0,bcs%my) 
      WRITE(iua,format2)(ix,(bcs%fsx(iqty,ix,iy),iy=0,bcs%my),          &
     &     ix=0,bcs%mx)                                                 
      WRITE(iua,format1)(iy,iy=0,bcs%my) 
      WRITE(iua,"(a/)")"fsx" 
!-----------------------------------------------------------------------
!     write fsy.                                                        
!-----------------------------------------------------------------------
      WRITE(iua,"(a)")"fsy" 
      WRITE(iua,format1)(iy,iy=0,bcs%my) 
      WRITE(iua,format2)(ix,(bcs%fsy(iqty,ix,iy),iy=0,bcs%my),          &
     &     ix=0,bcs%mx)                                                 
      WRITE(iua,format1)(iy,iy=0,bcs%my) 
      WRITE(iua,"(a/)")"fsy" 
!-----------------------------------------------------------------------
!     write fsxy.                                                       
!-----------------------------------------------------------------------
      WRITE(iua,"(a)")"fsxy" 
      WRITE(iua,format1)(iy,iy=0,bcs%my) 
      WRITE(iua,format2)(ix,(bcs%fsxy(iqty,ix,iy),iy=0,bcs%my),         &
     &     ix=0,bcs%mx)                                                 
      WRITE(iua,format1)(iy,iy=0,bcs%my) 
      WRITE(iua,"(a/)")"fsxy" 
!-----------------------------------------------------------------------
!     terminate routine.                                                
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE bicube_write_arrays 
!-----------------------------------------------------------------------
!     subprogram 12. bicube_assign_rsc.                                 
!     assign a real scalar value to a bicube structure.                 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE bicube_assign_rsc(bc,rscalar) 
                                                                        
      TYPE(bicube_type), INTENT(INOUT) :: bc 
      REAL(r8), INTENT(IN) :: rscalar 
                                                                        
!-----------------------------------------------------------------------
!     derivatives are 0.                                                
!-----------------------------------------------------------------------
      bc%fs=rscalar 
      bc%fsx=0 
      bc%fsy=0 
      bc%fsxy=0 
      IF (ASSOCIATED(bc%cmats)) DEALLOCATE(bc%cmats) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE bicube_assign_rsc 
!-----------------------------------------------------------------------
!     subprogram 13. bicube_assign_bc.                                  
!     set one bicube structure array values equal to another.           
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE bicube_assign_bc(bc1,bc2) 
                                                                        
      TYPE(bicube_type), INTENT(INOUT) :: bc1 
      TYPE(bicube_type), INTENT(IN) :: bc2 
                                                                        
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.                   
!-----------------------------------------------------------------------
      bc1%fs=bc2%fs 
      bc1%fsx=bc2%fsx 
      bc1%fsy=bc2%fsy 
      bc1%fsxy=bc2%fsxy 
      IF (ASSOCIATED(bc1%cmats)) DEALLOCATE(bc1%cmats) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE bicube_assign_bc 
!-----------------------------------------------------------------------
!     subprogram 14. bicube_assign_int.                                 
!     assign a integer value to a bicube structure.                     
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE bicube_assign_int(bc,int) 
                                                                        
      TYPE(bicube_type), INTENT(INOUT) :: bc 
      INTEGER(i4), INTENT(IN) :: int 
                                                                        
!-----------------------------------------------------------------------
!     derivatives are 0.                                                
!-----------------------------------------------------------------------
      bc%fs=int 
      bc%fsx=0 
      bc%fsy=0 
      bc%fsxy=0 
      IF (ASSOCIATED(bc%cmats)) DEALLOCATE(bc%cmats) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE bicube_assign_int 
!-----------------------------------------------------------------------
!     subprogram 15. dump_read_bicube.                                  
!-----------------------------------------------------------------------
      SUBROUTINE dump_read_bicube(rb_bc,lx,ly,name,title) 
      USE io 
                                                                        
      INTEGER(i4), INTENT(IN) :: lx,ly 
      CHARACTER(*), INTENT(IN) :: name 
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title 
      TYPE(bicube_type), INTENT(OUT) :: rb_bc 
                                                                        
      REAL(r8) :: iread 
                                                                        
      READ(rstrt_unit) iread 
      rb_bc%nqty=NINT(iread) 
      CALL bicube_alloc(rb_bc,lx,ly,rb_bc%nqty) 
      rb_bc%name=name 
      IF (SIZE(title)<SIZE(rb_bc%title)) THEN 
        rb_bc%title=title(1) 
      ELSE 
        rb_bc%title=title 
      ENDIF 
      READ(rstrt_unit) rb_bc%xs 
      READ(rstrt_unit) rb_bc%ys 
      READ(rstrt_unit) rb_bc%fs 
      READ(rstrt_unit) rb_bc%fsx 
      READ(rstrt_unit) rb_bc%fsy 
      READ(rstrt_unit) rb_bc%fsxy 
                                                                        
      RETURN 
      END SUBROUTINE dump_read_bicube 
!-----------------------------------------------------------------------
!     subprogram 16. dump_write_bicube.                                 
!-----------------------------------------------------------------------
      SUBROUTINE dump_write_bicube(rb_bc) 
      USE io 
                                                                        
      TYPE(bicube_type), INTENT(IN) :: rb_bc 
                                                                        
      WRITE(dump_unit) REAL(rb_bc%nqty,r8) 
      WRITE(dump_unit) rb_bc%xs 
      WRITE(dump_unit) rb_bc%ys 
      WRITE(dump_unit) rb_bc%fs 
      WRITE(dump_unit) rb_bc%fsx 
      WRITE(dump_unit) rb_bc%fsy 
      WRITE(dump_unit) rb_bc%fsxy 
                                                                        
      RETURN 
      END SUBROUTINE dump_write_bicube 
!-----------------------------------------------------------------------
!     close module.                                                     
!-----------------------------------------------------------------------
      END MODULE bicube                                        
