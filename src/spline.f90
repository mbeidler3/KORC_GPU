!-----------------------------------------------------------------------
!     $Id: spline.f90 6919 2019-10-08 20:19:41Z ehowell $
!     fits functions to cubic splines.                                  
!     Reference: H. Spaeth, "Spline Algorithms for Curves and Surfaces,"
!     Translated from the German by W. D. Hoskins and H. W. Sager.      
!     Utilitas Mathematica Publishing Inc., Winnepeg, 1974.             
!                                                                       
!     Spline boundary conditions supported:                             
!       1. Extrap                                                       
!       2. Not-a-Knot                                                   
!       3. Periodic                                                     
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.                                                
!-----------------------------------------------------------------------
!     0. spline_type definition.                                        
!     1. spline_alloc.                                                  
!     2. spline_dealloc.                                                
!     3. spline_fit.                                                    
!     5. spline_eval.                                                   
!     6. spline_all_eval.                                               
!     7. spline_write.                                                  
!     8. spline_int.                                                    
!     PRIVATE FUNCTIONS                                                 
!     4. spline_fac.                                                    
!     9. spline_triluf.                                                 
!     10. spline_trilus.                                                
!     11. spline_sherman.  (used for extrap & N-a-K B.C.s)              
!     12. spline_morrison. (used for periodic B.C.s)                    
!-----------------------------------------------------------------------
!     subprogram 0. spline_type definition.                             
!     defines spline_type.                                              
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      MODULE spline 
      USE local 
      IMPLICIT NONE 
                                                                        
      TYPE :: spline_type 
      INTEGER(i4) :: nodes,nqty,inode 
      ! Grid                                                            
      REAL(r8), DIMENSION(:), POINTER :: xs 
                                                                        
      ! Data, derivatives, and integrated quantities.                   
      REAL(r8), DIMENSION(:,:), POINTER :: fs,fs1,fsi 
                                                                        
      ! Data from evaluations                                           
      REAL(r8), DIMENSION(:), POINTER :: f,f1,f2,f3 
                                                                        
      ! Conveniences lables                                             
      CHARACTER(6), DIMENSION(:), POINTER :: title 
      CHARACTER(6) :: name 
      ! Allocated flag
      LOGICAL :: alloced =.FALSE.
      END TYPE spline_type 
                                                                        
      CONTAINS 
!-----------------------------------------------------------------------
!     subprogram 1. spline_alloc.                                       
!     allocates space for spline_type.                                  
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE spline_alloc(spl,nodes,nqty) 
                                                                        
      INTEGER(i4), INTENT(IN) :: nodes,nqty 
      TYPE(spline_type), INTENT(INOUT) :: spl 
!-----------------------------------------------------------------------
!     allocate space.                                                   
!-----------------------------------------------------------------------
      spl%nodes=nodes 
      spl%nqty=nqty 
      spl%inode=0 
      ALLOCATE(spl%xs(0:nodes)) 
      ALLOCATE(spl%f(nqty)) 
      ALLOCATE(spl%f1(nqty)) 
      ALLOCATE(spl%f2(nqty)) 
      ALLOCATE(spl%f3(nqty)) 
      ALLOCATE(spl%title(0:nqty)) 
      ALLOCATE(spl%fs(0:nodes,nqty)) 
      ALLOCATE(spl%fs1(0:nodes,nqty)) 
      NULLIFY(spl%fsi) 
      spl%alloced=.TRUE.
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE spline_alloc 
!-----------------------------------------------------------------------
!     subprogram 2. spline_dealloc.                                     
!     deallocates space for spline_type.                                
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE spline_dealloc(spl) 
                                                                        
      TYPE(spline_type), INTENT(INOUT) :: spl 
!-----------------------------------------------------------------------
!     deallocate space.                                                 
!-----------------------------------------------------------------------
      IF (spl%alloced) THEN
        DEALLOCATE(spl%xs) 
        DEALLOCATE(spl%f) 
        DEALLOCATE(spl%f1) 
        DEALLOCATE(spl%f2) 
        DEALLOCATE(spl%f3) 
        DEALLOCATE(spl%title) 
        DEALLOCATE(spl%fs) 
        DEALLOCATE(spl%fs1)
        IF (ASSOCIATED(spl%fsi)) DEALLOCATE(spl%fsi)
        spl%alloced=.FALSE.
      END IF
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE spline_dealloc 
!-----------------------------------------------------------------------
!     subprogram 3. spline_fit.                                         
!     fits REAL functions to cubic splines.                             
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE spline_fit(spl,endmode) 
                                                                        
      TYPE(spline_type), INTENT(INOUT) :: spl 
      CHARACTER(*), INTENT(IN) :: endmode 
                                                                        
      INTEGER(i4) :: iqty 
      REAL(r8), DIMENSION(-1:1,0:spl%nodes) :: a 
      REAL(r8), DIMENSION(spl%nodes) :: b 
      REAL(r8), DIMENSION(4) :: cl,cr 
!-----------------------------------------------------------------------
!     set up grid matrix.                                               
!-----------------------------------------------------------------------
      CALL spline_fac(spl,a,b,cl,cr,endmode) 
!-----------------------------------------------------------------------
!     compute first derivatives, interior.                              
!-----------------------------------------------------------------------
      DO iqty=1,spl%nqty 
         spl%fs1(1:spl%nodes-1,iqty)=                                   &
     &        3*((spl%fs(2:spl%nodes,iqty)-spl%fs(1:spl%nodes-1,iqty))  &
     &        *b(2:spl%nodes)                                           &
     &        +(spl%fs(1:spl%nodes-1,iqty)-spl%fs(0:spl%nodes-2,iqty))  &
     &        *b(1:spl%nodes-1))                                        
      ENDDO 
!-----------------------------------------------------------------------
!     extrapolation boundary conditions.                                
!-----------------------------------------------------------------------
      IF(endmode == "extrap")THEN 
         DO iqty=1,spl%nqty 
            spl%fs1(0,iqty)=SUM(cl(1:4)*spl%fs(0:3,iqty)) 
            spl%fs1(spl%nodes,iqty)=SUM(cr(1:4)                         &
     &           *spl%fs(spl%nodes:spl%nodes-3:-1,iqty))                
            spl%fs1(1,iqty)=spl%fs1(1,iqty)-spl%fs1(0,iqty)             &
     &           /(spl%xs(1)-spl%xs(0))                                 
            spl%fs1(spl%nodes-1,iqty)=                                  &
     &           spl%fs1(spl%nodes-1,iqty)-spl%fs1(spl%nodes,iqty)      &
     &           /(spl%xs(spl%nodes)-spl%xs(spl%nodes-1))               
            CALL spline_trilus(a(:,1:),spl%fs1(1:,iqty),                &
     &           spl%nodes-1_i4,1_i4)                                   
         ENDDO 
!-----------------------------------------------------------------------
!     not-a-knot boundary conditions.                                   
!-----------------------------------------------------------------------
      ELSE IF(endmode == "not-a-knot")THEN 
         spl%fs1(1,:)=spl%fs1(1,:)-(2*spl%fs(1,:)                       &
     &        -spl%fs(0,:)-spl%fs(2,:))*2*b(1)                          
         spl%fs1(spl%nodes-1,:)=spl%fs1(spl%nodes-1,:)                  &
     &        +(2*spl%fs(spl%nodes-1,:)-spl%fs(spl%nodes,:)             &
     &        -spl%fs(spl%nodes-2,:))*2*b(spl%nodes)                    
         DO iqty=1,spl%nqty 
            CALL spline_trilus(a(:,1:),spl%fs1(1:,iqty),                &
     &           spl%nodes-1_i4,1_i4)                                   
         ENDDO 
         spl%fs1(0,:)=(2*(2*spl%fs(1,:)-spl%fs(0,:)-spl%fs(2,:))        &
     &        +(spl%fs1(1,:)+spl%fs1(2,:))*(spl%xs(2)-spl%xs(1))        &
     &        -spl%fs1(1,:)*(spl%xs(1)-spl%xs(0)))/(spl%xs(1)-spl%xs(0))
         spl%fs1(spl%nodes,:)=                                          &
     &        (2*(spl%fs(spl%nodes-2,:)+spl%fs(spl%nodes,:)             &
     &        -2*spl%fs(spl%nodes-1,:))                                 &
     &        +(spl%fs1(spl%nodes-1,:)+spl%fs1(spl%nodes-2,:))          &
     &        *(spl%xs(spl%nodes-1)-spl%xs(spl%nodes-2))                &
     &        -spl%fs1(spl%nodes-1,:)                                   &
     &        *(spl%xs(spl%nodes)-spl%xs(spl%nodes-1)))                 &
     &        /(spl%xs(spl%nodes)-spl%xs(spl%nodes-1))                  
!-----------------------------------------------------------------------
!     periodic boudary conditions.                                      
!-----------------------------------------------------------------------
      ELSE IF(endmode == "periodic")THEN 
         spl%fs1(0,:)=3*((spl%fs(1,:)-spl%fs(0,:))*b(1)                 &
     &        +(spl%fs(0,:)-spl%fs(spl%nodes-1,:))*b(spl%nodes))        
         DO iqty=1,spl%nqty 
            CALL spline_morrison(a,spl%fs1(:,iqty),                     &
     &           spl%nodes,1_i4)                                        
         ENDDO 
         spl%fs1(spl%nodes,:)=spl%fs1(0,:) 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE spline_fit 
!-----------------------------------------------------------------------
!     subprogram 4. spline_eval.                                        
!     evaluates REAL cubic spline function.                             
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE spline_eval(spl,x,mode) 
                                                                        
      TYPE(spline_type), INTENT(INOUT) :: spl 
      REAL(r8), INTENT(IN) :: x 
      INTEGER(i4), INTENT(IN) :: mode 
                                                                        
      REAL(r8) :: d,z,z1 
!-----------------------------------------------------------------------
!     find cubic spline interval.                                       
!-----------------------------------------------------------------------
      spl%inode=MAX(spl%inode,0_i4) 
      spl%inode=MIN(spl%inode,spl%nodes-1_i4) 
      DO 
         IF(x >= spl%xs(spl%inode).OR.spl%inode <= 0)EXIT 
         spl%inode=spl%inode-1 
      ENDDO 
      DO 
         IF(x < spl%xs(spl%inode+1).OR.spl%inode >= spl%nodes-1)EXIT 
         spl%inode=spl%inode+1 
      ENDDO 
!-----------------------------------------------------------------------
!     evaluate ofspl%fset and related quantities.                       
!-----------------------------------------------------------------------
      d=spl%xs(spl%inode+1)-spl%xs(spl%inode) 
      z=(x-spl%xs(spl%inode))/d 
      z1=1-z 
!-----------------------------------------------------------------------
!     evaluate output values.                                           
!-----------------------------------------------------------------------
      spl%f(1:spl%nqty)=spl%fs(spl%inode,1:spl%nqty)*z1*z1*(3-2*z1)     &
     &     +spl%fs(spl%inode+1,1:spl%nqty)*z*z*(3-2*z)                  &
     &     +d*z*z1*(spl%fs1(spl%inode,1:spl%nqty)*z1                    &
     &     -spl%fs1(spl%inode+1,1:spl%nqty)*z)                          
      IF(mode == 0)RETURN 
!-----------------------------------------------------------------------
!     evaluate first derivatives.                                       
!-----------------------------------------------------------------------
      spl%f1(1:spl%nqty)=6*(spl%fs(spl%inode+1,1:spl%nqty)              &
     &     -spl%fs(spl%inode,1:spl%nqty))*z*z1/d                        &
     &     +spl%fs1(spl%inode,1:spl%nqty)*z1*(3*z1-2)                   &
     &     +spl%fs1(spl%inode+1,1:spl%nqty)*z*(3*z-2)                   
      IF(mode == 1)RETURN 
!-----------------------------------------------------------------------
!     evaluate second derivatives.                                      
!-----------------------------------------------------------------------
      spl%f2(1:spl%nqty)=(6*(spl%fs(spl%inode+1,1:spl%nqty)             &
     &     -spl%fs(spl%inode,1:spl%nqty))*(z1-z)/d                      &
     &     -spl%fs1(spl%inode,1:spl%nqty)*(6*z1-2)                      &
     &     +spl%fs1(spl%inode+1,1:spl%nqty)*(6*z-2))/d                  
      IF(mode == 2)RETURN 
!-----------------------------------------------------------------------
!     evaluate third derivatives.                                       
!-----------------------------------------------------------------------
      spl%f3(1:spl%nqty)=(12*(spl%fs(spl%inode,1:spl%nqty)              &
     &     -spl%fs(spl%inode+1,1:spl%nqty))/d                           &
     &     +6*(spl%fs1(spl%inode,1:spl%nqty)                            &
     &     +spl%fs1(spl%inode+1,1:spl%nqty)))/(d*d)                     
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE spline_eval 
!-----------------------------------------------------------------------
!     subprogram 5. spline_all_eval.                                    
!     evaluates cubic spline function.                                  
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE spline_all_eval(spl,z,f,f1,f2,f3,mode) 
                                                                        
      TYPE(spline_type), INTENT(INOUT) :: spl 
      REAL(r8), INTENT(IN) :: z 
      REAL(r8), DIMENSION(spl%nodes,spl%nqty), INTENT(OUT) :: f,f1,f2,f3 
      INTEGER(i4), INTENT(IN) :: mode 
                                                                        
      INTEGER(i4) :: iqty,nqty,n 
      REAL(r8) :: z1 
      REAL(r8), DIMENSION(spl%nodes) :: d 
!-----------------------------------------------------------------------
!     evaluate offset and related quantities.                           
!-----------------------------------------------------------------------
      n=spl%nodes 
      nqty=spl%nqty 
      z1=1-z 
      d=spl%xs(1:n)-spl%xs(0:n-1) 
!-----------------------------------------------------------------------
!     evaluate output values.                                           
!-----------------------------------------------------------------------
      DO iqty=1,nqty 
         f(:,iqty)=spl%fs(0:n-1,iqty)*z1*z1*(3-2*z1)                    &
     &        +spl%fs(1:n,iqty)*z*z*(3-2*z)                             &
     &        +d*z*z1*(spl%fs1(0:n-1,iqty)*z1-spl%fs1(1:n,iqty)*z)      
      ENDDO 
      IF(mode == 0)RETURN 
!-----------------------------------------------------------------------
!     evaluate first derivatives.                                       
!-----------------------------------------------------------------------
      DO iqty=1,nqty 
         f1(:,iqty)=6*(spl%fs(1:n,iqty)-spl%fs(0:n-1,iqty))*z*z1/d      &
     &        +spl%fs1(0:n-1,iqty)*z1*(3*z1-2)                          &
     &        +spl%fs1(1:n,iqty)*z*(3*z-2)                              
      ENDDO 
      IF(mode == 1)RETURN 
!-----------------------------------------------------------------------
!     evaluate second derivatives.                                      
!-----------------------------------------------------------------------
      DO iqty=1,nqty 
         f2(:,iqty)=(6*(spl%fs(1:n,iqty)-spl%fs(0:n-1,iqty))*(z1-z)/d   &
     &        -spl%fs1(0:n-1,iqty)*(6*z1-2)                             &
     &        +spl%fs1(1:n,iqty)*(6*z-2))/d                             
      ENDDO 
      IF(mode == 2)RETURN 
!-----------------------------------------------------------------------
!     evaluate third derivatives.                                       
!-----------------------------------------------------------------------
      DO iqty=1,nqty 
         f3(:,iqty)=(12*(spl%fs(0:n-1,iqty)-spl%fs(1:n,iqty))/d         &
     &        +6*(spl%fs1(0:n-1,iqty)+spl%fs1(1:n,iqty)))/(d*d)         
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE spline_all_eval 
!-----------------------------------------------------------------------
!     subprogram 6. spline_write.                                       
!     produces ascii and binary output for REAL cubic spline fits.      
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE spline_write(spl,out,bin,iua,iub) 
                                                                        
      TYPE(spline_type), INTENT(INOUT) :: spl 
      LOGICAL, INTENT(IN) :: out,bin 
      INTEGER, INTENT(IN) :: iua,iub 
                                                                        
      CHARACTER(30) :: format1,format2 
      INTEGER(i4) :: i,j 
      REAL(r8) :: x,dx 
!-----------------------------------------------------------------------
!     formats.                                                          
!-----------------------------------------------------------------------
   10 FORMAT('(/4x,''i'',',i2.2,'(4x,a6,1x)/)') 
   20 FORMAT('(i5,1p,',i2.2,'e11.3)') 
   30 FORMAT('(/4x,''i'',2x,''j'',',i2.2,'(4x,a6,1x)/)') 
!-----------------------------------------------------------------------
!     print ascii tables of node values and derivatives.                
!-----------------------------------------------------------------------
      IF(.NOT.out.AND..NOT.bin)RETURN 
      IF(out)THEN 
         WRITE(format1,10)spl%nqty+1 
         WRITE(format2,20)spl%nqty+1 
         WRITE(iua,'(/1x,a)')'input values:' 
         WRITE(iua,format1)spl%title(0:spl%nqty) 
         WRITE(iua,format2)(i,spl%xs(i),spl%fs(i,1:spl%nqty),           &
     &        i=0,spl%nodes)                                            
         WRITE(iua,format1)spl%title(0:spl%nqty) 
         WRITE(iua,'(/1x,a)')'derivatives:' 
         WRITE(iua,format1)spl%title(0:spl%nqty) 
         WRITE(iua,format2)(i,spl%xs(i),spl%fs1(i,1:spl%nqty),          &
     &        i=0,spl%nodes)                                            
         WRITE(iua,format1)spl%title(0:spl%nqty) 
      ENDIF 
!-----------------------------------------------------------------------
!     print binary table of node values.                                
!-----------------------------------------------------------------------
      IF(bin)THEN 
         DO i=0,spl%nodes 
            WRITE(iub)(/REAL(spl%xs(i),4),REAL(spl%fs(i,1:spl%nqty),4)/) 
         ENDDO 
         WRITE(iub) 
      ENDIF 
!-----------------------------------------------------------------------
!     print header for interpolated values.                             
!-----------------------------------------------------------------------
      IF(out)THEN 
         WRITE(iua,'(/1x,a)')'interpolated values:' 
         WRITE(iua,format1)spl%title(0:spl%nqty) 
      ENDIF 
!-----------------------------------------------------------------------
!     print interpolated values.                                        
!-----------------------------------------------------------------------
      DO i=0,spl%nodes-1 
         dx=(spl%xs(i+1)-spl%xs(i))/4 
         DO j=0,4 
            x=spl%xs(i)+j*dx 
            CALL spline_eval(spl,x,0_i4) 
            IF(out)WRITE(iua,format2)i,x,spl%f(1:spl%nqty) 
            IF(bin)WRITE(iub)(/REAL(x,4),REAL(spl%f(1:spl%nqty),4)/) 
         ENDDO 
         IF(out)WRITE(iua,'(1x)') 
      ENDDO 
!-----------------------------------------------------------------------
!     print final interpolated values.                                  
!-----------------------------------------------------------------------
      x=spl%xs(spl%nodes) 
      CALL spline_eval(spl,x,0_i4) 
      IF(out)THEN 
         WRITE(iua,format2)i,x,spl%f(1:spl%nqty) 
         WRITE(iua,format1)spl%title(0:spl%nqty) 
      ENDIF 
      IF(bin)THEN 
         WRITE(iub)(/REAL(x,4),REAL(spl%f(1:spl%nqty),4)/) 
         WRITE(iub) 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE spline_write 
!-----------------------------------------------------------------------
!     subprogram 7. spline_int.                                         
!     integrates REAL cubic splines.                                    
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE spline_int(spl) 
                                                                        
      TYPE(spline_type), INTENT(INOUT) :: spl 
                                                                        
      INTEGER(i4) :: i,iqty 
      REAL(r8) :: d 
!-----------------------------------------------------------------------
!     integrate cubic splines.                                          
!-----------------------------------------------------------------------
      IF(.NOT.ASSOCIATED(spl%fsi))                                      &
     &     ALLOCATE(spl%fsi(0:spl%nodes,spl%nqty))                      
      DO iqty=1,spl%nqty 
         spl%fsi(0,iqty)=0 
         DO i=1,spl%nodes 
            d=spl%xs(i)-spl%xs(i-1) 
            spl%fsi(i,iqty)=spl%fsi(i-1,iqty)                           &
     &           +d/12*(6*(spl%fs(i-1,iqty)+spl%fs(i,iqty))             &
     &           +d*(spl%fs1(i-1,iqty)-spl%fs1(i,iqty)))                
         ENDDO 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE spline_int 
!-----------------------------------------------------------------------
!     subprogram 9. spline_fac.                                         
!     sets up matrix for cubic spline fitting.                          
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE spline_fac(spl,a,b,cl,cr,endmode) 
                                                                        
      TYPE(spline_type), INTENT(IN) :: spl 
      REAL(r8), DIMENSION(-1:1,0:spl%nodes), INTENT(OUT) :: a 
      REAL(r8), DIMENSION(spl%nodes), INTENT(OUT) :: b 
      REAL(r8), DIMENSION(4), INTENT(OUT) :: cl,cr 
      CHARACTER(*), INTENT(IN) :: endmode 
                                                                        
      INTEGER(i4) :: j 
!-----------------------------------------------------------------------
!     compute interior matrix.                                          
!-----------------------------------------------------------------------
      b=1/(spl%xs(1:spl%nodes)-spl%xs(0:spl%nodes-1)) 
      DO j=1,spl%nodes-1 
         a(-1,j)=b(j) 
         a(0,j)=2*(b(j)+b(j+1)) 
         a(1,j)=b(j+1) 
      ENDDO 
!-----------------------------------------------------------------------
!     extrapolation boundary conditions.                                
!-----------------------------------------------------------------------
      IF(endmode == "extrap")THEN 
         b=b*b 
         cl(1)=(spl%xs(0)*(3*spl%xs(0)                                  &
     &        -2*(spl%xs(1)+spl%xs(2)+spl%xs(3)))                       &
     &        +spl%xs(1)*spl%xs(2)+spl%xs(1)*spl%xs(3)                  &
     &        +spl%xs(2)*spl%xs(3))                                     &
     &        /((spl%xs(0)-spl%xs(1))*(spl%xs(0)-spl%xs(2))             &
     &        *(spl%xs(0)-spl%xs(3)))                                   
         cl(2)=((spl%xs(2)-spl%xs(0))*(spl%xs(3)-spl%xs(0)))            &
     &        /((spl%xs(1)-spl%xs(0))*(spl%xs(1)-spl%xs(2))             &
     &        *(spl%xs(1)-spl%xs(3)))                                   
         cl(3)=((spl%xs(0)-spl%xs(1))*(spl%xs(3)-spl%xs(0)))            &
     &        /((spl%xs(0)-spl%xs(2))*(spl%xs(1)-spl%xs(2))             &
     &        *(spl%xs(3)-spl%xs(2)))                                   
         cl(4)=((spl%xs(1)-spl%xs(0))*(spl%xs(2)-spl%xs(0)))            &
     &        /((spl%xs(3)-spl%xs(0))*(spl%xs(3)-spl%xs(1))             &
     &        *(spl%xs(3)-spl%xs(2)))                                   
         cr(1)=(spl%xs(spl%nodes)*(3*spl%xs(spl%nodes)                  &
     &        -2*(spl%xs(spl%nodes-1)+spl%xs(spl%nodes-2)               &
     &        +spl%xs(spl%nodes-3)))+spl%xs(spl%nodes-1)                &
     &        *spl%xs(spl%nodes-2)+spl%xs(spl%nodes-1)*spl%xs(spl%nodes &
     &        -3)+spl%xs(spl%nodes-2)*spl%xs(spl%nodes-3))              &
     &        /((spl%xs(spl%nodes)-spl%xs(spl%nodes-1))                 &
     &        *(spl%xs(spl%nodes)-spl%xs(spl%nodes-2))*(spl%xs(spl%nodes&
     &        )-spl%xs(spl%nodes-3)))                                   
         cr(2)=((spl%xs(spl%nodes-2)-spl%xs(spl%nodes))                 &
     &        *(spl%xs(spl%nodes-3)-spl%xs(spl%nodes)))                 &
     &        /((spl%xs(spl%nodes-1)-spl%xs(spl%nodes))                 &
     &        *(spl%xs(spl%nodes-1)-spl%xs(spl%nodes-2))                &
     &        *(spl%xs(spl%nodes-1)-spl%xs(spl%nodes-3)))               
         cr(3)=((spl%xs(spl%nodes)-spl%xs(spl%nodes-1))                 &
     &        *(spl%xs(spl%nodes-3)-spl%xs(spl%nodes)))                 &
     &        /((spl%xs(spl%nodes)-spl%xs(spl%nodes-2))                 &
     &        *(spl%xs(spl%nodes-1)-spl%xs(spl%nodes-2))                &
     &        *(spl%xs(spl%nodes-3)-spl%xs(spl%nodes-2)))               
         cr(4)=((spl%xs(spl%nodes-1)-spl%xs(spl%nodes))                 &
     &        *(spl%xs(spl%nodes-2)-spl%xs(spl%nodes)))                 &
     &        /((spl%xs(spl%nodes-3)-spl%xs(spl%nodes))                 &
     &        *(spl%xs(spl%nodes-3)-spl%xs(spl%nodes-1))                &
     &        *(spl%xs(spl%nodes-3)-spl%xs(spl%nodes-2)))               
         CALL spline_triluf(a(:,1:),spl%nodes-1_i4) 
!-----------------------------------------------------------------------
!     not-a-knot boundary conditions.                                   
!-----------------------------------------------------------------------
      ELSE IF(endmode == "not-a-knot")THEN 
         b=b*b 
         a(0,1)=a(0,1)+(spl%xs(2)+spl%xs(0)-2*spl%xs(1))*b(1) 
         a(1,1)=a(1,1)+(spl%xs(2)-spl%xs(1))*b(1) 
         a(0,spl%nodes-1)=a(0,spl%nodes-1)                              &
     &        +(2*spl%xs(spl%nodes-1)-spl%xs(spl%nodes-2)               &
     &        -spl%xs(spl%nodes))*b(spl%nodes)                          
         a(-1,spl%nodes-1)=a(-1,spl%nodes-1)                            &
     &        +(spl%xs(spl%nodes-1)-spl%xs(spl%nodes-2))*b(spl%nodes)   
         CALL spline_triluf(a(:,1:),spl%nodes-1_i4) 
!-----------------------------------------------------------------------
!     periodic boundary conditions.                                     
!-----------------------------------------------------------------------
      ELSE IF(endmode == "periodic")THEN 
         a(0,0:spl%nodes:spl%nodes)=2*(b(spl%nodes)+b(1)) 
         a(1,0)=b(1) 
         a(-1,0)=b(spl%nodes) 
         b=b*b 
         CALL spline_sherman(a,spl%nodes) 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE spline_fac 
!-----------------------------------------------------------------------
!     subprogram 9. spline_triluf.                                      
!     performs tridiagonal LU factorization.                            
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE spline_triluf(a,n) 
                                                                        
      INTEGER(i4), INTENT(IN) :: n 
      REAL(r8), DIMENSION(-1:1,n), INTENT(INOUT) :: a 
                                                                        
      INTEGER(i4) :: i,j,k,jmin,jmax 
!-----------------------------------------------------------------------
!     begin loop over rows and define limits.                           
!-----------------------------------------------------------------------
      DO i=1,n 
         jmin=MAX(1_i4-i,-1_i4) 
         jmax=MIN(n-i,1_i4) 
!-----------------------------------------------------------------------
!     compute lower elements.                                           
!-----------------------------------------------------------------------
         DO j=jmin,-1 
            DO k=MAX(jmin,j-1_i4),j-1_i4 
               a(j,i)=a(j,i)-a(k,i)*a(j-k,i+k) 
            ENDDO 
            a(j,i)=a(j,i)*a(0,i+j) 
         ENDDO 
!-----------------------------------------------------------------------
!     compute diagonal element                                          
!-----------------------------------------------------------------------
         DO k=MAX(jmin,-1_i4),-1 
            a(0,i)=a(0,i)-a(k,i)*a(-k,i+k) 
         ENDDO 
         a(0,i)=1/a(0,i) 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE spline_triluf 
!-----------------------------------------------------------------------
!     subprogram 10. spline_trilus.                                     
!     performs tridiagonal LU solution.                                 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE spline_trilus(a,x,n,ns) 
                                                                        
      INTEGER(i4), INTENT(IN) :: n,ns 
      REAL(r8), DIMENSION(-1:1,n), INTENT(IN) :: a 
      REAL(r8), DIMENSION(n,ns), INTENT(INOUT) :: x 
                                                                        
      INTEGER(i4) :: i,j,is 
!-----------------------------------------------------------------------
!     down sweep.                                                       
!-----------------------------------------------------------------------
      DO is=1,ns 
         DO i=1,n 
            DO j=MAX(1_i4-i,-1_i4),-1 
               x(i,is)=x(i,is)-a(j,i)*x(i+j,is) 
            ENDDO 
         ENDDO 
!-----------------------------------------------------------------------
!     up sweep.                                                         
!-----------------------------------------------------------------------
         DO i=n,1,-1 
            DO j=1,MIN(n-i,1_i4) 
               x(i,is)=x(i,is)-a(j,i)*x(i+j,is) 
            ENDDO 
            x(i,is)=x(i,is)*a(0,i) 
         ENDDO 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE spline_trilus 
!-----------------------------------------------------------------------
!     subprogram 11. spline_sherman.                                    
!     uses Sherman-Morrison formula to factor periodic matrix.          
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE spline_sherman(a,n) 
                                                                        
      INTEGER(i4) :: n 
      REAL(r8), DIMENSION(-1:1,n), INTENT(INOUT) :: a 
                                                                        
      INTEGER(i4) :: j 
      REAL(r8), DIMENSION(n,1) :: u 
!-----------------------------------------------------------------------
!     prepare matrices.                                                 
!-----------------------------------------------------------------------
      a(0,1)=a(0,1)-a(-1,1) 
      a(0,n)=a(0,n)-a(-1,1) 
      u=RESHAPE((/1._r8,(0._r8,j=2,n-1),1._r8/),SHAPE(u)) 
      CALL spline_triluf(a,n) 
      CALL spline_trilus(a,u,n,1_i4) 
      a(-1,1)=a(-1,1)/(1+a(-1,1)*(u(1,1)+u(n,1))) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE spline_sherman 
!-----------------------------------------------------------------------
!     subprogram 12. spline_morrison.                                   
!     uses Sherman-Morrison formula to solve periodic matrix.           
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE spline_morrison(a,x,n,ns) 
                                                                        
      INTEGER(i4) :: n,ns 
      REAL(r8), DIMENSION(-1:1,n), INTENT(IN) :: a 
      REAL(r8), DIMENSION(n,ns), INTENT(INOUT) :: x 
                                                                        
      REAL(r8), DIMENSION(n,ns) :: y 
!-----------------------------------------------------------------------
!     solve for x.                                                      
!-----------------------------------------------------------------------
      y=x 
      CALL spline_trilus(a,y,n,ns) 
      x(1,:)=x(1,:)-a(-1,1)*(y(1,:)+y(n,:)) 
      x(n,:)=x(n,:)-a(-1,1)*(y(1,:)+y(n,:)) 
      CALL spline_trilus(a,x,n,ns) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE spline_morrison 
      END MODULE spline 
