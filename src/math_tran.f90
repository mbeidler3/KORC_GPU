!-----------------------------------------------------------------------
!     $Id: math_tran.f90 4517 2015-07-29 23:08:48Z jking $
!     module containing mathematical operations and coordinate          
!     transformations.                                                  
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.                                                
!-----------------------------------------------------------------------
!     1.  math_cart_cross_rr                                            
!     2.  math_cart_cross_rc                                            
!     3.  math_cart_cross_cr                                            
!     4.  math_cart_cross_cc                                            
!     5.  math_cadd_cross_rr                                            
!     6.  math_cadd_cross_rc                                            
!     7.  math_cadd_cross_cr                                            
!     8.  math_cadd_cross_cc                                            
!     10. math_inv_2x2                                                  
!     11. math_solve_3x3                                                
!     12. math_solve_nxn                                                
!     13. math_solve_q1_nxn                                             
!     14. math_solve_q1_cnxn                                            
!     15. math_solve_sym                                                
!     16. math_solve_q1_sym                                             
!     17. math_solve_herm                                               
!     18. math_solve_q1_herm                                            
!     19. math_grid                                                     
!     20. math_curl                                                     
!     21. math_grad                                                     
!-----------------------------------------------------------------------
!     delimit module.                                                   
!-----------------------------------------------------------------------
      MODULE math_tran 
      USE local 
      IMPLICIT NONE 
                                                                        
      INTERFACE math_cart_cross 
        MODULE PROCEDURE math_cart_cross_rr,math_cart_cross_rc,         &
     &                   math_cart_cross_cr,math_cart_cross_cc          
      END INTERFACE 
                                                                        
      INTERFACE math_cadd_cross 
        MODULE PROCEDURE math_cadd_cross_rr,math_cadd_cross_rc,         &
     &                   math_cadd_cross_cr,math_cadd_cross_cc          
      END INTERFACE 
                                                                        
      CONTAINS 
!-----------------------------------------------------------------------
!     subprogram 1. math_cart_cross_rr.                                 
!     evaluates the Cartesian cross product multiplied by a scalar      
!     over 2D arrays (both reals).                                      
!-----------------------------------------------------------------------
      SUBROUTINE math_cart_cross_rr(vprod,v1,v2,scal) 
                                                                        
        REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: vprod 
        REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: v1,v2 
        REAL(r8), INTENT(IN) :: scal 
!-----------------------------------------------------------------------
!     _x,_y, and _z or _r,_z, and _phi                                  
!-----------------------------------------------------------------------
      vprod(1,:,:)=scal*(v1(2,:,:)*v2(3,:,:)-v1(3,:,:)*v2(2,:,:)) 
      vprod(2,:,:)=scal*(v1(3,:,:)*v2(1,:,:)-v1(1,:,:)*v2(3,:,:)) 
      vprod(3,:,:)=scal*(v1(1,:,:)*v2(2,:,:)-v1(2,:,:)*v2(1,:,:)) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_cart_cross_rr 
!-----------------------------------------------------------------------
!     subprogram 2. math_cart_cross_rc.                                 
!     evaluates the Cartesian cross product multiplied by a real scalar 
!     over 2D arrays (one real and one complex).                        
!-----------------------------------------------------------------------
      SUBROUTINE math_cart_cross_rc(vprod,v1,v2,scal) 
                                                                        
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: vprod 
        REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: v1 
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: v2 
        REAL(r8), INTENT(IN) :: scal 
!-----------------------------------------------------------------------
!     _x,_y, and _z or _r,_z, and _phi                                  
!-----------------------------------------------------------------------
      vprod(1,:,:)=scal*(v1(2,:,:)*v2(3,:,:)-v1(3,:,:)*v2(2,:,:)) 
      vprod(2,:,:)=scal*(v1(3,:,:)*v2(1,:,:)-v1(1,:,:)*v2(3,:,:)) 
      vprod(3,:,:)=scal*(v1(1,:,:)*v2(2,:,:)-v1(2,:,:)*v2(1,:,:)) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_cart_cross_rc 
!-----------------------------------------------------------------------
!     subprogram 3. math_cart_cross_cr.                                 
!     evaluates the Cartesian cross product multiplied by a real scalar 
!     over 2D arrays (one real and one complex).                        
!-----------------------------------------------------------------------
      SUBROUTINE math_cart_cross_cr(vprod,v1,v2,scal) 
                                                                        
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: vprod 
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: v1 
        REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: v2 
        REAL(r8), INTENT(IN) :: scal 
!-----------------------------------------------------------------------
!     _x,_y, and _z or _r,_z, and _phi                                  
!-----------------------------------------------------------------------
      vprod(1,:,:)=scal*(v1(2,:,:)*v2(3,:,:)-v1(3,:,:)*v2(2,:,:)) 
      vprod(2,:,:)=scal*(v1(3,:,:)*v2(1,:,:)-v1(1,:,:)*v2(3,:,:)) 
      vprod(3,:,:)=scal*(v1(1,:,:)*v2(2,:,:)-v1(2,:,:)*v2(1,:,:)) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_cart_cross_cr 
!-----------------------------------------------------------------------
!     subprogram 4. math_cart_cross_cc.                                 
!     evaluates the Cartesian cross product multiplied by a real scalar 
!     over 2D arrays (both complex).                                    
!-----------------------------------------------------------------------
      SUBROUTINE math_cart_cross_cc(vprod,v1,v2,scal) 
                                                                        
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: vprod 
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: v1,v2 
        REAL(r8), INTENT(IN) :: scal 
!-----------------------------------------------------------------------
!     _x,_y, and _z or _r,_z, and _phi                                  
!-----------------------------------------------------------------------
      vprod(1,:,:)=scal*(v1(2,:,:)*v2(3,:,:)-v1(3,:,:)*v2(2,:,:)) 
      vprod(2,:,:)=scal*(v1(3,:,:)*v2(1,:,:)-v1(1,:,:)*v2(3,:,:)) 
      vprod(3,:,:)=scal*(v1(1,:,:)*v2(2,:,:)-v1(2,:,:)*v2(1,:,:)) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_cart_cross_cc 
!-----------------------------------------------------------------------
!     subprogram 5. math_cadd_cross_rr.                                 
!     adds the Cartesian cross product multiplied by a scalar           
!     over 2D arrays (both reals) to an array.                          
!-----------------------------------------------------------------------
      SUBROUTINE math_cadd_cross_rr(vprod,v1,v2,scal) 
                                                                        
        REAL(r8), DIMENSION(:,:,:), INTENT(INOUT) :: vprod 
        REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: v1,v2 
        REAL(r8), INTENT(IN) :: scal 
!-----------------------------------------------------------------------
!     _x,_y, and _z or _r,_z, and _phi                                  
!-----------------------------------------------------------------------
      vprod(1,:,:)=vprod(1,:,:)+                                        &
     &             scal*(v1(2,:,:)*v2(3,:,:)-v1(3,:,:)*v2(2,:,:))       
      vprod(2,:,:)=vprod(2,:,:)+                                        &
     &             scal*(v1(3,:,:)*v2(1,:,:)-v1(1,:,:)*v2(3,:,:))       
      vprod(3,:,:)=vprod(3,:,:)+                                        &
     &             scal*(v1(1,:,:)*v2(2,:,:)-v1(2,:,:)*v2(1,:,:))       
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_cadd_cross_rr 
!-----------------------------------------------------------------------
!     subprogram 6. math_cadd_cross_rc.                                 
!     adds the Cartesian cross product multiplied by a real scalar      
!     over 2D arrays (one real and one complex) to an array.            
!-----------------------------------------------------------------------
      SUBROUTINE math_cadd_cross_rc(vprod,v1,v2,scal) 
                                                                        
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(INOUT) :: vprod 
        REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: v1 
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: v2 
        REAL(r8), INTENT(IN) :: scal 
!-----------------------------------------------------------------------
!     _x,_y, and _z or _r,_z, and _phi                                  
!-----------------------------------------------------------------------
      vprod(1,:,:)=vprod(1,:,:)+                                        &
     &             scal*(v1(2,:,:)*v2(3,:,:)-v1(3,:,:)*v2(2,:,:))       
      vprod(2,:,:)=vprod(2,:,:)+                                        &
     &             scal*(v1(3,:,:)*v2(1,:,:)-v1(1,:,:)*v2(3,:,:))       
      vprod(3,:,:)=vprod(3,:,:)+                                        &
     &             scal*(v1(1,:,:)*v2(2,:,:)-v1(2,:,:)*v2(1,:,:))       
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_cadd_cross_rc 
!-----------------------------------------------------------------------
!     subprogram 7. math_cadd_cross_cr.                                 
!     adds the Cartesian cross product multiplied by a real scalar      
!     over 2D arrays (one real and one complex) to an array.            
!-----------------------------------------------------------------------
      SUBROUTINE math_cadd_cross_cr(vprod,v1,v2,scal) 
                                                                        
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(INOUT) :: vprod 
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: v1 
        REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: v2 
        REAL(r8), INTENT(IN) :: scal 
!-----------------------------------------------------------------------
!     _x,_y, and _z or _r,_z, and _phi                                  
!-----------------------------------------------------------------------
      vprod(1,:,:)=vprod(1,:,:)+                                        &
     &             scal*(v1(2,:,:)*v2(3,:,:)-v1(3,:,:)*v2(2,:,:))       
      vprod(2,:,:)=vprod(2,:,:)+                                        &
     &             scal*(v1(3,:,:)*v2(1,:,:)-v1(1,:,:)*v2(3,:,:))       
      vprod(3,:,:)=vprod(3,:,:)+                                        &
     &             scal*(v1(1,:,:)*v2(2,:,:)-v1(2,:,:)*v2(1,:,:))       
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_cadd_cross_cr 
!-----------------------------------------------------------------------
!     subprogram 8. math_cadd_cross_cc.                                 
!     adds the Cartesian cross product multiplied by a real scalar      
!     over 2D arrays (both complex) to an array.                        
!-----------------------------------------------------------------------
      SUBROUTINE math_cadd_cross_cc(vprod,v1,v2,scal) 
                                                                        
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(INOUT) :: vprod 
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: v1,v2 
        REAL(r8), INTENT(IN) :: scal 
!-----------------------------------------------------------------------
!     _x,_y, and _z or _r,_z, and _phi                                  
!-----------------------------------------------------------------------
      vprod(1,:,:)=vprod(1,:,:)+                                        &
     &             scal*(v1(2,:,:)*v2(3,:,:)-v1(3,:,:)*v2(2,:,:))       
      vprod(2,:,:)=vprod(2,:,:)+                                        &
     &             scal*(v1(3,:,:)*v2(1,:,:)-v1(1,:,:)*v2(3,:,:))       
      vprod(3,:,:)=vprod(3,:,:)+                                        &
     &             scal*(v1(1,:,:)*v2(2,:,:)-v1(2,:,:)*v2(1,:,:))       
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_cadd_cross_cc 
!-----------------------------------------------------------------------
!     subprogram 10.  math_inv_2x2                                      
!     invert a 2D array of 2x2 matrices.                                
!-----------------------------------------------------------------------
      SUBROUTINE math_inv_2x2(out11,out12,out21,out22,                  &
     &                         in11, in12, in21, in22)                  
                                                                        
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: out11,out12,out21,out22 
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: in11,in12,in21,in22 
                                                                        
      REAL(r8), DIMENSION(SIZE(out11,1),SIZE(out11,2)) :: det 
!-----------------------------------------------------------------------
!     find determinant.                                                 
!-----------------------------------------------------------------------
      det=in11*in22-in12*in21 
!-----------------------------------------------------------------------
!     invert in-matrix.                                                 
!-----------------------------------------------------------------------
      out11= in22/det 
      out12=-in12/det 
      out21=-in21/det 
      out22= in11/det 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_inv_2x2 
!-----------------------------------------------------------------------
!     subprogram 11. math_solve_3x3.                                    
!     use lu factorization to solve a 3x3 system.                       
!-----------------------------------------------------------------------
      SUBROUTINE math_solve_3x3(matrix,xx,bb,instruction) 
                                                                        
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: matrix 
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: bb 
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: xx 
      CHARACTER(6), INTENT(IN) :: instruction 
                                                                        
      REAL(r8), DIMENSION(SIZE(matrix,1),SIZE(matrix,2)) :: y2,y3,      &
     &          lu21,lu22,lu23,lu31,lu32,lu33                           
!-----------------------------------------------------------------------
!     factor matrix into L*U (first row is identical to matrix).        
!-----------------------------------------------------------------------
      IF (instruction(1:5)/='solve') THEN 
        lu21 = matrix(:,:,2,1)/matrix(:,:,1,1) 
        lu31 = matrix(:,:,3,1)/matrix(:,:,1,1) 
        lu22 = matrix(:,:,2,2)-lu21*matrix(:,:,1,2) 
        lu23 = matrix(:,:,2,3)-lu21*matrix(:,:,1,3) 
        lu32 =(matrix(:,:,3,2)-lu31*matrix(:,:,1,2))/lu22 
        lu33 = matrix(:,:,3,3)-lu31*matrix(:,:,1,3)-lu32*lu23 
      ENDIF 
      IF (instruction=='factor') THEN 
        matrix(:,:,2,1)=lu21 
        matrix(:,:,3,1)=lu31 
        matrix(:,:,2,2)=lu22 
        matrix(:,:,3,2)=lu32 
        matrix(:,:,2,3)=lu23 
        matrix(:,:,3,3)=lu33 
        RETURN 
      ENDIF 
!-----------------------------------------------------------------------
!     solve for xx by getting y first.                                  
!     {Ly=bb then Uxx=y}  Note that L has ones on the diagonal.         
!-----------------------------------------------------------------------
      IF (instruction(1:5)=='solve') THEN 
        lu21=matrix(:,:,2,1) 
        lu31=matrix(:,:,3,1) 
        lu22=matrix(:,:,2,2) 
        lu32=matrix(:,:,3,2) 
        lu23=matrix(:,:,2,3) 
        lu33=matrix(:,:,3,3) 
      ENDIF 
      y2 = bb(:,:,2)-lu21*bb(:,:,1) 
      y3 = bb(:,:,3)-lu31*bb(:,:,1)-lu32*y2 
      xx(:,:,3) = y3/lu33 
      xx(:,:,2) =(y2-lu23*xx(:,:,3))/lu22 
      xx(:,:,1) =(bb(:,:,1)-matrix(:,:,1,2)*xx(:,:,2)                   &
     &                     -matrix(:,:,1,3)*xx(:,:,3))/matrix(:,:,1,1)  
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_solve_3x3 
!-----------------------------------------------------------------------
!     subprogram 12. math_solve_nxn.                                    
!     use lu factorization to solve a nxn system.  note that L has ones 
!     on the diagonal.  if the instruction is 'factor', the             
!     decomposition is returned in matrix.  if the instruction is       
!     'solve', the incoming matrix is expected to be the decomposition. 
!     if the instruction is 'both', matrix is not changed, but the      
!     decomposition is not saved after the solution vector is found.    
!-----------------------------------------------------------------------
      SUBROUTINE math_solve_nxn(n,matrix,xx,bb,instruction) 
                                                                        
      INTEGER(i4), INTENT(IN) :: n 
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: matrix 
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: bb 
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: xx 
      CHARACTER(*), INTENT(IN) :: instruction 
                                                                        
      REAL(r8), DIMENSION(SIZE(matrix,1),SIZE(matrix,2),n) :: y 
      REAL(r8), DIMENSION(SIZE(matrix,1),SIZE(matrix,2),n,n) :: lu 
      INTEGER(i4) :: isub,j 
!-----------------------------------------------------------------------
!     factor matrix into L*U (first row and column are treated first).  
!-----------------------------------------------------------------------
      IF (instruction/='solve') THEN 
        lu(:,:,1,:)=matrix(:,:,1,:) 
        DO j=2,n 
          lu(:,:,j,1)=matrix(:,:,j,1)/lu(:,:,1,1) 
        ENDDO 
!-----------------------------------------------------------------------
!       submatrix loop: for each find the row of L then the column of U.
!-----------------------------------------------------------------------
        DO isub=2,n 
          DO j=2,isub-1 
            lu(:,:,isub,j)=(matrix(:,:,isub,j)                          &
     &               -SUM(lu(:,:,isub,1:j-1)*lu(:,:,1:j-1,j),3))        &
     &                   /lu(:,:,j,j)                                   
          ENDDO 
          DO j=2,isub 
            lu(:,:,j,isub)= matrix(:,:,j,isub)                          &
     &               -SUM(lu(:,:,j,1:j-1)*lu(:,:,1:j-1,isub),3)         
          ENDDO 
        ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     copy matrix and return if solution is not required.               
!-----------------------------------------------------------------------
      IF (instruction=='factor') THEN 
        matrix=lu 
        RETURN 
      ENDIF 
!-----------------------------------------------------------------------
!     copy LU decomposition if already factored.                        
!-----------------------------------------------------------------------
      IF (instruction=='solve') THEN 
        lu=matrix 
      ENDIF 
!-----------------------------------------------------------------------
!     solve for xx by getting y first.  {Ly=bb then Uxx=y}              
!-----------------------------------------------------------------------
      y(:,:,1)=bb(:,:,1) 
      DO j=2,n 
        y(:,:,j)=bb(:,:,j)-SUM(lu(:,:,j,1:j-1)*y(:,:,1:j-1),3) 
      ENDDO 
      xx(:,:,n)=y(:,:,n)/lu(:,:,n,n) 
      DO j=n-1,1,-1 
        xx(:,:,j)=(y(:,:,j)-SUM(lu(:,:,j,j+1:n)*xx(:,:,j+1:n),3))       &
     &                         /lu(:,:,j,j)                             
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_solve_nxn 
!-----------------------------------------------------------------------
!     subprogram 13. math_solve_q1_nxn.                                 
!     use lu factorization to solve a real nxn system.  note that L has 
!     ones on the diagonal.  if the instruction is 'factor', the        
!     decomposition is returned in matrix.  if the instruction is       
!     'solve', the incoming matrix is expected to be the decomposition. 
!     if the instruction is 'both', matrix is not changed, but the      
!     decomposition is not saved after the solution vector is found.    
!                                                                       
!     this version has the quantity indices first and an assumed size   
!     index last.  also, the column index is first and the row index    
!     is second.                                                        
!-----------------------------------------------------------------------
      SUBROUTINE math_solve_q1_nxn(n,na,matrix,xx,bb,instruction,       &
     &                             singular)                            
                                                                        
      INTEGER(i4), INTENT(IN) :: n,na 
      REAL(r8), DIMENSION(n,n,*), INTENT(INOUT) :: matrix 
      REAL(r8), DIMENSION(n,*), INTENT(IN) :: bb 
      REAL(r8), DIMENSION(n,*), INTENT(OUT) :: xx 
      CHARACTER(*), INTENT(IN) :: instruction 
      LOGICAL, INTENT(OUT) :: singular 
                                                                        
      REAL(r8), DIMENSION(n) :: y 
      REAL(r8), DIMENSION(n,n) :: lu 
      INTEGER(i4) :: isub,i,ia 
!-----------------------------------------------------------------------
!     factor matrix into L*U (first row and column are treated first).  
!-----------------------------------------------------------------------
      singular=.false. 
      DO ia=1,na 
        IF (instruction/='solve') THEN 
          lu=matrix(:,:,ia) 
          IF (lu(1,1)==0) THEN 
            singular=.true. 
            RETURN 
          ENDIF 
          DO i=2,n 
            lu(1,i)=lu(1,i)/lu(1,1) 
          ENDDO 
!-----------------------------------------------------------------------
!         submatrix loop: for each find the row of L then the column    
!         of U.                                                         
!-----------------------------------------------------------------------
          DO isub=2,n 
            DO i=2,isub-1 
              lu(i,isub)=(lu(i,isub)                                    &
     &                  -SUM(lu(1:i-1,isub)*lu(i,1:i-1)))/lu(i,i)       
            ENDDO 
            DO i=2,isub 
              lu(isub,i)=lu(isub,i)-SUM(lu(1:i-1,i)*lu(isub,1:i-1)) 
            ENDDO 
            IF (lu(isub,isub)==0) THEN 
              singular=.true. 
              RETURN 
            ENDIF 
          ENDDO 
        ENDIF 
!-----------------------------------------------------------------------
!       copy matrix and return if solution is not required.             
!-----------------------------------------------------------------------
        IF (instruction=='factor') THEN 
          matrix(:,:,ia)=lu 
          CYCLE 
        ENDIF 
!-----------------------------------------------------------------------
!       copy LU decomposition if already factored.                      
!-----------------------------------------------------------------------
        IF (instruction=='solve') THEN 
          lu=matrix(:,:,ia) 
        ENDIF 
!-----------------------------------------------------------------------
!       solve for xx by getting y first.  {Ly=bb then Uxx=y}            
!-----------------------------------------------------------------------
        y(1)=bb(1,ia) 
        DO i=2,n 
          y(i)=bb(i,ia)-SUM(lu(1:i-1,i)*y(1:i-1)) 
        ENDDO 
        xx(n,ia)=y(n)/lu(n,n) 
        DO i=n-1,1,-1 
          xx(i,ia)=(y(i)-SUM(lu(i+1:n,i)*xx(i+1:n,ia)))/lu(i,i) 
        ENDDO 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_solve_q1_nxn 
!-----------------------------------------------------------------------
!     subprogram 14. math_solve_q1_cnxn.                                
!     use lu factorization to solve a complex nxn system.  note that L  
!     has ones on the diagonal.  if the instruction is 'factor', the    
!     decomposition is returned in matrix.  if the instruction is       
!     'solve', the incoming matrix is expected to be the decomposition. 
!     if the instruction is 'both', matrix is not changed, but the      
!     decomposition is not saved after the solution vector is found.    
!                                                                       
!     this version has the quantity indices first and an assumed size   
!     index last.  also, the column index is first and the row index    
!     is second.  furthermore, the matrix is complex.                   
!-----------------------------------------------------------------------
      SUBROUTINE math_solve_q1_cnxn(n,na,matrix,xx,bb,instruction,      &
     &                              singular)                           
                                                                        
      INTEGER(i4), INTENT(IN) :: n,na 
      COMPLEX(r8), DIMENSION(n,n,*), INTENT(INOUT) :: matrix 
      COMPLEX(r8), DIMENSION(n,*), INTENT(IN) :: bb 
      COMPLEX(r8), DIMENSION(n,*), INTENT(OUT) :: xx 
      CHARACTER(*), INTENT(IN) :: instruction 
      LOGICAL, INTENT(OUT) :: singular 
                                                                        
      COMPLEX(r8), DIMENSION(n) :: y 
      COMPLEX(r8), DIMENSION(n,n) :: lu 
      INTEGER(i4) :: isub,i,ia 
!-----------------------------------------------------------------------
!     factor matrix into L*U (first row and column are treated first).  
!-----------------------------------------------------------------------
      singular=.false. 
      DO ia=1,na 
        IF (instruction/='solve') THEN 
          lu=matrix(:,:,ia) 
          IF (lu(1,1)==0) THEN 
            singular=.true. 
            RETURN 
          ENDIF 
          DO i=2,n 
            lu(1,i)=lu(1,i)/lu(1,1) 
          ENDDO 
!-----------------------------------------------------------------------
!         submatrix loop: for each find the row of L then the column    
!         of U.                                                         
!-----------------------------------------------------------------------
          DO isub=2,n 
            DO i=2,isub-1 
              lu(i,isub)=(lu(i,isub)                                    &
     &                  -SUM(lu(1:i-1,isub)*lu(i,1:i-1)))/lu(i,i)       
            ENDDO 
            DO i=2,isub 
              lu(isub,i)=lu(isub,i)-SUM(lu(1:i-1,i)*lu(isub,1:i-1)) 
            ENDDO 
            IF (lu(isub,isub)==0) THEN 
              singular=.true. 
              RETURN 
            ENDIF 
          ENDDO 
        ENDIF 
!-----------------------------------------------------------------------
!       copy matrix and return if solution is not required.             
!-----------------------------------------------------------------------
        IF (instruction=='factor') THEN 
          matrix(:,:,ia)=lu 
          CYCLE 
        ENDIF 
!-----------------------------------------------------------------------
!       copy LU decomposition if already factored.                      
!-----------------------------------------------------------------------
        IF (instruction=='solve') THEN 
          lu=matrix(:,:,ia) 
        ENDIF 
!-----------------------------------------------------------------------
!       solve for xx by getting y first.  {Ly=bb then Uxx=y}            
!-----------------------------------------------------------------------
        y(1)=bb(1,ia) 
        DO i=2,n 
          y(i)=bb(i,ia)-SUM(lu(1:i-1,i)*y(1:i-1)) 
        ENDDO 
        xx(n,ia)=y(n)/lu(n,n) 
        DO i=n-1,1,-1 
          xx(i,ia)=(y(i)-SUM(lu(i+1:n,i)*xx(i+1:n,ia)))/lu(i,i) 
        ENDDO 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_solve_q1_cnxn 
!-----------------------------------------------------------------------
!     subprogram 15. math_solve_sym.                                    
!     use Cholesky factorization to solve a symmetric nxn system.       
!     if the instruction is 'factor', the lower triangular part of the  
!     decomposition is returned in matrix.  if the instruction is       
!     'solve', the incoming matrix is expected to be the decomposition. 
!     if the instruction is 'both', matrix is not changed, but the      
!     decomposition is not saved after the solution vector is found.    
!-----------------------------------------------------------------------
      SUBROUTINE math_solve_sym(n,matrix,xx,bb,instruction,singular) 
                                                                        
      INTEGER(i4), INTENT(IN) :: n 
      REAL(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: matrix 
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: bb 
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: xx 
      CHARACTER(*), INTENT(IN) :: instruction 
      LOGICAL, INTENT(OUT) :: singular 
                                                                        
      REAL(r8), DIMENSION(SIZE(matrix,1),SIZE(matrix,2),n) :: y 
      REAL(r8), DIMENSION(SIZE(matrix,1),SIZE(matrix,2),n,n) :: lu 
      INTEGER(i4) :: jcol,i 
!-----------------------------------------------------------------------
!     factor matrix into L*transpose(L) (first column is treated first).
!     if a singularity is encountered, return the pivot in the iq=jq=1  
!     position of matrix.                                               
!-----------------------------------------------------------------------
      singular=.false. 
      IF (instruction/='solve') THEN 
        lu=0 
        IF (MINVAL(matrix(:,:,1,1))<=0) THEN 
          singular=.true. 
          RETURN 
        ENDIF 
        lu(:,:,1,1)=SQRT(matrix(:,:,1,1)) 
        DO i=2,n 
          lu(:,:,i,1)=matrix(:,:,i,1)/lu(:,:,1,1) 
        ENDDO 
!-----------------------------------------------------------------------
!       column loop for L.                                              
!-----------------------------------------------------------------------
        DO jcol=2,n 
          lu(:,:,jcol,jcol)=matrix(:,:,jcol,jcol)                       &
     &          -SUM(lu(:,:,jcol,1:jcol-1)*lu(:,:,jcol,1:jcol-1),3)     
          IF (MINVAL(lu(:,:,jcol,jcol))<=0) THEN 
            singular=.true. 
            matrix(:,:,1,1)=lu(:,:,jcol,jcol) 
            RETURN 
          ENDIF 
          lu(:,:,jcol,jcol)=SQRT(lu(:,:,jcol,jcol)) 
          DO i=jcol+1,n 
            lu(:,:,i,jcol)=(matrix(:,:,i,jcol)                          &
     &          -SUM(lu(:,:,   i,1:jcol-1)*lu(:,:,jcol,1:jcol-1),3))    &
     &              /lu(:,:,jcol,jcol)                                  
          ENDDO 
        ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     copy matrix and return if solution is not required.               
!-----------------------------------------------------------------------
      IF (instruction=='factor') THEN 
        matrix=lu 
        RETURN 
      ENDIF 
!-----------------------------------------------------------------------
!     copy decomposition if already factored.                           
!-----------------------------------------------------------------------
      IF (instruction=='solve') THEN 
        lu=matrix 
      ENDIF 
!-----------------------------------------------------------------------
!     solve for xx by getting y first.  {Ly=bb then transpose(L)xx=y}   
!-----------------------------------------------------------------------
      y(:,:,1)=bb(:,:,1)/lu(:,:,1,1) 
      DO i=2,n 
        y(:,:,i)=(bb(:,:,i)-SUM(lu(:,:,i,1:i-1)*y(:,:,1:i-1),3))        &
     &                         /lu(:,:,i,i)                             
      ENDDO 
      xx(:,:,n)=y(:,:,n)/lu(:,:,n,n) 
      DO i=n-1,1,-1 
        xx(:,:,i)=(y(:,:,i)-SUM(lu(:,:,i+1:n,i)*xx(:,:,i+1:n),3))       &
     &                         /lu(:,:,i,i)                             
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_solve_sym 
!-----------------------------------------------------------------------
!     subprogram 16. math_solve_q1_sym.                                 
!     use Cholesky factorization to solve a symmetric nxn system.       
!     if the instruction is 'factor', the lower triangular part of the  
!     decomposition is returned in matrix.  if the instruction is       
!     'solve', the incoming matrix is expected to be the decomposition. 
!     if the instruction is 'both', matrix is not changed, but the      
!     decomposition is not saved after the solution vector is found.    
!                                                                       
!     this version has the quantity indices first and an assumed size   
!     index last.  also, the column index is first and the row index    
!     is second.                                                        
!-----------------------------------------------------------------------
      SUBROUTINE math_solve_q1_sym(n,na,matrix,xx,bb,instruction,       &
     &                             singular)                            
                                                                        
      INTEGER(i4), INTENT(IN) :: n,na 
      REAL(r8), DIMENSION(n,n,*), INTENT(INOUT) :: matrix 
      REAL(r8), DIMENSION(n,*), INTENT(IN) :: bb 
      REAL(r8), DIMENSION(n,*), INTENT(OUT) :: xx 
      CHARACTER(*), INTENT(IN) :: instruction 
      LOGICAL, INTENT(OUT) :: singular 
                                                                        
      REAL(r8), DIMENSION(n) :: y 
      REAL(r8), DIMENSION(n,n) :: lu 
      INTEGER(i4) :: jcol,i,ia 
!-----------------------------------------------------------------------
!     factor matrix into L*transpose(L) (first column is treated first).
!     if a singularity is encountered, return the pivot in the iq=jq=1  
!     position of matrix.                                               
!-----------------------------------------------------------------------
      singular=.false. 
      DO ia=1,na 
        IF (instruction=='solve') THEN 
          lu=matrix(:,:,ia) 
        ELSE 
          lu=0 
          IF (matrix(1,1,ia)<=0) THEN 
            singular=.true. 
            RETURN 
          ENDIF 
          lu(1,1)=SQRT(matrix(1,1,ia)) 
          DO i=2,n 
            lu(1,i)=matrix(1,i,ia)/lu(1,1) 
          ENDDO 
!-----------------------------------------------------------------------
!         column loop for L.                                            
!-----------------------------------------------------------------------
          DO jcol=2,n 
            lu(jcol,jcol)=matrix(jcol,jcol,ia)                          &
     &            -SUM(lu(1:jcol-1,jcol)*lu(1:jcol-1,jcol))             
            IF (lu(jcol,jcol)<=0) THEN 
              singular=.true. 
              matrix(1,1,ia)=lu(jcol,jcol) 
              RETURN 
            ENDIF 
            lu(jcol,jcol)=SQRT(lu(jcol,jcol)) 
            DO i=jcol+1,n 
              lu(jcol,i)=(matrix(jcol,i,ia)                             &
     &            -SUM(lu(1:jcol-1,i)*lu(1:jcol-1,jcol)))               &
     &                /lu(jcol,jcol)                                    
            ENDDO 
          ENDDO 
        ENDIF 
!-----------------------------------------------------------------------
!       copy matrix if solution is not required.                        
!-----------------------------------------------------------------------
        IF (instruction=='factor') THEN 
          matrix(:,:,ia)=lu 
          CYCLE 
        ENDIF 
!-----------------------------------------------------------------------
!       solve for xx by getting y first.  {Ly=bb then transpose(L)xx=y} 
!-----------------------------------------------------------------------
        y(1)=bb(1,ia)/lu(1,1) 
        DO i=2,n 
          y(i)=(bb(i,ia)-SUM(lu(1:i-1,i)*y(1:i-1)))/lu(i,i) 
        ENDDO 
        xx(n,ia)=y(n)/lu(n,n) 
        DO i=n-1,1,-1 
          xx(i,ia)=(y(i)-SUM(lu(i,i+1:n)*xx(i+1:n,ia)))/lu(i,i) 
        ENDDO 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_solve_q1_sym 
!-----------------------------------------------------------------------
!     subprogram 17. math_solve_herm.                                   
!     use Cholesky factorization to solve an Hermitian nxn system.      
!     if the instruction is 'factor', the lower triangular part of the  
!     decomposition is returned in matrix.  if the instruction is       
!     'solve', the incoming matrix is expected to be the decomposition. 
!     if the instruction is 'both', matrix is not changed, but the      
!     decomposition is not saved after the solution vector is found.    
!-----------------------------------------------------------------------
      SUBROUTINE math_solve_herm(n,matrix,xx,bb,instruction,singular) 
                                                                        
      INTEGER(i4), INTENT(IN) :: n 
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(INOUT) :: matrix 
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: bb 
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: xx 
      CHARACTER(*), INTENT(IN) :: instruction 
      LOGICAL, INTENT(OUT) :: singular 
                                                                        
      COMPLEX(r8), DIMENSION(SIZE(matrix,1),SIZE(matrix,2),n) :: y 
      COMPLEX(r8), DIMENSION(SIZE(matrix,1),SIZE(matrix,2),n,n) :: lu 
      INTEGER(i4) :: jcol,i 
!-----------------------------------------------------------------------
!     factor matrix into L*transpose(L^*) (first column is treated      
!     first).  if a singularity is encountered, return the pivot in     
!     the iq=jq=1 position of matrix.                                   
!-----------------------------------------------------------------------
      singular=.false. 
      IF (instruction/='solve') THEN 
        lu=0 
        IF (MINVAL(REAL(matrix(:,:,1,1),r8))<=0) THEN 
          singular=.true. 
          RETURN 
        ENDIF 
        lu(:,:,1,1)=SQRT(REAL(matrix(:,:,1,1),r8)) 
        DO i=2,n 
          lu(:,:,i,1)=matrix(:,:,i,1)/lu(:,:,1,1) 
        ENDDO 
!-----------------------------------------------------------------------
!       column loop for L.                                              
!-----------------------------------------------------------------------
        DO jcol=2,n 
          lu(:,:,jcol,jcol)=matrix(:,:,jcol,jcol)                       &
     &      -SUM(lu(:,:,jcol,1:jcol-1)*CONJG(lu(:,:,jcol,1:jcol-1)),3)  
          IF (MINVAL(REAL(lu(:,:,jcol,jcol),r8))<=0) THEN 
            singular=.true. 
            matrix(:,:,1,1)=lu(:,:,jcol,jcol) 
            RETURN 
          ENDIF 
          lu(:,:,jcol,jcol)=SQRT(REAL(lu(:,:,jcol,jcol),r8)) 
          DO i=jcol+1,n 
            lu(:,:,i,jcol)=(matrix(:,:,i,jcol)                          &
     &          -SUM(       lu(:,:,   i,1:jcol-1)                       &
     &               *CONJG(lu(:,:,jcol,1:jcol-1)),3))                  &
     &              /lu(:,:,jcol,jcol)                                  
          ENDDO 
        ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     copy matrix and return if solution is not required.               
!-----------------------------------------------------------------------
      IF (instruction=='factor') THEN 
        matrix=lu 
        RETURN 
      ENDIF 
!-----------------------------------------------------------------------
!     copy decomposition if already factored.                           
!-----------------------------------------------------------------------
      IF (instruction=='solve') THEN 
        lu=matrix 
      ENDIF 
!-----------------------------------------------------------------------
!     solve for xx by getting y first.  {Ly=bb then transpose(L^*)xx=y} 
!-----------------------------------------------------------------------
      y(:,:,1)=bb(:,:,1)/lu(:,:,1,1) 
      DO i=2,n 
        y(:,:,i)=(bb(:,:,i)-SUM(lu(:,:,i,1:i-1)*y(:,:,1:i-1),3))        &
     &                         /lu(:,:,i,i)                             
      ENDDO 
      xx(:,:,n)=y(:,:,n)/lu(:,:,n,n) 
      DO i=n-1,1,-1 
        xx(:,:,i)=(y(:,:,i)-SUM(CONJG(lu(:,:,i+1:n,i))*xx(:,:,i+1:n),3))&
     &                               /lu(:,:,i,i)                       
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_solve_herm 
!-----------------------------------------------------------------------
!     subprogram 18. math_solve_q1_herm.                                
!     use Cholesky factorization to solve an Hermitian nxn system.      
!     if the instruction is 'factor', the lower triangular part of the  
!     decomposition is returned in matrix.  if the instruction is       
!     'solve', the incoming matrix is expected to be the decomposition. 
!     if the instruction is 'both', matrix is not changed, but the      
!     decomposition is not saved after the solution vector is found.    
!                                                                       
!     this version has the quantity indices first and an assumed size   
!     index last.  also, the column index is first and the row index    
!     is second.                                                        
!-----------------------------------------------------------------------
      SUBROUTINE math_solve_q1_herm(n,na,matrix,xx,bb,instruction,      &
     &                              singular)                           
                                                                        
      INTEGER(i4), INTENT(IN) :: n,na 
      COMPLEX(r8), DIMENSION(n,n,*), INTENT(INOUT) :: matrix 
      COMPLEX(r8), DIMENSION(n,*), INTENT(IN) :: bb 
      COMPLEX(r8), DIMENSION(n,*), INTENT(OUT) :: xx 
      CHARACTER(*), INTENT(IN) :: instruction 
      LOGICAL, INTENT(OUT) :: singular 
                                                                        
      COMPLEX(r8), DIMENSION(n) :: y 
      COMPLEX(r8), DIMENSION(n,n) :: lu 
      INTEGER(i4) :: jcol,i,ia 
!-----------------------------------------------------------------------
!     factor matrix into L*transpose(L^*) (first column is treated      
!     first).  if a singularity is encountered, return the pivot in     
!     the iq=jq=1 position of matrix.                                   
!-----------------------------------------------------------------------
      singular=.false. 
      DO ia=1,na 
        IF (instruction=='solve') THEN 
          lu=matrix(:,:,ia) 
        ELSE 
          lu=0 
          IF (REAL(matrix(1,1,ia),r8)<=0) THEN 
            singular=.true. 
            RETURN 
          ENDIF 
          lu(1,1)=SQRT(REAL(matrix(1,1,ia),r8)) 
          DO i=2,n 
            lu(1,i)=matrix(1,i,ia)/lu(1,1) 
          ENDDO 
!-----------------------------------------------------------------------
!         column loop for L.                                            
!-----------------------------------------------------------------------
          DO jcol=2,n 
            lu(jcol,jcol)=matrix(jcol,jcol,ia)                          &
     &        -SUM(lu(1:jcol-1,jcol)*CONJG(lu(1:jcol-1,jcol)))          
            IF (REAL(lu(jcol,jcol),r8)<=0) THEN 
              singular=.true. 
              matrix(1,1,ia)=lu(jcol,jcol) 
              RETURN 
            ENDIF 
            lu(jcol,jcol)=SQRT(REAL(lu(jcol,jcol),r8)) 
            DO i=jcol+1,n 
              lu(jcol,i)=(matrix(jcol,i,ia)                             &
     &            -SUM(lu(1:jcol-1,i)*CONJG(lu(1:jcol-1,jcol))))        &
     &                /lu(jcol,jcol)                                    
            ENDDO 
          ENDDO 
        ENDIF 
!-----------------------------------------------------------------------
!       copy matrix if solution is not required.                        
!-----------------------------------------------------------------------
        IF (instruction=='factor') THEN 
          matrix(:,:,ia)=lu 
          CYCLE 
        ENDIF 
!-----------------------------------------------------------------------
!       solve for xx by getting y first. {Ly=bb then transpose(L^*)xx=y}
!-----------------------------------------------------------------------
        y(1)=bb(1,ia)/lu(1,1) 
        DO i=2,n 
          y(i)=(bb(i,ia)-SUM(lu(1:i-1,i)*y(1:i-1)))/lu(i,i) 
        ENDDO 
        xx(n,ia)=y(n)/lu(n,n) 
        DO i=n-1,1,-1 
          xx(i,ia)=(y(i)-SUM(CONJG(lu(i,i+1:n))*xx(i+1:n,ia)))/lu(i,i) 
        ENDDO 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_solve_q1_herm 
!-----------------------------------------------------------------------
!     subprogram 19. math_grid.                                         
!     find metric-related derivatives of r and z.                       
!-----------------------------------------------------------------------
      SUBROUTINE math_grid(option,drzdx,drzdy,jac,dxdr,dxdz,dydr,dydz) 
                                                                        
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: drzdx,drzdy 
      CHARACTER(*), INTENT(IN) :: option 
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: jac,dxdr,dxdz,dydr,dydz 
                                                                        
!-----------------------------------------------------------------------
!     find the jacobian, then perform a matrix inversion.               
!-----------------------------------------------------------------------
      jac = drzdx(1,:,:)*drzdy(2,:,:)-drzdx(2,:,:)*drzdy(1,:,:) 
      IF (option=='jaco') RETURN 
      dxdr= drzdy(2,:,:)/jac 
      dydr=-drzdx(2,:,:)/jac 
      dxdz=-drzdy(1,:,:)/jac 
      dydz= drzdx(1,:,:)/jac 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_grid 
!-----------------------------------------------------------------------
!     subprogram 20. math_curl                                          
!     find the cylindrical curl of the incoming vector (with cylindrical
!     phi component) where derivatives with respect to r and z are      
!     supplied.  this version uses complex arrays.                      
!-----------------------------------------------------------------------
      SUBROUTINE math_curl(nmodes,keff,geom,bigr,vec,dvecr,dvecz,       &
     &                     curl,scal)                                   
                                                                        
      INTEGER(i4), INTENT(IN) :: nmodes 
      REAL(r8), INTENT(IN) :: scal 
      REAL(r8), DIMENSION(:), INTENT(IN) :: keff 
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr 
      CHARACTER(*), INTENT(IN) :: geom 
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(IN) :: vec,dvecr,dvecz 
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: curl 
                                                                        
      INTEGER(i4) :: im 
!-----------------------------------------------------------------------
!     note:  vector storage is (comp1,comp2,comp3)                      
!     where (1,2,3)=(x,y,z) or (r,z,phi).                               
!                                                                       
!     separate coding is now required for toroidal and linear           
!     geometries since phi component is cylindrical not covariant.      
!-----------------------------------------------------------------------
      IF (geom=='tor') THEN 
        DO im=1,nmodes 
          curl(1,:,:,im)= scal*(                                        &
     &       dvecz(3,:,:,im)-(0,1)*keff(im)*vec(2,:,:,im) /bigr )       
          curl(2,:,:,im)= scal*(-dvecr(3,:,:,im)                        &
     &      +(-vec(3,:,:,im)+(0,1)*keff(im)*vec(1,:,:,im))/bigr )       
          curl(3,:,:,im)= scal*( dvecr(2,:,:,im)-dvecz(1,:,:,im) ) 
        ENDDO 
      ELSE 
        DO im=1,nmodes 
          curl(1,:,:,im)= scal*(                                        &
     &       dvecz(3,:,:,im)-(0,1)*keff(im)*vec(2,:,:,im) )             
          curl(2,:,:,im)= scal*(                                        &
     &      -dvecr(3,:,:,im)+(0,1)*keff(im)*vec(1,:,:,im) )             
          curl(3,:,:,im)= scal*( dvecr(2,:,:,im)-dvecz(1,:,:,im) ) 
        ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_curl 
!-----------------------------------------------------------------------
!     subprogram 21. math_grad                                          
!     assemble a field and its r and z derivatives into storage for a   
!     gradient where the index cycles over the partial-derivative       
!     directions fastest.                                               
!                                                                       
!     additional terms for 3-vectors in cylindrical/toroidal geometry   
!     are added, and an n=0 component is optional.                      
!-----------------------------------------------------------------------
      SUBROUTINE math_grad(nmodes,keff,nq,geom,vec,dvecr,dvecz,         &
     &                     grad,bigr,vec0,dvec0r,dvec0z)                
                                                                        
      INTEGER(i4), INTENT(IN) :: nmodes,nq 
      REAL(r8), DIMENSION(:), INTENT(IN) :: keff 
      REAL(r8), DIMENSION(:,:), INTENT(IN), OPTIONAL :: bigr 
      REAL(r8), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: vec0,         &
     &          dvec0r,dvec0z                                           
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(IN) :: vec,dvecr,dvecz 
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: grad 
      CHARACTER(*), INTENT(IN) :: geom 
                                                                        
      INTEGER(i4) :: im,iq,igr 
      LOGICAL :: add_n0 
!-----------------------------------------------------------------------
!     note:  incoming vector storage is (comp1,comp2,comp3)             
!     where (1,2,3)=(x,y,z) or (r,z,phi).  the outgoing gradient        
!     storage is (dcomp1/dr,dcomp1/dz,...)                              
!-----------------------------------------------------------------------
      IF (PRESENT(vec0)) THEN 
        add_n0=.true. 
      ELSE 
        add_n0=.false. 
      ENDIF 
!-----------------------------------------------------------------------
!     nq-specific cases are intended for optimization.                  
!-----------------------------------------------------------------------
      IF (geom=='tor') THEN 
        SELECT CASE(nq) 
        CASE(1) 
          DO im=1,nmodes 
            IF (keff(im)==0.AND.add_n0) THEN 
              grad(1,:,:,im)=dvecr(1,:,:,im)+dvec0r(1,:,:) 
              grad(2,:,:,im)=dvecz(1,:,:,im)+dvec0z(1,:,:) 
              grad(3,:,:,im)=0 
            ELSE 
              grad(1,:,:,im)=dvecr(1,:,:,im) 
              grad(2,:,:,im)=dvecz(1,:,:,im) 
              grad(3,:,:,im)=  vec(1,:,:,im)*(0,1)*keff(im)/bigr 
            ENDIF 
          ENDDO 
        CASE(3) !  Including curvature terms.                          
          DO im=1,nmodes 
            IF (keff(im)==0.AND.add_n0) THEN 
              grad(1,:,:,im)=dvecr(1,:,:,im)+dvec0r(1,:,:) 
              grad(2,:,:,im)=dvecz(1,:,:,im)+dvec0z(1,:,:) 
              grad(3,:,:,im)=-(vec(3,:,:,im)+  vec0(3,:,:))/bigr 
              grad(4,:,:,im)=dvecr(2,:,:,im)+dvec0r(2,:,:) 
              grad(5,:,:,im)=dvecz(2,:,:,im)+dvec0z(2,:,:) 
              grad(6,:,:,im)=0 
              grad(7,:,:,im)=dvecr(3,:,:,im)+dvec0r(3,:,:) 
              grad(8,:,:,im)=dvecz(3,:,:,im)+dvec0z(3,:,:) 
              grad(9,:,:,im)= (vec(1,:,:,im)+  vec0(1,:,:))/bigr 
            ELSE 
              grad(1,:,:,im)=dvecr(1,:,:,im) 
              grad(2,:,:,im)=dvecz(1,:,:,im) 
              grad(3,:,:,im)= (vec(1,:,:,im)*(0,1)*keff(im)-            &
     &                         vec(3,:,:,im))/bigr                      
              grad(4,:,:,im)=dvecr(2,:,:,im) 
              grad(5,:,:,im)=dvecz(2,:,:,im) 
              grad(6,:,:,im)=  vec(2,:,:,im)*(0,1)*keff(im)/bigr 
              grad(7,:,:,im)=dvecr(3,:,:,im) 
              grad(8,:,:,im)=dvecz(3,:,:,im) 
              grad(9,:,:,im)= (vec(3,:,:,im)*(0,1)*keff(im)+            &
     &                         vec(1,:,:,im))/bigr                      
            ENDIF 
          ENDDO 
        CASE DEFAULT 
          DO im=1,nmodes 
            IF (keff(im)==0.AND.add_n0) THEN 
              igr=1 
              DO iq=1,nq 
                grad(igr  ,:,:,im)=dvecr(iq,:,:,im)+dvec0r(iq,:,:) 
                grad(igr+1,:,:,im)=dvecz(iq,:,:,im)+dvec0z(iq,:,:) 
                grad(igr+2,:,:,im)=0 
                igr=igr+3 
              ENDDO 
            ELSE 
              igr=1 
              DO iq=1,nq 
                grad(igr  ,:,:,im)=dvecr(iq,:,:,im) 
                grad(igr+1,:,:,im)=dvecz(iq,:,:,im) 
                grad(igr+2,:,:,im)=  vec(iq,:,:,im)*(0,1)*keff(im)/bigr 
                igr=igr+3 
              ENDDO 
            ENDIF 
          ENDDO 
        END SELECT 
!-----------------------------------------------------------------------
!     same for linear geometry.                                         
!-----------------------------------------------------------------------
      ELSE !  geom='lin'                                               
        SELECT CASE(nq) 
        CASE(1) 
          DO im=1,nmodes 
            IF (keff(im)==0.AND.add_n0) THEN 
              grad(1,:,:,im)=dvecr(1,:,:,im)+dvec0r(1,:,:) 
              grad(2,:,:,im)=dvecz(1,:,:,im)+dvec0z(1,:,:) 
              grad(3,:,:,im)=0 
            ELSE 
              grad(1,:,:,im)=dvecr(1,:,:,im) 
              grad(2,:,:,im)=dvecz(1,:,:,im) 
              grad(3,:,:,im)=  vec(1,:,:,im)*(0,1)*keff(im) 
            ENDIF 
          ENDDO 
        CASE(3) 
          DO im=1,nmodes 
            IF (keff(im)==0.AND.add_n0) THEN 
              grad(1,:,:,im)=dvecr(1,:,:,im)+dvec0r(1,:,:) 
              grad(2,:,:,im)=dvecz(1,:,:,im)+dvec0z(1,:,:) 
              grad(3,:,:,im)=0 
              grad(4,:,:,im)=dvecr(2,:,:,im)+dvec0r(2,:,:) 
              grad(5,:,:,im)=dvecz(2,:,:,im)+dvec0z(2,:,:) 
              grad(6,:,:,im)=0 
              grad(7,:,:,im)=dvecr(3,:,:,im)+dvec0r(3,:,:) 
              grad(8,:,:,im)=dvecz(3,:,:,im)+dvec0z(3,:,:) 
              grad(9,:,:,im)=0 
            ELSE 
              grad(1,:,:,im)=dvecr(1,:,:,im) 
              grad(2,:,:,im)=dvecz(1,:,:,im) 
              grad(3,:,:,im)=  vec(1,:,:,im)*(0,1)*keff(im) 
              grad(4,:,:,im)=dvecr(2,:,:,im) 
              grad(5,:,:,im)=dvecz(2,:,:,im) 
              grad(6,:,:,im)=  vec(2,:,:,im)*(0,1)*keff(im) 
              grad(7,:,:,im)=dvecr(3,:,:,im) 
              grad(8,:,:,im)=dvecz(3,:,:,im) 
              grad(9,:,:,im)=  vec(3,:,:,im)*(0,1)*keff(im) 
            ENDIF 
          ENDDO 
        CASE DEFAULT 
          DO im=1,nmodes 
            IF (keff(im)==0.AND.add_n0) THEN 
              igr=1 
              DO iq=1,nq 
                grad(igr  ,:,:,im)=dvecr(iq,:,:,im)+dvec0r(iq,:,:) 
                grad(igr+1,:,:,im)=dvecz(iq,:,:,im)+dvec0z(iq,:,:) 
                grad(igr+2,:,:,im)=0 
                igr=igr+3 
              ENDDO 
            ELSE 
              igr=1 
              DO iq=1,nq 
                grad(igr  ,:,:,im)=dvecr(iq,:,:,im) 
                grad(igr+1,:,:,im)=dvecz(iq,:,:,im) 
                grad(igr+2,:,:,im)=  vec(iq,:,:,im)*(0,1)*keff(im) 
                igr=igr+3 
              ENDDO 
            ENDIF 
          ENDDO 
        END SELECT 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE math_grad 
!-----------------------------------------------------------------------
!     close module.                                                     
!-----------------------------------------------------------------------
      END MODULE math_tran
