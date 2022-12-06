!-----------------------------------------------------------------------
!     $Id: tri_linear.f90 4517 2015-07-29 23:08:48Z jking $
!     fits functions to tri_linear splines.                             
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.                                                
!-----------------------------------------------------------------------
!     0.  tri_linear_type definition.                                   
!     1.  tri_linear_geom_alloc.                                        
!     2.  tri_linear_geom_dealloc.                                      
!     3.  tri_linear_get_areas.                                         
!     4.  tri_linear_3D_alloc.                                                                              
!     6.  tri_linear_3D_all_eval.                                       
!     7.  tri_linear_3D_eval.                                           
!     8.  tri_linear_2D_alloc.                                                                               
!     10. tri_linear_2D_all_eval.                                       
!     11. tri_linear_2D_eval.                                           
!     12. tri_linear_3D_assign_rsc.                                     
!     13. tri_linear_3D_assign_csc.                                     
!     14. tri_linear_3D_assign_tl3.                                     
!     15. tri_linear_3D_assign_int.                                     
!     16. tri_linear_2D_assign_rsc.                                     
!     17. tri_linear_2D_assign_csc.                                     
!     18. tri_linear_2D_assign_tl2.                                     
!     19. tri_linear_2D_assign_int.                                     
!     16. dump_write_tri_linear.                                        
!     17. dump_write_tri_linear_2D.                                     
!     18. dump_read_tri_linear.                                         
!     19. dump_read_tri_linear_2D.                                      
!-----------------------------------------------------------------------
!     subprogram 0. tri_linear_type definition.                         
!     defines tri_linear_type.                                          
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      MODULE tri_linear 
      USE local 
      IMPLICIT NONE 
                                                                        
      TYPE :: neighbor_type 
      INTEGER(i4), DIMENSION(:), POINTER :: vertex 
      END TYPE neighbor_type 
                                                                        
      TYPE :: tri_linear_geom_type 
        INTEGER(i4) :: mvert,mcell 
        INTEGER(i4), DIMENSION(:,:), POINTER :: vertex 
        INTEGER(i4), DIMENSION(:,:,:), POINTER :: segment 
        REAL(r8), DIMENSION(:), POINTER :: xs,ys,area 
        REAL(r8), DIMENSION(3,1:7) :: alpha,alphab,delta 
        REAL(r8), DIMENSION(:,:), POINTER :: wjac,bigr 
        REAL(r8), DIMENSION(:,:,:), POINTER :: alpha_x,alpha_y,         &
     &    alpha_arr,dalpdr,dalpdz,dalpdrc,dalpmdr,dalpmdz,alpham_arr,   &
     &    dalpmdrc                                                      
        TYPE(neighbor_type), DIMENSION(:), POINTER :: neighbor 
      END TYPE tri_linear_geom_type 
                                                                        
      TYPE :: tri_linear_type 
        CHARACTER(6) :: name 
        CHARACTER(6), DIMENSION(:), POINTER :: title 
        INTEGER(i4) :: mvert,nqty,nfour 
        COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: fs 
        COMPLEX(r8), DIMENSION(:,:), POINTER :: f,fx,fy 
      END TYPE tri_linear_type 
                                                                        
      TYPE :: tri_linear_2D_type 
        CHARACTER(6) :: name 
        CHARACTER(6), DIMENSION(:), POINTER :: title 
        INTEGER(i4) :: mvert,nqty 
        REAL(r8), DIMENSION(:,:,:), POINTER :: fs 
        REAL(r8), DIMENSION(:), POINTER :: f,fx,fy 
      END TYPE tri_linear_2D_type 
                                                                        
!-----------------------------------------------------------------------
!     subprogram name interfaces                                        
!-----------------------------------------------------------------------
      INTERFACE tri_linear_alloc 
        MODULE PROCEDURE tri_linear_2D_alloc,tri_linear_3D_alloc 
      END INTERFACE 
                                                                        
                                                                  
      INTERFACE tri_linear_all_eval 
        MODULE PROCEDURE tri_linear_2D_all_eval,tri_linear_3D_all_eval 
      END INTERFACE 
                                                                        
      INTERFACE tri_linear_eval 
        MODULE PROCEDURE tri_linear_2D_eval,tri_linear_3D_eval 
      END INTERFACE 
!-----------------------------------------------------------------------
!     overloaded assignment.                                            
!-----------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=) 
        MODULE PROCEDURE tri_linear_3D_assign_csc,                      &
     &    tri_linear_3D_assign_rsc,tri_linear_3D_assign_tl3,            &
     &    tri_linear_3D_assign_int,tri_linear_2D_assign_csc,            &
     &    tri_linear_2D_assign_rsc,tri_linear_2D_assign_tl2,            &
     &    tri_linear_2D_assign_int                                      
      END INTERFACE 
                                                                        
      CONTAINS 
!-----------------------------------------------------------------------
!     subprogram 1. tri_linear_geom_alloc.                              
!     allocates space for tri_linear geometry.                          
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_geom_alloc(v,mvert,mcell) 
                                                                        
      INTEGER(i4), INTENT(IN) :: mvert,mcell 
      TYPE(tri_linear_geom_type), INTENT(INOUT) :: v 
!-----------------------------------------------------------------------
!     allocate space.                                                   
!-----------------------------------------------------------------------
      v%mvert=mvert 
      v%mcell=mcell 
      ALLOCATE(v%xs(0:mvert)) 
      ALLOCATE(v%ys(0:mvert)) 
      ALLOCATE(v%vertex(mcell,3)) 
      ALLOCATE(v%alpha_x(mcell,1,3)) 
      ALLOCATE(v%alpha_y(mcell,1,3)) 
      ALLOCATE(v%neighbor(0:mvert)) 
      ALLOCATE(v%area(mcell)) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_geom_alloc 
!-----------------------------------------------------------------------
!     subprogram 2. tri_linear_geom_dealloc.                            
!     deallocates space for spline_type.                                
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_geom_dealloc(v) 
                                                                        
      TYPE(tri_linear_geom_type), INTENT(INOUT) :: v 
!-----------------------------------------------------------------------
!     allocate space.                                                   
!-----------------------------------------------------------------------
      DEALLOCATE(v%xs) 
      DEALLOCATE(v%ys) 
      DEALLOCATE(v%vertex) 
      DEALLOCATE(v%alpha_x) 
      DEALLOCATE(v%alpha_y) 
      DEALLOCATE(v%neighbor) 
      DEALLOCATE(v%area) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_geom_dealloc 
!-----------------------------------------------------------------------
!     subprogram 3. tri_linear_get_areas.                               
!     evaluates areas and basis functions of triangular cells.          
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_get_areas(v) 
                                                                        
      TYPE(tri_linear_geom_type), INTENT(INOUT) :: v 
                                                                        
      INTEGER(i4) :: icell,iqty 
      REAL(r8), PARAMETER ::                                            &
     &     alpha1=0.05971587178976982_r8,                               &
     &     alpha2=0.1012865073234563_r8,                                &
     &     alpha3=0.3333333333333333_r8,                                &
     &     alpha4=0.4701420641051151_r8,                                &
     &     alpha5=0.7974269853530873_r8                                 
!-----------------------------------------------------------------------
!     set weight functions.                                             
!-----------------------------------------------------------------------
      v%alpha=TRANSPOSE(RESHAPE((/                                      &
     &     alpha3,alpha1,alpha4,alpha4,alpha5,alpha2,alpha2,            &
     &     alpha3,alpha4,alpha1,alpha4,alpha2,alpha5,alpha2,            &
     &     alpha3,alpha4,alpha4,alpha1,alpha2,alpha2,alpha5             &
     &     /),(/7,3/)))                                                 
      v%alphab=alpha3 
      v%delta=TRANSPOSE(RESHAPE((/                                      &
     &     1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/),(/7,3/)))        
!-----------------------------------------------------------------------
!     compute areas.                                                    
!-----------------------------------------------------------------------
      DO icell=1,v%mcell 
         v%area(icell)=0.5_r8*                                          &
     &        (v%xs(v%vertex(icell,1))*v%ys(v%vertex(icell,2))          &
     &        -v%xs(v%vertex(icell,2))*v%ys(v%vertex(icell,1))          &
     &        +v%xs(v%vertex(icell,2))*v%ys(v%vertex(icell,3))          &
     &        -v%xs(v%vertex(icell,3))*v%ys(v%vertex(icell,2))          &
     &        +v%xs(v%vertex(icell,3))*v%ys(v%vertex(icell,1))          &
     &        -v%xs(v%vertex(icell,1))*v%ys(v%vertex(icell,3)))         
      ENDDO 
!-----------------------------------------------------------------------
!     compute factors for x derivatives.  the second array index is to  
!     match the basis array structure in quadrilateral blocks.          
!-----------------------------------------------------------------------
      DO icell=1,v%mcell 
         v%alpha_x(icell,1,1)=0.5_r8/v%area(icell)                      &
     &        *(v%ys(v%vertex(icell,2))-v%ys(v%vertex(icell,3)))        
         v%alpha_x(icell,1,2)=0.5_r8/v%area(icell)                      &
     &        *(v%ys(v%vertex(icell,3))-v%ys(v%vertex(icell,1)))        
         v%alpha_x(icell,1,3)=0.5_r8/v%area(icell)                      &
     &        *(v%ys(v%vertex(icell,1))-v%ys(v%vertex(icell,2)))        
      ENDDO 
!-----------------------------------------------------------------------
!     compute factors for y derivatives.                                
!-----------------------------------------------------------------------
      DO icell=1,v%mcell 
         v%alpha_y(icell,1,1)=-0.5_r8/v%area(icell)                     &
     &        *(v%xs(v%vertex(icell,2))-v%xs(v%vertex(icell,3)))        
         v%alpha_y(icell,1,2)=-0.5_r8/v%area(icell)                     &
     &        *(v%xs(v%vertex(icell,3))-v%xs(v%vertex(icell,1)))        
         v%alpha_y(icell,1,3)=-0.5_r8/v%area(icell)                     &
     &        *(v%xs(v%vertex(icell,1))-v%xs(v%vertex(icell,2)))        
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_get_areas 
!-----------------------------------------------------------------------
!     subprogram 4. tri_linear_3D_alloc.                                
!     allocates space for tri_linear dependent variables.               
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_alloc(u,mvert,nqty,nfour,name,title) 
                                                                        
      INTEGER(i4), INTENT(IN) :: mvert,nqty,nfour 
      TYPE(tri_linear_type), INTENT(INOUT) :: u 
      CHARACTER(*), INTENT(IN), OPTIONAL :: name 
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title 
!-----------------------------------------------------------------------
!     allocate space.                                                   
!-----------------------------------------------------------------------
      u%mvert=mvert 
      u%nqty=nqty 
      u%nfour=nfour 
      ALLOCATE(u%fs(nqty,0:mvert,0:0,nfour)) 
      ALLOCATE(u%title(nqty)) 
      ALLOCATE(u%f(nqty,nfour)) 
      ALLOCATE(u%fx(nqty,nfour)) 
      ALLOCATE(u%fy(nqty,nfour)) 
!-----------------------------------------------------------------------
!     character descriptors, if present in input.                       
!-----------------------------------------------------------------------
      IF (PRESENT(name)) u%name=name 
      IF (PRESENT(title)) THEN 
        IF (SIZE(title)==nqty) THEN 
          u%title=title 
        ELSE 
          u%title=title(1) 
        ENDIF 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_3D_alloc 
!-----------------------------------------------------------------------
!     subprogram 6. tri_linear_3D_all_eval.                             
!     evaluates bicubic splines in all intervals for equal spacing.     
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_all_eval(u,v,node,mode,f,fx,fy) 
      USE io
                                                                        
      TYPE(tri_linear_type), INTENT(IN) :: u 
      TYPE(tri_linear_geom_type), INTENT(IN) :: v 
      INTEGER(i4), INTENT(IN) :: node,mode 
      COMPLEX(r8), INTENT(OUT),DIMENSION(:,:,0:,:) :: f,fx,fy 
                                                                        
      INTEGER(i4) :: icell,iqty 
!-----------------------------------------------------------------------
!     verify consistency between u and v.                               
!-----------------------------------------------------------------------
      IF(u%mvert/=v%mvert)THEN 
         WRITE(nim_wr,'(a)')                                            &
     &      "tri_linear_all_eval: geometry and dependent",              &
     &      " variables have inconsistent sizes."                       
         STOP 
      ENDIF 
!-----------------------------------------------------------------------
!     interpolate functions.                                            
!-----------------------------------------------------------------------
      DO icell=1,v%mcell 
         f(:,icell,0,:)                                                 &
     &        =u%fs(:,v%vertex(icell,1),0,:)*v%alpha(1,node)            &
     &        +u%fs(:,v%vertex(icell,2),0,:)*v%alpha(2,node)            &
     &        +u%fs(:,v%vertex(icell,3),0,:)*v%alpha(3,node)            
      ENDDO 
      IF(mode==0)RETURN 
!-----------------------------------------------------------------------
!     compute derivatives.                                              
!-----------------------------------------------------------------------
      DO icell=1,v%mcell 
         fx(:,icell,0,:)                                                &
     &        =u%fs(:,v%vertex(icell,1),0,:)*v%alpha_x(icell,1,1)       &
     &        +u%fs(:,v%vertex(icell,2),0,:)*v%alpha_x(icell,1,2)       &
     &        +u%fs(:,v%vertex(icell,3),0,:)*v%alpha_x(icell,1,3)       
         fy(:,icell,0,:)                                                &
     &        =u%fs(:,v%vertex(icell,1),0,:)*v%alpha_y(icell,1,1)       &
     &        +u%fs(:,v%vertex(icell,2),0,:)*v%alpha_y(icell,1,2)       &
     &        +u%fs(:,v%vertex(icell,3),0,:)*v%alpha_y(icell,1,3)       
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_3D_all_eval 
!-----------------------------------------------------------------------
!     subprogram 7. tri_linear_3D_eval.                                 
!     linear interpolation at (x,y) of u base on the                    
!     three points associated with cell icell on grid v.                
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_eval(u,v,x,y,icell,mode) 
      USE io 
                                                                        
      REAL(r8), INTENT(IN) :: x,y 
      INTEGER(i4), INTENT(IN) :: icell,mode 
      TYPE(tri_linear_type), INTENT(INOUT) :: u 
      TYPE(tri_linear_geom_type), INTENT(IN) :: v 
                                                                        
      INTEGER(i4) :: iqty,ifour 
      TYPE :: position_type 
        REAL(r8) :: r,z 
      END TYPE position_type 
      TYPE(position_type) :: p,p1,p2,p3 
      REAL(r8) :: denom 
      COMPLEX(r8) :: fx,fy 
                                                                        
!-----------------------------------------------------------------------
!     verify consistency between u and v.                               
!-----------------------------------------------------------------------
      IF(u%mvert/=v%mvert)THEN 
         WRITE(nim_wr,'(a)')                                            &
     &     "tri_linear_eval: geometry and dependent",                   &
     &     " variables have inconsistent sizes."                        
         STOP 
      ENDIF 
!-----------------------------------------------------------------------
!     verify icell is a valid cell                                      
!-----------------------------------------------------------------------
      IF((icell < 0) .OR. (icell > v%mcell))THEN 
         CALL nim_stop("tri_linear_eval: invalid cell") 
      ENDIF 
!-----------------------------------------------------------------------
!     Define positions                                                  
!-----------------------------------------------------------------------
      p%r=x 
      p%z=y 
      p1%r=v%xs(v%vertex(icell,1)) 
      p1%z=v%ys(v%vertex(icell,1)) 
      p2%r=v%xs(v%vertex(icell,2)) 
      p2%z=v%ys(v%vertex(icell,2)) 
      p3%r=v%xs(v%vertex(icell,3)) 
      p3%z=v%ys(v%vertex(icell,3)) 
      denom=1.0_r8/((p2%r-p1%r)*(p3%z-p1%z)-(p3%r-p1%r)*(p2%z-p1%z)) 
!-----------------------------------------------------------------------
!     Evaluate derivatives and functions.                               
!-----------------------------------------------------------------------
      DO ifour=1,u%nfour 
        DO iqty=1,u%nqty 
          fx=denom*(                                                    &
     &            (u%fs(iqty,v%vertex(icell,2),0,ifour)-                &
     &             u%fs(iqty,v%vertex(icell,1),0,ifour))*(p3%z-p1%z)-   &
     &            (u%fs(iqty,v%vertex(icell,3),0,ifour)-                &
     &             u%fs(iqty,v%vertex(icell,1),0,ifour))*(p2%z-p1%z))   
          fy=denom*(                                                    &
     &            (u%fs(iqty,v%vertex(icell,3),0,ifour)-                &
     &             u%fs(iqty,v%vertex(icell,1),0,ifour))*(p2%r-p1%r)-   &
     &            (u%fs(iqty,v%vertex(icell,2),0,ifour)-                &
     &             u%fs(iqty,v%vertex(icell,1),0,ifour))*(p3%r-p1%r))   
          u%f(iqty,ifour)=                                              &
     &            u%fs(iqty,v%vertex(icell,1),0,ifour)+                 &
     &            fy*(p%z-p1%z)+                                        &
     &            fx*(p%r-p1%r)                                         
          IF(mode /= 0 )THEN 
            u%fx(iqty,ifour) = fx 
            u%fy(iqty,ifour) = fy 
          ENDIF 
        ENDDO 
      ENDDO 
                                                                        
      RETURN 
      END SUBROUTINE tri_linear_3D_eval 
!-----------------------------------------------------------------------
!     subprogram 8. tri_linear_2D_alloc.                                
!     allocates space for tri_linear dependent variables.               
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_alloc(u,mvert,nqty,name,title) 
                                                                        
      INTEGER(i4), INTENT(IN) :: mvert,nqty 
      TYPE(tri_linear_2D_type), INTENT(INOUT) :: u 
      CHARACTER(*), INTENT(IN), OPTIONAL :: name 
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title 
!-----------------------------------------------------------------------
!     allocate space.                                                   
!-----------------------------------------------------------------------
      u%mvert=mvert 
      u%nqty=nqty 
      ALLOCATE(u%fs(nqty,0:mvert,0:0)) 
      ALLOCATE(u%title(nqty)) 
      ALLOCATE(u%f(nqty)) 
      ALLOCATE(u%fx(nqty)) 
      ALLOCATE(u%fy(nqty)) 
!-----------------------------------------------------------------------
!     character descriptors, if present in input.                       
!-----------------------------------------------------------------------
      IF (PRESENT(name)) u%name=name 
      IF (PRESENT(title)) THEN 
        IF (SIZE(title)==nqty) THEN 
          u%title=title 
        ELSE 
          u%title=title(1) 
        ENDIF 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_2D_alloc 
!-----------------------------------------------------------------------
!     subprogram 10. tri_linear_2D_all_eval.                            
!     evaluates bicubic splines in all intervals for equal spacing.     
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_all_eval(u,v,node,mode,f,fx,fy) 
      USE io 
                                                                        
      TYPE(tri_linear_2D_type), INTENT(IN) :: u 
      TYPE(tri_linear_geom_type), INTENT(IN) :: v 
      INTEGER(i4), INTENT(IN) :: node,mode 
      REAL(r8), INTENT(OUT),DIMENSION(:,:,0:) :: f,fx,fy 
                                                                        
      INTEGER(i4) :: icell,iqty 
!-----------------------------------------------------------------------
!     verify consistency between u and v.                               
!-----------------------------------------------------------------------
      IF(u%mvert/=v%mvert)THEN 
         WRITE(nim_wr,'(a)')                                            &
     &     "tri_linear_all_eval: geometry and dependent",               &
     &     " variables have inconsistent sizes."                        
         STOP 
      ENDIF 
!-----------------------------------------------------------------------
!     interpolate functions.                                            
!-----------------------------------------------------------------------
      DO icell=1,v%mcell 
         f(:,icell,0)                                                   &
     &        =u%fs(:,v%vertex(icell,1),0)*v%alpha(1,node)              &
     &        +u%fs(:,v%vertex(icell,2),0)*v%alpha(2,node)              &
     &        +u%fs(:,v%vertex(icell,3),0)*v%alpha(3,node)              
      ENDDO 
      IF(mode==0)RETURN 
!-----------------------------------------------------------------------
!     compute derivatives.                                              
!-----------------------------------------------------------------------
      DO icell=1,v%mcell 
         fx(:,icell,0)                                                  &
     &        =u%fs(:,v%vertex(icell,1),0)*v%alpha_x(icell,1,1)         &
     &        +u%fs(:,v%vertex(icell,2),0)*v%alpha_x(icell,1,2)         &
     &        +u%fs(:,v%vertex(icell,3),0)*v%alpha_x(icell,1,3)         
         fy(:,icell,0)                                                  &
     &        =u%fs(:,v%vertex(icell,1),0)*v%alpha_y(icell,1,1)         &
     &        +u%fs(:,v%vertex(icell,2),0)*v%alpha_y(icell,1,2)         &
     &        +u%fs(:,v%vertex(icell,3),0)*v%alpha_y(icell,1,3)         
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_2D_all_eval 
!-----------------------------------------------------------------------
!     subprogram 11. tri_linear_2D_eval.                                
!     linear interpolation at (x,y) of u base on the                    
!     three points associated with cell icell on grid v.                
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_eval(u,v,x,y,icell,mode) 
      USE io 
                                                                        
      REAL(r8), INTENT(IN) :: x,y 
      INTEGER(i4), INTENT(IN) :: icell,mode 
      TYPE(tri_linear_2D_type), INTENT(INOUT) :: u 
      TYPE(tri_linear_geom_type), INTENT(IN) :: v 
                                                                        
      INTEGER(i4) :: iqty 
      TYPE :: position_type 
        REAL(r8) :: r,z 
      END TYPE position_type 
      TYPE(position_type) :: p,p1,p2,p3 
      REAL(r8) :: denom,fx,fy 
                                                                        
!-----------------------------------------------------------------------
!     verify consistency between u and v.                               
!-----------------------------------------------------------------------
      IF(u%mvert/=v%mvert)THEN 
         WRITE(nim_wr,'(a)')                                            &
     &     "tri_linear_eval: geometry and dependent",                   &
     &     " variables have inconsistent sizes."                        
         STOP 
      ENDIF 
!-----------------------------------------------------------------------
!     verify icell is a valid cell                                      
!-----------------------------------------------------------------------
      IF((icell < 0) .OR. (icell > v%mcell))THEN 
         CALL nim_stop("tri_linear_eval: invalid cell") 
      ENDIF 
!-----------------------------------------------------------------------
!     Define positions                                                  
!-----------------------------------------------------------------------
      p%r=x 
      p%z=y 
      p1%r=v%xs(v%vertex(icell,1)) 
      p1%z=v%ys(v%vertex(icell,1)) 
      p2%r=v%xs(v%vertex(icell,2)) 
      p2%z=v%ys(v%vertex(icell,2)) 
      p3%r=v%xs(v%vertex(icell,3)) 
      p3%z=v%ys(v%vertex(icell,3)) 
      denom=1.0_r8/((p2%r-p1%r)*(p3%z-p1%z)-(p3%r-p1%r)*(p2%z-p1%z)) 
!-----------------------------------------------------------------------
!     Evaluate derivatives and functions.                               
!-----------------------------------------------------------------------
      DO iqty=1,u%nqty 
        fx=denom*(                                                      &
     &          (u%fs(iqty,v%vertex(icell,2),0)-                        &
     &           u%fs(iqty,v%vertex(icell,1),0))*(p3%z-p1%z)-           &
     &          (u%fs(iqty,v%vertex(icell,3),0)-                        &
     &           u%fs(iqty,v%vertex(icell,1),0))*(p2%z-p1%z))           
        fy=denom*(                                                      &
     &          (u%fs(iqty,v%vertex(icell,3),0)-                        &
     &           u%fs(iqty,v%vertex(icell,1),0))*(p2%r-p1%r)-           &
     &          (u%fs(iqty,v%vertex(icell,2),0)-                        &
     &           u%fs(iqty,v%vertex(icell,1),0))*(p3%r-p1%r))           
        u%f(iqty)=                                                      &
     &          u%fs(iqty,v%vertex(icell,1),0)+                         &
     &          fy*(p%z-p1%z)+                                          &
     &          fx*(p%r-p1%r)                                           
        IF(mode /= 0 )THEN 
          u%fx(iqty) = fx 
          u%fy(iqty) = fy 
        ENDIF 
      ENDDO 
                                                                        
      RETURN 
      END SUBROUTINE tri_linear_2D_eval 
!-----------------------------------------------------------------------
!     subprogram 12. tri_linear_3D_assign_rsc.                          
!     assign a real scalar value to a 3D tri_linear structure.          
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_assign_rsc(tl,rscalar) 
                                                                        
      TYPE(tri_linear_type), INTENT(INOUT) :: tl 
      REAL(r8), INTENT(IN) :: rscalar 
                                                                        
!-----------------------------------------------------------------------
!     grid vertex nodes completely represent the scalar.                
!-----------------------------------------------------------------------
      tl%fs=rscalar 
!-PRE IF (ASSOCIATED(tl%fss)) THEN                                      
!       tl%fss=0                                                        
!       tl%fsi=0                                                        
!     ENDIF                                                             
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_3D_assign_rsc 
!-----------------------------------------------------------------------
!     subprogram 13. tri_linear_3D_assign_csc.                          
!     assign a complex scalar value to a 3D tri_linear structure.       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_assign_csc(tl,cscalar) 
                                                                        
      TYPE(tri_linear_type), INTENT(INOUT) :: tl 
      COMPLEX(r8), INTENT(IN) :: cscalar 
                                                                        
!-----------------------------------------------------------------------
!     grid vertex nodes completely represent the scalar.                
!-----------------------------------------------------------------------
      tl%fs=cscalar 
!-PRE IF (ASSOCIATED(tl%fss)) THEN                                      
!       tl%fss=0                                                        
!       tl%fsi=0                                                        
!     ENDIF                                                             
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_3D_assign_csc 
!-----------------------------------------------------------------------
!     subprogram 14. tri_linear_3D_assign_tl3.                          
!     set one 3D tri_linear structure equal to another.                 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_assign_tl3(tl1,tl2) 
                                                                        
      TYPE(tri_linear_type), INTENT(INOUT) :: tl1 
      TYPE(tri_linear_type), INTENT(IN) :: tl2 
                                                                        
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.                   
!-----------------------------------------------------------------------
      tl1%fs=tl2%fs 
!-PRE IF (ASSOCIATED(tl1%fss)) THEN                                     
!       tl1%fss=tl2%fss                                                 
!       tl1%fsi=tl2%fsi                                                 
!     ENDIF                                                             
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_3D_assign_tl3 
!-----------------------------------------------------------------------
!     subprogram 15. tri_linear_3D_assign_int.                          
!     assign a integer value to a 3D tri_linear structure.              
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_assign_int(tl,int) 
                                                                        
      TYPE(tri_linear_type), INTENT(INOUT) :: tl 
      INTEGER(i4), INTENT(IN) :: int 
                                                                        
!-----------------------------------------------------------------------
!     grid vertex nodes completely represent the scalar.                
!-----------------------------------------------------------------------
      tl%fs=int 
!-PRE IF (ASSOCIATED(tl%fss)) THEN                                      
!       tl%fss=0                                                        
!       tl%fsi=0                                                        
!     ENDIF                                                             
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_3D_assign_int 
!-----------------------------------------------------------------------
!     subprogram 12. tri_linear_2D_assign_rsc.                          
!     assign a real scalar value to a 2D tri_linear structure.          
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_assign_rsc(tl,rscalar) 
                                                                        
      TYPE(tri_linear_2D_type), INTENT(INOUT) :: tl 
      REAL(r8), INTENT(IN) :: rscalar 
                                                                        
!-----------------------------------------------------------------------
!     grid vertex nodes completely represent the scalar.                
!-----------------------------------------------------------------------
      tl%fs=rscalar 
!-PRE IF (ASSOCIATED(tl%fss)) THEN                                      
!       tl%fss=0                                                        
!       tl%fsi=0                                                        
!     ENDIF                                                             
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_2D_assign_rsc 
!-----------------------------------------------------------------------
!     subprogram 12. tri_linear_2D_assign_csc.                          
!     assign a complex scalar value to a 2D tri_linear structure.       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_assign_csc(tl,cscalar) 
                                                                        
      TYPE(tri_linear_2D_type), INTENT(INOUT) :: tl 
      COMPLEX(r8), INTENT(IN) :: cscalar 
                                                                        
!-----------------------------------------------------------------------
!     grid vertex nodes completely represent the scalar.                
!-----------------------------------------------------------------------
      tl%fs=cscalar 
!-PRE IF (ASSOCIATED(tl%fss)) THEN                                      
!       tl%fss=0                                                        
!       tl%fsi=0                                                        
!     ENDIF                                                             
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_2D_assign_csc 
!-----------------------------------------------------------------------
!     subprogram 14. tri_linear_2D_assign_tl2.                          
!     set one 2D tri_linear structure equal to another.                 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_assign_tl2(tl1,tl2) 
                                                                        
      TYPE(tri_linear_2D_type), INTENT(INOUT) :: tl1 
      TYPE(tri_linear_2D_type), INTENT(IN) :: tl2 
                                                                        
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.                   
!-----------------------------------------------------------------------
      tl1%fs=tl2%fs 
!-PRE IF (ASSOCIATED(tl1%fss)) THEN                                     
!       tl1%fss=tl2%fss                                                 
!       tl1%fsi=tl2%fsi                                                 
!     ENDIF                                                             
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_2D_assign_tl2 
!-----------------------------------------------------------------------
!     subprogram 15. tri_linear_2D_assign_int.                          
!     assign a integer value to a 2D tri_linear structure.              
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_assign_int(tl,int) 
                                                                        
      TYPE(tri_linear_2D_type), INTENT(INOUT) :: tl 
      INTEGER(i4), INTENT(IN) :: int 
                                                                        
!-----------------------------------------------------------------------
!     grid vertex nodes completely represent the scalar.                
!-----------------------------------------------------------------------
      tl%fs=int 
!-PRE IF (ASSOCIATED(tl%fss)) THEN                                      
!       tl%fss=0                                                        
!       tl%fsi=0                                                        
!     ENDIF                                                             
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tri_linear_2D_assign_int 
!-----------------------------------------------------------------------
!     subprogram 16. dump_write_tri_linear.                             
!     real and imaginary parts are witten separately for possible       
!     data conversion.                                                  
!-----------------------------------------------------------------------
      SUBROUTINE dump_write_tri_linear(tb_l) 
      USE io 
                                                                        
      TYPE(tri_linear_type), INTENT(IN) :: tb_l 
                                                                        
      WRITE(dump_unit) REAL(tb_l%nqty,r8) 
      WRITE(dump_unit) REAL(tb_l%nfour,r8) 
      WRITE(dump_unit) REAL(tb_l%fs,r8) 
      WRITE(dump_unit) REAL(-(0,1)*tb_l%fs,r8) 
                                                                        
      RETURN 
      END SUBROUTINE dump_write_tri_linear 
!-----------------------------------------------------------------------
!     subprogram 17. dump_write_tri_linear_2D.                          
!-----------------------------------------------------------------------
      SUBROUTINE dump_write_tri_linear_2D(tb_l) 
      USE io 
                                                                        
      TYPE(tri_linear_2D_type), INTENT(IN) :: tb_l 
                                                                        
      WRITE(dump_unit) REAL(tb_l%nqty,r8) 
      WRITE(dump_unit) tb_l%fs 
                                                                        
      RETURN 
      END SUBROUTINE dump_write_tri_linear_2D 
!-----------------------------------------------------------------------
!     subprogram 18. dump_read_tri_linear.                              
!     real and imaginary parts are read separately for possible         
!     data conversion.                                                  
!-----------------------------------------------------------------------
      SUBROUTINE dump_read_tri_linear(tb_l,lv,name,title) 
      USE io 
                                                                        
      INTEGER(i4), INTENT(IN) :: lv 
      CHARACTER(*), INTENT(IN) :: name 
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title 
      TYPE(tri_linear_type), INTENT(OUT) :: tb_l 
                                                                        
      REAL(r8) :: iread 
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: rread 
                                                                        
      READ(rstrt_unit) iread 
      tb_l%nqty=NINT(iread) 
      READ(rstrt_unit) iread 
      tb_l%nfour=NINT(iread) 
      CALL tri_linear_alloc(tb_l,lv,tb_l%nqty,tb_l%nfour) 
      tb_l%name=name 
      IF (SIZE(title)<SIZE(tb_l%title)) THEN 
        tb_l%title=title(1) 
      ELSE 
        tb_l%title=title 
      ENDIF 
                                                                        
      ALLOCATE(rread(tb_l%nqty,0:lv,0:0,tb_l%nfour)) 
      READ(rstrt_unit) rread 
      tb_l%fs=rread 
      READ(rstrt_unit) rread 
      tb_l%fs=tb_l%fs+(0,1)*rread 
      DEALLOCATE(rread) 
                                                                        
      RETURN 
      END SUBROUTINE dump_read_tri_linear 
!-----------------------------------------------------------------------
!     subprogram 19. dump_read_tri_linear_2D.                           
!-----------------------------------------------------------------------
      SUBROUTINE dump_read_tri_linear_2D(tb_l,lv,name,title) 
      USE io 
                                                                        
      INTEGER(i4), INTENT(IN) :: lv 
      CHARACTER(*), INTENT(IN) :: name 
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title 
      TYPE(tri_linear_2D_type), INTENT(OUT) :: tb_l 
                                                                        
      REAL(r8) :: iread 
                                                                        
      READ(rstrt_unit) iread 
      tb_l%nqty=NINT(iread) 
      CALL tri_linear_alloc(tb_l,lv,tb_l%nqty) 
      tb_l%name=name 
      IF (SIZE(title)<SIZE(tb_l%title)) THEN 
        tb_l%title=title(1) 
      ELSE 
        tb_l%title=title 
      ENDIF 
      READ(rstrt_unit) tb_l%fs 
                                                                        
      RETURN 
      END SUBROUTINE dump_read_tri_linear_2D 
!-----------------------------------------------------------------------
!     close module.                                                     
!-----------------------------------------------------------------------
      END MODULE tri_linear 
