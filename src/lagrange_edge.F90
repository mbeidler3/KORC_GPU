!-----------------------------------------------------------------------
!     $Id: lagrange_edge.F90 6998 2020-01-07 21:32:33Z ehowell $
!     routines for evaluating Lagrange finite elements on edge of the bl
!     structured quadrilaterals.  this package is based on lagr_quad.f  
!                                                                       
!     throughout this package, basis functions are catagorized into node
!     (grid vertex) and side for each element.                          
!-----------------------------------------------------------------------
#include "config.f"
!-----------------------------------------------------------------------
!     code organization.                                                
!-----------------------------------------------------------------------
!      0.  lagr_edge_mod.                                               
!      1.  lagr_edge_bases.                                             
!      2.  lagr_edge_2D_alloc.                                                                            
!      4.  lagr_edge_2D_eval.                                           
!      5.  lagr_edge_2D_all_eval.                                       
!      6.  lagr_edge_2D_assign_rsc.                                     
!      7.  lagr_edge_2D_assign_csc.                                     
!      8.  lagr_edge_2D_assign_laq.                                     
!      9.  lagr_edge_2D_assign_int.                                     
!      10. lagr_edge_2D_basis_assign_arr                                
!      11. lagr_edge_2D_basis_add_arr                                   
!      12. lagr_edge_2D_basis_assign_loc                                
!      13. lagr_edge_2D_basis_add_loc                                   
!      14. lagr_edge_1D_alloc.                                                                               
!      16. lagr_edge_1D_eval.                                           
!      17. lagr_edge_1D_all_eval.                                       
!      18. lagr_edge_1D_assign_rsc.                                     
!      19. lagr_edge_1D_assign_csc.                                     
!      20. lagr_edge_1D_assign_laq.                                     
!      21. lagr_edge_1D_assign_int.                                     
!      22. lagr_edge_1D_basis_assign_arr                                
!      23. lagr_edge_1D_basis_add_arr                                   
!      24. lagr_edge_1D_basis_assign_loc                                
!      25. lagr_edge_1D_basis_add_loc   
!      26. h5_dump_lagr_edge
!      27. h5_dump_lagr_edge_2D
!      28. h5_read_lagr_edge
!      29. h5_read_lagr_edge_2d                                
!-----------------------------------------------------------------------
!     subprogram 0. lagr_edge_mod definition.                           
!     defines the lagrangian quadrilateral type.                        
!-----------------------------------------------------------------------
      MODULE lagr_edge_mod 
      USE local 
      IMPLICIT NONE 
                                                                        
      TYPE :: lagr_edge_type 
        INTEGER(i4) :: mx,nqty,nfour,n_side,n_int 
        !INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ix0 
        REAL(r8), DIMENSION(:), ALLOCATABLE :: dx 
        COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: fs 
        COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: fsh 
        COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: f,fx,fy 
        CHARACTER(6), DIMENSION(:), ALLOCATABLE :: title 
        CHARACTER(6) :: name 
      END TYPE lagr_edge_type 
                                                                        
      TYPE :: lagr_edge_1D_type 
        INTEGER(i4) :: mx,nqty,n_side,n_int 
        !INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ix0 
        REAL(r8), DIMENSION(:), ALLOCATABLE :: dx 
        REAL(r8), DIMENSION(:,:), ALLOCATABLE :: fs 
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: fsh 
        REAL(r8), DIMENSION(:), ALLOCATABLE :: f,fx,fy 
        CHARACTER(6), DIMENSION(:), ALLOCATABLE :: title 
        CHARACTER(6) :: name 
      END TYPE lagr_edge_1D_type 
                                                                        
      INTEGER(i4), PARAMETER, PRIVATE :: npoly_max=20 
      REAL(r8), DIMENSION(0:npoly_max), PRIVATE :: alx,aly,dalx,daly 
!-----------------------------------------------------------------------
!     subprogram name interfaces                                        
!-----------------------------------------------------------------------
      INTERFACE lagr_edge_alloc 
        MODULE PROCEDURE lagr_edge_1D_alloc,lagr_edge_2D_alloc 
      END INTERFACE 
                                                                                                                                           
      INTERFACE lagr_edge_all_eval 
        MODULE PROCEDURE lagr_edge_1D_all_eval,lagr_edge_2D_all_eval 
      END INTERFACE 
                                                                        
      INTERFACE lagr_edge_eval 
        MODULE PROCEDURE lagr_edge_1D_eval,lagr_edge_2D_eval 
      END INTERFACE 
                                                                        
      INTERFACE lagr_edge_basis_assign_arr 
        MODULE PROCEDURE                                                &
     &    lagr_edge_1D_basis_assign_arr,lagr_edge_2D_basis_assign_arr   
      END INTERFACE 
                                                                        
      INTERFACE lagr_edge_basis_add_arr 
        MODULE PROCEDURE                                                &
     &    lagr_edge_1D_basis_add_arr,lagr_edge_2D_basis_add_arr         
      END INTERFACE 
                                                                        
      INTERFACE lagr_edge_basis_assign_loc 
        MODULE PROCEDURE                                                &
     &    lagr_edge_1D_basis_assign_loc,lagr_edge_2D_basis_assign_loc   
      END INTERFACE 
                                                                        
      INTERFACE lagr_edge_basis_add_loc 
        MODULE PROCEDURE                                                &
     &    lagr_edge_1D_basis_add_loc,lagr_edge_2D_basis_add_loc         
      END INTERFACE 
                                                                        
!-----------------------------------------------------------------------
!     overloaded assignment.                                            
!-----------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=) 
        MODULE PROCEDURE                                                &
     &    lagr_edge_2D_assign_csc,lagr_edge_2D_assign_rsc,              &
     &    lagr_edge_2D_assign_laq,lagr_edge_2D_assign_int,              &
     &    lagr_edge_1D_assign_csc,lagr_edge_1D_assign_rsc,              &
     &    lagr_edge_1D_assign_laq,lagr_edge_1D_assign_int               
      END INTERFACE 
                                                                        
      CONTAINS 
!-----------------------------------------------------------------------
!     subprogram 1. lagr_edge_bases.                                    
!     computes basis functions and their derivatives.                   
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_bases(x,alpha,alphax,dmode) 
                                                                        
      REAL(r8), INTENT(IN) :: x 
      REAL(r8), DIMENSION(:), INTENT(INOUT) :: alpha,alphax 
      INTEGER(i4), INTENT(IN) :: dmode 
                                                                        
      INTEGER(i4) :: pd,i,k 
!-----------------------------------------------------------------------
!     interface block for lagr_1D external routine.                     
!-----------------------------------------------------------------------
      INTERFACE 
        SUBROUTINE lagr_1D(pd,x,al,dal,dmode,d2al) 
        USE local 
        INTEGER(i4), INTENT(IN) :: pd,dmode 
        REAL(r8), INTENT(IN) :: x 
        REAL(r8), DIMENSION(0:), INTENT(OUT) :: al,dal 
        REAL(r8), DIMENSION(0:), INTENT(OUT), OPTIONAL :: d2al
        END SUBROUTINE lagr_1D 
      END INTERFACE 
!-----------------------------------------------------------------------
!     first determine the polynomial degree and set 1D basis functions. 
!-----------------------------------------------------------------------
      pd=NINT(REAL(SIZE(alpha)))-1 
      CALL lagr_1D(pd,x,alx,dalx,dmode) 
!-----------------------------------------------------------------------
!     if derivatives are not needed, limit the computations.            
!-----------------------------------------------------------------------
      IF (dmode==0) THEN 
!-----------------------------------------------------------------------
!       compute basis functions with respect to the logical coordinates 
!-----------------------------------------------------------------------
        alpha(1)=alx(0) 
        alpha(2)=alx(pd) 
!-----------------------------------------------------------------------
!       horizontal side centered bases if needed.                       
!-----------------------------------------------------------------------
        k=3 
        DO i=1,pd-1 
          alpha(k)= alx(i) 
          k=k+1 
        ENDDO 
      ELSE 
!-----------------------------------------------------------------------
!       compute basis functions and derivatives with respect to the     
!       logical coordinates of the unit square.  node (grid vertex-     
!       centered) bases first; start in lower left corner, work left to 
!       right, bottom to top.                                           
!-----------------------------------------------------------------------
        alpha(1)=alx(0) 
        alpha(2)=alx(pd) 
        alphax(1)=dalx(0) 
        alphax(2)=dalx(pd) 
!-----------------------------------------------------------------------
!       horizontal side centered bases if needed.  the ordering is      
!       bottom to top, then left to right in pairs.                     
!-----------------------------------------------------------------------
        k=3 
        DO i=1,pd-1 
          alpha (k)= alx(i) 
          alphax(k)=dalx(i) 
          k=k+1 
        ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_bases 
!-----------------------------------------------------------------------
!     subprogram 2. lagr_edge_2D_alloc.                                 
!     allocates space for lagr_edge_2D_type.                            
!     NOTE: the lagr_edge_type has only an fsh.  Because for periodic   
!      geometries (cylinder, flux, rect w/ periodicity=y) the "edge" is 
!      in the vertical direction, in a sense, we are matching the fsh   
!      arrays with fsv arrays. Perhaps confusing for the most common    
!      case.                                                            
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_2D_alloc(laq,mx,nqty,nfour,poly_degree,      &
     &                              name,title)                         
                                                                        
      INTEGER(i4), INTENT(IN) :: mx,nqty,nfour,poly_degree 
      CHARACTER(*), INTENT(IN), OPTIONAL :: name 
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title 
      TYPE(lagr_edge_type), INTENT(OUT) :: laq 
                                                                        
      INTEGER(i4) :: ix,ib,ip 
      REAL(r8), DIMENSION(0:poly_degree) :: x_node 
!-----------------------------------------------------------------------
!     store grid, vector, and fourier series dimensions, and set the    
!     number of side and interior basis functions.                      
!-----------------------------------------------------------------------
      laq%mx=mx 
      laq%nqty=nqty 
      laq%nfour=nfour 
      laq%n_side=poly_degree-1 
!-----------------------------------------------------------------------
!     allocate space.                                                   
!-----------------------------------------------------------------------
      ALLOCATE(laq%fs(nqty,0:mx,nfour)); laq%fs=0. 
      IF (poly_degree>1) THEN 
        ALLOCATE(laq%fsh(nqty,laq%n_side,1:mx,nfour)); laq%fsh=0. 
      ENDIF 
      ALLOCATE(laq%title(nqty)) 
      ALLOCATE(laq%f(nqty,nfour)); laq%f=0. 
      ALLOCATE(laq%fx(nqty,nfour)); laq%fx=0. 
!-----------------------------------------------------------------------
!     character descriptors, if present in input.                       
!-----------------------------------------------------------------------
      IF (PRESENT(name)) laq%name=name 
      IF (PRESENT(title)) THEN 
        IF (SIZE(title)==nqty) THEN 
          laq%title=title 
        ELSE 
          laq%title=title(1) 
        ENDIF 
      ENDIF 
!-----------------------------------------------------------------------
!     compute basis starting indices and centerings in logical space.   
!     note that the centerings are for standard lagrange polynomial     
!     bases.                                                            
!                                                                       
!     dx gives the relative positions within a cell for the             
!     different bases.  note that these quantities may be used to access
!     all data with loops like:                                         
!-----------------------------------------------------------------------
      !ALLOCATE(laq%ix0(poly_degree)) 
      ALLOCATE(laq%dx(poly_degree)) 
      CALL poly_nodes(poly_degree,x_node) 
      DO ib=1,poly_degree 
        IF (ib==1) THEN 
          !laq%ix0(ib)=0 
          laq%dx(ib)=0 
        ELSE 
          !laq%ix0(ib)=1 
          laq%dx(ib)=x_node(ib-1_i4) 
        ENDIF 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_2D_alloc 
!-----------------------------------------------------------------------
!     subprogram 4. lagr_edge_2D_eval.                                  
!     evaluates complex lagr_edge quantities at a single point within a 
!     grid block.                                                       
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_2D_eval(laq,x,dmode) 
                                                                        
      TYPE(lagr_edge_type), INTENT(INOUT) :: laq 
      REAL(r8), INTENT(IN) :: x 
      INTEGER(i4), INTENT(IN) :: dmode 
                                                                        
      INTEGER(i4) :: ix,pd,i,k,im 
!-----------------------------------------------------------------------
!     interface block for lagr_1D external routine.                     
!-----------------------------------------------------------------------
      INTERFACE 
        SUBROUTINE lagr_1D(pd,x,al,dal,dmode,d2al) 
        USE local 
        INTEGER(i4), INTENT(IN) :: pd,dmode 
        REAL(r8), INTENT(IN) :: x 
        REAL(r8), DIMENSION(0:), INTENT(OUT) :: al,dal 
        REAL(r8), DIMENSION(0:), INTENT(OUT), OPTIONAL :: d2al
        END SUBROUTINE lagr_1D 
      END INTERFACE 
!-----------------------------------------------------------------------
!     find the interval, and compute 1D basis coefficients.             
!-----------------------------------------------------------------------
      ix=MAX(MIN(INT(x),laq%mx-1),0_i4) 
      pd=laq%n_side+1 
      CALL lagr_1D(pd,x-ix,alx,dalx,dmode) 
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.  bilinear   
!     derivatives simplify.                                             
!-----------------------------------------------------------------------
      IF (pd==1) THEN 
        laq%f=laq%fs(:,ix  ,:)*alx(0)+laq%fs(:,ix+1,:)*alx(1) 
        IF (dmode==0) RETURN 
        laq%fx=laq%fs(:,ix+1,:)-laq%fs(:,ix,:) 
!-----------------------------------------------------------------------
!     other polynomials:                                                
!-----------------------------------------------------------------------
      ELSE 
        laq%f=laq%fs(:,ix,:)*alx(0)+laq%fs(:,ix+1,:)*alx(pd) 
        DO i=1,pd-1 
          laq%f=laq%f+laq%fsh(:,i,ix+1,:)*alx(i) 
        ENDDO 
        IF (dmode==0) RETURN 
        laq%fx=laq%fs(:,ix  ,:)*dalx(0)+laq%fs(:,ix+1,:)*dalx(pd) 
        DO i=1,pd-1 
          laq%fx=laq%fx+laq%fsh(:,i,ix+1,:)*dalx(i) 
        ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_2D_eval 
!-----------------------------------------------------------------------
!     subprogram 5. lagr_edge_2D_all_eval.                              
!     evaluates complex lagr_edge quantities in all elements in a grid  
!     block for equal spacing.                                          
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_2D_all_eval(laq,x,f,fx,dmode) 
                                                                        
      TYPE(lagr_edge_type), INTENT(IN) :: laq 
      REAL(r8), INTENT(IN) :: x 
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: f,fx 
      INTEGER(i4), INTENT(IN) :: dmode 
                                                                        
      INTEGER(i4) :: mx 
      INTEGER(i4) :: pd,i,im,ix,j 
      REAL(r8), DIMENSION(laq%n_side+2) :: alpha,dalpdx 
!-----------------------------------------------------------------------
!     compute index limits.                                             
!-----------------------------------------------------------------------
      mx=laq%mx 
      pd=laq%n_side+1 
      CALL lagr_edge_bases(x,alpha,dalpdx,dmode) 
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.  bilinear   
!     derivatives simplify.  [do-loops are for optimization, not        
!     aesthetics.]                                                      
!-----------------------------------------------------------------------
      SELECT CASE(pd) 
      CASE(1) 
        IF (dmode==0) THEN 
          DO im=1,laq%nfour 
            DO ix=1,mx 
              f(:,ix,im)=laq%fs(:,ix-1,im)*alpha(1)+                    &
     &                   laq%fs(:,ix  ,im)*alpha(2)                     
            ENDDO 
          ENDDO 
        ELSE 
          DO im=1,laq%nfour 
            DO ix=1,mx 
              f(:,ix,im)=laq%fs(:,ix-1,im)*alpha(1)+                    &
     &                   laq%fs(:,ix  ,im)*alpha(2)                     
              fx(:,ix,im)=laq%fs(:,ix  ,im)-laq%fs(:,ix-1,im) 
              ENDDO 
          ENDDO 
        ENDIF 
!-----------------------------------------------------------------------
!     biquadratics--no sums.                                            
!-----------------------------------------------------------------------
      CASE(2) 
        IF (dmode==0) THEN 
          DO im=1,laq%nfour 
            DO ix=1,mx 
              f(:,ix,im)=                                               &
     &               laq%fs(:,ix-1,im)*alpha(1)+                        &
     &               laq%fs(:,ix  ,im)*alpha(2)+                        &
     &               laq%fsh(:,1,ix,im)*alpha(3)                        
            ENDDO 
          ENDDO 
        ELSE 
          DO im=1,laq%nfour 
            DO ix=1,mx 
                f(:,ix,im)=                                             &
     &               laq%fs(:,ix-1,im)*alpha(1)+                        &
     &               laq%fs(:,ix  ,im)*alpha(2)+                        &
     &               laq%fsh(:,1,ix,im)*alpha(3)                        
                fx(:,ix,im)=                                            &
     &                laq%fs(:,ix-1,im)*dalpdx(1)+                      &
     &                laq%fs(:,ix  ,im)*dalpdx(2)+                      &
     &                laq%fsh(:,1,ix,im)*dalpdx(3)                      
            ENDDO 
          ENDDO 
        ENDIF 
!-----------------------------------------------------------------------
!     other polynomials:                                                
!-----------------------------------------------------------------------
      CASE DEFAULT 
        j=pd+1 
        IF (dmode==0) THEN 
          DO im=1,laq%nfour 
            DO ix=1,mx 
              DO i=1,laq%nqty 
                f(i,ix,im)=                                             &
     &                 laq%fs(i,ix-1,im)*alpha(1)+                      &
     &                 laq%fs(i,ix  ,im)*alpha(2)+                      &
     &             SUM(laq%fsh(i,:,ix,im)*alpha(3:j))                   
              ENDDO 
            ENDDO 
          ENDDO 
        ELSE 
          DO im=1,laq%nfour 
            DO ix=1,mx 
              DO i=1,laq%nqty 
                f(i,ix,im)=                                             &
     &                 laq%fs(i,ix-1,im)*alpha(1)+                      &
     &                 laq%fs(i,ix  ,im)*alpha(2)+                      &
     &             SUM(laq%fsh(i,:,ix,im)*alpha(3:j))                   
                  fx(i,ix,im)=                                          &
     &                  laq%fs(i,ix-1,im)*dalpdx(1)+                    &
     &                  laq%fs(i,ix  ,im)*dalpdx(2)+                    &
     &              SUM(laq%fsh(i,:,ix,im)*dalpdx(3:j))                 
              ENDDO 
            ENDDO 
          ENDDO 
        ENDIF 
      END SELECT 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_2D_all_eval 
!-----------------------------------------------------------------------
!     subprogram 6. lagr_edge_2D_assign_rsc.                            
!     assign a real scalar value to a complex lagrange edge structure.  
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_2D_assign_rsc(laq,rscalar) 
                                                                        
      TYPE(lagr_edge_type), INTENT(INOUT) :: laq 
      REAL(r8), INTENT(IN) :: rscalar 
                                                                        
      laq%fs=rscalar 
      IF (ALLOCATED(laq%fsh)) laq%fsh=rscalar 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_2D_assign_rsc 
!-----------------------------------------------------------------------
!     subprogram 7. lagr_edge_2D_assign_csc.                            
!     assign a complex scalar value to a complex lagrange edge          
!     structure.                                                        
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_2D_assign_csc(laq,cscalar) 
                                                                        
      TYPE(lagr_edge_type), INTENT(INOUT) :: laq 
      COMPLEX(r8), INTENT(IN) :: cscalar 
                                                                        
      laq%fs=cscalar 
      IF (ALLOCATED(laq%fsh)) laq%fsh=cscalar 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_2D_assign_csc 
!-----------------------------------------------------------------------
!     subprogram 8. lagr_edge_2D_assign_laq.                            
!     set one complex lagrange edge structure equal to another.         
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_2D_assign_laq(laq1,laq2) 
                                                                        
      TYPE(lagr_edge_type), INTENT(INOUT) :: laq1 
      TYPE(lagr_edge_type), INTENT(IN) :: laq2 
                                                                        
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.                   
!-----------------------------------------------------------------------
      laq1%fs=laq2%fs 
      IF (ALLOCATED(laq1%fsh)) laq1%fsh=laq2%fsh 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_2D_assign_laq 
!-----------------------------------------------------------------------
!     subprogram 9. lagr_edge_2D_assign_int.                            
!     assign a integer value to a complex lagrange edge structure.      
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_2D_assign_int(laq,int) 
                                                                        
      TYPE(lagr_edge_type), INTENT(INOUT) :: laq 
      INTEGER(i4), INTENT(IN) :: int 
                                                                        
      laq%fs=int 
      IF (ALLOCATED(laq%fsh)) laq%fsh=int 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_2D_assign_int 
!-----------------------------------------------------------------------
!     subprogram 10. lagr_edge_2D_basis_assign_arr                      
!     assign data into coefficient arrays for one basis function.       
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_2D_basis_assign_arr(laq,data,ibasis) 
                                                                        
      TYPE(lagr_edge_type), INTENT(INOUT) :: laq 
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: data 
      INTEGER(i4), INTENT(IN) :: ibasis 
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to     
!     the value of ibasis.                                              
!     grid vertex-centered data is first.                               
!-----------------------------------------------------------------------
      IF (ibasis==1) THEN 
        laq%fs=data 
        RETURN 
!-----------------------------------------------------------------------
!     horizontal sides.                                                 
!-----------------------------------------------------------------------
      ELSE 
        laq%fsh(:,ibasis-1,:,:)=data 
        RETURN 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_2D_basis_assign_arr 
!-----------------------------------------------------------------------
!     subprogram 11. lagr_edge_2D_basis_add_arr                         
!     add data into coefficient arrays for one basis function.          
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_2D_basis_add_arr(laq,data,ibasis) 
                                                                        
      TYPE(lagr_edge_type), INTENT(INOUT) :: laq 
      COMPLEX(r8), DIMENSION(:,:,:), INTENT(IN) :: data 
      INTEGER(i4), INTENT(IN) :: ibasis 
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to     
!     the value of ibasis.                                              
!     grid vertex-centered data is first.                               
!-----------------------------------------------------------------------
      IF (ibasis==1) THEN 
        laq%fs=laq%fs+data 
        RETURN 
!-----------------------------------------------------------------------
!     horizontal sides.                                                 
!-----------------------------------------------------------------------
      ELSE 
        laq%fsh(:,ibasis-1,:,:)=laq%fsh(:,ibasis-1,:,:)+data 
        RETURN 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_2D_basis_add_arr 
!-----------------------------------------------------------------------
!     subprogram 12. lagr_edge_2D_basis_assign_loc                      
!     assign data into coefficient arrays for one basis function.       
!                                                                       
!     this is a local version of lagr_edge_basis_assign_arr, where the  
!     the data is located at a given poloidal and fourier indices       
!     triplet, only.                                                    
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_2D_basis_assign_loc(laq,data,ibasis,ix,im) 
                                                                        
      TYPE(lagr_edge_type), INTENT(INOUT) :: laq 
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: data 
      INTEGER(i4), INTENT(IN) :: ibasis,ix,im 
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to     
!     the value of ibasis.                                              
!     grid vertex-centered data is first.                               
!-----------------------------------------------------------------------
      IF (ibasis==1) THEN 
        laq%fs(:,ix,im)=data 
        RETURN 
!-----------------------------------------------------------------------
!     horizontal sides.                                                 
!-----------------------------------------------------------------------
      ELSE 
        laq%fsh(:,ibasis-1,ix,im)=data 
        RETURN 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_2D_basis_assign_loc 
!-----------------------------------------------------------------------
!     subprogram 13. lagr_edge_2D_basis_add_loc                         
!     add data into coefficient arrays for one basis function.          
!                                                                       
!     this is a local version of lagr_edge_basis_add_arr, where the     
!     the data is located at a given poloidal and fourier indices       
!     triplet, only.                                                    
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_2D_basis_add_loc(laq,data,ibasis,ix,im) 
                                                                        
      TYPE(lagr_edge_type), INTENT(INOUT) :: laq 
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: data 
      INTEGER(i4), INTENT(IN) :: ibasis,ix,im 
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to     
!     the value of ibasis.                                              
!     grid vertex-centered data is first.                               
!-----------------------------------------------------------------------
      IF (ibasis==1) THEN 
        laq%fs(:,ix,im)=laq%fs(:,ix,im)+data 
        RETURN 
!-----------------------------------------------------------------------
!     horizontal sides.                                                 
!-----------------------------------------------------------------------
      ELSE 
        laq%fsh(:,ibasis-1,ix,im)=laq%fsh(:,ibasis-1,ix,im)+data 
        RETURN 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_2D_basis_add_loc 
!-----------------------------------------------------------------------
!     subprogram 14. lagr_edge_1D_alloc.                                
!     allocates space for lagr_edge_1D_type.                            
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_1D_alloc(laq,mx,nqty,poly_degree,name,title) 
                                                                        
      INTEGER(i4), INTENT(IN) :: mx,nqty,poly_degree 
      CHARACTER(*), INTENT(IN), OPTIONAL :: name 
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title 
      TYPE(lagr_edge_1D_type), INTENT(OUT) :: laq 
                                                                        
      INTEGER(i4) :: ix,ib,ip,pd1 
      REAL(r8), DIMENSION(0:poly_degree) :: x_node 
!-----------------------------------------------------------------------
!     store grid and vector dimensions, and set the                     
!     number of side and interior basis functions.                      
!-----------------------------------------------------------------------
      laq%mx=mx 
      laq%nqty=nqty 
      laq%n_side=poly_degree-1 
!-----------------------------------------------------------------------
!     allocate space.                                                   
!-----------------------------------------------------------------------
      ALLOCATE(laq%fs(nqty,0:mx)); laq%fs=0. 
      IF (poly_degree>1) THEN 
        ALLOCATE(laq%fsh(nqty,laq%n_side,1:mx)); laq%fsh=0. 
      ENDIF 
      ALLOCATE(laq%title(nqty)) 
      ALLOCATE(laq%f(nqty)); laq%f=0. 
      ALLOCATE(laq%fx(nqty)); laq%fx=0. 
      ALLOCATE(laq%fy(nqty)); laq%fy=0. 
!-----------------------------------------------------------------------
!     character descriptors, if present in input.                       
!-----------------------------------------------------------------------
      IF (PRESENT(name)) laq%name=name 
      IF (PRESENT(title)) THEN 
        IF (SIZE(title)==nqty) THEN 
          laq%title=title 
        ELSE 
          laq%title=title(1) 
        ENDIF 
      ENDIF 
!-----------------------------------------------------------------------
!     compute basis starting indices and centerings in logical space.   
!     note that the centerings are for standard lagrange polynomial     
!     bases.                                                            
!-----------------------------------------------------------------------
      pd1=poly_degree-1 
      !ALLOCATE(laq%ix0(poly_degree)) 
      ALLOCATE(laq%dx(poly_degree)) 
      CALL poly_nodes(poly_degree,x_node) 
      DO ib=1,poly_degree 
        IF (ib==1) THEN 
          !laq%ix0(ib)=0 
          laq%dx(ib)=0 
        ELSE IF (ib<=poly_degree) THEN 
          !laq%ix0(ib)=1 
          laq%dx(ib)=x_node(ib-1_i4) 
        ENDIF 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_1D_alloc 
!-----------------------------------------------------------------------
!     subprogram 16. lagr_edge_1D_eval.                                 
!     evaluates real lagr_edge quantities at a single point within a    
!     grid block.                                                       
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_1D_eval(laq,x,dmode) 
                                                                        
      TYPE(lagr_edge_1D_type), INTENT(INOUT) :: laq 
      REAL(r8), INTENT(IN) :: x 
      INTEGER(i4), INTENT(IN) :: dmode 
                                                                        
      INTEGER(i4) :: ix,pd,i,im 
!-----------------------------------------------------------------------
!     interface block for lagr_1D external routine.                     
!-----------------------------------------------------------------------
      INTERFACE 
        SUBROUTINE lagr_1D(pd,x,al,dal,dmode,d2al) 
        USE local 
        INTEGER(i4), INTENT(IN) :: pd,dmode 
        REAL(r8), INTENT(IN) :: x 
        REAL(r8), DIMENSION(0:), INTENT(OUT) :: al,dal 
        REAL(r8), DIMENSION(0:), INTENT(OUT), OPTIONAL :: d2al
        END SUBROUTINE lagr_1D 
      END INTERFACE 
!-----------------------------------------------------------------------
!     find the interval, and compute 1D basis coefficients.             
!-----------------------------------------------------------------------
      ix=MAX(MIN(INT(x),laq%mx-1),0_i4) 
      pd=laq%n_side+1 
      CALL lagr_1D(pd,x-ix,alx,dalx,dmode) 
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.  bilinear   
!     derivatives simplify.                                             
!-----------------------------------------------------------------------
      IF (pd==1) THEN 
        laq%f=laq%fs(:,ix)*alx(0)+laq%fs(:,ix+1)*alx(1) 
        IF (dmode==0) RETURN 
        laq%fx=laq%fs(:,ix+1)-laq%fs(:,ix) 
!-----------------------------------------------------------------------
!     other polynomials:                                                
!-----------------------------------------------------------------------
      ELSE 
        laq%f=laq%fs(:,ix  )*alx(0)+laq%fs(:,ix+1)*alx(pd) 
        DO i=1,pd-1 
          laq%f=laq%f+laq%fsh(:,i,ix+1)*alx(i) 
        ENDDO 
        IF (dmode==0) RETURN 
        laq%fx=laq%fs(:,ix  )*dalx(0)+laq%fs(:,ix+1)*dalx(pd) 
        DO i=1,pd-1 
          laq%fx=laq%fx+laq%fsh(:,i,ix+1)*dalx(i) 
        ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_1D_eval 
!-----------------------------------------------------------------------
!     subprogram 17. lagr_edge_1D_all_eval.                             
!     evaluates real lagr_edge quantities in all elements in a grid     
!     block for equal spacing.                                          
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_1D_all_eval(laq,x,f,fx,dmode) 
                                                                        
      TYPE(lagr_edge_1D_type), INTENT(IN) :: laq 
      REAL(r8), INTENT(IN) :: x 
      REAL(r8), INTENT(OUT), DIMENSION(:,:) :: f,fx 
      INTEGER(i4), INTENT(IN) :: dmode 
                                                                        
      INTEGER(i4) :: mx,pd,i,j,ix 
      REAL(r8), DIMENSION(laq%n_side+2) :: alpha,dalpdx 
!-----------------------------------------------------------------------
!     compute index limits.                                             
!-----------------------------------------------------------------------
      mx=laq%mx 
      pd=laq%n_side+1 
      CALL lagr_edge_bases(x,alpha,dalpdx,dmode) 
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.  bilinear   
!     derivatives simplify.  [do-loops are for optimization, not        
!     aesthetics.]                                                      
!-----------------------------------------------------------------------
      SELECT CASE(pd) 
      CASE(1) 
        IF (dmode==0) THEN 
          DO ix=1,mx 
            f(:,ix)=laq%fs(:,ix-1)*alpha(1)+laq%fs(:,ix)*alpha(2) 
          ENDDO 
        ELSE 
          DO ix=1,mx 
            f(:,ix)=laq%fs(:,ix-1)*alpha(1)+laq%fs(:,ix)*alpha(2) 
            fx(:,ix)=laq%fs(:,ix)-laq%fs(:,ix-1) 
          ENDDO 
        ENDIF 
!-----------------------------------------------------------------------
!     biquadratics--no sums.                                            
!-----------------------------------------------------------------------
      CASE(2) 
        IF (dmode==0) THEN 
          DO ix=1,mx 
            f(:,ix)=laq%fs(:,ix-1) *alpha(1)+laq%fs(:,ix)*alpha(2)+     &
     &              laq%fsh(:,1,ix)*alpha(3)                            
          ENDDO 
        ELSE 
          DO ix=1,mx 
            f(:,ix)=laq%fs(:,ix-1)*alpha(1)+laq%fs(:,ix  )*alpha(2)+    &
     &             laq%fsh(:,1,ix)*alpha(3)                             
            fx(:,ix)=laq%fs(:,ix-1) *dalpdx(1)+laq%fs(:,ix)*dalpdx(2)+  &
     &               laq%fsh(:,1,ix)*dalpdx(3)                          
          ENDDO 
        ENDIF 
!-----------------------------------------------------------------------
!     other polynomials:                                                
!-----------------------------------------------------------------------
      CASE DEFAULT 
        j=pd+1 
        IF (dmode==0) THEN 
          DO ix=1,mx 
            DO i=1,laq%nqty 
              f(i,ix)=laq%fs(i,ix-1)*alpha(1)+laq%fs(i,ix)*alpha(2)+    &
     &                SUM(laq%fsh(i,:,ix)*alpha(3:j))                   
            ENDDO 
          ENDDO 
        ELSE 
          DO ix=1,mx 
            DO i=1,laq%nqty 
              f(i,ix)=laq%fs(i,ix-1)*alpha(1)+laq%fs(i,ix)*alpha(2)+    &
     &                SUM(laq%fsh(i,:,ix)*alpha(3:j))                   
              fx(i,ix)=laq%fs(i,ix-1)*dalpdx(1)+laq%fs(i,ix)*dalpdx(2)+ &
     &             SUM(laq%fsh(i,:,ix)*dalpdx(3:j))                     
              ENDDO 
          ENDDO 
        ENDIF 
      END SELECT 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_1D_all_eval 
!-----------------------------------------------------------------------
!     subprogram 18. lagr_edge_1D_assign_rsc.                           
!     assign a real scalar value to a real lagrange edge structure.     
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_1D_assign_rsc(laq,rscalar) 
                                                                        
      TYPE(lagr_edge_1D_type), INTENT(INOUT) :: laq 
      REAL(r8), INTENT(IN) :: rscalar 
                                                                        
      laq%fs=rscalar 
      IF (ALLOCATED(laq%fsh)) laq%fsh=rscalar 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_1D_assign_rsc 
!-----------------------------------------------------------------------
!     subprogram 19. lagr_edge_1D_assign_csc.                           
!     assign a real scalar value to a real lagrange edge                
!     structure.                                                        
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_1D_assign_csc(laq,cscalar) 
                                                                        
      TYPE(lagr_edge_1D_type), INTENT(INOUT) :: laq 
      COMPLEX(r8), INTENT(IN) :: cscalar 
                                                                        
      laq%fs=cscalar 
      IF (ALLOCATED(laq%fsh)) laq%fsh=cscalar 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_1D_assign_csc 
!-----------------------------------------------------------------------
!     subprogram 20. lagr_edge_1D_assign_laq.                           
!     set one real lagrange edge structure equal to another.            
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_1D_assign_laq(laq1,laq2) 
                                                                        
      TYPE(lagr_edge_1D_type), INTENT(INOUT) :: laq1 
      TYPE(lagr_edge_1D_type), INTENT(IN) :: laq2 
                                                                        
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.                   
!-----------------------------------------------------------------------
      laq1%fs=laq2%fs 
      IF (ALLOCATED(laq1%fsh)) laq1%fsh=laq2%fsh 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_1D_assign_laq 
!-----------------------------------------------------------------------
!     subprogram 21. lagr_edge_1D_assign_int.                           
!     assign a integer value to a real lagrange edge structure.         
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_1D_assign_int(laq,int) 
                                                                        
      TYPE(lagr_edge_1D_type), INTENT(INOUT) :: laq 
      INTEGER(i4), INTENT(IN) :: int 
                                                                        
      laq%fs=int 
      IF (ALLOCATED(laq%fsh)) laq%fsh=int 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_1D_assign_int 
!-----------------------------------------------------------------------
!     subprogram 22. lagr_edge_1D_basis_assign_arr                      
!     assign data into coefficient arrays for one basis function.       
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_1D_basis_assign_arr(laq,data,ibasis) 
                                                                        
      TYPE(lagr_edge_1D_type), INTENT(INOUT) :: laq 
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: data 
      INTEGER(i4), INTENT(IN) :: ibasis 
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to     
!     the value of ibasis.                                              
!     grid vertex-centered data is first.                               
!-----------------------------------------------------------------------
      IF (ibasis==1) THEN 
        laq%fs=data 
        RETURN 
!-----------------------------------------------------------------------
!     horizontal sides.                                                 
!-----------------------------------------------------------------------
      ELSE 
        laq%fsh(:,ibasis-1,:)=data 
        RETURN 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_1D_basis_assign_arr 
!-----------------------------------------------------------------------
!     subprogram 23. lagr_edge_1D_basis_add_arr                         
!     add data into coefficient arrays for one basis function.          
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_1D_basis_add_arr(laq,data,ibasis) 
                                                                        
      TYPE(lagr_edge_1D_type), INTENT(INOUT) :: laq 
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: data 
      INTEGER(i4), INTENT(IN) :: ibasis 
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to     
!     the value of ibasis.                                              
!     grid vertex-centered data is first.                               
!-----------------------------------------------------------------------
      IF (ibasis==1) THEN 
        laq%fs=laq%fs+data 
        RETURN 
!-----------------------------------------------------------------------
!     horizontal sides.                                                 
!-----------------------------------------------------------------------
      ELSE 
        laq%fsh(:,ibasis-1,:)=laq%fsh(:,ibasis-1,:)+data 
        RETURN 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_1D_basis_add_arr 
!-----------------------------------------------------------------------
!     subprogram 24. lagr_edge_1D_basis_assign_loc                      
!     assign data into coefficient arrays for one basis function.       
!                                                                       
!     this is a local version of lagr_edge_basis_assign_arr, where the  
!     the data is located at given poloidal indices, only.              
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_1D_basis_assign_loc(laq,data,ibasis,ix) 
                                                                        
      TYPE(lagr_edge_1D_type), INTENT(INOUT) :: laq 
      REAL(r8), DIMENSION(:), INTENT(IN) :: data 
      INTEGER(i4), INTENT(IN) :: ibasis,ix 
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to     
!     the value of ibasis.                                              
!     grid vertex-centered data is first.                               
!-----------------------------------------------------------------------
      IF (ibasis==1) THEN 
        laq%fs(:,ix)=data 
        RETURN 
!-----------------------------------------------------------------------
!     horizontal sides.                                                 
!-----------------------------------------------------------------------
      ELSE 
        laq%fsh(:,ibasis-1,ix)=data 
        RETURN 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_1D_basis_assign_loc 
!-----------------------------------------------------------------------
!     subprogram 25. lagr_edge_1D_basis_add_loc                         
!     add data into coefficient arrays for one basis function.          
!                                                                       
!     this is a local version of lagr_edge_basis_add_arr, where the     
!     the data is located at given poloidal indices only.               
!-----------------------------------------------------------------------
      SUBROUTINE lagr_edge_1D_basis_add_loc(laq,data,ibasis,ix) 
                                                                        
      TYPE(lagr_edge_1D_type), INTENT(INOUT) :: laq 
      REAL(r8), DIMENSION(:), INTENT(IN) :: data 
      INTEGER(i4), INTENT(IN) :: ibasis,ix 
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to     
!     the value of ibasis.                                              
!     grid vertex-centered data is first.                               
!-----------------------------------------------------------------------
      IF (ibasis==1) THEN 
        laq%fs(:,ix)=laq%fs(:,ix)+data 
        RETURN 
!-----------------------------------------------------------------------
!     horizontal sides.                                                 
!-----------------------------------------------------------------------
      ELSE 
        laq%fsh(:,ibasis-1,ix)=laq%fsh(:,ibasis-1,ix)+data 
        RETURN 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_edge_1D_basis_add_loc
!-----------------------------------------------------------------------
!     subprogram 28. dump_read_lagr_edge.                               
!     real and imaginary parts are read separately for possible         
!     data conversion.                                                  
!-----------------------------------------------------------------------
      SUBROUTINE dump_read_lagr_edge(laq,lx,name,title) 
      USE io 
                                                                        
      INTEGER(i4), INTENT(IN) :: lx
      CHARACTER(*), INTENT(IN) :: name 
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title 
      TYPE(lagr_edge_type), INTENT(OUT) :: laq 
                                                
      INTEGER(i4) :: nqtmp,nftmp                        
      REAL(r8) :: iread 
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rread 
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: rread2 
                                                                        
      READ(rstrt_unit) iread 
      nqtmp=NINT(iread)
      READ(rstrt_unit) iread 
      nftmp=NINT(iread)
      READ(rstrt_unit) iread 
      laq%n_side=NINT(iread) 
      READ(rstrt_unit) iread 
      laq%n_int=NINT(iread) 
      CALL lagr_edge_alloc(laq,lx,nqtmp,nftmp,laq%n_side+1_i4)
      laq%name=name 
      IF (SIZE(title)<SIZE(laq%title)) THEN 
        laq%title=title(1) 
      ELSE 
        laq%title=title 
      ENDIF 
                                                                        
      ALLOCATE(rread(laq%nqty,0:lx,laq%nfour)) 
      READ(rstrt_unit) rread 
      laq%fs=rread 
      READ(rstrt_unit) rread 
      laq%fs=laq%fs+(0,1)*rread 
      DEALLOCATE(rread) 
                                                                        
      IF (laq%n_side>0) THEN 
        ALLOCATE(rread2(laq%nqty,laq%n_side,1:lx,laq%nfour)) 
        READ(rstrt_unit) rread2 
        laq%fsh=rread2 
        READ(rstrt_unit) rread2 
        laq%fsh=laq%fsh+(0,1)*rread2 
        DEALLOCATE(rread2) 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE dump_read_lagr_edge 
!-----------------------------------------------------------------------
!     subprogram 29. dump_read_lagr_edge_1D.                            
!     1D version.                                                       
!-----------------------------------------------------------------------
      SUBROUTINE dump_read_lagr_edge_1D(laq,lx,name,title) 
      USE io 
                                                                        
      INTEGER(i4), INTENT(IN) :: lx
      CHARACTER(*), INTENT(IN) :: name 
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title 
      TYPE(lagr_edge_1D_type), INTENT(OUT) :: laq 
        
      INTEGER(i4) :: nqtmp
      REAL(r8) :: iread 
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: rread 
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rread2 
                                                                        
      READ(rstrt_unit) iread
      nqtmp=NINT(iread)
      READ(rstrt_unit) iread 
      laq%n_side=NINT(iread) 
      READ(rstrt_unit) iread 
      laq%n_int=NINT(iread) 
      CALL lagr_edge_alloc(laq,lx,nqtmp,laq%n_side+1_i4) 
      laq%name=name 
      IF (SIZE(title)<SIZE(laq%title)) THEN 
        laq%title=title(1) 
      ELSE 
        laq%title=title 
      ENDIF 
                                                                        
      ALLOCATE(rread(laq%nqty,0:lx)) 
      READ(rstrt_unit) rread 
      laq%fs=rread 
      DEALLOCATE(rread) 
                                                                        
      IF (laq%n_side>0) THEN 
        ALLOCATE(rread2(laq%nqty,laq%n_side,1:lx)) 
        READ(rstrt_unit) rread2 
        laq%fsh=rread2 
        DEALLOCATE(rread2) 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE dump_read_lagr_edge_1D  
#ifdef HAVE_FC_HDF5
!-----------------------------------------------------------------------
!     subprogram 33. h5_dump_lagr_edge.                              
!     real and imaginary parts are witten separately for possible       
!     data conversion.                                                  
!-----------------------------------------------------------------------
      SUBROUTINE h5_dump_lagr_edge(laq,fname,gid,blid)
      USE io 
                                                                        
      TYPE(lagr_edge_type), INTENT(IN) :: laq 
      CHARACTER(*), INTENT(IN) :: fname,blid
      INTEGER(HID_T), INTENT(IN) :: gid

      INTEGER(i4) :: dn,nx,iq,ifr,ic,is
      REAL(r8), ALLOCATABLE :: tarr(:,:),ttarr(:,:)
      CHARACTER(64) :: mdname

      h5in%vsMD=TRIM(fname)
      mdname = TRIM(h5in%vsMD)

      dn=laq%n_side+1
      nx=laq%mx*dn
      ALLOCATE(tarr(0:nx,laq%nqty*laq%nfour),                      &
     &         ttarr(laq%nqty*laq%nfour,0:nx))
      DO iq=1,laq%nqty
        DO ifr=1,laq%nfour
          ic=iq+(ifr-1)*laq%nqty
          tarr(0:nx:dn,ic)=REAL(laq%fs(iq,:,ifr))
          IF (ALLOCATED(laq%fsh)) THEN
            DO is=1,laq%n_side
              tarr(is:nx:dn,ic)=REAL(laq%fsh(iq,is,:,ifr))
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      h5in%vsMD="Re_"//mdname
      DO iq=1,laq%nqty*laq%nfour
        ttarr(iq,0:nx)=tarr(0:nx,iq)
      ENDDO
      CALL dump_h5(gid,"re"//TRIM(fname)//TRIM(blid),ttarr,h5in,h5err)
      DO iq=1,laq%nqty
        DO ifr=1,laq%nfour
          ic=iq+(ifr-1)*laq%nqty
          tarr(0:nx:dn,ic)=AIMAG(laq%fs(iq,:,ifr))
          IF (ALLOCATED(laq%fsh)) THEN
            DO is=1,laq%n_side
              tarr(is:nx:dn,ic)=AIMAG(laq%fsh(iq,is,:,ifr))
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      h5in%vsMD="Im_"//mdname
      DO iq=1,laq%nqty*laq%nfour
        ttarr(iq,0:nx)=tarr(0:nx,iq)
      ENDDO
      CALL dump_h5(gid,"im"//TRIM(fname)//TRIM(blid),ttarr,h5in,h5err)
      DEALLOCATE(tarr,ttarr)
      h5in%vsMD=" "
                                                                        
      RETURN 
      END SUBROUTINE h5_dump_lagr_edge
!-----------------------------------------------------------------------
!     subprogram 34. h5_dump_lagr_quad_2D.                           
!     2D version.                                                       
!-----------------------------------------------------------------------
      SUBROUTINE h5_dump_lagr_edge_1D(laq,fname,gid,blid) 
      USE io 
                                                                        
      TYPE(lagr_edge_1D_type), INTENT(IN) :: laq 
      CHARACTER(*), INTENT(IN) :: fname,blid
      INTEGER(HID_T), INTENT(IN) :: gid

      INTEGER(i4) :: dn,nx,ny,iq,is
      REAL(r8), ALLOCATABLE :: tarr(:,:)

      h5in%vsMD=TRIM(fname)

      dn=laq%n_side+1
      nx=laq%mx*dn
      ALLOCATE(tarr(laq%nqty,0:nx))
      DO iq=1,laq%nqty
        tarr(iq,0:nx:dn)=laq%fs(iq,:)
        IF (ALLOCATED(laq%fsh)) THEN
          DO is=1,laq%n_side
            tarr(iq,is:nx:dn)=laq%fsh(iq,is,:)
          ENDDO
        ENDIF
      ENDDO

      CALL dump_h5(gid,TRIM(fname)//TRIM(blid),tarr,h5in,h5err)
      DEALLOCATE(tarr)
      h5in%vsMD=" "
                                                                        
      RETURN 
      END SUBROUTINE h5_dump_lagr_edge_1D 
!-----------------------------------------------------------------------
!     subprogram 35. h5_read_lagr_edge                               
!     real and imaginary parts are read separately for possible         
!     data conversion.                                                  
!-----------------------------------------------------------------------
      SUBROUTINE h5_read_lagr_edge(laq,lx,pd,nfour,nqty2,fname,gid,blid)
      USE io 
                                                                        
      INTEGER(i4), INTENT(IN) :: lx
      TYPE(lagr_edge_type), INTENT(OUT) :: laq 
      CHARACTER(*), INTENT(IN) :: fname,blid
      INTEGER(i4), INTENT(IN) :: pd,nfour,nqty2
      INTEGER(HID_T), INTENT(IN) :: gid
                                                                        
      INTEGER(i4) :: dn,nx,iq,ifr,ic,is,nqty,ii
      REAL(r8), ALLOCATABLE :: tarr(:,:),itarr(:,:),ttarr(:,:)
      INTEGER(HID_T) :: dsetid
      INTEGER :: error
      INTEGER(HSIZE_T), DIMENSION(2) :: dims

      IF (.NOT.obj_exists(gid,"re"//TRIM(fname)//TRIM(blid),h5err))     &
        THEN ! Allocate, set to zero and return
        CALL lagr_edge_alloc(laq,lx,nqty2,nfour,pd,TRIM(fname))  
        laq%fs=0._r8
        IF (ALLOCATED(laq%fsh)) laq%fsh=0._r8
        RETURN
      ENDIF

      CALL h5dopen_f(gid,"re"//TRIM(fname)//TRIM(blid),dsetid,error)
      CALL read_dims(dsetid,dims,h5err)
      CALL h5dclose_f(dsetid,error)
      nqty=dims(1)/nfour
      IF (nqty/=nqty2)                                                  &
     &  CALL nim_stop('unexpected lagr_edge nqty in '//fname)
      IF (dims(2)/=lx*pd+1.)                                            &
     &  CALL nim_stop('unexpected lagr_edge size in '//fname)

      CALL lagr_edge_alloc(laq,lx,nqty,nfour,pd,TRIM(fname))
      ii=len(fname)
      IF (ii>=6) THEN
        laq%name=TRIM(fname(1:6))
      ELSE
        laq%name=TRIM(fname)
      ENDIF
 
      dn=laq%n_side+1
      nx=laq%mx*dn
      ALLOCATE( tarr(0:nx,laq%nqty*laq%nfour),                     &
     &         itarr(0:nx,laq%nqty*laq%nfour),                     &
     &         ttarr(laq%nqty*laq%nfour,0:nx)) 
      CALL read_h5(gid,"re"//TRIM(fname)//TRIM(blid),ttarr,h5in,h5err)
      DO iq=1,laq%nqty*laq%nfour
        tarr(0:nx,iq)=ttarr(iq,0:nx)
      ENDDO
      CALL read_h5(gid,"im"//TRIM(fname)//TRIM(blid),ttarr,h5in,h5err)
      DO iq=1,laq%nqty*laq%nfour
        itarr(0:nx,iq)=ttarr(iq,0:nx)
      ENDDO
      DO iq=1,laq%nqty
        DO ifr=1,laq%nfour
          ic=iq+(ifr-1)*laq%nqty
          laq%fs(iq,:,ifr)=tarr(0:nx:dn,ic)                   &
     &                      +(0,1)*itarr(0:nx:dn,ic)
          IF (ALLOCATED(laq%fsh)) THEN
            DO is=1,laq%n_side
              laq%fsh(iq,is,:,ifr)=tarr(is:nx:dn,ic)          &
     &                              +(0,1)*itarr(is:nx:dn,ic)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE(tarr,itarr,ttarr)
                                                                        
      RETURN 
      END SUBROUTINE h5_read_lagr_edge
!-----------------------------------------------------------------------
!     subprogram 36. h5_read_lagr_edge_2D.                            
!     2D version.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE h5_read_lagr_edge_1D(laq,lx,pd,nqty2,fname,gid,blid)
      USE io 
                                                                        
      INTEGER(i4), INTENT(IN) :: lx
      TYPE(lagr_edge_1D_type), INTENT(OUT) :: laq 
      CHARACTER(*), INTENT(IN) :: fname,blid
      INTEGER(HID_T), INTENT(IN) :: gid
      INTEGER(i4), INTENT(IN) :: pd,nqty2
                                                                        
      INTEGER(i4) :: dn,nx,iq,is,ins,ine,ii,nqty
      REAL(r8), ALLOCATABLE :: tarr(:,:)
      INTEGER(HID_T) :: dsetid
      INTEGER :: error
      INTEGER(HSIZE_T), DIMENSION(2) :: dims

      IF (.NOT.obj_exists(gid,TRIM(fname)//TRIM(blid),h5err))           &
        THEN ! Allocate, set to zero and return
        CALL lagr_edge_alloc(laq,lx,nqty2,pd,TRIM(fname))  
        laq%fs=0._r8
        IF (ALLOCATED(laq%fsh)) laq%fsh=0._r8
        RETURN
      ENDIF

      CALL h5dopen_f(gid,TRIM(fname)//TRIM(blid),dsetid,error)
      CALL read_dims(dsetid,dims,h5err)
      CALL h5dclose_f(dsetid,error)
      nqty=dims(1)
      IF (nqty/=nqty2)                                                  &
     &  CALL nim_stop('unexpected 2D lagr_edge nqty in '//fname)
      IF (dims(2)/=lx*pd+1.)                         &
     &  CALL nim_stop('unexpected 2D lagr_edge size in '//fname)

      CALL lagr_edge_alloc(laq,lx,nqty,pd,TRIM(fname))
      ii=len(fname)
      IF (ii>=6) THEN
        laq%name=TRIM(fname(1:6))
      ELSE
        laq%name=TRIM(fname)
      ENDIF
       
      dn=laq%n_side+1
      nx=laq%mx*dn
      ALLOCATE(tarr(laq%nqty,0:nx)) 
      CALL read_h5(gid,TRIM(fname)//TRIM(blid),tarr,h5in,h5err)
      DO iq=1,laq%nqty
        laq%fs(iq,:)=tarr(iq,0:nx:dn)
        IF (ALLOCATED(laq%fsh)) THEN
          DO is=1,laq%n_side
            laq%fsh(iq,is,:)=tarr(iq,is:nx:dn)
          ENDDO
        ENDIF
      ENDDO
      DEALLOCATE(tarr)
                                                                        
      RETURN 
      END SUBROUTINE h5_read_lagr_edge_1D 
#endif /* HAVE_FC_HDF5 */
!-----------------------------------------------------------------------
!     close module.                                                     
!-----------------------------------------------------------------------
      END MODULE lagr_edge_mod
