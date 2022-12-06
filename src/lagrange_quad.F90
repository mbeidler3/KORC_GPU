!-----------------------------------------------------------------------
!     $Id: lagrange_quad.F90 6863 2019-08-19 23:55:22Z jking $
!     routines for evaluating Lagrange finite elements on blocks of     
!     structured quadrilaterals.  this package is based on the          
!     previous bilinear-only version.                                   
!                                                                       
!     throughout this package, basis functions are catagorized into node
!     (grid vertex), side, and interior for each element.               
!-----------------------------------------------------------------------
#include "config.f"
!-----------------------------------------------------------------------
!     code organization.                                                
!-----------------------------------------------------------------------
!     0. lagr_type_mod definition.
!     1. lagr_quad_mod definition.
!     2. lagr_quad_bases.                                    
!     3. lagr_quad_3D_alloc.                                                              
!     5. lagr_quad_3D_eval.                                  
!     5a. lagr_quad_3D_eval_ts. 
!     6. lagr_quad_3D_all_eval.                              
!     7. lagr_quad_3D_assign_rsc.                            
!     8. lagr_quad_3D_assign_csc.                            
!     9. lagr_quad_3D_assign_laq.                            
!     10. lagr_quad_3D_assign_int.                                                               
!     13. lagr_quad_3D_basis_assign_loc                                             
!     16. lagr_quad_2D_alloc.                                                             
!     18. lagr_quad_2D_eval.                                 
!     18a. lagr_quad_2D_eval_ts. 
!     19. lagr_quad_2D_all_eval.                             
!     20. lagr_quad_2D_assign_rsc.                           
!     21. lagr_quad_2D_assign_csc.                           
!     22. lagr_quad_2D_assign_laq.                           
!     23. lagr_quad_2D_assign_int.                                                                   
!     26. lagr_quad_2D_basis_assign_loc                                                                                            
!     31. dump_read_lagr_quad.                               
!     32. dump_read_lagr_quad_2D.                            
!     33. h5_dump_lagr_quad.                              
!     34. h5_dump_lagr_quad_2D.                           
!     35. h5_read_lagr_quad.                               
!     36. h5_read_lagr_quad_2D.   
!-----------------------------------------------------------------------
!     subprogram 0. lagr_type_mod definition.
!     defines the lagrangian quadrilateral types and may be used
!     separately from lagr_quad_mod.
!-----------------------------------------------------------------------
      MODULE lagr_type_mod 
      USE local 
      IMPLICIT NONE 
                                                                        
      TYPE :: lagr_quad_type 
        INTEGER(i4) :: mx,my,nqty,nfour,n_side,n_int 
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ix0,iy0 
        REAL(r8), DIMENSION(:), ALLOCATABLE :: dx,dy 
        COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: fs 
        COMPLEX(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: fsh,fsv,fsi 
        COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: f,fx,fy
        CHARACTER(6), DIMENSION(:), ALLOCATABLE :: title 
        CHARACTER(6) :: name 
        INTEGER(i4), POINTER :: mem_id
      END TYPE lagr_quad_type 
                                                                        
      TYPE :: lagr_quad_2D_type 
        INTEGER(i4) :: mx,my,nqty,n_side,n_int 
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ix0,iy0 
        REAL(r8), DIMENSION(:), ALLOCATABLE :: dx,dy 
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: fs 
        REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: fsh,fsv,fsi 
        REAL(r8), DIMENSION(:), ALLOCATABLE :: f,fx,fy
        CHARACTER(6), DIMENSION(:), ALLOCATABLE :: title 
        CHARACTER(6) :: name 
        INTEGER(i4), POINTER :: mem_id
      END TYPE lagr_quad_2D_type 
 
      END MODULE lagr_type_mod
!-----------------------------------------------------------------------
!     subprogram 1. lagr_quad_mod definition.
!     contains the subprograms for manipulating lagrange_quad expansions
!     in structured blocks of quadrlaters, where the bases are
!     continuous across element borders.
!-----------------------------------------------------------------------
      MODULE lagr_quad_mod
      USE lagr_type_mod
      IMPLICIT NONE

      INTEGER(i4), PARAMETER, PRIVATE :: npoly_max=20 
!-----------------------------------------------------------------------
!     subprogram name interfaces                                        
!-----------------------------------------------------------------------
      INTERFACE lagr_quad_alloc 
        MODULE PROCEDURE lagr_quad_2D_alloc,lagr_quad_3D_alloc 
      END INTERFACE 
                                                                                                                                                
      INTERFACE lagr_quad_all_eval 
        MODULE PROCEDURE lagr_quad_2D_all_eval,lagr_quad_3D_all_eval 
      END INTERFACE 
                                                                        
      INTERFACE lagr_quad_eval 
        MODULE PROCEDURE lagr_quad_2D_eval,lagr_quad_3D_eval,           &
     &    lagr_quad_2D_eval_ts,lagr_quad_3D_eval_ts
      END INTERFACE 
                                                                                                                                                                                                                      
      INTERFACE lagr_quad_basis_assign_loc 
        MODULE PROCEDURE                                                &
     &    lagr_quad_2D_basis_assign_loc,lagr_quad_3D_basis_assign_loc   
      END INTERFACE 
                                                                                                                                                
!-----------------------------------------------------------------------
!     overloaded assignment.                                            
!-----------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=) 
        MODULE PROCEDURE                                                &
     &    lagr_quad_3D_assign_csc,lagr_quad_3D_assign_rsc,              &
     &    lagr_quad_3D_assign_laq,lagr_quad_3D_assign_int,              &
     &    lagr_quad_2D_assign_csc,lagr_quad_2D_assign_rsc,              &
     &    lagr_quad_2D_assign_laq,lagr_quad_2D_assign_int               
      END INTERFACE 
                                                                        
      CONTAINS 
!-----------------------------------------------------------------------
!     subprogram 2. lagr_quad_bases.                                    
!     computes basis functions and their derivatives.                   
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_bases(x,y,alpha,alphax,alphay,dmode,alxo,    &
     &               alyo)
                                                                        
      REAL(r8), INTENT(IN) :: x,y 
      REAL(r8), DIMENSION(:), INTENT(INOUT) :: alpha,alphax,alphay 
      INTEGER(i4), INTENT(IN) :: dmode 
      REAL(r8), DIMENSION(:), OPTIONAL, INTENT(OUT) :: alxo,alyo
                                                                        
      INTEGER(i4) :: pd,i,j,k
      REAL(r8), DIMENSION(0:npoly_max) :: alx,aly,dalx,daly
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
      pd=NINT(SQRT(REAL(SIZE(alpha))))-1 
      CALL lagr_1D(pd,x,alx,dalx,dmode) 
      CALL lagr_1D(pd,y,aly,daly,dmode) 
      IF (PRESENT(alxo)) THEN
        alxo=alx
      ENDIF
      IF (PRESENT(alyo)) THEN
        alyo=aly
      ENDIF
!-----------------------------------------------------------------------
!     if derivatives are not needed, limit the computations.            
!-----------------------------------------------------------------------
      IF (dmode==0) THEN 
!-----------------------------------------------------------------------
!       compute basis functions with respect to the logical coordinates 
!       of the unit square.  node (grid vertex-centered) bases first;   
!       start in lower left corner, work left to right, bottom to top.  
!-----------------------------------------------------------------------
        alpha(1)=alx(0) *aly(0) 
        alpha(2)=alx(pd)*aly(0) 
        alpha(3)=alx(0) *aly(pd) 
        alpha(4)=alx(pd)*aly(pd) 
!-----------------------------------------------------------------------
!       horizontal side centered bases if needed.  the ordering is      
!       bottom to top, then left to right in pairs.                     
!-----------------------------------------------------------------------
        k=5 
        DO i=1,pd-1 
          alpha(k)= alx(i)* aly(0) 
          k=k+1 
          alpha(k)= alx(i)* aly(pd) 
          k=k+1 
        ENDDO 
!-----------------------------------------------------------------------
!       vertical side centered bases if needed.  the ordering is        
!       left to right, then bottom to top in pairs.                     
!-----------------------------------------------------------------------
        DO j=1,pd-1 
          alpha(k)= alx(0)* aly(j) 
          k=k+1 
          alpha(k)= alx(pd)* aly(j) 
          k=k+1 
        ENDDO 
!-----------------------------------------------------------------------
!       interior bases if needed.  the ordering is                      
!       left to right, bottom to top.                                   
!-----------------------------------------------------------------------
        DO j=1,pd-1 
          DO i=1,pd-1 
            alpha(k)= alx(i)* aly(j) 
            k=k+1 
          ENDDO 
        ENDDO 
      ELSE 
!-----------------------------------------------------------------------
!       compute basis functions and derivatives with respect to the     
!       logical coordinates of the unit square.  node (grid vertex-     
!       centered) bases first; start in lower left corner, work left to 
!       right, bottom to top.                                           
!-----------------------------------------------------------------------
        alpha(1)=alx(0) *aly(0) 
        alpha(2)=alx(pd)*aly(0) 
        alpha(3)=alx(0) *aly(pd) 
        alpha(4)=alx(pd)*aly(pd) 
        alphax(1)=dalx(0) *aly(0) 
        alphax(2)=dalx(pd)*aly(0) 
        alphax(3)=dalx(0) *aly(pd) 
        alphax(4)=dalx(pd)*aly(pd) 
        alphay(1)=alx(0) *daly(0) 
        alphay(2)=alx(pd)*daly(0) 
        alphay(3)=alx(0) *daly(pd) 
        alphay(4)=alx(pd)*daly(pd) 
!-----------------------------------------------------------------------
!       horizontal side centered bases if needed.  the ordering is      
!       bottom to top, then left to right in pairs.                     
!-----------------------------------------------------------------------
        k=5 
        DO i=1,pd-1 
          alpha (k)= alx(i)* aly(0) 
          alphax(k)=dalx(i)* aly(0) 
          alphay(k)= alx(i)*daly(0) 
          k=k+1 
          alpha (k)= alx(i)* aly(pd) 
          alphax(k)=dalx(i)* aly(pd) 
          alphay(k)= alx(i)*daly(pd) 
          k=k+1 
        ENDDO 
!-----------------------------------------------------------------------
!       vertical side centered bases if needed.  the ordering is        
!       left to right, then bottom to top in pairs.                     
!-----------------------------------------------------------------------
        DO j=1,pd-1 
          alpha (k)= alx(0)* aly(j) 
          alphax(k)=dalx(0)* aly(j) 
          alphay(k)= alx(0)*daly(j) 
          k=k+1 
          alpha (k)= alx(pd)* aly(j) 
          alphax(k)=dalx(pd)* aly(j) 
          alphay(k)= alx(pd)*daly(j) 
          k=k+1 
        ENDDO 
!-----------------------------------------------------------------------
!       interior bases if needed.  the ordering is                      
!       left to right, bottom to top.                                   
!-----------------------------------------------------------------------
        DO j=1,pd-1 
          DO i=1,pd-1 
            alpha (k)= alx(i)* aly(j) 
            alphax(k)=dalx(i)* aly(j) 
            alphay(k)= alx(i)*daly(j) 
            k=k+1 
          ENDDO 
        ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_bases 
!-----------------------------------------------------------------------
!     subprogram 3. lagr_quad_3D_alloc.                                 
!     allocates space for lagr_quad_3D_type.                            
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_alloc(laq,mx,my,nqty,nfour,poly_degree,   &
     &                              name,title) 
#ifdef OBJ_MEM_PROF
      USE memlog, only: memlogger
#endif                        
      INTEGER(i4), INTENT(IN) :: mx,my,nqty,nfour,poly_degree 
      CHARACTER(*), INTENT(IN), OPTIONAL :: name 
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title 
      TYPE(lagr_quad_type), INTENT(OUT) :: laq 
                                                                        
      INTEGER(i4) :: ix,iy,ib,ip,pd1,sz
      REAL(r8), DIMENSION(0:poly_degree) :: x_node 
      CHARACTER(16) :: obj_name
!-----------------------------------------------------------------------
!     throw an error if allocated
!-----------------------------------------------------------------------
      IF (ALLOCATED(laq%fs)) THEN
        CALL nim_stop("Tried to reallocate lagr_quad_2D_alloc "         &
     &                //TRIM(title(1)))
      ENDIF
!-----------------------------------------------------------------------
!     store grid, vector, and fourier series dimensions, and set the    
!     number of side and interior basis functions.                      
!-----------------------------------------------------------------------
      laq%mx=mx 
      laq%my=my 
      laq%nqty=nqty 
      laq%nfour=nfour 
      laq%n_side=poly_degree-1 
      laq%n_int=(poly_degree-1)**2 
!-----------------------------------------------------------------------
!     allocate space.                                                   
!-----------------------------------------------------------------------
      ALLOCATE(laq%fs(nqty,0:mx,0:my,nfour)) 
      IF (poly_degree>1) THEN 
        ALLOCATE(laq%fsh(nqty,laq%n_side,1:mx,0:my,nfour)) 
        ALLOCATE(laq%fsv(nqty,laq%n_side,0:mx,1:my,nfour)) 
        ALLOCATE(laq%fsi(nqty,laq%n_int,1:mx,1:my,nfour)) 
      ENDIF 
      ALLOCATE(laq%title(nqty)) 
      ALLOCATE(laq%f(nqty,nfour)) 
      ALLOCATE(laq%fx(nqty,nfour)) 
      ALLOCATE(laq%fy(nqty,nfour)) 
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
!     grid-vertex data has grid indices running (0:mx,0:my), so         
!       ix0=0 and iy0=0.                                                
!     horizontal side data has grid indices running (1:mx,0:my), so     
!       ix0=1 and iy0=0.                                                
!     vertical side data has grid indices running (0:mx,1:my), so       
!       ix0=0 and iy0=1.                                                
!     interior data has grid indices running (1:mx,1:my), so            
!       ix0=1 and iy0=1.                                                
!                                                                       
!     dx and dy give the relative positions within a cell for the       
!     different bases.  note that these quantities may be used to access
!     all data with loops like:                                         
!                                                                       
!     DO ibase=1,SIZE(laq%ix0)                                          
!       DO iy=laq%iy0(ibase),laq%my                                     
!         DO ix=laq%ix0(ibase),laq%mx                                   
!           dx_data=ix-laq%ix0(ibase)+laq%dx(ibase)                     
!           dy_data=iy-laq%iy0(ibase)+laq%dy(ibase)                     
!           CALL lagr_quad_eval(laq,dx_data,dy_data,1_i4)               
!           xxx=laq%f(yyy,zzz)                                          
!         ENDDO                                                         
!       ENDDO                                                           
!     ENDDO                                                             
!-----------------------------------------------------------------------
      pd1=poly_degree-1 
      ALLOCATE(laq%ix0(poly_degree**2)) 
      ALLOCATE(laq%iy0(poly_degree**2)) 
      ALLOCATE(laq%dx(poly_degree**2)) 
      ALLOCATE(laq%dy(poly_degree**2)) 
      CALL poly_nodes(poly_degree,x_node) 
      DO ib=1,poly_degree**2 
        IF (ib==1) THEN 
          laq%ix0(ib)=0 
          laq%iy0(ib)=0 
          laq%dx(ib)=0 
          laq%dy(ib)=0 
        ELSE IF (ib<=poly_degree) THEN 
          laq%ix0(ib)=1 
          laq%iy0(ib)=0 
          laq%dx(ib)=x_node(ib-1_i4) 
          laq%dy(ib)=0 
        ELSE IF (ib<2*poly_degree) THEN 
          laq%ix0(ib)=0 
          laq%iy0(ib)=1 
          laq%dx(ib)=0 
          laq%dy(ib)=x_node(ib-poly_degree) 
        ELSE 
          laq%ix0(ib)=1 
          laq%iy0(ib)=1 
          ip=ib-2*poly_degree 
          ix=MODULO(ip,pd1)+1 
          iy=ip/pd1+1 
          laq%dx(ib)=x_node(ix) 
          laq%dy(ib)=x_node(iy) 
        ENDIF 
      ENDDO 
!-----------------------------------------------------------------------
!     register this object.
!-----------------------------------------------------------------------
      NULLIFY(laq%mem_id)
#ifdef OBJ_MEM_PROF
      sz= SIZEOF(laq%fs)+SIZEOF(laq%f)+SIZEOF(laq%fx)+SIZEOF(laq%fy)    &
     &   +SIZEOF(laq%title)+SIZEOF(laq%ix0)+SIZEOF(laq%iy0)             &
     &   +SIZEOF(laq%dx)+SIZEOF(laq%dy)+SIZEOF(laq)
      IF(poly_degree>1) sz=sz+SIZEOF(laq%fsh)                           &
     &                       +SIZEOF(laq%fsv)+SIZEOF(laq%fsi)
      obj_name='unknown'
      IF (PRESENT(name)) obj_name=name
      CALL memlogger%update(laq%mem_id,'laq3D',obj_name,sz)
#endif
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_3D_alloc 
!-----------------------------------------------------------------------
!     subprogram 5. lagr_quad_3D_eval.                                  
!     evaluates complex lagr_quad quantities at a single point within a 
!     grid block.                                                       
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_eval(laq,x,y,dmode) 
                                                                        
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq 
      REAL(r8), INTENT(IN) :: x,y 
      INTEGER(i4), INTENT(IN) :: dmode 
                                                                        
      INTEGER(i4) :: ix,iy,pd,i,j,k,im 
      REAL(r8), DIMENSION(0:npoly_max) :: alx,aly,dalx,daly
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
      iy=MAX(MIN(INT(y),laq%my-1),0_i4) 
      pd=laq%n_side+1 
      CALL lagr_1D(pd,x-ix,alx,dalx,dmode) 
      CALL lagr_1D(pd,y-iy,aly,daly,dmode) 
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.  bilinear   
!     derivatives simplify.                                             
!-----------------------------------------------------------------------
      IF (pd==1) THEN 
        laq%f=(laq%fs(:,ix  ,iy  ,:)*alx(0)                             &
     &        +laq%fs(:,ix+1,iy  ,:)*alx(1))*aly(0)+                    &
     &        (laq%fs(:,ix  ,iy+1,:)*alx(0)                             &
     &        +laq%fs(:,ix+1,iy+1,:)*alx(1))*aly(1)                     
        IF (dmode==0) RETURN 
        laq%fx=(laq%fs(:,ix+1,iy  ,:)                                   &
     &         -laq%fs(:,ix  ,iy  ,:))*aly(0)+                          &
     &         (laq%fs(:,ix+1,iy+1,:)                                   &
     &         -laq%fs(:,ix  ,iy+1,:))*aly(1)                           
        laq%fy=(laq%fs(:,ix  ,iy+1,:)                                   &
     &         -laq%fs(:,ix  ,iy  ,:))*alx(0)+                          &
     &         (laq%fs(:,ix+1,iy+1,:)                                   &
     &         -laq%fs(:,ix+1,iy  ,:))*alx(1)                           
!-----------------------------------------------------------------------
!     other polynomials:                                                
!-----------------------------------------------------------------------
      ELSE 
        laq%f=(laq%fs(:,ix  ,iy  ,:)*alx(0)                             &
     &        +laq%fs(:,ix+1,iy  ,:)*alx(pd))*aly(0)+                   &
     &        (laq%fs(:,ix  ,iy+1,:)*alx(0)                             &
     &        +laq%fs(:,ix+1,iy+1,:)*alx(pd))*aly(pd)                   
        DO i=1,pd-1 
          laq%f=laq%f+                                                  &
     &          (laq%fsh(:,i,ix+1,iy  ,:)*aly(0)                        &
     &          +laq%fsh(:,i,ix+1,iy+1,:)*aly(pd))*alx(i)+              &
     &          (laq%fsv(:,i,ix  ,iy+1,:)*alx(0)                        &
     &          +laq%fsv(:,i,ix+1,iy+1,:)*alx(pd))*aly(i)               
        ENDDO 
        k=1 
        DO j=1,pd-1 
          DO i=1,pd-1 
            laq%f=laq%f+laq%fsi(:,k,ix+1,iy+1,:)*alx(i)*aly(j) 
            k=k+1 
          ENDDO 
        ENDDO 
        IF (dmode==0) RETURN 
        laq%fx=(laq%fs(:,ix  ,iy  ,:)*dalx(0)                           &
     &         +laq%fs(:,ix+1,iy  ,:)*dalx(pd))*aly(0)+                 &
     &         (laq%fs(:,ix  ,iy+1,:)*dalx(0)                           &
     &         +laq%fs(:,ix+1,iy+1,:)*dalx(pd))*aly(pd)                 
        laq%fy=(laq%fs(:,ix  ,iy  ,:)*alx(0)                            &
     &         +laq%fs(:,ix+1,iy  ,:)*alx(pd))*daly(0)+                 &
     &         (laq%fs(:,ix  ,iy+1,:)*alx(0)                            &
     &         +laq%fs(:,ix+1,iy+1,:)*alx(pd))*daly(pd)                 
        DO i=1,pd-1 
          laq%fx=laq%fx+                                                &
     &           (laq%fsh(:,i,ix+1,iy  ,:)* aly(0)                      &
     &           +laq%fsh(:,i,ix+1,iy+1,:)* aly(pd))*dalx(i)+           &
     &           (laq%fsv(:,i,ix  ,iy+1,:)*dalx(0)                      &
     &           +laq%fsv(:,i,ix+1,iy+1,:)*dalx(pd))*aly(i)             
          laq%fy=laq%fy+                                                &
     &           (laq%fsh(:,i,ix+1,iy  ,:)*daly(0)                      &
     &           +laq%fsh(:,i,ix+1,iy+1,:)*daly(pd))*alx(i)+            &
     &           (laq%fsv(:,i,ix  ,iy+1,:)* alx(0)                      &
     &           +laq%fsv(:,i,ix+1,iy+1,:)* alx(pd))*daly(i)            
        ENDDO 
        k=1 
        DO j=1,pd-1 
          DO i=1,pd-1 
            laq%fx=laq%fx+laq%fsi(:,k,ix+1,iy+1,:)*dalx(i)* aly(j) 
            laq%fy=laq%fy+laq%fsi(:,k,ix+1,iy+1,:)* alx(i)*daly(j) 
            k=k+1 
          ENDDO 
        ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_3D_eval 
!-----------------------------------------------------------------------
!     subprogram 5a. lagr_quad_3D_eval_ts. 
!     evaluates complex lagr_quad quantities at a single point within a 
!     grid block. This is a thread safe version.                              
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_eval_ts(laq,x,y,dmode,f,fx,fy,fxx,fxy,fyy) 
                                                                        
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq 
      REAL(r8), INTENT(IN) :: x,y 
      INTEGER(i4), INTENT(IN) :: dmode 
      COMPLEX(r8), DIMENSION(laq%nqty,laq%nfour), INTENT(OUT) :: f
      COMPLEX(r8), DIMENSION(laq%nqty,laq%nfour), INTENT(OUT),          &
     &             OPTIONAL :: fx,fy,fxx,fxy,fyy
                                                                        
      INTEGER(i4) :: ix,iy,pd,i,j,k,im 
      REAL(r8), DIMENSION(0:npoly_max) :: alx,aly,dalx,daly,d2alx,d2aly
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
      iy=MAX(MIN(INT(y),laq%my-1),0_i4) 
      pd=laq%n_side+1 
      CALL lagr_1D(pd,x-ix,alx,dalx,dmode,d2alx) 
      CALL lagr_1D(pd,y-iy,aly,daly,dmode,d2aly) 
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.  bilinear   
!     derivatives simplify.                                             
!-----------------------------------------------------------------------
      IF (pd==1) THEN
        f=(laq%fs(:,ix  ,iy  ,:)*alx(0)                                 &
     &        +laq%fs(:,ix+1,iy  ,:)*alx(1))*aly(0)+                    &
     &        (laq%fs(:,ix  ,iy+1,:)*alx(0)                             &
     &        +laq%fs(:,ix+1,iy+1,:)*alx(1))*aly(1)                     
        IF (dmode==0) RETURN 
        fx=(laq%fs(:,ix+1,iy  ,:)                                       &
     &         -laq%fs(:,ix  ,iy  ,:))*aly(0)+                          &
     &         (laq%fs(:,ix+1,iy+1,:)                                   &
     &         -laq%fs(:,ix  ,iy+1,:))*aly(1)                           
        fy=(laq%fs(:,ix  ,iy+1,:)                                       &
     &         -laq%fs(:,ix  ,iy  ,:))*alx(0)+                          &
     &         (laq%fs(:,ix+1,iy+1,:)                                   &
     &         -laq%fs(:,ix+1,iy  ,:))*alx(1)                           
!-----------------------------------------------------------------------
!     other polynomials:                                                
!-----------------------------------------------------------------------
      ELSE 
        f=(laq%fs(:,ix  ,iy  ,:)*alx(0)                                 &
     &    +laq%fs(:,ix+1,iy  ,:)*alx(pd))*aly(0)+                       &
     &    (laq%fs(:,ix  ,iy+1,:)*alx(0)                                 &
     &    +laq%fs(:,ix+1,iy+1,:)*alx(pd))*aly(pd)                   
        DO i=1,pd-1 
          f=f+                                                          &
     &       (laq%fsh(:,i,ix+1,iy  ,:)*aly(0)                           &
     &       +laq%fsh(:,i,ix+1,iy+1,:)*aly(pd))*alx(i)+                 &
     &       (laq%fsv(:,i,ix  ,iy+1,:)*alx(0)                           &
     &       +laq%fsv(:,i,ix+1,iy+1,:)*alx(pd))*aly(i)               
        ENDDO 
        k=1 
        DO j=1,pd-1 
          DO i=1,pd-1 
            f=f+laq%fsi(:,k,ix+1,iy+1,:)*alx(i)*aly(j) 
            k=k+1 
          ENDDO 
        ENDDO 
        IF (dmode==0) RETURN 
        fx=(laq%fs(:,ix  ,iy  ,:)*dalx(0)                               &
     &         +laq%fs(:,ix+1,iy  ,:)*dalx(pd))*aly(0)+                 &
     &         (laq%fs(:,ix  ,iy+1,:)*dalx(0)                           &
     &         +laq%fs(:,ix+1,iy+1,:)*dalx(pd))*aly(pd)                 
        fy=(laq%fs(:,ix  ,iy  ,:)*alx(0)                                &
     &         +laq%fs(:,ix+1,iy  ,:)*alx(pd))*daly(0)+                 &
     &         (laq%fs(:,ix  ,iy+1,:)*alx(0)                            &
     &         +laq%fs(:,ix+1,iy+1,:)*alx(pd))*daly(pd)                 
        DO i=1,pd-1 
          fx=fx+                                                        &
     &           (laq%fsh(:,i,ix+1,iy  ,:)* aly(0)                      &
     &           +laq%fsh(:,i,ix+1,iy+1,:)* aly(pd))*dalx(i)+           &
     &           (laq%fsv(:,i,ix  ,iy+1,:)*dalx(0)                      &
     &           +laq%fsv(:,i,ix+1,iy+1,:)*dalx(pd))*aly(i)             
          fy=fy+                                                        &
     &           (laq%fsh(:,i,ix+1,iy  ,:)*daly(0)                      &
     &           +laq%fsh(:,i,ix+1,iy+1,:)*daly(pd))*alx(i)+            &
     &           (laq%fsv(:,i,ix  ,iy+1,:)* alx(0)                      &
     &           +laq%fsv(:,i,ix+1,iy+1,:)* alx(pd))*daly(i)            
        ENDDO 
        k=1 
        DO j=1,pd-1 
          DO i=1,pd-1 
            fx=fx+laq%fsi(:,k,ix+1,iy+1,:)*dalx(i)* aly(j) 
            fy=fy+laq%fsi(:,k,ix+1,iy+1,:)* alx(i)*daly(j) 
            k=k+1 
          ENDDO 
        ENDDO 
        IF (dmode==1) RETURN 
        fxx=(laq%fs(:,ix  ,iy  ,:)*d2alx(0)                             &
     &      +laq%fs(:,ix+1,iy  ,:)*d2alx(pd))*aly(0)+                   &
     &      (laq%fs(:,ix  ,iy+1,:)*d2alx(0)                             &
     &      +laq%fs(:,ix+1,iy+1,:)*d2alx(pd))*aly(pd)               
        fxy=(laq%fs(:,ix  ,iy  ,:)*dalx(0)                              &
     &      +laq%fs(:,ix+1,iy  ,:)*dalx(pd))*daly(0)+                   &
     &      (laq%fs(:,ix  ,iy+1,:)*dalx(0)                              &
     &      +laq%fs(:,ix+1,iy+1,:)*dalx(pd))*daly(pd)              
        fyy=(laq%fs(:,ix  ,iy  ,:)*alx(0)                               &
     &      +laq%fs(:,ix+1,iy  ,:)*alx(pd))*d2aly(0)+                   &
     &      (laq%fs(:,ix  ,iy+1,:)*alx(0)                               &
     &      +laq%fs(:,ix+1,iy+1,:)*alx(pd))*d2aly(pd)               
        DO i=1,pd-1 
          fxx=fxx+                                                      &
     &        (laq%fsh(:,i,ix+1,iy  ,:)*  aly(0)                        &
     &        +laq%fsh(:,i,ix+1,iy+1,:)*  aly(pd))*d2alx(i)+            &
     &        (laq%fsv(:,i,ix  ,iy+1,:)*d2alx(0)                        &
     &        +laq%fsv(:,i,ix+1,iy+1,:)*d2alx(pd))*aly(i)          
          fxy=fxy+                                                      &
     &        (laq%fsh(:,i,ix+1,iy  ,:)*daly(0)                         &
     &        +laq%fsh(:,i,ix+1,iy+1,:)*daly(pd))*dalx(i)+              &
     &        (laq%fsv(:,i,ix  ,iy+1,:)*dalx(0)                         &
     &        +laq%fsv(:,i,ix+1,iy+1,:)*dalx(pd))*daly(i)          
          fyy=fyy+                                                      &
     &        (laq%fsh(:,i,ix+1,iy  ,:)*d2aly(0)                        &
     &        +laq%fsh(:,i,ix+1,iy+1,:)*d2aly(pd))*  alx(i)+            &
     &        (laq%fsv(:,i,ix  ,iy+1,:)*  alx(0)                        &
     &        +laq%fsv(:,i,ix+1,iy+1,:)*  alx(pd))*d2aly(i)         
        ENDDO 
        k=1 
        DO j=1,pd-1 
          DO i=1,pd-1 
            fxx=fxx+laq%fsi(:,k,ix+1,iy+1,:)*d2alx(i)*  aly(j) 
            fxy=fxy+laq%fsi(:,k,ix+1,iy+1,:)* dalx(i)* daly(j) 
            fyy=fyy+laq%fsi(:,k,ix+1,iy+1,:)*  alx(i)*d2aly(j) 
            k=k+1 
          ENDDO 
        ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_3D_eval_ts
!-----------------------------------------------------------------------
!     subprogram 6. lagr_quad_3D_all_eval.                              
!     evaluates complex lagr_quad quantities in all elements in a grid  
!     block for equal spacing.                                          
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_all_eval(laq,x,y,f,fx,fy,dmode) 
                                                                        
      TYPE(lagr_quad_type), INTENT(IN) :: laq 
      REAL(r8), INTENT(IN) :: x,y 
      COMPLEX(r8), INTENT(OUT), DIMENSION(:,:,:,:) :: f,fx,fy 
      INTEGER(i4), INTENT(IN) :: dmode 
                                                                        
      INTEGER(i4) :: mx,my,mx1,my1 
      INTEGER(i4) :: pd,i,j,k,im,ix,iy 
      REAL(r8), DIMENSION((laq%n_side+2)**2) :: alpha,dalpdx,dalpdy 
      REAL(r8), DIMENSION(0:npoly_max) :: alx,aly
!-----------------------------------------------------------------------
!     compute index limits.                                             
!-----------------------------------------------------------------------
      mx=laq%mx 
      mx1=mx-1 
      my=laq%my 
      my1=my-1 
      pd=laq%n_side+1 
      CALL lagr_quad_bases(x,y,alpha,dalpdx,dalpdy,dmode,alx,aly)
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.  bilinear   
!     derivatives simplify.  [do-loops are for optimization, not        
!     aesthetics.]                                                      
!-----------------------------------------------------------------------
      SELECT CASE(pd) 
      CASE(1) 
        IF (dmode==0) THEN 
          DO im=1,laq%nfour 
            DO iy=1,my 
              DO ix=1,mx 
                f(:,ix,iy,im)=laq%fs(:,ix-1,iy-1,im)*alpha(1)+          &
     &                        laq%fs(:,ix  ,iy-1,im)*alpha(2)+          &
     &                        laq%fs(:,ix-1,iy  ,im)*alpha(3)+          &
     &                        laq%fs(:,ix  ,iy  ,im)*alpha(4)           
              ENDDO 
            ENDDO 
          ENDDO 
        ELSE 
          DO im=1,laq%nfour 
            DO iy=1,my 
              DO ix=1,mx 
                f(:,ix,iy,im)=laq%fs(:,ix-1,iy-1,im)*alpha(1)+          &
     &                        laq%fs(:,ix  ,iy-1,im)*alpha(2)+          &
     &                        laq%fs(:,ix-1,iy  ,im)*alpha(3)+          &
     &                        laq%fs(:,ix  ,iy  ,im)*alpha(4)           
                fx(:,ix,iy,im)=(laq%fs(:,ix  ,iy-1,im)                  &
     &                         -laq%fs(:,ix-1,iy-1,im))*aly(0)+         &
     &                         (laq%fs(:,ix  ,iy  ,im)                  &
     &                         -laq%fs(:,ix-1,iy  ,im))*aly(1)          
                fy(:,ix,iy,im)=(laq%fs(:,ix-1,iy  ,im)                  &
     &                         -laq%fs(:,ix-1,iy-1,im))*alx(0)+         &
     &                         (laq%fs(:,ix  ,iy  ,im)                  &
     &                         -laq%fs(:,ix  ,iy-1,im))*alx(1)          
              ENDDO 
            ENDDO 
          ENDDO 
        ENDIF 
!-----------------------------------------------------------------------
!     biquadratics--no sums.                                            
!-----------------------------------------------------------------------
      CASE(2) 
        IF (dmode==0) THEN 
          DO im=1,laq%nfour 
            DO iy=1,my 
              DO ix=1,mx 
                f(:,ix,iy,im)=                                          &
     &               laq%fs(:,ix-1,iy-1,im)*alpha(1)+                   &
     &               laq%fs(:,ix  ,iy-1,im)*alpha(2)+                   &
     &               laq%fs(:,ix-1,iy  ,im)*alpha(3)+                   &
     &               laq%fs(:,ix  ,iy  ,im)*alpha(4)+                   &
     &               laq%fsh(:,1,ix,iy-1,im)*alpha(5)+                  &
     &               laq%fsh(:,1,ix,iy  ,im)*alpha(6)+                  &
     &               laq%fsv(:,1,ix-1,iy,im)*alpha(7)+                  &
     &               laq%fsv(:,1,ix  ,iy,im)*alpha(8)+                  &
     &               laq%fsi(:,1,ix,iy,im)*alpha(9)                     
              ENDDO 
            ENDDO 
          ENDDO 
        ELSE 
          DO im=1,laq%nfour 
            DO iy=1,my 
              DO ix=1,mx 
                f(:,ix,iy,im)=                                          &
     &               laq%fs(:,ix-1,iy-1,im)*alpha(1)+                   &
     &               laq%fs(:,ix  ,iy-1,im)*alpha(2)+                   &
     &               laq%fs(:,ix-1,iy  ,im)*alpha(3)+                   &
     &               laq%fs(:,ix  ,iy  ,im)*alpha(4)+                   &
     &               laq%fsh(:,1,ix,iy-1,im)*alpha(5)+                  &
     &               laq%fsh(:,1,ix,iy  ,im)*alpha(6)+                  &
     &               laq%fsv(:,1,ix-1,iy,im)*alpha(7)+                  &
     &               laq%fsv(:,1,ix  ,iy,im)*alpha(8)+                  &
     &               laq%fsi(:,1,ix,iy,im)*alpha(9)                     
                fx(:,ix,iy,im)=                                         &
     &                laq%fs(:,ix-1,iy-1,im)*dalpdx(1)+                 &
     &                laq%fs(:,ix  ,iy-1,im)*dalpdx(2)+                 &
     &                laq%fs(:,ix-1,iy  ,im)*dalpdx(3)+                 &
     &                laq%fs(:,ix  ,iy  ,im)*dalpdx(4)+                 &
     &                laq%fsh(:,1,ix,iy-1,im)*dalpdx(5)+                &
     &                laq%fsh(:,1,ix,iy  ,im)*dalpdx(6)+                &
     &                laq%fsv(:,1,ix-1,iy,im)*dalpdx(7)+                &
     &                laq%fsv(:,1,ix  ,iy,im)*dalpdx(8)+                &
     &                laq%fsi(:,1,ix,iy,im)*dalpdx(9)                   
                fy(:,ix,iy,im)=                                         &
     &                laq%fs(:,ix-1,iy-1,im)*dalpdy(1)+                 &
     &                laq%fs(:,ix  ,iy-1,im)*dalpdy(2)+                 &
     &                laq%fs(:,ix-1,iy  ,im)*dalpdy(3)+                 &
     &                laq%fs(:,ix  ,iy  ,im)*dalpdy(4)+                 &
     &                laq%fsh(:,1,ix,iy-1,im)*dalpdy(5)+                &
     &                laq%fsh(:,1,ix,iy  ,im)*dalpdy(6)+                &
     &                laq%fsv(:,1,ix-1,iy,im)*dalpdy(7)+                &
     &                laq%fsv(:,1,ix  ,iy,im)*dalpdy(8)+                &
     &                laq%fsi(:,1,ix,iy,im)*dalpdy(9)                   
              ENDDO 
            ENDDO 
          ENDDO 
        ENDIF 
!-----------------------------------------------------------------------
!     other polynomials:                                                
!-----------------------------------------------------------------------
      CASE DEFAULT 
        j=2*pd+2 
        k=4*pd 
        IF (dmode==0) THEN 
          DO im=1,laq%nfour 
            DO iy=1,my 
              DO ix=1,mx 
                DO i=1,laq%nqty 
                  f(i,ix,iy,im)=                                        &
     &                 laq%fs(i,ix-1,iy-1,im)*alpha(1)+                 &
     &                 laq%fs(i,ix  ,iy-1,im)*alpha(2)+                 &
     &                 laq%fs(i,ix-1,iy  ,im)*alpha(3)+                 &
     &                 laq%fs(i,ix  ,iy  ,im)*alpha(4)+                 &
     &             SUM(laq%fsh(i,:,ix,iy-1,im)*alpha(5:j-1:2)+          &
     &                 laq%fsh(i,:,ix,iy  ,im)*alpha(6:j  :2)+          &
     &                 laq%fsv(i,:,ix-1,iy,im)*alpha(j+1:k-1:2)+        &
     &                 laq%fsv(i,:,ix  ,iy,im)*alpha(j+2:k:2))+         &
     &             SUM(laq%fsi(i,:,ix,iy,im)*alpha(k+1:))               
                ENDDO 
              ENDDO 
            ENDDO 
          ENDDO 
        ELSE 
          DO im=1,laq%nfour 
            DO iy=1,my 
              DO ix=1,mx 
                DO i=1,laq%nqty 
                  f(i,ix,iy,im)=                                        &
     &                 laq%fs(i,ix-1,iy-1,im)*alpha(1)+                 &
     &                 laq%fs(i,ix  ,iy-1,im)*alpha(2)+                 &
     &                 laq%fs(i,ix-1,iy  ,im)*alpha(3)+                 &
     &                 laq%fs(i,ix  ,iy  ,im)*alpha(4)+                 &
     &             SUM(laq%fsh(i,:,ix,iy-1,im)*alpha(5:j-1:2)+          &
     &                 laq%fsh(i,:,ix,iy  ,im)*alpha(6:j  :2)+          &
     &                 laq%fsv(i,:,ix-1,iy,im)*alpha(j+1:k-1:2)+        &
     &                 laq%fsv(i,:,ix  ,iy,im)*alpha(j+2:k:2))+         &
     &             SUM(laq%fsi(i,:,ix,iy,im)*alpha(k+1:))               
                  fx(i,ix,iy,im)=                                       &
     &                  laq%fs(i,ix-1,iy-1,im)*dalpdx(1)+               &
     &                  laq%fs(i,ix  ,iy-1,im)*dalpdx(2)+               &
     &                  laq%fs(i,ix-1,iy  ,im)*dalpdx(3)+               &
     &                  laq%fs(i,ix  ,iy  ,im)*dalpdx(4)+               &
     &              SUM(laq%fsh(i,:,ix,iy-1,im)*dalpdx(5:j-1:2)+        &
     &                  laq%fsh(i,:,ix,iy  ,im)*dalpdx(6:j  :2)+        &
     &                  laq%fsv(i,:,ix-1,iy,im)*dalpdx(j+1:k-1:2)+      &
     &                  laq%fsv(i,:,ix  ,iy,im)*dalpdx(j+2:k:2))+       &
     &              SUM(laq%fsi(i,:,ix,iy,im)*dalpdx(k+1:))             
                  fy(i,ix,iy,im)=                                       &
     &                  laq%fs(i,ix-1,iy-1,im)*dalpdy(1)+               &
     &                  laq%fs(i,ix  ,iy-1,im)*dalpdy(2)+               &
     &                  laq%fs(i,ix-1,iy  ,im)*dalpdy(3)+               &
     &                  laq%fs(i,ix  ,iy  ,im)*dalpdy(4)+               &
     &              SUM(laq%fsh(i,:,ix,iy-1,im)*dalpdy(5:j-1:2)+        &
     &                  laq%fsh(i,:,ix,iy  ,im)*dalpdy(6:j  :2)+        &
     &                  laq%fsv(i,:,ix-1,iy,im)*dalpdy(j+1:k-1:2)+      &
     &                  laq%fsv(i,:,ix  ,iy,im)*dalpdy(j+2:k:2))+       &
     &              SUM(laq%fsi(i,:,ix,iy,im)*dalpdy(k+1:))             
                ENDDO 
              ENDDO 
            ENDDO 
          ENDDO 
        ENDIF 
      END SELECT 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_3D_all_eval 
!-----------------------------------------------------------------------
!     subprogram 7. lagr_quad_3D_assign_rsc.                            
!     assign a real scalar value to a complex lagrange quad structure.  
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_assign_rsc(laq,rscalar) 
                                                                        
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq 
      REAL(r8), INTENT(IN) :: rscalar 
                                                                        
      laq%fs=rscalar 
      IF (ALLOCATED(laq%fsh)) THEN 
        laq%fsh=rscalar 
        laq%fsv=rscalar 
        laq%fsi=rscalar 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_3D_assign_rsc 
!-----------------------------------------------------------------------
!     subprogram 8. lagr_quad_3D_assign_csc.                            
!     assign a complex scalar value to a complex lagrange quad          
!     structure.                                                        
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_assign_csc(laq,cscalar) 
                                                                        
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq 
      COMPLEX(r8), INTENT(IN) :: cscalar 
                                                                        
      laq%fs=cscalar 
      IF (ALLOCATED(laq%fsh)) THEN 
        laq%fsh=cscalar 
        laq%fsv=cscalar 
        laq%fsi=cscalar 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_3D_assign_csc 
!-----------------------------------------------------------------------
!     subprogram 9. lagr_quad_3D_assign_laq.                            
!     set one complex lagrange quad structure equal to another.         
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_assign_laq(laq1,laq2) 
                                                                        
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq1 
      TYPE(lagr_quad_type), INTENT(IN) :: laq2 
                                                                        
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.                   
!-----------------------------------------------------------------------
      laq1%fs=laq2%fs 
      IF (ALLOCATED(laq1%fsh)) THEN 
        laq1%fsh=laq2%fsh 
        laq1%fsv=laq2%fsv 
        laq1%fsi=laq2%fsi 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_3D_assign_laq 
!-----------------------------------------------------------------------
!     subprogram 10. lagr_quad_3D_assign_int.                            
!     assign a integer value to a complex lagrange quad structure.      
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_assign_int(laq,int) 
                                                                        
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq 
      INTEGER(i4), INTENT(IN) :: int 
                                                                        
      laq%fs=int 
      IF (ALLOCATED(laq%fsh)) THEN 
        laq%fsh=int 
        laq%fsv=int 
        laq%fsi=int 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_3D_assign_int 
!-----------------------------------------------------------------------
!     subprogram 13. lagr_quad_3D_basis_assign_loc                      
!     assign data into coefficient arrays for one basis function.       
!                                                                       
!     this is a local version of lagr_quad_basis_assign_arr, where the  
!     the data is located at a given poloidal and fourier indices       
!     triplet, only.                                                    
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_basis_assign_loc(laq,data,ibasis,ix,iy,im) 
                                                                        
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq 
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: data 
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy,im 
                                                                        
      INTEGER(i4) :: poly_degree 
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to     
!     the value of ibasis.                                              
!     grid vertex-centered data is first.                               
!-----------------------------------------------------------------------
      poly_degree=laq%n_side+1 
      IF (ibasis==1) THEN 
        laq%fs(:,ix,iy,im)=data 
        RETURN 
!-----------------------------------------------------------------------
!     horizontal sides.                                                 
!-----------------------------------------------------------------------
      ELSE IF (ibasis<=poly_degree) THEN 
        laq%fsh(:,ibasis-1,ix,iy,im)=data 
        RETURN 
!-----------------------------------------------------------------------
!     vertical sides.                                                   
!-----------------------------------------------------------------------
      ELSE IF (ibasis<2*poly_degree) THEN 
        laq%fsv(:,ibasis-poly_degree,ix,iy,im)=data 
        RETURN 
!-----------------------------------------------------------------------
!     interior bases.                                                   
!-----------------------------------------------------------------------
      ELSE 
        laq%fsi(:,ibasis-2*poly_degree+1,ix,iy,im)=data 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_3D_basis_assign_loc 
!-----------------------------------------------------------------------
!     subprogram 16. lagr_quad_2D_alloc.                                
!     allocates space for lagr_quad_2D_type.                            
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_alloc(laq,mx,my,nqty,poly_degree,         &
     &                              name,title)                         
#ifdef OBJ_MEM_PROF
      USE memlog, only: memlogger
#endif                        
                                                                        
      INTEGER(i4), INTENT(IN) :: mx,my,nqty,poly_degree 
      CHARACTER(*), INTENT(IN), OPTIONAL :: name 
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title 
      TYPE(lagr_quad_2D_type), INTENT(OUT) :: laq 
                                                                        
      INTEGER(i4) :: ix,iy,ib,ip,pd1,sz
      REAL(r8), DIMENSION(0:poly_degree) :: x_node 
      CHARACTER(16) :: obj_name
!-----------------------------------------------------------------------
!     throw an error if allocated
!-----------------------------------------------------------------------
      IF (ALLOCATED(laq%fs)) THEN
        CALL nim_stop("Tried to reallocate lagr_quad_2D_alloc "         &
     &                //TRIM(title(1)))
      ENDIF
!-----------------------------------------------------------------------
!     store grid and vector dimensions, and set the                     
!     number of side and interior basis functions.                      
!-----------------------------------------------------------------------
      laq%mx=mx 
      laq%my=my 
      laq%nqty=nqty 
      laq%n_side=poly_degree-1 
      laq%n_int=(poly_degree-1)**2 
!-----------------------------------------------------------------------
!     allocate space.                                                   
!-----------------------------------------------------------------------
      ALLOCATE(laq%fs(nqty,0:mx,0:my)) 
      IF (poly_degree>1) THEN 
        ALLOCATE(laq%fsh(nqty,laq%n_side,1:mx,0:my)) 
        ALLOCATE(laq%fsv(nqty,laq%n_side,0:mx,1:my)) 
        ALLOCATE(laq%fsi(nqty,laq%n_int,1:mx,1:my)) 
      ENDIF 
      ALLOCATE(laq%title(nqty)) 
      ALLOCATE(laq%f(nqty)) 
      ALLOCATE(laq%fx(nqty)) 
      ALLOCATE(laq%fy(nqty)) 
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
!     grid-vertex data has grid indices running (0:mx,0:my), so         
!       ix0=0 and iy0=0.                                                
!     horizontal side data has grid indices running (1:mx,0:my), so     
!       ix0=1 and iy0=0.                                                
!     vertical side data has grid indices running (0:mx,1:my), so       
!       ix0=0 and iy0=1.                                                
!     interior data has grid indices running (1:mx,1:my), so            
!       ix0=1 and iy0=1.                                                
!                                                                       
!     dx and dy give the relative positions within a cell for the       
!     different bases.  note that these quantities may be used to access
!     all data with loops like:                                         
!                                                                       
!     DO ibase=1,SIZE(laq%ix0)                                          
!       DO iy=laq%iy0(ibase),laq%my                                     
!         DO ix=laq%ix0(ibase),laq%mx                                   
!           dx_data=ix-laq%ix0(ibase)+laq%dx(ibase)                     
!           dy_data=iy-laq%iy0(ibase)+laq%dy(ibase)                     
!           CALL lagr_quad_eval(laq,dx_data,dy_data,1_i4)               
!           xxx=laq%f(yyy)                                              
!         ENDDO                                                         
!       ENDDO                                                           
!     ENDDO                                                             
!-----------------------------------------------------------------------
      pd1=poly_degree-1 
      ALLOCATE(laq%ix0(poly_degree**2)) 
      ALLOCATE(laq%iy0(poly_degree**2)) 
      ALLOCATE(laq%dx(poly_degree**2)) 
      ALLOCATE(laq%dy(poly_degree**2)) 
      CALL poly_nodes(poly_degree,x_node) 
      DO ib=1,poly_degree**2 
        IF (ib==1) THEN 
          laq%ix0(ib)=0 
          laq%iy0(ib)=0 
          laq%dx(ib)=0 
          laq%dy(ib)=0 
        ELSE IF (ib<=poly_degree) THEN 
          laq%ix0(ib)=1 
          laq%iy0(ib)=0 
          laq%dx(ib)=x_node(ib-1_i4) 
          laq%dy(ib)=0 
        ELSE IF (ib<2*poly_degree) THEN 
          laq%ix0(ib)=0 
          laq%iy0(ib)=1 
          laq%dx(ib)=0 
          laq%dy(ib)=x_node(ib-poly_degree) 
        ELSE 
          laq%ix0(ib)=1 
          laq%iy0(ib)=1 
          ip=ib-2*poly_degree 
          ix=MODULO(ip,pd1)+1 
          iy=ip/pd1+1 
          laq%dx(ib)=x_node(ix) 
          laq%dy(ib)=x_node(iy) 
        ENDIF 
      ENDDO 
!-----------------------------------------------------------------------
!     register this object.
!-----------------------------------------------------------------------
      NULLIFY(laq%mem_id)
#ifdef OBJ_MEM_PROF
      sz= SIZEOF(laq%fs)+SIZEOF(laq%f)+SIZEOF(laq%fx)+SIZEOF(laq%fy)    &
     &   +SIZEOF(laq%title)+SIZEOF(laq%ix0)+SIZEOF(laq%iy0)             &
     &   +SIZEOF(laq%dx)+SIZEOF(laq%dy)+SIZEOF(laq)
      IF(poly_degree>1) sz=sz+SIZEOF(laq%fsh)                           &
     &                       +SIZEOF(laq%fsv)+SIZEOF(laq%fsi)
      obj_name='unknown'
      IF (PRESENT(name)) obj_name=name
      CALL memlogger%update(laq%mem_id,'laq2D',obj_name,sz)
#endif
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_2D_alloc 
!-----------------------------------------------------------------------
!     subprogram 18. lagr_quad_2D_eval.                                 
!     evaluates real lagr_quad quantities at a single point within a    
!     grid block.                                                       
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_eval(laq,x,y,dmode) 
                                                                        
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq 
      REAL(r8), INTENT(IN) :: x,y 
      INTEGER(i4), INTENT(IN) :: dmode 
                                                                        
      INTEGER(i4) :: ix,iy,pd,i,j,k,im 
      REAL(r8), DIMENSION(0:npoly_max) :: alx,aly,dalx,daly
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
      iy=MAX(MIN(INT(y),laq%my-1),0_i4) 
      pd=laq%n_side+1 
      CALL lagr_1D(pd,x-ix,alx,dalx,dmode) 
      CALL lagr_1D(pd,y-iy,aly,daly,dmode) 
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.  bilinear   
!     derivatives simplify.                                             
!-----------------------------------------------------------------------
      IF (pd==1) THEN 
        laq%f=(laq%fs(:,ix  ,iy  )*alx(0)                               &
     &        +laq%fs(:,ix+1,iy  )*alx(1))*aly(0)+                      &
     &        (laq%fs(:,ix  ,iy+1)*alx(0)                               &
     &        +laq%fs(:,ix+1,iy+1)*alx(1))*aly(1)                       
        IF (dmode==0) RETURN 
        laq%fx=(laq%fs(:,ix+1,iy  )                                     &
     &         -laq%fs(:,ix  ,iy  ))*aly(0)+                            &
     &         (laq%fs(:,ix+1,iy+1)                                     &
     &         -laq%fs(:,ix  ,iy+1))*aly(1)                             
        laq%fy=(laq%fs(:,ix  ,iy+1)                                     &
     &         -laq%fs(:,ix  ,iy  ))*alx(0)+                            &
     &         (laq%fs(:,ix+1,iy+1)                                     &
     &         -laq%fs(:,ix+1,iy  ))*alx(1)                             
!-----------------------------------------------------------------------
!     other polynomials:                                                
!-----------------------------------------------------------------------
      ELSE 
        laq%f=(laq%fs(:,ix  ,iy  )*alx(0)                               &
     &        +laq%fs(:,ix+1,iy  )*alx(pd))*aly(0)+                     &
     &        (laq%fs(:,ix  ,iy+1)*alx(0)                               &
     &        +laq%fs(:,ix+1,iy+1)*alx(pd))*aly(pd)                     
        DO i=1,pd-1 
          laq%f=laq%f+                                                  &
     &          (laq%fsh(:,i,ix+1,iy  )*aly(0)                          &
     &          +laq%fsh(:,i,ix+1,iy+1)*aly(pd))*alx(i)+                &
     &          (laq%fsv(:,i,ix  ,iy+1)*alx(0)                          &
     &          +laq%fsv(:,i,ix+1,iy+1)*alx(pd))*aly(i)                 
        ENDDO 
        k=1 
        DO j=1,pd-1 
          DO i=1,pd-1 
            laq%f=laq%f+laq%fsi(:,k,ix+1,iy+1)*alx(i)*aly(j) 
            k=k+1 
          ENDDO 
        ENDDO 
        IF (dmode==0) RETURN 
        laq%fx=(laq%fs(:,ix  ,iy  )*dalx(0)                             &
     &         +laq%fs(:,ix+1,iy  )*dalx(pd))*aly(0)+                   &
     &         (laq%fs(:,ix  ,iy+1)*dalx(0)                             &
     &         +laq%fs(:,ix+1,iy+1)*dalx(pd))*aly(pd)                   
        laq%fy=(laq%fs(:,ix  ,iy  )*alx(0)                              &
     &         +laq%fs(:,ix+1,iy  )*alx(pd))*daly(0)+                   &
     &         (laq%fs(:,ix  ,iy+1)*alx(0)                              &
     &         +laq%fs(:,ix+1,iy+1)*alx(pd))*daly(pd)                   
        DO i=1,pd-1 
          laq%fx=laq%fx+                                                &
     &           (laq%fsh(:,i,ix+1,iy  )* aly(0)                        &
     &           +laq%fsh(:,i,ix+1,iy+1)* aly(pd))*dalx(i)+             &
     &           (laq%fsv(:,i,ix  ,iy+1)*dalx(0)                        &
     &           +laq%fsv(:,i,ix+1,iy+1)*dalx(pd))*aly(i)               
          laq%fy=laq%fy+                                                &
     &           (laq%fsh(:,i,ix+1,iy  )*daly(0)                        &
     &           +laq%fsh(:,i,ix+1,iy+1)*daly(pd))*alx(i)+              &
     &           (laq%fsv(:,i,ix  ,iy+1)* alx(0)                        &
     &           +laq%fsv(:,i,ix+1,iy+1)* alx(pd))*daly(i)              
        ENDDO 
        k=1 
        DO j=1,pd-1 
          DO i=1,pd-1 
            laq%fx=laq%fx+laq%fsi(:,k,ix+1,iy+1)*dalx(i)* aly(j) 
            laq%fy=laq%fy+laq%fsi(:,k,ix+1,iy+1)* alx(i)*daly(j) 
            k=k+1 
          ENDDO 
        ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_2D_eval
!-----------------------------------------------------------------------
!     subprogram 18a. lagr_quad_2D_eval_ts. 
!     evaluates real lagr_quad quantities at a single point within a    
!     grid block. This version is thread safe.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_eval_ts(laq,x,y,dmode,f,fx,fy,fxx,fxy,fyy) 
                                                                        
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq 
      REAL(r8), INTENT(IN) :: x,y 
      INTEGER(i4), INTENT(IN) :: dmode 
      REAL(r8), DIMENSION(laq%nqty), INTENT(OUT) :: f
      REAL(r8), DIMENSION(laq%nqty), INTENT(OUT), OPTIONAL :: fx,fy,    &
     &                                            fxx,fxy,fyy
                                                                        
      INTEGER(i4) :: ix,iy,pd,i,j,k,im 
      REAL(r8), DIMENSION(0:npoly_max) :: alx,aly,dalx,daly,d2alx,d2aly
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
      iy=MAX(MIN(INT(y),laq%my-1),0_i4) 
      pd=laq%n_side+1 
      CALL lagr_1D(pd,x-ix,alx,dalx,dmode,d2alx) 
      CALL lagr_1D(pd,y-iy,aly,daly,dmode,d2aly) 
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.  bilinear   
!     derivatives simplify.                                             
!-----------------------------------------------------------------------
      IF (pd==1) THEN 
        f=(laq%fs(:,ix  ,iy  )*alx(0)                                   &
     &        +laq%fs(:,ix+1,iy  )*alx(1))*aly(0)+                      &
     &        (laq%fs(:,ix  ,iy+1)*alx(0)                               &
     &        +laq%fs(:,ix+1,iy+1)*alx(1))*aly(1)                       
        IF (dmode==0) RETURN 
        fx=(laq%fs(:,ix+1,iy  )                                         &
     &         -laq%fs(:,ix  ,iy  ))*aly(0)+                            &
     &         (laq%fs(:,ix+1,iy+1)                                     &
     &         -laq%fs(:,ix  ,iy+1))*aly(1)                             
        fy=(laq%fs(:,ix  ,iy+1)                                         &
     &         -laq%fs(:,ix  ,iy  ))*alx(0)+                            &
     &         (laq%fs(:,ix+1,iy+1)                                     &
     &         -laq%fs(:,ix+1,iy  ))*alx(1)                             
!-----------------------------------------------------------------------
!     other polynomials:                                                
!-----------------------------------------------------------------------
      ELSE 
        f=(laq%fs(:,ix  ,iy  )*alx(0)                                   &
     &        +laq%fs(:,ix+1,iy  )*alx(pd))*aly(0)+                     &
     &        (laq%fs(:,ix  ,iy+1)*alx(0)                               &
     &        +laq%fs(:,ix+1,iy+1)*alx(pd))*aly(pd)                     
        DO i=1,pd-1 
          f=f+                                                          &
     &          (laq%fsh(:,i,ix+1,iy  )*aly(0)                          &
     &          +laq%fsh(:,i,ix+1,iy+1)*aly(pd))*alx(i)+                &
     &          (laq%fsv(:,i,ix  ,iy+1)*alx(0)                          &
     &          +laq%fsv(:,i,ix+1,iy+1)*alx(pd))*aly(i)                 
        ENDDO 
        k=1 
        DO j=1,pd-1 
          DO i=1,pd-1 
            f=f+laq%fsi(:,k,ix+1,iy+1)*alx(i)*aly(j) 
            k=k+1 
          ENDDO 
        ENDDO 
        IF (dmode==0) RETURN 
        fx=(laq%fs(:,ix  ,iy  )*dalx(0)                                 &
     &         +laq%fs(:,ix+1,iy  )*dalx(pd))*aly(0)+                   &
     &         (laq%fs(:,ix  ,iy+1)*dalx(0)                             &
     &         +laq%fs(:,ix+1,iy+1)*dalx(pd))*aly(pd)                   
        fy=(laq%fs(:,ix  ,iy  )*alx(0)                                  &
     &         +laq%fs(:,ix+1,iy  )*alx(pd))*daly(0)+                   &
     &         (laq%fs(:,ix  ,iy+1)*alx(0)                              &
     &         +laq%fs(:,ix+1,iy+1)*alx(pd))*daly(pd)                   
        DO i=1,pd-1 
          fx=fx+                                                        &
     &           (laq%fsh(:,i,ix+1,iy  )* aly(0)                        &
     &           +laq%fsh(:,i,ix+1,iy+1)* aly(pd))*dalx(i)+             &
     &           (laq%fsv(:,i,ix  ,iy+1)*dalx(0)                        &
     &           +laq%fsv(:,i,ix+1,iy+1)*dalx(pd))*aly(i)               
          fy=fy+                                                        &
     &           (laq%fsh(:,i,ix+1,iy  )*daly(0)                        &
     &           +laq%fsh(:,i,ix+1,iy+1)*daly(pd))*alx(i)+              &
     &           (laq%fsv(:,i,ix  ,iy+1)* alx(0)                        &
     &           +laq%fsv(:,i,ix+1,iy+1)* alx(pd))*daly(i)              
        ENDDO 
        k=1 
        DO j=1,pd-1 
          DO i=1,pd-1 
            fx=fx+laq%fsi(:,k,ix+1,iy+1)*dalx(i)* aly(j) 
            fy=fy+laq%fsi(:,k,ix+1,iy+1)* alx(i)*daly(j) 
            k=k+1 
          ENDDO 
        ENDDO 
        IF (dmode==1) RETURN 
        fxx=(laq%fs(:,ix  ,iy  )*d2alx(0)                               &
     &      +laq%fs(:,ix+1,iy  )*d2alx(pd))*aly(0)+                     &
     &      (laq%fs(:,ix  ,iy+1)*d2alx(0)                               &
     &      +laq%fs(:,ix+1,iy+1)*d2alx(pd))*aly(pd)                
        fxy=(laq%fs(:,ix  ,iy  )*dalx(0)                                &
     &      +laq%fs(:,ix+1,iy  )*dalx(pd))*daly(0)+                     &
     &      (laq%fs(:,ix  ,iy+1)*dalx(0)                                &
     &      +laq%fs(:,ix+1,iy+1)*dalx(pd))*daly(pd)              
        fyy=(laq%fs(:,ix  ,iy  )*alx(0)                                 &
     &      +laq%fs(:,ix+1,iy  )*alx(pd))*d2aly(0)+                     &
     &      (laq%fs(:,ix  ,iy+1)*alx(0)                                 &
     &      +laq%fs(:,ix+1,iy+1)*alx(pd))*d2aly(pd)               
        DO i=1,pd-1 
          fxx=fxx+                                                      &
     &        (laq%fsh(:,i,ix+1,iy  )*  aly(0)                          &
     &        +laq%fsh(:,i,ix+1,iy+1)*  aly(pd))*d2alx(i)+              &
     &        (laq%fsv(:,i,ix  ,iy+1)*d2alx(0)                          &
     &        +laq%fsv(:,i,ix+1,iy+1)*d2alx(pd))*aly(i)          
          fxy=fxy+                                                      &
     &        (laq%fsh(:,i,ix+1,iy  )*daly(0)                           &
     &        +laq%fsh(:,i,ix+1,iy+1)*daly(pd))*dalx(i)+                &
     &        (laq%fsv(:,i,ix  ,iy+1)*dalx(0)                           &
     &        +laq%fsv(:,i,ix+1,iy+1)*dalx(pd))*daly(i)          
          fyy=fyy+                                                      &
     &        (laq%fsh(:,i,ix+1,iy  )*d2aly(0)                          &
     &        +laq%fsh(:,i,ix+1,iy+1)*d2aly(pd))*  alx(i)+              &
     &        (laq%fsv(:,i,ix  ,iy+1)*  alx(0)                          &
     &        +laq%fsv(:,i,ix+1,iy+1)*  alx(pd))*d2aly(i)         
        ENDDO 
        k=1 
        DO j=1,pd-1 
          DO i=1,pd-1 
            fxx=fxx+laq%fsi(:,k,ix+1,iy+1)*d2alx(i)*  aly(j) 
            fxy=fxy+laq%fsi(:,k,ix+1,iy+1)* dalx(i)* daly(j) 
            fyy=fyy+laq%fsi(:,k,ix+1,iy+1)*  alx(i)*d2aly(j) 
            k=k+1 
          ENDDO 
        ENDDO 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_2D_eval_ts
!-----------------------------------------------------------------------
!     subprogram 19. lagr_quad_2D_all_eval.                             
!     evaluates real lagr_quad quantities in all elements in a grid     
!     block for equal spacing.                                          
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_all_eval(laq,x,y,f,fx,fy,dmode) 
                                                                        
      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq 
      REAL(r8), INTENT(IN) :: x,y 
      REAL(r8), INTENT(OUT), DIMENSION(:,:,:) :: f,fx,fy 
      INTEGER(i4), INTENT(IN) :: dmode 
                                                                        
      INTEGER(i4) :: mx,my,mx1,my1 
      INTEGER(i4) :: pd,i,j,k,ix,iy 
      REAL(r8), DIMENSION((laq%n_side+2)**2) :: alpha,dalpdx,dalpdy 
      REAL(r8), DIMENSION(0:npoly_max) :: alx,aly
!-----------------------------------------------------------------------
!     compute index limits.                                             
!-----------------------------------------------------------------------
      mx=laq%mx 
      mx1=mx-1 
      my=laq%my 
      my1=my-1 
      pd=laq%n_side+1 
      CALL lagr_quad_bases(x,y,alpha,dalpdx,dalpdy,dmode,alx,aly)
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.  bilinear   
!     derivatives simplify.  [do-loops are for optimization, not        
!     aesthetics.]                                                      
!-----------------------------------------------------------------------
      SELECT CASE(pd) 
      CASE(1) 
        IF (dmode==0) THEN 
          DO iy=1,my 
            DO ix=1,mx 
              f(:,ix,iy)=laq%fs(:,ix-1,iy-1)*alpha(1)+                  &
     &                   laq%fs(:,ix  ,iy-1)*alpha(2)+                  &
     &                   laq%fs(:,ix-1,iy  )*alpha(3)+                  &
     &                   laq%fs(:,ix  ,iy  )*alpha(4)                   
            ENDDO 
          ENDDO 
        ELSE 
          DO iy=1,my 
            DO ix=1,mx 
              f(:,ix,iy)=laq%fs(:,ix-1,iy-1)*alpha(1)+                  &
     &                   laq%fs(:,ix  ,iy-1)*alpha(2)+                  &
     &                   laq%fs(:,ix-1,iy  )*alpha(3)+                  &
     &                   laq%fs(:,ix  ,iy  )*alpha(4)                   
              fx(:,ix,iy)=(laq%fs(:,ix  ,iy-1)                          &
     &                    -laq%fs(:,ix-1,iy-1))*aly(0)+                 &
     &                    (laq%fs(:,ix  ,iy  )                          &
     &                    -laq%fs(:,ix-1,iy  ))*aly(1)                  
              fy(:,ix,iy)=(laq%fs(:,ix-1,iy  )                          &
     &                    -laq%fs(:,ix-1,iy-1))*alx(0)+                 &
     &                    (laq%fs(:,ix  ,iy  )                          &
     &                    -laq%fs(:,ix  ,iy-1))*alx(1)                  
            ENDDO 
          ENDDO 
        ENDIF 
!-----------------------------------------------------------------------
!     biquadratics--no sums.                                            
!-----------------------------------------------------------------------
      CASE(2) 
        IF (dmode==0) THEN 
          DO iy=1,my 
            DO ix=1,mx 
              f(:,ix,iy)=                                               &
     &             laq%fs(:,ix-1,iy-1)*alpha(1)+                        &
     &             laq%fs(:,ix  ,iy-1)*alpha(2)+                        &
     &             laq%fs(:,ix-1,iy  )*alpha(3)+                        &
     &             laq%fs(:,ix  ,iy  )*alpha(4)+                        &
     &             laq%fsh(:,1,ix,iy-1)*alpha(5)+                       &
     &             laq%fsh(:,1,ix,iy  )*alpha(6)+                       &
     &             laq%fsv(:,1,ix-1,iy)*alpha(7)+                       &
     &             laq%fsv(:,1,ix  ,iy)*alpha(8)+                       &
     &             laq%fsi(:,1,ix,iy)*alpha(9)                          
            ENDDO 
          ENDDO 
        ELSE 
          DO iy=1,my 
            DO ix=1,mx 
              f(:,ix,iy)=                                               &
     &             laq%fs(:,ix-1,iy-1)*alpha(1)+                        &
     &             laq%fs(:,ix  ,iy-1)*alpha(2)+                        &
     &             laq%fs(:,ix-1,iy  )*alpha(3)+                        &
     &             laq%fs(:,ix  ,iy  )*alpha(4)+                        &
     &             laq%fsh(:,1,ix,iy-1)*alpha(5)+                       &
     &             laq%fsh(:,1,ix,iy  )*alpha(6)+                       &
     &             laq%fsv(:,1,ix-1,iy)*alpha(7)+                       &
     &             laq%fsv(:,1,ix  ,iy)*alpha(8)+                       &
     &             laq%fsi(:,1,ix,iy)*alpha(9)                          
              fx(:,ix,iy)=                                              &
     &              laq%fs(:,ix-1,iy-1)*dalpdx(1)+                      &
     &              laq%fs(:,ix  ,iy-1)*dalpdx(2)+                      &
     &              laq%fs(:,ix-1,iy  )*dalpdx(3)+                      &
     &              laq%fs(:,ix  ,iy  )*dalpdx(4)+                      &
     &              laq%fsh(:,1,ix,iy-1)*dalpdx(5)+                     &
     &              laq%fsh(:,1,ix,iy  )*dalpdx(6)+                     &
     &              laq%fsv(:,1,ix-1,iy)*dalpdx(7)+                     &
     &              laq%fsv(:,1,ix  ,iy)*dalpdx(8)+                     &
     &              laq%fsi(:,1,ix,iy)*dalpdx(9)                        
              fy(:,ix,iy)=                                              &
     &              laq%fs(:,ix-1,iy-1)*dalpdy(1)+                      &
     &              laq%fs(:,ix  ,iy-1)*dalpdy(2)+                      &
     &              laq%fs(:,ix-1,iy  )*dalpdy(3)+                      &
     &              laq%fs(:,ix  ,iy  )*dalpdy(4)+                      &
     &              laq%fsh(:,1,ix,iy-1)*dalpdy(5)+                     &
     &              laq%fsh(:,1,ix,iy  )*dalpdy(6)+                     &
     &              laq%fsv(:,1,ix-1,iy)*dalpdy(7)+                     &
     &              laq%fsv(:,1,ix  ,iy)*dalpdy(8)+                     &
     &              laq%fsi(:,1,ix,iy)*dalpdy(9)                        
            ENDDO 
          ENDDO 
        ENDIF 
!-----------------------------------------------------------------------
!     other polynomials:                                                
!-----------------------------------------------------------------------
      CASE DEFAULT 
        j=2*pd+2 
        k=4*pd 
        IF (dmode==0) THEN 
          DO iy=1,my 
            DO ix=1,mx 
              DO i=1,laq%nqty 
                f(i,ix,iy)=                                             &
     &               laq%fs(i,ix-1,iy-1)*alpha(1)+                      &
     &               laq%fs(i,ix  ,iy-1)*alpha(2)+                      &
     &               laq%fs(i,ix-1,iy  )*alpha(3)+                      &
     &               laq%fs(i,ix  ,iy  )*alpha(4)+                      &
     &           SUM(laq%fsh(i,:,ix,iy-1)*alpha(5:j-1:2)+               &
     &               laq%fsh(i,:,ix,iy  )*alpha(6:j  :2)+               &
     &               laq%fsv(i,:,ix-1,iy)*alpha(j+1:k-1:2)+             &
     &               laq%fsv(i,:,ix  ,iy)*alpha(j+2:k:2))+              &
     &           SUM(laq%fsi(i,:,ix,iy)*alpha(k+1:))                    
              ENDDO 
            ENDDO 
          ENDDO 
        ELSE 
          DO iy=1,my 
            DO ix=1,mx 
              DO i=1,laq%nqty 
                f(i,ix,iy)=                                             &
     &               laq%fs(i,ix-1,iy-1)*alpha(1)+                      &
     &               laq%fs(i,ix  ,iy-1)*alpha(2)+                      &
     &               laq%fs(i,ix-1,iy  )*alpha(3)+                      &
     &               laq%fs(i,ix  ,iy  )*alpha(4)+                      &
     &           SUM(laq%fsh(i,:,ix,iy-1)*alpha(5:j-1:2)+               &
     &               laq%fsh(i,:,ix,iy  )*alpha(6:j  :2)+               &
     &               laq%fsv(i,:,ix-1,iy)*alpha(j+1:k-1:2)+             &
     &               laq%fsv(i,:,ix  ,iy)*alpha(j+2:k:2))+              &
     &           SUM(laq%fsi(i,:,ix,iy)*alpha(k+1:))                    
                fx(i,ix,iy)=                                            &
     &                laq%fs(i,ix-1,iy-1)*dalpdx(1)+                    &
     &                laq%fs(i,ix  ,iy-1)*dalpdx(2)+                    &
     &                laq%fs(i,ix-1,iy  )*dalpdx(3)+                    &
     &                laq%fs(i,ix  ,iy  )*dalpdx(4)+                    &
     &            SUM(laq%fsh(i,:,ix,iy-1)*dalpdx(5:j-1:2)+             &
     &                laq%fsh(i,:,ix,iy  )*dalpdx(6:j  :2)+             &
     &                laq%fsv(i,:,ix-1,iy)*dalpdx(j+1:k-1:2)+           &
     &                laq%fsv(i,:,ix  ,iy)*dalpdx(j+2:k:2))+            &
     &            SUM(laq%fsi(i,:,ix,iy)*dalpdx(k+1:))                  
                fy(i,ix,iy)=                                            &
     &                laq%fs(i,ix-1,iy-1)*dalpdy(1)+                    &
     &                laq%fs(i,ix  ,iy-1)*dalpdy(2)+                    &
     &                laq%fs(i,ix-1,iy  )*dalpdy(3)+                    &
     &                laq%fs(i,ix  ,iy  )*dalpdy(4)+                    &
     &            SUM(laq%fsh(i,:,ix,iy-1)*dalpdy(5:j-1:2)+             &
     &                laq%fsh(i,:,ix,iy  )*dalpdy(6:j  :2)+             &
     &                laq%fsv(i,:,ix-1,iy)*dalpdy(j+1:k-1:2)+           &
     &                laq%fsv(i,:,ix  ,iy)*dalpdy(j+2:k:2))+            &
     &            SUM(laq%fsi(i,:,ix,iy)*dalpdy(k+1:))                  
              ENDDO 
            ENDDO 
          ENDDO 
        ENDIF 
      END SELECT 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_2D_all_eval 
!-----------------------------------------------------------------------
!     subprogram 20. lagr_quad_2D_assign_rsc.                           
!     assign a real scalar value to a real lagrange quad structure.     
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_assign_rsc(laq,rscalar) 
                                                                        
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq 
      REAL(r8), INTENT(IN) :: rscalar 
                                                                        
      laq%fs=rscalar 
      IF (ALLOCATED(laq%fsh)) THEN 
        laq%fsh=rscalar 
        laq%fsv=rscalar 
        laq%fsi=rscalar 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_2D_assign_rsc 
!-----------------------------------------------------------------------
!     subprogram 21. lagr_quad_2D_assign_csc.                           
!     assign a real scalar value to a real lagrange quad                
!     structure.                                                        
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_assign_csc(laq,cscalar) 
                                                                        
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq 
      COMPLEX(r8), INTENT(IN) :: cscalar 
                                                                        
      laq%fs=cscalar 
      IF (ALLOCATED(laq%fsh)) THEN 
        laq%fsh=cscalar 
        laq%fsv=cscalar 
        laq%fsi=cscalar 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_2D_assign_csc 
!-----------------------------------------------------------------------
!     subprogram 22. lagr_quad_2D_assign_laq.                           
!     set one real lagrange quad structure equal to another.            
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_assign_laq(laq1,laq2) 
                                                                        
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq1 
      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq2 
                                                                        
!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.                   
!-----------------------------------------------------------------------
      laq1%fs=laq2%fs 
      IF (ALLOCATED(laq1%fsh)) THEN 
        laq1%fsh=laq2%fsh 
        laq1%fsv=laq2%fsv 
        laq1%fsi=laq2%fsi 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_2D_assign_laq 
!-----------------------------------------------------------------------
!     subprogram 23. lagr_quad_2D_assign_int.                           
!     assign a integer value to a real lagrange quad structure.         
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_assign_int(laq,int) 
                                                                        
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq 
      INTEGER(i4), INTENT(IN) :: int 
                                                                        
      laq%fs=int 
      IF (ALLOCATED(laq%fsh)) THEN 
        laq%fsh=int 
        laq%fsv=int 
        laq%fsi=int 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_2D_assign_int 
!-----------------------------------------------------------------------
!     subprogram 26. lagr_quad_2D_basis_assign_loc                      
!     assign data into coefficient arrays for one basis function.       
!                                                                       
!     this is a local version of lagr_quad_basis_assign_arr, where the  
!     the data is located at given poloidal indices, only.              
!-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_basis_assign_loc(laq,data,ibasis,ix,iy) 
                                                                        
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq 
      REAL(r8), DIMENSION(:), INTENT(IN) :: data 
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy 
                                                                        
      INTEGER(i4) :: poly_degree 
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to     
!     the value of ibasis.                                              
!     grid vertex-centered data is first.                               
!-----------------------------------------------------------------------
      poly_degree=laq%n_side+1 
      IF (ibasis==1) THEN 
        laq%fs(:,ix,iy)=data 
        RETURN 
!-----------------------------------------------------------------------
!     horizontal sides.                                                 
!-----------------------------------------------------------------------
      ELSE IF (ibasis<=poly_degree) THEN 
        laq%fsh(:,ibasis-1,ix,iy)=data 
        RETURN 
!-----------------------------------------------------------------------
!     vertical sides.                                                   
!-----------------------------------------------------------------------
      ELSE IF (ibasis<2*poly_degree) THEN 
        laq%fsv(:,ibasis-poly_degree,ix,iy)=data 
        RETURN 
!-----------------------------------------------------------------------
!     interior bases.                                                   
!-----------------------------------------------------------------------
      ELSE 
        laq%fsi(:,ibasis-2*poly_degree+1,ix,iy)=data 
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE lagr_quad_2D_basis_assign_loc 
!-----------------------------------------------------------------------
!     subprogram 31. dump_read_lagr_quad.                               
!     real and imaginary parts are read separately for possible         
!     data conversion.                                                  
!-----------------------------------------------------------------------
      SUBROUTINE dump_read_lagr_quad(laq,lx,ly,name,title) 
      USE io 
                                                                        
      INTEGER(i4), INTENT(IN) :: lx,ly 
      CHARACTER(*), INTENT(IN) :: name 
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title 
      TYPE(lagr_quad_type), INTENT(OUT) :: laq 
                                                
      INTEGER(i4) :: nqtmp,nftmp                        
      REAL(r8) :: iread 
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: rread 
      REAL(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: rread2 
                                                                        
      READ(rstrt_unit) iread 
      nqtmp=NINT(iread)
      READ(rstrt_unit) iread 
      nftmp=NINT(iread)
      READ(rstrt_unit) iread 
      laq%n_side=NINT(iread) 
      READ(rstrt_unit) iread 
      laq%n_int=NINT(iread) 
      CALL lagr_quad_alloc(laq,lx,ly,nqtmp,nftmp,laq%n_side+1_i4)
      laq%name=name 
      IF (SIZE(title)<SIZE(laq%title)) THEN 
        laq%title=title(1) 
      ELSE 
        laq%title=title 
      ENDIF 
                                                                        
      ALLOCATE(rread(laq%nqty,0:lx,0:ly,laq%nfour)) 
      READ(rstrt_unit) rread 
      laq%fs=rread 
      READ(rstrt_unit) rread 
      laq%fs=laq%fs+(0,1)*rread 
      DEALLOCATE(rread) 
                                                                        
      IF (laq%n_side>0) THEN 
        ALLOCATE(rread2(laq%nqty,laq%n_side,1:lx,0:ly,laq%nfour)) 
        READ(rstrt_unit) rread2 
        laq%fsh=rread2 
        READ(rstrt_unit) rread2 
        laq%fsh=laq%fsh+(0,1)*rread2 
        DEALLOCATE(rread2) 
        ALLOCATE(rread2(laq%nqty,laq%n_side,0:lx,1:ly,laq%nfour)) 
        READ(rstrt_unit) rread2 
        laq%fsv=rread2 
        READ(rstrt_unit) rread2 
        laq%fsv=laq%fsv+(0,1)*rread2 
        DEALLOCATE(rread2) 
        ALLOCATE(rread2(laq%nqty,laq%n_int,1:lx,1:ly,laq%nfour)) 
        READ(rstrt_unit) rread2 
        laq%fsi=rread2 
        READ(rstrt_unit) rread2 
        laq%fsi=laq%fsi+(0,1)*rread2 
        DEALLOCATE(rread2) 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE dump_read_lagr_quad 
!-----------------------------------------------------------------------
!     subprogram 32. dump_read_lagr_quad_2D.                            
!     2D version.                                                       
!-----------------------------------------------------------------------
      SUBROUTINE dump_read_lagr_quad_2D(laq,lx,ly,name,title) 
      USE io 
                                                                        
      INTEGER(i4), INTENT(IN) :: lx,ly 
      CHARACTER(*), INTENT(IN) :: name 
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title 
      TYPE(lagr_quad_2D_type), INTENT(OUT) :: laq 
        
      INTEGER(i4) :: nqtmp
      REAL(r8) :: iread 
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rread 
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: rread2 
                                                                        
      READ(rstrt_unit) iread
      nqtmp=NINT(iread)
      READ(rstrt_unit) iread 
      laq%n_side=NINT(iread) 
      READ(rstrt_unit) iread 
      laq%n_int=NINT(iread) 
      CALL lagr_quad_alloc(laq,lx,ly,nqtmp,laq%n_side+1_i4) 
      laq%name=name 
      IF (SIZE(title)<SIZE(laq%title)) THEN 
        laq%title=title(1) 
      ELSE 
        laq%title=title 
      ENDIF 
                                                                        
      ALLOCATE(rread(laq%nqty,0:lx,0:ly)) 
      READ(rstrt_unit) rread 
      laq%fs=rread 
      DEALLOCATE(rread) 
                                                                        
      IF (laq%n_side>0) THEN 
        ALLOCATE(rread2(laq%nqty,laq%n_side,1:lx,0:ly)) 
        READ(rstrt_unit) rread2 
        laq%fsh=rread2 
        DEALLOCATE(rread2) 
        ALLOCATE(rread2(laq%nqty,laq%n_side,0:lx,1:ly)) 
        READ(rstrt_unit) rread2 
        laq%fsv=rread2 
        DEALLOCATE(rread2) 
        ALLOCATE(rread2(laq%nqty,laq%n_int,1:lx,1:ly)) 
        READ(rstrt_unit) rread2 
        laq%fsi=rread2 
        DEALLOCATE(rread2) 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE dump_read_lagr_quad_2D 
#ifdef HAVE_FC_HDF5
!-----------------------------------------------------------------------
!     subprogram 33. h5_dump_lagr_quad.                              
!     real and imaginary parts are witten separately for possible       
!     data conversion.                                                  
!-----------------------------------------------------------------------
      SUBROUTINE h5_dump_lagr_quad(laq,fname,gid,blid)
      USE io 
                                                                        
      TYPE(lagr_quad_type), INTENT(IN) :: laq 
      CHARACTER(*), INTENT(IN) :: fname,blid
      INTEGER(HID_T), INTENT(IN) :: gid

      INTEGER(i4) :: dn,nx,ny,iq,ifr,ic,is,ins,ine,ii
      REAL(r8), ALLOCATABLE :: tarr(:,:,:),ttarr(:,:,:)
      CHARACTER(64) :: mdname

      h5in%vsMD=TRIM(fname)
      mdname = TRIM(h5in%vsMD)

      dn=laq%n_side+1
      nx=laq%mx*dn
      ny=laq%my*dn
      ALLOCATE(tarr(0:nx,0:ny,laq%nqty*laq%nfour),                      &
     &         ttarr(laq%nqty*laq%nfour,0:nx,0:ny))
      DO iq=1,laq%nqty
        DO ifr=1,laq%nfour
          ic=iq+(ifr-1)*laq%nqty
          tarr(0:nx:dn,0:ny:dn,ic)=REAL(laq%fs(iq,:,:,ifr))
          IF (ALLOCATED(laq%fsh)) THEN
            DO is=1,laq%n_side
              tarr(is:nx:dn,0:ny:dn,ic)=REAL(laq%fsh(iq,is,:,:,ifr))
              tarr(0:nx:dn,is:ny:dn,ic)=REAL(laq%fsv(iq,is,:,:,ifr))
              ins=(is-1)*laq%n_side+1; ine=is*laq%n_side
              DO ii=ins,ine
                tarr(ii-ins+1:nx:dn,is:ny:dn,ic)=                       &
     &            REAL(laq%fsi(iq,ii,:,:,ifr))
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      h5in%vsMD="Re_"//mdname
      DO iq=1,laq%nqty*laq%nfour
        ttarr(iq,0:nx,0:ny)=tarr(0:nx,0:ny,iq)
      ENDDO
      CALL dump_h5(gid,"re"//TRIM(fname)//TRIM(blid),ttarr,h5in,h5err)
      DO iq=1,laq%nqty
        DO ifr=1,laq%nfour
          ic=iq+(ifr-1)*laq%nqty
          tarr(0:nx:dn,0:ny:dn,ic)=AIMAG(laq%fs(iq,:,:,ifr))
          IF (ALLOCATED(laq%fsh)) THEN
            DO is=1,laq%n_side
              tarr(is:nx:dn,0:ny:dn,ic)=AIMAG(laq%fsh(iq,is,:,:,ifr))
              tarr(0:nx:dn,is:ny:dn,ic)=AIMAG(laq%fsv(iq,is,:,:,ifr))
              ins=(is-1)*laq%n_side+1; ine=is*laq%n_side
              DO ii=ins,ine
                tarr(ii-ins+1:nx:dn,is:ny:dn,ic)=                       &
     &            AIMAG(laq%fsi(iq,ii,:,:,ifr))
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      h5in%vsMD="Im_"//mdname
      DO iq=1,laq%nqty*laq%nfour
        ttarr(iq,0:nx,0:ny)=tarr(0:nx,0:ny,iq)
      ENDDO
      CALL dump_h5(gid,"im"//TRIM(fname)//TRIM(blid),ttarr,h5in,h5err)
      DEALLOCATE(tarr,ttarr)
      h5in%vsMD=" "
                                                                        
      RETURN 
      END SUBROUTINE h5_dump_lagr_quad 
!-----------------------------------------------------------------------
!     subprogram 34. h5_dump_lagr_quad_2D.                           
!     2D version.                                                       
!-----------------------------------------------------------------------
      SUBROUTINE h5_dump_lagr_quad_2D(laq,fname,gid,blid) 
      USE io 
                                                                        
      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq 
      CHARACTER(*), INTENT(IN) :: fname,blid
      INTEGER(HID_T), INTENT(IN) :: gid

      INTEGER(i4) :: dn,nx,ny,iq,is,ins,ine,ii
      REAL(r8), ALLOCATABLE :: tarr(:,:,:),ttarr(:,:,:)

      h5in%vsMD=TRIM(fname)

      dn=laq%n_side+1
      nx=laq%mx*dn
      ny=laq%my*dn
      ALLOCATE(tarr(0:nx,0:ny,laq%nqty)) 
      ALLOCATE(ttarr(laq%nqty,0:nx,0:ny))
      DO iq=1,laq%nqty
        tarr(0:nx:dn,0:ny:dn,iq)=laq%fs(iq,:,:)
        IF (ALLOCATED(laq%fsh)) THEN
          DO is=1,laq%n_side
            tarr(is:nx:dn,0:ny:dn,iq)=laq%fsh(iq,is,:,:)
            tarr(0:nx:dn,is:ny:dn,iq)=laq%fsv(iq,is,:,:)
            ins=(is-1)*laq%n_side+1; ine=is*laq%n_side
            DO ii=ins,ine
              tarr(ii-ins+1:nx:dn,is:ny:dn,iq)=laq%fsi(iq,ii,:,:)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
!     Reorder indices as compMajorF crashes visit
!     (really this reordering is compMinorF, but visit displays 
!     correctly with the default compMinorC).
!-----------------------------------------------------------------------
      DO iq=1,laq%nqty
        ttarr(iq,0:nx,0:ny)=tarr(0:nx,0:ny,iq) 
      ENDDO
      CALL dump_h5(gid,TRIM(fname)//TRIM(blid),ttarr,h5in,h5err)
      DEALLOCATE(tarr,ttarr)
      h5in%vsMD=" "
                                                                        
      RETURN 
      END SUBROUTINE h5_dump_lagr_quad_2D 
!-----------------------------------------------------------------------
!     subprogram 35. h5_read_lagr_quad.                               
!     real and imaginary parts are read separately for possible         
!     data conversion.                                                  
!-----------------------------------------------------------------------
      SUBROUTINE h5_read_lagr_quad(laq,lx,ly,pd,nfour,nqty2,            &
     &                             fname,gid,blid)
      USE io 
                                                                        
      INTEGER(i4), INTENT(IN) :: lx,ly 
      TYPE(lagr_quad_type), INTENT(OUT) :: laq 
      CHARACTER(*), INTENT(IN) :: fname,blid
      INTEGER(i4), INTENT(IN) :: pd,nfour,nqty2
      INTEGER(HID_T), INTENT(IN) :: gid
                                                                        
      INTEGER(i4) :: dn,nx,ny,iq,ifr,ic,is,ins,ine,ii,nqty
      REAL(r8), ALLOCATABLE :: tarr(:,:,:),itarr(:,:,:),ttarr(:,:,:)
      INTEGER(HID_T) :: dsetid
      INTEGER :: error
      INTEGER(HSIZE_T), DIMENSION(3) :: dims

      IF (.NOT.obj_exists(gid,"re"//TRIM(fname)//TRIM(blid),h5err))     &
        THEN ! Allocate, set to zero and return
        CALL lagr_quad_alloc(laq,lx,ly,nqty2,nfour,pd,TRIM(fname))  
        laq%fs=0._r8
        IF (ALLOCATED(laq%fsh)) laq%fsh=0._r8
        IF (ALLOCATED(laq%fsv)) laq%fsv=0._r8
        IF (ALLOCATED(laq%fsi)) laq%fsi=0._r8
        RETURN
      ENDIF

      CALL h5dopen_f(gid,"re"//TRIM(fname)//TRIM(blid),dsetid,error)
      CALL read_dims(dsetid,dims,h5err)
      CALL h5dclose_f(dsetid,error)
      nqty=dims(1)/nfour
      IF (nqty/=nqty2)                                                  &
     &  CALL nim_stop('unexpected lagr_quad nqty in '//fname)
      IF (dims(2)/=lx*pd+1.OR.dims(3)/=ly*pd+1)                         &
     &  CALL nim_stop('unexpected lagr_quad size in '//fname)

      CALL lagr_quad_alloc(laq,lx,ly,nqty,nfour,pd,TRIM(fname))
      ii=len(fname)
      IF (ii>=6) THEN
        laq%name=TRIM(fname(1:6))
      ELSE
        laq%name=TRIM(fname)
      ENDIF
 
      dn=laq%n_side+1
      nx=laq%mx*dn
      ny=laq%my*dn
      ALLOCATE( tarr(0:nx,0:ny,laq%nqty*laq%nfour),                     &
     &         itarr(0:nx,0:ny,laq%nqty*laq%nfour),                     &
     &         ttarr(laq%nqty*laq%nfour,0:nx,0:ny)) 
      CALL read_h5(gid,"re"//TRIM(fname)//TRIM(blid),ttarr,h5in,h5err)
      DO iq=1,laq%nqty*laq%nfour
        tarr(0:nx,0:ny,iq)=ttarr(iq,0:nx,0:ny)
      ENDDO
      CALL read_h5(gid,"im"//TRIM(fname)//TRIM(blid),ttarr,h5in,h5err)
      DO iq=1,laq%nqty*laq%nfour
        itarr(0:nx,0:ny,iq)=ttarr(iq,0:nx,0:ny)
      ENDDO
      DO iq=1,laq%nqty
        DO ifr=1,laq%nfour
          ic=iq+(ifr-1)*laq%nqty
          laq%fs(iq,:,:,ifr)=tarr(0:nx:dn,0:ny:dn,ic)                   &
     &                      +(0,1)*itarr(0:nx:dn,0:ny:dn,ic)
          IF (ALLOCATED(laq%fsh)) THEN
            DO is=1,laq%n_side
              laq%fsh(iq,is,:,:,ifr)=tarr(is:nx:dn,0:ny:dn,ic)          &
     &                              +(0,1)*itarr(is:nx:dn,0:ny:dn,ic)
              laq%fsv(iq,is,:,:,ifr)=tarr(0:nx:dn,is:ny:dn,ic)          &
     &                              +(0,1)*itarr(0:nx:dn,is:ny:dn,ic)
              ins=(is-1)*laq%n_side+1; ine=is*laq%n_side
              DO ii=ins,ine
                laq%fsi(iq,ii,:,:,ifr)=tarr(ii-ins+1:nx:dn,is:ny:dn,ic) &
     &                          +(0,1)*itarr(ii-ins+1:nx:dn,is:ny:dn,ic)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE(tarr,itarr,ttarr)
                                                                        
      RETURN 
      END SUBROUTINE h5_read_lagr_quad 
!-----------------------------------------------------------------------
!     subprogram 36. h5_read_lagr_quad_2D.                            
!     2D version.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE h5_read_lagr_quad_2D(laq,lx,ly,pd,nqty2,fname,gid,blid)
      USE io 
                                                                        
      INTEGER(i4), INTENT(IN) :: lx,ly 
      TYPE(lagr_quad_2D_type), INTENT(OUT) :: laq 
      CHARACTER(*), INTENT(IN) :: fname,blid
      INTEGER(HID_T), INTENT(IN) :: gid
      INTEGER(i4), INTENT(IN) :: pd,nqty2
                                                                        
      INTEGER(i4) :: dn,nx,ny,iq,is,ins,ine,ii,nqty
      REAL(r8), ALLOCATABLE :: tarr(:,:,:),ttarr(:,:,:)
      INTEGER(HID_T) :: dsetid
      INTEGER :: error
      INTEGER(HSIZE_T), DIMENSION(3) :: dims

      IF (.NOT.obj_exists(gid,TRIM(fname)//TRIM(blid),h5err))           &
        THEN ! Allocate, set to zero and return
        CALL lagr_quad_alloc(laq,lx,ly,nqty2,pd,TRIM(fname))  
        laq%fs=0._r8
        IF (ALLOCATED(laq%fsh)) laq%fsh=0._r8
        IF (ALLOCATED(laq%fsv)) laq%fsv=0._r8
        IF (ALLOCATED(laq%fsi)) laq%fsi=0._r8
        RETURN
      ENDIF

      CALL h5dopen_f(gid,TRIM(fname)//TRIM(blid),dsetid,error)
      CALL read_dims(dsetid,dims,h5err)
      CALL h5dclose_f(dsetid,error)
      nqty=dims(1)
      IF (nqty/=nqty2)                                                  &
     &  CALL nim_stop('unexpected 2D lagr_quad nqty in '//fname)
      IF (dims(2)/=lx*pd+1.OR.dims(3)/=ly*pd+1)                         &
     &  CALL nim_stop('unexpected 2D lagr_quad size in '//fname)

      CALL lagr_quad_alloc(laq,lx,ly,nqty,pd,TRIM(fname))
      ii=len(fname)
      IF (ii>=6) THEN
        laq%name=TRIM(fname(1:6))
      ELSE
        laq%name=TRIM(fname)
      ENDIF
       
      dn=laq%n_side+1
      nx=laq%mx*dn
      ny=laq%my*dn
      ALLOCATE(tarr(0:nx,0:ny,laq%nqty),ttarr(laq%nqty,0:nx,0:ny)) 
      CALL read_h5(gid,TRIM(fname)//TRIM(blid),ttarr,h5in,h5err)
      DO iq=1,laq%nqty
        tarr(0:nx,0:ny,iq)=ttarr(iq,0:nx,0:ny)
      ENDDO
      DO iq=1,laq%nqty
        laq%fs(iq,:,:)=tarr(0:nx:dn,0:ny:dn,iq)
        IF (ALLOCATED(laq%fsh)) THEN
          DO is=1,laq%n_side
            laq%fsh(iq,is,:,:)=tarr(is:nx:dn,0:ny:dn,iq)
            laq%fsv(iq,is,:,:)=tarr(0:nx:dn,is:ny:dn,iq)
            ins=(is-1)*laq%n_side+1; ine=is*laq%n_side
            DO ii=ins,ine
              laq%fsi(iq,ii,:,:)=tarr(ii-ins+1:nx:dn,is:ny:dn,iq)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      DEALLOCATE(tarr,ttarr)
                                                                        
      RETURN 
      END SUBROUTINE h5_read_lagr_quad_2D 
#endif /* HAVE_FC_HDF5 */

!-----------------------------------------------------------------------
!     close module.                                                     
!-----------------------------------------------------------------------
      END MODULE lagr_quad_mod
