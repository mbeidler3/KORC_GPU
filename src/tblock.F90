!-----------------------------------------------------------------------
!     $Id: tblock.F90 4620 2015-12-03 19:18:40Z jking $
!     contains the module for handling finite element computations      
!     on unstructured triangular grid-blocks.                           
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.                                                
!-----------------------------------------------------------------------
!     1.  tblock_set                                                    
!     2.  tblock_make_real_matrix.                                      
!     3.  tblock_make_comp_matrix.                                      
!     4.  tblock_get_real_rhs.                                          
!     5.  tblock_get_comp_rhs.                                          
!     6.  tblock_get_comp_rhs_q.                                        
!     7.  tblock_basis_set.                                             
!     8.  tblock_basis_dealloc.                                         
!     9.  tblock_real_qp_update.                                        
!     10. tblock_comp_qp_update.                                        
!     11. tblock_real_qp_alloc.                                         
!     12. tblock_comp_qp_alloc.                                         
!     13. tblock_real_qp_dealloc.                                       
!     14. tblock_comp_qp_dealloc.                                       
!-----------------------------------------------------------------------
!     subprogram 0. tblock.                                             
!     module definitions.                                               
!-----------------------------------------------------------------------
      MODULE tblock 
      USE local 
      USE tblock_type_mod 
      IMPLICIT NONE 
                                                                        
      INTERFACE tblock_make_matrix 
        MODULE PROCEDURE tblock_make_real_matrix,tblock_make_comp_matrix 
      END INTERFACE 
                                                                        
      INTERFACE tblock_get_rhs 
        MODULE PROCEDURE tblock_get_real_rhs,tblock_get_comp_rhs 
      END INTERFACE 
                                                                        
      INTERFACE tblock_qp_update 
        MODULE PROCEDURE tblock_real_qp_update,tblock_comp_qp_update 
      END INTERFACE 
                                                                        
      INTERFACE tblock_qp_alloc 
        MODULE PROCEDURE tblock_real_qp_alloc,tblock_comp_qp_alloc 
      END INTERFACE 
                                                                        
      INTERFACE tblock_qp_dealloc 
        MODULE PROCEDURE tblock_real_qp_dealloc,tblock_comp_qp_dealloc 
      END INTERFACE 
                                                                        
!-----------------------------------------------------------------------
!-PRE declarations for Gaussian quadrature weights.                     
!-----------------------------------------------------------------------
      INTEGER(i4), PARAMETER, PRIVATE :: ngauss=7 
      REAL(r8), PARAMETER, PRIVATE :: wq0=0.225_r8,                     &
     &     wq1=0.1259391805448271_r8,wq2=0.1323941527885062_r8          
      REAL(r8), DIMENSION(ngauss), PARAMETER, PRIVATE ::                &
     &     quad_weight=(/wq0,wq1,wq1,wq1,wq2,wq2,wq2/)                  
                                                                        
      CONTAINS 
!-----------------------------------------------------------------------
!     subprogram 1. tblock_set.                                         
!     set the weights for quadratures.                                  
!-----------------------------------------------------------------------
      SUBROUTINE tblock_set(tb) 
                                                                        
      TYPE(tblock_type), INTENT(INOUT) :: tb 
!-----------------------------------------------------------------------
!-PRE                                                                   
!     set number of quadrature points and weights according to input.   
!     the number of points is now adjusted automatically with           
!     poly_degree.                                                      
!-----------------------------------------------------------------------
      tb%ng=ngauss 
      ALLOCATE(tb%wg(tb%ng)) 
      tb%wg=quad_weight 
!-----------------------------------------------------------------------
!     terminate tblock_set.                                             
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tblock_set 
!-----------------------------------------------------------------------
!     subprogram 2. tblock_make_real_matrix.                            
!     computes a linear response matrix for a supplied integrand        
!     subprogram that fits the interface block at the beginning         
!     of this subprogram.                                               
!-----------------------------------------------------------------------
      SUBROUTINE tblock_make_real_matrix(tb,mat,integrand,nqty)
      USE rblock_type_mod 
      USE matrix_type_mod 
                                                                        
      TYPE(tblock_type), INTENT(INOUT) :: tb 
      TYPE(matrix_element_type3), DIMENSION(0:), INTENT(OUT) :: mat 
      INTEGER(i4), INTENT(IN) :: nqty 
                                                                        
      TYPE(rblock_type) :: rdum 
                                                                        
      REAL(r8), DIMENSION(nqty,nqty,tb%mcell,3,3) :: integr
      REAL(r8), DIMENSION(:,:), POINTER :: bigr 
                                                                        
      INTEGER(i4) :: icell,iv,jv,iqty,jqty,inbr,nnbr,node,iv1,jv1,      &
     &               ierror
!-----------------------------------------------------------------------
!     interface block for the integrand computation.                    
!-----------------------------------------------------------------------
#include "integrand_real.finc"
!-----------------------------------------------------------------------
!     zero matrix.  quadrature-point loop is now within integrands.     
!-----------------------------------------------------------------------
      DO iv=0,tb%mvert 
        mat(iv)%element=0._r8 
      ENDDO 
!-----------------------------------------------------------------------
!     evaluate the integrand at all quadrature point.                   
!-----------------------------------------------------------------------
      bigr=>tb%tgeom%bigr 
      CALL integrand(integr,bigr,rdum,tb,1_i4)
!-----------------------------------------------------------------------
!     quadrature-point looping is now done inside the integrand routine,
!     and integrand is summed for each basis function, element by       
!     element.                                                          
!                                                                       
!     factors of Jacobian and quadrature weight are already in the basis
!     function arrays.                                                  
!-----------------------------------------------------------------------
      DO icell=1,tb%mcell 
        DO iv=1,3 
          iv1=tb%tgeom%vertex(icell,iv) 
          DO jv=1,3 
            jv1=tb%tgeom%vertex(icell,jv) 
            DO inbr=0,SIZE(tb%tgeom%neighbor(iv1)%vertex)-1 
              IF (tb%tgeom%neighbor(iv1)%vertex(inbr)==jv1) THEN 
                mat(iv1)%element(1:nqty,1:nqty,inbr)                    &
     &            =mat(iv1)%element(1:nqty,1:nqty,inbr)                 &
     &            +integr(:,:,icell,jv,iv)
              ENDIF 
            ENDDO 
          ENDDO 
        ENDDO 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tblock_make_real_matrix 
!-----------------------------------------------------------------------
!     subprogram 3. tblock_make_comp_matrix.                            
!     computes a linear response matrix for a supplied integrand        
!     subprogram that fits the interface block at the beginning         
!     of this subprogram.                                               
!                                                                       
!     complex version                                                   
!-----------------------------------------------------------------------
      SUBROUTINE tblock_make_comp_matrix(tb,mat,integrand,nqty)
      USE rblock_type_mod 
      USE matrix_type_mod 
      USE time 
                                                                        
      TYPE(tblock_type), INTENT(INOUT) :: tb 
      TYPE(comp_matrix_element_type3), DIMENSION(0:), INTENT(OUT) :: mat 
      INTEGER(i4), INTENT(IN) :: nqty 
                                                                        
      TYPE(rblock_type) :: rdum 
                                                                        
      COMPLEX(r8), DIMENSION(nqty,nqty,tb%mcell,3,3) :: integr
      REAL(r8), DIMENSION(:,:), POINTER :: bigr 
                                                                        
      INTEGER(i4) :: icell,iv,jv,iqty,jqty,inbr,nnbr,node,iv1,jv1,ierror
!-----------------------------------------------------------------------
!     interface block for the integrand computation.                    
!-----------------------------------------------------------------------
#include "integrand_comp.finc"
!-----------------------------------------------------------------------
!     zero matrix.  quadrature-point loop is now within integrands.     
!-----------------------------------------------------------------------
      DO iv=0,tb%mvert 
        mat(iv)%element=0._r8 
      ENDDO 
!-----------------------------------------------------------------------
!     evaluate the integrand at all quadrature points.                  
!-----------------------------------------------------------------------
      bigr=>tb%tgeom%bigr 
      CALL integrand(integr,bigr,rdum,tb,1_i4)
!-----------------------------------------------------------------------
!     quadrature-point looping is now done inside the integrand routine,
!     and integrand is summed for each basis function, element by       
!     element.                                                          
!                                                                       
!     factors of Jacobian and quadrature weight are already in the basis
!     function arrays.                                                  
!-----------------------------------------------------------------------
      DO icell=1,tb%mcell 
        DO iv=1,3 
          iv1=tb%tgeom%vertex(icell,iv) 
          DO jv=1,3 
            jv1=tb%tgeom%vertex(icell,jv) 
            DO inbr=0,SIZE(tb%tgeom%neighbor(iv1)%vertex)-1 
              IF (tb%tgeom%neighbor(iv1)%vertex(inbr)==jv1) THEN 
                mat(iv1)%element(1:nqty,1:nqty,inbr)                    &
     &            =mat(iv1)%element(1:nqty,1:nqty,inbr)                 &
     &            +integr(:,:,icell,jv,iv)
              ENDIF 
            ENDDO 
          ENDDO 
        ENDDO 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tblock_make_comp_matrix 
!-----------------------------------------------------------------------
!     subprogram 4. tblock_get_real_rhs.                                
!     performs finite-element integrations for the rhs of an equation   
!     where the integrand is computed with a supplied subroutine.       
!-----------------------------------------------------------------------
      SUBROUTINE tblock_get_real_rhs(tb,rhs,integrand,nq)
      USE rblock_type_mod 
      USE vector_type_mod 
                                                                        
      INTEGER(i4), INTENT(IN) :: nq 
      TYPE(tblock_type), INTENT(INOUT) :: tb 
      TYPE(vector_type), INTENT(INOUT) :: rhs 
                                                                        
      TYPE(rblock_type) :: rdum 
                                                                        
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: integr
      REAL(r8), DIMENSION(:,:), POINTER :: bigr 
      REAL(r8), DIMENSION(:,:,:), POINTER :: rhsg 
      REAL(r8), DIMENSION(:,:,:,:), POINTER :: rhss,rhsi 
                                                                        
      INTEGER(i4) :: node,ivert,icell,ierror,iseg 
      INTEGER(i4) :: iq,iv,nv,ib,mcell
      INTEGER(i4) :: start_seg,start_int,n_grid,n_seg,n_int 
!-----------------------------------------------------------------------
!     interface block for the integrand computation.                    
!-----------------------------------------------------------------------
#include "integrand_real_rhs.finc"
!-----------------------------------------------------------------------
!     examine the vector to determine what bases are used.              
!-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN 
        n_grid=3 
        rhsg=>rhs%arr 
        rhsg(1:nq,:,:)=0 
      ELSE 
        n_grid=0 
      ENDIF 
      iv=n_grid+1 
      IF (ASSOCIATED(rhs%arrh)) THEN 
        start_seg=iv 
        n_seg=SIZE(rhs%arrh,4) 
        rhss=>rhs%arrh 
        rhss(1:nq,:,:,:)=0 
      ELSE 
        n_seg=0 
      ENDIF 
      iv=iv+3*n_seg 
      IF (ASSOCIATED(rhs%arri)) THEN 
        start_int=iv 
        n_int=SIZE(rhs%arri,4) 
        rhsi=>rhs%arri 
        rhsi(1:nq,:,:,:)=0 
      ELSE 
        n_int=0 
      ENDIF 
      iv=iv+n_int-1 
      ALLOCATE(integr(nq,tb%mcell,iv))
!-----------------------------------------------------------------------
!     evaluate the integrand at all quadrature points.                  
!-----------------------------------------------------------------------
      bigr=>tb%tgeom%bigr 
      CALL integrand(integr,bigr,rdum,tb,1_i4)
!-----------------------------------------------------------------------
!     accumulate element contributions.                                 
!     grid vertex-centered bases first.                                 
!                                                                       
!     factors of Jacobian and quadrature weight are already in the      
!     test function arrays.                                             
!-----------------------------------------------------------------------
      IF (n_grid==3) THEN 
        DO ib=1,3 
          DO icell=1,tb%mcell 
            ivert=tb%tgeom%vertex(icell,ib) 
            rhsg(1:nq,ivert,0)=rhsg(1:nq,ivert,0)+integr(:,icell,ib)
          ENDDO 
        ENDDO 
      ENDIF 
      DEALLOCATE(integr)
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tblock_get_real_rhs 
!-----------------------------------------------------------------------
!     subprogram 5. tblock_get_comp_rhs.                                
!     performs finite-element integrations for the rhs of an equation   
!     where the integrand is computed with a supplied subroutine.
!     Depreciated until needed.
!-----------------------------------------------------------------------
      SUBROUTINE tblock_get_comp_rhs(tb,rhs,integrand,nq,nfour)
      USE rblock_type_mod 
      USE vector_type_mod 
                                                                        
      INTEGER(i4), INTENT(IN) :: nq,nfour 
      TYPE(tblock_type), INTENT(INOUT) :: tb 
      TYPE(cvector_type), INTENT(INOUT) :: rhs 
                                                                        
      TYPE(rblock_type) :: rdum 
                                                                        
      REAL(r8), DIMENSION(:,:), POINTER :: bigr 
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: integr
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: rhsg 
      COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: rhss,rhsi 
                                                                        
      INTEGER(i4) :: node,ivert,icell,ierror,iseg 
      INTEGER(i4) :: iq,iv,nv,ib,mcell,jf
      INTEGER(i4) :: start_seg,start_int,n_grid,n_seg,n_int 
!-----------------------------------------------------------------------
!     interface block for the integrand computation.                    
!-----------------------------------------------------------------------
#include "integrand_comp_rhs.finc"
!-----------------------------------------------------------------------
!     examine the vector to determine what bases are used.              
!-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN 
        n_grid=3 
        rhsg=>rhs%arr 
        rhsg(1:nq,:,:,1:nfour)=0 
      ELSE 
        n_grid=0 
      ENDIF 
      iv=n_grid+1 
      IF (ASSOCIATED(rhs%arrh)) THEN 
        start_seg=iv 
        n_seg=SIZE(rhs%arrh,4) 
        rhss=>rhs%arrh 
        rhss(1:nq,:,:,:,1:nfour)=0 
      ELSE 
        n_seg=0 
      ENDIF 
      iv=iv+3*n_seg 
      IF (ASSOCIATED(rhs%arri)) THEN 
        start_int=iv 
        n_int=SIZE(rhs%arri,4) 
        rhsi=>rhs%arri 
        rhsi(1:nq,:,:,:,1:nfour)=0 
      ELSE 
        n_int=0 
      ENDIF 
      iv=iv+n_int-1 
      ALLOCATE(integr(nq,tb%mcell,iv,nfour))
!-----------------------------------------------------------------------
!     evaluate the integrand at the quadrature point.                   
!-----------------------------------------------------------------------
      bigr=>tb%tgeom%bigr 
      CALL integrand(integr,bigr,rdum,tb,1_i4)
!-----------------------------------------------------------------------
!     accumulate element contributions.                                 
!     grid vertex-centered bases first.                                 
!                                                                       
!     factors of Jacobian and quadrature weight are already in the      
!     test function arrays.                                             
!-----------------------------------------------------------------------
      IF (n_grid==3) THEN 
        DO jf=1,nfour 
          DO ib=1,3 
            DO icell=1,tb%mcell 
              ivert=tb%tgeom%vertex(icell,ib) 
              rhsg(1:nq,ivert,0,jf)=rhsg(1:nq,ivert,0,jf)               &
     &          +integr(:,icell,ib,jf)
            ENDDO 
          ENDDO 
        ENDDO 
      ENDIF 
      DEALLOCATE(integr)
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tblock_get_comp_rhs 
!-----------------------------------------------------------------------
!     subprogram 6. tblock_get_comp_rhs_q.                              
!     performs finite-element integrations for the rhs of an equation   
!     where the integrand is computed with a supplied subroutine.       
!     same as tblock_get_comp_rhs but the quantity and Fourier component
!     indices are assumed to have the correct dimension for efficiency. 
!-----------------------------------------------------------------------
      SUBROUTINE tblock_get_comp_rhs_q(tb,rhs,integrand)
      USE rblock_type_mod 
      USE vector_type_mod 
                                                                        
      TYPE(tblock_type), INTENT(INOUT) :: tb 
      TYPE(cvector_type), INTENT(INOUT) :: rhs 
                                                                        
      TYPE(rblock_type) :: rdum 
                                                                        
      REAL(r8), DIMENSION(:,:), POINTER :: bigr 
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: integr
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: rhsg 
      COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: rhss,rhsi 
                                                                        
      INTEGER(i4) :: node,ivert,icell,ierror,iseg 
      INTEGER(i4) :: iq,iv,nv,ib,mcell,jf 
      INTEGER(i4) :: start_seg,start_int,n_grid,n_seg,n_int,nq,nfour 
!-----------------------------------------------------------------------
!     interface block for the integrand computation.                    
!-----------------------------------------------------------------------
#include "integrand_comp_rhs.finc"
!-----------------------------------------------------------------------
!     examine the vector to determine what bases are used.              
!-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN 
        n_grid=3 
        rhsg=>rhs%arr 
        nq=SIZE(rhsg,1) 
        nfour=SIZE(rhsg,4) 
        rhsg=0 
      ELSE 
        n_grid=0 
      ENDIF 
      iv=n_grid+1 
      IF (ASSOCIATED(rhs%arrh)) THEN 
        start_seg=iv 
        n_seg=SIZE(rhs%arrh,4) 
        rhss=>rhs%arrh 
        nq=SIZE(rhss,1) 
        nfour=SIZE(rhss,5) 
        rhss=0 
      ELSE 
        n_seg=0 
      ENDIF 
      iv=iv+3*n_seg 
      IF (ASSOCIATED(rhs%arri)) THEN 
        start_int=iv 
        n_int=SIZE(rhs%arri,4) 
        rhsi=>rhs%arri 
        nq=SIZE(rhsi,1) 
        nfour=SIZE(rhsi,5) 
        rhsi=0 
      ELSE 
        n_int=0 
      ENDIF 
      iv=iv+n_int-1 
      ALLOCATE(integr(nq,tb%mcell,iv,nfour))
!-----------------------------------------------------------------------
!     evaluate the integrand at all quadrature points.                  
!-----------------------------------------------------------------------
      bigr=>tb%tgeom%bigr 
      CALL integrand(integr,bigr,rdum,tb,1_i4)
!-----------------------------------------------------------------------
!     accumulate element contributions.                                 
!     grid vertex-centered bases first.                                 
!                                                                       
!     factors of Jacobian and quadrature weight are already in the      
!     test function arrays.                                             
!-----------------------------------------------------------------------
      IF (n_grid==3) THEN 
        DO jf=1,nfour 
          DO ib=1,3 
            DO icell=1,tb%mcell 
              ivert=tb%tgeom%vertex(icell,ib) 
              rhsg(:,ivert,0,jf)=rhsg(:,ivert,0,jf)                     &
     &          +integr(:,icell,ib,jf)
            ENDDO 
          ENDDO 
        ENDDO 
      ENDIF 
      DEALLOCATE(integr)
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tblock_get_comp_rhs_q 
!-----------------------------------------------------------------------
!     subprogram 7. tblock_basis_set.                                   
!     set the locations and weights for quadratures.                    
!-----------------------------------------------------------------------
      SUBROUTINE tblock_basis_set(tb,geom) 
                                                                        
      TYPE(tblock_type), INTENT(INOUT) :: tb 
      CHARACTER(*), INTENT(IN) :: geom 
                                                                        
      INTEGER(i4) :: node,iv,icell 
!-----------------------------------------------------------------------
!     copy the basis function values into arrays with dimensions to     
!     match those in structured blocks (except wjac, which is the       
!     quadrature weight times the Jacobian).                            
!                                                                       
!     for toroidal geometry dalpdrc holds d(alpha)/dr + alpha/r, and    
!     for linear geometry it's just d(alpha)/dx.                        
!-----------------------------------------------------------------------
      ALLOCATE(tb%tgeom%alpha_arr(tb%ng,tb%tgeom%mcell,3)) 
      ALLOCATE(tb%tgeom%dalpdr(tb%ng,tb%tgeom%mcell,3)) 
      ALLOCATE(tb%tgeom%dalpdz(tb%ng,tb%tgeom%mcell,3)) 
      ALLOCATE(tb%tgeom%dalpdrc(tb%ng,tb%tgeom%mcell,3)) 
      ALLOCATE(tb%tgeom%alpham_arr(tb%ng,tb%tgeom%mcell,3)) 
      ALLOCATE(tb%tgeom%dalpmdr(tb%ng,tb%tgeom%mcell,3)) 
      ALLOCATE(tb%tgeom%dalpmdz(tb%ng,tb%tgeom%mcell,3)) 
      ALLOCATE(tb%tgeom%dalpmdrc(tb%ng,tb%tgeom%mcell,3)) 
      ALLOCATE(tb%tgeom%bigr(tb%ng,tb%tgeom%mcell)) 
      ALLOCATE(tb%tgeom%wjac(tb%ng,tb%tgeom%mcell)) 
      DO node=1,tb%ng 
        DO iv=1,3 
          tb%tgeom%alpha_arr(node,:,iv)=tb%tgeom%alpha(iv,node) 
        ENDDO 
        tb%tgeom%dalpdr(node,:,:)=tb%tgeom%alpha_x(:,1,:) 
        tb%tgeom%dalpdz(node,:,:)=tb%tgeom%alpha_y(:,1,:) 
        tb%tgeom%dalpdrc(node,:,:)=tb%tgeom%alpha_x(:,1,:) 
        IF (geom=='tor') THEN 
          DO icell=1,tb%tgeom%mcell 
            tb%tgeom%bigr(node,icell)                                   &
     &        =tb%tgeom%xs(tb%tgeom%vertex(icell,1))*                   &
     &         tb%tgeom%alpha(1,node)                                   &
     &        +tb%tgeom%xs(tb%tgeom%vertex(icell,2))*                   &
     &         tb%tgeom%alpha(2,node)                                   &
     &        +tb%tgeom%xs(tb%tgeom%vertex(icell,3))*                   &
     &         tb%tgeom%alpha(3,node)                                   
          ENDDO 
          DO iv=1,3 
            tb%tgeom%dalpdrc(node,:,iv)=                                &
     &        tb%tgeom%dalpdrc(node,:,iv)+                              &
     &        tb%tgeom%alpha_arr(node,:,iv)/tb%tgeom%bigr(node,:)       
          ENDDO 
        ENDIF 
        tb%tgeom%wjac(node,:)=tb%tgeom%area(:)*tb%wg(node) 
      ENDDO 
      IF (geom=='tor') THEN 
        tb%tgeom%wjac=tb%tgeom%wjac*tb%tgeom%bigr 
      ELSE 
        tb%tgeom%bigr=1 
      ENDIF 
!-----------------------------------------------------------------------
!     generate copies of the basis/test functions for matrices, where   
!     the square root of Jacobian times the quadrature weight is        
!     already multiplied.                                               
!                                                                       
!     test functions used for rhs computations are multiplied by a full 
!     factor of Jacobian time quadrature weight.                        
!-----------------------------------------------------------------------
      DO iv=1,3 
        tb%tgeom%alpham_arr(:,:,iv)=tb%tgeom%alpha_arr(:,:,iv)*         &
     &                              SQRT(tb%tgeom%wjac)                 
        tb%tgeom%dalpmdr(:,:,iv)=tb%tgeom%dalpdr(:,:,iv)*               &
     &                           SQRT(tb%tgeom%wjac)                    
        tb%tgeom%dalpmdz(:,:,iv)=tb%tgeom%dalpdz(:,:,iv)*               &
     &                           SQRT(tb%tgeom%wjac)                    
        tb%tgeom%dalpmdrc(:,:,iv)=tb%tgeom%dalpdrc(:,:,iv)*             &
     &                            SQRT(tb%tgeom%wjac)                   
        tb%tgeom%alpha_arr(:,:,iv)=tb%tgeom%alpha_arr(:,:,iv)*          &
     &                             tb%tgeom%wjac                        
        tb%tgeom%dalpdr(:,:,iv)=tb%tgeom%dalpdr(:,:,iv)*tb%tgeom%wjac 
        tb%tgeom%dalpdz(:,:,iv)=tb%tgeom%dalpdz(:,:,iv)*tb%tgeom%wjac 
        tb%tgeom%dalpdrc(:,:,iv)=tb%tgeom%dalpdrc(:,:,iv)*tb%tgeom%wjac 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tblock_basis_set 
!-----------------------------------------------------------------------
!     subprogram 8. tblock_basis_dealloc.                               
!     deallocates space for tblock bases.                               
!-----------------------------------------------------------------------
      SUBROUTINE tblock_basis_dealloc(tg) 
                                                                        
      TYPE(tri_linear_geom_type), INTENT(INOUT) :: tg 
                                                                        
      DEALLOCATE(tg%alpha_arr) 
      DEALLOCATE(tg%dalpdr) 
      DEALLOCATE(tg%dalpdz) 
      DEALLOCATE(tg%dalpdrc) 
      DEALLOCATE(tg%alpham_arr) 
      DEALLOCATE(tg%dalpmdr) 
      DEALLOCATE(tg%dalpmdz) 
      DEALLOCATE(tg%dalpmdrc) 
      DEALLOCATE(tg%bigr) 
      DEALLOCATE(tg%wjac) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tblock_basis_dealloc 
!-----------------------------------------------------------------------
!     subprogram 9. tblock_real_qp_update.                              
!     evaluate 2D tri_linear data at the gaussian quadrature points     
!     for this block.                                                   
!-----------------------------------------------------------------------
      SUBROUTINE tblock_real_qp_update(tl,qtl,tb) 
                                                                        
      TYPE(tri_linear_2D_type), INTENT(IN) :: tl 
      TYPE(tb_real_qp_type), INTENT(INOUT) :: qtl 
      TYPE(tblock_type), INTENT(IN) :: tb 
                                                                        
      REAL(r8), DIMENSION(tl%nqty,tb%mcell,1) :: f,fr,fz 
                                                                        
      INTEGER(i4) :: ig 
!-----------------------------------------------------------------------
!     loop over quadrature points, evaluate, and save.                  
!-----------------------------------------------------------------------
      DO ig=1,tb%ng 
        CALL tri_linear_all_eval(tl,tb%tgeom,ig,1_i4,f,fr,fz) 
        qtl%qpf(:,ig,:)=f(:,:,1) 
        qtl%qpfr(:,ig,:)=fr(:,:,1) 
        qtl%qpfz(:,ig,:)=fz(:,:,1) 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tblock_real_qp_update 
!-----------------------------------------------------------------------
!     subprogram 10. tblock_comp_qp_update.                             
!     evaluate 3D tri_linear data at the gaussian quadrature points     
!     for this block.                                                   
!-----------------------------------------------------------------------
      SUBROUTINE tblock_comp_qp_update(tl,qtl,tb) 
                                                                        
      TYPE(tri_linear_type), INTENT(IN) :: tl 
      TYPE(tb_comp_qp_type), INTENT(INOUT) :: qtl 
      TYPE(tblock_type), INTENT(IN) :: tb 
                                                                        
      COMPLEX(r8), DIMENSION(tl%nqty,tb%mcell,1,tl%nfour) :: f,fr,fz 
                                                                        
      INTEGER(i4) :: ig 
!-----------------------------------------------------------------------
!     loop over quadrature points, evaluate, and save.                  
!-----------------------------------------------------------------------
      DO ig=1,tb%ng 
        CALL tri_linear_all_eval(tl,tb%tgeom,ig,1_i4,f,fr,fz) 
        qtl%qpf(:,ig,:,:)=f(:,:,1,:) 
        qtl%qpfr(:,ig,:,:)=fr(:,:,1,:) 
        qtl%qpfz(:,ig,:,:)=fz(:,:,1,:) 
      ENDDO 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tblock_comp_qp_update 
!-----------------------------------------------------------------------
!     subprogram 11. tblock_real_qp_alloc.                              
!     allocate space for 2D data at the gaussian quadrature points for  
!     this block.                                                       
!-----------------------------------------------------------------------
      SUBROUTINE tblock_real_qp_alloc(qtl,tb,nqty) 
                                                                        
      TYPE(tb_real_qp_type), INTENT(OUT) :: qtl 
      TYPE(tblock_type), INTENT(IN) :: tb 
      INTEGER(i4), INTENT(IN) :: nqty 
                                                                        
      ALLOCATE(qtl%qpf(nqty,tb%ng,tb%mcell)) 
      ALLOCATE(qtl%qpfr(nqty,tb%ng,tb%mcell)) 
      ALLOCATE(qtl%qpfz(nqty,tb%ng,tb%mcell)) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tblock_real_qp_alloc 
!-----------------------------------------------------------------------
!     subprogram 12. tblock_comp_qp_alloc.                              
!     allocate space for 3D data at the gaussian quadrature points for  
!     this block.                                                       
!-----------------------------------------------------------------------
      SUBROUTINE tblock_comp_qp_alloc(qtl,tb,nqty,nfour) 
                                                                        
      TYPE(tb_comp_qp_type), INTENT(OUT) :: qtl 
      TYPE(tblock_type), INTENT(IN) :: tb 
      INTEGER(i4), INTENT(IN) :: nqty,nfour 
                                                                        
      ALLOCATE(qtl%qpf(nqty,tb%ng,tb%mcell,nfour)) 
      ALLOCATE(qtl%qpfr(nqty,tb%ng,tb%mcell,nfour)) 
      ALLOCATE(qtl%qpfz(nqty,tb%ng,tb%mcell,nfour)) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tblock_comp_qp_alloc 
!-----------------------------------------------------------------------
!     subprogram 13. tblock_real_qp_dealloc.                            
!     deallocate space for 2D data at the gaussian quadrature points for
!     this block.                                                       
!-----------------------------------------------------------------------
      SUBROUTINE tblock_real_qp_dealloc(qtl) 
                                                                        
      TYPE(tb_real_qp_type), INTENT(INOUT) :: qtl 
                                                                        
      DEALLOCATE(qtl%qpf) 
      IF (ASSOCIATED(qtl%qpfr)) DEALLOCATE(qtl%qpfr) 
      IF (ASSOCIATED(qtl%qpfz)) DEALLOCATE(qtl%qpfz) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tblock_real_qp_dealloc 
!-----------------------------------------------------------------------
!     subprogram 14. tblock_comp_qp_dealloc.                            
!     deallocate space for 3D data at the gaussian quadrature points for
!     this block.                                                       
!-----------------------------------------------------------------------
      SUBROUTINE tblock_comp_qp_dealloc(qtl) 
                                                                        
      TYPE(tb_comp_qp_type), INTENT(INOUT) :: qtl 
                                                                        
      DEALLOCATE(qtl%qpf) 
      IF (ASSOCIATED(qtl%qpfr)) DEALLOCATE(qtl%qpfr) 
      IF (ASSOCIATED(qtl%qpfz)) DEALLOCATE(qtl%qpfz) 
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE tblock_comp_qp_dealloc 
!-----------------------------------------------------------------------
!     close module.                                                     
!-----------------------------------------------------------------------
      END MODULE tblock 
