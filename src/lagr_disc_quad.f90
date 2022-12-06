!-----------------------------------------------------------------------
!     $Id: lagr_disc_quad.f90 6863 2019-08-19 23:55:22Z jking $
!     routines for evaluating Lagrange finite elements on blocks of
!     structured quadrilaterals, where basis functions are discontinuous
!     across element borders.
!
!     the lagr_type_mod data types are used; though, only 'interior'
!     arrays are allocated.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.
!-----------------------------------------------------------------------
!      0.  lagr_disc_mod.
!      1.  lagr_disc_bases.
!      2.  lagr_disc_3D_alloc.
!      3.  lagr_disc_3D_dealloc.
!      4.  lagr_disc_3D_eval.
!      5.  lagr_disc_3D_all_eval.
!      6.  lagr_disc_3D_assign_rsc.
!      7.  lagr_disc_3D_assign_csc.
!      8.  lagr_disc_3D_assign_laq.
!      9.  lagr_disc_3D_assign_int.
!      10. lagr_disc_3D_basis_assign_arr
!      11. lagr_disc_3D_basis_add_arr
!      12. lagr_disc_3D_basis_assign_loc
!      13. lagr_disc_3D_basis_add_loc
!      14. lagr_disc_2D_alloc.
!      15. lagr_disc_2D_dealloc.
!      16. lagr_disc_2D_eval.
!      17. lagr_disc_2D_all_eval.
!      18. lagr_disc_2D_assign_rsc.
!      19. lagr_disc_2D_assign_csc.
!      20. lagr_disc_2D_assign_laq.
!      21. lagr_disc_2D_assign_int.
!      22. lagr_disc_2D_basis_assign_arr
!      23. lagr_disc_2D_basis_add_arr
!      24. lagr_disc_2D_basis_assign_loc
!      25. lagr_disc_2D_basis_add_loc
!-----------------------------------------------------------------------
!     subprogram 0. lagr_disc_mod definition.
!     contains the subprograms and interfaces.
!-----------------------------------------------------------------------
      MODULE lagr_disc_mod
      USE lagr_type_mod
      IMPLICIT NONE

      INTEGER(i4), PARAMETER, PRIVATE :: npoly_max=20
      REAL(r8), DIMENSION(0:npoly_max), PRIVATE :: alx,aly,dalx,daly
!-----------------------------------------------------------------------
!     subprogram name interfaces
!-----------------------------------------------------------------------
      INTERFACE lagr_disc_alloc
        MODULE PROCEDURE lagr_disc_2D_alloc,lagr_disc_3D_alloc
      END INTERFACE

      INTERFACE lagr_disc_dealloc
        MODULE PROCEDURE lagr_disc_2D_dealloc,lagr_disc_3D_dealloc
      END INTERFACE

      INTERFACE lagr_disc_all_eval
        MODULE PROCEDURE lagr_disc_2D_all_eval,lagr_disc_3D_all_eval
      END INTERFACE

      INTERFACE lagr_disc_eval
        MODULE PROCEDURE lagr_disc_2D_eval,lagr_disc_3D_eval
      END INTERFACE

      INTERFACE lagr_disc_basis_assign_arr
        MODULE PROCEDURE                                                &
     &    lagr_disc_2D_basis_assign_arr,lagr_disc_3D_basis_assign_arr
      END INTERFACE

      INTERFACE lagr_disc_basis_add_arr
        MODULE PROCEDURE                                                &
     &    lagr_disc_2D_basis_add_arr,lagr_disc_3D_basis_add_arr
      END INTERFACE

      INTERFACE lagr_disc_basis_assign_loc
        MODULE PROCEDURE                                                &
     &    lagr_disc_2D_basis_assign_loc,lagr_disc_3D_basis_assign_loc
      END INTERFACE

      INTERFACE lagr_disc_basis_add_loc
        MODULE PROCEDURE                                                &
     &    lagr_disc_2D_basis_add_loc,lagr_disc_3D_basis_add_loc
      END INTERFACE

      INTERFACE lagr_disc_assignment
        MODULE PROCEDURE                                                &
     &    lagr_disc_3D_assign_csc,lagr_disc_3D_assign_rsc,              &
     &    lagr_disc_3D_assign_laq,lagr_disc_3D_assign_int,              &
     &    lagr_disc_2D_assign_csc,lagr_disc_2D_assign_rsc,              &
     &    lagr_disc_2D_assign_laq,lagr_disc_2D_assign_int
      END INTERFACE

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 1. lagr_disc_bases.
!     computes basis functions and their derivatives.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_bases(x,y,alpha,alphax,alphay,dmode)

      REAL(r8), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:), INTENT(INOUT) :: alpha,alphax,alphay
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: pd,i,j,k
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
!-----------------------------------------------------------------------
!     if derivatives are not needed, limit the computations.
!-----------------------------------------------------------------------
      IF (dmode==0) THEN
!-----------------------------------------------------------------------
!       compute basis functions with respect to the logical coordinates
!       of the unit square.  the ordering is left to right, bottom to
!       top.
!-----------------------------------------------------------------------
        k=1
        DO j=0,pd
          DO i=0,pd
            alpha(k)= alx(i)* aly(j) 
            k=k+1
          ENDDO
        ENDDO
      ELSE
!-----------------------------------------------------------------------
!       compute basis functions and derivatives.
!-----------------------------------------------------------------------
        k=1
        DO j=0,pd
          DO i=0,pd
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
      END SUBROUTINE lagr_disc_bases
!-----------------------------------------------------------------------
!     subprogram 2. lagr_disc_3D_alloc.
!     allocates space for lagr_disc_3D_type.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_3D_alloc(laq,mx,my,nqty,nfour,poly_degree,   &
     &                              name,title)

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,nfour,poly_degree
      CHARACTER(*), INTENT(IN), OPTIONAL :: name
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title
      TYPE(lagr_quad_type), INTENT(OUT) :: laq

      INTEGER(i4) :: ix,iy,ib,ip,pdp1
      REAL(r8), DIMENSION(0:poly_degree) :: x_node
!-----------------------------------------------------------------------
!     store grid, vector, and fourier series dimensions, and set the
!     number of interior basis functions.
!-----------------------------------------------------------------------
      laq%mx=mx
      laq%my=my
      laq%nqty=nqty
      laq%nfour=nfour
      laq%n_side=0
      laq%n_int=(poly_degree+1_i4)**2
!-----------------------------------------------------------------------
!     allocate space.
!-----------------------------------------------------------------------
      ALLOCATE(laq%fsi(nqty,laq%n_int,1:mx,1:my,nfour))
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
!     all discontinuous bases are element-centered, so ix0=iy0=1.
!     this descriptor information is completed for compatibility.
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
!           CALL lagr_disc_eval(laq,dx_data,dy_data,1_i4)
!           xxx=laq%f(yyy,zzz)
!         ENDDO
!       ENDDO
!     ENDDO
!-----------------------------------------------------------------------
      pdp1=poly_degree+1
      ALLOCATE(laq%ix0(laq%n_int))
      ALLOCATE(laq%iy0(laq%n_int))
      ALLOCATE(laq%dx(laq%n_int))
      ALLOCATE(laq%dy(laq%n_int))
      CALL poly_nodes(poly_degree,x_node)
      DO ib=1,laq%n_int
        laq%ix0(ib)=1
        laq%iy0(ib)=1
        ix=MODULO(ib-1_i4,pdp1)
        iy=(ib-1_i4)/pdp1
        laq%dx(ib)=x_node(ix)
        laq%dy(ib)=x_node(iy)
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_3D_alloc
!-----------------------------------------------------------------------
!     subprogram 3. lagr_disc_3D_dealloc.
!     deallocates space for lagr_disc_3D_type.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_3D_dealloc(laq)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
!-----------------------------------------------------------------------
!     deallocate space.
!-----------------------------------------------------------------------
      DEALLOCATE(laq%fsi)
      DEALLOCATE(laq%title)
      DEALLOCATE(laq%f)
      DEALLOCATE(laq%fx)
      DEALLOCATE(laq%fy)
      DEALLOCATE(laq%ix0)
      DEALLOCATE(laq%iy0)
      DEALLOCATE(laq%dx)
      DEALLOCATE(laq%dy)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_3D_dealloc
!-----------------------------------------------------------------------
!     subprogram 4. lagr_disc_3D_eval.
!     evaluates complex lagr_disc quantities at a single point within a
!     grid block.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_3D_eval(laq,x,y,dmode)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      REAL(r8), INTENT(IN) :: x,y
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: ix,iy,pd,i,j,k,im
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
      pd=NINT(SQRT(REAL(laq%n_int)))-1_i4
      CALL lagr_1D(pd,x-ix,alx,dalx,dmode)
      CALL lagr_1D(pd,y-iy,aly,daly,dmode)
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.
!-----------------------------------------------------------------------
      laq%f=0._r8
      k=1
      DO j=0,pd
        DO i=0,pd
          laq%f=laq%f+laq%fsi(:,k,ix+1,iy+1,:)*alx(i)*aly(j)
          k=k+1
        ENDDO
      ENDDO
      IF (dmode==0) RETURN
      laq%fx=0._r8
      laq%fy=0._r8
      k=1
      DO j=0,pd
        DO i=0,pd
          laq%fx=laq%fx+laq%fsi(:,k,ix+1,iy+1,:)*dalx(i)* aly(j)
          laq%fy=laq%fy+laq%fsi(:,k,ix+1,iy+1,:)* alx(i)*daly(j)
          k=k+1
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_3D_eval
!-----------------------------------------------------------------------
!     subprogram 5. lagr_disc_3D_all_eval.
!     evaluates complex lagr_disc quantities in all elements in a grid
!     block for equal spacing.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_3D_all_eval(laq,x,y,f,fx,fy,dmode)

      TYPE(lagr_quad_type), INTENT(IN) :: laq
      REAL(r8), INTENT(IN) :: x,y
      COMPLEX(r8), INTENT(OUT), DIMENSION(:,:,:,:) :: f,fx,fy
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: mx,my,mx1,my1
      INTEGER(i4) :: pd,i,j,k,im,ix,iy
      REAL(r8), DIMENSION(laq%n_int) :: alpha,dalpdx,dalpdy
!-----------------------------------------------------------------------
!     compute index limits.
!-----------------------------------------------------------------------
      mx=laq%mx
      mx1=mx-1
      my=laq%my
      my1=my-1
      pd=NINT(SQRT(REAL(laq%n_int)))-1_i4
      CALL lagr_disc_bases(x,y,alpha,dalpdx,dalpdy,dmode)
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.
!-----------------------------------------------------------------------
      IF (dmode==0) THEN
        DO im=1,laq%nfour
          DO iy=1,my
            DO ix=1,mx
              DO i=1,laq%nqty
                f(i,ix,iy,im)=SUM(laq%fsi(i,:,ix,iy,im)*alpha(:))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO im=1,laq%nfour
          DO iy=1,my
            DO ix=1,mx
              DO i=1,laq%nqty
                f(i,ix,iy,im)=SUM(laq%fsi(i,:,ix,iy,im)*alpha(:))
                fx(i,ix,iy,im)=SUM(laq%fsi(i,:,ix,iy,im)*dalpdx(:))
                fy(i,ix,iy,im)=SUM(laq%fsi(i,:,ix,iy,im)*dalpdy(:))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_3D_all_eval
!-----------------------------------------------------------------------
!     subprogram 6. lagr_disc_3D_assign_rsc.
!     assign a real scalar value to a complex lagrange quad structure.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_3D_assign_rsc(laq,rscalar)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      REAL(r8), INTENT(IN) :: rscalar

      laq%fsi=rscalar
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_3D_assign_rsc
!-----------------------------------------------------------------------
!     subprogram 7. lagr_disc_3D_assign_csc.
!     assign a complex scalar value to a complex lagrange quad
!     structure.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_3D_assign_csc(laq,cscalar)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      COMPLEX(r8), INTENT(IN) :: cscalar

      laq%fsi=cscalar
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_3D_assign_csc
!-----------------------------------------------------------------------
!     subprogram 8. lagr_disc_3D_assign_laq.
!     set one complex lagrange quad structure equal to another.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_3D_assign_laq(laq1,laq2)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq1
      TYPE(lagr_quad_type), INTENT(IN) :: laq2

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      laq1%fsi(:,:,:,:,:)=laq2%fsi(:,:,:,:,:)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_3D_assign_laq
!-----------------------------------------------------------------------
!     subprogram 9. lagr_disc_3D_assign_int.
!     assign a integer value to a complex lagrange quad structure.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_3D_assign_int(laq,int)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      INTEGER(i4), INTENT(IN) :: int

      laq%fsi=int
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_3D_assign_int
!-----------------------------------------------------------------------
!     subprogram 10. lagr_disc_3D_basis_assign_arr
!     assign data into coefficient arrays for one basis function.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_3D_basis_assign_arr(laq,data,ibasis)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      laq%fsi(:,ibasis,:,:,:)=data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_3D_basis_assign_arr
!-----------------------------------------------------------------------
!     subprogram 11. lagr_disc_3D_basis_add_arr
!     add data into coefficient arrays for one basis function.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_3D_basis_add_arr(laq,data,ibasis)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      laq%fsi(:,ibasis,:,:,:)=laq%fsi(:,ibasis,:,:,:)+data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_3D_basis_add_arr
!-----------------------------------------------------------------------
!     subprogram 12. lagr_disc_3D_basis_assign_loc
!     assign data into coefficient arrays for one basis function.
!
!     this is a local version of lagr_disc_basis_assign_arr, where the
!     the data is located at a given poloidal and fourier indices
!     triplet, only.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_3D_basis_assign_loc(laq,data,ibasis,ix,iy,im)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy,im
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      laq%fsi(:,ibasis,ix,iy,im)=data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_3D_basis_assign_loc
!-----------------------------------------------------------------------
!     subprogram 13. lagr_disc_3D_basis_add_loc
!     add data into coefficient arrays for one basis function.
!
!     this is a local version of lagr_disc_basis_add_arr, where the
!     the data is located at a given poloidal and fourier indices
!     triplet, only.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_3D_basis_add_loc(laq,data,ibasis,ix,iy,im)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy,im
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      laq%fsi(:,ibasis,ix,iy,im)=laq%fsi(:,ibasis,ix,iy,im)+data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_3D_basis_add_loc
!-----------------------------------------------------------------------
!     subprogram 14. lagr_disc_2D_alloc.
!     allocates space for lagr_disc_2D_type.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_2D_alloc(laq,mx,my,nqty,poly_degree,         &
     &                              name,title)

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,poly_degree
      CHARACTER(*), INTENT(IN), OPTIONAL :: name
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title
      TYPE(lagr_quad_2D_type), INTENT(OUT) :: laq

      INTEGER(i4) :: ix,iy,ib,ip,pdp1
      REAL(r8), DIMENSION(0:poly_degree) :: x_node
!-----------------------------------------------------------------------
!     store grid, vector, and fourier series dimensions, and set the
!     number of interior basis functions.
!-----------------------------------------------------------------------
      laq%mx=mx
      laq%my=my
      laq%nqty=nqty
      laq%n_side=0
      laq%n_int=(poly_degree+1_i4)**2
!-----------------------------------------------------------------------
!     allocate space.
!-----------------------------------------------------------------------
      ALLOCATE(laq%fsi(nqty,laq%n_int,1:mx,1:my))
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
!     all discontinuous bases are element-centered, so ix0=iy0=1.
!     this descriptor information is completed for compatibility.
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
!           CALL lagr_disc_eval(laq,dx_data,dy_data,1_i4)
!           xxx=laq%f(yyy,zzz)
!         ENDDO
!       ENDDO
!     ENDDO
!-----------------------------------------------------------------------
      pdp1=poly_degree+1
      ALLOCATE(laq%ix0(laq%n_int))
      ALLOCATE(laq%iy0(laq%n_int))
      ALLOCATE(laq%dx(laq%n_int))
      ALLOCATE(laq%dy(laq%n_int))
      CALL poly_nodes(poly_degree,x_node)
      DO ib=1,laq%n_int
        laq%ix0(ib)=1
        laq%iy0(ib)=1
        ix=MODULO(ib-1_i4,pdp1)
        iy=(ib-1_i4)/pdp1
        laq%dx(ib)=x_node(ix)
        laq%dy(ib)=x_node(iy)
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_2D_alloc
!-----------------------------------------------------------------------
!     subprogram 15. lagr_disc_2D_dealloc.
!     deallocates space for lagr_disc_2D_type.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_2D_dealloc(laq)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
!-----------------------------------------------------------------------
!     deallocate space.
!-----------------------------------------------------------------------
      DEALLOCATE(laq%fsi)
      DEALLOCATE(laq%title)
      DEALLOCATE(laq%f)
      DEALLOCATE(laq%fx)
      DEALLOCATE(laq%fy)
      DEALLOCATE(laq%ix0)
      DEALLOCATE(laq%iy0)
      DEALLOCATE(laq%dx)
      DEALLOCATE(laq%dy)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_2D_dealloc
!-----------------------------------------------------------------------
!     subprogram 16. lagr_disc_2D_eval.
!     evaluates real lagr_disc quantities at a single point within a
!     grid block.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_2D_eval(laq,x,y,dmode)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      REAL(r8), INTENT(IN) :: x,y
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: ix,iy,pd,i,j,k,im
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
      pd=NINT(SQRT(REAL(laq%n_int)))-1_i4
      CALL lagr_1D(pd,x-ix,alx,dalx,dmode)
      CALL lagr_1D(pd,y-iy,aly,daly,dmode)
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.
!-----------------------------------------------------------------------
      laq%f=0._r8
      k=1
      DO j=0,pd
        DO i=0,pd
          laq%f=laq%f+laq%fsi(:,k,ix+1,iy+1)*alx(i)*aly(j)
          k=k+1
        ENDDO
      ENDDO
      IF (dmode==0) RETURN
      laq%fx=0._r8
      laq%fy=0._r8
      k=1
      DO j=0,pd
        DO i=0,pd
          laq%fx=laq%fx+laq%fsi(:,k,ix+1,iy+1)*dalx(i)* aly(j)
          laq%fy=laq%fy+laq%fsi(:,k,ix+1,iy+1)* alx(i)*daly(j)
          k=k+1
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_2D_eval
!-----------------------------------------------------------------------
!     subprogram 17. lagr_disc_2D_all_eval.
!     evaluates real lagr_disc quantities in all elements in a grid
!     block for equal spacing.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_2D_all_eval(laq,x,y,f,fx,fy,dmode)

      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq
      REAL(r8), INTENT(IN) :: x,y
      REAL(r8), INTENT(OUT), DIMENSION(:,:,:) :: f,fx,fy
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: mx,my,mx1,my1
      INTEGER(i4) :: pd,i,j,k,ix,iy
      REAL(r8), DIMENSION(laq%n_int) :: alpha,dalpdx,dalpdy
!-----------------------------------------------------------------------
!     compute index limits.
!-----------------------------------------------------------------------
      mx=laq%mx
      mx1=mx-1
      my=laq%my
      my1=my-1
      pd=NINT(SQRT(REAL(laq%n_int)))-1_i4
      CALL lagr_disc_bases(x,y,alpha,dalpdx,dalpdy,dmode)
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.
!-----------------------------------------------------------------------
      IF (dmode==0) THEN
        DO iy=1,my
          DO ix=1,mx
            DO i=1,laq%nqty
              f(i,ix,iy)=SUM(laq%fsi(i,:,ix,iy)*alpha(:))
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO iy=1,my
          DO ix=1,mx
            DO i=1,laq%nqty
              f(i,ix,iy)=SUM(laq%fsi(i,:,ix,iy)*alpha(:))
              fx(i,ix,iy)=SUM(laq%fsi(i,:,ix,iy)*dalpdx(:))
              fy(i,ix,iy)=SUM(laq%fsi(i,:,ix,iy)*dalpdy(:))
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_2D_all_eval
!-----------------------------------------------------------------------
!     subprogram 18. lagr_disc_2D_assign_rsc.
!     assign a real scalar value to a real lagrange quad structure.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_2D_assign_rsc(laq,rscalar)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      REAL(r8), INTENT(IN) :: rscalar

      laq%fsi=rscalar
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_2D_assign_rsc
!-----------------------------------------------------------------------
!     subprogram 19. lagr_disc_2D_assign_csc.
!     assign a real scalar value to a real lagrange quad
!     structure.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_2D_assign_csc(laq,cscalar)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      COMPLEX(r8), INTENT(IN) :: cscalar

      laq%fsi=cscalar
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_2D_assign_csc
!-----------------------------------------------------------------------
!     subprogram 20. lagr_disc_2D_assign_laq.
!     set one real lagrange quad structure equal to another.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_2D_assign_laq(laq1,laq2)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq1
      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq2

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      laq1%fsi(:,:,:,:)=laq2%fsi(:,:,:,:)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_2D_assign_laq
!-----------------------------------------------------------------------
!     subprogram 21. lagr_disc_2D_assign_int.
!     assign a integer value to a real lagrange quad structure.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_2D_assign_int(laq,int)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      INTEGER(i4), INTENT(IN) :: int

      laq%fsi=int
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_2D_assign_int
!-----------------------------------------------------------------------
!     subprogram 22. lagr_disc_2D_basis_assign_arr
!     assign data into coefficient arrays for one basis function.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_2D_basis_assign_arr(laq,data,ibasis)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      laq%fsi(:,ibasis,:,:)=data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_2D_basis_assign_arr
!-----------------------------------------------------------------------
!     subprogram 23. lagr_disc_2D_basis_add_arr
!     add data into coefficient arrays for one basis function.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_2D_basis_add_arr(laq,data,ibasis)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      laq%fsi(:,ibasis,:,:)=laq%fsi(:,ibasis,:,:)+data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_2D_basis_add_arr
!-----------------------------------------------------------------------
!     subprogram 24. lagr_disc_2D_basis_assign_loc
!     assign data into coefficient arrays for one basis function.
!
!     this is a local version of lagr_disc_basis_assign_arr, where the
!     the data is located at given poloidal indices, only.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_2D_basis_assign_loc(laq,data,ibasis,ix,iy)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      REAL(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      laq%fsi(:,ibasis,ix,iy)=data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_2D_basis_assign_loc
!-----------------------------------------------------------------------
!     subprogram 25. lagr_disc_2D_basis_add_loc
!     add data into coefficient arrays for one basis function.
!
!     this is a local version of lagr_disc_basis_add_arr, where the
!     the data is located at given poloidal indices only.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_disc_2D_basis_add_loc(laq,data,ibasis,ix,iy)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      REAL(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      laq%fsi(:,ibasis,ix,iy)=laq%fsi(:,ibasis,ix,iy)+data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_disc_2D_basis_add_loc
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      END MODULE lagr_disc_mod
