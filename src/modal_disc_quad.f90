!-----------------------------------------------------------------------
!     $Id: modal_disc_quad.f90 6329 2018-12-18 19:29:48Z jking $
!     routines for evaluating discontinuous modal bases in blocks of
!     structured quadrilateral elements.
!
!     unlike lagr_disc_quad, these modal bases have their own type
!     definition.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.
!-----------------------------------------------------------------------
!      0.  modal_type_mod.
!      0.1 modal_disc_mod.
!      1.  modal_disc_bases.
!      2.  modal_disc_3D_alloc.
!      3.  modal_disc_3D_dealloc.
!      4.  modal_disc_3D_eval.
!      5.  modal_disc_3D_all_eval.
!      6.  modal_disc_3D_assign_rsc.
!      7.  modal_disc_3D_assign_csc.
!      8.  modal_disc_3D_assign_modq.
!      9.  modal_disc_3D_assign_int.
!      10. modal_disc_3D_basis_assign_arr
!      11. modal_disc_3D_basis_add_arr
!      12. modal_disc_3D_basis_assign_loc
!      13. modal_disc_3D_basis_add_loc
!      14. modal_disc_2D_alloc.
!      15. modal_disc_2D_dealloc.
!      16. modal_disc_2D_eval.
!      17. modal_disc_2D_all_eval.
!      18. modal_disc_2D_assign_rsc.
!      19. modal_disc_2D_assign_csc.
!      20. modal_disc_2D_assign_modq.
!      21. modal_disc_2D_assign_int.
!      22. modal_disc_2D_basis_assign_arr
!      23. modal_disc_2D_basis_add_arr
!      24. modal_disc_2D_basis_assign_loc
!      25. modal_disc_2D_basis_add_loc
!-----------------------------------------------------------------------
!     subprogram 0. modal_type_mod definition.
!     contains the F90 user-type definition for modal bases in
!     quadrilateral elements.  the bases only have arrays for element-
!     interior coefficients, and the starting indices are always 1.
!
!     a distinction from standard bases is that the polynomials here
!     can be incomplete representations.  the expansion can be described
!     as the outer product of Legendre polynomials of degrees pdmin to
!     pdmax in the logical coordinate x and the complete basis from
!     0 to pd in y, 'unioned' with the same with x and y swapped.
!     it is assumed that pdmax<=pd.
!
!     for example, if pdmin=pdmax=3 and pd=4, bases include the outer
!     product of the cubic Legendre polynomial in x with the full
!     quartic expansion in y, plus the cubic Legendre polynomial in y
!     with the full quartic expansion in x, without duplication.
!-----------------------------------------------------------------------
      MODULE modal_type_mod
      USE local
      IMPLICIT NONE

      TYPE :: modal_quad_type
        INTEGER(i4) :: mx,my,nqty,nfour,n_int,pd,pdmin,pdmax
        INTEGER(i4), DIMENSION(:), POINTER :: ix0,iy0
        COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: fsi
        COMPLEX(r8), DIMENSION(:,:), POINTER :: f,fx,fy
        CHARACTER(6), DIMENSION(:), POINTER :: title
        CHARACTER(6) :: name
      END TYPE modal_quad_type

      TYPE :: modal_quad_2D_type
        INTEGER(i4) :: mx,my,nqty,n_int,pd,pdmin,pdmax
        INTEGER(i4), DIMENSION(:), POINTER :: ix0,iy0
        REAL(r8), DIMENSION(:,:,:,:), POINTER :: fsi
        REAL(r8), DIMENSION(:), POINTER :: f,fx,fy
        CHARACTER(6), DIMENSION(:), POINTER :: title
        CHARACTER(6) :: name
      END TYPE modal_quad_2D_type

      END MODULE modal_type_mod
!-----------------------------------------------------------------------
!     subprogram 0.1 modal_disc_mod definition.
!     contains the subprograms and interfaces.
!-----------------------------------------------------------------------
      MODULE modal_disc_mod
      USE modal_type_mod
      IMPLICIT NONE

      INTEGER(i4), PARAMETER, PRIVATE :: npoly_max=20
!-----------------------------------------------------------------------
!     subprogram name interfaces
!-----------------------------------------------------------------------
      INTERFACE modal_disc_alloc
        MODULE PROCEDURE modal_disc_2D_alloc,modal_disc_3D_alloc
      END INTERFACE

      INTERFACE modal_disc_dealloc
        MODULE PROCEDURE modal_disc_2D_dealloc,modal_disc_3D_dealloc
      END INTERFACE

      INTERFACE modal_disc_all_eval
        MODULE PROCEDURE modal_disc_2D_all_eval,modal_disc_3D_all_eval
      END INTERFACE

      INTERFACE modal_disc_eval
        MODULE PROCEDURE modal_disc_2D_eval,modal_disc_3D_eval
      END INTERFACE

      INTERFACE modal_disc_basis_assign_arr
        MODULE PROCEDURE                                                &
     &    modal_disc_2D_basis_assign_arr,modal_disc_3D_basis_assign_arr
      END INTERFACE

      INTERFACE modal_disc_basis_add_arr
        MODULE PROCEDURE                                                &
     &    modal_disc_2D_basis_add_arr,modal_disc_3D_basis_add_arr
      END INTERFACE

      INTERFACE modal_disc_basis_assign_loc
        MODULE PROCEDURE                                                &
     &    modal_disc_2D_basis_assign_loc,modal_disc_3D_basis_assign_loc
      END INTERFACE

      INTERFACE modal_disc_basis_add_loc
        MODULE PROCEDURE                                                &
     &    modal_disc_2D_basis_add_loc,modal_disc_3D_basis_add_loc
      END INTERFACE

      INTERFACE modal_disc_assignment
        MODULE PROCEDURE                                                &
     &    modal_disc_3D_assign_csc,modal_disc_3D_assign_rsc,            &
     &    modal_disc_3D_assign_modq,modal_disc_3D_assign_int,           &
     &    modal_disc_2D_assign_csc,modal_disc_2D_assign_rsc,            &
     &    modal_disc_2D_assign_modq,modal_disc_2D_assign_int    
      END INTERFACE

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 1. modal_disc_bases.
!     computes modal basis functions and their derivatives at a given
!     logical position (x,y) within an element, based on nimrod's
!     0<=x,y<=1.  note that the minimum polynomial degree (pdmin)
!     and the maximum degree in each coordinate (pd) must be specified.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_bases(x,y,alpha,alphax,alphay,dmode,        &
     &                            pd,pdmin,pdmax)

      REAL(r8), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:), INTENT(INOUT) :: alpha,alphax,alphay
      INTEGER(i4), INTENT(IN) :: dmode,pd,pdmin,pdmax

      INTEGER(i4) :: i,j,k
      REAL(r8) :: xstd,ystd
      REAL(r8), DIMENSION(0:npoly_max) :: alx,aly,dalx,daly
!-----------------------------------------------------------------------
!     first determine the logical coordinates in the standard,
!     -1<=xstd,ystd,<=+1 range and set 1D basis functions.
!-----------------------------------------------------------------------
      IF (pdmax>pd) CALL nim_stop("Modal_disc_bases: pdmax>pd")
      xstd=2._r8*x-1._r8
      ystd=2._r8*y-1._r8
      DO i=0,pd
        CALL leg_poly1_sub(xstd,i,alx(i),dalx(i))
        CALL leg_poly1_sub(ystd,i,aly(i),daly(i))
      ENDDO
!-----------------------------------------------------------------------
!     if derivatives are not needed, limit the computations.
!-----------------------------------------------------------------------
      IF (dmode==0) THEN
!-----------------------------------------------------------------------
!       compute basis functions with respect to the logical coordinates
!       of the unit square.  the basis ordering starts from pdmin in y
!       with increasing degree in x from 0 to pd then completes from
!       pdmin in x and increasing degree in y.
!-----------------------------------------------------------------------
        k=1
        DO j=pdmin,pdmax
          DO i=0,pd
            alpha(k)=alx(i)*aly(j) 
            k=k+1
          ENDDO
        ENDDO
        DO i=pdmin,pdmax
          DO j=0,pdmin-1
            alpha(k)=alx(i)*aly(j) 
            k=k+1
          ENDDO
          DO j=pdmax+1,pd
            alpha(k)=alx(i)*aly(j) 
            k=k+1
          ENDDO
        ENDDO
      ELSE
!-----------------------------------------------------------------------
!       compute basis functions and derivatives.
!-----------------------------------------------------------------------
        k=1
        DO j=pdmin,pdmax
          DO i=0,pd
            alpha(k)=  alx(i)* aly(j) 
            alphax(k)=dalx(i)* aly(j) 
            alphay(k)= alx(i)*daly(j) 
            k=k+1
          ENDDO
        ENDDO
        DO i=pdmin,pdmax
          DO j=0,pdmin-1
            alpha(k)=  alx(i)* aly(j) 
            alphax(k)=dalx(i)* aly(j) 
            alphay(k)= alx(i)*daly(j) 
            k=k+1
          ENDDO
          DO j=pdmax+1,pd
            alpha(k)=  alx(i)* aly(j) 
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
      END SUBROUTINE modal_disc_bases
!-----------------------------------------------------------------------
!     subprogram 2. modal_disc_3D_alloc.
!     allocates space for modal_disc_3D_type.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_alloc(modq,mx,my,nqty,nfour,pd,pdmin,    &
     &                               pdmax,name,title)

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,nfour,pd,pdmin,pdmax
      CHARACTER(*), INTENT(IN), OPTIONAL :: name
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title
      TYPE(modal_quad_type), INTENT(OUT) :: modq
!-----------------------------------------------------------------------
!     store grid, vector, and fourier series dimensions, and set the
!     number of interior basis functions.
!-----------------------------------------------------------------------
      IF (pdmax>pd) CALL nim_stop("Modal_disc_alloc: pdmax>pd")
      modq%mx=mx
      modq%my=my
      modq%nqty=nqty
      modq%nfour=nfour
      modq%pd=pd
      modq%pdmin=pdmin
      modq%pdmax=pdmax
      modq%n_int=(pdmax-pdmin+1)*(2*pd+pdmin-pdmax+1)
!-----------------------------------------------------------------------
!     allocate space.
!-----------------------------------------------------------------------
      ALLOCATE(modq%fsi(nqty,modq%n_int,1:mx,1:my,nfour))
      ALLOCATE(modq%title(nqty))
      ALLOCATE(modq%f(nqty,nfour))
      ALLOCATE(modq%fx(nqty,nfour))
      ALLOCATE(modq%fy(nqty,nfour))
!-----------------------------------------------------------------------
!     character descriptors, if present in input.
!-----------------------------------------------------------------------
      IF (PRESENT(name)) modq%name=name
      IF (PRESENT(title)) THEN
        IF (SIZE(title)==nqty) THEN
          modq%title=title
        ELSE
          modq%title=title(1)
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     all discontinuous bases are element-centered, so ix0=iy0=1.
!     for modal elements, their are no nodal positions (dx and dy in
!     nimrod's nodal representations).
!-----------------------------------------------------------------------
      ALLOCATE(modq%ix0(modq%n_int))
      ALLOCATE(modq%iy0(modq%n_int))
      modq%ix0=1
      modq%iy0=1
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_alloc
!-----------------------------------------------------------------------
!     subprogram 3. modal_disc_3D_dealloc.
!     deallocates space for modal_disc_3D_type.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_dealloc(modq)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
!-----------------------------------------------------------------------
!     deallocate space.
!-----------------------------------------------------------------------
      DEALLOCATE(modq%fsi)
      DEALLOCATE(modq%title)
      DEALLOCATE(modq%f)
      DEALLOCATE(modq%fx)
      DEALLOCATE(modq%fy)
      DEALLOCATE(modq%ix0)
      DEALLOCATE(modq%iy0)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_dealloc
!-----------------------------------------------------------------------
!     subprogram 4. modal_disc_3D_eval.
!     evaluates complex modal_disc quantities at a single point within a
!     grid block.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_eval(modq,x,y,dmode)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      REAL(r8), INTENT(IN) :: x,y
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: ix,iy,i,j,k,im
      REAL(r8) :: xstd,ystd
      REAL(r8), DIMENSION(0:npoly_max) :: alx,aly,dalx,daly
!-----------------------------------------------------------------------
!     find the interval, and compute 1D basis coefficients with logical
!     coordinates for a standard polynomial in -1<=xstd,ystd<=+1.
!-----------------------------------------------------------------------
      ix=MAX(MIN(INT(x),modq%mx-1),0_i4)
      iy=MAX(MIN(INT(y),modq%my-1),0_i4)
      xstd=2._r8*(x-ix)-1._r8
      ystd=2._r8*(y-iy)-1._r8
      DO i=0,modq%pd
        CALL leg_poly1_sub(xstd,i,alx(i),dalx(i))
        CALL leg_poly1_sub(ystd,i,aly(i),daly(i))
      ENDDO
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.  it is
!     important that the basis-function order is the same as defined
!     in modal_disc_bases.
!-----------------------------------------------------------------------
      k=1
      modq%f=0._r8
      modq%fx=0._r8
      modq%fy=0._r8

      IF (dmode==0) THEN
        DO j=modq%pdmin,modq%pdmax
          DO i=0,modq%pd
            modq%f=modq%f+modq%fsi(:,k,ix+1,iy+1,:)*alx(i)*aly(j)
            k=k+1
          ENDDO
        ENDDO
        DO i=modq%pdmin,modq%pdmax
          DO j=0,modq%pdmin-1
            modq%f=modq%f+modq%fsi(:,k,ix+1,iy+1,:)*alx(i)*aly(j)
            k=k+1
          ENDDO
          DO j=modq%pdmax+1,modq%pd
            modq%f=modq%f+modq%fsi(:,k,ix+1,iy+1,:)*alx(i)*aly(j)
            k=k+1
          ENDDO
        ENDDO
      ELSE
        DO j=modq%pdmin,modq%pdmax
          DO i=0,modq%pd
            modq%f =modq%f +modq%fsi(:,k,ix+1,iy+1,:)* alx(i)* aly(j)
            modq%fx=modq%fx+modq%fsi(:,k,ix+1,iy+1,:)*dalx(i)* aly(j)
            modq%fy=modq%fy+modq%fsi(:,k,ix+1,iy+1,:)* alx(i)*daly(j)
            k=k+1
          ENDDO
        ENDDO
        DO i=modq%pdmin,modq%pdmax
          DO j=0,modq%pdmin-1
            modq%f =modq%f +modq%fsi(:,k,ix+1,iy+1,:)* alx(i)* aly(j)
            modq%fx=modq%fx+modq%fsi(:,k,ix+1,iy+1,:)*dalx(i)* aly(j)
            modq%fy=modq%fy+modq%fsi(:,k,ix+1,iy+1,:)* alx(i)*daly(j)
            k=k+1
          ENDDO
          DO j=modq%pdmax+1,modq%pd
            modq%f =modq%f +modq%fsi(:,k,ix+1,iy+1,:)* alx(i)* aly(j)
            modq%fx=modq%fx+modq%fsi(:,k,ix+1,iy+1,:)*dalx(i)* aly(j)
            modq%fy=modq%fy+modq%fsi(:,k,ix+1,iy+1,:)* alx(i)*daly(j)
            k=k+1
          ENDDO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_eval
!-----------------------------------------------------------------------
!     subprogram 5. modal_disc_3D_all_eval.
!     evaluates complex modal_disc quantities in all elements in a grid
!     block for equal spacing.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_all_eval(modq,x,y,f,fx,fy,dmode)

      TYPE(modal_quad_type), INTENT(IN) :: modq
      REAL(r8), INTENT(IN) :: x,y
      COMPLEX(r8), INTENT(OUT), DIMENSION(:,:,:,:) :: f,fx,fy
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: i,j,k,im,ix,iy
      REAL(r8), DIMENSION(modq%n_int) :: alpha,dalpdx,dalpdy
!-----------------------------------------------------------------------
!     evaluate the bases at the specified offset (x,y) within each
!     element.
!-----------------------------------------------------------------------
      CALL modal_disc_bases(x,y,alpha,dalpdx,dalpdy,dmode,              &
     &                      modq%pd,modq%pdmin,modq%pdmax)
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.
!-----------------------------------------------------------------------
      IF (dmode==0) THEN
        DO im=1,modq%nfour
          DO iy=1,modq%my
            DO ix=1,modq%mx
              DO i=1,modq%nqty
                f(i,ix,iy,im)=SUM(modq%fsi(i,:,ix,iy,im)*alpha(:))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO im=1,modq%nfour
          DO iy=1,modq%my
            DO ix=1,modq%mx
              DO i=1,modq%nqty
                f(i,ix,iy,im)=SUM(modq%fsi(i,:,ix,iy,im)*alpha(:))
                fx(i,ix,iy,im)=SUM(modq%fsi(i,:,ix,iy,im)*dalpdx(:))
                fy(i,ix,iy,im)=SUM(modq%fsi(i,:,ix,iy,im)*dalpdy(:))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_all_eval
!-----------------------------------------------------------------------
!     subprogram 6. modal_disc_3D_assign_rsc.
!     assign a real scalar value to a complex modalange quad structure.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_assign_rsc(modq,rscalar)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      REAL(r8), INTENT(IN) :: rscalar

      modq%fsi=rscalar
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_assign_rsc
!-----------------------------------------------------------------------
!     subprogram 7. modal_disc_3D_assign_csc.
!     assign a complex scalar value to a complex modalange quad
!     structure.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_assign_csc(modq,cscalar)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      COMPLEX(r8), INTENT(IN) :: cscalar

      modq%fsi=cscalar
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_assign_csc
!-----------------------------------------------------------------------
!     subprogram 8. modal_disc_3D_assign_modq.
!     set one complex modalange quad structure equal to another.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_assign_modq(modq1,modq2)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq1
      TYPE(modal_quad_type), INTENT(IN) :: modq2

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      modq1%fsi(:,:,:,:,:)=modq2%fsi(:,:,:,:,:)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_assign_modq
!-----------------------------------------------------------------------
!     subprogram 9. modal_disc_3D_assign_int.
!     assign a integer value to a complex modalange quad structure.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_assign_int(modq,int)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      INTEGER(i4), INTENT(IN) :: int

      modq%fsi=int
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_assign_int
!-----------------------------------------------------------------------
!     subprogram 10. modal_disc_3D_basis_assign_arr
!     assign data into coefficient arrays for one basis function.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_basis_assign_arr(modq,data,ibasis)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      modq%fsi(:,ibasis,:,:,:)=data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_basis_assign_arr
!-----------------------------------------------------------------------
!     subprogram 11. modal_disc_3D_basis_add_arr
!     add data into coefficient arrays for one basis function.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_basis_add_arr(modq,data,ibasis)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      modq%fsi(:,ibasis,:,:,:)=modq%fsi(:,ibasis,:,:,:)+data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_basis_add_arr
!-----------------------------------------------------------------------
!     subprogram 12. modal_disc_3D_basis_assign_loc
!     assign data into coefficient arrays for one basis function.
!
!     this is a local version of modal_disc_basis_assign_arr, where the
!     the data is located at a given poloidal and fourier indices
!     triplet, only.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_basis_assign_loc(modq,data,ibasis,ix,iy, &
     &                                          im)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy,im
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      modq%fsi(:,ibasis,ix,iy,im)=data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_basis_assign_loc
!-----------------------------------------------------------------------
!     subprogram 13. modal_disc_3D_basis_add_loc
!     add data into coefficient arrays for one basis function.
!
!     this is a local version of modal_disc_basis_add_arr, where the
!     the data is located at a given poloidal and fourier indices
!     triplet, only.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_basis_add_loc(modq,data,ibasis,ix,iy,im)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy,im
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      modq%fsi(:,ibasis,ix,iy,im)=modq%fsi(:,ibasis,ix,iy,im)+data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_basis_add_loc
!-----------------------------------------------------------------------
!     subprogram 14. modal_disc_2D_alloc.
!     allocates space for modal_disc_2D_type.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_alloc(modq,mx,my,nqty,pd,pdmin,pdmax,    &
     &                               name,title)

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,pd,pdmin,pdmax
      CHARACTER(*), INTENT(IN), OPTIONAL :: name
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title
      TYPE(modal_quad_2D_type), INTENT(OUT) :: modq
!-----------------------------------------------------------------------
!     store grid and vector dimensions, and set the
!     number of interior basis functions.
!-----------------------------------------------------------------------
      IF (pdmax>pd) CALL nim_stop("Modal_disc_alloc: pdmax>pd")
      modq%mx=mx
      modq%my=my
      modq%nqty=nqty
      modq%pd=pd
      modq%pdmin=pdmin
      modq%pdmax=pdmax
      modq%n_int=(pdmax-pdmin+1)*(2*pd+pdmin-pdmax+1)
!-----------------------------------------------------------------------
!     allocate space.
!-----------------------------------------------------------------------
      ALLOCATE(modq%fsi(nqty,modq%n_int,1:mx,1:my))
      ALLOCATE(modq%title(nqty))
      ALLOCATE(modq%f(nqty))
      ALLOCATE(modq%fx(nqty))
      ALLOCATE(modq%fy(nqty))
!-----------------------------------------------------------------------
!     character descriptors, if present in input.
!-----------------------------------------------------------------------
      IF (PRESENT(name)) modq%name=name
      IF (PRESENT(title)) THEN
        IF (SIZE(title)==nqty) THEN
          modq%title=title
        ELSE
          modq%title=title(1)
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     all discontinuous bases are element-centered, so ix0=iy0=1.
!     for modal elements, their are no nodal positions (dx and dy in
!     nimrod's nodal representations).
!-----------------------------------------------------------------------
      ALLOCATE(modq%ix0(modq%n_int))
      ALLOCATE(modq%iy0(modq%n_int))
      modq%ix0=1
      modq%iy0=1
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_alloc
!-----------------------------------------------------------------------
!     subprogram 15. modal_disc_2D_dealloc.
!     deallocates space for modal_disc_2D_type.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_dealloc(modq)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
!-----------------------------------------------------------------------
!     deallocate space.
!-----------------------------------------------------------------------
      DEALLOCATE(modq%fsi)
      DEALLOCATE(modq%title)
      DEALLOCATE(modq%f)
      DEALLOCATE(modq%fx)
      DEALLOCATE(modq%fy)
      DEALLOCATE(modq%ix0)
      DEALLOCATE(modq%iy0)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_dealloc
!-----------------------------------------------------------------------
!     subprogram 16. modal_disc_2D_eval.
!     evaluates complex modal_disc quantities at a single point within a
!     grid block.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_eval(modq,x,y,dmode)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      REAL(r8), INTENT(IN) :: x,y
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: ix,iy,i,j,k,im
      REAL(r8) :: xstd,ystd
      REAL(r8), DIMENSION(0:npoly_max) :: alx,aly,dalx,daly
!-----------------------------------------------------------------------
!     find the interval, and compute 1D basis coefficients with logical
!     coordinates for a standard polynomial in -1<=xstd,ystd<=+1.
!-----------------------------------------------------------------------
      ix=MAX(MIN(INT(x),modq%mx-1),0_i4)
      iy=MAX(MIN(INT(y),modq%my-1),0_i4)
      xstd=2._r8*(x-ix)-1._r8
      ystd=2._r8*(y-iy)-1._r8
      DO i=0,modq%pd
        CALL leg_poly1_sub(xstd,i,alx(i),dalx(i))
        CALL leg_poly1_sub(ystd,i,aly(i),daly(i))
      ENDDO
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.  it is
!     important that the basis-function order is the same as defined
!     in modal_disc_bases.
!-----------------------------------------------------------------------
      k=1
      modq%f=0._r8
      modq%fx=0._r8
      modq%fy=0._r8

      IF (dmode==0) THEN
        DO j=modq%pdmin,modq%pdmax
          DO i=0,modq%pd
            modq%f=modq%f+modq%fsi(:,k,ix+1,iy+1)*alx(i)*aly(j)
            k=k+1
          ENDDO
        ENDDO
        DO i=modq%pdmin,modq%pdmax
          DO j=0,modq%pdmin-1
            modq%f=modq%f+modq%fsi(:,k,ix+1,iy+1)*alx(i)*aly(j)
            k=k+1
          ENDDO
          DO j=modq%pdmax+1,modq%pd
            modq%f=modq%f+modq%fsi(:,k,ix+1,iy+1)*alx(i)*aly(j)
            k=k+1
          ENDDO
        ENDDO
      ELSE
        DO j=modq%pdmin,modq%pdmax
          DO i=0,modq%pd
            modq%f =modq%f +modq%fsi(:,k,ix+1,iy+1)* alx(i)* aly(j)
            modq%fx=modq%fx+modq%fsi(:,k,ix+1,iy+1)*dalx(i)* aly(j)
            modq%fy=modq%fy+modq%fsi(:,k,ix+1,iy+1)* alx(i)*daly(j)
            k=k+1
          ENDDO
        ENDDO
        DO i=modq%pdmin,modq%pdmax
          DO j=0,modq%pdmin-1
            modq%f =modq%f +modq%fsi(:,k,ix+1,iy+1)* alx(i)* aly(j)
            modq%fx=modq%fx+modq%fsi(:,k,ix+1,iy+1)*dalx(i)* aly(j)
            modq%fy=modq%fy+modq%fsi(:,k,ix+1,iy+1)* alx(i)*daly(j)
            k=k+1
          ENDDO
          DO j=modq%pdmax+1,modq%pd
            modq%f =modq%f +modq%fsi(:,k,ix+1,iy+1)* alx(i)* aly(j)
            modq%fx=modq%fx+modq%fsi(:,k,ix+1,iy+1)*dalx(i)* aly(j)
            modq%fy=modq%fy+modq%fsi(:,k,ix+1,iy+1)* alx(i)*daly(j)
            k=k+1
          ENDDO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_eval
!-----------------------------------------------------------------------
!     subprogram 17. modal_disc_2D_all_eval.
!     evaluates complex modal_disc quantities in all elements in a grid
!     block for equal spacing.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_all_eval(modq,x,y,f,fx,fy,dmode)

      TYPE(modal_quad_2D_type), INTENT(IN) :: modq
      REAL(r8), INTENT(IN) :: x,y
      REAL(r8), INTENT(OUT), DIMENSION(:,:,:) :: f,fx,fy
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: i,j,k,ix,iy
      REAL(r8), DIMENSION(modq%n_int) :: alpha,dalpdx,dalpdy
!-----------------------------------------------------------------------
!     evaluate the bases at the specified offset (x,y) within each
!     element.
!-----------------------------------------------------------------------
      CALL modal_disc_bases(x,y,alpha,dalpdx,dalpdy,dmode,              &
     &                      modq%pd,modq%pdmin,modq%pdmax)
!-----------------------------------------------------------------------
!     evaluate the function up to the requested derivative.
!-----------------------------------------------------------------------
      IF (dmode==0) THEN
        DO iy=1,modq%my
          DO ix=1,modq%mx
            DO i=1,modq%nqty
              f(i,ix,iy)=SUM(modq%fsi(i,:,ix,iy)*alpha(:))
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO iy=1,modq%my
          DO ix=1,modq%mx
            DO i=1,modq%nqty
              f(i,ix,iy)=SUM(modq%fsi(i,:,ix,iy)*alpha(:))
              fx(i,ix,iy)=SUM(modq%fsi(i,:,ix,iy)*dalpdx(:))
              fy(i,ix,iy)=SUM(modq%fsi(i,:,ix,iy)*dalpdy(:))
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_all_eval
!-----------------------------------------------------------------------
!     subprogram 18. modal_disc_2D_assign_rsc.
!     assign a real scalar value to a real modalange quad structure.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_assign_rsc(modq,rscalar)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      REAL(r8), INTENT(IN) :: rscalar

      modq%fsi=rscalar
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_assign_rsc
!-----------------------------------------------------------------------
!     subprogram 19. modal_disc_2D_assign_csc.
!     assign a real scalar value to a real modalange quad
!     structure.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_assign_csc(modq,cscalar)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      COMPLEX(r8), INTENT(IN) :: cscalar

      modq%fsi=cscalar
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_assign_csc
!-----------------------------------------------------------------------
!     subprogram 20. modal_disc_2D_assign_modq.
!     set one real modalange quad structure equal to another.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_assign_modq(modq1,modq2)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq1
      TYPE(modal_quad_2D_type), INTENT(IN) :: modq2

!-----------------------------------------------------------------------
!     warning:  there are no checks on compatibility.
!-----------------------------------------------------------------------
      modq1%fsi(:,:,:,:)=modq2%fsi(:,:,:,:)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_assign_modq
!-----------------------------------------------------------------------
!     subprogram 21. modal_disc_2D_assign_int.
!     assign a integer value to a real modalange quad structure.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_assign_int(modq,int)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      INTEGER(i4), INTENT(IN) :: int

      modq%fsi=int
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_assign_int
!-----------------------------------------------------------------------
!     subprogram 22. modal_disc_2D_basis_assign_arr
!     assign data into coefficient arrays for one basis function.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_basis_assign_arr(modq,data,ibasis)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      modq%fsi(:,ibasis,:,:)=data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_basis_assign_arr
!-----------------------------------------------------------------------
!     subprogram 23. modal_disc_2D_basis_add_arr
!     add data into coefficient arrays for one basis function.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_basis_add_arr(modq,data,ibasis)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      modq%fsi(:,ibasis,:,:)=modq%fsi(:,ibasis,:,:)+data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_basis_add_arr
!-----------------------------------------------------------------------
!     subprogram 24. modal_disc_2D_basis_assign_loc
!     assign data into coefficient arrays for one basis function.
!
!     this is a local version of modal_disc_basis_assign_arr, where the
!     the data is located at given poloidal indices, only.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_basis_assign_loc(modq,data,ibasis,ix,iy)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      REAL(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      modq%fsi(:,ibasis,ix,iy)=data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_basis_assign_loc
!-----------------------------------------------------------------------
!     subprogram 25. modal_disc_2D_basis_add_loc
!     add data into coefficient arrays for one basis function.
!
!     this is a local version of modal_disc_basis_add_arr, where the
!     the data is located at given poloidal indices only.
!-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_basis_add_loc(modq,data,ibasis,ix,iy)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      REAL(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy
!-----------------------------------------------------------------------
!     decide the proper storage location for this data according to
!     the value of ibasis.
!-----------------------------------------------------------------------
      modq%fsi(:,ibasis,ix,iy)=modq%fsi(:,ibasis,ix,iy)+data
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_basis_add_loc
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      END MODULE modal_disc_mod
