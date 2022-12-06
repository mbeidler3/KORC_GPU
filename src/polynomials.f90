!-----------------------------------------------------------------------
!     $Id: polynomials.f90 7440 2021-03-19 23:27:32Z held $
!     subprograms for the 1D polynomials used in the finite elements
!     and their numerical integration.  this is not a module.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.
!-----------------------------------------------------------------------
!     1.  poly_mod.
!     2.  poly_set.
!     3.  poly_inquire.
!     4.  poly_nodes.
!     5.  lagr_1D.
!     6.  gauleg.
!     7.  lobleg.
!     8.  radleg.
!     9.  polint.
!     10. gaulag.
!     11. gammln.
!     12. hierarch_poly.
!     13. leg_poly1_sub.
!     14. legendre_poly.
!     15. legendre_polyd.
!     16. legendre_polyd2.
!     17. legendre_eval.
!     18. legendre_deriv.
!     19. legendre_deriv2.
!     20. vp_init.
!     21. vp_init_quads.
!     22. vp_init_w_projection.
!-----------------------------------------------------------------------
!     subprogram 1. poly_mod.
!     holds information that influences the location of basis functions
!     nodes.
!-----------------------------------------------------------------------
      MODULE poly_mod
      USE local
      IMPLICIT NONE

      CHARACTER(7) :: node_dist='uniform'
!     CHARACTER(7) :: node_dist='gll'

!-----------------------------------------------------------------------
!     These save the node locations such that they are not recomputed
!     unless n_last and dist_last are changed.
!-----------------------------------------------------------------------
      INTEGER(i4), PARAMETER :: pd_max=24
      INTEGER(i4), SAVE :: n_last=-1
      CHARACTER(7), SAVE :: dist_last='none'
      REAL(r8), DIMENSION(0:pd_max), SAVE :: x_last
      REAL(r8), DIMENSION(0:pd_max,0:pd_max), SAVE :: cardcoefs

      END MODULE poly_mod

!-----------------------------------------------------------------------
!     subprogram 2. poly_set.
!     provides a mechanism to change subsequent basis node position
!     evaluations during the execution of a program.
!     if this is used, it should be called before any allocation of
!     "lagr_quad" data structures; otherwise, subsequent evaluations
!     will be erroneous!
!-----------------------------------------------------------------------
      SUBROUTINE poly_set(new_dist)
      USE poly_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: new_dist

      node_dist=new_dist

      RETURN
      END SUBROUTINE poly_set

!-----------------------------------------------------------------------
!     subprogram 3. poly_inquire.
!     returns the character variable that determines how nodes are
!     distributed.
!-----------------------------------------------------------------------
      SUBROUTINE poly_inquire(dist)
      USE poly_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(OUT) :: dist

      dist=node_dist

      RETURN
      END SUBROUTINE poly_inquire

!-----------------------------------------------------------------------
!     subprogram 4. poly_nodes.
!     finds the n+1 location of nodes for 1D lagrange polynomials in the
!     domain 0<=x<=1.
!
!     the distribution of nodes is set according to the module variable
!     node_dist.  The value 'uniform' gives uniform nodes, and 'gll'
!     gives a nonuniform distribution.  For the latter,
!     the location of nodes corresponds to the zeros of the
!     Gauss-Lobatto polynomials, (1+x)(1-x) * d L_n(x)/dx, where
!     L_n is the n-th order Legendre polynomial.  [See "Spectral/hp
!     Element Methods for CFD," Karniadakis and Sherwin, for example.]
!-----------------------------------------------------------------------
      SUBROUTINE poly_nodes(n,x)
      USE poly_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), DIMENSION(0:n), INTENT(OUT) :: x

      INTEGER(i4) :: i,j
      REAL(r8) :: wght
      REAL(r8), DIMENSION(0:n) :: wgll
      REAL(r8), EXTERNAL :: legendre_poly
      CHARACTER(128) :: msg

!-----------------------------------------------------------------------
!     don't repeat computation if the nodes haven't changed.
!     with OpenMP this requires a double check. Enter a openMP critical
!     section to update, and then do a check that a previous thread
!     has not already done the update while we are waiting.
!-----------------------------------------------------------------------
      IF (n==n_last.AND.node_dist==dist_last) THEN
        x=x_last(0:n)
      ELSE
        IF (n>pd_max) THEN
          WRITE(msg,'(a,i8,a)') "POLY_NODES: requested ",n,             &
     &      " exceeds pd_max. Increase pd_max and recompile."
          CALL nim_stop(TRIM(msg))
        ENDIF     
        !$omp critical
        IF (n==n_last.AND.node_dist==dist_last) THEN
          x=x_last(0:n) ! a previous thread has performed the update
        ELSE
          n_last=n
          dist_last=node_dist
!-----------------------------------------------------------------------
!         uniform distribution.
!-----------------------------------------------------------------------
          IF (node_dist=='uniform') THEN
            DO i=0,n
              x(i)=REAL(i,r8)/REAL(n,r8)
            ENDDO
          ELSE
!-----------------------------------------------------------------------
!           for the Gauss-Lobatto-Legendre nodes, use the lobleg routine
!           and generate the coefficients of the Legendre polynomial
!           for each cardinal function.  note that cardcoefs is ordered
!           with the expansion index first and the cardinal-function
!           index second.
!-----------------------------------------------------------------------
            CALL lobleg(-1._r8,1._r8,x,wgll,n+1_i4)
            DO j=0,n
              wght=0.5_r8*REAL(2_i4*j+1_i4,r8)
              IF (j==n) wght=0.5_r8*REAL(j,r8)
              DO i=0,n
                cardcoefs(j,i)=legendre_poly(x(i),j)*wgll(i)*wght
              ENDDO
            ENDDO
            x=0.5_r8*x+0.5_r8
          ENDIF
          x_last(0:n)=x
        ENDIF
        !$omp end critical
      ENDIF

      END SUBROUTINE poly_nodes

!-----------------------------------------------------------------------
!     subprogram 5. lagr_1D.
!     finds the coefficients for general 1D lagrange polynomials.
!     note:  the assumed shape arrays al and dal require an interface
!     block in the calling routines.
!-----------------------------------------------------------------------
      SUBROUTINE lagr_1D(pd,x,al,dal,dmode,d2al)
      USE poly_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: pd,dmode
      REAL(r8), INTENT(IN) :: x
      REAL(r8), DIMENSION(0:), INTENT(OUT) :: al,dal
      REAL(r8), DIMENSION(0:), INTENT(OUT), OPTIONAL :: d2al

      INTEGER(i4) :: i,j,k
      REAL(r8), DIMENSION(0:pd) :: c_norm,x_node
      REAL(r8) :: dxtmp

      REAL(r8), EXTERNAL :: legendre_eval,legendre_deriv,legendre_deriv2
!-----------------------------------------------------------------------
!     get the locations of the nodes (zeros of the basis functions).
!-----------------------------------------------------------------------
      CALL poly_nodes(pd,x_node)
!-----------------------------------------------------------------------
!     if the node distribution is GLL, use the coefficients of the
!     Legendre-polynomial expansion of each cardinal function.
!-----------------------------------------------------------------------
      IF (node_dist=='gll') THEN
        DO i=0,pd
          al(i)=legendre_eval(pd,x,0._r8,1._r8,cardcoefs(0:pd,i))
        ENDDO
        IF (dmode==1) THEN
          DO i=0,pd
            dal(i)=legendre_deriv(pd,x,0._r8,1._r8,cardcoefs(0:pd,i))
          ENDDO
        ELSEIF (dmode==2) THEN
          DO i=0,pd
            dal(i)=legendre_deriv(pd,x,0._r8,1._r8,cardcoefs(0:pd,i))
            d2al(i)=legendre_deriv2(pd,x,0._r8,1._r8,cardcoefs(0:pd,i))
          ENDDO
        ENDIF
        RETURN
      ENDIF
!-----------------------------------------------------------------------
!     for other distributions get the normalization constant for the
!     formal cardinal-function relation.
!-----------------------------------------------------------------------
      c_norm=1
      DO i=0,pd
        DO j=0,pd
          IF (j==i) CYCLE
          c_norm(i)=c_norm(i)/(x_node(i)-x_node(j))
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     compute 1D basis values.
!-----------------------------------------------------------------------
      DO i=0,pd
        al(i)=c_norm(i)
        DO j=0,pd
          IF (j==i) CYCLE
          al(i)=al(i)*(x-x_node(j))
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     compute first derivatives.
!-----------------------------------------------------------------------
      IF (dmode<1) RETURN
      DO i=0,pd
        dal(i)=0
        DO k=0,pd
          IF (k==i) CYCLE
          dxtmp=c_norm(i)
          DO j=0,pd
            IF (j==i) CYCLE
            IF (j==k) CYCLE
            dxtmp=dxtmp*(x-x_node(j))
          ENDDO
          dal(i)=dal(i)+dxtmp
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate execution of lagr_1D.
!-----------------------------------------------------------------------
      RETURN

      CONTAINS

!-----------------------------------------------------------------------
!       factorial function.
!-----------------------------------------------------------------------
        FUNCTION lagr_1D_fac(n) RESULT(nfac)

        INTEGER(i4) :: nfac
        INTEGER(i4), INTENT(IN) :: n

        INTEGER(i4) :: jj

        nfac=1
        DO jj=2,n
          nfac=nfac*jj
        ENDDO

        END FUNCTION lagr_1D_fac

      END SUBROUTINE lagr_1D

!-----------------------------------------------------------------------
!     subprogram 6. gauleg.
!     abscissas and weights for Gauss-Legendre integration, adapted
!     from Numerical Recipies, 2nd ed., Cambridge Press.
!-----------------------------------------------------------------------
      SUBROUTINE gauleg(x1,x2,x,w,n)
      USE local
      USE physdat
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: x1,x2
      REAL(r8), DIMENSION(n), INTENT(OUT) :: x,w

      INTEGER(i4) :: i,j,m
      REAL(r8) :: p1,p2,p3,pp,xl,xm,z,z1
      REAL(r8), PARAMETER :: eps=1.e-14

      m=(n+1)/2
      xm=0.5_r8*(x2+x1)
      xl=0.5_r8*(x2-x1)
      DO i=1,m
        z=COS(pi*(REAL(i,r8)-0.25_r8)/(REAL(n,r8)+0.5_r8))
        DO
          p1=1
          p2=0
          DO j=1,n
            p3=p2
            p2=p1
            p1=(REAL(2*j-1,r8)*z*p2-REAL(j-1,r8)*p3)/REAL(j,r8)
          ENDDO
          pp=REAL(n,r8)*(z*p1-p2)/(z**2-1._r8)
          z1=z
          z=z1-p1/pp
          IF (ABS(z-z1)<=eps) EXIT
        ENDDO
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2._r8*xl/((1._r8-z**2)*pp**2)
        w(n+1-i)=w(i)
      ENDDO

      RETURN
      END SUBROUTINE gauleg

!-----------------------------------------------------------------------
!     subprogram 7. lobleg.
!     abscissas and weights for Lobatto-Legendre integration, based on
!     Abromowitz and Stegun.
!
!     this version is now independent of poly_nodes and is used by
!     poly_nodes.  its computation has been refined to get rid of the
!     numerical differentiation.
!-----------------------------------------------------------------------
      SUBROUTINE lobleg(x1,x2,x,w,n)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: x1,x2
      REAL(r8), DIMENSION(n), INTENT(OUT) :: x,w

      INTEGER(i4) :: i,j,m,nleg
      REAL(r8) :: fac,xl,xm,z,z1,p1,p2,p3,gg,dg
      REAL(r8), PARAMETER :: eps=1.e-15
      REAL(r8), DIMENSION(n-1) :: xgau,wgau
!-----------------------------------------------------------------------
!     the zeros of the Legendre polynomial of degree n-1 are used to
!     start Newton's method for finding the zeros of the derivative of
!     the polynomial.
!-----------------------------------------------------------------------
      nleg=n-1_i4
      CALL gauleg(-1._r8,1._r8,xgau,wgau,nleg)
!-----------------------------------------------------------------------
!     the GLL node locations are the zeros of (x**2-1)*d(L_nleg)/dx for
!     the interval -1<=x<=1.  bisect the intervals between zeros of
!     L_nleg to start each Newton iteration and use recurrence to
!     evaluate
!         g(x)=(x**2-1)*d(L_nleg)/dx=nleg*(x*L_nleg-L_(nleg-1))
!     and
!         dg/dx=nleg*(nleg+1)*L_nleg
!-----------------------------------------------------------------------
      m=(n+1_i4)/2_i4
      xm=0.5_r8*(x2+x1)
      xl=0.5_r8*(x2-x1)
      fac=2._r8*xl/REAL(n*(n-1_i4),r8)

      w(1)=fac
      x(1)=x1
      w(n)=fac
      x(n)=x2
      DO i=2,m
        z=0.5_r8*(xgau(i)+xgau(i-1))
        DO
          p1=1._r8
          p2=0._r8
          DO j=1,nleg
            p3=p2
            p2=p1
            p1=(REAL(2*j-1,r8)*z*p2-REAL(j-1,r8)*p3)/REAL(j,r8)
          ENDDO
          gg=REAL(nleg,r8)*(z*p1-p2)
          dg=REAL(nleg*(nleg+1_i4),r8)*p1
          z1=z
          z=z1-gg/dg
          IF (ABS(z-z1)<=eps) EXIT
        ENDDO
        w(i)=fac/p1**2
        w(n+1-i)=w(i)
        x(i)=xm+xl*z
        x(n+1-i)=xm-xl*z
      ENDDO

      RETURN
      END SUBROUTINE lobleg

!-----------------------------------------------------------------------
!     subprogram 8. radleg.
!     abscissas and weights for Radau-Legendre integration, based on
!     Abromowitz and Stegun.
!-----------------------------------------------------------------------
      SUBROUTINE radleg(x1,x2,x,w,n)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: x1,x2
      REAL(r8), DIMENSION(n), INTENT(OUT) :: x,w

      INTEGER(i4) :: i,j,m,nleg
      REAL(r8) :: fac,xl,xm,z,z1,gg,dg,hh,dh,ff,df
      REAL(r8), PARAMETER :: eps=1.e-15
      REAL(r8), DIMENSION(n-1) :: xgau,wgau
      REAL(r8), DIMENSION(n) :: xgaun,wgaun
!-----------------------------------------------------------------------
!     the zeros of the Legendre polynomials of degree n-1 and n are used
!     to start Newton's method for finding the zeros of their sums.
!-----------------------------------------------------------------------
      nleg=n-1_i4
      CALL gauleg(-1._r8,1._r8,xgau,wgau,nleg)
      CALL gauleg(-1._r8,1._r8,xgaun,wgaun,n)
!-----------------------------------------------------------------------
!     the GRL node locations are the zeros of the polynomial
!
!      (L_nleg+L_(nleg+1))/(x+1) = L_nleg+(x-1)*[d(L_nleg)/dx]/(nleg+1)
!
!     the interval -1<=x<=1.  Use 1st-order Taylor expansion about zeros
!     of the two Legendre polynomials on the lhs for an initial guess.
!     to start each Newton iteration.  Use recurrence to
!     evaluate
!         g(x)=L_nleg+(x-1)*[d(L_nleg)/dx]/(nleg+1)
!     and
!         dg/dx={nleg*L_nleg+
!                [(1+x)*(nleg+1)+1-x]*[d(L_nleg)/dx]/(nleg+1)}/(1+x)
!-----------------------------------------------------------------------
      m=(n+1_i4)/2_i4
      xm=0.5_r8*(x2+x1)
      xl=0.5_r8*(x2-x1)
      fac=xl/REAL(n**2,r8)

      w(1)=2._r8*fac
      x(1)=x1
      DO i=2,n
        CALL leg_poly1_sub(xgau(i-1),nleg,gg,dg)
        CALL leg_poly1_sub(xgaun(i),n,ff,df)
        z=(xgau(i-1)*dg+xgaun(i)*df)/(dg+df)
        DO
          CALL leg_poly1_sub(z,nleg,ff,df)
          gg=ff+(z-1._r8)*df/REAL(n,r8)
          dg=(REAL(n,r8)*ff+(z+1._r8)*df-gg)/(z+1._r8)
          z1=z
          z=z1-gg/dg
          IF (ABS(z-z1)<=eps) EXIT
        ENDDO
        w(i)=fac*(1._r8-z)/ff**2
        x(i)=xm+xl*z
      ENDDO

      RETURN
      END SUBROUTINE radleg

!-----------------------------------------------------------------------
!     subprogram 9. polint
!     Perform polynomial interpolation given various segments.
!-----------------------------------------------------------------------
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      USE local
      implicit none
      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: x
      REAL(r8), INTENT(OUT) :: dy,y
      REAL(r8), DIMENSION(n), INTENT(IN) :: xa,ya
      INTEGER(i4), PARAMETER :: NMAX=20
      INTEGER(i4) i,m,ns
      REAL(r8) ::  den,dif,dift,ho,hp,w
      REAL(r8), DIMENSION(:), ALLOCATABLE ::  c,d

      ALLOCATE(c(n),d(n))
      c=ya; d=ya
      ns=1
      dif=abs(x-xa(1))
      DO i=1,n
        dift=abs(x-xa(i))
        IF (dift<dif) THEN
          ns=i
          dif=dift
        ENDIF
      ENDDO
      y=ya(ns)
      ns=ns-1
      DO m=1,n-1
        DO i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den==0.) CALL nim_stop('failure in polint')
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
        ENDDO
        IF (2*ns<n-m)THEN
          dy=c(ns+1)
        ELSE
          dy=d(ns)
          ns=ns-1
        ENDIF
        y=y+dy
      ENDDO
      DEALLOCATE(c,d)

      RETURN
      END SUBROUTINE polint

!-----------------------------------------------------------------------
!     subprogram 10. gaulag.
!     abscissas and weights for Gauss-Laguerre integration, adapted
!     from Numerical Recipies, 2nd ed., Cambridge Press.
!-----------------------------------------------------------------------
      SUBROUTINE gaulag(x,w,n,alf)
      USE local
      implicit none
      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), DIMENSION(n), INTENT(INOUT) :: w(n),x(n)
      REAL(r8), INTENT(IN) :: alf
!U    USES gammln
      INTEGER(i4), PARAMETER :: MAXIT=10
      REAL(r8), PARAMETER :: EPS=3.e-14
      INTEGER(i4) :: i,its,j
      REAL(r8) :: ai,gammln
      REAL(r8) :: p1,p2,p3,pp,z,z1
      do 13 i=1,n
        if(i==1)then
          z=(1.+alf)*(3.+.92*alf)/(1.+2.4*n+1.8*alf)
        else if(i==2)then
          z=z+(15.+6.25*alf)/(1.+.9*alf+2.5*n)
        else
          ai=i-2
          z=z+((1.+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.+3.5*ai))*          &
     &(z-x(i-2))/(1.+.3*alf)
        endif
        do 12 its=1,MAXIT
          p1=1.
          p2=0.
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
   11     continue
          pp=(n*p1-(n+alf)*p2)/z
          z1=z
          z=z1-p1/pp
          if(abs(z-z1)<=EPS)goto 1
   12   continue
!        pause 'too many iterations in gaulag'
        WRITE(*,*) 'too many iterations in gaulag'
    1   x(i)=z
        w(i)=-exp(gammln(alf+n)-gammln(1._r8*n))/(pp*n*p2)
   13 continue
      return
      END SUBROUTINE gaulag

!-----------------------------------------------------------------------
!     subprogram 11. gammln.
!     used in conjunction with gaulag above and
!     adapted from Numerical Recipies, 2nd ed., Cambridge Press.
!-----------------------------------------------------------------------
      FUNCTION gammln(xx)
      USE local
      implicit none
      REAL(r8) :: gammln,xx
      INTEGER(i4) :: j
      REAL ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,            &
     &24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,    &
     &-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
   11 continue
      gammln=tmp+log(stp*ser/x)
      return
      END

!-----------------------------------------------------------------------
!     subprogram 12. hierarch_poly.
!     compute value and derivative of n-th order hierarchical
!     polynomial defined on domain -1 <= x <= 1 as:
!     f_0 = (1 + x)/2
!     f_pd_xi = (1 - x)/2
!     f_n = (P_n - P_n-2)/SQRT(2(2n-1)) for 1 < n < pd_xi
!     where P_n(x) are Legendre polynomials.
!     form used here is from Commun. Numer. Meth. Engng 2007; 00:1-13
!     by Sprague and Geers.
!-----------------------------------------------------------------------
      SUBROUTINE hierarch_poly(pd_xi,x,h,dh)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: pd_xi
      REAL(r8), INTENT(OUT) :: h(0:pd_xi),dh(0:pd_xi)
      REAL(r8), INTENT(IN) :: x

      INTEGER(i4) :: n
      REAL(r8) :: p_n,p_nm1,p_nm2
      REAL(r8), EXTERNAL :: legendre_poly

      h(0) = (1._r8+x)/2._r8
      dh(0) =  1._r8/2._r8
      h(pd_xi) = (1._r8-x)/2._r8
      dh(pd_xi)= -1._r8/2._r8
      DO n=2,pd_xi
        p_n   = legendre_poly(x,n)
        p_nm1 = legendre_poly(x,n-1_i4)
        p_nm2 = legendre_poly(x,n-2_i4)
        h(n-1) = (p_n - p_nm2)/SQRT(2._r8*(2._r8*n-1._r8))
        dh(n-1) = p_nm1*SQRT(n-.5_r8)
      ENDDO

      END SUBROUTINE hierarch_poly

!-----------------------------------------------------------------------
!     subprogram 13. leg_poly1_sub.
!
!     This subroutine evaluates a single Legendre polynomial and its
!     derivative, using recurrence.  The coding for the value is
!     extracted from gauleg.
!
!     The argument list is:
!
!     xx [real] {input} -- The value of the independent variable in
!        the standard -1 <= xx <= +1 range.
!     nn [integer] {input} -- The index of the polynomial.
!     val [real] {output} -- The value of L_n(xx).
!     drv [real] {output} -- The derivative of L_n(xx).
!-----------------------------------------------------------------------
      SUBROUTINE leg_poly1_sub(xx,nn,val,drv)
      USE local
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: xx
      INTEGER(i4), INTENT(IN) :: nn
      REAL(r8), INTENT(OUT) :: val,drv

      INTEGER(i4) :: j
      REAL(r8) :: p2,p3

      val=1._r8
      drv=0._r8
      p2=0._r8
      DO j=1,nn
        p3=p2
        p2=val
        val=(REAL(2*j-1,r8)*xx*p2-REAL(j-1,r8)*p3)/REAL(j,r8)
        drv=xx*drv+REAL(j,r8)*p2
      ENDDO

      RETURN
      END SUBROUTINE leg_poly1_sub

!-----------------------------------------------------------------------
!     subprogram 14. legendre_poly.
!     compute the value of a legendre polynomial of non-negative
!     integer order at a specified position.
!
!     note that the abscissa is expected to be in the standard range
!     of -1 <= x <= +1.
!-----------------------------------------------------------------------
      FUNCTION legendre_poly(x,n) RESULT(ln)
      USE local
      IMPLICIT NONE

      REAL(r8) :: ln
      REAL(r8), INTENT(IN) :: x
      INTEGER(i4), INTENT(IN) :: n

      REAL(r8) :: lim1,lim2
      INTEGER(i4) :: i
!-----------------------------------------------------------------------
!     use standard recursion to generate the value of the desired
!     legendre polynomial.  [Schaum's outline, Mathematical Handbook,
!     by M. Spiegel, for example.]
!-----------------------------------------------------------------------
      IF (n==0) THEN
        ln=1._r8
        RETURN
      ELSE IF (n==1) THEN
        ln=x
        RETURN
      ELSE
        lim2=1._r8
        lim1=x
        DO i=2,n
          ln=(REAL(2*i-1,r8)*x*lim1-REAL(i-1,r8)*lim2)/REAL(i,r8)
          lim2=lim1
          lim1=ln
        ENDDO
      ENDIF

      END FUNCTION legendre_poly

!-----------------------------------------------------------------------
!     subprogram 15. legendre_polyd.
!
!     This is a function-version of leg_poly1_sub that returns just the
!     derivative of a Legendre polynomial of non-negative integer order.
!
!     The argument list is:
!
!     xx [real] {input} -- The value of the independent variable in
!        the standard -1 <= xx <= +1 range.
!     nn [integer] {input} -- The index of the polynomial.
!
!     The result of the function is:
!
!     drv [real]  -- The derivative of L_n(xx).
!-----------------------------------------------------------------------
      FUNCTION legendre_polyd(xx,nn) RESULT(drv)
      USE local
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: xx
      INTEGER(i4), INTENT(IN) :: nn
      REAL(r8) :: drv

      INTEGER(i4) :: j
      REAL(r8) :: p2,p3,val

      val=1._r8
      drv=0._r8
      p2=0._r8
      DO j=1,nn
        p3=p2
        p2=val
        val=(REAL(2*j-1,r8)*xx*p2-REAL(j-1,r8)*p3)/REAL(j,r8)
        drv=xx*drv+REAL(j,r8)*p2
      ENDDO

      RETURN
      END FUNCTION legendre_polyd

!-----------------------------------------------------------------------
!     subprogram 16. legendre_polyd2.
!
!     This function returns just the derivative of a Legendre polynomial
!     of non-negative integer order.
!
!     The argument list is:
!
!     x [real] {input} -- The value of the independent variable in
!        the standard -1 <= xx <= +1 range.
!     n [integer] {input} -- The index of the polynomial.
!
!     The result of the function is:
!
!     drv2 [real]  -- The second derivative of L_n(xx).
!-----------------------------------------------------------------------
      FUNCTION legendre_polyd2(x,n) RESULT(drv2)
      USE local
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: x
      INTEGER(i4), INTENT(IN) :: n
      REAL(r8) :: drv2

      REAL(r8), EXTERNAL :: legendre_polyd
      
      integer(i4) :: j
      REAL(r8) :: lp, lpp1, lpp2
  
      drv2 = 0._r8
      lpp1 = 0._r8
      lpp2 = 0._r8

      DO j=2,n
        lp = legendre_polyd(x,j-1)
        lpp2 = lpp1
        lpp1 = drv2
        drv2 = (REAL(2*j-1,r8)*(2._r8*lp+x*lpp1)-(REAL(j-1,r8))*lpp2) / &
   &            REAL(j, r8)
      ENDDO
  
      RETURN
      END FUNCTION legendre_polyd2

!-----------------------------------------------------------------------
!     subprogram 17. legendre_eval.
!
!     The function legendre_eval returns the value of a
!     Legendre series approximation with nmax+1 terms.  The routine
!     needs the nmax+1 expansion coefficients, and the independent
!     variable y in the domain ymin <= y <= ymax, which gets mapped
!     to the standard -1 <= x <= 1.
!
!     The argument list for the subroutine is:
!
!     nmax [integer] {input} - Maximum index for the expansion with
!	   (nmax+1) terms for indices 0 through nmax.
!     yy   [real] {input} - The value of the independent variable
!          for the evaluation.
!     ymin [real] {input} - The minimum value of the independent
!          variable for functions on an arbitrary domain.
!     ymax [real] {input} - The maximum value of the independent
!          variable for functions on an arbitrary domain.
!     coef [real(0:nmax)] {input} - Holds the coefficients of the
!	   Legendre series.
!
!     The real legendre_eval function then evaluates to the sum
!
!	   f(x) = sum_n{ c_n*L_n(x), 0<=n<=nmax }
!
!     where x = x(y) is a linear function of yy.
!
!     This routine is adapted from cheb_eval from cyl_spec.
!-----------------------------------------------------------------------

      FUNCTION legendre_eval(nmax,yy,ymin,ymax,coef) RESULT(approx)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nmax
      REAL(r8), INTENT(IN) :: yy,ymin,ymax
      REAL(r8), DIMENSION(0:nmax), INTENT(IN) :: coef
      REAL(r8) :: approx

      REAL(r8), EXTERNAL :: legendre_poly
!-----------------------------------------------------------------------
!     Local variables:
!-----------------------------------------------------------------------

      REAL(r8) :: xx
      INTEGER(i4) :: mm

!-----------------------------------------------------------------------
!     Evaluate the argument for the Legendre polynomials.
!-----------------------------------------------------------------------

      xx=(2._r8*yy-ymin-ymax)/(ymax-ymin)

!-----------------------------------------------------------------------
!     Use the legendre_poly function for each term in the series.
!-----------------------------------------------------------------------

      approx=coef(0)
      IF (nmax>0) approx=approx+coef(1)*xx

      DO mm=2,nmax
        approx=approx+coef(mm)*legendre_poly(xx,mm)
      ENDDO

      RETURN
      END FUNCTION legendre_eval

!-----------------------------------------------------------------------
!     subprogram 18. legendre_deriv.
!
!     The function legendre_deriv returns the derivative of a
!     Legendre series approximation with nmax+1 terms.  The routine
!     needs the nmax+1 expansion coefficients, and the independent
!     variable y in the domain ymin <= y <= ymax, which gets mapped
!     to the standard -1 <= x <= 1.
!
!     The argument list for the subroutine is:
!
!     nmax [integer] {input} - Maximum index for the expansion with
!	   (nmax+1) terms for indices 0 through nmax.
!     yy   [real] {input} - The value of the independent variable
!          for the evaluation.
!     ymin [real] {input} - The minimum value of the independent
!          variable for functions on an arbitrary domain.
!     ymax [real] {input} - The maximum value of the independent
!          variable for functions on an arbitrary domain.
!     coef [real(0:nmax)] {input} - Holds the coefficients of the
!	   Legendre series.
!
!     The real legendre_eval function then evaluates to the sum
!
!	   f(x) = sum_n{ c_n*L_n(x), 0<=n<=nmax }
!
!     where x = x(y) is a linear function of yy.
!
!     This routine is adapted from cheb_eval from cyl_spec.
!-----------------------------------------------------------------------

      FUNCTION legendre_deriv(nmax,yy,ymin,ymax,coef) RESULT(approx)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nmax
      REAL(r8), INTENT(IN) :: yy,ymin,ymax
      REAL(r8), DIMENSION(0:nmax), INTENT(IN) :: coef
      REAL(r8) :: approx

      REAL(r8), EXTERNAL :: legendre_polyd
!-----------------------------------------------------------------------
!     Local variables:
!-----------------------------------------------------------------------

      REAL(r8) :: xx
      INTEGER(i4) :: mm

!-----------------------------------------------------------------------
!     Evaluate the argument for the Legendre polynomials.
!-----------------------------------------------------------------------

      xx=(2._r8*yy-ymin-ymax)/(ymax-ymin)

!-----------------------------------------------------------------------
!     Use the legendre_polyd function for each term in the series.
!-----------------------------------------------------------------------

      approx=0._r8
      IF (nmax>0) approx=approx+coef(1)

      DO mm=2,nmax
        approx=approx+coef(mm)*legendre_polyd(xx,mm)
      ENDDO
      approx=approx*2._r8/(ymax-ymin)

      RETURN
      END FUNCTION legendre_deriv

!-----------------------------------------------------------------------
!     subprogram 19. legendre_deriv2.
!
!     The function legendre_deriv2 returns the second derivative of a
!     Legendre series approximation with nmax+1 terms.  The routine
!     needs the nmax+1 expansion coefficients, and the independent
!     variable y in the domain ymin <= y <= ymax, which gets mapped
!     to the standard -1 <= x <= 1.
!
!     The argument list for the subroutine is:
!
!     nmax [integer] {input} - Maximum index for the expansion with
!	   (nmax+1) terms for indices 0 through nmax.
!     yy   [real] {input} - The value of the independent variable
!          for the evaluation.
!     ymin [real] {input} - The minimum value of the independent
!          variable for functions on an arbitrary domain.
!     ymax [real] {input} - The maximum value of the independent
!          variable for functions on an arbitrary domain.
!     coef [real(0:nmax)] {input} - Holds the coefficients of the
!	   Legendre series.
!
!     The real legendre_eval function then evaluates to the sum
!
!	   f(x) = sum_n{ c_n*L_n(x), 0<=n<=nmax }
!
!     where x = x(y) is a linear function of yy.
!
!     This routine is adapted from cheb_eval from cyl_spec.
!-----------------------------------------------------------------------

      FUNCTION legendre_deriv2(nmax,yy,ymin,ymax,coef) RESULT(approx)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nmax
      REAL(r8), INTENT(IN) :: yy,ymin,ymax
      REAL(r8), DIMENSION(0:nmax), INTENT(IN) :: coef
      REAL(r8) :: approx

      REAL(r8), EXTERNAL :: legendre_polyd2
!-----------------------------------------------------------------------
!     Local variables:
!-----------------------------------------------------------------------

      REAL(r8) :: xx
      INTEGER(i4) :: mm

!-----------------------------------------------------------------------
!     Evaluate the argument for the Legendre polynomials.
!-----------------------------------------------------------------------

      xx=(2._r8*yy-ymin-ymax)/(ymax-ymin)

!-----------------------------------------------------------------------
!     Use the legendre_polyd function for each term in the series.
!-----------------------------------------------------------------------

      approx=0._r8
      IF (nmax>1) approx=approx+3._r8*coef(2)

      DO mm=3,nmax
        approx=approx+coef(mm)*legendre_polyd2(xx,mm)
      ENDDO
      approx=approx*(2._r8/(ymax-ymin))**2

      RETURN
      END FUNCTION legendre_deriv2
!-----------------------------------------------------------------------
!     subprogram 20. vp_init.                                 
!     generate grids related to pitch-angle: th=ACOS(xi), xi=vpll/v
!     works for the nodes and the quadrature points
!-----------------------------------------------------------------------
      SUBROUTINE vp_init(x,m_xi,nq_xi,basis,vp,ja,tpb,m_xit,m_xip,qps)
      USE local
      USE physdat, ONLY: pi
      IMPLICIT NONE 
                                                                        
      REAL(r8), DIMENSION(nq_xi), INTENT(IN) :: x
      INTEGER(i4), INTENT(IN) :: m_xi,nq_xi
      CHARACTER(3), INTENT(IN) :: basis
      REAL(r8), DIMENSION(m_xi,nq_xi), INTENT(OUT) :: vp,ja
      REAL(r8), INTENT(IN) :: tpb
      INTEGER(i4), INTENT(IN) :: m_xit,m_xip
      LOGICAL, INTENT(IN) :: qps

      INTEGER(i4) :: i,ig,iy,ip,m_xi_pp
      REAL(r8) :: th,th_tp,tht,xi,eta,thp
!-----------------------------------------------------------------------
!     set up velocity space grid.
!-----------------------------------------------------------------------
      IF (basis(1:3)=='Leg') THEN
!-----------------------------------------------------------------------
!       xi and jacobian at quadrature points for Legendre polynomials
!       which use theta=th=-ACOS(xi), dth_dxi=1./SQRT(1.-xi**2)
!-----------------------------------------------------------------------
        ja(1,:)= 1._r8
        IF (qps) THEN
          vp(1,:)= x               ! xi=vp/v used in quadrature
        ELSE
          DO iy=1,m_xi
            DO ig=1,nq_xi
              th=pi*(iy-1+.5_r8*(1._r8+x(ig)))/m_xi
              vp(iy,ig) = -cos(th) ! uniform grid in th for plotting
            ENDDO
          ENDDO
        ENDIF
      ELSE
!-----------------------------------------------------------------------
!       semi-circular vpll/v grid w/ nodes at trapped/passing boundaries
!       or nodes at user specified boundary by passing tpb = xi_node_bnd
!
!       3 domains in pitch-angle with either
!       gridshape_v(1:7) = 'circuni', xi_t=xi_node_bnd and
!       m_xip nodes in:                      -1 < xi < xi_node_bnd
!       m_xit nodes in:            -xi_node_bnd < xi < xi_node_bnd 
!       m_xi-m_xip-m_xit nodes in:  xi_node_bnd < xi < 1, or
!
!       gridshape_v(1:7) = 'circtpb', xi_t=xi_trapped and
!       m_xip nodes in:                      -1 < xi < xi_trapped
!       m_xit nodes in:            -xi_trapped  < xi < xi_trapped
!       m_xi-m_xip-m_xit nodes in:  xi_trapped  < xi < 1
!-----------------------------------------------------------------------
        th_tp = acos(tpb) 
        i = 0 
!-----------------------------------------------------------------------
!       negative passing space: -1 < xi < -xt_t: 
!       xi=-cos(th_tp*eta/m_xip)
!       0 < eta < m_xip
!-----------------------------------------------------------------------
        DO iy=0,m_xip-1 
          i = i+1
          DO ig=1,nq_xi
            eta=iy+x(ig)
            th=th_tp*eta/m_xip
            vp(i,ig) = -cos(th)
            ja(i,ig) = th_tp*sin(th)/m_xip
          ENDDO
        ENDDO
!-----------------------------------------------------------------------
!       trapped space: -xi_t < xi < xt_t 
!       xi=-cos(th_tp + tht*(eta-m_xip)/m_xit), tht=pi-2*th_tp
!       m_xip < eta < m_xip + m_xit 
!-----------------------------------------------------------------------
        tht=pi-2.*th_tp
        DO iy=m_xip,m_xip+m_xit-1 
          i = i+1
          DO ig=1,nq_xi
            eta=iy+x(ig)
            th=th_tp + tht*(eta-m_xip)/m_xit
            vp(i,ig) = -cos(th)
            ja(i,ig) = tht*sin(th)/m_xit
          ENDDO
        ENDDO
!-----------------------------------------------------------------------
!       positive passing space: xi_t < xi < 1  
!       xi = cos(th_tp + tht + th_tp*(eta-m_xip-m_xit)/m_xi_pp)
!       m_xip+m_xit < eta < m_xi 
!-----------------------------------------------------------------------
        m_xi_pp = m_xi - m_xip - m_xit ! # of cells in + passing space
        DO iy=m_xip+m_xit,m_xi-1 
          i = i+1
          DO ig=1,nq_xi
            eta=iy+x(ig)
            th=th_tp + tht + th_tp*(eta-m_xip-m_xit)/m_xi_pp
            vp(i,ig) = -cos(th)
            ja(i,ig) = th_tp*sin(th)/m_xi_pp
          ENDDO
        ENDDO
      ENDIF 
      RETURN
      END SUBROUTINE vp_init
!-----------------------------------------------------------------------
!     subprogram 21. vp_init_quads.                                 
!     generates additional pitch-angle grid information for computing
!     integrals with various pitch-angle dependence.
!-----------------------------------------------------------------------
      SUBROUTINE vp_init_quads(x,m_xi,nq_xi,basis,vp,ja,                &
     &                             tpb,m_xit,m_xip,dth_deta,deta_dth_tp)
      USE local
      USE physdat, ONLY: pi
      IMPLICIT NONE 
                                                                        
      REAL(r8), DIMENSION(nq_xi), INTENT(IN) :: x
      INTEGER(i4), INTENT(IN) :: m_xi,nq_xi
      CHARACTER(3), INTENT(IN) :: basis
      REAL(r8), DIMENSION(m_xi,nq_xi), INTENT(OUT) :: vp,ja
      REAL(r8), INTENT(IN) :: tpb
      INTEGER(i4), INTENT(IN) :: m_xit,m_xip
      REAL(r8), DIMENSION(m_xi,nq_xi), INTENT(OUT)::dth_deta,deta_dth_tp
      !REAL(r8), INTENT(IN), OPTIONAL :: th_x,th_y,Rx,Ry,Zx,Zy

      INTEGER(i4) :: i,ig,iy,ip,m_xi_pp
      REAL(r8) :: th,th_tp,tht,xi,eta
!-----------------------------------------------------------------------
!     set up velocity space grid.
!-----------------------------------------------------------------------
      CALL vp_init(x,m_xi,nq_xi,basis,vp,ja,tpb,m_xit,m_xip,.true.)
!-----------------------------------------------------------------------
!     construct additional terms.
!-----------------------------------------------------------------------
      IF (basis(1:3)=='Leg') THEN
!-----------------------------------------------------------------------
!       Legendre polynomials th=-ACOS(xi), dth_dxi=1./SQRT(1.-xi**2) 
!       where  xi=vpll/v in [-1,1].
!       put dth_dxi in dth_deta even if it is a misnomer.
!-----------------------------------------------------------------------
        deta_dth_tp=0._r8  ! no tpb dependence for vertex nodes
        dth_deta(1,:)=1._r8/SQRT(1._r8-x(:)**2)
      ELSE
!-----------------------------------------------------------------------
!       semi-circular vpll/v grid w/ nodes at trapped/passing boundaries
!       or nodes at user specified boundary by passing tpb = xi_node_bnd
!
!       3 domains in pitch-angle with either
!       gridshape_v(1:7) = 'circuni', xi_t=xi_node_bnd and
!       m_xip nodes in:                      -1 < xi < xi_node_bnd
!       m_xit nodes in:            -xi_node_bnd < xi < xi_node_bnd 
!       m_xi-m_xip-m_xit nodes in:  xi_node_bnd < xi < 1, or
!
!       gridshape_v(1:7) = 'circtpb', xi_t=xi_trapped and
!       m_xip nodes in:                      -1 < xi < xi_trapped
!       m_xit nodes in:            -xi_trapped  < xi < xi_trapped
!       m_xi-m_xip-m_xit nodes in:  xi_trapped  < xi < 1
!-----------------------------------------------------------------------
        th_tp = acos(tpb) 
        i = 0 
!-----------------------------------------------------------------------
!       negative passing space: -1 < xi < -xt_t: 
!       xi=-cos(th_tp*eta/m_xip)
!       0 < eta < m_xip
!-----------------------------------------------------------------------
        DO iy=0,m_xip-1 
          i = i+1
          DO ig=1,nq_xi
            eta=iy+x(ig)
            th=th_tp*eta/m_xip
            dth_deta(i,ig) = th_tp/m_xip
            deta_dth_tp(i,ig) = -m_xip*th/th_tp**2
          ENDDO
        ENDDO
!-----------------------------------------------------------------------
!       trapped space: -xi_t < xi < xt_t 
!       xi=-cos(th_tp + tht*(eta-m_xip)/m_xit), tht=pi-2*th_tp
!       m_xip < eta < m_xip + m_xit 
!-----------------------------------------------------------------------
        tht=pi-2.*th_tp
        DO iy=m_xip,m_xip+m_xit-1 
          i = i+1
          DO ig=1,nq_xi
            eta=iy+x(ig)
            th=th_tp + tht*(eta-m_xip)/m_xit
            dth_deta(i,ig) = tht/m_xit
            deta_dth_tp(i,ig) = (2.*th - pi)*m_xit/tht**2
          ENDDO
        ENDDO
!-----------------------------------------------------------------------
!       positive passing space: xi_t < xi < 1  
!       xi = cos(th_tp + tht + th_tp*(eta-m_xip-m_xit)/m_xi_pp)
!       m_xip+m_xit < eta < m_xi 
!-----------------------------------------------------------------------
        m_xi_pp = m_xi - m_xip - m_xit ! # of cells in + passing space
        DO iy=m_xip+m_xit,m_xi-1 
          i = i+1
          DO ig=1,nq_xi
            eta=iy+x(ig)
            th=th_tp + tht + th_tp*(eta-m_xip-m_xit)/m_xi_pp
            dth_deta(i,ig) = th_tp/m_xi_pp
            deta_dth_tp(i,ig) =-m_xi_pp*(th-pi)/th_tp**2
          ENDDO
        ENDDO
      ENDIF 
      RETURN
      END SUBROUTINE vp_init_quads
!-----------------------------------------------------------------------
!     subprogram 22. vp_init_w_projection.                                 
!     generates pitch-angle-like grid for numerical quadrature and
!     projects basis functions onto Legendre polynomials (P_Leg) and 
!     powers of xi (P_xi).
!-----------------------------------------------------------------------
      SUBROUTINE vp_init_w_projection(x,w,m_xi,nq_xi,basis,vp,ja,       &
     &                   tpb,m_xit,m_xip,dof,lmax,pd_xi,qLag,P_Leg,P_xi)
      USE local
      USE physdat, ONLY: pi
      IMPLICIT NONE 
                                                                        
      REAL(r8), DIMENSION(nq_xi), INTENT(IN) :: x,w
      INTEGER(i4), INTENT(IN) :: m_xi,nq_xi
      CHARACTER(3), INTENT(IN) :: basis
      REAL(r8), DIMENSION(m_xi,nq_xi), INTENT(OUT) :: vp,ja
      REAL(r8), INTENT(IN) :: tpb
      INTEGER(i4), INTENT(IN) :: m_xit,m_xip
      INTEGER(i4), INTENT(IN) :: dof,lmax,pd_xi
      REAL(r8), DIMENSION(nq_xi,0:pd_xi), INTENT(IN) :: qLag
      REAL(r8), INTENT(OUT) :: P_Leg(1:dof,0:lmax),P_xi(1:dof,0:lmax)

      INTEGER(i4) :: i,ig,iy,ip,il,l,n
      REAL(r8) :: th,th_tp,tht,xi,eta,wjac
      REAL(r8), EXTERNAL :: legendre_poly
!-----------------------------------------------------------------------
!     set up pitch-angle-like grid.
!-----------------------------------------------------------------------
      CALL vp_init(x,m_xi,nq_xi,basis,vp,ja,tpb,m_xit,m_xip,.true.)
!-----------------------------------------------------------------------
!     project onto Legendre polynomials and powers of xi
!-----------------------------------------------------------------------
      P_Leg = 0. ; P_xi = 0.
      DO iy=1,m_xi
       il = 1+pd_xi*(iy-1)
       DO l=0,lmax
        DO n=0,pd_xi
         DO ig=1,SIZE(w)
           xi = vp(iy,ig)
           wjac = ja(iy,ig)*w(ig)
           P_Leg(il+n,l) = P_Leg(il+n,l)+legendre_poly(xi,l)            &
     &                                        *wjac*qLag(ig,n)
           P_xi(il+n,l) = P_xi(il+n,l)  +xi**l*wjac*qLag(ig,n)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      RETURN
      END SUBROUTINE vp_init_w_projection
