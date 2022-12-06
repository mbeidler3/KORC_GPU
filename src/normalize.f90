!-----------------------------------------------------------------------
!     $Id: normalize.f90 5564 2018-01-03 20:59:37Z jking $
!     routines to normalize/unnormalize variables before and after
!     dump reads and writes.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.                                                
!-----------------------------------------------------------------------
!     1. normalize              
!     2. unnormalize
!-----------------------------------------------------------------------
      MODULE normalize_mod
      USE local, ONLY: i4,r8
      IMPLICIT NONE

      REAL(r8) :: b0,va,p0,n0,t0,ta,j0
      LOGICAL :: set_norms=.FALSE.
      LOGICAL :: applied_norms=.FALSE.

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 1. normalize_fields
!     apply the normlizations such that
!     n -> n/n0
!     v -> v/vA
!     B -> B/B0
!     p -> p/p0
!     T -> T/T0
!       where vA=B0/SQRT(n0*mu*mi); p0=n0*kb*T0
!     pll_diff -> pll_diff/tau_A (where pll_diff = par_visc and k_pll)
!     perp_diff -> perp_diff*tau_A (perp_diff = iso_visc, k_prp etc)
!     elecd -> elecd/vA L
!       where tau_A = L / vA and L=1m
!     mi -> 1
!     mu0 -> 1
!     e -> L/di (di = SQRT(mi/n0*mu0*e^2)
!     kb -> 1
!     me -> me/mi
!     
!     This transforms all times to t/tA. Also, L=1m is assumed.
!-----------------------------------------------------------------------
      SUBROUTINE normalize_fields(dump_time)
      USE fields, ONLY: nrbl,rb
      IMPLICIT NONE
      REAL(r8), INTENT(INOUT) :: dump_time
      INTEGER(i4) :: ibl

      IF (.NOT.set_norms) CALL set_normalizations
      dump_time=dump_time/ta
!-----------------------------------------------------------------------
!     normalize fields
!-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL norm_field_2D(rb(ibl)%be_eq,1._r8/b0)
        CALL norm_field(rb(ibl)%be,1._r8/b0)
        CALL norm_field_2D(rb(ibl)%ja_eq,1._r8/j0)
        CALL norm_field_2D(rb(ibl)%ve_eq,1._r8/va)
        CALL norm_field(rb(ibl)%ve,1._r8/va)
        CALL norm_field_2D(rb(ibl)%pres_eq,1._r8/p0)
        CALL norm_field(rb(ibl)%pres,1._r8/p0)
        CALL norm_field_2D(rb(ibl)%prese_eq,1._r8/p0)
        CALL norm_field(rb(ibl)%prese,1._r8/p0)
        CALL norm_field_2D(rb(ibl)%nd_eq,1._r8/n0)
        CALL norm_field(rb(ibl)%nd,1._r8/n0)
        CALL norm_field(rb(ibl)%tele,1._r8/t0)
        CALL norm_field(rb(ibl)%tion,1._r8/t0)
      ENDDO
      applied_norms=.TRUE.
      END SUBROUTINE normalize_fields
!-----------------------------------------------------------------------
!     subprogram 2. unnormalize_fields
!     reverse the field normlizations
!-----------------------------------------------------------------------
      SUBROUTINE unnormalize_fields(dump_time)
      USE fields, ONLY: nrbl,rb
      IMPLICIT NONE
      REAL(r8), INTENT(INOUT) :: dump_time
      INTEGER(i4) :: ibl

      IF (.NOT.set_norms) CALL set_normalizations
      dump_time=dump_time*ta
!-----------------------------------------------------------------------
!     normalize fields
!-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL norm_field_2D(rb(ibl)%be_eq,b0)
        CALL norm_field(rb(ibl)%be,b0)
        CALL norm_field_2D(rb(ibl)%ja_eq,j0)
        CALL norm_field_2D(rb(ibl)%ve_eq,va)
        CALL norm_field(rb(ibl)%ve,va)
        CALL norm_field_2D(rb(ibl)%pres_eq,p0)
        CALL norm_field(rb(ibl)%pres,p0)
        CALL norm_field_2D(rb(ibl)%prese_eq,p0)
        CALL norm_field(rb(ibl)%prese,p0)
        CALL norm_field_2D(rb(ibl)%nd_eq,n0)
        CALL norm_field(rb(ibl)%nd,n0)
        CALL norm_field(rb(ibl)%tele,t0)
        CALL norm_field(rb(ibl)%tion,t0)
      ENDDO
      applied_norms=.FALSE.
      END SUBROUTINE unnormalize_fields
!-----------------------------------------------------------------------
!     subprogram 3. set_normalizations
!     set the normalization parameters
!-----------------------------------------------------------------------
      SUBROUTINE set_normalizations
      USE physdat
      USE input
      USE io, ONLY: out_unit,ofname,nim_wr
      USE pardata, ONLY: node
      IMPLICIT NONE

      IF (k_prp_ref_b>0._r8) THEN
        b0=k_prp_ref_b
      ELSE
        b0=be0 ! defaults to 1 T
      ENDIF
      IF (k_ref_n>0._r8) THEN
        n0=k_ref_n
      ELSE
        n0=ndens ! defaults to 1e20 m^-3
      ENDIF
      va=b0/SQRT(mu0*n0*ms(2))
      ta=1._r8/va
      p0=b0**2/mu0
      T0=p0/(kboltz*n0)
      j0=b0/mu0
!-----------------------------------------------------------------------
!     write out normalization used
!-----------------------------------------------------------------------
      IF (node==0) THEN
        OPEN(UNIT=out_unit,FILE=TRIM(ofname),STATUS='UNKNOWN',          &
             POSITION='APPEND')
        WRITE(nim_wr,'(a,es12.5)') new_line('c')// &
          'Normalizing units using tau_a = ',ta
        WRITE(out_unit,'(a,es12.5)') new_line('c')// &
          'Normalizing units using tau_a = ',ta
        CLOSE(UNIT=out_unit)
      ENDIF
!-----------------------------------------------------------------------
!     normalize diffusivities
!-----------------------------------------------------------------------
      iso_visc=iso_visc*ta/(n0*mtot)
      kin_visc=kin_visc*ta/(n0*mtot)
      par_visc=par_visc*ta/(n0*mtot)
      k_perpi =k_perpi *ta/n0
      k_perpe =k_perpe *ta/n0
      k_perp  =k_perp  *ta/n0
      k_plli  =k_plli  *ta/n0
      k_plle  =k_plle  *ta/n0
      k_pll   =k_pll   *ta/n0
      k_pll_min=k_pll_min*ta/n0
      k_pll_max=k_pll_max*ta/n0
      elecd=elecd/va
      divbd=divbd/va
      hyp_eta=hyp_eta/va
      hyp_dbd=hyp_dbd/va
      nd_diff =nd_diff /va
      nd_hypd =nd_hypd /va
!-----------------------------------------------------------------------
!     normalize input times
!-----------------------------------------------------------------------
      dtm=dtm/ta
      dt_initial=dt_initial/ta
      tmax=tmax/ta
!-----------------------------------------------------------------------
!     normalize constants
!-----------------------------------------------------------------------
      elementary_q=SQRT(mu0*n0*elementary_q**2/ms(2))
      qs(1)=     -elementary_q
      qs(2)= zeff*elementary_q
      mu0=1._r8
      kboltz=1._r8
      ms(1)=ms(1)/ms(2)
      ms(2)=1._r8
      mtot=ms(1)+ms(2)/zeff
      !clight=1._r8
      !eps0=1._r8

      set_norms=.TRUE.

      END SUBROUTINE set_normalizations
!-----------------------------------------------------------------------
!     subprogram 4. norm_field
!     normalize the field by norm
!-----------------------------------------------------------------------
      SUBROUTINE norm_field(field,norm)
      USE lagr_type_mod, ONLY: lagr_quad_type
      IMPLICIT NONE
      TYPE(lagr_quad_type), INTENT(INOUT) :: field
      REAL(r8), INTENT(IN) :: norm

      field%fs =field%fs *norm
      field%fsh=field%fsh*norm
      field%fsv=field%fsv*norm
      field%fsi=field%fsi*norm

      END SUBROUTINE norm_field
!-----------------------------------------------------------------------
!     subprogram 5. norm_field_2D
!     normalize the field by norm
!-----------------------------------------------------------------------
      SUBROUTINE norm_field_2D(field,norm)
      USE lagr_type_mod, ONLY: lagr_quad_2D_type
      IMPLICIT NONE
      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: field
      REAL(r8), INTENT(IN) :: norm

      field%fs =field%fs *norm
      field%fsh=field%fsh*norm
      field%fsv=field%fsv*norm
      field%fsi=field%fsi*norm

      END SUBROUTINE norm_field_2D
!-----------------------------------------------------------------------
!     close module                                                      
!-----------------------------------------------------------------------
      END MODULE normalize_mod
