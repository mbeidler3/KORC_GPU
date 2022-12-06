!-----------------------------------------------------------------------
!     $Id: global.f90 6251 2018-11-05 13:27:58Z jking $
!     contains run-time information that is not directly associated     
!     with input.                                                       
!-----------------------------------------------------------------------
      MODULE global 
      USE local 
      IMPLICIT NONE 
                                                                        
      REAL(r8) :: t,dt=0,dt_old=0,delta2,kdivb_2,smallnum 
      REAL(r8) :: fl_cfl=0,lin_cfl=0,nl_cfl=0 
      REAL(r8) :: cross_section,total_volume,cross_s_overr 
      REAL(r8) :: i_eq,i_n0=0,flux_eq,flux_n0,ff,theta 
      REAL(r8) :: volt=0,volt_old=0,i_n0_old 
      REAL(r8) :: e_vert=0,e_vert_old=0,tmatscale=1._r8
      REAL(r8), DIMENSION(:), POINTER :: keff,keff_total,k2ef 
      INTEGER(i4) :: istep,mhdits=0,hallits=0,jaits=0,                  &
     &               vmhdits=0,viscits=0,teits=0,tiits=0,dbits=0,       &
     &               ndits=0,neoits=0,vnits=0,nstop,                    &
     &               jmode,nmodes,nmodes_total,nsym_pre_band,           &
     &               ncontb,ndiscb
      INTEGER(i4) :: n_dt_increase=0 
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: mps_block,nindex,       &
     &             nindex_total,ipolst_block,ipolen_block,mpsq_block,   &
     &             ipqst_block,ipqen_block                              
      CHARACTER(32) :: integrand_flag 
      LOGICAL :: b0_changed=.false.,p0_changed=.false.,                 &
     &           nl_changed=.false.,n0_changed=.false.,                 &
     &           eta_changed=.false.,kpll_changed=.false.,              &
     &           v0_changed=.false.,surface_e=.false.,                  &
     &           ti0_changed=.false.,te0_changed=.false.,               &
     &           nn0_changed=.false.,vn0_changed=.false.,               &
     &           tn0_changed=.false.,nln_changed=.false.,               &
     &           pn0_changed=.false.
      LOGICAL :: impladv=.false.,threedeta=.false.,closure_calc=.false.,&
     &           closure_n0_only=.true.,disc_nd=.false.,                &
     &           fgnimeq_mode=.false.
      INTEGER(i4) :: nprobe,ompcnt=0_i4
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: probe_ibl
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: probe_rzp,probe_xy
!-----------------------------------------------------------------------
!     close module                                                      
!-----------------------------------------------------------------------
      END MODULE global 
