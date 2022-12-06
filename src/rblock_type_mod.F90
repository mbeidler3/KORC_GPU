!-----------------------------------------------------------------------
!     $Id: rblock_type_mod.F90 7623 2022-01-18 03:55:24Z tbechtel $
!     module containing rblock_type definition.
!-----------------------------------------------------------------------
!     2. dump_read_rblock.
!     4. h5_dump_rblock.
!     5. h5_dump_rbdum.
!     6. h5_dump_rdum_2D
!     7. h5_dump_rdum
!     8. h5_dump_rdum_edge_1D
!     9. h5_dump_rdum_edge
!     10. h5_read_rblock.
!-----------------------------------------------------------------------
!     declarations.
!-----------------------------------------------------------------------
#include "config.f"
      MODULE rblock_type_mod
      USE lagr_quad_mod
      USE lagr_edge_mod
      USE lagr_disc_mod
      USE modal_disc_mod
      IMPLICIT NONE
!-----------------------------------------------------------------------
!     types for saving data at quadrature points
!-----------------------------------------------------------------------
      TYPE :: rb_real_qp_type
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: qpf,qpfr,qpfz
      END TYPE rb_real_qp_type

      TYPE :: rb_real_big_qp_type
        REAL(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: qpf,qpfr,qpfz
      END TYPE rb_real_big_qp_type

      TYPE :: rb_real_bigger_qp_type
        REAL(r8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: qpf,qpfr,qpfz
      END TYPE rb_real_bigger_qp_type

      TYPE :: rb_comp_qp_type
        COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: qpf,qpfr,qpfz
      END TYPE rb_comp_qp_type

      TYPE :: rb_comp_big_qp_type
        COMPLEX(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: qpf,qpfr,qpfz
      END TYPE rb_comp_big_qp_type

      ! For surface integrals.
      TYPE :: rb_real_qpe_type
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: qpf,qpfr,qpfz
      END TYPE rb_real_qpe_type

      TYPE :: rb_comp_qpe_type
        COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: qpf,qpfr,qpfz
      END TYPE rb_comp_qpe_type
!-----------------------------------------------------------------------
!     types for holding basis-function values and their derivatives
!     at quadrature points.
!-----------------------------------------------------------------------
      TYPE :: rb_basis_type
        INTEGER(i4) :: poly_deg_basis
        INTEGER(i4) :: poly_degmin_basis
        INTEGER(i4) :: poly_degmax_basis
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: alpha,dalpdr,dalpdz, &
     &            dalpdrc,alpham,dalpmdr,dalpmdz,dalpmdrc
      END TYPE rb_basis_type
!-----------------------------------------------------------------------
!     data saved in rb types
!     The *_sv variables are for the rt_transfer_eq option.  Tion_sv
!      and tele_sv aren't strictly necessary for dump file consistency,
!      but are used in the sources
!     mxi and myi for nonuniform domain decomposition
!-----------------------------------------------------------------------
      TYPE :: rblock_type
        CHARACTER(64) :: name
        INTEGER(i4) :: mx_total, my_total, nxbl_total, nybl_total
        INTEGER(i4) :: id
        INTEGER(i4) :: mx,my,mxe
        INTEGER(i4) :: ng,nge
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: mxi,myi
        LOGICAL :: degenerate,r0block,periodicy

        TYPE(lagr_edge_1D_type) :: rzedge,norm,tang,grdn,lhss,          &
     &                             beq_edge,neq_edge,veq_edge,peq_edge, &
     &                             nneq_edge,vneq_edge,lambda_edge 
        TYPE(lagr_edge_type)    :: e_surface,b_edge,displacement,b_rmp, &
     &                             chi,dchi,ework,n_edge,v_edge,p_edge, &
     &                             nn_edge,vn_edge,pn_edge

        TYPE(lagr_quad_2D_type) :: rz,be_eq,ja_eq,ve_eq,pres_eq         &
     &            ,prese_eq,nd_eq,diff_shape,pflux,vee_eq,trap,ndn_eq   &
     &            ,vn_eq,tn_eq,pn_eq,F_diff_shape
        TYPE(lagr_quad_2D_type) :: be_sv,ja_sv,ve_sv,pres_sv,prese_sv   &
     &            ,tion_sv,tele_sv,ndn_sv,vn_sv,tn_sv,psi_sv,rz_sv
        TYPE(lagr_quad_2D_type) :: be_n0,pres_n0,nd_n0                  &
     &            ,si_nl_pres,tele_eq,bn0mod,ndn_n0,vn_n0,tn_n0,pn_n0   &
     &            ,tion_eq,nd_eq2,rwork1,rwork2,rwork3,ve_n0,vee_n0     &
     &            ,psi_eq,ti_n0,te_n0,Fion_eq,Fele_eq,Fhot_eq,eef_eq    &
     &            ,pmom_n0,insep,minsol,thot_eq,nhot_eq                 &
     &            ,phot_eq,vhot_eq,G_eq,rz_ph,th_tp,partvar_eq          &
     &            ,bdglnB,ramos2,exb,piv,keflx,qi0,twist_eq,fsa_beq2
        TYPE(lagr_quad_type) :: be,ja,ve,pres,prese,nd,tele,tion        &
     &            ,conc,work1,work2,work3,work4,work5,work6,work7,work8 &
     &            ,qpe,dqpe,qpi,dqpi,ppe,ppi,rpe,rpi,phot,phot1,bmod,eef&
     &            ,curv,dppi,dppe                                       &
     &            ,nZzt,nZnt,nZrt,nZel,nZTt,nz,nZeff                    &
     &            ,qlosd,qlosb,qlosr,qlosl,qlosi                        &
     &            ,Fele,Fion,Fhot,Fwork3,vee                            &
     &            ,Fework1,Fework2,Fiwork1,Fiwork2,Fhwork1,Fhwork2      &
     &            ,qmfeme,qmfemi,ephi,jhot,gradB,w6v1,w6v2,qi,qe,divpi  &
     &            ,gradp,ndn,vn,tn,pn,poten,partvar,Fework1_1,Fiwork1_1 &
     &            ,deef,twist_ln,ephip,flux,gei,BdB1,fp
        TYPE(modal_quad_type) :: auxv,auxb,mwork1,mwork2,mwork3,mwork4

        TYPE(rb_real_qp_type) :: qbe_eq,qja_eq,qve_eq,qpres_eq          &
     &            ,qprese_eq,qnd_eq,qbe_n0,qpres_n0,qnd_n0              &
     &            ,qtion_tot,qtele_tot,qsi_nl_pres                      &
     &            ,qbb,qtele_eq,qtion_eq,qrwork1,qbe_tot                &
     &            ,qnd_tot,qelecd_phi,qelecd_n0,qelecd_eq,qdiff_shape   &
     &            ,qkaprpi_phi,qkaprpe_phi,qkappli_phi,qkapple_phi      &
     &            ,qkappli_n0,qkapple_n0,qkaprpi_n0,qkaprpe_n0          &
     &            ,qte_b2,qti_b2,qbcrgte,qbcrgti                        &
     &            ,qve_n0,qve_tot,qgrdv,qja_tot,qdvv_eq,qeq_force       &
     &            ,qdart,qupw_phi,qupw_n0,qvv,qupti_phi,qupti_n0        &
     &            ,qupte_phi,qupte_n0,qti_n0,qte_n0,qgrdveq,qpi_veq     &
     &            ,qpi_pareq,qpi_gyreq,qvee_eq,qgrdve,qgrdveeq          &
     &            ,qvee_tot,qvee_n0,qndn_eq,qndn_n0,qndn_tot            &
     &            ,qvn_eq,qvn_n0,qvn_tot,qtn_eq,qtn_n0,qtn_tot          &
     &            ,qpn_eq,qpn_n0,qdvvn_eq,qpin_eq                       &
     &            ,qFele_eq,qFion_eq,qFhot_eq,qeef_eq,qF_diff_shape     &
     &            ,qsion_tot,qsion_sym,qsrec_tot,qsrec_sym              &
     &            ,qscx_tot,qscx_sym,qsnn_eq,qivn_eq,qrvn_eq,qcvn_eq    &
     &            ,qcxf_eq,qsion_eq,qsrec_eq,qscx_eq,qstn_eq,qsti_eq    &
     &            ,qsvn_eq,qsve_eq,qpmom_n0,qgrdvn_tot,qgrdvn_eq        &
     &            ,qnZTt_tot,qnimp_tot,qthot_eq,qnhot_eq                &
     &            ,qgrdb,qpr_tot,qgrdp,qphot_eq                         &
     &            ,qnu_ei,qrpe_eq,qppe_eq,qqpe_eq,qrpi_eq,qppi_eq       &
     &            ,qqpi_eq,qDe0,qDe1,qmDe,qDi0,qDi1,qti_dtn0,qte_dtn0   &
     &            ,qupz_phi,qupz_n0                                     &
     &            ,qnZel_tot,qnZzt_tot,qnZnt_tot,qnZrt_tot,qfsa_beq2    &
     &            ,qnZrt_n0,qnZTt_n0,qnZeff_tot,qrgeom_tot,qneodiff     &
     &            ,qndiff_phi,qndiff_n0
        TYPE(rb_comp_qp_type) :: qbe,qja,qve,qpres,qprese,qnd,qtele     &
     &            ,qtion,qconc,qvisc,qbe_off,qja_off,qve_off,qnd_off    &
     &            ,qwork1,qwork2,qwork3,qwork4,qvee,qevisc,qephi,qephip &
     &            ,qtion_dt,qtele_dt,qndn,qvn,qtn,qpn,qbb1,qflux        &
     &            ,qkapple_off,qkaprpe_off,qkappli_off,qkaprpi_off      &
     &            ,qneut_vheat,qfp,qndiff,qndiffa
        TYPE(rb_basis_type), DIMENSION (:), ALLOCATABLE :: base_pd      &
     &            ,base_disc,base_modal
        ! This is used by nimset
        TYPE(rb_real_qp_type) :: qrz,qpsi_eq,qtrap
        TYPE(rb_comp_qp_type) :: qpsi
        ! Sources - and things for sources.
        TYPE(rb_real_qp_type) :: qja_sv,qve_sv,qtele_sv,qtion_sv,qnd_sv,&
     &                           qrho_tot,qrho_sym
        TYPE(rb_real_qp_type) :: qe_src_n0,qqi_src_n0,qqe_src_n0,       &
     &                           qnd_src_n0,qve_src_n0
        TYPE(rb_comp_qp_type) :: qe_source,qqi_source,qqe_source,       &
     &                           qnd_source,qve_source,qrho
        TYPE(rb_real_qp_type) :: qrelax_mask
        ! Closures
        TYPE(rb_comp_qp_type) :: qqpe,qdqpe,qqpi,qdqpi,qppe,qppi,qrpe   &
     &         ,qphot,qjhot,qnimp,qnz,qgei,qdppi,qdppe                  &
     &         ,qnZzt,qnZnt,qnZrt,qnZel,qnZTt,qnZeff,qrgeom             &
     &         ,qqlosd,qqloso,qqlosb,qqlosr,qqlosl,qqlosi               &
     &         ,qFele,qFion,qFhot,qFwork3,qrpi,qBdB1
        TYPE(rb_real_qp_type) :: qlmfpi_n0,qlmfpe_n0,qfmi_eq,qfme_eq

        TYPE(rb_real_bigger_qp_type) :: qcee,qcei,qcii,qcie,qddTcee,    &
     &                                  qddTcii,qcee_1,qcii_1
        TYPE(rb_real_bigger_qp_type) :: qchion,qchele
        TYPE(rb_real_big_qp_type) :: qcee_n0,qcei_n0,qcii_n0,qcie_n0
        TYPE(rb_real_big_qp_type) :: qcee_lk,qcei_lk,qcii_lk,qcie_lk
        TYPE(rb_real_big_qp_type) :: qcaee_n0,qcaii_n0 ! like only
        TYPE(rb_real_big_qp_type) :: qcaee_lk,qcaii_lk ! like only
        TYPE(rb_real_big_qp_type) :: qmei

        TYPE(rb_real_big_qp_type) :: qvde_eq ! bulk electrons
        TYPE(rb_real_big_qp_type) :: qvdi_eq ! bulk ions
        TYPE(rb_real_big_qp_type) :: qvdh_eq ! hot particles
        TYPE(rb_real_big_qp_type) :: qaxi_eq, qasi_eq ! ion acceleration
        TYPE(rb_real_big_qp_type) :: qaxe_eq, qase_eq ! ele acceleration
        TYPE(rb_real_big_qp_type) :: qaxh_eq, qash_eq ! hot acceleration

        TYPE(rb_real_qp_type) :: qcLe_n0,qcLi_n0
        TYPE(rb_real_qp_type) :: qcLee_n0,qcLei_n0,qcLei_t_n0
        TYPE(rb_real_big_qp_type) :: qdrift_n0
        TYPE(rb_comp_big_qp_type) :: qFe_mom,qFi_mom,qdrift,qdriftn
        TYPE(rb_real_qp_type) :: qbn0mod,qcurv_n0
        TYPE(rb_comp_qp_type) :: qdivbb,qbdglnB,qdeef
        TYPE(rb_real_qp_type) :: qcurv_eq,qdivbb_eq,qbdglnB_eq,qtwist_eq
        TYPE(rb_comp_qp_type) :: qbmod,qcurv,qeef,qgradB
        TYPE(rb_comp_qp_type) :: qcurv_ln,qdivbb_ln,qbdglnB_ln,qbhat_ln
        TYPE(rb_comp_qp_type) :: qglnB_ln,qelecd_ln,qtwist_ln
        TYPE(rb_comp_qp_type) :: qCe_mom,qCi_mom

        ! Surface quantities
        TYPE(rb_real_qpe_type) :: qrzedge,qnorm,qtang,qds,qbeq_edge,    &
     &            qneq_edge,qveq_edge,qnneq_edge,qvneq_edge,qpeq_edge
        TYPE(rb_real_qpe_type) :: qlambda
        TYPE(rb_comp_qpe_type) :: qb_edge,qe_surface,pn_surface,qchi,   &
     &            qn_edge,qv_edge,qnn_edge,qvn_edge,qp_edge,qpn_edge,   &
     &            qDchi,qb_rmp

        REAL(r8), DIMENSION(:), ALLOCATABLE :: xg,yg,wg,xge,wge
        REAL(r8), DIMENSION(:,:), ALLOCATABLE :: dxdr,dxdz,dydr,dydz,   &
     &            bigr,wjac,jac2d,wgrin2nim
        REAL(r8), DIMENSION(:,:), ALLOCATABLE :: cell_vol,r_cent
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: wall_matrix,         &
     &            deriv_matrix

      END TYPE rblock_type

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 2. dump_read_rblock.
!     reads rblock data from dump file.
!-----------------------------------------------------------------------
      SUBROUTINE dump_read_rblock(rb)
      USE input, ONLY:phdump,nF,zimp,mx,my,F_dvac,dump_psi_eq,          &
     &       nxbl,nybl,closure_model,neutral_model,dump_read_F,         &
     &       dump_eexp,dump_fsa_beq2,lagr_edge_rmp,use_fp
      USE io

      TYPE (rblock_type), INTENT(OUT) :: rb
      REAL(r8) :: iread
      REAL(r8) :: rread
!-----------------------------------------------------------------------
!     read block descriptors.
!-----------------------------------------------------------------------
      rb%name='rblock'
      READ(rstrt_unit) iread
      rb%id=NINT(iread)
      READ(rstrt_unit) iread
      rb%mx=NINT(iread)
      READ(rstrt_unit) iread
      rb%my=NINT(iread)
      READ(rstrt_unit) rread
      IF (rread>0) THEN
        rb%degenerate=.true.
      ELSE
        rb%degenerate=.false.
      ENDIF
!-----------------------------------------------------------------------
!     Give me information about the global arrays comming from input.f
!-----------------------------------------------------------------------
      rb%mx_total=mx
      rb%my_total=my
      rb%nxbl_total=nxbl
      rb%nybl_total=nybl
!-----------------------------------------------------------------------
!     read coordinates.
!-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%rz,rb%mx,rb%my,'rz',               &
     &                            (/'   r  ','   z  '/))
!-----------------------------------------------------------------------
!     read magnetic fields.
!-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%be_eq,rb%mx,rb%my,'bq',            &
     &                            (/' be_eq'/))
      CALL dump_read_lagr_quad(rb%be,rb%mx,rb%my,'be',(/'  be  '/))
!-----------------------------------------------------------------------
!     read equilibrium current density.
!-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%ja_eq,rb%mx,rb%my,'jq',            &
     &                            (/' ja_eq'/))
!-----------------------------------------------------------------------
!     read fluid velocities.
!-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%ve_eq,rb%mx,rb%my,'vq',            &
     &                            (/' ve_eq'/))
      CALL dump_read_lagr_quad(rb%ve,rb%mx,rb%my,'ve',(/'  ve  '/))
!-----------------------------------------------------------------------
!     read pressures.
!-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%pres_eq,rb%mx,rb%my,'pq',          &
     &                            (/' pr_eq'/))
      CALL dump_read_lagr_quad(rb%pres,rb%mx,rb%my,'pr',(/' pres '/))
      CALL dump_read_lagr_quad_2D(rb%prese_eq,rb%mx,rb%my,'pq',         &
     &                            (/'pre_eq'/))
      CALL dump_read_lagr_quad(rb%prese,rb%mx,rb%my,'pe',(/' prese'/))
!-----------------------------------------------------------------------
!     read number densities.
!-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%nd_eq,rb%mx,rb%my,'nq',            &
     &                            (/' nd_eq'/))
      CALL dump_read_lagr_quad(rb%nd,rb%mx,rb%my,'nd',(/'  nd  '/))
!-----------------------------------------------------------------------
!     read diffusivity shape factor and material concentration.
!-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%diff_shape,rb%mx,rb%my,'ds',       &
     &                            (/'dif_sh'/))
      CALL dump_read_lagr_quad(rb%conc,rb%mx,rb%my,'co',(/' conc '/))
!-----------------------------------------------------------------------
!     read temperatures.
!-----------------------------------------------------------------------
      CALL dump_read_lagr_quad(rb%tele,rb%mx,rb%my,'te',(/' tele '/))
      CALL dump_read_lagr_quad(rb%tion,rb%mx,rb%my,'ti',(/' tion '/))
!-----------------------------------------------------------------------
!     read hot particle pressure tensor or current.
!-----------------------------------------------------------------------
      IF (phdump) THEN
        CALL dump_read_lagr_quad(rb%phot,rb%mx,rb%my,'phot',            &
     &         (/'phot  '/))
      ENDIF
!-----------------------------------------------------------------------
!     read kinetic distortions and closures.
!-----------------------------------------------------------------------
      IF (nF(1)>0.AND.dump_read_F(1)) THEN
       CALL dump_read_lagr_quad_2D(rb%Fele_eq,rb%mx,rb%my,'Feq',        &
     &                         (/' Fe_eq'/))
       CALL dump_read_lagr_quad(rb%Fele,rb%mx,rb%my,'Fele',(/' Fele '/))
       CALL dump_read_lagr_quad(rb%qpe,rb%mx,rb%my,'qpe',(/' qpe  '/))
       CALL dump_read_lagr_quad(rb%ppe,rb%mx,rb%my,'ppe',(/' ppe  '/))
       CALL dump_read_lagr_quad(rb%rpe,rb%mx,rb%my,'rpe',(/' rpe  '/))
       CALL dump_read_lagr_quad(rb%gei,rb%mx,rb%my,'gei',(/' gei  '/))
      ENDIF
      IF (nF(2)>0.AND.dump_read_F(2)) THEN
       CALL dump_read_lagr_quad_2D(rb%Fion_eq,rb%mx,rb%my,'Fiq',        &
     &                         (/' Fi_eq'/))
       CALL dump_read_lagr_quad(rb%Fion,rb%mx,rb%my,'Fion',(/' Fion '/))
       CALL dump_read_lagr_quad(rb%qpi,rb%mx,rb%my,'qpi',(/' qpi  '/))
       CALL dump_read_lagr_quad(rb%ppi,rb%mx,rb%my,'ppi',(/' ppi  '/))
       CALL dump_read_lagr_quad(rb%rpi,rb%mx,rb%my,'rpi',(/' rpi  '/))
      ENDIF
      IF (nF(3)>0.AND.dump_read_F(3)) THEN
       CALL dump_read_lagr_quad_2D(rb%Fhot_eq,rb%mx,rb%my,'Fhq',        &
     &                         (/' Fh_eq'/))
       CALL dump_read_lagr_quad(rb%Fhot,rb%mx,rb%my,'Fhot',(/' Fhot '/))
       CALL dump_read_lagr_quad_2D(rb%nhot_eq,rb%mx,rb%my,'nhq',        &
     &                         (/' nh_eq'/))
       CALL dump_read_lagr_quad_2D(rb%thot_eq,rb%mx,rb%my,'thq',        &
     &                         (/' th_eq'/))
      ENDIF
      IF (F_dvac>0.AND.(dump_read_F(1).OR.dump_read_F(2)                &
     &                                .OR.dump_read_F(3))) THEN
       CALL dump_read_lagr_quad_2D(rb%F_diff_shape,rb%mx,rb%my,'Fds',   &
     &                             (/'Fdifsh'/))
      ENDIF
!-----------------------------------------------------------------------
      IF (zimp/=0) THEN
      CALL dump_read_lagr_quad(rb%nZzt,rb%mx,rb%my,'nZzt',(/' nZzt '/))
!      CALL dump_read_lagr_quad(rb%nZnt,rb%mx,rb%my,'nZnt',(/' nZnt '/))
!      CALL dump_read_lagr_quad(rb%nZrt,rb%mx,rb%my,'nZrt',(/' nZrt '/))
      CALL dump_read_lagr_quad(rb%nZel,rb%mx,rb%my,'nZel',(/' nZel '/))
      CALL dump_read_lagr_quad(rb%nz,rb%mx,rb%my,  'nz',  (/'  nz  '/))
      CALL dump_read_lagr_quad(rb%nZeff,rb%mx,rb%my,'nZeff',            &
     &                            (/' nZeff'/))
      CALL dump_read_lagr_quad(rb%qlosd,rb%mx,rb%my,'qlosd',            &
     &                            (/' qlosd'/))
      CALL dump_read_lagr_quad(rb%qlosb,rb%mx,rb%my,'qlosb',            &
     &                            (/' qlosb'/))
      CALL dump_read_lagr_quad(rb%qlosr,rb%mx,rb%my,'qlosr',            &
     &                            (/' qlosr'/))
      CALL dump_read_lagr_quad(rb%qlosl,rb%mx,rb%my,'qlosl',            &
     &                            (/' qlosl'/))
      CALL dump_read_lagr_quad(rb%qlosi,rb%mx,rb%my,'qlosi',            &
     &                            (/' qlosi'/))
      CALL dump_read_lagr_quad(rb%ndn,rb%mx,rb%my,'ndn',   (/' ndn  '/))
      IF (dump_eexp) &
     &  CALL dump_read_lagr_quad(rb%eef,rb%mx,rb%my,'ef',   (/' Efld '/))
!      IF (dump_ja) &
!     &  CALL dump_read_lagr_quad(rb%ja,rb%mx,rb%my,'ja',   (/' ja   '/))
      ENDIF
!-----------------------------------------------------------------------
!     read neutrals data.
!-----------------------------------------------------------------------
      IF (neutral_model/='none') THEN
        CALL dump_read_lagr_quad_2D(rb%ndn_eq,rb%mx,rb%my,'ndn_eq',     &
     &                                                (/'ndn_eq'/))
        CALL dump_read_lagr_quad(rb%ndn,rb%mx,rb%my,'ndn',(/' ndn  '/))
        CALL dump_read_lagr_quad_2D(rb%vn_eq,rb%mx,rb%my,'vn_eq',       &
     &                                                (/'vn_eq '/))
        CALL dump_read_lagr_quad(rb%vn,rb%mx,rb%my,'vn',(/'  vn  '/))
        CALL dump_read_lagr_quad_2D(rb%tn_eq,rb%mx,rb%my,' tn_eq',      &
     &                                                (/'tn_eq '/))
        CALL dump_read_lagr_quad(rb%tn,rb%mx,rb%my,'tn',(/'  tn  '/))
      ENDIF
!-----------------------------------------------------------------------
!     read psi_eq data.
!-----------------------------------------------------------------------
      IF (dump_psi_eq) THEN
        CALL dump_read_lagr_quad_2D(rb%psi_eq,rb%mx,rb%my,'psi_eq',     &
     &                              (/'psi_eq'/))
      ENDIF
!-----------------------------------------------------------------------
!     read fsa_beq2 data.
!-----------------------------------------------------------------------
      IF (dump_fsa_beq2) THEN
        CALL dump_read_lagr_quad_2D(rb%fsa_beq2,rb%mx,rb%my,'fsbeq2',   &
     &                              (/'fsbeq2'/))
      ENDIF
!-----------------------------------------------------------------------
!     write edge coordinates if lagr_edge_rmp
!-----------------------------------------------------------------------
      IF (lagr_edge_rmp) THEN
        READ(rstrt_unit) iread
        rb%mxe=NINT(iread)
        IF (rb%mxe>0) THEN
          CALL dump_read_lagr_edge_1D(rb%rzedge,rb%mxe,'rzedge',        &
     &                                 (/'rzedge'/))
!-----------------------------------------------------------------------
!     write RMP magnetic fields.
!-----------------------------------------------------------------------
          CALL dump_read_lagr_edge(rb%b_rmp,rb%mxe,'brmp',(/' brmp '/))
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     read ponderomotive force.
!-----------------------------------------------------------------------
      IF (use_fp)                                                       &
     &  CALL dump_read_lagr_quad(rb%fp,rb%mx,rb%my,'fp',(/'  fp  '/))
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dump_read_rblock
#ifdef HAVE_FC_HDF5
!-----------------------------------------------------------------------
!     subprogram 4. h5_dump_rblock.
!     writes rblock data to h5 dump file.
!-----------------------------------------------------------------------
      SUBROUTINE h5_dump_rblock(rb,rbid)
      USE input, ONLY: phdump,nxbl,nF,zimp,poly_degree,dump_psi_eq,     &
     &                 closure_model,neutral_model,dump_eexp,F_dvac,    &
     &                 dump_ja,dump_q,dump_divpi,dump_fsa_beq2,         & 
     &                 separate_pe,ph_pd,lagr_edge_rmp,use_fp
      USE global, ONLY: nmodes
      USE io

      TYPE (rblock_type), INTENT(IN) :: rb
      INTEGER(HID_T), INTENT(IN) :: rbid
      CHARACTER(4) :: blid
      INTEGER(HID_T) :: gid
!-----------------------------------------------------------------------
!     setup rblock group.
!-----------------------------------------------------------------------
      WRITE(blid,fmt='(i4.4)') rb%id
      CALL make_group(rbid,TRIM(blid),gid,h5in,h5err)
!-----------------------------------------------------------------------
!     write block descriptors.
!-----------------------------------------------------------------------
      CALL write_attribute(gid,"id",rb%id,h5in,h5err)
      CALL write_attribute(gid,"mx",rb%mx,h5in,h5err)
      CALL write_attribute(gid,"my",rb%my,h5in,h5err)
      CALL write_attribute(gid,"nfour",rb%be%nfour,h5in,h5err)
      CALL write_attribute(gid,"poly_degree",rb%be%n_side+1,h5in,h5err)
      IF (rb%degenerate) THEN
        CALL write_attribute(gid,"degenerate",1_i4,h5in,h5err)
      ELSE
        CALL write_attribute(gid,"degenerate",0_i4,h5in,h5err)
      ENDIF
      CALL write_attribute(gid,"mxe",rb%mxe,h5in,h5err)
!-----------------------------------------------------------------------
!     write coordinates.
!-----------------------------------------------------------------------
      h5in%mesh="mesh-structured"
      h5in%vsAxisLabels="R,Z"
      ! see comment in h5_dump_lagr_quad_2D on ordering
!      h5in%vsIndexOrder="compMajorF"
      CALL h5_dump_lagr_quad_2D(rb%rz,"rz",gid,blid)
!-----------------------------------------------------------------------
!     write magnetic fields.
!-----------------------------------------------------------------------
      h5in%mesh="/rblocks/"//TRIM(blid)//"/rz"//TRIM(blid)
      h5in%vsTimeGroup="/dumpTime"
      CALL h5_dump_lagr_quad_2D(rb%be_eq,"bq",gid,blid)
      CALL h5_dump_lagr_quad(rb%be,"be",gid,blid)
!-----------------------------------------------------------------------
!     write equilibrium current density.
!-----------------------------------------------------------------------
      CALL h5_dump_lagr_quad_2D(rb%ja_eq,"jq",gid,blid)
!-----------------------------------------------------------------------
!     write fluid velocities.
!-----------------------------------------------------------------------
      CALL h5_dump_lagr_quad_2D(rb%ve_eq,"vq",gid,blid)
      CALL h5_dump_lagr_quad(rb%ve,"ve",gid,blid)
!-----------------------------------------------------------------------
!     write pressures.
!-----------------------------------------------------------------------
      CALL h5_dump_lagr_quad_2D(rb%pres_eq,"prq",gid,blid)
      CALL h5_dump_lagr_quad(rb%pres,"pr",gid,blid)
      CALL h5_dump_lagr_quad_2D(rb%prese_eq,"peq",gid,blid)
      CALL h5_dump_lagr_quad(rb%prese,"pe",gid,blid)
!-----------------------------------------------------------------------
!     write number densities.
!-----------------------------------------------------------------------
      CALL h5_dump_lagr_quad_2D(rb%nd_eq,"nq",gid,blid)
      CALL h5_dump_lagr_quad(rb%nd,"nd",gid,blid)
!-----------------------------------------------------------------------
!     write diffusivity shape function and material concentration.
!-----------------------------------------------------------------------
      CALL h5_dump_lagr_quad_2D(rb%diff_shape,"diff",gid,blid)
      CALL h5_dump_lagr_quad(rb%conc,"conc",gid,blid)
!-----------------------------------------------------------------------
!     write temperatures.
!-----------------------------------------------------------------------
      CALL h5_dump_lagr_quad(rb%tele,"te",gid,blid)
      CALL h5_dump_lagr_quad(rb%tion,"ti",gid,blid)
!-----------------------------------------------------------------------
!     write hot particle pressure tensor or current.
!-----------------------------------------------------------------------
      IF (phdump) THEN 
        IF (ph_pd/=poly_degree.AND.ph_pd==1) THEN
          h5in%mesh="mesh-structured"
          CALL h5_dump_lagr_quad_2D(rb%rz_ph,"rz_ph",gid,blid)
          h5in%mesh="/rblocks/"//TRIM(blid)//"/rz_ph"//TRIM(blid)
        ENDIF
        CALL h5_dump_lagr_quad(rb%phot,"phot",gid,blid)
        h5in%mesh="/rblocks/"//TRIM(blid)//"/rz"//TRIM(blid)
      ENDIF
!-----------------------------------------------------------------------
!     write kinetic distortions and closures.
!-----------------------------------------------------------------------
      IF (nF(1)>0) THEN
        CALL h5_dump_lagr_quad_2D(rb%Fele_eq,"Fele_eq",gid,blid)
        CALL h5_dump_lagr_quad(rb%Fele,"Fele",gid,blid)
        CALL h5_dump_lagr_quad(rb%qpe,"qpe",gid,blid)
        CALL h5_dump_lagr_quad(rb%ppe,"ppe",gid,blid)
        CALL h5_dump_lagr_quad(rb%rpe,"rpe",gid,blid)
        CALL h5_dump_lagr_quad(rb%gei,"gei",gid,blid)
      ENDIF
      IF (nF(2)>0) THEN
        CALL h5_dump_lagr_quad_2D(rb%Fion_eq,"Fion_eq",gid,blid)
        CALL h5_dump_lagr_quad(rb%Fion,"Fion",gid,blid)
        CALL h5_dump_lagr_quad(rb%qpi,"qpi",gid,blid)
        CALL h5_dump_lagr_quad(rb%ppi,"ppi",gid,blid)
        CALL h5_dump_lagr_quad(rb%rpi,"rpi",gid,blid)
      ENDIF
      IF (nF(3)>0) THEN
        CALL h5_dump_lagr_quad_2D(rb%Fhot_eq,"Fhot_eq",gid,blid)
        CALL h5_dump_lagr_quad(rb%Fhot,"Fhot",gid,blid)
        CALL h5_dump_lagr_quad_2D(rb%nhot_eq,"nhot_eq",gid,blid)
        CALL h5_dump_lagr_quad_2D(rb%thot_eq,"thot_eq",gid,blid)
      ENDIF
      IF (F_dvac>0) CALL h5_dump_lagr_quad_2D(rb%F_diff_shape,"F_diff", &
     &                                        gid,blid)
!-----------------------------------------------------------------------
!     write kprad data.
!-----------------------------------------------------------------------
      IF (zimp/=0) THEN
        CALL h5_dump_lagr_quad(rb%nZzt,"nZzt",gid,blid)
!        CALL h5_dump_lagr_quad(rb%nZnt,"nZnt",gid,blid)
!        CALL h5_dump_lagr_quad(rb%nZrt,"nZrt",gid,blid)
!        CALL h5_dump_lagr_quad(rb%nZTt,"nZTt",gid,blid)
        CALL h5_dump_lagr_quad(rb%nZel,"nZel",gid,blid)
        CALL h5_dump_lagr_quad(rb%nz,"nz",gid,blid)
        CALL h5_dump_lagr_quad(rb%nZeff,"nZeff",gid,blid)
!CCK condense into a single array of dimension nqlz
        CALL h5_dump_lagr_quad(rb%qlosd,"qlosd",gid,blid)
        CALL h5_dump_lagr_quad(rb%qlosb,"qlosb",gid,blid)
        CALL h5_dump_lagr_quad(rb%qlosr,"qlosr",gid,blid)
        CALL h5_dump_lagr_quad(rb%qlosl,"qlosl",gid,blid)
        CALL h5_dump_lagr_quad(rb%qlosi,"qlosi",gid,blid)
        CALL h5_dump_lagr_quad(rb%ndn,"ndn",gid,blid)
      ENDIF
!-----------------------------------------------------------------------
!     write neutrals data.
!-----------------------------------------------------------------------
      IF (neutral_model/='none') THEN
        CALL h5_dump_lagr_quad_2D(rb%ndn_eq,"ndn_eq",gid,blid)
        CALL h5_dump_lagr_quad(rb%ndn,"ndn",gid,blid)
        CALL h5_dump_lagr_quad_2D(rb%vn_eq,"vn_eq",gid,blid)
        CALL h5_dump_lagr_quad(rb%vn,"vn",gid,blid)
        CALL h5_dump_lagr_quad_2D(rb%tn_eq,"tn_eq",gid,blid)
        CALL h5_dump_lagr_quad(rb%tn,"tn",gid,blid)
      ENDIF
!-----------------------------------------------------------------------
!     write current, electric field, psi_eq, and fsa_beq^2
!-----------------------------------------------------------------------
      IF (dump_ja) THEN
        CALL h5_dump_lagr_quad(rb%ja,"ja",gid,blid)
      ENDIF
      IF (dump_eexp) THEN
        CALL h5_dump_lagr_quad(rb%eef,"ee",gid,blid)
        CALL h5_dump_lagr_quad_2D(rb%eef_eq,"eq",gid,blid)
      ENDIF
      IF (dump_q) THEN
        CALL h5_dump_lagr_quad(rb%qi,"qi",gid,blid)
        IF (separate_pe) CALL h5_dump_lagr_quad(rb%qe,"qe",gid,blid)
      ENDIF
      IF (dump_divpi) THEN
        CALL h5_dump_lagr_quad(rb%divpi,"divpi",gid,blid)
      ENDIF
      IF (dump_psi_eq) THEN
        CALL h5_dump_lagr_quad_2D(rb%psi_eq,"psi_eq",gid,blid)
      ENDIF
      IF (dump_fsa_beq2) THEN 
        CALL h5_dump_lagr_quad_2D(rb%fsa_beq2,"fsabeq2",gid,blid)
      ENDIF
!-----------------------------------------------------------------------
!     write edge coordinates if lagr_edge_rmp
!-----------------------------------------------------------------------
      IF (lagr_edge_rmp .and. rb%mxe>0) THEN
        h5in%mesh="mesh-structured"
        h5in%vsAxisLabels="R,Z"
        CALL h5_dump_lagr_edge_1D(rb%rzedge,"rzedge",gid,blid)
!-----------------------------------------------------------------------
!     write RMP magnetic fields.
!-----------------------------------------------------------------------
        h5in%mesh="/rblocks/"//TRIM(blid)//"/rzedge"//TRIM(blid)
        CALL h5_dump_lagr_edge(rb%b_rmp,"brmp",gid,blid)
      ENDIF
!-----------------------------------------------------------------------
!     write ponderomotive force.
!-----------------------------------------------------------------------
      IF (use_fp) CALL h5_dump_lagr_quad(rb%fp,"fp",gid,blid)
!-----------------------------------------------------------------------
!     terminate
!-----------------------------------------------------------------------
      h5in%mesh=" "
      h5in%vsAxisLabels=" "
      h5in%vsIndexOrder=" "
      h5in%vsTimeGroup=" "
      CALL close_group(TRIM(blid),gid,h5err)
      RETURN
      END SUBROUTINE h5_dump_rblock
!-----------------------------------------------------------------------
!     subprogram 5. h5_dump_rbdum.
!     writes rblock data to h5 dump file.
!-----------------------------------------------------------------------
      SUBROUTINE h5_dump_rbdum(id,rbid)
      USE input, ONLY: phdump,phqty,nxbl,nF,zimp,poly_degree,           &
     &                 closure_model,neutral_model,dump_eexp,F_dvac,    &
     &                 dump_psi_eq,dump_ja,dump_q,dump_divpi,           &
     &                 dump_fsa_beq2,separate_pe,ph_pd,ds_nqty,         &
     &                 lagr_edge_rmp,use_fp
      USE global, ONLY: nmodes,nmodes_total
      USE pardata, ONLY: block_sizes
      USE io

      INTEGER(i4), INTENT(IN):: id
      INTEGER(HID_T), INTENT(IN) :: rbid
      CHARACTER(4) :: blid
      INTEGER(i4) :: mtmp,nfour
      INTEGER(HID_T) :: gid
!-----------------------------------------------------------------------
!     setup rblock group.
!-----------------------------------------------------------------------
      WRITE(blid,fmt='(i4.4)') id
      CALL make_group(rbid,TRIM(blid),gid,h5in,h5err)
      CALL write_attribute(gid,"id",id,h5in,h5err)
      CALL write_attribute(gid,"mx",block_sizes(3,id),h5in,h5err)
      CALL write_attribute(gid,"my",block_sizes(4,id),h5in,h5err)
      nfour=nmodes_total
      CALL write_attribute(gid,"nfour",nfour,h5in,h5err)
      CALL write_attribute(gid,"poly_degree",poly_degree,h5in,h5err)
      CALL write_attribute(gid,"degenerate",block_sizes(5,id),          &
     &                     h5in,h5err)
      CALL write_attribute(gid,"mxe",block_sizes(6,id),h5in,h5err)
!-----------------------------------------------------------------------
!     write coordinates.
!-----------------------------------------------------------------------
      h5in%mesh="mesh-structured"
      h5in%vsAxisLabels="R,Z"
      CALL h5_dump_rdum_2D("rz",2)
!-----------------------------------------------------------------------
!     write magnetic fields.
!-----------------------------------------------------------------------
      h5in%mesh="/rblocks/"//TRIM(blid)//"/rz"//TRIM(blid)
      h5in%vsTimeGroup="/dumpTime"
      CALL h5_dump_rdum_2D("bq",3)
      CALL h5_dump_rdum("be",nfour,3)
!-----------------------------------------------------------------------
!     write equilibrium current density.
!-----------------------------------------------------------------------
      CALL h5_dump_rdum_2D("jq",3)
!-----------------------------------------------------------------------
!     write fluid velocities.
!-----------------------------------------------------------------------
      CALL h5_dump_rdum_2D("vq",3)
      CALL h5_dump_rdum("ve",nfour,3)
!-----------------------------------------------------------------------
!     write pressures.
!-----------------------------------------------------------------------
      CALL h5_dump_rdum_2D("prq",1)
      CALL h5_dump_rdum("pr",nfour,1)
      CALL h5_dump_rdum_2D("peq",1)
      CALL h5_dump_rdum("pe",nfour,1)
!-----------------------------------------------------------------------
!     write number densities.
!-----------------------------------------------------------------------
      CALL h5_dump_rdum_2D("nq",1)
      CALL h5_dump_rdum("nd",nfour,1)
!-----------------------------------------------------------------------
!     write diffusivity shape function and material concentration.
!-----------------------------------------------------------------------
      CALL h5_dump_rdum_2D("diff",ds_nqty)
      CALL h5_dump_rdum("conc",nfour,1)
!-----------------------------------------------------------------------
!     write temperatures.
!-----------------------------------------------------------------------
      CALL h5_dump_rdum("te",nfour,1)
      CALL h5_dump_rdum("ti",nfour,1)
!-----------------------------------------------------------------------
!     write hot particle pressure tensor or current.
!-----------------------------------------------------------------------
      IF (phdump) THEN
        IF (ph_pd/=poly_degree.AND.ph_pd==1) THEN
          h5in%mesh="mesh-structured"
          CALL h5_dump_rdum("rz_ph",nfour,2)
          h5in%mesh="/rblocks/"//TRIM(blid)//"/rz_ph"//TRIM(blid)
        ENDIF
        CALL h5_dump_rdum("phot",nfour,phqty)
      ENDIF
!-----------------------------------------------------------------------
!     write kinetic distortions and closures.
!-----------------------------------------------------------------------
      IF (nF(1)>0) THEN
        CALL h5_dump_rdum_2D("Fele_eq",nF(1))
        CALL h5_dump_rdum("Fele",nfour,nF(1))
        CALL h5_dump_rdum("qpe",nfour,1)
        CALL h5_dump_rdum("ppe",nfour,1)
        CALL h5_dump_rdum("rpe",nfour,1)
        CALL h5_dump_rdum("gei",nfour,1)
      ENDIF
      IF (nF(2)>0) THEN
        CALL h5_dump_rdum_2D("Fion_eq",nF(2))
        CALL h5_dump_rdum("Fion",nfour,nF(2))
        CALL h5_dump_rdum("qpi",nfour,1)
        CALL h5_dump_rdum("ppi",nfour,1)
        CALL h5_dump_rdum("rpi",nfour,1)
      ENDIF
      IF (nF(3)>0) THEN
        CALL h5_dump_rdum_2D("Fhot_eq",nF(3))
        CALL h5_dump_rdum("Fhot",nfour,nF(3))
        CALL h5_dump_rdum_2D("nhot_eq",1)
        CALL h5_dump_rdum_2D("thot_eq",1)
      ENDIF
      IF (F_dvac>0)                                                     &
     &  CALL h5_dump_rdum_2D("F_diff",1)
!-----------------------------------------------------------------------
!     write kprad data.
!-----------------------------------------------------------------------
      IF (zimp/=0) THEN
        CALL h5_dump_rdum("nZzt",nfour,1)
!        CALL h5_dump_rdum("nZnt",nfour,1)
!        CALL h5_dump_rdum("nZrt",nfour,1)
!        CALL h5_dump_rdum("nZTt",nfour,1)
        CALL h5_dump_rdum("nZel",nfour,1)
        CALL h5_dump_rdum("nz",nfour,INT(zimp+1,i4))
        CALL h5_dump_rdum("nZeff",nfour,1)
        CALL h5_dump_rdum("qlosd",nfour,1)
        CALL h5_dump_rdum("qlosb",nfour,1)
        CALL h5_dump_rdum("qlosr",nfour,1)
        CALL h5_dump_rdum("qlosl",nfour,1)
        CALL h5_dump_rdum("qlosi",nfour,1)
        CALL h5_dump_rdum("ndn",nfour,1)
      ENDIF
!-----------------------------------------------------------------------
!     write neutrals data.
!-----------------------------------------------------------------------
      IF (neutral_model/='none') THEN
        CALL h5_dump_rdum_2D("ndn_eq",1)
        CALL h5_dump_rdum("ndn",nfour,1)
        CALL h5_dump_rdum_2D("vn_eq",3)
        CALL h5_dump_rdum("vn",nfour,3)
        CALL h5_dump_rdum_2D("tn_eq",1)
        CALL h5_dump_rdum("tn",nfour,1)
      ENDIF
!-----------------------------------------------------------------------
!     write current, electric field and psi_eq.
!-----------------------------------------------------------------------
      IF (dump_ja) THEN
        CALL h5_dump_rdum("ja",nfour,3)
      ENDIF
      IF (dump_eexp) THEN
        CALL h5_dump_rdum("ee",nfour,3)
        CALL h5_dump_rdum_2D("eq",3)
      ENDIF
      IF (dump_q) THEN
        CALL h5_dump_rdum("qi",nfour,3)
        IF (separate_pe) CALL h5_dump_rdum("qe",nfour,3)
      ENDIF
      IF (dump_divpi) THEN
        CALL h5_dump_rdum("divpi",nfour,3)
      ENDIF
      IF (dump_psi_eq) THEN
        CALL h5_dump_rdum_2D("psi_eq",1)
      ENDIF
      IF (dump_fsa_beq2) THEN
        CALL h5_dump_rdum_2D("fsabeq2",1)
      ENDIF
!-----------------------------------------------------------------------
!     write edge coordinates if lagr_edge_rmp
!-----------------------------------------------------------------------
      IF (lagr_edge_rmp .AND. block_sizes(6,id)>0) THEN
        h5in%mesh="mesh-structured"
        h5in%vsAxisLabels="R,Z"
        CALL h5_dump_rdum_edge_1D("rzedge",2)
!-----------------------------------------------------------------------
!     write RMP magnetic fields.
!-----------------------------------------------------------------------
        h5in%mesh="/rblocks/"//TRIM(blid)//"/rzedge"//TRIM(blid)
        CALL h5_dump_rdum_edge("brmp",nfour,3)
      ENDIF
!-----------------------------------------------------------------------
!     write ponderomotive force.
!-----------------------------------------------------------------------
      IF (use_fp) CALL h5_dump_rdum("fp",nfour,3)
!-----------------------------------------------------------------------
!     terminate
!-----------------------------------------------------------------------
      h5in%mesh=" "
      h5in%vsAxisLabels=" "
      h5in%vsIndexOrder=" "
      h5in%vsTimeGroup=" "
      CALL close_group(TRIM(blid),gid,h5err)
      RETURN

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 6. h5_dump_rdum_2D
!     dump a dummy block for parallel I/O
!-----------------------------------------------------------------------
      SUBROUTINE h5_dump_rdum_2D(fname,nqty)
      USE io
      USE input, ONLY: poly_degree
      USE pardata, ONLY: block_sizes
      CHARACTER(*), INTENT(IN) :: fname
      INTEGER(i4), INTENT(IN) :: nqty
      INTEGER(i4) :: dims(3),ntmp

      dims(1)=nqty
      dims(2:3)=block_sizes(3:4,id)*poly_degree+1_i4
      h5in%vsMD=TRIM(fname)
      CALL dump_h5_rl_dum(gid,TRIM(fname)//TRIM(blid),                  &
     &                    dims,h5in,h5err)
      h5in%vsMD=" "

      RETURN
      END SUBROUTINE h5_dump_rdum_2D
!-----------------------------------------------------------------------
!     subprogram 7. h5_dump_rdum
!     dump a dummy block for parallel I/O
!-----------------------------------------------------------------------
      SUBROUTINE h5_dump_rdum(fname,nfour,nqty)
      USE io
      USE input, ONLY: poly_degree
      USE pardata, ONLY: block_sizes
      CHARACTER(*), INTENT(IN) :: fname
      INTEGER(i4), INTENT(IN) :: nfour,nqty

      INTEGER(i4) :: dims(3),ntmp
      CHARACTER(64) :: mdname

      dims(1)=nfour*nqty
      dims(2:3)=block_sizes(3:4,id)*poly_degree+1_i4
      h5in%vsMD=TRIM(fname)
      mdname=TRIM(h5in%vsMD)
      h5in%vsMD="Re_"//mdname
      CALL dump_h5_rl_dum(gid,"re"//TRIM(fname)//TRIM(blid),            &
     &                    dims,h5in,h5err)
      h5in%vsMD="Im_"//mdname
      CALL dump_h5_rl_dum(gid,"im"//TRIM(fname)//TRIM(blid),            &
     &                    dims,h5in,h5err)
      h5in%vsMD=" "

      RETURN
      END SUBROUTINE h5_dump_rdum
!-----------------------------------------------------------------------
!     subprogram 8. h5_dump_rdum_edge_1D
!     dump a dummy lagr_edge_1D for parallel I/O
!-----------------------------------------------------------------------
      SUBROUTINE h5_dump_rdum_edge_1D(fname,nqty)
      USE io
      USE input, ONLY: poly_degree
      USE pardata, ONLY: block_sizes
      CHARACTER(*), INTENT(IN) :: fname
      INTEGER(i4), INTENT(IN) :: nqty
      INTEGER(i4) :: dims(2),ntmp

      dims(1)=nqty
      dims(2)=block_sizes(6,id)*poly_degree+1_i4
      h5in%vsMD=TRIM(fname)
      CALL dump_h5_rl_dum(gid,TRIM(fname)//TRIM(blid),                  &
     &                    dims,h5in,h5err)
      h5in%vsMD=" "

      RETURN
      END SUBROUTINE h5_dump_rdum_edge_1D
!-----------------------------------------------------------------------
!     subprogram 9. h5_dump_rdum_edge
!     dump a dummy lagr_edge for parallel I/O
!-----------------------------------------------------------------------
      SUBROUTINE h5_dump_rdum_edge(fname,nfour,nqty)
      USE io
      USE input, ONLY: poly_degree
      USE pardata, ONLY: block_sizes
      CHARACTER(*), INTENT(IN) :: fname
      INTEGER(i4), INTENT(IN) :: nfour,nqty

      INTEGER(i4) :: dims(2),ntmp
      CHARACTER(64) :: mdname

      dims(1)=nfour*nqty
      dims(2)=block_sizes(6,id)*poly_degree+1_i4
      h5in%vsMD=TRIM(fname)
      mdname=TRIM(h5in%vsMD)
      h5in%vsMD="Re_"//mdname
      CALL dump_h5_rl_dum(gid,"re"//TRIM(fname)//TRIM(blid),            &
     &                    dims,h5in,h5err)
      h5in%vsMD="Im_"//mdname
      CALL dump_h5_rl_dum(gid,"im"//TRIM(fname)//TRIM(blid),            &
     &                    dims,h5in,h5err)
      h5in%vsMD=" "

      RETURN
      END SUBROUTINE h5_dump_rdum_edge

      END SUBROUTINE h5_dump_rbdum
!-----------------------------------------------------------------------
!     subprogram 10. h5_read_rblock.
!     reads rblock data from h5 dump file.
!-----------------------------------------------------------------------
      SUBROUTINE h5_read_rblock(rb,id,rbid)
      USE input, ONLY:phdump,nF,zimp,mx,my,F_dvac,dump_psi_eq,          &
     &                nxbl,nybl,closure_model,neutral_model,            &
     &                ph_pd,phqty,dump_eexp,dump_ja,dump_divpi,dump_q,  &
     &                dump_fsa_beq2,separate_pe,ds_nqty,lagr_edge_rmp,  &
     &                use_fp
      USE global, ONLY: nmodes_total
      USE io

      TYPE (rblock_type), INTENT(INOUT) :: rb
      INTEGER(i4), INTENT(IN):: id
      INTEGER(HID_T), INTENT(IN) :: rbid

      CHARACTER(4) :: blid
      INTEGER(i4) :: degenerate,pd,nfour
      INTEGER(HID_T) :: gid
!-----------------------------------------------------------------------
!     setup rblock group.
!-----------------------------------------------------------------------
      WRITE(blid,fmt='(i4.4)') id
      CALL open_group(rbid,TRIM(blid),gid,h5err)
      rb%name='rblock'//TRIM(blid)
!-----------------------------------------------------------------------
!     read block descriptors.
!-----------------------------------------------------------------------
      CALL read_attribute(gid,"id",rb%id,h5in,h5err)
      CALL read_attribute(gid,"mx",rb%mx,h5in,h5err)
      CALL read_attribute(gid,"my",rb%my,h5in,h5err)
      IF (attr_exists(gid,".","mxe",h5err)) THEN
        CALL read_attribute(gid,"mxe",rb%mxe,h5in,h5err)
      ELSE
        rb%mxe=0
      ENDIF
      CALL read_attribute(gid,"nfour",nfour,h5in,h5err)
      CALL read_attribute(gid,"poly_degree",pd,h5in,h5err)
      CALL read_attribute(gid,"degenerate",degenerate,h5in,h5err)
      IF (degenerate>0) THEN
        rb%degenerate=.true.
      ELSE
        rb%degenerate=.false.
      ENDIF
!-----------------------------------------------------------------------
!     information about the global arrays from input.f
!-----------------------------------------------------------------------
      rb%mx_total=mx
      rb%my_total=my
      rb%nxbl_total=nxbl
      rb%nybl_total=nybl
!-----------------------------------------------------------------------
!     read coordinates.
!-----------------------------------------------------------------------
      CALL h5_read_lagr_quad_2D(rb%rz,rb%mx,rb%my,pd,2,                 &
     &                          'rz',gid,blid)
!-----------------------------------------------------------------------
!     read magnetic fields.
!-----------------------------------------------------------------------
      CALL h5_read_lagr_quad_2D(rb%be_eq,rb%mx,rb%my,pd,3,              &
     &                          'bq',gid,blid)
      CALL h5_read_lagr_quad(rb%be,rb%mx,rb%my,pd,nfour,3,              &
     &                       'be',gid,blid)
!-----------------------------------------------------------------------
!     read equilibrium current density.
!-----------------------------------------------------------------------
      CALL h5_read_lagr_quad_2D(rb%ja_eq,rb%mx,rb%my,pd,3,              &
     &                          'jq',gid,blid)
!-----------------------------------------------------------------------
!     read fluid velocities.
!-----------------------------------------------------------------------
      CALL h5_read_lagr_quad_2D(rb%ve_eq,rb%mx,rb%my,pd,3,              &
     &                          'vq',gid,blid)
      CALL h5_read_lagr_quad(rb%ve,rb%mx,rb%my,pd,nfour,3,              &
     &                       've',gid,blid)
!-----------------------------------------------------------------------
!     read pressures.
!-----------------------------------------------------------------------
      CALL h5_read_lagr_quad_2D(rb%pres_eq,rb%mx,rb%my,pd,1,            &
     &                          'prq',gid,blid)
      CALL h5_read_lagr_quad(rb%pres,rb%mx,rb%my,pd,nfour,1,            &
     &                       'pr',gid,blid)
      CALL h5_read_lagr_quad_2D(rb%prese_eq,rb%mx,rb%my,pd,1,           &
     &                          'peq',gid,blid)
      CALL h5_read_lagr_quad(rb%prese,rb%mx,rb%my,pd,nfour,1,           &
     &                       'pe',gid,blid)
!-----------------------------------------------------------------------
!     read number densities.
!-----------------------------------------------------------------------
      CALL h5_read_lagr_quad_2D(rb%nd_eq,rb%mx,rb%my,pd,1,              &
     &                          'nq',gid,blid)
      CALL h5_read_lagr_quad(rb%nd,rb%mx,rb%my,pd,nfour,1,              &
     &                       'nd',gid,blid)
!-----------------------------------------------------------------------
!     read diffusivity shape factor and material concentration.
!-----------------------------------------------------------------------
      CALL h5_read_lagr_quad_2D(rb%diff_shape,rb%mx,rb%my,pd,ds_nqty,   &
     &                          'diff',gid,blid)
      CALL h5_read_lagr_quad(rb%conc,rb%mx,rb%my,pd,nfour,1,            &
     &                       'conc',gid,blid) 
!-----------------------------------------------------------------------
!     read temperatures.
!-----------------------------------------------------------------------
      CALL h5_read_lagr_quad(rb%tele,rb%mx,rb%my,pd,nfour,1,            &
     &                       'te',gid,blid)
      CALL h5_read_lagr_quad(rb%tion,rb%mx,rb%my,pd,nfour,1,            &
     &                       'ti',gid,blid)
!----------------------------------------------------------------------
!     read hot particle pressure tensor or current.
!----------------------------------------------------------------------
      IF (phdump) THEN
        IF (ph_pd/=pd.AND.ph_pd==1) THEN
          CALL h5_read_lagr_quad_2D(rb%rz_ph,rb%mx,rb%my,ph_pd,2,       &
     &                           'rz_ph',gid,blid)
        ENDIF
        CALL h5_read_lagr_quad(rb%phot,rb%mx,rb%my,ph_pd,nfour,phqty,   &
     &                         'phot',gid,blid)
      ENDIF
!----------------------------------------------------------------------
!     read kinetic distortions and closures.
!----------------------------------------------------------------------
      IF (nF(1)>0) THEN
        CALL h5_read_lagr_quad_2D(rb%Fele_eq,rb%mx,rb%my,pd,nF(1),      &
     &                            'Fele_eq',gid,blid)
        CALL h5_read_lagr_quad(rb%Fele,rb%mx,rb%my,pd,nfour,nF(1),      &
     &                         'Fele',gid,blid)
        CALL h5_read_lagr_quad(rb%qpe,rb%mx,rb%my,pd,nfour,1,           &
     &                         'qpe',gid,blid)
        CALL h5_read_lagr_quad(rb%ppe,rb%mx,rb%my,pd,nfour,1,           &
     &                         'ppe',gid,blid)
        CALL h5_read_lagr_quad(rb%rpe,rb%mx,rb%my,pd,nfour,1,           &
     &                         'rpe',gid,blid)
        CALL h5_read_lagr_quad(rb%gei,rb%mx,rb%my,pd,nfour,1,           &
     &                         'gei',gid,blid)
      ENDIF
      IF (nF(2)>0) THEN
        CALL h5_read_lagr_quad_2D(rb%Fion_eq,rb%mx,rb%my,pd,nF(2),      &
     &                            'Fion_eq',gid,blid)
        CALL h5_read_lagr_quad(rb%Fion,rb%mx,rb%my,pd,nfour,nF(2),      &
     &                         'Fion',gid,blid)
        CALL h5_read_lagr_quad(rb%qpi,rb%mx,rb%my,pd,nfour,1,           &
     &                         'qpi',gid,blid)
        CALL h5_read_lagr_quad(rb%ppi,rb%mx,rb%my,pd,nfour,1,           &
     &                         'ppi',gid,blid)
        CALL h5_read_lagr_quad(rb%rpi,rb%mx,rb%my,pd,nfour,1,           &
     &                         'rpi',gid,blid)
      ENDIF
      IF (nF(3)>0) THEN
        CALL h5_read_lagr_quad_2D(rb%Fhot_eq,rb%mx,rb%my,pd,nF(3),      &
     &                            'Fhot_eq',gid,blid)
        CALL h5_read_lagr_quad(rb%Fhot,rb%mx,rb%my,pd,nfour,nF(3),      &
     &                         'Fhot',gid,blid)
        CALL h5_read_lagr_quad_2D(rb%nhot_eq,rb%mx,rb%my,ph_pd,1,       &
     &                            'nhot_eq',gid,blid)
        CALL h5_read_lagr_quad_2D(rb%thot_eq,rb%mx,rb%my,pd,1,          &
     &                            'thot_eq',gid,blid)
      ENDIF
      IF (F_dvac>0) THEN
        CALL h5_read_lagr_quad_2D(rb%F_diff_shape,rb%mx,rb%my,pd,1,     &
     &                            'F_diff',gid,blid)
      ENDIF
!----------------------------------------------------------------------
      IF (zimp/=0) THEN
        CALL h5_read_lagr_quad(rb%nZzt,rb%mx,rb%my,pd,nfour,1,          &
     &                         'nZzt',gid,blid)
!        CALL h5_read_lagr_quad(rb%nZnt,rb%mx,rb%my,pd,nfour,1,          &
!     &                         'nZnt',gid,blid)
!        CALL h5_read_lagr_quad(rb%nZrt,rb%mx,rb%my,pd,nfour,1,          &
!     &                         'nZrt',gid,blid)
!        CALL h5_read_lagr_quad(rb%nZTt,rb%mx,rb%my,pd,nfour,1,          &
!     &                         'nZTt',gid,blid)
        CALL h5_read_lagr_quad(rb%nZel,rb%mx,rb%my,pd,nfour,1,          &
     &                         'nZel',gid,blid)
        CALL h5_read_lagr_quad(rb%nz,rb%mx,rb%my,pd,nfour,              &
     &                         INT(zimp+1,i4),'nz',gid,blid)
        CALL h5_read_lagr_quad(rb%nZeff,rb%mx,rb%my,pd,nfour,1,         &
     &                         'nZeff',gid,blid)
        CALL h5_read_lagr_quad(rb%qlosd,rb%mx,rb%my,pd,nfour,1,         &
     &                         'qlosd',gid,blid)
        CALL h5_read_lagr_quad(rb%qlosb,rb%mx,rb%my,pd,nfour,1,         &
     &                         'qlosb',gid,blid)
        CALL h5_read_lagr_quad(rb%qlosr,rb%mx,rb%my,pd,nfour,1,         &
     &                         'qlosr',gid,blid)
        CALL h5_read_lagr_quad(rb%qlosl,rb%mx,rb%my,pd,nfour,1,         &
     &                         'qlosl',gid,blid)
        CALL h5_read_lagr_quad(rb%qlosi,rb%mx,rb%my,pd,nfour,1,         &
     &                         'qlosi',gid,blid)
        CALL h5_read_lagr_quad(rb%ndn,rb%mx,rb%my,pd,nfour,1,           &
     &                         'ndn',gid,blid)
      ENDIF
!-----------------------------------------------------------------------
!     read neutrals data.
!-----------------------------------------------------------------------
      IF (neutral_model/='none') THEN
        CALL h5_read_lagr_quad_2D(rb%ndn_eq,rb%mx,rb%my,pd,1,           &
     &                            'ndn_eq',gid,blid)
        CALL h5_read_lagr_quad(rb%ndn,rb%mx,rb%my,pd,nfour,1,           &
     &                         'ndn',gid,blid)
        CALL h5_read_lagr_quad_2D(rb%vn_eq,rb%mx,rb%my,pd,3,            &
     &                            'vn_eq',gid,blid)
        CALL h5_read_lagr_quad(rb%vn,rb%mx,rb%my,pd,nfour,3,            &
     &                         'vn',gid,blid)
        CALL h5_read_lagr_quad_2D(rb%tn_eq,rb%mx,rb%my,pd,1,            &
     &                            'tn_eq',gid,blid)
        CALL h5_read_lagr_quad(rb%tn,rb%mx,rb%my,pd,nfour,1,            &
     &                         'tn',gid,blid)
      ENDIF
!-----------------------------------------------------------------------
!     read ja data.
!-----------------------------------------------------------------------
      IF (dump_ja) THEN 
        CALL h5_read_lagr_quad(rb%ja,rb%mx,rb%my,pd,nfour,3,            &
     &                         'ja',gid,blid)
      ENDIF
!-----------------------------------------------------------------------
!     read eexp data.
!-----------------------------------------------------------------------
      IF (dump_eexp) THEN
        CALL h5_read_lagr_quad_2D(rb%eef_eq,rb%mx,rb%my,pd,3,           &
     &                            'eq',gid,blid)
        CALL h5_read_lagr_quad(rb%eef,rb%mx,rb%my,pd,nfour,3,           &
     &                         'ee',gid,blid)
      ENDIF
!-----------------------------------------------------------------------
!     read divpi/q
!-----------------------------------------------------------------------
      IF (dump_divpi) THEN 
        CALL h5_read_lagr_quad(rb%divpi,rb%mx,rb%my,pd,nfour,3,         &
     &                         'divpi',gid,blid)
      ENDIF
      IF (dump_q) THEN
        CALL h5_read_lagr_quad(rb%qi,rb%mx,rb%my,pd,nfour,3,            &
     &                         'qi',gid,blid)
        IF (separate_pe)                                                &
     &     CALL h5_read_lagr_quad(rb%qe,rb%mx,rb%my,pd,nfour,3,         &
     &                            'qe',gid,blid)
      ENDIF
!-----------------------------------------------------------------------
!     read psi_eq data.
!-----------------------------------------------------------------------
      IF (dump_psi_eq) THEN
        CALL h5_read_lagr_quad_2D(rb%psi_eq,rb%mx,rb%my,pd,1,           &
     &                         'psi_eq',gid,blid)
      ENDIF
!-----------------------------------------------------------------------
!     read fsa_beq2 data.
!-----------------------------------------------------------------------
      IF (dump_fsa_beq2) THEN
        CALL h5_read_lagr_quad_2D(rb%fsa_beq2,rb%mx,rb%my,pd,1,         &
     &                         'fsabeq2',gid,blid)
      ENDIF
!-----------------------------------------------------------------------
!    read edge coordinates if lagr_edge_rmp
!-----------------------------------------------------------------------
      IF (lagr_edge_rmp .AND. rb%mxe>0) THEN
        CALL h5_read_lagr_edge_1D(rb%rzedge,rb%mxe,pd,2,"rzedge",       &
     &                            gid,blid)
!-----------------------------------------------------------------------
!     read RMP magnetic fields.
!-----------------------------------------------------------------------
        CALL h5_read_lagr_edge(rb%b_rmp,rb%mxe,pd,nfour,3,              &
     &                         "brmp",gid,blid)
      ENDIF
!-----------------------------------------------------------------------
!     read ponderomotive force.
!-----------------------------------------------------------------------
      IF (use_fp) THEN 
        CALL h5_read_lagr_quad(rb%fp,rb%mx,rb%my,pd,nfour,3,            &
     &                         'fp',gid,blid)
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      CALL close_group(TRIM(blid),gid,h5err)
      RETURN
      END SUBROUTINE h5_read_rblock
#endif /* HAVE_FC_HDF5 */
!-----------------------------------------------------------------------
!     close module
!-----------------------------------------------------------------------
      END MODULE rblock_type_mod
