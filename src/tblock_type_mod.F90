!-----------------------------------------------------------------------
!     $Id: tblock_type_mod.f90 7623 2022-01-18 03:55:24Z tbechtel $
!     contains a module that defines tblock types.
!     subprogram 2. dump_write_tblock.
!     subprogram 3. dump_read_tblock.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.
!-----------------------------------------------------------------------
      MODULE tblock_type_mod
      USE local
      USE tri_linear
      IMPLICIT NONE

!-----------------------------------------------------------------------
!     types for saving data at quadrature points
!-----------------------------------------------------------------------
      TYPE :: tb_real_qp_type
        REAL(r8), DIMENSION(:,:,:), POINTER :: qpf,qpfr,qpfz
      END TYPE tb_real_qp_type

      TYPE :: tb_real_big_qp_type
        REAL(r8), DIMENSION(:,:,:,:), POINTER :: qpf,qpfr,qpfz
      END TYPE tb_real_big_qp_type

      TYPE :: tb_comp_qp_type
        COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: qpf,qpfr,qpfz
      END TYPE tb_comp_qp_type

      TYPE :: tb_comp_big_qp_type
        COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: qpf,qpfr,qpfz
      END TYPE tb_comp_big_qp_type
!-----------------------------------------------------------------------
!     data saved in tblock types
!-----------------------------------------------------------------------
      TYPE :: tblock_type
        CHARACTER(64) :: name
        INTEGER(i4) :: id
        INTEGER(i4) :: mvert,mcell
        INTEGER(i4) :: ng
        TYPE(tri_linear_geom_type) :: tgeom

        TYPE(tri_linear_2D_type) :: be_eq,ja_eq,ve_eq,pres_eq,prese_eq, &
     &            nd_eq,diff_shape,pflux,be_n0,pres_n0,nd_n0,           &
     &            e_applied,q_applied,si_nl_pres,tele_eq,psi_eq,        &
     &            tion_eq,nd_eq2,rwork1,rwork2,rwork3,ve_n0,ti_n0,te_n0,&
     &            Fele_eq,Fion_eq,Fhot_eq,F_diff_shape,                 &
     &            phot_eq,vhot_eq,G_eq,fsa_beq2
        TYPE(tri_linear_2D_type) :: be_sv,ja_sv,ve_sv,pres_sv,prese_sv, &
     &            tele_sv,tion_sv
        TYPE(tri_linear_type) :: be,ja,ve,pres,prese,nd,tele,tion,      &
     &            conc,work1,work2,work3,work4,work5,work6,work7,work8  &
     &            ,qpe,dqpe,qpi,dqpi,ppe,ppi,rpe,rpi,gei,phot,phot1     &
     &            ,nimp,nion,qloss,w6v1,w6v2                            &
     &            ,Fele,Fion,Fhot,Fwork3,vee                            &
     &            ,Fework1,Fework2,Fiwork1,Fiwork2,Fhwork1,Fhwork2      &
     &            ,poten,dppe,dppi
        TYPE(tri_linear_type) :: auxv,auxb,mwork1,mwork2,mwork3,mwork4

        TYPE(tb_real_qp_type) :: qbe_eq,qja_eq,qve_eq,qpres_eq,         &
     &            qprese_eq,qnd_eq,qbe_n0,qpres_n0,qnd_n0,qndn_eq,      &
     &            qe_applied,qq_applied,qsi_nl_pres,qndn_tot,qndn_n0,   &
     &            qbb,qtele_eq,qtion_eq,qrwork1,qbe_tot,                &
     &            qtion_tot,qtele_tot,qvn_eq,qvn_n0,qtn_eq,qvn_tot,     &
     &            qtn_tot,qtn_n0,qpn_eq,qpn_n0,qgrdvn_tot,qgrdvn_eq,    &
     &            qnd_tot,qelecd_phi,qelecd_n0,qelecd_eq,qdiff_shape,   &
     &            qkappli_phi,qkapple_phi,qkappli_n0,qkapple_n0,        &
     &            qkaprpi_phi,qkaprpe_phi,                              &
     &            qkaprpi_n0,qkaprpe_n0,qte_b2,qti_b2,qbcrgte,qbcrgti,  &
     &            qve_n0,qve_tot,qgrdv,qja_tot,qdvv_eq,qeq_force,       &
     &            qdart,qupw_phi,qupw_n0,qvv,qupti_phi,qupti_n0,        &
     &            qupte_phi,qupte_n0,qti_n0,qte_n0,qgrdveq,qpi_veq,     &
     &            qpi_pareq,qpi_gyreq,qpin_eq,                          &
     &            qbe_tot_cel,qkappli_phi_cel,qkapple_phi_cel,          &
     &            qsigmav_phi,qsigmav_n0,qsvion_phi,qsvrec_ph,          &
     &            qsion_tot,qsion_sym,qsrec_tot,qsrec_sym,              &
     &            qsnn_eq,qscx_tot,qscx_sym,qsion_eq,qsrec_eq,qscx_eq,  &
     &            qstn_eq,qsti_eq,qsvn_eq,qsve_eq,                      &
     &            qFele_eq,qFion_eq,qFhot_eq,qF_diff_shape,             &
     &            qgrdb,qpr_tot,qgrdp,qivn_eq,qrvn_eq,qcvn_eq,qcxf_eq,  &
     &            qdvvn_eq,qfsa_beq2,qneodiff
        TYPE(tb_comp_qp_type) :: qbe,qja,qve,qpres,qprese,qnd,qtele,    &
     &            qtion,qconc,qvisc,qbe_off,qja_off,qve_off,qnd_off,    &
     &            qwork1,qwork2,qwork3,qwork4                           &
     &            ,qkapple_off,qkaprpe_off,qkappli_off,qkaprpi_off      &
     &            ,qqpe,qdqpe,qqpi,qdqpi,qppe,qppi,qpreseh,qndh,qphot   &
     &            ,qe_surface,qnimp,qnion,qndn,qvn,qtn,                 &
     &            qqloss,qrho                                           &
     &           ,qFele,qFion,qFhot,qFwork3,qrpe                        &
     &           ,qvee,qrpi,qbb1,qneut_vheat,qgei,qdppe,qdppi

        TYPE(tb_real_big_qp_type) :: qcee_n0,qcei_n0,qcii_n0,qcie_n0
        TYPE(tb_comp_big_qp_type) :: qFe_mom,qFi_mom

        ! This is used by nimset
        TYPE(tb_real_qp_type) :: qrz,qpsi_eq
        TYPE(tb_comp_qp_type) :: qpsi
        ! Sources - and things for sources.
        TYPE(tb_real_qp_type) :: qe_src_n0,qqi_src_n0,qqe_src_n0,       &
     &                           qnd_src_n0,qve_src_n0
        TYPE(tb_comp_qp_type) :: qe_source,qqi_source,qqe_source,       &
     &                           qnd_source,qve_source
        TYPE(tb_real_qp_type) :: qja_sv,qve_sv,qtele_sv,qtion_sv,qnd_sv,&
     &                           qrho_sym,qrho_tot
        TYPE(tb_real_qp_type) :: qrelax_mask

        REAL(r8), DIMENSION(:), POINTER :: wg
        REAL(r8), DIMENSION(:), POINTER :: cell_vol,r_cent
      END TYPE tblock_type
!-----------------------------------------------------------------------
!     type used for defining the mass matrices.
!-----------------------------------------------------------------------
      TYPE :: matrix_element_type1
        REAL(r8), DIMENSION(:), POINTER :: element
      END TYPE matrix_element_type1

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 2. dump_write_tblock.
!     writes tblock data to dump file.
!-----------------------------------------------------------------------
      SUBROUTINE dump_write_tblock(tb)
      USE input, ONLY:nF,dump_fsa_beq2
      USE io 

      TYPE (tblock_type), INTENT(IN) :: tb
      INTEGER(i4) :: ivert,nnbr
!-----------------------------------------------------------------------
!     write geometry.
!-----------------------------------------------------------------------
      WRITE(dump_unit) REAL(tb%id,r8)
      WRITE(dump_unit) REAL(tb%tgeom%mvert,r8)
      WRITE(dump_unit) REAL(tb%tgeom%mcell,r8)
      WRITE(dump_unit)tb%tgeom%xs
      WRITE(dump_unit)tb%tgeom%ys
      WRITE(dump_unit) REAL(tb%tgeom%vertex,r8)
      DO ivert=0,tb%tgeom%mvert
         nnbr=SIZE(tb%tgeom%neighbor(ivert)%vertex)-1
         WRITE(dump_unit) REAL(nnbr,r8)
         WRITE(dump_unit) REAL(tb%tgeom%neighbor(ivert)%vertex,r8)
      ENDDO
!-----------------------------------------------------------------------
!     write magnetic fields.
!-----------------------------------------------------------------------
      CALL dump_write_tri_linear_2D(tb%be_eq)
      CALL dump_write_tri_linear(tb%be)
!-----------------------------------------------------------------------
!     write equilibrium current density.
!-----------------------------------------------------------------------
      CALL dump_write_tri_linear_2D(tb%ja_eq)
!-----------------------------------------------------------------------
!     write fluid velocities.
!-----------------------------------------------------------------------
      CALL dump_write_tri_linear_2D(tb%ve_eq)
      CALL dump_write_tri_linear(tb%ve)
!-----------------------------------------------------------------------
!     write pressures.
!-----------------------------------------------------------------------
      CALL dump_write_tri_linear_2D(tb%pres_eq)
      CALL dump_write_tri_linear(tb%pres)
      CALL dump_write_tri_linear_2D(tb%prese_eq)
      CALL dump_write_tri_linear(tb%prese)
!-----------------------------------------------------------------------
!     write number densities.
!-----------------------------------------------------------------------
      CALL dump_write_tri_linear_2D(tb%nd_eq)
      CALL dump_write_tri_linear(tb%nd)
!-----------------------------------------------------------------------
!     write diffusivity shape factor and material concentration.
!-----------------------------------------------------------------------
      CALL dump_write_tri_linear_2D(tb%diff_shape)
      CALL dump_write_tri_linear(tb%conc)
!-----------------------------------------------------------------------
!     write temperatures.
!-----------------------------------------------------------------------
      CALL dump_write_tri_linear(tb%tele)
      CALL dump_write_tri_linear(tb%tion)
!-----------------------------------------------------------------------
!     write kinetic distortions and closures.                                        
!-----------------------------------------------------------------------
      IF (nF(1)>0) THEN 
        CALL dump_write_tri_linear_2D(tb%Fele_eq) 
        CALL dump_write_tri_linear(tb%Fele) 
        CALL dump_write_tri_linear(tb%qpe) 
        CALL dump_write_tri_linear(tb%ppe) 
        CALL dump_write_tri_linear(tb%rpe) 
        CALL dump_write_tri_linear(tb%gei) 
      ENDIF
      IF (nF(2)>0) THEN 
        CALL dump_write_tri_linear_2D(tb%Fion_eq) 
        CALL dump_write_tri_linear(tb%Fion) 
        CALL dump_write_tri_linear(tb%qpi) 
        CALL dump_write_tri_linear(tb%ppi) 
        CALL dump_write_tri_linear(tb%rpi) 
      ENDIF
      IF (nF(3)>0) THEN 
        CALL dump_write_tri_linear_2D(tb%Fhot_eq) 
        CALL dump_write_tri_linear(tb%Fhot) 
      ENDIF
!-----------------------------------------------------------------------
!     write fsa beq2.
!-----------------------------------------------------------------------
      IF (dump_fsa_beq2) CALL dump_write_tri_linear_2D(tb%fsa_beq2) 
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dump_write_tblock
!-----------------------------------------------------------------------
!     subprogram 3. dump_read_tblock.
!     reads tblock data to dump file.
!-----------------------------------------------------------------------
      SUBROUTINE dump_read_tblock(tb)
      USE input, ONLY:nF,dump_read_F, dump_fsa_beq2
      USE io 

      TYPE (tblock_type), INTENT(OUT) :: tb
      INTEGER(i4) :: ivert,nnbr,nmodes
      REAL(r8) :: iread
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ia1read
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: ia2read
!-----------------------------------------------------------------------
!     read geometry.
!-----------------------------------------------------------------------
      tb%name='tblock'
      READ(rstrt_unit) iread
      tb%id=NINT(iread)
      READ(rstrt_unit) iread
      tb%tgeom%mvert=NINT(iread)
      READ(rstrt_unit) iread
      tb%tgeom%mcell=NINT(iread)
      CALL tri_linear_geom_alloc(tb%tgeom,tb%tgeom%mvert,tb%tgeom%mcell)
      READ(rstrt_unit) tb%tgeom%xs
      READ(rstrt_unit) tb%tgeom%ys
      ALLOCATE(ia2read(tb%tgeom%mcell,3))
      READ(rstrt_unit) ia2read
      tb%tgeom%vertex=NINT(ia2read)
      DEALLOCATE(ia2read)
      DO ivert=0,tb%tgeom%mvert
         READ(rstrt_unit) iread
         nnbr=NINT(iread)
         ALLOCATE(tb%tgeom%neighbor(ivert)%vertex(0:nnbr))
         ALLOCATE(ia1read(0:nnbr))
         READ(rstrt_unit) ia1read
         tb%tgeom%neighbor(ivert)%vertex=NINT(ia1read)
         DEALLOCATE(ia1read)
      ENDDO
      tb%mvert=tb%tgeom%mvert
      tb%mcell=tb%tgeom%mcell
!-----------------------------------------------------------------------
!     read magnetic fields.
!-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%be_eq,tb%mvert,'bq',(/' be_eq'/))
      CALL dump_read_tri_linear(tb%be,tb%mvert,'be',(/'  be  '/))
!-----------------------------------------------------------------------
!     read equilibrium current density.
!-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%ja_eq,tb%mvert,'jq',(/' ja_eq'/))
!-----------------------------------------------------------------------
!     read fluid velocities.
!-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%ve_eq,tb%mvert,'vq',(/' ve_eq'/))
      CALL dump_read_tri_linear(tb%ve,tb%mvert,'ve',(/'  ve  '/))
!-----------------------------------------------------------------------
!     read pressures.
!-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%pres_eq,tb%mvert,'pq',            &
     &                             (/' pr_eq'/))
      CALL dump_read_tri_linear(tb%pres,tb%mvert,'pr',(/' pres '/))
      CALL dump_read_tri_linear_2D(tb%prese_eq,tb%mvert,'pq',           &
     &                             (/'pre_eq'/))
      CALL dump_read_tri_linear(tb%prese,tb%mvert,'pe',(/' prese'/))
!-----------------------------------------------------------------------
!     read number densities.
!-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%nd_eq,tb%mvert,'nq',(/' nd_eq'/))
      CALL dump_read_tri_linear(tb%nd,tb%mvert,'nd',(/'  nd  '/))
!-----------------------------------------------------------------------
!     read diffusivity shape factor and material concentration.
!-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%diff_shape,tb%mvert,'ds',         &
     &                             (/'dif_sh'/))
      CALL dump_read_tri_linear(tb%conc,tb%mvert,'co',(/' conc '/))
!-----------------------------------------------------------------------
!     read temperatures.
!-----------------------------------------------------------------------
      CALL dump_read_tri_linear(tb%tele,tb%mvert,'te',(/' tele '/))
      CALL dump_read_tri_linear(tb%tion,tb%mvert,'ti',(/' tion '/))
!-----------------------------------------------------------------------
!     read kinetic distortions and closures.                                         
!-----------------------------------------------------------------------
      IF (nF(1)>0.AND.dump_read_F(1)) THEN 
       CALL dump_read_tri_linear_2D(tb%Fele_eq,tb%mvert,'Feq',          &
     &                                               (/' Fe_eq'/)) 
       CALL dump_read_tri_linear(tb%Fele,tb%mvert,'Fele',(/' Fele '/)) 
       CALL dump_read_tri_linear(tb%qpe,tb%mvert,'qpe',(/' qpe  '/)) 
       CALL dump_read_tri_linear(tb%ppe,tb%mvert,'ppe',(/' ppe  '/)) 
       CALL dump_read_tri_linear(tb%rpe,tb%mvert,'rpe',(/' rpe  '/)) 
       CALL dump_read_tri_linear(tb%gei,tb%mvert,'gei',(/' gei  '/)) 
      ENDIF
      IF (nF(2)>0.AND.dump_read_F(2)) THEN 
       CALL dump_read_tri_linear_2D(tb%Fion_eq,tb%mvert,'Fiq',          &
     &                                               (/' Fi_eq'/)) 
       CALL dump_read_tri_linear(tb%Fion,tb%mvert,'Fion',(/' Fion '/)) 
       CALL dump_read_tri_linear(tb%qpi,tb%mvert,'qpi',(/' qpi  '/)) 
       CALL dump_read_tri_linear(tb%ppi,tb%mvert,'ppi',(/' ppi  '/)) 
       CALL dump_read_tri_linear(tb%rpi,tb%mvert,'rpi',(/' rpi  '/)) 
      ENDIF
      IF (nF(3)>0.AND.dump_read_F(3)) THEN 
       CALL dump_read_tri_linear_2D(tb%Fhot_eq,tb%mvert,'Fhq',          &
     &                                               (/' Fh_eq'/)) 
       CALL dump_read_tri_linear(tb%Fhot,tb%mvert,'Fhot',(/' Fhot '/)) 
      ENDIF
!-----------------------------------------------------------------------
!     read fsa_beq2.
!-----------------------------------------------------------------------
      IF (dump_fsa_beq2) THEN
        CALL dump_read_tri_linear_2D(tb%fsa_beq2,tb%mvert,'fsbeq2',     &
     &  (/'fsbeq2'/))
      ENDIF 
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dump_read_tblock
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      END MODULE tblock_type_mod
