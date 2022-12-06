!-----------------------------------------------------------------------
!     $Id: fields.f90 7623 2022-01-18 03:55:24Z tbechtel $
!     module containing block and block-connection structures; all grid
!     and spatially dependent variables are within.
!     note:  nbl_total = total # of blocks in problem (on all
!     processors), and nbl = # of blocks stored on this processor.
!     However, for nimset, only nbl is used--nimset does not run in
!     parallel.
!     nrbl (nrbl_total) is the number of rectangular blocks numbered
!     from 1 to nrbl, and the block of unstructured triangles are
!     indexed from nrbl+1 to nbl (nrbl_total+1 to nbl_total).
!-----------------------------------------------------------------------
      MODULE fields
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE vector_type_mod
      IMPLICIT NONE

      INTEGER(i4) :: nbl_total,nrbl_total
      INTEGER(i4) :: nbl,nrbl
      TYPE(rblock_type), DIMENSION(:), POINTER :: rb
      TYPE(tblock_type), DIMENSION(:), POINTER :: tb
!-----------------------------------------------------------------------
!     the following are pointers to nodal data.  the form is the same
!     for all blocks.
!-----------------------------------------------------------------------
      TYPE(cvector_type), DIMENSION(:), POINTER ::                      &
     &   be,ja,ve,pres,prese,nd,tele,tion,conc,ndn,vn,tn                &
     &  ,work1,work2,work3,work4,work5,work6,work7,work8                &
     &  ,qpe,qpi,ppe,ppi,rpe,rpi,e_surface,poten,phot,nimp              &
     &  ,qlosd,qlosb,qlosr,qlosl,qlosi,nzp,rhomn,auxb,auxv              &
     &  ,w6v1,w6v2,Fele,Fion,Fhot,Fwork3,nd2n                           &
     &  ,bmod,eef,curv,vee,qmfeme,qmfemi,ephi,ephip,gradB,gradp,jhot    &
     &  ,qi,qe,divpi,qDhFe,qDhFi,BdB1                                   &
     &  ,Fework1,Fework2,Fiwork1,Fiwork2,Fhwork1,Fhwork2,deef,twist_ln  &
     &  ,dqpe,dqpi,dppe,dppi,gei,Fework1_1,Fiwork1_1                    &
     &  ,nZzt,nZnt,nZrt,nZel,nZTt,nZeff
      TYPE(cvector_type), DIMENSION(:), POINTER ::                      &
     &   be_old,ja_old,ve_old                                           &
     &  ,p_old,pe_old,nd_old,te_old,ti_old,ndn_old,vn_old,tn_old        &
     &  ,nimp_old,nzp_old,Fe_old,Fi_old,Fh_old,phot_old,phot1           &
     &  ,jhot_old,vee_old                                               &
     &  ,nZzt_old,nZnt_old,nZrt_old,nZel_old,nZTt_old,nZeff_old
      TYPE(vector_type), DIMENSION(:), POINTER ::                       &
     &   be_eq,ja_eq,ve_eq,pres_eq,prese_eq,nd_eq,diff_shape,rzpt       &
     &  ,be_n0,pres_n0,nd_n0,rwork1,rwork2,rwork3,e_applied,q_applied   &
     &  ,si_nl_pres,tele_eq,tion_eq,nd_eq2,ve_n0,elecd_n0               &
     &  ,kappli_n0,kapple_n0,kaprpi_n0,kaprpe_n0,pflux,psi_eq           &
     &  ,rho_sym,rho_tot,ti_n0,te_n0,pmom_n0,minsol                     &
     &  ,Fele_eq,Fion_eq,Fhot_eq,eef_eq,F_diff_shape                    &
     &  ,ndn_eq,vn_eq,tn_eq,pn_eq,ndn_n0,vn_n0,tn_n0,pn_n0,twist_eq
!-----------------------------------------------------------------------
!     close module
!-----------------------------------------------------------------------
      END MODULE fields
