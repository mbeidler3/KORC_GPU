!-----------------------------------------------------------------------
!     Routines useful for doing the mapping problem. The key structure
!     is the rb_cel structure which is setup using the exchange
!     subroutines.  This puts all of the information into 1 block.  The
!     key routine then is rz_to_xy which finds xy given rz.  Once that
!     information is know, one can use the get_fields subroutines to
!     evaluate the information needed.
!
!     Definitions:
!           the direct  map is: x(R,Z), y(R,Z)
!           the inverse map is: R(x,y), Z(x,y)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.
!-----------------------------------------------------------------------
!     module map_mod
!     1. map_type_init
!     6. rz_to_xy_hard.
!     end module map_mod
!-----------------------------------------------------------------------
!     subprogram 0. map_mod
!     module containing the type declaration and default values for
!     nimfl input variables.
!-----------------------------------------------------------------------
      MODULE map_mod
      USE bicube
      USE global
      USE rblock_type_mod
      USE pardata
      USE lagr_quad_mod
      USE exchange_mod
      IMPLICIT NONE
!-----------------------------------------------------------------------
!     global cel data saved in rb types
!-----------------------------------------------------------------------
      TYPE :: cel_rblock_type
        CHARACTER(64) :: name
        INTEGER(i4) :: id
        INTEGER(i4) :: mx,my
        INTEGER(i4) :: ng,nge
        LOGICAL :: degenerate,r0block
        TYPE(lagr_quad_2D_type) :: rz,be_eq,ja_eq,ve_eq,pres_eq,        &
     &            prese_eq,nd_eq,tion_eq,tele_eq,Fion_eq,Fele_eq,       &
     &            Fhot_eq,psi_eq,ramos2,eef_eq,exb,piv,keflx,qi0,bdglnB,&
     &            partvar_eq,phot_eq,rwork2,diff_shape,fsa_beq2
        TYPE(lagr_quad_type) :: be,ja,ve,pres,prese,nd,tele,tion,       &
     &            qpe,qpi,ppe,ppi,phot,qpll,qprp,conc,ne,nimp,Fion,Fele,&
     &            rpe,rpi,gei,Fhot,eef,poten,ephip,partvar,ndn,nz,      &
     &            qlosd,qlosb,qlosr,qlosl,qlosi,fp,nimpw,bmod,nZeff

        REAL(r8), DIMENSION(:), POINTER :: xg,yg,wg
        REAL(r8), DIMENSION(:,:,:), POINTER :: dxdr,dxdz,dydr,dydz,     &
     &            bigr,wjac,jac2d,cell_int,cell_min,cell_max
        REAL(r8), DIMENSION(:,:,:,:), POINTER :: alpha,dalpdr,dalpdz,   &
     &                                           dalpdrc

      END TYPE cel_rblock_type
!-----------------------------------
!     Key storage type for the map
!-----------------------------------
      TYPE(cel_rblock_type), DIMENSION(:), POINTER :: rb_cel
      TYPE(rblock_type) :: rbrecv
      REAL(r8), PUBLIC :: rmin,rmax,zmin,zmax
      INTEGER(i4), DIMENSION(:,:,:,:), ALLOCATABLE :: mapl2g
      INTEGER(i4), DIMENSION(:,:,:), ALLOCATABLE :: mapg2l
      LOGICAL :: mapmod_initialized=.FALSE.
!-----------------------------------
!     Stuff for the spline map guess
!-----------------------------------
      LOGICAL :: spline_guess=.FALSE.
      INTEGER(i4) :: mdmap=15 ! Size of dmap_rz
      TYPE(bicube_type) :: dmap_rz

      INTEGER(i4), PRIVATE :: mimap,mext,mint,mxm,mym
      REAL(r8), PRIVATE :: rnn,znn,rnx,znx,rxn,zxn,rxx,zxx,rm,zm
      TYPE(bicube_type), PRIVATE :: imap_rz

      INTERFACE rz_to_xy
        MODULE PROCEDURE rz_to_xy_hard
      END INTERFACE

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 1.  map_type_init
!     Initialize the main rb_cel structure as well as arrays for data
!     communication (exchange_init)
!-----------------------------------------------------------------------
      SUBROUTINE map_type_init(b,j,v,p,pe,n,ti,te,e,qpe,qpi,ppe,ppi,    &
     &                         rpe,rpi,gei,phot,conc,alp,ne,nimp,bdglnB,&
     &                         Fion,Fele,Fhot,Fion_eq,Fele_eq,Fhot_eq,  &
     &                         be_eq,psi_eq,ramos2,flx,poten,ephip,     &
     &                         partvar_eq,ja_eq,ve_eq,p_eq,pe_eq,nd_eq, &
     &                         ti_eq,te_eq,rwork2,ndn,los,f,nimpw,bmod, &
     &                         diff_shape,fsa_beq2)
      USE input
      USE physdat
      IMPLICIT NONE

      CHARACTER*(*), INTENT(IN), OPTIONAL :: b,j,v,p,pe,n,ti,te,e,bdglnB
      CHARACTER*(*), INTENT(IN), OPTIONAL :: qpe,qpi,ppe,ppi,phot,conc
      CHARACTER*(*), INTENT(IN), OPTIONAL :: Fele,Fion,Fhot,rpe,rpi,gei
      CHARACTER*(*), INTENT(IN), OPTIONAL :: Fele_eq,Fion_eq,Fhot_eq
      CHARACTER*(*), INTENT(IN), OPTIONAL :: alp,be_eq,psi_eq,ephip,nimpw
      CHARACTER*(*), INTENT(IN), OPTIONAL :: ne,nimp,ramos2,flx,poten
      CHARACTER*(*), INTENT(IN), OPTIONAL :: partvar_eq,ja_eq,ve_eq
      CHARACTER*(*), INTENT(IN), OPTIONAL :: p_eq,pe_eq,nd_eq,ti_eq
      CHARACTER*(*), INTENT(IN), OPTIONAL :: te_eq,rwork2,ndn,los,f,bmod
      CHARACTER*(*), INTENT(IN), OPTIONAL :: diff_shape,fsa_beq2
      INTEGER(i4) :: ierr,ibl,nmd,pd,is,ii,ins,ine,nny,nnx
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xm,ym,xnode
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rzall
      LOGICAL :: degenerate
!-----------------------------------------------------------------------
!     INTERFACE for polynomial interpolation and for getting nodes
!-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE poly_nodes(n,x)
        USE poly_mod
        IMPLICIT NONE
        INTEGER(i4), INTENT(IN) :: n
        REAL(r8), DIMENSION(0:n), INTENT(OUT) :: x
        END SUBROUTINE poly_nodes
      END INTERFACE

      mapmod_initialized=.TRUE.
!-----------------------------------------------------------------------
!     Initialize the arrays used for communication exchange
!-----------------------------------------------------------------------
      CALL exchange_init
!-----------------------------------------------------------------------
!     Set up rb_cel structure
!-----------------------------------------------------------------------
      ALLOCATE(rb_cel(1))
      rb_cel(1)%name = rb(1)%name
      rb_cel(1)%id = 1
      rb_cel(1)%mx = mx
      rb_cel(1)%my = my

      ! Useful vars
      nmd=nmodes_total
      IF (ALLOCATED(rb(1)%rz%fsh)) THEN
        pd=SIZE(rb(1)%rz%fsh,2)+1
      ELSE
        pd=1
      ENDIF

      ! Set degeneracy
      IF(node==0) degenerate = rb(1)%degenerate
      IF (nprocs>1)                                                     &
     &  CALL mpi_bcast(degenerate,1,mpi_nim_logical,0,comm_nimrod,ierr)
      rb_cel(1)%degenerate = degenerate
      NULLIFY(rb_cel(1)%xg,rb_cel(1)%yg,rb_cel(1)%wg,rb_cel(1)%dxdr,    &
     &        rb_cel(1)%dxdz,rb_cel(1)%dydr,rb_cel(1)%dydz,             &
     &        rb_cel(1)%bigr,rb_cel(1)%wjac,rb_cel(1)%jac2d,            &
     &        rb_cel(1)%cell_int,rb_cel(1)%cell_min,rb_cel(1)%cell_max, &
     &        rb_cel(1)%alpha,rb_cel(1)%dalpdr,rb_cel(1)%dalpdz,        &
     &        rb_cel(1)%dalpdrc)
!-----------------------------------------------------------------------
!     Exchange rz since we always need it.
!-----------------------------------------------------------------------
      CALL lagr_quad_alloc(rb_cel(1)%rz,mx,my,2_i4,pd)
      DO ibl=1,nrbl
        CALL exch_lagr(rb(ibl)%rz,rb_cel(1)%rz,rb(ibl)%id)
      ENDDO
!-----------------------------------------------------------------------
!     find rz min and max
!-----------------------------------------------------------------------
      IF (ALLOCATED(rb_cel(1)%rz%fsv)) THEN
        rmin=MIN(MINVAL(rb_cel(1)%rz%fs(1,:,:)),                        &
     &           MIN(MINVAL(rb_cel(1)%rz%fsv(1,:,:,:)),                 &
     &               MINVAL(rb_cel(1)%rz%fsh(1,:,:,:))))
        rmax=MAX(MAXVAL(rb_cel(1)%rz%fs(1,:,:)),                        &
     &           MAX(MAXVAL(rb_cel(1)%rz%fsv(1,:,:,:)),                 &
     &               MAXVAL(rb_cel(1)%rz%fsh(1,:,:,:))))
        zmin=MIN(MINVAL(rb_cel(1)%rz%fs(2,:,:)),                        &
     &           MIN(MINVAL(rb_cel(1)%rz%fsv(2,:,:,:)),                 &
     &               MINVAL(rb_cel(1)%rz%fsh(2,:,:,:))))
        zmax=MAX(MAXVAL(rb_cel(1)%rz%fs(2,:,:)),                        &
     &           MAX(MAXVAL(rb_cel(1)%rz%fsv(2,:,:,:)),                 &
     &               MAXVAL(rb_cel(1)%rz%fsh(2,:,:,:))))
      ELSE
        rmin=MINVAL(rb_cel(1)%rz%fs(1,:,:))
        rmax=MAXVAL(rb_cel(1)%rz%fs(1,:,:))
        zmin=MINVAL(rb_cel(1)%rz%fs(2,:,:))
        zmax=MAXVAL(rb_cel(1)%rz%fs(2,:,:))
      ENDIF
!-----------------------------------------------------------------------
!     Other variables get set if requested
!-----------------------------------------------------------------------
      IF(PRESENT(b).AND.b=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%be,mx,my,3_i4,nmd,pd)
        CALL lagr_quad_alloc(rb_cel(1)%be_eq,mx,my,3_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%be_eq,rb_cel(1)%be_eq,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%be,rb_cel(1)%be,rb(ibl)%id)
        ENDDO
      ELSEIF(PRESENT(be_eq).AND.be_eq=='y') THEN
        CALL lagr_quad_alloc(rb_cel(1)%be_eq,mx,my,3_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%be_eq,rb_cel(1)%be_eq,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(psi_eq).AND.psi_eq=='y') THEN
        CALL lagr_quad_alloc(rb_cel(1)%psi_eq,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%psi_eq,rb_cel(1)%psi_eq,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(bdglnB).AND.bdglnB=='y') THEN
        CALL lagr_quad_alloc(rb_cel(1)%bdglnB,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%bdglnB,rb_cel(1)%bdglnB,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(j).AND.j=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%ja,mx,my,3_i4,nmd,pd)
        CALL lagr_quad_alloc(rb_cel(1)%ja_eq,mx,my,3_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%ja_eq,rb_cel(1)%ja_eq,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%ja,rb_cel(1)%ja,rb(ibl)%id)
        ENDDO
      ELSEIF(PRESENT(ja_eq).AND.ja_eq=='y') THEN
        CALL lagr_quad_alloc(rb_cel(1)%ja_eq,mx,my,3_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%ja_eq,rb_cel(1)%ja_eq,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(v).AND.v=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%ve,mx,my,3_i4,nmd,pd)
        CALL lagr_quad_alloc(rb_cel(1)%ve_eq,mx,my,3_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%ve_eq,rb_cel(1)%ve_eq,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%ve,rb_cel(1)%ve,rb(ibl)%id)
        ENDDO
      ELSEIF(PRESENT(ve_eq).AND.ve_eq=='y') THEN
        CALL lagr_quad_alloc(rb_cel(1)%ve_eq,mx,my,3_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%ve_eq,rb_cel(1)%ve_eq,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(e).AND.e=='y')  THEN
        IF (.NOT.transfer_eq) THEN
          CALL lagr_quad_alloc(rb_cel(1)%eef_eq,mx,my,3_i4,pd)
          DO ibl=1,nrbl
            CALL exch_lagr(rb(ibl)%eef_eq,rb_cel(1)%eef_eq,rb(ibl)%id)
          ENDDO
        ENDIF
!
        CALL lagr_quad_alloc(rb_cel(1)%eef,mx,my,3_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%eef,rb_cel(1)%eef,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(p).AND.p=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%pres,mx,my,1_i4,nmd,pd)
        CALL lagr_quad_alloc(rb_cel(1)%pres_eq,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%pres_eq,rb_cel(1)%pres_eq,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%pres,rb_cel(1)%pres,rb(ibl)%id)
        ENDDO
      ELSEIF(PRESENT(p_eq).AND.p_eq=='y') THEN
        CALL lagr_quad_alloc(rb_cel(1)%pres_eq,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%pres_eq,rb_cel(1)%pres_eq,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(ephip).AND.ephip=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%ephip,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%ephip,rb_cel(1)%ephip,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(pe).AND.pe=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%prese,mx,my,1_i4,nmd,pd)
        CALL lagr_quad_alloc(rb_cel(1)%prese_eq,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%prese_eq,rb_cel(1)%prese_eq,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%prese,rb_cel(1)%prese,rb(ibl)%id)
        ENDDO
      ELSEIF(PRESENT(pe_eq).AND.pe_eq=='y') THEN
        CALL lagr_quad_alloc(rb_cel(1)%prese_eq,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%prese_eq,rb_cel(1)%prese_eq,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(n).AND.n=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%nd,mx,my,1_i4,nmd,pd)
        CALL lagr_quad_alloc(rb_cel(1)%nd_eq,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%nd,rb_cel(1)%nd,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%nd_eq,rb_cel(1)%nd_eq,rb(ibl)%id)
        ENDDO
      ELSEIF(PRESENT(nd_eq).AND.nd_eq=='y') THEN
        CALL lagr_quad_alloc(rb_cel(1)%nd_eq,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%nd_eq,rb_cel(1)%nd_eq,rb(ibl)%id)
        ENDDO
      ENDIF


      IF(PRESENT(ne).AND.ne=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%ne,mx,my,1_i4,nmd,pd)
        IF (.NOT. PRESENT(n))                                           &
     &    CALL lagr_quad_alloc(rb_cel(1)%nd_eq,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%nZel,rb_cel(1)%ne,rb(ibl)%id)
        IF (.NOT. PRESENT(n))                                           &
     &    CALL exch_lagr(rb(ibl)%nd_eq,rb_cel(1)%nd_eq,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(ndn).AND.ndn=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%ndn,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%ndn,rb_cel(1)%ndn,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(nimp).AND.nimp=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%nimp,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%nZzt,rb_cel(1)%nimp,rb(ibl)%id)
        ENDDO
        CALL lagr_quad_alloc(rb_cel(1)%nZeff,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%nZeff,rb_cel(1)%nZeff,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(nimpw).AND.nimpw=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%nimpw,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%nZrt,rb_cel(1)%nimpw,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(los).AND.los=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%nz,mx,my,zimp+1,nmd,pd)
        CALL lagr_quad_alloc(rb_cel(1)%qlosd,mx,my,1_i4,nmd,pd)
        CALL lagr_quad_alloc(rb_cel(1)%qlosb,mx,my,1_i4,nmd,pd)
        CALL lagr_quad_alloc(rb_cel(1)%qlosr,mx,my,1_i4,nmd,pd)
        CALL lagr_quad_alloc(rb_cel(1)%qlosl,mx,my,1_i4,nmd,pd)
        CALL lagr_quad_alloc(rb_cel(1)%qlosi,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%nz,rb_cel(1)%nz,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%qlosd,rb_cel(1)%qlosd,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%qlosb,rb_cel(1)%qlosb,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%qlosr,rb_cel(1)%qlosr,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%qlosl,rb_cel(1)%qlosl,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%qlosi,rb_cel(1)%qlosi,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(ti).AND.ti=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%tion,mx,my,1_i4,nmd,pd)
        CALL lagr_quad_alloc(rb_cel(1)%tion_eq,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%tion,rb_cel(1)%tion,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%tion_eq,rb_cel(1)%tion_eq,rb(ibl)%id)
        ENDDO
      ELSEIF(PRESENT(ti_eq).AND.ti_eq=='y') THEN
        CALL lagr_quad_alloc(rb_cel(1)%tion_eq,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%tion_eq,rb_cel(1)%tion_eq,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(te).AND.te=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%tele,mx,my,1_i4,nmd,pd)
        CALL lagr_quad_alloc(rb_cel(1)%tele_eq,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%tele,rb_cel(1)%tele,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%tele_eq,rb_cel(1)%tele_eq,rb(ibl)%id)
        ENDDO
      ELSEIF(PRESENT(te_eq).AND.te_eq=='y') THEN
        CALL lagr_quad_alloc(rb_cel(1)%tele_eq,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%tele_eq,rb_cel(1)%tele_eq,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(conc).AND.conc=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%conc,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%conc,rb_cel(1)%conc,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(qpe) .AND. ALLOCATED(rb(1)%qpe%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%qpe,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%qpe,rb_cel(1)%qpe,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(qpi) .AND. ALLOCATED(rb(1)%qpi%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%qpi,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%qpi,rb_cel(1)%qpi,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(ppe) .AND. ALLOCATED(rb(1)%ppe%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%ppe,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%ppe,rb_cel(1)%ppe,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(ppi) .AND. ALLOCATED(rb(1)%ppi%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%ppi,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%ppi,rb_cel(1)%ppi,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(rpe) .AND. ALLOCATED(rb(1)%rpe%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%rpe,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%rpe,rb_cel(1)%rpe,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(gei) .AND. ALLOCATED(rb(1)%gei%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%gei,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%gei,rb_cel(1)%gei,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(rpi) .AND. ALLOCATED(rb(1)%rpi%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%rpi,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%rpi,rb_cel(1)%rpi,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(phot) .AND. ALLOCATED(rb(1)%phot%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%phot,mx,my,phqty,nmd,ph_pd)
        CALL lagr_quad_alloc(rb_cel(1)%phot_eq,mx,my,phqty,ph_pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%phot,rb_cel(1)%phot,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(Fele) .AND. ALLOCATED(rb(1)%Fele%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%Fele,mx,my,nF(1),nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%Fele,rb_cel(1)%Fele,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(Fion) .AND. ALLOCATED(rb(1)%Fion%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%Fion,mx,my,nF(2),nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%Fion,rb_cel(1)%Fion,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(Fhot) .AND. ALLOCATED(rb(1)%Fhot%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%Fhot,mx,my,nF(3),nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%Fhot,rb_cel(1)%Fhot,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(Fele_eq) .AND. ALLOCATED(rb(1)%Fele_eq%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%Fele_eq,mx,my,nF(1),pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%Fele_eq,rb_cel(1)%Fele_eq,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(Fion_eq) .AND. ALLOCATED(rb(1)%Fion_eq%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%Fion_eq,mx,my,nF(2),pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%Fion_eq,rb_cel(1)%Fion_eq,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(Fhot_eq) .AND. ALLOCATED(rb(1)%Fhot_eq%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%Fhot_eq,mx,my,nF(3),pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%Fhot_eq,rb_cel(1)%Fhot_eq,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(bmod) .AND. ALLOCATED(rb(1)%bmod%fs)) THEN
        CALL lagr_quad_alloc(rb_cel(1)%bmod,mx,my,1_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%bmod,rb_cel(1)%bmod,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(ramos2).AND.ramos2=='y') THEN
        CALL lagr_quad_alloc(rb_cel(1)%ramos2,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%ramos2,rb_cel(1)%ramos2,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(flx).AND.flx=='y') THEN
        CALL lagr_quad_alloc(rb_cel(1)%exb,mx,my,3_i4,pd)
        CALL lagr_quad_alloc(rb_cel(1)%piv,mx,my,3_i4,pd)
        CALL lagr_quad_alloc(rb_cel(1)%keflx,mx,my,3_i4,pd)
        CALL lagr_quad_alloc(rb_cel(1)%qi0,mx,my,3_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%exb,rb_cel(1)%exb,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%piv,rb_cel(1)%piv,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%keflx,rb_cel(1)%keflx,rb(ibl)%id)
          CALL exch_lagr(rb(ibl)%qi0,rb_cel(1)%qi0,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(poten).AND.poten=='y') THEN
        CALL lagr_quad_alloc(rb_cel(1)%poten,mx,my,3_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%poten,rb_cel(1)%poten,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(partvar_eq).AND.partvar_eq=='y')  THEN
!        CALL lagr_quad_alloc(rb_cel(1)%partvar_eq,mx,my,8,phf_pd)
        CALL lagr_quad_alloc(rb_cel(1)%partvar_eq,mx,my,10,phf_pd)
        CALL lagr_quad_alloc(rb_cel(1)%partvar,mx,my,6,nmd,phf_pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%partvar_eq,rb_cel(1)%partvar_eq,       &
     &                   rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(rwork2).AND.rwork2=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%rwork2,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%rwork2,rb_cel(1)%rwork2,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(f).AND.f=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%fp,mx,my,3_i4,nmd,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%fp,rb_cel(1)%fp,rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(diff_shape).AND.diff_shape=='y')  THEN
        CALL lagr_quad_alloc(rb_cel(1)%diff_shape,mx,my,ds_nqty,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%diff_shape,rb_cel(1)%diff_shape,       &
    &                    rb(ibl)%id)
        ENDDO
      ENDIF

      IF(PRESENT(fsa_beq2).AND.fsa_beq2=='y'.AND.dump_fsa_beq2)  THEN
        CALL lagr_quad_alloc(rb_cel(1)%fsa_beq2,mx,my,1_i4,pd)
        DO ibl=1,nrbl
          CALL exch_lagr(rb(ibl)%fsa_beq2,rb_cel(1)%fsa_beq2,           &
    &                    rb(ibl)%id)
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE map_type_init
!-----------------------------------------------------------------------
!     subprogram 6. rz_to_xy_hard.
!     Mapping Problem: given (r,z) find (xc,yc) in grid coordinates
!     Moves from starting point to desired point along vector direction
!      working in grid coordinates.
!     xc,yc should come in with a good guess
!     Hard boundary:
!     In this routine, if the logical coordinate goes immediately beyond
!      the boundary, then it exits with a fail (i.e., it has a hard
!      boundary.  This is useful for things like field lines where one
!      knows that there is a possibility of exiting the domain.
!-----------------------------------------------------------------------
      SUBROUTINE rz_to_xy_hard(r,z,xc,yc,fail,err)
      USE local
      USE input
      USE physdat
      IMPLICIT NONE
      REAL(r8), INTENT(IN) :: r,z
      REAL(r8), INTENT(INOUT) :: xc,yc,err
      LOGICAL, INTENT(INOUT) :: fail

      INTEGER(i4), PARAMETER :: max_loop = 30000
      INTEGER(i4) :: iter_loop
      REAL(r8), PARAMETER :: etol = 1.e-10
      REAL(r8) :: dr,dz,dx,dy,jacobian,md,zt,delz
      REAL(r8) :: drdx,dzdx,drdy,dzdy,yt,rc,zc
      REAL(r8) :: slow_down,rzf(2),rzfx(2),rzfy(2)
      LOGICAL, PARAMETER :: diagnose=.FALSE.
      slow_down=0.2
!-----------------------------------------------------------------------
!     diagnose
!-----------------------------------------------------------------------
      IF (diagnose) THEN
        WRITE(*,*) 'Rdesire, Zdesire ', r,z
        OPEN(UNIT=92,FILE='path.dat',STATUS='unknown',POSITION='REWIND')
      ENDIF
!-----------------------------------------------------------------------
!     check that rz are reasonable.
!-----------------------------------------------------------------------
      IF (gridshape/='rect'.AND.periodicity/='y-dir'.OR.                &
     &    periodicity/='both') THEN
        IF (r>rmax.OR.r<rmin.OR.z>zmax.OR.z<zmin) THEN
          fail=.TRUE.
          RETURN
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     main loop
!-----------------------------------------------------------------------
      iter_loop = 0_i4
      DO
!-----------------------------------------------------------------------
!       Evaluate new r and z - don't allow y to fly off the grid
!       assume periodicity
!-----------------------------------------------------------------------
        yt = yc
        md = FLOOR(ABS(yt/my))
        IF(yt<0) yt=(md+1)*my+yt
        IF(yt>my) yt=-md*my+yt

        CALL lagr_quad_2D_eval_ts(rb_cel(1)%rz,xc,yt,1_i4,rzf,rzfx,rzfy)
        IF (diagnose) WRITE(92,*) rzf(1),rzf(2)
!-----------------------------------------------------------------------
!       test for convergence. dr and dz are difference
!-----------------------------------------------------------------------
        zc = z ; rc = r
        IF(gridshape == 'rect')THEN
         SELECT CASE(periodicity)
         CASE("y-dir")
           IF(zc>ymax) zc=ymin+MOD(zc-ymin,ymax-ymin)
           IF(zc<ymin) zc=ymax+MOD(zc-ymin,ymax-ymin)
         CASE("both")
           IF(zc>ymax) zc=ymin+MOD(zc-ymin,ymax-ymin)
           IF(zc<ymin) zc=ymax+MOD(zc-ymin,ymax-ymin)
           IF(rc>xmax) rc=xmin+MOD(rc-xmin,xmax-xmin)
           IF(rc<xmin) rc=xmax+MOD(rc-xmin,xmax-xmin)
         END SELECT
        ENDIF
        dr=rc-rzf(1)
        dz=zc-rzf(2)
        err=abs(dr)+abs(dz)
        IF(err<etol) THEN
         fail = .FALSE.
         EXIT
        ENDIF
!-----------------------------------------------------------------------
!       limit the number of iterations.
!-----------------------------------------------------------------------
        iter_loop=iter_loop+1
        IF(iter_loop>max_loop.or.xc>mx) THEN
          fail = .true.
          IF (diagnose) CLOSE(92)
          RETURN
        ENDIF
!-----------------------------------------------------------------------
!       Update x and y
!-----------------------------------------------------------------------
        drdx = rzfx(1)
        drdy = rzfy(1)
        dzdx = rzfx(2)
        dzdy = rzfy(2)
        ! Jac MUST be positive
        jacobian=ABS(1./(drdx*dzdy-drdy*dzdx))

        dx = (dzdy*dr - drdy*dz)*jacobian
        dy = (drdx*dz - dzdx*dr)*jacobian
!-----------------------------------------------------------------------
!       Slow down near the axis
!-----------------------------------------------------------------------
        IF(xc < 0.05*REAL(mx)) slow_down=0.02
        ! Jac too small
        IF(xc+dx*slow_down <=0) dx=-0.1*xc
        xc=xc+dx*slow_down
        yc=yc+dy*slow_down

      ENDDO
      md = FLOOR(ABS(yc/my))
      IF(yc<0.) yc=(md+1)*my+yc
      IF(yc>my) yc=-md*my+yc

      IF(diagnose) CLOSE(92)
      RETURN
      END SUBROUTINE rz_to_xy_hard 
!-----------------------------------------------------------------------
!     End module
!-----------------------------------------------------------------------
      END MODULE map_mod