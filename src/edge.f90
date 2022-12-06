!-----------------------------------------------------------------------
!     &Id: edge.f90 4620 2015-12-03 19:18:40Z jking $
!     module for handling communication across grid block borders.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.
!-----------------------------------------------------------------------
!     1.  edge_init.
!     2.  edge_segment_init.
!     3.  edge_network.
!     4.  edge_accumulate.
!     5.  edge_comp_accumulate.
!     6.  edge_seg_network.
!     7.  edge_seg_accumulate.
!     8.  edge_seg_comp_accumulate.
!     9.  edge_load_arr.
!     10. edge_load_carr.
!     11. edge_load_2D_carr.
!     12. edge_unload_arr.
!     13. edge_unload_carr.
!     14. edge_unload_2D_carr.
!     15. edge_load_save.
!     16. edge_load_csave.
!     17. edge_add_save.
!     18. edge_add_csave.
!     19. edge_zero_save.
!     20. edge_zero_csave.
!-----------------------------------------------------------------------
!     declarations.
!-----------------------------------------------------------------------
      MODULE edge
      USE local
      USE edge_type_mod
      USE seam_storage_mod
      USE vector_type_mod
      IMPLICIT NONE

      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 1. edge_init.
!     initializes arrays used for border communication.
!-----------------------------------------------------------------------
      SUBROUTINE edge_init(nq_alloc,nq_save)
      USE fields

      INTEGER(i4), INTENT(IN) :: nq_alloc,nq_save

      INTEGER(i4) :: ibl,iv,nv,np,npb,ipb,ip,n0,in,jbv,jvv
      INTEGER(i4), DIMENSION(1) :: loc
      LOGICAL :: match
      LOGICAL, DIMENSION(:), ALLOCATABLE :: ord_mask
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ord_arr
!-----------------------------------------------------------------------
!     create a second ptr array without references to seam0 and with
!     interior vertex label cross references for parallel coding.
!-----------------------------------------------------------------------
      max_imags=0
      DO ibl=1,nbl
        DO iv=1,seam(ibl)%nvert
          np=SIZE(seam(ibl)%vertex(iv)%ptr,2)
          n0=0
          DO ip=1,np
            IF (seam(ibl)%vertex(iv)%ptr(1,ip)==0) n0=n0+1
          ENDDO
          npb=np-n0
          seam(ibl)%vertex(iv)%nimage = npb
          max_imags=MAX(npb,max_imags)
          ALLOCATE(seam(ibl)%vertex(iv)%ptr2(2,npb))
          ipb=0
          DO ip=1,np
            IF (seam(ibl)%vertex(iv)%ptr(1,ip) == 0) CYCLE
            ipb=ipb+1
            seam(ibl)%vertex(iv)%ptr2(1:2,ipb) =                        &
     &           seam(ibl)%vertex(iv)%ptr(1:2,ip)
          ENDDO
          IF (ipb/=npb) CALL nim_stop('Error in edge_init.')
          np=seam(ibl)%vertex(iv)%nimage+1
          ALLOCATE(seam(ibl)%vertex(iv)%seam_in(nq_alloc))
          ALLOCATE(seam(ibl)%vertex(iv)%seam_cin(nq_alloc))
          ALLOCATE(seam(ibl)%vertex(iv)%seam_out(nq_alloc))
          ALLOCATE(seam(ibl)%vertex(iv)%seam_cout(nq_alloc))
          ALLOCATE(seam(ibl)%vertex(iv)%seam_hold(np,nq_alloc))
          ALLOCATE(seam(ibl)%vertex(iv)%seam_chold(np,nq_alloc))
          ALLOCATE(seam(ibl)%vertex(iv)%seam_save(nq_save))
          ALLOCATE(seam(ibl)%vertex(iv)%seam_csave(nq_save))
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     determine a unique order for summing seam vertices with more than
!     two images to prevent the occurrence of different round-off errors
!     on different blocks, which leads to an instability.
!
!     degenerate points have one unshared connection per block, which
!     requires special treatment.  ordering is according to sending
!     vertex, so the order of each contribution on the receiving side
!     must be the same as the order on the sending side, even when they
!     have different connections in this unique case of degens.
!-----------------------------------------------------------------------
      ALLOCATE(ord_mask(max_imags+1),ord_arr(max_imags+1))
      DO ibl=1,nbl
        nv=seam(ibl)%nvert
        DO iv=1,nv
          np=seam(ibl)%vertex(iv)%nimage
          ALLOCATE(seam(ibl)%vertex(iv)%order(0:np))
          ord_arr(1)=seam(ibl)%id+1.e-9*iv
          IF (ibl<=nrbl) THEN
            IF (rb(ibl)%degenerate.AND.iv==2*rb(ibl)%mx+rb(ibl)%my)     &
     &        ord_arr(1)=-ord_arr(1)
          ENDIF
          DO ip=1,np
            ord_arr(ip+1)=seam(ibl)%vertex(iv)%ptr2(1,ip)               &
     &             +1.e-9*seam(ibl)%vertex(iv)%ptr2(2,ip)
            IF (ibl<=nrbl) THEN
              IF (rb(ibl)%degenerate.AND.iv==nv) THEN
                jbv=seam(ibl)%vertex(iv)%ptr2(1,ip)
                jvv=seam(ibl)%vertex(iv)%ptr2(2,ip)
                IF (jvv==seam(ibl)%vertex(1)%ptr2(2,1)+1.AND.           &
     &              jbv==seam(ibl)%vertex(1)%ptr2(1,1))                 &
     &            ord_arr(ip+1)=-ord_arr(ip+1)
              ENDIF
            ENDIF
          ENDDO
          ord_mask(1:np+1)=.true.
          DO ip=1,np+1
            loc=MINLOC(ord_arr(1:np+1),ord_mask(1:np+1))
            in=loc(1)
            seam(ibl)%vertex(iv)%order(in-1)=ip
            ord_mask(in)=.false.
          ENDDO
          ord_mask(1:np+1)=.true.
          DO ip=0,np
            ord_mask(seam(ibl)%vertex(iv)%order(ip))=.false.
          ENDDO
          DO ip=1,np+1
            IF (ord_mask(ip))                                           &
     &        CALL nim_stop('Edge_init: vertex order not unique.')
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(ord_mask,ord_arr)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_init
!-----------------------------------------------------------------------
!     subprogram 2. edge_segment_init.
!     initializes the segment part of seam structures used for border
!     communication of matrix elements.  note that this routine requires
!     having expoint and r0point defined within the seam structure, and
!     those flags are set in other initialization routines.
!-----------------------------------------------------------------------
      SUBROUTINE edge_segment_init(nq,nqm,offmat)
      USE fields
      USE pardata

      INTEGER(i4), INTENT(IN) :: nq,nqm,offmat

      INTEGER(i4) :: iv,ip,np,ibl,ibv,npb,ipb,nv,ivp,icell,             &
     &               jbv,jvv,jbpr,jvpr,ippr,ivt,ixp,iyp,ixn,iyn
      INTEGER(i4), DIMENSION(3) :: vdiff1,vdiff2
      LOGICAL :: match
      CHARACTER(64) :: msg
!-----------------------------------------------------------------------
!     create internal references and communication pointers for
!     edge-segment centered quantities.  the latter uses the fact that 
!     seams running around two adjacent blocks always advance in 
!     opposite directions.
!-----------------------------------------------------------------------
      DO ibl=1,nbl
        nv=seam(ibl)%nvert
        ALLOCATE(seam(ibl)%segment(nv))
        DO iv=1,nv
          ivp=iv-1
          IF (ivp==0) ivp=nv
          seam(ibl)%segment(iv)%intxyn=seam(ibl)%vertex(iv)%intxy
          seam(ibl)%segment(iv)%intxyp=seam(ibl)%vertex(ivp)%intxy
          seam(ibl)%segment(iv)%ptr=0
          IF (ibl>nrbl) THEN
            match=.false.
            ivt=seam(ibl)%vertex(iv)%intxy(1)
            DO ip=1,SIZE(tb(ibl)%tgeom%neighbor(ivt)%vertex)-1
              np=tb(ibl)%tgeom%neighbor(ivt)%vertex(ip)
              IF (np==seam(ibl)%segment(iv)%intxyp(1)) THEN
                seam(ibl)%segment(iv)%intxyn(2)=ip
                match=.true.
                CYCLE
              ENDIF
            ENDDO
            IF (.NOT.match) THEN
              WRITE(msg,'(a,i3,a,i3)')                                  &
     &          'No tri 1 match in edge_segment_init,ibl=',ibl,'iv=',iv
              CALL nim_stop(msg)
            ENDIF
            match=.false.
            ivt=seam(ibl)%vertex(ivp)%intxy(1)
            DO ip=1,SIZE(tb(ibl)%tgeom%neighbor(ivt)%vertex)-1
              np=tb(ibl)%tgeom%neighbor(ivt)%vertex(ip)
              IF (np==seam(ibl)%segment(iv)%intxyn(1)) THEN
                seam(ibl)%segment(iv)%intxyp(2)=ip
                match=.true.
                CYCLE
              ENDIF
            ENDDO
            IF (.NOT.match) THEN
              WRITE(msg,'(a,i3,a,i3)')                                  &
     &          'No tri 2 match in edge_segment_init,ibl=',ibl,'iv=',iv
              CALL nim_stop(msg)
            ENDIF
          ELSE
            IF (rb(ibl)%degenerate.AND.iv>2*rb(ibl)%mx+rb(ibl)%my)      &
     &        CYCLE
          ENDIF
          IF ((seam(ibl)%expoint(iv).OR.seam(ibl)%r0point(iv)).AND.     &
     &        (seam(ibl)%expoint(ivp).OR.                               &
     &         seam(ibl)%r0point(ivp))) CYCLE
          match=.false.
          ptr_loop: DO ip=1,SIZE(seam(ibl)%vertex(iv)%ptr2,2)
            jbv=seam(ibl)%vertex(iv)%ptr2(1,ip)
            jvv=seam(ibl)%vertex(iv)%ptr2(2,ip)
            IF (jvv==block_sizes(1,jbv)) jvv=0
            DO ippr=1,SIZE(seam(ibl)%vertex(ivp)%ptr2,2)
              jbpr=seam(ibl)%vertex(ivp)%ptr2(1,ippr)
              jvpr=seam(ibl)%vertex(ivp)%ptr2(2,ippr)
              IF (jbpr==jbv.AND.jvv==jvpr-1) THEN
                match=.true.
                seam(ibl)%segment(iv)%ptr(1)=jbv
                seam(ibl)%segment(iv)%ptr(2)=jvpr
                EXIT ptr_loop
              ENDIF
            ENDDO
          ENDDO ptr_loop
          IF (.NOT.match) THEN
            WRITE(msg,'(a,i3,a,i3)')                                    &
     &        'No edge match in edge_segment_init, ibl=', ibl,' iv=',iv
            CALL nim_stop(msg)
          ENDIF
        ENDDO
!-----------------------------------------------------------------------
!       allocate arrays:  those used for communicating off-diagonal
!       elements of matrices should have two for vertex-vertex elements,
!       and possibly two for each off-diagonal vertex-side or side-side
!       pair, depending on whether the preconditioner needs this
!       communication.
!-----------------------------------------------------------------------
        DO iv=1,seam(ibl)%nvert
          ALLOCATE(seam(ibl)%segment(iv)%seam_in(nq))
          ALLOCATE(seam(ibl)%segment(iv)%seam_cin(nq))
          ALLOCATE(seam(ibl)%segment(iv)%seam_out(nq))
          ALLOCATE(seam(ibl)%segment(iv)%seam_cout(nq))
          ALLOCATE(seam(ibl)%segment(iv)%seam_save(nq))
          ALLOCATE(seam(ibl)%segment(iv)%seam_csave(nq))
          ALLOCATE(seam(ibl)%segment(iv)%seam_mat_in(nqm,nqm,offmat))
          ALLOCATE(seam(ibl)%segment(iv)%seam_mat_cin(nqm,nqm,offmat))
          ALLOCATE(seam(ibl)%segment(iv)%seam_mat_out(nqm,nqm,offmat))
          ALLOCATE(seam(ibl)%segment(iv)%seam_mat_cout(nqm,nqm,offmat))
          ALLOCATE(seam(ibl)%segment(iv)%seam_mat_save(nqm,nqm,offmat))
          ALLOCATE(seam(ibl)%segment(iv)%seam_mat_csave(nqm,nqm,offmat))
        ENDDO
!-----------------------------------------------------------------------
!       indices for side-centered arrays.
!-----------------------------------------------------------------------
        IF (ibl>nrbl) THEN
          tri_side: DO iv=1,seam(ibl)%nvert
            seam(ibl)%segment(iv)%h_side=.true.
            DO icell=1,tb(ibl)%tgeom%mcell
              vdiff1=tb(ibl)%tgeom%vertex(icell,:)-                     &
     &          seam(ibl)%segment(iv)%intxyn(1)
              vdiff2=tb(ibl)%tgeom%vertex(icell,:)-                     &
     &          seam(ibl)%segment(iv)%intxyp(1)
              IF (MINVAL(ABS(vdiff1))==0.AND.MINVAL(ABS(vdiff2))==0)    &
     &          THEN
                seam(ibl)%segment(iv)%intxys=(/icell,1_i4/)
                CYCLE tri_side
              ENDIF
            ENDDO
          ENDDO tri_side
        ELSE
          DO iv=1,rb(ibl)%mx
            seam(ibl)%segment(iv)%h_side=.true.
            seam(ibl)%segment(iv)%intxys=seam(ibl)%segment(iv)%intxyn
          ENDDO
          DO iv=rb(ibl)%mx+1,rb(ibl)%mx+rb(ibl)%my
            seam(ibl)%segment(iv)%h_side=.false.
            seam(ibl)%segment(iv)%intxys=seam(ibl)%segment(iv)%intxyn
          ENDDO
          DO iv=rb(ibl)%mx+rb(ibl)%my+1,2*rb(ibl)%mx+rb(ibl)%my
            seam(ibl)%segment(iv)%h_side=.true.
            seam(ibl)%segment(iv)%intxys=seam(ibl)%segment(iv)%intxyp
          ENDDO
          DO iv=2*rb(ibl)%mx+rb(ibl)%my+1,seam(ibl)%nvert
            seam(ibl)%segment(iv)%h_side=.false.
            seam(ibl)%segment(iv)%intxys=seam(ibl)%segment(iv)%intxyp
          ENDDO
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_segment_init
!-----------------------------------------------------------------------
!     subprogram 3. edge_network.
!     manages the communication between blocks.  note that the first
!     n_four Fourier components of the complex arrays are communicated
!     if n_four>0.  IF n_four=0, the real arrays are communicated.
!-----------------------------------------------------------------------
      SUBROUTINE edge_network(nqty,n_four,n_side,tx_rhs)
      USE fields
      USE mpi_nim
      USE pardata
      USE time
      USE computation_pointers

      INTEGER(i4), INTENT(IN) :: nqty,n_four,n_side
      LOGICAL, INTENT(IN) :: tx_rhs

      INTEGER(i4) :: ierror,ib,iv,nv,kv,nfnq,is,iq,im
      REAL(r8), DIMENSION(2*nqty,nrbl) :: dg_save
      COMPLEX(r8), DIMENSION(2*nqty*n_four,nrbl) :: dg_csave
!-----------------------------------------------------------------------
!     if tx_rhs is true, communicate the first nqty items in rb%rhs,
!     otherwise seam%vertex%seam_in should be loaded already.
!-----------------------------------------------------------------------
      nfnq=n_four*nqty
      IF (tx_rhs) THEN
        IF (n_four>0) THEN
          DO ib=1,nbl
            CALL edge_load_carr(crhs(ib),nqty,1_i4,n_four,n_side,       &
     &                          seam(ib))
          ENDDO
        ELSE
          DO ib=1,nbl
            CALL edge_load_arr(rhs(ib),nqty,n_side,seam(ib))
          ENDDO
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     communicate element side-centered data, if needed.
!-----------------------------------------------------------------------
      IF (n_side/=0) THEN
        IF (n_four>0) THEN
          CALL edge_seg_network(nfnq*n_side,0_i4,.true.)
        ELSE
          CALL edge_seg_network(nqty*n_side,0_i4,.false.)
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     start timer for vertex-centered data only.
!-----------------------------------------------------------------------
      CALL timer(timestart)
!-----------------------------------------------------------------------
!     for degenerate rblocks, save the seam_in information at the
!     corners of the degenerate boundary, then do an internal
!     accumulation there.
!-----------------------------------------------------------------------
      IF (n_four>0) THEN
        DO ib=1,nrbl
          IF (rb(ib)%degenerate) THEN
            kv=2*rb(ib)%mx+rb(ib)%my
            nv=seam(ib)%nvert
            dg_csave(1:nfnq,ib)=seam(ib)%vertex(kv)%seam_cin(1:nfnq)
            dg_csave(nfnq+1:2*nfnq,ib)=                                 &
     &        seam(ib)%vertex(nv)%seam_cin(1:nfnq)
            DO iv=kv,nv-1
              seam(ib)%vertex(nv)%seam_cin(1:nfnq)=                     &
     &                seam(ib)%vertex(nv)%seam_cin(1:nfnq)              &
     &               +seam(ib)%vertex(iv)%seam_cin(1:nfnq)
            ENDDO
            seam(ib)%vertex(kv)%seam_cin(1:nfnq)=0
            IF (n_side>0) THEN
              DO iv=kv+1,nv
                DO im=0,(n_four-1)*nqty,nqty
                  iq=1
                  DO is=1,n_side
                    seam(ib)%vertex(nv)%seam_cin(im+1:im+nqty)=         &
     &                 seam(ib)%vertex(nv)%seam_cin(im+1:im+nqty)       &
     &                +seam(ib)%segment(iv)%seam_cin                    &
     &                   (im*n_side+iq:im*n_side+nqty+iq-1)
                    iq=iq+nqty
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ELSE
        DO ib=1,nrbl
          IF (rb(ib)%degenerate) THEN
            kv=2*rb(ib)%mx+rb(ib)%my
            nv=seam(ib)%nvert
            dg_save(1:nqty,ib)=seam(ib)%vertex(kv)%seam_in(1:nqty)
            dg_save(nqty+1:2*nqty,ib)=                                  &
     &        seam(ib)%vertex(nv)%seam_in(1:nqty)
            DO iv=kv,nv-1
              seam(ib)%vertex(nv)%seam_in(1:nqty)=                      &
     &                seam(ib)%vertex(nv)%seam_in(1:nqty)               &
     &               +seam(ib)%vertex(iv)%seam_in(1:nqty)
            ENDDO
            seam(ib)%vertex(kv)%seam_in(1:nqty)=0
            IF (n_side>0) THEN
              DO iv=kv+1,nv
                iq=1
                DO is=1,n_side
                  seam(ib)%vertex(nv)%seam_in(1:nqty)=                  &
     &                    seam(ib)%vertex(nv)%seam_in(1:nqty)           &
     &                   +seam(ib)%segment(iv)%seam_in(iq:nqty+iq-1)
                  iq=iq+nqty
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     communicate vertex-centered data, including an entire block's
!     contribution to a degenerate point.
!-----------------------------------------------------------------------
      IF (nprocs == 1) THEN
        IF (n_four>0) THEN
          CALL edge_comp_accumulate(nfnq)
        ELSE
          CALL edge_accumulate(nqty)
        ENDIF
      ELSE
        IF (n_four>0) THEN
          CALL parallel_seam_comm_comp(nfnq)
        ELSE
          CALL parallel_seam_comm(nqty)
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     for degenerate rblocks, do an internal scatter and restore
!     seam_in at the corners of the degenerate boundary.
!-----------------------------------------------------------------------
      IF (n_four>0) THEN
        DO ib=1,nrbl
          IF (rb(ib)%degenerate) THEN
            kv=2*rb(ib)%mx+rb(ib)%my
            nv=seam(ib)%nvert
            DO iv=kv,nv-1
              seam(ib)%vertex(iv)%seam_cout(1:nfnq)=                    &
     &                seam(ib)%vertex(nv)%seam_cout(1:nfnq)
            ENDDO
            seam(ib)%vertex(kv)%seam_cin(1:nfnq)=dg_csave(1:nfnq,ib)
            seam(ib)%vertex(nv)%seam_cin(1:nfnq)=                       &
     &        dg_csave(nfnq+1:2*nfnq,ib)
            IF (n_side>0) THEN
              DO iv=kv+1,nv
                DO im=0,(n_four-1)*nqty,nqty
                  iq=1
                  DO is=1,n_side
                    seam(ib)%segment(iv)%seam_cout                      &
     &                (im*n_side+iq:im*n_side+nqty+iq-1)=               &
     &                seam(ib)%vertex(nv)%seam_cout(im+1:im+nqty)
                    iq=iq+nqty
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ELSE
        DO ib=1,nrbl
          IF (rb(ib)%degenerate) THEN
            kv=2*rb(ib)%mx+rb(ib)%my
            nv=seam(ib)%nvert
            DO iv=kv,nv-1
              seam(ib)%vertex(iv)%seam_out(1:nqty)=                     &
     &                seam(ib)%vertex(nv)%seam_out(1:nqty)
            ENDDO
            seam(ib)%vertex(kv)%seam_in(1:nqty)=dg_save(1:nqty,ib)
            seam(ib)%vertex(nv)%seam_in(1:nqty)=                        &
     &        dg_save(nqty+1:2*nqty,ib)
            IF (n_side>0) THEN
              DO iv=kv+1,nv
                iq=1
                DO is=1,n_side
                  seam(ib)%segment(iv)%seam_out(iq:nqty+iq-1)=          &
     &                    seam(ib)%vertex(nv)%seam_out(1:nqty)
                  iq=iq+nqty
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     complete communication timing.
!-----------------------------------------------------------------------
      CALL timer(timeend)
      time_seam = time_seam + timeend-timestart
      seamcount = seamcount + 1
!-----------------------------------------------------------------------
!     copy to the rhs arrays if needed.
!-----------------------------------------------------------------------
      IF (tx_rhs) THEN
        IF (n_four>0) THEN
          DO ib=1,nbl
            CALL edge_unload_carr(crhs(ib),nqty,1_i4,n_four,n_side,     &
     &                            seam(ib))
          ENDDO
        ELSE
          DO ib=1,nbl
            CALL edge_unload_arr(rhs(ib),nqty,n_side,seam(ib))
          ENDDO
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_network
!-----------------------------------------------------------------------
!     subprogram 4. edge_accumulate.
!     accumulate all representations of edge vertex information from
!     seam_in into seam_out.
!-----------------------------------------------------------------------
      SUBROUTINE edge_accumulate(nqty)

      INTEGER(i4), INTENT(IN) :: nqty

      INTEGER(i4) :: iv,ib,jv,jb,ip,ni,il
!-----------------------------------------------------------------------
!     accumulate the different representations.  place them in an
!     array in a specified order (same on all sides of the border)
!     before summing to ensure the same round-off error.
!-----------------------------------------------------------------------
      DO ib=1,SIZE(seam)
        DO iv=1,seam(ib)%nvert
          ni=seam(ib)%vertex(iv)%nimage
          SELECT CASE(ni)
          CASE(0)
            seam(ib)%vertex(iv)%seam_out(1:nqty)=                       &
     &         seam(ib)%vertex(iv)%seam_in(1:nqty)
          CASE(1)
            jb=seam(ib)%vertex(iv)%ptr2(1,1)
            jv=seam(ib)%vertex(iv)%ptr2(2,1)
            seam(ib)%vertex(iv)%seam_out(1:nqty)=                       &
     &         seam(ib)%vertex(iv)%seam_in(1:nqty)                      &
     &        +seam(jb)%vertex(jv)%seam_in(1:nqty)
          CASE DEFAULT
            il=seam(ib)%vertex(iv)%order(0)
            seam(ib)%vertex(iv)%seam_hold(il,1:nqty)=                   &
     &        seam(ib)%vertex(iv)%seam_in(1:nqty)
            DO ip=1,ni
              il=seam(ib)%vertex(iv)%order(ip)
              jb=seam(ib)%vertex(iv)%ptr2(1,ip)
              jv=seam(ib)%vertex(iv)%ptr2(2,ip)
              seam(ib)%vertex(iv)%seam_hold(il,1:nqty)=                 &
     &          seam(jb)%vertex(jv)%seam_in(1:nqty)
            ENDDO
            seam(ib)%vertex(iv)%seam_out(1:nqty)=                       &
     &        SUM(seam(ib)%vertex(iv)%seam_hold(:,1:nqty),1)
          END SELECT
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_accumulate
!-----------------------------------------------------------------------
!     subprogram 5. edge_comp_accumulate.
!     accumulate all representations of edge vertex information from
!     seam_cin into seam_cout.  complex version of edge_accumulate
!-----------------------------------------------------------------------
      SUBROUTINE edge_comp_accumulate(nfnq)

      INTEGER(i4), INTENT(IN) :: nfnq

      INTEGER(i4) :: iv,ib,jv,jb,ip,ni,il
!-----------------------------------------------------------------------
!     accumulate the different representations.  place them in an
!     array in a specified order (same on all sides of the border)
!     before summing to ensure the same round-off error.
!-----------------------------------------------------------------------
      DO ib=1,SIZE(seam)
        DO iv=1,seam(ib)%nvert
          ni=seam(ib)%vertex(iv)%nimage
          SELECT CASE(ni)
          CASE(0)
            seam(ib)%vertex(iv)%seam_cout(1:nfnq)=                      &
     &         seam(ib)%vertex(iv)%seam_cin(1:nfnq)
          CASE(1)
            jb=seam(ib)%vertex(iv)%ptr2(1,1)
            jv=seam(ib)%vertex(iv)%ptr2(2,1)
            seam(ib)%vertex(iv)%seam_cout(1:nfnq)=                      &
     &         seam(ib)%vertex(iv)%seam_cin(1:nfnq)                     &
     &        +seam(jb)%vertex(jv)%seam_cin(1:nfnq)
          CASE DEFAULT
            il=seam(ib)%vertex(iv)%order(0)
            seam(ib)%vertex(iv)%seam_chold(il,1:nfnq)=                  &
     &        seam(ib)%vertex(iv)%seam_cin(1:nfnq)
            DO ip=1,ni
              il=seam(ib)%vertex(iv)%order(ip)
              jb=seam(ib)%vertex(iv)%ptr2(1,ip)
              jv=seam(ib)%vertex(iv)%ptr2(2,ip)
              seam(ib)%vertex(iv)%seam_chold(il,1:nfnq)=                &
     &          seam(jb)%vertex(jv)%seam_cin(1:nfnq)
            ENDDO
            seam(ib)%vertex(iv)%seam_cout(1:nfnq)=                      &
     &        SUM(seam(ib)%vertex(iv)%seam_chold(:,1:nfnq),1)
          END SELECT
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_comp_accumulate
!-----------------------------------------------------------------------
!     subprogram 6. edge_seg_network.
!     manages the communication between blocks for segment-centered
!     information.  this may be matrix elements or side-centered
!     data.
!-----------------------------------------------------------------------
      SUBROUTINE edge_seg_network(nqty,nmat,comp)
      USE mpi_nim
      USE pardata
      USE time

      INTEGER(i4), INTENT(IN) :: nqty,nmat
      LOGICAL, INTENT(IN) :: comp

      INTEGER(i4) :: ierror,ib,iv,im
!-----------------------------------------------------------------------
!     if nmat>0, communicate nqtyXnqty matrix blocks, where opposite
!     direction connections are stored consecutively--symmetry is no
!     longer assumed.  if not, communicate a 1D array of data.
!-----------------------------------------------------------------------
      CALL timer(timestart)
      IF (MODULO(nmat,2_i4)/=0) CALL nim_stop                           &
     &  ('Edge_seg_network:  nmat must be even.')
!-----------------------------------------------------------------------
!     complex data.
!-----------------------------------------------------------------------
      IF (comp) THEN
        IF (nmat>0) THEN
          DO ib=1,SIZE(seam)
            DO iv=1,seam(ib)%nvert
              DO im=1,nmat
                seam(ib)%segment(iv)%seam_mat_cout(1:nqty,1:nqty,im)=   &
     &            seam(ib)%segment(iv)%seam_mat_cin(1:nqty,1:nqty,im)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO ib=1,SIZE(seam)
            DO iv=1,seam(ib)%nvert
              seam(ib)%segment(iv)%seam_cout(1:nqty)=                   &
     &          seam(ib)%segment(iv)%seam_cin(1:nqty)
            ENDDO
          ENDDO
        ENDIF
        IF (nprocs == 1) THEN
          CALL edge_seg_comp_accumulate(nqty,nmat)
        ELSE IF (nmat>0) THEN
          CALL parallel_matseg_comm_comp(nqty,nmat)
        ELSE
          CALL parallel_seg_comm_comp(nqty)
        ENDIF
!-----------------------------------------------------------------------
!     real data.
!-----------------------------------------------------------------------
      ELSE
        IF (nmat>0) THEN
          DO ib=1,SIZE(seam)
            DO iv=1,seam(ib)%nvert
              DO im=1,nmat
                seam(ib)%segment(iv)%seam_mat_out(1:nqty,1:nqty,im)=    &
     &            seam(ib)%segment(iv)%seam_mat_in(1:nqty,1:nqty,im)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO ib=1,SIZE(seam)
            DO iv=1,seam(ib)%nvert
              seam(ib)%segment(iv)%seam_out(1:nqty)=                    &
     &          seam(ib)%segment(iv)%seam_in(1:nqty)
            ENDDO
          ENDDO
        ENDIF
        IF (nprocs == 1) THEN
          CALL edge_seg_accumulate(nqty,nmat)
        ELSE IF (nmat>0) THEN
          CALL parallel_matseg_comm(nqty,nmat)
        ELSE
          CALL parallel_seg_comm(nqty)
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     increment time and terminate.
!-----------------------------------------------------------------------
      CALL timer(timeend)
      time_seg = time_seg + timeend-timestart
      segcount = segcount + 1

      RETURN
      END SUBROUTINE edge_seg_network
!-----------------------------------------------------------------------
!     subprogram 7. edge_seg_accumulate.
!     accumulate all representations of edge segment information from
!     seam_in into seam_out.
!-----------------------------------------------------------------------
      SUBROUTINE edge_seg_accumulate(nqty,nmat)

      INTEGER(i4), INTENT(IN) :: nqty,nmat

      INTEGER(i4) :: iv,ib,jv,jb,im
!-----------------------------------------------------------------------
!     accumulate vector data or off-diagonal matrix elements.  symmetry
!     is no longer assumed for the matrix elements, and the previous
!     and next labels on one side of the block border are opposite of
!     what they are on the other.
!-----------------------------------------------------------------------
      DO ib=1,SIZE(seam)
        DO iv=1,seam(ib)%nvert
          jb=seam(ib)%segment(iv)%ptr(1)
          jv=seam(ib)%segment(iv)%ptr(2)
          IF (jb/=0) THEN
            IF (nmat>0) THEN
              DO im=1,nmat,2
                seam(ib)%segment(iv)%seam_mat_out(1:nqty,1:nqty,im)     &
     &            =seam(ib)%segment(iv)%seam_mat_out(1:nqty,1:nqty,im)  &
     &            +seam(jb)%segment(jv)%seam_mat_in(1:nqty,1:nqty,im+1)
                seam(ib)%segment(iv)%seam_mat_out(1:nqty,1:nqty,im+1)   &
     &            =seam(ib)%segment(iv)%seam_mat_out(1:nqty,1:nqty,im+1)&
     &            +seam(jb)%segment(jv)%seam_mat_in(1:nqty,1:nqty,im)
              ENDDO
            ELSE
              seam(ib)%segment(iv)%seam_out(1:nqty)                     &
     &          =seam(ib)%segment(iv)%seam_out(1:nqty)                  &
     &          +seam(jb)%segment(jv)%seam_in(1:nqty)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_seg_accumulate
!-----------------------------------------------------------------------
!     subprogram 8. edge_seg_comp_accumulate.
!     accumulate all representations of edge segment information from
!     seam_cin into seam_cout.  complex version of edge_seg_accumulate
!-----------------------------------------------------------------------
      SUBROUTINE edge_seg_comp_accumulate(nqty,nmat)

      INTEGER(i4), INTENT(IN) :: nqty,nmat

      INTEGER(i4) :: iv,ib,jv,jb,im
!-----------------------------------------------------------------------
!     accumulate vector data or off-diagonal matrix elements.  symmetry
!     is no longer assumed for the matrix elements, and the previous
!     and next labels on one side of the block border are opposite of
!     what they are on the other.
!-----------------------------------------------------------------------
      DO ib=1,SIZE(seam)
        DO iv=1,seam(ib)%nvert
          jb=seam(ib)%segment(iv)%ptr(1)
          jv=seam(ib)%segment(iv)%ptr(2)
          IF (jb/=0) THEN
            IF (nmat>0) THEN
              DO im=1,nmat,2
                seam(ib)%segment(iv)%seam_mat_cout(1:nqty,1:nqty,im)    &
     &           =seam(ib)%segment(iv)%seam_mat_cout(1:nqty,1:nqty,im)  &
     &           +seam(jb)%segment(jv)%seam_mat_cin(1:nqty,1:nqty,im+1)
                seam(ib)%segment(iv)%seam_mat_cout(1:nqty,1:nqty,im+1)  &
     &           =seam(ib)%segment(iv)%seam_mat_cout(1:nqty,1:nqty,im+1)&
     &           +seam(jb)%segment(jv)%seam_mat_cin(1:nqty,1:nqty,im)
              ENDDO
            ELSE
              seam(ib)%segment(iv)%seam_cout(1:nqty)                    &
     &          =seam(ib)%segment(iv)%seam_cout(1:nqty)                 &
     &          +seam(jb)%segment(jv)%seam_cin(1:nqty)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_seg_comp_accumulate
!-----------------------------------------------------------------------
!     subprogram 9. edge_load_arr
!     copies real block data into seam_in before block-to-block
!     communication.
!-----------------------------------------------------------------------
      SUBROUTINE edge_load_arr(vec,nqty,n_side,seam)

      TYPE(vector_type), INTENT(IN) :: vec
      INTEGER(i4), INTENT(IN) :: nqty,n_side
      TYPE(edge_type), INTENT(INOUT) :: seam

      INTEGER(i4) :: ix,iy,iv,ib,iq,is
!-----------------------------------------------------------------------
!     copy block-internal data to seam storage.  vertex-centered data
!     first.
!-----------------------------------------------------------------------
      DO iv=1,seam%nvert
        seam%vertex(iv)%seam_in(1:nqty)=                                &
     &    vec%arr(1:nqty,seam%vertex(iv)%intxy(1),                      &
     &                   seam%vertex(iv)%intxy(2))
      ENDDO
!-----------------------------------------------------------------------
!     element side-centered data.
!-----------------------------------------------------------------------
      IF (n_side/=0) THEN
        DO iv=1,seam%nvert
          ix=seam%segment(iv)%intxys(1)
          iy=seam%segment(iv)%intxys(2)
!-PRE     this will be removed when higher-order triangles work.
          IF (ix<0) THEN
            seam%segment(iv)%seam_in=0
            CYCLE
          ENDIF
          iq=1
          IF (seam%segment(iv)%h_side) THEN
            DO is=1,n_side
              seam%segment(iv)%seam_in(iq:nqty+iq-1)=                   &
     &          vec%arrh(1:nqty,is,ix,iy)
              iq=iq+nqty
            ENDDO
          ELSE
            DO is=1,n_side
              seam%segment(iv)%seam_in(iq:nqty+iq-1)=                   &
     &          vec%arrv(1:nqty,is,ix,iy)
              iq=iq+nqty
            ENDDO
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_load_arr
!-----------------------------------------------------------------------
!     subprogram 10. edge_load_carr
!     copies complex block data into seam_cin before block-to-block
!     communication.
!-----------------------------------------------------------------------
      SUBROUTINE edge_load_carr(vec,nqty,ifs,ife,n_side,seam)

      TYPE(cvector_type), INTENT(IN) :: vec
      INTEGER(i4), INTENT(IN) :: nqty,n_side,ifs,ife
      TYPE(edge_type), INTENT(INOUT) :: seam

      INTEGER(i4) :: ix,iy,iv,ib,iq,is,im
!-----------------------------------------------------------------------
!     copy block-internal data to seam storage.  vertex-centered data
!     first.
!-----------------------------------------------------------------------
      DO iv=1,seam%nvert
        iq=1
        DO im=ifs,ife
          seam%vertex(iv)%seam_cin(iq:nqty+iq-1)=                       &
     &      vec%arr(1:nqty,seam%vertex(iv)%intxy(1),                    &
     &                     seam%vertex(iv)%intxy(2),im)
          iq=iq+nqty
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     element side-centered data.
!-----------------------------------------------------------------------
      IF (n_side/=0) THEN
        DO iv=1,seam%nvert
          ix=seam%segment(iv)%intxys(1)
          iy=seam%segment(iv)%intxys(2)
!-PRE     this will be removed when higher-order triangles work.
          IF (ix<0) THEN
            seam%segment(iv)%seam_cin=0
            CYCLE
          ENDIF
          iq=1
          IF (seam%segment(iv)%h_side) THEN
            DO im=ifs,ife
              DO is=1,n_side
                seam%segment(iv)%seam_cin(iq:nqty+iq-1)=                &
     &            vec%arrh(1:nqty,is,ix,iy,im)
                iq=iq+nqty
              ENDDO
            ENDDO
          ELSE
            DO im=ifs,ife
              DO is=1,n_side
                seam%segment(iv)%seam_cin(iq:nqty+iq-1)=                &
     &            vec%arrv(1:nqty,is,ix,iy,im)
                iq=iq+nqty
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_load_carr
!-----------------------------------------------------------------------
!     subprogram 11. edge_load_2D_carr
!     copies 2D complex block data into seam_cin before block-to-block
!     communication.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.
!-----------------------------------------------------------------------
      SUBROUTINE edge_load_2D_carr(vec,nqty,n_side,seam)

      TYPE(cvector_2D_type), INTENT(IN) :: vec
      INTEGER(i4), INTENT(IN) :: nqty,n_side
      TYPE(edge_type), INTENT(INOUT) :: seam

      INTEGER(i4) :: ix,iy,iv,ib,iq,is
!-----------------------------------------------------------------------
!     copy block-internal data to seam storage.  vertex-centered data
!     first.
!-----------------------------------------------------------------------
      DO iv=1,seam%nvert
        seam%vertex(iv)%seam_cin(1:nqty)=                               &
     &    vec%arr(1:nqty,seam%vertex(iv)%intxy(1),                      &
     &                   seam%vertex(iv)%intxy(2))
      ENDDO
!-----------------------------------------------------------------------
!     element side-centered data.
!-----------------------------------------------------------------------
      IF (n_side/=0) THEN
        DO iv=1,seam%nvert
          ix=seam%segment(iv)%intxys(1)
          iy=seam%segment(iv)%intxys(2)
!-PRE     this will be removed when higher-order triangles work.
          IF (ix<0) THEN
            seam%segment(iv)%seam_cin=0
            CYCLE
          ENDIF
          iq=1
          IF (seam%segment(iv)%h_side) THEN
            DO is=1,n_side
              seam%segment(iv)%seam_cin(iq:nqty+iq-1)=                  &
     &          vec%arrh(1:nqty,is,ix,iy)
              iq=iq+nqty
              ENDDO
          ELSE
            DO is=1,n_side
              seam%segment(iv)%seam_cin(iq:nqty+iq-1)=                  &
     &          vec%arrv(1:nqty,is,ix,iy)
              iq=iq+nqty
            ENDDO
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_load_2D_carr
!-----------------------------------------------------------------------
!     subprogram 12. edge_unload_arr
!     copies real data from seam_out to a block array after
!     block-to-block communication.
!-----------------------------------------------------------------------
      SUBROUTINE edge_unload_arr(vec,nqty,n_side,seam)

      TYPE(vector_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: nqty,n_side
      TYPE(edge_type), INTENT(IN) :: seam

      INTEGER(i4) :: ix,iy,iv,ib,iq,is
!-----------------------------------------------------------------------
!     copy from seam storage to block-internal data. vertex-centered
!     data first.
!-----------------------------------------------------------------------
      DO iv=1,seam%nvert
        vec%arr(1:nqty,seam%vertex(iv)%intxy(1),                        &
     &                 seam%vertex(iv)%intxy(2))=                       &
     &    seam%vertex(iv)%seam_out(1:nqty)
      ENDDO
!-----------------------------------------------------------------------
!     element side-centered data.
!-----------------------------------------------------------------------
      IF (n_side/=0) THEN
        DO iv=1,seam%nvert
          ix=seam%segment(iv)%intxys(1)
          iy=seam%segment(iv)%intxys(2)
!-PRE     this will be removed when higher-order triangles work.
          IF (ix<0) CYCLE
          iq=1
          IF (seam%segment(iv)%h_side) THEN
            DO is=1,n_side
              vec%arrh(1:nqty,is,ix,iy)=                                &
     &          seam%segment(iv)%seam_out(iq:nqty+iq-1)
              iq=iq+nqty
            ENDDO
          ELSE
            DO is=1,n_side
              vec%arrv(1:nqty,is,ix,iy)=                                &
     &          seam%segment(iv)%seam_out(iq:nqty+iq-1)
              iq=iq+nqty
            ENDDO
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_unload_arr
!-----------------------------------------------------------------------
!     subprogram 13. edge_unload_carr
!     copies complex data from seam_cout to a block array after
!     block-to-block communication.
!-----------------------------------------------------------------------
      SUBROUTINE edge_unload_carr(vec,nqty,ifs,ife,n_side,seam)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: nqty,n_side,ifs,ife
      TYPE(edge_type), INTENT(IN) :: seam

      INTEGER(i4) :: ix,iy,iv,ib,iq,is,im
!-----------------------------------------------------------------------
!     copy from seam storage to block-internal data. vertex-centered
!     data first.
!-----------------------------------------------------------------------
      DO iv=1,seam%nvert
        iq=1
        DO im=ifs,ife
          vec%arr(1:nqty,seam%vertex(iv)%intxy(1),                      &
     &                   seam%vertex(iv)%intxy(2),im)=                  &
     &      seam%vertex(iv)%seam_cout(iq:nqty+iq-1)
          iq=iq+nqty
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     element side-centered data.
!-----------------------------------------------------------------------
      IF (n_side/=0) THEN
        DO iv=1,seam%nvert
          ix=seam%segment(iv)%intxys(1)
          iy=seam%segment(iv)%intxys(2)
!-PRE     this will be removed when higher-order triangles work.
          IF (ix<0) CYCLE
          iq=1
          IF (seam%segment(iv)%h_side) THEN
            DO im=ifs,ife
              DO is=1,n_side
                vec%arrh(1:nqty,is,ix,iy,im)=                           &
     &            seam%segment(iv)%seam_cout(iq:nqty+iq-1)
                iq=iq+nqty
              ENDDO
            ENDDO
          ELSE
            DO im=ifs,ife
              DO is=1,n_side
                vec%arrv(1:nqty,is,ix,iy,im)=                           &
     &            seam%segment(iv)%seam_cout(iq:nqty+iq-1)
                iq=iq+nqty
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_unload_carr
!-----------------------------------------------------------------------
!     subprogram 14. edge_unload_2D_carr
!     copies complex data from seam_cout to a 2D block array after
!     block-to-block communication.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.
!-----------------------------------------------------------------------
      SUBROUTINE edge_unload_2D_carr(vec,nqty,n_side,seam)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: nqty,n_side
      TYPE(edge_type), INTENT(IN) :: seam

      INTEGER(i4) :: ix,iy,iv,ib,iq,is
!-----------------------------------------------------------------------
!     copy from seam storage to block-internal data. vertex-centered
!     data first.
!-----------------------------------------------------------------------
      DO iv=1,seam%nvert
        vec%arr(1:nqty,seam%vertex(iv)%intxy(1),                        &
     &                 seam%vertex(iv)%intxy(2))=                       &
     &    seam%vertex(iv)%seam_cout(1:nqty)
      ENDDO
!-----------------------------------------------------------------------
!     element side-centered data.
!-----------------------------------------------------------------------
      IF (n_side/=0) THEN
        DO iv=1,seam%nvert
          ix=seam%segment(iv)%intxys(1)
          iy=seam%segment(iv)%intxys(2)
!-PRE     this will be removed when higher-order triangles work.
          IF (ix<0) CYCLE
          iq=1
          IF (seam%segment(iv)%h_side) THEN
            DO is=1,n_side
              vec%arrh(1:nqty,is,ix,iy)=                                &
     &          seam%segment(iv)%seam_cout(iq:nqty+iq-1)
              iq=iq+nqty
            ENDDO
          ELSE
            DO is=1,n_side
              vec%arrv(1:nqty,is,ix,iy)=                                &
     &          seam%segment(iv)%seam_cout(iq:nqty+iq-1)
              iq=iq+nqty
            ENDDO
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_unload_2D_carr
!-----------------------------------------------------------------------
!     subprogram 15. edge_load_save
!     copies real block data into seam_save.  otherwise, this is
!     similar to edge_load_arr.
!-----------------------------------------------------------------------
      SUBROUTINE edge_load_save(vec,nqty,n_side,seam)

      TYPE(vector_type), INTENT(IN) :: vec
      INTEGER(i4), INTENT(IN) :: nqty,n_side
      TYPE(edge_type), INTENT(INOUT) :: seam

      INTEGER(i4) :: ix,iy,iv,ib,iq,is
!-----------------------------------------------------------------------
!     copy block-internal data to seam storage.  vertex-centered data
!     first.
!-----------------------------------------------------------------------
      DO iv=1,seam%nvert
        seam%vertex(iv)%seam_save(1:nqty)=                              &
     &    vec%arr(1:nqty,seam%vertex(iv)%intxy(1),                      &
     &                   seam%vertex(iv)%intxy(2))
      ENDDO
!-----------------------------------------------------------------------
!     element side-centered data.
!-----------------------------------------------------------------------
      IF (n_side/=0) THEN
        DO iv=1,seam%nvert
          ix=seam%segment(iv)%intxys(1)
          iy=seam%segment(iv)%intxys(2)
!-PRE     this will be removed when higher-order triangles work.
          IF (ix<0) THEN
            seam%segment(iv)%seam_save=0
            CYCLE
          ENDIF
          iq=1
          IF (seam%segment(iv)%h_side) THEN
            DO is=1,n_side
              seam%segment(iv)%seam_save(iq:nqty+iq-1)=                 &
     &          vec%arrh(1:nqty,is,ix,iy)
              iq=iq+nqty
            ENDDO
          ELSE
            DO is=1,n_side
              seam%segment(iv)%seam_save(iq:nqty+iq-1)=                 &
     &          vec%arrv(1:nqty,is,ix,iy)
              iq=iq+nqty
            ENDDO
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_load_save
!-----------------------------------------------------------------------
!     subprogram 16. edge_load_csave
!     copies complex block data into seam_csave.  otherwise, this is
!     similar to edge_load_carr.
!-----------------------------------------------------------------------
      SUBROUTINE edge_load_csave(vec,nqty,ifs,ife,n_side,seam)

      TYPE(cvector_type), INTENT(IN) :: vec
      INTEGER(i4), INTENT(IN) :: nqty,n_side,ifs,ife
      TYPE(edge_type), INTENT(INOUT) :: seam

      INTEGER(i4) :: ix,iy,iv,ib,iq,is,im
!-----------------------------------------------------------------------
!     copy block-internal data to seam storage.  vertex-centered data
!     first.
!-----------------------------------------------------------------------
      DO iv=1,seam%nvert
        iq=1
        DO im=ifs,ife
          seam%vertex(iv)%seam_csave(iq:nqty+iq-1)=                     &
     &      vec%arr(1:nqty,seam%vertex(iv)%intxy(1),                    &
     &                     seam%vertex(iv)%intxy(2),im)
          iq=iq+nqty
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     element side-centered data.
!-----------------------------------------------------------------------
      IF (n_side/=0) THEN
        DO iv=1,seam%nvert
          ix=seam%segment(iv)%intxys(1)
          iy=seam%segment(iv)%intxys(2)
!-PRE     this will be removed when higher-order triangles work.
          IF (ix<0) THEN
            seam%segment(iv)%seam_csave=0
            CYCLE
          ENDIF
          iq=1
          IF (seam%segment(iv)%h_side) THEN
            DO im=ifs,ife
              DO is=1,n_side
                seam%segment(iv)%seam_csave(iq:nqty+iq-1)=              &
     &            vec%arrh(1:nqty,is,ix,iy,im)
                iq=iq+nqty
              ENDDO
            ENDDO
          ELSE
            DO im=ifs,ife
              DO is=1,n_side
                seam%segment(iv)%seam_csave(iq:nqty+iq-1)=              &
     &            vec%arrv(1:nqty,is,ix,iy,im)
                iq=iq+nqty
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_load_csave
!-----------------------------------------------------------------------
!     subprogram 17. edge_add_save
!     adds real data from seam_save to a vector structure.  otherwise,
!     this is similar to edge_unload_arr.
!-----------------------------------------------------------------------
      SUBROUTINE edge_add_save(vec,nqty,n_side,seam)

      TYPE(vector_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: nqty,n_side
      TYPE(edge_type), INTENT(IN) :: seam

      INTEGER(i4) :: ix,iy,iv,ib,iq,is
!-----------------------------------------------------------------------
!     copy from seam storage to block-internal data. vertex-centered
!     data first.
!-----------------------------------------------------------------------
      DO iv=1,seam%nvert
        vec%arr(1:nqty,seam%vertex(iv)%intxy(1),                        &
     &                 seam%vertex(iv)%intxy(2))=                       &
     &    vec%arr(1:nqty,seam%vertex(iv)%intxy(1),                      &
     &                   seam%vertex(iv)%intxy(2))+                     &
     &    seam%vertex(iv)%seam_save(1:nqty)
      ENDDO
!-----------------------------------------------------------------------
!     element side-centered data.
!-----------------------------------------------------------------------
      IF (n_side/=0) THEN
        DO iv=1,seam%nvert
          ix=seam%segment(iv)%intxys(1)
          iy=seam%segment(iv)%intxys(2)
!-PRE     this will be removed when higher-order triangles work.
          IF (ix<0) CYCLE
          iq=1
          IF (seam%segment(iv)%h_side) THEN
            DO is=1,n_side
              vec%arrh(1:nqty,is,ix,iy)=vec%arrh(1:nqty,is,ix,iy)+      &
     &          seam%segment(iv)%seam_save(iq:nqty+iq-1)
              iq=iq+nqty
            ENDDO
          ELSE
            DO is=1,n_side
              vec%arrv(1:nqty,is,ix,iy)=vec%arrv(1:nqty,is,ix,iy)+      &
     &          seam%segment(iv)%seam_save(iq:nqty+iq-1)
              iq=iq+nqty
            ENDDO
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_add_save
!-----------------------------------------------------------------------
!     subprogram 18. edge_add_csave
!     adds complex data from seam_csave to a complex vector structure.
!     otherwise, this is similar to edge_unload_carr.
!-----------------------------------------------------------------------
      SUBROUTINE edge_add_csave(vec,nqty,ifs,ife,n_side,seam)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: nqty,n_side,ifs,ife
      TYPE(edge_type), INTENT(IN) :: seam

      INTEGER(i4) :: ix,iy,iv,ib,iq,is,im
!-----------------------------------------------------------------------
!     copy from seam storage to block-internal data. vertex-centered
!     data first.
!-----------------------------------------------------------------------
      DO iv=1,seam%nvert
        iq=1
        DO im=ifs,ife
          vec%arr(1:nqty,seam%vertex(iv)%intxy(1),                      &
     &                   seam%vertex(iv)%intxy(2),im)=                  &
     &      vec%arr(1:nqty,seam%vertex(iv)%intxy(1),                    &
     &                     seam%vertex(iv)%intxy(2),im)+                &
     &      seam%vertex(iv)%seam_csave(iq:nqty+iq-1)
          iq=iq+nqty
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     element side-centered data.
!-----------------------------------------------------------------------
      IF (n_side/=0) THEN
        DO iv=1,seam%nvert
          ix=seam%segment(iv)%intxys(1)
          iy=seam%segment(iv)%intxys(2)
!-PRE     this will be removed when higher-order triangles work.
          IF (ix<0) CYCLE
          iq=1
          IF (seam%segment(iv)%h_side) THEN
            DO im=ifs,ife
              DO is=1,n_side
                vec%arrh(1:nqty,is,ix,iy,im)=                           &
     &            vec%arrh(1:nqty,is,ix,iy,im)+                         &
     &            seam%segment(iv)%seam_csave(iq:nqty+iq-1)
                iq=iq+nqty
              ENDDO
            ENDDO
          ELSE
            DO im=ifs,ife
              DO is=1,n_side
                vec%arrv(1:nqty,is,ix,iy,im)=                           &
     &            vec%arrv(1:nqty,is,ix,iy,im)+                         &
     &            seam%segment(iv)%seam_csave(iq:nqty+iq-1)
                iq=iq+nqty
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edge_add_csave
!-----------------------------------------------------------------------
!     subprogram 19. edge_zero_save
!     set the seam_save data arrays to 0.
!-----------------------------------------------------------------------
      SUBROUTINE edge_zero_save(seam)

      TYPE(edge_type), INTENT(INOUT) :: seam

      INTEGER(i4) :: iv
!-----------------------------------------------------------------------
!     loop over vertex and segment structures to clear the seam_save
!     arrays.
!-----------------------------------------------------------------------
      DO iv=1,seam%nvert
        seam%vertex(iv)%seam_save=0._r8
        seam%segment(iv)%seam_save=0._r8
      ENDDO

      RETURN
      END SUBROUTINE edge_zero_save
!-----------------------------------------------------------------------
!     subprogram 20. edge_zero_csave
!     set the seam_csave data arrays to 0.
!-----------------------------------------------------------------------
      SUBROUTINE edge_zero_csave(seam)

      TYPE(edge_type), INTENT(INOUT) :: seam

      INTEGER(i4) :: iv
!-----------------------------------------------------------------------
!     loop over vertex and segment structures to clear the seam_csave
!     arrays.
!-----------------------------------------------------------------------
      DO iv=1,seam%nvert
        seam%vertex(iv)%seam_csave=0._r8
        seam%segment(iv)%seam_csave=0._r8
      ENDDO

      RETURN
      END SUBROUTINE edge_zero_csave
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      END MODULE edge
