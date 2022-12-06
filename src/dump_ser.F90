!-----------------------------------------------------------------------
!     $Id: dump_ser.F90 7592 2021-11-09 18:30:00Z held $
!     routines reading and writing data dumps.                          
!-----------------------------------------------------------------------
#include "config.f"
!-----------------------------------------------------------------------
!     code organization.                                                
!-----------------------------------------------------------------------                                               
!     2.  dump_read.                                                                                         
!-----------------------------------------------------------------------
      MODULE dump
      IMPLICIT NONE
      CONTAINS
!-----------------------------------------------------------------------
!     subprogram 2. dump_read.                                                                                        
!-----------------------------------------------------------------------
      SUBROUTINE dump_read(nmodes,nmodes_total,keff,keff_total,         &
     &                     dump_time,dump_step) 
      USE input 
      USE io
      USE fields 
      USE seam_storage_mod 
      USE time
      USE normalize_mod
      USE rt_swap_mod
      USE physdat, ONLY: twopi
      IMPLICIT NONE 
                                                                       
      INTEGER(i4), INTENT(IN) :: nmodes_total
      INTEGER(i4), INTENT(OUT) :: nmodes,dump_step
      REAL(r8), INTENT(OUT) :: dump_time 
      REAL(r8), DIMENSION(:), POINTER :: keff,keff_total
                                                                        
      INTEGER(i4) :: ixbl,ib,nmodes_dump,ierror,tmp
      REAL(r8) :: iread 
      LOGICAL :: file_stat 
#ifdef HAVE_FC_HDF5
      INTEGER(HID_T) :: sid,rbid
#endif
      IF (ASSOCIATED(keff)) DEALLOCATE(keff)
      IF (ASSOCIATED(keff_total)) DEALLOCATE(keff_total)
      ALLOCATE(keff(nmodes_total),keff_total(nmodes_total))
      CALL timer(timestart)
!-----------------------------------------------------------------------
!     use h5 format if requested. 
!-----------------------------------------------------------------------
      IF (h5dump) THEN
#ifdef HAVE_FC_HDF5
        CALL fch5init
        INQUIRE(FILE=TRIM(dump_file),EXIST=file_stat) 
        IF (.NOT.file_stat) CALL nim_stop                               &
     &    ('Dump file '//TRIM(dump_file)//' does not exist.')
        CALL open_oldh5file(TRIM(dump_file),fileid,rootgid,h5in,h5err)
!-----------------------------------------------------------------------
!       read global data.
!-----------------------------------------------------------------------
        CALL open_group(rootgid,"dumpTime",sid,h5err)
        CALL read_attribute(sid,"vsTime",dump_time,h5in,h5err)
        CALL read_attribute(sid,"vsStep",dump_step,h5in,h5err)
        CALL close_group("dumpTime",sid,h5err)
        CALL read_attribute(rootgid,"nbl",nbl_total,h5in,h5err)
        CALL read_attribute_int_sc(rootgid,"nrbl",nrbl_total,h5in,h5err)
        nbl=nbl_total; nrbl=nrbl_total
        CALL read_attribute(rootgid,"poly_degree",tmp,h5in,h5err)
        IF (tmp/=poly_degree) THEN 
           CALL nim_stop                                                &
     &       ('NIMROD input is not consistent with dumped poly_degree.')
        ENDIF 
        CALL read_attribute(rootgid,"nmodes",tmp,h5in,h5err)
        IF (tmp/=nmodes_total) THEN 
           CALL nim_stop                                                &
     &       ('NIMROD input is not consistent with dumped nmodes.')     
        ENDIF 
        CALL read_h5(rootgid,"keff",keff,h5in,h5err)
        IF (tmp>1) THEN
          IF (geom=='lin'.AND.                                          &
     &     NINT((keff(tmp)-keff(tmp-1))/twopi*per_length)/=zperiod) THEN 
             CALL nim_stop                                              &
     &         ('NIMROD input is not consistent with dumped zperiod.')
          ELSEIF (geom=='tor'.AND.                                      &
     &                        NINT(keff(tmp)-keff(tmp-1))/=zperiod) THEN 
             CALL nim_stop                                              &
     &         ('NIMROD input is not consistent with dumped zperiod.')
          ENDIF     
        ENDIF
        keff_total=keff
!-----------------------------------------------------------------------
!       read seam and rblock data.
!-----------------------------------------------------------------------
        CALL open_group(rootgid,"seams",sid,h5err)
        CALL open_group(rootgid,"rblocks",rbid,h5err)
        ALLOCATE(seam0,seam(nbl),rb(nrbl)) 
        CALL h5_read_seam(seam0,0_i4,sid)
        DO ib=1,nrbl
          CALL h5_read_seam(seam(ib),ib,sid)
          CALL h5_read_rblock(rb(ib),ib,rbid)
        ENDDO
        CALL close_group("seams",sid,h5err)
        CALL close_group("rblocks",rbid,h5err)
        CALL close_h5file(fileid,rootgid,h5err)
!-----------------------------------------------------------------------
!       PRE tblocks
!-----------------------------------------------------------------------
#else
        CALL nim_stop('FC_HDF5 not linked, set h5dump=F')
#endif /* HAVE_FC_HDF5 */
      ELSE
!-----------------------------------------------------------------------
!       open dump file.
!-----------------------------------------------------------------------
        INQUIRE(FILE=TRIM(dump_file),EXIST=file_stat) 
        IF (.NOT.file_stat) CALL nim_stop                               &
     &    ('Dump file '//TRIM(dump_file)//' does not exist.')           
!-----------------------------------------------------------------------
!       read global data.  integers are read into 64 bit reals then 
!       converted upon copying.                         
!-----------------------------------------------------------------------
        CALL open_bin(rstrt_unit,TRIM(dump_file),'OLD','REWIND',64_i4) 
        READ(rstrt_unit) dump_time 
        READ(rstrt_unit) iread 
        dump_step=NINT(iread) 
        READ(rstrt_unit) iread 
        nbl_total=NINT(iread) 
        READ(rstrt_unit) iread 
        nrbl_total=NINT(iread) 
        nbl=nbl_total; nrbl=nrbl_total
!-----------------------------------------------------------------------
!       check for basis consistency.                                    
!-----------------------------------------------------------------------
        READ(rstrt_unit) iread 
        IF (NINT(iread)/=poly_degree) THEN 
           CALL nim_stop                                                &
     &       ('NIMROD input is not consistent with dumped poly_degree.')
        ENDIF 
!-----------------------------------------------------------------------
!       check for fourier rep. consistency, then read wavenumber array. 
!-----------------------------------------------------------------------
        READ(rstrt_unit) iread 
        nmodes_dump=NINT(iread) 
        IF (nmodes_dump/=nmodes_total) THEN 
           CALL nim_stop                                                &
     &       ('NIMROD input is not consistent with dumped nmodes.')     
        ENDIF 
        READ(rstrt_unit) keff 
        IF (nmodes_dump>1.AND.geom=='tor') THEN
          IF (NINT(keff(nmodes_dump)-keff(nmodes_dump-1))/=zperiod) THEN 
             CALL nim_stop                                              &
     &         ('NIMROD input is not consistent with dumped zperiod.')
          ENDIF     
        ENDIF
        keff_total=keff 
!-----------------------------------------------------------------------
!       read seam data.                                                   
!-----------------------------------------------------------------------
        ALLOCATE(seam0,seam(nbl)) 
        CALL dump_read_seam(seam0) 
        DO ib=1,nbl 
           CALL dump_read_seam(seam(ib)) 
        ENDDO 
!-----------------------------------------------------------------------
!       read rblock data.                                                 
!-----------------------------------------------------------------------
        ALLOCATE(rb(nrbl)) 
        DO ib=1,nrbl 
           CALL dump_read_rblock(rb(ib)) 
        ENDDO 
!-----------------------------------------------------------------------
!       read tblock data.                                                 
!-----------------------------------------------------------------------
        ALLOCATE(tb(nrbl+1:nbl)) 
!-PRE   DO ib=nrbl+1,nbl 
!-PRE           CALL dump_read_tblock(tb(ib))                             
!-PRE   ENDDO 
!-----------------------------------------------------------------------
!       this is a place holder for detecting closure i/o type             
!-----------------------------------------------------------------------
        read(rstrt_unit,iostat=ierror) 
        IF (ierror==0) write(*,*)'Charlson will fix this' 
!-----------------------------------------------------------------------
!       close dump file.                                                  
!-----------------------------------------------------------------------
        CALL close_bin(rstrt_unit,TRIM(dump_file)) 
      ENDIF
!-----------------------------------------------------------------------
!     allocate mxi for nonuniform domain decompositions.
!-----------------------------------------------------------------------
      DO ib=1,nrbl
        ALLOCATE(rb(ib)%mxi(1:nxbl))
      ENDDO
      DO ixbl=1,nxbl
        DO ib=1,nrbl
          rb(ib)%mxi(ixbl)=INT(mx*ixbl/nxbl)
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!     allocate parallel arrays even though this is the serial interface
!-----------------------------------------------------------------------
      nlayers=1 ! serial programs use dump_read
      CALL parallel_block_init(nmodes,nmodes_total,nlayers,decompflag)
      CALL parallel_block_init2
!-----------------------------------------------------------------------
!     real time transfer eq
!-----------------------------------------------------------------------
      IF (rt_transfer_eq.AND.nonlinear) CALL eq_swap(nmodes,keff)
!-----------------------------------------------------------------------
!     normalize parameters
!-----------------------------------------------------------------------
      IF (normalize.AND.ofname=='nimrod.out')                           &
     &  CALL normalize_fields(dump_time)
!-----------------------------------------------------------------------
!     terminate.                                                        
!-----------------------------------------------------------------------
      CALL timer(timeend)
      time_io = time_io + timeend-timestart
      RETURN 
      END SUBROUTINE dump_read 
      END MODULE dump