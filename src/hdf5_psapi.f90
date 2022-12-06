!-----------------------------------------------------------------------
! file hdf_psapi
! hdf_api module
!  Very generic module meant for writing HDF5 files with particular
!  attributes.  In case we want to convert over to C, it should make
!  this easier.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! code organization for hdf_psapi.
!-----------------------------------------------------------------------
! 4.  vshdf5_initprovenance(vsprov)
! 1 . write_provenance(gInId,grName,prin,h5err)
!-----------------------------------------------------------------------
! module hdf_api
!-----------------------------------------------------------------------
  module hdf5_psapi
  use hdf5_api
!-----------------------------------------------------------------------
! Provenance data for the codes.  See:
!  https://ice.txcorp.com/trac/vizschema/wiki/ProvenanceMetaData
! Some of these are hard to do in fortran
!-----------------------------------------------------------------------
  type provenancedata
   !att vstype = "runinfo"
   character(len=30) :: software        ! code name
   character(len=30) :: swVersion       ! software version
   character(len=30) :: swRevision      ! revision number from subversion
   character(len=30) :: vsVersion       ! vs compliance version
   character(len=30) :: fcCompiler      ! fortran compiler
   character(len=30) :: fcCompilerVersion      ! fc version
   character(len=30) :: fcCompilerFlags        ! fc flags
   character(len=30) :: buildHost              ! build host
   character(len=30) :: buildHostType          ! host type
   character(len=30) :: buildDate              ! Build date
   character(len=30) :: runHost                ! run host
   character(len=30) :: runHostType            ! run host type
   character(len=30) :: user                   ! user
   character(len=30) :: runDate                ! run date
   character(len=30) :: commandLine            ! commandline
  end type
  contains
!-----------------------------------------------------------------------
! subprogram 0.1 vshdf5_initprovenance
! Initialize these variables to default values.
!-----------------------------------------------------------------------
  subroutine vshdf5_initprovenance(vsprov)
  type(provenancedata), intent(inout) :: vsprov
  
  vsprov%software=" "    ! code name
  vsprov%swVersion=" "   ! software version
  vsprov%swRevision=" "  ! revision number from subversion
  vsprov%vsVersion="3.0.0"   ! VS compliance version
  vsprov%fcCompiler=" "  ! fortran compiler
  vsprov%fcCompilerVersion=" "  ! FC version
  vsprov%fcCompilerFlags=" "    ! FC flags
  vsprov%buildHost=" "          ! Build host
  vsprov%buildHostType=" "      ! Host type
  vsprov%buildDate=" "          ! Build host
  vsprov%runHost=" "            ! Date of build
  vsprov%runHostType=" "        ! Run Host type
  vsprov%runDate=" "            ! run date
  vsprov%user=" "               ! user
  vsprov%commandLine=" "        ! commandline
  return
  end subroutine vshdf5_initprovenance
!-----------------------------------------------------------------------
! subprogram 20. write_provenance
! Make a group that contains the provenance data.  See:
!  https://ice.txcorp.com/trac/vizschema/wiki/ProvenanceMetaData
!-----------------------------------------------------------------------
  subroutine write_provenance(gInId,grName,prin,h5in,h5err)
  integer(HID_T), intent(in) :: gInId
  character*(*), intent(in) :: grName
  type(provenanceData), intent(in) :: prin
  type(hdf5InOpts), intent(in) :: h5in
  type(hdf5ErrorType), intent(inoUT) :: h5err
  integer,parameter :: FAIL=-1
  integer(HID_T) :: prId
  integer :: error
!-----------------------------------------------------------------------
! Open the group
!-----------------------------------------------------------------------
  call make_group(gInId, grName, prId,h5in,h5err)
!-----------------------------------------------------------------------
!     Add the VisSchema attributes for provenance data
!  Att vsType = "runInfo"                             // Required string
!  Att vsSoftware = "VORPAL"                          // Required string
!  Att vsSwVersion = "3.1.0"                          // Required string
!  Att vsVsVersion = "2.1.0"                          // Required string
!  Att vsSwRevision = "r10168"                        // Optional string
!  Att vsCxxCompiler = "g++"                          // Optional string
!  Att vsCxxCompilerVersion = "4.3.0"                 // Optional string
!  Att vsCxxCompilerFlags = "-O2 -g"                  // Optional string
!  Att vsBuildHost = "boron.txcorp.com"               // Optional string
!  Att vsBuildHostType = "x86_64-unknown-linux-gnu"   // Optional string
!  Att vsRunHost = "boron.txcorp.com"                 // Optional string
!  Att vsUser = "cary"                                // Optional string
!  Att vsRunDate = "Mon Nov 10 09:06:14 MST 2008"     // Optional string
!  Att vsCommandLine = "vorpal -i run.in"             // Optional string
!-----------------------------------------------------------------------
  call write_attribute(prId,'vsType',"runInfo",h5in,h5err)
  if(LEN_TRIM(prin%software)>0)  &
   call write_attribute(prId,'vsSoftware',prin%software,h5in,h5err)
  if(LEN_TRIM(prin%swVersion)>0)  &
    call write_attribute(prId,'vsSwVersion',prin%swVersion,h5in,h5err)
  if(LEN_TRIM(prin%swRevision)>0)  &
    call write_attribute(prId,'vsSwRevision',prin%swRevision,h5in,h5err)
  if(LEN_TRIM(prin%vsVersion)>0)  &
    call write_attribute(prId,'vsVsVersion',prin%vsVersion,h5in,h5err)
  if(LEN_TRIM(prin%fcCompiler)>0)  &
    call write_attribute(prId,'vsFcCompiler',prin%fcCompiler,h5in,h5err)
  if(LEN_TRIM(prin%fcCompilerVersion)>0)  &
    call write_attribute(prId,'vsFcCompilerVersion',prin%fcCompilerVersion,h5in,h5err)
  if(LEN_TRIM(prin%fcCompilerFlags)>0)  &
    call write_attribute(prId,'vsFcCompilerFlags',prin%fcCompilerFlags,h5in,h5err)
  if(LEN_TRIM(prin%buildHost)>0)  &
    call write_attribute(prId,'vsBuildHost',prin%buildHost,h5in,h5err)
  if(LEN_TRIM(prin%buildHostType)>0)  &
    call write_attribute(prId,'vsBuildHostType',prin%buildHostType,h5in,h5err)
  if(LEN_TRIM(prin%runHost)>0)  &
    call write_attribute(prId,'vsRunHost',prin%runHost,h5in,h5err)
  if(LEN_TRIM(prin%runHostType)>0)  &
    call write_attribute(prId,'vsRunHostType',prin%runHostType,h5in,h5err)
  if(LEN_TRIM(prin%user)>0)  &
    call write_attribute(prId,'vsUser',prin%user,h5in,h5err)
  if(LEN_TRIM(prin%runHost)>0)  &
    call write_attribute(prId,'vsRunHost',prin%runHost,h5in,h5err)
  if(LEN_TRIM(prin%runDate)>0)  &
    call write_attribute(prId,'vsRunDate',prin%runDate,h5in,h5err)
  if(LEN_TRIM(prin%commandLine)>0)  &
    call write_attribute(prId,'vsCommandLine',prin%commandLine,h5in,h5err)
  call close_group(grName,prId,h5err)
  return
  end subroutine write_provenance

  end module hdf5_psapi

