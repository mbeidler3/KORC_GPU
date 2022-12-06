!-----------------------------------------------------------------------
!     $Id: input.F90 7642 2022-03-08 19:08:41Z held $
!     declarations and default values of input parameters.
!     note that new parameters should be added to the lists in module
!     input and to the namelist declarations in subroutine
!     read_namelist.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     code organization.
!-----------------------------------------------------------------------
!     1. input
!     2. read_namelist
!     3. rmcomment
!     4. get_k2min.
!     5. broadcast_input
!     6. bcast_str
!     7. input_normalization
!     8. direct_check.
!     9. perturbation_set
!-----------------------------------------------------------------------
!     subprogram 1. input
!     module containing the type declaration and default values for
!     all input variables.
!-----------------------------------------------------------------------
#include "config.f"
      MODULE input
      USE local
      !USE physdat
      IMPLICIT NONE

!=======================================================================
!     grid parameters.
!=======================================================================
!     The type of grid is set by gridshape and supporting input:
!       'rect':  rectangular cross section with the periodicity of the
!           x and y directions also specified as 'none', 'y-dir',
!           or 'both'.
!           The dimensions are set by xmin, xmax, ymin and ymax.
!           This grid may be made nonorthogonal by setting skew/=0.
!           To get a nonuniform grid, set firstx and firsty to the
!           first delta-x and delta-y desired, respectively,
!           otherwise set these parameters to 0.
!        'circ':  circular or annular cross section with minimum and
!           maximum polar radii set by xmin and xmax.
!           firstx and firsty have some functionality here, too.
!           for toroidal geometries (see geom below), firstx sets
!           an offset for the center of the grid from the geometric
!           center of the cross section to match a Shafranov shift.
!           firsty stretches the grid in the radial direction from
!           its center.  Unlike the application to rect gridshapes,
!           firsty is nondimensional, relative to uniform radial
!           gridding.
!        'flux':  shape based on equilibrium magnetic field which is
!           controlled by the fluxgrid.in namelist.  See
!           fluxgrid/input.f for details
!        'rect_cir': a physically smooth, i.e. circular or elliptical,
!       	 but logically rectangular grid with maximum dimensions
!       	 determined by xmin, xmax, ymin, and ymax.  The npc
!       	 numerical parameter controls the number of smoothing
!       	 passes.

      !REQUIRED
      CHARACTER(8) :: gridshape="circ"
      !REQUIRED
      CHARACTER(8) :: periodicity="none"
      !IF GRIDSHAPE='rect' or 'circ'
      REAL(r8) :: xmin=0
      REAL(r8) :: xmax=1
      REAL(r8) :: ymin=0
      REAL(r8) :: ymax=1
      REAL(r8) :: skew=0
      REAL(r8) :: firstx=0
      REAL(r8) :: firsty=0
      !IF grishape='tor'
      CHARACTER(64) :: eqfile="ldx.dat"

!-----------------------------------------------------------------------
!     Select a periodic linear system ('lin') with length specified by
!     per_length or a toroidal system ('tor').  For circular grids the
!     major radius is set by xo, and yo controls the vertical offset.
!     For flux grids, this is determined by the data, and for
!     rectangular grids, it's set by xmin, xmax, ymin and ymax.

      CHARACTER(3) :: geom="tor"
      REAL(r8) :: per_length=1
      REAL(r8) :: xo=2
      REAL(r8) :: yo=0
      REAL(r8) :: x1=2
      REAL(r8) :: y1=0

!-----------------------------------------------------------------------
!     Set zperiod to simulate only a fraction of a periodic device
!     (similar to that used in BOUT++). zperiod=1 simulates the full
!     device and zperiod=n simulates 1/n of the device.  This factor
!     just multiplies the wavenumber array keff in nimset.

      INTEGER(i4) :: zperiod=1

!-----------------------------------------------------------------------
!     The following parameters control the mesh refinement.  The two
!     parameters mx and my set the number of cells in the central
!     part of the grid.  For rectangular grids, mx and my are the
!     number of cells across the x and y directions.  For circular or
!     flux grids, they are the number of radial and poloidal cells,
!     respectively, in the polar region gridded by rectangles
!     (logically rectangular).  This region is annular in shape and
!     extends from a small hole at the magnetic axis to the separatrix
!     (or some fraction of the separatrix, depending on parameter
!     psihigh set in fluxgrid.in).  Note that mesh refinement set by
!     mx and my is independent of the     block decomposition set by nxb
!     and nybl (later), but each block in the central rblock region
!     must have at least 3 cells in each direction.

      INTEGER(i4) :: mx=32
      INTEGER(i4) :: my=32

!-----------------------------------------------------------------------
!     poly_degree is the degree of the polynomials used in the
!     lagrange elements for dependent variables.  The amount of
!     information is then given by (mx*poly_degree)*(my*poly_degree)

      INTEGER(i4) :: poly_degree=1

!-----------------------------------------------------------------------
!     A mapped grid always uses the same polynomial degree as other
!     elements, but in constructing the grid, one can use a different
!     ordering of polynomials in helping to construct the mesh.  One
!     always wants poly_degree_mesh < poly_degree to have convergence

      INTEGER(i4) :: poly_degree_mesh=1

!-----------------------------------------------------------------------
!     The default high-order element type uses Lagrange elements, which
!     has uniformly spaced node points within a cell.  By changing the
!     distribution of the nodes, then different polynomial elements
!     can be chosen.  Currently support, one can also choose
!     poly_distribution='gll'.  The gll polynomials can lead to better
!     matrices and avoid problems at high poly_degree.
!     When poly_distribution='gll' is chosen, the string 'gll' is
!     appended to the "dump_name" parameter to name the dump files
!     (see dump_name description below); i.e., dump files will be of the
!     form, dumpgll.00000

      CHARACTER(16) :: poly_distribution='uniform' ! 'gll'

!-----------------------------------------------------------------------
!     For grids that are based on Grad_Shafranov equilibria, most of the
!     parameters controlling the grid are in the fluxgrid.in file and
!     need to be controlled there; however, one important parameter is
!     how many points to place in the "vacuum" region (see
!     fluxgrid/input.f for caveats) which is mvac.  Specifying vac_frac
!     the value of mvac is set to the value:
!                             mvac=ANINT(vac_frac_fg*mx)
!     For typical DIII-D, cases, a vac_frac_fg~0.25 will work

      REAL(r8) :: vac_frac=0._r8

!-----------------------------------------------------------------------
!     For gridshape="circ", one can pack around rational surfaces whose
!     value are specified by qpack(:).  The packing is Gaussian (i.e.,
!     a plot of dr vs. i shows a Gaussian) with the amplitude and width
!     of the Guassian by ampx(:) and wxpack(:).  The amplitude is roughl
!     dr0/dr where dr0 is the value far away from the Gaussian and wxpac
!     is the full-width given as a percentage of a = xmax-xmin.
!     For example, if mx=164 and you want to pack an extra 100 points
!     around the q=1 rational surface within a width of 0.1*a then:
!          qpack(1)=1.
!          wxpack(1)=0.1
!          ampx(1)=15.625
!        where ampx = dr0/dr = 1/(164-100) / (0.1*1/100) = 15.625
!     A value of zero for any quantity turns off the packing.
!
!     For gridshape="rect", one can pack around x or y positions whose
!       value are specified by xpack(:) or ypack(:).  The packing is
!       Gaussian (i.e., a plot of dx vs. i shows a Gaussian) with the
!       amplitude and width of the Guassian by ampx(:) and wxpack(:).
!       The amplitude is roughly dx0/dx where dx0 is the value far
!       away from the Gaussian and wxpack is the full-width given as
!       a percentage of a = xmax-xmin.
!       For example, if mx=164 and you want to pack an extra 100 points
!       around the q=1 rational surface within a width of 0.1*a then:
!       	xpack(1)=1.
!       	wxpack(1)=0.1
!       	ampx(1)=15.625
!        where amp = dr0/dr = 1/(164-100) / (0.1*1/100) = 15.625
!       A value of zero for any quantity turns off the packing.

      REAL(r8), DIMENSION(15) :: xpack=0._r8
      REAL(r8), DIMENSION(15) :: qpack=0._r8
      REAL(r8), DIMENSION(15) :: wxpack=0._r8
      REAL(r8), DIMENSION(15) :: ampx=0._r8
      REAL(r8), DIMENSION(15) :: ypack=0._r8
      REAL(r8), DIMENSION(15) :: wypack=0._r8
      REAL(r8), DIMENSION(15) :: ampy=0._r8

!     This also works for gridshape="rect" for packing in the
!       x-direction if firstx=0.  In these cases, qpack has nothing
!       to do with safety factor, its just the x-value for packing
!       locations.
!     The switch to get this form of packing has been revised (5/06) to
!       based on amp alone.

!-----------------------------------------------------------------------
!     The block decomposition of the region of rectangular cells is
!     controlled by nxbl and nybl [the number of rblocks in the x and
!     y (or radial and poloidal) directions, respectively].  The
!     distribution of blocks on parallel machines is controlled by
!     decompflag (0 = clumped, 1 = strided).

      INTEGER(i4) :: nxbl=1
      INTEGER(i4) :: nybl=1
      INTEGER(i4) :: decompflag=0

!-----------------------------------------------------------------------
!     For nonlinear runs without FFTW linked, set the number of cells
!     in the phi-direction through lphi.  The number of cells is 2**lphi,
!     and the number of toroidal modes after dealiasing is 2**lphi/3 + 1
!     (including n=0) using integer algebra.
!     For nonlinear runs with FFTW linked, set the number of cells
!     in the phi-direction directly through nphi. The number of toroidal
!     modes after dealiasing is nphi/dealiase + 1 (see below). Although
!     any nphi may be used, optimal values are the products of small
!     primes (e.g. 2^l * 3^n * 5^m * 7^p * 11 * 13).
!     For linear runs, the number of Fourier modes is set directly with
!     lin_nmodes.  One may specify the maximum toroidal mode number
!     for these cases to get the range
!          lin_nmax-(lin_nmodes-1) <= n <= lin_nmax,
!     otherwise the default gives
!          0 <= n <= lin_nmodes-1

      INTEGER(i4) :: lphi=2  ! used without FFTW linked
      INTEGER(i4) :: nphi=-1 ! used with FFTW linked
      INTEGER(i4) :: lin_nmodes=1
      INTEGER(i4) :: lin_nmax=0

!-----------------------------------------------------------------------
!     The following option allows nimrod to skip the dealiasing in
!       nonlinear computations (when dealiase is less than 3).  In this
!       case, there are 2**lphi/2 + 1 Fourier components. Otherwise,
!       there are nphi/dealiase + 1 Fourier components. The default
!       value, 3, dealiases for quadratic quantites. Values of
!       dealiase greater than 3 are only available if FFTW is linked.

      INTEGER(i4) :: dealiase=3

!-----------------------------------------------------------------------
!     For parallel computing, set the number of processor-layers into
!     which the Fourier components are distributed.

      INTEGER(i4) :: nlayers=1

!-----------------------------------------------------------------------
!     Parameters that control map_mod spline guess behavior. Good
!     defaults are chosen, however these can be tweaked as a last
!     resort if the spline guess mapping fails to initialize.
!     mm_mint/ext control the number of interior/exterior map_mod
!     points for the inverse map of a spline guess. The mm_imap_use_int
!     boolean determines whether nodes other than vertex nodes are
!     used to create the inverse map. It can be helpful to set it to
!     false for cases with large amounts of packing. mm_mdmap sets
!     the resolution of the direct map. mm_rw_dmap reads a dmap file
!     (if found), otherwise it writes a dmap file.

      INTEGER(i4) :: mm_mint=-1
      INTEGER(i4) :: mm_mext=-1
      LOGICAL :: mm_imap_use_int=.TRUE.
      INTEGER(i4) :: mm_mdmap=-1
      LOGICAL :: mm_rw_dmap=.FALSE.

!=======================================================================
!     physical constants.
!=======================================================================
!     This set of input parameters allow the user to change physical
!       constants at run-time.  Previously, these constants had be
!       fortran paramters that could only be changes by modifying
!       physdat and recompiling.  These input specifications allow
!       greater flexibility, including running normalized sets of
!       equations, but they also allow inconsistencies between pre-
!       processing and the calculation if the parameters are
!       inadvertently changed between the two steps.

!     The default values are MKS, except the Boltzmann constant is for
!       temperature in eV.  Note that the appearance of factors of c
!       are for MKS and not CGS.  If you really want CGS, the closest
!       possibility may be to consider the computed B as B/c
!       (Gauss-s/cm) and set mu0 to 4*pi/c**2.  Vacuum permittivity
!       is always computed from 1/(c**2*mu0).

!     Note that if the const_input namelist does not appear in the
!       nimrod.in file, the default values appearing in physdat.f
!       are used, not the input defaults below.  This is unique among
!       input parameters.

      REAL(r8) :: chrg_input=1.60217733e-19_r8 !  elementary charge
      REAL(r8) :: zeff_input=1._r8  !  relative ion charge
      REAL(r8) :: mi_input=3.3435860e-27_r8 !  ion mass
      REAL(r8) :: me_input=9.1093898e-31_r8 !  electron mass
      REAL(r8) :: gam_input=5._r8/3._r8  !  ratio of sp. heats
      REAL(r8) :: kblz_input=1.60217733e-19_r8 !  Boltzmann constant
      REAL(r8) :: mu0_input=3.1415926535897932385_r8*4.e-7_r8 !  vacuum
      REAL(r8) :: c_input=2.99792458e8_r8 !  speed of light
      REAL(r8) :: mc2_kT_input=1._r8 ! test particle mc**2/kT
      REAL(r8) :: mc2_kTb_input=1._r8 ! background particle mc**2/kT

!     The following is not an input parameter; do not put it in the
!     namelist.

      LOGICAL :: set_phys_constants=.false.

!=======================================================================
!     initialization parameters.
!=======================================================================
!     The initial perturbation is controlled by the following set.  bamp
!     is the magnitude of the perturbed B.  nx and ny set the number
!     of half-wavelengths in the x and y (radial and poloidal)
!     directions, respectively.  init_type selects what is being
!     initialized; though, you are not guaranteed an eigenfunction.
!     Note that "linear b" matches the DEBS perturbation and is only
!     suitable for linear geometry.
!
!     "sound" perturbations - one has the option of using theta_init and
!      phi_init to rotate the perturbation with respect to k_pll
!
!     For all wave-type perturbations, allow the ability to make it
!     smaller or larger than the box with linit_fac; e.g.,
!           cos(kx*(xmax-xmin)*linit_fac)
!
!     "eig" perturbations - for gridshape='flux' only.
!       These perturbations are either eigenfunctions coming directly
!       from an external file (currently the "chum" files which contain
!       GATO eigenfunctions) or by creating an eigenfunction (see below)
!
!     "stell" - field initialization for stellarator simulations.
!       For geom="lin" an analytic magnetic field is formed from helical
!         components defined by poloidal mode number, pert_m, toroidal
!         mode number, pert_n, with amplitude, pert_cos.
!       For geom="tor" the normal component of the magnetic field is
!         loaded along the domain boundary from the file bfield.txt
!         which contains the bi-Fourier decomposition.
!       "stell te" sets a uniform perturbed temperature across the
!         domain with magnitude odd_pert_str and allows temperature to
!         diffuse through the boundary when nimrod is run.
!       "stell nd" sets the perturbed density to have the same shape as
!         the temperature with size renormalized to ndens/odd_pert_str.
!
!     "coil" - use an external Biot-Savart routine to set magnetic field
!       If bs_type='oculus' then the routine bs00aa is called.
!         The number and specification of coils is set within oculus.in
!         (see input_oculus for coil options)
!       IF bs_type='bsc' then the bsc_mod routines are used.
!         The coils are read from the coils."cfext" file with standard
!         VMEC form.
!       The currents for each coil/group can be scaled separately by
!         setting the corresponding entries of cur_fac.
!       The full coil system can also be spatially shifted in the
!         periodic dimension by setting per_shift to be a fraction of
!         the periodic coordinate. (e.g. 0<per_shift<1)
!         *If more than nc_max coils/groups are needed then the
!            hard-coded value of this parameter must be increased or
!            seg faults will ensue.
!       "coil bnd" loads fields along the domain boundary only.
!       "coil fltr" and "coil bnd fltr" uses high toroidal resolution to
!         set fields but only loads a minimal number of modes which
!         reduces divergence error associated with coils
!       **this is not related to the applied external Bfield option
!         described below under Extra Parameters.
!
!
!     For init_type='reset', see below.

      REAL(r8) :: bamp=1.e-5
      INTEGER(i4) :: nx=1
      INTEGER(i4) :: ny=2
      REAL(r8) :: init_x, init_y
      REAL(r8) :: blob_amp, blob_width
      REAL(r8) :: theta_init=0._r8
      REAL(r8) :: phi_init=0.5_r8
      REAL(r8) :: linit_fac=1._r8
      CHARACTER(10) :: bs_type="bsc" ! "oculus"
      CHARACTER(100) :: cfext="in"
      INTEGER(i4), PARAMETER :: nc_max=20
      REAL(r8), DIMENSION(nc_max) :: cur_fac=1._r8
      REAL(r8) :: per_shift=0._r8

      CHARACTER(16) :: init_type="compr alf" ! "compr alf", "shear alf",
                                             ! "whistler", "linear b",
                                             ! "eig","blob","reset"
                                             ! "stell","coil"

!-----------------------------------------------------------------------
!     init_type="eig-v" or "eig-b"
!     Perturbations have a general form of:
!        S(r,T,P)= Sum_mn[ A_mn(r) COS(mT-nP) + B_mn(r) SIN(mT-nP) ]
!     where
!        S=psi [ B_poloidal = grad(phi) X grad(psi) ]
!        T=theta [ Straight-field-line poloidal angle ]
!        P=phi   [ Normal toroidal angle used in nimrod ]
!
!     The resonant components are chosen by pert_m and pert_n
!     The magnitudes are chose by pert_cos (A_mn above) and pert_sin
!     (B_mn above).
!
!     Perturbations that have A_mn(r=r_s)=0 are dubbed interchange
!     parity (odd parity) and those that have A_mn(r=r_s) != 0 are dubbe
!     tearing parity (even parity).  Controlling the shape of A_mn and
!     B_mn are controlled by the pert_type.
!
!     Resonant perturbations (where q_s=M/N) always have perturbations
!     of the form COS(MT-NP).  To enable perturbations of the opposite
!     helicity, COS(MT+NP), one just changes the sign of pert_n.  pert_m
!     should always be positive.
!
!     As a concrete example, the following NAMELIST:
!               pert_m   = 2   3   2
!               pert_n   = 1   1  -1
!               pert_cos = 1.  1.  1.
!               pert_sin = 1.  0   0
!     yields, ignoring the radial shape functions:
!             S(r,T,P) = cos(2T-P) + sin(2T-P) + cos(3T-P) + cos(2T+P)
!
!     These parameters are meant for setting perturbations in nimset.
!     They are also used as runtime parameters for setting perturbations
!     on the wall.  See section on Dirichlet conditions below.
!
!     The use of width is a parameter that controls the shape fo the
!     radial function depending upon pert_type.
!     pert_type='gaussian'
!                 radial_function=EXP(-(x/pert_width)**2)
!                 where x = distance from rational surface
!
!     A phase can also be added to the helical for each mode
!
      INTEGER(i4),   DIMENSION(20) :: pert_m= 99999999
      INTEGER(i4),   DIMENSION(20) :: pert_n=-99999999
      REAL(r8),      DIMENSION(20) :: pert_cos=0.
      REAL(r8),      DIMENSION(20) :: pert_sin=0.
      REAL(r8),      DIMENSION(20) :: pert_width=-99999999.
      REAL(r8),      DIMENSION(20) :: pert_phase=0.
      CHARACTER(10), DIMENSION(20) :: pert_type=''

!-----------------------------------------------------------------------
!     init_type="reset"
!     When this option is chosen, the initial perturbation comes from
!     another NIMROD dump file ("reset_file").  Things to keep in mind:
!     1) the equilibrium and mapping files are set in the same way.
!        This option affects ONLY the perturbation fields,
!     2) Only modes that are contained in the reset_file are
!        initialized, all others are set to zero.  If you are requesting
!        modes 0-5 and reset_file only contains mode 1, then modes 0,2-5
!        will be zero.
!     3) If the grid has not changed (that is if the RZ locations are
!        within a "reset_tol"), then the fields are just copied over.
!        This is very useful for changing the block decomposition or for
!        adding/deleting modes.  For these cases, make sure you change
!        dump_file appropriately.
!     4) Extreme care must be made when changing the grid size during a
!        nonlinear run since energy could be introduced in the system.
!        Changing the grid works best for linear scans.
!
!     Other options for init_type=reset:
!       reset_time and reset_step set the physical time and the
!         numerical time-step index in the new file.  Note, if either of
!         these values are negative, then the values contained in the
!         reset_file are used.
!
!     "rescale_factor" rescales the perturbations read in from a dump
!       file.   This is useful for initializing a linear scans where the
!       magnitudes can rapidly become very large.  If rescale_factor is
!       negative, then rescale_factor is set to
!                   bamp/MAXVAL(REAL(B_phi(lowest k)))
!
!     reset_ve_file is old velocity_grid_*.dat for resetting F's

      CHARACTER(128) :: reset_file="none"    ! For init_type="reset"
      REAL(r8) ::       reset_tol=1.e-7
      REAL(r8) ::       reset_time=0         ! Reset the time
      INTEGER(i4) ::    reset_step=0         ! Reset the step
      REAL(r8) ::       rescale_factor=1.
      CHARACTER(128) :: reset_ve_file="none" !path_to/velocity_grid_e.dat
      CHARACTER(128) :: reset_vi_file="none" !path_to/velocity_grid_i.dat
      CHARACTER(128) :: reset_vh_file="none" !path_to/velocity_grid_h.dat

!-----------------------------------------------------------------------
!    If a eq_dump_file is specified then that file is used to set the
!      equilibrium, instead of generation by standard means. The
!      number of modes in the eq_dump_file does not need to match that
!      of the file being created. However, the grid, poly_degree, and
!      block decomposition must match.  Note eq_nmodes must be set.

      CHARACTER(128) :: eq_dump_file='none'
      INTEGER(i4) :: eq_nmodes=1_i4

!-----------------------------------------------------------------------
!     The transfer_eq option allows one to transfer the equilibrium
!       field arrays (from an equilibrium file or an analytic model)
!       into the n=0 solution vector for use as an initial condition.
!       The transfer_n0 is the inverse of this operation. This is done
!       with nimset, and may be done during a dump file reset, but
!       nonlinear must be true.  Set transfer_eq_fields to transfer by
!       specifying the field as a string, up to 'bvnt'. Care should be
!       taken with density as continuity options are not preserved.
!       The transfer_overwrite option sets the eq or n0 to be overwritten
!       by the transfer, instead of incremented.  Finally, if
!       init_type_transfer_n0/='none', the perturbations will be
!       re-calculated with init_type set to init_type_transfer_n0.
!     See also the rt_transfer_eq option in physics parameters

      LOGICAL ::       transfer_eq=.false. ! logical for compatibility
      CHARACTER(4) ::  transfer_eq_fields='bvt' !'all','bvt'
      CHARACTER(4) ::  transfer_n0='none'       !'all','bvt'
      CHARACTER(16) :: init_type_transfer_n0='none' ! Like init_type
      LOGICAL :: transfer_overwrite=.false.

!=======================================================================
!     equilibrium parameters.
!=======================================================================

!-----------------------------------------------------------------------
!     Parameters for the 0-th order plasma:  The electron number
!     density, ndens, is in m**(-3).  In select cases, additional
!     profile flexibility is offered for density as a function
!     of some normlalized coordinate rho (=SQRT(psi_normal)).
!     For "n_profile"=
!      constant:     n(rho)=ndens
!      power:  n(rho)=(ndens-nedge)*(1-rho**npp(1))**npp(2)+nedge
!      poly:   n(rho)=ndens+SUM(npp(j)*rho**(j), [j,1,10])
!      tanh_pow: power method core blended with pedestal tanh
!      eqfile: Use n from the equilibrium data file if available.
!         Currently this is only available from some efit files.

      REAL(r8) :: ndens=1.e20                   ! density in core
      REAL(r8) :: nped=1.e20                    ! density in pedestal
      REAL(r8) :: nedge=1.e20                   ! density at edge
      CHARACTER(8) :: n_profile="constant"
      REAL(r8), DIMENSION(10) :: npp=0.         ! n profile parameter

!-----------------------------------------------------------------------
!     Parameters for the steady-state electron pressure: Like density,
!     the electron pressure is not normally considered as part of the
!     steady-state equations; therefore, one is often free to choose the
!     steady-state profile for the electron pressure.  Currently, only
!     two options are given for specifying the steady-state profile.
!
!     For "pe_profile"=
!      pe_frac: pe_0=pe_frac*p_0 (total steady-state pressure)
!         Note that this option is closely related to "separate_pe",
!         "pe_frac", and "tequil_rate" which control the dynamics of
!         the plasma evolution.  If separate pressures (now separate
!         temperatures) and temperature equilibration are used,
!         equilibrium electron and ion temperatures should be the same.
!         Thus, pe_frac should be zeff/(1+zeff).
!      eqfile: Use profiles coming from equilibrium file if available.
!      pfile: Use profiles from equilibrium pfile if available.
!         Note that for eqfile or pfile, separate_pe MUST be true, and
!         pe_frac has no real meaning since pe_0/p_0 can be a function
!         of radius now.
!      neTifile: Use profiles ne and Ti from equilibrium pfile if
!         available.  Electron pressure is then determined by
!         pe_0=p_0 - (ne from pfile/zeff)*(Ti from pfile)*kboltz*mu0
!         Note that n_profile must equal 'pfile' for this option

      CHARACTER(8) :: pe_profile="pe_frac"
      REAL(r8) :: te_cent=10000  ! in eV
      REAL(r8) :: te_edge=50     ! in eV
      REAL(r8), DIMENSION(10) :: tepp=0. ! T_e profile params

!-----------------------------------------------------------------------
!       For many of the lamprof profiles, it is specified using the
!        pres_2,pres_4 variables.  See those profiles before.
!        some of the profiles can be generalized:
!
!       For "p_profile"=
!         power: p(rho)=(p_0-p_edge)*(1-rho**pres_2)**pres_4+p_edge
!            where p_0   =ndens *te_cent* e * mu0
!            where p_edge=n_edge*te_edge* e * mu0
!         eqfile: Use profiles from equilibrium file if available.
!         pfile: Use profiles from equilibrium pfile if available.
!         power: p(rho) = (p_cent-pedge)*(1-rho**ppp(1))**ppp(2) + pedge
!                the option is only used in fgnimeq, and pe

!       Note that the power profile is only used in fg_nimeq. Here pedge
!        is calculated from the pressure gauge or the density and the
!        terperature gauges.

      CHARACTER(16) :: p_profile="eqfile  " ! power
      REAL(r8) :: pres_2=0
      REAL(r8) :: pres_4=0
      REAL(r8), DIMENSION(10) :: ppp=0.  ! P profile parameters
      REAL(r8) :: p_cent=1.0e5           ! in Pa

!     Parameters to add a pedestal to the pressure profile when using
!     fg_nimeq. This is added to the pressure profile in p_profile when

!     For "pped_profile"=
!       none: don't add perturbation
!       tanh: Add a tanh pedestal
!       p(psi) = pped_height (1-tanh ((psi-pped_psi0)/pped_width))

      CHARACTER(16) :: pped_profile="none  " ! tanh
      REAL(r8) :: pped_height = 0.0
      REAL(r8) :: pped_width = 0.0 ! in psi_norm
      REAL(r8) :: pped_psi0 = 0.90 !

!     For "pmod_profile"=
!       none: don't add perturbation
!       erf: Add error function to the pprofile

      CHARACTER(16) :: pmod_profile="none  " ! none
      REAL(r8) :: pmod_amp = 0.0
      REAL(r8) :: pmod_center = 0.0 ! in psi_norm
      REAL(r8) :: pmod_width = 0.0 ! in delta psi_norm
!-----------------------------------------------------------------------
!     Parameters for the steady-state f-profile. These inputs are inteded
!     for use with fg_nimeq. These inputs are set here for consistency
!     with the fg_nimeq profile inputs.

!     For "f_profile"=
!       eqfile: Use profiles coming from equilibrium file.

      CHARACTER(16) :: f_profile="eqfile" ! use eqfile

      REAL(r8) :: f_cent=3.247      ! in Tm
      REAL(r8) :: f_edge=3.243      ! in Tm
      REAL(r8), DIMENSION(10) :: fpp=0. ! f profile params

!     In addation we can add a a small perturbation to the above f_profile
!     The intent is to add a small localized current to the equilibrium
!     which helps drive a tearing mode unstable.

!     For "fmod_profile"=
!       none: don't add perturbation
!       erf: Add error function to the fprofile

      CHARACTER(16) :: fmod_profile="none  " ! none
      REAL(r8) :: fmod_amp = 0.0
      REAL(r8) :: fmod_center = 0.0 ! in psi_norm
      REAL(r8) :: fmod_width = 0.0 ! in delta psi_norm

!-----------------------------------------------------------------------
!     Parameters for the 0-th order magnetic field for cases not reading
!     an equilibrium file:  be0 is the magnitude (at the polar origin
!     for circular grids).  For rectangular grids (with lam0=0), the
!     field is uniform, and its direction is given in polar angles
!     thetab and phib.
!
!     For circular grids in linear geometry, these parameters specify
!     the J/B profile.  lam0 is the J/B ratio at the magnetic axis in
!     units of 1/(minor radius).  Either rbreak or alpha are used to
!     control the shape of the profile according to the choice of
!     lamprof.  For paramagnetic equilibria, just specify lam0 and
!     be0, and the parallel current will be determined as if E_axial
!     and resistivity were uniform.
!
!     The pressure (and perpendicular current) are specified with the
!     beta, pres_2, and pres_4 parameters described for the pitprs
!     lamprof option below.
!
!     These lambda profiles can be applied to rectangular
!     cases, too, if mx is even and the grid spacing is uniform in
!     x.  This produces a sheared slab model, which can be rotated
!     about the x axis by setting thetab/=0.
!
!     Profiles for rectangular slab cases:
!     lamprof="sheet", "zcos", "jchanx", and "jchany"
!     1) 'sheet' => From Mirnov et al. PoP 2004, Ahedo and Ramos NF 2009
!          2fl benchmark of Sovinec et al. JCP 2010
!          By~Tanh[x/LB] Bz~Guide field, epsilon_B=By(x->inf)/Bz,
!          lambda=epsilon_B/LB*Sech[x/LB]**2
!          With epsilon_B << 1 : Delta'=(2/k*LB**2)*(1-(k*LB)**2)
!          phib=LB, be0=Bz(0), lam0=By(x->inf)/be0
!          p~beta*be0**2/2mu0*(1 + pres_2*Tanh[x/LB]) (drift-tearing)
!          reflect_slab=True, reflects this profile
!     2) 'zcos' => From Fitzpatrick PoP p1782 2003 (Taylor reconn. probl
!          By=Bc*(Tanh[x/a]-(x/a)Sech[1]^2), be0=Bz~const., a=r_wall
!          beta=p(0)*2mu0/Bc**2, lam0=a*lambda(0)*be0*=Bc*(1-Sech[1]**2)
!          Pressure profile peaked at and symmetric about x=0
!          With lamprof='zcosp' p(0) is chosen such that p(a)=0

      REAL(r8) :: be0=1
      REAL(r8) :: thetab=0
      REAL(r8) :: phib=0.5

      CHARACTER(6) :: lamprof="zero"    ! "mbfm", "alpha", "para"
                                        !  also "pitprs", "qspec" see below
      REAL(r8) :: lam0=0
      REAL(r8) :: rbreak=0.2  ! break point (in a) for MBFM
      REAL(r8) :: alpha=3.    ! lam=lam0*( (a-r)/a )**alpha
      LOGICAL :: reflect_slab=.FALSE.

!-----------------------------------------------------------------------
!     For rectangular linear geometry, one can specify a "gauss" profile
!       initial condition which is useful for studying reconnection
!       problems similar to Biskamp.  The equilibrium is given by:
!       psi=pit_0*(exp(-r0**4/4) + exp(-r1**4/4)) where
!           r0=((x-x0)**2 +(y-y0)**2)**0.5
!       The magnetic field is given by B=Be0 zhat + zhat x grad(psi)
!
!     For linear geometry with a circular cross section, one may create
!       nonzero beta equilibria with the pitch & pressure profile or
!       with the q and pressure profile (more convenient for tokamaks).
!       To invoke these profiles, set lamprof="pitprs" or "qspec".
!
!       For lamprof="pitprs":
!       Pitch is defined as (r*Bz/a*Bth), i.e. q*R/a, and it's
!       specified through the input variables pit_0, pit_2, and pit_4:
!
!       pitch(r)=pit_0*(1+pit_2*x**2+pit_4*x**4) ,
!
!       where x=(r-xmin)/(xmax-xmin)
!       The pressure is specified through beta, pres_2, pres_4 as
!
!       p(r)=beta*be0**2/(2*mu0)*(1+pres_2*x**2+pres_4*x**4) .
!
!       For lamprof="qspec":
!               See Holmes et.al., Phys. Fluid 26 (1983) 2569
!       The q profile is specified through the input
!       variables pit_0, pit_2, and pit_4:
!
!       q(r)=pit_0*(1+(x/pit_2)**(2*pit_4))**(1/pit_4) .
!
!       The pressure is specified through beta as
!
!       p(r)=beta*be0**2/(2*mu0)*Integrate[ x/q**3 * (2q - xq') ,(x,1,r)
!
!       Note that xmin must be zero for the qspec option.

      REAL(r8) :: pit_0=1
      REAL(r8) :: pit_2=0
      REAL(r8) :: pit_4=0

!-----------------------------------------------------------------------
!     Select an equilibrium flow profile.  This is treated like all
!     other _eq quantities--bicubic and fixed in time.  The choices
!     are:
!     1) 'none'
!     2) 'uniform'=>This produces a uniform flow with magnitude
!           ve0, specified in MKS units (m/s).
!           The direction is set by thetav and phiv
!           -- see nimset/nimset_init.f.
!     3) 'gauss'=>This gives a Gaussian flow distribution centered
!           at the axis of polar grids (or the horizontal center
!           of rectangular grids) with width eqflow_width, relative
!           to the root of the flux function (approximately).
!           The peak of the profile has magnitude ve0, and the
!           direction with different gridshapes is the same as for
!           the uniform option.
!     4) 'pinch'=>An axial or toroidal electric field is computed
!           from the resistivity and equilibrium current density at
!           the magnetic axis.  ve_eq is then solved from the
!           perpendicular part of E_tor=ve_eqXbe_eq - eta*j_eq.
!           Caution:  E_tor is considered to be ~ 1/R, whether or
!           not this is appropriate for a given equilibrium.
!     5) 'eq_file'=>The equilibrium flow is taken from the
!           equilibrium solver file where possible.
!     6) 'diamagnetic'=>Compute the ion diamagnetic flow from the
!           equilibrium B and J, assuming J_perp is from grad(P)
!           and using P_i=(1-pe_frac)*P.
!
!     Note that one must set the advect variable to get equilibrium
!     flow in the momentum equation, regardless of the eq_flow
!     selection.

      CHARACTER(16) :: eq_flow='none'
      REAL(r8) :: ve0=0
      REAL(r8) :: thetav=0
      REAL(r8) :: phiv=0.5
      REAL(r8) :: eqflow_width=0.5
      REAL(r8) :: dfac=1

!-----------------------------------------------------------------------
!     Set the total plasma pressure based on the equilibrium magnetic
!     field (pres_eq and be_eq, respectively) at the magnetic axis for
!     equilibria not initialized from an equilibrium file.  In ALL cases
!     the perturbed pressure is advanced, and momentum is affected by
!     pressure gradients only when beta>0.

      REAL(r8) :: beta=0

!-----------------------------------------------------------------------
!     For slab equilibria:
!     Set the equilibrium length for the graviational mode equilibrium
!       The first vector component of the gravity parameter, declared
!       with physics input, also affects the profile.

      REAL(r8) :: glength=0. ! gravitational scale length in meters
!     For the ITG/ETG studies, these are useful:

      REAL(r8) :: telength=0 ! T_e scale length in meters
      REAL(r8) :: tilength=0 ! T_i scale length in meters

!     Normalized ion temperature profile can be given by
!         T_i0(x) = 1 + bb_tanh * TANH(x / tilength)
!      Requires bb_tanh<1
!      DS: 01/23/2012

      REAL(r8) :: bb_tanh=0.9

!     eta_i = L_n/L_Ti
!     eta_e = L_n/L_Te
!     DS: 3/30/2011

      REAL(r8) :: eta_i = 0.
      REAL(r8) :: eta_e = 0.

!-----------------------------------------------------------------------
!     Parameters for the neutral density:
!     Options shown are implemented as will corresponding plasma
!     variables

      CHARACTER(8) :: neut_dens_profile="constant"
      REAL(r8) :: neutral_dens=1.e20    ! m**(-3)
      CHARACTER(8) :: eq_neut_flow="none" ! or constant, gauss
      REAL(r8) :: vn0 = 0.0             ! m/s
      CHARACTER(8) :: neut_temp_profile="constant"
      REAL(r8) :: neutral_temp=1.0      !  eV

!-----------------------------------------------------------------------
!     Shaping parameters for both resistivity and viscosity,
!     each is the constant specified above times a function.
!     ds_function="orig":
!           (1+(SQRT(dvac)-1)*(r/rv)**dexp)**2
!     ds_function="tanh":
!           0.5*(dvac -1)*tanh((r/rv-1)*dexp) + 0.5*(dvac +1)
!     where r is the radial logical coordinate in the collection
!     of polar rblocks and rv = radius of vacuum (xvac).
!     In core triangles the factor is 1, and in blocks outside the
!     separatrix it is dvac.
!
!     The character variable, ds_use, controls which diffusivities are
!     shaped by this function.
!
!     For ds_function='spitzer', the function used is that when
!     eta_model='full' or 'n=0 only' only the function remains fixed:
!             D=MIN( elecd*(T/eta_ref_t)**(-1.5), elecd_max ).
!     This is useful for comparing linear and nonlinear tokamak cases.
!     See comments for eta_model and eta_ref_t.
!
!     For the orig and tanh functions, we also have the "-gq" options
!     which make the functions into nimrod runtime parameters by
!     modifying the values of the diffusivity shape function at the
!     gaussian quadrature points to have the following form:
!              WHERE(ds>0.5*(1+dvac))
!                      ds=dvac
!              ELSEWHERE
!                       ds=1
!              END WHERE
!     This creats a sharp function exactly at the transition region
!     while avoid numerical problems that those functions can provide.
!     This can be especially useful for matching with ideal MHD codes.
!     This option only affects the resistivity
!
!     If ds_nqty is >1 then the diff_shape nqty will be increased
!     such that the shaping varies by field-solve:
!     nq=1 for B, nq=2 for v, nq=3 for n, nq=4 for T
!     The velocity diff_shape is only used with ds_use="both" or
!     "kin_visc".

      REAL(r8), DIMENSION(4) :: dvac=1
      REAL(r8), DIMENSION(4) :: dexp=20
      CHARACTER(8) :: ds_use="both" ! or "elecd","kin_visc","neither"
      CHARACTER(16), DIMENSION(4) :: ds_function="orig" ! or "tanh","spitzer"
      INTEGER(i4) :: ds_nqty=1

!-----------------------------------------------------------------------
!     For circular gridshapes, the following allows one to put a
!     vacuum region from xvac < r < xmax.  The dimension of the
!     current channel is then scaled to xvac instead of xmax.  The
!     default value gives no vacuum region.

      REAL(r8) :: xvac=0

!-----------------------------------------------------------------------
!     The following parameter is meant as a last ditch resort for
!     diagnosing possible problems with the equilibrium.  It has been
!     successfully used with the "CDXU benchmark problem", but should
!     not be used routinely.  When activated, it will calculate the
!     parallel current using finite element calculations in
!     nimrod_init.f

      LOGICAL :: tor_eqja_fe=.false.

!-----------------------------------------------------------------------
!     The following options are used by nimset and are most useful for
!     experimental equilibria. If calc_eq_fe is set to either b or j
!     the poloidal magnetic field is computed from psi, and the current
!     is computed from b with a finite element solve. The parameter
!     calc_jeq_fe_nq determines the number of components of j to be
!     updated in this solve (2 for RZ, and 3 for RZPhi). If
!     calc_jeq_fe_gs is true, the Grad-Shafranov equation is used to
!     determine the phi component of current. This equation involves
!     single derivatives of psi and F, whereas the calculation from B
!     implies a second derivate of psi. However, it also requires
!     division by grad(psi)**2. Near x-points this may be unsuitable
!     as there is a saddle point and grad(psi)->0. calc_jeq_fe_bnd
!     may be used to specify a bound on grad(psi)**2 where if it is
!     smaller than the bound, the RHS is set by curl(Bpol).
!
!     Additionally, a diffusion equation may be advanced for the vector
!     fields by setting any combination of 'bjv' in diffuse_eq_var.
!     The _divfd, _diff, and _phif arrays hold a divergence cleaner
!     diffusion coefficient, a diffusivity (similar to elecd for B), and
!     a coefficient for the phi-component diffusivity for their
!     respective fields. The time-step size for this advance is set by
!     dtm, and the number of steps for each field by the _steps array,
!     respectively for each field. For the last _dvbsteps steps, the
!     diffusivity is set to zero and only divergence cleaning is
!     performed. diffuse_eq_report sets the step frequency of reporting
!     from diagnostics during the advances. The diffusivity uses the
!     diffusivity profile set by ds_function, unless specified
!     differently through diffuse_eq_{dsftn,dexp,dvac}.

      CHARACTER(2) :: calc_eq_fe='no' ! 'bj'
      INTEGER(i4) ::  calc_jeq_fe_nq=3
      LOGICAL ::      calc_jeq_fe_gs=.true.
      REAL(r8) ::     calc_jeq_fe_bnd=0     ! Try 1.e-3 to 1.e-6

      CHARACTER(4) :: diffuse_eq_var='none' !'bjvs'
      REAL(r8), DIMENSION(4) :: diffuse_eq_divfd   =(/0.,0.,0.,0./)
      REAL(r8), DIMENSION(4) :: diffuse_eq_diff    =(/0.,0.,0.,0./)
      REAL(r8), DIMENSION(4) :: diffuse_eq_phif    =(/1.,1.,1.,1./)
      INTEGER(i4), DIMENSION(4) :: diffuse_eq_steps   =(/5,5,5,5/)
      INTEGER(i4), DIMENSION(4) :: diffuse_eq_dvbsteps=(/0,0,0,0/)
      CHARACTER(10) :: diffuse_eq_dsftn  = 'same'
      REAL(r8) ::      diffuse_eq_dexp   = 1._r8
      REAL(r8) ::      diffuse_eq_dvac   = 1._r8
      INTEGER(i4) ::   diffuse_eq_report = 1

!-----------------------------------------------------------------------
!     Parameters used to enfore ideal-like tokamak conditions.
!     Comparison with ideal-MHD codes is difficult, often more difficult
!     than more realistic cases. For comparison, the vacuum region
!     should have low density and high resistivity with a discontinuous
!     jump to the plasma values across the separatix. Setting ideal_like
!     to true will adjust the following fluxgrid input parameters:
!       psihigh -> 0.9999 - Element boundary on the separatrix
!       ideal_like_fg -> .TRUE. - Use the plasma values on the
!                                 separatrix nodes
!     and the following NIMROD input parameters:
!       n_profile -> 'bump' - Discontinuous density profile
!       npp(1) -> 0.99999
!       ds_function -> 'psep' - Discontinuous diffusivity profile
!                               set vacuum diffusivity to dvac
!       scal_byparts -> .FALSE. - Use div(nv) directly in the
!                                 continuity equations.
!       p_computation -> 'at nodes'
!     Finally, ideal_like=.true. will force the linear nimrod advance
!     to use a pressure formulation instead of temperature. This advance
!     only supports p_model='isothermal', 'adiabatic' and 'isotropic'
!     where k_perpi is used as the diffusivity coefficient with unis
!     of m^2/s. It also ignores all heating and sources.
!
!     If set_vac_quadpts is true, the quadrature point data for density,
!     diff_shape and temperature is post-adjusted in the vacuum
!     region. This requires mvac_nim and pressure_gauge_nim to be set.
!     These latter two parameters correspond to their fluxgrid
!     equivalents

      LOGICAL :: ideal_like=.false.
      LOGICAL :: set_vac_quadpts=.false.
      INTEGER(i4) :: mvac_nim=-1_i4         ! Set if set_vac_quadpts
      REAL(r8) :: pressure_gauge_nim=-1._r8 ! Set if set_vac_quadpts

!-----------------------------------------------------------------------
!     Equilibrium loop voltage (not used in the fluid advance).
!     If using an EFIT k-file, reverse the sign.

      REAL(r8) :: eq_loop_voltage=0  ! Need dump_eexp=T

!-----------------------------------------------------------------------
!     Determine what angle is used for setting and analyzing
!     perturbations in toroidal geometry.  The choices are to use a
!     simple geometric theta, or with a straight-field-line angle.
!     A straight field line angle is more accurate for determining
!     resonant perturbations, but becomes degenerate near the separatrix
!     which can cause problems.  In these cases, using the geometric
!     angle works better.
!     The sfa_ratsurf option is a hybrid where the theta value is
!     defined as a constant on the grid, but the values are chosen
!     as the straight-field line angle at the last rational surface of
!     a perturbation given by pert_m(1) and pert_n(1)

      CHARACTER(12) :: tor_angle="sfa"  ! or "geom", "sfa_ratsurf"

!-----------------------------------------------------------------------
!     Plasma state parameters
!     For SWIM runs, these are useful to control the plasma state output
!     gfilename=name of eqdsk file based on NIMROD dump file
!     cur_state_file=name of plasma_state file to write out
!     The plasma state labels are defined like this:
!       ps%Tokamak_ID = trim(pstok)
!       ps%Shot_Number = istep   ! Reuse
!       ps%Global_label = trim(run_id)//'_'//trim(pstok)

      CHARACTER(15) :: gfilename='g000001.0000'
      CHARACTER(25) :: cur_state_file='ps.cdf'
      CHARACTER(25) :: pstok="diiid-like"
      CHARACTER(128) :: run_id='swim'

!=======================================================================
!     physics parameters.
!=======================================================================
!     If nonlinear is false, the equilibrium field is held fixed in
!     time, and electric fields and Lorentz forces use only linear
!     terms.

      LOGICAL :: nonlinear=.true.

!-----------------------------------------------------------------------
!     The rt_transfer_eq option allows one to transfer the equilibrium
!     field arrays (from an equilibrium file or an analytic model)
!     into the n=0 solution vector at RunTime (rt).  This is
!     functionally equivalent to the initial condition functionality in
!     nimset, but the dump files continue to have the equilibrium fields
!     separated from the n=0 fields so that one can change this
!     parameter on the fly.  It also allows the various *_source options
!     (see closures namelist) to work in "transfer eq mode".

      LOGICAL :: rt_transfer_eq=.false.

!-----------------------------------------------------------------------
!     Select the terms to appear in the generalized Ohm's law.  Valid
!       choices are: 1) 'mhd', 2) 'hall', 3) 'mhd&hall', or 4) '2fl'
!       Options 2 - 4 now include the electron pressure term
!       term when beta>0, and 4 adds electron inertia.

      CHARACTER(8) :: ohms="mhd"

!-----------------------------------------------------------------------
!     Potentially temporary variable to allow for Hall physics without
!       contributions of Grad(Pe) to the magnetic field advance.

      LOGICAL :: no_pe_hall=.false.

!-----------------------------------------------------------------------
!     The eta_model input allows the user to select whether to use
!       a temperature (electron) dependent resistivity, proportional
!       to T**-3/2.  Valid input options are:
!
!       Note that for linear computations, the "eta full" option will
!       include the linear resistivity perturbation according to the
!       linear electron temperature perturbation, while the
!       "eta n=0 only" option will not.  Both will use an unperturbed
!       resistivity profile based on the equilibrium temperature,
!       however.
!
!       1) "fixed"		skip the T-dependence
!       2) "eta n=0 only"	use only the symmetric part of T
!       3) "eta full"         use 3D eta(T)
!       4) "eta full ds"	same as 3) but multiply by diff_shape
!       5) "chodura"		uses the Chodura model--see below
!
!       For options 2 & 3, the electrical diffusivity has the form
!
!       D=MAX( MIN( elecd*(T/eta_ref_t)**(-1.5), elecd_max ), elecd_min)
!
!       The elecd coefficient is the electrical diffusivity in m**2/s
!       (MKS resistivity/mu0).
!
!       For eta_ref_t < 0, the central temperature is used
!           (corresponding to 0,0 on the logical grid)
!
!       When eta_model="chodura" the electrical diffusivity is computed
!       from the sum of the "eta full" T**(-1.5) model and the
!       phenomenological Chodura model [Nucl. Fusion 15, 55 (1975)]:
!       {actual ref. is Milroy and Brackbill Phys. of Fl. 25, 775 (1982)}
!
!       elecd_chodura*SQRT(ndens/n)*(1-EXP(-f_chodura*v_drift/v_ia))
!
!       Where v_drift is the electron drift speed, |J/ne|, and v_ia
!       is the ion-acoustic speed, SQRT(gamma*P/mtot*n), each computed
!       locally in 3D.  At present the drift speed computation is not
!       centered in the time-advance of B, so this contribution is
!       first-order accurate in dt.

      CHARACTER(16) :: eta_model="fixed"
      REAL(r8) :: elecd=0
      REAL(r8) :: elecd_max=1000
      REAL(r8) :: elecd_min=0
      REAL(r8) :: eta_ref_t=1
      REAL(r8) :: f_chodura=1
      REAL(r8) :: elecd_chodura=1

!-----------------------------------------------------------------------
!     Select the advection terms to appear in the the momentum equation
!     and the generalized Ohm's law.  Valid choices are: 1) 'none',
!     2) 'V only', or 3) 'all'.

      CHARACTER(8) :: advect="V only"

!-----------------------------------------------------------------------
!     Select how the continuity equation for number density is applied.
!     Valid choices are:
!     1) 'none' density perturbations are not evolved.
!     2) 'fix profile' evolve density but use only the equilibrium
!          density in the velocity advance,
!     3) 'n=0 only' evolve density but use only the equilibrium plus
!          n=0 part of density in the velocity advance.
!     4) 'full' evolve density and use the full 3D density profile in
!          the velocity advance.  This is the only self-consistent
!          option for nonlinear runs, but it's also the most
!          computationally intensive, since the velocity advance
!          requires a 3D matrix solve at each step.

      CHARACTER(16) :: continuity="none"

!-----------------------------------------------------------------------
!     Isotropic diffusivity used in the number density evolution
!       (continuity) equation that crudely represents particle transport
!       for long time-scale simulations.  Units are m**2/s.
!     The nd_hypd is a hyper-diffusivity for number density evolution
!       and is only applied to computations with implicit advection,
!       i.e. ohms/=mhd or mhdadv_alg=centered.  Units are m**4/s.

      REAL(r8) :: nd_diff=0
      REAL(r8) :: nd_hypd=0

!     When either is nonzero, it is mathematically sound to choose
!       (homogeneous) Dirichlet conditions.
!     You should set ndensity_bc="dirichlet" see below

!-----------------------------------------------------------------------
!     Viscous dissipation in the velocity equation:
!       different viscous stress tensors may now be added together.
!
!       kinematic:     Pi(V) = -rho*kin_visc*grad(V)
!       isotropic:     Pi(V) = -rho*iso_visc*W
!       parallel:      Pi(V) = -rho*3*par_visc/2*(b.W.b)*(bb-I/3)
!       gyroviscosity: Pi(V) = -eta3/2*gyr_visc*
!                               ((bxW).(3bb+I)-(3bb+I).(Wxb))
!
!       where b is the magnetic direction vector, eta0 is the ion
!       pressure times the ion collision time, and W is the rate of
!       strain tensor,
!
!         W=grad(V)+grad(V)^T-(2/3)*div(V)
!
!       and kin_visc, iso_visc, and par_visc have units of m**2/s.
!       the diffusivity shaping specified below is used in both the
!       kinematic and the isotropic coefficients.

      REAL(r8) :: kin_visc=0
      REAL(r8) :: iso_visc=0
      REAL(r8) :: par_visc=0
      REAL(r8) :: gyr_visc=0

!-----------------------------------------------------------------------
!     Filter for the Fourier direction.  The filter is applied during
!     the velocity advance as additiona viscosity
!
!       kinv_visc k^2 + hv_coeff (k/kmax)**(2 hv_power)
!
!     hv_coeff will be set to zero for lin geomery and when
!     total_nmodes=1
!
!     ap***** Old filter implmentation:
!     input varibales hv_exp, hv_time, hv_filter are obsolete
!     *****ap
!     At the end of each time
!     step, a "filter" is applied:
!                 V=V A(k)
!     where
!     A(k) = exp(.5*ln(1-hv_exp)*(k/kmax)**4 (dt/hv_time))
!     if hv_filter=1
!     and
!     A(k)=1 - hv_coeff * (k/kmax)**(2 hv_power)*(dt/hv_time)**2
!     otherwise.
!
!     In choosing values for these parameters, hv_coeff is the amount to
!     decrease the maximum Fourier number at each time step, and
!     hv_power controls how rapidly the hyperviscosity increases with
!     toroidal mode number.

      INTEGER(i4) :: hv_filter=0
      REAL(r8)    :: hv_coeff=0.
      REAL(r8)    :: hv_time=1.e-7
      INTEGER(i4) :: hv_power=10

!-----------------------------------------------------------------------
!     If gravity.ne.0, add gravitational force to simulate the effect of
!     the magnetic lines' curvature

      REAL(r8), DIMENSION(3) :: gravity = 0.

!-----------------------------------------------------------------------
!     Select the independent calculation of the electron pressure
!     Warning: separate_pe=.TRUE. with ohms/="2fl" and p_model/="aniso2"
!      introduces an additional (unphysical) heating source in the model

      LOGICAL :: separate_pe=.FALSE.

!-----------------------------------------------------------------------
!     Set the electron pressure as a fraction of total.  For cases
!     where the electron pressure is advanced separately
!     (separate_pe=T) AND the pressure equilibrium is not read in from
!     an equilibrium file, pe_frac sets the equilibrium and initial
!     electron pressures.
!     When the equilibrium is read in from a file, the equilibrium
!     electron pressure is read in if available, and pe_frac is set to
!     be the values at the axis (pe_frac=peo/po), otherwise, the input
!     value of pe_frac is used to set it.
!     For cases without a separate advance, electron pressure is pe_frac
!     times the total pressure for all time.  Note: set 0<=pe_frac<=1.
!
!     Note that with separate pressures (now separate temperatures)
!     and temperature equilibration, equilibrium electron and ion
!     temperatures should be the same.  Thus, pe_frac should be
!     zeff/(1+zeff).

      REAL(r8) :: pe_frac=0.5

!-----------------------------------------------------------------------
!     Apply a toroidal electric field to the n=0 'perturbation.'
!     This is used to self-consistently generate the 0-th order
!     part of the field for nonlinear runs when an equilibrium file
!     is not available.  Use be0 to set the vacuum toroidal field at
!     the center of the polar mesh and loop_volt to apply a loop
!     voltage (in Volts).  The 'equilibrium' arrays will remain 0,
!     except for the vacuum toroidal field, and the n=0 'perturbed'
!     part will serve as the equilibrium field - the vacuum toroidal B
!     + the n=0 perturbation.  tloopv0 & tloopv1 are the
!     beginning and ending times to ramp the applied voltage from 0
!     to the full value of loop_volt.
!     CAUTION:  If the the equilibrium field has current, loop_volt
!     will effectively add to it.  To avoid equilibrium current when
!     gridshape='circ', ensure that lam0=0.

      REAL(r8) :: loop_volt=0
      REAL(r8) :: tloopv0=-1
      REAL(r8) :: tloopv1=0
      LOGICAL  :: loop_volt_file=.FALSE.
      REAL(r8), DIMENSION(20000) :: loop_volt_time
      REAL(r8), DIMENSION(20000) :: loop_volt_volt
      INTEGER :: loop_volt_nrow

!-----------------------------------------------------------------------
!     The following allows the voltage to vary in order to achieve a
!     desired amount of current.  The voltage is adjusted according
!     to V = V_old + dt*loop_rate*(i_desired-i_n0)
!            - dt*loop_rate2*d(i_n0)/dt      ,
!     where i_n0 is just the perturbed part.  If the equilibrium
!     fields have current, i_desired is interpretted as the difference
!     between the actual desired current and the equilibrium field
!     current.  This is used instead of the voltage ramp whenever
!     i_desired is not 0.  Note that the di/dt part of the feedback
!     is not applied during the first time-step (leading to minor
!     discrepancies if a restart is inserted into a sequence of
!     time-steps), and V_old at the beginning of a simulation is
!     manually set with loop_volt.
!
!     Be aware that the behavior of the current
!     depends on the simulated plasma!

      REAL(r8) :: i_desired=0       !   in Amperes
      REAL(r8) :: loop_rate=1.e-4   !   in Ohms/s
      REAL(r8) :: loop_rate2=0      !   in Ohms

!-----------------------------------------------------------------------
!     The next set applies a vertical electric field, possibly for
!     electrostatic current drive in compact devices.  The
!     t_e_vert0-&-1 are beginning and ending times for ramping the
!     electric field.
!     t_e_vert2 is the ramp down time uses slope of (t_e_vert1-t_e_vert0)

      REAL(r8) :: e_vertical=0      !   in Volts/m
      REAL(r8) :: t_e_vert0=-1
      REAL(r8) :: t_e_vert1=0
      REAL(r8) :: t_e_vert2=1.e12

!-----------------------------------------------------------------------
!     adding gun_max for flux injection

      REAL(r8) :: gun_max=0.
      REAL(r8) :: gun_lim0=0.
      REAL(r8) :: gun_lim1=-1.
      CHARACTER(8) :: gun_func='ramp'

!-----------------------------------------------------------------------
!     Normally, homogeneous Dirichlet boundary conditions are applied to
!     the change in the normal component of the magnetic field.
!     To eliminate possible round-off errors, one can explicitly set the
!     normal magnetic field (perturbation, not equilibrium) to zero afte
!     returning from the matrix solve.  This is the normal procedure:
!          magnetic_bc="homogeneous"
!       If set to magnetic_bc="no enforce homo", homogeneous
!     dirichlet boundary conditions will naturally applied to delta-B.
!
!     Non-Homogeneous Dirichlet boundary conditions can be applied in
!     nimrod for select cases.
!     Taylor problem:
!      The Taylor problem (given by lamprof='cosh' in slab geometry)
!      is to study forced reconnection when a conducting wall is
!      perturbed according to:
!      x_w = F(t) G(y,z)
!        where F(t)=exp(-t/mag_tau) * ( 1- exp(-t/mag_tau) )
!              G(y,z)=SUM[ A_mn*COS(my-nz) + B_mn*SIN(my-nz) ]
!      The time dependence is chosen to make the perturbation and
!      derivative continuous at time t=0.
!      See discussion of perturbations in initialization parameters for
!      a complete description of the spatial dependence.
!
!     Resistive wall characterized by a wall velocity vwall, which
!     allows for direct comparison to the alfven velocity.  The "ef"
!     parameters determine the optional error field surface current
!     efkphi in the phi direction, and the radial location of the
!     current, efrad.
!     reswall: Use matrix-vector formalism for calculating the boundary
!              condition.
!        -mmax sets number of poloidal fourier modes included
!              in cyl_wall_matrix
!        -reswall_nmax: set's max torodial mode number for resistive wall
!              n>nmax will be treated as a conducting wall
!     reswall_rmp: applies an error field outside of resistive wall.
!              see rmp inputs.
!     cylreswall: Use analytic boundary conditions for resistive wall in
!            cylinder

      CHARACTER(16) :: magnetic_bc="homogeneous" !"no enforce homo"
      REAL(r8) :: mag_tau=1.
      REAL(r8) :: vwall=0.    !vwall=eta_wall/(mu0*delta_wall)
      REAL(r8) :: efkphi=0.
      REAL(r8) :: efrad=1.
      INTEGER(i4) :: mmax=6
      INTEGER(i4) :: reswall_nmax=100
      REAL(r8) :: almcon=1.   ! almcon=fourier transform normilization
                              ! for per_cyl vac response mat_create
!-----------------------------------------------------------------------
!     Variables for applying RMP boundary conditions.
!
!     RMP_type = "linear" (formerly IZZO) linearly ramps the RMP using
!       B(t) = BRMP (t )/(t_RMP) if 0<t<t_RMP
!     RMP_type = pulse (formerly BEID) applies aa pulse using
!       B(t) = BRMP (1.0 -(1.0 + t/t_rmp) exp(-t/T_rmp) )
!       The rmp has a max of BRMP(1+e**-2) at T_RMP, and then asymptotes
!       to BRMP at infty
!     RMP_type = "expdecay" ramps the profile using an exp "decay"
!       Bn(t) = BRMP (1 - exp( -(t-toffset)/T_RMP) )
!     RMP_type = bump uses a bump function pulse
!       x = ( (t-t_rmp_offset) / t_rmp) **2
!       Bn(t) = BRMP ( exp ( -1/(1-x) ) ) for x < 1
!       Bn(t) = 0 for x >= 1
!
!     The value of BRMP is saved as a lagr_edge_type.
!
!     For nonlinear cases the flag load_n0_rmp determines how the n=0 rmp
!     fields are treated
!
!     If rmp_mult_fac/=0 then multiply the RMP fields by this factor
!
!     rmp_omega is the angular frequence of the rmp
!
!     All rmp bcs now use lagr_edge_rmp, this is set to true if
!     magnetic_bc = rmp or reswall_rmp

      REAL(r8) :: t_rmp=0
      REAL(r8) :: t_rmp_offset =-1.0e-8
      REAL(r8) :: rmp_mult_fac = 0.0_r8
      CHARACTER(12) :: rmp_type="IZZO"
      CHARACTER(12) :: rmp_file_type="allnode" !surfnode
      LOGICAL :: b_nonsym=.false.
      LOGICAL :: reverse_ip=.false.
      LOGICAL :: load_n0_rmp =.false.
      LOGICAL :: lagr_edge_rmp = .FALSE.
      REAL(r8) :: rmp_omega = 0.0_r8
!-----------------------------------------------------------------------
!     We typically use no-slip boundary conditions, apart from a
!       possible specified surface normal from ExB.  The following
!       allows you to choose free-slip conditions independent of
!       viscosity.  However, no-slip should be used if equilibrium or
!       perturbed magnetic field penetrates the wall.

      CHARACTER(16) :: velocity_bc="no-slip"  !   or "free-slip"

!-----------------------------------------------------------------------
!     Density boundary conditions should at a minimum be consistent with
!     density diffusion.

      CHARACTER(16) :: ndensity_bc="none" ! or "dirichlet"

!-----------------------------------------------------------------------
!     Temperature boundary conditions.  Note that "none" implies the
!     default which is:
!       Enforce dirichlet when diffusion is used.
!       Do not enforce dirichlet when diffusion is not used.
!     Using the option temp_*_bc="no enforce homo" allows one to
!     not enforce dirichlet even when diffusion is used.  This is
!     equivalent to nimuw's insulate option.

      CHARACTER(16) :: temp_ion_bc="none" ! or "no enforce homo"
      CHARACTER(16) :: temp_ele_bc="none"

!-----------------------------------------------------------------------
!   n=0 parameters: if any of these are >0 then they will overwrite the
!     corresponding parameters (same names without the trailing zero)
!     for layer 0 only (requires nlayers=nmodes)

      REAL(r8) :: elecd0=-1
      REAL(r8) :: nd_diff0=-1
      REAL(r8) :: kin_visc0=-1
      REAL(r8) :: iso_visc0=-1

!-----------------------------------------------------------------------
!     Select if the ponderomotive force is added to momentum eq or not
!     NOTE: ponderomotive force comes from VSim simulations

      LOGICAL :: use_fp=.FALSE.

!=======================================================================
!     closure physics parameters.
!=======================================================================

!-----------------------------------------------------------------------
!     When performing transport studies, it is useful to control the
!     equations used; e.g., hold the magnetic field steady and evolve
!     only the temperature.  eqn_model is for these types of modeling.
!     eqn_model(6:9)="Tcor" determines if one does nonlinear correction
!     step for T advance.

      CHARACTER(12) :: eqn_model="all"  ! "T", "T,V", "T,V,n"

!-----------------------------------------------------------------------
!     THERMAL CONDUCTION CLOSURE MODEL
!
!     The variable qpi_model determines which closure model is
!       used to calculate thermal transport coefficients.
!       The available options are:
!       1) "std kprp n0"  !  simplified model type chosen by p_model
!       2) "standard"     !  p_model options with 3D perp conductivity
!       3) "braginskii"   !  Braginskii's magnetizing conductivity in 3D
!       4) "k2"           !  Ji's k2 model
!-----------------------------------------------------------------------
!     STANDARD MODEL
!     The variable p_model controls the level of complexity in these
!       simplified models, varying from no conduction up to the
!       parametric dependence of the Braginskii model in the
!       high-magnetization limit.  The available options are:
!
!     1) 'adiabat' uses advective and compressive terms but no
!           thermal conduction.
!     2) 'isothermal' excludes the compressive term in the
!           temperature equation(s).
!     3) 'isotropic' adds isotropic thermal conduction to the
!           temperature equation(s) with k_perpi coefficient
!           (and k_perpe for separate_pe=T)
!           This is like the former 'iso' option but renamed to
!           help distinguish the new isothermal option.
!     4) 'aniso1' uses anisotropic thermal conduction in the
!           temperature equation(s) with diffusivities k_perpi
!           and k_plle for single T and k_perpe and k_plli for
!           separate Ti and Te.
!     5) 'aniso2' This is the same as aniso1 except:
!           separate_pe MUST be true, the center-of-mass velocity is
!           used in the electron temperature advance instead of V_e.
!           The goal of this option is to allow for separate T_e and T_i
!           equilibrium profiles without running into the V_e CFL limit.
!           This is also more general than allowing for a pe_frac
!           profile.
!     6) 'aniso_plltdep' nonlinear-only option for running
!           temperature-dependent parallel thermal conductivity,
!           k_plls-> k_plls*(Ts/k_pll_ref_t)**2.5
!           where s = e or i (species).
!     7) 'aniso_tdep' nonlinear-only option for running
!           temperature-dependent parallel thermal conductivity
!           (see 6.) and temperature- and mod(B)-dependent
!           perpendicular thermal conductivity with
!           k_perps-> k_perps*(k_pll_ref_t/Ts)**0.5*(kprp_ref_b/B)**2
!           where s = e or i (species), B is in Tesla, and the fields
!           used in Ts and B are set by the qpi_model (eq+n=0 or 3D)
!     8) 'aniso_ntdep' option for running
!           temperature-dependent parallel thermal conductivity
!           (see 6.) and temperature-, density- and mod(B)-dependent
!           perpendicular thermal conductivity with
!           k_perps-> k_perps*(k_pll_ref_t/Ts)**0.5*
!                       (n/k_ref_n)**2*(kprp_ref_b/B)**2
!           where s = e or i (species), B is in Tesla, and the fields
!           used in Ts and B are set by the qpi_model (eq+n=0 or 3D)
!
!     With dimensions of m**2/s, all of the k_perps and k_plls would be
!       better labeled as chis. But for all anisotropic conduction
!       options the input diffusivities in are immediately converted to
!       conductivities by multiplying by ndens.
!
!     Note that there are also some subtleties with respect to the
!       influence of number density.
!     All anisotropic p_model options have the correct form of the
!       temperature equation
!         n*dT/dt=-(gamma-1)*[n*T*div(v)+div(n*chi.grad(T))+...]
!       where d/dt includes v.grad() and chi is the temperature
!       diffusivity.  Conductive heat flux is proportional to n
!       and spatial variations in n are included according to the
!       continuity input parameter.
!     For p_model='isotropic' the equation is simplified to
!         dT/dt=-(gamma-1)*[T*div(v)+div(chi*grad(T))]
!       ignoring the possible effects of density gradients.
!
!     Additionally, if anisotropic conduction is used, we assume
!       B_eq.grad(T_eq)=0 and J_eq.grad(n_eq)=0 (separate_pe=T),
!       which are stronger assumptions on the equilibrium fields
!       than the usual steady-state assumption.  If this is not
!       satisfied, consider using transfer_eq to make the equilibrium
!       fields part of the initial condition only.
!-----------------------------------------------------------------------
!     BRAGINSKII MODEL
!     Reference:  Braginskii. "Transport Processes in a Plasma."
!       Reviews of Plasma Physics. 1965.
!     Thermal diffusivity coefficients have the form:
!       chi_pll,s  = Ts*tau_s/ms * (gamma0_s/delta0_s)
!       chi_perp,s = Ts*tau_s/ms *
!           (gamma1_s*xs**2+gamma0_s)/(xs**4+delta1_s*xs**2+delta0_s).
!       where tau_i=SQRT(2)*tau_ii and tau_e=tau_ei
!
!     Unless artificially scaling or normalizing thermal diffusivity
!       coefficients, use k_plli=k_perpi=k_plle=k_perpe=k_pll_ref_t=1.
!
!     The default coefficients given below are for a Z = 1 plasma.
!       Ion coefficients for different Z plasma are given in the table
!       on page 251 of the reference.
!-----------------------------------------------------------------------
!     K2 MODEL
!     Reference:  Ji. "Closure and transport theory for
!       high-collisionality electron-ion plasma." 2011.
!     Thermal diffusivity coefficients have the form:
!       chi_pll,s  = T_s*tau_s/ms * f1(xs,zeta)
!       chi_perp,s = T_s*tau_s/ms * f2(xs,zeta)
!       where tau_s=tau_ss and zeta = sqrt(me/mi*ti/te)/Z
!     Unless artificially scaling thermal diffusivity coefficients,
!       use k_plli=k_perpi=k_plle=k_perpe=k_pll_ref_t=1.

      CHARACTER(12) :: qpi_model="std kprp n0"

!-----------------------------------------------------------------------
!     STANDARD CLOSURE MODEL:
!     The complexity of the model is determined by p_model as described
!       above.  Reference values for T, B, and n should correspond to
!       the values of k_plls and k_perps chosen for aniso_dep models.

      CHARACTER(16) :: p_model="adiabat"
      REAL(r8) :: k_prp_ref_b=-1._r8  ! unused if negative
      REAL(r8) :: k_ref_n=-1._r8      ! unused if negative

!-----------------------------------------------------------------------
!     BRAGINSKII CLOSURE MODEL:
!     The default coefficients given below are for a Z = 1 plasma.
!       Ion coefficients for different Z plasma are given in the table
!       on page 251 of the reference.

      REAL(r8) :: gamma1_ele=4.6640_r8
      REAL(r8) :: gamma0_ele=11.920_r8
      REAL(r8) :: delta1_ele=14.790_r8
      REAL(r8) :: delta0_ele=3.7703_r8

      REAL(r8) :: gamma1_ion=2.0000_r8
      REAL(r8) :: gamma0_ion=2.6450_r8
      REAL(r8) :: delta1_ion=2.7000_r8
      REAL(r8) :: delta0_ion=0.6770_r8

!-----------------------------------------------------------------------
!     BRAGINSKII AND K2 CLOSURE MODELS:
!     The magnetization coefficients used to calculate
!       chi_perp have the form:  xs=omega_cs*tau_s
!     Unless artifically scaling or normalizing magnetization
!       coefficients, leave magfac_ele=magfac_ion=1.
!
!     If tdep_coul_log=T, a temperature dependent coulomb logarithm
!       is computed for thermal diffusivity coefficients and
!       the thermal equilibration rate.  Otherwise, a constant value
!       of coulomb_logarithm is used for the calculations.
!
!     When using tdep_coul_log=T, leave coulomb_logarithm=1.

      REAL(r8) :: magfac_ele=1._r8
      REAL(r8) :: magfac_ion=1._r8

      LOGICAL :: tdep_coul_log=.false.
      REAL(r8) :: coulomb_logarithm=1._r8

!-----------------------------------------------------------------------
!     THE FOLLOWING CLOSURE PARAMETERS CAN APPLY TO ANY MODEL
!
!     k_perpe, k_plle, k_perpi, k_plli and k_pll_ref_t should be set
!       according to the qpi_model as described above.
!
!     k_perp replaces k_perpe and/or k_perpi if either equal 0.
!     k_pll replaces k_plle and/or k_plli if either equal 0.
!
!     k_pll_min and k_pll_max set limits on k_plle and k_plli in all
!       closure models where spatial variation is allowed.
!
!     Note that unlike in nimuw, the gamma-1 factors should not be
!       included in any of the input conduction coefficients, for any
!       model. This factor is properly accounted for in the code.

      REAL(r8) :: k_perp=0
      REAL(r8) :: k_perpe=0
      REAL(r8) :: k_perpi=0
      REAL(r8) :: k_pll=0
      REAL(r8) :: k_plle=0
      REAL(r8) :: k_plli=0
      REAL(r8) :: k_pll_ref_t=1._r8
      REAL(r8) :: k_pll_max=1.e14
      REAL(r8) :: k_pll_min=0

!-----------------------------------------------------------------------
!     for the nonlinear runs using k_pll and/or k_cross:
!       bhat(1:4)=='full' -> normalize bhat by full B
!       bhat(1:4)=='neq0' -> normalize bhat by n=0 + eq B

      CHARACTER(4) :: bhat="neq0"

!-----------------------------------------------------------------------
!     k_cross is a switch for bXgrad(T) Braginskii heat fluxes.
!       Unless artificially scaling or normalizing thermal diffusivity,
!        use k_cross=1. This term can only be used if separate_pe=.true.

      REAL(r8) :: k_cross=0

!-----------------------------------------------------------------------
!     To provide extra thermal diffusivity near boundaries
!       that represent a gap, the k_pll_ref_t/<Ts> ratio
!       has a lower bound of 1+(kprp_mnrat*(diff_shape-1)/dvac)**2
!       where diff_shape is the evalated diffusivity shape
!       profile determined by dvac and dexp.

      REAL(r8) :: kprp_mnrat=0

!-----------------------------------------------------------------------
!     When separate_pe=.true., tequil_rate sets heat flux density
!       from electrons to ions (and vice versa) as
!         Q_i=-Q_e=n_e*k*(Te-Ti)*tequil_rate/(gamma-1)
!       and n_s*k*dT_s/dt equation has the additional term, Q_s.
!     Equilibration is not used if tequil_rate=0, the default.

        REAL(r8) :: tequil_rate=0._r8

!-----------------------------------------------------------------------
!     Specify whether to use a temperature-dependent thermal
!       equilibration rate.  When tdep_tequil=T, use tequil_rate=1.
!       Coefficients are calculated within the thermal_equil routine
!       within field_comps.f.  This option only applies to nonlinear
!	computations; linear computations use the fixed tequil_rate.

      LOGICAL :: tdep_tequil=.false.

!-----------------------------------------------------------------------
!     closure_model controls the type of parallel closures that are
!       used.  The first three characters (1:3) are for electron
!       heat flow, the next three (4:6) are for ion heat flow, the next
!       three (7:9) are for electron stress and the final three (10:12)
!       are for ion stress.  The default is for standard Braginskii:
!         closure_model = 'stpstpstpstp'
!       stp stands for "standard temperature and pressure tensor"
!     To change any of the parallel closures to the mixed finite-element
!       (MFE) method, replace 'stp' with 'mfe'.
!     To change any of the parallel closures to the continuum solution
!       of the CEL-DKE, replace 'stp' with 'cel'.

      CHARACTER(12) :: closure_model="stpstpstpstp"

!-----------------------------------------------------------------------
!     kplle_mfe, kplli_mfe, and par_visc_mfe are scaling coefficients
!       for the MFE method which aid in convergence of the iterative
!       solver when the auxiliary equations are in play.
!     setting k_pll(e,i)_mfe = ndens*(k_pll(e,i)-k_prp(e,i))**.25 works
!       well for the MFE temperature advance.

      REAL(r8) :: kpll_mfe=0
      REAL(r8) :: kplle_mfe=0
      REAL(r8) :: kplli_mfe=0
      REAL(r8) :: par_visc_mfe=0

!-----------------------------------------------------------------------
!     divq_mfe_ibp controls the weak form of the mixed finite element
!     method divergence of q_parallel term. False means no integration
!     by parts is done, while true means the standard weak form is used.

      LOGICAL :: divq_mfe_ibp=.false.

!-----------------------------------------------------------------------
!     TFtheta sets how much kinetic and fluid closures get used.
!     TFtheta=0 is purely kinetic, while TFtheta=1 is purely fluid.

      REAL(r8) :: TFtheta = 0.d0

!-----------------------------------------------------------------------
!     closure_dump controls whether or not the parallel heat flows
!     and stresses are written out to the dump file.

      LOGICAL :: closure_dump=.false.

!-----------------------------------------------------------------------
!   n=0 parameters: if any of these are >0 then they will overwrite the
!     corresponding parameters (same names without the trailing zero)
!     for layer 0 only (requires nlayers=nmodes)

      REAL(r8) :: k_perpi0=-1
      REAL(r8) :: k_perpe0=-1
      REAL(r8) :: k_perp0=-1

!-----------------------------------------------------------------------
!     SOURCES
!
!     Understanding the sources is complicated by the separation of
!     variables into an "equilibrium component", and a "perturbed
!     component".  The upshot is that when running with separated
!     variables (most common mode of operation for the tokamak
!     simulations), NIMROD will use a "diffusive source" term which is a
!     source term of the form: D grad^2(Q_eq) for quantity Q and diffusi
!     coefficient D.
!
!     To give the full range of steady-state behavior, we have the
!     "diffusive source" and "zero ss" options.  See the presentation
!     by S. Kruger on this subject.  Note that the rt_transfer_eq
!     option is used
!
!     The self-similar heat source increases the n=0 component
!     linearly with time with gamma_heat specifying the rate of increase
!
!     For heat_ion_source option descriptions see comment block below.
!
!     For electric_source='RF', a localized source is set. See
!     comment block below.
!
!     All sources can be set to 'bump-eq' which relaxes the field
!     back to the EQ fields on the nu_relax time-scale. A bump
!     function mask is applied to make this a localized source:
!     mask = EXP(-1/(1-r/bump_a))/EXP(-1) where r=sqrt((R-xo)**2+(Z-yo)**2)
!     is less than bump_a and zero elsewhere.
!     For each field then the following time-centered term is added:
!     df/dt = ... - nu_relax * \tilde{f}_n0

      CHARACTER(20) :: momentum_source="none"
      CHARACTER(20) :: ndensity_source="none"
      CHARACTER(20) :: electric_source="none"    ! 'RF'
      CHARACTER(20) :: heat_ion_source="none"
      CHARACTER(20) :: heat_ele_source="none"
      REAL(r8) :: gamma_heat=0.
      REAL(r8) :: nu_relax=1.e4 ! in 1/s
      REAL(r8) :: bump_a=1. ! in m

!-----------------------------------------------------------------------
!     These are generic time dependence parameters used by several
!     source options

      REAL(r8) :: toffset=5._r8
      REAL(r8) :: tperiod=50._r8
      REAL(r8) :: tdown=1000000000._r8
      REAL(r8) :: tmodperiod=50._r8
      REAL(r8) :: tmodphase=0._r8
      REAL(r8) :: t_on=10000000._r8
      REAL(r8) :: t_off=10000000._r8

!-----------------------------------------------------------------------
!     Ion Heat Sources:
!     Only a few heating source options are described here.  For the
!     full list, users should review sources.F90.
!
!       'uniform':  heat is applied uniformly throughout the domain with
!           magnitude q_cent (in SI units: eV/m^3/s)
!       'n0':  heat is applied to the n=0 mode only with a bump
!           function profile in the poloidal plane.  the profile has
!           is centered on the point (q_x,q_y) in real space with radial
!           extent psi_lcfs and magnitude q_cent
!       'point':  heat is applied locally around (q_x,q_y,phi) in real
!           space with radial extent psi_lcfs and mag. q_cent, where
!           phi is the periodic coordinate corresponding to the
!           fractional coordinate odd_pert string (0<odd_pert_sting<1)
!       'temperature':  heat is applied following the existing
!           temperature profile, which can approximate flux surface
!           shapes if anisotropic condition is used.  a bump function
!           profile with magnitude q_cent is used so that heat is only
!           applied in regions that have T>psi_lcfs(%)*Tmax

      LOGICAL :: odd_pert=.false.
      REAL(r8) :: odd_pert_str=0.0
      REAL(r8) :: psi_lcfs=0.0
      REAL(r8) :: q_cent=0.0
      REAL(r8) :: q_x=0.0
      REAL(r8) :: q_y=0.0

!-----------------------------------------------------------------------
!     v_cent and psi_per are parameters for controlling stellarator
!     momentum sources
      REAL(r8) :: v_cent=0
      REAL(r8) :: psi_per=0.0

!-----------------------------------------------------------------------
!     RF SOURCES:
!     In general we will get RF data from the GENRAY code.  This can be
!     done directly (electric_source = 'rf_genray', more_phiplanes =
!     'false'), in which case the toroidal resolution needs to be very
!     good (perhaps prohibitively so, as the nphi toroidal planes
!     corresponding to the discrete toroidal Fourier transform need
!     to exist in sufficient number to resolve the details of the RF
!     deposition).
!     For this case, setting more_phiplanes = 'true' will generally
!     not make sense.
!
!     RF deposition calculations can also use GENRAY data for the R-Z
!     planes, and then average this data over phi and spread it
!     out over some characteristic toroidal width set by delta_tor.
!     For this case, one uses electric_source = 'rf_genray_aht' (ad hoc
!     toroidal) and more_phiplanes = 'true', with mynphi a power of 2
!     significantly exceeding nphi (factors of 8 or 16 to capture the
!     details of the deposition are probably necessary).
!     For this case, setting more_phiplanes = 'false' will generally
!     result in unphysical data, so the code will stop.
!
!     One can also put in the data in a fully ad hoc manner (with
!     electric_source = 'rf'), wherein the parameters lam0rf, rrf,
!     zrf, phirf, delta_pol, delta_tor parameterize the current.
!
!     See sources.F90 for a more detailed description.
!-----------------------------------------------------------------------
      REAL(r8) :: lam0rf=0._r8
      REAL(r8) :: rrf=0._r8
      REAL(r8) :: zrf=0._r8
      REAL(r8) :: delta_pol=1._r8
      REAL(r8) :: delta_tor=1._r8
      REAL(r8) :: phirf=3.1415926535897932385_r8

!-----------------------------------------------------------------------
!     GENRAY RF SOURCE:
!     Several parameters need to be set as nimrod input in addtion to
!     those set in genray input files.
!     The wave frequency (Hz), total RF source power (W), and parameters
!     for extrapolations must be set.
!
!     It is often useful to use a larger number of phi planes for
!     bi_interp_mod which can optionally be set.
!
!     For a more detailed description of genray useage, see sources.F90.

      REAL(r8) :: rf_power_ec=0._r8
      REAL(r8) :: rf_freq_ec=1._r8
      INTEGER(i4) :: discrete_extrap=4
      LOGICAL :: more_phiplanes=.false.
      INTEGER(i4) :: mynphi=8
      INTEGER(i4) :: rf_extrap=10
      CHARACTER(20) :: shepard_alg="cshep"

!-----------------------------------------------------------------------
!     If rf_modulation='on', multiply the RF by the factor
!     sin(2 pi t/tmodperiod + tmodphase), otherwise, no modulation.

      CHARACTER(20) :: control="none"
      INTEGER(i4), DIMENSION(2) :: pcsstatus=0
      CHARACTER(20) :: rf_modulation="off"

!-----------------------------------------------------------------------
!     The following parameters determine whether dissipated energy
!     is used as a source in the temperature equations.  If
!     separate_pe=F, they would both add to the single-fluid
!     temperature evolution.  If separate_pe=T, Ohmic heating
!     (eta*J**2) is used in the electron temperature evolution only,
!     and viscous heating rho*kin_visc*grad(v)^T.grad(v) is used in
!     the ion temperature evolution only.

      LOGICAL :: ohm_heat=.false.
      LOGICAL :: visc_heat=.false.

!-----------------------------------------------------------------------
!     Computation of the perpendicular viscosity coefficient is
!       controlled by perpvisc_model.  Its default behavior uses
!       perp_visc as a fixed viscous diffusivity. When set
!       to "prpdep" it uses the n=0 dependence part of the p_model
!       choices aniso_tdep and aniso_ntdep. "prpdep" is only valid when
!       p_model is one of these choices and the k_perp factor is
!       replaced by kin_visc or iso_visc.
!
!     This options overrides the use of the diffusivity profile through
!       ds_use='kin_visc' or 'both'

      CHARACTER(16) :: perpvisc_model="fixed" ! "prpdep"

!-----------------------------------------------------------------------
!     Computation of the parallel viscosity coefficient is controlled
!       by parvisc_model.  Its default behavior uses par_visc as a
!       fixed viscous diffusivity.  When it is set to "plltdep," the
!       Braginskii T**5/2 temperature dependence is used.  Since the
!       coefficient for parallel thermal diffusivity has exactly the
!       same dependence, the physics kernel computes the temperature-
!       dependent thermal diffusivity as described above for aniso_*
!       of p_model (using k_pll_ref_t, k_pll_max, and k_pll_min) and
!       multiplies the result by par_visc/k_plli to get the viscous
!       diffusivity.  This plltdep option for viscosity may be used
!       with or without the temperature-dependent thermal diffusivity.

      CHARACTER(16) :: parvisc_model="fixed" ! or "plltdep"
!-----------------------------------------------------------------------
!     HEURISTIC CLOSURE MODELS

!     neo_flag controls which heuristic closures are used for the
!     ion and electron stress tensors. Currently the options are
!     none and gianakon. The gianakon model corresponds to equation
!     12 in Gianakon, Kruger and Hegna POP 9, 2002.

!     The inputs mu_e and mu_i are the electron and ion viscous damping
!     frequencies.

!     These models require that <B_eq^2> is calculated in a preprocessing
!     step and stored in the dump file. dump_fsa_beq2 must be set to true.

      CHARACTER(8) :: neoe_flag="none" !gianakon
      CHARACTER(8) :: neoi_flag="none"
      REAL(r8) :: mu_e=0._r8    ! electon viscous damping frequency
      REAL(r8) :: mu_i=0._r8    ! ion viscous damping frequency

!     The heuristic closures for the neoclassical stress scale as
!     V_p,J_p \cdot B_p B_p/(B_p4) and are singular at the magenetic
!     axis. A bump function of the form
!                   psi = 1-amp*exp(-1/ (1-(r/r0)^2) ) for r/r0< 1
!                       = 0 r>r0
!     is used to zero out the neoclassic stresses near the axis
!     here r is the distance from the magnetic axis
!     calcualed using neo_axis_r,z

      REAL(r8) :: neo_bump_r0 = 1 ! in units of length
      REAL(r8) :: neo_bump_amp = 1.0 !bump function amp (set to e)
      REAL(r8) :: neo_axis_r = 0.5 !r location of axis
      REAL(r8) :: neo_axis_z = 0.0 !z location of axis
!-----------------------------------------------------------------------
!     KINETIC CLOSURE MODELS
!     kin_model controls kinetic equations that are solved.
!     kin_model(1:3)='d_f' or 'cel' for delta-f or CEL electrons
!     kin_model(4:6)='d_f' or 'cel' for delta-f or CEL ions
!     kin_model(7:9)='d_f'          for delta-f energetic particles
!     kin_model(10:12)='sim'        simultaneous delta-f e/i advance
!     kin_model(13:17)='onlyF' only advance F (w or wo limited fluid eqs)

      CHARACTER(17) :: kin_model="none"
!-----------------------------------------------------------------------
!     For continuum solution to CEL-DKE, nF = number coefficients
!     in the velocity space representation of the distribution
!     function, F = sum_l F_li(x,t,s_i) phi_l(vp/v), where
!     F_li's are known on a grid of ns speed grid points
!     used for the quadrature in the speed variable.
!     1=electrons, 2=ions and 3=hot particles
!-----------------------------------------------------------------------
!     s_quad# = 'gau_lags' generates speed grid points for weight
!     s**2*s**sp1*exp(-s**2) on s=[smin(#),smid(#)] and
!     s**2*s**sp2*exp(-s**2) on s=[smid(#),smax(#)] and
!     if smid(#)/=infinity then domain is split with ns-ns2 grid points in
!     s=[smin(#),smid(#)] and ns2 points in  s=[smid(#),smax(#)] for gau_lags.
!     s#_quad = 'legendre' uses speed grid points consistent with
!     Legendre quadrature on domain s=[smin,smax].
!     s#_quad = 'hot_part' uses speed grid points consistent with
!     slowing down distribution which has a speed weighting that goes
!     as s**4/(1 + s**4) on domain s=[smin,smax].
!
!     capability for adding assigned nodes in s grid.
!     can be used to put a node at s=0 or some other desired place
!     like a transition point to the tail of a hot particle dist.
!     note: each additional point reduces exactness of s-quadrature
!     scheme by one order.
!     ns#_fixed: # of assigned nodes is s grid in domain #; 0,1,2 are possible
!
!     whenever nF>0, the dump file contains the coefficients, F_ln.

      INTEGER(i4), DIMENSION(3) :: nF =(/0,0,0/) ! e,i,h
      INTEGER(i4), DIMENSION(3) :: ns =(/0,0,0/) ! e,i,h
      INTEGER(i4), DIMENSION(3) :: ns2=(/0,0,0/) ! e,i,h
      INTEGER(i4), DIMENSION(3) :: ns1_fixed=(/0,0,0/)       ! 0, 1 or 2
      REAL(r8), DIMENSION(6) ::     s1_fixed=(/0,0,0,0,0,0/) ! ee,ii,hh
      INTEGER(i4), DIMENSION(3) :: ns2_fixed=(/0,0,0/)       ! 0, 1 or 2
      REAL(r8), DIMENSION(6) ::     s2_fixed=(/0,0,0,0,0,0/) ! ee,ii,hh
      REAL(r8), DIMENSION(3) :: smin=(/0,0,0/) ! e,i,h
      REAL(r8), DIMENSION(3) :: smid=(/0,0,0/) ! e,i,h
      REAL(r8), DIMENSION(3) :: smax=(/0,0,0/) ! e,i,h
      CHARACTER(8), DIMENSION(3) :: s1_quad =                           &
     &                  (/'gau_lags','gau_lags','hot_part'/)
      CHARACTER(8), DIMENSION(3) :: s2_quad =                           &
     &                  (/'gau_lags','gau_lags','hot_tail'/)
      REAL(r8), DIMENSION(3) :: sp1=(/0,0,0/)
      REAL(r8), DIMENSION(3) :: sp2=(/0,0,0/)

!     controls if s collocation approach uses Lagrange (F) or associated
!     s polynomials (T) for interpolation and derivative computations.
      LOGICAL :: s_poly_exp = .true. ! .false. = Lagrange

!     old_nF, old_ns, old_ns2, and old_ph_pd are for reset dumpfile when
!     remapping F and are put into the nimrod.in used for the reset

      INTEGER(i4), DIMENSION(3) :: old_nF =(/0,0,0/) ! e,i,h
      INTEGER(i4), DIMENSION(3) :: old_ns =(/0,0,0/) ! e,i,h
      INTEGER(i4), DIMENSION(3) :: old_ns2=(/0,0,0/) ! e,i,h
      INTEGER(i4) :: old_ph_pd=1

!-----------------------------------------------------------------------
!     for shaping high-energy tail in pitch angle a la Belova

      REAL(r8) :: dxi = 0._r8
      REAL(r8) :: xi0 = 0._r8

!-----------------------------------------------------------------------
!     for shaping high-energy tail in speed; parameters below are ~
!     for matching orbit-RF's calculation of beam ions driven by RF.

      REAL(r8) :: tpow=3.5_r8
      REAL(r8) :: tcoe=1._r8

!-----------------------------------------------------------------------
!     gridshape_v(1:7) = 'circuni', 3 domains in pitch-angle with
!     m_xip nodes in:                      -1 < xi < xi_node_bnd
!     m_xit nodes in:            -xi_node_bnd < xi < xi_node_bnd
!     m_xi-m_xip-m_xit nodes in:  xi_node_bnd < xi < 1
!
!     gridshape_v(1:7) = 'circtpb'
!     m_xip nodes in:                      -1 < xi < xi_trapped
!     m_xit nodes in:            -xi_trapped  < xi < xi_trapped
!     m_xi-m_xip-m_xit nodes in:  xi_trapped  < xi < 1
!     minimum cell size for trapped space -tpb_min<=xi<=tpb_min
!
!     gridshape_v(5:9)='tpbap' xi_trapped(R,Z) at outboard midplane
!     gridshape_v(5:9)='tpbex' xi_trapped(R,Z) at each spatial location
!     'tpbex' requires more computation (especially for the collision
!     operator) but is generally worth it in terms of convergence in xi

      CHARACTER(9) :: gridshape_v = 'none' ! 'circuni', 'circtpb'
      INTEGER(i4) :: m_xi=3
      INTEGER(i4) :: m_xip=1
      INTEGER(i4) :: m_xit=1
      REAL(r8) :: xi_node_bnd=0.5
      REAL(r8) :: tpb_min=0.5

!-----------------------------------------------------------------------
!     pd_xi is the polynomial degree of the 1D FE in pitch angle, xi,
!     Total number of degrees of freedom in vp/v domain = m_xi*pd_xi+1
!     and nF = ns * (m_xi*pd_xi+1).
!     pa_basis controls basis type
!     Lagrange and Gauss-Lobatto nodal bases as well as a hierarchecal
!     model basis built from Legendre polynomials are possible

      INTEGER(i4) :: pd_xi=0
      CHARACTER(3) :: pa_basis='Lag' !'Leg','Lob'

!-----------------------------------------------------------------------
!     nq_xi is the number of quadrature points in pitch-angle
!     integration. nq_xi = 2*(pd_xi+1) seems to work well for most of
!     the terms in the drift kinetic equation.
!     nq_xi = pd_xi+1 makes qp's for GLL polynomials the same as nodes.
!     C_nq_xi is the number of quadrature points used in the projection
!     of pitch-angle FE's onto Legendre moments in the calculation of
!     the collision operator.  C_nq_xi=200 underdoes it in extremely
!     collisional regimes.  We're talking nu_star like 1e2 to 1e4 where
!     4000 or 8000 would be better.  C_nq_xi is not used when Legendre
!     basis functions are used in pitch-angle.

      INTEGER(i4) :: nq_xi=0
      INTEGER(i4) :: C_nq_xi=800

!-----------------------------------------------------------------------
!     fraction of hot particle pressure in slowing-down distribution.

      REAL(r8) :: betafrac_sdd=1.0_r8

!-----------------------------------------------------------------------
!     close_nprocs controls number of processors devoted to the
!     kinetic closure calculation;  must be at least half of processors
!     if kinetic closures are calculated.

      INTEGER(i4) :: close_nprocs=0

!-----------------------------------------------------------------------
!     analytic_rhs should be set to 'vector' in closure_input if one
!     needs the jfac and jmat structures to project an analytic rhs
!     onto the NIMROD basis; see adv_aux.f for further documentation

      CHARACTER(6) :: analytic_rhs="none"

!-----------------------------------------------------------------------
!     single_s -> logical variable for determining if DKE is solved
!                 separately for each s.
!     steady_state = .true. removes partial F/ partial t from CEL-DKE
!     cel_output -> choose to write out matrix elements and such.
!     fcel -> centering parameter for split F advance for free-streaming
!             and particle trapping terms
!     fcol -> centering parameter for collision terms.
!     fder -> centering parameter for s derivative terms and d / dxi
!             acceleration term involving b dot grad ln n
!     fint -> centering parameter for implicit treatment of moments of
!             f_NM in CEL-DKEs
!     fadv -> centering parameter for implicit treatment of Maxwellian
!             advective terms, such as v dot grad lnT f_M
!     coll_model(1:7)=moments ! uses moment expansion form
!     coll_model(1:7)=s_coloc ! uses consistent F expansion
!     coll_model(1:7)=lorentz ! pitch-angle scattering only
!     coll_model(1:7)=neither

      LOGICAL :: single_s = .false.
      LOGICAL :: steady_state = .false.
      LOGICAL :: cel_output = .false.
      REAL(r8) :: fcel=1.
      REAL(r8) :: fcol=0.
      REAL(r8) :: fder=0.
      REAL(r8) :: fint=0.
      REAL(r8) :: fadv=0.
      CHARACTER(7) :: coll_model = 'neither'
      REAL(r8) :: spitz_fac  ! corrects Ramos form of coll. friction

!-----------------------------------------------------------------------
!     it_calls -> number of calls to iterative solver
!     xi_prec = .true. does preconditioning in pitch-angle
!     ss_prec = .true. does preconditioning at single speed
!     ss_rhs = .true. compute rhs and dot routines at single speed
!     ssoff_prec = .true. combined with ss_prec = .true. does
!     preconditioning at single speed, but includes a small set of
!     off-diagonal components

      INTEGER(i4) :: it_calls=1
      LOGICAL :: xi_prec = .false.
      LOGICAL :: ss_prec = .false.
      LOGICAL :: ss_rhs = .false.
      LOGICAL :: ssoff_prec = .false. !Only has an effect for ss_prec=T

!-----------------------------------------------------------------------
!     cel_flow and cel_temp control whether flow and temperature are
!     from Chapman-Enskog-like approach with zero V and T moments in F
!     F_init_type is used to set up initial distribution functions for
!     kinetic calculations.
!     dump_read_F is part of reset capability for adding distribution
!     fuctions to dump files that don't have them.
!     Deff_s includes an isotropic spatial diffusion term in F eqs.
!     Dh_s includes an isotropic spatial hyper-diffusion in F eqs.
!     both provide spatial smoothing but one must be careful that they
!     do not affect the physical results.
!     Shaped versions are controlled using F_dvac, F_dexp to construct
!     the F_diff_shape 2D lagrange structure

      LOGICAL :: cel_flow=.false.
      LOGICAL :: cel_temp=.false.
      CHARACTER(16) :: F_init_type="homogeneous" ! 'phmix', 'spitz_the
      LOGICAL :: dump_read_F(3)=(/.false.,.false.,.false./)
      LOGICAL :: F_mass_op_flag=.false.
      LOGICAL :: kinetic_eq=.false.
      REAL(r8) :: Deff_e = 0._r8
      REAL(r8) :: Deff_i = 0._r8
      REAL(r8) :: Deff_h = 0._r8
      REAL(r8) :: Deff_epll = 0._r8
      REAL(r8) :: Deff_ipll = 0._r8
      REAL(r8) :: Deff_hpll = 0._r8
      REAL(r8) :: F_dvac=-1  ! 1.e5 is for large (hyper) diffusion at edge
      REAL(r8) :: F_dexp=100 ! ramps up in last cell
!-----------------------------------------------------------------------
!     B for synchrotron radiation time in a homogeneous plasma

      REAL(r8) :: B_synch

!-----------------------------------------------------------------------
!     input flag for controlling whether or not time is normalized to
!     the collision time in continuum relativistic electron calculations

      LOGICAL :: tau_c_rel_norm

!-----------------------------------------------------------------------
!     controls boundary conditions on F.  Simple options of zero change
!     (homogeneous) or solving on the boundary (none) exist at present.

      CHARACTER(16) :: Fele_bc="none" ! "homogeneous"
      CHARACTER(16) :: Fion_bc="none" ! "homogeneous"
      CHARACTER(16) :: Fhot_bc="none" ! "homogeneous"

!-----------------------------------------------------------------------
!     controls for terms in Ramos CEL-DKE
!     qpi_drive -> turns on Pi_|| drive
!     qpe_drive -> turns on Pe_|| drive
!     ppi_drive -> turns on Pi_|| drive
!     ppe_drive -> turns on Pe_|| drive
!     rpe_drive -> turns on Fe_|| drive
!     gei_drive -> turns on Gei drive
!     dr_drive ->  turns on drift drives
!     te_drive ->  turns on T_e drives
!     ti_drive ->  turns on T_i drives
!     ve_drive ->  turns on electron flow drives
!     vi_drive ->  turns on ion flow drives

      LOGICAL :: qpi_drive = .false.
      LOGICAL :: qpe_drive = .false.
      LOGICAL :: ppi_drive = .false.
      LOGICAL :: ppe_drive = .false.
      LOGICAL :: rpe_drive = .false.
      LOGICAL :: gei_drive = .false.
      LOGICAL :: dr_drive = .false.
      LOGICAL :: te_drive = .false.
      LOGICAL :: ti_drive = .false.
      LOGICAL :: ve_drive = .false.
      LOGICAL :: vi_drive = .false.

!-----------------------------------------------------------------------
!     T_sfac -> controls weight of parallel temperature gradient drive
!     (T_sfac - s^2) v_|| dot grad ln T f_M
!     2.5 is correct form when flow is evolved.
!     1.5 seems more correct when only T and F are evolved.
      REAL(r8) :: T_sfac = 2.5 ! 1.5

!-----------------------------------------------------------------------
!     FLUID NEUTRAL MODEL
!     As opposed to the massive gas injection or pellet injection
!     parameters below, this is meant more for the modeling of
!     neutrals coming in from the edge.
!     See: Meier, Glasser, Lukin, Shumlak JCP 231 (7), 2963  2012
!     Meier, Shumlak PoP 19, 072508 (2012)
!
!     neutral_model has four options with respect to the
!     coupling coefficients (ionization, recombination, charge exchage):
!       (1) none - skip neutral fluid equations
!       (2) const_coef - use constant coefficients
!       (3) n0_coef - use n=0 dependent coefficients
!       (3) full - full 3D coefficients
!
!     Neutral_rates for ionization and recombination can be set by the
!     default rates or using a bicubic spline fit to the degas2 tables.

      CHARACTER(16) :: neutral_model="none"
      REAL(r8) :: ion_eff=0.           ! Effective ionization potential
      REAL(r8) :: sion_fac=0.          ! Ionization factor
      REAL(r8) :: srec_fac=0.          ! Recombination factor
      REAL(r8) :: scx_fac=0.           ! Charge-exchange factor
      REAL(r8) :: ndn_diff=0.          ! Neutral density diffusivity m^2/s
      REAL(r8) :: vn_iso_visc=0.       ! Neutral viscosity m^2/s
      REAL(r8) :: vn_kin_visc=0.       ! Neutral viscosity m^2/s
      CHARACTER(16) :: neut_vel_bc='no-slip' !'free-slip'
      CHARACTER(16) :: neutral_rates="default" ! degas2
      REAL(r8) :: k_perpn=0            ! Set below, not input

!-----------------------------------------------------------------------
!   IZZO parameters used for the kprad gas injection routines

      INTEGER(i4) :: zimp=0_i4           ! impurity atomic number
                                         ! zimp=0 ==> dont use kprad
      LOGICAL :: re_orbits=.false.       ! integrate RE drift orbits

!=======================================================================
!     numerical parameters.
!=======================================================================
!     The following are general run-controlling parameters.
!     The cpu time is checked for each processor during parallel runs,
!     and it's not the total time across processors.  Note that
!     nim_stop is called at the end of the first time step that
!     crosses cpu_tmax, so allow a small margin when setting this
!     value and your batch limit.

      REAL(r8) :: dtm=1.e-8      ! maximum time step, seconds
      REAL(r8) :: dt_initial=0   ! initial time step if nonzero
      REAL(r8) :: dt_stop=1.e-4  ! code stops if dt<=dt_stop*dtm
      REAL(r8) :: tmax=10.0      ! simulation time limit, seconds
      REAL(r8) :: cpu_tmax=1.e9  ! approximate cpu time limit, second
      INTEGER(i4) :: nstep=2000  ! limit on number of time steps
      INTEGER(i4) :: npc=1       ! presently not operable

!-----------------------------------------------------------------------
!     The next set of parameters influence the finite element
!      discretization.
!
!     The following parameter controls whether 1) pressure used in the
!      velocity advance is computed from the nT product at the nodes
!      of the FE expansion and interpolated to the quad points or 2)
!      n and T are interpolated then combined.

      CHARACTER(16) :: p_computation='at nodes'  !  or 'at quads'

!-----------------------------------------------------------------------
!     ngr sets the number of Gaussian quadrature points used in each
!     direction of a rectangular cell.  The number of quad points
!     is set by the function:
!
!     #_rblock_quad_pts_per_dir=ngr+poly_degree-1
!
!     so that it automatically increases with the polynomial
!     degree of the expansions for the dependent variables.  (see
!     poly_degree above).  ngr<2 may produce singular matrices.
!
!     The integration_formula lets you choose the type of numerical
!     integration, either Gaussian (no end-point abscissas) or
!     Lobatto (abscissas located at the GLL nodes if ngr=2, including
!     end points).  However, degenerate rblocks and rblocks touching
!     the geometric axis always use Gaussian to avoid mass matrices
!     with zeros on the diagonal at end points where the Jacobian is
!     zero.  To get Lobatto in other blocks, set nxbl>1.
!
!     met_spl selects what form the partial derivatives of the grid
!     have within the rblocks (bicubic, bilinear, or piecewise
!     constant).

      INTEGER(i4) :: ngr=2
      CHARACTER(5) :: met_spl='iso' !  'iso', 'liner', or 'pcnst'
      CHARACTER(8) :: integration_formula='gaussian'  !  or 'lobatto'
      CHARACTER(8) :: surface_int_formula='gaussian'  !  or 'lobatto'

!-----------------------------------------------------------------------
!     rwinterp set the method used to interpolate chi from the GRIN
!     basis to the NIMROD basis. Options:
!     1) average - a uniformly weighted average of width 1/poly_degree
!                  around each NIMROD node.
!     2) spline  - use cubic splines fitted to the GRIN nodes to
!                  iterpolate to the NIMROD nodes.
!     3) cylgll  - for the cylindral cases only. compute the matrix at
!                  the exact NIMROD GLL node locations.
!     4) nimbnd  - use nimbnd matrix without interpolation
!     5) nimbnd-derivmat - use nimbnd with a derivative matrix

      CHARACTER(16) :: rwinterp='average'

!-----------------------------------------------------------------------
!     Note:  the conform and lumping option are disabled permanently.

      LOGICAL :: conform=.true.
      LOGICAL :: lump_all=.false.

!-----------------------------------------------------------------------
!     Time-discretization parameters:
!       fom is the time-centering on the Lorentz forces.
!       feta is the time centering on resistive diffusion.
!       fvsc is the centering on viscous dissipation.
!       fthc is the implicit centering for thermal conduction.

      REAL(r8) :: fom=1  ! presently not operable
      REAL(r8) :: feta=1
      REAL(r8) :: fvsc=1
      REAL(r8) :: fthc=1

!-----------------------------------------------------------------------
!     Time-centering parameters for advective terms:  fb_vxb is the
!     time-centering (fraction from the predictor step used in the
!     corrector step) of the b in vxb of Ohm's law. fv_vdgv, fn_vdgn,
!     and fp_vdgp are similar centerings for v in momentum advection,
!     n in continuity, and p in the pressure equation, respectively.
!     All of these parameters should be <=1 (there are possible
!     exceptions when there is dissipation), but the minimum value
!     depends on v_cfl, which limits the time step to
!     v_cfl*(cell size)/V.
!
!     1D planar analysis (no dissipation) demonstrates that numerical
!     stability associated with advection is satisfied when
!     1-2*f+f^2*cfl^2<=0, where f is the centering parameter.  Here
!     cfl is an effective cfl [v*dt/(length scale)].  For the Fourier
!     direction, 1/(length scale) -> k_max.  For the poloidal plane,
!     the distributed mass matrix introduces an extra factor of
!     sqrt(3): 1/(length scale) -> sqrt(3)/(cell size).  These
!     factors are now computed within nimrod for the appropriate
!     components.

      REAL(r8) :: fb_vxb=0.54
      REAL(r8) :: fv_vdgv=0.54
      REAL(r8) :: fp_vdgp=0.54
      REAL(r8) :: fn_vdgn=0.54
      REAL(r8) :: v_cfl=0.5

!-----------------------------------------------------------------------
!     The implicit leapfrog algorithm for HMHD requires time-centered
!     implicit advection, and it is always used for ohms/='mhd'.
!       However, the following allows a choice of predictor/corrector
!     advection (with the coefficients defined just above) or
!     time-centered implicit advection.  This allows us to
!     reproduce earlier computations with P/C advection, but it leads
!     to additional coding complications and may be removed at some
!     point.

      CHARACTER(8) :: mhdadv_alg='precor' ! or 'centered'

!-----------------------------------------------------------------------
!     Switch to advance B and T concurrently. Requires beta>0,
!     and mhdadv_alg='centered'.
!     Implementation is currently only nonlinear=F.
!     Switch to advance V and Fhot concurrently.  Requires nF(3)>0,
!     and mhdadv_alg='centered'.
!     Implementation is currently only nonlinear=F.

      LOGICAL :: concur_bt_adv=.false.
      LOGICAL :: concur_vF_adv=.false.

!-----------------------------------------------------------------------
!     The semi-implicit operator in the present formulation may not be
!     effective for stabilizing nonlinear activity.  nl_cfl_lim limits
!     the time step so that a cfl condition based on wave speeds
!     computed with nonlinear pressures is satisfied.  The default
!     is effectively no limit.

      REAL(r8) :: nl_cfl_lim=1.e10

!     Nonlinear velocity advances may use either the standard semi-
!     implicit operator that is described in the 2004 JCP paper or a
!     3D version.  The standard operator acts separately on each
!     Fourier component, which is convenient for the linear algebra.
!     The 3D operator more accurately represents toroidal deformations
!     and improves accuracy at large time-step when toroidal variation
!     is substantial.

      CHARACTER(8) :: siop_type='standard'    !   or '3D'

!-----------------------------------------------------------------------
!     Semi-implicit coefficients:  In general, unity gives the minimum
!       coefficient based on linear numerical stability analysis.
!       Factors>1 are necessary to stabilize nonlinear terms.
!       si_fac_mhd is a coefficient for the semi-implicit
!       operator used for the vxb term in Faraday's law.  si_fac_pres is
!       the coefficient for the pressure advance.  si_fac_j0 is a
!       coefficient for the parallel and perpendicular j0 terms in the
!       velocity equations (j0Xb_tilde and v.grad(P0)).  si_fac_nl is
!       a coefficient for an isotropic operator in velocity, whose
!       magnitude is based on the nonlinear pressures.

      REAL(r8) :: si_fac_mhd=1
      REAL(r8) :: si_fac_pres=1

      REAL(r8) :: si_fac_j0=1
      REAL(r8) :: si_fac_nl=1.1

!-----------------------------------------------------------------------
!     At present, si_fac_hall is only used in the preconditioner for
!       nonlinear computations with ohms='2fl' and multiplies the
!       partial(J)/partial(t) term to increase diagonal dominance when
!       its value is greater than unity.  It should not affect the
!       numerical result provided that the solver tolerance is
!       sufficiently small.

      REAL(r8) :: si_fac_hall=2.

!-----------------------------------------------------------------------
!     The parameter mhd_si_iso is used to control the degree of isotropy
!     for the linear semi-implicit mhd operator.  It must satisfy
!     0<=mhd_si_iso<=1.  The operator is fully isotropic
!     when this factor is 1 and fully anisotropic when this factor is
!     zero.  This adds to the nonlinear isotropic operator controlled
!     by si_fac_nl.

      REAL(r8) :: mhd_si_iso=0.

!-----------------------------------------------------------------------
!     Split resistive diffusion from the mhd part of the magnetic field
!     advance. This may improve numerical stability at larger time
!     steps at the possible expense of the order of convergence in dt.

      LOGICAL :: split_resist=.false.    !  this option is no longer
                                         !  functioning.  always = F

!-----------------------------------------------------------------------
!     Split the viscous forces from and advective & jxb in the velocity
!     advance?  This will improve numerical stability at larger time
!     steps at the possible expense of the order of convergence in dt.

      LOGICAL :: split_visc=.false.

!-----------------------------------------------------------------------
!     Divergence of B cleaning:  divbd is the diffusivity used with
!     the grad(div(B)) term, and fdivb is the time-centering parameter
!     for the implicit solve.  ndivb is the cycle frequency for
!     calling the divergence diffusing equation.  If split_divb is
!     true, the divb cleaning will be done in a separate time split
!     after the rest of the time advance.  Otherwise it is applied
!     in the mhd time split AND in the hall time split if there is
!     one.
!
!     [ndivb is only applicable when split_divb=true; if not,
!     divergence cleaning is applied every time step.]

      REAL(r8) :: divbd=0
      REAL(r8) :: fdivb=1
      INTEGER(i4) :: ndivb=1
      LOGICAL :: split_divb=.false.
      REAL(r8) :: divbdc=0

!-----------------------------------------------------------------------
!     Limit on the tolerated wavenumber**2.  This should be set to
!     1/L**2, approximately, where L is the length scale of the
!     domain.

      REAL(r8) :: kdivb_2_limit=1

!-----------------------------------------------------------------------
!     This version of nimrod allows a separate discontinuous expansion
!      for cleaning div(B), employing a mixed finite-element method
!      with parabolic correction.  The auxiliary
!      discontinuous scalar is of polynomial degree poly_divb and is
!      allocated if the parameter is set >=0.
!
!      The poly_divb_min and poly_divb_max parameters are used
!      to indicate that the modal expansion for the auxiliary field
!      is incomplete.  Where these limits are used, the bases are
!      the outer product of a set of Legendre polynomials from 0
!      to poly_divb in the x logical coordinate and the Legendre
!      polynomials from poly_divb_min to poly_divb_max in the
!      y logical coordinate.  The union of these 2D polynomials and
!      ones that have x and y swapped defines the incomplete expansion.
!
!      The parabolic correction method solves the system,
!        dB/dt = -curl(E) + disc_dbd*grad(auxb)
!         auxb = div(B)
!
!      The method uses fdivb for implicit centering.

      REAL(r8) :: disc_dbd=0._r8   !  in L**2/time
      INTEGER(i4) :: poly_divb=-1  !  set to poly_degree-1 if used
      INTEGER(i4) :: poly_divb_min=-1  !  best left at -1 for a
      INTEGER(i4) :: poly_divb_max=-1  !  complete lower-order expansion

!-----------------------------------------------------------------------
!     When hyp_eta is greater than zero, a numerical hyper-resistivity
!       is applied in the advance of magnetic field.  The coefficient
!       is uniform and time-independent.  It has units of L**4/time.
!       Note that the magnetic advance requires solution of a 6-vector,
!       so it takes much more time per solve.

!       The split_hypeta option allows a choice of time-splitting the
!       hyper-resistive diffusion as a separate solve.

      REAL(r8) :: hyp_eta=0._r8
      REAL(r8) :: fhyp_eta=1._r8  !  temporal centering parameter
      LOGICAL :: split_hypeta=.false.

!     The hyp_dbd is a hyper-diffusivity for divergence cleaning that
!       is similar to the hyper-resistivity, except that the operator
!       is grad(div(grad(div()))).  The split_hypeta input affects this
!       computation, too.

      REAL(r8) :: hyp_dbd=0._r8
      REAL(r8) :: fhyp_dbd=1._r8

!-----------------------------------------------------------------------
!     Computations can use auxiliary scalar fields with incomplete
!      discontinuous representations to alter convergence properties.
!      The two fields described here work like the correction method
!      for magnetic divergence described above, but they act on
!      the momentum-density evolution.  One helps stabilize flow
!      divergence, and the other helps stabilize parallel vorticity.
!
!      When used in a computation, the modified system is (loosely)
!
!         rho*dV/dt = F + sqrt[ddivv*(B^2/mu0+gamma*P)*dt]*grad(auxv)
!                       + sqrt[dpvrt*dt/mu0]*BXgrad(auxw)
!         auxv = sqrt[ddivv*(B^2/mu0+gamma*P)*dt]*div(V)
!         auxw = sqrt[dpvrt*dt/mu0]*B.curl(V)
!
!      where ddivv is a normalized viscosity for divergence error,
!      and dpvrt is a normalized viscosity for vorticity error.
!      This description is not precise, because the implementation is
!      defined directly in weak form and not derived from the strong
!      form shown above.
!
!      The representations need to be incomplete so that they do not
!      duplicate the physics responses to compression
!      and parallel vorticity in numerically resolved evolution.
!      Thus, these auxiliary fields are expanded in incomplete modal,
!      Legendre polynomial representations.  For each of the two
!      logical coordinates in a 2D element, the polynomials from
!      poly_divv_min to poly_divv_max in one coordinate and from 0
!      to poly_divv in the other (and vice versa).  The same holds
!      for the auxiliary field for parallel-vorticity.
!
!      As described in the Aug. 2013 team meeting presentation on
!      projection, the most useful incomplete expansion has
!      poly_divv=poly_divv_min=poly_divv_max=poly_degree.  If the
!      poly_divv_auto input is set to true, the three separate
!      parameters are set to poly_degree automatically.
!
!      Speeds for the artificial waves can be comparable to
!      characteristic perpendicular and shear waves, respectively.
!
!      fdivv and fpvrt are implicit centering paramters for these
!      terms.

      REAL(r8) :: ddivv=1._r8
      REAL(r8) :: dpvrt=0.1_r8
      REAL(r8) :: fdivv=1.0_r8
      REAL(r8) :: fpvrt=1.0_r8
      INTEGER(i4) :: poly_divv=-1
      INTEGER(i4) :: poly_divv_min=-1
      INTEGER(i4) :: poly_divv_max=-1
      LOGICAL :: poly_divv_auto=.false.

!      NOTE: it is assumed that the two auxiliary fields have the
!      same expansion, so they are combined into a single two-vector,
!      which simplifies coding.
!-----------------------------------------------------------------------
!     The following parameters affect how frequently matrices and their
!     preconditioning factors are computed in nonlinear cases.  first,
!     nearly all of the matrices depend on dt.  thus, dt_change_frac
!     is used to force the time step to decrease by this fraction when
!     CFL conditions are limiting it.  it also prevents increases
!     until a change of dt_change_frac is allowable by the CFL
!     conditions.  this prevents time step from changing every cycle.
!
!     The n_dt_release parameter is a release from waiting for the CFL-
!       allowable dt to go back up by dt_change_frac.  After
!       n_dt_release cycles, it is allowed to increase, anyway.  This
!       had been a parameter in subroutine new_dt.
!
!     When the time-step is allowed to increase, it is limited to
!       changing by a factor of dt_incr.  This helps keep the time-step
!       from repeatedly bumping against the CFL limit with subsequent
!       drops by a fraction of dt_change_frac.
!
!     The semi-implicit operators used for advancing B and P depend on
!     the equilibrium fields plus the n=0 part of the perturbed
!     solution.  these matrices are recomputed if the fractional
!     change in the sum of the two changes by more than
!     ave_change_limit.
!
!     For nonlinear computations with slowly changing fields, it is
!       useful to regularly update matrices regardless of the tests
!       on how much the fields have changed.  n_mat_update sets this
!       frequency in number of steps.

      REAL(r8) :: dt_change_frac=0.1
      REAL(r8) :: ave_change_limit=0.01
      REAL(r8) :: dt_incr=1.07
      INTEGER(i4) :: n_dt_release=10
      INTEGER(i4) :: n_mat_update=100000

!-----------------------------------------------------------------------
!     If normalize is true normalization of the primary field variables
!       is applied before and after dump reads/writes. Fields,
!       diffusivities, constants, and times are normalized such that
!       n -> n/n0, v -> v/vA, B -> B/B0, p -> p/p0, T -> T/T0 where
!       n0 and B0 are determined from k_prp_ref_b and k_ref_n if
!       available. See normalize.f90 for specifics.

      LOGICAL :: normalize=.FALSE.

!-----------------------------------------------------------------------
!     The following parameters enhance local particle diffusion in
!       regions where number density drops below a 'floor' value.  The
!       extra diffusivity is
!
!         nd_dart_fac*(1/dt)*[cross_section/(mx*my)]*
!           MAX(0._r8,TANH( (nd_floor*ndens-n)/(nd_exp*ndens) ) )
!
!       where the factor in brackets is an average element area.  The
!       computation uses the minimum n over toroidal angle, and the
!       diffusivity is just a function of position in the poloidal
!       plane.
!
!       The input parameters are:
!
!       nd_floor    -- the floor value is nd_floor*ndens
!       nd_exp      -- nd_exp*ndens determines the transition 'width'
!                      of the diffusivity enhancement in terms of n
!       nd_dart_fac -- sets the magnitude of the diffusivity enhancement
!                      in terms of cell area per dt.
!
!       The default nd_floor (=0) skips this enhancement, and it's only
!       active when nd_diff>0.

      REAL(r8) :: nd_floor=0._r8
      REAL(r8) :: nd_exp=0.02_r8
      REAL(r8) :: nd_dart_fac=1._r8

!-----------------------------------------------------------------------
!     This is a second artificial diffusivity that reproduces an
!       upwinding-like smoothing.  The coefficient nd_dart_upw is
!       dimensionless, and it is multiplied by
!
!       (jac/dt)*[ (dt*V.grad(n)/n_tot)**2 +
!                  MAX(0,(nd_floor_upw-n_tot)/nd_width_upw)) ] ,
!
!       where jac is the local 2D jacobian for the poloidal plane.  This
!       makes an effective diffusivity of jac/dt that becomes active
!       in places where (dt*V/L_n)**2 is significant and where
!       n_tot drops below nd_floor_upw.  The diffusion is parallel to
!       the flow only with upw_aniso=1 (see below), and it is only
!       applied to nonlinear cases with implicit advection.  Here, the
!       diffusivity is fully 3D.

      REAL(r8) :: nd_dart_upw=0._r8
      REAL(r8) :: nd_floor_upw=0._r8   !  [m**-3 unlike nd_floor]
      REAL(r8) :: nd_width_upw=1._r8   !  [m**-3 unlike nd_exp]

!-----------------------------------------------------------------------
!     Artificial thermal diffusivities for an upwinding-like smoothing
!       are also available for the temperature advance.  The coefficient
!       t_dart_upw is used in the same way as nd_dart_upw for number
!       density.  Here, the diffusivities
!
!       (jac/dt)*[ (dt*V.grad(T)/T_tot)**2 +
!                  MAX(0,(t_floor_upw-T_tot)/t_width_upw)) ] ,
!
!       use the species T and V if separate_pe is true; though, the VV
!       dyad used to create the preconditioner matrix is formed from
!       the COM flow only.  Again, this is only used in nonlinear
!       computations with implicit advection, and the diffusivity is 3D.

      REAL(r8) :: t_dart_upw=0._r8
      REAL(r8) :: t_floor_upw=0._r8   !  in eV
      REAL(r8) :: t_width_upw=1._r8   !  in eV

!-----------------------------------------------------------------------
!     The following parameters alter the upwinding-like diffusivities
!       for number density, temperature, etc.  upw_aniso controls
!       whether the diffusion is purely along the direction of flow.
!       It should have values ranging from 0 to 1 with 0 being
!       isotropic, and 1 being fully anisotropic.  upw_limit sets
!       a limit on (dt*V.grad(n)/n)**2 and (dt*V.grad(T)/T)**2.  Use a
!       value less than 1 and increase nd_dart_upw and t_dart_upw to
!       enlarge the region over which the artificial diffusivity is
!       significant.

      REAL(r8) :: n_upw_aniso=1.
      REAL(r8) :: n_upw_limit=1.
      REAL(r8) :: t_upw_aniso=1.
      REAL(r8) :: t_upw_limit=1.
!-----------------------------------------------------------------------
!     Either of the artificial particle diffusivities (nd_diff and
!	nd_hypd) leads to violation of conservation of momentum and
!	energy.  The following option subtracts V*S_n from the flow-
!	velocity equation and adds (mV^2/2-k_B*T/gamma-1)*Sn to the
!	temperature equation to correct these errors.  The effective
!	source/sink factor is
!	  Sn=div(nd_diff*grad(n)-nd_hypd*grad(grad^2(n))).

      LOGICAL :: nd_correrr=.false.
!-----------------------------------------------------------------------
!     The following parameters set minimum or "floor" values directly
!     on the values of density and temperature at the nodal locations
!     of the spectral-element/Fourier expansions.  This is only used
!     in nonlinear computations.  Unlike diffusion, this approach is
!     not conservative and should only be used to control the density
!     and/or internal energy where they are negligibly small, such
!     as artificial cold-plasma regions.
!
!     Note that dealiasing, and dropping the largest-possible Fourier
!     coefficient when dealiase=F, make our FFTs non-invertible.
!     Thus, a specified minimum and what will be realized after the
!     forward FFT will differ.  A test can be used to compare the
!     two (change the check_min parameter in set_nodal_min in
!     utilities.f), but it performs an extra FFT.
!
!     This option is only used when the paramters are set>=0.

      REAL(r8) :: nd_nodal_floor=-1._r8  !  in units of 1/m^3
      REAL(r8) :: ti_nodal_floor=-1._r8  !  in units of eV
      REAL(r8) :: te_nodal_floor=-1._r8  !  used when separate_pe=T

!-----------------------------------------------------------------------
!     For simulations that have vertices along r=0, regularity
!     conditions must be satisfied for the different Fourier
!     components.  One of these conditions is that d/dr of any
!     radial or toroidal vector component of n=1 should be 0.  This
!     is not satisfied through usual finite element natural boundary
!     conditions, since the effective surface area is 0.  Instead,
!     a volume constraint integral is added to the matrices to
!     inhibit the growth of d/dr for these components.  The
!     coefficient r0dr_weight is a factor for the integrand, which
!     is proportional to r**2*(d/dr)**2 along cells touching r=0.
!     The integral is also scaled to diagonal elements of each matrix
!     to which it is added.
!
!     If these components appear to violate the Neumann condition
!     excessively in a particular simulation, increase r0dr_weight.
!     Conversely, if the conditions distort the interior of the
!     solution, decrease r0dr_weight.  In any case, the influence
!     of the constraint integral decreases with increasing radial
!     resolution.

      REAL(r8) :: r0dr_weight=1        !  NO LONGER OPERATIONAL

!-----------------------------------------------------------------------
!     Use natural magnetic field boundary conditions so that divb
!       cleaner can be applied along the boundary.  This is intended as
!       a pre-solve procedure and must be run with nonlinear=.FALSE.
!       For this to be successful it should be run for a few short
!       time-steps (dtm=10^-13) with large divbd=elecd (10^13) and no
!       other physical effects present (ohms='mhd', eta_model='fixed'
!       hyp_eta=hyp_dif=0, etc).

      LOGICAL :: clean_bnd_divb=.FALSE.

!=======================================================================
!     algebraic solver parameters.
!=======================================================================
!     The following control the iterative solver.
!     The choice of preconditioner is
!     determined by solver.  Existing selections include:
!     1) "diagonal"  where the inverted diagonal elements are the
!           sub-matrices consisting of all quantities located at
!           each grid vertex, e.g., the real and imaginary vector
!           components of B at each vertex for the B-field solve.
!     2) "lapack"  a banded serial direct solve over all grid-blocks
!           (rblocks only, tblocks use diagonal) using Lapack
!           routines.  Here, the cg is only used for iterative
!           refinement on 2D matrices. The old "bl_drect" option
!           now reverts to "lapack".
!     3) "bl_ilu_n"  where n is "0" or "1".  two versions of
!           incomplete factorization within each grid-block (again
!           rblocks only).  When n=0, there is no fill-in; the
!           sparsity pattern in the LU factors is the same as the
!           original matrix.  When n=1, there is 'first-level' fill-
!           in; the LU factors have the sparsity pattern of the
!           product of the n=0 LU factors.  Unlike option 2), the
!           actual matrix structure of degenerate rblocks is used,
!           and elements of the factors are vertex sub-matrices
!           (like the diagonal option).  Note that this option
!           has not been updated for high-order elements and only
!           works for bilinear elements.
!     4) "bl_diagx" or "bl_diagy" this is a line-Jacobi iteration,
!           where a direct solve is computed along either the x-
!           or y-direction.  It's mainly intended for single-block
!           problems on vector machines.  Another related choice
!           is "bl_diaga", which is an average of the separate
!           x- and y- direction solves.
!     5) "gl_diaga" global line-Jacobi preconditioning, where the
!           lines extend over contiguous, conforming rblocks to the
!           greatest extent possible.  Parallel decomposition
!           swapping is performed with asynchronous mpi calls.
!           This is probably the best large-problem, multi-block
!           preconditioner that we have at this point.
!     6) "seq_slu"  a sparse serial direct solve over all rblocks
!           using the Sequential SuperLU library, which needs to
!           be linked at the load step.  As with "lapack,"
!           the cg is only used for iterative refinement on 2D mats.
!     7) "slu_dist"  a sparse parallel direct solve over all rblocks
!           using the SuperLU_DIST library, which needs to
!           be linked at the load step.  As with "lapack,"
!           the cg is only used for iterative refinement on 2D mats.

      CHARACTER(8) :: solver="diagonal"

!-----------------------------------------------------------------------
!     The following parameters allow separate choices for specific
!       equations. Any parameters left at their default value will
!       be set to the value of the solver input parameter.

      CHARACTER(8) :: vmhd_solver="none"
      CHARACTER(8) :: bmhd_solver="none"
      CHARACTER(8) :: temp_solver="none"

!-----------------------------------------------------------------------
!     The parameter maxit is the maximum number of CG iterations or
!       GMRES basis vectors allowed, and tol is the relative convergence
!       tolerance.
!     If maxit is set too small some simulations may face prohibitively
!       small timesteps, if it is too large then users may face other
!       problems and warnings will be issued.

      INTEGER(i4) :: maxit=50
      REAL(r8) :: tol=1.e-8

!-----------------------------------------------------------------------
!     The existing ilu and line-solve based (like bl_diagx
!       and bl_adi_y) preconditioners do not use pivoting and are
!       unstable for matrices with large condition numbers.
!       off_diag_fac is used to improve the condition number (of the
!       factored matrix only).  The off-diagonal sub-matrices are
!       multiplied by this factor prior to the incomplete factorization.
!       If you run into   problems with the factorization, reduce
!       off_diag_fac.  If you are adventurous, try increasing it.
!       [Both ilu and line factorizations now have a loop to find the
!       largest factor possible.  off_diag_fac is now used as a lower
!       limit.]

      REAL(r8) :: off_diag_fac=0.9

!-----------------------------------------------------------------------
!     Old solution vectors from each equation are saved to extrapolate
!     guesses.  The extrapolation is based on a polynomial that
!     passes through the old data.  extrap_order is the order of this
!     polynomial.  As its value is increased, there is potentially
!       more accuracy in the guesses at the expense of more memory.

      INTEGER(i4) :: extrap_order=1

!-----------------------------------------------------------------------
!     When nsym_*pre_band is greater than zero, nimrod generates
!       bands representing Fourier-component coupling above and below
!       the diagonal for use in block SOR preconditioning during the
!       nonsymmetric 3d solves.  The nsym_*pre_rpass sets the number of
!       passes, which use alternating (up and down in F-comp) direction
!       between even and odd numbered iterations.  The relaxation factor
!       is nsym_pre_rfac, and nsym_pre_rtype indicates whether to use
!       Jacobi- or Gauss-Seidel-style passes.
!
!       nsym_tpre_band corresponds to the (electron if split) temperature
!       advance, and nsym_bpre_band controls the magnetic-field advance.
!       If both are greater than zero, they must be equal.

      INTEGER(i4) :: nsym_bpre_band=0
      INTEGER(i4) :: nsym_tpre_band=0
      INTEGER(i4) :: nsym_pre_rpass=2
      REAL(r8) :: nsym_pre_rfac=2._r8/3._r8
      CHARACTER(12) :: nsym_pre_rtype="Gauss Seidel"  !  or "Jacobi"

!-----------------------------------------------------------------------
!     Note that the total number of GMRES iterations (the outer loop
!       over these preconditioning passes) is independent of layer
!       decomposition only when using Jacobi-style passes.  When using
!       Gauss-Seidel-style passes, only Fourier couplings within a layer
!       use the most current data; couplings outside the layer use data
!       from the previous iteration, like Jacobi.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     The following parameters control convergence on nonlinear
!       algebraic systems in the implicit leapfrog algorithm.  In
!       particular, the delta_V.grad(delta_V) term can be iterated in
!       the COM velocity advance, and the delta_Jxdelta_B nonlinear
!       Hall term can be iterated in the two-fluid B-advance.

!       The maxit_nl parameter is the maximum number of iterations
!       allowed, but setting it to 1 (the default) is a switch to
!       avoid nonlinear iteration.  The tol_nl specifies the relative
!       tolerance of the 2-norm of the residual.  Note that we
!       must have tol_nl>>tol.

      INTEGER(i4) :: maxit_nl=1
      REAL(r8) :: tol_nl=1.e-4

!-----------------------------------------------------------------------
!     Long listing from iterative solver.

      LOGICAL :: itflag=.false.

!-----------------------------------------------------------------------
!     Hyper/Mumps/Pardiso options (depend on the solver used). See the
!     bridge routines for usage.

      INTEGER(i4) :: iopts(solve_nopts)=-1
      REAL(r8) :: dopts(solve_nopts)=-1._r8
      INTEGER(i4) :: v_iopts(solve_nopts)=-1
      REAL(r8) :: v_dopts(solve_nopts)=-1._r8
      INTEGER(i4) :: b_iopts(solve_nopts)=-1
      REAL(r8) :: b_dopts(solve_nopts)=-1._r8
      INTEGER(i4) :: n_iopts(solve_nopts)=-1
      REAL(r8) :: n_dopts(solve_nopts)=-1._r8
      INTEGER(i4) :: t_iopts(solve_nopts)=-1
      REAL(r8) :: t_dopts(solve_nopts)=-1._r8

      REAL(r8) :: solver_thread_reduction=1._r8

!-----------------------------------------------------------------------
!     If false, -div(nv) is used in the density advance instead of
!       integrating this term by parts.
      LOGICAL :: scal_byparts=.true.

!=======================================================================
!     output parameters.
!     io  parameters.
!=======================================================================
!     These parameters control the restart dumps read and written by
!     nimrod.  The run is initiated by reading dump_file, which was
!     written by nimset or during the course of another nimrod run.
!     It writes dumps into the dump_dir directory with the name prefix
!     dump_name.  ndump controls the dumping frequency in terms of
!     number of time steps, and dump_over tells nimrod how to handle
!     pre-existing dump files with the same name. If h5dump is true
!     then a hdf5 format is used instead of binary (if FCIOWRAPPERS
!     is linked). tdump sets a dump frequency in time where ndump and
!     tdump can be combined.

      CHARACTER(128) :: dump_file="dump.00000"
      CHARACTER(64) :: dump_dir="."
      CHARACTER(64) :: dump_name="dump"
      INTEGER(i4) :: ndump=10000
      REAL(i4) :: tdump=-1._r8
      INTEGER(i4) :: dump_over=0 ! 0: overwrite; 1: append; 2: error
      LOGICAL :: h5dump=.FALSE.
      INTEGER(i4) :: h5io_block_proc_stride=-1_i4
      INTEGER(i4) :: h5io_node_proc_stride=1_i4


!-----------------------------------------------------------------------
!     Parameters for time history plots, created if nhist>0.  This
!     includes the single-point probe, and spatially-integrated
!     energy and div(b) diagnostics.
!
      INTEGER(i4) :: n_probe = 1    ! number of total probes
      INTEGER(i4) :: nhist=0          ! time-history stride
      LOGICAL :: hist_binary = .true. ! true=binary file (text file is always made)
      INTEGER(i4), DIMENSION(10) :: ihist=0   ! i-index init grid for time hists.
      INTEGER(i4), DIMENSION(10) :: jhist=0   ! j-index init grid for time hists.
      LOGICAL :: lin0eq_energy = .true. ! include eq in lin n=0 energy
      LOGICAL :: exn0_divb = .false. ! exclude n=0 divB, only for transfer_eq=T
!-----------------------------------------------------------------------
!     Parameters for xdraw plots--created if the *_stride
!     values are >0.  Note:  these should be set to zero for running
!     in parallel, since the associated routines do not have
!     processor-to-processsor communication.  Intead, use nimplot to
!     create the same files from dump files.

      CHARACTER(64) :: xdraw_dir="."  ! xdraw file directory
      INTEGER(i4) :: xy_stride=0 ! time stride for xy slice
      INTEGER(i4) :: xt_stride=0 ! time stride for xt slice
      INTEGER(i4) :: yt_stride=0 ! time stride for yt slice
      REAL(r8) :: x0fac=0.25 ! relative position of yt slice
      REAL(r8) :: y0fac=0.25 ! relative position of xt slice

!-----------------------------------------------------------------------
!     Dirctory name and time-step cycle frequency for IBM Data
!     Explorer output.  Note:  set these to zero for running in
!     parallel, too (see above).

      CHARACTER(64) :: dx_dir="." ! DX file directory
      INTEGER(i4) :: ndxout=0     ! stride for DX output

!-----------------------------------------------------------------------
!     To include the computed perturbed current, electric field, heat
!     flux, divergence of the stress tensor, equilibrium psi, or the
!     <Beq^2> in
!     the dump file, set these to true.
!     This will increase the size of the dump files, and make them
!     incompatible with those produced (with nimset or nimrod) without
!     dump_XXX=true (unless dump_h5=T).

      LOGICAL :: dump_ja=.false.
      LOGICAL :: dump_eexp=.false.
      LOGICAL :: dump_q=.false.
      LOGICAL :: dump_divpi=.false.
      LOGICAL :: dump_psi_eq=.false.
      LOGICAL :: dump_fsa_beq2=.false.
!-----------------------------------------------------------------------
!     Write the boundary node RZ locations to a text file
      LOGICAL :: write_bdry_nodes = .false.
!-----------------------------------------------------------------------
!     Debugging output for nimset.

      LOGICAL :: detflag=.false. ! print detflaged ascii output?

!=======================================================================
!     extra parameters.
!=======================================================================
!-applied external B-field from coils
!
!     For rectangular grids add magnetic fields due to symmetric
!     toroidal coils to the equilibrium poloidal magnetic field.  The
!     magnetic field components are based on the eqns. 5.40 of Jackson,
!     including second-order terms.

      INTEGER(i4), PARAMETER :: ncoil_max=20
      REAL(r8), DIMENSION(ncoil_max) :: coil_current=5e3
      REAL(r8), DIMENSION(ncoil_max) :: coil_r=0
      REAL(r8), DIMENSION(ncoil_max) :: coil_z=0
      REAL(r8), DIMENSION(ncoil_max) :: coil_r_length=0.05
      REAL(r8), DIMENSION(ncoil_max) :: coil_z_length=0.05
      INTEGER(i4), DIMENSION(ncoil_max) :: coil_r_windings=10
      INTEGER(i4), DIMENSION(ncoil_max) :: coil_z_windings=10
      INTEGER(i4) :: ncoil=0

!=======================================================================
!     particle parameters.
!=======================================================================
!     nm is the requested number of particles in the simulation.
!     However due to roundoffs in the integer arithmetic involved in
!     initializing the particles (can't load 5.352 particles into a cell),
!     the actual number may very by a few.  nm=0 skips all particle
!     subroutines
!
!     part_type(1:5) picks the equations of motion used to advance the
!     particles and can be either
!               drift
!               boris
!               cartb
!     the remain CHARACTERS can be anything
!
!     pnorm_field sets the spatial profile.
!               pres_eq     profile proportional to pres_eq
!                                scaled by betafrac
!               dens_eq     profile proportional to ndens
!                                scaled by betafrac
!               uniformp    uniform normalize to specified pressure
!                                set by nh0*1/2*mass*vhmx^2
!               uniformn    uniform normalize to density
!                                set by nh0
!               gauss       gaussian profile uses [nh0,rh0,hwdth]
!                                nh0*EXP(- (r-rh0)^2/(L*hwdth)
!                                   where L is domain size
!               tanh        profile for ITPA TAE benchmark
!                              nh0*EXP{-(hwdth/R0)*TANH[(s-rh0)/hwdth]}
!               none        this is for tracing nm<10K
!                             use sign of avw0(3)=vbn and psipmax to control
!                               vbn>0 psipmax>0 line in R
!                               vbn<0 psipmax>0 line in Z
!                               vbn>0 psipmax<0 line in phi
!                               vbn<0 psipmax<0 circle load
!                             R0 and psipmax avw0(1,2) control the line position
!                              see particle_interface.f and particle_spatial_load.f
!
!     hot_distro picks the velocity space distribution.  Associate
!     parameters are [vhmx,avu0,ecrit] current options
!     are (see details in particle_velocity_load.f)
!               monoenergetic - all particles have same energy set by vhmx
!               maxwellian - maxwellian distribution there vth=vhmx
!               slowing down - a slowing down distribution function
!                     ecrit is the slowing down critical energy
!                       units are KeV
!                     vhmx is the birth velocity in m/s
!                     Psi0 sets the scale length in the weight equation
!
!               slowing down iso - option for a slowing down distribution function that
!                   is isotropic in velocity space (no pitch angle dependence) both
!                   in the load (which is the assumed norm of the delta-f distribution) and
!                   in the weight equation.
!
!       the default is to load particles isotropically in velocity space
!         however a direction can be specified with a trailing letter
!             see particle_velocity_load.f
!          's'   takes the direction vector from avu0(1:3)
!          'S'   spread in v with MAX(avu0) setting the dominant direction
!                   avw0(3) sets the number of bins
!
!     orbits is the number of subcyles per MHD timestep
!
!     phqty = [2,6] set the number of independent components in pressure tensor
!               2 - CGL like pressure
!               6 - full symmetric pressure tensor
!
!     add_phot = 0 no hot particle coupling to fields
!              = 1 usual phot pressure coupling in vrhs
!              = 2 phot_add2 (particles), isotropic only (continuum)
!              = 3 J_hot current subtracted from total J
!                  in momentum equation and Ohm's law
!
!     pdbg - turns on verbose output useful for debugging with optimizations
!              = 1 for now only level
!
!     dump_part - logical flag to turn off part_dump output, useful for
!       short runs where no restart is expected and output is a large
!       overhead
!
!     restart = 0 does not read in a part_dump file and does a usual load
!     restart /= 0 reads in the part_dump*.bin files and restart is used
!       as the time index of the part_dump file.  It is possible to start
!       a run with inconsistent particle and field information
!     restart < 0 reads in old format of part_dump file if dump_part is .TRUE.
!           i.e. don't timestamp the part_dump filename
!     restart = -1 OVERWRITES part_dump file write if dump_part is .TRUE.
!           i.e. don't timestamp the part_dump filename
!
!     tpartdr is the directory name of the particle output files
!
!     iseed is used to pick the initial prime array used for bit revers
!       algorithm for particle loading
!
!     other parameters are advanced features and the code should be
!       referred to for details of its current incarnation.
!=======================================================================
!
      REAL(r8), DIMENSION(1:3) :: avu0=0.0_r8
      REAL(r8), DIMENSION(1:3) :: avw0=0.0_r8
      REAL(r8) :: mass0=1.6726e-27_r8
      REAL(r8) :: charge0=1.6022e-19_r8
      REAL(r8) :: vhmx=1.0e6_r8
      REAL(r8) :: ecrit=10.0_r8
      REAL(r8) :: psipmax=0.0_r8
      REAL(r8) :: Ppsi0=0.0_r8
      REAL(r8) :: R0=1.0_r8
      REAL(r8) :: betafrac=0.0_r8
      CHARACTER(8) :: reduce_p='both'     ! 'both' (1-betafrac) applied to tion_eq and tele_eq
      CHARACTER(6) :: find_type='exact '  ! 'biline','biquad','bicube','biquar'
      CHARACTER(6) :: field_type='exact ' ! 'biline','biquad','bicube','biquar'
      CHARACTER(6) :: phot_type='biline'  ! 'biquad','bicube','biquar'
      REAL(r8) :: orbits=1.0_r8
      REAL(r8) :: anisop=1.0_r8
      REAL(r8) :: nh0=1.0_r8
      REAL(r8) :: rh0=0.0_r8
      REAL(r8) :: hwdth=1.0_r8
      REAL(r8), DIMENSION(1:16) :: pphparam=0.0_r8
      INTEGER(i4) :: nm=0
      INTEGER(i4) :: eparallel=0          ! controls efld(1:3): 0=zero, 1=only axisym, 2=total
                                          !                    -1=no etaJ in efld(4:6)
      INTEGER(i4) :: restart=0
      INTEGER(i4) :: phqty=0
      INTEGER(i4) :: ph_pd=1
      INTEGER(i4) :: phf_pd=1
      INTEGER(i4) :: bmxs=0
      INTEGER(i4) :: emxs=0
      INTEGER(i4) :: nbins=0
      INTEGER(i4), DIMENSION(1:8) :: vdlms=9999
      INTEGER(i4) :: iseed=1
      INTEGER(i4) :: add_phot=0
      INTEGER(i4) :: pdbg=0
      CHARACTER(50) :: pnorm_field='none'
      CHARACTER(50) :: tpartdr='./'
      CHARACTER(50) :: part_type=' '
      CHARACTER(50) :: phase_use='all'
      CHARACTER(50) :: weight_use='all'
      CHARACTER(50) :: hot_distro=' '
      LOGICAL :: trace=.false.
      LOGICAL :: phdump=.false.
      LOGICAL :: ho_pfield=.false.
      LOGICAL :: dump_part=.true.
      LOGICAL :: scale4phot=.false.
      REAL(r8) :: sv_input=0._r8 ! flow for shifted relmax cases

!-----------------------------------------------------------------------
!     These are useful for development

      CHARACTER(10), DIMENSION(10) :: chrdbug
      LOGICAL, DIMENSION(10) ::       logdbug
      INTEGER(i4), DIMENSION(10) ::   intdbug
      REAL(r8), DIMENSION(10) ::      dbldbug

!-----------------------------------------------------------------------
!     if model_probes=T,
!     the magnetic field components will be written to a file
!     file for list of single-point probes (whose position is given in
!     the file probes.in).  the time-stride is given by nhist and
!     the file format is set by hist_binary.

!     the present implementation outputs the magnetic field components
!     at the probe locations, but additional fields may be included
!     in future updates.

      LOGICAL :: model_probes=.false.

!-----------------------------------------------------------------------
!     Define a meta-variable for deciding how the gridshape is set
!-----------------------------------------------------------------------
      LOGICAL :: experiment_grid = .FALSE.
!-----------------------------------------------------------------------
!     Perturbation data type for more easily handling the pert_* input
!     parameters.  See associated routine for filling it.
!-----------------------------------------------------------------------
      TYPE :: perturbation
        INTEGER(i4) :: n ! Mode number for pe
        INTEGER(i4) :: nuse ! Flag for usage
        INTEGER(i4), DIMENSION(:), POINTER :: m

        ! P: M theta - N phi       R: M theta + N phi
        REAL(r8), DIMENSION(:), POINTER :: Pcos, Psin, Rcos, Rsin
        REAL(r8), DIMENSION(:), POINTER :: Arr,Ari, Air, Aii
        REAL(r8), DIMENSION(:), POINTER :: width
        CHARACTER(10), DIMENSION(:), POINTER :: ptype
      END TYPE perturbation
      TYPE(perturbation), DIMENSION(:), POINTER :: pert
!-----------------------------------------------------------------------
!     close module.
!-----------------------------------------------------------------------
      END MODULE input

!-----------------------------------------------------------------------
!     subprogram 2. read_namelist
!     open and read the namelist input.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.
!-----------------------------------------------------------------------
      SUBROUTINE read_namelist(infile,echo_in,code)
      USE local
      USE input
      USE pardata, ONLY: nprocs,node
      USE io
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: infile, code
      LOGICAL, INTENT(IN) :: echo_in

      INTEGER :: read_stat,nc
      INTEGER(i4) :: number_of_namelists=0,il,inst,ii
      CHARACTER(14) :: tempfile
      CHARACTER(128) :: ctmp
      LOGICAL :: reading, nl_init_exist=.FALSE.
      INTEGER(i4) :: tiopts(solve_nopts)
      REAL(r8) :: k2min,tdopts(solve_nopts)
!-----------------------------------------------------------------------
!     Deprecated variables.  Defined just to allow them to be read in.
!-----------------------------------------------------------------------
      LOGICAL :: lump_b=.false.
      INTEGER(i4) :: n_gm=15 ! stack limit for QGMRES
      REAL(r8) :: vsink_rate=0
      LOGICAL :: advance_pvi=.FALSE.
      CHARACTER(16) :: heat_source="none" ! 'self-similar'
      INTEGER(i4) :: mxpie=1
      CHARACTER(8) :: pieflag="rblock"
      CHARACTER(8) :: rimflag="none"
      INTEGER(i4):: nbl_rim=2
      REAL(r8) :: rim_length=-1.0
      REAL(r8), DIMENSION(15) :: wpack=0._r8
      REAL(r8), DIMENSION(15) :: amp=0._r8
      LOGICAL :: T_only = .false.
      LOGICAL :: si_in_v=.true. ! option no longer functions;
      ! si_in_v is always true.
      LOGICAL :: eqhist = .true.  ! eq energy in nonlinear cases
!-----------------------------------------------------------------------
!     Neoclassical closures are deprecated and do nothing at present.
!      CHARACTER(8) :: neoe_flag="none"
!      CHARACTER(8) :: neoi_flag="none"
      LOGICAL :: neo_debug=.false.
      INTEGER(i4) :: ng_ft=10   ! trapezoidal quadrature for trapped par
!      REAL(r8) :: mu_e=0._r8    ! electon viscous damping frequency
!      REAL(r8) :: mu_i=0._r8    ! ion viscous damping frequency
      REAL(r8) :: neo_rad=0.0_r8 ! Smoothing radius near magnetic axis
      INTEGER(i4) :: m_neo = 2  ! poloidal mode number (Hole)
      INTEGER(i4) :: n_neo = 1  ! toroidal mode number (Hole)
      CHARACTER(8) :: neoe_det="all"
      CHARACTER(8) :: neoi_det="all"
!-----------------------------------------------------------------------
!     These are legacy and compatiblity variables. See discussion below.
      CHARACTER(16) :: nd_bc=""
      CHARACTER(16) :: flow_bc=""
      CHARACTER(16) :: norm_flow=""
      CHARACTER(16) :: magnetic_dirichlet=""
      CHARACTER(16) :: velocity_dirichlet=""
      CHARACTER(16) :: ndensity_dirichlet=""
      LOGICAL :: zero_bnorm=.true.
      LOGICAL :: insulate=.false.
      REAL(r8) :: divbeqd = 0
!-----------------------------------------------------------------------
!     Parameters for resistive wall model.
!     NIMROD uses the thin-shell resistive wall model which allows the
!     surface electric field to be written simply as:
!                 E_tangential = elecd_wall/wall_width normal x [B]
!                   where [B] is the jump in B across the resistive wall
!     elecd_wall=0 gives the perfectly-conducting wall boundary
!     conditions which are the natural boundary conditions if
!     e_vertical, loop_voltage, and i_desired are all 0.
!
!     Notes: currently this gives a constant current boundary condition.
!     It is not clear what would happen if this is combined with
!     parameters that change the current as well.
      REAL(r8) :: elecd_wall = 0
      REAL(r8) :: wall_width = 0
      REAL(r8) :: cw_distance = 0 ! only used in circ and rect
!-----------------------------------------------------------------------
!     Extra variables for time varying loop voltage
      INTEGER(i4) :: loop_volt_ncol, loop_volt_i
!-----------------------------------------------------------------------
!     namelist declaration.  if new namelists are added, increase the
!     number_of_namelists parameter.
!-----------------------------------------------------------------------
!=======================================================================
      NAMELIST/grid_input/gridshape,xmin,xmax,ymin,ymax,mx,my,vac_frac  &
     &        ,nxbl,nybl,firstx,firsty,skew,geom,per_length,nphi,zperiod&
     &        ,periodicity,xo,yo,x1,y1,lphi,lin_nmodes,lin_nmax,nlayers &
     &        ,poly_degree,poly_degree_mesh,poly_distribution,dealiase  &
     &        ,qpack,wpack,amp,xpack,wxpack,ampx,ypack,wypack,ampy      &
     &        ,mxpie,pieflag,rimflag,nbl_rim,rim_length,decompflag      &
     &        ,eqfile,mm_mint,mm_mext,mm_imap_use_int,mm_mdmap          &
     &        ,mm_rw_dmap
!=======================================================================
      NAMELIST/const_input/chrg_input,zeff_input,mi_input,me_input      &
     &        ,gam_input,kblz_input,mu0_input,c_input,mc2_kT_input      &
     &        ,mc2_kTb_input
!=======================================================================
      NAMELIST/init_input/nx,ny,init_type,bamp,reset_file,reset_tol     &
     &        ,rescale_factor,reset_time,reset_step,transfer_eq         &
     &        ,theta_init,phi_init,init_x,init_y,blob_amp,blob_width    &
     &        ,pert_m,pert_n,pert_cos,pert_sin,pert_type,pert_width     &
     &        ,pert_phase,transfer_n0,init_type_transfer_n0             &
     &        ,transfer_overwrite,eq_dump_file,eq_nmodes                &
     &        ,transfer_eq_fields,reset_ve_file,reset_vi_file           &
     &        ,reset_vh_file,bs_type,cfext,cur_fac,per_shift
!=======================================================================
      NAMELIST/equil_input/ndens,nedge,n_profile,npp,be0                &
     &        ,pe_profile,te_cent,te_edge,tepp,p_profile                &
     &        ,thetab,phib,lamprof,lam0,rbreak,alpha,beta               &
     &        ,dvac,dexp,ds_use,ds_function,xvac,dfac,ds_nqty           &
     &        ,eq_flow,ve0,thetav,phiv,eqflow_width                     &
     &        ,pit_0,pit_2,pit_4,pres_2,pres_4,divbeqd                  &
     &        ,glength,tor_eqja_fe,tor_angle,bb_tanh                    &
     &        ,gfilename,cur_state_file,pstok,run_id                    &
     &        ,eta_i,eta_e,telength,tilength,diffuse_eq_var             &
     &        ,diffuse_eq_divfd,diffuse_eq_steps,diffuse_eq_diff        &
     &        ,calc_eq_fe,diffuse_eq_dvbsteps,diffuse_eq_phif           &
     &        ,calc_jeq_fe_nq,calc_jeq_fe_gs,diffuse_eq_report          &
     &        ,diffuse_eq_dsftn,diffuse_eq_dvac,diffuse_eq_dexp         &
     &        ,calc_jeq_fe_bnd,ideal_like,pressure_gauge_nim,mvac_nim   &
     &        ,set_vac_quadpts,eq_loop_voltage                          &
     &        ,neut_dens_profile,neutral_dens                           &
     &        ,eq_neut_flow,vn0,neut_temp_profile,neutral_temp          &
     &        ,f_profile,f_cent,f_edge,fpp,fmod_profile,fmod_amp        &
     &        ,fmod_center,fmod_width,pmod_profile,pmod_amp             &
     &        ,pmod_center,pmod_width,ppp,p_cent,pped_profile           &
     &        ,pped_height,pped_width,pped_psi0,reflect_slab,nped
!=======================================================================
      NAMELIST/physics_input/ndens,nonlinear,ohms,advect                &
     &        ,elecd,kin_visc,no_pe_hall,gun_max,gun_lim0,gun_lim1      &
     &        ,gun_func                                                 &
     &        ,loop_volt,tloopv0,tloopv1,i_desired,loop_rate,loop_rate2 &
     &        ,loop_volt_file,e_vertical,t_e_vert0,t_e_vert1,t_e_vert2  &
     &        ,vsink_rate,zero_bnorm,insulate,separate_pe,pe_frac       &
     &        ,nd_diff,nd_hypd,continuity,eta_model,eta_ref_t,elecd_min &
     &        ,elecd_max,beta,dvac,dexp,ds_use,advance_pvi,xvac         &
     &        ,f_chodura,elecd_chodura                                  &
     &        ,lamprof,lam0,rbreak,alpha,be0,thetab,phib,bamp           &
     &        ,eq_flow,ve0,thetav,phiv,eqflow_width                     &
     &        ,pit_0,pit_2,pit_4,pres_2,pres_4,reset_file               &
     &        ,init_type,nx,ny,iso_visc,par_visc,gyr_visc,gravity       &
     &        ,rt_transfer_eq,mag_tau,norm_flow,flow_bc,nd_bc           &
     &        ,magnetic_dirichlet,velocity_dirichlet,ndensity_dirichlet &
     &        ,magnetic_bc,velocity_bc,ndensity_bc                      &
     &        ,temp_ion_bc,temp_ele_bc                                  &
     &        ,pert_m,pert_n,pert_cos,pert_sin,pert_type,pert_width     &
     &        ,pert_phase,toffset,tperiod,tdown,more_phiplanes          &
     &        ,rf_power_ec,rf_freq_ec,discrete_extrap,mynphi            &
     &        ,lam0rf,rrf,zrf,phirf,delta_pol,delta_tor                 &
     &        ,neutral_model,ion_eff,sion_fac,srec_fac,scx_fac,ndn_diff &
     &        ,vn_iso_visc,vn_kin_visc,neut_vel_bc                      &
     &        ,shepard_alg,hv_coeff,hv_power,hv_time,hv_filter          &
     &        ,vwall,efkphi,efrad,mmax,almcon,reswall_nmax              &
     &        ,neutral_rates,elecd0,nd_diff0,iso_visc0,kin_visc0,use_fp &
     &        ,t_rmp,t_rmp_offset,rmp_mult_fac,rmp_type,rmp_file_type   &
     &        ,b_nonsym,reverse_ip,load_n0_rmp,lagr_edge_rmp,rmp_omega
!=======================================================================
      NAMELIST/closure_input/                                           &
     &        p_model,k_perp,k_pll,k_perpe,k_plle,k_perpi,k_plli        &
     &        ,k_perpn,kpll_mfe,kplle_mfe,kplli_mfe,par_visc_mfe        &
     &        ,divq_mfe_ibp,TFtheta,tequil_rate,ohm_heat,visc_heat      &
     &        ,k_perpi0,k_perpe0,k_perp0                                &
     &        ,tdep_coul_log,coulomb_logarithm,magfac_ele,magfac_ion    &
     &        ,delta0_ele,delta1_ele,gamma0_ele,gamma1_ele              &
     &        ,delta0_ion,delta1_ion,gamma0_ion,gamma1_ion,tdep_tequil  &
     &        ,gamma_heat,heat_source,electric_source,ndensity_source   &
     &        ,heat_ion_source,heat_ele_source,momentum_source          &
     &        ,analytic_rhs,nF,close_nprocs,T_only,eqn_model,kin_model  &
     &        ,pd_xi,nq_xi,C_nq_xi,m_xi,m_xit,m_xip,tpb_min,gridshape_v &
     &        ,pa_basis,betafrac_sdd,bump_a,nu_relax                    &
     &        ,smin,smid,smax                                           &
     &        ,dxi,xi0,tpow,tcoe                                        &
     &        ,ns,sp1,ns1_fixed,s1_fixed,s1_quad                        &
     &        ,ns2,sp2,ns2_fixed,s2_fixed,s2_quad                       &
     &        ,F_init_type,dump_read_F,F_mass_op_flag                   &
     &        ,cel_flow,cel_temp,kinetic_eq                             &
     &        ,qpi_model,closure_model,closure_dump,parvisc_model       &
     &        ,bhat,perpvisc_model                                      &
     &        ,k_pll_max,k_pll_min,k_pll_ref_t,k_prp_ref_b,k_ref_n      &
     &        ,kprp_mnrat,k_cross                                       &
     &        ,pert_m,pert_n,pert_cos,pert_sin,pert_type,pert_width     &
     &        ,pert_phase,toffset,tperiod,tdown,more_phiplanes          &
     &        ,rf_power_ec,rf_freq_ec,discrete_extrap,mynphi            &
     &        ,lam0rf,rrf,zrf,phirf,delta_pol,delta_tor                 &
     &        ,tmodperiod,tmodphase,rf_modulation,rf_extrap             &
     &        ,shepard_alg,t_on,t_off,control,pcsstatus                 &
     &        ,zimp,q_cent,b_nonsym,t_RMP,RMP_type,reverse_Ip           &
     &        ,t_rmp_offset,rmp_file_type,load_n0_rmp,rmp_mult_fac      &
     &        ,odd_pert,odd_pert_str,psi_lcfs,q_x,q_y,v_cent,psi_per    &
     &        ,single_s,steady_state,Deff_i,Deff_e,Deff_h,cel_output    &
     &        ,fcel,fcol,coll_model,it_calls,xi_prec,ss_prec,ss_rhs     &
     &        ,ssoff_prec,F_dexp,F_dvac,re_orbits                       &
     &        ,fder,fint,fadv,qpi_drive,qpe_drive,ppi_drive,ppe_drive   &
     &        ,rpe_drive,dr_drive,te_drive,ti_drive,ve_drive,vi_drive   &
     &        ,Fele_bc,Fion_bc,Fhot_bc,gei_drive                        &
     &        ,old_nF,old_ns,old_ns2,old_ph_pd                          &
     &        ,Deff_ipll,Deff_epll,Deff_hpll,T_sfac,s_poly_exp          &
     &        ,neoe_flag,neoi_flag,mu_e,mu_i,neo_bump_r0,neo_axis_r     &
     &        ,neo_axis_z,neo_bump_amp,lagr_edge_rmp,rmp_omega,B_synch  &
     &        ,tau_c_rel_norm
!=======================================================================
      NAMELIST/numerical_input/dtm,tmax,nstep,npc,fom,si_fac_hall,      &
     &     conform,lump_b,ngr,feta,divbd,fdivb,ndivb,met_spl,           &
     &     si_fac_pres,si_fac_mhd,fvsc,mhd_si_iso,dt_stop,v_cfl,        &
     &     split_resist,fb_vxb,fv_vdgv,fp_vdgp,fn_vdgn,                 &
     &     kdivb_2_limit,dt_change_frac,ave_change_limit,dt_incr,       &
     &     split_divb,nl_cfl_lim,si_in_v,split_visc,cpu_tmax,dt_initial,&
     &     lump_all,si_fac_j0,r0dr_weight,si_fac_nl,transfer_eq,        &
     &     poly_degree,fthc,n_dt_release,mhdadv_alg,concur_bt_adv,      &
     &     nd_floor,nd_exp,nd_dart_fac,nd_dart_upw,t_dart_upw,          &
     &     n_upw_aniso,n_upw_limit,t_upw_aniso,t_upw_limit,nd_floor_upw,&
     &     nd_width_upw,t_floor_upw,t_width_upw,n_mat_update,           &
     &     nd_nodal_floor,ti_nodal_floor,te_nodal_floor,nd_correrr,     &
     &     p_computation,chrdbug,logdbug,intdbug,dbldbug,               &
     &     integration_formula,surface_int_formula,poly_divv_auto,      &
     &     poly_divb,disc_dbd,poly_divv,poly_divv_min,fdivv,fpvrt,      &
     &     poly_divv_max,ddivv,dpvrt,poly_divb_min,poly_divb_max,       &
     &     scal_byparts,normalize,hyp_eta,fhyp_eta,split_hypeta,hyp_dbd,&
     &     fhyp_dbd,siop_type,rwinterp,clean_bnd_divb,concur_vF_adv
!=======================================================================
      NAMELIST/solver_input/tol,maxit,solver,off_diag_fac,extrap_order, &
     &     vmhd_solver,bmhd_solver,temp_solver,nsym_bpre_band,          &
     &     nsym_tpre_band,nsym_pre_rpass,nsym_pre_rfac,nsym_pre_rtype,  &
     &     maxit_nl,tol_nl,b_iopts,b_dopts,iopts,dopts,v_iopts,v_dopts, &
     &     n_iopts,n_dopts,t_iopts,t_dopts,solver_thread_reduction
!=======================================================================
      NAMELIST/output_input/detflag,nhist,ihist,jhist,hist_binary,ndump,&
     &     dump_name,dump_dir,ndxout,dump_over,dump_file,xy_stride,     &
     &     xt_stride,yt_stride,x0fac,y0fac,xdraw_dir,dx_dir,itflag,     &
     &     reset_file,reset_time,reset_step,eqhist,dump_eexp,h5dump,    &
     &     dump_psi_eq,h5io_block_proc_stride,h5io_node_proc_stride,    &
     &     dump_ja,dump_q,dump_divpi,dump_fsa_beq2,tdump,n_probe,       &
     &     write_bdry_nodes,lin0eq_energy,exn0_divb
!=======================================================================
      NAMELIST/particle_input/nm,restart,phqty,phdump,vhmx,             &
     &      avu0,avw0,mass0,charge0,psipmax,Ppsi0,R0,betafrac,orbits,   &
     &      eparallel,anisop,tpartdr,part_type,add_phot,trace,          &
     &      ho_pfield,ph_pd,phf_pd,bmxs,emxs,pnorm_field,ecrit,         &
     &      phase_use,weight_use,hot_distro,nh0,rh0,hwdth,pdbg,         &
     &      pphparam,dump_part,scale4phot,nbins,iseed,vdlms,reduce_p,   &
     &      find_type,field_type,phot_type,sv_input
!=======================================================================
      NAMELIST/extra_input/ncoil,coil_current,coil_r,coil_z,            &
     &     coil_r_length,coil_z_length,coil_r_windings,coil_z_windings, &
     &     model_probes

!-----------------------------------------------------------------------
!     open input file.
!     Remove comments from input file and put into temporary file.
!-----------------------------------------------------------------------
      tempfile='tempnimrod.in'
      CALL rmcoment(infile,tempfile)
      OPEN(UNIT=in_unit,FILE=tempfile,STATUS='OLD',POSITION='REWIND')
!-----------------------------------------------------------------------
!     check namelist file for namelist order and number.
!-----------------------------------------------------------------------
      DO
        READ(UNIT=in_unit,FMT='(a)',IOSTAT=read_stat) ctmp
        IF (read_stat/=0) EXIT
        nc=LEN_TRIM(ctmp)
        IF (nc<1) CYCLE
        ctmp=ADJUSTL(ctmp)
        reading=.false.
        IF (ctmp(1:1)=='&') THEN
          number_of_namelists=number_of_namelists+1
!-----------------------------------------------------------------------
!         trim all but the namelist name.
!-----------------------------------------------------------------------
          DO il=2,nc+1
            IF (ctmp(il:il)/=' ') THEN
              IF (.NOT.reading) inst=il
              reading=.true.
              CYCLE
            ENDIF
            IF (ctmp(il:il)==' '.AND.reading) THEN
              ctmp=ctmp(inst:il-1)
              EXIT
            ENDIF
          ENDDO
          BACKSPACE(in_unit)
!-----------------------------------------------------------------------
!         select and read namelist.
!-----------------------------------------------------------------------
          SELECT CASE(TRIM(ctmp))
          CASE('grid_input')
            READ(UNIT=in_unit,NML=grid_input,IOSTAT=read_stat)
          CASE('const_input')
            set_phys_constants=.true.
            READ(UNIT=in_unit,NML=const_input,IOSTAT=read_stat)
          CASE('init_input')
            READ(UNIT=in_unit,NML=init_input,IOSTAT=read_stat)
            nl_init_exist=.TRUE.
          CASE('equil_input')
            READ(UNIT=in_unit,NML=equil_input,IOSTAT=read_stat)
          CASE('physics_input')
            READ(UNIT=in_unit,NML=physics_input,IOSTAT=read_stat)
          CASE('closure_input')
            READ(UNIT=in_unit,NML=closure_input,IOSTAT=read_stat)
          CASE('numerical_input')
            READ(UNIT=in_unit,NML=numerical_input,IOSTAT=read_stat)
          CASE('solver_input')
            READ(UNIT=in_unit,NML=solver_input,IOSTAT=read_stat)
          CASE('output_input')
            READ(UNIT=in_unit,NML=output_input,IOSTAT=read_stat)
          CASE('particle_input')
            READ(UNIT=in_unit,NML=particle_input,IOSTAT=read_stat)
          CASE('extra_input')
            READ(UNIT=in_unit,NML=extra_input,IOSTAT=read_stat)
          CASE DEFAULT
            CALL nim_stop                                               &
     &        (TRIM(ctmp)//' is an unrecognized namelist.')
          END SELECT
          IF (read_stat/=0) CALL nim_stop                               &
     &      ('Error reading namelist '//TRIM(ctmp)//'.')
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
!         Check to see if it is a 3.0 version namelist
!-----------------------------------------------------------------------
      IF (.NOT. nl_init_exist) THEN
        WRITE(*,*) 'I think you are using the 3.0 namelist format.'
        WRITE(*,*) 'I recommend using bin/nl_convert to update'
        WRITE(*,*) '  your namelist.'
      ENDIF
!-----------------------------------------------------------------------
!     close input file.
!       Delete it since it is the temporary file
!-----------------------------------------------------------------------
      CLOSE(in_unit,STATUS='DELETE')
!-----------------------------------------------------------------------
!     check for old input specifications, inconsistencies, and set
!     any unset defaults.
!-----------------------------------------------------------------------
      IF (.NOT.conform)                                                 &
     &  CALL nim_stop('The conform=F option is not available.')
      IF (.NOT.si_in_v)                                                 &
     &  CALL nim_stop('The si_in_v=F option is not available.')
      IF (lump_b)                                                       &
     &  CALL nim_stop('The lump_b=T option is not available.')
      IF (ngr<2)                                                        &
     &  CALL nim_stop('Setting ngr<2 may produce singular matrices.')
      IF (ndxout>0)                                                     &
     &  CALL nim_stop('DX output is no longer produced directly '//     &
     &  'from nimrod.  Set ndxout to 0.')
      IF (set_vac_quadpts.AND.(mvac_nim<=0.OR.pressure_gauge_nim<0))    &
     &  CALL nim_stop('mvac_nim and pressure_gauge_nim must be set.' // &
     &                'with set_vac_quadpts=.true.')
      IF (set_vac_quadpts.AND..NOT.ideal_like)                          &
     &  CALL nim_stop('set_vac_quadpts requires ideal_like=T.')
      IF (nonlinear.AND.ideal_like)                                     &
     &  CALL nim_stop('ideal_like is linear only.')
      IF (ideal_like) THEN
        n_profile='bump'
        npp(1)=0.99999
        ds_function(1)='psep'
        p_computation='at nodes'
        scal_byparts=.false.
      ENDIF
      IF (continuity/='none'.AND.continuity/='fix profile'.AND.         &
     &    continuity/='n=0 only'.AND.continuity/='full')                &
     &  CALL nim_stop('Continuity value '//TRIM(continuity)//           &
     &  ' is not recognized.')
      IF (continuity=='none'.AND.beta>0) CALL nim_stop('Continuity '    &
     &  //'must not be set to none if beta>0.')
      IF (p_model=='aniso2'.AND..NOT.separate_pe)                       &
     &  CALL nim_stop('Separate_pe MUST be true if p_model=aniso2.')
      IF (split_visc.AND.kin_visc<=0) CALL nim_stop('Set split_visc='   &
     &  //'F when kin_visc=0.')
      IF (split_visc.AND.(iso_visc>0.OR.par_visc>0)) CALL nim_stop      &
     &  ('Isotropic and parallel viscosities do not work with '         &
     &  //'split_visc.')
      IF (eta_model/='fixed'.AND.eta_model/='eta n=0 only'.AND.         &
     &    eta_model/='eta full'.AND.eta_model/='eta full ds'.AND.       &
     &    eta_model/='chodura'.AND.eta_model/='cql3d_sptz')             &
     &  CALL nim_stop('Eta_model value '//TRIM(eta_model)//             &
     &  ' is not recognized.')
      IF (split_resist) CALL nim_stop('The split_resist option is no '  &
     &  //'longer available.')
      IF (met_spl=='cubic'.OR.met_spl=='bicubic')                       &
     &  CALL nim_stop('Met_spl value '//TRIM(met_spl)//' is not '       &
     &  //'available with the improved mapping.')
      IF (transfer_eq.AND..NOT. nonlinear)                              &
     &  CALL nim_stop('Can only transfer_eq when nonlinear.')
      IF (rt_transfer_eq.AND..NOT. nonlinear) CALL nim_stop('Can only'  &
     &  //' rt_transfer_eq when nonlinear.')
      IF (gyr_visc>0.AND.beta<=0)                                       &
     &  CALL nim_stop('Gyroviscosity requires beta>0.')
      IF (k_cross>0.AND.separate_pe.AND..NOT.p_model(1:5)=='aniso')     &
     &  CALL nim_stop('Thermal drifts require the anisotropic p_model.')
      IF (k_cross>0.AND.nonlinear.AND.continuity/='full')               &
     &  CALL nim_stop('Nonlinear thermal drifts require '//             &
     &                'continuity=full.')
      IF (vsink_rate>0)                                                 &
     &  CALL nim_stop('Vsink_rate is no longer available.')
      IF (heat_source/='none')                                          &
     &  CALL nim_stop('heat_source is deprecated.'                      &
     &  //' Change to heat_ion_source')
      IF (concur_bt_adv.AND.mhdadv_alg/='centered')                     &
     &  CALL nim_stop('concur_bt_adv requires mhdadv_alg=centered')
      IF (concur_vF_adv.AND.mhdadv_alg/='centered')                     &
     &  CALL nim_stop('concur_vF_adv requires mhdadv_alg=centered')
      IF (concur_bt_adv.AND.beta<=0)                                    &
     &  CALL nim_stop('concur_bt_adv requires beta>0')
      IF (concur_bt_adv.AND.nonlinear)                                  &
     &  CALL nim_stop('concur_bt_adv not implemented for nonlinear=T')
      IF (concur_vF_adv.AND.nonlinear)                                  &
     &  CALL nim_stop('concur_vF_adv not implemented for nonlinear=T')
      IF(T_only) eqn_model="T"
      IF (MINVAL(pert_m)<0) CALL nim_stop('Cannot have negative pert_m')
      IF (p_model=='iso') p_model='isotropic'
      IF (k_perpi==0.AND.k_perp>0) k_perpi=k_perp
      IF (k_perpe==0.AND.k_perp>0) k_perpe=k_perp
      IF (k_plli==0.AND.k_pll>0) k_plli=k_pll
      IF (k_plle==0.AND.k_pll>0) k_plle=k_pll
      IF (.NOT.separate_pe) k_plli=k_plle
      IF (kplli_mfe==0.AND.kpll_mfe>0) kplli_mfe=kpll_mfe
      IF (kplle_mfe==0.AND.kpll_mfe>0) kplle_mfe=kpll_mfe
      IF (.NOT.separate_pe) kplli_mfe=kplle_mfe
      IF (closure_model(1:3)/="mfe".AND.closure_model(4:6)/="mfe".AND.  &
     &    closure_model(7:9)/="mfe".AND.closure_model(10:12)/="mfe")THEN
        kplli_mfe=0
        kplle_mfe=0
        kpll_mfe=0
      ENDIF
      IF (perpvisc_model=='prpdep'.AND.p_model/='aniso_tdep'.AND.       &
     &    p_model/='aniso_ntdep')                                       &
     &  CALL nim_stop('perpvisc_model=prpdep requires p_model='//       &
     &                'aniso_tdep or aniso_ntdep')
      IF ((p_model=='aniso_tdep'.OR.p_model=='aniso_ntdep').AND.        &
     &    qpi_model/='std kprp n0'.AND.bhat(1:4)/='full')               &
     &  CALL nim_stop('3D perp conductivity requires bhat=full')

      IF (vmhd_solver=='none') vmhd_solver=solver
      IF (bmhd_solver=='none') bmhd_solver=solver
      IF (temp_solver=='none') temp_solver=solver
!      IF (neoi_flag/='none'.OR.neoe_flag/='none')
!     &  CALL nim_stop('The neoclassical options are not available'
!     &  //' in this version.')
      IF ((neoi_flag/='none'.OR.neoe_flag/='none').AND.                 &
     &    (.NOT.dump_fsa_beq2))                                         &
     &  CALL nim_stop('The neoclassical options require '//             &
     &                'dump_fsa_beq2=TRUE')
      IF (solver=='bl_drect') solver='lapack'
      IF (poly_divv_auto) THEN
        poly_divv=poly_degree
        poly_divv_min=poly_degree
        poly_divv_max=poly_degree
      ENDIF
      IF (poly_divv_max<0) poly_divv_max=poly_divv
      IF (poly_divv_max>=0) poly_divv_min=MAX(poly_divv_min,0_i4)
      IF (poly_divb_max<0) poly_divb_max=poly_divb
      IF (poly_divb_max>=0) poly_divb_min=MAX(poly_divb_min,0_i4)
      IF ((poly_divb_max>=0.OR.poly_divv_max>=0)                        &
     &    .AND.mhdadv_alg=='precor')                                    &
     &  CALL nim_stop('discontinuous fields require '//                 &
     &                'mhdadv_alg==centered')
      IF (concur_bt_adv.AND.(hyp_eta>0._r8.OR.hyp_dbd>0._r8).AND.       &
     &    .NOT.split_hypeta)                                            &
     &   CALL nim_stop('unsplit hyper-resistivity is not compatible'//  &
     &                 'with concur_bt_adv=T')
      IF (nsym_bpre_band>0.AND.(hyp_eta>0._r8.OR.hyp_dbd>0._r8).AND.    &
     &    .NOT.split_hypeta)                                            &
     &   CALL nim_stop('unsplit hyper-resistivity is not compatible'//  &
     &                 'with nsym_bpre_band>0')
      IF (magnetic_bc=='reswall' .OR. magnetic_bc=='reswall_rmp') THEN
        IF (poly_distribution/='gll')                                   &
     &   CALL nim_stop('mag_bc=reswall requires poly_distribution=gll')
        IF (surface_int_formula/='lobatto')                             &
     &   CALL nim_stop('mag_bc=reswall requires surf_int_form=lobatto')
        IF (ngr/=2)                                                     &
     &   CALL nim_stop('mag_bc=reswall requires ngr=2')
        IF (geom/='lin'.AND.rwinterp=='cylgll')                         &
     &    CALL nim_stop('rwinterp=cylgll requires geom=lin')
      ENDIF
      IF (magnetic_bc=="reswall_rmp" .OR. magnetic_bc=='rmp') THEN
        IF (.NOT. lagr_edge_rmp) THEN
          WRITE(*,*) "RMP BCs requre lagr_edge_rmp"
          WRITE(*,*) "Setting lagr_edge_rmp"
        ENDIF
        lagr_edge_rmp=.TRUE.
      ENDIF
      IF (zimp/=0) THEN
        IF ((.not.transfer_eq).AND.(reset_file.EQ."none"))              &
     &    CALL nim_stop('zimp>0 requires transfer_eq=.TRUE.')
        IF (dump_eexp) dump_name=TRIM(dump_name)//'ZwE'
      ENDIF
      IF (.NOT.nonlinear.AND.neutral_model=='full')                     &
     &   CALL nim_stop('neutral_model=full is a nonlinear only option.')
      IF (neutral_model/='none') THEN
        IF (.NOT.nonlinear)                                             &
     &    CALL nim_stop('neutral_model is a nonlinear only option.')
        IF (zimp/=0)                                                    &
     &    CALL nim_stop('zimpurity and neutral models are not currently'&
     &                 //' compatible.')
        IF (nd_hypd/=0)                                                 &
     &    CALL nim_stop('neutral_model can not be currently used with'//&
     &                 ' hyperdiffusivity.')
        IF (continuity/='full')                                         &
     &    CALL nim_stop('neutral_model requires continuity=full')
        IF (advect=='none')                                             &
     &    CALL nim_stop('neutral_model requires advect/=none')
        IF (beta<=0)                                                    &
     &    CALL nim_stop('neutral_model requires beta>0')
        IF (ohms=='hall')                                               &
     &    CALL nim_stop('neutral_model requires ohms/=hall')
        IF (split_visc)                                                 &
     &    CALL nim_stop('neutral_model requires split_visc=F')
        IF (maxit_nl>1)                                                 &
     &    CALL nim_Stop('neutral_model *currently* requires maxit_nl=1')
      ENDIF
      IF (ho_pfield.AND.(poly_degree/=phf_pd))                          &
     &  CALL nim_Stop('ho_pfield=T requires poly_degree=phf_pd')
      IF(TRIM(code)=="nimset") THEN
        IF (MAXVAL(wpack)>0) THEN
          WRITE(*,*) 'wpack parameter name has changed.  The new'
          WRITE(*,*) 'name is wxpack to distinguish it from wypack'
          WRITE(*,*) 'I will transfer it into wxpack, but please'
          WRITE(*,*) 'change in the future'
          wxpack=wpack
        ENDIF
        IF (MAXVAL(amp)>0) THEN
          WRITE(*,*) 'amp parameter name has changed.  The new'
          WRITE(*,*) 'name is ampx to distinguish it from ampy'
          WRITE(*,*) 'I will transfer it into ampx, but please'
          WRITE(*,*) 'change in the future'
          ampx=amp
        ENDIF
      ENDIF
!      IF ((dump_q.OR.dump_divpi.OR.dump_eexp.OR.dump_ja).AND.           &
!     &    .NOT.h5dump) THEN
!         CALL nim_stop('dump_{q,divpi,eexp,ja} require h5dump=T')
!      ENDIF
      IF (zperiod<1) CALL nim_stop('illogic zperiod specified')
      IF (gridshape/='rect') periodicity='none'
!-----------------------------------------------------------------------
!     I'm trying to unify profile input between fg_nimeq and nimset
!     while maintaining backwards compatibility
!     adding checks to make sure that the input is not ambiguous
!-----------------------------------------------------------------------
      IF (p_profile=='power') THEN
         IF ((pres_2 == 0.0) .AND. (pres_4 == 0.0)) THEN
           WRITE(*,*) "P profile is set to power "
           WRITE(*,*) "pres_2 and pres_4 are 0.0."
           WRITE(*,*) "Using ppp values"
           pres_2 = ppp(1)
           pres_4 = ppp(2)
         ELSEIF ((ppp(1) == 0.0) .AND. (ppp(2) == 0.0)) THEN
           WRITE(*,*) "P profile is set to power"
           WRITE(*,*) "ppp(1) and ppp(2) are 0.0."
           WRITE(*,*) "Using pres_2 and pres_4 values"
           ppp(1) = pres_2
           ppp(2) = pres_4
         ELSE
            WRITE(*,*) "P profile is set to power, but "
            WRITE(*,*) "ppp(1) and ppp(2) are nonzero and "
            WRITE(*,*) "and pres_2 and pres_4 are nonzero"
            WRITE(*,*) "Please set on pair of values to zero"
            CALL nim_stop("Ambiguous input with p_profile='power'")
         ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     set Fourier mode information.
!-----------------------------------------------------------------------
      IF (nonlinear) THEN
#ifdef HAVE_FFTW3
        IF (nphi<0) nphi=2**lphi  ! Use lphi for input
#else
        nphi=2**lphi  ! Use lphi for input
#endif
      ENDIF
      IF (dealiase < 3) dealiase=0_i4
#ifndef HAVE_FFTW3
      IF (dealiase > 3) dealiase=3_i4
#endif
!-----------------------------------------------------------------------
!     check domain limits.
!-----------------------------------------------------------------------
      IF (xmax<=xmin) CALL nim_stop("xmax must be > xmin.")
      IF (ymax<=ymin) CALL nim_stop("ymax must be > ymin.")
!-----------------------------------------------------------------------
!     Handle what to do if poly_distribution is not uniform
!-----------------------------------------------------------------------
      IF(poly_distribution=='gll') THEN
        dump_name=TRIM(dump_name)//'gll'
        CALL poly_set(poly_distribution)
      ELSEIF (poly_distribution/='uniform') THEN
        CALL nim_stop('poly_distribution value not recognized')
      ENDIF
!-----------------------------------------------------------------------
!     nimset now uses CASE statements to thread its way through the
!     options.  This is a method for simplifying that logic
!-----------------------------------------------------------------------
      IF (trim(lamprof)=='qspec'.or.lamprof(1:5)=='gmode' .or.          &
     &    trim(lamprof)=='itg' .or. trim(lamprof)=='itg2' .or.         &
     &    trim(lamprof)=='pitprs') THEN
      ELSE
       IF(lam0==0.AND.(TRIM(lamprof)/='zero'                             &
     &    .AND.TRIM(lamprof)/='zerob0')) THEN
        WRITE(*,*) "Because you have lam0=0, changing lamprof to 'zero'"
        lamprof='zero'
       ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     Experimental grids are often one-off development.  Rather than
!-----------------------------------------------------------------------
      IF (gridshape(1:2)=='ss'                                          &
     &         .OR.gridshape(1:3)=='msp'                                &
     &         .OR.gridshape(1:3)=='ldx'                                &
     &         .OR.gridshape(1:3)=='frc'                                &
     &         .OR.gridshape(1:4)=='besl'                               &
     &         .OR.gridshape(1:4)=='wedg'                               &
     &         .OR.gridshape(1:4)=='gwed'                               &
     &         .OR.gridshape(1:4)=='chal'                               &
     &         .OR.gridshape(1:4)=='Xcha'                               &
     &         .OR.gridshape(1:4)=='CTXs'                               &
     &         .OR.gridshape(1:4)=='sort'                               &
     &         .OR.gridshape(1:4)=='gene'                               &
     &         .OR.gridshape(1:4)=='extm'                               &
     &         .OR.gridshape(1:4)=='bell') experiment_grid=.TRUE.
!-----------------------------------------------------------------------
!     Sorting out the right flags.  To make it give the same behavior
!     regardless of the status of transfer_eq.  This makes the logic in
!     sources.f much easier.
!-----------------------------------------------------------------------
      IF (heat_ion_source=='diffusive source') THEN
        IF (rt_transfer_eq) THEN
          heat_ion_source='diffusive-zero_ss'
        ELSE
          heat_ion_source='none'
        ENDIF
      ENDIF
      IF (heat_ion_source=='zero ss') THEN
        IF (rt_transfer_eq) THEN
          heat_ion_source='none'
        ELSE
          heat_ion_source='diffusive-zero_ss'
        ENDIF
      ENDIF
      IF (heat_ele_source=='diffusive source') THEN
        IF (rt_transfer_eq) THEN
          heat_ele_source='diffusive-zero_ss'
        ELSE
          heat_ele_source='none'
        ENDIF
      ENDIF
      IF (heat_ele_source=='zero ss') THEN
        IF (rt_transfer_eq) THEN
          heat_ele_source='none'
        ELSE
          heat_ele_source='diffusive-zero_ss'
        ENDIF
      ENDIF
      IF (electric_source=='diffusive source') THEN
        IF (rt_transfer_eq) THEN
          electric_source='diffusive-zero_ss'
        ELSE
          electric_source='none'
        ENDIF
      ENDIF
      IF (electric_source=='zero ss') THEN
        IF (rt_transfer_eq) THEN
          electric_source='none'
        ELSE
          electric_source='diffusive-zero_ss'
        ENDIF
      ENDIF
      IF (momentum_source=='diffusive source') THEN
        IF (rt_transfer_eq) THEN
          momentum_source='diffusive-zero_ss'
        ELSE
          momentum_source='none'
        ENDIF
      ENDIF
      IF (momentum_source=='zero ss') THEN
        IF (rt_transfer_eq) THEN
          momentum_source='none'
        ELSE
          momentum_source='diffusive-zero_ss'
        ENDIF
      ENDIF
      IF (ndensity_source=='diffusive source') THEN
        IF (rt_transfer_eq) THEN
          ndensity_source='diffusive-zero_ss'
          CALL nim_stop('ndensity does not transfer so'                 &
     &  //' ndensity_source="diffusive source" is not valid.')
        ELSE
          ndensity_source='none'
        ENDIF
      ENDIF
      IF (ndensity_source=='zero ss') THEN
        IF (rt_transfer_eq) THEN
          CALL nim_stop('ndensity does not transfer so'                 &
     &  //' ndensity_source="zero ss" is not valid.')
        ELSE
          ndensity_source='diffusive-zero_ss'
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     initialize adaptive divb diffusion coefficient
!-----------------------------------------------------------------------
      IF (divbd>=0) THEN
        divbdc=divbd
      ELSE
        CALL get_k2min(k2min)
        divbdc=-divbd/dtm/k2min
        divbd=divbdc
      ENDIF
!-----------------------------------------------------------------------
!     Set hv_coeff for linear geometry and if nmodes_total=1
!-----------------------------------------------------------------------
      IF (nonlinear .AND. nphi<2) THEN
        hv_coeff = 0.
      ELSEIF (.NOT.(nonlinear) .AND. lin_nmodes==1) THEN
        hv_coeff = 0.
      ENDIF
      IF (geom=='lin') hv_coeff=0.
!-----------------------------------------------------------------------
!     logic for hot particle pressure anisotropy.
!     for continuum kinetics, phot has same spatial representation as
!     other fields, i.e., ph_pd=poly_degree.
!-----------------------------------------------------------------------
      IF (nF(3)>0) THEN
!        ph_pd=poly_degree
        IF (dump_read_F(3)) phdump=.true.
      ENDIF
!CCK      IF (nm>0) phdump=.true.
!CCK      IF (.NOT.phdump) phqty=0
!-----------------------------------------------------------------------
!     Sovinec and Kruger independently came up with different boundary
!     conditions.  The variables chosen are standardized to match the
!     sources so that it is easier to remember, but these are needed to
!     make input files backward compatible
!-----------------------------------------------------------------------
      IF (nd_bc/="") ndensity_bc=nd_bc
      ! Standardize names.  Homogeneous more accurate description
      IF (ndensity_bc=="dirichlet") ndensity_bc="homogeneous"
      IF (.NOT. zero_bnorm) magnetic_bc="no enforce homo"
      IF (insulate) temp_ion_bc="no enforce homo"
      IF (insulate) temp_ele_bc="no enforce homo"
      IF (flow_bc/="") velocity_bc=flow_bc
      IF (velocity_bc=="dirichlet") velocity_bc="homogeneous" ! nimdevel
      IF (magnetic_dirichlet/="") magnetic_bc=magnetic_dirichlet
      IF (velocity_dirichlet/="") velocity_bc=velocity_dirichlet
      IF (ndensity_dirichlet/="") ndensity_bc=ndensity_dirichlet
      IF (t_RMP>0 .and. magnetic_bc/="reswall_rmp") magnetic_bc='rmp'
      IF (vwall>0 .and. magnetic_bc=="homogeneous") THEN
         magnetic_bc='cylreswall'
      ENDIF
      IF(divbeqd/=0) THEN
        WRITE(nim_wr,*) 'Setting diffuse_eq options from divbeqd'
        diffuse_eq_var='bxxx'
        !diffuse_eq_divfd   =(/divbeqd,0.,0.,0./)
        diffuse_eq_divfd   =(/0.0,0.,0.,0./)
        diffuse_eq_diff    =(/0.,0.,0.,0./)
        diffuse_eq_steps   =(/1,0,0,0/)
      ENDIF
      IF (nsym_tpre_band>0.AND.                                         &
          (temp_ion_bc=='no enforce homo'.AND..NOT.separate_pe          &
          .OR.temp_ele_bc=='no enforce homo'.AND.separate_pe))          &
     &   CALL nim_stop('temperature toroidal preconditioning not '//    &
     &     'implemented for temp_ele_bc==no enforce homo')
!-----------------------------------------------------------------------
!     For GMRES (mhdadv_alg='centered'), every new basis vector used
!       adds to the memory consumption. Depending on the simulation
!       parameters and compute system used this could cause problems.
!-----------------------------------------------------------------------
      IF (node==0.AND.TRIM(ofname)=='nimrod.out') THEN
        IF (mhdadv_alg=='centered') THEN
          IF (maxit>500) THEN
            WRITE(nim_wr,*)                                             &
     &              'WARNING: High maxit requires large memory usage.'
            WRITE(nim_wr,*)'         Proceed with caution!'
            WRITE(out_unit,*)                                           &
     &              'WARNING: High maxit requires large memory usage.'
            WRITE(out_unit,*)'         Proceed with caution!'
          ENDIF
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     set h5io_block_proc_stride to nprocs_layer if -1
!-----------------------------------------------------------------------
      IF (h5io_block_proc_stride==-1_i4) THEN
        h5io_block_proc_stride=nprocs/nlayers
        IF (h5io_block_proc_stride==0) h5io_block_proc_stride=1
      ENDIF
!-----------------------------------------------------------------------
!     set appropriate default solver inputs
!-----------------------------------------------------------------------
      tiopts=0_i4
      tdopts=0_i4
      IF (solver(1:5)=='slu_d') THEN
        ! 1 - report 0/1
        ! 2 - Equil 0/1
        ! 3 - RowPerm 0-NO, 1-LargeDiag
        ! 4 - ColPerm 0-NATURAL, 1-MMD_ATA, 2-MMD_AT_PLUS_A, 3-COLAMD
        !             4-METIS_AT_PLUS_A
        ! 5 - printStat 0/1
        ! 6 - 3d factorization Zgrid size
        ! 7 - write the matrix if 1
        tiopts=(/0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        tdopts=(/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
      ELSEIF (solver(1:5)=='mumps') THEN
        tiopts=(/0, 1, 4, 8, 1, 2, 20, 0, 0, 0, 0, 0, 0, 0, 0/)
        tdopts=(/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
      ELSEIF (solver(1:5)=='hypre') THEN
        tiopts=(/21, 1, 8, 3, 6, 25, 0, 1, 1, 0, 0, 0, 1, 1, 0/)
        !tdopts=(/0.6_r8,0.1_r8,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
        tdopts=(/0.6,0.1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
      ELSEIF (solver(1:7)=='pardiso') THEN
        tiopts=(/3, 13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        tdopts=(/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
      ENDIF
      DO ii=1,solve_nopts
        IF (iopts(ii)==-1)     iopts(ii)=tiopts(ii)
        IF (v_iopts(ii)==-1) v_iopts(ii)=iopts(ii)
        IF (n_iopts(ii)==-1) n_iopts(ii)=iopts(ii)
        IF (t_iopts(ii)==-1) t_iopts(ii)=iopts(ii)
        IF (b_iopts(ii)==-1) b_iopts(ii)=iopts(ii)
        IF (dopts(ii)==-1)     dopts(ii)=tdopts(ii)
        IF (v_dopts(ii)==-1) v_dopts(ii)=dopts(ii)
        IF (n_dopts(ii)==-1) n_dopts(ii)=dopts(ii)
        IF (t_dopts(ii)==-1) t_dopts(ii)=dopts(ii)
        IF (b_dopts(ii)==-1) b_dopts(ii)=dopts(ii)
      ENDDO
!-----------------------------------------------------------------------
!     read in loop voltage
!-----------------------------------------------------------------------
      IF (loop_volt_file) THEN
        OPEN(19,FILE='loop_volt_in',STATUS='OLD',POSITION='REWIND')
        READ(19,*) loop_volt_ncol, loop_volt_nrow
        DO loop_volt_i=1, loop_volt_nrow
          READ(19,*) loop_volt_time(loop_volt_i),                       &
     &               loop_volt_volt(loop_volt_i)
        END DO
      ENDIF
!-----------------------------------------------------------------------
!     check that clean_bnd_divb is being run linearly, with mhd ohms
!     law, fixed resistivity profile, and not hyper resistivity or divb
!     cleaning
!-----------------------------------------------------------------------
      IF (clean_bnd_divb) THEN
        IF (nonlinear) THEN
          CALL nim_stop('clean_bnd_divb requries nonlinear=.FALSE.')
        ENDIF
        IF (ohms/='mhd') THEN
          CALL nim_stop('clean_bnd_divb requries mhd ohms law')
        ENDIF
        IF (eta_model/='fixed') THEN
          CALL nim_stop('clean_bnd_divb requries fixed eta')
        ENDIF
        IF (hyp_eta>0._r8) THEN
          CALL nim_stop('clean_bnd_divb requries hyp_eta=0')
        ENDIF
        IF (hyp_dbd>0._r8) THEN
          CALL nim_stop('clean_bnd_divb requries hyp_dbd=0')
        ENDIF
      ENDIF
      IF ((exn0_divb).AND.(.NOT.transfer_eq)) THEN
        CALL nim_stop('exn0_divb only for transfer_eq=T')
      ENDIF
!-----------------------------------------------------------------------
!     echo the input parameters to the output file.
!-----------------------------------------------------------------------
      IF (echo_in) THEN
        WRITE(out_unit,'(a,/)') 'grid_input:'
        WRITE(UNIT=out_unit,NML=grid_input)
        IF (set_phys_constants) THEN
          WRITE(out_unit,'(/,a,/)') 'const_input:'
          WRITE(UNIT=out_unit,NML=const_input)
        ENDIF
        WRITE(out_unit,'(/,a,/)') 'equil_input:'
        WRITE(UNIT=out_unit,NML=equil_input)
        WRITE(out_unit,'(/,a,/)') 'init_input:'
        WRITE(UNIT=out_unit,NML=init_input)
        WRITE(out_unit,'(/,a,/)') 'physics_input:'
        WRITE(UNIT=out_unit,NML=physics_input)
        WRITE(out_unit,'(/,a,/)') 'closure_input:'
        WRITE(UNIT=out_unit,NML=closure_input)
        WRITE(out_unit,'(/,a,/)') 'numerical_input:'
        WRITE(UNIT=out_unit,NML=numerical_input)
        WRITE(out_unit,'(/,a,/)') 'solver_input:'
        WRITE(UNIT=out_unit,NML=solver_input)
        WRITE(out_unit,'(/,a,/)') 'output_input:'
        WRITE(UNIT=out_unit,NML=output_input)
        WRITE(out_unit,'(/,a,/)') 'extra_input:'
        WRITE(UNIT=out_unit,NML=extra_input)
        IF (nm.ne.0) THEN
          WRITE(out_unit,'(/,a,/)')'particle_input:'
          WRITE(UNIT=out_unit,NML=particle_input)
        ENDIF
        WRITE(out_unit,'(//)')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      END SUBROUTINE read_namelist
!-----------------------------------------------------------------------
!     subprogram 3. rmcomment
!     routine strips out comments beginning with an exclamation point
!-----------------------------------------------------------------------
      SUBROUTINE rmcoment(fileold,filenew)

      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: fileold,filenew
      CHARACTER(128) :: line
      INTEGER, PARAMETER :: nold=55,nnew=56
      INTEGER cmax, ios
      LOGICAL :: file_exist
!-----------------------------------------------------------------------
!     Open files, but make sure the old one exists first.
!-----------------------------------------------------------------------
      INQUIRE(FILE=fileold,EXIST=file_exist)
      IF(.NOT. file_exist) THEN
        CALL nim_stop('The file "',fileold,'" could not be found.')
      ENDIF

      OPEN(UNIT=nold,FILE=fileold,status="OLD")
      OPEN(UNIT=nnew,FILE=filenew,status='REPLACE')
!-----------------------------------------------------------------------
!     Strip comments. Note: line lengths limited to 127 characters
!-----------------------------------------------------------------------
      DO
        READ(UNIT=nold,FMT='(a)',IOSTAT=ios) line
        IF (ios /= 0) EXIT
        cmax=1
        DO WHILE(line(cmax:cmax)/='!' .AND. cmax <= 127)
           cmax=cmax+1
        ENDDO
        IF(cmax > 1) WRITE(nnew,'(a)') TRIM(line(1:cmax-1))
      ENDDO
!-----------------------------------------------------------------------
!     Close files and exit
!-----------------------------------------------------------------------
      CLOSE(nold)
      CLOSE(nnew)

      RETURN
      END SUBROUTINE rmcoment
!-----------------------------------------------------------------------
!     subprogram 4. get_k2min.
!     estimate max of k**2 based on mesh size
!-----------------------------------------------------------------------
      SUBROUTINE get_k2min(k2min)
      USE local
      USE input
      IMPLICIT NONE

      REAL(r8) :: k2min
!-----------------------------------------------------------------------
!     only rectangular initialization currently implemented
!     -DCB: 10/01/04
!-----------------------------------------------------------------------
      IF (gridshape=='rect') THEN
        k2min=(1/(xmax-xmin))**2+(1/(ymax-ymin))**2
      ELSE
        CALL nim_stop                                                   &
     &    ('get_k2min not yet implemented for this gridshape')
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE get_k2min
!-----------------------------------------------------------------------
!     subprogram 5. broadcast_input
!     broadcast info read by proc 0 out of nimrod.in (namelist reads)
!     to all processors. every quantity in input module is broadcast in
!     case it was read.
!     IMPORTANT: anytime a variable gets added to input module, it MUST
!     added to this routine as well
!-----------------------------------------------------------------------
      SUBROUTINE broadcast_input(lnode)
      USE input
      USE mpi_nim
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: lnode
      INTEGER(i4) :: ierr,ii,cw

      cw=mpi_comm_world
!-----------------------------------------------------------------------
!     physics specifications
!-----------------------------------------------------------------------
      CALL mpi_bcast(nonlinear,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(lin_nmax,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(ndens,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nedge,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(npp,10,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(n_profile,8,lnode)
      CALL bcast_str(eta_model,16,lnode)
      CALL mpi_bcast(elecd,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(eta_ref_t,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(elecd_max,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(elecd_min,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(elecd_chodura,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(f_chodura,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(gravity,3,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(gamma_heat,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(t_RMP,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(t_rmp_offset,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(rmp_mult_fac,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(rmp_omega,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(RMP_type,12,lnode)
      CALL bcast_str(rmp_file_type,12,lnode)
      CALL mpi_bcast(lagr_edge_rmp,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(load_n0_rmp,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(b_nonsym,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(reverse_Ip,1,mpi_nim_logical,0,cw,ierr)
      CALL bcast_str(heat_ion_source,20,lnode)
      CALL bcast_str(heat_ele_source,20,lnode)
      CALL bcast_str(electric_source,20,lnode)
      CALL bcast_str(momentum_source,20,lnode)
      CALL bcast_str(ndensity_source,20,lnode)
      CALL mpi_bcast(bump_a,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nu_relax,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(ndensity_bc,16,lnode)
      CALL bcast_str(magnetic_bc,16,lnode)
      CALL bcast_str(temp_ion_bc,16,lnode)
      CALL bcast_str(temp_ele_bc,16,lnode)
      CALL bcast_str(velocity_bc,16,lnode)
      CALL mpi_bcast(separate_pe,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(kin_visc,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(iso_visc,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(par_visc,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(gyr_visc,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(hv_filter,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(hv_time,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(hv_coeff,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(hv_power,1,mpi_nim_int,0,cw,ierr)
      DO ii=1,4
        CALL mpi_bcast(dvac(ii),1,mpi_nim_real,0,cw,ierr)
        CALL mpi_bcast(dexp(ii),1,mpi_nim_real,0,cw,ierr)
        CALL bcast_str(ds_function(ii),8,lnode)
      ENDDO
      CALL bcast_str(ds_use,8,lnode)
      CALL mpi_bcast(ds_nqty,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(xvac,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(beta,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(neut_dens_profile,8,lnode)
      CALL mpi_bcast(neutral_dens,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(eq_neut_flow,8,lnode)
      CALL mpi_bcast(vn0,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(neut_temp_profile,8,lnode)
      CALL mpi_bcast(neutral_temp,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(elecd0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nd_diff0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(iso_visc0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(kin_visc0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_perpi0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_perpe0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_perp0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(use_fp,1,mpi_nim_logical,0,cw,ierr)

!-----------------------------------------------------------------------
!     loop voltage specifications
!-----------------------------------------------------------------------
      CALL mpi_bcast(loop_volt_file,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(loop_volt_nrow,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(loop_volt_time,loop_volt_nrow,mpi_nim_real,0,      &
     &     cw,ierr)
      CALL mpi_bcast(loop_volt_volt,loop_volt_nrow,mpi_nim_real,0,      &
     &     cw,ierr)
      CALL mpi_bcast(eq_loop_voltage,1,mpi_nim_real,0,cw,ierr)

!-----------------------------------------------------------------------
!     constants specifications
!-----------------------------------------------------------------------
      CALL mpi_bcast(set_phys_constants,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(chrg_input,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(zeff_input,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(mi_input,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(me_input,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(gam_input,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(kblz_input,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(mu0_input,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(c_input,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(mc2_kT_input,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(mc2_kTb_input,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(pe_frac,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(advect,8,lnode)
      CALL bcast_str(continuity,16,lnode)
      CALL mpi_bcast(nd_diff,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nd_hypd,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nd_floor,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nd_exp,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nd_dart_fac,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nd_dart_upw,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(t_dart_upw,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(ohms,8,lnode)
      CALL mpi_bcast(n_upw_aniso,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(n_upw_limit,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(t_upw_aniso,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(t_upw_limit,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nd_floor_upw,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nd_width_upw,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(t_floor_upw,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(t_width_upw,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nd_nodal_floor,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ti_nodal_floor,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(te_nodal_floor,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nd_correrr,1,mpi_nim_logical,0,mpi_comm_world,ierr)
      CALL bcast_str(p_computation,16,lnode)
      CALL mpi_bcast(no_pe_hall,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(be0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(thetab,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(phib,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(theta_init,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(phi_init,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(lamprof,6,lnode)
      CALL bcast_str(f_profile,16,lnode)
      CALL bcast_str(fmod_profile,16,lnode)
      CALL bcast_str(pmod_profile,16,lnode)
      CALL bcast_str(p_profile,16,lnode)
      CALL bcast_str(pe_profile,8,lnode)
      CALL mpi_bcast(lam0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(rbreak,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(alpha,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(eq_flow,16,lnode)
      CALL mpi_bcast(dfac,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(pit_0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(pit_2,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(pit_4,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(pres_2,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(pres_4,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(f_cent,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(f_edge,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fpp,10,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fmod_amp,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fmod_center,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fmod_width,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(pmod_amp,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(pmod_center,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(pmod_width,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ppp,10,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(tepp,10,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(te_edge,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(te_cent,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ve0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(thetav,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(phiv,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(eqflow_width,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(loop_volt,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(tloopv0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(tloopv1,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(i_desired,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(loop_rate,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(loop_rate2,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(e_vertical,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(t_e_vert0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(t_e_vert1,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(t_e_vert2,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(gun_max,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(gun_lim0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(gun_lim1,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(gun_func,8,lnode)
      CALL mpi_bcast(mag_tau,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(pert_m,20,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(pert_n,20,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(pert_cos,20,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(pert_sin,20,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(pert_width,20,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(pert_phase,20,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(toffset,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(tdown,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(tperiod,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(lam0rf,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(rrf,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(zrf,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(phirf,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(rf_power_ec,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(rf_freq_ec,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(discrete_extrap,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(mynphi,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(more_phiplanes,1,mpi_nim_logical,0,cw,ierr)
      CALL bcast_str(shepard_alg,20,lnode)
      CALL mpi_bcast(delta_pol,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(delta_tor,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(tmodperiod,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(tmodphase,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(t_on,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(t_off,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(rf_extrap,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(pcsstatus,2,mpi_nim_int,0,cw,ierr)
      CALL bcast_str(rf_modulation,20,lnode)
      CALL bcast_str(control,20,lnode)
      CALL bcast_str(init_type,16,lnode)
      CALL mpi_bcast(bamp,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nx,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(ny,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(tor_eqja_fe,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(ideal_like,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(mvac_nim,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(pressure_gauge_nim,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(set_vac_quadpts,1,mpi_nim_logical,0,cw,ierr)
      DO ii=1,10
        CALL bcast_str(chrdbug(ii),10,lnode)
      ENDDO
      CALL mpi_bcast(logdbug,10,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(intdbug,10,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(dbldbug,10,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(neutral_model,16,lnode)
      CALL mpi_bcast(ion_eff,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(sion_fac,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(srec_fac,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(scx_fac,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ndn_diff,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(vn_iso_visc,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(vn_kin_visc,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(neut_vel_bc,16,lnode)
      CALL bcast_str(neutral_rates,16,lnode)
!-----------------------------------------------------------------------
!     closure specifications
!-----------------------------------------------------------------------
      CALL mpi_bcast(close_nprocs,1,mpi_nim_int,0,cw,ierr)
      CALL bcast_str(p_model,16,lnode)
      CALL bcast_str(closure_model,12,lnode)
      CALL bcast_str(bhat,4,lnode)
      CALL mpi_bcast(closure_dump,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(cel_flow,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(cel_temp,1,mpi_nim_logical,0,cw,ierr)
      CALL bcast_str(analytic_rhs,6,lnode)
      CALL bcast_str(eqn_model,12,lnode)
      CALL bcast_str(kin_model,17,lnode)
      CALL bcast_str(F_init_type,16,lnode)
      CALL bcast_str(Fele_bc,16,lnode)
      CALL bcast_str(Fion_bc,16,lnode)
      CALL bcast_str(Fhot_bc,16,lnode)
      CALL mpi_bcast(dump_read_F,3,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(F_mass_op_flag,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(kinetic_eq,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(nF,3,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(pd_xi,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(nq_xi,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(C_nq_xi,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(m_xi,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(m_xit,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(m_xip,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(tpb_min,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(gridshape_v,9,lnode)
      CALL bcast_str(pa_basis,3,lnode)
      CALL mpi_bcast(ns,3,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(ns2,3,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(smin,3,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(smid,3,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(smax,3,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(s1_quad(1),8,lnode)
      CALL bcast_str(s1_quad(2),8,lnode)
      CALL bcast_str(s1_quad(3),8,lnode)
      CALL bcast_str(s2_quad(1),8,lnode)
      CALL bcast_str(s2_quad(2),8,lnode)
      CALL bcast_str(s2_quad(3),8,lnode)
      CALL mpi_bcast(sp1,3,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(sp2,3,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ns1_fixed,3,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(ns2_fixed,3,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(s1_fixed,6,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(s2_fixed,6,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(dxi,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(xi0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(tpow,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(tcoe,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(betafrac_sdd,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(single_s,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(s_poly_exp,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(steady_state,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(Deff_e,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(Deff_i,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(Deff_h,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(Deff_epll,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(Deff_ipll,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(Deff_hpll,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(F_dexp,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(F_dvac,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(B_synch,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(tau_c_rel_norm,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(cel_output,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(fcel,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fcol,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fder,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fint,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fadv,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(coll_model,7,lnode)
      CALL mpi_bcast(xi_prec,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(ss_prec,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(ss_rhs,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(ssoff_prec,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(it_calls,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(kpll_mfe,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(kplle_mfe,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(kplli_mfe,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(divq_mfe_ibp,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(TFtheta,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(par_visc_mfe,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_perp,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_pll,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_perpe,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_plle,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_perpi,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_plli,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_perpn,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_pll_max,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_pll_min,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_pll_ref_t,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_prp_ref_b,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_ref_n,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(kprp_mnrat,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(k_cross,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(tequil_rate,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ohm_heat,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(visc_heat,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(zimp,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(re_orbits,1,mpi_nim_logical,0,cw,ierr)
      CALL bcast_str(parvisc_model,16,lnode)
      CALL bcast_str(perpvisc_model,16,lnode)
      CALL bcast_str(qpi_model,12,lnode)
      CALL mpi_bcast(odd_pert,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(odd_pert_str,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(psi_lcfs,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(q_cent,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(q_x,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(q_y,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(psi_per,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(v_cent,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(tdep_coul_log,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(coulomb_logarithm,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(magfac_ele,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(magfac_ion,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(delta0_ele,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(delta1_ele,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(gamma0_ele,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(gamma1_ele,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(delta0_ion,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(delta1_ion,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(gamma0_ion,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(gamma1_ion,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(tdep_tequil,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(qpi_drive,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(qpe_drive,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(ppi_drive,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(ppe_drive,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(rpe_drive,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(gei_drive,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(dr_drive,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(te_drive,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(ti_drive,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(ve_drive,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(vi_drive,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(old_nF,3,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(old_ns,3,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(old_ns2,3,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(old_ph_pd,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(T_sfac,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(neoe_flag,8,lnode)
      CALL bcast_str(neoi_flag,8,lnode)
      CALL mpi_bcast(mu_e,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(mu_i,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(neo_bump_r0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(neo_bump_amp,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(neo_axis_r,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(neo_axis_z,1,mpi_nim_real,0,cw,ierr)
!-----------------------------------------------------------------------
!     numerical specifications
!-----------------------------------------------------------------------
      CALL mpi_bcast(dtm,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(dt_initial,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(dt_stop,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(dt_incr,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(v_cfl,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(mhdadv_alg,8,lnode)
      CALL mpi_bcast(concur_bt_adv,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(concur_vF_adv,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(nl_cfl_lim,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(tmax,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(cpu_tmax,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nstep,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(npc,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(fom,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(feta,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(split_resist,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(fvsc,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(split_visc,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(fthc,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fb_vxb,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fv_vdgv,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fp_vdgp,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fn_vdgn,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(divbd,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(divbdc,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fdivb,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ndivb,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(kdivb_2_limit,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(split_divb,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(dt_change_frac,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ave_change_limit,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(n_dt_release,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(n_mat_update,1,mpi_nim_int,0,cw,ierr)
      CALL bcast_str(siop_type,8,lnode)
      CALL mpi_bcast(si_fac_mhd,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(mhd_si_iso,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(si_fac_hall,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(si_fac_pres,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(si_fac_j0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(si_fac_nl,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(conform,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(lump_all,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(ngr,1,mpi_nim_int,0,cw,ierr)
      CALL bcast_str(integration_formula,8,lnode)
      CALL bcast_str(surface_int_formula,8,lnode)
      CALL bcast_str(rwinterp,8,lnode)
      CALL mpi_bcast(poly_degree,1,mpi_nim_int,0,cw,ierr)
      CALL bcast_str(poly_distribution,16,lnode)
      CALL mpi_bcast(poly_divb,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(poly_divb_min,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(poly_divb_max,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(disc_dbd,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(hyp_eta,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fhyp_eta,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(split_hypeta,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(hyp_dbd,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fhyp_dbd,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fdivv,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(fpvrt,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ddivv,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(dpvrt,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(poly_divv,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(poly_divv_min,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(poly_divv_max,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(poly_divv_auto,1,mpi_nim_logical,0,cw,ierr)
      CALL bcast_str(met_spl,5,lnode)
      CALL mpi_bcast(r0dr_weight,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(transfer_eq,1,mpi_nim_logical,0,cw,ierr)
      CALL bcast_str(transfer_eq_fields,4,lnode)
      CALL bcast_str(transfer_n0,4,lnode)
      CALL bcast_str(init_type_transfer_n0,16,lnode)
      CALL mpi_bcast(transfer_overwrite,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(rt_transfer_eq,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(scal_byparts,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(normalize,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(clean_bnd_divb,1,mpi_nim_logical,0,cw,ierr)

!-----------------------------------------------------------------------
!     linear solver specifications
!-----------------------------------------------------------------------
      CALL mpi_bcast(tol,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(maxit,1,mpi_nim_int,0,cw,ierr)
      CALL bcast_str(solver,8,lnode)
      CALL bcast_str(vmhd_solver,8,lnode)
      CALL bcast_str(bmhd_solver,8,lnode)
      CALL bcast_str(temp_solver,8,lnode)
      CALL mpi_bcast(off_diag_fac,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(extrap_order,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(nsym_bpre_band,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(nsym_tpre_band,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(nsym_pre_rpass,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(nsym_pre_rfac,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(nsym_pre_rtype,12,lnode)
      CALL mpi_bcast(maxit_nl,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(tol_nl,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(iopts,solve_nopts,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(dopts,solve_nopts,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(v_iopts,solve_nopts,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(v_dopts,solve_nopts,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(b_iopts,solve_nopts,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(b_dopts,solve_nopts,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(n_iopts,solve_nopts,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(n_dopts,solve_nopts,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(t_iopts,solve_nopts,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(t_dopts,solve_nopts,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(solver_thread_reduction,1,mpi_nim_real,0,cw,ierr)

!-----------------------------------------------------------------------
!     grid specifications
!-----------------------------------------------------------------------
      CALL bcast_str(gridshape,8,lnode)
      CALL bcast_str(periodicity,8,lnode)
      CALL mpi_bcast(xmin,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(xmax,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ymin,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ymax,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(xo,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(yo,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(zperiod,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(mx,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(my,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(firstx,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(firsty,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(skew,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nxbl,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(nybl,1,mpi_nim_int,0,cw,ierr)
      CALL bcast_str(geom,3,lnode)
      CALL mpi_bcast(lphi,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(nphi,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(lin_nmodes,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(dealiase,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(per_length,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(decompflag,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(nlayers,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(qpack,15,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(xpack,15,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(wxpack,15,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ampx,15,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ypack,15,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(wypack,15,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ampy,15,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(mm_mext,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(mm_mint,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(mm_imap_use_int,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(mm_mdmap,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(mm_rw_dmap,1,mpi_nim_logical,0,cw,ierr)

!-----------------------------------------------------------------------
!     output specifications
!-----------------------------------------------------------------------
      CALL mpi_bcast(detflag,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(dump_ja,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(dump_eexp,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(dump_q,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(dump_divpi,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(dump_psi_eq,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(dump_fsa_beq2,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(itflag,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(nhist,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(n_probe,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(ihist,10,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(jhist,10,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(hist_binary,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(lin0eq_energy,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(exn0_divb,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(xy_stride,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(xt_stride,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(yt_stride,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(x0fac,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(y0fac,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(xdraw_dir,64,lnode)
      CALL mpi_bcast(ndxout,1,mpi_nim_int,0,cw,ierr)
      CALL bcast_str(dx_dir,64,lnode)
      CALL mpi_bcast(ndump,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(tdump,1,mpi_nim_real,0,cw,ierr)
      CALL bcast_str(dump_file,128,lnode)
      CALL bcast_str(dump_name,64,lnode)
      CALL bcast_str(dump_dir,64,lnode)
      CALL mpi_bcast(dump_over,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(h5dump,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(h5io_block_proc_stride,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(h5io_node_proc_stride,1,mpi_nim_int,0,cw,ierr)
      CALL bcast_str(reset_file,128,lnode)
      CALL bcast_str(reset_ve_file,128,lnode)
      CALL bcast_str(reset_vi_file,128,lnode)
      CALL bcast_str(reset_vh_file,128,lnode)
      CALL mpi_bcast(reset_time,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(reset_step,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(vwall,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(efkphi,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(efrad,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(reswall_nmax,1,mpi_nim_int,0,cw,ierr)
!-----------------------------------------------------------------------
!     particle paramaters
!-----------------------------------------------------------------------
      CALL mpi_bcast(nm,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(vhmx,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(avu0,3,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(avw0,3,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(mass0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(charge0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(ecrit,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(psipmax,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(Ppsi0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(R0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(betafrac,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(sv_input,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(reduce_p,8,mpi_nim_char,0,cw,ierr)
      CALL mpi_bcast(orbits,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(eparallel,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(anisop,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(nh0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(rh0,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(hwdth,1,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(pphparam,16,mpi_nim_real,0,cw,ierr)
      CALL mpi_bcast(restart,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(phqty,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(nbins,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(iseed,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(vdlms,8,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(ph_pd,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(phf_pd,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(bmxs,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(emxs,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(add_phot,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(pdbg,1,mpi_nim_int,0,cw,ierr)
      CALL mpi_bcast(pnorm_field,50,mpi_nim_char,0,cw,ierr)
      CALL mpi_bcast(tpartdr,50,mpi_nim_char,0,cw,ierr)
      CALL mpi_bcast(part_type,50,mpi_nim_char,0,cw,ierr)
      CALL mpi_bcast(phase_use,50,mpi_nim_char,0,cw,ierr)
      CALL mpi_bcast(weight_use,50,mpi_nim_char,0,cw,ierr)
      CALL mpi_bcast(hot_distro,50,mpi_nim_char,0,cw,ierr)
      CALL mpi_bcast(trace,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(ho_pfield,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(phdump,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(dump_part,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(scale4phot,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(model_probes,1,mpi_nim_logical,0,cw,ierr)
      CALL mpi_bcast(find_type,6,mpi_nim_char,0,cw,ierr)
      CALL mpi_bcast(field_type,6,mpi_nim_char,0,cw,ierr)
      CALL mpi_bcast(phot_type,6,mpi_nim_char,0,cw,ierr)
!-----------------------------------------------------------------------
!     Handle global fixing of parameters
!-----------------------------------------------------------------------
      CALL poly_set(poly_distribution)
      RETURN
      END SUBROUTINE broadcast_input
!-----------------------------------------------------------------------
!     subprogram 6. bcast_str
!     broadcast a string from one processor to another by
!     treating them as integer should not have to call these from
!     send_rblock and send_seam if I can eventually get MPI characters
!     data typeto work correctly !!
!
!     TODO: this can be removed
!-----------------------------------------------------------------------
      SUBROUTINE bcast_str(str,n,lnode)
      USE local
      USE mpi_nim
      implicit none
      character*(*) str
      INTEGER(i4),INTENT(IN) :: lnode
      INTEGER(i4) :: n,i,value,ierr

      do i = 1,n
        if (lnode == 0) value = ichar(str(i:i))
        CALL mpi_bcast(value,1,mpi_nim_int,0,mpi_comm_world,ierr)
        str(i:i) = char(value)
      enddo

      return
      END
!-----------------------------------------------------------------------
!     subprogram 7. input_normalization
!     routine for conversion of input parameters to ones used by the
!     code.
!-----------------------------------------------------------------------
      SUBROUTINE input_normalization
      USE physdat
      USE input
      USE pardata
      USE global
      USE io
      IMPLICIT NONE
      REAL(r8) :: ndnorm

!-----------------------------------------------------------------------
!     If n=0 diffusivities are specified set them here for layer 0
!-----------------------------------------------------------------------
      IF (elecd0>0._r8.OR.nd_diff0>0._r8.OR.iso_visc0>0._r8.OR.         &
     &    kin_visc0>0._r8.OR.k_perpi0>0._r8.OR.k_perpe0>0._r8.OR.       &
     &    k_perp0>0._r8) THEN
        IF (nlayers/=nmodes_total.AND.TRIM(ofname)=='nimrod.out')       &
     &  CALL nim_stop('Using n=0 diffusivities requires nlayers=nmodes')
        IF (ilayer==0) THEN
          IF (elecd0>0._r8) elecd=elecd0
          IF (nd_diff0>0._r8) nd_diff=nd_diff0
          IF (iso_visc0>0._r8) iso_visc=iso_visc0
          IF (kin_visc0>0._r8) kin_visc=kin_visc0
          IF (k_perpi0>0._r8) k_perpi=k_perpi0
          IF (k_perpe0>0._r8) k_perpe=k_perpe0
          IF (k_perp0>0._r8) k_perp=k_perp0
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     Set density normalization
!-----------------------------------------------------------------------
      IF (k_ref_n>0._r8) THEN
        ndnorm=k_ref_n
      ELSE
        ndnorm=ndens
      ENDIF
!-----------------------------------------------------------------------
!     To get the closures to match up properly independent of form,
!     we change the thermal diffusivities to thermal conductivities.
!     We do it here since it has to be after zeff is defined.
!     This is done differently than nimuw because of the closures issue.
!     Multiply conductivies by gamma-1 factor in T evolution equation;
!     nimrod.in now contains the correct thermal diffusivities.
!-----------------------------------------------------------------------
      IF (.NOT.ideal_like) THEN
        IF (qpi_model(1:2)=="st") THEN
          k_perpi   = gamm1*k_perpi   *ndnorm/zeff
          k_perpe   = gamm1*k_perpe   *ndnorm
          k_perpn   = gamm1*k_perp    *neutral_dens
          k_plli    = gamm1*k_plli    *ndnorm/zeff
          k_plle    = gamm1*k_plle    *ndnorm
!-----------------------------------------------------------------------
!         These aren't really used, but multiply by ndnorm anyway
!-----------------------------------------------------------------------
          k_pll     = gamm1*k_pll     *ndnorm
          k_perp    = gamm1*k_perp    *ndnorm
!-----------------------------------------------------------------------
!         This will result in slightly different results from nimuw
!         if zeff!=1 and separate_pe=true since the temperature-dependent
!         k_pll_e will have different min/max
!-----------------------------------------------------------------------
          k_pll_min = gamm1*k_pll_min *ndnorm/zeff
          k_pll_max = gamm1*k_pll_max *ndnorm/zeff
        ENDIF
!-----------------------------------------------------------------------
!       viscosity is especially tricky:
!       specify flow diffusivity, convert to flow conductivity by
!       multiplying by n_i, finally end up with momentum conductivity by
!       multiplying by mtot
!-----------------------------------------------------------------------
        IF (qpi_model(1:2)=="st".OR.perpvisc_model=='fixed') THEN
          iso_visc  = iso_visc *ndnorm/zeff *mtot
          kin_visc  = kin_visc *ndnorm/zeff *mtot
        ELSE
          iso_visc  = iso_visc *mtot
          kin_visc  = kin_visc *mtot
        ENDIF
        IF (qpi_model(1:2)=="st".OR.parvisc_model=='fixed') THEN
          par_visc  = par_visc *ndnorm/zeff *mtot
        ELSE
          par_visc  = par_visc *mtot
        ENDIF
        vn_iso_visc  = vn_iso_visc *neutral_dens *(ms(1)+ms(2))
        vn_kin_visc  = vn_kin_visc *neutral_dens *(ms(1)+ms(2))
!-----------------------------------------------------------------------
!       modify the coefficients with the reference density/mag field.
!-----------------------------------------------------------------------
        IF (qpi_model(1:2)=="st".AND..NOT.normalize) THEN
          IF (p_model=='aniso_ntdep'.AND.k_ref_n>0._r8) THEN
            k_perpi = k_perpi/k_ref_n**2
            k_perpe = k_perpe/k_ref_n**2
            k_perp  = k_perp /k_ref_n**2
          ENDIF
          IF (perpvisc_model=='prpdep'.AND.k_ref_n>0._r8) THEN
            kin_visc  = kin_visc /k_ref_n**2
            iso_visc  = iso_visc /k_ref_n**2
          ENDIF
          IF (p_model=='aniso_tdep'.OR.p_model=='aniso_ntdep'.AND.      &
     &        k_prp_ref_b>0._r8) THEN
            k_perpi = k_perpi*k_prp_ref_b**2
            k_perpe = k_perpe*k_prp_ref_b**2
            k_perp  = k_perp *k_prp_ref_b**2
          ENDIF
          IF (perpvisc_model=='prpdep'.AND.k_prp_ref_b>0._r8) THEN
            kin_visc  = kin_visc *k_prp_ref_b**2
            iso_visc  = iso_visc *k_prp_ref_b**2
          ENDIF
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE input_normalization
!-----------------------------------------------------------------------
!     subprogram 8. direct_check.
!     test if the preconditioner option is a direct solve over the
!     fe mesh.
!-----------------------------------------------------------------------
      LOGICAL PURE FUNCTION direct_check(precon)
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: precon

      IF (precon=='lapack'    .OR.precon=='seq_slu'.OR.                 &
     &    precon(1:5)=='slu_d'.OR.precon(1:6)=='pastix'.OR.             &
     &    precon(1:5)=='mumps'.OR.precon(1:5)=='hypre'.OR.              &
     &    precon(1:7)=='pardiso') THEN
        direct_check=.true.
      ELSE
        direct_check=.false.
      ENDIF
!-----------------------------------------------------------------------
!     terminate routine
!-----------------------------------------------------------------------
      RETURN
      END FUNCTION direct_check
!-----------------------------------------------------------------------
!     subprogram 9. perturbation_set
!     open and read the namelist input.
!-----------------------------------------------------------------------
      SUBROUTINE perturbation_set(nmodes,keff)
      USE local
      USE input
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nmodes
      REAL(r8), DIMENSION(nmodes), INTENT(IN) :: keff

      INTEGER(i4) :: im, ipert, mpert, jpert, jpn, jnn
      INTEGER(i4), DIMENSION(nmodes) :: num_m, npert
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: ipn,iipn, inn,iinn
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: num_pn, num_nn
      LOGICAL :: found_m
!-----------------------------------------------------------------------
!     Each mode on a processor will have a perturbation set with it
!-----------------------------------------------------------------------
      num_m=0
      ALLOCATE(pert(nmodes))
      ALLOCATE(ipn(nmodes,SIZE(pert_m))); ipn=0
      ALLOCATE(iipn(nmodes,SIZE(pert_m))); iipn=0
      ALLOCATE(inn(nmodes,SIZE(pert_m))); inn=0
      ALLOCATE(iinn(nmodes,SIZE(pert_m))); iinn=0
      ALLOCATE(num_pn(nmodes,SIZE(pert_m))); num_pn=0
      ALLOCATE(num_nn(nmodes,SIZE(pert_m))); num_nn=0
!-----------------------------------------------------------------------
!     First initialize the basic parts of the datastructure
!-----------------------------------------------------------------------
      IF(geom=='lin') THEN
         npert=NINT(keff*(per_length/2./3.1415926535897932385_r8))
      ELSE
         npert=keff
      ENDIF
      DO ipert=1,nmodes
          pert(ipert)%n=npert(ipert)
          pert(ipert)%nuse=0
      ENDDO
!-----------------------------------------------------------------------
!     Map the pert_* arrays into the datastructure.  Auxilary arrays
!     handle the sorting.  The negative n's are a pain.
!-----------------------------------------------------------------------
      DO ipert=1,SIZE(pert_m)
        IF(pert_n(ipert)>npert(nmodes)) CYCLE
        DO im=1,nmodes
           pert(im)%n=npert(im)
           pert(im)%nuse=0
           IF(npert(im)==pert_n(ipert)) THEN
              num_m(im)=num_m(im)+1
              num_pn(im,num_m(im))=num_pn(im,num_m(im))+1
              ipn(im,num_m(im))=ipert
              iipn(im,ipert)=num_m(im)
           ENDIF
        ENDDO
      ENDDO
      ! Now go and handle the negative modes
      DO ipert=1,SIZE(pert_m)
        IF(pert_n(ipert)>npert(nmodes)) CYCLE
        im_loop: DO im=1,nmodes
           IF(npert(im)==-pert_n(ipert).AND.keff(im)>0.) THEN
              ! See if we find a matching m
              found_m=.false.
              DO jpert=1,SIZE(pert_m)
                 IF(ipert==jpert) CYCLE
                 IF(pert_m(jpert)==pert_m(ipert)) THEN
                    found_m=.true.
                    inn(im,iipn(im,jpert))=ipert
                    iinn(im,ipert)=iipn(im,jpert)
                  num_nn(im,iipn(im,jpert))=num_nn(im,iipn(im,jpert))+1
                    EXIT im_loop
                 ENDIF
              ENDDO
              IF(.NOT. found_m) THEN
                 num_m(im)=num_m(im)+1
                 inn(im,num_m(im))=ipert;      iinn(im,ipert)=num_m(im)
                 num_nn(im,num_m(im))=num_pn(im,num_m(im))+1
                 EXIT im_loop
              ENDIF
           ENDIF
        ENDDO im_loop
      ENDDO
      DO im=1,nmodes
         ALLOCATE(pert(im)%m(num_m(im)))
         ALLOCATE(pert(im)%Pcos(num_m(im))); pert(im)%Pcos=0
         ALLOCATE(pert(im)%Psin(num_m(im))); pert(im)%Psin=0
         ALLOCATE(pert(im)%Rcos(num_m(im))); pert(im)%Rcos=0
         ALLOCATE(pert(im)%Rsin(num_m(im))); pert(im)%Rsin=0
         ALLOCATE(pert(im)%Arr(num_m(im))); pert(im)%Arr=0
         ALLOCATE(pert(im)%Ari(num_m(im))); pert(im)%Ari=0
         ALLOCATE(pert(im)%Air(num_m(im))); pert(im)%Air=0
         ALLOCATE(pert(im)%Aii(num_m(im))); pert(im)%Aii=0
         ALLOCATE(pert(im)%width(num_m(im))); pert(im)%width=0.
         ALLOCATE(pert(im)%ptype(num_m(im)))
         DO ipert=1,num_m(im)
           pert(im)%m(ipert)=pert_m(ipn(im,ipert))
           DO jpn=1,num_pn(im,ipert)
             pert(im)%Pcos(ipert)= pert_cos(ipn(im,ipert))
             pert(im)%Psin(ipert)= pert_sin(ipn(im,ipert))
             pert(im)%ptype(ipert)=pert_type(ipn(im,ipert))
             pert(im)%width(ipert)=pert_width(ipn(im,ipert))
           ENDDO
           DO jnn=1,num_nn(im,ipert)
              pert(im)%Rcos(ipert)= pert_cos(inn(im,ipert))
              pert(im)%Rsin(ipert)= pert_sin(inn(im,ipert))
              pert(im)%ptype(ipert)=pert_type(inn(im,ipert))
              pert(im)%width(ipert)=pert_width(inn(im,ipert))
           ENDDO
         ENDDO
         pert(im)%Arr(:)=(pert(im)%Pcos(:) + pert(im)%Rcos(:))/2.
         pert(im)%Aii(:)=(pert(im)%Pcos(:) - pert(im)%Rcos(:))/2.
         ! For n=0, no imaginary components
         IF (INT(keff(im))/=0) THEN
           pert(im)%Ari(:)=-(pert(im)%Psin(:) + pert(im)%Rsin(:))/2.
           pert(im)%Air(:)= (pert(im)%Psin(:) - pert(im)%Rsin(:))/2.
         ENDIF
         IF(MAXVAL(ABS(pert(im)%Arr(:))) > 0) pert(im)%nuse=1
         IF(MAXVAL(ABS(pert(im)%Ari(:))) > 0) pert(im)%nuse=1
         IF(MAXVAL(ABS(pert(im)%Air(:))) > 0) pert(im)%nuse=1
         IF(MAXVAL(ABS(pert(im)%Aii(:))) > 0) pert(im)%nuse=1
      ENDDO
      DEALLOCATE(ipn,iipn,inn,iinn,num_pn,num_nn)
!-----------------------------------------------------------------------
!     terminate.
!-----------------------------------------------------------------------
      END SUBROUTINE perturbation_set
