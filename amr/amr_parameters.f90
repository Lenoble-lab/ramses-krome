module amr_parameters

  ! Define real types
  integer,parameter::sp=kind(1.0E0)
#ifndef NPRE
  integer,parameter::dp=kind(1.0E0) ! default
#else
#if NPRE==4
  integer,parameter::dp=kind(1.0E0) ! real*4
#else
  integer,parameter::dp=kind(1.0D0) ! real*8
#endif
#endif

! Precision for uold, unew, fgrav
#ifndef NUPRE
  integer,parameter::udp=kind(1.0D0) ! default
#else
#if NUPRE==4
  integer,parameter::udp=kind(1.0E0) ! real*4
#else
  integer,parameter::udp=kind(1.0D0) ! real*8
#endif
#endif

! Precision for radiative transfer
#ifndef NRTPRE
  integer,parameter::rtdp=kind(1.0D0) ! default
#else
#if NRTPRE==4
  integer,parameter::rtdp=kind(1.0E0) ! real*4
#else
  integer,parameter::rtdp=kind(1.0D0) ! real*8
#endif
#endif

#ifdef QUADHILBERT
  integer,parameter::qdp=kind(1.0_16) ! real*16
#else
  integer,parameter::qdp=kind(1.0_8) ! real*8
#endif
  integer,parameter::MAXOUT=1000
  integer,parameter::MAXLEVEL=100

  ! Define integer types (for particle IDs mostly)
  ! Warning: compiler needs to accept fortran standard 2003.
  ! Specific numbers for fortran kinds are, in principle, implementation
  ! dependent, so "i8b=8" with "integer(i8b)" is implementation-dependent.
  ! See portability note for non-gcc compilers: (boud 2016-11-29) -
  ! https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
  ! The selected_int_kind approach below is more likely to be portable:
  integer,parameter::i4b=selected_int_kind(9) ! since log(2*10^9)/log(2)=30.9
#ifndef LONGINT
  !integer,parameter::i8b=4  ! default long int are 4-byte int
  integer,parameter::i8b=selected_int_kind(9) ! since log(2*10^9)/log(2)=30.9
#else
  integer, parameter :: i8b=selected_int_kind(18) ! since log(2*10^18)/log(2)=60.8
  !integer,parameter::i8b=8  ! long int are 8-byte int
#endif /* LONGINT */

  ! Number of dimensions
#ifndef NDIM
  integer,parameter::ndim=1
#else
  integer,parameter::ndim=NDIM
#endif
  integer,parameter::twotondim=2**ndim
  integer,parameter::threetondim=3**ndim
  integer,parameter::twondim=2*ndim

  ! Vectorization parameter
#ifndef NVECTOR
  integer,parameter::nvector=500  ! Size of vector sweeps
#else
  integer,parameter::nvector=NVECTOR
#endif

  integer, parameter :: nstride = 65536

  ! Run control
  logical::verbose =.false.   ! Write everything
  logical::hydro   =.false.   ! Hydro activated
  logical::pic     =.false.   ! Particle In Cell activated
  logical::poisson =.false.   ! Poisson solver activated
  logical::cosmo   =.false.   ! Cosmology activated
  logical::star    =.false.   ! Star formation activated
  logical::sink    =.false.   ! Sink particles activated
  logical::rt      =.false.   ! Radiative transfer activated
  logical::debug   =.false.   ! Debug mode activated
  logical::static  =.false.   ! Static mode activated
  logical::static_dm=.false.  ! Static mode for dm only activated
  logical::static_gas=.false. ! Static mode for gas only activated
  logical::static_stars=.false.! Static mode for stars only activated
  logical::tracer  =.false.   ! Tracer particles activated
  logical::lightcone=.false.  ! Enable lightcone generation
  logical::clumpfind=.false.  ! Enable clump finder
  logical::aton=.false.       ! Enable ATON coarse grid radiation transfer
  
  ! Mesh parameters
  integer::nx=1,ny=1,nz=1                  ! Number of coarse cells in each dimension
  integer::levelmin=1                      ! Full refinement up to levelmin
  integer::nlevelmax=1                     ! Maximum number of level
  integer::ngridmax=0                      ! Maximum number of grids
  integer,dimension(1:MAXLEVEL)::nexpand=1 ! Number of mesh expansion
  integer::nexpand_bound=1                 ! Number of mesh expansion for virtual boundaries
  real(dp)::boxlen=1.0D0                   ! Box length along x direction
  character(len=128)::ordering='hilbert'
  logical::cost_weighting=.true.           ! Activate load balancing according to cpu time
  ! Recursive bisection tree parameters
  integer::nbilevelmax=1                   ! Max steps of bisection partitioning
  integer::nbinodes=3                      ! Max number of internal nodes
  integer::nbileafnodes=2                  ! Max number of leaf (terminal) nodes
  real(dp)::bisec_tol=0.05d0               ! Tolerance for bisection load balancing

                                 ! Step parameters
  integer::nrestart=0            ! New run or backup file number
  integer::nrestart_quad=0       ! Restart with double precision Hilbert keys
  real(dp)::trestart=0           ! Restart time
  logical::restart_remap=.false. ! Force load balance on restart
  integer::nstepmax=1000000      ! Maximum number of time steps
  integer::ncontrol=1            ! Write control variables
  integer::nremap=0              ! Load balancing frequency (0: never)
  integer,allocatable,dimension(:)::remap_pscalar
  
  ! Output parameters
  integer::iout=1                ! Increment for output times
  integer::ifout=1               ! Increment for output files
  integer::iback=1               ! Increment for backup files
  integer::noutput=1             ! Total number of outputs
  integer::foutput=1000000       ! Frequency of outputs
  logical::gadget_output=.false. ! Output in gadget format
  logical::output_now=.false.    ! write output next step
  real(dp)::walltime_hrs=-1      ! Wallclock time for submitted job
  real(dp)::minutes_dump=1       ! Dump an output minutes before walltime ends

  ! Lightcone parameters
  real(dp)::thetay_cone=12.5d0
  real(dp)::thetaz_cone=12.5d0
  real(dp)::zmax_cone=2d0

  ! Cosmology and physical parameters
  real(dp)::boxlen_ini               ! Box size in h-1 Mpc
  real(dp)::omega_b=0.045d0          ! Omega Baryon
  real(dp)::omega_m=1                ! Omega Matter
  real(dp)::omega_l=0                ! Omega Lambda
  real(dp)::omega_k=0                ! Omega Curvature
  real(dp)::h0     =1                ! Hubble constant in km/s/Mpc
  real(dp)::aexp   =1                ! Current expansion factor
  real(dp)::hexp   =0                ! Current Hubble parameter
  real(dp)::texp   =0                ! Current proper time
  real(dp)::n_sink = -1              ! Sink particle density threshold in H/cc
  real(dp)::rho_sink = -1            ! Sink particle density threshold in g/cc
  real(dp)::d_sink = -1              ! Sink particle density threshold in user units
  real(dp)::m_star =-1               ! Star particle mass in units of mass_sph
  real(dp)::n_star =0.1d0            ! Star formation density threshold in H/cc
  real(dp)::eps_star=0               ! Global star formation efficiency
  real(dp)::T2_star=0                ! Typical ISM polytropic temperature
  real(dp)::g_star =1.6d0            ! Typical ISM polytropic index
  real(dp)::jeans_ncells=-1          ! Jeans polytropic EOS
  real(dp)::del_star=200             ! Minimum overdensity to define ISM
  real(dp)::eta_sn =0                ! Supernova mass fraction
  real(dp)::eta_ssn=0.95d0           ! Single supernova ejected mass fraction (sf_imf=.true. only)
  real(dp)::yield  =0                ! Supernova yield
  real(dp)::f_ek   =1                ! Supernovae kinetic energy fraction (only between 0 and 1)
  real(dp)::rbubble=0                ! Supernovae superbubble radius in pc
  real(dp)::f_w    =0                ! Supernovae mass loading factor
  integer ::ndebris=1                ! Supernovae debris particle number
  real(dp)::mass_gmc=-1              ! Stochastic exploding GMC mass
  real(dp)::z_ave  =0                ! Average metal abundance
  real(dp)::B_ave  =0                ! Average magnetic field
  real(dp)::z_reion=8.5d0            ! Reionization redshift
  real(dp)::T2_start                 ! Starting gas temperature
  real(dp)::T2max=huge(1._dp)        ! Temperature ceiling for cooling_fine
  real(dp)::t_diss =20               ! Dissipation timescale for feedback
  real(dp)::t_sne =10                ! Supernova blast time
  real(dp)::J21    =0                ! UV flux at threshold in 10^21 units
  real(dp)::a_spec =1                ! Slope of the UV spectrum
  real(dp)::beta_fix=0               ! Pressure fix parameter
  real(dp)::kappa_IR=0               ! IR dust opacity
  real(dp)::ind_rsink=4              ! Number of cells defining the radius of the sphere where AGN feedback is active
  real(dp)::ir_eff=0.75d0            ! efficiency of the IR feedback (only when ir_feedback=.true.)
  real(dp)::sf_trelax=0              ! Relaxation time for star formation (cosmo=.false. only)
  real(dp)::sf_tdiss=0               ! Dissipation timescale for subgrid turbulence in units of turbulent crossing time
  integer::sf_model=3                ! Virial star formation model
  integer::nlevel_collapse=3         ! Number of levels to follow initial dark matter collapse (cosmo=.true. only)
  real(dp)::mass_star_max=120        ! Maximum mass of a star in solar mass
  real(dp)::mass_sne_min=10          ! Minimum mass of a single supernova in solar mass
  integer::momentum_feedback=0       ! Use supernovae momentum feedback if cooling radius not resolved

  logical ::self_shielding=.false.
  logical ::pressure_fix=.false.
  logical ::nordlund_fix=.true.
  logical ::cooling=.false.
  logical ::neq_chem=.false.            ! Non-equilbrium chemistry activated
  logical ::isothermal=.false.
  logical ::metal=.false.
  logical ::haardt_madau=.false.
  logical ::delayed_cooling=.false.
  logical ::improved_h2_cooling=.false.  ! Improved collisional dissociation (H2) cooling
  logical ::smbh=.false.
  logical ::agn=.false.
  logical ::use_proper_time=.false.
  logical ::convert_birth_times=.false. ! Convert stellar birthtimes: conformal -> proper
  logical ::ir_feedback=.false.         ! Activate ir feedback from accreting sinks
  logical ::sf_virial=.false.           ! Activate SF Virial criterion
  logical ::sf_log_properties=.false.   ! Log in ascii files birth properties of stars and supernovae
  logical ::sf_imf=.false.              ! Activate IMF sampling for SN feedback when resolution allows it
  logical ::sf_compressive=.false.      ! Advect compressive and solenoidal turbulence terms separately

  ! Cosmic Rays parameters
  logical::conduction  =.false.  ! Conduction module activated
  logical::conduction_ion=.false.! Conduction of ions activated
  logical::isotrope_cond=.false. ! Activate isotropic conduction instead of anisotropic
  logical::semi_implicit=.false. ! Semi-implicit integrator for anisotropic conduction (DOES NOT WORK)
  logical::coupling    =.false.  ! Couple the two temperatures within the conduction module (old feature)
  logical::coupling_out_of_conduction=.true. ! Couple the two temperatures outside the conduction module
  logical::saturation=.true.     ! Saturation of the heat flux for large mean free path
  logical::cr_diffusion=.false.  ! Cosmic ray diffusion module activated
  logical::fix_temp_diff=.true.  ! Cosmic ray diffusion and temperature conduction fix
  logical::twotemp     =.false.  ! Two-temperatures (Te,Ti) model for conduction
  logical::alfven_diff_coeff=.false. ! CR diffusion coeffcient dependant on the Alfvenic Mach number
  real(dp)::Tfloor=1.d0          ! Temperature floor for conduction temperature fix in K
  real(dp)::TCRmax=-1.d0         ! Maximum CR Temperature for CR diffusion step in K
  real(dp)::TCRmin=-1.d0         ! Minimum CR Temperature for CR diffusion step in K
  real(dp)::RelVar=0.1d0         ! Relative temperature variation criterion for 2 Temp. coupling
  logical ::cooling_cr=.false.   ! Radiative losses as in Booth+13
  !!!!!logical ::cooling_scalar=.false.   ! passive scalar to turn off the cooling
  logical::streaming_diffusion=.false. ! CR streaming treated as a diffusion term
  logical::streaming_heating=.false.   ! CR streaming heating term
  logical::do_limiter_aniso=.false.    ! Activate slope limiter for anisotropic diffusion
  real(dp)::DCRmax=1d30          ! Maximum allowed CR streaming diffusion coefficient in cgs
  real(dp)::nlength_str=1d30     ! Maximum number of cell sizes for CR streaming diffusion length scale
  real(dp)::fudge_streamboost=1d0! Allow for the Boost of the streaming velocity
  logical::debug_bicg=.false. ! Debug mode activated for BiCGSTAB
  logical ::slopelim_cond=.false. ! Conduction slope limiter. TODO for the asymetric scheme
  real(dp)::fecr=0d0             ! SN fraction of CR energy
  real(dp)::dtc12,dtc23,dtc34,dtcneig,dtcfath
  real(dp)::alpha_limiter=1.0d0 ! Choose the alpha parameter in the limiter of the normal component (0<alpha<=1) => 1 switches to asymmetric scheme (other values are not good slope limiters)
  real(dp)::epsilon_restartbicg=0.0d0 ! Criterion to restart the value of rzero in BiCGSTAB when rho^2 below some threshold
  integer::slope_limiter_aniso=0 ! Choose your slope_limiter for anisotropic diffusion (0:none,1:MinMod, 2:MonCen)
  logical::explicit_diff=.false.    !Uses an explicit solver, only works for a run without streaming for now
  logical::autoswitch_diff=.false.  !Allows to switch between implicit and explicit solver. Uses explicit if dt_hydro/dt_diff < nmax_explicit, implicit otherwise
  logical::sts_diff=.false.
  integer::nmax_explicit=100        !Limit for hydro to diffusion timestep ratio to define if explicit or implicit solver will be used

  ! Output times
  real(dp),dimension(1:MAXOUT)::aout=1.1d0      ! Output expansion factors
  real(dp),dimension(1:MAXOUT)::tout=HUGE(1.0D0)! Output times

  ! Movie
  integer,parameter::NMOV=5
  integer::imovout=0             ! Increment for output times
  integer::imov=1                ! Initialize
  real(kind=8)::tstartmov=0,astartmov=0
  real(kind=8)::tendmov=0,aendmov=0
  real(kind=8),allocatable,dimension(:)::amovout,tmovout
  logical::movie=.false.
  integer::nw_frame=512 ! prev: nx_frame, width of frame in pixels
  integer::nh_frame=512 ! prev: ny_frame, height of frame in pixels
  integer::levelmax_frame=0
  real(kind=8),dimension(1:4*NMOV)::xcentre_frame=0d0
  real(kind=8),dimension(1:4*NMOV)::ycentre_frame=0d0
  real(kind=8),dimension(1:4*NMOV)::zcentre_frame=0d0
  real(kind=8),dimension(1:2*NMOV)::deltax_frame=0d0
  real(kind=8),dimension(1:2*NMOV)::deltay_frame=0d0
  real(kind=8),dimension(1:2*NMOV)::deltaz_frame=0d0
  real(kind=8),dimension(1:NMOV)::dtheta_camera=0d0
  real(kind=8),dimension(1:NMOV)::dphi_camera=0d0
  real(kind=8),dimension(1:NMOV)::theta_camera=0d0
  real(kind=8),dimension(1:NMOV)::phi_camera=0d0
  real(kind=8),dimension(1:NMOV)::tstart_theta_camera=0d0
  real(kind=8),dimension(1:NMOV)::tstart_phi_camera=0d0
  real(kind=8),dimension(1:NMOV)::tend_theta_camera=0d0
  real(kind=8),dimension(1:NMOV)::tend_phi_camera=0d0
  real(kind=8),dimension(1:NMOV)::focal_camera=0d0
  real(kind=8),dimension(1:NMOV)::dist_camera=0d0
  real(kind=8),dimension(1:NMOV)::ddist_camera=0d0
  real(kind=8),dimension(1:NMOV)::smooth_frame=1d0
  real(kind=8),dimension(1:NMOV)::varmin_frame=-1d60
  real(kind=8),dimension(1:NMOV)::varmax_frame=1d60
  integer,dimension(1:NMOV)::ivar_frame=0
  logical,dimension(1:NMOV)::perspective_camera=.false.
  logical,dimension(1:NMOV)::zoom_only_frame=.false.
  character(LEN=NMOV)::proj_axis='z' ! x->x, y->y, projection along z
  character(LEN=6),dimension(1:NMOV)::shader_frame='square'
  character(LEN=10),dimension(1:NMOV)::method_frame='mean_mass'
  character(len=10),dimension(1:50)::movie_vars_txt=''
  integer::n_movie_vars
  integer::i_mv_temp=-1,  i_mv_dens=-1,       i_mv_p=-1
  integer::i_mv_speed=-1, i_mv_metallicity=-1
  integer::i_mv_vx=-1,    i_mv_vy=-1,         i_mv_vz=-1
  integer::i_mv_dm=-1,    i_mv_stars=-1,      i_mv_lum=-1
  integer::i_mv_var=-1,   i_mv_xh2=-1,        i_mv_xhi=-1
  integer:: i_mv_xhii=-1, i_mv_xheii=-1,      i_mv_xheiii=-1
  integer::i_mv_fp=-1,    i_mv_pmag=-1
  integer,dimension(1:50)::movie_vars=-1
  integer,dimension(1:50)::movie_var_number=1


                                                 ! Refinement parameters for each level
  real(dp),dimension(1:MAXLEVEL)::m_refine =-1   ! Lagrangian threshold
  real(dp),dimension(1:MAXLEVEL)::r_refine =-1   ! Radius of refinement region
  real(dp),dimension(1:MAXLEVEL)::x_refine = 0   ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::y_refine = 0   ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::z_refine = 0   ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::exp_refine = 2 ! Exponent for distance
  real(dp),dimension(1:MAXLEVEL)::a_refine = 1   ! Ellipticity (Y/X)
  real(dp),dimension(1:MAXLEVEL)::b_refine = 1   ! Ellipticity (Z/X)
  real(dp)::var_cut_refine=-1                    ! Threshold for variable-based refinement
  real(dp)::mass_cut_refine=-1                   ! Mass threshold for particle-based refinement
  integer::ivar_refine=-1                        ! Variable index for refinement
  logical::sink_refine=.false.                   ! Fully refine on sink particles

  ! Initial condition files for each level
  logical::multiple=.false.
  character(LEN=80),dimension(1:MAXLEVEL)::initfile=' '
  character(LEN=20)::filetype='ascii'

  ! Initial condition regions parameters
  integer,parameter::MAXREGION=100
  integer                           ::nregion=0
  character(LEN=10),dimension(1:MAXREGION)::region_type='square'
  real(dp),dimension(1:MAXREGION)   ::x_center=0
  real(dp),dimension(1:MAXREGION)   ::y_center=0
  real(dp),dimension(1:MAXREGION)   ::z_center=0
  real(dp),dimension(1:MAXREGION)   ::length_x=1d10
  real(dp),dimension(1:MAXREGION)   ::length_y=1d10
  real(dp),dimension(1:MAXREGION)   ::length_z=1d10
  real(dp),dimension(1:MAXREGION)   ::exp_region=2

  ! Boundary conditions parameters
  integer,parameter::MAXBOUND=100
  logical                           ::simple_boundary=.false.
  integer                           ::nboundary=0
  integer                           ::icoarse_min=0
  integer                           ::icoarse_max=0
  integer                           ::jcoarse_min=0
  integer                           ::jcoarse_max=0
  integer                           ::kcoarse_min=0
  integer                           ::kcoarse_max=0
  integer ,dimension(1:MAXBOUND)    ::boundary_type=0
  integer ,dimension(1:MAXBOUND)    ::ibound_min=0
  integer ,dimension(1:MAXBOUND)    ::ibound_max=0
  integer ,dimension(1:MAXBOUND)    ::jbound_min=0
  integer ,dimension(1:MAXBOUND)    ::jbound_max=0
  integer ,dimension(1:MAXBOUND)    ::kbound_min=0
  integer ,dimension(1:MAXBOUND)    ::kbound_max=0
  logical                           ::no_inflow=.false.

  ! Number of processes sharing one token
  ! Only one process can write at a time in an I/O group
  integer::IOGROUPSIZE=0           ! Main snapshot
  integer::IOGROUPSIZECONE=0       ! Lightcone
  integer::IOGROUPSIZEREP=0        ! Subfolder size
  logical::withoutmkdir=.false.    ! If true mkdir should be done before the run
  logical::print_when_io=.false.   ! If true print when IO
  logical::synchro_when_io=.false. ! If true synchronize when IO

  ! Allow increasing refinement for constant physical resolution
  logical::constant_physical_resolution=.true.
  
  character(len=10)::feedback_model='thermal' ! 'mechanical', 'agertz', 'mom2', 'kinetic'

  ! Taysun Kimm's Mechanical feedback stuff:
  ! Efficiency of stellar feedback energy used to heat up/blow out the gas:
  real(dp)::E_SNII=1d51        ! different from ESN used in feedback.f90 
  ! Realistic time delay for individual star particle. t_delay is oldest age to consider if set:
  logical ::sn2_real_delay=.false.
  ! Activate mechanical feedback
  logical ::mechanical_geen=.false.  ! Geen boost
  logical ::log_mfb=.false.
  logical ::log_mfb_mega=.false.
  real(dp)::A_SN=2.5d5
  real(dp)::A_SN_Geen=5d5
  real(dp)::expN_SN=-2d0/17d0
  logical ::mechanical_bpass=.false.  ! use the SN rates based on bpass_v2_300 model

  ! Star formation stuff:
  character(len=10)::star_maker='density' ! density, federrath, FK2
  real(dp)::T2thres_SF=2d4    ! Temperature threshold
  real(dp)::fstar_min=1d0     ! Mstar,min = nH*dx_min^3*fstar_min
  real(dp)::M_SNII=10d0       ! Mean progenitor mass of the TypeII SNe
  real(dp)::sf_lam=-1.        ! Jeans length criterion for Thermo-turbulent SF
  logical::sf_lmax=.false.    ! Form stars only on maximum allowed level
  character(len=8 )::star_imf=''          ! salpeter,kroupa,(chabrier)
  real(dp)::n_gmc=1D-1        ! Density above which Hopkins and Cen SF routines are evaluated
  integer ::nsn2mass=-1       ! Star particle mass in units of the number of SN
  real(dp)::eps_star_loc=0.5d0! Local star formation efficiency, unresolved (for us) proto-stellar feedback
  ! Runaway stars:
  real(dp)::sf_kick_kms=-1    ! Upper value of the speed distribution. No runaway stars if negative.
  real(dp)::sf_kick_pow=1d0   ! Index of the power law (used only if negative).
  ! The dense gas core-to-star efficiency would be 1.0 without feedback (Federrath & Klessen 2012)
  integer ::nlevelmax_current=0           ! nlevelmax that is actually used
  ! POPIII stars ?  
  logical  :: pop3=.false.
  real(dp) :: pop3_mass = -1.
  real(dp) :: Zcrit_pop3 = 1e-4


#ifdef SOLVERmhd
  integer,parameter::nvarMHD=NVAR+3
#else
  integer,parameter::nvarMHD=NVAR
#endif
  
  ! ----------- Oscar additions --------------
  real(dp)::fh2_min=0.0       ! Minimum fraction of H_2 for SF
  real(dp)::fh2_rho=1.0d5          !Threshold where we assume we have reached 100% molecular gas (regardless of metals)
  real(dp)::vmaxFB = 1.0d3  !maximum feedback velocity
  real(dp)::eta_rap =1.0D0     ! Radiation fudge factor for Oscar's model
  real(dp)::tau_IR=-1          ! infrared optical depth
  logical ::metalscaling=.true. !scale Prad and winds with metallicity
  real(dp)::Nrcool=3.0   !Number of resolved cooling radii for thermal SNe
  real(dp)::mstar_dt=1.0     ! sampling mass
  real(dp)::maxadvfb=1.0d10  ! maximum velocity allowed in km/s
  logical ::KMT09=.false.  !KMT09 f_h2 star formation
  logical ::momST=.false.    !S-T momentum (Blondin et al. 1998)
  logical ::ysc_stats=.false.
  logical ::supernovae=.true.
  logical ::winds=.true.
  logical ::energy=.true.
  logical ::momentum=.true.
  logical ::radpressure=.false.  !Radiation pressure on dust from young stars
  logical ::efb_advection=.true.
  logical ::pnonthermal=.false.
  logical ::fbsafety=.false.
  logical ::is_oxygen=.false.
   

end module amr_parameters
