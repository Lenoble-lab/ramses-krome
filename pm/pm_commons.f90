module pm_commons

  use amr_parameters
  use pm_parameters
  use random

  implicit none

  ! Sink particle related arrays
  real(dp),allocatable,dimension(:)    ::msink,xmsink
  real(dp),allocatable,dimension(:)    ::msink_new,msink_all
  real(dp),allocatable,dimension(:)    ::msmbh,msmbh_new,msmbh_all
  real(dp),allocatable,dimension(:)    ::oksink_new,oksink_all
  real(dp),allocatable,dimension(:)    ::tsink,tsink_new,tsink_all
  real(dp),allocatable,dimension(:)    ::dMsink_overdt,dMBHoverdt
  real(dp),allocatable,dimension(:)    ::dMsmbh_overdt,dMBHoverdt_smbh
  real(dp),allocatable,dimension(:)    ::rho_gas,volume_gas,eps_sink,c2sink
  real(dp),allocatable,dimension(:,:)  ::vel_gas
  real(dp),allocatable,dimension(:)    ::delta_mass,delta_mass_new,delta_mass_all
  real(dp),allocatable,dimension(:)    ::wden,weth,wvol,wdiv,wden_new,weth_new,wvol_new,wdiv_new
  real(dp),allocatable,dimension(:,:)  ::wmom,wmom_new
  real(dp),allocatable,dimension(:,:)  ::vsink,vsink_new,vsink_all
  real(dp),allocatable,dimension(:,:)  ::fsink,fsink_new,fsink_all
  real(dp),allocatable,dimension(:,:,:)::vsnew,vsold
  real(dp),allocatable,dimension(:,:,:)::fsink_partial,sink_jump
  real(dp),allocatable,dimension(:,:)  ::lsink,lsink_new,lsink_all
  real(dp),allocatable,dimension(:,:)  ::xsink,xsink_new,xsink_all
  real(dp),allocatable,dimension(:,:)  ::weighted_density,weighted_volume,weighted_ethermal,weighted_divergence
  real(dp),allocatable,dimension(:,:,:)::weighted_momentum
  real(dp),allocatable,dimension(:)    ::rho_sink_tff
  real(dp),allocatable,dimension(:)    ::msum_overlap
  integer,allocatable,dimension(:)     ::idsink,idsink_new,idsink_old,idsink_all
  logical,allocatable,dimension(:)     ::ok_blast_agn,ok_blast_agn_all
  logical,allocatable,dimension(:)     ::direct_force_sink
  logical,allocatable,dimension(:)     ::new_born,new_born_all,new_born_new
  integer,allocatable,dimension(:)     ::idsink_sort
  integer::ncloud_sink,ncloud_sink_massive
  integer::nindsink=0
  integer::sinkint_level=0         ! maximum level currently active is where the global sink variables are updated
  real(dp)::ssoft                  ! sink softening lenght in code units

  ! Particles related arrays
  real(dp),allocatable,dimension(:,:)  ::xp       ! Positions
  real(dp),allocatable,dimension(:,:)  ::vp       ! Velocities
  real(dp),allocatable,dimension(:)    ::mp       ! Masses
  real(dp),allocatable,dimension(:)    ::mp0      ! Initial masses
#ifdef OUTPUT_PARTICLE_POTENTIAL
  real(dp),allocatable,dimension(:)    ::ptcl_phi ! Potential of particle added by AP for output purposes
#endif
#ifdef NTRACEGROUPS
  integer ,allocatable,dimension(:)    ::ptracegroup
#endif
  real(dp),allocatable,dimension(:)    ::tp       ! Birth epoch
  real(dp),allocatable,dimension(:)    ::zp       ! Birth metallicity
  real(dp),allocatable,dimension(:)    ::zp_ox    ! Birth oxygen metallicity (O Agertz)
  integer ,allocatable,dimension(:)    ::nextp    ! Next particle in list
  integer ,allocatable,dimension(:)    ::prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)    ::levelp   ! Current level of particle
  integer(i8b),allocatable,dimension(:)::idp    ! Identity of particle
  ! Tree related arrays
  integer ,allocatable,dimension(:)    ::headp    ! Head particle in grid
  integer ,allocatable,dimension(:)    ::tailp    ! Tail particle in grid
  integer ,allocatable,dimension(:)    ::numbp    ! Number of particles in grid
  ! Global particle linked lists
  integer::headp_free,tailp_free,numbp_free=0,numbp_free_tot=0
  ! Local and current seed for random number generator
  integer,dimension(IRandNumSize) :: localseed=-1

  ! Particle types
  integer, parameter   :: NFAMILIES=5
  integer(1),parameter :: FAM_DM=1, FAM_STAR=2, FAM_CLOUD=3, FAM_DEBRIS=4, FAM_OTHER=5, FAM_UNDEF=127
  integer(1),parameter :: FAM_TRACER_GAS=0
  integer(1),parameter :: FAM_TRACER_DM=-1, FAM_TRACER_STAR=-2, FAM_TRACER_CLOUD=-3, FAM_TRACER_DEBRIS=-4, FAM_TRACER_OTHER=-5

  ! Customize here for particle tags within particle types (e.g. different kind of stars).
  ! Note that the type should be integer(1) (1 byte integers) for memory concerns.
  ! Also don't forget to create a function is_<type>_<tag>. See the wiki for a more complete example.
  ! By default, the tag is always 0.

  ! Particle keys for outputing. They should match the above particle
  ! types, except for 'under' family
  character(len=13), dimension(-NFAMILIES:NFAMILIES), parameter :: particle_family_keys = (/ &
       ' other_tracer', 'debris_tracer', ' cloud_tracer', '  star_tracer', ' other_tracer', &
       '   gas_tracer', &
       '           DM', '         star', '        cloud', '       debris', '        other'/)

  type(part_t), allocatable, dimension(:) :: typep  ! Particle type array

contains
  function cross(a,b)
    use amr_parameters, only:dp
    real(dp),dimension(1:3)::a,b
    real(dp),dimension(1:3)::cross
    !computes the cross product c= a x b
    cross(1)=a(2)*b(3)-a(3)*b(2)
    cross(2)=a(3)*b(1)-a(1)*b(3)
    cross(3)=a(1)*b(2)-a(2)*b(1)
  end function cross

  elemental logical pure function is_DM(typep)
    type(part_t), intent(in) :: typep
    is_DM = typep%family == FAM_DM
  end function is_DM

  elemental logical pure function is_star(typep)
    type(part_t), intent(in) :: typep
    is_star = typep%family == FAM_STAR
  end function is_star

  elemental logical pure function is_cloud(typep)
    type(part_t), intent(in) :: typep
    is_cloud = typep%family == FAM_CLOUD
  end function is_cloud

  elemental logical pure function is_debris(typep)
    type(part_t), intent(in) :: typep
    is_debris = typep%family == FAM_DEBRIS
  end function is_debris

  elemental logical pure function is_tracer(typep)
    type(part_t), intent(in) :: typep
    is_tracer = typep%family <= 0
  end function is_tracer

  elemental logical pure function is_not_tracer(typep)
    type(part_t), intent(in) :: typep
    is_not_tracer = typep%family > 0
  end function is_not_tracer

  elemental logical pure function is_not_DM(typep)
    type(part_t), intent(in) :: typep
    is_not_DM = typep%family /= FAM_DM
  end function is_not_DM
  

  elemental function part2int (part)
    ! Convert a particle into an integer
    ! This saves some space e.g. when communicating
    integer :: part2int
    type(part_t), intent(in) :: part

    ! This is the largest value for integer(1)
    integer, parameter :: a = 128, b = 2*a

    part2int = (int(part%family) + a) * b + (int(part%tag) + a)
  end function part2int

  elemental function int2part(index)
    ! Convert from an index to particle type
    type(part_t) :: int2part
    integer, intent(in) :: index

    ! This is the largest value for integer(1)
    integer, parameter :: a = 128, b = 2*a

    int2part%family = int(index / b - a, 1)
    int2part%tag = int(mod(index, b) - a, 1)
  end function int2part

  function props2type(idpii, tpii, mpii)
    use amr_commons
    use pm_parameters, only : part_t

    ! Converts from "old" ramses to "new" ramses
    !
    ! Here's the match, add yours here for backward compatibility purposes
    ! DM     tpii == 0
    ! stars  tpii != 0 and idpii > 0
    ! sinks  tpii != 0 and idpii < 0
    !
    ! This is mostly for support of GRAFFIC I/O.
    ! The reason we use idpii instead of idp is to prevent name clashes
    real(dp), intent(in) :: tpii, mpii
    integer, intent(in)  :: idpii

    type(part_t) :: props2type

    if (tpii == 0) then
       props2type%family = FAM_DM
    else if (idpii > 0) then
       props2type%family = FAM_STAR
    else if (idpii < 0) then
       props2type%family = FAM_CLOUD
    else if (mpii == 0) then
       props2type%family = FAM_TRACER_GAS
    end if
    props2type%tag = 0
  end function props2type
end module pm_commons

!####################################################################
!####################################################################
!####################################################################
module mechanical_commons
   use amr_commons, ONLY:dp,ndim,ncpu
#ifndef NCHEM
   integer,parameter::nchem=0 ! number of chemical elements (max: 8)
#else
   integer,parameter::nchem=NCHEM
#endif
   integer(kind=4),parameter::nvarSN=13+nchem

   ! Important Note: SN stands for SN cell, not SN particle (SNp) 

   ! Array to define neighbors for SN
   integer, parameter::nSNnei=48  ! number of neighboring cells to deposit mass/momentum/energy
   real(dp),parameter::nSNcen=4   ! number of cells corresponding to the central cell to deposit mass
   real(dp),dimension(1:3,1:nSNnei)::xSNnei
   real(dp),dimension(1:3,1:nSNnei)::vSNnei
   real(dp)::f_LOAD,f_LEFT,f_ESN,f_PCAN
   ! SN cells that are needed to communicate across different CPUs
   ! note that SNe are processed on a cell-by-cell basis
   ! hence the central position is taken as the central leaf cell
   integer ::ncomm_SN   ! the number of cells to be communicated (specific to myid)

   ! momentum input
   ! p_sn = A_SN*nH**(alpha)*ESN**(beta)*ZpSN**(gamma)
   ! ex) Thornton et al.
   !     A_SN = 3e5, alphaN = -2/17, beta = 16/17, gamma = -0.14
   ! ex) Kim & Ostriker (2015) uniform case
   !     A_SN = 2.17e5, alpha = -0.13, beta = 0.93
   real(dp),parameter:: expE_SN=+16d0/17d0
   real(dp),parameter:: expZ_SN=-0.14
   real(dp),parameter:: expN_SN_boost=-0.15
   real(dp),parameter:: expE_SN_boost=0.90

   ! For PopIII
   real(dp),parameter:: A_pop3=2.5d5
   real(dp) :: expN_pop3 =  -2d0/17d0
   real(dp) :: expE_pop3 = +16d0/17d0
   real(dp) :: expZ_pop3 = -0.14

#ifndef WITHOUTMPI
  ! common lists for SNe across different cpus
   integer, parameter ::ncomm_max = 50000
   integer ,dimension(1:ncomm_max)::iSN_comm  ! cpu list
   real(dp),dimension(1:ncomm_max)::nSN_comm          ! number of SNe
   real(dp),dimension(1:ncomm_max)::mSN_comm          ! gas mass of SNe 
   real(dp),dimension(1:ncomm_max)::mloadSN_comm      ! ejecta mass + gas entrained
   real(dp),dimension(1:ncomm_max)::eloadSN_comm      ! kinetic energy of ejecta + gas entrained
   real(dp),dimension(1:ncomm_max)::mZloadSN_comm     ! metals ejected + entrained
   real(dp),dimension(1:3,1:ncomm_max)::xSN_comm      ! pos of SNe host cell (leaf)
   real(dp),dimension(1:3,1:ncomm_max)::ploadSN_comm  ! momentum from original star + gas entrained
   real(dp),dimension(1:ncomm_max)::floadSN_comm      ! fraction of gas to be loaded from the central cell
   integer,dimension(:,:),allocatable::icpuSN_comm,icpuSN_comm_mpi
   integer,dimension(:)  ,allocatable::ncomm_SN_cpu,ncomm_SN_mpi
   real(dp),dimension(1:ncomm_max)::rSt_comm          ! Stromgren radius in pc
   real(dp),dimension(1:ncomm_max,1:nchem)::mchloadSN_comm     ! chemical species ejected + entrained
#endif

   ! refinement for resolved feedback (used 'ring' instead of 'shell' to be more catchy)
   integer::ncshell3                            ! the maximum number of cells within a (1+2*nshell_re
   integer,dimension(:,:),allocatable::xrefnei  ! relative position of the neighboring cells
   integer,dimension(:),allocatable::irefnei    ! index of the nei cells
   integer,dimension(:),allocatable::lrefnei    ! level of the neighboring cells
   integer,dimension(:),allocatable::icellnei   ! cell index of the neighbors
   real(dp),dimension(:),allocatable::mrefnei   ! gas mass within each cell in Msun
   real(dp),dimension(:),allocatable::mzrefnei  ! metal mass within each cell in Msun
   integer,dimension(:),allocatable::nrefnei_ring  ! cumulative number of nei cells per ring - useful
   real(dp),dimension(:),allocatable::mrefnei_ring ! gas mass within each shell in Msun
   real(dp),dimension(:),allocatable::mzrefnei_ring! metal mass within each shell in Msun

   integer,dimension(:),allocatable::icommr     ! for communication

   ! parameters for refinement based on feedback
   integer ::nshell_resolve=3   ! r_shell will be resolved on 3 cells in one direction
   real(dp)::nsn_resolve=1d0    ! the number of SN that we want to resolve
   real(dp)::reduce_mass_tr=27  ! 27 means mass_tr/27 should be resolved

   ! For Binary stellar evolution; BPASS_v2
   integer,parameter::nZ_bpass2=11
   integer,parameter::nt_bpass2=41
   real(dp),dimension(1:nZ_bpass2)::Zgrid_bpass2   ! Logarithmic metallicity grid
   real(dp),dimension(1:nt_bpass2,1:nZ_bpass2)::snr_bpass2   ! Supernova Type II rates
   real(dp),dimension(1:nZ_bpass2)::snr_bpass2_sum ! Total Supernova Type II rates
   real(dp),dimension(1:nZ_bpass2)::snr_bpass2_max ! Maximum Supernova Type II rates
   real(dp)::binsize_yr_bpass2,logmin_yr_bpass2,logmax_yr_bpass2

   ! Stellar winds stuff: not used for EoR
   logical::stellar_winds=.false.
   character(LEN=256)::stellar_winds_file='/home/kimm/soft/lib/swind_krp_pagb.dat'
   character(LEN=2),dimension(1:8):: chem_list=(/'H ','O ','Fe','Mg','C ','N ','Si','S '/)

   ! SN Type Ia  stuff: not used for EoR
   logical ::snIa=.false.
   real(dp)::A_snIa=0.0013
   real(dp)::E_SNIa=1d51
   logical ::variable_yield_SNII=.false.  ! TypeII yields are computed according to metallicities base
   real(dp)::mass_loss_boost=1d0 ! 1.5 for Chabrier

   ! For chemical abundance due to SN II
   real(dp)::Zejecta_chem_II(1:nchem)

   ! For chemical abundance due to SN Ia
   real(dp)::mejecta_Ia
   real(dp)::Zejecta_chem_Ia(1:nchem)

end module
!####################################################################
!####################################################################
!####################################################################
subroutine init_mechanical
   use amr_commons
   use mechanical_commons
#if NCHEM>0
   use hydro_parameters, ONLY:ichem
#endif
   implicit none
   integer::i,j,k,ind,indall
   real(kind=dp)::x,y,z,r
   logical::ok
   integer,allocatable::ltmp(:)
   integer::ncshell,iring,i2,j2,k2
#if NCHEM>0
   character(len=2)::element_name
   integer::ich
#endif

   !------------------------------------------------
   ! Warning messages
   !------------------------------------------------
   ok=.false.
   if(.not.metal .and. trim(feedback_model).eq.'mechanical')then
      print *, '>>> mechanical Err: Please turn on metal'
      ok=.true.
   endif
   if(ok) call clean_stop

#ifndef WITHOUTMPI
   allocate(ncomm_SN_cpu(1:ncpu))
   allocate(ncomm_SN_mpi(1:ncpu))
#endif

   ! some parameters
   f_LOAD = nSNnei / dble(nSNcen + nSNnei)
   f_LEFT = nSNcen / dble(nSNcen + nSNnei)
   f_ESN  = 0.676   ! Blondin+(98) at t=trad
   f_PCAN = 0.9387  ! correction due to direct momentum cancellation 
                    ! due to the assumption of 48 neighboring cells 
                    ! even in the uniform case where there are 18 immediate neighbors

   ! Arrays to define neighbors (center=[0,0,0])
   ! normalized to dx = 1 = size of the central leaf cell in which a SN particle sits
   ! from -0.75 to 0.75 
   ind=0
   do k=1,4
   do j=1,4
   do i=1,4
      ok=.true.
      if((i==1.or.i==4).and.&
         (j==1.or.j==4).and.&
         (k==1.or.k==4)) ok=.false. ! edge
      if((i==2.or.i==3).and.&
         (j==2.or.j==3).and.&
         (k==2.or.k==3)) ok=.false. ! centre
      if(ok)then
         ind=ind+1
         x = (i-1)+0.5d0 - 2  
         y = (j-1)+0.5d0 - 2  
         z = (k-1)+0.5d0 - 2  
         r = dsqrt(dble(x*x+y*y+z*z))
         xSNnei(1,ind) = x/2d0
         xSNnei(2,ind) = y/2d0
         xSNnei(3,ind) = z/2d0
         vSNnei(1,ind) = x/r  
         vSNnei(2,ind) = y/r  
         vSNnei(3,ind) = z/r  
         !indall(i+(j-1)*4+(k-1)*4*4) = ind      
      endif
   enddo
   enddo
   enddo

   !=======================================================
   ! For careful refinement (for resolved feedback)
   !=======================================================
   ncshell  = (1+2*nshell_resolve)
   ncshell3 = ncshell**3

   allocate(xrefnei(1:3,1:ncshell3))
   allocate(irefnei(1:ncshell3))
   ! notice that xrefnei,irefnei have a different ordering than the rest of variables
   ! xrefnei <- irefnei <- mrefnei 
   allocate(mrefnei(1:ncshell3))
   allocate(mzrefnei(1:ncshell3))
   allocate(icellnei(1:ncshell3))
   allocate(lrefnei(1:ncshell3))
   allocate(ltmp   (1:ncshell3))

   allocate(nrefnei_ring (0:nshell_resolve))
   allocate(mrefnei_ring (0:nshell_resolve))
   allocate(mzrefnei_ring(0:nshell_resolve))

   allocate(icommr (1:ncshell3))

   xrefnei=0;irefnei=0;nrefnei_ring=1;ltmp=0
   mrefnei=0d0;mrefnei_ring=0d0;icellnei=0;lrefnei=0

   ind=1
   do iring=0,nshell_resolve
      do k=-iring,iring
      do j=-iring,iring
      do i=-iring,iring

         i2=i+nshell_resolve  ! [0,nshell_resolve*2]
         j2=j+nshell_resolve
         k2=k+nshell_resolve

         indall=i2+(j2+k2*ncshell)*ncshell+1

         if(ltmp(indall)==0)then
           ltmp(indall)=1
           irefnei(ind)=indall
           nrefnei_ring(iring)=ind   ! just keep track of the last index
           ind=ind+1
         endif

         xrefnei(1,indall)=i ! [-3,3]
         xrefnei(2,indall)=j ! [-3,3]
         xrefnei(3,indall)=k ! [-3,3]

      end do
      end do
      end do

   end do

   deallocate(ltmp)

#if NCHEM>0
   if(.not.variable_yield_SNII)then ! assume solar case
      do ich=1,nchem
         element_name=chem_list(ich)
         select case (element_name)
            case ('H ')
               Zejecta_chem_II(ich) = 10.**(-0.30967822)
            case ('He')
               Zejecta_chem_II(ich) = 10.**(-0.40330181)
            case ('C ')
               Zejecta_chem_II(ich) = 10.**(-1.9626259)
            case ('N ')
               Zejecta_chem_II(ich) = 10.**(-2.4260355)
            case ('O ')
               Zejecta_chem_II(ich) = 10.**(-1.1213435)
            case ('Mg')
               Zejecta_chem_II(ich) = 10.**(-2.3706062)
            case ('Si')
               Zejecta_chem_II(ich) = 10.**(-2.0431845)
            case ('S ')
               Zejecta_chem_II(ich) = 10.**(-2.2964020)
            case ('Fe')
               Zejecta_chem_II(ich) = 10.**(-2.2126987)
            case default
               Zejecta_chem_II(ich)=0
         end select

      end do
   endif
#endif
  if (snIa) call init_snIa_yield

  ! Binary stellar evolution
  if (mechanical_bpass) call init_SNII_bpass_v2_300

end subroutine init_mechanical
!####################################################################
!####################################################################
!####################################################################
subroutine init_snIa_yield
   use amr_commons
   use mechanical_commons
   implicit none
   real(kind=8)::yield_snIa(1:66)
#if NCHEM>0
   character(len=2)::element_name
   integer::ich
#endif
!----------------------------------------------------------------------------
!  Nomoto et al. (1997,1984) W7 (carbon-deflagration model)
!            (better fit with Tycho observation)
!----------------------------------------------------------------------------
!       12C ,     13C ,   14N ,    15N ,    16O ,    17O ,    18O ,    19F ,
!       20Ne,     21Ne,   22Ne,    23Na,    24Mg,    25Mg,    26Mg,    27Al,
!       28Si,     29Si,   30Si,    31P ,    32S ,    33S ,    34S ,    36S ,
!       35Cl,     37Cl,   36Ar,    38Ar,    40Ar,    39K ,    41K ,    40Ca,
!       42Ca,     43Ca,   44Ca,    46Ca,    48Ca,    45Sc,    46Ti,    47Ti,
!       48Ti,     49Ti,   50Ti,    50V ,    51V ,    50Cr,    52Cr,    53Cr,
!       54Cr,     55Mn,   54Fe,    56Fe,    57Fe,    58Fe,    59Co,    58Ni,
!       60Ni,     61Ni,   62Ni,    64Ni,    63Cu,    65Cu,    64Zn,    66Zn,
!       67Zn,     68Zn                                                   
!----------------------------------------------------------------------------

   yield_snIa = (/ &  ! Msun per SN
    &4.83E-02,1.40E-06,1.16E-06,1.32E-09,1.43E-01,3.54E-08,8.25E-10,5.67E-10,&
    &2.02E-03,8.46E-06,2.49E-03,6.32E-05,8.50E-03,4.05E-05,3.18E-05,9.86E-04,&
    &1.50E-01,8.61E-04,1.74E-03,4.18E-04,8.41E-02,4.50E-04,1.90E-03,3.15E-07,&
    &1.34E-04,3.98E-05,1.49E-02,1.06E-03,1.26E-08,8.52E-05,7.44E-06,1.23E-02,&
    &3.52E-05,1.03E-07,8.86E-06,1.99E-09,7.10E-12,2.47E-07,1.71E-05,6.04E-07,&
    &2.03E-04,1.69E-05,1.26E-05,8.28E-09,5.15E-05,2.71E-04,5.15E-03,7.85E-04,&
    &1.90E-04,8.23E-03,1.04E-01,6.13E-01,2.55E-02,9.63E-04,1.02E-03,1.28E-01,&
    &1.05E-02,2.51E-04,2.66E-03,1.31E-06,1.79E-06,6.83E-07,1.22E-05,2.12E-05,&
    &1.34E-08,1.02E-08/)

   mejecta_Ia=sum(yield_snIa)
#if NCHEM>0
   do ich=1,nchem
      element_name=chem_list(ich)
      select case (element_name)
         case ('C ')
            Zejecta_chem_Ia(ich)=sum(yield_snIa(1:2))/mejecta_Ia
         case ('N ')
            Zejecta_chem_Ia(ich)=sum(yield_snIa(3:4))/mejecta_Ia
         case ('O ')
            Zejecta_chem_Ia(ich)=sum(yield_snIa(5:7))/mejecta_Ia
         case ('Mg')
            Zejecta_chem_Ia(ich)=sum(yield_snIa(13:15))/mejecta_Ia
         case ('Si')
            Zejecta_chem_Ia(ich)=sum(yield_snIa(17:19))/mejecta_Ia
         case ('S ')
            Zejecta_chem_Ia(ich)=sum(yield_snIa(21:24))/mejecta_Ia
         case ('Fe')
            Zejecta_chem_Ia(ich)=sum(yield_snIa(51:54))/mejecta_Ia
         case default
            Zejecta_chem_Ia(ich)=0d0
      end select     
   enddo
#endif
end subroutine init_snIa_yield
!####################################################################
!####################################################################
!####################################################################
subroutine init_SNII_bpass_v2_300
   ! data for BPASSv2
   ! Data from BPASSv2_imf135_300/OUTPUT_POP/
   ! Notice that log_yr should be regularly spaced
   ! The number of SNe is normalised to 10^6 Msun (instantaenous burst)
   ! ~/soft/bpass2/sn/convert_format.pro
   use amr_commons, ONLY:dp
   use mechanical_commons
   implicit none
   integer::iz

   binsize_yr_bpass2=0.100
   logmin_yr_bpass2= 6.000
   logmax_yr_bpass2=10.000
   zgrid_bpass2=(/0.001,0.002,0.003,0.004,0.006,0.008,0.010,0.014,0.020,0.030,0.040/)
   zgrid_bpass2=log10(zgrid_bpass2)

   snr_bpass2=reshape( (/ &
             0.00,    0.00,    0.00,    0.46,   77.96,  224.41,  348.77, &
           502.84,  652.74,  539.69,  908.56,  772.26, 1265.41,  856.63, &
          1480.79, 1148.45,  998.62,  720.07, 1068.05,  281.12,   25.05, &
            18.67,   17.84,    6.13,  118.99,    1.16,  160.16,  195.13, &
             3.12,  222.66,    7.08,  780.25,    0.34,    1.17,  675.08, &
             3.45,    0.00,    0.00,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.53,   72.71,  230.42,  332.05, &
           476.21,  676.36,  534.01,  897.12,  757.50, 1273.40,  909.43, &
          1452.34, 1089.24, 1064.76, 1048.65,  838.29,  160.65,    9.37, &
            11.88,   15.68,    0.00,    0.00,    0.00,    2.07,    0.00, &
             0.62,    2.20,  177.84,  291.99,  660.05,    1.78,  653.03, &
             3.82,    1.24,    3.09,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.46,   73.13,  232.04,  338.47, &
           485.73,  679.15,  528.15,  880.84,  758.31, 1105.45,  969.84, &
          1548.76, 1103.21,  925.31, 1078.13,  419.69,  475.99,    0.03, &
            11.88,   24.01,    0.00,    2.14,    4.15,    2.49,    2.57, &
             2.08,  426.39,  272.02,    0.88,  817.98,    2.01,  469.56, &
             0.60,    0.89,    1.99,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.38,   72.79,  230.91,  343.56, &
           490.47,  679.35,  580.65,  802.86,  770.11, 1109.26,  919.54, &
          1467.59, 1007.99,  970.72, 1090.46,  361.96,  675.09,   30.52, &
             1.21,   17.58,   23.58,    1.10,    0.00,   91.81,    0.80, &
             0.00,  135.07,    3.45,    0.85,    1.20, 1692.05,    0.00, &
             0.18,    0.39,    1.47,    2.24,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,   72.79,  201.07,  290.78, &
           458.55,  609.32,  560.30,  718.60,  828.29, 1165.96,  764.98, &
          1374.80,  965.64,  850.00, 1035.70,   84.66,  566.90,   54.53, &
            39.59,    1.38,    0.25,    1.07,    1.09,    1.36,    2.13, &
             9.61,    0.00,    0.72,  511.86,    0.46, 1315.77,  162.42, &
             0.38,    0.70,    1.24,    3.13,    0.74,    0.24, &
             0.00,    0.00,    0.00,    0.00,   64.16,  218.46,  295.45, &
           447.63,  603.20,  561.27,  710.78,  830.33, 1222.37,  737.60, &
          1345.90,  924.51,  804.49,  627.80,   76.23,  481.26,   53.06, &
            63.19,   31.36,   21.50,    0.79,    0.00,    1.42,    2.50, &
             2.47,    0.00,  130.88,    1.39, 1052.48,  279.24,    1.24, &
           174.10,    2.48,    0.49,    0.31,    0.00,    0.44, &
             0.00,    0.00,    0.00,    0.00,   63.55,  217.12,  301.70, &
           435.42,  609.95,  561.43,  714.21,  827.49, 1263.88,  690.07, &
          1340.79,  734.93,  941.25,  410.05,  321.92,  155.41,  157.58, &
           110.70,   83.85,   36.90,    0.02,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,    1.78,  273.24,  328.59, &
           427.49,  601.57,  578.42,  692.81,  822.85, 1141.12,  797.29, &
          1235.43,  764.37,  914.09,  314.58,  520.01,  234.85,   39.67, &
           341.08,  544.32,   26.80,   13.47,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,    0.00,  282.05,  338.93, &
           434.61,  600.35,  681.23,  810.65,  926.88,  932.90,  696.16, &
          1174.05,  701.36,  889.97,  379.65,   65.21,   85.41,   93.25, &
            70.14,  104.58,   30.24,   27.56,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,    0.00,  287.40,  550.00, &
           340.22,  383.23,  799.91,  816.31,  782.90,  934.82, 1257.31, &
           572.55,  772.16,  349.50,  287.36,   80.06,   74.23,   62.94, &
            70.96,    7.40,    0.01,    2.66,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,    0.00,  285.90,  526.52, &
           404.28,  584.32,  761.42,  698.97, 1087.02,  680.53, 1117.37, &
           695.28,  746.51,  290.25,  144.03,   57.40,   60.99,   71.28, &
            46.82,    0.00,    6.03,    0.00,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,    0.00,    0.00,    0.00, &
             0.00,    0.00,    0.00,    0.00,    0.00,    0.00 /), (/41,11/) )

   ! normalise to 1Msun
   snr_bpass2(:,:) = snr_bpass2(:,:)/1e6

   do iz=1,11
      snr_bpass2_sum(iz) = sum(snr_bpass2(:,iz))
      snr_bpass2_max(iz) = maxval(snr_bpass2(:,iz))
   end do

end subroutine init_SNII_bpass_v2_300

