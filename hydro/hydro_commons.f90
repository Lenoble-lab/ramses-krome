module hydro_commons
  use amr_parameters
  use hydro_parameters
  real(udp),allocatable,dimension(:,:)::uold,unew ! State vector and its update
  real(udp),allocatable,dimension(:)::pstarold,pstarnew ! Stellar momentum and its update
  real(udp),allocatable,dimension(:)::divu, enew ! Non conservative variables
  real(dp)::mass_tot=0,mass_tot_0=0
  real(dp)::ana_xmi,ana_xma,ana_ymi,ana_yma,ana_zmi,ana_zma
  integer::nbins
#ifdef SOLVERmhd
  integer,allocatable,dimension(:)::liste_ind
  integer::nb_ind
  real(dp):: dt_imp                              ! Implicit timestep 
  integer :: niter=0                             ! Total number of iteration                 
  real(dp):: epsilon_diff=1d-6                   ! CG iteration break criteria                
  real(dp),allocatable,dimension(:)::apotx,apoty,apotz ! MHD ICs
#endif
end module hydro_commons

module const
  use amr_parameters
  real(dp)::bigreal = 1.0d+30
  real(dp)::zero = 0
  real(dp)::one = 1
  real(dp)::two = 2
  real(dp)::three = 3
  real(dp)::four = 4
  real(dp)::two3rd = 2/3d0
  real(dp)::half = 1/2d0
  real(dp)::third = 1/3d0
  real(dp)::forth = 1/4d0
  real(dp)::sixth = 1/6d0
end module const

