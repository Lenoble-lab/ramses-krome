module rt_hydro_commons
  use amr_parameters
  use rt_parameters
  real(rtdp),allocatable,dimension(:,:)::rtuold,rtunew ! State vector and its update
end module rt_hydro_commons

