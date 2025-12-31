!> Module containing all globally defined parameters
module param_mod
implicit none
complex, parameter :: i=(0.0,1.0)
real, parameter :: pi= 4*atan(1.0)
integer :: Nspecies
integer :: sign
real :: theta
integer, allocatable, dimension (:) :: mode
real, allocatable, dimension (:) :: mu
real, allocatable, dimension (:) :: q
real, allocatable, dimension (:) :: dens
real, allocatable, dimension (:) :: drift
real, allocatable, dimension (:) :: beta_para
real, allocatable, dimension (:):: beta_perp
real, allocatable, dimension (:) :: beta_ratio
real, allocatable, dimension (:,:,:) :: distribution
real, allocatable, dimension (:,:) :: vpara,vperp
integer, allocatable, dimension (:) :: npara, nperp
integer :: npara_max, nperp_max
integer :: narb
real :: delta 
real :: rf_error
real :: eps_error


! for fort-kv-ints
integer, parameter :: kv_nwds = 8 ! 108 digits; sets the working precision of fort-kv-ints
integer, parameter :: sigma_max_spank1 = 6
integer, parameter :: sigma_max_standard = 4
integer, parameter :: lam1_max = 4
integer, parameter :: lam2_max = 3
integer, parameter :: lam3_max = 2
integer, parameter :: q_min = -3
integer, parameter :: q_max = 7

integer, parameter :: q_minn = -7
integer, parameter :: q_maxx_spank1 = 22
integer, parameter :: q_maxx_standard = 19
integer, parameter :: n_max = 8
integer, parameter :: vperp_pow_max = 10
integer, parameter :: kv_root_lam_max = 3 
integer, parameter :: f_spl_degr = 6 ! quintic spline


end module param_mod
