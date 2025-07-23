function disp_deriv(k,sol,splcoeff1,splcoeff2)
  use param_mod
  implicit none

  real :: k
  complex :: sol
  complex :: rh_disp_val
  complex :: disp_val1, disp_val2
  real :: disp_deriv
  real :: delta_sol

  real, dimension(npara_max-1,nperp_max-1,4,3,narb) :: splcoeff1, splcoeff2

  external rh_disp_val

  delta_sol = 0.001
  disp_val1 = rh_disp_val(sol - delta_sol,k,splcoeff1,splcoeff2)/((sol-delta_sol)*delta)**2
  disp_val2 = rh_disp_val(sol + delta_sol,k,splcoeff1,splcoeff2)/((sol+delta_sol)*delta)**2

  disp_deriv = real(disp_val2 - disp_val1)/(2*delta_sol)

end function disp_deriv




