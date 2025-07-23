! adapted from disp_det.f90 and integrator.f90 in LEOPARD
!> Computes the dispersion constant (= dialectric constant - (ck/omega)^2 ) of right-handed parallel-propagating electromagnetic waves for a plasma with bi-Maxwellian particle species and/or particles with arbitrary gyrotropic velocity distribution
!! \param omega complex wave frequency
!! \param k wavenumber
!! \param splcoeff1 array of coefficients obtained from cubic spline interpolation of the provided velocity distribution data 
!! \param splcoeff2 array of coefficients obtained from cubic spline interpolation of the provided velocity distribution data 
function rh_disp_val(omega,k,splcoeff1,splcoeff2)
  use param_mod
  implicit none

  complex :: rh_disp_val

  real, dimension(npara_max-1,nperp_max-1,4,3,narb) :: splcoeff1, splcoeff2

  complex :: omega
  real :: k
  integer :: m

  complex :: intgrl_a
  complex :: intgrl_b
  complex :: intgrl_c

  complex, dimension(3) :: sum0
  real, dimension(3) :: sum1

  real, allocatable, dimension(:,:) :: Ivpe
  complex, allocatable, dimension(:,:) :: Kvpa
  integer :: iperp, ipara
  integer :: iarb

  complex :: zeta1
  complex(kind=16) :: zeta

  complex :: Z_func
  external Z_func

  complex(kind=16), dimension (6) :: Kvpa_dummy
  integer :: n, ndp
  real :: zeta_r
  real :: dvpa

  n=-1

  !note that the dielectric tensor is renormalized according to epsilon_ij -> epsilon_ij * (v_A ^2/c^2 * omega^2)  
  rh_disp_val = 0.0
  rh_disp_val = rh_disp_val + (delta*omega)**2

  iarb=0

  do m=1,Nspecies

     if(mode(m).eq.0) then

        !for bi-Maxwellian scenario
        rh_disp_val=rh_disp_val+(beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2
        zeta1=(omega-k*drift(m)-n*mu(m)*q(m))/sqrt(beta_para(m))/k/sqrt(mu(m)) *sqrt(dens(m))
        rh_disp_val=rh_disp_val + Z_func(zeta1) * mu(m)*dens(m)*q(m)**2 * &
          ( zeta1*beta_ratio(m)-1.0/sqrt(beta_para(m))/k*q(m)*sqrt(dens(m)*mu(m)) )

     else if(mode(m).eq.1) then
!    
       !for particles with arbitrary gyrotropic velocity distributions
! 
        iarb=iarb+1
        allocate(Ivpe(5,nperp(iarb)))
        allocate(Kvpa(6,npara(iarb)-1))

   do iperp=1,nperp(iarb)

      Ivpe(1,iperp)=(vperp(iperp,iarb)**3)/3.0
      Ivpe(2,iperp)=(vperp(iperp,iarb)**4)/4.0
      Ivpe(3,iperp)=(vperp(iperp,iarb)**5)/5.0
      Ivpe(4,iperp)=(vperp(iperp,iarb)**6)/6.0
      Ivpe(5,iperp)=(vperp(iperp,iarb)**7)/7.0

   enddo

  zeta=(omega-n*mu(m)*q(m))/k
  dvpa=abs(vpara(2,iarb)-vpara(1,iarb))
  zeta_r = real(zeta)

  do ipara=1,npara(iarb)-1
    call acc_Kvpa(dvpa,zeta_r,ndp)
! the number of words ( mpwds = int (ndp / log10(2^60) + 2)  ) for mpfun numbers must be greater than 3 
    if(ndp.le.36) then 
      call int_para(iarb,ipara,k,zeta,Kvpa_dummy)
    else
      call int_para_mpfun(iarb,ipara,ndp,k,zeta,Kvpa_dummy)
    endif
    Kvpa(1,ipara)=Kvpa_dummy(1)
    Kvpa(2,ipara)=Kvpa_dummy(2)
    Kvpa(3,ipara)=Kvpa_dummy(3)
    Kvpa(4,ipara)=Kvpa_dummy(4)
    Kvpa(5,ipara)=Kvpa_dummy(5)
    Kvpa(6,ipara)=Kvpa_dummy(6)
  enddo

  Kvpa=-Kvpa/k

  intgrl_a = (0.0,0.0)
  intgrl_b = (0.0,0.0)
  intgrl_c = (0.0,0.0)

  ! lines below adapted from integrator.f90 in LEOPARD
  do iperp=1,nperp(iarb)-1

       sum0=(0.0,0.0)

       do ipara=1,npara(iarb)-1

          sum0=sum0+splcoeff1(ipara,iperp,1,:,iarb)*Kvpa(4,ipara)+splcoeff1(ipara,iperp,2,:,iarb)*Kvpa(3,ipara)+&
               &    splcoeff1(ipara,iperp,3,:,iarb)*Kvpa(2,ipara)+splcoeff1(ipara,iperp,4,:,iarb)*Kvpa(  1,ipara)

       enddo

       intgrl_a = intgrl_a+&
            & (Ivpe(3,iperp+1)-Ivpe(3,iperp))*sum0(1)+&
            & (Ivpe(2,iperp+1)-Ivpe(2,iperp))*sum0(2)+&
            & (Ivpe(1,iperp+1)-Ivpe(1,iperp))*sum0(3)

       sum0=(0.0,0.0)

       do ipara=1,npara(iarb)-1

          sum0=sum0+splcoeff1(ipara,iperp,1,:,iarb)*Kvpa(5,ipara)+splcoeff1(ipara,iperp,2,:,iarb)*Kvpa(4,ipara)+&
               &    splcoeff1(ipara,iperp,3,:,iarb)*Kvpa(3,ipara)+splcoeff1(ipara,iperp,4,:,iarb)*Kvpa(  2,ipara)

       enddo

       intgrl_b = intgrl_b+&
            & (Ivpe(3,iperp+1)-Ivpe(3,iperp))*sum0(1)+&
            & (Ivpe(2,iperp+1)-Ivpe(2,iperp))*sum0(2)+&
            & (Ivpe(1,iperp+1)-Ivpe(1,iperp))*sum0(3)
    enddo

    do ipara=1,npara(iarb)-1
       sum1=0.0
       do iperp=1,nperp(iarb)-1
         sum1=sum1+&
           &    splcoeff2(ipara,iperp,1,:,iarb)*(Ivpe(5,iperp+1)-Ivpe(5,iperp))+&
           &    splcoeff2(ipara,iperp,2,:,iarb)*(Ivpe(4,iperp+1)-Ivpe(4,iperp))+&
           &    splcoeff2(ipara,iperp,3,:,iarb)*(Ivpe(3,iperp+1)-Ivpe(3,iperp))+&
           &    splcoeff2(ipara,iperp,4,:,iarb)*(Ivpe(2,iperp+1)-Ivpe(2,iperp))
      enddo

      intgrl_c = intgrl_c+&
          & Kvpa(3,ipara)*sum1(1)+&
          & Kvpa(2,ipara)*sum1(2)+&
          & Kvpa(1,ipara)*sum1(3)
    enddo

    rh_disp_val=rh_disp_val + pi * mu(m) * q(m)**2 * dens(m) * &
        & (  omega*intgrl_a-&
        &    k*intgrl_b+&
        &    k*intgrl_c)

    deallocate(Ivpe)
    deallocate(Kvpa)

     else
        write(*,*) 'Check dist_in!'
        stop
     endif

  enddo

  rh_disp_val = rh_disp_val - k**2

end function rh_disp_val
