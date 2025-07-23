!> LEOPARD: Initializes the setup, scans through the requested wavenumber interval, computes corresponding frequencies, and prints dispersion relation to output file
! fort-kv-ints: does part of the computation of the general form of the weak turbulence integral ( eq. (2) in @@ ) , and sums such integrals up into the induced scattering nonlinear growth rate (eq. (1) in @@)
program main
  use param_mod
  use kv_ints_mod
  use pppack_mod_mp, only : make_interp_spline_quad_mp, make_interp_spline_2d_nak_mp
  use omp_lib
  implicit none
  complex :: omega_start, increment
  integer :: nk,ik,ik2, iarb

  real :: kstart, kend, dk
  real, allocatable, dimension (:) :: disp_derivs
  real, allocatable, dimension (:) :: krange
  type(mp_real), allocatable, dimension(:) :: kknots
  type(mp_real), allocatable, dimension(:) :: krange_mp, solution_mp
  complex, allocatable, dimension (:) :: solution

  real, allocatable, dimension(:,:,:,:,:) :: splcoeff1, splcoeff2
  type(mp_real), allocatable, dimension(:,:,:,:,:) :: splcoeff4

  real :: start, finish
  real :: start2, finish2

  integer :: i_int
  integer :: ipara

  real, allocatable, dimension (:) :: Bksq

  type(mp_complex), allocatable, dimension (:,:,:) :: gam2_is

  type(mp_complex), allocatable, dimension(:,:,:,:) :: v_int_vub_lam1
  type(mp_complex), allocatable, dimension(:,:,:,:) :: v_int_vlb_lam1

  real :: disp_deriv
  external :: disp_deriv
  integer, dimension(416,15) :: all_int_params
  integer :: k_pow
  complex :: rh_disp_val
  external rh_disp_val
  type(mp_complex) :: t1_root

  type(mp_real) :: klb, kub, k1, om1, vub, vlb
  type(mp_real), dimension(3) :: om2splcoeffs ! pfdsplcoeffs
  type(mp_real), allocatable, dimension(:,:) :: om2splcoeffs_nk

  real :: tp1, tp2
  integer :: kpow

  type(mp_real) :: mppic

  ! parameter negk for fort-kv-ints:
  ! if a distribution is symmetric w.r.t. B0->-B0, waves may grow propagating along both B0 and -B0,
  ! in which case the program should be run once with negk = .false. and once with negk = .true. and the results summed.
  ! if waves only propagate only along one direction, then the input distribution can be 
  ! configured such that this direction is +k, and negk should then be set to .false.
  logical, parameter :: negk = .false.

  mppic=mppi(kv_nwds)


  open(unit=7,status='unknown',file='omega.dat')

  call cpu_time(start)

  write(*,*) 'Read input data'
  call read_data(omega_start, increment, kstart, kend, nk)
  write(*,*) '...done'

  write(*,*) 'Read velocity distributions from files'
  call read_distr
  write(*,*) '...done.'


  write(*,*) 'Read integral params'
  open(unit=72,status='old',file='all_int_params.txt')
  do i_int=1,size(all_int_params, 1)
    read(72,*) all_int_params(i_int,:)
  enddo
  close(72)
  write(*,*) '...done'
  
  allocate(krange(nk),solution(nk),Bksq(nk))
  allocate(om2splcoeffs_nk(3,nk))
  allocate(kknots(nk-1))
  allocate(krange_mp(nk),solution_mp(nk))
  allocate(disp_derivs(nk))
  allocate(gam2_is(nk-1,nk,0:3))
  dk=(kend-kstart)/(nk-1.0)
  do ik=1,nk
     krange(ik)=kstart+(ik-1)*dk
  enddo

  allocate(splcoeff1(npara_max-1,nperp_max-1,4,3,narb))
  allocate(splcoeff2(npara_max-1,nperp_max-1,4,3,narb))
  allocate(splcoeff4(npara_max-1,nperp_max-1,6,6,narb))
  allocate(v_int_vub_lam1(0:n_max,0:lam1_max,npara_max,narb))
  allocate(v_int_vlb_lam1(0:n_max,0:lam1_max,npara_max,narb))

  do iarb=1,narb
     call get_splinecoeff(iarb,splcoeff1(:,:,:,:,iarb),splcoeff2(:,:,:,:,iarb))
     call make_interp_spline_2d_nak_mp(iarb,splcoeff4(:,:,:,:,iarb),6)
  enddo


  write(*,*) 'Compute the dispersion relation'
  !scan through wavenumber interval
  do ik=1,nk

     write(*,*) ' '
     write(*,'(A7,I6,A10,F12.8)') '-------',ik,'------- k=', krange(ik)

     call cpu_time(start2)

     !use Muller method to iterate root of dispersion relation
     call muller(omega_start,krange(ik),solution(ik),splcoeff1,splcoeff2)
     disp_derivs(ik) = disp_deriv(krange(ik),solution(ik),splcoeff1,splcoeff2)

     call cpu_time(finish2)

     if(aimag(solution(ik)).gt.0.0) then
       Bksq(ik)=1e-6*aimag(solution(ik))
    else 
      Bksq(ik)=0.0
    endif

     write(*,'(A9,E20.10,A9,E20.10,A13,E20.10)')  &
       '   omega:', real(solution(ik)), &
       '   gamma:',aimag(solution(ik)), &
       '  derivative:',disp_derivs(ik)
     write(*,*) 'time elapsed:', finish2-start2


     if ((ik .ge. 3).and.(ik .lt. nk))  then

        !if three subsequent solutions omega(k) are found, use quadratic polynomial fit 
        !to guess next starting frequency for Muller iteration
        call polyfit(krange(ik-2:ik+1),solution(ik-2:ik),omega_start)

     else

        !for the first two solution omega(k) guess next starting frequency for Muller iteration
        !by raising the computed omega by an increment which is provided by the user
        omega_start=solution(ik)+increment

     end if
     write(7,'(F12.8,E20.10,E20.10,E20.10)') krange(ik), real(solution(ik)), aimag(solution(ik)),disp_derivs(ik)

  enddo
  close(7)

  call cpu_time(finish)
 
  write(*,*) 'Total time elapsed:', finish-start

  do ik=1,nk
    krange_mp(ik) = mpreald(krange(ik),kv_nwds)
    solution_mp(ik) = mpreald(real(solution(ik)),kv_nwds)
  enddo
  call make_interp_spline_quad_mp(krange_mp, solution_mp, om2splcoeffs_nk, kknots)

  do iarb=1,narb
  enddo

  gam2_is = mpcmplx((0.0,0.0),kv_nwds)

  start = omp_get_wtime()
  write(*,*) ' '
  write(*,*) 'Compute induced scattering wave coupling coefficients'
 call set_eps(1e-100)
 do ik=2,nk-1
  write(*,*) ' '
  write(*,'(A7,I6,A13,F12.8)') '-------',ik,'------- k₁=', krange(ik)
    start2 = omp_get_wtime()
    k1=krange_mp(ik)
    om1=solution_mp(ik)

    v_int_vlb_lam1 = mpcmplx(cmplx(0.0,0.0),kv_nwds)
    v_int_vub_lam1 = mpcmplx(cmplx(0.0,0.0),kv_nwds)


! compute 
!  vb
! ⌠       vⁿ
! ⎮ dv ────────
! ⌡          λ₁
!      (v-vₒ)
! for each bound vb in the parallel velocity grid.  Here, vₒ=(-1+ω₁+i*eps)/k₁, i.e., vₒ is a zero of the denominator term t₁ = ω₁ - v*k₁ - 1 + i*eps .
    do iarb=1,narb
      do ipara=1,npara(iarb)-1
        vlb = mpreald(vpara(ipara,iarb),kv_nwds)
        vub = mpreald(vpara(ipara+1,iarb),kv_nwds)
        t1_root = (-1.0 + om1 + mpcmplx(i,kv_nwds) * eps)/k1 
        call do_v_int(t1_root, vub, v_int_vub_lam1(:,:,ipara,iarb))
        call do_v_int(t1_root, vlb, v_int_vlb_lam1(:,:,ipara,iarb))
      enddo
    enddo

    !$omp parallel do private(ik2,klb,kub,om2splcoeffs,kpow)&
    !$omp shared(gam2_is)
    do ik2=1,nk-2

      if(negk) then
        klb=-kknots(ik2+1)
        kub=-kknots(ik2)
      else
        klb=kknots(ik2)
        kub=kknots(ik2+1)
      endif

      om2splcoeffs(:) = om2splcoeffs_nk(:,ik2)
      if(negk) then
        om2splcoeffs(2) = -om2splcoeffs(2)
      endif

      call gam2_is_for_ik2(om1,k1,splcoeff4,om2splcoeffs,all_int_params,&
        gam2_is(ik,ik2,:),klb,kub,v_int_vub_lam1,v_int_vlb_lam1,disp_derivs(ik))

    write(*,*) ' '
    write(*,'(A13,I6,A14,F12.8,A4,F12.8)') '-------------',ik2,'------- k₂ =', qreal(klb),' -->',qreal(kub)
    write(*,'(A18)') 'S(k1,klb,kub,kpow)'
    do kpow=0,2
      tp1 = aimag(gam2_is(ik,ik2,kpow))
      write(*,'(A7,I2,A7,E20.10)') 'kpow =',kpow,':  S = ',tp1
    enddo
  enddo
  !$omp end parallel do
  finish2 = omp_get_wtime()
  write(*,*) '     k₁ time elapsed:', finish2-start2
  enddo
  open(unit=7,status='unknown',file='omega2.dat')
  do ik=2,nk-1
    do ik2=1,nk-2
      do k_pow=0,3
        tp1 = aimag(gam2_is(ik,ik2,k_pow))
        tp2 = kknots(ik2)
        write(7,'(F12.8,F12.8,I5,E20.10)') krange(ik), tp2,k_pow, tp1
      enddo
    enddo
  enddo
  close(7)
  finish = omp_get_wtime()
  write(*,*) 'Total time elapsed:', finish-start
            
 
  deallocate(krange,solution,Bksq)
  deallocate(om2splcoeffs_nk)
  deallocate(kknots)
  deallocate(krange_mp,solution_mp)
  deallocate(disp_derivs)
  deallocate(gam2_is)
  deallocate(splcoeff1,splcoeff2)
  deallocate(splcoeff4)
  deallocate(v_int_vub_lam1,v_int_vlb_lam1)

  deallocate(mu,q)
  deallocate(beta_para,beta_perp,beta_ratio)
  deallocate(mode, dens, drift)
  deallocate(distribution,vpara,vperp,npara,nperp)

end program main
