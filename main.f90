!> LEOPARD: Initializes the setup, scans through the requested wavenumber interval, computes corresponding frequencies, and prints dispersion relation to output file
! fort-kv-ints: calls subroutines from kv_ints_mod that compute eq. (1) in www.github.com/brentfpage/fort-kv-ints-r/blob/main/preprint.pdf , which is a nonlinear growth rate attributable to induced scattering 
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
  type(mp_real), allocatable, dimension(:) :: sol_at_kknots
  complex, allocatable, dimension (:) :: solution

  real, allocatable, dimension(:,:,:,:,:) :: splcoeff1, splcoeff2
  type(mp_real), allocatable, dimension(:,:,:,:,:) :: splcoeff4

  real :: start, finish
  real :: start2, finish2

  integer :: i_int
  integer :: ipara

  real, allocatable, dimension (:) :: Bksq

  type(mp_complex), allocatable, dimension (:,:,:) :: gam2_is_wccs

  type(mp_complex), allocatable, dimension(:,:,:) :: v_int_lam1_diff

  real :: disp_deriv
  external :: disp_deriv
  integer, dimension(416,15) :: all_int_params_standard
  integer, dimension(653,15) :: all_int_params_spank1
  complex :: rh_disp_val
  external rh_disp_val
  type(mp_complex) :: t1_root

  type(mp_real) :: k1, om1, vub, vlb
  type(mp_real), allocatable, dimension(:,:) :: om2splcoeffs_nk, pfdsplcoeffs_nk

  integer :: k_pow

  type(mp_real) :: mppic

  ! parameter negk for fort-kv-ints:
  ! if a distribution is symmetric w.r.t. B0->-B0, waves may grow propagating along both B0 and -B0,
  ! in which case the program should be run once with negk = .false. and once with negk = .true. and the results summed.
  ! if waves only propagate only along one direction, then the input distribution can be 
  ! configured such that this direction is +k, and negk should then be set to .false.
  logical, parameter :: negk = .false.
  integer :: spank1_ik2

  type(mp_real), allocatable, dimension(:,:) :: Ivpe
  type(mp_complex), allocatable, dimension(:,:,:,:,:) :: int_kq_roots_diff
  type(mp_complex), allocatable, dimension(:,:,:,:) :: int_kq_roots_diff_spank1
  type(mp_complex), allocatable, dimension(:,:) :: kroots_nk
  integer :: iperp

  mppic=mppi(kv_nwds)

  open(unit=7,status='unknown',file='omega.dat')

  call cpu_time(start)

  write(*,*) 'Read input data'
  call read_data(omega_start, increment, kstart, kend, nk)
  write(*,*) '...done'

  write(*,*) 'Read velocity distributions from files'
  call read_distr
  write(*,*) '...done.'

  ! the headers of gam2_is_wccs_for_ik and sum_gam2_is_wccs_over_ints_and_vperp in kv_ints_mod.f90 describe the contents of these integral
  ! parameter text files.
  write(*,*) 'Read integral params'
  open(unit=72,status='old',file='all_int_params_standard.txt')
  do i_int=1,size(all_int_params_standard, 1)
    read(72,*) all_int_params_standard(i_int,:)
  enddo
  close(72)
  ! in the spank1 interval, (k2-k1) can be factored out of the denominator, t₃ = ω₁–ω₂ - v(k₁-k₂) + i*eps, of the quantity
  ! g⁰_{₁K₊-₂K₊} appearing in the expression in the header of kv_ints_mod.f90.  Also, after a lot of algebra, (k2-k1) can be
  ! factored out of the inner bracketed term in that expression.  The two factors of (k2-k1) then cancel, which is important because
  ! otherwise the integral over t₃ would be undefined in some circumstances.  The process for generating the integral list
  ! all_int_params_spank1.txt, referenced below, included pulling out and cancelling these factors of (k2-k1).
  open(unit=72,status='old',file='all_int_params_spank1.txt')
  do i_int=1,size(all_int_params_spank1, 1)
    read(72,*) all_int_params_spank1(i_int,:)
  enddo
  close(72)
  write(*,*) '...done'
  
  allocate(krange(nk),solution(nk),Bksq(nk))
  allocate(om2splcoeffs_nk(3,nk),pfdsplcoeffs_nk(3,nk), kroots_nk(3,nk))
  allocate(kknots(nk-1))
  allocate(sol_at_kknots(nk-1))
  allocate(krange_mp(nk),solution_mp(nk))
  allocate(disp_derivs(nk))
  allocate(gam2_is_wccs(nk-1,nk-2,0:2))
  dk=(kend-kstart)/(nk-1.0)
  do ik=1,nk
     krange(ik)=kstart+(ik-1)*dk
  enddo

  allocate(splcoeff1(npara_max-1,nperp_max-1,4,3,narb))
  allocate(splcoeff2(npara_max-1,nperp_max-1,4,3,narb))
  allocate(splcoeff4(npara_max-1,nperp_max-1,f_spl_degr,f_spl_degr,narb))
  allocate(v_int_lam1_diff(0:n_max,0:lam1_max,npara_max))

  allocate(Ivpe(0:vperp_pow_max,nperp_max))
  allocate(int_kq_roots_diff(&
    nk-2,&
    q_minn:q_maxx_standard,&
    0:n_max+lam3_max,&
    0:sigma_max_standard,&
    0:sigma_max_standard))
  allocate(int_kq_roots_diff_spank1(&
    q_minn:q_maxx_spank1,&
    0:n_max+lam3_max,&
    0:sigma_max_spank1,&
    0:sigma_max_spank1))


  do iarb=1,narb
     call get_splinecoeff(iarb,splcoeff1(:,:,:,:,iarb),splcoeff2(:,:,:,:,iarb))
     call make_interp_spline_2d_nak_mp(iarb,splcoeff4(:,:,:,:,iarb),f_spl_degr)
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
  if(negk) then
    do ik2=1,nk-2
      om2splcoeffs_nk(2,ik2) = -om2splcoeffs_nk(2,ik2)
    enddo
    do ik2=1,nk-1
      kknots(ik2) = -kknots(ik2)
    enddo
    kknots(:) = kknots(size(kknots):1:-1)
    om2splcoeffs_nk(:,1:nk-2) = om2splcoeffs_nk(:,nk-2:1:-1)
  endif

  gam2_is_wccs = mpcmplx((0.0,0.0),kv_nwds)

  start = omp_get_wtime()
  write(*,*) ' '
  write(*,*) 'Compute induced scattering wave coupling coefficients'
 call set_eps(1.e-100)
 do ik=2,nk-1
  write(*,*) ' '
  write(*,'(A7,I6,A13,F12.8)') '-------',ik,'------- k₁=', krange(ik)
    start2 = omp_get_wtime()
    k1=krange_mp(ik)
    om1=solution_mp(ik)

    v_int_lam1_diff = mpcmplx(cmplx(0.0,0.0),kv_nwds)

! compute 
!  vb
! ⌠       vⁿ
! ⎮ dv ────────
! ⌡          λ₁
!      (v-vₒ)
! for each bound vb in the parallel velocity grid.  Here, vₒ=(-1+ω₁+i*eps)/k₁, i.e., vₒ is a zero of the denominator term t₁ = ω₁ - v*k₁ - 1 + i*eps .
    iarb=1 ! only narb=1 is supported at present
    do ipara=1,npara(iarb)-1
      vlb = mpreald(vpara(ipara,iarb),kv_nwds)
      vub = mpreald(vpara(ipara+1,iarb),kv_nwds)
      t1_root = (-1.0 + om1 + mpcmplx(i,kv_nwds) * eps)/k1 
      v_int_lam1_diff(:,:,ipara) = reshape(flatsubtract(&
        do_v_int(t1_root,vub),&
        do_v_int(t1_root,vlb),&
        size(v_int_lam1_diff(:,:,ipara))),&
        shape(v_int_lam1_diff(:,:,ipara)))
    enddo
    do iperp=1,nperp(iarb)
      do i_int=0,vperp_pow_max
        Ivpe(i_int,iperp) = mpreald(vperp(iperp,iarb),kv_nwds)**(i_int+1)/(i_int+1)
      enddo
    enddo

    do ik2=1,nk-2
      if((kknots(ik2).lt.k1).and.(kknots(ik2+1).gt.k1)) then 
        spank1_ik2 = ik2
        exit
      endif
    enddo

    pfdsplcoeffs_nk = compute_pfdsplcoeffs_nk(om2splcoeffs_nk, k1, om1, spank1_ik2)
    kroots_nk(2:3,:) = compute_pfd_kroots_nk(pfdsplcoeffs_nk, spank1_ik2)
    kroots_nk(1,:) = k1 * mpcmplx((1.0,0.0), kv_nwds)

    !$omp parallel do private(ik2)
    do ik2=1,nk-2
      if(ik2.ne.spank1_ik2) then
        int_kq_roots_diff(ik2,:,:,:,:) = mpcmplx((0.0,0.0),kv_nwds)
        int_kq_roots_diff(ik2,:,:,:,:) = integrate_over_v_independent_k_roots(kroots_nk(:,ik2),&
          kknots(ik2), kknots(ik2+1), .false.)
      else
        int_kq_roots_diff_spank1(:,:,:,:) = mpcmplx((0.0,0.0),kv_nwds)
        int_kq_roots_diff_spank1(:,:,:,:) = integrate_over_v_independent_k_roots(kroots_nk(:,ik2),&
          kknots(ik2), kknots(ik2+1), .true.)
      endif
    enddo
    !$omp end parallel do

    gam2_is_wccs(ik,:,:) = gam2_is_wccs_for_ik(int_kq_roots_diff, int_kq_roots_diff_spank1, v_int_lam1_diff, Ivpe,&
      om2splcoeffs_nk, pfdsplcoeffs_nk, splcoeff4, kroots_nk, kknots, vpara(:npara(iarb),iarb), om1, k1,&
      all_int_params_spank1, all_int_params_standard, spank1_ik2)

    do ik2=1,nk-2
      do k_pow=0,2
        gam2_is_wccs(ik,ik2,k_pow) = gam2_is_wccs(ik,ik2,k_pow) * mppic / mpreald(delta**2 * disp_derivs(ik),kv_nwds)
      enddo

      write(*,*) ' '
      write(*,'(A13,I6,A14,F12.8,A4,F12.8)') '-------------',ik2,'------- k₂ =', qreal(kknots(ik2)),' -->',qreal(kknots(ik2+1))
      write(*,'(A20)') 'S(k1,klb,kub,k_pow)'
      do k_pow=0,2
        write(*,'(A7,I2,A7,E27.17)') 'k_pow =',k_pow,':  S = ',qreal(aimag(gam2_is_wccs(ik,ik2,k_pow)))
      enddo
    enddo
  finish2 = omp_get_wtime()
  write(*,*) '     k₁ time elapsed:', finish2-start2
  enddo
  open(unit=7,status='unknown',file='gam2_is_wccs.dat')
  do ik=2,nk-1
    do ik2=1,nk-2
      do k_pow=0,2
        write(7,'(F12.8,F12.8,I5,E20.10)') krange(ik), qreal(kknots(ik2)) ,k_pow,&
          qreal(aimag(gam2_is_wccs(ik,ik2,k_pow)))
      enddo
    enddo
  enddo
  close(7)
  finish = omp_get_wtime()
  write(*,*) 'Total time elapsed:', finish-start
            
 
  deallocate(krange,solution,Bksq)
  deallocate(om2splcoeffs_nk)
  deallocate(kknots)
  deallocate(sol_at_kknots)
  deallocate(krange_mp,solution_mp)
  deallocate(disp_derivs)
  deallocate(gam2_is_wccs)
  deallocate(splcoeff1,splcoeff2)
  deallocate(splcoeff4)
  deallocate(v_int_lam1_diff)

  deallocate(mu,q)
  deallocate(beta_para,beta_perp,beta_ratio)
  deallocate(mode, dens, drift)
  deallocate(distribution,vpara,vperp,npara,nperp)

end program main
