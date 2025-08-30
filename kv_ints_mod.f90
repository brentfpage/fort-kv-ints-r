! fort-kv-ints : a program for computing induced scattering wave coupling coefficients S(k₁,klb,kub,kpow) = 
!               ⎡  kub                                                                                                      ⎤
!       π       ⎢ ⌠      (kpow)⎛ω₂⎞²⌠                                                                                       ⎥
! – ──────────Im⎢i⎮ dk₂ k₂     ⎜──⎟ ⎮ dv dv  v² g⁻¹  T    g⁰      ⎡⎛T    - D    ⎞g⁻¹ T    + ⎛T   - D   ⎞g⁻¹  T    ⎤F(v ,v  )⎥ .
!       +  ,+   ⎣ ⌡            ⎝k₂⎠ ⌡   ⟂  || ⟂  ₁K₊  ₂K₊  ₁K₊-₂K₊⎣⎝ -₂K₊   -₂K₊⎠ ₁K₊ ₁K₊   ⎝ ₁K₊   ₁K₊⎠ -₂K₊ -₂K₊⎦   ⟂  || ⎦
!   δ² ω₁ Λ        klb
!          ₁K₊
! Here,
!        ⎛ω-kv          kv       ⎞
!        ⎜    ||  d       ⟂   d  ⎟
! T  = -i⎜────── ──── + ─── ─────⎟ ,
!  K     ⎜  ω    d v     ω  d v  ⎟
!        ⎝          ⟂          ||⎠
!
!       ω-kv
!           || 1
! D  = i────── ── ,
!  K      ω    v
!               ⟂
! 
!  M         1
! g  = ────────────── ,
!  K   ω-kv  +M+i*eps
!          ||
!                           +           +
! ₁K₊ is shorthand for (k₁,ω₁) , where ω₁ is the frequency of waves with wavenumber k₁ ,
!                           +           +
! ₂K₊ is shorthand for (k₂,ω₂) , where ω₂ is the frequency of waves with wavenumber k₂ , and
! ω₂ is a quadratic function of k₂ in the considered interval from k₂=klb to k₂=kub,
! ω₂ = om2splcoeffs(1) + om2splcoeffs(2) * k₂ + om2splcoeffs(3) * pow(k₂, 2).  Also, 
! the velocity distribution function F(v  , v ) is a quintic spline function of 
!                                       ||   ⟂                                  
! v   (parallel velocity) and v (perpendicular velocity)  ,
!  ||                          ⟂
!              5   5
!              ⎲   ⎲                                       i   j
! F(v  ,v ) =  ⎳   ⎳  splcoeff4(ipara,iperp,6-i,6-j,iarb) v   v   ,
!    ||  ⟂    i=0 j=0                                      ||  ⟂
! where the indices ipara and iperp depend on what grid intervals the evaluation points v   and v  lie in.  
!                                                                                        ||      ⟂
! The index iarb of splcoeff4 labels the species.  At present, only narb=1 is implemented for fort-kv-ints, so iarb should always be 1 .
! Further, δ is a program parameter specified by the user,
!  ,+                       
! Λ    = ⎛d Re ⎛Λ⁺     ⎞ ╱   ⎞ ⎢       ,
!  ₁K₊   ⎝     ⎝ (k₁,ω)⎠╱ dω ⎠ ⎢   +
!                              ⎢ω=ω₁
! and Λ⁺ is a quantity related to the dialectric constant described in www.github.com/brentfpage/fort-kv-ints-r/blob/main/preprint.pdf .
!      K
!  ,+ 
! Λ    has been calculated in main.f90 and is 
!  ₁K₊
! provided as the argument disp_deriv_ik to the subroutine gam2_is_for_ik2 below.
! the variables k and ω from www.github.com/brentfpage/fort-kv-ints-r/blob/main/preprint.pdf are referred to in this program as k₁ and ω₁, while the variables k' and ω' from that writeup are referred to as k₂ and ω₂.

module kv_ints_mod
    use param_mod, only : i, kv_nwds, sigma_max, lam1_max, lam2_max, lam3_max, q_min, q_max,&
      nperp, narb, npara_max, vpara, vperp, npara, nperp_max, q_minn, q_maxx, n_max, delta
    use mpmodule
    implicit none
    type(mp_real) :: eps ! small imaginary part for resonant denominators

! wrapper to accommodate log of negative real numbers
    interface logw
      module procedure mp_logw
      module procedure mp_clogw
    end interface

    interface combine_multi_roots_ints
      module procedure combine4_multi_roots_ints
      module procedure combine5_multi_roots_ints
      module procedure combine6_multi_roots_ints
    end interface

    interface do_rational_ints
        module procedure do5_rational_ints
        module procedure do6_rational_ints
    end interface
! 
    contains

! Computer algebra software has been used to expand the integral described in the module header into many sub-integrals, the parameters of which are collected in all_int_params.  The wavenumber and parallel velocity part of every sub-integral has the same general form and is computed by the subroutine int_driver.  For a given k₂ interval k₂=klb->kub , the subroutine below loops over the parallel velocity grid and for each relevant interval calls int_driver.  Then, for each given parallel velocity interval, sum_gam2_is_vperp is called, which loops over the list of integrals in all_int_params as well as over perpendicular velocity and at each loop iteration increments the result for the complete integral, gam2_is_ik_ik2.

! om1 : ω₁
! k1 : k₁
! v_int_v(l/u)b_lam1 : include some helper integrals relevant for all k₂ that are done in main.f90

subroutine gam2_is_for_ik2(om1, k1, splcoeff4, om2splcoeffs, all_int_params, &
  gam2_is_ik_ik2,klb,kub,v_int_vub_lam1,v_int_vlb_lam1,disp_deriv_ik)
  implicit none
  type(mp_real) :: om1, k1
  type(mp_real), dimension(:,:,:,:,:) :: splcoeff4
  type(mp_real), dimension(:) :: om2splcoeffs
  integer, dimension(:,:) :: all_int_params
  type(mp_complex), dimension(0:) :: gam2_is_ik_ik2

  type(mp_real) :: klb, kub
  real :: disp_deriv_ik
  type(mp_complex), dimension(0:,0:,:,:), intent(in) :: v_int_vub_lam1, v_int_vlb_lam1
  logical :: spank1

  type(mp_real) :: vlb, vub
  integer :: iarb, ipara, iperp
  integer, dimension(npara_max) :: vres_idxs
  integer :: n_vres, iv, i_int
  type(mp_complex), dimension(3) :: kroots

  type(mp_real), allocatable, dimension(:,:) :: Ivpe
  type(mp_complex), allocatable, dimension(:,:,:,:,:,:,:) :: t123_int ! stores integrals over k₂ and v_parallel
  type(mp_real), dimension(3) :: pfdsplcoeffs
  type(mp_real) :: mppic

  integer :: k_pow
  type(mp_real) :: to_mult

  integer, parameter :: kv_root_lam_max = 3 
  integer, parameter :: vperp_pow_max = 10

  type(mp_complex), allocatable, dimension(:,:,:,:,:,:) :: int_kq_roots_klb, int_kq_roots_kub

  allocate(int_kq_roots_klb(q_minn:q_maxx,0:n_max+lam3_max,0:sigma_max,0:sigma_max,0:kv_root_lam_max,0:kv_root_lam_max))
  allocate(int_kq_roots_kub(q_minn:q_maxx,0:n_max+lam3_max,0:sigma_max,0:sigma_max,0:kv_root_lam_max,0:kv_root_lam_max))
  allocate(t123_int(q_minn:q_maxx,0:n_max+lam3_max,0:sigma_max,0:n_max,0:lam1_max,0:lam2_max,0:lam3_max))
  allocate(Ivpe(0:vperp_pow_max,nperp_max))

  mppic = mppi(kv_nwds)

  if((klb.lt.k1).and.(kub.gt.k1)) then 
    spank1 = .true.
  else
    spank1 = .false.
  endif

  kroots = mpcmplx((0.0,0.0),kv_nwds)
  kroots(1) = k1
  ! computing helper integrals over k₂
  call k_int_driver(kub,klb,kroots,int_kq_roots_kub,int_kq_roots_klb,k1,om1,om2splcoeffs,pfdsplcoeffs,spank1)

  do iarb=1,narb
    do iperp=1,nperp(iarb)
      do i_int=0,vperp_pow_max
        Ivpe(i_int,iperp) = mpreald(vperp(iperp,iarb),kv_nwds)**(i_int+1)/(i_int+1)
      enddo
    enddo

    vres_idxs = 0
    n_vres = 0
    ! only parallel velocity intervals that include a resonance for the given k₂ interval yield a non-zero contribution to the integral in the module header.  the subroutine is_resonant_vec finds these parallel velocity intervals
    call is_resonant_vec(om2splcoeffs,om1,k1,klb,kub,vpara(:npara(iarb),iarb),vres_idxs,n_vres)

    do iv=1,n_vres
      ipara = vres_idxs(iv)
      t123_int = mpcmplx((0.0,0.0),kv_nwds)

      vlb=mpreald(vpara(ipara,iarb),kv_nwds)
      vub=mpreald(vpara(ipara+1,iarb),kv_nwds)

      call int_driver(kub,klb,kroots,t123_int,int_kq_roots_kub,int_kq_roots_klb,&
        v_int_vub_lam1(:,:,ipara,iarb),v_int_vlb_lam1(:,:,ipara,iarb),vlb,vub,k1,om1,om2splcoeffs,&
        pfdsplcoeffs,spank1)

      call sum_gam2_is_vperp(om1,k1, splcoeff4(ipara,:,:,:,iarb),om2splcoeffs,&
        all_int_params,Ivpe, t123_int, gam2_is_ik_ik2,iarb)
    enddo
  enddo

  to_mult = mppic / mpreald(delta**2 * disp_deriv_ik,kv_nwds)
  do k_pow=0,3
    gam2_is_ik_ik2(k_pow) = gam2_is_ik_ik2(k_pow) * to_mult 
  enddo

  deallocate(int_kq_roots_klb, int_kq_roots_kub,t123_int,Ivpe)
end subroutine gam2_is_for_ik2


!
! computes
! t123_int(q,τ,σ,n,λ₁, λ₂, λ₃) = 
!  kub       q+λ₂       λ₃
! ⌠         k₂   (k₂-k₁)      ⌠       vⁿ
! ⎮ dk₂ ───────────────────── ⎮ dv ───────   .
! ⌡           τ      σ      σ ⌡     ₃
!  klb  (k-₁k) (k-₂k) (k-₃k)       ┬─┬  λₐ
!                                  │ │ tₐ
!                                  ᵃ⁼¹
! Here, the velocity v is specifically parallel velocity
! Also,
! t₁ = ω₁ - v*k₁ - 1 + i*eps
! t₂ = –ω₂ + k₂v + 1 + i*eps
! t₃ = ω₁–ω₂ - v(k₁-k₂) + i*eps
! Further, ₘk = kroots(m)

subroutine int_driver(kub,klb,kroots,t123_int,int_kq_roots_kub,int_kq_roots_klb,&
    v_int_vub_lam1_ipara,v_int_vlb_lam1_ipara,vlb,vub,k1,om1,om2splcoeffs,pfdsplcoeffs,spank1)
    implicit none
    type(mp_complex), dimension(:) :: kroots
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:) :: int_kq_roots_kub, int_kq_roots_klb
    type(mp_complex), dimension(0:,0:), intent(in) :: v_int_vub_lam1_ipara, v_int_vlb_lam1_ipara
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:,0:) :: t123_int
    type(mp_real) :: k1, om1, kub, klb, vlb, vub
    type(mp_real), dimension(:) :: om2splcoeffs,pfdsplcoeffs
    logical :: spank1

    integer :: lam1, lam2, lam3, qq, ip, tau, sig, n

    integer :: t13_lam1_max, t13_lam3_max
    integer :: t123_lam1_max, t123_lam2_max, t123_lam3_max
    logical :: necessary
    integer, dimension(4,3) :: hardest_ints

    ! does the integral in the subroutine header for (λ₁>0, λ₂=0, λ₃=0)
    call assemble_t1_ints(int_kq_roots_klb,int_kq_roots_kub,&
        v_int_vlb_lam1_ipara,v_int_vub_lam1_ipara,&
        t123_int(:,:,:,:,:,0,0),spank1)

    ! does the integral in the subroutine header for (λ₁=0, λ₂>0, λ₃=0)
    call t23_int_driver(kub,klb,kroots,t123_int(:,:,:,:,0,0,:),&
      & int_kq_roots_kub,int_kq_roots_klb,vlb,vub,k1,om1,om2splcoeffs, .true.,spank1)

    ! does the integral in the subroutine header for (λ₁=0, λ₂=0, λ₃>0)
    call t23_int_driver(kub,klb,kroots,t123_int(:,:,:,:,0,:,0),&
      & int_kq_roots_kub,int_kq_roots_klb,vlb,vub,k1,om1,om2splcoeffs, .false.,spank1)


    t13_lam1_max = 4
    t13_lam3_max = 2

    t123_lam1_max = 4
    t123_lam2_max = 3
    t123_lam3_max = 2

    hardest_ints(1,:) = (/ 1, 2, 2 /)
    hardest_ints(2,:) = (/ 1, 3, 1 /)
    hardest_ints(3,:) = (/ 3, 0, 2 /)
    hardest_ints(4,:) = (/ 4, 0, 1 /)

! with the integral in the subroutine header having been computed for (λ₁>0, λ₂=0, λ₃=0), (λ₁=0, λ₂>0, λ₃=0), and (λ₁=0, λ₂=0, λ₃>0),
! the code below computes the general (λ₁>0, λ₂>0, λ₃>0) cases.

! the quadratic polynomial h(k)≡pfdsplcoeffs(ubound(pfdsplcoeffs,1)) * (k-₂k) * (k-₃k) gets introduced by partial fraction decomposition operations applied to 1/(t₁ t₂) , 1/(t₂ t₃) , and 1/(t₃ t₁) . This specific expression for h(k) applies for the non-spank1 case.  In the spank1 case, k₁ is a root of h(k)=pfdsplcoeffs(ubound(pfdsplcoeffs,1)-1) * (k-k₁) * (k-₂k), which introduces the need to use if(spank1) ... else ... blocks extensively, 

    do lam1=1,t123_lam1_max
      do lam2=1,t123_lam2_max

        necessary = .false.
        do ip=1,size(hardest_ints,1)
          if ((lam1.le.hardest_ints(ip,1)).and.(lam2.le.hardest_ints(ip,2))) then
              necessary = .true.
              exit
          endif
        enddo

        if (necessary) then
          if (spank1) then
            do qq=q_min,ubound(t123_int,1) - (lam1+lam2-1) 
              do tau=0,ubound(t123_int,2)-1
                do sig=0,sigma_max - (lam1+lam2-1) 
                do n=0,ubound(t123_int,4)
                  t123_int(qq,tau,sig,n,lam1,lam2,0) = k1/pfdsplcoeffs(ubound(pfdsplcoeffs,1)-1) * &
                      & (t123_int(qq+1,tau+1,sig+1,n,lam1,lam2-1,0) - &
                      &  t123_int(qq+1,tau+1,sig+1,n,lam1-1,lam2,0))
                    enddo
                enddo
              enddo
            enddo
          else
            do qq=q_min,ubound(t123_int,1) - (lam1+lam2-1)
              do tau=0,ubound(t123_int,2)
                do sig=0,sigma_max - (lam1+lam2-1)
                  do n=0,ubound(t123_int,4)
                    t123_int(qq,tau,sig,n,lam1,lam2,0)=k1/pfdsplcoeffs(ubound(pfdsplcoeffs,1)) * &
                        & (t123_int(qq+1,tau,sig+1,n,lam1,lam2-1,0) - &
                        &  t123_int(qq+1,tau,sig+1,n,lam1-1,lam2,0))
                  enddo
                enddo
              enddo
            enddo
          endif
        endif
      enddo
    enddo

    do lam1=1,t13_lam1_max
      do lam3=1,t13_lam3_max

        necessary = .false.
        do ip=1,size(hardest_ints,1)
        ! these cover all the cases needed for t123 integrals as well
          if ((lam1.le.hardest_ints(ip,1)).and.(lam3.le.hardest_ints(ip,3))) then
              necessary = .true.
              exit
          endif
        enddo

        if (necessary) then
          if(spank1) then  
            do qq=q_min,ubound(t123_int,1) - ((lam1+lam3)-1)
              do tau=0,ubound(t123_int,2)
                do sig=0,sigma_max - (lam1 + lam3 - 1)
                  do n=0,ubound(t123_int,4)
                  t123_int(qq,tau,sig,n,lam1,0,lam3) = &
                      & k1*(t123_int(qq,tau,sig+1,n,lam1,0,lam3-1) - &
                      &  t123_int(qq,tau,sig+1,n,lam1-1,0,lam3))/pfdsplcoeffs(ubound(pfdsplcoeffs,1)-1)
                  enddo
                enddo
              enddo
            enddo
          else
            do qq=q_min,ubound(t123_int,1) - (lam1+lam3-1)
              do sig=0,sigma_max - (lam1+lam3-1)
                do n=0,ubound(t123_int,4)

                  t123_int(qq,0,sig,n,lam1,0,lam3) = (&
                      k1 * t123_int(qq+1,0,sig+1,n,lam1,0,lam3-1) - &
                      k1 * t123_int(qq+1,0,sig+1,n,lam1-1,0,lam3) - &
                      k1**2 * t123_int(qq,0,sig+1,n,lam1,0,lam3-1) + &
                      k1**2 * t123_int(qq,0,sig+1,n,lam1-1,0,lam3))/pfdsplcoeffs(ubound(pfdsplcoeffs,1)) 
                  do tau=1,lam3_max
                      t123_int(qq,tau,sig,n,lam1,0,lam3) = (&
                          k1 * t123_int(qq,tau-1,sig+1,n,lam1,0,lam3-1) - &
                          k1 * t123_int(qq,tau-1,sig+1,n,lam1-1,0,lam3))/pfdsplcoeffs(ubound(pfdsplcoeffs,1)) 
                  enddo
                enddo
              enddo
            enddo
          endif
        endif
      enddo 
    enddo

    do lam2=1,t123_lam2_max
      do lam3=1,t123_lam3_max

        necessary = .false.
        do ip=1,size(hardest_ints,1)
          if ((lam2.le.hardest_ints(ip,2)).and.(lam3.le.hardest_ints(ip,3))) then
              necessary = .true.
              exit
          endif
        enddo

        if (necessary) then
          if(spank1) then
            do qq=q_min,ubound(t123_int,1) - (2*(lam2+lam3) - 1)
              do tau=0,ubound(t123_int,2)
                do sig=0,sigma_max - (lam2 + lam3 - 1) 
                  do n=0,ubound(t123_int,4)
                    do lam1=0,ubound(t123_int,5)
                      t123_int(qq,tau,sig,n,lam1,lam2,lam3) = &
                          & (t123_int(qq+1,tau,sig+1,n,lam1,lam2,lam3-1) - &
                          & t123_int(qq+1,tau,sig+1,n,lam1,lam2-1,lam3))/pfdsplcoeffs(ubound(pfdsplcoeffs,1)-1)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          else
            do qq=q_min,ubound(t123_int,1) - (2*(lam2+lam3) - 1)
              do sig=0,sigma_max - (lam2 + lam3 - 1)
                do n=0,ubound(t123_int,4)
                  do lam1=0,ubound(t123_int,5)
                    t123_int(qq,0,sig,n,lam1,lam2,lam3) = &
                        & (&
                        & t123_int(qq+2,0,sig+1,n,lam1,lam2,lam3-1) - &
                        & t123_int(qq+2,0,sig+1,n,lam1,lam2-1,lam3) - &
                            & k1 * t123_int(qq+1,0,sig+1,n,lam1,lam2,lam3-1) +&
                        & k1 * t123_int(qq+1,0,sig+1,n,lam1,lam2-1,lam3))/pfdsplcoeffs(ubound(pfdsplcoeffs,1))
                    do tau=1,lam3_max
                      t123_int(qq,tau,sig,n,lam1,lam2,lam3) = &
                          & (&
                          & t123_int(qq+1,tau-1,sig+1,n,lam1,lam2,lam3-1) - &
                          & t123_int(qq+1,tau-1,sig+1,n,lam1,lam2-1,lam3) &
                          & )/pfdsplcoeffs(ubound(pfdsplcoeffs,1))
                    enddo
                  enddo
                enddo
              enddo
            enddo
          endif
        endif
      enddo
    enddo


    do qq=q_min,q_max 
      do tau=0,ubound(t123_int,2)
        do sig=0,sigma_max 
          do n=0,ubound(t123_int,4)
            do lam1=0,lam1_max
              do lam2=0,ubound(t123_int,6)
                do lam3=0,ubound(t123_int,7)
                  t123_int(qq,tau,sig,n,lam1,lam2,lam3) = t123_int(qq,tau,sig,n,lam1,lam2,lam3)/(-k1)**lam1
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
end subroutine int_driver


! with
! int_kq_roots(q,p1,...,pₘ,...) =
!               q
! ⌠kb          k
! ⎮   ─────────────────── f(k) dk
! ⌡   (k-₁k)ᵖ¹...(k-ₘk)ᵖᵐ
! having already been computed for q=qₘᵢₙ...qₘₐₓ and p₁=0...max(p₁), ...,pₘ=0...max(pₘ),
! as well as
! int_kq_roots(q,0,...,0,pₘ₊₁,0,...,0) =
!         q
! ⌠kb    k
! ⎮   ────────        f(k) dk
! ⌡   (k-₍ₘ₊₁₎k)ᵖ⁽ᵐ⁺¹⁾

! having already been computed for q=qₘᵢₙ...qₘₐₓ and p₍ₘ₊₁₎=0...max(p₍ₘ₊₁₎),
! append_root computes
! int_kq_roots(q,p1,...,pₘ,pₘ₊₁,...) =
!                      q
! ⌠kb                 k
! ⎮   ─────────────────────────────────   f(k) dk
! ⌡   (k-₁k)ᵖ¹...(k-ₘk)ᵖᵐ(k-₍ₘ₊₁₎k)ᵖ⁽ᵐ⁺¹⁾
!
! for q=qₘᵢₙ...qₘₐₓ and p₁=0...max(p₁), ..., p₍ₘ₊₁₎=0...max(p₍ₘ₊₁₎)
! Here, ₘk = kroots(m) .
!                                                               ⎛┬─┬      ⎞
! Also, f(k) is any function. In practice, f(k) = 1 or f(k) = ln⎜│ │(k-kᵣ)⎟ for this program
!                                                               ⎝ ʳ       ⎠

    subroutine append_root(int_kq_roots,n,inshape,kroots,newdim)
        implicit none
        type(mp_complex), dimension(n) :: int_kq_roots
        type(mp_complex), dimension(:) :: kroots
        integer, dimension(:) :: inshape
        integer :: id, is, n, slicemax, uptostride, newstride
        integer :: uptodim,newdim
        integer :: ui, ni
        slicemax=1
!         slice all dimensions up to uptodim
        do uptodim=2,newdim-1
            slicemax=slicemax*inshape(uptodim-1)
            uptostride = slicemax ! increments the index of dimension uptodim
            newstride = slicemax ! increments the index of dimension newdim
            do id=uptodim,newdim-1
               newstride = newstride * inshape(id)
            enddo
              do ni=1,inshape(newdim)-1
                do ui=1,inshape(uptodim)-1
                  do is=1,slicemax
                    int_kq_roots(uptostride*ui + newstride*ni + is) = &
                         ( int_kq_roots(uptostride*ui + newstride*(ni-1) + is) - &
                          int_kq_roots(uptostride*(ui-1) + newstride*ni + is))/ &
                        (kroots(uptodim-1)-kroots(newdim-1))
                    enddo
                enddo
              enddo
        enddo
    end subroutine append_root
          

!  computes the integrals
! int_rational_ln_k_k0(q,p₁, ...,pₘ) = 
!            q
! ⌠ᵏᵇ       k         ┬─┬
! ⎮   ─────────────ln(│ │(k-ᵣk))dk
! ⌡   ┬─┬              ʳ
!     │ │(k-ₘk)ᵖ⁽ᵐ⁾
!      ᵐ
! for q=qₘᵢₙ...qₘₐₓ and p₁=0...max(p₁), ...,pₘ=0...max(pₘ)
! ᵣk = t_kroots(r)
! ₘk = kroots(m)
    subroutine drive_ln_k_k0_ints(int_rational_ln_k_k0,&
        int_kq_roots, kroots, kb, t_kroots, k_crossings,span_k_crossings, t3, spank1)
        implicit none
        type(mp_complex), dimension(q_minn:,0:,0:,0:) :: int_rational_ln_k_k0
        type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:) :: int_kq_roots
        type(mp_complex), dimension(2) :: t_kroots
        type(mp_complex), dimension(3) :: kroots
        type(mp_complex), dimension(size(kroots)+1) :: kroots_helper
        type(mp_complex) :: spence
        type(mp_real) :: kb
        integer ::corr1
        type(mp_real) :: angle_sum 
        logical :: t3, spank1
        external spence
        type(mp_real) :: mppic
        type(mp_real), dimension(:,:) :: k_crossings
        logical, dimension(:,:) :: span_k_crossings

        mppic=mppi(kv_nwds)

        ! possible split ln correction from ln((k-t_kroots(1))*(k-t_kroots(2))) = ln(k-kroots(1)) + ln(k-t_kroots(2)) +? 2*pi*i
        if(t3.and.spank1) then
          corr1=0 ! only one kv_root
        else
          angle_sum = arg(kb-t_kroots(1)) + arg(kb-t_kroots(2))
          call split_ln_corr(angle_sum, corr1)  ! corr1 might always be zero because the two t_kroots are nearly pure real or have oppositely signed imaginary parts
        endif

        kroots_helper(1) = mpcmplx(cmplx(0.0,0.0),kv_nwds)
        kroots_helper(2:2+size(kroots)-1) = kroots

        call do_case1_ln_k_k0_ints(int_rational_ln_k_k0, int_kq_roots, kroots_helper, kb, t_kroots, corr1, t3, spank1)

        ! dilogarithm/spence integrals
        call do_case2_ln_k_k0_ints(int_rational_ln_k_k0, kroots_helper, kb, t_kroots, corr1,&
            k_crossings, span_k_crossings, t3, spank1)

        call do_case3_ln_k_k0_ints(int_rational_ln_k_k0(0,0,0,0), kb, t_kroots, corr1, t3, spank1)

        call combine_multi_roots_ints(int_rational_ln_k_k0,kroots,1,spank1)

    end subroutine drive_ln_k_k0_ints

!  computes the integrals
! int_rational_ln_k_k0(...,pₘ,...) =  (all indices except pₘ are zero)
!             
! ⌠kb       1         ┬─┬
! ⎮   ─────────────ln(│ │(k-ᵣk))dk
! ⌡     (k-ₘk)ᵖ⁽ᵐ⁾     ʳ
! for each ₘk = kroots_helper(m) and in turn for a range of pₘ specified by the dimension bounds of int_rational_ln_k_k0 , except skips pₘ=1 and pₘ=0
! kroots_helper(1) = 0.0, kroots_helper(2:) = kroots
! ᵣk = t_kroots(r)
! corr1 : coefficient of a possible split ln correction 

subroutine do_case1_ln_k_k0_ints(int_rational_ln_k_k0, int_kq_roots, kroots_helper, kb, t_kroots, corr1, t3, spank1)
  implicit none

  type(mp_complex), dimension(q_minn:,0:,0:,0:) :: int_rational_ln_k_k0
  type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:) :: int_kq_roots
  type(mp_complex), dimension(:) :: t_kroots
  type(mp_complex), dimension(:) :: kroots_helper
  type(mp_real) :: kb
  integer :: corr1
  logical :: t3, spank1

  integer :: m, pn, pn2, ir
  integer, dimension(2) :: ids1
  integer, dimension(4) :: ids2, ids3
  type(mp_real) :: prefactor1
  type(mp_complex) :: prefactor2
  type(mp_real) :: mppic

  mppic=mppi(kv_nwds)
  do m=1,size(kroots_helper)
    if(spank1.and.(m.eq.(size(kroots_helper)))) then
        cycle 
    endif
    ids2(:) = 0
    ids3(:) = 0

    do pn=lbound(int_rational_ln_k_k0,m),ubound(int_rational_ln_k_k0,m)
      ids2(m) = pn
      if(m.eq.1) then
          pn2 = pn 
          ids3(m) = pn+1
          if(pn.eq.ubound(int_rational_ln_k_k0,m)) then
            cycle
          endif
      else
          pn2 = -pn
          ids3(m) = pn-1
          if(pn.eq.lbound(int_rational_ln_k_k0,m)) then
            cycle
          endif
      endif
      if((pn2.eq.-1).or.(pn2.eq.0)) then
        cycle
      endif

      prefactor1 = mpreal(1.0,kv_nwds)/(pn2+1.0)
      prefactor2 = (kb-kroots_helper(m))**(pn2+1)

      do ir=1,size(t_kroots)
        ids1(:) = 0
        ids1(ir) = 1
        if(t3.and.spank1.and.(ir.eq.size(t_kroots))) then
            cycle 
        endif

        ! from integration by parts
        int_rational_ln_k_k0(ids2(1),ids2(2),ids2(3),ids2(4)) = &
            int_rational_ln_k_k0(ids2(1),ids2(2),ids2(3),ids2(4)) + &
             prefactor1 * (  prefactor2*logw(kb-t_kroots(ir)) - &
            int_kq_roots(ids3(1),ids3(2),ids3(3),ids3(4),ids1(1),ids1(2))  )

      enddo
      ! split ln correction
      int_rational_ln_k_k0(ids2(1),ids2(2),ids2(3),ids2(4)) = &
          int_rational_ln_k_k0(ids2(1),ids2(2),ids2(3),ids2(4)) + &
          prefactor1 * prefactor2 * mpcmplx(2*corr1*i,kv_nwds)*mppic

    enddo
  enddo
end subroutine do_case1_ln_k_k0_ints

!  computes the integrals
! int_rational_ln_k_k0(...,pₘ=1,...) =  (all indices except pₘ are zero)
!
!  kb         
! ⌠     1       ┬─┬
! ⎮   ────── ln(│ │(k-ᵣk))dk
! ⌡   (k-ₘk)     ʳ
! for each ₘk = kroots_helper(m) using spence's/the dilogarithm function

! kroots_helper(1) = 0.0, kroots_helper(2:) = kroots
! ᵣk = t_kroots(r)
! corr1 : coefficient of a possible split ln correction 
! k_crossings(m,ir) : k values where a log function involved in the expression of this integral may have a branch cut, see www.github.com/brentfpage/fort-kv-ints-r/blob/main/preprint.pdf
! span_k_crossings(m,ir) : whether the considered interval k=klb->kub spans k_crossings(m,ir) 
subroutine do_case2_ln_k_k0_ints(int_rational_ln_k_k0, kroots_helper, kb, t_kroots, corr1,&
    k_crossings, span_k_crossings, t3, spank1)
  implicit none
  type(mp_complex), dimension(q_minn:,0:,0:,0:) :: int_rational_ln_k_k0
  type(mp_complex), dimension(:) :: t_kroots
  type(mp_complex), dimension(:) :: kroots_helper
  type(mp_real) :: kb
  integer :: corr1
  logical :: t3, spank1
  type(mp_real), dimension(:,:) :: k_crossings
  logical, dimension(:,:) :: span_k_crossings

  type(mp_complex) :: spence_arg
  type(mp_complex) :: spence_res
  type(mp_complex) :: spence_arg_at_k_crossing
  type(mp_complex) :: spence
  external spence

  type(mp_real) :: angle_sum 
  integer :: corr2

  integer :: m, ir
  integer, dimension(4) :: ids2
  type(mp_real) :: mppic


  mppic = mppi(kv_nwds)
  do m=1,size(kroots_helper)
    if(spank1.and.(m.eq.size(kroots_helper))) then
        cycle 
    endif
    ids2(:) = 0
    if(m.eq.1) then
        ids2(m) = -1
    else
        ids2(m) = 1
    endif

    do ir=1,size(t_kroots)
      if(t3.and.spank1.and.(ir.eq.size(t_kroots))) then
          cycle 
      endif

      spence_arg = (t_kroots(ir) - kb)/(t_kroots(ir)-kroots_helper(m))
      spence_res = spence(spence_arg)

      int_rational_ln_k_k0(ids2(1),ids2(2),ids2(3),ids2(4)) = &
        int_rational_ln_k_k0(ids2(1),ids2(2),ids2(3),ids2(4)) &
        - spence_res + logw(kb-kroots_helper(m))*logw(kroots_helper(m)-t_kroots(ir))

      angle_sum = arg(spence_arg) + arg(kroots_helper(m)-t_kroots(ir))
      call split_ln_corr(angle_sum, corr2)
      
      if(corr2.ne.0) then
        spence_arg_at_k_crossing = (t_kroots(ir)-k_crossings(m,ir))/(t_kroots(ir)-kroots_helper(m))
        if (span_k_crossings(m,ir).and.(mpreal(spence_arg_at_k_crossing,kv_nwds).lt.mpreal(0.0,kv_nwds))) then
          int_rational_ln_k_k0(ids2(1),ids2(2),ids2(3),ids2(4)) = &
            int_rational_ln_k_k0(ids2(1),ids2(2),ids2(3),ids2(4)) &
            + mpcmplx(2*i*corr2,kv_nwds)* mppic*&
            (logw(kb-kroots_helper(m)) - logw(k_crossings(m,ir)-kroots_helper(m)))

          int_rational_ln_k_k0(ids2(1),ids2(2),ids2(3),ids2(4)) = &
            int_rational_ln_k_k0(ids2(1),ids2(2),ids2(3),ids2(4)) &
            + mpcmplx(2*i*corr2,kv_nwds)* mppic* &
            logw(abs(mpreal(1.0,kv_nwds)-spence_arg_at_k_crossing))
        else
          int_rational_ln_k_k0(ids2(1),ids2(2),ids2(3),ids2(4)) = &
            int_rational_ln_k_k0(ids2(1),ids2(2),ids2(3),ids2(4)) &
              + mpcmplx(2*i*corr2,kv_nwds)* mppic*logw(kb-kroots_helper(m))
        endif
      endif
    enddo

    int_rational_ln_k_k0(ids2(1),ids2(2),ids2(3),ids2(4)) = &
      int_rational_ln_k_k0(ids2(1),ids2(2),ids2(3),ids2(4)) &
      + logw(kb-kroots_helper(m))*mpcmplx(2*i*corr1,kv_nwds)*mppic

  enddo
end subroutine do_case2_ln_k_k0_ints


!  computes
! int_rational_ln_k_k0_0000 = 
!     
! ⌠kb    ┬─┬
! ⎮   ln(│ │(k-ᵣk))dk
! ⌡       ʳ
!     
! ᵣk = t_kroots(r)
subroutine do_case3_ln_k_k0_ints(int_rational_ln_k_k0_0000, kb, t_kroots, corr1, t3, spank1)
  implicit none
  type(mp_complex) :: int_rational_ln_k_k0_0000
  type(mp_real) :: kb
  type(mp_complex), dimension(:) :: t_kroots
  integer :: corr1
  logical :: t3, spank1

  integer :: ir
  type(mp_real) :: mppic

  mppic=mppi(kv_nwds)
  do ir=1,size(t_kroots)
    if(t3.and.spank1.and.(ir.eq.size(t_kroots))) then
        cycle 
    endif
    int_rational_ln_k_k0_0000 = int_rational_ln_k_k0_0000 + &
        logw(kb-t_kroots(ir))*(kb-t_kroots(ir)) - (kb-t_kroots(ir))
  enddo

  ! possible split ln correction
  int_rational_ln_k_k0_0000 = int_rational_ln_k_k0_0000 + mpcmplx(2*i*corr1,kv_nwds)*mppic*kb
end subroutine do_case3_ln_k_k0_ints


! drives computation of int_kq_roots_k(l/u)b(q,p₁,p₂,p₃,0,0) = 
!  k(l/u)b     q
! ⌠           k
! ⎮ dk₂ ─────────────
! ⌡      ₃
!       ┬─┬
!       │ │(k-ₘk)ᵖ⁽ᵐ⁾
!       m=1
! where ₘk=kroots(m).  Minor exceptions to this description apply in the spank1 case.
! This function also computes ₂k and ₃k , which are roots of a quadratic polynomial h(k₂) = k₂*ω₁ - k₁*ω₂ + k₁- k₂ that is introduced by application of partial fraction decomposition operations w.r.t v_parallel to the integral described in the header of int_driver.
! In the spank1 case, k1=₁k is one of the roots of this polynomial, so the integral above is only computed for two total roots.
subroutine k_int_driver(kub,klb,kroots,int_kq_roots_kub,int_kq_roots_klb,&
    k1,om1,om2splcoeffs,pfdsplcoeffs,spank1)
    implicit none
    type(mp_complex), dimension(:), intent(inout) :: kroots
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(out) :: int_kq_roots_kub, int_kq_roots_klb
    type(mp_real), intent(in) :: k1, om1, kub, klb
    type(mp_real), dimension(:), intent(in) :: om2splcoeffs
    type(mp_real), dimension(:), intent(out) :: pfdsplcoeffs
    logical :: spank1

    integer :: ip
    type(mp_real) :: sqrt_arg
    type(mp_complex), dimension(size(kroots)+1) :: kroots_helper

! pfdsplcoeffs contains the spline coefficients of h(k₂) = k₂*ω₁ - k₁*ω₂ + k₁ - k₂, where
! ω₂ = om2splcoeffs(1) + om2splcoeffs(2) * k₂ + om2splcoeffs(3) * pow(k₂, 2) .
! If the k integral bounds span k1, i.e., klb < k1 < kub, then k1 is a root of h(k₂) .
! In this case, (k₂ - k₁) gets factored out of h(k₂) .
        if (spank1) then
          pfdsplcoeffs(1) = om1 - k1*om2splcoeffs(2)-k1**2*om2splcoeffs(3) - 1.0 ! k^0 term
          pfdsplcoeffs(2) = -k1*om2splcoeffs(3) ! k^1 term
          pfdsplcoeffs(3) = mpreal(0.0,kv_nwds)

          kroots(2) = -pfdsplcoeffs(1)/pfdsplcoeffs(2)
          kroots(3) = mpreal(0.0,kv_nwds) ! second pfd root is k1

        else
          do ip=1,size(pfdsplcoeffs)
            pfdsplcoeffs(ip) = -k1*om2splcoeffs(ip)
          enddo
          pfdsplcoeffs(2) = pfdsplcoeffs(2) + om1 - 1.0
          pfdsplcoeffs(1) = pfdsplcoeffs(1) + k1

          sqrt_arg = pfdsplcoeffs(2)**2-4.0*pfdsplcoeffs(3)*pfdsplcoeffs(1)
          if(sqrt_arg.ge.0.0) then
            kroots(2) =  (-pfdsplcoeffs(2)+sqrt(sqrt_arg)) / (2.0 * pfdsplcoeffs(3)) 
            kroots(3) =  (-pfdsplcoeffs(2)-sqrt(sqrt_arg)) / (2.0 * pfdsplcoeffs(3))
          else
            sqrt_arg = -sqrt_arg
            kroots(2) =  (-pfdsplcoeffs(2)+sqrt(sqrt_arg)*mpcmplx(i,kv_nwds)) / (2.0 * pfdsplcoeffs(3)) 
            kroots(3) =  (-pfdsplcoeffs(2)-sqrt(sqrt_arg)*mpcmplx(i,kv_nwds)) / (2.0 * pfdsplcoeffs(3))
          endif

        endif
        int_kq_roots_klb = mpcmplx((0.0,0.0),kv_nwds)
        int_kq_roots_kub = mpcmplx((0.0,0.0),kv_nwds)

        kroots_helper(1) = mpcmplx(cmplx(0.0,0.0),kv_nwds)
        kroots_helper(2:2+size(kroots)-1) = kroots

        call do_rational_ints(int_kq_roots_kub, kub, kroots_helper, 1, spank1)
        call do_rational_ints(int_kq_roots_klb, klb, kroots_helper, 1, spank1)

end subroutine k_int_driver


! with
!  v(l/u)b
! ⌠       vⁿ
! ⎮ dv ────────
! ⌡          λ₁
!      (v-vₒ)
! having been pre-computed for vₒ=(-1+ω₁+i*eps)/k₁ in v_int_v(l/u)b_lam1_ipara
! and
!  k(l/u)b        q      
! ⌠              k₂       
! ⎮    dk₂ ───────────── 
! ⌡         ₃            
!          ┬─┬           
!          │ │(k₂-ₘk)ᵖ⁽ᵐ⁾
!          ᵐ⁼¹
! have been pre-computed in int_kq_roots_k(l/u)b (q,p1,p2,p3), assemble_t1_ints computes the product, i.e.,
!                            q
! ⌠        vⁿ   ⌠           k₂
! ⎮ dv ─────────⎮ dk₂ ──────────────
! ⌡           λ₁⌡      ₃
!     (-t₁/k₁)        ┬─┬
!                     │ │(k₂-ₘk)ᵖ⁽ᵐ⁾
!                     ᵐ⁼¹
! where t₁ = ω₁ - v*k₁ - 1 + i*eps
subroutine assemble_t1_ints(int_kq_roots_klb,int_kq_roots_kub,v_int_vlb_lam1_ipara,v_int_vub_lam1_ipara,t_int,spank1)
  implicit none
  type(mp_complex), dimension(0:,0:), intent(in) :: v_int_vub_lam1_ipara
  type(mp_complex), dimension(0:,0:), intent(in) :: v_int_vlb_lam1_ipara
  type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:) :: int_kq_roots_klb
  type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:) :: int_kq_roots_kub
  type(mp_complex) :: k_int
  type(mp_complex), dimension(q_minn:,0:,0:,0:,0:) :: t_int
  integer :: qq, tau, sig, n, lam1
  logical :: spank1

  do qq=lbound(int_kq_roots_klb,1),ubound(int_kq_roots_klb,1)
    do tau=0,ubound(int_kq_roots_klb,2)!-ubound(int_kq_roots_klb,3) 
      do sig=0,ubound(int_kq_roots_klb,3)
        if(spank1) then ! one of the pfd kroots is k1
          k_int = (int_kq_roots_kub(qq,tau,sig,0,0,0) - int_kq_roots_klb(qq,tau,sig,0,0,0))
        else
          k_int = (int_kq_roots_kub(qq,tau,sig,sig,0,0) - int_kq_roots_klb(qq,tau,sig,sig,0,0))
        endif

        do n=0,ubound(t_int,4)
          do lam1=0,ubound(t_int,5)
            t_int(qq,tau,sig,n,lam1) = k_int *&
              & (v_int_vub_lam1_ipara(n,lam1)-v_int_vlb_lam1_ipara(n,lam1))
          enddo
        enddo
      enddo
    enddo
  enddo
end subroutine assemble_t1_ints



! performs the following integrals:
! t_int(q,τ,σ,n,λ) = 
!   if t3
!                 λ₃
! ⌠              k₂           ⌠           vⁿ
! ⎮ dk₂ ───────────────────── ⎮ dv ────────────────
! ⌡           τ      σ      σ ⌡                  λ₃
!       (k-₁k) (k-₂k) (k-₃k)       ⎛    ω₁-ω₂-iϵ⎞
!                                  ⎜v - ────────⎟
!                                  ⎝     k₁-k₂  ⎠
!   else
!                 λ₂
! ⌠              k₂           ⌠           vⁿ
! ⎮ dk₂ ───────────────────── ⎮ dv ────────────────
! ⌡           τ      σ      σ ⌡                  λ₂
!       (k-₁k) (k-₂k) (k-₃k)       ⎛    -ω₂+1+iϵ⎞
!                                  ⎜v - ────────⎟
!                                  ⎝       k₂   ⎠

! where ₘk = kroots(m).
! Here, ω₂ = om2splcoeffs(1) + om2splcoeffs(2) * k₂ + om2splcoeffs(3) * pow(k₂, 2) .

! In addition,
!  k(l/u)b        q      
! ⌠              k₂       
! ⎮    dk₂ ───────────── 
! ⌡         ₃            
!          ┬─┬           
!          │ │(k₂-ₘk)ᵖ⁽ᵐ⁾
!          ᵐ⁼¹
! have been pre-computed in int_kq_roots_k(l/u)b (q,p1,p2,p3)
    subroutine t23_int_driver(kub,klb,kroots,t_int,int_kq_roots_kub,int_kq_roots_klb,vlb,vub,k1,om1,om2splcoeffs, t3, spank1)
        implicit none
        type(mp_complex), dimension(:) :: kroots
        type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:) :: int_kq_roots_kub, int_kq_roots_klb
        type(mp_complex), dimension(q_minn:,0:,0:,0:,0:) :: t_int
        type(mp_real) :: k1, om1, kub, klb, vlb, vub
        type(mp_real), dimension(:) :: om2splcoeffs
        logical :: t3, spank1

        t_int = mpcmplx(cmplx(0.0,0.0),kv_nwds) ! necessary to zero out this slice of t123_int to avoid some double counting of pure n, q integrals

        call indef_t23_int_driver(klb,kub,t_int,int_kq_roots_kub,int_kq_roots_klb,kroots,vub,k1,om1,om2splcoeffs,1,t3,spank1)
        call indef_t23_int_driver(klb,kub,t_int,int_kq_roots_kub,int_kq_roots_klb,kroots,vlb,k1,om1,om2splcoeffs,-1,t3,spank1)

    end subroutine t23_int_driver

! does the integral described above "subroutine t23_int_driver" for one velocity bound.
    subroutine indef_t23_int_driver(klb,kub,t_int,int_kq_roots_kub,int_kq_roots_klb,&
      kroots_in,vb,k1,om1,om2splcoeffs,sgn,t3,spank1)
      implicit none
      type(mp_real), dimension(:), intent(in) :: om2splcoeffs
      type(mp_complex), dimension(:), intent(in) :: kroots_in
      type(mp_complex), dimension(size(kroots_in)) :: kroots
      type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(inout) :: int_kq_roots_kub
      type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(inout) :: int_kq_roots_klb
      type(mp_complex), dimension(q_minn:,0:,0:,0:,0:), intent(inout) :: t_int
      type(mp_complex), dimension(q_minn:ubound(t_int,1),0:ubound(t_int,2),0:ubound(t_int,3),&
          & 0:ubound(t_int,4),0:ubound(t_int,5)) :: t_int_at_vb ! dims: q, k1_root_pow, other_roots_pow, n, lam
      type(mp_real) :: vb, k1, klb, kub, om1
      integer :: n, lam, q,sig, sgn
      logical :: t3, spank1
      integer :: tau
      real :: fact
      external fact
      type(mp_real), dimension(size(kroots_in)+1,2) :: k_crossings
      logical, dimension(size(kroots_in)+1,2) :: span_k_crossings
      type(mp_complex), dimension(size(om2splcoeffs)) :: tsplcoeffs
      type(mp_complex), dimension(size(om2splcoeffs)-1) :: t_kroots

      t_int_at_vb = mpcmplx(cmplx(0.0,0.0),kv_nwds)
      
      kroots(:) = kroots_in(:)

      call compute_tsplcoeffs_and_roots(om2splcoeffs,tsplcoeffs,t_kroots,k1,om1,vb,spank1,t3)

      call find_cut_crossings(kroots, t_kroots, k_crossings, span_k_crossings, klb, kub)

      call inindef_t23_int_driver(kub,t_int_at_vb,int_kq_roots_kub,kroots,vb,k1,om1,&
          om2splcoeffs,tsplcoeffs,t_kroots,k_crossings,span_k_crossings,1,t3,spank1)
      call inindef_t23_int_driver(klb,t_int_at_vb,int_kq_roots_klb,kroots,vb,k1,om1,&
          om2splcoeffs,tsplcoeffs,t_kroots,k_crossings,span_k_crossings,-1,t3,spank1)

      do q=lbound(t_int_at_vb,1),ubound(t_int_at_vb,1)
        do tau=0,ubound(t_int_at_vb,2)
          do sig=0,ubound(t_int_at_vb,3)
            do n=0,ubound(t_int_at_vb,4)
              do lam=0,ubound(t_int_at_vb,5)
                t_int(q,tau,sig,n,lam) = t_int(q,tau,sig,n,lam) + (sgn*1.0)*t_int_at_vb(q,tau,sig,n,lam)
              enddo
            enddo
          enddo
        enddo
      enddo

    end subroutine indef_t23_int_driver

  ! tsplcoeffs: the k₂ spline coefficients of either 
  ! t₂(vb) = –ω₂ + k₂(vb) + 1 + i*eps or 
  ! t₃(vb) = ω₁–ω₂ - (vb)(k₁-k₂) + i*eps, depending on whether t3=.true. or .false.
  ! vb : a parallel velocity bound for the considered integral
  ! t_kroots: the zeros of t₂(vb) or t₃(vb). i.e., t(vb) = 0 when evaluated at a t_kroot
subroutine compute_tsplcoeffs_and_roots(om2splcoeffs, tsplcoeffs, t_kroots, k1,om1,vb , spank1, t3)
  implicit none
  type(mp_real), dimension(:), intent(in) :: om2splcoeffs
  type(mp_complex), dimension(:), intent(out) :: tsplcoeffs, t_kroots
  type(mp_real) :: k1, om1, vb
  logical :: spank1, t3

  if(t3.and.spank1) then
! If the k integral bounds span k1, i.e., klb < k1 < kub, then k1 is a root of t₃ .
! In this case, (k₂ - k₁) gets factored out of t₃ .
    tsplcoeffs(3) = -om2splcoeffs(2)-k1*om2splcoeffs(3)+vb + mpcmplx(i,kv_nwds) * eps
    tsplcoeffs(2) = -om2splcoeffs(3)
    tsplcoeffs(1) = tsplcoeffs(2) ! put the coefficient of the highest power in
    ! the lowest index to make later code work
    t_kroots(1) = -tsplcoeffs(3)/tsplcoeffs(2)
    t_kroots(2) = mpcmplxdc(cmplx(0.0,0.0),kv_nwds) ! no second root
  else
    if(t3) then
      tsplcoeffs(3) = -(-om1 + om2splcoeffs(1) + vb*k1 - mpcmplx(i,kv_nwds) * eps)
      tsplcoeffs(2) = -(om2splcoeffs(2) - vb)
      tsplcoeffs(1) = -om2splcoeffs(3)
    else
      tsplcoeffs(3) = -om2splcoeffs(1) + 1.0 + mpcmplx(i,kv_nwds) * eps
      tsplcoeffs(2) = -om2splcoeffs(2) + vb
      tsplcoeffs(1) = -om2splcoeffs(3)
    endif
      t_kroots(1) =  (-tsplcoeffs(2) + &
  & sqrt(tsplcoeffs(2)**2-4.0*tsplcoeffs(1)*tsplcoeffs(3)))/(2.0*tsplcoeffs(1))
      t_kroots(2) =  (-tsplcoeffs(2) - &
  & sqrt(tsplcoeffs(2)**2-4.0*tsplcoeffs(1)*tsplcoeffs(3)))/(2.0*tsplcoeffs(1))
  endif
end subroutine compute_tsplcoeffs_and_roots


! determine the k_crossings where the argument of a log function involved in the spence integral may cross a branch cut, see www.github.com/brentfpage/fort-kv-ints-r/blob/main/preprint.pdf
subroutine find_cut_crossings(kroots, t_kroots, k_crossings, span_k_crossings, klb, kub)
  implicit none
  type(mp_complex), dimension(:), intent(in) :: kroots, t_kroots
  type(mp_real), dimension(:,:), intent(out) :: k_crossings
  logical, dimension(:,:), intent(out) :: span_k_crossings
  type(mp_real) :: klb, kub

  integer :: ir, j
  type(mp_complex), dimension(size(kroots)+1) :: kroots_helper

  kroots_helper(1) = mpcmplx(cmplx(0.0,0.0),kv_nwds)
  kroots_helper(2:2+size(kroots)-1) = kroots

  do j=1,size(kroots_helper)
    do ir=1,2
      k_crossings(j,ir) = aimag(conjg(t_kroots(ir))*(kroots_helper(j)))/aimag( conjg(t_kroots(ir) - kroots_helper(j)) )
      if ((klb.lt.k_crossings(j,ir)).and.(kub.gt.k_crossings(j,ir))) then
        span_k_crossings(j,ir) = .true.
      else
        span_k_crossings(j,ir) = .false.
      endif
    enddo
  enddo
end subroutine find_cut_crossings

! does the integral described above "subroutine t23_int_driver" for one velocity bound and one k bound
! tsplcoeffs : contains k spline coefficients of either 
!     t₂ = –ω₂ + k₂(vb) + 1 + i*eps or 
!     t₃ = ω₁–ω₂ - (vb)(k₁-k₂) + i*eps, depending on whether t3=.true. or .false.
! vb : a velocity bound for the considered integral
! t_kroots : the zeros of the quadratic polynomial specified by tsplcoeffs
    subroutine inindef_t23_int_driver(kb,t_int_at_vb,int_kq_roots,kroots_in,vb,k1,om1,&
        om2splcoeffs,tsplcoeffs,t_kroots,k_crossings,span_k_crossings,sgn,t3,spank1)
        implicit none
        type(mp_real), dimension(:), intent(in) :: om2splcoeffs
        type(mp_complex), dimension(:), intent(in) :: kroots_in
        type(mp_complex), dimension(size(kroots_in)) :: kroots
        type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(inout) :: int_kq_roots
        type(mp_complex), dimension(&
            & q_minn:ubound(int_kq_roots,1), &
            & 0:ubound(int_kq_roots,2), &
            & 0:ubound(int_kq_roots,3), &
            & 0:ubound(int_kq_roots,4)) :: int_rational_ln_k_k0
        type(mp_complex), dimension(q_minn:,0:,0:,0:,0:), intent(inout) :: t_int_at_vb
        type(mp_complex), dimension(q_minn:ubound(t_int_at_vb,1),0:ubound(t_int_at_vb,2),&
          0:ubound(t_int_at_vb,3), 0:ubound(t_int_at_vb,4),0:ubound(t_int_at_vb,5)) :: t_int_at_vb_kb ! dims: q, k1_root_pow, other_roots_pow, n, lam
        type(mp_real) :: vb, k1, kb, om1
        type(mp_complex), dimension(:) :: tsplcoeffs
        type(mp_complex), dimension(:) :: t_kroots
        type(mp_complex), dimension(1 + size(kroots) + size(t_kroots)) :: kroots_full
        integer :: n, lam, q, sig, tau, sgn
        logical :: t3, spank1
        type(mp_real), dimension(:,:) :: k_crossings
        logical, dimension(:,:) :: span_k_crossings

        kroots(:) = kroots_in
        kroots_full(1) = mpcmplx((0.0,0.0),kv_nwds)
        kroots_full(2:size(kroots)+1) = kroots
        kroots_full(size(kroots)+2:) = t_kroots


! with int_kq_roots(q,p₁,p₂,p₃,y₁,y₂)
!                   q
! ⌠                k
! ⎮ dk ──────────────────────────
! ⌡     ₃            ₂
!      ┬─┬          ┬─┬
!      │ │(k-ₘk)ᵖ⁽ᵐ⁾│ │(k-ₙk)ʸ⁽ⁿ⁾
!      ᵐ⁼¹          ⁿ⁼¹

! having already been computed for y₁=0,y₂=0, doing the y₁>0,y₂>0 cases below.
! here, ₘk=kroots(m) , ₙk=t_kroots(n) , and kroots_full is the concatenation of 0.0, kroots, and t_kroots. 
! small exceptions to this description apply if spank1
  if(spank1) then
    call do_rational_ints(int_kq_roots(:,:,:,0,:,:),&
          kb,&
          kroots_full((/ 1, 2, 3, 5, 6 /)),&
          size(kroots)+1,&
          t3)
    else
      call do_rational_ints(int_kq_roots,&
          kb,&
          kroots_full,&
          size(kroots)+2,&
          .false.)
    endif


!   t_int_at_vb_kb is the integral described above t23_int_driver for a single v_parallel bound and single k bound.
  t_int_at_vb_kb(:,:,:,:,:) = mpcmplx(cmplx(0.0,0.0),kv_nwds)
  call int_kq_roots_2_t_int_lam0(t_int_at_vb_kb, int_kq_roots, vb, spank1)
  call int_kq_roots_2_t_int_lam_gt1(t_int_at_vb_kb, int_kq_roots, k1, tsplcoeffs, spank1, t3)

    ! do lam=1 for all q, sigma, n=0
    int_rational_ln_k_k0(:,:,:,:) = mpcmplx(cmplx(0.0,0.0),kv_nwds)
    call drive_ln_k_k0_ints(int_rational_ln_k_k0,int_kq_roots, kroots, kb, t_kroots,&
      k_crossings,span_k_crossings,t3, spank1)

    if(spank1) then
        do sig=0,ubound(t_int_at_vb_kb,3)
          t_int_at_vb_kb(:,:,sig,0,1) = int_rational_ln_k_k0(:,:,sig,0)
        enddo
    else
      do sig=0,ubound(t_int_at_vb_kb,3)
        t_int_at_vb_kb(:,:,sig,0,1) = int_rational_ln_k_k0(:,:,sig,sig)
      enddo
    endif

    call increment_n(t_int_at_vb_kb, k1, om1, om2splcoeffs, spank1, t3)

    do q=lbound(t_int_at_vb,1),ubound(t_int_at_vb,1)
      do tau=0,ubound(t_int_at_vb,2)
        do sig=0,ubound(t_int_at_vb,3)
          do n=0,ubound(t_int_at_vb,4)
            do lam=0,ubound(t_int_at_vb,5)
              t_int_at_vb(q,tau,sig,n,lam) = t_int_at_vb(q,tau,sig,n,lam) + (sgn*1.0)*t_int_at_vb_kb(q,tau,sig,n,lam)
            enddo
          enddo
        enddo
      enddo
    enddo
end subroutine inindef_t23_int_driver


! converts values of int_kq_roots to values of t_int_at_vb_kb(...,lam) for lam=0
subroutine int_kq_roots_2_t_int_lam0(t_int_at_vb_kb, int_kq_roots, vb, spank1)
  implicit none

  type(mp_complex), dimension(q_minn:,0:,0:, 0:,0:), intent(inout) :: t_int_at_vb_kb ! dims: q, k1_root_pow, other_roots_pow, n, lam
  type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(in) :: int_kq_roots
  type(mp_real) :: vb
  logical :: spank1

  integer :: q, sig, ir, tau, n

  if(spank1) then
    do sig=0,ubound(int_kq_roots,3)
      t_int_at_vb_kb(:,:,sig,:,0) = spread(int_kq_roots(:,:,sig,0,0,0), 3, size(t_int_at_vb_kb,4)) !one of the pdf kroots is k1
    enddo
  else 
    do ir=0,ubound(int_kq_roots,3)
      t_int_at_vb_kb(:,:,ir,:,0) = spread(int_kq_roots(:,:,ir,ir,0,0), 3, size(t_int_at_vb_kb,4))
    enddo
  endif

  do q=lbound(t_int_at_vb_kb,1),ubound(t_int_at_vb_kb,1)
    do tau=0,ubound(t_int_at_vb_kb,2)
      do sig=0,ubound(t_int_at_vb_kb,3)
        do n=0,ubound(t_int_at_vb_kb,4)
          t_int_at_vb_kb(q,tau,sig,n,0) = t_int_at_vb_kb(q,tau,sig,n,0) * vb**(n+1)/mpreal(n+1.0,kv_nwds)
        enddo
      enddo
    enddo
  enddo
end subroutine int_kq_roots_2_t_int_lam0

! converts values of int_kq_roots to values of t_int_at_vb_kb(...,lam) for lam>1 
subroutine int_kq_roots_2_t_int_lam_gt1(t_int_at_vb_kb, int_kq_roots, k1, tsplcoeffs, spank1, t3)
  implicit none
  type(mp_complex), dimension(q_minn:,0:,0:,0:,0:), intent(inout) :: t_int_at_vb_kb ! dims: q, k1_root_pow, other_roots_pow, n, lam
  type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(in) :: int_kq_roots
  type(mp_real) :: k1
  type(mp_complex), dimension(:) :: tsplcoeffs
  logical :: spank1, t3

  integer :: q, sig, tau, lam, lam_prime, lamlam
  real :: fact
  external :: fact
  if(spank1) then
    if (t3) then
      do q=q_minn,ubound(t_int_at_vb_kb,1) 
        do lam=2,ubound(t_int_at_vb_kb,5)
          do tau=0,ubound(t_int_at_vb_kb,2)
            do sig=0,ubound(t_int_at_vb_kb,3)
              t_int_at_vb_kb(q,tau,sig,0,lam) = int_kq_roots(q,tau,sig,0,lam-1,0)/&
                mpreal(-lam+1.0,kv_nwds)/tsplcoeffs(1)**(lam-1)
            enddo
          enddo
        enddo
      enddo
    else
      do lam=2,ubound(t_int_at_vb_kb,5)
        do q=q_minn,ubound(t_int_at_vb_kb,1) - (lam-1)
          do tau=0,ubound(t_int_at_vb_kb,2)
            do sig=0,ubound(t_int_at_vb_kb,3)
              t_int_at_vb_kb(q,tau,sig,0,lam) = int_kq_roots(q+lam-1,tau,sig,0,lam-1,lam-1)/&
                &mpreal(-lam+1.0,kv_nwds)/tsplcoeffs(1)**(lam-1)
            enddo
          enddo
        enddo
      enddo
    endif
  else
    if (t3) then 
      do lam=2,ubound(t_int_at_vb_kb,5)
        do q=q_minn,ubound(t_int_at_vb_kb,1)-(lam-1)
          do tau=0,ubound(t_int_at_vb_kb,2) 
            lam_prime = lam - tau - 1
            do sig=0,ubound(int_kq_roots,3)
              if (lam_prime.gt.0) then
                do lamlam=0,lam_prime
                  t_int_at_vb_kb(q,tau,sig,0,lam) = t_int_at_vb_kb(q,tau,sig,0,lam) + &
      & ( fact(lam_prime)/mpreald(fact(lamlam)*fact(lam_prime-lamlam),kv_nwds) * &
      & (-k1)**(lam_prime-lamlam) / (-lam + 1) / tsplcoeffs(1)**(lam-1) ) &
      & * int_kq_roots(q+lamlam,0,sig,sig,lam-1,lam-1)
                enddo
              else
                t_int_at_vb_kb(q,tau,sig,0,lam) = (mpreal(1.0,kv_nwds)/ mpreal(-lam+1.0,kv_nwds) / tsplcoeffs(1)**(lam-1) )&
                  &*int_kq_roots(q,-lam_prime,sig,sig,lam-1,lam-1)
              endif
            enddo
          enddo
        enddo
      enddo
    else
      do lam=2,ubound(t_int_at_vb_kb,5)
        do q=q_minn,ubound(t_int_at_vb_kb,1) - 1 * (lam-1) 
          do tau=0,ubound(int_kq_roots,2)
            do sig=0,ubound(int_kq_roots,3)
              t_int_at_vb_kb(q,tau,sig,0,lam) = int_kq_roots(q+lam-1,tau,sig,sig,lam-1,lam-1)&
                  /( (-lam+1.0)*tsplcoeffs(1)**(lam-1) )
            enddo
          enddo
        enddo
      enddo
    endif
  endif
end subroutine int_kq_roots_2_t_int_lam_gt1


! assemble computations of t_int_at_vb_kb for n>0 from those already done for n=0 by applying
!    v         a
!   ─── = 1 + ───
!   v-a       v-a
subroutine increment_n(t_int_at_vb_kb, k1, om1, om2splcoeffs, spank1, t3)
  implicit none
  type(mp_complex), dimension(q_minn:,0:,0:,0:,0:), intent(inout) :: t_int_at_vb_kb ! dims: q, k1_root_pow, other_roots_pow, n, lam
  type(mp_real) :: k1, om1
  type(mp_real), dimension(:), intent(in) :: om2splcoeffs
  logical :: spank1, t3

  type(mp_complex), dimension(3) :: helper_splcoeffs
  integer :: n, tau, sig, lam, q, iq, q2

  if (t3) then
    if(spank1) then

      helper_splcoeffs(1) = om2splcoeffs(2)+k1*om2splcoeffs(3)
      helper_splcoeffs(2) = om2splcoeffs(3)
      helper_splcoeffs(3) = mpreal(0.0,kv_nwds)

      do n=1,ubound(t_int_at_vb_kb,4)
        do tau=0,ubound(t_int_at_vb_kb,2)
          do sig=0,ubound(t_int_at_vb_kb,3)
            do lam=1,ubound(t_int_at_vb_kb,5)
              do q=lbound(t_int_at_vb_kb,1),ubound(t_int_at_vb_kb,1)-2*n
                t_int_at_vb_kb(q,tau,sig,n,lam) = t_int_at_vb_kb(q,tau,sig,n-1,lam-1)
                iq = 1
                do q2=q,q+size(helper_splcoeffs)-2 
                  t_int_at_vb_kb(q,tau,sig,n,lam) = t_int_at_vb_kb(q,tau,sig,n,lam) + &
                      t_int_at_vb_kb(q2,tau,sig,n-1,lam)*helper_splcoeffs(iq)
                  iq = iq + 1
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    else
      helper_splcoeffs(1) = om1-om2splcoeffs(1) + mpcmplx(i,kv_nwds)*eps 
      helper_splcoeffs(2) = -om2splcoeffs(2) 
      helper_splcoeffs(3) = -om2splcoeffs(3) 


      do n=1,ubound(t_int_at_vb_kb,4)
        do tau=0,ubound(t_int_at_vb_kb,2)-n
          do lam=1,ubound(t_int_at_vb_kb,5)
            do q=lbound(t_int_at_vb_kb,1),ubound(t_int_at_vb_kb,1)-2*n
              do sig=0,ubound(t_int_at_vb_kb,3)
                t_int_at_vb_kb(q,tau,sig,n,lam) = t_int_at_vb_kb(q,tau,sig,n-1,lam-1)
                iq = 1
                do q2=q,q+size(helper_splcoeffs)-1 !q -> q+1
                  t_int_at_vb_kb(q,tau,sig,n,lam) = t_int_at_vb_kb(q,tau,sig,n,lam) - &
                      t_int_at_vb_kb(q2,tau+1,sig,n-1,lam)*helper_splcoeffs(iq)
                  iq = iq+1
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    endif
  else
    helper_splcoeffs(1) = om2splcoeffs(1) - 1.0 - mpcmplx(i,kv_nwds)*eps 
    helper_splcoeffs(2) = om2splcoeffs(2) 
    helper_splcoeffs(3) = om2splcoeffs(3) 
    do n=1,ubound(t_int_at_vb_kb,4)
      do tau=0,ubound(t_int_at_vb_kb,2)
        do sig=0,ubound(t_int_at_vb_kb,3)
          do lam=1,ubound(t_int_at_vb_kb,5)
            do q=lbound(t_int_at_vb_kb,1) + n,ubound(t_int_at_vb_kb,1)-2*n
              t_int_at_vb_kb(q,tau,sig,n,lam) = t_int_at_vb_kb(q,tau,sig,n-1,lam-1) 
              iq = 1
              do q2=q-1,q+size(helper_splcoeffs)-2 
                t_int_at_vb_kb(q,tau,sig,n,lam) = t_int_at_vb_kb(q,tau,sig,n,lam) + &
                    t_int_at_vb_kb(q2,tau,sig,n-1,lam)*helper_splcoeffs(iq)
                iq=iq+1
              enddo 
            enddo
          enddo
        enddo
      enddo
    enddo
  endif
end subroutine increment_n


! with 
! ⌠ᵏᵇ  q
! ⎮   k  f(k) dk
! ⌡
! and
! ⌠ᵏᵇ    1
! ⎮   ──────── f(k) dk
! ⌡   (k-ₘk)ᵖᵐ
! having already been computed for qₘᵢₙ...qₘₐₓ and pₘ=0...max(pₘ), combineN_single_root_ints computes
! int_kq_roots(q,0,...,0,pₘ,0,...,0) =
!         q
! ⌠ᵏᵇ    k
! ⎮   ──────── f(k) dk
! ⌡   (k-ₘk)ᵖᵐ
! where m=opdimnum, and ₘk = kroot.  
!                                                               ⎛┬─┬      ⎞
! Also, f(k) is any function. In practice, f(k) = 1 or f(k) = ln⎜│ │(k-kᵣ)⎟ for this program
!                                                               ⎝ ʳ       ⎠
subroutine combine6_single_root_ints(int_kq_roots,kroot,opdimnum)
    implicit none
    type(mp_complex) :: kroot
    integer :: opdimnum
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(inout) :: int_kq_roots
    integer, dimension(5) :: ids1, ids2
    integer :: ir,q
    ids1(:) = 0
    ids2(:) = 0
    do ir=1,size(int_kq_roots,opdimnum)-1 
      ids1(opdimnum-1) = ir
      ids2(opdimnum-1) = ir-1
      do q=-1,q_minn,-1
        int_kq_roots(q,ids1(1),ids1(2),ids1(3),ids1(4),ids1(5)) = mpreal(1.0,kv_nwds)/kroot * &
          ( int_kq_roots(q+1,ids1(1),ids1(2),ids1(3),ids1(4),ids1(5)) - &
            int_kq_roots(q, ids2(1),ids2(2),ids2(3),ids2(4),ids2(5)))
      enddo
      do q=1,ubound(int_kq_roots,1)
        int_kq_roots(q,ids1(1),ids1(2),ids1(3),ids1(4),ids1(5)) = &
          kroot * int_kq_roots(q-1,ids1(1),ids1(2),ids1(3),ids1(4),ids1(5)) + &
            int_kq_roots(q-1, ids2(1),ids2(2),ids2(3),ids2(4),ids2(5))
      enddo
    enddo
end subroutine combine6_single_root_ints

subroutine combine5_single_root_ints(int_kq_roots,kroot,opdimnum)
    implicit none
    type(mp_complex) :: kroot
    integer :: opdimnum
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:), intent(inout) :: int_kq_roots
    integer, dimension(4) :: ids1, ids2
    integer :: ir,q
    ids1(:) = 0
    ids2(:) = 0
    do ir=1,size(int_kq_roots,opdimnum)-1 
      ids1(opdimnum-1) = ir
      ids2(opdimnum-1) = ir-1
      do q=-1,q_minn,-1
        int_kq_roots(q,ids1(1),ids1(2),ids1(3),ids1(4)) = mpreal(1.0,kv_nwds)/kroot * &
          ( int_kq_roots(q+1,ids1(1),ids1(2),ids1(3),ids1(4)) - &
            int_kq_roots(q, ids2(1),ids2(2),ids2(3),ids2(4)))
      enddo
      do q=1,ubound(int_kq_roots,1)
        int_kq_roots(q,ids1(1),ids1(2),ids1(3),ids1(4)) = &
          kroot * int_kq_roots(q-1,ids1(1),ids1(2),ids1(3),ids1(4)) + &
            int_kq_roots(q-1, ids2(1),ids2(2),ids2(3),ids2(4))
      enddo
    enddo
end subroutine combine5_single_root_ints

subroutine combine4_single_root_ints(int_kq_roots,kroot,opdimnum)
    implicit none
    type(mp_complex) :: kroot
    integer :: opdimnum
    type(mp_complex), dimension(q_minn:,0:,0:,0:), intent(inout) :: int_kq_roots
    integer, dimension(3) :: ids1, ids2
    integer :: ir,q
    ids1(:) = 0
    ids2(:) = 0
    do ir=1,size(int_kq_roots,opdimnum)-1 
      ids1(opdimnum-1) = ir
      ids2(opdimnum-1) = ir-1
      do q=-1,q_minn,-1
        int_kq_roots(q,ids1(1),ids1(2),ids1(3)) = &
            ( int_kq_roots(q+1,ids1(1),ids1(2),ids1(3)) - &
              int_kq_roots(q, ids2(1),ids2(2),ids2(3)))/kroot
      enddo
      do q=1,ubound(int_kq_roots,1)
        int_kq_roots(q,ids1(1),ids1(2),ids1(3)) = &
            kroot * int_kq_roots(q-1,ids1(1),ids1(2),ids1(3)) + &
              int_kq_roots(q-1, ids2(1),ids2(2),ids2(3))
      enddo
    enddo
end subroutine combine4_single_root_ints

! with 
! int_kq_roots(q,p1,...,pₘ,...) =
!          q
! ⌠       k
! ⎮ ───────────── f(k) dk
! ⌡  s
!   ┬─┬
!   │ │(k-ₘk)ᵖ⁽ᵐ⁾
!   ᵐ⁼¹
! having already been computed for qₘᵢₙ...qₘₐₓ , pₘ=0...max(pₘ), and m=1...s=start_idx-1, combineN_q_all_roots_ints drives computation of
!          q
! ⌠       k
! ⎮ ───────────── f(k) dk
! ⌡  M
!   ┬─┬
!   │ │(k-ₘk)ᵖ⁽ᵐ⁾
!   ᵐ⁼¹
! where M=size(sel_kroots) and ₘk=sel_kroots(m).
!                                                               ⎛┬─┬      ⎞
! Also, f(k) is any function. In practice, f(k) = 1 or f(k) = ln⎜│ │(k-kᵣ)⎟ for this program
!                                                               ⎝ ʳ       ⎠
subroutine combine4_multi_roots_ints(int_kq_roots_fk, sel_kroots, start_idx, skip_last)
  implicit none

  type(mp_complex), dimension(:) :: sel_kroots
  type(mp_complex), dimension(q_minn:,0:,0:,0:), intent(inout) :: int_kq_roots_fk
  integer :: start_idx
  logical :: skip_last

  integer :: ir

  do ir=start_idx,size(sel_kroots)
     if(skip_last.and.(ir.eq.size(sel_kroots))) then
         cycle
     endif
     call combine4_single_root_ints(int_kq_roots_fk,sel_kroots(ir),ir+1)
     if(ir.gt.1) then
         call append_root(int_kq_roots_fk,size(int_kq_roots_fk),shape(int_kq_roots_fk), sel_kroots, ir+1)
     endif
  enddo
end subroutine combine4_multi_roots_ints

subroutine combine5_multi_roots_ints(int_kq_roots_fk, sel_kroots, start_idx, skip_last)
  implicit none

  type(mp_complex), dimension(:) :: sel_kroots
  type(mp_complex), dimension(q_minn:,0:,0:,0:,0:), intent(inout) :: int_kq_roots_fk
  logical :: skip_last

  integer :: start_idx
  integer :: ir

  do ir=start_idx,size(sel_kroots)
     if(skip_last.and.(ir.eq.size(sel_kroots))) then
         cycle
     endif
     call combine5_single_root_ints(int_kq_roots_fk,sel_kroots(ir),ir+1)
     if(ir.gt.1) then
         call append_root(int_kq_roots_fk,size(int_kq_roots_fk),shape(int_kq_roots_fk), sel_kroots, ir+1)
     endif
  enddo
end subroutine combine5_multi_roots_ints

subroutine combine6_multi_roots_ints(int_kq_roots_fk, sel_kroots, start_idx, skip_last)
  implicit none

  type(mp_complex), dimension(:) :: sel_kroots
  type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(inout) :: int_kq_roots_fk
  integer :: start_idx
  logical :: skip_last

  integer :: ir

  do ir=start_idx,size(sel_kroots)
     if(skip_last.and.(ir.eq.size(sel_kroots))) then
         cycle
     endif
     call combine6_single_root_ints(int_kq_roots_fk,sel_kroots(ir),ir+1)
     if(ir.gt.1) then
         call append_root(int_kq_roots_fk,size(int_kq_roots_fk),shape(int_kq_roots_fk), sel_kroots, ir+1)
     endif
  enddo
end subroutine combine6_multi_roots_ints


! performs the following integrals:
! v_int(n,λ) = 
! ⌠    vⁿ
! ⎮ ─────── dv
! ⌡       λ
!   (v-vₒ)  ,
! where vₒ=kv_root
    subroutine do_v_int(kv_root,vel,v_int)
        implicit none
        type(mp_complex), dimension(0:,0:), intent(inout) :: v_int
        type(mp_complex) :: kv_root
        type(mp_real) :: vel
        integer :: n, lam

          do n=0,ubound(v_int,1)
            v_int(n,0) = ( mpreal(1.0,kv_nwds)/(n+1.0) ) * vel**(n+1)
          enddo

        v_int(0,1) = logw(vel-kv_root)
        do lam=2,ubound(v_int,2)
          v_int(0,lam) = mpreal(1.0,kv_nwds)/(-lam+1.0)*(vel-kv_root)**(-lam+1)
        enddo

          do lam=1,ubound(v_int,2)
            do n=1,ubound(v_int,1)
              v_int(n,lam) = kv_root * v_int(n-1,lam) + v_int(n-1,lam-1)
            enddo
          enddo
    end subroutine do_v_int


! determine which of the denominator terms
! t₁ = ω₁ - v*k₁ - 1 + i*eps
! t₂ = –ω₂ + k₂v + 1 + i*eps
! t₃ = ω₁–ω₂ - v(k₁-k₂) + i*eps,
! pass through t = i*eps for the given k integral bounds klb->kub and for each of the intervals on the velocity grid.

! vres_idxs : intervals of the velocity grid for which t₁, t₂, or t₃ pass through t = i*eps
! n_vres : number intervals in the velocity grid for which t₁, t₂, or t₃ pass through t = i*eps
subroutine is_resonant_vec(om2splcoeffs, om1, k1, klb, kub, vels, vres_idxs, n_vres)
    implicit none
    type(mp_real) :: t1_res_vel, t2_min_res_vel, t2_max_res_vel, t3_min_res_vel, t3_max_res_vel
    type(mp_real), dimension(:) :: om2splcoeffs

    type(mp_real) :: ca, cb, cc, om1, k1, klb, kub, vlb, vub
    type(mp_real) :: k2_vres_maxmin1, k2_vres_maxmin2, k2_vres_arg
    logical :: t1_resonant, t2_resonant, t3_resonant, spank1
    real, dimension(:) :: vels
    integer, dimension(:) :: vres_idxs
    integer :: iv, n_vres

    spank1 = (klb.lt.k1).and.(kub.gt.k1)

    ca = om2splcoeffs(1)
    cb = om2splcoeffs(2)
    cc = om2splcoeffs(3)

    t1_res_vel = (om1-1)/k1
    t2_min_res_vel = min(-(1-ca-cb*klb-cc*klb**2)/klb,-(1-ca-cb*kub-cc*kub**2)/kub)
    t2_max_res_vel = max(-(1-ca-cb*klb-cc*klb**2)/klb,-(1-ca-cb*kub-cc*kub**2)/kub)

    t3_min_res_vel = min((om1-ca-cb*klb-cc*klb**2)/(k1-klb),(om1-ca-cb*kub-cc*kub**2)/(k1-kub))
    t3_max_res_vel = max((om1-ca-cb*klb-cc*klb**2)/(k1-klb),(om1-ca-cb*kub-cc*kub**2)/(k1-kub))

    k2_vres_arg = cc**2*k1**2-cc*(om1-ca-k1*cb)
    if((.not.spank1).and.((k2_vres_arg).gt.0.0)) then
      k2_vres_maxmin1 = k1+1.0/cc*sqrt(k2_vres_arg)
      if((klb.lt.k2_vres_maxmin1).and.(kub.gt.k2_vres_maxmin1)) then
        t3_min_res_vel=min(t3_min_res_vel,&
          (om1-ca-cb*k2_vres_maxmin1-cc*k2_vres_maxmin1**2)/(k1-k2_vres_maxmin1))
        t3_max_res_vel=max(t3_max_res_vel,&
          (om1-ca-cb*k2_vres_maxmin1-cc*k2_vres_maxmin1**2)/(k1-k2_vres_maxmin1))
      endif
      k2_vres_maxmin2 = k1-1.0/cc*sqrt(k2_vres_arg)
      if((klb.lt.k2_vres_maxmin2).and.(kub.gt.k2_vres_maxmin2)) then
        t3_min_res_vel=min(t3_min_res_vel,&
          (om1-ca-cb*k2_vres_maxmin2-cc*k2_vres_maxmin2**2)/(k1-k2_vres_maxmin2))
        t3_max_res_vel=max(t3_max_res_vel,&
          (om1-ca-cb*k2_vres_maxmin2-cc*k2_vres_maxmin2**2)/(k1-k2_vres_maxmin2))
      endif
    endif

    n_vres = 0
    do iv=1,size(vels)-1
      vlb = mpreald(vels(iv),kv_nwds)
      vub = mpreald(vels(iv+1),kv_nwds)
    
      t1_resonant = (vlb<t1_res_vel).neqv.(vub<t1_res_vel)
      t2_resonant = .not.((vub.lt.t2_min_res_vel).or.(vlb.gt.t2_max_res_vel))
      t3_resonant = .not.((vub.lt.t3_min_res_vel).or.(vlb.gt.t3_max_res_vel))

    if( spank1 ) then
      if(t1_resonant.or.t2_resonant) then
        n_vres = n_vres + 1
        vres_idxs(n_vres) = iv
      endif
    else
      if(t1_resonant.or.t2_resonant.or.t3_resonant) then
        n_vres = n_vres + 1
        vres_idxs(n_vres) = iv
      endif
    endif
  enddo
    end subroutine is_resonant_vec


  ! with the integral described in the header of int_driver having been computed and stored in t123_int, as well as elementary perpendicular velocity integrals having been computed and stored in Ivpe, sum_gam2_is_vperp computes the full integral described in the module header.
  ! specifically, sum_gam2_is_vperp sums a list of sub-integrals whose parameters are stored in all_int_params and that all can be computed easily given results for t123_int and Ivpe.

  ! in terms of the column labels
  ! [a,b,g,f,α,β,γ,q,m,n,--,λ₁,λ₂,λ₃,X]
  ! the integrand corresponding to a given row in all_int_params is
  !       f    g    α    β    γ    q        1           m  n   dᵃ    dᵇ
  ! X (k₁) (ω₁) (c₁) (c₂) (c₃) (k₂) ────────────────── v  v   ──── ─────F(v ,v  )
  !                                     λ₁    λ₂    λ₃  ⟂  || d vᵃ d vᵇ    ⟂  ||
  !                                 (t₁)  (t₂)  (t₃)             ⟂    ||
  ! where c₁, c₂, and c₃ are spline coefficients of ω₂=c₁ + c₂ * k₂ + c₃ * pow(k₂,2).
  ! Also, as mentioned previously,
  !              5   5
  !              ⎲   ⎲                                       i   j
  ! F(v  ,v ) =  ⎳   ⎳  splcoeff4(ipara,iperp,6-i,6-j,iarb) v   v   ,
  !    ||  ⟂    i=0 j=0                                      ||  ⟂
  ! where the indices ipara and iperp depend on what grid intervals the evaluation points v   and v  lie in.  
  !                                                                                        ||      ⟂

subroutine sum_gam2_is_vperp(om1, k1, splcoeff4, om2splcoeffs, all_int_params, Ivpe, t123_int,gam2_is_ik_ik2,iarb)
  implicit none
  type(mp_real) :: om1, k1
  type(mp_real), dimension(:,:,:) :: splcoeff4
  type(mp_real), dimension(:) :: om2splcoeffs
  integer, dimension(:,:) :: all_int_params
  type(mp_real), dimension(0:,:) :: Ivpe
  type(mp_complex), dimension(0:) :: gam2_is_ik_ik2
  type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:,0:) :: t123_int
  integer :: iarb

  integer :: i_int_max, i_int
  integer :: q_int, p_int, vperp_pow
  integer :: iperp
  type(mp_complex), dimension(0:3) :: new_sum
  integer, dimension(15) :: int_params
  integer :: k_pow
  type(mp_complex) :: to_mult
  real :: fact
  external :: fact

  i_int_max = size(all_int_params,1)

  do i_int=1,i_int_max
    new_sum = mpcmplx((0.0,0.0),kv_nwds)
    int_params = all_int_params(i_int,:)

    do q_int=int_params(1),5
      do p_int=int_params(2),5
        vperp_pow = (2+q_int-int_params(1) + int_params(9))

        do iperp=1,nperp(iarb)-1  
          to_mult = &
            splcoeff4(iperp,6-p_int,6-q_int)*&
            (Ivpe(vperp_pow,iperp+1)-Ivpe(vperp_pow,iperp))*&
                  (nint(fact(q_int)/fact(q_int-int_params(1)))*&
                  nint(fact(p_int)/fact(p_int-int_params(2))))
                  
            do k_pow=0,3
              new_sum(k_pow) = new_sum(k_pow) + to_mult *&
                  t123_int(int_params(8) + k_pow - int_params(13),&
                          int_params(14),&
                          0,&
                          p_int-int_params(2)+int_params(10),&
                          int_params(12),&
                          int_params(13),&
                          int_params(14)+int_params(11))
            enddo
        enddo
      enddo
    enddo
                            
    to_mult = int_params(15) * om1**(int_params(3)-1)*k1**(int_params(4))*&
                 om2splcoeffs(1)**int_params(5)*&
                 om2splcoeffs(2)**int_params(6)*&
                 om2splcoeffs(3)**int_params(7)
    do k_pow=0,3
      gam2_is_ik_ik2(k_pow) = gam2_is_ik_ik2(k_pow) - &
        new_sum(k_pow) * to_mult
    enddo

  enddo
  end subroutine sum_gam2_is_vperp


! doM_rational_ints computes
! int_kq_roots(...,pₘ,...) = (all indices except pₘ are zero)
!     ⌠     1
!     ⎮ ────────── dk
!     ⌡ (k-ₘk)ᵖ⁽ᵐ⁾
! for each ₘk = some_kroots(m) from m=start_idx to M=size(some_kroots) and in turn for a range of pₘ specified by the dimension bounds of int_kq_roots.  Also, ₁k=some_kroots(m=1) = 0.0.
! doN_rational_ints then calls another subroutine to compute
! int_kq_roots(q,p₂,p₃,...,p )  =
!                           M
!          q
! ⌠       k
! ⎮ ───────────── dk
! ⌡  M
!   ┬─┬
!   │ │(k-ₘk)ᵖ⁽ᵐ⁾
!   ᵐ⁼²

subroutine do6_rational_ints(int_kq_roots, kb, some_kroots, start_idx, skip_last)
    implicit none
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:) :: int_kq_roots
    type(mp_real) :: kb
    type(mp_complex), dimension(:) :: some_kroots
    integer :: start_idx
    logical :: skip_last

    integer :: ir, ip, ip2
    integer, dimension(6) :: idxs


    do ir=start_idx,size(some_kroots)
      if(skip_last.and.(ir.eq.size(some_kroots))) then
          cycle
      endif
      idxs=0
      do ip=lbound(int_kq_roots,ir),ubound(int_kq_roots,ir)
        if(ir.eq.1) then
          ip2=-ip
        else
          ip2=ip
        endif
        idxs(ir)=ip
        if(ip2.eq.1) then
          int_kq_roots(idxs(1),idxs(2),idxs(3),idxs(4),idxs(5),idxs(6)) = logw(kb-some_kroots(ir))
        else
          int_kq_roots(idxs(1),idxs(2),idxs(3),idxs(4),idxs(5),idxs(6)) = mpreal(1.0,kv_nwds)/(-ip2+1.0) * (kb-some_kroots(ir))**(-ip2+1)
        endif
      enddo
    enddo

    call combine_multi_roots_ints(int_kq_roots,some_kroots(2:),max(start_idx-1,1),skip_last)

end subroutine do6_rational_ints

subroutine do5_rational_ints(int_kq_roots, kb, some_kroots, start_idx, skip_last)
    implicit none
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:) :: int_kq_roots
    type(mp_real) :: kb
    type(mp_complex), dimension(:) :: some_kroots
    integer :: start_idx
    logical :: skip_last

    integer :: ir, ip, ip2
    integer, dimension(5) :: idxs

    do ir=start_idx,size(some_kroots)
      if(skip_last.and.(ir.eq.size(some_kroots))) then
          cycle
      endif
      idxs=0
      do ip=lbound(int_kq_roots,ir),ubound(int_kq_roots,ir)
        if(ir.eq.1) then
          ip2=-ip
        else
          ip2=ip
        endif
        idxs(ir)=ip
        if(ip2.eq.1) then
          int_kq_roots(idxs(1),idxs(2),idxs(3),idxs(4),idxs(5)) = logw(kb-some_kroots(ir))
        else
          int_kq_roots(idxs(1),idxs(2),idxs(3),idxs(4),idxs(5)) = mpreal(1.0,kv_nwds)/(-ip2+1.0) * (kb-some_kroots(ir))**(-ip2+1)
        endif
      enddo
    enddo

    call combine_multi_roots_ints(int_kq_roots,some_kroots(2:),max(start_idx-1,1),skip_last)

end subroutine do5_rational_ints


    subroutine split_ln_corr(angle, corr)
        implicit none
        type(mp_real) :: angle
        integer :: corr
        type(mp_real) :: mppic

        mppic=mppi(kv_nwds)
          if( angle.gt.mppic) then
              corr = -1
          else if ( angle.le.-mppic) then
              corr = 1
          else 
              corr = 0
          endif
    end subroutine split_ln_corr


  subroutine set_eps(new_eps)
    implicit none
    type(real) :: new_eps
    eps=mpreald(new_eps,kv_nwds)
  end subroutine set_eps

type(mp_real) function arg(z)
  implicit none
  type(mp_complex) :: z
  arg = atan2w(aimag(z),mpreal(z,kv_nwds))
end function arg

! think this wrapper is needed so that the interface logw (defined above) will have a log
! fn for real numbers in addition to the one below for complex numbers
function mp_logw(ra)
  implicit none
  type(mp_real) :: mp_logw
  type(mp_real), intent(in) :: ra
  mp_logw = log(ra)
  return
end function mp_logw

! wrapper to accommodate log of negative real numbers with complex types
function mp_clogw (za)
  implicit none
  type (mp_complex):: mp_clogw
  type (mp_complex), intent (in):: za
  type(mp_real) :: mppic
  mppic=mppi(kv_nwds)
  if((aimag(za).eq.mpreal(0.0_8,kv_nwds)).and.mpreal(za,kv_nwds).lt.mpreal(0.0,kv_nwds)) then
    mp_clogw = log(abs(za)) + mppic*mpcmplx(i,kv_nwds)
  else
    mp_clogw = log(za)
  endif
  return
end function mp_clogw

! wrapper to accommodate atan2 of negative real numbers
function atan2w(ra, rb)
  implicit none
  type (mp_real):: atan2w
  type (mp_real), intent (in):: ra, rb
  type(mp_real) :: mppic

  mppic=mppi(kv_nwds)
  if(ra.eq.0) then
    if(rb.gt.0) then
      atan2w = mpreal(0.0,kv_nwds)
    else if (rb.lt.0) then
      atan2w = mppic
    else
      write(*,*) 'error: taking atan2 of 0,0'
      atan2w = mpreal(-99.0,kv_nwds)
    endif
  else
    atan2w = atan2(ra,rb)
  endif
end function atan2w

end module kv_ints_mod
