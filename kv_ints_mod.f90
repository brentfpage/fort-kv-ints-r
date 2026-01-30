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
! provided as the argument disp_deriv_ik to the subroutine gam2_is_wccs_for_ik below.
! the variables k and ω from www.github.com/brentfpage/fort-kv-ints-r/blob/main/preprint.pdf are referred to in this program as k₁ and ω₁, while the variables k' and ω' from that writeup are referred to as k₂ and ω₂.

module kv_ints_mod
  use param_mod, only : i, kv_nwds, lam1_max, lam2_max, lam3_max, q_min, q_max,&
    nperp, narb, npara_max, vpara, vperp, npara, nperp_max, q_minn, n_max, delta,&
    vperp_pow_max, kv_root_lam_max, f_spl_degr, sigma_max_standard, sigma_max_spank1,&
    q_maxx_spank1, q_maxx_standard
  use mpmodule
  implicit none
  type(mp_real) :: eps ! small imaginary part for resonant denominators

! wrapper to accommodate log of negative real numbers
  interface logw
    module procedure mp_logw
    module procedure mp_clogw
  end interface

  interface do_single_kroot_integrals
    module procedure do4_single_kroot_integrals
    module procedure do6_single_kroot_integrals
  end interface

  interface combine_multi_roots_ints
    module procedure combine3_multi_roots_ints
    module procedure combine4_multi_roots_ints
    module procedure combine5_multi_roots_ints
    module procedure combine6_multi_roots_ints
  end interface

  contains

  ! compute the spline coefficients of h(k₂) = k₂*ω₁ - k₁*ω₂ + k₁ - k₂, which is introduced by application of partial fraction
  ! decomposition operations w.r.t v_parallel to the integral described in the header of t123_int_driver. Here,
  ! ω₂ = om2splcoeffs(1) + om2splcoeffs(2) * k₂ + om2splcoeffs(3) * pow(k₂, 2).
  ! If the k integral bounds span k1, i.e., klb < k1 < kub, then k1 is a root of h(k₂) .
  ! In this case, (k₂ - k₁) gets factored out of h(k₂) .
  function compute_pfdsplcoeffs_nk(om2splcoeffs_nk, k1, om1, spank1_ik2) result(pfdsplcoeffs_nk)
    implicit none
    type(mp_real), dimension(:,:), intent(in) :: om2splcoeffs_nk
    type(mp_real), dimension(3,size(om2splcoeffs_nk,2)) :: pfdsplcoeffs_nk
    type(mp_real), intent(in) :: k1, om1
    integer :: spank1_ik2

    integer :: ip, ik2

    do ik2=1,size(om2splcoeffs_nk,2)-2
      if (ik2.eq.spank1_ik2) then
        pfdsplcoeffs_nk(1, ik2) = om1 - k1*om2splcoeffs_nk(2, ik2)-k1**2*om2splcoeffs_nk(3, ik2) - 1.0 ! k^0 term
        pfdsplcoeffs_nk(2, ik2) = -k1*om2splcoeffs_nk(3, ik2) ! k^1 term
        pfdsplcoeffs_nk(3, ik2) = mpreal(0.0,kv_nwds)
      else
        do ip=1,size(pfdsplcoeffs_nk,1)
          pfdsplcoeffs_nk(ip, ik2) = -k1*om2splcoeffs_nk(ip, ik2)
        enddo
        pfdsplcoeffs_nk(2, ik2) = pfdsplcoeffs_nk(2, ik2) + om1 - 1.0
        pfdsplcoeffs_nk(1, ik2) = pfdsplcoeffs_nk(1, ik2) + k1
      endif
    enddo
  end function compute_pfdsplcoeffs_nk

  function compute_pfd_kroots_nk(pfdsplcoeffs_nk, spank1_ik2) result(kroots_nk)
    implicit none
    type(mp_real), dimension(:,:), intent(in) :: pfdsplcoeffs_nk
    type(mp_complex), dimension(2,size(pfdsplcoeffs_nk,2)) :: kroots_nk
    type(mp_real) :: sqrt_arg
    integer :: ik2
    integer :: spank1_ik2
    do ik2=1,size(pfdsplcoeffs_nk,2)-2

      if(ik2.eq.spank1_ik2) then
        kroots_nk(1, ik2) = -pfdsplcoeffs_nk(1, ik2)/pfdsplcoeffs_nk(2, ik2)
        kroots_nk(2, ik2) = mpreal(0.0,kv_nwds) ! second pfd root is k1
      else
        sqrt_arg = pfdsplcoeffs_nk(2, ik2)**2-4.0*pfdsplcoeffs_nk(3, ik2)*pfdsplcoeffs_nk(1, ik2)
        if(sqrt_arg.ge.0.0) then
          kroots_nk(1, ik2) =  (-pfdsplcoeffs_nk(2, ik2)+sqrt(sqrt_arg)) / (2.0 * pfdsplcoeffs_nk(3, ik2)) 
          kroots_nk(2, ik2) =  (-pfdsplcoeffs_nk(2, ik2)-sqrt(sqrt_arg)) / (2.0 * pfdsplcoeffs_nk(3, ik2))
        else
          sqrt_arg = -sqrt_arg
          kroots_nk(1, ik2) =  (-pfdsplcoeffs_nk(2, ik2)+sqrt(sqrt_arg)*mpcmplx(i,kv_nwds)) / (2.0 * pfdsplcoeffs_nk(3, ik2)) 
          kroots_nk(2, ik2) =  (-pfdsplcoeffs_nk(2, ik2)-sqrt(sqrt_arg)*mpcmplx(i,kv_nwds)) / (2.0 * pfdsplcoeffs_nk(3, ik2))
        endif
      endif
    enddo
  end function compute_pfd_kroots_nk


  ! drive computation of int_kq_roots_diff(q,p₁,p₂,p₃) = 
  !              q
  ! ⌠kub        k
  ! ⎮     ───────────── dk₂
  ! ⌡klb     ₃
  !         ┬─┬
  !         │ │(k-ₘk)ᵖ⁽ᵐ⁾
  !         m=1
  ! where ₘk=kroots(m).  Minor exceptions to this description apply in the spank1 case.
  function integrate_over_v_independent_k_roots(kroots,klb,kub,spank1) result(int_kq_roots_diff)
    implicit none
    type(mp_complex), dimension(:), intent(in) :: kroots
    type(mp_complex), dimension(size(kroots)+1) :: kroots_helper
    type(mp_complex), dimension(:,:,:,:), allocatable :: int_kq_roots_diff
    type(mp_real) :: mppic
    type(mp_real) :: klb,kub
    logical :: spank1
    integer :: q_maxx, sigma_max
    integer :: ir

    if(spank1) then
      sigma_max = sigma_max_spank1
      q_maxx = q_maxx_spank1
    else
      sigma_max = sigma_max_standard
      q_maxx = q_maxx_standard
    endif

    allocate(int_kq_roots_diff(q_minn:q_maxx,0:n_max+lam3_max,0:sigma_max,0:sigma_max))

    mppic = mppi(kv_nwds)

    kroots_helper(1) = mpcmplx(cmplx(0.0,0.0),kv_nwds)
    kroots_helper(2:2+size(kroots)-1) = kroots

    do ir=1,size(kroots_helper)
      if(spank1.and.(ir.eq.size(kroots_helper))) then
        cycle
      endif
      call do_single_kroot_integrals(klb,kub,kroots_helper(ir),int_kq_roots_diff,ir)
    enddo

    if(spank1) then
  ! this integral is undefined and should be treated in principal value.
  ! this is accomplished by subtracting out this contribution from the pole at k=k1.
      int_kq_roots_diff(0,1,0,0) = int_kq_roots_diff(0,1,0,0) + mpcmplx(i,kv_nwds)*mppic
    endif
    if(spank1) then
      call combine_multi_roots_ints(int_kq_roots_diff(:,:,:,0),kroots([1,2]),1)
    else
      call combine_multi_roots_ints(int_kq_roots_diff,kroots,1)
    endif
  end function integrate_over_v_independent_k_roots


  ! computes
  ! t123_int(q,τ,σ,n,λ₁,λ₂,λ₃) = 
  !  kub     q+λ₂       λ₃     λ₁
  ! ⌠       k₂   (k₂-k₁)  (-k1)    ⌠       vⁿ
  ! ⎮ dk₂ ───────────────────────  ⎮ dv ───────   .
  ! ⌡            τ       σ       σ ⌡     ₃
  !  klb  (k₂-₁k) (k₂-₂k) (k₂-₃k)       ┬─┬  λₐ
  !                                     │ │ tₐ
  !                                     ᵃ⁼¹
  ! Here, the velocity v is specifically parallel velocity
  ! Also,
  ! t₁ = ω₁ - v*k₁ - 1 + i*eps
  ! t₂ = –ω₂ + k₂v + 1 + i*eps
  ! t₃ = ω₁–ω₂ - v(k₁-k₂) + i*eps
  ! Further, ₘk = kroots(m)

  function t123_int_driver(klb,kub,kroots,int_kq_roots_diff,&
    v_int_lam1_diff,vlb,vub,k1,om1,om2splcoeffs,pfdsplcoeffs,spank1,t2_is_res,t3_is_res) result(t123_int)
    implicit none
    type(mp_complex), dimension(:) :: kroots
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:) :: int_kq_roots_diff
    type(mp_complex), dimension(0:,0:), intent(in) :: v_int_lam1_diff
    type(mp_complex), dimension(:,:,:,:,:,:,:), allocatable :: t123_int
    type(mp_real) :: k1, om1, kub, klb, vlb, vub
    type(mp_real), dimension(:) :: om2splcoeffs,pfdsplcoeffs
    logical :: spank1, t2_is_res, t3_is_res

    integer :: lam1, lam2, lam3, qq, ip, tau, sig, n

    integer :: t13_lam1_max, t13_lam3_max
    integer :: t123_lam1_max, t123_lam2_max, t123_lam3_max
    logical :: necessary
    integer, dimension(4,3) :: hardest_ints
    integer :: sigma_max, q_maxx

    if(spank1) then
      sigma_max = sigma_max_spank1
      q_maxx = q_maxx_spank1
    else
      sigma_max = sigma_max_standard
      q_maxx = q_maxx_standard
    endif

    allocate(t123_int(q_minn:q_maxx,0:n_max+lam3_max,0:sigma_max,0:n_max,0:lam1_max,0:lam2_max,0:lam3_max))
    t123_int = mpcmplx((0.0,0.0),kv_nwds)

    ! does the integral in the function header for (λ₁>=0, λ₂=0, λ₃=0) apart from factors of k₁
    t123_int(:,:,:,:,:,0,0) = def_t1_and_t0_int_driver(int_kq_roots_diff(:,:,:,:,0,0),&
      v_int_lam1_diff, spank1)

    ! does the integral in the function header for (λ₁=0, λ₂=0, λ₃>=1)
    if(t3_is_res) then
      t123_int(:,:,:,0,0,0,1:) = reshape(flatsubtract(&
        indef_t23_int_driver(klb,kub,int_kq_roots_diff,kroots,vub,k1,om1,om2splcoeffs,.true.,spank1,lam3_max),&
        indef_t23_int_driver(klb,kub,int_kq_roots_diff,kroots,vlb,k1,om1,om2splcoeffs,.true.,spank1,lam3_max),&
        size(t123_int(:,:,:,0,0,0,1:))),&
        shape(t123_int(:,:,:,0,0,0,1:)))
      call increment_n(t123_int(:,:,:,:,0,0,:),k1,om1,om2splcoeffs,spank1,.true.)
    endif

  ! does the integral in the function header for (λ₁=0, λ₂>=1, λ₃=0)
    if(t2_is_res) then
      t123_int(:,:,:,0,0,1:,0) = reshape(flatsubtract(&
        indef_t23_int_driver(klb,kub,int_kq_roots_diff,kroots,vub,k1,om1,om2splcoeffs,.false.,spank1,lam2_max),&
        indef_t23_int_driver(klb,kub,int_kq_roots_diff,kroots,vlb,k1,om1,om2splcoeffs,.false.,spank1,lam2_max),&
        size(t123_int(:,:,:,0,0,1:,0))),&
        shape(t123_int(:,:,:,0,0,1:,0)))
      call increment_n(t123_int(:,:,:,:,0,:,0),k1,om1,om2splcoeffs,spank1,.false.)
    endif

    t13_lam1_max = 4
    t13_lam3_max = 2

    t123_lam1_max = 4
    t123_lam2_max = 3
    t123_lam3_max = 2

    hardest_ints = 0
    if(spank1) then
       hardest_ints(1,:) = [3, 2, 2]
       hardest_ints(2,:) = [3, 3, 1]
       hardest_ints(3,:) = [4, 2, 1]
    else
       hardest_ints(1,:) = [1, 2, 2]
       hardest_ints(2,:) = [1, 3, 1]
       hardest_ints(3,:) = [3, 0, 2]
       hardest_ints(4,:) = [4, 0, 1]
    endif


  ! with the integral in the function header having been computed for (λ₁>0, λ₂=0, λ₃=0), (λ₁=0, λ₂>0, λ₃=0), and (λ₁=0, λ₂=0, λ₃>0),
  ! the code below computes the general (λ₁>0, λ₂>0, λ₃>0) cases.

  ! the quadratic polynomial h(k)≡pfdsplcoeffs(ubound(pfdsplcoeffs,1)) * (k-₂k) * (k-₃k) gets introduced by partial fraction
  ! decomposition operations applied to 1/(t₁ t₂) , 1/(t₂ t₃), and 1/(t₃ t₁) . This specific expression for h(k) applies for the
  ! non-spank1 case.  In the spank1 case, k₁ is a root of h(k)=pfdsplcoeffs(ubound(pfdsplcoeffs,1)-1) * (k-k₁) * (k-₂k), which
  ! introduces the need to use if(spank1) ... else ... blocks extensively, 

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
                do sig=0,ubound(t123_int,3) - (lam1+lam2-1) 
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
                do sig=0,ubound(t123_int,3) - (lam1+lam2-1)
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
                do sig=0,ubound(t123_int,3) - (lam1 + lam3 - 1)
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
              do sig=0,ubound(t123_int,3) - (lam1+lam3-1)
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
                do sig=0,ubound(t123_int,3) - (lam2 + lam3 - 1) 
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
              do sig=0,ubound(t123_int,3) - (lam2 + lam3 - 1)
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

  end function t123_int_driver

  ! with
  !  
  ! ⌠vub       vⁿ
  ! ⎮    dv ────────
  ! ⌡vlb          λ₁
  !         (v-vₒ)
  ! having been pre-computed for vₒ=(-1+ω₁+i*eps)/k₁ in v_int_lam1_diff(n,λ₁)
  ! and
  !                 q      
  ! ⌠kub           k₂       
  ! ⎮    dk₂ ───────────── 
  ! ⌡klb      ₃            
  !          ┬─┬           
  !          │ │(k₂-ₘk)ᵖ⁽ᵐ⁾
  !          ᵐ⁼¹
  ! having been pre-computed in int_kq_roots_diff(q,p1,p2,p3), this function computes
  ! t_int(q,τ,σ,n,λ₁) = 
  !                     q                         
  ! ⌠kub               k₂              ⌠vub      vⁿ 
  ! ⎮    dk₂ ────────────────────────  ⎮   dv ──────────
  ! ⌡klb           τ        σ       σ  ⌡vlb           λ₁         
  !          (k₂-₁k) (k₂-₂k) (k₂-₃k)          (-t₁/k1)
  ! where t₁ is defined in the header of t123_int_driver,
  !  -t₁       ω₁-1+iϵ 
  !  ─── = v - ────────
  !   k₁          k₁   
  function def_t1_and_t0_int_driver(int_kq_roots_diff,v_int_lam1_diff,spank1) result(t_int)
    implicit none
    type(mp_complex), dimension(q_minn:,0:,0:,0:) :: int_kq_roots_diff
    type(mp_complex), dimension(&
      q_minn:ubound(int_kq_roots_diff,1),&
      0:ubound(int_kq_roots_diff,2),&
      0:ubound(int_kq_roots_diff,3)&
      ) :: int_kq_roots_helper

    type(mp_complex), dimension(0:,0:), intent(in) :: v_int_lam1_diff
    type(mp_complex), dimension(&
      q_minn:ubound(int_kq_roots_diff,1),&
      0:ubound(int_kq_roots_diff,2),&
      0:ubound(int_kq_roots_diff,3),&
      0:ubound(v_int_lam1_diff,1),&
      0:ubound(v_int_lam1_diff,2)&
      ) :: t_int
    integer :: sig
    logical :: spank1

    if(spank1) then
        int_kq_roots_helper = int_kq_roots_diff(:,:,:,0)
    else 
      do sig=0,ubound(int_kq_roots_diff,3)
        int_kq_roots_helper(:,:,sig) = int_kq_roots_diff(:,:,sig,sig)
      enddo
    endif

!     t_int(qq,tau,sig,n,lam1) = int_kq_roots_helper(qq,tau,sig) * v_int_lam1_diff(n,lam1)
    t_int = reshape(flatmultiply(&
      spread(spread(int_kq_roots_helper,4,size(t_int,4)),5,size(t_int,5)),&
      spread(spread(spread(v_int_lam1_diff,1,size(t_int,3)),1,size(t_int,2)),1,size(t_int,1)),&
      size(t_int)),&
      shape(t_int))

  end function def_t1_and_t0_int_driver


  ! performs the following integrals for one velocity bound
  ! t_int_at_vb(q,τ,σ,n,λ) = 
  !  if t3
  !  kub       q          λ        vb
  ! ⌠         k₂   (k₂-k₁)         ⌠     vⁿ
  ! ⎮ dk₂ ───────────────────────  ⎮ dv ──── 
  ! ⌡            τ       σ       σ ⌡      λ 
  !  klb  (k₂-₁k) (k₂-₂k) (k₂-₃k)        t₃
  !                                    
  !  else
  ! kub             q+λ           vb
  ! ⌠              k₂             ⌠     vⁿ
  ! ⎮ dk₂ ─────────────────────   ⎮ dv ────
  ! ⌡            τ       σ      σ ⌡      λ 
  ! klb   (k₂-₁k) (k₂-₂k) (k₂-₃k)       t₂
  !                                
  !                                
  ! where ₘk = kroots(m). Also, tₐ are defined in the header of t123_int_driver.
  ! Further, ω₂ = om2splcoeffs(1) + om2splcoeffs(2) * k₂ + om2splcoeffs(3) * pow(k₂, 2) .

  ! In addition,
  !                 q      
  ! ⌠kub           k₂       
  ! ⎮    dk₂ ───────────── 
  ! ⌡klb      ₃            
  !          ┬─┬           
  !          │ │(k₂-ₘk)ᵖ⁽ᵐ⁾
  !          ᵐ⁼¹
  ! have been pre-computed in int_kq_roots_diff(q,p1,p2,p3,0,0)
  function indef_t23_int_driver(klb,kub,int_kq_roots_diff,&
    kroots_in,vb,k1,om1,om2splcoeffs,t3,spank1,lam_2or3_max) result(t_int_at_vb)
    implicit none
    type(mp_real), dimension(:), intent(in) :: om2splcoeffs
    type(mp_complex), dimension(:), intent(in) :: kroots_in
    type(mp_complex), dimension(size(kroots_in)) :: kroots
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(inout) :: int_kq_roots_diff
    type(mp_complex), dimension(q_minn:ubound(int_kq_roots_diff,1),&
                                0:ubound(int_kq_roots_diff,2),&
                                0:ubound(int_kq_roots_diff,3),&
                                1:lam_2or3_max) :: t_int_at_vb
    type(mp_real) :: vb, k1, klb, kub, om1
    type(mp_complex), dimension(&
        & q_minn:ubound(int_kq_roots_diff,1), &
        & 0:ubound(int_kq_roots_diff,2), &
        & 0:ubound(int_kq_roots_diff,3), &
        & 0:ubound(int_kq_roots_diff,4)) :: int_rational_ln_k_k0_diff
    integer :: sig
    logical :: t3, spank1
    real :: fact
    external fact
    type(mp_complex), dimension(size(om2splcoeffs)) :: tsplcoeffs
    type(mp_complex), dimension(size(om2splcoeffs)-1) :: t_kroots
    type(mp_complex), dimension(1 + size(kroots) + size(t_kroots)) :: kroots_full
    integer :: lam_2or3_max

    t_int_at_vb = mpcmplx(cmplx(0.0,0.0),kv_nwds)

    kroots(:) = kroots_in(:)

    tsplcoeffs = compute_tsplcoeffs(om2splcoeffs,k1,om1,vb,t3)
    t_kroots = compute_t_kroots(tsplcoeffs,om2splcoeffs,k1,spank1,t3)

    kroots_full(1) = mpcmplx((0.0,0.0),kv_nwds)
    kroots_full(2:size(kroots)+1) = kroots
    kroots_full(size(kroots)+2:) = t_kroots

    call finish_integrals_over_k_roots(int_kq_roots_diff,klb,kub,kroots_full,size(kroots)+2,t3,spank1)

    t_int_at_vb(:,:,:,2:) = int_kq_roots_2_t_int_lam_gt1(int_kq_roots_diff, k1, tsplcoeffs, lam_2or3_max, spank1, t3)

    int_rational_ln_k_k0_diff = compute_int_rational_ln_k_k0_diff(klb,kub,int_kq_roots_diff,&
      kroots_full(:size(kroots)+1),t_kroots,t3,spank1)

    if(spank1) then ! the pfd denom h(k) only has one root
      if(t3) then ! the k=k1 root typically associated with index 2 of these arrays is not needed
        t_int_at_vb(:,0,:,1) = int_rational_ln_k_k0_diff(:,0,:,0)
      else
        t_int_at_vb(:,:,:,1) = int_rational_ln_k_k0_diff(:,:,:,0)
      endif
    else
      do sig=0,ubound(t_int_at_vb,3)
        t_int_at_vb(:,:,sig,1) = int_rational_ln_k_k0_diff(:,:,sig,sig)
      enddo
    endif
  end function indef_t23_int_driver

  ! tsplcoeffs: the k₂ spline coefficients of either 
  ! t₂(vb) = –ω₂ + k₂(vb) + 1 + i*eps or 
  ! t₃(vb) = ω₁–ω₂ - (vb)(k₁-k₂) + i*eps, depending on whether t3=.true. or .false.
  ! vb : a parallel velocity bound for the considered integral
  function compute_tsplcoeffs(om2splcoeffs,k1,om1,vb,t3) result(tsplcoeffs)
    implicit none
    type(mp_real), dimension(:), intent(in) :: om2splcoeffs
    type(mp_complex), dimension(3) :: tsplcoeffs
    type(mp_real) :: k1, om1, vb
    logical :: t3

    if(t3) then
      tsplcoeffs(3) = om1 - om2splcoeffs(1) - vb*k1 + mpcmplx(i,kv_nwds) * eps
      tsplcoeffs(2) = -om2splcoeffs(2) + vb
      tsplcoeffs(1) = -om2splcoeffs(3)
    else
      tsplcoeffs(3) = -om2splcoeffs(1) + 1.0 + mpcmplx(i,kv_nwds) * eps
      tsplcoeffs(2) = -om2splcoeffs(2) + vb
      tsplcoeffs(1) = -om2splcoeffs(3)
    endif
  end function compute_tsplcoeffs


  ! compute the roots of the spline representations of 
  ! t₂(vb) = –ω₂ + k₂(vb) + 1 + i*eps or 
  ! t₃(vb) = ω₁–ω₂ - (vb)(k₁-k₂) + i*eps, depending on whether t3=.true. or .false.
  function compute_t_kroots(tsplcoeffs,om2splcoeffs,k1,spank1,t3) result(t_kroots)
    implicit none
    type(mp_complex), dimension(:) :: tsplcoeffs
    type(mp_real), dimension(:) :: om2splcoeffs
    type(mp_complex), dimension(2) :: t_kroots
    type(mp_real) :: k1
    integer :: sign1, sign2
    type(mp_complex) :: sqrt_val
    logical :: spank1, t3

    sqrt_val = sqrt(tsplcoeffs(2)**2-4.0*tsplcoeffs(1)*tsplcoeffs(3))
    if(t3.and.spank1) then
        if (mpreal(tsplcoeffs(2)-2*om2splcoeffs(3)*k1,kv_nwds)>0.0) then
            sign1=1
        else
            sign1=-1
        endif
        if (mpreal(sqrt_val,kv_nwds)>0.0) then
            sign2=1
        else
            sign2=-1
        endif
        if(sign1.ne.sign2) then
            sqrt_val = -sqrt_val
        endif
    endif

  ! in the t3.and.spank1 case, the root k=k1-i*eps/(vb-om2splcoeffs(2)-2*om2splcoeffs(3)*k1) will be in the t_kroots(2) slot
    t_kroots(1) =  (-tsplcoeffs(2) - sqrt_val)/(2.0*tsplcoeffs(1))
    t_kroots(2) =  (-tsplcoeffs(2) + sqrt_val)/(2.0*tsplcoeffs(1))

  end function compute_t_kroots


  ! with int_kq_roots_diff(q,p₁,p₂,p₃,y₁,y₂)
  !                   q
  ! ⌠kub             k
  ! ⎮    ────────────────────────── dk
  ! ⌡klb  ₃            ₂
  !      ┬─┬          ┬─┬
  !      │ │(k-ₘk)ᵖ⁽ᵐ⁾│ │(k-ₙk)ʸ⁽ⁿ⁾
  !      ᵐ⁼¹          ⁿ⁼¹

  ! having already been computed for y₁=0,y₂=0, doing the y₁>0,y₂>0 cases below.
  ! here, ₘk=kroots(m) and ₙk=t_kroots(n) .
  ! small exceptions to this description apply if spank1
  subroutine finish_integrals_over_k_roots(int_kq_roots_diff,klb,kub,kroots_full,start_idx,t3,spank1)
    implicit none
    type(mp_complex), dimension(:), intent(in) :: kroots_full
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(inout) :: int_kq_roots_diff
    type(mp_real) :: klb, kub
    integer :: start_idx
    logical :: t3,spank1
    integer :: ir
    do ir=start_idx,size(kroots_full)
      if(spank1.and.t3.and.(ir.eq.size(kroots_full))) then
        call do_single_kroot_integrals(klb,kub,kroots_full(ir),int_kq_roots_diff(:,:,:,:,:,:1),ir)
      else
        call do_single_kroot_integrals(klb,kub,kroots_full(ir),int_kq_roots_diff,ir)
      endif
    enddo
    if(spank1) then
      if(t3) then
        call combine_multi_roots_ints(int_kq_roots_diff(:,0,:,0,:,:1),kroots_full([3,5,6]),2)
      else
        call combine_multi_roots_ints(int_kq_roots_diff(:,:,:,0,:,:),kroots_full([2,3,5,6]),3)
      endif
    else
      call combine_multi_roots_ints(int_kq_roots_diff,kroots_full(2:),4)
    endif
  end subroutine finish_integrals_over_k_roots

  ! converts values of int_kq_roots_diff to values of t_int_at_vb(...,lam) for lam>1 
  function int_kq_roots_2_t_int_lam_gt1(int_kq_roots_diff, k1, tsplcoeffs, lam_2or3_max, spank1, t3) result(t_int_at_vb)
    implicit none
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(in) :: int_kq_roots_diff
    type(mp_complex), dimension(&
      & q_minn:ubound(int_kq_roots_diff,1),&
      & 0:ubound(int_kq_roots_diff,2),&
      & 0:ubound(int_kq_roots_diff,3),&
      & 2:lam_2or3_max&
      ) :: t_int_at_vb ! dims: q, k1_root_pow, other_roots_pow, n
    type(mp_real) :: k1
    type(mp_complex), dimension(:) :: tsplcoeffs
    logical :: spank1, t3
    integer :: lam_2or3_max

    integer :: q, sig, tau, lam, lam_prime, lamlam
    real :: fact
    external :: fact

    t_int_at_vb = mpcmplx((0.0,0.0),kv_nwds)

    if(spank1) then
      if (t3) then
        do q=q_minn,ubound(t_int_at_vb,1) 
          do lam=2,ubound(t_int_at_vb,4)
            do sig=0,ubound(t_int_at_vb,3)
              t_int_at_vb(q,0,sig,lam) = int_kq_roots_diff(q,0,sig,0,lam-1,0)/&
                mpreal(-lam+1.0,kv_nwds)/tsplcoeffs(1)**(lam-1)
            enddo
          enddo
        enddo
      else
        do lam=2,ubound(t_int_at_vb,4)
          do q=q_minn,ubound(t_int_at_vb,1) - (lam-1)
            do tau=0,ubound(t_int_at_vb,2)
              do sig=0,ubound(t_int_at_vb,3)
                t_int_at_vb(q,tau,sig,lam) = int_kq_roots_diff(q+lam-1,tau,sig,0,lam-1,lam-1)/&
                  &mpreal(-lam+1.0,kv_nwds)/tsplcoeffs(1)**(lam-1)
              enddo
            enddo
          enddo
        enddo
      endif
    else
      if (t3) then 
        do lam=2,ubound(t_int_at_vb,4)
          do q=q_minn,ubound(t_int_at_vb,1)-(lam-1)
            do tau=0,ubound(t_int_at_vb,2) 
              lam_prime = lam - tau - 1
              do sig=0,ubound(int_kq_roots_diff,3)
                if (lam_prime.gt.0) then
                  do lamlam=0,lam_prime
                    t_int_at_vb(q,tau,sig,lam) = t_int_at_vb(q,tau,sig,lam) + &
        & ( fact(lam_prime)/mpreald(fact(lamlam)*fact(lam_prime-lamlam),kv_nwds) * &
        & (-k1)**(lam_prime-lamlam) / (-lam + 1) / tsplcoeffs(1)**(lam-1) ) &
        & * int_kq_roots_diff(q+lamlam,0,sig,sig,lam-1,lam-1)
                  enddo
                else
                  t_int_at_vb(q,tau,sig,lam) = (mpreal(1.0,kv_nwds)/ mpreal(-lam+1.0,kv_nwds) / tsplcoeffs(1)**(lam-1) )&
                    &*int_kq_roots_diff(q,-lam_prime,sig,sig,lam-1,lam-1)
                endif
              enddo
            enddo
          enddo
        enddo
      else
        do lam=2,ubound(t_int_at_vb,4)
          do q=q_minn,ubound(t_int_at_vb,1) - 1 * (lam-1) 
            do tau=0,ubound(int_kq_roots_diff,2)
              do sig=0,ubound(int_kq_roots_diff,3)
                t_int_at_vb(q,tau,sig,lam) = int_kq_roots_diff(q+lam-1,tau,sig,sig,lam-1,lam-1)&
                    /( (-lam+1.0)*tsplcoeffs(1)**(lam-1) )
              enddo
            enddo
          enddo
        enddo
      endif
    endif
  end function int_kq_roots_2_t_int_lam_gt1

  ! computes the integrals
  ! int_rational_ln_k_k0_diff(q,p₁, ...,pₘ) = 
  !             q
  ! ⌠kub       k          ┬─┬
  ! ⎮    ───────────── ln(│ │(k-ᵤk))dk
  ! ⌡klb ┬─┬               ᵤ
  !      │ │(k-ₘk)ᵖ⁽ᵐ⁾
  !       ᵐ
  ! for q=qₘᵢₙ...qₘₐₓ and p₁=0...max(p₁), ...,pₘ=0...max(pₘ)
  ! ᵤk = t_kroots(u)
  ! ₘk = kroots(m)
  function compute_int_rational_ln_k_k0_diff(klb,kub,int_kq_roots_diff,&
      kroots_helper,t_kroots,t3,spank1) result(int_rational_ln_k_k0_diff)
    implicit none
    type(mp_real) :: klb, kub
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(inout) :: int_kq_roots_diff
    type(mp_complex), dimension(&
        & q_minn:ubound(int_kq_roots_diff,1), &
        & 0:ubound(int_kq_roots_diff,2), &
        & 0:ubound(int_kq_roots_diff,3), &
        & 0:ubound(int_kq_roots_diff,4)) :: int_rational_ln_k_k0_diff
    type(mp_complex), dimension(:) :: kroots_helper
    type(mp_complex), dimension(:) :: t_kroots
    type(mp_real), dimension(size(kroots_helper),2) :: k_crossings
    logical :: t3,spank1
    integer :: ir

    k_crossings = compute_cut_crossings(kroots_helper,t_kroots)

    int_rational_ln_k_k0_diff = mpcmplx((0.0,0.0),kv_nwds)
    do ir=1,size(kroots_helper)
      if(&
        (spank1.and.(ir.eq.(size(kroots_helper)))).or.&
        (spank1.and.t3.and.(ir.eq.2))&
      ) then
        cycle
      endif
      call do_ln_k_k0_case1_int(int_kq_roots_diff,klb,kub,kroots_helper(ir),t_kroots,int_rational_ln_k_k0_diff,ir)
      call do_ln_k_k0_case2_int(klb,kub,kroots_helper(ir),t_kroots,&
        k_crossings(ir,:),int_rational_ln_k_k0_diff,ir,spank1)
    enddo
    int_rational_ln_k_k0_diff(0,0,0,0) = ln_k_k0_case3_int(klb,kub,t_kroots)

    if(spank1) then
      if(t3) then 
        call combine2_single_root_ints(int_rational_ln_k_k0_diff(:,0,:,0), kroots_helper(3), 2)
      else 
        call combine_multi_roots_ints(int_rational_ln_k_k0_diff(:,:,:,0), kroots_helper([2,3]), 1)
      endif
    else
      call combine_multi_roots_ints(int_rational_ln_k_k0_diff, kroots_helper(2:), 1)
    endif
  end function compute_int_rational_ln_k_k0_diff


  ! compute the integrals
  ! int_rational_ln_k_k0_diff(...,ip,...) =
  ! ⌠kub      1           ┬─┬
  ! ⎮   ────────────── ln(│ │(k-ᵤk))dk
  ! ⌡klb  (k-kroot)ⁱᵖ      ᵤ 
  ! for a range of ip specified by the dimension bounds of int_rational_ln_k_k0_diff, except skip ip=1 and ip=0.  Also, ip is in
  ! index ir of int_rational_ln_k_k0_diff, and all other indices are 0.
  ! ᵤk = t_kroots(u)
  ! int_kq_roots_diff: integrals computed by 
  subroutine do_ln_k_k0_case1_int(int_kq_roots_diff,klb,kub,kroot,t_kroots,int_rational_ln_k_k0_diff,ir)
    implicit none
    type(mp_complex) :: kroot
    type(mp_complex), dimension(2) :: t_kroots
    type(mp_real) :: klb, kub
    type(mp_real) :: mppic
    type(mp_real) :: prefactor1
    type(mp_complex) :: prefactor2_kub
    type(mp_complex) :: prefactor2_klb
    type(mp_complex), dimension(q_minn:,0:,0:,0:), intent(out) :: int_rational_ln_k_k0_diff
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(in) :: int_kq_roots_diff
    integer :: ir, iu, ip, ip2
    integer, dimension(4) :: ids2, ids3
    integer, dimension(2) :: ids1
    integer :: corr1_klb, corr1_kub
    type(mp_real) :: angle_sum
    mppic=mppi(kv_nwds)

    ! corr1 might always be zero because the two t_kroots are nearly pure real or have oppositely signed imaginary parts
    angle_sum = arg(klb-t_kroots(1)) + arg(klb-t_kroots(2))
    corr1_klb = split_ln_corr(angle_sum)  
    angle_sum = arg(kub-t_kroots(1)) + arg(kub-t_kroots(2))
    corr1_kub = split_ln_corr(angle_sum)

    ids2 = 0
    ids3 = 0
    do ip=lbound(int_rational_ln_k_k0_diff,ir),ubound(int_rational_ln_k_k0_diff,ir)
      ids2(ir) = ip
      if(ir.eq.1) then
          ip2 = ip 
          ids3(ir) = ip+1
          if(ip.eq.ubound(int_rational_ln_k_k0_diff,ir)) then
            cycle
          endif
      else
          ip2 = -ip
          ids3(ir) = ip-1
          if(ip.eq.lbound(int_rational_ln_k_k0_diff,ir)) then
            cycle
          endif
      endif
      if((ip2.eq.-1).or.(ip2.eq.0)) then
        cycle
      endif

      prefactor1 = mpreal(1.0,kv_nwds)/(ip2+1.0)
      prefactor2_kub = (kub-kroot)**(ip2+1)
      prefactor2_klb = (klb-kroot)**(ip2+1)

      do iu=1,size(t_kroots)
        ! from integration by parts
        ids1(:) = 0
        ids1(iu) = 1
        int_rational_ln_k_k0_diff(ids2(1),ids2(2),ids2(3),ids2(4)) = &
            int_rational_ln_k_k0_diff(ids2(1),ids2(2),ids2(3),ids2(4)) + &
             prefactor1 * (&
             prefactor2_kub*logw(kub-t_kroots(iu)) - prefactor2_klb*logw(klb-t_kroots(iu)) - &
             int_kq_roots_diff(ids3(1),ids3(2),ids3(3),ids3(4),ids1(1),ids1(2))&
             )
      enddo
      ! split ln correction
      int_rational_ln_k_k0_diff(ids2(1),ids2(2),ids2(3),ids2(4)) = &
          int_rational_ln_k_k0_diff(ids2(1),ids2(2),ids2(3),ids2(4)) + &
          prefactor1 *  (&
            prefactor2_kub * mpcmplx(2*corr1_kub*i,kv_nwds) -&
            prefactor2_klb * mpcmplx(2*corr1_klb*i,kv_nwds))*mppic
    enddo
  end subroutine do_ln_k_k0_case1_int

  ! compute the integrals
  ! int_rational_ln_k_k0_diff(...,1,...) =
  ! ⌠kub      1          ┬─┬
  ! ⎮   ───────────── ln(│ │(k-ᵤk))dk
  ! ⌡klb  (k-kroot)       ᵤ
  ! using spence's function. Here, 1 is in index ir of int_rational_ln_k_k0_diff, and all other indices are 0.
  ! ᵤk = t_kroots(u)
  ! k_crossings : k values where the utilized indefinite integral expression may have a brunch cut, see www.github.com/brentfpage/fort-kv-ints-r/blob/main/preprint.pdf
  subroutine do_ln_k_k0_case2_int(klb,kub,kroot,t_kroots,k_crossings,int_rational_ln_k_k0_diff,ir,spank1)
    implicit none
    type(mp_complex) :: kroot
    type(mp_complex), dimension(2) :: t_kroots
    type(mp_real), dimension(2) :: k_crossings
    type(mp_real) :: klb, kub
    type(mp_real) :: mppic
    type(mp_complex) :: spence_arg_kub
    type(mp_complex) :: spence_arg_klb
    type(mp_complex) :: spence_eval_kub
    type(mp_complex) :: spence_eval_klb
    type(mp_complex) :: spence_arg_at_k_crossing
    type(mp_complex) :: spence
    type(mp_complex), dimension(q_minn:,0:,0:,0:), intent(out) :: int_rational_ln_k_k0_diff
    integer :: ir, iu
    integer, dimension(4) :: ids2
    logical :: spank1
    integer :: sign_help
    mppic=mppi(kv_nwds)

    ids2 = 0
    if(ir.eq.1) then
      ids2(ir) = -1
    else
      ids2(ir) = 1
    endif

    do iu=1,size(t_kroots)
      spence_arg_kub = (kub - kroot)/(t_kroots(iu)-kroot)
      spence_eval_kub = spence(spence_arg_kub)
      spence_arg_klb = (klb - kroot)/(t_kroots(iu)-kroot)
      spence_eval_klb = spence(spence_arg_klb)
      if(aimag(spence_eval_kub).gt.0) then
        sign_help = -1
      else
        sign_help = 1
      endif

      int_rational_ln_k_k0_diff(ids2(1),ids2(2),ids2(3),ids2(4)) = &
        int_rational_ln_k_k0_diff(ids2(1),ids2(2),ids2(3),ids2(4)) + &
        (spence_eval_kub + logw(kub-t_kroots(iu))*logw(spence_arg_kub)) - &
        (spence_eval_klb + logw(klb-t_kroots(iu))*logw(spence_arg_klb))

      if(spank1.and.(ir.eq.2)) then
        int_rational_ln_k_k0_diff(ids2(1),ids2(2),ids2(3),ids2(4)) = &
          int_rational_ln_k_k0_diff(ids2(1),ids2(2),ids2(3),ids2(4)) &
          - mpcmplx(sign_help *i,kv_nwds)* mppic*logw(k_crossings(iu)-t_kroots(iu))
      else if ((klb.lt.k_crossings(iu)).and.(kub.gt.k_crossings(iu))) then
        spence_arg_at_k_crossing = (k_crossings(iu) - kroot)/(t_kroots(iu)-kroot)

        if(mpreal(spence_arg_at_k_crossing,kv_nwds).lt.0.0) then
          int_rational_ln_k_k0_diff(ids2(1),ids2(2),ids2(3),ids2(4)) = &
            int_rational_ln_k_k0_diff(ids2(1),ids2(2),ids2(3),ids2(4)) &
            - mpcmplx(sign_help*2*i,kv_nwds)* mppic*&
            (logw(mpreal(1.0,kv_nwds)-spence_arg_at_k_crossing) + logw(k_crossings(iu)-t_kroots(iu)))
        endif
      endif
    enddo
  end subroutine do_ln_k_k0_case2_int

  ! compute 
  ! ⌠kub    ┬─┬
  ! ⎮    ln(│ │(k-ᵤk))dk
  ! ⌡klb     ᵤ
  ! ᵤk = t_kroots(u)
  function ln_k_k0_case3_int(klb, kub, t_kroots) result(int_rational_ln_k_k0_0000)
    implicit none
    type(mp_complex) :: int_rational_ln_k_k0_0000
    type(mp_real) :: klb, kub
    type(mp_complex), dimension(:) :: t_kroots
    integer :: corr1_klb, corr1_kub
    type(mp_real) :: angle_sum

    integer :: iu
    type(mp_real) :: mppic

    ! corr1 might always be zero because the two t_kroots are nearly pure real or have oppositely signed imaginary parts
    angle_sum = arg(klb-t_kroots(1)) + arg(klb-t_kroots(2))
    corr1_klb = split_ln_corr(angle_sum)  
    angle_sum = arg(kub-t_kroots(1)) + arg(kub-t_kroots(2))
    corr1_kub = split_ln_corr(angle_sum)

    mppic=mppi(kv_nwds)
    int_rational_ln_k_k0_0000  = mpcmplx((0.0,0.0),kv_nwds)
    do iu=1,size(t_kroots)
      int_rational_ln_k_k0_0000 = int_rational_ln_k_k0_0000 + &
          (logw(kub-t_kroots(iu))*(kub-t_kroots(iu)) - (kub-t_kroots(iu))) - &
          (logw(klb-t_kroots(iu))*(klb-t_kroots(iu)) - (klb-t_kroots(iu)))
    enddo

    ! possible split ln correction
    int_rational_ln_k_k0_0000 = int_rational_ln_k_k0_0000 + mppic*(&
      kub * mpcmplx(2*i*corr1_kub,kv_nwds) -&
      klb * mpcmplx(2*i*corr1_klb,kv_nwds))
  end function ln_k_k0_case3_int

  ! determine the k_crossings where the indefinite integral involving spence's function may have
  ! a branch cut, see www.github.com/brentfpage/fort-kv-ints-r/blob/main/preprint.pdf
  function compute_cut_crossings(kroots_helper, t_kroots) result(k_crossings)
    implicit none
    type(mp_complex), dimension(:), intent(in) :: kroots_helper, t_kroots
    type(mp_real), dimension(size(kroots_helper),size(t_kroots)) :: k_crossings

    integer :: ir, j

    do j=1,size(kroots_helper)
      do ir=1,size(t_kroots)
        k_crossings(j,ir) = aimag(conjg(t_kroots(ir))*(kroots_helper(j)))/aimag( conjg(t_kroots(ir) - kroots_helper(j)) )
      enddo
    enddo
  end function compute_cut_crossings

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
        helper_splcoeffs(1) = om1-om2splcoeffs(1)
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
      helper_splcoeffs(1) = om2splcoeffs(1) - 1.0
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

  ! with the integral described in the header of t123_int_driver having been computed and stored in t123_int, as well as elementary
  ! perpendicular velocity integrals having been computed and stored in Ivpe, sum_gam2_is_wccs_over_ints_and_vperp computes the full integral described in the module header.
  ! specifically, sum_gam2_is_wccs_over_ints_and_vperp sums a list of sub-integrals whose parameters are stored in all_int_params and that all can be computed easily given results for t123_int and Ivpe.

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
  function sum_gam2_is_wccs_over_ints_and_vperp(om1, k1, splcoeff4, om2splcoeffs,&
      all_int_params, Ivpe, t123_int,iarb) result(gam2_is_wccs_ik_ik2_ipara)
    implicit none
    type(mp_real) :: om1, k1
    type(mp_real), dimension(:,:,:) :: splcoeff4
    type(mp_real), dimension(:) :: om2splcoeffs
    integer, dimension(:,:) :: all_int_params
    type(mp_real), dimension(0:,:) :: Ivpe
    type(mp_complex), dimension(0:2) :: gam2_is_wccs_ik_ik2_ipara
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:,0:) :: t123_int
    integer :: iarb

    integer :: i_int_max, i_int
    integer :: q_int, p_int, vperp_pow
    integer :: iperp
    type(mp_complex), dimension(0:2) :: new_sum
    integer, dimension(15) :: int_params
    integer :: k_pow
    type(mp_real) :: new_sum2
    type(mp_complex) :: to_mult
    real :: fact
    external :: fact

    i_int_max = size(all_int_params,1)
    gam2_is_wccs_ik_ik2_ipara = mpcmplx((0.0,0.0),kv_nwds)
    do i_int=1,i_int_max
      new_sum = mpcmplx((0.0,0.0),kv_nwds)
      int_params = all_int_params(i_int,:)

      do p_int=int_params(2),f_spl_degr-1
        new_sum2 = mpreal(0.0,kv_nwds)
        do q_int=int_params(1),f_spl_degr-1
          vperp_pow = (2+q_int-int_params(1) + int_params(9))
          do iperp=1,nperp(iarb)-1  
            new_sum2 = new_sum2 + &
              splcoeff4(iperp,f_spl_degr-p_int,f_spl_degr-q_int)*&
              (Ivpe(vperp_pow,iperp+1)-Ivpe(vperp_pow,iperp))*&
                    (nint(fact(q_int)/fact(q_int-int_params(1)))*&
                    nint(fact(p_int)/fact(p_int-int_params(2))))
          enddo
        enddo
        do k_pow=0,2
          new_sum(k_pow) = new_sum(k_pow) + new_sum2 *&
              t123_int(int_params(8) + k_pow - int_params(13),&
                      int_params(14),&
                      0,&
                      p_int-int_params(2)+int_params(10),&
                      int_params(12),&
                      int_params(13),&
                      int_params(14)+int_params(11))
        enddo
      enddo
                              
      to_mult = int_params(15)*om1**(int_params(3)-1)*k1**(int_params(4)-int_params(12))*(-1)**(int_params(12))*&
        om2splcoeffs(1)**int_params(5)*&
        om2splcoeffs(2)**int_params(6)*&
        om2splcoeffs(3)**int_params(7)
      do k_pow=0,2
        gam2_is_wccs_ik_ik2_ipara(k_pow) = gam2_is_wccs_ik_ik2_ipara(k_pow) - new_sum(k_pow) * to_mult
      enddo
    enddo
  end function sum_gam2_is_wccs_over_ints_and_vperp

  function eval_om2(om2splcoeffs, k) result(om2_eval)
    implicit none
    type(mp_real) :: om2_eval, k
    type(mp_real), dimension(3) :: om2splcoeffs
    om2_eval = om2splcoeffs(1) + om2splcoeffs(2)*k + om2splcoeffs(3)*k**2
  end function eval_om2

! compute the resonant velocity for the t2 denominator for a given value of k2
  function t2_res_vel_at_k(om2splcoeffs, om1, k1, k)
    implicit none
    type(mp_real) :: om2_eval, k, om1, k1
    type(mp_real), dimension(3) :: om2splcoeffs
    type(mp_real) :: t2_res_vel_at_k
    om2_eval = eval_om2(om2splcoeffs, k)
    t2_res_vel_at_k = -(1-om2_eval)/k
  end function t2_res_vel_at_k

! compute the resonant velocity for the t3 denominator for a given value of k2
  function t3_res_vel_at_k(om2splcoeffs, om1, k1, k2)
    implicit none
    type(mp_real) :: om2_eval, k2, om1, k1
    type(mp_real), dimension(3) :: om2splcoeffs
    type(mp_real) :: t3_res_vel_at_k
    om2_eval = eval_om2(om2splcoeffs, k2)
    t3_res_vel_at_k = (om1-om2_eval)/(k1-k2)
  end function t3_res_vel_at_k

! determine the values of k2 at which the resonant velocities for the t2 denominator reach a min or max
  function t2_k2_args_vres_extrema(om2splcoeffs, om1, k1) result(k2_vres_extrema_args)
    implicit none
    type(mp_real), dimension(3) :: om2splcoeffs
    type(mp_real) :: om1, k1
    type(mp_real), allocatable, dimension(:) :: k2_vres_extrema_args
    type(mp_real) :: k2_vres_helper
    k2_vres_helper = (om2splcoeffs(1) - 1)/om2splcoeffs(3)
    if((k2_vres_helper).gt.0.0) then
      k2_vres_extrema_args = &
        [&
        k1+1.0/om2splcoeffs(3)*sqrt(k2_vres_helper),&
        k1-1.0/om2splcoeffs(3)*sqrt(k2_vres_helper)&
        ]
    else
      allocate(k2_vres_extrema_args(0))
    endif
  end function t2_k2_args_vres_extrema

! determine the values of k2 at which the resonant velocities for the t3 denominator reach a min or max
  function t3_k2_args_vres_extrema(om2splcoeffs, om1, k1) result(k2_vres_extrema_args)
    implicit none
    type(mp_real), dimension(3) :: om2splcoeffs
    type(mp_real) :: om1, k1
    type(mp_real), allocatable, dimension(:) :: k2_vres_extrema_args
    type(mp_real) :: k2_vres_helper
    k2_vres_helper = om2splcoeffs(3)**2*k1**2 - om2splcoeffs(3)*(om1-om2splcoeffs(1)-k1*om2splcoeffs(2))
    if((k2_vres_helper).gt.0.0) then
      k2_vres_extrema_args = &
        [&
        sqrt(k2_vres_helper),&
        -sqrt(k2_vres_helper)&
        ]
    else
      allocate(k2_vres_extrema_args(0))
    endif
  end function t3_k2_args_vres_extrema

  ! determine which of the denominator terms
  ! t₁ = ω₁ - v*k₁ - 1 + i*eps
  ! t₂ = –ω₂ + k₂v + 1 + i*eps
  ! t₃ = ω₁–ω₂ - v(k₁-k₂) + i*eps,
  ! pass through t = i*eps for the given k integral bounds klb->kub and for each of the intervals on the velocity grid.

  ! vres_idxs : intervals of the velocity grid for which t₁, t₂, or t₃ pass through t = i*eps
  function is_resonant_npara(om2splcoeffs, om1, k1, klb, kub, vels) result(vres_idxs_out)
      implicit none
      type(mp_real) :: t1_res_vel, t2_min_res_vel, t2_max_res_vel, t3_min_res_vel, t3_max_res_vel
      type(mp_real), dimension(3) :: om2splcoeffs

      type(mp_real) :: om1, k1, klb, kub, vlb, vub
      logical :: t1_resonant, t2_resonant, t3_resonant, spank1
      real, dimension(:) :: vels
      integer, dimension(npara_max,4) :: vres_idxs
      integer, dimension(:,:), allocatable :: vres_idxs_out
      integer :: ipara, n_vres
      type(mp_real), dimension(2) :: min_max
      procedure (t2_res_vel_at_k), pointer :: eval_res_vel_at_k
      procedure (t2_k2_args_vres_extrema), pointer :: eval_k2_args_vres_extrema

      spank1 = (klb.lt.k1).and.(kub.gt.k1)

      t1_res_vel = (om1-1)/k1

      eval_res_vel_at_k => t2_res_vel_at_k
      eval_k2_args_vres_extrema => t2_k2_args_vres_extrema
      min_max = get_t_2or3_min_max_res_vels_for_k2(&
        om2splcoeffs, om1, k1, klb, kub, .true.,&
        eval_res_vel_at_k, eval_k2_args_vres_extrema)
      t2_min_res_vel = min_max(1)
      t2_max_res_vel = min_max(2)

      eval_res_vel_at_k => t3_res_vel_at_k
      eval_k2_args_vres_extrema => t3_k2_args_vres_extrema
      min_max = get_t_2or3_min_max_res_vels_for_k2(&
        om2splcoeffs, om1, k1, klb, kub, .not.spank1,&
        eval_res_vel_at_k, eval_k2_args_vres_extrema)
      t3_min_res_vel = min_max(1)
      t3_max_res_vel = min_max(2)

      n_vres = 0
      do ipara=1,size(vels)-1
        vlb = mpreald(vels(ipara),kv_nwds)
        vub = mpreald(vels(ipara+1),kv_nwds)
      
        t1_resonant = (vlb<t1_res_vel).neqv.(vub<t1_res_vel)
        t2_resonant = .not.((vub.lt.t2_min_res_vel).or.(vlb.gt.t2_max_res_vel))
        t3_resonant = .not.((vub.lt.t3_min_res_vel).or.(vlb.gt.t3_max_res_vel))

        if(t1_resonant.or.t2_resonant.or.t3_resonant) then
          n_vres = n_vres + 1
          vres_idxs(n_vres,1) = ipara
          if(t1_resonant) then
            vres_idxs(n_vres,2) = 1
          else
            vres_idxs(n_vres,2) = 0
          endif
          if(t2_resonant) then
            vres_idxs(n_vres,3) = 1
          else
            vres_idxs(n_vres,3) = 0
          endif
          if(t3_resonant) then
            vres_idxs(n_vres,4) = 1
          else
            vres_idxs(n_vres,4) = 0
          endif
        endif
      enddo
    allocate(vres_idxs_out(n_vres,4))
    vres_idxs_out(:,:) = vres_idxs(:n_vres,:)
  end function is_resonant_npara

  function get_t_2or3_min_max_res_vels_for_k2(om2splcoeffs, om1, k1, klb, kub, check_extrema,&
      eval_res_vel_at_k, eval_k2_args_vres_extrema) result(min_max)
    implicit none
    type(mp_real), dimension(3) :: om2splcoeffs
    type(mp_real) :: om1, k1, klb, kub
    procedure (t2_res_vel_at_k), pointer :: eval_res_vel_at_k
    procedure (t2_k2_args_vres_extrema), pointer :: eval_k2_args_vres_extrema
    type(mp_real) :: min_res_vel, max_res_vel
    type(mp_real) :: res_vel_extremum
    type(mp_real), dimension(2) :: min_max
    integer :: ie
    logical :: check_extrema

    type(mp_real), allocatable, dimension(:) :: vres_extrema_args

    min_res_vel = min(eval_res_vel_at_k(om2splcoeffs, om1, k1, klb), eval_res_vel_at_k(om2splcoeffs, om1, k1, kub))
    max_res_vel = max(eval_res_vel_at_k(om2splcoeffs, om1, k1, klb), eval_res_vel_at_k(om2splcoeffs, om1, k1, kub))

    if(check_extrema) then
      vres_extrema_args = eval_k2_args_vres_extrema(om2splcoeffs, om1, k1)

      do ie=1,size(vres_extrema_args)
        if((klb.lt.vres_extrema_args(ie)).and.(kub.gt.vres_extrema_args(ie))) then
          res_vel_extremum = eval_res_vel_at_k(om2splcoeffs, om1, k1, vres_extrema_args(ie))
          min_res_vel = min(min_res_vel, res_vel_extremum)
          max_res_vel = min(max_res_vel, res_vel_extremum)
        endif
      enddo
      deallocate(vres_extrema_args)
    endif
    min_max = [min_res_vel, max_res_vel]
  end function get_t_2or3_min_max_res_vels_for_k2

  function is_resonant_npara_nk(om2splcoeffs_nk, om1, k1, kknots, vels) result(vres_idxs_out)
    implicit none
    integer, allocatable, dimension(:,:) :: vres_idxs_out
    type(mp_real), dimension(:,:) :: om2splcoeffs_nk
    real, dimension(:) :: vels
    type(mp_real) :: om1, k1
    integer, dimension(:, :), allocatable :: vres_idxs_for_k2
    type(mp_real), dimension(:) :: kknots
    integer, dimension(size(kknots)*size(vels),5) :: vres_idxs

    integer :: ik2,n, idn
    
    n=0
    do ik2=1,size(kknots)-1
      vres_idxs_for_k2 = is_resonant_npara(om2splcoeffs_nk(:,ik2), om1, k1, kknots(ik2), kknots(ik2+1), vels)

      do idn=1,size(vres_idxs_for_k2,1)
        vres_idxs(n+idn,1) = ik2
        vres_idxs(n+idn,2) = vres_idxs_for_k2(idn,1)
        vres_idxs(n+idn,3:5) = vres_idxs_for_k2(idn,2:4)
      enddo
      n = n+size(vres_idxs_for_k2,1)
      deallocate(vres_idxs_for_k2)
    enddo

    allocate(vres_idxs_out(n,5))
    vres_idxs_out(:,:) = vres_idxs(1:n,:)
  end function is_resonant_npara_nk

  ! The function below drives computation of the integral in the module header for k₁=krange(ik) as a function of ik2, with
  ! klb=kknots(ik2) and kub=kknots(ik2+1). Computer algebra software has been used to expand the integral into many sub-integrals, the
  ! parameters of which are collected in all_int_params.  The wavenumber (k₂) and parallel velocity (vpara) part of every
  ! sub-integral has the same generic form and is computed by the function t123_int_driver.  For a given k₂ and vpara interval,
  ! integration over vperp and summation over the sub-integrals is handled by sum_gam2_is_wccs_over_ints_and_vperp.

  ! om1 : ω₁
  ! k1 : k₁
  ! v_int_lam1_diff : include some helper integrals relevant for all k₂ that are done in main.f90
  ! remaining arguments are described in the module header
  function gam2_is_wccs_for_ik(int_kq_roots_diff, int_kq_roots_diff_spank1, v_int_lam1_diff, Ivpe,&
      om2splcoeffs_nk, pfdsplcoeffs_nk, splcoeff4, kroots_nk, kknots, vpara, om1, k1, &
      all_int_params_spank1, all_int_params_standard, spank1_ik2)
    implicit none

    integer :: ik2, ipara, ii
    integer :: q_maxx, sigma_max
    type(mp_complex), dimension(:,:,:,:,:), intent(in) :: int_kq_roots_diff
    type(mp_complex), dimension(:,:,:,:), intent(in) :: int_kq_roots_diff_spank1
    type(mp_complex), dimension(:,:,:) :: v_int_lam1_diff
    type(mp_real), dimension(:,:) :: Ivpe
    type(mp_complex), dimension(:,:), intent(in) :: kroots_nk
    type(mp_real), dimension(:,:) :: om2splcoeffs_nk, pfdsplcoeffs_nk
    type(mp_real), dimension(:,:,:,:,:) :: splcoeff4
    integer, dimension(:,:), allocatable :: k2_vpara_res_idxs
    type(mp_real), dimension(:), intent(in) :: kknots
    real, dimension(:), intent(in) :: vpara
    type(mp_complex), allocatable, dimension(:,:,:,:,:,:) :: int_kq_roots_diff_ipara
    type(mp_complex), dimension(size(kknots)-1,0:2) :: gam2_is_wccs_for_ik
    type(mp_complex), allocatable, dimension(:,:) :: gam2_is_wccs_ik_ik2_ipara
    type(mp_real) :: vlb, vub, om1, k1
    type(mp_complex), allocatable, dimension(:,:,:,:,:,:,:) :: t123_int ! stores integrals over k₂ and v_parallel
    integer, dimension(:,:) :: all_int_params_spank1, all_int_params_standard
    integer, allocatable, dimension(:,:) :: all_int_params
    integer :: spank1_ik2
    integer :: iarb

    iarb=1
  ! only parallel velocity intervals that include a resonance yield a non-zero contribution to the integral in the module header.
  ! the subroutine is_resonant_npara_nk finds these intervals.
    k2_vpara_res_idxs = is_resonant_npara_nk(om2splcoeffs_nk, om1, k1, kknots, vpara)
    allocate(gam2_is_wccs_ik_ik2_ipara(size(k2_vpara_res_idxs,1),0:2))
    gam2_is_wccs_ik_ik2_ipara = mpcmplx((0.0,0.0),kv_nwds)

    !$omp parallel do private(ii,ik2,ipara,vlb,vub,t123_int,int_kq_roots_diff_ipara,all_int_params, q_maxx, sigma_max)
    do ii=1,size(k2_vpara_res_idxs,1)
      ik2 = k2_vpara_res_idxs(ii,1)
      ipara = k2_vpara_res_idxs(ii,2)

      vlb=mpreald(vpara(ipara),kv_nwds)
      vub=mpreald(vpara(ipara+1),kv_nwds)

      if(ik2.eq.spank1_ik2) then
        q_maxx = q_maxx_spank1
        sigma_max = sigma_max_spank1
      else
        q_maxx = q_maxx_standard
        sigma_max = sigma_max_standard
      endif

      allocate(int_kq_roots_diff_ipara(&
        q_minn:q_maxx,&
        0:n_max+lam3_max,&
        0:sigma_max,&
        0:sigma_max,&
        0:kv_root_lam_max-1,&
        0:kv_root_lam_max-1))

      int_kq_roots_diff_ipara = mpcmplx((0.0,0.0),kv_nwds)
      if(ik2.eq.spank1_ik2) then
        int_kq_roots_diff_ipara(:,:,:,:,0,0) = int_kq_roots_diff_spank1(:,:,:,:)
      else
        int_kq_roots_diff_ipara(:,:,:,:,0,0) = int_kq_roots_diff(ik2,:,:,:,:)
      endif

      t123_int = t123_int_driver(kknots(ik2),kknots(ik2+1),kroots_nk(:,ik2),int_kq_roots_diff_ipara,&
        v_int_lam1_diff(:,:,ipara),vlb,vub,k1,om1,om2splcoeffs_nk(:,ik2),pfdsplcoeffs_nk(:,ik2),&
        ik2.eq.spank1_ik2, k2_vpara_res_idxs(ii,4).eq.1, k2_vpara_res_idxs(ii,5).eq.1&
        )
      deallocate(int_kq_roots_diff_ipara)

      if(ik2.eq.spank1_ik2) then
        all_int_params = all_int_params_spank1
      else
        all_int_params = all_int_params_standard
      endif
      gam2_is_wccs_ik_ik2_ipara(ii,:) = sum_gam2_is_wccs_over_ints_and_vperp(om1,k1,&
        splcoeff4(ipara,:,:,:,iarb),om2splcoeffs_nk(:,ik2),all_int_params,Ivpe, t123_int,iarb)
      deallocate(all_int_params, t123_int)
    enddo
    !$omp end parallel do
    gam2_is_wccs_for_ik = mpcmplx((0.0,0.0),kv_nwds)
    do ii=1,size(k2_vpara_res_idxs,1)
      ik2 = k2_vpara_res_idxs(ii,1)
      gam2_is_wccs_for_ik(ik2,:) = flatadd(&
        gam2_is_wccs_for_ik(ik2,:),&
        gam2_is_wccs_ik_ik2_ipara(ii,:),size(gam2_is_wccs_for_ik,2))
    enddo

    deallocate(gam2_is_wccs_ik_ik2_ipara, k2_vpara_res_idxs)
  end function gam2_is_wccs_for_ik


  ! performs the following integrals:
  ! v_int(n,λ) = 
  ! ⌠    vⁿ
  ! ⎮ ─────── dv
  ! ⌡       λ
  !   (v-vₒ)  ,
  ! where vₒ=kv_root
  function do_v_int(kv_root,vel) result(v_int)
    implicit none
    type(mp_complex), dimension(0:n_max,0:lam1_max) :: v_int
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
  end function do_v_int

  function split_ln_corr(angle) result(corr)
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
  end function split_ln_corr

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

  function flatadd(vec1, vec2, n) result (vecsum)
    implicit none
    type(mp_complex), dimension(n) :: vec1
    type(mp_complex), dimension(n) :: vec2
    type(mp_complex), dimension(n) :: vecsum

    integer :: n
    integer :: i

    do i=1,size(vec1)
      vecsum(i) = vec1(i) + vec2(i)
    enddo
  end function

  function flatsubtract(vec1, vec2, n) result (vecdiff)
    implicit none
    type(mp_complex), dimension(n) :: vec1
    type(mp_complex), dimension(n) :: vec2
    integer :: n

    integer :: i

    type(mp_complex), dimension(n) :: vecdiff

    do i=1,n
      vecdiff(i) = vec1(i) - vec2(i)
    enddo
  end function

  function flatmultiply(vec1, vec2, n) result (vecprod)
    implicit none
    type(mp_complex), dimension(n) :: vec1
    type(mp_complex), dimension(n) :: vec2
    integer :: n

    integer :: i

    type(mp_complex), dimension(n) :: vecprod

    do i=1,n
      vecprod(i) = vec1(i) * vec2(i)
    enddo
  end function

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

  subroutine combine3_single_root_ints(int_kq_roots,kroot,opdimnum)
      implicit none
      type(mp_complex) :: kroot
      integer :: opdimnum
      type(mp_complex), dimension(q_minn:,0:,0:), intent(inout) :: int_kq_roots
      integer, dimension(2) :: ids1, ids2
      integer :: ir,q
      ids1(:) = 0
      ids2(:) = 0
      do ir=1,size(int_kq_roots,opdimnum)-1 
        ids1(opdimnum-1) = ir
        ids2(opdimnum-1) = ir-1
        do q=-1,q_minn,-1
          int_kq_roots(q,ids1(1),ids1(2)) = &
              ( int_kq_roots(q+1,ids1(1),ids1(2)) - &
                int_kq_roots(q, ids2(1),ids2(2)))/kroot
        enddo
        do q=1,ubound(int_kq_roots,1)
          int_kq_roots(q,ids1(1),ids1(2)) = &
              kroot * int_kq_roots(q-1,ids1(1),ids1(2)) + &
                int_kq_roots(q-1, ids2(1),ids2(2))
        enddo
      enddo
  end subroutine combine3_single_root_ints

  subroutine combine2_single_root_ints(int_kq_roots,kroot,opdimnum)
      implicit none
      type(mp_complex) :: kroot
      integer :: opdimnum
      type(mp_complex), dimension(q_minn:,0:), intent(inout) :: int_kq_roots
      integer, dimension(1) :: ids1, ids2
      integer :: ir,q
      ids1(:) = 0
      ids2(:) = 0
      do ir=1,size(int_kq_roots,opdimnum)-1 
        ids1(opdimnum-1) = ir
        ids2(opdimnum-1) = ir-1
        do q=-1,q_minn,-1
          int_kq_roots(q,ids1(1)) = &
              ( int_kq_roots(q+1,ids1(1)) - &
                int_kq_roots(q, ids2(1)))/kroot
        enddo
        do q=1,ubound(int_kq_roots,1)
          int_kq_roots(q,ids1(1)) = &
              kroot * int_kq_roots(q-1,ids1(1)) + &
                int_kq_roots(q-1, ids2(1))
        enddo
      enddo
  end subroutine combine2_single_root_ints

  ! with 
  ! int_kq_roots(q,p1,...,pₘ,...) =
  !          q
  ! ⌠       k
  ! ⎮ ───────────── f(k) dk
  ! ⌡  s
  !   ┬─┬
  !   │ │(k-ₘk)ᵖ⁽ᵐ⁾
  !   ᵐ⁼¹
  ! having already been computed for qₘᵢₙ...qₘₐₓ , pₘ=0...max(pₘ), and m=1...s=start_idx-1, combineN_multi_roots_ints drives
  ! computation of
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

  subroutine combine3_multi_roots_ints(int_kq_roots_fk, sel_kroots, start_idx)
    implicit none

    type(mp_complex), dimension(:) :: sel_kroots
    type(mp_complex), dimension(q_minn:,0:,0:), intent(inout) :: int_kq_roots_fk
    integer :: start_idx

    integer :: ir

    do ir=start_idx,size(sel_kroots)
       call combine3_single_root_ints(int_kq_roots_fk,sel_kroots(ir),ir+1)
       if(ir.gt.1) then
           call append_root(int_kq_roots_fk,size(int_kq_roots_fk),shape(int_kq_roots_fk), sel_kroots, ir+1)
       endif
    enddo
  end subroutine combine3_multi_roots_ints

  subroutine combine4_multi_roots_ints(int_kq_roots_fk, sel_kroots, start_idx)
    implicit none

    type(mp_complex), dimension(:) :: sel_kroots
    type(mp_complex), dimension(q_minn:,0:,0:,0:), intent(inout) :: int_kq_roots_fk
    integer :: start_idx

    integer :: ir

    do ir=start_idx,size(sel_kroots)
       call combine4_single_root_ints(int_kq_roots_fk,sel_kroots(ir),ir+1)
       if(ir.gt.1) then
           call append_root(int_kq_roots_fk,size(int_kq_roots_fk),shape(int_kq_roots_fk), sel_kroots, ir+1)
       endif
    enddo
  end subroutine combine4_multi_roots_ints

  subroutine combine5_multi_roots_ints(int_kq_roots_fk, sel_kroots, start_idx)
    implicit none

    type(mp_complex), dimension(:) :: sel_kroots
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:), intent(inout) :: int_kq_roots_fk

    integer :: start_idx
    integer :: ir

    do ir=start_idx,size(sel_kroots)
       call combine5_single_root_ints(int_kq_roots_fk,sel_kroots(ir),ir+1)
       if(ir.gt.1) then
           call append_root(int_kq_roots_fk,size(int_kq_roots_fk),shape(int_kq_roots_fk), sel_kroots, ir+1)
       endif
    enddo
  end subroutine combine5_multi_roots_ints

  subroutine combine6_multi_roots_ints(int_kq_roots_fk, sel_kroots, start_idx)
    implicit none

    type(mp_complex), dimension(:) :: sel_kroots
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(inout) :: int_kq_roots_fk
    integer :: start_idx

    integer :: ir

    do ir=start_idx,size(sel_kroots)
       call combine6_single_root_ints(int_kq_roots_fk,sel_kroots(ir),ir+1)
       if(ir.gt.1) then
           call append_root(int_kq_roots_fk,size(int_kq_roots_fk),shape(int_kq_roots_fk), sel_kroots, ir+1)
       endif
    enddo
  end subroutine combine6_multi_roots_ints

  ! doM_single_kroot_integrals computes
  ! int_kq_roots_diff(...,ip,...) = 
  !     ⌠kub       1
  !     ⎮    ──────────── dk
  !     ⌡klb (k-kroot)ⁱᵖ
  ! for a range of ip specified by the dimension bounds of int_kq_roots_diff.  Also, ip is in
  ! index ir of int_kq_roots_diff, and all other indices are 0.
  subroutine do4_single_kroot_integrals(klb,kub,kroot,int_kq_roots_diff,ir)
    implicit none
    type(mp_complex) :: kroot
    type(mp_real) :: klb, kub
    type(mp_complex), dimension(q_minn:,0:,0:,0:), intent(out) :: int_kq_roots_diff
    integer :: ir, ip, ip2
    integer, dimension(rank(int_kq_roots_diff)) :: idxer

    idxer=0
    do ip=lbound(int_kq_roots_diff,ir),ubound(int_kq_roots_diff,ir)
      idxer(ir) = ip
      if(kroot.eq.mpcmplx((0.0,0.0),kv_nwds)) then
        ip2=-ip
      else
        ip2=ip
      endif
      if(ip2.eq.1) then
        int_kq_roots_diff(idxer(1),idxer(2),idxer(3),idxer(4))&
          = logw(kub-kroot) - logw(klb-kroot)
      else
        int_kq_roots_diff(idxer(1),idxer(2),idxer(3),idxer(4))&
          = mpreal(1.0,kv_nwds)/(-ip2+1.0) *((kub-kroot)**(-ip2+1) - (klb-kroot)**(-ip2+1))
      endif
    enddo
  end subroutine do4_single_kroot_integrals

  subroutine do6_single_kroot_integrals(klb,kub,kroot,int_kq_roots_diff,ir)
    implicit none
    type(mp_complex) :: kroot
    type(mp_real) :: klb, kub
    type(mp_complex), dimension(q_minn:,0:,0:,0:,0:,0:), intent(out) :: int_kq_roots_diff
    integer :: ir, ip, ip2
    integer, dimension(rank(int_kq_roots_diff)) :: idxer

    idxer=0
    do ip=lbound(int_kq_roots_diff,ir),ubound(int_kq_roots_diff,ir)
      idxer(ir) = ip
      if(kroot.eq.mpcmplx((0.0,0.0),kv_nwds)) then
        ip2=-ip
      else
        ip2=ip
      endif
      if(ip2.eq.1) then
        int_kq_roots_diff(idxer(1),idxer(2),idxer(3),idxer(4),idxer(5),idxer(6))&
          = logw(kub-kroot) - logw(klb-kroot)
      else
        int_kq_roots_diff(idxer(1),idxer(2),idxer(3),idxer(4),idxer(5),idxer(6))&
          = mpreal(1.0,kv_nwds)/(-ip2+1.0) *((kub-kroot)**(-ip2+1) -  (klb-kroot)**(-ip2+1))
      endif
    enddo
  end subroutine do6_single_kroot_integrals

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
