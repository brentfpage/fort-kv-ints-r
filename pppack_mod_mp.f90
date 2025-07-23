!> computes splines using arbitrary precision arithmetic.
module pppack_mod_mp
  use mpmodule
  use param_mod, only : distribution,nperp,npara,vperp,vpara,kv_nwds
  implicit none

  contains

    ! computes spline coefficients for a 2d velocity distribution function.  the pppack algorithm for a 1d spline is applied twice for generation of a 2d spline.  the approach used is the same as that used in get_splinecoeff.f90.  The purpose of this new code is to allow higher spline orders.  The number of continuous derivatives that a spline has depends on its order.  Orders as high as six may be necessary to ensure that the weak turbulence computations are smooth as functions of wavenumber.  For the pppack definition of order, adopted here, the highest power of the abscissa in a spline of order six is five.

  subroutine make_interp_spline_2d_nak_mp(iarb,splcoeff4d, k)
    implicit none
    integer :: k ! the spline order as defined in pppack
    integer :: iarb
    type(mp_real), dimension(npara(iarb),&
                             nperp(iarb)-k+1,&
                             k) :: splcoeff4a 
    type(mp_real), dimension(npara(iarb)-k+1,&
                             k,&
                             k) :: splcoeff4b
    type(mp_real), dimension(npara(iarb),&
                             k) :: splcoeff4b_helper
    type(mp_real), dimension(npara(iarb)-k+1,&
                             nperp(iarb)-k+1,&
                             k,&
                             k) :: splcoeff4c
    type(mp_real), dimension(npara(iarb)-1,&
                             nperp(iarb)-1,&
                             k,&
                             k) :: splcoeff4d


    type(mp_real), dimension(npara(iarb)+k) :: t_para
    type(mp_real), dimension(nperp(iarb)+k) :: t_perp

    type(mp_real), dimension((2*k-1)*npara(iarb)) :: q_para
    type(mp_real), dimension((2*k-1)*nperp(iarb)) :: q_perp
    type(mp_real), dimension(npara(iarb)) :: bcoef_para
    type(mp_real), dimension(nperp(iarb)) :: bcoef_perp
    type(mp_real), dimension(k,k) :: scrtch
    type(mp_real), dimension(npara(iarb)-k+2) :: break_para
    type(mp_real), dimension(nperp(iarb)-k+2) :: break_perp
    type(mp_real), dimension(k, nperp(iarb)) :: coef_perp
    type(mp_real), dimension(k, npara(iarb)) :: coef_para

    type(mp_real), dimension(npara(iarb)) :: vpara_mp
    type(mp_real), dimension(nperp(iarb)) :: vperp_mp
    type(mp_real), dimension(npara(iarb),nperp(iarb)) :: distribution_mp

    integer (kind=4) :: iflag, L
    integer (kind=4) :: j, p, p1, p2, m
    integer :: ipara, iperp
    real(kind=8) :: fact
    external fact

    do ipara=1,npara(iarb)
      vpara_mp(ipara) = mpreald(vpara(ipara,iarb),kv_nwds)
    enddo
    do iperp=1,nperp(iarb)
      vperp_mp(iperp) = mpreald(vperp(iperp,iarb),kv_nwds)
    enddo

    do ipara=1,npara(iarb)
      do iperp=1,nperp(iarb)
        distribution_mp(ipara,iperp) = mpreald(distribution(ipara,iperp,iarb),kv_nwds)
      enddo
    enddo


    ! for not-a-knot b.c.s
    t_perp(1:k) = vperp_mp(1)
    t_perp(size(t_perp)-(k-1):) = vperp_mp(nperp(iarb))
    t_perp(k+1:size(t_perp)-k) = vperp_mp(k/2+1:nperp(iarb)-k/2)

    t_para(1:k) = vpara_mp(1)
    t_para(size(t_para)-(k-1):) = vpara_mp(npara(iarb))
    t_para((k+1):size(t_para)-k) = vpara_mp(k/2+1:npara(iarb)-k/2)

    splcoeff4a=mpreald(0.0,kv_nwds)
    splcoeff4c=mpreald(0.0,kv_nwds)
    splcoeff4d=mpreald(0.0,kv_nwds)
    do ipara=1,npara(iarb)

      bcoef_perp=mpreald(0.0,kv_nwds)
      call splint_mp(vperp_mp, distribution_mp(ipara,:),&
        t_perp, nperp(iarb), k, q_perp, bcoef_perp, iflag)
      scrtch=mpreald(0.0,kv_nwds)
      call bsplpp_mp(t_perp, bcoef_perp, size(bcoef_perp), k, scrtch, break_perp, &
        & coef_perp,L)
      do p1=1,k
        splcoeff4a(ipara,:,k+1-p1) = coef_perp(p1,:size(splcoeff4a,2))
      enddo
    enddo
      ! highest power of splcoeff4a in the lowest index

    do iperp=1,size(splcoeff4a,2)
      splcoeff4b_helper = mpreald(0.0,kv_nwds)
      splcoeff4b = mpreald(0.0,kv_nwds)
      do ipara=1,npara(iarb)
        do p=0,k-1 
          do m=p,k-1
            splcoeff4b_helper(ipara,k-p) = splcoeff4b_helper(ipara,k-p) + &
               & splcoeff4a(ipara,iperp,k-m) * (-break_perp(iperp))**(m-p)/&
               & nint(fact(p)*fact(m-p))! a factor of m! is built into the def. of splcoeff4a
          enddo
        enddo
      enddo
      ! highest power of splcoeff4b_helper in the lowest index

      do p1=1,k
        bcoef_para=mpreald(0.0,kv_nwds)
        call splint_mp(vpara_mp,splcoeff4b_helper(:,p1), t_para, npara(iarb), k, q_para,&
          bcoef_para, iflag)
        scrtch=mpreald(0.0,kv_nwds)
        call bsplpp_mp(t_para, bcoef_para, size(bcoef_para), k, scrtch, break_para, &
        & coef_para,L)
        do p2=1,k
          splcoeff4b(:,k+1-p2,p1) = coef_para(p2,:size(splcoeff4b,1))

        enddo
      enddo

      do ipara=1,size(splcoeff4b,1)
        do p1=1,k
          do p=0,k-1
            do m=p,k-1
              splcoeff4c(ipara,iperp,k-p,p1) = splcoeff4c(ipara,iperp,k-p,p1) + &
                 & splcoeff4b(ipara,k-m,p1) * (-break_para(ipara))**(m-p)/&
                 & nint(fact(p)*fact(m-p))! a factor of m! is built into the def. of splcoeff4b
            enddo
          enddo
        enddo
      enddo

    enddo

! not-a-knot b.c.s: for both vperp and vpara, first two intervals are the same polynomial, and same with the last two intervals
    do p1=1,k
      do m=1,k
        do p=1,k/2-1
          do j=1,k/2-1
            splcoeff4d(p,j,p1,m) = splcoeff4c(1,1,p1,m)
            splcoeff4d(size(splcoeff4d,1)-p+1,size(splcoeff4d,2)-j+1,p1,m) = splcoeff4c(size(splcoeff4c,1),size(splcoeff4c,2),p1,m)
            splcoeff4d(p,size(splcoeff4d,2)-j+1,p1,m) = splcoeff4c(1,size(splcoeff4c,2),p1,m)
            splcoeff4d(size(splcoeff4d,1)-p+1,j,p1,m) = splcoeff4c(size(splcoeff4c,1),1,p1,m)
          enddo
          splcoeff4d(k/2:size(splcoeff4d,1)-k/2+1,p,p1,m) = splcoeff4c(:,1,p1,m)
          splcoeff4d(k/2:size(splcoeff4d,1)-k/2+1,size(splcoeff4d,2)-p+1,p1,m) = splcoeff4c(:,size(splcoeff4c,2),p1,m)

          splcoeff4d(p,k/2:size(splcoeff4d,2)-k/2+1,p1,m) = splcoeff4c(1,:,p1,m)
          splcoeff4d(size(splcoeff4d,1)-p+1,k/2:size(splcoeff4d,2)-k/2+1,p1,m) = splcoeff4c(size(splcoeff4c,1),:,p1,m)
        enddo
        ! not padding here, just copying
        splcoeff4d(k/2:size(splcoeff4d,1)-k/2+1,&
                   k/2:size(splcoeff4d,2)-k/2+1,&
                   p1,m) = splcoeff4c(:,:,p1,m)
                
      enddo
    enddo
    

  end subroutine make_interp_spline_2d_nak_mp

! computes an order 1 (quadratic) spline of y with not-a-knot boundary conditions
  subroutine make_interp_spline_quad_mp(x, y, coefc, break)
    implicit none
    integer, parameter :: k=3
    type(mp_real), dimension(:) :: x
    type(mp_real), dimension(:) :: y
    type(mp_real), dimension(size(x)+k) :: t
    type(mp_real), dimension((2*k-1)*size(x)) :: q
    type(mp_real), dimension(size(x)) :: bcoef
    type(mp_real), dimension(k,k) :: scrtch
    type(mp_real), dimension(k, size(x)) :: coef
    type(mp_real), dimension(:, :) :: coefc
    type(mp_real), dimension(:) :: break
    real :: tp1
    integer (kind=4) :: iflag, L
    integer (kind=4) :: j

    t(1:3) = x(1)
    t(size(t)-2:size(t)) = x(size(x))
    do j=4,size(t)-3
      t(j) = (x(j-1) + x(j-2))/mpreal(2.0,kv_nwds)
      tp1 = t(j)

    enddo

    call splint_mp(x, y, t, size(x), k, q, bcoef, iflag)
    call bsplpp_mp(t, bcoef, size(bcoef), k, scrtch, break, coef, L)
  ! with qi=scipy.interpolate.PPoly.from_spline(scipy.interpolate.make_interp_spline(x, y, k=2)),
  ! coef is similar to qi.c
  do j=1,L
    coefc(3,j) = coef(3,j)/mpreal(2.0,kv_nwds)
    coefc(2,j) = coef(2,j) - break(j)*coef(3,j)
    coefc(1,j) = coef(1,j) - coef(2,j)*break(j) + break(j)**2 * coef(3,j)/mpreal(2.0,kv_nwds)
  enddo

  end subroutine

  ! some subroutines from https://people.math.sc.edu/Burkardt/f_src/pppack/pppack.html
  subroutine bsplpp_mp ( t, bcoef, n, k, scrtch, break, coef, l )

  !*****************************************************************************80
  !
  !! BSPLPP converts from B-spline to piecewise polynomial form.
  !
  !  Discussion:
  !
  !    The B-spline representation of a spline is 
  !      ( T, BCOEF, N, K ),
  !    while the piecewise polynomial representation is 
  !      ( BREAK, COEF, L, K ).
  !
  !    For each breakpoint interval, the K relevant B-spline coefficients 
  !    of the spline are found and then differenced repeatedly to get the 
  !    B-spline coefficients of all the derivatives of the spline on that 
  !    interval. 
  !
  !    The spline and its first K-1 derivatives are then evaluated at the 
  !    left end point of that interval, using BSPLVB repeatedly to obtain 
  !    the values of all B-splines of the appropriate order at that point.
  !
  !  Modified:
  !
  !    14 February 2007
  !    27 April 2025, made break assumed shape, changed reals to multiple precision, hard-coded dot_product
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Carl de Boor.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Carl de Boor,
  !    A Practical Guide to Splines,
  !    Springer, 2001,
  !    ISBN: 0387953663,
  !    LC: QA1.A647.v27.
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) T(N+K), the knot sequence.
  ! 
  !    Input, real ( kind = 8 ) BCOEF(N), the B spline coefficient sequence.
  ! 
  !    Input, integer ( kind = 4 ) N, the number of B spline coefficients.
  ! 
  !    Input, integer ( kind = 4 ) K, the order of the spline.
  ! 
  !    Work array, real ( kind = 8 ) SCRTCH(K,K).
  ! 
  !    Output, real ( kind = 8 ) BREAK(L+1), the piecewise polynomial breakpoint 
  !    sequence.  BREAK contains the distinct points in the sequence T(K:N+1)
  ! 
  !    Output, real ( kind = 8 ) COEF(K,N), with COEF(I,J) = (I-1)st derivative 
  !    of the spline at BREAK(J) from the right.
  ! 
  !    Output, integer ( kind = 4 ) L, the number of polynomial pieces which 
  !    make up the spline in the interval ( T(K), T(N+1) ).
  !
    implicit none

    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    integer ( kind = 4 ) n

    type(mp_real) bcoef(n)
    type(mp_real) biatx(k)
    type(mp_real) break(:)
    type(mp_real) coef(k,n)
    type(mp_real) diff
    integer ( kind = 4 ) i
    integer ( kind = 4 ) j
    integer ( kind = 4 ) jp1
    integer ( kind = 4 ) left
    integer ( kind = 4 ) lsofar
    type(mp_real) scrtch(k,k)      
    type(mp_real) t(n+k)

    lsofar = 0
    break(1) = t(k)
    
    do left = k, n
  !
  !  Find the next nontrivial knot interval.
  !
      if ( t(left+1) == t(left) ) then
        cycle
      end if

      lsofar = lsofar + 1
      break(lsofar+1) = t(left+1)

      if ( k <= 1 ) then
        coef(1,lsofar) = bcoef(left)
        cycle
      end if
  !
  !  Store the K B-spline coefficients relevant to current knot 
  !  interval in SCRTCH(*,1).
  !
      do i = 1, k
        scrtch(i,1) = bcoef(left-k+i)
      end do
  !
  !  For J=1,...,K-1, compute the  K-J  B-spline coefficients relevant to
  !  the current knot interval for the J-th derivative by differencing
  !  those for the (J-1)st derivative, and store in SCRTCH(.,J+1).
  !
      do jp1 = 2, k
        j = jp1 - 1
        do i = 1, k - j
          diff = t(left+i) - t(left+i-(k-j))
          if ( mpreal(0.0,kv_nwds) < diff ) then
            scrtch(i,jp1) = ( ( scrtch(i+1,j) - scrtch(i,j) ) / diff ) &
              * real ( k - j, kind = 8 )
          end if
        end do
      end do
  !
  !  For J=0, ..., K-1, find the values at T(left) of the J+1
  !  B-splines of order J+1 whose support contains the current
  !  knot interval from those of order J (in  BIATX ), then combine
  !  with the B-spline coefficients (in SCRTCH(.,K-J) ) found earlier
  !  to compute the (K-J-1)st derivative at  T(LEFT) of the given
  !  spline.
  !
      call bsplvb_mp ( t, 1, 1, t(left), left, biatx )

      coef(k,lsofar) = scrtch(1,k)
      
      do jp1 = 2, k
      
        call bsplvb_mp ( t, jp1, 2, t(left), left, biatx )
        coef(k+1-jp1,lsofar) = mpreal(0.0,kv_nwds)

        do i=1,jp1
          coef(k+1-jp1,lsofar) = coef(k+1-jp1,lsofar) + biatx(i) * scrtch(i,k+1-jp1)
        enddo
        
      end do

    end do
     
    l = lsofar

    return
  end
  subroutine splint_mp ( tau, gtau, t, n, k, q, bcoef, iflag )

  !*****************************************************************************80
  !
  !! SPLINT produces the B-spline coefficients BCOEF of an interpolating spline.
  !
  !  Discussion:
  !
  !    The spline is of order K with knots T(1:N+K), and takes on the 
  !    value GTAU(I) at TAU(I), for I = 1 to N.
  !
  !    The I-th equation of the linear system 
  !      A * BCOEF = B 
  !    for the B-spline coefficients of the interpolant enforces interpolation
  !    at TAU(1:N).
  !
  !    Hence, B(I) = GTAU(I), for all I, and A is a band matrix with 2*K-1
  !    bands, if it is invertible.
  !
  !    The matrix A is generated row by row and stored, diagonal by diagonal,
  !    in the rows of the array Q, with the main diagonal going
  !    into row K.  See comments in the program.
  !
  !    The banded system is then solved by a call to BANFAC, which 
  !    constructs the triangular factorization for A and stores it again in
  !    Q, followed by a call to BANSLV, which then obtains the solution
  !    BCOEF by substitution.
  !
  !    BANFAC does no pivoting, since the total positivity of the matrix
  !    A makes this unnecessary.
  !
  !    The linear system to be solved is (theoretically) invertible if
  !    and only if
  !      T(I) < TAU(I) < TAU(I+K), for all I.
  !    Violation of this condition is certain to lead to IFLAG = 2.
  !
  !  Modified:
  !
  !    14 February 2007
  !    27 April 2025, changed reals to multiple precision
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Carl de Boor.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Carl de Boor,
  !    A Practical Guide to Splines,
  !    Springer, 2001,
  !    ISBN: 0387953663,
  !    LC: QA1.A647.v27.
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) TAU(N), the data point abscissas.  The entries in
  !    TAU should be strictly increasing.
  !
  !    Input, real ( kind = 8 ) GTAU(N), the data ordinates.
  !
  !    Input, real ( kind = 8 ) T(N+K), the knot sequence.
  !
  !    Input, integer ( kind = 4 ) N, the number of data points.
  !
  !    Input, integer ( kind = 4 ) K, the order of the spline.
  !
  !    Output, real ( kind = 8 ) Q((2*K-1)*N), the triangular factorization
  !    of the coefficient matrix of the linear system for the B-coefficients 
  !    of the spline interpolant.  The B-coefficients for the interpolant 
  !    of an additional data set can be obtained without going through all 
  !    the calculations in this routine, simply by loading HTAU into BCOEF 
  !    and then executing the call:
  !      call banslv ( q, 2*k-1, n, k-1, k-1, bcoef )
  !
  !    Output, real ( kind = 8 ) BCOEF(N), the B-spline coefficients of 
  !    the interpolant.
  !
  !    Output, integer ( kind = 4 ) IFLAG, error flag.
  !    1, = success.
  !    2, = failure.
  !
    implicit none

    integer ( kind = 4 ) n

    type(mp_real) bcoef(n)
    type(mp_real) gtau(n)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) iflag
    integer ( kind = 4 ) ilp1mx
    integer ( kind = 4 ) j
    integer ( kind = 4 ) jj
    integer ( kind = 4 ) k
    integer ( kind = 4 ) kpkm2
    integer ( kind = 4 ) left
    type(mp_real) q((2*k-1)*n)
    type(mp_real) t(n+k)
    type(mp_real) tau(n)
    type(mp_real) taui

    kpkm2 = 2 * ( k - 1 )
    left = k
    q(1:(2*k-1)*n) = mpreal(0.0,kv_nwds)
  !
  !  Loop over I to construct the N interpolation equations.
  !
    do i = 1, n
    
      taui = tau(i)
      ilp1mx = min ( i + k, n + 1 )
  !
  !  Find LEFT in the closed interval (I,I+K-1) such that
  !
  !    T(LEFT) <= TAU(I) < T(LEFT+1)
  !
  !  The matrix is singular if this is not possible.
  !
      left = max ( left, i )

      if ( taui < t(left) ) then
        iflag = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINT - Fatal Error!'
        write ( *, '(a)' ) '  The linear system is not invertible!'
        return
      end if

      do while ( t(left+1) <= taui )

        left = left + 1

        if ( left < ilp1mx ) then
          cycle
        end if

        left = left - 1

        if ( t(left+1) < taui ) then
          iflag = 2
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SPLINT - Fatal Error!'
          write ( *, '(a)' ) '  The linear system is not invertible!'
          return
        end if

        exit

      end do
  !
  !  The I-th equation enforces interpolation at TAUI, hence for all J,
  !    A(I,J) = B(J,K,T)(TAUI).
  !
  !  Only the K entries with J = LEFT-K+1,...,LEFT actually might be nonzero.
  !
  !  These K numbers are returned, in BCOEF (used for temporary storage here),
  !  by the following.
  !
      call bsplvb_mp ( t, k, 1, taui, left, bcoef )
  !
  !  We therefore want BCOEF(J) = B(LEFT-K+J)(TAUI) to go into
  !  A(I,LEFT-K+J), that is, into Q(I-(LEFT+J)+2*K,(LEFT+J)-K) since
  !  A(I+J,J) is to go into Q(I+K,J), for all I, J, if we consider Q
  !  as a two-dimensional array, with  2*K-1 rows.  See comments in
  !  BANFAC.
  !
  !  In the present program, we treat Q as an equivalent
  !  one-dimensional array, because of fortran restrictions on
  !  dimension statements.
  !
  !  We therefore want  BCOEF(J) to go into the entry of Q with index:
  !
  !    I -(LEFT+J)+2*K + ((LEFT+J)-K-1)*(2*K-1)
  !   = I-LEFT+1+(LEFT -K)*(2*K-1) + (2*K-2)*J
  !
      jj = i - left + 1 + ( left - k ) * ( k + k - 1 )

      do j = 1, k
        jj = jj + kpkm2
        q(jj) = bcoef(j)
      end do
      
    end do
  !
  !  Obtain factorization of A, stored again in Q.
  !
    call banfac_mp ( q, k+k-1, n, k-1, k-1, iflag )
    
    if ( iflag == 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINT - Fatal Error!'
      write ( *, '(a)' ) '  The linear system is not invertible!'
      return
    end if
  !
  !  Solve 
  !
  !    A * BCOEF = GTAU
  !
  !  by back substitution.
  !
    bcoef(1:n) = gtau(1:n)

    call banslv_mp ( q, k+k-1, n, k-1, k-1, bcoef )

    return
  end
  subroutine banfac_mp ( w, nroww, nrow, nbandl, nbandu, iflag )

  !*****************************************************************************80
  !
  !! BANFAC factors a banded matrix without pivoting.
  !
  !  Discussion:
  !
  !    BANFAC returns in W the LU-factorization, without pivoting, of 
  !    the banded matrix A of order NROW with (NBANDL+1+NBANDU) bands 
  !    or diagonals in the work array W.
  ! 
  !    Gauss elimination without pivoting is used.  The routine is 
  !    intended for use with matrices A which do not require row 
  !    interchanges during factorization, especially for the totally 
  !    positive matrices which occur in spline calculations.
  !
  !    The matrix storage mode used is the same one used by LINPACK 
  !    and LAPACK, and results in efficient innermost loops.
  ! 
  !    Explicitly, A has 
  ! 
  !      NBANDL bands below the diagonal
  !      1     main diagonal
  !      NBANDU bands above the diagonal
  !
  !    and thus, with MIDDLE=NBANDU+1,
  !    A(I+J,J) is in W(I+MIDDLE,J) for I=-NBANDU,...,NBANDL, J=1,...,NROW.
  !
  !    For example, the interesting entries of a banded matrix
  !    matrix of order 9, with NBANDL=1, NBANDU=2:
  !
  !      11 12 13  0  0  0  0  0  0
  !      21 22 23 24  0  0  0  0  0
  !       0 32 33 34 35  0  0  0  0
  !       0  0 43 44 45 46  0  0  0
  !       0  0  0 54 55 56 57  0  0
  !       0  0  0  0 65 66 67 68  0
  !       0  0  0  0  0 76 77 78 79
  !       0  0  0  0  0  0 87 88 89
  !       0  0  0  0  0  0  0 98 99
  !
  !    would appear in the first 1+1+2=4 rows of W as follows:
  !
  !       0  0 13 24 35 46 57 68 79
  !       0 12 23 34 45 56 67 78 89
  !      11 22 33 44 55 66 77 88 99
  !      21 32 43 54 65 76 87 98  0
  ! 
  !    All other entries of W not identified in this way with an
  !    entry of A are never referenced.
  ! 
  !    This routine makes it possible to solve any particular linear system 
  !    A*X=B for X by the call
  !
  !      call banslv ( w, nroww, nrow, nbandl, nbandu, b )
  !
  !    with the solution X contained in B on return.
  ! 
  !    If IFLAG=2, then one of NROW-1, NBANDL, NBANDU failed to be nonnegative, 
  !    or else one of the potential pivots was found to be zero 
  !    indicating that A does not have an LU-factorization.  This 
  !    implies that A is singular in case it is totally positive.
  !
  !  Modified:
  !
  !    14 February 2007
  !    27 April 2025, changed reals to multiple precision
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Carl de Boor.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Carl de Boor,
  !    A Practical Guide to Splines,
  !    Springer, 2001,
  !    ISBN: 0387953663,
  !    LC: QA1.A647.v27.
  !
  !  Parameters:
  ! 
  !    Input/output, real ( kind = 8 ) W(NROWW,NROW).
  !    On input, W contains the "interesting" part of a banded 
  !    matrix A, with the diagonals or bands of A stored in the
  !    rows of W, while columns of A correspond to columns of W. 
  !    On output, W contains the LU-factorization of A into a unit 
  !    lower triangular matrix L and an upper triangular matrix U 
  !    (both banded) and stored in customary fashion over the 
  !    corresponding entries of A.  
  !
  !    Input, integer ( kind = 4 ) NROWW, the row dimension of the work array W.
  !    NROWW must be at least NBANDL+1 + NBANDU.
  ! 
  !    Input, integer ( kind = 4 ) NROW, the number of rows in A.
  !
  !    Input, integer ( kind = 4 ) NBANDL, the number of bands of A below 
  !    the main diagonal.
  ! 
  !    Input, integer ( kind = 4 ) NBANDU, the number of bands of A above 
  !    the main diagonal.
  ! 
  !    Output, integer ( kind = 4 ) IFLAG, error flag.
  !    1, success.
  !    2, failure, the matrix was not factored.
  !
    implicit none

    integer ( kind = 4 ) nrow
    integer ( kind = 4 ) nroww

    type(mp_real) factor
    integer ( kind = 4 ) i
    integer ( kind = 4 ) iflag
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) middle
    integer ( kind = 4 ) nbandl
    integer ( kind = 4 ) nbandu
    type(mp_real) pivot
    type(mp_real) w(nroww,nrow)

    iflag = 1

    if ( nrow < 1 ) then
      iflag = 2
      return
    end if
  !
  !  W(MIDDLE,*) contains the main diagonal of A.
  !
    middle = nbandu + 1
    
    if ( nrow == 1 ) then
      if ( w(middle,nrow) == mpreal(0.0,kv_nwds) ) then
        iflag = 2
      end if
      return
    end if
  !
  !  A is upper triangular.  Check that the diagonal is nonzero.
  !
    if ( nbandl <= 0 ) then

      do i = 1, nrow-1
        if ( w(middle,i) == mpreal(0.0,kv_nwds) ) then
          iflag = 2
          return
        end if
      end do

      if ( w(middle,nrow) == mpreal(0.0,kv_nwds) ) then
        iflag = 2
      end if

      return
  !
  !  A is lower triangular.  Check that the diagonal is nonzero and
  !  divide each column by its diagonal.
  !
    else if ( nbandu <= 0 ) then

      do i = 1, nrow - 1

        pivot = w(middle,i)

        if ( pivot == mpreal(0.0,kv_nwds) ) then
          iflag = 2
          return
        end if

        do j = 1, min ( nbandl, nrow-i )
          w(middle+j,i) = w(middle+j,i) / pivot
        end do

      end do

      return

    end if
  !
  !  A is not just a triangular matrix.  
  !  Construct the LU factorization.
  !
    do i = 1, nrow - 1
  !
  !  W(MIDDLE,I) is the pivot for the I-th step.
  !
      if ( w(middle,i) == mpreal(0.0,kv_nwds) ) then
        iflag = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BANFAC - Fatal error!'
        write ( *, '(a,i8)' ) '  Zero pivot encountered in column ', i
        stop 1
      end if
  !
  !  Divide each entry in column I below the diagonal by PIVOT.
  !
      do j = 1, min ( nbandl, nrow-i )
        w(middle+j,i) = w(middle+j,i) / w(middle,i)
      end do
  !
  !  Subtract A(I,I+K)*(I-th column) from (I+K)-th column (below row I).
  !
      do k = 1, min ( nbandu, nrow-i )
        factor = w(middle-k,i+k)
        do j = 1, min ( nbandl, nrow-i )
          w(middle-k+j,i+k) = w(middle-k+j,i+k) - w(middle+j,i) * factor
        end do
      end do
   
    end do
  !
  !  Check the last diagonal entry.
  !
    if ( w(middle,nrow) == mpreal(0.0,kv_nwds) ) then
      iflag = 2
    end if

    return
    end
  subroutine banslv_mp ( w, nroww, nrow, nbandl, nbandu, b )

  !*****************************************************************************80
  !
  !! BANSLV solves a banded linear system A * X = B factored by BANFAC.
  !
  !  Modified:
  !
  !    14 February 2007
  !    27 April 2025, changed reals to multiple precision
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Carl de Boor.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Carl de Boor,
  !    A Practical Guide to Splines,
  !    Springer, 2001,
  !    ISBN: 0387953663,
  !    LC: QA1.A647.v27.
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) W(NROWW,NROW).  W contains the banded matrix,
  !    after it has been factored by BANFAC.
  !
  !    Input, integer ( kind = 4 ) NROWW, the row dimension of the work array W.
  !    NROWW must be at least NBANDL+1 + NBANDU.
  ! 
  !    Input, integer ( kind = 4 ) NROW, the number of rows in A.
  !
  !    Input, integer ( kind = 4 ) NBANDL, the number of bands of A below the 
  !    main diagonal.
  ! 
  !    Input, integer ( kind = 4 ) NBANDU, the number of bands of A above the 
  !    main diagonal.
  ! 
  !    Input/output, real ( kind = 8 ) B(NROW).
  !    On input, B contains the right hand side of the system to be solved.
  !    On output, B contains the solution, X.
  !
    implicit none

    integer ( kind = 4 ) nrow
    integer ( kind = 4 ) nroww

    type(mp_real) b(nrow)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) j
    integer ( kind = 4 ) jmax
    integer ( kind = 4 ) middle
    integer ( kind = 4 ) nbandl
    integer ( kind = 4 ) nbandu
    type(mp_real) w(nroww,nrow)

    middle = nbandu + 1

    if ( nrow == 1 ) then
      b(1) = b(1) / w(middle,1)
      return
    end if
  !
  !  Forward pass:
  !
  !  For I = 1, 2, ..., NROW-1, subtract RHS(I)*(I-th column of L) 
  !  from the right hand side, below the I-th row.
  !
    if ( 0 < nbandl ) then
      do i = 1, nrow - 1
        jmax = min ( nbandl, nrow-i )
        do j = 1, jmax
          b(i+j) = b(i+j) - b(i) * w(middle+j,i)
        end do
      end do
    end if
  !
  !  Backward pass:
  !
  !  For I=NROW, NROW-1,...,1, divide RHS(I) by 
  !  the I-th diagonal entry of U, then subtract 
  !  RHS(I)*(I-th column of U) from right hand side, above the I-th row.
  !
    do i = nrow, 2, -1
     
      b(i) = b(i) / w(middle,i)

      do j = 1, min ( nbandu, i - 1 )
        b(i-j) = b(i-j) - b(i) * w(middle-j,i)
      end do

    end do

    b(1) = b(1) / w(middle,1)

    return
  end
  subroutine bsplvb_mp ( t, jhigh, index, x, left, biatx )

  !*****************************************************************************80
  !
  !! BSPLVB evaluates B-splines at a point X with a given knot sequence.
  !
  !  Discusion:
  !
  !    BSPLVB evaluates all possibly nonzero B-splines at X of order
  !
  !      JOUT = MAX ( JHIGH, (J+1)*(INDEX-1) ) 
  !  
  !    with knot sequence T.
  ! 
  !    The recurrence relation
  ! 
  !                     X - T(I)               T(I+J+1) - X
  !    B(I,J+1)(X) = ----------- * B(I,J)(X) + --------------- * B(I+1,J)(X)
  !                  T(I+J)-T(I)               T(I+J+1)-T(I+1)
  ! 
  !    is used to generate B(LEFT-J:LEFT,J+1)(X) from B(LEFT-J+1:LEFT,J)(X)
  !    storing the new values in BIATX over the old. 
  !
  !    The facts that 
  !
  !      B(I,1)(X) = 1  if  T(I) <= X < T(I+1)
  !
  !    and that 
  !
  !      B(I,J)(X) = 0  unless  T(I) <= X < T(I+J)
  !
  !    are used. 
  !
  !    The particular organization of the calculations follows 
  !    algorithm 8 in chapter X of the text.
  !
  !  Modified:
  !
  !    14 February 2007
  !    27 April 2025, changed reals to multiple precision
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Carl de Boor.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Carl de Boor,
  !    A Practical Guide to Splines,
  !    Springer, 2001,
  !    ISBN: 0387953663,
  !    LC: QA1.A647.v27.
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) T(LEFT+JOUT), the knot sequence.  T is assumed to 
  !    be nondecreasing, and also, T(LEFT) must be strictly less than 
  !    T(LEFT+1).
  ! 
  !    Input, integer ( kind = 4 ) JHIGH, INDEX, determine the order 
  !    JOUT = max ( JHIGH, (J+1)*(INDEX-1) )  
  !    of the B-splines whose values at X are to be returned.  
  !    INDEX is used to avoid recalculations when several 
  !    columns of the triangular array of B-spline values are
  !    needed, for example, in BVALUE or in BSPLVD.
  !    If INDEX = 1, the calculation starts from scratch and the entire 
  !    triangular array of B-spline values of orders
  !    1, 2, ...,JHIGH is generated order by order, that is, 
  !    column by column.
  !    If INDEX = 2, only the B-spline values of order J+1, J+2, ..., JOUT  
  !    are generated, the assumption being that BIATX, J, 
  !    DELTAL, DELTAR are, on entry, as they were on exit 
  !    at the previous call.  In particular, if JHIGH = 0, 
  !    then JOUT = J+1, that is, just the next column of B-spline 
  !    values is generated.
  !    Warning: the restriction  JOUT <= JMAX (= 20) is
  !    imposed arbitrarily by the dimension statement for DELTAL
  !    and DELTAR, but is nowhere checked for.
  ! 
  !    Input, real ( kind = 8 ) X, the point at which the B-splines 
  !    are to be evaluated.
  ! 
  !    Input, integer ( kind = 4 ) LEFT, an integer chosen so that 
  !    T(LEFT) <= X <= T(LEFT+1).
  ! 
  !    Output, real ( kind = 8 ) BIATX(JOUT), with BIATX(I) containing the
  !    value at X of the polynomial of order JOUT which agrees 
  !    with the B-spline B(LEFT-JOUT+I,JOUT,T) on the interval 
  !    (T(LEFT),T(LEFT+1)).
  !
    implicit none

    integer ( kind = 4 ), parameter :: jmax = 20

    integer ( kind = 4 ) jhigh

    type(mp_real) biatx(jhigh)
    type(mp_real), save, dimension ( jmax ) :: deltal
    type(mp_real), save, dimension ( jmax ) :: deltar
    integer ( kind = 4 ) i
    integer ( kind = 4 ) index
    integer ( kind = 4 ), save :: j = 1
    integer ( kind = 4 ) left
    type(mp_real) saved
    type(mp_real) t(left+jhigh)
    type(mp_real) term
    type(mp_real) x

    if ( index == 1 ) then 
      j = 1
      biatx(1) = mpreal(1.0,kv_nwds)
      if ( jhigh <= j ) then
        return
      end if
    end if

    if ( t(left+1) <= t(left) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BSPLVB - Fatal error!'
      write ( *, '(a)' ) '  It is required that T(LEFT) < T(LEFT+1).'
      write ( *, '(a,i8)' ) '  But LEFT = ', left
      write ( *, '(a,g14.6)' ) '  T(LEFT) =   ', t(left)
      write ( *, '(a,g14.6)' ) '  T(LEFT+1) = ', t(left+1)
      stop 1
    end if

    do
     
      deltar(j) = t(left+j) - x
      deltal(j) = x - t(left+1-j)

      saved = mpreal(0.0,kv_nwds)
      do i = 1, j
        term = biatx(i) / ( deltar(i) + deltal(j+1-i) )
        biatx(i) = saved + deltar(i) * term
        saved = deltal(j+1-i) * term
      end do

      biatx(j+1) = saved
      j = j + 1

      if ( jhigh <= j ) then
        exit
      end if

    end do

    return
  end
end module pppack_mod_mp
