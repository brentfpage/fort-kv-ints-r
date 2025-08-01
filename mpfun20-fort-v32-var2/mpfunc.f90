!*****************************************************************************

!  MPFUN20-Fort: A thread-safe arbitrary precision package with special functions
!  Binary-decimal, decimal-binary and I/O functions (module MPFUNC)

!  Updated: 28 Jun 2024

!  AUTHOR:
!    David H. Bailey
!    Lawrence Berkeley National Lab (retired)
!    Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2024 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs.  All basic arithmetic
!    operations and transcendental functions are supported, together with numerous
!    special functions.

!    In addition to fast execution times, one key feature of this package is a
!    100% THREAD-SAFE design, which means that user-level applications can be
!    easily converted for parallel execution, say using a threaded parallel
!    environment such as OpenMP.

!  DOCUMENTATION:
!    A detailed description of this package, and instructions for compiling
!    and testing this program on various specific systems are included in the
!    README file accompanying this package, and, in more detail, in the
!    following technical paper:

!    David H. Bailey, "MPFUN2020: A new thread-safe arbitrary precision package,"
!    available at http://www.davidhbailey.com/dhbpapers/mpfun2020.pdf.

!  DESCRIPTION OF THIS MODULE (MPFUNC):
!    This module contains subroutines for binary-decimal and decimal-binary
!    conversion, together with low-level E-format and F-format conversion, and
!    basic input and output.

module mpfunc
use mpfuna
use mpfunb

contains

subroutine mpcinp (iu, a, mpnw)

!  This inputs the MPC variable A from Fortran unit iu.

implicit none
integer, intent(in):: iu, mpnw
integer (mpiknd), intent(out):: a(0:)
integer la

la = a(0)
call mpinp (iu, a(0:), mpnw)
call mpinp (iu, a(la:), mpnw)
return
end subroutine mpcinp

subroutine mpcout (iu, n1, n2, a, mpnw)

!  This outputs the MPC variable A to Fortran unit iu, in En1.n2 format.

implicit none
integer, intent(in):: iu, n1, n2, mpnw
integer (mpiknd), intent(in):: a(0:)
integer la

la = a(0)
call mpout (iu, n1, n2, a(0:), mpnw)
call mpout (iu, n1, n2, a(la:), mpnw)
return
end

subroutine mpctomp (a, n, b, mpnw)

!  Converts the character(1) array A of length N into the MPR number B.
!  Restrictions: (a) no embedded blanks; (b) a leading digit (possibly
!  zero) must be present; and (c) a period must be present.  An exponent
!  (with "d" or "e") may optionally follow the numeric value.

implicit none
character(1), intent(in):: a(n)
integer, intent(in):: n, mpnw
integer (mpiknd), intent(out):: b(0:)
integer, parameter:: lexpmx = 9
character(10), parameter:: digits = '0123456789'
real (mpdknd), parameter::  d10w = 10.d0**mpndpw
integer i, iexp, ix, j, kde, kend, kexpend, kexpst, kexpsgn, knumend1, &
  knumend2, knumst1, knumst2, kper, ksgn, kstart, lexp, lnum, &
  lnum1, lnum2, mpnw1, n1, n2
real (mpdknd) t1
character(32) ca
integer (mpiknd) f(0:8), s0(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. b(0) < mpnw + 6) then
 write (mpldb, 1)
1 format ('*** MPCTOMP: uninitialized or inadequately sized arrays')
  call mpabrt ( 301)
endif

s0(0) = mpnw + 7
s1(0) = mpnw + 7
s2(0) = mpnw + 7
f(0) = 9
f(1) = mpnw

do i = 2, 8
  f(i) = 0
enddo

mpnw1 = mpnw + 1
kde = 0
kend = 0
kexpend = 0
kexpst = 0
kexpsgn = 0
knumend1 = 0
knumend2 = 0
knumst1 = 0
knumst2 = 0
kper = 0
ksgn = 0
kstart = 0

!   Locate:
!     kstart = index of first nonblank character.
!     kend = index of last nonblank character.

do i = 1, n
  if (a(i) /= ' ') goto 100
enddo

!   Input is completely blank.

write (6, 2) 1
2 format ('*** MPCTOMP: Syntax error in input string; code =',i4/ &
  'Restrictions: (a) no embedded blanks; (b) a leading digit (possibly'/ &
  'zero) must be present; and (c) a period must be present.  An exponent'/ &
  '(with "d" or "e") may optionally follow the numeric value.')
call mpabrt ( 302)

100 continue

kstart = i

do i = n, kstart, -1
  if (a(i) /= ' ') goto 110
enddo

i = kstart

110 continue

kend = i

!   Scan input for:
!     kde = index of 'd' or 'e'.
!     kexpend = index of end of exponent.
!     kexpst = index of start of exponent.
!     kespsgn = index of sign of exponent.
!     knumend1 = index of end of numeric part prior to period.
!     knumend2 = index of end of numeric part after period.
!     knumst1 = index of start of numeric part prior to period.
!     knumst2 = index of start of numeric part after period.
!     kper = index of period.
!     ksgn = index of sign of number.

do i = kstart, kend
  if (a(i) == ' ') then
    write (6, 2) 2
    call mpabrt ( 303)
  elseif (a(i) == '+' .or. a(i) == '-') then
    if (i == kstart) then
      ksgn = i
    elseif (kde > 0 .and. kexpsgn == 0 .and. kexpst == 0 .and. i < kend) then
      kexpsgn = i
    else
      write (6, 2) 3
      call mpabrt ( 304)
    endif
  elseif (a(i) == 'e' .or. a(i) == 'E' .or. a(i) == 'd' .or. a(i) == 'D') then
    if (kde == 0 .and. kper > 0 .and. i < kend) then
      kde = i
      knumend2 = i - 1
    else
      write (6, 2) 4
      call mpabrt ( 305)
    endif
  elseif (a(i) == '.') then
    if (kper == 0 .and. kde == 0 .and. knumst1 > 0 .and. knumst2 == 0) then
      kper = i
      knumend1 = i - 1
    else
      write (6, 2) 5
      call mpabrt ( 306)
    endif
  elseif (index (digits, a(i)) > 0) then
    if (knumst1 == 0) then
      knumst1 = i
    elseif (kper > 0 .and. knumst2 == 0 .and. kde ==  0) then
      knumst2 = i
    elseif (kde > 0 .and. kexpst == 0) then
      kexpst = i
    endif
    if (i == kend) then
      if (knumst2 > 0 .and. kde == 0) then
        knumend2 = i
      elseif (kexpst > 0) then
        kexpend = i
      else
        write (6, 2) 6
        call mpabrt ( 307)
      endif
    endif
  else
    write (6, 2) 7
    call mpabrt ( 308)
  endif
enddo

!   Decode exponent.

if (kexpst > 0) then
  lexp = kexpend - kexpst + 1
  if (lexp > lexpmx) then
    write (6, 3)
3   format ('*** MPCTOMP: exponent string is too long.')
    call mpabrt ( 309)
  endif

  do i = 1, lexp
    ca(i:i) = a(i+kexpst-1)
  enddo

  iexp = mpdigin (ca, lexp)
  if (kexpsgn > 0) then
    if (a(kexpsgn) == '-') iexp = -iexp
  endif
else
  iexp = 0
endif

!   Determine lengths of two sections of number.

lnum1 = knumend1 - knumst1 + 1
if (knumst2 > 0) then
  lnum2 = knumend2 - knumst2 + 1
else
  lnum2 = 0
endif
lnum = lnum1 + lnum2

!   Determine the number of chunks of digits and the left-over.

n1 = lnum / mpndpw
n2 = mod (lnum, mpndpw)

!   Construct first (left-over) portion, right-justified in CA.

ca(1:mpndpw - n2) = ' '
ix = knumst1 - 1

do i = 1, n2
  ix = ix + 1
  if (ix == kper) ix = ix + 1
  ca(i+mpndpw-n2:i+mpndpw-n2) = a(ix)
enddo

t1 = mpdigin (ca, mpndpw)
if (t1 > 0) then
  f(2) = 1
  f(3) = 0
  f(4) = t1
else
  f(2) = 0
  f(3) = 0
  f(4) = 0
endif
call mpeq (f, s0, mpnw1)

!   Process remaining chunks of digits.

do j = 1, n1
  do i = 1, mpndpw
    ix = ix + 1
    if (ix == kper) ix = ix + 1
    ca(i:i) = a(ix)
  enddo

  t1 = mpdigin (ca, mpndpw)
  if (t1 > 0) then
    f(2) = 1
    f(3) = 0
    f(4) = t1
  else
    f(2) = 0
    f(3) = 0
    f(4) = 0
  endif

  call mpmuld (s0, d10w, s1, mpnw1)
  call mpadd (s1, f, s0, mpnw1)
enddo

!  Correct exponent.

iexp = iexp - lnum2
f(2) = 1
f(3) = 0
f(4) = 10
call mpnpwr (f, iexp, s1, mpnw1)
call mpmul (s0, s1, s2, mpnw1)
if (ksgn > 0) then
  if (a(ksgn) == '-') s2(2) = -s2(2)
endif

!   Restore original precision and exit.

call mproun (s2, mpnw)
call mpeq (s2, b, mpnw)

return
end subroutine mpctomp

real (mpdknd) function mpdigin (ca, n)

!   This converts the string CA of nonblank length N to double precision.
!   CA may only be modest length and may only contain digits.  Blanks are ignored.
!   This is intended for internal use only.

  implicit none
  integer, intent(in):: n
  character(*), intent(in):: ca
  character(10), parameter:: digits = '0123456789'
  integer i, k
  real (mpdknd) d1

! End of declaration

  d1 = 0.d0

  do i = 1, n
    if (ca(i:i) /= ' ') then
      k = index (digits, ca(i:i)) - 1
      if (k < 0) then
        write (mpldb, 1) ca(i:i)
1       format ('*** MPDIGIN: non-digit in character string = ',a)
        call mpabrt ( 310)
      elseif (k <= 9) then
        d1 = 10.d0 * d1 + k
      endif
    endif
  enddo

  mpdigin = d1
end function mpdigin

character(32) function mpdigout (a, n)

!   This converts the double precision input A to a character(32) string of
!   nonblank length N.  A must be a whole number, and N must be sufficient
!   to hold it.  This is intended for internal use only.

  implicit none
  integer, intent(in):: n
  real (mpdknd), intent(in):: a
  character(10), parameter:: digits = '0123456789'
  integer i, k
  real (mpdknd) d1, d2
  character(32) ca

! End of declaration

  ca = ' '
  d1 = abs (a)

  do i = n, 1, -1
    d2 = aint (d1 / 10.d0)
    k = 1.d0 + (d1 - 10.d0 * d2)
    d1 = d2
    ca(i:i) = digits(k:k)
  enddo

  mpdigout = ca
  return
end function mpdigout

subroutine mpeformat (a, nb, nd, b, mpnw)

!   Converts the MPR number A into character form in the character(1) array B.
!   NB (input) is the length of the output string, and ND (input) is the
!   number of digits after the decimal point.  The format is analogous to
!   Fortran E format.  The result is left-justified among the NB cells of B.
!   The condition NB >= ND + 10 must hold or an error message will result.
!   NB cells must be available in array B.

implicit none
integer (mpiknd), intent(in):: a(0:)
integer, intent(in):: mpnw, nb, nd
integer i, ia, ix, ixp, i1, i2, j, k, mpnw1, na, nexp, nl
character(1), intent(out):: b(nb)
character(1) b2(nb+50)
character(10), parameter:: digits = '0123456789'
character(32) ca
real (mpdknd) aa, an, t1
real (mpdknd), parameter:: d10w = 10.d0**mpndpw
integer (mpiknd) f(0:8), s0(0:mpnw+6), s1(0:mpnw+6)

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. nb < nd + 10) then
  write (mpldb, 1)
1 format ('*** MPEFORMAT: uninitialized or inadequately sized arrays')
  call mpabrt ( 311)
endif

ia = sign (int (1, mpiknd), a(2))
na = min (int (abs (a(2))), mpnw)
s0(0) = mpnw + 7
s1(0) = mpnw + 7
mpnw1 = mpnw + 1

!   Set f = 10.

f(0) = 9
f(1) = mpnw1
f(2) = 1
f(3) = 0
f(4) = 10
f(5) = 0
f(6) = 0

!   Determine power of ten for exponent, and scale input to within 1 and 10.

if (na > 0) then
  aa = a(4)
  if (na >= 2) aa = aa + dble (a(5)) / mpbdx
  t1 = log10 (2.d0) * mpnbt * a(3) + log10 (aa)

  if (t1 >= 0.d0) then
    nexp = t1
  else
    nexp = t1 - 1.d0
  endif

  if (nexp == 0) then
    call mpeq (a, s1, mpnw1)
  elseif (nexp > 0) then
    call mpnpwr (f, nexp, s0, mpnw1)
    call mpdiv (a, s0, s1, mpnw1)
  elseif (nexp < 0) then
    call mpnpwr (f, - nexp, s0, mpnw1)
    call mpmul (a, s0, s1, mpnw1)
  endif

!   If we didn't quite get it exactly right, multiply or divide by 10 to fix.

100 continue

  if (s1(3) < 0) then
    nexp = nexp - 1
    call mpmuld (s1, 10.d0, s0, mpnw1)
    call mpeq (s0, s1, mpnw1)
    goto 100
  elseif (s1(4) >= 10) then
    nexp = nexp + 1
    call mpdivd (s1, 10.d0, s0, mpnw1)
    call mpeq (s0, s1, mpnw1)
    goto 100
  endif

  s1(2) = abs (s1(2))
else
  nexp = 0
  call mpeq (a, s1, mpnw1)
endif

!   Insert sign and first digit.

ix = 0
if (ia == -1) then
  ix = ix + 1
  b2(ix) = '-'
endif
if (na > 0) then
  an = s1(4)
else
  an = 0.d0
endif
ca = mpdigout (an, 1)
ix = ix + 1
b2(ix) = ca(1:1)
ix = ix + 1
b2(ix) = '.'
ixp = ix

!   Set f = an.

f(0) = 9
f(1) = mpnw1
f(2) = 1
f(3) = 0
f(4) = an
f(5) = 0
f(6) = 0
call mpsub (s1, f, s0, mpnw1)
call mpmuld (s0, d10w, s1, mpnw1)

!   Calculate the number of remaining chunks.

nl = nd / mpndpw + 1

!   Insert the digits of the remaining words.

do j = 1, nl
  if (s1(2) /= 0 .and. s1(3) == 0) then
    an = s1(4)
    f(2) = 1
    f(3) = 0
    f(4) = an
  else
    f(2) = 0
    f(3) = 0
    f(4) = 0
    an = 0.d0
  endif

  ca = mpdigout (an, mpndpw)

  do i = 1, mpndpw
    ix = ix + 1
    if (ix > nb + 50) then
      write (6, 2)
2     format ('MPEFORMAT: Insufficient space in B2 array.')
      call mpabrt ( 312)
    endif
    b2(ix) = ca(i:i)
  enddo

  call mpsub (s1, f, s0, mpnw1)
  call mpmuld (s0, d10w, s1, mpnw1)
enddo

!   Round the result.

if (ix >= nd + 1) then
  i1 = index (digits, b2(nd+1)) - 1
  if (i1 >= 5) then

!   Perform rounding, beginning at the last digit (position IX).  If the rounded
!   digit is 9, set to 0, then repeat at position one digit to left.  Continue
!   rounding if necessary until the decimal point is reached.

    do i = ix, ixp + 1, -1
      i2 = index (digits, b2(i)) - 1
      if (i2 <= 8) then
        b2(i) = digits(i2+2:i2+2)
        goto 180
      else
        b2(i) = '0'
      endif
    enddo

!   We have rounded up all digits to the right of the decimal point.  If the
!   digit to the left of the decimal point is a 9, then set that digit to 1
!   and increase the exponent by one; otherwise increase that digit by one.

    if (b2(ixp-1) == '9') then
      b2(ixp-1) = '1'
      nexp = nexp + 1
    else
      i1 = index (digits, b2(ixp-1)) - 1
      b2(ixp-1) = digits(i1+2:i1+2)
    endif
  endif
endif

180 continue

!   Done with mantissa.  Insert exponent.

ix = nd + 2
if (ia < 0) ix = ix + 1
b2(ix) = 'e'
if (nexp < 0) then
  ix = ix + 1
  b2(ix) = '-'
endif
ca = mpdigout (dble (abs (nexp)), 10)

do k = 1, 10
  if (ca(k:k) /= '0') goto 190
enddo

k = 10

190 continue

do i = k, 10
  ix = ix + 1
  b2(ix) = ca(i:i)
enddo

do i = ix + 1, nb
  b2(i) = ' '
enddo

!   Copy entire b2 array to B.

do i = 1, nb
  b(i) = b2(i)
enddo

return
end subroutine mpeformat

subroutine mpfformat (a, nb, nd, b, mpnw)

!   Converts the MPR number A into character form in the character(1) array B.
!   NB (input) is the length of the output string, and ND (input) is the
!   number of digits after the decimal point.  The format is analogous to
!   Fortran F format; the result is right-justified among the NB cells of B.
!   The condition NB >= ND + 10 must hold or an error message will result.
!   However, if it is found during execution that there is not sufficient space,
!   to hold all digits, the entire output field will be filled with asterisks.
!   NB cells of type character(1) must be available in B.

implicit none
integer (mpiknd), intent(in):: a(0:)
integer, intent(in):: mpnw, nb, nd
character(1), intent(out):: b(nb)
integer i, ia, ixp, i1, i2, i3, j, k, na, nb2, nb3, nexp
character(1) b2(nb+20)
character(16) ca
real (mpdknd) aa, t1

! End of declaration

if (mpnw < 4 .or. a(0) < abs (a(2)) + 4 .or. nb < nd + 10) then
  write (mpldb, 1)
1 format ('*** MPFFORMAT: uninitialized or inadequately sized arrays')
  call mpabrt ( 313)
endif

ia = sign (int (1, mpiknd), a(2))
if (a(2) == 0) ia = 0
na = min (int (abs (a(2))), mpnw)

if (ia == 0) then
  nb2 = nd + 13
else
  aa = a(4)
  if (na >= 2) aa = aa + dble (a(5)) / mpbdx
  nb2 = int (log10 (2.d0) * mpnbt * a(3) + log10 (aa)) + nd + 13
endif

nb3 = nb2 - 10
call mpeformat (a, nb2, nb3, b2, mpnw+1)

!   Trim off trailing blanks.

do i = nb2, 1, -1
  if (b2(i) /= ' ') goto 90
enddo

90 continue

nb2 = i

!   Look for the 'e' in B2.

do k = 1, nb2
  if (b2(k) == 'e') goto 100
enddo

write (6, 2)
2 format ('*** MPFFORMAT: Syntax error in output of mpeformat')
call mpabrt ( 314)

100 continue

!   Check the sign of the exponent.

k = k + 1
if (b2(k) == '-') then
  ixp = -1
  k = k + 1
else
  ixp = 1
endif
j = 0
ca = ' '

!   Copy the exponent into CA.

do i = k, nb2
  j = j + 1
  if (j <= 16) ca(j:j) = b2(i)
enddo

t1 = mpdigin (ca, j)

!   Check if there is enough space in the output array for all digits.

if (t1 + nd + 3 > nb) then
  do i = 1, nb
    b(i) = '*'
  enddo

  goto 210
endif
nexp = ixp * t1

!   Insert the sign of the number, if any.

i1 = 0
i2 = 0
if (b2(1) == '-') then
  i1 = i1 + 1
  b(i1) = '-'
  i2 = i2 + 1
endif

if (nexp == 0) then

!   Exponent is zero.  Copy first digit, period and ND more digits.

  do i = 1, nd + 2
    i1 = i1 + 1
    i2 = i2 + 1
    b(i1) = b2(i2)
  enddo

  goto 200
elseif (nexp > 0) then

!   Exponent is positive.  Copy first digit, skip the period, then copy
!   nexp digits.

  i1 = i1 + 1
  i2 = i2 + 1
  b(i1) = b2(i2)
  i2 = i2 + 1

  do i = 1, nexp
    i1 = i1 + 1
    i2 = i2 + 1
    b(i1) = b2(i2)
  enddo

!   Insert the period.

  i1 = i1 + 1
  b(i1) = '.'

!   Copy nd more digits.

  do i = 1, nd
    i1 = i1 + 1
    i2 = i2 + 1
    b(i1) = b2(i2)
  enddo

  goto 200
else

!   Exponent is negative.  Insert a zero, then a period, then nexp - 1
!   zeroes, then the first digit, then the remaining digits up to ND total
!   fractional digits.

  i1 = i1 + 1
  b(i1) = '0'
  i1 = i1 + 1
  b(i1) = '.'
  i3 = min (- nexp - 1, nd - 1)

  do i = 1, i3
    i1 = i1 + 1
    b(i1) = '0'
  enddo

  if (- nexp - 1 < nd) then
    i1 = i1 + 1
    i2 = i2 + 1
    b(i1) = b2(i2)
    i2 = i2 + 1
  endif

  do i = i3 + 2, nd
    i1 = i1 + 1
    i2 = i2 + 1
    b(i1) = b2(i2)
  enddo
endif

200 continue

!   Right-justify in field.

k = nb - i1

do i = 1, i1
  b(nb-i+1) = b(nb-i-k+1)
enddo

do i = 1, k
  b(i) = ' '
enddo

210 continue

return
end subroutine mpfformat

subroutine mpinp (iu, a, mpnw)

!   This routine reads the MPR number A from logical unit IU.  The digits of A
!   may span more than one line, provided that a "\" appears at the end of
!   a line to be continued (any characters after the "\" on the same line
!   are ignored).  Individual input lines may not exceed 2048 characters in
!   length, although this limit can be changed in the system parameters
!   (parameter mpnstr) in module MPFUNA.  Embedded blanks are allowed anywhere.
!   An exponent with "e" or "d" may optionally follow the numeric value.

!   A scratch array below (CHR1) holds character data for input to mpctomp.
!   It is dimensioned MPNW * (MPNDPW + 1) + 1000 (see below).  If more nonblank
!   input characters than this are input, they are ignored.

implicit none
integer, intent(in):: iu, mpnw
integer (mpiknd), intent(out):: a(0:)
integer i, i1, lnc1, lncx, ln1
character(mpnstr) line1
character(18), parameter:: validc = ' 0123456789+-.dDeE'
character(1) chr1(mpnw*int(mpdpw+1)+1000)

! End of declaration

if (mpnw < 4 .or. a(0) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPINP: uninitialized or inadequately sized arrays')
  call mpabrt ( 315)
endif

lnc1 = 0
lncx = mpnw * int (mpdpw + 1) + 1000

100 continue

read (iu, '(a)', end = 200) line1

!   Find the last nonblank character.

do i = mpnstr, 1, -1
  if (line1(i:i) /= ' ') goto 110
enddo

!   Input line is blank -- ignore.

goto 100

110 continue

ln1 = i

!   Scan input line, looking for valid characters.

do i = 1, ln1
  if (line1(i:i) == '\') goto 100
  i1 = index (validc, line1(i:i))
  if (i1 == 0 .and. line1(i:i) /= ' ') then
      write (6, 2) line1(i:i)
2     format ('*** MPINP: Invalid input character = ',a)
      call mpabrt ( 316)
  elseif (line1(i:i) /= ' ') then
    if (lnc1 < lncx) then
      lnc1 = lnc1 + 1
      chr1(lnc1) = line1(i:i)
    endif
  endif
enddo

call mpctomp (chr1, lnc1, a, mpnw)
goto 300

200  continue

write (mpldb, 4)
4 format ('*** MPINP: End-of-file encountered.')
call mpabrt ( 317)

300 return
end subroutine mpinp

subroutine mpout (iu, ln, nd, a, mpnw)

!   This routine writes MPR number A to logical unit IU in E(LN,ND) format.
!   This is output on MPOUTL characters per line.  The value of MPOUTL is set
!   in the system parameters at the start of module MPFUNA.

implicit none
integer, intent(in):: iu, ln, nd, mpnw
integer (mpiknd), intent(in):: a(0:)
integer i, ln1
character(1) chr1(ln)
character(32) cform1, cform2

! End of declaration

call mpeformat (a, ln, nd, chr1, mpnw)

write (cform1, 1) mpoutl
1 format ('(',i8,'a1)')
write (cform2, 2) mpoutl
2 format ('(',i8,'a1,"\")')

if (ln <= mpoutl) then
  write (iu, fmt = cform1) (chr1(i), i = 1, ln)
elseif (mod (ln, mpoutl) == 0) then
  ln1 = mpoutl * (ln / mpoutl) - mpoutl
  write (iu, fmt = cform2) (chr1(i), i = 1, ln1)
  write (iu, fmt = cform1) (chr1(i), i = ln1 + 1, ln)
else
  ln1 = mpoutl * (ln / mpoutl)
  write (iu, fmt = cform2) (chr1(i), i = 1, ln1)
  write (iu, fmt = cform1) (chr1(i), i = ln1 + 1, ln)
endif

return
end subroutine mpout

subroutine mprealch (aa, b, mpnw)

!   This converts the character string AA of arbitrary length to MPR.

implicit none
character (*), intent(in):: aa
integer (mpiknd), intent(out):: b(0:)
integer, intent(in):: mpnw
character (1) chr1(len(aa))
integer i, ln1

ln1 = len (aa)

do i = 1, ln1
  chr1(i) = aa(i:i)
enddo

call mpctomp (chr1, ln1, b, mpnw)
return
end subroutine mprealch


end module mpfunc

