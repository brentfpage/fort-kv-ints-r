!*****************************************************************************

!  MPFUN20-Fort: A thread-safe arbitrary precision package with special functions
!  Precision level declaration module (module MPFUNF)

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
!    minor changes to existing Fortran-90 programs. All basic arithmetic
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

!  DESCRIPTION OF THIS MODULE (MPFUNF):
!    This module defines the default standard precision level in digits (mpipl)
!    and default medium precision level in digits (mpiplm), and the equivalent
!    precision levels in words (mpwds and mpwdsm), which are calculated below as:
!       mpwds = int (mpipl / mpdpw + 2)
!       mpwdsm = int (mpiplm / mpdpw + 2)
!    (mpdpw is the approx. number of digits per word, set in module MPFUNA).
!    These precision levels are the maximum working precision levels for all
!    operations that use module MPFUNG and MPFUNH.

module mpfunf
use mpfuna

implicit none

!  *** Set the default standard and medium precision levels (in digits) here.

integer, public, parameter:: mpipl = 109, mpiplm = 109

!----------------------------------------------------------------------------

!  *** Do not change the following code (in normal usage).

integer, public, parameter:: mpwds = int (mpipl / mpdpw + 2.d0), mpwds6 = mpwds + 6, &
  mpwdsm = int (mpiplm / mpdpw + 2.d0), mpwdsm6 = mpwdsm + 6

end module mpfunf
