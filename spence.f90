! 
! Copyright © 2001, 2002 Enthought, Inc.
! All rights reserved.
! 
! Copyright © 2003-2019 SciPy Developers.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 
!     Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 
!     Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 
!     Neither the name of Enthought nor the names of the SciPy Developers may be used to endorse or promote products derived from this software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! source https://github.com/scipy/scipy/blob/main/scipy/special/_spence.pxd
!
! # Implement Spence's function, a.k.a. the dilogarithm, for complex
! # arguments. Note that our definition differs from that in the sources
! # by the mapping z -> 1 - z.
! #
! # Sources
! # [1] Zagier, "The Dilogarithm Function"
! # [2] functions.wolfram.com
! # [3] Ginsberg, Zaborowski, "The Dilogarithm Function of a Real Argument"
! #
! # Author: Josh Wilson
! #
! # Released under the same license as Scipy.
!
!           ⌠z  log(t)
! Spence(z)=⎮   ────── dt
!           ⌡1   1-t
!


function spence(zin)
  use mpmodule
  use param_mod, only:kv_nwds
  implicit none
type(mp_complex) :: zin
type(mp_complex) :: z
type(mp_complex) :: spence
type(mp_complex) :: cspence_series0
type(mp_complex) :: cspence_series1
type(mp_complex) :: zlog1
type(mp_real) :: mppic, PISQ_6
external cspence_series0
external cspence_series1
external zlog1
! type(mp_complex) :: logw

mppic=mppi(kv_nwds)
PISQ_6 = mppic**2 / 6.0


!     """
!     Compute Spence's function for complex arguments. The strategy is:
!     - If z is close to 0, use a series centered at 0.
!     - If z is far away from 1, use the reflection formula
! 
!     spence(z) = -spence(z/(z - 1)) - pi**2/6 - ln(z - 1)**2/2
! 
!     to move close to 1. See [1].
!     - If z is close to 1, use a series centered at 1.
! 
!     """

   z = zin 

    if (abs(z).lt.mpreal(0.5,kv_nwds)) then
!         # This step isn't necessary, but this series converges faster.
        spence = cspence_series0(z)
    else if (abs(mpreal(1.0,kv_nwds) - z).gt.mpreal(1.0,kv_nwds)) then
!         spence = -cspence_series1(z/(z - 1.0)) - PISQ_6 - 0.5*logw(z - 1.0)**2
        spence = -cspence_series1(z/(z - mpreal(1.0,kv_nwds))) - PISQ_6 - mpreal(0.5,kv_nwds)*zlog1(z - mpreal(1.0,kv_nwds))**2
    else
        spence = cspence_series1(z)
    endif

end function spence
