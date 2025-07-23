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
!
! source https://github.com/scipy/scipy/blob/main/scipy/special/_spence.pxd
! """
! A series centered at z = 1 which enjoys faster convergence than
! the Taylor series. See [3]. The number of terms used comes from
! bounding the absolute tolerance at the edge of the radius of
! convergence where the sum is O(1).
! 
! """
! # Author: Josh Wilson
! #
! # Released under the same license as Scipy.

! 7/21/2025 : changed the number of terms in the series to get more precision than given by 64 bit floats

function cspence_series1(z)
  use mpmodule
  use param_mod, only:kv_nwds
    implicit none
    integer(kind=8) :: lub
    integer(kind=8) :: n
    integer(kind=8) :: ten
    type(mp_real) :: n_mp
    type(mp_real) :: n1
    type(mp_real) :: n2
    type(mp_real) :: n3
    type(mp_complex) :: zfac 
    type(mp_complex) :: res 
    type(mp_complex) :: term, z, zz, mzm1
    type(mp_complex) :: cspence_series1
    type(mp_complex) :: zlog1
    type(mp_real) :: TOL 
    external zlog1

    TOL = mpreald(10.0,kv_nwds)**(-100) ! 100 ~ number of digits of working precision of fort_kv_ints.  100 should be updated if the parameter kv_nwds gets changed
    ten=10
    lub=ten**15

    zfac = mpcmplx(cmplx(1.0,0.0),kv_nwds)
    res = mpcmplx(cmplx(0.0,0.0),kv_nwds)

    if (z.eq.mpreal(1.0,kv_nwds)) then
        cspence_series1 = mpcmplx((0.0,0.0),kv_nwds)
    else
        mzm1 = mpreal(1.0,kv_nwds) - z
        zz = mzm1**2
        do n=1,lub
            zfac = mzm1 * zfac
            n_mp = mpreald(n*1.0_8,kv_nwds)
            n1 = n_mp**2
            n2 = (n_mp+1)**2
            n3 = (n_mp+2)**2
            term = ((zfac/n1)/n2)/n3
            res = term + res
            if(abs(term).le.(TOL*abs(res))) then
                exit
            endif
        enddo
        if(n.eq.(lub+1)) then
          write(*,*) 'spence_series1 lub reached'
        endif
        res = 4.0*zz * res
        res = res + 4.0*mzm1 + 5.75*zz + 3.0*(mpreal(1.0,kv_nwds) - zz)*zlog1(z)
        res = res/(mpreal(1.0,kv_nwds) + 4.0*mzm1 + zz)
    endif
    cspence_series1 = res
end function cspence_series1
