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
!
! source https://github.com/scipy/scipy/blob/main/scipy/special/_complexstuff.pxd
!
!     """
!     Compute log, paying special attention to accuracy around 1. We
!     implement this ourselves because some systems (most notably the
!     Travis CI machines) are weak in this regime.
! 
!     """
!
! 7/21/2025 : changed the number of terms in the series to get more precision than given by 64 bit floats
function zlog1(z)
  use mpmodule
  use kv_ints_mod, only: logw
  use param_mod, only: kv_nwds
        implicit none
        integer :: n
        integer :: lub
        type(mp_complex) :: coeff
        type(mp_complex) :: res
        type(mp_real) :: tol 
        type(mp_complex) :: z, zm1
        type(mp_complex) :: zlog1
        type(mp_complex) :: logw

        tol = mpreald(10.0,kv_nwds)**(-100) ! 100 ~ number of digits of working precision of fort_kv_ints.  100 should be updated if the parameter kv_nwds gets changed
        coeff = mpcmplx(cmplx(-1.0,0.0),kv_nwds)
        res = mpcmplx(cmplx(0.0,0.0),kv_nwds)

        lub = int(100.0*1.1) ! 100 ~ number of digits of working precision of fort_kv_ints.  100 should be updated if the parameter kv_nwds gets changed
    if (abs(z - mpreal(1.0,kv_nwds)).gt.mpreald(0.1,kv_nwds)) then
        zlog1 = logw(z)
    else
        zm1 = z - mpreal(1.0,kv_nwds)
        if (zm1.eq.mpreal(0.0,kv_nwds)) then
            zlog1 = mpcmplx((0.0,0.0),kv_nwds)
        else
          do n = 1,lub 
                coeff = -coeff * zm1 
                res = res + coeff/mpreal(n*1.0,kv_nwds)
                if (abs(coeff/res).lt.tol) then ! 7/21/2025, changed from res/coeff to coeff/res
                    exit
                endif
            enddo
            zlog1 = res
        endif
        if(n.eq.(lub+1)) then
          write(*,*) 'zlog1 lub reached'
        endif
    endif
end function zlog1
