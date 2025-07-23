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
! source https://github.com/scipy/scipy/blob/main/scipy/special/_spence.pxd
!
!     """
!     A series centered at z = 0; see
! 
!     http://functions.wolfram.com/10.07.06.0005.02
! 
!     """
! # Author: Josh Wilson
! #
! # Released under the same license as Scipy.

! 7/21/2025 : changed the number of terms in the series to get more precision than given by 64 bit floats

function cspence_series0(z)
  use mpmodule
  use param_mod, only:kv_nwds
        implicit none
type(mp_real) :: PISQ_6 
type(mp_real) :: TOL 
integer :: n, lub
type(mp_complex) :: zfac
type(mp_complex) :: sum1 
type(mp_complex) :: sum2 
type(mp_complex) :: term1, term2
type(mp_complex) :: z
type(mp_complex) :: cspence_series0
type(mp_complex) :: zlog1
type(mp_real) :: mppic
external zlog1

mppic=mppi(kv_nwds)
PISQ_6 = mppic**2 / 6.0
TOL = mpreald(10.0,kv_nwds)**(-100) ! 100 ~ number of digits of working precision of fort_kv_ints.  100 should be updated if the parameter kv_nwds gets changed

zfac = mpcmplx(cmplx(1.0,0.0),kv_nwds)
sum1 = mpcmplx(cmplx(0.0,0.0),kv_nwds)
sum2 = mpcmplx(cmplx(0.0,0.0),kv_nwds)

      ! expect convergence to mpipl digits before this loop upper bound is reached
    lub = int(1.1 * 100.0/log10(2.0) ) ! 100 ~ number of digits of working precision of fort_kv_ints.  100 should be updated if the parameter kv_nwds gets changed
    if (z == mpreal(0.0,kv_nwds)) then
        cspence_series0 = PISQ_6
    else
        do n=1,lub
!         do n=1,350
            zfac = z * zfac
            term1 = zfac/mpreal((n*1.0)**2,kv_nwds)
            sum1 = term1 + sum1
            term2 = zfac/mpreal(n*1.0,kv_nwds)
            sum2 = term2 + sum2
            if ((abs(term1).le.(TOL*abs(sum1))).and.(abs(term2).le.(TOL*abs(sum2)))) then
                exit
            endif
        enddo
        if(n.eq.(lub+1)) then
          write(*,*) 'spence series 0 end reached'
        endif
        cspence_series0 = PISQ_6 - sum1 + zlog1(z)*sum2
    endif

end function cspence_series0
