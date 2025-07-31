fort-kv-ints is an add-on to LEOPARD (www.github.com/pastfalk/LEOPARD) that partially computes the nonlinear growth rate of the weak turbulence theory for parallel-propagating electromagnetic waves in magnetized plasma.
The complete nonlinear growth rate includes contributions from three-wave interactions and induced scattering, and fort-kv-ints computes part of the induced scattering contribution.
Some details of the computation can be found in github.com/brentfpage/fort-kv-ints-r/blob/main/preprint.pdf .

File descriptions:
    The following files in /fort_kv_ints_r are revised versions of files from LEOPARD
        - all mp*f90 files (mpfuna.f90 mpfunbq.f90 mpfunc.f90 mpfund.f90 mpfune.f90 mpfunf.f90 mpfungq2.f90 mpmodule.f90)
        - main.f90
        - muller.f90
        - param_mod.f90
        - int_para_mpfun.f90
        - rh_disp_val.f90
        - input.dat
        - the directory /distribution/
        - makefile

        Descriptions of modifications:
            -- The mp*f90 files in /fort_kv_ints_r come from a more recent version (MPFUN20-Fort (v32) by David H. Bailey) of the arbitrary precision arithmetic package used in LEOPARD. 
            -- The new main.f90 retains much of the same content as in /LEOPARD but now also drives the computation of the induced scattering growth rate.
            -- Relative to the LEOPARD version, muller.f90 in /fort_kv_ints_r features the replacement 'disp_det'->'rh_disp_val'.
            -- fort-kv-ints parameters have been appended to the new param_mod.f90
            -- In the /fort_kv_ints_r version of int_para_mpfun.f90, the arguments to the mp type constructors have been changed.
            -- rh_disp_val.f90 in /fort_kv_ints_r uses code from the LEOPARD files disp_det.f90 and integrator.f90 to compute the dialectric constant for right-handed parallel-propagating waves
            -- /distribution/distribution1.dat in /fort_kv_ints_r is an electron vdf derived from measurements by the ARTEMIS electrostatic analyzer during a whistler wave event.  For details, see github.com/brentfpage/fort-kv-ints-r/blob/main/preprint.pdf .
            -- compilation blocks for fort-kv-ints files have been appended to the new makefile

    New files in /fort_kv_ints_r:
        - kv_ints_mod.f90: computes the induced scattering integrals
        - pppack_mod_mp.f90: includes select functions from PPPACK by Carl de Boor, specifically FORTRAN 90 versions made available by John Burkardt, as well as driver functions to create quintic 2d splines and quadratic 1d splines.
        - spence.f90 , cspence_series0.f90, cspence_series1.f90, zlog1.f90: functions for computing Spence's function, which is nearly equivalent to the dilogarithm.  These were originally written in Python/Cython/C++ by Josh Wilson for SciPy
        - disp_deriv.f90: evaluates the derivative of the real part of (dialectric_constant - (c*k/omega)^2 ) with respect to omega using LEOPARD computations
        - fact.f90: factorial function


Building the program:
    First, download LEOPARD and copy the following files into this directory (/fort_kv_ints_r).
        - acc_Kvpa.f90
        - cerror.f90
        - cont_frac.f90
        - get_splinecoeff.f90
        - int_para.f90
        - polyfit.f90
        - read_data.f90
        - read_distr.f90
        - spline_interpol.f90
        - Z_func.f90
    On the command line, navigate to this directory (/fort_kv_ints_r) , type "make" (without quotes), and hit return.  The program should then go through a compilation process.  If the build is successful, the command "./dsolve" will start the program.  The default input parameters and particle velocity distribution function (vdf) are those that were used in github.com/brentfpage/fort-kv-ints-r/blob/main/preprint.pdf.  If 32 threads are used, this default program run takes about 7 hours.


Scope:
    The weak turbulence computations in kv_ints_mod.f90 have only been tested for a one-species plasma, but generalizing them for multi-species plasma would not be complicated.  Also, the program currently runs for right-handed waves, but it would be straightforward to adapt rh_disp_val.f90, which computes the dialectric constant, and kv_ints_mod.f90 to make it instead run for left-handed waves.

fort-kv-ints output:
    The program output is the file 'omega2.dat'.  The columns in this file are
        k    k'   k_pow    gam2_is(k, k', k_pow)

    The k grid is specified by the user in input.dat, as described in the LEOPARD Readme.  The k' grid is staggered relative to the k grid.

    For a wave kinetic equation simulation, the magnetic field spectral density should be defined on the k grid and splined.  At present, a degree 2 spline with not-a-knot boundary conditions is a sensible choice.  In python, scipy.interpolate.make_interp_spline(x, y, 2) produces such a spline whose knots are the same as the k' grid values in omega2.dat.

    If the magnetic field spectral density B between grid points k'_1 and k'_2 has the spline representation
        B = Bsp0 + Bsp1*k' + Bsp2*(k')^2
    then the contribution to the nonlinear growth rate at some grid point k=k_eval from the region of wavenumber space k'_1 -> k'_2 is
        Bsp0*gam2_is(k_eval, k'_1, 0) + Bsp1*gam2_is(k_eval, k'_1, 1) + Bsp2*gam2_is(k_eval, k'_1, 2)
    where it has been assumed that the knots of the spline of B are the same as the k' grid values of gam2_is.


Re-use:
    Redistributing modified versions of kv_ints_mod.f90 is permitted, but please include a description of the modifications.


Contact:
    Questions about fort-kv-ints can be directed to brentfpage@gmail.com .
