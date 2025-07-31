objects=param_mod.o main.o Z_func.o muller.o polyfit.o read_data.o read_distr.o spline_interpol.o cerror.o cont_frac.o int_para.o int_para_mpfun.o acc_Kvpa.o get_splinecoeff.o
kv_ints_obj = kv_ints_mod.o fact.o zlog1.o cspence_series0.o cspence_series1.o spence.o rh_disp_val.o disp_deriv.o pppack_mod_mp.o
mpfun_obj = mpfuna.o mpfunb.o mpfunc.o mpfund.o mpfune.o mpfunf.o mpfung2.o mpfunh2.o mpmask13.o mpmodule.o
mpfun_mod = mpfuna.mod  mpfunb.mod mpfunc.mod mpfund.mod mpfune.mod mpfunf.mod mpfung.mod mpfunh.mod mpmodule.mod
f90comp = gfortran
options =  -fdefault-real-8 -g -O3 -fopenmp -fbounds-check -ffpe-trap=invalid -ffpe-trap=overflow -fimplicit-none # -ffpe-trap=denormal -Wall -Wcompare-reals
options_mp = -O3 -g
options_mp_fast = -Ofast -g
mpfun_dir = mpfun20-fort-v32-var2

dsolve: $(objects) $(kv_ints_obj) $(mpfun_obj)
	$(f90comp) -o dsolve $(options) $(objects) $(kv_ints_obj) $(mpfun_obj)

param_mod.mod: param_mod.o param_mod.f90
	$(f90comp) -c $(options) param_mod.f90

param_mod.o: param_mod.f90
	$(f90comp) -c $(options) param_mod.f90

mpmodule.mod: mpmodule.o ./$(mpfun_dir)/mpmodule.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpmodule.f90

mpmodule.o: ./$(mpfun_dir)/mpmodule.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpmodule.f90

mpfuna.mod: mpfuna.o ./$(mpfun_dir)/mpfuna.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfuna.f90

mpfuna.o: ./$(mpfun_dir)/mpfuna.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfuna.f90

mpfunb.mod: mpmask13.o mpfunb.o ./$(mpfun_dir)/mpfunb.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfunb.f90

mpfunb.o: ./$(mpfun_dir)/mpfunb.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfunb.f90

mpfunc.mod: mpfunc.o ./$(mpfun_dir)/mpfunc.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfunc.f90

mpfunc.o: ./$(mpfun_dir)/mpfunc.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfunc.f90

mpfund.mod: mpfund.o ./$(mpfun_dir)/mpfund.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfund.f90

mpfund.o: ./$(mpfun_dir)/mpfund.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfund.f90

mpfune.mod: mpfune.o ./$(mpfun_dir)/mpfune.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfune.f90

mpfune.o: ./$(mpfun_dir)/mpfune.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfune.f90

mpfunf.mod: mpfunf.o ./$(mpfun_dir)/mpfunf.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfunf.f90

mpfunf.o: ./$(mpfun_dir)/mpfunf.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfunf.f90

mpfung.mod: mpfung2.o ./$(mpfun_dir)/mpfung2.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfung2.f90

mpfung2.o: ./$(mpfun_dir)/mpfung2.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfung2.f90

mpfunh.mod: mpfunh2.o ./$(mpfun_dir)/mpfunh2.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfunh2.f90

mpfunh2.o: ./$(mpfun_dir)/mpfunh2.f90
	$(f90comp) -c $(options_mp_fast) ./$(mpfun_dir)/mpfunh2.f90

mpmask13.o: ./$(mpfun_dir)/mpmask13.f90
	$(f90comp) -c $(options_mp) ./$(mpfun_dir)/mpmask13.f90

main.o: param_mod.mod kv_ints_mod.mod pppack_mod_mp.mod main.f90
	$(f90comp) -c $(options) main.f90

disp_det.o: param_mod.mod disp_det.f90
	$(f90comp) -c $(options) disp_det.f90

Z_func.o: param_mod.mod Z_func.f90
	$(f90comp) -c $(options)  Z_func.f90

cerror.o: param_mod.mod cerror.f90
	$(f90comp) -c $(options)  cerror.f90

cont_frac.o: cont_frac.f90
	$(f90comp) -c $(options)  cont_frac.f90

muller.o: param_mod.mod muller.f90
	$(f90comp) -c $(options)  muller.f90

polyfit.o: polyfit.f90
	$(f90comp) -c $(options)  polyfit.f90

read_data.o: param_mod.mod read_data.f90
	$(f90comp) -c $(options)  read_data.f90

read_distr.o: param_mod.mod read_distr.f90
	$(f90comp) -c $(options)  read_distr.f90

spline_interpol.o: spline_interpol.f90
	$(f90comp) -c $(options)  spline_interpol.f90

get_splinecoeff.o: param_mod.mod get_splinecoeff.f90
	$(f90comp) -c $(options)  get_splinecoeff.f90

acc_Kvpa.o: acc_Kvpa.f90
	$(f90comp) -c $(options)  acc_Kvpa.f90

int_para.o: param_mod.mod int_para.f90
	$(f90comp) -c $(options)  int_para.f90

int_para_mpfun.o: param_mod.mod $(mpfun_mod) int_para_mpfun.f90
	$(f90comp) -c $(options)  int_para_mpfun.f90

kv_ints_mod.mod: $(mpfun_obj) kv_ints_mod.o kv_ints_mod.f90
	$(f90comp) -c $(options) kv_ints_mod.f90

kv_ints_mod.o: $(mpfun_obj) spence.o kv_ints_mod.f90
	$(f90comp) -c $(options) kv_ints_mod.f90

fact.o: fact.f90
	$(f90comp) -c $(options) fact.f90

zlog1.o: $(mpfun_obj) param_mod.mod zlog1.f90
	$(f90comp) -c $(options) zlog1.f90

cspence_series0.o: $(mpfun_obj) param_mod.mod cspence_series0.f90 zlog1.f90
	$(f90comp) -c $(options) cspence_series0.f90

cspence_series1.o: $(mpfun_obj) param_mod.mod cspence_series1.f90 zlog1.f90
	$(f90comp) -c $(options) cspence_series1.f90

spence.o: $(mpfun_obj) param_mod.mod spence.f90 cspence_series0.f90 cspence_series1.f90 zlog1.f90
	$(f90comp) -c $(options) spence.f90

rh_disp_val.o: param_mod.mod rh_disp_val.f90
	$(f90comp) -c $(options) rh_disp_val.f90

disp_deriv.o: param_mod.mod disp_deriv.f90
	$(f90comp) -c $(options) disp_deriv.f90

pppack_mod_mp.mod: pppack_mod_mp.o param_mod.o pppack_mod_mp.f90
	$(f90comp) -c $(options) pppack_mod_mp.f90

pppack_mod_mp.o: param_mod.mod pppack_mod_mp.f90
	$(f90comp) -c $(options) pppack_mod_mp.f90

clean:
	rm dsolve 
	rm param_mod.mod
	rm kv_ints_mod.mod
	rm pppack_mod_mp.mod
	rm $(objects)
	rm $(kv_ints_obj)
	rm $(mpfun_mod)
	rm $(mpfun_obj)
