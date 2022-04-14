all: mcta_1.plot mcta_2.plot mcta_3.plot mcta_ov.plot mctf_1.plot w2ta_1.plot mcta_d1.plot

clean:
	rm -f *.plot

# ports:
#  - 8084: f4_2019
#  - 8095: f4_2021


#### downloading data
mcta_1.dat:
	wget "http://slazav.xyz:8084/get_range?name=wire_mc_cont&t1=1577276001&t2=1577276288" -O $@

mcta_2.dat:
	wget "http://slazav.xyz:8084/get_range?name=wire_mc_cont&t1=1577300767&t2=1577301054" -O $@

mcta_3.dat:
	wget "http://slazav.xyz:8095/get_range?name=mcta_sweeps&t1=1620034365&t2=1620034479" -O $@

mcta_d1.dat:
	wget "http://slazav.xyz:8095/get_range?name=mcta_sweeps&t1=1649929050&t2=1649929275" -O $@

mcta_ov.dat:
	wget "http://slazav.xyz:8084/get_range?name=wire_mc_cont&t1=1577309281&t2=1577309575" -O $@

mctf_1.dat:
	wget "http://slazav.xyz:8095/get_range?name=mctf_sweeps&t1=1617873768&t2=1617873812" -O $@

w2ta_1.dat:
	wget "http://slazav.xyz:8095/get_range?name=w2ta_sweeps&t1=1617659809&t2=1617662075" -O $@

#### fit_res parameters
PARS=
mcta_1.plot: PARS=--pars 6 --overload 0 --coord 0
mcta_2.plot: PARS=--pars 6 --overload 0 --coord 0
mcta_3.plot: PARS=--pars 6 --overload 0
mcta_d1.plot: PARS=--pars 10

%.plot: %.dat make_plot Makefile
	./make_plot $< $@ $(PARS)