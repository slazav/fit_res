all: mcta_1.plot mcta_2.plot mcta_3.plot mcta_ov.plot mctf_1.plot

clean:
	rm -f *.plot

# ports:
#  - 8084: f4_2019
#  - 8095: f4_2021

mcta_1.dat:
	wget "http://slazav.xyz:8084/get_range?name=wire_mc_cont&t1=1577276001&t2=1577276288" -O $@

mcta_2.dat:
	wget "http://slazav.xyz:8084/get_range?name=wire_mc_cont&t1=1577300767&t2=1577301054" -O $@

mcta_3.dat:
	wget "http://slazav.xyz:8095/get_range?name=mcta_sweeps&t1=1620034365&t2=1620034479" -O $@

mcta_ov.dat:
	wget "http://slazav.xyz:8084/get_range?name=wire_mc_cont&t1=1577309281&t2=1577309575" -O $@

mctf_1.dat:
	wget "http://slazav.xyz:8095/get_range?name=mctf_sweeps&t1=1617873768&t2=1617873812" -O $@


%.plot: %.dat make_plot
	./make_plot $< $@