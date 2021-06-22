#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include "math.h"

#include "fit.h"

/*
 Program reads resonance data (time, freq, x, y) from stdin,
 fits data with a Lorentzian function:
   (X + iY) = (A + iB) + (C + iD)/(w0^2 - w^2 - iw*dw)

 and prints a line with 14 values to stdout:
   - time -- center of the time range
   - function error -- mean square difference between data and fit
   - A, A_error, B, B_error -- shift in x and y component
   - C, C_error, D, D_error -- driving force (not amplitude!
   - w0, w0_error == resonance frequency (in Hz or rad/s depending on input data)
   - dw, dw_error -- width at 1/2 height of amplitude curve,
     or (approximately) distance between dispersion minimum
     and maximum. If frequency data is in rad/s then df=2/tau,
     if data in Hz then df = pi/tau

The program can be used as a filter in
graphene_filter script.

*/


void
print_help() {
  std::cerr <<
  "Usage: fit_res [options] < file\n"
  "Options\n"
  " --do_fit (1|0)     -- do fitting or stop after initial guess for parameters, default 1\n"
  " --overload (1|0)   -- use overload detection, default 1\n"
  " --pars (6|8)       -- number of parameters, default 8\n"
  " --show_zeros (1|0) -- write trailing zeros for unused parameters, default 0\n"
  ;
}

int
main (int argc, char *argv[]) {

  // default parameters
  bool do_fit = true;
  bool overload_detection = true;
  bool show_zeros = false;
  size_t p = 8;


  // parse command-line options
  if (argc%2 != 1) {
    print_help(); return 1;
  }
  for (int i=1; i<argc-1; i+=2) {
    if (strcasecmp(argv[i], "--do_fit") == 0)
      do_fit = atoi(argv[i+1]);
    else
    if (strcasecmp(argv[i], "--overload") == 0)
      overload_detection = atoi(argv[i+1]);
    else
    if (strcasecmp(argv[i], "--pars") == 0)
      p = atoi(argv[i+1]);
    else
    if (strcasecmp(argv[i], "--show_zeros") == 0)
      show_zeros = atoi(argv[i+1]);
    else {
      print_help(); return 1;
    }
  }

  if (p != 6 && p != 8) {
    print_help(); return 1;
  }


  std::vector<double> freq, real, imag, time;
  std::vector<double> pars(MAXPARS), pars_e(MAXPARS);

  double maxx=0, maxy=0;
  while (!std::cin.eof()){
    std::string l;
    getline(std::cin, l);

    std::istringstream ss(l);
    double t,f,x,y;
    ss >> t >> f >> x >> y;
    time.push_back(t);
    freq.push_back(f);
    real.push_back(x);
    imag.push_back(y);

    // find maximum - for overload detection
    if (fabs(x)>maxx) maxx=fabs(x);
    if (fabs(y)>maxy) maxy=fabs(y);
  }

  if (freq.size()<p) return 0;

  // initial guess:
  fit_res_init(freq.size(), p,
     freq.data(), real.data(), imag.data(),
     pars.data());

  // fit
  double func_e = 0;
  if (do_fit) {
    func_e = fit_res(freq.size(), p,
       freq.data(), real.data(), imag.data(),
       pars.data(), pars_e.data());

    // overload detection (remove largest values and compare result)
    if (overload_detection) {
      std::vector<double> freq1, real1, imag1;
      std::vector<double> pars1(pars), pars_e1(MAXPARS);
      for (int i=0; i<freq.size(); i++){
        if (fabs(real[i]) > maxx*0.95 || fabs(imag[i]) > maxy*0.95) continue;
        freq1.push_back(freq[i]);
        real1.push_back(real[i]);
        imag1.push_back(imag[i]);
      }
      if (freq1.size() >= p) {
        double func_e1 = fit_res(freq1.size(), p,
           freq1.data(), real1.data(), imag1.data(),
           pars1.data(), pars_e1.data());
        if (func_e1 < func_e) {
          pars.swap(pars1);
          pars_e.swap(pars_e1);
          func_e = func_e1;
        }
      }
    }
  }

  double t = (*time.begin() + *time.rbegin())/2;

  std::cout << std::setprecision(14)
            << std::fixed << " " << t
            << std::scientific
            << " " << func_e;
  for (size_t i = 0; i<p; i++) {
    std::cout << " " << pars[i]
              << " " << pars_e[i];
  }
  if (show_zeros) {
    for (size_t i = p; i<MAXPARS; i++) std::cout << " 0 0";
  }

  std::cout << "\n";
  return 0;
}
