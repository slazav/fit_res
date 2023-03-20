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
   coordinate:
   (X + iY) = (A + iB) + (C + iD)/(w0^2 - w^2 - iw*dw)
   speed:
   (X + iY) = (A + iB) + i*w*(C + iD)/(w0^2 - w^2 - iw*dw)

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
  " --coord (1|0)      -- do coordinate or speed fitting, default 1\n"
  " --pars (6|8)       -- number of parameters, default 8\n"
  " --show_zeros (1|0) -- write trailing zeros for unused parameters, default 0\n"
  " --fmt_out (1|0)    -- write <name>=<value> lines instead of table, default 0\n"
  ;
}

int
main (int argc, char *argv[]) {

  // default parameters
  bool do_fit = true;
  bool overload_detection = true;
  bool coord  = true;
  size_t p = 8;
  bool show_zeros = false;
  bool fmt_out = 0;


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
    if (strcasecmp(argv[i], "--coord") == 0)
      coord = atoi(argv[i+1]);
    else
    if (strcasecmp(argv[i], "--pars") == 0)
      p = atoi(argv[i+1]);
    else
    if (strcasecmp(argv[i], "--show_zeros") == 0)
      show_zeros = atoi(argv[i+1]);
    else
    if (strcasecmp(argv[i], "--fmt_out") == 0)
      fmt_out = atoi(argv[i+1]);
    else {
      print_help(); return 1;
    }
  }

  fit_func_t fit_func;
  if      (p==6 && coord==1) fit_func = OSCX_COFFS;
  else if (p==8 && coord==1) fit_func = OSCX_LOFFS;
  else if (p==6 && coord==0) fit_func = OSCV_COFFS;
  else if (p==8 && coord==0) fit_func = OSCV_LOFFS;
  else if (p==10 && coord==1) fit_func = DOSCX_COFFS;
  else if (p==10 && coord==0) fit_func = DOSCV_COFFS;
  else {
    print_help(); return 1;
  }

  std::vector<double> freq, real, imag, time;
  std::vector<double> pars(MAXPARS), pars_e(MAXPARS);

  // Read data (t,f,x,y), find max/min values
  double maxx=-INFINITY, maxy=-INFINITY, maxf=-INFINITY;
  double minx=INFINITY, miny=INFINITY, minf=INFINITY;
  while (!std::cin.eof()){
    std::string l;
    getline(std::cin, l);

    std::istringstream ss(l);
    double t,f,x,y;
    ss >> t >> f >> x >> y;
    if (ss.fail()) continue;
    time.push_back(t);
    freq.push_back(f);
    real.push_back(x);
    imag.push_back(y);

    // find max/min values
    if (x>maxx) maxx=x;
    if (y>maxy) maxy=y;
    if (f>maxf) maxf=f;
    if (x<minx) minx=x;
    if (y<miny) miny=y;
    if (f<minf) minf=f;
  }
  // for overload detection
  double maxax=std::max(fabs(maxx),fabs(minx));
  double maxay=std::max(fabs(maxx),fabs(miny));

  // too few data points
  if (freq.size()<p) return 0;

  // shift/scale data
  double x0 = (maxx+minx)/2;
  double y0 = (maxy+miny)/2;
  double sa = std::min(maxx-minx, maxy-miny);
  double sf = (maxf+minf)/2;
  for (size_t i=0; i<freq.size(); i++){
    real[i] = (real[i]-x0)/sa;
    imag[i] = (imag[i]-y0)/sa;
    freq[i] = freq[i]/sf;
  }

  // initial guess:
  fit_res_init(freq.size(), p,
     freq.data(), real.data(), imag.data(),
     pars.data(), fit_func);

  // avoid zero values in init.cond
  if (fabs(pars[0]) < 1e-6) pars[0] = 1e-6;
  if (fabs(pars[1]) < 1e-6) pars[1] = 1e-6;
  if (fabs(pars[6]) < 1e-6) pars[6] = 1e-6;
  if (fabs(pars[7]) < 1e-6) pars[7] = 1e-6;

  // fit
  double func_e = 0;
  if (do_fit) {
    func_e = fit_res(freq.size(), p,
       freq.data(), real.data(), imag.data(),
       pars.data(), pars_e.data(), fit_func);


    // overload detection (remove largest values and compare result)
    if (overload_detection) {
      std::vector<double> freq1, real1, imag1;
      std::vector<double> pars1(pars), pars_e1(MAXPARS);
      for (int i=0; i<freq.size(); i++){
        if (fabs(real[i]*sa+x0) > maxax*0.95 ||
            fabs(imag[i]*sa+y0) > maxay*0.95) continue;
        freq1.push_back(freq[i]);
        real1.push_back(real[i]);
        imag1.push_back(imag[i]);
      }
      if (freq1.size() >= p) {
        double func_e1 = fit_res(freq1.size(), p,
           freq1.data(), real1.data(), imag1.data(),
           pars1.data(), pars_e1.data(), fit_func);
        if (func_e1 < func_e) {
          pars.swap(pars1);
          pars_e.swap(pars_e1);
          func_e = func_e1;
        }
      }
    }
  }

  // shift/scale back
  func_e *= sa;
  pars[0] = (pars[0]*sa)+x0;  pars_e[0] *= sa;
  pars[1] = (pars[1]*sa)+y0;  pars_e[1] *= sa;
  if (coord) {
    pars[2] *= sa*sf*sf; pars_e[2] *= sa*sf*sf;
    pars[3] *= sa*sf*sf; pars_e[3] *= sa*sf*sf;
  }
  else {
    pars[2] *= sa*sf; pars_e[2] *= sa*sf;
    pars[3] *= sa*sf; pars_e[3] *= sa*sf;
  }
  pars[4] *= sf; pars_e[4] *= sf;
  pars[5] *= sf; pars_e[5] *= sf;

  if (p==8){
    pars[6] *= sa/sf; pars_e[6] *= sa/sf;
    pars[7] *= sa/sf; pars_e[7] *= sa/sf;
  }
  if (p==10){
    if (coord) {
      pars[6] *= sa*sf*sf; pars_e[6] *= sa*sf*sf;
      pars[7] *= sa*sf*sf; pars_e[7] *= sa*sf*sf;
    }
    else {
      pars[6] *= sa*sf; pars_e[6] *= sa*sf;
      pars[7] *= sa*sf; pars_e[7] *= sa*sf;
    }
    pars[8] *= sf; pars_e[8] *= sf;
    pars[9] *= sf; pars_e[9] *= sf;
  }

  double t = (*time.begin() + *time.rbegin())/2;

  if (fmt_out==0) {
    std::cout << std::setprecision(14)
              << std::fixed << " " << t
              << std::scientific
              << " " << func_e;
    for (size_t i = 0; i<p; i++) {
      std::cout << " " << pars[i]
                << " " << pars_e[i];
    }
    if (show_zeros && p==6) {
      std::cout << " 0 0 0 0";
    }
  }

  if (fmt_out==1) {
    std::cout << std::setprecision(14)
              << std::fixed << "t0=" << t << "\n"
              << std::scientific
              << "err=" << func_e << "\n";

    std::cout << "A="  << pars[0] << "\nA_err=" << pars_e[0] << "\n";
    std::cout << "B="  << pars[1] << "\nB_err=" << pars_e[1] << "\n";
    std::cout << "C="  << pars[2] << "\nC_err=" << pars_e[2] << "\n";
    std::cout << "D="  << pars[3] << "\nC_err=" << pars_e[3] << "\n";
    std::cout << "f0=" << pars[4] << "\nf0_err=" << pars_e[4] << "\n";
    std::cout << "df=" << pars[5] << "\ndf_err=" << pars_e[5] << "\n";
    if (p==8){
      std::cout << "E="  << pars[6] << "\nE_err=" << pars_e[6] << "\n";
      std::cout << "F="  << pars[7] << "\nF_err=" << pars_e[7] << "\n";
    }
    if (p==10){
      std::cout << "C2="  << pars[6] << "\nC2_err=" << pars_e[6] << "\n";
      std::cout << "D2="  << pars[7] << "\nC2_err=" << pars_e[7] << "\n";
      std::cout << "f02=" << pars[8] << "\nf02_err=" << pars_e[8] << "\n";
      std::cout << "df2=" << pars[9] << "\ndf2_err=" << pars_e[9] << "\n";
    }
    std::cout << "fit_func=" << (int)fit_func << "\n";
  }


  std::cout << "\n";
  return 0;
}
