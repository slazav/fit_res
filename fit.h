#ifndef FIT_RES_H
#define FIT_RES_H

#define MAXPARS 8

enum fit_func_t {
  // Coordinate response of lineal oscillator (Lorentzian function) with constant offset:
  // (X + iY) = (A + iB) + (C + iD)/(w0^2 - w^2 - iw*dw) [ + (E + iF)*(w-w0) ]
  // X(w) = A + (C*(w0^2-w^2) + D*w*dw) / ((w0^2-w^2)^2 + (w*dw)^2) [ + E*(w-w0) ]
  // Y(w) = B + (D*(w0^2-w^2) - C*w*dw) / ((w0^2-w^2)^2 + (w*dw)^2) [ + F*(w-w0) ]
  // 6 parameters
  OSCX_COFFS,

  // Same with linear offset, additional term
  // (X + iY) = ... + (E + iF)*(w-w0)
  // 8 parameters
  OSCX_LOFFS,

  // Velocity response of lineal oscillator with constant offset:
  // (X + iY) = (A + iB) + i*w*(C + iD)/(w0^2 - w^2 - iw*dw) [ + (E + iF)*(w-w0) ]
  // X(w) = A - w*(D*(w0^2-w^2) - C*w*dw) / ((w0^2-w^2)^2 + (w*dw)^2) [ + E*(w-w0) ]
  // Y(w) = B + w*(C*(w0^2-w^2) + D*w*dw) / ((w0^2-w^2)^2 + (w*dw)^2) [ + F*(w-w0) ]
  // 6 parameters
  OSCV_COFFS,

  // Same with linear offset, additional term
  // (X + iY) = ... + (E + iF)*(w-w0)
  // 8 parameters
  OSCV_LOFFS,
};

/*
Find initial conditions by some trivial assumptions.
Arguments:
  n   - number of points
  p   - number of parameters (6 or 8)
  freq - frequency data [0..n-1]
  real - X (real part) of the data [0..n-1]
  imag - Y (imag part) of the data [0..n-1]
  pars - array of size 8, fit parameters to be returned:
         A,B,C,D,w,dw,E,F
*/
void fit_res_init (const size_t n, const size_t p,
         double * freq, double * real, double * imag,
         double pars[MAXPARS], fit_func_t fit_func);


/*
Fit resonance with Lorentzian curve
Arguments:
  n   - number of points
  p   - number of parameters (6 or 8)
  freq - frequency data [0..n-1]
  real - X (real part) of the data [0..n-1]
  imag - Y (imag part) of the data [0..n-1]
  pars - array of size 8, fit parameters:
         A,B,C,D,w,dw,E,F
         On input initial values (e.g from fit_res_init) should be provided,
         On output parameters are changed to new values.
  pars_e -- On output: parameter errors
  coord -- coord/speed function (1|0)
Return value:
  mean square difference of the fitted function
*/

double fit_res (const size_t n, const size_t p,
                double * freq, double * real, double * imag,
                double pars[MAXPARS], double pars_e[MAXPARS],
                fit_func_t fit_func);

#endif
