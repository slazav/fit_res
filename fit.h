#ifndef FIT_RES_H
#define FIT_RES_H

#define MAXPARS 8

/*
Model function:

"coordinate":
(X + iY) = (A + iB) + (C + iD)/(w0^2 - w^2 - iw*dw) [ + (E + iF)*(w-w0) ]
X(w) = A + (C*(w0^2-w^2) + D*w*dw) / ((w0^2-w^2)^2 + (w*dw)^2) [ + E*(w-w0) ]
Y(w) = B + (D*(w0^2-w^2) - C*w*dw) / ((w0^2-w^2)^2 + (w*dw)^2) [ + F*(w-w0) ]

"speed":
(X + iY) = (A + iB) + i*w*(C + iD)/(w0^2 - w^2 - iw*dw) [ + (E + iF)*(w-w0) ]
X(w) = A - w*(D*(w0^2-w^2) - C*w*dw) / ((w0^2-w^2)^2 + (w*dw)^2) [ + E*(w-w0) ]
Y(w) = B + w*(C*(w0^2-w^2) + D*w*dw) / ((w0^2-w^2)^2 + (w*dw)^2) [ + F*(w-w0) ]

*/

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
         double pars[MAXPARS], bool coord);


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
                double pars[MAXPARS], double pars_e[MAXPARS], bool coord);

#endif
