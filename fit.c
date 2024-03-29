#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "fit.h"
/********************************************************************/

// modified Gaussian example from
//https://www.gnu.org/software/gsl/doc/html/nls.html#c.gsl_multifit_nlinear_fdf

struct data {
  double *w;
  double *x;
  double *y;
  size_t n;
  fit_func_t fit_func;
};

int
func_f (const gsl_vector * x, void *params, gsl_vector * f) {
  struct data *d = (struct data *) params;
  double A = gsl_vector_get(x, 0);
  double B = gsl_vector_get(x, 1);
  double C = gsl_vector_get(x, 2);
  double D = gsl_vector_get(x, 3);
  double w0 = gsl_vector_get(x, 4);
  double dw = gsl_vector_get(x, 5);
  double E = (x->size == 8) ? gsl_vector_get(x, 6) : 0.0;
  double F = (x->size == 8) ? gsl_vector_get(x, 7) : 0.0;
  double C2 = (x->size == 10) ? gsl_vector_get(x, 6) : 0.0;
  double D2 = (x->size == 10) ? gsl_vector_get(x, 7) : 0.0;
  double w02 = (x->size == 10) ? gsl_vector_get(x, 8) : 0.0;
  double dw2 = (x->size == 10) ? gsl_vector_get(x, 9) : 0.0;

  size_t i;
  for (i = 0; i < d->n; ++i) {
    double wi = d->w[i];
    double Xi = d->x[i];
    double Yi = d->y[i];

    double wa = w0*w0 - wi*wi;
    double wb = wi*dw;
    double z = wa*wa + wb*wb;

    double wa2 = w02*w02 - wi*wi;
    double wb2 = wi*dw2;
    double z2 = wa2*wa2 + wb2*wb2;

    double X=0, Y=0;
    switch (d->fit_func) {
     case OSCX_COFFS:
        X = A + (C*wa + D*wb)/z;
        Y = B + (D*wa - C*wb)/z;
        break;
     case OSCX_LOFFS:
        X = A + (C*wa + D*wb)/z + E*(wi-w0);
        Y = B + (D*wa - C*wb)/z + F*(wi-w0);
        break;
     case OSCV_COFFS:
        X = A - wi*(D*wa - C*wb)/z;
        Y = B + wi*(C*wa + D*wb)/z;
        break;
     case OSCV_LOFFS:
        X = A - wi*(D*wa - C*wb)/z + E*(wi-w0);
        Y = B + wi*(C*wa + D*wb)/z + F*(wi-w0);
        break;
     case DOSCX_COFFS:
        X = A + (C*wa + D*wb)/z + (C2*wa2 + D2*wb2)/z2;
        Y = B + (D*wa - C*wb)/z + (D2*wa2 - C2*wb2)/z2;
        break;
     case DOSCV_COFFS:
        X = A - wi*(D*wa - C*wb)/z - wi*(D2*wa2 - C2*wb2)/z2;
        Y = B + wi*(C*wa + D*wb)/z + wi*(C2*wa2 + D2*wb2)/z2;
        break;
    }
    gsl_vector_set(f, 2*i,   Xi - X);
    gsl_vector_set(f, 2*i+1, Yi - Y);
  }

  return GSL_SUCCESS;
}


// function derivatives
int
func_df (const gsl_vector * x, void *params, gsl_matrix * J) {
  struct data *d = (struct data *) params;
  double A = gsl_vector_get(x, 0);
  double B = gsl_vector_get(x, 1);
  double C = gsl_vector_get(x, 2);
  double D = gsl_vector_get(x, 3);
  double w0 = gsl_vector_get(x, 4);
  double dw = gsl_vector_get(x, 5);
  double E = (x->size >= 8) ? gsl_vector_get(x, 6) : 0.0;
  double F = (x->size >= 8) ? gsl_vector_get(x, 7) : 0.0;
  double C2 = (x->size == 10) ? gsl_vector_get(x, 6) : 0.0;
  double D2 = (x->size == 10) ? gsl_vector_get(x, 7) : 0.0;
  double w02 = (x->size == 10) ? gsl_vector_get(x, 8) : 0.0;
  double dw2 = (x->size == 10) ? gsl_vector_get(x, 9) : 0.0;
  size_t i;

  for (i = 0; i < d->n; ++i) {
    double wi = d->w[i];
    double Xi = d->x[i];
    double Yi = d->y[i];

    double wa = w0*w0 - wi*wi;
    double wb = wi*dw;
    double z = wa*wa + wb*wb;

    double wa2 = w02*w02 - wi*wi;
    double wb2 = wi*dw2;
    double z2 = wa2*wa2 + wb2*wb2;

    gsl_matrix_set(J, 2*i, 0, -1);  // -dX/dA
    gsl_matrix_set(J, 2*i, 1, 0);   // -dX/dB
    gsl_matrix_set(J, 2*i+1, 0, 0);   // -dY/dA
    gsl_matrix_set(J, 2*i+1, 1, -1);  // -dY/dB

    switch (d->fit_func) {
     case OSCX_COFFS:
        gsl_matrix_set(J, 2*i, 2, -wa/z); // -dX/dC
        gsl_matrix_set(J, 2*i, 3, -wb/z); // -dX/dD
        gsl_matrix_set(J, 2*i, 4, -2*C*w0/z + (C*wa+D*wb)/z/z * 4*wa*w0); // -dX/d(f0)
        gsl_matrix_set(J, 2*i, 5, -D*wi/z   + (C*wa+D*wb)/z/z * 2*wb*wi); // -dX/d(df)

        gsl_matrix_set(J, 2*i+1, 2, +wb/z);  // -dY/dC
        gsl_matrix_set(J, 2*i+1, 3, -wa/z); // -dY/dD
        gsl_matrix_set(J, 2*i+1, 4, -2*D*w0/z + (D*wa-C*wb)/z/z * 4*wa*w0); // -dY/d(f0)
        gsl_matrix_set(J, 2*i+1, 5, +C*wi/z   + (D*wa-C*wb)/z/z * 2*wb*wi); // -dY/d(df)
        break;
     case OSCX_LOFFS:
        gsl_matrix_set(J, 2*i, 2, -wa/z); // -dX/dC
        gsl_matrix_set(J, 2*i, 3, -wb/z); // -dX/dD
        gsl_matrix_set(J, 2*i, 4, -2*C*w0/z + (C*wa+D*wb)/z/z * 4*wa*w0 - E); // -dX/d(f0)
        gsl_matrix_set(J, 2*i, 5, -D*wi/z   + (C*wa+D*wb)/z/z * 2*wb*wi);     // -dX/d(df)
        gsl_matrix_set(J, 2*i, 6, w0-wi); // -dX/dE
        gsl_matrix_set(J, 2*i, 7, 0);     // -dX/dF

        gsl_matrix_set(J, 2*i+1, 2, +wb/z);  // -dY/dC
        gsl_matrix_set(J, 2*i+1, 3, -wa/z); // -dY/dD
        gsl_matrix_set(J, 2*i+1, 4, -2*D*w0/z + (D*wa-C*wb)/z/z * 4*wa*w0 - F); // -dY/d(f0)
        gsl_matrix_set(J, 2*i+1, 5, +C*wi/z   + (D*wa-C*wb)/z/z * 2*wb*wi);     // -dY/d(df)
        gsl_matrix_set(J, 2*i+1, 6, 0);     // dY/dE
        gsl_matrix_set(J, 2*i+1, 7, w0-wi); // dY/dF
        break;
     case OSCV_COFFS:
        gsl_matrix_set(J, 2*i, 2, -wi*wb/z); // -dX/dC
        gsl_matrix_set(J, 2*i, 3, +wi*wa/z); // -dX/dD
        gsl_matrix_set(J, 2*i, 4, -wi*(-2*D*w0/z + (D*wa-C*wb)/z/z * 4*wa*w0)); // -dX/d(f0)
        gsl_matrix_set(J, 2*i, 5, -wi*(+C*wi/z   + (D*wa-C*wb)/z/z * 2*wb*wi));     // -dX/d(df)

        gsl_matrix_set(J, 2*i+1, 2, -wi*wa/z); // -dY/dC
        gsl_matrix_set(J, 2*i+1, 3, -wi*wb/z); // -dY/dD
        gsl_matrix_set(J, 2*i+1, 4, wi*(-2*C*w0/z + (C*wa+D*wb)/z/z * 4*wa*w0)); // -dY/d(f0)
        gsl_matrix_set(J, 2*i+1, 5, wi*(-D*wi/z   + (C*wa+D*wb)/z/z * 2*wb*wi));     // -dY/d(df)
        break;
     case OSCV_LOFFS:
        gsl_matrix_set(J, 2*i, 2, -wi*wb/z); // -dX/dC
        gsl_matrix_set(J, 2*i, 3, +wi*wa/z); // -dX/dD
        gsl_matrix_set(J, 2*i, 4, -wi*(-2*D*w0/z + (D*wa-C*wb)/z/z * 4*wa*w0) - E); // -dX/d(f0)
        gsl_matrix_set(J, 2*i, 5, -wi*(+C*wi/z   + (D*wa-C*wb)/z/z * 2*wb*wi));     // -dX/d(df)
        gsl_matrix_set(J, 2*i, 6, w0-wi); // -dX/dE
        gsl_matrix_set(J, 2*i, 7, 0);     // -dX/dF

        gsl_matrix_set(J, 2*i+1, 2, -wi*wa/z); // -dY/dC
        gsl_matrix_set(J, 2*i+1, 3, -wi*wb/z); // -dY/dD
        gsl_matrix_set(J, 2*i+1, 4, wi*(-2*C*w0/z + (C*wa+D*wb)/z/z * 4*wa*w0) - F); // -dY/d(f0)
        gsl_matrix_set(J, 2*i+1, 5, wi*(-D*wi/z   + (C*wa+D*wb)/z/z * 2*wb*wi));     // -dY/d(df)
        gsl_matrix_set(J, 2*i+1, 6, 0);     // dY/dE
        gsl_matrix_set(J, 2*i+1, 7, w0-wi); // dY/dF
        break;
     case DOSCX_COFFS:
        gsl_matrix_set(J, 2*i, 2, -wa/z); // -dX/dC
        gsl_matrix_set(J, 2*i, 3, -wb/z); // -dX/dD
        gsl_matrix_set(J, 2*i, 4, -2*C*w0/z + (C*wa+D*wb)/z/z * 4*wa*w0); // -dX/d(f0)
        gsl_matrix_set(J, 2*i, 5, -D*wi/z   + (C*wa+D*wb)/z/z * 2*wb*wi);     // -dX/d(df)
        gsl_matrix_set(J, 2*i, 6, -wa2/z2); // -dX/dC2
        gsl_matrix_set(J, 2*i, 7, -wb2/z2); // -dX/dD2
        gsl_matrix_set(J, 2*i, 8, -2*C2*w02/z2 + (C2*wa2+D2*wb2)/z2/z2 * 4*wa2*w02); // -dX/d(f02)
        gsl_matrix_set(J, 2*i, 9, -D2*wi/z2   + (C2*wa2+D2*wb2)/z2/z2 * 2*wb2*wi);   // -dX/d(df2)

        gsl_matrix_set(J, 2*i+1, 2, +wb/z); // -dY/dC
        gsl_matrix_set(J, 2*i+1, 3, -wa/z); // -dY/dD
        gsl_matrix_set(J, 2*i+1, 4, -2*D*w0/z + (D*wa-C*wb)/z/z * 4*wa*w0); // -dY/d(f0)
        gsl_matrix_set(J, 2*i+1, 5, +C*wi/z   + (D*wa-C*wb)/z/z * 2*wb*wi); // -dY/d(df)
        gsl_matrix_set(J, 2*i+1, 6, +wb2/z2); // -dY/dC2
        gsl_matrix_set(J, 2*i+1, 7, -wa2/z2); // -dY/dD2
        gsl_matrix_set(J, 2*i+1, 8, -2*D2*w02/z2 + (D2*wa2-C2*wb2)/z2/z2 * 4*wa2*w02); // -dY/d(f02)
        gsl_matrix_set(J, 2*i+1, 9, +C2*wi/z2   + (D2*wa2-C2*wb2)/z2/z2 * 2*wb2*wi);   // -dY/d(df2)
        break;
     case DOSCV_COFFS:
        gsl_matrix_set(J, 2*i, 2, -wi*wb/z); // -dX/dC
        gsl_matrix_set(J, 2*i, 3, +wi*wa/z); // -dX/dD
        gsl_matrix_set(J, 2*i, 4, -wi*(-2*D*w0/z + (D*wa-C*wb)/z/z * 4*wa*w0)); // -dX/d(f0)
        gsl_matrix_set(J, 2*i, 5, -wi*(+C*wi/z   + (D*wa-C*wb)/z/z * 2*wb*wi));     // -dX/d(df)

        gsl_matrix_set(J, 2*i, 6, -wi*wb2/z2); // -dX/dC2
        gsl_matrix_set(J, 2*i, 7, +wi*wa2/z2); // -dX/dD2
        gsl_matrix_set(J, 2*i, 8, -wi*(-2*D2*w02/z2 + (D2*wa2-C2*wb2)/z2/z2 * 4*wa2*w02)); // -dX/d(f02)
        gsl_matrix_set(J, 2*i, 9, -wi*(+C2*wi/z2   + (D2*wa2-C2*wb2)/z2/z2 * 2*wb2*wi));   // -dX/d(df2)

        gsl_matrix_set(J, 2*i+1, 2, -wi*wa/z); // -dY/dC
        gsl_matrix_set(J, 2*i+1, 3, -wi*wb/z); // -dY/dD
        gsl_matrix_set(J, 2*i+1, 4, wi*(-2*C*w0/z + (C*wa+D*wb)/z/z * 4*wa*w0)); // -dY/d(f0)
        gsl_matrix_set(J, 2*i+1, 5, wi*(-D*wi/z   + (C*wa+D*wb)/z/z * 2*wb*wi));     // -dY/d(df)
        gsl_matrix_set(J, 2*i+1, 6, -wi*wa2/z2); // -dY/dC2
        gsl_matrix_set(J, 2*i+1, 7, -wi*wb2/z2); // -dY/dD2
        gsl_matrix_set(J, 2*i+1, 8, wi*(-2*C2*w02/z2 + (C2*wa2+D2*wb2)/z2/z2 * 4*wa2*w02)); // -dY/d(f02)
        gsl_matrix_set(J, 2*i+1, 9, wi*(-D2*wi/z2   + (C2*wa2+D2*wb2)/z2/z2 * 2*wb2*wi));     // -dY/d(df2)
        break;
    }
  }

  return GSL_SUCCESS;
}

/*
// Additional derivatives for the accelerated method
// (see fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel)
int
func_fvv (const gsl_vector * x, const gsl_vector * v,
          void *params, gsl_vector * fvv)
{
  struct data *d = (struct data *) params;
  double a = gsl_vector_get(x, 0);
  double b = gsl_vector_get(x, 1);
  double c = gsl_vector_get(x, 2);
  double va = gsl_vector_get(v, 0);
  double vb = gsl_vector_get(v, 1);
  double vc = gsl_vector_get(v, 2);
  size_t i;

  for (i = 0; i < d->n; ++i) {
      double ti = d->t[i];
      double zi = (ti - b) / c;
      double ei = exp(-0.5 * zi * zi);
      double Dab = -zi * ei / c;
      double Dac = -zi * zi * ei / c;
      double Dbb = a * ei / (c * c) * (1.0 - zi*zi);
      double Dbc = a * zi * ei / (c * c) * (2.0 - zi*zi);
      double Dcc = a * zi * zi * ei / (c * c) * (3.0 - zi*zi);
      double sum;

      sum = 2.0 * va * vb * Dab +
            2.0 * va * vc * Dac +
                  vb * vb * Dbb +
            2.0 * vb * vc * Dbc +
                  vc * vc * Dcc;

      gsl_vector_set(fvv, i, sum);
    }

  return GSL_SUCCESS;
}
*/

/********************************************************************/

void
callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w) {

  /*

  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double avratio = gsl_multifit_nlinear_avratio(w);
  double rcond;

  (void) params; // not used

  // compute reciprocal condition number of J(x)
  gsl_multifit_nlinear_rcond(&rcond, w);

  fprintf(stderr, "iter %2zu: a = %.4f, b = %.4f, c = %.4f, |a|/|v| = %.4f cond(J) = %8.4f, |f(x)| = %.4f\n",
          iter,
          gsl_vector_get(x, 0),
          gsl_vector_get(x, 1),
          gsl_vector_get(x, 2),
          avratio,
          1.0 / rcond,
          gsl_blas_dnrm2(f));
  */
}

double
solve_system(gsl_vector *x, gsl_vector *xe, gsl_multifit_nlinear_fdf *fdf,
             gsl_multifit_nlinear_parameters *params) {

  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;

  const size_t max_iter = 200;
  const double xtol = 1.0e-10;
  const double gtol = 1.0e-10;
  const double ftol = 1.0e-10;

  const size_t n = fdf->n;
  const size_t p = fdf->p;

  gsl_multifit_nlinear_workspace *work =
    gsl_multifit_nlinear_alloc(T, params, n, p);
  gsl_vector * f = gsl_multifit_nlinear_residual(work);
  gsl_vector * y = gsl_multifit_nlinear_position(work);


  int info;
  double chisq0, chisq, rcond;
  size_t i;


  /* initialize solver */
  gsl_multifit_nlinear_init(x, fdf, work);

  /* store initial cost */
  gsl_blas_ddot(f, f, &chisq0);

  /* iterate until convergence */
  gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,
                              callback, NULL, &info, work);

  /* store final cost */
  gsl_blas_ddot(f, f, &chisq);

  /* store cond(J(x)) */
  gsl_multifit_nlinear_rcond(&rcond, work);

  gsl_vector_memcpy(x, y);

  /* compute parameter errors (see first example in
     https://www.gnu.org/software/gsl/doc/html/nls.html ) */
  {
    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    gsl_matrix *J = gsl_multifit_nlinear_jac(work);
    double c = sqrt(chisq / (n-p));
    gsl_multifit_nlinear_covar (J, 0.0, covar);

    for (i=0; i<p; i++)
      gsl_vector_set(xe, i, c*sqrt(gsl_matrix_get(covar,i,i)));

    gsl_matrix_free (covar);
  }

  /* print summary */
  /*
  fprintf(stderr, "summary from method '%s/%s'\n",
          gsl_multifit_nlinear_name(work),
          gsl_multifit_nlinear_trs_name(work));
  fprintf(stderr, "number of iterations: %zu\n",
          gsl_multifit_nlinear_niter(work));
  fprintf(stderr, "function evaluations: %zu\n", fdf->nevalf);
  fprintf(stderr, "Jacobian evaluations: %zu\n", fdf->nevaldf);
  fprintf(stderr, "reason for stopping: %s\n",
          (info == 1) ? "small step size" : "small gradient");
  fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
  fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq));
  */

  gsl_multifit_nlinear_free(work);
  return sqrt(chisq/n);
}

/********************************************************************/
// Find initial conditions by some trivial assumptions.

void
fit_res_init (const size_t n, const size_t p,
         double * freq, double * real, double * imag,
         double pars[MAXPARS], fit_func_t fit_func) {

  // points with min/max freq
  size_t ifmin=0, ifmax = 0;
  for (size_t i = 0; i<n; i++) {
    if (freq[i] < freq[ifmin]) ifmin = i;
    if (freq[i] > freq[ifmax]) ifmax = i;
  }

  // A,B - in the middle between these points:
  double A = (real[ifmin] + real[ifmax])/2;
  double B = (imag[ifmin] + imag[ifmax])/2;

  // E,F - slope of the line connecting first and list points
  double E = (real[ifmax] - real[ifmin])/(freq[ifmax] - freq[ifmin]);
  double F = (imag[ifmax] - imag[ifmin])/(freq[ifmax] - freq[ifmin]);

  // Find furthest point from the line connecting these points,
  // It should be the resonance.
  double dmax=0;
  size_t imax=0;
  for (size_t i = 0; i<n; i++) {
    double d = hypot(real[i] - real[ifmin] - (freq[i]-freq[ifmin])*E,
                     imag[i] - imag[ifmin] - (freq[i]-freq[ifmin])*F);
    if (d>dmax) {dmax=d; imax=i;}
  }
  double w0 = freq[imax];

  // Find min/max freq where distance > dmax/sqrt(2).
  // This is resonance width.
  size_t idmin=imax, idmax=imax;
  double d0 = dmax/sqrt(2.0);
  for (size_t i = 0; i<n; i++) {
    double d = hypot(real[i] - real[ifmin] - (freq[i]-freq[ifmin])*E,
                     imag[i] - imag[ifmin] - (freq[i]-freq[ifmin])*F);
    if (d>d0 && freq[i] < freq[idmin]) idmin = i;
    if (d>d0 && freq[i] > freq[idmax]) idmax = i;
  }
  if (idmin == idmax) {
    if (idmin>0)  idmin--;
    if (idmax<n-1) idmax++;
  }
  double dw = freq[idmax]-freq[idmin];

  // amplitudes (coord):
  double C,D;
  C = -freq[imax]*dw*(imag[imax]-B);
  D =  freq[imax]*dw*(real[imax]-A);

  // fill parameters
  pars[0] = A; pars[1] = B;

  switch (fit_func) {
    case OSCX_COFFS:
      pars[2] = C;  pars[3] = D;
      pars[4] = w0; pars[5] = dw;
      break;
   case OSCX_LOFFS:
      pars[2] = C;  pars[3] = D;
      pars[4] = w0; pars[5] = dw;
      pars[6] = E;  pars[7] = F;
      break;
   case OSCV_COFFS:
      pars[2] = D/w0; pars[3] = -C/w0;
      pars[4] = w0;   pars[5] = dw;
      break;
   case OSCV_LOFFS:
      pars[2] = D/w0; pars[3] = -C/w0;
      pars[4] = w0;   pars[5] = dw;
      pars[6] = E;    pars[7] = F;
      break;
   case DOSCX_COFFS:
      pars[2] = C;     pars[3] = D;
      pars[4] = w0+dw; pars[5] = dw;
      pars[6] = C;     pars[7] = D;
      pars[8] = w0-dw; pars[9] = dw;
      break;
   case DOSCV_COFFS:
      pars[2] = D/w0;  pars[3] = -C/w0;
      pars[4] = w0+dw; pars[5] = dw;
      pars[6] = D/w0;  pars[7] = -C/w0;
      pars[8] = w0-dw; pars[9] = dw;
      break;
  }
}

/********************************************************************/
// Fit resonance with Lorentzian curve
double
fit_res (const size_t n, const size_t p,
         double * freq, double * real, double * imag,
         double pars[MAXPARS], double pars_e[MAXPARS],
         fit_func_t fit_func) {

  gsl_vector *f = gsl_vector_alloc(2*n);
  gsl_vector *x  = gsl_vector_alloc(p);
  gsl_vector *xe = gsl_vector_alloc(p);
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();
  size_t i;
  struct data fit_data;

  fit_data.fit_func = fit_func;
  fit_data.n = n;
  fit_data.w = freq;
  fit_data.x = real;
  fit_data.y = imag;

  /* define function to be minimized */
  fdf.f = func_f;
  fdf.df = func_df; // NULL;
  fdf.fvv = NULL; //func_fvv;
  fdf.n = 2*n;
  fdf.p = p;
  fdf.params = & fit_data;

  /* starting point */
  for (i=0; i<p; i++) gsl_vector_set(x, i, pars[i]);

//  fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;
  fdf_params.trs = gsl_multifit_nlinear_trs_lm;
  double res = solve_system(x, xe, &fdf, &fdf_params);

  for (i=0; i<p; i++) pars[i]  = gsl_vector_get(x, i);
  for (i=0; i<p; i++) pars_e[i] = gsl_vector_get(xe, i);
  for (i=p; i<MAXPARS; i++) pars[i] = 0;
  for (i=p; i<MAXPARS; i++) pars_e[i] = 0;

  gsl_vector_free(f);
  gsl_vector_free(x);
  gsl_vector_free(xe);

  return res;
}

