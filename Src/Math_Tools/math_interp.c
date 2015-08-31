/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Miscellaneous functions for handling interpolation.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Jan 2, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void MonotoneSplineCoeffs (double *x, double *y, double *dydx, int n,
                           double *a, double *b, double *c, double *d)
/*!
 * Compute the cubic spline coefficients for interpolating a monotonic
 * dataset known at node values  f[k] = f(x[k]) and derivative 
 * dfdx(x[k]).
 * The resulting spline will be the smoothest curve that passes through 
 * the control points while preserving monotonocity of the dataset 
 * but not necessarily continuity of the second derivative.
 * Coefficients will be computed in each interval   x[k] < x < x[k+1].
 * 
 * The x[k] may not be equally spaced.
 *
 * \param [in]  x      1D array of abscissas
 * \param [in]  f      1D array of function values 
 * \param [in]  dfdx   1D array of function derivative 
 * \param [in]  n      the number of function values.
 * \param [out] a      the x^3 cubic coefficient
 * \param [out] b      the x^2 cubic coefficient
 * \param [out] c      the x   cubic coefficient
 * \param [out] d      the (1) cubic coefficient
 * 
 * \b Reference:
 *    - "Monotonic cubic spline interpolation", 
 *      G. Wolberg, I. Alfy, 
 *      Proceedings of Computer Graphics International (1999)
 *    - "An energy-minimization framework for monotonic cubic
 *       spline interpolation", Wolberg \& Alfy,
 *       J. of comput. and applied math. (2002) 143, 145
 * 
 *********************************************************************** */
{
  int    k, nfail = 0;
  char   cm, c1, c2, c2a, c2b, c2c;
  double m, alpha, beta, dx;

  for (k = 0; k < n-1; k++){
    m     = (y[k+1] - y[k])/(x[k+1] - x[k]);
    alpha = dydx[k]/m;
    beta  = dydx[k+1]/m;

  /* -----------------------------------------------------
      Check monotonicity conditions. Quit if not verified
     ----------------------------------------------------- */

    cm  = (DSIGN(dydx[k]) == DSIGN(dydx[k+1])) && (DSIGN(dydx[k]) == DSIGN(m));
    c1  = (alpha + beta - 2.0)     <= 0.0; 
    c2  = (alpha + 2.0*beta - 2.0)  > 0.0;
    c2a = (2.0*alpha + beta - 3.0) <= 0.0;
    c2b = (alpha + 2.0*beta - 3.0) <= 0.0;
    c2c = (alpha*alpha + alpha*(beta-6.0) + (beta-3.0)*(beta-3.0)) < 0.0;
    if ( (cm && c1) || (c2 && cm && (c2a || c2b || c2c))){
      dx = x[k+1] - x[k];
      a[k] = dx*(-2.0*m + dydx[k]     + dydx[k+1]);
      b[k] = dx*( 3.0*m - 2.0*dydx[k] - dydx[k+1]);
      c[k] = dx*dydx[k];
      d[k] = y[k];
    }else{
      print1 ("! MonotoneSplineCoeffs(): monotonicity condition not ");
      print1 ("satisifed in SPLINE1 \n");
      print1 ("cm = %d, c1 = %d, c2 = %d, c2a = %d\n",cm,c1,c2,c2a);
      QUIT_PLUTO(1);
    }

  /* ---------------------------------------------------------
      Reset a = b = 0 if the interpolant does not have
      enough curvature.
      This is done by checking the size of |y''/y| at x = 0.5
      (midpoint).
     --------------------------------------------------------- */

    alpha = 6.0*a[k]*0.5 + 2.0*b[k];
    beta  = d[k] + 0.5*(c[k] + 0.5*(b[k] + 0.5*a[k]));
    m     = fabs(alpha/beta);
    if (m < 1.e-16){
      nfail++;
      a[k] = b[k] = 0.0;
    }
  }

}
/* ********************************************************************* */
void SplineCoeffs (double *x, double *f, double dfL, double dfR, int n,
                   double *a, double *b, double *c, double *d)
/*!
 * Compute the cubic spline coefficients for interpolating a
 * dataset (not necessarily monotonic) known at node values  
 * f[k] = f(x[k]).
 * The resulting spline will be continuous up to the second derivative.
 * Coefficients will be computed in each interval   x[k] < x < x[k+1].
 * 
 * The x[k] may not be equally spaced.
 *
 * \param [in]  x      1D array of abscissas
 * \param [in]  f      1D array of function values 
 * \param [in]  dfL    derivative at the leftmost node 
 * \param [in]  dfR    derivative at the rightmost node 
 * \param [in]  n      the number of function values.
 * \param [out] a      the x^3 cubic coefficient
 * \param [out] b      the x^2 cubic coefficient
 * \param [out] c      the x   cubic coefficient
 * \param [out] d      the (1) cubic coefficient
 * 
 * \b Reference:
 *    - "Numerical Recipe in C, second edition", Section 3.3.
 * 
 *********************************************************************** */
{
  int    i, nfail=0;
  double p, sig, qn, un, dx;
  double *u, *d2f;
  double xm, xp, det, aa,bb,cc; /* coefficients of the derivative (parabola) */
  
  if (u == NULL) {
    u   = ARRAY_1D(n, double);
    d2f = ARRAY_1D(n, double);
  }

  d2f[0] = -0.5;
  u[0]   = (3.0/(x[1] - x[0])*(f[1] - f[0])/(x[1]-x[0]) - dfL);
    
/* ---------------------------------------------------------------
    Start decomposition loop of the tridiagonal algorithm.
    d2f[] and u[] are used for temporary storage of the 
    decomposed factors.
   --------------------------------------------------------------- */ 
   
  for (i = 1; i <= n-2;i++) { 
    sig    = (x[i] - x[i-1])/(x[i+1] - x[i-1]);
    p      = sig*d2f[i-1] + 2.0;
    d2f[i] = (sig - 1.0)/p;
    u[i]   = (f[i+1] - f[i])/(x[i+1]-x[i]) - (f[i] - f[i-1])/(x[i]-x[i-1]);
    u[i]   = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p;
  }
  qn = 0.5;
  un = (3.0/(x[n-1]-x[n-2]))*(dfR - (f[n-1] - f[n-2])/(x[n-1]-x[n-2]));
  
  d2f[n-1] = (un - qn*u[n-2])/(qn*d2f[n-2] + 1.0);
  for (i = n-2; i >= 0; i--) d2f[i] = d2f[i]*d2f[i+1] + u[i];

/* -----------------------------------------------------
    Construct spline coefficients, Eq. [3.3.3], [3.3.4] 
   ----------------------------------------------------- */

  for (i = 0; i < n-1; i++){
    dx   = x[i+1] - x[i]; 
    a[i] = dx*dx/6.0*(d2f[i+1] - d2f[i]);
    b[i] = dx*dx*0.5*d2f[i];
    c[i] = f[i+1] - f[i] - dx*dx/6.0*(2.0*d2f[i] + d2f[i+1]);
    d[i] = f[i];

  /* -- check monotonicity (quit if not satisfied) -- */
      
    aa = 3.0*a[i];
    bb = 2.0*b[i];
    cc = 1.0*c[i];
    det = bb*bb - 4.0*aa*cc;
    if (det >= 0.0){
      xm = (-bb + sqrt(det))/(2.0*aa);
      xp = (-bb - sqrt(det))/(2.0*aa);

      if ((xm > 0.0 && xm < 1.0) || (xp > 0.0 && xp < 1.0)){
        printf ("! SplineCoeffs(): cubic is not monotonic, in [%8.3e, %8.3e]\n",
                 x[i], x[i+1]);
        printf ("! xm = %8.3e, xp = %8.3e\n",xm,xp);
        QUIT_PLUTO(1);
      }   
    }

  /* ---------------------------------------------------------
      Reset a = b = 0 if the interpolant does not have
      enough curvature.
      This is done by checking the size of |y''/y| at x = 0.5
      (midpoint).
     --------------------------------------------------------- */

    xm = 6.0*a[i]*0.5 + 2.0*b[i];
    xp = d[i] + 0.5*(c[i] + 0.5*(b[i] + 0.5*a[i]));
    det     = fabs(xm/xp);
    if (det < 1.e-16){
      nfail++;
      a[i] = b[i] = 0.0;
    }
  }
}

