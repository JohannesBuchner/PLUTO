/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Collection of handy numerical math tools.

  This file provides a number of standard numerical routines 
  to achieve simple basic tasks such as

  - LU decomposition functions;
  - Numerical quadrature;
  - Ordinary differential equation solver;
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   March 28, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
double GaussQuadrature(double (*func)(double, void *), void * par, double xb, double xe,
                       int nstep, int order)
/*!
 * Perform numerical quadrature of the function f(x) between 
 * the lower bound xb and upper bound xe by subdividing the interval 
 * into 'nstep' steps.
 * A 3 or 5-point Gaussian quadrature rule is used depending on the 
 * input variable order (=3 or =5)
 * 
 * \param [in] *func   a pointer to the function func(x) (returning double)
 *                     to be integrated
 * \param [in] xb      the lower interval bound
 * \param [in] xe      the upper interval bound
 * \param [in] nstep   the number of sub-intervals into which the 
 *                     original interval [xb,xe] has to be divided
 * \param [in] order   the number of Gaussian points (only 3 or 5)    
 *                 
 *********************************************************************** */
{
  int    i, n;
  double w[8], z[8], x;
  double I, Isub;
  double xb0, xe0, dx;
  
  if (order == 3){
    double s3 = sqrt(3.0/5.0);

    z[0] = -s3;  w[0] = 5.0/9.0;
    z[1] = 0.0;  w[1] = 8.0/9.0;
    z[2] =  s3;  w[2] = 5.0/9.0;
  }else if (order == 5){
    double s1, s7;
    
    s1 = sqrt(10.0/7.0);
    s7 = sqrt(70.0);
    z[0] = -1.0/3.0*sqrt(5.0 - 2.0*s1); w[0] = (322.0 + 13.0*s7)/900.0;
    z[1] = -z[0]; w[1] = w[0];

    z[2] = -1.0/3.0*sqrt(5.0 + 2.0*s1); w[2] = (322.0 - 13.0*s7)/900.0;
    z[3] = -z[2]; w[3] = w[2];

    z[4] = 0.0; w[4] = 128.0/225.0;
  }else{
    print ("! GaussQuadrature: order must be either 3 or 5\n");
    QUIT_PLUTO(1);
  }
  
  if (nstep <= 0){
    print ("! GaussQuadrature: nstep must be > 0\n");
    QUIT_PLUTO(1);
  }
  
  xb0 = xb; xe0 = xe; /* save original interval endpoints */
  dx = (xe - xb)/(double)nstep; /* sub-interval length */
  I  = 0.0;
  for (i = 0; i < nstep; i++){
    xb = xb0 + i*dx;
    xe = xb + dx;
    
    Isub = 0.0;  /* intgrate sub-interval */
    for (n = 0; n < order; n++){
      x     = 0.5*(xe - xb)*z[n] + (xe + xb)*0.5;
      Isub += w[n]*func(x, par);
    }    
    Isub *= 0.5*(xe - xb);
    I    += Isub;
  }
  
  return I;
}

