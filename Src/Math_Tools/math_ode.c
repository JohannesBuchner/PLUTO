/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Ordinary differential equation solvers.

  \author A. Mignone (mignone@ph.unito.it)
  \date   June 20, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
    
/* ---------------------------------------------------------------------
     Standard Ordinary differential equation (ODE) solvers
   --------------------------------------------------------------------- */
   
#define ODE_NVAR_MAX 256
static int RK4Solve (double *, int, double, double, 
                     void (*rhs)(double, double *, double *));
static int CK45Solve (double *, int, double, double *, 
                      void (*rhs)(double, double *, double *), double);

/* ********************************************************************* */
void ODE_Solve(double *y0, int nvar, double xbeg, double xend,
               double dx, void (*rhs)(double, double *, double *), int method)
/*!
 *  Main driver for the numerical integration of a system of standard 
 *  ordinary differential equations (ODEs) of the type
 *  \f[
 *      \frac{dy}{dx} = f(x,y)
 *  \f]
 * where x is the independent variable and y[] is a vector of unknowns.
 *  
 * \param [in,out] y0  on input, y0 is an array containing the initial 
 *                     condition at x=xbeg; on output it is replaced by 
 *                     the value of the function at x=xend
 * \param [in] nvar    an integer number specifying the number of variable
 *                     (i.e. the dimension of y0[])
 * \param [in] xbeg    the initial integration point
 * \param [in] xend    the final integration point
 * \param [in] dx      the initial (suggested) step size
 * \param [in] rhs     a pointer to a function of the form rhs(x,*y,*f)
 *                     giving the right hand side \f$f(x,y)\f$ of the ODE. 
 * \param [in] method  an integer specifying the integration algorithm.
 *                     Possible choices are ODE_RK4 (fixed step size, 4th 
 *                     order Runge-Kutta), ODE_CK45 (adaptive step size,
 *                     5th order Cash-Karp algorithm).
 *                                        
 *********************************************************************** */
{
  int    nv, i, nsteps;
  double x0, x1;

  nsteps = (int)((xend - xbeg)/dx);
  
  x0 = xbeg;

/*printf ("[ODE_Solve] xbeg = %12.6e, xend = %12.6e\n",xbeg, xend); */

  while (x0 < xend){

  /* -- shorten step size to match final integration point -- */

    if (method == ODE_RK4){

      RK4Solve(y0, nvar, x0, dx, rhs);
      x0 += dx;

    }else if (method == ODE_CK45){
      int    err, k=0;
      double dx0;

      dx0 = dx; /* -- save initial step  size -- */

      while (++k) {
        if ((x0+dx) > xend) {
          dx = xend - x0;
        }
        x1  = x0 + dx;
        err = CK45Solve (y0, nvar, x0, &dx, rhs, 1.e-9);
        if (err == 0) {
          x0 = x1;
          break;
        }
        if (fabs(dx/dx0) < 1.e-4) {
          print ("! ODE_Solve: step size too small \n");
          QUIT_PLUTO(1);
        }
        if (k > 16384){
          print ("! ODE_Solve: too many steps\n");
          print ("! xbeg = %f, xend = %f\n",xbeg, xend);
          QUIT_PLUTO(1);
        }
      }
          
    }else{
      print ("! ODE_Solve: integration method unknown\n");
      QUIT_PLUTO(1);
    }
  }
}

/* ********************************************************************* */
int RK4Solve (double *y0, int nvar, double x0, 
              double dx, void (*rhs)(double, double *, double *))
/*
 *
 *
 *********************************************************************** */
{
  int    nv;
  double y1[ODE_NVAR_MAX];
  double k1[ODE_NVAR_MAX], k2[ODE_NVAR_MAX];
  double k3[ODE_NVAR_MAX], k4[ODE_NVAR_MAX];

/* -- step 1 -- */

  rhs (x0, y0, k1);
  for (nv = 0; nv < nvar; nv++) y1[nv] = y0[nv] + 0.5*dx*k1[nv];

/* -- step 2 -- */

  rhs (x0 + 0.5*dx, y1, k2);
  for (nv = 0; nv < nvar; nv++) y1[nv] = y0[nv] + 0.5*dx*k2[nv];

/* -- step 3 -- */

  rhs (x0 + 0.5*dx, y1, k3);
  for (nv = 0; nv < nvar; nv++) y1[nv] = y0[nv] + dx*k3[nv];

/* -- step 4 -- */

  rhs (x0 + dx, y1, k4);
  for (nv = 0; nv < nvar; nv++)
    y0[nv] += dx*(k1[nv] + 2.0*(k2[nv] + k3[nv]) + k4[nv])/6.0;

  return (0);
}

/* ********************************************************************* */
int CK45Solve (double *y0, int nvar, double x0, 
               double *dx0, void (*rhs)(double, double *, double *), double tol)
/*
 *     use explicit Cash-Karp 4-5 integrator
 *
 * On output (*dx0) is replaced with the suggested new time step.
 *********************************************************************** */
{
  int   i, n, rhs_err, isgn, useZ=0;

  double scrh, err, dx_shrink, dx_grow;
  double Yscal, dx;

  const double a2 = 0.2;
  const double a3 = 0.3;
  const double a4 = 0.6;
  const double a5 = 1.0;
  const double a6 = 0.875;

  const double c1 = 37.0/378.0;
  const double c2 = 0.0;
  const double c3 = 250.0/621.0;
  const double c4 = 125.0/594.0;
  const double c5 = 0.0;
  const double c6 = 512.0/1771.0; 

  const double cs1 = 2825.0/27648.0;
  const double cs2 = 0.0;
  const double cs3 = 18575.0/48384.0;
  const double cs4 = 13525.0/55296.0;
  const double cs5 = 277.0/14336.0;
  const double cs6 = 0.25;

  const double b21 =  0.2;
  const double b31 =  3.0/40.0 , b32 = 9.0/40.0;
  const double b41 =  0.3      , b42 = -0.9, b43 =  1.2;
  const double b51 = -11.0/54.0, b52 =  2.5, b53 = -70.0/27.0, b54 = 35.0/27.0;
  const double b61 = 1631.0/55296.0, b62 = 175.0/512.0, b63 = 575.0/13824.0, b64 = 44275.0/110592.0, b65 = 253.0/4096.0;
    
  double y1[ODE_NVAR_MAX], ys[ODE_NVAR_MAX];
  double y4th[ODE_NVAR_MAX], y5th[ODE_NVAR_MAX];
  double k1[ODE_NVAR_MAX], k2[ODE_NVAR_MAX], k3[ODE_NVAR_MAX]; 
  double k4[ODE_NVAR_MAX], k5[ODE_NVAR_MAX], k6[ODE_NVAR_MAX];

  dx = (*dx0);

/* -- check increment sign -- */

  if (dx < 0.0) isgn = -1;
  else          isgn =  1;

/* -- make a safety copy of the initial condition -- */

  for (n = 0; n < nvar; n++) y1[n] = y0[n];


  rhs (x0, y0, k1);
  for (n = 0; n < nvar; n++) ys[n] = y1[n] + dx*b21*k1[n];

  rhs (x0 + a2*dx, ys, k2);
  for (n = 0; n < nvar; n++) ys[n] = y1[n] + dx*(b31*k1[n] + b32*k2[n]);

  rhs (x0 + a3*dx, ys, k3);
  for (n = 0; n < nvar; n++) {
    ys[n] = y1[n] + dx*(b41*k1[n] + b42*k2[n] + b43*k3[n]);
  }

  rhs (x0 + a4*dx, ys, k4);
  for (n = 0; n < nvar; n++) {
    ys[n] = y1[n] + dx*(b51*k1[n] + b52*k2[n] + b53*k3[n] + b54*k4[n]);
  }

  rhs (x0 + a5*dx, ys, k5);
  for (n = 0; n < nvar; n++) {
    ys[n] = y1[n] + dx*(b61*k1[n] + b62*k2[n] + b63*k3[n] 
                                  + b64*k4[n] + b65*k5[n]);
  }

  rhs (x0 + a6*dx, ys, k6);

/* -- compute 5th and 4th order solutions -- */

  for (n = 0; n < nvar; n++) {
    y5th[n] = y1[n] + dx*(c1*k1[n] + c2*k2[n] + c3*k3[n] 
                        + c4*k4[n] + c5*k5[n] + c6*k6[n]);
    y4th[n] = y1[n] + dx*(cs1*k1[n] + cs2*k2[n] + cs3*k3[n] 
                        + cs4*k4[n] + cs5*k5[n] + cs6*k6[n]);
  }

/* --------------------------------------------------------
    compute error

    the following should give an absolute error if the
    solution is close to zero and a relative one in the 
    limit of large |Y|
   -------------------------------------------------------- */

  err   = 0.0;
  for (n = 0; n < nvar; n++) {
    scrh = fabs(y5th[n] - y4th[n])/(1.0 + fabs(y1[n]));
    err  = MAX(err, scrh);
  }

  err /= tol;

  if (err < 1.0){  /* -- ok, accept step -- */

    err = MAX(err, 1.e-18);

  /* -- provide an estimate for next dx -- */

    dx_grow  = 0.9*fabs(dx)*pow(err, -0.2);
    dx_grow  = MIN(dx_grow, 5.0*fabs(dx)); /* -- do not increase more than 5 -- */
    *dx0     = isgn*dx_grow;

    for (n = 0; n < nvar; n++) y0[n] = y1[n] = y5th[n];
    return 0;

  }else{   /* -- shrink dx and redo time step -- */

    dx_shrink = 0.9*fabs(dx)*pow(err, -0.25);
    dx        = MAX(dx_shrink, 0.05*fabs(dx)); /* -- do not decrease more than 20 -- */
    *dx0      = isgn*dx;
    return 1;
  } 
}    
#undef ODE_NVAR_MAX
