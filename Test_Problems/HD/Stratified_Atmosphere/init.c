/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Hydrostatic atmosphere in a smoothed gravitational potential.

  Set initial conditions for an hydrostatic atmosphere in cylindrical
  coordinates.
  Gravity is given as
  \f[
     g = \left\{\begin{array}{ll}
          -1/R^2              & \quad{\rm for} \quad R > 1 \\ \noalign{\medskip}
          aR + bR^2 + cR^3 & \quad{\rm for} \quad R < 1
     \end{array}\right.\,,
  \f]
  where \f$R = \sqrt{r^2 + z^2}\f$ is the spherical radius. The coefficients
  \f$a\f$, \f$b\f$ and \f$c\f$ are chosen to guarantee continuity of \f$g\f$,
  its first and second derivative (optionally).
  Density and pressure are tied by the isothermal condition
  \f$P = \rho/a\f$ so that the hydrostatic condition is
  \f[
    \frac{1}{\rho} \frac{d\rho}{dr} = ag\,,
  \f]
  with the normalization \f$\rho = 1\f$ at \f$R = 1\f$.

  The runtime parameters that are read from \c pluto.ini are 
  - <tt>g_inputParam[ALPHA]</tt>: sets the value of \f$a\f$;

  Configurations:
  - #1-4: 2D cylindrical
  - #5: 3D cartesian

  \b References: 

  \author A. Mignone (mignone@ph.unito.it)
  \date   Fen 28, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static real acf, bcf, ccf; /*  polynomial coefficient for g when r < 1 */

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double rs, scrh;

/* ------------------------------------------
    with this choice g, g' and g''
    will be continuous at R = 1
   ----------------------------------------- */

  acf = -10.0;
  bcf =  15.0;
  ccf = -6.0;

/* -----------------------------------
     with this choice g and g' 
     will be continuous
   ---------------------------------- */

  acf = -3.0;
  bcf =  2.0;
  ccf = 0.0;

  #if GEOMETRY == CARTESIAN
   rs = sqrt(x1*x1 + x2*x2 + x3*x3);
  #elif GEOMETRY == CYLINDRICAL 
   rs = sqrt(x1*x1 + x2*x2);
  #elif GEOMETRY == SPHERICAL
   rs = sqrt(x1*x1);
  #endif

  if (rs > 1.0){
    us[RHO] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0));
  }else{
    scrh   = 0.5*acf*(rs*rs - 1.0) + 1.0/3.0*bcf*(rs*rs*rs - 1.0)
              + 0.25*ccf*(rs*rs*rs*rs - 1.0);
    us[RHO] = exp(scrh*g_inputParam[ALPHA]);
  }

  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = 0.0;
  us[PRS] = us[RHO]/g_inputParam[ALPHA];
  us[TRC] = 0.0;
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *********************************************************************** */
{
  int   i, j, k;
  double *x1, *x2, *x3;
  double rs;

  x1 = grid->xgc[IDIR];
  x2 = grid->xgc[JDIR];
  x3 = grid->xgc[KDIR];

  if (side == X1_END) {

    X1_END_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN
       rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
      #elif GEOMETRY == CYLINDRICAL
       rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j]);
      #elif GEOMETRY == SPHERICAL
       rs = sqrt(x1[i]*x1[i]);
      #endif
      d->Vc[RHO][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0));
      EXPAND(d->Vc[VX1][k][j][i] = 0.0;   ,
             d->Vc[VX2][k][j][i] = 0.0;   ,
             d->Vc[VX3][k][j][i] = 0.0;)
      d->Vc[PRS][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0))/g_inputParam[ALPHA];
    }

  } else if (side == X2_END) {

    X2_END_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN
       rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
      #elif GEOMETRY == CYLINDRICAL
       rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j]);
      #elif GEOMETRY == SPHERICAL
       rs = sqrt(x1[i]*x1[i]);
      #endif
      d->Vc[RHO][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0));
      EXPAND(d->Vc[VX1][k][j][i] = 0.0;   ,
             d->Vc[VX2][k][j][i] = 0.0;   ,
             d->Vc[VX3][k][j][i] = 0.0;)
      d->Vc[PRS][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0))/g_inputParam[ALPHA];
    }

  } else if (side == X3_END) {   /* Only Cartesian */

    X3_END_LOOP(k,j,i){  
      rs = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
      d->Vc[RHO][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0));
      EXPAND(d->Vc[VX1][k][j][i] = 0.0;   ,
             d->Vc[VX2][k][j][i] = 0.0;   ,
             d->Vc[VX3][k][j][i] = 0.0;)
      d->Vc[PRS][k][j][i] = exp(g_inputParam[ALPHA]*(1.0/rs - 1.0))/g_inputParam[ALPHA];
    }
  }
}

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  double gs, rs;
  double acf, bcf, ccf;

  acf = -3.0;
  bcf =  2.0;
  ccf =  0.0;
  #if GEOMETRY == CARTESIAN
   rs = sqrt(x1*x1 + x2*x2 + x3*x3);
  #elif GEOMETRY == CYLINDRICAL
   rs = sqrt(x1*x1 + x2*x2);
  #endif

  if (rs > 1.0) gs = -1.0/rs/rs;
  else          gs = rs*(acf + rs*(bcf + rs*ccf));

  #if GEOMETRY == CARTESIAN
   g[IDIR] = gs*x1/rs;
   g[JDIR] = gs*x2/rs;
   g[KDIR] = gs*x3/rs;
  #elif GEOMETRY == CYLINDRICAL
   g[IDIR] = gs*x1/rs;
   g[JDIR] = gs*x2/rs;
   g[KDIR] = 0.0;
  #endif

}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  double rs, phi;
  double acf, bcf, ccf, C;

  acf = -3.0;
  bcf =  2.0;
  ccf =  0.0;

  #if GEOMETRY == CARTESIAN
   rs = sqrt(x1*x1 + x2*x2 + x3*x3);
  #elif GEOMETRY == CYLINDRICAL
   rs = sqrt(x1*x1 + x2*x2);
  #endif

  C = (0.5*acf + bcf/3.0 + ccf*0.25);  /* integration constant to make phi continuous */
  if (rs > 1.0) phi = -1.0/rs;
  else          phi = -rs*rs*(0.5*acf + rs*(bcf/3.0 + rs*ccf*0.25)) + C - 1.0;

  return phi;
}
#endif
