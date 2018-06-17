/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief X-point test particle acceleration.

  Test-particle acceleration near an X-type magnetic reconnection region.
  (see sect. 4.6 of Mignone et al. 2018).


  \author G. Mattia, 
          A. Mignone (mignone@ph.unito.it)
  \date   May 13, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x, double y, double z)
/*! 
 *
 *
 *********************************************************************** */
{
  double alpha = g_inputParam[ALPHA];
  double beta  = g_inputParam[BETA];
  double E[3], B[3];
  double coor[3] = {x,y,z};

  Particles_CR_EMFields(coor, E, B);

  v[RHO] = 1.0;

  v[BX1] = B[IDIR];
  v[BX2] = B[JDIR];
  v[BX3] = 0.0;

  v[PRS] = 1.0;

  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = alpha*y*y/2.0 - beta*x*x/2.0;

  g_smallPressure = 1.e-5;
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 *
 *
 *********************************************************************** */
{
  
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *
 *********************************************************************** */
{
  
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif

/* ********************************************************************* */
void Particles_CR_EMFields(double *coor, double *E, double *B)
/*
 * Compute Electromagnetic fields as a function of the coordinate.
 *********************************************************************** */
{
  double alpha = g_inputParam[ALPHA];
  double beta  = g_inputParam[BETA];

  double Bz = g_inputParam[BMAG_Z];
  double Ez = g_inputParam[EMAG_Z];

  double x = coor[IDIR];
  double y = coor[JDIR];
  
/* -- Experiments -- */
/*
  B[IDIR] = 1.0;
  B[JDIR] = 0.0;
  B[KDIR] = 0.0;

  E[IDIR] = 0.0;
  E[JDIR] = 0.0;
  E[KDIR] = 0.3*sin(x*y);
*/
/* -- Mori et al 1998 -- */

  B[IDIR] = alpha*y;
  B[JDIR] = beta*x;
  B[KDIR] = Bz;
  
  E[IDIR] = 0.0;
  E[JDIR] = 0.0;
  E[KDIR] = Ez;
  
/* -- Zharkova 2011 (page 386-387, exepct for a minus sign) -- */  
/*
  B[IDIR] = tanh(-alpha*y);
  B[JDIR] = beta*x;
  B[KDIR] = Bz;

  E[IDIR] = 0.0;
  E[JDIR] = 0.0;
  E[KDIR] = Ez;
*/
}
