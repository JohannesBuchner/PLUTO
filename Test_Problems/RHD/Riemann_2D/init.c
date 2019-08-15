/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Two-dimensional Riemann problem.

  Sets the initial condition for the 2D Riemann problem 
  described in Mignone et al. (2005).
  The computational domain is initially divided into four 
  states and the outcoming wave pattern involves the formtation of
  2 shocks and two contact waves.

  \image html rhd_riemann2D.01.jpg "Final state for configuration #01."

  \authors A. Mignone (mignone@ph.unito.it)
  \date   July 09, 2014

  \b Reference:
     - "The Piecewise Parabolic Method for Multidimensional Relativistic 
        Fluid Dynamics", Mignone, Plewa & Bodo, ApJS (2005)
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double x, y, scrh;

  x = x1;
  y = x2;

  us[VX1] = us[VX2] = 0.0;
  if (x > 0.0 && y > 0.0){
    us[RHO] = 5.477875e-3;
    us[PRS] = 2.762987e-3;
    us[VX1] = 0.0;
    us[VX2] = 0.0;
  }else if(x < 0.0 && y > 0.0){
    us[RHO] = 0.1;
    us[PRS] = 1.0;
    us[VX1] = 0.99;
    us[VX2] = 0.0;
  }else if(x < 0.0 && y < 0.0){
    us[RHO] = 0.5;
    us[PRS] = 1.0;
    us[VX1] = 0.0;
    us[VX2] = 0.0;
  }else if(x > 0.0 && y < 0.0){
    us[RHO] = 0.1;
    us[PRS] = 1.0;
    us[VX1] = 0.0;
    us[VX2] = 0.99;
  }   

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
 *
 *********************************************************************** */
{ }

