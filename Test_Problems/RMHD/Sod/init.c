#include "pluto.h"
/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Collection of relativistic MHD shock-tube problems.

  This directory contains several 1D shock-tube configurations 
  commonly used for benchmarking of numerical methods.
  There are 16  input parameters that can be used to specify entirely
  the left and right states with respect to the discontinuity initially
  placed at \c x=1/2.
  For the left state one has:

  - <tt>g_inputParam[RHO_LEFT]</tt>: density for the left state
  - <tt>g_inputParam[VX_LEFT]</tt>: normal (x) velocity for the left state
  - <tt>g_inputParam[VY_LEFT]</tt>: tangential (y) velocity for the left state
  - <tt>g_inputParam[VZ_LEFT]</tt>: tangential (z) velocity for the left state
  - <tt>g_inputParam[BY_LEFT]</tt>: tangential (y) magnetic field for the left state
  - <tt>g_inputParam[BZ_LEFT]</tt>: tangential (z) magnetic field for the left state
  - <tt>g_inputParam[PR_LEFT]</tt>: pressure in the left state
  
  Similar parameters are used to specify the right states.
  In addition, 
  
  - <tt>g_inputParam[BX_CONST]</tt>: normal component of magnetic field (this 
    cannot have jump)
  - <tt>g_inputParam[GAMMA_EOS]</tt>: the specific heat ratio

  The IDEAL equation of state is used.
  The available configurations are taken from 

  <CENTER>
  Conf | Reference
  -----|-----------------------
  01-02| 1st problem in [MB06]
  03-04| 2nd problem in [MB06]
  05-06| 3rd problem in [MB06]
  07-08| 4th problem in [MB06]
  09-10| 2nd problem in [MUB09]
  11-12| 4th problem (generic Alfven test) in [MUB09]
  13-14| [MMB05]
  </CENTER>
  
  \image html rmhd_sod.11.jpg "Results for the generic Alfven test at t = 0.5 on 1600 grid zones for conf. #11."
  	 
  \authors A. Mignone (mignone@ph.unito.it)\n
  \date    Oct 6, 2014 16, 2014

  \b References 
     - [MB06]: "An HLLC Riemann solver for relativistic flows -- 
       II. Magnetohydrodynamics", Mignone & Bodo MNRAS(2006) 368,1040 
     - [MUB09]: "A five-wave Harten-Lax-van Leer Riemann solver for 
       relativistic magnetohydrodynamics", Mignone et al., MNRAS(2009) 393,1141 
     - [MMB05]: "Relativistic MHD simulations of jets with toroidal 
       magnetic fields", Mignone, Massaglia & Bodo, SSRv (2005), 121,21
     - Balsara, ApJS (2001), 132, 83
     - Del Zanna, Bucciantini, Londrillo, A&A (2003) 
*/
/* ///////////////////////////////////////////////////////////////////// */

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  g_gamma = g_inputParam[GAMMA_EOS];

  if (x1 < 0.5){
    v[RHO] = g_inputParam[RHO_LEFT];
    v[VX1] = g_inputParam[VX_LEFT];
    v[VX2] = g_inputParam[VY_LEFT];
    v[VX3] = g_inputParam[VZ_LEFT];
    v[BX1] = g_inputParam[BX_CONST];
    v[BX2] = g_inputParam[BY_LEFT];
    v[BX3] = g_inputParam[BZ_LEFT];
    v[PRS] = g_inputParam[PR_LEFT];
  }else{
    v[RHO] = g_inputParam[RHO_RIGHT];
    v[VX1] = g_inputParam[VX_RIGHT];
    v[VX2] = g_inputParam[VY_RIGHT];
    v[VX3] = g_inputParam[VZ_RIGHT];    
    v[BX1] = g_inputParam[BX_CONST];
    v[BX2] = g_inputParam[BY_RIGHT];
    v[BX3] = g_inputParam[BZ_RIGHT];    
    v[PRS] = g_inputParam[PR_RIGHT]; 
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
{
}
