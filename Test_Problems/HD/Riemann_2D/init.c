/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Two-Dimensional Riemann problem.

  Sets the initial condition for a 2D Riemann problem in terms of 4 four
  constant states at each corner of the computational domain.

  For example, the top right corner, \f$x > 0\f$, \f$y > 0\f$, is denoted
  with \c PP and its initial values are:

  \f[
    \rho = \rho_{\rm PP}
    \,,\quad\quad
    P = P_{\rm PP}
    \,,\quad\quad
    v_x = v_{x\,{\rm PP}}
    \,,\quad\quad
    v_y = v_{y\,{\rm PP}}\,.
  \f]

  Similar conditions hold for 
  the bottom right corner (\f$x > 0\f$, \f$y < 0\f$; \c PM),
  the top left (\f$x < 0\f$, \f$y > 0\f$; \c MP), and
  the bottom left (\f$x < 0\f$, \f$y < 0\f$; \c MM).
  The 4 flow quantities of each corner are defined by the 16 parameters
  that are read from pluto.ini.
  - <tt>DN_PP, PR_PP, VX_PP, VY_PP</tt>
  - <tt>DN_MP, PR_MP, VX_MP, VY_MP</tt>
  - <tt>DN_PM, PR_PM, VX_PM, VY_PM</tt>
  - <tt>DN_MM, PR_MM, VX_MM, VY_MM</tt>


  - Configurations #01, #02 and #05 involve the interaction between four contact 
    discontinuities (see fig. below).
  - Configurations #03 and #04 (with AMR) involve the interaction of shocks.

  \image html hd_riemann2D.05.jpg "Final state for configuration #05."

  \author A. Mignone (mignone@ph.unito.it)
  \date    Sept 15, 2014
  
  \b References
     - C. W. Schulz-Rinne, "Classification of the Riemann problem for
        two-dimensional gas dynamics", SIAM J. Math. Anal., 24 
        (1993), pp. 76-88.
     - C. W. Schulz-Rinne, J. P. Collins and H. M. Glaz, 
       "Numerical solution of the Riemann problem for two-dimensional gas
        dynamics", SIAM J. Sci. Comp., 14 (1993), pp. 1394-1414.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *********************************************************************** */
{
  double x,y;

  x = x1;
  y = x2;
  #if EOS == IDEAL
   g_gamma = 1.4;
  #elif EOS == ISOTHERMAL
   g_isoSoundSpeed = 1.0;
 g_isoSoundSpeed = 2.0;
  #endif

  if (x > 0.0 && y > 0.0){
    us[RHO] = g_inputParam[DN_PP];
    #if EOS == IDEAL
     us[PRS] = g_inputParam[PR_PP];
    #endif
    us[VX1] = g_inputParam[VX_PP];
    us[VX2] = g_inputParam[VY_PP];
  }else if(x < 0.0 && y > 0.0){
    us[RHO] = g_inputParam[DN_MP];
    #if EOS == IDEAL
     us[PRS] = g_inputParam[PR_MP];
    #endif
    us[VX1] = g_inputParam[VX_MP];
    us[VX2] = g_inputParam[VY_MP];
  }else if(x < 0.0 && y < 0.0){
    us[RHO] = g_inputParam[DN_MM];
    #if EOS == IDEAL
     us[PRS] = g_inputParam[PR_MM];
    #endif
    us[VX1] = g_inputParam[VX_MM];
    us[VX2] = g_inputParam[VY_MM];
  }else if(x > 0.0 && y < 0.0){
    us[RHO] = g_inputParam[DN_PM];
    #if EOS == IDEAL
     us[PRS] = g_inputParam[PR_PM];
    #endif
    us[VX1] = g_inputParam[VX_PM];
    us[VX2] = g_inputParam[VY_PM];
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
 *********************************************************************** */
{ }

