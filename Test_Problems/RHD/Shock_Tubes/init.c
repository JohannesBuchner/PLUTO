/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Relativistic shock tube problems.

  Shock tube problems are widely used for code benchmarking since they
  directly probe the solver's ability in resolving elementary waves and their
  propagation.
  In addition, numerical results can be compared with analytical solutions.
  
  Here we provide a simple suite of shock tube for the relativistic hydro
  module, by following Mignone \& Bodo (2005).  
  \f[
     \left(\rho,\, v_x,\, p\right) = \left\{\begin{array}{ll}
       \left(\rho,\, v_x,\, p\right)_L & \quad\mathrm{for}\quad x < 0.5
       \\ \noalign{\medskip}
       \left(\rho,\, v_x,\, p\right)_R & \quad\mathrm{otherwise} 
      \end{array}\right.
  \f]  
  The four available configurations correspond to Problem 1-4 of Mignone
  \& Bodo (2005):
    
  <CENTER>
  Test | rho  | vx  | p   |   Gamma
  -----|------|-----|-----|--------------
   1L  |  1   | 0.9 | 1   |    4/3
   1R  |  1   | 0.0 | 10  |    -
   2L  |  1   |-0.6 | 10  |    5/3
   2R  |  10  | 0.5 | 20  |    - 
   3L  |  10  |-0.6 | 40/3|    5/3
   3R  |  1   | 0.0 | 0   |    -
   4L  |  1   | 0.0 | 1.e3|    5/3
   4R  |  1   | 0.0 |1.e-2|    -
  </CENTER>

  Results are shown in the four figures below:

  \image html rhd_shock_tube.01.jpg "Flow profiles for the 1st shock tube problem"
  \image html rhd_shock_tube.02.jpg "Flow profiles for the 2nd shock tube problem"
  \image html rhd_shock_tube.03.jpg "Flow profiles for the 3rd shock tube problem"
  \image html rhd_shock_tube.04.jpg "Flow profiles for the 4th shock tube problem"

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 10, 2015

  \b References
     - "An HLLC Riemann solver for relativistic flows - I. Hydrodynamics"
        Mignone \& Bodo, MNRAS (2005) 364, 126
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x, double y, double z)
/*
 *
 *
 *
 *********************************************************************** */
{
  g_gamma = g_inputParam[GAMMA_EOS];
  if (x < 0.5){
    us[RHO] = g_inputParam[DN_L];
    us[VX1] = g_inputParam[VX_L];
    us[VX2] = g_inputParam[VY_L];
    us[PRS] = g_inputParam[PR_L];
  }else{
    us[RHO] = g_inputParam[DN_R];
    us[VX1] = g_inputParam[VX_R];
    us[VX2] = g_inputParam[VY_R];
    us[PRS] = g_inputParam[PR_R];
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
