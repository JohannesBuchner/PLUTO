/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Relativistic magnetized blast wave.

  Set the initial condition for the relativistic magnetized blast wave 
  problem in 2D or 3D. 
  It consists of a highly pressurized region inside a circle (in 2D) or
  a sphere (in 3D) embeddd in a static uniform medium with lower pressure.
  The magnetic field is constant and threads the whole computational domain.
  \f[
     (\rho,\, p) = \left\{\begin{array}{lcl}
        (\rho_{\rm in},\, p_{\rm in})  & \quad\mathrm{or}\quad & r < r_c
        \\ \noalign{\medskip}
        (\rho_{\rm out},\, p_{\rm out})  & \quad\mathrm{or}\quad & r \ge r_c
     \end{array}\right.
     \,,\qquad
      |\vec{B}| = B_0 
   \f]
  In 3D, a linear smoothing is applied in the region \f$ r_c<r<1\f$.
  The input parameters used in this problem are:
  
  -# <tt>g_inputParam[PRS_IN]</tt>:  pressure inside the initial circular (2D) 
                                    or spherical (3D) region.
  -# <tt>g_inputParam[PRS_OUT]</tt>: ambient pressure
  -# <tt>g_inputParam[RHO_OUT]</tt>: ambient density
  -# <tt>g_inputParam[BMAG]</tt>:    magnetic field intensity
  -# <tt>g_inputParam[THETA]</tt>:   angle between mag. field and z-axis 
                                    (Cartesian only )
  -# <tt>g_inputParam[PHI]</tt>:     angle between mag. field and xy-plane 
                                    (Cartesian only)
  -# <tt>g_inputParam[RADIUS]</tt>:  radius of the initial over-pressurized region.
 
  Note that a given choice of parameters can be re-scaled by an arbitrary 
  factor \f$\eta\f$  by letting
  \f$ \{\rho,\, p\} \to \eta^2\{\rho,\,p\},\, B \to \eta B \f$.

  The different configurations are: 
  - #01 and #02 are taken from  Del Zanna et al, A&A (2003) 400,397
  - #03 and #04 are taken from  Mignone et al, ApJS (2007), 170, 228

  - #05 and #06 are taken from Beckwith & Stone, ApJS (2011), 193, 6
    Strongly magnetized case, sec. 4.6 (Fig. 14)

  Strongly magnetized configurations can pass this test only by taking 
  some precautions (e.g. correcting total energy with staggered magnetic 
  field).

  \image html rmhd_blast.02.jpg "Density map (in log scale) for configuration #02"

  \authors A. Mignone (mignone@ph.unito.it)\n
  \date    Sept 16, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  double r, rc, theta, phi;
  double dc, pc, de, pe;

  g_gamma = 4./3.;
                        
  #if DIMENSIONS == 2
   r = sqrt(x1*x1 + x2*x2);
  #elif DIMENSIONS == 3
   r = sqrt(x1*x1 + x2*x2 + x3*x3);
  #endif 
  rc = g_inputParam[RADIUS];

  dc = g_inputParam[RHO_IN];
  pc = g_inputParam[PRS_IN];
  de = g_inputParam[RHO_OUT];
  pe = g_inputParam[PRS_OUT];

  if (r <= rc) {
    us[RHO] = dc;
    us[PRS] = pc;
  #if DIMENSIONS == 3
   }else if (r > rc && r < 1.0){
     us[RHO] = de*(r - rc)/(1.0 - rc) + dc*(r - 1.0)/(rc - 1.0);
     us[PRS] = pe*(r - rc)/(1.0 - rc) + pc*(r - 1.0)/(rc - 1.0);
  #endif
  }else{
    us[RHO] = de;
    us[PRS] = pe;
  }

  us[VX1] = us[VX2] = us[VX3] = 0.0;
  us[AX1] = us[AX2] = us[AX3] = 0.0;

  theta = g_inputParam[THETA]*CONST_PI/180.0;
  phi   =   g_inputParam[PHI]*CONST_PI/180.0;

  #if GEOMETRY == CARTESIAN
   us[BX1]  = g_inputParam[BMAG]*sin(theta)*cos(phi);
   us[BX2]  = g_inputParam[BMAG]*sin(theta)*sin(phi);
   us[BX3]  = g_inputParam[BMAG]*cos(theta);

   us[AX1] = 0.0;
   us[AX2] = us[BX3]*x1;
   us[AX3] = -us[BX2]*x1 + us[BX1]*x2;
  #elif GEOMETRY == CYLINDRICAL
   us[BX1]  = 0.0;
   us[BX2]  = g_inputParam[BMAG];
   us[BX3]  = 0.0;

   us[AX1] = us[AX2] = 0.0;
   us[AX3] = 0.5*us[BX2]*x1;
  #endif

  g_smallPressure = 1.e-6;  
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
 *
 *********************************************************************** */
{ }

