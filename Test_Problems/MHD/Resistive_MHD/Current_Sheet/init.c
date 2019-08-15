/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Reconnection test (Harris sheet) in 2D.
 
  In this setup - see also Section 5.3 of Mignone et al., ApJS (2012) 198:7 -
  we reproduce a 2D Harris current sheet with magnetic field profile given by 
  \f[
    B_x(y) = B_0 \tanh(y/l)
  \f] 
  where \c l is the half thickness of the layer.
  The density profile is given by
  \f[ 
    \rho(y) = \rho_0 \cosh^{-2}(y/l) + \rho_{\infty}
  \f]
  We use \f$ \rho_0=1 \f$ and \f$  \rho_{\infty}  = 0.2 \f$,
  following the guidelines of Birn et al., 2001, while
  \c l is user supplied. \n
  In order to achieve equilibrium with the magnetic pressure, 
  the thermal pressure is chosen to be \f$ p = c_s^2 \rho \f$, where
  \f$ c_s^2 = \frac{B_0^2}{2\rho_0} \f$.
  The initial equilibrium is pertubed by an additional magnetic field 
  defined as
  \f[
     \begin{array}{lcl}
      B_x(x,y) &=& \DS -\Psi_0\frac{\pi}{L_y}\cos\left(\frac{2\pi x}{L_x}\right)
                    \sin\left(\frac{\pi y}{L_y}\right),  \\ \noalign{\medskip}
     B_y(x,y) &=& \DS +\Psi_0 \frac{2\pi}{L_x}\sin\left(\frac{2\pi x}{L_x}\right)
                    \cos\left(\frac{\pi y}{L_y}\right).      
     \end{array}
  \f]  
  
  The Lundquist number \f$ S \f$ of a plasma is defined as 
  \f[
    S = \frac{v_A L}{\eta}
  \f]
  where \f$ v_A \f$ is the Alfvén velocity, 
  \f$ v_A = \DS \frac{B}{\sqrt{\rho}}\f$, \f$ L \f$ is a typical lenght scale, 
  and \f$ \eta \f$ the plasma resistivity.
  The reconnection rate \f$\mathcal{E} = \DS \frac{v_{in}}{v_{out}}\f$, with
  \f$ v_{in} \f$ and \f$ v_{out}\f$ the plasma inflow and outflow velocities,
  follows the Sweet-Parker scaling law 
  \f$\mathcal{E} \sim \frac{\delta}{L} \sim \frac{1}{\sqrt{S}}\f$.
  In this example several values of the resitivity \f$ \eta \f$, 
  that correspond to different values of the Lundquist number \f$ S \f$, 
  are provided. 
  The reconnection rate, calculated as the ratio \f$ \frac{\delta}{L} \f$
  (see Mignone et al., 2012) verifies the Sweet-Parker scaling in the range 
  \f$ \eta = 10^{-2} - 10^{-4} \f$ (see the first figure below).

  The input parameters (read from \c pluto.ini) for this test problem are:

  - <tt> g_inputParam[ETA]</tt>:   sets the value of resistivity \f$ \eta \f$;
  - <tt> g_inputParam[WIDTH]</tt>: sets the layer width \c l;
  - <tt> g_inputParam[PSI0]</tt>:  sets the amplitude of perturbation \f$\Psi_0\f$.

  \note
  - Configuration #02 employs a small width (l -\> 0, current-sheet)
    large resistivity (test passes only with the new implementation of 
    the resistive-CT module in PLUTO 4.1. Crash with PLUTO 4.0).
  - Configuratation #09 employs adaptive mesh refinement as in the original
    PLUTO-Chombo paper (ApJS 2012).
   

  \image html fig1.png "Computed Sweet-Parker scaling for different values of eta with a resolution of 512x256."
  \image html vx_rho_plot__137.png "Density map and magnetic field lines for eta = 2.e-3 at t = 137."

  \authors E. Striani (edoardo.striani@iaps.inaf.it)\n
           A. Mignone (mignone@ph.unito.it)
  \date    March 02, 2017
  
  \b Reference
     - "The PLUTO Code for Adaptive Mesh Computations in
        Astrophysical FLuid Dynamics" Mignone et al., ApJS (2012) 198,7
     - "Geospace Environmental Modeling (GEM) magnetic
        reconnection challenge" Birn et al., JGR (2001) 106, 3715
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x, double y, double x3)
/*
 *********************************************************************** */
{
  double cs2 = 0.5, b0 = 1.0, l, Psi0;
  double Lx, Ly, kx, ky;

  l = g_inputParam[WIDTH];
  v[RHO] = 0.2 + 1.0/(cosh(y/l)*(cosh(y/l)));
  v[PRS] = cs2*v[RHO];  /* v[PRS] = 0.5 in the original PLUTO (2012) paper. */
  v[VX1] = 0.0;
  v[VX2] = 0.0;

#if PHYSICS == MHD || PHYSICS == RMHD
  v[BX1] = b0*tanh(y/l);
  v[BX2] = 0.0;

  Lx = g_domEnd[IDIR] - g_domBeg[IDIR]; kx = CONST_PI/Lx;
  Ly = g_domEnd[JDIR] - g_domBeg[JDIR]; ky = CONST_PI/Ly;

  Psi0    = g_inputParam[PSI0];
  v[BX1] += -Psi0*ky*sin(ky*y)*cos(2.0*kx*x);
  v[BX2] +=  Psi0*2.0*kx*sin(2.0*kx*x)*cos(ky*y);
  v[BX3]  = 0.0;

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = Psi0*cos(ky*y)*cos(2.0*kx*x);
#endif
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
/* *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/* *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* *********************************************************************** */
{
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/* *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/* *********************************************************************** */
{
  return 0.0;
}
#endif


