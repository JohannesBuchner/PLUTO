/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Taylor-Couette Flow in 2D cylindrical coordinates.

  This problem considers a fluid rotating between two concentric 
  cylinders situated at \f$ R_{\rm int}\f$ and \f$ R_{\rm ext}\f$
  (fixed by the computational domain in \c pluto.ini). 
  The outer cylinder is not rotating while the inner one rotates with 
  anular velocity \f$\omega\f$.
  Viscous effects are controlled by the Reynolds number defined in 
  \c visc_nu.c as
  \f[
     {\rm Re} = \Omega R_{int}\left(R_{\rm ext} - R_{\rm int}\right) 
                \frac{\rho}{\nu_1}
  \f]
  For Reynolds numbers \f${\rm Re} > {\rm Re}_{\rm critical}\sim 130\f$ 
  vortices are formed, the axial distribution of which is controlled by 
  the wave-number of the initial perturbation \f$\kappa = 2 \pi / \lambda\f$.
  For small Reynolds numbers, viscosity suppresses the vortex formation.
  
  The input parameters for this problem are:

  - <tt>g_inputParam[OMEGA]</tt>: the angular velocity of the inner cylinder.
  - <tt>g_inputParam[REYN]</tt>:  the Reynolds number.

  \author P. Tzeferacos (petros.tzeferacos@ph.unito.it)\n
          A. Mignone (mignone@ph.unito.it)
  \date   March 3, 2017

  \b References
     - "Stability of a viscous liquid contained between two rotating cylinders",
         Taylor G.I.,  1923, Philos. Trans. R. Soc. London, Ser. A, 223, 289
     - "Taylor vortex formation in axial through-flow: Linear and 
        weakly nonlinear analysis", Recktenwald A., L\"ucke M., M\"uller H.W.,  
        1993, Phys. Rev. E, 48, 4444 
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  double rmin = g_domBeg[IDIR];
  double rmax = g_domEnd[IDIR];
  double eta = rmin/rmax;
  double lambda = (rmax-rmin); 
  double A = -g_inputParam[OMEGA]*eta*eta/(1.0 - eta*eta);
  double B =  g_inputParam[OMEGA]*rmin*rmin/(1.0 - eta*eta);
  double SPER = sin(2.*CONST_PI*x2/lambda);
  double CPER = cos(2.*CONST_PI*x2/lambda);
  
  v[RHO] = 1.0; 
  v[PRS] = 10.+ 0.5*(A*A*x1*x1 +4.*A*B*log(x1) - B*B/x1/x1) 
           + 0.01*0.5*(A*A*x1*x1 +4.*A*B*log(x1) - B*B/x1/x1)*CPER; 
  v[VX1] = 0.01*(A*x1 + B/x1)*CPER;
  v[VX2] = 0.01*(A*x1 + B/x1)*SPER;
  v[VX3] = A*x1 + B/x1 + 0.01*(A*x1 + B/x1)*CPER;
  v[TRC] = 1.0;
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
/*! 
 *  Assign user-defined boundary conditions at inner and outer 
 *  radial boundaries. Reflective conditions are applied except for 
 *  the azimuthal velocity which is fixed.
 *
 *********************************************************************** */
{
  int     i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){};
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary / Cylinder @ R = R_INT -- */
    X1_BEG_LOOP(k,j,i){
      d->Vc[RHO][k][j][i] =   d->Vc[RHO][k][j][2*IBEG - i - 1];
      d->Vc[PRS][k][j][i] =   d->Vc[PRS][k][j][2*IBEG - i - 1];
      d->Vc[VX1][k][j][i] = - d->Vc[VX1][k][j][2*IBEG - i - 1];
      d->Vc[VX2][k][j][i] =   d->Vc[VX2][k][j][2*IBEG - i - 1];
      d->Vc[VX3][k][j][i] =   g_inputParam[OMEGA]*x1[i];	    	    
    }
  }

  if (side == X1_END){  /* -- X1_END boundary / Cylinder @ R = R_EXT -- */
    X1_END_LOOP(k,j,i){
      d->Vc[RHO][k][j][i] =   d->Vc[RHO][k][j][2*IEND - i + 1];
      d->Vc[PRS][k][j][i] =   d->Vc[PRS][k][j][2*IEND - i + 1];
      d->Vc[VX1][k][j][i] = - d->Vc[VX1][k][j][2*IEND - i + 1];
      d->Vc[VX2][k][j][i] =   d->Vc[VX2][k][j][2*IEND - i + 1];
      d->Vc[VX3][k][j][i] =   0.0;	  	    
    }
  }
}

