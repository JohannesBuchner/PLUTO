/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Hydrodynamic jet propagation in 2D cylindrical coordinates.

  This problem considers the propagation of a hydrodynamic jet into a
  static uniform medium with constant density and pressure.
  The ambient density, in units of the jet density, is prescribed to be
  \f$ \rho_a = \eta \f$ where \f$\eta\f$ is the ambient/jet density ratio.
  The ambient pressure is \f$ p = 1/\Gamma \f$ when an \c IDEAL EoS is used
  or it is computed from temperature as \f$ p = p(\rho_a, T_a) \f$ for the
   \c PVTE_LAW EoS (here Ta is the ambient temperature).
  These values are set in Init() while the jet inflow is set through 
  a user-defined boundary condition at the lower z-boundary.
  A simple top-hat injection nozzle is used.
  
  The configuration is defined in terms of the following parameters:

  - <tt>g_inputParam[ETA]</tt>:   density ratio between ambient and jet;
  - <tt>g_inputParam[MACH]</tt>:  jet Mach number;
  - <tt>g_inputParam[TJET]</tt>:  jet temperature (only for \c PVTE_LAW EoS).

  defined in \c pluto.ini.
  The reference density and length are given by the jet density and radius
  while the reference velocity is 1 Km/s.
  The actual numerical values are needed only when using the \c PVTE_LAW EoS.
     
  - Configuration #01 uses an \c IDEAL EoS
  - Configurations #02 and #03 set a highly supersonic molecular jet
    evolving with the \c PVTE_LAW EoS.
    The first one adopts the root-finder version while the second one
    adopts the tabulated version.

  \image html hd_jet.jpg "Pressure (left) and density (right) maps for configuration #01 at t=15" 
  \author A. Mignone (mignone@ph.unito.it)
  \date   March 2, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  double Ta = 50.0;  /* Ambient temperature */

  v[RHO] = g_inputParam[ETA];
  v[VX1] = 0.0;
  v[VX2] = 0.0;

  #if EOS == IDEAL
   g_gamma = 5.0/3.0;
   v[PRS]  = 1.0/g_gamma;
  #elif EOS == PVTE_LAW
   v[PRS] = Pressure(v,Ta);
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

/* **************************************************************** */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 **************************************************************** */
{ }
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions in the lower boundary ghost
 *  zones.  The profile is top-hat: 
 *  \f[
 *     V_{ij} = \left\{\begin{array}{ll}
 *     V_{\rm jet} & \quad\mathrm{for}\quad r_i < 1 \\ \noalign{\medskip}
 *     \mathrm{Reflect}(V)  & \quad\mathrm{otherwise}
 *    \end{array}\right.
 *  \f]
 * where \f$ V_{\rm jet} = (\rho,v,p)_{\rm jet} = (1,M,1/\Gamma)\f$ and
 * \c M is the flow Mach number (the unit velocity is the jet sound speed, 
 * so \f$ v = M\f$).
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  double cs, Tj, vj[NVAR];

  x1 = grid->xgc[IDIR];  /* -- array pointer to x1 coordinate -- */
  x2 = grid->xgc[JDIR];  /* -- array pointer to x2 coordinate -- */
  x3 = grid->xgc[KDIR];  /* -- array pointer to x3 coordinate -- */

  vj[RHO] = 1.0;
  #if EOS == IDEAL
   vj[PRS] = 1.0/g_gamma;         /* -- Pressure-matched jet -- */
   vj[VX2] = g_inputParam[MACH];  /* -- Sound speed is one in this case -- */
  #elif EOS == PVTE_LAW
   Tj      = g_inputParam[TJET];
   vj[PRS] = Pressure(vj,Tj);
   cs      = sqrt(vj[PRS]/vj[RHO]); /* -- For simplicity, isothermal cs is used -- */
   vj[VX2] = g_inputParam[MACH]*cs;
  #endif

  if (side == X2_BEG){     /* -- select the boundary side -- */
    BOX_LOOP(box,k,j,i){   /* -- Loop over boundary zones -- */
      if (x1[i] <= 1.0){   /* -- set jet values for r <= 1 -- */
        d->Vc[RHO][k][j][i] = vj[RHO];
        d->Vc[VX1][k][j][i] = 0.0;
        d->Vc[VX2][k][j][i] = vj[VX2];
        d->Vc[PRS][k][j][i] = vj[PRS]; 
      }else{                /* -- reflective boundary for r > 1 --*/
        VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][2*JBEG-j-1][i];
        d->Vc[VX2][k][j][i] *= -1.0;
/*
        d->Vc[RHO][k][j][i] =  d->Vc[RHO][k][2*JBEG - j - 1][i];
        d->Vc[VX1][k][j][i] =  d->Vc[VX1][k][2*JBEG - j - 1][i];
        d->Vc[VX2][k][j][i] = -d->Vc[VX2][k][2*JBEG - j - 1][i];
        d->Vc[PRS][k][j][i] =  d->Vc[PRS][k][2*JBEG - j - 1][i];
*/
      }
    }
  }
}

#if (BODY_FORCE & VECTOR)
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the graviational potential as function of the coordinates.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
