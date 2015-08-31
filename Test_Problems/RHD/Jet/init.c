/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Relativistic jet propagation.

  This problem sets initial and boundary conditions for a jet
  propagating into a uniform medium with constant density and pressure.
  The computation can be carried using \c CYLINDRICAL or \c SPHERICAL
  coordinates.
  In the first case, the jet enters at the lower z-boundary with speed
  \f$ v_z = \beta\f$ while in the second case, a conical beam with a
  small aperture (\f$\theta = 5^\circ\f$) is injected with the same
  speed from the lower radial boundary.
  At the border of the nozzle, jet values are smoothly joined with
  ambient values using the Profile() function
  (in the current setting the profile is a sharp transition).
  The jet is pressure-matched so that the beam and ambient pressure
  coincide.
 
  The configuration is defined in terms of the following parameters:

  -# <tt>g_inputParam[BETA]</tt>:     the jet velocity;
  -# <tt>g_inputParam[RHO_IN]</tt>:   the jet density;
  -# <tt>g_inputParam[RHO_OUT]</tt>:  the ambient density.
  -# <tt>g_inputParam[PRESS_IN]</tt>: the jet pressure (also equal to ambient
                                     pressure)

  defined in \c pluto.ini.
  The \c TAUB equation of state is used.

  - Configurations #01 and #02 use \c CYLINDRICAL coordinates;
  - Configuration #03 employs \c SPHERICAL coordinates (see snapshot below)

  \image html rhd_jet.03.jpg "Density (log) for configuration #03 at t=200"

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 18, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static double Profile(double r, int nv);
static void GetJetValues (double *vjet);

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double scrh;
  #if EOS == IDEAL
   g_gamma = 5.0/3.0;
  #endif
  
  v[RHO] = g_inputParam[RHO_OUT];
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[PRS] = g_inputParam[PRESS_IN];

  g_smallPressure = v[PRS]/500.0;
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
  int   nv, i, j, k;
  double  r, vjet[NVAR], vout[NVAR];
 
  #if GEOMETRY == SPHERICAL 
  if (side == X1_BEG){
    GetJetValues(vjet);
    X1_BEG_LOOP(k,j,i){
      VAR_LOOP(nv) vout[nv] = d->Vc[nv][k][j][2*IBEG-i-1];
      vout[VX2] *= -1.0;

      r = grid[JDIR].x[j];
      VAR_LOOP(nv)
        d->Vc[nv][k][j][i] = vout[nv] + (vjet[nv] - vout[nv])*Profile(r,nv);
    }
  }
  #endif

  #if GEOMETRY == CYLINDRICAL || GEOMETRY == CARTESIAN
  if (side == X2_BEG){
    GetJetValues(vjet);
    X2_BEG_LOOP(k,j,i){
      VAR_LOOP(nv) vout[nv] = d->Vc[nv][k][2*JBEG-j-1][i];
      vout[VX2] *= -1.0;

      r = grid[IDIR].x[i];
      for (nv = 0; nv < NVAR; nv++) 
        d->Vc[nv][k][j][i] = vout[nv] + (vjet[nv] - vout[nv])*Profile(r,nv);
    }
  }
  #endif
}

/* ********************************************************************* */
void GetJetValues (double *vjet)
/*
 *
 *
 *********************************************************************** */
{
  vjet[RHO] = g_inputParam[RHO_IN];
  #if GEOMETRY == CYLINDRICAL || GEOMETRY == CARTESIAN
   vjet[VX1] = 0.0;
   vjet[VX2] = g_inputParam[BETA];
  #elif GEOMETRY == SPHERICAL
   vjet[VX1] = g_inputParam[BETA];
   vjet[VX2] = 0.0;
  #endif
  vjet[PRS] = g_inputParam[PRESS_IN];
}
 
/* ********************************************************************* */
double Profile(double r, int nv)
/* 
 *
 *
 *********************************************************************** */
{
  int xn = 14;
  double r0 = 1.0;

  if (nv == RHO) r0 = 1.1;

  #if GEOMETRY == SPHERICAL
   r0 = 5.0/180.0*CONST_PI;
  #endif
  return 1.0/cosh(pow(r/r0,xn));
}
