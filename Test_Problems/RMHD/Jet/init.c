/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Relativistic magnetized jet in axisymmetric coordinates.

  The jet problem is set up in axisymmetric cylindrical coordinates
  \f$ (r,z) \f$ as described in section 4.2.3 of 
  Mignone, Ugliano \& Bodo (2009, [MUB09] hereafter).
  We label ambient and jet values with the suffix "a" and "j", respectively.
  The ambient medium is at rest and is characterized by uniform density
  and pressure, \f$ \rho_a \f$ and \f$ p_a \f$.
  The beam enters from the lower-z boundary from a circular nozzle of
  radius 1, carries a constant poloidal field \f$ B_z \f$
  and a (radially-varying) toroidal component \f$ B_\phi(R) \f$.
  Flow variables are prescribed as
  \f[ \left\{\begin{array}{lcl}
     \rho(R)   &=& \rho_j  \\ \noalign{\medskip}
     v_z(R)    &=& v_j     \\ \noalign{\medskip}
     B_\phi(R) &=& \DS \left\{\begin{array}{ll}
             -\gamma_j b_m R/a & \quad{\rm for}\quad r<a \\ \noalign{\medskip}
             -\gamma_j b_m a/R & \quad{\rm for}\quad r>a
             \end{array}\right. \\ \noalign{\medskip}
     B_z(R)    &=& B_{z0}\quad\mathrm{(const)}
  \end{array}\right. \f]
  Here \c a is the magnetization radius while \f$ v_r = B_r = 0 \f$.
  The pressure distribution is found by solving the radial momentum
  balance between thermal, centrifugal and magnetic forces. 
  Neglecting rotation and assuming \c Bz to be constant the solution 
  is given by
  \f[
    p(R) = p_j + b_m^2\left[1-\min\left(\frac{r^2}{a^2},1\right)\right]
  \f]
  Here \f$p_j=p_a\f$ is the jet/ambient pressure at \c r=1 and its
  value depends on the Mach number, see [MUB09].

  The parameters controlling the problem is
  -# <tt>g_inputParam[MACH]</tt>: the jet Mach number
     \f$ M = v_j/cs \f$ where \f$c_s=\sqrt{\Gamma p_j/(\rho h)}\f$ is 
     the sound speed and \f$\rho h = \rho + \Gamma p_j/(\Gamma-1)\f$.
     This is used to recover \f$p_j\f$;
  -# <tt>g_inputParam[LORENTZ]</tt>:  the jet Lorentz factor;
  -# <tt>g_inpurParam[RHOJ]</tt>: the jet density;
  -# <tt>g_inpurParam[SIGMA_POL]</tt>: magnetization strength for poloidal 
     magnetic field component (see [MUB09], Eq [65]);
  -# <tt>g_inpurParam[SIGMA_TOR]</tt>: magnetization strength for toroidal
     magnetic field component (see [MUB09], Eq [65]);

  The different configurations are:

  - configurations #01-02 employ a purely poloidal field and are similar 
    to section 4.3.3 of Mignone \& Bodo (2006);
  - configuration #03

  The following image show the (log) density map at the end of
  simulation for setup #01.

  \image html rmhd_jet.01.jpg "Density map at the end of the computation for configuration #01"

  \b References
     - "A five-wave Harten-Lax-van Leer Riemann solver for relativistic
        magnetohydrodynamics", Mignone, Ugliano \& Bodo, MNARS (200) 393, 1141
     - "An HLLC Rieman solver for relativistic flows - II. Magnetohydrodynamics"
        Mignone \& Bodo, MNRAS (2006) 368, 1040   

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 17, 2014

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static void   GetJetValues (double x1, double x2, double x3, double *vj);
static double pj;

#if EOS == TAUB
 static double g_gamma;
#endif

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  int    nv;
  double r, z, p0, vjet[256], vamb[256];

  g_gamma = 5./3.;
  GetJetValues (x1, x2, x3, vjet);

/* -- set ambient values -- */

  v[RHO] = g_inputParam[RHOA];
  EXPAND(v[VX1] = 0.0;   ,
         v[VX2] = 0.0;   ,
         v[VX3] = 0.0;)
  v[PRS] = pj;
  #if PHYSICS == RMHD
   EXPAND(v[BX1] = 0.0;        ,
          v[BX2] = vjet[BX2];  ,
          v[BX3] = 0.0;)
   v[AX1] = 0.0;
   v[AX2] = 0.0;
   v[AX3] = vjet[AX3];
  #endif
  #if NTRACER > 0
   v[TRC] = 0.0;
  #endif
  #ifdef PSI_GLM 
   v[PSI_GLM] = 0.0;
  #endif

  g_smallPressure = 1.e-5;
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
  int   i, j, k, nv;
  double *x1, *x2, *x3;
  double ***bxs, ***bys, ***bzs;
  double prof, vjet[256], vout[NVAR];

  #ifdef STAGGERED_MHD
   D_EXPAND(bxs = d->Vs[BX1s];  , 
            bys = d->Vs[BX2s];  , 
            bzs = d->Vs[BX3s];)
  #endif

  x1  = grid[IDIR].xgc;  
  x2  = grid[JDIR].xgc;  
  x3  = grid[KDIR].xgc;

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER){    /* -- cell-centered boundary conditions -- */
      BOX_LOOP(box,k,j,i){
        GetJetValues (x1[i], x2[j], x3[k], vjet); /* Jet Values */
        VAR_LOOP(nv) vout[nv] = d->Vc[nv][k][2*JBEG - j - 1][i]; /* Ambient */

        vout[VX2] *= -1.0; 
        EXPAND(vout[BX1] *= -1.0;  ,  
                                ;  , 
               vout[BX3] *= -1.0;)
        #ifdef PSI_GLM
         vjet[PSI_GLM]  = d->Vc[PSI_GLM][k][JBEG][i] - d->Vc[BX2][k][JBEG][i];
         vout[PSI_GLM] *= -1.0;
        #endif

        prof = (fabs(x1[i]) <= 1.0); 
        VAR_LOOP(nv) d->Vc[nv][k][j][i] = vout[nv] - (vout[nv] - vjet[nv])*prof;
      }

    }else if (box->vpos == X1FACE){  /* -- staggered fields -- */
      #ifdef STAGGERED_MHD
       x1 = grid[IDIR].xr;
       vjet[BX1] = 0.0; 
       BOX_LOOP(box,k,j,i){
          vout[BX1] = -bxs[k][2*JBEG - j - 1][i];
          prof = (fabs(x1[i]) <= 1.0);
          bxs[k][j][i] = vout[BX1] - (vout[BX1] - vjet[BX1])*prof;
       }
      #endif
    }
  }
}

/* **************************************************************** */
void GetJetValues (double x1, double x2, double x3, double *vj)
/*
 *
 *
 *
 * 
 ****************************************************************** */
{
  static int  first_call = 1;
  double bm, b2_av, bphi, lor, a = 0.5;
  double r, x, scrh, sig_z, sig_phi;

  r   = x1;
  lor = g_inputParam[LORENTZ];
  sig_z   = g_inputParam[SIGMA_POL];
  sig_phi = g_inputParam[SIGMA_TOR];
  if (fabs(r) < 1.e-9) r = 1.e-9;

  x = r/a;
  vj[RHO] = g_inputParam[RHOJ];

  EXPAND(vj[VX1] = 0.0;                        ,
         vj[VX2] = sqrt(1.0 - 1.0/(lor*lor));  , /* 3-vel */
         vj[VX3] = 0.0;)

  scrh = g_inputParam[MACH]/vj[VX2];
  pj   = vj[RHO]/(g_gamma*scrh*scrh - g_gamma/(g_gamma - 1.0));

  bm  = 4.0*pj*sig_phi;
  bm /= a*a*(1.0 - 4.0*log(a) - 2.0*sig_phi);
  bm  = sqrt(bm);

  scrh    = MIN(x*x, 1.0);
  vj[PRS] = pj + bm*bm*(1.0 - scrh);

  bphi  = bm*(fabs(x) < 1.0 ? x: 1.0/x);

  #if GEOMETRY == CYLINDRICAL
   EXPAND(vj[iBR]   = 0.0;                                 ,
          vj[iBZ]   = sqrt(sig_z*(bm*bm*a*a + 2.0*pj));  ,
          vj[iBPHI] = lor*bphi;)
   vj[AX1] = 0.0;
   vj[AX2] = 0.0;
   vj[AX3] = 0.5*r*vj[iBZ];
  #endif

/* -- set tracer value -- */

  #if NTRACER > 0
   vj[TRC] = 1.0;
  #endif
  #ifdef PSI_GLM 
   vj[PSI_GLM] = 0.0;
  #endif

}
