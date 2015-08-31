/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Diffusion of a linear force-free magnetic field in 
         cylindrical coordinates.

  Solve the diffusion equation for a linear force-free magnetic field 
  in cylindrical coordinates.
  A force-free field satisfies
  \f[
      \vec{J} \times \vec{B} = 0 
     \qquad\Longrightarrow\qquad 
     \vec{J} = \nabla\times\vec{B} = \mu \vec{B}
  \f]
  which means that \f$\vec{J}\f$ and \f$\vec{B}\f$ must be parallel.
  If \f$ \mu \f$ is constant then the solution is given by the Bessel 
  functions:
  \f[
     B_z(r)    = B_0 J_1(r) \,,\qquad
     B_\phi(r) = B_0 J_0(r)
  \f]
  For constant resistivity and zero velocity the induction equation 
  simplifies as follows:
  \f[
   \frac{d\vec{B}}{dt} = - \nabla\times (\eta \vec{J}) 
   = - \nabla\times (\eta \mu \vec{B}) = -\eta \mu^2 \vec{B}
  \f]
  which admits the exact analytical solution 
  \f$\vec{B}(r,t) = \vec{B}(r,0) \exp(-\eta\mu^2 t)\f$ meaning that 
  the field remains force free also at subsequent times.
   
  Note that if pressure and density are initially constant and the
  velocity is also initially zero everywhere, the previous solution
  is also an exact solution of the isothermal MHD equations but not 
  of the adiabatic MHD equations because of the Ohmic dissipation term 
  (magnetic energy transforms into heat).
  However, using a large density makes pressure effects negligible.

  - Configurations #01-03 solve the problem on a cylindrical domain
    away from the origin.
  - Configurations #04-06 extends the integration up to the origin
    and tests the quality of the discretization at r = 0

  \author T. Matsakos\n
          A. Mignone (mignone@ph.unito.it)
  \date   April 13, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  us[RHO] = g_inputParam[DENSITY];
  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = 0.0;
  #if HAVE_ENERGY
   us[PRS] = 1.0;
  #endif
  us[TRC] = 0.0;

  us[BX1] = 0.0;
  us[BX2] = BesselJ0(x1);
  us[BX3] = BesselJ1(x1);
  #ifdef GLM_MHD
   us[PSI_GLM] = 0.0;
  #endif
}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 *********************************************************************** */
{
  int  i, j, k, nv;
  double *r, e_t;
  double t, dt, vin[256];
  static float *j0, *j1; /* store bessel function to avoid overhead */

  t   = g_time;
  r   = grid[IDIR].x;
  e_t = exp(-t*1.0);

/* -------------------------------------------------
   Pre-compute Bessel functions to avoid overhead
   ------------------------------------------------- */

  if (j0 == NULL){
    j0 = ARRAY_1D(NX1_TOT, float);
    j1 = ARRAY_1D(NX1_TOT, float);
    ITOT_LOOP(i) {
      j0[i] = BesselJ0(r[i]);
      j1[i] = BesselJ1(r[i]);
    }
  }

/* ----------------------------------------------------
    Set default boundary values
   ---------------------------------------------------- */

  vin[RHO] = g_inputParam[DENSITY];
  vin[VX1] = vin[VX2] = vin[VX3] = 0.0;
  vin[BX1] = vin[BX2] = vin[BX3] = 0.0;
  #if HAVE_ENERGY
   vin[PRS] = 1.0;
  #endif
  #ifdef GLM_MHD
   vin[PSI_GLM] = 0.0;
  #endif

/* ---------------------------------------------------
    Set boundary conditions
   --------------------------------------------------- */

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        VAR_LOOP(nv) d->Vc[nv][k][j][i] = vin[nv];
        d->Vc[BX2][k][j][i] = j0[i]*e_t;
        d->Vc[BX3][k][j][i] = j1[i]*e_t;
      }
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = j0[i]*e_t; 
      #endif
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        VAR_LOOP(nv)  d->Vc[nv][k][j][i] = vin[nv];
        d->Vc[BX2][k][j][i] = j0[i]*e_t;
        d->Vc[BX3][k][j][i] = j1[i]*e_t;
      }
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = j0[i]*e_t; 
      #endif
    }
  }
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{
}
