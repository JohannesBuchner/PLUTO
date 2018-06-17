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
     B_z(r)    = B_0 J_1(\mu r) \,,\qquad
     B_\phi(r) = B_0 J_0(\mu r)
  \f]
  with vector potential given by
  \f[
     A_z(r)    = - \int J_1(\mu r)\,dr = \frac{1}{\mu}J_0(\mu r) \,,\qquad
     A_\phi(r) = - \frac{1}{r}\int rJ_0(\mu r)\,dr
               =  \frac{1}{\mu}J_0(\mu r)
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
  
  In the test we set \f$\mu = \eta = 1\f$.

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
  - Configurations #08-09 are done in spherical coordinates.

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 20, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  double Br, Bz, Bphi;
  double Ar, Az, Aphi;
#if GEOMETRY == CYLINDRICAL
  double r = x1;
#elif GEOMETRY == SPHERICAL
  double r = x1*sin(x2);
#endif

  v[RHO] = g_inputParam[DENSITY];
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
#if HAVE_ENERGY
  v[PRS] = 1.0;
 #endif
  v[TRC] = 0.0;

  Br   = 0.0;
  Bz   = BesselJ0(r);
  Bphi = BesselJ1(r);

  Ar   = 0.0;
  Az   = BesselJ0(r);
  Aphi = BesselJ1(r);

#if GEOMETRY == CYLINDRICAL
  v[BX1] = Br;
  v[BX2] = Bz;
  v[BX3] = Bphi;

  v[AX1] = 0.0;
  v[AX2] = Az;
  v[AX3] = Aphi;
#elif GEOMETRY == SPHERICAL
  v[BX1] = Br*sin(x2) + Bz*cos(x2);
  v[BX2] = Br*cos(x2) - Bz*sin(x2);
  v[BX3] = Bphi;

  v[AX1] = Ar*sin(x2) + Az*cos(x2);
  v[AX2] = Ar*cos(x2) - Az*sin(x2);
  v[AX3] = Aphi;
#endif

#ifdef GLM_MHD
  v[PSI_GLM] = 0.0;
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
{
  int  i, j, k, nv;
  double *r, e_t;
  double t, vin[256];
  double *x1 = grid->x[IDIR], *x1r = grid->xr[IDIR];
  double *x2 = grid->x[JDIR], *x2r = grid->xr[JDIR];
  double *x3 = grid->x[KDIR], *x3r = grid->xr[KDIR];

  t   = g_time;
#if RESISTIVITY != NO  
  e_t = exp(-t*1.0);
#else
  e_t = 1.0;  /* Means no decay */
#endif

/* ---------------------------------------------------
    Set boundary conditions
   --------------------------------------------------- */

  if (side == X1_BEG || side == X1_END){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        Init (vin, x1[i], x2[j], x3[k]);
        NVAR_LOOP(nv) d->Vc[nv][k][j][i] = vin[nv];

        d->Vc[BX1][k][j][i] *= e_t;
        d->Vc[BX2][k][j][i] *= e_t;
        d->Vc[BX3][k][j][i] *= e_t;
      }
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i) {
        Init (vin, x1[i], x2r[j], x3[k]);
        d->Vs[BX2s][k][j][i] = vin[BX2]*e_t;
      }
      #endif
    }
  }

}
