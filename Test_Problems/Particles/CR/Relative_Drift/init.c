/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Fluid-CR relative drift test.

  Set the fluid configuration for the Fluid-particles relative drift test,
  as described in section 4.3 of [MVBM18].
  
  The six configuration corresponds to the five cases plotted in Fig4
  using a MUSCL-Hancock scheme plus one case with RK2:
  
  - #01: Nsub = 1, np
  - #02: Nsub = 1, wp1
  - #03: Nsub = 5, np
  - #04: Nsub = 5, wp1
  - #05: Nsub = 4, wp2
  - #06: Nsub = 4, wp2 (RK) 
  
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   April 2, 2018
  \b References: \n
   - [MVBM18] "A PARTICLE MODULE FOR THE PLUTO CODE: I - AN IMPLEMENTATION OF
               THE MHD-PIC EQUATIONS", Mignone etal.ApJS (2018)  [ Sec. 4.3 ]
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 *
 *********************************************************************** */
{
  int np_cell  = RuntimeGet()->Nparticles_cell;
  double rho_p = 1.e-2*g_inputParam[RHO_GAS]/np_cell; /* Particle mass */

  double rho_g = g_inputParam[RHO_GAS];
  double vpx0  = g_inputParam[VPX1];
  double vpy0  = g_inputParam[VPX2];
  double B0    = 2.0*CONST_PI;
  
/* --------------------------------------------------------
    Assign initial condition in a frame of reference where
    total momentum is initially 0.
   -------------------------------------------------------- */
     
  v[RHO] = rho_g;
  v[VX1] = -rho_p*vpx0/rho_g;
  v[VX2] = -rho_p*vpy0/rho_g;
  v[VX3] = 0.0;
  v[PRS] = 1.0;
  v[TRC] = 0.0;

  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = B0;

double L = g_domEnd[IDIR] - g_domBeg[IDIR];
//v[BX3] += 0.1*B0*sin(2.0*CONST_PI*x1/L)*cos(4.0*CONST_PI*x2/L);
  
  #if DIMENSIONS == 3
  v[BX1] = B0*sin(CONST_PI/3.0);
  v[BX3] = B0*cos(CONST_PI/3.0);
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
/*! 
 * Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
  static int first_call = 1;
  int np_cell  = RuntimeGet()->Nparticles_cell;
  double rho_p = 1.e-2*g_inputParam[RHO_GAS]/np_cell; /* Particle mass */
  double Bz0  = 2.0*CONST_PI;
  double rho  = g_inputParam[RHO_GAS];
  double vpx0, vgx0, vgx_ex, vpx_ex, vgx, dvx, dvx_ex;
  double vpy0, vgy0, vgy_ex, vpy_ex, vgy, dvy, dvy_ex;

  double Omega_g, Omega_p, Omega, R, err;
  double v[256];
  double c, s;
  Particle *p;
  FILE *fp;
  
/* --------------------------------------------------------
   0. Set initial conditions and constants
   -------------------------------------------------------- */

  Init(v, 0.0, 0.0, 0.0);
  vgx0 = v[VX1];
  vgy0 = v[VX2];
  vpx0 = g_inputParam[VPX1];
  vpy0 = g_inputParam[VPX2];

//  p = &(d->PHead->p);
  p = Particles_Select(d->PHead, 1);

/* -- Compute exact solution -- */

#if PARTICLES_CR_FEEDBACK == YES
  R  = PARTICLES_CR_E_MC*rho_p;
  R /= PARTICLES_CR_E_MC*rho_p + PARTICLES_CR_E_MC_GAS*rho;
#else
  R  = 0.0;
#endif

  Omega_g = PARTICLES_CR_E_MC_GAS*R*Bz0;
  Omega_p = PARTICLES_CR_E_MC*(1.0 - R)*Bz0;
  Omega   = Omega_g + Omega_p;
  c       = cos(Omega*g_time);
  s       = sin(Omega*g_time);

/* --------------------------------------------------------
   1. Compute exact analytical solution & errors
   -------------------------------------------------------- */

/* -- 1a. Fluid exact solution -- */

  vgx_ex =   Omega_g*(vgx0 - vpx0)*c + Omega_g*(vgy0 - vpy0)*s
           + Omega_g*vpx0 + Omega_p*vgx0;

  vgy_ex = - Omega_g*(vgx0 - vpx0)*s + Omega_g*(vgy0 - vpy0)*c
           + Omega_g*vpy0 + Omega_p*vgy0;

/* -- 1b. Particle exact solution -- */

  vpx_ex =   Omega_p*(vpx0 - vgx0)*c + Omega_p*(vpy0 - vgy0)*s
           + Omega_p*vgx0 + Omega_g*vpx0;

  vpy_ex = - Omega_p*(vpx0 - vgx0)*s + Omega_p*(vpy0 - vgy0)*c
           + Omega_p*vgy0 + Omega_g*vpy0;

  vgx_ex /= Omega;
  vgy_ex /= Omega;
  vpx_ex /= Omega;
  vpy_ex /= Omega;

/* -- 1c. Fluid velocity -- */

  vgx = d->Vc[VX1][KBEG][JBEG][IBEG];
  vgy = d->Vc[VX2][KBEG][JBEG][IBEG];

/* -- 1d. Compute relative velocity dv = vg-vp -- */

  dvx = vgx - p->speed[IDIR];
  dvy = vgy - p->speed[JDIR];

/* -- 1e. Compute exact relative velocity -- */

  dvx_ex = vgx_ex - vpx_ex;
  dvy_ex = vgy_ex - vpy_ex;

/* --------------------------------------------------------
   2. Open file for writing / appending 
   -------------------------------------------------------- */

  if (first_call) {
    fp = fopen("error.dat","w");
    fprintf (fp, "#  R       = %12.6e\n",R);
    fprintf (fp, "#  Omega_g = %12.6e\n",Omega_g);
    fprintf (fp, "#  Omega_p = %12.6e\n",Omega_p);
    fprintf (fp, "#  Omega   = %12.6e  (%12.6e, %12.6e)\n",Omega,
                      R*PARTICLES_CR_E_MC_GAS, (1.0-R)*PARTICLES_CR_E_MC);
    fprintf (fp, "# -------------------------------------------------");
    fprintf (fp, "---------------------------------------------------\n");
  } else {
    fp = fopen("error.dat","a");
  }

/* -- Compute L1-norm error and write it to disk -- */

  err =   fabs(vgx_ex - vgx)            + fabs(vgy_ex - vgy)
        + fabs(vpx_ex - p->speed[IDIR]) + fabs(vpy_ex - p->speed[JDIR]);

  fprintf (fp,"%12.6e  %12.6e  %12.6e\n", g_time, g_dt, err);
/*
  fprintf (fp,"%12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n",
           g_time, g_dt, fabs(p->speed[IDIR]-vp_ex[IDIR]),
                         fabs(p->speed[JDIR]-vp_ex[JDIR]),
                         fabs(p->speed[KDIR]-vp_ex[KDIR]),
                         fabs(vg[IDIR] - vg_ex[IDIR]),
                         fabs(vg[JDIR] - vg_ex[JDIR]),
                         fabs(vg[KDIR] - vg_ex[KDIR]));
*/
  fclose(fp);

/* -- Write a particle trajectory to file -- */

  if (first_call) Particles_WriteTrajectory (p, 'w');
  else            Particles_WriteTrajectory (p, 'a');

  first_call = 0;

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){};
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_END){  /* -- X3_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
