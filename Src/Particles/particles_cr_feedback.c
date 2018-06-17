/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief File containing functions for feedback from CR to the fluid.
  
  \authors A. Mignone (mignone@ph.unito.it)\n

  \b References
    - "MAGNETOHYDRODYNAMIC-PARTICLE-IN-CELL METHOD FOR COUPLING COSMIC RAYS 
       WITH A THERMAL PLASMA: APPLICATION TO NON-RELATIVISTIC SHOCKS"\n
       Bai et al., ApJ (2015) 809, 55

  \date   June 15, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if PARTICLES_TYPE == COSMIC_RAYS

#define CROSS_PRODUCT_X1(a,b)  ( (a)[JDIR]*(b)[KDIR] - (a)[KDIR]*(a)[JDIR] )
#define CROSS_PRODUCT_X2(a,b)  ( (a)[KDIR]*(b)[IDIR] - (a)[IDIR]*(a)[KDIR] )
#define CROSS_PRODUCT_X3(a,b)  ( (a)[IDIR]*(b)[JDIR] - (a)[JDIR]*(a)[IDIR] )

/* ********************************************************************* */
void Particles_CR_ConservativeFeedback(Data_Arr U, Data_Arr Fcr,
                                       double dt, RBox *box)
/*!
 *  Add source term to the right hand side of the conservative MHD
 *  equations.
 *
 * \param [out]  U   An array of conserved quantities, U[k][j][i][nv]
 * \param [in]  Fcr  An array containing the CR feedback, Fcr[4][k][j][i]
 * \param [in]  dt   The time step increment
 * \param [in] box   A pointer to a RBox structure.
 * 
 *********************************************************************** */
{
#if PARTICLES_CR_FEEDBACK == YES
  int    i,j,k;

  BOX_LOOP(box,k,j,i){
    U[k][j][i][MX1] -= dt*Fcr[IDIR][k][j][i];
    U[k][j][i][MX2] -= dt*Fcr[JDIR][k][j][i];
    U[k][j][i][MX3] -= dt*Fcr[KDIR][k][j][i];
    #if HAVE_ENERGY
    U[k][j][i][ENG] -= dt*Fcr[3][k][j][i];
    #endif
  }
#endif  /* PARTICLES_CR_FEEDBACK == YES */
}

/* ********************************************************************* */
void Particles_CR_Flux(const State *state, int beg, int end)
/*!
 *  Add CR flux contribution to MHD fluxes.
 *  On input, it requires \c state->Fcr to be already computed.
 *  On output, it returns \c state->fluxCR.
 * 
 *********************************************************************** */
{
#if PARTICLES_CR_FEEDBACK == YES
  int    nv, i, dir;
  double *v, *B, *Fcr;
  double B2, FcrB, inv_qg;
  double **fluxCR = state->fluxCR;

  DEBUG_FUNC_BEG ("Particles_CR_Flux");

  if (g_dir == IDIR){
    for (i = beg; i <= end; i++) {
      v   = state->v[i];
      B   = v + BX1;
      Fcr = state->Fcr[i];
      
    /* -- Clean parallel component (Fcr.B = 0) -- */
    
      FcrB = DOT_PRODUCT(Fcr, B);
      B2   = DOT_PRODUCT(B,B);
      for (dir = 0; dir < 3; dir++) Fcr[dir] -= FcrB*B[dir]/(B2 + 1.e-16);
      
      NVAR_LOOP(nv) fluxCR[i][nv] = 0.0;

      inv_qg  = 1.0/(PARTICLES_CR_E_MC_GAS*v[RHO]);
      fluxCR[i][BX2] += Fcr[KDIR]*inv_qg;
      fluxCR[i][BX3] -= Fcr[JDIR]*inv_qg;
      #if HAVE_ENERGY
      fluxCR[i][ENG] -= (Fcr[JDIR]*v[BX3] - Fcr[KDIR]*v[BX2])*inv_qg;
      #endif
    }
  }

  if (g_dir == JDIR){
    for (i = beg; i <= end; i++) {
      v   = state->v[i];
      B   = v + BX1;
      Fcr = state->Fcr[i];

    /* -- Clean parallel component (Fcr.B = 0) -- */
   
      FcrB = DOT_PRODUCT(Fcr, B);
      B2   = DOT_PRODUCT(B,B);
      for (dir = 0; dir < 3; dir++) Fcr[dir] -= FcrB*B[dir]/(B2 + 1.e-16);

      NVAR_LOOP(nv) fluxCR[i][nv] = 0.0;

      inv_qg  = 1.0/(PARTICLES_CR_E_MC_GAS*v[RHO]);
      fluxCR[i][BX1] -= Fcr[KDIR]*inv_qg;
      fluxCR[i][BX3] += Fcr[IDIR]*inv_qg;
      #if HAVE_ENERGY
      fluxCR[i][ENG] -= (Fcr[KDIR]*v[BX1] - Fcr[IDIR]*v[BX3])*inv_qg;
      #endif
    }
  }

  if (g_dir == KDIR){
    for (i = beg; i <= end; i++) {
      v   = state->v[i];
      B   = v + BX1;
      Fcr = state->Fcr[i];

    /* -- Clean parallel component (Fcr.B = 0) -- */
    
      FcrB = DOT_PRODUCT(Fcr, B);
      B2   = DOT_PRODUCT(B,B);
      for (dir = 0; dir < 3; dir++) Fcr[dir] -= FcrB*B[dir]/(B2 + 1.e-16);

      NVAR_LOOP(nv) fluxCR[i][nv] = 0.0;

      inv_qg  = 1.0/(PARTICLES_CR_E_MC_GAS*v[RHO]);
      fluxCR[i][BX1] += Fcr[JDIR]*inv_qg;
      fluxCR[i][BX2] -= Fcr[IDIR]*inv_qg;
      #if HAVE_ENERGY
      fluxCR[i][ENG] -= (Fcr[IDIR]*v[BX2] - Fcr[JDIR]*v[BX1])*inv_qg;
      #endif
    }
  }

  DEBUG_FUNC_END ("Particles_CR_Flux");
#endif
}

/* ********************************************************************* */
void Particles_CR_States1DCopy(const Data *d, const Sweep *sweep, int beg, int end)
/*!
 * Load Fcr from 3D data into 1D state structures.
 * This is needed later to compute the flux (Particles_CR_Flux).
 *
 *********************************************************************** */
{
#if PARTICLES_CR_FEEDBACK == YES
  int    i, j, k, dir;
  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);
  double dFcr_lim;
  static double **dFcr;
  
  DEBUG_FUNC_BEG ("Particles_CR_States1DCopy");
  
/* --------------------------------------------------------
   0. Allocate memory.
   -------------------------------------------------------- */
  
  if (dFcr == NULL){
    dFcr = ARRAY_2D(NMAX_POINT, 3, double);  
  }
  
/* --------------------------------------------------------
   1. Load Fcr[0..3] into 1D state, cell centers.
   -------------------------------------------------------- */

  if (g_dir == IDIR){

    j = g_j; k = g_k;
    for (i = 0; i < NX1_TOT; i++){
      for (dir = 0; dir < 4; dir++) stateC->Fcr[i][dir] = d->Fcr[dir][k][j][i];
    }

  }else if (g_dir == JDIR){

    i = g_i; k = g_k;
    for (j = 0; j < NX2_TOT; j++){
      for (dir = 0; dir < 4; dir++) stateC->Fcr[j][dir] = d->Fcr[dir][k][j][i];
    }

  } else if (g_dir == KDIR){

    i = g_i; j = g_j;
    for (k = 0; k < NX3_TOT; k++){
      for (dir = 0; dir < 4; dir++) stateC->Fcr[k][dir] = d->Fcr[dir][k][j][i];
    }
  }

/* --------------------------------------------------------
   2. Reconstruct Fcr[0..2] at cell interfaces.
   -------------------------------------------------------- */

  for (i = beg; i <= end+1; i++){
    for (dir = 0; dir < 3; dir++){
      dFcr[i][dir] = stateC->Fcr[i][dir] - stateC->Fcr[i-1][dir];
    }
  }
  
  for (i = beg; i <= end; i++){
    for (dir = 0; dir < 3; dir++){
      dFcr_lim = MC(dFcr[i+1][dir], dFcr[i][dir]);
      stateL->Fcr[i][dir]   = stateC->Fcr[i][dir] + 0.5*dFcr_lim;
      stateR->Fcr[i-1][dir] = stateC->Fcr[i][dir] - 0.5*dFcr_lim;
//stateL->Fcr[i][dir]   = 0.5*(stateC->Fcr[i][dir] + stateC->Fcr[i+1][dir]);
//stateR->Fcr[i-1][dir] = 0.5*(stateC->Fcr[i][dir] + stateC->Fcr[i-1][dir]);
    }  
  }
  
  DEBUG_FUNC_END ("Particles_CR_States1DCopy");
#endif  /* PARTICLES_CR_FEEDBACK == YES */
}

/* ********************************************************************* */
void Particles_CR_StatesSource(const Sweep *sweep, double dt,
                                 int beg, int end, Grid *grid)
/*!
 *  Add CR flux-difference and source term contributions to L/R
 *  primitive states during the predictor step of CTU algorithm.
 *  Note that, although the update should be done in conservative variables,
 *  we choose to directly update primitive variables by performing, in place,
 *  a conservative-to-primitive recovery so that:
 *  \f[
 *     \left\{\begin{array}{lcl}
 *       \vec{m}_\pm &\leftarrow& \DS
 *       \vec{m}_\pm - \frac{\Delta t}{2}\vec{F}_{CR}
 *       \\ \noalign{\medskip}
 *       \vec{B}_\pm &\leftarrow& \DS
 *       \vec{B}_\pm + \frac{\Delta t}{2}\nabla_x\times
 *                     \left(\frac{\vec{F}_{CR}}{q_g}\right)
 *       \\ \noalign{\medskip}
 *       {\cal E}_\pm &\leftarrow& \DS
 *       {\cal E}_\pm - \frac{\Delta t}{2}\nabla_x\cdot
 *       \left(\frac{\vec{F}_{CR}\times\vec{B}}{q_g} + \vec{F}_{CR}\cdot\vec{v}\right)
 *     \end{array}\right.
 *     \quad\Longrightarrow\quad
 *     \left\{\begin{array}{lcl}
 *       \vec{v}_1  &=& \vec{v}_0 + \Delta\vec{v}
 *       \\ \noalign{\medskip}
 *       \vec{B}_1  &=& \vec{B}_0 + \Delta\vec{B}
 *       \\ \noalign{\medskip}
 *       (\rho e)_1 &=& \DS (\rho e)_0 + \frac{\rho}{2}(\vec{v}_0^2 - \vec{v}_1^2)
 *       + \frac{1}{2}\left(\vec{B}_0^2-\vec{B}_1^2\right) + \Delta {\cal E};
 *     \end{array}\right.
 *  \f]
 *
 *  Note that curl and divergence operators are discretized using finite
 *  differences of \c fluxCR computed in Particles_CR_Flux().
 *
 *  This step is done only for the primitive Hancock / Char. Trac. schemes
 *  since these are carried out in primitive variables.
 *  This step is not needed for the conservative Hancock scheme since
 *  these terms are included during the predictor. 
 *
 *  For source terms not involving derivatives, we include *ALL* of the
 *  components here since no other partcles-related source term will be
 *  added during the transverse predictors.
 *  For source terms involving derivatives, we include only the normal
 *  components since the transverse predictor will add the other
 *  components through Riemann solver.
 *
 *********************************************************************** */
{
#if PARTICLES_CR_FEEDBACK == YES
#if   (TIME_STEPPING == HANCOCK && PRIMITIVE_HANCOCK == YES) \
   || (TIME_STEPPING == CHARACTERISTIC_TRACING)
  int   i, nv;

  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);

  double scrh, dtdx, src[NVAR];
  double *v, *vp, *vm, *Fcr;
  double dv2_p, dv2_m, dB2_p, dB2_m;

  double **fluxCR_p = stateL->fluxCR;
  double **fluxCR_m = stateR->fluxCR - 1;

  intList var_list = {6, VX1, VX2, VX3, BX1, BX2, BX3};


  for (i = beg; i <= end; i++){
    Fcr = stateC->Fcr[i];
    v   = stateC->v[i];
    vp  = stateL->v[i];
    vm  = stateR->v[i-1];

    scrh = dt/v[RHO];
    
    src[VX1] = -scrh*Fcr[IDIR]; 
    src[VX2] = -scrh*Fcr[JDIR]; 
    src[VX3] = -scrh*Fcr[KDIR]; 

    dtdx     = dt/grid->dx[g_dir][i];
    src[BX1] = -dtdx*(fluxCR_p[i][BX1] - fluxCR_m[i][BX1]);
    src[BX2] = -dtdx*(fluxCR_p[i][BX2] - fluxCR_m[i][BX2]);
    src[BX3] = -dtdx*(fluxCR_p[i][BX3] - fluxCR_m[i][BX3]);

    dv2_p = vp[VX1]*vp[VX1] + vp[VX2]*vp[VX2] + vp[VX3]*vp[VX3];
    dv2_m = vm[VX1]*vm[VX1] + vm[VX2]*vm[VX2] + vm[VX3]*vm[VX3];

    dB2_p = vp[BX1]*vp[BX1] + vp[BX2]*vp[BX2] + vp[BX3]*vp[BX3];
    dB2_m = vm[BX1]*vm[BX1] + vm[BX2]*vm[BX2] + vm[BX3]*vm[BX3];

    FOR_EACH(nv, &var_list){
      vp[nv] += src[nv];
      vm[nv] += src[nv];
    }

    dv2_p -= vp[VX1]*vp[VX1] + vp[VX2]*vp[VX2] + vp[VX3]*vp[VX3];
    dv2_m -= vm[VX1]*vm[VX1] + vm[VX2]*vm[VX2] + vm[VX3]*vm[VX3];

    dB2_p -= vp[BX1]*vp[BX1] + vp[BX2]*vp[BX2] + vp[BX3]*vp[BX3];
    dB2_m -= vm[BX1]*vm[BX1] + vm[BX2]*vm[BX2] + vm[BX3]*vm[BX3];

    src[ENG] = -dtdx*(fluxCR_p[i][ENG] - fluxCR_m[i][ENG]) - dt*Fcr[3];

    vp[PRS] += (g_gamma - 1.0)*(0.5*vp[RHO]*dv2_p + 0.5*dB2_p + src[ENG]);
    vm[PRS] += (g_gamma - 1.0)*(0.5*vm[RHO]*dv2_m + 0.5*dB2_m + src[ENG]);
  }

#endif /* TIME_STEPPING == HANCOCK */
#endif /* PARTICLES_FEEDBACK == YES */
}

/* ********************************************************************* */
void Particles_CR_StatesSourceOld(const Sweep *sweep, double dt,
                                int beg, int end, Grid *grid)
/*!
 *  Add CR flux-difference and source term contributions to L/R
 *  conservative states during the predictor step of CTU algorithm.
 *
 *  This step is done only for the primitive Hancock / Char. Trac. schemes
 *  since these are carried out in primitive variables.
 *  This step is not needed for the conservative Hancock scheme since
 *  these terms are included during the predictor. 
 *
 *  Note that, unlike gravity, we include *ALL* components of the
 *  force here since no other partcles-related source term will be
 *  added during the transverse predictors.
 *
 *********************************************************************** */
{
#if PARTICLES_CR_FEEDBACK == YES
#if   (TIME_STEPPING == HANCOCK && PRIMITIVE_HANCOCK == YES) \
   || (TIME_STEPPING == CHARACTERISTIC_TRACING)
  int   i, nv;
  static unsigned char *flag;

  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);

  double  **up = stateL->u;
  double  **um = stateR->u-1;

  double **fluxCR   = stateC->fluxCR;
  double **fluxCR_p = stateL->fluxCR;
  double **fluxCR_m = stateR->fluxCR - 1;
  double dtdx, du, dup, dum;
  #if HAVE_ENERGY
  intList var_list = {4, BX1, BX2, BX3, ENG};
  #else
  intList var_list = {3, BX1, BX2, BX3};
  #endif

  DEBUG_FUNC_BEG("Particles_CR_StatesSource");
  
/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */

  if (flag == NULL) {
    flag    = ARRAY_1D(NMAX_POINT, unsigned char);
  }
  
/* --------------------------------------------------------
   2. Loop over one row of zones and modify
      conservative L/R 1D states
   -------------------------------------------------------- */
  
  for (i = beg; i <= end; i++){

  /* --------------------------------------------
     2b. Add flux-difference contribution to
         magnetic field and energy equations
     -------------------------------------------- */

    dtdx = dt/grid->dx[g_dir][i];
    FOR_EACH(nv, &var_list){
      du = dtdx*(fluxCR_p[i][nv] - fluxCR_m[i][nv]);
      um[i][nv] -= du;
      up[i][nv] -= du;
    }

  /* --------------------------------------------
     2c. Add source term contribution to
         momentum and energy
     -------------------------------------------- */

    um[i][MX1] -= stateC->Fcr[i][IDIR]*dt;  
    up[i][MX1] -= stateC->Fcr[i][IDIR]*dt;  
    
    um[i][MX2] -= stateC->Fcr[i][JDIR]*dt;  
    up[i][MX2] -= stateC->Fcr[i][JDIR]*dt;  
    
    um[i][MX3] -= stateC->Fcr[i][KDIR]*dt;
    up[i][MX3] -= stateC->Fcr[i][KDIR]*dt;  
    
    #if HAVE_ENERGY
    um[i][ENG] -= stateC->Fcr[i][3]*dt;
    up[i][ENG] -= stateC->Fcr[i][3]*dt;  
    #endif
  }

/* --------------------------------------------------------
   3. Convert conservative to primitive states
   -------------------------------------------------------- */

  ConsToPrim (stateL->u, stateL->v, beg, end, flag);
  ConsToPrim (stateR->u, stateR->v, beg-1, end-1, flag);

  DEBUG_FUNC_END("Particles_CR_StatesSource");
#endif /* TIME_STEPPING == HANCOCK */
#endif /* PARTICLES_FEEDBACK == YES */
}


#endif  /* PARTICLES_TYPE == COSMIC_RAYS */

