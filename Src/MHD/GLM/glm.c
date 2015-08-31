/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief GLM module implementation.

  Contains functions for the GLM module.

  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date   Aug 16, 2012

  \b Reference
     - "Hyperbolic Divergence Cleaning for the MHD equations" 
        Dedner et al., JCP (2002) 175, 645 
     - "A second-order unsplit Godunov scheme for cell-centered MHD:
        the CTU-GLM scheme" Mignone \& Tzeferacos, JCP (2010) 229, 2117
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

double glm_ch = -1.0;

/* ********************************************************************* */
void GLM_Solve (const State_1D *state, double **VL, double **VR,
                int beg, int end, Grid *grid)
/*!
 * Solve the 2x2 linear hyperbolic GLM-MHD system given by the divergence 
 * cleaning approach.
 * Build new states VL and VR for Riemann problem.
 * We use Eq. (42) of Dedner et al (2002)
 *
 * \param [in,out] state pointer to a State_1D structure
 * \param [out]    VL    left-interface state to be passed to the
 *                       Riemann solver
 * \param [out]    VR    right-interface state to be passed to the
 *                       Riemann solver
 * \param [in]    beg    starting index of computation
 * \param [in]    end    final index of computation
 * \param [in]    grid   pointer to array of Grid structures
 *
 * The purpose of this function is two-fold:
 *
 * 1. assign a unique value to the normal component of magnetic field 
 *     and to te scalar psi \e before the actual Riemann solver is called.
 * 2. compute GLM fluxes in Bn and psi.
 *
 * The following MAPLE script has been used
 * \code
 *  restart;
 *  with(linalg);
 *  A := matrix(2,2,[ 0, 1, c^2, 0]);
 *
 *  R := matrix(2,2,[ 1, 1, c, -c]);
 *  L := inverse(R);
 *  E := matrix(2,2,[c,0,0,-c]);
 *  multiply(R,multiply(E,L));
 * \endcode
 *
 *********************************************************************** */
{
  int    i, nv, nflag=0;
  double Bm, psim;
  double dB, dpsi;
  double **v;
  static double *wp, *wm, *Bxp, *Bxm;

  if (wp == NULL){
    wp = ARRAY_1D(NMAX_POINT, double);
    wm = ARRAY_1D(NMAX_POINT, double);
    Bxp = ARRAY_1D(NMAX_POINT, double);
    Bxm = ARRAY_1D(NMAX_POINT, double);
  }

  #ifdef CTU 
   #if PARABOLIC_FLUX != EXPLICIT
    if (g_intStage == 1) nflag = 1;
   #endif

  /* ------------------------------------------------------
      The nflag = 1 is not necessary but it improves the 
      results when dimensionally unsplit CTU is used. 
      In practice, it redefines the input normal field 
      Bn to be at time level n rather than n+1/2.
     ------------------------------------------------------- */

   if (nflag){
     v = state->v;
     for (i = beg-1; i <= end+1; i++){
       dB   = v[i+1][BXn] - v[i][BXn];
       dpsi = v[i+1][PSI_GLM] - v[i][PSI_GLM];
       wm[i] =  0.5*(dB - dpsi/glm_ch);
       wp[i] =  0.5*(dB + dpsi/glm_ch);
     }
     for (i = beg; i <= end+1; i++){
       dB  = MC(wm[i], wm[i-1]);
       dB += MC(wp[i], wp[i-1]);
       Bxp[i] = v[i][BXn] + 0.5*dB;
       Bxm[i] = v[i][BXn] - 0.5*dB;
     }
   }  
  #endif

/* -------------------------------------------------
    Solve the Riemann problem for the 2x2 linear 
    GLM system. Use the solution to assign a single 
    value to both left and right (Bn,psi).
   ------------------------------------------------- */

  for (i = beg; i <= end; i++){
    for (nv = NVAR; nv--;   ){
      VL[i][nv] = state->vL[i][nv];
      VR[i][nv] = state->vR[i][nv];
    }

    #ifdef CTU
     if (nflag){
       VL[i][BXn] = Bxp[i];
       VR[i][BXn] = Bxm[i+1];
     }
    #endif

    dB   = VR[i][BXn]     - VL[i][BXn];
    dpsi = VR[i][PSI_GLM] - VL[i][PSI_GLM];

    Bm   = 0.5*(VL[i][BXn]     + VR[i][BXn])     - 0.5*dpsi/glm_ch;
    psim = 0.5*(VL[i][PSI_GLM] + VR[i][PSI_GLM]) - 0.5*glm_ch*dB;

    state->bn[i] = VL[i][BXn] = VR[i][BXn] = Bm;
    VL[i][PSI_GLM] = VR[i][PSI_GLM] = psim;
  }

  #if COMPUTE_DIVB == YES
   if (g_intStage == 1) GLM_ComputeDivB(state, grid);  
  #endif

}

/* ********************************************************************* */
void GLM_Source (const Data_Arr Q, double dt, Grid *grid)
/*!
 * Include the parabolic source term of the Lagrangian multiplier
 * equation in a split fashion for the mixed GLM formulation. 
 * Ref. Mignone & Tzeferacos, JCP (2010) 229, 2117, Equation (27).
 *   
 *
 *********************************************************************** */
{
  int    i,j,k;
  double cr, cp, scrh;
  double dx, dy, dz, dtdx;

  #ifdef CHOMBO
   dtdx = dt;
  #else
   dx   = grid[IDIR].dx[IBEG];
   dtdx = dt/dx;
  #endif

  cp   = sqrt(dx*glm_ch/GLM_ALPHA);
  scrh = dtdx*glm_ch*GLM_ALPHA; /* CFL*g_inputParam[ALPHA]; */

  scrh = exp(-scrh);
  DOM_LOOP(k,j,i) Q[PSI_GLM][k][j][i] *= scrh;
  return;
}

/* ********************************************************************* */
void GLM_ExtendedSource (const State_1D *state, double dt,
                         int beg, int end, Grid *grid)
/*!
 * Add source terms to the right hand side of the conservative equations, 
 * momentum and energy equations only. 
 * This yields the extended GLM equations given by Eq. (24a)--(24c) in 
 *
 * "Hyperbolic Divergence cleaning for the MHD Equations"
 * Dedner et al. (2002), JcP, 175, 645
 *
 *********************************************************************** */
{
  int    i;
  double bx, by, bz;
  double ch2, r, s;
  double *Ar, *Ath, **bgf;
  double *rhs, *v, *Bm, *pm;
  static double *divB, *dpsi, *Bp, *pp;
  Grid   *GG;
/* 
  #if PHYSICS == RMHD
   print1 ("! EGLM not working for the RMHD module\n");
   QUIT_PLUTO(1);
  #endif
*/
  if (divB == NULL){
    divB = ARRAY_1D(NMAX_POINT, double);
    dpsi = ARRAY_1D(NMAX_POINT, double);
    Bp   = ARRAY_1D(NMAX_POINT, double);
    pp   = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------- ---------------------
    compute magnetic field normal component 
    interface value by arithmetic averaging
   -------------------------------------------- */

  ch2 = glm_ch*glm_ch;
  for (i = beg - 1; i <= end; i++) {
    Bp[i] = state->flux[i][PSI_GLM]/ch2;
    pp[i] = state->flux[i][BXn];
  }
  Bm  = Bp - 1;
  pm  = pp - 1;
  GG  = grid + g_dir;
  
/* --------------------------------------------
    Compute div.B contribution from the normal 
    direction (1) in different geometries 
   -------------------------------------------- */
  
  #if GEOMETRY == CARTESIAN

   for (i = beg; i <= end; i++) {
     divB[i] = (Bp[i] - Bm[i])/GG->dx[i];
     dpsi[i] = (pp[i] - pm[i])/GG->dx[i];
   }

  #elif GEOMETRY == CYLINDRICAL

   if (g_dir == IDIR){   /* -- r -- */
     Ar  = grid[IDIR].A;
     for (i = beg; i <= end; i++) {
       divB[i] = (Bp[i]*Ar[i] - Bm[i]*Ar[i - 1])/GG->dV[i];
       dpsi[i] = (pp[i] - pm[i])/GG->dx[i];
     }
   }else if (g_dir == JDIR){  /* -- z -- */
     for (i = beg; i <= end; i++) {
       divB[i] = (Bp[i] - Bm[i])/GG->dx[i];
       dpsi[i] = (pp[i] - pm[i])/GG->dx[i];
     }
   }

  #elif GEOMETRY == POLAR

   if (g_dir == IDIR){  /* -- r -- */
     Ar  = grid[IDIR].A;
     for (i = beg; i <= end; i++) {
       divB[i] = (Bp[i]*Ar[i] - Bm[i]*Ar[i - 1])/GG->dV[i];
       dpsi[i] = (pp[i] - pm[i])/GG->dx[i];
     }
   }else if (g_dir == JDIR){  /* -- phi -- */
     r = grid[IDIR].x[g_i];
     for (i = beg; i <= end; i++) {
       divB[i] = (Bp[i] - Bm[i])/(r*GG->dx[i]);
       dpsi[i] = (pp[i] - pm[i])/(r*GG->dx[i]);
     }
   }else if (g_dir == KDIR){  /* -- z -- */
     for (i = beg; i <= end; i++) {
       divB[i] = (Bp[i] - Bm[i])/GG->dx[i];
       dpsi[i] = (pp[i] - pm[i])/GG->dx[i];
     }
   }

  #elif GEOMETRY == SPHERICAL

   if (g_dir == IDIR){  /* -- r -- */
     Ar  = grid[IDIR].A;
     for (i = beg; i <= end; i++) {
       divB[i] = (Bp[i]*Ar[i] - Bm[i]*Ar[i - 1])/GG->dV[i];
       dpsi[i] = (pp[i] - pm[i])/GG->dx[i];
     }
   }else if (g_dir == JDIR){  /* -- theta -- */
     Ath = grid[JDIR].A;
     r   = grid[IDIR].x[g_i];
     for (i = beg; i <= end; i++) {
       divB[i] = (Bp[i]*Ath[i] - Bm[i]*Ath[i - 1]) /
                 (r*GG->dV[i]);
       dpsi[i] = (pp[i] - pm[i])/(r*GG->dx[i]);

     }
   }else if (g_dir == KDIR){  /* -- phi -- */
     r = grid[IDIR].x[g_i];
     s = sin(grid[JDIR].x[g_j]);
     for (i = beg; i <= end; i++) {
       divB[i] = (Bp[i] - Bm[i])/(r*s*GG->dx[i]);
       dpsi[i] = (pp[i] - pm[i])/(r*s*GG->dx[i]);
     }
   }

  #endif

  #if BACKGROUND_FIELD == YES
   bgf = GetBackgroundField (beg - 1, end, CELL_CENTER, grid);
  #endif

/* ---------------------
     Add source terms
   -------------------- */

  for (i = beg; i <= end; i++) {

    v   = state->vh[i];
    rhs = state->rhs[i];

    EXPAND(bx = v[BX1];  ,
           by = v[BX2];  ,
           bz = v[BX3];)

    #if BACKGROUND_FIELD == YES
     EXPAND(bx += bgf[i][BX1];  ,
            by += bgf[i][BX2];  ,
            bz += bgf[i][BX3];)
    #endif

    EXPAND (rhs[MX1] -= dt*divB[i]*bx;  ,
            rhs[MX2] -= dt*divB[i]*by;  ,
            rhs[MX3] -= dt*divB[i]*bz;)

    #if HAVE_ENERGY
     rhs[ENG] -= dt*v[BXn]*dpsi[i];
    #endif
  }
}

/* ********************************************************************* */
void GLM_Init (const Data *d, const Time_Step *Dts, Grid *grid)
/*!
 * Initialize the maximum propagation speed ::glm_ch.
 *
 *
 *********************************************************************** */
{
  int i,j,k,nv;
  double dxmin, gmaxc;

  #if PHYSICS == MHD 
   if (glm_ch < 0.0){
     double        *x1, *x2, *x3;
     static double **v, **lambda, *cs2;
     Index indx;

     if (v == NULL) {
       v      = ARRAY_2D(NMAX_POINT, NVAR, double);
       lambda = ARRAY_2D(NMAX_POINT, NVAR, double);
       cs2    = ARRAY_1D(NMAX_POINT, double);
     }
     x1 = grid[IDIR].x;
     x2 = grid[JDIR].x;
     x3 = grid[KDIR].x;

     glm_ch = 0.0;
     
  /* -- X1 - g_dir -- */

     g_dir = IDIR;  SetIndexes (&indx, grid);
     KDOM_LOOP(k) JDOM_LOOP(j){
       g_j = j; g_k = k;
       IDOM_LOOP(i) for (nv = NVAR; nv--;  ) v[i][nv] = d->Vc[nv][k][j][i];
       SoundSpeed2  (v, cs2, NULL, IBEG, IEND, CELL_CENTER, grid);
       Eigenvalues  (v, cs2, lambda, IBEG, IEND);
       IDOM_LOOP(i){
         glm_ch = MAX(glm_ch, fabs(lambda[i][KFASTP]));
         glm_ch = MAX(glm_ch, fabs(lambda[i][KFASTM]));
       }
     }

  /* -- X2 - g_dir -- */

     #if DIMENSIONS >= 2
     g_dir = JDIR;  SetIndexes (&indx, grid);
     KDOM_LOOP(k) IDOM_LOOP(i){
       g_i = i; g_k = k;
       JDOM_LOOP(j) for (nv = NVAR; nv--;  ) v[j][nv] = d->Vc[nv][k][j][i];
       SoundSpeed2  (v, cs2, NULL, JBEG, JEND, CELL_CENTER, grid);
       Eigenvalues  (v, cs2, lambda, JBEG, JEND);
       JDOM_LOOP(j){
         glm_ch = MAX(glm_ch, fabs(lambda[j][KFASTP]));
         glm_ch = MAX(glm_ch, fabs(lambda[j][KFASTM]));
       }
     }
     #endif

  /* -- X3 - g_dir -- */

     #if DIMENSIONS == 3
     g_dir = KDIR;  SetIndexes (&indx, grid);
     JDOM_LOOP(j) IDOM_LOOP(i){
       g_i = i; g_j = j;
       KDOM_LOOP(k) for (nv = NVAR; nv--;  ) v[k][nv] = d->Vc[nv][k][j][i];
       SoundSpeed2  (v, cs2, NULL, KBEG, KEND, CELL_CENTER, grid);
       Eigenvalues  (v, cs2, lambda, KBEG, KEND);
       KDOM_LOOP(k){
         glm_ch = MAX(glm_ch, fabs(lambda[k][KFASTP]));
         glm_ch = MAX(glm_ch, fabs(lambda[k][KFASTM]));
       }
     }
     #endif

   }
  #elif PHYSICS == RMHD 
   if (glm_ch < 0.0) glm_ch = 1.0;
  #endif

  #ifdef PARALLEL
   MPI_Allreduce (&glm_ch, &gmaxc, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   glm_ch = gmaxc;
  #endif
}

#if COMPUTE_DIVB == YES
static double ***divB;
/* ********************************************************* */
void GLM_ComputeDivB(const State_1D *state, Grid *grid)
/*
 *
 *
 * 
 *
 *********************************************************** */
{
  int    i,j,k;
  static int old_nstep = -1;
  double   *B, *A, *dV;

  if (divB == NULL) divB = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

  if (old_nstep != g_stepNumber){
    old_nstep = g_stepNumber;
    DOM_LOOP(k,j,i) divB[k][j][i] = 0.0;
  }

  B  = state->bn;
  A  = grid[g_dir].A;
  dV = grid[g_dir].dV;
  if (g_dir == IDIR){
    k = g_k; j = g_j;
    IDOM_LOOP(i) divB[k][j][i] += (A[i]*B[i] - A[i-1]*B[i-1])/dV[i];
  }else if (g_dir == JDIR){
    k = g_k; i = g_i;
    JDOM_LOOP(j) divB[k][j][i] += (A[j]*B[j] - A[j-1]*B[j-1])/dV[j];
  }else if (g_dir == KDIR){
    j = g_j; i = g_i;
    KDOM_LOOP(k) divB[k][j][i] += (A[k]*B[k] - A[k-1]*B[k-1])/dV[k];
  } 
}

double ***GLM_GetDivB(void)
{
  return divB;
}
#endif
