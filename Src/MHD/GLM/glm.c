/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief GLM module implementation.

  Contains functions for the GLM module.

  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date    Dec 04, 2017

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
void GLM_Solve (const Sweep *sweep, int beg, int end, Grid *grid)
/*!
 * Solve the 2x2 linear hyperbolic GLM-MHD system given by the divergence 
 * cleaning approach.
 * Modify inteface states (Bx and psi components) for input to full
 * Riemann problem.
 * We use Eq. (42) of Dedner et al (2002)
 *
 * \param [in,out] sweep pointer to a Sweep structure
 * \param [in]     beg    starting index of computation
 * \param [in]     end    final index of computation
 * \param [in]     grid   pointer to Grid structure
 *
 * The purpose of this function is two-fold:
 *
 * 1. assign a unique value to the normal component of magnetic field 
 *    and to te scalar psi \e before the actual Riemann solver is called.
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
  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);
  double Bm, psim;
  double dB, dpsi;
  double **v, *vL, *vR;
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
    v = sweep->vn;
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
    vL = stateL->v[i];
    vR = stateR->v[i];


    #ifdef CTU
    if (nflag){
      vL[BXn] = Bxp[i];
      vR[BXn] = Bxm[i+1];
    }
    #endif

    dB   = vR[BXn]     - vL[BXn];
    dpsi = vR[PSI_GLM] - vL[PSI_GLM];

    Bm   = 0.5*(vL[BXn]     + vR[BXn])     - 0.5*dpsi/glm_ch;
    psim = 0.5*(vL[PSI_GLM] + vR[PSI_GLM]) - 0.5*glm_ch*dB;

    sweep->bn[i] = vL[BXn] = vR[BXn] = Bm;
    vL[PSI_GLM] = vR[PSI_GLM] = psim;

    #if PHYSICS == RMHD && RESISTIVITY != NO
    {
      double dE, dphi, Em, phim;

      dE   = vR[EXn]     - vL[EXn];
      dpsi = vR[PHI_GLM] - vL[PHI_GLM];

      Em   = 0.5*(vL[EXn]     + vR[EXn])     - 0.5*dphi/glm_ch;
      phim = 0.5*(vL[PHI_GLM] + vR[PHI_GLM]) - 0.5*glm_ch*dE;

      vL[EXn] = vR[EXn] = Em;
      vL[PHI_GLM] = vR[PHI_GLM] = phim;
    }
    #endif
  }

  PrimToCons(stateL->v, stateL->u, beg, end);
  PrimToCons(stateR->v, stateR->u, beg, end);

  #if GLM_COMPUTE_DIVB == YES
  if (g_intStage == 1) GLM_ComputeDivB(sweep, grid);  
  #endif
  #if GLM_COMPUTE_DIVE == YES
  if (g_intStage == 1) GLM_ComputeDivE(sweep, grid);  
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
  dx   = grid->dx[IDIR][IBEG];
  dtdx = dt/dx;
#endif

  cp   = sqrt(dx*glm_ch/GLM_ALPHA);
  scrh = dtdx*glm_ch*GLM_ALPHA; 

  scrh = exp(-scrh);
#ifdef CHOMBO
  DOM_LOOP(k,j,i) Q[k][j][i][PSI_GLM] *= scrh;
#else
  DOM_LOOP(k,j,i) {
    #if PHYSICS == MHD || (PHYSICS == RMHD && RESISTIVITY == NO)
    Q[PSI_GLM][k][j][i] *= scrh;
    #elif (PHYSICS == RMHD) && (RESISTIVITY != NO)
    scrh = exp(-dt);
    Q[PSI_GLM][k][j][i] *= scrh;
    Q[PHI_GLM][k][j][i] *= scrh;
    #endif
  }  
#endif
  return;
}

/* ********************************************************************* */
void GLM_ExtendedSource (const Sweep *sweep, double dt,
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
  int    i,j,k;
  const State *stateC = &(sweep->stateC);
  double bx, by, bz, ch2;
  double *rhs, *v, *Bm, *pm;
  double ***A, *dx, *r, *th;
  static double *divB, *dpsi, *Bp, *pp;
/* 
  #if PHYSICS == RMHD
   print ("! EGLM not working for the RMHD module\n");
   QUIT_PLUTO(1);
  #endif
*/

/* ----------------------------------------------------
   0. Allocate memory and create pointer shortcuts
   ---------------------------------------------------- */

  if (divB == NULL){
    divB = ARRAY_1D(NMAX_POINT, double);
    dpsi = ARRAY_1D(NMAX_POINT, double);
    Bp   = ARRAY_1D(NMAX_POINT, double);
    pp   = ARRAY_1D(NMAX_POINT, double);
  }
  r  = grid->x[IDIR];
  th = grid->x[JDIR];
  dx = grid->dx[g_dir];
  A  = grid->A[g_dir];

/* ---------------------------------------------
   1. Compute magnetic field normal component 
      interface value by arithmetic averaging
   -------------------------------------------- */

  ch2 = glm_ch*glm_ch;
  for (i = beg - 1; i <= end; i++) {
    Bp[i] = sweep->flux[i][PSI_GLM]/ch2;
    pp[i] = sweep->flux[i][BXn];
  }
  Bm  = Bp - 1;
  pm  = pp - 1;
  
/* -----------------------------------------------
   2. Compute div.B contribution from the normal 
      direction (1) in different geometries 
   ----------------------------------------------- */
  
    
#if GEOMETRY == CARTESIAN
  for (i = beg; i <= end; i++) {
    divB[i] = (Bp[i] - Bm[i])/dx[i];
    dpsi[i] = (pp[i] - pm[i])/dx[i];
  }
#else
  if (g_dir == IDIR){
    j = g_j;
    k = g_k;
    for (i = beg; i <= end; i++) {
      divB[i] = (A[k][j][i]*Bp[i] - A[k][j][i-1]*Bm[i])/grid->dV[k][j][i];
      dpsi[i] = (pp[i] - pm[i])/dx[i];
    }
  }else if (g_dir == JDIR){
    i = g_i;
    k = g_k;
    for (j = beg; j <= end; j++) {
      divB[j] = (A[k][j][i]*Bp[j] - A[k][j-1][i]*Bm[j])/grid->dV[k][j][i];
      dpsi[j] = (pp[j] - pm[j])/dx[j];
      #if GEOMETRY == POLAR || GEOMETRY == SPHERICAL
      dpsi[j] /= r[i];
      #endif
    }
  }else if (g_dir == KDIR){
    i = g_i;
    j = g_j;
    for (k = beg; k <= end; k++) {
      divB[k] = (A[k][j][i]*Bp[k] - A[k-1][j][i]*Bm[k])/grid->dV[k][j][i];
      dpsi[k] = (pp[k] - pm[k])/dx[k];
      #if GEOMETRY == SPHERICAL
      dpsi[k] /= r[i]*sin(th[j]);
      #endif
    }
  }
#endif

/* ---------------------------------------------------
   3. Add source terms
   --------------------------------------------------- */

#if BACKGROUND_FIELD == YES
  GetBackgroundField (stateC, beg - 1, end, CELL_CENTER, grid);
#endif

  for (i = beg; i <= end; i++) {
    v   = stateC->v[i];
    rhs = sweep->rhs[i];

    EXPAND(bx = v[BX1];  ,
           by = v[BX2];  ,
           bz = v[BX3];)

    #if BACKGROUND_FIELD == YES
    EXPAND(bx += stateC->Bbck[i][IDIR];  ,
           by += stateC->Bbck[i][JDIR];  ,
           bz += stateC->Bbck[i][KDIR];)
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
void GLM_Init (const Data *d, const timeStep *Dts, Grid *grid)
/*!
 * Initialize the maximum propagation speed ::glm_ch at the beginning
 * of integration cycle. 
 *
 *********************************************************************** */
{
  int    i,j,k,nv;
  double dxmin, gmaxc;

#if PHYSICS == MHD 
  if (glm_ch < 0.0){
    int nbeg, nend, *in;
    double *x1, *x2, *x3;
    RBox   sweepBox;
    static double **v, **lambda, *cs2;
    State state;

    if (v == NULL) {
      v      = ARRAY_2D(NMAX_POINT, NVAR, double);
      lambda = ARRAY_2D(NMAX_POINT, NVAR, double);
      cs2    = ARRAY_1D(NMAX_POINT, double);
    }
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

    state.v  = v;
    state.a2 = cs2;
     
    for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){

      RBoxDefine(IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &sweepBox);
      RBoxSetDirections (&sweepBox, g_dir);
      SetVectorIndices (g_dir);

      nbeg = *sweepBox.nbeg;
      nend = *sweepBox.nend;

      BOX_TRANSVERSE_LOOP(&sweepBox, k,j,i){
        in  = sweepBox.n;
        g_i = i; g_j = j; g_k = k;
        for (*in = nbeg; *in <= nend; (*in)++){
          NVAR_LOOP(nv) v[*in][nv] = d->Vc[nv][k][j][i];
        }
        SoundSpeed2  (&state, nbeg, nend, CELL_CENTER, grid);
        Eigenvalues  (v, cs2, lambda, nbeg, nend);
        for (*in = nbeg; *in <= nend; (*in)++){
          glm_ch = MAX(glm_ch, fabs(lambda[*in][KFASTP]));
          glm_ch = MAX(glm_ch, fabs(lambda[*in][KFASTM]));
        }    
      }
    }
  }
#elif PHYSICS == RMHD 
  if (glm_ch < 0.0) glm_ch = 1.0;
#endif

#ifdef PARALLEL
  MPI_Allreduce (&glm_ch, &gmaxc, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  glm_ch = gmaxc;
#endif
}

#if GLM_COMPUTE_DIVB == YES
static double ***divB;
/* ********************************************************* */
void GLM_ComputeDivB(const Sweep *sweep, Grid *grid)
/*
 *
 *
 * 
 *
 *********************************************************** */
{
  int    i,j,k;
  static int old_nstep = -1;
  double *B    = sweep->bn;
  double ***A  = grid->A[g_dir];
  double ***dV = grid->dV;

  if (divB == NULL) divB = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

  if (old_nstep != g_stepNumber){
    old_nstep = g_stepNumber;
    DOM_LOOP(k,j,i) divB[k][j][i] = 0.0;
  }

  if (g_dir == IDIR){
    k = g_k; j = g_j;
    IDOM_LOOP(i) {
      divB[k][j][i] += (A[k][j][i]*B[i] - A[k][j][i-1]*B[i-1])/dV[k][j][i];
    }
  }else if (g_dir == JDIR){
    k = g_k; i = g_i;
    JDOM_LOOP(j) {
      divB[k][j][i] += (A[k][j][i]*B[j] - A[k][j-1][i]*B[j-1])/dV[k][j][i];
    }
  }else if (g_dir == KDIR){
    j = g_j; i = g_i;
    KDOM_LOOP(k) {
      divB[k][j][i] += (A[k][j][i]*B[k] - A[k-1][j][i]*B[k-1])/dV[k][j][i];
    }
  } 
}

double ***GLM_GetDivB(void)
{
  return divB;
}
#endif

#if (GLM_COMPUTE_DIVE == YES) && (PHYSICS == RMHD && RESISTIVITY != NO)
static double ***divE;
/* ********************************************************* */
void GLM_ComputeDivE(const Sweep *sweep, Grid *grid)
/*!
 *  Compute the divergence of E using Godunov fluxes
 *  previously obtained at cell interfaces.
 *  This function may be used in Resistive RMHD.
 *
 *********************************************************** */
{
  int    i,j,k;
  static int old_nstep = -1;
  double **flux = sweep->flux;
  double ***A   = grid->A[g_dir];
  double ***dV  = grid->dV;

  if (divE == NULL) divE = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

  if (old_nstep != g_stepNumber){
    old_nstep = g_stepNumber;
    DOM_LOOP(k,j,i) divE[k][j][i] = 0.0;
  }

  if (g_dir == IDIR){
    k = g_k; j = g_j;
    IDOM_LOOP(i) {
      divE[k][j][i] += (  A[k][j][i]*flux[i][PHI_GLM]
                        - A[k][j][i-1]*flux[i-1][PHI_GLM])/dV[k][j][i];
    }
  }else if (g_dir == JDIR){
    k = g_k; i = g_i;
    JDOM_LOOP(j) {
      divE[k][j][i] += (  A[k][j][i]*flux[j][PHI_GLM]
                        - A[k][j-1][i]*flux[j-1][PHI_GLM])/dV[k][j][i];
    }
  }else if (g_dir == KDIR){
    j = g_j; i = g_i;
    KDOM_LOOP(k) {
      divE[k][j][i] += (  A[k][j][i]*flux[k][PHI_GLM]
                        - A[k-1][j][i]*flux[k-1][PHI_GLM])/dV[k][j][i];
    }
  }
}

double ***GLM_GetDivE(void)
{
  return divE;
}
#endif
