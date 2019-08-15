/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute source terms for Powell formulation.

  This file contains two implementations of Powell's source term 
  in the 8-wave formulation for MHD:
  \f[
    Q\nabla\cdot\vec{B} \approx Q_i\, 
    \frac{A_{i\HALF}B_{i+\HALF} - A_{i-\HALF}B_{i-\HALF}}{\Delta{\cal V}_i}
  \f]
  where \c B is the magnetic field in the normal direction (::g_dir),
  \c Q is a fluid quantity, \c A is the area and the denominator is
  the cell volume.
  The first function, Roe_DivBSource() is called by Roe_Solver() and 
  TVDLF_Solver() and computes the normal component using arithmetic average 
  of the left and right states: \f$ B_{i+\HALF} = (B^+_i + B^-_{i+1})/2 \f$.
  The second implementation contained in HLL_DivBSource() computes the term 
  using upwinding:
  \f[
     B_{i+\HALF} = \left\{\begin{array}{ll}
      B^+_i      & \qquad{\rm if}\quad F_{\rho,i+\HALF} > 0\\ \noalign{\medskip}
      B^-_{i+1}  & \qquad{\rm otherwise}
      \end{array}\right.
  \f]      
  
  \b Reference:
    - "A positive conservative method for MHD based on HLL and Roe methods"
       P. Janhunen, JCP (2000), 160, 649.
       
  \authors A. Mignone (mignone@ph.unito.it)
  \date    March 3, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if DIVB_CONTROL == EIGHT_WAVES
/* ********************************************************************* */
void Roe_DivBSource (const Sweep *sweep, int beg, int end, Grid *grid)
/*!
 *  Include Powell div.B source term to momentum, induction
 *  and energy equation for Roe and TVDLF solvers.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;

  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);

  double btx, bty, btz, bx, by, bz, vx, vy, vz;
  double r, s;
  double *vc;
  double **bgf, *src;
  double *dx   = grid->dx[g_dir];
  double ***A  = grid->A[g_dir];
  double ***dV = grid->dV;
  static double *divB, *Bn;

/* --------------------------------------------
   0. Allocate memory
   -------------------------------------------- */

  if (divB == NULL){
    divB = ARRAY_1D(NMAX_POINT, double);
    Bn   = ARRAY_1D(NMAX_POINT, double);
  }

/* --------------------------------------------
   1. Compute magnetic field normal component 
      interface value by arithmetic averaging
   -------------------------------------------- */

  for (i = beg - 1; i <= end; i++) {
    Bn[i] = 0.5*(stateL->v[i][BXn] + stateR->v[i][BXn]);
  }
  
/* --------------------------------------------
   2. Compute div.B contribution from the normal 
      direction in different geometries 
   -------------------------------------------- */
  
#if GEOMETRY == CARTESIAN

  for (i = beg; i <= end; i++) divB[i] = (Bn[i] - Bn[i-1])/dx[i];

#else

  if (g_dir == IDIR){
    j = g_j;
    k = g_k;
    for (i = beg; i <= end; i++) {
      divB[i] = (A[k][j][i]*Bn[i] - A[k][j][i-1]*Bn[i-1])/dV[k][j][i];
    }
  }else if (g_dir == JDIR){
    i = g_i;
    k = g_k;
    for (j = beg; j <= end; j++) {
      divB[j] = (A[k][j][i]*Bn[j] - A[k][j-1][i]*Bn[j-1])/dV[k][j][i];
    }
  }else if (g_dir == KDIR){
    i = g_i;
    j = g_j;
    for (k = beg; k <= end; k++) {
      divB[k] = (A[k][j][i]*Bn[k] - A[k-1][j][i]*Bn[k-1])/dV[k][j][i];
    }
  }

#endif

#if BACKGROUND_FIELD == YES
  GetBackgroundField (stateC, beg - 1, end, CELL_CENTER, grid);
#endif

/* -------------------------------------------
   3. Compute Powell's source term
   ------------------------------------------- */

  for (i = beg; i <= end; i++) {
  
    vc  = stateC->v[i];
    src = sweep->src[i];
    bgf = stateC->Bbck;

    EXPAND (vx = vc[VX1];  ,
            vy = vc[VX2];  ,
            vz = vc[VX3];)

    EXPAND (bx = btx = vc[BX1];  ,
            by = bty = vc[BX2];  ,
            bz = btz = vc[BX3];)

    #if BACKGROUND_FIELD == YES
    btx += bgf[i][BX1];
    bty += bgf[i][BX2];
    btz += bgf[i][BX3];
    #endif

    src[RHO] = 0.0;
    EXPAND (src[MX1] = -divB[i]*btx;  ,
            src[MX2] = -divB[i]*bty;  ,
            src[MX3] = -divB[i]*btz;)

    #if HAVE_ENERGY
    src[ENG] = -divB[i]*(EXPAND(vx*bx, +vy*by, +vz*bz));
    #endif
    EXPAND (src[BX1] = -divB[i]*vx;   ,
            src[BX2] = -divB[i]*vy;   ,
            src[BX3] = -divB[i]*vz;)
  }
}

/* *********************************************************************  */
void HLL_DivBSource (const Sweep *sweep, double **Uhll, 
                     int beg, int end, Grid *grid)
/*! 
 *  Include div.B source term to momentum, induction
 *  and energy equation. Used in conjunction with 
 *  an HLL-type Riemann solver. 
 *********************************************************************** */
{
  int i, j, k, nv;

  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);

  double *vc, *src, vB;
  double *dx   = grid->dx[g_dir];
  double ***A  = grid->A[g_dir];
  double ***dV = grid->dV;
  static double *divB, *Bn;

/* --------------------------------------------
   0. Allocate memory
   -------------------------------------------- */

  if (divB == NULL){
    divB = ARRAY_1D(NMAX_POINT, double);
    Bn   = ARRAY_1D(NMAX_POINT, double);
  }

/* --------------------------------------------
   1. Compute normal component of the field
      using upwinding
   -------------------------------------------- */

  for (i = beg - 1; i <= end; i++) {
    Bn[i] = sweep->flux[i][RHO] < 0.0 ? stateR->v[i][BXn]: stateL->v[i][BXn];
  }

/* --------------------------------------------
   2. Compute div.B contribution from the normal 
      direction in different geometries 
   -------------------------------------------- */

#if GEOMETRY == CARTESIAN

  for (i = beg; i <= end; i++) divB[i] = (Bn[i] - Bn[i-1])/dx[i];

#else

  if (g_dir == IDIR){
    j = g_j;
    k = g_k;
    for (i = beg; i <= end; i++) {
      divB[i] = (A[k][j][i]*Bn[i] - A[k][j][i-1]*Bn[i-1])/dV[k][j][i];
    }
  }else if (g_dir == JDIR){
    i = g_i;
    k = g_k;
    for (j = beg; j <= end; j++) {
      divB[j] = (A[k][j][i]*Bn[j] - A[k][j-1][i]*Bn[j-1])/dV[k][j][i];
    }
  }else if (g_dir == KDIR){
    i = g_i;
    j = g_j;
    for (k = beg; k <= end; k++) {
      divB[k] = (A[k][j][i]*Bn[k] - A[k-1][j][i]*Bn[k-1])/dV[k][j][i];
    }
  }
#endif

/* -----------------------------------------
   3. Compute divB source terms
   ----------------------------------------- */
     
  for (i = beg; i <= end; i++) {

    vc  = stateC->v[i];
    src = sweep->src[i];

    vB = EXPAND(vc[VX1]*vc[BX1], + vc[VX2]*vc[BX2], + vc[VX3]*vc[BX3]);
    src[RHO] = 0.0;
    EXPAND(src[MX1] = -vc[BX1]*divB[i];  ,
           src[MX2] = -vc[BX2]*divB[i];  ,
           src[MX3] = -vc[BX3]*divB[i];)

    EXPAND(src[BX1] = -vc[VX1]*divB[i];  ,
           src[BX2] = -vc[VX2]*divB[i];  ,
   	       src[BX3] = -vc[VX3]*divB[i];)

    #if HAVE_ENERGY 
    src[ENG] = -vB*divB[i];
    #endif
  }  
}
#endif  /* DIVB_CONTROL == EIGHT_WAVES */
