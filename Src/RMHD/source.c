#include "pluto.h"

#if DIVB_CONTROL == EIGHT_WAVES
/* *************************************************************************  */
void POWELL_DIVB_SOURCE(const Sweep *sweep, int beg, int end, Grid *grid)
/*
 *
 *
 ************************************************************************** */
{
  int    i,j,k;

  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);

  double g_2, vB, *vc;
  double *dx   = grid->dx[g_dir];
  double ***A  = grid->A[g_dir];
  double ***dV = grid->dV;
  static double *divB, *Bn;
  
  if (divB == NULL) {
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

/* ---------------------------------------------------------
   2. Compute normal component of div(B)
   --------------------------------------------------------- */

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

/* ------------------------------------------------------------
              POWELL's divB monopole source term
   ------------------------------------------------------------ */

  for (i = beg; i <= end; i++) {

    vc = stateC->v[i];

    vB  = EXPAND(vc[VX1]*vc[BX1], + vc[VX2]*vc[BX2], + vc[VX3]*vc[BX3]);
    g_2 = EXPAND(vc[VX1]*vc[VX1], + vc[VX2]*vc[VX2], + vc[VX3]*vc[VX3]);
    g_2 = 1.0 - g_2;

    EXPAND (sweep->src[i][MX1] = -divB[i]*(vc[BX1]*g_2 + vB*vc[VX1]);  ,
            sweep->src[i][MX2] = -divB[i]*(vc[BX2]*g_2 + vB*vc[VX2]);  ,
            sweep->src[i][MX3] = -divB[i]*(vc[BX3]*g_2 + vB*vc[VX3]);)

    sweep->src[i][ENG] = -divB[i]*vB;

    EXPAND (sweep->src[i][BX1] = -divB[i]*vc[VX1];  ,
            sweep->src[i][BX2] = -divB[i]*vc[VX2];  ,
            sweep->src[i][BX3] = -divB[i]*vc[VX3];)

  }
}

/* *********************************************************************  */
void HLL_DIVB_SOURCE (const Sweep *sweep, double **Uhll, int beg, 
                      int end, Grid *grid)
/*!
 * Include div.B source term to momentum, induction
 * and energy equation. Used in conjunction with 
 * an HLL-type Riemann solver. 
 * 
 * LAST_MODIFIED
 *
 *   Setp 18th 2008, by Andrea Mignone  (mignone@to.astro.it)
 *
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
    EXPAND(src[MX1] = 0.0;  ,
           src[MX2] = 0.0;  ,
           src[MX3] = 0.0;)

    EXPAND(src[BX1] = -vc[VX1]*divB[i];  ,
           src[BX2] = -vc[VX2]*divB[i];  ,
   	       src[BX3] = -vc[VX3]*divB[i];)

    #if HAVE_ENERGY 
    src[ENG] = 0.0;
    #endif
  }  
}
#endif  /* DIVB_CONTROL == EIGHT_WAVES */
