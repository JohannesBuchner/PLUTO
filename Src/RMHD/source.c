#include "pluto.h"

#if DIVB_CONTROL == EIGHT_WAVES
/* *************************************************************************  */
void POWELL_DIVB_SOURCE(const State_1D *state, int is, int ie, Grid *grid)
/*
 *
 *
 ************************************************************************** */
{
  int    i;
  double divB;
  double g_2, vB, *vc;
  Grid   *GG;

  GG = grid + g_dir;
  
/* ------------------------------------------------------------
              POWELL's divB monopole source term
   ------------------------------------------------------------ */

  for (i = is; i <= ie; i++) {

    vc = state->vh[i];

/*     Make (divB) contribution from the normal direction (1)     */

    divB = (state->vR[i][BXn]     + state->vL[i][BXn])*GG->A[i] -
           (state->vR[i - 1][BXn] + state->vL[i - 1][BXn])*GG->A[i - 1];
/*
    divB = state->bn[i]*GG->A[i] - state->bn[i - 1]*GG->A[i - 1];
*/
    divB /= 2.0*GG->dV[i];

    vB  = EXPAND(vc[VX1]*vc[BX1], + vc[VX2]*vc[BX2], + vc[VX3]*vc[BX3]);
    g_2 = EXPAND(vc[VX1]*vc[VX1], + vc[VX2]*vc[VX2], + vc[VX3]*vc[VX3]);
    g_2 = 1.0 - g_2;

    EXPAND (state->src[i][MX1] = -divB*(vc[BX1]*g_2 + vB*vc[VX1]);  ,
            state->src[i][MX2] = -divB*(vc[BX2]*g_2 + vB*vc[VX2]);  ,
            state->src[i][MX3] = -divB*(vc[BX3]*g_2 + vB*vc[VX3]);)

    state->src[i][ENG] = -divB*vB;

    EXPAND (state->src[i][BX1] = -divB*vc[VX1];  ,
            state->src[i][BX2] = -divB*vc[VX2];  ,
            state->src[i][BX3] = -divB*vc[VX3];)

  }
}

/* *********************************************************************  */
void HLL_DIVB_SOURCE (const State_1D *state, double **Uhll, int beg, 
                      int end, Grid *grid)
/* 
 *
 * PURPOSE
 *
 *   Include div.B source term to momentum, induction
 *   and energy equation. Used in conjunction with 
 *   an HLL-type Riemann solver. 
 * 
 * LAST_MODIFIED
 *
 *   Setp 18th 2008, by Andrea Mignone  (mignone@to.astro.it)
 *
 *********************************************************************** */
{
  int i, nv;
  double vc[NVAR], *A, *src, *vm;
  double r, s, vB;
  static double *divB, *vp;
  Grid *GG;

  if (divB == NULL){
    divB = ARRAY_1D(NMAX_POINT, double);
    vp   = ARRAY_1D(NMAX_POINT, double);
  }

  GG = grid + g_dir;
  vm = vp - 1;
  A  = grid[g_dir].A;

/* --------------------------------------------
    Compute normal component of the field 
   -------------------------------------------- */

  for (i = beg - 1; i <= end; i++) {
    vp[i] = state->flux[i][RHO] < 0.0 ? state->vR[i][BXn]: state->vL[i][BXn];
  }

/* --------------------------------------------
    Compute div.B contribution from the normal 
    direction (1) in different geometries 
   -------------------------------------------- */

  
  #if GEOMETRY == CARTESIAN

   for (i = beg; i <= end; i++) {
     divB[i] = (vp[i] - vm[i])/GG->dx[i];
   }

  #elif GEOMETRY == CYLINDRICAL

   if (g_dir == IDIR){   /* -- r -- */
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i]*A[i] - vm[i]*A[i - 1])/GG->dV[i];
     }
   }else if (g_dir == JDIR){  /* -- z -- */
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i] - vm[i])/GG->dx[i];
     }
   }

  #elif GEOMETRY == POLAR

   if (g_dir == IDIR){  /* -- r -- */
     Ar  = grid[IDIR].A;
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i]*A[i] - vm[i]*A[i - 1])/GG->dV[i];
     }
   }else if (g_dir == JDIR){  /* -- phi -- */
     r = grid[IDIR].x[*g_i];
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i] - vm[i])/(r*GG->dx[i]);
     }
   }else if (g_dir == KDIR){  /* -- z -- */
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i] - vm[i])/GG->dx[i];
     }
   }

  #elif GEOMETRY == SPHERICAL

   if (g_dir == IDIR){  /* -- r -- */
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i]*A[i] - vm[i]*A[i - 1])/GG->dV[i];
     }
   }else if (g_dir == JDIR){  /* -- theta -- */
     r   = grid[IDIR].x[*g_i];
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i]*A[i] - vm[i]*A[i - 1])/(r*GG->dV[i]);
     }
   }else if (g_dir == KDIR){  /* -- phi -- */
     r = grid[IDIR].x[*g_i];
     s = sin(grid[JDIR].x[*g_j]);
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i] - vm[i])/(r*s*GG->dx[i]);
     }
   }

  #endif

  /* -----------------------------------------
          compute total source terms
     ----------------------------------------- */
     
  for (i = beg; i <= end; i++) {

    src = state->src[i];
    for (nv = NFLX; nv--;  ) vc[nv] = state->vh[i][nv];

    vB = EXPAND(vc[VX1]*vc[BX1], + vc[VX2]*vc[BX2], + vc[VX3]*vc[BX3]);
    src[RHO] = 0.0;
    EXPAND(src[MX1] = 0.0;  ,
           src[MX2] = 0.0;  ,
           src[MX3] = 0.0;)

    EXPAND(src[BX1] = -vc[VX1]*divB[i];  ,
           src[BX2] = -vc[VX2]*divB[i];  ,
           src[BX3] = -vc[VX3]*divB[i];)

    #if EOS != ISOTHERMAL
     src[ENG] = 0.0;
    #endif

  }  
}
#endif
