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
  \date    June 8, 2007
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if DIVB_CONTROL == EIGHT_WAVES
/* ********************************************************************* */
void Roe_DivBSource (const State_1D *state, int is, int ie, Grid *grid)
/*!
 *  Include Powell div.B source term to momentum, induction
 *  and energy equation for Roe and TVDLF solvers.
 *
 *********************************************************************** */
{
  int    i;
  double btx, bty, btz, bx, by, bz, vx, vy, vz;
  double r, s;
  double *Ar, *Ath;
  double *vm, **bgf;
  double *src, *v;
  static double *divB, *vp;
  Grid   *GG;

  if (divB == NULL){
    divB = ARRAY_1D(NMAX_POINT, double);
    vp   = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------- ---------------------
    compute magnetic field normal component 
    interface value by arithmetic averaging
   -------------------------------------------- */

  for (i = is - 1; i <= ie; i++) {
    vp[i] = 0.5*(state->vL[i][BXn] + state->vR[i][BXn]);
  }
  vm  = vp - 1;
  GG  = grid + g_dir;
  
/* --------------------------------------------
    Compute div.B contribution from the normal 
    direction (1) in different geometries 
   -------------------------------------------- */
  
  #if GEOMETRY == CARTESIAN

   for (i = is; i <= ie; i++) {
     divB[i] = (vp[i] - vm[i])/GG->dx[i];
   }

  #elif GEOMETRY == CYLINDRICAL

   if (g_dir == IDIR){   /* -- r -- */
     Ar  = grid[IDIR].A;
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i]*Ar[i] - vm[i]*Ar[i - 1])/GG->dV[i];
     }
   }else if (g_dir == JDIR){  /* -- z -- */
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i] - vm[i])/GG->dx[i];
     }
   }

  #elif GEOMETRY == POLAR

   if (g_dir == IDIR){  /* -- r -- */
     Ar  = grid[IDIR].A;
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i]*Ar[i] - vm[i]*Ar[i - 1])/GG->dV[i];
     }
   }else if (g_dir == JDIR){  /* -- phi -- */
     r = grid[IDIR].x[g_i];
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i] - vm[i])/(r*GG->dx[i]);
     }
   }else if (g_dir == KDIR){  /* -- z -- */
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i] - vm[i])/GG->dx[i];
     }
   }

  #elif GEOMETRY == SPHERICAL

   if (g_dir == IDIR){  /* -- r -- */
     Ar  = grid[IDIR].A;
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i]*Ar[i] - vm[i]*Ar[i - 1])/GG->dV[i];
     }
   }else if (g_dir == JDIR){  /* -- theta -- */
     Ath = grid[JDIR].A;
     r   = grid[IDIR].x[g_i];
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i]*Ath[i] - vm[i]*Ath[i - 1]) /
                 (r*GG->dV[i]);
     }
   }else if (g_dir == KDIR){  /* -- phi -- */
     r = grid[IDIR].x[g_i];
     s = sin(grid[JDIR].x[g_j]);
     for (i = is; i <= ie; i++) {
       divB[i] = (vp[i] - vm[i])/(r*s*GG->dx[i]);
     }
   }

  #endif

  #if BACKGROUND_FIELD == YES
   bgf = GetBackgroundField (is - 1, ie, CELL_CENTER, grid);
  #endif

/* -------------------------------------------
          Compute Powell's source term
   ------------------------------------------- */

  for (i = is; i <= ie; i++) {

    v   = state->vh[i];
    src = state->src[i];

    EXPAND (vx = v[VX1];  ,
            vy = v[VX2];  ,
            vz = v[VX3];)

    EXPAND (bx = btx = v[BX1];  ,
            by = bty = v[BX2];  ,
            bz = btz = v[BX3];)

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
void HLL_DivBSource (const State_1D *state, double **Uhll, 
                     int beg, int end, Grid *grid)
/*! 
 *  Include div.B source term to momentum, induction
 *  and energy equation. Used in conjunction with 
 *  an HLL-type Riemann solver. 
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

/* -------------------------------------------------------
    Compute normal component of the field using upwinding
   ------------------------------------------------------- */

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

   if (g_dir == IDIR){        /* -- r -- */
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i]*A[i] - vm[i]*A[i - 1])/GG->dV[i];
     }
   }else if (g_dir == JDIR){  /* -- phi -- */
     r = grid[IDIR].x[g_i];
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
     r   = grid[IDIR].x[g_i];
     for (i = beg; i <= end; i++) {
       divB[i] = (vp[i]*A[i] - vm[i]*A[i - 1])/(r*GG->dV[i]);
     }
   }else if (g_dir == KDIR){  /* -- phi -- */
     r = grid[IDIR].x[g_i];
     s = sin(grid[JDIR].x[g_j]);
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
#endif
