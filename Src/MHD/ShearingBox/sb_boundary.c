/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Wrapper function to assign shearing-box boundary conditions.

  SB_Boundary() is a wrapper function used to assign shearing-box 
  boundary conditions on both cell-centered and staggered variables at
  the X1_BEG and X1_END boundaries.
  It is called as a regular boundary condition after periodic boundary
  conditions have already been imposed on the solution arrays.\n
  The function defines, for each type of variable (centered or staggered),
  the corresponding boundary region where shearing conditions must be
  applied using the RBox structure.
  The actual boundary condition is imposed by calling SB_SetBoundaryVar() 
  with the desired array and its box layout.

  \authors A. Mignone (mignone@ph.unito.it)\n
           G. Muscianisi (g.muscianisi@cineca.it)
  \date Jan 31, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

double sb_vy;

/* ********************************************************************* */
void SB_Boundary (const Data *d, int side, Grid *grid) 
/*! 
 * Main wrapper function used to assign shearing-box boundary conditions
 * on flow variables.
 *
 * \param d    pointer to the PLUTO Data structure
 * \param side the side of the computational domain (X1_BEG or X1_END) 
 * \param grid pointer to an array of Grid structures
 *
 * \return This function has no return value.
 *
 * \todo Check if sb_vy needs to be global.
 *********************************************************************** */
{
  int    i, j, k, nv;
  double t, Lx;
  RBox   box;

  Lx    = g_domEnd[IDIR] - g_domBeg[IDIR];
  sb_vy = fabs(2.0*SB_A*Lx);

{
t = g_time;
#if TIME_STEPPING == RK2
 if (g_intStage == 2) t += g_dt;
#elif TIME_STEPPING == RK3
 if (g_intStage == 2) t += 0.5*g_dt;
 if (g_intStage == 3) t += g_dt;
#endif
#ifdef CTU
 if (g_intStage == 2) t += 0.5*g_dt;
#endif
}

/* -------------------------------------------------
                  X1 Beg Boundary
   ------------------------------------------------- */

  if (side == X1_BEG){

  /* ---- loop on cell-centered variables ---- */

    box.ib = 0; box.ie = IBEG-1;
    box.jb = 0; box.je = NX2_TOT-1;
    box.kb = 0; box.ke = NX3_TOT-1;
    for (nv = 0; nv < NVAR; nv++){

      #ifdef STAGGERED_MHD
       D_EXPAND(if (nv == BX) continue;  ,
                if (nv == BY) continue;  ,
                if (nv == BZ) continue;)
      #endif

      SB_SetBoundaryVar(d->Vc[nv], &box, side, t, grid);
      if (nv == VY) {
        X1_BEG_LOOP(k,j,i) d->Vc[nv][k][j][i] += sb_vy;
      }
    }  /* -- end loop on cell-centered variables -- */

    #ifdef STAGGERED_MHD
      box.ib =  0; box.ie = IBEG-1;
      box.jb = -1; box.je = NX2_TOT-1;
      box.kb =  0; box.ke = NX3_TOT-1;
      SB_SetBoundaryVar(d->Vs[BX2s], &box, side, t, grid);
      #if DIMENSIONS == 3
       box.ib =  0; box.ie = IBEG-1;
       box.jb =  0; box.je = NX2_TOT-1;
       box.kb = -1; box.ke = NX3_TOT-1;
       SB_SetBoundaryVar(d->Vs[BX3s], &box, side, t, grid);
      #endif
    #endif /* STAGGERED_MHD */
  } /* -- END side X1_BEG -- */

/* -------------------------------------------------
                  X1 End Boundary
   ------------------------------------------------- */

  if (side == X1_END){

  /* ---- loop on cell-centered variables ---- */

    box.ib = IEND+1; box.ie = NX1_TOT-1;
    box.jb =      0; box.je = NX2_TOT-1;
    box.kb =      0; box.ke = NX3_TOT-1;
    for (nv = 0; nv < NVAR; nv++){

      #ifdef STAGGERED_MHD
       D_EXPAND(if (nv == BX) continue;  ,
                if (nv == BY) continue;  ,
                if (nv == BZ) continue;)
      #endif

      SB_SetBoundaryVar(d->Vc[nv], &box, side, t, grid);
      if (nv == VY){
        X1_END_LOOP(k,j,i) d->Vc[nv][k][j][i] -= sb_vy;
      }
    }  /* -- end loop on cell-centered variables -- */

    #ifdef STAGGERED_MHD
     box.ib = IEND+1; box.ie = NX1_TOT-1;
     box.jb =     -1; box.je = NX2_TOT-1;
     box.kb =      0; box.ke = NX3_TOT-1;
     SB_SetBoundaryVar(d->Vs[BX2s], &box, side, t, grid);
     #if DIMENSIONS == 3
      box.ib = IEND+1; box.ie = NX1_TOT-1;
      box.jb =      0; box.je = NX2_TOT-1;
      box.kb =     -1; box.ke = NX3_TOT-1;
      SB_SetBoundaryVar(d->Vs[BX3s], &box, side, t, grid);
     #endif /* DIMENSIONS == 3 */
    #endif /* STAGGERED_MHD */
  }
}
