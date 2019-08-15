/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Perform index permutation and set domain integration indexes.

  The function SetIndexes() performs two different tasks, depending on
  the value of ::g_dir and the time-stepping algorithm:
  - perform a cyclic permutation of the array indexes corresponding to
    vector components (velocity, momentum and magnetic field);
  - set the starting and final grid indexes (IBEG, IEND,...) of the 
    computational domain before commencing integration.

  The indexes coincide with the usual ones (IBEG, IEND,...) most of the time,
  but can be expanded one further zone either in the transverse or normal
  directions or both. This depends on the integration algorithm.

  With cell-centered fields, the following table holds:
  \verbatim
                        RK        CTU
                   +-------------------------------------------
    Transverse++   |    NO        YES, at predictor step
                   |
    Normal++       |    NO        NO
  \endverbatim
 
  If constrained transport is enabled, then
  \verbatim
                        RK        CTU
                   +-------------------------------------------
    Transverse++   |    YES       YES
                   |
    Normal++       |    NO        YES, at predictor step
  \endverbatim
  
  Predictor step (g_intStage = 1) in CTU requires extra transverse loop
  to obtain transverse predictor in any case.
  Also, with CTU+CT, one needs to expand the grid of one zone in the
  \e normal direction as well.
  This allows to computed fully corner coupled sweeps in the boundary to
  get electric field components during the constrained transport algorithm. 

  \author A. Mignone (mignone@ph.unito.it)\n
  \date   Nov 27, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void SetVectorIndices (int dir)
/*!
 * Set vector indices and integration index range.
 *
 * \param [out] indx pointer to an Index structure
 * \param [in]  grid pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  if (dir == IDIR) {   /* -- Order: X-Y-Z  {in,t1,t2 = i,j,k} -- */

    EXPAND(VXn = MXn = VX1; , 
           VXt = MXt = VX2; , 
           VXb = MXb = VX3;)
    #if PHYSICS == MHD || PHYSICS == RMHD
    EXPAND(BXn = BX1; , 
           BXt = BX2; ,  
           BXb = BX3;)
    #endif
    #if (PHYSICS == RMHD) && (RESISTIVITY != NO)
    EXPAND(EXn = EX1; , 
           EXt = EX2; ,  
           EXb = EX3;)
    #endif

#if DUST_FLUID == YES
    EXPAND(VXn_D = MXn_D = VX1_D; , 
           VXt_D = MXt_D = VX2_D; ,  
           VXb_D = MXb_D = VX3_D;)
#endif
 
  }else if (dir == JDIR){ /* -- Order: Y-X-Z  {in,t1,t2 = j,i,k} -- */

    EXPAND(VXn = MXn = VX2;  , 
           VXt = MXt = VX1;  , 
           VXb = MXb = VX3;)
    #if PHYSICS == MHD || PHYSICS == RMHD
    EXPAND(BXn = BX2; , 
           BXt = BX1; ,  
           BXb = BX3;)
    #endif
    #if (PHYSICS == RMHD) && (RESISTIVITY != NO)
    EXPAND(EXn = EX2; , 
           EXt = EX1; ,  
           EXb = EX3;)
    #endif
#if DUST_FLUID == YES
    EXPAND(VXn_D = MXn_D = VX2_D; , 
           VXt_D = MXt_D = VX1_D; ,  
           VXb_D = MXb_D = VX3_D;)
#endif
    
  }else if (dir == KDIR){ /* -- Order: Z-X-Y  {in,t1,t2 = k,i,j} -- */

    VXn = MXn = VX3;
    VXt = MXt = VX1;
    VXb = MXb = VX2;
    #if PHYSICS == MHD || PHYSICS == RMHD
     BXn = BX3;
     BXt = BX1;
     BXb = BX2;
    #endif
    #if (PHYSICS == RMHD) && (RESISTIVITY != NO)
    EXPAND(EXn = EX3; , 
           EXt = EX1; ,  
           EXb = EX2;)
    #endif
#if DUST_FLUID == YES
    EXPAND(VXn_D = MXn_D = VX3_D; , 
           VXt_D = MXt_D = VX1_D; ,  
           VXb_D = MXb_D = VX2_D;)
#endif

  }
}

/* ********************************************************************* */
void ResetState (const Data *d, Sweep *sweep, Grid *grid)
/*!
 * Initialize some of the elements of the Sweep structure to zero
 * in order to speed up computations. These includes:
 *
 *    - source term
 *    - left and right eigenvectors
 *    - the maximum eigenvalue ???
 *
 * \param [in] d  pointer to Data structure
 * \param [out] sweep pointer to a Sweep structure
 * \param [in]  grid pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int i, j, k, nv;
  double v[NVAR], lambda[NVAR];
  double a;

  memset ((void *)sweep->src[0],   '\0', NMAX_POINT*NVAR*sizeof(double));

  memset ((void *)sweep->stateC.Rp[0][0], '\0', NMAX_POINT*NFLX*NFLX*sizeof(double));
  memset ((void *)sweep->stateC.Lp[0][0], '\0', NMAX_POINT*NFLX*NFLX*sizeof(double));
  memset ((void *)sweep->stateL.Rp[0][0], '\0', NMAX_POINT*NFLX*NFLX*sizeof(double));
  memset ((void *)sweep->stateL.Lp[0][0], '\0', NMAX_POINT*NFLX*NFLX*sizeof(double));
  memset ((void *)sweep->stateR.Rp[-1][0], '\0', NMAX_POINT*NFLX*NFLX*sizeof(double));
  memset ((void *)sweep->stateR.Lp[-1][0], '\0', NMAX_POINT*NFLX*NFLX*sizeof(double));
/*
  for (i = 0; i < NMAX_POINT; i++){
    for (j = NVAR; j--;  ) sweep->src[i][j] = 0.0;

    for (j = NFLX; j--;  ){
    for (k = NFLX; k--;  ){
      sweep->Lp[i][j][k] = sweep->Rp[i][j][k] = 0.0;
    }}
  }  
*/
/* ---------------------------------------------------------
     When using Finite Difference methods, we need to find,
     for each characteristic k, its maximum over the 
     whole computational domain (LF global splitting).
   --------------------------------------------------------- */

  #ifdef FINITE_DIFFERENCE
   FD_GetMaxEigenvalues (d, sweep, grid);
  #endif
}
