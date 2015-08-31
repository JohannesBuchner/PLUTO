/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Include source terms using operator splitting.

  The SplitSource() function handles source terms in a separate
  step using operator splitting.
  It is called from Integrate() between standard hydro advances.
  At present these source terms are one or more of the following:

  - optically thin radiative losses (cooling)
  - Diffusion operators: 
    - resistivity 
    - Thermal conduction
    - Viscosity
  - additional user-defined terms may also be included here.

  \authors A. Mignone (mignone@ph.unito.it)
  \date    May 10, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void SplitSource (const Data *d, double dt, Time_Step *Dts, Grid *grid)
/*! 
 *  Take one step on operator-split source terms.
 *
 *  \param [in,out] d   pointer to PLUTO Data structure containing 
 *                      the solution array updated from the most 
 *                      recent call
 *  \param[in]      dt  the time step used to integrate the source 
 *                      terms
 *  \param[in]     Dts  pointer to the time step structure
 *  \param[in]    grid  pointer to an array of grid structures
 *
 *********************************************************************** */
{
/*  ---- GLM source term treated in main ----  */
/*
  #ifdef GLM_MHD
   GLM_SOURCE (d->Vc, dt, grid);
  #endif
*/
/*  ---------------------------------------------
             Cooling/Heating losses
    ---------------------------------------------  */

  #if COOLING != NO
   #if COOLING == POWER_LAW  /* -- solve exactly -- */
    PowerLawCooling (d->Vc, dt, Dts, grid);
   #else
    CoolingSource (d, dt, Dts, grid);
   #endif
  #endif

/* ----------------------------------------------
    Parabolic terms using STS:

    - resistivity 
    - thermal conduction
    - viscosity 
   ---------------------------------------------- */

  #if (PARABOLIC_FLUX & SUPER_TIME_STEPPING)
   STS (d, Dts, grid);
  #endif

  #if (PARABOLIC_FLUX & RK_CHEBYSHEV)
   RKC (d, Dts, grid);
  #endif
}
