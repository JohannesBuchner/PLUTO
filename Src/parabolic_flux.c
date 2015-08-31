/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute diffusion fluxes for explicit time stepping.

  Compute parabolic fluxes and corresponding source terms for explicit 
  time stepping only and add them to upwind fluxes.
  For ::STS and ::RKC integration see the ParabolicRHS() function.
  Note that source terms are only included for viscous terms, since 
  resistivity and thermal conduction do not need any.
  
  \note For explicit resistive MHD, the EMF is comprised of 2 terms:
      EMF = E(hyp) + E(par)
   The correct sequence of steps for building the EMF are:
   - E(hyp) is the flux computed with Riemann solver 
   - for STAGGERED_MHD E(hyp) is stored at appropriate 
     location for later re-use by calling ::CT_StoreEMF
   - Parabolic fluxes are added to the hyperbolic flux 
     (useful only for cell-centered MHD).
    
   
  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date    Sept 4, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if THERMAL_CONDUCTION != NO && (EOS == ISOTHERMAL || EOS == BAROTROPIC) 
 #error ! No Energy Equation: Thermal Conduction cannot be included
#endif

/* ********************************************************************* */
void ParabolicFlux (Data_Arr V, Data_Arr J, double ***T,
                    const State_1D *state,
                    double **dcoeff, int beg, int end, Grid *grid)
/*! 
 * Add the diffusion fluxes to the upwind fluxes for explicit time
 * integration.
 *
 * \param [in] V   pointer to the 3D array of cell-centered primitive
 *                 variables
 * \param [in,out] state pointer to a State_1D structure
 * \param [out]    dcoeff the diffusion coefficients 1D array
 * \param [in]      beg   initial index of computation
 * \param [in]      end   final   index of computation
 * \param [in]      grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int i, j, k, nv;
  
/* -------------------------------------------------
   1. Viscosity 
   ------------------------------------------------- */
    
  #if VISCOSITY == EXPLICIT
   #ifdef FARGO
    print ("! ParabolicFlux: FARGO incompatible with explicit viscosity.\n");
    print ("                 Try STS or RKC instead\n");
    QUIT_PLUTO(1);
   #endif
   ViscousFlux (V, state->visc_flux, state->visc_src, dcoeff, beg, end, grid);
   for (i = beg; i <= end; i++){
     EXPAND(state->flux[i][MX1] -= state->visc_flux[i][MX1];  ,
            state->flux[i][MX2] -= state->visc_flux[i][MX2];  ,
            state->flux[i][MX3] -= state->visc_flux[i][MX3];  ) 
     #if HAVE_ENERGY
      state->flux[i][ENG] -= state->visc_flux[i][ENG];
     #endif

   }
  #endif 

/* -------------------------------------------------
   2. Thermal conduction
   ------------------------------------------------- */

  #if THERMAL_CONDUCTION == EXPLICIT
   TC_Flux (T, state, dcoeff, beg, end, grid);  
   for (i = beg; i <= end; i++) state->flux[i][ENG] -= state->tc_flux[i][ENG];
   /* !!!! tc_flux can be redefined as a 1D array since only the energy
           component is required !!! */
  #endif 

/* ----------------------------------------------------------------
   3. Resistivity.
      Note for the entropy flux: differently from viscosity and
      thermal conduction, the contribution to the entropy due to
      currents is given by J^2 and cannot be expressed as a simple
      two-point flux difference.
      This is done in a separate step (...)
   ---------------------------------------------------------------- */

  #if (RESISTIVITY == EXPLICIT) 
   ResistiveFlux (V, J, state->res_flux, dcoeff, beg, end, grid);
   for (i = beg; i <= end; i++){

   /* ------------------------------------------
       normal component of magnetic field does 
       not evolve during the current sweep. 
      ------------------------------------------ */

     state->res_flux[i][BXn] = 0.0;  

  /* ---------------------------------------------------
      add the parabolic part of the EMF, only for 
      cell-centered MHD. CT is handled in a truly
      multidimensional way.
     --------------------------------------------- */

     EXPAND(state->flux[i][BX1] += state->res_flux[i][BX1];  ,
            state->flux[i][BX2] += state->res_flux[i][BX2];  ,
            state->flux[i][BX3] += state->res_flux[i][BX3]; )
     #if HAVE_ENERGY
      state->flux[i][ENG] += state->res_flux[i][ENG];
     #endif
   }
  #endif 

}
