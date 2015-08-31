/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute flux for passive scalars.                            

  Compute the interface upwind flux for passive scalars \c q obeying 
  advection equations of the form:
  \f[ \partial_tq + v\cdot \nabla q = 0 
      \qquad\Longleftrightarrow\qquad
      \partial_t(\rho q) + \nabla\cdot(\rho\vec{v}q) = 0 
  \f]  
  Fluxes are computed using an upwind selection rule based on the 
  density flux, already computed during a previous Riemann solver:
  \f[ 
    (\rho vq)_{i+\HALF} = \left\{\begin{array}{ll}
    (\rho v)_{i+\HALF}q_L & \;\textrm{if} 
           \quad (\rho v)_{i+\HALF} \ge 0 \\ \noalign{\medskip}
    (\rho v)_{i+\HALF}q_R & \; \textrm{otherwise}  
   \end{array}\right.
  \f]
  where \f$ (\rho v)_{i+\HALF}\f$ is the density flux computed with 
  the employed Riemann solver.
  
  When ionization fractions are present, we employ a technique similar
  to the CMA (Consistent multi-fluid advection method) to normalize
  the sum of mass fractions to one.
  
  The CMA can also be switched on for standard tracers (<tt> 
  #define USE_CMA  YES</tt>)
  
  \author A. Mignone (mignone@ph.unito.it)\n
          O. Tesileanu
  \date   02 April, 2015

  \b Reference\n
     "The consistent multi-fluid advection method"
      Plewa and Muller, A&A (1999) 342, 179.
*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#ifndef USE_CMA
  #define USE_CMA NO
#endif

/* ********************************************************************* */
void AdvectFlux (const State_1D *state, int beg, int end, Grid *grid)
/*! 
 *
 * \param [in,out] state
 * \param [in]      beg    initial index of computation 
 * \param [in]      end    final   index of computation
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int    i, nv;
  double *ts, *flux, *vL, *vR;
  double s, rho;
  double phi;
  static double *sigma, **vi;
  
/* -- compute scalar's fluxes -- */

  for (i = beg; i <= end; i++){

    flux = state->flux[i];
    vL   = state->vL[i];
    vR   = state->vR[i];

    ts = flux[RHO] > 0.0 ? vL:vR;

    NSCL_LOOP(nv) flux[nv] = flux[RHO]*ts[nv];
    
    #if COOLING == MINEq

     /* -----   He   -----  */
                         
     phi = ts[X_HeI] + ts[X_HeII];
     for (nv = X_HeI; nv <= X_HeII; nv++) flux[nv] /= phi;

     /* -----   C   -----  */

     phi = 0.0;
     for (nv = X_CI; nv < X_CI + C_IONS; nv++) phi += ts[nv]; 
     for (nv = X_CI; nv < X_CI + C_IONS; nv++) flux[nv] /= phi;

     /* -----   N   -----  */

     phi = 0.0;
     for (nv = X_NI; nv < X_NI+N_IONS; nv++) phi += ts[nv]; 
     for (nv = X_NI; nv < X_NI+N_IONS; nv++) flux[nv] /= phi;

     /* -----   O   -----  */

     phi = 0.0;
     for (nv = X_OI; nv < X_OI+O_IONS; nv++) phi += ts[nv]; 
     for (nv = X_OI; nv < X_OI+O_IONS; nv++) flux[nv] /= phi;

     /* -----   Ne   -----  */

     phi = 0.0;
     for (nv = X_NeI; nv < X_NeI+Ne_IONS; nv++) phi += ts[nv];
     for (nv = X_NeI; nv < X_NeI+Ne_IONS; nv++) flux[nv] /= phi;

     /* -----   S   -----  */

     phi = 0.0;
     for (nv = X_SI; nv < X_SI+S_IONS; nv++) phi += ts[nv]; 
     for (nv = X_SI; nv < X_SI+S_IONS; nv++) flux[nv] /= phi;

    #endif

    #if COOLING == H2_COOL
     phi = ts[X_HI] + 2.0*ts[X_H2] + ts[X_HII];
     for (nv = X_HI; nv < X_HI + NIONS; nv++) flux[nv] /= phi; 
    #endif

    #if USE_CMA == YES  /* -- only for tracers -- */
     phi = 0.0;
     NTRACER_LOOP(nv) phi += ts[nv];
     NTRACER_LOOP(nv) flux[nv] /= phi;
    #endif

    #if ENTROPY_SWITCH
     if (flux[RHO] >= 0.0) flux[ENTR] = state->vL[i][ENTR]*flux[RHO];
     else                  flux[ENTR] = state->vR[i][ENTR]*flux[RHO];
    #endif
  }
}
