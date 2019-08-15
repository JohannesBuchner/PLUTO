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
  
  \author  A. Mignone (mignone@ph.unito.it)\n
           O. Tesileanu
  \date    Feb 28, 2017

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
void AdvectFlux (const Sweep *sweep, int beg, int end, Grid *grid)
/*! 
 *
 * \param [in,out] sweep
 * \param [in]      beg    initial index of computation 
 * \param [in]      end    final   index of computation
 * \param [in]      grid   Pointer to Grid structure
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

    flux = sweep->flux[i];
    vL   = sweep->stateL.v[i];
    vR   = sweep->stateR.v[i];

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

    #if COOLING == KROME
     phi = 0.0;
     for (nv = X_H; nv < X_H + NIONS; nv++) phi += ts[nv];
     for (nv = X_H; nv < X_H + NIONS; nv++) flux[nv] /= phi;
    #endif

    #if USE_CMA == YES  /* -- only for tracers -- */
     phi = 0.0;
     NTRACER_LOOP(nv) phi += ts[nv];
     NTRACER_LOOP(nv) flux[nv] /= phi;
    #endif

    #if ENTROPY_SWITCH
    if (flux[RHO] >= 0.0) flux[ENTR] = vL[ENTR]*flux[RHO];
    else                  flux[ENTR] = vR[ENTR]*flux[RHO];
    #endif
  }
}

/* ********************************************************************* */
void StoreAMRFlux (double **flux, double **aflux, int sign,
                    int nvar_beg, int nvar_end, int beg, int end, Grid *grid)
/*!
 * \param [in]  flux      pointer to a 1D flux array
 * \param [out] aflux     pointer to a 1D flux array for AMR refluxing
 *                        operation
 * \param [in]  sign      an integer equal to 0, +1 or -1.
 *                        When equal to 0, flux is initialized, otherwise
 *                        it is added (+1) or subtracted (-1)
 * \param [in]  nvar_beg  the starting variable index
 * \param [in]  nvar_end  the final variable index
 * \param [in]  grid      a pointer to the grid structure.
 *
 *********************************************************************** */
{
  int i,j,k,nv;
  int nxf, nyf, nzf;
  int nxb, nyb, nzb;
  int *in;
  long int indf, ind1;
  double w;  

#ifdef CTU                      /* With CTU, fluxes are saved at the */
  if (g_intStage == 1) return;  /* corrector step.                   */ 
#endif

  nxf = grid->np_int[IDIR] + (g_dir == IDIR);
  nyf = grid->np_int[JDIR] + (g_dir == JDIR);
  nzf = grid->np_int[KDIR] + (g_dir == KDIR);

  nxb = grid->lbeg[IDIR] - (g_dir == IDIR);
  nyb = grid->lbeg[JDIR] - (g_dir == JDIR);
  nzb = grid->lbeg[KDIR] - (g_dir == KDIR);

  i = g_i; j = g_j; k = g_k;
  if (g_dir == IDIR) in = &i;
  if (g_dir == JDIR) in = &j;
  if (g_dir == KDIR) in = &k;

#if TIME_STEPPING == RK2
  w = 0.5;
#else 
  w = 1.0;
#endif

  if (sign == 0){

    for ((*in) = beg; (*in) <= end; (*in)++) {
      ind1 = (k - nzb)*nyf*nxf + (j - nyb)*nxf + (i - nxb);
      for (nv = nvar_beg; nv <= nvar_end; nv++){
        indf = nv*nzf*nyf*nxf + ind1;
        aflux[g_dir][indf] = w*flux[*in][nv];
      }
    }

  }else{ 

    for ((*in) = beg; (*in) <= end; (*in)++) {
      ind1 = (k - nzb)*nyf*nxf + (j - nyb)*nxf + (i - nxb);
      for (nv = nvar_beg; nv <= nvar_end; nv++){
        indf = nv*nzf*nyf*nxf + ind1;
        aflux[g_dir][indf] += sign*w*flux[*in][nv];
      }
    }

  }
}
