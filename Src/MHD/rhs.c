/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute the right hand side of the conservative 
         HD/MHD equations.

  This function constructs the one-dimensional right hand side of 
  the conservative MHD or HD equations in the direction given by ::g_dir 
  in different geometries.
  The right hand side in the \c d direction is computed as a two-point
  flux difference term plus a source term:  
  \f[
      {\cal R}_{\ivec}^{(d)} = 
      -\frac{\Delta t}{\Delta{\cal V}_{\ivec}}\
      \Big[  (A{\cal F})_{\ivec+\HALF\hvec{e}_d} 
            -(A{\cal F})_{\ivec-\HALF\hvec{e}_d} \Big]
         + \Delta t{\cal S}_{\ivec}
  \f] 
   where \f$\ivec = (i,j,k)\f$ while   
   - \f$ A        \f$: interface areas
   - \f$ {\cal V} \f$: cell volume
   - \f$ \F   \f$: interface fluxes
   - \f$ \Delta t \f$: time step
   - \f$ {\cal S} \f$: source term including geometrical terms and 
                       body forces.

  See also the \ref RHS_page.
  The right hand side is assembled through the following steps:
 
  - If either one of FARGO, ROTATION or gravitational potential is used, 
    fluxes are combined to enforce conservation of total angular momentum
    and/or energy, see TotalFlux();
  - initialize rhs with flux differences;
  - add dissipative effects (viscosity and thermal conduction) to
    entropy equation. Ohmic dissipation is included in a separate step.

  Source terms are added later in RightHandSideSource().

  For the entropy equation, the dissipative contributions can be recovered
  from the internal energy equation, for which (Boyd, Eq. [3.27]): 
  \f[
      \pd{p}{t} + \vec{v}\cdot\nabla p + \Gamma p \nabla\cdot\vec{v}
    = (\Gamma-1)\left[\nabla\cdot\left(\kappa\nabla T\right)
      + \eta\vec{J}\cdot\vec{J} + \mu\left(\pd{v_i}{x_j}+\pd{v_j}{x_i}
        - \frac{2}{3}\delta_{ij}\nabla\cdot\vec{v}\right)\pd{v_i}{x_j}\right]
  \f] 
  To obtain the corresponding contribution to the (conservative form of)
  entropy equation, just divide the previous one by 
  \f$\rho^{\Gamma-1}\f$: 
  \f[
     \pd{(\rho s)}{t} + \nabla\cdot(\rho s\vec{v}) = 
     (\Gamma-1)\rho^{1-\Gamma}\left[\nabla\cdot\left(\kappa\nabla T\right)
       + \eta\vec{J}^2 + \mu \Pi_{ij}\pd{v_i}{x_j}\right]
  \f]
  See the books by
  - Goedbloed & Poedts, "Principle of MHD", page 165-166
  - Boyd, "The Physics of Plasmas", page 57, Eq. 3.27

  \b References
     - "PLUTO: A Numerical Code for Computational Astrophysics."
       Mignone et al, ApJS (2007) 170, 228
     - "The PLUTO Code for Adaptive Mesh Computations in Astrophysical Fluid
        Dynamics" 
       Mignone et al, ApJS (2012) 198, 7M
   
  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 20, 2015

  \todo
    - merge GetAreaFlux and TotalFlux ?
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static void TotalFlux (const State_1D *, double *, int, int, Grid *);

#ifdef CH_SPACEDIM   /*  implies Chombo is being used  */
 #define USE_PR_GRADIENT  YES   
#else
 #ifdef FINITE_DIFFERENCE 
  #define USE_PR_GRADIENT  NO   /* -- default for Finite Difference schemes -- */
 #else
  #define USE_PR_GRADIENT  YES   /* -- default, do not change!! -- */
 #endif
#endif

/* *********************************************************************** */
void RightHandSide (const State_1D *state, Time_Step *Dts, 
                    int beg, int end, double dt, Grid *grid)
/*! 
 *
 * \param [in,out]  state  pointer to State_1D structure
 * \param [in]      Dts    pointer to time step structure
 * \param [in]      beg    initial index of computation
 * \param [in]      end    final   index of computation
 * \param [in]      dt     time increment
 * \param [in]      grid  pointer to Grid structure
 *
 * \return This function has no return value.
 * \note    --
 * \todo    --
 ************************************************************************* */
{
  int    i, j, k, nv;
  double dtdx, dtdV, scrh, rhog;
  double *x1  = grid[IDIR].x,  *x2  = grid[JDIR].x,  *x3  = grid[KDIR].x;
  double *x1p = grid[IDIR].xr, *x2p = grid[JDIR].xr, *x3p = grid[KDIR].xr;
  double *dx1 = grid[IDIR].dx, *dx2 = grid[JDIR].dx, *dx3 = grid[KDIR].dx;
  double *dV1 = grid[IDIR].dV, *dV2 = grid[JDIR].dV, *dV3 = grid[KDIR].dV;
  double **rhs  = state->rhs;
  double **flux = state->flux, *p = state->press;
  double **vh   = state->vh, **vp = state->vp, **vm = state->vm;
  double *A     = grid[g_dir].A, *dV    = grid[g_dir].dV;
  
  double cl;
  double w, wp, vphi, phi_c;
  double g[3];
  static double **fA, *phi_p;
  static double **fvA;
#if ENTROPY_SWITCH
  double rhs_entr;
  double **visc_flux = state->visc_flux; 
  double **visc_src  = state->visc_src; 
  double **tc_flux   = state->tc_flux; 
  double **res_flux  = state->res_flux;   
  if (fvA == NULL) fvA = ARRAY_2D(NMAX_POINT, NVAR, double);
#endif

  #if GEOMETRY != CARTESIAN
   if (fA == NULL) fA = ARRAY_2D(NMAX_POINT, NVAR, double);
  #endif
  if (phi_p == NULL) phi_p = ARRAY_1D(NMAX_POINT, double);

/* --------------------------------------------------
   1. Compute fluxes for passive scalars and dust
   -------------------------------------------------- */

#if NSCL > 0
  AdvectFlux (state, beg - 1, end, grid);
#endif

#if DUST == YES
  Dust_Solver(state, beg - 1, end, Dts->cmax, grid);
#endif

  i = g_i;  /* will be redefined during x1-sweep */
  j = g_j;  /* will be redefined during x2-sweep */
  k = g_k;  /* will be redefined during x3-sweep */
  
/* ------------------------------------------------
     Add pressure to normal component of 
     momentum flux if necessary.
   ------------------------------------------------ */

#if USE_PR_GRADIENT == NO
  for (i = beg - 1; i <= end; i++) flux[i][MXn] += p[i];
#endif

/* -----------------------------------------------------
     Compute gravitational potential at cell interfaces
   -----------------------------------------------------  */

#if (BODY_FORCE & POTENTIAL)
  if (g_dir == IDIR) {
    for (i = beg-1; i <= end; i++) {
      phi_p[i] = BodyForcePotential(x1p[i], x2[j], x3[k]);
    }
  }else if (g_dir == JDIR){
    for (j = beg-1; j <= end; j++) {
      phi_p[j] = BodyForcePotential(x1[i], x2p[j], x3[k]);
    }
  }else if (g_dir == KDIR){
    for (k = beg-1; k <= end; k++) {
      phi_p[k] = BodyForcePotential(x1[i], x2[j], x3p[k]);
    }
  }
#endif
  
/* -----------------------------------------------------
    Compute total flux
   -----------------------------------------------------  */

  #if (defined FARGO && !defined SHEARINGBOX) ||\
      (ROTATING_FRAME == YES) || (BODY_FORCE & POTENTIAL)
   TotalFlux(state, phi_p, beg-1, end, grid);
  #endif

#if GEOMETRY == CARTESIAN

  if (g_dir == IDIR){

    for (i = beg; i <= end; i++) {
      dtdx = dt/dx1[i];

    /* -- I1. initialize rhs with flux difference -- */

      NVAR_LOOP(nv) rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i-1][nv]);
      #if USE_PR_GRADIENT == YES
       rhs[i][MX1] -= dtdx*(p[i] - p[i-1]);
      #endif

    /* -- I5. Add dissipative terms to entropy equation -- */

      #if (ENTROPY_SWITCH) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vh[i][RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (visc_flux[i][ENG] - visc_flux[i-1][ENG]);
        rhs_entr -= EXPAND(  vh[i][VX1]*(visc_flux[i][MX1] - visc_flux[i-1][MX1])  ,
                           + vh[i][VX2]*(visc_flux[i][MX2] - visc_flux[i-1][MX2])  ,
                           + vh[i][VX3]*(visc_flux[i][MX3] - visc_flux[i-1][MX3]));
       #endif

       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (tc_flux[i][ENG] - tc_flux[i-1][ENG]);       
       #endif
       rhs[i][ENTR] += rhs_entr*dtdx*rhog;              
      #endif

    }
  } else if (g_dir == JDIR){

    for (j = beg; j <= end; j++) {
      dtdx = dt/dx2[j];

    /* -- J1. initialize rhs with flux difference -- */

      NVAR_LOOP(nv) rhs[j][nv] = -dtdx*(flux[j][nv] - flux[j-1][nv]);
      #if USE_PR_GRADIENT == YES
       rhs[j][MX2] -= dtdx*(p[j] - p[j-1]);
      #endif

    /* -- J5. Add dissipative terms to entropy equation -- */

      #if (ENTROPY_SWITCH) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vh[j][RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (visc_flux[j][ENG] - visc_flux[j-1][ENG]);
        rhs_entr -= EXPAND(  vh[j][VX1]*(visc_flux[j][MX1] - visc_flux[j-1][MX1])  ,
                           + vh[j][VX2]*(visc_flux[j][MX2] - visc_flux[j-1][MX2])  ,
                           + vh[j][VX3]*(visc_flux[j][MX3] - visc_flux[j-1][MX3]));
       #endif

       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (tc_flux[j][ENG] - tc_flux[j-1][ENG]);       
       #endif
       rhs[j][ENTR] += rhs_entr*dtdx*rhog;              
      #endif
      
    }

  }else if (g_dir == KDIR){

    for (k = beg; k <= end; k++) {
      dtdx = dt/dx3[k];

    /* -- K1. initialize rhs with flux difference -- */

      NVAR_LOOP(nv) rhs[k][nv] = -dtdx*(flux[k][nv] - flux[k-1][nv]);
      #if USE_PR_GRADIENT == YES
       rhs[k][MX3] -= dtdx*(p[k] - p[k-1]);
      #endif

    /* -- K5. Add dissipative terms to entropy equation -- */

      #if (ENTROPY_SWITCH) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vh[k][RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma); 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (visc_flux[k][ENG] - visc_flux[k-1][ENG]);
        rhs_entr -= EXPAND(  vh[k][VX1]*(visc_flux[k][MX1] - visc_flux[k-1][MX1])  ,
                           + vh[k][VX2]*(visc_flux[k][MX2] - visc_flux[k-1][MX2])  ,
                           + vh[k][VX3]*(visc_flux[k][MX3] - visc_flux[k-1][MX3]));
       #endif

       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (tc_flux[k][ENG] - tc_flux[k-1][ENG]);
       #endif
       rhs[k][ENTR] += rhs_entr*dtdx*rhog;              
      #endif
    }
  }

#elif GEOMETRY == CYLINDRICAL

{
  double R, z, phi, R_1; 

#if DUST == YES
  #error "DUST not implemented in CYLINDRICAL coordinates"
#endif

  if (g_dir == IDIR) {  
    double vc[NVAR];

    GetAreaFlux (state, fA, fvA, beg, end, grid);
    for (i = beg; i <= end; i++){ 
      R    = x1[i];
      dtdV = dt/dV1[i];
      dtdx = dt/dx1[i];
      R_1  = 1.0/R;

    /* -- I1. initialize rhs with flux difference -- */

      rhs[i][RHO] = -dtdV*(fA[i][RHO] - fA[i-1][RHO]);
      EXPAND(rhs[i][iMR]   = - dtdV*(fA[i][iMR]   - fA[i-1][iMR]);  ,
             rhs[i][iMZ]   = - dtdV*(fA[i][iMZ]   - fA[i-1][iMZ]);  ,
             rhs[i][iMPHI] = - dtdV*(fA[i][iMPHI] - fA[i-1][iMPHI])*fabs(R_1);)
      #if USE_PR_GRADIENT == YES
       rhs[i][iMR] -= dtdx*(p[i] - p[i-1]);  
      #endif
      #if PHYSICS == MHD
       EXPAND(rhs[i][iBR]   = - dtdV*(fA[i][iBR]   - fA[i-1][iBR]);  ,
              rhs[i][iBZ]   = - dtdV*(fA[i][iBZ]   - fA[i-1][iBZ]);  ,
              rhs[i][iBPHI] = - dtdx*(flux[i][iBPHI] - flux[i-1][iBPHI]);)
       #ifdef GLM_MHD
        rhs[i][iBR]     = - dtdx*(flux[i][iBR]   - flux[i-1][iBR]);
        rhs[i][PSI_GLM] = - dtdV*(fA[i][PSI_GLM] - fA[i-1][PSI_GLM]);
       #endif
      #endif
      IF_ENERGY(rhs[i][ENG] = -dtdV*(fA[i][ENG] - fA[i-1][ENG]);)
      NSCL_LOOP(nv)  rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);
      
    /* -- I5. Add dissipative terms to entropy equation -- */

      #if (ENTROPY_SWITCH) && (PARABOLIC_FLUX & EXPLICIT)
       NVAR_LOOP(nv) vc[nv] = 0.5*(vp[i][nv] + vm[i][nv]); 
       rhog = vc[RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (fvA[i][ENG] - fvA[i-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(fvA[i][MX1] - fvA[i-1][MX1])      ,
                           + vc[VX2]*(fvA[i][MX2] - fvA[i-1][MX2])      ,
                           + vc[VX3]*(fvA[i][MX3] - fvA[i-1][MX3])*fabs(R_1));
                    
        rhs_entr -= EXPAND(  vc[VX1]*visc_src[i][MX1], 
                           + vc[VX2]*visc_src[i][MX2],
                           + vc[VX3]*visc_src[i][MX3]);
       #endif

       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (A[i]*tc_flux[i][ENG] - A[i-1]*tc_flux[i-1][ENG]);       
       #endif
       rhs[i][ENTR] += rhs_entr*dtdV*rhog;              
      #endif
    }
     
  } else if (g_dir == JDIR) { 

    for (j = beg; j <= end; j++){ 
      dtdx = dt/dx2[j];

    /* -- J1. initialize rhs with flux difference -- */

      NVAR_LOOP(nv) rhs[j][nv] = -dtdx*(flux[j][nv] - flux[j-1][nv]);
      #if USE_PR_GRADIENT == YES
       rhs[j][iMZ] += - dtdx*(p[j] - p[j-1]);
      #endif

    /* -- J5. Add dissipative terms to entropy equation -- */

      #if (ENTROPY_SWITCH) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vh[j][RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (visc_flux[j][ENG] - visc_flux[j-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(visc_flux[j][MX1] - visc_flux[j-1][MX1])  ,
                           + vc[VX2]*(visc_flux[j][MX2] - visc_flux[j-1][MX2])  ,
                           + vc[VX3]*(visc_flux[j][MX3] - visc_flux[j-1][MX3]));
       #endif
       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (tc_flux[j][ENG] - tc_flux[j-1][ENG]);
       #endif
       rhs[j][ENTR] += rhs_entr*dtdx*rhog;              
      #endif
    }
  }
}

#elif GEOMETRY == POLAR
{
  double R, phi, z; 
  double R_1;
   
  if (g_dir == IDIR) { 
    double vc[NVAR];

    GetAreaFlux (state, fA, fvA, beg, end, grid);
    for (i = beg; i <= end; i++) {
      R    = x1[i];
      dtdV = dt/dV1[i];
      dtdx = dt/dx1[i];
      R_1  = grid[IDIR].r_1[i];

    /* -- I1. initialize rhs with flux difference -- */

      rhs[i][RHO] = -dtdV*(fA[i][RHO] - fA[i-1][RHO]);
      EXPAND(rhs[i][iMR]   = - dtdV*(fA[i][iMR]   - fA[i-1][iMR])
                             - dtdx*(p[i] - p[i-1]);                      ,      
             rhs[i][iMPHI] = - dtdV*(fA[i][iMPHI] - fA[i-1][iMPHI])*R_1;  ,
             rhs[i][iMZ]   = - dtdV*(fA[i][iMZ]   - fA[i-1][iMZ]);)
  #if PHYSICS == MHD 
      EXPAND(rhs[i][iBR]   = -dtdV*(fA[i][iBR]     - fA[i-1][iBR]);      ,
             rhs[i][iBPHI] = -dtdx*(flux[i][iBPHI] - flux[i-1][iBPHI]);  ,
             rhs[i][iBZ]   = -dtdV*(fA[i][iBZ]     - fA[i-1][iBZ]);)
    #ifdef GLM_MHD
      rhs[i][iBR]     = -dtdx*(flux[i][iBR]   - flux[i-1][iBR]);
      rhs[i][PSI_GLM] = -dtdV*(fA[i][PSI_GLM] - fA[i-1][PSI_GLM]);
    #endif
  #endif
      IF_ENERGY(rhs[i][ENG] = -dtdV*(fA[i][ENG] - fA[i-1][ENG]);)

      NSCL_LOOP(nv)  rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);

      IF_DUST(NDUST_LOOP(nv) rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);  
              rhs[i][MX2_D] *= R_1;)

    /* -- I5. Add dissipative terms to entropy equation -- */

      #if (ENTROPY_SWITCH) && (PARABOLIC_FLUX & EXPLICIT)
       for (nv = NVAR; nv--;  ) vc[nv] = 0.5*(vp[i][nv] + vm[i][nv]); 
       rhog = vc[RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (fvA[i][ENG] - fvA[i-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(fvA[i][iMR]   - fvA[i-1][iMR])        ,
                           + vc[VX2]*(fvA[i][iMPHI] - fvA[i-1][iMPHI])*R_1  ,
                           + vc[VX3]*(fvA[i][iMZ]   - fvA[i-1][iMZ]));
                    
        rhs_entr -= EXPAND(  vc[VX1]*visc_src[i][MX1], 
                           + vc[VX2]*visc_src[i][MX2],
                           + vc[VX3]*visc_src[i][MX3]);
       #endif

       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (A[i]*tc_flux[i][ENG] - A[i-1]*tc_flux[i-1][ENG]);       
       #endif
       rhs[i][ENTR] += rhs_entr*dtdV*rhog;              
      #endif

    }
     
  } else if (g_dir == JDIR) {

    scrh = dt/x1[i];
    for (j = beg; j <= end; j++){ 
      dtdx = scrh/dx2[j];

    /* -- J1. Compute equations rhs for phi-contributions -- */

      NVAR_LOOP(nv) rhs[j][nv] = -dtdx*(flux[j][nv] - flux[j-1][nv]);
      rhs[j][iMPHI] -= dtdx*(p[j] - p[j-1]);

    /* -- J5. Add dissipative terms to entropy equation -- */

      #if (ENTROPY_SWITCH) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vh[j][RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (visc_flux[j][ENG] - visc_flux[j-1][ENG]);
        rhs_entr -= EXPAND(  vh[j][VX1]*(visc_flux[j][MX1] - visc_flux[j-1][MX1]) ,
                           + vh[j][VX2]*(visc_flux[j][MX2] - visc_flux[j-1][MX2]) ,
                           + vh[j][VX3]*(visc_flux[j][MX3] - visc_flux[j-1][MX3]));
                    
        rhs_entr -= EXPAND(  vh[j][VX1]*visc_src[j][MX1], 
                           + vh[j][VX2]*visc_src[j][MX2],
                           + vh[j][VX3]*visc_src[j][MX3]);
       #endif
       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (tc_flux[j][ENG] - tc_flux[j-1][ENG]);
       #endif
       rhs[j][ENTR] += rhs_entr*dtdx*rhog;              
      #endif
      
    }

  } else if (g_dir == KDIR) { 

    for (k = beg; k <= end; k++){ 
      dtdx = dt/dx3[k];

    /* -- K1. initialize rhs with flux difference -- */

      VAR_LOOP(nv) rhs[k][nv] = -dtdx*(flux[k][nv] - flux[k-1][nv]);
      rhs[k][MX3] -= dtdx*(p[k] - p[k-1]);

    /* -- K5. Add dissipative terms to entropy equation -- */

      #if (ENTROPY_SWITCH) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vh[k][RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (visc_flux[k][ENG] - visc_flux[k-1][ENG]);
        rhs_entr -= EXPAND(  vh[j][VX1]*(visc_flux[k][MX1] - visc_flux[k-1][MX1])  ,
                           + vh[j][VX2]*(visc_flux[k][MX2] - visc_flux[k-1][MX2])  ,
                           + vh[j][VX3]*(visc_flux[k][MX3] - visc_flux[k-1][MX3]));
       #endif

       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (tc_flux[k][ENG] - tc_flux[k-1][ENG]);
       #endif
       rhs[k][ENTR] += rhs_entr*dtdx*rhog;              
      #endif
 
    }
  }
}
#elif GEOMETRY == SPHERICAL
{
  double r_1, th, s, s_1;

  if (g_dir == IDIR) { 
    double vc[NVAR], dVdx;

    GetAreaFlux (state, fA, fvA, beg, end, grid);
    for (i = beg; i <= end; i++) { 
      dtdV = dt/dV1[i];
      dtdx = dt/dx1[i];
      r_1  = grid[IDIR].r_1[i];

    /* -- I1. initialize rhs with flux difference -- */

/* Alternative sequence 
dVdx = dV1[i]/dx1[i];
NVAR_LOOP(nv) rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);
rhs[i][MX1] -= dtdx*(p[i] - p[i-1]);
rhs[i][MX3] *= r_1;
#if PHYSICS == MHD
 EXPAND(             
                               ,
      rhs[i][iBTH]  *= dVrdx;  ,
      rhs[i][iBPHI] *= dVrdx;
  )
 #ifdef GLM_MHD
   rhs[i][iBR]     *= dVdx;
 #endif
#endif
IF_DUST(rhs[i][MX3_D] *= r_1;)
*/

      rhs[i][RHO] = -dtdV*(fA[i][RHO] - fA[i-1][RHO]);
      EXPAND(
        rhs[i][iMR]   = - dtdV*(fA[i][iMR] - fA[i-1][iMR])
                        - dtdx*(p[i] - p[i-1]);                    ,
        rhs[i][iMTH]  = -dtdV*(fA[i][iMTH]  - fA[i-1][iMTH]);      ,
        rhs[i][iMPHI] = -dtdV*(fA[i][iMPHI] - fA[i-1][iMPHI])*r_1; 
      )
      #if PHYSICS == MHD
       EXPAND(                                                     
         rhs[i][iBR]   = -dtdV*(fA[i][iBR]   - fA[i-1][iBR]);       ,
         rhs[i][iBTH]  = -dtdx*(fA[i][iBTH]  - fA[i-1][iBTH])*r_1;  ,
         rhs[i][iBPHI] = -dtdx*(fA[i][iBPHI] - fA[i-1][iBPHI])*r_1;
       )
       #ifdef GLM_MHD
        rhs[i][iBR]     = -dtdx*(flux[i][iBR]   - flux[i-1][iBR]);
        rhs[i][PSI_GLM] = -dtdV*(fA[i][PSI_GLM] - fA[i-1][PSI_GLM]);
       #endif
      #endif
      IF_ENERGY(rhs[i][ENG] = -dtdV*(fA[i][ENG] - fA[i-1][ENG]);)

      NSCL_LOOP(nv) rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);
      IF_DUST(NDUST_LOOP(nv) rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);  
              rhs[i][MX3_D] *= r_1;)

    /* -- I5. Add dissipative terms to entropy equation -- */

      #if (ENTROPY_SWITCH) && (PARABOLIC_FLUX & EXPLICIT)
       NVAR_LOOP(nv) vc[nv] = 0.5*(vp[i][nv] + vm[i][nv]);
     
       rhog = vc[RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (fvA[i][ENG] - fvA[i-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(fvA[i][iMR]   - fvA[i-1][iMR])        ,
                           + vc[VX2]*(fvA[i][iMTH]  - fvA[i-1][iMTH])       ,
                           + vc[VX3]*(fvA[i][iMPHI] - fvA[i-1][iMPHI])*r_1);
                    
        rhs_entr -= EXPAND(  vc[VX1]*visc_src[i][MX1], 
                           + vc[VX2]*visc_src[i][MX2],
                           + vc[VX3]*visc_src[i][MX3]);
       #endif
       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (A[i]*tc_flux[i][ENG] - A[i-1]*tc_flux[i-1][ENG]);       
       #endif
       rhs[i][ENTR] += rhs_entr*dtdV*rhog;              
      #endif
    }

  } else if (g_dir == JDIR) {

    GetAreaFlux (state, fA, fvA, beg, end, grid);
    r_1 = 0.5*(x1p[i]*x1p[i] - x1p[i-1]*x1p[i-1])/dV1[i];
    scrh = dt*r_1;
    for (j = beg; j <= end; j++){
      th   = x2[j];
      dtdV = scrh/dV2[j];
      dtdx = scrh/dx2[j];
      s    = sin(th);
      s_1  = 1.0/s;   

    /* -- J1. initialize rhs with flux difference -- */

      rhs[j][RHO] = -dtdV*(fA[j][RHO] - fA[j-1][RHO]);
      EXPAND(
        rhs[j][iMR]   = - dtdV*(fA[j][iMR] - fA[j-1][iMR]);  , 
        rhs[j][iMTH]  = - dtdV*(fA[j][iMTH] - fA[j-1][iMTH])
                        - dtdx*(p[j] - p[j-1]);              , 
        rhs[j][iMPHI] = - dtdV*(fA[j][iMPHI] - fA[j-1][iMPHI])*fabs(s_1);
      )       
      #if PHYSICS == MHD
       EXPAND(                                                     
         rhs[j][iBR]   = -dtdV*(fA[j][iBR]   - fA[j-1][iBR]);   ,
         rhs[j][iBTH]  = -dtdV*(fA[j][iBTH]  - fA[j-1][iBTH]);  ,
         rhs[j][iBPHI] = -dtdx*(flux[j][iBPHI] - flux[j-1][iBPHI]);
       )
       #ifdef GLM_MHD
        rhs[j][iBTH]    = -dtdx*(flux[j][iBTH] - flux[j-1][iBTH]);
        rhs[j][PSI_GLM] = -dtdV*(fA[j][PSI_GLM] - fA[j-1][PSI_GLM]);
       #endif
      #endif
      IF_ENERGY(rhs[j][ENG] = -dtdV*(fA[j][ENG] - fA[j-1][ENG]);)
      
      NSCL_LOOP(nv) rhs[j][nv] = -dtdV*(fA[j][nv] - fA[j-1][nv]);

      IF_DUST(NDUST_LOOP(nv) rhs[j][nv] = -dtdV*(fA[j][nv] - fA[j-1][nv]);
              rhs[j][MX3_D] *= fabs(s_1);)
      
    /* -- J5. Add TC dissipative term to entropy equation -- */

      #if (ENTROPY_SWITCH) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vh[j][RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (fvA[j][ENG] - fvA[j-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(fvA[j][iMR]   - fvA[j-1][iMR])        ,
                           + vc[VX2]*(fvA[j][iMTH]  - fvA[j-1][iMTH])       ,
                           + vc[VX3]*(fvA[j][iMPHI] - fvA[j-1][iMPHI])*fabs(s_1));
                    
        rhs_entr -= EXPAND(  vc[VX1]*visc_src[j][MX1], 
                           + vc[VX2]*visc_src[j][MX2],
                           + vc[VX3]*visc_src[j][MX3]);
       #endif
       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (A[j]*tc_flux[j][ENG] - A[j-1]*tc_flux[j-1][ENG]);
       #endif
       rhs[j][ENTR] += rhs_entr*dtdV*rhog;              
      #endif

    }

  } else if (g_dir == KDIR) {

    r_1  = 0.5*(x1p[i]*x1p[i] - x1p[i-1]*x1p[i-1])/dV1[i];
    scrh = dt*r_1*dx2[j]/dV2[j];
    for (k = beg; k <= end; k++) {
      dtdx = scrh/dx3[k];

    /* -- K1.  initialize rhs with flux difference -- */

      VAR_LOOP(nv) rhs[k][nv] = -dtdx*(flux[k][nv] - flux[k-1][nv]);
      rhs[k][iMPHI] -= dtdx*(p[k] - p[k-1]); 

    /* -- K5. Add dissipative terms to entropy equation -- */

      #if (ENTROPY_SWITCH) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vh[k][RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (visc_flux[k][ENG] - visc_flux[k-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(visc_flux[k][MX1] - visc_flux[k-1][MX1])  ,
                           + vc[VX2]*(visc_flux[k][MX2] - visc_flux[k-1][MX2])  ,
                           + vc[VX3]*(visc_flux[k][MX3] - visc_flux[k-1][MX3]));
       #endif
       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (tc_flux[k][ENG] - tc_flux[k-1][ENG]);
       #endif
       rhs[k][ENTR] += rhs_entr*dtdx*rhog;              
      #endif      
    }
  }
}
#endif  /* GEOMETRY == SPHERICAL */

/* --------------------------------------------------------------
    Add source terms
   -------------------------------------------------------------- */

  RightHandSideSource (state, Dts, beg, end, dt, phi_p, grid);

/* --------------------------------------------------
    Reset right hand side in internal boundary zones
   -------------------------------------------------- */
   
  #if INTERNAL_BOUNDARY == YES
   InternalBoundaryReset(state, Dts, beg, end, grid);
  #endif
  
/* --------------------------------------------------
           Time step determination
   -------------------------------------------------- */

#if !GET_MAX_DT
return;
#endif

  cl = 0.0;
  for (i = beg-1; i <= end; i++) {
    scrh = Dts->cmax[i]*grid[g_dir].inv_dxi[i];
    cl = MAX(cl, scrh);   
  }
  #if GEOMETRY == POLAR || GEOMETRY == SPHERICAL
   if (g_dir == JDIR) {
     cl /= fabs(grid[IDIR].xgc[g_i]);
   }
   #if GEOMETRY == SPHERICAL
    if (g_dir == KDIR){
      cl /= fabs(grid[IDIR].xgc[g_i])*sin(grid[JDIR].xgc[g_j]);
    }
   #endif
  #endif
  Dts->inv_dta = MAX(cl, Dts->inv_dta);  
}

/* ********************************************************************* */
void TotalFlux (const State_1D *state, double *phi_p,
                int beg, int end, Grid *grid)
/*!
 * Compute the total flux in order to enforce conservation of 
 * angular momentum and energy in presence of FARGO source 
 * terms, rotation or gravitational potential.
 *
 * In Cartesian coordinates, for instance, we modify the energy and
 * the y-momentum flux as:
 * \f{eqnarray*}{
 *    \F^{E}_{i+\HALF} &\quad \to\quad&
 *    \F^{E}_{i+\HALF} + w_{i+\HALF}\left(\HALF w_{i+\HALF}\F^{\rho}_{i+\HALF}
 *     + \F^{m_y}_{i+\HALF}\right) +  \Phi_{i+\HALF}\F^{\rho}_{i+\HALF}                           
 *    \\ \noalign{\medskip}
 *    \F^{m_y}_{i+\HALF} &\quad \to\quad&
 *    \F^{m_y}_{i+\HALF} + w_{i+\HALF}\F^{\rho}_{i+\HALF}
 * \f}
 * where \f$w\f$ is the average orbital velocity computed during the
 * FARGO algorithm and \f$\Phi\f$ is the gravitational potential.
 * This is described in Mignone et al. (A\&A 2012), see appendix there.
 *
 * \note This function does not change the fluxes when the Shearing-box
 *       and FARGO modules are both enabled.
 *       Instead, we incorporate the source terms elsewhere.
 * 
 *  
 * \param [in]     state  pointer to State_1D structure;
 * \param [in,out] phi_p  1D array defining the gravitational potential
 *                        at grid interfaces;
 * \param [in]     beg    initial index of computation; 
 * \param [in]     end    final   index of computation;
 * \param [in]     grid   pointer to Grid structure;
 *
 * \b Reference
 *    - "A conservative orbital advection scheme for simulations of magnetized
 *        shear flows with the PLUTO code"
 *        Mignone et al., A\&A (2012) 545A, 152M
 *
 *********************************************************************** */
#ifndef iMPHI
 #define iMPHI MX2  /* -- for Cartesian coordinates -- */
#endif
{
  int i,j,k;
  double wp, R;
  double **flux, *vp;
  double *x1,  *x2,  *x3;
  double *x1p, *x2p, *x3p;
#ifdef FARGO
  double **wA;
  wA = FARGO_GetVelocity();
#endif

  flux = state->flux;
  x1  = grid[IDIR].x;  x1p = grid[IDIR].xr;
  x2  = grid[JDIR].x;  x2p = grid[JDIR].xr;
  x3  = grid[KDIR].x;  x3p = grid[KDIR].xr;

  i = g_i;  /* will be redefined during x1-sweep */
  j = g_j;  /* will be redefined during x2-sweep */
  k = g_k;  /* will be redefined during x3-sweep */

  if (g_dir == IDIR){ 
    for (i = beg; i <= end; i++){

    /* ----------------------------------------------------
        include flux contributions from FARGO or Rotation 
        Note: ShearingBox terms are not included here but
              only within the BodyForce function.
       ---------------------------------------------------- */

    #if (defined FARGO && !defined SHEARINGBOX) || (ROTATING_FRAME == YES)
      wp = 0.0;
      #if GEOMETRY == SPHERICAL
      IF_FARGO(wp = 0.5*(wA[j][i] + wA[j][i+1]);)
      R = x1p[i]*sin(x2[j]);  /* -- cylindrical radius -- */
      #else
      IF_FARGO(wp = 0.5*(wA[k][i] + wA[k][i+1]);)
      R = x1p[i];                   /* -- cylindrical radius -- */
      #endif
      IF_ROTATING_FRAME(wp += g_OmegaZ*R;)
      IF_ENERGY  (flux[i][ENG] += wp*(0.5*wp*flux[i][RHO] + flux[i][iMPHI]);)
      flux[i][iMPHI] += wp*flux[i][RHO];
    #endif

    /* -- gravitational potential -- */

    #if (BODY_FORCE & POTENTIAL) && HAVE_ENERGY
      flux[i][ENG] += flux[i][RHO]*phi_p[i];                          
    #endif
    }
  }else if (g_dir == JDIR){
    for (j = beg; j <= end; j++){ 

    /* ----------------------------------------------------
        include flux contributions from FARGO and Rotation 
       ---------------------------------------------------- */

    #if GEOMETRY == SPHERICAL
      #if (defined FARGO) || (ROTATING_FRAME == YES)
      wp = 0.0;
      R  = x1[i]*sin(x2p[j]);
      IF_FARGO   (wp += 0.5*(wA[j][i] + wA[j+1][i]);)
      IF_ROTATING_FRAME(wp += g_OmegaZ*R;)
      IF_ENERGY  (flux[j][ENG] += wp*(0.5*wp*flux[j][RHO] + flux[j][iMPHI]);)
      flux[j][iMPHI] += wp*flux[j][RHO];
      #endif
    #endif

    #if (BODY_FORCE & POTENTIAL) && HAVE_ENERGY
      flux[j][ENG] += flux[j][RHO]*phi_p[j];
    #endif      

    }
  }else if (g_dir == KDIR){
    for (k = beg; k <= end; k++) {

    /* ----------------------------------------------------
        include flux contributions from FARGO
        (polar/cartesian geometries only)
       ---------------------------------------------------- */

    #if (GEOMETRY != SPHERICAL) && (defined FARGO && !defined SHEARINGBOX)
      wp = 0.5*(wA[k][i] + wA[k+1][i]);
      IF_ENERGY(flux[k][ENG] += wp*(0.5*wp*flux[k][RHO] + flux[k][iMPHI]);)
      flux[k][iMPHI] += wp*flux[k][RHO];
    #endif

    #if (BODY_FORCE & POTENTIAL) && HAVE_ENERGY
      flux[k][ENG] += flux[k][RHO]*phi_p[k];
    #endif
    }
  }
}
