/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute the thermal conduction flux.

  Compute the thermal conduction flux along one row of computational 
  zones for the HD and MHD  modules according to Spitzer (1962):
  \f[
     \vec{F}_c = \frac{q}{|\vec{F}_{\rm class}| + q}\vec{F}_{\rm class}
  \f] 
  where \f$ \vec{F}_{\rm class} = \kappa_\parallel \nabla T\f$ is the
  classical (hydrodynamic) thermal conduction flux,
  \f$ q = F_{\rm sat} = 5\phi \rho c_{\rm iso}^3\f$ is the saturated flux.
  Temperature is, at present, simply computed as \f$p/\rho\f$.
  
  The classical flux is purely parabolic, and it is discretized using
  standard finite differences. 
  The saturated flux is included only when the macro ::TC_SATURATED_FLUX
  is enabled and it is treated in an upwind manner following the guidelines
  given in the appendix of Mignone et al. 2012 (see also Balsara (2008)
  for alternative discretization methods).

  In MHD, the classical flux further splits into 2 components, along
  and across the magnetic field lines (see Eq. [7] of Mignone et al.).

  This function also computes the inverse of the time
  step and return its maximum over the current sweep.
  
  \b References
     - "The PLUTO Code for Adaptive Mesh Computations in Astrophysical 
        Fluid Dynamics" \n
       Mignone et al, ApJS (2012) 198, 7M
     - "Simulating anisotropic thermal conduction in supernova remnants
        - I. Numerical methods",
        Balsara, Tilley, Howk, MNRAS (2008) 386, 627 
       
  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos
  \date    June 20, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#ifndef TC_SATURATED_FLUX
/*! When set to YES (default), include saturation effects when
    computing the  conduction flux.                                */ 
  #define TC_SATURATED_FLUX  YES
#endif

/* ********************************************************************* */
void TC_Flux (double ***T, const Sweep *sweep,  
              double *dcoeff, int beg, int end, Grid *grid)
/*! 
 * Compute the thermal conduction flux, sweep->par_flx.
 *
 * \param [in]     T       3D array containing the dimensionless 
 *                         temperature
 * \param [in,out] sweep   pointer to a Sweep structure
 * \param [out]    dcoeff  the diffusion coefficient needed for computing
 *                         the time step.
 * \param [in]     beg     initial index of computation
 * \param [in]     end     final   index of computation
 * \param [in]     grid    pointer to an array of Grid structures
 *
 * \return This function has no return value.                       
 *********************************************************************** */
{
  int  i, j, k, nv;
  double bgradT, Bmag, Fc_mag, dT_mag;
  double Fc, Fsat, g1 = g_gamma - 1.0;
  double x1, x2, x3;
  double dvp, dvm, dvl;
  double alpha, uL, uR, bn;
  double vi[NVAR], kpar=0.0, knor=0.0, phi;
  double bck_fld[3];
  double **vc = sweep->vn;
  static double **gradT, *pp, *pm;

/* -----------------------------------------------------------
   1. Allocate memory, compute temperature gradient in the
      required direction and set 2 out of the 3 coordinates.
   ----------------------------------------------------------- */

  if (gradT == NULL) {
    gradT = ARRAY_2D(NMAX_POINT, 3, double);
    pp    = ARRAY_1D(NMAX_POINT, double);
    pm    = ARRAY_1D(NMAX_POINT, double);
  }

  GetGradient (T, gradT, beg, end, grid);
  if (g_dir == JDIR || g_dir == KDIR) x1 = grid->x[IDIR][g_i];
  if (g_dir == IDIR || g_dir == KDIR) x2 = grid->x[JDIR][g_j];
  if (g_dir == IDIR || g_dir == JDIR) x3 = grid->x[KDIR][g_k];

/* --------------------------------------------------------------
   2. Compute limited slopes in pressure (only for saturated TC)
   -------------------------------------------------------------- */

#if TC_SATURATED_FLUX == YES  
  #if (defined CTU) && (THERMAL_CONDUCTION == EXPLICIT)
  if (g_intStage == 1){             /* Predictor stage of CTU requires */
    for (i = beg; i <= end+1; i++){ /* 1st order states                */
      pp[i] = vc[i][PRS];
      pm[i] = vc[i][PRS];
    }
  }else
  #endif
  for (i = beg; i <= end+1; i++){
    dvp = vc[i+1][PRS] - vc[i][PRS];
    dvm = vc[i][PRS] - vc[i-1][PRS];
    dvl = VAN_LEER(dvp, dvm);
    pp[i] = vc[i][PRS] + 0.5*dvl;
    pm[i] = vc[i][PRS] - 0.5*dvl;
  }
#endif

/* ----------------------------------------------- 
   3. Compute Thermal Conduction Flux (tcflx).
   ----------------------------------------------- */

  for (i = beg; i <= end; i++){
    
  /* -- 3a. Compute interface values -- */

    NVAR_LOOP(nv) vi[nv]  = 0.5*(vc[i][nv] + vc[i+1][nv]);    

  /* -- 3b. Get TC coefficients at cell interfaces -- */

    if (g_dir == IDIR) x1 = grid->xr[IDIR][i];
    if (g_dir == JDIR) x2 = grid->xr[JDIR][i];
    if (g_dir == KDIR) x3 = grid->xr[KDIR][i];

#if PHYSICS == MHD && BACKGROUND_FIELD == YES
    BackgroundField(x1,x2,x3, bck_fld);
    EXPAND(vi[BX1] += bck_fld[0];  ,
           vi[BX2] += bck_fld[1];  ,
           vi[BX3] += bck_fld[2];)
#endif

    TC_kappa(vi, x1, x2, x3, &kpar, &knor, &phi);
    dT_mag  = D_EXPAND(  gradT[i][0]*gradT[i][0], 
                       + gradT[i][1]*gradT[i][1], 
                       + gradT[i][2]*gradT[i][2]);
    dT_mag = sqrt(dT_mag) + 1.e-12;


  /* -- 3c. Compute the classical TC flux and diffusion coeff -- */

#if PHYSICS == HD

    Fc        = kpar*gradT[i][g_dir];  /* Classical thermal conduction flux */
    Fc_mag    = kpar*dT_mag;          /* Flux magnitude (include all dirs) */
    dcoeff[i] = fabs(kpar/vi[RHO])*g1; /* Diffusion coefficient */

#elif PHYSICS == MHD

    Bmag = EXPAND(vi[BX1]*vi[BX1], + vi[BX2]*vi[BX2], + vi[BX3]*vi[BX3]);
    Bmag = sqrt(Bmag) + 1.e-12;

    bgradT = D_EXPAND(  vi[BX1]*gradT[i][IDIR], 
                      + vi[BX2]*gradT[i][JDIR], 
                      + vi[BX3]*gradT[i][KDIR]);
    bgradT /= Bmag;

    bn     = vi[BX1 + g_dir]/Bmag; /* -- unit vector component -- */
    Fc     = kpar*bgradT*bn + knor*(gradT[i][g_dir] - bn*bgradT);

    Fc_mag = sqrt(  (kpar*kpar - knor*knor)*bgradT*bgradT
                   + knor*knor*dT_mag*dT_mag); 

    dcoeff[i] = fabs(Fc_mag/dT_mag*bn/vi[RHO])*g1;

#endif

  /* -- 3d. Compute saturated flux using upwinding -- */

#if TC_SATURATED_FLUX == YES
    uL       = pp[i];
    uR       = pm[i+1];
    if      (Fc > 0.0) Fsat = 5.0*phi/sqrt(vi[RHO])*uR*sqrt(uR);
    else if (Fc < 0.0) Fsat = 5.0*phi/sqrt(vi[RHO])*uL*sqrt(uL);
    else               Fsat = 5.0*phi*vi[PRS]*sqrt(vi[PRS]/vi[RHO]);

    alpha      = Fsat/(Fsat + Fc_mag);
    Fc        *= alpha;
    dcoeff[i] *= alpha;
#endif

    sweep->tc_flux[i][ENG] = Fc;
  }
}
