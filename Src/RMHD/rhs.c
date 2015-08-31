#include "pluto.h"

#ifdef CH_SPACEDIM   /*  implies Chombo is being used  */
 #define USE_PR_GRADIENT  YES   
#else
 #ifdef FINITE_DIFFERENCE 
  #define USE_PR_GRADIENT  NO   /* -- default for Finite Difference schemes -- */
 #else
  #define USE_PR_GRADIENT  YES   /* -- default, do not change!! -- */
 #endif
#endif

/* *********************************************************** */
void RightHandSide (const State_1D *state, Time_Step *Dts, 
                    int beg, int end, double dt, Grid *grid)
/* 
 *   Compute right hand side of the MHD equations in different geometries,
 *   taking contributions from one direction at a time. 
 *   The right hand side (rhs) consists of the following contributions:
 * 
 *    rhs = dt/dV * (Ap*Fp - Am*Fm) + dt*S
 * 
 *   where 
 *
 *    Ap, Am:  interface areas, 
 *    Fp, Fm:  interface fluxes,
 *    dt:      time step
 *    S:       source term including 
 *             * geometrical source terms
 *             * gravity
 *  
 *   In order to compute rhs, this function takes the following steps:
 *
 *    - initialize rhs with flux differences  (#1)
 *    - add geometrical source terms          (#2)
 *    - add gravity                           (#4)
 *
 * LAST MODIFIED
 * 
 *   July, 31  2012 by A. Mignone (mignone@ph.unito.it)
 *
 ************************************************************* */
{
  int    i, j, k, nv;
  double dtdx, dtdV, scrh;
  double *x1, *x1p, *dx1, *dV1;
  double *x2, *x2p, *dx2, *dV2;
  double *x3, *x3p, *dx3, *dV3;
  double **flux, **rhs, *p, *v;
  double cl;
  double g[3];
  static double **fA, *h;
  
  #if GEOMETRY != CARTESIAN
   if (fA == NULL) {
     fA = ARRAY_2D(NMAX_POINT, NVAR, double);
     h  = ARRAY_1D(NMAX_POINT, double);
   }
  #endif

/* --------------------------------------------------
             Compute passive scalar fluxes
   -------------------------------------------------- */

  #if NSCL > 0
   AdvectFlux (state, beg - 1, end, grid);
  #endif

/* --------------------------
      pointer shortcuts
   -------------------------- */

  rhs  = state->rhs;
  flux = state->flux;
  p    = state->press;
  
  x1  = grid[IDIR].x;  x2  = grid[JDIR].x;  x3  = grid[KDIR].x;
  x1p = grid[IDIR].xr; x2p = grid[JDIR].xr; x3p = grid[KDIR].xr;
  dx1 = grid[IDIR].dx; dx2 = grid[JDIR].dx; dx3 = grid[KDIR].dx;
  dV1 = grid[IDIR].dV; dV2 = grid[JDIR].dV; dV3 = grid[KDIR].dV;

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

#if GEOMETRY == CARTESIAN
/* ***********************************************************
   
        Compute right-hand side of the RMHD/RHD 
        equations in CARTESIAN geometry.
    
   *********************************************************** */
{
  double x, y, z;

  if (g_dir == IDIR){

  /* ****************************************************
      Cartesian x-direction,

       - initialize rhs with flux differences (I1)
       - enforce conservation of total angular
         momentum and/or energy               (I3)
       - add gravity                          (I4)
     **************************************************** */

    y = x2[j];
    z = x3[k];
    for (i = beg; i <= end; i++) {
      x    = x1[i];
      dtdx = dt/dx1[i];

    /* -----------------------------------------------
       I1. initialize rhs with flux difference
       ----------------------------------------------- */

      NVAR_LOOP(nv) rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i-1][nv]);
      #if USE_PR_GRADIENT == YES
       rhs[i][MX1] -= dtdx*(p[i] - p[i-1]);
      #endif

    /* ----------------------------------------------------
       I4. Include gravity
       ---------------------------------------------------- */

      v = state->vh[i];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(v, g, x1[i], x2[j], x3[k]);
       rhs[i][MX1] += dt*v[RHO]*g[IDIR];
       #if DIMENSIONS == 1 && COMPONENTS > 1 /* In 1D and when COMPONENTS == 2  */
                                             /* or 3, add gravity contributions */
                                             /* along non-existing dimensions. */
        EXPAND(                                     ,
               rhs[i][MX2] += dt*vc[RHO]*g[JDIR];   ,
               rhs[i][MX3] += dt*vc[RHO]*g[KDIR];)
       #endif

       #if HAVE_ENERGY
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[IDIR];
        #if DIMENSIONS == 1 && COMPONENTS > 1
         rhs[i][ENG] += dt*(EXPAND(0.0, + vc[RHO]*vc[VX2]*g[JDIR],
                                        + vc[RHO]*vc[VX3]*g[KDIR]));
        #endif
       #endif
      #endif
    }
  } else if (g_dir == JDIR){

  /* ****************************************************
      Cartesian y-direction,

       - initialize rhs with flux differences (J1)
       - add gravity                          (J4)
     **************************************************** */

    x = x1[i];
    z = x3[k];
    for (j = beg; j <= end; j++) {
      y    = x2[j];
      dtdx = dt/dx2[j];

    /* -----------------------------------------------
       J1. initialize rhs with flux difference
       ----------------------------------------------- */

      NVAR_LOOP(nv) rhs[j][nv] = -dtdx*(flux[j][nv] - flux[j-1][nv]);
      #if USE_PR_GRADIENT == YES
       rhs[j][MX2] -= dtdx*(p[j] - p[j-1]);
      #endif

    /* ----------------------------------------------------
       J4. Include gravity
       ---------------------------------------------------- */

      v = state->vh[j];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(v, g, x1[i], x2[j], x3[k]);
       rhs[j][MX2] += dt*v[RHO]*g[JDIR];
       #if DIMENSIONS == 2 && COMPONENTS == 3
        rhs[j][MX3] += dt*vc[RHO]*g[KDIR];
       #endif
       #if HAVE_ENERGY
        rhs[j][ENG] += dt*0.5*(flux[j][RHO] + flux[j-1][RHO])*g[JDIR];
        #if DIMENSIONS == 2 && COMPONENTS == 3
         rhs[j][ENG] += dt*vc[RHO]*vc[VX3]*g[KDIR];
        #endif
       #endif
      #endif
    }

  }else if (g_dir == KDIR){

  /* ****************************************************
      Cartesian z-direction,

       - initialize rhs with flux differences (K1)
       - enforce conservation of total angular
         momentum and/or energy               (K3)
       - add gravity                          (K4)
     **************************************************** */

    x = x1[i];
    y = x2[j];
    for (k = beg; k <= end; k++) {
      z    = x3[k];
      dtdx = dt/dx3[k];

    /* -----------------------------------------------
       K1. initialize rhs with flux difference
       ----------------------------------------------- */

      NVAR_LOOP(nv) rhs[k][nv] = -dtdx*(flux[k][nv] - flux[k-1][nv]);
      #if USE_PR_GRADIENT == YES
       rhs[k][MX3] -= dtdx*(p[k] - p[k-1]);
      #endif

    /* ----------------------------------------------------
       K4. Include gravity
       ---------------------------------------------------- */

      v = state->vh[k];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(v, g, x1[i], x2[j], x3[k]);
       rhs[k][MX3] += dt*v[RHO]*g[KDIR];
       #if HAVE_ENERGY
        rhs[k][ENG] += dt*0.5*(flux[k][RHO] + flux[k-1][RHO])*g[KDIR];
       #endif
      #endif
    }
  }
}
#elif GEOMETRY == CYLINDRICAL

/* ***********************************************************
   
        Compute right-hand side of the RMHD/RHD 
        equations in CYLINDRICAL geometry.
    
   *********************************************************** */
{
  double R, z, phi; 
  double R1, R_1;
  double vB, vel2, lor2, wt, mphi, B2;

  if (g_dir == IDIR) {  

  /* ****************************************************
      Cylindrical radial direction:
      multiply fluxes times interface area
     **************************************************** */

    z   = x2[g_j];
    phi = 0.0;
    for (i = beg - 1; i <= end; i++){ 
      R = grid[IDIR].A[i];

      fA[i][RHO] = flux[i][RHO]*R;
      EXPAND(fA[i][iMR]   = flux[i][iMR]*R;     ,
             fA[i][iMZ]   = flux[i][iMZ]*R;     ,
             fA[i][iMPHI] = flux[i][iMPHI]*R*R;)       
      #if PHYSICS == RMHD
       EXPAND(fA[i][iBR]   = flux[i][iBR]*R;  ,
              fA[i][iBZ]   = flux[i][iBZ]*R;  ,
              fA[i][iBPHI] = flux[i][iBPHI]*R;)
      #endif
      #if HAVE_ENERGY
       fA[i][ENG] = flux[i][ENG]*R;
      #endif
      #ifdef GLM_MHD
       fA[i][PSI_GLM] = flux[i][PSI_GLM]*R;
      #endif
       NSCL_LOOP(nv) fA[i][nv] = flux[i][nv]*R;
    }

  /* ****************************************************
      Cylindrical radial direction,

       - initialize rhs with flux differences (I1)
       - add source terms                     (I2)
       - add gravity                          (I4)
     **************************************************** */

    #if COMPONENTS == 3
     Enthalpy (state->v, h, beg, end);
    #endif
    for (i = beg; i <= end; i++){ 
      R    = x1[i];
      dtdV = dt/dV1[i];
      dtdx = dt/dx1[i];
      R_1  = 1.0/R;

    /* -----------------------------------------------
       I1. initialize rhs with flux difference
       ----------------------------------------------- */

      rhs[i][RHO] = -dtdV*(fA[i][RHO] - fA[i-1][RHO]);
      EXPAND(rhs[i][iMR]   = - dtdV*(fA[i][iMR]   - fA[i-1][iMR]);  ,
             rhs[i][iMZ]   = - dtdV*(fA[i][iMZ]   - fA[i-1][iMZ]);  ,
             rhs[i][iMPHI] = - dtdV*(fA[i][iMPHI] - fA[i-1][iMPHI])*fabs(R_1);)
      #if USE_PR_GRADIENT == YES
       rhs[i][iMR] -= dtdx*(p[i] - p[i-1]);  
      #endif
      #if PHYSICS == RMHD
       EXPAND(rhs[i][iBR]   = - dtdV*(fA[i][iBR]   - fA[i-1][iBR]);  ,
              rhs[i][iBZ]   = - dtdV*(fA[i][iBZ]   - fA[i-1][iBZ]);  ,
              rhs[i][iBPHI] = - dtdV*(fA[i][iBPHI] - fA[i-1][iBPHI]);)
       #ifdef GLM_MHD
        rhs[i][iBR]     = - dtdx*(flux[i][iBR]   - flux[i-1][iBR]);
        rhs[i][PSI_GLM] = - dtdV*(fA[i][PSI_GLM] - fA[i-1][PSI_GLM]);
       #endif
      #endif
      #if HAVE_ENERGY
       rhs[i][ENG] = -dtdV*(fA[i][ENG] - fA[i-1][ENG]);
      #endif
       NSCL_LOOP(nv) rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);
      

    /* -------------------------------------------------------
       I2. Add source terms
           [for Bphi we use the non-conservative formulation
            with the source term since it has been found to 
            be more stable at low resolution in toroidal jet
            simulations]
       ------------------------------------------------------- */

      v     = state->vh[i];
      #if COMPONENTS == 3
       vel2  = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
       lor2  = 1.0/(1.0 - vel2);
       #if PHYSICS == RMHD
        vB   = EXPAND(v[VX1]*v[BX1], + v[VX2]*v[BX2], + v[VX3]*v[BX3]);
        B2   = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
        wt   = v[RHO]*h[i]*lor2 + B2;
        mphi = wt*v[iVPHI] - vB*v[iBPHI]; 
       #elif PHYSICS == RHD
        wt   = v[RHO]*h[i]*lor2;
        mphi = wt*v[iVPHI]; 
       #endif       
      
       rhs[i][iMR] += dt*mphi*v[iVPHI]*R_1;
       #if PHYSICS == RMHD
        rhs[i][iMR]   -= dt*(v[iBPHI]/lor2 + vB*v[iVPHI])*v[iBPHI]*R_1;
        rhs[i][iBPHI] -= dt*(v[iVPHI]*v[iBR] - v[iBPHI]*v[iVR])*R_1;   
       #endif
      #endif  /* COMPONENTS == 3 */

    /* ----------------------------------------------------
       I4. Include gravity
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
       BodyForceVector(v, g, x1[i], x2[g_j], x3[g_k]);
       rhs[i][iMR] += dt*v[RHO]*g[g_dir];
       #if HAVE_ENERGY
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[g_dir];
       #endif
      #endif
    }
     
  } else if (g_dir == JDIR) { 

  /* ****************************************************
      Cylindrical vertical direction:

       - initialize rhs with flux differences (J1)
       - add gravity                          (J4)
     **************************************************** */

    R   = x1[i];
    phi = 0.0;
    for (j = beg; j <= end; j++){ 
      z    = x2[j];   
      dtdx = dt/dx2[j];

    /* -----------------------------------------------
       J1. initialize rhs with flux difference
       ----------------------------------------------- */

      NVAR_LOOP(nv) rhs[j][nv] = -dtdx*(flux[j][nv] - flux[j-1][nv]);
      #if USE_PR_GRADIENT == YES
       rhs[j][iMZ] += - dtdx*(p[j] - p[j-1]);
      #endif

    /* ----------------------------------------------------
       J4. Include gravity
       ---------------------------------------------------- */

      v = state->vh[j];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(v, g, x1[i], x2[j], x3[k]);
       rhs[j][iMZ] += dt*v[RHO]*g[JDIR];
       #if HAVE_ENERGY
        rhs[j][ENG] += dt*0.5*(flux[j][RHO] + flux[j-1][RHO])*g[JDIR];
       #endif
      #endif
    }
  }
}

#elif GEOMETRY == POLAR

/* ***********************************************************
   
        Compute right-hand side of the RMHD/RHD 
        equations in POLAR geometry.
    
   *********************************************************** */
{
  double R, phi, z; 
  double R_1;
  double vB, vel2, g_2, wt, mphi, B2;
   
  if (g_dir == IDIR) { 

  /* ****************************************************
      Polar radial direction:
      multiply fluxes times interface area
     **************************************************** */

    phi = x2[j];
    z   = x3[k];
    for (i = beg - 1; i <= end; i++) { 
      R = grid[IDIR].A[i];

      fA[i][RHO] = flux[i][RHO]*R;
      EXPAND(fA[i][iMR]   = flux[i][iMR]*R;      ,
             fA[i][iMPHI] = flux[i][iMPHI]*R*R;  ,
             fA[i][iMZ]   = flux[i][iMZ]*R;)       
      #if PHYSICS == RMHD
       EXPAND(fA[i][iBR]   = flux[i][iBR]*R;    ,
              fA[i][iBPHI] = flux[i][iBPHI];    ,
              fA[i][iBZ]   = flux[i][iBZ]*R;)
      #endif
      #if HAVE_ENERGY
       fA[i][ENG] = flux[i][ENG]*R;
      #endif
      #ifdef GLM_MHD
       fA[i][PSI_GLM] = flux[i][PSI_GLM]*R;
      #endif
       NSCL_LOOP(nv) fA[i][nv] = flux[i][nv]*R;
    }

  /* ****************************************************
      Polar radial direction,

       - initialize rhs with flux differences (I1)
       - add source terms                     (I2)
       - enforce conservation of total angular
         momentum and/or energy               (I3)
       - add gravity                          (I4)
     **************************************************** */

    #if COMPONENTS >= 2  /* -- need enthalpy for source term computation -- */
     Enthalpy (state->v, h, beg, end);
    #endif
    for (i = beg; i <= end; i++) {
      R    = x1[i];
      dtdV = dt/dV1[i];
      dtdx = dt/dx1[i];
      R_1  = grid[IDIR].r_1[i];

    /* -----------------------------------------------
       I1. initialize rhs with flux difference
       ----------------------------------------------- */

      rhs[i][RHO] = -dtdV*(fA[i][RHO] - fA[i-1][RHO]);
      EXPAND(rhs[i][iMR]   = - dtdV*(fA[i][iMR]   - fA[i-1][iMR])
                             - dtdx*(p[i] - p[i-1]);                      ,      
             rhs[i][iMPHI] = - dtdV*(fA[i][iMPHI] - fA[i-1][iMPHI])*R_1;  ,
             rhs[i][iMZ]   = - dtdV*(fA[i][iMZ]   - fA[i-1][iMZ]);)
      #if PHYSICS == RMHD
       EXPAND(rhs[i][iBR]   = -dtdV*(fA[i][iBR]   - fA[i-1][iBR]);    ,
              rhs[i][iBPHI] = -dtdx*(fA[i][iBPHI] - fA[i-1][iBPHI]);  ,
              rhs[i][iBZ]   = -dtdV*(fA[i][iBZ]   - fA[i-1][iBZ]);)
       #ifdef GLM_MHD
        rhs[i][iBR]     = -dtdx*(flux[i][iBR]   - flux[i-1][iBR]);
        rhs[i][PSI_GLM] = -dtdV*(fA[i][PSI_GLM] - fA[i-1][PSI_GLM]);
       #endif
      #endif
      #if HAVE_ENERGY
       rhs[i][ENG] = -dtdV*(fA[i][ENG] - fA[i-1][ENG]);
      #endif
      
       NSCL_LOOP(nv)  rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);
      

    /* ----------------------------------------------------
       I2. Add source terms to the radial momentum eqn.
           S = w*gamma^2*v(phi)^2 - B(phi)^2 - E(phi)^2
       ---------------------------------------------------- */

      v     = state->vh[i];
      #if COMPONENTS >= 2
       vel2  = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
       g_2   = 1.0 - vel2;
       #if PHYSICS == RMHD
        vB   = EXPAND(v[VX1]*v[BX1], + v[VX2]*v[BX2], + v[VX3]*v[BX3]);
        B2   = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
        wt   = v[RHO]*h[i]/g_2 + B2;
        mphi = wt*v[iVPHI] - vB*v[iBPHI]; 
       #elif PHYSICS == RHD
        wt   = v[RHO]*h[i]/g_2;
        mphi = wt*v[iVPHI]; 
       #endif       
      
       rhs[i][iMR] += dt*mphi*v[iVPHI]*R_1;
       #if PHYSICS == RMHD
        rhs[i][iMR] -= dt*(v[iBPHI]*g_2 + vB*v[iVPHI])*v[iBPHI]*R_1;
       #endif
      #endif  /* COMPONENTS >= 2 */
 
    /* ----------------------------------------------------
       I4. Include gravity
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
       BodyForceVector(v, g, x1[i], x2[j], x3[k]);
       rhs[i][iMR] += dt*v[RHO]*g[IDIR];
       #if HAVE_ENERGY
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[JDIR];
       #endif
      #endif
    }
     
  } else if (g_dir == JDIR) {

  /* ****************************************************
      Polar azimuthal direction:

       - initialize rhs with flux differences (J1)
       - add gravity                          (J4)
     **************************************************** */

    R = x1[i];
    z = x3[k];
    scrh = dt/R;
    for (j = beg; j <= end; j++){ 
      phi  = x2[j];
      dtdx = scrh/dx2[j];

    /* ------------------------------------------------
       J1. Compute equations rhs for phi-contributions
       ------------------------------------------------ */

      NVAR_LOOP(nv) rhs[j][nv] = -dtdx*(flux[j][nv] - flux[j-1][nv]);
      rhs[j][iMPHI] -= dtdx*(p[j] - p[j-1]);

    /* -------------------------------------------------------
       J4. Include gravity
       ------------------------------------------------------- */

      v = state->vh[j];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(v, g, x1[i], x2[j], x3[k]);
       rhs[j][iMPHI] += dt*v[RHO]*g[JDIR];
       #if HAVE_ENERGY
        rhs[j][ENG] += dt*0.5*(flux[j][RHO] + flux[j-1][RHO])*g[JDIR];
       #endif
      #endif
    }

  } else if (g_dir == KDIR) { 

  /* ****************************************************
      Polar vertical direction:

       - initialize rhs with flux differences (K1)
       - enforce conservation of total angular
         momentum and/or energy               (K3)
       - add gravity                          (K4)
     **************************************************** */

    R   = x1[i];
    phi = x2[j];
    for (k = beg; k <= end; k++){ 
      z    = x3[k];
      dtdx = dt/dx3[k];

    /* -----------------------------------------------
       K1. initialize rhs with flux difference
       ----------------------------------------------- */

      NVAR_LOOP(nv) rhs[k][nv] = -dtdx*(flux[k][nv] - flux[k-1][nv]);
      rhs[k][iMZ] -= dtdx*(p[k] - p[k-1]);

    /* ----------------------------------------------------
       K4. Include gravity
       ---------------------------------------------------- */

      v = state->vh[k];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(v, g, x1[i], x2[j], x3[k]); 
       rhs[k][iMZ] += dt*v[RHO]*g[KDIR];
       #if HAVE_ENERGY
        rhs[k][ENG] += dt*0.5*(flux[k][RHO] + flux[k-1][RHO])*g[KDIR];
       #endif
      #endif
    }
  }
}
#elif GEOMETRY == SPHERICAL

/* ***********************************************************
   
        Compute right-hand side of the RMHD/RHD 
        equations in SPHERICAL geometry.
    
   *********************************************************** */
{
  double r, th, phi;
  double r2, r3, r_1;
  double s, s2, ct, s_1;
  double vB, vel2, lor2, wt, B2;
  double mth = 0.0, mphi = 0.0;

  if (g_dir == IDIR) { 
    double Sm;

  /* ****************************************************
      Spherical radial direction: 
      multiply fluxes by interface area 
     **************************************************** */

    th  = x2[j]; s = sin(th);
    phi = x3[k];
    for (i = beg - 1; i <= end; i++){
      r  = x1p[i];
      r2 = r*r; 
      r3 = r2*r;

      fA[i][RHO] = flux[i][RHO]*r2;
      EXPAND(fA[i][iMR]   = flux[i][iMR]*r2;   ,
             fA[i][iMTH]  = flux[i][iMTH]*r2;  ,
             fA[i][iMPHI] = flux[i][iMPHI]*r3;)
      #if PHYSICS == RMHD
       EXPAND(fA[i][iBR]   = flux[i][iBR]*r2;   ,
              fA[i][iBTH]  = flux[i][iBTH]*r;  ,
              fA[i][iBPHI] = flux[i][iBPHI]*r;)
      #endif
      #if HAVE_ENERGY
       fA[i][ENG] = flux[i][ENG]*r2;
      #endif
      #ifdef GLM_MHD
       fA[i][PSI_GLM] = flux[i][PSI_GLM]*r2;
      #endif
       NSCL_LOOP(nv) fA[i][nv] = flux[i][nv]*r2;
    } 

  /* ****************************************************
      Spherical radial direction:

       - initialize rhs with flux differences (I1)
       - add source terms                     (I2)
       - enforce conservation of total angular 
         momentum and/or energy               (I3)
       - add gravity                          (I4)
     **************************************************** */

    #if COMPONENTS >= 2  /* -- need enthalpy for source term computation -- */
     Enthalpy (state->v, h, beg, end);
    #endif
    for (i = beg; i <= end; i++) { 
      r    = x1[i];
      dtdV = dt/dV1[i];
      dtdx = dt/dx1[i];
      r_1  = grid[IDIR].r_1[i];

    /* -----------------------------------------------
       I1. initialize rhs with flux difference
       ----------------------------------------------- */

      rhs[i][RHO] = -dtdV*(fA[i][RHO] - fA[i-1][RHO]);
      EXPAND(
        rhs[i][iMR]   = - dtdV*(fA[i][iMR] - fA[i-1][iMR])
                        - dtdx*(p[i] - p[i-1]);                    ,
        rhs[i][iMTH]  = -dtdV*(fA[i][iMTH]  - fA[i-1][iMTH]);      ,
        rhs[i][iMPHI] = -dtdV*(fA[i][iMPHI] - fA[i-1][iMPHI])*r_1; 
      )
      #if PHYSICS == RMHD
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
      #if HAVE_ENERGY
       rhs[i][ENG] = -dtdV*(fA[i][ENG] - fA[i-1][ENG]);
      #endif

      
       NSCL_LOOP(nv)  rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);
      

    /* ----------------------------------------------------
       I2. Add source terms 
       ---------------------------------------------------- */
  
      v = state->vh[i];

      vel2 = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
      lor2 = 1.0/(1.0 - vel2);
      #if PHYSICS == RMHD
       vB   = EXPAND(v[VX1]*v[BX1], + v[VX2]*v[BX2], + v[VX3]*v[BX3]);
       B2   = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
       wt   = v[RHO]*h[i]*lor2 + B2;
       EXPAND(                                  ,
              mth  = wt*v[iVTH]  - vB*v[iBTH];  ,
              mphi = wt*v[iVPHI] - vB*v[iBPHI];)
      #elif PHYSICS == RHD
       wt   = v[RHO]*h[i]*lor2;
       EXPAND(                 ;     ,
              mth  = wt*v[iVTH];     ,
              mphi = wt*v[iVPHI];)
      #endif

      Sm = EXPAND(  0.0, + mth*v[iVTH], + mphi*v[iVPHI]);
      #if PHYSICS == RMHD 
       Sm += EXPAND(  0.0, - (v[iBTH]/lor2  + vB*v[iVTH])*v[iBTH], 
                           - (v[iBPHI]/lor2 + vB*v[iVPHI])*v[iBPHI]);
      #endif
      rhs[i][iMR] += dt*Sm*r_1;

    /* ----------------------------------------------------
       I4. Include gravity
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
       BodyForceVector(v, g, x1[i], x2[j], x3[k]); 
       rhs[i][iMR] += dt*v[RHO]*g[IDIR];
       #if HAVE_ENERGY
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[IDIR]; 
       #endif
      #endif
    }

  } else if (g_dir == JDIR) {
    double Sm;

  /* ****************************************************
      Spherical meridional direction:
      multiply fluxes by zone-interface area
     **************************************************** */

    r   = x1[i];
    phi = x3[k];
    for (j = beg - 1; j <= end; j++){ 
      s  = grid[JDIR].A[j];
      s2 = s*s;

      fA[j][RHO] = flux[j][RHO]*s;
      EXPAND(fA[j][iMR]   = flux[j][iMR]*s;   ,
             fA[j][iMTH]  = flux[j][iMTH]*s;  ,
             fA[j][iMPHI] = flux[j][iMPHI]*s2;) 
      #if PHYSICS == RMHD
       EXPAND(fA[j][iBR]   = flux[j][iBR]*s;   ,
              fA[j][iBTH]  = flux[j][iBTH]*s;  ,
              fA[j][iBPHI] = flux[j][iBPHI];)
      #endif
      #if HAVE_ENERGY
       fA[j][ENG] = flux[j][ENG]*s;
      #endif
      #ifdef GLM_MHD
       fA[j][PSI_GLM] = flux[j][PSI_GLM]*s;
      #endif
       NSCL_LOOP(nv) fA[j][nv] = flux[j][nv]*s;
    }

  /* ****************************************************
      Spherical meridional direction:

       - initialize rhs with flux differences (J1)
       - add source terms                     (J2)
       - enforce conservation of total angular
         momentum and/or energy               (J3)
       - add gravity                          (J4)
     **************************************************** */
    
    r_1 = 0.5*(x1p[i]*x1p[i] - x1p[i-1]*x1p[i-1])/dV1[i];

    #if COMPONENTS >= 2  /* -- need enthalpy for source term computation -- */
     Enthalpy (state->v, h, beg, end);
    #endif
    for (j = beg; j <= end; j++){
      th   = x2[j];
      dtdV = dt/dV2[j]*r_1;
      dtdx = dt/dx2[j]*r_1;      
      s    = sin(th);
      s_1  = 1.0/s;   
      ct   = grid[JDIR].ct[j];         /* = cot(theta)  */

    /* -----------------------------------------------
       J1. initialize rhs with flux difference
       ----------------------------------------------- */

      rhs[j][RHO] = -dtdV*(fA[j][RHO] - fA[j-1][RHO]);
      EXPAND(
        rhs[j][iMR]   = - dtdV*(fA[j][iMR] - fA[j-1][iMR]);  , 
        rhs[j][iMTH]  = - dtdV*(fA[j][iMTH] - fA[j-1][iMTH])
                        - dtdx*(p[j] - p[j-1]);              , 
        rhs[j][iMPHI] = - dtdV*(fA[j][iMPHI] - fA[j-1][iMPHI])*fabs(s_1);
      )       
      #if PHYSICS == RMHD
       EXPAND(                                                     
         rhs[j][iBR]   = -dtdV*(fA[j][iBR]   - fA[j-1][iBR]);  ,
         rhs[j][iBTH]  = -dtdV*(fA[j][iBTH]  - fA[j-1][iBTH]);  ,
         rhs[j][iBPHI] = -dtdx*(fA[j][iBPHI] - fA[j-1][iBPHI]);
       )
       #ifdef GLM_MHD
        rhs[j][iBTH]    = -dtdx*(flux[j][iBTH] - flux[j-1][iBTH]);
        rhs[j][PSI_GLM] = -dtdV*(fA[j][PSI_GLM] - fA[j-1][PSI_GLM]);
       #endif
      #endif
      #if HAVE_ENERGY
       rhs[j][ENG] = -dtdV*(fA[j][ENG] - fA[j-1][ENG]);
      #endif
      
       NSCL_LOOP(nv)  rhs[j][nv] = -dtdV*(fA[j][nv] - fA[j-1][nv]);
      

    /* ----------------------------------------------------
       J2. Add source terms
       ---------------------------------------------------- */
       
      v    = state->vh[j];
      vel2 = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
      lor2 = 1.0/(1.0 - vel2);
      #if PHYSICS == RMHD
       vB  = EXPAND(v[VX1]*v[BX1], + v[VX2]*v[BX2], + v[VX3]*v[BX3]);
       B2  = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
       wt  = v[RHO]*h[j]*lor2 + B2;
       EXPAND(                                  ,
              mth  = wt*v[iVTH]  - vB*v[iBTH];  ,
              mphi = wt*v[iVPHI] - vB*v[iBPHI];)
      #elif PHYSICS == RHD
       wt   = v[RHO]*h[j]*lor2;
       EXPAND(                 ;      ,
              mth  = wt*v[iVTH];      ,
              mphi = wt*v[iVPHI];)
      #endif

      Sm = EXPAND(  0.0, - mth*v[iVR], + ct*mphi*v[iVPHI]);
      #if PHYSICS == RMHD
      Sm += EXPAND(  0.0, +    (v[iBTH]/lor2  + vB*v[iVTH])*v[iBR],
                          - ct*(v[iBPHI]/lor2 + vB*v[iVPHI])*v[iBPHI]);
      #endif
      rhs[j][iMTH] += dt*Sm*r_1;

    /* ----------------------------------------------------
       J4. Include gravity
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
       BodyForceVector(v, g, x1[i], x2[j], x3[k]);
       rhs[j][iMTH] += dt*v[RHO]*g[JDIR];
       #if HAVE_ENERGY
        rhs[j][ENG] += dt*0.5*(flux[j][RHO] + flux[j-1][RHO])*g[JDIR];
       #endif
      #endif
    }

  } else if (g_dir == KDIR) {

  /* ****************************************************
      Spherical azimuthal direction:

       - initialize rhs with flux differences (K1)
       - add gravity                          (K4)
     **************************************************** */

    r  = x1[i];
    th = x2[j];
    r_1  = 0.5*(x1p[i]*x1p[i] - x1p[i-1]*x1p[i-1])/dV1[i];
    scrh = dt*r_1*dx2[j]/dV2[j];

    for (k = beg; k <= end; k++) {
      phi  = x3[k];
      dtdx = scrh/dx3[k];

    /* ------------------------------------------------
       K1.  initialize rhs with flux difference
       ------------------------------------------------ */

      NVAR_LOOP(nv) rhs[k][nv] = -dtdx*(flux[k][nv] - flux[k-1][nv]);
      rhs[k][iMPHI] -= dtdx*(p[k] - p[k-1]); 

    /* -------------------------------------------------------
       K4. Include gravity
       ------------------------------------------------------- */

      v = state->vh[k];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(v, g, x1[i], x2[j], x3[k]);
       rhs[k][iMPHI] += dt*v[RHO]*g[KDIR];
       #if HAVE_ENERGY
        rhs[k][ENG] += dt*0.5*(flux[k][RHO] + flux[k-1][RHO])*g[KDIR];
       #endif
      #endif
    }
  }
}
#endif  /* GEOMETRY == SPHERICAL */

/* --------------------------------------------------
              Powell's source terms
   -------------------------------------------------- */

  #if PHYSICS == RMHD
   #if DIVB_CONTROL == EIGHT_WAVES
    for (i = beg; i <= end; i++) {
      EXPAND(rhs[i][MX1] += dt*state->src[i][MX1];  ,
             rhs[i][MX2] += dt*state->src[i][MX2];  ,
             rhs[i][MX3] += dt*state->src[i][MX3];)

      EXPAND(rhs[i][BX1] += dt*state->src[i][BX1];  ,
             rhs[i][BX2] += dt*state->src[i][BX2];  ,
             rhs[i][BX3] += dt*state->src[i][BX3];)
      #if HAVE_ENERGY
       rhs[i][ENG] += dt*state->src[i][ENG];
      #endif
    }
   #endif
  #endif

/* -------------------------------------------------
            Extended GLM source terms
   ------------------------------------------------- */

  #if (defined GLM_MHD) && (GLM_EXTENDED == YES)
   print1 ("! RightHandSide(): Extended GLM source terms not defined for RMHD\n");
   QUIT_PLUTO(1);
  #endif

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
