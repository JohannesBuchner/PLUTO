/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Loop on the computational cells to assign initial conditions.

  This function is called anyway, even if restart from file is enabled.
  This is useful to initialized a number of global variables and/or 
  user-defined parameters by calling Init().

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Startup (Data *d, Grid *G)
/*! 
 *
 *
 *
 *
 *********************************************************************** */
{
  int i, j, k;
  int isub, jsub, ksub, nsub = 5;
  int nv,  l_convert;
  static double **ucons, **uprim;
  double x1,  x2,  x3;
  double x1p, x2p, x3p;
  double x1s, x2s, x3s;
  double dx1, dx2, dx3;
  double us[256], u_av[256], b[3]; 
  double scrh;
  struct GRID *GX, *GY, *GZ;

  GX = G;
  GY = G + 1;
  GZ = G + 2;

  if (uprim == NULL){
    uprim = ARRAY_2D(NMAX_POINT, NVAR, double);
    ucons = ARRAY_2D(NMAX_POINT, NVAR, double); 
  }

/* ----  set labels  ---- */

  EXPAND(MXn = VXn = VX1;  ,
         MXt = VXt = VX2;  ,
         MXb = VXb = VX3;)

  #if PHYSICS == MHD || PHYSICS == RMHD
   EXPAND(BXn = BX1;  ,
          BXt = BX2;  ,
          BXb = BX3;)
  #endif

/* ----------------------------------------------------------
    If the -input-data-file command line argument is given, 
    try to read initial conditions from external files.
   ---------------------------------------------------------- */
    
  print1 ("> Assigning initial conditions (Startup) ...\n");

/* --------------------------------------------------------------
                    Assign initial conditions   
   -------------------------------------------------------------- */

  KTOT_LOOP(k) { 
  JTOT_LOOP(j) { 
  ITOT_LOOP(i) { 

    #if GEOMETRY == CYLINDRICAL
     x1 = GX->xgc[i]; x1p = GX->xr[i]; dx1 = GX->dx[i];
     x2 = GY->xgc[j]; x2p = GY->xr[j]; dx2 = GY->dx[j];
     x3 = GZ->xgc[k]; x3p = GZ->xr[k]; dx3 = GZ->dx[k];
    #else
     x1 = GX->x[i]; x1p = GX->xr[i]; dx1 = GX->dx[i];
     x2 = GY->x[j]; x2p = GY->xr[j]; dx2 = GY->dx[j];
     x3 = GZ->x[k]; x3p = GZ->xr[k]; dx3 = GZ->dx[k];
    #endif
    
    for (nv = NVAR; nv--; ) d->Vc[nv][k][j][i] = u_av[nv] = 0.0;

/*  ----------------------------------------------------------------
                Compute volume averages      
    ---------------------------------------------------------------- */

    #ifdef PSI_GLM
     u_av[PSI_GLM] = us[PSI_GLM] = 0.0;
    #endif

    #if INITIAL_SMOOTHING == YES

     for (ksub = 0; ksub < nsub; ksub++){ 
     for (jsub = 0; jsub < nsub; jsub++){ 
     for (isub = 0; isub < nsub; isub++){ 
            
       x1s = x1 + (double)(1.0 - nsub + 2.0*isub)/(double)(2.0*nsub)*dx1;
       x2s = x2 + (double)(1.0 - nsub + 2.0*jsub)/(double)(2.0*nsub)*dx2;
       x3s = x3 + (double)(1.0 - nsub + 2.0*ksub)/(double)(2.0*nsub)*dx3;

       Init (us, x1s, x2s, x3s);
       for (nv = 0; nv < NVAR; nv++) {
         u_av[nv] += us[nv]/(double)(nsub*nsub*nsub);
       }
     }}}

    #else

     Init (u_av, x1, x2, x3);
  
    #endif

    for (nv = NVAR; nv--;  ) d->Vc[nv][k][j][i] = u_av[nv];

/* -----------------------------------------------------
        Initialize cell-centered vector potential 
        (only for output purposes)
   ----------------------------------------------------- */

    #if PHYSICS == MHD || PHYSICS == RMHD
     #if UPDATE_VECTOR_POTENTIAL == YES
      D_EXPAND(                     ,
        d->Ax3[k][j][i] = u_av[AX3];  ,
        d->Ax1[k][j][i] = u_av[AX1];  
        d->Ax2[k][j][i] = u_av[AX2];)
     #endif
    #endif

 /* -------------------------------------------------------------
     Assign staggered components;
     If a vector potential is used (ASSIGN_VECTOR_POTENTIAL == YES),
     use the STAGGERED_INIT routine;
     otherwise assign staggered components directly from 
     the init.c and ignore the vector potential.

     NOTE: in N dimensions only N components are assigned 
             through this call.
    ------------------------------------------------------------- */    

    #if PHYSICS == MHD || PHYSICS == RMHD
     #if ASSIGN_VECTOR_POTENTIAL == YES
      VectorPotentialDiff(b, i, j, k, G);

      #ifdef STAGGERED_MHD
       for (nv = 0; nv < DIMENSIONS; nv++) {
         d->Vs[nv][k][j][i] = b[nv];
       }
      #else
       for (nv = 0; nv < DIMENSIONS; nv++) {
         d->Vc[BX1+nv][k][j][i] = b[nv];
       }
      #endif

     #else

      #ifdef STAGGERED_MHD
       D_EXPAND(
         Init (u_av, x1p, x2, x3);
         d->Vs[BX1s][k][j][i] = u_av[BX1];       ,

         Init (u_av, x1, x2p, x3);
         d->Vs[BX2s][k][j][i] = u_av[BX2];       ,

         Init (u_av, x1, x2, x3p);
         d->Vs[BX3s][k][j][i] = u_av[BX3];
       )
      #endif
     #endif  /* ASSIGN_VECTOR_POTENTIAL */
    #endif /* PHYSICS == MHD || PHYSICS == RMHD */

  }}}

  #ifdef STAGGERED_MHD

  /* ---------------------------------------------------
       Compute cell-centered magnetic field  
       by simple arithmetic averaging. This 
       is useful only for saving the first 
       output, since the average will be 
       re-computed at the beginning of the computation.
     --------------------------------------------------- */

   DOM_LOOP(k,j,i){
     D_EXPAND( 
       d->Vc[BX1][k][j][i] = 0.5*(d->Vs[BX1s][k][j][i] + d->Vs[BX1s][k][j][i-1]); , 
       d->Vc[BX2][k][j][i] = 0.5*(d->Vs[BX2s][k][j][i] + d->Vs[BX2s][k][j-1][i]); ,
       d->Vc[BX3][k][j][i] = 0.5*(d->Vs[BX3s][k][j][i] + d->Vs[BX3s][k-1][j][i]);
     )
   }

  #endif

/* --------------------------------------------------------------------
         Check if values have physical meaning...
   -------------------------------------------------------------------- */

#if PHYSICS != ADVECTION
  KDOM_LOOP(k){ x3 = GZ->x[k];
  JDOM_LOOP(j){ x2 = GY->x[j];
  IDOM_LOOP(i){ x1 = GX->x[i];

    for (nv = NVAR; nv--;  ) us[nv] = d->Vc[nv][k][j][i];

    if (us[RHO] <= 0.0) {
      print ("! Startup(): negative density, zone [%f, %f, %f]\n", x1,x2,x3);
      QUIT_PLUTO(1);
    }
    #if HAVE_ENERGY
     if (us[PRS] <= 0.0) {
       print ("! Startup(): negative pressure, zone [%f, %f, %f]\n",x1,x2,x3);
       QUIT_PLUTO(1);
     }
    #endif

    #if (PHYSICS == RHD || PHYSICS == RMHD) 
     scrh = EXPAND(us[VX1]*us[VX1], + us[VX2]*us[VX2], + us[VX3]*us[VX3]);
     if (scrh >= 1.0){
       print ("! Startup(): total velocity exceeds 1\n"); 
       QUIT_PLUTO(1);
     }
    #endif
  }}}
#endif

/* --------------------------------------------------------------------
     Convert primitive variables to conservative ones
   -------------------------------------------------------------------- */

  Boundary (d, -1, G);
  
/* --------------------------------------------------------------------
                    Convert to conservative
   -------------------------------------------------------------------- */
/*
  KDOM_LOOP(k) {
  JDOM_LOOP(j){
    IDOM_LOOP(i) {
    for (nv = NVAR; nv--;  ) {
      uprim[i][nv] = d->Vc[nv][k][j][i];
    }}
    PrimToCons(uprim, d->Uc[k][j], IBEG, IEND);
  }}
*/
}

