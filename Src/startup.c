/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Loop on the computational cells to assign initial conditions.

  This function is called anyway, even if restart from file is enabled.
  This is useful to initialized a number of global variables and/or 
  user-defined parameters by calling Init().

  \author A. Mignone (mignone@ph.unito.it)
          B. Vaidya

  \Last Change  D. Mukherjee
  \date         May 25, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
double FieldAverage(double, double, Grid *, int, int);

/* ********************************************************************* */
void Startup (Data *d, Grid *grid)
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
  double us[256], u_av[256], bp[3], bm[3]; 
  double scrh;

/* ------------------------------------------------------
   0. Initialize function by allocating memory
      and setting labels.
   ------------------------------------------------------ */

  print ("> Assigning initial conditions (Startup) ...\n");

  if (uprim == NULL){
    uprim = ARRAY_2D(NMAX_POINT, NVAR, double);
    ucons = ARRAY_2D(NMAX_POINT, NVAR, double); 
  }

  EXPAND(MXn = VXn = VX1;  ,
         MXt = VXt = VX2;  ,
         MXb = VXb = VX3;)

  #if PHYSICS == MHD || PHYSICS == RMHD
   EXPAND(BXn = BX1;  ,
          BXt = BX2;  ,
          BXb = BX3;)
  #endif

/* ------------------------------------------------------
   1. Assign initial conditions point by point using
      primitive variables.
   ------------------------------------------------------ */

  KTOT_LOOP(k) { 
  JTOT_LOOP(j) { 
  ITOT_LOOP(i) { 

  /* -- Define coordinates -- */

    #if GEOMETRY == CYLINDRICAL
    x1 = grid->xgc[IDIR][i]; x1p = grid->xr[IDIR][i]; dx1 = grid->dx[IDIR][i];
    x2 = grid->xgc[JDIR][j]; x2p = grid->xr[JDIR][j]; dx2 = grid->dx[JDIR][j];
    x3 = grid->xgc[KDIR][k]; x3p = grid->xr[KDIR][k]; dx3 = grid->dx[KDIR][k];
    #else
    x1 = grid->x[IDIR][i]; x1p = grid->xr[IDIR][i]; dx1 = grid->dx[IDIR][i];
    x2 = grid->x[JDIR][j]; x2p = grid->xr[JDIR][j]; dx2 = grid->dx[JDIR][j];
    x3 = grid->x[KDIR][k]; x3p = grid->xr[KDIR][k]; dx3 = grid->dx[KDIR][k];
    #endif
    
  /* -- Reset arrays -- */

    NVAR_LOOP(nv) d->Vc[nv][k][j][i] = u_av[nv] = 0.0;

    #ifdef GLM_MHD
    u_av[PSI_GLM] = us[PSI_GLM] = 0.0;
    #ifdef PHI_GLM
    u_av[PHI_GLM] = us[PHI_GLM] = 0.0;
    #endif
    #endif

    #if INITIAL_SMOOTHING == YES

    for (ksub = 0; ksub < nsub; ksub++){ 
    for (jsub = 0; jsub < nsub; jsub++){ 
    for (isub = 0; isub < nsub; isub++){ 
            
      x1s = x1 + (double)(1.0 - nsub + 2.0*isub)/(double)(2.0*nsub)*dx1;
      x2s = x2 + (double)(1.0 - nsub + 2.0*jsub)/(double)(2.0*nsub)*dx2;
      x3s = x3 + (double)(1.0 - nsub + 2.0*ksub)/(double)(2.0*nsub)*dx3;

      Init (us, x1s, x2s, x3s);
      NVAR_LOOP(nv) u_av[nv] += us[nv]/(double)(nsub*nsub*nsub);
    }}}

    #else

    Init (u_av, x1, x2, x3);
  
    #endif

    for (nv = NVAR; nv--;  ) d->Vc[nv][k][j][i] = u_av[nv];

  /* -- Initialize cell-centered vector potential 
        (only for output purposes)                -- */

    #if (PHYSICS == MHD || PHYSICS == RMHD) && (ASSIGN_VECTOR_POTENTIAL == YES)
  /* ------------------------------------------------------
     Store vector potential into array for later
     differentiation. 
     - Note that staggered MHD requires A[]
     to be staggered (Ax1[0++], Ax2[+0+], Ax3[++0]).
     - All 3 components must be given here. 
     - Vector potential is always staggered.
     ------------------------------------------------------ */
    Init (u_av, x1, x2p, x3p);
    d->Ax1[k][j][i] = u_av[AX1];  

    Init (u_av, x1p, x2, x3p);
    d->Ax2[k][j][i] = u_av[AX2];

    Init (u_av, x1p, x2p, x3);
    d->Ax3[k][j][i] = u_av[AX3];

    #endif

 /* -------------------------------------------------------
     Assign staggered components directly from 
     the init.c  if ASSIGN_VECTOR_POTENTIAL not defined 
    NOTE: in N dimensions only N components are assigned 
             through this call.
    ------------------------------------------------------- */    

    #if (defined STAGGERED_MHD) && (ASSIGN_VECTOR_POTENTIAL == NO)
    D_EXPAND(
             Init (u_av, x1p, x2, x3);
             d->Vs[BX1s][k][j][i] = u_av[BX1];       ,

             Init (u_av, x1, x2p, x3);
             d->Vs[BX2s][k][j][i] = u_av[BX2];       ,

             Init (u_av, x1, x2, x3p);
             d->Vs[BX3s][k][j][i] = u_av[BX3];)
    #endif

  }}}

/* --------------------------------------------------------
   2. Call Init_Domain() to assign primitive variables
      by looping over computational cells.
      This is new in PLUTO 4.3 and provides an
      alternative to the pointwise initialization.
   -------------------------------------------------------- */

  InitDomain(d, grid);

 /*----------------------------------------------------------------- 
   3. Assign magnetic field from vector potential.
      Assignments run from array elements 1 to NX#_TOT-2,
      because the magnetic field needs extra cells for derivatives.
      The extremal points are to be assigned at the 
      boundary conditions or domain swap.
  ------------------------------------------------------------------*/

  #if (PHYSICS ==  MHD || PHYSICS == RMHD) && ASSIGN_VECTOR_POTENTIAL == YES
  #ifdef STAGGERED_MHD
  for (k = KOFFSET; k < NX3_TOT-KOFFSET; k++) { 
  for (j = JOFFSET; j < NX2_TOT-JOFFSET; j++){
  for (i = IOFFSET; i < NX1_TOT-IOFFSET; i++){ 
      VectorPotentialDiff(bp, d, i, j, k, grid);
      for (nv = 0; nv < DIMENSIONS; nv++) d->Vs[nv][k][j][i] = bp[nv];
  }}}
  #else

/*Note: Start from 2 as we call bm, which passes i-1 
       to vec_pot_diff which calls i-2 */
  for (k = 2*KOFFSET; k < NX3_TOT-KOFFSET; k++){  
  for (j = 2*JOFFSET; j < NX2_TOT-JOFFSET; j++){
  for (i = 2*IOFFSET; i < NX1_TOT-IOFFSET; i++){ 

      VectorPotentialDiff(bp, d, i, j, k, grid);
      VectorPotentialDiff(bm, d, i-1, j, k, grid);
      d->Vc[BX1][k][j][i] = FieldAverage(bp[0], bm[0], grid, IDIR, i);

      #if DIMENSIONS == 2
      VectorPotentialDiff(bm, d, i, j-1, k, grid);
      d->Vc[BX2][k][j][i] = FieldAverage(bp[1], bm[1], grid, JDIR, j);
      #endif

      #if DIMENSIONS == 3    
      VectorPotentialDiff(bm, d, i, j, k-1, grid);
      d->Vc[BX3][k][j][i] = FieldAverage(bp[2], bm[2], grid, KDIR, k);
      #endif
  }}}
  #endif /* STAGGERED_MHD */
  #endif

/* --------------------------------------------------------
   4. Initialize Forced-Turbulence module
   -------------------------------------------------------- */

  #if FORCED_TURB == YES
    struct ForcedTurb *Ft;
    Ft = d->Ft;
    ForcedTurb_Init(Ft);
    print("> Initialized Forced Turbulence %d Modes\n",Ft->NModes);
    ForcedTurb_OUNoiseInit(Ft->OUPhases, 6*(Ft->NModes), Ft->OUVar);
  #endif

/* --------------------------------------------------------
   5. Compute cell-centered magnetic field by simple
      arithmetic averaging.
      Useful only for saving the first output, since the
      average will be re-computed anyway at the beginning
      of the computation.
   -------------------------------------------------------- */

#ifdef STAGGERED_MHD
  DOM_LOOP(k,j,i){
    D_EXPAND( 
      d->Vc[BX1][k][j][i] = 0.5*(d->Vs[BX1s][k][j][i] + d->Vs[BX1s][k][j][i-1]); , 
      d->Vc[BX2][k][j][i] = 0.5*(d->Vs[BX2s][k][j][i] + d->Vs[BX2s][k][j-1][i]); ,
      d->Vc[BX3][k][j][i] = 0.5*(d->Vs[BX3s][k][j][i] + d->Vs[BX3s][k-1][j][i]);
    )
  /*  Use the new fieldAverage routine to be consistent? */
/*   
    D_EXPAND(
     d->Vc[BX1][k][j][i] = FieldAverage(d->Vs[BX1s][k][j][i], d->Vs[BX1s][k][j][i-1], 
                                        grid, IDIR, i); ,
     d->Vc[BX2][k][j][i] = FieldAverage(d->Vs[BX2s][k][j][i], d->Vs[BX2s][k][j-1][1], 
                                        grid, JDIR, j); ,
     d->Vc[BX3][k][j][i] = FieldAverage(d->Vs[BX3s][k][j][i], d->Vs[BX3s][k-1][j][1], 
                                        grid, KDIR, k); )
*/                                        
  }
#endif

/* ------------------------------------------------------
   6.  Check if values have physical meaning.
   ------------------------------------------------------ */

#if PHYSICS != ADVECTION
  KDOM_LOOP(k){ x3 = grid->x[KDIR][k];
  JDOM_LOOP(j){ x2 = grid->x[JDIR][j];
  IDOM_LOOP(i){ x1 = grid->x[IDIR][i];

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

/* ------------------------------------------------------
   7. Set boundary conditions.
   ------------------------------------------------------ */

  Boundary (d, -1, grid);
}

/* ********************************************************************* */
double FieldAverage(double bp, double bm, Grid *grid, int dir, int l)
/*!
 *  Computes the volume average of staggered magnetic fields to return
 *  cell-centred values. 
 *  bp: B component at right staggered face (e.g. Bx[i, j, k])
 *  bm: B component at left staggered face  (e.g. Bx[i-1, j, k])
 *  dir: direction of staggered. options: IDIR, JDIR, KDIR
 *  l:  integer for cell along direction of staggered.
 *  Grid: PLUTO grid structure.
 *  return: average field.

 *********************************************************************** */
{
  double bavg;
 #if (GEOMETRY == CYLINDRICAL) || (GEOMETRY == POLAR) || (GEOMETRY == SPHERICAL)  
  double rp, rm;
 #endif
 #if GEOMETRY == SPHERICAL 
  double thp, thm, sp, sm, dx, dV;
 #endif
 
  if ((dir != IDIR) && (dir != JDIR) && (dir != KDIR)) {
     print("! Wrong dir in FieldAvergae. Abort! \n");
     QUIT_PLUTO(1);
     } 

 #if   GEOMETRY == CARTESIAN
  bavg = 0.5*(bp + bm);

 #elif GEOMETRY == CYLINDRICAL 
  if (dir == IDIR) {
     rp   = grid->xr[dir][l];
     rm   = grid->xl[dir][l];
     bavg = (rp*bp + rm*bm)/(rp + rm);              
     } else bavg = 0.5*(bp + bm);

 #elif GEOMETRY == POLAR /* Note: Polar is same as cylindrical.
                                  Creating separate entries if
                                  changes are needed in future. */
  if (dir == IDIR) {
     rp = grid->xr[dir][l];
     rm = grid->xl[dir][l];
     bavg = (rp*bp + rm*bm)/(rp + rm);              
     } else bavg = 0.5*(bp + bm);

 #elif GEOMETRY == SPHERICAL
  if (dir == IDIR) {
     rp   = grid->xr[dir][l];
     rm   = grid->xl[dir][l];
     dx   = grid->dx[dir][l];
     dV   = (rp*rp*rp - rm*rm*rm)/3.;
     bavg = 0.5*(rp*rp*bp + rm*rm*bm)*dx/dV;
     } 
     else if (dir == JDIR){
     thp  = grid->xr[dir][l];     
     thm  = grid->xl[dir][l];
     sp   = fabs(sin(thp));
     sm   = fabs(sin(thm));
     dx   = grid->dx[dir][l];
     dV   = fabs(cos(thm) - cos(thp)); 
     bavg = 0.5*(sp*bp + sm*bm)*dx/dV;
     } else bavg = 0.5*(bp + bm);
  #endif

 return bavg;

}

