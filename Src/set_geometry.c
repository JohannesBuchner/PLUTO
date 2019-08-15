/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Initialize geometry-dependent grid quantities.

  Compute grid quantities (such as interface areas, volumes, 
  centroid of volume, etc..) that depend on the geometry.

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 13, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void MakeGeometry (Grid *grid)
/*!
 *
 * \param [in,out] grid  Pointer to an array of Grid structures;
 *
 *********************************************************************** */
{
  int     i, j, k, idim;
  int     iend = grid->lend[IDIR] + grid->nghost[IDIR];
  int     jend = grid->lend[JDIR] + grid->nghost[JDIR];
  int     kend = grid->lend[KDIR] + grid->nghost[KDIR];
  int     nx1_tot = grid->np_tot[IDIR];
  int     nx2_tot = grid->np_tot[JDIR];
  int     nx3_tot = grid->np_tot[KDIR];
  double  dVr, dmu, xL, xR;

  double *x1 = grid->x[IDIR], *dx1 = grid->dx[IDIR];
  double *x2 = grid->x[JDIR], *dx2 = grid->dx[JDIR];
  double *x3 = grid->x[KDIR], *dx3 = grid->dx[KDIR];

  double *x1p = grid->xr[IDIR], *x1m = grid->xl[IDIR];
  double *x2p = grid->xr[JDIR], *x2m = grid->xl[JDIR];
  double *x3p = grid->xr[KDIR], *x3m = grid->xl[KDIR];

  double ***Ax1, ***Ax2, ***Ax3;

/* --------------------------------------------------------------
   0. Allocate memory for the Grid structue.
      All values are defined at the cell center 
      with the exception of the area element which is defined on a
      staggered mesh and therefore starts at [-1].
   ----------------------------------------------------------- */

  grid->dV      = ARRAY_3D(nx3_tot, nx2_tot, nx1_tot, double);
  grid->A[IDIR] = ArrayBox( 0, nx3_tot-1,  0, nx2_tot-1, -1, nx1_tot-1);
  grid->A[JDIR] = ArrayBox( 0, nx3_tot-1, -1, nx2_tot-1,  0, nx1_tot-1);
  grid->A[KDIR] = ArrayBox(-1, nx3_tot-1,  0, nx2_tot-1,  0, nx1_tot-1);

  grid->dx_dl[IDIR] = ARRAY_2D(nx2_tot, nx1_tot, double);
  grid->dx_dl[JDIR] = ARRAY_2D(nx2_tot, nx1_tot, double);
  grid->dx_dl[KDIR] = ARRAY_2D(nx2_tot, nx1_tot, double);
  
  grid->rt  = ARRAY_1D(grid->np_tot[IDIR], double);
  grid->sp  = ARRAY_1D(grid->np_tot[JDIR], double);
  grid->s   = ARRAY_1D(grid->np_tot[JDIR], double);
  grid->dmu = ARRAY_1D(grid->np_tot[JDIR], double);
  for (idim = 0; idim < 3; idim++) {
    grid->xgc[idim]     = ARRAY_1D(grid->np_tot[idim], double);
    grid->inv_dx[idim]  = ARRAY_1D(grid->np_tot[idim], double);
    grid->inv_dxi[idim] = ARRAY_1D(grid->np_tot[idim], double);
  }

/* ----------------------------------------------------------
   1a. Compute positions arrays in the X1 (IDIR) direction
   ---------------------------------------------------------- */

  for (i = 0; i <= iend; i++){
    xL = x1m[i];
    xR = x1p[i];
    #if GEOMETRY == CARTESIAN
    grid->xgc[IDIR][i] = x1[i];
    grid->rt[i]        = x1[i];
    #elif GEOMETRY == CYLINDRICAL || GEOMETRY == POLAR
    grid->xgc[IDIR][i] = x1[i] + dx1[i]*dx1[i]/(12.0*x1[i]); 
    grid->rt[i]        = x1[i];
    #elif GEOMETRY == SPHERICAL
    grid->xgc[IDIR][i] = x1[i] + 2.0*x1[i]*dx1[i]*dx1[i]/
                                (12.0*x1[i]*x1[i] + dx1[i]*dx1[i]);
    grid->rt[i]  = (xR*xR*xR - xL*xL*xL)/(xR*xR - xL*xL)/1.5;
    #endif
  }

/* ----------------------------------------------------------
   1b. Compute positions arrays in the X2 (JDIR) direction
   ---------------------------------------------------------- */

  for (j = 0; j <= jend; j++){
    xL = x2m[j];
    xR = x2p[j];

    #if GEOMETRY != SPHERICAL
    grid->xgc[JDIR][j] = x2[j];
    #else
    grid->xgc[JDIR][j]  = sin(xR) - sin(xL)+ xL*cos(xL) - xR*cos(xR);
    grid->xgc[JDIR][j] /= cos(xL) - cos(xR);
    grid->sp[j]         = fabs(sin(xR));
    grid->s[j]          = fabs(sin(x2[j]));
    grid->dmu[j]        = fabs(cos(xL) - cos(xR));
    #endif
  }

/* ------------------------------------------------------------
   1c. Compute positions arrays in the X3 (KDIR) direction
   ------------------------------------------------------------ */

  for (k = 0; k <= kend; k++){
    grid->xgc[KDIR][k] = x3[k];
  }

/* ------------------------------------------------------------
   2. Compute volumes 
   ------------------------------------------------------------ */

  for (k = 0; k <= kend; k++){
  for (j = 0; j <= jend; j++){
  for (i = 0; i <= iend; i++){
    #if GEOMETRY == CARTESIAN
    grid->dV[k][j][i]  = D_EXPAND(dx1[i], *dx2[j], *dx3[k]);  /* = dx*dy*dz */
    #elif GEOMETRY == CYLINDRICAL
    dVr = fabs(x1[i])*dx1[i];
    grid->dV[k][j][i]  = D_EXPAND(dVr, *dx2[j], *1.0);        /* = |r|*dr*dz */
    #elif GEOMETRY == POLAR
    dVr = fabs(x1[i])*dx1[i];
    grid->dV[k][j][i]  = D_EXPAND(dVr, *dx2[j], *dx3[k]);     /* = |r|*dr*dphi*dz */
    #elif GEOMETRY == SPHERICAL
    dVr = fabs(x1p[i]*x1p[i]*x1p[i] - x1m[i]*x1m[i]*x1m[i])/3.0;
    dmu = fabs(cos(x2m[j]) - cos(x2p[j]));
    grid->dV[k][j][i]  = D_EXPAND(dVr, *dmu, *dx3[k]);        /* = dVr*dmu*dphi*/
    #endif
  }}}

/* ------------------------------------------------------------
   3a. Compute area in the x1-direction
   ------------------------------------------------------------ */

  Ax1 = grid->A[IDIR];
  Ax2 = grid->A[JDIR];
  Ax3 = grid->A[KDIR];

  for (k =  0; k <= kend; k++){
  for (j =  0; j <= jend; j++){
  for (i = -1; i <= iend; i++){
    #if GEOMETRY == CARTESIAN
    Ax1[k][j][i] = D_EXPAND(1.0, *dx2[j], *dx3[k]);         /* = dy*dz */ 
    #elif GEOMETRY == CYLINDRICAL
    if (i == -1) {
      Ax1[k][j][i] = D_EXPAND(fabs(x1m[0]), *dx2[j], *1.0); /* = rp*dz */
    }else{
      Ax1[k][j][i] = D_EXPAND(fabs(x1p[i]), *dx2[j], *1.0); /* = rp*dz */
    }
    #elif GEOMETRY == POLAR
    if (i == -1) {
      Ax1[k][j][i] = D_EXPAND(fabs(x1m[0]), *dx2[j], *dx3[k]); /* = rp*dphi*dz */
    }else{
      Ax1[k][j][i] = D_EXPAND(fabs(x1p[i]), *dx2[j], *dx3[k]); /* = rp*dphi*dz */
    }
    #elif GEOMETRY == SPHERICAL
    dmu = fabs(cos(x2m[j]) - cos(x2p[j]));
    if (i == -1) {
      Ax1[k][j][i] = D_EXPAND(x1m[0]*x1m[0], *dmu, *dx3[k]); /* = rp^2*dmu*dphi */
    }else{
      Ax1[k][j][i] = D_EXPAND(x1p[i]*x1p[i], *dmu, *dx3[k]); /* = rp^2*dmu*dphi */
    }
    #endif
  }}}

/* ------------------------------------------------------------
   3b. Compute area in the x2-direction
   ------------------------------------------------------------ */

  for (k =  0; k <= kend; k++){
  for (j = -1; j <= jend; j++){
  for (i =  0; i <= iend; i++){
    #if GEOMETRY == CARTESIAN
    Ax2[k][j][i] = D_EXPAND(dx1[i], *1.0, *dx3[k]);        /* = dx*dz */
    #elif GEOMETRY == CYLINDRICAL
    Ax2[k][j][i] = D_EXPAND(fabs(x1[i]), *dx1[i], *1.0);   /* = r*dr */
    #elif GEOMETRY == POLAR
    Ax2[k][j][i] = D_EXPAND(dx1[i], *1.0, *dx3[k]);        /* = dr*dz */    
    #elif GEOMETRY == SPHERICAL
    if (j == -1){
      Ax2[k][j][i] = D_EXPAND(x1[i]*dx1[i], *fabs(sin(x2m[0])), *dx3[k]); /* = r*dr*sin(thp)*dphi */
    }else{
      Ax2[k][j][i] = D_EXPAND(x1[i]*dx1[i], *fabs(sin(x2p[j])), *dx3[k]); /* = r*dr*sin(thp)*dphi */
    }
    #endif
  }}}

/* ------------------------------------------------------------
   3c. Compute area in the x3-direction
   ------------------------------------------------------------ */

  for (k = -1; k <= kend; k++){
  for (j =  0; j <= jend; j++){
  for (i =  0; i <= iend; i++){
    #if GEOMETRY == CARTESIAN
    Ax3[k][j][i] = D_EXPAND(dx1[i], *dx2[j], *1.0);          /* = dx*dy */ 
    #elif GEOMETRY == CYLINDRICAL
    Ax3[k][j][i] = 1.0;   /* No 3rd direction in cylindrical coords */
    #elif GEOMETRY == POLAR
    Ax3[k][j][i] = D_EXPAND(x1[i]*dx1[i], *dx2[j], *1.0);   /* = r*dr*dphi */        
    #elif GEOMETRY == SPHERICAL
    Ax3[k][j][i] = D_EXPAND(x1[i]*dx1[i], *dx2[j], *1.0);   /* = r*dr*dth */        
    #endif
  }}}

/* ------------------------------------------------------------
   4a. Compute dx/dl in the x1-direction
   ------------------------------------------------------------ */

  for (k =  0; k <= kend; k++){
  for (j =  0; j <= jend; j++){
  for (i =  0; i <= iend; i++){
    #if GEOMETRY == CARTESIAN
    grid->dx_dl[IDIR][j][i] = 1.0;
    grid->dx_dl[JDIR][j][i] = 1.0;
    grid->dx_dl[KDIR][j][i] = 1.0;
    #elif GEOMETRY == CYLINDRICAL
    grid->dx_dl[IDIR][j][i] = 1.0;
    grid->dx_dl[JDIR][j][i] = 1.0;
    #elif GEOMETRY == POLAR
    grid->dx_dl[IDIR][j][i] = 1.0;
    grid->dx_dl[JDIR][j][i] = 1.0/x1[i];
    grid->dx_dl[KDIR][j][i] = 1.0;
    #elif GEOMETRY == SPHERICAL
    grid->dx_dl[IDIR][j][i] = 1.0;
    grid->dx_dl[JDIR][j][i] = 1.0/grid->rt[i];
    grid->dx_dl[KDIR][j][i] = dx2[j]/(grid->rt[i]*grid->dmu[j]);
    #endif
  }}}
  
/* ---------------------------------------------------------
   5. Compute and store the reciprocal of cell spacing
      between interfaces (inv_dx) and cell centers (inv_dxi)
   --------------------------------------------------------- */

  for (idim = 0; idim < DIMENSIONS; idim++){
    for (i = 0; i < grid->np_tot[idim]; i++) {
      grid->inv_dx[idim][i] = 1.0/(grid->dx[idim][i]);
    }

    for (i = 0; i < grid->np_tot[idim]-1; i++) {
      grid->inv_dxi[idim][i] = 2.0/(grid->dx[idim][i] + grid->dx[idim][i+1]);
    }
  }  
}

/* ********************************************************************** */
double Length_1 (int i, int j, int k, Grid *grid)
/*
 *
 *
 *
 *
 *
 ************************************************************************ */
{
  return (grid->dx[IDIR][i]);
}

/* ********************************************************************** */
double Length_2 (int i, int j, int k, Grid *grid)
/*
 *
 *
 *
 *
 *
 ************************************************************************ */
{
#if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
  return (grid->dx[JDIR][j]);
#elif GEOMETRY == POLAR ||  GEOMETRY == SPHERICAL
  return (fabs(grid->xgc[IDIR][i])*grid->dx[JDIR][j]);
#endif
}

/* ********************************************************************** */
double Length_3 (int i, int j, int k, Grid *grid)
/*
 *
 *
 *
 *
 *
 ************************************************************************ */
{
#if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
  return grid->dx[KDIR][k];
#elif GEOMETRY == CYLINDRICAL
  return fabs(grid->xgc[IDIR][i])*grid->dx[KDIR][k];
#elif GEOMETRY == SPHERICAL
  return fabs(grid->xgc[IDIR][i]*sin(grid->xgc[JDIR][j]))*grid->dx[KDIR][k];
#endif
}

/* ******************************************************************* */
double *GetInverse_dl (const Grid *grid)
/*
 *
 *  Return an array of (inverse) physical cell lengths in the
 *  direction given by g_dir.
 *  For spherical coordinates, for instance this will be
 *
 *    {dr}_i                      if g_dir == IDIR
 *    {r_i*dtheta}_j              if g_dir == JDIR
 *    {r_i*sin(theta_j)*dphi}_k   if g_dir == KDIR
 *   
 *
 ********************************************************************* */
{
#if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL

  return grid->inv_dx[g_dir];

#elif GEOMETRY == POLAR

  if (g_dir != JDIR){
    return grid->inv_dx[g_dir];
  }else{
    int    j;
    double r_1;
    static double *inv_dl;
   
    if (inv_dl == NULL) {
     #ifdef CHOMBO
      inv_dl = ARRAY_1D(NX2_MAX, double);
     #else
      inv_dl = ARRAY_1D(NX2_TOT, double);
     #endif
    }
    r_1 = 1.0/grid->x[IDIR][g_i];
    JTOT_LOOP(j) inv_dl[j] = grid->inv_dx[JDIR][j]*r_1;
    return inv_dl;
  }

#elif GEOMETRY == SPHERICAL

  int    j, k;
  double r_1, s;
  static double *inv_dl2, *inv_dl3;

  if (inv_dl2 == NULL) {
   #ifdef CHOMBO
    inv_dl2 = ARRAY_1D(NX2_MAX, double);
    inv_dl3 = ARRAY_1D(NX3_MAX, double);
   #else
    inv_dl2 = ARRAY_1D(NX2_TOT, double);
    inv_dl3 = ARRAY_1D(NX3_TOT, double);
   #endif
  }

  if (g_dir == IDIR){
    return grid->inv_dx[IDIR];
  }else if (g_dir == JDIR) {
    r_1 = 1.0/grid->x[IDIR][g_i];
    JTOT_LOOP(j) inv_dl2[j] = grid->inv_dx[JDIR][j]*r_1;
    return inv_dl2;
  }else if (g_dir == KDIR){
    r_1 = 1.0/grid->x[IDIR][g_i];
    s   = grid->x[JDIR][g_j];
    s   = sin(s);
    KTOT_LOOP(k) inv_dl3[k] = grid->inv_dx[KDIR][k]*r_1/s;
    return inv_dl3;
  }

#endif

  return NULL;
}

