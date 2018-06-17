/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Update staggered magnetic field.

  Update face-centered magnetic field in the constrained transport 
  formulation using a discrete version of Stoke's theorem.
  The update consists of a single Euler step:
  \f[
    \mathtt{d->Vs} = \mathtt{Bs} + \Delta t R
  \f]
  where \c d->Vs is the main staggered array used by PLUTO, 
  \c Bs is the magnetic field to be updated and \c R is the 
  right hand side already computed during the unsplit integrator.
  \c d->Vs and \c Bs may be the same array or may be different.
  
  \b References
   - "A staggered mesh algorithm using high-order Godunov fluxes to 
      ensure solenoidal magnetic field in MHD simulations"\n
      Balsara \& Spicer, JCP (1999) 149, 270
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   Sep 09, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void CT_Update(const Data *d, Data_Arr Bs, double dt, Grid *grid)
/*!
 * Update staggered magnetic field using discrete version of 
 * Stoke's theorem.
 * Only \c d->Vs is updated, while \c Bs is the original array:\n
 * \c d->Vs = \c Bs + \c dt * \c R, where R = curl(E) is the electric field.
 *
 * \param [in,out] d     pointer to PLUTO Data structure. 
 *                       d->Vs will be updated.
 * \param [in]     Bs    the original array (not updated).
 * \param [in]     dt    step size
 * \param [in]     grid  pointer to Grid structure
 *
 *********************************************************************** */
{
  int  i, j, k, nv;
  int  ibeg, jbeg, kbeg;
  int  iend, jend, kend;
  double rhs_x, rhs_y, rhs_z;

  double *dx1 = grid->dx[IDIR], *x1p = grid->xr[IDIR], *x1m = grid->xl[IDIR];
  double *dx2 = grid->dx[JDIR], *x2p = grid->xr[JDIR], *x2m = grid->xl[JDIR];
  double *dx3 = grid->dx[KDIR], *x3p = grid->xr[KDIR], *x3m = grid->xl[KDIR];

  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];

  double Ax1, Ax2, Ax3;
  double dV1, dV2 , dV3;

  EMF *emf = d->emf;
  double ***Ex1 = emf->ex;  /* X1-comp. of emf at cell edges (i,j+1/2,k+1/2) */
  double ***Ex2 = emf->ey;  /* X2-comp. of emf at cell edges (i+1/2,j,k+1/2) */
  double ***Ex3 = emf->ez;  /* X3-comp. of emf at cell edges (i+1/2,j+1/2,k) */

/* ---- check div.B ---- */

#if CHECK_DIVB_CONDITION == YES
  if (g_intStage == 1) CT_CheckDivB (d->Vs, grid);
#endif
  
#if UPDATE_VECTOR_POTENTIAL == YES
  VectorPotentialUpdate (d, emf, NULL, grid);
#endif

/* --------------------------------------------------------
   2a. Update Bx1 at (i+1/2, j, k) faces
   -------------------------------------------------------- */

  for (k = emf->kbeg + KOFFSET; k <= emf->kend; k++){
  for (j = emf->jbeg + 1      ; j <= emf->jend; j++){
  for (i = emf->ibeg          ; i <= emf->iend; i++){
     
    Ax1 = grid->A[IDIR][k][j][i];

    #if GEOMETRY == CARTESIAN

    rhs_x = D_EXPAND(  0.0                                          , 
                     - dt/dx2[j]*(Ex3[k][j][i] - Ex3[k][j-1][i])    ,
                     + dt/dx3[k]*(Ex2[k][j][i] - Ex2[k-1][j][i]) ); 

    #elif GEOMETRY == CYLINDRICAL
    /* ------------------------------------------------
        Note that Ex3 = -E_\phi since cylindrical
        coordinates (R,z) are not right-handed
       ------------------------------------------------- */

    rhs_x = - dt/dx2[j]*(Ex3[k][j][i] - Ex3[k][j-1][i]); 

    #elif GEOMETRY == POLAR 
    Ax1 = fabs(x1p[i]);

    rhs_x = D_EXPAND(   0.0                                               ,
                     - dt/(Ax1*dx2[j])*(Ex3[k][j][i] - Ex3[k][j - 1][i])  , 
                     + dt/dx3[k]      *(Ex2[k][j][i] - Ex2[k - 1][j][i])); 

/*
    rhs_x = D_EXPAND(    0.0                                              ,
                     - dt*(Ex3[k][j][i] - Ex3[k][j-1][i])*dx3[k]          ,
                     + dt/Ax1*(Ex2[k][j][i] - Ex2[k-1][j][i])*x1p[i]*dx2[j] );
*/
    #elif GEOMETRY == SPHERICAL 

    double Ax2p = fabs(sin(x2p[j]));
    double Ax2m = fabs(sin(x2m[j]));

    dV2 = fabs(cos(x2m[j]) - cos(x2p[j]));
    rhs_x = D_EXPAND( 
            0.0                                                               ,
          - dt/(x1p[i]*dV2)*(Ax2p*Ex3[k][j][i] - Ax2m*Ex3[k][j - 1][i]) ,
          + dt*dx2[j]/(x1p[i]*dV2*dx3[k])*(Ex2[k][j][i] - Ex2[k - 1][j][i]));

    #endif    
      
    d->Vs[BX1s][k][j][i] = Bs[BX1s][k][j][i] + rhs_x;

  }}}

/* --------------------------------------------------------
   2b. Update Bx2 at (i, j+1/2, k) faces
   -------------------------------------------------------- */

  for (k = emf->kbeg + KOFFSET; k <= emf->kend; k++){
  for (j = emf->jbeg          ; j <= emf->jend; j++){
  for (i = emf->ibeg + 1      ; i <= emf->iend; i++){

    Ax2 = grid->A[JDIR][k][j][i];
     
    #if GEOMETRY == CARTESIAN

    rhs_y = D_EXPAND(  dt/dx1[i]*(Ex3[k][j][i] - Ex3[k][j][i-1])   ,
                                                                   ,   
                     - dt/dx3[k]*(Ex1[k][j][i] - Ex1[k-1][j][i]) ); 

    #elif GEOMETRY == CYLINDRICAL
    /* ------------------------------------------------
        Note that Ex3 = -E_\phi since cylindrical
        coordinates (R,z) are not right-handed
       ------------------------------------------------- */

double scrh = dt/(fabs(x1[i]*dx1[i]));
    rhs_y = scrh*(fabs(x1p[i])*Ex3[k][j][i] - fabs(x1m[i])*Ex3[k][j][i-1]);
/*
    rhs_y = dt/(r[i]*dx1[i])*(   fabs(rp[i])*Ex3[k][j][i]
                               - fabs(rp[i-1])*Ex3[k][j][i-1]);
*/

/*
double dV1 = fabs(grid->x[IDIR][i])*dx1[i];
double Ap   = grid->x[IDIR][i] + 0.5*dx1[i];
double Am   = grid->x[IDIR][i] - 0.5*dx1[i];
rhs_y =  dt/dV1*(fabs(Ap)*Ex3[k][j][i] - fabs(Am)*Ex3[k][j][i - 1]);
*/
    #elif GEOMETRY == POLAR 

     rhs_y =  D_EXPAND(  dt/dx1[i]*(Ex3[k][j][i] - Ex3[k][j][i - 1])    , 
                                                                        ,
                       - dt/dx3[k]*(Ex1[k][j][i] - Ex1[k - 1][j][i]));
 
    #elif GEOMETRY == SPHERICAL 

    Ax2 = fabs(sin(x2p[j]));
    rhs_y = D_EXPAND( 
            + dt/(x1[i]*dx1[i])*(x1p[i]*Ex3[k][j][i] - x1p[i-1]*Ex3[k][j][i - 1])    ,
                                                                               ,
            - dt/(x1[i]*Ax2*dx3[k])*(Ex1[k][j][i] - Ex1[k-1][j][i]));

    #endif    
      
    d->Vs[BX2s][k][j][i] = Bs[BX2s][k][j][i] + rhs_y;  

  }}}

/* --------------------------------------------------------
   2c. Update Bx3 at (i, j, k+1/2) faces
   -------------------------------------------------------- */

#if DIMENSIONS == 3
  for (k = emf->kbeg    ; k <= emf->kend; k++){
  for (j = emf->jbeg + 1; j <= emf->jend; j++){
  for (i = emf->ibeg + 1; i <= emf->iend; i++){
     
    Ax3 = grid->A[KDIR][k][j][i];
    #if GEOMETRY == CARTESIAN

    rhs_z = - dt/dx1[i]*(Ex2[k][j][i] - Ex2[k][j][i-1])
            + dt/dx2[j]*(Ex1[k][j][i] - Ex1[k][j-1][i]); 

    #elif GEOMETRY == POLAR 
    Ax1 = fabs(x1p[i]);
    dV1 = fabs(x1[i])*dx1[i];
    rhs_z = - dt/dV1*(fabs(x1p[i])*Ex2[k][j][i] - fabs(x1m[i])*Ex2[k][j][i-1])
            + dt/(x1[i]*dx2[j])*(Ex1[k][j][i] - Ex1[k][j - 1][i]);

/*
    rhs_z = - dt/Ax3*(fabs(x1p[i])*dx2[j]*Ex2[k][j][i]
                               fabs(x1p[i-1)*dx2[j]*Ex2[k][j][i - 1])
            + dt/Ax3*(Ex1[k][j][i]*dx1[i] - Ex1[k][j - 1][i]*dx1[j]);
*/
    #elif GEOMETRY == SPHERICAL 

    rhs_z = - dt/(x1[i]*dx1[i])*(x1p[i]*Ex2[k][j][i] - x1p[i-1]*Ex2[k][j][i-1])
            + dt/(x1[i]*dx2[j])*(Ex1[k][j][i] - Ex1[k][j-1][i]);
    #endif
      
    d->Vs[BX3s][k][j][i] = Bs[BX3s][k][j][i] + rhs_z;

  }}}
#endif

}

/* ********************************************************************* */
void CT_CheckDivB (double ***bf[], Grid *grid)
/*!
 * Check the divergence-free condition of magnetic field in the 
 * constrained transport formalism.
 * The solenoidal condition is discretized in a finite-volume sense:
 * \f[
 *      \int \nabla\cdot\vec{B}\, d^3x
 *    = \int \vec{B}\cdot d\hvec{S}
 *    = \sum_d \Big(\bar{B}_{d,+}S_{d,+} - \bar{B}_{d,-}S_{d,-} \Big)
 * \f]
 * where \f$S_{d,\pm}\f$ denotes the right (+) or left (-) interface surface
 * areas orthogonal to the \c d direction.
 * Thus in Cartesian coordinates one has
 * \f[
 *      \int \nabla\cdot\vec{B}\, d^3x
 *    =  \Delta y\Delta z \left(B_{x,+} - B_{x,-}\right)
 *     + \Delta x\Delta z \left(B_{y,+} - B_{y,-}\right)
 *     + \Delta x\Delta y \left(B_{z,+} - B_{z,-}\right)
 * \f]
 * while in spherical coordinates (\f$ dS_1 = r^2\,\sin\theta\,d\theta\,d\phi,
 * \, dS_2 = r\sin\theta\, dr\,d\phi,\, dS_3 = r\,dr\,d\theta\f$) we have
 * \f[
 *      \int \nabla\cdot\vec{B}\, d^3x
 *    =  \Delta \mu\Delta\phi
 *         \Big(r^2_+B_{r,+} - r^2_-B_{r,-}\Big)
 *     + \Delta\left(\frac{r^2}{2}\right) \Delta\phi
 *         \Big(  \sin\theta_+B_{\theta,+} - \sin\theta_-B_{\theta,-}\Big)
 *     + \Delta\left(\frac{r^2}{2}\right)\Delta\theta
 *         \Big(B_{\phi,+} - B_{\phi,-}\Big)
 * \f]
 * where \f$\mu = 1-\cos\theta\f$.
 * Notice also that \f$\Delta (r^2/2) = r\Delta r\f$ where \c r is the cell
 * center.
 * 
 * \param [in]  bf    an array of staggered magnetic field components
 * \param [in]  grid  a pointer to Grid structure
 *********************************************************************** */
{
  int i,j,k;
  double divB;
  D_EXPAND(double ***Bx1 = bf[0];  ,
           double ***Bx2 = bf[1];  ,
           double ***Bx3 = bf[2];)
  double dBx1, dBx2, dBx3, dbmax=0.0;
  double ***Ax1 = grid->A[IDIR];
  double ***Ax2 = grid->A[JDIR];
  double ***Ax3 = grid->A[KDIR];

/* ---------------------------------------
    Loop over computational domain
   --------------------------------------- */

/* 
double *dx = grid->dx[IDIR];
double *dy = grid->dx[JDIR];
double Ap, Am;
double *r  = grid->x[IDIR];
*/
  TOT_LOOP(k,j,i){
    D_EXPAND(dBx1 = Ax1[k][j][i]*Bx1[k][j][i] - Ax1[k][j][i-1]*Bx1[k][j][i-1];  ,
             dBx2 = Ax2[k][j][i]*Bx2[k][j][i] - Ax2[k][j-1][i]*Bx2[k][j-1][i];  ,
             dBx3 = Ax3[k][j][i]*Bx3[k][j][i] - Ax3[k-1][j][i]*Bx3[k-1][j][i];)
    divB = D_EXPAND(dBx1, + dBx2, + dBx3);

/*
Ap = fabs(r[i] + 0.5*dx[i]);
Am = fabs(r[i] - 0.5*dx[i]);
double divB_old;
double dBx1_old;
double dBx2_old;
dBx1_old = dy[j]*(Bx1[k][j][i]*Ap - Bx1[k][j][i-1]*Am);
dBx2_old = + fabs(r[i])*dx[i]*(Bx2[k][j][i] - Bx2[k][j-1][i]);
divB_old = dBx1_old + dBx2_old;

if (fabs(divB-divB_old) > 1.e-6){
  print ("! Divb not equal. (1) = %12.6e; (2) = %12.6e\n", divB, divB_old);
  print ("! dB[New] = %12.6e  %12.6e\n", dBx1, dBx2);
  print ("! dB[Old] = %12.6e  %12.6e\n", dBx1_old, dBx2_old);
  QUIT_PLUTO(1);
}
*/
    dbmax = MAX(dbmax,fabs(divB));
    if (fabs(divB) > 1.e-6) {
      print ("! CT_CheckDivB: div(B) = %12.6e, rank = %d, ijk = %d %d %d \n",
              divB, prank, i,j,k);
      D_EXPAND( print ("  Bx1: %12.6e  %12.6e\n",Bx1[k][j][i],Bx1[k][j][i-1]); ,
                print ("  Bx2: %12.6e  %12.6e\n",Bx2[k][j][i],Bx2[k][j-1][i]); ,
                print ("  Bx3: %12.6e  %12.6e\n",Bx3[k][j][i],Bx3[k-1][j][i]); )

      QUIT_PLUTO(1);
    }
  }
/*
#ifdef PARALLEL
  MPI_Allreduce (&dbmax, &divB, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  dbmax = divB;
#endif
  print ("       [divB] = %8.3e\n",dbmax);
*/
  
}
