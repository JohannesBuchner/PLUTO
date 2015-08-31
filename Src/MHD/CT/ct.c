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
  \date   April 1, 2014
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
 * \param [in]     grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int  i, j, k, nv;
  int  ibeg, jbeg, kbeg;
  int  iend, jend, kend;
  double rhs_x, rhs_y, rhs_z;
  double *dx, *dy, *dz, *A1, *A2, *dV1, *dV2;
  double *r, *rp, *th, r_2;
  double ***ex, ***ey, ***ez;
  EMF *emf;

/* ---- check div.B ---- */

  #if CHECK_DIVB_CONDITION == YES
   if (g_intStage == 1) CT_CheckDivB (d->Vs, grid);
  #endif

  emf = CT_GetEMF (d, grid);

  #if UPDATE_VECTOR_POTENTIAL == YES
   VectorPotentialUpdate (d, emf, NULL, grid);
  #endif

  ex = emf->ex;
  ey = emf->ey;
  ez = emf->ez;

/* ------------------------------------------------------------
    Correct EMF if the ShearingBox module is enabled.
    For CTU, this step is carried only during the final stage.
   ------------------------------------------------------------ */
   
  #ifdef SHEARINGBOX
   #ifdef CTU 
    if (g_intStage == 2) SB_CorrectEMF (emf, Bs, grid);
   #else
    SB_CorrectEMF (emf, d->Vs, grid);
   #endif
  #endif

  dx = grid[IDIR].dx; 
  dy = grid[JDIR].dx;
  dz = grid[KDIR].dx;

  r   = grid[IDIR].x;
  rp  = grid[IDIR].xr;
  A1  = grid[IDIR].A;
  dV1 = grid[IDIR].dV;

  th  = grid[JDIR].x;
  A2  = grid[JDIR].A;
  dV2 = grid[JDIR].dV;

/* --------------------------------------------------------
          Update Bx1 at (i+1/2, j, k)  faces
   -------------------------------------------------------- */

  for (k = emf->kbeg + KOFFSET; k <= emf->kend; k++){
  for (j = emf->jbeg + 1      ; j <= emf->jend; j++){
  for (i = emf->ibeg          ; i <= emf->iend; i++){
     
    #if GEOMETRY == CARTESIAN

     rhs_x = D_EXPAND(  0.0                                       , 
                      - dt/dy[j]*(ez[k][j][i] - ez[k][j - 1][i])  ,
                      + dt/dz[k]*(ey[k][j][i] - ey[k - 1][j][i]) ); 

    #elif GEOMETRY == CYLINDRICAL

     rhs_x = - dt/dy[j]*(ez[k][j][i] - ez[k][j - 1][i]);

    #elif GEOMETRY == POLAR 

     rhs_x = D_EXPAND(   0.0                                               ,
                       - dt/(A1[i]*dy[j])*(ez[k][j][i] - ez[k][j - 1][i])  , 
                       + dt/dz[k]        *(ey[k][j][i] - ey[k - 1][j][i])); 

    #elif GEOMETRY == SPHERICAL 

     rhs_x = D_EXPAND( 
             0.0                                                               ,
           - dt/(rp[i]*dV2[j])*(A2[j]*ez[k][j][i] - A2[j - 1]*ez[k][j - 1][i]) ,
           + dt*dy[j]/(rp[i]*dV2[j]*dz[k])*(ey[k][j][i] - ey[k - 1][j][i]));

    #endif    
      
    d->Vs[BX1s][k][j][i] = Bs[BX1s][k][j][i] + rhs_x;

  }}}

/* --------------------------------------------------------
          Update Bx2 at (i, j+1/2, k)  faces
   -------------------------------------------------------- */

  for (k = emf->kbeg + KOFFSET; k <= emf->kend; k++){
  for (j = emf->jbeg          ; j <= emf->jend; j++){
  for (i = emf->ibeg + 1      ; i <= emf->iend; i++){
     
    #if GEOMETRY == CARTESIAN

     rhs_y = D_EXPAND(  dt/dx[i]*(ez[k][j][i] - ez[k][j][i - 1])   ,  
                                                                   ,   
                      - dt/dz[k]*(ex[k][j][i] - ex[k - 1][j][i]) ); 

    #elif GEOMETRY == CYLINDRICAL

     rhs_y =   dt/dV1[i]*(A1[i]*ez[k][j][i] - A1[i - 1]*ez[k][j][i - 1]);

    #elif GEOMETRY == POLAR 

     rhs_y =  D_EXPAND(  dt/dx[i]*(ez[k][j][i] - ez[k][j][i - 1])    , 
                                                                     ,
                       - dt/dz[k]*(ex[k][j][i] - ex[k - 1][j][i]));
 
    #elif GEOMETRY == SPHERICAL 

     rhs_y = D_EXPAND( 
            + dt/(r[i]*dx[i])*(rp[i]*ez[k][j][i] - rp[i - 1]*ez[k][j][i - 1])    ,
                                                                               ,            - dt/(r[i]*A2[j]*dz[k])*(ex[k][j][i] - ex[k - 1][j][i]));

    #endif    
      
    d->Vs[BX2s][k][j][i] = Bs[BX2s][k][j][i] + rhs_y;  

  }}}

/* --------------------------------------------------------
          Update Bx3 at (i, j, k+1/2)  faces
   -------------------------------------------------------- */

  #if DIMENSIONS == 3
  for (k = emf->kbeg    ; k <= emf->kend; k++){
  for (j = emf->jbeg + 1; j <= emf->jend; j++){
  for (i = emf->ibeg + 1; i <= emf->iend; i++){
     
     #if GEOMETRY == CARTESIAN

      rhs_z = - dt/dx[i]*(ey[k][j][i] - ey[k][j][i - 1])
              + dt/dy[j]*(ex[k][j][i] - ex[k][j - 1][i]); 

     #elif GEOMETRY == POLAR 

      rhs_z = - dt/dV1[i]*(A1[i]*ey[k][j][i] - A1[i - 1]*ey[k][j][i - 1])
              + dt/(r[i]*dy[j])*(ex[k][j][i] - ex[k][j - 1][i]);

     #elif GEOMETRY == SPHERICAL 

      rhs_z = - dt/(r[i]*dx[i])*(rp[i]*ey[k][j][i] - rp[i - 1]*ey[k][j][i - 1])
              + dt/(r[i]*dy[j])*(ex[k][j][i] - ex[k][j - 1][i]);
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
 *
 *********************************************************************** */
{
  int i,j,k;
  double divB;
  double ***bx, ***by, ***bz;
  double *dx, *dy, *dz;
  double *Ar, *s, *dmu, *r, *rp, *rm;
  double dbmax=0.0;

  D_EXPAND( bx = bf[0];  ,
            by = bf[1];  ,
            bz = bf[2]; )
    
  rp  = grid[IDIR].xr;
  rm  = grid[IDIR].xl;

  Ar  = grid[IDIR].A;
  s   = grid[JDIR].A;
  dmu = grid[JDIR].dV;

  r    = grid[IDIR].x;
  dx   = grid[IDIR].dx;
  dy   = grid[JDIR].dx;
  dz   = grid[KDIR].dx;
 

/*
for (i = IBEG-1; i >= -1; i--){
for (j=0; j < NX2_TOT; j++){
  divB = bx[0][j][i]-bx[0][j][i+NX];
  if (fabs(divB) > 1.e-12){ 
    printf ("!! periodicity not respected (bx), zone = %d, %d\n",i,j);
    exit(1);
  }
}}
for (i = IBEG-1; i >= 0; i--){
for (j = -1; j < NX2_TOT; j++){
  divB = by[0][j][i]-by[0][j][i+NX];
  if (fabs(divB) > 1.e-12){ 
    printf ("!! periodicity not respected (by)\n");
    exit(1);
  }
}}
*/
  for (k = 0; k < NX3_TOT; k++){

    #if DIMENSIONS < 3
     dz[k] = 1.0;
    #endif

    for (j = 0; j < NX2_TOT; j++){
    for (i = 0; i < NX1_TOT; i++){

    #if GEOMETRY == CARTESIAN

     divB = D_EXPAND(   (bx[k][j][i] - bx[k][j][i-1])*dy[j]*dz[k]  ,
                      + (by[k][j][i] - by[k][j-1][i])*dx[i]*dz[k]  , 
                      + (bz[k][j][i] - bz[k-1][j][i])*dx[i]*dy[j] );

    #elif GEOMETRY == CYLINDRICAL

     divB =        dy[j]*(bx[k][j][i]*Ar[i] - bx[k][j][i-1]*Ar[i-1])
            + fabs(r[i])*dx[i]*(by[k][j][i] - by[k][j-1][i]);

/*
     divB =        dy[j]*(bx[k][j][i]*fabs(rp[i]) - bx[k][j][i-1]*fabs(rm[i]))
            + fabs(r[i])*dx[i]*(by[k][j][i] - by[k][j-1][i]);
*/
    #elif GEOMETRY == POLAR

     divB = D_EXPAND(
              dz[k]*dy[j]*(bx[k][j][i]*Ar[i] - bx[k][j][i-1]*Ar[i-1]) ,
            + dx[i]*dz[k]*(by[k][j][i] - by[k][j-1][i])               ,
            + r[i]*dx[i]*dy[j]*(bz[k][j][i] - bz[k-1][j][i])
           );

    #elif GEOMETRY == SPHERICAL
         
     divB = D_EXPAND(  
               dmu[j] *dz[k]*(bx[k][j][i]*Ar[i]  - bx[k][j][i-1]*Ar[i-1])   ,
             + r[i]*dx[i]*dz[k]*(by[k][j][i]*s[j] - by[k][j-1][i]*s[j-1])  ,
             + r[i]*dx[i]*dy[j]*(bz[k][j][i]        - bz[k-1][j][i])
            );

    #endif

     dbmax = MAX(dbmax,fabs(divB));
     if (fabs(divB) > 1.e-6) {
       print ("! CHECK_DIVB: div(B) = %12.6e != 0, rank = %d, ijk = %d %d %d \n",
               divB, prank, i,j,k);
       print ("! rank =%d, g_intStage = %d, Beg, End = [%d, %d]  [%d, %d]  [%d, %d]\n",
               prank, g_intStage, IBEG, IEND,
               JBEG, JEND, KBEG, KEND);
       D_EXPAND( print ("b1: %12.6e  %12.6e\n",bx[k][j][i],bx[k][j][i-1]); ,
                 print ("b2: %12.6e  %12.6e\n",by[k][j][i],by[k][j-1][i]); ,
                 print ("b3: %12.6e  %12.6e\n",bz[k][j][i],bz[k-1][j][i]); )

       QUIT_PLUTO(1);
     }
    }}
  }

}
