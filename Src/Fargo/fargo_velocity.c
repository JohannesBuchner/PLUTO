/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Functions for computing/retrieving the mean aziumthal velocity.

  Collects various functions for computing, adding, subtracting and
  retrieving the azimuthally-averaged velocity.
  Depending on the macro ::FARGO_AVERAGE_VELOCITY (YES/NO) this is
  done either numerically or analytically.

  \authors A. Mignone (mignone@ph.unito.it)\n
           G. Muscianisi (g.muscianisi@cineca.it)
         
  \date    June 22, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/*! Defines a 2D array containing the azimuthally-averaged velocity.
    When the average is done along the X2 direction (Cartesian and polar
    geometries) the array should have size (NX3_TOT, NX1_TOT).
    When the average is performed along the X3 direction (spherical geometry)
    the array has dimensions (NX2_TOT, NX1_TOT)                           */
static double **wA;  
/*static double **wAx1, **wAx2, **wAx3;*/  /* -- average orbital speed at
                                            x1, x2 and x3 faces -- */
                                            
/* ********************************************************************* */
void FARGO_ComputeVelocity(const Data *d, Grid *grid)
/*!
 * Compute the average orbital velocity as a 2D array by averaging
 * for each pair (x,z)/(r,z)/(r,theta)  (in Cartesian/polar/spherical)
 * along the orbital coordinate y/phi/phi.
 *
 *********************************************************************** */
#if GEOMETRY != SPHERICAL
 #define s j   /* -- orbital direction index -- */
#else
 #define s k
#endif
{
  int i,j,k;
  static int first_call = 1;
  double *x1, *x2, *x3, th, w;
  double ***vphi;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
   if (wA == NULL){
     wA   = ARRAY_2D(NX3_TOT, NX1_TOT, double);
/*
     wAx1 = ARRAY_2D(NX3_TOT, NX1_TOT, double);  
     wAx3 = ARRAY_2D(NX3_TOT, NX1_TOT, double);
*/
   }
  #else
   if (wA == NULL){
     wA   = ARRAY_2D(NX2_TOT, NX1_TOT, double);
/*
     wAx1 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
     wAx2 = ARRAY_2D(NX2_TOT, NX1_TOT, double);
*/
   }
  #endif

#if GEOMETRY == CARTESIAN
  vphi = d->Vc[VX2];
#else
  vphi = d->Vc[iVPHI];
#endif

#if FARGO_AVERAGE_VELOCITY == YES

/* ---------------------------------------------------------------------
      average velocity along the transport direction every 10 steps
   --------------------------------------------------------------------- */

  if (g_stepNumber%FARGO_NSTEP_AVERAGE == 0 || first_call){
  #ifdef PARALLEL
    double w_recv=0.0, w_sum=0.0;
    int count, coords[3], src, dst;
    MPI_Comm cartcomm;
    MPI_Status status;
  #endif
     
  /* -- fill ghost zones -- */

    Boundary(d, ALL_DIR, grid); 

  /* ---- get ranks of the upper and lower processsors ---- */

  #ifdef PARALLEL
    AL_Get_cart_comm(SZ, &cartcomm);
    for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
    coords[SDIR] += 1;
    MPI_Cart_rank(cartcomm, coords, &dst);

    for (i = 0; i < DIMENSIONS; i++) coords[i] = grid->rank_coord[i];
    coords[SDIR] -= 1;
    MPI_Cart_rank(cartcomm, coords, &src);
  #endif

  /* ---- compute average orbital speed ---- */

  #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    KTOT_LOOP(k) ITOT_LOOP(i) {
  #else
    JTOT_LOOP(j) ITOT_LOOP(i) {
  #endif
      w = 0.0;
      SDOM_LOOP(s) w += vphi[k][j][i];
      if (grid->nproc[SDIR] > 1){
      #ifdef PARALLEL
        w_sum = w;
        for (count=1; count < grid->nproc[SDIR]; count++ ){
          MPI_Sendrecv(&w, 1, MPI_DOUBLE, dst, 0,  &w_recv, 1,
                              MPI_DOUBLE, src, 0, cartcomm, &status);
          w      = w_recv;
          w_sum += w;
        }
        w = w_sum/(double)(grid->np_int_glob[SDIR]);
      #endif
      }else{
         w = w/(double)NS;
      }
      #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
        wA[k][i] = w;
      #elif GEOMETRY == SPHERICAL
        wA[j][i] = w;
      #endif
    }
  }

  /* -- compute interface velocity -- */
/*
   #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR

    KTOT_LOOP(k) for (i = 0; i < NX1_TOT-1; i++){
      wAx1[k][i] = 0.5*(wA[k][i] + wA[k][i+1]);
    }
    for (k = 0; k < NX3_TOT-1; k++) ITOT_LOOP(i){
      wAx3[k][i] = 0.5*(wA[k][i] + wA[k+1][i]);
    }

   #elif GEOMETRY == SPHERICAL

    JTOT_LOOP(j) for (i = 0; i < NX1_TOT-1; i++){
      wAx1[j][i] = 0.5*(wA[j][i] + wA[j][i+1]);
    }
    for (j = 0; j < NX2_TOT-1; j++) ITOT_LOOP(i){
      wAx2[j][i] = 0.5*(wA[j][i] + wA[j+1][i]);
    }
   #endif
*/

#else

/* ---------------------------------------------------------------------
              define velocity analytically
   --------------------------------------------------------------------- */

  if (first_call){
  #if GEOMETRY == CARTESIAN 
    KTOT_LOOP(k) ITOT_LOOP(i) wA[k][i] = FARGO_SetVelocity(x1[i], x3[k]);
  #elif GEOMETRY == POLAR
    KTOT_LOOP(k) ITOT_LOOP(i) wA[k][i] = x1[i]*FARGO_SetVelocity(x1[i], x3[k]);
  #elif GEOMETRY == SPHERICAL
    JTOT_LOOP(j) ITOT_LOOP(i) {
      th = x2[j];
      wA[j][i] = x1[i]*sin(th)*FARGO_SetVelocity(x1[i], x2[j]);
    }
  #else
    print ("! FARGO not supported in this geometry\n");
    QUIT_PLUTO(1);
   #endif
  }
#endif  /* FARGO_AVERAGE_VELOCITY */
  
  first_call = 0;
}

static int tot_vphi=1; /* A flag that tracks whether the input velocity field
                          contains the residual (tot_vphi=0) or the total
                          velocity (tot_vphi=1).
                          On first call it is assumed that tot_vphi=1.
/* ********************************************************************* */
void FARGO_SubtractVelocity(const Data *d, Grid *grid)
/*!
 * Subtract the mean background contribution from the total 
 * velocity. In other words, compute the residual velocity.
 *
 * \param [in,out] d  pointer to PLUTO Data structure
 * \param [in]     grid pointer to array of Grid structures
 *
 * \return This function has no return value.
 *********************************************************************** */
{
  int i,j,k;
  double ***vphi;

/* -------------------------------------------------------
   0. Check that input velocity does not contain residual
   ------------------------------------------------------- */

  if (tot_vphi == 0) {
    print ("! FARGO_SubtractVelocity(): input velocity already contains residual\n");
    QUIT_PLUTO(1);
  }

/* -----------------------------------------------------
   1. Assign pointer shortcut
   ----------------------------------------------------- */

#if GEOMETRY == CARTESIAN
  vphi = d->Vc[VX2];
#else
  vphi = d->Vc[iVPHI];
#endif

/* ------------------------------------------------------
   2. Subtract velocity
   ------------------------------------------------------ */

  TOT_LOOP(k,j,i){
    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    vphi[k][j][i] -= wA[k][i];
    #elif GEOMETRY == SPHERICAL
    vphi[k][j][i] -= wA[j][i];
    #else
    print ("! FARGO not supported in this geometry\n");
    QUIT_PLUTO(1);
    #endif
  }


#ifdef PARTICLES
  int i1,j1,k1;
  particleNode *curNode;
  Particle *p;
  double v;
  static double ***wp;
  if (wp == NULL){
    wp = ArrayBox (-1, 1, -1, 1, -1, 1);
  }

  PARTICLES_LOOP(curNode, d->PHead){
    p = &(curNode->p);

  /* -- Interpolate fargo shift velocity at particle position -- */

    Particles_GetWeights(p, p->cell, wp, grid);  
    i = p->cell[IDIR];
    j = p->cell[JDIR];
    k = p->cell[KDIR];

    v = 0.0;
    for (k1 = -KOFFSET; k1 <= KOFFSET; k1++){
    for (j1 = -JOFFSET; j1 <= JOFFSET; j1++){
    for (i1 = -IOFFSET; i1 <= IOFFSET; i1++){
      v += wp[k1][j1][i1]*wA[k + k1][i + i1];
    }}}


// v = -SB_Q*SB_OMEGA*p->coord[IDIR];
//print ("vy' = %12.6e\n",p->speed[JDIR]);
    p->speed[JDIR] -= v;
//print ("vy' = %12.6e\n",p->speed[JDIR]);
  }
#endif

  tot_vphi = 0;
}
/* ************************************************************** */
void FARGO_AddVelocity(const Data *d, Grid *grid)
/*!
 * Add the mean backgroun contribution to the residual
*  velocity in order to obtain the total velocity.
 *
 * \param [in,out] d  pointer to PLUTO Data structure
 * \param [in]     grid pointer to array of Grid structures
 *
 * \return This function has no return value.
 *********************************************************************** */
{
  int i,j,k;
  double ***vphi;

/* -----------------------------------------------------------
    Check that input velocity does not contain total velocity
   ----------------------------------------------------------- */

  if (tot_vphi == 1) {
    print ("! FARGO_AddVelocity(): input velocity already contains total velocity\n");
    QUIT_PLUTO(1);
  }

#if GEOMETRY == CARTESIAN
  vphi = d->Vc[VX2];
#else
  vphi = d->Vc[iVPHI];
#endif

  TOT_LOOP(k,j,i){
    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
    vphi[k][j][i] += wA[k][i];
    #elif GEOMETRY == SPHERICAL
    vphi[k][j][i] += wA[j][i];
    #else
    print ("! FARGO not supported in this geometry\n");
    QUIT_PLUTO(1);
    #endif 
  }

#ifdef PARTICLES
  int i1,j1,k1;
  particleNode *curNode;
  Particle *p;
  double v;
  static double ***wp;
  if (wp == NULL){
    wp = ArrayBox (-1, 1, -1, 1, -1, 1);
  }

  PARTICLES_LOOP(curNode, d->PHead){
    p = &(curNode->p);

  /* -- Interpolate fargo shift velocity at particle position -- */

    Particles_GetWeights(p, p->cell, wp, grid);  
    i = p->cell[IDIR];
    j = p->cell[JDIR];
    k = p->cell[KDIR];

    v = 0.0;
    for (k1 = -KOFFSET; k1 <= KOFFSET; k1++){
    for (j1 = -JOFFSET; j1 <= JOFFSET; j1++){
    for (i1 = -IOFFSET; i1 <= IOFFSET; i1++){
      v += wp[k1][j1][i1]*wA[k + k1][i + i1];
    }}}
// v = -SB_Q*SB_OMEGA*p->coord[IDIR];
    p->speed[JDIR] += v;
  }
#endif

  tot_vphi = 1;
}

/* ********************************************************************* */
double **FARGO_GetVelocity(void)
/*!
 * Return a pointer to the background orbital velocity ::wA.
 *
 *********************************************************************** */
{
  return wA;
/*
  int i, j, k;

  if (where == CENTER){

    return wA;

  }else if (where == X1FACE){

     return wAx1;

  }else if (where == X2FACE){

    #if GEOMETRY == SPHERICAL
     return wAx2;
    #else
     print ("! FARGO_GetVelocity: incorrect direction\n");
     QUIT_PLUTO(1);
    #endif
    
  }else if (where == X3FACE){

    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
     return wAx3;
    #else
     print ("! FARGO_GetVelocity: incorrect direction\n");
     QUIT_PLUTO(1);
    #endif

  }

*/
}
/* ********************************************************************* */
int FARGO_TotalVelocityIsSet(void)
/*!
 * Return 1 if velocity is the azimuthal direction contains
 * the total velocity. Return 0 is velocity contains the residual.
 *
 *********************************************************************** */
{
  return tot_vphi;
}

