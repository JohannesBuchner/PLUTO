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
         
  \date   Sept 12, 2012
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
 * Compute the mean orbital velocity as a 2D array by averaging
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

  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

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
      for (i = 0; i < DIMENSIONS; i++) coords[i] = grid[i].rank_coord;
      coords[SDIR] += 1;
      MPI_Cart_rank(cartcomm, coords, &dst);

      for (i = 0; i < DIMENSIONS; i++) coords[i] = grid[i].rank_coord;
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
       if (grid[SDIR].nproc > 1){
         #ifdef PARALLEL
          w_sum = w;
          for (count=1; count < grid[SDIR].nproc; count++ ){
            MPI_Sendrecv(&w, 1, MPI_DOUBLE, dst, 0,  &w_recv, 1,
                                MPI_DOUBLE, src, 0, cartcomm, &status);
            w      = w_recv;
            w_sum += w;
          }
          w = w_sum/(double)(grid[SDIR].np_int_glob);
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
      print1 ("! FARGO not supported in this geometry\n");
      QUIT_PLUTO(1);
     #endif
   }
  #endif  /* FARGO_AVERAGE_VELOCITY */
  
  first_call = 0;
}

static int vphi_is_total_velocity = YES; /* Used to make sure either operation
                                            is not repeated twice in a row */
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

  #if GEOMETRY == CARTESIAN
   vphi = d->Vc[VX2];
  #else
   vphi = d->Vc[iVPHI];
  #endif

/* ----------------------------------------------------------------------
                      subtract velocity
   ---------------------------------------------------------------------- */

  TOT_LOOP(k,j,i){
    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
     vphi[k][j][i] -= wA[k][i];
    #elif GEOMETRY == SPHERICAL
     vphi[k][j][i] -= wA[j][i];
    #else
     print1 ("! FARGO not supported in this geometry\n");
     QUIT_PLUTO(1);
    #endif
  }

  vphi_is_total_velocity = NO;
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
     print1 ("! FARGO not supported in this geometry\n");
     QUIT_PLUTO(1);
    #endif 
  }

  vphi_is_total_velocity = YES;
}

/* ********************************************************************* */
int FARGO_HasTotalVelocity ()
/*!
 *  Return 1 if the velocity array contains the total velocity.
 *  Return 0 otherwise.
 *********************************************************************** */
{
  return vphi_is_total_velocity;
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
     print1 ("! FARGO_GetVelocity: incorrect direction\n");
     QUIT_PLUTO(1);
    #endif
    
  }else if (where == X3FACE){

    #if GEOMETRY == CARTESIAN || GEOMETRY == POLAR
     return wAx3;
    #else
     print1 ("! FARGO_GetVelocity: incorrect direction\n");
     QUIT_PLUTO(1);
    #endif

  }

*/
}
