#include "pluto.h"

/* ------------------------------------------------------------------
           set boundary conditions on electric field 
   ------------------------------------------------------------------ */

/* ******************************************************************* */
void EMF_BOUNDARY (EMF *emf, Grid *grid)
/*
 *
 * PURPOSE:
 *
 *  Set boundary conditions on electric field. 
 *  Quantities are evaluated at zone face, NOT edge.
 *  
 *  The electric field components are given by 
 *  the fluxes of the induction equation:
 *
 *   +----+     
 *   |  0 |
 *   |-Ezi|  at x - faces
 *   | Eyi|
 *   +----+
 *
 *   +----+     
 *   | Ezj|
 *   | 0  |  at y - faces
 *   |-Exj|
 *   +----+
 *
 *   +----+     
 *   |-Eyk|
 *   | Exk|  at z - faces
 *   | 0  |
 *   +----+
 *
 *  Only the components required to update the normal megnetic
 *  field are required.
 *  For example, at the lower x-boundary, the necessary 
 *  electric field components to update bx are Ez and Ey (no Ex).
 *  Thus one needs to specify 
 *
 *     Ez at y-faces   (ezj)
 *     Ey at z-faces   (eyk)
 *  
 *
 * LAST MODIFIED: 
 *
 *   April 27 2007 by A. Mignone (email:mignone@to.astro.it)
 *
 ******************************************************************* */
{
  int side[6] = {X1_BEG, X1_END, X2_BEG, X2_END, X3_BEG, X3_END};
  int type[6], is;
  int par_dim[3] = {0, 0, 0};
  double ***e1, ***e2, vsign;

  D_EXPAND(par_dim[0] = grid[IDIR].nproc > 1;  ,
           par_dim[1] = grid[JDIR].nproc > 1;  ,
           par_dim[2] = grid[KDIR].nproc > 1;)

  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);

   AL_Exchange_dim ((char *)emf->ezj[0][0], par_dim, SZ);
   AL_Exchange_dim ((char *)emf->ezi[0][0], par_dim, SZ);
   #if DIMENSIONS == 3
    AL_Exchange_dim ((char *)emf->exj[0][0], par_dim, SZ);
    AL_Exchange_dim ((char *)emf->exk[0][0], par_dim, SZ);
    AL_Exchange_dim ((char *)emf->eyi[0][0], par_dim, SZ);
    AL_Exchange_dim ((char *)emf->eyk[0][0], par_dim, SZ);
   #endif
/*
   AL_Exchange (emf->ezj[0][0], SZ);
   AL_Exchange (emf->ezi[0][0], SZ);
   #if DIMENSIONS == 3
    AL_Exchange (emf->exj[0][0], SZ);
    AL_Exchange (emf->exk[0][0], SZ);
    AL_Exchange (emf->eyi[0][0], SZ);
    AL_Exchange (emf->eyk[0][0], SZ);
   #endif
*/
   MPI_Barrier (MPI_COMM_WORLD);

  #endif

return;


  type[0] = grid[IDIR].lbound;
  type[1] = grid[IDIR].rbound;

  type[2] = grid[JDIR].lbound;
  type[3] = grid[JDIR].rbound;

  type[4] = grid[KDIR].lbound;
  type[5] = grid[KDIR].rbound;

/* --------------------------------------------------------
                  Loop on directions       
   -------------------------------------------------------- */

  for (is = 0; is < 2*DIMENSIONS; is++){

    if (type[is] == 0) continue;  /* no physical boundary: skip */

  /* -- set pointers to the emf components 
        that need boundaries               -- */

    if (side[is] == X1_BEG || side[is] == X1_END){
      D_EXPAND(   , e1 = emf->ezj;, e2 = emf->eyk;)
    } else if (side[is] == X2_BEG || side[is] == X2_END){
      D_EXPAND(   , e1 = emf->ezi;, e2 = emf->exk;)
    } else if (side[is] == X3_BEG || side[is] == X3_END){
      D_EXPAND(   , e1 = emf->exj;, e2 = emf->eyi;)
    }

    if (type[is] == OUTFLOW){
/* --> replace in OutflowBound(..)
      D_EXPAND(                                          ,
        OUTFLOW_BOUND(e1, grid + (is/2), -1, side[is]);  ,
        OUTFLOW_BOUND(e2, grid + (is/2), -1, side[is]);)
*/
    }else if (   (type[is] == REFLECTIVE) 
              || (type[is] == AXISYMMETRIC)
              || (type[is] == EQTSYMMETRIC) ){

      vsign = -1.0;
      if (type[is] == EQTSYMMETRIC) vsign = 1.0;
/* -- > replace in ReflectiveBound(...)
      D_EXPAND(                                 ,
        REFLECTIVE_BOUND(e1, -1, side[is], vsign);  ,
        REFLECTIVE_BOUND(e2, -1, side[is], vsign);)
*/    
    }else if (type[is] == USERDEF) {

      /* -------------------------------------------------
          for userdef boundaries, ghost EMF must be 
          computed from the time marching algorithm by 
          extending integration in the corresponding 
          (transverse) boundary zones.
         ------------------------------------------------- */
          
      EMF_USERDEF_BOUNDARY (emf, side[is], FACE_EMF, grid);  

    }else if (type[is] == PERIODIC){


      if (!par_dim[is/2]) {

/* -- replace in PeriodicBound
        D_EXPAND(                        ,
          PERIODIC_BOUND(e1, -1, side[is]);  ,
          PERIODIC_BOUND(e2, -1, side[is]);)
*/
      } 
/*
      #ifndef PARALLEL   
       D_EXPAND(                        ,
         PERIODIC_BOUND(e1, side[is]);  ,
         PERIODIC_BOUND(e2, side[is]);)
      #endif
*/
    } 
  }
}

