/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Placing the particles on the grid and interpolating dervied quantites 
        at particle positions.
 
 \authors   A. Mignone (mignone@ph.unito.it)\n
            B. Vaidya (bvaidya@unito.it)\n
 
 \date    June 29, 2017
 */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Particles_Set(Data *d, Grid *grid)
/*!
 * Initialize particles on the grid and distribute them among
 * processors.
 *
 *  \param [in]    d     Pointer to the PLUTO data structure.
 *  \param [in]    grid  Pointer to the PLUTO grid structure.
 *
 *********************************************************************** */
{ 
#if PARTICLES_TYPE == LAGRANGIAN
  char particles_type[32] = "LAGRANGIAN";
#elif PARTICLES_TYPE == COSMIC_RAYS
  char particles_type[32] = "COSMIC_RAYS";
#elif PARTICLES_TYPE == DUST
  char particles_type[32] = "DUST";  
#endif
  Particle p1;

  print ("> Particles_Set():\n");

/* ------------------------------------------------------------
   0. Do some checks before continuing...
   ------------------------------------------------------------ */

#if PARTICLES_TYPE == COSMIC_RAYS || PARTICLES_TYPE == DUST
  #if COMPONENTS != 3
  print ("! Particles_Set(): number of COMPONENTS must be 3\n");
  QUIT_PLUTO(1);
  #endif
#endif

/* ----------------------------------------------------------
   1. Define MPI data type (not strictly necessary) 
   ---------------------------------------------------------- */ 

#ifdef PARALLEL
  Particles_StructDatatype();
#endif  
  p_nrestart = 0;
  
/* ----------------------------------------------------------
   2. Initialize particles on grid
   ---------------------------------------------------------- */ 

  p_nparticles   = 0;
  Particles_Init(d, grid);

/* ----------------------------------------------------------
   3. Do some printing
   ---------------------------------------------------------- */ 

  print ("  particles type: %s\n",particles_type); 
  print ("  number of particles created [local]: %d\n\n",p_nparticles);  

#if PARTICLES_TYPE == LAGRANGIAN
  particleNode *CurNode = d->PHead;
  Particle *pl;
  #if PARTICLES_LP_SPECTRA == YES
  PARTICLES_LOOP(CurNode, d->PHead){
    pl = &(CurNode->p);
    Particles_LP_FixValue(pl, d, grid);
  }
  #endif
#endif

#if PARTICLES_USE_ARRAY == YES                                    
    Particles_ListToArray(d);                                     
 #endif

#ifdef PARALLEL 
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}
