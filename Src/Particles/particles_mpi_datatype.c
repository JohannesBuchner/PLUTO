/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Define MPI data-type structure for MPI communications

  Define the MPI Datatype of the Particle structure so that it can be 
  passed among processors. This has to be done for each PARTICLES_TYPE
  depending on the quantities present in the structure.
 
  \authors   A. Mignone (mignone@ph.unito.it)\n
             B. Vaidya (bvaidya@unito.it)\n
  
  \date      March 31, 2018
  
  \b References
  http://mpi.deino.net/mpi_functions/MPI_Type_create_struct.html
  https://stackoverflow.com/questions/33618937/trouble-understanding-mpi-type-create-struct
  https://stackoverflow.com/questions/9864510/struct-serialization-in-c-and-transfer-over-mpi
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include <stddef.h> /* Required for using offsetof() function */

#ifdef PARALLEL
/* ************************************************************************** */
void Particles_StructDatatype()
/*!
 *
 **************************************************************************** */
{
  int i = 0;
  Particle pl[1];

#if (PARTICLES_TYPE == COSMIC_RAYS) || (PARTICLES_TYPE == DUST)
  #define N_ELEMENTS  9
  MPI_Datatype types[N_ELEMENTS] = {MPI_DOUBLE,    /* coord[3]      */
                                    MPI_DOUBLE,    /* speed[3]      */
                                    MPI_DOUBLE,    /* coord_old[3]  */
                                    MPI_DOUBLE,    /* speed_old[3] */
                                    MPI_DOUBLE,    /* rho          */
                                    MPI_FLOAT,     /* tinj         */ 
                                    MPI_FLOAT,     /* color        */ 
                                    MPI_INT,       /* cell[3]      */
                                    MPI_UINT32_T}; /* id           */

  int blocklengths[N_ELEMENTS]  = {3,3,3,3,1,1,1,3,1};
  MPI_Aint offsets[N_ELEMENTS];

  i=0;
  offsets[i++] = offsetof(Particle, coord);
  offsets[i++] = offsetof(Particle, speed);
  offsets[i++] = offsetof(Particle, coord_old);
  offsets[i++] = offsetof(Particle, speed_old);
  offsets[i++] = offsetof(Particle, rho);
  offsets[i++] = offsetof(Particle, tinj);
  offsets[i++] = offsetof(Particle, color);
  offsets[i++] = offsetof(Particle, cell);
  offsets[i++] = offsetof(Particle, id);
  
  MPI_Type_create_struct(N_ELEMENTS, blocklengths, offsets, types, &MPI_PARTICLE);
  MPI_Type_commit(&MPI_PARTICLE);
                                   
#elif (PARTICLES_TYPE == LAGRANGIAN) && (PARTICLES_LP_SPECTRA == NO)
  #define N_ELEMENTS  9

  MPI_Datatype types[N_ELEMENTS] = {MPI_DOUBLE,     /* coord[3]      */
                                    MPI_DOUBLE,     /* speed[3]     */
                                    MPI_DOUBLE,     /* coord_old[3]  */
                                    MPI_DOUBLE,     /* speed_old[3] */
                                    MPI_DOUBLE,     /* density      */
                                    MPI_FLOAT,      /* tinj         */ 
                                    MPI_FLOAT,      /* color        */ 
                                    MPI_INT,	      /* cell[3]      */
                                    MPI_UINT32_T};  /* id           */

  int blocklengths[N_ELEMENTS]  = {3,3,3,3,1,1,1,3,1};
  MPI_Aint offsets[N_ELEMENTS];

  i=0;
  offsets[i++] = offsetof(Particle, coord);
  offsets[i++] = offsetof(Particle, speed);
  offsets[i++] = offsetof(Particle, coord_old);
  offsets[i++] = offsetof(Particle, speed_old);
  offsets[i++] = offsetof(Particle, density);
  offsets[i++] = offsetof(Particle, tinj);
  offsets[i++] = offsetof(Particle, color);
  offsets[i++] = offsetof(Particle, cell);
  offsets[i++] = offsetof(Particle, id);
  
  MPI_Type_create_struct(N_ELEMENTS, blocklengths, offsets, types, &MPI_PARTICLE);
  MPI_Type_commit(&MPI_PARTICLE);
  
#elif (PARTICLES_TYPE == LAGRANGIAN) && (PARTICLES_LP_SPECTRA == YES)
  #define N_ELEMENTS  33
  MPI_Datatype types[N_ELEMENTS] = {MPI_DOUBLE,   /* coord[3]     */
                                    MPI_DOUBLE,   /* speed[3]     */
                                    MPI_DOUBLE,   /* coord_old[3] */
                                    MPI_DOUBLE,   /* speed_old[3] */
                                    MPI_DOUBLE,   /* density      */
                                    
                                    MPI_DOUBLE,   /* fourvel[4]   */
                                    MPI_DOUBLE,   /* density_old  */
                                    MPI_DOUBLE,   /* pressure     */
                                    MPI_DOUBLE,   /* pressure_old */
                                    MPI_DOUBLE,   /* shk_gradp    */
                                    MPI_DOUBLE,   /* divv         */
                                    
                                    MPI_DOUBLE,   /* ca_old       */
                                    MPI_DOUBLE,   /* ca           */
                                    MPI_DOUBLE,   /* cr_old       */
                                    MPI_DOUBLE,   /* cr           */
                                    MPI_DOUBLE,   /* cmp_ratio    */
                                    MPI_DOUBLE,   /* lorG         */
                                    MPI_DOUBLE,   /* lorG_old     */
                                    
                                    MPI_DOUBLE,   /* nmicro       */
                                    MPI_DOUBLE,   /* mag[3]       */
                                    MPI_DOUBLE,   /* gradp[3]     */
                                    MPI_DOUBLE,   /* shk_vL[NVAR] */
                                    MPI_DOUBLE,   /* shk_vR[NVAR] */
                                    
                                    MPI_DOUBLE,   /* eng [EBINS] */
                                    MPI_DOUBLE,   /* chi [EBINS] */
                                    
                                    MPI_DOUBLE,   /* eng_old [EBINS] */
                                    MPI_DOUBLE,   /* chi_old [EBINS] */
                                    
                                    MPI_FLOAT,    /* tinj         */ 
                                    MPI_FLOAT,    /* color        */ 
                                    MPI_INT,	  /* cell[3]      */
                                    MPI_UINT32_T, /* id           */

                                    MPI_CHAR,     /*   shkflag    */
                                    MPI_CHAR,     /*   prev_shkflag    */
   				    };

  int blocklengths[N_ELEMENTS]  = {3, 3, 3, 3, 1, 4,
                                   1, 1, 1, 1, 1, 1, 
                                   1, 1, 1, 1, 1, 1, 1,
                                    3, 3, NVAR, NVAR,
                                   PARTICLES_LP_NEBINS, PARTICLES_LP_NEBINS,
                                   PARTICLES_LP_NEBINS, PARTICLES_LP_NEBINS,
                                   1, 1, 3, 1, 1, 1};

  MPI_Aint offsets[N_ELEMENTS];

  i = 0;
  offsets[i++] = offsetof(Particle, coord);
  offsets[i++] = offsetof(Particle, speed);
  offsets[i++] = offsetof(Particle, coord_old);
  offsets[i++] = offsetof(Particle, speed_old);
  offsets[i++] = offsetof(Particle, density);
  offsets[i++] = offsetof(Particle, fourvel);
  offsets[i++] = offsetof(Particle, density_old);
  offsets[i++] = offsetof(Particle, pressure);
  offsets[i++] = offsetof(Particle, pressure_old);
  offsets[i++] = offsetof(Particle, shk_gradp);
  offsets[i++] = offsetof(Particle, divv);
  offsets[i++] = offsetof(Particle, ca_old);
  offsets[i++] = offsetof(Particle, ca);
  offsets[i++] = offsetof(Particle, cr_old);
  offsets[i++] = offsetof(Particle, cr);
  offsets[i++] = offsetof(Particle, cmp_ratio);
  offsets[i++] = offsetof(Particle, lorG);
  offsets[i++] = offsetof(Particle, lorG_old);
  offsets[i++] = offsetof(Particle, nmicro);
  offsets[i++] = offsetof(Particle, mag);
  offsets[i++] = offsetof(Particle, gradp);
  offsets[i++] = offsetof(Particle, shk_vL);
  offsets[i++] = offsetof(Particle, shk_vR);
  offsets[i++] = offsetof(Particle, eng);
  offsets[i++] = offsetof(Particle, chi);
  offsets[i++] = offsetof(Particle, eng_old);
  offsets[i++] = offsetof(Particle, chi_old);
  offsets[i++] = offsetof(Particle, tinj);
  offsets[i++] = offsetof(Particle, color);
  offsets[i++] = offsetof(Particle, cell);
  offsets[i++] = offsetof(Particle, id);
  offsets[i++] = offsetof(Particle, shkflag);
  offsets[i++] = offsetof(Particle, prev_shkflag);
    
  MPI_Type_create_struct(N_ELEMENTS, blocklengths, offsets, types, &MPI_PARTICLE);
  MPI_Type_commit(&MPI_PARTICLE);
#endif
}
#endif

/* ********************************************************************* */
void Particles_SetID (particleNode *PHead)
/*!
 *
 * Assign id to newly injected particles.
 * The id is incrementally assigned starting from the global counter
 * p_idCounter.
 * 
 *     |--------------|--------------|--------------|--------------|
 * [A]      1,2,3           4,5            6             7,8        = N0
 *
 * [B]   +dn(0)=2        +dn(1)=5        +dn(2)=3      +dn(3)=4   
 *
 * The new id must take the value:
 *
 * rank #0:    1,2,3, [N0+1, N0+dn(0)]
 * rank #1:    4,5,   [N0+dn(0)+1, N0+dn(0)+dn(1)+1]
 * rank #2:    6,     [N0+dn(0)+dn(1)+1, N0+dn(0)+dn(1)+dn(2)+1]
 * ...
 *
 *********************************************************************** */
{
  long int Np0_loc, Np0_glob;
  int np, dn;
  int count, offset;
  static int *dn_arr;
  particleNode *curNode;
  Particle *p;

#ifdef PARALLEL
//print (">> [Particles_ID]\n");

/* --------------------------------------------------------
   0. Get number of processors and initialize arrays
   -------------------------------------------------------- */

  if (dn_arr == NULL) dn_arr = ARRAY_1D(g_nprocs, int);
  
/* --------------------------------------------------------
   1. Count how many particles have been injected (id = -1)
   -------------------------------------------------------- */

  dn = 0;
  PARTICLES_LOOP(curNode, PHead){
    p = &(curNode->p);
    if (p->id == -1)  dn++;
  }
/*
print ("  global idCounter           [before]: %d\n", p_idCounter);
print ("  local  number of particles [before]: %d\n", p_nparticles-dn);
print ("  local  increment                   : %d\n", dn);
print ("  local  number of particles [after]:  %d\n", p_nparticles);
*/
/* --------------------------------------------------------
   2. In parallel, do an MPI_Gather operation so the
      increments of all processors are gathered.
   -------------------------------------------------------- */
  
  MPI_Allgather (&dn, 1, MPI_INT, dn_arr, 1, MPI_INT, MPI_COMM_WORLD);
//  for (np = 0; np < g_nprocs; np++) print ("dn[%d] = %d\n",np, dn_arr[np]);
  
/* --------------------------------------------------------
   3. Set id to newly added particles. 
   -------------------------------------------------------- */

  offset = p_idCounter;
  for (np = 0; np < prank; np++) offset += dn_arr[np];
  count  = offset;
  PARTICLES_LOOP(curNode, PHead){
    p = &(curNode->p);
    if (p->id == -1){
      p->id = ++count;  
//print ("  assigning new id = %d\n",count);
    }
  }

/* --------------------------------------------------------
   4. Increment counter
   -------------------------------------------------------- */

  for (np = 0; np < g_nprocs; np++){
    p_idCounter += dn_arr[np];
  }
//print ("  global idCounter           [after]:  %d\n", p_idCounter);

#else

  count  = p_idCounter;
  PARTICLES_LOOP(curNode, PHead){
    p = &(curNode->p);
    if (p->id == -1){
      p->id = ++count;  
    }
  }
  p_idCounter = count;

#endif
}
