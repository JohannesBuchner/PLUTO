/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Particle to grid deposition functions.

  Transfer particle quantities to the grid using deposition:
  \f[
     Q_{ijk} = \sum_p W_p(x_{ijk}-x_p)q_p
  \f]
  where \c qp is a particle-related quantities (e.g. energy or mass).
  For deposition to be consistently done, particles can be everywhere in
  the domain, including ghost zones with the important exception of
  the very first and very last one.
  An error will occur in Particles_GetWeights() otherwise.
  
  Alternatively, integer deposition can be used to verify the
  correctness of serial/parallel implementation:
  \f[
     Q_{ijk} = \frac{q_{\max}}{C} \sum_p {\rm int}\left(C
               W_p(x_{ijk}-x_p) \frac{q_p}{q_{\max}}\right)
  \f]
  where \c C is a large number (\c 1.e6 by default) while
  \c q_max is the maximum value of \c qp over all particles.

  Deposition is always followed by a boundary condition call,
  since quantities deposited in the ghost zones must be exchanged
  with neighbour processors. 

  \authors A. Mignone (mignone@ph.unito.it)\n

  \b References

  \date  March 31, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define NELEM_MAX     64
/* ********************************************************************* */
void Particles_Deposit(particleNode *PHead, void (*Func)(Particle *, double *),
                       Data_Arr Q, int nelem, Grid *grid)
/*!
 *
 *********************************************************************** */
{
  int    i, j, k, n, dir;
  int    i1, j1, k1;
  long int np_tot = NX1_TOT*NX2_TOT*NX3_TOT;
  double qd[NELEM_MAX];
  static double ***w;
  particleNode *CurNode;
  Particle     *p;

/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */

  if (w == NULL) {
    w  = ArrayBox (-1, 1, -1, 1, -1, 1);
  }

  if (nelem > NELEM_MAX){
    print ("! Particles_Deposit: exceeded max number of elements (%d)\n",
           NELEM_MAX);
    QUIT_PLUTO(1);
  }

/* --------------------------------------------------------
   1. Initialize deposition: call boundary and
      set array elements to zero
   -------------------------------------------------------- */

  for (n = 0; n < nelem; n++) {
    memset ((void *)Q[n][0][0], '\0', np_tot*sizeof(double));
  }

/* --------------------------------------------------------
   2. Loop on particles and deposit to grid zones
   -------------------------------------------------------- */

#if PARTICLES_DEPOSIT == REAL
  PARTICLES_LOOP(CurNode, PHead){
    p = &(CurNode->p);

    Particles_GetWeights(p, p->cell, w, grid);

    Func(p, qd);    /* Compute quantities to be deposited */
    i = p->cell[IDIR];
    j = p->cell[JDIR];
    k = p->cell[KDIR];

    for (n = 0; n < nelem; n++){
      for (k1 = -KOFFSET; k1 <= KOFFSET; k1++){
      for (j1 = -JOFFSET; j1 <= JOFFSET; j1++){
      for (i1 = -IOFFSET; i1 <= IOFFSET; i1++){
        Q[n][k+k1][j+j1][i+i1] += qd[n]*w[k1][j1][i1];
      }}}
    }
  }
#elif PARTICLES_DEPOSIT == INTEGER
  long int wq;
  double C = 1.e6, max_qd[16], glob_max_qd[16];

/* --------------------------------------------------------
   2a. Compute the maximum value, among all particles, of
       the quantities to be deposited.
   -------------------------------------------------------- */
      
  for (n = 0; n < nelem; n++)  max_qd[n] = 0.0;
  PARTICLES_LOOP(CurNode, PHead){
    p = &(CurNode->p);

    Func(p, qd);    
    for (n = 0; n < nelem; n++) max_qd[n] = MAX(max_qd[n], fabs(qd[n]));
  }

  #ifdef PARALLEL
  MPI_Allreduce (max_qd, glob_max_qd, nelem, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  for (n = 0; n < nelem; n++)  max_qd[n] = glob_max_qd[n];
  #endif

/* --------------------------------------------------------
   2b. Deposit normalized values (qd/max(qd)) by
       transforming to integer, so that addition is
       associative.
   -------------------------------------------------------- */
  
  PARTICLES_LOOP(CurNode, PHead){
    p = &(CurNode->p);

    Particles_GetWeights(p, p->cell, w, grid);

    Func(p, qd);    /* Compute quantities to be deposited */    
    i = p->cell[IDIR];
    j = p->cell[JDIR];
    k = p->cell[KDIR];

    for (n = 0; n < nelem; n++){
      qd[n] /= max_qd[n];
      for (k1 = -KOFFSET; k1 <= KOFFSET; k1++){
      for (j1 = -JOFFSET; j1 <= JOFFSET; j1++){
      for (i1 = -IOFFSET; i1 <= IOFFSET; i1++){
        wq = (long)(C*qd[n]*w[k1][j1][i1]);
        Q[n][k+k1][j+j1][i+i1] += wq;
      }}}
    }
  }
#endif

/* --------------------------------------------------------
   3. Exchange ghost zones between processors, or
      impose periodicity (if required).
   -------------------------------------------------------- */

  Particles_DepositBoundaryExchange(Q, nelem, grid);

/* --------------------------------------------------------
   3a. Normalize back to real. 
   -------------------------------------------------------- */
  
#if PARTICLES_DEPOSIT == INTEGER
  for (n = 0; n < nelem; n++) TOT_LOOP(k,j,i) Q[n][k][j][i] *= max_qd[n]/C;
#endif

//print ("<<[Particles_Deposit()]\n");

}

/* ********************************************************************* */
void Particles_DepositBoundaryExchange(Data_Arr Q, int nelem, Grid *grid)
/*!
 * Exchange boundary values between processors to complete the
 * deposition process.
 * If Q contains quantities in the boundaries, we first copy these
 * values into auxiliary buffers \c snd_bufL and \c snd_bufR and then
 * transfer them through an MPI Send/Recv call:
 *
 * \code
 *   Rank #0        A0|B0------C0|D0
 *                             ^   |
 *                             |   v 
 *   Rank #1                   A1|B1------C1|D1
 *                                        ^   |
 *                                        |   v
 *   Rank #2                              A2|B2------C2|D2
 * \endcode
 *
 * Seen from rank # 1, we have
 * - \c A1 = snd_bufL
 * - \c D1 = snd_bufR
 * - \c D0 = rcv_bufL
 * - \c A2 = rcv_bufR
 *
 * Here 'L' or 'R' refer to the processor to the left (rank #0) and to
 * the right (rank #2), respectively.
 * In 2D or 3D we iterate over dimensions.
 *********************************************************************** */
{
  int i,j,k,n, dir;
  int ngh, procL, procR;
  int ibeg, jbeg, kbeg;
  int iend, jend, kend;
  static int **neigh;
  long int countL, countR, nsndL, nsndR, nrcvL, nrcvR;
  static double  *snd_bufL, *snd_bufR;
  static double  *rcv_bufL, *rcv_bufR;
  RBox boxL, boxR;
  
/* ---------------------------------------------
   0. Allocate memory
   --------------------------------------------- */

#ifdef PARALLEL
  if (snd_bufL == NULL) {
    long size_max;
    size_max  = MAX(NX1_TOT*NX2_TOT, NX1_TOT*NX3_TOT);
    size_max  = MAX(size_max, NX2_TOT*NX3_TOT);
    size_max *= grid->nghost[IDIR];
    
    snd_bufL = ARRAY_1D(size_max, double);
    snd_bufR = ARRAY_1D(size_max, double);
    rcv_bufL = ARRAY_1D(size_max, double);
    rcv_bufR = ARRAY_1D(size_max, double);
    neigh    = ARRAY_2D(3, 2, int);    

    GetNeighbourRanks (grid, neigh);
  }
  MPI_Status status;
  
  for (n = 0; n < nelem; n++){
    for (dir = 0; dir < DIMENSIONS; dir++){

      procR = neigh[dir][1]; /* Rank of processor to the right */
      procL = neigh[dir][0]; /* Rank of processor to the left */

      ngh = grid->nghost[dir];
      
    /* -- Define left box -- */
    
      ibeg = 0; iend = (dir == IDIR ? IBEG-1:NX1_TOT-1);
      jbeg = 0; jend = (dir == JDIR ? JBEG-1:NX2_TOT-1);
      kbeg = 0; kend = (dir == KDIR ? KBEG-1:NX3_TOT-1);

      RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, CENTER, &boxL);

    /* -- Define right box -- */

      ibeg = (dir == IDIR ? IEND+1:0); iend = NX1_TOT-1;
      jbeg = (dir == JDIR ? JEND+1:0); jend = NX2_TOT-1;
      kbeg = (dir == KDIR ? KEND+1:0); kend = NX3_TOT-1;
          
      RBoxDefine(ibeg, iend, jbeg, jend, kbeg, kend, CENTER, &boxR);

    /* -- Copy buffer -- */
      
      countL = countR = 0;
      BOX_LOOP(&boxL, k,j,i) snd_bufL[countL++] = Q[n][k][j][i];
      BOX_LOOP(&boxR, k,j,i) snd_bufR[countR++] = Q[n][k][j][i];
/*
if (countL != countR){
  print ("! err, countL = %d != countR = %d\n", countL, countR);
  QUIT_PLUTO(1);
}
*/
    /* -- Send left buffer, recv right buffer -- */

      nsndL = nrcvL = countL;
      nsndR = nrcvR = countR;
      if (procL < 0){
        nsndL = 0;
        nrcvL = 0;
        procL = MPI_PROC_NULL;
      }
      
      if (procR < 0){
        nsndR = 0;
        nrcvR = 0;
        procR = MPI_PROC_NULL;
      }
   
      MPI_Sendrecv (snd_bufL, nsndL, MPI_DOUBLE, procL, 1,
                    rcv_bufR, nrcvR, MPI_DOUBLE, procR, 1, MPI_COMM_WORLD, &status);

      MPI_Sendrecv (snd_bufR, nsndR, MPI_DOUBLE, procR, 2,
                    rcv_bufL, nrcvL, MPI_DOUBLE, procL, 2, MPI_COMM_WORLD, &status);
      
    /* -- Shift left buffer to active domain, add data -- */

      if (nrcvL > 0){
        boxL.ibeg += (dir == IDIR ? ngh:0); boxL.iend += (dir == IDIR ? ngh:0);
        boxL.jbeg += (dir == JDIR ? ngh:0); boxL.jend += (dir == JDIR ? ngh:0);
        boxL.kbeg += (dir == KDIR ? ngh:0); boxL.kend += (dir == KDIR ? ngh:0);
        countL = 0;                 
        BOX_LOOP(&boxL, k,j,i)  Q[n][k][j][i] += rcv_bufL[countL++];
      }  
      
    /* -- Shift right buffer to active domain -- */

      if (nrcvR > 0){
        boxR.ibeg -= (dir == IDIR ? ngh:0); boxR.iend -= (dir == IDIR ? ngh:0);
        boxR.jbeg -= (dir == JDIR ? ngh:0); boxR.jend -= (dir == JDIR ? ngh:0);
        boxR.kbeg -= (dir == KDIR ? ngh:0); boxR.kend -= (dir == KDIR ? ngh:0);
        countR = 0;                 
        BOX_LOOP(&boxR, k,j,i)  Q[n][k][j][i] += rcv_bufR[countR++];
      }  

    /* -- Add buffer to active zones values -- */

      int pdir[3];
      pdir[IDIR] = (dir == IDIR ? 1:0);
      pdir[JDIR] = (dir == JDIR ? 1:0);
      pdir[KDIR] = (dir == KDIR ? 1:0);
      AL_Exchange_dim ((char *)Q[n][0][0], pdir, SZ);
    }
  }  
  
#else

  #if DIMENSIONS >= 1
  if (   (grid->lbound[IDIR] == PERIODIC && grid->rbound[IDIR] == PERIODIC)
      || (grid->lbound[IDIR] == SHEARING && grid->rbound[IDIR] == SHEARING)){
    for (n = 0; n < nelem; n++){
      X1_BEG_LOOP(k,j,i) Q[n][k][j][i+NX1] += Q[n][k][j][i];
      X1_END_LOOP(k,j,i) Q[n][k][j][i-NX1] += Q[n][k][j][i];

      X1_BEG_LOOP(k,j,i) Q[n][k][j][i] = Q[n][k][j][i+NX1];
      X1_END_LOOP(k,j,i) Q[n][k][j][i] = Q[n][k][j][i-NX1];
    }
  }
  #endif

  #if DIMENSIONS >= 2
  if (grid->lbound[JDIR] == PERIODIC && grid->rbound[JDIR] == PERIODIC){
    for (n = 0; n < nelem; n++){
      X2_BEG_LOOP(k,j,i) Q[n][k][j+NX2][i] += Q[n][k][j][i];
      X2_END_LOOP(k,j,i) Q[n][k][j-NX2][i] += Q[n][k][j][i];

      X2_BEG_LOOP(k,j,i) Q[n][k][j][i] = Q[n][k][j+NX2][i];
      X2_END_LOOP(k,j,i) Q[n][k][j][i] = Q[n][k][j-NX2][i];
    }
  }
  #endif

  #if DIMENSIONS == 3
  if (grid->lbound[KDIR] == PERIODIC && grid->rbound[KDIR] == PERIODIC){
    for (n = 0; n < nelem; n++){
      X3_BEG_LOOP(k,j,i) Q[n][k+NX3][j][i] += Q[n][k][j][i];
      X3_END_LOOP(k,j,i) Q[n][k-NX3][j][i] += Q[n][k][j][i];

      X3_BEG_LOOP(k,j,i) Q[n][k][j][i] = Q[n][k+NX3][j][i];
      X3_END_LOOP(k,j,i) Q[n][k][j][i] = Q[n][k-NX3][j][i];
    }
  }
  #endif
  
#endif  /* PARALLEL */
}

/* ********************************************************************* */
void Particles_Density(Particle *p, double *qd)
/*!
 *  Simple deposition function used to compute particle density
 *  on the grid.
 *********************************************************************** */
{
#if PARTICLES_TYPE == DUST
  qd[0] = p->rho;
#elif PARTICLES_TYPE == COSMIC_RAYS
  qd[0] = p->rho;
#endif
}

