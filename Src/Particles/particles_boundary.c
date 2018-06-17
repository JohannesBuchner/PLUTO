/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Set boundary conditions on particles.
 
  Set boundary conditions at physical boundaries and exchange
  particles between neighbour processors.
 
  The function Particles_Boundary() sets boundary conditions at
  physical boundaries only. The function Particles_BoundaryExchange()
  exchanges particles between adjacent processors.
 
  \authors   A. Mignone (mignone@ph.unito.it)\n
             B. Vaidya (bvaidya@unito.it)\n
  
  \date      June 26, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Particles_Boundary(Data *d, Grid *grid)
/*!
 * Sets physical boundary conditions on particles.
 *
 * \param [in]  d         Pointer to PLUTO data structure.
 * \param [in]  grid      Pointer to PLUTO grid structure.        
 *********************************************************************** */
{
  int      dir, ncopyL, ncopyR;
  double   xbeg, xend, L[3];
  particleNode *curr, *next;
  Particle *p;

  DEBUG_FUNC_BEG ("Particles_Boundary");

/* --------------------------------------------------------
   0. Do some compatbility checking...
   -------------------------------------------------------- */

#ifdef SHEARINGBOX
  if (grid->nproc[JDIR] != 1) {
    print ("! Particles_Boundary(): ShearingBox+Particles requires -no-x2par\n");
    QUIT_PLUTO(1);
  }
#endif

  
  DIM_LOOP(dir) L[dir] = grid->xend_glob[dir] - grid->xbeg_glob[dir];

/* --------------------------------------------------------
   1. Loop over physical boudnaries
   -------------------------------------------------------- */

  for (dir = 0; dir < DIMENSIONS; dir++){
    curr   = d->PHead;
    ncopyL = ncopyR = 0;
    xbeg   = grid->xbeg_glob[dir];
    xend   = grid->xend_glob[dir];

  /* ----------------------------------------------
     1a.Process left physical boundary
     ---------------------------------------------- */
  
    switch (grid->lbound[dir]){
    case 0:  /* Not a physical boundary */
      break;
    case AXISYMMETRIC: case EQTSYMMETRIC: case REFLECTIVE:
      PARTICLES_LOOP(curr, d->PHead) {  
        p    = &(curr->p);        
        if (p->coord[dir] < xbeg) {
          p->coord[dir]   = 2.0*xbeg - p->coord[dir];
          p->speed[dir] *= -1.0;          
        }
      }
      break;
    case PERIODIC:
      #ifndef PARALLEL
      PARTICLES_LOOP(curr, d->PHead) {  /* Loop on particles */
        p    = &(curr->p);  
        if (p->coord[dir] < xbeg) p->coord[dir]  += L[dir];
      }
      #endif
      break;
    case SHEARING:
      #ifdef SHEARINGBOX
      PARTICLES_LOOP(curr, d->PHead) {  /* Loop on particles */
        p    = &(curr->p);  
        if (p->coord[dir] < xbeg) {
          double delta_Ly = fmod(sb_vy*g_time, L[JDIR]);
          #ifndef PARALLEL
          p->coord[IDIR]  += L[IDIR];   /* Shift to right side */
          #endif
          #ifndef FARGO  
          p->coord[JDIR]  -= delta_Ly;
          p->speed[JDIR] -= sb_vy;
          #endif   
        }
      }
      #endif
      break;
    case USERDEF:
      if      (dir == IDIR) Particles_UserDefBoundary (d, X1_BEG, grid);
      else if (dir == JDIR) Particles_UserDefBoundary (d, X2_BEG, grid);
      else if (dir == KDIR) Particles_UserDefBoundary (d, X3_BEG, grid);
      break;
    default:
      curr = d->PHead;
      while (curr != NULL) {  /* Do not use macro looping here  */
        p    = &(curr->p);  
        next = curr->next; 
        if (p->coord[dir] < xbeg) Particles_Destroy (curr,d);
        curr = next;
      }
      break;
    }          

  /* ----------------------------------------------
     1b.Process right physical boundary
     ---------------------------------------------- */
  
    switch (grid->rbound[dir]){
    case 0:  /* Not a physical boundary */
      break;
    case AXISYMMETRIC: case EQTSYMMETRIC: case REFLECTIVE:
      PARTICLES_LOOP(curr, d->PHead) {  
        p    = &(curr->p);

        if (p->coord[dir] >= xend) {
          p->coord[dir]   = 2.0*xend - p->coord[dir];
          p->speed[dir] *= -1.0;          
        }
      }
      break;
    case PERIODIC:
      #ifndef PARALLEL
      PARTICLES_LOOP(curr, d->PHead) {  
        p = &(curr->p);
        if (p->coord[dir] >=  xend) p->coord[dir]  -= L[dir];
      }
      #endif
      break;
    case SHEARING:
      #ifdef SHEARINGBOX
      PARTICLES_LOOP(curr, d->PHead) {  /* Loop on particles */
        p    = &(curr->p);

        if (p->coord[dir] >=  xend) {
          double delta_Ly = fmod(sb_vy*g_time, L[JDIR]);
          #ifndef PARALLEL
          p->coord[IDIR]  -= L[IDIR];    /* Shift to left side of domain */
          #endif
          #ifndef FARGO   
          p->coord[JDIR]  += delta_Ly;
          p->speed[JDIR] += sb_vy;
          #endif  
        }
      }
      #endif
      break;
    case USERDEF:
      if      (dir == IDIR) Particles_UserDefBoundary (d, X1_END, grid);
      else if (dir == JDIR) Particles_UserDefBoundary (d, X2_END, grid);
      else if (dir == KDIR) Particles_UserDefBoundary (d, X3_END, grid);
      break;
    default:
      curr = d->PHead;
      while (curr != NULL) {  /* Do not use macro looping here  */
        p    = &(curr->p);  
        next = curr->next; 
        if (p->coord[dir] >= xend) Particles_Destroy (curr,d);
        curr = next;
      }
      break;
    }          
  } /* End loop on dimensions */

  
#if DEBUG      
if (d_condition)print("Particles_Check = %d, %d, %d; p_nparticles = %d\n",
        Particles_CheckAll(d->PHead, 0, grid),
        Particles_CheckAll(d->PHead, 1, grid),
        Particles_CheckAll(d->PHead, 2, grid), p_nparticles);
#endif

  DEBUG_FUNC_END ("Particles_Boundary");
}

/* ********************************************************************* */
void Particles_BoundaryExchange(Data *d, Grid *grid)
/*!
 * In parallel mode, exchange/duplicate particles between neighbour
 * processor.
 * The function is done direction-wise: if a particle escapes from the
 * left (right) boundary it is sent to the left (right) processor.
 * The exchange is repeated direction by direction.
 *
 * \param [in]   d        Pointer to PLUTO data structure.
 * \param [in]   grid     Pointer to the PLUTO grid structure.
 *********************************************************************** */
{
#ifdef PARALLEL
  int dir, i, beg, end;
  int nsendL, nsendR;
  int nrecvL, nrecvR;
  int procL, procR;     /* Rank of left and right neighbour processors */
  int success;
  static int first_call = 1, **neigh;

  double xbeg, xend;

  particleNode *curr, *next;
  Particle *p;         /* Just a shortcut */
  Particle *recv_bufL, *recv_bufR;  /* Array of received particles (runtime) */
  Particle *send_bufL, *send_bufR;  /* Array collecting particles to be sent (runtime) */ 
  particleNode **pstackL, **pstackR; /* Array collecting particle pointers to
                                        be trasnferred (runtime) */
  MPI_Status status;

  DEBUG_FUNC_BEG("Particles_BoundaryExchange");

/* --------------------------------------------------------------
   1. On first call, finds processor's neighbour ranks from
      Cartesian topology and the corresponding domain coordinates
      for each direction.
   -------------------------------------------------------------- */       

  if (first_call){
    neigh = ARRAY_2D(3, 2, int);    
    GetNeighbourRanks (grid, neigh);
    first_call = 0;       
  } 
  MPI_Barrier (MPI_COMM_WORLD);

/* -------------------------------------------------------------
   2. Loop on directions and set boundary conditions
      direction-wise.      
   ------------------------------------------------------------- */
        
// long int nsend_tot=0, nrecv_tot=0; 
  DIM_LOOP(dir){

  /* ------------------------------------------
     2a. Allocate stack buffer
         Note: p_nparticles may change from
               one direction to the next.
     ------------------------------------------ */

    pstackL = ARRAY_1D(p_nparticles, particleNode *);
    pstackR = ARRAY_1D(p_nparticles, particleNode *);
    
    procR = neigh[dir][1]; /* Rank of processor to the right */
    procL = neigh[dir][0]; /* Rank of processor to the left */

    nsendL = nsendR = 0;

  /* -------------------------------
     2b. Set boundary region
     ------------------------------- */

    beg  = grid->lbeg[dir];
    end  = grid->lend[dir];

    xbeg = grid->xl[dir][beg]; 
    xend = grid->xr[dir][end];
    
  /* -------------------------------------------------------------
     2c. Count how many particles are leaving the local processor
         domain from either the left or right side and stack
         pointers into the arrays pstackL[] and pstackR[].
     ------------------------------------------------------------- */

    curr = d->PHead;
    while (curr != NULL){   /* Loop on particles */
       
      next = curr->next;
      p    = &(curr->p);

      if      (p->coord[dir] <  xbeg) pstackL[nsendL++] = curr;
      else if (p->coord[dir] >= xend) pstackR[nsendR++] = curr;

      curr = next;
    }  /* End loop on particles */

  /* -------------------------------------------------------------
     2d. Allocate memory for send buffers and copy particles
         from stacks 
     ------------------------------------------------------------- */

    send_bufL = ARRAY_1D(nsendL, Particle);
    send_bufR = ARRAY_1D(nsendR, Particle);

    for (i = 0; i < nsendL; i++) send_bufL[i] = pstackL[i]->p; 
    for (i = 0; i < nsendR; i++) send_bufR[i] = pstackR[i]->p; 

  /* ---------------------------------------------------------------
     2e. Exchange particles between neighbouring processors. 
         Data is transferred only between neighbours and only if 
         one or more particles are to be sent/received.
         No communication should be done if a processor touches
         a non-periodic physical boundary.
     --------------------------------------------------------------- */
     
    nrecvL = 0;
    nrecvR = 0;
    if (procR < 0) {
      if (nsendR != 0){
        print ("! Particles_BoundaryExchange(): nsendR != 0\n");
        QUIT_PLUTO(1);
      }
      nsendR = 0;
      procR  = MPI_PROC_NULL;
    }
    if (procL < 0) {
      if (nsendL != 0){
        print ("! Particles_BoundaryExchange(): nsendL != 0\n");
        QUIT_PLUTO(1);
      }
      nsendL = 0;
      procL  = MPI_PROC_NULL;
    }

    MPI_Sendrecv(&nsendR, 1, MPI_INT, procR, 1, 
                 &nrecvL, 1, MPI_INT, procL, 1, MPI_COMM_WORLD, &status);

    MPI_Sendrecv(&nsendL, 1, MPI_INT, procL, 5, 
                 &nrecvR, 1, MPI_INT, procR, 5, MPI_COMM_WORLD, &status);

    recv_bufL = ARRAY_1D(nrecvL, Particle);
    recv_bufR = ARRAY_1D(nrecvR, Particle);
   
#if PARTICLES_USE_MPI_DATATYPE == YES   
    MPI_Sendrecv(send_bufR, nsendR, MPI_PARTICLE, procR, 7, 
                 recv_bufL, nrecvL, MPI_PARTICLE, procL, 7, MPI_COMM_WORLD, &status);
#else
    MPI_Sendrecv(send_bufR, sizeof(Particle)*nsendR, MPI_BYTE, procR, 7, 
                 recv_bufL, sizeof(Particle)*nrecvL, MPI_BYTE, procL, 7, MPI_COMM_WORLD, &status);
#endif

    for(i = 0; i < nrecvL; i++){
      p = recv_bufL + i;
      if (grid->lbound[dir] == PERIODIC || grid->lbound[dir] == SHEARING){
        p->coord[dir] -= g_domEnd[dir] - g_domBeg[dir]; /* Apply Periodicity. */
      }
      success = Particles_Insert(p, d, PARTICLES_TRANSFER, grid);
    }
     
#if PARTICLES_USE_MPI_DATATYPE == YES
    MPI_Sendrecv(send_bufL, nsendL, MPI_PARTICLE, procL , 9, 
                 recv_bufR, nrecvR, MPI_PARTICLE, procR,  9, MPI_COMM_WORLD, &status);
#else
    MPI_Sendrecv(send_bufL, sizeof(Particle)*nsendL, MPI_BYTE, procL , 9, 
                 recv_bufR, sizeof(Particle)*nrecvR, MPI_BYTE, procR,  9, MPI_COMM_WORLD, &status);
#endif

    for(i = 0; i < nrecvR; i++){
      p = recv_bufR + i;
      if (grid->rbound[dir] == PERIODIC || grid->rbound[dir] == SHEARING){
        p->coord[dir] += g_domEnd[dir] - g_domBeg[dir]; /* Apply Periodicity. */
      }
      success = Particles_Insert(p, d, PARTICLES_TRANSFER, grid);
    }     

    MPI_Barrier(MPI_COMM_WORLD); /* Synchronize after each sweep */

/* -------------------------------------------------------
    Destroy particles that have been transferred
   ------------------------------------------------------- */

    for (i = 0; i < nsendL; i++) Particles_Destroy (pstackL[i], d); 
    for (i = 0; i < nsendR; i++) Particles_Destroy (pstackR[i], d); 

/* --------------------------------------------------------
     Free memory
   -------------------------------------------------------- */
      
    FreeArray1D(pstackL);
    FreeArray1D(pstackR);
    FreeArray1D(recv_bufL);
    FreeArray1D(recv_bufR);
    FreeArray1D(send_bufL);
    FreeArray1D(send_bufR);

    g_usedMemory -= 2*p_nparticles*sizeof(particleNode *);
    g_usedMemory -= (nsendL+nsendR)*sizeof(Particle);
    g_usedMemory -= (nrecvL+nrecvR)*sizeof(Particle);

  }  /* End loop on dimensions */

#if DEBUG    
if (d_condition)print("Particles_Check = %d, %d, %d; p_nparticles = %d\n",
        Particles_CheckAll(d->PHead, 0, grid),
        Particles_CheckAll(d->PHead, 1, grid),
        Particles_CheckAll(d->PHead, 2, grid), p_nparticles);
#endif
   
#endif /* PARALLEL */
/*
  int long count = Particles_CheckAll(d->PHead, 1, grid);
  if (count != p_nparticles){
    print ("! Particles_BoundaryExchange(): %d particles out of %d are\n",
           count,p_nparticles);
    print ("! inside the active domain. \n");
    QUIT_PLUTO(1);
  }
*/
#if PARTICLES_USE_ARRAY == YES
  Particles_ListToArray(d);
#endif

  DEBUG_FUNC_END("Particles_BoundaryExchange");
}

/* ********************************************************************* */
void Particles_ListToArray (Data *d)
/*!
 * Convert a linked list into an array of pointers.
 * This function must be called only after particles have been 
 * redistributed.
 * 
 *********************************************************************** */
{
  long int count;
  particleNode *curr, *next;

/* ------------------------------------------------
   0. Destroy current array and create a new one
   ------------------------------------------------ */

  if (d->pstr != NULL) FreeArray1D(d->pstr);
  d->pstr = ARRAY_1D(p_nparticles+1, Particle *); 

  count = 0L;
  curr  = d->PHead;
  while (curr != NULL){   /* Loop on particles */
    d->pstr[count++] = &(curr->p);
    next = curr->next;
    curr = next;
  }  /* End loop on particles */
  d->pstr[count] = NULL; /* -- Set last element to NULL -- */

  if (count != p_nparticles){
    print ("! Particles_ListToArray(): number of particles does not match\n");
    QUIT_PLUTO(1);
  }

  return;
/* -- Verify -- */

  int i;
  Particle *p;
  Particle **node = d->pstr;
  
  PARTICLES_LOOP(curr, d->PHead){
    p = &(curr->p);
    Particles_Display(p);
  }

/*
  for (i = 0; i < p_nparticles; i++){
    p = d->pstr[i];
    Particles_Display(p);
  }
*/
/*
  for (p = d->pstr[0], i=0; i < p_nparticles; i++, p=d->pstr[i]){
    Particles_Display(p);
  }
*/
  for (node = d->pstr, p = node[0]; p != NULL; node++, p = node[0]){
//  for (p = d->pstr[0]; p != NULL; (*p)++){
    Particles_Display(p);
  }




  print ("Done in Particles_ListToArray\n");
  exit(1);
}

/* ********************************************************************* */
int Particles_BoundaryCheck(Particle *p, Grid *grid)
/*!
 * Check if a particle is in the ghost zones or outside the whole
 * computational domain.
 *
 *
 *      <--------->                               <--------->
 *     | last |      |  .... Active Domain ....    |      | last |
 *     
 *    xb2    xb1    xb0                           xe0    xe1    xe2
 *
 *  Here 'last' marks the very last (or first) ghost zones.
 *  If a particle lies here, then deposition may not be done
 *  correctly.
 *  
 *  \return This function returns:
 *          - 0  if the particle is in the active domain
 *          - 1  if the particle is in the ghost zones, but in the last one
 *          - 2  if the particle is in the last ghost zone
 *          - 3  if the particle is outside the whole comput. domain
 *********************************************************************** */
{
  int dir;
  int beg0, beg1, beg2;
  int end0, end1, end2, np_tot;
  double xb0, xb1, xb2;
  double xe0, xe1, xe2;
  
  for (dir = 0; dir < DIMENSIONS; dir++){
    np_tot = grid->np_tot[dir];
    
    beg2 = 0;
    end2 = np_tot-1;

    beg1 = beg2 + 1;
    end1 = end2 - 1;

    beg0 = grid->lbeg[dir];
    end0 = grid->lend[dir];
    
    xb0 = grid->xl[dir][beg0];
    xe0 = grid->xr[dir][end0];

    xb1 = grid->xl[dir][beg1];
    xe1 = grid->xr[dir][end1];

    xb2 = grid->xl[dir][beg2];
    xe2 = grid->xr[dir][end2];
    
    if (p->coord[dir] < xb2 || p->coord[dir] > xe2){
      return 3;
    }else if (p->coord[dir] < xb1 || p->coord[dir] > xe1){
      return 2;
    }else if (p->coord[dir] < xb0 || p->coord[dir] > xe0){
      return 1;
    }
  }

  return 0;
}

