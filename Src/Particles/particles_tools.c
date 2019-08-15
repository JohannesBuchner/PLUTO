/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Tools required to define the Particle MPI struct and 
        Interpolating quantities from grid to particles.
 
 \authors   A. Mignone (mignone@ph.unito.it)\n
            B. Vaidya (bvaidya@unito.it)\n
  
 \date     March 30, 2018
 */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
long Particles_CheckAll (particleNode *PHead, int mode, Grid *grid)
/*!
 *  Count and return the number of particles inside a given region
 *  of the computational domain.
 *  The region is specified using mode ( = 0,1,2),
 *  see Particles_CheckSingle()
 *********************************************************************** */
{
  long int count=0, check;
  Particle     *p;
  particleNode *curr, *next;
  
  curr  = PHead;
  while (curr != NULL) {  /* Loop on particles */
    p    = &(curr->p);

    check = Particles_CheckSingle(p, mode, grid);
    if (check) count++;
    next = curr->next; 
    curr = next;
  }

  return count;
}

/* ********************************************************************* */
int Particles_CheckSingle(Particle *p, int mode, Grid *grid)
/*!
 * Check if the particle belongs to the local processor domain
 *
 *  \param [in]   p      pointer to Particle structure.
 *  \param [in]   mode   operation mode:
 *                       mode = 0   check if particle is in active+ghost domain;
 *                       mode = 1   check if particle is in active zones;
 *                       mode = 2   check just ghost zones;
 *
 * \return TRUE if the particle is inside the specified region.
 *         FALSE otherwise.
 *********************************************************************** */
{
  int    dir, bbeg, bend, abeg, aend, ngh;
  int    cond0, cond1;
  double xbeg0, xend0, xbeg1, xend1;
  
  if (mode < 0 || mode > 2){
    print ("! Particles_CheckSingle(): invalid mode %d\n",mode);
    QUIT_PLUTO(1);
  }

  cond0 = cond1 = 1;
  DIM_LOOP(dir){
    abeg = grid->lbeg[dir];
    aend = grid->lend[dir];
    ngh  = grid->nghost[dir];
    
    bbeg = abeg - ngh;
    bend = aend + ngh;

    xbeg0 = grid->xl[dir][bbeg];
    xend0 = grid->xr[dir][bend];

    xbeg1 = grid->xl[dir][abeg];
    xend1 = grid->xr[dir][aend];
    
    cond0 *= (p->coord[dir] <= xend0);
    cond0 *= (p->coord[dir] >  xbeg0);
    
    cond1 *= (p->coord[dir] <= xend1);
    cond1 *= (p->coord[dir] >  xbeg1);
    
  }
  
  if (mode == 0) return cond0;
  if (mode == 1) return cond1;
  if (mode == 2) return (cond1 == FALSE ? cond0: FALSE);
    
}

/* ********************************************************************* */
void Particles_GetWeights (Particle *p, int *cell, double ***w, Grid *grid)
/*! 
 * Compute particle weights such that 
 * \f[
 *      q_p = \sum_{ij} w_{ij} Q_{ij}
 * \f]
 * 
 * \param [in]  p       pointer to particle structure
 * \param [out] cell    a 3-element array containing the indices (i,j,k) of 
 *                      the grid cell that is closer and to the left of 
 *                      the particle.
 * \param [out] w       a 3x3x3 array containing the weights
 * \param [in]  grid    a pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int    err, i, j, k, dir;
  double *xg, *xp, *inv_dx;
  double w1[3][3], delta;

  err = Particles_LocateCell(p->coord, cell, grid);
  if (err){
    print ("! Particles_GetWeights(): index out of range  [%d, %d, %d\n",
             cell[IDIR], cell[JDIR], cell[KDIR]); 
    Particles_Display(p);
    QUIT_PLUTO(1);
  }
  
/* -- Set default value for 1D weights  -- */

  for (dir = 0; dir < 3; dir++) {
    w1[dir][0] = 0.0;
    w1[dir][1] = 0.0;
    w1[dir][2] = 0.0;
  }  
    
  DIM_LOOP(dir){

    xg     = grid->x[dir];
    xp     = p->coord;
    inv_dx = grid->inv_dx[dir];
    i      = cell[dir];

  /* ------------------------------------------------------
      Weights can be computed properly only if particles
      do not lie in the very first or very last ghost zone.
     ------------------------------------------------------ */
    
    if (i > 0 && i < grid->np_tot[dir]-1){
      double eps = 1.e-12;

      delta = (p->coord[dir] - xg[i])*inv_dx[i];  /* -1/2 <= delta < 1/2 */
      
      if ( delta < -(0.5+eps) || delta > (0.5+eps)){
        print ("! Particles_GetWeights(): delta+1/2 = %12.6e; deltaL-1/2= %12.6e\n",
                delta+0.5, delta-0.5);
        QUIT_PLUTO(1);
      }

#if PARTICLES_SHAPE == 1

/* -- Compute 1D weights:  "Nearest Grid Point" (NGP) method -- */

      w1[dir][0] = 0.0;
      w1[dir][1] = 1.0;
      w1[dir][2] = 0.0;

#elif PARTICLES_SHAPE == 2

/* -- Compute 1D weights:  "Cloud-In-Cell" (CIC) method -- */

      w1[dir][0] = MAX(0.0,- delta);
      w1[dir][1] = 1.0 - fabs(delta);
      w1[dir][2] = MAX(0.0,+ delta);

#elif PARTICLES_SHAPE == 3

  /* -- Compute 1D weights:  "Triangular-Shaped Cloud" (TSC) method -- */

      w1[dir][0] = 0.125*(1.0 - 2.0*delta)*(1.0 - 2.0*delta);
      w1[dir][1] = 0.75 - delta*delta;
      w1[dir][2] = 0.125*(1.0 + 2.0*delta)*(1.0 + 2.0*delta);

      w1[dir][1] = 1.0 - w1[dir][0] - w1[dir][2];

#else
      print ("! Particles_GetWeights(): invalid shape\n");
      QUIT_PLUTO(1);
#endif

    }else{
      print ("! Particles_GetWeights(): particle (%d)", p->id);
      print (" outside computational domain [dir = %d]\n",dir);
      print ("!                         pcell = (%d, %d, %d)\n",
              cell[IDIR], cell[JDIR], cell[KDIR]);
      
      print ("!                         xp = %f, Active = [%f, %f], full = [%f, %f]\n",
              p->coord[dir], grid->xbeg[dir], grid->xend[dir],
                            grid->xl[dir][0], grid->xr[dir][grid->np_tot[dir]-1]);
      Particles_Display(p);     
      QUIT_PLUTO(1); 
    }
  }  /* End loop on directions */

  for (k = -KOFFSET; k <= KOFFSET; k++) {
  for (j = -JOFFSET; j <= JOFFSET; j++) { 
  for (i = -IOFFSET; i <= IOFFSET; i++) {
    w[k][j][i] = D_EXPAND(w1[IDIR][i+1], *w1[JDIR][j+1], *w1[KDIR][k+1]);
  }}}   
    
}

/* ********************************************************************* */
double Particles_Interpolate(double ***V, double ***w, int *indx)
/*! 
 * Interpolate a grid quantity V to particle position x.
 *
 * \param [in]   V    a 3D array defined on the fluid grid->
 * \param [in]   w    a 3x3x3 array containing the weights
 * \param [in]  indx  a 3 element array giving the starting indices 
 *                    (i,j,k). 
 *
 * \return The interpolated value of V at the desired location.
 *********************************************************************** */
{
  int    i, j, k;
  int    i1, j1, k1;
  double v = 0.0;
  
  i = indx[IDIR];
  j = indx[JDIR];
  k = indx[KDIR];
 
/* --------------------------------------------------------
    Safety check: ensure particle lies inside the active 
    computational zone or, at least, in the 1st ghost zone 
   -------------------------------------------------------- */
/*      
  if (i < 1 || i > NX1_TOT-2){
    print ("! Particles_Interpolate(): i = %d out of bound\n", i);
    print ("  p_intStage = %d\n",p_intStage);
    QUIT_PLUTO(1);
  }
  if (j < 1 || j > NX2_TOT-2){
    print ("! Particles_Interpolate(): j = %d out of bound\n", j);
    print ("  p_intStage = %d\n",p_intStage);
    QUIT_PLUTO(1);
  }
#if DIMENSIONS == 3  
  if (k < 1 || k > NX3_TOT-2){
    print ("! Particles_Interpolate(): k = %d out of bound\n", k);
    print ("  p_intStage = %d\n",p_intStage);
    QUIT_PLUTO(1);  
  }
#endif
*/
/* ----------------------------------------------------
    Interpolate 
   ---------------------------------------------------- */
                 
  for (k1 = -KOFFSET; k1 <= KOFFSET; k1++){
  for (j1 = -JOFFSET; j1 <= JOFFSET; j1++){
  for (i1 = -IOFFSET; i1 <= IOFFSET; i1++){
    v += w[k1][j1][i1]*V[k+k1][j+j1][i+i1];
  }}}

  return v; 
}

/* ********************************************************************* */
int Particles_LocateCell(double *xp, int *indx, Grid *grid) 
/*! 
 * Determine the index of the computational zone hosting the particle,
 * \f$ x_{i-\HALF} <= x_p < x_{i+\HALF} \f$.
 *    
 * This function works for both uniform and stretched grids and employs
 * a binary search algorithm.
 * 
 * \param [in]  xp    a 3-element array specifying the particle position
 * \param [out] indx  a 3-element array giving the cell coordinate (i,j,k)
 * \param [in]  grid  a pointer to an array of grid structures.
 *
 * Return 0 on success, 1 on error.
 *********************************************************************** */
{
  int dir, ngh;
  int l_ind, r_ind, m_ind;  
  double xL;
   
  indx[IDIR] = indx[JDIR] = indx[KDIR] = 0;


#if UNIFORM_CARTESIAN_GRID == YES

  /* Fast search for uniform grid */

  ngh = GetNghost();
  DIM_LOOP(dir) {
    xL        = grid->xbeg[dir] - grid->dx[dir][0]*ngh;
    indx[dir] = (int)((xp[dir] - xL)*grid->inv_dx[dir][IBEG]);
    
    if (indx[dir] < 0 || indx[dir] > grid->np_tot[dir]){
      return 1;
    }  
        
    if (xp[dir] <  grid->xl[dir][indx[dir]]) indx[dir]--; /* Prevent round-off */
    if (xp[dir] >= grid->xr[dir][indx[dir]]) indx[dir]++; /* Prevent round-off */
  }

#else

/* Binary search (any grid)  */

  DIM_LOOP(dir) {
    l_ind = 0;
    r_ind = grid->lend[dir] + grid->nghost[dir];
    while (l_ind != r_ind) {
      m_ind = l_ind + (r_ind - l_ind) / 2;
   
      if (xp[dir] < grid->xr[dir][m_ind]) {
        r_ind = m_ind;
      } else {
        l_ind = m_ind + 1;
      }
    }   
    indx[dir] = r_ind;

  }
#endif

  return 0;
////////////////////////////////////////////////////////
// Check that the particle lies indeed  cell i,j,k
{
  int i = indx[IDIR];
  int j = indx[JDIR];
  int k = indx[KDIR];
  int check = TRUE;
  double xl[3], xr[3];
  
  
  for (dir = 0; dir < DIMENSIONS; dir++) {
    xl[dir] = grid->xl[dir][indx[dir]];
    xr[dir] = grid->xr[dir][indx[dir]];
//    check *= ( (xp[dir] > xr[dir] || xp[dir] < xl[dir]) ? FALSE:TRUE);
    check *= ( (xp[dir] >= xl[dir] && xp[dir] < xr[dir]) ? TRUE:FALSE);
  }  
  
  if ( !check ){
    print ("\n");
    print ("! Particles_LocateCell(): xp = [%12.6e, %12.6e, %12.6e]\n", 
           xp[IDIR],xp[JDIR],xp[KDIR]);
    print ("! cell                       = [%12.6e, %12.6e], [%12.6e, %12.6e], [%12.6e, %12.6e]\n",
              grid->xl[IDIR][i], grid->xr[IDIR][i],    
              grid->xl[JDIR][j], grid->xr[JDIR][j],    
              grid->xl[KDIR][k], grid->xr[KDIR][k]);
    print ("diff(x) =  %12.6e, %12.6e\n",xp[IDIR]-xl[IDIR],-xp[IDIR]+xr[IDIR]);
    print ("diff(y) =  %12.6e, %12.6e\n",xp[JDIR]-xl[JDIR],-xp[JDIR]+xr[JDIR]);

    print ("! indx = %d %d %d\n",indx[IDIR], indx[JDIR], indx[KDIR]);
    print ("! grid extent = [%d, %d]  [%d, %d]  [%d, %d]\n",IBEG,IEND,JBEG,JEND,KBEG,KEND); 
    print ("  outside cell [%d, %d, %d]; stage = %d\n",i,j,k,p_intStage);           
    QUIT_PLUTO(1);
  } 
  
}
////////////////////////////////////////////////////////
  
}


/* ********************************************************************* */
Particle *Particles_Select(particleNode *PHead, int id)
/*!
 *  Loop over particle and return the one with specified id.
 *
 * \param [in] PHead       pointer to the head node of the particle
 *                         linked list.
 * \param [in] id          the particle id
 *
 * Return the particle (if any), or NULL
 *********************************************************************** */
{
  particleNode *curNode;
  Particle     *p;

  PARTICLES_LOOP(curNode, PHead){
    p = &(curNode->p);
    if (p->id == id)  return p;
  }
  return NULL;
}

