/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Adjust the effective domain size when the -xNjet option is given.

  The SetJetDomain() function shortens the domain integration index in the 
  direction specified by either '\c -x1jet', '\c -x2jet' or '\c -x3jet' to save 
  computational time.
  This is done by setting the final integration index in the direction
  of jet propagation to that of the first zone (counting from the top) with 
  non-zero pressure gradient (GetRightmostIndex()).
  The function UnsetJetDomain() restores the initial computational domain
  size.
  Useful for problems involving jet propagation.

  \note In parallel, the domain is \e not decomposed along the
        propagation direction (see ParseCmdLineArgs()).
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   Feb 28, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static int jd_nbeg, jd_nend, jd_npt, jd_ntot, rbound;
static int GetRightmostIndex(int, double ***);

/* ********************************************************************* */
void SetJetDomain (const Data *d, int dir, int log_freq, Grid *grid)
/*!
 *
 * Adjust the size of computational domain by reducing the final
 * computationa index in the direction of propagation.
 *
 * \param [in]         d   pointer to Data structure
 * \param [in]       dir   the direction of propagation
 * \param [in]  log_freq   the output log frequency
 * \param [in,out]  grid   pointer to Grid structure
 * 
 *********************************************************************** */
{
  int i, j, k, ngh;
  int n, n_glob;
  static int first_call = 1;
  double ***pr, ***dn, dp;

  dn = d->Vc[RHO];
  #if HAVE_ENERGY
   pr = d->Vc[PRS];
  #else
   pr = d->Vc[RHO];
  #endif

  ngh = grid->nghost[dir];

/* -- save original domain offsets 
      and return for the first time -- */

  if (first_call){
    if (dir == IDIR) {
      jd_nbeg = IBEG; jd_nend = IEND; 
      jd_npt  = NX1;  jd_ntot = NX1_TOT;
    }else if (dir == JDIR){
      jd_nbeg = JBEG; jd_nend = JEND; 
      jd_npt  = NX2;  jd_ntot = NX2_TOT;
    }else if (dir == KDIR){
      jd_nbeg = KBEG; jd_nend = KEND; 
      jd_npt  = NX3;  jd_ntot = NX3_TOT;
    }

    rbound = grid->rbound[dir];
    first_call = 0;
    return;
  }

/* --------------------------------------
    find where grad(p) is not zero and 
    add safety guard cells 
   -------------------------------------- */

  n = GetRightmostIndex(dir, pr) + 2*ngh;

  if (n < jd_nbeg + ngh) n = jd_nbeg + ngh;
  if (n > jd_nend)       n = jd_nend;

  #ifdef PARALLEL
   MPI_Allreduce (&n, &n_glob, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   n = n_glob;
  #endif

  if (g_stepNumber%log_freq==0){
/*    print ("- SetJetDomain: index %d / %d\n",n,jd_nend); */
  }

/* ------------------------------------------------------------------
    Change global integer variables giving the rightmost integration 
    indexes. Also, suppress the rightmost physics boundary (rbound) 
    if solution has not yet reached it. 
   ------------------------------------------------------------------ */

  if (dir == IDIR) {
    IEND = grid->lend[IDIR] = n;
    NX1_TOT = IEND + ngh;
    NX1     = IEND - ngh;

    grid->rbound[IDIR] = rbound;
    if (jd_nend != IEND) grid->rbound[IDIR] = 0;
    
  } else if (dir == JDIR){

    JEND = grid->lend[JDIR] = n;
    NX2_TOT = JEND + ngh;
    NX2     = JEND - ngh;
    grid->rbound[JDIR] = rbound;
    if (jd_nend != JEND) grid->rbound[JDIR] = 0;

  }else if (dir == KDIR){

    KEND = grid->lend[KDIR] = n;
    NX3_TOT = KEND + ngh;
    NX3     = KEND - ngh;

    grid->rbound[KDIR] = rbound;
    if (jd_nend != KEND) grid->rbound[KDIR] = 0;  
  }
  
}

/* ********************************************************************* */
void UnsetJetDomain (const Data *d, int dir, Grid *grid)
/*!
 *  Restore original (full domain) indexes.
 *
 *********************************************************************** */
{

  if (dir == IDIR){
    IBEG = grid->lbeg[IDIR] = jd_nbeg;
    IEND = grid->lend[IDIR] = jd_nend;
    NX1_TOT = jd_ntot;
    NX1     = jd_npt;
    grid->rbound[IDIR] = rbound;
  }

  if (dir == JDIR){
    JBEG = grid->lbeg[JDIR] = jd_nbeg;
    JEND = grid->lend[JDIR] = jd_nend;
    NX2_TOT = jd_ntot;
    NX2     = jd_npt;
    grid->rbound[JDIR] = rbound;
  }

  if (dir == KDIR){
    KBEG = grid->lbeg[KDIR] = jd_nbeg;
    KEND = grid->lend[KDIR] = jd_nend;
    NX3_TOT = jd_ntot;
    NX3     = jd_npt;
    grid->rbound[KDIR] = rbound;
  }

}

/* ********************************************************************* */
int GetRightmostIndex (int dir, double ***q)
/*!
 *
 * Find the local index j where grad(p) exceeds
 * a ceratin threshold
 *
 *********************************************************************** */
{
  int    i, j, k;
  double dp;

/* -- backward sweep on x-axis -- */

  if (dir == IDIR){
    for (i = jd_nend; i >= jd_nbeg; i--){
    KDOM_LOOP(k){
    JDOM_LOOP(j){

      dp  = q[k][j][i+1] - q[k][j][i-1];
      dp /= q[k][j][i+1] + q[k][j][i-1];

      if (fabs(dp) > 1.e-5) return(i);
    }}}
  }

/* -- backward sweep on y-axis -- */

  if (dir == JDIR){
    for (j = jd_nend; j >= jd_nbeg; j--){
    KDOM_LOOP(k){
    IDOM_LOOP(i){

      dp  = q[k][j+1][i] - q[k][j-1][i];
      dp /= q[k][j+1][i] + q[k][j-1][i];

      if (fabs(dp) > 1.e-5) return(j);
    }}}
  }

/* -- backward sweep on z-axis -- */

  if (dir == KDIR){
    for (k = jd_nend; k >= jd_nbeg; k--){
    JDOM_LOOP(j){
    IDOM_LOOP(i){

      dp  = q[k+1][j][i] - q[k-1][j][i];
      dp /= q[k+1][j][i] + q[k-1][j][i];

      if (fabs(dp) > 1.e-5) return(k);
    }}}
  }

  return(-1);
}
