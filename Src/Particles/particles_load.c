/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Initial particle loading.

  \authors A. Mignone (mignone@ph.unito.it)\n

  \b References
     - "Title" \n
       Authors, Journal (year) vol, page

  \date  May 18, 2016
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Particles_LoadRandom(double *qbeg, double *qend,
                          double (*DistribFunc)(double, double, double),
                          double *q)
/*!
 *  Compute particle coordinate (space or velocity) assuming a
 *  a given distribution function using the acceptance-rejection method. 
 *
 *  \param [in]  qbeg         an array specifying the coordinate lower bound 
 *  \param [in]  qend         an array specifying the coordinate upper bound 
 *  \param [in]  DistribFunc  Particle Distribution function.
 *  \param [out] q            the output coordinate
 *
 *********************************************************************** */
{
  int    dir;
  double rnd, prb, dq[3];

  for (dir = 0; dir < 3; dir++) dq[dir] = qend[dir] - qbeg[dir];

/* ------------------------------------------------------------------
    Generate random coordinates q[] and probability rnd.
    If rnd < DistribFunc() then accept.

    Note: this function works in parallel inasmuch the same
    seed is employed on all processors.
   ------------------------------------------------------------------ */

  int success = 0;
  while (!success){

    for (dir = 0; dir < 3; dir++) {
      q[dir] = qbeg[dir] + RandomNumber(0,1)*dq[dir];
    }
    rnd = RandomNumber(0,1);

  /* -- Throw a dice and see if we can accept this particle position */
  
    prb = DistribFunc(q[IDIR], q[JDIR], q[KDIR]);
    if (rnd < prb){  /* Accept */
      success = 1;
    } else {
      success = 0;
    }
  }  
}

/* ********************************************************************* */
void Particles_LoadUniform(int i, int ntot, double *xbeg, double *xend,
                           double *coor)
/*! 
 *  Compute particles coordinate assuming a regular spacing 
 *  given by l = 1/(n)^(1/d) where d is the number of spatial 
 *  dimensions.  
 *
 * \param [in]   i      an integer giving the particle counter inside
 *                      the region [xbeg, xend]
 * \param [in]   ntot   the total number of particles inside the region
 *                      [xbeg, xend]
 * \param [in]   xbeg   an array specifying the coordinates of the lower
 *                      bound in the domain to be filled
 * \param [in]   xbeg   an array specifying the coordinates of the upper 
 *                      bound in the domain to be filled
 * \param [out]  coor   an array containing the particle position
 *  
 *********************************************************************** */
{  
  double Lx = xend[IDIR] - xbeg[IDIR];
  double Ly = xend[JDIR] - xbeg[JDIR];
  double Lz = xend[KDIR] - xbeg[KDIR];

#if DIMENSIONS == 1
  double l = Lx/(double)ntot; 
#elif DIMENSIONS == 2  
  double l = sqrt(Lx*Ly/(double)ntot); 
#elif DIMENSIONS == 3
  double l = pow(Lx*Ly*Lz/(double)ntot,1.0/3.0); 
#endif  
   
  double x0 = xbeg[IDIR] + l*0.5;  /* x-coordinate of 1st particles */
  double y0 = xbeg[JDIR] + l*0.5;  /* y-coordinate of 1st particles */
  double z0 = xbeg[KDIR] + l*0.5;  /* z-coordinate of 1st particles */
  double x1, y1,z1;
      
  x1 = x0 + i*l;
  coor[IDIR] = xbeg[IDIR] + fmod(x1 - xbeg[IDIR], Lx);

  y1 = y0 + (int)((x1 - xbeg[IDIR])/Lx)*l;
  coor[JDIR] = xbeg[JDIR] + fmod(y1 - xbeg[JDIR], Ly);

  z1 = z0 + (int)((y1 - xbeg[JDIR])/Ly)*l;
  coor[KDIR] = xbeg[KDIR] + fmod(z1 - xbeg[KDIR], Lz);
}

