/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Write trajectoy of a particle
  
  \authors A. Mignone (mignone@ph.unito.it)\n

  \date   Sep 23, 2016
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Particles_WriteTrajectory (Particle *p, char mode)
/*!
 *  Write particle coordinates to an ASCII file, named
 *  "particle.<br>.<id>.dat" where <br> and <id> are the birth rank
 *  and id of the particle.
 *  The file used a multiple column format,
 *
 *  <t>  <x1>  <x2>  <x3>
 *
 * where <t> is the time column and the other contain coordinates.
 *
 * \param [in]   p     A pointer to a particle structure
 * \param [in]   mode  Either "w" or "a" to write or append to the file.
 * 
 *********************************************************************** */
{
  FILE *fp;
  char fname[64];

  sprintf (fname,"particle.%04d.dat",p->id);
  if      (mode == 'w') fp = fopen (fname,"w");
  else if (mode == 'a') fp = fopen (fname,"a");
  else{
    print ("! ParticlesWriteTrajectory: invalid mode\n");
    QUIT_PLUTO(1);
  }

  fprintf (fp, "%12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n",
           g_time, p->coord[IDIR], p->coord[JDIR], p->coord[KDIR],
                   p->speed[IDIR], p->speed[JDIR], p->speed[KDIR]);
  fclose(fp);
}



