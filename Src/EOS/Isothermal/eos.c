/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file
  \brief Implementation of the isothermal EOS.
                    
  \author A. Mignone (mignone@ph.unito.it)
  \date   April 14, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void SoundSpeed2 (double **v, double *cs2, double *h, int beg, int end,
                  int pos, Grid *grid)
/*!
 * Define the square of the sound speed.
 * 
 * \param [in]    v   1D array of primitive quantities
 * \param [out] cs2   1D array containing the square of the sound speed
 * \param [in]    h   1D array of enthalpy values
 * \param [in]  beg   initial index of computation 
 * \param [in]  end   final   index of computation
 * \param [in]  pos   an integer specifying the spatial position 
 *                    inside the cell (only for spatially-dependent EOS)
 * \param [in]  grid  pointer to an array of Grid structures
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int  i, j,k;  /* -- used as multidimensional indices -- */
  double *x1, *x2, *x3;

  #if PHYSICS == HD || PHYSICS == MHD
   x1 = grid[IDIR].x;
   x2 = grid[JDIR].x;
   x3 = grid[KDIR].x;

   i = g_i; j = g_j; k = g_k;
    
   if (g_dir == IDIR) {

     x1 = (pos == FACE_CENTER ? grid[IDIR].xr : grid[IDIR].x);
     for (i = beg; i <= end; i++) cs2[i] = g_isoSoundSpeed*g_isoSoundSpeed;

   }else if (g_dir == JDIR){

     x2 = (pos == FACE_CENTER ? grid[JDIR].xr : grid[JDIR].x);
     for (j = beg; j <= end; j++) cs2[j] = g_isoSoundSpeed*g_isoSoundSpeed;

   }else if (g_dir == KDIR){

     x3 = (pos == FACE_CENTER ? grid[KDIR].xr : grid[KDIR].x);
     for (k = beg; k <= end; k++) cs2[k] = g_isoSoundSpeed*g_isoSoundSpeed;

  }
  #else
   print ("! SoundSpeed2: not defined for this EoS\n");
   QUIT_PLUTO(1);
  #endif
}

/* ********************************************************************* */
void Enthalpy (double **v, real *h, int beg, int end)
/*!
 * Compute the enthalpy.
 *
 * \param [in]    v   1D array of primitive quantities
 * \param [in]    h   1D array of enthalpy values
 * \param [in]  beg   initial index of computation 
 * \param [in]  end   final   index of computation
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  print ("! Enthalpy: enthalpy not defined in isothermal EOS\n");
  QUIT_PLUTO(1);
}

/* ********************************************************************* */
void Entropy (double **v, double *s, int beg, int end)
/*!
 * Compute the entropy.
 * 
 * \param [in]    v   1D array of primitive quantities
 * \param [in]    s   1D array of entropy values
 * \param [in]   is   initial index of computation 
 * \param [in]   ie   final   index of computation
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  print ("Entropy: entropy not defined for isothermal EOS\n");
  QUIT_PLUTO(1);
}
