/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file
  \brief Taub-Matthews (TM) EOS for relativistic hydro and MHD.
                    
  \author A. Mignone (mignone@ph.unito.it)
  \date   Oct 13, 2016
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void SoundSpeed2 (const State *q, int beg, int end, int pos, Grid *grid)
/*!
 * Define the square of the sound speed.
 * 
 * \param [in]   p    pointer to a state structure
 * \param [in]  beg   initial index of computation 
 * \param [in]  end   final   index of computation
 * \param [in]  pos   an integer specifying the spatial position 
 *                    inside the cell (only for spatially-dependent EOS)
 * \param [in]  grid  pointer to an array of Grid structures
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int  i;
  double theta;

  #if PHYSICS == RHD || PHYSICS == RMHD
   Enthalpy(q->v, q->h, beg, end);
   for (i = beg; i <= end; i++) {
     theta = q->v[i][PRS]/q->v[i][RHO];
     #if EOS == IDEAL
     q->a2[i] = g_gamma*theta/q->h[i];
     #elif EOS == TAUB
     q->a2[i] = theta/(3.0*q->h[i])*(5.0*q->h[i] - 8.0*theta)/(q->h[i] - theta);
     #endif
   }
  #else
   print ("! SoundSpeed2: Taub EOS not defined for this physics module.\n");
   QUIT_PLUTO(1);
  #endif
}

/* ********************************************************************* */
void Enthalpy (double **v, double *h, int beg, int end)
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
  int i;
  double gmmr, theta;

  #if PHYSICS == RHD || PHYSICS == RMHD
   for (i = beg; i <= end; i++) {
     theta = v[i][PRS]/v[i][RHO];
     #if EOS == IDEAL
      h[i] = 1.0 + gmmr*theta;
     #elif EOS == TAUB
      h[i] = 2.5*theta + sqrt(2.25*theta*theta + 1.0);
    #endif
   }
  #endif
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
  int i;
  double rho, th;

  #if PHYSICS == RHD || PHYSICS == RMHD
   for (i = beg; i <= end; i++) {
     rho = v[i][RHO];
     #if EOS == IDEAL
      s[i] = v[i][PRS]/pow(rho,g_gamma);
     #elif EOS == TAUB
      th   = v[i][PRS]/rho;
      s[i] = v[i][PRS]/pow(rho,5./3.)*(1.5*th + sqrt(2.25*th*th + 1.0));
     #endif
   }
  #endif
}
