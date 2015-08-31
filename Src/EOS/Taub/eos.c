/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file
  \brief Taub-Matthews (TM) EOS for relativistic hydro and MHD.
                    
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
  int  i;
  double theta;

  #if PHYSICS == RHD || PHYSICS == RMHD
   Enthalpy(v, h, beg, end);
   for (i = beg; i <= end; i++) {
     theta = v[i][PRS]/v[i][RHO];
     #if EOS == IDEAL
      cs2[i] = g_gamma*theta/h[i];
     #elif EOS == TAUB
      cs2[i] = theta/(3.0*h[i])*(5.0*h[i] - 8.0*theta)/(h[i] - theta);
     #endif
   }
  #else
   print ("! SoundSpeed2: Taub EOS not defined for this physics module.\n");
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
