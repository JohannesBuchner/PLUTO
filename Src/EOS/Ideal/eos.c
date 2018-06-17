/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file
  \brief Ideal (constant gamma) EOS.
  
  Implements essential functions for the constant-gamma law for
  hydro/MHD, relativistic hydro/relativistic MHD.
                    
  \author A. Mignone (mignone@ph.unito.it)
  \date   Oct 12, 2016
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void SoundSpeed2 (const State *p, int beg, int end, int pos, Grid *grid)
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

  #if PHYSICS == HD || PHYSICS == MHD
  for (i = beg; i <= end; i++) p->a2[i] = g_gamma*p->v[i][PRS]/p->v[i][RHO];
  #elif PHYSICS == RHD || PHYSICS == RMHD
  double theta;
  Enthalpy(p->v, p->h, beg, end);
  for (i = beg; i <= end; i++) {
    theta = p->v[i][PRS]/p->v[i][RHO];
    p->a2[i] = g_gamma*theta/p->h[i];
  }
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

  gmmr = g_gamma/(g_gamma - 1.0);

/* ---------------------------------------------------------------
              Classical equations of state
   --------------------------------------------------------------- */

  #if PHYSICS == HD || PHYSICS == MHD
   for (i = beg; i <= end; i++) h[i] = gmmr*v[i][PRS]/v[i][RHO];
  #elif PHYSICS == RHD || PHYSICS == RMHD
   for (i = beg; i <= end; i++) {
     theta = v[i][PRS]/v[i][RHO];
     h[i] = 1.0 + gmmr*theta;
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

  #if PHYSICS == HD || PHYSICS == MHD
   for (i = beg; i <= end; i++){
     rho  = v[i][RHO];
     s[i] = v[i][PRS]/pow(rho,g_gamma);
   }
  #elif PHYSICS == RHD || PHYSICS == RMHD
   for (i = beg; i <= end; i++) {
     rho = v[i][RHO];
     s[i] = v[i][PRS]/pow(rho,g_gamma);
   }
  #endif
}
