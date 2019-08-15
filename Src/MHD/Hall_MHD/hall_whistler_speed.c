/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute the whistler speed for the Hall MHD equations.

  Compute whistler wave speed in and save it to \c state->cw:
  \f[
    c_w =   \left|\frac{B}{2en_e\Delta l}\right|
          + \sqrt{\left(\frac{B}{2en_e\Delta l}\right)^2 + \frac{B^2}{\rho}}
  \f]
  where \f$\Delta l\f$ is the grid spacing.

  \authors E. Striani (edoardo.striani@iaps.inaf.it) \n
           A. Mignone (mignone@ph.unito.it) \n
           B. Vaidya  (bhargav.vaidya@unito.it) 

  \b References
     - "Thanatology in Protoplanetary Discs"
       Lesur, Kunz \& Fromang, A&A (2014) 

  \date   May 10, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void HallMHD_WhistlerSpeed (const State *p, int beg, int end, Grid *grid)
/*!
 * Define the whistler speed defined from the equation in the
 * Appendix of Lesur et al. 2014.
 * 
 * \param [in]   p    pointer to a state structure
 * \param [in]  beg   initial index of computation 
 * \param [in]  end   final   index of computation
 * \param [in]  grid  pointer to an array of Grid structures
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int    i, j, k;
  double *x1  = grid->x[IDIR], *dx1 = grid->dx[IDIR];
  double *x2  = grid->x[JDIR], *dx2 = grid->dx[JDIR];
  double *x3  = grid->x[KDIR], *dx3 = grid->dx[KDIR];
  double *rt  = grid->rt;
  double *dmu = grid->dmu;

  double dl, ne, Bsq, B, cA2, B_ne, hallterm, scrh;
  double *v;

/* ---------------------------------------------------
    Compute the whistler speed direction by direction
    since geometrical factors may enter in polar
    and spherical geometries.
   --------------------------------------------------- */

  i = g_i; 
  j = g_j;
  k = g_k;
  if (g_dir == IDIR){
    for (i = beg; i <= end; i++) {
      dl   = dx1[i];
      v    = p->v[i];
      ne   = HallMHD_ne(v);
      Bsq  = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
      B    = sqrt(Bsq);
      cA2  = Bsq/v[RHO];  /* Alfven velocity */
      B_ne = B/ne;
      hallterm = B_ne/(2.0*dl);
      p->cw[i] = fabs(hallterm) + sqrt(hallterm*hallterm + cA2); 
    }
  }else if (g_dir == JDIR){
    for (j = beg; j <= end; j++) {
      dl   = dx2[j];
      #if GEOMETRY == POLAR
      dl   *= x1[i];
      #elif GEOMETRY == SPHERICAL
      dl   *= rt[i];
      #endif
      v    = p->v[j];
      ne   = HallMHD_ne(v);
      Bsq  = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
      B    = sqrt(Bsq);
      cA2  = Bsq/v[RHO];  /* Alfven velocity */
      B_ne = B/ne;
      hallterm = B_ne/(2.0*dl);
      p->cw[j] = fabs(hallterm) + sqrt(hallterm*hallterm + cA2); 
    }
  }else if (g_dir == KDIR){
    scrh = rt[i]*dmu[j]/dx2[j];  /* Multiplicative factor to get length */
    for (k = beg; k <= end; k++) {
      dl   = dx3[k];
      #if GEOMETRY == SPHERICAL
      dl  *= scrh;
      #endif
      v    = p->v[k];
      ne   = HallMHD_ne(v);
      Bsq  = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
      B    = sqrt(Bsq);
      cA2  = Bsq/v[RHO];  /* Alfven velocity */
      B_ne = B/ne;
      hallterm = B_ne/(2.0*dl);
      p->cw[k] = fabs(hallterm) + sqrt(hallterm*hallterm + cA2); 
    }
  }



for (i = beg; i <= end; i++) {
//  print ("dir = %d, i = %d, cw = %12.6e\n",g_dir,i,p->cw[i]);
}


}
