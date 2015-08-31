/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Add source terms before the FARGO advection step.
  
  This function is called prior to the FARGO advection algorithm to add
  source terms.
  At present, we use it only for the energy equation in the
  ShearingBox module:
  \f[ 
     \pd{E'}{t} + w\pd{E'}{y} = (B_xB_y - \rho v_xv'_y)\pd{w}{x}
  \f]
  where \f$w = -q\Omega x\f$.
  The discretization follows the algorithm of [GS10], see
  Eq. (51) and (63) of that paper.

  \b Reference
     - [GS10] "Implementation of the shearing box approximation in Athena",
       Stone & Gardiner, ApJS (2010) 189, 142.
 
  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 26, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void FARGO_Source(Data_Arr UU, double dt, Grid *grid)
/*!
 *
 * \param [in,out] UU     an array of conserved variables
 * \param [in]     dt     the current time increment
 * \param [in]     grid   pointer to an array of Grid structures 
 *********************************************************************** */
{
#if (HAVE_ENERGY) && (defined SHEARINGBOX)
  int    i,j,k;
  double scrh, rho, mx, my, Bx, By;

  scrh = SB_Q*SB_OMEGA;
  DOM_LOOP(k,j,i){
    rho = UU[k][j][i][RHO];
    mx  = UU[k][j][i][MX1];
    my  = UU[k][j][i][MX2];
    Bx  = UU[k][j][i][BX1];
    By  = UU[k][j][i][BX2];

    UU[k][j][i][ENG] += - dt*scrh*Bx*(By - 0.5*dt*Bx*scrh)
                        + dt*scrh*mx*my/rho;
  }
#endif
}
