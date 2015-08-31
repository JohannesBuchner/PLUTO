/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Hydrodynamics blast wave problem.

  Set up an isothermal blast wave in a non-uniform, randomnly
  perturbed density field.
  The background perturbation is interpolated read from an input 
  file generated externally.

  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date   May 06, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 *
 *********************************************************************** */
{
  static int first_call = 1;
  double r;

  #if GEOMETRY == CARTESIAN
   r = D_EXPAND(x1*x1, + x2*x2, + x3*x3);
   r = sqrt(r);
  #elif GEOMETRY == POLAR
   r = x1;
  #endif

  v[RHO] = 1.0;

  #if ADD_TURBULENCE == YES
  if (first_call){
    int k, input_var[256];
    
    input_var[0] = RHO;
    input_var[1] = -1;
    InputDataSet ("./grid0.out",input_var);
    InputDataRead("./rho0.dbl"," ");
    first_call = 0;
  }
  InputDataInterpolate(v, x1, x2, x3);  /* -- interpolate density from 
                                              input data file -- */
  #endif

  #if EOS == IDEAL
   v[PRS] = 1.0/g_gamma;  /* ambient pressure */
  #elif EOS == PVTE_LAW
   v[PRS] = 1.0;
  #endif

/* -- set a high pressure region inside r = 1 -- */

  if (r < 0.1) {
    v[RHO] = g_inputParam[RHO_IN];
    #if HAVE_ENERGY
     v[PRS] = g_inputParam[PRS_IN];
    #endif
  }

  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

}
/* **************************************************************** */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 * PURPOSE
 *  
 *   Perform some pre-processing data
 *
 * ARGUMENTS
 *
 *   d:      the PLUTO Data structure.
 *   grid:   pointer to array of GRID structures  
 *
 **************************************************************** */
{

}

/* ************************************************************** */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 *
 **************************************************************** */
{
  int   i, j, k, nv;
  real  *x1, *x2, *x3;

  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){};
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    X1_BEG_LOOP(k,j,i){}
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    X1_END_LOOP(k,j,i){}
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    X2_BEG_LOOP(k,j,i){}
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    X2_END_LOOP(k,j,i){}
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    X3_BEG_LOOP(k,j,i){}
  }
 
  if (side == X3_END) {  /* -- X3_END boundary -- */
    X3_END_LOOP(k,j,i){}
  }
}
