/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute background magnetic field at the desired location.

  This function is normally called from Riemann solvers to compute the
  background magnetic field into the *state structure.

  \authors A. Mignone (mignone@ph.unito.it)\n


  \date   March 05, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if BACKGROUND_FIELD == YES
/* ********************************************************************* */
void GetBackgroundField (const State *state, int beg, int end, int where, Grid *grid)
/*!
 *
 * \param [out] *state   pointer to state structure.
 * \param [in]  beg      starting index of computation.
 * \param [in]  end      final  index of computation.
 * \param [in]  where    the grid location (CELL_CENTER/FACE_CENTER)
 * \param [in]  *grid    pointer to Grid structure.
 *
 *********************************************************************** */
{
  int i;
  double *x1, *x2, *x3;

/* ----------------------------------------------------
         Check for incompatibilities 
   ---------------------------------------------------- */
  
#if (TIME_STEPPING != RK2) && (TIME_STEPPING != RK3)
  print ("! Background field splitting works with RK integrators ONLY \n");
  QUIT_PLUTO(1);
#endif
#if DIVB_CONTROL == EIGHT_WAVES
  print ("! Background field splitting works with CT or GLM ONLY \n");
  QUIT_PLUTO(1);
#elif DIVB_CONTROL == CONSTRAINED_TRANSPORT
  #if (CT_EMF_AVERAGE != ARITHMETIC) && (CT_EMF_AVERAGE != UCT_HLL)
  print ("! Background field splitting works with ARITHMETIC or");
  print (" UCT_HLL averages only\n");
  QUIT_PLUTO(1);
  #endif
#endif

/* ----------------------------------
     get pointers to coordinates 
   ---------------------------------- */
   
  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (g_dir == IDIR){
    if (where == FACE_CENTER) x1 = grid->xr[IDIR];
    for (i = beg; i <= end; i++){
      BackgroundField (x1[i],x2[g_j],x3[g_k], state->Bbck[i]);
    }
  }else if (g_dir == JDIR){
    if (where == FACE_CENTER) x2 = grid->xr[JDIR];
    for (i = beg; i <= end; i++){
      BackgroundField (x1[g_i],x2[i],x3[g_k], state->Bbck[i]);
    }
  }else{
    if (where == FACE_CENTER) x3 = grid->xr[KDIR];
    for (i = beg; i <= end; i++){
      BackgroundField (x1[g_i],x2[g_j],x3[i], state->Bbck[i]);
    }
  }    

}
#endif
