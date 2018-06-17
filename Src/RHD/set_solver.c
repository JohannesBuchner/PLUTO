#include"pluto.h"

/* ****************************************************************** */
Riemann_Solver *SetSolver (const char *solver)
/*
 *
 * PURPOSE
 *
 *   return a pointer to a riemann solver function
 *
 *
 *
 *
 ********************************************************************* */
{
  
/* ------------------------------------------------------
       Set Pointers for SOLVERS 
   ------------------------------------------------------ */

  if ( !strcmp(solver, "two_shock"))  {
  
    return(&TwoShock_Solver);

  }else if (!strcmp(solver, "tvdlf")) {

    return (&LF_Solver);

  }else if (!strcmp(solver, "hlle") || 
            !strcmp(solver, "hll")) {

    return (&HLL_Solver);

  }else if (!strcmp(solver, "hllc")) {

    return (&HLLC_Solver);

/*
  }else if (!strcmp(solver, "musta")) {
    return (&MUSTA_FLUX);
*/

  }

  print ("\n! SetSolver: '%s' is not available.\n", solver);
  QUIT_PLUTO(1);
    
}
