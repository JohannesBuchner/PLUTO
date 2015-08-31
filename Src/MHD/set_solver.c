/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Return a pointer to a Riemann solver function.

  \author A. Mignone (mignone@ph.unito.it)
  \date   June 5, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/* ********************************************************************* */
Riemann_Solver *SetSolver (const char *solver)
/*!
 *  Depending on the choice of the Riemann solver specified in
 *  pluto.ini, return a pointer to the corresponding Riemann solver
 *  function
 *
 *********************************************************************** */
{

  #ifdef FINITE_DIFFERENCE
   return (&FD_Flux);
  #endif

/* ------------------------------------------------------
       Set Pointers for SOLVERS
   ------------------------------------------------------ */

  #if EOS == IDEAL
  
   if      (!strcmp(solver, "tvdlf"))   return (&LF_Solver);
   else if (!strcmp(solver, "roe"))     return (&Roe_Solver);
   else if (!strcmp(solver, "hlle") ||
            !strcmp(solver, "hll"))     return (&HLL_Solver);
   else if (!strcmp(solver, "hllc"))    return (&HLLC_Solver);
   else if (!strcmp(solver, "hlld"))    return (&HLLD_Solver);
/*
   else if (!strcmp(solver, "musta"))   return (&MUSTA_Solver);
*/

   #elif EOS == PVTE_LAW
   if (!strcmp(solver, "tvdlf"))       return (&LF_Solver);
   else if (!strcmp(solver, "hlle") || 
            !strcmp(solver, "hll"))    return (&HLL_Solver);
   else if (!strcmp(solver, "hllc"))   return (&HLLC_Solver);
   else if (!strcmp(solver, "hlld"))   return (&HLLD_Solver);

  #elif EOS == ISOTHERMAL

   if      (!strcmp(solver, "tvdlf"))   return (&LF_Solver);
   else if (!strcmp(solver, "roe"))     return (&Roe_Solver);
   else if (!strcmp(solver, "hlle") ||
            !strcmp(solver, "hll"))     return (&HLL_Solver);
   else if (!strcmp(solver, "hlld"))    return (&HLLD_Solver);
/*
   else if (!strcmp(solver, "musta"))   return (&MUSTA_Solver);
*/

  #elif EOS == BAROTROPIC

   if (!strcmp(solver, "tvdlf"))        return (&LF_Solver);
   else if (!strcmp(solver, "hlle") ||
             !strcmp(solver, "hll"))    return (&HLL_Solver);

  #endif

  print1 ("\n! SetSolver: '%s' is not available.\n", solver);
  QUIT_PLUTO(1);

}

