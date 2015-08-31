#include "pluto.h"

/* ****************************************************************** */
void SetDefaultVarNames(Output *output)
/*
 *
 *  PURPOSE
 *
 *    Set file names for I/O
 *
 *
 ******************************************************************** */
{
  int nv;

/* ----------------------------------------------
    Physics module file names; 
    these pertain to the physics module ONLY
   ---------------------------------------------- */

  output->var_name[RHO] = "rho";
  EXPAND(output->var_name[VX1] = "vx1";  ,
         output->var_name[VX2] = "vx2";  ,
         output->var_name[VX3] = "vx3";)
#if HAVE_ENERGY
  output->var_name[PRS] = "prs";
#endif
  
#if PHYSICS == MHD || PHYSICS == RMHD
  EXPAND(output->var_name[BX1] = "bx1";  ,
         output->var_name[BX2] = "bx2";  ,
         output->var_name[BX3] = "bx3";)
#endif
  
  /* (staggered field names are set in SetOutput) */

#ifdef GLM_MHD
  output->var_name[PSI_GLM] = "psi_glm";
#endif
 
/* ------------------------------------------------
    Dust
   ------------------------------------------------ */

#if DUST == YES
  output->var_name[RHO_D] = "rho_d";
  EXPAND(output->var_name[VX1_D] = "vx1_d";  ,
         output->var_name[VX2_D] = "vx2_d";  ,
         output->var_name[VX3_D] = "vx3_d";)
#endif

/* ------------------------------------------------
                   Tracers 
   ------------------------------------------------ */

  NTRACER_LOOP(nv) sprintf (output->var_name[nv],"tr%d",nv - TRC + 1);

  #if ENTROPY_SWITCH
   sprintf (output->var_name[ENTR],"entropy");
  #endif

/* ------------------------------------------------
               Cooling vars
   ------------------------------------------------ */

#if COOLING == MINEq
  {
    static char *ion_name[] = {"X_HI", "X_HeI", "X_HeII" 
                        C_EXPAND("X_CI","X_CII", "X_CIII", "X_CIV", "X_CV")
                        N_EXPAND("X_NI","X_NII", "X_NIII", "X_NIV", "X_NV")
                        O_EXPAND("X_OI","X_OII", "X_OIII", "X_OIV", "X_OV")
                       Ne_EXPAND("X_NeI","X_NeII", "X_NeIII", "X_NeIV", "X_NeV")
                        S_EXPAND("X_SI","X_SI I", "X_SIII", "X_SIV", "X_SV")
                       Fe_EXPAND("X_FeI", "X_FeII", "X_FeIII")};

    NIONS_LOOP(nv) output->var_name[nv] = ion_name[nv-NFLX];  
  }
#elif COOLING == SNEq

   output->var_name[X_HI] = "X_HI";

#elif COOLING == H2_COOL
  {
    static char *molnames[] = {"X_HI", "X_H2", "X_HII"};
    NIONS_LOOP(nv) output->var_name[nv] = molnames[nv-NFLX];
  }  

#endif

}
