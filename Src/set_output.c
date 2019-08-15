/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Set/retrieve output data attributes.
  
  The function SetOutput() sets, for each output data type (DBL, FLT, 
  VTK etc..) the default attributes of the corresponding ::Output structures.
  These include the variable name, a pointer to the actual 
  3D array, the centering of the variable (center/staggered), a 
  conditional inclusion flag (telling if the corresponding variable has
  to be written in the specified format), and so on.
  
  The function SetOutputVar() can be used to include or exclude a given 
  variable to be written using a particular output format.
  
  The function GetUserVar() returns the memory address to a 
  user-defined 3D array.

  \note Starting with PLUTO 4.1 velocity and magnetic field components 
        will be saved as scalars when writing VTK output. 
        If this is not what you want and prefer to save them as vector 
        fields (VTK VECTOR attribute), set VTK_VECTOR_DUMP to YES
        in your definitions.h.        
  
  \authors A. Mignone (mignone@ph.unito.it)
  \date    March 29, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#ifndef VTK_VECTOR_DUMP  
 #define VTK_VECTOR_DUMP NO
#endif

static Output *all_outputs;
/* ********************************************************************* */
void SetOutput (Data *d, Runtime *runtime)
/*!
 *  Set default attributes (variable names, pointers to data structures, 
 *  filename extensions, etc...) of the output structures.
 *
 * \param [in] d        pointer to Data structure
 * \param [in] runtime  pointer to Runtime structure
 *
 *********************************************************************** */
{
  int nv, i, k;
  Output *output;

  if (runtime->user_var > 0)
    d->Vuser = ARRAY_4D(runtime->user_var, NX3_TOT, NX2_TOT, NX1_TOT, double);
  else
    d->Vuser = NULL;

  all_outputs = runtime->output;

/* ---------------------------------------------
          Loop on output types 
   --------------------------------------------- */

  for (k = 0; k < MAX_OUTPUT_TYPES; k++){
    
    output = runtime->output + k;
    output->var_name = ARRAY_2D(MAX_OUTPUT_VARS, 128,char);
//    output->stag_var = ARRAY_1D(64, int);
//    output->dump_var = ARRAY_1D(64, int);
    strcpy(output->dir, runtime->output_dir); /* output directory is the same     */
                                              /* for all outputs (easy to change) */
    output->nfile    = -1;  

  /* -- Set particles filename extensions -- */

/* -- Exclude particle outputs -- */

    #ifdef PARTICLES
    if (output->type == PARTICLES_DBL_OUTPUT ||
        output->type == PARTICLES_FLT_OUTPUT ||
        output->type == PARTICLES_VTK_OUTPUT ||
        output->type == PARTICLES_TAB_OUTPUT ||
        output->type == PARTICLES_HDF5_OUTPUT) continue;
    #endif

  /* -- set variables names -- */

    SetDefaultVarNames(output);

  /* -- Set array pointers -- */

    for (nv = 0; nv < NVAR; nv++){
      output->V[nv]        = d->Vc[nv];
      output->stag_var[nv] = -1; /* -- means cell centered -- */ 
    }
    nv = NVAR;
    #ifdef STAGGERED_MHD
     D_EXPAND(
       output->var_name[nv]   = "Bx1s"; 
       output->V[nv]          = d->Vs[BX1s]; 
       output->stag_var[nv++] = 0;           ,

       output->var_name[nv]   = "Bx2s"; 
       output->V[nv]          = d->Vs[BX2s]; 
       output->stag_var[nv++] = 1;           ,

       output->var_name[nv]   = "Bx3s"; 
       output->V[nv]          = d->Vs[BX3s]; 
       output->stag_var[nv++] = 2; 
     )
    #endif

    #if UPDATE_VECTOR_POTENTIAL == YES
     #if DIMENSIONS == 3   
      output->var_name[nv]   = "Ax1";
      output->V[nv]          = d->Ax1;
      output->stag_var[nv++] = -1;  /* -- vector potential is 
                                          computed at cell center -- */

      output->var_name[nv]   = "Ax2";
      output->V[nv]          = d->Ax2;
      output->stag_var[nv++] = -1;  
     #endif
     output->var_name[nv]   = "Ax3";
     output->V[nv]          = d->Ax3;
     output->stag_var[nv++] = -1;  
    #endif
    output->nvar = nv;

  /* -- repeat for user defined vars -- */

    for (i = 0; i < runtime->user_var; i++){
      sprintf (output->var_name[i + nv], "%s", runtime->user_var_name[i]);
      output->V[i + nv] = d->Vuser[i];
      output->stag_var[i + nv] = -1; /* -- assume cell-centered -- */
    }

  /* -- add user vars to total number of variables -- */

    output->nvar += runtime->user_var;
 
  /* -- select which variables are going to be dumped to disk  -- */

    for (nv = output->nvar; nv--; ) output->dump_var[nv] = YES;
    #if ENTROPY_SWITCH
     output->dump_var[ENTR] = NO;
    #endif

    switch (output->type){
      case DBL_OUTPUT:   /* -- dump ALL variables -- */
        sprintf (output->ext,"dbl");
        break;
      case FLT_OUTPUT:   /* -- do not dump staggered fields (below)-- */
        sprintf (output->ext,"flt");
        break;
      case DBL_H5_OUTPUT:   /* -- dump ALL variables -- */
        sprintf (output->ext,"dbl.h5");
        break;
      case FLT_H5_OUTPUT:   /* -- do not dump staggered fields (below)-- */
        sprintf (output->ext,"flt.h5");
        break;
      case VTK_OUTPUT:   /* -- do not dump staggered fields (below) -- */
        sprintf (output->ext,"vtk");
        #if VTK_VECTOR_DUMP == YES
         D_EXPAND(output->dump_var[VX1] = VTK_VECTOR;  ,
                  output->dump_var[VX2] = NO;          ,
                  output->dump_var[VX3] = NO;)
         #if PHYSICS == MHD || PHYSICS == RMHD
          D_EXPAND(output->dump_var[BX1] = VTK_VECTOR;  ,
                   output->dump_var[BX2] = NO;          ,
                   output->dump_var[BX3] = NO;)
         #endif
        #endif
        break;
      case TAB_OUTPUT:   /* -- do not dump staggered fields -- */
        sprintf (output->ext,"tab");
        break;
      case PPM_OUTPUT:   /* -- dump density only  -- */
        sprintf (output->ext,"ppm");
        for (nv = output->nvar; nv--; ) output->dump_var[nv] = NO;
        break;
      case PNG_OUTPUT:   /* -- dump density only  -- */
        sprintf (output->ext,"png");
        for (nv = output->nvar; nv--; ) output->dump_var[nv] = NO;
        break;
    }
    
  /* ---------------------------------------------------------------
      for divergence cleaning never dump the scalar psi unless
      the output type can be potentially used for restart
     --------------------------------------------------------------- */
   
    #ifdef GLM_MHD
     if (output->type == DBL_OUTPUT || output->type == DBL_H5_OUTPUT) {
       output->dump_var[PSI_GLM] = YES;
       #ifdef PHI_GLM
       output->dump_var[PHI_GLM] = YES;
       #endif       
     } else {
       output->dump_var[PSI_GLM] = NO;
       #ifdef PHI_GLM
       output->dump_var[PHI_GLM] = NO;
       #endif       
     }
    #endif
  }
 
/* -- Exclude staggered components from all output except .dbl and .h5.dbl -- */

  #ifdef STAGGERED_MHD
  D_EXPAND( SetOutputVar ("Bx1s", VTK_OUTPUT, NO);  ,
            SetOutputVar ("Bx2s", VTK_OUTPUT, NO);  ,
            SetOutputVar ("Bx3s", VTK_OUTPUT, NO);)
  D_EXPAND( SetOutputVar ("Bx1s", FLT_OUTPUT, NO);  ,
            SetOutputVar ("Bx2s", FLT_OUTPUT, NO);  ,
            SetOutputVar ("Bx3s", FLT_OUTPUT, NO);)
  D_EXPAND( SetOutputVar ("Bx1s", FLT_H5_OUTPUT, NO);  ,
            SetOutputVar ("Bx2s", FLT_H5_OUTPUT, NO);  ,
            SetOutputVar ("Bx3s", FLT_H5_OUTPUT, NO);)
  D_EXPAND( SetOutputVar ("Bx1s", TAB_OUTPUT, NO);  ,
            SetOutputVar ("Bx2s", TAB_OUTPUT, NO);  ,
            SetOutputVar ("Bx3s", TAB_OUTPUT, NO);)
  #endif

/* -- defaults: dump density only in ppm and png formats -- */

  SetOutputVar ("rho", PPM_OUTPUT, YES);
  SetOutputVar ("rho", PNG_OUTPUT, YES);
  
}

/* ********************************************************************* */
int SetOutputVar (char *var_name, int output_type, int flag)
/*!
 *  Include ('flag == YES') or exclude ('flag == NO') the variable
 *  corresponding to 'var_name' in or from the output type 'out_type'. 
 *  If 'out_type' corresponds to an image (ppm or png), create a
 *  correspdonding Image structure.
 *
 * \param [in] var_name     the name of the variable (e.g. "rho", "vx1",...)
 * \param [in] output_type  select the output type (e.g., DBL_OUTPUT, 
 *                          VTK_OUTPUT, and so forth)
 * \param [in] flag         an integer values (YES/NO).      
 *********************************************************************** */
{
  int k, nv;
  Output *output;

  for (k = 0; k < MAX_OUTPUT_TYPES; k++){ 
    output = all_outputs + k;
    if (output->type == output_type) break;
  }

  for (nv = output->nvar; nv--; ) { 
    if (strcmp(output->var_name[nv], var_name) == 0) {
      output->dump_var[nv] = flag;
      if (flag == YES){
        if (output_type == PPM_OUTPUT || output_type == PNG_OUTPUT){
          CreateImage (var_name);
        }
      }
      return(0);
    }
  }

  print ("! var_name '%s' cannot be set/unset for writing\n",var_name);
  return(1);
}

/* ********************************************************************* */
int GetOutputVarNames(int output_type, char *var_name[NVAR])
/*!
 *  Return the number and the names of the variables being written to
 *  disk with a particular output type.
 *
 * \param [in]  output_type    The output type (e.g. FLT_OUTPUT)
 * \param [out] var_names      An array of strings containing the actual
 *                              variable names.
 *
 * \return An integer giving the number of variables.
 *
 *********************************************************************** */
{
  int k, nv, count = 0;
  Output *output;

  for (k = 0; k < MAX_OUTPUT_TYPES; k++){ 
    output = all_outputs + k;
    if (output->type == output_type) break;
  }

  for (nv = 0; nv < output->nvar; nv++) {
    if (output->dump_var[nv]) sprintf (var_name[count++], "%s", output->var_name[nv]);
  }
  return count;
}


/* ********************************************************************* */
double ***GetUserVar (char *var_name)
/*! 
 *  return a pointer to the 3D array associated with the 
 *  variable named 'var_name'.
 *
 *********************************************************************** */
{
  int indx = -1;
  
  while (strcmp(all_outputs->var_name[++indx], var_name)){
    if (all_outputs->V[indx] == NULL){
      print ("! Error: uservar '%s' is not allocated\n"); 
      QUIT_PLUTO(1);
    }
  }
  return (all_outputs->V[indx]);
}
