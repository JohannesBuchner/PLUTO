/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Read runtime information from pluto.ini.

  Parse and read runtime information from the initialization file 
  pluto.ini (default) and sets value of the runtime structure.

  \authors A. Mignone (mignone@ph.unito.it)
  \date    Aug 24, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define NOPT      32              /*  # of possible options in a menu */
#define NLEN      128              /*  # default string length         */

#define COMPARE(s1,s2,ii) \
        for ( (ii) = 1 ; (ii) < NOPT && !(strcmp ( (s1), (s2)) == 0); (ii)++);

/* ********************************************************************* */
int RuntimeSetup (Runtime *runtime, Cmd_Line *cmd_line, char *ini_file)
/*!
 * Open and parse the runtime initialization file. 
 * Assign values to the runtime structure.
 *
 * \param [out]  runtime     pointer to a Runtime structure
 * \param [in]   cmd_line  pointer to a Cmd_Line structure (useful, e.g.,
 *                         to resize the domain using the \c -xres option)
 * \param [in]   ini_file  the name of the initialization file (default
 *                         is "pluto.ini") specified with the \c -i option.
 *
 *********************************************************************** */
{
  int    idim, ip, ipos, itype, nlines;
  char   *bound_opt[NOPT], str_var[512], *str;
  char  *glabel[]     = {"X1-grid", "X2-grid","X3-grid"};
  char  *bbeg_label[] = {"X1-beg", "X2-beg","X3-beg"};
  char  *bend_label[] = {"X1-end", "X2-end","X3-end"};
  double dbl_var, rx;
  Output *output;
  FILE *fp;

  for (itype = 0; itype < NOPT; itype++) {
    bound_opt[itype] = "0000";
  }

/*  ---------------------------------------------------
      available options are given as two set of names;
      This facilitates when updating the code and
      people are too lazy to read the manual !
    --------------------------------------------------- */

  bound_opt[OUTFLOW]      = "outflow";
  bound_opt[REFLECTIVE]   = "reflective";
  bound_opt[AXISYMMETRIC] = "axisymmetric";
  bound_opt[EQTSYMMETRIC] = "eqtsymmetric";
  bound_opt[PERIODIC]     = "periodic";
  bound_opt[SHEARING]     = "shearingbox";
  bound_opt[USERDEF]      = "userdef";

  runtime->log_freq = 1; /* -- default -- */
 
  nlines = ParamFileRead(ini_file);

/* ------------------------------------------------------------
   [Grid] Section 
   ------------------------------------------------------------ */

  for (idim = 0; idim < 3; idim++){
    runtime->npatch[idim] = atoi(ParamFileGet(glabel[idim], 1));
    runtime->npoint[idim] = 0;

    ipos = 1;
    for (ip = 1; ip <= runtime->npatch[idim]; ip++) {

      runtime->patch_left_node[idim][ip] = atof(ParamFileGet(glabel[idim], ++ipos));
      runtime->patch_npoint[idim][ip]    = atoi(ParamFileGet(glabel[idim], ++ipos));
      runtime->npoint[idim]             += runtime->patch_npoint[idim][ip];
      runtime->grid_is_uniform[idim]     = 0;

      strcpy (str_var, ParamFileGet(glabel[idim], ++ipos)); 
/*
printf ("%f  %d %s\n",runtime->patch_left_node[idim][ip],runtime->patch_npoint[idim][ip],str_var);
*/
      if (strcmp(str_var,"u") == 0 || strcmp(str_var,"uniform") == 0) {
        runtime->patch_type[idim][ip] = UNIFORM_GRID;
        if (runtime->npatch[idim] == 1) runtime->grid_is_uniform[idim] = 1;        
      }else if (strcmp(str_var,"s") == 0 || strcmp(str_var,"strecthed") == 0) { 
        runtime->patch_type[idim][ip] = STRETCHED_GRID;
      }else if (strcmp(str_var,"l+") == 0){
        runtime->patch_type[idim][ip] = LOGARITHMIC_INC_GRID;
      }else if (strcmp(str_var,"l-") == 0){
        runtime->patch_type[idim][ip] = LOGARITHMIC_DEC_GRID;
      }else{ 
        printf ("\nSetup: You must specify either 'u', 's', 'l+' or 'l-' as grid-type in %s\n",
                ini_file);
        QUIT_PLUTO(1);
      }
    }
    
    runtime->patch_left_node[idim][ip] = atof(ParamFileGet(glabel[idim], ++ipos));

    if ( (ipos+1) != (runtime->npatch[idim]*3 + 3)) {
      printf ("! Setup: domain #%d setup is not properly defined \n", idim);
      QUIT_PLUTO(1);
    }
    if (idim >= DIMENSIONS && runtime->npoint[idim] != 1) {
      printf ("! Setup: %d point(s) on dim. %d is NOT valid, resetting to 1\n",
              runtime->npoint[idim],idim+1);
      runtime->npoint[idim]          = 1;
      runtime->npatch[idim]          = 1;
      runtime->patch_npoint[idim][1] = 1;
    }
  }

/* ------------------------------------------------------------
      Change the resolution if cmd_line->xres has been given
   ------------------------------------------------------------ */

  if (cmd_line->xres > 1) {
    rx =  (double)cmd_line->xres/(double)runtime->patch_npoint[IDIR][1];
    for (idim = 0; idim < DIMENSIONS; idim++){
      if (runtime->npatch[idim] > 1){  
        printf ("! Setup: -xres option works on uniform, single patch grid\n");
        QUIT_PLUTO(1);
      }
      
      dbl_var = (double)runtime->patch_npoint[idim][1];
      runtime->patch_npoint[idim][1] = MAX( (int)(dbl_var*rx), 1);
      dbl_var = (double)runtime->npoint[idim];
      runtime->npoint[idim] = MAX( (int)(dbl_var*rx), 1); 
    }  
  }

/* ------------------------------------------------------------
   [Time] Section 
   ------------------------------------------------------------ */

  runtime->cfl         = atof(ParamFileGet("CFL", 1));

  if (ParamExist ("CFL_par")) runtime->cfl_par = atof(ParamFileGet("CFL_par", 1));
  else                        runtime->cfl_par = 0.8/(double)DIMENSIONS;

  if (ParamExist ("rmax_par")) runtime->rmax_par = atof(ParamFileGet("rmax_par", 1));
  else                         runtime->rmax_par = 100.0;

  runtime->cfl_max_var = atof(ParamFileGet("CFL_max_var", 1));
  runtime->tstop       = atof(ParamFileGet("tstop", 1));
  runtime->first_dt    = atof(ParamFileGet("first_dt", 1));

/* ------------------------------------------------------------
   [Solver] Section 
   ------------------------------------------------------------ */

  sprintf (runtime->solv_type,"%s",ParamFileGet("Solver",1));

/* ------------------------------------------------------------
   [Boundary] Section 
   ------------------------------------------------------------ */

  for (idim = 0; idim < 3; idim++){

    str = ParamFileGet(bbeg_label[idim], 1);
    COMPARE (str, bound_opt[itype], itype);
    if (itype == NOPT) {
      printf ("! Setup: don't know how to put left boundary '%s'  \n", str);
      QUIT_PLUTO(1);
    }
    runtime->left_bound[idim] = itype;
  }

  for (idim = 0; idim < 3; idim++){

    str = ParamFileGet(bend_label[idim], 1);
    COMPARE (str, bound_opt[itype], itype);
    if (itype == NOPT) {
      printf ("! Setup: don't know how to put left boundary '%s'  \n", str);
      QUIT_PLUTO(1);
    }
    runtime->right_bound[idim] = itype;
  }

/* ------------------------------------------------------------
   [Static Grid Output] Section 
   ------------------------------------------------------------ */

#ifndef CHOMBO
  runtime->user_var = atoi(ParamFileGet("uservar", 1));
  for (ip = 0; ip < runtime->user_var; ip++){
    if ( (str = ParamFileGet("uservar", 2 + ip)) != NULL){
      sprintf (runtime->user_var_name[ip], "%s", str);
    }else{
      printf ("! Setup: missing name after user var name '%s'\n", 
              runtime->user_var_name[ip-1]);
      QUIT_PLUTO(1);
    } 
  }

/* ---- set output directory ---- */

  sprintf (runtime->output_dir, "%s","./");  /* default value is current directory */
  if (ParamExist("output_dir")){
    str = ParamFileGet("output_dir",1);
    sprintf (runtime->output_dir, "%s",str);
  }

/* -- check if we have write access and if the directory exists -- */

  sprintf (str_var,"%s/tmp09123.txt",runtime->output_dir); /* -- test file -- */
  fp = fopen(str_var,"w");  /* -- open test file for writing -- */
  if (fp == NULL){
    printf ("! Setup: cannot access directory '%s'.\n", runtime->output_dir);
    printf ("!        Please check that the directory exists\n");
    printf ("!        and you have write permission.\n");
    QUIT_PLUTO(1);
  }else{
    fclose(fp);
    remove(str_var);  /* -- remove file -- */
  }

/* ---- dbl output ---- */

  ipos = 0;
  output = runtime->output + (ipos++);
  output->type  = DBL_OUTPUT;
  output->cgs   = 0;  /* cannot write .dbl using cgs units */
  GetOutputFrequency(output, "dbl");

  sprintf (output->mode,"%s",ParamFileGet("dbl",3));
  #ifdef USE_ASYNC_IO
   if (    strcmp(output->mode,"single_file") 
        && strcmp(output->mode,"single_file_async")
        && strcmp(output->mode,"multiple_files")){
      printf ("! Setup: expecting 'single_file', 'single_file_async' ");
      printf ("or 'multiple_files' in dbl output\n");
      QUIT_PLUTO(1);
   }
  #else
   if (   strcmp(output->mode,"single_file")
       && strcmp(output->mode,"multiple_files")){
      printf (
      "! Setup: expecting 'single_file' or 'multiple_files' in dbl output\n");
      QUIT_PLUTO(1);
   }     
  #endif

 /* ---- flt output ---- */

  if (ParamExist("flt")){
    output = runtime->output + (ipos++);
    output->type  = FLT_OUTPUT;
    GetOutputFrequency(output, "flt");

    sprintf (output->mode,"%s",ParamFileGet("flt",3));  
    #ifdef USE_ASYNC_IO
     if (    strcmp(output->mode,"single_file") 
          && strcmp(output->mode,"single_file_async")
          && strcmp(output->mode,"multiple_files")){
        printf ("! Setup: expecting 'single_file', 'single_file_async' ");
        printf ("or 'multiple_files' in flt output\n");
        QUIT_PLUTO(1);
     }
    #else
     if (    strcmp(output->mode,"single_file") 
          && strcmp(output->mode,"multiple_files")){
        printf (
        "! Setup: expecting 'single_file' or 'multiple_files' in flt output\n");
        QUIT_PLUTO(1);
     }  
    #endif
    if (ParamFileHasBoth ("flt","cgs")) output->cgs = 1;
    else                                output->cgs = 0;
  }

 /* -- hdf5 output -- */

  if (ParamExist("dbl.h5")){
    output = runtime->output + (ipos++);
    output->type  = DBL_H5_OUTPUT;
    output->cgs   = 0;  /* cannot write .h5 using cgs units */
    GetOutputFrequency(output, "dbl.h5");
  }
  if (ParamExist("flt.h5")){
    output = runtime->output + (ipos++);
    output->type  = FLT_H5_OUTPUT;
    output->cgs   = 0;  /* cannot write .h5 using cgs units */
    GetOutputFrequency(output, "flt.h5");
  }

 /* -- vtk output -- */

  if (ParamExist ("vtk")){
    output = runtime->output + (ipos++);
    output->type  = VTK_OUTPUT;
    GetOutputFrequency(output, "vtk");

    if (ParamFileGet("vtk",3) == NULL){
      printf ("! Setup: extra field missing in vtk output\n");
      QUIT_PLUTO(1);
    }
    sprintf (output->mode,"%s",ParamFileGet("vtk",3));
    if (   strcmp(output->mode,"single_file")
        && strcmp(output->mode,"multiple_files")){
       printf ("! Setup: expecting 'single_file' or 'multiple_files' in\n");
       printf ("         vtk output\n");
       QUIT_PLUTO(1);
    }
    if (ParamFileHasBoth ("vtk","cgs")) output->cgs = 1;
    else                                output->cgs = 0;
  }

 /* -- tab output -- */

  if (ParamExist ("tab")){
    output = runtime->output + (ipos++);
    output->type  = TAB_OUTPUT;
    GetOutputFrequency(output, "tab");
    if (ParamFileHasBoth ("tab","cgs")) output->cgs = 1;
    else                                output->cgs = 0;
  }

 /* -- ppm output -- */

  if (ParamExist ("ppm")){
    output = runtime->output + (ipos++);
    output->type  = PPM_OUTPUT;
    output->cgs   = 0;   /* Cannot write ppm in cgs units */
    GetOutputFrequency(output, "ppm");
  }

 /* -- png output -- */

  if (ParamExist ("png")){
    output = runtime->output + (ipos++);
    output->type  = PNG_OUTPUT;
    output->cgs   = 0;   /* Cannot write png in cgs units */
    GetOutputFrequency(output, "png");
  }


 /* -- log frequency -- */

  runtime->log_freq = atoi(ParamFileGet("log", 1));
  runtime->log_freq = MAX(runtime->log_freq, 1);
  
 /* -- set default for remaining output type -- */

  while (ipos < MAX_OUTPUT_TYPES){
    output = runtime->output + ipos;
    output->type   = -1;
    output->cgs    =  0;
    output->dt     = -1.0;
    output->dn     = -1;
    output->dclock = -1.0;
    ipos++;
  }

 /* -- analysis -- */

  if (ParamExist ("analysis")){
    runtime->anl_dt = atof(ParamFileGet("analysis", 1));
    runtime->anl_dn = atoi(ParamFileGet("analysis", 2));
  }else{
    runtime->anl_dt = -1.0;   /* -- defaults -- */
    runtime->anl_dn = -1;
  }
#endif

#ifdef CHOMBO

/* ------------------------------------------------------------
   [Chombo HDF5 output] section
   ------------------------------------------------------------ */

/* ---- set output directory ---- */

  sprintf (runtime->output_dir, "%s","./");  /* default value is current directory */
  if (ParamExist("Output_dir")){
    str = ParamFileGet("Output_dir",1);
    sprintf (runtime->output_dir, "%s",str);
  }

/* -- check if we have write access and if the directory exists -- */

  sprintf (str_var,"%s/tmp09123.txt",runtime->output_dir); /* -- test file -- */
  fp = fopen(str_var,"w");  /* -- open test file for writing -- */
  if (fp == NULL){
    printf ("! Setup: cannot access directory '%s'.\n", runtime->output_dir);
    printf ("!        Please check that the directory exists\n");
    printf ("!        and you have write permission.\n");
    QUIT_PLUTO(1);
  }else{
    fclose(fp);
    remove(str_var);  /* -- remove file -- */
  }
#endif

/* ------------------------------------------------------------
   [Parameters] Section 
   ------------------------------------------------------------ */

  fp = fopen(ini_file,"r");
  
/* -- find position at "[Parameters" -- */

  for (ipos = 0; ipos <= nlines; ipos++){ 
    fgets(str_var, 512, fp);
    
    if (strlen(str_var) > 0) {
      str = strtok (str_var,"]");
      if (strcmp(str,"[Parameters") == 0) break;
    }
  }

  fgets(str_var, 512, fp); 
  
  for (ip = 0; ip < USER_DEF_PARAMETERS; ip++){
    fscanf (fp,"%s \n", str_var);
    dbl_var = atof(ParamFileGet(str_var,1));
    runtime->aux[ip] = dbl_var;
    fgets(str_var, sizeof(str_var), fp); /* use fgets to advance to next line */
  }
  fclose(fp);

  return(0);
}
#undef COMPARE
#undef NOPT
#undef NLEN

/* ********************************************************************* */
void GetOutputFrequency(Output *output, const char *output_format)
/*!
 *  Set the intervals between output files. 
 *  This can be done in three different ways:
 *
 *  - dt: time interval in code units
 *  - dn: step interval
 *  - dclock: actual clock time (in hours)
 * 
 * However, dn and dclock are mutually exclusive.
 *
 *********************************************************************** */
{
  char *str;
  int len, nhrs, nmin, nsec;

/* -- time interval in code units (dt) -- */

  output->dt = atof(ParamFileGet(output_format, 1));

/* -- check the 2nd field and decide to set "dn" or "dclock" -- */

  str = ParamFileGet(output_format,2);
  len = strlen(str);
  if (str[len-1] == 'h'){
    output->dclock = atof(str);    /* clock interval in hours */
    nhrs = (int)output->dclock;    /* integer part */
    nmin = (int)((output->dclock - nhrs)*100.0); /* remainder in minutes */
    if (nmin >= 60){
      printf ("! OutputFrequency: number of minutes exceeds 60 in %s output\n",
               output_format);
      QUIT_PLUTO(1);
    }
    output->dclock = nhrs*3600.0 + nmin*60;  /* convert to seconds */
    output->dn     = -1;
  }else if (str[len-1] == 'm'){
    output->dclock = atof(str);      /* clock interval in minutes */
    nmin = (int)output->dclock;      /* integer part */
    nsec = (int)((output->dclock - nmin)*100.0); /* remainder in seconds */
    if (nsec >= 60){
      printf ("! OutputFrequency: number of seconds exceeds 60 in %s output\n",
               output_format);
      QUIT_PLUTO(1);
    }
    output->dclock = nmin*60.0 + nsec;
    output->dn     = -1;
  }else if (str[len-1] == 's'){
    output->dclock = atof(str);           /* clock interval in seconds */
    output->dn     = -1;
  }else{
    output->dclock = -1.0;    
    output->dn     = atoi(ParamFileGet(output_format, 2));
  }

}

/* ********************************************************************* */
static Runtime q;
void RuntimeSet(Runtime *runtime)
/*!
 *  Store a static copy of the runtime structure for later access.
 *********************************************************************** */
{
  q = *runtime;
}
Runtime *RuntimeGet(void)
{
  return &q;
}
