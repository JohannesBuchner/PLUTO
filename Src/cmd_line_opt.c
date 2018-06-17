/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Parse command line options.

  Parse command line options at runtime and set values in
  the Cmd_Line structure.

  \authors A. Mignone (mignone@ph.unito.it)
  \date    March 16, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static void PrintUsage();

/* ********************************************************************* */
void ParseCmdLineArgs (int argc, char *argv[], char *ini_file, 
                       Cmd_Line *cmd)
/*!
 *
 * Parse command line options. Error messages will be output to
 * stdout.
 *
 * \param[in]  argc       argument count
 * \param[in]  argv       argument vector
 * \param[in]  ini_file   a pointer to a string containing the name of
 *                        the initialization file (default: "pluto.ini")
 * \param[out] cmd        the command-line structure.
 *
 *********************************************************************** */
{
  int i,j;

/* -----------------------------------------------
           Set default values first
   ----------------------------------------------- */

  cmd->restart   = NO;
  cmd->h5restart = NO;
  cmd->maxsteps  = -1;
  cmd->write     = YES;
  cmd->makegrid  = NO; 
  cmd->jet       = -1; /* -- means no direction -- */
  cmd->xres      = -1; /* -- means no grid resizing -- */

  cmd->nproc[IDIR] = -1; /* means autodecomp will be used */
  cmd->nproc[JDIR] = -1;
  cmd->nproc[KDIR] = -1;

  #ifdef PARALLEL
   cmd->parallel_dim[IDIR] = YES;  /* by default, we parallelize */
   cmd->parallel_dim[JDIR] = YES;  /* all directions             */
   cmd->parallel_dim[KDIR] = YES;
  #endif

/* --------------------------------------------------------------------
             Parse Command Line Options.

    Note: Since at this time the output directory is now known and
          "pluto.log" has not been opened yet, we use printf to issue
          error messages.
   -------------------------------------------------------------------- */

  for (i = 1; i < argc ; i++){

    if (!strcmp(argv[i],"-dec")) {

    /* -- start reading integers at i+1 -- */

      for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){

        if ((++i) >= argc){
          if (prank == 0){
            D_SELECT(printf ("! You must specify -dec n1\n");  ,
                     printf ("! You must specify -dec n1  n2\n");  ,
                     printf ("! You must specify -dec n1  n2  n3\n");)
          }
          QUIT_PLUTO(1);
        }
        cmd->nproc[g_dir] = atoi(argv[i]);
        if (cmd->nproc[g_dir] == 0){
            if (prank == 0) {
              printf ("! Incorrect number of processor for g_dir = %d \n", g_dir);
            }
            QUIT_PLUTO(0);
        }
      }

    }else if (!strcmp(argv[i],"-i")) {

      sprintf (ini_file,"%s",argv[++i]);
    
    } else if (!strcmp(argv[i],"-makegrid")) {

      cmd->makegrid = YES;

    }else if (!strcmp(argv[i],"-maxsteps")){

      if ((++i) >= argc){
        if (prank == 0) printf ("! You must specify -maxsteps nn\n");
        QUIT_PLUTO(1);
      }else{
        cmd->maxsteps = atoi(argv[i]);
        if (cmd->maxsteps < 0) {
          if (prank == 0)
            printf ("! You must specify -maxsteps nn, with nn >= 0 \n");
            QUIT_PLUTO(0);
        }
      }

    }else if (!strcmp(argv[i],"-no-write")) {

      cmd->write = NO;

    }else if (!strcmp(argv[i],"-no-x1par")) {

      cmd->parallel_dim[IDIR] = NO;

    }else if (!strcmp(argv[i],"-no-x2par")) {

      cmd->parallel_dim[JDIR] = NO;

    }else if (!strcmp(argv[i],"-no-x3par")) {

      cmd->parallel_dim[KDIR] = NO;

    } else if (!strcmp(argv[i],"-x1jet")) {

      cmd->jet = IDIR;

    } else if (!strcmp(argv[i],"-x2jet")) {

      cmd->jet = JDIR;

    } else if (!strcmp(argv[i],"-x3jet")) {

      cmd->jet = KDIR;

    }else if (!strcmp(argv[i],"-xres")){

      if ((++i) >= argc){
        if (prank == 0) printf ("! You must specify -xres nn\n");
        QUIT_PLUTO(1);
      }else{
        cmd->xres = atoi(argv[i]);
        if (cmd->xres <= 1) {
          if (prank == 0) printf ("! You must specify -xres nn, with nn > 1 \n");
          QUIT_PLUTO(0)
        }
      }

    }else  if (!strcmp(argv[i],"-restart") || !strcmp(argv[i],"-h5restart")) {

     /* --------------------------------------------- 
          default restart is last written file (-1)
        --------------------------------------------- */

      if (!strcmp(argv[i], "-restart")) cmd->restart   = YES;  /* can only take YES/NO values */
      else                              cmd->h5restart = YES;
      cmd->nrestart = -1;   /* the file number to restart from */

      if ((++i) < argc){
        char *endptr;
/*         cmd->restart = atoi(argv[i]); */
        cmd->nrestart = (int)strtol(argv[i], &endptr, 10);

      /* ----------------------------------------------
          if a non-numerical character is encountered, 
          cmd->nrestart should reset to -1
         ---------------------------------------------- */

        if (endptr == argv[i]){
          i--;
          cmd->nrestart = -1;
        }
      }

    }else if (!strcmp(argv[i],"--help")){

      PrintUsage();
      QUIT_PLUTO(1);

    }else{
      if (prank == 0) printf ("! Unknown option '%s'\n",argv[i]);
      QUIT_PLUTO(1);
    }
  }

/* -- disable domain decomposition in the 
      direction specified by -xnjet  -- */

  if      (cmd->jet == IDIR) cmd->parallel_dim[IDIR] = NO;
  else if (cmd->jet == JDIR) cmd->parallel_dim[JDIR] = NO;
  else if (cmd->jet == KDIR) cmd->parallel_dim[KDIR] = NO;

}
/* ******************************************************************* */
void PrintUsage()
/*
 *
 *
 ********************************************************************* */
{

  if (prank != 0) return;
  
  printf ("Usage: pluto [options]\n\n");
  printf ("           or \n\n");
  printf ("       mpirun -np NP ./pluto [options]\n\n");
  printf ("[options] are:\n\n");
  printf (" -dec n1 [n2] [n3]\n");  
  printf ("    Enable user-defined parallel decomposition mode. The integers\n");
  printf ("    n1, n2 and n3 specify the number of processors along the x1,\n");
  printf ("    x2, and x3 directions. There must be as many integers as the\n");
  printf ("    number of dimensions and their product must equal the total\n");
  printf ("    number of processors used by mpirun or an error will occurr.\n\n"); 

  printf (" -i <name>\n");
  printf ("    Use <name> as initialization file instead of pluto.ini.\n\n");

  printf (" --help\n");
  printf ("    Show this option summary.\n\n");

  printf (" -h5restart n\n");
  printf ("    Restart computations from the n-th output file in HDF5\n");
  printf ("    double precision format (.dbl.h5).\n\n");

  printf (" -makegrid\n");
  printf ("    Generate grid only, do not start computations.\n\n");

  printf (" -maxsteps n\n");
  printf ("    Stop computations after n steps.\n\n");

  printf (" -no-write\n");
  printf ("    Do not write data to disk.\n\n");
  
  printf (" -no-x1par, -no-x2par, -no-x3par\n");
  printf ("    Do not perform parallel domain decomposition along the x1, x2\n");
  printf ("    or x3 direction, respectively.\n\n");

  printf (" -restart n\n");
  printf ("    Restart computations from the n-th output file in double in\n");
  printf ("    precision format (.dbl).\n\n");

  printf (" -show-dec\n");
  printf ("    Show domain decomposition when running in parallel mode.\n\n");
  
  printf (" -x1jet, -x2jet, -x3jet\n");
  printf ("    Exclude from integration regions of zero pressure gradient\n");
  printf ("    that extends up to the end of the domain in x1, x2 or x3\n");
  printf ("    direction, respectively. This option is specifically\n");
  printf ("    designed for jets propagating along one of the coordinate\n");
  printf ("    axis. In parallel mode, parallel decomposition is not\n");
  printf ("    performed along the selected direction.\n\n");

  printf (" -xres n1\n");
  printf ("    Set the grid resolution in the x1 direction to n1 zones\n");
  printf ("    by overriding pluto.ini. Cell aspect ratio is preserved by\n");
  printf ("    modifying the grid resolution in the other coordinate\n"); 
  printf ("    directions accordingly.\n");

}
