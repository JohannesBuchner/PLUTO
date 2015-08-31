/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Collection of general-purpose functions.

  This file contains some miscellaneous tools for

  - Debugging / printing / error functions such as Trace(), Show(), 
    Where(), PlutoError
  - function for swapping/detecting endianity
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 1, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
int CheckNaN (double **u, int is, int ie, int id)
/*!
 * Cheeck whether the array \c u contains Not-a-Number
 *  (NaN)
 *
 *********************************************************************** */
{
  int ii, nv, i, j;

  for (ii = is; ii <= ie; ii++) {
  for (nv = 0; nv < NVAR; nv++) {
    if (u[ii][nv] != u[ii][nv]) {
      print (" > NaN found (%d), |", id);
      Show (u, ii);
      QUIT_PLUTO(1);
    }
  }}
  return (0);
}

/* ********************************************************************* */
int IsLittleEndian (void) 
/*!
 * Return 1 if the current architecture has little endian order
 *
 *********************************************************************** */
{
  int TestEndian = 1;
  return *(char*)&TestEndian;
}

/* ********************************************************************* */
void MakeState (State_1D *state)
/*!
 *
 * Allocate memory areas for arrays inside the state
 * structure.
 *
 *
 *********************************************************************** */
{
  state->v       = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->vp      = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->vm      = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->up      = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->um      = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->flux    = ARRAY_2D(NMAX_POINT, NVAR, double);

  state->src     = ARRAY_2D(NMAX_POINT, NVAR, double);

  state->visc_flux = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->visc_src  = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->tc_flux   = ARRAY_2D(NMAX_POINT, NVAR, double);    
  state->res_flux  = ARRAY_2D(NMAX_POINT, NVAR, double); 

  state->rhs     = ARRAY_2D(NMAX_POINT, NVAR, double);
  state->press   = ARRAY_1D(NMAX_POINT, double);
  state->bn      = ARRAY_1D(NMAX_POINT, double);
  state->SL      = ARRAY_1D(NMAX_POINT, double);
  state->SR      = ARRAY_1D(NMAX_POINT, double);

/* -- eigenvectors -- */

  state->Lp      = ARRAY_3D(NMAX_POINT, NFLX, NFLX, double);
  state->Rp      = ARRAY_3D(NMAX_POINT, NFLX, NFLX, double);
  state->lambda  = ARRAY_2D(NMAX_POINT, NFLX, double);
  state->lmax    = ARRAY_1D(NVAR, double);

  state->a2   = ARRAY_1D(NMAX_POINT, double);
  state->h    = ARRAY_1D(NMAX_POINT, double);

/*  state->dwlim   = ARRAY_2D(NMAX_POINT, NVAR, double);*/

  state->flag    = ARRAY_1D(NMAX_POINT, unsigned char);

/* --------------------------------------
     define shortcut pointers for
     left and right values with respect
     to the cell center
   -------------------------------------- */
   
  state->vL = state->vp;
  state->vR = state->vm + 1;

  state->uL = state->up;
  state->uR = state->um + 1;

  #if (TIME_STEPPING == HANCOCK) || (TIME_STEPPING == CHARACTERISTIC_TRACING)
   state->vh = ARRAY_2D(NMAX_POINT, NVAR, double);
  #else
   state->vh = state->v;
  #endif

}

/* ********************************************************************* */
void PlutoError (int condition, char *str)
/*!
 * If condition is true, issue an error and quit the code.
 *
 *********************************************************************** */
{
  char *str_err="! Error: ";

  if (condition) {
    print (str_err);
    print (str);
    print ("\n");
    QUIT_PLUTO(1);
  }
}

/* ********************************************************************* */
void Show (double **a, int ip)
/*! 
 * Print the component of the array \c a at grid index \c ip  
 *
 *********************************************************************** */
{
  int nv, ix, iy, iz;

  if (g_dir == IDIR) {
    print ("X-sweep");
    ix = ip;
    iy = g_j;
    iz = g_k;
  } else if (g_dir == JDIR) {
    print ("Y-sweep");
    ix = g_i;
    iy = ip;
    iz = g_k;
  } else if (g_dir == KDIR) {
    print ("Z-sweep");
    ix = g_i;
    iy = g_j;
    iz = ip;
  }

  D_SELECT( print (" (%d)> ", ix);     ,
            print (" (%d,%d)> ", ix, iy);  ,
            print (" (%d,%d,%d)> ", ix, iy, iz);  )


  for (nv = 0; nv < NVAR; nv++) {
    print ("%12.6e  ", a[ip][nv]);
  }
  print ("\n");
}

/* ********************************************************************* */
void ShowVector (double *v, int n)
/*! 
 * Print the component of the array \c a at grid index \c ip  
 *
 *********************************************************************** */
{
  int k;

  for (k = 0; k < n; k++)  print ("%12.6e  ", v[k]);
  print ("\n");
}

/* ********************************************************************* */
void ShowMatrix(double **A, int n, double eps)
/*!
 * Make a nice printing of a 2D matrix \c A[0..n-1][0..n-1]
 * Entries with values below eps will display "0.0"
 *
 *
 *********************************************************************** */
{
  int k1,k2;

  print ("----------------------------------------------------------------\n");
  for (k1 = 0; k1 < n; k1++){
    for (k2 = 0; k2 < n; k2++){
      print ("%12.3e   ", fabs(A[k1][k2]) > eps ? A[k1][k2]:0.0);
    }
    printf ("\n");
  }
  print ("----------------------------------------------------------------\n");
}

/* ********************************************************************* */
void SwapEndian (void *x, const int nbytes) 
/*!
 * Swap the byte order of x.
 * 
 * \param [in] x      pointer to the variable being swapped
 * \param [in] nbytes data type size
 * \return This function has no return value.
 *********************************************************************** */
{
  int k;
  static char Swapped[16];
  char *c;

  c = (char *) x;

  for (k = nbytes; k--; ) Swapped[k] = *(c + nbytes - 1 - k);
  for (k = nbytes; k--; ) c[k] = Swapped[k];
}

/* ********************************************************************* */
void Trace (double xx)
/*!
 * Print a number xx and the number of times it has been called.
 *
 *********************************************************************** */
{
  static int ik;

  printf ("Trace ------> %f ,  %d\n", xx, ++ik);
}

/* ********************************************************************* */
void Where (int i, Grid *grid)
/*!
 *  Print the location of a particular zone (i,j,k)
 *  in the computational domain.
 *  \note This function must be initialized before using it 
 *        to store grid information. This is done  by calling 
 *        Where(i, grid) the very first time.
 *        Subsequent calls can be then done by simply using 
 *        Where(i,NULL). 
 *
 *********************************************************************** */
{
  int    ii=0, jj=0, kk=0;
  double x1, x2, x3;
  static Grid *grid1, *grid2, *grid3;

/* --------------------------------------------------
    Keep a local copy of grid for subsequent calls
   -------------------------------------------------- */
 
  if (grid != NULL){
    grid1 = grid + IDIR;
    grid2 = grid + JDIR;
    grid3 = grid + KDIR;
    return;
  }

  #ifdef CH_SPACEDIM
   if (g_intStage < 0) return; /* HOT FIX used by CHOMBO
                             (g_intStage = -1) when writing HDF5 file */
  #endif

/* -- ok, proceed normally -- */
  
  if (g_dir == IDIR){
    D_EXPAND(ii = i;, jj = g_j;, kk = g_k;)
  }else if (g_dir == JDIR){
    D_EXPAND(ii = g_i;, jj = i;, kk = g_k;)
  }else if (g_dir == KDIR){
    D_EXPAND(ii = g_i;, jj = g_j;, kk = i;)
  }

  D_EXPAND(
    x1 = grid1->x[ii];  ,
    x2 = grid2->x[jj];  ,
    x3 = grid3->x[kk];
  )

  D_SELECT(
    print ("zone [x1(%d) = %f]",
            ii, grid1->x[ii]);  ,

    print ("zone [x1(%d) = %f, x2(%d) = %f]",
            ii, grid1->x[ii], 
            jj, grid2->x[jj]);  ,

    print ("zone [x1(%d) = %f, x2(%d) = %f, x3(%d) = %f]",
            ii, grid1->x[ii], 
            jj, grid2->x[jj],
            kk, grid3->x[kk]);
  )

  #ifdef CHOMBO
   print (", Level = %d\n", grid1->level);
   return;
  #endif
  #ifdef PARALLEL
   print (", proc %d\n", prank);
   return;
  #else
   print ("\n");
   return;
  #endif  
}

/* /////////////////////////////////////////////////////////////////////
    The next set of functions provides basic functionalities to
     
     - set the log file
     - formatted output to the log file through the print() and print1()
       functions
   ///////////////////////////////////////////////////////////////////// */

static char log_file_name[512];

/* ********************************************************************* */
int SetLogFile(char *output_dir, Cmd_Line *cmd)
/*!
 * Set the name of the log file and open in write or append mode
 * depending on whether restart is enabled or not.
 *
 * \param [in] output_dir  the name of the output directory
 * \param [in] cmd         pointer to cmd line option structure.
 *                           
 *********************************************************************** */
{
#if PRINT_TO_FILE == YES
  FILE *fl;

/* ------------------------------------------------
    All processors set log file name
   ------------------------------------------------ */

  sprintf (log_file_name, "%s/pluto.log",output_dir);    

/* ------------------------------------------------
    Proc. #0 opens log file for writing if 
    -restart or -h5restart have not been given.
    Otherwise, open in append mode.
   ------------------------------------------------ */
  
  if (prank == 0){
    if (cmd->restart == NO && cmd->h5restart == NO){
      fl = fopen(log_file_name,"w");
    }else{
      fl = fopen(log_file_name,"aw");
    } 

  /* -- check that we have a valid directory name -- */

    if (fl == NULL){
      printf ("! SetLogFile: pluto.log cannot be written.\n");
      QUIT_PLUTO(1);
    }
    fprintf(fl,"\n");
    fclose(fl);
  }
#endif
}

#ifndef CH_SPACEDIM
/* ********************************************************************* */
void print (const char *fmt, ...)
/*!
 * Define print function for the static grid version
 * of PLUTO. The Chombo version is defined in Chombo/amrPLUTO.cpp
 *
 *********************************************************************** */
{
  FILE *fl;

  va_list args;
  va_start(args, fmt);

#if PRINT_TO_FILE == YES
  fl = fopen(log_file_name,"a");
  vfprintf(fl, fmt, args);
  fclose (fl);
#else
  vprintf(fmt, args);
#endif

  va_end(args);
}
/* ********************************************************************* */
void print1 (const char *fmt, ...)
/*!
 *
 *   Define print1 function
 *
 *********************************************************************** */
{
  FILE *fl;

  va_list args;
  va_start(args, fmt);

  #if PRINT_TO_FILE == YES
   if (prank == 0){
     fl = fopen(log_file_name,"a");
     vfprintf(fl,fmt, args);
     fclose(fl);
   }
  #else
   if (prank == 0) vprintf(fmt, args);
  #endif
  va_end(args);
}
#endif

/* ********************************************************************* */
void WriteAsciiFile (char *fname, double *q, int nvar)
/*! 
 *  Write one row a multi-column ascii file 
 *
 *********************************************************************** */
{
  int n;
  static char old_file[64] = "\0";
  FILE *fp;
  
/* --------------------------------------------------------
    If the old file name matches the new one, then open 
    in append mode. Otherwise, open in write mode.
   -------------------------------------------------------- */
   
  if (strcmp (fname,old_file) == 0){  
    fp = fopen(fname,"a");
  }else{
    fp = fopen(fname,"w");
    sprintf (old_file,"%s",fname);
  }
  
  for (n = 0; n < nvar; n++) fprintf (fp,"%12.6e  ",q[n]);
  fprintf (fp,"\n");
  fclose(fp);
  
}
