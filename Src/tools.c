/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Collection of general-purpose functions.

  This file contains some miscellaneous tools for

  - Debugging / printing / error functions such as Trace(), Show(), 
    Where(), PlutoError
  - function for swapping/detecting endianity
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   July 15, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

   
/* ********************************************************************* */
int CheckNaN (double **u, int is, int ie, int id)
/*!
 * Check whether the array \c u contains Not-a-Number
 *  (NaN). QUIT if true.
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
void GetNeighbourRanks (Grid *grid, int **nranks)
/*! 
 *  Find the ranks of neighbour processors direction by direction.
 *  Store them in the array \c nrank[dir][s] where \c dir is the
 *  direction and <tt> s = 0,1 </tt> stands for the left (0) or right
 *  (1) processor.
 *  If a processor touches the physical boundary that is not periodic,
 *  the corresponding neighbour rank will be set to -1.
 *
 * \param [in]   grid    pointer to an array of grid structures.
 * \param [out]  nranks  an array of integers containing the ranks of
 *                       the neighbouring processors
 *********************************************************************** */
{
#ifdef PARALLEL
  int dir, coords[3], rnk;
  int lbound, rbound;
  MPI_Comm cartcomm;
  

  AL_Get_cart_comm(SZ, &cartcomm);
    
/* ------------------------------------------------------------
    Neighbour exists when there's no physical boundary, or
    when PERIODIC or SHEARING boundary conditions are imposed
    along that direction.
   ------------------------------------------------------------ */

  for (dir = 0; dir < 3;dir++){

  /* -- Obtain local processor coordinates -- */

    D_EXPAND(coords[IDIR] = grid->rank_coord[IDIR];  ,
             coords[JDIR] = grid->rank_coord[JDIR];  ,
             coords[KDIR] = grid->rank_coord[KDIR];)

    nranks[dir][0] = nranks[dir][1] = -1;

    lbound = grid->lbound[dir];
    if (lbound == 0 || lbound == PERIODIC || lbound == SHEARING){
      coords[dir] = grid->rank_coord[dir] - 1;
      MPI_Cart_rank(cartcomm, coords, &rnk);
      nranks[dir][0] = rnk;
    }

    rbound = grid->rbound[dir];

    if (rbound == 0 || rbound == PERIODIC || rbound == SHEARING){
      coords[dir] = grid->rank_coord[dir] + 1;
      MPI_Cart_rank(cartcomm, coords, &rnk);
      nranks[dir][1] = rnk;
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  return;
#endif  /* PARALLEL */
}

#if !(defined CHOMBO) && (defined GNUPLOT)
/* ********************************************************************* */
void GnuplotSetting(Grid *grid)
/*
 * Set-up a gnuplot script containing grid info, variable names, etc...
 *********************************************************************** */
{
  int  nv, count;
  char *var_name[256];
  FILE *fp;

/* --------------------------------------------------------
   0. Allocate memory for an array of strings containing
      the names of the variables.
      Get the variable names being written to disk.
   -------------------------------------------------------- */

  for (nv = 0; nv < 255; nv++) var_name[nv] = ARRAY_1D(32,char);
  count = GetOutputVarNames(FLT_OUTPUT, var_name);

/* --------------------------------------------------------
   1. Open file and write relevant info.
   -------------------------------------------------------- */

  fp = fopen("pluto.gp","w");

  fprintf (fp, "# ***********************************************************\n");
  fprintf (fp, "# Initialize grid quantities for reading binary data with \n");
  fprintf (fp, "# Gnuplot.\n");
  fprintf (fp, "# ***********************************************************\n");

  fprintf (fp, "nx = %d\n", grid->np_int_glob[IDIR]);
  fprintf (fp, "ny = %d\n", grid->np_int_glob[JDIR]);
  fprintf (fp, "nz = %d\n", grid->np_int_glob[KDIR]);
  
  fprintf (fp, "xbeg = %f; xend = %f\n", grid->xbeg_glob[IDIR],
                                         grid->xend_glob[IDIR]);

  fprintf (fp, "ybeg = %f; yend = %f\n", grid->xbeg_glob[JDIR],
                                         grid->xend_glob[JDIR]);
    
  fprintf (fp, "\n# -- Set variables indices --\n\n");
  for (nv = 0; nv < count; nv++) fprintf (fp, "%s = %d\n",var_name[nv], nv);
  
  fclose(fp);
}
#endif

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
void MakeState (Sweep *sweep)
/*!
 * Allocate memory areas for arrays inside the sweep
 * structure.
 *********************************************************************** */
{
  State *stateC = &(sweep->stateC);
  State *stateL = &(sweep->stateL);
  State *stateR = &(sweep->stateR);

/* --------------------------------------------------------
   0. Allocate memory for sweep structure members
   -------------------------------------------------------- */

  sweep->vn      = ARRAY_2D(NMAX_POINT, NVAR, double);
  sweep->flux    = ARRAY_2D(NMAX_POINT, NVAR, double);
  sweep->src     = ARRAY_2D(NMAX_POINT, NVAR, double);

  sweep->tc_flux = ARRAY_2D(NMAX_POINT, NVAR, double);    

  sweep->rhs     = ARRAY_2D(NMAX_POINT, NVAR, double);
  sweep->press   = ARRAY_1D(NMAX_POINT, double);
  sweep->bn      = ARRAY_1D(NMAX_POINT, double);
  sweep->SL      = ARRAY_1D(NMAX_POINT, double);
  sweep->SR      = ARRAY_1D(NMAX_POINT, double);

/* -- eigenvectors -- */

  sweep->lmax    = ARRAY_1D(NVAR, double);
  sweep->flag    = ARRAY_1D(NMAX_POINT, unsigned char);

/* --------------------------------------------------------
   1. Allocate memory for state structure members.
      C/L stand, respectively, for cell center and
      left interfaces (i+1/2, L).
   -------------------------------------------------------- */
  
  StateStructAllocate (stateC);
  StateStructAllocate (stateL);

/* --------------------------------------------------------
   2. Allocate memory for the right state structure.
      Note that we add an offset +1 in order to access
      stateR->v[i-1].
   -------------------------------------------------------- */
  
  stateR->v      = ARRAY_2D(NMAX_POINT, NVAR, double)+1;
  stateR->u      = ARRAY_2D(NMAX_POINT, NVAR, double)+1;
  stateR->flux   = ARRAY_2D(NMAX_POINT, NVAR, double)+1;
  stateR->lambda = ARRAY_2D(NMAX_POINT, NVAR, double)+1;
  stateR->prs    = ARRAY_1D(NMAX_POINT, double)+1;
  stateR->a2     = ARRAY_1D(NMAX_POINT, double)+1;
  stateR->cw     = ARRAY_1D(NMAX_POINT, double)+1;
  stateR->h      = ARRAY_1D(NMAX_POINT, double)+1;
  stateR->Lp     = ARRAY_3D(NMAX_POINT, NFLX, NFLX, double)+1;
  stateR->Rp     = ARRAY_3D(NMAX_POINT , NFLX, NFLX, double)+1;
  stateR->J      = stateL->J;
  stateR->cCR    = stateL->cCR;
  stateR->Fcr    = ARRAY_2D(NMAX_POINT, 4, double) + 1; 
  stateR->fluxCR = ARRAY_2D(NMAX_POINT, NVAR, double)+1;
  stateR->Bbck   = stateL->Bbck;

//  StateStructAllocate (stateR);

/* -- Allocate shared memory areas -- */

//  stateL->Bbck = ARRAY_2D(NMAX_POINT, 3, double);
//  stateR->Bbck = stateL->Bbck;
//  stateC->Bbck = ARRAY_2D(NMAX_POINT, 3, double);
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
void StateStructAllocate (State *p)
/*!
 * Allocate memory areas for arrays inside the State structure.
 *********************************************************************** */
{
  p->v      = ARRAY_2D(NMAX_POINT, NVAR, double);
  p->u      = ARRAY_2D(NMAX_POINT, NVAR, double);
  p->flux   = ARRAY_2D(NMAX_POINT, NVAR, double);
  p->lambda = ARRAY_2D(NMAX_POINT, NVAR, double);
  p->prs    = ARRAY_1D(NMAX_POINT, double);
  p->a2     = ARRAY_1D(NMAX_POINT, double);
  p->cw     = ARRAY_1D(NMAX_POINT, double);
  p->h      = ARRAY_1D(NMAX_POINT, double);
  p->Lp     = ARRAY_3D(NMAX_POINT, NFLX, NFLX, double);
  p->Rp     = ARRAY_3D(NMAX_POINT, NFLX, NFLX, double);
  p->J      = ARRAY_2D(NMAX_POINT, 3, double);
  p->cCR    = ARRAY_2D(NMAX_POINT, 3, double);
  p->Fcr    = ARRAY_2D(NMAX_POINT, 4, double); 
  p->fluxCR = ARRAY_2D(NMAX_POINT, NVAR, double);
  p->Bbck   = ARRAY_2D(NMAX_POINT, 3, double);
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
    print ("%8.3e  ", a[ip][nv]);
  }
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
    print ("\n");
  }
  print ("----------------------------------------------------------------\n");
}

/* ********************************************************************* */
void ShowVector (double *v, int n)
/*! 
 * Print the first n components of the vector v[]  
 *
 *********************************************************************** */
{
  int k;

  for (k = 0; k < n; k++)  print ("%12.6e  ", v[k]);
  print ("\n");
}

/* ********************************************************************* */
int StringArrayIndex (char *str_arr[], int size, char *str)
/*!
 *  Find the index of an array of strings whose value matches the
 *  input string \c *str.
 *
 * \param [in]  str_arr  the array of strings
 * \param [in]  size     the size of the array
 * \param [in]  str      the string to search for
 *********************************************************************** */
{
  int i;

  for (i = 0; i < size; i++){
    if (!strcmp(str_arr[i],str)) return i;
  }
  return -1;

}

/* ********************************************************************* */
void SymmetryCheck (Data_Arr V, int where, RBox *box)
/*!
 * Check if vectors preserve symmetry / antisymmetry properties
 * at a symmetry axis (r = 0 for cylindrical or theta = 0 for spherical)
 * Vectors can only be cell-centered (meaning that *all* components must
 * be cell-centered) or face-centered (meaning that *all* components must
 * be face-centered).
 * Staggered fields are not allowed.
 *********************************************************************** */
{
  int i,j,k;
  int ip, jp;
  double ***Vx1 = V[0];
  double ***Vx2 = V[1];
  double ***Vx3 = V[2];
  double scrh1, scrh2, scrh3;

  BOX_LOOP(box,k,j,i){

    #if GEOMETRY == CYLINDRICAL
    if (i >= IBEG) continue;

    if (where == X1FACE){   /* Vector is defined at radial i+1/2 interfaces */
      ip = 2*IBEG-i-2;
      jp = j;

    /* At IBEG-1/2 interface, check that there're no Inf or Nan values */
      if (i == IBEG-1){
        if (   fabs(Vx1[k][j][i]) > 1.e8 || Vx1[k][j][i] != Vx1[k][j][i]
            || fabs(Vx2[k][j][i]) > 1.e8 || Vx2[k][j][i] != Vx2[k][j][i]
            || fabs(Vx3[k][j][i]) > 1.e8 || Vx3[k][j][i] != Vx3[k][j][i]){
          print ("! SymmetryCheck(): invalid value\n");
          QUIT_PLUTO(1);
        }
        continue; /* Skip this point since i = ip would be the same */
      }
    }else {                /* Vector is defined at center or j+1/2 interface */
      ip = 2*IBEG-i-1;
      jp = j;
    }
    scrh1 = fabs(Vx1[k][j][i] + Vx1[k][jp][ip]); /* r-component:   antisymmetric */
    scrh2 = fabs(Vx2[k][j][i] - Vx2[k][jp][ip]); /* z-component:   symmetric */
    scrh3 = fabs(Vx3[k][j][i] + Vx3[k][jp][ip]); /* phi-component: antisymmetric */
    #elif GEOMETRY == SPHERICAL
    if (j >= JBEG) continue;

    if (where == X2FACE){  /* Vector is define at meridional j+1/2 interfaces */
      ip = i;
      jp = 2*JBEG-j-2;
    /* At JBEG-1/2 interface, check that there're no Inf or Nan values */
      if (j == JBEG-1){
        if (   fabs(Vx1[k][j][i]) > 1.e8 || Vx1[k][j][i] != Vx1[k][j][i]
            || fabs(Vx2[k][j][i]) > 1.e8 || Vx2[k][j][i] != Vx2[k][j][i]
            || fabs(Vx3[k][j][i]) > 1.e8 || Vx3[k][j][i] != Vx3[k][j][i]){
          print ("! SymmetryCheck(): invalid value\n");
          print ("! V = (%12.6e, %12.6e, %12.6e\n",
                    Vx1[k][j][i],Vx2[k][j][i],Vx3[k][j][i]);
          QUIT_PLUTO(1);
        }
        continue; /* Skip this point since j = jp would be the same */
      }
    }else {             /* Vector is define at center or radial interfaces */     
      ip = i;
      jp = 2*JBEG-j-1;
    }
    scrh1 = fabs(Vx1[k][j][i] - Vx1[k][jp][ip]); /* r-component:   symmetric */
    scrh2 = fabs(Vx2[k][j][i] + Vx2[k][jp][ip]); /* th-component:  antisymmetric */
    scrh3 = fabs(Vx3[k][j][i] + Vx3[k][jp][ip]); /* phi-component: antisymmetric */
    #endif


    if (scrh1 > 1.e-14 || scrh2 > 1.e-14 || scrh3 > 1.e-14){
      #if GEOMETRY == CYLINDRICAL
      if (where == X1FACE){
        print ("! SymmetryCheck(): Vector not symmetric at i+1/2,j = %d+1/2, %d\n",i,j);
      }else if (where == X2FACE){
        print ("! SymmetryCheck(): Vector not symmetric at i,j+1/2 = %d, %d+1/2\n",i,j);
      }else {
        print ("! SymmetryCheck(): Vector not symmetric at i,j = %d, %d\n",i,j);
      }
      #elif GEOMETRY == SPHERICAL
      if (where == X1FACE){
        print ("! SymmetryCheck(): Vector not symmetric at i+1/2,j = %d+1/2, %d\n",i,j);
      }else if (where == X2FACE){
        print ("! SymmetryCheck(): Vector not symmetric at i,j+1/2 = %d, %d+1/2\n",i,j);
      }else{
        print ("! SymmetryCheck(): Vector not symmetric at i,j = %d, %d\n",i,j);
      }
      #endif
      print ("! i, j, k = %d, %d, %d --> V = (%12.6e, %12.6e, %12.6e) \n",
                   i,j,k,Vx1[k][j][i],Vx2[k][j][i], Vx3[k][j][i]);
      print ("! ip,jp,k = %d, %d, %d --> V = (%12.6e, %12.6e, %12.6e)\n",
                 ip,jp,k, Vx1[k][jp][ip],Vx2[k][jp][ip], Vx3[k][jp][ip]);

      QUIT_PLUTO(1);
    } 
  }
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

  print ("Trace ------> %f ,  %d\n", xx, ++ik);
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
  static Grid *grid_copy;

/* --------------------------------------------------
    Keep a local copy of grid for subsequent calls
   -------------------------------------------------- */
 
  if (grid != NULL){
    grid_copy = grid;
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
    x1 = grid_copy->x[IDIR][ii];  ,
    x2 = grid_copy->x[JDIR][jj];  ,
    x3 = grid_copy->x[KDIR][kk];
  )

  D_SELECT(
    print ("  [i = %d], [x1 = %f]", ii, x1);  ,

    print ("  [i,j = %d, %d], [x1,x2 =  %f, %f]", ii, jj, x1, x2);  ,

    print ("  [i,j,k = %d, %d, %d], [x1,x2,x3 = %f, %f, %f]", ii, jj, kk,
            x1, x2, x3);
  )

  #ifdef CHOMBO
  print (", Level = %d\n", grid_copy->level);
  return;
  #endif
  print ("\n");
}

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

