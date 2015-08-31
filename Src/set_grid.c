/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Set up the global and local grids in the three coordinate 
         directions.

  Collects functions for allocating memory and defining grid 
  coordinates.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 24, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static void InitializeGrid (Runtime *, Grid *);
static void MakeGrid       (int, Runtime *, double *, double *, double *);
static void stretch_fun (double, double *, double *, double *);

/* ********************************************************************* */
void SetGrid (Runtime *rtime, Grid *GXYZ)
/*! 
 *
 * \param [in]  rtime   pointer to a Runtime  structure
 * \param [out] GXYZ    pointer to array of Grid structures
 * 
 *********************************************************************** */
{
  int  i, idim;
  int  iL, iR, ngh;
  char fname[512];
  double *dx, *xlft, *xrgt;
  Grid *G;
  FILE *fg;
  double xpatch_lft, xpatch_rgt;

  InitializeGrid(rtime, GXYZ);

  i    = MAX(rtime->npoint[0], MAX(rtime->npoint[1], rtime->npoint[2]));
  dx   = ARRAY_1D(i + 2, double);
  xlft = ARRAY_1D(i + 2, double);
  xrgt = ARRAY_1D(i + 2, double);

  for (idim = 0; idim < 3; idim++) {
     
    G   = GXYZ + idim;
    ngh = G->nghost;

    iL = ngh - 1;
    MakeGrid (idim, rtime, xlft, xrgt, dx);
/*MakeGrid (idim, rtime, G->xl_glob + iL, G->xr_glob + iL, G->dx_glob + iL); */

  /* ---- Assign values to grid structure members ----  */

    for (i = 1; i <= rtime->npoint[idim]; i++) {
      G->dx_glob[i + ngh - 1] = dx[i];
      G->xl_glob[i + ngh - 1] = xlft[i];
      G->xr_glob[i + ngh - 1] = xrgt[i];
    }

    iL = ngh;
    iR = rtime->npoint[idim] + ngh - 1;
    if (idim < DIMENSIONS){
      G->xr_glob[iL - 1] = G->xl_glob[iL];
      G->xl_glob[iR + 1] = G->xr_glob[iR];
    }
/*  ----  fill boundary values by copying adjacent cells  ---- */

    for (i = 0; i < ngh; i++) {
        
     /*  ---- left boundary  ----  */        
     
      G->dx_glob[i] = G->dx_glob[iL];
      G->xl_glob[i] = G->xl_glob[iL] - (ngh - i)*G->dx_glob[iL];
      G->xr_glob[i] = G->xl_glob[i] + G->dx_glob[iL];
      
     /*  ---- right boundary  ----  */   
          
      G->dx_glob[iR + i + 1] = G->dx_glob[iR];
      G->xl_glob[iR + i + 1] = G->xl_glob[iR] + (i + 1.0)*G->dx_glob[iR];
      G->xr_glob[iR + i + 1] = G->xl_glob[iR] + (i + 2.0)*G->dx_glob[iR];
    }

/*  ----  define geometrical cell center  ----  */

    for (i = 0; i <= iR + ngh; i++) {
      G->x_glob[i] = 0.5*(G->xl_glob[i] + G->xr_glob[i]);
    }
    
/*  ---- define leftmost and rightmost domain extrema  ---- */
    
    G->xi = G->xl_glob[iL];
    G->xf = G->xr_glob[iR];

    g_domBeg[idim] = G->xl_glob[iL];
    g_domEnd[idim] = G->xr_glob[iR];
  }

/*  ----  free memory  ----  */

  FreeArray1D(dx);
  FreeArray1D(xlft);
  FreeArray1D(xrgt);

/* ---------------------------------------------------------------------
                       Write grid file 
   --------------------------------------------------------------------- */

  sprintf (fname,"%s/grid.out", rtime->output_dir);
#ifdef PLUTO3_Grid
  fg = fopen(fname,"w");
  for (idim = 0; idim < 3; idim++) {
    G   = GXYZ + idim;
    ngh = G->nghost;
    iL  = G->gbeg;
    iR  = G->gend;
    
    fprintf (fg, "%d \n", iR - iL + 1);
    for (i = iL; i <= iR; i++) {
      fprintf (fg, " %d   %12.6e    %12.6e  %12.6e  %12.6e\n", 
        i-ngh+1, G->xl_glob[i], G->x_glob[i], G->xr_glob[i], G->dx_glob[i]);
    }
  }
  fclose(fg);
#else
  if (prank == 0){
    time_t time_now;
    time(&time_now);

    fg = fopen(fname,"w");
    fprintf (fg, "# ******************************************************\n");
  #ifdef CHOMBO  
    fprintf (fg, "# PLUTO-Chombo %s (Base) Grid File\n",PLUTO_VERSION);
  #else
    fprintf (fg, "# PLUTO %s Grid File\n",PLUTO_VERSION);
  #endif

    fprintf (fg, "# Generated on  %s",asctime(localtime(&time_now)));

/*    fprintf (fg, "# %s\n", getenv("PWD")); */

    fprintf (fg, "# \n");
    fprintf (fg,"# DIMENSIONS: %d\n",DIMENSIONS);
    #if GEOMETRY == CARTESIAN 
     fprintf (fg, "# GEOMETRY:   CARTESIAN\n");
    #elif GEOMETRY == CYLINDRICAL
     fprintf (fg, "# GEOMETRY:   CYLINDRICAL\n");
    #elif GEOMETRY == POLAR
     fprintf (fg, "# GEOMETRY:   POLAR\n");
    #elif GEOMETRY == SPHERICAL
     fprintf (fg, "# GEOMETRY:   SPHERICAL\n");
    #endif
    for (idim = 0; idim < DIMENSIONS; idim++){
     fprintf (fg, "# X%d: [% f, % f], %d point(s), %d ghosts\n", idim+1,
              g_domBeg[idim], g_domEnd[idim], 
              GXYZ[idim].np_int_glob, GXYZ[idim].nghost);
    }
    fprintf (fg, "# ******************************************************\n");
    for (idim = 0; idim < 3; idim++) {
      G   = GXYZ + idim;
      ngh = G->nghost;
      iL  = G->gbeg;
      iR  = G->gend;
    
      fprintf (fg, "%d \n", iR - iL + 1);
      for (i = iL; i <= iR; i++) {
       fprintf (fg, " %d   %18.12e    %18.12e\n", 
          i-ngh+1, G->xl_glob[i], G->xr_glob[i]);
      }
    }
    fclose(fg);
  }
#endif

/*  ----  define geometry factors, vol, area, etc...  ----  */

  MakeGeometry(GXYZ);

/* ----------------------------------------
         print domain specifications
   ---------------------------------------- */

  for (idim = 0; idim < DIMENSIONS; idim++){
   print1 ("  X%d: [% f, % f], %d point(s), %d ghosts\n", idim+1,
            g_domBeg[idim], g_domEnd[idim], 
            GXYZ[idim].np_int_glob, GXYZ[idim].nghost);
  }
}

/* ********************************************************************* */
void FreeGrid (Grid *grid)
/*!
 * Free array memory allocated previously.
 *
 *********************************************************************** */
{
  int dir;
  
  for (dir = 0; dir < 3; dir++){
    FreeArray1D(grid[dir].x);
    FreeArray1D(grid[dir].xl);
    FreeArray1D(grid[dir].xr);
    FreeArray1D(grid[dir].dx);
    FreeArray1D(grid[dir].xgc);
    FreeArray1D(grid[dir].dV);
    FreeArray1D(grid[dir].A-1);
    FreeArray1D(grid[dir].r_1);
    FreeArray1D(grid[dir].ct);
    FreeArray1D(grid[dir].inv_dx);
    FreeArray1D(grid[dir].inv_dxi);
  }
}

/* ********************************************************************* */
void InitializeGrid (Runtime *rtime, Grid *GXYZ)
/*!
 * Allocate memory for grid arrays, set shortcut pointers 
 * for local grids.
 * 
 *
 *********************************************************************** */
{
  int    idim, ngh;
  int    np_int, np_tot, np_int_glob, np_tot_glob;
  struct GRID *G;

  for (idim = 0; idim < 3; idim++) {

    G = GXYZ + idim;
    if (GEOMETRY == CARTESIAN){
      G->uniform = rtime->grid_is_uniform[idim];
    }else{ 
      G->uniform = 0;
    }
   
/* ----------------------------------------------------------------
    Dimensions that are not used will have 1 grid point, no bound
   ---------------------------------------------------------------- */

    if (idim >= DIMENSIONS) {
      G->np_int = G->np_tot = G->np_int_glob = G->np_tot_glob = 1;
      ngh = G->nghost = G->beg = G->end = G->gbeg = G->gend = 
            G->lbeg   = G->lend = 0;
      G->nproc = 1;
    } else {
      ngh = G->nghost;
    }

    np_tot_glob = G->np_tot_glob;
    np_int_glob = G->np_int_glob;
    np_tot = G->np_tot;
    np_int = G->np_int;

/*  -----------------------------------------------------------
                 Memory allocation. 
    ----------------------------------------------------------- */
     
    G->x_glob  = ARRAY_1D(np_tot_glob, double);
    G->xr_glob = ARRAY_1D(np_tot_glob, double);
    G->xl_glob = ARRAY_1D(np_tot_glob, double);
    G->dx_glob = ARRAY_1D(np_tot_glob, double);

/*  ----  define shortcuts for local grids  ----  */

    G->x  = G->x_glob  + G->beg - ngh;
    G->xr = G->xr_glob + G->beg - ngh;
    G->xl = G->xl_glob + G->beg - ngh;
    G->dx = G->dx_glob + G->beg - ngh;
  }
}

/* ********************************************************************* */
void MakeGrid (int idim, Runtime *rtime, double *xlft, double *xrgt, double *dx)
/*! 
 *
 *  Build grid nodes as defined by pluto.ini.
 *
 *  Options are:
 *  
 *  'u' = uniform grid, simply defined as
 *
 *    dx = (xR - xL)/npoint, xleft(i) = xl + i*dx, xright(i) = xleft(i) + dx
 *  
 *  's' = stretched grid; solve
 *
 *          dx*( 1 + r + r^2 + r^3 + ... r^(N-1)) = xR - xL
 *
 *        in the stretching ratio r, provided dx, N, xR and xL are known.
 *        dx is taken from the closest uniform grid.
 *
 *  'l+' = logarithmic grid, mesh size increases with x; it is defined as
 *
 *
 *               x + |xL| - xL 
 *      y = Log[ ------------- ] , with uniform spacing  y(i+1/2) - y(i-1/2) = dy
 *                  |xL|
 *  
 *      dy    = (yR - yL)/N  and dx(i) becomes
 *      dx(i) = (x(i-) + fabs(xL) - xL)*(10^dy - 1.0);
 *
 *     NOTE: xR must be positive and xL can take any value different from 0
 *
 *  'l-' = logarithmic grid, mesh size decreases with x; it is defined as
 *
 *
 *               xR + |xL| - x 
 *      y = Log[ -------------- ] , with uniform spacing  y(i+1/2) - y(i-1/2) = dy
 *                    |xL|
 *
 *      dy    = -(yR - yL)/N   
 *      dx(i) = (x(i-) - fabs(xL) - xR)*(10^dy - 1.0);
 *     
 *********************************************************************** */
#define MAX_ITER   50
{
  int    i, iseg, n, nstart, nseg, npoint;
  int    log_inc, log_dec, next_seg_is;
  int    done_with_segment[16], all_segments_done;
  int    iL, iR;
  int    i_patch_lft[16], i_patch_rgt[16];
  double x_patch_lft[16], x_patch_rgt[16];
  double xR, xL;
  double par[4];
  double dalpha, alpha, f, df;
  double dy;

  nseg = rtime->npatch[idim];

/* ---------------------------------------------------
       for each patch, find the leftmost and 
       rightmost indexes ilft and irgt          
   --------------------------------------------------- */

  i_patch_lft[1] = 1;
  i_patch_rgt[1] = rtime->patch_npoint[idim][1];
  x_patch_lft[1] = rtime->patch_left_node[idim][1];
  x_patch_rgt[1] = rtime->patch_left_node[idim][2];
 
  for (iseg = 2; iseg <= nseg; iseg++) {
    i_patch_lft[iseg] = i_patch_rgt[iseg - 1] + 1;
    i_patch_rgt[iseg] = i_patch_lft[iseg] + rtime->patch_npoint[idim][iseg] - 1;
    x_patch_lft[iseg] = rtime->patch_left_node[idim][iseg];
    x_patch_rgt[iseg] = rtime->patch_left_node[idim][iseg + 1];
  }

  done_with_segment[0] = 0;
  done_with_segment[nseg + 1] = 0;

  for (iseg = 1; iseg <= nseg; iseg++) {
   
    done_with_segment[iseg] = 0;
    xL     = x_patch_lft[iseg];
    xR     = x_patch_rgt[iseg];
    iL     = i_patch_lft[iseg];
    iR     = i_patch_rgt[iseg];

    npoint = rtime->patch_npoint[idim][iseg];

/*  ----  first process only uniform or logarithmic grids  ----  */

    if (rtime->patch_type[idim][iseg] == UNIFORM_GRID) {

      for (i = iL; i <= iR; i++) {
        dx[i]   = (xR - xL)/(double)(npoint);
        xlft[i] = xL + (double)(i - iL)*dx[i];
        xrgt[i] = xlft[i] + dx[i];
      }
      done_with_segment[iseg] = 1;

    } else if ( rtime->patch_type[idim][iseg] == LOGARITHMIC_INC_GRID ||
                rtime->patch_type[idim][iseg] == LOGARITHMIC_DEC_GRID) {
 
      log_inc = rtime->patch_type[idim][iseg] == LOGARITHMIC_INC_GRID;
      log_dec = rtime->patch_type[idim][iseg] == LOGARITHMIC_DEC_GRID;

      dy  = log10( (xR + fabs(xL) - xL)/fabs(xL));
      dy /= (double) (npoint);

      if (log_dec) dy *= -1.0;

      xlft[iL] = xL;

      for (i = iL; i <= iR; i++) {
        if (log_inc) {
          dx[i] = (xlft[i] + fabs(xL) - xL)*(pow(10.0,dy) - 1.0);
        }else{
          dx[i] = (xlft[i] - fabs(xL) - xR)*(pow(10.0,dy) - 1.0);
        }
        xrgt[i]     = xlft[i] + dx[i];
        xlft[i + 1] = xrgt[i];
      }        
      done_with_segment[iseg] = 1;
   
    } else {
      continue;
    }
  }

  all_segments_done = 1;
  for (iseg = 1; iseg <= nseg; iseg++) { 
    all_segments_done = all_segments_done && done_with_segment[iseg];
  }

/* --------------------------------------------------------------
                  now do stretched grids ...      
   -------------------------------------------------------------- */
        
  while (!all_segments_done) {  /* loop until all segments are processed... */

/* ---------------------------------------------------
    scan the entire grid and process only those 
    segments close to a segment that is already done
    (i.e. done_with_segment[iseg] = 1)
   --------------------------------------------------- */

    for (iseg = 1; iseg <= nseg; iseg++) {

      xL     = x_patch_lft[iseg];
      xR     = x_patch_rgt[iseg];
      iL     = i_patch_lft[iseg];
      iR     = i_patch_rgt[iseg];
      npoint = rtime->patch_npoint[idim][iseg];

/*  -----------------------------------------------------
     now find whether the segment to right (iseg+1) or 
     to the left (iseg-1) is already done;     
    -----------------------------------------------------  */

      next_seg_is = 0;
      if (done_with_segment[iseg + 1]) next_seg_is = iseg + 1;
      if (done_with_segment[iseg - 1]) next_seg_is = iseg - 1;

/*  -----------------------------------------------------
     if none of them is done, skip this segment and 
      move forward,  otherwise process segment iseg,
      otherwise process the grid    
    -----------------------------------------------------  */

      if (rtime->patch_type[idim][iseg] == STRETCHED_GRID && 
          next_seg_is != 0 && !done_with_segment[iseg]) {

/*  ----------------------------------------------------------
     nstart is:
     
      *  the rightmost point if the next grid is done,
         i.e., next_seg_is > iseg;
      *  the leftmost  point if the previous grid is done,
         i.e., next_seg_is < iseg;  
    ----------------------------------------------------------  */

        nstart = (next_seg_is > iseg ? iR:iL);
        
        alpha  = 1.01;           /* provide a first guess */

/*  --------------------------------------------------------------------
     par[0] : contains the number of points  
     par[1] : Ratio L/dx : the first dx is the one that belongs to the 
              previous (iseg>1) or next (iseg = 1) uniform segment
    --------------------------------------------------------------------- */

        par[0] = (double)npoint;      /* Number of points */

        if (next_seg_is > iseg)
          par[1] = (xR - xL)/dx[nstart + 1];
        else
          par[1] = (xR - xL)/dx[nstart - 1];
  
/*  ----  Use NEWTON algorithm to find the stretch factor  ----  */

        for (n = 0; n <= MAX_ITER; n++) {
          if (n == MAX_ITER) {
            print1 ("Too many iterations during grid (%d) generation!\n",idim);
            QUIT_PLUTO(1);
          }

          stretch_fun (alpha, par, &f, &df);

          dalpha = f / df;
          alpha -= dalpha;

          if (fabs (dalpha) < 1.e-14*alpha)
            break;
        }

        if (alpha > 1.2) {
          print1 (" ! WARNING:  alpha=%12.6e > 1.05 , dimension %d\n", alpha, idim);
          print1 (" ! While stretching segment %d\n", iseg);
          QUIT_PLUTO(1);
        }

        print1 ("     - Stretched grid on dim %d (ratio: %f)\n",idim, alpha);

        if (next_seg_is > iseg) {
          xrgt[nstart] = xR;
          for (i = nstart; i >= iL; i--) {
            dx[i]       = pow (alpha, nstart - i + 1)*dx[nstart + 1];
            xlft[i]     = xrgt[i] - dx[i];
            xrgt[i - 1] = xlft[i]; 
          }
        } else {
          xlft[nstart] = xL;
          for (i = nstart; i <= iR; i++) {
            dx[i]       = pow (alpha, i - nstart + 1)*dx[nstart - 1];
            xlft[i + 1] = xrgt[i] = xlft[i] + dx[i];
          }
        }

/*  ----  ok, this segment has been processed  ----  */

        done_with_segment[iseg] = 1;

/*  ----  check whether there are other segments to process  ----  */

        all_segments_done = 1;
        for (i = 1; i <= nseg; i++) 
          all_segments_done = all_segments_done && done_with_segment[i];

/*  ----  exit from  segment loop (iseg) and rescan the grid from the beginning  ----   */

        break;
      }
    }
  }
}
#undef MAX_ITER

/* ############################################################## */
void stretch_fun (double x, double *par, double *f, double *dfdx)
/* 
 #
 # S
 #
 ################################################################ */
{
  double xN, L_dx, scrh;

  xN   = par[0];
  L_dx = par[1];

  scrh  = x*(pow (x, xN) - 1.0)/(x - 1.0);
  *f    = log (scrh) - log (L_dx);
  *dfdx = 1.0/x + xN*pow (x, xN - 1.0)/(pow (x, xN) - 1.0) -
          1.0/(x - 1.0);
}
