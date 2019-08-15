#include"pluto.h"

/* ********************************************************************** */
void Flatten (const Sweep *sweep, int beg, int end, Grid *grid)
/*
 *
 * Flatten distribution using 
 *
 *  - multidimensional flattening (revert to minmod limiter)
 *  - oned flattening (decrease slopes as in the original 
 *    PPM algorithm)
 *
 *
 *
 ************************************************************************ */
#if SHOCK_FLATTENING == ONED

/* -----------------------------------------
    One dimensional shock flattening:
    Requires 7 point stencil:
  
         |---|---|---|---|---|---|---|
                       i
   f_t             o   o   o
           x   x   x   x   x   x   x

    It is kept mostly for backward 
    compatibility reasons.
   ----------------------------------------- */

  /* ---- PPM PARAMETERS ---- */

#if PHYSICS == HD

 #define   EPS2     0.33
 #define   OME1     0.75
 #define   OME2     10.0

#endif

#if PHYSICS == MHD

 #define   EPS2     0.33
 #define   OME1     0.75
 #define   OME2     10.0

#endif

#if PHYSICS == RHD || PHYSICS == RMHD

 #define   EPS2     1.0
 #define   OME1     0.52
 #define   OME2     10.0

#endif

{
  int    i, nv, sj;

  double scrh, dp, d2p, min_p, vf, fj;
  double **v, **vp, **vm;
  static double *f_t;
   
#if EOS == ISOTHERMAL 
  int PRS = RHO;
#endif
   
  if (f_t == NULL){
    f_t   = ARRAY_1D(NMAX_POINT, double);    
  }

  v  = sweep->stateC.v;
  vp = sweep->stateL.v;
  vm = sweep->stateR.v-1;

/* ---------------------------------------------------------
     the following constraints is necessary for
     a particular combinations of algorithms:
 
      - CTU + CT + MHD + userdef boundary 
      - shock flattening.

     This is necessary when user defined boundary conditions 
     are used in the fully corner coupled unsplit schemes 
     (HANCOCK and CHARACTERISTIC_TRACING), since the 
     grid is expanded by one point in the userdef boundary 
     to get the correct electric field 
     (see the EXPAND_CTU_GRID in unsplit_ctu.c).
   -------------------------------------------------------- */

  beg = MAX(beg, 3);
  end = MIN(end, grid->np_tot[g_dir] - 4);

  for (i = beg - 1; i <= end + 1; i++) {
    dp    = v[i + 1][PRS] - v[i - 1][PRS];
    min_p = MIN(v[i + 1][PRS], v[i - 1][PRS]);
    d2p   = v[i + 2][PRS] - v[i - 2][PRS];
    scrh = fabs(dp)/min_p;
    if (scrh < EPS2 || (v[i + 1][VXn] > v[i - 1][VXn])){
      f_t[i] = 0.0;
    }else{ 
      scrh   = OME2*(fabs(dp/d2p) - OME1);
      scrh   = MIN(1.0, scrh);
      f_t[i] = MAX(0.0, scrh);
    }
  }

  for (i = beg; i <= end; i++) {
    sj = (v[i + 1][PRS] < v[i - 1][PRS] ?  1 : -1);
    fj = MAX(f_t[i], f_t[i + sj]);
    for (nv = 0; nv < NVAR; nv++){
      vf   = v[i][nv]*fj;
      scrh = 1.0 - fj;
      vm[i][nv] = vf + vm[i][nv]*scrh;
      vp[i][nv] = vf + vp[i][nv]*scrh;
    }
  }
}
#undef   EPS2
#undef   OME1
#undef   OME2
#else
{

}
#endif

