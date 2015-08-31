#include"pluto.h"

/* ********************************************************************** */
void Flatten (const State_1D *state, int beg, int end, Grid *grid)
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
  real **v, **vp, **vm;
  static real *f_t;
   
  #if EOS == ISOTHERMAL 
   int PRS = RHO;
  #endif
   
  if (f_t == NULL){
    f_t   = ARRAY_1D(NMAX_POINT, double);    
  }

  v  = state->v;
  vp = state->vp;
  vm = state->vm;

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
  end = MIN(end, grid[g_dir].np_tot - 4);

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

/*
  int   beg, end, i, nv, sj;
  real   scrh1, scrh2, scrh3;
  real **a, **ap, **am;
  static real  *f_t, *fj, *dp, *d2p, *min_p;
   
  #if EOS == ISOTHERMAL 
   int PR = DN;
  #endif
     
  if (dp == NULL){
    f_t   = ARRAY_1D(NMAX_POINT, double);    
    fj    = ARRAY_1D(NMAX_POINT, double);
    dp    = ARRAY_1D(NMAX_POINT, double);
    d2p   = ARRAY_1D(NMAX_POINT, double);
    min_p = ARRAY_1D(NMAX_POINT, double);
  }

  a  = state->v;
  ap = state->vp;
  am = state->vm;

  beg = grid[g_dir].lbeg - 1;
  end = grid[g_dir].lend + 1;

  beg = MAX(beg, 3);
  end = MIN(end, grid[g_dir].np_tot - 3);

  for (i = beg - 2; i <= end + 2; i++) {
    dp[i]    = a[i + 1][PRS] - a[i - 1][PRS];
    min_p[i] = MIN(a[i + 1][PRS], a[i - 1][PRS]);
  }
  for (i = beg - 1; i <= end + 1; i++) {
    d2p[i]   = a[i + 2][PRS] - a[i - 2][PRS];
  }

  for (i = beg - 1; i <= end + 1; i++) {
    scrh1 = fabs(dp[i]) / min_p[i];
    scrh2 = a[i + 1][VXn] - a[i - 1][VXn];
    if (scrh1 < EPS2 || scrh2 > 0.0){
      f_t[i] = 0.0;
    }else{ 
      scrh3  = OME2*(fabs(dp[i]/d2p[i]) - OME1);
      scrh3  = MIN(1.0, scrh3);
      f_t[i] = MAX(0.0, scrh3);
    }
  }

  for (i = beg; i <= end; i++) {
    sj = (dp[i] < 0.0 ?  1 : -1);
    fj[i] = MAX(f_t[i], f_t[i + sj]);
  }

  for (i = beg; i <= end; i++) {
  for (nv = 0; nv < NVAR; nv++){
    scrh1 = a[i][nv]*fj[i];
    scrh2 = 1.0 - fj[i];
    am[i][nv] = scrh1 + am[i][nv]*scrh2;
    ap[i][nv] = scrh1 + ap[i][nv]*scrh2;
  }}
*/

}
#undef   EPS2
#undef   OME1
#undef   OME2
#else
{

}
#endif

