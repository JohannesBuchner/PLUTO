#include "pluto.h"

#if BACKGROUND_FIELD == YES
/* ********************************************************************* */
double **GetBackgroundField (int beg, int end, int where, Grid *grid)
/*
 *
 *
 *
 *
 *
 *********************************************************************** */
{
  int i;
  double *x1, *x2, *x3;
  static double **bck_fld;

/* ----------------------------------------------------
         Check for incompatibilities 
   ---------------------------------------------------- */
  
  #if (TIME_STEPPING != RK2) && (TIME_STEPPING != RK3)
   print1 ("! Background field splitting works with RK integrators ONLY \n");
   QUIT_PLUTO(1);
  #endif
  #if DIVB_CONTROL == EIGHT_WAVES
   print1 ("! Background field splitting works with CT or GLM ONLY \n");
   QUIT_PLUTO(1);
  #elif DIVB_CONTROL == CONSTRAINED_TRANSPORT
   #if (CT_EMF_AVERAGE != ARITHMETIC) && (CT_EMF_AVERAGE != UCT_HLL)
    print1 ("! Background field splitting works with ARITHMETIC or");
    print1 (" UCT_HLL averages only\n");
    QUIT_PLUTO(1);
   #endif
  #endif

  if (bck_fld == NULL) {
    bck_fld = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

/* ----------------------------------
     get pointers to coordinates 
   ---------------------------------- */
   
  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  if (g_dir == IDIR){
    if (where == FACE_CENTER) x1 = grid[IDIR].xr;
    for (i = beg; i <= end; i++){
      BackgroundField (x1[i],x2[g_j],x3[g_k], bck_fld[i] + BX1);
    }
  }else if (g_dir == JDIR){
    if (where == FACE_CENTER) x2 = grid[JDIR].xr;
    for (i = beg; i <= end; i++){
      BackgroundField (x1[g_i],x2[i],x3[g_k], bck_fld[i] + BX1);
    }
  }else{
    if (where == FACE_CENTER) x3 = grid[KDIR].xr;
    for (i = beg; i <= end; i++){
      BackgroundField (x1[g_i],x2[g_j],x3[i], bck_fld[i] + BX1);
    }
  }    

  return(bck_fld);
}
#endif
