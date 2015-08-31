#include "pluto.h"

/* **************************************************************** */
void SoundSpeed2 (double **u, double *cs2, double *h, int beg, int end,
                  int pos, Grid *grid)
/*
 *
 *    Define the square of the sound speed for different EOS
 *
 ****************************************************************** */
{
  int  i;

  #if EOS == IDEAL
   for (i = beg; i <= end; i++) cs2[i] = g_gamma*u[i][PRS]/u[i][RHO];
  #elif EOS == ISOTHERMAL
  {
    int    j,k;  /* -- used as multidimensional indices -- */
    double *x1, *x2, *x3;

    x1 = grid[IDIR].x;
    x2 = grid[JDIR].x;
    x3 = grid[KDIR].x;

    i = g_i;
    j = g_j;
    k = g_k;

    if (g_dir == IDIR) {
      double R;

      x1 = (pos == FACE_CENTER ? grid[IDIR].xr : grid[IDIR].x);
      for (i = beg; i <= end; i++){
        #if GEOMETRY == POLAR
         R = x1[i];
        #elif GEOMETRY == SPHERICAL
         R = x1[i]*sin(x2[j]);
        #endif
        cs2[i] = g_isoSoundSpeed*g_isoSoundSpeed/R;
      }

    }else if (g_dir == JDIR){
      double R;

      x2 = (pos == FACE_CENTER ? grid[JDIR].xr : grid[JDIR].x);
      for (j = beg; j <= end; j++) {
        #if GEOMETRY == POLAR
         R = x1[i];
        #elif GEOMETRY == SPHERICAL
         R = x1[i]*sin(x2[j]);
        #endif
        cs2[j] = g_isoSoundSpeed*g_isoSoundSpeed/R;
      }
    }else if (g_dir == KDIR){
      double R;

      x3 = (pos == FACE_CENTER ? grid[KDIR].xr : grid[KDIR].x);
      for (k = beg; k <= end; k++){
        #if GEOMETRY == POLAR
         R = x1[i];
        #elif GEOMETRY == SPHERICAL
         R = x1[i]*sin(x2[j]);
        #endif
        cs2[k] = g_isoSoundSpeed*g_isoSoundSpeed/R;
      }
    }
  }
  #else
   print ("! SoundSpeed2: not defined for this EoS\n");
   QUIT_PLUTO(1);
  #endif
}


/* *************************************************************** */
void Enthalpy (real **uprim, real *h, int beg, int end)
/*
 *
 *
 *
 ***************************************************************** */
{
  int i;
  double g_gammar;

  #if EOS == IDEAL
   g_gammar = g_gamma/(g_gamma - 1.0);
   for (i = beg; i <= end; i++){
     h[i] = g_gammar*uprim[i][PRS]/uprim[i][RHO];
   }
  #elif EOS == ISOTHERMAL 
   print (" Enthalpy not defined for isothermal EoS\n");
   QUIT_PLUTO(1);
  #endif
}
/* *************************************************************** */
void ENTROPY (real **v, real *s, int is, int ie)
/*
 *
 *
 *
 ***************************************************************** */
{
  int i;
  double rho;

  #if EOS == IDEAL
   for (i = is; i <= ie; i++){
     rho  = v[i][RHO];
     s[i] = v[i][PRS]/pow(rho,g_gamma);
   }
  #elif EOS == ISOTHERMAL || EOS == BAROTROPIC
   print (" Entropy not defined in isothermal or barotropic MHD\n");
   QUIT_PLUTO(1);
  #endif
}

