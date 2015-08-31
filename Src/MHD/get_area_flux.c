#include "pluto.h"

/* ********************************************************************* */
void GetAreaFlux (const State_1D *state, double **fA,
                  double **fvA, int beg, int end, Grid *grid)
/*
 *
 *********************************************************************** */
{
  int    i,j,k, nv;
  double A, rp, rp2, rp3, s, s2;
  double **flux = state->flux;
#if ENTROPY_SWITCH
  double **visc_flux = state->visc_flux; 
  double **visc_src  = state->visc_src; 
  double **tc_flux   = state->tc_flux; 
  double **res_flux  = state->res_flux;
#endif 

  if (g_dir == IDIR) { 

    for (i = beg - 1; i <= end; i++){ 
      A  = grid[IDIR].A[i];
      rp = grid[IDIR].xr[i];
#if GEOMETRY == CYLINDRICAL
      NVAR_LOOP(nv) fA[i][nv] = flux[i][nv]*A;
  #if COMPONENTS == 3    
      fA[i][MX3] = flux[i][MX3]*A*A;
      IF_DUST(fA[i][MX3_D] = flux[i][MX3_D]*A*A;)
  #endif    
#elif GEOMETRY == POLAR
      NVAR_LOOP(nv) fA[i][nv] = flux[i][nv]*A;
  #if COMPONENTS >= 2
      fA[i][MX2] *= A;
      IF_DUST(fA[i][MX2_D] *= A;)
  #endif    
#elif GEOMETRY == SPHERICAL
      rp2 = rp*rp;
      NVAR_LOOP(nv) fA[i][nv] = flux[i][nv]*rp2;
  #if COMPONENTS == 3
      rp3 = rp2*rp;
      fA[i][MX3] = flux[i][MX3]*rp3;
      IF_DUST(fA[i][MX3_D] = flux[i][MX3_D]*rp3;)
  #endif
  #if PHYSICS == MHD
      EXPAND(                            ;   ,
             fA[i][BX2] = flux[i][BX2]*rp;  ,
             fA[i][BX3] = flux[i][BX3]*rp;)
  #endif
#endif

      #if (ENTROPY_SWITCH) && (VISCOSITY == EXPLICIT)
       EXPAND(fvA[i][MX1] = visc_flux[i][MX1]*A; ,
              fvA[i][MX2] = visc_flux[i][MX2]*A; ,
              fvA[i][MX3] = visc_flux[i][MX3]*A;)
       #if GEOMETRY == CYLINDRICAL && COMPONENTS == 3
        fvA[i][MX3] *= rp;
       #elif GEOMETRY == POLAR && COMPONENTS >= 2
        fvA[i][MX2] *= rp; 
       #elif GEOMETRY == SPHERICAL && COMPONENTS == 3
        fvA[i][MX3] *= rp; 
       #endif
       #if HAVE_ENERGY
        fvA[i][ENG] = visc_flux[i][ENG]*A;
       #endif
      #endif
    }

  }else if (g_dir == JDIR){

    for (j = beg - 1; j <= end; j++){

#if GEOMETRY == SPHERICAL  
      s  = grid[JDIR].A[j];
      s2 = s*s; 

      NVAR_LOOP(nv) fA[j][nv] = flux[j][nv]*s;
  #if COMPONENTS == 3
      fA[j][MX3] = flux[j][MX3]*s2;
      IF_DUST(fA[j][MX3_D] = flux[j][MX3_D]*s2;)
  #endif

      #if (ENTROPY_SWITCH) && (VISCOSITY == EXPLICIT)
       EXPAND(fvA[j][iMR]   = visc_flux[j][iMR]*s;    ,
              fvA[j][iMTH]  = visc_flux[j][iMTH]*s;   ,
              fvA[j][iMPHI] = visc_flux[j][iMPHI]*s2;)
       #if HAVE_ENERGY
        fvA[j][ENG] = visc_flux[j][ENG]*s;
       #endif
      #endif
#endif

    }

  }else if (g_dir == KDIR){

    print1 ("! GetAreaFlux(): should not be called during x3 sweep\n");
    QUIT_PLUTO(1);

  }

}
  