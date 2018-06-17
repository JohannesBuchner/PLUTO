#include "pluto.h"

/* ********************************************************************* */
void HallMHD_Flux (const Sweep *sweep, Data_Arr curlB, double **hall_flux,
                   double **lambdaH, int beg, int end, Grid *grid)
/*! 
 *  Compute Hall MHD fluxes for the induction and energy
 *  equations. Also, compute the diffusion coefficient.
 *  Called by either ParabolicFlux or ParabolicRHS.
 *
 * \param [in]   sweep      Pointer to the sweep structure
 * \param [out]  hall_flux  The hall mhd fluxes
 * \param [out]  lambdaH    The diffusion coefficients evaluated at 
 *                          cell interfaces
 * \param [in]   beg        initial index of computation
 * \param [in]   end        final   index of computation
 * \param [in]  grid        pointer to an array of Grid structures.
 *
 \author A. Mignone (mignone@ph.unito.it)
 \date   March 21, 2017
 *********************************************************************** */
{
  int  i, j, k, nv;
  double *x1, *x2, *x3, scrh;
  double eta, Jp[NVAR], v[NVAR]; 
  double ***Jx1, ***Jx2, ***Jx3;
  double *J, Vh, VhB, Bsq, Vhx, Vhy, Vhz;
  double JxB[3], **vc;

  /* !!!!!!!!
  1D Primitive Vector at the Left and Right Interface is needed 
  in vector v[NVAR]. 
  Should we use v = state->v[i] and provide state as Input to the function.
  Or should we use v = sweep->vn[i] and then interpolate in place.
   */

  #ifdef STAGGERED_MHD
   Jx1 = curlB[IDIR];
   Jx2 = curlB[JDIR]; 
   Jx3 = curlB[KDIR]; 
  #endif

  x1 = grid->x[IDIR];  
  x2 = grid->x[JDIR];  
  x3 = grid->x[KDIR]; 
  
  vc = sweep->vn; /* Cell Centered Variables to compute v */

/* ----------------------------------------------- 
     Compute resistive flux for the induction
     and energy equations
   ----------------------------------------------- */
 
  if (g_dir == IDIR){
  
    j  = g_j; k = g_k;    
    x1 = grid->xr[IDIR];
    /*
    You do not require here to compute Whistler Speed! 
    HallMHD_WhistlerSpeed (stateL, beg, end, grid);
    HallMHD_WhistlerSpeed (stateR, beg, end, grid);
    Dcoeffs are obtained from Drift speed.! 
    */

    for (i = beg; i <= end; i++){
    /* -- interface values at (i+1/2, j, k) -- */
   
      for (nv = 0; nv < NVAR; nv++){
        v[nv] = vc[i][nv];
      }

      eta = HallMHD_ne(v);

      #ifdef STAGGERED_MHD
       EXPAND(                                 ;   ,
              Jp[KDIR] = AVERAGE_Y(Jx3,k,j-1,i);   ,
              Jp[JDIR] = AVERAGE_Z(Jx2,k-1,j,i);)

      #else
       for (nv = 0; nv < 3; nv++) J[nv] = curlB[nv][k][j][i];
      #endif
      
      JxB[JDIR] = -J[IDIR]*v[BX3] + J[KDIR]*v[BX1];
      JxB[KDIR] =  J[IDIR]*v[BX2] - J[JDIR]*v[BX1];     
      
      hall_flux[i][BX1] += 0.0;
      hall_flux[i][BX2] += eta*JxB[KDIR];
      hall_flux[i][BX3] += -eta*JxB[JDIR];


      #if HAVE_ENERGY
     
      Vh =  -J[g_dir]*eta; 
      Vhx = -J[IDIR]*eta;
      Vhy = -J[JDIR]*eta;  
      Vhz = -J[KDIR]*eta;

      VhB   = EXPAND(Vhx*v[BX1] , + Vhy*v[BX2], + Vhz*v[BX3]);
      Bsq = EXPAND(v[BX1]*v[BX1] , + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
      hall_flux[i][ENG] += 0.5*Vh*(Bsq) - (VhB)*v[BXn]; 
      #endif

    /* -- store diffusion coefficient -- */

      EXPAND(lambdaH[i][BX1] = 0.0; ,
	     lambdaH[i][BX2] = eta*J[JDIR]; , 
	     lambdaH[i][BX3] = eta*J[KDIR]; )
	     
    }

  }else if (g_dir == JDIR){
  
    i = g_i; k = g_k;    
    x2 = grid->xr[JDIR];
   
    for (j = beg; j <= end; j++){

    /* -- interface value -- */

      for (nv = 0; nv < NVAR; nv++) {
        v[nv] = vc[j][nv];
      }
      eta = HallMHD_ne(v);

      #ifdef STAGGERED_MHD
       EXPAND(Jp[KDIR] = AVERAGE_X(Jx3,k,j,i-1);   ,
                                                    ;   ,
              Jp[IDIR] = AVERAGE_Z(Jx1,k-1,j,i);)
                                                 
      #else
       for (nv = 0; nv < 3; nv++) Jp[nv] = curlB[nv][k][j][i];
      #endif

      JxB[IDIR] =  J[JDIR]*v[BX3] - J[KDIR]*v[BX2];
      JxB[KDIR] =  J[IDIR]*v[BX2] - J[JDIR]*v[BX1];
      
      hall_flux[j][BX1] += -eta*JxB[KDIR];
      hall_flux[j][BX2] += 0.0;
      hall_flux[j][BX3] += eta*JxB[IDIR];

      #if HAVE_ENERGY
     
      Vh =  -J[g_dir]*eta; 
      Vhx = -J[IDIR]*eta;
      Vhy = -J[JDIR]*eta;  
      Vhz = -J[KDIR]*eta;

      VhB   = EXPAND(Vhx*v[BX1] , + Vhy*v[BX2], + Vhz*v[BX3]);
      Bsq = EXPAND(v[BX1]*v[BX1] , + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
      hall_flux[j][ENG] += 0.5*Vh*(Bsq) - (VhB)*v[BXn]; 
      #endif

     
     
     /* -- store diffusion coefficient -- */
      EXPAND(lambdaH[i][BX1] = eta*J[IDIR]; ,
	     lambdaH[i][BX2] = 0.0; , 
	     lambdaH[i][BX3] = eta*J[KDIR]; )
    }

  }else if (g_dir == KDIR){
  
    i = g_i; j = g_j;    
    x3 = grid->xr[KDIR];

    for (k = beg; k <= end; k++){

    /* -- interface value -- */

      for (nv = 0; nv < NVAR; nv++) {
         v[nv] = vc[k][nv];
      }

      eta = HallMHD_ne(v);

      #ifdef STAGGERED_MHD
       Jp[JDIR] = AVERAGE_X(Jx2,k,j,i-1);
       Jp[IDIR] = AVERAGE_Y(Jx1,k,j-1,i);
      
      #else
       for (nv = 0; nv < 3; nv++) Jp[nv] = curlB[nv][k][j][i];
      #endif
      JxB[IDIR] =  J[JDIR]*v[BX3] - J[KDIR]*v[BX2];
      JxB[JDIR] = -J[IDIR]*v[BX3] + J[KDIR]*v[BX1];
      
      hall_flux[k][BX1] += eta*JxB[JDIR];
      hall_flux[k][BX2] += -eta*JxB[IDIR];
      hall_flux[k][BX3] += 0.0;
  
      
      #if HAVE_ENERGY
      Vh =  -J[g_dir]*eta; 
      Vhx = -J[IDIR]*eta;
      Vhy = -J[JDIR]*eta;  
      Vhz = -J[KDIR]*eta;

      VhB   = EXPAND(Vhx*v[BX1] , + Vhy*v[BX2], + Vhz*v[BX3]);
      Bsq = EXPAND(v[BX1]*v[BX1] , + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
      hall_flux[k][ENG] += 0.5*Vh*(Bsq) - (VhB)*v[BXn]; 
      #endif
 
   /* -- store diffusion coefficient -- */
      EXPAND(lambdaH[i][BX1] = eta*J[IDIR]; ,
	    lambdaH[i][BX2] = eta*J[JDIR]; , 
	    lambdaH[i][BX3] = 0.0; )

    }
  }
}

