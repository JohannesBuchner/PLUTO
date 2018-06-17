#include "pluto.h"
#include "cooling_defs.h"

#define MAX_NLEVEL 16

/* ********************************************************************* */
double CompEquil (double N, double T, double *v0)
/*!
 *  Compute the equilibrium ionization balance for (rho,T)   
 *
 *
 *********************************************************************** */
{
  int nlev, ind_start, i,j, k, nv, nrt;
  double d, n_el, n_el_old, sT;
  double v[NVAR];
  static double **M;
  static double *rhs;
  static int *p;
   
/* ----  start index in the ions matrix for each element ---- */

  int ind_starts[] = {X_HI, X_HeI, X_CI, X_NI, X_OI, X_NeI, X_SI, X_SI+S_IONS};

/* ----  number of ions for each element  ---- */

  int nlevs[]      = {1, 2, C_IONS, N_IONS, O_IONS, Ne_IONS, S_IONS, Fe_IONS};

/* -- Copy solution vector before any modification -- */

  NVAR_LOOP(nv) v[nv] = v0[nv];  

  sT = sqrt(T);
/*   N  = find_N(rho);   -- Total number density -- */
   
  if (M == NULL) {     /* -- Allocate matrices for the linear eq system -- */
    M   = ARRAY_2D(MAX_NLEVEL, MAX_NLEVEL, double);
    rhs = ARRAY_1D(MAX_NLEVEL, double);
    p   = ARRAY_1D(MAX_NLEVEL, int);
  }
  for (nv = 0; nv < 8; nv++) ind_starts[nv] -= NFLX; /* -- offset index to start from zero -- */
    
/* ---------------------------------------------------------------
     Initialize the relative abundances.
     Note: at the first iteration step below the charge-transfer
     terms will be 0, as the HI and HeII abundances are set to 0.
   --------------------------------------------------------------- */
    
  for (j = 0; j < NIONS; j++) v[NFLX+j] = 0.0;

/* --------------------------------------------------------------------------
     Initially all ion fractions are 0, so also n_el will be zero. As first
     guess for the iterations though, if T>5000K, let's use a higher value.
   -------------------------------------------------------------------------- */

  if (T > 5000.) n_el = N*0.1;
  else n_el = N*1.e-3;
 
  n_el_old = n_el/2.;
   
/* ----------------------------------------------------------------------
    Start the iteration for finding the equilibrium ionization balance. 
    Stop if n_el converges or is very small, or the number of iterations 
    becomes too large.
   ---------------------------------------------------------------------- */

  nrt = 0;
  while ( fabs(n_el - n_el_old) / n_el_old > 1.e-5  &&  nrt<500 && n_el > 1.0e-14) {
   
    n_el_old = n_el;     
 
    Find_Rates(T, n_el, N, v);
    for (k = 1; k < 8; k++) {     /* Main loop on elements */
      nlev = nlevs[k];            /* Number of ionization states for element k  */
      ind_start = ind_starts[k];  /* Start index in the ions matrix for element k */

      if (nlev == 0) continue; /* Skip element with zero ions */

   /* -----------------------------------------
            compute coefficient matrix
      -----------------------------------------    */

      for (i = 0; i < MAX_NLEVEL; i++) {
      for (j = 0; j < MAX_NLEVEL; j++) {
        M[i][j] = 0.0;
      }}
      for (j = 0 ; j < nlev-1 ; j++) {
        M[j][j+1] +=  CoolCoeffs.Rrate[ind_start+j];
        M[j][j]   += -CoolCoeffs.Crate[ind_start+j];
        if (j>0) M[j][j-1] += CoolCoeffs.Lrate[ind_start+j];
      }

  /*  -------------------------------------
       Replace 1st eq. with normalization
       condition and define rhs
      -------------------------------------  */

      for (j = 0; j < nlev; j++) {
        M[nlev-1][j] = 1.0;
        rhs[j] = 0.0;
      }
      rhs[nlev-1] = 1.0;  

     /* ----------------------------------
           Solve equations by LU decomp
        ---------------------------------- */
    
      LUDecompose (M, nlev, p, &d);
      LUBackSubst (M, nlev, p, rhs); 

      for (i = 0 ; i < nlev ; i++) {
        if (rhs[i] >= 0. && rhs[i] <= 1.0) v[NFLX+i+ind_start] = rhs[i];
        else if (rhs[i] < 0.0) v[NFLX + i + ind_start] = 0.0;  
        else v[NFLX + i + ind_start] = 1.0;
      }
   
    }  /* end main loop on elements  */

  /* ------------------------------------------------------
         ... and now compute the ionization balance for H.
     ------------------------------------------------------  */
 
    v[X_HI] = CoolCoeffs.Rrate[0] / (CoolCoeffs.Crate[0] + CoolCoeffs.Rrate[0]);
     
    if (v[X_HI] > 1.0){
      print ("\n!CompEquil: GR(H) = %12.6e\nDR(H) = %12.6e\n\n", 
              CoolCoeffs.Rrate[0], CoolCoeffs.Crate[0]);
      QUIT_PLUTO(1);
    }

  /* -----------------------------------
       compute electron number density 
     ----------------------------------- */

    n_el = (1.0 - v[X_HI])*elem_ab[el_H]*N; /* contribution from ionized H  */
    for (nv = 0; nv < NIONS; nv++) {
       n_el += v[NFLX + nv]*(rad_rec_z[nv] - 1.)*elem_ab[elem_part[nv]]*N;    
    }
     
    if (n_el != n_el) {
      print ("! CompEquil: error!! n_el NaN\n\n");
      QUIT_PLUTO(1);
    }             
    nrt++;  /* increment the iteration counter  */
  }  /* --  end main iteration  -- */
 
/* -- Update input solution vector -- */

  for (nv = 0; nv < NIONS; nv++) v0[NFLX+nv] = v[NFLX+nv];
  
  return (n_el);
}
