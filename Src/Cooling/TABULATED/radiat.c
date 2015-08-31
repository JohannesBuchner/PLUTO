#include "pluto.h"
/* ***************************************************************** */
void Radiat (double *v, double *rhs)
/*!
 *   Provide r.h.s. for tabulated cooling.
 * 
 ******************************************************************* */
{
  int    klo, khi, kmid;
  static int ntab;
  double  mu, T, Tmid, scrh, dT, prs;
  static double *L_tab, *T_tab, E_cost;
  
  FILE *fcool;

/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */

  if (T_tab == NULL){
    print1 (" > Reading table from disk...\n");
    fcool = fopen("cooltable.dat","r");
    if (fcool == NULL){
      print1 ("! Radiat: cooltable.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab = ARRAY_1D(20000, double);
    T_tab = ARRAY_1D(20000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", T_tab + ntab, 
                                       L_tab + ntab)!=EOF) {
      ntab++;
    }
    E_cost = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);
  }

/* ---------------------------------------------
            Get pressure and temperature 
   --------------------------------------------- */

  prs = v[RHOE]*(g_gamma-1.0);
  if (prs < 0.0) {
    prs     = g_smallPressure;
    v[RHOE] = prs/(g_gamma - 1.0);
  }

  mu  = MeanMolecularWeight(v);
  T   = prs/v[RHO]*KELVIN*mu;

  if (T != T){
    printf (" ! Nan found in radiat \n");
    printf (" ! rho = %12.6e, prs = %12.6e\n",v[RHO], prs);
    QUIT_PLUTO(1);
  }

  if (T < g_minCoolingTemp) { 
    rhs[RHOE] = 0.0;
    return;
  }

/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */

  klo = 0;
  khi = ntab - 1;

  if (T > T_tab[khi] || T < T_tab[klo]){
    print (" ! T out of range   %12.6e\n",T);
    QUIT_PLUTO(1);
  }

  while (klo != (khi - 1)){
    kmid = (klo + khi)/2;
    Tmid = T_tab[kmid];
    if (T <= Tmid){
      khi = kmid;
    }else if (T > Tmid){
      klo = kmid;
    }
  }

/* -----------------------------------------------
    Compute r.h.s
   ----------------------------------------------- */

  dT       = T_tab[khi] - T_tab[klo];
  scrh     = L_tab[klo]*(T_tab[khi] - T)/dT + L_tab[khi]*(T - T_tab[klo])/dT;
  rhs[RHOE] = -scrh*v[RHO]*v[RHO];
  
  scrh       = UNIT_DENSITY/(CONST_amu*mu);  
  rhs[RHOE] *= E_cost*scrh*scrh;
}
