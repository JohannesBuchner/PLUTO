#include "pluto.h"

#define X 0.711
#define Y 0.2741

#define A_Z      30.0    /*   mean atomic weight of heavy elements  */
#define A_He     4.004   /*   atomic weight of Helium  */
#define A_H      1.008   /*   atomic weight of Hydrogen  */

/* ********************************************************************* */
void Radiat (double *v, double *rhs)
/*!
 *  Cooling for neutral or singly ionized gas: good up to about 35,000 K
 *  in equilibrium or shocks in neutral gas up to about 80 km/s.
 *  Assumed abundances in ab
 *  Uses t : Kelvin
 *       dene :  electron density cm*-3
 *       fneut : hydrogen neutral fraction    (adimensionale)
 *       ci,cr : H ionization and recombination rate coefficients
 *
 *
 * em(1) = TOTAL EMISSIVITY : (ergs cm**3 s**-1)
 * em(2) = Ly alpha  + two photon continuum: Aggarwal MNRAS 202,
 *         10**4.3 K
 * em(3) = H alpha: Aggarwal, Case B
 * em(4) = He I 584 + two photon + 623 (all n=2 excitations): Berrington
 *         &Kingston,JPB 20
 * em(5) = C I 9850 + 9823: Mendoza, IAU 103, 5000 K
 * em(6) = C II, 156 micron: Mendoza, 10,000 K
 * em(7) = C II] 2325 A: Mendoza, 15,000 K
 * em(8) = N I 5200 A: Mendoza, 7500 K
 * em(9) = N II 6584 + 6548 A: Mendoza
 * em(10) = O I 63 micron: Mendoza,2500 K
 * em(11) = O I 6300 A + 6363 A: Mendoza, 7500 K
 * em(12) = O II 3727: Mendoza
 * em(13) = Mg II 2800: Mendoza
 * em(14) = Si II 35 micron: Dufton&Kingston, MNRAS 248
 * em(15) = S II 6717+6727: Mendoza
 * em(16) = Fe II 25 micron: Nussbaumer&Storey
 * em(17) = Fe II 1.6 micron
 * em(18) = thermal energy lost by ionization
 * em(19) = thermal energy lost by recombination (2/3 kT per
 *          recombination.
 *          The ionization energy lost is not included here.
 *
 *********************************************************************** */
{

  int   ii, k;
  double  T, mu, st, rho, pr, fn;
  double  N_H, n_el, cr, ci, rlosst, em[20], src_pr;

  static double t00[18] = {0.0   , 0.0   , 1.18e5, 1.40e5, 2.46e5, 1.46e4,
                           92.1  , 6.18e4, 2.76e4, 2.19e4, 228.0 , 2.28e4,
                           3.86e4, 5.13e4, 410.0 , 2.13e4, 575.0 , 8980.0};

  static double  ep[18] = {0.0   , 0.0 , 10.2  , 1.89, 21.2  , 1.26,
                           0.0079, 5.33, 2.38  , 1.89, 0.0197, 1.96,
                           3.33  , 4.43, 0.0354, 1.85, 0.0495, 0.775};

  static double critn[18] = {0.0  , 0.0   , 1.e10, 1.e10, 1.e10,  312.0,
                             0.849, 1.93e7, 124.0, 865.0, 1090.0, 3950.0,
                             177.0, 1.e10 ,  16.8,  96.0,  580.0, 1130.0};

  static double om[18] = {0.0 , 0.0 , 0.90, 0.35, 0.15  , 0.067,
                          0.63, 0.52, 0.90, 0.30, 0.0055, 0.19,
                          0.33, 8.0 , 2.85, 1.75, 0.3   ,0.39};

  static double ab[18] = {0.0   , 0.0    , 1.0    , 1.0    , 0.1    , 0.0003,
                          0.0003, 0.0003 , 0.0001 , 0.0001 , 0.0006 , 0.0006,
                          0.0006, 0.00002, 0.00004, 0.00004, 0.00004, 0.00004};

  static double fn1[20] = {0.0, 0.0, 0.0, 0.0, 0.0,
                           0.1, 1.0, 1.0, 0.0, 1.0,
                           0.0, 0.0, 1.0, 1.0, 1.0,
                           1.0, 1.0, 1.0, 0.0, 1.0};

  static double fn2[20] = {0.0, 0.0, 1.0, 1.0, 1.0,
                           0.0, 0.0, 0.0, 1.0,-1.0,
                           1.0, 1.0,-1.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 1.0,-1.0};
  static int first_call = 1;
  static double E_cost, Unit_Time, N_H_rho, frac_Z, frac_He;

  if (first_call) {

    E_cost    = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);
    Unit_Time = UNIT_LENGTH/UNIT_VELOCITY;

  /*  ------------------------------------------------------
       conversion factor from total density to hydrogen
        number density, i.e.    nH = N_H_rho * rho
      ------------------------------------------------------   */

    N_H_rho  = (UNIT_DENSITY/CONST_amu)*(X/A_H);  /* n(H)  */
    frac_He  = (Y/A_He)*(A_H/X);
    frac_Z   = ((1 - X - Y)/A_Z)*(A_H/X);
    
     first_call = 0;
  }

  rho = v[RHO];
  pr  = v[PRS];

/* -----------------------------------
    Force fneut to stay between [0,1]
   ----------------------------------- */
  
  v[X_HI] = MAX(v[X_HI], 0.0);
  v[X_HI] = MIN(v[X_HI], 1.0);

  if (v[PRS] < 0.0) v[PRS] = g_smallPressure;

  fn  = v[X_HI];
  mu  = MeanMolecularWeight(v); 
  T   = pr/rho*KELVIN*mu;

 /*  --------------------------------------------------------
       Set source terms equal to zero when the temperature
       falls below g_minCoolingTemp 
    -------------------------------------------------------- */

/*    if (T < g_minCoolingTemp){
      a = b = src_pr = 0.0;
      continue;
    }
*/
  if (mu < 0.0){
    printf ("Negative mu in radiat \n");
    exit(1);
  }

  st    = sqrt(T);
  N_H   = N_H_rho*rho;  /* -- number density of hydrogen N_H = N(HI) + N(HII)  -- */

  /*  ----  coeff di ionizz. e ricomb. della particella fluida i,j  ----  */

  cr = 2.6e-11/st;
  ci = 1.08e-8*st*exp(-157890.0/T)/(13.6*13.6);

  n_el = N_H*(1.0 - fn + frac_Z);  /* -- electron number density, in cm^{-3} -- */
  rhs[X_HI] = Unit_Time*n_el*(-(ci + cr)*fn + cr); 

  em[1] = 0.0;
  for (k = 2; k <= 17; k++){
    em[k] = 1.6e-12*8.63e-6*om[k]*ep[k]*exp(-t00[k]/T)/st;
    em[k] = em[k]*critn[k]*st/(n_el + critn[k]*st);
    em[k] = em[k]*ab[k]*(fn1[k] + fn*fn2[k]);
    em[1] = em[1] + em[k];
  }

  em[18] = ci*13.6*1.6e-12*fn;
  em[19] = cr*0.67*1.6e-12*(1.0 - fn)*T/11590.0;
  em[1]  = em[1] + em[18] + em[19];

  /* ---------------------------------------------------
      rlosst is the energy loss in units of erg/cm^3/s;
      it must be multiplied by cost_E in order to match 
      non-dimensional units.
      Source term for the neutral fraction scales with 
      UNIT_TIME
     --------------------------------------------------- */
  double g_gamma;
#if EOS == IDEAL 
  g_gamma = 1.4; 
#elif EOS == PVTE_LAW
  g_gamma = Gamma1(v); 
#endif

  rlosst  =  em[1]*n_el*N_H;
  rhs[PRS] = -E_cost*rlosst*(g_gamma - 1.0);
  rhs[PRS] *= 1.0/(1.0 + exp( -(T - g_minCoolingTemp)/100.0)); /* -- cutoff -- */

} 

    
/* ********************************************************************* */
double MeanMolecularWeight (double *V)
/*
 *
 *
 * PURPOSE
 * 
 *   Compute the mean molecular weight as function of the 
 *   composition of the gas.
 *   The definitiion of the mean molecular weight \mu is 
 *   the standard one:
 *
 *     1     \sum_k f_k n_k
 *    --- = ----------------     (Clayton, pag 82-83)
 *    \mu    \sum_k f_k A_k
 * 
 *   where 
 *
 *    f_k   : is the fractional abundance (by number) with
 *            respect to hydrogen, f_k = N_k/N_H
 *
 *    A_K   : is the atomic weight
 *
 *    n_k   : is the number of free particles 
 *            contributed to the gas by element k
 *
 *   The mean molecular weight satifies 
 *
 *               \rho = \mu m_{amu} N_{tot}
 *   
 *   where N_{tot} is the total number of particles
 *
 *   For the ``Raymond'' cooling module \mu is calculated
 *   as follows:
 *
 *            A_H + f_He*A_He + f_Z*A_z
 *    \mu =  ---------------------------
 *             2 - fn + f_He + 2*f_Z
 *
 * 
 * ARGUMENTS
 *
 *   V:   a set of primitive variables
 *
 *********************************************************************** */
{
  int nv;

  for (nv = NFLX;nv <(NFLX+NIONS);nv++){
    if (V[nv] < 0.0) V[nv] = 0.0;
    if (V[nv] > 1.0) V[nv] = 1.0;
  }

  double N_H  = V[RHO]*(UNIT_DENSITY/CONST_amu)*(X/A_H);  /* n(H) */
  double N_He = V[RHO]*(UNIT_DENSITY/CONST_amu)*(Y/A_He); /* n(He) */
  double N_Z  = V[RHO]*(UNIT_DENSITY/CONST_amu)*((1.0-X-Y)/A_Z);  /*/ n(metals) */

  double fracHe = N_He/N_H;
  double fracZ = N_Z/N_H;
  
  double fn = V[X_HI];
 
  double munum = 1.0 + A_He*(fracHe) + A_Z*(fracZ);
  double muden = 2 - fn + fracHe + fracZ + 0.5*A_Z*(fracZ);

  return munum/muden;


  /*
  return  ( (A_H + frac_He*A_He + frac_Z*A_Z) /
         (2.0 + frac_He + 2.0*frac_Z - V[X_HI]));
  */
}

/* /\* ********************************************************************* *\/ */
/* double H_MassFrac (void) */
/* /\* */
/*  * */
/*  * */
/*  * PURPOSE */
/*  *  */
/*  *   Compute the mass fraction X of Hydrogen as function of the  */
/*  *   composition of the gas. */
/*  * */
/*  *             f_H A_H */
/*  *    X = ----------------    */
/*  *         \sum_k f_k A_k */
/*  *  */
/*  *   where  */
/*  * */
/*  *    f_k   : is the fractional abundance (by number),  */
/*  *               f_k = N_k/N_tot */
/*  *            of atomic species (no matter ionization degree). */
/*  * */
/*  *    A_K   : is the atomic weight */
/*  * */
/*  *    Note: In this module, f_H = 1.0 */
/*  * */
/*  *   where N_{tot} is the total number density of particles */
/*  * */
/*  * ARGUMENTS */
/*  * */
/*  *   none */
/*  * */
/*  *********************************************************************** *\/ */
/* { */
/*   return A_H/(A_H + frac_He*A_He + frac_Z*A_Z); */
/* } */

/* /\* ********************************************************************* *\/ */
/* double CompEquil (double N, double T, double *v) */
/* /\* */
/*  * */
/*  *      compute the equilibrium ionization balance for (rho,T)    */
/*  * */
/*  * */
/*  *********************************************************************** *\/ */
/* { */
/*   double cr, ci, st, n_e; */
/*   st = sqrt(T); */
/*   cr = 2.6e-11/st;  /\* recombination / ionization coefficients  *\/ */
/*   ci = 1.08e-8*st*exp(-157890.0/T)/(13.6*13.6); */
/*   v[X_HI] = cr / (ci + cr);   /\* compute fraction of neutrals at equilibrium *\/ */
/*   n_e = (1.0 - v[X_HI] + frac_Z)*N; */
/*   return n_e; */
/* } */

