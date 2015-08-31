#include "pluto.h"

#define POW10(x)  (pow(10.0, x))
#define TABLE_NPT   2048

/* ********************************************************************* */
void H2RateTables(double T, double *krvals)
/*!
 *  On first call, compute and store rate coefficient arrays 
 *  needed for the H2_COOL module.
 *  On subsequent calls, use linear interpolation between adjacent
 *  tabulated values to obtain the coefficients at the desired 
 *  temperature T. 
 *********************************************************************** */
{
  int i, indx_lo, indx_hi;
  double st, t3, tev, lnT, lnTmin, lnTmax, dlnT;
  double scrh, dum_lo, dum_hi;
  double Tmin, Tmax;
  static double *lnTarr, *Rate_arr1, *Rate_arr2, *Rate_arr3, *Rate_arr4;
  static double *Rate_arr5, *Rate_arr6;
  static double *Emiss_arr1, *Emiss_arr2, *Emiss_arr3;
  double Tdust = 15.0;

  Tmin = 1.0e-2;
  Tmax = 1.0e8;
  lnTmin = log10(Tmin);
  lnTmax = log10(Tmax);
  dlnT   = (lnTmax - lnTmin)/((double)TABLE_NPT - 1.0);

  T = MAX(T, Tmin);  /* We constrain local T to lie between [Tmin, Tmax] so */
  T = MIN(T, Tmax);  /* that coefficients will be constant for values of
                        temperature overflow / underflow.  */

  if (lnTarr == NULL){
    lnTarr     = ARRAY_1D(TABLE_NPT, double);
    Rate_arr1  = ARRAY_1D(TABLE_NPT, double);
    Rate_arr2  = ARRAY_1D(TABLE_NPT, double);
    Rate_arr3  = ARRAY_1D(TABLE_NPT, double);
    Rate_arr4  = ARRAY_1D(TABLE_NPT, double);
    Rate_arr5  = ARRAY_1D(TABLE_NPT, double);
    Rate_arr6  = ARRAY_1D(TABLE_NPT, double);
    Emiss_arr1 = ARRAY_1D(TABLE_NPT, double);
    Emiss_arr2 = ARRAY_1D(TABLE_NPT, double);
    Emiss_arr3 = ARRAY_1D(TABLE_NPT, double);
    
    for(i = 0 ; i < TABLE_NPT; i++){
      lnTarr[i] = lnTmin + i*dlnT;
      T     = POW10(lnTarr[i]);
      st    = sqrt(T);
      t3    = T/1.e3;  
      tev   = T*CONST_kB/CONST_eV;
      
      Rate_arr1[i] = 3.e-18*st/(1.0 + 0.04*st + 2.e-3*T + 8.e-6*T*T); /* fa ~ 1 and Tg << T  */
      Rate_arr2[i] = 1.067e-10*pow(tev,2.012)*exp(-(4.463/tev)*pow((1.0 + 0.2472),3.512));
      Rate_arr3[i] = 1.e-8*exp(-84100.0/T); 
      Rate_arr4[i] = 4.4e-10*pow(T,0.35)*exp(-102000.0/T);
      Rate_arr5[i] = 2.6e-11/st;
      Rate_arr6[i] = 1.08e-8*st*exp(-157890.0/T)/(13.6*13.6);
      Emiss_arr1[i] = 6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3);
      Emiss_arr2[i] = (9.5e-22*pow(t3,3.76)/(1. + 0.12*pow(t3,2.1)))*exp(pow((-0.13/t3),3.0))
	+ 3.0e-24*exp(-0.51/t3); 
      Emiss_arr3[i] = 3.8e-33*st*(T-Tdust)*(1.0 - 0.8*exp(-75.0/T));
    }
  }else{   /* Perform linear interpolation */
    lnT = log10(T);

    scrh    = (lnT - lnTmin)/dlnT;
    indx_lo = INT_FLOOR(scrh);
    indx_hi = indx_lo + 1;

    dum_lo  = (lnTarr[indx_hi] - lnT)/dlnT;
    dum_hi  = (lnT - lnTarr[indx_lo])/dlnT;

    if (indx_lo < 0){
      print ("! Make_RateTable: Cannot Locate Index for T = %12.6e\n",T);
      QUIT_PLUTO(1);
    }else{
      krvals[0] = Rate_arr1[indx_lo]*dum_lo + Rate_arr1[indx_hi]*dum_hi; /* kr1 */
      krvals[1] = Rate_arr2[indx_lo]*dum_lo + Rate_arr2[indx_hi]*dum_hi; /* kr2 */
      krvals[2] = Rate_arr3[indx_lo]*dum_lo + Rate_arr3[indx_hi]*dum_hi; /* kr3 */
      krvals[3] = Rate_arr4[indx_lo]*dum_lo + Rate_arr4[indx_hi]*dum_hi; /* kr4 */
      krvals[4] = Rate_arr5[indx_lo]*dum_lo + Rate_arr5[indx_hi]*dum_hi; /* cr */
      krvals[5] = Rate_arr6[indx_lo]*dum_lo + Rate_arr6[indx_hi]*dum_hi; /* ci */
      krvals[6] = Emiss_arr1[indx_lo]*dum_lo + Emiss_arr1[indx_hi]*dum_hi; /* em[1] */
      krvals[7] = Emiss_arr2[indx_lo]*dum_lo + Emiss_arr2[indx_hi]*dum_hi; /* em[2] */
      krvals[8] = Emiss_arr3[indx_lo]*dum_lo + Emiss_arr3[indx_hi]*dum_hi; /* em[11] */
    }
  }
  return;
}

/* ***************************************************************** */
void Radiat (double *v, double *rhs)
/*
 *  Cooling for molecular+atomic+ionized hydrogen gas: 
 * 
 * em(7) = Cooling due to H2 Rotational Vibration (Glover, 2008).
 * em(8) = Cooling due to Hydrogen Ionization. 
 * em(9) = Cooling due to recombination
 * em(10) = Cooling due to H2 Dissociation
 * em(11) = Gas - Grain cooling in molecular hydrogen.
 * em(13) = Total Cooling Function from all of the above.  
 *
 ******************************************************************* */
{
  int   ii, k, nv, status;
  double  T, mu, t3, rho, prs, fn, gn, hn, x;
  double  N_H, n_el, cr, ci, rlosst, src_pr, crold,crnew;
  double kr1, kr2, kr3, kr4, em[20], krvalues[9];
  
  static int first_call = 1;
  static real E_cost, Unit_Time, N_H_rho, frac_He, frac_Z, N_tot;

  if (first_call) {
    E_cost    = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);
    Unit_Time = UNIT_LENGTH/UNIT_VELOCITY;
        
    N_H_rho  = (UNIT_DENSITY/CONST_amu)*(H_MASS_FRAC/CONST_AH);
    frac_He  = (He_MASS_FRAC/CONST_AHe)*(CONST_AH/H_MASS_FRAC);
    frac_Z   = ((1 - H_MASS_FRAC - He_MASS_FRAC)/CONST_AZ)*(CONST_AH/H_MASS_FRAC);
    N_tot  =  N_H_rho*(1.0 + frac_He + frac_Z);
    H2RateTables(100.0, krvalues);
    first_call = 0;
  }

/* ---------------------------------------------
    Force fneut and fmol to stay between [0,1]
   --------------------------------------------- */

  NIONS_LOOP(nv){
    v[nv] = MAX(v[nv], 0.0);
    v[nv] = MIN(v[nv], 1.0);
  }
  v[X_H2] = MIN(v[X_H2], 0.5);

/* ---------------------------------------------------
    Compute temperature from rho, rhoe and fractions
   --------------------------------------------------- */
   
  rho = v[RHO];
  mu  = MeanMolecularWeight(v); 
  if (mu < 0.0){
    print1 ("! Radiat: mu = %f < 0 \n",mu);
    QUIT_PLUTO(1);
  }
  #if EOS == IDEAL
   if (v[RHOE] < 0.0) v[RHOE] = g_smallPressure/(g_gamma-1.0);
   prs = v[RHOE]*(g_gamma-1.0);
   T   = prs/rho*KELVIN*mu;
  #else
   status = GetEV_Temperature(v[RHOE], v, &T);
   if (status != 0) {
     T = T_CUT_RHOE;
     v[RHOE] = InternalEnergy(v, T);
   }
  #endif

  N_H = N_H_rho*rho;  /* Hydrogen number density
                         N_H = N(X_HI) + 2N(X_H2)+  N(X_HII)  */
  fn  = v[X_HI];
  gn  = v[X_H2];
  hn  = v[X_HII];

/* Recombination and ionization coefficients */ 

  n_el = N_H*(hn + 0.5*CONST_AZ*frac_Z);  /* -- electron number density, in cm^{-3} -- */
  
  if (n_el < 0){
    print1 ("! Radiat: negative electron density\n");
    QUIT_PLUTO(1);
  }

  t3    = T/1.e3;  
 
/* The Units of kr1, kr2, kr3, kr4, cr, ci = cm^{3}s^{-1}*/
  
  H2RateTables(T, krvalues);
  kr1 = krvalues[0];
  kr2 = krvalues[1];
  kr3 = krvalues[2];
  kr4 = krvalues[3];
  cr  = krvalues[4];
  ci  = krvalues[5];

  /*
  kr1 = 3.e-18*st/(1.0 + 0.04*st + 2.e-3*T + 8.e-6*T*T);
  kr2 = 1.067e-10*pow(tev,2.012)*exp(-(4.463/tev)*pow((1.0 + 0.2472),3.512));
  kr3 = 1.e-8*exp(-84100.0/T); 
  kr4 = 4.4e-10*pow(T,0.35)*exp(-102000.0/T);
 */
/* ************ ENSURING SUM IS UNITY ********************** 
  * This is done to ensure that the sum of speicies 
  * equal to 1. Total Hydrogen number N_H composed of
  * N_H = nHI + 2.*nH2 + nHII. 
  * while fraction gn = nH2/N_H, thus the rhs[X_H2] has 
  * a pre-factor of 0.5. 
  * And the fact the sum = 1 is ensured if the
  * corresponding sum of RHS = 0. 
  * i.e., rhs[X_H1] + 2.0*rhs[X_H2] + rhs[X_HII] = 0.0.
  ********************************************************** */

  x = n_el/N_H;
  rhs[X_HI] = Unit_Time*N_H*(  gn*(kr2*fn + kr3*gn + kr4*x) + cr*hn*x 
                              -fn*(ci*x + kr1*fn));
  rhs[X_H2]  = 0.5*Unit_Time*N_H*(kr1*fn*fn - gn*(kr2*fn + kr3*gn + kr4*x));
  rhs[X_HII] = Unit_Time*N_H*(ci*fn - cr*hn)*x;
/*
  double sum;
  sum = 0.0;
  sum=fabs(rhs[X_HI] + 2.0*rhs[X_H2] + rhs[X_HII]);
  if (sum > 1.e-9){
    print("%le > Sum(RHS) != 0 for H!!\n",sum);
    QUIT_PLUTO(1);
  }
*/

  double logt3 = log(t3);

/* The unit of cooling function Lambda in erg cm^{-3} s^{-1}*/
  
  /* em[1] = 6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3);
  em[2] = (9.5e-22*pow(t3,3.76)/(1. + 0.12*pow(t3,2.1)))*exp(pow((-0.13/t3),3.0))
           + 3.0e-24*exp(-0.51/t3);*/
  

  em[1] = krvalues[6];
  em[2] = krvalues[7];
  em[3] = em[1] + em[2];

  if (t3 < 0.1){ /* Following the Table 8 from Glover & Abel, MNRAS (2008) */
    em[4] = -16.818342 + logt3*(37.383713 + logt3*(58.145166 + logt3*(48.656103 + logt3*(20.159831 + logt3*(3.8479610)))));
  }else if(t3 < 1.0){
    em[4] = -24.311209 + logt3*(3.5692468 + logt3*(-11.332860 + logt3*(-27.850082 + logt3*(-21.328264 + logt3*(-4.2519023)))));

  }else{
    em[4] = -24.311209 + logt3*(4.6450521 + logt3*(-3.7209846 + logt3*(5.9369081 + logt3*(-5.5108047 + logt3*(1.5538288))))); 

  }
  em[4] = POW10(em[4]);

  em[5] = -23.962112 + logt3*(2.09433740 + logt3*(-0.77151436 + logt3*(0.43693353 + logt3*(-0.14913216 + logt3*(-0.033638326))))); 

  em[5] = POW10(em[5]);
  em[6] = (fn*em[4] + gn*em[5])*N_H; 
  
/* Cooling due to X_H2 rotational vibration.  */

  em[7] = gn*N_H*em[3]/(1.0 + (em[3]/em[6])); 

/* Cooling due to ionization */
  
  em[8] = ci*13.6*1.6e-12*fn;           
  
/* Cooling due to radiative recombination */

  em[9] = cr*0.67*1.6e-12*hn*T/11590.0; 
  
/* Cooling due to H2 disassociation */

  em[10] = 7.18e-12*(kr2*fn*gn + kr3*gn*gn + kr4*gn*(n_el/N_H))*N_H*N_H; 

/* Cooling due to gas grain cooling */

  em[11] = krvalues[8]*N_H*N_H*fn*fn; //3.8e-33*st*(T-Tdust)*(1.0 - 0.8*exp(-75.0/T))*N_H*N_H;

  em[13]  = em[11] + em[10] + em[7] + (em[8] + em[9])*n_el*N_H;    
     
/* ---------------------------------------------------
    rlosst is the energy loss in units of erg/cm^3/s;
    it must be multiplied by cost_E in order to match 
    non-dimensional units.
    Source term for the neutral fraction scales with 
    UNIT_TIME
   --------------------------------------------------- */

  rlosst  =  em[13];
  rhs[RHOE] = -E_cost*rlosst;
  rhs[RHOE] *= 1.0/(1.0 + exp(-(T - g_minCoolingTemp)/100.0)); /* -- lower cutoff -- */
}


#ifdef CHOMBO
/* ********************************************************************* */
void NormalizeIons (double *u)
/*!
 *  This function re-normalize the set of conservative variables u in
 *  such a way that the sum of all ion fractions is equal to 1. 
 *  Since u is a vector of conserved quantities, the normalization 
 *  condition is:
 *
 *      \sum_j (\rho*X_j) = \rho
 * 
 *  where, for instance, X_j = {OI, OII, OIII, OIV, OV}.
 *  In order to fulfill the previous equation, each conservative
 *  variable is redefined according to 
 *
 *                               (\rho*X_j)
 *     (\rho*X_j) -->  \rho * -------------------
 *                             \sum_j (\rho*X_j) 
 *
 *  This function is necessary only when used with Chombo AMR since
 *  the coarse-to-fine interpolation fails to conserve the correct
 *  normalization. Re-fluxing and fine-to-coarse do not need this
 *  treatment since the normalization is preserved.
 *
 *********************************************************************** */
{
  int nv;
  double phi;
  
  phi = 2.0*u[X_H2] + u[X_HI] + u[X_HII];
  for (nv = X_HI; nv <= X_HII; nv++) u[nv] *= u[RHO]/phi;
}
#endif
