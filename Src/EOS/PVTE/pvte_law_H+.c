/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief PVTE_LAW for a partially ionized gas.

  Compute the internal energy for a partially ionized hydrogen gas.

  \author A. Mignone (mignone@ph.unito.it)\n
          B. Vaidya
  \date   29 Aug 2014

  \b Reference 
     - PLUTO Users' guide.
*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h" 

static double SahaXFrac(double T, double rho);
#define DIFF_ORDER  4   /* = 2 or 4 for 2nd or 4th accurate approximations
                           to derivatives in Gamma1()  */
#define INTE_EXACT 1    /* Compute internal energy exactly in Gamma1() */

/* ***************************************************************** */
double InternalEnergyFunc(double *v, double T)
/*!
 *  Compute the gas internal energy as a function of temperature
 *  and fractions (or density):
 *  - <tt> rhoe = rhoe(T,rho) </tt> in LTE or CIE;
 *  - <tt> rhoe = rhoe(T,X) </tt> in non-equilibrium chemistry.
 *
 *  \param [in]   v   1D Array of primitive variables containing
 *                    density and species. Other variables are
 *                    ignored.
 *  \param [in]   T   Gas temperature
 *
 *  \return The gas internal energy (\c rhoe) in code units.
 ******************************************************************* */
{
  double x, rho, e, mu, kT;
  double rhoe, chi = 13.6*CONST_eV;
  double p0 = UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY;

  kT   = CONST_kB*T;
  x    = SahaXFrac(T, v[RHO]);
  GetMu(T,v[RHO],&mu); 
  rho  = v[RHO]*UNIT_DENSITY;
  e    = 1.5*kT/(mu*CONST_amu) + chi*x/CONST_mH;

/*
FILE *fp;
double p, lnT;
fp = fopen("gamma.dat","w");
v[RHO] = 1.0;
for (lnT = 2; lnT < 6; lnT += 0.01){
  T = pow(10.0, lnT);
  kT   = CONST_kB*T;
  x    = SahaXFrac(T, v[RHO]);
  GetMu(T,v[RHO],&mu); 
  rho  = v[RHO]*UNIT_DENSITY;
  e    = 1.5*kT/(mu*CONST_amu) + chi*x/CONST_mH;
p = kT/(mu*CONST_amu)*rho;
fprintf (fp,"%12.6e   %12.6e   %12.6e\n",T,p/(rho*e) + 1.0, x);
}
fclose(fp);
exit(1);
*/

  return rho*e/p0; /* Convert to code units */ 
}

/* ********************************************************************* */
double SahaXFrac(double T, double rho)
/*!
 * Use Saha equation to compute the degree of ionization.
 *
 *
 *********************************************************************** */
{
  double me, kT, h3, c, x, n;
  double chi = 13.6*CONST_eV;

  rho *= UNIT_DENSITY; 

  me = 2.0*CONST_PI*CONST_me;
  kT = CONST_kB*T;
  h3 = CONST_h*CONST_h*CONST_h;

  n  = rho/CONST_mp; /* = n(protons) + n(neutrals)   not   n(total) */
  c  = me*kT*sqrt(me*kT)/(h3*n)*exp(-chi/kT);
  
  x = 2.0/(sqrt(1.0 + 4.0/c) + 1.0);
  return x;
}

/* ********************************************************************* */
void GetMu(double T, double rho, double *mu)
/*!
 *  Calculate the mean molecular weight for the case in which 
 *  hydrogen fractions are estimated using Saha Equations.
 *
 *  \param [in]   T    Gas temperature in Kelvin.
 *  \param [in]   rho  Gas density (code units)
 *  \param [out]  mu   Mean molecular weight
 *
 *********************************************************************** */
{
  double x;
  x = SahaXFrac(T, rho);
  *mu = 1.0/(1.0 + x);
}

/* ********************************************************************* */
double Gamma1(double *v)
/*!
 *  Calculate the value of the first adiabatic index:
 *  \f[
 *     \Gamma_1 = \frac{1}{c_V}\frac{p}{\rho T} \chi_T^2 + \chi_\rho^2
 *     \qquad{\rm where}\quad
 *     \chi_T = \left(\pd{\log p}{\log T}\right)_{\rho} = 
 *              1 - \pd{\log{\mu}}{\log T} \,;\quad
 *     \chi_\rho = \left(\pd{\log p}{\log\rho}\right)_{T} = 
 *                  1 - \pd{\log{\mu}}{\log\rho} \,;\quad
 *  \f]
 *  where \c p and \c rho are in c.g.s units.
 *  Note that if species are evolved explicitly (non-equilibrium chemistry),
 *  we set \c chi=1.
 *
 *  The heat capacity at constant volume, \c cV, is defined as the 
 *  derivative of specific internal energy with respect to temperature:
 *  \f[
 *      c_V = \left.\pd{e}{T}\right|_V
 *  \f]
 *  and it is computed numerically using a centered derivative.
 *
 *  This function is needed (at present) only when computing the 
 *  sound speed in the Riemann solver.
 *  Since this is only needed for an approximated value, 5/3 (upper
 *  bound) should be ok.
 *
 * \param [in]    v       1D array of primitive quantities
 * 
 * \return Value of first adiabatic index Gamma1.
 *********************************************************************** */
{
  double gmm1, T;
  double cv, mu, chirho = 1.0, chiT = 1.0, rho, rhoe;
  double epp, emm, ep, em, delta = 1.e-3;
  double Tpp, Tmm, Tp, Tm, mupp, mumm, mup, mum, rhopp, rhomm, rhop, rhom;
  double dmu_dT, dmu_drho, de_dT;

  return 1.666667; 
  
/* ---------------------------------------------
    Obtain temperature and fractions.
   --------------------------------------------- */
   
  #if NIONS == 0
   GetPV_Temperature(v, &T);
   GetMu(T, v[RHO], &mu);
  #else
   mu  = MeanMolecularWeight(v); 
   T   = v[PRS]/v[RHO]*KELVIN*mu;
  #endif

/* ---------------------------------------------------
    Compute cV (Specific heat capacity for constant 
    volume) using centered derivative. 
    cV will be in code units. The corresponding cgs 
    units will be the same velocity^2/Kelvin.
   --------------------------------------------------- */
  
  Tp = T*(1.0 + delta);
  Tm = T*(1.0 - delta);
  Tmm = T*(1.0 - 2.0*delta);
  Tpp = T*(1.0 + 2.0*delta);

  #if INTE_EXACT
   emm = InternalEnergyFunc(v, Tmm)/v[RHO]; /* in code units */
   em  = InternalEnergyFunc(v, Tm)/v[RHO]; /* in code units */
   ep  = InternalEnergyFunc(v, Tp)/v[RHO]; /* in code units */
   epp = InternalEnergyFunc(v, Tpp)/v[RHO]; /* in code units */
  #else
   emm = InternalEnergy(v, Tmm)/v[RHO]; /* in code units */
   em  = InternalEnergy(v, Tm)/v[RHO]; /* in code units */
   ep  = InternalEnergy(v, Tp)/v[RHO]; /* in code units */
   epp = InternalEnergy(v, Tpp)/v[RHO]; /* in code units */
  #endif
 
  #if DIFF_ORDER == 2
   de_dT = (ep - em)/(2.0*delta*T);
  #elif DIFF_ORDER == 4  
   de_dT = (emm - 8.0*em + 8.0*ep - epp)/(12.0*delta*T);
  #else
   #error Order not allowed in Gamma1()  
  #endif
 
  cv    = de_dT;  /* this is code units. */

  #if NIONS == 0
   rho = v[RHO];

   GetMu(Tp, rho, &mup);
   GetMu(Tm, rho, &mum);
   GetMu(Tpp, rho, &mupp);
   GetMu(Tmm, rho, &mumm);

   #if DIFF_ORDER == 2
    dmu_dT = (mup - mum)/(2.0*delta*T);      
   #elif DIFF_ORDER == 4
    dmu_dT = (mumm - 8.0*mum + 8.0*mup - mupp)/(12.0*delta*T);
   #else
    #error Order not allowed in Gamma1()     
   #endif   

   chiT  = 1.0 - T/mu*dmu_dT;

   rhop = rho*(1.0 + delta);
   rhom = rho*(1.0 - delta);
   rhopp = rho*(1.0 + 2.0*delta);
   rhomm = rho*(1.0 - 2.0*delta);

   GetMu(T, rhop, &mup);
   GetMu(T, rhom, &mum);
   GetMu(T, rhopp, &mupp);
   GetMu(T, rhomm, &mumm);

   #if DIFF_ORDER == 2
    dmu_drho = (mup - mum)/(2.0*delta*rho);      
   #elif DIFF_ORDER == 4
    dmu_drho = (mumm - 8.0*mum + 8.0*mup - mupp)/(12.0*delta*rho); 
   #else
    #error Order not allowed in Gamma1()     
   #endif   
   chirho   = 1.0 - rho/mu*dmu_drho;
  #endif

/* --------------------------------------------
    Compute first adiabatic index
   -------------------------------------------- */
/*  
printf ("gmm1:  prs = %12.6e, T = %12.6e, chiT = %12.6e, chirho = %12.6e\n",
        v[PRS], T, chiT, chirho);
*/
  gmm1   = v[PRS]/(cv*T*v[RHO])*chiT*chiT  + chirho;

//printf ("gmm1 = %f\n",gmm1);
  return gmm1;
  
}
