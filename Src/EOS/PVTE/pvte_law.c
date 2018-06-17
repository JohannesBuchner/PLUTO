/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Basic interface definition for the \c PVTE_LAW EoS.
  
  Collect the basic set of functions required by the 
  \c PVTE_LAW equation of state:

  - InternalEnergyFunc() defines the gas internal energy as a function
    temperature and ionization fractions (for non-equilibrium chemistry) 
    or temperature and density (in Local Thermodynamic Equilibrium - LTE -
    or Collisional Ionization Equilibrium - CIE).
  - GetMu() computes the mean molecular weight (in LTE or CIE).
  
  Auxiliary functions (needed by the the previous ones) are also
  included here.
  
  The current version implements the EoS given by D'Angelo (2013).

  \author A. Mignone (mignone@ph.unito.it)\n
          B. Vaidya
  \date   20 June, 2014

  \b Reference 
     - D'Angelo et. al  ApJ 778, 2013 (Eq. [26-27])

*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h" 

#define HELIUM_IONIZATION   NO

#define DEG_x      0  /* hydrogen ionization degree */
#define DEG_y      1  /* molecular hydrogen dissociation degree */
#define DEG_z1     2  /* helium single ionization degree */
#define DEG_z2     3  /* helium double ionization degree */


void        GetFuncDum(double, double *);
#if NIONS == 0
 static void GetSahaHFracs (double, double, double *);
#else
 static void GetHFracs (double *, double *);
#endif

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
  double eH, eHe, eHpH, eHplus, eH2, rhoe, eHeplus, eHeplusplus;
  double func_zetaR;
  double f[4];

  #if NIONS == 0
   GetSahaHFracs(T, v[RHO], f);
  #else
   GetHFracs (v, f);
   if (f[DEG_x] > 1.001 || f[DEG_y] > 1.001){
     printf ("! InternalEnergy: x or y are too large %f, %f\n",
              f[DEG_x],f[DEG_y]);
     QUIT_PLUTO(1);
   }
  #endif

  GetFuncDum(T, &func_zetaR);   /* = (1.5 + e(rot) + e(vib)) */
/*   func_zetaR = 1.5;   to recover ideal EoS  */ 
  
/* -- Estimate contributions to Egas -- */

  eH   = 1.5*H_MASS_FRAC*(1.0 + f[DEG_x])*f[DEG_y]; 
  eHe  = 3.0*He_MASS_FRAC/8.0;
  eH2  = 0.5*H_MASS_FRAC*(1.0 - f[DEG_y])*func_zetaR; 
 
/* -- constant terms -- */
  
  #if COOLING == NO
   eHpH   = 4.48*CONST_eV*H_MASS_FRAC*f[DEG_y]/(2.0*CONST_kB*T);
   eHplus = 13.60*CONST_eV*H_MASS_FRAC*f[DEG_x]*f[DEG_y]/(CONST_kB*T);
   
   #if HELIUM_IONIZATION == YES
    eHeplus = 24.59*CONST_eV*He_MASS_FRAC*f[DEG_z1]*(1.0 - f[DEG_z2])/(4.0*CONST_kB*T);
    eHeplusplus = 54.42*CONST_eV*He_MASS_FRAC*f[DEG_z1]*f[DEG_z2]/(4.0*CONST_kB*T);
   #else
    eHeplus = eHeplusplus = 0.0;
   #endif
  
  #else
   eHpH = eHplus = eHeplus = eHeplusplus = 0.0;
  #endif

/* ----------------------------------------------------------------
    Compute rhoe in cgs except for density which is left in 
    code units (to avoid extra multiplication and division later 
   ---------------------------------------------------------------- */
      
  rhoe = (eH2 + eH + eHe + eHpH + eHplus + eHeplus + eHeplusplus)*(CONST_kB*T*v[RHO]/CONST_mH);
  
  return rhoe/(UNIT_VELOCITY*UNIT_VELOCITY);  /* convert to code units */
}

#if NIONS == 0
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
 *  \b Reference 
 *     - Black and Bodenheimer, ApJ 199, 1975 (Eq. 8) 
 *     - D'Angelo, G. et al ApJ 778, 2013 (Eqs. 15, 16.) 
 *********************************************************************** */
{
  double f[4];

  GetSahaHFracs(T, rho, f);
  #if HELIUM_IONIZATION == YES
  *mu = 4.0/(  2.*H_MASS_FRAC*(1.0 + f[DEG_y] + 2.0*f[DEG_x]*f[DEG_y]) 
               + He_MASS_FRAC*(1.0 + f[DEG_z1]*(1.0 + f[DEG_z2])));
  #else
  *mu = 4.0/(2.*H_MASS_FRAC*(1.0 + f[DEG_y] + 2.0*f[DEG_x]*f[DEG_y]) 
             + He_MASS_FRAC);
  #endif
}

/* ********************************************************************* */
void GetSahaHFracs(double T, double rho, double *fdeg)
/*!
 * Compute degree of ionization and dissociation using Saha Equations. 
 * The quadratic equation ay^2 + by + c = 0  is solved using 
 * a = 1, c = -b (for hydrogen). A similar expression is used for Helium.
 *
 * \param [in]   T     Gas temperature in Kelvin.
 * \param [in]   rho   Gas density (in code units)
 * \param [out]  fdeg  array of ionization/dissociation degrees.
 *
 * \b  Reference 
 *     - D'Angelo, G. et al ApJ 778, 2013 (Eqs. 15, 16.) 
 *
 *********************************************************************** */
{
  double rhs1, rhs2, rhs3;
  double b,c, scrh;
  double kT   = CONST_kB*T;
  double hbar = CONST_h/(2.0*CONST_PI);

  rho *= UNIT_DENSITY;

  rhs1 = CONST_mH/(2.0*H_MASS_FRAC*rho);
  rhs2 = CONST_mH*kT/(4.0*CONST_PI*hbar*hbar);
  rhs3 = exp(-4.48*CONST_eV/kT);
  b    = rhs1*rhs2*sqrt(rhs2)*rhs3;

  fdeg[DEG_y] = 2.0/(1.0 + sqrt(1.0 + 4.0/b)); /* solution of quadratic equation */

  rhs1 = CONST_mH/(H_MASS_FRAC*rho);
  rhs2 = CONST_me*kT/(2.0*CONST_PI*hbar*hbar);
  rhs3 = exp(-13.60*CONST_eV/kT);
  b    = rhs1*rhs2*sqrt(rhs2)*rhs3;

  fdeg[DEG_x] = 2.0/(1.0 + sqrt(1.0 + 4.0/b)); /* solution of quadratic equation */
  
  #if HELIUM_IONIZATION == YES
   rhs3  = 4.0*CONST_mH/rho*rhs2*sqrt(rhs2)*exp(-24.59*CONST_eV/kT);
   b     = 4.0/He_MASS_FRAC*(H_MASS_FRAC + rhs3);
   c     = -rhs3*4.0/He_MASS_FRAC;
   fdeg[DEG_z1] = -2.0*c/(b + sqrt(b*b - 4.0*c));
   
   rhs3  = CONST_mH/rho*rhs2*sqrt(rhs2)*exp(-54.42*CONST_eV/kT);
   b     = 4.0/He_MASS_FRAC*(H_MASS_FRAC + 0.25*He_MASS_FRAC + rhs3);
   c     = -rhs3*4.0/He_MASS_FRAC;
   fdeg[DEG_z2] = -2.0*c/(b + sqrt(b*b - 4.0*c));
  #endif
}

#else  

/* ********************************************************************* */
void GetHFracs(double *v, double *f)
/*!
 * Compute degree of ionization and dissociation 
 * from non-equillibrium evolution of various hydrogen fractions 
 * using the cooling module. 
 *
 * \param [in]   v    1D Array of primitive variables.
 * \param [out]  x    Deree of Ionization of hydrogen.
 * \param [out]  y    Degree of Dissociation of molecular hydrogen.
 *
 *********************************************************************** */
{
  #if COOLING == SNEq
   f[DEG_x] = 1.0 - v[X_HI];
   f[DEG_y] = 1.0;
  #elif COOLING == H2_COOL
   /* Ionization Fraction. */
   f[DEG_x] = v[X_HII]/(v[X_HII] + v[X_HI] + 1.e-40); 

 /* ----------------------------------------------------------
     Dissociation Fraction y = 1.0 - mol. fraction,
     the small epsilon = 1.0e-12 in mol.fraction is used
     to avoid the NaN in case of very high temperature limit
     where X_H2 = X_HI = 0.0 and fHII = 1.0.
     Further it satisfies :  y -> 1.0 as T -> Infinity. 
    ---------------------------------------------------------- */
      
   f[DEG_y] = 1.0 - 2.0*v[X_H2]/(2.0*v[X_H2] + v[X_HI] + 1.0e-40);

  #else
   print ("! GetHFracs: not defined for this cooling\n");
   QUIT_PLUTO(1);
  #endif
}
#endif  /* NIONS == 0 */

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
 *  \b Reference 
 *     - D'Angelo et. al  ApJ 778, 2013 (Eq. [26-27])
 *
 * \param [in]    v       1D array of primitive quantities
 * 
 * \return Value of first adiabatic index Gamma1.
 *********************************************************************** */
{
  double gmm1, T;
  double cv, mu, chirho = 1.0, chiT = 1.0, rho, rhoe;
  double ep, em, delta = 1.e-2;
  double Tp, Tm, mup, mum, rhop, rhom;
  double dmu_dT, dmu_drho, de_dT;

  return 1.6667; /* !! For the current implementation, we always return
                       5/3 rather then computing the actual Gamma1 index
                       (which is quite expensive).
                        We do not think this is a major problem since
                       Gamma1 only appears in the estimate of the sound speed
                       needed by the Riemann solver. An upper bound will do. */
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
  em = InternalEnergy(v, Tm)/v[RHO]; /* in code units */
  ep = InternalEnergy(v, Tp)/v[RHO]; /* in code units */
  
  de_dT = (ep - em)/(2.0*delta*T);
  cv    = de_dT;  /* this is code units. */

  #if NIONS == 0
   rho = v[RHO];

   GetMu(Tp, rho, &mup);
   GetMu(Tm, rho, &mum);
   dmu_dT = (mup - mum)/(2.0*delta*T);      
   chiT  = 1.0 - T/mu*dmu_dT;

   rhop = rho*(1.0 + delta);
   rhom = rho*(1.0 - delta);
   GetMu(T, rhop, &mup);
   GetMu(T, rhom, &mum);
   dmu_drho = (mup - mum)/(2.0*delta*rho);      
   chirho   = 1.0 - rho/mu*dmu_drho;
  #endif

/* --------------------------------------------
    Compute first adiabatic index
   -------------------------------------------- */
   
  gmm1   = v[PRS]/(cv*T*v[RHO])*chiT*chiT  + chirho;

  return gmm1;
}

#if NIONS > 0
/* ********************************************************************* */
void InternalEnergyBracket (double rhoe, double *v, double *Tlo, double *Thi)
/*
 *
 *********************************************************************** */
{
  double rho = v[RHO];
  const double norm=CONST_kB*rho/CONST_mH;
  const double a = 1.5*H_MASS_FRAC;
  const double b = 3.0*He_MASS_FRAC/8.0;
  const double c = 0.5*H_MASS_FRAC;
  double q = rhoe*(UNIT_VELOCITY*UNIT_VELOCITY)/norm;
  double f[4];

  GetHFracs (v, f);
   
  (*Thi) = q/(a*(1.0 + f[DEG_x])*f[DEG_y] + b + c*(1.0 - f[DEG_y])*1.49);
  (*Tlo) = q/(a*(1.0 + f[DEG_x])*f[DEG_y] + b + c*(1.0 - f[DEG_y])*3.7);

/* -- enlarge the interval by a little to avoid problems when y = 1 -- */

  (*Tlo) *= 0.99;  /* avoid negative temperature */
  (*Thi) += 1.0;
}
#endif

