/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Templeate file for the \c PVTE_LAW EoS.
  
  Collect the basic set of functions required by the 
  \c PVTE_LAW equation of state:

  - InternalEnergyFunc() defines the gas internal energy as a function
    temperature and ionization fractions (for non-equilibrium chemistry) 
    or temperature and density (in Local Thermodynamic Equilibrium - LTE -
    or Collisional Ionization Equilibrium - CIE).
  - GetMu() computes the mean molecular weight (in LTE or CIE).
  
  \author A. Mignone (mignone@ph.unito.it)\n
          B. Vaidya
  \date   20 June, 2014

*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h" 

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
 *  \b Reference 
 *     - D'Angelo et. al  ApJ 778, 2013 (Eq. [26-27])
 *
 * \param [in]    v       1D array of primitive quantities
 * 
 * \return Value of first adiabatic index Gamma1.
 *********************************************************************** */
{
  return 1.6667; 
}
