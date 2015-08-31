/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute the mean molecular weight.

  Compute and return the mean molecular weight as a function of the gas
  composition under non-equilibrium conditions.
  The mean molecular weight is usually needed to compute the temperature
  or mass density from the particle number density:
  \f[
    T = \frac{p}{n_{\rm tot}k_B} = \frac{p}{\rho}\frac{m_u\mu}{k_B}
    \,,\qquad
    \rho = \mu m_u n_{\rm tot}
  \f]
  where \f$ m_u \f$ is the atomic mass unit while 
  \f$ n_{\rm tot} \f$ is the number density of all particles. 
  The Mean molecular weight is defined as the average mass of a particle of 
  gas in terms of the atomic mass unit and it is expressed by the weighted
  sum of the mass of particles in atomic mass unit divided by total number
  of particles (see book by Ryan & Norton [Eq. 1.7])
  \f[ 
    \mu = \frac{\DS\sum_k N_k \frac{m_k}{m_u}}{\DS\sum_k N_k}
  \f]  
  where 
  - \f$ N_k \f$ is the number of particles in the gas of element k and it
    can be related to mass fraction \f$X_k\f$ as 
   \f$\DS \frac{N_k}{V} = n_k = \frac{\rho}{m_u}\frac{X_k}{A_k}\f$, where
   \f$A_k\f$ is the atomic mass number.
  - \f$ m_k \f$ is the mass of each particle of element k.
   
  \b H2_COOL: to compute \f$\mu\f$ we proceed as follows:
  -# using standard solar composition \f$ X_{\odot},\, Y_{\odot},\, Z_{\odot}\f$
     we derive 
     \f[
       \frac{N_{He}}{N_H} = \frac{Y_{\odot}}{A_{He}}\frac{A_H}{X_\odot}
       ,\,\qquad
       \frac{N_Z}{N_H} = \frac{Z_{\odot}}{A_{Z}}\frac{A_H}{X_\odot}
     \f]
  -# The weighted sum of the mass of particles in the numerator of \f$\mu\f$
     is given by 
     \f[
       \sum_k N_k \frac{m_k}{m_u} =   N_H   \frac{m_H}{m_u}
                                    + N_{He}\frac{m_{He}}{m_u}
                                    + N_Z   \frac{m_Z}{m_u}
     \f]
     Note that since mass of electron is negligible, so is its contribution to
     the previous summation.
  -# For the denominator we have
     \f[
       \sum_k N_k = N_{HI} + N_{HII} + N_{H2} + N_e + N_{He} + N_Z
                   + \frac{A_ZN_Z}{2}
     \f]
     where two sources of electrons considered here: \f$ N_e \f$ electrons
     corresponding to \f$ N_{HII} \f$ protons due to ionization of hydrogen and
     \f$ A_ZN_Z/2 \f$ number of electrons due to metals. 
     Note that now the electrons contribute to the total number of particles
     and cannot be  neglected.
  -# Next define the total number of hydrogen,
     \f$ N_H = N_{HI} + N_{HII} + 2N_{H2} \f$ as the sum of number of atomic
     hydrogen (HI), ionized hydrogen (HII, or protons) and twice the number 
     of molecular hydrogen (H2) and the corresponding number fractions:
     \f[
        f_{HI}  = \frac{N_{HI}}{N_H},\quad
        f_{HII} = \frac{N_{HII}}{N_H},\quad
        f_{H2}  = \frac{N_{H2}}{N_H},\quad
     \f]
     so that \f$ f_{HI} + 2f_{H2} + f_{HII} = 1 \f$.

  Putting it all together:
  \f[
     \mu = \frac{A_H + A_{He}f_{He} + A_Zf_Z}
                {f_{HI} + f_{H2} + 2f_{HII} + f_{He} + f_Z + A_Z f_Z/2}
  \f]
  where   
  - \f$A_{He}, A_Z \f$ are atomic mass numbers of helium and metals respectively.
  - \f$f_{He} = N_{He}/N_H\f$ is the fixed number fraction of helium with
    respect to hydrogen;
  - \f$f_Z = N_Z/N_H\f$ is the fixed number fraction of metals with respect to
    hydrogen.
    

  \b MINEq: please see Eq. [12] of Tesileanu (2008)

  \b SNEq: the derivation is similar to H2_COOL with \f$ f_{H2} = 0\f$ yielding 
     \f[
        \mu = \frac{A_H + A_{He}f_{He} + A_Zf_Z}
                   {2 - f_{HI} + f_{He} + 2f_Z}
     \f]
     where one electron from metals is assumed. 
  
  \b No \b Chemistry: in case where chemical reaction are not incuded,
     the mean molecular weight is computed from the mass fractions assuming
     a fully ionized gas:
     \f[
        \mu = \frac{A_H + A_{He}f_{He} + A_Zf_Z}{2  + f_{He} + f_Z(1 + A_Z/2)}
     \f]


  \b References
     - "Stellar Evolution and Nucleosynthesis"
        Sean G. Ryan and Andrew J. Norton.
        The University of Chicago Press 
     - "Simulating radiative astrophysical flows with the PLUTO code:
        a non-equilibrium, multi-species cooling function",
        Tesileanu, Mignone \& Massaglia, A\&A (2008) 488, 429

  \authors A. Mignone (mignone@ph.unito.it)\n
           B. Vaidya  

  \date   Aug 11, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
//#if (COOLING == SNEq) || (COOLING == MINEq) || (COOLING == H2_COOL)
#if (COOLING == MINEq)
  #include "cooling_defs.h"
#endif

/* ********************************************************************* */
double MeanMolecularWeight(double *v)
/*!
 *
 * Return the mean molecular weight.
 * 
 * \param [in]  v  array of primitive variables (including ions)
 *
 *********************************************************************** */
{
  int nv;
  double mu;
  
#if COOLING == NO
  
  mu =  (CONST_AH + FRAC_He*CONST_AHe + FRAC_Z*CONST_AZ) /
        (2.0 + FRAC_He + FRAC_Z*(1.0 + CONST_AZ*0.5));
    
#elif COOLING == TABULATED

  mu =  (CONST_AH + FRAC_He*CONST_AHe + FRAC_Z*CONST_AZ) /
        (2.0 + FRAC_He + FRAC_Z*(1.0 + CONST_AZ*0.5));

#elif COOLING == SNEq

  mu = (CONST_AH + FRAC_He*CONST_AHe + FRAC_Z*CONST_AZ) /
       (2.0 + FRAC_He + 2.0*FRAC_Z - v[X_HI]);
    
/*
  return  ( (CONST_AH + frac_He*CONST_AHe + frac_Z*CONST_AZ) /
            (2.0 + frac_He + 2.0*frac_Z - v[X_HI]));
*/            

#elif COOLING == H2_COOL

  double munum, muden;
 
  NIONS_LOOP(nv){
    v[nv] = MAX(v[nv], 0.0);
    v[nv] = MIN(v[nv], 1.0);
  }

  v[X_H2] = MIN(v[X_H2], 0.5);

  double fn = v[X_HI];
  double gn = v[X_H2];
  double hn = v[X_HII];

  mu  = (CONST_AH + CONST_AHe*FRAC_He + CONST_AZ*FRAC_Z) /
        (fn + gn + 2*hn + FRAC_He + FRAC_Z + 0.5*CONST_AZ*FRAC_Z);
  
/*
  double N_H  = (H_MASS_FRAC/CONST_AH);  
  double N_He = (He_MASS_FRAC/CONST_AHe); 
  double N_Z  = ((1.0 - H_MASS_FRAC - He_MASS_FRAC)/CONST_AZ);  

  double fracHe = N_He/N_H;
  double fracZ  = N_Z/N_H;
  
  double fn = v[X_HI];
  double gn = v[X_H2];
  double hn = v[X_HII];
 
  munum = 1.0 + CONST_AHe*(fracHe) + CONST_AZ*(fracZ);
  muden = fn + gn + 2*hn + fracHe + fracZ + 0.5*CONST_AZ*(fracZ);

  return munum/muden;
*/
#elif COOLING == MINEq

  double mmw1, mmw2;
  int    i, j;
  
  mmw1 = mmw2 = 0.0;
  for (i = 0; i < NIONS; i++) {
    if (v[NFLX+i] < 0.0) v[NFLX+i] = 0.0;
    if (v[NFLX+i] > 1.0) v[NFLX+i] = 1.0;
    CoolCoeffs.dmuN_dX[i] = elem_mass[elem_part[i]]*elem_ab[elem_part[i]];
    CoolCoeffs.dmuD_dX[i] = elem_ab[elem_part[i]]  *rad_rec_z[i];
    mmw1 += CoolCoeffs.dmuN_dX[i]*v[NFLX+i];  /*    Numerator part of mu    */
    mmw2 += CoolCoeffs.dmuD_dX[i]*v[NFLX+i];  /*    Denominator part of mu  */    
  }

/* -- Add contributions from ionized H --  */

  CoolCoeffs.dmuN_dX[0] += -elem_mass[0]*elem_ab[el_H];
  CoolCoeffs.dmuD_dX[0] += -2.0*elem_ab[el_H];

  mmw1 += elem_mass[0]*elem_ab[el_H]*(1.0 - v[X_HI]); 
  mmw2 += elem_ab[el_H]*(1.0 - v[X_HI])*2.;

  CoolCoeffs.muN = mmw1;
  CoolCoeffs.muD = mmw2;

  if (mmw1 != mmw1) {
    print(">>> Error!  MMW1  NaN! %ld\n",g_stepNumber);
    for (i = 0; i < NIONS; i++) {
       print ("%d   %10.4e\n",i,v[NFLX+i]);
    }
    QUIT_PLUTO(1);
  }
  if (mmw2 != mmw2) {
    print(">>> Error!  MMW2  NaN!\n");
    for (i = 0; i < NIONS; i++) {
       print ("%d   %10.4e\n",i,v[NFLX+i]);
    }
    QUIT_PLUTO(1);
  }
  
  mu = mmw1/mmw2;

#endif
    
  return mu;
}