/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute equilibrium fractions for the H2_COOL module.

  Compute the equilibrium fractions of atomic and molecular hydrogen 
  for a given density and temperature (Collisional Excitation 
  Equlibrium or CIE).
  This corresponds to the simultaneous solution of the three
  rate equations

      R_HI (f,h,g) = 0 
      R_HII(f,h,g) = 0 
      R_H2 (f,h,g) = 0 

  where f = X(HI), h = X(HII), g = X(H2) are the number fraction of
  atomic, ionized and molecular hydrogen, respectively.
  The functions R_HI, R_HII and R_H2 are given in radiat.c.
  The previous system is reduced to a single scalar equation in either
  f or h by using  h = ci/cr*f and the normalization condition
  g = (1-h-f)/2. 
  A maple script can be found at the end of this file.
  The final equation is a quadratic,
  
      a*f^2 + b*f + c = 0
         
  and we have verified that the solution with the - sign is the physical 
  acceptable one since 0 <= f <= 1 for any T.
  
  Note that the solution for g at very large temperatures may suffer
  from machine precision when g becomes very small (g < 1.e-19). 
  Using long double for the coefficients of the quadratic helps 
  a little bit.  
 
  \authors A. Mignone (mignone@ph.unito.it)\n
           B. Vaidya 

  \date    March 1, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void CompEquil(double n, double T, double *v)
/*!
 *  \param [in]      n  the particle number density (not needed, but
 *                      kept for compatibility)
 *  \param [in]      T  the temperature (in K) for which equilibrium must 
 *                      be found.
 *  \param [in,out] v  an array of primitive variables. 
 *                      On input, only density needs to be defined.
 *                      On output, fractions will be updated to the 
 *                      equilibrium values.
 *
 *********************************************************************** */
{
  double cr, ci, kr1, kr2, kr3, kr4;
  double frac_Z, st, tev;
  long double a, b, c, qc, scrh, az;
  long double f, h, g;
  double krv[9];

  H2RateTables(T, krv);
  kr1 = krv[0];
  kr2 = krv[1];
  kr3 = krv[2];
  kr4 = krv[3];
  cr  = krv[4];
  ci  = krv[5];

  kr2 /= kr1;
  kr3 /= kr1;
  kr4 /= kr1;
  kr1  = 1.0;
  frac_Z = ((1.0 - H_MASS_FRAC - He_MASS_FRAC)/CONST_AZ)*(CONST_AH/H_MASS_FRAC);

  az   = 0.5*CONST_AZ*frac_Z;

/* ---------------------------------------------------
    Solve quadratic equation for fn or hn.
    The solution depends only on temperature and the 
    physical solution is the one with - sign.
   --------------------------------------------------- */
  
  if (T < 1.e16) {   /* Solve for fn */
    qc   = ci/cr;
    scrh = -0.5*(1.0 + qc);
    a    = (long double)(scrh*(kr2 + kr3*scrh + kr4*qc) - kr1);
    b    = (long double)( 0.5*(kr2 + kr3*scrh + kr4*qc + scrh*(kr3 + 2.0*kr4*az)));
    c    = (long double)(0.5*(0.5*kr3 + kr4*az));

  /* ------------------------------------------
      Set fn = 0.0 below T = 300 K to prevent
      machine accuracy problems
     ------------------------------------------ */

    if (T > 300.0){
      if   (b >= 0.0) f = -(b + sqrtl(b*b - 4.0*a*c))/(2.0*a);
      else            f = 2.0*c/(sqrtl(b*b - 4.0*a*c) - b);
    }else{
      f = 0.0;
    }  
    h = qc*f;
    g = 0.5*(1.0 - f*(1.0 + qc));

  /* -- Set g = 0 if we are close to machine accuracy -- */

/*    if ( fabs(1.0 - h - f) < 1.e-15) g = 0.0;  */
 
  }else{             /* Solve for hn */

    qc   = cr/ci;
    scrh = -0.5*(1.0 + qc);
    a    = scrh*(kr2*qc + kr3*scrh + kr4) - kr1*qc*qc;
    b    =  0.5*(kr2*qc + kr3*scrh + kr4 + scrh*(kr3 + 2.0*kr4*az));
    c    =  0.5*(0.5*kr3 + kr4*az);
    h    = -(b + sqrtl(b*b - 4.0*a*c))/(2.0*a);
    f    = qc*h;
    g    = 0.5*(1.0 - h*(1.0 + qc));
  }

  v[X_HI]  = f;
  v[X_HII] = h;
  v[X_H2]  = g;
}

/*
########################################################################
# CIE solution for the H2COOL cooling module.
# Reduce the three rate equations to one scalar equation by using
# the third one (h = ci/cr)*f and the normalization condition.
# Solve for f
########################################################################
restart;
h   := ci*f/cr;
x   := h + az;
g   := (1 - h - f)/2;

RI  := (g*(k2*f + k3*g + k4*x) + cr*h*x - f*(ci*x + k1*f));

coeff(RI,f,2);
coeff(RI,f,1);
coeff(RI,f,0);

########################################################################
# Same as before but solve for h
########################################################################
restart;
f   := cr*h/ci;
x   := h + az;
g   := (1 - h - f)/2;

RI  := (g*(k2*f + k3*g + k4*x) + cr*h*x - f*(ci*x + k1*f));

coeff(RI,h,2);
coeff(RI,h,1);
coeff(RI,h,0);

########################################################################
# Same as before but solve for g
########################################################################
restart;
f   := (1 - 2*g)/alpha; #(1 + ci/cr);
h   := ci*f/cr;
x   := h + az;

RI  := (g*(k2*f + k3*g + k4*x) + cr*h*x - f*(ci*x + k1*f));

coeff(RI,g,2);
coeff(RI,g,1);
coeff(RI,g,0);
*/



