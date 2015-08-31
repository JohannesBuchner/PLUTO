/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Propagation of a conduction front.

  The problem sets the initial condition for a propagating conduction front
  for which an analytical solution exists (Reale 1995).

  The equation
  \f[
    \frac{\partial T}{\partial t} = a\frac{\partial}{\partial s}
    \left(T^n\frac{\partial T}{\partial s}\right)
  \f]
  has the solution
  \f[
    T = T_c\left(1 - \frac{s^2}{s_f^2}\right)^{1/n}
  \f]
  where
  \f[
    T_c = \left(\frac{Q^2}{at}\right)^{1/(n+2)}
      \left(\frac{n}{2(n+2)}\zeta_0^2\right)^{1/n}\,,
    \quad\quad\quad
    s_f = (aQ^nt)^{1/(n+2)}\zeta_0\,,
    \quad\quad\quad
    \zeta_0 = \left[\frac{(n+2)^{1+n}2^{1-n}}{n\pi^{n/2}}\right]^{1/(n+2)}
      \left[\frac{\Gamma(1/2 + 1/n)}{\Gamma(1/n)}\right]^{n/(n+2)}\,,
  \f]
  \f$Q\f$ is the integral over the whole space, and \f$\Gamma\f$ is the
  gamma function.

  The setup is built to compare the numerical solution with the analytical one.

  \image html 1D-hd-exp.jpg "Evolution of the TC front in the 1D case; numerical solution (points) and analytical (lines)".

  Configurations (\c EXPLICIT and \c STS):
  - #01-#04: 1D cartesian
  - #05-#06: 2D cartesian
  - #07-#09: 3D cartesian
  - #10-#11: 3D polar
  - #12-#13: 3D spherical

  \author A. Mignone (mignone@ph.unito.it)
  \date   July 09, 2014 

  \b References:
    - Section 4 of Reale, F. Computer Physics Communications (1995), 86, 13

*/
/* ///////////////////////////////////////////////////////////////////// */


#include "pluto.h"
/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *********************************************************************** */
{
#if DIMENSIONS == 1
  double n0, k0, Q, t0, j0, xf, Tc, T0;
  double G1_2p1_n, G1_n;

  n0 = 5./2.;
  k0 = 4.412;
  Q = 1.2e15;
  t0 = 0.1;

  G1_2p1_n = 1.068628702119319;
  G1_n = 2.218159543757688;

  j0 = pow((pow(n0+2.,1.+n0)*pow(2.,1.-n0))/(n0*pow(CONST_PI,n0/2.)),1./(n0+2.))*pow(G1_2p1_n/G1_n,n0/(n0+2.));
  xf = pow(k0*pow(Q,n0)*t0, 1./(n0+2.))*j0/1.0e8;
  Tc = pow(Q*Q/(k0*t0), 1./(n0+2.))*pow(n0*j0*j0/(2*(n0+2.)), 1./n0);
  T0 = 6.05736872274111e7;

  us[RHO] = 1.0;
  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = 0.0;

  if ((1.0 - x1*x1/xf/xf) < 1.0e-20){
     us[PRS] = us[RHO]*Tc/T0*pow(1.0e-20,0.4);
  }else{
     us[PRS] = us[RHO]*Tc/T0*pow( (1.0 - x1*x1/xf/xf) ,0.4);
  }
  us[TRC] = 0.0;
  
  #if PHYSICS == MHD
   us[BX1] = 1.0;
   us[BX2] = 0.0;
   us[BX3] = 0.0;
  #endif

#else

  double n0, k0, Q, t0, j0, xf, Tc, T0;
  double G5_2p1_n, G1p1_n, G3_2, x;

  n0 = 5./2.;
  k0 = 4.412;
  Q = 1.2e30;
  t0 = 0.1;

  G5_2p1_n = 1.827355080624036;
  G1p1_n = 0.887263817503075;
  G3_2 = 0.886226925452758;

  j0 = pow((3.*n0+2.)/(pow(2., n0-1.)*n0*pow(CONST_PI,n0)),1./(3.*n0+2.))*pow(G5_2p1_n/G1p1_n/G3_2,n0/(3.*n0+2.));
  xf = pow(k0*pow(Q,n0)*t0, 1./(3.*n0+2.))*j0/1.e8;
  Tc = j0*j0*j0/xf/xf/xf*1.e-24*Q*pow(n0*j0*j0/(2.*(3.*n0+2.)), 1./n0);
  T0 = 6.05736872274111e7;

  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
   x = sqrt(D_EXPAND(x1*x1, + x2*x2, + x3*x3));
  #elif GEOMETRY == POLAR
   x = sqrt(D_EXPAND(x1*x1, + 0.0, + x3*x3));
  #elif GEOMETRY == SPHERICAL
   x = x1;
  #endif


  us[RHO] = 1.0;
  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = 0.0;
  if ((1.0 - x*x/xf/xf) < 1.0e-20){
    #if GEOMETRY == SPHERICAL
     us[PRS] = us[RHO]*Tc/T0*pow(1.0e-8,0.4);
    #else
     us[PRS] = us[RHO]*Tc/T0*pow(1.0e-20,0.4);
    #endif
  }else{
     us[PRS] = us[RHO]*Tc/T0*pow( (1.0 - x*x/xf/xf) ,0.4);
  }
  us[TRC] = 0.0;
#endif
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{
}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  x1, x2, x3;

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){
      EXPAND(d->Vc[VX1][k][j][i] = 0.0;  ,
             d->Vc[VX2][k][j][i] = 0.0;  ,
             d->Vc[VX3][k][j][i] = 0.0;)
    }
    #if ENTROPY_SWITCH == YES
     TOT_LOOP(k,j,i) d->flag[k][j][i] |= FLAG_ENTROPY;
    #endif
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        x1 = grid[IDIR].x[i];
        d->Vc[RHO][k][j][i] =  d->Vc[RHO][k][j][2*IBEG - i - 1];
        EXPAND(
          d->Vc[VX1][k][j][i] = -d->Vc[VX1][k][j][2*IBEG - i - 1];  ,
          d->Vc[VX2][k][j][i] =  d->Vc[VX2][k][j][2*IBEG - i - 1];  ,
          d->Vc[VX3][k][j][i] =  d->Vc[VX3][k][j][2*IBEG - i - 1];  
        )
        d->Vc[PRS][k][j][i] =  d->Vc[PRS][k][j][2*IBEG - i - 1];
        #if PHYSICS == MHD
         EXPAND(
           d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][2*IBEG - i - 1];  ,
           d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][2*IBEG - i - 1];  ,
           d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][2*IBEG - i - 1];  
         )
        #endif
      }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    X1_END_LOOP(k,j,i){}
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    X2_BEG_LOOP(k,j,i){}
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    X2_END_LOOP(k,j,i){}
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    X3_BEG_LOOP(k,j,i){}
  }
 
  if (side == X3_END) {  /* -- X3_END boundary -- */
    X3_END_LOOP(k,j,i){}
  }
}

