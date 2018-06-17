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

  In order to solve just the internal energy equation, we force the velocity
  to be zero using the internal boundary.
  We then deal with the internal energy equation (in cgs units) in presence of
  conduction,
  \f[
  \label{eq:internal_energy}
  \pd{(\rho\epsilon)}{t} = \pd{}{x}\left[\kappa T^{5/2}\pd{T}{x}\right]
  \f]
  where \f$\rho\epsilon\f$ is the gas internal energy (in cgs units), \f$T\f$ is the
  temperature (in K), \f$\kappa = 9.22\cdot 10^{-7}\f$ is the plasma thermal
  conductivity (in cgs units):
  The previous equation can be written in terms of the temperature variable
  only using
  \f[
  \rho\epsilon = \frac{p}{\gamma-1} = \frac{2nk_BT}{\gamma-1}
  \qquad\qquad\Longrightarrow\qquad\qquad
  \pd{T}{x} = a \pd{}{x}\left(T^{5/2}\pd{T}{x}\right)
  \f]
  where \f$n=10^{10} \,{\rm cm}^{-3}\f$ (\f$10^{9}\f$ is employed in 1D)
  is the plasma number density, \f$\gamma=5/3\f$ while
  \f[
   a = \frac{(\gamma-1)\kappa}{2nk_B} \approx 0.226
   \qquad (Titos\quad uses\quad 4.212 \quad )
  \f]

  In code (non-dimensional) units, we adopt \f$t = t_0 \tilde{t}\f$,
  \f$x = L_0\tilde{x}\f$ and therefore Eq. (\ref{eq:internal_energy}) is solved
  as

  \f[\label{eq:internal_energy_code}
    \pd{(\widetilde{\rho\epsilon})}{\tilde{t}} 
   = \left(\kappa \frac{t_0}{L_0^2}\frac{T_0^{5/2}}{2nk_B}  \right)
     \pd{}{\tilde{x}}\left[\tilde{T}^{5/2}
                           \pd{\tilde{T}}{\tilde{x}}\right]
  \f]

  Values with a tilde are dimensionless.
  Note that the reference velocity and temperature are computed directly from
  \f$ t_0 = 1\,{\rm s}$ and $L_0 = 10^8\,{\rm cm} \f$:
  \f[
  v_0 = \frac{L_0}{t_0}\,,\qquad
  T_0 = \frac{m}{2k_B}v_0^2
  \f]

  \image html 1D-hd-exp.jpg "Evolution of the TC front in the 1D case; numerical solution (points) and analytical (lines)".

  Configurations (\c EXPLICIT, \c STS, \c RKL):
  - #01-#03: 1D Cartesian
  - #04-#06: 2D Cylindrical
  - #07-#09: 3D Cartesian
  - #10-#12: 3D Polar
  - #13-#15: 2D Spherical
  - #16      2D cylindrical + chombo

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 11, 2018

  \b References:
    - Section 4 of Reale, F. Computer Physics Communications (1995), 86, 13

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static double TProfile (double, double, double, double);

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *********************************************************************** */
{
  double tbeg = g_inputParam[TBEG];
  double T  = TProfile (x1, x2, x3, tbeg);
  double T0 = TREF;
  
  v[RHO] = 1.0;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  v[PRS] = v[RHO]*T/T0;
//print ("Init: (x1,x2) = %f, %f, p = %12.6e\n",x1,x2, v[PRS]);
  v[TRC] = 0.0;
  
  #if PHYSICS == MHD
  v[BX1] = 1.0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;
  #endif
}

/* ********************************************************************* */
double TProfile (double x1, double x2, double x3, double t)
/*
 *
 *********************************************************************** */
{
  double n, p, a, Q, t0, Tc, T, arg;
  double kpar = g_inputParam[KPAR];
  double nref = g_inputParam[NREF];

  a  = (g_gamma - 1.0)*kpar/(2.0*nref*CONST_kB);
#if DIMENSIONS == 1
  double zeta0, sf;
  double G1_2p1_n, G1_n;
  
  n = 2.5;
  p = 1.0/(n + 2.0);
  Q = 1.2e15;

  G1_2p1_n = 1.068628702119319;
  G1_n     = 2.218159543757688;

  zeta0 = pow(n+2,1+n)*pow(2.,1-n)/(n*pow(CONST_PI,n/2.0));
  zeta0 = pow(zeta0, p)*pow(G1_2p1_n/G1_n, n*p);
  sf    = pow(a*pow(Q,n)*t, p)*zeta0;
  Tc    = pow(Q*Q/(a*t), p)*pow(0.5*n*zeta0*zeta0*p, 1.0/n);

  sf  /= UNIT_LENGTH;
  arg  = 1.0 - x1*x1/(sf*sf);
  T    = Tc*pow(MAX(arg,1.e-20), 1.0/n);
  
  return T;
#else  /* The following is valid in 3D space */
  double zeta1, rf, r;
  double G5_2p1_n, G1p1_n, G3_2, x;

  n = 2.5;
  p = 1.0/(3.0*n + 2.0);
  Q = 1.e30;

  G5_2p1_n = 1.827355080624036;
  G1p1_n   = 0.887263817503075;
  G3_2     = 0.886226925452758;

  zeta1  = (3.0*n + 2.0)/(pow(2.0, n - 1.0)*n*pow(CONST_PI,n));
  zeta1 *= pow(G5_2p1_n/G1p1_n/G3_2,n);
  zeta1  = pow(zeta1,p);
  
  rf = pow(a*pow(Q,n)*t, p)*zeta1;
  Tc = zeta1*zeta1*zeta1/(rf*rf*rf)*Q*pow(n*zeta1*zeta1*p/2.0, 1.0/n);

  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
  r = sqrt(D_EXPAND(x1*x1, + x2*x2, + x3*x3));
  #elif GEOMETRY == POLAR
  r = sqrt(D_EXPAND(x1*x1, + 0.0, + x3*x3));
  #elif GEOMETRY == SPHERICAL
  r = x1;
  #endif

  rf  /= UNIT_LENGTH;
  arg  = 1.0 - r*r/(rf*rf);
  T    = Tc*pow(MAX(arg,1.e-18), 1.0/n);

//print ("x1,x2 = (%f, %f), Tc = %12.6e, arg = %12.6e, T = %12.6e\n",x1,x2,Tc, arg,T);
  return T;
#endif
}  
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 * Generate Analytical solution 
 *
 *********************************************************************** */
{
  static int count = 0; 
  char fname[32];
  double tbeg = g_inputParam[TBEG];
  double T, s;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];
  FILE *fp;
  
  if (prank == 0){
    sprintf (fname,"tcfront%dD_%02d.dat",DIMENSIONS, count);
    fp = fopen(fname,"w");
    for (s = 0.01; s < 5.0; s += 0.01){
      T = TProfile (s, 0.0, 0.0, g_time + tbeg);
      fprintf (fp,"%f  %12.6e\n", s, T/TREF);
    }
    fclose(fp);
    count++;
  }  
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
        x1 = grid->x[IDIR][i];
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

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}


