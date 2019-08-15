#include "pluto.h"
#include "cooling_defs.h"

/* ********************************************************************* */
void Solve_System (Ion *X, double Ne, double T)
/*!
 *  Solve sytem  A*x = rhs ( adapted from the GSL library ).
 *  This is used to compute the level populations for the ions.
 *
 *********************************************************************** */
{
  int    i, j, k, nlev;
  double scrh,tmpx, d;
  double **q;
  double **M, *rhs;
  int *p;  
  
  nlev = X->nlev;
  M    = ARRAY_2D(nlev, nlev, double);
  rhs  = ARRAY_1D(nlev, double);
  p    = ARRAY_1D(nlev, int);
  q    = ARRAY_2D (nlev, nlev, double);

  for (i = 0; i < nlev; i++)  {
    rhs[i] = 0.0;
    p[i] = 0.0;
    for (j = 0; j < nlev; j++)  {
      M[i][j] = 0.0;
      q[i][j] = 0.0;
    }
  }

/* -------------------------------------------
                Define qij
   -------------------------------------------  */

  for (i = 0; i < nlev; i++) { 
  for (j = i + 1; j < nlev; j++)  {
    if ( (!X->isMAP) && (!X->isCV) && (!X->isH) && (!X->isCHEB) )
      scrh = lagrange(X->Tom, X->omega[i][j], T, X->nTom, i, j); /* remove i and j after testing */
    if (X->isCV)  scrh = X->omega[i][j][0] * pow( (T / pow(10.,X->omega[i][j][2])), X->omega[i][j][1]);
    if (X->isMAP) scrh = X->omega[i][j][0] * pow( T/10000., X->omega[i][j][1]);
    if (X->isH)  {
      if (T<55000) scrh = X->omega[i][j][0] + T*X->omega[i][j][1] + T*T*X->omega[i][j][2] + T*T*T*X->omega[i][j][3];
      else scrh = X->omega[i][j][4] + T*X->omega[i][j][5] + T*T*X->omega[i][j][6] + T*T*T*X->omega[i][j][7];
    }
    if (X->isCHEB) {
      tmpx = 0.6199646 * log(T) - 6.2803580;
      scrh = 0.5 * X->omega[i][j][0] + X->omega[i][j][1] * tmpx + X->omega[i][j][2] * (2.*tmpx*tmpx - 1) + X->omega[i][j][3] * (4.*tmpx*tmpx*tmpx - 3*tmpx);
      scrh = exp(scrh);
    }

    q[j][i]  = scrh/X->wght[j];
    q[j][i] *= 8.629e-6/sqrt(T);
    q[i][j]  = (X->wght[j]/X->wght[i])*q[j][i]*exp( - X->dE[i][j]/(kB*T)); 
  }}               
  for (i = 0; i < nlev; i++) q[i][i] = 0.0;

/* -----------------------------------------
       compute coefficient matrix 
   ----------------------------------------- */

  for (i = 0 ; i < nlev ; i++) {
  for (j = 0 ; j < nlev ; j++) {
    scrh = 0.0; 
    if (j == i) {
      for (k = 0 ; k < nlev ; k++) {
          if (k != i) scrh += Ne*q[i][k];
          if (k < i)  scrh += X->A[i][k];
      }
      M[i][j] = -scrh;
    }
    if (j < i) M[i][j] = Ne*q[j][i];
    if (j > i) M[i][j] = Ne*q[j][i] + X->A[j][i];
  }}

/* -------------------------------------
    Replace 1st eq. with normalization
    condition and define rhs 
   ------------------------------------- */

  for (j = 0; j < nlev; j++) {
    M[nlev-1][j] = 1.0;
    rhs[j] = 0.0;
  } 
  rhs[nlev-1] = 1.0;

/* ----------------------------------
    Solve equations by LU decomp
   ---------------------------------- */
  LUDecompose (M, nlev, p, &d);
  LUBackSubst (M, nlev, p, rhs);
  for (i = 0 ; i < nlev ; i++) X->Ni[i] = rhs[i];
}

/* ********************************************************************* */
void Symmetrize_Coeff (Ion *X)
/*
 *
 *
 * 
 *********************************************************************** */
{
  int i,j,n;
  
/*  -------------------------------------------------- 
         symmetrize dE, A, and Omega   
         in order to prevent swapped index in the
         ion routines  
    -------------------------------------------------- */

  for (i = 0; i < X->nlev; i++) {
  for (j = 0; j < X->nlev; j++) {
    if (X->dE[i][j] > 0.0) {
      X->dE[j][i] = X->dE[i][j];  
      X->A[j][i]  = X->A[i][j];     
      for (n = 0; n < X->nTom; n++){          
        if (   (X->omega[i][j][n] > 0.0) 
            || (X->isH    && X->omega[i][j][n] != 0.0) 
            || (X->isCHEB && X->omega[i][j][n] != 0.0) )
          X->omega[j][i][n]  = X->omega[i][j][n];
      }
    }
  }}

/* ----------------------------------------------------
    2. convert energy to correct units
   ---------------------------------------------------- */

  for (i = 0; i < X->nlev; i++){ 
  for (j = 0; j < X->nlev; j++) {
    if (X->dE[i][j] < 0.0) { X->dE[i][j] = 0.0; X->A[i][j] = 0.0; for (n = 0; n < X->nTom; n++) X->omega[i][j][n] = 0.000; }
    if (X->A[i][j] < 0.0) { X->A[i][j] = 0.0; for (n = 0; n < X->nTom; n++) X->omega[i][j][n] = 0.000; }
    if (X->dE[i][j] > 0.0) X->dE[i][j] = 1.2399e-4/(1.e-8*X->dE[i][j])  ;  /* multply by hc to get energy  */
    if ( (X->omega[i][j][0] < 0.0)&&(!X->isCV)&&(!X->isH)&&(!X->isCHEB)&&(!X->isMAP) ) for (n = 0; n < X->nTom; n++) X->omega[i][j][n] = 0.000; 
  }}


/* -------------------------------------------------- 
     3. set A(i,j) = 0  when j >= i
   -------------------------------------------------- */

  for (i = 0; i < X->nlev; i++) {
  for (j = i; j < X->nlev; j++) {
    X->A[i][j] = 0.0;                     
  }}

/* --------------------------------------------------- 
    4. set diagonal elements of dE and omega to zero
   --------------------------------------------------- */

  for (i = 0; i < X->nlev; i++) {
    for (n = 0; n < X->nTom; n++) X->omega[i][i][n] = 0.000;
    X->dE[i][i]    = 0.0;
    X->A[i][i]     = 0.0;
  }

}

/* ********************************************************************* */
double lagrange (double *x, double *y, double xp, int n, int ii, int jj)
/*
 * PURPOSE:
 *
 *   Return lagrange interpolation for a tabulated function.
 *
 *  x : a vector of n-components with the tabulated values
 *      of the independent variables x;
 *
 *  y : a vecotr of n-components with the tabulated values
 *      y[n] = y(x[n]);
 *
 * xp: the point at which interpolation is desired
 *
 * n : the number of points in the table, also
 *     the degree of the intrpolant
 *
 *
 ********************************************************** */
{
  int    i, k, j;
  double scrh, yp;

  yp = 0.0;
  if (xp < x[0])    return (y[0]);
  if (xp > x[n-1])  return (y[n-1]);

  for (i = 0; i < n; i++){
    scrh = 1.0;
    for (k = 0; k < n; k++){
      if (k == i) continue;
      scrh *= (xp - x[k])/(x[i] - x[k]);
    }
    yp += scrh*y[i];
  }
  return(yp);
}

/* ********************************************************************* */
int Create_Ion_Coeff_Tables(double ***tbl)
/*!
 *   Compute and save to memory the coefficients depending only
 *   on temperature used at runtime to compute the ionization balance.
 *
 *********************************************************************** */
#define ZERO_5X  0.0, 0.0, 0.0, 0.0, 0.0
{
  double T, t4, tmprec, lam, f1, f2, x1, x2, P, Q, F[2], Tn[9];

  long int nv, i, j, ti1, ti2, cindex;
  
  static double **coll_ion, **rad_rec, **diel_rec, **chtr_hp, **chtr_h, **chtr_he, **tot_recs; 

/* Collisional ionization rates , fit parameters - Voronov 1997, ATOMIC DATA AND NUCLEAR DATA TABLES 65, 1-35  */

double coll_ion_P[] = { 0., 0., 1.
                        C_EXPAND(0., 1., 1., 1., 1.)
                        N_EXPAND(0., 0., 1., 1., 0.) 
                        O_EXPAND(0., 1., 1., 0., 1.) 
                       Ne_EXPAND(1., 0., 1., 1., 1.) 
                        S_EXPAND(1., 1., 1., 1., 1.)
                       Fe_EXPAND(0., 1., 0.) }; 
double coll_ion_A[] = { 0.291e-7, 0.175e-7, 0.0000 
                        C_EXPAND(0.685e-7, 0.185e-7, 0.635e-8, 0.150e-8, 0.0000)
                        N_EXPAND(0.482e-7, 0.298e-7, 0.810e-8, 0.371e-8, 0.0000) 
                        O_EXPAND(0.359e-7, 0.139e-7, 0.931e-8, 0.102e-7, 0.0000)
                       Ne_EXPAND(0.150e-7, 0.198e-7, 0.703e-8, 0.424e-8, 0.0000) 
                        S_EXPAND(0.549e-7, 0.681e-7, 0.214e-7, 0.166e-7, 0.0000)
                       Fe_EXPAND(0.252e-6, 0.221e-7, 0.0000) };  /*  in cm^3/s  */
double coll_ion_X[] = { 0.232, 0.180, 0.265 
                        C_EXPAND( 0.193, 0.286,  0.427,  0.416, 0.0000) 
                        N_EXPAND(0.0652, 0.310,  0.350,  0.549, 0.0167) 
                        O_EXPAND( 0.073, 0.212,  0.270,  0.614, 0.630)
                       Ne_EXPAND(0.0329, 0.295, 0.0677, 0.0482, 0.0000) 
                        S_EXPAND( 0.100, 0.693,  0.353,  1.030, 0.0000)
                       Fe_EXPAND(0.7010, 0.033, 0.0000) };
double coll_ion_K[] = { 0.39, 0.35, 0.25 
                        C_EXPAND(0.25, 0.24, 0.21, 0.13, 0.02) 
                        N_EXPAND(0.42, 0.30, 0.24, 0.18, 0.74) 
                        O_EXPAND(0.34, 0.22, 0.27, 0.27, 0.17)
                       Ne_EXPAND(0.43, 0.20, 0.39, 0.58, 0.00) 
                        S_EXPAND(0.25, 0.21, 0.24, 0.14, 0.00)
                       Fe_EXPAND(0.25, 0.45, 0.17) };
double coll_ion_Tmin[] = { 0., 0., 0.
                           C_EXPAND(0., 0., 0., 0., 0.)
                           N_EXPAND(0., 0., 0., 0., 0.) 
                           O_EXPAND(0., 0., 0., 0., 0.) 
                          Ne_EXPAND(0., 0., 0., 0., 0.) 
                           S_EXPAND(0., 0., 0., 0., 0.)
                          Fe_EXPAND(0., 0., 0.) };   /*  in eV   */

/*
double coll_ion_Tmin[28] = { 1., 1., 3., 1., 1., 2., 3., 20., 
                        1., 2., 2., 4., 4., 1., 2., 3., 4., 5., 
                        1., 3., 3., 5., 7., 1., 1., 2., 2., 3. };   */


/* Charge transfer with H II ionization rates, fit parameters - Kingdon & Ferland 1996, ApJSS  */

double chtrH_ion_a[] = { 0.00, 0.00, 0.00 
                         C_EXPAND(   0.00, 0.00, 0.00, 0.00, 0.00) 
                         N_EXPAND(4.55e-3, 0.00, 0.00, 0.00, 0.00) 
                         O_EXPAND( 7.4e-2, 0.00, 0.00, 0.00, 0.00) 
                        Ne_EXPAND(   0.00, 0.00, 0.00, 0.00, 0.00) 
                         S_EXPAND(   0.00, 0.00, 0.00, 0.00, 0.00)
                        Fe_EXPAND(   5.40, 2.10, 0.00) };  /*  in 10^(-9) cm^3/s */
double chtrH_ion_b[] = { 0.00, 0.00, 0.00 
                         C_EXPAND( 0.00, 0.00, 0.00, 0.00, 0.00) 
                         N_EXPAND(-0.29, 0.00, 0.00, 0.00, 0.00) 
                         O_EXPAND( 0.47, 0.00, 0.00, 0.00, 0.00) 
                        Ne_EXPAND( 0.00, 0.00, 0.00, 0.00, 0.00) 
                         S_EXPAND( 0.00, 0.00, 0.00, 0.00, 0.00)
                        Fe_EXPAND( 0.00, 7.72e-2, 0.00) };
double chtrH_ion_c[] = { 0.00, 0.00, 0.00 
                         C_EXPAND( 0.00, 0.00, 0.00, 0.00, 0.00) 
                         N_EXPAND(-0.92, 0.00, 0.00, 0.00, 0.00) 
                         O_EXPAND(24.37, 0.00, 0.00, 0.00, 0.00) 
                        Ne_EXPAND( 0.00, 0.00, 0.00, 0.00, 0.00) 
                         S_EXPAND( 0.00, 0.00, 0.00, 0.00, 0.00)
                        Fe_EXPAND( 0.00, -0.41, 0.00)};
double chtrH_ion_d[] = { 0.00, 0.00, 0.00 
                         C_EXPAND( 0.00, 0.00, 0.00, 0.00, 0.00) 
                         N_EXPAND(-8.38, 0.00, 0.00, 0.00, 0.00) 
                         O_EXPAND(-0.74, 0.00, 0.00, 0.00, 0.00) 
                        Ne_EXPAND( 0.00, 0.00, 0.00, 0.00, 0.00) 
                         S_EXPAND( 0.00, 0.00, 0.00, 0.00, 0.00)
                        Fe_EXPAND( 0.00, -7.31, 0.00) };
double chtrH_ion_upT[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(1.0,  0.00, 0.00, 0.00, 0.00) 
                           N_EXPAND(5.0,  0.00, 0.00, 0.00, 0.00) 
                           O_EXPAND(1.0,  0.00, 0.00, 0.00, 0.00) 
                          Ne_EXPAND(0.00, 0.00, 0.00, 0.00, 0.00) 
                           S_EXPAND(1.0,  0.00, 0.00, 0.00, 0.00)
                          Fe_EXPAND(1.0,  10.0, 0.0) };  /*   in 10^4 K   */

/* Charge transfer with H recombination rates, fit parameters - Kingdon & Ferland 1996, ApJSS  */

double chtrH_rec_a[] = { 0.00, 7.46e-6, 0.00 
                         C_EXPAND(1.76e-9, 1.67e-4,  3.25, 332.46, 0.00) 
                         N_EXPAND(1.01e-3, 3.05e-1,  4.54,   3.28, 0.00) 
                         O_EXPAND(   1.04,    1.04,  3.98,  0.252, 0.00)
                        Ne_EXPAND(   0.00, 1.00e-5, 14.73,   6.47, 0.00) 
                         S_EXPAND(3.82e-7, 1.00e-5,  2.29,   6.44, 0.00)
                        Fe_EXPAND(   0.00,    1.26, 0.00) }; 
        /*  in  10^(-9) cm^3/s , per Ion - 1 (e.g., the coeffs for He II recomb to He I are written to the He I position */
double chtrH_rec_b[] = { 0.00, 2.06, 0.00 
                         C_EXPAND(   8.33,    2.79,    0.21, -0.11, 0.00)
                         N_EXPAND(  -0.29,    0.60,    0.57,  0.52, 0.00) 
                         O_EXPAND(3.15e-2,    0.27,    0.26,  0.63, 0.00) 
                        Ne_EXPAND(   0.00,    0.00, 4.52e-2,  0.54, 0.00) 
                         S_EXPAND(  11.10,    0.00, 4.02e-2,  0.13, 0.00)
                        Fe_EXPAND(   0.00, 7.72e-2, 0.00) };
double chtrH_rec_c[] = { 0.00, 9.93, 0.00 
                         C_EXPAND(4278.78, 304.72,  0.19, -0.995, 0.00)
                         N_EXPAND(  -0.92,   2.65, -0.65,  -0.52, 0.00) 
                         O_EXPAND(  -0.61,   2.02,  0.56,   2.08, 0.00) 
                        Ne_EXPAND(   0.00,   0.00, -0.84,   3.59, 0.00) 
                         S_EXPAND(2.57e+4,   0.00,  1.59,   2.69, 0.00)
                        Fe_EXPAND(   0.00,  -0.41,  0.00) };
double chtrH_rec_d[] = { 0.00, -3.89, 0.00 
                         C_EXPAND(-6.41, -4.07, -3.29, -1.58e-3, 0.00)
                         N_EXPAND(-8.38, -0.93, -0.89,    -0.19, 0.00) 
                         O_EXPAND(-9.73, -5.92, -2.62,    -4.16, 0.00) 
                        Ne_EXPAND( 0.00,  0.00, -0.31,    -5.22, 0.00) 
                         S_EXPAND(-8.22,  0.00, -6.06,    -5.69, 0.00)
                        Fe_EXPAND( 0.00, -7.31,  0.00) };
double chtrH_rec_upT[] = { 0.00, 1.e+5, 1.e+7 
                           C_EXPAND(1.e+4, 5.e+4, 1.e+5, 1.e+5, 0.00)
                           N_EXPAND(5.e+4, 1.e+5, 1.e+5, 3.e+5, 0.00) 
                           O_EXPAND(1.e+4, 1.e+5, 5.e+4, 3.e+4, 0.00)
                          Ne_EXPAND( 0.00, 5.e+4, 5.e+4, 3.e+4, 0.00) 
                           S_EXPAND(1.e+4, 3.e+4, 3.e+4, 3.e+4, 0.00)
                          Fe_EXPAND( 0.00, 1.e+5, 0.00) };  /*  in K   */

/* Charge transfer of O ions with H recombination rates, fit parameters - 
   M. Rakovic', J. G. Wang, D. R. Schultz, and P. C. Stancil (2001)  
   ORNL Charge Transfer Database                                                          */

double chtrH_rec_o1[8] = {3.348706E+00, 1.281922E+00, -1.870663E+00,-7.195090E-01,-2.473806E-01, 8.577410E-02, 1.721933E-02, 9.213735E-03};
double chtrH_rec_o2[8] = {3.266558E+00, 2.264383E+00, -1.019642E+00,-1.090796E+00,-4.403829E-01,-1.369024E-01, 1.606288E-01, 1.925501E-01};
double chtrH_rec_o3[8] = {4.860561E+00, 1.680463E+00, -1.110986E+00,-1.346895E+00,-2.697744E-01, 1.329470E-01,-4.511763E-02, 4.614949E-02};
double chtrH_rec_o4[8] = {4.680553E+00, 2.456278E+00, -1.026270E+00,-1.587952E+00,-1.318074E-01, 1.451866E-01,-1.526627E-01, 7.596889E-02};
double chtrH_rec_o5[8] = {5.723788E+00, 1.771316E+00, -1.254652E+00,-9.009546E-01,-3.558194E-01 -7.138940E-02, 2.846941E-02, 1.313952E-01};
double chtrH_rec_Tmin = 1.0e+2, chtrH_rec_Tmax = 1.0e+10;

/* Charge transfer with He  -  ORNL Charge Transfer Database 
                               http://www-cfadc.phy.ornl.gov/astro/ps/data/cx/helium/rates/fits.data   */ 

/* 1st stage - up to 5000K */
double chtrHe_rec_a1[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00,    0.00,    1.12, 3.12e-7, 0.00) 
                           N_EXPAND(0.00, 4.84e-1,    2.05, 1.26e-2, 0.00) 
                           O_EXPAND(0.00, 7.10e-3,    1.12,   0.997, 0.00)
                          Ne_EXPAND(0.00, 1.00e-5, 1.00e-5,    1.77, 0.00) 
                           S_EXPAND(0.00,    0.00,    3.58, 7.44e-4, 0.00)
                          Fe_EXPAND(0.00,    0.00, 0.00) };/*  in 10^(-9) cm^3/s */
double chtrHe_rec_b1[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00, 0.00,    0.42, -7.37e-2, 0.00)
                           N_EXPAND(0.00, 0.92,    0.23,     1.55, 0.00) 
                           O_EXPAND(0.00, 2.60,    0.42,     0.40, 0.00)
                          Ne_EXPAND(0.00, 0.00,    0.00,     0.14, 0.00) 
                           S_EXPAND(0.00, 0.00, 7.77e-3,     0.34, 0.00)
                          Fe_EXPAND(0.00, 0.00,   0.00) };
double chtrHe_rec_c1[] = { 0.00, 0.00, 0.00
                           C_EXPAND(0.00, 0.00, -0.69, 3.50e+1, 0.00)
                           N_EXPAND(0.00, 2.37, -0.72,    11.2, 0.00) 
                           O_EXPAND(0.00, 8.99, -0.71,   -0.46, 0.00)
                          Ne_EXPAND(0.00, 0.00,  0.00, 4.88e-2, 0.00) 
                           S_EXPAND(0.00, 0.00, -0.94,    3.74, 0.00)
                          Fe_EXPAND(0.00, 0.00,  0.00) };
double chtrHe_rec_d1[] = { 0.00, 0.00, 0.00  
                           C_EXPAND(0.00,  0.00,    -0.34,  2.40, 0.00)
                           N_EXPAND(0.00, -10.2,    -0.19, -7.82, 0.00) 
                           O_EXPAND(0.00, -0.78, -1.98e-2, -0.35, 0.00)
                          Ne_EXPAND(0.00,  0.00,     0.00, -3.35, 0.00) 
                           S_EXPAND(0.00,  0.00,    -0.30, -5.18, 0.00)
                          Fe_EXPAND(0.00,  0.00,     0.00) };
/* 2nd stage - 5000K to 10000K */
double chtrHe_rec_a2[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00,    0.00,    1.12, 3.12e-7, 0.00) 
                           N_EXPAND(0.00, 4.84e-1,    2.05, 1.26e-2, 0.00) 
                           O_EXPAND(0.00, 7.10e-3,    1.12,   0.997, 0.00)
                          Ne_EXPAND(0.00, 8.48e-3, 1.00e-5,    1.77, 0.00) 
                           S_EXPAND(0.00,    0.00,    3.58, 7.44e-4, 0.00)
                          Fe_EXPAND(0.00,    0.00, 0.00) };
double chtrHe_rec_b2[] = { 0.00, 0.00, 0.00
                           C_EXPAND(0.00, 0.00,    0.42, -7.37e-2, 0.00)
                           N_EXPAND(0.00, 0.92,    0.23,     1.55, 0.00) 
                           O_EXPAND(0.00, 2.60,    0.42,     0.40, 0.00)
                          Ne_EXPAND(0.00, 3.35,    0.00,     0.14, 0.00) 
                           S_EXPAND(0.00, 0.00, 7.77e-3,     0.34, 0.00)
                          Fe_EXPAND(0.00, 0.00, 0.00) };
double chtrHe_rec_c2[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00,  0.00, -0.69, 3.50e+1, 0.00)
                           N_EXPAND(0.00,  2.37, -0.72,    11.2, 0.00) 
                           O_EXPAND(0.00,  8.99, -0.71,   -0.46, 0.00)
                          Ne_EXPAND(0.00, -1.92,  0.00, 4.88e-2, 0.00) 
                           S_EXPAND(0.00,  0.00, -0.94,    3.74, 0.00)
                          Fe_EXPAND(0.00,  0.00,  0.00) };
double chtrHe_rec_d2[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00,  0.00,    -0.34,  2.40, 0.00)
                           N_EXPAND(0.00, -10.2,    -0.19, -7.82, 0.00) 
                           O_EXPAND(0.00, -0.78, -1.98e-2, -0.35, 0.00) 
                          Ne_EXPAND(0.00, -1.50,     0.00, -3.35, 0.00) 
                           S_EXPAND(0.00,  0.00,    -0.30, -5.18, 0.00)
                          Fe_EXPAND(0.00,  0.00, 0.00) };
/* 3rd stage - 10,000K to 50,000K */
double chtrHe_rec_a3[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00,    0.00,    1.12, 1.49e-5, 0.00) 
                           N_EXPAND(0.00, 4.84e-1,    2.05, 1.26e-2, 0.00) 
                           O_EXPAND(0.00, 7.10e-3,    1.12,   0.997, 0.00)
                          Ne_EXPAND(0.00, 2.52e-2, 1.34e-4,    1.77, 0.00) 
                           S_EXPAND(0.00,    0.00,    3.58, 7.44e-4, 0.00)
                          Fe_EXPAND(0.00,    0.00, 0.00) };
double chtrHe_rec_b3[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00, 0.00,    0.42, 2.73, 0.00)
                           N_EXPAND(0.00, 0.92,    0.23, 1.55, 0.00) 
                           O_EXPAND(0.00, 2.60,    0.42, 0.40, 0.00)
                          Ne_EXPAND(0.00, 0.14,    2.33, 0.14, 0.00) 
                           S_EXPAND(0.00, 0.00, 7.77e-3, 0.34, 0.00)
                          Fe_EXPAND(0.00, 0.00, 0.00) };
double chtrHe_rec_c3[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00,  0.00, -0.69,    5.93, 0.00)
                           N_EXPAND(0.00,  2.37, -0.72,    11.2, 0.00) 
                           O_EXPAND(0.00,  8.99, -0.71,   -0.46, 0.00)
                          Ne_EXPAND(0.00, -1.99, -2.55, 4.88e-2, 0.00) 
                           S_EXPAND(0.00,  0.00, -0.94,    3.74, 0.00)
                          Fe_EXPAND(0.00, 0.00, 0.00) };
double chtrHe_rec_d3[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00,  0.00,    -0.34, -8.74e-2, 0.00)
                           N_EXPAND(0.00, -10.2,    -0.19,    -7.82, 0.00) 
                           O_EXPAND(0.00, -0.78, -1.98e-2,    -0.35, 0.00) 
                          Ne_EXPAND(0.00, -0.91,    -0.37,    -3.35, 0.00) 
                           S_EXPAND(0.00,  0.00,    -0.30,    -5.18, 0.00)
                          Fe_EXPAND(0.00, 0.00, 0.00) };
/* 4th stage - 50,000K to 100,000K */
double chtrHe_rec_a4[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00,    0.00,    1.12, 1.49e-5, 0.00) 
                           N_EXPAND(0.00,    3.17,    2.05, 1.26e-2, 0.00) 
                           O_EXPAND(0.00, 6.21e-1,    1.12,   0.997, 0.00)
                          Ne_EXPAND(0.00, 2.52e-2, 1.34e-4, 2.67e-1, 0.00) 
                           S_EXPAND(0.00,    0.00,    0.00,    0.00, 0.00)
                          Fe_EXPAND(0.00,    0.00, 0.00) };
double chtrHe_rec_b4[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00, 0.00, 0.42, 2.73, 0.00)
                           N_EXPAND(0.00, 0.20, 0.23, 1.55, 0.00) 
                           O_EXPAND(0.00, 0.53, 0.42, 0.40, 0.00)
                          Ne_EXPAND(0.00, 0.14, 2.33, 0.54, 0.00) 
                           S_EXPAND(0.00, 0.00, 0.00, 0.00, 0.00)
                          Fe_EXPAND(0.00, 0.00, 0.00) };
double chtrHe_rec_c4[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00,  0.00, -0.69,  5.93, 0.00)
                           N_EXPAND(0.00, -0.72, -0.72,  11.2, 0.00) 
                           O_EXPAND(0.00, -0.66, -0.71, -0.46, 0.00)
                          Ne_EXPAND(0.00, -1.99, -2.55,  0.91, 0.00) 
                           S_EXPAND(0.00,  0.00,  0.00,  0.00, 0.00)
                          Fe_EXPAND(0.00, 0.00, 0.00) };
double chtrHe_rec_d4[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00,     0.00,    -0.34, -8.74e-2, 0.00)
                           N_EXPAND(0.00, -4.81e-2,    -0.19,    -7.82, 0.00) 
                           O_EXPAND(0.00, -2.22e-2, -1.98e-2,    -0.35, 0.00) 
                          Ne_EXPAND(0.00,    -0.91,    -0.37, -1.88e-2, 0.00) 
                           S_EXPAND(0.00,     0.00,     0.00,     0.00, 0.00)
                          Fe_EXPAND(0.00,     0.00, 0.00) };
/* 5th stage - more than 100,000K */
double chtrHe_rec_a5[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00,    0.00, 1.12, 1.49e-5, 0.00) 
                           N_EXPAND(0.00,    3.17, 2.05, 3.75e-1, 0.00) 
                           O_EXPAND(0.00, 6.21e-1, 1.12,   0.997, 0.00)
                          Ne_EXPAND(0.00, 2.52e-2,  0.1, 2.67e-1, 0.00) 
                           S_EXPAND(0.00,    0.00, 0.00,    0.00, 0.00)
                          Fe_EXPAND(0.00, 0.00, 0.00) };
double chtrHe_rec_b5[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00, 0.00, 0.42, 2.73, 0.00)
                           N_EXPAND(0.00, 0.20, 0.23, 0.54, 0.00) 
                           O_EXPAND(0.00, 0.53, 0.42, 0.40, 0.00)
                          Ne_EXPAND(0.00, 0.14, 0.24, 0.54, 0.00) 
                           S_EXPAND(0.00, 0.00, 0.00, 0.00, 0.00)
                          Fe_EXPAND(0.00, 0.00, 0.00) };
double chtrHe_rec_c5[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00,  0.00, -0.69,  5.93, 0.00)
                           N_EXPAND(0.00, -0.72, -0.72, -0.82, 0.00) 
                           O_EXPAND(0.00, -0.66, -0.71, -0.46, 0.00)
                          Ne_EXPAND(0.00, -1.99, -1.09,  0.91, 0.00) 
                           S_EXPAND(0.00,  0.00,  0.00,  0.00, 0.00)
                          Fe_EXPAND(0.00, 0.00, 0.00) };
double chtrHe_rec_d5[] = { 0.00, 0.00, 0.00 
                           C_EXPAND(0.00,     0.00,    -0.34, -8.74e-2, 0.00)
                           N_EXPAND(0.00, -4.81e-2,    -0.19, -2.07e-2, 0.00) 
                           O_EXPAND(0.00, -2.22e-2, -1.98e-2,    -0.35, 0.00) 
                          Ne_EXPAND(0.00,    -0.91, -2.47e-2, -1.88e-2, 0.00) 
                           S_EXPAND(0.00,     0.00,     0.00,     0.00, 0.00)
                          Fe_EXPAND(0.00, 0.00, 0.00) };

/* Radiative recombination - Pequignot & al 1991, A&A   */
double rad_rec_a[] = { 5.596, 0.0000, 0.0000 
                       C_EXPAND(5.068, 5.434, 4.742,  4.051, 0.00)
                       N_EXPAND(3.874, 4.974, 4.750,  4.626, 0.00) 
                       O_EXPAND(3.201, 4.092, 4.890, 14.665, 0.00)
                      Ne_EXPAND(0.000, 0.000, 0.000, 0.0000, 0.00) 
                       S_EXPAND(0.000, 0.000, 0.000, 0.0000, 0.00) };  /* H value NOT used */
                        /* Ne values: 11.800, 5.841, 15.550, 7.538, */
                        /* He I and II values: 8.295, 5.596,   */
double rad_rec_b[] = { -0.6038, 0.0000, 0.0000 
                       C_EXPAND(-0.6192, -0.6116, -0.6167, -0.6270, 0.00)
                       N_EXPAND(-0.6487, -0.6209, -0.5942, -0.9521, 0.00) 
                       O_EXPAND(-0.6880, -0.6413, -0.6213, -0.5140, 0.00)
                      Ne_EXPAND( 0.0000,  0.0000,  0.0000,  0.0000, 0.00) 
                       S_EXPAND( 0.0000,  0.0000,  0.0000,  0.0000, 0.00)};  /*   S, H and Ne to add separately ; prev. H value: -0.6038 */
                        /* Ne values: -0.5404, -0.5921, -0.4825, -0.5540, */
                        /* He I and II values:  -0.6193, -0.6038, */
double rad_rec_c[] = { 0.3436, 0.0000, 0.0000 
                       C_EXPAND(-0.0815, 0.0694, 0.2960, 0.5054, 0.00) 
                       N_EXPAND( 0.0000, 0.0000, 0.8452, 0.4729, 0.00) 
                       O_EXPAND(-0.0174, 0.0000, 0.0184, 2.7300, 0.00)
                      Ne_EXPAND( 0.0000, 0.0000, 0.0000, 0.0000, 0.00) 
                       S_EXPAND( 0.0000, 0.0000, 0.0000, 0.0000, 0.00)};/* prev. H value: 0.3436*/
                        /* Ne values: 3.0300, 0.4498, 3.2740, 1.2960, */
                        /* He I and II values:  0.9164, 0.3436,  */
double rad_rec_d[] = { 0.4479, 0.0000, 0.0000 
                       C_EXPAND(1.2910, 0.7866, 0.6167,  0.6692, 0.00) 
                       N_EXPAND(1.0000, 1.0000, 2.8450, -0.4477, 0.00) 
                       O_EXPAND(1.7070, 1.0000, 1.5550,  0.2328, 0.00)
                      Ne_EXPAND(0.0000, 0.0000, 0.0000,  0.0000, 0.00) 
                       S_EXPAND(0.0000, 0.0000, 0.0000,  0.0000, 0.00) };/* prev. H value: 0.4479*/
                        /* Ne values: 0.2050, 0.6395, 0.3030, 0.3472, */
                        /* He I and II values:  0.2667, 0.4479, */

/* Dielectronic recombination  - Nussbaumer, Storey, 1983*/
double diel_rec_a[] = { 0.0000, 0.0000, 0.0000 
                        C_EXPAND(0.0108,  1.8267,  2.3196, 0.0000, 0.00)
                        N_EXPAND(0.0000,  0.0320, -0.8806, 0.4134, 0.00) 
                        O_EXPAND(0.0000, -0.0036,  0.0000, 0.0061, 0.00)
                       Ne_EXPAND(0.0000,  0.0000,  0.0000, 0.0000, 0.00) 
                        S_EXPAND(0.0000,  0.0000,  0.0000, 0.0000, 0.00) };
double diel_rec_b[] = { 0.0000, 0.0000, 0.0000 
                        C_EXPAND(-0.1075,  4.1012, 10.7328,  0.0000, 0.00)
                        N_EXPAND( 0.6310, -0.6624, 11.2406, -4.6319, 0.00) 
                        O_EXPAND( 0.0238,  0.7519, 21.8790,  0.2269, 0.00)
                       Ne_EXPAND( 0.0000,  0.0000,  0.0000,  0.0000, 0.00) 
                        S_EXPAND( 0.0000,  0.0000,  0.0000,  0.0000, 0.00) };
double diel_rec_c[] = { 0.0000, 0.0000, 0.0000 
                        C_EXPAND(0.2810, 4.8443,  6.8830,  0.0000, 0.00)
                        N_EXPAND(0.1990, 4.3191, 30.7066, 25.9172, 0.00) 
                        O_EXPAND(0.0659, 1.5252, 16.2730, 32.1419, 0.00)
                       Ne_EXPAND(0.0000, 0.0000,  0.0000,  0.0000, 0.00) 
                        S_EXPAND(0.0000, 0.0000,  0.0000,  0.0000, 0.00) };
double diel_rec_d[] = { 0.0000, 0.0000, 0.0000 
                        C_EXPAND(-0.0193,  0.2261, -0.1824,  0.0000, 0.00)
                        N_EXPAND(-0.0197,  0.0003, -1.1721, -2.2290, 0.00) 
                        O_EXPAND( 0.0349, -0.0838, -0.7020,  1.9939, 0.00)
                       Ne_EXPAND( 0.0000,  0.0000,  0.0000,  0.0000, 0.00) 
                        S_EXPAND( 0.0000,  0.0000,  0.0000,  0.0000, 0.00) };
double diel_rec_f[] = { 0.0000, 0.0000, 0.0000 
                        C_EXPAND(-0.1127, 0.5960, 0.4101,  0.0000, 0.00)
                        N_EXPAND( 0.4398, 0.5946, 0.6127,  0.2360, 0.00) 
                        O_EXPAND( 0.5334, 0.2769, 1.1899, -0.0646, 0.00)
                       Ne_EXPAND( 0.0000, 0.0000, 0.0000,  0.0000, 0.00) 
                        S_EXPAND( 0.0000, 0.0000, 0.0000,  0.0000, 0.00) };


/*  Total recombination rate for H, Ne and S  -  NIFS-DATA-54 - Kato & Asano 1999 */
/* recombination coefficient to the X^n ion (from the X^(n+1) is written in array in the position of X^n */
double tot_rec_T[5] = { 5.e+3, 1.1e+4, 5.5e+4, 1.1e+5, 2.2e+5 };  /*  in K  */
static int tot_rec_ions[] = { 0, 1, 18, 19, 20, 21, 23, 24, 25, 26};
double tot_rec[][5] = { { -12.5, -13.0, -13.2, -13.4, -13.65 },    /*  H     */
                          { -12.5, -12.6, -12.8, -12.1, -11.75 },  /*  He I  */
                /*        { -11.4, -11.7, -12.2, -12.45, -12.8 },   */
                          { -12.3, -12.5, -12.6, -11.8, -11.7 },     /*  Ne I  */
                          { -11.7, -11.9, -11.6, -11.3, -11.3 },    
                          { -11.2, -11.4, -11.6, -10.9, -11.1 },
                          { -11.0, -11.1, -11.3, -10.6, -10.7 },
                          { -12.5, -12.4, -10.9, -10.95, -11.4 },    /*  S I    */
                          { -11.7, -11.8, -10.6, -10.3, -10.6 },
                          { -11.4, -11.6, -10.0,  -9.8, -10.1 },
                          { -11.2, -11.3,  -9.9,  -9.7, -10.0 } };   /*  log (Alpha_tot)  [cm^3/s]*/

double tot_recH_He[2][5]  = {{ -12.5, -13.0, -13.2, -13.4, -13.65 },   /* -- H -- */
                             { -12.5, -12.6, -12.8, -12.1, -11.75 }};  /* -- He -- */
double tot_recNe[][5] = { { -12.3, -12.5, -12.6, -11.8, -11.7 },     
                          { -11.7, -11.9, -11.6, -11.3, -11.3 },    
                          { -11.2, -11.4, -11.6, -10.9, -11.1 },
                          { -11.0, -11.1, -11.3, -10.6, -10.7 }};
double tot_recS[][5]  = { { -12.5, -12.4, -10.9, -10.95, -11.4 },   
                          { -11.7, -11.8, -10.6, -10.3, -10.6 },
                          { -11.4, -11.6, -10.0,  -9.8, -10.1 },
                          { -11.2, -11.3,  -9.9,  -9.7, -10.0 }};


/* Atomic data for Fe */
/*  Radiative recombination, Arnaud & Raymond 1992, ApJ 398, 394     */
double fe_rr_A  [3] = { 1.42e-13, 1.02e-12, 0.00 };
double fe_rr_eta [3] = { 0.891, 0.843, 0.00 };
/*  Dielectronic recombination, Arnaud & Raymond 1992, ApJ 398, 394  */
double fe_dr_e1 [3] = { 1.670, 2.860, 0.000};
double fe_dr_e2 [3] = { 31.40, 52.1, 0.00};
double fe_dr_e3 [3] = { 0.00, 0.00, 0.00};
double fe_dr_e4 [3] = { 0.00, 0.00, 0.00};
double fe_dr_c1 [3] = { 2.30e-3, 1.50e-2, 0.00};
double fe_dr_c2 [3] = { 2.7e-3, 4.7e-3, 0.00};
double fe_dr_c3 [3] = { 0.00, 0.00, 0.00};
double fe_dr_c4 [3] = { 0.00, 0.00, 0.00};
/* Charge transfer with H    - above   */
/* Charge transfer with H+   - above   */
/* Charge transfer with He   - not relevant */

  print("> MINeq: creating ionization coefficients tables...\n");
  print("  * collisional ionization\n"); 

  if (coll_ion == NULL) coll_ion = ARRAY_2D(I_g_stepNumber, NIONS,double);

  for (i = 0; i < I_g_stepNumber; i++ ) {
    T = I_TBEG + (float)i*I_TSTEP;
    for (nv = 0; nv < NIONS; nv++) {
      if ( coll_ion_Tmin[nv] < to_ev(T) ) 
        coll_ion[i][nv] = coll_ion_A[nv] * ( 1. + coll_ion_P[nv] * sqrt(coll_ion_dE[nv]/to_ev(T)) ) 
                                         / ( coll_ion_X[nv] + ( coll_ion_dE[nv]/to_ev(T) ) ) 
                                         * pow( coll_ion_dE[nv]/to_ev(T), coll_ion_K[nv] ) 
                                         * exp ( -coll_ion_dE[nv]/to_ev(T) );
      else coll_ion[i][nv] = coll_ion_A[nv] * ( 1. + coll_ion_P[nv] * sqrt(coll_ion_dE[nv]/to_ev(T)) ) 
                                         / ( coll_ion_X[nv] + ( coll_ion_dE[nv]/to_ev(T) ) ) 
                                         * pow( coll_ion_dE[nv]/to_ev(T), coll_ion_K[nv] ) 
                                         * exp ( -coll_ion_dE[nv]/to_ev(T) ) * exp( -(coll_ion_Tmin[nv]-to_ev(T))*1.e+2 );
    }
  }

  print("  * radiative recombination\n"); 
  if (rad_rec == NULL) rad_rec = ARRAY_2D(I_g_stepNumber, NIONS,double);
  
  for (i = 0; i < I_g_stepNumber; i++ ) {
    T = I_TBEG + i*I_TSTEP;
    for (nv = 0; nv < NIONS - Fe_IONS; nv++) {
      if ( (T*1.e-4/(rad_rec_z[nv]*rad_rec_z[nv])>0.004) && (T*1.e-4/(rad_rec_z[nv]*rad_rec_z[nv])<4.0) ) 
          rad_rec[i][nv] = 1.e-13 * rad_rec_z[nv] * rad_rec_a[nv] * pow( 1.e-4*T/pow(rad_rec_z[nv],2.), rad_rec_b[nv] ) 
                                  / ( 1. + rad_rec_c[nv] * pow( 1.e-4*T/pow(rad_rec_z[nv],2.), rad_rec_d[nv] ) ) ;
      else {
        if (T*1.e-4/(rad_rec_z[nv]*rad_rec_z[nv])>=4.0 && (T*1.e-4/(rad_rec_z[nv]*rad_rec_z[nv])<5.0)) {
          rad_rec[i][nv] = 1.e-13 * rad_rec_z[nv] * rad_rec_a[nv] * pow( 1.e-4*T/pow(rad_rec_z[nv],2.), rad_rec_b[nv] ) 
                                  / ( 1. + rad_rec_c[nv] * pow( 1.e-4*T/pow(rad_rec_z[nv],2.), rad_rec_d[nv] ) ) 
                                  * exp((4.0 - T*1.e-4/(rad_rec_z[nv]*rad_rec_z[nv]))*5.0);
        }
        else rad_rec[i][nv] = 0.0;
      }
    }
    lam = 157890. / T;
    rad_rec[i][0] = 5.197e-14 * sqrt(lam) * (0.4288 + 0.5 * log(lam) + 0.469 * pow(lam,-1./3.) ); 
#if Fe_IONS > 0
    for (nv = NIONS-Fe_IONS; nv < NIONS; nv++)  /* data for Fe ions */
      rad_rec[i][nv] = fe_rr_A[nv - NIONS + 3] * pow( (T*1.e-4), -fe_rr_eta[nv-NIONS+3]);
#endif
  }

  print("  * dielectronic recombination\n"); 
  if (diel_rec == NULL) diel_rec = ARRAY_2D(I_g_stepNumber,NIONS,double);
  
  for (i = 0; i < I_g_stepNumber; i++ ) {
    T = I_TBEG + i*I_TSTEP;
    t4 = T/10000.0;
    for (nv = 0; nv < NIONS - Fe_IONS; nv++) {
      diel_rec[i][nv] = 1.e-12 * ( diel_rec_a[nv] / t4 + diel_rec_b[nv] + diel_rec_c[nv]*t4 + diel_rec_d[nv] * pow(t4,2.) ) 
                               * pow( t4, -3./2. ) * exp( -diel_rec_f[nv]/t4 );
    }
#if Fe_IONS > 0
    for (nv = NIONS-Fe_IONS; nv < NIONS; nv++) {
      diel_rec[i][nv] = 1.6e-12 * pow(T, -1.5) * (  fe_dr_c1[nv-NIONS+3]*exp(-fe_dr_e1[nv-NIONS+3]/to_ev(T)) 
                                                  + fe_dr_c2[nv-NIONS+3]*exp(-fe_dr_e2[nv-NIONS+3]/to_ev(T)) 
                                                  + fe_dr_c3[nv-NIONS+3]*exp(-fe_dr_e3[nv-NIONS+3]/to_ev(T)) 
                                                  + fe_dr_c4[nv-NIONS+3]*exp(-fe_dr_e4[nv-NIONS+3]/to_ev(T)) );
    }
#endif
  }

  print("  * charge transfer with H+\n");
  if (chtr_hp == NULL) chtr_hp = ARRAY_2D(I_g_stepNumber,NIONS,double);
  
  for (i = 0; i < I_g_stepNumber; i++ ) {
    T = I_TBEG + i*I_TSTEP;
    t4 = T/10000.0;
    for (nv = 0; nv < NIONS; nv++) {
      if (t4<=chtrH_ion_upT[nv])
        chtr_hp[i][nv] = 1.e-9 * chtrH_ion_a[nv] * pow( t4, chtrH_ion_b[nv] ) 
                               * ( 1. +  chtrH_ion_c[nv] * exp ( chtrH_ion_d[nv]*t4 ));
      else chtr_hp[i][nv] = 1.e-9 * chtrH_ion_a[nv] * pow( t4, chtrH_ion_b[nv] ) 
                               * ( 1. +  chtrH_ion_c[nv] * exp ( chtrH_ion_d[nv]*t4 ))
                               * exp( -(t4-chtrH_ion_upT[nv])/(5.e-2*chtrH_ion_upT[nv])); /*  smooth to 0  */
    }
  }

  print("  * charge transfer with H\n"); 
  if (chtr_h == NULL) chtr_h = ARRAY_2D(I_g_stepNumber,NIONS,double);
  
  for (i = 0; i < I_g_stepNumber; i++ ) {
    T = I_TBEG + i*I_TSTEP;
    t4 = T/10000.0;
    for (nv = 0; nv < NIONS; nv++) {
      if (T<=chtrH_rec_upT[nv]) 
        chtr_h[i][nv] = 1.e-9 * chtrH_rec_a[nv] * pow( t4, chtrH_rec_b[nv] ) 
                              * ( 1. +  chtrH_rec_c[nv] * exp ( chtrH_rec_d[nv]*t4 ));
      else chtr_h[i][nv] = 1.e-9 * chtrH_rec_a[nv] * pow( t4, chtrH_rec_b[nv] ) 
                              * ( 1. +  chtrH_rec_c[nv] * exp ( chtrH_rec_d[nv]*t4 ))
                              * exp( -(T - chtrH_rec_upT[nv])/(5.e-2*chtrH_rec_upT[nv])); /*  smooth to 0  */
    }
  }

  print("  * charge transfer with He\n"); 
  if (chtr_he == NULL) chtr_he = ARRAY_2D(I_g_stepNumber,NIONS,double);
  
  for (i = 0; i < I_g_stepNumber; i++ ) {
    T = I_TBEG + i*I_TSTEP;
    t4 = T/10000.0;
    for (nv = 0; nv < NIONS; nv++) {
/*
      if (T<5000.) 
        chtr_he[i][nv] = elem_ab[el_He] * 1.e-9 * chtrHe_rec_a1[nv] * pow( t4, chtrHe_rec_b1[nv] ) 
                                    * ( 1. +  chtrHe_rec_c1[nv] * exp ( chtrHe_rec_d1[nv]*t4 ));
      if ( (T>=5000.) && (T<10000.) )
        chtr_he[i][nv] = elem_ab[el_He] * 1.e-9 * chtrHe_rec_a2[nv] * pow( t4, chtrHe_rec_b2[nv] ) 
                                    * ( 1. +  chtrHe_rec_c2[nv] * exp ( chtrHe_rec_d2[nv]*t4 ));
      if ( (T>=10000.) && (T<50000.) )
        chtr_he[i][nv] = elem_ab[el_He] * 1.e-9 * chtrHe_rec_a3[nv] * pow( t4, chtrHe_rec_b3[nv] ) 
                                    * ( 1. +  chtrHe_rec_c3[nv] * exp ( chtrHe_rec_d3[nv]*t4 ));
      if ( (T>=50000.) && (T<100000.) )
        chtr_he[i][nv] = elem_ab[el_He] * 1.e-9 * chtrHe_rec_a4[nv] * pow( t4, chtrHe_rec_b4[nv] ) 
                                    * ( 1. +  chtrHe_rec_c4[nv] * exp ( chtrHe_rec_d4[nv]*t4 ));
      if (T>=100000.)
        chtr_he[i][nv] = elem_ab[el_He] * 1.e-9 * chtrHe_rec_a5[nv] * pow( t4, chtrHe_rec_b5[nv] ) 
                                    * ( 1. +  chtrHe_rec_c5[nv] * exp ( chtrHe_rec_d5[nv]*t4 ));
*/
/* Replaced temperature intervals with the smooth formula below  */
      chtr_he[i][nv] = elem_ab[el_He] * 1.e-9 * chtrHe_rec_a1[nv] * pow( t4, chtrHe_rec_b1[nv] ) 
                                    * ( 1. +  chtrHe_rec_c1[nv] * exp ( chtrHe_rec_d1[nv]*t4 )) 
                                    * (T>5000.0?exp( -(T-5000.0)/1.e2):1.0)
                     + elem_ab[el_He] * 1.e-9 * chtrHe_rec_a2[nv] * pow( t4, chtrHe_rec_b2[nv] ) 
                                    * ( 1. +  chtrHe_rec_c2[nv] * exp ( chtrHe_rec_d2[nv]*t4 )) 
                                    * (T>10000.0?exp( -(T-10000.0)/2.e2):1.0) * (T<5000.0?exp( -(5000.0-T)/1.e2):1.0)
                     + elem_ab[el_He] * 1.e-9 * chtrHe_rec_a3[nv] * pow( t4, chtrHe_rec_b3[nv] ) 
                                    * ( 1. +  chtrHe_rec_c3[nv] * exp ( chtrHe_rec_d3[nv]*t4 ))
                                    * (T>50000.0?exp( -(T-50000.0)/1.e3):1.0) * (T<10000.0?exp( -(10000.0-T)/2.e2):1.0)
                     + elem_ab[el_He] * 1.e-9 * chtrHe_rec_a4[nv] * pow( t4, chtrHe_rec_b4[nv] ) 
                                    * ( 1. +  chtrHe_rec_c4[nv] * exp ( chtrHe_rec_d4[nv]*t4 ))
                                    * (T>100000.0?exp( -(T-100000.0)/2.e3):1.0) * (T<50000.0?exp( -(50000.0-T)/1.e3):1.0)
                     + elem_ab[el_He] * 1.e-9 * chtrHe_rec_a5[nv] * pow( t4, chtrHe_rec_b5[nv] ) 
                                    * ( 1. +  chtrHe_rec_c5[nv] * exp ( chtrHe_rec_d5[nv]*t4 ))
                                    * (T<100000.0?exp( -(100000.0-T)/2.e3):1.0);
    }
  }
  
  /* And now, for ions for which we only have the total electron-ion recombination coefficient:   */

  print("  * total recombination coefficients\n");
  if (tot_recs == NULL) tot_recs = ARRAY_2D(I_g_stepNumber,NIONS,double);
  for (i = 0; i < I_g_stepNumber; i++ ) 
    for (nv = 0; nv < NIONS; nv++ ) tot_recs[i][nv] = 0.0;

  for (i = 0; i < I_g_stepNumber; i++ ) {
    T = (double) (I_TBEG + i*I_TSTEP);
    if (T <= tot_rec_T[0]) { ti1=-1; ti2=0; }   /*  Find temperature vector indexes to interpolate -  ti1 and ti2 */
    if (T >  tot_rec_T[4]) { ti1=4;  ti2=-1; }
    for (j = 0; j < 4; j++) {
      if ( (tot_rec_T[j] < T) && (tot_rec_T[j+1]>=T) ) { ti1=j; ti2=j+1; }
    }

  /* -- H/He recombination -- */

    for (nv = 0; nv <= 1; nv++){
      if      (ti1 == -1) tmprec = tot_recH_He[nv][0];
      else if (ti2 == -1) tmprec = tot_recH_He[nv][4];
      else tmprec =  (T-tot_rec_T[ti1])/(tot_rec_T[ti2]-tot_rec_T[ti1])*tot_recH_He[nv][ti2] 
                   + (tot_rec_T[ti2]-T)/(tot_rec_T[ti2]-tot_rec_T[ti1])*tot_recH_He[nv][ti1];
      tot_recs[i][nv] = pow(10.0,tmprec);  
    }

  /* -- Ne recombination -- */

    for (nv = 0; nv < Ne_IONS-1; nv++){
      if      (ti1 == -1) tmprec = tot_recNe[nv][0];
      else if (ti2 == -1) tmprec = tot_recNe[nv][4];
      else tmprec =  (T-tot_rec_T[ti1])/(tot_rec_T[ti2]-tot_rec_T[ti1])*tot_recNe[nv][ti2] 
                   + (tot_rec_T[ti2]-T)/(tot_rec_T[ti2]-tot_rec_T[ti1])*tot_recNe[nv][ti1];
      tot_recs[i][X_NeI+nv-NFLX] = pow(10.0,tmprec);  
    }

  /* -- S recombination -- */

    for (nv = 0; nv < S_IONS-1; nv++){
      if      (ti1 == -1) tmprec = tot_recS[nv][0];
      else if (ti2 == -1) tmprec = tot_recS[nv][4];
      else tmprec =  (T-tot_rec_T[ti1])/(tot_rec_T[ti2]-tot_rec_T[ti1])*tot_recS[nv][ti2] 
                   + (tot_rec_T[ti2]-T)/(tot_rec_T[ti2]-tot_rec_T[ti1])*tot_recS[nv][ti1];
      tot_recs[i][X_SI+nv-NFLX] = pow(10.0,tmprec);  
    }
continue;
  
    for (nv = 0; nv < NIONS; nv++) {
      cindex=-1;
      for (j = 0; j < 10; j++) 
        if (tot_rec_ions[j] == nv) cindex = j;  /*  Find the ion index in the total rec. array   */

      if (cindex >= 0) {
        tmprec=0.;   /*    make a linear interpolation...   */
        if      (ti1 == -1) tmprec = tot_rec[cindex][0];
        else if (ti2 == -1) tmprec = tot_rec[cindex][4];
        else tmprec =  (T-tot_rec_T[ti1])/(tot_rec_T[ti2]-tot_rec_T[ti1])*tot_rec[cindex][ti2] 
                     + (tot_rec_T[ti2]-T)/(tot_rec_T[ti2]-tot_rec_T[ti1])*tot_rec[cindex][ti1];
        tot_recs[i][nv] = pow(10.0,tmprec);   /*   and obtain the rate * N(x+1)   */
/*
        if (tmprec) tot_recs[i][nv] = pow(10.0,tmprec);  
        else        tot_recs[i][nv] = 0.0;
*/
      }
    }
  }

  /* -- generate table -- */
  
  for (i = 0; i < I_g_stepNumber; i++){
    for (nv = 0; nv < NIONS; nv++) {
      tbl[i][0][nv] = coll_ion[i][nv];
      tbl[i][1][nv] =  rad_rec[i][nv];
      tbl[i][2][nv] = diel_rec[i][nv];
      tbl[i][3][nv] =  chtr_hp[i][nv];
      tbl[i][4][nv] =   chtr_h[i][nv];
      tbl[i][5][nv] =  chtr_he[i][nv];
      tbl[i][6][nv] = tot_recs[i][nv];
    }

    for (j = 0; j < 7; j++) {
      tbl[i][j][X_CI  +  C_IONS-1-NFLX] = 0.0;
      tbl[i][j][X_NI  +  N_IONS-1-NFLX] = 0.0;
      tbl[i][j][X_OI  +  O_IONS-1-NFLX] = 0.0;
      tbl[i][j][X_NeI + Ne_IONS-1-NFLX] = 0.0;
      tbl[i][j][X_SI  +  S_IONS-1-NFLX] = 0.0;
    }      
  }
  print("  * Done.\n");
  return 0;
}
#undef ZERO_5X

/* ********************************************************************* */
int Create_Losses_Tables(double ***losstables, int *nT, int *nN)
/*!
 *  Compute and save to disk the radiative 
 *  losses tables function of Ne and T
 *
 *********************************************************************** */
{
  int    ib, ie, first_time, atom_id, i, j;
  double Ne, T, erg;
  double tmpN, tmpT;
  double line_1, line_2, d;
  char fname[100],fnn[10];
  Ion atoms[NIONS], *X;   
  
  erg = 1.602177e-12;

  Ne = C_NeMIN;
  *nN = 0;
  while (Ne < C_NeMAX) {
    T   = C_TMIN;
    *nT = 0;
    while (T  < C_TMAX)  {
      tmpT = T;
      T    = T*exp(C_TSTEP);   /* should be *exp(0.02)  */
      *nT  = *nT + 1;
    }
    tmpN = Ne;
    Ne   = Ne*exp(C_NeSTEP);   /* should be *exp(0.06)  */
    *nN  = *nN + 1;
  }

  for (i = 0; i < NIONS; i++) atoms[i].dE = NULL;

  first_time = 1; 
  for (atom_id = 0; atom_id < NIONS; atom_id++) {
    X = &atoms[atom_id];
    INIT_ATOM(X,atom_id);
    Ne = C_NeMIN;
    *nN = 0;
    while (Ne < C_NeMAX) {
      T  = C_TMIN;
      *nT = 0;
      while (T < C_TMAX)  {
      
        line_2 = 0.0;
        Solve_System (X, Ne, T);
        for (ib = 1; ib < X->nlev; ib++) { 
        for (ie = 0; ie < ib; ie++) {
          if ( (X->Ni[ib] != X->Ni[ib]) || (X->A[ib][ie] != X->A[ib][ie]) || 
               (X->dE[ib][ie] != X->dE[ib][ie]) || 
               (X->Ni[ib] * X->A[ib][ie] * X->dE[ib][ie] > 1.0e+200) ){
            printf (" Nan found  atom = %d   ib,ie = %d, %d\n", atom_id,ib,ie);
            printf (" T, ne = %12.6e  %12.6e   %12.6e   %12.6e   %12.6e  %12.6e\n",T, Ne,X->Ni[ib],X->A[ib][ie],X->dE[ib][ie], X->omega[ib][ie][0] );
            QUIT_PLUTO(1);
          }      

          if (X->A[ib][ie] > 0.0 && X->Ni[ib] > 0.0) 
            line_2 += X->Ni[ib] * X->A[ib][ie] * X->dE[ib][ie] / Ne ;
                  /* cooling / ion in eV, should be multiplied by ion abundance and Ne to obtain energy/cmc/s*/
        }}
        line_2 = line_2*erg;

        losstables[atom_id][*nN][*nT] = line_2;
        tmpT = T;
        T    = T*exp(C_TSTEP);   /* should be *exp(0.02)  */
        *nT  = *nT + 1;
      }
      first_time = 0;
      tmpN = Ne;
      Ne   = Ne*exp(C_NeSTEP);   /* should be *exp(0.06)  */
      *nN  = *nN + 1;
    }
  }
  
  print("> MINEq radiative losses tables generated and saved to memory.\n");
  return(0);
}
