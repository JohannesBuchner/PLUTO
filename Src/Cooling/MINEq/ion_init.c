#include "pluto.h"
#include "cooling_defs.h"

/* ---------------------------------------------------
    Initialize energy and transition 
    probabilities between levels i,j
  
    Note: matrix are assumed to be symmetric
          ONLY for definition purposes.

          you can initialize indifferently 
 
           A(i,j) or A(j,i) to the same number.
           
          Do not initialize both.
          Elements not defined must be 
          initialized to a negative number.
                                                
     The correct symmetrization is taken care
     of in the Symmetrize_Coeff routine.  
   --------------------------------------------------- */
    

/* ********************************************************************* */
void INIT_ATOM(Ion *at, int nat)
/*
 *  
 *********************************************************************** */
{
  int n=0;
  void (*X_INIT[])(Ion *at) = {HI_INIT, HeI_INIT, HeII_INIT
                               C_EXPAND( CI_INIT,  CII_INIT,  CIII_INIT,  CIV_INIT, CV_INIT)
                               N_EXPAND( NI_INIT,  NII_INIT,  NIII_INIT,  NIV_INIT, NV_INIT)
                               O_EXPAND( OI_INIT,  OII_INIT,  OIII_INIT,  OIV_INIT, OV_INIT)
                              Ne_EXPAND(NeI_INIT, NeII_INIT, NeIII_INIT, NeIV_INIT, NeV_INIT)
                               S_EXPAND( SI_INIT,  SII_INIT,  SIII_INIT,  SIV_INIT, SV_INIT)
                              Fe_EXPAND(FeI_INIT, FeII_INIT, FeIII_INIT)};
  at->isMAP  = 0;
  at->isCV   = 0;
  at->isH    = 0;
  at->isCHEB = 0;
 
  (X_INIT[nat])(at);
  
/*
  switch (nat) {
    case HI:    HI_INIT(at); break;
    case HeI:   HeI_INIT(at); break;
    case HeII:  HeII_INIT(at); break;
    case CI:    CI_INIT(at); break;
    case CII:   CII_INIT(at); break;
    case CIII:  CIII_INIT(at); break;
    case CIV:   CIV_INIT(at); break;
    case CV:    CV_INIT(at); break;
    case NI:  NI_INIT(at); break;
    case NII:  NII_INIT(at); break;
    case NIII: NIII_INIT(at); break;
    case NIV: NIV_INIT(at); break;
    case NV: NV_INIT(at); break;
    case OI: OI_INIT(at); break;
    case OII: OII_INIT(at); break;
    case OIII: OIII_INIT(at); break;
    case OIV: OIV_INIT(at); break;
    case OV: OV_INIT(at); break;
    case NeI: NeI_INIT(at); break;  
    case NeII: NeII_INIT(at); break;
    case NeIII: NeIII_INIT(at); break;
    case NeIV: NeIV_INIT(at); break;
    case NeV: NeV_INIT(at); break;
    case 23: SI_INIT(at); break;
    case 24: SII_INIT(at); break;
    case 25: SIII_INIT(at); break;
    case 26: SIV_INIT(at); break;
    case 27: SV_INIT(at); break;
    case 28: FeI_INIT(at); break;
    case 29: FeII_INIT(at); break;
    case 30: FeIII_INIT(at); break;
  }
*/
}

/* ********************************************************************* */
void HI_INIT (Ion *HIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int S1, S2, P2, S3, P3;  /* new data from Giovanardi et al, 1987 A&AS  */
/*  int S1, S2, P;  */  /* Commented is the old data (3 levels Pradhan)  */ 
                                                                                                    
  HIv->N = 1.0 ;    /* abundance with respect to H (Raymond) */  
  HIv->nlev = 5;
  
  HIv->isH = 1;
  
/*  HIv->nlev = 3;
 
  S1 = 0;
  P = 1;
  S2 = 2;
 
  HIv->nTom  = 4;
  HIv->Tom[0] = 0.5e4;
  HIv->Tom[1] = 1.0e4;
  HIv->Tom[2] = 1.5e4;
  HIv->Tom[3] = 2.0e4;  */

 
  S1 = 0;
  S2 = 1;
  P2 = 2;
  P3 = 3;
  S3 = 4;
 
  HIv->nTom  = 8;
  HIv->Tom[0] = 0.5e4;
  HIv->Tom[1] = 1.0e4;
  HIv->Tom[2] = 1.5e4;
  HIv->Tom[3] = 2.0e4;  
  HIv->Tom[4] = 2.5e4;
  HIv->Tom[5] = 3.0e4;
  HIv->Tom[6] = 3.5e4;
  HIv->Tom[7] = 4.0e4;  
  

  if (HIv->dE == NULL){
    HIv->dE    = ARRAY_2D(HIv->nlev, HIv->nlev, double);
    HIv->A     = ARRAY_2D(HIv->nlev, HIv->nlev, double);
    HIv->omega = ARRAY_3D(HIv->nlev, HIv->nlev, HIv->nTom, double);
  }
  

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < HIv->nlev ; i++) {
  for (j = 0 ; j < HIv->nlev ; j++) {
    HIv->dE[i][j] = -1.0;
    HIv->A[i][j]  = -1.0;
    for (n = 0; n < HIv->nTom; n++){
      HIv->omega[i][j][n]  = -1.0;
    }
  }}

   
  HIv->wght[S1] = 2.; /* 2 J + 1 */
  HIv->wght[S2] = 2.; /*         */
  HIv->wght[P2] = 4.;  /*         */
  HIv->wght[P3] = 2.;  /*         */
  HIv->wght[S3] = 2.;  /*         */


/*  HIv->omega[S1][S2][0] = 2.55e-1;
  HIv->omega[S1][S2][1] = 2.74e-1;
  HIv->omega[S1][S2][2] = 2.81e-1;
  HIv->omega[S1][S2][3] = 2.84e-1;
  
  HIv->omega[S1][P][0] = 4.16e-1;
  HIv->omega[S1][P][1] = 4.72e-1;
  HIv->omega[S1][P][2] = 5.28e-1;
  HIv->omega[S1][P][3] = 5.85e-1;

  HIv->omega[S2][P][0] = 4.16e-30;
  HIv->omega[S2][P][1] = 4.16e-30;
  HIv->omega[S2][P][2] = 4.16e-30;
  HIv->omega[S2][P][3] = 4.16e-30; */
  
  HIv->omega[S2][S1][0] = 2.297e-1;
  HIv->omega[S2][S1][1] = 5.318e-6;
  HIv->omega[S2][S1][2] = -1.180e-10;
  HIv->omega[S2][S1][3] = 8.636e-16;
  HIv->omega[S2][S1][4] = 2.694e-1;
  HIv->omega[S2][S1][5] = 7.883e-7;
  HIv->omega[S2][S1][6] = -1.394e-12;
  HIv->omega[S2][S1][7] = 1.451e-18;
  
  HIv->omega[P2][S1][0] = 3.435e-1;
  HIv->omega[P2][S1][1] = 1.297e-5;
  HIv->omega[P2][S1][2] = 2.178e-12;
  HIv->omega[P2][S1][3] = 7.928e-17;
  HIv->omega[P2][S1][4] = 3.162e-1;
  HIv->omega[P2][S1][5] = 1.472e-5;
  HIv->omega[P2][S1][6] = -8.275e-12;
  HIv->omega[P2][S1][7] = -8.794e-19;
/*
  HIv->omega[P2][S2][0] = 0.000000;
  HIv->omega[P2][S2][1] = 0.000000;
  HIv->omega[P2][S2][2] = 0.000000;
  HIv->omega[P2][S2][3] = 0.000000;
  HIv->omega[P2][S2][4] = 0.000000;
  HIv->omega[P2][S2][5] = 0.000000;
  HIv->omega[P2][S2][6] = 0.000000;
  HIv->omega[P2][S2][7] = 0.000000;
*/
/*  HIv->omega[S3][S1][0] = 6.250e-2;
  HIv->omega[S3][S1][1] = -1.299e-6;
  HIv->omega[S3][S1][2] = 2.666e-11;
  HIv->omega[S3][S1][3] = -1.596e-16;
  HIv->omega[S3][S1][4] = 3.337e-2;
  HIv->omega[S3][S1][5] = 2.223e-7;
  HIv->omega[S3][S1][6] = -2.794e-13;
  HIv->omega[S3][S1][7] = 1.516e-19; */

  HIv->omega[P3][S1][0] = 9.941e-2;
  HIv->omega[P3][S1][1] = -3.714e-7;
  HIv->omega[P3][S1][2] = 6.134e-11;
  HIv->omega[P3][S1][3] = -3.973e-16;
  HIv->omega[P3][S1][4] = 6.985e-2;
  HIv->omega[P3][S1][5] = 2.538e-6;
  HIv->omega[P3][S1][6] = -8.729e-13;
  HIv->omega[P3][S1][7] = -1.291e-18;
/*
  HIv->omega[S3][S2][0] = 0.00000;
  HIv->omega[S3][S2][1] = 0.00000;
  HIv->omega[S3][S2][2] = 0.00000;
  HIv->omega[S3][S2][3] = 0.00000;
  HIv->omega[S3][S2][4] = 0.00000;
  HIv->omega[S3][S2][5] = 0.00000;
  HIv->omega[S3][S2][6] = 0.00000;
  HIv->omega[S3][S2][7] = 0.00000;  
*/
/*  HIv->omega[S3][S2][0] = 1.326;
  HIv->omega[S3][S2][1] = -1.727e-5;
  HIv->omega[S3][S2][2] = 8.914e-10;
  HIv->omega[S3][S2][3] = -6.101e-15;
  HIv->omega[S3][S2][4] = 6.311e-1;
  HIv->omega[S3][S2][5] = 2.881e-5;
  HIv->omega[S3][S2][6] = -5.372e-11;
  HIv->omega[S3][S2][7] = 4.095e-17;  */
 
  HIv->omega[P3][S2][0] = 2.040;
  HIv->omega[P3][S2][1] = -1.580e-5;
  HIv->omega[P3][S2][2] = 1.908e-9;
  HIv->omega[P3][S2][3] = -1.027e-14;
  HIv->omega[P3][S2][4] = -1.334;
  HIv->omega[P3][S2][5] = 1.229e-4;
  HIv->omega[P3][S2][6] = -9.676e-11;
  HIv->omega[P3][S2][7] = 2.842e-17; 
 
  HIv->omega[S3][P2][0] = 1.690;
  HIv->omega[S3][P2][1] = 4.563e-5;
  HIv->omega[S3][P2][2] = -6.605e-10;
  HIv->omega[S3][P2][3] = 4.445e-15;
  HIv->omega[S3][P2][4] = 2.325;
  HIv->omega[S3][P2][5] = 1.361e-5;
  HIv->omega[S3][P2][6] = -2.702e-11;
  HIv->omega[S3][P2][7] = 2.397e-17; 

/*  HIv->omega[P3][P2][0] = 4.923;
  HIv->omega[P3][P2][1] = 1.525e-4;
  HIv->omega[P3][P3][2] = 4.370e-11;
  HIv->omega[P3][P2][3] = -3.914e-15;
  HIv->omega[P3][P2][4] = 6.984;
  HIv->omega[P3][P2][5] = 1.260e-4;
  HIv->omega[P3][P2][6] = -3.014e-10;
  HIv->omega[P3][P2][7] = 2.655e-16;  */

  HIv->omega[P3][P2][0] = 0.000000;
  HIv->omega[P3][P2][1] = 0.000000;
  HIv->omega[P3][P3][2] = 0.000000;
  HIv->omega[P3][P2][3] = 0.000000;
  HIv->omega[P3][P2][4] = 0.000000;
  HIv->omega[P3][P2][5] = 0.000000;
  HIv->omega[P3][P2][6] = 0.000000;
  HIv->omega[P3][P2][7] = 0.000000;  


  HIv->dE[S2][S1]  = 1215.6740;  HIv->A[S2][S1] = 6.265e+8;  /* dE [ HIGHER ][LOWER ]  */
  HIv->dE[P2][S1]  = 1215.6680;  HIv->A[P2][S1]  = 6.265e+8;
/*  HI->dE[P2][S2]  = 2.5e+16;    HI->A[P2][S2]  = 0.000000;
  HI->dE[S3][S1]  = 2.5e+16;    HI->A[S3][S1]  = 0.000000;
*/
  HIv->dE[P3][S1]  = 1025.72;    HIv->A[P3][S1]  = 1.672e+8;
/*  HI->dE[S3][S2]  = 2.5e+16;    HI->A[S3][S2]  = 0.000000;
*/
  HIv->dE[P3][S2]  = 6562.77;    HIv->A[P3][S2]  = 2.245e+7;
  HIv->dE[S3][P2]  = 6562.91;    HIv->A[S3][P2]  = 4.209e+6;
/*  HIv->dE[P3][P2]  = 2.5e+16;    HIv->A[P3][P2]  = 0.000000;
  HIv->dE[S3][P3]  = 2.5e+16;    HIv->A[S3][P3]  = 0.000000;
*/                                                                                                                
  Symmetrize_Coeff (HIv);
}

/* ********************************************************************* */
void HeI_INIT (Ion *HeIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int S11, S23, S21, P23, P21;
                                                                                                    
  HeIv->N = 0.1;    /* abundance of He with respect to H (Raymond) */  
  HeIv->nlev = 5;
 
  S11 = 0;
  S23 = 1;
  S21 = 2;
  P23 = 3;
  P21 = 4;
 
  HeIv->nTom  = 4;
  HeIv->Tom[0] = 0.5e4;
  HeIv->Tom[1] = 1.0e4;
  HeIv->Tom[2] = 1.5e4;
  HeIv->Tom[3] = 2.0e4;

  if (HeIv->dE == NULL){
    HeIv->dE    = ARRAY_2D(HeIv->nlev, HeIv->nlev, double);
    HeIv->A     = ARRAY_2D(HeIv->nlev, HeIv->nlev, double);
    HeIv->omega = ARRAY_3D(HeIv->nlev, HeIv->nlev, HeIv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < HeIv->nlev ; i++) {
  for (j = 0 ; j < HeIv->nlev ; j++) {
    HeIv->dE[i][j] = -1.0;
    HeIv->A[i][j]  = -1.0;
    for (n = 0; n < HeIv->nTom; n++){
      HeIv->omega[i][j][n]  = -1.0;
    }
  }}

   
  HeIv->wght[S11] = 1.; /*    2 J  +  1     */
  HeIv->wght[S23] = 3.; /*         */
  HeIv->wght[S21] = 1.; /*         */
  HeIv->wght[P21] = 3.; /*         */
  HeIv->wght[P23] = 1.; /*        */


  HeIv->omega[S11][S23][0] = 6.50e-2;
  HeIv->omega[S11][S23][1] = 6.87e-2;
  HeIv->omega[S11][S23][2] = 6.81e-2;
  HeIv->omega[S11][S23][3] = 6.72e-2;
  
  HeIv->omega[S11][S21][0] = 3.11e-2;
  HeIv->omega[S11][S21][1] = 3.61e-2;
  HeIv->omega[S11][S21][2] = 3.84e-2;
  HeIv->omega[S11][S21][3] = 4.01e-2;

  HeIv->omega[S11][P23][0] = 1.60e-2;
  HeIv->omega[S11][P23][1] = 2.27e-2;
  HeIv->omega[S11][P23][2] = 2.71e-2;
  HeIv->omega[S11][P23][3] = 3.07e-2;
  
  HeIv->omega[S11][P21][0] = 9.92e-3;
  HeIv->omega[S11][P21][1] = 1.54e-2;
  HeIv->omega[S11][P21][2] = 1.98e-2;
  HeIv->omega[S11][P21][3] = 2.40e-2;

  HeIv->omega[S23][S21][0] = 2.24;
  HeIv->omega[S23][S21][1] = 2.40;
  HeIv->omega[S23][S21][2] = 2.32;
  HeIv->omega[S23][S21][3] = 2.20;

  HeIv->omega[S23][P23][0] = 1.50e+1;
  HeIv->omega[S23][P23][1] = 2.69e+1;
  HeIv->omega[S23][P23][2] = 3.74e+1;
  HeIv->omega[S23][P23][3] = 4.66e+1;
  
  HeIv->omega[S23][P21][0] = 7.70e-1;
  HeIv->omega[S23][P21][1] = 9.75e-1;
  HeIv->omega[S23][P21][2] = 1.05e+0;
  HeIv->omega[S23][P21][3] = 1.08e+0;

  HeIv->omega[S21][P23][0] = 1.50;
  HeIv->omega[S21][P23][1] = 1.70;
  HeIv->omega[S21][P23][2] = 1.74;
  HeIv->omega[S21][P23][3] = 1.72;

  HeIv->omega[S21][P21][0] = 9.73e+0;
  HeIv->omega[S21][P21][1] = 1.86e+1;
  HeIv->omega[S21][P21][2] = 2.58e+1;
  HeIv->omega[S21][P21][3] = 3.32e+1;
/*
  HeIv->omega[P23][P21][0] = 1.45e-30;
  HeIv->omega[P23][P21][1] = 1.45e-30;
  HeIv->omega[P23][P21][2] = 1.45e-30;
  HeIv->omega[P23][P21][3] = 1.45e-30;
*/
  HeIv->dE[S11][S23]   = 625.48;  HeIv->A[S11][S23]  = 1.13e-4;  /* dE [ HIGHER ][LOWER ]  */
  HeIv->dE[S11][S21]   = 601.30;  HeIv->A[S11][S21]  = 5.13e+1;
  HeIv->dE[S11][P23]   = 591.29;  HeIv->A[S11][P23]  = 1.76e+2;
  HeIv->dE[S11][P21]   = 584.21;  HeIv->A[S11][P21]  = 1.80e+9;
  HeIv->dE[S23][S21]   = 15553.7;  HeIv->A[S23][S21]  = 1.51e-7;
  HeIv->dE[S23][P23]   = 10817.0;  HeIv->A[S23][P23]  = 1.02e+7;
  HeIv->dE[S23][P21]   = 8854.5;  HeIv->A[S23][P21]  = 1.29;
  HeIv->dE[S21][P23]   = 35519.5;  HeIv->A[S21][P23]  = 2.70e-2;
  HeIv->dE[S21][P21]   = 20557.7;  HeIv->A[S21][P21]  = 1.98e+6;
/*  HeIv->dE[P23][P21]   = 1e+11;   HeIv->A[P23][P21]  = 0.0000;    
  */                                                                                                                
  Symmetrize_Coeff (HeIv);
}

/* ********************************************************************* */
void HeII_INIT (Ion *HeIIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int S1, S2, P;
                                                                                                    
  HeIIv->N = 0.1;    /* abundance with respect to H (Raymond) */  
  HeIIv->nlev = 3;
 
  S1 = 0;
  S2 = 2;
  P = 1;
 
  HeIIv->nTom  = 4;
  HeIIv->Tom[0] = 0.5e4;
  HeIIv->Tom[1] = 1.0e4;
  HeIIv->Tom[2] = 1.5e4;
  HeIIv->Tom[3] = 2.0e4;

  if (HeIIv->dE == NULL){
    HeIIv->dE    = ARRAY_2D(HeIIv->nlev, HeIIv->nlev, double);
    HeIIv->A     = ARRAY_2D(HeIIv->nlev, HeIIv->nlev, double);
    HeIIv->omega = ARRAY_3D(HeIIv->nlev, HeIIv->nlev, HeIIv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < HeIIv->nlev ; i++) {
  for (j = 0 ; j < HeIIv->nlev ; j++) {
    HeIIv->dE[i][j] = -1.0;
    HeIIv->A[i][j]  = -1.0;
    for (n = 0; n < HeIIv->nTom; n++){
      HeIIv->omega[i][j][n]  = -1.0;
    }
  }}

   
  HeIIv->wght[S1] = 2.; /*    2 J  +  1     */
  HeIIv->wght[S2] = 2.; /*         */
  HeIIv->wght[P] = 2.; /*         */


  HeIIv->omega[S1][S2][0] = 1.60e-1;
  HeIIv->omega[S1][S2][1] = 1.59e-1;
  HeIIv->omega[S1][S2][2] = 1.57e-1;
  HeIIv->omega[S1][S2][3] = 1.56e-1;
  
  HeIIv->omega[S1][P][0] = 3.40e-1;
  HeIIv->omega[S1][P][1] = 3.53e-1;
  HeIIv->omega[S1][P][2] = 3.63e-1;
  HeIIv->omega[S1][P][3] = 3.73e-1;
/*  
  HeIIv->omega[S2][P][0] = 3.40e-30;  
  HeIIv->omega[S2][P][1] = 3.40e-30;
  HeIIv->omega[S2][P][2] = 3.40e-30;
  HeIIv->omega[S2][P][3] = 3.40e-30;
*/
  HeIIv->dE[S1][S2]   = 303.92;  HeIIv->A[S1][S2]  = 5.66e+2;  /* dE [ HIGHER ][LOWER ]  */
  HeIIv->dE[S1][P]   = 303.92;  HeIIv->A[S1][P]  = 1.0e+10;
/*  HeIIv->dE[S2][P]   = 2.0e+30;  HeIIv->A[S2][P]  = 0.00000;   Fake data */
                                                                                                                  
  Symmetrize_Coeff (HeIIv);
}

/* ********************************************************************* */
void CI_INIT (Ion *CIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int _3P0, _3P1, _3P2, _1D2, _1S0, _5S2;
            
  /* Old data, commented, is from Pradhan. New T-dep CS from MAPPINGS  */
                                                                                                    
  CIv->N = 5.e-4;    /* abundance of C with respect to H (Raymond) */  
  CIv->nlev = 5;
 
  CIv->isMAP = 0;
  
  _3P0 = 0;
  _3P1 = 1;
  _3P2 = 2;
  _1D2 = 3;
  _1S0 = 4;
/*  _5S2 = 5;  */
 
  CIv->nTom  = 4;
  CIv->Tom[0] = 0.5e4;
  CIv->Tom[1] = 1.0e4;
  CIv->Tom[2] = 1.5e4;
  CIv->Tom[3] = 2.0e4;

  if (CIv->dE == NULL){
    CIv->dE    = ARRAY_2D(CIv->nlev, CIv->nlev, double);
    CIv->A     = ARRAY_2D(CIv->nlev, CIv->nlev, double);
    CIv->omega = ARRAY_3D(CIv->nlev, CIv->nlev, CIv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < CIv->nlev ; i++) {
  for (j = 0 ; j < CIv->nlev ; j++) {
    CIv->dE[i][j] = -1.0;
    CIv->A[i][j]  = -1.0;
    for (n = 0; n < CIv->nTom; n++){
      CIv->omega[i][j][n]  = -1.0;
    }
  }}

   
  CIv->wght[_3P0] = 1.; /*    2 J  +  1     */
  CIv->wght[_3P1] = 3.; /*         */
  CIv->wght[_3P2] = 5.; /*         */
  CIv->wght[_1D2] = 5.; /*         */
  CIv->wght[_1S0] = 1.; /*        */
/*  CIv->wght[_5S2] = 5.;         */


  CIv->omega[_1D2][_3P0][0] = 6.03e-1 * (1.0)/(3.0*3.0);
  CIv->omega[_1D2][_3P0][1] = 1.14 * (1.0)/(3.0*3.0);
  CIv->omega[_1D2][_3P0][2] = 1.60 * (1.0)/(3.0*3.0);
  CIv->omega[_1D2][_3P0][3] = 1.96 * (1.0)/(3.0*3.0);

  CIv->omega[_1D2][_3P1][0] = CIv->omega[_1D2][_3P0][0] * (3.0)/(3.0*3.0);
  CIv->omega[_1D2][_3P1][1] = CIv->omega[_1D2][_3P0][1] * (3.0)/(3.0*3.0);
  CIv->omega[_1D2][_3P1][2] = CIv->omega[_1D2][_3P0][2] * (3.0)/(3.0*3.0);
  CIv->omega[_1D2][_3P1][3] = CIv->omega[_1D2][_3P0][3] * (3.0)/(3.0*3.0);

  CIv->omega[_1D2][_3P2][0] = CIv->omega[_1D2][_3P0][0] * (5.0)/(3.0*3.0);
  CIv->omega[_1D2][_3P2][1] = CIv->omega[_1D2][_3P0][0] * (5.0)/(3.0*3.0);
  CIv->omega[_1D2][_3P2][2] = CIv->omega[_1D2][_3P0][0] * (5.0)/(3.0*3.0);
  CIv->omega[_1D2][_3P2][3] = CIv->omega[_1D2][_3P0][0] * (5.0)/(3.0*3.0);

  CIv->omega[_1S0][_3P1][0] = 1.49e-1 * (3.0)/(3.0*3.0);
  CIv->omega[_1S0][_3P1][1] = 2.52e-1 * (3.0)/(3.0*3.0);
  CIv->omega[_1S0][_3P1][2] = 3.20e-1 * (3.0)/(3.0*3.0);
  CIv->omega[_1S0][_3P1][3] = 3.65e-1 * (3.0)/(3.0*3.0);

  CIv->omega[_1S0][_3P2][0] = CIv->omega[_1S0][_3P1][0] * (5.0)/(3.0*3.0);
  CIv->omega[_1S0][_3P2][1] = CIv->omega[_1S0][_3P1][1] * (5.0)/(3.0*3.0);
  CIv->omega[_1S0][_3P2][2] = CIv->omega[_1S0][_3P1][2] * (5.0)/(3.0*3.0);
  CIv->omega[_1S0][_3P2][3] = CIv->omega[_1S0][_3P1][3] * (5.0)/(3.0*3.0);

  CIv->omega[_1S0][_1D2][0] = 1.96e-1;
  CIv->omega[_1S0][_1D2][1] = 2.77e-1;
  CIv->omega[_1S0][_1D2][2] = 3.40e-1;
  CIv->omega[_1S0][_1D2][3] = 3.92e-1;

  CIv->omega[_3P1][_3P0][0] = 2.43e-1;
  CIv->omega[_3P1][_3P0][1] = 3.71e-1;
  CIv->omega[_3P1][_3P0][2] = 3.71e-1; 
  CIv->omega[_3P1][_3P0][3] = 3.71e-1; 

  CIv->omega[_3P2][_3P0][0] = 1.82e-1;
  CIv->omega[_3P2][_3P0][1] = 2.46e-1;
  CIv->omega[_3P2][_3P0][2] = 2.46e-1; 
  CIv->omega[_3P2][_3P0][3] = 2.46e-1; 
 
  CIv->omega[_3P2][_3P1][0] = 7.14e-1;
  CIv->omega[_3P2][_3P1][1] = 1.02;
  CIv->omega[_3P2][_3P1][2] = 1.02; 
  CIv->omega[_3P2][_3P1][3] = 1.02; 
/*
  CIv->omega[_5S2][_3P1][0] = 4.75e-1;
  CIv->omega[_5S2][_3P1][1] = 6.71e-1;
  CIv->omega[_5S2][_3P1][2] = 8.22e-1;
  CIv->omega[_5S2][_3P1][3] = 9.50e-1;

  CIv->omega[_5S2][_3P2][0] = CIv->omega[_5S2][_3P1][0] * (5.0)/(3.0*3.0);
  CIv->omega[_5S2][_3P2][1] = CIv->omega[_5S2][_3P1][1] * (5.0)/(3.0*3.0);
  CIv->omega[_5S2][_3P2][2] = CIv->omega[_5S2][_3P1][2] * (5.0)/(3.0*3.0);
  CIv->omega[_5S2][_3P2][3] = CIv->omega[_5S2][_3P1][3] * (5.0)/(3.0*3.0);
*/
/* 
   And now some fake data to make things work 
   We will add A=0 for these transitions, this is only to initialize the matrix and avoid NaNs
*/
/*
  CIv->omega[_1S0][_3P0][0] = 1.00e-30;
  CIv->omega[_1S0][_3P0][1] = 1.00e-30;
  CIv->omega[_1S0][_3P0][2] = 1.00e-30;
  CIv->omega[_1S0][_3P0][3] = 1.00e-30;

  CIv->omega[_5S2][_1S0][0] = 1.00e-30;
  CIv->omega[_5S2][_1S0][1] = 1.00e-30;
  CIv->omega[_5S2][_1S0][2] = 1.00e-30;
  CIv->omega[_5S2][_1S0][3] = 1.00e-30;

  CIv->omega[_5S2][_3P0][0] = 1.00e-30;
  CIv->omega[_5S2][_3P0][1] = 1.00e-30;
  CIv->omega[_5S2][_3P0][2] = 1.00e-30;
  CIv->omega[_5S2][_3P0][3] = 1.00e-30;

  CIv->omega[_5S2][_1D2][0] = 1.00e-30;
  CIv->omega[_5S2][_1D2][1] = 1.00e-30;
  CIv->omega[_5S2][_1D2][2] = 1.00e-30;
  CIv->omega[_5S2][_1D2][3] = 1.00e-30;
*/

/* End fake data */

/*
  CIv->omega[_1D2][_3P0][0] = 1.267e-1;
  CIv->omega[_1D2][_3P0][1] = 0.9190;

  CIv->omega[_1D2][_3P1][0] = 3.80e-1;
  CIv->omega[_1D2][_3P1][1] = 0.9190;

  CIv->omega[_1D2][_3P2][0] = 6.333e-1;
  CIv->omega[_1D2][_3P2][1] = 0.9190;

  CIv->omega[_1S0][_3P1][0] = 8.40e-2;
  CIv->omega[_1S0][_3P1][1] = 0.7580;

  CIv->omega[_1S0][_3P2][0] = 1.40e-1;
  CIv->omega[_1S0][_3P2][1] = 0.7580;

  CIv->omega[_1S0][_1D2][0] = 2.77e-1;
  CIv->omega[_1S0][_1D2][1] = 0.4990;

  CIv->omega[_3P1][_3P0][0] = 3.71e-1;
  CIv->omega[_3P1][_3P0][1] = 0.0000; 

  CIv->omega[_3P2][_3P0][0] = 2.46e-1; 
  CIv->omega[_3P2][_3P0][1] = 0.0000;
 
  CIv->omega[_3P2][_3P1][0] = 1.00;
  CIv->omega[_3P2][_3P1][1] = 0.0000;
*/

  CIv->dE[_1D2][_3P0]   = 9808.13;    CIv->A[_1D2][_3P0]  = 8.62e-8;  
  CIv->dE[_1D2][_3P1]   = 9824.12;   CIv->A[_1D2][_3P1]  = 6.05e-5;
  CIv->dE[_1D2][_3P2]   = 9850.26;    CIv->A[_1D2][_3P2]  = 1.80e-4;  
  CIv->dE[_1S0][_3P1]   = 4621.57;    CIv->A[_1S0][_3P1]  = 2.1e-3;
  CIv->dE[_1S0][_3P2]   = 4627.35;    CIv->A[_1S0][_3P2]  = 1.79e-5;
  CIv->dE[_1S0][_1D2]   = 8727.18;    CIv->A[_1S0][_1D2]  = 6.34e-1;
  CIv->dE[_3P1][_3P0]   = 6.094e+6;   CIv->A[_3P1][_3P0]  = 7.88e-8;
  CIv->dE[_3P2][_3P0]   = 2304147.;   CIv->A[_3P2][_3P0]  = 1.81e-14;
  CIv->dE[_3P2][_3P1]   = 3704140.;   CIv->A[_3P2][_3P1]  = 2.65e-7;
/*  CIv->dE[_5S2][_3P1]   = 2965.70;    CIv->A[_5S2][_3P1]  = 6.94e+0;
  CIv->dE[_5S2][_3P2]   = 2968.08;    CIv->A[_5S2][_3P2]  = 1.56e+1;
  CIv->dE[_5S2][_1S0]   = 1.e+30;    CIv->A[_5S2][_1S0]  = 0.00000; 
  CIv->dE[_5S2][_3P0]   = 1.e+30;    CIv->A[_5S2][_3P0]  = 0.00000;
  CIv->dE[_5S2][_1D2]   = 1.e+30;    CIv->A[_5S2][_1D2]  = 0.00000;  
  CIv->dE[_1S0][_3P0]   = 1.e+30;    CIv->A[_1S0][_3P0]  = 0.00000; */


  Symmetrize_Coeff (CIv);
}

/* ********************************************************************* */
void CII_INIT (Ion *CIIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int _2P12, _2P32, _4P12, _4P32, _4P52;
            
  /* Old data from Pradhan commented. Enlarged T-dep from Blum & Pradhan 1992, ApJSS  */
                                                                                                    
  CIIv->N = 5.e-4;    /* abundance of C with respect to H (Raymond) */  
  CIIv->nlev = 5;
 
  _2P12 = 0;
  _2P32 = 1;
  _4P12 = 2;
  _4P32 = 3;
  _4P52 = 4;
 
  CIIv->nTom  = 4;
  CIIv->Tom[0] = 0.3e4;
  CIIv->Tom[1] = 1.0e4;
  CIIv->Tom[2] = 2.2e4;
  CIIv->Tom[3] = 4.0e4;

  if (CIIv->dE == NULL){
    CIIv->dE    = ARRAY_2D(CIIv->nlev, CIIv->nlev, double);
    CIIv->A     = ARRAY_2D(CIIv->nlev, CIIv->nlev, double);
    CIIv->omega = ARRAY_3D(CIIv->nlev, CIIv->nlev, CIIv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < CIIv->nlev ; i++) {
  for (j = 0 ; j < CIIv->nlev ; j++) {
    CIIv->dE[i][j] = -1.0;
    CIIv->A[i][j]  = -1.0;
    for (n = 0; n < CIIv->nTom; n++){
      CIIv->omega[i][j][n]  = -1.0;
    }
  }}

   
  CIIv->wght[_2P12] = 2.; /*    2 J  +  1     */
  CIIv->wght[_2P32] = 4.; /*         */
  CIIv->wght[_4P12] = 2.; /*         */
  CIIv->wght[_4P32] = 4.; /*         */
  CIIv->wght[_4P52] = 6.; /*        */

/*  CIIv->omega[_2P32][_2P12][0] = 1.89;
  CIIv->omega[_2P32][_2P12][1] = 2.15;
  CIIv->omega[_2P32][_2P12][2] = 2.26;
  CIIv->omega[_2P32][_2P12][3] = 2.28;

  CIIv->omega[_4P12][_2P12][0] = 2.43e-1;
  CIIv->omega[_4P12][_2P12][1] = 2.42e-1;
  CIIv->omega[_4P12][_2P12][2] = 2.46e-1;
  CIIv->omega[_4P12][_2P12][3] = 2.48e-1;

  CIIv->omega[_4P12][_2P32][0] = 1.74e-1;
  CIIv->omega[_4P12][_2P32][1] = 1.77e-1;
  CIIv->omega[_4P12][_2P32][2] = 1.82e-1;
  CIIv->omega[_4P12][_2P32][3] = 1.84e-1;

  CIIv->omega[_4P32][_2P12][0] = 3.61e-1;
  CIIv->omega[_4P32][_2P12][1] = 3.62e-1;
  CIIv->omega[_4P32][_2P12][2] = 3.68e-1;
  CIIv->omega[_4P32][_2P12][3] = 3.70e-1;

  CIIv->omega[_4P32][_2P32][0] = 4.72e-1;
  CIIv->omega[_4P32][_2P32][1] = 4.77e-1;
  CIIv->omega[_4P32][_2P32][2] = 4.88e-1;
  CIIv->omega[_4P32][_2P32][3] = 4.93e-1;

  CIIv->omega[_4P32][_4P12][0] = 6.60e-1;
  CIIv->omega[_4P32][_4P12][1] = 8.24e-1;
  CIIv->omega[_4P32][_4P12][2] = 9.64e-1;
  CIIv->omega[_4P32][_4P12][3] = 1.06;

  CIIv->omega[_4P52][_2P12][0] = 2.29e-1;
  CIIv->omega[_4P52][_2P12][1] = 2.34e-1;
  CIIv->omega[_4P52][_2P12][2] = 2.42e-1;
  CIIv->omega[_4P52][_2P12][3] = 2.45e-1;

  CIIv->omega[_4P52][_2P32][0] = 1.02;
  CIIv->omega[_4P52][_2P32][1] = 1.02;
  CIIv->omega[_4P52][_2P32][2] = 1.04;
  CIIv->omega[_4P52][_2P32][3] = 1.05;

  CIIv->omega[_4P52][_4P12][0] = 7.30e-1;
  CIIv->omega[_4P52][_4P12][1] = 8.53e-1;
  CIIv->omega[_4P52][_4P12][2] = 9.32e-1;
  CIIv->omega[_4P52][_4P12][3] = 9.71e-1;

  CIIv->omega[_4P52][_4P32][0] = 1.65;
  CIIv->omega[_4P52][_4P32][1] = 1.98;
  CIIv->omega[_4P52][_4P32][2] = 2.23;
  CIIv->omega[_4P52][_4P32][3] = 2.39;
*/

  CIIv->omega[_2P32][_2P12][0] = 1.7157;
  CIIv->omega[_2P32][_2P12][1] = 2.1519;
  CIIv->omega[_2P32][_2P12][2] = 2.2816;
  CIIv->omega[_2P32][_2P12][3] = 2.2460;

  CIIv->omega[_4P12][_2P12][0] = 2.507e-1;
  CIIv->omega[_4P12][_2P12][1] = 2.425e-1;
  CIIv->omega[_4P12][_2P12][2] = 2.476e-1;
  CIIv->omega[_4P12][_2P12][3] = 2.412e-1;

  CIIv->omega[_4P12][_2P32][0] = 1.776e-1;
  CIIv->omega[_4P12][_2P32][1] = 1.771e-1;
  CIIv->omega[_4P12][_2P32][2] = 1.845e-1;
  CIIv->omega[_4P12][_2P32][3] = 1.799e-1;

  CIIv->omega[_4P32][_2P12][0] = 3.719e-1;
  CIIv->omega[_4P32][_2P12][1] = 3.618e-1;
  CIIv->omega[_4P32][_2P12][2] = 3.701e-1;
  CIIv->omega[_4P32][_2P12][3] = 3.614e-1;

  CIIv->omega[_4P32][_2P32][0] = 4.847e-1;
  CIIv->omega[_4P32][_2P32][1] = 4.774e-1;
  CIIv->omega[_4P32][_2P32][2] = 4.934e-1;
  CIIv->omega[_4P32][_2P32][3] = 4.807e-1;

  CIIv->omega[_4P32][_4P12][0] = 6.20e-1;
  CIIv->omega[_4P32][_4P12][1] = 8.237e-1;
  CIIv->omega[_4P32][_4P12][2] = 1.0898;
  CIIv->omega[_4P32][_4P12][3] = 1.2207;

  CIIv->omega[_4P52][_2P12][0] = 2.340e-1;
  CIIv->omega[_4P52][_2P12][1] = 2.349e-1;
  CIIv->omega[_4P52][_2P12][2] = 2.457e-1;
  CIIv->omega[_4P52][_2P12][3] = 2.395e-1;

  CIIv->omega[_4P52][_2P32][0] = 1.0508;
  CIIv->omega[_4P52][_2P32][1] = 1.0238;
  CIIv->omega[_4P52][_2P32][2] = 1.0508;
  CIIv->omega[_4P52][_2P32][3] = 1.0237;

  CIIv->omega[_4P52][_4P12][0] = 6.895e-1;
  CIIv->omega[_4P52][_4P12][1] = 8.533e-1;
  CIIv->omega[_4P52][_4P12][2] = 9.789e-1;
  CIIv->omega[_4P52][_4P12][3] = 9.924e-1;

  CIIv->omega[_4P52][_4P32][0] = 1.5521;
  CIIv->omega[_4P52][_4P32][1] = 1.9818;
  CIIv->omega[_4P52][_4P32][2] = 2.4300;
  CIIv->omega[_4P52][_4P32][3] = 2.5885;

  CIIv->dE[_2P32][_2P12]   = 1.5774e+5;CIIv->A[_2P32][_2P12]  = 2.29e-6;  
  CIIv->dE[_4P12][_2P12]   = 2325.;    CIIv->A[_4P12][_2P12]  = 7.0e+1;
  CIIv->dE[_4P12][_2P32]   = 2329.;    CIIv->A[_4P12][_2P32]  = 6.3e+1;  
  CIIv->dE[_4P32][_2P12]   = 2324.;    CIIv->A[_4P32][_2P12]  = 1.4e+0;
  CIIv->dE[_4P32][_2P32]   = 2328.;    CIIv->A[_4P32][_2P32]  = 9.4e+0;
  CIIv->dE[_4P32][_4P12]   = 4.55e+6;  CIIv->A[_4P32][_4P12]  = 2.39e-7;
/*  CIIv->dE[_4P52][_2P12]   = 2323.;    CIIv->A[_4P52][_2P12]  = 0.00000;
*/
  CIIv->dE[_4P52][_2P32]   = 2326.;    CIIv->A[_4P52][_2P32]  = 5.1e+1;
  CIIv->dE[_4P52][_4P12]   = 1.99e+6;  CIIv->A[_4P52][_4P12]  = 3.49e-14;
  CIIv->dE[_4P52][_4P32]   = 3.53e+6;  CIIv->A[_4P52][_4P32]  = 3.67e-7;

  Symmetrize_Coeff (CIIv);
}


/* ********************************************************************* */
void CIII_INIT (Ion *CIIIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int _3P2, _3P1, _3P0, _1P1, _1S0;
                                                                                                    
  CIIIv->N = 5.e-4;    /* abundance of C with respect to H (Raymond) */  
  CIIIv->nlev = 5;
 
  _1S0 = 0;
  _3P0 = 1;
  _3P1 = 2;
  _3P2 = 3;
  _1P1 = 4;
 
  CIIIv->nTom  = 4;
  CIIIv->Tom[0] = 0.5e4;
  CIIIv->Tom[1] = 1.0e4;
  CIIIv->Tom[2] = 1.5e4;
  CIIIv->Tom[3] = 2.0e4;

  if (CIIIv->dE == NULL){
    CIIIv->dE    = ARRAY_2D(CIIIv->nlev, CIIIv->nlev, double);
    CIIIv->A     = ARRAY_2D(CIIIv->nlev, CIIIv->nlev, double);
    CIIIv->omega = ARRAY_3D(CIIIv->nlev, CIIIv->nlev, CIIIv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < CIIIv->nlev ; i++) {
  for (j = 0 ; j < CIIIv->nlev ; j++) {
    CIIIv->dE[i][j] = -1.0;
    CIIIv->A[i][j]  = -1.0;
    for (n = 0; n < CIIIv->nTom; n++){
      CIIIv->omega[i][j][n]  = -1.0;
    }
  }}

   
  CIIIv->wght[_1S0] = 1.; /*    2 J  +  1     */
  CIIIv->wght[_3P0] = 1.; /*         */
  CIIIv->wght[_3P1] = 3.; /*         */
  CIIIv->wght[_3P2] = 5.; /*         */
  CIIIv->wght[_1P1] = 3.; /*        */

  CIIIv->omega[_3P2][_1S0][0] = 1.12 * (5.0)/(3.0 * 3.0);
  CIIIv->omega[_3P2][_1S0][1] = 1.01 * (5.0)/(3.0 * 3.0);
  CIIIv->omega[_3P2][_1S0][2] = 9.90e-1 * (5.0)/(3.0 * 3.0);
  CIIIv->omega[_3P2][_1S0][3] = 9.96e-1 * (5.0)/(3.0 * 3.0);

  CIIIv->omega[_3P1][_1S0][0] = 1.12 * (3.0)/(3.0 * 3.0);
  CIIIv->omega[_3P1][_1S0][1] = 1.01 * (3.0)/(3.0 * 3.0);
  CIIIv->omega[_3P1][_1S0][2] = 9.90e-1 * (3.0)/(3.0 * 3.0);
  CIIIv->omega[_3P1][_1S0][3] = 9.96e-1 * (3.0)/(3.0 * 3.0);

  CIIIv->omega[_1P1][_1S0][0] = 3.85;
  CIIIv->omega[_1P1][_1S0][1] = 4.34;
  CIIIv->omega[_1P1][_1S0][2] = 4.56;
  CIIIv->omega[_1P1][_1S0][3] = 4.69;
  
  CIIIv->omega[_3P1][_3P0][0] = 8.48e-1;
  CIIIv->omega[_3P1][_3P0][1] = 9.11e-1;
  CIIIv->omega[_3P1][_3P0][2] = 9.75e-1;
  CIIIv->omega[_3P1][_3P0][3] = 1.03;
  
  CIIIv->omega[_3P2][_3P1][0] = 2.36;
  CIIIv->omega[_3P2][_3P1][1] = 2.66;
  CIIIv->omega[_3P2][_3P1][2] = 2.97;
  CIIIv->omega[_3P2][_3P1][3] = 3.23;

  CIIIv->omega[_3P2][_3P0][0] = 0.579;
  CIIIv->omega[_3P2][_3P0][1] = 0.677;
  CIIIv->omega[_3P2][_3P0][2] = 0.776;
  CIIIv->omega[_3P2][_3P0][3] = 0.867;

  /*   fake data   */
/*
  CIIIv->omega[_3P2][_3P0][0] = 8.48e-31;
  CIIIv->omega[_3P2][_3P0][1] = 9.11e-31;
  CIIIv->omega[_3P2][_3P0][2] = 9.75e-31;
  CIIIv->omega[_3P2][_3P0][3] = 1.03e-30;

  CIIIv->omega[_1P1][_3P0][0] = 8.48e-31;
  CIIIv->omega[_1P1][_3P0][1] = 9.11e-31;
  CIIIv->omega[_1P1][_3P0][2] = 9.75e-31;
  CIIIv->omega[_1P1][_3P0][3] = 1.03e-30;

  CIIIv->omega[_1P1][_3P1][0] = 8.48e-31;
  CIIIv->omega[_1P1][_3P1][1] = 9.11e-31;
  CIIIv->omega[_1P1][_3P1][2] = 9.75e-31;
  CIIIv->omega[_1P1][_3P1][3] = 1.03e-30;

  CIIIv->omega[_1P1][_3P2][0] = 8.48e-31;
  CIIIv->omega[_1P1][_3P2][1] = 9.11e-31;
  CIIIv->omega[_1P1][_3P2][2] = 9.75e-31;
  CIIIv->omega[_1P1][_3P2][3] = 1.03e-30;

  CIIIv->omega[_3P0][_1S0][0] = 8.48e-31;
  CIIIv->omega[_3P0][_1S0][1] = 9.11e-31;
  CIIIv->omega[_3P0][_1S0][2] = 9.75e-31;
  CIIIv->omega[_3P0][_1S0][3] = 1.03e-30;
*/

  CIIIv->dE[_3P2][_1S0]   = 1907.  ;  CIIIv->A[_3P2][_1S0]  = 5.19e-3;  
  CIIIv->dE[_3P1][_1S0]   = 1909.  ;  CIIIv->A[_3P1][_1S0]  = 1.14e+2;
  CIIIv->dE[_1P1][_1S0]   = 977.02;  CIIIv->A[_1P1][_1S0]  = 1.767e+9;  
  CIIIv->dE[_3P1][_3P0]   = 4.22e+6;  CIIIv->A[_3P1][_3P0]  = 1.4e+0;
  CIIIv->dE[_3P2][_3P1]   = 1.774e+6; CIIIv->A[_3P2][_3P1]  = 9.4e+0;
  CIIIv->dE[_3P2][_3P0]   = 1.25e+6 ;  CIIIv->A[_3P2][_3P0]  = 1.5e-13;
/*  CIIIv->dE[_1P1][_3P0]   = 2323.e+13;CIIIv->A[_1P1][_3P0]  = 1.e-31;
  CIIIv->dE[_1P1][_3P1]   = 2326.e+13;CIIIv->A[_1P1][_3P1]  = 1.e-31;
  CIIIv->dE[_1P1][_3P2]   = 1.99e+16; CIIIv->A[_1P1][_3P2]  = 1.e-31;
  CIIIv->dE[_3P0][_1S0]   = 3.53e+16; CIIIv->A[_3P0][_1S0]  = 1.e-31;
*/
  Symmetrize_Coeff (CIIIv);
}

/* ********************************************************************* */
void CIV_INIT (Ion *CIVv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int S12, P12, P32;
                                                                                                    
  CIVv->N = 5.e-4;    /* abundance of C with respect to H (Raymond) */  
  CIVv->nlev = 3;
 
  S12 = 0;
  P12 = 1;
  P32 = 2;
 
  CIVv->nTom  = 2;
  CIVv->Tom[0] = 1.0e4;
  CIVv->Tom[1] = 2.0e4;

  if (CIVv->dE == NULL){
    CIVv->dE    = ARRAY_2D(CIVv->nlev, CIVv->nlev, double);
    CIVv->A     = ARRAY_2D(CIVv->nlev, CIVv->nlev, double);
    CIVv->omega = ARRAY_3D(CIVv->nlev, CIVv->nlev, CIVv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < CIVv->nlev ; i++) {
  for (j = 0 ; j < CIVv->nlev ; j++) {
    CIVv->dE[i][j] = -1.0;
    CIVv->A[i][j]  = -1.0;
    for (n = 0; n < CIVv->nTom; n++){
      CIVv->omega[i][j][n]  = -1.0;
    }
  }}

   
  CIVv->wght[S12] = 2.; /*    2 J  +  1     */
  CIVv->wght[P12] = 2.; /*         */
  CIVv->wght[P32] = 4.; /*         */

  CIVv->omega[P32][S12][0] = 8.88*(4.0)/(2.0 * 3.0);
  CIVv->omega[P32][S12][1] = 8.95*(4.0)/(2.0 * 3.0);

  CIVv->omega[P12][S12][0] = CIVv->omega[P32][S12][0]*(2.0)/(2.0 * 3.0);
  CIVv->omega[P12][S12][1] = CIVv->omega[P32][S12][1]*(2.0)/(2.0 * 3.0);
/*
  CIVv->omega[P32][P12][0] = 8.88e-31;
  CIVv->omega[P32][P12][1] = 8.88e-31;
*/

  CIVv->dE[P32][S12]   = 1548.2  ;  CIVv->A[P32][S12]  = 2.65e+8;  
  CIVv->dE[P12][S12]   = 1550.8  ;  CIVv->A[P12][S12]  = 2.63e+8;
/*  CIVv->dE[P32][P12]   = 1.e+16  ;  CIVv->A[P32][P12]  = 1.e-31;  
*/

  Symmetrize_Coeff (CIVv);
}


/* ********************************************************************* */
void CV_INIT (Ion *CVv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int S0, P2, P2_1, P3_1;
  
  CVv->isCV = 1;
                                                                                                      
  CVv->N = 5.e-4;    /* abundance of C with respect to H (Raymond) */  
  CVv->nlev = 4;
 
  S0   = 0;
  P2   = 1;
  P2_1 = 2;
  P3_1 = 3;
 
  CVv->nTom  = 3;
  CVv->Tom[0] = 1.0e4;
  CVv->Tom[1] = 2.0e4;
  CVv->Tom[2] = 3.0e4;

  if (CVv->dE == NULL){
    CVv->dE    = ARRAY_2D(CVv->nlev, CVv->nlev, double);
    CVv->A     = ARRAY_2D(CVv->nlev, CVv->nlev, double);
    CVv->omega = ARRAY_3D(CVv->nlev, CVv->nlev, CVv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < CVv->nlev ; i++) {
  for (j = 0 ; j < CVv->nlev ; j++) {
    CVv->dE[i][j] = -1.0;
    CVv->A[i][j]  = -1.0;
    for (n = 0; n < CVv->nTom; n++){
      CVv->omega[i][j][n]  = -1.0;
    }
  }}

   
  CVv->wght[S0]   = 1.; /*    2 J  +  1     */
  CVv->wght[P2]   = 5.; /*         */
  CVv->wght[P2_1] = 3.; /*         */
  CVv->wght[P3_1] = 3.; /*         */

  CVv->omega[P2][S0][0] = 2.06e-2;  /* alpha      */
  CVv->omega[P2][S0][1] = 0.064;    /* beta       */
  CVv->omega[P2][S0][2] = 5.50;     /* log T_m    */
  
  CVv->omega[P2_1][S0][0] = 2.96e-2;  
  CVv->omega[P2_1][S0][1] = 0.282;  
  CVv->omega[P2_1][S0][2] = 5.50;  

  CVv->omega[P3_1][S0][0] = 1.09e-2;
  CVv->omega[P3_1][S0][1] = 0.094;
  CVv->omega[P3_1][S0][2] = 5.50;

  CVv->omega[P3_1][P2][0] = 0.00000;
  CVv->omega[P3_1][P2][1] = 0.094;
  CVv->omega[P3_1][P2][2] = 5.50;

  CVv->omega[P3_1][P2_1][0] = 0.00000;
  CVv->omega[P3_1][P2_1][1] = 0.094;
  CVv->omega[P3_1][P2_1][2] = 5.50;

  CVv->omega[P2_1][P2][0] = 0.00000;
  CVv->omega[P2_1][P2][1] = 0.094;
  CVv->omega[P2_1][P2][2] = 5.50;

  
  CVv->dE[P2][S0]   = 34.97;  CVv->A[P2][S0]   = 2.554e+11;  
  CVv->dE[P2_1][S0] = 40.27;  CVv->A[P2_1][S0] = 8.873e+11;
  CVv->dE[P3_1][S0] = 40.73;  CVv->A[P3_1][S0] = 2.62e+4;  
/*  
  CVv->dE[P3_1][P2]   = 1.0e+30; CVv->A[P3_1][P2]   = 0.00000;  
  CVv->dE[P3_1][P2_1] = 1.0e+30; CVv->A[P3_1][P2_1] = 0.00000;  
  CVv->dE[P2_1][P2]   = 1.0e+30; CVv->A[P2_1][P2]   = 0.00000;  
*/

  Symmetrize_Coeff (CVv);
}

/* ********************************************************************* */
void NI_INIT (Ion *NIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int S32, D52, D32, P32, P12;
                              
                              
  /* Data from T-dependent CS from MAPPINGS commented */
  
  NIv->isMAP = 0;
                                                                                                    
  NIv->N = 0.9e-4 ;    /* abundance of N with respect to H (Raymond) */  
  NIv->nlev = 5;
 
  S32 = 0;
  D52 = 1;
  D32 = 2;
  P12 = 3;
  P32 = 4;
 
  NIv->nTom  = 3;
  NIv->Tom[0] = 0.5e4;
  NIv->Tom[1] = 1.0e4;
  NIv->Tom[2] = 2.0e4;

  if (NIv->dE == NULL){
    NIv->dE    = ARRAY_2D(NIv->nlev, NIv->nlev, double);
    NIv->A     = ARRAY_2D(NIv->nlev, NIv->nlev, double);
    NIv->omega = ARRAY_3D(NIv->nlev, NIv->nlev, NIv->nTom, double);
  }


/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < NIv->nlev ; i++) {
  for (j = 0 ; j < NIv->nlev ; j++) {
    NIv->dE[i][j] = -1.0;
    NIv->A[i][j]  = -1.0;
    for (n = 0; n < NIv->nTom; n++){
      NIv->omega[i][j][n]  = -1.0;
    }
  }}

   
  NIv->wght[S32] = 4.; /*    2 J  +  1     */
  NIv->wght[D52] = 6.; /*         */
  NIv->wght[D32] = 4.; /*         */
  NIv->wght[P32] = 4.; /*         */
  NIv->wght[P12] = 2.; /*        */


  NIv->omega[D52][S32][0] = 1.55e-1;
  NIv->omega[D52][S32][1] = 2.90e-1;
  NIv->omega[D52][S32][2] = 4.76e-1;
 
  NIv->omega[D32][S32][0] = 1.03-1;
  NIv->omega[D32][S32][1] = 1.94-1;
  NIv->omega[D32][S32][2] = 3.18-1;
  
  NIv->omega[P32][S32][0] = 5.97e-2;
  NIv->omega[P32][S32][1] = 1.13e-1;
  NIv->omega[P32][S32][2] = 1.89e-1;
  
  NIv->omega[P12][S32][0] = 2.98e-2;
  NIv->omega[P12][S32][1] = 5.67e-2;
  NIv->omega[P12][S32][2] = 9.47e-2;
  
  NIv->omega[D52][D32][0] = 1.28e-1;
  NIv->omega[D52][D32][1] = 2.69e-1;
  NIv->omega[D52][D32][2] = 4.65e-1;

  NIv->omega[P32][P12][0] = 3.29e-2;
  NIv->omega[P32][P12][1] = 7.10e-2;
  NIv->omega[P32][P12][2] = 1.53e-1;

  NIv->omega[P32][D52][0] = 1.62e-1;
  NIv->omega[P32][D52][1] = 2.66e-1;
  NIv->omega[P32][D52][2] = 4.38e-1;

  NIv->omega[P32][D32][0] = 8.56e-2;
  NIv->omega[P32][D32][1] = 1.47e-1;
  NIv->omega[P32][D32][2] = 2.52e-1;

  NIv->omega[P12][D52][0] = 6.26e-2;
  NIv->omega[P12][D52][1] = 1.09e-1;
  NIv->omega[P12][D52][2] = 1.90e-1;

  NIv->omega[P12][D32][0] = 6.01e-2;
  NIv->omega[P12][D32][1] = 9.70e-2;
  NIv->omega[P12][D32][2] = 1.57e-1;

/*
  NIv->omega[D52][S32][0] = 2.90e-1;
  NIv->omega[D52][S32][1] = 9.10e-1;
 
  NIv->omega[D32][S32][0] = 1.94-1;
  NIv->omega[D32][S32][1] = 9.10-1;
  
  NIv->omega[P32][S32][0] = 1.13e-1;
  NIv->omega[P32][S32][1] = 9.30e-1;
  
  NIv->omega[P12][S32][0] = 5.67e-2;
  NIv->omega[P12][S32][1] = 9.30e-1;
  
  NIv->omega[D52][D32][0] = 2.69e-1;
  NIv->omega[D52][D32][1] = 1.08;

  NIv->omega[P32][P12][0] = 7.10e-2;
  NIv->omega[P32][P12][1] = 1.11;

  NIv->omega[P32][D52][0] = 2.66e-1;
  NIv->omega[P32][D52][1] = 0.72;

  NIv->omega[P32][D32][0] = 1.47e-1;
  NIv->omega[P32][D32][1] = 7.8e-1;

  NIv->omega[P12][D52][0] = 1.09e-1;
  NIv->omega[P12][D52][1] = 8.0e-1;

  NIv->omega[P12][D32][0] = 9.70e-2;
  NIv->omega[P12][D32][1] = 0.69;
*/

  NIv->dE[D52][S32]   = 5200.4;   NIv->A[D52][S32]  = 6.13e-6 ;  
  NIv->dE[D32][S32]   = 5197.9;   NIv->A[D32][S32]  = 2.28e-5;
  NIv->dE[P32][S32]   = 3466.5;   NIv->A[P32][S32]  = 6.60e-3;  
  NIv->dE[P12][S32]   = 3466.5;   NIv->A[P12][S32]  = 2.72e-3;
  NIv->dE[D52][D32]   = 1.148e+7; NIv->A[D52][D32]  = 1.24e-8;
  NIv->dE[P32][P12]   = 2.59e+8;  NIv->A[P32][P12]  = 5.17e-13;
  NIv->dE[P32][D52]   = 10397.7;  NIv->A[P32][D52]  = 5.59e-2;
  NIv->dE[P32][D32]   = 10407.2;  NIv->A[P32][D32]  = 2.52e-2;
  NIv->dE[P12][D52]   = 10398.2;  NIv->A[P12][D52]  = 3.14e-2;
  NIv->dE[P12][D32]   = 10407.6;  NIv->A[P12][D32]  = 4.80e-2; 

                                                                                                                  
  Symmetrize_Coeff (NIv);
}

/* ********************************************************************* */
void NII_INIT (Ion *NIIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int S0, P0, P1, P2, D2;
                 
  /*  Old data from Pradhan, commented. New T-dep data (Chebyshev pol fit) from Stafford et al. 1994  */

  NIIv->isCHEB = 0;  /* Chebyshev fit indicator  */
                                                                                     
  NIIv->N = 0.9e-4 ;    /* abundance of N with respect to H (Raymond) */  
  NIIv->nlev = 5;
 
  P0 = 0;
  P1 = 1;
  P2 = 2;
  D2 = 3;
  S0 = 4;
 
  NIIv->nTom  = 4;
  NIIv->Tom[0] = 0.5e4;
  NIIv->Tom[1] = 1.0e4;
  NIIv->Tom[2] = 1.5e4;
  NIIv->Tom[3] = 2.0e4;

  if (NIIv->dE == NULL){
    NIIv->dE    = ARRAY_2D(NIIv->nlev, NIIv->nlev, double);
    NIIv->A     = ARRAY_2D(NIIv->nlev, NIIv->nlev, double);
    NIIv->omega = ARRAY_3D(NIIv->nlev, NIIv->nlev, NIIv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < NIIv->nlev ; i++) {
  for (j = 0 ; j < NIIv->nlev ; j++) {
    NIIv->dE[i][j] = -1.0;
    NIIv->A[i][j]  = -1.0;
    for (n = 0; n < NIIv->nTom; n++){
      if (NIIv->isCHEB) NIIv->omega[i][j][n]  = 0.0;
      else NIIv->omega[i][j][n]  = -1.0;
    }
  }}

   
  NIIv->wght[S0] = 1.; /*    2 J  +  1     */
  NIIv->wght[P0] = 1.; /*         */
  NIIv->wght[P1] = 3.; /*         */
  NIIv->wght[P2] = 5.; /*         */
  NIIv->wght[D2] = 5.; /*        */

if (!NIIv->isCHEB) {
  NIIv->omega[D2][P0][0] = 2.57 * (1.0 / ( 3.0 * 3.0) );
  NIIv->omega[D2][P0][1] = 2.64 * (1.0 / ( 3.0 * 3.0) );
  NIIv->omega[D2][P0][2] = 2.70 * (1.0 / ( 3.0 * 3.0) );
  NIIv->omega[D2][P0][3] = 2.73 * (1.0 / ( 3.0 * 3.0) );

  NIIv->omega[D2][P1][0] = 2.57 * (3.0)/(3.0 * 3.0);
  NIIv->omega[D2][P1][1] = 2.64 * (3.0)/(3.0 * 3.0);
  NIIv->omega[D2][P1][2] = 2.70 * (3.0)/(3.0 * 3.0);
  NIIv->omega[D2][P1][3] = 2.73 * (3.0)/(3.0 * 3.0);

  NIIv->omega[D2][P2][0] = 2.57 * (5.0)/(3.0 * 3.0);
  NIIv->omega[D2][P2][1] = 2.64 * (5.0)/(3.0 * 3.0);
  NIIv->omega[D2][P2][2] = 2.70 * (5.0)/(3.0 * 3.0);
  NIIv->omega[D2][P2][3] = 2.73 * (5.0)/(3.0 * 3.0);

  NIIv->omega[S0][P1][0] = .287 * (3.0)/(3.0*3.0);
  NIIv->omega[S0][P1][1] = .293 * (3.0)/(3.0*3.0);
  NIIv->omega[S0][P1][2] = .300 * (3.0)/(3.0*3.0);
  NIIv->omega[S0][P1][3] = .305 * (3.0)/(3.0*3.0);

  NIIv->omega[S0][P2][0] = .287 * (5.0)/(3.0*3.0);
  NIIv->omega[S0][P2][1] = .293 * (5.0)/(3.0*3.0);
  NIIv->omega[S0][P2][2] = .300 * (5.0)/(3.0*3.0);
  NIIv->omega[S0][P2][3] = .305 * (5.0)/(3.0*3.0);

  NIIv->omega[S0][D2][0] = .959;
  NIIv->omega[S0][D2][1] = .834;
  NIIv->omega[S0][D2][2] = .761;
  NIIv->omega[S0][D2][3] = .734;

  NIIv->omega[P1][P0][0] = .371;
  NIIv->omega[P1][P0][1] = .408;
  NIIv->omega[P1][P0][2] = .429;
  NIIv->omega[P1][P0][3] = .443;
                                                                                                                        
  NIIv->omega[P2][P0][0] = .243;
  NIIv->omega[P2][P0][1] = .272;
  NIIv->omega[P2][P0][2] = .301;
  NIIv->omega[P2][P0][3] = .316;
                                                                                                                        
  NIIv->omega[P2][P1][0] = 1.01;
  NIIv->omega[P2][P1][1] = 1.12;
  NIIv->omega[P2][P1][2] = 1.21;
  NIIv->omega[P2][P1][3] = 1.26;
/*
  NIIv->omega[S0][P0][0]  = 1.e-31;
  NIIv->omega[S0][P0][1]  = 1.e-31;
  NIIv->omega[S0][P0][2]  = 1.e-31;
  NIIv->omega[S0][P0][3]  = 1.e-31;
*/
}
                                                                                                                        
if (NIIv->isCHEB) {
  NIIv->omega[D2][P0][0] = -2.1882;
  NIIv->omega[D2][P0][1] = -0.0023;
  NIIv->omega[D2][P0][2] = 0.0082;
  NIIv->omega[D2][P0][3] = -0.0005;

  NIIv->omega[D2][P1][0] = 0.0090;
  NIIv->omega[D2][P1][1] = -0.0023;
  NIIv->omega[D2][P1][2] = 0.0082;
  NIIv->omega[D2][P1][3] = -0.0006;

  NIIv->omega[D2][P2][0] = 1.0307;
  NIIv->omega[D2][P2][1] = -0.0023;
  NIIv->omega[D2][P2][2] = 0.0082;
  NIIv->omega[D2][P2][3] = -0.0006;

  NIIv->omega[S0][P1][0] = -4.1820;
  NIIv->omega[S0][P1][1] = -0.0076;
  NIIv->omega[S0][P1][2] = 0.0190;
  NIIv->omega[S0][P1][3] = -0.0013;

  NIIv->omega[S0][P2][0] = -3.1604;
  NIIv->omega[S0][P2][1] = -0.0076;
  NIIv->omega[S0][P2][2] = 0.0190;
  NIIv->omega[S0][P2][3] = -0.0013;

  NIIv->omega[S0][D2][0] = -1.2463;
  NIIv->omega[S0][D2][1] = 0.0961;
  NIIv->omega[S0][D2][2] = 0.0291;
  NIIv->omega[S0][D2][3] = -0.0024;

  NIIv->omega[P1][P0][0] = -1.5340;
  NIIv->omega[P1][P0][1] = 0.1436;
  NIIv->omega[P1][P0][2] = 0.0107;
  NIIv->omega[P1][P0][3] = -0.0076;
                                                                                                                        
  NIIv->omega[P2][P0][0] = -2.4004;
  NIIv->omega[P2][P0][1] = 0.1962;
  NIIv->omega[P2][P0][2] = -0.0440;
  NIIv->omega[P2][P0][3] = -0.0296;
                                                                                                                        
  NIIv->omega[P2][P1][0] = 0.4597;
  NIIv->omega[P2][P1][1] = 0.1717;
  NIIv->omega[P2][P1][2] = -0.0189;
  NIIv->omega[P2][P1][3] = -0.0198;

  NIIv->omega[S0][P0][0]  = -6.3792;
  NIIv->omega[S0][P0][1]  = -0.0076;
  NIIv->omega[S0][P0][2]  = 0.0190;
  NIIv->omega[S0][P0][3]  = -0.0013;

}

  NIIv->dE[D2][P2]   = 6583.4;  NIIv->A[D2][P2]  = 2.72e-3 ;  
  NIIv->dE[D2][P1]   = 6548.1;  NIIv->A[D2][P1]  = 9.19e-4;
  NIIv->dE[D2][P0]   = 6527.1;  NIIv->A[D2][P0]  = 5.4e-7;  
  NIIv->dE[P1][P0]   = 204.e4;  NIIv->A[P1][P0]  = 2.1e-6;
  NIIv->dE[P2][P0]   = 76.e4 ;  NIIv->A[P2][P0]  = 1.2e-12;
  NIIv->dE[P2][P1]   = 122.e4;  NIIv->A[P2][P1]  = 7.5e-6;
  NIIv->dE[S0][P1]   = 3062.9;  NIIv->A[S0][P1]  = 3.38e-2;
  NIIv->dE[S0][P2]   = 3071.4;  NIIv->A[S0][P2]  = 1.51e-4;
  NIIv->dE[S0][D2]   = 5754.6;  NIIv->A[S0][D2]  = 1.12;
/*  NIIv->dE[S0][P0]   = 2.e+21;  NIIv->A[S0][P0]  = 0.0000;    come S2-P2 (Pradhan)  */

                                                                                                                  
  Symmetrize_Coeff (NIIv);
}

/* ********************************************************************* */
void NIII_INIT (Ion *NIIIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int _2P12, _2P32, _4P12, _4P32, _4P52;
  
  /*  Old data from Pradhan, commented. New T-dep data (Chebyshev pol fit) from Stafford et al. 1994  */
                                                                                                    
  NIIIv->N = 0.9e-4;    /* abundance of N with respect to H (Raymond) */  
  NIIIv->nlev = 5;
 
  _2P12 = 0;
  _2P32 = 1;
  _4P12 = 2;
  _4P32 = 3;
  _4P52 = 4;

  NIIIv->isCHEB = 1; /* Use Chebyshev interpolation  ( log(T)<5.1 ) */ 
  NIIIv->isMAP  = 0; /* Use Mappings T-dep CS data   T < 10^7 K     */


  if (NIIIv->isCHEB) {
    NIIIv->nTom  = 4;
    NIIIv->Tom[0] = 0.5e4;
    NIIIv->Tom[1] = 1.0e4;
    NIIIv->Tom[2] = 1.5e4;
    NIIIv->Tom[3] = 2.0e4;
  }
  if (NIIIv->isMAP) {
    NIIIv->nTom  = 2;
    NIIIv->Tom[0] = 0.5e4;
    NIIIv->Tom[1] = 1.0e4;
  }

  if (NIIIv->dE == NULL){
    NIIIv->dE    = ARRAY_2D(NIIIv->nlev, NIIIv->nlev, double);
    NIIIv->A     = ARRAY_2D(NIIIv->nlev, NIIIv->nlev, double);
    NIIIv->omega = ARRAY_3D(NIIIv->nlev, NIIIv->nlev, NIIIv->nTom, double);
  }
 
  for (i = 0 ; i < NIIIv->nlev ; i++) {
  for (j = 0 ; j < NIIIv->nlev ; j++) {
    NIIIv->dE[i][j] = -1.0;
    NIIIv->A[i][j]  = -1.0;
    for (n = 0; n < NIIIv->nTom; n++){
      NIIIv->omega[i][j][n]  = 0.0;
    }
  }}

   
  NIIIv->wght[_2P12] = 2.; /*    2 J  +  1     */
  NIIIv->wght[_2P32] = 4.; /*         */
  NIIIv->wght[_4P12] = 2.; /*         */
  NIIIv->wght[_4P32] = 4.; /*         */
  NIIIv->wght[_4P52] = 6.; /*        */

/*
  NIIIv->omega[_2P32][_2P12][0] = 1.32;
  NIIIv->omega[_2P32][_2P12][1] = 1.45;
  NIIIv->omega[_2P32][_2P12][2] = 1.55;
  NIIIv->omega[_2P32][_2P12][3] = 1.64;

  NIIIv->omega[_4P12][_2P12][0] = 1.89e-1;
  NIIIv->omega[_4P12][_2P12][1] = 1.98e-1;
  NIIIv->omega[_4P12][_2P12][2] = 2.04e-1;
  NIIIv->omega[_4P12][_2P12][3] = 2.07e-1;

  NIIIv->omega[_4P12][_2P32][0] = 1.35e-1;
  NIIIv->omega[_4P12][_2P32][1] = 1.51e-1;
  NIIIv->omega[_4P12][_2P32][2] = 1.62e-1;
  NIIIv->omega[_4P12][_2P32][3] = 1.68e-1;
  
  NIIIv->omega[_4P32][_2P12][0] = 2.81e-1;
  NIIIv->omega[_4P32][_2P12][1] = 2.98e-1;
  NIIIv->omega[_4P32][_2P12][2] = 3.09e-1;
  NIIIv->omega[_4P32][_2P12][3] = 3.16e-1;
  
  NIIIv->omega[_4P32][_2P32][0] = 3.67e-1;
  NIIIv->omega[_4P32][_2P32][1] = 3.99e-1;
  NIIIv->omega[_4P32][_2P32][2] = 4.23e-1;
  NIIIv->omega[_4P32][_2P32][3] = 4.35e-1;

  NIIIv->omega[_4P32][_4P12][0] = 1.01;
  NIIIv->omega[_4P32][_4P12][1] = 1.10;
  NIIIv->omega[_4P32][_4P12][2] = 1.14;
  NIIIv->omega[_4P32][_4P12][3] = 1.16;

  NIIIv->omega[_4P52][_2P12][0] = 1.78e-1;
  NIIIv->omega[_4P52][_2P12][1] = 2.01e-1;
  NIIIv->omega[_4P52][_2P12][2] = 2.19e-1;
  NIIIv->omega[_4P52][_2P12][3] = 2.29e-1;

  NIIIv->omega[_4P52][_2P32][0] = 7.93e-1;
  NIIIv->omega[_4P52][_2P32][1] = 8.44e-1;
  NIIIv->omega[_4P52][_2P32][2] = 8.80e-1;
  NIIIv->omega[_4P52][_2P32][3] = 8.98e-1;

  NIIIv->omega[_4P52][_4P12][0] = 6.12e-1;
  NIIIv->omega[_4P52][_4P12][1] = 6.67e-1;
  NIIIv->omega[_4P52][_4P12][2] = 6.95e-1;
  NIIIv->omega[_4P52][_4P12][3] = 7.11e-1;

  NIIIv->omega[_4P52][_4P32][0] = 1.88;
  NIIIv->omega[_4P52][_4P32][0] = 2.04;
  NIIIv->omega[_4P52][_4P32][0] = 2.12;
  NIIIv->omega[_4P52][_4P32][0] = 2.16;
*/

if (NIIIv->isCHEB) {
  NIIIv->omega[_2P32][_2P12][0] = 0.7240;
  NIIIv->omega[_2P32][_2P12][1] = 0.2349;
  NIIIv->omega[_2P32][_2P12][2] = -0.0682;
  NIIIv->omega[_2P32][_2P12][3] = -0.0517;

  NIIIv->omega[_4P12][_2P12][0] = -3.6087;
  NIIIv->omega[_4P12][_2P12][1] = -0.0656;
  NIIIv->omega[_4P12][_2P12][2] = -0.0831;
  NIIIv->omega[_4P12][_2P12][3] = -0.0165;

  NIIIv->omega[_4P12][_2P32][0] = -4.1572;
  NIIIv->omega[_4P12][_2P32][1] = 0.0119;
  NIIIv->omega[_4P12][_2P32][2] = -0.1247;
  NIIIv->omega[_4P12][_2P32][3] = -0.0311;
  
  NIIIv->omega[_4P32][_2P12][0] = -2.7917;
  NIIIv->omega[_4P32][_2P12][1] = -0.0501;
  NIIIv->omega[_4P32][_2P12][2] = -0.0914;
  NIIIv->omega[_4P32][_2P12][3] = -0.0197;
  
  NIIIv->omega[_4P32][_2P32][0] = -2.2053;
  NIIIv->omega[_4P32][_2P32][1] = -0.0191;
  NIIIv->omega[_4P32][_2P32][2] = -0.1080;
  NIIIv->omega[_4P32][_2P32][3] = -0.0257;

  NIIIv->omega[_4P32][_4P12][0] = -0.0975;
  NIIIv->omega[_4P32][_4P12][1] = 0.0837;
  NIIIv->omega[_4P32][_4P12][2] = -0.0437;
  NIIIv->omega[_4P32][_4P12][3] = -0.0195;

  NIIIv->omega[_4P52][_2P12][0] = -3.5774;
  NIIIv->omega[_4P52][_2P12][1] = 0.0274;
  NIIIv->omega[_4P52][_2P12][2] = -0.1330;
  NIIIv->omega[_4P52][_2P12][3] = -0.0337;

  NIIIv->omega[_4P52][_2P32][0] = -0.7075;
  NIIIv->omega[_4P52][_2P32][1] = -0.0465;
  NIIIv->omega[_4P52][_2P32][2] = -0.0934;
  NIIIv->omega[_4P52][_2P32][3] = -0.0204;

  NIIIv->omega[_4P52][_4P12][0] = -0.7808;
  NIIIv->omega[_4P52][_4P12][1] = 0.1475;
  NIIIv->omega[_4P52][_4P12][2] = -0.0418;
  NIIIv->omega[_4P52][_4P12][3] = -0.0211;

  NIIIv->omega[_4P52][_4P32][0] = 1.2787;
  NIIIv->omega[_4P52][_4P32][0] = 0.1128;
  NIIIv->omega[_4P52][_4P32][0] = -0.0426;
  NIIIv->omega[_4P52][_4P32][0] = -0.0203;
}

if (NIIIv->isMAP) {  
  NIIIv->omega[_2P32][_2P12][0] = 1.4450;
  NIIIv->omega[_2P32][_2P12][1] = 0.1840;

  NIIIv->omega[_4P12][_2P12][0] = 0.1980;
  NIIIv->omega[_4P12][_2P12][1] = 0.0630;

  NIIIv->omega[_4P12][_2P32][0] = 0.1510;
  NIIIv->omega[_4P12][_2P32][1] = 0.1640;
  
  NIIIv->omega[_4P32][_2P12][0] = 0.2980;
  NIIIv->omega[_4P32][_2P12][1] = 0.0840;
  
  NIIIv->omega[_4P32][_2P32][0] = 0.3990;
  NIIIv->omega[_4P32][_2P32][1] = 0.1250;

  NIIIv->omega[_4P52][_2P32][0] = 0.8440;
  NIIIv->omega[_4P52][_2P32][1] = 0.0890;

  NIIIv->omega[_4P32][_4P12][0] = 1.102;
  NIIIv->omega[_4P32][_4P12][1] = 0.0720;

  NIIIv->omega[_4P52][_2P12][0] = 0.2010;
  NIIIv->omega[_4P52][_2P12][1] = 0.1830;

  NIIIv->omega[_4P52][_4P12][0] = 0.6680;
  NIIIv->omega[_4P52][_4P12][1] = 0.0910;

  NIIIv->omega[_4P52][_4P32][0] = 2.0440;
  NIIIv->omega[_4P52][_4P32][0] = 0.0800;

}
  NIIIv->dE[_2P32][_2P12]   = 5.73e+5;  NIIIv->A[_2P32][_2P12]  = 4.77e-5;  
  NIIIv->dE[_4P12][_2P12]   = 1748.;    NIIIv->A[_4P12][_2P12]  = 3.08e+2;
  NIIIv->dE[_4P12][_2P32]   = 1754.;    NIIIv->A[_4P12][_2P32]  = 5.22e+2;  
  NIIIv->dE[_4P32][_2P12]   = 1747.;    NIIIv->A[_4P32][_2P12]  = 9.49;
  NIIIv->dE[_4P32][_2P32]   = 1752.;    NIIIv->A[_4P32][_2P32]  = 6.50e+1;
  NIIIv->dE[_4P52][_2P32]   = 1747.;    NIIIv->A[_4P52][_2P32]  = 3.08e+2;
  NIIIv->dE[_4P32][_4P12]   = 1.68e+6;  NIIIv->A[_4P32][_4P12]  = 0.00000;
  NIIIv->dE[_4P52][_2P12]   = 1744.4;   NIIIv->A[_4P52][_2P12]  = 0.00000;
  NIIIv->dE[_4P52][_4P12]   = 7.10e+5;  NIIIv->A[_4P52][_4P12]  = 0.00000;
  NIIIv->dE[_4P52][_4P32]   = 1.23e+6; NIIIv->A[_4P52][_4P32]  = 0.00000;


  Symmetrize_Coeff (NIIIv);
}

/* ********************************************************************* */
void NIV_INIT (Ion *NIVv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int S0, P1, P2, P1_1;
                                                                                                    
  NIVv->N = 0.9e-4 ;    /* abundance of N with respect to H (Raymond) */  
  NIVv->nlev = 4;
 
  S0 = 0;
  P1 = 1;
  P2 = 2;
  P1_1 = 3;
 
  NIVv->nTom  = 4;
  NIVv->Tom[0] = 0.5e4;
  NIVv->Tom[1] = 1.0e4;
  NIVv->Tom[2] = 1.5e4;
  NIVv->Tom[3] = 2.0e4;

  if (NIVv->dE == NULL){
    NIVv->dE    = ARRAY_2D(NIVv->nlev, NIVv->nlev, double);
    NIVv->A     = ARRAY_2D(NIVv->nlev, NIVv->nlev, double);
    NIVv->omega = ARRAY_3D(NIVv->nlev, NIVv->nlev, NIVv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < NIVv->nlev ; i++) {
  for (j = 0 ; j < NIVv->nlev ; j++) {
    NIVv->dE[i][j] = -1.0;
    NIVv->A[i][j]  = -1.0;
    for (n = 0; n < NIVv->nTom; n++){
      NIVv->omega[i][j][n]  = -1.0;
    }
  }}

   
  NIVv->wght[S0] = 1.; /*    2 J  +  1     */
  NIVv->wght[P1] = 3.; /*         */
  NIVv->wght[P2] = 5.; /*         */
  NIVv->wght[P1_1] = 3.; /*        */


  NIVv->omega[P2][S0][0] = 9.37e-1 * (5.0)/(3.0*3.0);
  NIVv->omega[P2][S0][1] = 9.05e-1 * (5.0)/(3.0*3.0);
  NIVv->omega[P2][S0][2] = 8.79e-1 * (5.0)/(3.0*3.0);
  NIVv->omega[P2][S0][3] = 8.58e-1 * (5.0)/(3.0*3.0);
  
  NIVv->omega[P1][S0][0] = 9.37e-1 *(3.0)/(3.0 * 3.0);
  NIVv->omega[P1][S0][1] = 9.05e-1 *(3.0)/(3.0 * 3.0);
  NIVv->omega[P1][S0][2] = 9.05e-1 *(3.0)/(3.0 * 3.0);
  NIVv->omega[P1][S0][3] = 8.58e-1 *(3.0)/(3.0 * 3.0);
  
  NIVv->omega[P1_1][S0][0] = 3.84;
  NIVv->omega[P1_1][S0][1] = 3.53;
  NIVv->omega[P1_1][S0][2] = 3.41;
  NIVv->omega[P1_1][S0][3] = 3.36;

  /* fake data */
/*  
  NIVv->omega[P2][P1][0] = 9.37e-31;
  NIVv->omega[P2][P1][1] = 9.05e-31;
  NIVv->omega[P2][P1][2] = 8.79e-31;
  NIVv->omega[P2][P1][3] = 8.58e-31;

  NIVv->omega[P1_1][P1][0] = 9.37e-31;
  NIVv->omega[P1_1][P1][1] = 9.05e-31;
  NIVv->omega[P1_1][P1][2] = 8.79e-31;
  NIVv->omega[P1_1][P1][3] = 8.58e-31;

  NIVv->omega[P1_1][P2][0] = 9.37e-31;
  NIVv->omega[P1_1][P2][1] = 9.05e-31;
  NIVv->omega[P1_1][P2][2] = 8.79e-31;
  NIVv->omega[P1_1][P2][3] = 8.58e-31;
*/

  NIVv->dE[P2][S0]   = 1483.3;  NIVv->A[P2][S0]  = 3.e-3 ;  
  NIVv->dE[P1][S0]   = 1486.4;  NIVv->A[P1][S0]  = 1.e-3;
  NIVv->dE[P1_1][S0] = 765.15;  NIVv->A[P1_1][S0]  = 5.4e-7;  
/*  NIVv->dE[P2][P1]   = 6.94e+5;  NIVv->A[P2][P1]  = 2.1e-31;
  NIVv->dE[P1_1][P1] = 6.94e+5;  NIVv->A[P1_1][P1]  = 1.2e-31;
  NIVv->dE[P1_1][P2] = 6.94e+5;  NIVv->A[P1_1][P2]  = 7.5e-31;
*/
                                                                                                                  
  Symmetrize_Coeff (NIVv);
}

/* ********************************************************************* */
void NV_INIT (Ion *NVv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int S, P;
                                                                                                    
  NVv->N = 0.9e-4 ;    /* abundance of N with respect to H (Raymond) */  
  NVv->nlev = 2;
 
  S = 0;
  P = 1;
 
  NVv->nTom  = 4;
  NVv->Tom[0] = 0.5e4;
  NVv->Tom[1] = 1.0e4;
  NVv->Tom[2] = 1.5e4;
  NVv->Tom[3] = 2.0e4;

  if (NVv->dE == NULL){
    NVv->dE    = ARRAY_2D(NVv->nlev, NVv->nlev, double);
    NVv->A     = ARRAY_2D(NVv->nlev, NVv->nlev, double);
    NVv->omega = ARRAY_3D(NVv->nlev, NVv->nlev, NVv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < NVv->nlev ; i++) {
  for (j = 0 ; j < NVv->nlev ; j++) {
    NVv->dE[i][j] = -1.0;
    NVv->A[i][j]  = -1.0;
    for (n = 0; n < NVv->nTom; n++){
      NVv->omega[i][j][n]  = -1.0;
    }
  }}

   
  NVv->wght[S] = 2.; /*    2 J  +  1     */
  NVv->wght[P] = 4.; /*         */
/*  NVv->wght[P2] = 5.; 
  NVv->wght[P1_1] = 3.;         */


  NVv->omega[P][S][0] = 6.61;
  NVv->omega[P][S][1] = 6.65;
  NVv->omega[P][S][2] = 6.69;
  NVv->omega[P][S][3] = 6.72;

  NVv->dE[P][S]   = 1238.8;  NVv->A[P][S]  = 3.41e+8 ;  
                                                                                                                  
  Symmetrize_Coeff (NVv);
}

/* ********************************************************************* */
void OI_INIT (Ion *OIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int S0, P0, P1, P2, D2;
                                                                                                    
  OIv->N = 6.e-4 ;    /* abundance of O with respect to H (Raymond) */  
  OIv->nlev = 5;
  
  OIv->isMAP = 1;
 
  P0 = 2;
  P1 = 1;
  P2 = 0;
  D2 = 3;
  S0 = 4;
 
  OIv->nTom  = 2;
  OIv->Tom[0] = 0.5e4;
  OIv->Tom[1] = 1.0e4;

  if (OIv->dE == NULL){
    OIv->dE    = ARRAY_2D(OIv->nlev, OIv->nlev, double);
    OIv->A     = ARRAY_2D(OIv->nlev, OIv->nlev, double);
    OIv->omega = ARRAY_3D(OIv->nlev, OIv->nlev, OIv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < OIv->nlev ; i++) {
  for (j = 0 ; j < OIv->nlev ; j++) {
    OIv->dE[i][j] = -1.0;
    OIv->A[i][j]  = -1.0;
    for (n = 0; n < OIv->nTom; n++){
      OIv->omega[i][j][n]  = 0.0;
    }
  }}

   
  OIv->wght[S0] = 1.; /*    2 J  +  1     */
  OIv->wght[P0] = 1.; /*         */
  OIv->wght[P1] = 3.; /*         */
  OIv->wght[P2] = 5.; /*         */
  OIv->wght[D2] = 5.; /*        */


/* OLD data from Pradhan commented  -  T-dep data from MAPPINGS added  */

/*
  OIv->omega[D2][P0][0] = 1.24e-1;
  OIv->omega[D2][P0][1] = 2.66e-1;
  OIv->omega[D2][P0][2] = 5.01e-1;

  OIv->omega[D2][P1][0] = OIv->omega[D2][P0][0]*(3.0)/(3.0 * 3.0);
  OIv->omega[D2][P1][1] = OIv->omega[D2][P0][1]*(3.0)/(3.0 * 3.0);
  OIv->omega[D2][P1][2] = OIv->omega[D2][P0][2]*(3.0)/(3.0 * 3.0);

  OIv->omega[D2][P2][0] = OIv->omega[D2][P0][0]*(5.0)/(3.0 * 3.0);
  OIv->omega[D2][P2][1] = OIv->omega[D2][P0][1]*(5.0)/(3.0 * 3.0);
  OIv->omega[D2][P2][2] = OIv->omega[D2][P0][2]*(5.0)/(3.0 * 3.0);

  OIv->omega[S0][P1][0] = 1.53e-2;
  OIv->omega[S0][P1][1] = 3.24e-2;
  OIv->omega[S0][P1][2] = 6.07e-2;

  OIv->omega[S0][P2][0] = OIv->omega[S0][P1][0]*(5.0)/(3.0*3.0);
  OIv->omega[S0][P2][1] = OIv->omega[S0][P1][1]*(5.0)/(3.0*3.0);
  OIv->omega[S0][P2][2] = OIv->omega[S0][P1][2]*(5.0)/(3.0*3.0);

  OIv->omega[S0][D2][0] = 7.32e-2;
  OIv->omega[S0][D2][1] = 1.05e-1;
  OIv->omega[S0][D2][2] = 1.48e-1;

  OIv->omega[P0][P1][0] = 1.12e-2;
  OIv->omega[P0][P1][1] = 2.65e-2;
  OIv->omega[P0][P1][2] = 6.93e-2;

  OIv->omega[P0][P2][0] = 1.48e-2;
  OIv->omega[P0][P2][1] = 2.92e-2;
  OIv->omega[P0][P2][2] = 5.36e-2;

  OIv->omega[P1][P2][0] = 4.74e-2;
  OIv->omega[P1][P2][1] = 9.87e-2;
  OIv->omega[P1][P2][2] = 2.07e-1;

  OIv->omega[S0][P0][0]  = 1.0e-30;
  OIv->omega[S0][P0][1]  = 1.0e-30;
  OIv->omega[S0][P0][2]  = 1.0e-30;
*/

  OIv->omega[D2][P0][0] = 2.955e-2;
  OIv->omega[D2][P0][1] = 1.1;
  
  OIv->omega[D2][P1][0] = 8.866e-2;
  OIv->omega[D2][P1][1] = 1.1;
  
  OIv->omega[D2][P2][0] = 1.4778e-1;
  OIv->omega[D2][P2][1] = 1.1;
  
  OIv->omega[S0][P1][0] = 1.08e-2;
  OIv->omega[S0][P1][1] = 1.08;
  
  OIv->omega[S0][P2][0] =  1.80e-2;
  OIv->omega[S0][P2][1] =  1.08;
    
  OIv->omega[S0][D2][0] = 0.10500;
  OIv->omega[S0][D2][1] = 0.52;
  
  OIv->omega[P0][P1][0] = 3.49e-2;
  OIv->omega[P0][P1][1] = 1.2;
  
  OIv->omega[P0][P2][0] = 0.02710;
  OIv->omega[P0][P2][1] = 0.94;
  
  OIv->omega[P1][P2][0] = 0.10500;
  OIv->omega[P1][P2][1] = 1.02;
  
/*
  OIv->omega[S0][P0][0]  = 0.00360;
  OIv->omega[S0][P0][1]  = 1.08;
*/

  OIv->dE[D2][P2]   = 6300.3;  OIv->A[D2][P2]  = 6.34e-3;  /* dE [ HIGHER ][LOWER ]  */
  OIv->dE[D2][P1]   = 6363.8;  OIv->A[D2][P1]  = 2.11e-3;
  OIv->dE[D2][P0]   = 6393.5;  OIv->A[D2][P0]  = 7.23e-7;
  OIv->dE[P1][P0]   = 1.455e+6;OIv->A[P1][P0]  = 1.74e-5;
  OIv->dE[P2][P0]   = 4.406e+5;OIv->A[P2][P0]  = 1.0e-10;
  OIv->dE[P2][P1]   = 6.318e+5;OIv->A[P2][P1]  = 8.92e-5;
  OIv->dE[S0][P1]   = 2972.3;  OIv->A[S0][P1]  = 7.32e-2;
  OIv->dE[S0][P2]   = 2959.2;  OIv->A[S0][P2]  = 2.88e-4;
  OIv->dE[S0][D2]   = 5577.3;  OIv->A[S0][D2]  = 1.22e+0;
/*  OIv->dE[S0][P0]   = 2.e+30;  OIv->A[S0][P0]  = 0.00000;      fittizio  */
                                                                                                                  
  Symmetrize_Coeff (OIv);
}

/* ********************************************************************* */
void OII_INIT (Ion *OIIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int S32, D52, D32, P32, P12;
                                                                                                    
  OIIv->N = 6.e-4 ;    /* abundance of O with respect to H (Raymond) */  
  OIIv->nlev = 5;
 
  S32 = 0;
  D52 = 1;
  D32 = 2;
  P32 = 3;
  P12 = 4;
 
  OIIv->nTom  = 4;
  OIIv->Tom[0] = 0.5e4;
  OIIv->Tom[1] = 1.0e4;
  OIIv->Tom[2] = 1.5e4;
  OIIv->Tom[3] = 2.0e4;

  if (OIIv->dE == NULL){
    OIIv->dE    = ARRAY_2D(OIIv->nlev, OIIv->nlev, double);
    OIIv->A     = ARRAY_2D(OIIv->nlev, OIIv->nlev, double);
    OIIv->omega = ARRAY_3D(OIIv->nlev, OIIv->nlev, OIIv->nTom, double);
  }
  
/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < OIIv->nlev ; i++) {
  for (j = 0 ; j < OIIv->nlev ; j++) {
    OIIv->dE[i][j] = -1.0;
    OIIv->A[i][j]  = -1.0;
    for (n = 0; n < OIIv->nTom; n++){
      OIIv->omega[i][j][n]  = -1.0;
    }
  }}

   
  OIIv->wght[S32] = 4.; /*    2 J  +  1     */
  OIIv->wght[D52] = 6.; /*         */
  OIIv->wght[D32] = 4.; /*         */
  OIIv->wght[P32] = 4.; /*         */
  OIIv->wght[P12] = 2.; /*        */


  OIIv->omega[D52][S32][0] = 7.95e-1;
  OIIv->omega[D52][S32][1] = 8.01e-1;
  OIIv->omega[D52][S32][2] = 8.10e-1;
  OIIv->omega[D52][S32][3] = 8.18e-1;

  OIIv->omega[D32][S32][0] = 5.30e-1;
  OIIv->omega[D32][S32][1] = 5.34e-1;
  OIIv->omega[D32][S32][2] = 5.41e-1;
  OIIv->omega[D32][S32][3] = 5.45e-1;

  OIIv->omega[P32][S32][0] = 2.65e-1;
  OIIv->omega[P32][S32][1] = 2.70e-1;
  OIIv->omega[P32][S32][2] = 2.75e-1;
  OIIv->omega[P32][S32][3] = 2.80e-1;
  
  OIIv->omega[P12][S32][0] = 1.33e-1;
  OIIv->omega[P12][S32][1] = 1.35e-1;
  OIIv->omega[P12][S32][2] = 1.37e-1;
  OIIv->omega[P12][S32][3] = 1.40e-1;

  OIIv->omega[D52][D32][0] = 1.22;
  OIIv->omega[D52][D32][1] = 1.17;
  OIIv->omega[D52][D32][2] = 1.14;
  OIIv->omega[D52][D32][3] = 1.11;

  OIIv->omega[P32][P12][0] = 2.80e-1;
  OIIv->omega[P32][P12][1] = 2.87e-1;
  OIIv->omega[P32][P12][2] = 2.93e-1;
  OIIv->omega[P32][P12][3] = 3.00e-1;

  OIIv->omega[P32][D52][0] = 7.18e-1;
  OIIv->omega[P32][D52][1] = 7.30e-1;
  OIIv->omega[P32][D52][2] = 7.41e-1;
  OIIv->omega[P32][D52][3] = 7.55e-1;
  
  OIIv->omega[P32][D32][0] = 4.01e-1;
  OIIv->omega[P32][D32][1] = 4.08e-1;
  OIIv->omega[P32][D32][2] = 4.14e-1;
  OIIv->omega[P32][D32][3] = 4.22e-1;

  OIIv->omega[P12][D52][0] = 2.90e-1;
  OIIv->omega[P12][D52][1] = 2.95e-1;
  OIIv->omega[P12][D52][2] = 3.00e-1;
  OIIv->omega[P12][D52][3] = 3.05e-1;

  OIIv->omega[P12][D32][0] = 2.70e-1;
  OIIv->omega[P12][D32][1] = 2.75e-1;
  OIIv->omega[P12][D32][2] = 2.81e-1;
  OIIv->omega[P12][D32][3] = 2.84e-1;

  OIIv->dE[D52][S32]   = 3728.8;  OIIv->A[D52][S32]  = 3.50e-5;  /* dE [ HIGHER ][LOWER ]  */
  OIIv->dE[D32][S32]   = 3726.0;  OIIv->A[D32][S32]  = 1.79e-4;
  OIIv->dE[P32][S32]   = 2470.3;  OIIv->A[P32][S32]  = 5.70e-2;
  OIIv->dE[P12][S32]   = 2470.2;  OIIv->A[P12][S32]  = 2.34e-2;
  OIIv->dE[D52][D32]   = 4.97e+6;  OIIv->A[D52][D32]  = 1.30e-7;
  OIIv->dE[P32][P12]   = 5.00e+6;  OIIv->A[P32][P12]  = 2.08e-11;
  OIIv->dE[P32][D52]   = 7319.9;  OIIv->A[P32][D52]  = 1.07e-1;
  OIIv->dE[P32][D32]   = 7330.7;  OIIv->A[P32][D32]  = 5.78e-2;
  OIIv->dE[P12][D52]   = 7321.8;  OIIv->A[P12][D52]  = 6.15e-2;
  OIIv->dE[P12][D32]   = 7329.6;   OIIv->A[P12][D32]  = 1.02e-1;    
                                                                                                                  
  Symmetrize_Coeff (OIIv);
}

/* ********************************************************************* */
void OIII_INIT (Ion *OIIIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int P0, P1, P2, D2, S0;

  /*  MAPPINGS data used for this one  */
                                                                                                    
  OIIIv->N = 6.e-4;    /* abundance of O with respect to H (Raymond) */  
  OIIIv->nlev = 5;
 
  P0 = 0;
  P1 = 1;
  P2 = 2;
  D2 = 3;
  S0 = 4;

  OIIIv->isMAP = 1;
 
  OIIIv->nTom  = 2;
  OIIIv->Tom[0] = 0.5e4;
  OIIIv->Tom[1] = 1.5e4;


  if (OIIIv->dE == NULL){
    OIIIv->dE    = ARRAY_2D(OIIIv->nlev, OIIIv->nlev, double);
    OIIIv->A     = ARRAY_2D(OIIIv->nlev, OIIIv->nlev, double);
    OIIIv->omega = ARRAY_3D(OIIIv->nlev, OIIIv->nlev, OIIIv->nTom, double);
  }
  
/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < OIIIv->nlev ; i++) {
  for (j = 0 ; j < OIIIv->nlev ; j++) {
    OIIIv->dE[i][j] = -1.0;
    OIIIv->A[i][j]  = -1.0;
    for (n = 0; n < OIIIv->nTom; n++){
      OIIIv->omega[i][j][n]  = 0.0;
    }
  }}
   
  OIIIv->wght[P0] = 1.; /*    2 J  +  1     */
  OIIIv->wght[P1] = 3.; /*         */
  OIIIv->wght[P2] = 5.; /*         */
  OIIIv->wght[D2] = 5.; /*         */
  OIIIv->wght[S0] = 1.; /*         */


  OIIIv->omega[D2][P0][0] = 0.2411;
  OIIIv->omega[D2][P0][1] = 0.1400;
  
  OIIIv->omega[D2][P1][0] = 0.7233;
  OIIIv->omega[D2][P1][1] = 0.1400;

  OIIIv->omega[D2][P2][0] = 1.2056;
  OIIIv->omega[D2][P2][1] = 0.1400;

  OIIIv->omega[S0][P1][0] = 0.0902;
  OIIIv->omega[S0][P1][1] = 0.1860;

  OIIIv->omega[S0][P2][0] = 0.1530;
  OIIIv->omega[S0][P2][1] = 0.1860;

  OIIIv->omega[S0][D2][0] = 0.6200;
  OIIIv->omega[S0][D2][1] = 0.0000;

  OIIIv->omega[P1][P0][0] = 0.5420;
  OIIIv->omega[P1][P0][1] = 0.0000;

  OIIIv->omega[P2][P0][0] = 0.2710;
  OIIIv->omega[P2][P0][1] = 0.0000;

  OIIIv->omega[P2][P1][0] = 1.2900;
  OIIIv->omega[P2][P1][1] = 0.0000;

  OIIIv->omega[S0][P0][0] = 0.0307;
  OIIIv->omega[S0][P0][1] = 0.1860;


 /*   Aggarwal '93 CS data  (up to 40000 K)
  OIIIv->omega[D][P][0] = 2.039;
  OIIIv->omega[D][P][1] = 2.191;
  OIIIv->omega[D][P][2] = 2.519;
  OIIIv->omega[D][P][3] = 2.582;
  
  OIIIv->omega[S][P][0] = 2.732e-1;
  OIIIv->omega[S][P][1] = 2.885e-1;
  OIIIv->omega[S][P][2] = 3.404e-1;
  OIIIv->omega[S][P][3] = 3.550e-1;

  OIIIv->omega[S][D][0] = 4.312e-1;
  OIIIv->omega[S][D][1] = 5.227e-1;
  OIIIv->omega[S][D][2] = 5.812e-1;
  OIIIv->omega[S][D][3] = 5.745e-1;  
*/

  OIIIv->dE[D2][P0] = 4932.6;  OIIIv->A[D2][P0] = 2.74e-6;  /* dE [ HIGHER ][LOWER ]  */
  OIIIv->dE[D2][P1] = 4958.9;  OIIIv->A[D2][P1] = 6.74e-3;
  OIIIv->dE[D2][P2] = 5006.7;  OIIIv->A[D2][P2] = 1.96e-2;
  OIIIv->dE[S0][P1] = 2321.0;  OIIIv->A[S0][P1] = 2.17e-1;
  OIIIv->dE[S0][P2] = 2332.1;  OIIIv->A[S0][P2] = 7.85e-4;
  OIIIv->dE[S0][D2] = 4363.2;  OIIIv->A[S0][D2] = 1.78;
  OIIIv->dE[P1][P0] = 883562.; OIIIv->A[P1][P0] = 2.62e-5;
  OIIIv->dE[P2][P0] = 326611.; OIIIv->A[P2][P0] = 3.02e-11;
  OIIIv->dE[P2][P1] = 518145.; OIIIv->A[P2][P1] = 9.76e-5;
/*  OIIIv->dE[S0][P0] = 1.e+20;  OIIIv->A[S0][P0] = 0.0000;
*/ 
                                                                                                                  
  Symmetrize_Coeff (OIIIv);
}

/* ********************************************************************* */
void OIV_INIT (Ion *OIVv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int _2P12, _2P32, _4P12, _4P32, _4P52;
   
   /* MAPPINGS CS data here  */     

  OIVv->isMAP = 1;
                                                                                                    
  OIVv->N = 6.e-4 ;    /* abundance of O with respect to H (Raymond) */  
  OIVv->nlev = 5;
 
  _2P12 = 0;
  _2P32 = 1;
  _4P12 = 2;
  _4P32 = 3;
  _4P52 = 4;
 
  OIVv->nTom  = 2;
  OIVv->Tom[0] = 0.5e4;
  OIVv->Tom[1] = 1.0e4;

  if (OIVv->dE == NULL){
    OIVv->dE    = ARRAY_2D(OIVv->nlev, OIVv->nlev, double);
    OIVv->A     = ARRAY_2D(OIVv->nlev, OIVv->nlev, double);
    OIVv->omega = ARRAY_3D(OIVv->nlev, OIVv->nlev, OIVv->nTom, double);
  }
  
/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < OIVv->nlev ; i++) {
  for (j = 0 ; j < OIVv->nlev ; j++) {
    OIVv->dE[i][j] = -1.0;
    OIVv->A[i][j]  = -1.0;
    for (n = 0; n < OIVv->nTom; n++){
      OIVv->omega[i][j][n]  = 0.0;
    }
  }}

   
  OIVv->wght[_2P12] = 2.; /*    2 J  +  1     */
  OIVv->wght[_2P32] = 4.; /*         */
  OIVv->wght[_4P12] = 2.; /*         */
  OIVv->wght[_4P32] = 4.; /*         */
  OIVv->wght[_4P52] = 6.; /*        */


  OIVv->omega[_2P32][_2P12][0] = 2.4200;
  OIVv->omega[_2P32][_2P12][1] = 0.0920;

  OIVv->omega[_4P12][_2P12][0] = 0.1330;
  OIVv->omega[_4P12][_2P12][1] = 0.1610;

  OIVv->omega[_4P12][_2P32][0] = 0.1020;
  OIVv->omega[_4P12][_2P32][1] = 0.2940;

  OIVv->omega[_4P32][_2P12][0] = 0.2010;
  OIVv->omega[_4P32][_2P12][1] = 0.1890;

  OIVv->omega[_4P32][_2P32][0] = 0.2700;
  OIVv->omega[_4P32][_2P32][1] = 0.2430;

  OIVv->omega[_4P32][_4P12][0] = 1.0880;
  OIVv->omega[_4P32][_4P12][1] = 0.1060;

  OIVv->omega[_4P52][_2P12][0] = 0.1370;
  OIVv->omega[_4P52][_2P12][1] = 0.3200;

  OIVv->omega[_4P52][_2P32][0] = 0.5690;
  OIVv->omega[_4P52][_2P32][1] = 0.1950;

  OIVv->omega[_4P52][_4P12][0] = 0.6870;
  OIVv->omega[_4P52][_4P12][1] = 0.1070;

  OIVv->omega[_4P52][_4P32][0] = 2.0540;
  OIVv->omega[_4P52][_4P32][1] = 0.1070;

  
  OIVv->dE[_2P32][_2P12] = 2.587e+5; OIVv->A[_2P32][_2P12] = 5.18e-4;
  OIVv->dE[_4P12][_2P12] = 1426.46;  OIVv->A[_4P12][_2P12] = 1.54e+3;
  OIVv->dE[_4P12][_2P32] = 1434.07;  OIVv->A[_4P12][_2P32] = 1.56e+3;
  OIVv->dE[_4P32][_2P12] = 1423.84;  OIVv->A[_4P32][_2P12] = 4.50e+1;
  OIVv->dE[_4P32][_2P32] = 1431.42;  OIVv->A[_4P32][_2P32] = 2.56e+2;
  OIVv->dE[_4P32][_4P12] = 1.68e+6;  OIVv->A[_4P32][_4P12] = 5.070e-5;
  OIVv->dE[_4P52][_2P12] = 1420.19;  OIVv->A[_4P52][_2P12] = 0.00000;
  OIVv->dE[_4P52][_2P32] = 1427.78;  OIVv->A[_4P52][_2P32] = 1.24e+3;
  OIVv->dE[_4P52][_4P12] = 3.26e+5;  OIVv->A[_4P52][_4P12] = 0.00000;
  OIVv->dE[_4P52][_4P32] = 5.62e+5;  OIVv->A[_4P52][_4P32] = 1.02e-4;
                                                                                                                  
  Symmetrize_Coeff (OIVv);
}

/* ********************************************************************* */
void OV_INIT (Ion *OVv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int S0, P0, P1, P2, P11;
                         
  /*  Pradhan data here  */                      
                                                                                                    
  OVv->N = 6.e-4 ;    /* abundance of O with respect to H (Raymond) */  
  OVv->nlev = 5;
 
  S0 = 0;
  P0 = 1;
  P1 = 2;
  P2 = 3;
  P11 = 4;
 
  OVv->nTom  = 4;
  OVv->Tom[0] = 0.5e4;
  OVv->Tom[1] = 1.0e4;
  OVv->Tom[2] = 1.5e4;
  OVv->Tom[3] = 2.0e4;

  if (OVv->dE == NULL){
    OVv->dE    = ARRAY_2D(OVv->nlev, OVv->nlev, double);
    OVv->A     = ARRAY_2D(OVv->nlev, OVv->nlev, double);
    OVv->omega = ARRAY_3D(OVv->nlev, OVv->nlev, OVv->nTom, double);
  }
  
/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < OVv->nlev ; i++) {
  for (j = 0 ; j < OVv->nlev ; j++) {
    OVv->dE[i][j] = -1.0;
    OVv->A[i][j]  = -1.0;
    for (n = 0; n < OVv->nTom; n++){
      OVv->omega[i][j][n]  = -1.0;
    }
  }}

   
  OVv->wght[S0] = 1.; /*    2 J  +  1     */
  OVv->wght[P0] = 1.; /*         */
  OVv->wght[P1] = 3.; /*         */
  OVv->wght[P2] = 5.; /*         */
  OVv->wght[P11] = 3.; /*        */


  OVv->omega[P2][S0][0] = 7.33e-1 * (5.0)/(3.0*3.0);
  OVv->omega[P2][S0][1] = 7.21e-1 * (5.0)/(3.0*3.0);
  OVv->omega[P2][S0][2] = 6.74e-1 * (5.0)/(3.0*3.0);
  OVv->omega[P2][S0][3] = 6.39e-1 * (5.0)/(3.0*3.0);

  OVv->omega[P1][S0][0] = 7.33e-1 * (3.0)/(3.0*3.0);
  OVv->omega[P1][S0][1] = 7.21e-1 * (3.0)/(3.0*3.0);
  OVv->omega[P1][S0][2] = 6.74e-1 * (3.0)/(3.0*3.0);
  OVv->omega[P1][S0][3] = 6.39e-1 * (3.0)/(3.0*3.0);

  OVv->omega[P0][S0][0] = 7.33e-1 * (1.0)/(3.0*3.0);
  OVv->omega[P0][S0][1] = 7.21e-1 * (1.0)/(3.0*3.0);
  OVv->omega[P0][S0][2] = 6.74e-1 * (1.0)/(3.0*3.0);
  OVv->omega[P0][S0][3] = 6.39e-1 * (1.0)/(3.0*3.0);

  OVv->omega[P11][S0][0] = 2.66;
  OVv->omega[P11][S0][1] = 2.76;
  OVv->omega[P11][S0][2] = 2.82;
  OVv->omega[P11][S0][3] = 2.85;

  OVv->omega[P1][P0][0] = 7.26e-1;
  OVv->omega[P1][P0][1] = 8.39e-1;
  OVv->omega[P1][P0][2] = 8.65e-1;
  OVv->omega[P1][P0][3] = 8.66e-1;

  OVv->omega[P2][P0][0] = 2.74e-1;
  OVv->omega[P2][P0][1] = 6.02e-1;
  OVv->omega[P2][P0][2] = 7.51e-1;
  OVv->omega[P2][P0][3] = 8.16e-1;

  OVv->omega[P2][P1][0] = 3.19;
  OVv->omega[P2][P1][1] = 2.86;
  OVv->omega[P2][P1][2] = 2.80;
  OVv->omega[P2][P1][3] = 2.77;
/*
  OVv->omega[P11][P0][0] = 1.e-31;
  OVv->omega[P11][P0][1] = 1.e-31;
  OVv->omega[P11][P0][2] = 1.e-31;
  OVv->omega[P11][P0][3] = 1.e-31;

  OVv->omega[P11][P1][0] = 1.e-31;
  OVv->omega[P11][P1][1] = 1.e-31;
  OVv->omega[P11][P1][2] = 1.e-31;
  OVv->omega[P11][P1][3] = 1.e-31;

  OVv->omega[P11][P2][0] = 1.e-31;
  OVv->omega[P11][P2][1] = 1.e-31;
  OVv->omega[P11][P2][2] = 1.e-31;
  OVv->omega[P11][P2][3] = 1.e-31;
*/

  OVv->dE[P2][S0] = 1213.8;   OVv->A[P2][S0] = 2.16e-2;
  OVv->dE[P1][S0] = 1218.3;   OVv->A[P1][S0] = 2.25e+3;
  OVv->dE[P0][S0] = 1220.4;   OVv->A[P0][S0] = 0.00000;
  OVv->dE[P11][S0] = 629.7;   OVv->A[P11][S0] = 2.80e+9;
  OVv->dE[P1][P0] = 7.35e+5;  OVv->A[P1][P0] = 5.81e-5;
  OVv->dE[P2][P0] = 2.26e+5;  OVv->A[P2][P0] = 0.00000;
  OVv->dE[P2][P1] = 3.26e+5;  OVv->A[P2][P1] = 3.55e-4;
/*  OVv->dE[P11][P0] = 1.e+20;  OVv->A[P11][P0] = 0.00000;
  OVv->dE[P11][P1] = 1.e+20;  OVv->A[P11][P1] = 0.00000;
  OVv->dE[P11][P2] = 1.e+20;  OVv->A[P11][P2] = 0.00000;
*/                                                                                                                  
  Symmetrize_Coeff (OVv);
}

/* ********************************************************************* */
void NeI_INIT (Ion *NeIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int lev1, lev2;
                                                                                                    
  NeIv->N = 8.e-5;    /* abundance of Ne with respect to H (Raymond) */  
  NeIv->nlev = 2;
 
  lev1 = 0;
  lev2 = 1;
 
  NeIv->nTom  = 4;
  NeIv->Tom[0] = 0.5e4;
  NeIv->Tom[1] = 1.0e4;
  NeIv->Tom[2] = 1.5e4;
  NeIv->Tom[3] = 2.0e4;

  if (NeIv->dE == NULL){
    NeIv->dE    = ARRAY_2D(NeIv->nlev, NeIv->nlev, double);
    NeIv->A     = ARRAY_2D(NeIv->nlev, NeIv->nlev, double);
    NeIv->omega = ARRAY_3D(NeIv->nlev, NeIv->nlev, NeIv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < NeIv->nlev ; i++) {
  for (j = 0 ; j < NeIv->nlev ; j++) {
    NeIv->dE[i][j] = -1.0;
    NeIv->A[i][j]  = -1.0;
    for (n = 0; n < NeIv->nTom; n++){
      NeIv->omega[i][j][n]  = -1.0;
    }
  }}

   
  NeIv->wght[lev1] = 1.; /*    2 J  +  1     */
  NeIv->wght[lev2] = 3.; /*         */

  NeIv->omega[lev1][lev2][0] = 0.0;
  NeIv->omega[lev1][lev2][1] = 0.0;
  NeIv->omega[lev1][lev2][2] = 0.0;
  NeIv->omega[lev1][lev2][3] = 0.0;


  NeIv->dE[lev1][lev2] = 743.7;  NeIv->A[lev1][lev2]  = 4.76e+7;  

  Symmetrize_Coeff (NeIv);
}

/* ********************************************************************* */
void NeII_INIT (Ion *NeIIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int _2P12, _2P32;
                                                                                                    
  NeIIv->N = 8.e-5;    /* abundance of Ne with respect to H (Raymond) */  
  NeIIv->nlev = 2;
 
  _2P32 = 0;
  _2P12 = 1;
 
  NeIIv->nTom  = 4;
  NeIIv->Tom[0] = 0.5e4;
  NeIIv->Tom[1] = 1.0e4;
  NeIIv->Tom[2] = 1.5e4;
  NeIIv->Tom[3] = 2.0e4;

  if (NeIIv->dE == NULL){
    NeIIv->dE    = ARRAY_2D(NeIIv->nlev, NeIIv->nlev, double);
    NeIIv->A     = ARRAY_2D(NeIIv->nlev, NeIIv->nlev, double);
    NeIIv->omega = ARRAY_3D(NeIIv->nlev, NeIIv->nlev, NeIIv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < NeIIv->nlev ; i++) {
  for (j = 0 ; j < NeIIv->nlev ; j++) {
    NeIIv->dE[i][j] = -1.0;
    NeIIv->A[i][j]  = -1.0;
    for (n = 0; n < NeIIv->nTom; n++){
      NeIIv->omega[i][j][n]  = -1.0;
    }
  }}

   
  NeIIv->wght[_2P12] = 2.; /*    2 J  +  1     */
  NeIIv->wght[_2P32] = 4.; /*         */

  NeIIv->omega[_2P12][_2P32][0] = 2.96e-1;
  NeIIv->omega[_2P12][_2P32][1] = 3.03e-1;
  NeIIv->omega[_2P12][_2P32][2] = 3.10e-1;
  NeIIv->omega[_2P12][_2P32][3] = 3.17e-1;


  NeIIv->dE[_2P12][_2P32] = 1.28e+5;  NeIIv->A[_2P12][_2P32]  = 8.55e-3;  

  Symmetrize_Coeff (NeIIv);
}

/* ********************************************************************* */
void NeIII_INIT (Ion *NeIIIv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int _3P2, _3P1, _3P0, _1D2, _1S0;
                                                                                                    
  NeIIIv->N = 8.e-5;    /* abundance of Ne with respect to H (Raymond) */  
  NeIIIv->nlev = 5;
 
  _3P2 = 0;
  _3P1 = 1;
  _3P0 = 2;
  _1D2 = 3;
  _1S0 = 4;
 
  NeIIIv->nTom  = 4;
  NeIIIv->Tom[0] = 0.5e4;
  NeIIIv->Tom[1] = 1.0e4;
  NeIIIv->Tom[2] = 1.5e4;
  NeIIIv->Tom[3] = 2.0e4;

  if (NeIIIv->dE == NULL){
    NeIIIv->dE    = ARRAY_2D(NeIIIv->nlev, NeIIIv->nlev, double);
    NeIIIv->A     = ARRAY_2D(NeIIIv->nlev, NeIIIv->nlev, double);
    NeIIIv->omega = ARRAY_3D(NeIIIv->nlev, NeIIIv->nlev, NeIIIv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < NeIIIv->nlev ; i++) {
  for (j = 0 ; j < NeIIIv->nlev ; j++) {
    NeIIIv->dE[i][j] = -1.0;
    NeIIIv->A[i][j]  = -1.0;
    for (n = 0; n < NeIIIv->nTom; n++){
      NeIIIv->omega[i][j][n]  = -1.0;
    }
  }}

   
  NeIIIv->wght[_3P2] = 5.; /*    2 J  +  1     */
  NeIIIv->wght[_3P1] = 3.; /*         */
  NeIIIv->wght[_3P0] = 1.; /*         */
  NeIIIv->wght[_1D2] = 5.; /*         */
  NeIIIv->wght[_1S0] = 1.; /*         */

  NeIIIv->omega[_1D2][_3P0][0] = 1.63 * (1.0) / (3.0*3.0);
  NeIIIv->omega[_1D2][_3P0][1] = 1.65 * (1.0) / (3.0*3.0);
  NeIIIv->omega[_1D2][_3P0][2] = 1.65 * (1.0) / (3.0*3.0);
  NeIIIv->omega[_1D2][_3P0][3] = 1.64 * (1.0) / (3.0*3.0);

  NeIIIv->omega[_1D2][_3P1][0] = 1.63 * (3.0) / (3.0*3.0);
  NeIIIv->omega[_1D2][_3P1][1] = 1.65 * (3.0) / (3.0*3.0);
  NeIIIv->omega[_1D2][_3P1][2] = 1.65 * (3.0) / (3.0*3.0);
  NeIIIv->omega[_1D2][_3P1][3] = 1.64 * (3.0) / (3.0*3.0);

  NeIIIv->omega[_1D2][_3P2][0] = 1.63 * (5.0) / (3.0*3.0);
  NeIIIv->omega[_1D2][_3P2][1] = 1.65 * (5.0) / (3.0*3.0);
  NeIIIv->omega[_1D2][_3P2][2] = 1.65 * (5.0) / (3.0*3.0);
  NeIIIv->omega[_1D2][_3P2][3] = 1.64 * (5.0) / (3.0*3.0);

  NeIIIv->omega[_1S0][_3P1][0] = 1.51e-1 * (3.0) / (3.0*3.0);
  NeIIIv->omega[_1S0][_3P1][1] = 1.69e-1 * (3.0) / (3.0*3.0);
  NeIIIv->omega[_1S0][_3P1][2] = 1.75e-1 * (3.0) / (3.0*3.0);
  NeIIIv->omega[_1S0][_3P1][3] = 1.79e-1 * (3.0) / (3.0*3.0);

  NeIIIv->omega[_1S0][_3P2][0] = 1.51e-1 * (5.0) / (3.0*3.0);
  NeIIIv->omega[_1S0][_3P2][1] = 1.69e-1 * (5.0) / (3.0*3.0);
  NeIIIv->omega[_1S0][_3P2][2] = 1.75e-1 * (5.0) / (3.0*3.0);
  NeIIIv->omega[_1S0][_3P2][3] = 1.79e-1 * (5.0) / (3.0*3.0);

  NeIIIv->omega[_1S0][_1D2][0] = 2.00e-1;
  NeIIIv->omega[_1S0][_1D2][1] = 2.26e-1;
  NeIIIv->omega[_1S0][_1D2][2] = 2.43e-1;
  NeIIIv->omega[_1S0][_1D2][3] = 2.60e-1;

  NeIIIv->omega[_3P0][_3P1][0] = 3.31e-1;
  NeIIIv->omega[_3P0][_3P1][1] = 3.50e-1;
  NeIIIv->omega[_3P0][_3P1][2] = 3.51e-1;
  NeIIIv->omega[_3P0][_3P1][3] = 3.50e-1;

  NeIIIv->omega[_3P0][_3P2][0] = 3.00e-1;
  NeIIIv->omega[_3P0][_3P2][1] = 3.07e-1;
  NeIIIv->omega[_3P0][_3P2][2] = 3.03e-1;
  NeIIIv->omega[_3P0][_3P2][3] = 2.98e-1;

  NeIIIv->omega[_3P1][_3P2][0] = 1.09;
  NeIIIv->omega[_3P1][_3P2][1] = 1.65;
  NeIIIv->omega[_3P1][_3P2][2] = 1.65;
  NeIIIv->omega[_3P1][_3P2][3] = 1.64;

  NeIIIv->dE[_1D2][_3P0] = 4012.8;  NeIIIv->A[_1D2][_3P0]  = 8.51e-6;  
  NeIIIv->dE[_1D2][_3P1] = 3967.5;  NeIIIv->A[_1D2][_3P1]  = 5.42e-2;  
  NeIIIv->dE[_1D2][_3P2] = 3868.8;  NeIIIv->A[_1D2][_3P2]  = 1.71e-1;  
  NeIIIv->dE[_1S0][_3P1] = 1814.6;  NeIIIv->A[_1S0][_3P1]  = 2.00;  
  NeIIIv->dE[_1S0][_3P2] = 1793.7;  NeIIIv->A[_1S0][_3P2]  = 3.94e-3;  
  NeIIIv->dE[_1S0][_1D2] = 3342.5;  NeIIIv->A[_1S0][_1D2]  = 2.71;  
  NeIIIv->dE[_3P0][_3P1] = 3.60e+5; NeIIIv->A[_3P0][_3P1]  = 1.15e-3;  
  NeIIIv->dE[_3P0][_3P2] = 1.07e+5; NeIIIv->A[_3P0][_3P2]  = 2.18e-8;  
  NeIIIv->dE[_3P1][_3P2] = 1.56e+5; NeIIIv->A[_3P1][_3P2]  = 5.97e-3;  
/*  NeIIIv->dE[_1S0][_3P0] = 1827.6; NeIIIv->A[_1S0][_3P0]  = 0.00000;  fake */  


  Symmetrize_Coeff (NeIIIv);
}

/* ********************************************************************* */
void NeIV_INIT (Ion *NeIVv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int S32, D52, D32, P12, P32;
                                                                                                    
  NeIVv->N = 8.e-5;    /* abundance of Ne with respect to H (Raymond) */  
  NeIVv->nlev = 5;
 
  S32 = 0;
  D52 = 1;
  D32 = 2;
  P12 = 3;
  P32 = 4;
 
  NeIVv->nTom  = 4;
  NeIVv->Tom[0] = 0.5e4;
  NeIVv->Tom[1] = 1.0e4;
  NeIVv->Tom[2] = 1.5e4;
  NeIVv->Tom[3] = 2.0e4;

  if (NeIVv->dE == NULL){
    NeIVv->dE    = ARRAY_2D(NeIVv->nlev, NeIVv->nlev, double);
    NeIVv->A     = ARRAY_2D(NeIVv->nlev, NeIVv->nlev, double);
    NeIVv->omega = ARRAY_3D(NeIVv->nlev, NeIVv->nlev, NeIVv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < NeIVv->nlev ; i++) {
  for (j = 0 ; j < NeIVv->nlev ; j++) {
    NeIVv->dE[i][j] = -1.0;
    NeIVv->A[i][j]  = -1.0;
    for (n = 0; n < NeIVv->nTom; n++){
      NeIVv->omega[i][j][n]  = -1.0;
    }
  }}

   
  NeIVv->wght[S32] = 4.; /*    2 J  +  1     */
  NeIVv->wght[D52] = 6.; /*         */
  NeIVv->wght[D32] = 4.; /*         */
  NeIVv->wght[P12] = 2.; /*         */
  NeIVv->wght[P32] = 4.; /*         */

  NeIVv->omega[D52][S32][0] = 8.45e-1;
  NeIVv->omega[D52][S32][1] = 8.43e-1;
  NeIVv->omega[D52][S32][2] = 8.32e-1;
  NeIVv->omega[D52][S32][3] = 8.24e-1;
  
  NeIVv->omega[D32][S32][0] = 5.63e-1;
  NeIVv->omega[D32][S32][1] = 5.59e-1;
  NeIVv->omega[D32][S32][2] = 5.55e-1;
  NeIVv->omega[D32][S32][3] = 5.50e-1;

  NeIVv->omega[P32][S32][0] = 3.07e-1;
  NeIVv->omega[P32][S32][1] = 3.13e-1;
  NeIVv->omega[P32][S32][2] = 3.12e-1;
  NeIVv->omega[P32][S32][3] = 3.09e-1;

  NeIVv->omega[P12][S32][0] = 1.53e-1;
  NeIVv->omega[P12][S32][1] = 1.56e-1;
  NeIVv->omega[P12][S32][2] = 1.56e-1;
  NeIVv->omega[P12][S32][3] = 1.55e-1;

  NeIVv->omega[D52][D32][0] = 1.37;
  NeIVv->omega[D52][D32][1] = 1.36;
  NeIVv->omega[D52][D32][2] = 1.35;
  NeIVv->omega[D52][D32][3] = 1.33;

  NeIVv->omega[P32][P12][0] = 3.17e-1;
  NeIVv->omega[P32][P12][1] = 3.43e-1;
  NeIVv->omega[P32][P12][2] = 3.58e-1;
  NeIVv->omega[P32][P12][3] = 3.70e-1;

  NeIVv->omega[P32][D52][0] = 8.56e-1;
  NeIVv->omega[P32][D52][1] = 9.00e-1;
  NeIVv->omega[P32][D52][2] = 9.08e-1;
  NeIVv->omega[P32][D52][3] = 9.09e-1;

  NeIVv->omega[P32][D32][0] = 4.73e-1;
  NeIVv->omega[P32][D32][1] = 5.09e-1;
  NeIVv->omega[P32][D32][2] = 5.15e-1;
  NeIVv->omega[P32][D32][3] = 5.16e-1;

  NeIVv->omega[P12][D52][0] = 3.40e-1;
  NeIVv->omega[P12][D52][1] = 3.68e-1;
  NeIVv->omega[P12][D52][2] = 3.73e-1;
  NeIVv->omega[P12][D52][3] = 3.74e-1;

  NeIVv->omega[P12][D32][0] = 3.24e-1;
  NeIVv->omega[P12][D32][1] = 3.36e-1;
  NeIVv->omega[P12][D32][2] = 3.39e-1;
  NeIVv->omega[P12][D32][3] = 3.39e-1;


  NeIVv->dE[D52][S32] = 2420.9;   NeIVv->A[D52][S32]  = 4.58e-4;
  NeIVv->dE[D32][S32] = 2418.2;   NeIVv->A[D32][S32]  = 5.77e-3;  
  NeIVv->dE[P32][S32] = 1601.5;   NeIVv->A[P32][S32]  = 1.27;  
  NeIVv->dE[P12][S32] = 1601.7;   NeIVv->A[P12][S32]  = 5.21e-1;  
  NeIVv->dE[D52][D32] = 2.237e+6; NeIVv->A[D52][D32]  = 1.48e-6;  
  NeIVv->dE[P32][P12] = 1.56e+7;  NeIVv->A[P32][P12]  = 2.82e-9;
  NeIVv->dE[P32][D52] = 4714.3;   NeIVv->A[P32][D52]  = 3.88e-1;  
  NeIVv->dE[P32][D32] = 4724.2;   NeIVv->A[P32][D32]  = 4.37e-1;  
  NeIVv->dE[P12][D52] = 4717.0;   NeIVv->A[P12][D52]  = 1.15e-2;  
  NeIVv->dE[P12][D32] = 4725.6;   NeIVv->A[P12][D32]  = 3.93e-1;


  Symmetrize_Coeff (NeIVv);
}

/* ********************************************************************* */
void NeV_INIT (Ion *NeVv)
/*
 *  
 *********************************************************************** */
{
  int i, j, n;
  int _3P2, _3P1, _3P0, _1D2, _1S0;

  NeVv->N = 8.e-5;    /* abundance of Ne with respect to H (Raymond) */  
  NeVv->nlev = 5;
 
  _3P0 = 0;
  _3P1 = 1;
  _3P2 = 2;
  _1D2 = 3;
  _1S0 = 4;
 
  NeVv->nTom  = 4;
  NeVv->Tom[0] = 0.5e4;
  NeVv->Tom[1] = 1.0e4;
  NeVv->Tom[2] = 1.5e4;
  NeVv->Tom[3] = 2.0e4;

  if (NeVv->dE == NULL){
    NeVv->dE    = ARRAY_2D(NeVv->nlev, NeVv->nlev, double);
    NeVv->A     = ARRAY_2D(NeVv->nlev, NeVv->nlev, double);
    NeVv->omega = ARRAY_3D(NeVv->nlev, NeVv->nlev, NeVv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < NeVv->nlev ; i++) {
  for (j = 0 ; j < NeVv->nlev ; j++) {
    NeVv->dE[i][j] = -1.0;
    NeVv->A[i][j]  = -1.0;
    for (n = 0; n < NeVv->nTom; n++){
      NeVv->omega[i][j][n]  = -1.0;
    }
  }}

   
  NeVv->wght[_3P0] = 1.; /*    2 J  +  1     */
  NeVv->wght[_3P1] = 3.; /*         */
  NeVv->wght[_3P2] = 5.; /*         */
  NeVv->wght[_1D2] = 5.; /*         */
  NeVv->wght[_1S0] = 1.; /*         */

  NeVv->omega[_1D2][_3P0][0] = 2.13 * (1.0) / (3.0*3.0);
  NeVv->omega[_1D2][_3P0][1] = 2.09 * (1.0) / (3.0*3.0);
  NeVv->omega[_1D2][_3P0][2] = 2.11 * (1.0) / (3.0*3.0);
  NeVv->omega[_1D2][_3P0][3] = 2.14 * (1.0) / (3.0*3.0);

  NeVv->omega[_1D2][_3P1][0] = 2.13 * (3.0) / (3.0*3.0);
  NeVv->omega[_1D2][_3P1][1] = 2.09 * (3.0) / (3.0*3.0);
  NeVv->omega[_1D2][_3P1][2] = 2.11 * (3.0) / (3.0*3.0);
  NeVv->omega[_1D2][_3P1][3] = 2.14 * (3.0) / (3.0*3.0);

  NeVv->omega[_1D2][_3P2][0] = 2.13 * (5.0) / (3.0*3.0);
  NeVv->omega[_1D2][_3P2][1] = 2.09 * (5.0) / (3.0*3.0);
  NeVv->omega[_1D2][_3P2][2] = 2.11 * (5.0) / (3.0*3.0);
  NeVv->omega[_1D2][_3P2][3] = 2.14 * (5.0) / (3.0*3.0);

  NeVv->omega[_1S0][_3P1][0] = 2.54e-1 * (3.0) / (3.0*3.0);
  NeVv->omega[_1S0][_3P1][1] = 2.46e-1 * (3.0) / (3.0*3.0);
  NeVv->omega[_1S0][_3P1][2] = 2.49e-1 * (3.0) / (3.0*3.0);
  NeVv->omega[_1S0][_3P1][3] = 2.51e-1 * (3.0) / (3.0*3.0);

  NeVv->omega[_1S0][_3P2][0] = 2.54e-1 * (5.0) / (3.0*3.0);
  NeVv->omega[_1S0][_3P2][1] = 2.46e-1 * (5.0) / (3.0*3.0);
  NeVv->omega[_1S0][_3P2][2] = 2.49e-1 * (5.0) / (3.0*3.0);
  NeVv->omega[_1S0][_3P2][3] = 2.51e-1 * (5.0) / (3.0*3.0);

  NeVv->omega[_1S0][_1D2][0] = 6.63e-1;
  NeVv->omega[_1S0][_1D2][1] = 5.77e-1;
  NeVv->omega[_1S0][_1D2][2] = 6.10e-1;
  NeVv->omega[_1S0][_1D2][3] = 6.49e-1;

  NeVv->omega[_3P1][_3P0][0] = 1.68;
  NeVv->omega[_3P1][_3P0][1] = 1.41;
  NeVv->omega[_3P1][_3P0][2] = 1.19;
  NeVv->omega[_3P1][_3P0][3] = 1.10;

  NeVv->omega[_3P2][_3P0][0] = 2.44;
  NeVv->omega[_3P2][_3P0][1] = 1.81;
  NeVv->omega[_3P2][_3P0][2] = 1.42;
  NeVv->omega[_3P2][_3P0][3] = 1.26;

  NeVv->omega[_3P2][_3P1][0] = 7.59;
  NeVv->omega[_3P2][_3P1][1] = 5.82;
  NeVv->omega[_3P2][_3P1][2] = 4.68;
  NeVv->omega[_3P2][_3P1][3] = 4.20;


  NeVv->dE[_1D2][_3P0] = 3301.3;  NeVv->A[_1D2][_3P0]  = 2.37e-5;
  NeVv->dE[_1D2][_3P1] = 3345.8;  NeVv->A[_1D2][_3P1]  = 1.31e-1;
  NeVv->dE[_1D2][_3P2] = 3425.9;  NeVv->A[_1D2][_3P2]  = 3.65e-1;
  NeVv->dE[_1S0][_3P1] = 1574.8;  NeVv->A[_1S0][_3P1]  = 4.21;  
  NeVv->dE[_1S0][_3P2] = 1592.3;  NeVv->A[_1S0][_3P2]  = 6.69e-3;  
  NeVv->dE[_1S0][_1D2] = 2972.8;  NeVv->A[_1S0][_1D2]  = 2.85;  
  NeVv->dE[_3P1][_3P0] = 2.428e+5;NeVv->A[_3P1][_3P0]  = 1.28e-3;  
  NeVv->dE[_3P2][_3P0] = 90082.;  NeVv->A[_3P2][_3P0]  = 5.08e-9;  
  NeVv->dE[_3P2][_3P1] = 1.432e+5;NeVv->A[_3P2][_3P1]  = 4.59e-3;  
/*  NeVv->dE[_1S0][_3P0] = 1.6e+16; NeVv->A[_1S0][_3P0]  = 0.00000;  fake */  


  Symmetrize_Coeff (NeVv);
}

/* ********************************************************************* */
void SI_INIT (Ion *SIv)
/*
 * Data from Tayal 2004, ApJSS   
 *********************************************************************** */

{
  int i, j, n;
  int P2, P1, P0, D2, S0;

  SIv->N    = 4.e-5;              /* abundance of S with respect to H (Raymond) */
  SIv->nlev = 5;

  P2 = 0;
  P1 = 1;
  P0 = 2;
  D2 = 3;
  S0 = 4;
   
/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER, 
    O GET PROPER ORDER LATER  */

  SIv->nTom  = 4;
  SIv->Tom[0] = 1.995e3;
  SIv->Tom[1] = 1.0e4;
  SIv->Tom[2] = 1.0e5;
  SIv->Tom[3] = 3.981e5;

  if (SIv->dE == NULL){
    SIv->dE    = ARRAY_2D(SIv->nlev, SIv->nlev, double);
    SIv->A     = ARRAY_2D(SIv->nlev, SIv->nlev, double);
    SIv->omega = ARRAY_3D(SIv->nlev, SIv->nlev, SIv->nTom, double);
  }

  for (i = 0 ; i < SIv->nlev ; i++) {
  for (j = 0 ; j < SIv->nlev ; j++) {
    SIv->dE[i][j] = -1.0;
    SIv->A[i][j]  = -1.0;
    for (n = 0; n < SIv->nTom; n++){
      SIv->omega[i][j][n]  = -1.0;
    }
  }}
  
  SIv->wght[P2] = 5.; /*    3P_2       */
  SIv->wght[P1] = 3.; /*    3P_1    */
  SIv->wght[P0] = 1.; /*    3P_0     */
  SIv->wght[D2] = 5.; /*    1D_2      */
  SIv->wght[S0] = 1.; /*    1S_0      */

/* ---------------------------------------------------
     define collision strengths as function of
     temperature.
   --------------------------------------------------- */
  
 
  SIv->omega[P0][P1][0] = 0.056;
  SIv->omega[P0][P1][1] = 0.283;
  SIv->omega[P0][P1][2] = 1.862;
  SIv->omega[P0][P1][3] = 1.897;
 
  SIv->omega[P0][P2][0] = 0.011;
  SIv->omega[P0][P2][1] = 0.054;
  SIv->omega[P0][P2][2] = 0.415;
  SIv->omega[P0][P2][3] = 0.507;

  SIv->omega[D2][P0][0] = 0.058;
  SIv->omega[D2][P0][1] = 0.273;
  SIv->omega[D2][P0][2] = 1.609;
  SIv->omega[D2][P0][3] = 1.454;

  SIv->omega[S0][P0][0] = 0.007;
  SIv->omega[S0][P0][1] = 0.031;
  SIv->omega[S0][P0][2] = 0.195;
  SIv->omega[S0][P0][3] = 0.167;
/*
  SIv->omega[S0][D2][0] = 1.0e-31;
  SIv->omega[S0][D2][1] = 1.0e-31;
  SIv->omega[S0][D2][2] = 1.0e-31;
  SIv->omega[S0][D2][3] = 1.0e-31;

  SIv->omega[S0][P1][0] = 1.0e-31;
  SIv->omega[S0][P1][1] = 1.0e-31;
  SIv->omega[S0][P1][2] = 1.0e-31;
  SIv->omega[S0][P1][3] = 1.0e-31;

  SIv->omega[S0][P2][0] = 1.0e-31;
  SIv->omega[S0][P2][1] = 1.0e-31;
  SIv->omega[S0][P2][2] = 1.0e-31;
  SIv->omega[S0][P2][3] = 1.0e-31;

  SIv->omega[D2][P1][0] = 1.0e-31;
  SIv->omega[D2][P1][1] = 1.0e-31;
  SIv->omega[D2][P1][2] = 1.0e-31;
  SIv->omega[D2][P1][3] = 1.0e-31;

  SIv->omega[D2][P2][0] = 1.0e-31;
  SIv->omega[D2][P2][1] = 1.0e-31;
  SIv->omega[D2][P2][2] = 1.0e-31;
  SIv->omega[D2][P2][3] = 1.0e-31;

  SIv->omega[P1][P2][0] = 1.0e-31;
  SIv->omega[P1][P2][1] = 1.0e-31;
  SIv->omega[P1][P2][2] = 1.0e-31;
  SIv->omega[P1][P2][3] = 1.0e-31;
*/

  SIv->dE[P0][P1] = 563111. ;  SIv->A[P0][P1] = 3.02e-4; /* dE [ HIGHER ][LOWER ]  */
  SIv->dE[P0][P2] = 174325.;   SIv->A[P0][P2] = 7.1e-8;
  SIv->dE[D2][P0] = 11540.72;  SIv->A[D2][P0] = 5.0e-6;
  SIv->dE[S0][P0] = 4628.28;   SIv->A[S0][P0] = 3.5e-1;
/*  SIv->dE[S0][D2] = 1.0e+16;   SIv->A[S0][D2] = 0.00000;
  SIv->dE[S0][P1] = 1.0e+16;   SIv->A[S0][P1] = 0.00000;
  SIv->dE[S0][P2] = 1.0e+16;   SIv->A[S0][P2] = 0.00000;
  SIv->dE[D2][P1] = 1.0e+16;   SIv->A[D2][P1] = 0.00000;
  SIv->dE[D2][P2] = 1.0e+16;   SIv->A[D2][P2] = 0.00000;
  SIv->dE[P1][P2] = 1.0e+16;   SIv->A[P1][P2] = 0.00000;
*/
  Symmetrize_Coeff (SIv);
}



/* ********************************************************************* */
void SII_INIT (Ion *SIIv)
/*
 *
 *
 *      P32
 *
 *      P12
 *
 *      D52
 *
 *      D32
 *
 *      S
 *
 *********************************************************************** */
{
  int i, j, n;
  int S, D52, D32, P32, P12;

  /*  T-dep CS data from MAPPINGS. Old Pradhan data commented.   */

  SIIv->isMAP = 1;

  SIIv->N    = 4.e-5;              /* abundance of S with respect to H (Raymond) */
  SIIv->nlev = 5;

  S   = 0;
  D32 = 1;
  D52 = 2;
  P12 = 3;
  P32 = 4;
   
/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER, 
    O GET PROPER ORDER LATER  */

  SIIv->nTom  = 2;
  SIIv->Tom[0] = 0.5e4;
  SIIv->Tom[1] = 1.0e4;

  if (SIIv->dE == NULL){
    SIIv->dE    = ARRAY_2D(SIIv->nlev, SIIv->nlev, double);
    SIIv->A     = ARRAY_2D(SIIv->nlev, SIIv->nlev, double);
    SIIv->omega = ARRAY_3D(SIIv->nlev, SIIv->nlev, SIIv->nTom, double);
  }

  for (i = 0 ; i < SIIv->nlev ; i++) {
  for (j = 0 ; j < SIIv->nlev ; j++) {
    SIIv->dE[i][j] = -1.0;
    SIIv->A[i][j]  = -1.0;
    for (n = 0; n < SIIv->nTom; n++){
      SIIv->omega[i][j][n]  = 0.0;
    }
  }}
  
  SIIv->wght[S]   = 4.; /*    4S_3/2       */
  SIIv->wght[D52] = 6.; /*    2D_5/2    */
  SIIv->wght[D32] = 4.; /*    2D_3/2     */
  SIIv->wght[P32] = 4.; /*    4P_3/2      */
  SIIv->wght[P12] = 2.; /*    4P_1/2      */

/* ---------------------------------------------------
     define collision strengths as function of
     temperature.
   --------------------------------------------------- */
  
/*
  SIIv->omega[S][D52][0] = 4.90;
  SIIv->omega[S][D52][1] = 4.66;
  SIIv->omega[S][D52][2] = 4.44;
  SIIv->omega[S][D52][3] = 4.26;

  SIIv->omega[S][D32][0] = 3.27;
  SIIv->omega[S][D32][1] = 3.11;
  SIIv->omega[S][D32][2] = 2.97;
  SIIv->omega[S][D32][3] = 2.84;

  SIIv->omega[P32][S][0] = 1.67;
  SIIv->omega[P32][S][1] = 2.07;
  SIIv->omega[P32][S][2] = 1.98;
  SIIv->omega[P32][S][3] = 2.07;

  SIIv->omega[P12][S][0] = 0.831;
  SIIv->omega[P12][S][1] = 0.897;
  SIIv->omega[P12][S][2] = 0.987;
  SIIv->omega[P12][S][3] = 1.03;

  SIIv->omega[D52][D32][0] = 7.9;
  SIIv->omega[D52][D32][1] = 7.46;
  SIIv->omega[D52][D32][2] = 7.11;
  SIIv->omega[D52][D32][3] = 8.65;

  SIIv->omega[P32][P12][0] = 2.02;
  SIIv->omega[P32][P12][1] = 2.54;
  SIIv->omega[P32][P12][2] = 2.13;
  SIIv->omega[P32][P12][3] = 2.22;

  SIIv->omega[P32][D52][0] = 5.93;
  SIIv->omega[P32][D52][1] = 4.77;
  SIIv->omega[P32][D52][2] = 4.75;
  SIIv->omega[P32][D52][3] = 4.68;

  SIIv->omega[P32][D32][0] = 3.41;
  SIIv->omega[P32][D32][1] = 2.74;
  SIIv->omega[P32][D32][2] = 2.74;
  SIIv->omega[P32][D32][3] = 2.71;

  SIIv->omega[P12][D52][0] = 2.57;
  SIIv->omega[P12][D52][1] = 1.99;
  SIIv->omega[P12][D52][2] = 1.99;
  SIIv->omega[P12][D52][3] = 1.97;

  SIIv->omega[P12][D32][0] = 2.20;
  SIIv->omega[P12][D32][1] = 1.76;
  SIIv->omega[P12][D32][2] = 1.76;
  SIIv->omega[P12][D32][3] = 1.73;
*/

  SIIv->omega[S][D52][0] = 4.1900;
  SIIv->omega[S][D52][1] = -.0800;

  SIIv->omega[S][D32][0] = 2.7900;
  SIIv->omega[S][D32][1] = -.0800;

  SIIv->omega[P32][S][0] = 1.5200;
  SIIv->omega[P32][S][1] = 0.0000;

  SIIv->omega[P12][S][0] = 0.7590;
  SIIv->omega[P12][S][1] = 0.0000;

  SIIv->omega[D52][D32][0] = 7.5900;
  SIIv->omega[D52][D32][1] = -.1400;

  SIIv->omega[P32][P12][0] = 2.3800;
  SIIv->omega[P32][P12][1] = 0.0000;

  SIIv->omega[P32][D52][0] = 4.7900;
  SIIv->omega[P32][D52][1] = 0.0000;

  SIIv->omega[P32][D32][0] = 3.3800;
  SIIv->omega[P32][D32][1] = 0.0000;

  SIIv->omega[P12][D52][0] = 2.5600;
  SIIv->omega[P12][D52][1] = 0.0000;

  SIIv->omega[P12][D32][0] = 1.5200;
  SIIv->omega[P12][D32][1] = 0.0000;


  SIIv->dE[P32][P12] = 2.14e6 ;  SIIv->A[P32][P12] = 1.03e-6; /* dE [ HIGHER ][LOWER ]  */
  SIIv->dE[P32][D52] = 10320.5;  SIIv->A[D52][P32] = 0.18;
  SIIv->dE[P32][D32] = 10286.7;  SIIv->A[D32][P32] = 0.13;
  SIIv->dE[P12][D52] = 10370.5;  SIIv->A[D52][P12] = 7.8e-2;
  SIIv->dE[P12][D32] = 10336.4;  SIIv->A[P12][D32] = 0.16;
  SIIv->dE[P32][S]   = 4068.6;   SIIv->A[P32][S]   = 0.22;
  SIIv->dE[P12][S]   = 4076.4;   SIIv->A[P12][S]   = 9.1e-2;
  SIIv->dE[D52][D32] = 3.145e6;  SIIv->A[D52][D32] = 3.3e-7;
  SIIv->dE[S][D52]   = 6716.4;   SIIv->A[S][D52]   = 2.6e-4;
  SIIv->dE[S][D32]   = 6730.8;   SIIv->A[S][D32]   = 8.8e-4;

  Symmetrize_Coeff (SIIv);
}

/* ********************************************************************* */
void SIII_INIT (Ion *SIIIv)
/*
 *
 *********************************************************************** */
{
  int i, j, n;
  int _3P2, _3P1, _3P0, _1D2, _1S0;
                                                                                                    
  SIIIv->N = 8.e-5;    /* abundance of Ne with respect to H (Raymond) */  
  SIIIv->nlev = 5;
 
  _3P0 = 0;
  _3P1 = 1;
  _3P2 = 2;
  _1D2 = 3;
  _1S0 = 4;
 
  SIIIv->nTom  = 4;
  SIIIv->Tom[0] = 0.5e4;
  SIIIv->Tom[1] = 1.0e4;
  SIIIv->Tom[2] = 1.5e4;
  SIIIv->Tom[3] = 2.0e4;

  if (SIIIv->dE == NULL){
    SIIIv->dE    = ARRAY_2D(SIIIv->nlev, SIIIv->nlev, double);
    SIIIv->A     = ARRAY_2D(SIIIv->nlev, SIIIv->nlev, double);
    SIIIv->omega = ARRAY_3D(SIIIv->nlev, SIIIv->nlev, SIIIv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < SIIIv->nlev ; i++) {
  for (j = 0 ; j < SIIIv->nlev ; j++) {
    SIIIv->dE[i][j] = -1.0;
    SIIIv->A[i][j]  = -1.0;
    for (n = 0; n < SIIIv->nTom; n++){
      SIIIv->omega[i][j][n]  = -1.0;
    }
  }}

   
  SIIIv->wght[_3P0] = 1.; /*    2 J  +  1     */
  SIIIv->wght[_3P1] = 3.; /*         */
  SIIIv->wght[_3P2] = 5.; /*         */
  SIIIv->wght[_1D2] = 5.; /*         */
  SIIIv->wght[_1S0] = 1.; /*         */

  SIIIv->omega[_1D2][_3P0][0] = 9.07 * (1.0) / (3.0*3.0);
  SIIIv->omega[_1D2][_3P0][1] = 8.39 * (1.0) / (3.0*3.0);
  SIIIv->omega[_1D2][_3P0][2] = 8.29 * (1.0) / (3.0*3.0);
  SIIIv->omega[_1D2][_3P0][3] = 8.20 * (1.0) / (3.0*3.0);

  SIIIv->omega[_1D2][_3P1][0] = 9.07 * (3.0) / (3.0*3.0);
  SIIIv->omega[_1D2][_3P1][1] = 8.39 * (3.0) / (3.0*3.0);
  SIIIv->omega[_1D2][_3P1][2] = 8.29 * (3.0) / (3.0*3.0);
  SIIIv->omega[_1D2][_3P1][3] = 8.20 * (3.0) / (3.0*3.0);

  SIIIv->omega[_1D2][_3P2][0] = 9.07 * (5.0) / (3.0*3.0);
  SIIIv->omega[_1D2][_3P2][1] = 8.39 * (5.0) / (3.0*3.0);
  SIIIv->omega[_1D2][_3P2][2] = 8.29 * (5.0) / (3.0*3.0);
  SIIIv->omega[_1D2][_3P2][3] = 8.20 * (5.0) / (3.0*3.0);

  SIIIv->omega[_1S0][_3P1][0] = 1.16 * (3.0) / (3.0*3.0);
  SIIIv->omega[_1S0][_3P1][1] = 1.19 * (3.0) / (3.0*3.0);
  SIIIv->omega[_1S0][_3P1][2] = 1.21 * (3.0) / (3.0*3.0);
  SIIIv->omega[_1S0][_3P1][3] = 1.24 * (3.0) / (3.0*3.0);

  SIIIv->omega[_1S0][_3P2][0] = 1.16 * (5.0) / (3.0*3.0);
  SIIIv->omega[_1S0][_3P2][1] = 1.19 * (5.0) / (3.0*3.0);
  SIIIv->omega[_1S0][_3P2][2] = 1.21 * (5.0) / (3.0*3.0);
  SIIIv->omega[_1S0][_3P2][3] = 1.24 * (5.0) / (3.0*3.0);

  SIIIv->omega[_1S0][_1D2][0] = 1.42;
  SIIIv->omega[_1S0][_1D2][1] = 1.88;
  SIIIv->omega[_1S0][_1D2][2] = 2.02;
  SIIIv->omega[_1S0][_1D2][3] = 2.08;

  SIIIv->omega[_3P1][_3P0][0] = 2.64;
  SIIIv->omega[_3P1][_3P0][1] = 2.59;
  SIIIv->omega[_3P1][_3P0][2] = 2.38;
  SIIIv->omega[_3P1][_3P0][3] = 2.20;

  SIIIv->omega[_3P2][_3P0][0] = 1.11;
  SIIIv->omega[_3P2][_3P0][1] = 1.15;
  SIIIv->omega[_3P2][_3P0][2] = 1.15;
  SIIIv->omega[_3P2][_3P0][3] = 1.14;

  SIIIv->omega[_3P2][_3P1][0] = 5.79;
  SIIIv->omega[_3P2][_3P1][1] = 5.81;
  SIIIv->omega[_3P2][_3P1][2] = 5.56;
  SIIIv->omega[_3P2][_3P1][3] = 5.32;


  SIIIv->dE[_1D2][_3P0] = 8833.9;  SIIIv->A[_1D2][_3P0]  = 5.82e-6;
  SIIIv->dE[_1D2][_3P1] = 9068.9;  SIIIv->A[_1D2][_3P1]  = 2.21e-2;
  SIIIv->dE[_1D2][_3P2] = 9531.0;  SIIIv->A[_1D2][_3P2]  = 5.76e-2;
  SIIIv->dE[_1S0][_3P1] = 3721.7;  SIIIv->A[_1S0][_3P1]  = 7.96e-1;  
  SIIIv->dE[_1S0][_3P2] = 3797.8;  SIIIv->A[_1S0][_3P2]  = 1.05e-2;  
  SIIIv->dE[_1S0][_1D2] = 6312.1;  SIIIv->A[_1S0][_1D2]  = 2.22;  
  SIIIv->dE[_3P1][_3P0] = 3.347e+5;SIIIv->A[_3P1][_3P0]  = 4.72e-4;  
  SIIIv->dE[_3P2][_3P0] = 1.20e+5; SIIIv->A[_3P2][_3P0]  = 4.61e-8;  
  SIIIv->dE[_3P2][_3P1] = 187129.; SIIIv->A[_3P2][_3P1]  = 2.07e-3;  
/*  SIIIv->dE[_1S0][_3P0] = 1.6e+16; SIIIv->A[_1S0][_3P0]  = 0.00000;  fake */  


  Symmetrize_Coeff (SIIIv);
}

/* ********************************************************************* */
void SIV_INIT (Ion *SIVv)
/*
 *
 *********************************************************************** */
{
  int i, j, n;
  int _2P12, _2P32, _4P12, _4P32, _4P52;
  
  
  /*  T-dep CS data from MAPPINGS. Old Pradhan data commented.   */

  SIVv->isMAP = 1;
  
  SIVv->N    = 4.e-5;              /* abundance of S with respect to H (Raymond) */
  SIVv->nlev = 5;

  _2P12 = 0;
  _2P32 = 1;
  _4P12 = 2;
  _4P32 = 3;
  _4P52 = 4;
   
/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER, 
    O GET PROPER ORDER LATER  */

  SIVv->nTom  = 2;
  SIVv->Tom[0] = 1.0e4;
  SIVv->Tom[1] = 1.5e4;

  if (SIVv->dE == NULL){
    SIVv->dE    = ARRAY_2D(SIVv->nlev, SIVv->nlev, double);
    SIVv->A     = ARRAY_2D(SIVv->nlev, SIVv->nlev, double);
    SIVv->omega = ARRAY_3D(SIVv->nlev, SIVv->nlev, SIVv->nTom, double);
  }

  for (i = 0 ; i < SIVv->nlev ; i++) {
  for (j = 0 ; j < SIVv->nlev ; j++) {
    SIVv->dE[i][j] = -1.0;
    SIVv->A[i][j]  = -1.0;
    for (n = 0; n < SIVv->nTom; n++){
      SIVv->omega[i][j][n]  = 0.0;
    }
  }}
  
  SIVv->wght[_2P12] = 2.; /*    2P_1/2       */
  SIVv->wght[_2P32] = 4.; /*    2P_3/2    */
  SIVv->wght[_4P12] = 2.; /*    4P_1/2     */
  SIVv->wght[_4P32] = 4.; /*    4P_3/2      */
  SIVv->wght[_4P52] = 6.; /*    4P_5/2      */

/* ---------------------------------------------------
     define collision strengths as function of
     temperature.
   --------------------------------------------------- */
  
/* 
  SIVv->omega[_2P32][_2P12][0] = 6.42;
  SIVv->omega[_2P32][_2P12][1] = 6.41;
  
  SIVv->omega[_4P12][_2P12][0] = 5.50e-1;
  SIVv->omega[_4P12][_2P12][1] = 4.80e-1;

  SIVv->omega[_4P12][_2P32][0] = 6.60e-1;
  SIVv->omega[_4P12][_2P32][1] = 6.30e-1;

  SIVv->omega[_4P32][_2P12][0] = 8.70e-1;
  SIVv->omega[_4P32][_2P12][1] = 8.30e-1;

  SIVv->omega[_4P32][_2P32][0] = 1.47;
  SIVv->omega[_4P32][_2P32][1] = 1.40;

  SIVv->omega[_4P32][_4P12][0] = 3.04;
  SIVv->omega[_4P32][_4P12][1] = 2.85;

  SIVv->omega[_4P52][_2P12][0] = 9.5e-1;
  SIVv->omega[_4P52][_2P12][1] = 9.1e-1;

  SIVv->omega[_4P52][_2P32][0] = 2.53;
  SIVv->omega[_4P52][_2P32][1] = 2.41;

  SIVv->omega[_4P52][_4P12][0] = 2.92;
  SIVv->omega[_4P52][_4P12][1] = 2.71;

  SIVv->omega[_4P52][_4P32][0] = 7.01;
  SIVv->omega[_4P52][_4P32][1] = 6.57;
*/

  SIVv->omega[_2P32][_2P12][0] = 6.4200;
  SIVv->omega[_2P32][_2P12][1] = -.0070;
  
  SIVv->omega[_4P12][_2P12][0] = 0.5100;
  SIVv->omega[_4P12][_2P12][1] = -.1600;

  SIVv->omega[_4P12][_2P32][0] = 6.60e-1;
  SIVv->omega[_4P12][_2P32][1] = -.1220;

  SIVv->omega[_4P32][_2P12][0] = 8.70e-1;
  SIVv->omega[_4P32][_2P12][1] = -.1330;

  SIVv->omega[_4P32][_2P32][0] = 1.4700;
  SIVv->omega[_4P32][_2P32][1] = -.1330;

  SIVv->omega[_4P32][_4P12][0] = 3.0400;
  SIVv->omega[_4P32][_4P12][1] = -.1740;

  SIVv->omega[_4P52][_2P12][0] = 9.5e-1;
  SIVv->omega[_4P52][_2P12][1] = -.1210;

  SIVv->omega[_4P52][_2P32][0] = 2.53;
  SIVv->omega[_4P52][_2P32][1] = -.1230;

  SIVv->omega[_4P52][_4P12][0] = 2.9200;
  SIVv->omega[_4P52][_4P12][1] = -.2130;

  SIVv->omega[_4P52][_4P32][0] = 7.0100;
  SIVv->omega[_4P52][_4P32][1] = -.1960;

  SIVv->dE[_2P32][_2P12] = 1.05e+5; SIVv->A[_2P32][_2P12] = 7.73e-3; /* dE [ HIGHER ][LOWER ]  */
  SIVv->dE[_4P12][_2P12] = 1404.9;  SIVv->A[_4P12][_2P12] = 5.50e+4;
  SIVv->dE[_4P12][_2P32] = 1423.9;  SIVv->A[_4P12][_2P32] = 3.39e+4;
  SIVv->dE[_4P32][_2P12] = 1398.1;  SIVv->A[_4P32][_2P12] = 1.40e+2;
  SIVv->dE[_4P32][_2P32] = 1017.0;  SIVv->A[_4P32][_2P32] = 1.95e+4;
  SIVv->dE[_4P32][_4P12] = 2.91e+5; SIVv->A[_4P32][_4P12] = 0.00000;
  SIVv->dE[_4P52][_2P12] = 1387.5;  SIVv->A[_4P52][_2P12] = 0.00000;
  SIVv->dE[_4P52][_2P32] = 1406.1;  SIVv->A[_4P52][_2P32] = 3.95e+4;
  SIVv->dE[_4P52][_4P12] = 1.12e+5; SIVv->A[_4P52][_4P12] = 0.00000;
  SIVv->dE[_4P52][_4P32] = 1.85e+5; SIVv->A[_4P52][_4P32] = 0.00000;

  Symmetrize_Coeff (SIVv);
}

/* ********************************************************************* */
void SV_INIT (Ion *SVv)
/*
 *
 *********************************************************************** */
{
  int i, j, n;
  int _1S0, _3P0, _3P1, _3P2, _1P1;
                                                                                                    
  SVv->N = 8.e-5;    /* abundance of Ne with respect to H (Raymond) */  
  SVv->nlev = 5;
 
  _1S0 = 0;
  _3P0 = 1;
  _3P1 = 2;
  _3P2 = 3;
  _1P1 = 4;
 
  SVv->nTom  = 4;
  SVv->Tom[0] = 0.5e4;
  SVv->Tom[1] = 1.0e4;
  SVv->Tom[2] = 1.5e4;
  SVv->Tom[3] = 2.0e4;

  if (SVv->dE == NULL){
    SVv->dE    = ARRAY_2D(SVv->nlev, SVv->nlev, double);
    SVv->A     = ARRAY_2D(SVv->nlev, SVv->nlev, double);
    SVv->omega = ARRAY_3D(SVv->nlev, SVv->nlev, SVv->nTom, double);
  }

/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < SVv->nlev ; i++) {
  for (j = 0 ; j < SVv->nlev ; j++) {
    SVv->dE[i][j] = -1.0;
    SVv->A[i][j]  = -1.0;
    for (n = 0; n < SVv->nTom; n++){
      SVv->omega[i][j][n]  = -1.0;
    }
  }}

   
  SVv->wght[_1S0] = 1.; /*    2 J  +  1     */
  SVv->wght[_3P0] = 1.; /*         */
  SVv->wght[_3P1] = 3.; /*         */
  SVv->wght[_3P2] = 5.; /*         */
  SVv->wght[_1P1] = 3.; /*         */

  SVv->omega[_3P2][_1S0][0] = 9.11e-1 * (5.0) / (3.0*3.0);
  SVv->omega[_3P2][_1S0][1] = 9.10e-1 * (5.0) / (3.0*3.0);
  SVv->omega[_3P2][_1S0][2] = 9.14e-1 * (5.0) / (3.0*3.0);
  SVv->omega[_3P2][_1S0][3] = 9.05e-1 * (5.0) / (3.0*3.0);
  
  SVv->omega[_3P1][_1S0][0] = 9.11e-1 * (3.0) / (3.0*3.0);
  SVv->omega[_3P1][_1S0][1] = 9.10e-1 * (3.0) / (3.0*3.0);
  SVv->omega[_3P1][_1S0][2] = 9.14e-1 * (3.0) / (3.0*3.0);
  SVv->omega[_3P1][_1S0][3] = 9.05e-1 * (3.0) / (3.0*3.0);

  SVv->omega[_3P0][_1S0][0] = 9.11e-1 * (1.0) / (3.0*3.0);
  SVv->omega[_3P0][_1S0][1] = 9.10e-1 * (1.0) / (3.0*3.0);
  SVv->omega[_3P0][_1S0][2] = 9.14e-1 * (1.0) / (3.0*3.0);
  SVv->omega[_3P0][_1S0][3] = 9.05e-1 * (1.0) / (3.0*3.0);

  SVv->omega[_1P1][_1S0][0] = 7.30;
  SVv->omega[_1P1][_1S0][1] = 7.30;
  SVv->omega[_1P1][_1S0][2] = 7.29;
  SVv->omega[_1P1][_1S0][3] = 7.27;
  
  SVv->omega[_1P1][_3P0][0] = 2.72e-1;
  SVv->omega[_1P1][_3P0][1] = 2.72e-1;
  SVv->omega[_1P1][_3P0][2] = 2.72e-1;
  SVv->omega[_1P1][_3P0][3] = 2.72e-1;

  SVv->omega[_3P2][_3P0][0] = 4.00e-1;
  SVv->omega[_3P2][_3P0][1] = 4.00e-1;
  SVv->omega[_3P2][_3P0][2] = 4.00e-1;
  SVv->omega[_3P2][_3P0][3] = 4.00e-1;

  SVv->omega[_3P2][_3P1][0] = 1.24;
  SVv->omega[_3P2][_3P1][1] = 1.24;
  SVv->omega[_3P2][_3P1][2] = 1.24;
  SVv->omega[_3P2][_3P1][3] = 1.24;
/*
  SVv->omega[_1P1][_3P1][0] = 2.72e-31;
  SVv->omega[_1P1][_3P1][1] = 2.72e-31;
  SVv->omega[_1P1][_3P1][2] = 2.72e-31;
  SVv->omega[_1P1][_3P1][3] = 2.72e-31;

  SVv->omega[_1P1][_3P2][0] = 2.72e-31;
  SVv->omega[_1P1][_3P2][1] = 2.72e-31;
  SVv->omega[_1P1][_3P2][2] = 2.72e-31;
  SVv->omega[_1P1][_3P2][3] = 2.72e-31;

  SVv->omega[_3P1][_3P0][0] = 2.72e-31;
  SVv->omega[_3P1][_3P0][1] = 2.72e-31;
  SVv->omega[_3P1][_3P0][2] = 2.72e-31;
  SVv->omega[_3P1][_3P0][3] = 2.72e-31;
*/
  
  SVv->dE[_3P2][_1S0] = 1188.3;   SVv->A[_3P2][_1S0]  = 6.59e-2;
  SVv->dE[_3P1][_1S0] = 1199.1;   SVv->A[_3P1][_1S0]  = 1.26e+5;
  SVv->dE[_3P0][_1S0] = 1204.5;   SVv->A[_3P0][_1S0]  = 0.00000;
  SVv->dE[_1P1][_1S0] = 786.48;   SVv->A[_1P1][_1S0]  = 5.25e+9;  
  SVv->dE[_1P1][_3P0] = 2.71e+5;  SVv->A[_1P1][_3P0]  = 9.16e-4;  
  SVv->dE[_3P2][_3P0] = 88401.7;  SVv->A[_3P2][_3P0]  = 0.00000;  
  SVv->dE[_3P2][_3P1] = 1.312e+5; SVv->A[_3P2][_3P1]  = 5.49e-3;  
/*  SVv->dE[_1P1][_3P1] = 1.6e+16;  SVv->A[_1P1][_3P1]  = 0.00000;  
  SVv->dE[_1P1][_3P2] = 1.6e+16;  SVv->A[_1P1][_3P2]  = 0.00000;  
  SVv->dE[_3P1][_3P0] = 1.6e+16;  SVv->A[_3P1][_3P0]  = 0.00000;  
*/

  Symmetrize_Coeff (SVv);
}

/* ********************************************************************* */
void FeI_INIT (Ion *FeIv)
/*
 *
 *********************************************************************** */
{
  int i, j, n;
  int L1, L2;
                                                                                                    
  FeIv->N = 3.e-5;    /* abundance of Fe with respect to H (Raymond) */  
  FeIv->nlev = 2;
 
  L1 = 0;
  L2 = 1;
 
  FeIv->nTom  = 4;
  FeIv->Tom[0] = 0.5e4;
  FeIv->Tom[1] = 1.0e4;
  FeIv->Tom[2] = 1.5e4;
  FeIv->Tom[3] = 2.0e4;

  if (FeIv->dE == NULL){
    FeIv->dE    = ARRAY_2D(FeIv->nlev, FeIv->nlev, double);
    FeIv->A     = ARRAY_2D(FeIv->nlev, FeIv->nlev, double);
    FeIv->omega = ARRAY_3D(FeIv->nlev, FeIv->nlev, FeIv->nTom, double);
  }
  
/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < FeIv->nlev ; i++) {
  for (j = 0 ; j < FeIv->nlev ; j++) {
    FeIv->dE[i][j] = -1.0;
    FeIv->A[i][j]  = -1.0;
    for (n = 0; n < FeIv->nTom; n++){
      FeIv->omega[i][j][n]  = -1.0;
    }
  }}


  FeIv->wght[L1] = 3.; /*    2 J  +  1     */
  FeIv->wght[L2] = 3.; /*         */

/*
  FeIv->omega[L2][L1][0] = 0.00;
  FeIv->omega[L2][L1][1] = 0.00;
  FeIv->omega[L2][L1][2] = 0.00;
  FeIv->omega[L2][L1][3] = 0.00; 
*/
  FeIv->dE[L2][L1] = 2.598e+5;  FeIv->A[L2][L1]  = 2.13e-3;  

  Symmetrize_Coeff (FeIv);
}

/* ********************************************************************* */
void FeII_INIT (Ion *FeIIv)
/*
 *
 *********************************************************************** */
{
  int i, j, n;
  int _6D12, _6D32, _6D52, _6D72, _6D92, _4F92, _4D72, _4D52, _4D12, _4P52, _4F72, _4D32;

                                                                                                    
  FeIIv->N = 3.e-5;    /* abundance of Fe with respect to H (Raymond) */  
  FeIIv->nlev = 6;

 
  _6D92 = 0;
  _6D72 = 1;
  _4F92 = 2;
  _4F72 = 3;
  _4D32 = 4;
  _4P52 = 5; 

  FeIIv->nTom  = 3;
  FeIIv->Tom[0] = 2.37e4;
  FeIIv->Tom[1] = 3.15e4;
  FeIIv->Tom[2] = 3.94e4;


  if (FeIIv->dE == NULL){
    FeIIv->dE    = ARRAY_2D(FeIIv->nlev, FeIIv->nlev, double);
    FeIIv->A     = ARRAY_2D(FeIIv->nlev, FeIIv->nlev, double);
    FeIIv->omega = ARRAY_3D(FeIIv->nlev, FeIIv->nlev, FeIIv->nTom, double);
  }
  
/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < FeIIv->nlev ; i++) {
  for (j = 0 ; j < FeIIv->nlev ; j++) {
    FeIIv->dE[i][j] = -1.0;
    FeIIv->A[i][j]  = -1.0;
    for (n = 0; n < FeIIv->nTom; n++)
      FeIIv->omega[i][j][n]  = -1.0;
  }}


  FeIIv->wght[_6D92] = 10.; /*    2 J  +  1     */
  FeIIv->wght[_6D72] = 8.; /*         */
  FeIIv->wght[_4F92] = 10.; /*         */
  FeIIv->wght[_4F72] = 8.; /*         */
  FeIIv->wght[_4D32] = 4.; /*         */
  FeIIv->wght[_4P52] = 6.; /*         */


  FeIIv->omega[_6D72][_6D92][0] = 3.043;
  FeIIv->omega[_6D72][_6D92][1] = 3.046;
  FeIIv->omega[_6D72][_6D92][2] = 3.049;

  FeIIv->omega[_4F92][_6D92][0] = 1.811;
  FeIIv->omega[_4F92][_6D92][1] = 1.742;
  FeIIv->omega[_4F92][_6D92][2] = 1.671;

  FeIIv->omega[_4P52][_4F92][0] = 1.967;
  FeIIv->omega[_4P52][_4F92][1] = 2.296;
  FeIIv->omega[_4P52][_4F92][2] = 3.612;

  FeIIv->omega[_4D32][_4F72][0] = 0.457;  
  FeIIv->omega[_4D32][_4F72][1] = 0.447;  
  FeIIv->omega[_4D32][_4F72][2] = 0.429;  

  FeIIv->omega[_4F72][_4F92][0] = 3.521;  
  FeIIv->omega[_4F72][_4F92][1] = 3.405;  
  FeIIv->omega[_4F72][_4F92][2] = 3.330;  

  FeIIv->omega[_4F72][_6D72][0] = 0.645;  
  FeIIv->omega[_4F72][_6D72][1] = 0.620;  
  FeIIv->omega[_4F72][_6D72][2] = 0.595;  

  FeIIv->omega[_4F92][_6D72][0] = 0.742;  
  FeIIv->omega[_4F92][_6D72][1] = 0.714;  
  FeIIv->omega[_4F92][_6D72][2] = 0.685;  

  FeIIv->omega[_4F72][_6D92][0] = 0.663;  
  FeIIv->omega[_4F72][_6D92][1] = 0.638;  
  FeIIv->omega[_4F72][_6D92][2] = 0.612;  

  FeIIv->omega[_4D32][_6D92][0] = 0.132;  
  FeIIv->omega[_4D32][_6D92][1] = 0.132;  
  FeIIv->omega[_4D32][_6D92][2] = 0.131;  

  FeIIv->omega[_4P52][_6D72][0] = 0.622;  
  FeIIv->omega[_4P52][_6D72][1] = 0.594;  
  FeIIv->omega[_4P52][_6D72][2] = 0.567;  

  FeIIv->omega[_4P52][_4F72][0] = 1.224;  
  FeIIv->omega[_4P52][_4F72][1] = 1.554;  
  FeIIv->omega[_4P52][_4F72][2] = 2.211;  

  FeIIv->omega[_4P52][_4D32][0] = 0.241;  
  FeIIv->omega[_4P52][_4D32][1] = 0.220;  
  FeIIv->omega[_4P52][_4D32][2] = 0.216;   

/* add now some fake data for not defined transitions - always keep A=0 */
/*  for (i=1; i<FeIIv->nlev; i++)
    for (j=0; j<i; j++)
      if ( (FeIIv->omega[i][j][0] < 0.0) && (FeIIv->omega[j][i][0] < 0.0) ) {
        for (n=0; n<FeIIv->nTom; n++)
          FeIIv->omega[i][j][n]=1.0e-30;
        FeIIv->A[i][j]  = 0.0;
        FeIIv->dE[i][j] = 1.0e+30;      
      }
*/

  FeIIv->dE[_6D72][_6D92] = 2218.26;  FeIIv->A[_6D72][_6D92]  = 1.9e+8;  
  FeIIv->dE[_4F92][_6D92] = 2260.9;  FeIIv->A[_4F92][_6D92]  = 4.9e+6;  
  FeIIv->dE[_4P52][_4F92] = 8617.0;    FeIIv->A[_4P52][_4F92]  = 1.9e-2; 
  FeIIv->dE[_4D32][_4F72] = 1.5994e+4; FeIIv->A[_4D32][_4F72]  = 1.40e-3;
  FeIIv->dE[_4F72][_6D92] = 4382.7; FeIIv->A[_4F72][_6D92] = 5.5e-2;
  FeIIv->dE[_4D32][_6D92] = dEtoA(1.076239); FeIIv->A[_4D32][_6D92] = 0.0;
  FeIIv->dE[_4P52][_6D92] = dEtoA(1.6706139); FeIIv->A[_4P52][_6D92] = 0.0;
  FeIIv->dE[_4F92][_6D72] = 2279.9; FeIIv->A[_4F92][_6D72] = 3.9e+6;
  FeIIv->dE[_4F72][_6D72] = 2253.1; FeIIv->A[_4F72][_6D72] = 5.1e+6;
  FeIIv->dE[_4D32][_6D72] = dEtoA(1.029); FeIIv->A[_4D32][_6D72] = 0.0;
  FeIIv->dE[_4P52][_6D72] = dEtoA(0.99); FeIIv->A[_4P52][_6D72] = 1.5e-3;
  FeIIv->dE[_4F72][_4F92] = 2331.3; FeIIv->A[_4F72][_4F92] = 2.9e+7;
  FeIIv->dE[_4D32][_4F92] = dEtoA(0.84); FeIIv->A[_4D32][_4F92] = 0.0;
  FeIIv->dE[_4P52][_4F72] = dEtoA(0.69); FeIIv->A[_4P52][_4F72] = 4.94e-3;
  FeIIv->dE[_4P52][_4D32] = dEtoA(0.592); FeIIv->A[_4P52][_4D32] = 7.71e-5;

/*
  FeIIv->dE[_6D52][_6D92] = dEtoA(0.0827821); FeIIv->A[_6D52][_6D92] = 0.0;
  FeIIv->dE[_6D32][_6D92] = dEtoA(0.106950); FeIIv->A[_6D32][_6D92] = 0.0;
  FeIIv->dE[_6D12][_6D92] = dEtoA(0.121139); FeIIv->A[_6D12][_6D92] = 0.0;
  FeIIv->dE[_4F72][_6D92] = dEtoA(0.3012936); FeIIv->A[_4F72][_6D92] = 0.0;
  FeIIv->dE[_4D52][_6D92] = dEtoA(1.040468); FeIIv->A[_4D52][_6D92] = 0.0;
  FeIIv->dE[_4D32][_6D92] = dEtoA(1.076239); FeIIv->A[_4D32][_6D92] = 0.0;
  FeIIv->dE[_4D12][_6D92] = dEtoA(1.096859); FeIIv->A[_4D12][_6D92] = 0.0;
  FeIIv->dE[_4P52][_6D92] = dEtoA(1.6706139); FeIIv->A[_4P52][_6D92] = 0.0;
  FeIIv->dE[][_6D72] = ; FeIIv->A[][] = ;
  FeIIv->dE[][] = ; FeIIv->A[][] = ;
  FeIIv->dE[][] = ; FeIIv->A[][] = ;
  FeIIv->dE[][] = ; FeIIv->A[][] = ;
  FeIIv->dE[][] = ; FeIIv->A[][] = ;

  _6D92 = 0;
  _6D72 = 1;
  _4F92 = 2;
  _4F72 = 3;
  _4D32 = 4;
  _4P52 = 5; 

*/
  Symmetrize_Coeff (FeIIv);
}

/* ********************************************************************* */
void FeIIb_INIT (Ion *FeIIv)
/*
 *
 *********************************************************************** */
{
  int i, j, n;
  int _6D12, _6D32, _6D52, _6D72, _6D92, _4F92, _4D72, _4D52, _4D12, _4P52, _4F72, _4D32;
  double lev_en[12] = {0.0, 0.047708, 0.0837821, 0.106950, 0.121139, 0.2321687, 0.3012936, 0.9863313, 1.040468, 1.076239, 1.096859, 1.6706139};

                                                                                                    
  FeIIv->N = 3.e-5;    /* abundance of Fe with respect to H (Raymond) */  
  FeIIv->nlev = 12;

 
  _6D92 = 0;
  _6D72 = 1;
  _6D52 = 2;
  _6D32 = 3;
  _6D12 = 4;
  _4F92 = 5;
  _4F72 = 6;
  _4D72 = 7;
  _4D52 = 8;
  _4D32 = 9;
  _4D12 = 10;
  _4P52 = 11; 

  FeIIv->nTom  = 4;
  FeIIv->Tom[0] = 0.5e4;
  FeIIv->Tom[1] = 1.0e4;
  FeIIv->Tom[2] = 1.5e4;
  FeIIv->Tom[3] = 2.0e4;

  if (FeIIv->dE == NULL){
    FeIIv->dE    = ARRAY_2D(FeIIv->nlev, FeIIv->nlev, double);
    FeIIv->A     = ARRAY_2D(FeIIv->nlev, FeIIv->nlev, double);
    FeIIv->omega = ARRAY_3D(FeIIv->nlev, FeIIv->nlev, FeIIv->nTom, double);
  }
  
/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < FeIIv->nlev ; i++) {
  for (j = 0 ; j < FeIIv->nlev ; j++) {
    if (j <= i) {
      FeIIv->dE[i][j] = -1.0;
      FeIIv->A[i][j]  = -1.0;
      for (n = 0; n < FeIIv->nTom; n++)
        FeIIv->omega[i][j][n]  = -1.0;
    } else {
      FeIIv->dE[i][j] = dEtoA(lev_en[j] - lev_en[i]);
      FeIIv->A[i][j]  = 0.0;
      for (n = 0; n < FeIIv->nTom; n++)
        FeIIv->omega[i][j][n]  = 0.0;
    }
  }}


  FeIIv->wght[_6D92] = 10.; /*    2 J  +  1     */
  FeIIv->wght[_6D72] = 8.; /*         */
  FeIIv->wght[_6D52] = 6.; /*         */
  FeIIv->wght[_6D32] = 4.; /*         */
  FeIIv->wght[_6D12] = 4.; /*         */
  FeIIv->wght[_4F92] = 10.; /*         */
  FeIIv->wght[_4F72] = 8.; /*         */
  FeIIv->wght[_4D72] = 8.; /*         */
  FeIIv->wght[_4D52] = 6.; /*         */
  FeIIv->wght[_4D32] = 4.; /*         */
  FeIIv->wght[_4D12] = 2.;
  FeIIv->wght[_4P52] = 6.; /*         */


  FeIIv->omega[_6D72][_6D92][0] = 7.41e-3;
  FeIIv->omega[_6D72][_6D92][1] = 5.52;
  FeIIv->omega[_6D72][_6D92][2] = 5.46;
  FeIIv->omega[_6D72][_6D92][3] = 5.48; 

  FeIIv->omega[_4F92][_6D92][0] = 4.77e-2;
  FeIIv->omega[_4F92][_6D92][1] = 3.60;
  FeIIv->omega[_4F92][_6D92][2] = 3.19;
  FeIIv->omega[_4F92][_6D92][3] = 2.89;

  FeIIv->omega[_4D72][_6D92][0] = 3.10;
  FeIIv->omega[_4D72][_6D92][1] = 1.10e+1;
  FeIIv->omega[_4D72][_6D92][2] = 9.65;
  FeIIv->omega[_4D72][_6D92][3] = 8.77;

  FeIIv->omega[_4D52][_6D52][0] = 3.46;
  FeIIv->omega[_4D52][_6D52][1] = 2.85;
  FeIIv->omega[_4D52][_6D52][2] = 2.56;
  FeIIv->omega[_4D52][_6D52][3] = 2.36;

  FeIIv->omega[_4D52][_6D32][0] = 7.11e-1;
  FeIIv->omega[_4D52][_6D32][0] = 6.26e-1;
  FeIIv->omega[_4D52][_6D32][0] = 6.05e-1;
  FeIIv->omega[_4D52][_6D32][0] = 5.87e-1;

  FeIIv->omega[_4D12][_6D12][0] = 2.06;
  FeIIv->omega[_4D12][_6D12][1] = 1.65;
  FeIIv->omega[_4D12][_6D12][2] = 1.44;
  FeIIv->omega[_4D12][_6D12][3] = 1.30;  

  FeIIv->omega[_4D72][_4F92][0] = 2.28;
  FeIIv->omega[_4D72][_4F92][1] = 2.11;
  FeIIv->omega[_4D72][_4F92][2] = 2.14;
  FeIIv->omega[_4D72][_4F92][3] = 2.18;

  FeIIv->omega[_4D52][_4F92][0] = 7.77e-1;
  FeIIv->omega[_4D52][_4F92][1] = 7.25e-1;
  FeIIv->omega[_4D52][_4F92][2] = 7.54e-1;
  FeIIv->omega[_4D52][_4F92][3] = 7.83e-1;
  
  FeIIv->omega[_4P52][_4F92][0] = 1.22;
  FeIIv->omega[_4P52][_4F92][1] = 1.23;
  FeIIv->omega[_4P52][_4F92][2] = 1.27;
  FeIIv->omega[_4P52][_4F92][3] = 1.29;

  FeIIv->omega[_4D32][_4F72][0] = 0.457;  
  FeIIv->omega[_4D32][_4F72][1] = 0.447;  
  FeIIv->omega[_4D32][_4F72][2] = 0.429;  
  FeIIv->omega[_4D32][_4F72][3] = 0.429;  

/* add now some fake data for not defined transitions - always keep A=0 */
  for (i=1; i<FeIIv->nlev; i++)
    for (j=0; j<i; j++)
      if ( (FeIIv->omega[i][j][0] < 0.0) && (FeIIv->omega[j][i][0] < 0.0) ) {
        for (n=0; n<FeIIv->nTom; n++)
          FeIIv->omega[i][j][n]=1.0e-30;
        FeIIv->A[i][j]  = 0.0;
        FeIIv->dE[i][j] = 1.0e+30;      
      }


  FeIIv->dE[_6D72][_6D92] = 2.598e+5;  FeIIv->A[_6D72][_6D92]  = 2.13e-3;  
  FeIIv->dE[_4F92][_6D92] = 5.339e+4;  FeIIv->A[_4F92][_6D92]  = 4.17e-5;  
  FeIIv->dE[_4D72][_6D92] = 1.257e+4;  FeIIv->A[_4D72][_6D92]  = 4.83e-3;  
  FeIIv->dE[_4D52][_6D52] = 1.294e+4;  FeIIv->A[_4D52][_6D52]  = 1.94e-3;  
  FeIIv->dE[_4D52][_6D32] = 1.328e+4;  FeIIv->A[_4D52][_6D32]  = 1.21e-3;  
  FeIIv->dE[_4D12][_6D12] = 1.270e+4;  FeIIv->A[_4D12][_6D12]  = 2.91e-3;  
  FeIIv->dE[_4D72][_4F92] = 1.644e+4;  FeIIv->A[_4D72][_4F92]  = 4.65e-3;  
  FeIIv->dE[_4D52][_4F92] = 1.533e+4;  FeIIv->A[_4D52][_4F92]  = 2.44e-3;  
  FeIIv->dE[_4P52][_4F92] = 8617.0;    FeIIv->A[_4P52][_4F92]  = 2.73e-2; 
  FeIIv->dE[_4D32][_4F72] = 1.5994e+4; FeIIv->A[_4D32][_4F72]  = 1.40e-3;
/*
  FeIIv->dE[_6D52][_6D92] = dEtoA(0.0827821); FeIIv->A[_6D52][_6D92] = 0.0;
  FeIIv->dE[_6D32][_6D92] = dEtoA(0.106950); FeIIv->A[_6D32][_6D92] = 0.0;
  FeIIv->dE[_6D12][_6D92] = dEtoA(0.121139); FeIIv->A[_6D12][_6D92] = 0.0;
  FeIIv->dE[_4F72][_6D92] = dEtoA(0.3012936); FeIIv->A[_4F72][_6D92] = 0.0;
  FeIIv->dE[_4D52][_6D92] = dEtoA(1.040468); FeIIv->A[_4D52][_6D92] = 0.0;
  FeIIv->dE[_4D32][_6D92] = dEtoA(1.076239); FeIIv->A[_4D32][_6D92] = 0.0;
  FeIIv->dE[_4D12][_6D92] = dEtoA(1.096859); FeIIv->A[_4D12][_6D92] = 0.0;
  FeIIv->dE[_4P52][_6D92] = dEtoA(1.6706139); FeIIv->A[_4P52][_6D92] = 0.0;
  FeIIv->dE[][_6D72] = ; FeIIv->A[][] = ;
  FeIIv->dE[][] = ; FeIIv->A[][] = ;
  FeIIv->dE[][] = ; FeIIv->A[][] = ;
  FeIIv->dE[][] = ; FeIIv->A[][] = ;
  FeIIv->dE[][] = ; FeIIv->A[][] = ;

  _6D92 = 0;
  _6D72 = 1;
  _6D52 = 2;
  _6D32 = 3;
  _6D12 = 4;
  _4F92 = 5;
  _4F72 = 6;
  _4D72 = 7;
  _4D52 = 8;
  _4D32 = 9;
  _4D12 = 10;
  _4P52 = 11;
*/
  Symmetrize_Coeff (FeIIv);
}

/* ********************************************************************* */
void FeIII_INIT (Ion *FeIIIv)
/*
 *
 *********************************************************************** */
{
  int i, j, n;
  int L1, L2;
                                                                                                    
  FeIIIv->N = 3.e-5;    /* abundance of Fe with respect to H (Raymond) */  
  FeIIIv->nlev = 2;
 
  L1 = 0;
  L2 = 1;
 
  FeIIIv->nTom  = 4;
  FeIIIv->Tom[0] = 0.5e4;
  FeIIIv->Tom[1] = 1.0e4;
  FeIIIv->Tom[2] = 1.5e4;
  FeIIIv->Tom[3] = 2.0e4;

  if (FeIIIv->dE == NULL){
    FeIIIv->dE    = ARRAY_2D(FeIIIv->nlev, FeIIIv->nlev, double);
    FeIIIv->A     = ARRAY_2D(FeIIIv->nlev, FeIIIv->nlev, double);
    FeIIIv->omega = ARRAY_3D(FeIIIv->nlev, FeIIIv->nlev, FeIIIv->nTom, double);
  }
  
/*  INITIALIZE MATRIX WITH NEGATIVE NUMBER,
    O GET PROPER ORDER LATER  */
 
  for (i = 0 ; i < FeIIIv->nlev ; i++) {
  for (j = 0 ; j < FeIIIv->nlev ; j++) {
    FeIIIv->dE[i][j] = -1.0;
    FeIIIv->A[i][j]  = -1.0;
    for (n = 0; n < FeIIIv->nTom; n++){
      FeIIIv->omega[i][j][n]  = -1.0;
    }
  }}


  FeIIIv->wght[L1] = 3.; /*    2 J  +  1     */
  FeIIIv->wght[L2] = 3.; /*         */


  FeIIIv->omega[L2][L1][0] = 2.87;
  FeIIIv->omega[L2][L1][1] = 2.92;
  FeIIIv->omega[L2][L1][2] = 2.95;
  FeIIIv->omega[L2][L1][3] = 2.95; 

  FeIIIv->dE[L2][L1] = 3239.74;  FeIIIv->A[L2][L1]  = 2.3e-1;  

  Symmetrize_Coeff (FeIIIv);
}
