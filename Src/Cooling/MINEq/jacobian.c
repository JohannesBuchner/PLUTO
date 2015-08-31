#include "pluto.h"
#include "cooling_defs.h"

/* ********************************************************************* */
void Jacobian (double *v, double *rhs, double **dfdy)
/*!
 * Compute the jacobian J(k,l) = dfdy
 *
 * k = row index
 * l = column index
 *
 *    J(0,0)   J(0,1) ...  J(0, n-1)
 *    J(1,0)   J(1,1) .... J(1, n-1)
 *      .         .           .
 *      .         .           .
 *      .         .           .
 *    J(n-1,0)  ....      J(n-1, n-1)
 *
 *
 *   or, 
 *
 *   +-----------------------+
 *   +              |        |
 *   +              |        |
 *   +              |        |
 *   +    dX'/dX    | dX'/dp |
 *   +     (JXX)    |  (JXp) |
 *   +              |        |
 *   +              |        |
 *   +--------------+--------+
 *   +   dp'/dX     | dp'/dp | 
 *   +    (JpX)     |  Jpp   |
 *   +-----------------------+
 *
 *
 *
 *********************************************************************** */
{
  int    k, l, nv, n;
  double   dmu_dX[NIONS], N, mu;
  double   T, dLdX, dCdX, dRdX, eps, scrh;
  double   vp[NVAR], vm[NVAR], rhs_m[NVAR], rhs_p[NVAR];
  double   dfdp[NIONS];
  double   *L, *La, *Lb, *Lc;
  double   *C, *Ca, *Cb, *Cc;
  double   *R, *Ra, *Rb, *Rc;
  double   *X, *dnel_dX, *de;
  double   Unit_Time, E_cost;
  static double **delta;

  n = NIONS + 1;

  Unit_Time = UNIT_LENGTH/UNIT_VELOCITY;

  E_cost    = UNIT_LENGTH/UNIT_DENSITY/
              (UNIT_VELOCITY*UNIT_VELOCITY*UNIT_VELOCITY);

  if (delta == NULL) {
    delta = ARRAY_2D(NIONS, NIONS, double);
    for (k = 0; k < NIONS; k++){
    for (l = 0; l < NIONS; l++){
       delta[k][l] = (k == l ? 1.0:0.0);
    }}
  }

  Radiat (v, rhs);

  L = CoolCoeffs.Lrate; La = CoolCoeffs.La; Lb = CoolCoeffs.Lb; Lc = CoolCoeffs.Lc;
  C = CoolCoeffs.Crate; Ca = CoolCoeffs.Ca; Cb = CoolCoeffs.Cb; Cc = CoolCoeffs.Cc;
  R = CoolCoeffs.Rrate; Ra = CoolCoeffs.Ra; Rb = CoolCoeffs.Rb; Rc = CoolCoeffs.Rc;

  de = CoolCoeffs.de;

  X       = v + NFLX;
  dnel_dX = CoolCoeffs.dnel_dX;

  N  = v[RHO]*find_N_rho();       /* -- Total number density -- */
  mu = MeanMolecularWeight(v); 
  T  = v[PRS]/v[RHO]*KELVIN*mu;

/* --  compute the vector grad_X (mu)  --  */

  scrh = 1.0/CoolCoeffs.muD;
  for (k = 0; k < n - 1; k++){
    dmu_dX[k] =  (CoolCoeffs.dmuN_dX[k] 
                - CoolCoeffs.dmuD_dX[k]*mu)*scrh;
  }
  
/* -------------------------------------------------
         Compute dX'/dX, dp'/dX with 
 
      k = 0....n - 2
      l = 0....n - 2
   ------------------------------------------------- */

  /* -- compute d\dot{X}_0/dXl (H right hand side) separately -- */
  
  for (l = 0; l < n - 1; l++) {  /* -- first row -- */
    dfdy[0][l] = - (R[0] + C[0])*delta[0][l] + 
                 (1.0 - X[0])*CoolCoeffs.fRH*dnel_dX[l]
                      - X[0] *CoolCoeffs.fCH*dnel_dX[l];
  }

  for (k = 1; k < n - 1; k++) {  /* -- C's and L's contribution -- */
  for (l = 0; l < n - 1; l++) {
    dLdX       = La[k]*dnel_dX[l]      + Lb[k]*delta[0][l] + Lc[k]*delta[1][l];
    dCdX       = Ca[k]*dnel_dX[l]      + Cb[k]*delta[0][l] + Cc[k]*delta[1][l];
    dfdy[k][l] =  L[k]*delta[k - 1][l] + dLdX*X[k - 1] - C[k]*delta[k][l] - dCdX*X[k];
  }}

  for (k = 1; k < n - 2; k++) { /* -- R's contribution -- */
  for (l = 0; l < n - 1; l++) {
    dRdX        = Ra[k]*dnel_dX[l] + Rb[k]*delta[0][l] + Rc[k]*delta[1][l];
    dfdy[k][l] +=  R[k]*delta[k + 1][l] + dRdX*X[k + 1];
  }}

  for (k = 0; k < n - 1; k++) { /* -- Get correct dimensions -- */
  for (l = 0; l < n - 1; l++) {
    dfdy[k][l] *= Unit_Time;
  }}

/* -------------------------------------------------------
    Compute PDEs of cooling function with respect to X's
   ------------------------------------------------------- */

  for (l = 0; l < n - 1; l++){

    dfdy[n - 1][l] = dnel_dX[l]*rhs[PRS]/CoolCoeffs.Ne;

    scrh = 0.0; 
    for (nv = 0; nv < NIONS; nv++){ 
      scrh += X[nv]*CoolCoeffs.de_dne[nv]*elem_ab[elem_part[nv]];
    }
    scrh *= dnel_dX[l];
    scrh += de[l]*elem_ab[elem_part[l]];
    scrh += CoolCoeffs.dLIR_dX[0]*delta[0][l];
    scrh += CoolCoeffs.dLIR_dX[2]*delta[2][l];
    scrh *= N*CoolCoeffs.Ne*E_cost;

    dfdy[n - 1][l] -= (g_gamma-1.0)*scrh;
  }

/* --------------------------------------------------------
      Compute partial derivatives with respect to T
      using numerical differentiation.
   -------------------------------------------------------- */

  eps = 1.e-4;
  
  for (nv = 0; nv < NVAR; nv++){
    vm[nv] = vp[nv] = v[nv];
  }
  vp[PRS] = v[PRS]*(1.0 + eps);
  vm[PRS] = v[PRS]*(1.0 - eps);

  Radiat (vp, rhs_p);
  Radiat (vm, rhs_m);

/* -- Compute last column (Jxp and Jpp)  drhs/dp -- */

  for (k = 0; k < n - 1; k++){
    dfdp[k]        = (rhs_p[k + NFLX] - rhs_m[k + NFLX])/(2.0*eps*v[PRS]);
    dfdy[k][n - 1] = dfdp[k]; 
  }
  dfdp[n - 1]        = (rhs_p[PRS] - rhs_m[PRS])/(2.0*eps*v[PRS]);
  dfdy[n - 1][n - 1] = dfdp[n - 1];
	
/* -- Add df/dT*dT/dX term to all columns except last one -- */

  scrh = v[PRS]/mu;
  for (k = 0; k < n    ; k++) {
  for (l = 0; l < n - 1; l++) {
    dfdy[k][l] += dfdp[k]*scrh*dmu_dX[l];
  }}

}


