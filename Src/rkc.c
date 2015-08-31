/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Runge-Kutta-Chebyshev driver for integration of diffusion terms.

  Take one step in the solution of the diffusion equation \f$ dQ/dt = F\f$
  using Runge-Kutta-Chebyshev (RKC) method.
  \f[ \begin{array}{lcl}
      Y_0  & = & Q^n                       \\ \noalign{\medskip}
      Y_1  & = & Y_0 + \tilde{\mu}_1\tau F_0   \\ \noalign{\medskip}
      Y_j  & = & (1 - \mu_j - \nu_j) Y_0 + \mu_j Y_{j-1} + \nu_j Y_{j-2}
             + \tilde{\mu}_j \tau F_{j-1} + \tilde{\gamma}_j\tau F_0  
            \quad \rm{for} \quad j = 2,...,s    \\ \noalign{\medskip}
    Q^{n+1} & = & Y_s
      \end{array}
  \f]
  where 
  \f[\tau = \Delta t,\, 
     \tilde{\mu}_1 = b_1 w_1,\,
     \mu_j = 2 b_j \frac{w_0}{b_{j-1}}, \,
     \nu_j = -\frac{b_j}{b_{j-2}},\,
     \tilde{\mu}_j = 2 b_j \frac{w_1}{b_{j-1}},\,
     \tilde{\gamma}_j = -\alpha_{j-1} \tilde{\mu}_j,\,
  \f]
  \f[   
     \alpha_j = 1 - b_j T_j(w_0),\,
     w_0 = 1 + \epsilon/s^2,\,
     w_1 = T'_s(w_0)/T''_s(w_0),\,
     b_j = T''_j(w_0)/(T'_j(w_0))^2,\,
     b0=b1=b2,\,
     \epsilon = 2/13     
   \f]  
   Here the \e T's are the Chebyshev polynomial of the first kind:
   \f[
     T_j(x) = \cos(j {\rm arccos(x)})  ,\,
     T_0(x) = 1,\,
     T_1(x)=x,\,
     T_j(x)=2*x*T_jm1(x) - T_jm2(x), \quad \rm{for}\quad j = 2,...,s  
   \f]
    while sprad=Spectral_Radius-> 4,8,12 kappa/dx/dx (check)

  Only the parabolic terms of the equations are advanced by a 
  single step ::g_dt equal to the current integration time step and the 
  number of RKC steps (Nrkc) is computed...
  The explicit parabolic time step is computed in the same manner as in
  ::STS.
  
 \b References

  \authors P. Tzeferacos (petros.tzeferacos@ph.unito.it)\n
           A. Mignone (mignone@ph.unito.it)
  \date    Oct 29, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void RKC (const Data *d, Time_Step *Dts, Grid *grid)
/*!
 * Solve diffusion equation using Runge-Kutta-Chebyshev (RKC) method
 *
 * \param [in,out]  d    pointer to Data structure
 * \param [in,out]  Dts  pointer to Time_Step structure  
 * \param [in]     grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int    i, j, k, nv, s, s_RKC = 0;
  int    nv_indx, var_list[NVAR], nvar_rkc;
  double kappa_dl2, eta_dl2, nu_dl2;
  double absh, sprad, Y;
  double epsilon = 2./13.;
  double w_0, w_1;
  double scrh, scr1, scr2;
  double mu_j, nu_j, mu_tilde,a_jm1;
  double Tj, dTj, d2Tj, Tjm1, Tjm2, dTjm1, dTjm2, d2Tjm1, d2Tjm2;
  double b_j, b_jm1, b_jm2, arg;
  double dt_par; 
  static Data_Arr Y_jm1, Y_jm2, UU_0, F_0, F_jm1;
  static double **v;
  static unsigned char *flag;
  
  if (Y_jm1 == NULL) { 
    Y_jm1  = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);  
    Y_jm2  = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    UU_0   = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    F_0    = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    F_jm1  = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    flag   = ARRAY_1D(NMAX_POINT, unsigned char);
    v      = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

/* -------------------------------------------------------
     set the number of variables to be evolved with RKC
   ------------------------------------------------------- */
   
  i = 0;
  for (nv = 0; nv < NVAR; nv++) var_list[nv] = 0;
  #if VISCOSITY == RK_CHEBYSHEV
   EXPAND(var_list[i++] = MX1;  ,
          var_list[i++] = MX2;  ,
          var_list[i++] = MX3;)
  #endif

  #if RESISTIVITY == RK_CHEBYSHEV
   EXPAND(var_list[i++] = BX1;  ,
          var_list[i++] = BX2;  ,
          var_list[i++] = BX3;)
  #endif
  #if    (THERMAL_CONDUCTION == RK_CHEBYSHEV)  \
      || (VISCOSITY     == RK_CHEBYSHEV && EOS == IDEAL) \
      || (RESISTIVITY == RK_CHEBYSHEV && EOS == IDEAL)
   
   var_list[i++] = ENG;
  #endif
  nvar_rkc = i;
  
/* -------------------------------------------------
          get conservative vector UU 
   ------------------------------------------------- */

  KDOM_LOOP(k){
  JDOM_LOOP(j){
    IDOM_LOOP(i){
    for (nv = NVAR; nv--;   ) {
      v[i][nv] = d->Vc[nv][k][j][i];
    }}
    PrimToCons (v, UU_0[k][j], IBEG, IEND);
  }}   
  
/* -------------------------------------------------
    Compute the parabolic time step by calling 
    the RHS function once outside the main loop.
   ------------------------------------------------- */

  g_intStage = 1;
  Boundary(d, ALL_DIR, grid);
  #if SHOCK_FLATTENING == MULTID
   FindShock (d, grid);
  #endif

  Dts->inv_dtp  = ParabolicRHS(d, F_0, 1.0, grid);
  Dts->inv_dtp /= (double) DIMENSIONS;  
  #ifdef PARALLEL
   MPI_Allreduce (&Dts->inv_dtp, &scrh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   Dts->inv_dtp = scrh;
  #endif
  dt_par = Dts->cfl_par/(2.0*Dts->inv_dtp); /* -- explicit parabolic time step -- */   

/* -----------------------------------------------------
      Compute spectral radius (placeholder for now...)
   -----------------------------------------------------*/  

  D_EXPAND( sprad = 4.*Dts->inv_dtp;,
            sprad = 8.*Dts->inv_dtp;,
            sprad = 12.*Dts->inv_dtp;)

D_EXPAND( sprad = 2./dt_par; ,
          sprad = 4./dt_par; ,
          sprad = 6./dt_par;)

/* ----------------------------------------------------- 
    Compute number of RKC steps from sprad
   -----------------------------------------------------*/   

  absh = g_dt;   
  s_RKC = 1 + (int)(sqrt(1.54*absh*sprad + 1));
  Dts->Nrkc = s_RKC;    
/* limit dt */
  if(absh > 0.653*s_RKC*s_RKC/sprad){
    print1 ("! RKC: outside parameter range\n");
    QUIT_PLUTO(1);
  }    

/* ------------------------------------------------------------------ 
                     RKC main stage loop 
   ------------------------------------------------------------------ */

/* -- define coefficients -- */

  w_0 = 1. + epsilon/(1.*s_RKC)/(1.*s_RKC);      
  scr1 = w_0*w_0 - 1.;
  scr2 = sqrt(scr1);
  arg = s_RKC*log(w_0 + scr2);
  w_1 = sinh(arg)*scr1/(cosh(arg)*s_RKC*scr2 - w_0*sinh(arg));
  b_jm1 = 1./(2.*w_0)/(2.*w_0);
  b_jm2 = b_jm1;
  mu_tilde = w_1*b_jm1;

/* -------------------------------------------------------------
     First RKC stage here:
     - take one step in conservative variables,
     - convert to primitive 
   ------------------------------------------------------------- */

  DOM_LOOP (k,j,i){
    for (nv = 0; nv < NVAR; nv++) {
      Y_jm1[k][j][i][nv] = Y_jm2[k][j][i][nv] = UU_0[k][j][i][nv];
    }
    for (nv_indx = 0; nv_indx < nvar_rkc; nv_indx++){
      nv = var_list[nv_indx];
      Y_jm1[k][j][i][nv] = Y_jm2[k][j][i][nv] + mu_tilde*absh*F_0[k][j][i][nv];
    }
  }
	
  KDOM_LOOP(k){
  JDOM_LOOP(j){
    ConsToPrim (Y_jm1[k][j], v, IBEG, IEND, flag);
    IDOM_LOOP(i){
    for (nv = NVAR; nv--;   ){
      d->Vc[nv][k][j][i] = v[i][nv];
    }}
  }}
  
/* ----------------------------------------------------------------- 
    Loop over remaining RKC stages s = 2,..,s_RKC steps 
   ----------------------------------------------------------------- */

/* -- Initialize Chebyshev polynomials (recursive) -- */

  Tjm1   = w_0; Tjm2   = 1.0; dTjm1  = 1.0;
  dTjm2  = 0.0; d2Tjm1 = 0.0; d2Tjm2 = 0.0;
  
  for (s = 2; s <= s_RKC; s++){
  
  /* ---- compute coefficients ---- */
  
    Tj    =   2.*w_0*Tjm1 - Tjm2;
    dTj   =   2.*w_0*dTjm1 - dTjm2 + 2.*Tjm1;
    d2Tj  =   2.*w_0*d2Tjm1 - d2Tjm2 + 4.*dTjm1;
    b_j   =   d2Tj/dTj/dTj;
    a_jm1 =   1. - Tjm1*b_jm1;
    mu_j  =   2.*w_0*b_j/b_jm1;
    nu_j  = - b_j/b_jm2;
    mu_tilde   =   mu_j*w_1/w_0;
    
  /* -- call again boundary and take a new step -- */
  
    g_intStage = s;
    Boundary(d, ALL_DIR, grid);
    ParabolicRHS (d, F_jm1, 1.0, grid);
    DOM_LOOP (k,j,i){
      for (nv_indx = 0; nv_indx < nvar_rkc; nv_indx++){
        nv = var_list[nv_indx];
        Y                  = mu_j*Y_jm1[k][j][i][nv] + nu_j*Y_jm2[k][j][i][nv];
        Y_jm2[k][j][i][nv] = Y_jm1[k][j][i][nv];
        Y_jm1[k][j][i][nv] = Y + (1. - mu_j - nu_j)*UU_0[k][j][i][nv]
                               + absh*mu_tilde*(F_jm1[k][j][i][nv] 
                                              - a_jm1*F_0[k][j][i][nv]);
      }                                              
    } /* END DOM_LOOP  */    
    
  /* -- Put parameters of s -> s-1 (outside of dom loop) -- */
    
    b_jm2  = b_jm1;    b_jm1 = b_j;    Tjm2 = Tjm1;
    Tjm1   = Tj;       dTjm2 = dTjm1;  dTjm1 = dTj;    
    d2Tjm2 = d2Tjm1;  d2Tjm1 = d2Tj;

  /* ----------------------------------------------------------
      convert conservative variables to primitive for the next 
      s steps  (or timestep if s == s_RKC)
     ---------------------------------------------------------- */

    KDOM_LOOP(k){
    JDOM_LOOP(j){
      ConsToPrim (Y_jm1[k][j], v, IBEG, IEND, flag);
      IDOM_LOOP(i){
      for (nv = NVAR; nv--;   ){
        d->Vc[nv][k][j][i] = v[i][nv];
      }}
    }}   
  }/* s loop */
}/* func */
