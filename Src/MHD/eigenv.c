/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Wave-speeds and characteristic decomposition for the MHD equations.

  This file contains various functions containing Jacobian-related information
  such as characteristic signal speeds, eigenvalues and eigenvector 
  decomposition for the MHD module.
  
  The function MaxSignalSpeed() computes the maximum and minimum 
  characteristic signal velocity for the MHD equations.
  
  The function Eigenvalues() computes the 7 characteristic waves.
  
  The function PrimEigenvectors() calculates left and right eigenvectors
  and the corresponding eigenvalues for the \e primitive form the the
  MHD equations with an adiabatic or isothermal EoS.

  The function ConsEigenvectors() provides the characteristic decomposition
  of the convervative MHD equations.  
  
  The function PrimToChar()  compute the matrix-vector multiplcation 
  between the L matrix (containing primitive left eigenvectors) 
  and the vector v. The result containing the characteristic
  variables is stored in the vector w = L.v

  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date    Feb 23, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void MaxSignalSpeed (const State *state, double *cmin, double *cmax,
                     int beg, int end)
/*!
 * Compute the maximum and minimum characteristic velocities for the 
 * MHD equations, cmin= v - cf, cmax = v + cf
 *
 * \param [in]  state  pointer to a State structure
 * \param [out] cmin  1-D array containing the leftmost characteristic
 * \param [out] cmin  1-D array containing the rightmost characteristic
 * \param [in]  beg   starting index of computation
 * \param [in]  end   final   index of computation
 *
 *********************************************************************** */
{
  int  i;
  double gpr, Bmag2, Btmag2;
  double cf, b1, b2, b3, cw;
  double *v, *bgf, cs2;

/* ---------------------------------------------------------
   With CR, we need to solve a quartic equation.
   --------------------------------------------------------- */

  b1 = b2 = b3 = 0.0;
  for (i = beg; i <= end; i++) {

    v   = state->v[i];
    cs2 = state->a2[i];

  /* ----------------------------------------------------------
      the following are equivalent, but round-off may prevent
      a couple of tests from reproducing the same log files
     ---------------------------------------------------------- */

    #if EOS == IDEAL
    gpr = g_gamma*v[PRS];
    #else
    gpr = cs2*v[RHO];
    #endif

  /* -- get total field -- */

    EXPAND (b1 = v[BXn];  ,
            b2 = v[BXt];  ,
            b3 = v[BXb];)

    #if BACKGROUND_FIELD == YES
    bgf = state->Bbck[i] - BX1;
    EXPAND (b1 += bgf[BXn]; ,
            b2 += bgf[BXt]; ,
            b3 += bgf[BXb];)
    #endif
    Btmag2 = b2*b2 + b3*b3;
    Bmag2  = b1*b1 + Btmag2;

    cf = gpr - Bmag2;
    cf = gpr + Bmag2 + sqrt(cf*cf + 4.0*gpr*Btmag2);
    cf = sqrt(0.5*cf/v[RHO]);

    #if HALL_MHD    
    cw = state->cw[i];
  
    cmin[i] = v[VXn] - MAX(cw, cf);
    cmax[i] = v[VXn] + MAX(cw, cf);
    #else
    cmin[i] = v[VXn] - cf;
    cmax[i] = v[VXn] + cf;
    #endif
  }
}

/* ********************************************************************* */
void Eigenvalues(double **v, double *csound2, double **lambda, 
                 int beg, int end)
/*!
 * Compute eigenvalues for the MHD equations
 *
 * \param [in]  v        1-D array of primitive variables
 * \param [out] csound2  1-D array containing the square of sound speed
 * \param [out] lambda    1-D array [i][nv] containing the eigenvalues
 * \param [in]  beg   starting index of computation
 * \param [in]  end   final   index of computation
 *
 *********************************************************************** */
{
  int  i, k;
  double *q;
  double scrh0, scrh1, scrh2, scrh3, scrh4;
  double u, a, a2, ca2, cf2, cs2;
  double cs, ca, cf, b2, A2, At2;
  double tau;
  double sqrt_rho;

  for (i = beg; i <= end; i++){
    q   = v[i];

    u   = q[VXn];
    tau = 1.0/q[RHO];
    sqrt_rho = sqrt(q[RHO]);

    a2    = csound2[i];

    scrh2 = q[BXn]*q[BXn];            /* > 0 */
    scrh3 = EXPAND(0.0, + q[BXt]*q[BXt],  + q[BXb]*q[BXb]);   /* > 0 */

    b2  = scrh2 + scrh3;     /*  >0             */
    ca2 = scrh2*tau;         /*  >0  if tau >0  */
    A2  = b2*tau;            /*  >0  ''   ''    */
    At2 = scrh3*tau;

    scrh1 = a2 - A2;
    scrh0 = sqrt(scrh1*scrh1 + 4.0*a2*At2);      /*   >0   */

/*   Now get fast and slow speeds   */
    
    cf2 = 0.5*(a2 + A2 + scrh0);
    cs2 = a2*ca2/cf2;

    cf = sqrt(cf2);
    cs = sqrt(cs2);
    ca = sqrt(ca2);
    a  = sqrt(a2);

    lambda[i][KFASTM] = u - cf;
    lambda[i][KFASTP] = u + cf;

#if EOS == IDEAL 
    lambda[i][KENTRP] = u;
#endif

#ifndef GLM_MHD
  #if DIVB_CONTROL == EIGHT_WAVES
    lambda[i][KDIVB] = u;
  #else  
    lambda[i][KDIVB] = 0.0;
  #endif
#endif

#if COMPONENTS > 1
    lambda[i][KSLOWM] = u - cs;
    lambda[i][KSLOWP] = u + cs;
#endif

#if COMPONENTS == 3
    lambda[i][KALFVM] = u - ca;
    lambda[i][KALFVP] = u + ca;
#endif

#ifdef GLM_MHD
    lambda[i][KPSI_GLMM] = -glm_ch;
    lambda[i][KPSI_GLMP] =  glm_ch;
#endif
  }
}
/* ********************************************************************* */
void PrimEigenvectors(const State *state, int beg, int end)
/*!
 * Provide left and right eigenvectors and corresponding 
 * eigenvalues for the PRIMITIVE form of the MHD equations 
 * (adiabatic & isothermal cases).
 * 
 * \param [in,out]  state   Pointer to State structure
 * \param [in]      beg     starting grid index
 * \param [in]      end     final    grid index
 *
 *  
 * \note Eigenvectors \c state->Lp and \c state->Rp must be initialized to 
 *       zero \em before since only non-zero entries are treated here. 
 *
 *  Wave names and their order are defined as enumeration constants in 
 *  mod_defs.h. 
 *  Notice that, the characteristic decomposition may differ 
 *  depending on the way div.B is treated.
 * 
 *  Advection modes associated with passive scalars are simple cases 
 *  for which lambda = u (entropy mode) and l = r = (0, ... , 1, 0, ...). 
 *  For this reason they are NOT defined here and must be treated 
 *  seperately.
 *
 * \b References:
 *    - "Notes on the Eigensystem of Magnetohydrodynamics", \n
 *      P.L. Roe, D.S. Balsara,
 *      SIAM Journal on Applied Mathematics, 56, 57 (1996)
 *
 *    - "A solution adaptive upwind scheme for ideal MHD", \n
 *      K. G. Powell, P. L. Roe, and T. J. Linde,
 *      Journal of Computational Physics, 154, 284-309, (1999).
 *
 *    - "A second-order unsplit Godunov scheme for cell-centered MHD:
 *       the CTU-GLM scheme"\n
 *      Mignone \& Tzeferacos, JCP (2010) 229, 2117
 *
 *    - "ATHENA: A new code for astrophysical MHD",
 *      J. Stone, T. Gardiner, 
 *      ApJS, 178, 137 (2008)
 *
 * 
 * The derivation of the isothermal eigenvectors follows the
 * consideration given in roe.c
 *
 *********************************************************************** */
#define  sqrt_1_2  (0.70710678118654752440)
{
  int    i, k;
  double *q, a2, h;
  double **LL, **RR, *lambda;
  double scrh0, scrh1, scrh2, scrh3, scrh4;
  double u, a, ca2, cf2, cs2;
  double cs, ca, cf, b2, A2, At2;
  double tau, S;
  double alpha_f, alpha_s, beta_y, beta_z;
  double sqrt_rho;

#if BACKGROUND_FIELD == YES
  print ("! Background field does not support characteristic limiting\n");
  QUIT_PLUTO(1);
#endif

  for (i = beg; i <= end; i++){
    q      = state->v[i];
    a2     = state->a2[i];
    h      = state->h[i];
    LL     = state->Lp[i];
    RR     = state->Rp[i];
    lambda = state->lambda[i];

    u   = q[VXn];
    tau = 1.0/q[RHO];
    sqrt_rho = sqrt(q[RHO]);
  
    scrh2 = q[BXn]*q[BXn];                                 /*  Bx^2 */
    scrh3 = EXPAND(0.0, + q[BXt]*q[BXt],  + q[BXb]*q[BXb]);  /*  Bt^2 */
  
    b2  = scrh2 + scrh3;     /*  B^2 = Bx^2 + Bt^2 */
    ca2 = scrh2*tau;         /*  Bx^2/rho          */
    A2  = b2*tau;            /*  B^2/rho           */
    At2 = scrh3*tau;         /*  Bt^2/rho           */
  
    scrh1 = a2 - A2;
    scrh0 = sqrt(scrh1*scrh1 + 4.0*a2*At2);  /* sqrt( (g*p/rho-B^2/rho)^2 
                                                     + 4*g*p/rho*Bt^2/rho)  */
  /* --  Obtain fast and slow speeds -- */
      
    cf2 = 0.5*(a2 + A2 + scrh0);
    cs2 = a2*ca2/cf2;
  
    cf = sqrt(cf2);
    cs = sqrt(cs2);
    ca = sqrt(ca2);
    a  = sqrt(a2);
  
    if (cf == cs) {
      alpha_f = 1.0;
      alpha_s = 0.0;
    }else{
      scrh0   = 1.0/scrh0;
      alpha_f = (a2 - cs2)*scrh0;
      alpha_s = (cf2 - a2)*scrh0;
  
      alpha_f = MAX(0.0, alpha_f);
      alpha_s = MAX(0.0, alpha_s);
  
      alpha_f = sqrt(alpha_f);
      alpha_s = sqrt(alpha_s);
    }
  
    scrh0 = sqrt(scrh3);
    if (scrh0 > 1.e-9) {
      SELECT (                        , 
              beta_y = DSIGN(q[BXt]);  ,
              beta_y = q[BXt] / scrh0; 
              beta_z = q[BXb] / scrh0;)
    } else {
      SELECT (                  , 
              beta_y = 1.0;     ,
              beta_z = beta_y = sqrt_1_2;)
    }
  
    S = (q[BXn] >= 0.0 ? 1.0 : -1.0);
  
  /*  ------------------------------------------------------------
       define primitive right and left eigenvectors (RR and LL),
       for all of the 8 waves;
       left eigenvectors for fast & slow waves can be defined
       in terms of right eigenvectors (see page 296)  
      ------------------------------------------------------------  */
  
    /* -------------------------
        FAST WAVE,  (u - c_f)
       -------------------------  */
  
    k = KFASTM;  
    lambda[k] = u - cf;
    scrh0 = alpha_s*cs*S;
    scrh1 = alpha_s*sqrt_rho*a;
    scrh2 = 0.5 / a2;
    scrh3 = scrh2*tau;
  
    RR[RHO][k] = q[RHO]*alpha_f;        
    EXPAND(RR[VXn][k] = -cf*alpha_f;   ,
           RR[VXt][k] = scrh0*beta_y;  ,
           RR[VXb][k] = scrh0*beta_z;)
    EXPAND(                        ;  , 
           RR[BXt][k] = scrh1*beta_y;  ,
           RR[BXb][k] = scrh1*beta_z;)
  
  #if HAVE_ENERGY
  //  scrh4 = alpha_f*g_gamma*q[PRS];
    scrh4 = alpha_f*a2*q[RHO];
    RR[PRS][k] = scrh4;
  #endif
  
  #if EOS == ISOTHERMAL
    LL[k][RHO] = 0.5*alpha_f/q[RHO];
  #endif
    EXPAND(LL[k][VXn] = RR[VXn][k]*scrh2; ,
           LL[k][VXt] = RR[VXt][k]*scrh2; ,
           LL[k][VXb] = RR[VXb][k]*scrh2;) 
    EXPAND(                             ; , 
           LL[k][BXt] = RR[BXt][k]*scrh3; ,
           LL[k][BXb] = RR[BXb][k]*scrh3;)
  #if HAVE_ENERGY
    LL[k][PRS] = alpha_f*scrh3;
  #endif
  
    /* -------------------------
        FAST WAVE,  (u + c_f)
       -------------------------  */
  
    k = KFASTP; 
    lambda[k] = u + cf;
    RR[RHO][k] = RR[RHO][KFASTM];
    EXPAND(RR[VXn][k] = -RR[VXn][KFASTM];  ,
           RR[VXt][k] = -RR[VXt][KFASTM];  ,
           RR[VXb][k] = -RR[VXb][KFASTM];)
    EXPAND(                            ;   ,
           RR[BXt][k] = RR[BXt][KFASTM];   ,
           RR[BXb][k] = RR[BXb][KFASTM];)
  #if HAVE_ENERGY
    RR[PRS][k] = RR[PRS][KFASTM];
  #endif
  
  #if EOS == ISOTHERMAL
    LL[k][RHO] = LL[KFASTM][RHO];
  #endif
    EXPAND(LL[k][VXn] = -LL[KFASTM][VXn];  ,
           LL[k][VXt] = -LL[KFASTM][VXt];  ,
           LL[k][VXb] = -LL[KFASTM][VXb];) 
    EXPAND(                            ;   ,                         
           LL[k][BXt] = LL[KFASTM][BXt];   ,
           LL[k][BXb] = LL[KFASTM][BXb];)
  #if HAVE_ENERGY
    LL[k][PRS] = LL[KFASTM][PRS]; 
  #endif
  
    /* -------------------------
        entropy wave,  (u) only
        in ideal MHD
       -------------------------  */
  
  #if HAVE_ENERGY 
    k = KENTRP;
    lambda[k] = u;
    RR[RHO][k] =   1.0;
    LL[k][RHO] =   1.0; 
    LL[k][PRS] = - 1.0/a2;
  #endif
  
    /* -------------------------
          magnetic flux, (u)
       -------------------------  */
  
  #ifndef GLM_MHD
    k = KDIVB;
    #if DIVB_CONTROL == EIGHT_WAVES
    lambda[k] = u;
    RR[BXn][k] = 1.0;
    LL[k][BXn] = 1.0;
    #else  
    lambda[k] = 0.0;
    #endif
  #endif
  
  #if COMPONENTS > 1
  
    /* -------------------------
        SLOW WAVE,  (u - c_s)
       -------------------------  */
  
    k = KSLOWM;
    lambda[k] = u - cs;
    scrh0 = alpha_f*cf*S;
    scrh1 = alpha_f*sqrt_rho*a;
  
    RR[RHO][k] = q[RHO]*alpha_s;
    EXPAND(RR[VXn][k] = -cs*alpha_s;     ,
           RR[VXt][k] = -scrh0*beta_y;   ,
           RR[VXb][k] = -scrh0*beta_z;)
    EXPAND(                         ;  ,
           RR[BXt][k] = -scrh1*beta_y;  ,
           RR[BXb][k] = -scrh1*beta_z;)
    #if HAVE_ENERGY
  //  scrh4 = alpha_s*g_gamma*q[PRS]; 
    scrh4 = alpha_s*a2*q[RHO]; 
    RR[PRS][k] = scrh4;
    #endif
  
    #if EOS == ISOTHERMAL
    LL[k][RHO] = 0.5*alpha_s/q[RHO];
    #endif
    EXPAND(LL[k][VXn] = RR[VXn][k]*scrh2; ,
           LL[k][VXt] = RR[VXt][k]*scrh2; ,
           LL[k][VXb] = RR[VXb][k]*scrh2;) 
    EXPAND(                           ; ,
           LL[k][BXt] = RR[BXt][k]*scrh3; ,
           LL[k][BXb] = RR[BXb][k]*scrh3;) 
  
    #if HAVE_ENERGY
    LL[k][PRS] = alpha_s*scrh3;
    #endif
  
    /* -------------------------
        SLOW WAVE,  (u + c_s)
       -------------------------  */
  
    k = KSLOWP;
    lambda[k] = u + cs;
  
    RR[RHO][k] = RR[RHO][KSLOWM];
    EXPAND(RR[VXn][k] = -RR[VXn][KSLOWM];   ,
           RR[VXt][k] = -RR[VXt][KSLOWM];   ,
           RR[VXb][k] = -RR[VXb][KSLOWM];)
  
    EXPAND(                          ;  ,
           RR[BXt][k] = RR[BXt][KSLOWM];  ,
           RR[BXb][k] = RR[BXb][KSLOWM];)
    #if HAVE_ENERGY
    RR[PRS][k] = scrh4;
    #endif
  
    #if EOS == ISOTHERMAL
    LL[k][RHO] = LL[KSLOWM][RHO];
    #endif
    EXPAND(LL[k][VXn] = -LL[KSLOWM][VXn];   ,
           LL[k][VXt] = -LL[KSLOWM][VXt];   ,
           LL[k][VXb] = -LL[KSLOWM][VXb];) 
    EXPAND(                          ;   ,
           LL[k][BXt] = LL[KSLOWM][BXt];   ,
           LL[k][BXb] = LL[KSLOWM][BXb];) 
  
    #if HAVE_ENERGY
    LL[k][PRS] = LL[KSLOWM][PRS];
    #endif
  #endif
  
  #if COMPONENTS == 3
  
    /* -------------------------
        Alfven WAVE,  (u - c_a)
       -------------------------  */
     
    k = KALFVM;
    lambda[k] = u - ca;
    scrh2 = beta_y*sqrt_1_2;
    scrh3 = beta_z*sqrt_1_2;
  
    RR[VXt][k] = -scrh3;  
    RR[VXb][k] =  scrh2;
    RR[BXt][k] = -scrh3*sqrt_rho*S;   
    RR[BXb][k] =  scrh2*sqrt_rho*S;   
  
    LL[k][VXt] = RR[VXt][k]; 
    LL[k][VXb] = RR[VXb][k]; 
    LL[k][BXt] = RR[BXt][k]*tau;
    LL[k][BXb] = RR[BXb][k]*tau; 
  
    /* -------------------------
        Alfven WAVE,  (u + c_a)
        (-R, -L of the eigenv 
        defined by Gardiner Stone)
       -------------------------  */
  
    k = KALFVP;
    lambda[k] = u + ca;
    RR[VXt][k] =   RR[VXt][KALFVM]; 
    RR[VXb][k] =   RR[VXb][KALFVM]; 
    RR[BXt][k] = - RR[BXt][KALFVM]; 
    RR[BXb][k] = - RR[BXb][KALFVM]; 
  
    LL[k][VXt] =   LL[KALFVM][VXt];
    LL[k][VXb] =   LL[KALFVM][VXb];
    LL[k][BXt] = - LL[KALFVM][BXt];
    LL[k][BXb] = - LL[KALFVM][BXb];
  
  /* -------------------------------------------------------------------
      HotFix: when B=0 in a 2D plane and only Bz != 0, enforce zero
      jumps in BXt. This is useful with CharTracing since 2 equal and
      opposite jumps may accumulate due to round-off.
      It also perfectly consistent, since no (Bx,By) can be generated.  
     ------------------------------------------------------------------- */
  
    #if DIMENSIONS <= 2
    if (q[BXn] == 0 && q[BXt] == 0.0) for (k = 0; k < NFLX; k++) RR[BXt][k] = 0.0;
    #endif
  #endif
  
  
  #ifdef GLM_MHD
  
  /* -------------------------
      GLM wave,  -glm_ch
     -------------------------  */
  
    k = KPSI_GLMM;
    lambda[k] = -glm_ch;
    RR[BXn][k]      =  1.0;
    RR[PSI_GLM][k] = -glm_ch;
  
    LL[k][BXn]      =  0.5;
    LL[k][PSI_GLM] = -0.5/glm_ch;
  
  /* -------------------------
      GLM wave,  +glm_ch
     -------------------------  */
  
    k = KPSI_GLMP;
    lambda[k] = glm_ch;
    RR[BXn][k]      =  1.0;
    RR[PSI_GLM][k] =  glm_ch;
  
    LL[k][BXn]      =  0.5;
    LL[k][PSI_GLM] =  0.5/glm_ch;
  
  #endif
  
  /* ------------------------------------------------------------
      Verify eigenvectors consistency by
  
      1) checking that A = L.Lambda.R, where A is
         the Jacobian dF/dU
      2) verify orthonormality, L.R = R.L = I
     ------------------------------------------------------------ */
  
  #if CHECK_EIGENVECTORS == YES
  {
    int ip,jp,kp;
    double dA;
    static double **A, **ALR;
  
    if (A == NULL){
      A   = ARRAY_2D(NFLX, NFLX, double);
      ALR = ARRAY_2D(NFLX, NFLX, double);
    #if COMPONENTS != 3
      print ("! PrimEigenvectors: eigenvector check requires 3 components\n");
      return;
    #endif
    }
  
   /* --------------------------------------
       Construct the Jacobian analytically
      -------------------------------------- */

    for (ip = 0; ip < NFLX; ip++){
    for (jp = 0; jp < NFLX; jp++){
      A[ip][jp] = ALR[ip][jp] = 0.0;
    }}
  
    #if HAVE_ENERGY
  
     A[RHO][RHO] = q[VXn]; A[RHO][VXn] = q[RHO];
     A[VXn][VXn] = q[VXn]; A[VXn][BXt] =  q[BXt]*tau; A[VXn][BXb] = q[BXb]*tau; A[VXn][PRS] = tau;
     A[VXt][VXt] = q[VXn]; A[VXt][BXt] = -q[BXn]*tau; 
     A[VXb][VXb] = q[VXn]; A[VXb][BXb] = -q[BXn]*tau; 
     A[BXt][VXn] = q[BXt]; A[BXt][VXt] = -q[BXn]; A[BXt][BXn] = -q[VXt]; A[BXt][BXt] = q[VXn];
     A[BXb][VXn] = q[BXb]; A[BXb][VXb] = -q[BXn]; A[BXb][BXn] = -q[VXb]; A[BXb][BXb] = q[VXn];
     A[PRS][VXn] = a2*q[RHO]; A[PRS][PRS] =  q[VXn];
  
    #elif EOS == ISOTHERMAL
  
     A[RHO][RHO] = q[VXn] ; A[RHO][VXn] = q[RHO];
     A[VXn][VXn] = q[VXn] ; A[VXn][BXt] =  q[BXt]*tau; A[VXn][BXb] = q[BXb]*tau; 
     A[VXn][RHO] = tau*a2;
     A[VXt][VXt] = q[VXn] ; A[VXt][BXt] = -q[BXn]*tau; 
     A[VXb][VXb] = q[VXn] ; A[VXb][BXb] = -q[BXn]*tau; 
     A[BXt][VXn] = q[BXt] ; A[BXt][VXt] = -q[BXn]; A[BXt][BXn] = -q[VXt]; A[BXt][BXt] = q[VXn];
     A[BXb][VXn] = q[BXb] ; A[BXb][VXb] = -q[BXn]; A[BXb][BXn] = -q[VXb]; A[BXb][BXb] = q[VXn];
  
    #endif
  
    #ifdef GLM_MHD
     A[BXn][PSI_GLM] = 1.0;
     A[PSI_GLM][BXn] = glm_ch*glm_ch;
    #endif
  
    for (ip = 0; ip < NFLX; ip++){
    for (jp = 0; jp < NFLX; jp++){
      ALR[ip][jp] = 0.0;
      for (kp = 0; kp < NFLX; kp++){
        ALR[ip][jp] += RR[ip][kp]*lambda[kp]*LL[kp][jp];
      }
    }}
  
    for (ip = 0; ip < NFLX; ip++){
    for (jp = 0; jp < NFLX; jp++){
  
    /* --------------------------------------------------------
       NOTE: if the standard 7-wave formulation is adopted, 
             the column and the row corresponding to B(normal) 
             do not exist. 
             However, PLUTO uses a primitive matrix with 8 
             entries  and the B(normal) column contains two 
             entries (-vy, -vz) which cannot be recovered using 
             a 7x7 system.
       -------------------------------------------------------- */    
  
      if (jp == BXn) continue;
      dA = ALR[ip][jp] - A[ip][jp];
      if (fabs(dA) > 1.e-8){
        print ("! PrimEigenvectors: eigenvectors not consistent\n");
        print ("! g_dir = %d\n",g_dir);
        print ("! A[%d][%d] = %16.9e, R.Lambda.L[%d][%d] = %16.9e\n",
                  ip,jp, A[ip][jp], ip,jp,ALR[ip][jp]);
        print ("\n\n A = \n");
        ShowMatrix(A, NFLX, 1.e-8);
        print ("\n\n R.Lambda.L = \n");
        ShowMatrix(ALR, NFLX, 1.e-8);
        QUIT_PLUTO(1);
      }
    }}  
  
  /* -- check orthornomality -- */
  
    for (ip = 0; ip < NFLX; ip++){
    for (jp = 0; jp < NFLX; jp++){
      #if (DIVB_CONTROL == NO) || (DIVB_CONTROL == CONSTRAINED_TRANSPORT)
       if (ip == KDIVB || jp == KDIVB) continue;
      #endif
      a = 0.0;
      for (kp = 0; kp < NFLX; kp++) a += LL[ip][kp]*RR[kp][jp];
      if ( (ip == jp && fabs(a-1.0) > 1.e-8) ||
           (ip != jp && fabs(a)>1.e-8) ) {
        print ("! PrimEigenvectors: Eigenvectors not orthogonal\n");
        print ("!   i,j = %d, %d  %12.6e \n",ip,jp,a);
        print ("!   g_dir: %d\n",g_dir);
        QUIT_PLUTO(1);
      }
    }}
  } /* End block CHECK_EIGENVECTORS */
  #endif
  }  /* end loop i=beg, end */
}
#undef sqrt_1_2

/* ********************************************************************* */
void ConsEigenvectors (double *u, double *v, double a2, 
                       double **Lc, double **Rc, double *lambda)
/*!
 * Provide conservative eigenvectors for MHD equations.
 * 
 * \param [in]   u   array of conservative variables
 * \param [in]   v   array of primitive variables
 * \param [in]  a2   square of sound speed
 * \param [out] Lc   left conservative eigenvectors
 * \param [out] Rc   right  conservative eigenvectors
 * \param [out] lambda eigenvalues
 *
 * \b Reference
 *    - "High-order conservative finite difference GLM-MHD schemes for
 *       cell-centered MHD"\n
 *       Mignone, Tzeferacos \& Bodo, JCP (2010) 229, 5896
 *    - "A High-order WENO Finite Difference Scheme for the Equations
 *       of Ideal MHD"\n
 *       Jiang,Wu, JCP 150,561 (1999)
 *
 *  With the following corrections:
 *
 *  The components (By, Bz) of  L_{1,7}  page 571 should be 
 *  +\sqrt{rho}a\alpha_s instead of -.
 *  Also, since Bx is not considered as a parameter, one must also add a 
 *  component in the fast and slow waves.
 *  This can be seen by forming the conservative eigenvectors from the 
 *  primitive ones, see the paper from Powell, Roe Linde.
 *
 *  The Lc_{1,3,5,7}[Bx] components, according to Serna 09
 *  are -(1-g_gamma)*a_{f,s}*Bx and not (1-g_gamma)*a_{f,s}*Bx. Both are 
 *  orthonormal though. She is WRONG!
 *  -Petros-
 *********************************************************************** */
{
  int    nv, i,j,k;
  double rho;
  double vx, vy, vz;
  double Bx, By, Bz;
  double bx, by, bz;
  double beta_y, beta_z;
  double one_sqrho;
  double cf2, cs2, ca2, v2;
  double cf, cs, ca, a;
  double alpha_f, alpha_s;
  double g1, g2, tau, Gf, Ga, Gs;
  double scrh0, scrh1, S, bt2, b2, Btmag;
  double one_gmm;
  double LBX = 0.0;  /* setting LBX to 1.0 will include Powell's */
                     /* eight wave. Set it to 0 to recover the   */ 
                     /* standard 7 wave formulation              */  

#if BACKGROUND_FIELD == YES
  print ("! ConsEigenvectors: Background field does not support\n");
  print ("                    characteristic limiting.\n");
  QUIT_PLUTO(1);
#endif

#if EOS == PVTE_LAW
  print( "! ConsEigenvectors: cannot be used presently with PVTE_LAW EoS\n");
  QUIT_PLUTO(1);  
#endif

/* --------------------------------------------------------------
    If eigenvector check is required, we make sure  that 
    U = U(V) is exactly converted from V. 
    This is not the case, with the present version of the code, 
    with FD scheme where V = 0.5*(VL + VR) and U = 0.5*(UL + UR) 
    hence U != U(V).
   --------------------------------------------------------------- */
    
  #if CHECK_EIGENVECTORS == YES
   u[RHO] = v[RHO];
   u[MX1] = v[RHO]*v[VX1];
   u[MX2] = v[RHO]*v[VX2];
   u[MX3] = v[RHO]*v[VX3];
   u[BX1] = v[BX1];
   u[BX2] = v[BX2];
   u[BX3] = v[BX3];
   #if EOS == IDEAL
    u[ENG] =   v[PRS]/(g_gamma-1.0)
            + 0.5*v[RHO]*(v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3])
            + 0.5*(v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3]);
   #endif
  #endif

/* ----------------------------------------------------------
            compute speed and eigenvalues  
   ---------------------------------------------------------- */

  #if EOS == BAROTROPIC
   print ("! ConsEigenvectors: not defined for barotropic MHD\n");
   QUIT_PLUTO(1);
  #endif

  rho = v[RHO];

  EXPAND (vx = v[VXn];  ,
          vy = v[VXt];  ,
          vz = v[VXb];)

  EXPAND (Bx = v[BXn];  ,
          By = v[BXt];  ,
          Bz = v[BXb];)
  
  one_sqrho = 1.0/sqrt(rho);
  
  S   = (Bx >= 0.0 ? 1.0 : -1.0);
 
  EXPAND(bx = Bx*one_sqrho;  ,
         by = By*one_sqrho;  ,
         bz = Bz*one_sqrho;)
    
  bt2   = EXPAND(0.0  , + by*by, + bz*bz);
  b2    = bx*bx + bt2;
  Btmag = sqrt(bt2*rho);
  v2    = EXPAND(vx*vx , + vy*vy, + vz*vz);

/* ------------------------------------------------------------
    Compute fast and slow magnetosonic speeds.

    The following expression appearing in the definitions
    of the fast magnetosonic speed 
    
     (a^2 - b^2)^2 + 4*a^2*bt^2 = (a^2 + b^2)^2 - 4*a^2*bx^2

    is always positive and avoids round-off errors.
   ------------------------------------------------------------ */
        
  scrh0 = a2 - b2;
  ca2   = bx*bx;
  scrh0 = scrh0*scrh0 + 4.0*bt2*a2;    
  scrh0 = sqrt(scrh0);    

  cf2 = 0.5*(a2 + b2 + scrh0); 
  cs2 = a2*ca2/cf2;   /* -- same as 0.5*(a2 + b2 - scrh0) -- */
    
  cf = sqrt(cf2);
  cs = sqrt(cs2);
  ca = sqrt(ca2);
  a  = sqrt(a2); 

  if (Btmag > 1.e-9) {
    SELECT(                     , 
           beta_y = DSIGN(By);  ,
           beta_y = By/Btmag; 
           beta_z = Bz/Btmag;)
  } else {
    SELECT(                       , 
           beta_y = 1.0;          ,
           beta_z = beta_y = sqrt(0.5);)
  }
    
  if (cf == cs) {
    alpha_f = 1.0;
    alpha_s = 0.0;
  }else if (a <= cs) {
    alpha_f = 0.0;
    alpha_s = 1.0;
  }else if (cf <= a){
    alpha_f = 1.0;
    alpha_s = 0.0;
  }else{ 
    scrh0   = 1.0/(cf2 - cs2);
    alpha_f = (a2  - cs2)*scrh0;
    alpha_s = (cf2 -  a2)*scrh0;
    alpha_f = MAX(0.0, alpha_f);
    alpha_s = MAX(0.0, alpha_s);
    alpha_f = sqrt(alpha_f);
    alpha_s = sqrt(alpha_s);
  }

/* --------------------------------------------------------
    Compute non-zero entries of conservative
    eigenvectors (Rc, Lc).
   -------------------------------------------------------- */
 
  #if EOS == IDEAL
   one_gmm = 1.0 - g_gamma;
   g1  = 0.5*(g_gamma - 1.0);
   g2  = (g_gamma - 2.0)/(g_gamma - 1.0);
   tau = (g_gamma - 1.0)/a2;
  #elif EOS == ISOTHERMAL 
   one_gmm = g1 = tau = 0.0;
  #endif   

  scrh0 = EXPAND(0.0, + beta_y*vy, + beta_z*vz);
  Gf = alpha_f*cf*vx - alpha_s*cs*S*scrh0;
  Ga = EXPAND(0.0,    , + S*(beta_z*vy - beta_y*vz));
  Gs = alpha_s*cs*vx + alpha_f*cf*S*scrh0;
  
/* -----------------------
    FAST WAVE  (u - c_f) 
   ----------------------- */

  k = KFASTM;
  lambda[k] = vx - cf;

  scrh0 = alpha_s*cs*S;
  scrh1 = alpha_s*a*one_sqrho;
  Rc[RHO][k] = alpha_f;
  EXPAND( Rc[MXn][k] = alpha_f*lambda[k];          ,
          Rc[MXt][k] = alpha_f*vy + scrh0*beta_y;  ,
          Rc[MXb][k] = alpha_f*vz + scrh0*beta_z; ) 
  EXPAND(                         ;  ,  
          Rc[BXt][k] = scrh1*beta_y;  ,
          Rc[BXb][k] = scrh1*beta_z; )
  #if EOS == IDEAL
   Rc[ENG][k] = alpha_f*(0.5*v2 + cf2 - g2*a2) - Gf;
  #endif

  Lc[k][RHO] = (g1*alpha_f*v2 + Gf)*0.5/a2; 
  #if EOS == ISOTHERMAL
   Lc[k][RHO] += alpha_f*0.5;
  #endif
  EXPAND( Lc[k][MXn] = (one_gmm*alpha_f*vx - alpha_f*cf)  *0.5/a2;  ,
          Lc[k][MXt] = (one_gmm*alpha_f*vy + scrh0*beta_y)*0.5/a2;  ,
          Lc[k][MXb] = (one_gmm*alpha_f*vz + scrh0*beta_z)*0.5/a2;) 
  EXPAND( Lc[k][BXn] =  LBX*one_gmm*alpha_f*Bx*0.5/a2;                  , 
          Lc[k][BXt] = (one_gmm*alpha_f*By + scrh1*rho*beta_y)*0.5/a2;  ,
          Lc[k][BXb] = (one_gmm*alpha_f*Bz + scrh1*rho*beta_z)*0.5/a2; )
  #if EOS == IDEAL
   Lc[k][ENG] = alpha_f*(g_gamma - 1.0)*0.5/a2;
  #endif
  
  /* -----------------------
      FAST WAVE  (u + c_f) 
     ----------------------- */

  k = KFASTP;
  lambda[k] = vx + cf;

  Rc[RHO][k] = alpha_f;
  EXPAND( Rc[MXn][k] = alpha_f*lambda[k];          ,
          Rc[MXt][k] = alpha_f*vy - scrh0*beta_y;  ,
          Rc[MXb][k] = alpha_f*vz - scrh0*beta_z; ) 
  EXPAND(                         ;  , 
          Rc[BXt][k] = scrh1*beta_y;  ,
          Rc[BXb][k] = scrh1*beta_z; )
  #if EOS == IDEAL
   Rc[ENG][k] = alpha_f*(0.5*v2 + cf2 - g2*a2) + Gf;
  #endif

  Lc[k][RHO] = (g1*alpha_f*v2 - Gf)*0.5/a2; 
  #if EOS == ISOTHERMAL
   Lc[k][RHO] += alpha_f*0.5;
  #endif
  EXPAND( Lc[k][MXn] = (one_gmm*alpha_f*vx + alpha_f*cf)  *0.5/a2;  ,
          Lc[k][MXt] = (one_gmm*alpha_f*vy - scrh0*beta_y)*0.5/a2;  ,
          Lc[k][MXb] = (one_gmm*alpha_f*vz - scrh0*beta_z)*0.5/a2;) 
  EXPAND( Lc[k][BXn] =  LBX*one_gmm*alpha_f*Bx*0.5/a2;     , 
          Lc[k][BXt] = (one_gmm*alpha_f*By + sqrt(rho)*a*alpha_s*beta_y)*0.5/a2;  ,
          Lc[k][BXb] = (one_gmm*alpha_f*Bz + sqrt(rho)*a*alpha_s*beta_z)*0.5/a2; )
  #if EOS == IDEAL
   Lc[k][ENG] = alpha_f*(g_gamma - 1.0)*0.5/a2;
  #endif
  
  /* -----------------------
      Entropy wave  (u) 
     ----------------------- */

  #if EOS == IDEAL
   k = KENTRP;
   lambda[k] = vx;

   Rc[RHO][k] = 1.0;
   EXPAND( Rc[MXn][k] = vx; ,
           Rc[MXt][k] = vy; ,
           Rc[MXb][k] = vz; )
   Rc[ENG][k] = 0.5*v2;

   Lc[k][RHO] = 1.0 - 0.5*tau*v2;
   EXPAND(Lc[k][VXn] = tau*vx;  ,
          Lc[k][VXt] = tau*vy;  ,
          Lc[k][VXb] = tau*vz;)
   EXPAND(Lc[k][BXn] = LBX*tau*Bx;  ,
          Lc[k][BXt] = tau*By;  ,
          Lc[k][BXb] = tau*Bz;)
   Lc[k][ENG] = -tau;
  #endif

  /* --------------------------------
            div.B wave  (u) 
     -------------------------------- */

  #ifndef GLM_MHD
   k = KDIVB;
   #if DIVB_CONTROL == EIGHT_WAVES
    lambda[k] = vx;
    Rc[BXn][k] = 1.0;
    Lc[k][BXn] = 1.0;
   #else
    lambda[k] = 0.0;
    Rc[BXn][k] = Lc[k][BXn] = 0.0;
   #endif
  #endif
    
  #if COMPONENTS > 1    

 /* -----------------------
     SLOW WAVE  (u - c_s) 
    ----------------------- */

   k = KSLOWM;
   lambda[k] = vx - cs;
   scrh0 = alpha_f*cf*S;

   Rc[RHO][k] = alpha_s;
   EXPAND( Rc[MXn][k] = alpha_s*lambda[k];          ,
           Rc[MXt][k] = alpha_s*vy - scrh0*beta_y;  ,
           Rc[MXb][k] = alpha_s*vz - scrh0*beta_z; ) 
   EXPAND(                                         ;  , 
           Rc[BXt][k] = - alpha_f*a*beta_y*one_sqrho;  ,
           Rc[BXb][k] = - alpha_f*a*beta_z*one_sqrho; )

   #if EOS == IDEAL
    Rc[ENG][k] = alpha_s*(0.5*v2 + cs2 - g2*a2) - Gs;
   #endif

   Lc[k][RHO] = (g1*alpha_s*v2 + Gs)*0.5/a2; 
   
   #if EOS == ISOTHERMAL
    Lc[k][RHO] += alpha_s*0.5;
   #endif
   
   EXPAND( Lc[k][MXn] = (one_gmm*alpha_s*vx - alpha_s*cs)  *0.5/a2;  ,
           Lc[k][MXt] = (one_gmm*alpha_s*vy - scrh0*beta_y)*0.5/a2;  ,
           Lc[k][MXb] = (one_gmm*alpha_s*vz - scrh0*beta_z)*0.5/a2;) 
   EXPAND( Lc[k][BXn] = LBX*one_gmm*alpha_s*Bx*0.5/a2;                                ,                                
           Lc[k][BXt] = (one_gmm*alpha_s*By - sqrt(rho)*a*alpha_f*beta_y)*0.5/a2;  ,
           Lc[k][BXb] = (one_gmm*alpha_s*Bz - sqrt(rho)*a*alpha_f*beta_z)*0.5/a2; )
   #if EOS == IDEAL
    Lc[k][ENG] = alpha_s*(g_gamma - 1.0)*0.5/a2;
   #endif
   
   /* -----------------------
       SLOW WAVE  (u + c_s) 
      ----------------------- */

   k = KSLOWP;
   lambda[k] = vx + cs; 

   Rc[RHO][k] = alpha_s;
   EXPAND( Rc[MXn][k] = alpha_s*lambda[k];          ,
           Rc[MXt][k] = alpha_s*vy + scrh0*beta_y;  ,
           Rc[MXb][k] = alpha_s*vz + scrh0*beta_z; ) 
   EXPAND(                                         ;  , 
           Rc[BXt][k] = - alpha_f*a*beta_y*one_sqrho;  ,
           Rc[BXb][k] = - alpha_f*a*beta_z*one_sqrho; )

   #if EOS == IDEAL
    Rc[ENG][k] = alpha_s*(0.5*v2 + cs2 - g2*a2) + Gs;
   #endif

   Lc[k][RHO] = (g1*alpha_s*v2 - Gs)*0.5/a2; 

   #if EOS == ISOTHERMAL
    Lc[k][RHO] += alpha_s*0.5;
   #endif
   
   EXPAND( Lc[k][MXn] = (one_gmm*alpha_s*vx + alpha_s*cs)  *0.5/a2;  ,
           Lc[k][MXt] = (one_gmm*alpha_s*vy + scrh0*beta_y)*0.5/a2;  ,
           Lc[k][MXb] = (one_gmm*alpha_s*vz + scrh0*beta_z)*0.5/a2;) 
   EXPAND( Lc[k][BXn] = LBX*one_gmm*alpha_s*Bx*0.5/a2;                                 ,               
           Lc[k][BXt] = (one_gmm*alpha_s*By - sqrt(rho)*a*alpha_f*beta_y)*0.5/a2;  ,
           Lc[k][BXb] = (one_gmm*alpha_s*Bz - sqrt(rho)*a*alpha_f*beta_z)*0.5/a2; )
   #if EOS == IDEAL
    Lc[k][ENG] = alpha_s*(g_gamma - 1.0)*0.5/a2;
   #endif
  #endif

  #if COMPONENTS == 3

 /* ------------------------
     Alfven WAVE  (u - c_a) 
    ------------------------ */

   k = KALFVM;
   lambda[k] = vx - ca;

   Rc[MXt][k] = - beta_z*S;  
   Rc[MXb][k] = + beta_y*S;
   Rc[BXt][k] = - beta_z*one_sqrho;
   Rc[BXb][k] =   beta_y*one_sqrho;
   #if EOS == IDEAL
    Rc[ENG][k] = - Ga;
   #endif

   Lc[k][RHO] = 0.5*Ga;
   Lc[k][MXt] = - 0.5*beta_z*S;
   Lc[k][MXb] =   0.5*beta_y*S;
   Lc[k][BXt] = - 0.5*sqrt(rho)*beta_z;
   Lc[k][BXb] =   0.5*sqrt(rho)*beta_y;

 /* -----------------------
     Alfven WAVE  (u + c_a) 
    ----------------------- */

   k = KALFVP;
   lambda[k] = vx + ca;

   Rc[MXt][k] = Rc[MXt][KALFVM];  
   Rc[MXb][k] = Rc[MXb][KALFVM];
   Rc[BXt][k] = - Rc[BXt][KALFVM];   
   Rc[BXb][k] = - Rc[BXb][KALFVM];
   #if EOS == IDEAL
    Rc[ENG][k] = Rc[ENG][KALFVM];
   #endif

   Lc[k][RHO] =   Lc[KALFVM][RHO];
   Lc[k][MXt] =   Lc[KALFVM][MXt];
   Lc[k][MXb] =   Lc[KALFVM][MXb];
   Lc[k][BXt] = - Lc[KALFVM][BXt];
   Lc[k][BXb] = - Lc[KALFVM][BXb];
   #if EOS == IDEAL
    Lc[k][ENG] =   Lc[KALFVM][ENG];
   #endif
  #endif

  #ifdef GLM_MHD

  /* -------------------------
      GLM wave,  -glm_ch
     -------------------------  */

   k = KPSI_GLMM;
   lambda[k] = -glm_ch;
   Rc[BXn][k]     =  1.0;
   Rc[PSI_GLM][k] = -glm_ch;

   Lc[k][BXn]     =  0.5;
   Lc[k][PSI_GLM] = -0.5/glm_ch;

  /* -------------------------
      GLM wave,  +glm_ch
     -------------------------  */

   k = KPSI_GLMP;
   lambda[k] = glm_ch;
   Rc[BXn][k]     = 1.0;
   Rc[PSI_GLM][k] = glm_ch;

   Lc[k][BXn]     = 0.5;
   Lc[k][PSI_GLM] = 0.5/glm_ch;

  #endif

/* --------------------------------------------------------------
    Verify eigenvectors consistency by

    1) checking that A = L.Lambda.R, where A is 
       the Jacobian dF/dU
    2) verify orthonormality, L.R = R.L = I

    IMPORTANT: condition 1) can be satisfied only if U and
    V are consistently defined, that is, U = U(V). 
    This is why we perform an extra conversion at the 
    beginning of this function only when this check is enabled.
   -------------------------------------------------------------- */

#if CHECK_EIGENVECTORS == YES
{
  static double **A, **ALR;
  double dA, vel2, Bmag2, vB;

  if (A == NULL){
    A   = ARRAY_2D(NFLX, NFLX, double);
    ALR = ARRAY_2D(NFLX, NFLX, double);
    #if COMPONENTS != 3
     print ("! ConsEigenvectors: eigenvector check requires 3 components\n");
    #endif
  }
  #if COMPONENTS != 3
   return;
  #endif

 /* --------------------------------------
     Construct the Jacobian analytically
    -------------------------------------- */

  for (i = 0; i < NFLX; i++){
  for (j = 0; j < NFLX; j++){
    A[i][j] = ALR[i][j] = 0.0;
  }}

  vel2  = v[VXn]*v[VXn] + v[VXt]*v[VXt] + v[VXb]*v[VXb];
  Bmag2 = v[BXn]*v[BXn] + v[BXt]*v[BXt] + v[BXb]*v[BXb];
  vB    = v[VXn]*v[BXn] + v[VXt]*v[BXt] + v[VXb]*v[BXb];
  #if EOS == IDEAL
   A[RHO][MXn] = 1.0;

   A[MXn][RHO] = -v[VXn]*v[VXn] + (g_gamma-1.0)*vel2*0.5;  
   A[MXn][MXn] = (3.0 - g_gamma)*v[VXn];
   A[MXn][MXt] = (1.0 - g_gamma)*v[VXt];
   A[MXn][MXb] = (1.0 - g_gamma)*v[VXb];
   A[MXn][BXt] = (2.0 - g_gamma)*v[BXt];
   A[MXn][BXb] = (2.0 - g_gamma)*v[BXb];
   A[MXn][ENG] = g_gamma - 1.0;

   A[MXt][RHO] = - v[VXn]*v[VXt];
   A[MXt][MXn] =   v[VXt];
   A[MXt][MXt] =   v[VXn];
   A[MXt][BXt] = - v[BXn];

   A[MXb][RHO] = - v[VXn]*v[VXb];
   A[MXb][MXn] =   v[VXb];
   A[MXb][MXb] =   v[VXn];
   A[MXb][BXb] = - v[BXn];

   #ifdef GLM_MHD
    A[BXn][PSI_GLM] = 1.0;
   #endif

   A[BXt][RHO] =  (v[BXn]*v[VXt] - v[BXt]*v[VXn])/v[RHO];
   A[BXt][MXn] =  v[BXt]/v[RHO];
   A[BXt][MXt] = -v[BXn]/v[RHO];
   A[BXt][BXt] =  v[VXn];
   
   A[BXb][RHO] =  (v[BXn]*v[VXb] - v[BXb]*v[VXn])/v[RHO];
   A[BXb][MXn] =  v[BXb]/v[RHO];
   A[BXb][MXb] = -v[BXn]/v[RHO];
   A[BXb][BXb] =  v[VXn];

   A[ENG][RHO] =   0.5*(g_gamma - 1.0)*vel2*v[VXn] 
               - (u[ENG] + v[PRS] + 0.5*Bmag2)/v[RHO]*v[VXn] + vB*v[BXn]/v[RHO];

   A[ENG][MXn] =   (1.0 - g_gamma)*v[VXn]*v[VXn] 
               + (u[ENG] + v[PRS] + 0.5*Bmag2)/v[RHO] - v[BXn]*v[BXn]/v[RHO];

   A[ENG][MXt] = (1.0 - g_gamma)*v[VXn]*v[VXt] - v[BXt]*v[BXn]/v[RHO];
   A[ENG][MXb] = (1.0 - g_gamma)*v[VXn]*v[VXb] - v[BXb]*v[BXn]/v[RHO];

   A[ENG][BXt] = (2.0 - g_gamma)*v[BXt]*v[VXn] - v[VXt]*v[BXn];
   A[ENG][BXb] = (2.0 - g_gamma)*v[BXb]*v[VXn] - v[VXb]*v[BXn];
   A[ENG][ENG] = g_gamma*v[VXn];
   #ifdef PSI_GLM 
    A[PSI_GLM][BXn] = glm_ch;
   #endif

  #elif EOS == ISOTHERMAL
   A[RHO][MXn] = 1.0;

   A[MXn][RHO] = -v[VXn]*v[VXn] + a2;
   A[MXn][MXn] = 2.0*v[VXn];
   A[MXn][BXt] =     v[BXt];
   A[MXn][BXb] =     v[BXb];

   A[MXt][RHO] = - v[VXn]*v[VXt];
   A[MXt][MXn] =   v[VXt];
   A[MXt][MXt] =   v[VXn];
   A[MXt][BXt] = - v[BXn];

   A[MXb][RHO] = - v[VXn]*v[VXb];
   A[MXb][MXn] =   v[VXb];
   A[MXb][MXb] =   v[VXn];
   A[MXb][BXb] = - v[BXn];

   #ifdef GLM_MHD
    A[BXn][PSI_GLM] = 1.0;
   #endif

   A[BXt][RHO] =  (v[BXn]*v[VXt] - v[BXt]*v[VXn])/v[RHO];
   A[BXt][MXn] =  v[BXt]/v[RHO];
   A[BXt][MXt] = -v[BXn]/v[RHO];
   A[BXt][BXt] =  v[VXn];
   
   A[BXb][RHO] =  (v[BXn]*v[VXb] - v[BXb]*v[VXn])/v[RHO];
   A[BXb][MXn] =  v[BXb]/v[RHO];
   A[BXb][MXb] = -v[BXn]/v[RHO];
   A[BXb][BXb] =  v[VXn];

   #ifdef PSI_GLM 
    A[PSI_GLM][BXn] = glm_ch;
   #endif

  #endif

  for (i = 0; i < NFLX; i++){
  for (j = 0; j < NFLX; j++){
    ALR[i][j] = 0.0;
    for (k = 0; k < NFLX; k++){
      ALR[i][j] += Rc[i][k]*lambda[k]*Lc[k][j];
    }
  }}

  for (i = 0; i < NFLX; i++){
  for (j = 0; j < NFLX; j++){
    if (j == BXn) continue;
    if (fabs(ALR[i][j] - A[i][j]) > 1.e-6){
      print ("! ConsEigenvectors: eigenvectors not consistent\n");
      print ("! g_dir = %d\n",g_dir);
      print ("! A[%d][%d] = %16.9e, R.Lambda.L[%d][%d] = %16.9e\n",
                i,j, A[i][j], i,j,ALR[i][j]);
      print ("\n\n A = \n");   ShowMatrix(A, NFLX, 1.e-8);
      print ("\n\n R.Lambda.L = \n"); ShowMatrix(ALR, NFLX, 1.e-8);
      QUIT_PLUTO(1);
    }
  }}

/* -- check orthornomality -- */

  for (i = 0; i < NFLX; i++){
  for (j = 0; j < NFLX; j++){
    #if (DIVB_CONTROL == NO) || (DIVB_CONTROL == CONSTRAINED_TRANSPORT)
     if (i == KDIVB || j == KDIVB) continue;
    #endif
    a = 0.0;
    for (k = 0; k < NFLX; k++) a += Lc[i][k]*Rc[k][j];
    if ( (i == j && fabs(a-1.0) > 1.e-8) || 
         (i != j && fabs(a)>1.e-8) ) {
      print ("! ConsEigenvectors: Eigenvectors not orthogonal\n");
      print ("!   i,j = %d, %d  %12.6e \n",i,j,a);
      print ("!   g_dir: %d\n",g_dir);
      QUIT_PLUTO(1);
    }
  }}
}
#endif
}

/* ********************************************************************* */
void PrimToChar (double **Lp, double *v, double *w)
/*!
 *  Compute the matrix-vector multiplcation between the
 *  the L matrix (containing primitive left eigenvectors) 
 *  and the vector v. The result containing the characteristic
 *  variables is stored in the vector w = L.v
 *
 *  For efficiency purpose, multiplication is done 
 *  explicitly, so that only nonzero entries
 *  of the left primitive eigenvectors are considered.
 *  
 * \param [in]  Lp   Left eigenvectors
 * \param [in]   v   (difference of) primitive variables 
 * \param [out]  w   (difference of) characteristic variables
 *
 *********************************************************************** */
{
  int    k;
  double wv, wB, *L;

/* -- fast waves -- */

  L = Lp[KFASTM];
  wv = EXPAND(L[VXn]*v[VXn], + L[VXt]*v[VXt], + L[VXb]*v[VXb]); 
  #if HAVE_ENERGY
   wB = EXPAND(L[PRS]*v[PRS], + L[BXt]*v[BXt], + L[BXb]*v[BXb]); 
  #elif EOS == ISOTHERMAL
   wB = EXPAND(L[RHO]*v[RHO], + L[BXt]*v[BXt], + L[BXb]*v[BXb]); 
  #endif
  w[KFASTM] =  wv + wB;
  w[KFASTP] = -wv + wB;

/* -- entropy -- */

  #if HAVE_ENERGY
   L = Lp[KENTRP];
   w[KENTRP] = L[RHO]*v[RHO] + L[PRS]*v[PRS];
  #endif

  #ifndef GLM_MHD
   L = Lp[KDIVB];
   #if DIVB_CONTROL == EIGHT_WAVES
    w[KDIVB] = v[BXn];
   #else
    w[KDIVB] = 0.0;
   #endif
  #endif

  #if COMPONENTS > 1
   L = Lp[KSLOWM];
   wv = EXPAND(L[VXn]*v[VXn], + L[VXt]*v[VXt], + L[VXb]*v[VXb]); 

   #if HAVE_ENERGY
    wB = EXPAND(L[PRS]*v[PRS], + L[BXt]*v[BXt], + L[BXb]*v[BXb]); 
   #elif EOS == ISOTHERMAL
    wB = EXPAND(L[RHO]*v[RHO], + L[BXt]*v[BXt], + L[BXb]*v[BXb]); 
   #endif

   w[KSLOWM] =  wv + wB;
   w[KSLOWP] = -wv + wB;

   #if COMPONENTS == 3
    L = Lp[KALFVM];
    wv = L[VXt]*v[VXt] + L[VXb]*v[VXb]; 
    wB = L[BXt]*v[BXt] + L[BXb]*v[BXb]; 
    w[KALFVM] = wv + wB;
    w[KALFVP] = wv - wB;
   #endif
  #endif
  
  #ifdef GLM_MHD
   L = Lp[KPSI_GLMP];
   w[KPSI_GLMM] = L[BXn]*v[BXn] - L[PSI_GLM]*v[PSI_GLM];
   w[KPSI_GLMP] = L[BXn]*v[BXn] + L[PSI_GLM]*v[PSI_GLM];
  #endif

/* ----------------------------------------------- 
     For passive scalars, the characteristic 
     variable is equal to the primitive one, 
     since  l = r = (0,..., 1 , 0 ,....)
   ----------------------------------------------- */		   

#if NSCL > 0 
  NSCL_LOOP(k) w[k] = v[k];
#endif   

/* -------------------------------------------------
    verify that the previous simplified expressions 
    are indeed  w = Lp.v
   ------------------------------------------------- */

#if CHECK_EIGENVECTORS == YES
{
  int    i;
  double w2[NVAR];

  for (k = 0; k < NVAR; k++){
    w2[k] = 0.0;
    for (i = 0; i < NVAR; i++) w2[k] += Lp[k][i]*v[i];
    if (fabs(w[k] - w2[k]) > 1.e-8){
      printf ("! PrimToChar: projection not correct, k = %d\n",k);
      QUIT_PLUTO(1);
    }
  }
}
#endif
}
