/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Wave-speeds and characteristic decomposition for the HD equations.

  This file contains various functions containing Jacobian-related information
  such as characteristic signal speeds, eigenvalues and eigenvector
  decomposition for the HD module.

  The function MaxSignalSpeed() computes the maximum and minimum
  characteristic signal velocity for the HD equations.

  The function Eigenvalues() computes the 3 characteristic wave speed

  The function PrimEigenvectors() calculates left and right eigenvectors
  and the corresponding eigenvalues for the \e primitive form the the
  HD equations with an adiabatic, general convex or isothermal EoS.

  The function ConsEigenvectors() provides the characteristic decomposition
  of the convervative HD equations.

  The function PrimToChar()  compute the matrix-vector multiplcation
  between the L matrix (containing primitive left eigenvectors)
  and the vector v. The result containing the characteristic
  variables is stored in the vector w = L.v

  \author A. Mignone (mignone@ph.unito.it)
  \date   Feb 13, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void MaxSignalSpeed (const State *state, double *cmin, double *cmax,
                     int beg, int end)
/*!
 * Compute the maximum and minimum characteristic velocities for the
 * HD equations from the vector of primitive variables state->v.
 *
 * \param [in]  state  pointer to a State structure
 * \param [out] cmin  1-D array containing the leftmost characteristic
 * \param [out] cmin  1-D array containing the rightmost characteristic
 * \param [in]  beg   starting index of computation
 * \param [in]  end   final   index of computation
 *
 *********************************************************************** */
{
  int    i;
  double a;
  double **v  = state->v;
  double *cs2 = state->a2;

  for (i = beg; i <= end; i++) {
    #if HAVE_ENERGY
/*    a = sqrt(g_gamma*v[i][PRS]/v[i][RHO]);  */
     a = sqrt(cs2[i]);
    #elif EOS == ISOTHERMAL
/*     a = g_isoSoundSpeed;   */
     a = sqrt(cs2[i]);
    #endif

    cmin[i] = v[i][VXn] - a;
    cmax[i] = v[i][VXn] + a;
  }
}

/* ********************************************************************* */
void Eigenvalues(double **v, double *csound2, double **lambda, int beg, int end)
/*!
 * Compute eigenvalues for the HD equations
 *
 * \param [in]  v        1-D array of primitive variables
 * \param [out] csound2  1-D array containing the square of sound speed
 * \param [out] lambda   1-D array [i][nv] containing the eigenvalues
 * \param [in]  beg      starting index of computation
 * \param [in]  end      final    index of computation
 *
 *********************************************************************** */
{
  int    i, k;
  double cs;

  for (i = beg; i <= end; i++){
    cs = sqrt(csound2[i]);
    lambda[i][0] = v[i][VXn] - cs;  
    lambda[i][1] = v[i][VXn] + cs;  
    for (k = 2; k < NFLX; k++) lambda[i][k] = v[i][VXn];
  } 
}

/* ********************************************************************* */
void  PrimEigenvectors (const State *state, int beg, int end)
/*!
 * Provide left and right eigenvectors and corresponding
 * eigenvalues for the primitive form of the HD equations
 * (adiabatic, pvte & isothermal cases).
 *
 * \param [in]  state  pointer to a State structure
 * \param [in]  beg    starting index of computation
 * \param [in]  end    final   index of computation
 *
 * \note It is highly recommended that state->LL and state->RR 
 *       be initialized to zero *BEFORE* since only non-zero entries 
 *       are treated here.
 *
 *  Wave names and their order are defined as enumeration constants in
 *  mod_defs.h.
 *
 *  Advection modes associated with passive scalars are simple cases
 *  for which lambda = u (entropy mode) and l = r = (0, ... , 1, 0, ...).
 *  For this reason they are NOT defined here and must be treated
 *  seperately.
 *
 * \b References:
 *    - "Riemann Solvers and Numerical Methods for Fluid Dynamics", \n
 *      Toro, 1997 Springer-Verlag, Eq. [3.18]
 *
 *********************************************************************** */
{
  int     i, j, k;
  double  rhocs, rho_cs, cs;
  double  cs2, *q, **LL, **RR, *lambda;
#if CHECK_EIGENVECTORS == YES
  double Aw1[NFLX], Aw0[NFLX], AA[NFLX][NFLX], a;
#endif

  for (i = beg; i <= end; i++){
    q      = state->v[i];
    cs2    = state->a2[i];
    LL     = state->Lp[i];
    RR     = state->Rp[i];
    lambda = state->lambda[i];

    cs = sqrt(cs2);
  
    rhocs  = q[RHO]*cs;
    rho_cs = q[RHO]/cs;
  
  /* ------------------------------------------------------
     1. Compute RIGHT eigenvectors 
     ------------------------------------------------------ */
  
    lambda[0]  = q[VXn] - cs;  /*  lambda = u - c   */
    RR[RHO][0] =  0.5*rho_cs;
    RR[VXn][0] = -0.5;
  #if HAVE_ENERGY
    RR[PRS][0] =  0.5*rhocs;
  #endif
  
    lambda[1]  = q[VXn] + cs;  /*  lambda = u + c   */  
    RR[RHO][1] = 0.5*rho_cs;
    RR[VXn][1] = 0.5; 
  #if HAVE_ENERGY
    RR[PRS][1] = 0.5*rhocs;
  #endif
 
    for (k = 2; k < NFLX; k++) lambda[k] = q[VXn];  /*  lambda = u   */
  
  #if HAVE_ENERGY
    EXPAND(RR[RHO][2] = 1.0;  ,
           RR[VXt][3] = 1.0;  ,
           RR[VXb][4] = 1.0;)
  #elif EOS == ISOTHERMAL
    EXPAND(                   ,
           RR[VXt][2] = 1.0;  ,
           RR[VXb][3] = 1.0;)
  #endif
  
  /* ------------------------------------------------------
     2. Compute LEFT eigenvectors 
     ------------------------------------------------------ */
  
    LL[0][VXn] = -1.0;  
  #if HAVE_ENERGY
    LL[0][PRS] =  1.0/rhocs;
  #elif EOS == ISOTHERMAL
    LL[0][RHO] =  1.0/rho_cs;
  #endif
  
    LL[1][VXn] = 1.0;  
  #if HAVE_ENERGY
    LL[1][PRS] = 1.0/rhocs;
  #elif EOS == ISOTHERMAL
    LL[1][RHO] = 1.0/rho_cs;
  #endif
  
  #if HAVE_ENERGY
    EXPAND(LL[2][RHO] = 1.0; ,
           LL[3][VXt] = 1.0; ,
           LL[4][VXb] = 1.0;)
  
    LL[2][PRS] = -1.0/cs2;
  #elif EOS == ISOTHERMAL
    EXPAND(                  ,
           LL[2][VXt] = 1.0; ,
           LL[3][VXb] = 1.0;)
  #endif
  
  /* ------------------------------------------------------------
     3. If required, verify eigenvectors consistency by
  
      1) checking that A = L.Lambda.R, where A is
         the Jacobian dF/dU
      2) verify orthonormality, L.R = R.L = I
     ------------------------------------------------------------ */
  
  #if CHECK_EIGENVECTORS == YES
  {
    static double **A, **ALR;
    double dA;
  
    if (A == NULL){
      A   = ARRAY_2D(NFLX, NFLX, double);
      ALR = ARRAY_2D(NFLX, NFLX, double);
      #if COMPONENTS != 3
       print ("! PrimEigenvectors(): eigenvector check requires 3 components\n");
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
  
    #if HAVE_ENERGY
     A[RHO][RHO] = q[VXn]    ; A[RHO][VXn] = q[RHO];
     A[VXn][VXn] = q[VXn]    ; A[VXn][PRS] = 1.0/q[RHO];
     A[VXt][VXt] = q[VXn]    ;  
     A[VXb][VXb] = q[VXn]    ;  
     A[PRS][VXn] = cs2*q[RHO]; A[PRS][PRS] =  q[VXn];
    #elif EOS == ISOTHERMAL
     A[RHO][RHO] = q[VXn]; A[RHO][VXn] = q[RHO];
     A[VXn][VXn] = q[VXn];
     A[VXn][RHO] = cs2/q[RHO];
     A[VXt][VXt] = q[VXn];
     A[VXb][VXb] = q[VXn];
    #endif
  
    for (i = 0; i < NFLX; i++){
    for (j = 0; j < NFLX; j++){
      ALR[i][j] = 0.0;
      for (k = 0; k < NFLX; k++) ALR[i][j] += RR[i][k]*lambda[k]*LL[k][j];
    }}
  
    for (i = 0; i < NFLX; i++){
    for (j = 0; j < NFLX; j++){
      dA = ALR[i][j] - A[i][j];
      if (fabs(dA) > 1.e-8){
        print ("! PrimEigenvectors: eigenvectors not consistent\n");
        print ("! A[%d][%d] = %16.9e, R.Lambda.L[%d][%d] = %16.9e\n",
                  i,j, A[i][j], i,j,ALR[i][j]);
        print ("! cs2 = %12.6e\n",cs2);
        print ("\n\n A = \n");   ShowMatrix(A, NFLX, 1.e-8);
        print ("\n\n R.Lambda.L = \n"); ShowMatrix(ALR, NFLX, 1.e-8);
        QUIT_PLUTO(1);
      }
    }}  
  
  /* -- check orthornomality -- */
  
    for (i = 0; i < NFLX; i++){
    for (j = 0; j < NFLX; j++){
      a = 0.0;
      for (k = 0; k < NFLX; k++) a += LL[i][k]*RR[k][j];
      if ( (i == j && fabs(a-1.0) > 1.e-8) ||
           (i != j && fabs(a)>1.e-8) ) {
        print ("! PrimEigenvectors: Eigenvectors not orthogonal\n");
        print ("!   i,j = %d, %d  %12.6e \n",i,j,a);
        print ("!   g_dir: %d\n",g_dir);
        QUIT_PLUTO(1);
      }
    }}
  }
#endif
  } /* End loop i = beg, end */
}
/* ********************************************************************* */
void ConsEigenvectors (double *u, double *v, double a2, 
                       double **LL, double **RR, double *lambda)
/*!
 * Provide conservative eigenvectors for HD equations.
 *
 * \param [in]   u   array of conservative variables
 * \param [in]   v   array of primitive variables
 * \param [in]  a2   square of sound speed
 * \param [out] LL   left conservative eigenvectors
 * \param [out] RR   right  conservative eigenvectors
 * \param [out] lambda eigenvalues
 *
 * \b References:
 *    - "Riemann Solvers and Numerical Methods for Fluid Dynamics", \n
 *      Toro, 1997 Springer-Verlag  (Page 107)
 *
 *********************************************************************** */
{
  int    i, j, k, nv;
  double H, a, gt_1, vmag2;

  #if EOS == PVTE_LAW
   print( "! ConsEigenvectors: cannot be used presently\n");
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
   #if EOS == IDEAL
    u[ENG] =   v[PRS]/(g_gamma-1.0)
            + 0.5*v[RHO]*(v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3]);
   #endif
  #endif

  #if EOS == IDEAL
   gt_1 = 1.0/(g_gamma - 1.0);
/*
   a2   = g_gamma*v[PRS]/v[RHO];
   a    = sqrt(a2);
*/
  #elif EOS == ISOTHERMAL
/*
   a2 = g_isoSoundSpeed2;
   a  = g_isoSoundSpeed;
*/
  #endif
  a = sqrt(a2);

  vmag2 = EXPAND(v[VXn]*v[VXn], + v[VXt]*v[VXt], + v[VXb]*v[VXb]);

/* -----------------------------------------------
    to compute H we use primitive variable since
    U and V may have been averaged in WENO_RF and
    are not the map of each other.
    Otherwise left and right eigenvectors wouldn't
    be orthonormal.
   ----------------------------------------------- */
/*
  H  = (u[ENG] + v[PRS])/v[RHO];
*/
  #if EOS == IDEAL
   H = 0.5*vmag2 + a2*gt_1;
  #endif

/* ======================================================
                RIGHT EIGENVECTORS  
   ====================================================== */

/*       lambda = u - c       */

  k = 0;
  lambda[k] = v[VXn] - a;

  RR[RHO][k] = 1.0;
  EXPAND(RR[MXn][k] = lambda[k];  ,
         RR[MXt][k] = v[VXt];      ,
         RR[MXb][k] = v[VXb]; )
  #if EOS == IDEAL
   RR[ENG][k] = H - a*v[VXn];
  #endif
     
/*       lambda = u + c       */
 
  k = 1;
  lambda[k] = v[VXn] + a;

  RR[RHO][k] = 1.0;
  EXPAND(RR[MXn][k] = lambda[k];  ,
         RR[MXt][k] = v[VXt];      ,
         RR[MXb][k] = v[VXb];)
  #if EOS == IDEAL
   RR[ENG][k] = H + a*v[VXn];
  #endif

/*       lambda = u        */

  #if EOS == IDEAL
   k = 2;
   lambda[k] = v[VXn];

   RR[RHO][k] = 1.0;  
   EXPAND(RR[MXn][k] = v[VXn];  ,
          RR[MXt][k] = v[VXt];  ,
          RR[MXb][k] = v[VXb];)
   RR[ENG][k] = 0.5*vmag2;  

   #if COMPONENTS > 1
    k = 3;
    lambda[k] = v[VXn];
    RR[MXt][k] = 1.0;  
    RR[ENG][k] = v[VXt];  
   #endif

   #if COMPONENTS == 3
    k = 4;
    lambda[k] = v[VXn];
    RR[MXb][k] = 1.0;  
    RR[ENG][k] = v[VXb];  
   #endif
  #elif EOS == ISOTHERMAL 
   EXPAND(                                      ,
          lambda[2] = v[VXn]; RR[MXt][2] = 1.0;   ,
          lambda[3] = v[VXn]; RR[MXb][3] = 1.0;)  
  #endif

/* ======================================================
                LEFT EIGENVECTORS  
   ====================================================== */

/*       lambda = u - c       */

  k = 0;
  #if EOS == IDEAL
   LL[k][RHO] = H + a*gt_1*(v[VXn] - a); 
   EXPAND(LL[k][MXn] = -(v[VXn] + a*gt_1); ,
          LL[k][MXt] = -v[VXt];            ,
          LL[k][MXb] = -v[VXb];)
   LL[k][ENG] = 1.0;
  #elif EOS == ISOTHERMAL
   LL[k][RHO] = 0.5*(1.0 + v[VXn]/a);
   LL[k][MXn] = -0.5/a; 
  #endif 

/*       lambda = u + c       */

  k = 1;
  #if EOS == IDEAL
   LL[k][RHO] = H - a*gt_1*(v[VXn] + a); 
   EXPAND(LL[k][MXn] = -v[VXn] + a*gt_1;  ,
          LL[k][MXt] = -v[VXt];         ,
          LL[k][MXb] = -v[VXb];)
    LL[k][ENG] = 1.0;
  #elif EOS == ISOTHERMAL
   LL[k][RHO] = 0.5*(1.0 - v[VXn]/a);
   LL[k][MXn] = 0.5/a; 
  #endif

/*       lambda = u       */

  #if EOS == IDEAL
   k = 2;
   LL[k][RHO] = -2.0*H + 4.0*gt_1*a2; 
   EXPAND(LL[k][MXn] = 2.0*v[VXn];   ,
          LL[k][MXt] = 2.0*v[VXt];   ,
          LL[k][MXb] = 2.0*v[VXb];)
   LL[k][ENG] = -2.0;

   #if COMPONENTS > 1
    k = 3;
    LL[k][RHO] = -2.0*v[VXt]*a2*gt_1; 
    LL[k][MXt] = 2.0*a2*gt_1;    
    LL[k][ENG] = 0.0;
   #endif

   #if COMPONENTS == 3
    k = 4; 
    LL[k][RHO] = -2.0*v[VXb]*a2*gt_1; 
    LL[k][MXb] = 2.0*a2*gt_1;
   #endif

   for (k = 0; k < NFLX; k++){   /* normalization */
   for (nv = 0; nv < NFLX; nv++){
     LL[k][nv] *= (g_gamma - 1.0)/(2.0*a2);
   }}

  #elif EOS == ISOTHERMAL
   EXPAND(                                              ,
          k = 2; LL[k][RHO] = -v[VXt]; LL[k][VXt] =  1.0;  ,
          k = 3; LL[k][RHO] = -v[VXb]; LL[k][VXb] =  1.0; )
  #endif

/* -----------------------------------------
         Check eigenvectors consistency
   ----------------------------------------- */

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
  #if EOS == IDEAL
   A[RHO][MXn] = 1.0;

   A[MXn][RHO] = -v[VXn]*v[VXn] + (g_gamma - 1.0)*vel2*0.5;  
   A[MXn][MXn] = (3.0 - g_gamma)*v[VXn];
   A[MXn][MXt] = (1.0 - g_gamma)*v[VXt];
   A[MXn][MXb] = (1.0 - g_gamma)*v[VXb];
   A[MXn][ENG] = g_gamma - 1.0;

   A[MXt][RHO] = - v[VXn]*v[VXt];
   A[MXt][MXn] =   v[VXt];
   A[MXt][MXt] =   v[VXn];

   A[MXb][RHO] = - v[VXn]*v[VXb];
   A[MXb][MXn] =   v[VXb];
   A[MXb][MXb] =   v[VXn];

   A[ENG][RHO] =  0.5*(g_gamma - 1.0)*vel2*v[VXn] - (u[ENG] + v[PRS])/v[RHO]*v[VXn];

   A[ENG][MXn] =  (1.0 - g_gamma)*v[VXn]*v[VXn] + (u[ENG] + v[PRS])/v[RHO];

   A[ENG][MXt] = (1.0 - g_gamma)*v[VXn]*v[VXt];
   A[ENG][MXb] = (1.0 - g_gamma)*v[VXn]*v[VXb];

   A[ENG][ENG] = g_gamma*v[VXn];
  #elif EOS == ISOTHERMAL
   A[RHO][MXn] = 1.0;

   A[MXn][RHO] = -v[VXn]*v[VXn] + a2;
   A[MXn][MXn] = 2.0*v[VXn];

   A[MXt][RHO] = - v[VXn]*v[VXt];
   A[MXt][MXn] =   v[VXt];
   A[MXt][MXt] =   v[VXn];

   A[MXb][RHO] = - v[VXn]*v[VXb];
   A[MXb][MXn] =   v[VXb];
   A[MXb][MXb] =   v[VXn];
  #endif

  for (i = 0; i < NFLX; i++){
  for (j = 0; j < NFLX; j++){
    ALR[i][j] = 0.0;
    for (k = 0; k < NFLX; k++){
      ALR[i][j] += RR[i][k]*lambda[k]*LL[k][j];
    }
  }}

  for (i = 0; i < NFLX; i++){
  for (j = 0; j < NFLX; j++){
    if (fabs(ALR[i][j] - A[i][j]) > 1.e-6){
      print ("! ConsEigenvectors: eigenvectors not consistent\n");
      print ("! g_dir = %d\n",g_dir);
      print ("! A[%d][%d] = %16.9e, R.Lambda.L[%d][%d] = %16.9e\n",
                i,j, A[i][j], i,j,ALR[i][j]);
      print ("\n\n A   = \n"); ShowMatrix(A, NFLX, 1.e-8);
      print ("\n\n R.Lambda.L = \n"); ShowMatrix(ALR, NFLX, 1.e-8);
      
      print ("\n\n RR   = \n"); ShowMatrix(RR, NFLX, 1.e-8);
      print ("\n\n LL   = \n"); ShowMatrix(LL, NFLX, 1.e-8);
      QUIT_PLUTO(1);
    }
  }}

/* -- check orthornomality -- */

  for (i = 0; i < NFLX; i++){
  for (j = 0; j < NFLX; j++){
    dA = 0.0;
    for (k = 0; k < NFLX; k++) dA += LL[i][k]*RR[k][j];
    if ( (i == j && fabs(dA-1.0) > 1.e-8) || 
         (i != j && fabs(dA)>1.e-8) ) {
      print ("! ConsEigenvectors: Eigenvectors not orthogonal\n");
      print ("!   i,j = %d, %d  %12.6e \n",i,j,dA);
      print ("!   g_dir: %d\n",g_dir);
      QUIT_PLUTO(1);
    }
  }}
}
#endif
}
/* ********************************************************************* */
void PrimToChar (double **L, double *v, double *w)
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
 * \param [in]   L   Left eigenvectors
 * \param [in]   v   (difference of) primitive variables
 * \param [out]  w   (difference of) characteristic variables
 *
 *********************************************************************** */
{
  int nv;

  #if HAVE_ENERGY
   w[0] = L[0][VXn]*v[VXn] + L[0][PRS]*v[PRS];
   w[1] = L[1][VXn]*v[VXn] + L[1][PRS]*v[PRS];
   EXPAND( w[2] = v[RHO] + L[2][PRS]*v[PRS];  ,
           w[3] = v[VXt];                   ,
           w[4] = v[VXb];)
  #elif EOS == ISOTHERMAL
    w[0] = L[0][RHO]*v[RHO] + L[0][VXn]*v[VXn];
    w[1] = L[1][RHO]*v[RHO] + L[1][VXn]*v[VXn];
    EXPAND(                     ,
            w[2] = v[VXt];       ,
            w[3] = v[VXb];)
  #endif

/* ----------------------------------------------- 
     For passive scalars, the characteristic 
     variable is equal to the primitive one, 
     since  l = r = (0,..., 1 , 0 ,....)
   ----------------------------------------------- */		   

#if NSCL > 0
  NSCL_LOOP(nv)  w[nv] = v[nv];
#endif   
}
