/* *************************************************************** */
/*!
  \file
  \brief Runge-Kutta-Legendre driver for integration of diffusion terms.
  
  Take one step in the solution of the diffusion equation \f$ dQ/dT = F\f$
  using Runge-Kutta-Legendre (RKL) method.
  Equations are evolved in time by taking \c s sub-steps including only the
  diffusion (parabolic) terms.   
  The full step \c tau equals the current integration time step (from hydro
  or a fixed one).

  This driver can be used with a 1st or 2nd order algorithm (default is 2).

  For the 1st order method (\c RKL_ORDER = 1):
  \f[ \begin{array}{lcl}
      Y_0 & = & U^n                            \\ \noalign{\medskip}
      Y_1 & = & Y_0 + \tilde{\mu}_1\tau M(Y_0) \\ \noalign{\medskip}
      Y_j & = & \mu_j Y_{j-1} + \nu_j Y_{j-2} + \tilde{\mu}_j \tau M(Y_{j-1}) 
            \quad \rm{for} \quad j = 2,...,s   \\ \noalign{\medskip}
      U^{n+1} & = & Y_s
      \end{array}
  \f]
  where
  \f[ \begin{array}{lclclcl}
      b_0 = b_1 = b_2 & = & 1 & &
      b_j & = & 1 \\ \noalign{\medskip}
      a_j & = & 0.0 & & w_1 & = & \frac{2}{s^2 + s} \\ \noalign{\medskip}
      \mu_j & = & \frac{2j - 1}{j} w_1 & &
      \nu_j & = & -\frac{j-1}{j} \\ \noalign{\medskip}
      \tilde{\mu}_j & = & w_1 \mu_j & & 
              \quad \rm{for} \quad j = 2,...,s   
      \end{array}
  \f]


  For the 2nd order method (\c RKL_ORDER = 2):
  \f[ \begin{array}{lcl}
      Y_0 & = & U^n                            \\ \noalign{\medskip}
      Y_1 & = & Y_0 + \tilde{\mu}_1\tau M(Y_0) \\ \noalign{\medskip}
      Y_j & = & \mu_j Y_{j-1} + \nu_j Y_{j-2} + (1 - \mu_j - \nu_j) Y_0
            + \tilde{\mu}_j \tau M(Y_{j-1}) + \tilde{\gamma}_j \tau M(Y_0)
            \quad \rm{for} \quad j = 2,...,s   \\ \noalign{\medskip}
      U^{n+1} & = & Y_s
      \end{array}
  \f]
  where
  \f[ \begin{array}{lclclcl}
      b_0 = b_1 = b_2 & = & \frac{1}{3} & &
              b_j & = & \frac{j^2 + j - 2}{2j(j + 1)} \\ \noalign{\medskip}
      a_j & = & 1 - b_j & & w_1 & = & \frac{4}{s^2 + s - 2} \\ \noalign{\medskip}
      \mu_j & = & \frac{2j - 1}{j} \frac{b_j}{b_{j-1}} & &
              \nu_j & = & -\frac{j-1}{j} \frac{b_j}{b_{j-2}} \\ \noalign{\medskip}
      \tilde{\mu}_j & = & w_1 \mu_j & & \tilde{\gamma}_j = -a_{j-1}\tilde{\mu}_j
              \quad \rm{for} \quad j = 2,...,s   
      \end{array}
  \f]
 
  
  \todo
    - selective update on variables
    - remove tau from multiplication inside inner loop (put it outside)

  \b References
     - Meyer, C., Balsara, D., \& Aslam T., 2012, 
       Mon. Not. R. Astron. Soc., 422

  \authors L. Rickler   (luca.rickler@edu.unito.it)\n
           A. Mignone (mignone@ph.unito.it)
    
  \date    March 6, 2018
*/
/* /////////////////////////////////////////////////////////////// */
#include "pluto.h"

#ifndef RKL_ORDER
 #define RKL_ORDER 2
#endif

/* *************************************************************** */
void RKL (const Data *d, double dt, timeStep *Dts, Grid *grid)
/*!
 *Solve diffusion equation using Runge-Kutta-Legendre (RKL) method
 *
 * \param [in,out] d     pointer to Data structure
 * \param [in]     dt    the time step increment
 * \param [in,out] Dts   pointer to timeStep structure
 * \param [in]     grid  pointer to an array of Grid structures
 *
 ***************************************************************** */
{
  int i, j, k, nv, s, s_RKL = 0;
  int nv_indx, var_list[NVAR], nvar_rkl;
  double mu_j, nu_j, mu_tilde_j, gamma_j, Y;
  double a_jm1, b_j, b_jm1, b_jm2, w1;
  double tau = dt, t0 = g_time;
  double dt_par, scrh;
  static Data_Arr Y_jm1, Y_jm2, MY_jm1, MY_0;
  static double **v;
  double s_str;
  RBox box;

/* --------------------------------------------------------
   0. Initialize arrays
   -------------------------------------------------------- */

  RBoxDefine (IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &box);

  if (Y_jm1 == NULL) { 
    Y_jm1  = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);  
    Y_jm2  = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    MY_0   = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    MY_jm1 = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    v      = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

/* --------------------------------------------------------
   1. Set the number of variables to be evolved with RKL
   -------------------------------------------------------- */
  
  i = 0;
  for (nv = 0; nv < NVAR; nv++) var_list[nv] = 0;
  #if VISCOSITY == RK_LEGENDRE
  EXPAND(var_list[i++] = MX1;  ,
         var_list[i++] = MX2;  ,
         var_list[i++] = MX3;)
  #endif

  #if RESISTIVITY == RK_LEGENDRE
  EXPAND(var_list[i++] = BX1;  ,
         var_list[i++] = BX2;  ,
         var_list[i++] = BX3;)
  #endif
  #if    (THERMAL_CONDUCTION == RK_LEGENDRE)  \
      || (VISCOSITY     == RK_LEGENDRE && EOS == IDEAL) \
      || (RESISTIVITY   == RK_LEGENDRE && EOS == IDEAL)
   
  var_list[i++] = ENG;
  #endif
  nvar_rkl = i;

/* --------------------------------------------------------
   2. Obtain conservative vector Uc
   -------------------------------------------------------- */

  PrimToCons3D(d->Vc, d->Uc, &box);

/* --------------------------------------------------------
   3. Compute the parabolic time step by calling 
      the RHS function once outside the main loop.
   -------------------------------------------------------- */
  
  g_intStage = 1;
  Boundary(d, ALL_DIR, grid);

  Dts->invDt_par  = ParabolicRHS(d, MY_0, &box, NULL, RK_LEGENDRE, 1.0, grid);
  Dts->invDt_par /= (double) DIMENSIONS;  
#ifdef PARALLEL
  MPI_Allreduce (&Dts->invDt_par, &scrh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  Dts->invDt_par = scrh;
#endif

  dt_par = Dts->cfl_par/(2.0*Dts->invDt_par); /* -- explicit parabolic time step -- */

/* --------------------------------------------------------
   4. Compute number of RKL steps from Dts->inv_dta
   -------------------------------------------------------- */   

#if RKL_ORDER == 1
  scrh  = tau/dt_par;                            /* Solution of quadratic Eq.  */
  s_str = 4.0*scrh/(1.0 + sqrt(1.0 + 8.0*scrh)); /* 2*tau/dt_exp = s^2 + s     */
#elif RKL_ORDER == 2
  scrh  = tau/dt_par;                      /*  Solution of quadratic Eq.   */
  s_str =   4.0*(1.0 + 2.0*scrh)           /*  4*tau/dt_exp = s^2 + s - 2  */
          /(1.0 + sqrt(9.0 + 16.0*scrh)); 
#else
  #error Invalid RKL_ORDER
#endif
  s_RKL = 1 + INT_FLOOR(s_str);
  Dts->Nrkl = s_RKL;

/* --------------------------------------------------------
   5. Unflag zone tagged with the ENTROPY_SWITCH since
      only total energy can be evolved using STS
   -------------------------------------------------------- */

  #if ENTROPY_SWITCH
  TOT_LOOP(k,j,i) d->flag[k][j][i] &= ~FLAG_ENTROPY;
  #endif

/* --------------------------------------------------------
   6. Compute coefficients
   -------------------------------------------------------- */

#if RKL_ORDER == 1
  w1 = 2.0/(s_RKL*s_RKL + s_RKL);
  mu_tilde_j = w1;
#elif RKL_ORDER == 2
  w1 = 4.0/(s_RKL*s_RKL + s_RKL - 2.0); 
  mu_tilde_j = w1/3.0;

  b_j = b_jm1 = b_jm2 = 1.0/3.0; 
  a_jm1 = 1.0 - b_jm1;
#endif

  DOM_LOOP (k,j,i) {
    NVAR_LOOP(nv) Y_jm1[k][j][i][nv] = Y_jm2[k][j][i][nv] = d->Uc[k][j][i][nv];
    for (nv_indx = 0; nv_indx < nvar_rkl; nv_indx++) {
      nv = var_list[nv_indx];
      Y_jm1[k][j][i][nv] = Y_jm2[k][j][i][nv] + mu_tilde_j*tau*MY_0[k][j][i][nv];
    } 
  }
 
  ConsToPrim3D(Y_jm1, d->Vc, d->flag, &box);
  /* s loop */
  s = 1;
#if RKL_ORDER == 1  
  g_time = t0 + 0.5*tau*(s*s+s)*w1;
#elif RKL_ORDER == 2
  g_time = t0 + 0.25*tau*(s*s+s-2)*w1;
#endif

/* --------------------------------------------------------
   7. RKL main stage loop 
   -------------------------------------------------------- */

  for (s = 2; s <= s_RKL; s++) {
  
    #if RKL_ORDER == 1
    mu_j       = (2.*s -1.)/s; /* Eq. [17] */
    mu_tilde_j = w1*mu_j;                                     
    nu_j       = -(s -1.)/s;
    #elif RKL_ORDER == 2
    mu_j       = (2.*s -1.)/s * b_j/b_jm1;   /* Eq. [17] */
    mu_tilde_j = w1*mu_j;
    gamma_j    = -a_jm1*mu_tilde_j;
    nu_j       = -(s -1.)*b_j/(s*b_jm2);

    b_jm2 = b_jm1;    /* Eq. [16] */
    b_jm1 = b_j;
    a_jm1 = 1.0 - b_jm1;
    b_j   = 0.5*(s*s+3.0*s)/(s*s+3.0*s+2);
    #endif
    
    g_intStage = s;
    Boundary(d, ALL_DIR, grid);
    ParabolicRHS (d, MY_jm1, &box, NULL, RK_LEGENDRE, 1.0, grid);
    DOM_LOOP (k,j,i){ 
      for (nv_indx = 0; nv_indx < nvar_rkl; nv_indx++) {  
        nv = var_list[nv_indx];

     /* Eq. [15] */

        Y                  = mu_j*Y_jm1[k][j][i][nv] + nu_j*Y_jm2[k][j][i][nv];
        Y_jm2[k][j][i][nv] = Y_jm1[k][j][i][nv];
        #if RKL_ORDER == 1 
        Y_jm1[k][j][i][nv] = Y + tau*mu_tilde_j*MY_jm1[k][j][i][nv];
        #elif RKL_ORDER == 2
        Y_jm1[k][j][i][nv] = Y + (1.0 - mu_j - nu_j)*d->Uc[k][j][i][nv]
                               + tau*mu_tilde_j*MY_jm1[k][j][i][nv] 
                               + gamma_j*tau*MY_0[k][j][i][nv];
        #endif                  
      }                                              
    } /* END DOM_LOOP  */

  /* -- Update staggered magnetic field -- */

    #if (defined STAGGERED_MHD) && (RESISTIVITY == RK_LEGENDRE)
    #error RKL and CONSTRAINED_TRANSPORT not compatible yet.
    #endif

    #if RKL_ORDER == 1
    g_time = t0 + 0.5*tau*(s*s+s)*w1;
    #elif RKL_ORDER == 2
    g_time = t0 + 0.25*tau*(s*s+s-2)*w1;
    #endif
 
  /* ------------------------------------------------------
     7b. Convert conservative variables to primitive for
         the next stage  (or timestep if s == s_RKL)
     ---------------------------------------------------------- */
        
    ConsToPrim3D(Y_jm1, d->Vc, d->flag, &box);
  }/* s loop */

  DOM_LOOP (k,j,i) NVAR_LOOP(nv) d->Uc[k][j][i][nv] = Y_jm1[k][j][i][nv];
   
  g_time = t0;
}
