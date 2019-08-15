/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Convert between primitive and conservative variables.

  The PrimToCons() converts an array of primitive quantities to 
  an array of conservative variables for the RMHD equations.
  
  The ConsToPrim() converts an array of conservative quantities to 
  an array of primitive quantities.
  During the conversion, pressure is normally recovered from total 
  energy using the algorithm outlined in
  - "Equation of sweep in relativistic magnetohydrodynamics: variable versus
     constant adiabatic index"\n
     Mignone \& Mc Kinney, MNRAS (2007) 378, 1118.

  However, if the zone has been tagged with FLAG_ENTROPY, primitive
  variables are recovered by using the conserved entropy rather than
  total energy.
  
  In other words:
  
      if (FLAG_ENTROPY is TRUE)  --> p = p(S)
      else                       --> p = p(E)
  
  If the inversion scheme fails and p cannot be obtained the
  PressureFix() function is used.

  \author A. Mignone (mignone@ph.unito.it)
  \date   June 25, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"


#if RMHD_RESISTIVE_CHECK
static int CheckConversionQuartic (double rho, double *m, double *B, double *E, 
                             double Enr, double glor, double *F_gamma, double *dF_gamma);
#endif

/* ********************************************************************* */
void PrimToCons (double **uprim, double **ucons, int beg, int end)
/*!
 * Convert primitive variables to conservative variables. 
 *
 * \param [in]  uprim array of primitive variables
 * \param [out] ucons array of conservative variables
 * \param [in]  beg   starting index of computation
 * \param [in]  end   final index of computation
 *
 *********************************************************************** */
{
  int   i, nv;
  double  vel2, usq, vB, Bmag2;
  double  vx1, vx2, vx3;
  double  g, g2, wt;
  double  *u, *v;
  static double *h;
  #if EOS == IDEAL
   double gmmr = g_gamma/(g_gamma - 1.0);
  #endif

  if (h == NULL) h = ARRAY_1D(NMAX_POINT, double);

  Enthalpy(uprim, h, beg, end);

  for (i = beg; i <= end; i++) {

    v = uprim[i];
    u = ucons[i];

    EXPAND(vx1 = v[VX1];  ,
           vx2 = v[VX2];  ,
           vx3 = v[VX3];)   
    vel2  = EXPAND(vx1*vx1, + vx2*vx2, + vx3*vx3);
    g2 = 1.0/(1.0 - vel2);
    g  = sqrt(g2);

    vB    = EXPAND(vx1*v[BX1], + vx2*v[BX2], + vx3*v[BX3]);
    Bmag2 = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
    wt  = v[RHO]*h[i]*g2 + Bmag2;

  /* -------------------------------------------------------
       Convert from primitive (v) to conservative (u)   
     ------------------------------------------------------- */

    u[RHO] = g*v[RHO];
    EXPAND (u[MX1] = wt*vx1 - vB*v[BX1];  ,
            u[MX2] = wt*vx2 - vB*v[BX2];  ,
            u[MX3] = wt*vx3 - vB*v[BX3];)

    EXPAND (u[BX1] = v[BX1];  ,
            u[BX2] = v[BX2];  ,
            u[BX3] = v[BX3];)

    #if RMHD_REDUCED_ENERGY == YES
     #if EOS == IDEAL
      u[ENG]  = v[PRS]*(g2*gmmr - 1.0) + u[RHO]*g2*vel2/(g + 1.0)
               + 0.5*(Bmag2*(1.0 + vel2) - vB*vB);
     #elif EOS == TAUB
      wt    = v[PRS]/v[RHO];
      u[ENG] = v[PRS]*(g2*2.5 - 1.0) 
               + u[RHO]*g2*(2.25*wt*wt + vel2)/(g*sqrt(1.0 + 2.25*wt*wt) + 1.0)
               + 0.5*(Bmag2*(1.0 + vel2) - vB*vB);
     #endif
    #else
     u[ENG]  = v[RHO]*h[i]*g2 - v[PRS] + 0.5*(Bmag2*(1.0 + vel2) - vB*vB);
    #endif

#if NSCL > 0
    NSCL_LOOP(nv) u[nv] = u[RHO]*v[nv];
#endif    
    
    #ifdef GLM_MHD
     u[PSI_GLM] = v[PSI_GLM]; 
    #endif
  }
}

/* ********************************************************************* */
int ConsToPrim (double **ucons, double **uprim, int beg, int end, 
                unsigned char *flag)
/*!
 * Convert from conservative to primitive variables.
 *
 * \param [in]  ucons      array of conservative variables
 * \param [out] uprim      array of primitive variables
 * \param [in]  beg        starting index of computation
 * \param [in]  end        final index of computation
 * \param [in,out] flag    array of flags tagging, in input, zones
 *                         where entropy must be used to recover pressure
 *                         and, on output, zones where conversion was
 *                         not successful.
 * 
 * \return Return 0 if conversion was successful in all zones in the 
 *         range [ibeg,iend].
 *         Return 1 if one or more zones could not be converted correctly
 *         and either pressure, density or energy took non-physical values. 
 *
 *********************************************************************** */
{
  int    i, nv, err, ifail;
  int    use_entropy, use_energy=1;
  double *u, *v, scrh, w_1;
  Map_param par;

  ifail = 0;
  for (i = beg; i <= end; i++) {

    u = ucons[i];
    v = uprim[i];

  /* -----------------------------------------------------------
      Define the input parameters of the parameter structure
     ----------------------------------------------------------- */

    par.D  = u[RHO];
    par.E  = u[ENG];
    par.S  = EXPAND(u[MX1]*u[BX1], + u[MX2]*u[BX2], + u[MX3]*u[BX3]);
    par.m2 = EXPAND(u[MX1]*u[MX1], + u[MX2]*u[MX2], + u[MX3]*u[MX3]);
    par.B2 = EXPAND(u[BX1]*u[BX1], + u[BX2]*u[BX2], + u[BX3]*u[BX3]); 
    par.S2 = par.S*par.S;

  /* -------------------------------------------
        Check density and energy positivity 
     ------------------------------------------- */
  
    if (u[RHO] < 0.0) {
      print("! ConsToPrim(): negative density (%8.2e), ", u[RHO]);
      Where (i, NULL);
      u[RHO]   = g_smallDensity;
      flag[i] |= FLAG_CONS2PRIM_FAIL;
      ifail    = 1;
    }

    if (u[ENG] < 0.0) {
      WARNING(
        print("! ConsToPrim(): negative energy (%8.2e), ", u[ENG]);
        Where (i, NULL);
      )
      u[ENG]   = 1.e-5;
      flag[i] |= FLAG_CONS2PRIM_FAIL;
      ifail    = 1;
    }

  /* ---------------------------------------------------------
      Attempt to recover pressure and velocity from energy or
      entropy.
      If an error occurs, use the PressureFix() function
     ------------------------------------------------------- */

    #if ENTROPY_SWITCH
    use_entropy = (flag[i] & FLAG_ENTROPY);
    use_energy  = !use_entropy;
    par.sigma_c = u[ENTR];
    if (use_entropy) {
      err = EntropySolve(&par);      
      if (err) {
        WARNING(Where (i, NULL);)
        err = PressureFix(&par);
        if (err){
          Where(i,NULL);
          QUIT_PLUTO(1);
        }
        flag[i] |= FLAG_CONS2PRIM_FAIL;
        ifail    = 1;
      }
      u[ENG] = par.E;  /* Redefine energy */
    } 
    #endif
 
    if (use_energy){
      err = EnergySolve(&par);
      if (err){
        WARNING(Where(i,NULL);)
        err = PressureFix(&par);
        if (err){
          Where(i,NULL);
          QUIT_PLUTO(1);
        }
        u[ENG]   = par.E;
        flag[i] |= FLAG_CONS2PRIM_FAIL;
        ifail    = 1;
      }
      #if ENTROPY_SWITCH     
      u[ENTR] = par.sigma_c;  /* Redefine entropy */
      #endif      
    }

 /* -- W, p and lor have been found. Now complete conversion --  */

    v[RHO] = u[RHO]/par.lor;
    v[PRS] = par.prs;

    w_1  = 1.0/(par.W + par.B2);
    scrh = par.S/par.W;
    EXPAND(v[VX1] = w_1*(u[MX1] + scrh*u[BX1]);  ,
           v[VX2] = w_1*(u[MX2] + scrh*u[BX2]);  ,
           v[VX3] = w_1*(u[MX3] + scrh*u[BX3]);)

    scrh = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
    if (scrh >= 1.0){
      print ("! ConsToPrim(): v^2 = %f > 1  (p = %12.6e); ", scrh, par.prs);
      Where (i, NULL);
      print ("!               Flag_Entropy = %d\n", (flag[i] & FLAG_ENTROPY)); 
      QUIT_PLUTO(1);
    }

    EXPAND(v[BX1] = u[BX1];  ,
           v[BX2] = u[BX2];  ,
           v[BX3] = u[BX3];)

    #if NSCL > 0 
    NSCL_LOOP(nv) v[nv] = u[nv]/u[RHO];
    #endif

    #ifdef GLM_MHD
    v[PSI_GLM] = u[PSI_GLM]; 
    #endif


    #if RMHD_RESISTIVE_CHECK
  /* -------------------------------------------------------------
        Now check the quartic that will be used in RRMHD 
     ------------------------------------------------------------- */

   /* -- Compute ideal electric field -- */

    double emf[3];
    emf[IDIR] = v[VX3]*v[BX2] - v[VX2]*v[BX3];
    emf[JDIR] = v[VX1]*v[BX3] - v[VX3]*v[BX1];
    emf[KDIR] = v[VX2]*v[BX1] - v[VX1]*v[BX2]; 

    double g, F, G, dF, g_lornew, g_lorold, p_new, v_new[3], rho_new, omega_new, E2, B2, g1 = g_gamma/(g_gamma-1.0);
    int fun, check1=0, check2=0;
    double glor_1;

    for (g = 0.95; g <= 30.0; g += 5){
      fun = CheckConversionQuartic (u[RHO], u + MX1, u + BX1, emf, u[ENG] + u[RHO], g, &F, &dF);
      fun = CheckConversionQuartic (u[RHO], u + MX1, u + BX1, emf, u[ENG] + u[RHO], g + 5, &G, &dF);
      if(F*G<=0){
        check1=1;
 
        g_lornew=g+5;
        g_lorold=g;
        while(fabs(g_lorold-g_lornew)>1.e-10){
          check2++;
          g_lorold=g_lornew;
          fun = CheckConversionQuartic (u[RHO], u + MX1, u + BX1, emf, u[ENG] + u[RHO], g_lorold, &F, &G);
          g_lornew=g_lornew-F/G;
     
          if(check2>20){
            print ("Error, Newton-Rhapson does not converge\n");
            QUIT_PLUTO(1);
          }
        }
        if(fabs(g_lornew-par.lor)>1.e-6){
          print ("Error, zero found not accurate\n");
          QUIT_PLUTO(1);
        }
      }
    }
    if(check1==0){
      print ("Error, zero not found!\n");
      QUIT_PLUTO(1);
    }

    /*Now we compute and check all the other primitive variables*/

    E2 = emf[IDIR]*emf[IDIR] + emf[JDIR]*emf[JDIR] + emf[KDIR]*emf[KDIR];
    B2 = u[BX1]*u[BX1] + u[BX2]*u[BX2] + u[BX3]*u[BX3];    

    rho_new = u[RHO]/g_lornew;
    if(fabs(v[RHO] - rho_new)>1.e-6) {
      print ("Density is not accurate!\n");
      QUIT_PLUTO(1);
    }

    p_new = (u[ENG] +u[RHO] - u[RHO]*g_lornew - 0.5*(E2 + B2))/(g_lornew*g_lornew*g1 - 1);
    if(fabs(v[PRS] - p_new)/v[PRS]/v[PRS]>1.e-6) {
      print ("Pressure is not accurate!\n");
      QUIT_PLUTO(1);
    }

    omega_new = rho_new+g1*p_new;

    glor_1 = 1/(omega_new*g_lornew*g_lornew);

    v_new[0] = (u[VX1] - emf[JDIR]*u[BX3] + emf[KDIR]*u[BX2])*glor_1;
    if(fabs(v[VX1] - v_new[0])>1.e-6) {
      print ("Velocity (X) is not accurate!\n");
      QUIT_PLUTO(1);
    }

    v_new[1] = (u[VX2] - emf[KDIR]*u[BX1] + emf[IDIR]*u[BX3])*glor_1;
    if(fabs(v[VX2] - v_new[1])>1.e-6) {
      print ("Velocity (Y) is not accurate!\n");
      QUIT_PLUTO(1);
    }

    v_new[2] = (u[VX3] - emf[IDIR]*u[BX2] + emf[JDIR]*u[BX1])*glor_1;
    if(fabs(v[VX3] - v_new[2])>1.e-6) {
      print ("Velocity (Z) is not accurate!\n");
      QUIT_PLUTO(1);
    }  
     
    #endif
  }
  return ifail;
}

#if RMHD_RESISTIVE_CHECK
/* ********************************************************************* */
int CheckConversionQuartic (double D, double *m, double *B, double *E, 
                             double tau, double glor, double *F_gamma, double *dF_gamma)
/*
 *
 *********************************************************************** */
{
  double EvecB[3], EB_2, m2, m_EB, E2, B2, C_1, C_2, g1 = ((g_gamma-1.0)/g_gamma);
  double A4, A3, A2, A1, A0;
  
  //input

  EvecB[IDIR] = E[JDIR]*B[KDIR] - E[KDIR]*B[JDIR];
  EvecB[JDIR] = E[KDIR]*B[IDIR] - E[IDIR]*B[KDIR];
  EvecB[KDIR] = E[IDIR]*B[JDIR] - E[JDIR]*B[IDIR];

  EB_2 = EvecB[IDIR]*EvecB[IDIR] + EvecB[JDIR]*EvecB[JDIR] + EvecB[KDIR]*EvecB[KDIR];
  m2 = m[IDIR]*m[IDIR] + m[JDIR]*m[JDIR] + m[KDIR]*m[KDIR];
  m_EB = -2*m[IDIR]*EvecB[IDIR] - 2*m[JDIR]*EvecB[JDIR] - 2*m[KDIR]*EvecB[KDIR];
  E2 = E[IDIR]*E[IDIR] + E[JDIR]*E[JDIR] + E[KDIR]*E[KDIR];
  B2 = B[IDIR]*B[IDIR] + B[JDIR]*B[JDIR] + B[KDIR]*B[KDIR];

  C_1 = m2 + EB_2 + m_EB;
  C_2 = tau - 0.5*E2 - 0.5*B2;

  A4 = C_1 - C_2*C_2;
  A3 = 2*C_2*g1*D;
  A2 = C_2*C_2 - 2*C_1*g1 - g1*g1*D*D;
  A1 = -2*C_2*g1*D;
  A0 = g1*g1*(C_1 + D*D);

  *F_gamma = A0 + glor*(A1 + glor*(A2 + glor*(A3 + glor*A4)));

  *dF_gamma = A1+glor*(2.0*A2+glor*(3.0*A3+4.0*glor*A4));

  return 0;
 
}
    
#endif

