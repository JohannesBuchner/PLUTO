/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Convert between primitive and conservative variables.

  The PrimToCons() converts an array of primitive quantities to 
  an array of conservative variables for the RHD equations.
  
  The ConsToPrim() converts an array of conservative quantities to 
  an array of primitive quantities.
  During the conversion, pressure is normally recovered from total 
  energy unless zone has been tagged with FLAG_ENTROPY.
  In this case we recover pressure from conserved entropy:
  
      if (FLAG_ENTROPY is TRUE)  --> p = p(S)
      else                       --> p = p(E)
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   June 25, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void PrimToCons  (double *uprim[], double *ucons[],
                 int ibeg, int iend)
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
  int     nv, ii;
  double  rhoh_g2, scrh, g;
  double  beta_fix = 0.9999;
  double  *u, *v;
  static double  *h;

  if (h == NULL){
    h = ARRAY_1D(NMAX_POINT, double);
  }

  Enthalpy (uprim, h, ibeg, iend);

  for (ii = ibeg; ii <= iend; ii++) {
   
    u = ucons[ii];
    v = uprim[ii];

    g = EXPAND(v[VX1]*v[VX1], +v[VX2]*v[VX2], +v[VX3]*v[VX3]);

    if (g >= 1.0){
      WARNING( 
        print ("! u^2 > 1 (%f) in PrimToCons\n", scrh);
        Where (ii, NULL);
      )

      g = beta_fix/sqrt(g);
      EXPAND(v[VX1] *= g;  ,
             v[VX2] *= g;  ,
             v[VX3] *= g;)

      g = beta_fix*beta_fix;
    }
    g    = 1.0/(1.0 - g);
    scrh = rhoh_g2 = v[RHO]*h[ii]*g;
    g    = sqrt(g);

    u[RHO] = v[RHO]*g;
    EXPAND(u[MX1] = scrh*v[VX1];  ,
           u[MX2] = scrh*v[VX2];  ,
           u[MX3] = scrh*v[VX3];)

    ucons[ii][ENG] = rhoh_g2 - v[PRS];

#if NSCL > 0    
    NSCL_LOOP(nv) u[nv] = v[nv]*u[RHO];
#endif    

    #if CHECK_CONSERVATIVE_VAR == YES
     m2 = EXPAND(u[MX1]*u[MX1], + u[MX2]*u[MX2], + u[MX3]*u[MX3]);
     g  = u[ENG] - sqrt(m2);
     if (g <= 0.0){
       printf ("! E - m < 0 in PrimToCons (%12.6e)\n",g);
       QUIT_PLUTO(1);
     }
     g = u[RHO]/u[ENG] - sqrt(1.0 - m2/u[ENG]/u[ENG]);
     if (g >=0){
        print ("! g > 0.0 in PrimToCons(2) (%12.6e)\n",g);
        Show(ucons,ii);
        QUIT_PLUTO(1);
     }
    #endif
  }
}

/* ********************************************************************* */
int ConsToPrim (double **ucons, double **uprim, int ibeg, int iend, 
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
  int    nv, i, err, ifail;
  int    use_entropy, use_energy=1;
  double  scrh, m, g;
  double *u, *v;
  Map_param par;

  ifail = 0;
  for (i = ibeg; i <= iend; i++) {

    u = ucons[i];
    v = uprim[i];

  /* -----------------------------------------------------------
      Define the input parameters of the parameter structure
     ----------------------------------------------------------- */

    par.D  = u[RHO];
    par.E  = u[ENG];
    par.m2 = EXPAND(u[MX1]*u[MX1], + u[MX2]*u[MX2], + u[MX3]*u[MX3]);

  /* -------------------------------------------
        Check density and energy positivity 
     ------------------------------------------- */
  
    if (u[RHO] < 0.0) {
      print("! ConsToPrim: negative density (%8.2e), ", u[RHO]);
      Where (i, NULL);
      u[RHO]   = g_smallDensity;
      flag[i] |= FLAG_CONS2PRIM_FAIL;
      ifail    = 1;
    }

    if (u[ENG] < 0.0) {
      WARNING(
        print("! ConsToPrim: negative energy (%8.2e), ", u[ENG]);
        Where (i, NULL);
      )
      u[ENG]   = sqrt(1.e-8 + par.m2 + u[RHO]*u[RHO]);
      flag[i] |= FLAG_CONS2PRIM_FAIL;
      ifail    = 1;
    }

  /* --------------------------------
      recover pressure by inverting 
      the energy or entropy equation.
     -------------------------------- */

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

  /*  ------------------------------------------
              complete conversion 
      ------------------------------------------ */

    v[PRS] = par.prs;
    scrh   = 1.0/(u[ENG] + v[PRS]);  /* = 1 / W */

    EXPAND(v[VX1] = u[MX1]*scrh; ,
           v[VX2] = u[MX2]*scrh; ,
           v[VX3] = u[MX3]*scrh;)

    g = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
    g = 1.0/sqrt(1.0 - g);
    v[RHO] = u[RHO]/g;

#if NSCL > 0     
    NSCL_LOOP(nv)  v[nv] = u[nv]/u[RHO];
#endif

  }

  return ifail;
}
