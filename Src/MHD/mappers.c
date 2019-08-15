/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Convert between primitive and conservative variables.

  The PrimToCons() converts an array of primitive quantities to 
  an array of conservative variables for the MHD equations.
  
  The ConsToPrim() converts an array of conservative quantities to 
  an array of primitive quantities.
  During the conversion, pressure is normally recovered from total 
  energy unless zone has been tagged with FLAG_ENTROPY.
  In this case we recover pressure from conserved entropy:
  
      if (FLAG_ENTROPY is TRUE)  --> p = p(S)
      else                       --> p = p(E)
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   June 24, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void PrimToCons (double **uprim, double **ucons, int ibeg, int iend)
/*!
 * Convert primitive variables in conservative variables. 
 *
 * \param [in]  uprim array of primitive variables
 * \param [out] ucons array of conservative variables
 * \param [in]  beg   starting index of computation
 * \param [in]  end   final index of computation
 *
 *********************************************************************** */
{
  int  i, nv, status;
  double *v, *u;
  double rhoe, kinb2, T, gmm1;

#if EOS == IDEAL
  gmm1 = g_gamma - 1.0;
#endif
  for (i = ibeg; i <= iend; i++) {
  
    v = uprim[i];
    u = ucons[i];

    u[RHO] = v[RHO];
        
    EXPAND (u[MX1] = v[RHO]*v[VX1];  ,
            u[MX2] = v[RHO]*v[VX2];  ,
            u[MX3] = v[RHO]*v[VX3];)

    EXPAND (u[BX1] = v[BX1];  ,
            u[BX2] = v[BX2];  ,
            u[BX3] = v[BX3];)

    kinb2   = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
    kinb2   = v[RHO]*kinb2 + EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
    kinb2  *= 0.5;

    #if EOS == IDEAL
    u[ENG] = kinb2 + v[PRS]/gmm1;
    #elif EOS == PVTE_LAW
    status = GetPV_Temperature(v, &T);
    if (status != 0){
      T      = T_CUT_RHOE;
      v[PRS] = Pressure(v, T);
    }
    rhoe   = InternalEnergy(v, T);
    u[ENG] = rhoe + kinb2;

    if (u[ENG] != u[ENG]){
      print("! PrimToCons: KE:%12.6e uRHO : %12.6e, m2 : %12.6e \n",rhoe,v[RHO],u[ENG]);
      QUIT_PLUTO(1);
    }
    #endif
    
    #ifdef GLM_MHD
    u[PSI_GLM] = v[PSI_GLM]; 
    #endif
    #if NSCL > 0 
    NSCL_LOOP(nv) u[nv] = v[RHO]*v[nv];
    #endif    
  }
}
/* ********************************************************************* */
int ConsToPrim (double **ucons, double **uprim, int ibeg, int iend, 
                unsigned char *flag)
/*!
 * Convert from conservative to primitive variables.
 *
 * \param [in]     ucons   array of conservative variables
 * \param [out]    uprim   array of primitive variables
 * \param [in]     beg     starting index of computation
 * \param [in]     end     final index of computation
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
  int  i, nv, err, ifail;
  int  use_entropy, use_energy=1;
  double tau, rho, gmm1, rhoe, T;
  double b2, m2, kinb2, rhog1;
  double *u, *v;

#if EOS == IDEAL
   gmm1 = g_gamma - 1.0;
#endif

  ifail = 0;
  for (i = ibeg; i <= iend; i++) {

    u = ucons[i];
    v = uprim[i];

    m2 = EXPAND(u[MX1]*u[MX1], + u[MX2]*u[MX2], + u[MX3]*u[MX3]);
    b2 = EXPAND(u[BX1]*u[BX1], + u[BX2]*u[BX2], + u[BX3]*u[BX3]);

  /* -- Check density positivity -- */
  
    if (u[RHO] < 0.0) {
      print("! ConsToPrim: negative density (%8.2e), ", u[RHO]);
      Where (i, NULL);
      u[RHO]   = g_smallDensity;
      flag[i] |= FLAG_CONS2PRIM_FAIL;
      ifail    = 1;
    }

  /* -- Compute density, velocity, mag. field and scalars -- */

    v[RHO] = rho = u[RHO];
    tau = 1.0/u[RHO];
    EXPAND(v[VX1] = u[MX1]*tau;  ,
           v[VX2] = u[MX2]*tau;  ,
           v[VX3] = u[MX3]*tau;)

    EXPAND(v[BX1] = u[BX1];  ,
           v[BX2] = u[BX2];  ,
           v[BX3] = u[BX3];)

    kinb2 = 0.5*(m2*tau + b2);    

  /* -- Check energy positivity -- */

#if HAVE_ENERGY
    if (u[ENG] < 0.0) {
      WARNING(
        print("! ConsToPrim: negative energy (%8.2e), ", u[ENG]);
        Where (i, NULL);
      )
      u[ENG]    = g_smallPressure/gmm1 + kinb2;
      flag[i] |= FLAG_CONS2PRIM_FAIL;
      ifail    = 1;
    }
#endif

  /* -- Compute pressure from total energy or entropy -- */

#if EOS == IDEAL
  #if ENTROPY_SWITCH
    use_entropy = (flag[i] & FLAG_ENTROPY);
    use_energy  = !use_entropy;
    if (use_entropy){
      rhog1 = pow(rho, gmm1);
      v[PRS] = u[ENTR]*rhog1; 
      if (v[PRS] < 0.0){
        WARNING(
               print("! ConsToPrim: negative p(S) (%8.2e, %8.2e), ", v[PRS], u[ENTR]);
               Where (i, NULL);
           )
        v[PRS]   = g_smallPressure;
        flag[i] |= FLAG_CONS2PRIM_FAIL;
        ifail    = 1;
      }
      u[ENG] = v[PRS]/gmm1 + kinb2; /* -- recompute energy -- */
    }
  #endif  /* ENTROPY_SWITCH */

    if (use_energy){
      v[PRS] = gmm1*(u[ENG] - kinb2);
      if (v[PRS] < 0.0){
        WARNING(
          print("! ConsToPrim: negative p(E) (%8.2e), ", v[PRS]);
Show(ucons,i);                
          Where (i, NULL);
        )
        v[PRS]    = g_smallPressure;
        u[ENG]    = v[PRS]/gmm1 + kinb2; /* -- recompute energy -- */
        flag[i] |= FLAG_CONS2PRIM_FAIL;
        ifail    = 1;
      }
      #if ENTROPY_SWITCH
      u[ENTR] = v[PRS]/pow(rho,gmm1);  /* -- Recompute entropy -- */
      #endif
    }

    #if NSCL > 0
    NSCL_LOOP(nv) v[nv] = u[nv]*tau;
    #endif    

#elif EOS == ISOTHERMAL

    #if NSCL > 0
    NSCL_LOOP(nv) v[nv] = u[nv]*tau;
    #endif    
 
#elif EOS == PVTE_LAW

  /* -- Convert scalars here since EoS may need ion fractions -- */

    #if NSCL > 0                       
    NSCL_LOOP(nv) v[nv] = u[nv]*tau;
    #endif    

    if (u[ENG] != u[ENG]){
      print("! ConsToPrim: NaN found\n");
      Show(ucons,i);
      QUIT_PLUTO(1);
    }
    rhoe  = u[ENG] - kinb2; 

    err = GetEV_Temperature (rhoe, v, &T);
    if (err){  /* If something went wrong while retrieving the  */
               /* temperature, we floor \c T to \c T_CUT_RHOE,  */
               /* recompute internal and total energies.        */
      T = T_CUT_RHOE;
      WARNING(  
        print ("! ConsToPrim: rhoe < 0 or T < T_CUT_RHOE; ");
        Where(i,NULL);
      )
      rhoe     = InternalEnergy(v, T);
      u[ENG]   = rhoe + kinb2; /* -- redefine total energy -- */
      flag[i] |= FLAG_CONS2PRIM_FAIL;
      ifail    = 1;
    }
    v[PRS] = Pressure(v, T);

#endif  /* EOS  */

    #ifdef GLM_MHD
    v[PSI_GLM] = u[PSI_GLM]; 
    #endif
  }
  return ifail;
}
