/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief PVTE_LAW for a partially ionized gas.

  Compute the internal energy for a partially ionized hydrogen gas.

  \author A. Mignone (mignone@ph.unito.it)\n
          B. Vaidya
  \date   29 Aug 2014

  \b Reference 
     - PLUTO Users' guide.
*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h" 

#define INTE_EXACT   (TV_ENERGY_TABLE == 1 ? 0:1)
#define SEC_ORDER 0

/* ****************************************************************** */
double FundamentalDerivative(double *v, double T)
/*!
 * To compute the fundamental derivative 
 * for each point. 
 *
 *
 ******************************************************************** */
{  
  int nv;
  double FDeriv, cs, Gam1, Gam2, Gam3;
  double prs, rho, e,  dT, drho, rho2, rho3, rho4;
  double dPdr, d2Pdr2, dPdt, d2Pdt2, d2Pdrdt, dEdt, d2Edt2;
  
  double prs_rmm, prs_rpp, prs_rm, prs_rp;
  double prs_tmm, prs_tpp, prs_tm, prs_tp; 
  double prs_rm_tm, prs_rp_tm, prs_rm_tp, prs_rp_tp;
  double prs_rmm_tmm, prs_rpp_tpp, prs_rmm_tpp, prs_rpp_tmm;
  double prs_rp_tmm, prs_rpp_tm, prs_rmm_tp, prs_rm_tpp;
  double prs_rm_tmm, prs_rmm_tm, prs_rp_tpp, prs_rpp_tp;
  
  double epp, emm, ep, em;
  double vm[NVAR], vp[NVAR], vpp[NVAR], vmm[NVAR];
 
  VAR_LOOP(nv) vm[nv] = vp[nv] = vpp[nv] = vmm[nv] = v[nv];

  rho  = v[RHO];
  rho2 = rho*rho; 
  rho3 = rho2*rho;
  rho4 = rho2*rho2;
  
#if INTE_EXACT
  e = InternalEnergyFunc(v, T)/rho;
#else
  e =  InternalEnergy(v, T)/rho;
#endif

  prs = Pressure(v, T);
  v[PRS] = prs;
  cs  = Gamma1(v)*prs/rho;

  dT = 1.0e-3;
  drho = 1.0e-3;

/* Compute fist and second derivative of Pressure w.r.t Density */ 


  vmm[RHO] = rho*(1.0 - 2.0*drho);
  vm[RHO]  = rho*(1.0 - 1.0*drho);
  vp[RHO]  = rho*(1.0 + 1.0*drho);
  vpp[RHO] = rho*(1.0 + 2.0*drho);  
  
  prs_rmm = Pressure(vmm, T);
  prs_rm  = Pressure(vm, T);
  prs_rp  = Pressure(vp, T);
  prs_rpp = Pressure(vpp, T);
  
#if SEC_ORDER
  dPdr = prs_rp - prs_rm;
  dPdr /= 2.0*rho*drho;
  
  d2Pdr2 = prs_rp - 2.0*prs + prs_rm;
  d2Pdr2 /= rho2*drho*drho;
#else
  dPdr =  prs_rmm - 8.0*prs_rm + 8.0*prs_rp - prs_rpp ;
  dPdr /= 12.0*rho*drho;
  
  d2Pdr2 = -prs_rmm + 16.0*prs_rm - 30.0*prs + 16.0*prs_rp - prs_rpp;
  d2Pdr2 /= 12.0*rho2*drho*drho;
#endif

 /* Compute fist and second derivative of Pressure w.r.t Temperature */ 
  
  prs_tmm = Pressure(v, T*(1.0 - 2.0*dT));
  prs_tm  = Pressure(v, T*(1.0 - 1.0*dT));
  prs_tp  = Pressure(v, T*(1.0 + 1.0*dT));
  prs_tpp = Pressure(v, T*(1.0 + 2.0*dT));

#if SEC_ORDER
  dPdt = prs_tp - prs_tm;
  dPdt /= 2.0*T*dT;
  
  d2Pdt2 = prs_tp - 2.0*prs + prs_tm;
  d2Pdt2 /= T*T*dT*dT;
#else
  dPdt = prs_tmm - 8.0*prs_tm + 8.0*prs_tp - prs_tpp;
  dPdt /= 12.0*T*dT;
  
  d2Pdt2 = -prs_tmm + 16.0*prs_tm - 30.0*prs + 16.0*prs_tp - prs_tpp;
  d2Pdt2 /= 12.0*T*T*dT*dT;
#endif


  /* Compute fist and second derivative of Internal energy w.r.t Temperature */ 

#if INTE_EXACT
  emm  = InternalEnergyFunc(v, T*(1.0 - 2.0*dT))/rho;
  em   = InternalEnergyFunc(v, T*(1.0 - 1.0*dT))/rho;
  ep   = InternalEnergyFunc(v, T*(1.0 + 1.0*dT))/rho;
  epp  = InternalEnergyFunc(v, T*(1.0 + 2.0*dT))/rho;
#else
  emm  = InternalEnergy(v, T*(1.0 - 2.0*dT))/rho;
  em   = InternalEnergy(v, T*(1.0 - 1.0*dT))/rho;
  ep   = InternalEnergy(v, T*(1.0 + 1.0*dT))/rho;
  epp  = InternalEnergy(v, T*(1.0 + 2.0*dT))/rho;
#endif

#if SEC_ORDER
  dEdt  = ep - em;
  dEdt /= 2.0*T*dT;
  
  d2Edt2 = ep - 2.0*e + em;
  d2Edt2 /= T*T*dT*dT;
#else
  dEdt  = emm - 8.0*em + 8.0*ep - epp;
  dEdt /= 12.0*T*dT;
  
  d2Edt2 = -emm + 16.0*em -30.0*e + 16.0*ep - epp;
  d2Edt2 /= 12.0*T*T*dT*dT;
#endif

 /* Compute mixed derivative of Pressure w.r.t Temperature and Density */ 

  prs_rm_tm  = Pressure(vm, T*(1.0 - 1.0*dT));
  prs_rp_tm  = Pressure(vp, T*(1.0 - 1.0*dT));
  prs_rm_tp  = Pressure(vm, T*(1.0 + 1.0*dT));
  prs_rp_tp  = Pressure(vp, T*(1.0 + 1.0*dT));
 
  prs_rmm_tmm  = Pressure(vmm, T*(1.0 - 2.0*dT));
  prs_rpp_tmm  = Pressure(vpp, T*(1.0 - 2.0*dT));
  prs_rmm_tpp  = Pressure(vmm, T*(1.0 + 2.0*dT));
  prs_rpp_tpp  = Pressure(vpp, T*(1.0 + 2.0*dT));

  prs_rp_tmm  = Pressure(vp, T*(1.0 - 2.0*dT));
  prs_rpp_tm  = Pressure(vpp, T*(1.0 - 1.0*dT));
  prs_rmm_tp  = Pressure(vmm, T*(1.0 + 1.0*dT));
  prs_rm_tpp  = Pressure(vm, T*(1.0 + 2.0*dT));

  prs_rm_tmm  = Pressure(vm, T*(1.0 - 2.0*dT));
  prs_rmm_tm  = Pressure(vmm, T*(1.0 - 1.0*dT));
  prs_rp_tpp  = Pressure(vp, T*(1.0 + 2.0*dT));
  prs_rpp_tp  = Pressure(vpp, T*(1.0 + 1.0*dT));



#if SEC_ORDER
  d2Pdrdt = (prs_rp_tp + prs_rm_tm - prs_rp_tm - prs_rm_tp);
  d2Pdrdt /= 4.0*T*dT*rho*drho;
#else
  d2Pdrdt = 74.0*(prs_rp_tp + prs_rm_tm - prs_rp_tm - prs_rm_tp)	\
    + 44.0*(prs_rpp_tpp + prs_rmm_tmm - prs_rpp_tmm - prs_rmm_tpp)	\
    - 63.0*(prs_rp_tmm + prs_rpp_tm + prs_rmm_tp + prs_rm_tpp)		\
    + 63.0*(prs_rm_tmm + prs_rmm_tm + prs_rp_tpp + prs_rpp_tp);
  
  d2Pdrdt /= 600.0*T*dT*rho*drho;
#endif

/* Compute three terms for the fundamental Gas derivative */ 

  Gam1 = rho4*d2Pdr2 + 2.0*rho2*dPdr;
  Gam2 = 3.0*T*(rho2/dEdt)*dPdt*d2Pdrdt;
  Gam3 = (T/dEdt)*(T/dEdt)*dPdt*dPdt;
  Gam3 *= 3.0*d2Pdt2 + (1./T)*dPdt*(1.0 - (T/dEdt)*d2Edt2);

  FDeriv = (Gam1 + Gam2 + Gam3)/(2.0*cs*cs*rho3);

  return FDeriv;
}
