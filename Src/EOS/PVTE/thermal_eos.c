/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Thermodynamic relations for the thermal EoS, \c p=nkT and
         <tt> T=p/(nK) </tt>.
  
  Collect miscellaneous functions involving computations with the 
  thermal equation of state,
  \f[
       p = nk_BT \quad {\rm(in\; cgs)} 
       \qquad \Longrightarrow\qquad
       p = \frac{\rho T}{K\mu(\vec{X})} \quad {\rm(in\; code\quad units)} 
  \f]
  where, in the second expression, \c rho and \c p are density and 
  pressure (in code units), \c T is the temperature (in Kelvin), \c K 
  is the \c ::KELVIN dimensionalization constant, \c mu(\b X) is the mean 
  molecular weight and \b X is an array containing the ionization fractions 
  of different elements.
  The previous relation is typically used to compute
  - <tt> p=p(T,rho,\b X) </tt> by calling Pressure(); 
  - <tt> T=T(p,rho,\b X) </tt> by calling GetPV_Temperature().
    This relation can be either:
    - \e Explicit: ion fractions are evolved independently using a 
      chemical reaction network (e.g. when  \c H2_COOL or \c SNEq is enabled).
      Then one has \f$ T=(p/\rho) K\mu(\vec{X})\f$.  
    - \e Implicit: ions are not evolved explicitly but are computed 
      from equilibrium considerations (e.g. Saha of collisional ionization 
      equilibrium) so that \c X=X(T,rho) and the mean molecular 
      weight now depends on temperature.
      Then one has \f$T = (p/\rho) K \mu(T,\rho)\f$.

  Nonlinear equation can be solved with a root finder (typically Brent's 
  method) or a tabulated approach (default).
  In both case, a 2D table (::Ttab) of the temperature \c T(i,j) is
  pre-calculated for fixed values of \c p/rho and \c rho by
  MakePV_TemperatureTable() and stored into memory.
  The table is then used to bracket the root (in the case of a root 
  finder) or to replace the runtime computation with a simpler 
  array indexing operation followed by a combination of lookup table 
  and bilinear (direct or inverse) interpolation.
  
  \author A. Mignone (mignone@ph.unito.it)\n
          B. Vaidya
  \date   31 Aug, 2014
*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h" 

#ifndef PV_TEMPERATURE_TABLE_NX
 #define PV_TEMPERATURE_TABLE_NX     768
#endif
#ifndef PV_TEMPERATURE_TABLE_NY
 #define PV_TEMPERATURE_TABLE_NY     768
#endif

static Table2D  Ttab; /**< A 2D table containing pre-computed values of 
                           temperature stored at equally spaced node values 
                           of <tt> Log(p/rho) </tt> and \c Log(rho)  */
#if NIONS == 0
static double TFunc(double T, void *par);

/* ********************************************************************* */
void MakePV_TemperatureTable()
/*!
 * Pre-calculate a 2D table of temperature <tt>T(p,rho)</tt> by 
 * inverting, at specified values of \c p and \c rho, the nonlinear 
 * equation (only in LTE or CIE)
 * \f[
 *   f(T) = T - \frac{p}{\rho} K\mu(T,\rho) = 0
 *   \qquad\Longrightarrow\qquad
 *   T_{ij} = T(x_i, y_j) \qquad
 *   \left(x=\log_{10}\frac{p}{\rho}\,,\quad y=\log_{10}\rho\right)
 * \f]
 * For convenience, the table is constructed using equally 
 * spaced bins of <tt>x=log10(p/rho)</tt> and <tt>y=log10(rho)</tt> 
 * where \c p and \c rho are in code units. 
 * This function must be called only once to initialize the 2D table 
 * \c ::Ttab.
 *
 *********************************************************************** */
{
  int i,j, status;
  double prs, rho, T1, Tlo, Thi, T, logK = log10(KELVIN);
  double mu_lo, mu_hi;
  struct func_param par;

  print1 ("> MakePV_TemperatureTable: Generating table...\n");

/* --------------------------------------------------------------
    Initialize table. The two table axis are given by 
    ln(x) = log(p/rho) and ln(y) = log(rho). 
    The lower and upper x-bounds correspond to T/mu = 1 K and 
    T/mu = 10^7 K, respectively.
   -------------------------------------------------------------- */

  InitializeTable2D(&Ttab,1.0/KELVIN, 1.e7/KELVIN, PV_TEMPERATURE_TABLE_NX, 
                          1.e-7     , 1.e7, PV_TEMPERATURE_TABLE_NY);

/* -----------------------------------------------------------------
    Guess the smallest and largest value of \mu. 
    This is likely to happen at very low and very high temperatures  
   ----------------------------------------------------------------- */

  GetMu(1.0,   1.0, &mu_lo);
  GetMu(1.e12, 1.0, &mu_hi);

  for (j = 0; j < Ttab.ny; j++){
  for (i = 0; i < Ttab.nx; i++){

    T1  = Ttab.x[i]*KELVIN;  /* = p/rho*KELVIN */
    rho = Ttab.y[j];

    par.T1     = T1;
    par.v[RHO] = rho; 

  /* -- set interval endpoints for root finder -- */ 
    
    Tlo = par.T1*mu_lo;  /*  T > T1*mu_\min  */
    Thi = par.T1*mu_hi;  /*  T < T1*mu_\max  */
     
/*    status = Ridder(TFunc, &par, Tlo, Thi, -1, 1.0e-12, &T); */
    status = Brent(TFunc, &par, Tlo, Thi, -1, 1.0e-12, &T);    
    if (status != 0) {
      print1 ("! MakePV_TemperatureTable: ");
      print1 ("root could not be found [%d]\n",status);
      QUIT_PLUTO(1);
    }
    Ttab.f[j][i] = T;
  }}
  
  FinalizeTable2D(&Ttab);
  if (prank == 0)  WriteBinaryTable2D("T_tab.bin",&Ttab);  

#if 0 /* Table-to-Table conversion: attempt to interpolate from another table */
{
  int status;
  double T, p, mu = 1.2;
  Table2D rhotab;  /* test table rho = rho(T,p);  */

  InitializeTable2D(&rhotab, 1.e2, 1.e7, 64, 1.e-2, 1.e4, 80);

  for (j = 0; j < rhotab.ny; j++){
  for (i = 0; i < rhotab.nx; i++){
    T = rhotab.x[i];
    p = rhotab.y[j];
mu = (T/5.e3);
mu = 1.2 + 0.0/cosh(mu*mu);
    rhotab.f[j][i] = p*KELVIN*mu/T;   /* rho */
  }}

  FinalizeTable2D(&rhotab);
  if (prank == 0)  WriteBinaryTable2D("rhotab.bin",&rhotab);

/* -- Build T(p/rho, rho) table -- */

  status = 0;
  for (j = 0; j < Ttab.ny; j++){
  for (i = 0; i < Ttab.nx; i++){
    rho = Ttab.y[j];
    p   = Ttab.x[i]*rho;
    status = InverseLookupTable2D(&rhotab, p, rho, &T);
    if (status != 0){
      printf ("! Error in inverting table\n");
      T = -1.0;
    }
    Ttab.f[j][i] = T;
  }}

  if (prank == 0)  WriteBinaryTable2D("T_tab2.bin",&Ttab);

  printf ("KELVIN = %12.6e\n",KELVIN);
  exit(1); 
}
#endif /* 0 */
}
#endif  /* NIONS == 0 */

/* ********************************************************************* */
double Pressure(double *v, double T)
/*!
 * Return pressure as a function of temperature, density and fractions.
 * If <tt>T(p/rho,rho)</tt> has been tabulated, this requires 
 * inverting the table by calling InverseLookupTable2D().
 *
 * \param [in]  v     a pointer to a vector of primitive quantities 
 *                    (only density and fractions will be actually used).
 * \param [in]  T     The temperature (in Kelvin);
 *
 * \return Pressure in code units.                   
 *********************************************************************** */
{
#if PV_TEMPERATURE_TABLE == YES 

  int    status;
  double rho, T1;

  rho = v[RHO];
  status = InverseLookupTable2D(&Ttab, rho, T, &T1);
  if (status != 0){
    print ("! Pressure: table interpolation failure [rho = %12.6e]\n", rho);
    QUIT_PLUTO(1);
  }
  
  return T1*rho;

#else

  double mu;

  #if NIONS == 0
   GetMu(T, v[RHO], &mu);
  #else
   mu = MeanMolecularWeight(v);
  #endif
  return v[RHO]*T/(KELVIN*mu);

#endif
}

/* ********************************************************************* */
int GetPV_Temperature (double *v, double *T)
/*!
 *  Compute gas temperature (in Kelvin) from density, pressure 
 *  and fractions:
 *  \f[
 *     T = \left\{\begin{array}{ll}
 *         \DS \frac{p}{\rho} K \mu(\vec{X}) &\quad{(\rm with\; chemistry)}
 *            \\ \noalign{\medskip}
 *         \DS \frac{p}{\rho}K\mu(\vec{X}(T,\rho)) & \quad{(\rm without\; chemistry)}
 *         \end{array}\right.
 *  \f]
 *  where \c K is the \c ::KELVIN macro.
 *
 *  The first relation is used when ion fractions are explicitly evolved
 *  by the code (\c NIONS>0) while the second one is inverted numerically
 *  using root finder or lookup table when either LTE or CIE is assumed.
 *
 *  \param [in]    v   1D array of primitive quantities.
 *  \param [out]   T   Temperature in K
 *
 *  \return 0 on success, 1 on failure.
 *********************************************************************** */
{
  int status;

#if NIONS == 0 && PV_TEMPERATURE_TABLE == YES

  double rho, T1;
  
  if (v[PRS] < 0.0) return 1;

  rho    = v[RHO];
  T1     = v[PRS]/rho;
  status = Table2DInterpolate(&Ttab, T1, rho, T);

  if (status != 0){
    print ("! GetPV_Temperature: table interpolation failure ");
    print ("[T1 = %8.3e, rho = %8.3e]\n",T1,rho);
    
    if (status == -1) return 1;      /* hit lower x range (easy fix) */
    else              QUIT_PLUTO(1); /* not so easy to fix */
  }

  if (*T < T_CUT_RHOE) return 1;
  else                 return 0;

#elif NIONS == 0 && PV_TEMPERATURE_TABLE == NO 

  int    i,j;
  double Tmin, Tmax, lnx, lny, rho, prs;
  struct func_param par;

  if (v[PRS] < 0.0) return 1;

  par.v[RHO] = v[RHO];
  par.T1     = v[PRS]/v[RHO]*KELVIN;

/* -- Locate indices (i,j) in the table -- */

  rho = v[RHO];
  prs = v[PRS];
  lnx = log10(prs/rho);
  lny = log10(rho);
  i   = INT_FLOOR((lnx - Ttab.lnxmin)*Ttab.dlnx_1);
  j   = INT_FLOOR((lny - Ttab.lnymin)*Ttab.dlny_1);

/* -------------------------------------------------------------------
    Since Ttab[j][i] increases with i and j, the min and max between
    four adjacent node values are the two opposite points on the main
    diagonal given, respectively, by (i,j) and (i+1,j+1).
   ------------------------------------------------------------------- */

  Tmin = MIN(Ttab.f[j][i], Ttab.f[j+1][i+1]);
  Tmax = MAX(Ttab.f[j][i], Ttab.f[j+1][i+1]);

/*  status = Ridder(TFunc, &p, Tmin, Tmax, -1, 1.0e-7, &T);  */
  status = Brent(TFunc, &par, Tmin, Tmax, -1, 1.0e-7, T);
  if (status != 0){ 
    print ("! PrimitiveTemperature: could not find root, ");
    print ("! [Tmin, Tmax] = [%12.6e, %12.6e]\n",Tmin, Tmax);
    QUIT_PLUTO(1);
  }

  if (*T < T_CUT_RHOE || status != 0) return 1;
  return 0;
  
#else

  *T = (v[PRS]/v[RHO])*KELVIN*MeanMolecularWeight(v);
  if (*T < T_CUT_RHOE) return 1;
  return 0;
   
#endif 

}

#if NIONS == 0
/* ********************************************************************* */
double TFunc(double T, void *par)
/*!
 * Return the nonlinear function <tt> f(T) = T - T1 mu(T,rho)</tt> used 
 * by the root finder to find temperature \c T.
 *
 * \param [in] T    The temperature in Kelvin
 * \param [in] par  A pointer to a \c func_param structure containing 
 *                  density required for computing internal energy.
 * \return f(T)
 *********************************************************************** */
{
  double  T1, rho, mu;
  struct func_param *p = (struct func_param *) par;
  
  T1   = p->T1;
  rho  = p->v[RHO];
  GetMu(T, rho, &mu);  
  return T - T1*mu;
}
#endif

