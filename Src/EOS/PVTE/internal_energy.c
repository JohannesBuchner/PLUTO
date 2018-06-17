/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Thermodynamic relations for the caloric EoS, \c e=e(T,rho) and 
         \c T=T(e,rho).
  
  Collect various functions involving computations between internal 
  energy, density and temperature,
  \f[
       e = e(T,\vec{X})\,,\qquad
       T = T(e,\vec{X})
  \f]
  where \c e is the specific internal energy, \c T is the temperature
  and \c X are the ionization fractions.
  The first relation is used to retrieve \c e as a function of 
  \c T and \c X by calling InternalEnergy() while the second one computes
  \c T from internal energy and fractions and it is handled by the
  function GetEV_Temperature():
  
  If no chemistry is used (\c NIONS==0 ), ionization fractions 
  are computed assuming equilibrium conditions and \c X=X(T,rho). 
  This requires the numerical inversion of the second equation
  <tt> T=(e,X(T,rho)) </tt> which can be carried out using either a root 
  finder algorithm (typically Brent's method) or, alternatively, lookup 
  table together with cubic/linear (default) or bilinear interpolation.

  For the root-finder version, a 2D table (::Trhoe_tab) giving
  <tt>T=T(rhoe,rho)</tt> is pre-computed in the function
  MakeEV_TemperatureTable()  and lookup table is used to narrow
  down the root interval using the values stored inside this table.

  In the tabulated EOS version, a different table (::rhoe_tab)
  <tt>rhoe=rhoe(T,rho)</tt> giving the internal energy at predefined values
  of temperature and density is pre-calculated in the function
  MakeInternalEnergyTable() so that internal energy and temperature can 
  be found by a combination of lookup table and bilinear (direct or inverse)
  interpolation.
  
  \author A. Mignone (mignone@ph.unito.it)\n
          B. Vaidya
  \date   7 Jan, 2015
*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h" 

#ifndef TV_ENERGY_TABLE_NX   
 #define TV_ENERGY_TABLE_NX    1024
#endif
#ifndef TV_ENERGY_TABLE_NY   
 #define TV_ENERGY_TABLE_NY    512
#endif

static double rhoeFunc(double, void *);
static double InternalEnergyDeriv (double rho, double T);

#if TV_ENERGY_TABLE == YES
static Table2D rhoe_tab; /**< A 2D table containing pre-computed values of 
                              internal energy stored at equally spaced node 
                              values of \c Log(T) and \c Log(rho) .*/

/* ********************************************************************* */
void MakeInternalEnergyTable()
/*!
 * Compute a 2D table of the internal energy <tt>rhoe(T,rho)</tt> as
 * function of temperature and density.
 * The grid is equally spaced in <tt>log10(T) (=x)</tt> and 
 *  <tt>log10(rho) (=y)</tt>.
 * Values are stored in the static table ::rhoe_tab.
 *
 *********************************************************************** */
{
  int i,j;
  double x, y, q;
  double T, rho, v[NVAR];

  print ("> MakeInternalEnergyTable(): Generating table (%d x %d points)\n",
           TV_ENERGY_TABLE_NX, TV_ENERGY_TABLE_NY);

  InitializeTable2D(&rhoe_tab, 1.0, 1.e8, TV_ENERGY_TABLE_NX, 
                               1.e-6, 1.e6, TV_ENERGY_TABLE_NY);  
  for (j = 0; j < rhoe_tab.ny; j++){
  for (i = 0; i < rhoe_tab.nx; i++){
    T   = rhoe_tab.x[i];
    rho = rhoe_tab.y[j];
    v[RHO] = rho;
    rhoe_tab.f[j][i] = InternalEnergyFunc(v,T);
  }}

  rhoe_tab.interpolation = SPLINE1;

/* ----------------------------------------------------------------
    Compute cubic spline coefficients
   ---------------------------------------------------------------- */

  static double *dfdx;

  if (dfdx == NULL) dfdx = ARRAY_1D(TV_ENERGY_TABLE_NX, double);

  for (j = 0; j < rhoe_tab.ny; j++){
    rho = rhoe_tab.y[j]; 
    for (i = 0; i < rhoe_tab.nx; i++){
      T       = rhoe_tab.x[i]; 
      dfdx[i] = InternalEnergyDeriv(rho, T);
    }
    if (rhoe_tab.interpolation == SPLINE1){
      MonotoneSplineCoeffs(rhoe_tab.x, rhoe_tab.f[j], dfdx, rhoe_tab.nx,
                           rhoe_tab.a[j], rhoe_tab.b[j], rhoe_tab.c[j],
                           rhoe_tab.d[j]);
    }else if (rhoe_tab.interpolation == SPLINE2){
      SplineCoeffs(rhoe_tab.x, rhoe_tab.f[j],
                   dfdx[0], dfdx[rhoe_tab.nx-1], rhoe_tab.nx,
                   rhoe_tab.a[j], rhoe_tab.b[j], rhoe_tab.c[j], rhoe_tab.d[j]);
    }
  }

/* let's plot the exact solutin and the spline using 20 points between
   adjacent nodes */
#if 0 
int ii;
double G, lnT, rhoe, rhoe_ex;
FILE *fp;

j = TV_ENERGY_TABLE_NY/2;

/* -- write exact internal energy, rhoe(T) -- */

fp  = fopen("rhoe.dat","w");
v[RHO] = rho = rhoe_tab.y[j];
for (i = 10; i < rhoe_tab.nx-11; i++){
  for (ii = 0; ii < 20; ii++){   // Sub grid 
    T   = rhoe_tab.x[i] + ii*rhoe_tab.dx[i]/20.0;
    rhoe_ex = InternalEnergyFunc(v,T);
    G       = FundamentalDerivative(v,T); 
    Table2DInterpolate(&rhoe_tab, T, rho, &rhoe);
    fprintf (fp, "%12.6e  %12.6e  %12.6e  %12.6e\n",T,rhoe_ex, rhoe, G);
  }
}
fclose(fp);
#endif


  FinalizeTable2D(&rhoe_tab);
  if (prank == 0) WriteBinaryTable2D("rhoe_tab.bin",&rhoe_tab);  
}

/* ********************************************************************* */
double InternalEnergyDeriv (double rho, double T)
/*!
 * Compute derivative of internal energy using numerical
 * differentiation
 *
 *********************************************************************** */
{
  double Epp, Emm, Ep, Em, dEdT;
  double dT = 1.e-3, v[256];

  v[RHO] = rho;
  Epp   = InternalEnergyFunc(v,T*(1.0 + 2.0*dT));
  Ep    = InternalEnergyFunc(v,T*(1.0 + 1.0*dT));
  Em    = InternalEnergyFunc(v,T*(1.0 - 1.0*dT));
  Emm   = InternalEnergyFunc(v,T*(1.0 - 2.0*dT));

  dEdT = (-Epp + 8.0*Ep - 8.0*Em + Emm)/(12.0*T*dT);
  return dEdT;
}


#endif /* TV_ENERGY_TABLE == YES */

#if NIONS == 0 && TV_ENERGY_TABLE == NO
static Table2D Trhoe_tab; /**< A 2D table containing pre-computed values of 
                               temperature by inverting e=e(T,rho) at
                               equally spaced node values of 
                               \c Log(rhoe) and \c Log(rho). This is used 
                               only in conjunction with a root finder to 
                               bracket the zero in a narrower interval. */
/* ********************************************************************* */
void MakeEV_TemperatureTable()
/*!
 * Pre-calculate a 2D table of temperature <tt>T(rhoe,rho)</tt> by 
 * inverting, at specified values of \c p and \c rho, the nonlinear 
 * equation 
 * \f[
 *   f(T) = \rho e - \rho e(T,\rho) = 0
 *   \qquad\Longrightarrow\qquad
 *   T_{ij} = T(x_i, y_j) \qquad
 *   \left(x=\log_{10}\rho e\,,\quad y=\log_{10}\rho\right)
 * \f]
 * For convenience, the table is constructed using equally 
 * spaced bins of <tt>x=Log(rhoe)</tt> and \c y=Log(rho)
 * where \c rhoe and \c rho are in code units. 
 * This function must be called only once to initialize the table 
 * \c ::Trhoe_tab which is private to this file only.
 *
 * This table is required only by the root-finder methods and it is
 * useless for the tabulated version (default).
 *
 *********************************************************************** */
{
  int i,j, status;
  double rhoe, rho, Tlo, Thi, T;
  struct func_param par;

  print ("> MakeEV_TemperatureTable(): Generating table...\n");

  InitializeTable2D(&Trhoe_tab, 1.e-9, 1.e9, 1200, 1.e-12, 1.e12, 1200);
  Tlo = 1.0;
  Thi = 1.e12; 

  for (j = 0; j < Trhoe_tab.ny; j++){
  for (i = 0; i < Trhoe_tab.nx; i++){
    
    par.rhoe   = rhoe = Trhoe_tab.x[i];
    par.v[RHO] = rho  = Trhoe_tab.y[j];

/*  status = Ridder(rhoeFunc, &param, Tlo, Thi, -1, 1.0e-12, &T);  */
    status = Brent(rhoeFunc, &par, Tlo, Thi, -1, 1.e-12, &T);
      
    if (status != 0) T = -1.0;

    Trhoe_tab.f[j][i] = T;
  }}


/* let's plot the exact solution and the fundamental derivative  */
#if 0 
int ii;
double G, lnT, v[256];
FILE *fp;

/* -- write exact internal energy, rhoe(T) -- */
fp  = fopen("rhoe.dat","w");
v[RHO] = rho = 1.245197084735e+00;
for (lnT = 2.0; lnT < 4.5; lnT += 1.e-3){
  T    = pow(10.0, lnT);
  rhoe = InternalEnergyFunc(v,T);
  G    = FundamentalDerivative(v,T); 
  fprintf (fp,"%12.6e  %12.6e  %12.6e\n",T,rhoe,G);
}
fclose(fp);
#endif

/* -------------------------------------------------------
    Compute forward differences 
   ------------------------------------------------------- */

  FinalizeTable2D(&Trhoe_tab);
  if (prank == 0)  WriteBinaryTable2D("Trhoe_tab.bin",&Trhoe_tab);  
}
#undef NRHO
#undef NRHOE
#endif /* NIONS == 0 && TV_ENERGY_TABLE == NO */

/* ********************************************************************* */
double InternalEnergy(double *v, double T)
/*!
 *  Compute and return the gas internal energy as a function of
 *  temperature and density or fractions:
 *  <tt> rhoe = rhoe(T,rho) </tt> or <tt> rhoe = rhoe(T,X) </tt>.
 *
 *  \param [in]   v     1D Array of primitive variables containing
 *                      density and species. Other variables are
 *                      ignored.
 *  \param [in]   T     Gas temperature
 *
 *  \return The gas internal energy (rhoe) in code units.
 *********************************************************************** */
{
#if TV_ENERGY_TABLE == YES

  int    status;
  double rho, rhoe;

  rho    = v[RHO];
  status = Table2DInterpolate(&rhoe_tab, T, rho, &rhoe);
  if (status != 0){
    print ("! InternalEnergy(): table interpolation failure (bound exceeded)\n");
    QUIT_PLUTO(1);
  }
  return rhoe;

#else 

  return InternalEnergyFunc (v,T);

#endif /* TV_ENERGY_TABLE == NO */
}

/* ********************************************************************* */
int GetEV_Temperature(double rhoe, double *v, double *T)
/*!
 *  Compute temperature from internal energy, density and fractions by
 *  solving either
 *  \f[
 *     \begin{array}{ll}
 *      \rho e - \rho e(T,\vec{X}) = 0  & \quad{(\rm with\; chemistry)}
 *            \\ \noalign{\medskip}
 *      \rho e - \rho e(T,\vec{X}(T,\rho)) = 0 & \quad{(\rm without\; chemistry)}
 *         \end{array}
 *  \f]
 *  
 *  The previous equations are solved using a root finder algorithm 
 *  (default is Brent or Ridder) or by lookup table / interpolation methods.
 *
 *  \param [in]  rhoe   Internal energy of gas in code units.
 *  \param [in]   v     1D array of primitive variables (only density
 *                      and species need to be set).
 *  \param [out]   T   Temperature in K
 *
 *  \return 0 on success, 1 on failure.
 *********************************************************************** */
{
#if TV_ENERGY_TABLE ==  YES
  int status;
  double rho;
  
  if (rhoe < 0.0) return 1;

  rho    = v[RHO];
  status = InverseLookupTable2D(&rhoe_tab, rho, rhoe, T);
  if (status != 0){
    WARNING(
      print ("! GetEV_Temperature(): table interpolation failure");
      print ("  [rho = %12.6e, rhoe = %12.6e]\n",rho,rhoe);
    )
    return 1;
  }  

  if (*T < T_CUT_RHOE) return 1;
  return 0;
#else

  int    nv, status, i, j;
  double Tmin, Tmax,rho, lnx, lny;
  struct func_param param;

  if (rhoe < 0.0) return 1;

/* ---------------------------------------------------
    Constrain species to lie between [0,1]
   --------------------------------------------------- */
 
#if NIONS > 0
  NIONS_LOOP(nv) {
    v[nv] = MAX(v[nv],0.0);
    v[nv] = MIN(v[nv],1.0);
  }
  #if COOLING == H2_COOL
    v[X_H2] = MIN(v[X_H2], 0.5);
  #endif
#endif
 
  param.v[RHO] = v[RHO];  /* maybe set param.v to point to v ?? */
  param.rhoe   = rhoe;
  NSCL_LOOP(nv) param.v[nv] = v[nv];

  #if NIONS == 0   /* No chemistry */

   rho = v[RHO];
   lnx = log10(rhoe);
   lny = log10(rho);
   i   = INT_FLOOR((lnx - Trhoe_tab.lnxmin)*Trhoe_tab.dlnx_1);
   j   = INT_FLOOR((lny - Trhoe_tab.lnymin)*Trhoe_tab.dlny_1);

/* ------------------------------------------------------------------
    Since Trhoe_tab[j][i] increases with i but decrease with j,
    the min and max between four adjacent node values are the two
    opposite points of the secondary diagonal (i,j+1) and (i+1,j),
    respectively.
   ------------------------------------------------------------------ */

   Tmin = Trhoe_tab.f[j+1][i];
   Tmax = Trhoe_tab.f[j][i+1];

  #else           /* Chemistry */
   InternalEnergyBracket(rhoe, v, &Tmin, &Tmax);
  #endif

  if (Tmin < 0.0 || Tmax < 0.0) return 1;

/*  status = Ridder(rhoeFunc, &param, Tmin, Tmax,-1, 1.0e-7, &T);  */
  status = Brent(rhoeFunc, &param, Tmin, Tmax, -1,   1.e-7,  T);

  if (*T < T_CUT_RHOE || status != 0) return 1;
  return 0;

#endif /* TV_ENERGY_TABLE == NO */
}
/* ********************************************************************* */
double rhoeFunc(double T, void *par)
/*! 
 * Return the nonlinear function <tt>f(T)=rhoe(T,rho)-rhoe </tt> used 
 * by the root finder to find temperature \c T.
 *
 * \param [in] T    The temperature in Kelvin
 * \param [in] par  A pointer to a \c func_param structure 
 *                  containing quantities required for the computation.
 * \return f(T)      
 *********************************************************************** */
{
  double func, rhoe;
  struct func_param *p = (struct func_param *) par;

  if (T < 0.0){
    printf ("! rhoeFunc: T = %12.6e\n",T);
    exit(1);
  }
  rhoe = InternalEnergyFunc(p->v, T);
  func = rhoe - p->rhoe;
  return func;
}

/* ********************************************************************* */
void SoundSpeed2 (const State *p, int beg, int end,int pos, Grid *grid)
/*!
 * Define the square of the sound speed.
 * 
 * \param [in]    v   1D array of primitive quantities
 * \param [out] cs2   1D array containing the square of the sound speed
 * \param [in]    h   1D array of enthalpy values
 * \param [in]  beg   initial index of computation 
 * \param [in]  end   final   index of computation
 * \param [in]  pos   an integer specifying the spatial position 
 *                    inside the cell (only for spatially-dependent EOS)
 * \param [in]  grid  pointer to an array of Grid structures
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int  i;

  for (i = beg; i <= end; i++){
    p->a2[i] = Gamma1(p->v[i])*p->v[i][PRS]/p->v[i][RHO];
/*     cs2[i] = 1.6667*v[i][PRS]/v[i][RHO];  */
  }
}
