/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief PVTE EoS header file.

  Contains basic macro definitions, function prototypes and structure 
  definition used by the PVTE_LAW equation of state.
  In this module, two equations of state are used:
  \f[
  \left\{\begin{array}{lcll}
  p &=&\DS \frac{\rho T}{{\rm K}\mu(\vec{X})} &\qquad\rm{(thermal)}\,,
    \\ \noalign{\bigskip}
  e &=&    e(T,\vec{X}) &\qquad\rm{(caloric)}\,.
  \end{array}\right.
  \f]

  Since these two equations can be nonlinear functions and PLUTO performs
  conversions between primitive and conservative variables quite often
  during a single update step, the employment of this EoS is more
  computationally intensive:

  -# *Non equilibrium case*: 
     when converting from primitive to conservative and
     viceversa, the following operations are performed:
     \f[ \begin{array}{lll}
         (p,\rho,\vec{X}) &\quad\Longrightarrow\quad \DS T = \frac{p}{\rho}K\mu(\vec{X}) 
           &\quad\Longrightarrow\quad e = e(T,\vec{X})   
           \\ \noalign{\medskip}
       (e,\rho,\vec{X}) &\quad\Longrightarrow\quad T = T(e,\vec{X})^* 
           &\quad\Longrightarrow\quad \DS p = \frac{T\rho}{K\mu(\vec{X})} 
           \\ \noalign{\medskip}
         \end{array}
     \f]
     The \* symbol means that a nonlinear equation must be inverted
     numerically using a root-finder.
  -# *Equilibrium case*: \b X \c (T,rho) and all the
     dependencies on \b X become dependencies on  \c (T,rho).
     In this case the operations carried out during primitive to conservative
     and viceversa are:
     \f[ \begin{array}{lll}
     (p,\rho) &\quad\Longrightarrow\quad \DS T = \frac{p}{\rho}K\mu(T,\rho)^*
             &\quad\Longrightarrow\quad e = e(T,\rho)
           \\ \noalign{\medskip}
  (e,\rho) &\quad\Longrightarrow\quad T = T(e,\rho)^*
           &\quad\Longrightarrow\quad \DS p = \frac{T\rho}{K\mu(T,\rho)} 
           \\ \noalign{\medskip}
 \end{array}
   \f]
       requiring thus two numerical inversions of nonlinear equations \*.

  The sequence of operation is handled by the following functions:

  <CENTER>
  Main variables | Calling Function | T(p,rho,X)          | e = e(rho,T,X) 
  -------------- | ---------------- | ------------------- | --------------
  (p,rho,X)      | PrimToCons()     | GetPV_Temperature() | InternalEnergy()

  Main variables | Calling Function | T(rho,e,X)          | p = (rho,T,X)
  -------------- | ---------------- | ------------------- | ------------
  (rho,e,X)      | ConsToPrim()     | GetEV_Temperature() | Pressure()
  </CENTER>

  \authors B. Vaidya (bhargav.vaidya@ph.unito.it)\n
           A. Mignone (mignone@ph.unito.it)
  \date    June 25, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#ifndef T_CUT_RHOE
 #define T_CUT_RHOE  10.0  /**< Sets the lowest cut-off temperature used in 
                                the PVTE_LAW equation of state.
                                Zones with temperature below T_CUT_RHOE
                                will be reset to this value and the internal
                                energy will be redefined accordingly. */
#endif


#if ENTROPY_SWITCH
 #error ! PVTE_LAW not working with ENTROPY_SWITCH
#endif


/* ***********************************************************
    \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES 
    Function prototyping
   *********************************************************** */

double FundamentalDerivative(double *, double T);
double Gamma1(double *);   /* User supplied */
int    GetEV_Temperature (double, double *, double *);
void   GetMu     (double, double, double *);     /* User supplied */
int    GetPV_Temperature (double *, double *);
double InternalEnergy (double *, double);     
double InternalEnergyFunc (double *, double);    /* User supplied */
void   InternalEnergyBracket (double, double *, double *, double *); 
void   MakeInternalEnergyTable();
void   MakePV_TemperatureTable();
void   MakeEV_TemperatureTable();
double Pressure(double *, double);
/* \endcond */

struct func_param {
  double v[NVAR];
  double rhoe;
  double T1;
};

/* The following set of switches should be enabled ONLY WITHOUT CHEMISTRY */
/* Use YES/NO to switch between the version with root finder in both   
   temperature and internal energy or faster version with lookup table   */

#if NIONS == 0
 #ifndef PV_TEMPERATURE_TABLE
  #define PV_TEMPERATURE_TABLE   YES
 #endif
 #ifndef TV_ENERGY_TABLE
  #define TV_ENERGY_TABLE        YES
 #endif
#else
 #undef PV_TEMPERATURE_TABLE
 #undef TV_ENERGY_TABLE
 #define PV_TEMPERATURE_TABLE   NO
 #define TV_ENERGY_TABLE        NO
#endif

