/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Reconstruction coefficients header file

  Define some useful macros for the computation of slope limiters 
  in the piecewise linear method.
  
  The macro ::UNIFORM_CARTESIAN_GRID can be set to YES to 
  enable faster computation when the grid is uniform and Cartesian.
  In the general case (non-uniform and/or non-Cartesian), set it to NO. 
  
  \authors A. Mignone (mignone@ph.unito.it)
  \date    May 19, 2014
  
  \b References
     - "High-order conservative reconstruction schemes for finite
        volume methods in cylindrical and spherical coordinates"
        A. Mignone, JCP (2014), 270, 784.
*/
/* ///////////////////////////////////////////////////////////////////// */

#ifndef UNIFORM_CARTESIAN_GRID
 #if GEOMETRY == CARTESIAN 
  #define UNIFORM_CARTESIAN_GRID   YES
 #else
  #define UNIFORM_CARTESIAN_GRID   NO
 #endif
#endif

#define CHECK_MONOTONICITY       NO

/* ********************************************************************* */
/*! Simple structure used to retrieve 1D reconstruction weights 
    (c, w, d) used by piecewise linear interpolation (see states_plm.c)
   ********************************************************************* */

typedef struct PLM_COEFFS{
  double *cp;
  double *cm;
  double *wp;
  double *wm;
  double *dp;
  double *dm;
} PLM_Coeffs;

void PLM_CoefficientsSet(Grid *grid);
void PLM_CoefficientsGet(PLM_Coeffs*, int);

/* ---------------------------------------------------------------------
     Set macro for limiter functions
   --------------------------------------------------------------------- */

/*! \name Limiter macros.
    The following macros define a few limiter functions that can be 
    used for piecewise linear reconstruction.
    Each limiter has the form \c SET_XX_LIMITER(dv, dvp, dvm, cp, cm) 
    where \c dvp and \c dvm are the geometrically weighted forward and 
    backward slopes (\f$\Delta Q^{F}\f$ and \f$\Delta Q^{B}\f$ in Eq. [29]),
    \c cp and \c cm are weight coefficients depending on the distance 
    between the cell interfaces and the centroid of volume (Eq. [33] in 
    the reference paper). 
    The limited slope (used in Eq. [30]) is returned in \c dv.
    
    Usually \c cp and \c cm are used only for a few limiters and 
    when the grid is either non-uniform or the geometry is not Cartesian.
    For this reason, the OSPRE (OS), van Leer (VL) and monotonized central
    (MC) are 
*/

/*! Set flat (zero slope) reconstruction (also non-uniform). */
#define SET_FL_LIMITER(dv, dvp, dvm, cp, cm) \
  dv = 0.0                                  

/*! Set Minmod limiter. */
#define SET_MM_LIMITER(dv, dvp, dvm, cp, cm) \
  dv = (dvp*dvm > 0.0 ? ABS_MIN(dvp, dvm): 0.0)

/*! Van Albada limiter */
#define SET_VA_LIMITER(dv, dvp, dvm, cp, cm)\
  if (dvp*dvm > 0.0) { \
    double _dpp= dvp*dvp, _dmm = dvm*dvm; \
    dv = (dvp*(_dmm + 1.e-18) + dvm*(_dpp + 1.e-18))/(_dpp + _dmm + 1.e-18); \
  }else dv = 0.0;

/*! Umist limiter */
#define SET_UM_LIMITER(dv, dvp, dvm, cp, cm)\
  if (dvp*dvm > 0.0){ \
    double _ddp = 0.25*(dvp + 3.0*dvm), _ddm = 0.25*(dvm + 3.0*dvp); \
    double _d2  = 2.0*ABS_MIN(dvp, dvm);\
    _d2  = ABS_MIN(_d2, _ddp);\
    dv   = ABS_MIN(_d2, _ddm);\
  }else dv = 0.0;

/*! Generalised minmod limiter */
#define SET_GM_LIMITER(dv, dvp, dvm, cp, cm) \
   if (dvp*dvm > 0.0) { \
     double _qc  = 0.5*(dvm + dvp), _scrh = ABS_MIN((dvp)*(cp), (dvm)*(cm)); \
     dv   = ABS_MIN(_qc, _scrh);  \
   }else dv = 0.0; 

/* -------------------------------------------------------------
           Limiters on uniform Cartesian grid
   ------------------------------------------------------------- */
   
#if UNIFORM_CARTESIAN_GRID == YES 

/*! OSPRE limiter (uniform Cart. grid) */
 #define SET_OS_LIMITER(dv, dvp, dvm, cp, cm)\
  dv = (dvp*dvm > 0.0? \
       dv = 1.5*dvp*dvm*(dvm + dvp)/(dvp*dvp + dvm*dvm + dvp*dvm): 0.0);

/*! Van Leer limiter (uniform Cartesian grid) */
 #define SET_VL_LIMITER(dv, dvp, dvm, cp, cm)\
   dv = (dvp*dvm > 0.0 ? 2.0*dvp*dvm/(dvp + dvm) :0.0)

/*! Monotonized central limiter (uniform cart. grid). 
    Here \c cp and \c cm are useless. */
 #define SET_MC_LIMITER(dv, dvp, dvm, cp, cm) \
   if (dvp*dvm > 0.0) { \
     double _qc  = 0.5*(dvm + dvp), _scrh = 2.0*ABS_MIN(dvp, dvm); \
     dv   = ABS_MIN(_qc, _scrh);  \
   }else dv = 0.0; 

/* -------------------------------------------------------------
          Limiters on irregular or non-Cartesian grids
   ------------------------------------------------------------- */
   
#else 

/*! OSPRE limiter (general grid case) */
 #define SET_OS_LIMITER(dv, dvp, dvm, cp, cm)\
  if (dvp*dvm > 0.0){  \
    double _den = 2.0*dvp*dvp + 2.0*dvm*dvm + (cp + cm - 2.0)*dvp*dvm;\
    dv = dvp*dvm*((1.0+cp)*dvm + (1.0+cm)*dvp)/_den; \
  }else dv = 0.0;

/* -- van Leer limiter (general grid) -- */

 #define SET_VL_LIMITER(dv, dvp, dvm, cp, cm)\
   dv = (dvp*dvm > 0.0 ? dvp*dvm*(cp*dvm + cm*dvp) \
                       /(dvp*dvp + dvm*dvm + (cp + cm - 2.0)*dvp*dvm) :0.0)

/* -- monotonized central (general grid) -- */

 #define SET_MC_LIMITER(dv, dvp, dvm, cp, cm) \
   if (dvp*dvm > 0.0) { \
     double _qc  = 0.5*(dvm + dvp), _scrh = ABS_MIN(dvp*cp, dvm*cm); \
     dv   = ABS_MIN(_qc, _scrh);  \
   }else dv = 0.0; 

#endif

/* -------------------------------------------------------------------
    when a single limiter is specified, use SET_LIMITER as
    a shortcut to the actual limiter macro
   ------------------------------------------------------------------- */

#ifdef LIMITER  /* May not be defined when using finite difference schemes */
 #if LIMITER == FLAT_LIM
  #define SET_LIMITER SET_FL_LIMITER
 #elif LIMITER == MINMOD_LIM
  #define SET_LIMITER  SET_MM_LIMITER
 #elif LIMITER == VANALBADA_LIM
  #define SET_LIMITER  SET_VA_LIMITER
 #elif LIMITER == OSPRE_LIM
  #define SET_LIMITER  SET_OS_LIMITER
 #elif LIMITER == UMIST_LIM
  #define SET_LIMITER  SET_UM_LIMITER
 #elif LIMITER == VANLEER_LIM
  #define SET_LIMITER  SET_VL_LIMITER
 #elif LIMITER == MC_LIM
  #define SET_LIMITER  SET_MC_LIMITER
 #endif
#endif

