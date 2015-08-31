/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Shearing-Box module header file

  The Shearing-Box module header file contains basic macro definitions, 
  function prototypes and declaration of global variables used by the 
  sheraring-box module.
  The variable ::sb_q and ::sb_Omega are the most important ones and 
  \b must be defined and initialized in your init.c in order to 
  configure your shearing-box problem.

  Optionally, the order of interpolation (default is 2) at physical 
  boundaries may be changed using the ::SB_ORDER macro.

  The additional macros ::SB_SYMMETRIZE_HYDRO, ::SB_SYMMETRIZE_EY and
  ::SB_SYMMETRIZE_EZ may be set to YES/NO to enable/disable enforcement
  of conservation at the radial (x) boundaries.

  \authors A. Mignone (mignone@ph.unito.it)\n
           G. Muscianisi (g.muscianisi@cineca.it)
  \date   Aug 26, 2015
  \todo   Check if sb_vy and sb_Ly are really needed as global variables.
*/
/* ///////////////////////////////////////////////////////////////////// */

/*! Sets the order of interpolation at physical boundaries (1, 2 or 3).*/
#if RECONSTRUCTION == LINEAR || RECONSTRUCTION == WENO3 || RECONSTRUCTION == LimO3
 #define SB_ORDER              2    
#else
 #define SB_ORDER              3    
#endif                                 

#ifndef FARGO
 #ifndef SB_SYMMETRIZE_HYDRO
  #define SB_SYMMETRIZE_HYDRO  YES  /**< Symmetrize the hydrodynamical fluxes 
          at the left and right x-boundaries in order to enforce conservation of
          hydrodynamic variables like density, momentum and energy
          (no magnetic field). Default is YES.*/
 #endif
 
 #ifndef SB_SYMMETRIZE_EY
  #define SB_SYMMETRIZE_EY    (YES  && (DIMENSIONS == 3)) /**< Symmetrize the
          y-component of the electric field at the left and right x-boundaries
          to enforce conservation of magnetic field (only in 3D). */
 #endif

 #ifndef SB_SYMMETRIZE_EZ
  #define SB_SYMMETRIZE_EZ     YES  /**< Symmetrize the z-component of electric
          field at the left and right x-boundaries to enforce conservation of
          magnetic field. */
 #endif

 #define SB_FORCE_EMF_PERIODS NO  /**< Force periodicity at y- and z-
                                       boundaries. */
#else
 #define SB_SYMMETRIZE_HYDRO  NO  
 #define SB_SYMMETRIZE_EY    (NO  && (DIMENSIONS == 3)) 
 #define SB_SYMMETRIZE_EZ     NO 
 #define SB_FORCE_EMF_PERIODS NO 
#endif

/* ----------------------------
    Global variables
   ---------------------------- */

#ifndef SB_Q
  #define SB_Q  1.5 /**< The shear parameter, \f$\DS q = -\HALF\frac{d\log
                         \Omega^2}{d\log R} \f$. */
#endif
#ifndef SB_OMEGA
  #define SB_OMEGA   1.0 /**< Disk local orbital frequency \f$ \Omega_0 = 
                              \Omega(R_0)\f$.  */
#endif

//extern double sb_q;  /**< The shear parameter, \f$\DS q = -\HALF\frac{d\log
//              \Omega^2}{d\log R} \f$. The explicit numerical value and the
//              variable definition should be set inside your Init() function. */

//extern double sb_Omega; /**< Disk local orbital frequency \f$ \Omega_0 = 
//              \Omega(R_0)\f$. The explicit numerical value and the variable
 //             definition should be set inside your Init() function.  */

extern double sb_vy;  /**< Velocity offset (>0), in SB_Boundary(). */ 

#define SB_A (-0.5*SB_OMEGA*SB_Q)  /**< Short-hand definition for the Oort
                                        constant \f$ A = -q\Omega_0/2 \f$. */

/* ***********************************************************
    \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES 
    Function prototyping
   *********************************************************** */

void SB_Boundary (const Data *, int, Grid *) ;
void SB_ShearingInterp (double *, double *, double, int, Grid *);
void SB_CorrectFluxes  (Data_Arr, double, double, Grid *);
#ifdef STAGGERED_MHD
void SB_CorrectEMF (EMF *, Data_Arr, Grid *);
#endif
int  SB_JSHIFT (int);
void SB_SaveFluxes (State_1D *, Grid *);
void SB_SetBoundaryVar(double ***, RBox *, int, double, Grid *);
void SB_ShiftBoundaryVar(double ***, RBox *, int, double, Grid *);
void SB_FillBoundaryGhost(double ***, RBox *, int, int, Grid *);

#ifdef PARALLEL
 void ExchangeX (real *, real *, int, Grid *);
#endif
 
/* \endcond */
