/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Collects global variables definitions.

  This file contains definitions for all global variables 
  (visible anywhere in the code) used by PLUTO. 
  Global variables names, by convention, are prefixed with a "g_" 
  unless they're used as constants throughout the code in which case 
  they keep the full-capitalized notation typical of macros.
  
  For modules, global variables are prefixed with the initial letters
  of the module name, e.g., \c sb_vy or \c glm_ch.
  
  In the following "local" means "for the local processor".
  "Interior" means inside the computational domain.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sep 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

 int SZ;
 int SZ_stagx;
 int SZ_stagy;
 int SZ_stagz;
 int SZ_float;
 int SZ_char;
 int SZ_Float_Vect;
 int SZ_rgb;
 int SZ_short;

int prank; /**< Processor rank. In serial mode it is defined to be 0. */

long int IBEG; /**< Lower grid index of the computational domain in the 
                    the X1 direction for the local processor. */
long int IEND; /**< Upper grid index of the computational domain in the 
                    the X1 direction for the local processor. */
long int JBEG; /**< Lower grid index of the computational domain in the 
                    the X2 direction for the local processor. */
long int JEND; /**< Upper grid index of the computational domain in the 
                    the X2 direction for the local processor. */
long int KBEG; /**< Lower grid index of the computational domain in the 
                    the X3 direction for the local processor. */
long int KEND; /**< Upper grid index of the computational domain in the 
                    the X3 direction for the local processor. */

long int NX1; /**< Number of interior zones in the X1 directions (boundaries
                   \e excluded) for the local processor. */
long int NX2; /**< Number of interior zones in the X2 directions (boundaries
                   \e excluded) for the local processor. */
long int NX3; /**< Number of interior zones in the X3 directions (boundaries
                   \e excluded) for the local processor. */

long int NX1_TOT; /**< Total number of zones in the X1 direction (boundaries
                       \e included) for the local processor.*/
long int NX2_TOT; /**< Total number of zones in the X2 direction (boundaries
                       \e included) for the local processor.*/
long int NX3_TOT; /**< Total number of zones in the X3 direction (boundaries
                       \e included) for the local processor.*/

long int NMAX_POINT;  /**< Maximum number of points among the three 
                           directions, boundaries \a excluded.*/

/*! \name Direction-dependent Vector Labels
    Vector indices permuted during sweeps are used to distinguish between 
    normal ("n"), tangent ("t") and bi-tangent ("b") directions.
    In vector notations, \f$ \hvec{b} = \hvec{n} \times \hvec{t} \f$, they
    form a right-handed triad.
    Values are set in the SetIndex() function before commencing 
    integration.                                                        */
/**@{ */
int VXn, VXt, VXb; 
int MXn, MXt, MXb;
int BXn, BXt, BXb;
#if DUST == YES
  int VXn_D, VXt_D, VXb_D;
  int MXn_D, MXt_D, MXb_D;
#endif
/**@} */
               
int g_i; /**<  x1 grid index when sweeping along the x2 or x3 direction. */
int g_j; /**<  x2 grid index when sweeping along the x1 or x3 direction. */
int g_k; /**<  x3 grid index when sweeping along the x1 or x2 direction. */

int g_dir; /**< Specifies the current sweep or direction of integration. 
                Its value is set usually in the time stepping functions 
                and can take the values
                - IDIR, for integration in the X1 dir; 
                - JDIR, for integration in the X2 dir; 
                - KDIR, for integration in the X3 dir; */

int      g_maxRiemannIter;  /**< Maximum number of iterations for 
                                 iterative Riemann Solver.       */
int      g_maxRootIter;  /**< Maximum number of iterations for root finder */
long int g_usedMemory;   /**< Amount of used memory in bytes. */
long int g_stepNumber;  /**< Gives the current integration step number. */
int      g_intStage;    /**< Gives the current integration stage of the time
                             stepping method (predictor = 0, 1st
                             corrector = 1, and so on). */
int      g_operatorStep; /**< Gives the current operator step. 
                              (HYPERBOLIC_STEP/PARABOLIC_STEP/SOURCE_STEP). */

double g_maxCoolingRate = 0.1;  /**< The maximum fractional variation due to 
                                     cooling from one step to the next. */
double g_minCoolingTemp = 50.0; /**< The minimum temperature (in K) below which
                                     cooling is suppressed. */
                                    
double g_smallDensity  = 1.e-12; /**< Small value for density fix. */
double g_smallPressure = 1.e-12; /**< Small value for pressure fix. */
#if EOS == IDEAL 
 double g_gamma = 5./3.;
#elif EOS == ISOTHERMAL
 double g_isoSoundSpeed = 1.0; /* g_isoSoundSpeed */
#endif

double g_time;     /**< The current integration time. */
double g_dt;       /**< The current integration time step. */
double g_maxMach;  /**< The maximum Mach number computed during integration. */
#if ROTATING_FRAME
 double g_OmegaZ;  /**< The angular rotation frequency when rotation is
                        included. */
#endif

double g_domBeg[3];  /**< Lower limits of the computational domain. */
double g_domEnd[3];  /**< Upper limits of the computational domain. */

/* g_inputParam is an array containing the user-defined
   parameters     */

double g_inputParam[32]; /**< Array containing the user-defined parameters. 
                              The index names of this array are defined in
                              definitions.h through the python interface. */
#ifdef CH_SPACEDIM
 double glm_ch_max, glm_ch_max_loc, g_coeff_dl_min; 
        /**< Variables used by Chombo to compute glm_ch. They must be 
             always declared since Chombo does not define GLM_MHD */

 double g_level_dx; /**< Global variable to pass m_dx ivalue 
                       (C++ Chombo library value) to C Pluto functions */
 double g_x2stretch, g_x3stretch; /**< Grid spacing stretch factor */
 #ifdef GLM_MHD
  int glm_is_defined = 1;  /**< Integer used to tell Chombo to compute glm_ch
                                (glm_is_defined = 1) or not (glm_is_not_defined = 0) */
 #else
  int glm_is_defined = 0;
 #endif
 #if GEOMETRY == CARTESIAN
  double g_stretch_fact;
 #endif
#endif

