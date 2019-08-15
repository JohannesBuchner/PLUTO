/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Module header file for relativistic MHD (RMHD).

  Set label, indexes and basic prototyping for the relativistic 
  MHD module.

  \authors A. Mignone (mignone@ph.unito.it)\n
           C. Zanni   (zanni@oato.inaf.it)\n
           G. Mattia
  \date    Nov 27, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */

/* *********************************************************
    Set flow variable indices.
    Extra vector components, when not needed, point to the
    last element (255) of the array stored by startup.c.  
   ********************************************************* */

#if RESISTIVITY != NO
  #define RMHD_NVEC 3  /* Number of vectors: velocity, B, E */
#else
  #define RMHD_NVEC 2  /* Number of vectors: velocity, B    */
#endif


/* Use the following switch to enable/disable chekcing of the
   inversion process used with resistive MHD (but in ideal case)  */

#ifndef RMHD_RESISTIVE_CHECK
  #define RMHD_RESISTIVE_CHECK  NO   
#endif


#define  RHO 0
#define  MX1 1
#define  MX2 (COMPONENTS >= 2 ? 2: 255)
#define  MX3 (COMPONENTS == 3 ? 3: 255)
#define  BX1 (COMPONENTS + 1)
#define  BX2 (COMPONENTS >= 2 ? (BX1+1): 255)
#define  BX3 (COMPONENTS == 3 ? (BX1+2): 255)

#if HAVE_ENERGY
  #define ENG  (2*COMPONENTS + 1)
  #define PRS  ENG
#endif

#if RESISTIVITY 
  #if COMPONENTS != 3
    #error "Must use 3 COMPONENTS with resistive RMHD"
  #endif
  #define  EX1 (ENG + 1)
  #define  EX2 (ENG + 2)
  #define  EX3 (ENG + 3)
  #define  CRG (ENG + 4)
#endif

#if DIVB_CONTROL == DIV_CLEANING
  #if RESISTIVITY
    #define PSI_GLM  (CRG + 1)
    #define PHI_GLM  (CRG + 2)
    #define DIV_COMP 2
  #else
    #define PSI_GLM  (ENG + 1)
    #define DIV_COMP 1
  #endif
#else
  #define DIV_COMP 0
#endif

#define VX1   MX1
#define VX2   MX2
#define VX3   MX3

#define NFLX (1 + RMHD_NVEC*COMPONENTS + DIV_COMP + HAVE_ENERGY + (RESISTIVITY != NO))

/* *********************************************************
    Label the different waves in increasing order 
    following the number of vector components.

    IMPORTANT: the KPSI_GLMM & KPSI_GLMP modes are 
               present only in the MHD-GLM formulation.
               We keep them at the END of the enumeration
               so we can skip them in unnecessary loops.
               Please do NOT change them !
   ********************************************************* */

enum KWAVES {
 KFASTM, KFASTP, KENTRP

 #if DIVB_CONTROL != DIV_CLEANING
  , KDIVB
 #endif

 #if COMPONENTS >= 2
  , KSLOWM, KSLOWP
  #if COMPONENTS == 3
   , KALFVM, KALFVP
  #endif
 #endif

 #if DIVB_CONTROL == DIV_CLEANING  
  , KPSI_GLMM, KPSI_GLMP 
 #endif
};

                                     
/* ********************************************************************* */
/*! The Map_param structure is used to pass input/output arguments 
    during the conversion from conservative to primitive variables 
    operated by the ConsToPrim() function in the relativistic modules
    (RHD and RMHD).
    The output parameter, rho, W, lor and p, must be set at the end
    of every root-finder routine (EnergySolve(), EntropySolve() and
    PressureFix()).
    Additionally, some of the input parameters must be re-computed in
    EntropySolve() and PressureFix().
   ********************************************************************* */
#if RESISTIVITY == NO
typedef struct MAP_PARAM{
 double D;       /**< Lab density       (input). */
 double sigma_c; /**< Conserved entropy (input). */
 double E;       /**< Total energy      (input). */
 double m2;      /**< Square of total momentum (input). */
 double S;       /**< m<dot>B                  (input). */
 double S2;      /**< Square of S              (input). */
 double B2;      /**< Square of magnetic field (input). */

 double rho;     /**< proper density     (output)  */
 double W;       /**< D*h*lor            (output). */
 double lor;     /**< Lorentz factor     (output). */
 double prs;     /**< Thermal pressure   (output). */
} Map_param;
#else
typedef struct MAP_PARAM{
 double D;       /**< Lab density       (input). */
 double sigma_c; /**< Conserved entropy (input). */
 double E;       /**< Total energy      (input). */
 double m2;      /**< Square of total momentum (input). */
 double S;       /**< m<dot>B                  (input). */
 double S2;      /**< Square of S              (input). */
 double B2;      /**< Square of magnetic field (input). */

 double rho;     /**< proper density     (output)  */
 double W;       /**< D*h*lor            (output). */
 double lor;     /**< Lorentz factor     (output). */
 double prs;     /**< Thermal pressure   (output). */
} Map_param;
#endif

/* ******************************************************
     Vector potential: these labels are and MUST only
     be used in the STARTUP / INIT  functions;
     they're convenient in obtaining a discretization 
     that preserve divB since the beginning.   
   ****************************************************** */

#define   AX1  (NVAR + 1)
#define   AX2  (NVAR + 2)
#define   AX3  (NVAR + 3)

#define AX  AX1  /* backward compatibility */
#define AY  AX2
#define AZ  AX3

/* *************************************************
     Now define more convenient and user-friendly 
     pointer labels for geometry setting      
   ************************************************* */

#if GEOMETRY == CYLINDRICAL
                               
 #define iVR    VX1
 #define iVZ    VX2
 #define iVPHI  VX3
                               
 #define iMR    MX1
 #define iMZ    MX2
 #define iMPHI  MX3
                               
 #define iBR    BX1
 #define iBZ    BX2
 #define iBPHI  BX3
                               
#endif

#if GEOMETRY == POLAR
                               
 #define iVR    VX1
 #define iVPHI  VX2
 #define iVZ    VX3
                               
 #define iMR    MX1
 #define iMPHI  MX2
 #define iMZ    MX3
                               
 #define iBR    BX1
 #define iBPHI  BX2
 #define iBZ    BX3
                               
#endif

#if GEOMETRY == SPHERICAL
                    
 #define iVR     VX1
 #define iVTH    VX2
 #define iVPHI   VX3
                             
 #define iMR    MX1
 #define iMTH   MX2
 #define iMPHI  MX3
                               
 #define iBR    BX1
 #define iBTH   BX2
 #define iBPHI  BX3
 
#endif

/* ******************************************************************
    Module-specific symbolic constants (switches)
   ****************************************************************** */

#ifndef RMHD_FAST_EIGENVALUES
  #define RMHD_FAST_EIGENVALUES  NO /**< If set to YES, use approximate (and
                                         faster) expressions when computing the
                                         fast magnetosonic speed, see Sect. 3.3
                                         of Del Zanna,  A&A (2006), 473.
                                         Solutions of quartic equation is avoided
                                         and replace with upper bounds provided by
                                         quadratic equation. */
#endif  

#ifndef RMHD_REDUCED_ENERGY
  #define RMHD_REDUCED_ENERGY    YES  /**< By turning RMHD_REDUCED_ENERGY to YES, 
                                            we let PLUTO evolve the total energy 
                                            minus the mass density contribution. */
#endif                                            

/* ---- Function prototyping ----  */

int  ConsToPrim   (double **, double **, int, int, unsigned char *);
void ConvertTo4vel (double **, int, int);
void ConvertTo3vel (double **, int, int);
void PRIM_EIGENVECTORS (double *, double, double, double *, double **, double **);
int  ApproximateFastWaves  (double *, double, double, double *);
int  EntropySolve (Map_param *);
int  EnergySolve  (Map_param *);
#if RESISTIVITY == NO
int  PressureFix  (Map_param *);
#else
int  PressureFix  (double *, double *, double *);
#endif
 
void Flux      (const State *, int, int);
void HLL_Speed (const State *, const State *, double *, double *, int, int);
int  MaxSignalSpeed (const State *, double *, double *, int, int);

void PrimToCons   (double **, double **, int, int);
void VelocityLimiter (double *, double *, double *);

int  Magnetosonic (double *vp, double cs2, double h, double *lambda);

Riemann_Solver LF_Solver, HLL_Solver, HLLC_Solver, HLLD_Solver, HLL_Linde_Solver; 
Riemann_Solver GMUSTA1_Solver;
Riemann_Solver Blended_HLLX_Solver;

#if RESISTIVITY != NO
 #include "Resistivity/resistivity.h"
#endif

#if DIVB_CONTROL == EIGHT_WAVES
 void POWELL_DIVB_SOURCE(const Sweep *, int, int, Grid *);
 void HLL_DIVB_SOURCE (const Sweep *, double **, int, int, Grid *);
#elif DIVB_CONTROL == DIV_CLEANING
 #include "MHD/GLM/glm.h"
#elif DIVB_CONTROL == CONSTRAINED_TRANSPORT
  #include "MHD/CT/ct.h"
#endif

