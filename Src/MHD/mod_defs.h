/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Set labels, indexes and prototypes for the MHD module.

  Contains basic macro definitions, structure definitions and global
  variable declarations used by the MHD module.

  \author A. Mignone (mignone@ph.unito.it)
  \date   June 19, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */

/* *********************************************************
    Set flow variable indices.
    Extra vector components, when not needed, point to the
    last element (255) of the array stored by startup.c.  
   ********************************************************* */

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
#if DIVB_CONTROL == DIV_CLEANING
  #define PSI_GLM  (2*COMPONENTS + 1 + HAVE_ENERGY)
#endif

#define VX1   MX1
#define VX2   MX2
#define VX3   MX3

#define NFLX (1 + 2*COMPONENTS + HAVE_ENERGY + (DIVB_CONTROL == DIV_CLEANING))

/* ********************************************************************* */
/*! Label the different waves in increasing order 
    following the number of vector components.

    \b IMPORTANT: the KPSI_GLMM & KPSI_GLMP modes are 
                  present only in the MHD-GLM formulation.
                  We keep them at the END of the enumeration
                  so we can skip them in unnecessary loops.
                  Please do NOT change them !
   ********************************************************************* */

enum KWAVES {
 KFASTM, KFASTP
 #if HAVE_ENERGY
  , KENTRP
 #endif

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

/*! \name Vector Potential Labels 
    These may only be used in the STARTUP / INIT  functions.
    They're convenient in obtaining a discretization that preserve 
    the divergence-free condition (for staggered field) or if you simply
    wish to initialize the magnetic field from the vector potential.     */
/**@{ */
#define   AX1  (NVAR + 1)
#define   AX2  (NVAR + 2)
#define   AX3  (NVAR + 3)
/**@} */

#define AX  AX1  
#define AY  AX2
#define AZ  AX3

/* *************************************************
     Now define more convenient and user-friendly 
     pointer labels for geometry setting      
   ************************************************* */

#if GEOMETRY == CARTESIAN
  #define VX    VX1
  #define VY    VX2
  #define VZ    VX3

  #define MX    MX1
  #define MY    MX2
  #define MZ    MX3

  #define BX    BX1
  #define BY    BX2
  #define BZ    BX3
#endif

#if GEOMETRY == CYLINDRICAL
  #if COMPONENTS >= 1
    #define iVR    VX1
    #define iMR    MX1
    #define iBR    BX1
  #endif

  #if COMPONENTS >= 2
    #define iVZ    VX2
    #define iMZ    MX2
    #define iBZ    BX2
  #endif

  #if COMPONENTS >= 3 
    #define iVPHI  VX3
    #define iMPHI  MX3
    #define iBPHI  BX3
  #endif  
#endif

#if GEOMETRY == POLAR
  #if COMPONENTS >= 1
    #define iVR    VX1
    #define iMR    MX1
    #define iBR    BX1
  #endif

  #if COMPONENTS >= 2 
    #define iVPHI  VX2
    #define iMPHI  MX2
    #define iBPHI  BX2
  #endif  

  #if COMPONENTS == 3
    #define iVZ    VX3
    #define iMZ    MX3
    #define iBZ    BX3
  #endif
#endif

#if GEOMETRY == SPHERICAL
  #if COMPONENTS >= 1
    #define iVR    VX1
    #define iMR    MX1
    #define iBR    BX1
  #endif
  
  #if COMPONENTS >= 2
    #define iVTH   VX2
    #define iMTH   MX2
    #define iBTH   BX2
  #endif
 
  #if COMPONENTS == 3    
    #define iVPHI  VX3
    #define iMPHI  MX3
    #define iBPHI  BX3
  #endif
#endif

/* ***********************************************************
    \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES 
    Function prototyping
   *********************************************************** */

void BackgroundField (double x1, double x2, double x3, double *B0);

void ConsEigenvectors (double *, double *, double,
                       double **, double **, double *);
int  ConsToPrim   (double **, double **, int , int, unsigned char *);
void Eigenvalues (double **, double *, double **, int, int);

void Flux (const State *, int, int);
#if BACKGROUND_FIELD == YES
void GetBackgroundField (const State *, int, int, int, Grid *);
#endif
void GetCurrent (const Data *, Grid *);
void HLL_Speed (const State *, const State *, double *, double *, int, int);

void MaxSignalSpeed (const State *, double *, double *, int, int);

void PrimEigenvectors(const State *, int, int);
void PrimRHS     (double *, double *, double, double, double *);
void PrimSource  (const State *, double **, int, int, Grid *);
void PrimToCons  (double **, double **, int, int);

#if DIVB_CONTROL == EIGHT_WAVES
 void Roe_DivBSource (const Sweep *, int, int, Grid *);
 void HLL_DivBSource (const Sweep *, double **, int, int, Grid *);
#elif DIVB_CONTROL == DIV_CLEANING

 #include "MHD/GLM/glm.h"

#elif DIVB_CONTROL == CONSTRAINED_TRANSPORT

 #include "MHD/CT/ct.h"

#endif

Riemann_Solver HLL_Solver, HLLC_Solver, HLLD_Solver;
Riemann_Solver LF_Solver, Roe_Solver;
Riemann_Solver HLL_Linde_Solver;

#if AMBIPOLAR_DIFFUSION != NO
 #include "Ambipolar_Diffusion/ad.h"
#endif
#if RESISTIVITY != NO
 #include "Resistivity/res.h"
#endif

#ifdef SHEARINGBOX
 #include "MHD/ShearingBox/shearingbox.h"
#endif
/* \endcond */
