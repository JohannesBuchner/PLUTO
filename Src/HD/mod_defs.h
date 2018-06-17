/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Set labels, indexes and prototypes for the HD module.

  Contains variable names and prototypes for the HD module

  \author A. Mignone (mignone@ph.unito.it)
  \date   Oct 11, 2016
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
#if HAVE_ENERGY
  #define ENG  (COMPONENTS + 1)
  #define PRS  ENG
#endif

#define VX1   MX1
#define VX2   MX2
#define VX3   MX3

#define NFLX (1 + COMPONENTS + HAVE_ENERGY)

/* *************************************************
     Now define more convenient and user-friendly 
     pointer labels for geometry setting      
   ************************************************* */

#if GEOMETRY == CYLINDRICAL 
  #if COMPONENTS >= 1
    #define iVR    VX1
    #define iMR    MX1
  #endif

  #if COMPONENTS >= 2
    #define iVZ    VX2
    #define iMZ    MX2
  #endif

  #if COMPONENTS == 3
    #define iVPHI  VX3
    #define iMPHI  MX3
  #endif
#endif

#if GEOMETRY == POLAR 
  #if COMPONENTS >= 1
    #define iVR    VX1
    #define iMR    MX1
  #endif

  #if COMPONENTS >= 2
    #define iVPHI  VX2
    #define iMPHI  MX2
  #endif

  #if COMPONENTS == 3
    #define iVZ    VX3
    #define iMZ    MX3
  #endif
#endif

#if GEOMETRY == SPHERICAL 
  #if COMPONENTS >= 1
    #define iVR    VX1
    #define iMR    MX1
  #endif

  #if COMPONENTS >= 2
    #define iVTH   VX2
    #define iMTH   MX2
  #endif

  #if COMPONENTS == 3
    #define iVPHI  VX3
    #define iMPHI  MX3
  #endif
#endif

/* *************************************************
     Label the different waves in increasing order 
     following the number of vector components.
   ************************************************* */

enum KWAVES {
 KSOUNDM, KSOUNDP
 #if HAVE_ENERGY
  , KENTRP
 #endif
};

/* ***********************************************************
                   Prototyping goes here          
   *********************************************************** */

int  ConsToPrim   (double **, double **, int, int, unsigned char *);
void Eigenvalues (double **, double *, double **, int, int);
void PrimEigenvectors (const State *, int, int);
void ConsEigenvectors (double *, double *, double, 
                       double **, double **, double *);

void Flux      (const State *, int, int);
void HLL_Speed (const State *, const State *, double *, double *, int, int);
void MaxSignalSpeed (const State *, double *, double *, int, int);
void PrimToCons   (double **, double **, int, int);
void PrimRHS    (double *, double *, double, double, double *);
void PrimSource (const State *, double **, int, int, Grid *);

Riemann_Solver TwoShock_Solver, LF_Solver, Roe_Solver, HLL_Solver, 
               HLLC_Solver, RusanovDW_Solver;
Riemann_Solver AUSMp_Solver;

