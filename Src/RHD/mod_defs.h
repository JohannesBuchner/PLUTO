/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Set labels, indexes and prototypes for the RHD module.

  Contains variable names and prototypes for the RHD module

  \author A. Mignone (mignone@ph.unito.it)
  \date   April, 2, 2015
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

#define NFLX (2 + COMPONENTS)

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

#endif

#if GEOMETRY == POLAR 

 #define iVR    VX1
 #define iVPHI  VX2
 #define iVZ    VX3

 #define iMR    MX1
 #define iMPHI  MX2
 #define iMZ    MX3

#endif

#if GEOMETRY == SPHERICAL 

 #define iVR    VX1
 #define iVTH   VX2
 #define iVPHI  VX3

 #define iMR    MX1
 #define iMTH   MX2
 #define iMPHI  MX3

#endif

/* ********************************************************************* */
/*! The Map_param structure is used to pass input/output arguments 
    during the conversion from conservative to primitive variables 
    operated by the ConsToPrim() function.
    The output parameter, rho, W, lor and p, must be set at the end
    of every root-finder routine (EnergySolve(), EntropySolve() and
    PressureFix()).
    Additionally, some of the input parameters must be re-computed in
    EntropySolve() and PressureFix().
   ********************************************************************* */
typedef struct MAP_PARAM{
 double D;       /**< Lab density       (input). */
 double sigma_c; /**< Conserved entropy (input). */
 double E;       /**< Total energy      (input). */
 double m2;      /**< Square of total momentum (input). */

 double rho;     /**< proper density     (output)  */
 double W;       /**< D*h*lor            (output). */
 double lor;     /**< Lorentz factor     (output). */
 double prs;     /**< Thermal pressure   (output). */
} Map_param;

/* ---- Function prototyping ---- */

int  ConsToPrim   (double **, double **, int, int, unsigned char *);
void PrimEigenvectors(double *, double, double, double *, double **, double **);

int EnergySolve  (Map_param *);
int EntropySolve (Map_param *);
int PressureFix  (Map_param *);

void Flux      (double **, double **, double *, double **, double *, int, int);
void HLL_Speed (double **, double **, double *, double *,
                double *, double *, int, int);
void MaxSignalSpeed (double **, double *, double *, double *, int, int);
void PrimToCons   (double **, double **, int, int);
void PrimRHS    (double *, double *, double, double, double *);
void PrimSource (const State_1D *, int, int, 
                 double *, double *, double **, Grid *);
void VelocityLimiter(double *, double *, double *);
void ConvertTo4vel (double **, int, int);
void ConvertTo3vel (double **, int, int);


Riemann_Solver TwoShock_Solver, LF_Solver, HLL_Solver, HLLC_Solver;

