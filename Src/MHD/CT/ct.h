/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Header file for Constrained-Transport (CT) module.

  Provides macros, function prototypes and structure definitions for 
  the constrained transport (CT) MHD module to control the divergence-free 
  condition.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Feb 21, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#define STAGGERED_MHD

/* ---- set labels for CT_EMF_AVERAGE ---- */

#define ARITHMETIC     1
#define UCT0           2
#define UCT_CONTACT    3
#define UCT_HLL        4
#define RIEMANN_2D     5

/* ---- staggered component labels ---- */

#define BX1s  0
#define BX2s  1
#define BX3s  2

/* ---- backward compatibility ---- */

#define BXs BX1s
#define BYs BX2s
#define BZs BX3s

#define FACE_EMF   11
#define EDGE_EMF   12

/* *********************************************************************
   Default values
   ********************************************************************* */

#ifndef CT_EMF_AVERAGE
 #define  CT_EMF_AVERAGE            UCT_HLL
#endif

#ifndef  CT_EN_CORRECTION
 #define  CT_EN_CORRECTION          NO
#endif

/* *********************************************************************
    The following switch can be used to select where, in the CT
    formalism, current should be computed.
    When GET_CURRENT_AT_EDGES == YES, current will be computed
    at edges (like the emf). Otherwise at cell-interfaces as it
    is usually done for cell-centered methods.
    For Hall MHD, we see better stability when it is set to NO.
   ********************************************************************* */

#ifndef GET_CURRENT_AT_EDGES
  #if HALL_MHD
  #define GET_CURRENT_AT_EDGES  NO
  #else
  #define GET_CURRENT_AT_EDGES  YES
  #endif
#endif

/*       Now define more convenient and user-friendly 
         pointer labels for geometry setting      */

#if GEOMETRY == CYLINDRICAL 
 #define iBRs    BX1s
 #define iBZs    BX2s
 #define iBPHIs  BX3s
#endif

#if GEOMETRY == POLAR 
 #define iBRs    BX1s
 #define iBPHIs  BX2s
 #define iBZs    BX3s
#endif

#if GEOMETRY == SPHERICAL 
 #define iBRs    BX1s
 #define iBTHs   BX2s
 #define iBPHIs  BX3s
#endif

/* ***********************************************************
    \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES
    Function prototyping
   *********************************************************** */

void CT_Allocate (EMF *);   
void CT_EMF_ArithmeticAverage (const EMF *, const double);
void CT_EMF_IntegrateToCorner (const Data *, const EMF *, Grid *);
void CT_AverageMagneticField (double ****bf, double ***UU[], Grid *);
void CT_AverageNormalMagField (const Data *, int, Grid *);
void CT_AverageTransverseMagField (const Data *, int, Grid *);
void CT_CheckDivB (double ***b[], Grid *);
void CT_ComputeCenterEMF(const Data *);
void CT_ComputeEMF (const Data *, Grid *);
void CT_EMF_HLL_Solver (const Data *, const EMF *, Grid *);
void CT_GetStagSlopes (const Data_Arr, EMF *, Grid *);
void CT_GetEMF (const Data *, Grid *);
void CT_ResistiveEMF (const Data *, int, Grid *);

void FillMagneticField (const Data *, int, Grid *); 

void CT_StoreUpwindEMF    (const Sweep *, EMF *, int, int, Grid *);
void CT_StoreVelSlopes (EMF *, const Sweep *, int, int);
void CT_Update(const Data *, Data_Arr, double, Grid *);

/* \endcond */
 
