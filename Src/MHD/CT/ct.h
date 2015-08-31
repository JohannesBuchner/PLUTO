/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Header file for Constrained-Transport (CT) module.

  Provides macros, function prototypes and structure definitions for 
  the constrained transport (CT) MHD module to control the divergence-free 
  condition.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 27, 2014
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

  
/* ********************************************************************* */
/*! The EMF structure is used to pull together all the information 
    necessary to build / use the electromotive force used to update
    the staggered components of magnetic field.
   ********************************************************************* */
typedef struct ElectroMotiveForce{

/*! \name Face-centered electric field components.
    Three-dimensional arrays storing the emf components computed
    at cell faces during the dimensional sweeps.     
*/
/**@{ */
  double ***exj; /**< Ex available at y-faces (j+1/2); */
  double ***exk; /**< Ex available at z-faces (k+1/2); */
  double ***eyi; /**< Ey available at x-faces (i+1/2); */
  double ***eyk; /**< Ey available at z-faces (k+1/2); */
  double ***ezi; /**< Ez available at x-faces (i+1/2); */
  double ***ezj; /**< Ez available at y-faces (j+1/2); */
/**@} */

  signed char ***svx, ***svy, ***svz;

/*! \name Range of existence */
/**@{ */
  int  ibeg, jbeg, kbeg;
  int  iend, jend, kend;
/**@} */

/*! \name Signal velocities  */
/**@{ */
  double ***SxL;
  double ***SxR;
  double ***SyL;
  double ***SyR;
  double ***SzL;
  double ***SzR;
/**@} */

/*! \name Edge-centered fields   */
/**@{ */
  double ***ex;
  double ***ey;
  double ***ez;
/**@} */

/*! \name Staggered magnetic field and velocity slopes */
/**@{ */
  double ***dbx_dy, ***dbx_dz;  
  double ***dby_dx, ***dby_dz;
  double ***dbz_dx, ***dbz_dy;

  double ***dvx_dx, ***dvy_dx, ***dvz_dx;
  double ***dvx_dy, ***dvy_dy, ***dvz_dy;
  double ***dvx_dz, ***dvy_dz, ***dvz_dz;
/**@} */

} EMF;

/* ***********************************************************
    \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES
    Function prototyping
   *********************************************************** */
   
void CT_EMF_ArithmeticAverage (const EMF *, const double);
void CT_EMF_IntegrateToCorner (const Data *, const EMF *, Grid *);
void CT_AverageMagneticField (double ****bf, double ***UU[], Grid *);
void CT_AverageNormalMagField (const Data *, int, Grid *);
void CT_AverageTransverseMagFiled (const Data *, int, Grid *);
void CT_Update(const Data *, Data_Arr, double, Grid *);
void CT_CheckDivB (double ***b[], Grid *);
void CT_GetStagSlopes (const Data_Arr, EMF *, Grid *);
void CT_StoreEMF (const State_1D *, int, int, Grid *);

EMF *CT_GetEMF (const Data *, Grid *);
void CT_AddResistiveEMF (const Data *d, Grid *grid);
void FillMagneticField (const Data *, int, Grid *); 

void CT_StoreVelSlopes (EMF *, const State_1D *, int, int);
void CT_EMF_HLL_Solver (const Data *, const EMF *, Grid *);
void CT_EMF_CMUSCL_Average (const Data *, const EMF *, Grid *);

void EMF_BOUNDARY (EMF *, Grid *);
void EMF_USERDEF_BOUNDARY (EMF *, int, int, Grid *);

/* \endcond */
 