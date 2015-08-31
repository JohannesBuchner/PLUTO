#define    PHYSICS   RMHD
#define    DIMENSIONS   3
#define    COMPONENTS   3
#define    GEOMETRY   CARTESIAN
#define    INCLUDE_BODY_FORCE   NO
#define    INCLUDE_COOLING   NO
#define    INCLUDE_PARTICLES   NO
#define    INTERPOLATION   LINEAR
#define    TIME_STEPPING   HANCOCK
#define    DIMENSIONAL_SPLITTING   YES
#define    NTRACER   1
#define    USER_DEF_PARAMETERS   6

/* -- physics dependent declarations -- */

#define    EOS   TAUB
#define    MHD_FORMULATION   DIV_CLEANING

/* -- pointers to user-def parameters -- */

#define  MACH   0
#define  LORENTZ   1
#define  ETA   2
#define  SIGMA_TOR   3
#define  SIGMA_POL   4
#define  EPS   5

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING     NO
#define  WARNING_MESSAGES      NO
#define  PRINT_TO_FILE         NO
#define  SHOCK_FLATTENING      MULTID
#define  ARTIFICIAL_VISCOSITY  NO
#define  CHAR_LIMITING         NO
#define  LIMITER               DEFAULT
#define  SAVE_VEC_POT          NO
#define  PRIMITIVE_HANCOCK     NO
