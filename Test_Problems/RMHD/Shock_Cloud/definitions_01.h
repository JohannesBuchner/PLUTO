#define  PHYSICS                 RMHD
#define  DIMENSIONS              2
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           HANCOCK
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     0

/* -- physics dependent declarations -- */

#define  EOS                     TAUB
#define  ENTROPY_SWITCH          NO
#define  DIVB_CONTROL            NO

/* -- user-defined parameters (labels) -- */


/* [Beg] user-defined constants (do not change this line) */

#define  CHOMBO_REF_VAR          RHO

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          YES
#define  PRINT_TO_FILE             NO
#define  INTERNAL_BOUNDARY         NO
#define  SHOCK_FLATTENING          MULTID
#define  CHAR_LIMITING             NO
#define  LIMITER                   VANLEER_LIM
#define  ASSIGN_VECTOR_POTENTIAL   NO
#define  PRIMITIVE_HANCOCK         NO
