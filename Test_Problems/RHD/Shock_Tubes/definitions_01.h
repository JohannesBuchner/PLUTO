#define  PHYSICS                 RHD
#define  DIMENSIONS              1
#define  COMPONENTS              1
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           HANCOCK
#define  DIMENSIONAL_SPLITTING   YES
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     9

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO

/* -- user-defined parameters (labels) -- */

#define  GAMMA_EOS               0
#define  DN_L                    1
#define  VX_L                    2
#define  VY_L                    3
#define  PR_L                    4
#define  DN_R                    5
#define  VX_R                    6
#define  VY_R                    7
#define  PR_R                    8

/* [Beg] user-defined constants (do not change this line) */


/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING   NO
#define  WARNING_MESSAGES    YES
#define  PRINT_TO_FILE       NO
#define  INTERNAL_BOUNDARY   NO
#define  SHOCK_FLATTENING    YES
#define  CHAR_LIMITING       NO
#define  LIMITER             FOURTH_ORDER_LIM
#define  PRIMITIVE_HANCOCK   YES
