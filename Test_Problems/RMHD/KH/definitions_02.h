#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  COMPONENTS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  HANCOCK
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        1
#define  USER_DEF_PARAMETERS            4

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT
#define  RESISTIVITY                    NO

/* -- user-defined parameters (labels) -- */

#define  SIGMA_TOR                      0
#define  SIGMA_POL                      1
#define  VEL0                           2
#define  MACH                           3

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        MC_LIM
#define  CT_EMF_AVERAGE                 ARITHMETIC
#define  ASSIGN_VECTOR_POTENTIAL        YES

/* [End] user-defined constants (do not change this line) */
