#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  COMPONENTS                     2
#define  GEOMETRY                       CYLINDRICAL
#define  BODY_FORCE                     VECTOR
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  HANCOCK
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            1

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  ALPHA                          0

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        VANLEER_LIM

/* [End] user-defined constants (do not change this line) */
