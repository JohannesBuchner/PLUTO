#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  COMPONENTS                     2
#define  GEOMETRY                       CYLINDRICAL
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  CHARACTERISTIC_TRACING
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            4

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  CS_WIND                        0
#define  RHO_AMB                        1
#define  CS_AMB                         2
#define  V_CSM                          3

/* [Beg] user-defined constants (do not change this line) */

#define  INITIAL_SMOOTHING              YES
#define  INTERNAL_BOUNDARY              YES
#define  CHAR_LIMITING                  YES

/* [End] user-defined constants (do not change this line) */
