#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  COMPONENTS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            2

/* -- physics dependent declarations -- */

#define  EOS                            PVTE_LAW
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  RHO_IN                         0
#define  PRS_IN                         1

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        MC_LIM
#define  ADD_TURBULENCE                 NO
#define  UNIT_LENGTH                    2.5e15
#define  UNIT_DENSITY                   (1.e5*CONST_amu)
#define  UNIT_VELOCITY                  2.5e5

/* [End] user-defined constants (do not change this line) */
