#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  COMPONENTS                     2
#define  GEOMETRY                       CYLINDRICAL
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            3

/* -- physics dependent declarations -- */

#define  EOS                            PVTE_LAW
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  ETA                            0
#define  MACH                           1
#define  TJET                           2

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        OSPRE_LIM
#define  PV_TEMPERATURE_TABLE           NO
#define  TV_ENERGY_TABLE                NO
#define  UNIT_DENSITY                   (1.e3*CONST_amu)
#define  UNIT_LENGTH                    2.5e15
#define  UNIT_VELOCITY                  1.e5

/* [End] user-defined constants (do not change this line) */
