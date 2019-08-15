#define  PHYSICS                        MHD
#define  DIMENSIONS                     1
#define  COMPONENTS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          YES
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            1

/* -- physics dependent declarations -- */

#define  EOS                            ISOTHERMAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   NO
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       EXPLICIT
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  PMODE                          0

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        MC_LIM
#define  UNIT_DENSITY                   (1.0e12*CONST_amu)
#define  UNIT_LENGTH                    (20.0)
#define  UNIT_VELOCITY                  (5487322.14)
#define  SETUP                          1

/* [End] user-defined constants (do not change this line) */
