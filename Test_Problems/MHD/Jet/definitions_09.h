#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  COMPONENTS                     3
#define  GEOMETRY                       CYLINDRICAL
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        H2_COOL
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  HANCOCK
#define  DIMENSIONAL_SPLITTING          YES
#define  NTRACER                        1
#define  USER_DEF_PARAMETERS            6

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   NO
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  ETA                            0
#define  JET_VEL                        1
#define  SIGMA_Z                        2
#define  SIGMA_PHI                      3
#define  PERT_AMPLITUDE                 4
#define  PERT_PERIOD                    5

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        VANLEER_LIM
#define  UNIT_DENSITY                   (CONST_amu*200.0)
#define  UNIT_LENGTH                    (2.5e15)
#define  UNIT_VELOCITY                  (1.e5)
#define  WARNING_MESSAGES               NO

/* [End] user-defined constants (do not change this line) */
