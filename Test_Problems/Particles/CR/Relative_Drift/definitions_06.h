#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  COMPONENTS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            3

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  RHO_GAS                        0
#define  VPX1                           1
#define  VPX2                           2

/* [Beg] user-defined constants (do not change this line) */

#define  WARNING_MESSAGES               NO
#define  CT_EMF_AVERAGE                 UCT_CONTACT
#define  PARTICLES_CR_C                 1.e40
#define  PARTICLES_CR_FEEDBACK          YES
#define  PARTICLES_CR_PREDICTOR         2
#define  PARTICLES_CR_E_MC              1.0
#define  PARTICLES_CR_NSUB              -4
#define  PARTICLES_SHAPE                3
#define  PARTICLES_TYPE                 COSMIC_RAYS
#define  SHOW_TIME_STEPS                FALSE
#define  PRIMITIVE_HANCOCK              TRUE

/* [End] user-defined constants (do not change this line) */
