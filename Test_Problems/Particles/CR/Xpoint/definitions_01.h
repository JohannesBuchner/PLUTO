#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  COMPONENTS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  HANCOCK
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            4

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

#define  ALPHA                          0
#define  BETA                           1
#define  BMAG_Z                         2
#define  EMAG_Z                         3

/* [Beg] user-defined constants (do not change this line) */

#define  WARNING_MESSAGES               NO
#define  CT_EMF_AVERAGE                 UCT_CONTACT
#define  ASSIGN_VECTOR_POTENTIAL        YES
#define  UPDATE_VECTOR_POTENTIAL        YES
#define  UNIT_DENSITY                   (1.0)
#define  UNIT_LENGTH                    (1.0)
#define  UNIT_VELOCITY                  (CONST_c/100.0)
#define  PARTICLES_CR_C                 (CONST_c/UNIT_VELOCITY)
#define  PARTICLES_CR_E_MC              1.0
#define  PARTICLES_CR_FEEDBACK          NO
#define  PARTICLES_CR_LARMOR_EPS        0.5
#define  PARTICLES_CR_NSUB              -1
#define  PARTICLES_CR_PREDICTOR         NO
#define  PARTICLES_TEST                 YES
#define  PARTICLES_TYPE                 COSMIC_RAYS

/* [End] user-defined constants (do not change this line) */
