#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  COMPONENTS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  CHARACTERISTIC_TRACING
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            2

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

#define  EPSILON                        0
#define  BPERP_AMPL                     1

/* [Beg] user-defined constants (do not change this line) */

#define  LINEAR_SETUP                   TRUE
#define  PARTICLES_TYPE                 COSMIC_RAYS
#define  PARTICLES_DEPOSIT              REAL
#define  PARTICLES_CR_C                 1.e3
#define  PARTICLES_CR_E_MC              (1.e-6*2.0*CONST_PI)
#define  PARTICLES_CR_E_MC_GAS          1.e6
#define  PARTICLES_CR_RHO               (2.e6*g_inputParam[EPSILON])
#define  PARTICLES_PREDICTOR            1
#define  PRIMITIVE_HANCOCK              TRUE
#define  VTK_TIME_INFO                  TRUE
#define  SHOW_TIME_STEPS                FALSE

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  INTERNAL_BOUNDARY         YES
#define  SHOCK_FLATTENING          NO
#define  CHAR_LIMITING             NO
#define  LIMITER                   MC_LIM
#define  CT_EMF_AVERAGE            UCT_CONTACT
#define  CT_EN_CORRECTION          NO
#define  ASSIGN_VECTOR_POTENTIAL   NO
#define  UPDATE_VECTOR_POTENTIAL   NO
