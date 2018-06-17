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

#define  INTERNAL_BOUNDARY              YES
#define  LIMITER                        MC_LIM
#define  CT_EMF_AVERAGE                 ARITHMETIC
#define  LINEAR_SETUP                   TRUE
#define  PARTICLES_TYPE                 COSMIC_RAYS
#define  PARTICLES_DEPOSIT              INTEGER
#define  PARTICLES_CR_C                 1.e3
#define  PARTICLES_CR_E_MC              (1.e-6*2.0*CONST_PI)
#define  PARTICLES_CR_E_MC_GAS          1.e6
#define  VTK_TIME_INFO                  TRUE
#define  SHOW_TIME_STEPS                FALSE

/* [End] user-defined constants (do not change this line) */
