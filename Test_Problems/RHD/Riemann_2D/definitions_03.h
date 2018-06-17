#define  PHYSICS                        RHD
#define  DIMENSIONS                     2
#define  COMPONENTS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  HANCOCK
#define  DIMENSIONAL_SPLITTING          YES
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            1

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO

/* -- user-defined parameters (labels) -- */

#define  SCRH                           0

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        MINMOD_LIM
#define  RECONSTRUCT_4VEL               YES

/* [End] user-defined constants (do not change this line) */
