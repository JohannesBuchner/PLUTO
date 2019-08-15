#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  COMPONENTS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 WENO3
#define  TIME_STEPPING                  RK3
#define  DIMENSIONAL_SPLITTING          YES
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            16

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  DN_PP                          0
#define  VX_PP                          1
#define  VY_PP                          2
#define  PR_PP                          3
#define  DN_MP                          4
#define  VX_MP                          5
#define  VY_MP                          6
#define  PR_MP                          7
#define  DN_MM                          8
#define  VX_MM                          9
#define  VY_MM                          10
#define  PR_MM                          11
#define  DN_PM                          12
#define  VX_PM                          13
#define  VY_PM                          14
#define  PR_PM                          15

/* [Beg] user-defined constants (do not change this line) */

#define  CHAR_LIMITING                  YES
#define  LIMITER                        VANLEER_LIM

/* [End] user-defined constants (do not change this line) */
