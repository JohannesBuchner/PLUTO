/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Particle module header file

  Contains basic macro definitions, structure definitions and global
  variable declarations used by the Particle module.

  \author A. Mignone (mignone@ph.unito.it)\n
          B. Vaidya
  \date   June 29, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */

/* Set default particle type */

#ifndef PARTICLES_TYPE
  #define PARTICLES_TYPE  LAGRANGIAN
#endif

#define INTEGER       1  /* Used for particle deposition */
#define REAL          2

#define PARTICLES_CREATE    1
#define PARTICLES_TRANSFER  2
#define PARTICLES_RESTART   3

#ifndef PARTICLES_LP_SPECTRA
  #define PARTICLES_LP_SPECTRA   NO 
#endif

#ifndef PARTICLES_SHAPE
  #define PARTICLES_SHAPE  3  /* Weight functions: 1 = Nearestg Grid Point (NGP)
                                                   2 = Cloud-In-Cell (CIC)
                                                   3 = Triangular Shape Cloud (TSC)  */
#endif

/* *********************************************************************
    The PARTICLES_TEST directive can be used to evolve just the
    particles and not the fluid.
   ********************************************************************* */

#ifndef PARTICLES_TEST
  #define PARTICLES_TEST  NO
#endif

/* *********************************************************************
    Test/Debug directive
   ********************************************************************* */
 
#ifndef PARTICLES_DEPOSIT
  #define PARTICLES_DEPOSIT   REAL
#endif

/* *********************************************************************
    MPI Data type
    Set PARTICLES_USE_MPI_DATATYPE to YES to define an appropriate
    datatype for MPI send/recv of structure. 
    Set PARTICLES_USE_MPI_DATATYPE to NO to pass the structure as a
    sequence of BYTEs (requires no specific datatype) 
   ********************************************************************* */

#ifndef PARTICLES_USE_MPI_DATATYPE
  #if PARTICLES_TYPE == LAGRANGIAN
    #define PARTICLES_USE_MPI_DATATYPE  YES
  #else
    #define PARTICLES_USE_MPI_DATATYPE  YES
  #endif
#endif

/* *********************************************************************
    Particle type: COSMIC RAYS
   ********************************************************************* */

#if PARTICLES_TYPE == COSMIC_RAYS

  #ifndef PARTICLES_CR_C                /* Speed of light */
    #if PHYSICS == RMHD
      #define PARTICLES_CR_C       1.0
    #else
      #define PARTICLES_CR_C         10000.0
    #endif
  #endif

  #ifndef PARTICLES_CR_FEEDBACK
    #define PARTICLES_CR_FEEDBACK  YES  /* Enable/disable particle feedback on the fluid */
  #endif

  #ifndef PARTICLES_CR_UPWIND_FLUX
    #define PARTICLES_CR_UPWIND_FLUX  NO
  #endif
  
  #if PARTICLES_CR_FEEDBACK == NO
    #undef PARTICLES_CR_UPWIND_FLUX
    #define PARTICLES_CR_UPWIND_FLUX  NO
  #endif    

  #ifndef PARTICLES_CR_E_MC             /* Charge to mass ratio (CR) */
    #define PARTICLES_CR_E_MC      1.0    
  #endif

  #ifndef PARTICLES_CR_E_MC_GAS         
    #define PARTICLES_CR_E_MC_GAS  1.0   /* Charge to mass ratio (fluid) */  
  #endif

  #ifndef PARTICLES_CR_LARMOR_EPS
    #define PARTICLES_CR_LARMOR_EPS  0.3 /* Fraction of Larmor time (1/OmegaL) covered */
                                         /* during one hydro time step */
  #endif

  #ifndef PARTICLES_CR_NCELLS_EPS         /* Maximum number of zones crossed  */
    #define PARTICLES_CR_NCELLS_EPS  1.8  /* during one sub-cycle             */
  #endif

  #ifndef PARTICLES_CR_NSUB
    #define PARTICLES_CR_NSUB      5   /* Max number of particles time steps */
                                       /* per hydro step. When set to a negative */
                                       /* integer, force the sub-stepping to be  */
                                       /* the same. */
  #endif

  #ifndef PARTICLES_CR_PREDICTOR         
    #define PARTICLES_CR_PREDICTOR  2  /* Predictor step in mover.        */
                                         /* This makes the scheme 2nd order */ 
  #endif

  #ifndef PARTICLES_CR_WRITE_4VEL                                                                                                                                    
    #define PARTICLES_CR_WRITE_4VEL  NO                                                                                                                              
  #endif

#endif

/* *********************************************************************
    Particle type: DUST
   ********************************************************************* */

#if PARTICLES_TYPE == DUST

  #ifndef PARTICLES_DUST_FEEDBACK
    #define PARTICLES_DUST_FEEDBACK   YES
  #endif

  #ifndef DUST_RHO
    #define DUST_RHO          1.0   /* Dust/gas density ratio */
  #endif

  #ifndef DUST_SB_ETA_VK
    #define DUST_SB_ETA_VK    0.0
  #endif

  #ifndef DUST_TSTOP
    #define DUST_TSTOP      0.5  /* Stopping time of the particles */
  #endif

  #ifndef DUST_TIME_STEPPING
    #define DUST_TIME_STEPPING  SEMI_IMPLICIT
  #endif
#endif

/* *********************************************************************
    Particle type: LAGRANGIAN
   ********************************************************************* */

#if PARTICLES_TYPE == LAGRANGIAN && PARTICLES_LP_SPECTRA == YES
  #ifndef PARTICLES_LP_NEBINS        /* --> PARTICLES_LP_NEBINS */
    #define PARTICLES_LP_NEBINS 100
  #endif

  #ifndef PARTICLES_LP_SHK_THRESHOLD
    #define PARTICLES_LP_SHK_THRESHOLD   50
  #endif

  #ifndef PARTICLES_LP_SPEC_ENERGY
    #define PARTICLES_LP_SPEC_ENERGY    0.01   /* Unit energy [ergs] to scale spectra */
  #endif

  #ifndef PARTICLES_LP_SHK_GRADP
    #define PARTICLES_LP_SHK_GRADP   0.1 /* To ensure the particles dont cross the same shock twice */
  #endif

  #ifndef PARTICLES_LP_NONTH_FRACN
    #define  PARTICLES_LP_NONTH_FRACN  0.01
  #endif

  #ifndef PARTICLES_LP_NONTH_FRACE
    #define  PARTICLES_LP_NONTH_FRACE  0.5
  #endif

  #ifndef PARTICLES_LP_MICROETA
    #define PARTICLES_LP_MICROETA     4.25
  #endif

  #ifndef PARTICLES_LP_ICCMBZ
   #define PARTICLES_LP_ICCMBZ		0.0
  #endif

  #define UNIT_MAGFIELD       (UNIT_VELOCITY*sqrt(4.0*CONST_PI*UNIT_DENSITY))
  #define UNIT_TIME           (UNIT_LENGTH/UNIT_VELOCITY)
  #define SYNCHROTRON_CONST   (0.0015829*UNIT_MAGFIELD*UNIT_MAGFIELD           \
                               *PARTICLES_LP_SPEC_ENERGY*UNIT_TIME)
#endif

/* *********************************************************************
    Global variables (see globals.h for doc)
   ********************************************************************* */

extern long int p_nparticles;
extern long int p_idCounter;
extern int p_nrestart;
extern int p_intStage;

/* *********************************************************************
    Structure definitions
   ********************************************************************* */

#if PARTICLES_TYPE == COSMIC_RAYS || \
    PARTICLES_TYPE == DUST
typedef struct Particle_{
  double   coord[3];     /**< Particle coordinates */
  double   speed[3];     /**< Particle velocity */
  double   coord_old[3]; /**< Particle coordinates at previous time level */
  double   speed_old[3]; /**< Particles velocity at previous time level */
  double   rho;      /**< Mass density of a single particle         */
  float    tinj;     /**< Particle injection time */  
  float    color;    /**< User-supplied real number to distinguish particles */
  int      cell[3];  /**< Indices (i,j,k) of the cell hosting the particle */
  uint32_t id;       /**< Particle id. Numbering is sequential for each rank */
} Particle;
#endif

#if (PARTICLES_TYPE == LAGRANGIAN) && (PARTICLES_LP_SPECTRA == NO) 
typedef struct Particle_{
  double coord[3];     /**< Particle coordinates */
  double speed[3];     /**< Particle velocity */
  double coord_old[3]; /**< Particle coordinates at previous time level */
  double speed_old[3]; /**< Particles velocity at previous time level */
  double density;
  float  tinj;         /**< Particle injection time */  
  float  color;        /**< User-supplied real number to distinguish particles */
  int    cell[3];      /**< Indices (i,j,k) of the cell hosting the particle */
  uint32_t id;         /**< Particle id. Numbering is sequential for each rank */
} Particle;
#endif

#if (PARTICLES_TYPE == LAGRANGIAN) &&  (PARTICLES_LP_SPECTRA == YES)
typedef struct Particle_{
  double coord[3];     /**< Particle coordinates */
  double speed[3];     /**< Particle velocity */
  double coord_old[3]; /**< Particle coordinates at previous time level */
  double speed_old[3]; /**< Particles velocity at previous time level */
  double density;

  double fourvel[4];
  double density_old;
  double pressure;
  double pressure_old;
  double shk_gradp;
  double divv;
  double ca_old;
  double ca;
  double cr_old;
  double cr;
  double cmp_ratio;
  double lorG;
  double lorG_old;
  double nmicro;
  double mag[3];
  double gradp[3];
  double shk_vL[NVAR];
  double shk_vR[NVAR];
  double eng[PARTICLES_LP_NEBINS];
  double chi[PARTICLES_LP_NEBINS];
  double eng_old[PARTICLES_LP_NEBINS];
  double chi_old[PARTICLES_LP_NEBINS];

  float  tinj;         /**< Particle injection time */  
  float  color;        /**< User-supplied real number to distinguish particles */
  int    cell[3];      /**< Indices (i,j,k) of the cell hosting the particle */
  uint32_t id;         /**< Particle id. Numbering is sequential for each rank */

  char shkflag;
  char prev_shkflag;

} Particle;
#endif

typedef struct particleNode_{
  struct Particle_     p;
  struct particleNode_ *next;
  struct particleNode_ *prev;  
} particleNode;


/* ***************************************************************
    Macro definitions
   *************************************************************** */

#ifndef PARTICLES_USE_ARRAY
  #define PARTICLES_USE_ARRAY   NO
#endif

/* Useful macro to loop over particles (ATTENTION: do *NOT* use this
  macro inside loops involving creation / destruction of particles) */
#define PARTICLES_LOOP(a,b)   for (a = b; a != NULL; a = a->next)

/* ***************************************************************
    Function Prototypes
   *************************************************************** */

void    Particles_Boundary(Data *, Grid *);
void    Particles_BoundaryExchange(Data *, Grid *);
int     Particles_BoundaryCheck(Particle *p, Grid *grid);
long    Particles_CheckAll   (particleNode *, int, Grid *);
int     Particles_CheckSingle(Particle *, int, Grid *);

#if PARTICLES_TYPE == COSMIC_RAYS
void    Particles_CR_ComputeCurrent(const Data *, Grid *);
void    Particles_CR_ComputeForce(Data_Arr, const Data *, Grid *);
void    Particles_CR_ConservativeFeedback(Data_Arr, Data_Arr, double, RBox *);
void    Particles_CR_EMFields(double *, double *, double *);
void    Particles_CR_Flux(const State *, int, int);
void    Particles_CR_States1DCopy(const Data *, const Sweep *, int, int);
void    Particles_CR_StatesSource(const Sweep *, double, int, int, Grid *);
void Particles_CR_StatesSourceOld(const Sweep *sweep, double dt,
                                int beg, int end, Grid *grid);

void    Particles_CR_Update(Data *, timeStep *, double, Grid *);
#endif

void    Particles_Density(Particle *, double *);
void    Particles_Deposit(particleNode *, void (*Func)(Particle *, double *),
                         Data_Arr, int, Grid *);
void    Particles_DepositBoundaryExchange(Data_Arr, int, Grid *);
void    Particles_Destroy(particleNode *, Data *);
void    Particles_Display(Particle*);

#if PARTICLES_TYPE == DUST
void    Particles_Dust_ComputeForce (Data_Arr, Data *, Grid *);
void    Particles_Dust_StatesSource(Sweep *, Data_Arr, Data_Arr, double, int, int, Grid *);

void    Particles_Dust_ConservativeFeedback(Data_Arr, Data_Arr, Data_Arr, double, RBox *);
void    Particles_Dust_Update(Data *, timeStep *, double, Grid *);
#endif

void    Particles_GetWeights (Particle *, int *, double ***, Grid *);

void    Particles_Init(Data *, Grid *);
void    Particles_Inject(Data *d, Grid *grid);
int     Particles_Insert(Particle *, Data *, char, Grid *);
double  Particles_Interpolate(double ***, double ***, int *);

void    Particles_ListToArray (Data *d);

void    Particles_LoadRandom(double *xbeg, double *xend,
                          double (*DistribFunc)(double, double, double),
                          double coor[]);
void    Particles_LoadUniform(int, int, double *, double *, double *);

int     Particles_LocateCell(double *, int *, Grid *);

//Spectra Routines
#if PARTICLES_TYPE == LAGRANGIAN
void    Particles_LP_Predictor(Data *, timeStep *, double, Grid *);
void    Particles_LP_Corrector(Data *, timeStep *, double, Grid *);
#if PARTICLES_LP_SPECTRA == YES
void    Particles_LP_FixValue(Particle*, Data *, Grid *);
void    Particles_LP_FlagShock(Data *, float ***, Grid *);  
void    Particles_LP_GradP(double*,int, double ***u[], Grid *, int indici[]);
void    Particles_LP_IC_Emissivity(Particle *, double, double, double *);
void    Particles_LP_InitSpectra(Particle*);
void    Particles_LP_Spectra(Data *, float ***, Particle*, Grid *, double);  // Useless
void    Particles_LP_UpdateSpectra(Data *, double, Grid *);

void    Particles_LP_Sync_Emissivity(Particle*, double, double, double *, double *, double *);
#endif
#endif

int     Particles_Number(particleNode *);
void    Particles_Restart(Data *, int, Grid *);
Particle *Particles_Select(particleNode *, int);
void    Particles_Set(Data *, Grid *);
void    Particles_SetID(particleNode *);
void    Particles_SetOutput (Data *, Runtime *);

void    Particles_ShowList(particleNode *, int);

void    Particles_UserDefBoundary(Data *d, int, Grid *);

void    Particles_WriteBinary(particleNode *, double, Output *, char *);
void    Particles_WriteData(Data *d, Output *, Grid *);
void    Particles_WriteTab   (particleNode*, char filename[]);
void    Particles_WriteTrajectory (Particle *, char);
void    Particles_WriteVTK   (particleNode*, Output *, char filename[]);

#ifdef PARALLEL
 extern MPI_Datatype MPI_PARTICLE;
 extern MPI_Datatype PartOutputType;
 void Particles_StructDatatype();
#endif
