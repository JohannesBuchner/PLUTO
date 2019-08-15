/* *********************************************************************
   PLUTO function prototypes
   ********************************************************************* */

int    AdvanceStep(Data *, Riemann_Solver *, timeStep *, Grid *);
void   AdvectFlux (const Sweep *, int, int, Grid *);
void   AMR_StoreFlux (double **, double **ViF, int,
                      int, int, int, int, Grid *);

void   Analysis (const Data *, Grid *);

void   Boundary    (const Data *, int, Grid *);

char     *Array1D (int, size_t);
char    **Array2D (int, int, size_t);
char   ***Array3D (int, int, int, size_t);
char  ****Array4D (int, int, int, int, size_t);
double ***ArrayBox(long int, long int, long int, long int, long int, long int);
double ***ArrayBoxMap (int, int, int, int, int, int, double *);
double ***ArrayMap (int, int, int, double *);
unsigned char ***ArrayCharMap (int, int, int, unsigned char *);

void BodyForceVectorGet (double **v, double **g,
                         double *, double *, double *, int beg, int end, 
                         Grid *grid);
void BodyForcePotentialGet (double **v, double **gphi,
                            double *phi_p, double *phi_c,
                            int beg, int end, Grid *grid);


double BodyForcePotential(double, double, double);
void   BodyForceVector(double *, double *, double, double, double);

void   ChangeOutputVar (void);
void   CharTracingStep(const Sweep *, int, int, Grid *);
void   CheckPrimStates (double **, double **, double **, int, int);
int    CheckNaN (double **, int, int, int);
void   ComputeUserVar (const Data *, Grid *);
float  ***Convert_dbl2flt (double ***, double, int);

void  ConsToPrim3D(Data_Arr, Data_Arr, unsigned char ***, RBox *);
void  CreateImage (char *);

void  ComputeEntropy      (const Data *, Grid *);
void  EntropySwitch       (const Data *, Grid *);
void  EntropyOhmicHeating (const Data *, Data_Arr, double, Grid *);

int    FileClose  (FILE *, int);
int    FileDelete (char *);
FILE  *FileOpen  (char *, int, char *);
void   FileReadData  (void *, size_t, int, FILE *, int, int);
void   FileWriteData (void *, size_t, int, FILE *, int);
void   FileWriteHeader(char *buffer, char fname[], int mode);
void   FileWriteArray(void *, long int, int, size_t, char *);

void  FreeArray1D (void *);
void  FreeArray2D (void **);
void  FreeArray3D (void ***);
void  FreeArray4D (void ****);
void  FreeArrayBox(double ***, long, long, long);
void  FreeArrayBoxMap (double ***, int, int, int, int, int, int);
void  FreeArrayMap (double ***);

void  FreeArrayCharMap(unsigned char ***);

#ifdef FINITE_DIFFERENCE  
 Riemann_Solver FD_Flux;
 Reconstruct MP5_Reconstruct, PPM_Reconstruct, LIMO3_Reconstruct,
             WENOZ_Reconstruct, WENO3_Reconstruct;
 void FD_GetMaxEigenvalues (const Data *d, Sweep *sweep, Grid *grid);
#endif

void FindShock (const Data *, Grid *);
void FlagShock (const Data *, Grid *);
void Flatten (const Sweep *, int, int, Grid *);
void FreeGrid (Grid *);

void     StateStructAllocate (State *);

void     GetAreaFlux (const Sweep *, double **, double **, int, int, Grid *);
void     GetCGSUnits (double *);
Image   *GetImage (char *);
double  *GetInverse_dl (const Grid *);
int      GetNghost (void);
void     GetNeighbourRanks (Grid *, int **);
void     GetOutputFrequency(Output *, const char *);
int      GetOutputVarNames(int, char *var_names[NVAR]);

double  ***GetUserVar (char *);

void    GnuplotSetting(Grid *);

void HancockStep    (const Sweep *, int, int, Grid *);

int LocateIndex(double *, int, int, double);

char  *IndentString();

void Init (double *, double, double, double);
void InitDomain (Data *, Grid *);

void Initialize(int argc, char *argv[], Data *, Runtime *, Grid *, Cmd_Line *);

void   InternalBoundaryReset (const Sweep *, timeStep *, int, int, Grid *);

void   InputDataClose(int);
double InputDataInterpolate (int, double, double, double);
int    InputDataOpen(char *, char *, char *, int);
void   InputDataReadSlice(int, int);

int    IsLittleEndian (void);

void   MakeState (Sweep *);
void   MakeGeometry (Grid *);
double MeanMolecularWeight(double *);
double Median (double a, double b, double c);

void   OutputLogPre  (Data *, timeStep *, Runtime *, Grid *);
void   OutputLogPost (Data *, timeStep *, Runtime *, Grid *);

void   ParabolicArrays(const Data *, int *, int, Grid *);
void   ParabolicFlux(Data_Arr, Data_Arr, double ***, const Sweep *,
                     double **, int, int, Grid *);
double ParabolicRHS   (const Data *, Data_Arr, RBox *, double **, int, double, Grid *);
void   ParabolicUpdate(const Data *, Data_Arr, RBox *, double **, double, timeStep *, Grid *);

void   ParseCmdLineArgs (int, char *argv[], char *, Cmd_Line *);

int    ParamFileRead    (char *);
char  *ParamFileGet     (const char *, int );
int    ParamExist       (const char *);
int    ParamFileHasBoth (const char *, const char *);

void   PrimToChar (double **, double *, double *); 
void   PrimToCons3D(Data_Arr, Data_Arr, RBox *);

void   RBoxDefine(int, int, int, int, int, int, int, RBox *);
void   RBoxSetDirections(RBox *, int);
void   RBoxShow(RBox *);


void   ReadHDF5 (Output *output, Grid *grid);
void   ResetState (const Data *, Sweep *, Grid *);
void   RestartFromFile (Runtime *, int, int, Grid *);
void   RestartDump     (Runtime *);
void   RestartGet      (Runtime *, int, int, int);

void   RightHandSide (const Sweep *, timeStep *, int, int, double, Grid *);
void   RightHandSideSource (const Sweep *, timeStep *, int, int, double,
                          double *, Grid *);
void     RKC (const Data *d, double, timeStep *, Grid *);
void     RKL (const Data *d, double, timeStep *, Grid *);

Runtime *RuntimeGet(void);
int      RuntimeSetup  (Runtime *, Cmd_Line *, char *);
void     RuntimeSet(Runtime *runtime);

void StoreAMRFlux (double **flux, double **aflux, int sign,
                    int nvar_beg, int nvar_end, int beg, int end, Grid *grid);


void SetColorMap (unsigned char *, unsigned char *, unsigned char *, char *);
void SetDefaultVarNames(Output *);
int  SetOutputVar (char *, int, int);

int  SetLogFile(char *, Cmd_Line *);
Riemann_Solver *SetSolver (const char *);
void SetGrid (Runtime *, Grid *);
void SetJetDomain   (const Data *, int, int, Grid *);
void SetOutput (Data *d, Runtime *input);
void SetVectorIndices (int);

void Show (double **, int);
void ShowMatrix(double **, int, double);
void ShowVector (double *, int);
int  StringArrayIndex (char *str_arr[], int, char *);

void SymmetryCheck (Data_Arr, int, RBox *);
void ShowConfig(int, char *a[], char *);
void ShowUnits ();
void SplitSource (const Data *, double, timeStep *, Grid *);
void STS (const Data *d, double, timeStep *, Grid *);


void Startup    (Data *, Grid *);
void States     (const Sweep *, int, int, Grid *);
void SwapEndian (void *, const int); 

void UnsetJetDomain (const Data *, int, Grid *);
void UpdateStage(Data *, Data_Arr, double **, Riemann_Solver *,
                 double, timeStep *, Grid *);
void UserDefBoundary (const Data *, RBox *, int,  Grid *); 

//void VectorPotentialDiff (double *, int, int, int, Grid *);
void VectorPotentialDiff (double *, Data *, int, int, int, Grid *);

void Where (int, Grid *);
void WriteData (const Data *, Output *, Grid *);
void WriteHDF5        (Output *output, Grid *grid);
void WriteVTK_Header (FILE *, Grid *);
void WriteVTK_Vector (FILE *, Data_Arr, double, char *, Grid *);
void WriteVTK_Scalar (FILE *, double ***, double, char *, Grid *);
void WriteVTKProcFile (double ***, int, int, int, char *);

void WriteTabArray (Output *, char *, Grid *);
void WritePPM (double ***, char *, char *, Grid *);
void WritePNG (double ***, char *, char *, Grid *);

#define ARRAY_1D(nx,type)          (type    *)Array1D(nx,sizeof(type))
#define ARRAY_2D(nx,ny,type)       (type   **)Array2D(nx,ny,sizeof(type))
#define ARRAY_3D(nx,ny,nz,type)    (type  ***)Array3D(nx,ny,nz,sizeof(type))
#define ARRAY_4D(nx,ny,nz,nv,type) (type ****)Array4D(nx,ny,nz,nv,sizeof(type))

/* ---------------------------------------------------------------------
            Prototyping for standard output/debugging
   --------------------------------------------------------------------- */

void print  (const char *fmt, ...);
void Trace (double);

/* --------------------------------------------------------------------- 
          Prototype for cooling functions
   --------------------------------------------------------------------- */
   
#if COOLING != NO
 void Numerical_Jacobian (double *v, double **J);
 void Jacobian (double *v, double *rhs, double **dfdy);
 void CoolingSource (const Data *, double, timeStep *, Grid *);
 #if COOLING == POWER_LAW
  void  PowerLawCooling (Data_Arr, double, timeStep *, Grid *);
 #endif
 /* move the following elsewhere ?  */
/*
double SolveODE_CK45  (double *, double *, double *, double, double);
double SolveODE_RKF23 (double *, double *, double *, double, double);
double SolveODE_RKF12 (double *, double *, double *, double, double);
*/
 double SolveODE_CK45  (double *, double *, double *, double, double, intList *);
 double SolveODE_RKF23 (double *, double *, double *, double, double, intList *);
 double SolveODE_RKF12 (double *, double *, double *, double, double, intList *);

 double SolveODE_ROS34 (double *, double *, double *, double, double);
 double SolveODE_RK4   (double *, double *, double *, double, intList *);
 double SolveODE_RK2   (double *, double *, double *, double);
 double SolveODE_Euler (double *, double *, double *, double);
#endif

/* ----------------------------------------------
           functions in tools.c 
   ---------------------------------------------- */

void PlutoError (int, char *);

double Length_1 (int i, int j, int k, Grid *);
double Length_2 (int i, int j, int k, Grid *);
double Length_3 (int i, int j, int k, Grid *);

#if UPDATE_VECTOR_POTENTIAL == YES
 void VectorPotentialUpdate (const Data *d, const void *vp, 
                             const Sweep *sweep, const Grid *grid);
#endif

/* --------------- EOS ------------------------- */

void Enthalpy    (double **, double *, int, int);
void Entropy     (double **, double *, int, int);
void SoundSpeed2 (const State *, int, int, int, Grid *);

double GetEntropy (double x);

/* --- New stuffs -- */

void   WriteAsciiFile (char *fname, double *q, int nvar);


