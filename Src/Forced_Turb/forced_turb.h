#define ForcedTurb_MAXMODES     1000
  
#ifndef FORCED_TURB_FREQ
  #define FORCED_TURB_FREQ      1
#endif

#ifndef FORCED_TURB_DECAY
  #define FORCED_TURB_DECAY     0.5
#endif

#ifndef FORCED_TURB_KMIN
  #define FORCED_TURB_KMIN      6.283
#endif

#ifndef FORCED_TURB_KMAX
  #define FORCED_TURB_KMAX      18.95
#endif
  
#ifndef  FORCED_TURB_ENERGY
  #define FORCED_TURB_ENERGY    2.0e-3
#endif

#ifndef  FORCED_TURB_WEIGHT
  #define FORCED_TURB_WEIGHT	1.0
#endif
  
typedef struct ForcedTurb{
    int NModes;
    int SpectForm;
    int StirFreq;
    double StirDecay;
    double StirEnergy;
    double StirKMin;
    double StirKMax;
    double SolveWt;
    double SolveWtNorm;
    double OUVar;
    double *OUPhases;
    double *StirAmpl;
    double **Mode;
    double **aka;
    double **akb;
    double ****Acc;
}ForcedTurb;

void   ForcedTurb_ComputeAcceleration(ForcedTurb *, Grid *);

void   ForcedTurb_Init(ForcedTurb *);
double ForcedTurb_GenRandNum(void);
void   ForcedTurb_OUNoiseInit(double *, int, double);
void   ForcedTurb_OUNoiseUpdate(double *, int, double, double, double);
void   ForcedTurb_CalcPhases(ForcedTurb *);
void   ForcedTurb_CorrectRHS(const Data *, const Sweep *, int, int, double, Grid *);



