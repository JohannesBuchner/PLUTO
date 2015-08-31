/* ############################################################
      
     FILE:     cooling_defs.h

     PURPOSE:  contains shared definitions with scope
               limited to the cooling module ONLY

   ############################################################ */

/* ##############################################################
     
                   P R O T O T Y P I N G 

   ############################################################## */

void   Check_Species (double *);
void   Find_Rates(double T, double Ne, double N, double *v);
double find_N_rho ();
int    Create_Ion_Coeff_Tables(double ***);
int    Create_Losses_Tables(double ***, int *, int *);

/* ############################################################################

                            New structures

   ############################################################################ */

typedef struct COOLING_COEFFICIENTS
{
  double Lrate[NIONS];  /*  Abundance increase rates. This is to be multiplied with N(x-1) - it's the ionization from x-1 to x */
  double Crate[NIONS];  /*  Abundance decrease rates. This is to be multiplied with N(x) - it's the ionization/recombination from x  */
  double Rrate[NIONS];  /*  Abundance increase rates. This is to be multiplied with N(x+1) - it's the recombination from x+1 to x  */
  double N;                    /*  Total number density                 */
  double Ne;                   /*  Electron number density ( cm^{-3} )  */
  double T;                    /*  Temperature (K)                      */
  double La[NIONS];     /*  The coefficient of N_el in Lrate     */
  double Lb[NIONS];     /*  The coefficient of X_H in Lrate      */
  double Lc[NIONS];     /*  The coefficient of X_He in Lrate     */
  double Ca[NIONS];     /*  The coefficient of N_el in Crate     */
  double Cb[NIONS];     /*  The coefficient of X_H in Crate      */
  double Cc[NIONS];     /*  The coefficient of X_He in Crate     */
  double Ra[NIONS];     /*  The coefficient of N_el in Rrate     */
  double Rb[NIONS];     /*  The coefficient of X_H in Rrate      */
  double Rc[NIONS];     /*  The coefficient of X_He in Rrate     */
  double muN, muD;        /*  The numerator and denominator of the mean molecular weight - used for d\mu computations  */
  double de[NIONS];     /*  The radiative losses read from cooling tables, interpolated for T and N_el               */
  double de_dne[NIONS]; /*  The slope coefficient in the radiative losses interpolation function of N_el             */

/* -------- NEW ENTRIES -------- */
  
  double dnel_dX[NIONS];
  double fCH, fRH;
  double dmuN_dX[NIONS];
  double dmuD_dX[NIONS];
  double dLIR_dX[NIONS];

} COOL_COEFF;


/* **********************************************************************
          Macros
   ********************************************************************** */

#define dEtoA(x) (12375./(x))
#define to_ev(T) (1.38/1.602*1.e-4*(T))

/* **********************************************************************
                label elements and ions
   ********************************************************************** */

#define el_H   0
#define el_He  1
#define el_C   2
#define el_N   3
#define el_O   4
#define el_Ne  5
#define el_S   6
#define el_Fe  7

/* ############################################################################

        Global variable declaration.  They are defined in radiat.c

   ############################################################################ */

extern const double elem_ab[8];   /* Number fractions   */
extern double elem_ab_sol[7];   /* Number fractions, Solar    */
extern double elem_ab_uni[7];   /* Number fractions, Universe */
extern const double elem_mass[8];  /*   Atomic mass, in a.m.u.   */
extern const int elem_part[31];
extern const double rad_rec_z[31];
extern const double coll_ion_dE[31];
extern COOL_COEFF CoolCoeffs;

/* ****************************************************
   First, the definitions for the 
   Radiative Losses Tables computation
   **************************************************** */


/* ****************************************************
   nlev  = number of levels involved 
   A     = transition probabilities in sec^(-1)
   wgth  = statistical weights
   dE    = energy of the transition in eV 
   N     = Abundance of the ion;
   Ni    = populations of the level;
           notice that we must have \sum Ni = N
   **************************************************** */

typedef struct ION {
   int      nlev, nTom, isMAP, isCV, isH, isCHEB;                 
   double   N;  
   double   wght[16];    /* max number of levels in a ion is 16 */    
   double   Ni[16];
   double   **A, **dE;  
   double   ***omega, Tom[8];
} Ion;

/*  Lagrange interpolation function */
double lagrange (double *x, double *y, double xp, int n, int ii, int jj); 

/* Linear system solver  */
void Solve_System (Ion *X, double Ne, double T);
void Symmetrize_Coeff (Ion *X);

/* Emission lines definitions  */
void INIT_ATOM(Ion *, int);
void HI_INIT (Ion *HIv);
void HeI_INIT (Ion *HeIv);
void HeII_INIT (Ion *HeIIv);
void CI_INIT (Ion *CIv);
void CII_INIT (Ion *CIIv);
void CIII_INIT (Ion *CIIv);
void CIV_INIT (Ion *CIIv);
void CV_INIT (Ion *CIIv);
void NI_INIT (Ion *NIv);
void NII_INIT (Ion *NIIv);
void NIII_INIT (Ion *NIIIv);
void NIV_INIT (Ion *NIIIv);
void NV_INIT (Ion *NIIIv);
void OI_INIT (Ion *OIv);
void OII_INIT (Ion *OIIv);
void OIII_INIT (Ion *OIIIv);
void OIV_INIT (Ion *OIVv);
void OV_INIT (Ion *OVv);
void NeI_INIT (Ion *NeIv);
void NeII_INIT (Ion *NeIIv);
void NeIII_INIT (Ion *NeIIIv);
void NeIV_INIT (Ion *NeIIIv);
void NeV_INIT (Ion *NeIIIv);
void SI_INIT (Ion *SIIv);
void SII_INIT (Ion *SIIv);
void SIII_INIT (Ion *SIIv);
void SIV_INIT (Ion *SIIv);
void SV_INIT (Ion *SIIv);
void FeI_INIT (Ion *FeIv);
void FeII_INIT (Ion *FeIIv);
void FeIII_INIT (Ion *FeIIIv);

/* -----------------------------------------------
    Configure ionization tables parameters here:
    start temperature, temperature step,
    number of steps, number of ion species
   ----------------------------------------------- */

#define  kB        8.617343e-5   

#define I_g_stepNumber    40000
#define I_TSTEP     10.0
#define I_TBEG     500.0
#define I_TEND   (I_TBEG+(I_g_stepNumber-1)*I_TSTEP)

/* --------------------------------------------------------
    Configure cooling tables parameters here:
    electron number density range (C_NeMIN to C_NeMAX),
   -------------------------------------------------------- */

#define C_NeMIN   1.0e-2
#define C_NeMAX   1.0e+8
#define C_NeSTEP   0.12

#define C_TMIN    2000.0
#define C_TMAX    2.0e+5
#define C_TSTEP    0.025






