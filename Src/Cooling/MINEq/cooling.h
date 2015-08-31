/* ############################################################
      
     FILE:     cooling.h

     PURPOSE:  contains common definitions for the 
               whole CODE

     Notice: the order is absolutely important and MUST NOT
             be changed !!!
   ############################################################ */


#define C_IONS  3   /* in [1,5] */
#define N_IONS  3   /* in [1,5] */
#define O_IONS  3   /* in [1,5] */
#define Ne_IONS 3   /* in [1,5] */
#define S_IONS  3   /* in [1,5] */
#define Fe_IONS 0   /* in [0,3] */

#if C_IONS == 0
 #define C_EXPAND(a,b,c,d,e)  
#elif C_IONS == 1       
 #define C_EXPAND(a,b,c,d,e)  ,a
#elif C_IONS == 2     
 #define C_EXPAND(a,b,c,d,e)  ,a,b
#elif C_IONS == 3       
 #define C_EXPAND(a,b,c,d,e)  ,a,b,c 
#elif C_IONS == 4
 #define C_EXPAND(a,b,c,d,e)  ,a,b,c,d
#elif C_IONS == 5      
 #define C_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
#endif

#if N_IONS == 0
 #define N_EXPAND(a,b,c,d,e)  
#elif N_IONS == 1      
 #define N_EXPAND(a,b,c,d,e)  ,a
#elif N_IONS == 2     
 #define N_EXPAND(a,b,c,d,e)  ,a,b
#elif N_IONS == 3       
 #define N_EXPAND(a,b,c,d,e)  ,a,b,c 
#elif N_IONS == 4
 #define N_EXPAND(a,b,c,d,e)  ,a,b,c,d
#elif N_IONS == 5      
 #define N_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
#endif

#if O_IONS == 0
 #define O_EXPAND(a,b,c,d,e)  
#elif O_IONS == 1      
 #define O_EXPAND(a,b,c,d,e)  ,a
#elif O_IONS == 2     
 #define O_EXPAND(a,b,c,d,e)  ,a,b
#elif O_IONS == 3       
 #define O_EXPAND(a,b,c,d,e)  ,a,b,c 
#elif O_IONS == 4
 #define O_EXPAND(a,b,c,d,e)  ,a,b,c,d
#elif O_IONS == 5      
 #define O_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
#endif

#if Ne_IONS == 0
 #define Ne_EXPAND(a,b,c,d,e)  
#elif Ne_IONS == 1      
 #define Ne_EXPAND(a,b,c,d,e)  ,a
#elif Ne_IONS == 2     
 #define Ne_EXPAND(a,b,c,d,e)  ,a,b
#elif Ne_IONS == 3       
 #define Ne_EXPAND(a,b,c,d,e)  ,a,b,c 
#elif Ne_IONS == 4
 #define Ne_EXPAND(a,b,c,d,e)  ,a,b,c,d
#elif Ne_IONS == 5      
 #define Ne_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
#endif

#if S_IONS == 0
 #define S_EXPAND(a,b,c,d,e)  
#elif S_IONS == 1      
 #define S_EXPAND(a,b,c,d,e)  ,a
#elif S_IONS == 2     
 #define S_EXPAND(a,b,c,d,e)  ,a,b
#elif S_IONS == 3       
 #define S_EXPAND(a,b,c,d,e)  ,a,b,c 
#elif S_IONS == 4
 #define S_EXPAND(a,b,c,d,e)  ,a,b,c,d
#elif S_IONS == 5      
 #define S_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
#endif

#if Fe_IONS == 0
 #define Fe_EXPAND(a,b,c)  
#elif Fe_IONS == 1      
 #define Fe_EXPAND(a,b,c)  ,a
#elif Fe_IONS == 2     
 #define Fe_EXPAND(a,b,c)  ,a,b
#elif Fe_IONS == 3       
 #define Fe_EXPAND(a,b,c)  ,a,b,c 
#endif

/* **********************************************************************
     Ions are labeled progressively, depending on how many ionization 
     stages are effectively included in the network through the previous 
     X_EXPAND macros. 
     Elements are ordered as {H, He, C, N, O, Ne, S, Fe} and must be 
     carefully respected everywhere in the code. 
     Hydrogen and Helium are always included.
   ********************************************************************** */

enum {
  X_HI = NFLX, X_HeI, X_HeII
  C_EXPAND(X_CI, X_CII, X_CIII, X_CIV, X_CV)
  N_EXPAND(X_NI, X_NII, X_NIII, X_NIV, X_NV)
  O_EXPAND(X_OI, X_OII, X_OIII, X_OIV, X_OV)
  Ne_EXPAND(X_NeI, X_NeII, X_NeIII, X_NeIV, X_NeV)
  S_EXPAND(X_SI, X_SII, X_SIII, X_SIV, X_SV)
  Fe_EXPAND(X_FeI, X_FeII, X_FeIII)
};

#define NIONS  (3+C_IONS+N_IONS+O_IONS+Ne_IONS+S_IONS+Fe_IONS)

double GetMaxRate (double *, double *, double);
double H_MassFrac (void);
double CompEquil  (double, double, double *);
double find_N_rho ();
void Radiat (double *, double *);
void NormalizeIons (double *);

