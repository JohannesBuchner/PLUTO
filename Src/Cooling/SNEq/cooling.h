/* ############################################################
      
     FILE:     cooling.h

     PURPOSE:  contains common definitions for the 
               whole CODE

   ############################################################ */

#define NIONS  1
#define X_HI   NFLX   

/* These are used in radiat.c and MeanMolecularWeight() function */

#define frac_Z   1.e-3   /*   = N(Z) / N(H), fractional number density of metals (Z)
                                with respect to hydrogen (H) */ 
#define frac_He  0.082   /*   = N(He) / N(H), fractional number density of helium (He)
                                with respect to hydrogen (H) */ 

double GetMaxRate (double *, double *, double);
double H_MassFrac (void);
double CompEquil (double, double, double *);
void   Radiat (double *, double *);










