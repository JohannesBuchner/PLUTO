/* ############################################################

     FILE:     cooling.h

     PURPOSE:  contains shared definitions with scope
               limited to the cooling module ONLY

   ############################################################ */

/* ##############################################################

                   P R O T O T Y P I N G

   ############################################################## */

void   CompEquil(double n, double T, double *v);
double GetMaxRate (double *, double *, double);
void   Radiat (double *, double *);
void   NormalizeIons (double *);
void   H2RateTables(double, double *);

/* ############################################################

         Fractions and Atomic Wts.

 ############################################################## */

#define NIONS    3
#define X_HI     (NFLX)
#define X_H2     (NFLX + 1)
#define X_HII    (NFLX + 2)
