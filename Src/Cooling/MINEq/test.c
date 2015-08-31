#include<stdio.h>

#define C_IONS  3
#define N_IONS  0


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

enum {
  HI = 10
  C_EXPAND(CI, CII, CIII, CIV, CV)
  N_EXPAND(NI, NII, NIII, NIV, NV)
};

int main()
{
  double coll_ion_P[] = { 0., 0., 1.
                          C_EXPAND(0., 1., 1., 1., 1.)
                          N_EXPAND(0., 0., 1., 1., 0.), 0., 1., 0. }; 
  double ab[] = {1 C_EXPAND(0.1,0.2,0.3,0.4, 0.5 ) N_EXPAND(10.2, 10.3, 10.4, 10.5, 10.6)};
  int a;
  printf ("CIII = %d, %d\n",CIII, CI);

  a = (1,0);
  return (0);

}


