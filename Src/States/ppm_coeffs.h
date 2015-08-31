#ifndef PPM_ORDER
  #define PPM_ORDER  4
#endif

#define POLY_2(a0,a1,a2,x)     \
        ( a0 + x*(a1 + x*a2) )
#define POLY_3(a0,a1,a2,a3,x)  \
        ( a0 + x*(a1 + x*(a2 + x*a3)) )
#define POLY_4(a0,a1,a2,a3,a4,x) \
        ( a0 + x*(a1 + x*(a2 + x*(a3 + x*a4))) )
#define POLY_5(a0,a1,a2,a3,a4,a5,x)  \
        ( a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*a5)))) )
#define POLY_6(a0,a1,a2,a3,a4,a5,a6,x) \
        ( a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*(a5 + x*a6)))))  )
#define POLY_7(a0,a1,a2,a3,a4,a5,a6,a7,x) \
        ( a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*(a5 + x*(a6 + x*a7)))))) )
 
#define POLY_8(a0,a1,a2,a3,a4,a5,a6,a7,a8,x) \
        ( a0 + x*(a1 + x*(a2 + x*(a3 + x*POLY_4(a4, a5, a6, a7, a8, x)))) )

#define POLY_10(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,x) \
( a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*POLY_5(a5, a6, a7, a8, a9, a10, x))))) )

/* ********************************************************************* */
/*! Simple structure used to retrieve 1D reconstruction weights 
    (c, w, d) used by piecewise linear interpolation (see states_plm.c)
   ********************************************************************* */

typedef struct PPM_COEFFS{
  double **wp;
  double **wm;
  double *hp;
  double *hm;
} PPM_Coeffs;

void PPM_CoefficientsSet(Grid *grid);
void PPM_CoefficientsGet(PPM_Coeffs*, int);


PPM_Coeffs* PPM_Coefficients(int action, Grid *grid);
void PPM_Q6_Coeffs(double *hp, double *hm, int dir, Grid *grid);
