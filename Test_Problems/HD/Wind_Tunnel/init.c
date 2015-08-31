#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
   
   g_gamma = 1.4;

   us[VX1] = (x1>0.6 && x2<0.2 ? 0.0:3.0);
   us[VX2] = 0.0;
   us[VX3] = 0.0;
   us[TRC] = 0.0;

   us[RHO] = 1.4;
   us[PRS] = 1.0;

}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, i1, j1;
  double  *x, *y;

  x = grid[IDIR].x;
  y = grid[JDIR].x;

  if (side == 0 && g_dir == IDIR){
    KDOM_LOOP(k) { 
    JDOM_LOOP(j) { 
      if (y[j] > 0.2) continue;
      IDOM_LOOP(i){ 
        if (x[i] > 0.6 && x[i-1] < 0.6){
          for (i1 = i; i1 <= i + 3; i1++){
            d->Vc[RHO][k][j][i1] =  d->Vc[RHO][k][j][2*i - i1 - 1];
            d->Vc[VX1][k][j][i1] = -d->Vc[VX1][k][j][2*i - i1 - 1];
            d->Vc[VX2][k][j][i1] =  d->Vc[VX2][k][j][2*i - i1 - 1];
            d->Vc[PRS][k][j][i1] =  d->Vc[PRS][k][j][2*i - i1 - 1];
          }
        }
      }
    }}
  }

  if (side == 0 && g_dir == JDIR){
    KDOM_LOOP(k) {
    IDOM_LOOP(i) {
      if (x[i] < 0.6) continue;
      JDOM_LOOP(j){ 
        if (y[j] > 0.2 && y[j-1] < 0.2){
          for (j1 = j - 1; j1 >= j - 4; j1--){
            d->Vc[RHO][k][j1][i] =  d->Vc[RHO][k][2*j - j1 - 1][i];
            d->Vc[VX1][k][j1][i] =  d->Vc[VX1][k][2*j - j1 - 1][i];
            d->Vc[VX2][k][j1][i] = -d->Vc[VX2][k][2*j - j1 - 1][i];
            d->Vc[PRS][k][j1][i] =  d->Vc[PRS][k][2*j - j1 - 1][i];
          }
        }
      }
    }}
  }
       

  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i){ 
      d->Vc[RHO][k][j][i] = 1.4;
      d->Vc[VX1][k][j][i] = 3.0;
      d->Vc[VX2][k][j][i] = 0.0;
      d->Vc[PRS][k][j][i] = 1.0;
    }
  } 

}


