#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double x,y;
  double c,s,alpha;

  g_gamma = g_inputParam[GAMMA_EOS];

  alpha = atan(2.0);
  
  x = x1; 
  y = x2;

  c = cos(alpha);
  s = sin(alpha);
  
  if ( (s*y + c*(x - 0.5)) < 1.e-6){

    us[RHO] = g_inputParam[RHO_LEFT];
    us[VX1] = g_inputParam[VX_LEFT]*c - g_inputParam[VY_LEFT]*s;
    us[VX2] = g_inputParam[VY_LEFT]*c + g_inputParam[VX_LEFT]*s;
    us[VX3] = g_inputParam[VZ_LEFT];
    us[BX1] = g_inputParam[BX_CONST]*c - g_inputParam[BY_LEFT]*s;
    us[BX2] = g_inputParam[BY_LEFT]*c  + g_inputParam[BX_CONST]*s;
    us[BX3] = g_inputParam[BZ_LEFT];
    us[PRS] = g_inputParam[PR_LEFT];

  }else{

    us[RHO] = g_inputParam[RHO_RIGHT];
    us[VX1] = g_inputParam[VX_RIGHT]*c - g_inputParam[VY_RIGHT]*s;
    us[VX2] = g_inputParam[VY_RIGHT]*c + g_inputParam[VX_RIGHT]*s;
    us[VX3] = g_inputParam[VZ_RIGHT];
    us[BX1] = g_inputParam[BX_CONST]*c - g_inputParam[BY_RIGHT]*s;
    us[BX2] = g_inputParam[BY_RIGHT]*c + g_inputParam[BX_CONST]*s;
    us[BX3] = g_inputParam[BZ_RIGHT];
    us[PRS] = g_inputParam[PR_RIGHT];
  }

  us[BX1] /= sqrt(4.0*CONST_PI);
  us[BX2] /= sqrt(4.0*CONST_PI);
  us[BX3] /= sqrt(4.0*CONST_PI);

  us[AX1] = us[AX2] = 0.0;
  us[AX3] = y*us[BX1] - x*us[BX2];
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
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
 * \param [in/out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies on which side boundary conditions need 
 *                    to be assigned. side can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  x1, x2, x3;

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    k = 0;
    for (j = JBEG;  j--;  ) {
    for (i = 1; i < NX1_TOT; i++){
      for (nv = NVAR; nv--; ){
        d->Vc[nv][k][j][i] = d->Vc[nv][k][j + 1][i - 1];
      }
      #ifdef STAGGERED_MHD
       d->Vs[BX1s][k][j][i] = d->Vs[BX1s][k][j + 1][i - 1];
      #endif
    }}
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    k = 0;
    for (j = JEND + 1; j < NX2_TOT; j++) {
    for (i = 0; i < NX1_TOT - 1; i++){
      for (nv = NVAR; nv--; ){
        d->Vc[nv][k][j][i] = d->Vc[nv][k][j - 1][i + 1];
      }
      #ifdef STAGGERED_MHD
       d->Vs[BX1s][k][j][i] = d->Vs[BX1s][k][j - 1][i + 1];
      #endif
    }}
  }
}

