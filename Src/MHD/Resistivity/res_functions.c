/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Compute resistivity.

  Compute resistivity using the same staggering of the current:

  -  eta_x1 at   (i, j+1/2, k+1/2)
  -  eta_x2 at   (i+1/2, j, k+1/2)
  -  eta_x3 at   (i+1/2, j+1/2, k)
 
  Since each eta_xn may depend on the total current, we also need
  to interpolate different current components at the same place.

  \c eta is stored in this function as a static array.

  \authors A. Mignone (mignone@ph.unito.it)\n
           
  \date   Sept 7, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static Data_Arr eta;

/* ********************************************************************* */
void ComputeStaggeredEta(const Data *d, Grid *grid)
/*!
 *
 * \param [in,out]  d       pointer to the PLUTO data structure
 * \param [in]      grid    pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int    i, j, k, nv;
  double Js[3], etas[3], vs[NVAR];
  double ***Jx1, ***Jx2, ***Jx3;
  double *x1,  *x2,  *x3;
  double *x1r, *x2r, *x3r;

  if (eta == NULL){
    eta = ARRAY_4D(3, NX3_TOT,NX2_TOT,NX1_TOT,double);
  }

  Jx1 = d->J[IDIR]; x1 = grid->x[IDIR]; x1r = grid->xr[IDIR];
  Jx2 = d->J[JDIR]; x2 = grid->x[JDIR]; x2r = grid->xr[JDIR];
  Jx3 = d->J[KDIR]; x3 = grid->x[KDIR]; x3r = grid->xr[KDIR];

  for (k = 0; k < NX3_TOT-KOFFSET; k++){
  for (j = 0; j < NX2_TOT-JOFFSET; j++){
  for (i = 0; i < NX1_TOT-IOFFSET; i++){

    #if COMPONENTS == 2

  /* -- Compute J and eta_x3 at  i+1/2, j+1/2, k -- */

    Js[IDIR] = 0.0;
    Js[JDIR] = 0.0;
    Js[KDIR] = Jx3[k][j][i];
    VAR_LOOP(nv) vs[nv] = AVERAGE_XY(d->Vc[nv],k,j,i);
    Resistive_eta (vs, x1r[i], x2r[j], x3[k], Js, etas);
    eta[KDIR][k][j][i] = etas[KDIR];

    #elif COMPONENTS == 3

  /* -- Compute J and eta_x1 at  i, j+1/2, k+1/2 -- */

    if (i > 0) {
      Js[IDIR] = Jx1[k][j][i];
      Js[JDIR] = AVERAGE_XY(Jx2,k,j,i-1);
      Js[KDIR] = AVERAGE_XZ(Jx3,k,j,i-1);
      VAR_LOOP(nv) vs[nv] = AVERAGE_YZ(d->Vc[nv],k,j,i);
      Resistive_eta (vs, x1[i], x2r[j], x3r[k], Js, etas);
      eta[IDIR][k][j][i] = etas[IDIR];
    }

  /* -- Compute J and eta_x2 at  i+1/2, j, k+1/2 -- */

    if (j > 0 || DIMENSIONS == 1) {
      Js[IDIR] = AVERAGE_XY(Jx1,k,j-1,i);
      Js[JDIR] = Jx2[k][j][i];
      Js[KDIR] = AVERAGE_YZ(Jx3,k,j-1,i);
      VAR_LOOP(nv) vs[nv] = AVERAGE_XZ(d->Vc[nv],k,j,i);
      Resistive_eta (vs, x1r[i], x2[j], x3r[k], Js, etas);
      eta[JDIR][k][j][i] = etas[JDIR];
    }

  /* -- Compute J and eta_x3 at  i+1/2, j+1/2, k -- */

    if (k > 0 || DIMENSIONS <= 2){
      Js[IDIR] = AVERAGE_XZ(Jx1,k-1,j,i);
      Js[JDIR] = AVERAGE_YZ(Jx2,k-1,j,i);
      Js[KDIR] = Jx3[k][j][i];
      VAR_LOOP(nv) vs[nv] = AVERAGE_XY(d->Vc[nv],k,j,i);
      Resistive_eta (vs, x1r[i], x2r[j], x3[k], Js, etas);
      eta[KDIR][k][j][i] = etas[KDIR];
    }

    #endif

  }}}

}

/* ********************************************************************* */
Data_Arr GetStaggeredEta()
/*
 *********************************************************************** */
{
  return eta;
}
