#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k, nv;  
  double ***tmp, T, v[NVAR], mu = 0.5;

  tmp = GetUserVar("tmp");
  DOM_LOOP(k,j,i){
    VAR_LOOP(nv) v[nv] = d->Vc[nv][k][j][i];
#if EOS == IDEAL
    T = v[PRS]/v[RHO]*KELVIN*mu;
#elif EOS == PVTE_LAW
    GetPV_Temperature(v, &T);
#endif
    tmp[k][j][i] = T;
  }
}
/* ************************************************************* */
void ChangeDumpVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

}





