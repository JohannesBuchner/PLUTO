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
  #if EOS == PVTE_LAW 
   int i, j, k, nv;  
   double ***tmp, T, v[NVAR];

   tmp = GetUserVar("tmp");
   DOM_LOOP(k,j,i){
     VAR_LOOP(nv) v[nv] = d->Vc[nv][k][j][i];
     GetPV_Temperature(v, &T);
     tmp[k][j][i] = T;
   }
  #endif
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

}





