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
  int i, j, k;  
  
  DOM_LOOP(k,j,i){
  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

  SetOutputVar("rho", FLT_OUTPUT, NO);
  SetOutputVar("prs", FLT_OUTPUT, NO);
  SetOutputVar("vx1", FLT_OUTPUT, NO);
  SetOutputVar("vx2", FLT_OUTPUT, NO);
  SetOutputVar("vx3", FLT_OUTPUT, NO);
}





