#include "pluto.h"

/* ********************************************************************* */
void ConvertTo4vel (double **v, int beg, int end)
/*
 *
 *
 *********************************************************************** */
{
  int i;
  double lor;

  for (i = beg; i <= end; i++){
    lor = EXPAND(  v[i][VX1]*v[i][VX1],
                 + v[i][VX2]*v[i][VX2],
                 + v[i][VX3]*v[i][VX3]);
    lor = 1.0/sqrt(1.0 - lor);
  
    EXPAND(v[i][VX1] *= lor;  ,
           v[i][VX2] *= lor;  ,
           v[i][VX3] *= lor;)
  }
}
/* ********************************************************************* */
void ConvertTo3vel (double **v, int beg, int end)
/*
 *
 *
 *********************************************************************** */
{
  int i;
  double lor;

  for (i = beg; i <= end; i++){
    lor = EXPAND(  v[i][VX1]*v[i][VX1],
                 + v[i][VX2]*v[i][VX2],
                 + v[i][VX3]*v[i][VX3]);
    lor = sqrt(1.0 + lor);
  
    EXPAND(v[i][VX1] /= lor;  ,
           v[i][VX2] /= lor;  ,
           v[i][VX3] /= lor;)
  }
}
