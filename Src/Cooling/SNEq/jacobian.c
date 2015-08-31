#include "pluto.h"

/* ******************************************************** */
void Jacobian (double *v, double *rhs, double **dfdy)
/*
 *
 *   Compute the jacobian J(k,l) = dfdy
 *
 *   k = row index
 *   l = column index
 *
 *    J(0,0)   J(0,1) ...  J(0, n-1)
 *    J(1,0)   J(1,1) .... J(1, n-1)
 *      .         .           .
 *      .         .           .
 *      .         .           .
 *    J(n-1,0)  ....      J(n-1, n-1)
 *
 *
 *   or, 
 *
 *   +-----------------------+
 *   +              |        |
 *   +              |        |
 *   +              |        |
 *   +    dX'/dX    | dX'/dp |
 *   +     (JXX)    |  (JXp) |
 *   +              |        |
 *   +              |        |
 *   +--------------+--------+
 *   +   dp'/dX     | dp'/dp | 
 *   +    (JpX)     |  Jpp   |
 *   +-----------------------+
 *
 *
 *
 ********************************************************** */
{
  print (" ! Jacobian not defined \n");
  QUIT_PLUTO(1);
}

