#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  static int first_call = 1;
  double Bx, By, Bz, B0;
  double vx, vy, vz, vA;
  double eta, w, kx, scrh;
  double c,s; 

  g_gamma = 4.0/3.0;

  B0 = eta = 1.0;
  us[RHO] = 1.0;
  us[PRS] = 1.0;
  
  w  = us[RHO] + g_gamma*us[PRS]/(g_gamma - 1.0);
  vA  = B0*B0/(w + B0*B0*(1.0 + eta*eta));
  vA /= 0.5*(1.0 + sqrt(1.0 - 4.0*eta*eta*vA*vA));
  vA  = sqrt(vA);

  #if DIMENSIONS == 1
   c = 1.0;
   s = 0.0;
   kx = 2.0*CONST_PI*x1;
  #elif DIMENSIONS == 2
   c = 1.0/sqrt(2.0);
   s = 1.0/sqrt(2.0);
   kx = 2.0*CONST_PI*sqrt(2.0)*(c*x1 + s*x2);
  #endif

/* -------------------------------------------------
             define 1D solution   
   ------------------------------------------------- */

  Bx = B0;
  By = eta*B0*cos(kx);
  Bz = eta*B0*sin(kx);
  
  vx = 0.0;
  vy = -vA*By/B0;
  vz = -vA*Bz/B0;

/* -------------------------------------------------
               rotate solution
   ------------------------------------------------- */

  us[VX1] = vx*c - vy*s;
  us[VX2] = vx*s + vy*c;
  us[VX3] = vz;
 
  us[BX1] = Bx*c - By*s;
  us[BX2] = Bx*s + By*c;
  us[BX3] = Bz;

  us[AX1] = 0.0;
  us[AX2] = 0.0;
  us[AX3] = B0*((-x1*s + x2*c) - eta*sin(kx)/(2.0*CONST_PI*sqrt(2.0)));

  if (first_call == 1){
    print1 ("vA = %18.12e\n",vA);
    print1 ("T  = %18.12e\n",c/vA);
    first_call = 0;
  }
/*
restart;
B[t] := eta*B[0]*cos(k*(c*x+s*y));
B[x] := B[0]*c - B[t]*s;
B[y] := B[0]*s + B[t]*c;

A[z]:= B[0]*((-x*s+y*c) - eta*sin(k*(c*x+s*y))/k);
simplify( diff(A[z],y) - B[x]);
simplify(-diff(A[z],x) - B[y]);
*/


}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ************************************************************** */
void BACKGROUND_FIELD (real x1, real x2, real x3, real *B0)
/* 
 *
 * PURPOSE
 *
 *   Define the component of a static, curl-free background 
 *   magnetic field.
 *
 *
 * ARGUMENTS
 *
 *   x1, x2, x3  (IN)    coordinates
 *
 *   B0         (OUT)    vector component of the background field.
 *
 *
 **************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

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
{ }

