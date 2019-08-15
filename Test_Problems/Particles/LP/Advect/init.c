/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief A simple advection test for Particles as Lagrangian tracers.
 
  Set the initial condition for a advection test along with Lagrangian tracers.
 
  The input parameters read from pluto.ini are labeled as:
  - <tt>g_inputParam[VEL]/tt>: Advection velocity magnitude
  - <tt>g_inputParam[THETA]/tt>: Spherical co-ordinate Angle : Theta
  - <tt>g_inputParam[PHI]/tt>: Spherical co-ordinate Angle : Phi
 
  The available configuration refer to:
   -#01 Cartesian (2D, HD, VEL=0.25, THETA = PI/2, PHI = 0)
   -#02 Cartesian (3D, HD, VEL=0.25, THETA = PI/2, PHI = 0)
   -#03 Cartesian (3D, HD, VEL=0.25, THETA = PI/2, PHI = PI/4)
   -#04 Cartesian (3D, HD, VEL=0.25, THETA = PI/4, PHI = PI/4)
 
  One particle per cell is initialized uniformly in the numerical
  domain. Spectral evolution is not considered for this test
  problem.
 
 \author B. Vaidya (bvaidya@iiti.ac.in)
         A. Mignone (mignone@ph.unito.it)
 
 \date   June 02, 2018

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  double V_0 = g_inputParam[VEL];
  double phi    = g_inputParam[PHI];
  double theta  = g_inputParam[THETA];
  v[RHO] = 1.0;

  #if HAVE_ENERGY
   v[PRS] = 1.0/g_gamma;
  #endif

 EXPAND( v[VX1] = V_0*sin(theta)*cos(phi); ,
         v[VX2] = V_0*sin(theta)*sin(phi); ,
         v[VX3] = V_0*cos(theta);)

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
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
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
{ }
