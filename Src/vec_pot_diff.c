/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file 
  \brief Compute magnetic field from vector potential.
  
  The function VectorPotentialDiff() computes either staggered or 
  cell-center magnetic field components by suitable central differencing
  of the vector potential, \f$ \vec{B} = \nabla\times\vec{A} \f$.
  - In Cartesian geometry:
    \f[
    B_x = \pd{A_z}{y} - \pd{A_y}{z} \,,\quad
    B_y = \pd{A_x}{z} - \pd{A_z}{x} \,,\quad
    B_z = \pd{A_y}{x} - \pd{A_x}{y}
    \f]
 
  - In cylindrical geometry:
    \f[
     B_R = \frac{1}{R}\pd{A_z}{\phi}    - \pd{A_\phi}{z}      \,,\quad
     B_z = \frac{1}{R}\left(\pd{(RA_\phi)}{R} - \pd{A_R}{\phi}\right)
    \f]
 
  - In polar geometry:
    \f[
    B_R    = \frac{1}{R}\pd{A_z}{\phi} - \pd{A_\phi}{z}  \,,\quad
    B_\phi = \pd{A_R}{z} - \pd{A_z}{R}             \,,\quad
    B_z    = \frac{1}{R}\pd{(RA_\phi)}{R} - \frac{1}{R}\pd{A_R}{\phi}
    \f]
 
  - In spherical geometry:
    \f[ 
    B_r      =   \frac{1}{r\sin\theta}\pd{(\sin\theta A_\phi)}{\theta} 
               - \frac{1}{r\sin\theta}\pd{A_\theta}{\phi}           \,,\quad
    B_\theta =   \frac{1}{r\sin\theta}\pd{A_r}{\phi} 
               - \frac{1}{r}\pd{(rA_\phi)}{r}                         \,,\quad
    B_\phi   =   \frac{1}{r}\pd{(rA_\theta)}{r}                  
               - \frac{1}{r}\pd{A_r}{\theta}
    \f]
 
  For cell-centered MHD vector potential is compute at the cell-center.
  In the case of staggered MHD, the position of A is edge-centered and
  it is shown below:
  \verbatim
            ______________________
           /                     /|  
          /                     / |  
         /       z face        /  |
        /                     Ax  |
       /                     /    |
      /                     /     Az
     ----------Ay-----------      |
     |                     |      |
     |                     |   y  |
     |                     | face |
     |                     |     / 
     Az       x face      Az    /  
     |                     |   Ax  
     |                     |  /     
     |                     | /
     |                     |/
     ----------Ay-----------
  \endverbatim
 
 
  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date   Sep 24, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if PHYSICS == MHD || PHYSICS == RMHD
/* ********************************************************************* */
void VectorPotentialDiff (double *b, int i, int j, int k, Grid *grid)
/*!
 * Assign face- or cell-centered magnetic field by differentiating
 * the vector potential.
 *
 * \param [out] b  array of magnetic field starting at 0, 
 *               \f$ B_{x_1} = b[0]\,, B_{x_2} = b[1]\,, B_{x_3} = b[2] \f$
 * \param [in]  i  the cell index in the first coordinate direction
 * \param [in]  j  the cell index in the first coordinate direction
 * \param [in]  k  the cell index in the first coordinate direction
 * \param [in] grid pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int    l_convert;
  double dx1, x1, x1p, x1m;
  double dx2, x2, x2p, x2m;
  double dx3, x3, x3p, x3m;
  double bx, by, bz, br, bphi, bth;
  double us_p[256], us_m[256];
  double r_2;
  double x1f, x2f, x3f; /* point at which magnetic field is desired */

  x1 = grid[IDIR].x[i]; 
  x2 = grid[JDIR].x[j]; 
  x3 = grid[KDIR].x[k]; 

  dx1 = grid[IDIR].dx[i];
  dx2 = grid[JDIR].dx[j];
  dx3 = grid[KDIR].dx[k];

  #ifdef STAGGERED_MHD   
   x1p = grid[IDIR].xr[i]; x1m = grid[IDIR].xl[i];
   x2p = grid[JDIR].xr[j]; x2m = grid[JDIR].xl[j];
   x3p = grid[KDIR].xr[k]; x3m = grid[KDIR].xl[k];
   x1f = grid[IDIR].xr[i]; /* for staggered MHD, we compute magnetic */
   x2f = grid[JDIR].xr[j]; /* field at face centers                  */
   x3f = grid[KDIR].xr[k];
   
  #else
   x1p = grid[IDIR].x[i] + dx1; x1m = grid[IDIR].x[i] - dx1;
   x2p = grid[JDIR].x[j] + dx2; x2m = grid[JDIR].x[j] - dx2;
   x3p = grid[KDIR].x[k] + dx3; x3m = grid[KDIR].x[k] - dx3;
   x1f = grid[IDIR].x[i]; /* for cell-centered MHD, we compute magnetic */
   x2f = grid[JDIR].x[j]; /* field at cell-centers                      */
   x3f = grid[KDIR].x[k];
   
   dx1 = 2.*grid[IDIR].dx[i]; /* redefine the spacing between cell-centers */
   dx2 = 2.*grid[JDIR].dx[j];
   dx3 = 2.*grid[KDIR].dx[k];
  #endif 
  
  #if GEOMETRY == CARTESIAN

  /* ---- assign bx at i_f, j, k  ---- */

   Init (us_p, x1f, x2p, x3);
   Init (us_m, x1f, x2m, x3);

   bx = (us_p[AX3] - us_m[AX3])/dx2;

   #if DIMENSIONS == 3
    Init (us_p, x1f, x2, x3p);
    Init (us_m, x1f, x2, x3m);

    bx -= (us_p[AX2] - us_m[AX2])/dx3;
   #endif

  /* ---- assign by at i, j_f, k  ---- */

   Init (us_p, x1p, x2f, x3);
   Init (us_m, x1m, x2f, x3);

   by = -(us_p[AX3] - us_m[AX3])/dx1;

   #if DIMENSIONS == 3
    Init (us_p, x1, x2f, x3p);
    Init (us_m, x1, x2f, x3m);
   
    by += (us_p[AX1] - us_m[AX1])/dx3;
   #endif

  /* ---- assign bz at i, j, k_f  ---- */

   Init (us_p, x1p, x2, x3f);
   Init (us_m, x1m, x2, x3f);
   
   bz = (us_p[AX2] - us_m[AX2])/dx1;

   Init (us_p, x1, x2p, x3f);
   Init (us_m, x1, x2m, x3f);
   
   bz -= (us_p[AX1] - us_m[AX1])/dx2;

   b[0] = bx; b[1] = by; b[2] = bz;

  #elif GEOMETRY == CYLINDRICAL  /* -- only 2D -- */
  
  /* ---- assign br at i_f, j, k  ---- */
  
   Init (us_p, x1f, x2p, x3);
   Init (us_m, x1f, x2m, x3);
   
   br = - (us_p[AX3] - us_m[AX3])/dx2;

  /* ---- assign bz at i, j_f, k  ---- */
     
   Init (us_p, x1p, x2f, x3);
   Init (us_m, x1m, x2f, x3);
   
   bz = (x1p*us_p[AX3] - x1m*us_m[AX3])/(x1*dx1);

  /* ---- assign bphi at i, j, k ---- */  /* -- Only non STAG-- */
  
   #ifdef STAGGERED_MHD   
    bphi = 0.0;
   #else
    Init (us_p, x1, x2p, x3);
    Init (us_m, x1, x2m, x3);
    
    bphi = (us_p[AX1] - us_m[AX1])/dx2;
    
    Init (us_p, x1p, x2, x3);
    Init (us_m, x1m, x2, x3);
    
    bphi -= (us_p[AX2] - us_m[AX2])/dx1;
   #endif
   
   b[0] = br; b[1] = bz; b[2] = bphi;

  #elif GEOMETRY == POLAR 

   /* ---- assign br at i_f, j, k ---- */

   Init (us_p, x1f, x2p, x3);
   Init (us_m, x1f, x2m, x3);
   
   br = (us_p[AX3] - us_m[AX3])/(x1f*dx2);

   Init (us_p, x1f, x2, x3p);
   Init (us_m, x1f, x2, x3m);

   br -= (us_p[AX2] - us_m[AX2])/dx3;

  /* ---- assign bphi at i, j_f, k  ---- */

   Init (us_p, x1p, x2f, x3);
   Init (us_m, x1m, x2f, x3);

   bphi = -(us_p[AX3] - us_m[AX3])/dx1;
     
   #if DIMENSIONS == 3
    Init (us_p, x1, x2f, x3p);
    Init (us_m, x1, x2f, x3m);

    bphi += (us_p[AX1] - us_m[AX1])/dx3;
   #endif

  /* ---- assign bz at i, j, k_f  ---- */

   Init (us_p, x1p, x2, x3f);
   Init (us_m, x1m, x2, x3f);

   bz = (x1p*us_p[AX2] - x1m*us_m[AX2])/(x1*dx1);

   Init (us_p, x1, x2p, x3f);
   Init (us_m, x1, x2m, x3f);

   bz -= (us_p[AX1] - us_m[AX1])/(x1*dx2);

   b[0] = br; b[1] = bphi; b[2] = bz;
 
  #elif GEOMETRY == SPHERICAL

   /* ---- assign br at i_f, j, k ---- */

   Init (us_p, x1f, x2p, x3);
   Init (us_m, x1f, x2m, x3);

   br = (sin(x2p)*us_p[AX3] - sin(x2m)*us_m[AX3])/(x1f*(cos(x2m) - cos(x2p)));

   #if DIMENSIONS == 3
    Init (us_p, x1f, x2, x3p);
    Init (us_m, x1f, x2, x3m);

    br -= 1.0/(x1f*sin(x2)*dx3)*(us_p[AX2] - us_m[AX2]);
   #endif
   
  /* ---- assign btheta at i, j_f, k  ---- */

   Init (us_p, x1p, x2f, x3);
   Init (us_m, x1m, x2f, x3);

   bth = - (x1p*us_p[AX3] - x1m*us_m[AX3])/(x1*dx1);

   #if DIMENSIONS == 3
    Init (us_p, x1, x2f, x3p);
    Init (us_m, x1, x2f, x3m);
   
    bth += (us_p[AX1] - us_m[AX1])/(x1*sin(x2f)*dx3);
   #endif

  /* ---- assign bphi at i, j, k_f  ---- */

   bphi = 0.0;

   Init (us_p, x1p, x2, x3f);
   Init (us_m, x1m, x2, x3f);

   bphi = (x1p*us_p[AX2] - x1m*us_m[AX2])/(x1*dx1);
 
   Init (us_p, x1, x2p, x3f);
   Init (us_m, x1, x2m, x3f);

   bphi -= (us_p[AX1] - us_m[AX1])/(x1*dx2);
   
   b[0] = br; b[1] = bth; b[2] = bphi;

  #endif
  
}
#endif /* PHYSICS == MHD || PHYSICS == RMHD */
