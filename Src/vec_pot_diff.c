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

  For staggered MHD the convention is the following:
  
  - <tt> Ax1[i,j,k] --> Ax1(0++) </tt> \f$= A_{x_1, i, j+\HALF, k+\HALF}\f$
  - <tt> Ax2[i,j,k] --> Ax2(+0+) </tt> \f$= A_{x_2, i+\HALF, j, k+\HALF}\f$
  - <tt> Ax3[i,j,k] --> Ax3(++0) </tt> \f$= A_{x_3, i+\HALF, j+\HALF, k}\f$

 
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
 
 
  \last change  D. Mukherjee, A. Mignone (dipanjan.mukherjee@unito.it) \n
  \date         May 21, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if PHYSICS == MHD || PHYSICS == RMHD
/* ********************************************************************* */
void VectorPotentialDiff (double *b, Data *d, int i, int j, int k, Grid *grid)
/*!
 * Assign face- or cell-centered magnetic field by differentiating
 * the vector potential.
 *
 * \param [out]    b  array of magnetic field starting at 0, 
 *                 \f$ B_{x_1} = b[0]\,, B_{x_2} = b[1]\,, B_{x_3} = b[2] \f$
 * \param [in/out] Data structure                
 * \param [in]  i  the cell index in the first coordinate direction
 * \param [in]  j  the cell index in the first coordinate direction
 * \param [in]  k  the cell index in the first coordinate direction
 * \param [in] grid pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  double dx1,dx2,dx3;
  double x1,x2,x3;
  double x1p,x2p,x3p;
  double x1m,x2m,x3m;
  double x1f, x2f, x3f; 
  double bx, by, bz, br, bphi, bth;

  double ***A1, ***A2, ***A3;
  double A1_x2p, A1_x3p;
  double A1_x2m, A1_x3m;
  double A2_x1p, A2_x3p;
  double A2_x1m, A2_x3m;
  double A3_x1p, A3_x2p;
  double A3_x1m, A3_x2m;

  x1  = grid->x[IDIR][i]; 
  x2  = grid->x[JDIR][j]; 
  x3  = grid->x[KDIR][k];
 
  dx1 = grid->dx[IDIR][i];
  dx2 = grid->dx[JDIR][j];
  dx3 = grid->dx[KDIR][k];

/* -- Define pointers to A[] -- */
  A1 = d->Ax1;  A2 = d->Ax2; A3 = d->Ax3; 


  x1p = grid->xr[IDIR][i]; x1m = grid->xl[IDIR][i];
  x2p = grid->xr[JDIR][j]; x2m = grid->xl[JDIR][j];
  x3p = grid->xr[KDIR][k]; x3m = grid->xl[KDIR][k];
  x1f = grid->xr[IDIR][i]; /* for staggered MHD, we compute magnetic */
  x2f = grid->xr[JDIR][j]; /* field at face centers                  */
  x3f = grid->xr[KDIR][k];
  
  A1_x2p = A1[k][j][i]; A1_x2m = A1[k][j-1][i];
  A1_x3p = A1[k][j][i]; A1_x3m = A1[k-1][j][i];

  A2_x1p = A2[k][j][i]; A2_x1m = A2[k][j][i-1];
  A2_x3p = A2[k][j][i]; A2_x3m = A2[k-1][j][i];

  A3_x1p = A3[k][j][i]; A3_x1m = A3[k][j][i-1];
  A3_x2p = A3[k][j][i]; A3_x2m = A3[k][j-1][i];
 


#if GEOMETRY == CARTESIAN
/* -----------------------------------------------------
    Compute Bx
    - For cell centred: Bx = dely Az   - delz Ay
    - For Staggered:    Bx =  (Az[++0] - Az[+-0])/dy
                             -(Ay[+0+] - Ay[+0-])/dz
   ----------------------------------------------------- */
  bx = (A3_x2p - A3_x2m)/dx2;

  #if DIMENSIONS == 3
  bx -= (A2_x3p - A2_x3m)/dx3;
  #endif 
  
/* -----------------------------------------------------
    Compute By
    - For cell centred: By = -delx Az  + delz Ax
    - For Staggered:    By = -(Az[++0] - Az[-+0])/dx
                             +(Ax[0++] - Ax[0+-])/dz
   ----------------------------------------------------- */
  by = -(A3_x1p - A3_x1m)/dx1;
  #if DIMENSIONS == 3
  by += (A1_x3p - A1_x3m)/dx3;
  #endif 

/* -----------------------------------------------------
    Compute Bz
    - For cell centred: Bz =  delx Ay  - dely Ax
    - For Staggered:    Bz =  (Ay[+0+] - Ay[-0+])/dx
                             -(Ax[0++] - Ax[0-+])/dy
   ----------------------------------------------------- */
  bz =  (A2_x1p - A2_x1m)/dx1 - (A1_x2p - A1_x2m)/dx2;

  b[0] = bx; b[1] = by; b[2] = bz;

#elif GEOMETRY == CYLINDRICAL  /* -- only 2D -- */
/* -----------------------------------------------------
   Compute Br   
   - For cell centred: Br = - delz Aphi

   Convention for  Components: 
   - Field:      (Br, Bz)
   - Potential:  (Ar, Az, -Aphi) = (AX1, AX2, AX3)

   Note that A3 = -A_\phi since in cylindrical
   coordinates (R,z) are not right-handed
   -----------------------------------------------------*/
   br = -(A3_x2p   - A3_x2m)/dx2;

/* ----------------------------------------------------------
   Compute Bz  
   - For cell centred: Bz = 1/r delr (r*Aphi) 
   - For Staggered:    Bz = (r+ A3[0,0] - r- A3[0,-])/(r dr)
  -----------------------------------------------------------*/
  bz =  (x1p*A3_x1p   - x1m*A3_x1m)/(x1*dx1);

/* ----------------------------------------------
    Compute Bphi  
    - For cell centred: Bphi = delz Ar - delr Az  
    - For Staggered:    Bphi = 0.
   ----------------------------------------------*/
  bphi = 0.0;
  b[0] = br; b[1] = bz; b[2] = bphi;

#elif GEOMETRY == POLAR 
/* ------------------------------------------------------------
   Compute Br 
     Convention for  Components: 
    - Field:      (Br, Bphi, Bz)
    - Potential:  (Ar, Aphi, Az)
 
    - For cell centred: Br =  1/r delphi Az - delz Aphi
    - For Staggered:    Br =  (Az[++0]   - Az[+-0]  )/(r dphi)
                             -(Aphi[+0+] - Aphi[+0-])/dz
  --------------------------------------------------------------*/
  br = (A3_x2p - A3_x2m)/(x1f*dx2);

  #if DIMENSIONS == 3
  br -= (A2_x3p   - A2_x3m)/dx3;
  #endif 

/* -----------------------------------------------------------
    Compute Bphi 
    - For cell centred: Bphi = - delr Az + delz Ar
    - For Staggered:    Bphi =  (Az[++0] - Az[-+0])/dr 
		               -(Ar[0++] - Ar[0+-])/dz
  -------------------------------------------------------------*/
  bphi = -(A3_x1p   - A3_x1m)/dx1;

  #if DIMENSIONS == 3
  bphi += (A1_x3p - A1_x3m)/dx3;
  #endif 

/* ------------------------------------------------------------------
   Compute Bz  
    - For cell centred: Bz = 1/r( delr (r Aphi) - delphi Ar)
    - For Staggered:    Bz = ((r+) Aphi[+0+] - (r-) Aphi[+0-])/(r dr)
		            -(     Ar[0++]   -      Ar[0-+]  )/(r dphi) 
   --------------------------------------------------------------------*/
  bz =  (x1p*A2_x1p - x1m*A2_x1m)/(x1*dx1) - (A1_x2p - A1_x2m)/(x1*dx2);
  b[0] = br; b[1] = bphi; b[2] = bz;

#elif GEOMETRY == SPHERICAL
/* ----------------------------------------------------------------------
    Compute Br 
     Convention for  Components: 
    - Field:      (Br, Bth, Bphi)
    - Potential:  (Ar, Ath, Aphi)

    - For cell centred: 
      Br = 1/(r sin_theta) (deltheta (sin_theta Aphi) - delphi Ath)
	 =  (1/(-r del_costheta) (sin_theta Aphi)) 
           -1/(r sin_theta) delphi Ath 
    - For Staggered: 
      Br = (sin(th+)Aphi[++0] - sin(th-)Aphi[+-0])/(r (cos(th-)-cos(th+)) ) 
          -(Ath[+0+] - Ath[-0+])/(r sin(th) dphi)  
   -------------------------------------------------------------------------*/
  br = (A3_x2p*sin(x2p)   - A3_x2m*sin(x2m))/(x1f*(cos(x2m)-cos(x2p)));

  #if DIMENSIONS == 3
  br -= 1./(x1p*sin(x2)*dx3)*(A2_x3p - A2_x3m);
  #endif 

/* -------------------------------------------------------------
   Compute Bth
    - For cell centred: 
      Bth = 1/(r) (-delr (r Aphi) +1/(sin_theta) delphi Ar)
    - For Staggered: 
      Bth = -( (r+) Aphi[++0] - (r-) Aphi[-+0])/(r dr)
            +(Ar[0++] - Ar[0+-])/(r+ sin(th) dphi)
   -------------------------------------------------------------*/
  bth = -(x1p*A3_x1p   - x1m*A3_x1m)/(x1*dx1);

  #if DIMENSIONS == 3
  bth += 1./(x1p*sin(x2)*dx3)*(A1_x3p   - A1_x3m);
  #endif 

/* -------------------------------------------------------------
    Compute Bphi
    - For cell centred: 
      Bphi = 1/(r) (delr (r Ath) - deltheta Ar)
    - For Staggered:
      Bphi = ((r+) Ath[+0+] - (r-) Ath[-0+])/(r dr)
            -(Ar[0++] - Ar[0-+])/(r dth) 
  -------------------------------------------------------------*/
  bphi = (x1p*A2_x1p - x1m*A2_x1m)/(x1*dx1) - (A1_x2p     - A1_x2m)/(x1*dx2);
  b[0] = br; b[1] = bth; b[2] = bphi;

#endif /* GEOMETRY == SPHERICAL */
  
}
#endif /* PHYSICS == MHD || PHYSICS == RMHD */

