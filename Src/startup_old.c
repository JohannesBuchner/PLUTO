/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Loop on the computational cells to assign initial conditions.

  This function is called anyway, even if restart from file is enabled.
  This is useful to initialized a number of global variables and/or 
  user-defined parameters by calling Init().

  \author A. Mignone (mignone@ph.unito.it)
          B. Vaidya
  \date   May 11, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static void VectorPotentialDiffOld (double *b, int i, int j, int k, Grid *grid);

/* ********************************************************************* */
void Startup (Data *d, Grid *grid)
/*! 
 *
 *
 *
 *
 *********************************************************************** */
{
  int i, j, k;
  int isub, jsub, ksub, nsub = 5;
  int nv,  l_convert;
  static double **ucons, **uprim;
  double x1,  x2,  x3;
  double x1p, x2p, x3p;
  double x1s, x2s, x3s;
  double dx1, dx2, dx3;
  double us[256], u_av[256], b[3]; 
  double scrh;

/* ------------------------------------------------------
   0. Initialize function by allocating memory
      and setting labels.
   ------------------------------------------------------ */

  print ("> Assigning initial conditions (Startup) ...\n");

  if (uprim == NULL){
    uprim = ARRAY_2D(NMAX_POINT, NVAR, double);
    ucons = ARRAY_2D(NMAX_POINT, NVAR, double); 
  }

  EXPAND(MXn = VXn = VX1;  ,
         MXt = VXt = VX2;  ,
         MXb = VXb = VX3;)

  #if PHYSICS == MHD || PHYSICS == RMHD
   EXPAND(BXn = BX1;  ,
          BXt = BX2;  ,
          BXb = BX3;)
  #endif

/* ------------------------------------------------------
   1. Assign initial conditions point by point using
      primitive variables.
   ------------------------------------------------------ */

  KTOT_LOOP(k) { 
  JTOT_LOOP(j) { 
  ITOT_LOOP(i) { 

  /* -- Define coordinates -- */

    #if GEOMETRY == CYLINDRICAL
    x1 = grid->xgc[IDIR][i]; x1p = grid->xr[IDIR][i]; dx1 = grid->dx[IDIR][i];
    x2 = grid->xgc[JDIR][j]; x2p = grid->xr[JDIR][j]; dx2 = grid->dx[JDIR][j];
    x3 = grid->xgc[KDIR][k]; x3p = grid->xr[KDIR][k]; dx3 = grid->dx[KDIR][k];
    #else
    x1 = grid->x[IDIR][i]; x1p = grid->xr[IDIR][i]; dx1 = grid->dx[IDIR][i];
    x2 = grid->x[JDIR][j]; x2p = grid->xr[JDIR][j]; dx2 = grid->dx[JDIR][j];
    x3 = grid->x[KDIR][k]; x3p = grid->xr[KDIR][k]; dx3 = grid->dx[KDIR][k];
    #endif
    
  /* -- Reset arrays -- */

    NVAR_LOOP(nv) d->Vc[nv][k][j][i] = u_av[nv] = 0.0;

    #ifdef GLM_MHD
    u_av[PSI_GLM] = us[PSI_GLM] = 0.0;
    #ifdef PHI_GLM
    u_av[PHI_GLM] = us[PHI_GLM] = 0.0;
    #endif
    #endif

    #if INITIAL_SMOOTHING == YES

    for (ksub = 0; ksub < nsub; ksub++){ 
    for (jsub = 0; jsub < nsub; jsub++){ 
    for (isub = 0; isub < nsub; isub++){ 
            
      x1s = x1 + (double)(1.0 - nsub + 2.0*isub)/(double)(2.0*nsub)*dx1;
      x2s = x2 + (double)(1.0 - nsub + 2.0*jsub)/(double)(2.0*nsub)*dx2;
      x3s = x3 + (double)(1.0 - nsub + 2.0*ksub)/(double)(2.0*nsub)*dx3;

      Init (us, x1s, x2s, x3s);
      NVAR_LOOP(nv) u_av[nv] += us[nv]/(double)(nsub*nsub*nsub);
    }}}

    #else

    Init (u_av, x1, x2, x3);
  
    #endif

    for (nv = NVAR; nv--;  ) d->Vc[nv][k][j][i] = u_av[nv];

  /* -- Initialize cell-centered vector potential 
        (only for output purposes)                -- */

    #if PHYSICS == MHD || PHYSICS == RMHD
    #if UPDATE_VECTOR_POTENTIAL == YES
    D_EXPAND(                              ,
             d->Ax3[k][j][i] = u_av[AX3];  ,
             d->Ax1[k][j][i] = u_av[AX1];  
             d->Ax2[k][j][i] = u_av[AX2];)
    #endif
    #endif

 /* -------------------------------------------------------
     Assign staggered components.
     If a vector potential is used
     (ASSIGN_VECTOR_POTENTIAL == YES), use the
     VectorPotentialDiff() function.
     Otherwise assign staggered components directly from 
     the init.c and ignore the vector potential.

     NOTE: in N dimensions only N components are assigned 
             through this call.
    ------------------------------------------------------- */    

    #if PHYSICS == MHD || PHYSICS == RMHD
    #if ASSIGN_VECTOR_POTENTIAL == YES
    VectorPotentialDiffOld(b, i, j, k, grid);

    #ifdef STAGGERED_MHD
    for (nv = 0; nv < DIMENSIONS; nv++) d->Vs[nv][k][j][i] = b[nv];
    #else
    for (nv = 0; nv < DIMENSIONS; nv++) d->Vc[BX1+nv][k][j][i] = b[nv];
    #endif

    #else

    #ifdef STAGGERED_MHD
    D_EXPAND(
             Init (u_av, x1p, x2, x3);
             d->Vs[BX1s][k][j][i] = u_av[BX1];       ,

             Init (u_av, x1, x2p, x3);
             d->Vs[BX2s][k][j][i] = u_av[BX2];       ,

             Init (u_av, x1, x2, x3p);
             d->Vs[BX3s][k][j][i] = u_av[BX3];)
    #endif
    #endif  /* ASSIGN_VECTOR_POTENTIAL */
    #endif /* PHYSICS == MHD || PHYSICS == RMHD */

  }}}

/* --------------------------------------------------------
   2. Call Init_Domain() to assign primitive variables
      by looping over computational cells.
      This is new in PLUTO 4.3 and provides an
      alternative to the pointwise initialization.
   -------------------------------------------------------- */

  InitDomain(d, grid);

  #if FORCED_TURB == YES
    struct ForcedTurb *Ft;
    Ft = d->Ft;
    ForcedTurb_Init(Ft);
    print("> Initialized Forced Turbulence %d Modes\n",Ft->NModes);
    ForcedTurb_OUNoiseInit(Ft->OUPhases, 6*(Ft->NModes), Ft->OUVar);
  #endif

/* --------------------------------------------------------
   3. Compute cell-centered magnetic field by simple
      arithmetic averaging.
      Useful only for saving the first output, since the
      average will be re-computed anyway at the beginning
      of the computation.
   -------------------------------------------------------- */

#ifdef STAGGERED_MHD
  DOM_LOOP(k,j,i){
    D_EXPAND( 
      d->Vc[BX1][k][j][i] = 0.5*(d->Vs[BX1s][k][j][i] + d->Vs[BX1s][k][j][i-1]); , 
      d->Vc[BX2][k][j][i] = 0.5*(d->Vs[BX2s][k][j][i] + d->Vs[BX2s][k][j-1][i]); ,
      d->Vc[BX3][k][j][i] = 0.5*(d->Vs[BX3s][k][j][i] + d->Vs[BX3s][k-1][j][i]);
    )
  }
#endif

/* ------------------------------------------------------
   4.  Check if values have physical meaning.
   ------------------------------------------------------ */

#if PHYSICS != ADVECTION
  KDOM_LOOP(k){ x3 = grid->x[KDIR][k];
  JDOM_LOOP(j){ x2 = grid->x[JDIR][j];
  IDOM_LOOP(i){ x1 = grid->x[IDIR][i];

    for (nv = NVAR; nv--;  ) us[nv] = d->Vc[nv][k][j][i];

    if (us[RHO] <= 0.0) {
      print ("! Startup(): negative density, zone [%f, %f, %f]\n", x1,x2,x3);
      QUIT_PLUTO(1);
    }
    #if HAVE_ENERGY
     if (us[PRS] <= 0.0) {
       print ("! Startup(): negative pressure, zone [%f, %f, %f]\n",x1,x2,x3);
       QUIT_PLUTO(1);
     }
    #endif

    #if (PHYSICS == RHD || PHYSICS == RMHD) 
     scrh = EXPAND(us[VX1]*us[VX1], + us[VX2]*us[VX2], + us[VX3]*us[VX3]);
     if (scrh >= 1.0){
       print ("! Startup(): total velocity exceeds 1\n"); 
       QUIT_PLUTO(1);
     }
    #endif
  }}}
#endif

/* ------------------------------------------------------
   5. Set boundary conditions.
   ------------------------------------------------------ */

  Boundary (d, -1, grid);
}
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
void VectorPotentialDiffOld (double *b, int i, int j, int k, Grid *grid)
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

  x1 = grid->x[IDIR][i]; 
  x2 = grid->x[JDIR][j]; 
  x3 = grid->x[KDIR][k]; 

  dx1 = grid->dx[IDIR][i];
  dx2 = grid->dx[JDIR][j];
  dx3 = grid->dx[KDIR][k];

  #ifdef STAGGERED_MHD   
   x1p = grid->xr[IDIR][i]; x1m = grid->xl[IDIR][i];
   x2p = grid->xr[JDIR][j]; x2m = grid->xl[JDIR][j];
   x3p = grid->xr[KDIR][k]; x3m = grid->xl[KDIR][k];
   x1f = grid->xr[IDIR][i]; /* for staggered MHD, we compute magnetic */
   x2f = grid->xr[JDIR][j]; /* field at face centers                  */
   x3f = grid->xr[KDIR][k];
   
  #else
   x1p = grid->x[IDIR][i] + dx1; x1m = grid->x[IDIR][i] - dx1;
   x2p = grid->x[JDIR][j] + dx2; x2m = grid->x[JDIR][j] - dx2;
   x3p = grid->x[KDIR][k] + dx3; x3m = grid->x[KDIR][k] - dx3;
   x1f = grid->x[IDIR][i]; /* for cell-centered MHD, we compute magnetic */
   x2f = grid->x[JDIR][j]; /* field at cell-centers                      */
   x3f = grid->x[KDIR][k];
   
   dx1 = 2.*grid->dx[IDIR][i]; /* redefine the spacing between cell-centers */
   dx2 = 2.*grid->dx[JDIR][j];
   dx3 = 2.*grid->dx[KDIR][k];
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
