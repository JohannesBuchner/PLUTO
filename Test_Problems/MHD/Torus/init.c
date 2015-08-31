/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Magnetized accretion torus.

  Magnetized accretion torus in 2D spherical coordinates in a 
  radially stratified atmosphere. 
  The accretion torus is determined by the equilibrium condition
  \f[
     -\frac{1}{R} + \frac{1}{2(1-a)}\frac{L_k^2}{r^{2-2a}} 
     +(n+1)\frac{p}{\rho} = const
  \f]
  where \f$p = K\rho^\gamma \f$, \c R is the spherical radius, \c r 
  is the cylindrical radius.
  Since \f$ n = 1/(\gamma-1) \f$ then \f$ 1+n = \gamma/(\gamma-1) \f$
  we can also write using \f$ a = 0 \f$ (constant angular momentum 
  distribution, Kuwabara Eq [8]):
  \f[
     -\frac{1}{R} + \frac{1}{2}\frac{L_k^2}{r^2} 
     + \frac{\gamma}{\gamma-1}K\rho^{\gamma-1} = c
     = -\frac{1}{r_{\min}} + \frac{1}{2}\frac{L_k^2}{r_{\min}^2} 
  \f]
  where the constant on the right hand side  is determined at the 
  zero pressure surface.
  The previous equation is used to compute \c K at the point where
  the density is maximum and then again to obtain the density 
  distribution:
  \f[
     K = \frac{\gamma-1}{\gamma}\left[
          c + \frac{1}{r_{\max}} - \frac{1}{2}\frac{L_k^2}{r_{\max}^2}
         \right]\frac{1}{\rho_{\max}^{\gamma-1}}\,;\qquad
     \rho =  \left[\frac{\gamma-1}{K\gamma}\left(c + \frac{1}{R} 
           - \frac{1}{2}\frac{L_k^2}{r^2}\right)\right]^{1/(\gamma-1)}
  \f]
  Thus, specifying \f$ r_{\min} \f$, \f$ r_{\max} \f$ and \f$ \rho_{\max}\f$
  determine the torus structure completely.
  The torus azimuthal velocity is compute from \f$v_\phi = L_k/r\f$.

  The atmosphere is assumed to be isothermal and it is given by the
  condition of hydrostatic balance:
  \f[
    \rho_a(R) = \eta\rho_{\max}
                \exp\left[\left(\frac{1}{R}-\frac{1}{r_{\min}}\right)
                          \frac{GM}{a^2}\right]
  \f]
  where \f$\eta\f$ is the density contrast between the (maximum) torus
  density and the atmosphere density at \f$(r_{\min},0)\f$.

  The torus surface is defined by the condition \f$p_t > p_a\f$.

  The dimensionless form of the equation employs on the following units:

  - \f$ [\rho] = \rho_{\max}\f$: reference density;
  - \f$ [L]    = R_{*} \f$: star radius;
  - \f$ [v^2] = GM/L \f$: reference (squared) velocity.

  The user defined parameters of this problem, as they appear 
  in pluto.ini, allows for control in the Torus shape and the 
  contrast of physical values: 

  -# <tt>g_inputParam[RMIN]</tt>: minimum cylindrical radius for the torus
                                  (inner rim)
  -# <tt>g_inputParam[RMAX]</tt>: radius of the Torus where pressure is maximum
  -# <tt>g_inputParam[RHO_CUT]</tt>:  minimum density to define the last
                                      contour of the magnetic vec. pot. 
  -# <tt>g_inputParam[BETA]</tt>: plasma beta 
  -# <tt>g_inputParam[ETA]</tt>   density contrast between atmosphere and Torus
  -# <tt>g_inputParam[SCALE_HEIGHT]</tt>:  atmosphere scale height \f$ H=GM/a^2\f$
 
  The magnetic field can be specified to be inside the torus or as a
  large-scale dipole field.
  In the first case, we give the vector potential:
  \f[
    A_\phi = B_0(\rho_t - \rho_{\rm cut})
  \f]
  In the second case, a dipole field is specified in the function
  ::DipoleField().
  The user-defined constant \c USE_DIPOLE can be used to select between the
  two configurations.
  
  \author A. Mignone (mignone@ph.unito.it) \n
          T. Matsakos (titos@oddjob.uchicago.edu)
  \date   Aug 16, 2012

  \b References
     - "Global Magnetohydrodynamical Simulations of Accretion Tori", 
       Hawley, J.F., ApJ (2000) 528 462.
     - "The dynamical stability of differentially rotating discs with
        constant specific angular momentum"
        Papaloizou, J.C.B, Pringle, J.E, MNRAS (1984) 208 721
     - "The Acceleration Mechanism of Resistive Magnetohydrodynamic Jets
        Launched from Accretion Disks." 
        Kuwabara, T., Shibata, K., Kudoh, T., Matsumoto, R. ApJ (2005) 621 921   
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static void DipoleField(double x1, double x2, double x3, 
                        double *Bx1, double *Bx2, double *A);

static double kappa;

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{ 
  double rmin, rmax, beta;
  double r, R, z, A_phi, gmmr, H;
  double a,c, l_k, l_k2;
  double Bo, Vo;
  double prt, pra, phi;
  double rhot, rhoa, rhocut;
  
  gmmr = g_gamma/(g_gamma - 1.0);
  #if GEOMETRY == CYLINDRICAL
   r  = x1;   /* cylindrical radius */
   z  = x2;   
   R  = sqrt(x1*x1 + x2*x2);           /* spherical radius */
  #elif GEOMETRY == SPHERICAL
   r  = x1*sin(x2);   /* cylindrical radius */
   z  = x1*cos(x2);   
   R  = x1;           /* spherical radius */
  #endif

  phi  = 0.0; 
  
  rmax    = g_inputParam[RMAX];
  rmin    = g_inputParam[RMIN];
  beta    = g_inputParam[BETA];
  rhocut  = g_inputParam[RHO_CUT];  
  H       = g_inputParam[SCALE_HEIGHT];  

  l_k     = sqrt(rmax);
  l_k2    = l_k*l_k;
  Vo      = l_k;

/* ----------------------------------------
   Solve for kappa and c assuming pra=0
   ---------------------------------------- */  
 
  c     =    - 1.0/rmin + 0.5*l_k2/(rmin*rmin);
  kappa = (c + 1.0/rmax - 0.5*l_k2/(rmax*rmax))/gmmr;

/* ----------------------------
   Torus density + pressure
   ---------------------------- */  

  a    = c + 1.0/R - 0.5*l_k2/(r*r);
  a    = MAX(a, 0.0);
  rhot = pow(a/(gmmr*kappa), 1.0/(g_gamma-1));
  prt  = kappa*pow(rhot,g_gamma);

/* -------------------------------------
    Atmospheric density + pressure 
    (force equilibrium Fg and \nabla P)
   ------------------------------------- */

  rhoa = g_inputParam[ETA]*exp((1./R - 1./rmin)/H); 
  pra  = rhoa*H;

  Bo = sqrt(2.0*kappa/beta);

  if (prt > pra && r > 2.0) {    /* Torus */
    v[RHO] = rhot;
    v[PRS] = prt;
    v[VX1] = 0.0;
    v[VX2] = 0.0;
    v[VX3] = Vo/r;
    v[TRC] = 1.0;
  }else{                          /* Atmosphere */
    v[RHO] = rhoa;
    v[PRS] = pra;
    v[VX1] = 0.0;
    v[VX2] = 0.0;
    v[VX3] = 0.0;
    v[TRC] = 0.0;
  }
  
  #if PHYSICS == MHD || PHYSICS == RMHD
   v[BX1] = 0.0;
   v[BX2] = 0.0;
   v[BX3] = 0.0;

   A_phi  = 0.0;
   #if USE_DIPOLE == NO
    if (rhot > rhocut && r > 2.0) A_phi = Bo*(rhot - rhocut);
   #endif
   
   v[AX1] = 0.0;
   v[AX2] = 0.0;
   v[AX3] = A_phi;
  #endif

  #if (USE_DIPOLE == YES) && (BACKGROUND_FIELD == NO)
   DipoleField (x1,x2,x3,v+BX1, v+BX2, v+AX3);
  #endif

  g_smallPressure = 1.e-8;
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
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background
 * magnetic field.
 *
 *********************************************************************** */
{
   double A;

   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
   #if USE_DIPOLE == YES
    DipoleField (x1,x2,x3,B0, B0+1, &A);
   #endif
}
#endif
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *  At the inner boundary we use outflow conditions, except for
 *  velocity which we reset to zero when there's an inflow
 *
 *********************************************************************** */
{
  int     i, j, k, nv;
  double  *x1, *x2, *x3;
  static double vR1[NVAR];
  double R;

  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;
 
  if (side == X1_BEG){
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        for (nv = 0; nv < NVAR; nv++){
          d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IBEG];
        }
        d->Vc[VX1][k][j][i] = MIN(d->Vc[VX1][k][j][i],0.0); 
      }
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
//       BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IBEG];
      #endif
    }
  }
}

#if (BODY_FORCE & VECTOR)
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *********************************************************************** */
{
  double  R;

  #if GEOMETRY == CYLINDRICAL
   R = sqrt(x1*x1 + x2*x2);
   g[IDIR] = -1.0/(R*R*R)*x1;
   g[JDIR] = -1.0/(R*R*R)*x2;
   g[KDIR] =  0.0; 
  #elif GEOMETRY == SPHERICAL
   R = x1;
   g[IDIR] = -1.0/(R*R);
   g[JDIR] =  0.0;
   g[KDIR] =  0.0; 
  #endif
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif

/* ********************************************************************* */
void DipoleField(double x1, double x2, double x3, 
                 double *Bx1, double *Bx2, double *A)
/*!
 *  Assign the polodial component of magnetic fields and vector
 *  potential for a dipole.
 *
 *********************************************************************** */
{
  double R, z, r, r2, r3, Bmag;
  double a;

  Bmag = sqrt(2.0*kappa/g_inputParam[BETA]);

  #if GEOMETRY == CYLINDRICAL
   R  = x1; z = x2; r = sqrt(x1*x1 + x2*x2);
   r2 = r*r;
   r3 = r2*r;
   *Bx1 = 3.0*Bmag*z*R/(r2*r3);
   *Bx2 =    -Bmag*(R*R - 3.0*z*z)/(r2*r3);
   *A   = Bmag*R/r3;
  #elif GEOMETRY == SPHERICAL
   r  = x1;
   r2 = r*r;
   r3 = r2*r;
   *Bx1 = 2.0*Bmag*cos(x2)/r3;
   *Bx2 =     Bmag*sin(x2)/r3;
   *A   = Bmag*sin(x2)/r2;
  #endif
}
