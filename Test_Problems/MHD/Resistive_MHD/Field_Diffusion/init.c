/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Magnetic field diffusion in 2D and 3D.

  Sets the initial conditions for a magnetic field diffusion problem
  in 2 or 3 dimensions.
  This is a useful test to check the ability of the code to solve 
  standard diffusion problems.
  The magnetic field has initially a Gaussian profile, and an anisotropic
  resistivity is possible.
  This problem has an analytical solution given, in 2D, by
  \f[
    B_x(y,t) = \exp(-y^2/4\eta_zt)/\sqrt{t}
    \quad\quad\quad
    B_y(x,t) = \exp(-x^2/4\eta_zt)/\sqrt{t}
    \quad\quad\quad
    B_z(x,y,t) = \exp(-x^2/4\eta_yt)\exp(-y^2/4\eta_xt)/t
  \f]
  and in 3D by
  \f[ 
    \begin{array}{lcl}
    B_x(y,z,t) &=& \exp(-y^2/4\eta_zt)\exp(-z^2/4\eta_yt)/t \\ \noalign{\medskip}
    B_y(x,z,t) &=& \exp(-x^2/4\eta_zt)\exp(-z^2/4\eta_xt)/t \\ \noalign{\medskip}
    B_z(x,y,t) &=& \exp(-x^2/4\eta_yt)\exp(-y^2/4\eta_xt)/t
    \end{array}
  \f]
  The initial condition is simply set using the previous profiles 
  with \c t=1.
  In order to solve only the parabolic term in the induction equation
  (\f$ \nabla\times\vec{J}=\nabla\times(\eta\nabla\times\vec{B})\f$) 
  we give to the fluid a very large intertia using a high density 
  value.
  Moreover, to avoid any fluid motion, the velocity is reset to zero 
  at each time step by using the \c INTERNAL_BOUNDARY.

  The runtime parameters that are read from \c pluto.ini are 
  - <tt>g_inputParam[ETAX]</tt>: sets the resistivity along the \f$x\f$ direction;
  - <tt>g_inputParam[ETAY]</tt>: sets the resistivity along the \f$y\f$ direction;
  - <tt>g_inputParam[ETAZ]</tt>: sets the resistivity along the \f$z\f$ direction;

  \image html 3D-cart-0.jpg "Initial and final profiles of the numerical (points) and analytical (lines) solutions for a component of the magnetic field."

  The configurations use both \c EXPLICIT and \c STS time integrators and
  different geometries are explored:

  <CENTER>
  Conf.|GEOMETRY  |DIM|T.STEPPING|divB|RESISTIVITY
  -----|----------|---|----------|----| ----------
   #01 |CARTESIAN | 3 | RK2      | 8W | EXPLICIT
   #02 |CARTESIAN | 3 | HANCOCK  | GLM| STS
   #03 |CARTESIAN | 3 | RK3      | CT | EXPLICIT
   #04 |CARTESIAN | 3 | HANCOCK  | GLM| EXPLICIT
   #05 |POLAR     | 3 | RK2      | 8W | EXPLICIT
   #06 |POLAR     | 3 | RK2      | 8W | STS
   #07 |SPHERICAL | 3 | RK2      | 8W | EXPLICIT
   #08 |SPHERICAL | 3 | RK2      | 8W | STS
   #09 |SPHERICAL | 3 | RK3      | CT | EXPLICIT
   #10 |CARTESIAN | 2 | HANCOCK  | CT | EXPLICIT
   #11 |CARTESIAN | 3 | HANCOCK  | CT | EXPLICIT
   #12 |CARTESIAN | 2 | HANCOCK  | CT | STS
   #13 |SPHERICAL | 3 | RK2      | CT | STS
   #14 |CARTESIAN | 3 | HANCOCK  | GLM| RKL
   #13 |SPHERICAL | 3 | RK2      | GLM| RKL
  </CENTER>

  \authors A. Mignone (mignone@ph.unito.it) \n
           T. Matsakos (titos@oddjob.uchicago.edu)
  \date    March 03, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

void BoundValues (double *v, double x1, double x2, double x3, double t);

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  BoundValues(us, x1, x2, x3, 1.0);
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
/* 
 *
 *********************************************************************** */
{}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *********************************************************************** */
{
  int  i, j, k, nv;
  double *x, *y, *z;
  double *xp, *yp, *zp;
  double t, dt, vb[256];

  if (side == 0){  
    #if GEOMETRY == CARTESIAN
    TOT_LOOP(k,j,i){
      EXPAND(d->Vc[VX1][k][j][i] = 0.0;  ,
             d->Vc[VX2][k][j][i] = 0.0;  ,
             d->Vc[VX3][k][j][i] = 0.0;)
      #if ENTROPY_SWITCH == YES
       d->flag[k][j][i] |= FLAG_ENTROPY;
      #endif
    }
    #endif
  }

  t = g_time + 1.0;

  x = grid->x[IDIR]; xp = grid->xr[IDIR];
  y = grid->x[JDIR]; yp = grid->xr[JDIR];
  z = grid->x[KDIR]; zp = grid->xr[KDIR];

  if (side == X1_BEG || side == X2_BEG || side == X3_BEG ||
      side == X1_END || side == X2_END || side == X3_END){  /* -- All Boundaries -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
        BoundValues (vb, x[i], y[j], z[k], t);
        for (nv = NVAR; nv--;  ) d->Vc[nv][k][j][i] = vb[nv];
      }
    }else if (box->vpos == X1FACE){
      #ifdef STAGGERED_MHD
       if (side == X1_BEG || side == X1_END) return;
       BOX_LOOP(box,k,j,i){
         BoundValues (vb, xp[i], y[j], z[k], t);
         d->Vs[BX1s][k][j][i] = vb[BX1];
       }
      #endif
    }else if (box->vpos == X2FACE){
      if (side == X2_BEG || side == X2_END) return;
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i){
         BoundValues (vb, x[i], yp[j], z[k], t);
         d->Vs[BX2s][k][j][i] = vb[BX2];
       }
      #endif

    }else if (box->vpos == X3FACE){
      if (side == X3_BEG || side == X3_END) return;
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) {
         BoundValues (vb, x[i], y[j], zp[k], t);
         d->Vs[BX3s][k][j][i] = vb[BX3];
       }
      #endif

    }
  }
}

/* ******************************************************************* */
void BoundValues (double *v, double x1, double x2, double x3, double t)
/*
 *
 * Time-dependent exact solution
 * 
 ********************************************************************* */
{
  double Bx, By, Bz;
  double x, y, z;
  double eta_x, eta_y, eta_z, st;

  eta_x = g_inputParam[ETAX];
  eta_y = g_inputParam[ETAY];
  eta_z = g_inputParam[ETAZ];

/* -- find Cartesian coordinates from (x1,x2,x3) -- */

  #if GEOMETRY == CARTESIAN
   x = x1; y = x2; z = x3;
  #elif GEOMETRY == POLAR
   x = x1*cos(x2) - 5.0; 
   y = x1*sin(x2) - 5.0; 
   z = x3 - 5.0;
  #elif GEOMETRY == SPHERICAL
   x = x1*cos(x3)*sin(x2) - 5.0;
   y = x1*sin(x3)*sin(x2) - 5.0;
   z = x1*cos(x2) - 5.0;
  #else
   #error geometry not valid
  #endif

  v[RHO] = 1.0e9;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  #if EOS != ISOTHERMAL
   v[PRS] = 1.0;
  #endif
  #ifdef GLM_MHD
   v[PSI_GLM] = 0.0;
  #endif

/* -- find Cartesian components of magnetic field -- */
/*
  Bx = exp(-0.25*(y*y)/g_inputParam[ETAZ]/t)*exp(-0.25*(z*z)/g_inputParam[ETAY]/t)/t;
  By = exp(-0.25*(x*x)/g_inputParam[ETAZ]/t)*exp(-0.25*(z*z)/g_inputParam[ETAX]/t)/t;
  Bz = exp(-0.25*(x*x)/g_inputParam[ETAY]/t)*exp(-0.25*(y*y)/g_inputParam[ETAX]/t)/t;
*/
 
  #if DIMENSIONS == 2 
   st = sqrt(t);

   Bx = exp(-0.25*y*y/(eta_z*t))/st;
   By = exp(-0.25*x*x/(eta_z*t))/st;
   Bz = exp(-0.25*x*x/(eta_y*t))*exp(-0.25*y*y/(eta_x*t))/t;
  #elif DIMENSIONS == 3
   Bx = exp(-0.25*y*y/(eta_z*t))*exp(-0.25*z*z/(eta_y*t))/t;
   By = exp(-0.25*x*x/(eta_z*t))*exp(-0.25*z*z/(eta_x*t))/t;
   Bz = exp(-0.25*x*x/(eta_y*t))*exp(-0.25*y*y/(eta_x*t))/t;
  #endif

/* -----------------------------------------------------------
     Transform back to original coordinate system
   ----------------------------------------------------------- */
  
  #if GEOMETRY == CARTESIAN
   v[BX1] = Bx;
   v[BX2] = By;
   v[BX3] = Bz;
  #elif GEOMETRY == POLAR
   v[BX1] =  cos(x2)*Bx + sin(x2)*By;
   v[BX2] = -sin(x2)*Bx + cos(x2)*By;
   v[BX3] =  Bz;
  #elif GEOMETRY == SPHERICAL
   v[BX1] =  cos(x3)*sin(x2)*Bx + sin(x3)*sin(x2)*By + cos(x2)*Bz;
   v[BX2] =  cos(x3)*cos(x2)*Bx + sin(x3)*cos(x2)*By - sin(x2)*Bz;
   v[BX3] = -sin(x3)*Bx + cos(x3)*By;
  #else
   print1 ("! BoundValues: GEOMETRY not defined\n");
   QUIT_PLUTO(1);
  #endif  

}
