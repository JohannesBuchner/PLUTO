/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file 
  \brief MHD Shock Cloud interaction.
    
  Set the initial condition for the 2D or 3D shock-cloud problem.
  The shock propagates in the x-direction and it is initially located
  at \c x=0.6 with post- and -pre shock value equal to 
  \f[
    \left(\begin{array}{l}
      \rho \\ \noalign{\medskip}
       v_x  \\ \noalign{\medskip}
       p   \\ \noalign{\medskip}
       B_y  \\ \noalign{\medskip}
       B_z  
    \end{array}\right) 
     =
    \left(\begin{array}{l}
        3.86859      \\ \noalign{\medskip}
          0          \\ \noalign{\medskip}
        167.345      \\ \noalign{\medskip}
        B_{\rm post} \\ \noalign{\medskip}
       -B_{\rm post} \\ \noalign{\medskip}
    \end{array}\right)
     \quad\mathrm{for}\quad x < 0.6\,,\qquad\qquad
    \left(\begin{array}{l}
      \rho \\ \noalign{\medskip}
       v_x  \\ \noalign{\medskip}
       p   \\ \noalign{\medskip}
       B_y  \\ \noalign{\medskip}
       B_z  
    \end{array}\right) 
     =
    \left(\begin{array}{l}
          1         \\ \noalign{\medskip}
        -11.2536    \\ \noalign{\medskip}
         1          \\ \noalign{\medskip}
        B_{\rm pre} \\ \noalign{\medskip}
       -B_{\rm pre} \\ \noalign{\medskip}
    \end{array}\right)
     \quad\mathrm{for}\quad x > 0.6\,,\qquad
  \f]
  while the remaining vector components are 0.
  The cloud has radius \f$R\f$ centered at
  \f$r_0 = (x_0,y_0,z_0) = (0.8, 0.5, 0.5)\f$ and larger density 
  \f$\rho = 10\f$.

  The interaction may be divided into two phases: 1) the collapse stage
  where the front of the cloud is strongly compressed and two fast shocks
  are generated and 2) the reexpansion phase which begins when the 
  transmitted fast shock overtake the back of the cloud.
  
  The runtime parameters that are read from \c pluto.ini are 
  - <tt>g_inputParam[B_POST]</tt>: sets the post-shock magnetic field;
  - <tt>g_inputParam[B_PRE]</tt>:  sets the pre-shock magnetic field;
  - <tt>g_inputParam[RADIUS]</tt>: sets the radius of the cloud;

  A list of the available configurations is given in the following table:

  <CENTER>
  Conf.| GEOMETRY  |DIM| T. STEPPING|RECONSTRUCTION|divB|  AMR
  -----|-----------|---|------------|--------------| ---|-------
   #01 |CARTESIAN  | 2 |  RK2       |  LINEAR      | CT | NO
   #02 |CARTESIAN  | 2 |  ChTr      |  LINEAR      | CT | NO    
   #03 |CARTESIAN  | 2 |  HANCOCK   |  LINEAR      | 8W | NO    
   #04 |CARTESIAN  | 3 |  ChTr      |  LINEAR      | CT | NO
   #05 |CARTESIAN  | 3 |  RK3       |  LINEAR      | CT | NO    
   #06 |CARTESIAN  | 3 | ChTr       |  LINEAR      | GLM| NO    
   #07 |CARTESIAN  | 3 |  RK3       |  WENO3_FD    | GLM| NO    
   #08 |CARTESIAN  | 3 | HANCOK     |  LINEAR      | GLM| YES
   #09 |CARTESIAN  | 2 |  RK2       |  LINEAR      | CT | NO 
  </CENTER>

  \image html mhd_shock_cloud.02.jpg "Density map with overplotted field lines for configuration #02"

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 17, 2014

  \b References:
     - "Interactions between magnetohydrodynamical shocks and denser clouds"
        W. Dai an P.R. Woodward, ApJ (1994) 436, 776.
     - "The divB = 0 Constraint in Shock-Capturing MHD codes"
        G. Toth, JCP 161, 605-652 (2000)
     - "A central constrained transport scheme for ideal magnetohydrodynamics" 
        U. Ziegler, JCP 196 (2004), 393
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x, double y, double z)
/*
 *
 *********************************************************************** */
{
 double r, x0, y0, z0;

  x0 = 0.8;
  y0 = 0.5;
  z0 = 0.5;

  if (x < 0.6) {
    v[RHO] = 3.86859;
    v[VX1] = 0.0;
    v[VX2] = 0.0;
    v[VX3] = 0.0;
    v[BX1] = 0.0;
    v[BX2] =  g_inputParam[B_POST];
    v[BX3] = -g_inputParam[B_POST];
    v[PRS] = 167.345;
  }else{
    v[RHO] = 1.0;
    v[VX1] = -11.2536;
    v[VX2] = 0.0;
    v[VX3] = 0.0;
    v[BX1] = 0.0;
    v[BX2] = g_inputParam[B_PRE];
    v[BX3] = g_inputParam[B_PRE];
    v[PRS] = 1.0;
  } 

  /*  ----  CLOUD  ----  */

  r = D_EXPAND((x-x0)*(x-x0), + (y-y0)*(y-y0) , + (z-z0)*(z-z0));

  if (sqrt(r) < g_inputParam[RADIUS]) v[RHO] = 10.0;

/* no need for potential vector, 
   since CT_VEC_POT_INIT is set to NO

  v[AX1] = x3*v[BX2] - x2*v[BX3];
  v[AX2] = 0.0;
  v[AX3] = 0.0;
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

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *********************************************************************** */
{
  int   i, j, k;

  if (side == X1_END){          /* -- select the boundary side -- */
    if (box->vpos == CENTER){   /* -- select the variable position -- */
      BOX_LOOP(box,k,j,i){      /* -- Loop over boundary zones -- */
        d->Vc[RHO][k][j][i] = 1.0;
        EXPAND(d->Vc[VX1][k][j][i] = -11.2536;  ,
               d->Vc[VX2][k][j][i] = 0.0;       ,
               d->Vc[VX3][k][j][i] = 0.0;)
        d->Vc[PRS][k][j][i] = 1.0;
        EXPAND(d->Vc[BX1][k][j][i] = 0.0;        ,
               d->Vc[BX2][k][j][i] = g_inputParam[B_PRE]; ,
               d->Vc[BX3][k][j][i] = g_inputParam[B_PRE];)
      }
    }else if (box->vpos == X2FACE){  /* -- y staggered field -- */
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = g_inputParam[B_PRE];
      #endif
    }else if (box->vpos == X3FACE){  /* -- z staggered field -- */
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX3s][k][j][i] = g_inputParam[B_PRE];
      #endif
    }
  }
}
