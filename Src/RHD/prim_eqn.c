/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute the right hand side of the relativistic 
         hydro (RHD) equations in primitive form.
  
  Implements the right hand side of the quasi-linear form of the 
  relativistic hydro equations. 
  In 1D this may be written as
  \f[ 
      \partial_t{\mathbf{V}} = - A\cdot\partial_x\mathbf{V} + \mathbf{S}
  \f]
  where \f$ A \f$ is the matrix of the primitive form of the equations,
  \f$ S \f$ is the source term.

  \b Reference:
    - "The Piecewise Parabolic Method for  Multidimensional Relativistic 
       Fluid Dynamics", Mignone, Plewa and Bodo, ApJS (2005) 160,199.
 
  The function PrimRHS() implements the first term while PrimSource() 
  implements the source term part.

  \author A. Mignone (mignone@ph.unito.it)
  \date   July 03, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void PrimRHS (double *w, double *dw, double cs2, double h, double *Adw)
/*!
 * Compute the matrix-vector multiplication \f$ A(\mathbf{v})\cdot 
 * \Delta\mathbf{v} \f$ where \c A is the matrix of the quasi-linear form 
 * of the RHD equations.
 *
 * \param [in]  w    vector of primitive variables;
 * \param [in]  dw   limited (linear) slopes;
 * \param [in]  cs2  local sound speed;
 * \param [in]  h    local enthalpy;
 * \param [out] Adw  matrix-vector product.
 *
 * \note 
 * \return This function has no return value.
 ************************************************************************** */
{
  int nv;
  double rho, u1, u2, u3;
  double d_2, g2, scrh;
  double g, v1, v2, v3, gx2;

#if RECONSTRUCT_4VEL
  rho = w[RHO];
  EXPAND(u1  = w[VXn];  ,
         u2  = w[VXt];  ,
         u3  = w[VXb];)

  scrh = EXPAND(u1*u1, + u2*u2, + u3*u3);
  g2   = 1.0 + scrh;
  g    = sqrt(g2);
  d_2  = g/(g2 - scrh*cs2);

/* Get 3vel  */

  EXPAND(v1 = u1/g; ,
         v2 = u2/g; ,
         v3 = u3/g;)

  gx2 = 1.0/(1.0 - v1*v1);

  scrh = EXPAND(v1*dw[VXn], + v2*dw[VXt], + v3*dw[VXb]);

  Adw[PRS] =  d_2*(rho*h*cs2*(dw[VXn] - v1*scrh)
                      + u1*(1.0 - cs2)*dw[PRS]);

  Adw[RHO] = v1*dw[RHO] - (v1*dw[PRS] - Adw[PRS])/(h*cs2);
 
  scrh = 1.0/(g*rho*h);
  d_2  = u1*dw[PRS] - g*Adw[PRS];

  EXPAND(Adw[VXn] = v1*dw[VXn] + scrh*(dw[PRS] + u1*d_2);  ,
         Adw[VXt] = v1*dw[VXt] + scrh*u2*d_2;                ,
         Adw[VXb] = v1*dw[VXb] + scrh*u3*d_2;)
  
#else
  rho = w[RHO];
  EXPAND(v1 = w[VXn];  ,
         v2 = w[VXt];  ,
         v3 = w[VXb];)

  g2  = EXPAND(v1*v1, + v2*v2, + v3*v3);
  d_2 = 1.0/(1.0 - g2*cs2);
  g2  = 1.0/(1.0 - g2);

  Adw[PRS] = d_2*(cs2*rho*h*dw[VXn] 
              + v1*(1.0 - cs2)*dw[PRS]);

  Adw[RHO] = v1*dw[RHO] - (v1*dw[PRS] - Adw[PRS])/(h*cs2);
  
  scrh = 1.0/(g2*rho*h);
  EXPAND(Adw[VXn] =  v1*dw[VXn] 
                  + scrh*(dw[PRS] - v1*Adw[PRS]); ,
         Adw[VXt] =  v1*dw[VXt] -  scrh*v2*Adw[PRS];  ,
         Adw[VXb] =  v1*dw[VXb] -  scrh*v3*Adw[PRS];)
#endif

#if NSCL > 0 
  NSCL_LOOP(nv)  Adw[nv] = v1*dw[nv];
#endif

}

/* *********************************************************************  */
void PrimSource (const State_1D *state, int beg, int end, 
                 double *a2, double *h, double **src, Grid *grid)
/*!
 * Compute source terms of the RHD equations in primitive variables.
 *
 *  - Geometrical sources;
 *  - Gravity;
 *
 *  The rationale for choosing during which sweep a particular source 
 *  term has to be incorporated should match the same criterion used 
 *  during the conservative update. 
 *  For instance, in polar or cylindrical coordinates, curvilinear source
 *  terms are included during the radial sweep only.
 * 
 * \param [in]  state pointer to a State_1D structure;
 * \param [in]  beg   initial index of computation;
 * \param [in]  end   final   index of computation;
 * \param [in]  a2    array of sound speed; 
 * \param [in]  h     array of enthalpies (not needed in MHD);
 * \param [out] src   array of source terms;
 * \param [in]  grid  pointer to a Grid structure.
 * \return This function has no return value.
 *
 *********************************************************************** */
{
  int    nv, i;
  double r_1, scrh, alpha;
  double vel2, delta;
  double *q, *x1, *x2, *x3, g[3];

  #if GEOMETRY == CYLINDRICAL
   x1 = grid[IDIR].xgc;
   x2 = grid[JDIR].xgc;
   x3 = grid[KDIR].xgc;
  #else  
   x1 = grid[IDIR].x; 
   x2 = grid[JDIR].x; 
   x3 = grid[KDIR].x; 
  #endif

  for (i = beg; i <= end; i++){
  for (nv = NVAR; nv--;  ){
    src[i][nv] = 0.0;
  }}

#if GEOMETRY == CARTESIAN

#elif GEOMETRY == CYLINDRICAL 
  
  if (g_dir == IDIR){
    for (i = beg; i <= end; i++) {
 
      r_1  = 1.0/x1[i];
      q    = state->v[i];
      vel2 = EXPAND(q[VX1]*q[VX1], +q[VX2]*q[VX2], +q[VX3]*q[VX3]);

      #if RECONSTRUCT_4VEL
      scrh    = sqrt(1.0 + vel2);
      alpha   = q[VXn]*r_1*scrh/(1.0 + vel2*(1.0 - a2[i]));
      scrh    = a2[i]*alpha;
      print1 ("! Primitive source terms not yet implemented\n");
      print1 ("! with 4-vel. Please try 3-vel\n");
      QUIT_PLUTO(1);
      #else
      alpha = q[VXn]*r_1/(1.0 - a2[i]*vel2);
      scrh  = a2[i]*(1.0 - vel2)*alpha;
      #endif

      src[i][RHO] = -q[RHO]*alpha;
      EXPAND (src[i][VX1] = scrh*q[VX1];  ,
              src[i][VX2] = scrh*q[VX2];  ,
              src[i][VX3] = scrh*q[VX3];)

      #if COMPONENTS == 3
      EXPAND(src[i][iVR]   +=  q[iVPHI]*q[iVPHI]*r_1;   ,
                                                        ,
             src[i][iVPHI] += -q[iVPHI]*q[iVR]*r_1;)
      #endif

      src[i][PRS] = -a2[i]*q[RHO]*h[i]*alpha;

    }
  }

#elif GEOMETRY == SPHERICAL && DIMENSIONS == 1 && COMPONENTS == 1

  if (g_dir == IDIR){
    for (i = beg; i <= end; i++) {
 
      r_1  = 1.0/x1[i];
      q    = state->v[i];
      vel2 = EXPAND(q[VX1]*q[VX1], +q[VX2]*q[VX2], +q[VX3]*u[VX3]);

      #if RECONSTRUCT_4VEL 
      print1 ("! PrimRHSD(): primitive source terms not yet implemented\n");
      print1 ("!            with 4-vel. Please try 3-vel\n");
      QUIT_PLUTO(1);
      #else
      delta = 1.0 - vel2*a2[i];
      scrh  = 2.0*q[iVR]/delta*r_1;
      #endif

      src[i][RHO] = -q[RHO]*scrh;
      EXPAND (src[i][VX1] = scrh*q[VX1]*a2[i]*(1.0 - vel2);  ,
              src[i][VX2] = scrh*q[VX2]*a2[i]*(1.0 - vel2);  ,
              src[i][VX3] = scrh*q[VX3]*a2[i]*(1.0 - vel2);)

      src[i][PRS] = src[i][RHO]*h[i]*a2[i];
    }
  }

#else 

  print1 ("! PrimRHS(): primitive source terms not available for this geometry\n");
  print1 ("!            Use RK integrators\n");
  QUIT_PLUTO(1);

#endif   
  
/* -----------------------------------------------------------
                   Add body force
   ----------------------------------------------------------- */

  #if (BODY_FORCE != NO)
   i = beg - 1;
   if (g_dir == IDIR) {
     #if BODY_FORCE & POTENTIAL
      gPhi[i] = BodyForcePotential(x1p[i], x2[g_j], x3[g_k]);
     #endif
     for (i = beg; i <= end; i++){
       #if BODY_FORCE & VECTOR
        v = state->v[i];
        BodyForceVector(v, g, x1[i], x2[g_j], x3[g_k]);
        src[i][VX1] += g[g_dir];
       #endif
       #if BODY_FORCE & POTENTIAL
        gPhi[i]     = BodyForcePotential(x1p[i], x2[g_j], x3[g_k]); 
        src[i][VX1] -= (gPhi[i] - gPhi[i-1])/(hscale*dx1[i]);
       #endif

       #if DIMENSIONS == 1
        EXPAND(                        ,
               src[i][VX2] += g[JDIR];  ,
               src[i][VX3] += g[KDIR];)
       #endif
     }
   }else if (g_dir == JDIR){
     #if BODY_FORCE & POTENTIAL
      gPhi[i] = BodyForcePotential(x1[g_i], x2p[i], x3[g_k]);
     #endif
     for (i = beg; i <= end; i++){
       #if BODY_FORCE & VECTOR
        v = state->v[i];
        BodyForceVector(v, g, x1[g_i], x2[i], x3[g_k]);
        src[i][VX2] += g[g_dir];
       #endif
       #if BODY_FORCE & POTENTIAL
        gPhi[i]     = BodyForcePotential(x1[g_i], x2p[i], x3[g_k]);
        src[i][VX2] -= (gPhi[i] - gPhi[i-1])/(hscale*dx2[i]);
       #endif

       #if DIMENSIONS == 2 && COMPONENTS == 3
        src[i][VX3] += g[KDIR];
       #endif
     }
   }else if (g_dir == KDIR){
     #if BODY_FORCE & POTENTIAL
      gPhi[i] = BodyForcePotential(x1[g_i], x2[g_j], x3p[i]);
     #endif
     for (i = beg; i <= end; i++){
       #if BODY_FORCE & VECTOR
        v = state->v[i];
        BodyForceVector(v, g, x1[g_i], x2[g_j], x3[i]);
        src[i][VX3] += g[g_dir];
       #endif
       #if BODY_FORCE & POTENTIAL
        gPhi[i]     = BodyForcePotential(x1[g_i], x2[g_j], x3p[i]); 
        src[i][VX3] -=  (gPhi[i] - gPhi[i-1])/(hscale*dx3[i]);
       #endif
     }
   }
  #endif
}
