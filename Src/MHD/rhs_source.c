/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Add source terms to the right hand side of HD/MHD eqns.

  Add source terms to the right hand side of hydro/MHD equations.
  These include

  -# Shearing box source terms;
  -# Geometrical source term in curvilinear coordinates;
     (these include centrifugal-like terms as well as conservative
      discretizations to achieve angular momentum and B-phi conservation)
  -# FARGO or ROTATION source terms (included to enforce conservation of
     total energy or angular mometum);
  -# Body forces;
  -# Powell's 8-waves source terms;
  -# Extended GLM source terms.
  
  \note
  <B> Shearing-Box source terms: </B>
      With the shearing-box module, the user is required to provide only
      gravity while Coriolis source terms are automatically added
      here. 
      Without FARGO, the source terms of the momentum and energy equations are
      \f{eqnarray*}{ 
        \vec{S}_{\vec{m}} &=&
           \Big[ 2\Omega_0^2 q x \rho + 2\Omega_0\rho v_y\Big]\hvec{i}
         + \Big[-2\Omega_0\rho v_x\Big]\hvec{j}
         + \Big[-\rho\Omega_0^2z\Big]\hvec{k} \\ 
         S_E &=& \rho \vec{v} \cdot\vec{g}
         = \rho v_x (2\Omega_0^2 qx) + \rho v_z (-\Omega_0^2 z)
      \f}
      When the FARGO module is used, they modify to 
      \f{eqnarray*}{ 
       \vec{S}_{\vec{m}} &=&
             \Big[ 2\Omega_0\rho v'_y\Big]\hvec{i}
           + \Big[(q-2)\Omega_0\rho v_x\Big]\hvec{j}
           + \Big[-\rho\Omega_0^2z\Big]\hvec{k} \\ \noalign{\medskip}
         S_E &=& (\rho v'_yv_x - B_yB_x)q\Omega_0 + \rho v_z(-\Omega_0^2 z)
      \f}
      The energy source term is included during the Fargo transport step,
      see FARGO_Source(). 

  \b Reference
    - "A conservative orbital advection scheme for simulations 
       of magnetized shear flows with the PLUTO code"\n
       Mignone et al, A&A (2012) 545, A152 (--> Appendix A.1.1)

  \author A. Mignone (mignone@ph.unito.it)\
          B. Vaidya 
  \date   May 13, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if PHYSICS == MHD
  #if BACKGROUND_FIELD == YES
    #define TotBB(v, b0, a, b) (v[a]*(v[b] + b0[b-BX1]) + b0[a-BX1]*v[b])
  #else
    #define TotBB(v, b0, a, b) (v[a]*v[b])
  #endif
#else
  #define TotBB(v, b0, a, b)  0.0
#endif

#ifndef iMPHI
 #define iMPHI MX2  /* -- for Cartesian coordinates -- */
#endif

/* *********************************************************************** */
void RightHandSideSource (const Sweep *sweep, timeStep *Dts,
                          int beg, int end, double dt, double *phi_p, Grid *grid)
/*! 
 *
 * \param [in,out]  state  pointer to State_1D structure
 * \param [in]      Dts    pointer to time step structure
 * \param [in]      beg    initial index of computation
 * \param [in]      end    final   index of computation
 * \param [in]      dt     time increment
 * \param [in]      phi_p  force potential at interfaces
 * \param [in]      grid   pointer to Grid structure
 *
 * \return This function has no return value.
 ************************************************************************* */
{
  int    i, j, k, nv;

  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);

  double r_1, g[3];

  double dtdx, scrh, ct;
  double Sm;
  double *x1   = grid->x[IDIR],  *x2  = grid->x[JDIR],  *x3  = grid->x[KDIR];
  double *x1p  = grid->xr[IDIR], *x2p = grid->xr[JDIR], *x3p = grid->xr[KDIR];
  double *x1m  = grid->xl[IDIR], *x2m = grid->xl[JDIR], *x3m = grid->xl[KDIR];
  double *dx1  = grid->dx[IDIR], *dx2 = grid->dx[JDIR], *dx3 = grid->dx[KDIR];
#if GEOMETRY == SPHERICAL
  double *rt   = grid->rt;
  double *sp   = grid->sp;
  double *sm   = grid->sp-1;
  double *s    = grid->s;
  double *dmu  = grid->dmu;
#endif
  double ***dV = grid->dV;
  double **rhs  = sweep->rhs;
  double **flux = sweep->flux;
  double **vp   = stateL->v;
  double **vm   = stateR->v-1;
  double *p;
  double cl;
  double **Bg0, **wA, w, wp, vphi, phi_c;
  double vc[NVAR], *vg;

#ifdef FARGO
  wA = FARGO_GetVelocity();
#endif

#if (PHYSICS == MHD) && (BACKGROUND_FIELD == YES)
  GetBackgroundField (stateC, beg, end, CELL_CENTER, grid);
  Bg0 = stateC->Bbck;
#endif

/* --------------------------
      pointer shortcuts
   -------------------------- */

  p  = sweep->press;
  
  i = g_i;  /* will be redefined during x1-sweep */
  j = g_j;  /* will be redefined during x2-sweep */
  k = g_k;  /* will be redefined during x3-sweep */
  
  if (g_dir == IDIR){

    for (i = beg; i <= end; i++) {
      dtdx = dt/dx1[i];

    /* --------------------------------------------
       I1. Add geometrical source term
       -------------------------------------------- */

#if GEOMETRY == CARTESIAN

      vg = stateC->v[i];
      #ifdef SHEARINGBOX  /* Include Coriolis term for x-momentum */
      rhs[i][MX1] += dt*2.0*vg[RHO]*vg[VX2]*SB_OMEGA;
      #endif

      #if (defined FARGO && !defined SHEARINGBOX) 
      w = wA[k][i];
      #endif
  
#elif GEOMETRY == CYLINDRICAL

      r_1 = 1.0/x1[i];
      vg  = vc;
      NFLX_LOOP(nv) vc[nv] = 0.5*(vp[i][nv] + vm[i][nv]);       
      #if COMPONENTS == 3
      vphi = vc[iVPHI];
      IF_ROTATING_FRAME(w     = g_OmegaZ*x1[i];
                        vphi += w;)
      rhs[i][MX1]   += dt*(vc[RHO]*vphi*vphi - TotBB(vc, Bg0[i], iBPHI, iBPHI))*r_1;
      #endif  /* COMPONENTS == 3 */

#elif GEOMETRY == POLAR

      r_1 = 1.0/x1[i];
      vg  = vc;
      NFLX_LOOP(nv) vc[nv] = 0.5*(vp[i][nv] + vm[i][nv]); 
      vphi = vc[iVPHI];
      #if (defined FARGO) || (ROTATING_FRAME == YES)
      w = 0.0;
      IF_FARGO         (w += wA[k][i];)
      IF_ROTATING_FRAME(w += g_OmegaZ*x1[i];)
      vphi += w;
      #endif
      rhs[i][MX1] += dt*(  vc[RHO]*vphi*vphi 
                         - TotBB(vc, Bg0[i], iBPHI, iBPHI))*r_1;

#elif GEOMETRY == SPHERICAL

      r_1 = 1.0/x1[i];
      vg  = vc;
      NFLX_LOOP(nv) vc[nv] = 0.5*(vp[i][nv] + vm[i][nv]);
      vphi = SELECT(0.0, 0.0, vc[iVPHI]);
      #if (defined FARGO) || (ROTATING_FRAME == YES)
      w = 0.0;
      IF_FARGO         (w += wA[j][i];)
      IF_ROTATING_FRAME(w += g_OmegaZ*x1[i]*s[j];)
      vphi += w;
      #endif
      Sm = vc[RHO]*(EXPAND(  0.0, + vc[VX2]*vc[VX2], + vphi*vphi));
      Sm += EXPAND(  0.0, - TotBB(vc, Bg0[i], iBTH, iBTH), 
                          - TotBB(vc, Bg0[i], iBPHI,iBPHI));
      rhs[i][MX1] += dt*Sm*r_1;

#endif

    /* ----------------------------------------------------
       I2. Modify rhs to enforce conservation
       ---------------------------------------------------- */

#if ((defined FARGO && !defined SHEARINGBOX) || (ROTATING_FRAME == YES)) 
      rhs[i][iMPHI] -= w*rhs[i][RHO];
      IF_ENERGY(rhs[i][ENG] -= w*(rhs[i][iMPHI] + 0.5*w*rhs[i][RHO]);)
#endif

    /* ----------------------------------------------------
       I3. Include body forces
       ---------------------------------------------------- */

#if (BODY_FORCE & VECTOR)
      BodyForceVector(vg, g, x1[i], x2[j], x3[k]);
      rhs[i][MX1]   += dt*vg[RHO]*g[IDIR];
      IF_ENERGY(rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[IDIR];)
      
      /* -- Add tangential components in 1D -- */
      
      #if DIMENSIONS == 1   
      EXPAND(                                    ,
             rhs[i][MX2] += dt*vg[RHO]*g[JDIR];  ,
             rhs[i][MX3] += dt*vg[RHO]*g[KDIR];)
      #if HAVE_ENERGY
      rhs[i][ENG] += dt*(EXPAND(0.0, + vg[RHO]*vg[VX2]*g[JDIR],
                                     + vg[RHO]*vg[VX3]*g[KDIR] );
      #endif                         
      #endif
#endif
      
#if (BODY_FORCE & POTENTIAL)
      rhs[i][MX1]   -= dtdx*vg[RHO]*(phi_p[i] - phi_p[i-1]);
      IF_ENERGY(phi_c       = BodyForcePotential(x1[i], x2[j], x3[k]); 
                rhs[i][ENG] -= phi_c*rhs[i][RHO];)
#endif

    }
    
  } else if (g_dir == JDIR){

    scrh = dt;
#if GEOMETRY == POLAR
    scrh /= x1[i];
    r_1   = 1.0/x1[i];
#elif GEOMETRY == SPHERICAL
    scrh /= rt[i];
    r_1   = 1.0/rt[i];
#endif
    for (j = beg; j <= end; j++) {
      dtdx = scrh/dx2[j];

    /* --------------------------------------------
       J1. Add geometrical source terms
       -------------------------------------------- */

#if GEOMETRY != SPHERICAL

      vg = stateC->v[j];

#ifdef SHEARINGBOX /* Include Coriolis term for y-momentum */
    #ifdef FARGO 
      rhs[j][MX2] += dt*(SB_Q-2.0)*vg[RHO]*vg[VX1]*SB_OMEGA;
    #else
      rhs[j][MX2] -= dt*2.0*vg[RHO]*vg[VX1]*SB_OMEGA;
    #endif
#endif

#elif GEOMETRY == SPHERICAL

      ct = 1.0/tan(x2[j]); 
      vg = vc;
      NFLX_LOOP(nv) vc[nv] = stateC->v[j][nv];
      vphi = SELECT(0.0, 0.0, vc[iVPHI]);
      #if (defined FARGO) || (ROTATING_FRAME == YES)
      w = 0.0; 
      IF_FARGO         (w += wA[j][i];)
//      IF_ROTATING_FRAME(w += g_OmegaZ*x1[i]*sin(x2[j]);)
      IF_ROTATING_FRAME(w += g_OmegaZ*x1[i]*s[j];)
      vphi += w;
      #endif
      Sm = vc[RHO]*(EXPAND(  0.0, - vc[iVTH]*vc[iVR], + ct*vphi*vphi));
      Sm += EXPAND(0.0, +    TotBB(vc, Bg0[j], iBTH, iBR), 
                        - ct*TotBB(vc, Bg0[j], iBPHI, iBPHI));
      rhs[j][MX2] += dt*Sm*r_1;

    /* ----------------------------------------------------
       J2. Modify rhs to enforce conservation
       ---------------------------------------------------- */

  #if (defined FARGO) || (ROTATING_FRAME == YES)
      rhs[j][MX3] -= w*rhs[j][RHO];
      IF_ENERGY(rhs[j][ENG]  -= w*(rhs[j][MX3] + 0.5*w*rhs[j][RHO]);)
  #endif

#endif  /* GEOMETRY == SPHERICAL */

    /* ----------------------------------------------------
       J3. Include Body force
       ---------------------------------------------------- */

#if (BODY_FORCE & VECTOR)
      BodyForceVector(vg, g, x1[i], x2[j], x3[k]);
      rhs[j][MX2]   += dt*vg[RHO]*g[JDIR];
      IF_ENERGY(rhs[j][ENG] += dt*0.5*(flux[j][RHO] + flux[j-1][RHO])*g[JDIR];)

      /* -- Add tangential components in 2D -- */
      
      #if DIMENSIONS == 2 && COMPONENTS == 3  
      rhs[j][MX3] += dt*vg[RHO]*g[KDIR];
      IF_ENERGY(rhs[j][ENG] += dt*vg[RHO]*vg[VX3]*g[KDIR];)
      #endif
#endif

#if (BODY_FORCE & POTENTIAL)
      rhs[j][MX2]   -= dtdx*vg[RHO]*(phi_p[j] - phi_p[j-1]);
      IF_ENERGY(phi_c        = BodyForcePotential(x1[i], x2[j], x3[k]); 
                rhs[j][ENG] -= phi_c*rhs[j][RHO];)
#endif
    }

  }else if (g_dir == KDIR){

    scrh  = dt;
#if GEOMETRY == SPHERICAL
    scrh *= dx2[j]/(rt[i]*dmu[j]);
#endif

    for (k = beg; k <= end; k++) {
      dtdx = scrh/dx3[k];
      vg   = stateC->v[k];

#if (GEOMETRY == CARTESIAN) || (GEOMETRY == POLAR)
    /* ------------------------------------------------------
       K2. modify rhs to enforce conservation (FARGO only)
           (solid body rotations are not included the
            velocity depends on the cylindrical radius only)
       ------------------------------------------------------ */

  #if (defined FARGO && !defined SHEARINGBOX) 
       w = wA[k][i];
       rhs[k][MX2] -= w*rhs[k][RHO];
       IF_ENERGY(rhs[k][ENG] -= w*(rhs[k][MX2] + 0.5*w*rhs[k][RHO]);)
  #endif
#endif

    /* ----------------------------------------------------
       K3. Include body forces
       ---------------------------------------------------- */

#if (BODY_FORCE & VECTOR)
      BodyForceVector(vg, g, x1[i], x2[j], x3[k]);
      rhs[k][MX3]   += dt*vg[RHO]*g[KDIR];
      IF_ENERGY(rhs[k][ENG] += dt*0.5*(flux[k][RHO] + flux[k-1][RHO])*g[KDIR];)
#endif

#if (BODY_FORCE & POTENTIAL)
      rhs[k][MX3]   -= dtdx*vg[RHO]*(phi_p[k] - phi_p[k-1]);
      IF_ENERGY(phi_c        = BodyForcePotential(x1[i], x2[j], x3[k]);
                rhs[k][ENG] -= phi_c*rhs[k][RHO];)
#endif
    }
  }

/* --------------------------------------------------
              Powell's source terms
   -------------------------------------------------- */

#if (PHYSICS == MHD) && (DIVB_CONTROL == EIGHT_WAVES)
  for (i = beg; i <= end; i++) {
    EXPAND(rhs[i][MX1] += dt*sweep->src[i][MX1];  ,
           rhs[i][MX2] += dt*sweep->src[i][MX2];  ,
           rhs[i][MX3] += dt*sweep->src[i][MX3];)

    EXPAND(rhs[i][BX1] += dt*sweep->src[i][BX1];  ,
           rhs[i][BX2] += dt*sweep->src[i][BX2];  ,
           rhs[i][BX3] += dt*sweep->src[i][BX3];)
  #if HAVE_ENERGY
    rhs[i][ENG] += dt*sweep->src[i][ENG];
  #endif
  }
#endif

/* -------------------------------------------------
            Extended GLM source terms
   ------------------------------------------------- */

#if (defined GLM_MHD) && (GLM_EXTENDED == YES)
  GLM_ExtendedSource (sweep, dt, beg, end, grid);
#endif

}
#undef TotBB
