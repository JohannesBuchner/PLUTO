/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute the right hand side of the primitive MHD equations.
  
  Implements the right hand side of the quasi-linear form of the MHD 
  equations. 
  In 1D this may be written as
  \f[ 
      \partial_t{\mathbf{V}} = - A\cdot\partial_x\mathbf{V} + \mathbf{S}
  \f]
  where \f$ A \f$ is the matrix of the primitive form of the equations,
  \f$ S \f$ is the source term.

  \b Reference
 
  - "A solution adaptive upwind scheme for ideal MHD",
    Powell et al., JCP (1999) 154, 284.
 
  The function PrimRHS() implements the first term while PrimSource() 
  implements the source term part.
 
  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 26, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* *********************************************************************** */
void PrimRHS (double *v, double *dv, double cs2, double h, double *Adv)
/*!
 * Compute the matrix-vector multiplication \f$ A(\mathbf{v})\cdot 
 * d\mathbf{v} \f$ where A is the matrix of the quasi-linear form 
 * of the MHD equations.
 *
 *  \b References
 *
 *  - "A solution adaptive upwind scheme for ideal MHD",
 *    Powell et al., JCP (1999) 154, 284
 *
 *  - "An unsplit Godunov method for ideal MHD via constrained transport"
 *    Gardiner \& Stone, JCP (2005) 205, 509
 *
 * \param [in]  v    vector of primitive variables
 * \param [in]  dv   limited (linear) slopes
 * \param [in]  cs2  local sound speed
 * \param [in]  h    local enthalpy
 * \param [out] AdV  matrix-vector product
 *
 * \return This function has no return value.
 ************************************************************************* */
{
  int nv;
  double tau, scrh;
  double ch2;

  tau = 1.0/v[RHO];

/* ---------------------------------------------
         Adv[k]  Contains A[k][*]*dv[*]  
   ---------------------------------------------  */

  Adv[RHO] = v[VXn]*dv[RHO] + v[RHO]*dv[VXn];
  scrh = EXPAND(0.0, + v[BXt]*dv[BXt], + v[BXb]*dv[BXb]);

  #if EOS == IDEAL || EOS == PVTE_LAW
   Adv[VXn] = v[VXn]*dv[VXn] + tau*(dv[PRS] + scrh);
  #elif EOS == ISOTHERMAL
   Adv[VXn] = v[VXn]*dv[VXn] + tau*(cs2*dv[RHO] + scrh);
  #else 
   print ("! PRIM_RHS: not defined for this EoS\n");
   QUIT_PLUTO(1);
  #endif

  EXPAND(                                         ;    ,
         Adv[VXt] = v[VXn]*dv[VXt] - tau*v[BXn]*dv[BXt];    ,
         Adv[VXb] = v[VXn]*dv[VXb] - tau*v[BXn]*dv[BXb]; ) 

  #if DIVB_CONTROL == EIGHT_WAVES
   Adv[BXn] = v[VXn]*dv[BXn];
  #elif DIVB_CONTROL == DIV_CLEANING 
   ch2 = glm_ch*glm_ch;
   Adv[BXn]      = dv[PSI_GLM];             
   Adv[PSI_GLM] = dv[BXn]*ch2; 
  #else
   Adv[BXn] = 0.0;
  #endif
   
/* ------------------------------------------------------------------- */
/*! \note 
    In the 7-wave and 8-wave formulations we use the the same matrix 
    being decomposed into right and left eigenvectors during 
    the Characteristic Tracing step.
    Note, however, that it DOES NOT include two additional  terms 
    (-vy*dV[BX] for By, -vz*dv[BX] for Bz) that are needed in the 
    7-wave form and are added using source terms.
   --------------------------------------------------------------------- */

  EXPAND(                                                    ;          ,
         Adv[BXt] = v[BXt]*dv[VXn] - v[BXn]*dv[VXt] + v[VXn]*dv[BXt];   ,
         Adv[BXb] = v[BXb]*dv[VXn] - v[BXn]*dv[VXb] + v[VXn]*dv[BXb];)

  #if EOS == IDEAL
   Adv[PRS] = g_gamma*v[PRS]*dv[VXn] + v[VXn]*dv[PRS];
  #elif EOS == PVTE_LAW
   Adv[PRS] = cs2*v[RHO]*dv[VXn] + v[VXn]*dv[PRS];
  #endif

/*  -------------------------------------------------------------
                         Now Define Tracers
    ------------------------------------------------------------- */

#if NSCL > 0
  NSCL_LOOP(nv) Adv[nv] = v[VXn]*dv[nv];
#endif

}

/* ********************************************************************* */
void PrimSource (const State_1D *state, int beg, int end, double *a2, 
                  double *h, double **src, Grid *grid)
/*!
 * Compute source terms of the MHD equations in primitive variables.
 * These include:
 *
 *  - Geometrical sources;
 *  - Shearing-box terms 
 *  - Gravity;
 *  - terms related to divergence of B control (Powell eight wave and GLM);
 *  - FARGO source terms.
 *
 *  The rationale for choosing during which sweep a particular source 
 *  term has to be incorporated should match the same criterion used 
 *  during the conservative update. 
 *  For instance, in polar or cylindrical coordinates, curvilinear source
 *  terms are included during the radial sweep only.
 * 
 * \param [in]  state pointer to a State_1D structure
 * \param [in]  beg   initial index of computation
 * \param [in]  end   final   index of computation
 * \param [in]  a2    array of sound speed
 * \param [in]  h     array of enthalpies (not needed in MHD)
 * \param [out] src   array of source terms
 * \param [in]  grid  pointer to a Grid structure
 *
 * \note This function does not work in spherical coordinates yet. 
 *   For future implementations we annotate hereafter the induction 
 *   equation in spherical coordinates:
 *
 *  \f[ \partial_tB_r + \frac{1}{r}\partial_\theta E_\phi
 *    - \frac{1}{r\sin\theta}\partial_\phi E_\theta = -E_\phi\cot\theta/r \f]
 *  \f[ \partial_t B_\theta + \frac{1}{r\sin\theta}\partial_\phi E_r
 *    - \partial_rE_\phi =   E_\phi/r \f]
 *  \f[ \partial_t B_\phi + \partial_r E_\theta 
 *    - \frac{1}{r}\partial_\theta E_r = - E_\theta/r\f]
 *
 * where 
 *  \f[ E_\phi   = -(v \times B)_\phi   = - (v_r B_\theta - v_\theta B_r) 
 *     \,,\qquad
 *      E_\theta = -(v \times B)_\theta = - (v_\phi B_r    - v_r B_\phi) \f]
 *********************************************************************** */
{
  int    nv, i, j, k;
  double tau, dA_dV, th;
  double hscale; /* scale factor */
  double *v, *vp,  *A, *dV, r_inv, ct;
  double *x1,  *x2,  *x3;
  double *x1p, *x2p, *x3p;
  double *dx1, *dx2, *dx3;
  static double *phi_p;
  double g[3], ch2, db, scrh;

  #if ROTATING_FRAME == YES
   print1 ("! PrimSource: does not work with rotations\n");
   QUIT_PLUTO(1);
  #endif

/* ----------------------------------------------------------
   1. Memory allocation and pointer shortcuts 
   ---------------------------------------------------------- */

  if (phi_p == NULL) phi_p = ARRAY_1D(NMAX_POINT, double);

  #if GEOMETRY == CYLINDRICAL
   x1 = grid[IDIR].xgc; x1p = grid[IDIR].xr; dx1 = grid[IDIR].dx;
   x2 = grid[JDIR].xgc; x2p = grid[JDIR].xr; dx2 = grid[JDIR].dx;
   x3 = grid[KDIR].xgc; x3p = grid[KDIR].xr; dx3 = grid[KDIR].dx;
  #else  
   x1 = grid[IDIR].x; x1p = grid[IDIR].xr; dx1 = grid[IDIR].dx;
   x2 = grid[JDIR].x; x2p = grid[JDIR].xr; dx2 = grid[JDIR].dx;
   x3 = grid[KDIR].x; x3p = grid[KDIR].xr; dx3 = grid[KDIR].dx;
  #endif
  
  #ifdef GLM_MHD
   ch2 = glm_ch*glm_ch;
  #endif

  A  = grid[g_dir].A;
  dV = grid[g_dir].dV;
  hscale  = 1.0;

  i = g_i; j = g_j; k = g_k;
  
/* ----------------------------------------------------------
     initialize all elements of src to zero
   ---------------------------------------------------------- */

  memset((void *)src[0], '\0',NMAX_POINT*NVAR*sizeof(double));
  
/* ----------------------------------------------------------
   2. Compute geometrical source terms
   ---------------------------------------------------------- */

#if GEOMETRY == CYLINDRICAL

  if (g_dir == IDIR) {
    for (i = beg; i <= end; i++){
      v = state->v[i]; 

      tau   = 1.0/v[RHO];
      dA_dV = 1.0/x1[i];

      src[i][RHO] = -v[RHO]*v[VXn]*dA_dV;
      EXPAND(                                                                  ,
             src[i][iBZ]   = (v[iBR]*v[iVZ] - v[iVR]*v[iBZ])*dA_dV;            ,
             src[i][iVR]   = (v[iVPHI]*v[iVPHI] - v[iBPHI]*v[iBPHI]*tau)*dA_dV; 
             src[i][iVPHI] = (-v[iVR]*v[iVPHI] + v[iBR]*v[iBPHI]*tau)*dA_dV;)
   
      #if EOS == IDEAL
       src[i][PRS] = a2[i]*src[i][RHO];
      #endif

      #ifdef GLM_MHD
       src[i][PSI_GLM] = -v[iBR]*dA_dV*ch2;
      #endif
    }
  }

#elif GEOMETRY == POLAR

  if (g_dir == IDIR) {
    for (i = beg; i <= end; i++){
      v = state->v[i]; 

      tau   = 1.0/v[RHO];
      dA_dV = 1.0/x1[i];
      src[i][RHO]  = -v[RHO]*v[VXn]*dA_dV;

      EXPAND(                                                                ,
         src[i][iVR]   = (v[iVPHI]*v[iVPHI] - v[iBPHI]*v[iBPHI]*tau)*dA_dV;  
         src[i][iVPHI] = (-v[iVR]*v[iVPHI]  + v[iBR]*v[iBPHI]*tau)*dA_dV;    ,
         src[i][iBZ]   = ( v[iBR]*v[iVZ]    - v[iVR]*v[iBZ])*dA_dV;)

      #if EOS == IDEAL
       src[i][PRS] = a2[i]*src[i][RHO];
      #endif
      #ifdef GLM_MHD
       src[i][PSI_GLM] = -v[iBR]*dA_dV*ch2;
      #endif

    }
  }

#elif GEOMETRY == SPHERICAL 

  print1 ("! PrimSource: not implemented in Spherical geometry\n");
  QUIT_PLUTO(1);
  
#endif

/* ----------------------------------------------------------
   3.  Add body forces. This includes:
       - Coriolis terms for the shearing box module
       - Body forces
   ---------------------------------------------------------- */

#ifdef SHEARINGBOX

  if (g_dir == IDIR){
    for (i = beg; i <= end; i++) {
      src[i][VX1] =  2.0*state->v[i][VX2]*SB_OMEGA;
    } 
  }else if (g_dir == JDIR){
    for (j = beg; j <= end; j++) {
    #ifdef FARGO 
      src[j][VX2] = (SB_Q - 2.0)*state->v[j][VX1]*SB_OMEGA;
    #else
      src[j][VX2] = -2.0*state->v[j][VX1]*SB_OMEGA;
    #endif
    }
  }

#endif
  

#if (BODY_FORCE != NO)
  if (g_dir == IDIR) {

    i = beg-1;
    j = g_j;
    k = g_k;
  #if BODY_FORCE & POTENTIAL
    phi_p[i] = BodyForcePotential(x1p[i], x2[j], x3[k]);
  #endif
    for (i = beg; i <= end; i++){
    #if BODY_FORCE & VECTOR
      v = state->v[i];
      BodyForceVector(v, g, x1[i], x2[j], x3[k]);
      src[i][VX1] += g[IDIR];
    #endif
    #if BODY_FORCE & POTENTIAL
      phi_p[i]     = BodyForcePotential(x1p[i], x2[j], x3[k]); 
      src[i][VX1] -= (phi_p[i] - phi_p[i-1])/(hscale*dx1[i]);
    #endif

    /* -- Add tangential components in 1D -- */
    
    #if DIMENSIONS == 1
      EXPAND(                         , 
             src[i][VX2] += g[JDIR];  ,
             src[i][VX3] += g[KDIR];)
    #endif
    }

  }else if (g_dir == JDIR){

    i = g_i;
    j = beg - 1;
    k = g_k;
  #if BODY_FORCE & POTENTIAL
    phi_p[j] = BodyForcePotential(x1[i], x2p[j], x3[k]);
  #endif
  #if GEOMETRY == POLAR || GEOMETRY == SPHERICAL
    hscale = x1[i];
  #endif
    for (j = beg; j <= end; j++){
    #if BODY_FORCE & VECTOR
      v = state->v[j];
      BodyForceVector(v, g, x1[i], x2[j], x3[k]);
      src[j][VX2] += g[JDIR];
    #endif
    #if BODY_FORCE & POTENTIAL
      phi_p[j]     = BodyForcePotential(x1[i], x2p[j], x3[k]);
      src[j][VX2] -= (phi_p[j] - phi_p[j-1])/(hscale*dx2[j]);
    #endif

    /* -- Add 3rd component in 2D -- */

    #if DIMENSIONS == 2 && COMPONENTS == 3
      src[j][VX3] += g[KDIR];
    #endif
    }

  }else if (g_dir == KDIR){

    i = g_i;
    j = g_j;
    k = beg - 1;
  #if BODY_FORCE & POTENTIAL
    phi_p[k] = BodyForcePotential(x1[i], x2[j], x3p[k]);
  #endif
  #if GEOMETRY == SPHERICAL
    th     = x2[j];
    hscale = x1[i]*sin(th);
  #endif
    for (k = beg; k <= end; k++){
    #if BODY_FORCE & VECTOR
      v = state->v[k];
      BodyForceVector(v, g, x1[i], x2[j], x3[k]);
      src[k][VX3] += g[KDIR];
    #endif
    #if BODY_FORCE & POTENTIAL
      phi_p[k]     = BodyForcePotential(x1[i], x2[j], x3p[k]); 
      src[k][VX3] -= (phi_p[k] - phi_p[k-1])/(hscale*dx3[k]);
    #endif
    }
  }
#endif

/* -----------------------------------------------------------
   4. MHD, div.B related source terms
   ----------------------------------------------------------- */

  #ifdef GLM_MHD
   #if GEOMETRY == POLAR || GEOMETRY == SPHERICAL
    print1 ("! Error: Div. Cleaning does not work in this configuration.\n");
    print1 ("!        Try RK integrator instead\n");
    QUIT_PLUTO(1);
   #endif
   for (i = beg; i <= end; i++){
     v = state->v[i];

     tau = 1.0/v[RHO]; 
     db  = 0.5*(  A[i]  *(state->v[i+1][BXn] + state->v[i][BXn]) 
                - A[i-1]*(state->v[i-1][BXn] + state->v[i][BXn]))/dV[i];
     #if GLM_EXTENDED == NO
      EXPAND(src[i][VXn] += v[BXn]*tau*db;  ,
             src[i][VXt] += v[BXt]*tau*db;  ,
             src[i][VXb] += v[BXb]*tau*db;)
     #endif
     EXPAND(                        ,
            src[i][BXt] += v[VXt]*db; ,
            src[i][BXb] += v[VXb]*db;)
     
     #if EOS == IDEAL
      scrh = EXPAND(v[VXn]*v[BXn], + v[VXt]*v[BXt], + v[VXb]*v[BXb]);
      src[i][PRS] += (1.0 - g_gamma)*scrh*db;
      #if GLM_EXTENDED == NO
       scrh = 0.5*(state->v[i+1][PSI_GLM] - state->v[i-1][PSI_GLM])/grid[g_dir].dx[i];
       src[i][PRS] += (g_gamma - 1.0)*v[BXn]*scrh;
      #endif        
     #endif
   }
  #endif

/* ---------------------------------------------------------------
   5. Add FARGO source terms (except for SHEARINGBOX).
      When SHEARINGBOX is also enabled, we do not include
      these source terms since they're all provided by body_force)
   --------------------------------------------------------------- */
  
  #if (defined FARGO && !defined SHEARINGBOX)
   #if GEOMETRY == POLAR || GEOMETRY == SPHERICAL
    print1 ("! Time Stepping works only in Cartesian or cylindrical coords\n");
    print1 ("! Use RK instead\n");
    QUIT_PLUTO(1);
   #endif

   double **wA, *dx, *dz;
   wA = FARGO_GetVelocity();
   if (g_dir == IDIR){
     k  = g_k;
     dx = grid[IDIR].dx;
     for (i = beg; i <= end; i++){
       v = state->v[i];
       src[i][VX2] -= 0.5*v[VX1]*(wA[k][i+1] - wA[k][i-1])/dx[i];
     }
   }else if (g_dir == KDIR){
     i  = g_i;
     dz = grid[KDIR].dx;
     for (k = beg; k <= end; k++){
       v = state->v[k];
       src[k][VX2] -= 0.5*v[VX3]*(wA[k+1][i] - wA[k-1][i])/dz[k];
     }
   }
  #endif
}
