/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute the resistive MHD flux for cell-centered field.

  Compute the resistive fluxes for the induction and energy equations
  to update cell-centered fields.
  Cell-centered fields are updated also in the CT formalism.

  In the induction equation, fluxes are computed by explicitly writing
  the curl operator in components. In Cartesian components, for instance,
  one has
  \f[
     \pd{\vec{B}}{t} = - \nabla\times{\vec{E}_{\rm res}} =
      \pd{}{x}\left(\begin{array}{c}
       0           \\ \noalign{\medskip}
       \eta_z J_z  \\ \noalign{\medskip}
     - \eta_y J_y \end{array}\right)
       +
      \pd{}{y}\left(\begin{array}{c}
     - \eta_z J_z  \\ \noalign{\medskip}
           0       \\ \noalign{\medskip}
       \eta_x J_x \end{array}\right)
        +
      \pd{}{z}\left(\begin{array}{c}
        \eta_y J_y \\ \noalign{\medskip}
       -\eta_x J_x \\ \noalign{\medskip}
            0    \end{array}\right)      
       \,,\qquad
     \left(\vec{E}_{\rm res} = \tens{\eta}\cdot\vec{J}\right)
  \f]
  where \f$\tens{\eta}\f$ is the resistive diagonal tensor and 
  \f$\vec{J} = \nabla\times\vec{B}\f$ is the current density.
  The corresponding contribution to the energy equation is
  \f[
    \pd{E}{t} = -\nabla\cdot\Big[(\tens{\eta}\cdot\vec{J})\times\vec{B}\Big]
              = - \pd{}{x}\left(\eta_yJ_yB_z - \eta_zJ_zB_y\right)
                - \pd{}{y}\left(\eta_zJ_zB_x - \eta_xJ_xB_z\right)
                - \pd{}{z}\left(\eta_xJ_xB_y - \eta_yJ_yB_x\right)
  \f] 
  
  The sign of the flux is given by writing each equation as
  \f[
     \pd{B_i}{t} = \nabla\cdot\vec{F}_i
  \f]
  (flux is positive on the right hand side).
  
  \b References
     - "The PLUTO Code for Adaptive Mesh Computations in Astrophysical 
        Fluid Dynamics" \n
       Mignone et al, ApJS (2012) 198, 7M
       
  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos
  \date    Sept 7, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void ResistiveFlux (Data_Arr V, Data_Arr curlB, double **res_flux,
                    double **dcoeff, int beg, int end, Grid *grid)
/*! 
 *
 * \param [in]   V      3D data array of primitive variables
 * \param [out]  curlB  the current
 * \param [out]  dcoeff the diffusion coefficients evaluated at 
 *                      cell interfaces
 * \param [in]   beg    initial index of computation
 * \param [in]   end    final   index of computation
 * \param [in]  grid    pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int  i, j, k, nv;
  double *x1, *x2, *x3, scrh;
  double eta[3], vp[NVAR], Jp[NVAR];
  Data_Arr etas;  
  double ***Jx1, ***Jx2, ***Jx3;
  double ***eta_x1, ***eta_x2, ***eta_x3;

  etas = GetStaggeredEta();
  Jx1 = curlB[IDIR]; eta_x1 = etas[IDIR];
  Jx2 = curlB[JDIR]; eta_x2 = etas[JDIR];
  Jx3 = curlB[KDIR]; eta_x3 = etas[KDIR];

  x1 = grid->x[IDIR];  
  x2 = grid->x[JDIR];  
  x3 = grid->x[KDIR]; 

/* -------------------------------------------------------- 
   1.  Compute resistive flux for the induction and energy
       equations at X1, X2 or X3 faces.
   -------------------------------------------------------- */

  if (g_dir == IDIR){
  
    j  = g_j; k = g_k;    
    x1 = grid->xr[IDIR];
    for (i = beg; i <= end; i++){

    /* -- interface values at (i+1/2, j, k) -- */

      for (nv = 0; nv < NVAR; nv++){
        vp[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j][i+1]);
      }

    /* ----------------------------------------------------
        Compute transverse components of J at X1-faces as 
        arithmetic averages of the staggered arrays.
       ---------------------------------------------------- */

      EXPAND(                                 ;   ,
             Jp[KDIR] = AVERAGE_Y(Jx3,k,j-1,i);   ,
             Jp[JDIR] = AVERAGE_Z(Jx2,k-1,j,i);)

      EXPAND(                                     ;   ,
             eta[KDIR] = AVERAGE_Y(eta_x3,k,j-1,i);   ,
             eta[JDIR] = AVERAGE_Z(eta_x2,k-1,j,i);)

      EXPAND(res_flux[i][BX1] =   0.0;                  ,
             res_flux[i][BX2] =   eta[KDIR]*Jp[KDIR];   ,
             res_flux[i][BX3] = - eta[JDIR]*Jp[JDIR];)
      #if HAVE_ENERGY
      res_flux[i][ENG] = EXPAND(0.0, + vp[BX2]*res_flux[i][BX2], 
                                     + vp[BX3]*res_flux[i][BX3]);
      #endif

    /* -- store diffusion coefficient -- */

      EXPAND(dcoeff[0][i] = 0.0;        ,
             dcoeff[1][i] = eta[KDIR];  ,
             dcoeff[2][i] = eta[JDIR];)
//print ("J = %12.6e  %12.6e  %12.6e\n",Jx1[k][j][i], Jx2[k][j][i], Jx3[k][j][i]);
//ShowVector(Jp,3);
//ShowVector(eta,3);
//ShowVector(res_flux[i]+BX1,3);
    }

  }else if (g_dir == JDIR){
  
    i  = g_i; k = g_k;    
    x2 = grid->xr[JDIR];
    for (j = beg; j <= end; j++){

    /* -- interface values at (i, j+1/2, k) -- */

      for (nv = 0; nv < NVAR; nv++) {
        vp[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j+1][i]);
      }

    /* ----------------------------------------------------
        Compute transverse components of J at X2-faces as 
        arithmetic averages of the staggered arrays.
       ---------------------------------------------------- */

      EXPAND(Jp[KDIR] = AVERAGE_X(Jx3,k,j,i-1);   ,
                                              ;   ,
             Jp[IDIR] = AVERAGE_Z(Jx1,k-1,j,i);)

      EXPAND(eta[KDIR] = AVERAGE_X(eta_x3,k,j,i-1);   ,
                                                  ;   ,
             eta[IDIR] = AVERAGE_Z(eta_x1,k-1,j,i);)

      EXPAND(res_flux[j][BX1] = - eta[KDIR]*Jp[KDIR];    ,
             res_flux[j][BX2] =   0.0;                   ,
             res_flux[j][BX3] =   eta[IDIR]*Jp[IDIR];)

      #if HAVE_ENERGY
      res_flux[j][ENG] = EXPAND(  vp[BX1]*res_flux[j][BX1],
                                + 0.0, 
                                + vp[BX3]*res_flux[j][BX3]);
      #endif

     /* -- store diffusion coefficient -- */

      EXPAND(dcoeff[0][j] = eta[KDIR];   ,
             dcoeff[1][j] = 0.0;         ,
             dcoeff[2][j] = eta[IDIR];)
    }

  }else if (g_dir == KDIR){
  
    i = g_i; j = g_j;    
    x3 = grid->xr[KDIR];
    for (k = beg; k <= end; k++){

    /* -- interface values at (i, j, k+1/2) -- */

      for (nv = 0; nv < NVAR; nv++) {
        vp[nv] = 0.5*(V[nv][k][j][i] + V[nv][k+1][j][i]);
      }

    /* ----------------------------------------------------
        Compute transverse components of J at X3-faces as 
        arithmetic averages of the staggered arrays.
       ---------------------------------------------------- */

      Jp[JDIR] = AVERAGE_X(Jx2,k,j,i-1);
      Jp[IDIR] = AVERAGE_Y(Jx1,k,j-1,i);

      eta[JDIR] = AVERAGE_X(eta_x2,k,j,i-1);
      eta[IDIR] = AVERAGE_Y(eta_x1,k,j-1,i);

      res_flux[k][BX1] =   eta[JDIR]*Jp[JDIR];
      res_flux[k][BX2] = - eta[IDIR]*Jp[IDIR];
      res_flux[k][BX3] =   0.0;

      #if HAVE_ENERGY
      res_flux[k][ENG] = vp[BX1]*res_flux[k][BX1] + vp[BX2]*res_flux[k][BX2];
      #endif

    /* -- store diffusion coefficient -- */

      EXPAND(dcoeff[0][k] = eta[JDIR];  ,
             dcoeff[1][k] = eta[IDIR];  ,
             dcoeff[2][k] = 0.0;)
    }
  }
}
