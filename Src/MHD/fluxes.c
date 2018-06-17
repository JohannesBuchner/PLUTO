/* ///////////////////////////////////////////////////////////////////// */
/*! 
 \file
 \brief Compute the MHD flux.                                             

  Compute the flux of the conservative MHD equations in the direction 
  given by ::g_dir.
  A useful formulary is

  - Divergence of a tensor:
    \f[
      \nabla\cdot(\vec{A}\vec{B}) =
      \pd{}{x}\left(\begin{array}{l}
        A_xB_x \\ \noalign{\medskip}
        A_xB_y \\ \noalign{\medskip}
        A_xB_z
      \end{array}\right)
      +
      \pd{}{y}\left(\begin{array}{l}
        A_yB_x \\ \noalign{\medskip}
        A_yB_y \\ \noalign{\medskip}
        A_yB_z
      \end{array}\right)
      +
      \pd{}{z}\left(\begin{array}{l}
        A_zB_x \\ \noalign{\medskip}
        A_zB_y \\ \noalign{\medskip}
        A_zB_z
      \end{array}\right)
      \qquad\Longrightarrow\qquad
      \sum_n\pd{}{x_n}\left(\begin{array}{l}
        A_nB_x \\ \noalign{\medskip}
        A_nB_y \\ \noalign{\medskip}
        A_nB_z
      \end{array}\right)
    \f]
 
  - Curl of a vector:
    \f[
      \nabla\times\vec{A} =
      \pd{}{x}\left(\begin{array}{l}
         0  \\ \noalign{\medskip}
        -A_z \\ \noalign{\medskip}
         A_y
      \end{array}\right)
      +
      \pd{}{y}\left(\begin{array}{l}
        A_z \\ \noalign{\medskip}
        0   \\ \noalign{\medskip}
        -A_x
      \end{array}\right)
      +
      \pd{}{z}\left(\begin{array}{l}
        -A_y\\ \noalign{\medskip}
         A_x \\ \noalign{\medskip}
         0
      \end{array}\right)
      \qquad\Longrightarrow\qquad
      \sum_n\pd{}{x_n} \left(\hvec{n}\times\vec{A}\right)
    \f]
 
  - Curl of external product
    \f[
       \nabla\times\left(\vec{A}\times\vec{B}\right) =
       \nabla\cdot \left(\vec{B}\vec{A} - \vec{A}\vec{B}\right) =
      \pd{}{x}\left(\begin{array}{l}
        0      \\ \noalign{\medskip}
        B_xA_y - B_yA_x \\ \noalign{\medskip}
        B_xA_z - B_zA_x
      \end{array}\right)
      +
      \pd{}{y}\left(\begin{array}{l}
        B_yA_x - B_xA_y \\ \noalign{\medskip}
        0               \\ \noalign{\medskip}
        B_yA_z - B_zA_y
      \end{array}\right)
      +
      \pd{}{z}\left(\begin{array}{l}
        B_zA_x - B_xA_z \\ \noalign{\medskip}
        B_zA_y - B_yA_z\\ \noalign{\medskip}
        0 
      \end{array}\right)
      \qquad\Longrightarrow\qquad
      \sum_n\pd{}{x_n}\left(\begin{array}{l}
        A_nB_x - A_xB_n\\ \noalign{\medskip}
        A_nB_y - A_yB_n\\ \noalign{\medskip}
        A_nB_z - A_zB_n
      \end{array}\right)
    \f]
    

  This function defines the component of the hyperbolic flux tensor 
  of the standard MHD equations.\n
  In what follows:
  - \c VXn, \c MXn, \c BXn are the velocity, momentum and magnetic field 
    components in the direction given by ::g_dir (normal, \c "n")
  - \c VXt, \c MXt, \c BXt and \c VXb, \c MXb, \c BXb are the transverse 
    components (tangent \c "t" and bi-tangent \c "b").

  Normal, tangent and bitangent indices are set in the
  SetVectorIndices() function.

 \author A. Mignone (mignone@ph.unito.it)
 \date   April 14, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Flux (const State *state, int beg, int end)
/*!
 * Compute the component of the flux tensor normal to the direction of
 * integration.
 *                
 * \param [in,out]  state   Pointer to a state structure
 * \param [in]      beg     initial index of computation 
 * \param [in]      end     final   index of computation
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int nv, i;
  double vB, ptot;
  double bt1, bt2, bt3;
  double *v, *u, *fx, *Bbck, B0B1;
  double Bmag2, cB, *cCR;
#if HALL_MHD
  double *J, vH[3], Eh[3], Bsq, ne, vHB;
  #if COMPONENTS != 3
    #error Hall_MHD requires COMPONENTS to be 3
  #endif
#endif

/* --------------------------------------------------------
   0. Set normal, tangent and bi-tanget indices
   -------------------------------------------------------- */

  int nDIR, tDIR, bDIR;
  double *Fcr, qg, st, sb;
  if (g_dir == IDIR){
    st = -1.0;
    sb = +1.0;
    
  }else if (g_dir == JDIR){
    st = +1.0;
    sb = -1.0;

  }else if (g_dir == KDIR){
    st = -1.0;
    sb = +1.0;
  }
  nDIR = VXn-VX1; tDIR = VXt-VX1; bDIR = VXb-VX1;


/* --------------------------------------------------------
   1. Compute MHD flux
   -------------------------------------------------------- */

  for (i = beg; i <= end; i++) {

    v  = state->v[i];
    u  = state->u[i];
    fx = state->flux[i];

    Bmag2 = EXPAND(v[BX1]*v[BX1] , + v[BX2]*v[BX2], + v[BX3]*v[BX3]);

    #if HAVE_ENERGY
    ptot  = v[PRS] + 0.5*Bmag2;
    #elif EOS == ISOTHERMAL
    ptot  = state->a2[i]*v[RHO] + 0.5*Bmag2;
    #else
    print ("! Flux(): not defined for this EoS\n");
    QUIT_PLUTO(1);
    #endif

  /* ------------------------------------------------------
     1a. Compute ideal MHD fluxes
     ------------------------------------------------------ */

    vB    = EXPAND(v[VX1]*v[BX1] , + v[VX2]*v[BX2], + v[VX3]*v[BX3]);

    fx[RHO] = u[MXn];
    EXPAND(fx[MX1] = v[VXn]*u[MX1] - v[BXn]*v[BX1];  ,
           fx[MX2] = v[VXn]*u[MX2] - v[BXn]*v[BX2];  ,
           fx[MX3] = v[VXn]*u[MX3] - v[BXn]*v[BX3]; ) 

    EXPAND(fx[BXn] = 0.0;                             ,
           fx[BXt] = v[VXn]*v[BXt] - v[BXn]*v[VXt];   ,
           fx[BXb] = v[VXn]*v[BXb] - v[BXn]*v[VXb]; )
    #if HAVE_ENERGY
    fx[ENG] = (u[ENG] + ptot)*v[VXn] - v[BXn]*vB;
    #endif

  /* ------------------------------------------------------
     1b. Add background field
     ------------------------------------------------------ */

    #if BACKGROUND_FIELD == YES
    Bbck  = state->Bbck[i] - BX1;
    B0B1  = EXPAND(Bbck[BX1]*v[BX1], + Bbck[BX2]*v[BX2], + Bbck[BX3]*v[BX3]);
    ptot += B0B1;

    EXPAND(fx[MX1] -= Bbck[BXn]*v[BX1] + v[BXn]*Bbck[BX1];  ,
           fx[MX2] -= Bbck[BXn]*v[BX2] + v[BXn]*Bbck[BX2];  ,
           fx[MX3] -= Bbck[BXn]*v[BX3] + v[BXn]*Bbck[BX3];)

    EXPAND(                                                 ,
           fx[BXt] += v[VXn]*Bbck[BXt] - Bbck[BXn]*v[VXt];  ,
           fx[BXb] += v[VXn]*Bbck[BXb] - Bbck[BXn]*v[VXb]; )
    #if HAVE_ENERGY
    fx[ENG] += B0B1*v[VXn] - Bbck[BXn]*vB;
    #endif
    #endif
  
    state->prs[i] = ptot;

    #ifdef GLM_MHD
    fx[BXn]     = v[PSI_GLM];
    fx[PSI_GLM] = glm_ch*glm_ch*v[BXn];
    #endif

  /* ------------------------------------------------------
     1c. Add Hall-MHD terms to the induction equation
     ------------------------------------------------------ */

    #if HALL_MHD == EXPLICIT
    J  = state->J[i];  
    ne = HallMHD_ne(v);
    vH[IDIR] = -J[IDIR]/ne;
    vH[JDIR] = -J[JDIR]/ne;
    vH[KDIR] = -J[KDIR]/ne;

    Eh[IDIR] = -(vH[JDIR]*v[BX3] - vH[KDIR]*v[BX2]);
    Eh[JDIR] = -(vH[KDIR]*v[BX1] - vH[IDIR]*v[BX3]);
    Eh[KDIR] = -(vH[IDIR]*v[BX2] - vH[JDIR]*v[BX1]);
      
    fx[BXt] += st*Eh[bDIR];
    fx[BXb] += sb*Eh[tDIR];

    #if HAVE_ENERGY
    vHB = EXPAND(vH[IDIR]*v[BX1], + vH[JDIR]*v[BX2], + vH[KDIR]*v[BX3]);
    Bsq = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
    fx[ENG] += vH[nDIR]*Bsq - vHB*v[BXn]; 
    #endif
    #endif  /* HALL_MHD */

  }

}
