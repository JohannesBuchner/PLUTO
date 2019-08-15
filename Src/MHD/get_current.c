/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Compute the curl of magnetic field.

  Compute the electric current (defined as J = curl(B)) at cell edges
  for both staggered MHD and cell-centered MHD.
  That is, the three components of J = (Jx, Jy, Jz) are always placed
  at different locations inside the cell, no matter what scheme is
  used to control the divergence condition:

  - \c Jx  at  <tt> (i, j+1/2, k+1/2) </tt>
  - \c Jy  at  <tt> (i+1/2, j, k+1/2) </tt>
  - \c Jz  at  <tt> (i+1/2, j+1/2, k) </tt>

  The same rule apply to the components of resistivity \c eta which
  are computed and stored inside this function.

  For a compact implementation, we note that the curl of a vector in 
  the three system of coordinates normally adopted may be written as 
  \f[
     \begin{array}{lcl}
      \left(\nabla\times\vec{B}\right)_{x_1} &=&  \DS
        \frac{1}{d_{23}}\pd{}{x_2}(a_{23}B_{x_3})
      - \frac{1}{d_{32}}\pd{B_{x_2}}{x_3}
   \\ \noalign{\bigskip}
      \left(\nabla\times\vec{B}\right)_{x_2} &=&  \DS
        \frac{1}{d_{31}}\pd{B_{x_1}}{x_3}
      - \frac{1}{d_{13}}\pd{}{x_1}(a_{13}B_{x_3})
   \\ \noalign{\bigskip}
      \left(\nabla\times\vec{B}\right)_{x_3} &=&  \DS
        \frac{1}{d_{12}}\pd{}{x_1}(a_{12}B_{x_2})
      - \frac{1}{d_{21}}\pd{B_{x_1}}{x_2}
     \end{array}
  \f]
  where the coefficients \f$d_{nm}=1\f$ and \f$a_{nm}=1\f$ except:
  \f[
     \begin{array}{lll}
            d_{nm} = 1; & 
     \quad  a_{nm} = 1; & 
     \qquad ({\rm Cartesian})
     \\ \noalign{\medskip}
           d_{13} = r\,,\qquad & 
           d_{32} = d_{31} = \infty;  
     \quad a_{13} = r;  & 
     \qquad ({\rm Cylindrical})
     \\ \noalign{\medskip}
           d_{23} = d_{12} = d_{21} = r; & 
     \quad a_{12} = r;  & 
     \qquad ({\rm Polar})
     \\ \noalign{\medskip}
       d_{23} = d_{32} = d_{31} = r\sin\theta \,,\qquad
       d_{13} = d_{12} = d_{21} = r; &  
     \quad  a_{23} = \sin\theta\,,\qquad a_{13} = a_{12} = r; &
     \qquad ({\rm Spherical})
     \end{array}
  \f]
  In the actual implementation we use \f$d_{nm} \to 1/(d_{nm}\Delta x_n)\f$.

  \note In cylindrical (or polar) geometry, there's a potential
        singular behavior with the definition of \f$ d_{13} (d_{12})
        \quad\mathrm{at}\quad r=0\f$.
        An equivalent formulation can be obtained by using absolute
        values in the geometrical coefficient and thus replacing
        \f[
           d_{13} = r_{i+\HALF}\Delta r \quad\to\quad
                    \frac{r_{i+1}|r_{i+1}| - r_i |r_i|}{2}\,,\qquad

           a_{13} = r_i \quad\to\quad  |r_i|.                
        \f]
        This is obviously the desired derivative away from the axis 
        (where r = |r|) and it gives the correct behavior at the 
        axis, where Bphi --> -Bphi.          
        Bphi = 0  --> Bphi \sim alpha*r --> d(r Bphi)/(r dr) = 2*alpha 
           

  \attention Do NOT use this function to compute terms like (curl V)^2 !!

  \authors A. Mignone (mignone@ph.unito.it)
           
  \date    Sep 07, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void GetCurrent (const Data *d, Grid *grid)
/*!
 *
 * Compute the current inside \c d->J.
 *
 * \param [in,out]  d       pointer to the PLUTO data structure
 * \param [in]      grid    pointer to an array of Grid structures
 *   
 *********************************************************************** */
{
  int  i, j, k;
  static int first_call = 1;
  double *dx1, *dx2, *dx3;
  double *r, *rp, *th, *thp;
  double s, sp;
  double ***Jx1, ***Jx2, ***Jx3;
  double dx1_Bx2 = 0.0, dx2_Bx1 = 0.0;
  double dx2_Bx3 = 0.0, dx3_Bx2 = 0.0;
  double dx1_Bx3 = 0.0, dx3_Bx1 = 0.0;
  double d12, d21, d13, d31, d23, d32;
  double h3;
  static double ***Bx1, ***Bx2, ***Bx3;
  static double ***a23Bx3, ***a13Bx3, ***a12Bx2;
  RBox box;

/* --------------------------------------------------------
   1. Set pointer shortcuts.
      The primary magnetic field will be the staggered one
      or the cell-centered field.
      
      Note: in 2.5 dimensions, we also need Jx1 and Jx2 which
            contribute to the total energy equation but not
            to the induction equation.
            For this reason, we set  Bx3 to be equal to the \b
            cell centered value.
   ---------------------------------------------------------------- */

#if (defined STAGGERED_MHD)
  D_EXPAND(Bx1 = d->Vs[BX1s];  , 
           Bx2 = d->Vs[BX2s];  ,
           Bx3 = d->Vs[BX3s];)
  #if DIMENSIONS == 2 && COMPONENTS == 3
  Bx3 = d->Vc[BX3];
  #endif
#else
  if (first_call){
    Bx1 = ARRAY_3D(NX3_MAX,NX2_MAX,NX1_MAX,double);
    Bx2 = ARRAY_3D(NX3_MAX,NX2_MAX,NX1_MAX,double);
    Bx3 = ARRAY_3D(NX3_MAX,NX2_MAX,NX1_MAX,double);
  }

/* -- Construct interface values -- */

  for (k = 0; k < NX3_TOT; k++){
  for (j = 0; j < NX2_TOT; j++){
  for (i = 0; i < NX1_TOT-IOFFSET; i++){
    Bx1[k][j][i] = AVERAGE_X(d->Vc[BX1],k,j,i);
  }}}

  #if COMPONENTS >= 2
  for (k = 0; k < NX3_TOT; k++){
  for (j = 0; j < NX2_TOT-JOFFSET; j++){
  for (i = 0; i < NX1_TOT; i++){
    Bx2[k][j][i] = AVERAGE_Y(d->Vc[BX2],k,j,i);
  }}}
  #endif

  #if COMPONENTS == 3
  for (k = 0; k < NX3_TOT-KOFFSET; k++){
  for (j = 0; j < NX2_TOT; j++){
  for (i = 0; i < NX1_TOT; i++){ 
    Bx3[k][j][i] = AVERAGE_Z(d->Vc[BX3],k,j,i);
  }}}
  #endif

#endif

  Jx1 = d->J[IDIR];
  Jx2 = d->J[JDIR];
  Jx3 = d->J[KDIR];

  dx1 = grid->dx[IDIR];
  dx2 = grid->dx[JDIR];
  dx3 = grid->dx[KDIR];

/* ----------------------------------------------------------------
   2. Allocate static memory areas for the three arrays
      a12Bx2, a23Bx3, a13Bx3 if necessary. Otherwise, they will
      just be shortcuts to the default magnetic field arrays.
   ---------------------------------------------------------------- */

  if (first_call){
    a12Bx2 = Bx2;  
    a23Bx3 = Bx3;
    a13Bx3 = Bx3;
    #if GEOMETRY == CYLINDRICAL && COMPONENTS == 3
    a13Bx3 = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    #elif GEOMETRY == POLAR && COMPONENTS >=  2
    a12Bx2 = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    #elif GEOMETRY == SPHERICAL
    EXPAND(                                                    ;   ,
           a12Bx2 = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);   ,
           a13Bx3 = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
           a23Bx3 = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);)
    #endif
  }

/* --------------------------------------------------------------
   3. Construct goemetrical coefficients and arrays
   -------------------------------------------------------------- */

#if GEOMETRY == CYLINDRICAL
  r   = grid->x[IDIR]; rp  = grid->xr[IDIR];
  #if COMPONENTS == 3
  TOT_LOOP(k,j,i) a13Bx3[k][j][i] = Bx3[k][j][i]*fabs(r[i]);
  #endif
#elif GEOMETRY == POLAR
  r   = grid->x[IDIR]; rp  = grid->xr[IDIR];
  #if COMPONENTS >= 2
  TOT_LOOP(k,j,i) a12Bx2[k][j][i] = Bx2[k][j][i]*r[i];
  #endif
#elif GEOMETRY == SPHERICAL
  r   = grid->x[IDIR]; rp  = grid->xr[IDIR];
  th  = grid->x[JDIR]; thp = grid->xr[JDIR];
  TOT_LOOP(k,j,i) {
    EXPAND(                                   ;        ,
           a12Bx2[k][j][i] = Bx2[k][j][i]*r[i];        ,
           a13Bx3[k][j][i] = Bx3[k][j][i]*r[i];
           a23Bx3[k][j][i] = Bx3[k][j][i]*fabs(sin(th[j]));  )
  }
#endif

/* -------------------------------------------------------------
   4a. Compute the three components of currents at cell edges,
       so that the different components of J have different
       spatial locations:

       Jx  at  (i    , j+1/2, k+1/2)
       Jy  at  (i+1/2, j    , k+1/2)
       Jz  at  (i+1/2, j+1/2, k    )
   -------------------------------------------------------------- */

  for (k = 0; k < NX3_TOT-KOFFSET; k++){
  for (j = 0; j < NX2_TOT-JOFFSET; j++){
  for (i = 0; i < NX1_TOT-IOFFSET; i++){

    D_EXPAND(d12 = d13 = 1.0/dx1[i];  ,
             d21 = d23 = 1.0/dx2[j];  ,
             d31 = d32 = 1.0/dx3[k];)

    #if GEOMETRY == CYLINDRICAL
    d32 = d31 = 0.0;
    d13 = 2.0/(fabs(r[i+1])*r[i+1] - fabs(r[i])*r[i]);  /* = 1.0/(rp*dr) */
    #elif GEOMETRY == POLAR
    d23 /= r[i];
    d21 /= rp[i];
    d12  = 2.0/(fabs(r[i+1])*r[i+1] - fabs(r[i])*r[i]);  /* = 1.0/(rp*dr) */
    #elif GEOMETRY == SPHERICAL
    s  = fabs(sin(th[j]));
    sp = 0.5*(fabs(sin(th[j])) + fabs(sin(th[j+1]))); 
    D_EXPAND(d12 /= rp[i];   d13 /= rp[i];     ,
             d21 /= rp[i];   d23 /= r[i]*sp;   ,
             d32 /= r[i]*sp; d31 /= rp[i]*s;)
    #endif

    #if COMPONENTS == 3
    Jx1[k][j][i] = FDIFF_X2(a23Bx3,k,j,i)*d23 - FDIFF_X3(   Bx2,k,j,i)*d32;
    Jx2[k][j][i] = FDIFF_X3(   Bx1,k,j,i)*d31 - FDIFF_X1(a13Bx3,k,j,i)*d13;
    #endif
    Jx3[k][j][i] = FDIFF_X1(a12Bx2,k,j,i)*d12 - FDIFF_X2(Bx1,k,j,i)*d21;
  }}}

  #if RESISTIVITY != NO
  ComputeStaggeredEta(d, grid);  
  #endif

  first_call = 0;

}
