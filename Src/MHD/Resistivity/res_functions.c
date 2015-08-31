/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Compute the curl of magnetic field.

  Compute the electric current (defined as J = curl(B))
  for the induction and the total energy equations.

  For constrained transport MHD, J has the same staggered location of
  the electric field and the three components (Jx, Jy, Jz) are placed
  at different locations inside the cell:

  - Jx  at  (i, j+1/2, k+1/2) 
  - Jy  at  (i+1/2, j, k+1/2)
  - Jz  at  (i+1/2, j+1/2, k)

  The same rule apply to the components of resistivity \c eta which
  are computed and stored inside this function.

  For cell-centered MHD, the three components of J are computed
  during each sweep direction at cell interfaces, that is,
  
  - IDIR:  (Jx, Jy, Jz)  at  (i+1/2, j, k)
  - JDIR:  (Jx, Jy, Jz)  at  (i, j+1/2, k)
  - KDIR:  (Jx, Jy, Jz)  at  (i, j, k+1/2)

  Although only the two transverse components of J are actually needed
  during the update step, we compute also the normal component 
  since the resistivity coefficients eta may depend on the total current.

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
        singular behavior with the definition of d13 (d12) at r=0.
        An equivalent formulation can be obtained by using absolute
        values in the geometrical coefficient and thus replacing
        d13 = r[i+1/2]*dr with (r[i+1] |r[i+1]| - r[i] |r[i]|)/2 and
        a13 = r[i] with abs(r[i]).                
        This is obviously the desired derivative away from the axis 
        (where r = |r|) and it gives the correct behavior at the 
        axis, where Bphi --> -Bphi.          
        Bphi = 0  --> Bphi \sim alpha*r --> d(r Bphi)/(r dr) = 2*alpha 
           

  \attention Do NOT use this function to compute terms like (curl V)^2 !!

  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos\n
           
  \date   March 10, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static double ***eta[3];

/* ********************************************************************* */
void GetCurrent (const Data *d, int dir, Grid *grid)
/*!
 * Compute the curl of magnetic field for constrained transport MHD or
 * cell-centered MHD.
 *
 * \param  [in,out]  d    pointer to the PLUTO data structure
 * \param  [in]    dir    sweep direction (useless for constrained
 *                        transport MHD)
 * \param [in]    grid    pointer to an array of Grid structures
 *   
 *********************************************************************** */
{
  int  i, j, k;
  static int first_call = 1;
  double *dx1, *dx2, *dx3;
  double *r, *rp, *th, *thp;
  double ***Bx1, ***Bx2, ***Bx3;
  double ***Jx1, ***Jx2, ***Jx3;
  double dx1_Bx2 = 0.0, dx2_Bx1 = 0.0;
  double dx2_Bx3 = 0.0, dx3_Bx2 = 0.0;
  double dx1_Bx3 = 0.0, dx3_Bx1 = 0.0;
  double d12, d21, d13, d31, d23, d32;
  double h3;
  static double ***a23Bx3, ***a13Bx3, ***a12Bx2;
  RBox box;

/* ------------------------------------------------------------------
   1. Set pointer shortcuts.
      The primary magnetic field will be the staggered one
      (if CONSTRAINTED_TRANSPORT is used) or the cell-centered field
      otherwise.
      
      Note: in 2+1/2 dimensions, we also need Jx1 and Jx2 which
            contribute to the total energy equation but not
            to the induction equation.
            For this reason, we set  Bx3 to be equal to the \b
            cell centered value.
   ---------------------------------------------------------------- */

  #ifdef STAGGERED_MHD
   D_EXPAND(Bx1 = d->Vs[BX1s];  , 
            Bx2 = d->Vs[BX2s];  ,
            Bx3 = d->Vs[BX3s];)
   #if DIMENSIONS == 2 && COMPONENTS == 3
    Bx3 = d->Vc[BX3];
   #endif
  #else
   EXPAND(Bx1 = d->Vc[BX1];    ,
          Bx2 = d->Vc[BX2];    ,
          Bx3 = d->Vc[BX3];)
  #endif

  Jx1 = d->J[IDIR];
  Jx2 = d->J[JDIR];
  Jx3 = d->J[KDIR];

  dx1 = grid[IDIR].dx;
  dx2 = grid[JDIR].dx;
  dx3 = grid[KDIR].dx;

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
    first_call = 0;
  }

/* --------------------------------------------------------------
   3. Construct goemetrical coefficients and arrays
   -------------------------------------------------------------- */

  #if GEOMETRY == CYLINDRICAL
   r   = grid[IDIR].x; rp  = grid[IDIR].xr;
   #if COMPONENTS == 3
    TOT_LOOP(k,j,i) a13Bx3[k][j][i] = Bx3[k][j][i]*fabs(r[i]);
   #endif
  #elif GEOMETRY == POLAR
   r   = grid[IDIR].x; rp  = grid[IDIR].xr;
   #if COMPONENTS >= 2
    TOT_LOOP(k,j,i) a12Bx2[k][j][i] = Bx2[k][j][i]*r[i];
   #endif
  #elif GEOMETRY == SPHERICAL
   r   = grid[IDIR].x; rp  = grid[IDIR].xr;
   th  = grid[JDIR].x; thp = grid[JDIR].xr;
   TOT_LOOP(k,j,i) {
     EXPAND(                                   ;        ,
            a12Bx2[k][j][i] = Bx2[k][j][i]*r[i];        ,
            a13Bx3[k][j][i] = Bx3[k][j][i]*r[i];
            a23Bx3[k][j][i] = Bx3[k][j][i]*sin(th[j]);  )
   }
  #endif

/* -------------------------------------------------------------
   4a. For staggered MHD, we compute the three components of
       currents at cell edges.
       In a staggered formulation the different components of J
       have different spatial location (same as electric field)
       and one needs

       Jx  at  (i    , j+1/2, k+1/2)
       Jy  at  (i+1/2, j    , k+1/2)
       Jz  at  (i+1/2, j+1/2, k    )
   -------------------------------------------------------------- */

#ifdef STAGGERED_MHD

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
     D_EXPAND(d12 /= rp[i];             d13 /= rp[i];              ,
              d21 /= rp[i];             d23 /= r[i]*sin(thp[j]);   ,
              d32 /= r[i]*sin(thp[j]);  d31 /= rp[i]*sin(th[j]);)
    #endif

    #if COMPONENTS == 3
     Jx1[k][j][i] = FDIFF_X2(a23Bx3,k,j,i)*d23 - FDIFF_X3(   Bx2,k,j,i)*d32;
     Jx2[k][j][i] = FDIFF_X3(   Bx1,k,j,i)*d31 - FDIFF_X1(a13Bx3,k,j,i)*d13;
    #endif
    Jx3[k][j][i] = FDIFF_X1(a12Bx2,k,j,i)*d12 - FDIFF_X2(Bx1,k,j,i)*d21;
  }}}
  ComputeStaggeredEta(d, grid);

#else

/* ------------------------------------------------------------- */
/*!
   4b. For cell-centered MHD, we compute the three components of
        currents at a cell interface.
  -------------------------------------------------------------- */

  if (dir == IDIR){

  /* ----------------------------------------------------
      IDIR: Compute {Jx1, Jx2, Jx3} at i+1/2,j,k faces.
     ---------------------------------------------------- */

    box.ib =       0; box.ie = NX1_TOT-1-IOFFSET;
    box.jb = JOFFSET; box.je = NX2_TOT-1-JOFFSET;
    box.kb = KOFFSET; box.ke = NX3_TOT-1-KOFFSET;

    BOX_LOOP(&box,k,j,i){

      d12 = d13 = 1.0/dx1[i];
      d21 = d23 = 1.0/dx2[j]; 
      d31 = d32 = 1.0/dx3[k];    
      
      #if GEOMETRY == CYLINDRICAL
       d32 = d31 = 0.0;
       d13 = 2.0/(fabs(r[i+1])*r[i+1] - fabs(r[i])*r[i]);  
      #elif GEOMETRY == POLAR
       d23 = d21 = 1.0/(rp[i]*dx2[j]); 
       d12 = 2.0/(fabs(r[i+1])*r[i+1] - fabs(r[i])*r[i]);  
      #elif GEOMETRY == SPHERICAL
       h3  = rp[i]*sin(th[j]);
       
       D_EXPAND(d12 = 1.0/(rp[i]*dx1[i]); d13 = 1.0/(rp[i]*dx1[i]);   ,
                d21 = 1.0/(rp[i]*dx2[j]); d23 = 1.0/(h3*dx2[j]);      ,
                d32 = 1.0/(h3*dx3[k]);    d31 = 1.0/(h3*dx3[k]);)
      #endif

      #if COMPONENTS == 3
       dx2_Bx3 = 0.5*(CDIFF_X2(a23Bx3,k,j,i) + CDIFF_X2(a23Bx3,k,j,i+1))*d23;
       dx3_Bx2 = 0.5*(CDIFF_X3(Bx2,k,j,i)    + CDIFF_X3(Bx2,k,j,i+1)  )*d32;
       Jx1[k][j][i] = (dx2_Bx3 - dx3_Bx2);

       dx3_Bx1 = 0.5*(CDIFF_X3(Bx1,k,j,i) + CDIFF_X3(Bx1,k,j,i+1))*d31;
       dx1_Bx3 = FDIFF_X1(a13Bx3,k,j,i)*d13;
       Jx2[k][j][i] = (dx3_Bx1 - dx1_Bx3);
      #endif

      dx1_Bx2 = FDIFF_X1(a12Bx2,k,j,i)*d12;
      dx2_Bx1 = 0.5*(CDIFF_X2(Bx1,k,j,i) + CDIFF_X2(Bx1,k,j,i+1))*d21;
      Jx3[k][j][i] = (dx1_Bx2 - dx2_Bx1);
    }
    
  }else if (dir == JDIR){

  /* ----------------------------------------------------
      JDIR: Compute {Jx1, Jx2, Jx3} at i,j+1/2,k faces.
     ---------------------------------------------------- */

    box.ib = IOFFSET; box.ie = NX1_TOT-1-IOFFSET;
    box.jb =       0; box.je = NX2_TOT-1-JOFFSET;
    box.kb = KOFFSET; box.ke = NX3_TOT-1-KOFFSET;

    BOX_LOOP(&box,k,j,i){
      d12 = d13 = 1.0/dx1[i];
      d21 = d23 = 1.0/dx2[j];
      d31 = d32 = 1.0/dx3[k];

      #if GEOMETRY == CYLINDRICAL
       d32 = d31 = 0.0; /* axisymmetry */
       d13 = 1.0/(r[i]*dx1[i]); 
      #elif GEOMETRY == POLAR
       d23 = d21 = 1.0/(r[i]*dx2[j]); 
       d12 = 1.0/(r[i]*dx1[i]); 
      #elif GEOMETRY == SPHERICAL
       h3  = r[i]*sin(thp[j]);

/* !!! check singularity !!! */

       D_EXPAND(d12 = 1.0/(r[i]*dx1[i]); d13 = 1.0/(r[i]*dx1[i]);   ,
                d21 = 1.0/(r[i]*dx2[j]); d23 = 1.0/(h3*dx2[j]);     ,
                d32 = 1.0/(h3*dx3[k]);   d31 = 1.0/(h3*dx3[k]);)
      #endif
      
      #if COMPONENTS == 3
       dx2_Bx3 = FDIFF_X2(a23Bx3,k,j,i)*d23;
       dx3_Bx2 = 0.5*(CDIFF_X3(Bx2,k,j,i) + CDIFF_X3(Bx2,k,j+1,i))*d32;
       Jx1[k][j][i] = (dx2_Bx3 - dx3_Bx2);

       dx3_Bx1 = 0.5*(CDIFF_X3(Bx1,k,j,i)    + CDIFF_X3(Bx1,k,j+1,i))*d31;
       dx1_Bx3 = 0.5*(CDIFF_X1(a13Bx3,k,j,i) + CDIFF_X1(a13Bx3,k,j+1,i))*d13;
       Jx2[k][j][i] = (dx3_Bx1 - dx1_Bx3);
      #endif
      dx1_Bx2 = 0.5*(CDIFF_X1(a12Bx2,k,j,i) + CDIFF_X1(a12Bx2,k,j+1,i))*d12;
      dx2_Bx1 = FDIFF_X2(Bx1,k,j,i)*d21;
      Jx3[k][j][i] = (dx1_Bx2 - dx2_Bx1);
    }

  }else if (dir == KDIR) {

  /* ----------------------------------------------------
      KDIR: Compute {Jx1, Jx2, Jx3} at i,j,k+1/2 faces.
     ---------------------------------------------------- */

    box.ib = IOFFSET; box.ie = NX1_TOT-1-IOFFSET;
    box.jb = JOFFSET; box.je = NX2_TOT-1-JOFFSET;
    box.kb =       0; box.ke = NX3_TOT-1-KOFFSET;

    BOX_LOOP(&box,k,j,i){
      d12 = d13 = 1.0/dx1[i];
      d21 = d23 = 1.0/dx2[j];
      d31 = d32 = 1.0/dx3[k];

      #if GEOMETRY == POLAR
       d23 = d21 = 1.0/(rp[i]*dx2[j]); 
       d12 = 1.0/(rp[i]*dx1[i]); 
      #elif GEOMETRY == SPHERICAL
       h3  = r[i]*sin(th[j]);
       
       D_EXPAND(d12 = 1.0/(r[i]*dx1[i]); d13 = 1.0/(r[i]*dx1[i]);   ,
                d21 = 1.0/(r[i]*dx2[j]); d23 = 1.0/(h3*dx2[j]);     ,
                d32 = 1.0/(h3*dx3[k]);   d31 = 1.0/(h3*dx3[k]);)
      #endif
      
      #if COMPONENTS == 3
       dx2_Bx3 = 0.5*(CDIFF_X2(a23Bx3,k,j,i) + CDIFF_X2(a23Bx3,k+1,j,i))*d23;
       dx3_Bx2 = FDIFF_X3(Bx2,k,j,i)*d32;
       Jx1[k][j][i] = (dx2_Bx3 - dx3_Bx2);

       dx3_Bx1 = FDIFF_X3(Bx1,k,j,i)*d31;
       dx1_Bx3 = 0.5*(CDIFF_X1(a13Bx3,k,j,i) + CDIFF_X1(a13Bx3,k+1,j,i))*d13;
       Jx2[k][j][i] = (dx3_Bx1 - dx1_Bx3);
      #endif
      dx1_Bx2 = 0.5*(CDIFF_X1(a12Bx2,k,j,i) + CDIFF_X1(a12Bx2,k+1,j,i))*d12;
      dx2_Bx1 = 0.5*(CDIFF_X2(Bx1,k,j,i)    + CDIFF_X2(Bx1,k+1,j,i))*d21;
      Jx3[k][j][i] = (dx1_Bx2 - dx2_Bx1);
    }
  }

#endif  /* STAGGERED_MHD */
}

#ifdef STAGGERED_MHD
/* ********************************************************************* */
void ComputeStaggeredEta(const Data *d, Grid *grid)
/*!
 *  Compute eta coefficients using the same staggering of the current:
 *
 *   -  eta_x1 at   (i, j+1/2, k+1/2)
 *   -  eta_x2 at   (i+1/2, j, k+1/2)
 *   -  eta_x3 at   (i+1/2, j+1/2, k)
 * 
 *  Since each eta_x may depend on the total current, we also need
 *  to interpolate different current components at the same place.
 * 
 *********************************************************************** */
{
  int    i, j, k, nv;
  double Js[3], etas[3], vs[NVAR];
  double ***Jx1, ***Jx2, ***Jx3;
  double *x1,  *x2,  *x3;
  double *x1r, *x2r, *x3r;

  if (eta[0] == NULL){
    eta[IDIR] = ARRAY_3D(NX3_TOT,NX2_TOT,NX1_TOT,double);
    eta[JDIR] = ARRAY_3D(NX3_TOT,NX2_TOT,NX1_TOT,double);
    eta[KDIR] = ARRAY_3D(NX3_TOT,NX2_TOT,NX1_TOT,double);
  }

  Jx1 = d->J[IDIR]; x1 = grid[IDIR].x; x1r = grid[IDIR].xr;
  Jx2 = d->J[JDIR]; x2 = grid[JDIR].x; x2r = grid[JDIR].xr;
  Jx3 = d->J[KDIR]; x3 = grid[KDIR].x; x3r = grid[KDIR].xr;

  for (k = 0; k < NX3_TOT-KOFFSET; k++){
  for (j = 0; j < NX2_TOT-JOFFSET; j++){
  for (i = 0; i < NX1_TOT-IOFFSET; i++){

    #if COMPONENTS == 2

  /* ---- Compute J and eta_x3 at  i+1/2, j+1/2, k ---- */

     Js[IDIR] = 0.0;
     Js[JDIR] = 0.0;
     Js[KDIR] = Jx3[k][j][i];
     VAR_LOOP(nv) vs[nv] = AVERAGE_XY(d->Vc[nv],k,j,i);
     Resistive_eta (vs, x1r[i], x2r[j], x3[k], Js, etas);
     eta[KDIR][k][j][i] = etas[KDIR];

    #elif COMPONENTS == 3

  /* ---- Compute J and eta_x1 at  i, j+1/2, k+1/2 ---- */

     if (i > 0) {
       Js[IDIR] = Jx1[k][j][i];
       Js[JDIR] = AVERAGE_XY(Jx2,k,j,i-1);
       Js[KDIR] = AVERAGE_XZ(Jx3,k,j,i-1);
       VAR_LOOP(nv) vs[nv] = AVERAGE_YZ(d->Vc[nv],k,j,i);
       Resistive_eta (vs, x1[i], x2r[j], x3r[k], Js, etas);
       eta[IDIR][k][j][i] = etas[IDIR];
     }

  /* ---- Compute J and eta_x2 at  i+1/2, j, k+1/2 ---- */

     if (j > 0) {
       Js[IDIR] = AVERAGE_XY(Jx1,k,j-1,i);
       Js[JDIR] = Jx2[k][j][i];
       Js[KDIR] = AVERAGE_YZ(Jx3,k,j-1,i);
       VAR_LOOP(nv) vs[nv] = AVERAGE_XZ(d->Vc[nv],k,j,i);
       Resistive_eta (vs, x1r[i], x2[j], x3r[k], Js, etas);
       eta[JDIR][k][j][i] = etas[JDIR];
     }

  /* ---- Compute J and eta_x3 at  i+1/2, j+1/2, k ---- */

     if (k > 0 || KOFFSET == 0){
       Js[IDIR] = AVERAGE_XZ(Jx1,k-1,j,i);
       Js[JDIR] = AVERAGE_YZ(Jx2,k-1,j,i);
       Js[KDIR] = Jx3[k][j][i];
       VAR_LOOP(nv) vs[nv] = AVERAGE_XY(d->Vc[nv],k,j,i);
       Resistive_eta (vs, x1r[i], x2r[j], x3[k], Js, etas);
       eta[KDIR][k][j][i] = etas[KDIR];
     }

    #endif

  }}}

}

/* ********************************************************************* */
Data_Arr GetStaggeredEta()
/*
 *********************************************************************** */
{
  return eta;
}
#endif
