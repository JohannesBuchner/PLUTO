/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Set boundary conditions.

  The Boundary() function sets both internal and physical boundary 
  conditions on one or more sides of the computational domain.
  It is used to fill ghost zones of both cell-centered and face-centered 
  data arrays.\n
  The type of boundary conditions at the leftmost or rightmost side of a 
  given grid is specified by the integers <c> grid[dir].lbound </c> or 
  <c> grid[dir].rbound </c>, respectively.
  When this value is different from zero, the local processor borders 
  the physical boundary and the admissible values for \c lbound or \c 
  rbound are OUTFLOW, REFLECTIVE, AXISYMMETRIC, EQTSYMMETRIC, PERIODIC, 
  SHEARING or USERDEF.
  Conversely, when this value is zero (internal boundary), two neighboring
  processors that share the same side need to fill ghost zones by exchanging 
  data values. 
  This step is done here only for parallel computations on static grids.
  
  Predefined physical boundary conditions are handled by the 
  following functions:
  
  - OutflowBound()
  - ReflectiveBound()
  - PeriodicBound()

  \author A. Mignone (mignone@ph.unito.it)
  \date   Dec 18, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"
                           
/* ********************************************************************* */
void Boundary (const Data *d, int idim, Grid *grid)
/*!
 * Set boundary conditions on one or more sides of the computational
 * domain.
 *
 * \param [in,out] d  pointer to PLUTO Data structure containing the 
 *                    solution array (including centered and staggered
 *                    fields)
 * \param [in]   idim specifies on which side(s) of the computational
 *               domain boundary conditions must be set. Possible values
 *               are  
 *        - idim = IDIR   first dimension (x1)
 *        - idim = JDIR   second dimenson (x2)
 *        - idim = KDIR   third dimension (x3)
 *        - idim = ALL_DIR all dimensions
 *
 * \param [in]  grid   pointer to an array of grid structures.
 *********************************************************************** */
{
  int  is, nv, fargo_velocity_has_changed;
  int  side[6] = {X1_BEG, X1_END, X2_BEG, X2_END, X3_BEG, X3_END};
  int  type[6], sbeg, send, vsign[NVAR];
  int  par_dim[3] = {0, 0, 0};
  double ***q;

/* ---------------------------------------------------
    Check the number of processors in each direction
   --------------------------------------------------- */

  D_EXPAND(par_dim[0] = grid[IDIR].nproc > 1;  ,
           par_dim[1] = grid[JDIR].nproc > 1;  ,
           par_dim[2] = grid[KDIR].nproc > 1;)

/* -------------------------------------------------
    With FARGO, boundary conditions must be set on 
    total velocity. 
    We temporarily add the background motion if 
    the input array contains the deviation only.
   ------------------------------------------------- */
   
   #ifdef FARGO
    fargo_velocity_has_changed = NO;
    if (FARGO_HasTotalVelocity() == NO) {
      FARGO_AddVelocity (d,grid);
      fargo_velocity_has_changed = YES;
    }
   #endif

/* -------------------------------------------------
    Call userdef internal boundary with side == 0
   -------------------------------------------------  */

  #if INTERNAL_BOUNDARY == YES
   UserDefBoundary (d, NULL, 0, grid);
  #endif
  
/* -------------------------------------
     Exchange data between processors 
   ------------------------------------- */
   
  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
   for (nv = 0; nv < NVAR; nv++) {
     AL_Exchange_dim ((char *)d->Vc[nv][0][0], par_dim, SZ);
   }
   #ifdef STAGGERED_MHD 
    D_EXPAND(
      AL_Exchange_dim ((char *)(d->Vs[BX1s][0][0] - 1), par_dim, SZ_stagx);  ,
      AL_Exchange_dim ((char *)d->Vs[BX2s][0][-1]     , par_dim, SZ_stagy);  ,
      AL_Exchange_dim ((char *)d->Vs[BX3s][-1][0]     , par_dim, SZ_stagz);)
   #endif
   MPI_Barrier (MPI_COMM_WORLD);
  #endif

/* ----------------------------------------------------------------
     When idim == ALL_DIR boundaries are imposed on ALL sides:
     a loop from sbeg = 0 to send = 2*DIMENSIONS - 1 is performed. 
     When idim = n, boundaries are imposed at the beginning and 
     the end of the i-th direction.
   ---------------------------------------------------------------- */ 

  if (idim == ALL_DIR) {
    sbeg = 0;
    send = 2*DIMENSIONS - 1;
  } else {
    sbeg = 2*idim;
    send = 2*idim + 1;
  }

/* --------------------------------------------------------
        Main loop on computational domain sides
   -------------------------------------------------------- */

  type[0] = grid[IDIR].lbound; type[1] = grid[IDIR].rbound;
  type[2] = grid[JDIR].lbound; type[3] = grid[JDIR].rbound;
  type[4] = grid[KDIR].lbound; type[5] = grid[KDIR].rbound;

  for (is = sbeg; is <= send; is++){

    if (type[is] == 0) continue;  /* no physical boundary: skip */

    if (type[is] == OUTFLOW) {

    /* -------------------------------
         OUTFLOW boundary condition 
       ------------------------------- */
 
      for (nv = 0; nv < NVAR; nv++){
        OutflowBound (d->Vc[nv], side[is], CENTER, grid + (is/2));   
      }
      #ifdef STAGGERED_MHD 
       D_EXPAND(OutflowBound (d->Vs[BX1s], side[is], X1FACE, grid+(is/2)); ,
                OutflowBound (d->Vs[BX2s], side[is], X2FACE, grid+(is/2)); ,
                OutflowBound (d->Vs[BX3s], side[is], X3FACE, grid+(is/2));)

    /* ---------------------------------------------------------
        average normal field only since transverse components
        are assigned consistently with cell-centered quantities.
       --------------------------------------------------------- */
            
       CT_AverageNormalMagField (d, side[is], grid);
      #endif

    }else if (  (type[is] == REFLECTIVE) 
             || (type[is] == AXISYMMETRIC)
             || (type[is] == EQTSYMMETRIC)){

    /* -------------------------------------
        REFLECTIVE-type boundary conditions 
       ------------------------------------- */

      FlipSign (side[is], type[is], vsign);
      for (nv = 0; nv < NVAR; nv++){
        ReflectiveBound (d->Vc[nv], vsign[nv], side[is], CENTER);
      }
      #ifdef STAGGERED_MHD  
       D_EXPAND(ReflectiveBound(d->Vs[BX1s], vsign[BX1], side[is], X1FACE);  ,
                ReflectiveBound(d->Vs[BX2s], vsign[BX2], side[is], X2FACE);  ,
                ReflectiveBound(d->Vs[BX3s], vsign[BX3], side[is], X3FACE);)
      #endif

    }else if (type[is] == PERIODIC) {

    /* ----------------------------------------
        PERIODIC boundary condition (only for 
        one processor in the direction) 
       ---------------------------------------- */

      if (!par_dim[is/2]) {
        for (nv = 0; nv < NVAR; nv++){
          PeriodicBound (d->Vc[nv], side[is], CENTER);
        }
        #ifdef STAGGERED_MHD
         D_EXPAND(PeriodicBound (d->Vs[BX1s], side[is], X1FACE);  ,
                  PeriodicBound (d->Vs[BX2s], side[is], X2FACE);  ,
                  PeriodicBound (d->Vs[BX3s], side[is], X3FACE);)
        #endif
 
      }

    }else if (type[is] == SHEARING) {  /* -- shearingbox boundary -- */

    /* ---------------------------------------------------------
         SHEARING-BOX boundary condition is implemented as

        1) apply periodic boundary conditions for all variables
          (except staggered BX)
        2) Perform spatial shift in the y-direction
       --------------------------------------------------------- */

      #ifdef SHEARINGBOX
       if (side[is] != X1_BEG && side[is] != X1_END){
         print1 ("! BOUNDARY: shearingbox can only be assigned at an X1 boundary\n");
         QUIT_PLUTO(1);
       }
       if (grid[IDIR].nproc == 1){
         for (nv = 0; nv < NVAR; nv++) {
           PeriodicBound (d->Vc[nv], side[is], CENTER);
         }
         #ifdef STAGGERED_MHD
          D_EXPAND(                                             ;  ,
                   PeriodicBound (d->Vs[BX2s], side[is], X2FACE);  ,
                   PeriodicBound (d->Vs[BX3s], side[is], X3FACE);)
         #endif
       }
       SB_Boundary (d, side[is], grid);

     /* -- assign normal component of staggered B 
           using the div.B = 0 condition           -- */

       #ifdef STAGGERED_MHD
        FillMagneticField (d, side[is], grid);  
        CT_AverageNormalMagField (d, side[is], grid);
        CT_AverageTransverseMagField (d, side[is], grid);
       #endif
      #endif  /* def SHEARINGBOX */

    }else if (type[is] == USERDEF) {

    /* --------------------------------------------------------------
             USER-DEFINED boundary condition
       -------------------------------------------------------------- */

      #ifdef GLM_MHD
      {int i,j,k;
       switch(side[is]){
         case X1_BEG: X1_BEG_LOOP(k,j,i) d->Vc[PSI_GLM][k][j][i] = 0.0; break;
         case X1_END: X1_END_LOOP(k,j,i) d->Vc[PSI_GLM][k][j][i] = 0.0; break;
         case X2_BEG: X2_BEG_LOOP(k,j,i) d->Vc[PSI_GLM][k][j][i] = 0.0; break;
         case X2_END: X2_END_LOOP(k,j,i) d->Vc[PSI_GLM][k][j][i] = 0.0; break;
         case X3_BEG: X3_BEG_LOOP(k,j,i) d->Vc[PSI_GLM][k][j][i] = 0.0; break;
         case X3_END: X3_END_LOOP(k,j,i) d->Vc[PSI_GLM][k][j][i] = 0.0; break;
       }
      }
      #endif

      RBox *box = GetRBox(side[is], CENTER);
      UserDefBoundary (d, box, side[is], grid);
      #ifdef STAGGERED_MHD
        D_EXPAND(box = GetRBox(side[is], X1FACE);
                 UserDefBoundary (d, box, side[is], grid);  ,
                 box = GetRBox(side[is], X2FACE);
                 UserDefBoundary (d, box, side[is], grid);  ,
                 box = GetRBox(side[is], X3FACE);
                 UserDefBoundary (d, box, side[is], grid);)

       /* -- assign normal component of staggered B 
             using the div.B = 0 condition           -- */

       FillMagneticField (d, side[is], grid);  
       CT_AverageNormalMagField (d, side[is], grid);
       CT_AverageTransverseMagField (d, side[is], grid);
      #endif
    }
  }

/* --------------------------------------------------
    With FARGO, restore velocity deviation if the 
    original input array to Boundary() did not 
    contain the background velocity.
   -------------------------------------------------- */
   
   #ifdef FARGO
    if (fargo_velocity_has_changed){
      if (FARGO_HasTotalVelocity() == YES) FARGO_SubtractVelocity(d,grid);
    }
   #endif

/* -------------------------------------------------------
    Compute entropy for the next time level
   ------------------------------------------------------- */

  #if ENTROPY_SWITCH
   ComputeEntropy (d, grid);
  #endif
    
}

/* ********************************************************************* */
void OutflowBound (double ***q, int side, int vpos, Grid *grid)
/*! 
 * Impose zero-gradient boundary conditions on 'q' on 
 * the boundary side specified by box->side.
 * The staggered component of magnetic field is assigned
 * using the div.B = 0 condition in a different loop.
 *
 * \param [in,out] q   a 3D array representing a flow variable
 * \param [in]    box pointer to a RBox structure defining the extent
 *                    of the boundary region and the variable position
 *                    inside the cell
 * \param [in]    grid  pointer to an array of Grid structures
 *********************************************************************** */
{
  int  i, j, k;
  double dB, *A, *dV;
  RBox *box = GetRBox (side, vpos);

  A  = grid->A;
  dV = grid->dV;

  if (side == X1_BEG) { 
    if (box->vpos != X1FACE){
      BOX_LOOP(box,k,j,i) q[k][j][i] = q[k][j][IBEG];   
    }else{
      BOX_LOOP(box,k,j,i){
        #if GEOMETRY == CARTESIAN
         q[k][j][i] = 2.0*q[k][j][i+1] - q[k][j][i+2];
        #else
         dB = (A[i+2]*q[k][j][i+2] - A[i+1]*q[k][j][i+1])/dV[i+2];
         q[k][j][i] = (q[k][j][i+1]*A[i+1] - dV[i+1]*dB)/A[i];
        #endif
      }
    }

  }else if (side == X1_END){

    if (box->vpos != X1FACE){
      BOX_LOOP(box,k,j,i) q[k][j][i] = q[k][j][IEND];   
    }else{
      BOX_LOOP(box,k,j,i){
        #if GEOMETRY == CARTESIAN
         q[k][j][i] = 2.0*q[k][j][i-1] - q[k][j][i-2];
        #else
         dB = (A[i-1]*q[k][j][i-1] - A[i-2]*q[k][j][i-2])/dV[i-1];
         q[k][j][i] = (q[k][j][i-1]*A[i-1] + dV[i]*dB)/A[i];
        #endif
      }
    }

  }else if (side == X2_BEG){

    if (box->vpos != X2FACE) {
      BOX_LOOP(box,k,j,i) q[k][j][i] = q[k][JBEG][i];
    }else{
      BOX_LOOP(box,k,j,i) {
        #if GEOMETRY == CARTESIAN
         q[k][j][i] = 2.0*q[k][j+1][i] - q[k][j+2][i];
        #else
         dB = (A[j+2]*q[k][j+2][i] - A[j+1]*q[k][j+1][i])/dV[j+2];
         q[k][j][i] = (q[k][j+1][i]*A[j+1] - dV[j+1]*dB)/A[j];
        #endif
      }
    }

  }else if (side == X2_END){

    if (box->vpos != X2FACE) {
      BOX_LOOP(box,k,j,i) q[k][j][i] = q[k][JEND][i];
    }else{
      BOX_LOOP(box,k,j,i){
        #if GEOMETRY == CARTESIAN
         q[k][j][i] = 2.0*q[k][j-1][i] - q[k][j-2][i];
        #else
         dB = (A[j-1]*q[k][j-1][i] - A[j-2]*q[k][j-2][i])/dV[j-1];
         q[k][j][i] = (q[k][j-1][i]*A[j-1] + dV[j]*dB)/A[j];
        #endif
      }
    }

  }else if (side == X3_BEG){

    if (box->vpos != X3FACE) {
      BOX_LOOP(box,k,j,i) q[k][j][i] = q[KBEG][j][i];
    }else{
      BOX_LOOP(box,k,j,i) {
        #if GEOMETRY == CARTESIAN
         q[k][j][i] = 2.0*q[k+1][j][i] - q[k+2][j][i];
        #else
         dB = (A[k+2]*q[k+2][j][i] - A[k+1]*q[k+1][j][i])/dV[k+2];
         q[k][j][i] = (q[k+1][j][i]*A[k+1] - dV[k+1]*dB)/A[k];
        #endif
      }
    }

  }else if (side == X3_END){

    if (box->vpos != X3FACE) {
      BOX_LOOP(box,k,j,i) q[k][j][i] = q[KEND][j][i];
    }else{
      BOX_LOOP(box,k,j,i) {
        #if GEOMETRY == CARTESIAN
         q[k][j][i] = 2.0*q[k-1][j][i] - q[k-2][j][i];  
        #else
         dB = (A[k-1]*q[k-1][j][i] - A[k-2]*q[k-2][j][i])/dV[k-1];
         q[k][j][i] = (q[k-1][j][i]*A[k-1] + dV[k]*dB)/A[k];
        #endif
      }
    }
  }
}

/* ********************************************************************* */
void FlipSign (int side, int type, int *vsign)
/*!
 * Reverse the sign of vector components with respect to axis side. 
 * Depending on type, one needs to symmetrize or anti-symmetrize:
 *
 * - REFLECTIVE:    \n
 *   o   Vn -> -Vn,  Bn -> -Bn  \n
 *   o   Vp ->  Vp,  Bp ->  Bp  \n 
 *   o   Vt ->  Vt,  Bt ->  Bt
 *
 * - AXISYMMETRIC:  \n
 *   o   Vn -> -Vn,  Bn -> -Bn \n
 *   o   Vp ->  Vp,  Bp ->  Bp \n
 *   o   Vt -> -Vt,  Bt -> -Bt
 *
 * - EQTSYMMETRIC:  \n 
 *   o   Vn -> -Vn,  Bn ->  Bn  \n
 *   o   Vp ->  Vp,  Bp -> -Bp  \n
 *   o   Vt ->  Vt,  Bt -> -Bt
 *
 * where (n) is the normal components, (p) and (t)
 * are the transverse (or poloidal and toroidal for
 * cylindrical and spherical coordinates) components.
 * 
 * \param [in]  side  boundary side
 * \param [in]  type boundary condition type
 * \param [out] vsign an array of values (+1 or -1) giving the sign
 *********************************************************************** */
{
  int nv;
  int Vn, Vp, Vt;
  int Bn, Bp, Bt;

  #if PHYSICS == ADVECTION
   for (nv = 0; nv < NVAR; nv++) vsign[nv] = 1.0;
   return;
  #endif

/* ----------------------------------------------------------
    get normal (n), poloidal (p) and toroidal (t) vector 
    components. The ordering is the same as in SetIndexes()
   ---------------------------------------------------------- */

  if (side == X1_BEG || side == X1_END){
    Vn = VX1; Vp = VX2; Vt = VX3;
    #if PHYSICS == MHD || PHYSICS == RMHD
     Bn = BX1; Bp = BX2; Bt = BX3;
    #endif
  }else if (side == X2_BEG || side == X2_END){
    Vn = VX2; Vp = VX1; Vt = VX3;
    #if PHYSICS == MHD || PHYSICS == RMHD
     Bn = BX2; Bp = BX1; Bt = BX3;
    #endif
  }else if (side == X3_BEG || side == X3_END){
    Vn = VX3; Vp = VX1; Vt = VX2;
    #if PHYSICS == MHD || PHYSICS == RMHD
     Bn = BX3; Bp = BX1; Bt = BX2;
    #endif
  }

/* ---------------------------------------
     decide which variable flips sign
   --------------------------------------- */

  for (nv = 0; nv < NVAR; nv++) vsign[nv] = 1.0;

  vsign[Vn] = -1.0;
  #if PHYSICS == MHD || PHYSICS == RMHD
   vsign[Bn] = -1.0;
  #endif

  #if COMPONENTS == 3
   if (type == AXISYMMETRIC){
     vsign[Vt] = -1.0;
     #if PHYSICS == MHD || PHYSICS == RMHD
      vsign[Bt] = -1.0;
     #endif
   }
  #endif

  #if PHYSICS == MHD || PHYSICS == RMHD
   if (type == EQTSYMMETRIC){
     EXPAND(vsign[Bn] =  1.0; ,
            vsign[Bp] = -1.0; ,
            vsign[Bt] = -1.0;)
     #ifdef GLM_MHD 
      vsign[PSI_GLM] = -1.0;
     #endif
   }
  #endif
}

/* ********************************************************************* */
void ReflectiveBound (double ***q, int s, int side, int vpos)
/*!
 * Make symmetric (s = 1) or anti-symmetric (s=-1) profiles 
 * with respect to the boundary plane specified by box->side.
 * The sign is set by the FlipSign() function. 
 *
 * \param [in,out] q   a 3D flow quantity
 * \param [in] box pointer to box structure
 * \param [in] s   an integer taking only the values +1 (symmetric 
 *                 profile) or -1 (antisymmetric profile)
 *   
 *********************************************************************** */
{
  int   i, j, k;
  RBox *box = GetRBox(side, vpos);

  if (side == X1_BEG) {   

    if (box->vpos != X1FACE){
      BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[k][j][2*IBEG-i-1];
    }else{
      BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[k][j][2*IBEG-i-2];
    }

  }else if (side == X1_END){  

    if (box->vpos != X1FACE){
      BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[k][j][2*IEND-i+1];
    }else{
      BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[k][j][2*IEND-i];
    }

  }else if (side == X2_BEG){  

    if (box->vpos != X2FACE){
      BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[k][2*JBEG-j-1][i];
    }else{
      BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[k][2*JBEG-j-2][i];
    }

  }else if (side == X2_END){

    if (box->vpos != X2FACE){
      BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[k][2*JEND-j+1][i];
    }else{
      BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[k][2*JEND-j][i];
    }
  }else if (side == X3_BEG){

    if (box->vpos != X3FACE){
      BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[2*KBEG-k-1][j][i];
    }else{
      BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[2*KBEG-k-2][j][i];
    }

  }else if (side == X3_END){

    if (box->vpos != X3FACE){
      BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[2*KEND-k+1][j][i];
    }else{
      BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[2*KEND-k][j][i];
    }

  }
}

/* ********************************************************************* */
void PeriodicBound (double ***q, int side, int vpos)
/*!
 * Implements periodic boundary conditions in serial mode.
 *
 *********************************************************************** */
{
  int  i, j, k;
  RBox *box = GetRBox(side, vpos);

  if (side == X1_BEG){

    BOX_LOOP(box,k,j,i) q[k][j][i] = q[k][j][i + NX1];

  }else if (side == X1_END){

    BOX_LOOP(box,k,j,i) q[k][j][i] = q[k][j][i - NX1];
    
  }else if (side == X2_BEG){

    BOX_LOOP(box,k,j,i) q[k][j][i] = q[k][j + NX2][i];

  }else if (side == X2_END){

    BOX_LOOP(box,k,j,i) q[k][j][i] = q[k][j - NX2][i];

  }else if (side == X3_BEG){

    BOX_LOOP(box,k,j,i) q[k][j][i] = q[k + NX3][j][i];

  }else if (side == X3_END){

    BOX_LOOP(box,k,j,i) q[k][j][i] = q[k - NX3][j][i];

  }
}

