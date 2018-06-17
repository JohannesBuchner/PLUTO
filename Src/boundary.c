/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Set boundary conditions.

  The Boundary() function sets both internal and physical boundary 
  conditions on one or more sides of the computational domain.
  It is used to fill ghost zones of both cell-centered and face-centered 
  data arrays.\n
  The type of boundary conditions at the leftmost or rightmost side of a 
  given grid is specified by the integers <c> grid.lbound[dir] </c> or 
  <c> grid.rbound[dir] </c>, respectively.
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
  \date   Nov 23, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

static void OutflowBound   (double ***, RBox *, int);
static void ReflectiveBound(double ***, int, RBox *, int);
static void PeriodicBound  (double ***, RBox *, int);
static void FlipSign       (int, int, int *);
                           
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
 * \param [in]  grid   pointer to grid structure.
 *********************************************************************** */
{
  int  is, nv;
  int  side[6] = {X1_BEG, X1_END, X2_BEG, X2_END, X3_BEG, X3_END};
  int  type[6], sbeg, send, vsign[NVAR];
  int  par_dim[3] = {0, 0, 0};
  int  ib,ie,jb,je,kb,ke;
  RBox center_box, x1face_box, x2face_box, x3face_box;

/* --------------------------------------------------------
   0. Check the number of processors in each direction
   -------------------------------------------------------- */

  D_EXPAND(par_dim[0] = grid->nproc[IDIR] > 1;  ,
           par_dim[1] = grid->nproc[JDIR] > 1;  ,
           par_dim[2] = grid->nproc[KDIR] > 1;)

/* --------------------------------------------------------
   1. When FARGO is on, boundary conditions must be set 
      on total velocity. 
   -------------------------------------------------------- */

#ifdef FARGO
  if (g_hydroStep && FARGO_TotalVelocityIsSet() == 0){
    FARGO_AddVelocity (d,grid);
  }
#endif

/* --------------------------------------------------------
   2. Call userdef internal boundary with side == 0
   --------------------------------------------------------  */

#if INTERNAL_BOUNDARY == YES
  UserDefBoundary (d, NULL, 0, grid);
#endif
  
/* --------------------------------------------------------
   3.  Exchange data between processors 
   -------------------------------------------------------- */
   
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

/* ---------------------------------------------------------
   4. When idim == ALL_DIR boundaries are imposed on ALL 
      sides: a loop from sbeg = 0 to send = 2*DIMENSIONS - 1 
      is performed. 
     
      When idim = n, boundaries are imposed at the 
      beginning and the end of the i-th direction.
   -------------------------------------------------------- */ 

  if (idim == ALL_DIR) {
    sbeg = 0;
    send = 2*DIMENSIONS - 1;
  } else {
    sbeg = 2*idim;
    send = 2*idim + 1;
  }

/* --------------------------------------------------------
   5. Main loop on computational domain sides
   -------------------------------------------------------- */

  type[0] = grid->lbound[IDIR]; type[1] = grid->rbound[IDIR];
  type[2] = grid->lbound[JDIR]; type[3] = grid->rbound[JDIR];
  type[4] = grid->lbound[KDIR]; type[5] = grid->rbound[KDIR];

  for (is = sbeg; is <= send; is++){

    if (type[is] == 0) continue;  /* no physical boundary: skip */

  /* ----------------------------------------------
     5a. Define boundary boxes, sweeping direction. 
     ---------------------------------------------- */

    ib = 0; ie = NX1_TOT-1;
    jb = 0; je = NX2_TOT-1;
    kb = 0; ke = NX3_TOT-1;

    if      (side[is] == X1_BEG) {ib = IBEG-1; ie = 0;}         /* Backward */
    else if (side[is] == X1_END) {ib = IEND+1; ie = NX1_TOT-1;} /* Forward */
    else if (side[is] == X2_BEG) {jb = JBEG-1; je = 0;}         /* Backward */
    else if (side[is] == X2_END) {jb = JEND+1; je = NX2_TOT-1;} /* Forward */
    else if (side[is] == X3_BEG) {kb = KBEG-1; ke = 0;}         /* Backward */
    else if (side[is] == X3_END) {kb = KEND+1; ke = NX3_TOT-1;} /* Forward */

    RBoxDefine (ib, ie, jb, je, kb, ke, CENTER, &center_box);
#ifdef STAGGERED_MHD
    /* -- Define RBoxes for staggered field (! Note the backward/forward order) -- */
    if (ib < ie) RBoxDefine (ib-1, ie  , jb  , je  , kb  , ke  , X1FACE, &x1face_box);
    else         RBoxDefine (ib  , ie-1, jb  , je  , kb  , ke  , X1FACE, &x1face_box);
    if (jb < je) RBoxDefine (ib  , ie  , jb-1, je  , kb  , ke  , X2FACE, &x2face_box);
    else         RBoxDefine (ib  , ie  , jb  , je-1, kb  , ke  , X2FACE, &x2face_box);
    if (kb < ke) RBoxDefine (ib  , ie  , jb  , je  , kb-1, ke  , X3FACE, &x3face_box);
    else         RBoxDefine (ib  , ie  , jb  , je  , kb  , ke-1, X3FACE, &x3face_box);
#endif

  /* --------------------------------------------------------
     5b. Apply boundary conditions.
     -------------------------------------------------------- */

    if (type[is] == OUTFLOW) {    /* ---- Outflow B.C. ---- */

      NVAR_LOOP(nv) OutflowBound (d->Vc[nv], &center_box, side[is]);
#ifdef STAGGERED_MHD
      D_EXPAND(if (side[is] != X1_BEG && side[is] != X1_END)
                 OutflowBound (d->Vs[BX1s], &x1face_box, side[is]); ,
               if (side[is] != X2_BEG && side[is] != X2_END)
                 OutflowBound (d->Vs[BX2s], &x2face_box, side[is]); ,
               if (side[is] != X3_BEG && side[is] != X3_END)                 
                 OutflowBound (d->Vs[BX3s], &x3face_box, side[is]);)

    /* ---------------------------------------------------------
        Recover normal magnetic field and average it to cell
        center. Transverse components are assigned consistently
        with cell-centered quantities.
       --------------------------------------------------------- */

      FillMagneticField (d, side[is], grid);
      CT_AverageNormalMagField (d, side[is], grid);
#endif

    }else if (  (type[is] == REFLECTIVE) 
             || (type[is] == AXISYMMETRIC)
             || (type[is] == EQTSYMMETRIC)){ /* ---- Reflective B.C. ---- */

      FlipSign (side[is], type[is], vsign);
      NVAR_LOOP(nv) ReflectiveBound (d->Vc[nv], vsign[nv], &center_box, side[is]);
#ifdef STAGGERED_MHD
      D_EXPAND(if (side[is] != X1_BEG && side[is] != X1_END)
                 ReflectiveBound(d->Vs[BX1s], vsign[BX1], &x1face_box, side[is]);  ,
               if (side[is] != X2_BEG && side[is] != X2_END)  
                 ReflectiveBound(d->Vs[BX2s], vsign[BX2], &x2face_box, side[is]);  ,
               if (side[is] != X3_BEG && side[is] != X3_END)
                 ReflectiveBound(d->Vs[BX3s], vsign[BX3], &x3face_box, side[is]);)
                  
      FillMagneticField (d, side[is], grid);
#endif

    }else if (type[is] == PERIODIC){  /* -- Periodic B.C. (serial or 1 proc) -- */

      if (!par_dim[is/2]) {
        NVAR_LOOP(nv)  PeriodicBound(d->Vc[nv], &center_box, side[is]);
#ifdef STAGGERED_MHD
        D_EXPAND(PeriodicBound(d->Vs[BX1s], &x1face_box, side[is]);  ,
                 PeriodicBound(d->Vs[BX2s], &x2face_box, side[is]);  ,
                 PeriodicBound(d->Vs[BX3s], &x3face_box, side[is]);)
#endif
      }

    }else if (type[is] == SHEARING) {  /* -- Shearingbox B.C. -- */

    /* ---------------------------------------------------------
         SHEARING-BOX boundary condition is implemented as

        1) apply periodic boundary conditions for all variables
          (except staggered BX)
        2) Perform spatial shift in the y-direction
       --------------------------------------------------------- */

#ifdef SHEARINGBOX
      if (side[is] != X1_BEG && side[is] != X1_END){
        print ("! Boundary(): shearingbox can only be assigned at an X1 boundary\n");
        QUIT_PLUTO(1);
      }
      if (grid->nproc[IDIR] == 1){
        NVAR_LOOP(nv) PeriodicBound(d->Vc[nv], &center_box, side[is]);
      #ifdef STAGGERED_MHD
        D_EXPAND(                                                 ;  ,
                 PeriodicBound(d->Vs[BX2s], &x2face_box, side[is]);  ,
                 PeriodicBound(d->Vs[BX3s], &x3face_box, side[is]);)
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
#else
      print ("! Boundary(): shearingbox module not loaded\n");
      QUIT_PLUTO(1);
#endif  /* #ifdef SHEARINGBOX */

    }else if (type[is] == USERDEF) { /* ---- User-defined B.C. ---- */

      #ifdef GLM_MHD
      {int i,j,k;
/*
       switch(side[is]){
         case X1_BEG: X1_BEG_LOOP(k,j,i) d->Vc[PSI_GLM][k][j][i] = 0.0; break;
         case X1_END: X1_END_LOOP(k,j,i) d->Vc[PSI_GLM][k][j][i] = 0.0; break;
         case X2_BEG: X2_BEG_LOOP(k,j,i) d->Vc[PSI_GLM][k][j][i] = 0.0; break;
         case X2_END: X2_END_LOOP(k,j,i) d->Vc[PSI_GLM][k][j][i] = 0.0; break;
         case X3_BEG: X3_BEG_LOOP(k,j,i) d->Vc[PSI_GLM][k][j][i] = 0.0; break;
         case X3_END: X3_END_LOOP(k,j,i) d->Vc[PSI_GLM][k][j][i] = 0.0; break;
       }
*/
       BOX_LOOP(&center_box, k, j, i) d->Vc[PSI_GLM][k][j][i] = 0.0;
       #ifdef PHI_GLM
       BOX_LOOP(&center_box, k, j, i) d->Vc[PHI_GLM][k][j][i] = 0.0;
       #endif
      }
      #endif

      UserDefBoundary (d, &center_box, side[is], grid);
#ifdef STAGGERED_MHD
      D_EXPAND(UserDefBoundary (d, &x1face_box, side[is], grid);  , /* Only 2 needed */
               UserDefBoundary (d, &x2face_box, side[is], grid);  ,
               UserDefBoundary (d, &x3face_box, side[is], grid);)

       /* -- assign normal component of staggered B 
             using the div.B = 0 condition           -- */

      FillMagneticField (d, side[is], grid);  
      CT_AverageNormalMagField (d, side[is], grid);
      CT_AverageTransverseMagField (d, side[is], grid);
#endif
    }
  }

/* --------------------------------------------------------
   6. Compute entropy for the next time level
   -------------------------------------------------------- */

#if ENTROPY_SWITCH
  ComputeEntropy (d, grid);
#endif

/* --------------------------------------------------------
   7. Subtract velocity if we're doing an hydro step
      (but not if we're doing an operator-split parabolic
       step)
   -------------------------------------------------------- */

#ifdef FARGO
  if (g_hydroStep && FARGO_TotalVelocityIsSet() == 1) {
    FARGO_SubtractVelocity(d, grid);
  }  
#endif

}

/* ********************************************************************* */
void FlipSign (int side, int type, int *vsign)
/*!
 * Reverse the sign of vector components with respect to axis side. 
 * Depending on type, one needs to symmetrize or anti-symmetrize:
 *
 * - REFLECTIVE:    \n
 *   o   Vn -> -Vn,  Bn -> -Bn,  En ->  En  \n
 *   o   Vp ->  Vp,  Bp ->  Bp,  Ep -> -Ep  \n 
 *   o   Vt ->  Vt,  Bt ->  Bt,  Et -> -Et
 *
 * - AXISYMMETRIC:  \n
 *   o   Vn -> -Vn,  Bn -> -Bn,  En -> -En \n
 *   o   Vp ->  Vp,  Bp ->  Bp,  Ep ->  Ep \n
 *   o   Vt -> -Vt,  Bt -> -Bt,  Et -> -Et
 *
 * - EQTSYMMETRIC:  \n 
 *   o   Vn -> -Vn,  Bn ->  Bn,  En -> -En  \n
 *   o   Vp ->  Vp,  Bp -> -Bp,  Ep ->  Ep  \n
 *   o   Vt ->  Vt,  Bt -> -Bt,  Et ->  Et
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
#if PHYSICS == MHD || PHYSICS == RMHD
  int Bn, Bp, Bt;
#endif
#if PHYSICS == RMHD || RESISTIVITY != NO
  int En, Ep, Et;
#endif

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
    #if PHYSICS == RMHD && RESISTIVITY != NO
    En = EX1; Ep = EX2; Et = EX3;
    #endif 
  }else if (side == X2_BEG || side == X2_END){
    Vn = VX2; Vp = VX1; Vt = VX3;
    #if PHYSICS == MHD || PHYSICS == RMHD
    Bn = BX2; Bp = BX1; Bt = BX3;
    #endif
    #if PHYSICS == RMHD && RESISTIVITY != NO
    En = EX2; Ep = EX1; Et = EX3;
    #endif 
  }else if (side == X3_BEG || side == X3_END){
    Vn = VX3; Vp = VX1; Vt = VX2;
    #if PHYSICS == MHD || PHYSICS == RMHD
    Bn = BX3; Bp = BX1; Bt = BX2;
    #endif
    #if PHYSICS == RMHD && RESISTIVITY != NO
    En = EX3; Ep = EX1; Et = EX2;
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
  
  if (type == REFLECTIVE){
    #if PHYSICS == RMHD && RESISTIVITY != NO
    vsign[Ep] = -1.0;
    vsign[Et] = -1.0;
    #endif 
  }  

#if COMPONENTS == 3
  if (type == AXISYMMETRIC){
    vsign[Vt] = -1.0;
    #if PHYSICS == MHD || PHYSICS == RMHD
    vsign[Bt] = -1.0;
    #endif
    #if PHYSICS == RMHD && RESISTIVITY != NO
    vsign[En] = -1.0;
    vsign[Et] = -1.0;
    #endif 
  }
#endif

  if (type == EQTSYMMETRIC){
    #if PHYSICS == MHD || PHYSICS == RMHD
    EXPAND(vsign[Bn] =  1.0; ,
           vsign[Bp] = -1.0; ,
           vsign[Bt] = -1.0;)
    #ifdef GLM_MHD 
    vsign[PSI_GLM] = -1.0;
    #endif
    #if PHYSICS == RMHD && RESISTIVITY != NO
    vsign[En] = -1.0;
    #endif 
    #endif
  }
}

/* ********************************************************************* */
void OutflowBound (double ***q, RBox *box, int side)
/*! 
 * Impose zero-gradient boundary conditions on 'q' on 
 * the boundary side specified by 'side'.
 * The input array 'q' must not represent the normal component
 * of a staggered magnetic fied.
 *
 * \param [in,out] q     a 3D array requiring ghost zone filling
 * \param [in]     box   pointer to a RBox structure defining the
 *                       extent of the boundary region
 * \param [in]     side  the side of the computational domain.
 *********************************************************************** */
{
  int  i, j, k;

  if (side == X1_BEG) { 

    BOX_LOOP(box,k,j,i) q[k][j][i] = q[k][j][IBEG];   

  }else if (side == X1_END){

    BOX_LOOP(box,k,j,i) q[k][j][i] = q[k][j][IEND];   

  }else if (side == X2_BEG){

    BOX_LOOP(box,k,j,i) q[k][j][i] = q[k][JBEG][i];

  }else if (side == X2_END){

    BOX_LOOP(box,k,j,i) q[k][j][i] = q[k][JEND][i];

  }else if (side == X3_BEG){

    BOX_LOOP(box,k,j,i) q[k][j][i] = q[KBEG][j][i];

  }else if (side == X3_END){

    BOX_LOOP(box,k,j,i) q[k][j][i] = q[KEND][j][i];

  }
}

/* ********************************************************************* */
void PeriodicBound (double ***q, RBox *box, int side)
/*!
 * Implements periodic boundary conditions in serial mode or when 
 * one processor only handle the periodic direction.
 *
 * \param [in,out] q     a 3D array requiring ghost zone filling
 * \param [in]     box   pointer to a RBox structure defining the
 *                       extent of the boundary region
 * \param [in]     side  the side of the computational domain.
 *********************************************************************** */
{
  int  i, j, k;

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

/* ********************************************************************* */
void ReflectiveBound (double ***q, int s, RBox *box, int side)
/*!
 * Make symmetric (s = 1) or anti-symmetric (s=-1) profiles 
 * with respect to the boundary plane specified by box->side.
 * The sign is set by the FlipSign() function. 
 *
 * \param [in,out] q     a 3D array requiring ghost zone filling
 * \param [in]     s     an integer taking only the values +1 (symmetric 
 *                       profile) or -1 (antisymmetric profile)
 * \param [in]     box   pointer to a RBox structure defining the
 *                       extent of the boundary region
 * \param [in]     side  the side of the computational domain.
 *********************************************************************** */
{
  int   i, j, k;

  if (side == X1_BEG) {   

    BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[k][j][2*IBEG-i-1];

  }else if (side == X1_END){  

    BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[k][j][2*IEND-i+1];

  }else if (side == X2_BEG){  

    BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[k][2*JBEG-j-1][i];

  }else if (side == X2_END){

    BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[k][2*JEND-j+1][i];

  }else if (side == X3_BEG){

    BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[2*KBEG-k-1][j][i];

  }else if (side == X3_END){

    BOX_LOOP(box,k,j,i) q[k][j][i] = s*q[2*KEND-k+1][j][i];

  }
}

