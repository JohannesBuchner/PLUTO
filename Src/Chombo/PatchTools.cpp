/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Contains some implementation of the PatchPluto class.

   PatchTools.cpp contains a number of function used by PatchUnsplit.
   Most of them are simple wrappers to other functions employed by the 
   static version of the code.

  \authors C. Zanni   (zanni@oato.inaf.it)\n
           A. Mignone (mignone@ph.unito.it)
  \date    June 25, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"

/* ********************************************************************* */
void PatchPluto::saveFluxes (double **aflux, Grid *grid)
/*! 
 *  Rebuild fluxes in a way suitable for AMR operation
 *  by and multiplying by area.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  int   nxf, nyf, nzf;
  int   nxb, nyb, nzb;
  long int indf, ind1;

  RBox  fluxBox;
  double A;
  double mdx_dim = D_EXPAND(1.0, *m_dx, *m_dx);
  double ***Ax1 = grid->A[IDIR], *x1 = grid->x[IDIR];
  double ***Ax2 = grid->A[JDIR], *x2 = grid->x[JDIR];
  double ***Ax3 = grid->A[KDIR], *x3 = grid->x[KDIR];

#if GEOMETRY == CARTESIAN
  if (g_stretch_fact == 1.0) return;
#endif

/* ------------------------------------------------------------
   1. Compute flux*area for the X1 direction
   ------------------------------------------------------------ */
  g_dir = IDIR;
  nxf = grid->np_int[IDIR] + (g_dir == IDIR);
  nyf = grid->np_int[JDIR] + (g_dir == JDIR);
  nzf = grid->np_int[KDIR] + (g_dir == KDIR);

  nxb = grid->lbeg[IDIR] - (g_dir == IDIR);
  nyb = grid->lbeg[JDIR] - (g_dir == JDIR);
  nzb = grid->lbeg[KDIR] - (g_dir == KDIR);
    
  RBoxDefine (IBEG-1, IEND, JBEG, JEND, KBEG, KEND, CENTER, &fluxBox);

  BOX_LOOP(&fluxBox,k,j,i){
    ind1 = (k - nzb)*nyf*nxf + (j - nyb)*nxf + (i - nxb);
    A = Ax1[k][j][i]/mdx_dim;

    for (nv = 0; nv < NVAR; nv++){
      indf = nv*nzf*nyf*nxf + ind1;
      aflux[IDIR][indf] *= A;
    }

  /* -- Modify flux for angular momentum conservation -- */

    #if GEOMETRY == POLAR && (COMPONENTS > 1) && (ENTROPY_SWITCH)
    indf = (iMPHI)*nzf*nyf*nxf + ind1;
    aflux[IDIR][indf] *= grid->xr[IDIR][i];
    #elif (GEOMETRY == SPHERICAL) && (COMPONENTS == 3) && (ENTROPY_SWITCH)
    indf = (iMPHI)*nzf*nyf*nxf + ind1;
    aflux[IDIR][indf] *= grid->xr[IDIR][i]*sin(x2[j]);
    #endif
  }

/* ------------------------------------------------------------
   2. Compute flux*area for the X2 direction
   ------------------------------------------------------------ */

  g_dir = JDIR;
  nxf = grid->np_int[IDIR] + (g_dir == IDIR);
  nyf = grid->np_int[JDIR] + (g_dir == JDIR);
  nzf = grid->np_int[KDIR] + (g_dir == KDIR);

  nxb = grid->lbeg[IDIR] - (g_dir == IDIR);
  nyb = grid->lbeg[JDIR] - (g_dir == JDIR);
  nzb = grid->lbeg[KDIR] - (g_dir == KDIR);
    
  RBoxDefine (IBEG, IEND, JBEG-1, JEND, KBEG, KEND, CENTER, &fluxBox);

  BOX_LOOP(&fluxBox,k,j,i){
    ind1 = (k - nzb)*nyf*nxf + (j - nyb)*nxf + (i - nxb);

    A = Ax2[k][j][i]/mdx_dim;
    for (nv = 0; nv < NVAR; nv++){
      indf = nv*nzf*nyf*nxf + ind1;
      aflux[JDIR][indf] *= A;
    }

  /* -- Modify flux for angular momentum conservation -- */

    #if GEOMETRY == POLAR && (COMPONENTS > 1) && (ENTROPY_SWITCH)
    indf = (iMPHI)*nzf*nyf*nxf + ind1;
    aflux[JDIR][indf] *= x1[i];
    #elif (GEOMETRY == SPHERICAL) && (COMPONENTS == 3) && (ENTROPY_SWITCH)
    indf = (iMPHI)*nzf*nyf*nxf + ind1;
    aflux[JDIR][indf] *= x1[i]*sin(grid->xr[JDIR][j]);
    #endif
  }
  
/* ------------------------------------------------------------
   3. Compute flux*area for the X3 direction
   ------------------------------------------------------------ */

#if DIMENSIONS == 3
  g_dir = KDIR;
  nxf = grid->np_int[IDIR] + (g_dir == IDIR);
  nyf = grid->np_int[JDIR] + (g_dir == JDIR);
  nzf = grid->np_int[KDIR] + (g_dir == KDIR);

  nxb = grid->lbeg[IDIR] - (g_dir == IDIR);
  nyb = grid->lbeg[JDIR] - (g_dir == JDIR);
  nzb = grid->lbeg[KDIR] - (g_dir == KDIR);
    
  RBoxDefine (IBEG, IEND, JBEG, JEND, KBEG-1, KEND, CENTER, &fluxBox);

  BOX_LOOP(&fluxBox,k,j,i){
    ind1 = (k - nzb)*nyf*nxf + (j - nyb)*nxf + (i - nxb);
    A = Ax3[k][j][i]/mdx_dim;
    for (nv = 0; nv < NVAR; nv++){
      indf = nv*nzf*nyf*nxf + ind1;
      aflux[KDIR][indf] *= A;
    }
    #if GEOMETRY == POLAR && (COMPONENTS > 1) && (ENTROPY_SWITCH)
    indf = (iMPHI)*nzf*nyf*nxf + ind1;
    aflux[KDIR][indf] *= x1[i]*sin(x2[j]);
    #elif (GEOMETRY == SPHERICAL) && (COMPONENTS == 3) && (ENTROPY_SWITCH)
    indf = (iMPHI)*nzf*nyf*nxf + ind1;
    aflux[KDIR][indf] *= x1[i]*sin(grid->xr[JDIR][j]);
    #endif
  }
#endif /* DIMENSIONS == 3 */  
}

/* ********************************************************************* */
void PatchPluto::getPrimitiveVars (Data_Arr U, Data *d, Grid *grid)
/*!
 * - Recover primitive variables from the input conservative array \c U
 * - Set physical boundary conditions and convert
 *
 * \date June 25, 2015
 *********************************************************************** */
{
  int   i,j,k;
  int   dir, err;
  int   nx, ny, nz;
  int   lft_side[3] = {0,0,0}, rgt_side[3]={0,0,0};
  static unsigned char ***flagEntr;
  RBox  tbox, box; 

  nx = grid->np_tot[IDIR];
  ny = grid->np_tot[JDIR];
  nz = grid->np_tot[KDIR];

/* ------------------------------------------------------- 
     Check whether the patch touches a physical boundary
   ------------------------------------------------------- */

  for (dir = 0; dir < DIMENSIONS; dir++){
    lft_side[dir] = (grid->lbound[dir - IDIR] != 0);
    rgt_side[dir] = (grid->rbound[dir - IDIR] != 0);
  }

/* ----------------------------------------------------------
    Extract the portion of the domain where U 
    is defined (i.e. NOT in the physical boundary).
    Here tbox is the  "total box" defining the whole 
    computational patch, while 'box' is defined from 
    time to time.
   --------------------------------------------------------- */

  RBoxDefine (0, nx-1, 0, ny-1, 0, nz-1, CENTER, &tbox);
  box = tbox;

/* -------------------------------------------------
    Exclude physical boundaries since the 
    conservative vector U is not yet defined.    
   ------------------------------------------------- */

  D_EXPAND(if (lft_side[IDIR]) box.ibeg = IBEG;  ,
           if (lft_side[JDIR]) box.jbeg = JBEG;  ,
           if (lft_side[KDIR]) box.kbeg = KBEG;)

  D_EXPAND(if (rgt_side[IDIR]) box.iend = IEND;  ,
           if (rgt_side[JDIR]) box.jend = JEND;  ,
           if (rgt_side[KDIR]) box.kend = KEND;)

/* ----------------------------------------------------------
    Convert conservative variables into primitive variables.
    Normally this operation is performed by using total
    energy density.
    However, when the ENTROPY_SWITCH is enabled, we force
    conversion to be done from the entropy by artificially
    setting flagEntr.
   ---------------------------------------------------------- */

#if    ENTROPY_SWITCH == SELECTIVE \
    || ENTROPY_SWITCH == ALWAYS    \
    || ENTROPY_SWITCH == CHOMBO_REGRID
  if (flagEntr == NULL) {
    flagEntr = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, unsigned char);
    for (k = 0; k < NX3_MAX; k++){
    for (j = 0; j < NX2_MAX; j++){
    for (i = 0; i < NX1_MAX; i++){
      flagEntr[k][j][i] = 0;
      flagEntr[k][j][i] |= FLAG_ENTROPY;
    }}}
  }
/*
  BOX_LOOP(&cbox,k,j,i){
    flagEntr[k][j][i]  = 0;
    flagEntr[k][j][i] |= FLAG_ENTROPY;
  }
*/
  ConsToPrim3D(U, d->Vc, flagEntr, &box);
#else
  ConsToPrim3D(U, d->Vc, d->flag, &box);
#endif
  Boundary (d, ALL_DIR, grid);

/* --------------------------------------------------------------
    Convert primitive variables to conservative in the ghost 
    zones.
   -------------------------------------------------------------- */

#if INTERNAL_BOUNDARY == YES
  PrimToCons3D(d->Vc, U, &tbox);
#else
  if (lft_side[IDIR]) {
    box      = tbox;
    box.iend = IBEG-1;
    PrimToCons3D(d->Vc, U, &box);
  }

  if (lft_side[JDIR]) {
    box      = tbox;
    box.jend = JBEG-1;
    PrimToCons3D(d->Vc, U, &box);
  }

  if (lft_side[KDIR]) {
    box      = tbox;
    box.kend = KBEG-1;
    PrimToCons3D(d->Vc, U, &box);
  }

  if (rgt_side[IDIR]) {
    box      = tbox;
    box.ibeg = IEND+1;
    PrimToCons3D(d->Vc, U, &box);
  }

  if (rgt_side[JDIR]) {
    box      = tbox;
    box.jbeg = JEND+1;
    PrimToCons3D(d->Vc, U, &box);
  }

  if (rgt_side[KDIR]) {
    box      = tbox;
    box.kbeg = KEND+1;
    PrimToCons3D(d->Vc, U, &box);
  }

#endif
}

/* ************************************************************ */
void PatchPluto::showPatch (Grid *grid)
/*
 *
 *
 *
 ************************************************************** */
{
  char pb[4]="*";  /* -- physical boundary -- */
  char ib[4]="o";  /* -- internal boundary -- */

  print ("+-----------------------------------------------------------+\n");
  print ("| Level = %d \n", grid->level);
  print ("| ib,ie = [%d, %d], xb,xe = %s[%f, %f]%s\n", 
             IBEG, IEND, (grid->lbound[IDIR] != 0 ? pb:ib), 
              grid->xr[IDIR][IBEG-1], grid->xr[IDIR][IEND],
             (grid->rbound[IDIR] != 0 ? pb:ib));
  print ("| jb,je = [%d, %d], yb,ye = %s[%f, %f]%s\n", 
            JBEG, JEND, (grid->lbound[JDIR] != 0 ? pb:ib), 
            grid->xr[JDIR][JBEG-1], grid->xr[JDIR][JEND],
            (grid->rbound[JDIR] != 0 ? pb:ib));
  print ("+-----------------------------------------------------------+\n");

}

/* ********************************************************************* */
void PatchPluto::convertFArrayBox(FArrayBox&  U)
/*!
 *  Convert a conservative array to primitive for 
 *  plotfile data. 
 *  Called from AMRLevelPluto::writePlotLevel
 *
 *
 *********************************************************************** */
{
  int ibeg, jbeg, kbeg;
  int iend, jend, kend;
  int i,j,k, nv;
  double ***UU[NVAR];
  static unsigned char *flag;
  static double **u, **v;
  RBox box;
  
  if (u == NULL){
    u    = ARRAY_2D(NMAX_POINT, NVAR, double);
    v    = ARRAY_2D(NMAX_POINT, NVAR, double);
    flag = ARRAY_1D(NMAX_POINT, unsigned char);
  }

  jbeg = jend = kbeg = kend = 0;
  D_EXPAND(ibeg = U.loVect()[IDIR]; 
           iend = U.hiVect()[IDIR]; ,
           jbeg = U.loVect()[JDIR];   
           jend = U.hiVect()[JDIR]; ,
           kbeg = U.loVect()[KDIR];   
           kend = U.hiVect()[KDIR]; );

  for (nv=0; nv<NVAR; nv++) {
    UU[nv] = ArrayBoxMap(kbeg,kend,jbeg,jend,ibeg,iend,U.dataPtr(nv));
  }

  box.ibeg = ibeg+IOFFSET; box.iend = iend-IOFFSET;
  box.jbeg = jbeg+JOFFSET; box.jend = jend-JOFFSET;
  box.kbeg = kbeg+KOFFSET; box.kend = kend-KOFFSET;

/* --------------------------------------------------------
    Conversion is done in the interior points only. 
    For this reason we exclude from (ibeg, iend, ...) one 
    boundary zone in each direction.
   -------------------------------------------------------- */

  for (k = kbeg + KOFFSET; k <= kend - KOFFSET; k++){
  for (j = jbeg + JOFFSET; j <= jend - JOFFSET; j++){
    for (i = ibeg + IOFFSET; i <= iend - IOFFSET; i++){
      flag[i-ibeg] = 0;
      NVAR_LOOP(nv) u[i-ibeg][nv] = UU[nv][k][j][i];
    }
    ConsToPrim (u, v, IOFFSET, iend-ibeg-IOFFSET, flag);
    for (i = ibeg; i <= iend; i++){
      NVAR_LOOP(nv) UU[nv][k][j][i] = v[i-ibeg][nv];
    }
  }}

  for (nv = 0; nv < NVAR; nv++) 
    FreeArrayBoxMap(UU[nv], kbeg, kend, jbeg, jend, ibeg, iend);
}

