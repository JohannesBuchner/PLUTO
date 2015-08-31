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
void PatchPluto::saveFluxes (const State_1D *state, int beg, int end, 
                             Grid *grid)
/*! 
 *  Rebuild fluxes in a way suitable for AMR operation
 *  by adding pressure and multiplying by area.
 *
 *********************************************************************** */
{
  int  i, nv;
  double **f, *p, r, area;

  f = state->flux;
  p = state->press;

  for (i = beg; i <= end; i++) f[i][MXn] += p[i];

  #if (GEOMETRY == CARTESIAN) && (CH_SPACEDIM > 1)
    if ((g_dir == IDIR) && (g_stretch_fact != 1.)) {
      for (i = beg; i <= end; i++) {
        VAR_LOOP(nv) f[i][nv] *= g_stretch_fact;
      }
    }
   #if (CH_SPACEDIM == 3)
    if ((g_dir == JDIR) && (g_x3stretch != 1.)) {
      for (i = beg; i <= end; i++) {
        VAR_LOOP(nv) f[i][nv] *= g_x3stretch;
      }
    }
    if ((g_dir == KDIR) && (g_x2stretch != 1.)) {
      for (i = beg; i <= end; i++) {
        VAR_LOOP(nv) f[i][nv] *= g_x2stretch;
      }
    }   
   #endif
  #endif

  #if GEOMETRY == CYLINDRICAL
   if (g_dir == IDIR){
     for (i = beg; i <= end; i++) {
     for (nv = 0; nv < NVAR; nv++) {
       f[i][nv] *= grid[IDIR].A[i];
       #if CH_SPACEDIM > 1
        f[i][nv] *= g_x2stretch;
       #endif
     }}
   }else{
     area = fabs(grid[IDIR].x[g_i]);
     for (i = beg; i <= end; i++) {
     for (nv = 0; nv < NVAR; nv++) {
       f[i][nv] *= area;
     }}
   }
  #endif

  #if GEOMETRY == SPHERICAL
   if (g_dir == IDIR){
    #if CH_SPACEDIM > 1
     area = grid[JDIR].dV[g_j]/m_dx;
     #if CH_SPACEDIM == 3
      area *= g_x3stretch;
     #endif 
    #endif
     for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NVAR; nv++) {
        f[i][nv] *= grid[IDIR].A[i];
       #if CH_SPACEDIM > 1
        f[i][nv] *= area;
       #endif 
      }
      #if (COMPONENTS == 3) && (ENTROPY_SWITCH)
       f[i][iMPHI] *= grid[IDIR].xr[i]*sin(grid[JDIR].x[g_j]);
      #endif
     }
   }
   if (g_dir == JDIR){
     area = fabs(grid[IDIR].x[g_i]);
     #if CHOMBO_LOGR == YES
      area *= grid[IDIR].dx[g_i]/m_dx;
     #endif
     #if CH_SPACEDIM == 3
      area *= g_x3stretch;
     #endif
     for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NVAR; nv++) {
        f[i][nv] *= grid[JDIR].A[i]*area;
      }
      #if (COMPONENTS == 3) && (ENTROPY_SWITCH)
       f[i][iMPHI] *= grid[IDIR].x[g_i]*sin(grid[JDIR].xr[i]);
      #endif
     }
   }
   if (g_dir == KDIR){
     area = g_x2stretch*fabs(grid[IDIR].x[g_i]);
    #if CHOMBO_LOGR == YES 
     area *= grid[IDIR].dx[g_i]/m_dx;
    #endif
     for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NVAR; nv++) {
        f[i][nv] *= area;
      }
      #if (COMPONENTS == 3) && (ENTROPY_SWITCH)
       f[i][iMPHI] *= grid[IDIR].x[g_i]*sin(grid[JDIR].x[g_j]);
      #endif
     }
   }
  #endif

  #if GEOMETRY == POLAR
   if (g_dir == IDIR){
     for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NVAR; nv++) {
        f[i][nv] *= grid[IDIR].A[i];
       #if CH_SPACEDIM > 1
        f[i][nv] *= g_x2stretch;
       #endif
       #if CH_SPACEDIM == 3
        f[i][nv] *= g_x3stretch;
       #endif
      }
      #if (COMPONENTS > 1) && (ENTROPY_SWITCH)
       f[i][iMPHI] *= grid[IDIR].xr[i];
      #endif
     }
   }
   if (g_dir == JDIR) {
     area = g_x3stretch;
     #if CHOMBO_LOGR == YES
      area *= grid[IDIR].dx[g_i]/m_dx;
     #endif
     if (area != 1.) {
      for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NVAR; nv++) {
        f[i][nv] *= area;
     }}}
     #if (COMPONENTS > 1) && (ENTROPY_SWITCH)
      for (i = beg; i <= end; i++) f[i][iMPHI] *= grid[IDIR].x[g_i];
     #endif
   }
   if (g_dir == KDIR) {
     area = g_x2stretch*fabs(grid[IDIR].x[g_i]);
    #if CHOMBO_LOGR == YES
     area *= grid[IDIR].dx[g_i]/m_dx;
    #endif
     for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NVAR; nv++) {
        f[i][nv] *= area;
      }
      #if (COMPONENTS > 1) && (ENTROPY_SWITCH)
       f[i][iMPHI] *= grid[IDIR].x[g_i];
      #endif
     }
   }
  #endif
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
  int i,j,k;
  int dir, err;
  int nx, ny, nz;
  int lft_side[3] = {0,0,0}, rgt_side[3]={0,0,0};
  static unsigned char ***flagEntr;
  RBox cbox, *box;

  nx = grid[IDIR].np_tot;
  ny = grid[JDIR].np_tot;
  nz = grid[KDIR].np_tot;

/* ------------------------------------------------------- 
     Check whether the patch touches a physical boundary
   ------------------------------------------------------- */

  for (dir = 0; dir < DIMENSIONS; dir++){
    lft_side[dir] = (grid[dir - IDIR].lbound != 0);
    rgt_side[dir] = (grid[dir - IDIR].rbound != 0);
  }

/* ---------------------------------------------------
    Extract the portion of the domain where U 
    is defined (i.e. NOT in the physical boundary).
   --------------------------------------------------- */

  cbox.ib = 0; cbox.ie = nx - 1;
  cbox.jb = 0; cbox.je = ny - 1;
  cbox.kb = 0; cbox.ke = nz - 1;

/* -------------------------------------------------
    Exclude physical boundaries since the 
    conservative vector U is not yet defined.    
   ------------------------------------------------- */

  D_EXPAND(if (lft_side[IDIR]) cbox.ib = IBEG;  ,
           if (lft_side[JDIR]) cbox.jb = JBEG;  ,
           if (lft_side[KDIR]) cbox.kb = KBEG;)

  D_EXPAND(if (rgt_side[IDIR]) cbox.ie = IEND;  ,
           if (rgt_side[JDIR]) cbox.je = JEND;  ,
           if (rgt_side[KDIR]) cbox.ke = KEND;)

/* ----------------------------------------------------------
    Convert conservative variables into primitive variables.
    Normally this operation is performed by using total
    energy density.
    However, when the ENTROPY_SWITCH is enabled, we force
    conversion to be done from the entropy by artificially
    setting flagEntr.
   ---------------------------------------------------------- */

#if ENTROPY_SWITCH
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
  ConsToPrim3D(U, d->Vc, flagEntr, &cbox);
#else
  ConsToPrim3D(U, d->Vc, d->flag, &cbox);
#endif
  Boundary (d, ALL_DIR, grid);

/* --------------------------------------------------------------
    Convert primitive variables to conservative in the ghost 
    zones.
   -------------------------------------------------------------- */

#if INTERNAL_BOUNDARY == YES
  box = GetRBox(TOT, CENTER);
  PrimToCons3D(d->Vc, U, box);
#else
  if (lft_side[IDIR]) {
    box = GetRBox(X1_BEG, CENTER);
    PrimToCons3D(d->Vc, U, box);
  }

  if (lft_side[JDIR]) {
    box = GetRBox(X2_BEG, CENTER);
    PrimToCons3D(d->Vc, U, box);
  }

  if (lft_side[KDIR]) {
    box = GetRBox(X3_BEG, CENTER);
    PrimToCons3D(d->Vc, U, box);
  }

  if (rgt_side[IDIR]) {
    box = GetRBox(X1_END, CENTER);
    PrimToCons3D(d->Vc, U, box);
  }

  if (rgt_side[JDIR]) {
    box = GetRBox(X2_END, CENTER);
    PrimToCons3D(d->Vc, U, box);
  }

  if (rgt_side[KDIR]) {
    box = GetRBox(X3_END, CENTER);
    PrimToCons3D(d->Vc, U, box);
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

  Grid *Gx, *Gy, *Gz;

  D_EXPAND(Gx = grid + IDIR;  ,
           Gy = grid + JDIR;  ,
           Gz = grid + KDIR;)

  print ("+-----------------------------------------------------------+\n");
  print ("| Level = %d \n", grid[IDIR].level);
  print ("| ib,ie = [%d, %d], xb,xe = %s[%f, %f]%s\n", 
             IBEG, IEND, (Gx->lbound != 0 ? pb:ib), 
              Gx->xr[IBEG-1], Gx->xr[IEND],(Gx->rbound != 0 ? pb:ib));
  print ("| jb,je = [%d, %d], yb,ye = %s[%f, %f]%s\n", 
            JBEG, JEND, (Gy->lbound != 0 ? pb:ib), 
            Gy->xr[JBEG-1], Gy->xr[JEND], (Gy->rbound != 0 ? pb:ib));
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

  box.ib = ibeg+IOFFSET; box.ie = iend-IOFFSET;
  box.jb = jbeg+JOFFSET; box.je = jend-JOFFSET;
  box.kb = kbeg+KOFFSET; box.ke = kend-KOFFSET;

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

