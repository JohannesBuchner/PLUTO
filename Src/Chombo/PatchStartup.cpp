#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"

void PatchPluto::startup(FArrayBox& a_U)
{
  int i, j, k;
  int isub, jsub, ksub, nsub = 5;
  int nxtot, nytot, nztot;
  int nv,  l_convert;
  int ibg, ieg, jbg, jeg, kbg, keg;
  static double **ucons, **uprim;
  double x1,  x2,  x3;
  double x1s, x2s, x3s;
  double dx1, dx2, dx3;
  double us[256], u_av[256], b[3];
  double scrh, dx, cylr;

  double ***UU[NVAR];
  Grid   grid[3];
  Box curBox = a_U.box();
  RBox tbox;

  FArrayBox dV(curBox,CHOMBO_NDV);

  curBox.grow(-m_numGhost);

  setGrid(curBox,grid,dV);  /* -- define grid -- */

  jbg = jeg = kbg = keg = 0;

  D_EXPAND(ibg = a_U.loVect()[0]; ieg = a_U.hiVect()[0];  ,
           jbg = a_U.loVect()[1]; jeg = a_U.hiVect()[1];  ,
           kbg = a_U.loVect()[2]; keg = a_U.hiVect()[2];)

  NX1_TOT = nxtot = ieg - ibg + 1;
  NX2_TOT = nytot = jeg - jbg + 1;
  NX3_TOT = nztot = keg - kbg + 1;

 if (uprim == NULL){
   uprim = ARRAY_2D(NMAX_POINT, NVAR, double);
   ucons = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

 for (nv=0 ; nv<NVAR ; nv ++)
   UU[nv] = ArrayMap(nztot,nytot,nxtot,a_U.dataPtr(nv));

 /* ----  set labels  ---- */

  EXPAND(MXn = VXn = VX1;  ,
         MXt = VXt = VX2;  ,
         MXb = VXb = VX3;)

  #if PHYSICS == MHD || PHYSICS == RMHD
   EXPAND(BXn = BX1;  ,
          BXt = BX2;  ,
          BXb = BX3;)
  #endif

  tbox.ib = 0; tbox.ie = nxtot-1;
  tbox.jb = 0; tbox.je = nytot-1;
  tbox.kb = 0; tbox.ke = nztot-1;
  
/* --------------------------------------------------------------
                    Assign initial conditions
   -------------------------------------------------------------- */

  BOX_LOOP(&tbox,k,j,i){
    x3 = grid[KDIR].x[k];
    x2 = grid[JDIR].x[j];
    x1 = grid[IDIR].x[i];

    for (nv = 0; nv < NVAR; nv++)  UU[nv][k][j][i] = u_av[nv] = 0.0;
    
#ifdef GLM_MHD
    u_av[PSI_GLM] = us[PSI_GLM] = 0.0;
#endif

/*  ----------------------------------------------------------------
                Compute volume averages
    ---------------------------------------------------------------- */

    #if INITIAL_SMOOTHING == YES

     for (ksub = 0; ksub < nsub; ksub++){
     for (jsub = 0; jsub < nsub; jsub++){
     for (isub = 0; isub < nsub; isub++){

       x1s = x1 + (double)(1.0 - nsub + 2.0*isub)/(double)(2.0*nsub)*dx1;
       x2s = x2 + (double)(1.0 - nsub + 2.0*jsub)/(double)(2.0*nsub)*dx2;
       x3s = x3 + (double)(1.0 - nsub + 2.0*ksub)/(double)(2.0*nsub)*dx3;

       Init (us, x1s, x2s, x3s);
       for (nv = 0; nv < NVAR; nv++) {
         u_av[nv] += us[nv]/(double)(nsub*nsub*nsub);
       }
     }}}

    #else

     Init (u_av, x1, x2, x3);
    
    #endif

    for (nv = 0; nv < NVAR; nv++) UU[nv][k][j][i] = u_av[nv];

    #if (PHYSICS == MHD || PHYSICS == RMHD) 
     #if ASSIGN_VECTOR_POTENTIAL == YES
      VectorPotentialDiff(b, i, j, k, grid);
      for (nv = 0; nv < DIMENSIONS; nv++) UU[BX1+nv][k][j][i] = b[nv];
     #endif  /* ASSIGN_VECTOR_POTENTIAL */
    #endif /* PHYSICS == MHD || PHYSICS == RMHD */

  }

/* --------------------------------------------------------------------
     Convert primitive variables to conservative ones
   -------------------------------------------------------------------- */

  for (k = 0; k < nztot ; k++) {
  for (j = 0; j < nytot ; j++) {

    for (i = 0; i < nxtot; i++) {
      for (nv = 0 ; nv < NVAR ; nv++) {
        us[nv] = uprim[i][nv] = UU[nv][k][j][i];
      }

  /* -- check if primitive values are physically ok -- */

      if (us[RHO] <= 0.0) {
        print ("! startup: density is negative\n");
        QUIT_PLUTO(1);
      }
      #if EOS != ISOTHERMAL
       if (us[PRS] <= 0.0) {
         print ("! startup: pressure is negative\n");
         QUIT_PLUTO(1);
       }
      #endif
      #if (PHYSICS == RHD) || (PHYSICS == RMHD)
       scrh = EXPAND(us[VX1]*us[VX1], + us[VX2]*us[VX2], + us[VX3]*us[VX3]);
       if (scrh >= 1.0){
         print ("! startup: total velocity exceeds 1\n"); 
         QUIT_PLUTO(1);
       }
      #endif
    }

    PrimToCons (uprim, ucons, 0, nxtot-1);
    for (i = 0; i < nxtot; i++) {
    for (nv = 0; nv < NVAR; nv++) {
      UU[nv][k][j][i] = ucons[i][nv];
    }}

#if ENTROPY_SWITCH
    Entropy(uprim, UU[ENTR][k][j], 0, nxtot-1);  /* -- primitive: s -- */
    for (i = 0; i < nxtot; i++) {
      UU[ENTR][k][j][i] *= UU[RHO][k][j][i];   /* -- conservative: s*D -- */
    }
#endif

  }}

/* --------------------------------------------------
     Pass U*dV/m_dx^3 to the library
   -------------------------------------------------- */

  #if GEOMETRY != CARTESIAN
   #if CHOMBO_CONS_AM == YES
    #if ROTATING_FRAME == YES
     Box aBox = a_U.box();
     for(BoxIterator bit(aBox); bit.ok(); ++bit) {
       const IntVect& iv = bit();
       a_U(iv,iMPHI) += a_U(iv,RHO)*dV(iv,1)*g_OmegaZ;
       a_U(iv,iMPHI) *= dV(iv,1);
     }
    #else
     a_U.mult(dV,1,iMPHI);
    #endif
   #endif
   for (nv = 0; nv < a_U.nComp(); nv++) a_U.mult(dV,0,nv);
  #else
   if (g_stretch_fact != 1.) a_U *= g_stretch_fact;
  #endif
 
  for (nv = 0; nv < NVAR; nv++) FreeArrayMap(UU[nv]);
  FreeGrid(grid);
}
