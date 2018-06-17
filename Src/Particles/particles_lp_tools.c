/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Collection of tools for the LAGRANGIAN particle module
        at particle positions.
 
 \authors  B. Vaidya (bvaidya@unito.it)\n
           A. Mignone (mignone@ph.unito.it)
            
 \date   May 31, 2017
 */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if PARTICLES_TYPE == LAGRANGIAN
/* ********************************************************************* */
void Particles_LP_FixValue(Particle *pl, Data *d, Grid *grid) 
/*!
 * This function finds the cell of the particles and fix all the phisical
 * variables, interpolating them on the particle coordinates.
 * Call this routine only after pl->coor have been correctly initialized or 
 * calculated.
 *
 *  \param [in]    pl         Pointer to the Particle.
 *  \param [in]    d          Pointer to the Data structure.
 *  \param [in]    grid       Pointer to the PLUTO grid structure.
 *********************************************************************** */ 
{
  int  i,j,k, dir, indx[3];
  static double ***wts;

  if (wts == NULL) wts = ArrayBox (-1, 1, -1, 1, -1, 1);
  
  Particles_GetWeights(pl, pl->cell, wts, grid);  

  EXPAND(pl->speed[IDIR] = Particles_Interpolate(d->Vc[VX1],  wts, pl->cell); ,
         pl->speed[JDIR] = Particles_Interpolate(d->Vc[VX2],  wts, pl->cell); ,
         pl->speed[KDIR] = Particles_Interpolate(d->Vc[VX3],  wts, pl->cell);)
  
  pl->density = Particles_Interpolate(d->Vc[RHO], wts, pl->cell);

#if PARTICLES_LP_SPECTRA == YES     
  #if HAVE_ENERGY
    pl->pressure = Particles_Interpolate (d->Vc[PRS], wts, pl->cell); 
  #endif
  pl->lorG = 1.0;

  #if PHYSICS == MHD || PHYSICS == RMHD
  double Btot2, Utot2, udotB, sin2alpha;

    /*double lab_density;
    static double ****Uvel; 
    if (Uvel == NULL) Uvel = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
    
    lab_density = Particles_Interpolate(d->Uc[RHO], wts, pl->cell); 
    pl->lorG = lab_density/pl->density;
    pl->fourvel[0] = pl->lorG;
*/
    /* Compute Cell by cell in the fluid grid the Four Velocity from 3 Velocity */
  /*  TOT_LOOP(k,j,i){
      EXPAND(Uvel[IDIR][k][j][i] = pl->lorG*d->Vc[VX1][k][j][i]; , 
             Uvel[JDIR][k][j][i] = pl->lorG*d->Vc[VX2][k][j][i];, 
             Uvel[KDIR][k][j][i] = pl->lorG*d->Vc[VX3][k][j][i]; )
          }
  
    pl->fourvel[1] = Particles_Interpolate(Uvel[IDIR],  wts, pl->cell); 
    pl->fourvel[2] = Particles_Interpolate(Uvel[JDIR],  wts, pl->cell);
    pl->fourvel[3] = Particles_Interpolate(Uvel[KDIR],  wts, pl->cell);
  
   for(dir=0;dir<3;dir++) pl->speed[dir] = pl->fourvel[dir+1]/pl->lorG;
   */
  EXPAND(pl->mag[IDIR] = Particles_Interpolate(d->Vc[BX1], wts, pl->cell); , 
         pl->mag[JDIR] = Particles_Interpolate(d->Vc[BX2], wts, pl->cell); , 
         pl->mag[KDIR] = Particles_Interpolate(d->Vc[BX3], wts, pl->cell);)
    

  Btot2 = EXPAND(pl->mag[IDIR]*pl->mag[IDIR] ,
               + pl->mag[JDIR]*pl->mag[JDIR] ,
               + pl->mag[KDIR]*pl->mag[KDIR]);

  Utot2 = EXPAND(  pl->speed[IDIR]*pl->speed[IDIR] ,
                 + pl->speed[JDIR]*pl->speed[JDIR] ,
                 + pl->speed[KDIR]*pl->speed[KDIR]);
  udotB = EXPAND(  pl->mag[IDIR]*pl->speed[IDIR] ,
                 + pl->mag[JDIR]*pl->speed[JDIR] ,
                 + pl->mag[KDIR]*pl->speed[KDIR]);    
     
  #if PHYSICS == RMHD
  pl->lorG = 1.0/sqrt(1.0 - Utot2);
  Btot2 = Btot2/(pl->lorG*pl->lorG) + udotB*udotB;
  #endif  
   
  /* -- Here the B2 total B field is used as particles are assumed isotropic -- */
  
  double Urad, Umag, zeta, Tcmb = 2.728;  
  Urad = 4.0*(CONST_sigma/CONST_c)*(pow(Tcmb*(1. + PARTICLES_LP_ICCMBZ), 4.0)); 
  Umag = Btot2*UNIT_MAGFIELD*UNIT_MAGFIELD;                                                
  zeta = MIN(Urad/Umag, 5000.); /* This is to avoid in regions with B = 0.0, !! 5000 is arbitary high number ?? !!*/
  pl->cr = SYNCHROTRON_CONST*Btot2*(1.0 + zeta);
                
  #endif /* PHYSICS == MHD || PHYSICS == RMHD */
  
  #if HAVE_ENERGY
  Particles_LP_GradP(pl->gradp, PRS, d->Vc, grid, pl->cell);
  #else
  Particles_LP_GradP(pl->gradp, RHO, d->Vc, grid, pl->cell);    
  #endif
#endif /* PARTICLES_LP_SPECTRA */
}

#if PARTICLES_LP_SPECTRA == YES
/* ********************************************************************* */
void Particles_LP_FlagShock(Data *d, float ***flag2, Grid *grid)
/*!
 * Flag cells that are shocked. The condition used is same as that
 * used for the fluid. i.e., divV < 0.0, |GradP|/P > Threshold.
 *
 * \param [in] d      Pointer to the data structure
 * \param [in] grid   Pointer to the grid structure
 * \param [out] flag2 3D Array of the flagged zones.
 *
 *********************************************************************** */
{
  int  i, j, k, nv;
  int  ip, jp, kp;
  double divv, gradp, pt_min;
  double dpx1, pt_min1, dvx1;
  double dpx2, pt_min2, dvx2;
  double dpx3, pt_min3, dvx3;
  D_EXPAND(double ***vx1 = d->Vc[VX1];  ,
                 double ***vx2 = d->Vc[VX2];  , 
                 double ***vx3 = d->Vc[VX3];)
   
  double *dx1 = grid->dx[IDIR];
  double *dx2 = grid->dx[JDIR];
  double *dx3 = grid->dx[KDIR];
 
  double ***Ax1 = grid->A[IDIR];
  double ***Ax2 = grid->A[JDIR];
  double ***Ax3 = grid->A[KDIR];

  static double ***pt;
  int i1, j1, k1, itag, jtag, ktag;
  int nbuf = 4;

/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */

  if (pt == NULL) pt = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
  
/* --------------------------------------------------------
   1. Compute total pressure and flag, initially all zones
      with ENTROPY_SWITCH.
   -------------------------------------------------------- */

  TOT_LOOP(k,j,i){
    flag2[k][j][i] = 0.0;
    #if EOS == ISOTHERMAL
    pt[k][j][i] = d->Vc[RHO][k][j][i]*g_isoSoundSpeed*g_isoSoundSpeed;
    #else
    #if HAVE_ENERGY 
    pt[k][j][i] = d->Vc[PRS][k][j][i];
    #endif
    #endif
  }
 
/* --------------------------------------------------------
   2. Track zones lying in a shock
   -------------------------------------------------------- */

  for (k = KOFFSET; k < NX3_TOT-KOFFSET; k++){
  for (j = JOFFSET; j < NX2_TOT-JOFFSET; j++){
  for (i = IOFFSET; i < NX1_TOT-IOFFSET; i++){
  
  /* -- Compute divergence of velocity -- */

    #if GEOMETRY == CARTESIAN
    D_EXPAND(dvx1 = (vx1[k][j][i+1] - vx1[k][j][i-1])/dx1[i];   ,
             dvx2 = (vx2[k][j+1][i] - vx2[k][j-1][i])/dx2[j];   ,
             dvx3 = (vx3[k+1][j][i] - vx3[k-1][j][i])/dx3[k];)

    divv = D_EXPAND(dvx1, + dvx2, + dvx3);
    #else
    D_EXPAND(dvx1 =   Ax1[k][j][i]*  (vx1[k][j][i+1] + vx1[k][j][i])
                    - Ax1[k][j][i-1]*(vx1[k][j][i-1] + vx1[k][j][i]);  ,

             dvx2 =   Ax2[k][j][i]*  (vx2[k][j+1][i] + vx2[k][j][i])
                    - Ax2[k][j-1][i]*(vx2[k][j-1][i] + vx2[k][j][i]);  ,

             dvx3 =   Ax3[k][j][i]*  (vx3[k+1][j][i] + vx3[k][j][i])
                    - Ax3[k-1][j][i]*(vx3[k-1][j][i] + vx3[k][j][i]))

    divv = (D_EXPAND(dvx1, + dvx2, + dvx3))/grid->dV[k][j][i];
    #endif
 
    if (divv < 0.0){
    
    /* -----------------------------------------------
        Compute undivided difference of the total
        pressure and minimum value in neighbour zones
       ----------------------------------------------- */
       
      pt_min = pt[k][j][i];
      D_EXPAND(pt_min1 = MIN(pt[k][j][i+1], pt[k][j][i-1]); ,
               pt_min2 = MIN(pt[k][j+1][i], pt[k][j-1][i]);  ,
               pt_min3 = MIN(pt[k+1][j][i], pt[k-1][j][i]); )

      D_EXPAND(pt_min = MIN(pt_min, pt_min1);  ,
               pt_min = MIN(pt_min, pt_min2);  ,
               pt_min = MIN(pt_min, pt_min3);)
      
      D_EXPAND(dpx1 = fabs(pt[k][j][i+1] - pt[k][j][i-1]);  ,
               dpx2 = fabs(pt[k][j+1][i] - pt[k][j-1][i]);  ,
               dpx3 = fabs(pt[k+1][j][i] - pt[k-1][j][i]);)
                
      gradp = D_EXPAND(dpx1, + dpx2, + dpx3);
      if (gradp > PARTICLES_LP_SHK_THRESHOLD*pt_min) {
        for(i1 = -nbuf; i1<=nbuf;i1++){    
        for(j1 = -nbuf; j1<=nbuf;j1++){
        for(k1 = -nbuf; k1<=nbuf;k1++){
          itag = MAX(i+i1, 0.0);
          itag = MIN(itag, NX1_TOT-1);
          jtag = MAX(j+j1, 0.0);
          jtag = MIN(jtag, NX2_TOT-1);
          ktag = MAX(k+k1, 0.0);
          ktag = MIN(ktag, NX3_TOT-1);
          flag2[ktag][jtag][itag]  = 1.0;
        }}}
      }
    }  /* end if (divv < 0) */
  }}}
}

/* ********************************************************************* */
void Particles_LP_GradP(double *gradval, int VAR, double ***uu[], 
                     Grid *grid, int indici[]) 
/*
 * This function computes the normal gradient of the variable  
 * VAR for each particles based on its position given by indici.
 *
 *********************************************************************** */
{
  int il, ic, ir, jl, jc, jr, kl, kc, kr;
  double dx1, dx2, dy1, dy2, dz1, dz2;
  double ax, bx, cx, ay, by, cy, az, bz, cz;
  double DelPx, DelPy, DelPz; 

#if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)
#if DIMENSIONS == 1

  ic = indici[0];
  il = ic - 1;
  ir = ic + 1;

  dx1 = grid->x[IDIR][ic] - grid->x[IDIR][il];
  dx2 = grid->x[IDIR][ir] - grid->x[IDIR][ic];

  ax = (dx2 * dx2 - dx1 * dx1) / (dx1 * dx2 * (dx1 + dx2));
  bx = -dx2 / (dx1 * (dx1 + dx2));
  cx = dx1 / (dx2 * (dx1 + dx2));
  DelPx = (ax * uu[VAR][0][0][ic] + bx * uu[VAR][0][0][il] + cx * uu[VAR][0][0][ir]);

  gradval[0] = DelPx 

#elif DIMENSIONS == 2

  ic = indici[0];
  il = ic - 1;
  ir = ic + 1;

  jc = indici[1];
  jl = jc - 1;
  jr = jc + 1;

  dx1 = grid->x[IDIR][ic] - grid->x[IDIR][il];
  dx2 = grid->x[IDIR][ir] - grid->x[IDIR][ic];

  dy1 = grid->x[JDIR][jc] - grid->x[JDIR][jl];
  dy2 = grid->x[JDIR][jr] - grid->x[JDIR][jc];


  ax = (dx2 * dx2 - dx1 * dx1) / (dx1 * dx2 * (dx1 + dx2));
  bx = -dx2 / (dx1 * (dx1 + dx2));
  cx = dx1 / (dx2 * (dx1 + dx2));

  ay = (dy2 * dy2 - dy1 * dy1) / (dy1 * dy2 * (dy1 + dy2));
  by = -dy2 / (dy1 * (dy1 + dy2));
  cy = dy1 / (dy2 * (dy1 + dy2));
  
  DelPx = ax * uu[VAR][0][jc][ic] + bx * uu[VAR][0][jc][il] + cx * uu[VAR][0][jc][ir];
  DelPy = ay * uu[VAR][0][jc][ic] + by * uu[VAR][0][jl][ic] + cy * uu[VAR][0][jr][ic];

  gradval[0] = DelPx;
  gradval[1] = DelPy;

#elif DIMENSIONS == 3

  ic = indici[0];
  il = ic - 1;
  ir = ic + 1;

  jc = indici[1];
  jl = jc - 1;
  jr = jc + 1;

  kc = indici[2];
  kl = kc - 1;
  kr = kc + 1;

  dx1 = grid->x[IDIR][ic] - grid->x[IDIR][il];
  dx2 = grid->x[IDIR][ir] - grid->x[IDIR][ic];

  dy1 = grid->x[JDIR][jc] - grid->x[JDIR][jl];
  dy2 = grid->x[JDIR][jr] - grid->x[JDIR][jc];

  dz1 = grid->x[KDIR][kc] - grid->x[KDIR][kl];
  dz2 = grid->x[KDIR][kr] - grid->x[KDIR][kc];

  ax = (dx2 * dx2 - dx1 * dx1) / (dx1 * dx2 * (dx1 + dx2));
  bx = -dx2 / (dx1 * (dx1 + dx2));
  cx = dx1 / (dx2 * (dx1 + dx2));

  ay = (dy2 * dy2 - dy1 * dy1) / (dy1 * dy2 * (dy1 + dy2));
  by = -dy2 / (dy1 * (dy1 + dy2));
  cy = dy1 / (dy2 * (dy1 + dy2));

  az = (dz2 * dz2 - dz1 * dz1) / (dz1 * dz2 * (dz1 + dz2));
  bz = -dz2 / (dz1 * (dz1 + dz2));
  cz = dz1 / (dz2 * (dz1 + dz2));

  DelPx = ax * uu[VAR][kc][jc][ic] + bx * uu[VAR][kc][jc][il] + cx * uu[VAR][kc][jc][ir];
  DelPy = ay * uu[VAR][kc][jc][ic] + by * uu[VAR][kc][jl][ic]+ cy * uu[VAR][kc][jr][ic];
  DelPz = az * uu[VAR][kc][jc][ic] + bz * uu[VAR][kl][jc][ic] + cz * uu[VAR][kr][jc][ic];

  gradval[0] = DelPx;
  gradval[1] = DelPy;
  gradval[2] = DelPz;

#endif  /* DIMENSIONS == 3*/
#endif /*  (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL) */

}
#endif  /* PARTICLES_LP_SPECTRA == YES */
#endif /* PARTICLES_TYPE == LAGRANGIAN */
