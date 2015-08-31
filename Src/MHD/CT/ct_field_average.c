/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Performs various magnetic field averaging operations.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void CT_AverageMagneticField (double ****bf, double ****UU, Grid *grid)
/*!
 * Average staggered magnetic field components to form a zone-centered 
 * field, e.g., \f$ <B_i> = (b_{i+1/2} + b_{i-1/2})/2\f$
 * 
 * The averaging procedure is performed inside the computational 
 * domain augmented by an additional row/plane of ghost zones in the 
 * boundary regions.
 * Inclusion of boundary regions is useful only during the predictor 
 * step of the CT-CTU algorithm, whereas is useless for RK time stepping.
 *
 * When the CT_EN_CORRECTION flag is enabled, we also redefine the
 * zone total energy using the newly formed cell-centered field, i.e.,
 *
 *  \f[ 
 *   E_i \to E_i - \frac{\mathbf{B}_i^2}{2} + \frac{<\mathbf{B}_i^2>}{2}
 *  \f]
 *
 * \param [in]   bf   array of staggered fields
 * \param [out]  UU   array of conservative variables
 * \param [in]   grid  pointer to Grid structure
 ************************************************************************ */
{
  int i, j, k;
  double b2_old, b2_new, bx_ave, by_ave, bz_ave;
  double rp, rm;
  double ***bx, ***by, ***bz;
  double *dx, *dy, *A1, *A2, *dV1, *dV2;
  double *r;
  
  D_EXPAND(bx = bf[BX1s];  ,
           by = bf[BX2s];  ,
           bz = bf[BX3s]; )
  
  dx = grid[IDIR].dx; 
  dy = grid[JDIR].dx;
  A1  = grid[IDIR].A;
  A2  = grid[JDIR].A;
  dV1 = grid[IDIR].dV;
  dV2 = grid[JDIR].dV;
  r   = grid[IDIR].x;

/* ---------------------------------------------------------
    Loop over all zones in the domain. 
    We include also one set of boundary zones which 
    is useful only during the predictor step of the CT-CTU 
    algorithm, whereas is useless for RK time stepping.
   --------------------------------------------------------- */

  for (k = KBEG - KOFFSET; k <= KEND + KOFFSET; k++){
  for (j = JBEG - JOFFSET; j <= JEND + JOFFSET; j++){
  for (i = IBEG - IOFFSET; i <= IEND + IOFFSET; i++){

    #if GEOMETRY == CARTESIAN 

     D_EXPAND( bx_ave = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  ,
               by_ave = 0.5*(by[k][j][i] + by[k][j - 1][i]);  ,
               bz_ave = 0.5*(bz[k][j][i] + bz[k - 1][j][i]); )

    #elif GEOMETRY == CYLINDRICAL
/*
     bx_ave =     (A1[i]*bx[k][j][i] + A1[i-1]*bx[k][j][i - 1])/(A1[i] + A1[i-1]);
     by_ave = 0.5*(   by[k][j][i] +    by[k][j - 1][i]);
*/

     bx_ave = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);
     by_ave = 0.5*(by[k][j][i] + by[k][j - 1][i]);

    #elif GEOMETRY == POLAR

     rp = grid[IDIR].xr[i];
     rm = grid[IDIR].xl[i];

     D_EXPAND(bx_ave = (rp*bx[k][j][i] + rm*bx[k][j][i - 1])/(rp + rm); ,
              by_ave = 0.5*(by[k][j][i] + by[k][j - 1][i]);             ,
              bz_ave = 0.5*(bz[k][j][i] + bz[k - 1][j][i]);)

    #elif GEOMETRY == SPHERICAL

     D_EXPAND(bx_ave  = 0.5*(A1[i]*bx[k][j][i] + A1[i - 1]*bx[k][j][i - 1]);
              bx_ave *= dx[i]/dV1[i];                                        ,

              by_ave  = 0.5*(A2[j]*by[k][j][i] + A2[j - 1]*by[k][j - 1][i]);
              by_ave *= dy[j]/dV2[j];                                        ,

              bz_ave  = 0.5*(bz[k][j][i] + bz[k - 1][j][i]);)
    #endif
   
   /* ---------------------------------------------
          apply energy correction if necessary 
       --------------------------------------------- */
 
    #if CT_EN_CORRECTION == YES && HAVE_ENERGY
     b2_old = D_EXPAND(  UU[k][j][i][BX1]*UU[k][j][i][BX1], 
                       + UU[k][j][i][BX2]*UU[k][j][i][BX2],  
                       + UU[k][j][i][BX3]*UU[k][j][i][BX3]);
    #endif   

    D_EXPAND( UU[k][j][i][BX1] = bx_ave;  ,
              UU[k][j][i][BX2] = by_ave;  ,
              UU[k][j][i][BX3] = bz_ave; )

    #if CT_EN_CORRECTION == YES && HAVE_ENERGY
     b2_new = D_EXPAND(bx_ave*bx_ave, + by_ave*by_ave, + bz_ave*bz_ave);
     UU[k][j][i][ENG] += 0.5*(b2_new - b2_old);
    #endif
  }}}

}
/* ********************************************************************* */
void CT_AverageNormalMagField (const Data *d, int side, Grid *grid)
/*!
 * \brief Compute the normal component of the volume-average magnetic 
 *        field from the staggered components in the ghost zones.
 * \details
 *
 *   For a given "side" of the boundary, average the staggered magnetic  
 *   field components normal to that boundary into a cell-centered field. 
 *   For instance, at X3_BEG we only average Bz.
 *   Transverse components should be averaged using a slightly different
 *   stencil.
 *   This function is automatically called from Boundary() for 
 *   consistency between staggered and zone-centered fields.
 * 
 * \param [in,out]    d  pointer to PLUTO Data structure
 * \param [in]     side  the side 
 * \param [in]     grid  pointer to PLUTO Grid structure
 *
 * \todo   replace the loops with more compact macro, such as
 *         X1_BEG_LOOP()...
 *********************************************************************** */
{
  int    i, j, k;
  double *Ap, *Am;
  double ***Bx, ***By, ***Bz;
  double ***bx, ***by, ***bz;

/* ------------------------------------------------------
                   X1 boundaries
   ------------------------------------------------------ */

  if (side == X1_BEG){
    Ap = grid[IDIR].A;  Am = Ap - 1;
    Bx = d->Vc[BX1]; bx = d->Vs[BX1s];
    KTOT_LOOP(k) JTOT_LOOP(j) for (i = 0; i < IBEG; i++) {
      #if GEOMETRY == CARTESIAN 
       Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  
      #elif GEOMETRY == CYLINDRICAL
       Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  
      #elif GEOMETRY == POLAR
       Bx[k][j][i] = (Ap[i]*bx[k][j][i] + Am[i]*bx[k][j][i - 1])/(Ap[i] + Am[i]); 
      #elif GEOMETRY == SPHERICAL
       Bx[k][j][i] = (Ap[i]*bx[k][j][i] + Am[i]*bx[k][j][i - 1])/(Ap[i] + Am[i]);   
      #endif
    }    
  }else if (side == X1_END){
    Ap = grid[IDIR].A; Am = Ap - 1;
    Bx = d->Vc[BX1]; bx = d->Vs[BX1s];
    KTOT_LOOP(k) JTOT_LOOP(j) for (i = IEND+1; i < NX1_TOT; i++) {
      #if GEOMETRY == CARTESIAN 
       Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  
      #elif GEOMETRY == CYLINDRICAL
       Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  
      #elif GEOMETRY == POLAR
       Bx[k][j][i] = (Ap[i]*bx[k][j][i] + Am[i]*bx[k][j][i - 1])/(Ap[i] + Am[i]); 
      #elif GEOMETRY == SPHERICAL
       Bx[k][j][i] = (Ap[i]*bx[k][j][i] + Am[i]*bx[k][j][i - 1])/(Ap[i] + Am[i]);   
      #endif
    }
  }

/* ------------------------------------------------------
                   X2 boundaries
   ------------------------------------------------------ */

  if (side == X2_BEG){
    Ap = grid[JDIR].A; Am = Ap - 1;
    By = d->Vc[BX2]; by = d->Vs[BX2s];
    KTOT_LOOP(k) for (j = 0; j < JBEG; j++) ITOT_LOOP(i){
      #if GEOMETRY == CARTESIAN 
       By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);  
      #elif GEOMETRY == CYLINDRICAL
       By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);  
      #elif GEOMETRY == POLAR
       By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);                
      #elif GEOMETRY == SPHERICAL
       By[k][j][i] = (Ap[j]*by[k][j][i] + Am[j]*by[k][j - 1][i])/(Ap[j] + Am[j]); 
      #endif
    }    
  }else if (side == X2_END){
    Ap = grid[JDIR].A; Am = Ap - 1;
    By = d->Vc[BX2]; by = d->Vs[BX2s];
    KTOT_LOOP(k) for (j = JEND+1; j < NX2_TOT; j++) ITOT_LOOP(i){
      #if GEOMETRY == CARTESIAN 
       By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);  
      #elif GEOMETRY == CYLINDRICAL
       By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);  
      #elif GEOMETRY == POLAR
       By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);   
      #elif GEOMETRY == SPHERICAL
       By[k][j][i] = (Ap[j]*by[k][j][i] + Am[j]*by[k][j - 1][i])/(Ap[j] + Am[j]); 
      #endif
    }    
  }

/* ------------------------------------------------------
                   X3 boundaries
   ------------------------------------------------------ */

  if (side == X3_BEG){
    Bz = d->Vc[BX3]; bz = d->Vs[BX3s];
    for (k = 0; k < KBEG; k++) JTOT_LOOP(j) ITOT_LOOP(i){
      #if GEOMETRY == CARTESIAN 
       Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]); 
      #elif GEOMETRY == CYLINDRICAL
       Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]); 
      #elif GEOMETRY == POLAR
       Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]);
      #elif GEOMETRY == SPHERICAL
       Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]);
      #endif
    }    
  }else if (side == X3_END){
    Bz = d->Vc[BX3]; bz = d->Vs[BX3s];
    for (k = KEND+1; k < NX3_TOT; k++) JTOT_LOOP(j) ITOT_LOOP(i){
      #if GEOMETRY == CARTESIAN 
       Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]); 
      #elif GEOMETRY == CYLINDRICAL
       Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]); 
      #elif GEOMETRY == POLAR
       Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]);
      #elif GEOMETRY == SPHERICAL
       Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]);
      #endif
    }    
  }
}

/* ********************************************************************* */
void CT_AverageTransverseMagField (const Data *d, int side, Grid *grid)
/*!
 * \brief Compute the transverse component of the volume-average magnetic 
 *        field from the staggered components in the ghost zones.
 * \details
 *
 *  For a given "side" of the boundary, average the staggered magnetic
 *  field components transverse to that boundary into a cell-centered  
 *  field.
 *  For instance, at X2_BEG we average Bx and Bz.
 *  This call is necessary only if one requires consistency between
 *  staggered and cell-centered magnetic field components.
 * 
 * \param [in,out]    d  pointer to PLUTO Data structure
 * \param [in]     side  the side 
 * \param [in]     grid  pointer to PLUTO Grid structure
 *
 * \todo   replace the loops using more compact notations, such as
 *         X1_BEG_LOOP()...
 *********************************************************************** */
#define IEXT_LOOP(i) for (i = IOFFSET; i < NX1_TOT - IOFFSET; i++)
#define JEXT_LOOP(j) for (j = JOFFSET; j < NX2_TOT - JOFFSET; j++)
#define KEXT_LOOP(k) for (k = KOFFSET; k < NX3_TOT - KOFFSET; k++)
{
  int    i, j, k;
  double *A1, *A2, *A3;
  double ***Bx, ***By, ***Bz;
  double ***bx, ***by, ***bz;

  A1 = grid[IDIR].A;
  A2 = grid[JDIR].A;
  A3 = grid[KDIR].A;

  D_EXPAND(Bx = d->Vc[BX1]; bx = d->Vs[BX1s];  ,
           By = d->Vc[BX2]; by = d->Vs[BX2s];  ,
           Bz = d->Vc[BX3]; bz = d->Vs[BX3s]; )

/* ------------------------------------------------------
                   X1 boundaries
   ------------------------------------------------------ */

  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i){
      #if GEOMETRY != SPHERICAL
       D_EXPAND(                                                    ,
                By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);  ,
                Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]); )
      #else
       D_EXPAND(                                                              ,
                By[k][j][i] = (A2[j]*by[k][j][i] + A2[j-1]*by[k][j - 1][i])/
                              (A2[j] + A2[j-1]);                              ,
                Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]);)
      #endif
    }    
  }else if (side == X1_END){
    X1_END_LOOP(k,j,i){
      #if GEOMETRY != SPHERICAL
       D_EXPAND(                                                    ,
                By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);  ,
                Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]); )
      #else
       D_EXPAND(                                                              ,
                By[k][j][i] = (A2[j]*by[k][j][i] + A2[j-1]*by[k][j - 1][i])/
                              (A2[j] + A2[j-1]);                              ,
                Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]);)
      #endif
    }    
  }

/* ------------------------------------------------------
                   X2 boundaries
   ------------------------------------------------------ */

  if (side == X2_BEG){
    X2_BEG_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN 
       D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  ,
                                                                     ,
                 Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]); )
      #elif GEOMETRY == CYLINDRICAL
       D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  ,
                                                                     ,
                 Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]); )
      #elif GEOMETRY == POLAR || GEOMETRY == SPHERICAL
       D_EXPAND(Bx[k][j][i] = (A1[i]*bx[k][j][i] + A1[i-1]*bx[k][j][i - 1])/
                              (A1[i] + A1[i-1]);                              ,
                                                                              ,
                Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]);)
      #endif
    }
  }else if (side == X2_END){
    X2_END_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN 
       D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  ,
                                                                     ,
                 Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]); )
      #elif GEOMETRY == CYLINDRICAL
       D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  ,
                                                                     ,
                 Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]); )
      #elif GEOMETRY == POLAR || GEOMETRY == SPHERICAL
       D_EXPAND(Bx[k][j][i] = (A1[i]*bx[k][j][i] + A1[i-1]*bx[k][j][i - 1])/
                              (A1[i] + A1[i-1]);                              ,
                                                                              ,
                Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]);)
      #endif
    }    
  }

/* ------------------------------------------------------
                   X3 boundaries
   ------------------------------------------------------ */

  if (side == X3_BEG){
    X3_BEG_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN 
       D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  ,
                 By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);  ,
                                                                    )

      #elif GEOMETRY == CYLINDRICAL
       D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  ,
                 By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);  ,
                                                                  )
      #elif GEOMETRY == POLAR
       D_EXPAND(Bx[k][j][i] = (A1[i]*bx[k][j][i] + A1[i-1]*bx[k][j][i - 1])/
                              (A1[i] + A1[i-1]);                             ,
                By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);           ,
                                                                 )
      #elif GEOMETRY == SPHERICAL

       D_EXPAND(Bx[k][j][i] = (A1[i]*bx[k][j][i] + A1[i-1]*bx[k][j][i - 1])/ 
                              (A1[i] + A1[i-1]);                              ,
                By[k][j][i] = (A2[j]*by[k][j][i] + A2[j-1]*by[k][j - 1][i])/
                              (A2[j] + A2[j-1]);                              ,
                                                             )
      #endif
    }
  }else if (side == X3_END){
    X3_END_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN 
       D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  ,
                 By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);  ,
                                                                    )

      #elif GEOMETRY == CYLINDRICAL
       D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  ,
                 By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);  ,
                                                                  )
      #elif GEOMETRY == POLAR
       D_EXPAND(Bx[k][j][i] = (A1[i]*bx[k][j][i] + A1[i-1]*bx[k][j][i - 1])/
                              (A1[i] + A1[i-1]);                             ,
                By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);           ,
                                                                 )
      #elif GEOMETRY == SPHERICAL

       D_EXPAND(Bx[k][j][i] = (A1[i]*bx[k][j][i] + A1[i-1]*bx[k][j][i - 1])/ 
                              (A1[i] + A1[i-1]);                              ,
                By[k][j][i] = (A2[j]*by[k][j][i] + A2[j-1]*by[k][j - 1][i])/
                              (A2[j] + A2[j-1]);                              ,
                                                             )
      #endif
    }    
  }
}
#undef IEXT_LOOP
#undef JEXT_LOOP
#undef KEXT_LOOP
