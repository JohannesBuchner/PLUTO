/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Performs various magnetic field averaging operations.

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 04, 2017
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
 * \param [in]   bf    array of staggered fields
 * \param [out]  UU    array of conservative variables
 * \param [in]   grid  pointer to Grid structure
 ************************************************************************ */
{
  int i, j, k;
  double b2_old, b2_new, bx_ave, by_ave, bz_ave;
  double rp, rm;
  double ***bx, ***by, ***bz;
  double *dx, *dy;
  double *r;
  
  D_EXPAND(bx = bf[BX1s];  ,
           by = bf[BX2s];  ,
           bz = bf[BX3s]; )
  
  dx  = grid->dx[IDIR]; 
  dy  = grid->dx[JDIR];
  r   = grid->x[IDIR];

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

     D_EXPAND( bx_ave = 0.5*(bx[k][j][i] + bx[k][j][i-1]);  ,
               by_ave = 0.5*(by[k][j][i] + by[k][j-1][i]);  ,
               bz_ave = 0.5*(bz[k][j][i] + bz[k-1][j][i]); )

    #elif GEOMETRY == CYLINDRICAL
/*
     bx_ave =     (A1[i]*bx[k][j][i] + A1[i-1]*bx[k][j][i - 1])/(A1[i] + A1[i-1]);
     by_ave = 0.5*(   by[k][j][i] +    by[k][j - 1][i]);
*/

     bx_ave = 0.5*(bx[k][j][i] + bx[k][j][i-1]);
     by_ave = 0.5*(by[k][j][i] + by[k][j-1][i]);

    #elif GEOMETRY == POLAR

     rp = grid->xr[IDIR][i];
     rm = grid->xl[IDIR][i];

     D_EXPAND(bx_ave = (rp*bx[k][j][i] + rm*bx[k][j][i-1])/(rp + rm); ,
              by_ave = 0.5*(by[k][j][i] + by[k][j-1][i]);             ,
              bz_ave = 0.5*(bz[k][j][i] + bz[k-1][j][i]);)

    #elif GEOMETRY == SPHERICAL
    {
     double Ap, Am, dV1;
     double thp, thm, sp, sm, dV2;
     rp  = grid->xr[IDIR][i];
     rm  = grid->xl[IDIR][i];
     Ap  = rp*rp;
     Am  = rm*rm;
     dV1 = (rp*rp*rp - rm*rm*rm)/3.0;
     
     thp = grid->xr[JDIR][j];
     thm = grid->xl[JDIR][j];
     sp  = fabs(sin(thp));
     sm  = fabs(sin(thm));
     dV2 = fabs(cos(thm) - cos(thp));

     D_EXPAND(bx_ave  = 0.5*(Ap*bx[k][j][i] + Am*bx[k][j][i - 1]);
              bx_ave *= dx[i]/dV1;                                            ,

              by_ave  = 0.5*(sp*by[k][j][i] + sm*by[k][j - 1][i]);
              by_ave *= dy[j]/dV2;                                        ,

              bz_ave  = 0.5*(bz[k][j][i] + bz[k - 1][j][i]);)
    }
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
  double ***Bx, ***By, ***Bz;
  double ***bx, ***by, ***bz;
  double *rp = grid->xr[IDIR];
  double *rm = grid->xl[IDIR];
  double Ap, Am, thp, thm;

/* ------------------------------------------------------
                   X1 boundaries
   ------------------------------------------------------ */

  if (side == X1_BEG){
    Bx = d->Vc[BX1]; bx = d->Vs[BX1s];
    KTOT_LOOP(k) JTOT_LOOP(j) for (i = 0; i < IBEG; i++) {
      #if GEOMETRY == CARTESIAN 
      Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i-1]);  
      #elif GEOMETRY == CYLINDRICAL
      Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i-1]);  
      #elif GEOMETRY == POLAR
      Bx[k][j][i] = (rp[i]*bx[k][j][i] + rm[i]*bx[k][j][i-1])/(rp[i] + rm[i]); 
      #elif GEOMETRY == SPHERICAL
      Ap = rp[i]*rp[i];
      Am = rm[i]*rm[i];

      Bx[k][j][i] = (Ap*bx[k][j][i] + Am*bx[k][j][i-1])/(Ap + Am);   
      #endif
    }    
  }else if (side == X1_END){
    Bx = d->Vc[BX1]; bx = d->Vs[BX1s];
    KTOT_LOOP(k) JTOT_LOOP(j) for (i = IEND+1; i < NX1_TOT; i++) {
      #if GEOMETRY == CARTESIAN 
      Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i-1]);  
      #elif GEOMETRY == CYLINDRICAL
      Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i-1]);
      #elif GEOMETRY == POLAR
      Bx[k][j][i] = (rp[i]*bx[k][j][i] + rm[i]*bx[k][j][i-1])/(rp[i] + rm[i]); 
      #elif GEOMETRY == SPHERICAL
      Ap = rp[i]*rp[i];
      Am = rm[i]*rm[i];
      Bx[k][j][i] = (Ap*bx[k][j][i] + Am*bx[k][j][i-1])/(Ap + Am);   
      #endif
    }
  }

/* ------------------------------------------------------
                   X2 boundaries
   ------------------------------------------------------ */

  if (side == X2_BEG){
    By = d->Vc[BX2]; by = d->Vs[BX2s];
    KTOT_LOOP(k) for (j = 0; j < JBEG; j++) ITOT_LOOP(i){
      #if GEOMETRY == CARTESIAN 
      By[k][j][i] = 0.5*(by[k][j][i] + by[k][j-1][i]);  
      #elif GEOMETRY == CYLINDRICAL
      By[k][j][i] = 0.5*(by[k][j][i] + by[k][j-1][i]);  
      #elif GEOMETRY == POLAR
      By[k][j][i] = 0.5*(by[k][j][i] + by[k][j-1][i]);                
      #elif GEOMETRY == SPHERICAL
      thp = grid->xr[JDIR][j];
      thm = grid->xl[JDIR][j];
      Ap  = fabs(sin(thp));
      Am  = fabs(sin(thm));

      By[k][j][i] = (Ap*by[k][j][i] + Am*by[k][j-1][i])/(Ap + Am); 
      #endif
    }    
  }else if (side == X2_END){
    By = d->Vc[BX2]; by = d->Vs[BX2s];
    KTOT_LOOP(k) for (j = JEND+1; j < NX2_TOT; j++) ITOT_LOOP(i){
      #if GEOMETRY == CARTESIAN 
      By[k][j][i] = 0.5*(by[k][j][i] + by[k][j-1][i]);  
      #elif GEOMETRY == CYLINDRICAL
      By[k][j][i] = 0.5*(by[k][j][i] + by[k][j-1][i]);  
      #elif GEOMETRY == POLAR
      By[k][j][i] = 0.5*(by[k][j][i] + by[k][j-1][i]);   
      #elif GEOMETRY == SPHERICAL
      thp = grid->xr[JDIR][j];
      thm = grid->xl[JDIR][j];
      Ap  = fabs(sin(thp));
      Am  = fabs(sin(thm));
      By[k][j][i] = (Ap*by[k][j][i] + Am*by[k][j-1][i])/(Ap + Am); 
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
      Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]); 
      #elif GEOMETRY == CYLINDRICAL
      Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]); 
      #elif GEOMETRY == POLAR
      Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]);
      #elif GEOMETRY == SPHERICAL
      Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]);
      #endif
    }    
  }else if (side == X3_END){
    Bz = d->Vc[BX3]; bz = d->Vs[BX3s];
    for (k = KEND+1; k < NX3_TOT; k++) JTOT_LOOP(j) ITOT_LOOP(i){
      #if GEOMETRY == CARTESIAN 
      Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]); 
      #elif GEOMETRY == CYLINDRICAL
      Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]); 
      #elif GEOMETRY == POLAR
      Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]);
      #elif GEOMETRY == SPHERICAL
      Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]);
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
  double ***Bx, ***By, ***Bz;
  double ***bx, ***by, ***bz;
  double *rp = grid->xr[IDIR];
  double *rm = grid->xl[IDIR];
  double thp, thm, A1p, A1m, A2p, A2m;

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
               By[k][j][i] = 0.5*(by[k][j][i] + by[k][j-1][i]);  ,
               Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]); )
      #else
      thp = grid->xr[JDIR][j];
      thm = grid->xl[JDIR][j];
      A2p = fabs(sin(thp));
      A2m = fabs(sin(thm));
      D_EXPAND(                                                      ,
               By[k][j][i] = (A2p*by[k][j][i] + A2m*by[k][j - 1][i])/
                             (A2p + A2m);                              ,
               Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]);)
      #endif
    }    
  }else if (side == X1_END){
    X1_END_LOOP(k,j,i){
      #if GEOMETRY != SPHERICAL
      D_EXPAND(                                                    ,
               By[k][j][i] = 0.5*(by[k][j][i] + by[k][j-1][i]);  ,
               Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]); )
      #else
      thp = grid->xr[JDIR][j];
      thm = grid->xl[JDIR][j];
      A2p = fabs(sin(thp));
      A2m = fabs(sin(thm));
      D_EXPAND(                                                        ,
               By[k][j][i] = (A2p*by[k][j][i] + A2m*by[k][j - 1][i])/
                             (A2p + A2m);                              ,
               Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]);)
      #endif
    }    
  }

/* ------------------------------------------------------
                   X2 boundaries
   ------------------------------------------------------ */

  if (side == X2_BEG){
    X2_BEG_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN 
      D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i-1]);  ,
                                                                    ,
                Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]); )
      #elif GEOMETRY == CYLINDRICAL
      D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i-1]);  ,
                                                                    ,
                Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]); )
      #elif GEOMETRY == POLAR || GEOMETRY == SPHERICAL
      D_EXPAND(Bx[k][j][i] = (rp[i]*bx[k][j][i] + rm[i]*bx[k][j][i - 1])/
                             (rp[i] + rm[i]);                                ,
                                                                             ,
               Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]);)
      #endif
    }
  }else if (side == X2_END){
    X2_END_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN 
      D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i-1]);  ,
                                                                    ,
                Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]); )
      #elif GEOMETRY == CYLINDRICAL
      D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i-1]);  ,
                                                                    ,
                Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]); )
      #elif GEOMETRY == POLAR || GEOMETRY == SPHERICAL
      D_EXPAND(Bx[k][j][i] = (rp[i]*bx[k][j][i] + rm[i]*bx[k][j][i-1])/
                             (rp[i] + rm[i]);                            ,
                                                                          ,
               Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k-1][j][i]);)
      #endif
    }    
  }

/* ------------------------------------------------------
                   X3 boundaries
   ------------------------------------------------------ */

  if (side == X3_BEG){
    X3_BEG_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN 
      D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i-1]);  ,
                By[k][j][i] = 0.5*(by[k][j][i] + by[k][j-1][i]);  ,
                                                                   )

      #elif GEOMETRY == CYLINDRICAL
      D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i-1]);  ,
                By[k][j][i] = 0.5*(by[k][j][i] + by[k][j-1][i]);  ,
                                                                  )
      #elif GEOMETRY == POLAR
      D_EXPAND(Bx[k][j][i] = (rp[i]*bx[k][j][i] + rm[i]*bx[k][j][i-1])/
                             (rp[i] + rm[i]);                             ,
               By[k][j][i] = 0.5*(by[k][j][i] + by[k][j-1][i]);           ,
                                                                )
      #elif GEOMETRY == SPHERICAL
      A1p = rp[i]*rp[i];
      A1m = rm[i]*rm[i];
      thp = grid->xr[JDIR][j];
      thm = grid->xl[JDIR][j];
      A2p = fabs(sin(thp));
      A2m = fabs(sin(thm));
      D_EXPAND(Bx[k][j][i] = (A1p*bx[k][j][i] + A1m*bx[k][j][i-1])/ 
                             (A1p + A1m);                              ,
               By[k][j][i] = (A2p*by[k][j][i] + A2m*by[k][j-1][i])/
                             (A2p + A2m);                              ,
                                                             )
      #endif
    }
  }else if (side == X3_END){
    X3_END_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN 
      D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i-1]);  ,
                By[k][j][i] = 0.5*(by[k][j][i] + by[k][j-1][i]);  ,
                                                                    )
      #elif GEOMETRY == CYLINDRICAL
      D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i-1]);  ,
                By[k][j][i] = 0.5*(by[k][j][i] + by[k][j-1][i]);  ,
                                                                  )
      #elif GEOMETRY == POLAR
      D_EXPAND(Bx[k][j][i] = (rp[i]*bx[k][j][i] + rm[i]*bx[k][j][i-1])/
                             (rp[i] + rm[i]);                             ,
               By[k][j][i] = 0.5*(by[k][j][i] + by[k][j-1][i]);           ,
                                                               )
      #elif GEOMETRY == SPHERICAL
      A1p = rp[i]*rp[i];
      A1m = rm[i]*rm[i];
      thp = grid->xr[JDIR][j];
      thm = grid->xl[JDIR][j];
      A2p = fabs(sin(thp));
      A2m = fabs(sin(thm));
      D_EXPAND(Bx[k][j][i] = (A1p*bx[k][j][i] + A1m*bx[k][j][i-1])/ 
                             (A1p + A1m);                              ,
               By[k][j][i] = (A2p*by[k][j][i] + A2m*by[k][j-1][i])/
                             (A2p + A2m);                              ,
                                                            )
      #endif
    }    
  }
}
#undef IEXT_LOOP
#undef JEXT_LOOP
#undef KEXT_LOOP
