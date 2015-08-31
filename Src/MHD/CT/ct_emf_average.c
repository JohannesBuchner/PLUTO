/* ///////////////////////////////////////////////////////////////////// */
/*! 
 * \file  
 * \brief  Collects different EMF averaging schemes.                     */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define EX(k,j,i)  (vz[k][j][i]*By[k][j][i] - vy[k][j][i]*Bz[k][j][i])
#define EY(k,j,i)  (vx[k][j][i]*Bz[k][j][i] - vz[k][j][i]*Bx[k][j][i])
#define EZ(k,j,i)  (vy[k][j][i]*Bx[k][j][i] - vx[k][j][i]*By[k][j][i])

/* ********************************************************************* */
void CT_EMF_ArithmeticAverage (const EMF *Z1, const double w)
/*!
 * \brief Compute arithmetic average of EMF at cell edges.
 * \details
 *
 *  Combine the four electric field values computed at zone faces 
 *  as upwind Godunov fluxes into an edge-centered value.\n
 *  The face-centered EMF should have been stored by previous calls
 *  to CT_StoreEMF() during the one-dimensional sweeps.\n
 *  This function employs a simple arithmetic averaging of the 
 *  face-centered electric field.
 *
 * \b References:
 *    - "A Staggered Mesh Algorithm Using High Order Godunov Fluxes to 
 *       Ensure Solenoidal Magnetic Fields in Magnetohydrodynamic 
 *       Simulations"\n
 *       Balsara \& Spicer, JCP (1999) 149, 270.
 *
 * \param [in]      Z1    pointer to EMF structure
 * \param [in]      w     weighting factor
 *
 * \return  This function has no return value.
 *
 *
 * \author A. Mignone (mignone@ph.unito.it)
 * \date   Aug 16, 31, 2012
 *********************************************************************** */
{
  int i, j, k;

  for (k = Z1->kbeg; k <= Z1->kend; k++){
  for (j = Z1->jbeg; j <= Z1->jend; j++){
  for (i = Z1->ibeg; i <= Z1->iend; i++){      
    #if DIMENSIONS == 3
     Z1->ex[k][j][i] = w*(  Z1->exk[k][j][i] + Z1->exk[k][j + 1][i] 
                          + Z1->exj[k][j][i] + Z1->exj[k + 1][j][i]);
     Z1->ey[k][j][i] = w*(  Z1->eyi[k][j][i] + Z1->eyi[k + 1][j][i] 
                          + Z1->eyk[k][j][i] + Z1->eyk[k][j][i + 1]);
    #endif 
    Z1->ez[k][j][i] = w*(  Z1->ezi[k][j][i] + Z1->ezi[k][j + 1][i] 
                         + Z1->ezj[k][j][i] + Z1->ezj[k][j][i + 1]);
  }}}
}

/* ********************************************************************* */
void CT_EMF_IntegrateToCorner (const Data *d, const EMF *emf, Grid *grid)
/*!
 * \details
 *
 *  Add derivatives to the 4-point arithmetic average of magnetic fields.
 *  Obtain the electric field at corners. 
 *
 * \b References: 
 *  - "An unsplit Godunov method for ideal MHD via constrained transport"\n
 *    Gardiner & Stone, JCP (2005) 205, 509.
 *    See Eq. (41), (45) and (50).
 *
 * \param [in]      d     pointer to PLUTO Data structure
 * \param [in]      emf   pointer to EMF structure
 * \param [in]      grid  pointer to Grid structure
 *
 *
 * \return  This function has no return value.
 *
 *
 * \author A. Mignone (mignone@ph.unito.it)
 * \date   Aug 16, 31, 2012
 *********************************************************************** */
#define dEx_dyp(k,j,i) (emf->exj[k][j][i] - EX(k,j,i))
#define dEx_dzp(k,j,i) (emf->exk[k][j][i] - EX(k,j,i))

#define dEy_dxp(k,j,i) (emf->eyi[k][j][i] - EY(k,j,i))
#define dEy_dzp(k,j,i) (emf->eyk[k][j][i] - EY(k,j,i))

#define dEz_dxp(k,j,i) (emf->ezi[k][j][i] - EZ(k,j,i))
#define dEz_dyp(k,j,i) (emf->ezj[k][j][i] - EZ(k,j,i))

#define dEx_dym(k,j,i) (EX(k,j,i) - emf->exj[k][j-1][i])
#define dEx_dzm(k,j,i) (EX(k,j,i) - emf->exk[k-1][j][i])

#define dEy_dxm(k,j,i) (EY(k,j,i) - emf->eyi[k][j][i-1])
#define dEy_dzm(k,j,i) (EY(k,j,i) - emf->eyk[k-1][j][i])

#define dEz_dxm(k,j,i) (EZ(k,j,i) - emf->ezi[k][j][i-1])
#define dEz_dym(k,j,i) (EZ(k,j,i) - emf->ezj[k][j-1][i])
{
  int    i, j, k;
  int    iu, ju, ku;
  signed char  sx, sy, sz;
  double ***vx, ***vy, ***vz;
  double ***Bx, ***By, ***Bz;

  #ifdef CTU
   if (g_intStage == 1) return; /* -- not needed in predictor step of CTU -- */
  #endif

  EXPAND(vx = d->Vc[VX1]; Bx = d->Vc[BX1];   ,
         vy = d->Vc[VX2]; By = d->Vc[BX2];   ,
         vz = d->Vc[VX3]; Bz = d->Vc[BX3];)

  for (k = emf->kbeg; k <= emf->kend + KOFFSET; k++){
  for (j = emf->jbeg; j <= emf->jend + 1      ; j++){
  for (i = emf->ibeg; i <= emf->iend + 1      ; i++){      

    D_EXPAND(sx = emf->svx[k][j][i];  ,
             sy = emf->svy[k][j][i];  ,
             sz = emf->svz[k][j][i];)

    D_EXPAND(iu = sx > 0 ? i:i+1;  ,  /* -- upwind index -- */
             ju = sy > 0 ? j:j+1;  ,
             ku = sz > 0 ? k:k+1;)

 /* ----------------------------------------
      Span X - Faces:    dEz/dy, dEy/dz
    ---------------------------------------- */

    if (sx == 0) {
      emf->ez[k][j][i]   += 0.5*(dEz_dyp(k,j,i) + dEz_dyp(k,j,i+1));
      emf->ez[k][j-1][i] -= 0.5*(dEz_dym(k,j,i) + dEz_dym(k,j,i+1));
      #if DIMENSIONS == 3
       emf->ey[k][j][i]   += 0.5*(dEy_dzp(k,j,i) + dEy_dzp(k,j,i+1));
       emf->ey[k-1][j][i] -= 0.5*(dEy_dzm(k,j,i) + dEy_dzm(k,j,i+1));
      #endif
    }else{
      emf->ez[k][j][i]   += dEz_dyp(k,j,iu);
      emf->ez[k][j-1][i] -= dEz_dym(k,j,iu);
      #if DIMENSIONS == 3
       emf->ey[k][j][i]   += dEy_dzp(k,j,iu);
       emf->ey[k-1][j][i] -= dEy_dzm(k,j,iu);
      #endif
    }

 /* ----------------------------------------
      Span Y - Faces:    dEz/dx, dEx/dz
    ---------------------------------------- */

    if (sy == 0) {
      emf->ez[k][j][i]   += 0.5*(dEz_dxp(k,j,i) + dEz_dxp(k,j+1,i));
      emf->ez[k][j][i-1] -= 0.5*(dEz_dxm(k,j,i) + dEz_dxm(k,j+1,i));
      #if DIMENSIONS == 3
       emf->ex[k][j][i]   += 0.5*(dEx_dzp(k,j,i) + dEx_dzp(k,j+1,i));
       emf->ex[k-1][j][i] -= 0.5*(dEx_dzm(k,j,i) + dEx_dzm(k,j+1,i));
      #endif
    }else{
      emf->ez[k][j][i]   += dEz_dxp(k,ju,i);
      emf->ez[k][j][i-1] -= dEz_dxm(k,ju,i);
      #if DIMENSIONS == 3
       emf->ex[k][j][i]   += dEx_dzp(k,ju,i);
       emf->ex[k-1][j][i] -= dEx_dzm(k,ju,i);
      #endif
    }

 /* ----------------------------------------
      Span Z - Faces:    dEx/dy, dEy/dx
    ---------------------------------------- */

    #if DIMENSIONS == 3

    if (sz == 0) {
      emf->ex[k][j][i]   += 0.5*(dEx_dyp(k,j,i) + dEx_dyp(k+1,j,i));
      emf->ex[k][j-1][i] -= 0.5*(dEx_dym(k,j,i) + dEx_dym(k+1,j,i));
      emf->ey[k][j][i]   += 0.5*(dEy_dxp(k,j,i) + dEy_dxp(k+1,j,i));
      emf->ey[k][j][i-1] -= 0.5*(dEy_dxm(k,j,i) + dEy_dxm(k+1,j,i));
    }else{
      emf->ex[k][j][i]   += dEx_dyp(ku,j,i);
      emf->ex[k][j-1][i] -= dEx_dym(ku,j,i);
      emf->ey[k][j][i]   += dEy_dxp(ku,j,i);
      emf->ey[k][j][i-1] -= dEy_dxm(ku,j,i);
    }

    #endif
  }}}
}
#undef EX
#undef EY
#undef EZ
#undef dEx_dyp
#undef dEx_dzp
#undef dEy_dxp
#undef dEy_dzp
#undef dEz_dxp
#undef dEz_dyp
#undef dEx_dym
#undef dEx_dzm
#undef dEy_dxm
#undef dEy_dzm
#undef dEz_dxm
#undef dEz_dym

/* ********************************************************************* */
void CT_EMF_HLL_Solver (const Data *d, const EMF *emf, Grid *grid)
/*!
 * \brief   Solve 2D Riemann problem for induction equation.
 * \details
 *
 *  Solve 2-D Riemann problem using the 2D HLL Riemann flux 
 *  formula of Londrillo & Del Zanna (2004), JCP, eq. 56.
 *  Here N, W, E, S refer to the following configuration:
 *  \verbatim
 *                    |
 *                    N
 *                 NW | NE
 *                    | 
 *               --W--+--E--
 *                    | 
 *                 SW | SE
 *                    S
 *                    |
 * \endverbatim
 *
 * \param [in]      d     pointer to PLUTO Data structure
 * \param [in]      grid  pointer to Grid structure;
 *
 * \return  This function has no return value.
 * \author A. Mignone (mignone@ph.unito.it)
 * \date   Aug 16, 2012
 *********************************************************************** */

/* -- 2D interpolation macro for integrating the velocity 
      from the cell center to the edge                     -- */
#define dPP(a,x,y,i,j,k) (a[k][j][i] + 0.5*(d##a##_##d##x[k][j][i] + \
                                            d##a##_##d##y[k][j][i]))

#define dPM(a,x,y,i,j,k) (a[k][j][i] + 0.5*(d##a##_##d##x[k][j][i] - \
                                            d##a##_##d##y[k][j][i]))

#define dMM(a,x,y,i,j,k) (a[k][j][i] - 0.5*(d##a##_##d##x[k][j][i] + \
                                            d##a##_##d##y[k][j][i]))

#define dMP(a,x,y,i,j,k) (a[k][j][i] - 0.5*(d##a##_##d##x[k][j][i] - \
                                            d##a##_##d##y[k][j][i]))

/* -- 1D interpolation macro for integrating the staggered  
      magnetic field from the face to the edge              -- */
#define dP(a,x,i,j,k) (a[k][j][i] + 0.5*(d##a##_##d##x[k][j][i]))
#define dM(a,x,i,j,k) (a[k][j][i] - 0.5*(d##a##_##d##x[k][j][i]))

{
  int i,j,k;
  int ip, jp, kp;
  double a_xp, a_yp, a_zp;
  double a_xm, a_ym, a_zm;

  double eSE, eSW, eNW, eNE;
  double bS, bN, bW, bE;

  double *x,  *y,  *z;
  double *xr, *yr, *zr;
  double B0[3];

  double ***bx, ***by, ***bz; 
  double ***vx, ***vy, ***vz;
  double   ***dvx_dx, ***dvx_dy, ***dvx_dz;
  double   ***dvy_dx, ***dvy_dy, ***dvy_dz;
  double   ***dvz_dx, ***dvz_dy, ***dvz_dz;

  double   ***dbx_dy, ***dbx_dz;
  double   ***dby_dx, ***dby_dz;
  double   ***dbz_dx, ***dbz_dy;

  D_EXPAND(vx = d->Vc[VX1]; bx = d->Vs[BX1s];  ,
           vy = d->Vc[VX2]; by = d->Vs[BX2s];  ,
           vz = d->Vc[VX3]; bz = d->Vs[BX3s];)

  dvx_dx = emf->dvx_dx; dvx_dy = emf->dvx_dy; dvx_dz = emf->dvx_dz;
  dvy_dx = emf->dvy_dx; dvy_dy = emf->dvy_dy; dvy_dz = emf->dvy_dz;
  dvz_dx = emf->dvz_dx; dvz_dy = emf->dvz_dy; dvz_dz = emf->dvz_dz;
 
                        dbx_dy = emf->dbx_dy; dbx_dz = emf->dbx_dz;
  dby_dx = emf->dby_dx;                       dby_dz = emf->dby_dz;
  dbz_dx = emf->dbz_dx; dbz_dy = emf->dbz_dy;                      

  x  = grid[IDIR].x;  y  = grid[JDIR].x;  z  = grid[KDIR].x;
  xr = grid[IDIR].xr; yr = grid[JDIR].xr; zr = grid[KDIR].xr;

  for (k = emf->kbeg; k <= emf->kend; k++){  kp = k + 1;
  for (j = emf->jbeg; j <= emf->jend; j++){  jp = j + 1;
  for (i = emf->ibeg; i <= emf->iend; i++){  ip = i + 1;

  /* ------------------------------------------- 
        EMF: Z component at (i+1/2, j+1/2, k)
     ------------------------------------------- */

    a_xp = MAX(emf->SxR[k][j][i], emf->SxR[k][jp][i]);
    a_xm = MAX(emf->SxL[k][j][i], emf->SxL[k][jp][i]);
    a_yp = MAX(emf->SyR[k][j][i], emf->SyR[k][j][ip]);
    a_ym = MAX(emf->SyL[k][j][i], emf->SyL[k][j][ip]);

    bS = dP(bx, y, i , j , k); bW = dP(by, x, i , j , k);
    bN = dM(bx, y, i , jp, k); bE = dM(by, x, ip, j , k);
 
    #if BACKGROUND_FIELD == YES
     BackgroundField (xr[i], yr[j], z[k], B0);
     bS += B0[0]; bW += B0[1];
     bN += B0[0]; bE += B0[1];
    #endif

    eSW = dPP(vy,x,y,i ,j ,k)*bS - dPP(vx,x,y,i ,j ,k)*bW;
    eSE = dMP(vy,x,y,ip,j ,k)*bS - dMP(vx,x,y,ip,j ,k)*bE;
    eNE = dMM(vy,x,y,ip,jp,k)*bN - dMM(vx,x,y,ip,jp,k)*bE;
    eNW = dPM(vy,x,y,i ,jp,k)*bN - dPM(vx,x,y,i ,jp,k)*bW;

    emf->ez[k][j][i] =   a_xp*a_yp*eSW + a_xm*a_yp*eSE 
                       + a_xm*a_ym*eNE + a_xp*a_ym*eNW;
    emf->ez[k][j][i] /= (a_xp + a_xm)*(a_yp + a_ym);
    emf->ez[k][j][i] -= a_yp*a_ym*(bN - bS)/(a_yp + a_ym);
    emf->ez[k][j][i] += a_xp*a_xm*(bE - bW)/(a_xp + a_xm);

  /* ------------------------------------------- 
        EMF: X component at (i, j+1/2, k+1/2)
     ------------------------------------------- */

    #if DIMENSIONS == 3
     a_xp = MAX(emf->SyR[k][j][i], emf->SyR[kp][j][i]);
     a_xm = MAX(emf->SyL[k][j][i], emf->SyL[kp][j][i]);
     a_yp = MAX(emf->SzR[k][j][i], emf->SzR[k][jp][i]);
     a_ym = MAX(emf->SzL[k][j][i], emf->SzL[k][jp][i]);

     bS = dP(by, z, i, j, k ); bW = dP(bz, y, i, j , k);
     bN = dM(by, z, i, j, kp); bE = dM(bz, y, i, jp, k);

     #if BACKGROUND_FIELD == YES
      BackgroundField (x[i], yr[j], zr[k], B0);
      bS += B0[1]; bW += B0[2];
      bN += B0[1]; bE += B0[2];
     #endif
 
     eSW = dPP(vz,y,z,i,j ,k )*bS - dPP(vy,y,z,i ,j ,k )*bW;
     eSE = dMP(vz,y,z,i,jp,k )*bS - dMP(vy,y,z,i ,jp,k )*bE;
     eNE = dMM(vz,y,z,i,jp,kp)*bN - dMM(vy,y,z,i ,jp,kp)*bE;
     eNW = dPM(vz,y,z,i,j ,kp)*bN - dPM(vy,y,z,i ,j ,kp)*bW;  

     emf->ex[k][j][i] =   a_xp*a_yp*eSW + a_xm*a_yp*eSE 
                        + a_xm*a_ym*eNE + a_xp*a_ym*eNW;
     emf->ex[k][j][i] /= (a_xp + a_xm)*(a_yp + a_ym);
     emf->ex[k][j][i] -= a_yp*a_ym*(bN - bS)/(a_yp + a_ym);
     emf->ex[k][j][i] += a_xp*a_xm*(bE - bW)/(a_xp + a_xm);

  /* ------------------------------------------- 
        EMF: Y component at (i+1/2, j, k+1/2)
     ------------------------------------------- */

     a_xp = MAX(emf->SzR[k][j][i], emf->SzR[k][j][ip]);
     a_xm = MAX(emf->SzL[k][j][i], emf->SzL[k][j][ip]);
     a_yp = MAX(emf->SxR[k][j][i], emf->SxR[kp][j][i]);
     a_ym = MAX(emf->SxL[k][j][i], emf->SxL[kp][j][i]);

     bS = dP(bz, x, i , j, k); bW = dP(bx, z, i, j, k);
     bN = dM(bz, x, ip, j, k); bE = dM(bx, z, i, j, kp);

     #if BACKGROUND_FIELD == YES
      BackgroundField (xr[i], y[j], zr[k], B0);
      bS += B0[2]; bW += B0[0];
      bN += B0[2]; bE += B0[0];
     #endif
 
     eSW = dPP(vx,z,x,i, j,k )*bS - dPP(vz,z,x,i ,j ,k )*bW;
     eSE = dMP(vx,z,x,i, j,kp)*bS - dMP(vz,z,x,i ,j ,kp)*bE;
     eNE = dMM(vx,z,x,ip,j,kp)*bN - dMM(vz,z,x,ip,j ,kp)*bE;
     eNW = dPM(vx,z,x,ip,j,k )*bN - dPM(vz,z,x,ip,j ,k )*bW;

     emf->ey[k][j][i] =   a_xp*a_yp*eSW + a_xm*a_yp*eSE 
                        + a_xm*a_ym*eNE + a_xp*a_ym*eNW;
     emf->ey[k][j][i] /= (a_xp + a_xm)*(a_yp + a_ym);
     emf->ey[k][j][i] -= a_yp*a_ym*(bN - bS)/(a_yp + a_ym);
     emf->ey[k][j][i] += a_xp*a_xm*(bE - bW)/(a_xp + a_xm);
    #endif  /* DIMENSIONS == 3 */
  }}}

}
#undef dPP
#undef dPM
#undef dMM
#undef dMP
#undef dP
#undef dM

/* ********************************************************************* */
void CT_EMF_CMUSCL_Average (const Data *d, const EMF *emf, Grid *grid)
/*!
 * \brief   --
 * \details
 *
 *   Used in the predictor scheme of CTU scheme in conjunction 
 *   with UCT_HLL average. 
 *   No riemann solver actually used, but only a simple
 *   pointwise average.
 *   
 * \b References:
 *   - "A High order Godunov scheme with constrained transport and 
 *      adaptie mesh refinement for astrophysical magnetohydrodynamics"\n
 *      Fromang, Hennebelle and Teyssier, A&A (2006) 457, 371
 *
 * \return  This function has no return value.
 * \author A. Mignone (mignone@ph.unito.it)
 * \date   Aug 16, 2012
 *********************************************************************** */
{
  int i,j,k;
  int ip, jp, kp;
  double bx2, by2, bz2;
  double vx4, vy4, vz4;
  double ***vx, ***vy, ***vz;
  double ***bx, ***by, ***bz;

  D_EXPAND(vx = d->Vc[VX1]; bx = d->Vs[BX1s];  ,
           vy = d->Vc[VX2]; by = d->Vs[BX2s];  ,
           vz = d->Vc[VX3]; bz = d->Vs[BX3s];)

  for (k = emf->kbeg; k <= emf->kend; k++){  kp = k + 1;
  for (j = emf->jbeg; j <= emf->jend; j++){  jp = j + 1;
  for (i = emf->ibeg; i <= emf->iend; i++){  ip = i + 1;

  /* ------------------------------------------- 
                 EMF: Z component
     ------------------------------------------- */

    vx4 = 0.25*(vx[k][j][i] + vx[k][j][ip] + vx[k][jp][i] + vx[k][jp][ip]); 
    vy4 = 0.25*(vy[k][j][i] + vy[k][j][ip] + vy[k][jp][i] + vy[k][jp][ip]); 

    bx2 = 0.5*(bx[k][j][i] + bx[k][jp][i]);
    by2 = 0.5*(by[k][j][i] + by[k][j][ip]);

    emf->ez[k][j][i] = vy4*bx2 - vx4*by2;

  /* ------------------------------------------- 
                 EMF: X component
     ------------------------------------------- */

    #if DIMENSIONS == 3
     vy4 = 0.25*(vy[k][j][i] + vy[k][jp][i] + vy[kp][j][i] + vy[kp][jp][i]); 
     vz4 = 0.25*(vz[k][j][i] + vz[k][jp][i] + vz[kp][j][i] + vz[kp][jp][i]); 

     by2 = 0.5*(by[k][j][i] + by[kp][j][i]);
     bz2 = 0.5*(bz[k][j][i] + bz[k][jp][i]);

     emf->ex[k][j][i] = vz4*by2 - vy4*bz2;

  /* ------------------------------------------- 
                 EMF: Y component
     ------------------------------------------- */

     vz4 = 0.25*(vz[k][j][i] + vz[kp][j][i] + vz[k][j][ip] + vz[kp][j][ip]); 
     vx4 = 0.25*(vx[k][j][i] + vx[kp][j][i] + vx[k][j][ip] + vx[kp][j][ip]); 

     bz2 = 0.5*(bz[k][j][i] + bz[k][j][ip]);
     bx2 = 0.5*(bx[k][j][i] + bx[kp][j][i]);

     emf->ey[k][j][i] = vx4*bz2 - vz4*bx2;

    #endif
  }}}

}
