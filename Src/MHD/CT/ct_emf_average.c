/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief  Collects different EMF averaging schemes.

  \author A. Mignone
  \date   March 1, 2017
*/  
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void CT_EMF_ArithmeticAverage (const EMF *Z1, const double w)
/*!
 *  Compute arithmetic average of EMF at cell edges by combining
 *  the four electric field values computed at zone faces 
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
 * \return  This function has no return value.
 *********************************************************************** */
#define DEX_DYP(k,j,i) (emf->exj[k][j][i] - d->Ex1[k][j][i])
#define DEX_DZP(k,j,i) (emf->exk[k][j][i] - d->Ex1[k][j][i])

#define DEY_DXP(k,j,i) (emf->eyi[k][j][i] - d->Ex2[k][j][i])
#define DEY_DZP(k,j,i) (emf->eyk[k][j][i] - d->Ex2[k][j][i])

#define DEZ_DXP(k,j,i) (emf->ezi[k][j][i] - d->Ex3[k][j][i])
#define DEZ_DYP(k,j,i) (emf->ezj[k][j][i] - d->Ex3[k][j][i])

#define DEX_DYM(k,j,i) (d->Ex1[k][j][i] - emf->exj[k][j-1][i])
#define DEX_DZM(k,j,i) (d->Ex1[k][j][i] - emf->exk[k-1][j][i])

#define DEY_DXM(k,j,i) (d->Ex2[k][j][i] - emf->eyi[k][j][i-1])
#define DEY_DZM(k,j,i) (d->Ex2[k][j][i] - emf->eyk[k-1][j][i])

#define DEZ_DXM(k,j,i) (d->Ex3[k][j][i] - emf->ezi[k][j][i-1])
#define DEZ_DYM(k,j,i) (d->Ex3[k][j][i] - emf->ezj[k][j-1][i])
{
  int    i, j, k;
  int    iu, ju, ku;
  signed char  sx, sy, sz;

#ifdef CTU
  if (g_intStage == 1) return; /* -- not needed in predictor step of CTU -- */
#endif

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
      emf->ez[k][j][i]   += 0.5*(DEZ_DYP(k,j,i) + DEZ_DYP(k,j,i+1));
      emf->ez[k][j-1][i] -= 0.5*(DEZ_DYM(k,j,i) + DEZ_DYM(k,j,i+1));
#if DIMENSIONS == 3
      emf->ey[k][j][i]   += 0.5*(DEY_DZP(k,j,i) + DEY_DZP(k,j,i+1));
      emf->ey[k-1][j][i] -= 0.5*(DEY_DZM(k,j,i) + DEY_DZM(k,j,i+1));
#endif
    }else{
      emf->ez[k][j][i]   += DEZ_DYP(k,j,iu);
      emf->ez[k][j-1][i] -= DEZ_DYM(k,j,iu);
#if DIMENSIONS == 3
      emf->ey[k][j][i]   += DEY_DZP(k,j,iu);
      emf->ey[k-1][j][i] -= DEY_DZM(k,j,iu);
#endif
    }

 /* ----------------------------------------
      Span Y - Faces:    dEz/dx, dEx/dz
    ---------------------------------------- */

    if (sy == 0) {
      emf->ez[k][j][i]   += 0.5*(DEZ_DXP(k,j,i) + DEZ_DXP(k,j+1,i));
      emf->ez[k][j][i-1] -= 0.5*(DEZ_DXM(k,j,i) + DEZ_DXM(k,j+1,i));
#if DIMENSIONS == 3
      emf->ex[k][j][i]   += 0.5*(DEX_DZP(k,j,i) + DEX_DZP(k,j+1,i));
      emf->ex[k-1][j][i] -= 0.5*(DEX_DZM(k,j,i) + DEX_DZM(k,j+1,i));
#endif
    }else{
      emf->ez[k][j][i]   += DEZ_DXP(k,ju,i);
      emf->ez[k][j][i-1] -= DEZ_DXM(k,ju,i);
#if DIMENSIONS == 3
      emf->ex[k][j][i]   += DEX_DZP(k,ju,i);
      emf->ex[k-1][j][i] -= DEX_DZM(k,ju,i);
#endif
    }

 /* ----------------------------------------
      Span Z - Faces:    dEx/dy, dEy/dx
    ---------------------------------------- */

#if DIMENSIONS == 3
    if (sz == 0) {
      emf->ex[k][j][i]   += 0.5*(DEX_DYP(k,j,i) + DEX_DYP(k+1,j,i));
      emf->ex[k][j-1][i] -= 0.5*(DEX_DYM(k,j,i) + DEX_DYM(k+1,j,i));
      emf->ey[k][j][i]   += 0.5*(DEY_DXP(k,j,i) + DEY_DXP(k+1,j,i));
      emf->ey[k][j][i-1] -= 0.5*(DEY_DXM(k,j,i) + DEY_DXM(k+1,j,i));
    }else{
      emf->ex[k][j][i]   += DEX_DYP(ku,j,i);
      emf->ex[k][j-1][i] -= DEX_DYM(ku,j,i);
      emf->ey[k][j][i]   += DEY_DXP(ku,j,i);
      emf->ey[k][j][i-1] -= DEY_DXM(ku,j,i);
    }
#endif
  }}}
}
#undef DEX_DYP
#undef DEX_DZP
#undef DEY_DXP
#undef DEY_DZP
#undef DEZ_DXP
#undef DEZ_DYP
#undef DEX_DYM
#undef DEX_DZM
#undef DEY_DXM
#undef DEY_DZM
#undef DEZ_DXM
#undef DEZ_DYM

/* ********************************************************************* */
void CT_EMF_HLL_Solver (const Data *d, const EMF *emf, Grid *grid)
/*!
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

  x  = grid->x[IDIR];  y  = grid->x[JDIR];  z  = grid->x[KDIR];
  xr = grid->xr[IDIR]; yr = grid->xr[JDIR]; zr = grid->xr[KDIR];

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

