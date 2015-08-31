/* ///////////////////////////////////////////////////////////////////// */
/*! \file
 *  \brief Fill normal staggered component in the ghost zones.           */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************** */
void FillMagneticField (const Data *d, int side, Grid *grid)
/*!
 * \brief Assign the normal component of the staggered magnetic field
 *        in the ghost zone-faces.
 *
 * \details
 *  
 *  Using the div.B = 0 condition in staggered MHD, this function 
 *  computes the staggered component of magnetic field lying on the 
 *  zone-face parallel to the boundary specified by "side".
 *  This is preformed by solving the div.B = 0 condition for one 
 *  variable only which in 2-D requires the knowledge of the other 
 *  3 components while in 3-D required the knowledge of the other 5 
 *  staggered components.\n
 * 
 *  Note that this operation is performed in the outermost ghost zones 
 *  only since the face value at IBEG-1 or IEND is supposed to be part 
 *  of the solution and is not changed during this function.
 *  Therefore, only nghost-1 faces are assigned:
 *
 *  \verbatim
 *  +-----+-----+-----+-----+-----+--
 *  |     |     |     |     |     |
 *  |     X     X     |     |     |
 *  |     |     |     |     |     |
 *  +-----+-----+-----+-----+-----+--
 *                    |
 *  <-----------------> BEG
 *   Physical boundary     
 *
 *   X = components assigned in this function.
 * \endverbatim
 *
 * \param [in,out]    d  pointer to PLUTO Data structure
 * \param [in]     side  the side
 * \param [in]     grid  pointer to PLUTO Grid structure
 *
 * \todo   replace the loops with more compact macro, such as
 *         X1_BEG_LOOP()...
 * \author A. Mignone (mignone@ph.unito.it)
 * \date   Aug 16, 2012
 *********************************************************************** */
{
  int  ibeg, jbeg, kbeg;
  int  iend, jend, kend;
  int  i,j,k;
  int  di,dj,dk;

  real r, dmu;
  real bxp, byp, bzp;
  real bxm, bym, bzm;
  real dBx = 0.0, dBy = 0.0, dBz = 0.0;
  real dVx, dVy, dVz;
  real  dx,  dy,  dz;
  real rp, sp, Axp, Ayp, Azp;
  real rm, sm, Axm, Aym, Azm;

  double ***bx, ***by, ***bz;  /* -- staggered magnetic field -- */
  double ***Bx, ***By, ***Bz;  /* -- cell centered mag. field -- */

  D_EXPAND(bx = d->Vs[BX1s]; Bx = d->Vc[BX1];  ,
           by = d->Vs[BX2s]; By = d->Vc[BX2];  ,
           bz = d->Vs[BX3s]; Bz = d->Vc[BX3];)

  ibeg = 0; iend = NX1_TOT-1; di = 1;
  jbeg = 0; jend = NX2_TOT-1; dj = 1;
  #if DIMENSIONS == 3
   kbeg = 0; kend = NX3_TOT-1; dk = 1;
  #else
   kbeg = kend = 0, dk = 1;
  #endif

  if (side == X1_BEG) {ibeg = IBEG-1; iend = 0; di = -1;}
  if (side == X1_END)  ibeg = IEND+1;

  if (side == X2_BEG) {jbeg = JBEG-1; jend = 0; dj = -1;}
  if (side == X2_END)  jbeg = JEND+1;

  if (side == X3_BEG) {kbeg = KBEG-1; kend = 0; dk = -1;}
  if (side == X3_END)  kbeg = KEND+1;

  for (k = kbeg; dk*k <= dk*kend; k += dk){
  for (j = jbeg; dj*j <= dj*jend; j += dj){
  for (i = ibeg; di*i <= di*iend; i += di){

    r  = grid[IDIR].x[i];
    rp = fabs(grid[IDIR].xr[i]);
    rm = fabs(grid[IDIR].xl[i]);

    dx = grid[IDIR].dx[i];
    dy = grid[JDIR].dx[j];
    dz = grid[KDIR].dx[k];

    #if GEOMETRY == SPHERICAL
     dmu = grid[JDIR].dV[j];
     sp  = grid[JDIR].A[j];
     sm  = grid[JDIR].A[j - 1];
    #endif

    #if DIMENSIONS == 2
     dz = 1.0;
    #endif

    D_EXPAND(bxp = bx[k][j][i]; bxm = bx[k][j][i - 1];  ,
             byp = by[k][j][i]; bym = by[k][j - 1][i];  ,
             bzp = bz[k][j][i]; bzm = bz[k - 1][j][i];)

  /* -------------------------------------------------------
       Divergence is written as 

         dVx*(Axp*bxp - Axm*dxm) + 
         dVy*(Ayp*byp - Aym*dym) + 
         dVz*(Azp*bzp - Azm*dzm) = 0

       so that the k-th component can be 
       recovered as
 
      bkp = bkm*Akm/Akp + 
            sum_(j != k) (Ajp*bjp - Ajm*bjm)*dVj/(dVk*Akp)
    ------------------------------------------------------- */

    #if GEOMETRY == CARTESIAN

     dVx = dy*dz; Axp = Axm = 1.0;
     dVy = dx*dz; Ayp = Aym = 1.0;
     dVz = dx*dy; Azp = Azm = 1.0;

    #elif GEOMETRY == CYLINDRICAL

     dVx =   dy*dz; Axp = rp ; Axm = rm;
     dVy = r*dx*dz; Ayp = 1.0; Aym = 1.0;
     dVz =   dx*dy; Azp = 1.0; Azm = 1.0;

    #elif GEOMETRY == POLAR

     dVx =   dy*dz; Axp = rp ; Axm = rm;
     dVy =   dx*dz; Ayp = 1.0; Aym = 1.0;
     dVz = r*dx*dy; Azp = 1.0; Azm = 1.0;

    #elif GEOMETRY == SPHERICAL

     dVx =  dmu*dz; Axp = rp*rp; Axm = rm*rm;
     dVy = r*dx*dz; Ayp = sp   ; Aym = sm;
     dVz = r*dx*dy; Azp = 1.0  ; Azm = 1.0;

    #endif

    D_EXPAND(dBx = dVx*(Axp*bxp - Axm*bxm);  ,
             dBy = dVy*(Ayp*byp - Aym*bym);  ,
             dBz = dVz*(Azp*bzp - Azm*bzm); )

/* -------------------------------------------------
      Assign a single face magnetic field 
   ------------------------------------------------- */

    if (side == X1_BEG){

      bxm = (bxp*Axp + (dBy + dBz)/dVx)/Axm;
      bx[k][j][i - 1] = bxm;

    }else if (side == X1_END){

      bxp = (bxm*Axm - (dBy + dBz)/dVx)/Axp;
      bx[k][j][i] = bxp;

    }else if (side == X2_BEG){

      bym = (byp*Ayp + (dBx + dBz)/dVy)/Aym;
      by[k][j - 1][i] = bym;

    }else if (side == X2_END){
  
      byp = (bym*Aym - (dBx + dBz)/dVy)/Ayp;
      by[k][j][i] = byp;

    #if DIMENSIONS == 3
     }else if (side == X3_BEG){

       bzm = (bzp*Azp + (dBx + dBy)/dVz)/Azm;
       bz[k - 1][j][i] = bzm;

     }else if (side == X3_END){

       bzp = (bzm*Azm - (dBx + dBy)/dVz)/Azp;
       bz[k][j][i] = bzp;

    #endif
    }

  /* -------------------------------------------------------
        Now redefine the cell-centered magnetic field 
     ------------------------------------------------------- */
/*
    #if GEOMETRY == CARTESIAN 
     D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  ,
               By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);  ,
               Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]); )
    #elif GEOMETRY == CYLINDRICAL
     D_EXPAND( Bx[k][j][i] = 0.5*(bx[k][j][i] + bx[k][j][i - 1]);  ,
               By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);  ,
               Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]); )
    #elif GEOMETRY == POLAR
     D_EXPAND(Bx[k][j][i] = (Axp*bx[k][j][i] + Axm*bx[k][j][i - 1])/(Axp + Axm);  ,
              By[k][j][i] = 0.5*(by[k][j][i] + by[k][j - 1][i]);                  ,
              Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]);)
    #elif GEOMETRY == SPHERICAL
     D_EXPAND(Bx[k][j][i] = (Axp*bx[k][j][i] + Axm*bx[k][j][i - 1])/(Axp + Axm);   ,
              By[k][j][i] = (Ayp*by[k][j][i] + Aym*by[k][j - 1][i])/(Ayp + Aym);  ,
              Bz[k][j][i] = 0.5*(bz[k][j][i] + bz[k - 1][j][i]);)
    #endif
*/
  }}}

}


