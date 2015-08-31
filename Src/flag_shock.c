/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Shock finding algorithm.
  
  Search and flag computational zones lying in a shock wave.
  The flagging strategy is based on two switches designed to detect 
  the presence of compressive motion or shock waves in the fluid:
  \f[
     \nabla\cdot\vec{v} < 0 \qquad{\rm and}\qquad
     \Delta x\frac{|\nabla p|}{p} > \epsilon_p
  \f]  
  where \f$\epsilon_p\f$ sets the shock strength.
  At the discrete level we replace the two conditions by 
  \f[
    \sum_d \frac{ A_{\vec{i}+\HALF\hvec{e}_d}v_{d,\vec{i}+\HALF\hvec{e}_d}
                 -A_{\vec{i}-\HALF\hvec{e}_d}v_{d,\vec{i}-\HALF\hvec{e}_d} }
                {\Delta{\cal V}_{d,\vec{i}}}  < 0
                \qquad{\rm and}\qquad
    \sum_{d} \left|p_{\vec{i}+\hvec{e}_d} - p_{\vec{i}-\hvec{e}_d}\right|
             < 
     \epsilon_p \min_d\left(p_{\vec{i}+\hvec{e}_d},
                            p_{\vec{i}-\hvec{e}_d},p_{\vec{i}}\right)
  \f]  
  where \f$\hvec{i} = (i,j,k)\f$ is a vector of integer numbers 
  giving the position of a computational zone, while \f$\hvec{e}_d =
  (\delta_{1d},\delta_{2d},\delta_{3d})\f$ is a unit vector in the direction
  given by \c d.
  Once a zone has been tagged as lying in a shock, different flags may be
  switched on or off to control the update strategy in these critical regions.

  This function can be called called when:
  - \c SHOCK_FLATTENING has been set to \c MULTID: in this case shocked zones
    are tagged with \c FLAG_MINMOD and \c FLAG_HLL that will later
    be used to force the reconstruction with the minmod limiter
    and the Riemann solver to HLL.
  - \c ENTROPY_SWITCH has been turned on: this flag will be checked later in
    the ConsToPrim() functions in order to recover pressure from the
    entropy density rather than from the total energy density.
    The update process is:
 
    - start with a vector of primitive variables  <tt> {V,s} </tt>
      where \c s is the entropy;
    - set boundary condition on \c {V}; 
    -  compute \c {s} from \c {V};
    - flag zones where entropy may be used (flag = 1);
    - evolve equations for one time step;
    - convert <tt> {U,S} </tt> to primitive:
      \code
        if (flag == 0) {  // Use energy
          p = p(E)
          s = s(p)
        }else{            // use entropy
          p = p(S)
          E = E(p)
        }  
        \endcode

  \b Reference
     - "Maintaining Pressure Positivity in Magnetohydrodynamics Simulations"
       Balsara \& Spicer, JCP (1999) 148, 133
  
  \authors A. Mignone (mignone@ph.unito.it)
  \date    July 28, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#ifndef EPS_PSHOCK_FLATTEN
 #define EPS_PSHOCK_FLATTEN  5.0
#endif

#ifndef EPS_PSHOCK_ENTROPY
 #define EPS_PSHOCK_ENTROPY  0.05
#endif

#if (SHOCK_FLATTENING == MULTID) || ENTROPY_SWITCH
/* *************************************************************** */
void FlagShock (const Data *d, Grid *grid)
/*!
 *
 *  
 ***************************************************************** */
{
  int  i, j, k, nv;
  int  ip, jp, kp;
  double divv, gradp, pt_min;
  double dpx1, *dx1, *dV1, ***vx1, pt_min1, dvx1;
  double dpx2, *dx2, *dV2, ***vx2, pt_min2, dvx2;
  double dpx3, *dx3, *dV3, ***vx3, pt_min3, dvx3;
  
  double *dVx, *dVy, *dVz;
  double *Ar, *Ath, *r, *th, s;
  static double ***pt;

/* ----------------------------------------------------
   0. Define pointers to variables and allocate memory
   ---------------------------------------------------- */

  dx1 = grid[IDIR].dx; dV1 = grid[IDIR].dV;
  dx2 = grid[JDIR].dx; dV2 = grid[JDIR].dV; 
  dx3 = grid[KDIR].dx; dV3 = grid[KDIR].dV; 
  
  Ar  = grid[IDIR].A;
  r   = grid[IDIR].x;
  Ath = grid[JDIR].A;
  th  = grid[JDIR].x;

  EXPAND(vx1 = d->Vc[VX1];  ,
         vx2 = d->Vc[VX2];  ,
         vx3 = d->Vc[VX3];)

  if (pt == NULL) pt = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    
/* --------------------------------------------------------
   1. Compute total pressure and flag, initially all zones
      with ENTROPY_SWITCH.
   -------------------------------------------------------- */

  TOT_LOOP(k,j,i){  
#if EOS == ISOTHERMAL
    pt[k][j][i] = d->Vc[RHO][k][j][i]*g_isoSoundSpeed*g_isoSoundSpeed;
#else
  #if HAVE_ENERGY 
     pt[k][j][i] = d->Vc[PRS][k][j][i];
  #endif    
#endif

#if (ENTROPY_SWITCH == SELECTIVE) || (ENTROPY_SWITCH == ALWAYS)
    d->flag[k][j][i] |= FLAG_ENTROPY;
#endif
  }
 
/* ----------------------------------------------
   2. Track zones lying in a shock
   ---------------------------------------------- */

  for (k = KOFFSET; k < NX3_TOT-KOFFSET; k++){ 
  for (j = JOFFSET; j < NX2_TOT-JOFFSET; j++){ 
  for (i = IOFFSET; i < NX1_TOT-IOFFSET; i++){

  /* -- Compute divergence of velocity -- */
     
    #if GEOMETRY == CARTESIAN

     D_EXPAND(dvx1 = vx1[k][j][i + 1] - vx1[k][j][i - 1];   ,
              dvx2 = vx2[k][j + 1][i] - vx2[k][j - 1][i];   ,
              dvx3 = vx3[k + 1][j][i] - vx3[k - 1][j][i];)

    #elif GEOMETRY == CYLINDRICAL

     D_EXPAND(dvx1 =   Ar[i]  *(vx1[k][j][i + 1] + vx1[k][j][i])
                     - Ar[i-1]*(vx1[k][j][i - 1] + vx1[k][j][i]);   ,
              dvx2 = vx2[k][j + 1][i] - vx2[k][j - 1][i];           , 
              dvx3 = vx3[k + 1][j][i] - vx3[k - 1][j][i];)

    #elif GEOMETRY == POLAR

     D_EXPAND(dvx1 =  Ar[i]  *(vx1[k][j][i + 1] + vx1[k][j][i])
                    - Ar[i-1]*(vx1[k][j][i - 1] + vx1[k][j][i]);  ,
              dvx2 = (vx2[k][j + 1][i] - vx2[k][j - 1][i])/r[i];      ,
              dvx3 =  vx3[k + 1][j][i] - vx3[k - 1][j][i];)

    #elif GEOMETRY == SPHERICAL

     D_EXPAND(dvx1 =  Ar[i]  *(vx1[k][j][i + 1] + vx1[k][j][i])
                    - Ar[i-1]*(vx1[k][j][i - 1] + vx1[k][j][i]);               ,
              dvx2 = (  Ath[j] *(vx2[k][j + 1][i] + vx2[k][j][i])
                     - Ath[j-1]*(vx2[k][j - 1][i] + vx2[k][j][i]))/fabs(r[i]); ,
              s    = th[j];
              dvx3 = (vx3[k + 1][j][i] - vx3[k - 1][j][i])/(r[i]*sin(s));)

    #endif

    divv = D_EXPAND(dvx1/dV1[i], + dvx2/dV2[j], + dvx3/dV3[k]);
 
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

      #if SHOCK_FLATTENING == MULTID
       if (gradp > EPS_PSHOCK_FLATTEN*pt_min) {
         d->flag[k][j][i]   |= FLAG_HLL;
        
         d->flag[k][j][i]   |= FLAG_MINMOD;
         D_EXPAND(
           d->flag[k][j][i+1] |= FLAG_MINMOD;
           d->flag[k][j][i-1] |= FLAG_MINMOD;  ,
           d->flag[k][j-1][i] |= FLAG_MINMOD;  
           d->flag[k][j+1][i] |= FLAG_MINMOD;  ,
           d->flag[k-1][j][i] |= FLAG_MINMOD;
           d->flag[k+1][j][i] |= FLAG_MINMOD;)
       }
      #endif

  /* -----------------------------------------------------
      When using entropy, we unflag those zones lying in 
      a shock as well as one neighbour cells to the left
      and to the right for each dimension.
     ----------------------------------------------------- */

      #if ENTROPY_SWITCH == SELECTIVE
       if (gradp > EPS_PSHOCK_ENTROPY*pt_min) { /* -- unflag zone -- */
         d->flag[k][j][i] &= ~(FLAG_ENTROPY);
         D_EXPAND(
           d->flag[k][j][i-1] &= ~(FLAG_ENTROPY);
           d->flag[k][j][i+1] &= ~(FLAG_ENTROPY); ,
           d->flag[k][j+1][i] &= ~(FLAG_ENTROPY);
           d->flag[k][j-1][i] &= ~(FLAG_ENTROPY); ,
           d->flag[k+1][j][i] &= ~(FLAG_ENTROPY);
           d->flag[k-1][j][i] &= ~(FLAG_ENTROPY);
         )
       }
      #endif
    } 
      
  }}}

  #ifdef PARALLEL
   AL_Exchange (d->flag[0][0], SZ_char);
  #endif
}

#undef  EPS_PSHOCK_FLATTEN 
#undef  EPS_PSHOCK_ENTROPY
#endif
