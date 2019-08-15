/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Riemann solver for dust (pressuless gas)

  Solve the Riemann problem for dust.
  When DUST_FLUID_SOLVER is set to 1, a standard Lax-Frierichs solver is used.
  Otherwise the exact solver from LeVeque is employed. 

  \b Reference:
   - "THE DYNAMICS OF PRESSURELESS DUST_FLUID CLOUDS AND DELTA WAVES"
      R. J. LeVeque, 
      J. Hyper. Differential Equations, 01, 315 (2004).

  \authors A. Mignone (mignone@ph.unito.it)
  \date    April 4, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#define DUST_FLUID_SOLVER  2

/* ********************************************************************* */
void Dust_Solver (const Sweep *sweep, int beg, int end, 
                  double *cmax, Grid *grid)
/*!
 * Solve Riemann problem for the adiabatic/isothermal MHD equations 
 * using the HLL Riemann solver.
 *
 * \param[in,out] sweep   pointer to Sweep structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int    nv, i;
  double *uR, *uL, *vL, *vR;
  double  vmax, fL[NVAR], fR[NVAR];


#if DUST_FLUID_SOLVER == 1

/* -------------------------------------------------------------
    Solver #1: Lax Friedrichs
   ------------------------------------------------------------- */

  for (i = beg; i <= end; i++) {
    uL   = sweep->uL[i]; vL = sweep->vL[i];
    uR   = sweep->uR[i]; vR = sweep->vR[i];
    
    fL[RHO_D] = uL[MXn_D];
    EXPAND(fL[MX1_D] = uL[MX1_D]*vL[VXn_D];  ,
           fL[MX2_D] = uL[MX2_D]*vL[VXn_D];  ,
           fL[MX3_D] = uL[MX3_D]*vL[VXn_D];)

    fR[RHO_D] = uR[MXn_D];
    EXPAND(fR[MX1_D] = uR[MX1_D]*vR[VXn_D];  ,
           fR[MX2_D] = uR[MX2_D]*vR[VXn_D];  ,
           fR[MX3_D] = uR[MX3_D]*vR[VXn_D];)
    
    vmax = MAX(fabs(vL[VXn_D]), fabs(vR[VXn_D]));
    NDUST_FLUID_LOOP(nv){
      sweep->flux[i][nv] = 0.5*(fL[nv] + fR[nv] -  vmax*(uR[nv] - uL[nv]));
    }  
    cmax[i] = MAX(cmax[i], vmax);
  }
#elif DUST_FLUID_SOLVER == 2  

/* -------------------------------------------------------------
    Solver #2: Exact solver (Le Veque 2003)
   ------------------------------------------------------------- */

  for (i = beg; i <= end; i++) {
    vR = sweep->vR[i]; uR = sweep->uR[i];
    vL = sweep->vL[i]; uL = sweep->uL[i];

    if (vL[VXn_D] < 0.0 && vR[VXn_D] > 0.0){
      sweep->flux[i][RHO_D] = 0.0;
      EXPAND(sweep->flux[i][MX1_D] = 0.0;  ,
             sweep->flux[i][MX2_D] = 0.0;  ,
             sweep->flux[i][MX3_D] = 0.0;)
    } else {
      double lambda;
      lambda  = sqrt(vL[RHO_D])*vL[VXn_D] + sqrt(vR[RHO_D])*vR[VXn_D];
      lambda /= sqrt(vL[RHO_D]) + sqrt(vR[RHO_D]);
      if (lambda > 0.0){
        sweep->flux[i][RHO_D] = uL[MXn_D];
        EXPAND(sweep->flux[i][MX1_D] = uL[MX1_D]*vL[VXn_D];  ,
               sweep->flux[i][MX2_D] = uL[MX2_D]*vL[VXn_D];  ,
               sweep->flux[i][MX3_D] = uL[MX3_D]*vL[VXn_D];)
      }else if (lambda < 0.0){
        sweep->flux[i][RHO_D] = uR[MXn_D];
        EXPAND(sweep->flux[i][MX1_D] = uR[MX1_D]*vR[VXn_D];  ,
               sweep->flux[i][MX2_D] = uR[MX2_D]*vR[VXn_D];  ,
               sweep->flux[i][MX3_D] = uR[MX3_D]*vR[VXn_D];)
      }else{
        sweep->flux[i][RHO_D] = 0.5*(uL[MXn_D] + uR[MXn_D]);
        EXPAND(sweep->flux[i][MX1_D] = 0.5*(uL[MX1_D]*vL[VXn_D] + uR[MX1_D]*vR[VXn_D]);  ,
               sweep->flux[i][MX2_D] = 0.5*(uL[MX2_D]*vL[VXn_D] + uR[MX2_D]*vR[VXn_D]);  ,
               sweep->flux[i][MX3_D] = 0.5*(uL[MX3_D]*vL[VXn_D] + uR[MX3_D]*vR[VXn_D]);)
      }
    }
    vmax    = MAX(fabs(vL[VXn_D]), fabs(vR[VXn_D]));
    cmax[i] = MAX(cmax[i], vmax);
  }
#endif
}
#undef DUST_FLUID_SOLVER

/* ********************************************************************* */
void Dust_DragForce(const Sweep *sweep, int beg, int end, double dt, Grid *grid)

/*!
 *  Dust-gas coupling coefficients (Drug force)
 * ********************************************************************* */
{
  int i;
  double fd[COMPONENTS], k=1.0, r;
  double **rhs = sweep->rhs;
  double **v   = sweep->v;
  
  for (i = beg; i <= end; i++){
    r = grid[IDIR].x[i];
    k = v[i][RHO_D]/1000.0;
    fd[g_dir] = -k*(v[i][VXn] - v[i][VXn_D]);
    rhs[i][MXn]   += dt*fd[g_dir]; 
    rhs[i][MXn_D] -= dt*fd[g_dir];  
  }

/*    
  for (i = beg; i <= end; i++){
    k = v[i][RHO_D]*g_isoSoundSpeed/10.0;  
    EXPAND(fd[IDIR] = -k*(v[i][VX1] - v[i][VX1_D]);  ,
           fd[JDIR] = -k*(v[i][VX2] - v[i][VX2_D]);  ,
           fd[KDIR] = -k*(v[i][VX3] - v[i][VX3_D]);)

    EXPAND(rhs[i][MX1] += dt*fd[IDIR];  ,
           rhs[i][MX2] += dt*fd[JDIR];  ,
           rhs[i][MX3] += dt*fd[KDIR];)

    EXPAND(rhs[i][MX1_D] -= dt*fd[IDIR];  ,
           rhs[i][MX2_D] -= dt*fd[JDIR];  ,
           rhs[i][MX3_D] -= dt*fd[KDIR];)
  }
*/
}
/* ********************************************************************* */
void Dust_DragForceImpliciUpdate()
/*
 *
 *  U -->
 *
 *********************************************************************** */
{
  int i,j,k,nv;

/* -----------------------------------------------------
   1. Compute drag matrix
   ----------------------------------------------------- */


/* -----------------------------------------------------
   2. Update solution
   ----------------------------------------------------- */
/*
  DOM_LOOP(k,j,i){
    d->Vc[nv][k][j][i] = exp(IK)*d->Vc[nv][k][j][i]; (something like that)
  }
*/

}
