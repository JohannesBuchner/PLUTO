/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Integrate cooling and reaction source terms.

  Solve the system of ordinary differential equations describing
  optically thin cooling and ionization network (if any).
  We use an adaptive step size integration method and follow the 
  strategy outlined in section 5.1 of Tesileanu et al. (2008) for handling
  stiff cells.

  On output, a time-step estimate for the next time level is computed using
  the relative or absolute variation obtained during the integration of the ODE 
  (from t(n) --> t(n+1))  
  \f[
     \Delta t_c = \min_{ijk}\left[\frac{\Delta t^n M_R}{\epsilon}\right]
     \,\quad\rm{where}\quad
     \epsilon = \max\left(\left|\frac{p^{n+1}}{p^n} - 1\right|,\,
                |X^{n+1}-X^n|\right)
  \f]
  where \f$ M_R \f$ is the maximum cooling rate (defined by the global variable  
  ::g_maxCoolingRate) and X are the chemical species.
  
  \b References
     - "Simulating radiative astrophysical flows with the PLUTO code:
        a non-equilibrium, multi-species cooling function" \n
       Tesileanu, Mignone and Massaglia, A&A (2008) 488, 429

  \authors A. Mignone (mignone@ph.unito.it)\n
           O. Tesileanu
           B. Vaidya
  \date    April 22, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void CoolingSource (const Data *d, double dt, Time_Step *Dts, Grid *GXYZ)
/*!
 * Integrate cooling and reaction source terms.
 *
 * \param [in,out]  d   pointer to Data structure
 * \param [in]     dt   the time step to be taken
 * \param [out]    Dts  pointer to the Time_Step structure
 * \param [in]    GXYZ  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int  nv, k, j, i, stiff, status;
  double err, scrh, min_tol = 2.e-5;
  double mu0, T0, T1, mu1, prs;
  double v0[NVAR], v1[NVAR], k1[NVAR];
  double maxrate;
  intList var_list;
  
/* --------------------------------------------------------
    Set number and indices of the time-dependent variables
   -------------------------------------------------------- */

  var_list.nvar    = NIONS+1;
  var_list.indx[0] = PRS;
  for (nv = 0; nv < NIONS; nv++) var_list.indx[nv+1] = NFLX+nv;

/* -----------------------------------------------------
    Zero k1
   ----------------------------------------------------- */
  
  NVAR_LOOP(nv) k1[nv] = 0.0;  

/*  ----------------------------------------------------------- 
                   Begin Integration 
    -----------------------------------------------------------  */

  DOM_LOOP(k,j,i){  /* -- span the computational domain -- */

  /* --------------------------------------------------
      Skip integration if cell has been tagged with 
      FLAG_INTERNAL_BOUNDARY or FLAG_SPLIT_CELL 
      (only for AMR)
     -------------------------------------------------- */ 

    #if INTERNAL_BOUNDARY == YES
     if (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY) continue;
    #endif
    if (d->flag[k][j][i] & FLAG_SPLIT_CELL) continue;
    
  /* ----------------------------------------------
      Compute temperature and internal energy from
      density, pressure and concentrations.
     ---------------------------------------------- */
    
    NVAR_LOOP(nv) v0[nv] = v1[nv] = d->Vc[nv][k][j][i];
    prs = v0[PRS];
    mu0 = MeanMolecularWeight(v0);
    T0  = v0[PRS]/v0[RHO]*KELVIN*mu0;
    #if EOS == IDEAL
     v0[RHOE] = v1[RHOE] = prs/(g_gamma-1.0);
    #else
     v1[RHOE] = v0[RHOE] = InternalEnergy(v0, T0);
    #endif

    if (T0 <= 0.0){
      print ("! CoolingSource: negative initial temperature\n");
      print (" %12.6e  %12.6e\n",v0[RHOE], v0[RHO]);
      print (" at: %f %f\n",GXYZ[IDIR].x[i], GXYZ[JDIR].x[j]);
      QUIT_PLUTO(1);
    }

  /* -------------------------------------------
      Get estimated time step based on 
      the max rate of the reaction network.
     ------------------------------------------- */

    Radiat(v0, k1);

    maxrate = GetMaxRate (v0, k1, T0);
    stiff = (dt*maxrate > 1.0 ? 1:0);

  /* ---------------------------------------------------
      If the system is not stiff, then try to advance
      with an explicit 2-nd order midpoint rule
     --------------------------------------------------- */

    if (!stiff){ /* -- no stiffness: try explicit -- */
 
      scrh = SolveODE_RKF12 (v0, k1, v1, dt, min_tol, &var_list); 
/*    scrh = SolveODE_RKF23 (v0, k1, v1, dt, min_tol);  */

/* -- error is too big ? --> use some other integrator -- */

      if (scrh < 0.0) SolveODE_CK45 (v0, k1, v1, dt, min_tol, &var_list);

    } else { /* -- if stiff = 1 use more sophisticated integrators -- */

/*  SolveODE_ROS34 (v0, k1, v1, dtsub, min_tol);  */

      int nsub, k;
      double dtsub, dtnew, t;

      nsub  = ceil(dt*maxrate);
      dtsub = dt/(double)nsub;

    /* ----------------------------------
            use sub-time stepping
       ---------------------------------- */

      t = 0.0;
      for (k = 1; 1; k++){
        dtnew = SolveODE_CK45 (v0, k1, v1, dtsub, min_tol, &var_list);
     /* dtnew = SolveODE_ROS34(v0, k1, v1, dtsub, min_tol);  */
        t    += dtsub;
     /* printf ("%d / %d   %12.6e  %12.6e, t/dt = %f\n",
        k,nsub,dtsub,dtnew, t/dt);*/

        if (fabs(t/dt - 1.0) < 1.e-9) break;

        v0[RHOE] = v1[RHOE];
        NIONS_LOOP(nv) v0[nv] = v1[nv];

        Radiat(v0, k1);
        dtsub = MIN (dtnew, dt - t);
      }

      if (k > 100) print ("! CoolingSource: Number of substeps exceeded 100 (%d)\n",k);
      if (fabs(t/dt - 1.0) > 1.e-12) {
        print ("! CoolingSource: dt mismatch\n");
        QUIT_PLUTO(1);
      }
    }  /* -- end if (stiff) -- */

  /* -- Constrain ions to lie between [0,1] -- */

    NIONS_LOOP(nv){
      v1[nv] = MAX(v1[nv], 0.0);
      v1[nv] = MIN(v1[nv], 1.0);
    }
    #if COOLING == H2_COOL
     v1[X_H2] = MIN(v1[X_H2], 0.5);
    #endif
    
  /* -- pressure must be positive -- */

    mu1 = MeanMolecularWeight(v1);
    #if EOS == IDEAL
     prs = v1[RHOE]*(g_gamma - 1.0);
     T1  = prs/v1[RHO]*KELVIN*mu1;
    #elif EOS == PVTE_LAW
     status = GetEV_Temperature(v1[RHOE], v1, &T1);
     prs    = v1[RHO]*T1/(KELVIN*mu1);
    #endif

    if (prs < 0.0) prs = g_smallPressure;

  /* -- Check final temperature -- */


    if (T1 < g_minCoolingTemp && T0 > g_minCoolingTemp)
      prs = g_minCoolingTemp*v1[RHO]/(KELVIN*mu1);

  /* ------------------------------------------
      Suggest next time step based on 
      fractional variaton.
     ------------------------------------------ */

    err = fabs(prs/d->Vc[PRS][k][j][i] - 1.0);

    #if COOLING == MINEq
     for (nv = NFLX; nv < NFLX + NIONS - Fe_IONS; nv++) 
    #else
       NIONS_LOOP(nv)
    #endif
      err = MAX(err, fabs(d->Vc[nv][k][j][i] - v1[nv]));

    scrh = dt*g_maxCoolingRate/err;

    Dts->dt_cool = MIN(Dts->dt_cool, scrh);

  /* ---- Update solution array ---- */

    d->Vc[PRS][k][j][i] = prs;
    NIONS_LOOP(nv) d->Vc[nv][k][j][i] = v1[nv];

  } /* -- end loop on points -- */
}
 
/* ********************************************************************* */
void Numerical_Jacobian (double *v, double **J)
/*!
 *  Compute Jacobian matrix numerically and 
 *  get eigenvector and eigenvalue decomposition
 *
 *  The purpose of this function is to detect
 *  stiffness. 
 *
 *  Note: This function is EXTREMELY time consuming.
 *
 *********************************************************************** */
{
  int  n, nv, k, l;
  double eps;
  double vpp[NVAR], vp[NVAR], vm[NVAR], vmm[NVAR];
  double Rpp[NVAR], Rp[NVAR], Rm[NVAR], Rmm[NVAR];

  eps = 1.e-5;
  n = NIONS + 1;

/* -- partial derivs with respect to pressure -- */

  for (nv = 0; nv < NVAR; nv++){
    vp[nv] = vm[nv] = v[nv];
  }
  vp[PRS] *= 1.0 + eps;
  vm[PRS] *= 1.0 - eps;
  
  Radiat (vp, Rp);
  Radiat (vm, Rm);

  for (k = 0; k < n - 1; k++){
    J[k][n - 1] = (Rp[k + NFLX] - Rm[k + NFLX])/(2.0*v[PRS]*eps);
  }
  J[n - 1][n - 1] = (Rp[PRS] - Rm[PRS])/(2.0*v[PRS]*eps);

/* -- partial deriv with respect to ions -- */

  for (l = 0; l < n - 1; l++){

    for (nv = 0; nv < NVAR; nv++){
      vp[nv] = vm[nv] = v[nv];
    }
    vp [l + NFLX] = v[l + NFLX] + eps;
    vm [l + NFLX] = v[l + NFLX] - eps;
                      
    vp [l + NFLX] = MIN(vp [l + NFLX], 1.0);
    vm [l + NFLX] = MAX(vm [l + NFLX], 0.0);
                         
    Radiat (vp , Rp);
    Radiat (vm , Rm);

    for (k = 0; k < n - 1; k++){
      J[k][l] = (Rp[k + NFLX] - Rm[k + NFLX])/(vp[l + NFLX] - vm[l + NFLX]); 
    }
    J[n - 1][l] = (Rp[PRS] - Rm[PRS])/(vp[l + NFLX] - vm[l + NFLX]);
  }
}
