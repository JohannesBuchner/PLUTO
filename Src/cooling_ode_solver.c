#include "pluto.h"
#if COOLING == MINEq
 #include "cooling_defs.h"
#endif

/* ********************************************************************* */
double SolveODE_CK45 (double *v0, double *k1, double *v5th, 
                      double dt, double tol, intList *vars)
/*
 *
 *     use explicit Cash-Karp 4-5 integrator
 *
 *
 * 
 *********************************************************************** */
{
  int   nv, ksub, kfail;

  const double c1 = 37.0/378.0;
  const double c2 = 0.0;
  const double c3 = 250.0/621.0;
  const double c4 = 125.0/594.0;
  const double c5 = 0.0;
  const double c6 = 512.0/1771.0; 

  const double cs1 = 2825.0/27648.0;
  const double cs2 = 0.0;
  const double cs3 = 18575.0/48384.0;
  const double cs4 = 13525.0/55296.0;
  const double cs5 = 277.0/14336.0;
  const double cs6 = 0.25;

  const double b21 = 0.2;
  const double b31 = 3.0/40.0      , b32 = 9.0/40.0;
  const double b41 = 0.3           , b42 = -0.9       , b43 = 1.2;
  const double b51 = -11.0/54.0    , b52 = 2.5        , b53 = -70.0/27.0   , b54 = 35.0/27.0;
  const double b61 = 1631.0/55296.0, b62 = 175.0/512.0, b63 = 575.0/13824.0, b64 = 44275.0/110592.0, b65 = 253.0/4096.0;
    
  double scrh, err;
  double t, tend, dt_shrink, dt_grow;
  double v1[NVAR], v4th[NVAR];
  double k2[NVAR], k3[NVAR], k4[NVAR], k5[NVAR], k6[NVAR];
  double vscal[NVAR];
  double tsub[4096];

/* -------------------------------------------
    copy ALL variables so that density is 
    defined when calling Radiat.
   ------------------------------------------- */

  for (nv = 0; nv < NVAR; nv++) v1[nv] = v0[nv];

  t    = tsub[0] = 0.0;
  tend = dt;

  ksub = kfail = 0;
  for (;;){  

    FOR_EACH(nv,0,vars) v1[nv] = v0[nv] + dt*b21*k1[nv];

  /* -- get K2 -- */

    Radiat(v1, k2);
    FOR_EACH(nv,0,vars) v1[nv] = v0[nv] + dt*(b31*k1[nv] + b32*k2[nv]);

  /* -- get K3 -- */

    Radiat(v1, k3);
    FOR_EACH(nv, 0, vars){
      v1[nv] = v0[nv] + dt*(b41*k1[nv] + b42*k2[nv] + b43*k3[nv]);
    }

  /* -- get K4 -- */

    Radiat(v1, k4);
    FOR_EACH(nv,0, vars) {
      v1[nv] = v0[nv] + dt*(b51*k1[nv] + b52*k2[nv] + b53*k3[nv]  
                          + b54*k4[nv]);
    }

  /* -- get K5 -- */

    Radiat(v1, k5);
    FOR_EACH(nv,0,vars){
      v1[nv] = v0[nv] + dt*(b61*k1[nv] + b62*k2[nv] + b63*k3[nv] 
                          + b64*k4[nv] + b65*k5[nv]);
    }

  /* -- get K6 -- */

    Radiat(v1, k6);

  /* -- 5th order solution -- */

    FOR_EACH(nv,0,vars){
      v5th[nv] = v0[nv] + dt*(c1*k1[nv] + c2*k2[nv] + c3*k3[nv] 
                            + c4*k4[nv] + c5*k5[nv] + c6*k6[nv]);
    }

  /* -- 4th order embedded solution -- */

    FOR_EACH(nv,0,vars){
      v4th[nv] = v0[nv] + dt*(cs1*k1[nv] + cs2*k2[nv] + cs3*k3[nv] 
                            + cs4*k4[nv] + cs5*k5[nv] + cs6*k6[nv]);
    }

  /* -----------------------------------
            compute error
     ----------------------------------- */

    vscal[PRS] = fabs(v0[PRS]) + dt*fabs(k1[PRS]);
    
    NIONS_LOOP(nv)  vscal[nv] = 1.0;  /* fabs(v0[nv]) + dt*fabs(k1[nv]) + 1.e-6; */    
    

   /* -- do not take error on Fe -- */

    #if COOLING == MINEq && Fe_IONS > 0
      vscal[FeI] = vscal[FeII] = vscal[FeIII] = 1.e18;
    #endif

    err = 0.0;
    FOR_EACH(nv,0,vars){
      scrh = fabs(v5th[nv] - v4th[nv])/fabs(vscal[nv]);
      err  = MAX(err, scrh);
    }
    err /= tol;

    if (err < 1.0){  /* -- ok, accept step -- */

      ksub++;
      err = MAX(0.0, 1.e-18);

      t          += dt;
      tsub[ksub]  = t;

    /* -- provide an estimate for next dt -- */

      dt_grow = 0.9*dt*pow(err, -0.2);
      dt_grow = MIN(dt_grow, 5.0*dt); /* -- do not increase more than 5 -- */
      dt      = dt_grow;

    /* -- exit loop if final time has been reached -- */

      if (fabs(t/tend - 1.0) < 1.e-9) break;

    /* -- adjust dt not to exceed exactly tend -- */

      if (dt > (tend - t)) dt = tend - t;

    /* -- initialize solution vector before continuing -- */
    
      FOR_EACH(nv,0,vars) v0[nv] = v5th[nv];
      Radiat(v0, k1);

      if (ksub > 1000) {
        print ("! SolveODE_CK45: Number of substeps too large (%d)\n",ksub);
        QUIT_PLUTO(1);
      }

    }else{   /* -- shrink dt and redo time step -- */

      kfail++;
      dt_shrink = 0.9*dt*pow(err, -0.25);
      dt        = MAX(dt_shrink, 0.05*dt); /* -- do not decrease more than 20 -- */
    }

  } /* -- end while -- */

  if (ksub > 100) {
    int i;
    print ("! SolveODE_CK45: number of substeps is %d\n", ksub);
  }

  return (dt);
}    

/* ********************************************************************* */
double SolveODE_RK4 (double *v0, double *k1, double *v4th,
                     double dt, intList *var_list)
/*!
 *  Solve the system of ODE with a standard RK4 integrator (no
 *  adaptive step).
 *  
 *
 *
 * 
 *********************************************************************** */
{
  int   nv;
  double v1[NVAR];
  double k2[NVAR], k3[NVAR], k4[NVAR];

/* -------------------------------------------
    copy ALL variables so that density is 
    defined when calling Radiat.
   ------------------------------------------- */

  for (nv = 0; nv < NVAR; nv++) v1[nv] = v0[nv];

/* -- step 1 -- */

  FOR_EACH(nv,0,var_list) v1[nv] = v0[nv] + 0.5*dt*k1[nv];

/* -- step 2 -- */

  Radiat(v1, k2);
  FOR_EACH(nv,0,var_list) v1[nv] = v0[nv] + 0.5*dt*k2[nv];

/* -- step 3 -- */

  Radiat(v1, k3);
  FOR_EACH(nv,0,var_list) v1[nv] = v0[nv] + 0.5*dt*k3[nv];

/* -- step 4 -- */

  Radiat(v1, k4);
  FOR_EACH(nv,0,var_list) {
    v4th[nv] = v0[nv] + dt*(k1[nv] + 2.0*(k2[nv] + k3[nv]) + k4[nv])/6.0;
  }

  return (0.0);
}

/* ********************************************************************* */
double SolveODE_RKF12 (double *v0, double *k1, double *v2nd, 
                       double dt, double tol, intList *vars)
/*
 *
 *
 * 
 *********************************************************************** */
{
  int  nv;
  double err, scrh, dt_grow;
  double v1[NVAR], v1st[NVAR], vscal[NVAR];
  double k2[NVAR];
  
/* -------------------------------------------
    copy ALL variables so that density is 
    defined when calling Radiat.
   ------------------------------------------- */

  for (nv = 0; nv < NVAR; nv++) v1[nv] = v0[nv];

/* -- Get K2 -- */

  FOR_EACH(nv, 0, vars)  v1[nv] = v0[nv] + 0.5*dt*k1[nv];

  Radiat(v1, k2);

/* -- 2nd order solution -- */

  FOR_EACH(nv, 0, vars) v2nd[nv] = v0[nv] + dt*k2[nv];

/* -- 1st order solution -- */

  FOR_EACH(nv, 0, vars) v1st[nv] = v0[nv] + dt*k1[nv];

/* -----------------------------------
          compute error
   ----------------------------------- */

  vscal[PRS] = fabs(v0[PRS]) + dt*fabs(k1[PRS]);
 
  NIONS_LOOP(nv) vscal[nv] = 1.0;  /* fabs(v0[nv]) + dt*fabs(k1[nv]) + 1.e-6; */    
 

  err = 0.0;
  FOR_EACH(nv,0,vars){
    scrh = fabs(v2nd[nv] - v1st[nv])/fabs(vscal[nv]);
    err  = MAX(err, scrh);
  }
  err /= tol;

  if (err < 1.0) return ( 1.0);
  else           return (-1.0);
}

/* ********************************************************************* */
double SolveODE_RKF23 (double *v0, double *k1, double *v3rd, 
                       double dt, double tol, intList *vars)
/*
 *
 *
 * y(3rd) = y(0) + dt*(k1 + k2 + 4*k3)/6
 *
 *  k1 = RHS(y(0))
 *  k2 = RHS(y(0) + dt*k1)
 *  k3 = RHS(y(0) + 0.25*dt*(k1 + k2))
 *
 * y(2nd) = y(0) + dt*(k1 + k2)/2
 *
 *********************************************************************** */
{
  int  nv;
  double err, scrh, dt_4, dt_6;
  double v1[NVAR], v2nd[NVAR], vscal[NVAR];
  double k2[NVAR], k3[NVAR];

  dt_4 = 0.25*dt;
  dt_6 = dt/6.0;

/* -------------------------------------------
    copy ALL variables so that density is 
    defined when calling Radiat.
   ------------------------------------------- */

  for (nv = 0; nv < NVAR; nv++) v1[nv] = v0[nv];

/* -- Get K2 -- */

  FOR_EACH(nv,0,vars) v1[nv] = v0[nv] + dt*k1[nv]; 

  Radiat(v1, k2);

/* -- Get K3 -- */

  FOR_EACH(nv,0,vars) v1[nv] = v0[nv] + dt_4*(k1[nv] + k2[nv]);

  Radiat(v1, k3);

/* -- 3rd order solution -- */

  FOR_EACH(nv,0,vars) v3rd[nv] = v0[nv] + dt_6*(k1[nv] + k2[nv] + 4.0*k3[nv]);

/* -- 2nd order solution -- */

  FOR_EACH(nv,0,vars) v2nd[nv] = v0[nv] + 0.5*dt*(k1[nv] + k2[nv]);

/* -----------------------------------
          compute error
   ----------------------------------- */

  vscal[PRS] = fabs(v0[PRS]) + dt*fabs(k1[PRS]);
  
  NIONS_LOOP(nv) vscal[nv] = 1.0;  /* fabs(v0[nv]) + dt*fabs(k1[nv]) + 1.e-6; */    
  

  err = 0.0;
  FOR_EACH(nv,0,vars){
    scrh = fabs(v3rd[nv] - v2nd[nv])/fabs(vscal[nv]);
    err  = MAX(err, scrh);
  }
  err /= tol;

  if (err < 1.0) return ( 1.0);
  else           return (-1.0);
}

#define GAM (1.0/2.0)
#define A21 2.0
#define A31 (48.0/25.0)
#define A32 (6.0/25.0)
#define C21 -8.0
#define C31 (372.0/25.0)
#define C32 (12.0/5.0)
#define C41 (-112.0/125.0)
#define C42 (-54.0/125.0)
#define C43 (-2.0/5.0)
#define B1 (19.0/9.0)
#define B2 (1.0/2.0)
#define B3 (25.0/108.0)
#define B4 (125.0/108.0)
#define E1 (17.0/54.0)
#define E2 (7.0/36.0)
#define E3 0.0
#define E4 (125.0/108.0)

#define C1X (1.0/2.0)
#define C2X (-3.0/2.0)
#define C3X (121.0/50.0)
#define C4X (29.0/250.0)
#define A2X 1.0
#define A3X (3.0/5.0)

/* ********************************************************************* */
double SolveODE_ROS34 (double *v0, double *k1, double *v4th, 
                       double dt, double tol)
/*
 * 
 *   Solve using Semi-Implicit Rosenbrock Method
 *
 *********************************************************************** */
{
  int    i, j, n, nv, ksub, kfail;
  static int *indx;
  double   err, scrh, t, tend;
  double   v1[NVAR], vscal[NVAR], k2[NVAR], k3[NVAR];
  double   tsub[4096], dt_grow, dt_shrink;
  static double  **a, **J, **J2;
  static double  *g1, *g2, *g3, *g4;

double vbeg[NVAR];

/* -------------------------------------------
    copy ALL variables so that density is
    defined when calling Radiat.
   ------------------------------------------- */
                                                                                                                             
  for (nv = 0; nv < NVAR; nv++) v1[nv] = vbeg[nv] = v0[nv];
  
  n = NIONS + 1;
  
  if (indx == NULL){
    indx  = ARRAY_1D (n, int);
    a     = ARRAY_2D (n, n, double);
    J     = ARRAY_2D (n, n, double);
    J2    = ARRAY_2D (n, n, double);
    g1    = ARRAY_1D (n, double);
    g2    = ARRAY_1D (n, double);
    g3    = ARRAY_1D (n, double);
    g4    = ARRAY_1D (n, double);
  }

  for (i = 0; i < n - 1; i++) vscal[i + NFLX] = 1.0;

 /* -- do not take error on Fe -- */

  #if COOLING == MINEq && Fe_IONS > 0
    vscal[FeI] = vscal[FeII] = vscal[FeIII] = 1.e18;
  #endif

  vscal[PRS] = fabs(v0[PRS]);

  Jacobian (v0, k1, J);   
/*
  Numerical_Jacobian (v0, J2); 
{
 int k,l;
 double T0;

 T0 = v0[PRS]/v0[RHO]*KELVIN*MeanMolecularWeight(v0);
 printf (" T0 = %12.6e\n",T0);

 for (k = 0; k < n; k++){

 for (l = 0; l < n; l++){
   printf ("%4.1f ", 100.*fabs(J2[k][l] - J[k][l])/(fabs(J[k][l]) + 1.e-12));
  }

//  for (l = 0; l < 1; l++){
//    printf ("%12.6e   %12.6e", J2[k][l], J[k][l]);
//  }
  printf ("\n");
 }
 exit(1);
}
*/
  t    = 0.0;
  tend = dt;

  ksub = kfail = 0;
  for (;;){

  /* --  Compute (I - hJ)  -- */

    for (i = 0; i < n; i++) {   
      for (j = 0; j < n; j++) a[i][j] = -J[i][j];
      a[i][i] += 1.0/(GAM*dt);
    }
    LUDecompose (a, n, indx, &scrh);    /*    LU decomposition of the matrix. */

  /* -- set right hand side for g1 -- */

    for (i = 0; i < n - 1; i++) {  
      g1[i] = k1[i + NFLX];
    }
    g1[n - 1] = k1[PRS];

  /* -- solve for g1 -- */

    LUBackSubst (a, n, indx, g1); 
    for (i = 0; i < n - 1; i++) {  
      v1[i + NFLX] = v0[i + NFLX] + A21*g1[i];
    }
    v1[PRS] = v0[PRS] + A21*g1[n - 1];

    Radiat (v1, k2);    

  /* -- set right hand side for g2 -- */

    for (i = 0; i < n - 1; i++) { 
      g2[i] = k2[i + NFLX] + C21*g1[i]/dt;
    }
    g2[n - 1] = k2[PRS] + C21*g1[n - 1]/dt;

  /* -- solve for g2 -- */

    LUBackSubst (a, n, indx, g2);  
    for (i = 0; i < n - 1; i++) {    
      v1[i + NFLX] = v0[i + NFLX] + A31*g1[i] + A32*g2[i];
    }
    v1[PRS] = v0[PRS] + A31*g1[n - 1] + A32*g2[n - 1];

    Radiat (v1, k3);
    
  /* -- set right hand side for g3 -- */

    for (i = 0; i < n - 1; i++) { 
      g3[i] = k3[i + NFLX] + (C31*g1[i] + C32*g2[i])/dt;
    }
    g3[n - 1] = k3[PRS] + (C31*g1[n - 1] + C32*g2[n - 1])/dt;

  /* -- solve for g3 -- */

    LUBackSubst (a, n, indx, g3);  

  /* -- set right hand side for g4 -- */

    for (i = 0; i < n - 1; i++) { 
      g4[i] = k3[i + NFLX] + (C41*g1[i] + C42*g2[i] + C43*g3[i])/dt;
    }
    g4[n - 1] = k3[PRS] + (C41*g1[n - 1] + C42*g2[n - 1] + C43*g3[n - 1])/dt;

  /* -- solve for g4  -- */ 

    LUBackSubst (a, n, indx, g4);    

  /* --  4th order solution & error estimation -- */

    i = n - 1;
    v4th[PRS] = v0[PRS] + B1*g1[i] + B2*g2[i] + B3*g3[i] + B4*g4[i];
    err      = fabs(E1*g1[i] + E2*g2[i] + E3*g3[i] + E4*g4[i])/vscal[PRS];

    for (i = 0; i < n - 1; i++) {   
      v4th[i + NFLX] = v0[i + NFLX] + B1*g1[i] + B2*g2[i] + B3*g3[i] + B4*g4[i];
      scrh           = E1*g1[i] + E2*g2[i] + E3*g3[i] + E4*g4[i];
      err            = MAX(err, fabs(scrh)/vscal[i + NFLX]);
    }
    
    err /= tol;

    if (err < 1.0) {

      ksub++;      
      err = MAX(0.0, 1.e-18);

      t          += dt;
      tsub[ksub]  = t;

    /* -- provide an estimate for next dt -- */

      dt_grow = 0.9*dt*pow(err, -0.25);
      dt_grow = MIN(dt_grow, 5.0*dt);
      dt      = dt_grow;

    /* -- exit loop if final time has been reached -- */

      if (fabs(t/tend - 1.0) < 1.e-9) break;

    /* -- adjust dt not to exceed tend -- */

      if (dt > (tend - t)) dt = tend - t;

    /* -- initialize solution vector continuing -- */
   
      v0[PRS]  = v4th[PRS];
      NIONS_LOOP(nv) v0[nv] = v4th[nv];
      Radiat (v0, k1);
      Jacobian (v0, k1, J);   

      if (ksub > 1000){
        print ("! SolveODE_ROS34: Number of substeps too large (%d)\n",ksub);
        QUIT_PLUTO(1);
      }
                                                                                                                             
    }else{   /* -- shrink dt and redo time step -- */

      kfail++;
      dt_shrink = 0.9*dt*pow(err, -1.0/3.0);
      dt        = MAX(dt_shrink, 0.05*dt);                                                                                                                        
    }
  }

  if (ksub > 2) {
    int i;
    print ("! SolveODE_ROS34: number of substeps is %d\n", ksub);
/*
    for (i = 1; i <= ksub; i++){
      printf ("# %d, dt = %12.6e, t/dt = %f,  tend = %12.6e\n",
           i, tsub[i]-tsub[i-1],tsub[i]/tend, tend);
    }
    printf ("kfail = %d\n",kfail);
    QUIT_PLUTO(1);
*/
  }

  return (dt);
}

#undef GAM 
#undef A21 
#undef A31 
#undef A32 
#undef C21 
#undef C31 
#undef C32 
#undef C41 
#undef C42 
#undef C43 
#undef B1 
#undef B2 
#undef B3 
#undef B4 
#undef E1 
#undef E2 
#undef E3 
#undef E4 
#undef C1X 
#undef C2X 
#undef C3X 
#undef C4X 
#undef A2X 
#undef A3X 


