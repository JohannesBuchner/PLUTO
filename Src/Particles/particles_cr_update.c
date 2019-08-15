/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Update CR particles.
 
  Boris pusher for updating Cosmic ray particles.
  \code
   dt = dt_g/Nsub
    for n = 0,Nsub-1 {
     Compute Fcr
     Loop on particles{
       x(n+1/2) = x(n)     + (dt/2)*v(n)     (drift)
       v(n+1)   = v(n)     +  dt*a(n+1/2)    (kick+rotate+kick)
       x(n+1)   = x(n+1/2) + (dt/2)*v(n+1)   (drift)
        save dL += dt*(v(n) + v(n+1))/2
       Compute time step:
         dL < Nz*dx --> dt < Nz*dx/[v(0)/2 + v(1) + .. V(Nsub-1)/2)]
     }
   }
  \endcode
 
  Time step restriction is computed by requiring that no particle
  travels more than \c Nmax = \c CR_NCELLS_EPS zones and that
  the Larmor scale is resolved with more than 1 cycle:
  \f[
   \left\{
     \begin{array}{lcl}
       \Delta s_d &=& \DS\Delta t\max_{n,p}\left(
                      \frac{|v_{p,d}^n + v_{p,d}^{n+1}|}{2}\right)
                   < \epsilon_s N_{\max}\Delta x_d \\ \noalign{\medskip}
       \Delta t  &<& \epsilon_L \Omega_L^{-1}
     \end{array}
   \right.
  \f]
  where the maximum extends to all sub-steps (\c n) and particles
  (p), \f$ \Omega_L = qB_\perp/(\gamma m c) \f$ is the Larmor frequency while
  \f$\epsilon_s (\sim 0.9)\f$ and \f$\epsilon_L\sim 0.3\f$ are safety factors
  while 
  \f[
   \vec{B}^2_\perp = \vec{B}^2 - \frac{(\vec{v}\cdot\vec{B})^2}
                                      {\vec{v}\cdot\vec{v}}                                
  \f]
  
  \authors A. Mignone (mignone@ph.unito.it)\n

  \b References
    - "MAGNETOHYDRODYNAMIC-PARTICLE-IN-CELL METHOD FOR COUPLING COSMIC RAYS 
       WITH A THERMAL PLASMA: APPLICATION TO NON-RELATIVISTIC SHOCKS"\n
       Bai et al., ApJ (2015) 809, 55

  \date   April 10, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static void Particles_CR_Predictor(Data *, timeStep *, double, Grid *);

#if (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_TEST == NO)
/* ********************************************************************* */
void Particles_CR_Update(Data *data, timeStep *Dts, double dt, Grid *grid)
/*!
 * Advance particle by a step dt.
 * 
 * \param [in,out]  d     Data structure (contains paeticles list)
 * \param [in,out]  Dts   timeStep structure
 * \param [in]      dt    time time increment
 * \param [in]      grid  pointer to Grid structure
 *********************************************************************** */
{
  int    i,j,k, dir;
  int    kcycle;
  long int dsize = NX1_TOT*NX2_TOT*NX3_TOT;
  double C = 1.e8;  /* Used only when PARTICLES_DEPOSIT == INTEGER */
  
  double pcoord_old[3], pspeed_old[3];
  double vg[3], B[3], E[3];
  
  double um[3], up[3], u1[3], u[3], b[3], u_old[3];
  double gamma, gamma_old, Bperp;
  double b2, Bmag2, u2, u2_old, v2; /* Squares of vector */
  double qd[4], scrh, omL;
  double dt0, dt_half, h_half, inv_dt, qg;
  double inv_dtL=0.0, inv_dts=0.0;
  double wF, wM;
  const double c2 = PARTICLES_CR_C*PARTICLES_CR_C;
  static double ***w, ****dM, ****dM_tot, ****Fcr0;
  static double ****emf0, ****emf_res;

  particleNode *CurNode;
  Particle *p;

  DEBUG_FUNC_BEG ("CR_Update");

#if SHOW_TIMING
  clock_t clock_beg = clock(), clock0;
#endif

  Boundary (data, ALL_DIR, grid);

/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */
  
  if (w == NULL){
    w      = ArrayBox (-1, 1, -1, 1, -1, 1);
    emf0   = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
    #if RESISTIVITY != NO
    emf_res = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
    #endif
    dM_tot = ARRAY_4D(4, NX3_TOT, NX2_TOT, NX1_TOT, double);
    dM     = ARRAY_4D(4, NX3_TOT, NX2_TOT, NX1_TOT, double);
    Fcr0   = ARRAY_4D(4, NX3_TOT, NX2_TOT, NX1_TOT, double);
  }

#if PARTICLES_CR_FEEDBACK == YES
  for (dir = 0; dir < 4; dir++) {
    TOT_LOOP(k,j,i) {
      dM_tot[dir][k][j][i] = 0.0;
      #if TIME_STEPPING == RK2                   /* Save Fcr at t^n */
      Fcr0[dir][k][j][i] = data->Fcr[dir][k][j][i]; /* for time-interpolation */
      #endif                                     
    }
  }
#endif

/* --------------------------------------------------------
   1. Initialize sub-cycling
   -------------------------------------------------------- */

  inv_dt  = 1.e-18;
  dt0     = dt;  /* Save initial dt */
  dt      = dt/(double)Dts->Nsub_particles;
  dt_half = 0.5*dt;
  h_half  = 0.5*dt*PARTICLES_CR_E_MC; 

/* --------------------------------------------------------
   2a. Compute convective electric field, cE0 = -v x B 
   -------------------------------------------------------- */

  TOT_LOOP(k,j,i){
    vg[IDIR] = data->Vc[VX1][k][j][i]; B[IDIR] = data->Vc[BX1][k][j][i];
    vg[JDIR] = data->Vc[VX2][k][j][i]; B[JDIR] = data->Vc[BX2][k][j][i];
    vg[KDIR] = data->Vc[VX3][k][j][i]; B[KDIR] = data->Vc[BX3][k][j][i];
  
    emf0[IDIR][k][j][i] = -(vg[JDIR]*B[KDIR] - vg[KDIR]*B[JDIR]);
    emf0[JDIR][k][j][i] = -(vg[KDIR]*B[IDIR] - vg[IDIR]*B[KDIR]);
    emf0[KDIR][k][j][i] = -(vg[IDIR]*B[JDIR] - vg[JDIR]*B[IDIR]);
  }

#if PARTICLES_CR_FEEDBACK == NO
  for (dir = 0; dir < 3; dir++) TOT_LOOP(k,j,i){       
     data->Ecr[dir][k][j][i] = emf0[dir][k][j][i];
  }
  
  /* -----------------------------------------------------
      Compute resistive electric field at cell centers
     ----------------------------------------------------- */
  
  #if RESISTIVITY != NO
  for (k = KOFFSET; k < NX3_TOT-KOFFSET; k++){
  for (j = JOFFSET; j < NX2_TOT-JOFFSET; j++){
  for (i = IOFFSET; i < NX1_TOT-IOFFSET; i++){
    int nv;
    double Vfluid[NVAR], J[3], eta[3];
    double *x1 = grid->x[IDIR], *dx1 = grid->dx[IDIR];
    double *x2 = grid->x[JDIR], *dx2 = grid->dx[JDIR];
    double *x3 = grid->x[KDIR], *dx3 = grid->dx[KDIR];
    
    #if PHYSICS == MHD
    J[IDIR] =  0.5*(data->Vc[BX3][k][j+1][i] - data->Vc[BX3][k][j-1][i])/dx2[j];
    J[JDIR] = -0.5*(data->Vc[BX3][k][j][i+1] - data->Vc[BX3][k][j][i-1])/dx1[i];
    J[KDIR] =  0.5*(data->Vc[BX2][k][j][i+1] - data->Vc[BX2][k][j][i-1])/dx1[i]
              -0.5*(data->Vc[BX1][k][j+1][i] - data->Vc[BX1][k][j-1][i])/dx2[j];
    
    NVAR_LOOP(nv) Vfluid[nv] = data->Vc[nv][k][j][i];
    
    /* -- Compute current and resistivity at cell center -- */
  
    Resistive_eta (Vfluid, x1[i], x2[j], x3[k], J, eta);
    emf_res[IDIR][k][j][i] = eta[IDIR]*J[IDIR];
    emf_res[JDIR][k][j][i] = eta[JDIR]*J[JDIR];
    emf_res[KDIR][k][j][i] = eta[KDIR]*J[KDIR];
    #elif PHYSICS == RMHD
    emf_res[IDIR][k][j][i] = data->Vc[EX1][k][j][i] - emf0[IDIR][k][j][i];
    emf_res[JDIR][k][j][i] = data->Vc[EX2][k][j][i] - emf0[JDIR][k][j][i];
    emf_res[KDIR][k][j][i] = data->Vc[EX3][k][j][i] - emf0[KDIR][k][j][i];
    #endif
    
  }}}
  
  /* ------------------------------------------------------
      In parallel, we need to fill ghost zone values for
      the resistive electric field since it is computed
      using central differences
     ------------------------------------------------------ */
  
  #ifdef PARALLEL
  MPI_Barrier (MPI_COMM_WORLD);
  AL_Exchange ((char *)emf_res[IDIR][0][0], SZ);
  AL_Exchange ((char *)emf_res[JDIR][0][0], SZ);
  AL_Exchange ((char *)emf_res[KDIR][0][0], SZ);
  #endif
  #endif /* RESISTIVITY != NO */
  
#else

  #if PARTICLES_CR_PREDICTOR == NO

  /* -- Recompute Ecr since vXB has been evolved to t^{n+1/2} -- */

  TOT_LOOP(k,j,i){       
    qg = PARTICLES_CR_E_MC_GAS*data->Vc[RHO][k][j][i];
    for (dir = 0; dir < 3; dir++){
      data->Ecr[dir][k][j][i] = emf0[dir][k][j][i] - data->Fcr[dir][k][j][i]/qg;
    }
  }
  #endif
#endif

/* --------------------------------------------------------
   3. Start sub-cycling
   -------------------------------------------------------- */

  for (kcycle = 0; kcycle < Dts->Nsub_particles; kcycle++){

    int correct_emf = 0;

  /* --------------------------------------------
     3a. Predictor step 
     -------------------------------------------- */

    #if PARTICLES_CR_FEEDBACK == YES
    #if PARTICLES_CR_PREDICTOR == 1
    if (kcycle == 0){
      Particles_CR_Predictor (data, Dts, 0.5*dt, grid);
    }else{
    
    /* -- Compute Fcr(x^{n+1/2}, v^n), extrapolate in time -- */
    
      Particles_CR_ComputeForce(data->Vc, data, grid); 
      for (dir = 0; dir < 3; dir++) TOT_LOOP(k,j,i){
        data->Fcr[dir][k][j][i] = 2.0*data->Fcr[dir][k][j][i] - dM[dir][k][j][i]/dt;
      }
    } /* end if (ncycle != 0) */
    correct_emf = 1;
    #elif PARTICLES_CR_PREDICTOR == 2
    if (Dts->Nsub_particles == 1){
      Particles_CR_Predictor (data, Dts, 0.5*dt, grid);
      correct_emf = 1;
    }else {
      if (kcycle == 0)  {          /* Pre-push particles */
        Particles_CR_Predictor (data, Dts, dt, grid);
        correct_emf = 1;
      }else  if (kcycle%2 == 0){   /* Extrapolate */
        Particles_CR_ComputeForce(data->Vc, data, grid); 
        wF = (kcycle + 2.0)/kcycle;
        wM = 2.0/kcycle;
        for (dir = 0; dir < 3; dir++) TOT_LOOP(k,j,i){
          data->Fcr[dir][k][j][i] =   wF*data->Fcr[dir][k][j][i] 
                                 - wM*dM_tot[dir][k][j][i]/(kcycle*dt);
        }
        correct_emf = 1;
      }
    } /* end if odd cycle */
    #endif  /* PARTICLES_CR_PREDICTOR == 2 */

  /* ------------------------------------------------
     3b. Correct total electric field
         E = E0 - Fcr/qg at x^{n+(k+1/2)/Nsub}
     ------------------------------------------------ */

    if (correct_emf){
      TOT_LOOP(k,j,i){
        qg = PARTICLES_CR_E_MC_GAS*data->Vc[RHO][k][j][i];
        for (dir = 0; dir < 3; dir++){
          data->Ecr[dir][k][j][i] = emf0[dir][k][j][i] - data->Fcr[dir][k][j][i]/qg;
        }
      }
    }

    for (dir = 0; dir < 4; dir++){
      memset ((void *)dM[dir][0][0], '\0', dsize*sizeof(double));
    }
    #endif  /* PARTICLES_CR_FEEDBACK == YES */

  /* ---------------------------------------------
     3c. Boris pusher:
         - Drift by dt/2
         - Kick, rotate by dt
         - Drift by dt/2
        
        At the beginning of this step, coordinate
        and velocity must be set as follows:

        p->coord = x^n
        p->speed = v^n
     -------------------------------------------- */

    PARTICLES_LOOP(CurNode, data->PHead){
      p = &(CurNode->p);

    /* -- A. Get particle 4-velocity and compute \gamma^n -- */

      v2    = DOT_PRODUCT(p->speed,p->speed);
      gamma = 1.0/sqrt(1.0 - v2/c2);

      for (dir = 0; dir < 3; dir++) {
        pcoord_old[dir] = p->coord[dir];        /* Save x^n      */
        pspeed_old[dir] = p->speed[dir];       /* Save v^n      */

        p->coord[dir]       += 0.5*dt*p->speed[dir];  /* Drift by dt/2  */
        u[dir] = u_old[dir] = p->speed[dir]*gamma; /* Compute 4-vel  */
      }
      gamma_old = gamma;
      u2_old    = v2*gamma*gamma;

    /* -- B. Compute weights and indices at x^{n+1/2} for later deposition -- */

      Particles_GetWeights(p, p->cell, w, grid);  
      i = p->cell[IDIR];
      j = p->cell[JDIR];
      k = p->cell[KDIR];

    /* -- C. Interpolate electromagnetic fields at x^{n+1/2} -- */
  
      B[IDIR] = Particles_Interpolate(data->Vc[BX1], w, p->cell); 
      B[JDIR] = Particles_Interpolate(data->Vc[BX2], w, p->cell); 
      B[KDIR] = Particles_Interpolate(data->Vc[BX3], w, p->cell);

      E[IDIR] = Particles_Interpolate(data->Ecr[IDIR], w, p->cell); 
      E[JDIR] = Particles_Interpolate(data->Ecr[JDIR], w, p->cell); 
      E[KDIR] = Particles_Interpolate(data->Ecr[KDIR], w, p->cell);

    /* -- D1. Clean parallel component of E(perp), so that E.B = 0 -- */

      Bmag2 = DOT_PRODUCT(B,B);
      scrh  = DOT_PRODUCT(E,B)/(Bmag2 + 1.e-18);

      E[IDIR] -= scrh*B[IDIR];
      E[JDIR] -= scrh*B[JDIR];
      E[KDIR] -= scrh*B[KDIR];

    /* -- D2. Add resistive field after cleaning step -- */
    
      #if RESISTIVITY != NO
      E[IDIR] += Particles_Interpolate(emf_res[IDIR], w, p->cell); 
      E[JDIR] += Particles_Interpolate(emf_res[JDIR], w, p->cell); 
      E[KDIR] += Particles_Interpolate(emf_res[KDIR], w, p->cell);
      #endif
      
    /* -- E. Boris pusher (kick + rotate + kick) -- */
     
      um[IDIR] = u[IDIR] + h_half*E[IDIR];
      um[JDIR] = u[JDIR] + h_half*E[JDIR];
      um[KDIR] = u[KDIR] + h_half*E[KDIR];

      scrh  = DOT_PRODUCT(um,um);
      gamma = sqrt(1.0 + scrh/c2);  /* Compute time-centered \gamma */

      scrh    = h_half/gamma;
      b[IDIR] = scrh*B[IDIR];
      b[JDIR] = scrh*B[JDIR];
      b[KDIR] = scrh*B[KDIR];
    
      b2      = DOT_PRODUCT(b,b);    

      u1[IDIR] = um[IDIR] + (um[JDIR]*b[KDIR] - um[KDIR]*b[JDIR]);
      u1[JDIR] = um[JDIR] + (um[KDIR]*b[IDIR] - um[IDIR]*b[KDIR]);
      u1[KDIR] = um[KDIR] + (um[IDIR]*b[JDIR] - um[JDIR]*b[IDIR]);

      scrh = 2.0/(1.0 + b2);

      up[IDIR] = um[IDIR] + scrh*(u1[JDIR]*b[KDIR] - u1[KDIR]*b[JDIR]);
      up[JDIR] = um[JDIR] + scrh*(u1[KDIR]*b[IDIR] - u1[IDIR]*b[KDIR]);
      up[KDIR] = um[KDIR] + scrh*(u1[IDIR]*b[JDIR] - u1[JDIR]*b[IDIR]);
    
    /* -- F. Update velocity by another half step -- */

      u[IDIR] = up[IDIR] + h_half*E[IDIR];
      u[JDIR] = up[JDIR] + h_half*E[JDIR];
      u[KDIR] = up[KDIR] + h_half*E[KDIR];
    
      u2    = DOT_PRODUCT(u,u);
      gamma = sqrt(1.0 + u2/c2);
      scrh  = 1.0/gamma;

      p->speed[IDIR] = u[IDIR]*scrh;
      p->speed[JDIR] = u[JDIR]*scrh;
      p->speed[KDIR] = u[KDIR]*scrh;

    /* -- G. Deposit momentum and energy variations at x^{n+1/2} -- */

      #if PARTICLES_CR_FEEDBACK == YES
      int i1,j1,k1,n,nelem=4;

      qd[IDIR] = p->rho*(u[IDIR] - u_old[IDIR]);
      qd[JDIR] = p->rho*(u[JDIR] - u_old[JDIR]);
      qd[KDIR] = p->rho*(u[KDIR] - u_old[KDIR]);
      qd[3]    = p->rho*(u2  /(gamma + 1.0) - u2_old/(gamma_old + 1.0));

      for (n = 0; n < nelem; n++){
        for (k1 = -KOFFSET; k1 <= KOFFSET; k1++){
        for (j1 = -JOFFSET; j1 <= JOFFSET; j1++){
        for (i1 = -IOFFSET; i1 <= IOFFSET; i1++){
          #if PARTICLES_DEPOSIT == INTEGER
          long a;
          scrh = qd[n]*w[k1][j1][i1]*C;
          a    = (long)(scrh);
          scrh = (double)(a);
          dM[n][k+k1][j+j1][i+i1] += scrh;
          #else
          dM[n][k+k1][j+j1][i+i1] += qd[n]*w[k1][j1][i1];
          #endif
        }}}
      }
      #endif

    /* -- H. Update spatial coordinate to obtain x^{n+1} -- */

      p->coord[IDIR] += dt_half*p->speed[IDIR];
      p->coord[JDIR] += dt_half*p->speed[JDIR];
      p->coord[KDIR] += dt_half*p->speed[KDIR];

    /* ---------------------------------------------------------
        I1. Compute time step restriction based on the maximum
            allowed distance that a particle can travel at
            its current speed:

            1/dt_1 = v^{n+1/2}/(eps * dx)

            where eps = CR_NCELL_EPS.
       --------------------------------------------------------- */
  
      for (dir = 0; dir < DIMENSIONS; dir++) {
        scrh   = 0.5*(fabs(p->speed[dir]) + fabs(p->speed_old[dir]));
        scrh  /= PARTICLES_CR_NCELLS_EPS*grid->dx[dir][p->cell[dir]];
        inv_dt = MAX(inv_dt,scrh); 
      }

    /* ---------------------------------------------------------
        I2. Compute time step restriction based on the maximum
            allowed fraction of Larmor time.
       --------------------------------------------------------- */

      scrh  = DOT_PRODUCT(p->speed,B);
      scrh *= scrh; 
      scrh /= u2/(gamma*gamma) + 1.e-12;
      scrh  = Bmag2 - scrh;
      Bperp = sqrt(MAX(scrh,0.0));
      omL   = Bperp*PARTICLES_CR_E_MC/(PARTICLES_CR_LARMOR_EPS*gamma);

      inv_dt = MAX(inv_dt, omL);  /* Larmor  */

    /* -- J. Check that particle has not travelled more than one cell -- */
/*
      int indx_old[3], dindx, dindx_max;

      for (dir = 0; dir < DIMENSIONS; dir++) indx_old[dir] = p->cell[dir]; 

      Particles_LocateCell (p->coord_old, indx_old, grid);
      Particles_LocateCell (p->coord, p->cell, grid);

      dindx_max = 0;
      for (dir = 0; dir < DIMENSIONS; dir++) {
        dindx = fabs(p->cell[dir] - indx_old[dir]);
        dindx_max = MAX(dindx,dindx_max);    
      }
      if (dindx_max > Nmax){//(int)max_zones){
        print ("! Particles_Update(): particle has travelled %d zones\n",
                dindx_max);
        print ("! nparticles = %d\n",p_nparticles);
        print ("  indx_old = %d %d %d\n",indx_old[IDIR], indx_old[JDIR], indx_old[KDIR]);
        print ("  indx_new = %d %d %d\n",p->cell[IDIR], p->cell[JDIR], p->cell[KDIR]);
        QUIT_PLUTO(1);
      }
*/    
    }  /* End loop on particles */

  /* ---------------------------------------------------
     3d. Synchronize single-cycle deposit array dM and
         add contribution to cumulative array dM_tot
     --------------------------------------------------- */

    #if PARTICLES_CR_FEEDBACK == YES
    Particles_DepositBoundaryExchange (dM, 4, grid);
    for (dir = 0; dir < 4; dir++) {
      #if PARTICLES_DEPOSIT == INTEGER
      TOT_LOOP(k,j,i) {
        dM[dir][k][j][i]     /= C;
        dM_tot[dir][k][j][i] += dM[dir][k][j][i];
      }
      #else
      TOT_LOOP(k,j,i) dM_tot[dir][k][j][i] += dM[dir][k][j][i];
      #endif
    }
    #endif

  /* ----------------------------------------------------
     3e. Set boundary condition after deposition at
         x^{n+1/2} has been done
     ---------------------------------------------------- */

    Particles_Boundary(data, grid);
    Particles_BoundaryExchange(data, grid);

    #if SHOW_TIMING
    clock0 = clock();
    #endif
  }  /* End loop on sub-cycles */

  Dts->invDt_particles = inv_dt;

/* ----------------------------------------------------------
   4. Compute feedback array Fcr at t(n+1/2) needed in the
      corrector step.
   ---------------------------------------------------------- */

#if PARTICLES_CR_FEEDBACK == YES
  scrh = 1.0/dt0;
  for (dir = 0; dir < 4; dir++) {
    TOT_LOOP(k,j,i) {
      #ifdef CTU
      data->Fcr[dir][k][j][i] = dM_tot[dir][k][j][i]*scrh;
//if (isnan(data->Fcr[dir][k][j][i])){
//  print ("! Particles_CR_Update(): Fcr at end is nan\n");
//  QUIT_PLUTO(1);
//}
      #else  /* Use this for RK2 time stepping */
      data->Fcr[dir][k][j][i] = 2.0*dM_tot[dir][k][j][i]/dt0 - Fcr0[dir][k][j][i];
      #endif
    }
  }
#endif


#if SHOW_TIMING
{
  double dclock_tot;

  Dts->clock_particles = (double)(clock() - clock_beg)/CLOCKS_PER_SEC;
  dclock_tot = (double)(clock() - clock_beg)/CLOCKS_PER_SEC;
  
  print ("  Total: %f, [%8.3e per particle]\n",dclock_tot,
                                               dclock_tot/p_nparticles);
  
}
#endif
  DEBUG_FUNC_END ("CR_Update");
}

/* ********************************************************************* */
void Particles_CR_Predictor(Data *data, timeStep *Dts, double dt, Grid *grid)
/*!
 * Advance particle velocity and speed by dt using an
 * implicit-explicit 1st order scheme.
 * Compute Fcr at F(t+dt), then restore particles into their original
 * position.
 *
 *********************************************************************** */
{
  static int nparticles = -1;
  int    i, j, k, dir;
  double a[3], Fcr[3], vg[3], vp[3], B[3], E[3];
  double v2, gamma, u2, um[3], up[3];
  double h, qg, den;
  const double c2 = PARTICLES_CR_C*PARTICLES_CR_C;
  static double ***w, **Ainv, **pcoord_old, **pspeed_old;
  particleNode *curNode;
  Particle *p;

/* --------------------------------------------------------
   0. Allocate memory.
      This includes also the two arrays pcoord_old[][]
      and pspeed_old[][] used to store position and
      velocity at the current time level.
      The two arrays are re-allocated when the current
      number of particles (p_nparticles) exceeds the
      previous one (nparticles).
   -------------------------------------------------------- */

  if (w == NULL){
    w    = ArrayBox (-1, 1, -1, 1, -1, 1);
    Ainv = ARRAY_2D(3,3,double);
  }
/*
  if (nparticles < p_nparticles){
print ("Reallocating.....old = %d, new = %d\n",nparticles, p_nparticles);
    if (pcoord_old != NULL){
      FreeArray2D((void **) pcoord_old);
      FreeArray2D((void **) pspeed_old);
    }
    nparticles = p_nparticles;
    pcoord_old = ARRAY_2D(nparticles,3, double);
    pspeed_old = ARRAY_2D(nparticles,3, double);
  }
*/

/* --------------------------------------------------------
   1. Advance particles with 1st order method
   -------------------------------------------------------- */

  h = dt*PARTICLES_CR_E_MC; 
  PARTICLES_LOOP(curNode, data->PHead){

    p = &(curNode->p);

  /* --------------------------------------------
     1a. Get weights and electric field at x^n.
         Here data->Ex1, data->Ex2, data->Ex3 are
         cell-centered electric fields computed
         at t^n.
         For CT, this is done in
         CT_ComputeCenterEMF())
     -------------------------------------------- */

    Particles_GetWeights(p, p->cell, w, grid);

    E[IDIR] = Particles_Interpolate(data->Ecr[IDIR], w, p->cell); 
    E[JDIR] = Particles_Interpolate(data->Ecr[JDIR], w, p->cell); 
    E[KDIR] = Particles_Interpolate(data->Ecr[KDIR], w, p->cell);
/*
//    E[IDIR] = Particles_Interpolate(data->Ex1, w, p->cell); 
//    E[JDIR] = Particles_Interpolate(data->Ex2, w, p->cell); 
//    E[KDIR] = Particles_Interpolate(data->Ex3, w, p->cell);
*/
  /* --------------------------------------------
     1b. Kick step: evolve four-velocity using
         electric field only:

         um^* = u^n + h/2*(cE)^n
     -------------------------------------------- */
   
    v2    = DOT_PRODUCT(p->speed, p->speed);
    gamma = 1.0/sqrt(1.0 - v2/c2);
    for (dir = 0; dir < 3; dir++) {
      um[dir] = gamma*p->speed[dir] + h*E[dir];
    }
    u2    = DOT_PRODUCT(um, um);
    gamma = sqrt(1.0 + u2/c2);

  /* -- 1c. Save initial coordinate and position

            p->coord_old = x^n
            p->speed_old = v^n

            and update spatial position: x^* = x^n + dt*vm^ -- */
  
    for (dir = 0; dir < 3; dir++){
      p->coord_old[dir] = p->coord[dir];
      p->speed_old[dir] = p->speed[dir];
      p->coord[dir]     += dt*um[dir]/gamma;
    }

  /* -- 1d. Get weights and magnetic field at x^* -- */

    Particles_GetWeights(p, p->cell, w, grid);  

    B[IDIR] = Particles_Interpolate(data->Vc[BX1], w, p->cell); 
    B[JDIR] = Particles_Interpolate(data->Vc[BX2], w, p->cell); 
    B[KDIR] = Particles_Interpolate(data->Vc[BX3], w, p->cell);

  /* -- 1e. Rotate: add magnetic field contribution -- */

    B[IDIR] *= h/gamma;
    B[JDIR] *= h/gamma;
    B[KDIR] *= h/gamma;

    den = 1.0 + B[IDIR]*B[IDIR] + B[JDIR]*B[JDIR] + B[KDIR]*B[KDIR];

    Ainv[IDIR][IDIR] =      1.0 + B[IDIR]*B[IDIR];
    Ainv[IDIR][JDIR] =  B[KDIR] + B[IDIR]*B[JDIR];
    Ainv[IDIR][KDIR] = -B[JDIR] + B[IDIR]*B[KDIR];

    Ainv[JDIR][IDIR] = -B[KDIR] + B[JDIR]*B[IDIR];
    Ainv[JDIR][JDIR] =      1.0 + B[JDIR]*B[JDIR];
    Ainv[JDIR][KDIR] =  B[IDIR] + B[JDIR]*B[KDIR];

    Ainv[KDIR][IDIR] =  B[JDIR] + B[KDIR]*B[IDIR];
    Ainv[KDIR][JDIR] = -B[IDIR] + B[KDIR]*B[JDIR];
    Ainv[KDIR][KDIR] =      1.0 + B[KDIR]*B[KDIR];

  /* -- 1f. Advance coordinate and velocity by dt: -- */

    for (dir = 0; dir < 3; dir++){
      up[dir]  = DOT_PRODUCT(Ainv[dir],um);
      up[dir] /= den;
    }
    u2    = DOT_PRODUCT(up,up);
    gamma = sqrt(1.0 + u2/c2);

    for (dir = 0; dir < 3; dir++){
      p->speed[dir] = up[dir]/gamma;
      p->coord[dir]  = p->coord_old[dir] + dt*p->speed_old[dir];
    }
  }
  
/* --------------------------------------------------------
   2. Compute Fcr(x^{n+1/2}, v^{n+1/2})
   -------------------------------------------------------- */
  
  Particles_CR_ComputeForce(data->Vc, data, grid);

/* --------------------------------------------------------
   3. Restore coordinates and velocities
   -------------------------------------------------------- */

  PARTICLES_LOOP(curNode, data->PHead){
    p = &(curNode->p);
    for(dir = 0; dir < 3; dir++) {
      p->speed[dir] = p->speed_old[dir];
      p->coord[dir] = p->coord_old[dir];
    }
  }

}
#endif  /* (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_TEST == NO) */
