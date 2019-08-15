/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Update CR particles (test particles only)
 
  This function provides a simple implementation of the Boris pusher in
  which electromagnetic fields are supplied externally through the
  user-defined function Particles_CR_EMFields().
  Fluid quanatities may not be defined at all.
  
  \authors A. Mignone (mignone@ph.unito.it)\n

  \date   Dec 28, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_TEST == YES)
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
  static double ***w;

  particleNode *CurNode;
  Particle *p;

#if SHOW_TIMING
  clock_t clock_beg = clock(), clock0;
#endif

/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */
  
  if (w == NULL)  w = ArrayBox (-1, 1, -1, 1, -1, 1);

/* --------------------------------------------------------
   1. Initialize sub-cycling
   -------------------------------------------------------- */

  inv_dt  = 1.e-18;
  dt0     = dt;  /* Save initial dt */
  dt      = dt/(double)Dts->Nsub_particles;
  dt_half = 0.5*dt;
  h_half  = 0.5*dt*PARTICLES_CR_E_MC; 

/* --------------------------------------------------------
   2. Start sub-cycling
   -------------------------------------------------------- */

  for (kcycle = 0; kcycle < Dts->Nsub_particles; kcycle++){

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

    /* -- B. Compute electromagnetic fields at x^{n+1/2} -- */
  
      Particles_GetWeights(p, p->cell, w, grid); /* Need cell indices to
                                                    compute max distance */
      Particles_CR_EMFields(p->coord, E, B);
      Bmag2 = DOT_PRODUCT(B,B);
  
    /* -- C. Boris pusher (kick + rotate + kick) -- */
     
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
    
    /* -- D. Update velocity by another half step -- */

      u[IDIR] = up[IDIR] + h_half*E[IDIR];
      u[JDIR] = up[JDIR] + h_half*E[JDIR];
      u[KDIR] = up[KDIR] + h_half*E[KDIR];
    
      u2    = DOT_PRODUCT(u,u);
      gamma = sqrt(1.0 + u2/c2);
      scrh  = 1.0/gamma;

      p->speed[IDIR] = u[IDIR]*scrh;
      p->speed[JDIR] = u[JDIR]*scrh;
      p->speed[KDIR] = u[KDIR]*scrh;

    /* -- E. Update spatial coordinate to obtain x^{n+1} -- */

      p->coord[IDIR] += dt_half*p->speed[IDIR];
      p->coord[JDIR] += dt_half*p->speed[JDIR];
      p->coord[KDIR] += dt_half*p->speed[KDIR];

    /* ---------------------------------------------------------
        F1. Compute time step restriction based on the maximum
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
        F2. Compute time step restriction based on the maximum
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

  /* ----------------------------------------------------
     Set boundary conditions.
     ---------------------------------------------------- */

    Particles_Boundary(data, grid);
    Particles_BoundaryExchange(data, grid);

    #if SHOW_TIMING
    clock0 = clock();
    #endif
  }  /* End loop on sub-cycles */

  Dts->invDt_particles = inv_dt;


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

#endif  /* (PARTICLES_TYPE == COSMIC_RAYS) && (PARTICLES_TEST == YES) */
