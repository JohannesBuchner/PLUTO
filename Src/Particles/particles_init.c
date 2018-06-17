/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Initialize particle distrbution, set B.C. and injection.
 
 This file contains routines to initialize particles on the grid,
 assign user-defined boundary conditions and inject particles.
 Particles attributes that can be set are: position, velocity, color.
 In case of evolution with spectra, the initial spectral profile is also
 prescribed here.
 
 \authors A. Mignone (mignone@ph.unito.it)\n
          B. Vaidya (bvaidya@unito.it)\n
 
 \date    May 13, 2018
 */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Particles_Init(Data *d, Grid *grid)
/*!
 *  Sets initial conditions on particles.
 *
 *  \param [in]    d       Pointer to the PLUTO data structure.
 *  \param [in]    grid    Pointer to the PLUTO grid structure.
 *
 *********************************************************************** */
{
  int i,j,k, np, dir;
  int np_glob = RuntimeGet()->Nparticles_glob;
  int np_cell = RuntimeGet()->Nparticles_cell;
  static int first_call = 1;
  double xbeg[3], xend[3];
  Particle p;


  if (first_call){
    RandomSeed(time(NULL),0);
    first_call = 0;
  }

/* --------------------------------------------------------------
   1. Global initialization
   -------------------------------------------------------------- */

  if (np_glob > 0){

    for (dir = 0; dir < 3; dir++){
      xbeg[dir] = grid->xbeg_glob[dir];
      xend[dir] = grid->xend_glob[dir];
    }

    for (np = 0; np < np_glob; np++){
      Particles_LoadUniform(np, np_glob, xbeg, xend, p.coord);
      #if (PARTICLES_TYPE == LAGRANGIAN) && (PARTICLES_LP_SPECTRA == YES)
      Particles_LP_InitSpectra(&p);  
      #endif
      #if PARTICLES_TYPE == COSMIC_RAYS
      p.speed[IDIR] = 0.0;
      p.speed[JDIR] = 0.0;
      p.speed[KDIR] = 0.0;
      p.rho   = 1.e-3;
      #endif
      p.color = 0.0; 
      Particles_Insert (&p, d, PARTICLES_CREATE, grid);
    }
  }

/* ------------------------------------------------------------------
   2. Cell by cell initialization.
      Note: You may use Particles_LoadRandom() to initialize
            velocity components but not spatial coordinates.
   ------------------------------------------------------------------ */

  if (np_cell > 0){

    DOM_LOOP(k,j,i){
      xbeg[IDIR] = grid->xl[IDIR][i]; xend[IDIR] = grid->xr[IDIR][i];
      xbeg[JDIR] = grid->xl[JDIR][j]; xend[JDIR] = grid->xr[JDIR][j];
      xbeg[KDIR] = grid->xl[KDIR][k]; xend[KDIR] = grid->xr[KDIR][k];

      for (np = 0; np < np_cell; np++){

      /* -- Spatial distribution -- */

        Particles_LoadUniform(np, np_cell, xbeg, xend, p.coord);

        #if (PARTICLES_TYPE == LAGRANGIAN) && (PARTICLES_LP_SPECTRA == YES)
        Particles_LP_InitSpectra(&p);
        #endif

     /* -- Velocity distribution -- */

        #if (PARTICLES_TYPE == COSMIC_RAYS)
        p.speed[IDIR] = RandomNumber(-1,1);
        p.speed[JDIR] = RandomNumber(-1,1);
        p.speed[KDIR] = RandomNumber(-1,1);
        p.rho         = 1.e-3/np_cell;
        #endif
        p.color = 0.0;
        Particles_Insert (&p, d, PARTICLES_CREATE, grid);
      }
    }
  }  
  Particles_SetID(d->PHead);
}

#if PARTICLES_LP_SPECTRA == YES
/* ********************************************************************** */
void Particles_LP_InitSpectra(Particle* pl)
/*!
 *  Initialize spectra for each particle (only for LAGRANGIAN).
 *  Specify here the initial distribution of N(E) with E for each particle
 *
 *  \param [in]      pl      Pointer to the Particle structure.
 * 
 ********************************************************************** */
{
  int i;
  double Emin, Emax, DeltaE, N_0, alpha1, alpha2, Eb, N1, N2;
  double lnEmin, lnEmax, dlnE, scrh;
  Emin = 1.0e-2;
  Emax = 1.0e4;
  lnEmin = log10(Emin);
  lnEmax = log10(Emax);
  dlnE = (lnEmax - lnEmin)/((double)PARTICLES_LP_NEBINS-1);
    
  pl->nmicro = 0.001; /* The number density of micro-particles */
  for (i=0; i< PARTICLES_LP_NEBINS; i++){
    scrh  = lnEmin + i*dlnE;
    pl->eng[i] = pow(10.0, scrh);
          
    /* /\*Single Power Law*\/ */
    alpha1 = 3.0;
    N_0 = (pl->nmicro)*(1.0-alpha1)/(pow(Emax,1.0-alpha1)-pow(Emin,1.0-alpha1));
    pl->chi[i] = N_0 * pow(pl->eng[i],-alpha1);
    pl->cmp_ratio = 1.0;
    pl->shkflag = 0;
    pl->shk_vL[RHO] = -1.0;
    pl->shk_vR[RHO] = -1.0;
    pl->ca = 0.0;
    pl->cr = 0.0;

    /* Relativistic Maxwellian Distribution*/
    /* double mu = 2.0; // mu = m_0c^2/kT_e */
    /* double k2_mu = BesselKn(2,mu); */
    /* pl->chi[i] = (pl->nmicro)*(mu/k2_mu)*(pl->eng[i] + 1.0)*sqrt(pl->eng[i]*(pl->eng[i] + 2.0))*exp(-mu*(pl->eng[i] + 1.0)); */
    /* pl->cmp_ratio = 1.0; */
    /* pl->shkflag = 0;*/
    /* pl->rho_min = -1.0;*/
    /*pl->rho_max = -1.0;*/
    /*pl->ca = 0.0; */
    /*pl->cr = 0.0; */
  }
}
#endif

/* ********************************************************************* */
void Particles_Inject(Data *data, Grid *grid)
/*!
 *  Inject particles as you wish.
 *
 *  \param [in]  data    Pointer to the PLUTO data structure.
 *  \param [in]  grid    Pointer to the PLUTO grid structure.
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Particles_UserDefBoundary(Data *d, int side, Grid *grid)
/*
 *
 *********************************************************************** */
{
  int    dir;
  double xbeg[3], xend[3];
  particleNode *curr = d->PHead, *next;
  Particle *p;
  
  for (dir = 0; dir < 3; dir++) {
    xbeg[dir] = grid->xbeg_glob[dir];
    xend[dir] = grid->xend_glob[dir];
  }

  if (side == X1_BEG){
    dir = IDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X1_END){
    dir = IDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X2_BEG){
    dir = JDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X2_END){
    dir = JDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X3_BEG){
    dir = KDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X3_END){
    dir = KDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }
}
