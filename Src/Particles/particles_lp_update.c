/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Routines for the temporal update of particles position and/or velocity.
 
 \authors  A. Mignone (mignone@ph.unito.it)\n
           B. Vaidya (bvaidya@unito.it)\n
 
 \date   Nov 21, 2017
 */
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if PARTICLES_TYPE == LAGRANGIAN
/* ********************************************************************* */
void Particles_LP_Predictor(Data *d, timeStep *Dts, double dt, Grid *grid)
/*!
 * Predictor step for time g_dt before fluid update.
 * \f[
 *    xstar = xcurr + dt*vparticle
 * \f]
 * Where vparticle is velocity of fluid at xcurr. 
 *
 * \param [in]   d      Pointer to PLUTO data structure.
 * \param [in]   dt     Time step dt.
 * \param [in]   grid   Pointer to the PLUTO grid structure.
 *********************************************************************** */
{  
  int dir;
  particleNode *curNode;
  Particle *p;
  static double ***w;

  Particles_Inject (d, grid);

  p_intStage = 1;
    
  if (w == NULL) w = ArrayBox (-1, 1, -1, 1, -1, 1);

  PARTICLES_LOOP(curNode, d->PHead){
    p = &(curNode->p);

    #if PARTICLES_LP_SPECTRA == YES
    int e;
    Particles_LP_FixValue(p, d, grid); /* Required for Particles that are Injected */
    p->pressure_old = p->pressure;
    p->ca_old       = p->ca;
    p->prev_shkflag = p->shkflag;
    p->density_old  = p->density;
    p->lorG_old     = p->lorG;
    
    for(e=0; e < PARTICLES_LP_NEBINS; e++){
      p->chi_old[e] = p->chi[e];
      p->eng_old[e]    = p->eng[e];
    }
    #if PHYSICS == MHD || PHYSICS == RMHD
    p->cr_old = p->cr;
    #endif

    #else
    Particles_GetWeights(p, p->cell, w, grid);  
    D_EXPAND(p->speed[IDIR] = Particles_Interpolate(d->Vc[VX1], w, p->cell); ,
             p->speed[JDIR] = Particles_Interpolate(d->Vc[VX2], w, p->cell); ,
             p->speed[KDIR] = Particles_Interpolate(d->Vc[VX3], w, p->cell);)
    #endif

    for (dir = 0; dir < DIMENSIONS; dir++) p->coord_old[dir] = p->coord[dir];

    #if GEOMETRY == CARTESIAN
    for (dir = 0; dir < DIMENSIONS; dir++) {
      p->coord[dir] = p->coord[dir] + dt*p->speed[dir];
    }
/*
print ("[Pred]: p= (%d, %d) p->coord_old = [%f, %f]; p->speed = [%f, %f]; p->coor = [%f, %f]\n",
        p->id, p->birth_rank,
        p->coord_old[IDIR],p->coord_old[JDIR],
        p->speed[IDIR],   p->speed[JDIR],
        p->coord[IDIR],    p->coord[JDIR]);
*/
    #elif GEOMETRY == POLAR 
    p->coord[IDIR] = p->coord[IDIR] + dt*p->speed[IDIR];
    p->coord[JDIR] = p->coord[JDIR] + dt*(p->speed[JDIR]/p->coord_old[IDIR]);
    p->coord[KDIR] = (DIMENSIONS == 3 ? p->coord[KDIR]+dt*p->speed[KDIR]:0.0);
    #elif GEOMETRY == SPHERICAL
    p->coord[IDIR] = p->coord[IDIR] + dt*p->speed[IDIR];
    p->coord[JDIR] = p->coord[JDIR] + dt*(p->speed[JDIR]/p->coord_old[IDIR]);
    p->coord[KDIR] = (DIMENSIONS == 3 ? p->coord[KDIR]+dt*(p->speed[KDIR]/(p->coord_old[IDIR]*sin(p->coord_old[JDIR]))):0.0);
    #elif GEOMETRY == CYLINDRICAL
    p->coord[IDIR] = p->coord[IDIR] + dt*p->speed[IDIR];
    p->coord[JDIR] = p->coord[JDIR] + dt*p->speed[JDIR];
    p->coord[KDIR] = (DIMENSIONS == 3 ? p->coord[KDIR]+dt*(p->speed[KDIR]/p->coord_old[IDIR]):0.0);
    #endif
           
  }  /* End loop on particles */
  
}

/* ********************************************************************* */
void Particles_LP_Corrector(Data *d, timeStep *Dts, double dt, Grid *grid)
/*!
 * Corrector step after fluid update
 * Locate Particle at xcurr and interpolate velocity at that position.
 * \f[
 *     xnew = xcurr + 0.5*dt*(vinterp(xstar) + vparticle)
 * \f]
 *
 * \param [in]   d      Pointer to PLUTO data structure.
 * \param [in]   dt     Time step dt.
 * \param [in]   grid   Pointer to the PLUTO grid structure.
 *********************************************************************** */
{
  int          dir;
  particleNode *curNode = d->PHead;
  Particle *p;
  double dl, scrh1, scrh2, scrh3;
  static double ***w;

  p_intStage = 2;
  
  if (w == NULL)  w = ArrayBox (-1, 1, -1, 1, -1, 1);

#if SHOW_TIMING
  clock_t clock_beg = clock();
#endif

  while (curNode != NULL){
    p = &(curNode->p);

    Particles_GetWeights(p, p->cell, w, grid);

    D_EXPAND(p->speed[IDIR] = Particles_Interpolate(d->Vc[VX1], w, p->cell); ,
             p->speed[JDIR] = Particles_Interpolate(d->Vc[VX2], w, p->cell); ,
             p->speed[KDIR] = Particles_Interpolate(d->Vc[VX3], w, p->cell);)

    #if GEOMETRY == CARTESIAN
    for (dir = 0; dir < DIMENSIONS; dir++) {
      #if TIME_STEPPING == RK2
      p->coord[dir] = 0.5*(p->coord[dir] + p->coord_old[dir] + dt*p->speed[dir]);
      #elif (TIME_STEPPING == HANCOCK) || (TIME_STEPPING == CHARACTERISTIC_TRACING)
      p->coord[dir] = p->coord_old[dir] + dt*p->speed[dir];
      #else
        #error Integrator not defined in Particles_Corrector()
      #endif
    }

    #elif GEOMETRY == POLAR
    scrh1 = p->coord[IDIR];
    p->coord[IDIR] = 0.5*(p->coord[IDIR] + p->coord_old[IDIR] + dt*p->speed[IDIR]);
    p->coord[JDIR] = 0.5*(p->coord[JDIR] + p->coord_old[JDIR] + dt*(p->speed[JDIR]/scrh1));
    scrh2 = 0.5*(p->coord[KDIR] + p->coord_old[KDIR] + dt*p->speed[KDIR]);
    p->coord[KDIR] = (DIMENSIONS == 3 ? scrh2:0.0);
    #elif GEOMETRY == SPHERICAL
    scrh1 = p->coord[IDIR];
    scrh2 = p->coord[JDIR];
    p->coord[IDIR] = 0.5*(p->coord[IDIR] + p->coord_old[IDIR] + dt*p->speed[IDIR]);
    p->coord[JDIR] = 0.5*(p->coord[JDIR] + p->coord_old[JDIR] + dt*p->speed[JDIR]/scrh1);
    scrh3 = 0.5*(p->coord[KDIR] + p->coord_old[KDIR] + dt*(p->speed[KDIR]/(scrh1*sin(scrh2))));
    p->coord[KDIR] = (DIMENSIONS == 3 ? scrh3:0.0);
    #elif GEOMETRY == CYLINDRICAL
    scrh1 = p->coord[IDIR];
    p->coord[IDIR] = 0.5*(p->coord[IDIR] + p->coord_old[IDIR] + dt*p->speed[IDIR]);
    p->coord[JDIR] = 0.5*(p->coord[JDIR] + p->coord_old[JDIR] + dt*p->speed[JDIR]);
    scrh2 = 0.5*(p->coord[KDIR] + p->coord_old[KDIR] + dt*(p->speed[KDIR]/scrh1));
    p->coord[KDIR] = (DIMENSIONS == 3 ? scrh2:0.0);
    #endif

    curNode = curNode->next;
  }
#if SHOW_TIMING
  Dts->clock_particles = (double)(clock() - clock_beg)/CLOCKS_PER_SEC;
#endif 
  Particles_Boundary(d, grid);
  Particles_BoundaryExchange(d, grid);

#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}
#endif /* PARTICLES_TYPE == LAGRANGIAN */

