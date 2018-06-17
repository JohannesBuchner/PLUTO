/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Initialize particle for gyration test.
 
 \authors A. Mignone (mignone@ph.unito.it)\n
          B. Vaidya (bvaidya@unito.it)\n
 
 \date    March 20, 2016
 
 \b References: \n
   - "A PARTICLE MODULE FOR THE PLUTO CODE: I - AN IMPLEMENTATION OF THE
      MHD-PIC EQUATIONS", Mignone etal.ApJS (2018)  [ Sec. 4.1 ]

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
  double xbeg[3], xend[3], vbeg[3], vend[3];
  double gamma, alpha, u1, v1;
  double vx1, vy1, vz1, vg;
  double c2 = PARTICLES_CR_C*PARTICLES_CR_C;
  Particle p;

/* --------------------------------------------------------------
    Global initialization
   -------------------------------------------------------------- */

  vg = g_inputParam[VFLUID_X];
  if (np_glob > 0){

    for (np = 0; np < np_glob; np++){

      p.coord[0] = 0.0;
      p.coord[1] = 0.0;
      p.coord[2] = 0.0;

      p.color = 0.0;

    /* -- Compute particle velocity in fluid co-moving frame -- */

      u1    = g_inputParam[PARTICLE_4VEL];  /* Particle 4-vel in co-moving frame */
      gamma = sqrt(1.0 + u1*u1/c2);
      v1    = u1/gamma;

      alpha = CONST_PI*g_inputParam[ALPHA]/180.0;
      vx1   = v1*cos(alpha);
      vy1   = v1*sin(alpha);
      vz1   = 0.0;

    /* -- Compute particle velocity in Lab frame -- */

      gamma = 1.0/sqrt(1.0 - vg*vg/c2);
      p.speed[IDIR] = (vx1 + vg)/(1.0 + vx1*vg/c2);
      p.speed[JDIR] = vy1/(1.0 + vx1*vg/c2)/gamma;
      p.speed[KDIR] = vz1/(1.0 + vx1*vg/c2)/gamma;

      p.rho = 1.0;
     
      Particles_Insert (&p, d, PARTICLES_CREATE, grid);
    }
  }

  Particles_SetID(d->PHead);
}

/* ********************************************************************* */
void Particles_Inject(Data *data, Grid *grid)
/*!
 *  Sets user-defined boundary conditions on particles.
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
