/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Initialize particles for the Bell instability test-
 
 \authors A. Mignone (mignone@ph.unito.it)\n
          B. Vaidya (bvaidya@unito.it)\n
 
 \date    March 20, 2018
 
 \b References: \n
    - [MVBM18] "A PARTICLE MODULE FOR THE PLUTO CODE: I - AN IMPLEMENTATION OF
                THE MHD-PIC EQUATIONS", Mignone etal.ApJS (2018)  [ Sec. 4.4 ]
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

void VectorRotate(double *v, int s);
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
  Particle p;

/* --------------------------------------------------------------
    Global initialization
   -------------------------------------------------------------- */

  if (np_glob > 0){
    print ("! Use cell-by-cell initialization\n");
    QUIT_PLUTO(1);
  }

/* --------------------------------------------------------------
    Cell by cell initialization
   -------------------------------------------------------------- */

  if (np_cell > 0){
    static int first_call = 1;
    double rho_p = 2.e6*g_inputParam[EPSILON]/np_cell;
    double uCR = 1.0/g_inputParam[EPSILON];  /* = vA/epsilon */
    double Jcr = PARTICLES_CR_E_MC*rho_p*uCR;   /* = Jcr/c      */ 

    DOM_LOOP(k,j,i){

      xbeg[IDIR] = grid->xl[IDIR][i]; xend[IDIR] = grid->xr[IDIR][i];
      xbeg[JDIR] = grid->xl[JDIR][j]; xend[JDIR] = grid->xr[JDIR][j];
      xbeg[KDIR] = grid->xl[KDIR][k]; xend[KDIR] = grid->xr[KDIR][k];

  /* -- Loop on particles -- */
  
      for (np = 0; np < np_cell; np++){
        Particles_LoadUniform(np, np_cell, xbeg, xend, p.coord);
        p.speed[IDIR] = uCR;
        p.speed[JDIR] = 0.0;
        p.speed[KDIR] = 0.0;
        p.rho         = rho_p;
        p.color       = 0.0;
        VectorRotate(p.speed, 1);

#if PARTICLES_TYPE == COSMIC_RAYS || PARTICLES_TYPE == DUST
//        Particles_LoadRandom(vbeg, vend, Particles_VelocityDistrib, p.speed);
#endif
        Particles_Insert (&p, d, PARTICLES_CREATE, grid);
      }
    }

    if (first_call){
      double lambda = 4.0*CONST_PI/Jcr;
      double R;

      R = PARTICLES_CR_E_MC*rho_p;
      R = R/(R + PARTICLES_CR_E_MC_GAS);
      print ("> Particles_Init():\n");
      print ("  ------------------------------------ \n");
      print ("  normal               = [%f, %f, %f]\n", p.speed[IDIR]/uCR,
                                                        p.speed[JDIR]/uCR,
                                                        p.speed[KDIR]/uCR);
      print ("  gyration radius      = %8.3e\n", uCR/PARTICLES_CR_E_MC);
      print ("  most unstable lambda = %f  (k = %f)\n", lambda, 2.0*CONST_PI/lambda);
      print ("  R                    = %8.3e\n", R);
      print ("  Lambda               = %8.3e\n", R/g_inputParam[EPSILON]);
      print ("  ------------------------------------ \n");
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
