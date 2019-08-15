#include "pluto.h"

static void MyDeposit (Particle *p, double *qd);

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k;  
  double ***Ep = GetUserVar("Ep");
  
  Particles_Deposit(d->PHead, MyDeposit, &Ep, 1, grid);
}

#ifdef PARTICLES
/* ************************************************************* */
void MyDeposit (Particle *p, double *qd)
/*
 *************************************************************** */
{
  qd[0] = (  p->speed[0]*p->speed[0]
           + p->speed[1]*p->speed[1]
           + p->speed[2]*p->speed[2] );
  qd[0] = p->speed[2]; // sqrt(qd[0]); 
}
#endif

/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

}





