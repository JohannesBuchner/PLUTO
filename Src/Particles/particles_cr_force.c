/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute CR-related quantities (like force) from the particles to the grid.

  This routine is used to compute 3D arrays of cell-centered quantities:

  - data->Fcr (Lorentz force felt by the fluid)
  - data->Jcr (CR current density)
  - data->qcr (CR charge density)
  - data->Ecr (Total electric field)

  Current and charge density are first deposited from particles to the grid.
  With these, we compute Fcr and Ecr.

  \authors A. Mignone (mignone@ph.unito.it)\n

  \b References
    - "MAGNETOHYDRODYNAMIC-PARTICLE-IN-CELL METHOD FOR COUPLING COSMIC RAYS 
       WITH A THERMAL PLASMA: APPLICATION TO NON-RELATIVISTIC SHOCKS"\n
       Bai et al., ApJ (2015) 809, 55

  \date   March 31, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

void Particles_CR_Current (Particle *p, double *qd);

#if PARTICLES_TYPE == COSMIC_RAYS 
/* ********************************************************************* */
void Particles_CR_ComputeForce(Data_Arr V, const Data *data, Grid *grid)
/*!
 *
 * \param [in]         V   3D array of primitive variables, V[nv][k][j][i]
 * \param [in,out]  data   PLUTO Data structure, containing d->Fcr[]
 * \param [in]      grid   array of grid structures
 *  
 *********************************************************************** */
{
  int    i,j,k,dir;
  double R1, qg, inv_qg, emf0[3], Jcr_X_B[3];
  double vx1, vx2, vx3;
  double Bx1, Bx2, Bx3;
  double ***Jmu[4];   /* Pointer shortcut to density and currents */

#if PARTICLES_CR_FEEDBACK == NO
  return;
#endif

//print (">>[Particles_CR_Force], step = %d\n", g_stepNumber);

/* -------------------------------------------------------------
   1. Compute charge and current density from particles to grid 
   ------------------------------------------------------------- */

  Jmu[0] = data->qcr;     
  Jmu[1] = data->Jcr[0];
  Jmu[2] = data->Jcr[1];
  Jmu[3] = data->Jcr[2];

  Particles_Deposit (data->PHead, Particles_CR_Current, Jmu, 4, grid);

/* --------------------------------------------------------------
   2. Loop over domain to compute total electric field and force
   -------------------------------------------------------------- */
  
  TOT_LOOP(k,j,i){

    vx1 = V[VX1][k][j][i]; Bx1 = V[BX1][k][j][i];
    vx2 = V[VX2][k][j][i]; Bx2 = V[BX2][k][j][i];
    vx3 = V[VX3][k][j][i]; Bx3 = V[BX3][k][j][i];

    emf0[IDIR] = -(vx2*Bx3 - vx3*Bx2);
    emf0[JDIR] = -(vx3*Bx1 - vx1*Bx3);
    emf0[KDIR] = -(vx1*Bx2 - vx2*Bx1);    

  /* -- Compute the J X B term -- */
    
    Jcr_X_B[IDIR] = data->Jcr[JDIR][k][j][i]*Bx3 - data->Jcr[KDIR][k][j][i]*Bx2;
    Jcr_X_B[JDIR] = data->Jcr[KDIR][k][j][i]*Bx1 - data->Jcr[IDIR][k][j][i]*Bx3;
    Jcr_X_B[KDIR] = data->Jcr[IDIR][k][j][i]*Bx2 - data->Jcr[JDIR][k][j][i]*Bx1;

  /* -- Compute Force on the grid -- */
    
    qg = PARTICLES_CR_E_MC_GAS*V[RHO][k][j][i];
    R1 = qg/(qg + data->qcr[k][j][i]);            /* This is (1 - R) */

    inv_qg = 1.0/qg;
    for (dir = 0; dir < 3; dir++) {   
      data->Fcr[dir][k][j][i] = R1*(data->qcr[k][j][i]*emf0[dir] + Jcr_X_B[dir]);
      data->Ecr[dir][k][j][i] = emf0[dir] - data->Fcr[dir][k][j][i]*inv_qg;
    }

  /* -- Compute source term to energy equation, vg*Fcr -- */

    data->Fcr[3][k][j][i] =   vx1*data->Fcr[IDIR][k][j][i]
                            + vx2*data->Fcr[JDIR][k][j][i]
                            + vx3*data->Fcr[KDIR][k][j][i];
  }

//print ("<<[Particles_CR_Force]\n");
}

/* ********************************************************************* */
void Particles_CR_ComputeCurrent(const Data *d, Grid *grid)
/*!
 * Compute Fcr (Lorentz force felt by particles and fluid) at cell centers.
 * Force is compute by depositing charge and current from individual
 * particles to the grid.
 *
 * \param [in]      V   3D array of primitive variables, V[nv][k][j][i]
 * \param [in,out]  d   PLUTO Data structure, containing d->Fcr[]
 * \param [in]   grid   array of grid structures
 *  
 *********************************************************************** */
{
  double ***Jmu[4];   /* Pointer shortcut to density and currents */

#if PARTICLES_CR_FEEDBACK == NO
  return;
#endif

//print (">>[Particles_CR_Force], step = %d\n", g_stepNumber);

/* -------------------------------------------------------------
   1. Compute charge and current density from particles to grid 
   ------------------------------------------------------------- */

  Jmu[0] = d->qcr;     
  Jmu[1] = d->Jcr[0];
  Jmu[2] = d->Jcr[1];
  Jmu[3] = d->Jcr[2];

  Particles_Deposit (d->PHead, Particles_CR_Current, Jmu, 4, grid);

}
 
/* ************************************************************* */
void Particles_CR_Current (Particle *p, double *qd)
/*!
 *  Compute the quantities to be deposited to the grid. 
 *************************************************************** */
{
  const double q = PARTICLES_CR_E_MC * p->rho;

  qd[0] = q;
  qd[1] = q*p->speed[IDIR];
  qd[2] = q*p->speed[JDIR];
  qd[3] = q*p->speed[KDIR];
}

#endif  /* PARTICLES_TYPE == COSMIC_RAYS */
