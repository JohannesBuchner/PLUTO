#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  static int first_call = 1;
  double rnd;
  double x, y, z = 0.0, r2; 
  double kp, mu;
  double dvx, dvy, dBx, dBy, arg;
  
  g_gamma = 5./3.;

  if (first_call == 1){
    srand(time(NULL) + prank);
    first_call = 0;
  }
  rnd = (double)(rand())/((double)RAND_MAX + 1.0);
  
  D_EXPAND(x = x1;, y = x2;, z = x3;)

  kp  = g_inputParam[KPAR];
  mu  = g_inputParam[MUPAR];
  r2  = x*x + y*y + z*z;
  arg = (1.0 - r2);

  dvx = -y*kp*exp(0.5*arg);
  dvy =  x*kp*exp(0.5*arg);

  dBx = -y*mu*exp(0.5*arg);
  dBy =  x*mu*exp(0.5*arg);
  
  us[RHO] = 1.0; // mu*mu/(kp*kp);

  us[VX1] = dvx + g_inputParam[PERT_AMPL]*(rnd-0.5)*exp(-x1*x1/4.0);
  us[VX2] = dvy + 0.5*g_inputParam[MACH]*tanh(x1);
  us[VX3] = 0.0;
  #if DIMENSIONS == 3
   us[VX3] = 0.1*sin(2.0*CONST_PI*x3/10.0);
  #endif

  us[PRS] = 1.0/g_gamma + exp(arg)*(mu*mu*(1.0 - r2 + z*z) - us[RHO]*kp*kp)*0.5;

  #if PHYSICS == MHD
   us[BX1] = dBx;
   us[BX2] = dBy;  
   us[BX3] = 0.0;

   #ifdef STAGGERED_MHD
    us[AX1] = us[AX2] = 0.0;
    us[AX3] = mu*exp(0.5*arg);
   #endif
  #endif

}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{
  int    i,j,k;
  double rhoe, pm, kin, my, E, BxBy, vol;
  double pm0, kinx, kiny;
  FILE *fp;

/* --------------------------------------------------------
              compute volume averages 
   -------------------------------------------------------- */

  rhoe = pm = kin = my = E = BxBy = 0.0;
  DOM_LOOP(k,j,i){
    #if PHYSICS == MHD
     pm0  = 0.5*(  d->Vc[BX1][k][j][i]*d->Vc[BX1][k][j][i]
                 + d->Vc[BX2][k][j][i]*d->Vc[BX2][k][j][i]);
    #endif
    kinx = 0.5*d->Vc[RHO][k][j][i]*d->Vc[VX1][k][j][i]*d->Vc[VX1][k][j][i];
    kiny = 0.5*d->Vc[RHO][k][j][i]*d->Vc[VX2][k][j][i]*d->Vc[VX2][k][j][i];

    pm  += pm0;
    kin += kinx;

    rhoe += d->Vc[PRS][k][j][i]/(g_gamma - 1.0);
    my   += d->Vc[RHO][k][j][i]*d->Vc[VX2][k][j][i];
    E    += d->Vc[PRS][k][j][i]/(g_gamma - 1.0) + kinx + kiny + pm0;
    #if PHYSICS == MHD
     BxBy += d->Vc[BX1][k][j][i]*d->Vc[BX2][k][j][i];
    #endif
  }

  vol  = grid[IDIR].dx[IBEG]*grid[JDIR].dx[JBEG];
  vol /= (g_domEnd[IDIR] - g_domBeg[IDIR])*(g_domEnd[JDIR] - g_domBeg[JDIR]);

  rhoe *= vol;
  pm   *= vol;
  kin  *= vol;
  my   *= vol;
  E    *= vol;
  BxBy *= vol;

  #ifdef PARALLEL
   MPI_Allreduce (&rhoe, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   rhoe = vol;

   MPI_Allreduce (&pm, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   pm = vol;

   MPI_Allreduce (&kin, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   kin = vol;

   MPI_Allreduce (&my, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   my = vol;

   MPI_Allreduce (&E, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   E = vol;

   MPI_Allreduce (&BxBy, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   BxBy = vol;

   MPI_Barrier (MPI_COMM_WORLD);
  #endif

/* -----------------------------------------------------
         processor #0 does the actual writing 
   ----------------------------------------------------- */

  if (prank == 0){
    static int first_call = 1;
    if (first_call){
      fp = fopen("averages.dat","w");
      fprintf (fp,"#   t\t\t  dt\t      <rhoe>\t   <B^2>/2\t  <rho*vx^2/2>\t  <rho*vy>\t  <E>\t    <BxBy>\n");
//      fprintf (fp,"# %12c  %12c  %12c\n",'t','delta_t', '<kin>');
    }else{
      fp = fopen("averages.dat","a");
    }
    fprintf (fp, "%12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", 
             g_time, g_dt, rhoe, pm, kin, my, E, BxBy);
    fclose(fp);
    first_call = 0;
  }
}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in/out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies on which side boundary conditions need 
 *                    to be assigned. side can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{ }


