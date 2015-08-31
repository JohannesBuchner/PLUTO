/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Super Time Stepping driver for integration of diffusion terms.

  Take one step in the solution of the diffusion equation 
  \f$ dU/dt = R(\partial^2U)  \f$ where R is a nonlinear right hand side
  involving second derivatives. 
  The super step is taken to be equal to the current time step
  ::g_dt and the number of substeps Nsts is given by solving the following
  nonlinear equation:
  \f[
    \frac{\Delta t}{\Delta t_{\rm par}}
     =
    \frac{N}{2\sqrt{\nu}}\frac{(1+\sqrt{\nu})^{2N} - (1-\sqrt{\nu})^{2N}}
                              {(1+\sqrt{\nu})^{2N} + (1-\sqrt{\nu})^{2N}}
  \f]
  where \f$\nu\f$ is set by the macro ::STS_NU (default equal to 0.01).
  The previous relation is given by Eq. [2.10] of Alexiades et al. with the
  explicit parabolic time step \f$\Delta t_{\rm par}\f$ being computed from
  \f[
    \frac{2}{N_d} \max\left[  \frac{\eta_x}{\Delta x^2} 
                            + \frac{\eta_y}{\Delta y^2} 
                            + \frac{\eta_z}{\Delta z^2} \right] 
    \Delta t_{\rm par} = C_p < \frac{1}{N_d}
  \f]
  where \f$C_p\f$ is the parabolic Courant number, \f$ N_d \f$ is the number
  of spatial dimensions and the maximum of the square
  bracket is computed during the call to ::ParabolicRHS.
  
  This function is called in an operator-split way before/after advection has 
  been carried out.
 
  \b References
     - Alexiades, V., Amiez, A., \& Gremaud E.-A. 1996, 
       Com. Num. Meth. Eng., 12, 31

  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos
  \date    Aug 27, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#ifndef STS_NU     
 #define STS_NU  0.01  /**< Sets the nu parameter used in STS. */
#endif

#define STS_MAX_STEPS 1024

static void   STS_ComputeSubSteps(double, double tau[], int);
static double STS_FindRoot(double, double, double);
static double STS_CorrectTimeStep(int, double);
/* ********************************************************************* */
void STS (const Data *d, Time_Step *Dts, Grid *grid)
/*!
 * Solve diffusion equation(s) using Super-Time-Stepping.
 *
 * \param [in,out]  d    pointer to Data structure
 * \param [in,out]  Dts  pointer to Time_Step structure  
 * \param [in]     grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int    i, j, k, nv, n, m;
  double N, ts[STS_MAX_STEPS];
  double dt_par, tau, tsave, inv_dtp;
  static Data_Arr rhs;
  RBox *box = GetRBox(DOM, CENTER);

  if (rhs == NULL) rhs  = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double); 

/* -------------------------------------------------------------
    Compute the conservative vector in order to start the cycle. 
    This step will be useless if the data structure 
    contains Vc as well as Uc (for future improvements).
   ------------------------------------------------------------- */

  PrimToCons3D(d->Vc, d->Uc, box);
  tsave = g_time;

/* ------------------------------------------------------------
               Main STS Loop starts here
   ------------------------------------------------------------ */

  m = 0;
  n = STS_MAX_STEPS;
  while (m < n){

    g_intStage = m + 1;
    Boundary(d, ALL_DIR, grid);
    inv_dtp = ParabolicRHS(d, rhs, 1.0, grid); 
    
  /* --------------------------------------------------------------
      At the first step (m=0) compute (explicit) parabolic time 
      step. Restriction on explicit time step should be
      
       [2*eta_x/dx^2 + 2*eta_y/dy^2 + 2*eta_z/dz^2]*dt < 1

      which in PLUTO is normally written as 
      
       [2/dtp_x + 2/dtp_y + 2/dtp_z]*dt/Ndim = 1/Ndim = cfl_par

      where Ndim is the number of spatial dimensions.
     -------------------------------------------------------------- */

    if (m == 0){
      Dts->inv_dtp = inv_dtp/(double)DIMENSIONS;  
      #ifdef PARALLEL
       MPI_Allreduce (&Dts->inv_dtp, &tau, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
       Dts->inv_dtp = tau;
      #endif
      Dts->inv_dtp = MAX(Dts->inv_dtp, 1.e-18);
      dt_par = Dts->cfl_par/(2.0*Dts->inv_dtp); /* explicit parabolic 
                                                   time step          */
    /* -------------------------------------------
        Compute the number of steps needed to fit 
        the supertimestep with the advection one.
       ------------------------------------------- */

      N = STS_FindRoot(1.0, dt_par, g_dt);
      N = floor(N+1.0);
      n = (int)N;

      Dts->Nsts = n;
      if (n > STS_MAX_STEPS){
        print1 ("! STS: the number of substeps (%d) is > %d\n",n, STS_MAX_STEPS); 
        QUIT_PLUTO(1);
      }

    /* ------------------------------------------
        Adjust explicit timestep to make Nsts an 
        integer number and compute the substeps.
       ------------------------------------------ */

      if (Dts->Nsts > 1){
        dt_par = STS_CorrectTimeStep(Dts->Nsts, g_dt);
        STS_ComputeSubSteps(dt_par, ts, Dts->Nsts);
      }
    
      if (Dts->Nsts == 1) ts[0] = g_dt;

    }   /* end if m == 0 */

    tau = ts[n-m-1];

  /* -----------------------------------------------
      Update cell-centered conservative variables
     ----------------------------------------------- */
     
    DOM_LOOP (k,j,i){
      #if VISCOSITY == SUPER_TIME_STEPPING
       EXPAND(d->Uc[k][j][i][MX1] += tau*rhs[k][j][i][MX1];  ,
              d->Uc[k][j][i][MX2] += tau*rhs[k][j][i][MX2];  ,
              d->Uc[k][j][i][MX3] += tau*rhs[k][j][i][MX3];)
      #endif
      #if (RESISTIVITY == SUPER_TIME_STEPPING)
       EXPAND(d->Uc[k][j][i][BX1] += tau*rhs[k][j][i][BX1];  ,
              d->Uc[k][j][i][BX2] += tau*rhs[k][j][i][BX2];  ,
              d->Uc[k][j][i][BX3] += tau*rhs[k][j][i][BX3];)
      #endif
      #if HAVE_ENERGY
       #if (THERMAL_CONDUCTION == SUPER_TIME_STEPPING) || \
           (RESISTIVITY      == SUPER_TIME_STEPPING) || \
           (VISCOSITY          == SUPER_TIME_STEPPING) 
        d->Uc[k][j][i][ENG] += tau*rhs[k][j][i][ENG]; 
       #endif
      #endif
    }
    #if (defined STAGGERED_MHD) && (RESISTIVITY == SUPER_TIME_STEPPING)

  /* -----------------------------------------------
      Update staggered magnetic field variables
     ----------------------------------------------- */

/* This piece of code should be used for simple Cartesian 2D
   configurations with constant resistivity.
   It shows that the stability limit of the resistive part of the
   induction equation is larger if the right hand side is written 
   as Laplacian(B) rather than curl(J) 
   
double ***Bx, ***By, ***Bz, Jz, dx, dy, dt;
static double ***Ez, ***Rx, ***Ry;
if (Ez == NULL) {
  Ez = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  Rx = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  Ry = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
}
Bx = d->Vs[BX1s];
By = d->Vs[BX2s];
dx = grid[IDIR].dx[IBEG];
dy = grid[JDIR].dx[JBEG];

dt = tau;
k=0;
for (j = 1; j < NX2_TOT-1; j++) for (i = 1; i < NX1_TOT-1; i++){
 Jz = (By[k][j][i+1] - By[k][j][i])/dx - (Bx[k][j+1][i] - Bx[k][j][i])/dy;
 Ez[k][j][i] = Jz*g_inputParam[ETAZ];
 Rx[k][j][i] =  (Bx[k][j][i+1] - 2.0*Bx[k][j][i] + Bx[k][j][i-1])/(dx*dx)
              + (Bx[k][j+1][i] - 2.0*Bx[k][j][i] + Bx[k][j-1][i])/(dy*dy);

 Ry[k][j][i] =   (By[k][j][i+1] - 2.0*By[k][j][i] + By[k][j][i-1])/(dx*dx)
               + (By[k][j+1][i] - 2.0*By[k][j][i] + By[k][j-1][i])/(dy*dy);
}

JDOM_LOOP(j) for (i = IBEG-1; i <= IEND; i++){
// Bx[k][j][i] += - dt/dy*(Ez[k][j][i] - Ez[k][j-1][i]);
 Bx[k][j][i] +=   dt*g_inputParam[ETAZ]*Rx[k][j][i];
}

for (j = JBEG-1; j <= JEND; j++) IDOM_LOOP(i){
// By[k][j][i] += dt/dx*(Ez[k][j][i] - Ez[k][j][i-1]);
 By[k][j][i] +=   dt*g_inputParam[ETAZ]*Ry[k][j][i];
}
*/
     CT_Update  (d, d->Vs, tau, grid);
     CT_AverageMagneticField (d->Vs, d->Uc, grid);
    #endif

  /* --------------------------------------------------
      Unflag zone tagged with the ENTROPY_SWITCH since
      only total energy can be evolved using STS
     -------------------------------------------------- */

    #if ENTROPY_SWITCH
     TOT_LOOP(k,j,i) d->flag[k][j][i] &= ~FLAG_ENTROPY;
    #endif

  /* ----------------------------------------------
      convert conservative variables to primitive 
      for next iteration. Increment loop index.
     ---------------------------------------------- */

    ConsToPrim3D(d->Uc, d->Vc, d->flag, box);
    g_time += ts[n-m-1];
    m++;
  }
  g_time = tsave;  /* restore initial time step */
}

/* ********************************************************************* */
void STS_ComputeSubSteps(double dtex, double tau[], int ssorder)
/*
 *
 *********************************************************************** */
{
  int i;
  double S=0.0;

  for (i = 0; i < ssorder; i++) {
    tau[i] = dtex / ((-1.0 + STS_NU)*cos(((2.0*i+1.0)*CONST_PI)/(2.0*ssorder)) 
                     + 1.0 + STS_NU);
    S += tau[i];
  }
}

/* ********************************************************************* */
double STS_FindRoot(double x0, double dtr, double dta)
/*
 *
 *********************************************************************** */
{
  double a,b,c;
  double Ns, Ns1;
  int n;

  n = 0;

  Ns  = x0+1.0;
  Ns1 = x0;
   
  while(fabs(Ns-Ns1) >= 1.0e-5){
    Ns = Ns1;
    a = (1.0-sqrt(STS_NU))/(1.0+sqrt(STS_NU));
    b = pow(a,2.0*Ns);
    c = (1.0-b)/(1.0+b);
    Ns1 = Ns + (dta - dtr*Ns/(2.0*sqrt(STS_NU))*c)
              /(dtr/(2.0*sqrt(STS_NU))*(c-2.0*Ns*b*log(a)*(1.0+c)/(1.0+b)));
    n += 1;
    if (n == 128){
      print1 ("! STS_FindRoot: max number of iterations exceeded");
      QUIT_PLUTO(1);
    }
  }
  return(Ns);
}

/* ********************************************************************* */
double STS_CorrectTimeStep(int n0, double dta)
/*
 *
 *********************************************************************** */
{
  double a,b,c;
  double dtr;

  a = (1.0-sqrt(STS_NU))/(1.0+sqrt(STS_NU));
  b = pow(a,2.0*n0);
  c = (1.0-b)/(1.0+b);

  dtr = dta*2.0*sqrt(STS_NU)/(n0*c);
  return(dtr);
}
#undef STS_MAX_STEPS 
