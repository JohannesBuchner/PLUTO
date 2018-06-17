/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Super Time Stepping driver for integration of diffusion terms.

  Take one step in the solution of the diffusion equation 
  \f$ dU/dt = R(\partial^2U)  \f$ where R is a nonlinear right hand side
  involving second derivatives. 
  The super step is taken to be equal to the current time step
  \c dt and the number of substeps Nsts is given by solving the following
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
  \date    May 15, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#ifndef STS_NU     
 #define STS_NU  0.01  /**< Sets the nu parameter used in STS. */
#endif

#define STS_MAX_STEPS 1024

static void   STS_ComputeSubSteps(double, double tau[], int);
static double STS_FindRoot(double, double);
static double STS_CorrectTimeStep(int, double);
/* ********************************************************************* */
void STS (const Data *d, double dt, timeStep *Dts, Grid *grid)
/*!
 * Solve diffusion equation(s) using Super-Time-Stepping.
 *
 * \param [in,out]  d    pointer to Data structure
 * \param [in]      dt   the time step increment
 * \param [in,out]  Dts  pointer to timeStep structure  
 * \param [in]     grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int    i, j, k, nv, n, m;
  double N, ts[STS_MAX_STEPS];
  double dt_par, tau, tsave, invDt_par;
  static Data_Arr rhs;
  RBox  box;

  if (rhs == NULL) rhs  = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double); 
  RBoxDefine (IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &box);

/* --------------------------------------------------------
   0. Compute the conservative vector in order to start
      the cycle. 
      This step will be useless if the data structure 
      contains Vc as well as Uc (for future improvements).
   -------------------------------------------------------- */

  PrimToCons3D(d->Vc, d->Uc, &box);
  tsave = g_time;

/* --------------------------------------------------------
   1. Main STS Loop starts here
   -------------------------------------------------------- */

  m = 0;
  n = STS_MAX_STEPS;
  while (m < n){

    g_intStage = m + 1;
    Boundary(d, ALL_DIR, grid);
    invDt_par = ParabolicRHS(d, rhs, &box, NULL, SUPER_TIME_STEPPING, 1.0, grid);

  /* -----------------------------------------------------------
     1a. At the first step (m=0) compute (explicit) parabolic
         time step. Restriction on explicit time step should be
      
       [2*eta_x/dx^2 + 2*eta_y/dy^2 + 2*eta_z/dz^2]*dt < 1

      which in PLUTO is normally written as 
      
       [2/dtp_x + 2/dtp_y + 2/dtp_z]*dt/Ndim = 1/Ndim = cfl_par

      where Ndim is the number of spatial dimensions.
     ----------------------------------------------------------- */

    if (m == 0){
      Dts->invDt_par = invDt_par/(double)DIMENSIONS;  
      #ifdef PARALLEL
      MPI_Allreduce (&Dts->invDt_par, &tau, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      Dts->invDt_par = tau;
      #endif
      Dts->invDt_par = MAX(Dts->invDt_par, 1.e-18);
      dt_par = Dts->cfl_par/(2.0*Dts->invDt_par); /* explicit parabolic 
                                                   time step          */
    /* -------------------------------------------
        Compute the number of steps needed to fit 
        the supertimestep with the advection one.
       ------------------------------------------- */

      N = STS_FindRoot(dt_par, dt);
      N = floor(N+1.0);
      n = (int)N;

      Dts->Nsts = n;
      if (n > STS_MAX_STEPS){
        print ("! STS: the number of substeps (%d) is > %d\n",n, STS_MAX_STEPS); 
        QUIT_PLUTO(1);
      }

    /* ------------------------------------------
        Adjust explicit timestep to make Nsts an 
        integer number and compute the substeps.
       ------------------------------------------ */

      if (Dts->Nsts > 1){
        dt_par = STS_CorrectTimeStep(Dts->Nsts, dt);
        STS_ComputeSubSteps(dt_par, ts, Dts->Nsts);
      }
    
      if (Dts->Nsts == 1) ts[0] = dt;
    }   /* end if m == 0 */

    tau = ts[n-m-1];

  /* ------------------------------------------------------
     1b. Update cell-centered conservative variables
     ------------------------------------------------------ */
     
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
           (RESISTIVITY        == SUPER_TIME_STEPPING) || \
           (VISCOSITY          == SUPER_TIME_STEPPING) 
       d->Uc[k][j][i][ENG] += tau*rhs[k][j][i][ENG]; 
       #endif
      #endif
    }

  /* ------------------------------------------------------
     1c. Update staggered magnetic field variables
     ------------------------------------------------------ */

    #if (defined STAGGERED_MHD) && (RESISTIVITY == SUPER_TIME_STEPPING)
    CT_ResistiveEMF(d, 0, grid);
    CT_Update  (d, d->Vs, tau, grid);
    CT_AverageMagneticField (d->Vs, d->Uc, grid);
    #endif

  /* ------------------------------------------------------
     1d. Unflag zone tagged with the ENTROPY_SWITCH since
         only total energy can be evolved using STS
     ------------------------------------------------------ */

    #if ENTROPY_SWITCH
    TOT_LOOP(k,j,i) d->flag[k][j][i] &= ~FLAG_ENTROPY;
    #endif

  /* ----------------------------------------------
     1e. Convert conservative variables to primitive 
         for next iteration. Increment loop index.
     ---------------------------------------------- */

    ConsToPrim3D(d->Uc, d->Vc, d->flag, &box);
    g_time += ts[n-m-1];
    m++;
  }
  g_time = tsave;  /* restore initial time step */
}

/* ********************************************************************* */
void STS_ComputeSubSteps(double dtex, double tau[], int N)
/*!
 * Compute the single sub-step sequence (Eq. [2.9]).
 * N must be an integer by now.
 *
 *********************************************************************** */
{
  int i;
  double S=0.0;

  for (i = 0; i < N; i++) {
    tau[i] = dtex / ((-1.0 + STS_NU)*cos(((2.0*i+1.0)*CONST_PI)/(2.0*N)) 
                     + 1.0 + STS_NU);
    S += tau[i];
  }
}

/* ********************************************************************* */
double STS_FindRoot(double dt_exp, double dT)
/*!
 * Find the number of sub-steps N by solving Eq. (2.10) of AAG using a
 * Newton-Raphson scheme. 
 * Input to the function are:
 * 
 * \param [in]  dt_exp   the explicit time step
 * \param [in]  dt       the super-step.
 *
 *********************************************************************** */
{
  int k;  /* Iteration number */
  double a,b,c, scrh;
  double fN, N, dN, dfN;
  double db_dN, sqrt_nu = sqrt(STS_NU);

  k = 0;
  N = 1.0;
  a = (1.0 - sqrt_nu)/(1.0 + sqrt_nu);
  while(k < 128){
    b     = pow(a,2.0*N);
    c     = (1.0-b)/(1.0+b);    /* round bracket in Eq. [10] in AAG */
    db_dN = 2.0*log(a)*b;
    scrh  = c - N*2.0/((1.0+b)*(1.0+b))*db_dN;

    fN  = dT - 0.5*dt_exp/sqrt_nu*N*c;
    dfN =    - 0.5*dt_exp/sqrt_nu*scrh;
    dN  = fN/dfN; 
    
    N -= dN;
    k++;

    if (fabs(dN) < 1.e-5) return N;
  }

  print ("! STS_FindRoot: too many iterations\n");
  QUIT_PLUTO(1);
  return -1.0;
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
