/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Shearing Box setup.

  Set initial conditions for a shearing box model in two or 
  three dimensions:
  \f[ 
    \rho = \left\{\begin{array}{ll}
      1              & \quad{\rm(unstratified)} \\ \noalign{\medskip}
     \exp(-z^2/H^2) & \quad{\rm(stratified)}
   \end{array}\right.
    ;\quad 
    v_y = -q\Omega x\,; \quad
    p = \rho c_s^2\,; \quad
    \vec{B} = \left\{\begin{array}{ll}
     (0,0,B_z)          & \quad{\rm(net\,flux)} \\ \noalign{\medskip}
     (0,0,B_z\sin 2\pi x/L_x) &  \quad{\rm(no\,net\,flux)} 
    \end{array}\right.
   \f]
  where, by default, the angular rotation frequency is \f$\Omega = 1\f$ 
  while the shear parameter is \f$ q = 3/2 \f$ (Keplerian rotation).
  In the stratified case the scale height is given by 
  \f$ H = \sqrt{2}c_s/\Omega \f$. 
  The parameters controlling the model are the elements of the 
  array \c g_inputParam:
 
  - <tt>g_inputParam[BETA]</tt>: sets the plasma beta (midplane),
     \f$ B_z = \sqrt{2c_s^2/\beta} \f$;
  - <tt>g_inputParam[CSOUND]</tt>: sets the speed of sound \f$ c_s \f$.
 
  The additional constants \c NET_FLUX and \c STRATIFICATION are 
  defined inside \c definitions.h and are used to:
 
  - \c NET_FLUX (\c YES/\c NO): set the magnetic field configuration 
     corresponding to a net magnetic flux (\c YES) or zero-net flux
     (\c NO);
  - \c STRATIFICATION (\c YES/\c NO): enable or disable stratified shearing 
    box models. With stratification, density and vertical gravity are changed.
    
  The \c Shearing_Box/ directory contains several configurations 
  for different purposes:

  <CENTER>
  Conf.|DIM|NET_FLUX|STRAT.|T. STEPPING|RECONSTR. | EOS      |Ref
  -----|---|--------|------|---------- |----------|----------|-------
   #01 | 2 | YES    | NO   | HANCOCK   |LINEAR    |ISOTHERMAL|  -
   #02 | 2 | YES    | NO   | ChTr      |PARABOLIC |ISOTHERMAL|  -
   #03 | 3 | YES    | NO   | RK2       |LINEAR    |ISOTHERMAL|[Bod08]
   #04 | 3 | YES    | NO   | RK2       |LINEAR    |ISOTHERMAL|[Bod08]
   #05 | 3 | YES    | NO   | ChTr      |PARABOLIC |ISOTHERMAL|[Bod08]
   #06 | 3 | YES    | NO   | ChTr      |PARABOLIC |ISOTHERMAL|[Bod08]
   #07 | 3 | NO     | NO   | RK2       |LINEAR    |ISOTHERMAL|[Mig12]
   #08 | 3 | NO     | NO   | RK2       |LINEAR    |ISOTHERMAL|[Mig12] (*)
   #09 | 3 | NO     | NO   | ChTr      |PARABOLIC |ISOTHERMAL|[Mig12]
   #10 | 3 | NO     | NO   | ChTr      |PARABOLIC |ISOTHERMAL|[Mig12] (*)
   #11 | 3 | NO     | YES  | HANCOCK   |LINEAR    |ISOTHERMAL|[Bod14] 
   #12 | 3 | NO     | YES  | HANCOCK   |LINEAR    |IDEAL     |[Bod12] (^)
   #13 | 3 | NO     | NO   | HANCOCK   |LINEAR    |IDEAL     |(**)
   #14 | 3 | YES    | NO   | RK3       |WENO3     |IDEAL     |(**)
  </CENTER>

  (*) used with FARGO to be compared to the previous configurations.
  
  (**) as conf. #9 but with IDEAL eos and FARGO.
  
  (^) This a simplified version of the [Bod12] setup.
      The actual configuration was also modifying the Riemann solver flux
      to guarantee zero mass flux across top and bottom boundary.

  \image html sb.06.png  "Density and magnetic field lines for configuration #06 at t = 50 using a grid resolution of 64 x 256 x 64."

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 2, 2017

  \b References:
     - [Bod08]: "Aspect ratio dependence in magnetorotational instability 
                 shearing box simulations" Bodo et al. A&A (2008) 487, 1
     - [Mig12]: "A conservative orbital advection scheme for simulations of
                 magnetized shear flows with the PLUTO code",
                 Mignone et al., A&A (2012) 545, Sec. 3.3.2
     - [Bod12]: "Magnetorotational turbulence in stratified shearing boxes
                 with perfect gas equation of state and finite thermal diffusivity",
                 Bodo et al., ApJ (2012) 761, 116
     - [Bod14]: "On the Convergence of Magnetorotational Turbulence in 
                 Stratified Isothermal Shearing Boxes" Bodo et al. (2014) 787 L13
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

//double sb_Omega = 1.0; /* Orbital frequency (global variable) */
//double sb_q     = 1.5; /* Shear parameter   (global variable) */


/* ********************************************************************* */
void Init (double *v, double x, double y, double z)
/*
 *
 *********************************************************************** */
{
  static int first_call = 1;
  double Bz0, rnd, dvy, cs, H;
  double Lx, Ly, Lz;
  double kx, ky, kz;

#ifndef SHEARINGBOX
  print ("! ShearingBox module has not been included.\n");
  print ("! Cannot continue.\n");
  QUIT_PLUTO(1);
#endif

/* -- compute domain sizes -- */

  Lx = g_domEnd[IDIR] - g_domBeg[IDIR]; kx = 2.0*CONST_PI/Lx;
  Ly = g_domEnd[JDIR] - g_domBeg[JDIR]; ky = 2.0*CONST_PI/Ly;
  Lz = g_domEnd[KDIR] - g_domBeg[KDIR]; kz = 2.0*CONST_PI/Lz;

/* -- get sound speed and pressure scale height -- */

  cs = g_inputParam[CSOUND];    /* sound speed */
//H  = sqrt(2.0)*cs/sb_Omega;   /* pressure scale height */
  H  = sqrt(2.0)*cs/SB_OMEGA;   /* pressure scale height */

/* -- seed random numbers differently for each processor -- */

  if (first_call == 1){
    srand(time(NULL) + prank);
    first_call = 0;
  }
  rnd = (double)(rand())/((double)RAND_MAX + 1.0);

/* -- set velocity perturbation [1]: random noise -- */

/*  dvy = 0.01*cs*rnd; */

/* -- set velocity perturbation [2], two harmonics for each direction -- */

  dvy  = sin(kx*x + 0.20) + sin(2.0*kx*x - 0.37);
  dvy *= sin(ky*y + 0.13) + sin(2.0*ky*y + 0.04);
  dvy *= sin(kz*z + 0.56) + sin(2.0*kz*z + 0.62);
  dvy *= 0.01*cs/8.0;

/* -- in 2D we don't use any perturbation -- */

  #if DIMENSIONS == 2
   dvy = 0.0;
  #endif

/* -- set initial condition -- */

  #if STRATIFICATION == YES
   v[RHO] = exp(-z*z/(H*H));
  #else
   v[RHO] = 1.0;
  #endif

  v[VX1] = 0.0;
//v[VX2] = -sb_q*sb_Omega*x + dvy;
  v[VX2] = -SB_Q*SB_OMEGA*x + dvy;
  v[VX3] = 0.0;
  #if EOS == IDEAL
   v[PRS] = cs*cs*v[RHO];
  #elif EOS == ISOTHERMAL
   g_isoSoundSpeed = cs;
  #endif

  v[TRC] = 0.0;

/* ----------------------------------------------------------------
    The magnetic field amplitude is set by the parameter
    beta = 2p/B^2 = 2*rho*c^2/b0^2   -->   b0 = c*sqrt(2/beta)
    where it is assumed that rho = 1 (midplane).
   ---------------------------------------------------------------- */

  #if PHYSICS == MHD 
   Bz0 = cs*sqrt(2.0/g_inputParam[BETA]);

   #if  NET_FLUX  == YES  /* -- Net flux, constant vertical field Bz = B0 -- */

    v[BX1] = 0.0;
    v[BX2] = 0.0;
    v[BX3] = Bz0;

    v[AX1] = 0.0;
    v[AX2] = Bz0*x;
    v[AX3] = 0.0;

   #else /* -- Zero net flux, Bz = B0*sin(2*pi*x/Lx) -- */

   v[BX1] = 0.0;
   v[BX2] = 0.0;
   v[BX3] = Bz0*sin(kx*x);

   v[AX1] = 0.0;
   v[AX2] = -Bz0*cos(kx*x)/kx;
   v[AX3] = 0.0;

   #endif

   #if DIMENSIONS == 2  /* 2D Case only for testing */
    v[BX1] = Bz0*sin(ky*y);
    v[BX2] = 0.0;
    v[BX3] = 0.0;

    v[AX1] = 0.0;
    v[AX2] = 0.0;
    v[AX3] = -Bz0*cos(ky*y)/ky;
   #endif
  #endif
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 * Compute volume-integrated magnetic pressure, Maxwell and
 * Reynolds stresses. Save them to "averages.dat"
 *
 *********************************************************************** */
{
  int    i,j,k;
  static int first_call=1;
  double *dx, *dy, *dz;
  double *x, *y, *z;
  double Lx, Ly, Lz, scrh;
  double pm, Mxy, Rxy, tot_vol, aM, aR, dV;
  FILE *fp;

/* -- compute domain sizes and box volume -- */

  Lx = g_domEnd[IDIR] - g_domBeg[IDIR];
  Ly = g_domEnd[JDIR] - g_domBeg[JDIR];
  Lz = g_domEnd[KDIR] - g_domBeg[KDIR];

  tot_vol = Lx*Ly*Lz;

/* -- pointers to mesh spacing and coordinates -- */

  dx = grid->dx[IDIR]; dy = grid->dx[JDIR]; dz = grid->dx[KDIR];
   x = grid->x[IDIR];   y = grid->x[JDIR];   z = grid->x[KDIR];

/* ------------------------------------------------------------
    Main analysis loop.
    Compute volume-integrated magnetic pressure and stresses
   ------------------------------------------------------------ */

  pm = Mxy = Rxy = 0.0;
  DOM_LOOP(k,j,i){
    dV  = dx[i]*dy[j]*dz[k];

  /* -- magnetic pressure -- */

    pm += 0.5*(EXPAND(  d->Vc[BX1][k][j][i]*d->Vc[BX1][k][j][i],
                      + d->Vc[BX2][k][j][i]*d->Vc[BX2][k][j][i],
                      + d->Vc[BX3][k][j][i]*d->Vc[BX3][k][j][i]))*dV;

  /* -- Maxwell and Reynolds stresses -- */

    aM = d->Vc[BX1][k][j][i]*d->Vc[BX2][k][j][i];
    aR = d->Vc[RHO][k][j][i]*d->Vc[VX1][k][j][i]*
        (d->Vc[VX2][k][j][i] + SB_Q*SB_OMEGA*x[i]);

    Mxy   += aM*dV;
    Rxy   += aR*dV;
  }

/* -- divide summations by total volume -- */

  pm  /= tot_vol;
  Mxy /= tot_vol;
  Rxy /= tot_vol;

/* -- parallel reduce -- */

  #ifdef PARALLEL
   MPI_Allreduce (&pm, &scrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   pm = scrh;

   MPI_Allreduce (&Mxy, &scrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   Mxy = scrh;

   MPI_Allreduce (&Rxy, &scrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   Rxy = scrh;

   MPI_Barrier (MPI_COMM_WORLD);
  #endif

/* -- only proc #0 writes ascii data-file to disk -- */

  if (prank == 0){
    static double tpos = -1.0;
    if (g_stepNumber == 0){   /* -- open for writing if initial step -- */
      char fname[64];
      sprintf (fname, "%s/averages.dat",RuntimeGet()->output_dir); 
      fp = fopen(fname,"w");
      fprintf (fp,"# %4s  %12s  %12s  %12s  %12s\n",
                   "time","  step  "," <B^2/2> "," <Bx*By> ","<rho*ux*duy>");
      first_call = 0;
    }else{
      if (tpos < 0.0){  /* obtain time coordinate of last written line */
        char   sline[512];
        fp = fopen("averages.dat","r");
        if (fp == NULL){
          print ("! Analysis(): file averages.dat not found\n");
          QUIT_PLUTO(1);
        }
        while (fgets(sline, 512, fp))  {
        }
        sscanf(sline, "%lf\n",&tpos);
        fclose(fp);
      }
      fp = fopen("averages.dat","a");
    }
    if (g_time > tpos){ /* -- write -- */
      fprintf (fp, "%12.6e  %4d  %12.6e  %12.6e  %12.6e\n", 
               g_time, g_stepNumber, pm, Mxy, Rxy);
    }
    fclose(fp);
  }
}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  The user-defined boundary is used to impose stress-free boundary
 *  and purely vertical magnetic field at the top and bottom boundaries,
 *  as done in Bodo et al. (2012).
 *  In addition, constant temperature and hydrostatic balance are imposed.
 *  For instance, at the bottom boundary, one has:
 *  \f[
 *     \left\{\begin{array}{lcl}
 *      \DS  \frac{dp}{dz} &=& \rho g_z  \\ \noalign{\medskip}
 *            p            &=& \rho c_s^2
 *     \end{array}\right.
 *     \qquad\Longrightarrow\qquad
 *     \frac{p_{k+1}-p_k}{\Delta z} = \frac{p_{k} + p_{k+1}}{2c_s^2}g_z
 *  \f]
 *  where \f$g_z\f$ is the value of gravity at the lower boundary.
 *  Solving for \f$p_k\f$ at the bottom boundary where \f$k=k_b-1\f$
 *  gives:
 *  \f[
 *     \left\{\begin{array}{lcl}
 *       p_k   &=& \DS  p_{k+1} \frac{1-a}{1+a}  \\ \noalign{\medskip}
 *      \rho_k &=& \DS \frac{p_k}{c_s^2}
 *     \end{array}\right.
 *     \qquad\mathrm{where}\qquad
 *     a = \frac{\Delta z g_z}{2c_s^2} > 0 
 *  \f]
 *  where, for simplicity, we keep constant temperature in the ghost
 *  zones rather than at the boundary interface (this seems to give a
 *  more stable behavior and avoids negative densities).
 *  A similar treatment holds at the top boundary.
 *********************************************************************** */
{
  int   i, j, k, nv;
  double *z, *dz;
  double ***rho = d->Vc[RHO];
#if HAVE_ENERGY
  double ***prs = d->Vc[PRS];
#endif
  double cs = g_inputParam[CSOUND], cs2 = cs*cs, H;
  EXPAND(double ***vx = d->Vc[VX1];  ,
         double ***vy = d->Vc[VX2];  ,
         double ***vz = d->Vc[VX3];)
  EXPAND(double ***Bx = d->Vc[BX1];  ,
         double ***By = d->Vc[BX2];  ,
         double ***Bz = d->Vc[BX3];)
#ifdef STAGGERED_MHD
  EXPAND(double ***Bxs = d->Vs[BX1s];  ,
         double ***Bys = d->Vs[BX2s];  ,
         double ***Bzs = d->Vs[BX3s];)
#endif

  H   = sqrt(2.0)*cs/SB_OMEGA;   /* pressure scale height */

  z  = grid->x[KDIR];
  dz = grid->dx[KDIR];

  if (side == 0) {    /* -- Density threshold -- */
    DOM_LOOP(k,j,i){
      if (d->Vc[RHO][k][j][i] < 1.e-5) d->Vc[RHO][k][j][i] = 1.e-5;
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    double gz, a;
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        gz = -g_domBeg[KDIR]*SB_OMEGA*SB_OMEGA;  /* > 0 */
        a  = 0.5*dz[k]*gz/cs2;
#if EOS == IDEAL
        prs[k][j][i] = prs[KBEG][j][i]*(1.0 - a)/(1.0 + a);
        rho[k][j][i] = prs[k][j][i]/cs2; 
#elif EOS == ISOTHERMAL
        rho[k][j][i] = rho[KBEG][j][i]*(1.0 - a)/(1.0 + a);
#endif
        EXPAND(vx[k][j][i]  =  vx[2*KBEG-k-1][j][i];  ,
               vy[k][j][i]  =  vy[2*KBEG-k-1][j][i];  ,
               vz[k][j][i]  = -vz[2*KBEG-k-1][j][i];)

        EXPAND(Bx[k][j][i]  = -Bx[2*KBEG-k-1][j][i];  ,
               By[k][j][i]  = -By[2*KBEG-k-1][j][i];  ,
               Bz[k][j][i]  =  Bz[2*KBEG-k-1][j][i];)
      }
    }
#ifdef STAGGERED_MHD
    if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i) Bxs[k][j][i] = -Bxs[2*KBEG-k-1][j][i];
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i) Bys[k][j][i] = -Bys[2*KBEG-k-1][j][i];
    }
#endif
  }

  if (side == X3_END) {  /* -- X3_END boundary -- */
    double gz, a;
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        gz = -g_domEnd[KDIR]*SB_OMEGA*SB_OMEGA;  /* < 0 */
        a  = 0.5*dz[k]*gz/cs2;
#if EOS == IDEAL
        prs[k][j][i] = prs[KEND][j][i]*(1.0 + a)/(1.0 - a);
        rho[k][j][i] = prs[k][j][i]/cs2; 
#elif EOS == ISOTHERMAL
        rho[k][j][i] = rho[KEND][j][i]*(1.0 + a)/(1.0 - a);
#endif
        EXPAND(vx[k][j][i]  =  vx[2*KEND-k+1][j][i];  ,
               vy[k][j][i]  =  vy[2*KEND-k+1][j][i];  ,
               vz[k][j][i]  = -vz[2*KEND-k+1][j][i];)

        EXPAND(Bx[k][j][i]  = -Bx[2*KEND-k+1][j][i];  ,
               By[k][j][i]  = -By[2*KEND-k+1][j][i];  ,
               Bz[k][j][i]  =  Bz[2*KEND-k+1][j][i];)
      }
    }
#ifdef STAGGERED_MHD
    if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i) Bxs[k][j][i] = -Bxs[2*KEND-k+1][j][i];
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i) Bys[k][j][i] = -Bys[2*KEND-k+1][j][i];
    }
#endif
  }
	
}

/* ********************************************************************* */
void BodyForceVector (double *v, double *g, double x, double y, double z)
/*!
 *  Include gravitational force in the shearing box module.
 *  Coriolis terms are included elsewhere.
 *
 *  Note: with FARGO, gravity in the x-direction must not be included.
 *  
 *********************************************************************** */
{
  double om2 = SB_OMEGA*SB_OMEGA;
  double zb  = g_domBeg[KDIR], ze  = g_domEnd[KDIR];

#ifdef FARGO
  g[IDIR] = 0.0;
#else
  g[IDIR] = om2*2.0*SB_Q*x;
#endif
  g[JDIR] = 0.0;

#if STRATIFICATION == YES
  g[KDIR] = -om2*z;

/* -----------------------------------------------------------
    With CTU and a reflective vertical boundary, gravity must
    be symmetrized since the predictor step of CTU does not
    not preserve the symmetry.
   ----------------------------------------------------------- */

  #ifdef CTU   
  if      (z > ze) g[KDIR] += 2.0*om2*ze;
  else if (z < zb) g[KDIR] += 2.0*om2*zb;
  #endif

#else
  g[KDIR] = 0.0; 
#endif
}
/* ********************************************************************* */
double BodyForcePotential(double x, double y, double z)
/*
 * Return the gravitational potential as function of the coordinates.
 *
 *********************************************************************** */
{
  double om2 = SB_OMEGA*SB_OMEGA;
  double zb  = g_domBeg[KDIR], ze = g_domEnd[KDIR];
  double psi;
 
#if STRATIFICATION == YES
  psi = om2*(0.5*z*z - SB_Q*x*x);
  #ifdef CTU /* see note in BodyForceVector() */
  if (z > ze) psi -= 2.0*om2*ze*(z - ze);
  if (z < zb) psi -= 2.0*om2*zb*(z - zb); 
  #endif
  return psi;
#else
  return -om2*SB_Q*x*x;
#endif
}

/* ************************************************************** */
double FARGO_SetVelocity(double x1, double x2)
/*!
 *   Compute the shear angular velocity to be subtracted from 
 *   the HD or MHD equations.
 * 
 **************************************************************** */
{
  return -SB_Q*SB_OMEGA*x1;
}
