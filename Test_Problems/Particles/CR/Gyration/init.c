/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file 
  \brief Single CR particle gyration (and drift) test.

  \author A. Mignone (mignone@ph.unito.it)

  \date   Jan 10, 2017

  \b References: \n
   - "A PARTICLE MODULE FOR THE PLUTO CODE: I - AN IMPLEMENTATION OF THE
      MHD-PIC EQUATIONS", Mignone et al. ApJS (2018)  [ Sec. 4.1 ]
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 *  Test particle gyration in EM fields.
 *
 *  The time steps for partciles and fluid scale as
 *
 *    dt_pa = Na*dx/vp        (if limited by translation)
 *    dt_pg = 1/(Ng*Omega)     (if limited by gyration)
 *    dt_f  = Ca*dx/lambda
 *
 *  where lambda = B/sqrt(rho)  (when p << 1), Omega = q*B/(m*c), Ng
 *  and Na are some integers, Ca is the Courant number.
 * 
 *  We set B = H*dx so the fluid time step is approximately the same
 *  at any resolution:
 *
 *    dt_pa = Na*dx/vp
 *    dt_pg = 1/(Ng*H*dx*q/mc)
 *    dt_f  = Ca*sqrt(rho)/H
 *
 *  The ratio between particle and fluid time steps is:
 *
 *    dt_pa / dt_f = Na*dx*H/(vp*Ca*sqrt(rho))
 *    dt_pg / dt_f = mc/(Ng*q*Ca*sqrt(rho)*dx)
 *
 *********************************************************************** */
{
  v[RHO] = 1.0;
  v[VX1] = g_inputParam[VFLUID_X]; 
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  v[PRS] = 1.0/g_gamma;
  v[TRC] = 0.0;

  #if PHYSICS == MHD || PHYSICS == RMHD
   v[BX1] = 0.0;
   v[BX2] = 0.0;
   v[BX3] = 1.0;

   v[AX1] = 0.0;
   v[AX2] = 0.0;
   v[AX3] = 0.0;
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
 *
 * - Set  {x, v}(t=0) in the Lab frame.
 * - Apply a Lorentz boost, so that E' = 0  (v/c = E/B = vgx/c)
 * - Solve particle motion in the primed system
 *   o transform initial condition (x,v,t) -> (x,v,t)'
 *   o solve equations
 * - Transform back into original system.
 *
 *********************************************************************** */
{
  int    dir;
  static int first_call = 1;
  char file_name[32];
  double Bz0, Bz1, vg, v1, u1;
  double c2 = PARTICLES_CR_C*PARTICLES_CR_C;
  double gamma_g, gamma, gamma1, Omega, Omega1; 
  double  x,  y,  z,  t,  vx,  vy, vz;
  double x1, y1, z1, t1, vx1, vy1, vz1, Ekin1;
  Particle *p;
  particleNode *curr;
  FILE *fp;

/* -- Compute important quantities -- */

  vg      = g_inputParam[VFLUID_X];      /* Fluid velocity in the x direction */
  gamma_g = 1.0/sqrt(1.0 - vg*vg/c2);    /* Lorentz boost (fluid velocity) */
  Bz0     = 1.0;                         /* Magnetic field in lab frame */
  Bz1     = Bz0/gamma_g;

  u1      = g_inputParam[PARTICLE_4VEL];  /* Particle 4-vel in co-moving frame */
  gamma1  = sqrt(1.0 + u1*u1/c2);  
  v1      = u1/gamma1;                    /* Particle 3-vel in co-moving frame */
  Omega1  = PARTICLES_CR_E_MC*Bz1/gamma1; /* Larmor Frequency in co-moving frame */

/* -- Loop over particles -- */

  PARTICLES_LOOP(curr, d->PHead){
    p = &(curr->p);

  /* -- Get particle coordinates (x and v) in Lab frame -- */
    
    x  = p->coord[IDIR];  y = p->coord[JDIR];  z = p->coord[KDIR];
    vx = p->speed[IDIR]; vy = p->speed[JDIR]; vz = p->speed[KDIR];
    t  = g_time;

  /* -- Compute coordinates in co-moving frame -- */

    x1 = gamma_g*(x - vg*t);
    y1 = y;
    z1 = z;
    t1 = gamma_g*(t - vg*x/c2);  

    vx1 = (vx - vg)/(1.0 - vx*vg/c2);
    vy1 = (vy     )/(1.0 - vx*vg/c2)/gamma_g;
    vz1 = (vz     )/(1.0 - vx*vg/c2)/gamma_g;

    gamma1 = 1.0/sqrt(1.0 - (vx1*vx1 + vy1*vy1 + vz1*vz1)/c2);
    Ekin1  = (gamma1 - 1.0)*c2;

    sprintf (file_name,"particle.%02d.dat", p->id);
    if (first_call) { 
      fp = fopen(file_name,"w");
      fprintf (fp,"# Omega1 = %12.6e\n",Omega1);
      fprintf (fp,"# R1     = %12.6e\n",v1/Omega1);
      fprintf (fp,"#\n");
      fprintf (fp,"#  %9s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n",
                   "t'","x'","y'","z'","vx'","vy'","vz'","E'");
      fprintf (fp,"# ---------------------------------------");
      fprintf (fp,"-----------------------------------------");
      fprintf (fp,"-----------------------------------------\n");  
    }else     fp = fopen(file_name,"a");
 
    fprintf (fp,"%12.6e  %12.6e  %12.6e  %12.6e  ",t1, x1, y1, z1);
    fprintf (fp,"%12.6e  %12.6e  %12.6e  %12.6e\n", vx1, vy1, vz1, Ekin1);

    fclose(fp);
  }

/* -------------------------------------------------------------
    Write exact solution in the comoving frame
   ------------------------------------------------------------- */

  if (first_call){
    double alpha;  
    double vx_ic, vy_ic, vz_ic, vx1_ic, vy1_ic, vz1_ic;

    gamma = 1.0/sqrt(1.0 - (vx*vx + vy*vy + vz*vz)/c2);
    Omega = PARTICLES_CR_E_MC*Bz0/gamma;

  /* -- Compute initial velocity components in Lab frame -- */

    alpha  = CONST_PI*g_inputParam[ALPHA]/180.0;
    vx1_ic = v1*cos(alpha);  /* Particle initial x-velocity, co-moving frame */
    vy1_ic = v1*sin(alpha);  /* Particle initial y-velocity, co-moving frame */
    vz1_ic = 0.0;

    print ("> Analysis(): Omega  = %f, Period  = %f, E  = %12.6e\n",
              Omega, 2.0*CONST_PI/Omega, (gamma - 1.0)*c2);
    print (">             Omega' = %f, Period' = %f, E' = %12.6e\n",
              Omega1, 2.0*CONST_PI/Omega1, (gamma1 - 1.0)*c2);

    fp = fopen("exact.dat","w");
    fprintf (fp,"# Omega1 = %12.6e\n",Omega1);
    fprintf (fp,"# R1     = %12.6e\n",v1/Omega1);
    fprintf (fp,"#\n");
    fprintf (fp,"#  %9s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n",
                   "t'","x'","y'","z'","vx'","vy'","vz'","E'");
    fprintf (fp,"# ---------------------------------------");
    fprintf (fp,"-----------------------------------------");
    fprintf (fp,"-----------------------------------------\n");

    double tstop = RuntimeGet()->tstop;
    for (t1 = 0.0; t1 < tstop; t1 += (tstop/1000.0)){
      alpha = Omega1*t1;
      x1  = (vy1_ic*(1.0 - cos(alpha)) + vx1_ic*sin(alpha))/Omega1;  
      y1  = (vx1_ic*(cos(alpha) - 1.0) + vy1_ic*sin(alpha))/Omega1;  
      z1  = 0.0;
      vx1 = (vy1_ic*sin(alpha) + vx1_ic*cos(alpha));
      vy1 = (vy1_ic*cos(alpha) - vx1_ic*sin(alpha));
      vz1 = 0.0;

      Ekin1  = 1.0/sqrt(1.0 - (vx1*vx1 + vy1*vy1 + vz1*vz1)/c2) - 1.0;
      Ekin1 *= c2;

      fprintf (fp,"%12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n",
                   t1, x1, y1, z1, vx1, vy1, vz1, Ekin1);
    }
    fclose(fp);
  }

  first_call = 0;
}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *
 *********************************************************************** */
{
  
}
