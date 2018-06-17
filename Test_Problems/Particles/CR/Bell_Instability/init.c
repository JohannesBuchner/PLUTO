/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Bell instability setup.

  Initialize fluid variables using Eqns. of [MVBM18], sect. 4.4. using a
  1D configuration.
  The configuration is then rotated around the y- and z-axis in a way similar
  to Mignone, Tzeferacos \& Bodo, JCP (2010).
  The initial phase is computed as \f$ \phi = k_0x_0 = \vec{k}\cdot\vec{x} \f$
  (which is invariant under rotations),  where \f$ \vec{k} \f$ and \f$ \vec{x} \f$
  are the wavevector and coordinate
  vector in the rotated system, while \f$ k_0 \f$ and \f$ x_0 \f$ are
  the original (unrotated) 1D vectors.

  The amount of rotation depends on the dimensionality of the problem
  and is uniquely specified by the domain size in the three directions
  Lx, Ly, Lz:
  \f[
    \tan(\alpha) = \frac{L_x}{L_y} = \frac{k_y}{k_x} \,,\quad
    \tan(\beta)  = \frac{L_x}{L_z} = \frac{k_z}{k_x}
  \f]
  where \f$ k_x = 2\pi/L_x,\, k_y = 2\pi/L_y\f$ and \f$ k_z = 2\pi/L_z \f$.
  In such a way the modulus of the wavevector is

  \f[
  |\vec{k}| = \sqrt{k_x^2 + k_y^2 + k_z^2}
            = 2\pi\sqrt{\frac{1}{L_x^2} + \frac{1}{L_y^2} + \frac{1}{L_z^2}}
  \f]

  In order to have one wavelength, Ly and Lz must be chosen so that
  \f$ |\vec{k}| = 2\pi \f$.
  Thus in 1D, 2D and 3D we choose:

  <CENTER>
  Dim |        Lx         |         Ly         |      Lz
  ----|-------------------|--------------------|------------
   1  |           1       |          -         |      -
   2  | \f$ \sqrt{5} \f$  | \f$ \sqrt{5}/2 \f$ |      -
   3  |           3       |        3/2         |     3/2
  </CENTER>

  \author A. Mignone (mignone@ph.unito.it)

  \date   March 20, 2018
  
  \b References: \n
   - [MVBM18] "A PARTICLE MODULE FOR THE PLUTO CODE: I - AN IMPLEMENTATION OF
               THE MHD-PIC EQUATIONS", Mignone etal.ApJS (2018)  [ Sec. 4.4 ]
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

void VectorRotate(double *v, int s);

/* ********************************************************************* */
void Init (double *v, double x, double y, double z)
/*! 
 *
 *********************************************************************** */
{
  static int first_call = 1;
  double xvec[3] = {x, y, z};
  double kvec[3];
  double Lx = g_domEnd[IDIR] - g_domBeg[IDIR];
  double Ly = g_domEnd[JDIR] - g_domBeg[JDIR];
  double Lz = g_domEnd[KDIR] - g_domBeg[KDIR];
  double theta   = asin(g_inputParam[EPSILON]);
  double k0, kx, ky=0.0, kz=0.0, phi;
  double Bperp   = g_inputParam[BPERP_AMPL];

  if (first_call){
    srand48(prank);
  }

/* -- Set wave vectorys and phase -- */

  D_EXPAND(kx = 2.0*CONST_PI/Lx;  ,
           ky = 2.0*CONST_PI/Ly;  ,
           kz = 2.0*CONST_PI/Lz;)
//kz = 0.0;
  if (fabs(Lx-1.0) < 1.e-9) ky = kz = 0.0; /** Special Case */

  phi = kx*x + ky*y + kz*z;
  k0  = sqrt(kx*kx + ky*ky + kz*kz);

/* -- Assign scalar quantities -- */

  v[RHO] = 1.0;
  v[PRS] = 1.0;

#if LINEAR_SETUP == TRUE
  if (first_call){
    print ("  ---------------------------------------------\n");
    print ("  normal   = [%f, %f, %f]\n",kx/k0,ky/k0,kz/k0);
    print ("  k0/(2pi) = %f\n",k0/(2.0*CONST_PI));
    print ("  lambda   = %f\n",2.0*CONST_PI/k0);
    print ("  ---------------------------------------------\n");
  }

/* --  Assign vector components in 1D, then rotate -- */

  v[VX1] =  0.0;
  v[VX2] =  Bperp*sin(phi - theta);
  v[VX3] = -Bperp*cos(phi - theta);
  VectorRotate(v+VX1, 1);

  v[BX1] = 1.0;
  v[BX2] = Bperp*cos(phi);
  v[BX3] = Bperp*sin(phi);
  VectorRotate(v+BX1, 1);

/* -- Vector potential is best assigned in 1D system using 1D coordinates -- */

  xvec[IDIR] = x; kvec[IDIR] = kx; 
  xvec[JDIR] = y; kvec[JDIR] = ky;
  xvec[KDIR] = z; kvec[KDIR] = kz;
  VectorRotate(xvec,-1);  /* Apply inverse rotation */
  VectorRotate(kvec,-1);  /* Apply inverse rotation */

  v[AX1] = 0.0;
  v[AX2] = -Bperp*cos(xvec[0]*kvec[0])/kvec[0];
  v[AX3] = -Bperp*sin(xvec[0]*kvec[0])/kvec[0] + xvec[1];
  VectorRotate(v+AX1, 1);

#else
  v[RHO] = 1.0;
  #if HAVE_ENERGY
   v[PRS] = 1.0;
  #endif

  lambda = (g_domEnd[IDIR]-g_domBeg[IDIR])/5.0;
  k      = 2.0*CONST_PI/lambda;

  v[VX1] = Bperp*RandomNumber(-1,1);
  v[VX2] = Bperp*RandomNumber(-1,1);
  v[VX3] = Bperp*RandomNumber(-1,1);

  v[TRC] = 0.0;

  v[BX1] = 1.0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;
#endif

  first_call = 0;
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
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
  int i, j, k;
  time_t tnow;
  double Emag_max, Emag_perp_ave, scrh, dvol;
  double By_max, By_min;
  double v[3];
  double *dx = grid->dx[IDIR];
  double *dy = grid->dx[JDIR];
  double *dz = grid->dx[KDIR];
  double vol =  (g_domEnd[IDIR]-g_domBeg[IDIR])
               *(g_domEnd[JDIR]-g_domBeg[JDIR])
               *(g_domEnd[KDIR]-g_domBeg[KDIR]);
  
  Emag_max = Emag_perp_ave = 0.0;
  By_max   = By_min   = 0.0;

  DOM_LOOP(k,j,i){
    dvol = dx[i]*dy[j]*dz[k];

    v[0] = d->Vc[BX1][k][j][i];
    v[1] = d->Vc[BX2][k][j][i];
    v[2] = d->Vc[BX3][k][j][i];

    VectorRotate(v,-1);

    scrh = 0.5*(v[1]*v[1] + v[2]*v[2]);
    Emag_perp_ave += scrh*dvol;
/*
    scrh = 0.5*(  d->Vc[BX1][k][j][i]*d->Vc[BX1][k][j][i]
                + d->Vc[BX2][k][j][i]*d->Vc[BX2][k][j][i]
                + d->Vc[BX3][k][j][i]*d->Vc[BX3][k][j][i]);
    
    Emag_perp_ave += (scrh-0.5)*dvol;
*/
    Emag_max  = MAX(Emag_max, scrh);
    By_max    = MAX(d->Vc[BX2][k][j][i], By_max);
    By_min    = MIN(d->Vc[BX2][k][j][i], By_min);
  }
  Emag_perp_ave /= vol;

#ifdef PARALLEL
  MPI_Allreduce (&Emag_perp_ave, &scrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Emag_perp_ave = scrh;

  MPI_Allreduce (&Emag_max, &scrh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  Emag_max = scrh;

  MPI_Allreduce (&By_max, &scrh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  By_max = scrh;

  MPI_Allreduce (&By_min, &scrh, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  By_min = scrh;
#endif

/* ---- Write ascii file "averages.dat" to disk ---- */

  if (prank == 0){
    char fname[512];
    static double tpos = -1.0;
    FILE *fp;
    sprintf (fname, "%s/averages.dat",RuntimeGet()->output_dir);
    if (g_stepNumber == 0){

    /* -------------------------------------------------------
        Open for writing only when we’re starting from
        beginning of computation. Write Header.
       ------------------------------------------------------- */

      fp = fopen(fname,"w");
      time(&tnow);
      fprintf (fp,"# File created on %s", asctime(localtime(&tnow)));
      fprintf (fp,"# Epsilon = %f\n",g_inputParam[EPSILON]);
      fprintf (fp,"# \n");
      fprintf (fp,"# Column legend:\n");
      fprintf (fp,"#  [0] t\n");
      fprintf (fp,"#  [1] <Emag_p>\n");
      fprintf (fp,"#  [2] Max(Emag)\n");
      fprintf (fp,"#  [3] delta By  = Max(By) - Min(By)\n");
      fprintf (fp,"# \n");
      fprintf (fp,"# %7s %13s %13s %13s\n", 
                   "[0]", "[1]", "[2]", "[3]");
      fprintf (fp,"# -------------------------------------------------------");
      fprintf (fp,"------------------------------------------------------- \n");

    }else{

    /* Append if this is not step 0 */

      if (tpos < 0.0){ /* Obtain time coordinate of to last written row */
        char
        sline[512];
        fp = fopen(fname,"r");
        while (fgets(sline, 512, fp)) {}
        sscanf(sline, "%lf\n",&tpos); /* tpos = time of the last written row */
        fclose(fp);
      }
      fp = fopen(fname,"a");
    }

    if (g_time > tpos){     /* Write if current time if > tpos */
      fprintf (fp, "%12.6e  %12.6e  %12.6e  %12.6e\n",g_time,
                    Emag_perp_ave, Emag_max, By_max - By_min);
    }
    fclose(fp);
  }

/* --------------------------------------------------
    Write particle trajectory
   -------------------------------------------------- */
#if 0
  particleNode *cur;
  Particle *p;

  PARTICLES_LOOP(cur, d->PHead){
    p = &(cur->p);
    if (p->id == 10 && p->birth_rank == 0){
      Particles_WriteTrajectory(p, g_stepNumber == 0 ? 'w':'a');
    }
  }
#endif
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
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */
    double Bm2, vA, vAmax = 100.0;
    static Data_Arr Uc;
    RBox   domBox;

    DOM_LOOP(k,j,i){
      Bm2 =   d->Vc[BX1][k][j][i]*d->Vc[BX1][k][j][i]
            + d->Vc[BX2][k][j][i]*d->Vc[BX2][k][j][i]
            + d->Vc[BX3][k][j][i]*d->Vc[BX3][k][j][i];
      vA = sqrt(Bm2/d->Vc[RHO][k][j][i]);
      if (vA > vAmax) {
        d->Vc[RHO][k][j][i] = sqrt(Bm2)/vAmax;
      }
    }

    g_smallPressure = 1.e-3;
    DOM_LOOP(k,j,i){
      if (d->Vc[PRS][k][j][i] < g_smallPressure) {
        d->Vc[PRS][k][j][i] = g_smallPressure;
      }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif

/* ********************************************************************* */
void VectorRotate(double *v, int s)
/*!
 *  Rotate a vector <v[0], v[1], v[2]> based on the domain aspect ratio.
 *  s =  1 use normal rotation
 *  s = -1 use inverse rotation
 *
 *********************************************************************** */
{
  double vx, vy, vz;
  double Lx = g_domEnd[IDIR] - g_domBeg[IDIR];
  double Ly = g_domEnd[JDIR] - g_domBeg[JDIR];
  double Lz = g_domEnd[KDIR] - g_domBeg[KDIR];
  double ta = 0.0, tb  = 0.0, sa, sb, ca, cb, tg, cg, sg;
  double k0, kx, ky, kz, phi;

/* -- Compute tangent, sin() and cos() of alpha and beta -- */

  D_EXPAND(           ,
           ta = Lx/Ly;,
           tb = Lx/Lz;)

  if (fabs(Lx-1.0) < 1.e-9) ta = tb = 0.0; /** Special case */

//tb = 0.0;

  ca = 1.0/sqrt(ta*ta + 1.0);
  cb = 1.0/sqrt(tb*tb + 1.0);
  sa = ta*ca;
  sb = tb*cb;

  tg = tb*ca;
  cg = 1.0/sqrt(tg*tg + 1.0);
  sg = cg*tg;

/* -- Save initial components -- */

  vx = v[0]; 
  vy = v[1];
  vz = v[2];

/* -- Rotate -- */

  if (s == 1){
    v[0] = vx*ca*cg - vy*sa - vz*ca*sg;
    v[1] = vx*sa*cg + vy*ca - vz*sa*sg;
    v[2] = vx*sg            + vz*cg;
  }else if (s == -1){
    v[0] =  vx*ca*cg + vy*sa*cg + vz*sg;
    v[1] = -vx*sa    + vy*ca;
    v[2] = -vx*ca*sg - vy*sa*sg + vz*cg;
  }else{
    print ("! VectorRotate(): invalid value of s = %d\n", s);
    QUIT_PLUTO(1);
  }

}
