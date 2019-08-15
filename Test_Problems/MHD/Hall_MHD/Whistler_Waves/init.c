/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Whistler Waves.


  The dispersion relation for whistler waves is given by (Huba, 2003):
  \f[
   \omega = \frac{k^2_x B}{e n_e}
  \f]
  where\f$k_x = 2\pi m/L_x\f$ is the wave number along the \c x direction,
  \f$ L_x \f$ is the domain lenght in the horizontal direction,
  \f$ e \f$ and \f$ n_e \f$ are the electron charge and number
  density \c B is the magnetic field strength.

  <B> Whistler waves 1 and 2D (\c SETUP \c == \c 1)</B>

  In configurations #01 and #02 we test the propagation of whistler waves
  in 1D and 2D, respectively, and the results can be compared with 
  the analytical dispersion relation, in a manner similar to Huba 2003.
  In 2D, the setup is rotated by an amount proportional to the ratio
  between the vertical and horizontal axis extents.
  The magnetic field is in the x-direction,  \f$\vec{B} = B_0 \hvec{e}_x \f$,
  and the system is perturbed with
  \f$ \delta B_y = +\delta B \cos(2\pi m x/Lx) \f$ and
  \f$ \delta B_z = -\delta B \sin(2\pi m x/Lx) \f$, where
  \f$ m = 1, 2, ..., \f$ is the mode number, \f$ B_0 = 100 \f$ G and
  \f$ \delta B = 10^{-3} \f$ G.
  We also set \f$ \rho = en_e = 1 \f$.

  We show in Fig. 1 the temporal evolution of the perpendicular component
  of the magnetic field (\f$B_y\f$) and the its fourier composition.
  In Fig. 2 we show the comparison between the analytical and numerical
  values of the whistler wave frequency as a function of \f$ m \f$.
  <br>
  \image html whistler_waves_fig1.png "Fig. 1 - Temporal evolution of the perpendicular component of the magnetic field (B_y) at x_0 = 0.5 (left), and its fourier decomposition (real amplitude) as a function of frequency for the run with mode number m = 6."
  <br>
  \image html whistler_waves_fig2.png "Fig. 2 - Pseudocolor rendering of the current density (color) with magnetic field lines at  t = 57 for the Hall MHD simulation of the Whistler waves."
  <br>

  <B> Whistler waves 2D (\c SETUP \c == \c 2)</B>

  In this section (configuration #02) we test the correct propagation of
  the whistler waves in 2D, in a manner similar to Viganò et al., 2012. 
  The initial magnetic field is \f$B_x = B_0 + B_1 \cos(k_x y) \cos(k_x x)\f$,
  \f$B_y = B_1 \sin(k_x y) \sin(k_x x)\f$ and
  \f$ B_z = \sqrt{2} B_1 \sin(k_x y) \cos(k_x x)\f$,
  where \f$ B_0 = 10^3 \f$ and \f$B_1 = 1\f$.
  Again we set \f$ \rho = en_e = 1 \f$.

  The equations above admit wave solution travelling along \f$ x \f$
  with speed 
  \f[
  v_w = -\frac{c}{4 e n_0} \frac{\sqrt{2} }{L_y} m B_0
  \f] 
  where \f$ L_y\f$ is the domain length in the vertical direction,
  and again \f$ m \f$ is the number of wave modes.

  In Fig. 3 we show the result of our simulation (\f$ B_y \f$) at
  \f$ t = 0 \f$. We checked the scaling of the whistler speed with
  \f$ m \f$, shown in Fig. 4.

  \image html whistler_waves_fig3.png "Fig. 3 - Pseudocolor rendering of B_y at  t = 0  for the simulation of the 2D whistler waves defined by the above equations  (positive values in black and negative values in red)."
<br>
  \image html whistler_waves_fig4.png "Fig. 4 - The whistler speed for different values of m (red circles) compared with the analitical values from the expression of v_w"


  \author E. Striani (edoardo.striani@iaps.inaf.it) \n
          A. Mignone (mignone@ph.unito.it) \n
          B. Vaidya  (bhargav.vaidya@unito.it) 
  \date   May 17, 2017 

  \b References: \n
     - Huba 2003, "Hall Magnetohydrodynamics - A Tutorial",
       Huba, J.D.\ 2003, Space Plasma Simulation, 615, 166 
     - Viganò et al., 2012, "A new code for the Hall-driven magnetic evolution of neutron stars",
       Viganò, D., Pons, J.A., \& Miralles, J.A.\ 2012, Computer Physics Communications, 183, 2042 

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
void VectorRotate(double *v, int s);

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *********************************************************************** */
{ 
  int m;
  double Lx, Ly, B0, dB, k;

  m  = g_inputParam[PMODE];
  Lx = g_domEnd[IDIR] - g_domBeg[IDIR];
  Ly = g_domEnd[JDIR] - g_domBeg[JDIR];
  k  = 2.0*CONST_PI*m/Lx;
  dB = 1.0;
  B0 = 1.0e3;

  v[RHO] = 1.0;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

#if SETUP == 1  
  double bref;
  dB *= 1.e-3; /* field in Gauss */
  B0 *= 0.1;   /* field in Gauss */
  bref = sqrt(4.0*CONST_PI*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);

  v[BX1] = B0/bref;
  v[BX2] =  (dB/bref)*cos(k*x1);
  v[BX3] = -(dB/bref)*sin(k*x1);

/* ------------------------------------
   In 2D, rotate vectors.
   ------------------------------------ */

  #if DIMENSIONS > 1
  m  = g_inputParam[PMODE];
  double kx  = 2.0*CONST_PI*m/Lx;
  double ky  = 2.0*CONST_PI/Ly;
  double phi = kx*x1 + ky*x2;

  v[BX1] = B0/bref;
  v[BX2] =  (dB/bref)*cos(phi);
  v[BX3] = -(dB/bref)*sin(phi);
  VectorRotate(v+BX1,1);
  #endif

#elif SETUP == 2
   v[BX1] = B0 + dB*cos(k*x2)*cos(k*x1);
   v[BX2] = dB*sin(k*x2)*sin(k*x1);
   v[BX3] = sqrt(2.0)*dB*sin(k*x2)*cos(k*x1);
#endif

}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/* *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* ********************************************************************** */
{
  char fname[32];
  long i;
  double Bz;
  FILE *fp;

  if (prank != 0) return;  /* Only proc. 0 does the job */

/* ---------------------------------------------------
   1. Open files for writing/appending:
      - "whistler_waves.m.dat" contains Bz(x0) as a
        function of time;
      - "whistler_waves_max.m.dat" contains the maxima
        of Bz(x0) and the computed frequency;
   --------------------------------------------------- */

  sprintf (fname,"whistler_waves.%d.dat",(int)g_inputParam[PMODE]);
  if (g_stepNumber == 0) {
    fp = fopen (fname,"w");
    fprintf (fp,"# --------------------------------------------\n");
    fprintf (fp,"#  File: %s\n",fname);
    fprintf (fp,"#  Column legend:\n");
    fprintf (fp,"#  [0] = t\n");
    fprintf (fp,"#  [1] = Bz(x0)\n");
    fprintf (fp,"# --------------------------------------------\n");
  } else fp = fopen (fname,"a");


/* ----------------------------------------------------
   2. Write Bz at a given point.
   ---------------------------------------------------- */

  Bz = d->Vc[BX3][KBEG][JBEG][IBEG];
  fprintf (fp, "%13.6e  %13.6e\n", g_time, Bz);
  
/* -----------------------------------------------
   3. Close file
   ----------------------------------------------- */

  fclose (fp);

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/* ********************************************************************** */
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

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/* ********************************************************************* */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/* ********************************************************************** */
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

