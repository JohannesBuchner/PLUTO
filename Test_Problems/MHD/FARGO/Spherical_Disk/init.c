/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Accretion disk in 3D spherical coordinates.

  Setup a magnetized disk in 3D spherical coordinates \f$(r,\theta,\phi\f$).
  The disk model is taken form [Flo11],[Mig12] and it consists of an initial
  equilibrium in a point mass gravity \f$ 1/r^2\f$ with equilibrium profiles
  given by
  \f[
     \rho = \frac{1}{R^{3/2}}\exp\left[\frac{\sin(\theta) - 1}{c_0^2}\right]
      \,,\qquad
     p    = c_s^2\rho\
      ,,\qquad
     v_\phi = \frac{1}{\sqrt{r}}\left(1 - \frac{5}{2\sin\theta}c_0^2\right)
  \f]
  where \f$R = r\sin\theta\f$ is the cylindrical radius while
  \f$c_0 = H/R\f$ and \f$c_s = H/R^{3/2}\f$ is the sound speed.
  Here \c H is the scale height and it is proportional to the cylindrical
  radius.
  The constant of proportionality is the ratio \f$ H/R \f$ defined by the
  input parameter <tt>g_inputParam[H_R]</tt>.
  With the current setting, one rotation period of the inner orbit (\c r=1)
  is \f$\Delta T = 2\pi\f$.

  Differently from [Flo11], a polodial magnetic field is used here, with
  vector potential prescribed as follows:
  \f[
     A_\phi = A_0\frac{\sin(2\pi R) - R\cos(2\pi R)}{R}
              \exp\left[-\left(\frac{z}{H}\right)^4\right]
              \exp\left[-\left(\frac{R-6}{2}\right)^4\right]
  \f]
  where \c A0 is a constant chosen to prescribe a given plasma beta.
  The exponential terms on the right hand side confine the magnetic field
  in the midplane around \c R=6.
  We also make the vector potential vanish identically for \f$ z/H > 2.5\f$
  or \f$ R < 2.5\f$.

  The boundary conditions are all set to \c userdef:

  - at the inner radial boundary, zero-gradient is imposed except for
    the azimuthal velocity (which is fixed the equilibrium value) and the
    radial velocity which is zero if velocity is directed inside the domain
    or copied otherwise (this is referred to as the "diode" b.c.)
  - at the outer radial boundary we impose zero-gradient except for the
    azimuthal velocity (fixed to equil. value) and radial velocity
    which is reflected.
  - at the vertical (theta) boundaries, ghost zones are filled with
    the initial equilibrium solution.
   
  For the sake of simplicity only a quarter of a disk is considered.
  This test problems can be run normally or using the FARGO module,
  giving a saving factor of (roughly) 5.
  The figure below was produced by running the current setup at twice the
  resolution and for 90 orbits using FARGO.

  \image html mhd_spherical_disk.02.png "Density map of the disk after 90 innter orbits at twice the resolution."
  \authors A. Mignone (mignone@ph.unito.it)\n
  \date   Aug 16, 2012
  
  \b References
     - [Flo11]: "Turbulence and Steady Flows in 3D Global Stratified
       MHD Simulations of Accretion Disks", Flock et al., ApJ (2011), 735:122
     - [Mig12]: "A Conservative orbital advection scheme for simulations of
       magnetized shear flows with the PLUTO code",
       Mignone et al, A&A (2012) 545, A152
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  static int first_call = 1;
  double rnd, s, c, R, r, z;
  double H, c0, cs, Bphi;

  g_gamma = 1.0001;
  r   = x1;
  s   = sin(x2);  
  c   = cos(x2);

  R  = r*s;   /* Cylindrical radius */ 
  z  = r*c;
  c0 = g_inputParam[H_R];
  cs = c0/sqrt(R);
  H  = g_inputParam[H_R]*R;

  us[RHO]   = exp((s-1.0)/(c0*c0))/(R*sqrt(R));
  us[iVR]   = 0.0;
  us[iVTH]  = 1.e-4*cs*cos(8.0*x3)*exp(-2.0*z*z/(H*H));  /* Perturbation */
  us[iVPHI] = 1.0/sqrt(r)*(1.0 - 2.5*c0*c0/s);
  us[PRS]   = cs*cs*us[RHO];

  #if PHYSICS == MHD
   us[BX1] = us[BX2] = us[BX3] = 0.0;
   us[AX1] = us[AX2] = us[AX3] = 0.0;

   us[AX3]  = (sin(2.0*CONST_PI*R) - R*cos(2.0*CONST_PI*R))/R;
   us[AX3] *= 3.5e-4*exp(-pow(z/H,4.0))*exp(-pow((R-6.0)/2.0,4.0));
   if (fabs(z/H) > 2.5 || R < 2.5) us[AX3] = 0.0;
  #endif
  g_smallPressure = 1.e-9;
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
/* 
 *
 *
 *********************************************************************** */
{
  int    i,j,k,nv;
  static int first_call = 1;
  double v[NVAR];
  double *r, *th, vol;
  double pm, kin;
  double pmtot, Tmtot;
  static double voltot;
  FILE *fp;

  r   = grid->x[IDIR];
  th  = grid->x[JDIR];

/* ---------------------------------------------------------
         compute total volume once at the beginning
   --------------------------------------------------------- */

  if (first_call){
    voltot = 0.0;
    DOM_LOOP(k,j,i) voltot += grid->dV[k][j][i];
    #ifdef PARALLEL
     MPI_Allreduce (&voltot, &kin, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
     MPI_Barrier (MPI_COMM_WORLD);
     voltot = kin;
    #endif
  }

/* --------------------------------------------------------
              compute volume averages 
   -------------------------------------------------------- */

  pmtot = Tmtot = 0.0;
  DOM_LOOP(k,j,i){
    vol = grid->dV[k][j][i]/voltot;

    for (nv = NVAR; nv--;   ) v[nv] = d->Vc[nv][k][j][i];

    #if PHYSICS == MHD
     pm     = 0.5*(v[BX1]*v[BX1] + v[BX2]*v[BX2] + v[BX3]*v[BX3]);
     pmtot += pm*vol;
     Tmtot += v[BX1]*v[BX3]*vol;
    #endif
  }
 
  #ifdef PARALLEL
   MPI_Allreduce (&pmtot, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   pmtot = vol;

   MPI_Allreduce (&Tmtot, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   Tmtot = vol;

   MPI_Barrier (MPI_COMM_WORLD);
  #endif

/* -----------------------------------------------------
    Write "averages.dat" to disk. 
    The file is created as new if this is the very
    initial time step. Data is appended otherwise.
    Processor #0 does the actual writing.
   ----------------------------------------------------- */

  if (prank == 0){
    static double tpos = -1.0;
    if (g_stepNumber == 0){
      fp = fopen("averages.dat","w");
      fprintf (fp,"#%s %12s  %12s  %12s\n", 
               "      t","dt","<B^2>/2","<Tm>");
    }else{
      if (tpos < 0.0){  /* obtain time coordinate of to last entry */
        char   sline[512];
        fp = fopen("averages.dat","r");
        while (fgets(sline, 512, fp))  {}
        sscanf(sline, "%lf\n",&tpos);
        fclose(fp);
      }
      fp = fopen("averages.dat","a");
    }
    if (g_time > tpos){
      fprintf (fp, "%12.6e  %12.6e  %12.6e  %12.6e\n",
               g_time, g_dt, pmtot, Tmtot);
    }
    fclose(fp);
  }
  first_call = 0; 
}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 *
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double c0, cs, qx, qy, qz, s, v[256];
  double *r, *th, *phi;
 
  r   = grid->xgc[IDIR];
  th  = grid->xgc[JDIR];
  phi = grid->xgc[KDIR];
  c0  = g_inputParam[H_R];
  
  if (side == 0) {    /* -- check solution inside domain -- */
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i]= d->Vc[nv][k][j][IBEG];
        s = sin(th[j]);
        d->Vc[iVPHI][k][j][i] = 1.0/sqrt(r[i])*(1.0 - 2.5*c0*c0/s);
        if (d->Vc[iVR][k][j][i] > 0.0) d->Vc[iVR][k][j][i] = 0.0;
      }
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IBEG];  
      #endif
    }else if (box->vpos == X3FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX3s][k][j][i] = d->Vs[BX3s][k][j][IBEG];  
      #endif
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER){
      BOX_LOOP(box, k, j, i){
        for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i]= d->Vc[nv][k][j][IEND];
        s = sin(th[j]);
        d->Vc[iVPHI][k][j][i] = 1.0/sqrt(r[i])*(1.0 - 2.5*c0*c0/s);
        d->Vc[iVR][k][j][i]   = -d->Vc[iVR][k][j][2*IEND - i + 1];
      }
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i)  d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IEND];
      #endif
    }else if (box->vpos == X3FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i)  d->Vs[BX3s][k][j][i] = d->Vs[BX3s][k][j][IEND];
      #endif
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER){
      BOX_LOOP(box, k, j, i){
        Init (v, r[i], th[j], phi[k]);
        cs      = c0/sqrt(r[i]*sin(th[j]));
        v[iVTH] = 0.0;
        v[PRS]  = v[RHO]*cs*cs; 
        for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i]= v[nv];
      }
    }else if (box->vpos == X1FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box, k, j, i) d->Vs[BX1s][k][j][i] = 0.0;
      #endif
   }else if (box->vpos == X3FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box, k, j, i) d->Vs[BX3s][k][j][i] = 0.0;
      #endif
    }

  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    if (box->vpos == CENTER){
      BOX_LOOP(box, k, j, i){
        Init (v, r[i], th[j], phi[k]);
        cs      = c0/sqrt(r[i]*sin(th[j]));
        v[iVTH] = 0.0;
        v[PRS]  = v[RHO]*cs*cs; 
        for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i]= v[nv];
      }
    }else if (box->vpos == X1FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box, k, j, i) d->Vs[BX1s][k][j][i] = 0.0;
      #endif
    }else if (box->vpos == X3FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box, k, j, i) d->Vs[BX3s][k][j][i] = 0.0;
      #endif
    }
  }
}

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  g[IDIR] = -1.0/(x1*x1);
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  return 0.0;
}
#endif
