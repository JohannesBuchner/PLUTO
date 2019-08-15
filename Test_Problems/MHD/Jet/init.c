/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Axisymmetric jet propagation.
 
  \b Description.\n
  The jet problem is set up in axisymmetric cylindrical coordinates 
  \f$ (R,z) \f$.
  We label ambient and jet values with the suffix "a" and "j", respectively.
  The ambient medium is at rest and is characterized by uniform density 
  and pressure, \f$ \rho_a \f$ and \f$ p_a \f$.
  The beam enters from the lower-z boundary from a circular nozzle of 
  radius \f$ R_j \f$ and carries a constant poloidal field \f$ B_z \f$
  and a (radially-varying) toroidal component \f$ B_\phi(R) \f$.
  Flow variables are prescribed as 
  \f[ \left\{\begin{array}{lcl}
     \rho(R)   &=& \rho_j  \\ \noalign{\medskip}
     v_z(R)    &=& v_j     \\ \noalign{\medskip}
     B_\phi(R) &=& \DS \left\{\begin{array}{ll}
                        -B_m R/a & \quad{\rm for}\quad R<a \\ \noalign{\medskip}
                        -B_m a/R & \quad{\rm for}\quad R>a 
                       \end{array}\right. \\ \noalign{\medskip}
     B_z(R)    &=& B_{z0}\quad\mathrm{(const)}
  \end{array}\right. \f]
  Here \c a is the magnetization radius while \f$ v_R = B_R = 0 \f$.
  These profiles are similar to the ones used by Tesileanu et al. (2008).
  The pressure distribution is found by solving the radial momentum 
  balance between thermal, centrifugal and magnetic forces:
  \f[
    \frac{dp}{dR} = \frac{\rho v_\phi^2}{R} -
                    \frac{1}{2}\left[\frac{1}{R^2}\frac{d(R^2B^2_\phi)}{dR} 
                                    +\frac{dB_z^2}{dR}\right]
  \f]
  Neglecting rotation and assuming \c Bz to be constant the solution of 
  radial momentum balance becomes
  \f[
    p(R) = p_a + B_m^2\left[1-\min\left(\frac{R^2}{a^2},1\right)\right]
  \f]
  and the jet on-axis pressure increases for increasing toroidal field:
  \f[
    p(R=0) \equiv p_j = p_a + B_m^2
  \f]
  where \f$p_a\f$ is the ambient pressure.

  \b Normalization.\n
  Since the MHD equations are scale-invariant we have the freedom to 
  specify a reference length, density and velocity. Here we choose
  - Length: jet radius \f$R_j=1\f$;
  - Density: ambient density \f$\rho_a=1\f$; 
  - Velocity:
    - adiabtic setup (no cooling):
      ambient sound speed: \f$c_a = \Gamma p_a/\rho_a = 1\f$
      (from which it follows \f$p_a = 1/\Gamma\f$).
    - radiative setup (with cooling): 1 Km/s (set by \c UNIT_VELOCITY).
      In this case, the ambient pressure is computed from the ambient 
      temperature <tt> Ta = 2500 K </tt>.
       
  In this way the number of parameters is reduced to 4:
  -# <tt>g_inputParam[ETA]</tt>: density contrast 
     \f$ \eta = \rho_j/\rho_a
         \qquad\Longrightarrow\qquad \rho_j = \eta
     \f$
  -# <tt>g_inputParam[JET_VEL]</tt>:  Jet velocity  \f$ v_j \f$.
     This is also the Mach number of the beam with respect to the ambient
     in the adiabatic setup;
  -# <tt>g_inpurParam[SIGMA_Z]</tt>: poloidal magnetization;
  -# <tt>g_inpurParam[SIGMA_PHI]</tt>: toroidal magnetization;
  Magnetization are defined as follows:
  \f[ \left\{\begin{array}{ll}
       \sigma_z    &= \DS \frac{B_z^2}{2p_a}      \\ \noalign{\medskip}
       \sigma_\phi &= \DS \frac{<B_\phi^2>}{2p_a} 
     \end{array}\right.
   \qquad\Longrightarrow\qquad
     \left\{\begin{array}{ll}
        B_z^2 &= \DS 2\sigma_zp_a       \\ \noalign{\medskip}
        B_m^2 &= \DS \frac{2\sigma_\phi}{a^2(1/2-2\log a)}p_a
      \end{array}\right.    
  \f]
  Here the average value of \f$B^2_\phi\f$ is simply found as
  \f[
     <B^2_\phi> = \frac{\int_0^1 B^2_\phi R\,dR}{\int_0^1 R\,dR}
                = B_m^2a^2\left(\frac{1}{2} - 2\log a\right)
  \f]
  
  The following MAPLE code can be used to verify the solution:
  \code
   restart;   # Solve radial equilibrium balance equation
   assume (a < 1, a > 0);
   B[phi] := -B[m]*R/a;
   ode    := diff(p(R),R) = -1/(2*R^2)*diff(R^2*B[phi]^2,R);
   dsolve({ode,p(a)=pa});   # only for R < a
  \endcode

  When cooling is enabled, two additional parameters controlling the
  amplitude and frequency of perturbation can be used:
  -# <tt>g_inputParam[PERT_AMPLITUDE]</tt>: perturbation amplitude (in 
     units of jet velocity);
  -# <tt>g_inputParam[PERT_PERIOD]</tt>: perturbation period (in years).

  The jet problem has been tested using a variety of configurations, namely:
  <CENTER>
  Conf.| GEOMETRY  |DIM| T. STEPPING|RECONSTRUCTION|divB| EOS    | COOLING
  -----|-----------|---|----------- |--------------| ---|--------|--------
   #01 |CYLINDRICAL| 2 |  RK2       |  LINEAR      | CT | IDEAL  |   NO
   #02 |CYLINDRICAL| 2 |  HANCOCK   |  LINEAR      | CT | IDEAL  |   NO 
   #03 |CYLINDRICAL| 2 |  ChTr      |  PARABOLIC   | CT | IDEAL  |   NO 
   #04 |CYLINDRICAL|2.5|  RK2       |  LINEAR      |NONE| IDEAL  |   NO
   #05 |CYLINDRICAL|2.5|  RK2       |  LINEAR      | CT | IDEAL  |   NO 
   #06 |CYLINDRICAL|2.5|  RK2       |  LINEAR      | CT |PVTE_LAW|   NO 
   #07 |CYLINDRICAL|2.5|  ChTr      |  PARABOLIC   |NONE| IDEAL  |   SNEq    
   #08 |CYLINDRICAL|2.5|  ChTr      |  PARABOLIC   |NONE| IDEAL  |   MINEq
   #09 |CYLINDRICAL|2.5|  ChTr      |  PARABOLIC   |NONE| IDEAL  |   H2_COOL
  </CENTER>
  Here 2.5 is a short-cut for to 2 dimensions and 3 components.

  The following image show the (log) density map at the end of 
  simulation for setup #01.

  \image html mhd_jet.01.jpg "Density map at the end of the computation for configuration #01"

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 02, 2017

  \b References
     - "Simulating radiative astrophysical flows with the PLUTO code:
       a non-equilibrium, multi-species cooling function"
       Tesileanu, Mignone and Massaglia, A&A (2008) 488,429.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h" 

void GetJetValues    (double, double *);
static double Profile (double, double);

static double a = 0.8;      /* -- Magnetization radius -- */
static double Ta, pa,mua;   /* -- Ambient pressure -- */

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 * Define ambient values and units.
 *********************************************************************** */
{
  int    nv;
  static int first_call = 1;
  double nH = 200.0;
  static double veq[NVAR];
  
  v[RHO] = 1.0;
  v[VX1] = v[VX2] = v[VX3] = 0.0;
#if COOLING == NO
   v[PRS] = pa = 0.6;
   /* v[PRS] = Pressure(v, Ta); */
#else
  #if COOLING == H2_COOL
   Ta = 100.0;  /* Ambient temperature for molecular cooling  */
  #else
   Ta = 2500.0; /* Ambient temperature for atomic cooling  */
  #endif
  if (first_call) CompEquil(nH, Ta, veq);
  for (nv = NFLX; nv < NFLX + NIONS; nv++) v[nv] = veq[nv];
  mua    = MeanMolecularWeight(v);
  v[PRS] = pa = Ta*v[RHO]/(KELVIN*mua);
#endif

#if PHYSICS == MHD
  EXPAND(v[BX1] = 0.0;                                ,   
         v[BX2] = sqrt(2.0*g_inputParam[SIGMA_Z]*pa); ,
         v[BX3] = 0.0;)

  v[AX1] = v[AX2] = 0.0;
  v[AX3] = 0.5*x1*v[BX2];
#endif

  v[TRC] = 0.0;

#if COOLING == H2_COOL
  g_minCoolingTemp = 10.0;
  g_maxCoolingRate = 0.1;
#else
  g_minCoolingTemp = 2500.0;
  g_maxCoolingRate = 0.2;
#endif
  g_smallPressure  = 1.e-3;

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
/* 
 *
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 * Set the injection boundary condition at the lower z-boundary
 * (\c X2-beg must be set to \c userdef in \c pluto.ini).
 * For \f$ R \le 1 \f$ we set constant input values (given by the
 * ::GetJetValues() function while for $ R > 1 $ the solution has
 * equatorial symmetry with respect to the \c z=0 plane.
 * To avoid numerical problems with a "top-hat" discontinuous jump,
 * we smoothly merge the inlet and reflected value using a profile 
 * function ::Profile().
 *
 *********************************************************************** */
{
  int    i, j, k, nv;
  double *x1, *x2, *x3;
  double r, vjet[256], vout[NVAR], q;
  double t, omega; /* pulsed jet frequency */

  x1  = grid->xgc[IDIR];
  x2  = grid->xgc[JDIR];
  x3  = grid->xgc[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){}
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER){    /* -- cell-centered boundary conditions -- */
      BOX_LOOP(box,k,j,i){
        GetJetValues (x1[i], vjet);

      /* -- add pulsed perturbation -- */

        #if COOLING != NO
         t = g_time*UNIT_LENGTH/UNIT_VELOCITY;  /* time in seconds */
         omega = 2.0*CONST_PI/(86400.*365.*g_inputParam[PERT_PERIOD]);
         vjet[VX2] *= 1.0 + g_inputParam[PERT_AMPLITUDE]*sin(omega*t);
        #endif

      /* -- copy and reflect ambient medium values -- */

        VAR_LOOP(nv) vout[nv] = d->Vc[nv][k][2*JBEG - j - 1][i];
        vout[VX2] *= -1.0;
        #if PHYSICS == MHD
        EXPAND(vout[BX1] *= -1.0; ,  
                                ; , 
               vout[BX3] *= -1.0;)
        #endif
        VAR_LOOP(nv){
          d->Vc[nv][k][j][i] = vout[nv]-(vout[nv]-vjet[nv])*Profile(x1[i], nv);
        }
      }

    }else if (box->vpos == X1FACE){  /* -- staggered fields -- */
      #ifdef STAGGERED_MHD
       x1 = grid->xr[IDIR];
       BOX_LOOP(box,k,j,i){
         vout[BX1] = -d->Vs[BX1s][k][2*JBEG - j - 1][i];
         d->Vs[BX1s][k][j][i] =    vout[BX1] 
                                - (vout[BX1] - vjet[BX1])*Profile(fabs(x1[i]), BX1);
       }
      #endif
    }else if (box->vpos == X3FACE){

    }
  }
}

/* ********************************************************************* */
void GetJetValues (double R, double *vj)
/*! 
 * Define jet value as functions of the (cylindrical) radial
 * coordinate.
 * For a periodic jet configuration, density and verical velocity
 * smoothly join the ambient values. This has no effect on equilibrium
 * as long as rotation profile depends on the chosen density profile.
 * 
 *********************************************************************** */
{
  int nv;
  static int first_call = 1;
  static double veq[NVAR];
  double  Bz, Bm = 0.0, x;

  if (fabs(R) < 1.e-9) R = 1.e-9;
  x  = R/a;

  vj[RHO] = 1.0 + (g_inputParam[ETA] - 1.0);
  Bz = sqrt(2.0*g_inputParam[SIGMA_Z]*pa);
  Bm = sqrt(2.0*g_inputParam[SIGMA_PHI]*pa/(a*a*(0.5 - 2.0*log(a))));

#if PHYSICS == MHD
  EXPAND( vj[BX1] =  0.0;                      ,
          vj[BX2] =  Bz;                       ,
          vj[BX3] = -Bm*(x < 1.0 ? x:1.0/x);)
  vj[AX1] = vj[AX2] = 0.0;
  vj[AX3] = 0.5*R*Bz;
#endif

  EXPAND( vj[VX1] = 0.0;                     ,
          vj[VX2] = g_inputParam[JET_VEL];   ,
          vj[VX3] = 0.0; )
  
  vj[PRS] = pa + Bm*Bm*(1.0 - MIN(x*x,1.0));
  vj[TRC] = 1.0;

  #if COOLING != NO
   static double Tj,muj;
   #if COOLING == H2_COOL
    Tj = Ta*(1.0 + Bm*Bm/pa)/g_inputParam[ETA]; /* Guess Tj assuming muj = mua */
    if (first_call){
      CompEquil(200.0*g_inputParam[ETA], Tj, veq);
      muj = MeanMolecularWeight(veq);
      Tj *= muj/mua;   /* Improve guess */
    }
   #else
    if (first_call) CompEquil(200.0*g_inputParam[ETA], Ta, veq);
   #endif
    for (nv = NFLX; nv < NFLX+NIONS; nv++) vj[nv] = veq[nv];
  #endif
  first_call = 0;

}

/* ********************************************************************* */
double Profile (double R, double nv)
/*
 *
 *
 *********************************************************************** */
{
  double R4 = R*R*R*R, R8 = R4*R4;

#if PHYSICS == MHD && COMPONENTS == 3    /* Start with a smoother profile */
  if (g_time < 0.1)  return 1.0/cosh(R4); /* with toroidal magnetic fields */
#endif
  return 1.0/cosh(R8);
}
