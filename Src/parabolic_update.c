/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Update solution array with parabolic (diffusion) terms.

  The ParabolicUpdate() function updates a solution array of conservative
  variables using contributions coming from diffusion operators only.
  It is called by EXPLICIT methods only immediately after the hyerbolic
  flux contributions have been added to \c U. 
  Note that this is \e not an operator split formalism since contributions
  are added to the base time level:
  \f[
    \begin{array}{llcll}
      (1)\quad & \vec{U}^{n} &\to& \vec{U}^*     & = \vec{U}^n
                              - \Delta t\nabla\cdot\vec{F}_h^n
      \\ \noalign{\medskip}
      (2)\quad & \vec{U}^{*} &\to& \vec{U}^{n+1} & = \vec{U}^*
                                + \Delta t\nabla\cdot\vec{F}_p^n
    \end{array}
  \f]
  Step (2) is the one performed by ParabolicUpdate().
  By convention, parabolic fluxes are written with the plus sign
  when they are on the right hand side.

  The ParabolicRHS() function computes the right hand side, in 
  divergence form, of the parabolic (diffusion) operators only:
  \f[
    \pd{U}{t} = \nabla\cdot\Big(D\nabla U\Big) + S
  \f]  
  Here \c U is a generic cell-centered quantity (like momentum, magnetic
  field or total energy) and \f$\vec{F} = D\nabla U\f$ is the corresponding
  flux and \c S is the source terms (including geometrical source terms).
  Contributions may simultaneously come from ambipolar diffusion, resistivity,
  thermal conduction and viscosity computed, respectively, by the functions
  AD_Flux(), ResistiveFlux(), TC_Flux() and ViscousFlux().
 
  For resistivity / ambipolar diffusion / Hall MHD, we use the divergence
  form of the induction equation when updating cell-centered magnetic fields:
  \f[
    \pd{\vec{B}}{t} = - \nabla\times\vec{E} \equiv -\nabla\cdot\tens{M}
    \qquad\mathrm{where}\qquad
    \tens{M} = \left(\begin{array}{ccc}
        0    &  E_z  & -E_y  \\ \noalign{\medskip}
       -E_z  &  0    &  E_x  \\ \noalign{\medskip}
        E_y  & -E_x  & 0
    \end{array}\right) \qquad\mathrm{or\quad also}\qquad
    \tens{M}_{ij} = (\vec{E}\times\hvec{e}_i)\cdot\hvec{e}_j
  \f]
  which is valid for any electric vector while \f$\tens{M}\f$ is an
  antisymmetric tensor.
  Note that in curvilinear geometries, the right hand side contains
  source terms.
  In cylindrical coordinates, for instance:
  \f[
     \begin{array}{lcl}
     (\nabla\cdot\tens{M})_r &=&\DS
                                   \frac{1}{r}\pd{}{r}(rM_{rr})
                                 + \frac{1}{r}\pd{M_{\phi r}}{\phi}
                                 +            \pd{M_{zr}}{z}
                                 - \frac{M_{\phi\phi}}{r}
     \\ \noalign{\medskip}
     (\nabla\cdot\tens{M})_\phi &=&\DS
                                   \frac{1}{r}\pd{}{r}(rM_{r\phi})
                                 + \frac{1}{r}\pd{M_{\phi\phi}}{\phi}
                                 +            \pd{M_{z\phi}}{z}
                                 + \frac{M_{\phi r}}{r}
     \\ \noalign{\medskip}
     (\nabla\cdot\tens{M})_z &=&\DS
                                   \frac{1}{r}\pd{}{r}(rM_{rz})
                                 + \frac{1}{r}\pd{M_{\phi z}}{\phi}
                                 +            \pd{M_{zz}}{z}
     \end{array}
  \f]

  This function is intended for operator split algorithms (STS, RKL) as well
  as explicit time stepping.
  
  This function also computes the inverse diffusion time step by adding, 
  for each diffusion equation, contributions coming from different directions:
  \f[ \Delta t_{p}^{-1} = 
      \max\left[\frac{D_{i+\HALF} + D_{i-\HALF}}{2\Delta x^2} +
                \frac{D_{j+\HALF} + D_{j-\HALF}}{2\Delta y^2} +
                \frac{D_{k+\HALF} + D_{k-\HALF}}{2\Delta z^2}\right]
  \f]
  where \f$\eta_{x,y,z}\f$ are the diffusion coefficients available at 
  cell interfaces in the three directions, respectively, and the maximum
  is taken over the local processor grid.
  
       
  \date    May 13, 2018
  \authors A. Mignone (mignone@ph.unito.it)\n
           B. Vaidya
           Z. Ahmane
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"
#if PARABOLIC_FLUX != NO

#define MAX_OP   8   /* Maximum number of diffusion operators */

/* Define diffusion operator labels, in increasing order */
enum PARABOLIC_OPERATORS{
  AMB_DIFF_OP,      /* AMBIPOLAR DIFFUSION */
  HALL_OP,          /* HALL_MHD            */
  RES_OP,           /* RESISITIVITY  (3 operators since it's a tensor) */
  TC_OP = RES_OP+3, /* THERMAL CONDUCTION */
  VISC_OP,          /* VISCOSITY */
};

/* ********************************************************************* */
void ParabolicUpdate(const Data *d, Data_Arr UU, RBox *domBox, double **aflux,
                     double dt, timeStep *Dts, Grid *grid)
/*!
 * Take one step in updating the conservative vector \c UU:
 * \f[
 *      U^1 = U^0 + \Delta t R(V^0)
 * \f]
 * where \f$ U^0 =\f$ \c UU is a conservative array while
 * \f$ V^0 = \f$ \c d->Vc is the array of primitive variabels.
 * Note that \c d->Vc and \c UU may not necessarily be the map of each
 * other.
 * It is called only by explicit schemes.
 *
 *
 * \note When the entropy switch is enabled, we actually do not compute the
 *       right hand side of the entropy equation but, instead, we recover
 *       entropy variation from total energy:
 *       \f[ \begin{array}{lcl}
 *             E^*     &=& E(S^*)  \\ \noalign{\medskip}
 *             E^{n+1} &=& E^* + \Delta t (\nabla\cdot\vec{F}_{E,p}) \\ \noalign{\medskip}
 *             S^{n+1} &=& S(E^{n+1})
 *           \end{array}             
 *       \f]
 *       where \f$S^*\f$ is the value of entropy after the hyperbolic update, 
 *       \f$ \vec{F}_{E,p} \f$ is the energy diffusion flux and starred
 *       quantities are coming from the previous step.
 *
 * \param [in]     d        Pointer to the PLUTO data structure.
 *                          When set to \c NULL, the right hand side is not
 *                          recomputed, and the most recent computed value
 *                          of rhs[] is employed (useful for CTU algorithms).
 * \param [in,out] UU       Array of conservative variables to be updated
 * \param [in]     domBox   A pointer to an RBox structure defining
 *                          the zones of the domain to be updated
 * \param [in,out] aflux    A 2D pointer to store fluxes (needed by Chombo)
 * \param [in]     dt       The time step increment
 * \param [out]    Dts      Pointer to the timeStep structure
 * \param [in]     grid     Pointer to the grid structure
 *********************************************************************** */
{
  int    i,j,k,nv;
  int    beg_dir, end_dir;
  static unsigned char ***flag; 
  double invDt_par, Ehyp, Es, rho, p, *u;
  double Emag = 0.0, Ekin;
#if EOS == IDEAL
  double g1 = g_gamma - 1.0;
#endif
  static double ****rhs;
  
/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */

  if (rhs == NULL){
    rhs = ARRAY_4D(NX3_MAX, NX2_MAX, NX1_MAX, NVAR, double);
  }

/* --------------------------------------------------------
   1. Compute parabolic RHS when is d is not a NULL pointer
   -------------------------------------------------------- */

  if (d != NULL){
    invDt_par = ParabolicRHS(d, rhs, domBox, aflux, EXPLICIT,  1.0, grid);

    if (g_intStage == 1){
      #if (defined CTU) || (DIMENSIONAL_SPLITTING == YES)  
      Dts->invDt_par = MAX(Dts->invDt_par, invDt_par);
      #else
      invDt_par /= (double) DIMENSIONS;
      Dts->invDt_par = MAX(Dts->invDt_par, invDt_par);
      #endif
    }
    flag = d->flag;  /* Take the address of d->flag for later re-use */
  }

/* --------------------------------------------------------
   2. Loop over computational Box and update conservative 
      variables.
   -------------------------------------------------------- */

  BOX_LOOP(domBox, k,j,i){

  /* -- Entropy switch: recompute energy from entropy --  */
    #if ENTROPY_SWITCH && (EOS == IDEAL)
    if (flag[k][j][i] & FLAG_ENTROPY){
      u    = UU[k][j][i];
      rho  = u[RHO];      
      p    = u[ENTR]*pow(rho,g1);

      Ekin = 0.5*(EXPAND(u[MX1]*u[MX1], + u[MX2]*u[MX2], + u[MX3]*u[MX3]))/u[RHO];
      #if PHYSICS == MHD
      Emag = 0.5*(EXPAND(u[BX1]*u[BX1], + u[BX2]*u[BX2], + u[BX3]*u[BX3]));
      #endif
      u[ENG] = p/g1 + Ekin + Emag;
    }
    #endif

    #if VISCOSITY == EXPLICIT
    EXPAND(UU[k][j][i][MX1] += dt*rhs[k][j][i][MX1];  ,
           UU[k][j][i][MX2] += dt*rhs[k][j][i][MX2];  ,
           UU[k][j][i][MX3] += dt*rhs[k][j][i][MX3];)
    #endif
    #if (AMBIPOLAR_DIFFUSION == EXPLICIT) || (RESISTIVITY == EXPLICIT)
    EXPAND(UU[k][j][i][BX1] += dt*rhs[k][j][i][BX1];  ,
           UU[k][j][i][BX2] += dt*rhs[k][j][i][BX2];  ,
           UU[k][j][i][BX3] += dt*rhs[k][j][i][BX3];)
    #endif
    
    #if HAVE_ENERGY
    #if (AMBIPOLAR_DIFFUSION == EXPLICIT) ||\
        (RESISTIVITY        == EXPLICIT)  || \
        (THERMAL_CONDUCTION == EXPLICIT)  || \
        (VISCOSITY          == EXPLICIT) 

    UU[k][j][i][ENG] += dt*rhs[k][j][i][ENG];

    #endif
    #endif

  /* -- Entropy switch: recompute entropy from energy --  */
    #if ENTROPY_SWITCH && (EOS == IDEAL)
    if (flag[k][j][i] & FLAG_ENTROPY){
      Ekin = 0.5*(EXPAND(u[MX1]*u[MX1], + u[MX2]*u[MX2], + u[MX3]*u[MX3]))/u[RHO];
      #if PHYSICS == MHD
      Emag = 0.5*(EXPAND(u[BX1]*u[BX1], + u[BX2]*u[BX2], + u[BX3]*u[BX3]));
      #endif
      p = g1*(u[ENG] - Ekin - Emag); 
      u[ENTR] = p/pow(rho,g1);
    }
    #endif
  } /* End BOX_LOOP() */

}

/* ********************************************************************* */
double ParabolicRHS (const Data *d, Data_Arr dU, RBox *domBox, double **aflux,
                     int timeStepping, double dt, Grid *grid)
/*!  
 * \param [in]  d             Pointer to the PLUTO data structure.
 * \param [out] dU            3D array containing the conservative right hand sides 
 * \param [in]  domBox        A pointer to an RBox structure defining
 *                            the zones of the domain to be updated
 * \param [in]  aflux         A 2D pointer to store fluxes (needed by Chombo)
 * \param [in]  timeStepping  an integer specifying the time stepping
 *                            method (= EXPLICIT / STS / RKL)
 * \param [in]  dt            the time step
 * \param [in]  grid          pointer to Grid structure
 *
 * \return On output it returns the maximum diffusion coefficients 
 *         among all dissipative term over the local processor grid.
 *********************************************************************** */
{
  int     i, j, k, nv;
  int     nbeg, nend;
  int     includeDir[3], include[8];
  double  scrh;
  double  max_invDt_par = 0.0, invDt_par;
  static  double ***C_dtp[MAX_OP], *dcoeff, **dcoeff_res;
  
/* --------------------------------------------------------
   0. Allocate storage memory for sweep structure,
      area-weighted flux and diffusion coefficients.
     
      We use C_dt[AMB_DIFF_OP] for ambipolar diffusion, 
             C_dt[RES_OP+IDIR/JDIR/KDIR] for resistivity (eta_x/y/z),
             C_dt[TC_OP] for thermal conduction, etc...
/* -------------------------------------------------------- */

  if (dcoeff == NULL) {
    dcoeff  = ARRAY_1D(NMAX_POINT, double);
    dcoeff_res  = ARRAY_2D(3, NMAX_POINT, double);

    if (AMBIPOLAR_DIFFUSION) {
      C_dtp[AMB_DIFF_OP] = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    }
    if (HALL_MHD){
      C_dtp[HALL_OP] = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    }
    if (RESISTIVITY) {
      C_dtp[RES_OP+0] = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
      C_dtp[RES_OP+1] = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
      C_dtp[RES_OP+2] = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    }  
    if (THERMAL_CONDUCTION){
      C_dtp[TC_OP] = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    }
    if (VISCOSITY){
      C_dtp[VISC_OP] = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
    }
  }
 
  for (nv = 0; nv < MAX_OP; nv++) {
    if (C_dtp[nv] != NULL) TOT_LOOP(k,j,i) C_dtp[nv][k][j][i] = 0.0;
  }

#if DIMENSIONAL_SPLITTING == YES
  TOT_LOOP(k,j,i)  NVAR_LOOP(nv) dU[k][j][i][nv] = 0.0;
#endif

/* --------------------------------------------------------
   1. Select which operator(s) should be included during
      this call.
      When a diffusion operator is enabled, it may be
      called explicitly or using STS methods.
   -------------------------------------------------------- */

  include[AMB_DIFF_OP] = (AMBIPOLAR_DIFFUSION == timeStepping);
  include[RES_OP]      = (RESISTIVITY         == timeStepping);
  include[TC_OP]       = (THERMAL_CONDUCTION  == timeStepping);
  include[VISC_OP]     = (VISCOSITY           == timeStepping);

  if (timeStepping == EXPLICIT && DIMENSIONAL_SPLITTING == YES){
    includeDir[IDIR] = (g_dir == IDIR);
    includeDir[JDIR] = (g_dir == JDIR);
    includeDir[KDIR] = (g_dir == KDIR);
  }else{
    includeDir[IDIR] = 1;
    includeDir[JDIR] = (DIMENSIONS >= 2);
    includeDir[KDIR] = (DIMENSIONS == 3);
  }

  i = j = k = 0;

/* --------------------------------------------------------
   3. Compute current at cell edges before sweeping.
   -------------------------------------------------------- */
  
#if PHYSICS == MHD
  if (include[RES_OP] || include[AMB_DIFF_OP]) GetCurrent (d, grid);
#endif

#if THERMAL_CONDUCTION
  if (include[TC_OP]){   
    TOT_LOOP(k,j,i) d->Tc[k][j][i] = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i];
  }
#endif

/* --------------------------------------------------------
   4.  X1-Sweep (g_dir == IDIR)
   -------------------------------------------------------- */

  if (includeDir[IDIR]){

    g_dir = IDIR;

    KBOX_LOOP (domBox, k){ 
    JBOX_LOOP (domBox, j){ 

      g_j = j; g_k = k;
      nbeg = domBox->ibeg;
      nend = domBox->iend;

  /* -- Compute parabolic fluxes -- */

      #ifdef SHEARINGBOX
      /* With the shearingbox, diffusion terms are computed normally using
         grid and boudnary values as they are and no symmetrization is
         implemented at present. */
      #endif

  /* -- Start main X1-sweep -- */
  
      ITOT_LOOP(i) NVAR_LOOP(nv) dU[k][j][i][nv] = 0.0;

  /* -- Compute total parabolic flux -- */

      #if RESISTIVITY
      if (include[RES_OP]){
        ResistiveRHS (d, dU, dcoeff_res, aflux, dt, nbeg, nend, grid);
        if (g_intStage == 1){
          double *inv_dl = GetInverse_dl(grid);
          IBOX_LOOP (domBox, i){  
            double inv_dl2 = inv_dl[i]*inv_dl[i];

            EXPAND(C_dtp[RES_OP+0][k][j][i] += 0.5*( dcoeff_res[0][i-1]
                                                   + dcoeff_res[0][i])*inv_dl2;  ,
                   C_dtp[RES_OP+1][k][j][i] += 0.5*( dcoeff_res[1][i-1]
                                                   + dcoeff_res[1][i])*inv_dl2;  ,
                   C_dtp[RES_OP+2][k][j][i] += 0.5*(  dcoeff_res[2][i-1]
                                                    + dcoeff_res[2][i])*inv_dl2;)

            EXPAND(invDt_par     = dcoeff_res[0][i]*inv_dl2;
                   max_invDt_par = MAX(max_invDt_par, invDt_par);         , 

                   invDt_par     = dcoeff_res[1][i]*inv_dl2;
                   max_invDt_par = MAX(max_invDt_par, invDt_par);         , 

                   invDt_par     = dcoeff_res[2][i]*inv_dl2;
                   max_invDt_par = MAX(max_invDt_par, invDt_par);)

          }  
        }  
      }
      #endif

      #if THERMAL_CONDUCTION 
      if (include[TC_OP]){
        TC_RHS (d, dU, dcoeff, aflux, dt, nbeg, nend, grid);
        if (g_intStage == 1){
          double *inv_dl = GetInverse_dl(grid);
          IBOX_LOOP (domBox, i){  
            double inv_dl2 = inv_dl[i]*inv_dl[i];
            C_dtp[TC_OP][k][j][i] += 0.5*(dcoeff[i-1] + dcoeff[i])*inv_dl2;
            invDt_par     = dcoeff[i]*inv_dl2;
            max_invDt_par = MAX(max_invDt_par, invDt_par);
          }  
        }  
      }
      #endif

      #if VISCOSITY
      if (include[VISC_OP]) { 
        ViscousRHS (d, dU, dcoeff, aflux, dt, nbeg, nend, grid);
        if (g_intStage == 1){
          double *inv_dl = GetInverse_dl(grid);
          IBOX_LOOP (domBox, i){  
            double inv_dl2 = inv_dl[i]*inv_dl[i];
            C_dtp[VISC_OP][k][j][i] += 0.5*(dcoeff[i-1] + dcoeff[i])*inv_dl2;
            invDt_par     = dcoeff[i]*inv_dl2;
            max_invDt_par = MAX(max_invDt_par, invDt_par);
          }  
        }  
      }
      #endif /* VISCOSITY */
    }}
  } /* end if (includeDir(IDIR)) */

/* --------------------------------------------------------
   5.  X2-Sweep (g_dir == JDIR)
   -------------------------------------------------------- */

  if (includeDir[JDIR]){

    g_dir = JDIR;

    KBOX_LOOP (domBox,k){
    IBOX_LOOP (domBox,i){

      g_i = i; g_k = k;
      nbeg = domBox->jbeg;
      nend = domBox->jend;

  /* -- Compute total parabolic flux -- */

      #if RESISTIVITY
      if (include[RES_OP]){
        ResistiveRHS (d, dU, dcoeff_res, aflux, dt, nbeg, nend, grid);
        if (g_intStage == 1){
          double *inv_dl = GetInverse_dl(grid);
          JBOX_LOOP (domBox, j){  
            double inv_dl2 = inv_dl[j]*inv_dl[j];

            EXPAND(C_dtp[RES_OP+0][k][j][i] += 0.5*( dcoeff_res[0][j-1]
                                                  + dcoeff_res[0][j])*inv_dl2;  ,
                   C_dtp[RES_OP+1][k][j][i] += 0.5*( dcoeff_res[1][j-1]
                                                  + dcoeff_res[1][j])*inv_dl2;  ,
                   C_dtp[RES_OP+2][k][j][i] += 0.5*(  dcoeff_res[2][j-1]
                                                   + dcoeff_res[2][j])*inv_dl2;)

            EXPAND(invDt_par     = dcoeff_res[0][j]*inv_dl2;
                   max_invDt_par = MAX(max_invDt_par, invDt_par);         , 

                   invDt_par     = dcoeff_res[1][j]*inv_dl2;
                   max_invDt_par = MAX(max_invDt_par, invDt_par);         , 

                   invDt_par     = dcoeff_res[2][j]*inv_dl2;
                   max_invDt_par = MAX(max_invDt_par, invDt_par);)

          }  
        }  
      }
      #endif
  
      #if THERMAL_CONDUCTION 
      if (include[TC_OP]){
        TC_RHS (d, dU, dcoeff, aflux, dt, nbeg, nend, grid);
        if (g_intStage == 1){  
          double *inv_dl = GetInverse_dl(grid);
          JBOX_LOOP (domBox, j){  
            double inv_dl2 = inv_dl[j]*inv_dl[j];
            C_dtp[TC_OP][k][j][i] += 0.5*(dcoeff[j-1] + dcoeff[j])*inv_dl2;
            invDt_par     = dcoeff[j]*inv_dl2;
            max_invDt_par = MAX(max_invDt_par, invDt_par);
          }  
        }  
      }
      #endif

      #if VISCOSITY
      if (include[VISC_OP]) { 
        ViscousRHS (d, dU, dcoeff, aflux, dt, nbeg, nend, grid);
        if (g_intStage == 1){  
          double *inv_dl = GetInverse_dl(grid);
          JBOX_LOOP (domBox, j){  
            double inv_dl2 = inv_dl[j]*inv_dl[j];
            C_dtp[VISC_OP][k][j][i] += 0.5*(dcoeff[j-1] + dcoeff[j])*inv_dl2;
            invDt_par     = dcoeff[j]*inv_dl2;
            max_invDt_par = MAX(max_invDt_par, invDt_par);
          }  
        }  
      }  
      #endif /* VISCOSITY */
    }}
  }  /* end JDIR   */

/* --------------------------------------------------------
   6.  X3-Sweep (g_dir == KDIR)
   -------------------------------------------------------- */

  if (includeDir[KDIR]){

    g_dir = KDIR;

    JBOX_LOOP (domBox, j){
    IBOX_LOOP (domBox, i){

      g_i = i; g_j = j;
      nbeg = domBox->kbeg;
      nend = domBox->kend;

  /* -- Compute total parabolic flux -- */

      #if RESISTIVITY
      if (include[RES_OP]){
        ResistiveRHS (d, dU, dcoeff_res, aflux, dt, nbeg, nend, grid);
        if (g_intStage == 1){
          double *inv_dl = GetInverse_dl(grid);
          KBOX_LOOP (domBox, k){  
            double inv_dl2 = inv_dl[k]*inv_dl[k];

            EXPAND(C_dtp[RES_OP+0][k][j][i] += 0.5*( dcoeff_res[0][k-1]
                                                  + dcoeff_res[0][k])*inv_dl2;  ,
                   C_dtp[RES_OP+1][k][j][i] += 0.5*( dcoeff_res[1][k-1]
                                                  + dcoeff_res[1][k])*inv_dl2;  ,
                   C_dtp[RES_OP+2][k][j][i] += 0.5*(  dcoeff_res[2][k-1]
                                                   + dcoeff_res[2][k])*inv_dl2;)

            EXPAND(invDt_par     = dcoeff_res[0][k]*inv_dl2;
                   max_invDt_par = MAX(max_invDt_par, invDt_par);     , 

                   invDt_par     = dcoeff_res[1][k]*inv_dl2;
                   max_invDt_par = MAX(max_invDt_par, invDt_par);     , 

                   invDt_par     = dcoeff_res[2][k]*inv_dl2;
                   max_invDt_par = MAX(max_invDt_par, invDt_par);)
          }  
        }  
      }
      #endif
  
      #if THERMAL_CONDUCTION 
      if (include[TC_OP]){
        TC_RHS (d, dU, dcoeff, aflux, dt, nbeg, nend, grid);
        if (g_intStage == 1){  
          double *inv_dl = GetInverse_dl(grid);
          KBOX_LOOP (domBox, k){  
            double inv_dl2 = inv_dl[k]*inv_dl[k];
            C_dtp[TC_OP][k][j][i] += 0.5*(dcoeff[k-1] + dcoeff[k])*inv_dl2;
            invDt_par     = dcoeff[k]*inv_dl2;
            max_invDt_par = MAX(max_invDt_par, invDt_par);
          }  
        }  
      }
      #endif

      #if VISCOSITY
      if (include[VISC_OP]) { 
        ViscousRHS (d, dU, dcoeff, aflux, dt, nbeg, nend, grid);
        if (g_intStage == 1){  
          double *inv_dl = GetInverse_dl(grid);
          KBOX_LOOP (domBox, k){  
            double inv_dl2 = inv_dl[k]*inv_dl[k];
            C_dtp[VISC_OP][k][j][i] += 0.5*(dcoeff[k-1] + dcoeff[k])*inv_dl2;
            invDt_par     = dcoeff[k]*inv_dl2;
            max_invDt_par = MAX(max_invDt_par, invDt_par);
          }  
        }  
      }
      #endif /* VISCOSITY */
    }}
  }  /* end KDIR */
  
/* --------------------------------------------------------
   7. Take the maximum of inverse dt over domain and zero 
      right hand side in the internal boundary zones.
      Note: for CTU or dimensionally split methods, we
            take the maximum at the cell interface
            (the same is done for hyperbolic terms).
   -------------------------------------------------------- */

  if (timeStepping == EXPLICIT){
    #if (defined CTU) || DIMENSIONAL_SPLITTING == YES
    return max_invDt_par;
    #endif
  }

  scrh = 0.0;
  BOX_LOOP(domBox, k,j,i){
    #if AMBIPOLAR_DIFFUSION
    if (include[AMB_DIFF_OP]){
      scrh = MAX(scrh, C_dtp[AMB_DIFF_OP][k][j][i]);
    }
    #endif

    #if RESISTIVITY
    if (include[RES_OP]){
      EXPAND(scrh = MAX(scrh, C_dtp[RES_OP+0][k][j][i]);   ,  
             scrh = MAX(scrh, C_dtp[RES_OP+1][k][j][i]);   ,
             scrh = MAX(scrh, C_dtp[RES_OP+2][k][j][i]);)
    }
    #endif

    #if THERMAL_CONDUCTION
    if (include[TC_OP]) scrh = MAX(scrh, C_dtp[TC_OP][k][j][i]);
    #endif

    #if VISCOSITY
    if (include[VISC_OP]) scrh = MAX(scrh, C_dtp[VISC_OP][k][j][i]);
    #endif

    #if INTERNAL_BOUNDARY == YES
    if (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY) {
      NVAR_LOOP(nv) dU[k][j][i][nv] = 0.0;
    }
    #endif
  }
  return scrh;
}

#endif /* PARABOLIC_FLUX != NO */
