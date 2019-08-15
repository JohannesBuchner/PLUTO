/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief PLUTO main function.
 
  The file main.c contains the PLUTO main function and several other 
  top-level routines.
  main() provides basic code initialization, handles the the principal 
  integration loop and calls the output driver write_data.c.
  Other useful functions contained in this file are Integrate() which does
  the actual integration, NextTimeStep() responsible for computing the
  next time step based on the information available at the last time
  level.
  
  The standard integration loop consists of the following steps:
 
   - Check for last step & adjust dt if necessary
   - Dump log information, n, t(n), dt(n), MAX_MACH(n-1), etc..
   - Check output/analysis:  t(n) < tout < t(n)+dt(n)
   - write to disk/call analysis using {U(n), t(n), dt(n)}
   - Advance solution using dt(n): U(n) --> U(n+1)
   - Increment t(n+1) = t(n) + dt(n)
   - [MPI] Show dominant time step (n)
   - Get next time step dt(n+1)
   - [MPI] reduction operations (n)
   - Increment n --> n+1
 
  \author A. Mignone (mignone@ph.unito.it)
  \date   April 9, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "globals.h"

#ifndef SHOW_TIME_STEPS
  #define SHOW_TIME_STEPS  NO  /* Show time steps due to different processes */
#endif
#ifndef SHOW_TIMING
  #define SHOW_TIMING      NO  /* Compute CPU timing between steps */
#endif

static double NextTimeStep (timeStep *, Runtime *, Grid *);
static char *TotalExecutionTime (double);
static int Integrate (Data *, Riemann_Solver *, timeStep *, Grid *);
static void CheckForOutput (Data *, Runtime *, time_t, Grid *);
static void CheckForAnalysis (Data *, Runtime *, Grid *);

/* ********************************************************************* */
int main (int argc, char *argv[])
/*!
 * Start PLUTO, initialize functions, define data structures and 
 * handle the main integration loop.
 *
 * \param [in] argc Argument counts.
 * \param [in] argv Array of pointers to the strings.
 * \return This function return 0 on normal exit.
 *
 *********************************************************************** */
{
  int    nv, idim, err;
  char   first_step=1, last_step = 0;
  double scrh;
  Data   data;
  clock_t clock_beg, clock_end; /* measures CPU time */
  time_t  tbeg, tend;           /* measures the real time */
  Riemann_Solver *Solver;
  Grid      grd[3];
  timeStep Dts;
  Cmd_Line cmd_line;
  Runtime  ini;
  Output *output;

/* --------------------------------------------------------
   0. Initialize environment
   -------------------------------------------------------- */

#ifdef PARALLEL
  AL_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &prank);
#endif

  time (&tbeg);
  Initialize (argc, argv, &data, &ini, grd, &cmd_line);

  data.Dts = &Dts;

/* -- 0a. Initialize members of timeStep structure -- */

  Dts.cmax      = ARRAY_1D(NMAX_POINT, double);
  Dts.invDt_hyp = 0.0;
  Dts.invDt_par = 0.5e-38; 
  Dts.dt_cool   = 1.e38;
  Dts.cfl       = ini.cfl;
  Dts.cfl_par   = ini.cfl_par;
  Dts.rmax_par  = ini.rmax_par;
  Dts.Nsts      = Dts.Nrkc = Dts.Nrkl = 0;
#if (defined PARTICLES) && (PARTICLES_TYPE == COSMIC_RAYS)
  Dts.Nsub_particles  = MAX(1, -PARTICLES_CR_NSUB);
#else
  Dts.Nsub_particles = 1;
#endif

  Dts.invDt_particles = 1.0/ini.first_dt;
  
  Solver = SetSolver (ini.solv_type);
  
  g_stepNumber = 0;
  
/* --------------------------------------------------------
   1. Check if restart is necessary. 
      If not, write initial condition to disk.
   ------------------------------------------------------- */
  
  if (cmd_line.restart == YES) {
    RestartFromFile (&ini, cmd_line.nrestart, DBL_OUTPUT, grd);
    #ifdef PARTICLES
    Particles_Restart(&data, cmd_line.nrestart, grd);
    #endif
  }else if (cmd_line.h5restart == YES){
    RestartFromFile (&ini, cmd_line.nrestart, DBL_H5_OUTPUT, grd);
  }else if (cmd_line.write){
    CheckForOutput (&data, &ini, tbeg, grd);
    CheckForAnalysis (&data, &ini, grd);
  }

  if (cmd_line.maxsteps == 0) last_step = 1;
  print ("> Starting computation... \n\n");  

/* =====================================================================
   2.  M A I N      L O O P      S T A R T S      H E R E
   ===================================================================== */

  while (!last_step){

    #if SHOW_TIMING
    clock_beg = clock();
    #endif

  /* ------------------------------------------------------
     2a. Check if this is the last integration step:
         - final tstop has been reached: adjust time step 
         - or max number of steps has been reached
     ------------------------------------------------------ */

    if ((g_time + g_dt) >= ini.tstop*(1.0 - 1.e-8)) {
      g_dt      = (ini.tstop - g_time);
      last_step = 1;
    }
    if (g_stepNumber == cmd_line.maxsteps && cmd_line.maxsteps >= 0) {
      last_step = 1;
    }

  /* ------------------------------------------------------
     2b. Update log file
     ------------------------------------------------------ */

    if (g_stepNumber%ini.log_freq == 0) OutputLogPre(&data, &Dts, &ini, grd);

  /* ----------------------------------------------------
     2c. Check if it's time to write or perform analysis
     ---------------------------------------------------- */
     
    if (!first_step && !last_step && cmd_line.write) {
      CheckForOutput  (&data, &ini, tbeg, grd);
      CheckForAnalysis(&data, &ini, grd);
    }

  /* ----------------------------------------------------
     2d. Advance solution array by a single time step
         g_dt = dt(n). After this step U^n -> U^{n+1}
     ---------------------------------------------------- */

    if (cmd_line.jet != -1) SetJetDomain (&data, cmd_line.jet, ini.log_freq, grd); 
    err = Integrate (&data, Solver, &Dts, grd);
    if (cmd_line.jet != -1) UnsetJetDomain (&data, cmd_line.jet, grd); 

  /* ----------------------------------------------------
     2e. Integration didn't go through. Step must
         be redone from previously saved solution.
     ---------------------------------------------------- */
/*
    if (err != 0){
      print ("! Step failed. Re-trying\n");
      zones with problems must be tagged with MINMOD_FLAG and HLL_FLAG
      time step should be halved
      GET_SOL(&data);
    }
*/

  /* ------------------------------------------------------
     2h. Global MPI reduction operations
     ------------------------------------------------------ */
  
    #ifdef PARALLEL
    MPI_Allreduce (&g_maxMach, &scrh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    g_maxMach = scrh;

    MPI_Allreduce (&g_maxRiemannIter, &nv, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    g_maxRiemannIter = nv;
    
    #if (PHYSICS == RMHD) && (RESISTIVITY != NO)
    MPI_Allreduce (&g_maxIMEXIter, &nv, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    g_maxIMEXIter = nv;
    #endif
    #endif

    if (g_stepNumber%ini.log_freq == 0) {
      OutputLogPost(&data, &Dts, &ini, grd);
    }

  /* ----------------------------------------------------
     2f. Increment time, t(n+1) = t(n) + dt(n)
     ---------------------------------------------------- */ 

    g_time += g_dt;

    #if SHOW_TIMING
    clock_end = clock();
    if (g_stepNumber%ini.log_freq == 0) {
      scrh = (double)(clock_end - clock_beg)/CLOCKS_PER_SEC;
      print("%s [clock (total)         = %f (s)]\n",IndentString(), scrh);
      print("%s [clock (AdvanceStep()) = %f (s)]\n",IndentString(),Dts.clock_hyp);
      #ifdef PARTICLES 
      print("%s [clock (particles)     = %f (s)]\n",IndentString(),Dts.clock_particles);
      #endif
    }
    #endif

  /* ------------------------------------------------------
     2g. Get next time step dt(n+1).
         Do it every two steps if cooling or dimensional
         splitting are used.
     ------------------------------------------------------ */

    #if (COOLING == NO) && ((DIMENSIONS == 1) || (DIMENSIONAL_SPLITTING == NO))
    g_dt = NextTimeStep(&Dts, &ini, grd);
    #else
    if (g_stepNumber%2 == 1) g_dt = NextTimeStep(&Dts, &ini, grd);
    #endif

    g_stepNumber++;
    first_step = 0;
  }

/* =====================================================================
          M A I N       L O O P      E N D S       H E R E 
   ===================================================================== */

  /*  Prevent double output when maxsteps == 0 */
  if ((cmd_line.write) && !(cmd_line.maxsteps == 0)){
    CheckForOutput (&data, &ini, tbeg, grd);
    CheckForAnalysis (&data, &ini, grd);
  }

  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
   print  ("\n> Total allocated memory  %6.2f Mb (proc #%d)\n",
             (float)g_usedMemory/1.e6,prank);
   MPI_Barrier (MPI_COMM_WORLD);
  #else
   print  ("\n> Total allocated memory  %6.2f Mb\n",(float)g_usedMemory/1.e6);
  #endif

  time(&tend);
  g_dt = difftime(tend, tbeg);
  print("> Elapsed time             %s\n", TotalExecutionTime(g_dt));

  /*  Check if stepNumber = 0. 
   *  Prevent NaN print at maxsteps = 0 */
  if (g_stepNumber > 0)
      print("> Average time/step       %10.2e  (sec)  \n", 
           difftime(tend,tbeg)/(double)g_stepNumber);
  else print("> Average time/step       %10.2e  (sec)  \n",difftime(tend,tbeg));

  print("> Local time                %s",asctime(localtime(&tend)));
  print("> Done\n");

  FreeArray4D ((void *) data.Vc);
  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
   AL_Finalize ();
  #endif

  return (0);
}

/* ********************************************************************* */
int Integrate (Data *d, Riemann_Solver *Solver, timeStep *Dts, Grid *grid)
/*!
 * Advance equations by a single time-step.
 *
 * \param  d      pointer to PLUTO Data structure;
 * \param  Solver pointer to a Riemann solver function;
 * \param  Dts    pointer to time Step structure;
 * \param  grid   pointer to grid structure.
 * 
 * \return An integer giving success / failure (development).
 * 
 *********************************************************************** */
{
  int    idim, err = 0;
  int    i,j,k;

  g_maxMach        = 0.0;
  g_maxRiemannIter = 0;
  g_maxRootIter    = 0;
	
#ifdef PARTICLES 
  #if (PARTICLES_TEST == YES)
  Particles_CR_Update (d, Dts, g_dt, grid);
  return 0;
  #endif
#endif

/* --------------------------------------------------------
   1. Initialize max propagation speed in Dedner's approach
   -------------------------------------------------------- */

#ifdef GLM_MHD  /* -- initialize glm_ch -- */
  GLM_Init (d, Dts, grid);   
  GLM_Source (d->Vc, 0.5*g_dt, grid);
#endif

/* --------------------------------------------------------
   2. Perform Strang Splitting on directions (if necessary)
      and sources 
   -------------------------------------------------------- */

  TOT_LOOP(k,j,i) d->flag[k][j][i] = 0;

#ifdef FARGO
  FARGO_ComputeVelocity(d, grid);
#endif
  
  if ((g_stepNumber%2) == 0){
    g_hydroStep = TRUE;
    #if DIMENSIONAL_SPLITTING == YES
    for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){
      if (AdvanceStep (d, Solver, Dts, grid) != 0) return (1);
    }
    #else
    if (AdvanceStep (d, Solver, Dts, grid) != 0) return(1);
    #endif
    g_hydroStep = FALSE;
    SplitSource (d, g_dt, Dts, grid);
  }else{
    g_hydroStep = FALSE;
    SplitSource (d, g_dt, Dts, grid);
    g_hydroStep = TRUE;
    #if DIMENSIONAL_SPLITTING == YES
    for (g_dir = DIMENSIONS - 1; g_dir >= 0; g_dir--){
      if (AdvanceStep(d, Solver, Dts, grid) != 0) return (1);
    }
    #else
    if (AdvanceStep (d, Solver, Dts, grid) != 0) return(1);
    #endif
  }       

#ifdef GLM_MHD  /* -- GLM source for dt/2 -- */
  GLM_Source (d->Vc, 0.5*g_dt, grid);
#endif
  g_hydroStep = 0;

  return (0); /* -- ok, step achieved -- */
}

/* ********************************************************************* */
char *TotalExecutionTime (double dt)
/*!
 *
 *   convert a floating-point variable (dt, in seconds) to a string 
 *   displaying days:hours:minutes:seconds
 *   
 *********************************************************************** */
{
  static char c[128];
  int days, hours, mins, secs;

  days  = (int) (dt/86400.0);
  hours = (int) ((dt - 86400.0*days)/3600.0);
  mins  = (int) ((dt - 86400.0*days - 3600.0*hours)/60.);
  secs  = (int) (dt - 86400.0*days - 3600.0*hours - 60.0*mins);

  sprintf (c, " %dd:%dh:%dm:%ds", days,hours, mins, secs);
  return (c);
}

/* ********************************************************************* */
double NextTimeStep (timeStep *Dts, Runtime *ini, Grid *grid)
/*!
 * Compute and return the time step for the next time level
 * using the information from the previous integration
 * (Dts->invDt_hyp and Dts->invDt_par).
 *
 * \param [in] Dts    pointer to the timeStep structure
 * \param [in] ini    pointer to the Runtime structure
 * \param [in] grid   pointer to array of Grid structures
 *
 * \return The time step for next time level
 *********************************************************************** */
{
  int    idim;
  double dt_hyp, dt_par, dt_particles, dtnext;
  double scrh;
  double dxmin;
  double xloc, xglob;

/* --------------------------------------------------------
   1. Take the maximum of invDt_hyp, invDt_par, etc...
      across all processors
   -------------------------------------------------------- */

#ifdef PARALLEL
  xloc = Dts->invDt_hyp;
  MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  Dts->invDt_hyp = xglob;
  #if (PARABOLIC_FLUX != NO)
  xloc = Dts->invDt_par;
  MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  Dts->invDt_par = xglob;
  #endif
  #if COOLING != NO
  xloc = Dts->dt_cool;
  MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  Dts->dt_cool = xglob;
  #endif
  #ifdef PARTICLES
  xloc = Dts->invDt_particles;
  MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  Dts->invDt_particles = xglob;
  #endif
#endif

/* --------------------------------------------------------
   2. Show the time step ratios between the actual g_dt
      and the advection, diffusion and cooling time scales.
   -------------------------------------------------------- */

#if SHOW_TIME_STEPS == YES
  if (g_stepNumber%ini->log_freq == 0) {
    char *str = IndentString();
    print("%s [dt(adv)       = cfl x %10.4e]\n",str, 1.0/Dts->invDt_hyp);
    #if PARABOLIC_FLUX != NO
    print("%s [dt(par)       = cfl x %10.4e]\n",str, 1.0/(2.0*Dts->invDt_par));
    #endif
    #if COOLING != NO
    print("%s [dt(cool)      =       %10.4e]\n",str, Dts->dt_cool);
    #endif
    #ifdef PARTICLES 
    print("%s [dt(particles) =       %10.4e]\n",str, 1.0/Dts->invDt_particles);
    #endif
  }
#endif

/* --------------------------------------------------------
   3. Compute hyperbolic (htydro) time step
   -------------------------------------------------------- */

#if (PARABOLIC_FLUX & EXPLICIT)
  dt_hyp  = 1.0/(Dts->invDt_hyp + 2.0*Dts->invDt_par);
#else
  dt_hyp  = 1.0/Dts->invDt_hyp;
#endif
  dt_hyp *= ini->cfl;
  dtnext  = dt_hyp;

/* --------------------------------------------------------
   4. For GLM-MHD recompute maximum propagation speed from
      dtnext and update glm_ch.
      However, for the first time, glm_ch is recomputed
      using GLM_Init().
      Note: in the case of relativistic flows, glm_ch
      can become less than 1 (when cfast << 1) allowing
      larger time steps to be taken.
   -------------------------------------------------------- */

#ifdef GLM_MHD
  dxmin = grid->dl_min[IDIR];
  for (idim = 1; idim < DIMENSIONS; idim++){ /*  Min cell length   */
    dxmin = MIN(dxmin, grid->dl_min[idim]);
  }
  glm_ch = ini->cfl*dxmin/dtnext;
/*
scrh = glm_ch;
glm_ch = dxmin*Dts->invDt_hyp;
if (fabs(scrh-glm_ch) > 1.e-12){
  print ("! Different values of glm_ch  = %12.6e,  %12.6e\n",glm_ch,scrh);
  QUIT_PLUTO(1);
}
*/
#endif

/* --------------------------------------------------------
   5. With STS-like methods, the ratio between advection
      (full) and parabolic time steps should not exceed
      ini->rmax_par.
   -------------------------------------------------------- */
      
#if (PARABOLIC_FLUX & SUPER_TIME_STEPPING) || \
    (PARABOLIC_FLUX & RK_CHEBYSHEV)        || \
    (PARABOLIC_FLUX & RK_LEGENDRE)
   dt_par  = ini->cfl_par/(2.0*Dts->invDt_par);
   dtnext *= MIN(1.0, ini->rmax_par/(dt_hyp/dt_par));
#endif

/* --------------------------------------------------------
   6. Compute Cooling time step
   -------------------------------------------------------- */

#if COOLING != NO
  dtnext = MIN(dtnext, Dts->dt_cool);
#endif

/* --------------------------------------------------------
   7. Limit time step due to particles
      (only if PARTICLES_CR_NSUB > 0)
   -------------------------------------------------------- */

#if    (defined PARTICLES) && (PARTICLES_TYPE == COSMIC_RAYS) \
                           && (PARTICLES_CR_NSUB > 0)
  dt_particles = PARTICLES_CR_NSUB/Dts->invDt_particles;
  dtnext       = MIN(dtnext, dt_particles);
#endif

/* --------------------------------------------------------
   8. Allow time step to increase at most by a factor
      ini->cfl_max_var.
      Quit if dt gets too small, issue a warning if
      first_dt has been overestimated.
   -------------------------------------------------------- */

  dtnext = MIN(dtnext, ini->cfl_max_var*g_dt);

  if (dtnext < ini->first_dt*1.e-9){
    char *str = IndentString();
    print("! NextTimeStep(): dt is too small (%12.6e).\n", dtnext);
    print("! %s [dt(adv)       = cfl x %10.4e]\n",str, 1.0/Dts->invDt_hyp);
    print("! %s [dt(par)       = cfl x %10.4e]\n",str, 1.0/(2.0*Dts->invDt_par));
    print("! %s [dt(cool)      =       %10.4e]\n",str, Dts->dt_cool);
    print("! %s [dt(particles) =       %10.4e]\n",str, 1.0/Dts->invDt_particles);
    print("! Cannot continue.\n");
    QUIT_PLUTO(1);
  }

  if (g_stepNumber <= 1 && (ini->first_dt > dtnext/ini->cfl)){
    print ("! NextTimeStep(): initial dt exceeds stability limit\n");
  }

/* --------------------------------------------------------
   9. Reset time step coefficients
   -------------------------------------------------------- */

#if (defined PARTICLES) && (PARTICLES_TYPE == COSMIC_RAYS)
  #if PARTICLES_CR_NSUB > 0
  Dts->Nsub_particles = MIN(ceil(dtnext*Dts->invDt_particles), PARTICLES_CR_NSUB);
  #if PARTICLES_CR_PREDICTOR == 2
  if (Dts->Nsub_particles > 2){
    Dts->Nsub_particles += (Dts->Nsub_particles%2 == 1);
  }
  #endif
  #elif PARTICLES_CR_NSUB < 0                  
  Dts->Nsub_particles = abs(PARTICLES_CR_NSUB);  /* Fix Nsub no matter what */
  #else
  #error ! NextTimeStep(): PARTICLES_CR_NSUB cannot be = 0
  #endif
#endif

  DIM_LOOP(idim) Dts->cmax[idim] = 0.0;
  Dts->invDt_hyp       = 0.0;
  Dts->invDt_par       = 0.5e-38;
  Dts->invDt_particles = 1.e-38;
  Dts->dt_cool         = 1.e38;

  return(dtnext);
}

/* ********************************************************************* */
void CheckForOutput (Data *d, Runtime *ini, time_t t0, Grid *grid)
/*!
 *  Check if file output has to be performed.
 *  
 *********************************************************************** */
{
  static int first_call = 1;
  int    n, check_dt, check_dn, check_dclock;
  int    first_step     = (g_stepNumber == 0 ? 1:0);
  int    last_step      = (fabs(g_time-ini->tstop) < 1.e-12 ? 1:0);
  int    restart_update = 0;
  double dtime[MAX_OUTPUT_TYPES];
  static time_t tbeg[MAX_OUTPUT_TYPES], tend;
  Output *output;

/* ----------------------------------------------
   0. Initialize
   ---------------------------------------------- */

  if (first_call) for (n = 0; n < MAX_OUTPUT_TYPES; n++) tbeg[n] = t0;

  time (&tend);

  for (n = 0; n < MAX_OUTPUT_TYPES; n++) dtime[n] = difftime(tend, tbeg[n]);
#ifdef PARALLEL
  MPI_Bcast(dtime, MAX_OUTPUT_TYPES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  
/* -------------------------------------------------------
   1. Start main loop on outputs
   ------------------------------------------------------- */

  for (n = 0; n < MAX_OUTPUT_TYPES; n++){
    output = ini->output + n;    
    check_dt = check_dn = check_dclock = 0; 

  /* ----------------------------------------------------
     1a. Enable writing if  t(n) < t(write) < t(n+1)
     ---------------------------------------------------- */

    if (output->dt > 0.0){
      check_dt = (int) ((g_time+g_dt)/output->dt) - (int)(g_time/output->dt);
      check_dt = check_dt || first_step || last_step;
    }   

  /* ----------------------------------------------------
     1b. Enable writing if current step is an integer
         multiple of output interval
     ---------------------------------------------------- */

    if (output->dn > 0){
      check_dn = (g_stepNumber%output->dn) == 0;
      check_dn = check_dn || first_step || last_step;
    }

  /* ----------------------------------------------------
     1c. Enable writing if clock time has been reached
     ---------------------------------------------------- */

    if (output->dclock > 0.0){
      if (dtime[n] >= output->dclock) {
        check_dclock = 1;
        tbeg[n] = tend;
      }else{ 
        check_dclock = 0;
      }
      check_dclock = check_dclock || first_step || last_step;
    }

  /* -- if any of the previous is true dump data to disk -- */

    if (check_dt || check_dn || check_dclock) { 

#ifdef PARTICLES
      if (  output->type == PARTICLES_DBL_OUTPUT
         || output->type == PARTICLES_FLT_OUTPUT
         || output->type == PARTICLES_VTK_OUTPUT
         || output->type == PARTICLES_TAB_OUTPUT
         || output->type == PARTICLES_HDF5_OUTPUT){
        Particles_WriteData(d, output, grid);
      } else
#endif
      WriteData(d, output, grid);  

    /* ----------------------------------------------------------
        save the file number of the dbl and dbl.h5 output format
        for writing restart.out once we exit the loop.
       ---------------------------------------------------------- */

      if ((output->type == DBL_OUTPUT) ||
          (output->type == DBL_H5_OUTPUT)) restart_update = 1;
    }
  }

/* -------------------------------------------------------
   2. Dump restart information if required 

      Note that if both dbl and dbl.h5 formats are used,
      bookeeping is done using dbl format.
   ------------------------------------------------------- */

  if (restart_update) RestartDump (ini);

  first_call = 0;
}

/* ******************************************************************** */
void CheckForAnalysis (Data *d, Runtime *ini, Grid *grid)
/*
 *
 * PURPOSE 
 *
 *   Check if Analysis needs to be called
 *
 ********************************************************************** */
{
  int check_dt, check_dn;
  double t, tnext;

  t     = g_time;
  tnext = t + g_dt;
  check_dt = (int) (tnext/ini->anl_dt) - (int)(t/ini->anl_dt);
  check_dt = check_dt || g_stepNumber == 0 || fabs(t - ini->tstop) < 1.e-9; 
  check_dt = check_dt && (ini->anl_dt > 0.0);

  check_dn = (g_stepNumber%ini->anl_dn) == 0;
  check_dn = check_dn && (ini->anl_dn > 0);

  if (check_dt || check_dn) Analysis (d, grid);
}

