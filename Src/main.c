/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief PLUTO main function.
 
  The file main.c contains the PLUTO main function and several other 
  top-level routines.
  main() provides basic code initialization, handles the the principal 
  integration loop and calls the output driver write_data.c.
  Other useful functions contained in this file are Integrate() which does
  the actual integration, GetNextTimeStep() responsible for computing the
  next time step based on the information available at the last time
  level.
  
  We use two slightly different integration loops depending on whether
  asnchrounous I/O has to be performed (macro USE_ASYNC_IO).
  If the macro USE_ASYNC_IO is not defined, the standard 
  integration loop consists of the following steps:
 
   - Check for last step & adjust dt if necessary
   - Dump log information, n, t(n), dt(n), MAX_MACH(n-1), etc..
   - Check output/analysis:  t(n) < tout < t(n)+dt(n)
   - write to disk/call analysis using {U(n), t(n), dt(n)}
   - Advance solution using dt(n): U(n) --> U(n+1)
   - Increment t(n+1) = t(n) + dt(n)
   - [MPI] Show dominant time step (n)
   - [MPI] Get next time step dt(n+1)
   - [MPI] reduction operations (n)
   - Increment n --> n+1
 
  Otherwise, using Asynchrounous I/O:
 
   - Check for last step & adjust dt
   - check for output/analysis:   t(n) < tout < t(n+1)
     - Write data/call analysis using {U(n), t(n), dt(n)}
   - [MPI] Show dominant time step (n-1)
   - [MPI] reduction operations (n-1)
   - [AIO]: finish writing
   - Dump log information, n, t(n), dt(n), MAX_MACH(n-1), etc..
   - Advance solution using dt(n), U(n) --> U(n+1)
   - Increment t(n+1) = t(n) + dt(n)
   - [MPI] Get next time step dt(n)
   - Increment n --> n+1

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "globals.h"

#ifndef SHOW_TIME_STEPS
  #define SHOW_TIME_STEPS  NO   /* -- show time steps due to advection,
                                       diffusion and cooling */
#endif

static double NextTimeStep (Time_Step *, Runtime *, Grid *);
static char *TotalExecutionTime (double);
static int Integrate (Data *, Riemann_Solver *, Time_Step *, Grid *);
static void CheckForOutput (Data *, Runtime *, Grid *);
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
  time_t  tbeg, tend;
  Riemann_Solver *Solver;
  Grid      grd[3];
  Time_Step Dts;
  Cmd_Line cmd_line;
  Runtime  ini;
  Output *output;
  double *dbl_pnt;
  int    *int_pnt;

  #ifdef PARALLEL
   AL_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &prank);
  #endif

  Initialize (argc, argv, &data, &ini, grd, &cmd_line);

  print1 ("> Basic data type:\n");
  print1 ("  sizeof (char)     = %d\n", sizeof(char));
  print1 ("  sizeof (uchar)    = %d\n", sizeof(unsigned char));
  print1 ("  sizeof (short)    = %d\n", sizeof(short));
  print1 ("  sizeof (ushort)   = %d\n", sizeof(unsigned short));
  print1 ("  sizeof (int)      = %d\n", sizeof(int));
  print1 ("  sizeof (*int)     = %d\n", sizeof(int_pnt));
  print1 ("  sizeof (float)    = %d\n", sizeof(float));
  print1 ("  sizeof (double)   = %d\n", sizeof(double));
  print1 ("  sizeof (*double)  = %d\n", sizeof(dbl_pnt));
  
/*
  print1 ("\n> Structure data type:\n");
  print1 ("  sizeof (CMD_LINE)   = %d\n", sizeof(Cmd_Line));
  print1 ("  sizeof (DATA)       = %d\n", sizeof(Data));
  print1 ("  sizeof (STATE_1D)   = %d\n", sizeof(State_1D));
  print1 ("  sizeof (GRID)       = %d\n", sizeof(Grid));
  print1 ("  sizeof (TIME_STEP)  = %d\n", sizeof(Time_Step));
  print1 ("  sizeof (OUTPUT)     = %d\n", sizeof(Output));
  print1 ("  sizeof (RUNTIME)    = %d\n", sizeof(Runtime));
  print1 ("  sizeof (RESTART)    = %d\n", sizeof(Restart));
  print1 ("  sizeof (RGB)        = %d\n", sizeof(RGB));
  print1 ("  sizeof (IMAGE)      = %d\n", sizeof(Image));
  print1 ("  sizeof (FLOAT_VECT) = %d\n", sizeof(Float_Vect));
  print1 ("  sizeof (INDEX)      = %d\n", sizeof(Index));
  print1 ("  sizeof (RBOX)       = %d\n", sizeof(RBox));
*/

/* -- initialize members of Time_Step structure -- */

  Dts.cmax     = ARRAY_1D(NMAX_POINT, double);
  Dts.inv_dta  = 0.0;
  Dts.inv_dtp  = 0.5e-38; 
  Dts.dt_cool  = 1.e38;
  Dts.cfl      = ini.cfl;
  Dts.cfl_par  = ini.cfl_par;
  Dts.rmax_par = ini.rmax_par;
  Dts.Nsts     = Dts.Nrkc = 0;
  
  Solver = SetSolver (ini.solv_type);
  
  time (&tbeg);
  g_stepNumber = 0;
  
/* --------------------------------------------------------
    Check if restart is necessary. 
    If not, write initial condition to disk.
   ------------------------------------------------------- */
   
  if (cmd_line.restart == YES) {
    RestartFromFile (&ini, cmd_line.nrestart, DBL_OUTPUT, grd);
  }else if (cmd_line.h5restart == YES){
    RestartFromFile (&ini, cmd_line.nrestart, DBL_H5_OUTPUT, grd);
  }else if (cmd_line.write){
    CheckForOutput (&data, &ini, grd);
    CheckForAnalysis (&data, &ini, grd);
    #ifdef USE_ASYNC_IO
     Async_EndWriteData (&ini);
    #endif
  }

  print1 ("> Starting computation... \n\n");

/* =====================================================================
          M A I N      L O O P      S T A R T S      H E R E
   ===================================================================== */

#ifndef USE_ASYNC_IO  /* -- Standard loop, don't use Asynchrouns I/O -- */

  while (!last_step){

  /* ------------------------------------------------------
      Check if this is the last integration step:
      - final tstop has been reached: adjust time step 
      - or max number of steps has been reached
     ------------------------------------------------------ */

    if ((g_time + g_dt) >= ini.tstop*(1.0 - 1.e-8)) {
      g_dt   = (ini.tstop - g_time);
      last_step = 1;
    }
    if (g_stepNumber == cmd_line.maxsteps && cmd_line.maxsteps > 0) {
      last_step = 1;
    }

  /* ------------------------------------------------------
                Dump log information
     ------------------------------------------------------ */

    if (g_stepNumber%ini.log_freq == 0) {
      print1 ("step:%d ; t = %10.4e ; dt = %10.4e ; %d %% ; [%f, %d",
               g_stepNumber, g_time, g_dt, (int)(100.0*g_time/ini.tstop), 
               g_maxMach, g_maxRiemannIter);
/*      if (g_maxRootIter > 0) print1 (", root it. # = %d",g_maxRootIter);  */
      #if (PARABOLIC_FLUX & SUPER_TIME_STEPPING)
       print1 (", Nsts = %d",Dts.Nsts);
      #endif
      #if (PARABOLIC_FLUX & RK_CHEBYSHEV)
       print1 (", Nrkc = %d",Dts.Nrkc);
      #endif
      print1 ("]\n");      
    }

  /* ------------------------------------------------------
       check if it's time to write or perform analysis
     ------------------------------------------------------ */

    if (!first_step && !last_step && cmd_line.write) {
      CheckForOutput  (&data, &ini, grd);
      CheckForAnalysis(&data, &ini, grd);
    }

  /* ------------------------------------------------------
      Advance solution array by a single time step
      g_dt = dt(n)
     ------------------------------------------------------ */

    if (cmd_line.jet != -1) SetJetDomain (&data, cmd_line.jet, ini.log_freq, grd); 
    err = Integrate (&data, Solver, &Dts, grd);
    if (cmd_line.jet != -1) UnsetJetDomain (&data, cmd_line.jet, grd); 

  /* ------------------------------------------------------
       Integration didn't go through. Step must
       be redone from previously saved solution.
     ------------------------------------------------------ */
/*
    if (err != 0){
      print1 ("! Step failed. Re-trying\n");
      zones with problems must be tagged with MINMOD_FLAG and HLL_FLAG
      time step should be halved
      GET_SOL(&data);
    }
*/

  /* ------------------------------------------------------
      Increment time, t(n+1) = t(n) + dt(n)
     ------------------------------------------------------ */

    g_time += g_dt;

  /* ------------------------------------------------------
      Show the time step ratios between the actual g_dt
      and the advection, diffusion and cooling time scales.
     ------------------------------------------------------ */

    #if SHOW_TIME_STEPS == YES
     if (g_stepNumber%ini.log_freq == 0) {
       double cg, dta, dtp, dtc;
       dta = 1.0/Dts.inv_dta;
       dtp = 0.5/Dts.inv_dtp;
       dtc = Dts.dt_cool;
       #ifdef PARALLEL
        MPI_Allreduce (&dta, &cg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dta = cg;

        MPI_Allreduce (&dtp, &cg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dtp = cg;

        MPI_Allreduce (&dtc, &cg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dtc = cg;
       #endif
       /*
   print1 ("[dt/dt(adv) = %10.4e, dt/dt(par) = %10.4e, dt/dt(cool) = %10.4e]\n",
                g_dt/dta, g_dt/dtp, g_dt/dtc);
                */
       print1 ("  dt(adv)  = cfl x %10.4e;\n",dta);
       print1 ("  dt(par)  = cfl x %10.4e;\n",dtp);
       print1 ("  dt(cool) =       %10.4e;\n",dtc);
     }
    #endif

  /* ------------------------------------------------------
      Get next time step dt(n+1).
      Do it every two steps if cooling or dimensional
      splitting are used.
     ------------------------------------------------------ */

    #if (COOLING == NO) && ((DIMENSIONS == 1) || (DIMENSIONAL_SPLITTING == NO))
      g_dt = NextTimeStep(&Dts, &ini, grd);
    #else
     if (g_stepNumber%2 == 1) g_dt = NextTimeStep(&Dts, &ini, grd);
    #endif

  /* ------------------------------------------------------
          Global MPI reduction operations
     ------------------------------------------------------ */
  
    #ifdef PARALLEL
     MPI_Allreduce (&g_maxMach, &scrh, 1, 
                    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     g_maxMach = scrh;

     MPI_Allreduce (&g_maxRiemannIter, &nv, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
     g_maxRiemannIter = nv;
    #endif

    g_stepNumber++;
    
    first_step = 0;
  }

#else /* Use Asynchrounous I/O */

  while (!last_step){

  /* ------------------------------------------------------
      Check if this is the last integration step:
      - final tstop has been reached: adjust time step 
      - or max number of steps has been reached
     ------------------------------------------------------ */

    if ((g_time + g_dt) >= ini.tstop*(1.0 - 1.e-8)) {
      g_dt   = (ini.tstop - g_time);
      last_step = 1;
    }
    if (g_stepNumber == cmd_line.maxsteps && cmd_line.maxsteps > 0) {
      last_step = 1;
    }

  /* ------------------------------------------------------
       check if it's time to write or perform analysis
     ------------------------------------------------------ */

    if (!first_step && !last_step && cmd_line.write) {
      CheckForOutput  (&data, &ini, grd);
      CheckForAnalysis(&data, &ini, grd);
    }

  /* ------------------------------------------------------
      Show the time step ratios between the actual g_dt
      and the advection, diffusion and cooling time scales.
     ------------------------------------------------------ */

    #if SHOW_TIME_STEPS == YES
     if (!first_step && g_stepNumber%ini.log_freq == 0) {
       double cg, dta, dtp, dtc;
       dta = 1.0/Dts.inv_dta;
       dtp = 0.5/Dts.inv_dtp;
       dtc = Dts.dt_cool;
       #ifdef PARALLEL
        MPI_Allreduce (&dta, &cg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dta = cg;

        MPI_Allreduce (&dtp, &cg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dtp = cg;

        MPI_Allreduce (&dtc, &cg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dtc = cg;
       #endif
       print1 ("\t[dt/dta = %10.4e, dt/dtp = %10.4e, dt/dtc = %10.4e \n",
                g_dt/dta, g_dt/dtp, g_dt/dtc);
     }
    #endif

  /* ------------------------------------------------------
          Global MPI reduction operations
     ------------------------------------------------------ */
  
    #ifdef PARALLEL
     MPI_Allreduce (&g_maxMach, &scrh, 1, 
                    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     g_maxMach = scrh;

     MPI_Allreduce (&g_maxRiemannIter, &nv, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
     g_maxRiemannIter = nv;
    #endif

  /* ------------------------------------------------------
             Finish writing using Async I/O
     ------------------------------------------------------ */

    #ifdef USE_ASYNC_IO
     Async_EndWriteData (&ini);
    #endif

  /* ------------------------------------------------------
                Dump log information
     ------------------------------------------------------ */

    if (g_stepNumber%ini.log_freq == 0) {
      print1 ("step:%d ; t = %10.4e ; dt = %10.4e ; %d %% ; [%f, %d",
               g_stepNumber, g_time, g_dt, (int)(100.0*g_time/ini.tstop), 
               g_maxMach, g_maxRiemannIter);
      #if (PARABOLIC_FLUX & SUPER_TIME_STEPPING)
       print1 (", Nsts = %d",Dts.Nsts);
      #endif
      #if (PARABOLIC_FLUX & RK_CHEBYSHEV)
       print1 (", Nrkc = %d",Dts.Nrkc);
      #endif
      print1 ("]\n");      
    }
    
  /* ------------------------------------------------------
      Advance solution array by a single time step
      g_dt = dt(n)
     ------------------------------------------------------ */

    if (cmd_line.jet != -1) SetJetDomain (&data, cmd_line.jet, ini.log_freq, grd); 
    err = Integrate (&data, Solver, &Dts, grd);
    if (cmd_line.jet != -1) UnsetJetDomain (&data, cmd_line.jet, grd); 

  /* ------------------------------------------------------
       Integration didn't go through. Step must
       be redone from previously saved solution.
     ------------------------------------------------------ */
/*
    if (err != 0){
      print1 ("! Step failed. Re-trying\n");
      zones with problems must be tagged with MINMOD_FLAG and HLL_FLAG
      time step should be halved
      GET_SOL(&data);
    }
*/
  /* ------------------------------------------------------
      Increment time, t(n+1) = t(n) + dt(n)
     ------------------------------------------------------ */

    g_time += g_dt;

  /* ------------------------------------------------------
                Get next time step dt(n+1)
     ------------------------------------------------------ */

    g_dt = NextTimeStep(&Dts, &ini, grd);
    g_stepNumber++;
    first_step = 0;
  }
#endif /* USE_ASYNC_IO */

/* =====================================================================
          M A I N       L O O P      E N D S       H E R E 
   ===================================================================== */

  if (cmd_line.write){
    CheckForOutput (&data, &ini, grd);
    CheckForAnalysis (&data, &ini, grd);
    #ifdef USE_ASYNC_IO
     Async_EndWriteData (&ini);
    #endif
  }

  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
   print1  ("\n> Total allocated memory  %6.2f Mb (proc #%d)\n",
             (float)g_usedMemory/1.e6,prank);
   MPI_Barrier (MPI_COMM_WORLD);
  #else
   print1  ("\n> Total allocated memory  %6.2f Mb\n",(float)g_usedMemory/1.e6);
  #endif

  time(&tend);
  g_dt = difftime(tend, tbeg);
  print1("> Elapsed time             %s\n", TotalExecutionTime(g_dt));
  print1("> Average time/step       %10.2e  (sec)  \n", 
          difftime(tend,tbeg)/(double)g_stepNumber);
  print1("> Local time                %s",asctime(localtime(&tend)));
  print1("> Done\n");

  FreeArray4D ((void *) data.Vc);
  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
   AL_Finalize ();
  #endif

  return (0);
}
#undef SHOW_TIME_STEPS
/* ******************************************************************** */
int Integrate (Data *d, Riemann_Solver *Solver, Time_Step *Dts, Grid *grid)
/*!
 * Advance equations by a single time-step.

 * \param  d      pointer to PLUTO Data structure;
 * \param  Solver pointer to a Riemann solver function;
 * \param  Dts    pointer to time Step structure;
 * \param  grid   pointer to grid structure.
 * 
 * \return An integer giving success / failure (development).
 * 
 ********************************************************************** */
{
  int idim, err = 0;
  int i,j,k;
  
  g_maxMach = 0.0;
  g_maxRiemannIter = 0;
  g_maxRootIter    = 0;

/* -------------------------------------------------------
    Initialize max propagation speed in Dedner's approach
   ------------------------------------------------------- */

  #ifdef GLM_MHD  /* -- initialize glm_ch -- */
   GLM_Init (d, Dts, grid);   
   GLM_Source (d->Vc, 0.5*g_dt, grid);
  #endif

  /* ---------------------------------------------
        perform Strang Splitting on directions 
        (if necessary) and sources 
     --------------------------------------------- */

  TOT_LOOP(k,j,i) d->flag[k][j][i] = 0;

  #ifdef FARGO
   FARGO_ComputeVelocity(d, grid);
  #endif
  if ((g_stepNumber%2) == 0){
    g_operatorStep = HYPERBOLIC_STEP;
    #if DIMENSIONAL_SPLITTING == YES
     for (g_dir = 0; g_dir < DIMENSIONS; g_dir++){
       if (AdvanceStep (d, Solver, Dts, grid) != 0) return (1);
     }
    #else
     if (AdvanceStep (d, Solver, Dts, grid) != 0) return(1);
    #endif
    g_operatorStep = PARABOLIC_STEP;
    SplitSource (d, g_dt, Dts, grid);
  }else{
    g_operatorStep = PARABOLIC_STEP;
    SplitSource (d, g_dt, Dts, grid);
    g_operatorStep = HYPERBOLIC_STEP;
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
double NextTimeStep (Time_Step *Dts, Runtime *ini, Grid *grid)
/*!
 * Compute and return the time step for the next time level
 * using the information from the previous integration
 * (Dts->inv_dta and Dts->inv_dp).
 *
 * \param [in] Dts    pointer to the Time_Step structure
 * \param [in] ini    pointer to the Runtime structure
 * \param [in] grid   pointer to array of Grid structures
 *
 * \return The time step for next time level
 *********************************************************************** */
{
  int idim;
  double dt_adv, dt_par, dtnext;
  double dxmin;
  double xloc, xglob;

/* ---------------------------------------------------
   1. Take the maximum of inv_dt across all processors
   --------------------------------------------------- */

  #ifdef PARALLEL
   xloc = Dts->inv_dta;
   MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   Dts->inv_dta = xglob;
   #if (PARABOLIC_FLUX != NO)
    xloc = Dts->inv_dtp;
    MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    Dts->inv_dtp = xglob;
   #endif
   #if COOLING != NO
    xloc = Dts->dt_cool;
    MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    Dts->dt_cool = xglob;
   #endif
  #endif

/* ----------------------------------
   2. Compute time step
   ---------------------------------- */

  #if (PARABOLIC_FLUX & EXPLICIT)
   dt_adv  = 1.0/(Dts->inv_dta + 2.0*Dts->inv_dtp);
  #else
   dt_adv  = 1.0/Dts->inv_dta;
  #endif
  dt_adv *= ini->cfl;
  dtnext  = dt_adv;

/* -------------------------------------------------------
   3. Maximum propagation speed for the local processor.
      Global glm_ch will be computed later in GLM_Init.
   ------------------------------------------------------- */

  #ifdef GLM_MHD
   dxmin = grid[IDIR].dl_min;
   for (idim = 1; idim < DIMENSIONS; idim++){ /*  Min cell length   */
     dxmin = MIN(dxmin, grid[idim].dl_min);
   }
   glm_ch = ini->cfl*dxmin/dtnext;
  #endif

/* ---------------------------------------------------------
   4. With STS, the ratio between advection (full) and 
      parabolic time steps should not exceed ini->rmax_par.
   --------------------------------------------------------- */
      
  #if (PARABOLIC_FLUX & SUPER_TIME_STEPPING) || (PARABOLIC_FLUX & RK_CHEBYSHEV)
   dt_par  = ini->cfl_par/(2.0*Dts->inv_dtp);
   dtnext *= MIN(1.0, ini->rmax_par/(dt_adv/dt_par));
  #endif

/* ----------------------------------
   5. Compute Cooling time step
   ---------------------------------- */

  #if COOLING != NO
   dtnext = MIN(dtnext, Dts->dt_cool);
  #endif
   
/* --------------------------------------------------------------
    6. Allow time step to vary at most by a factor 
       ini->cfl_max_var.
       Quit if dt gets too small, issue a warning if first_dt has
       been overestimated.
   -------------------------------------------------------------- */

  dtnext = MIN(dtnext, ini->cfl_max_var*g_dt);

  if (dtnext < ini->first_dt*1.e-9){
    print1 ("! NextTimeStep(): dt is too small (%12.6e). Cannot continue.\n", dtnext);
    QUIT_PLUTO(1);
  }

  if (g_stepNumber <= 1 && (ini->first_dt > dtnext/ini->cfl)){
    print1 ("! NextTimeStep(): initial dt exceeds stability limit\n");
  }

/* --------------------------------------------
   7. Reset time step coefficients
   -------------------------------------------- */

  DIM_LOOP(idim) Dts->cmax[idim] = 0.0;
  Dts->inv_dta = 0.0;
  Dts->inv_dtp = 0.5e-38;
  Dts->dt_cool = 1.e38;

  return(dtnext);
}

/* ********************************************************************* */
void CheckForOutput (Data *d, Runtime *ini, Grid *grid)
/*!
 *  Check if file output has to be performed.
 *  
 *********************************************************************** */
{
  static int first_call = 1;
  int  n, check_dt, check_dn, check_dclock;
  int  restart_update, last_step;
  double t, tnext;
  Output *output;
  static time_t clock_beg[MAX_OUTPUT_TYPES], clock_end;
  static double tbeg[MAX_OUTPUT_TYPES], tend;
  double dclock;

  restart_update = 0;
  t     = g_time;
  tnext = t + g_dt;
  
  last_step = (fabs(t-ini->tstop) < 1.e-12 ? 1:0);

/* -- on first execution initialize
      current beginning time for all output types -- */

  if (first_call){
    #ifdef PARALLEL
     if (prank == 0){
       double tstart;
       tstart = MPI_Wtime();
       for (n = 0; n < MAX_OUTPUT_TYPES; n++) tbeg[n] = tstart;
     }
     MPI_Bcast(tbeg, MAX_OUTPUT_TYPES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #else
     for (n = 0; n < MAX_OUTPUT_TYPES; n++) time(clock_beg + n);
    #endif
  } 

/* -- get current time -- */

  #ifdef PARALLEL  
   if (prank == 0) tend = MPI_Wtime();
   MPI_Bcast(&tend, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #else
   time(&clock_end);
  #endif
  
/* -------------------------------------------------------
          start main loop on outputs
   ------------------------------------------------------- */
   
  for (n = 0; n < MAX_OUTPUT_TYPES; n++){
    output = ini->output + n;    
    check_dt = check_dn = check_dclock = 0; 

  /* -- check time interval in code units (dt) -- */

    if (output->dt > 0.0){
      check_dt = (int) (tnext/output->dt) - (int)(t/output->dt);
      check_dt = check_dt || g_stepNumber == 0 || last_step;
    }   

  /* -- check time interval in number of steps (dn) -- */

    if (output->dn > 0){
      check_dn = (g_stepNumber%output->dn) == 0;
      check_dn = check_dn || g_stepNumber == 0 || last_step;
    }

  /* -- check time interval in clock time (dclock) -- */

    if (output->dclock > 0.0){
      #ifdef PARALLEL
       dclock = tend - tbeg[n];
      #else
       dclock = difftime(clock_end, clock_beg[n]);
      #endif
      if (dclock >= output->dclock) {
        check_dclock = 1;
        #ifdef PARALLEL
         tbeg[n] = tend;
        #else
         time(clock_beg + n);
        #endif
      }else{ 
        check_dclock = 0;
      }
      check_dclock = check_dclock || g_stepNumber == 0 || last_step;
    }

  /* -- if any of the previous is true dump data to disk -- */

    if (check_dt || check_dn || check_dclock) { 

      #ifdef USE_ASYNC_IO
       if (!strcmp(output->mode,"single_file_async")){
         Async_BegWriteData (d, output, grid);
       }else{
         WriteData(d, output, grid);
       }
      #else     
       WriteData(d, output, grid);
      #endif   

    /* ----------------------------------------------------------
        save the file number of the dbl and dbl.h5 output format
        for writing restart.out once we exit the loop.
       ---------------------------------------------------------- */

      if ((output->type == DBL_OUTPUT) ||
          (output->type == DBL_H5_OUTPUT)) restart_update = 1;
    }
  }

/* -------------------------------------------------------
    Dump restart information if required 

    Note that if both dbl and dbl.h5 formats are used,
    bookkeeping is done using dbl format.
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
