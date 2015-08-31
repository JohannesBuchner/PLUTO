/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Print useful information about the current computations.

  

  \author A. Mignone (mignone@ph.unito.it)
  \date   May 15, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
static void CheckConfig();

/* ********************************************************************* */
void ShowConfig (int argc, char *argv[], char *ini_file)
/*!
 *  Write a summary of the selected options
 *
      ___  __   __  ____________ 
     / _ \/ /  / / / /_  __/ __ \
    / ___/ /__/ /_/ / / / / /_/ /
   / /  /____/\____/ /_/  \____/ 
   ============================== 
                                  
 *
 * 
 *                        
 *
 *********************************************************************** */
{
  int  n;
  FILE *fp;
  time_t time_now;
  char  str1[128], str2[128], str3[128], sline[512];

  CheckConfig();

  print1 ("\n");
  print1("   ___  __   __  ____________   \n");
  print1("  / _ \\/ /  / / / /_  __/ __ \\ \n");
  print1(" / ___/ /__/ /_/ / / / / /_/ /  \n");
  print1("/_/  /____/\\____/ /_/  \\____/   \n");
  print1("=============================    v. %s  \n", PLUTO_VERSION);
  
  print1 ("\n> System:\n\n");

  if ( (fp = fopen("sysconf.out","r")) != NULL){

    while (fscanf (fp, "%s %s %s\n", str1, str2, str3) != EOF) {
      if (!strcmp(str1,"USER")) 
        print1 ("  %s:             %s\n",str1, str3);
      else if (!strcmp(str1,"WORKING_g_dir"))
        print1 ("  %s:      %s\n",str1, str3);
      else if (!strcmp(str1,"SYSTEM_NAME"))
        print1 ("  %s:      %s\n",str1, str3);
      else if (!strcmp(str1,"NODE_NAME"))
        print1 ("  %s:        %s\n",str1, str3);
      else if (!strcmp(str1,"ARCH"))
        print1 ("  %s:             %s\n",str1, str3);
      else if (!strcmp(str1,"BYTE_ORDER"))
        print1 ("  %s:       %s\n\n",str1, str3);
    }
    fclose(fp);

  }else{
    print1 ("! sysconf.out file not found \n\n");
  }
 
  time(&time_now);
  print1("> Local time:       %s\n",asctime(localtime(&time_now)));
      
  if (COMPONENTS < DIMENSIONS) {
    print1 ("Sorry, but the number of vector components can\n");
    print1 ("not be less than the dimension number.\n");
    print1 ("Please edit definitions.h and fix this.\n");
    QUIT_PLUTO(1);
  }

/* -- print command line arguments -- */

  print1 ("> Cmd line args:    ");
  for (n = 1; n < argc; n++) print1 ("%s ",argv[n]);
  print1 ("\n\n");

/* -- print problem configuration -- */

  print1 ("> Header configuration:\n\n");

  if (PHYSICS == ADVECTION) print1 ("  PHYSICS:          ADVECTION\n");
  if (PHYSICS == HD)        print1 ("  PHYSICS:          HD\n");
  if (PHYSICS == RHD)       print1 ("  PHYSICS:          RHD\n");
  if (PHYSICS == MHD)       print1 ("  PHYSICS:          MHD [div.B: ");
  if (PHYSICS == RMHD)      print1 ("  PHYSICS:          RMHD [div.B: ");
#if PHYSICS == MHD || PHYSICS == RMHD
  #if DIVB_CONTROL == NO
  print1 ("None]\n");
  #elif DIVB_CONTROL == EIGHT_WAVES
    print1 ("Powell's 8wave]\n");
  #elif DIVB_CONTROL == DIV_CLEANING
    #if GLM_EXTENDED == NO 
    print1 ("Divergence Cleaning (GLM)]\n");
    #elif GLM_EXTENDED == YES
    print1 ("Divergence Cleaning (Extended GLM)]\n");
    #endif
    
  #elif DIVB_CONTROL == CONSTRAINED_TRANSPORT
    print1 ("CT/");
    #if CT_EMF_AVERAGE == ARITHMETIC
    print1 ("Ar. average]\n");
    #elif CT_EMF_AVERAGE == UCT_CONTACT
    print1 ("UCT_CONTACT]\n");
    #elif CT_EMF_AVERAGE == UCT0
    print1 ("UCT0]\n");
    #elif CT_EMF_AVERAGE == UCT_HLL
    print1 ("UCT_HLL]\n");
    #endif
  #endif
#endif

  print1 ("  DIMENSIONS:       %d\n", DIMENSIONS);
  print1 ("  COMPONENTS:       %d\n", COMPONENTS);

  print1 ("  GEOMETRY:         ");
  if (GEOMETRY == CARTESIAN)    print1 ("Cartesian\n");
  if (GEOMETRY == CYLINDRICAL)  print1 ("Cylindrical\n");
  if (GEOMETRY == POLAR)        print1 ("Polar\n");
  if (GEOMETRY == SPHERICAL)    print1 ("Spherical\n");

  print1 ("  BODY_FORCE:       ");
  print1 (BODY_FORCE == NO ? "NO\n":"EXPLICIT\n");

  print1 ("  RECONSTRUCTION:   ");
#ifndef FINITE_DIFFERENCE
   if (RECONSTRUCTION == FLAT)          print1 ("Flat");
   if (RECONSTRUCTION == LINEAR)        print1 ("Linear TVD");
   if (RECONSTRUCTION == LINEAR_MULTID) print1 ("Linear_Multid");
   if (RECONSTRUCTION == LimO3)         print1 ("LimO3");
   if (RECONSTRUCTION == WENO3)         print1 ("WENO 3rd order");
   if (RECONSTRUCTION == PARABOLIC)     print1 ("Parabolic");
 #ifdef CHAR_LIMITING
   if (CHAR_LIMITING == YES) print1 (" (Characteristic lim)\n");
   else                      print1 (" (Primitive lim)\n");
 #endif
#endif

#ifdef FINITE_DIFFERENCE
  if (RECONSTRUCTION == LIMO3_FD)     print1 ("LimO3 (FD), 3rd order\n");
  if (RECONSTRUCTION == WENO3_FD)     print1 ("WENO3 (FD), 3rd order\n");
  if (RECONSTRUCTION == WENOZ_FD)     print1 ("WENOZ (FD) 5th order\n");
  if (RECONSTRUCTION == MP5_FD)       print1 ("MP5 (FD), 5th order\n");
#endif

  print1 ("  TRACERS:          %d\n", NTRACER);
  print1 ("  VARIABLES:        %d\n", NVAR);
  print1 ("  ENTROPY_SWITCH:   %s\n",(ENTROPY_SWITCH != NO ? "ENABLED":"NO"));
#if PHYSICS == MHD 
  print1 ("  BACKGROUND_FIELD: %s\n",(BACKGROUND_FIELD == YES ? "YES":"NO"));
#endif

  print1 ("  LOADED MODULES:\n");
  #if PHYSICS == MHD
   #ifdef SHEARINGBOX
    print1 ("\n  o [SHEARINGBOX]\n");
    print1 ("     - Order:             %d\n", SB_ORDER);
    print1 ("     - Sym Hydro Flux:    %s\n", 
             (SB_SYMMETRIZE_HYDRO == YES ? "YES":"NO"));
    print1 ("     - Sym Ey:            %s\n", 
             (SB_SYMMETRIZE_EY == YES ? "YES":"NO"));
    print1 ("     - Sym Ez:            %s\n", 
             (SB_SYMMETRIZE_EZ == YES ? "YES":"NO"));
    print1 ("     - Force EMF periods: %s\n", 
             (SB_FORCE_EMF_PERIODS == YES ? "YES":"NO"));
   #endif
  #endif
  #ifdef FARGO
   print1 ("\n  o [FARGO]\n");
   print1 ("     - Order:         %d\n", FARGO_ORDER);
   print1 ("     - Average Speed: %s\n", 
            (FARGO_AVERAGE_VELOCITY == YES ? "YES":"NO"));
   print1 ("     - Av. Frequency: %d\n", FARGO_NSTEP_AVERAGE);

  #endif
  print1 ("\n");

  print1 ("  ROTATION:         ");
  print1(ROTATING_FRAME == YES ? "YES\n":"NO\n");

  print1 ("  EOS:              ");
  if      (EOS == IDEAL)        print1 ("Ideal\n");
  else if (EOS == PVTE_LAW)     print1 ("PVTE_LAW\n");
  else if (EOS == BAROTROPIC)   print1 ("Barotropic\n");
  else if (EOS == ISOTHERMAL)   print1 ("Isothermal\n");
  else if (EOS == TAUB)         print1 ("Taub - TM\n");
  else                          print1 ("None\n");

  print1 ("  TIME INTEGRATOR:  ");
  if (TIME_STEPPING == EULER)            print1 ("Euler\n");
  if (TIME_STEPPING == RK2)              print1 ("Runga-Kutta II\n");
  if (TIME_STEPPING == RK3)              print1 ("Runga_Kutta III\n");
  if (TIME_STEPPING == CHARACTERISTIC_TRACING)
                                         print1 ("Characteristic Tracing\n");
  if (TIME_STEPPING == HANCOCK)          print1 ("Hancock\n");

  print1 ("  DIM. SPLITTING:   ");
  if (DIMENSIONAL_SPLITTING == YES)  print1 ("Yes\n");
  else                               print1 ("No\n");
  

  #if PARABOLIC_FLUX != NO
   print1 ("  DIFFUSION TERMS:");
   #if (RESISTIVITY == EXPLICIT) 
    print1 ("  Resistivity  [EXPLICIT]\n");  
   #elif (RESISTIVITY == SUPER_TIME_STEPPING)
    print1 ("  Resistivity  [STS]\n");  
   #endif

   #if (THERMAL_CONDUCTION == EXPLICIT) 
    print1 ("  Thermal Conduction [EXPLICIT]\n");  
   #elif (THERMAL_CONDUCTION == SUPER_TIME_STEPPING)
    print1 ("  Thermal Conduction [STS]\n");  
   #endif

   #if (VISCOSITY == EXPLICIT) 
    print1 ("  Viscosity [EXPLICIT]\n");  
   #elif (VISCOSITY == SUPER_TIME_STEPPING)
    print1 ("  Viscosity [STS]\n");  
   #endif
  #endif

  print1 ("\n");

/* -----------------------------------------------------
    Print runtime configuration info (definitions.h 
    and from pluto.ini)
   ----------------------------------------------------- */
/*   
  print1 ("> Header file configuration (definitions.h):\n\n");
  print1 ("  +----------------------------------------------------------\n");
  fp = fopen("definitions.h","r");
  while ( fgets(sline, 512, fp) != NULL ) {
    print1 ("  | %s",sline);
  }
  fclose(fp);
  print1 ("  +---------------------------------------------------------\n\n");
*/
  print1 ("> Runtime configuration (%s):\n\n", ini_file);
  print1 ("  +----------------------------------------------------------\n");
  fp = fopen(ini_file,"r");
  while ( fgets(sline, 512, fp) != NULL ) {
    print1 ("  | %s",sline);
  }
  fclose(fp);
  print1 ("  +---------------------------------------------------------\n");


}

/* ********************************************************************* */
void ShowUnits ()
/*!
 *  Show units when cooling is enabled.
 *
 *
 *********************************************************************** */
{

#if COOLING != NO
  print1 ("> Cooling Module:    ");
  if (COOLING == SNEq)  print1 (" SNEq\n");
  if (COOLING == MINEq) print1 (" MINEq\n");
  if (COOLING == TABULATED) print1 (" TABULATED\n");
  if (COOLING == H2_COOL) print1 (" H2_COOL \n");
#endif

  print1 ("> Normalization Units:\n\n");
  print1 ("  [Density]:      %8.3e (gr/cm^3), %8.3e (1/cm^3)\n",
          UNIT_DENSITY,UNIT_DENSITY/CONST_mp);
  print1 ("  [Pressure]:     %8.3e (dyne/cm^2)\n",
          UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
  print1 ("  [Velocity]:     %8.3e (cm/s)\n",UNIT_VELOCITY);
  print1 ("  [Length]:       %8.3e (cm)\n",UNIT_LENGTH);
  print1 ("  [Temperature]:  %8.3e X (p/rho*mu) (K)\n",KELVIN);
  print1 ("  [Time]:         %8.3e (sec), %8.3e (yrs) \n",
       UNIT_LENGTH/UNIT_VELOCITY, UNIT_LENGTH/UNIT_VELOCITY/86400./365.);
#if PHYSICS == MHD || PHYSICS == RMHD
  print1 ("  [Mag Field]:    %8.3e (Gauss)\n",
           UNIT_VELOCITY*sqrt(4.0*CONST_PI*UNIT_DENSITY));
#endif

  print1 (" \n");
}

/* ********************************************************************* */
void ShowDomainDecomposition(int nprocs, Grid *GXYZ)
/*!
 * Show the parallel domain decomposition by having each processor print
 * its own computational domain.
 * This is activated with the -show-dec command line argument.
 * It may be long for several thousand processors. 
 *
 * \param [in] nprocs   the total number of processors
 * \param [in] GXYZ     a pointer to an array of grid structures
 *********************************************************************** */
{
  int p;
  int ib, ie, jb, je, kb, ke;
  int nxp, nyp, nzp;
  int nghx, nghy, nghz;

  static int *ib_proc, *ie_proc;
  static int *jb_proc, *je_proc;
  static int *kb_proc, *ke_proc;

  double xb, xe, yb, ye, zb, ze;

  double *xb_proc, *xe_proc;
  double *yb_proc, *ye_proc;
  double *zb_proc, *ze_proc;

  Grid *Gx, *Gy, *Gz;
  
  Gx = GXYZ;
  Gy = GXYZ + 1;
  Gz = GXYZ + 2;

/* ---- Allocate memory ---- */

  ib_proc = ARRAY_1D(nprocs, int); ie_proc = ARRAY_1D(nprocs, int);
  jb_proc = ARRAY_1D(nprocs, int); je_proc = ARRAY_1D(nprocs, int);
  kb_proc = ARRAY_1D(nprocs, int); ke_proc = ARRAY_1D(nprocs, int);
  xb_proc = ARRAY_1D(nprocs, double); xe_proc = ARRAY_1D(nprocs, double);
  yb_proc = ARRAY_1D(nprocs, double); ye_proc = ARRAY_1D(nprocs, double);
  zb_proc = ARRAY_1D(nprocs, double); ze_proc = ARRAY_1D(nprocs, double);

#ifdef PARALLEL  
  nxp = Gx->np_tot;
  nyp = Gy->np_tot;
  nzp = Gz->np_tot;

/* -- Local beg and end indices -- */

  ib = nghx = Gx->nghost; ie = ib + Gx->np_int - 1;  
  jb = nghy = Gy->nghost; je = jb + Gy->np_int - 1;
  kb = nghz = Gz->nghost; ke = kb + Gz->np_int - 1;

/* -- Leftmost and rightmost processor coordinates -- */

  xb = Gx->xl[ib]; xe = Gx->xr[ie];
  yb = Gy->xl[jb]; ye = Gy->xr[je];
  zb = Gz->xl[kb]; ze = Gz->xr[ke];

  D_EXPAND(
    MPI_Gather (&xb, 1, MPI_DOUBLE, xb_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&xe, 1, MPI_DOUBLE, xe_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  ,
    
    MPI_Gather (&yb, 1, MPI_DOUBLE, yb_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&ye, 1, MPI_DOUBLE, ye_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  ,
    
    MPI_Gather (&zb, 1, MPI_DOUBLE, zb_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather (&ze, 1, MPI_DOUBLE, ze_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  )

/* -- Global beg and end indices -- */

  ib  = Gx->beg; ie += Gx->beg - nghx;
  jb  = Gy->beg; je += Gy->beg - nghy;
  kb  = Gz->beg; ke += Gz->beg - nghz;

  D_EXPAND(
    MPI_Gather (&ib, 1, MPI_INT, ib_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather (&ie, 1, MPI_INT, ie_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);  ,
    
    MPI_Gather (&jb, 1, MPI_INT, jb_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather (&je, 1, MPI_INT, je_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);  ,
    
    MPI_Gather (&kb, 1, MPI_INT, kb_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather (&ke, 1, MPI_INT, ke_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);
  )

  print1 ("> Domain Decomposition (%d procs):\n\n", nprocs);

  for (p = 0; p < nprocs; p++){ 
    D_EXPAND(
      print1 ("  - Proc # %d, X1: [%f, %f], i: [%d, %d], NX1: %d\n",
              p, xb_proc[p], xe_proc[p], 
                 ib_proc[p], ie_proc[p], 
                 ie_proc[p]-ib_proc[p]+1);    ,
      print1 ("              X2: [%f, %f], j: [%d, %d]; NX2: %d\n",
                 yb_proc[p], ye_proc[p], 
                 jb_proc[p], je_proc[p],
                 je_proc[p]-jb_proc[p]+1);    ,
      print1 ("              X3: [%f, %f], k: [%d, %d], NX3: %d\n\n",  
                 zb_proc[p], ze_proc[p], 
                 kb_proc[p], ke_proc[p],
                 ke_proc[p]-kb_proc[p]+1);
    )
  }

  MPI_Barrier (MPI_COMM_WORLD);
#endif

/* ---- Free memory ---- */

  FreeArray1D((void *) ib_proc); FreeArray1D((void *) ie_proc);
  FreeArray1D((void *) jb_proc); FreeArray1D((void *) je_proc);
  FreeArray1D((void *) kb_proc); FreeArray1D((void *) ke_proc);
  FreeArray1D((void *) xb_proc); FreeArray1D((void *) xe_proc);
  FreeArray1D((void *) yb_proc); FreeArray1D((void *) ye_proc);
  FreeArray1D((void *) zb_proc); FreeArray1D((void *) ze_proc);

  print1 ("\n");
}

/* ********************************************************************* */
void CheckConfig()
/*
 *
 *
 * Check if the selected configuration is 
 * allowed.
 *
 *
 *********************************************************************** */
{
  #if DIMENSIONS == 3 

   #if GEOMETRY  == CYLINDRICAL 
    print1 ("\n! Cylindrical coordinates are only 2D.\n");
    print1 ("! Use polar instead.\n");
    QUIT_PLUTO(1);
   #endif

   #if GEOMETRY == SPHERICAL 
    #if (TIME_STEPPING == HANCOCK) || (TIME_STEPPING == CHARACTERISTIC_TRACING)
     print1 ("\n ! Spherical 3D only works with RK integrators\n");
     QUIT_PLUTO(1);
    #endif
   #endif

  #endif

  #if DIMENSIONAL_SPLITTING == NO && DIMENSIONS == 1
   #ifndef CH_SPACEDIM
    print1 ("! CheckConfig(): Cannot integrate a 1-D problem with an unsplit method \n");
    QUIT_PLUTO(1);
   #endif
  #endif

#if (defined STAGGERED_MHD) && (DIMENSIONAL_SPLITTING == YES)
  print1 ("! CheckConfig(): CT requires dimensional unsplit scheme.\n");
  QUIT_PLUTO(1); 
#endif


/*
  #if GEOMETRY == SPHERICAL || GEOMETRY == POLAR
   #if TIME_STEPPING != RK2 || TIME_STEPPING != RK3
    print1 (" ! Spherical and Polar geometries work with RK integrators\");
    QUIT_PLUTO(1);
   #endif
  #endif
*/

   
}
