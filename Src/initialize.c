/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Initialize PLUTO.

  Initialize() performs a number of initialization tasks before 
  starting the main computation loop.\n
  More precisely, it completes the following sequence of steps:

  - parse command line options 
  - runtime initialization (i.e. pluto.ini)
  - parallel domain decomposition
  - grid generation
  - memory allocation
  - initial conditions
  - set output attributes

  The function GetDecompMode() sets the parallel domain decomposition 
  mode which can be equal to
  
  - AL_AUTO_DECOMP   default;
  - AL_USER_DECOMP   if the -dec n1 [n2] [n3] command line 
                     argument has been given;
                     In this case only, procs[] contains the number
                     of processors in the three directions;
  - AL_MPI_DECOMP   [todo]

  \author A. Mignone (mignone@ph.unito.it)
          B. Vaidya
  \date   May 21, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static int GetDecompMode (Cmd_Line *cmd_line, int procs[]);

/* ********************************************************************* */
void Initialize(int argc, char *argv[], Data *data, 
                Runtime *runtime, Grid *grid, Cmd_Line *cmd_line)
/*!
 * Initialize computational grid, domain decomposition and memory
 * allocation. Also, set initial conditions and output attributes.
 *
 * \param [in]     argc       the number of command-line argument passed to 
 *                            the code 
 * \param [in]     argv       the argument value as a 1D array of char
 * \param [in,out] data       a pointer to the main PLUTO data structure
 * \param [in,out] runtime    a pointer to the Runtime structure
 * \param [in]      grid      pointer to an array of Grid structures
 * \param [in]     cmd_line   pointer to the Cmd_Line structure
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int  nprocs, decomp_mode;
  int  i, j, k, idim, nv;
  int  nx, ny, nz, nghost, status;
  int  gsize[DIMENSIONS], lsize[DIMENSIONS];
  int  beg[DIMENSIONS], end[DIMENSIONS];
  int  gbeg[DIMENSIONS], gend[DIMENSIONS];
  int  lbeg[DIMENSIONS], lend[DIMENSIONS];
  int  is_gbeg[DIMENSIONS], is_gend[DIMENSIONS];
  int  ghosts[DIMENSIONS];
  int  periods[DIMENSIONS];
  int  pardim[DIMENSIONS], stagdim[DIMENSIONS];
  int  procs[DIMENSIONS];
  char ini_file[128];
  double scrh, dxmin[3], dxming[3];
  Output *output;
  #ifdef PARALLEL
   MPI_Datatype rgb_type;
   MPI_Datatype Float_Vect_type;
  #endif

/* -- set default input file name -- */

  sprintf (ini_file,"pluto.ini");

/* -- parse command line options -- */

  ParseCmdLineArgs (argc, argv, ini_file, cmd_line);

#ifdef PARALLEL

/* -- get number of processors -- */

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

/* -- read initialization file -- */

  if (prank == 0) RuntimeSetup (runtime, cmd_line, ini_file);
  MPI_Bcast (runtime,  sizeof (Runtime) , MPI_BYTE, 0, MPI_COMM_WORLD);

/* -- get number of ghost zones and set periodic boundaries -- */

  nghost = GetNghost();
  MPI_Allreduce (&nghost, &idim, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  nghost = idim;
   
  for (idim = 0; idim < DIMENSIONS; idim++) {
    gsize[idim]   = runtime->npoint[idim];
    ghosts[idim]  = nghost;
    periods[idim] =    (runtime->left_bound[idim] == PERIODIC ? 1:0)
                    || (runtime->left_bound[idim] == SHEARING ? 1:0);
    pardim[idim]  = cmd_line->parallel_dim[idim];
  }

/* -- find parallel decomposition mode and number of processors -- */

  decomp_mode = GetDecompMode(cmd_line, procs);

/* ---- double distributed array descriptor ---- */

/* SetDistributedArray (SZ, type, gsize, periods, stagdim); 
  return args: beg, end, lsize, lbeg, lend, gbeg, gend, is_gbeg, is_gend
*/

  AL_Sz_init (MPI_COMM_WORLD, &SZ);
  AL_Set_type (MPI_DOUBLE, 1, SZ);
  AL_Set_dimensions (DIMENSIONS, SZ);
  AL_Set_global_dim (gsize, SZ);
  AL_Set_ghosts (ghosts, SZ);
  AL_Set_periodic_dim (periods, SZ);
  AL_Set_parallel_dim (pardim, SZ);

  AL_Decompose (SZ, procs, decomp_mode);
  AL_Get_local_dim (SZ, lsize);
  AL_Get_bounds (SZ, beg, end, ghosts, AL_C_INDEXES);
  AL_Get_lbounds (SZ, lbeg, lend, ghosts, AL_C_INDEXES);
  AL_Get_gbounds (SZ, gbeg, gend, ghosts, AL_C_INDEXES);
  AL_Is_boundary (SZ, is_gbeg, is_gend);

/* ---- float distributed array descriptor ---- */

  AL_Sz_init (MPI_COMM_WORLD, &SZ_float);
  AL_Set_type (MPI_FLOAT, 1, SZ_float);
  AL_Set_dimensions (DIMENSIONS, SZ_float);
  AL_Set_global_dim (gsize, SZ_float);
  AL_Set_ghosts (ghosts, SZ_float);
  AL_Set_periodic_dim (periods, SZ_float);
  AL_Set_parallel_dim (pardim, SZ_float);

  AL_Decompose (SZ_float, procs, decomp_mode);
  AL_Get_local_dim (SZ_float, lsize);
  AL_Get_bounds (SZ_float, beg, end, ghosts, AL_C_INDEXES);
  AL_Get_lbounds (SZ_float, lbeg, lend, ghosts, AL_C_INDEXES);
  AL_Get_gbounds (SZ_float, gbeg, gend, ghosts, AL_C_INDEXES);
  AL_Is_boundary (SZ_float, is_gbeg, is_gend);

/* ---- char distributed array descriptor ---- */

  AL_Sz_init (MPI_COMM_WORLD, &SZ_char);
  AL_Set_type (MPI_CHAR, 1, SZ_char);
  AL_Set_dimensions (DIMENSIONS, SZ_char);
  AL_Set_global_dim (gsize, SZ_char);
  AL_Set_ghosts (ghosts, SZ_char);
  AL_Set_periodic_dim (periods, SZ_char);
  AL_Set_parallel_dim (pardim, SZ_char);

  AL_Decompose (SZ_char, procs, decomp_mode);
  AL_Get_local_dim (SZ_char, lsize);
  AL_Get_bounds (SZ_char, beg, end, ghosts, AL_C_INDEXES);
  AL_Get_lbounds (SZ_char, lbeg, lend, ghosts, AL_C_INDEXES);
  AL_Get_gbounds (SZ_char, gbeg, gend, ghosts, AL_C_INDEXES);
  AL_Is_boundary (SZ_char, is_gbeg, is_gend);

/* ---- Float_Vec distributed array descriptor ---- */

  MPI_Type_contiguous (3, MPI_FLOAT, &Float_Vect_type);
  MPI_Type_commit (&Float_Vect_type);
 
  AL_Sz_init (MPI_COMM_WORLD, &SZ_Float_Vect);
  AL_Set_type (MPI_FLOAT, 3, SZ_Float_Vect);
  AL_Set_dimensions (DIMENSIONS, SZ_Float_Vect);
  AL_Set_global_dim (gsize, SZ_Float_Vect);
  AL_Set_ghosts (ghosts, SZ_Float_Vect);
  AL_Set_periodic_dim (periods, SZ_Float_Vect);
  AL_Set_parallel_dim (pardim, SZ_Float_Vect);

  AL_Decompose (SZ_Float_Vect, procs, decomp_mode);
  AL_Get_local_dim (SZ_Float_Vect, lsize);
  AL_Get_bounds  (SZ_Float_Vect, beg, end, ghosts, AL_C_INDEXES);
  AL_Get_lbounds (SZ_Float_Vect, lbeg, lend, ghosts, AL_C_INDEXES);
  AL_Get_gbounds (SZ_Float_Vect, gbeg, gend, ghosts, AL_C_INDEXES);
  AL_Is_boundary (SZ_Float_Vect, is_gbeg, is_gend);

  for (idim = 0; idim < DIMENSIONS; idim++) {
    grid->nghost[idim]      = nghost;
    grid->np_tot[idim]      = lsize[idim] + 2*ghosts[idim];
    grid->np_int[idim]      = lsize[idim];
    grid->np_tot_glob[idim] = runtime->npoint[idim] + 2*ghosts[idim];
    grid->np_int_glob[idim] = runtime->npoint[idim];
    grid->beg[idim]         = beg[idim];
    grid->end[idim]         = end[idim];
    grid->gbeg[idim]        = gbeg[idim];
    grid->gend[idim]        = gend[idim];
    grid->lbeg[idim]        = lbeg[idim];
    grid->lend[idim]        = lend[idim];
    grid->lbound[idim] = runtime->left_bound[idim]*is_gbeg[idim];
    grid->rbound[idim] = runtime->right_bound[idim]*is_gend[idim];
    grid->nproc[idim]  = procs[idim];
  }

/* -- Find total number of processors & decomposition mode -- */

/*    AL_Get_size(SZ, &nprocs);  Replaced by MPI_Comm_size above */

  #ifdef STAGGERED_MHD

/* ---- x-staggered array descriptor (double) ---- */

    #if DIMENSIONS >= 1
    D_EXPAND(stagdim[IDIR] = AL_TRUE;  ,
             stagdim[JDIR] = AL_FALSE; , 
             stagdim[KDIR] = AL_FALSE;)

    DIM_LOOP(idim) gsize[idim] = runtime->npoint[idim];
    gsize[IDIR] += 1;

    #ifdef SHEARINGBOX
    periods[IDIR] = 0;
    #endif

    AL_Sz_init (MPI_COMM_WORLD, &SZ_stagx);
    AL_Set_type (MPI_DOUBLE, 1, SZ_stagx);
    AL_Set_dimensions (DIMENSIONS, SZ_stagx);
    AL_Set_global_dim (gsize, SZ_stagx);
    AL_Set_ghosts (ghosts, SZ_stagx);
    AL_Set_staggered_dim(stagdim, SZ_stagx);
    AL_Set_periodic_dim (periods, SZ_stagx);
    AL_Set_parallel_dim (pardim, SZ_stagx);

    AL_Decompose (SZ_stagx, procs, decomp_mode);
    AL_Get_local_dim (SZ_stagx, lsize);
    AL_Get_bounds (SZ_stagx, beg, end, ghosts, AL_C_INDEXES);
    AL_Get_lbounds (SZ_stagx, lbeg, lend, ghosts, AL_C_INDEXES);
    AL_Get_gbounds (SZ_stagx, gbeg, gend, ghosts, AL_C_INDEXES);
    AL_Is_boundary (SZ_stagx, is_gbeg, is_gend);

    #ifdef SHEARINGBOX
     periods[IDIR] = 1;
    #endif
    #endif

    #if DIMENSIONS >= 2
    D_EXPAND(stagdim[IDIR] = AL_FALSE;  ,
             stagdim[JDIR] = AL_TRUE;   , 
             stagdim[KDIR] = AL_FALSE;)

    DIM_LOOP(idim) gsize[idim] = runtime->npoint[idim];
    gsize[JDIR] += 1;

    AL_Sz_init (MPI_COMM_WORLD, &SZ_stagy);
    AL_Set_type (MPI_DOUBLE, 1, SZ_stagy);
    AL_Set_dimensions (DIMENSIONS, SZ_stagy);
    AL_Set_global_dim (gsize, SZ_stagy);
    AL_Set_ghosts (ghosts, SZ_stagy);
    AL_Set_staggered_dim(stagdim, SZ_stagy);
    AL_Set_periodic_dim (periods, SZ_stagy);
    AL_Set_parallel_dim (pardim, SZ_stagy);

    AL_Decompose (SZ_stagy, procs, decomp_mode);
    AL_Get_local_dim (SZ_stagy, lsize);
    AL_Get_bounds  (SZ_stagy, beg, end, ghosts, AL_C_INDEXES);
    AL_Get_lbounds (SZ_stagy, lbeg, lend, ghosts, AL_C_INDEXES);
    AL_Get_gbounds (SZ_stagy, gbeg, gend, ghosts, AL_C_INDEXES);
    AL_Is_boundary (SZ_stagy, is_gbeg, is_gend);
    #endif

    #if DIMENSIONS == 3
    D_EXPAND(stagdim[IDIR] = AL_FALSE;  ,
             stagdim[JDIR] = AL_FALSE;  , 
             stagdim[KDIR] = AL_TRUE;)

    DIM_LOOP(idim) gsize[idim] = runtime->npoint[idim];
    gsize[KDIR] += 1;

    AL_Sz_init (MPI_COMM_WORLD, &SZ_stagz);
    AL_Set_type (MPI_DOUBLE, 1, SZ_stagz);
    AL_Set_dimensions (DIMENSIONS, SZ_stagz);
    AL_Set_global_dim (gsize, SZ_stagz);
    AL_Set_ghosts (ghosts, SZ_stagz);
    AL_Set_staggered_dim(stagdim, SZ_stagz);
    AL_Set_periodic_dim (periods, SZ_stagz);
    AL_Set_parallel_dim (pardim, SZ_stagz);

    AL_Decompose (SZ_stagz, procs, decomp_mode);
    AL_Get_local_dim (SZ_stagz, lsize);
    AL_Get_bounds  (SZ_stagz, beg, end, ghosts, AL_C_INDEXES);
    AL_Get_lbounds (SZ_stagz, lbeg, lend, ghosts, AL_C_INDEXES);
    AL_Get_gbounds (SZ_stagz, gbeg, gend, ghosts, AL_C_INDEXES);
    AL_Is_boundary (SZ_stagz, is_gbeg, is_gend);
    #endif
  #endif /* STAGGERED_MHD */

  /* -- find processors coordinates in a Cartesian topology -- */

  {
    int coords[3] = {0,0,0};
    int rank;
    MPI_Comm cartcomm;
    AL_Get_cart_comm(SZ, &cartcomm);
    MPI_Cart_get(cartcomm, DIMENSIONS, procs, periods, coords);
    MPI_Cart_rank(cartcomm, coords, &rank);
    if (rank != prank) {
      printf ("! Initialize: rank and prank are different\n");
      QUIT_PLUTO(1);
    }
    for (idim = 0; idim < DIMENSIONS; idim++) {
      grid->rank_coord[idim] = coords[idim];
    }
  }

#else  /* if NOT PARALLEL */

/* -----------------------------------------------------
               Serial Initialization
   ----------------------------------------------------- */

  RuntimeSetup (runtime, cmd_line, ini_file);
  nghost = GetNghost();

  for (idim = 0; idim < DIMENSIONS; idim++) {
    grid->nghost[idim]  = nghost;
    grid->np_int[idim]  = grid->np_int_glob[idim] = runtime->npoint[idim];
    grid->np_tot[idim]  = grid->np_tot_glob[idim] = runtime->npoint[idim] + 2*nghost;
    grid->beg[idim]     = grid->gbeg[idim] = grid->lbeg[idim] = nghost;
    grid->end[idim]     = grid->gend[idim] = grid->lend[idim] 
                        = (grid->lbeg[idim] - 1) + grid->np_int[idim];
    grid->lbound[idim]  = runtime->left_bound[idim];
    grid->rbound[idim]  = runtime->right_bound[idim];
    grid->nproc[idim]   = 1;
  }
  nprocs = 1;

#endif

  RuntimeSet (runtime);

/* ----------------------------------------------------
    Set output directory and log file
   ---------------------------------------------------- */

  SetLogFile  (runtime->log_dir, cmd_line);
  ShowConfig  (argc, argv, ini_file);

/* ---------------------------------------------------
                Grid Generation
   --------------------------------------------------- */

  print ("\n> Generating grid...\n\n");
  SetGrid (runtime, grid);
  Where (-1, grid);     /* -- store grid inside the "Where" 
                              function for subsequent calls -- */
  if (cmd_line->makegrid == YES) {
    print ("\n> Done < \n");
    QUIT_PLUTO(0);
  }

/* ------------------------------------------------
     Initialize global variables
   ------------------------------------------------ */

  g_dt             = runtime->first_dt;
  g_time           = 0.0;
  g_maxMach        = 0.0;
  g_maxRiemannIter = 0;
  g_nprocs         = nprocs;
  g_usedMemory     = 0;
  
  IBEG = grid->lbeg[IDIR]; IEND = grid->lend[IDIR];
  JBEG = grid->lbeg[JDIR]; JEND = grid->lend[JDIR];
  KBEG = grid->lbeg[KDIR]; KEND = grid->lend[KDIR];

  NX1 = grid->np_int[IDIR];
  NX2 = grid->np_int[JDIR];
  NX3 = grid->np_int[KDIR];

  NX1_TOT = grid->np_tot[IDIR]; 
  NX2_TOT = grid->np_tot[JDIR];
  NX3_TOT = grid->np_tot[KDIR];

/* ---------------------------------------
     Get the maximum number of points 
     among all directions
   --------------------------------------- */

  NMAX_POINT = MAX(NX1_TOT, NX2_TOT);
  NMAX_POINT = MAX(NMAX_POINT, NX3_TOT);

/* --------------------------------------------------------------------
     Find the minimum physical cell length for each direction
   -------------------------------------------------------------------- */

  for (idim = 0; idim < DIMENSIONS; idim++)  dxmin[idim] = 1.e30;

  for (i = IBEG; i <= IEND; i++) {
  for (j = JBEG; j <= JEND; j++) {
  for (k = KBEG; k <= KEND; k++) {

    scrh = Length_1(i, j, k, grid);
    dxmin[IDIR] = MIN (dxmin[IDIR], scrh);

    scrh = Length_2(i, j, k, grid);
    dxmin[JDIR] = MIN (dxmin[JDIR], scrh);
 
    scrh = Length_3(i, j, k, grid); 
    dxmin[KDIR] = MIN (dxmin[KDIR], scrh);
     
  }}}

#ifdef PARALLEL
  MPI_Allreduce (dxmin, dxming, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  dxmin[IDIR] = dxming[IDIR];
  dxmin[JDIR] = dxming[JDIR];
  dxmin[KDIR] = dxming[KDIR];
#endif

  grid->dl_min[IDIR] = dxmin[IDIR];
  grid->dl_min[JDIR] = dxmin[JDIR];
  grid->dl_min[KDIR] = dxmin[KDIR];

/* --------------------------------------------------------------
    Define geometrical coefficients for linear interpolation
   -------------------------------------------------------------- */
   
  PLM_CoefficientsSet (grid);   /* -- these may be needed by
                                      shock flattening algorithms */
#if RECONSTRUCTION == PARABOLIC
  PPM_CoefficientsSet (grid);  
#endif

/* --------------------------------------------------------------
    Copy user defined parameters into global array g_inputParam 
   -------------------------------------------------------------- */

  for (nv = 0; nv < USER_DEF_PARAMETERS; nv++) g_inputParam[nv] = runtime->aux[nv];

/* ------------------------------------------------------------
          Allocate memory for 3D data arrays
   ------------------------------------------------------------ */

  print ("\n> Memory allocation\n");
  data->Vc = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
  data->Uc = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double); 

#ifdef STAGGERED_MHD
  data->Vs = ARRAY_1D(DIMENSIONS, double ***);
  D_EXPAND(
    data->Vs[BX1s] = ArrayBox( 0, NX3_TOT-1, 0, NX2_TOT-1,-1, NX1_TOT-1); ,
    data->Vs[BX2s] = ArrayBox( 0, NX3_TOT-1,-1, NX2_TOT-1, 0, NX1_TOT-1); ,
    data->Vs[BX3s] = ArrayBox(-1, NX3_TOT-1, 0, NX2_TOT-1, 0, NX1_TOT-1);
  )

  data->Ex1 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  data->Ex2 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  data->Ex3 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

  data->emf = malloc(sizeof(EMF));
  CT_Allocate (data->emf);
#endif  

/* ------------------------------------------------------------
    Allocate memory for vector potential.
   ------------------------------------------------------------ */

#if UPDATE_VECTOR_POTENTIAL == YES || ASSIGN_VECTOR_POTENTIAL == YES
  data->Ax3 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);   
  data->Ax1 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  
  data->Ax2 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  
#endif

#if RESISTIVITY || HALL_MHD
  data->J = ARRAY_4D(3,NX3_TOT, NX2_TOT, NX1_TOT, double);
#endif

#if THERMAL_CONDUCTION
  data->Tc = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
#endif

#if FORCED_TURB == YES
  data->Ft = malloc(sizeof(ForcedTurb));
#endif

#ifdef PARTICLES
  data->PHead = NULL;
  data->pstr  = NULL;
  #if PARTICLES_TYPE == COSMIC_RAYS
  data->Ecr = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
  data->Fcr = ARRAY_4D(4, NX3_TOT, NX2_TOT, NX1_TOT, double);
  data->Jcr = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
  data->qcr = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  #endif
  #if PARTICLES_TYPE == DUST
  data->Vdust = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
  data->Fdust = ARRAY_4D(3, NX3_TOT, NX2_TOT, NX1_TOT, double);
  #endif
#endif

  data->flag = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, unsigned char);

/* ------------------------------------------------------------
    Initialize tables needed for EOS 
   ------------------------------------------------------------ */
  
#if EOS == PVTE_LAW && NIONS == 0
  #if TV_ENERGY_TABLE == YES
  MakeInternalEnergyTable();
  #else
  MakeEV_TemperatureTable();
  #endif
  MakePV_TemperatureTable();
#endif

/* ------------------------------------------------------------
    Assign initial conditions.
    We also compute the average orbital velocity in case
    output requires writing the residual.
   ------------------------------------------------------------ */

  Startup (data, grid);
#ifdef PARTICLES
  if (cmd_line->restart == NO) Particles_Set(data, grid); 
#endif
#ifdef FARGO
  FARGO_ComputeVelocity(data, grid);
#endif

/* --------------------------------------------------------
    Set output attributes (names, images, number of
    outputs...)
   -------------------------------------------------------- */

  SetOutput (data, runtime);
#ifdef PARTICLES
  Particles_SetOutput(data, runtime);
#endif
  ChangeOutputVar();

/* --------------------------------------------------------
      print normalization units
   -------------------------------------------------------- */

  ShowUnits();  

  print ("> Number of processors: %d\n",nprocs);
  D_EXPAND(
           print ("> Proc size:            %d",grid->np_int[IDIR]);  ,
           print (" X %d", grid->np_int[JDIR]);                      ,
           print (" X %d", grid->np_int[KDIR]);)
  print ("\n");
#ifdef PARALLEL
  print ("> Parallel Directions: ");
  D_EXPAND(if (pardim[IDIR]) print (" X1");  ,
           if (pardim[JDIR]) print ("/X2");  ,
           if (pardim[KDIR]) print ("/X3");)
  print ("\n");
#endif

/* --------------------------------------------------------
   Enable writing of grid information to be read by
   gnuplot. This is done for educational purposes only.
   -------------------------------------------------------- */

#if !(defined CHOMBO) && (defined GNUPLOT)
  GnuplotSetting(grid);
#endif

}

#ifdef PARALLEL
/* ********************************************************************* */
int GetDecompMode (Cmd_Line *cmd_line, int procs[])
/*!
 * Returns the parallel domain decomposition mode.
 *
 * \param [in]  cmd_line  pointer to the Cmd_Line structure
 * \param [out] procs     an array of integers giving the number
 *                        of processors in each direction only if
 *                        the -dec command line option has been given
 *
 *  \return  The decomposition mode:
 *  
 *   - AL_AUTO_DECOMP   defaults
 *   - AL_USER_DECOMP   if the -dec n1 [n2] [n3] command line 
 *                      argument has been given;
 *                      In this case only, procs[] contains the number
 *                      of processors in the three directions;
 *   -  AL_MPI_DECOMP   if [todo]
 * \todo AL_MPI_DECOMP mode
 *********************************************************************** */ 
{
  int  nprocs;
  long int npx = 1, npy = 1, npz = 1;

  D_EXPAND(npx = cmd_line->nproc[IDIR];  ,
           npy = cmd_line->nproc[JDIR];  ,
           npz = cmd_line->nproc[KDIR];)

/* ------------------------------------------------
    if -dec has not been given
    decomposition will be set to be AL_AUTO_DECOMP
   ------------------------------------------------ */

  if (npx == -1 || npy == -1 || npz == -1){
    return AL_AUTO_DECOMP;
  }

/* -----------------------------------------------
    enter in user decomp mode if the number of
    processors along all directions has been given
   ------------------------------------------------ */

  if (npx > 0 && npy > 0 && npz > 0){
    procs[IDIR] = npx;
    procs[JDIR] = npy;
    procs[KDIR] = npz;
    
  /* -- check if decomposition is correct -- */
  
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (procs[IDIR]*procs[JDIR]*procs[KDIR] != nprocs){
      printf ("! The specified parallel decomposition (%d, %d, %d) is not\n",
               procs[IDIR],procs[JDIR],procs[KDIR]);
      printf ("! consistent with the number of processors (%d).\n",nprocs);
      QUIT_PLUTO(1);
    }
    return AL_USER_DECOMP;
  }

  printf ("! GetDecompMode: invalid decomposition mode");
  QUIT_PLUTO(1);
}
#endif
