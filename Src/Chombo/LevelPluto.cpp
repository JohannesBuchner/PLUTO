#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>

#include "LevelPluto.H"
#include "BoxIterator.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "LoHiSide.H"

#include "CH_Timer.H"

#include "NamespaceHeader.H"

// Constructor - set up some defaults
LevelPluto::LevelPluto()
{
  m_dx           = 0.0;
  m_dl_min       = 1.e30;
  m_refineCoarse = 0;
  m_patchPluto   = NULL;
  m_isDefined    = false;
}

// Destructor - free up storage
LevelPluto::~LevelPluto()
{
  if (m_patchPluto != NULL)
    {
      delete m_patchPluto;
    }
}

// Define the object so that time stepping can begin
void LevelPluto::define(const DisjointBoxLayout&  a_thisDisjointBoxLayout,
                        const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
                        const ProblemDomain&      a_domain,
                        const int&                a_refineCoarse,
                        const int&                a_level,
                        const Real&               a_dx,
                        const PatchPluto*         a_patchPlutoFactory,
                        const bool&               a_hasCoarser,
                        const bool&               a_hasFiner)
{
  CH_TIME("LevelPluto::define");

  // Sanity checks
  CH_assert(a_refineCoarse > 0);
  CH_assert(a_dx > 0.0);

  // Make a copy of the current grids
  m_grids  = a_thisDisjointBoxLayout;

  // Cache data
  m_dx = a_dx;
  m_level = a_level;
  m_domain = a_domain;
  m_refineCoarse = a_refineCoarse;
  m_hasCoarser = a_hasCoarser;
  m_hasFiner = a_hasFiner;

  // Remove old patch integrator (if any), create a new one, and initialize
  if (m_patchPluto != NULL)
    {
      delete m_patchPluto;
    }

 // Determing the number of ghost cells necessary here
  m_numGhost = GetNghost();

  m_patchPluto = a_patchPlutoFactory->new_patchPluto();
  m_patchPluto->define(m_domain,m_dx,m_level,m_numGhost);
 
  // Set the grid for the entire level
  setGridLevel();

  // Get the number of conserved variable and face centered fluxes
  m_numCons   = m_patchPluto->numConserved();
  m_numFluxes = m_patchPluto->numFluxes();

  m_exchangeCopier.exchangeDefine(a_thisDisjointBoxLayout,
                                  m_numGhost*IntVect::Unit);

  // Setup an interval corresponding to the conserved variables
  Interval UInterval(0,m_numCons-1);

 #if (TIME_STEPPING == RK2)
  // Create temporary storage with a layer of "m_numGhost" ghost cells
  // for the flags passing from predictor to corrector (RK2 only)
  m_Flags.define(m_grids,1,m_numGhost*IntVect::Unit);
  m_Utmp.define(m_grids,m_numCons,m_numGhost*IntVect::Unit);
 #endif

  // Create temporary storage with a layer of "m_numGhost" ghost cells
  {
    CH_TIME("setup::Udefine");
    m_U.define(m_grids,m_numCons,m_numGhost*IntVect::Unit);
  }

  // Set up the interpolator if there is a coarser level
  if (m_hasCoarser)
    {
      m_patcher.define(a_thisDisjointBoxLayout,
                       a_coarserDisjointBoxLayout,
                       m_numCons,
                       coarsen(a_domain,a_refineCoarse),
                       a_refineCoarse,
                       m_dx,
                       m_numGhost);
    }

  // Everything is defined
  m_isDefined = true;
}

// Advance the solution by "a_dt" by using an unsplit method.
// "a_finerFluxRegister" is the flux register with the next finer level.
// "a_coarseFluxRegister" is flux register with the next coarser level.
// If source terms do not exist, "a_S" should be null constructed and not
// defined (i.e. its define() should not be called).
Real LevelPluto::step(LevelData<FArrayBox>&       a_U,
                      LevelData<FArrayBox>        a_flux[CH_SPACEDIM],
                      LevelFluxRegister&          a_finerFluxRegister,
                      LevelFluxRegister&          a_coarserFluxRegister,
                      LevelData<FArrayBox>&       a_split_tags,
                      const LevelData<FArrayBox>& a_UCoarseOld,
                      const Real&                 a_TCoarseOld,
                      const LevelData<FArrayBox>& a_UCoarseNew,
                      const Real&                 a_TCoarseNew,
                      const Real&                 a_time,
                      const Real&                 a_dt,
                      const Real&                 a_cfl)
{
  CH_TIMERS("LevelPluto::step");

  CH_TIMER("LevelPluto::step::setup"   ,timeSetup);
  CH_TIMER("LevelPluto::step::update"  ,timeUpdate);
  CH_TIMER("LevelPluto::step::reflux"  ,timeReflux);
  CH_TIMER("LevelPluto::step::conclude",timeConclude);

  // Make sure everything is defined
  CH_assert(m_isDefined);

  CH_START(timeSetup);

  // Clear flux registers with next finer level
  if (m_hasFiner && (g_intStage == 1))
    {
      a_finerFluxRegister.setToZero();
    }

  // Setup an interval corresponding to the conserved variables
  Interval UInterval(0,m_numCons-1);

  {
    CH_TIME("setup::localU");
    for (DataIterator dit = m_U.dataIterator(); dit.ok(); ++dit)
    {
      m_U[dit].setVal(0.0); // Gets rid of denormalized crap.
      m_U[dit].copy(a_U[dit]);
    }

    m_U.exchange(m_exchangeCopier);
  }

  // Fill m_U's ghost cells using fillInterp
  if (m_hasCoarser)
    {
      // Fraction "a_time" falls between the old and the new coarse times
      Real alpha = (a_time - a_TCoarseOld) / (a_TCoarseNew - a_TCoarseOld);

      if (alpha > 1.0) alpha = 1.0;
      if (alpha < 0.0) alpha = 0.0;

      // Interpolate ghost cells from next coarser level using both space
      // and time interpolation
      m_patcher.fillInterp(m_U,
                           a_UCoarseOld,
                           a_UCoarseNew,
                           alpha,
                           0,0,m_numCons);
    }

  // Potentially used in boundary conditions
  m_patchPluto->setCurrentTime(a_time);

  // Use to restrict maximum wave speed away from zero
  Real maxWaveSpeed = 1.e-12;
  Real minDtCool    = 1.e38;

  // The grid structure
  Grid *grid;
  static timeStep Dts;
  Real inv_dt;
  
  #ifdef GLM_MHD
   glm_ch = g_coeff_dl_min*m_dx/(a_dt + 1.e-16)*a_cfl;
//   glm_ch = g_coeff_dl_min/(a_dt + 1.e-16)*a_cfl; /* If subcycling is turned off */
   glm_ch = MIN(glm_ch,glm_ch_max*g_coeff_dl_min);
  #endif

  CH_STOP(timeSetup);
  g_level_dx = m_dx;

  // Beginning of loop through patches/grids.
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    CH_START(timeUpdate);

    // The current box
    Box curBox = m_grids.get(dit());

    // The current grid of conserved variables
    FArrayBox& curU = m_U[dit];

    // The current grid of volumes
    #if GEOMETRY != CARTESIAN
     const FArrayBox& curdV = m_dV[dit()];
    #else
     const FArrayBox  curdV;
    #endif
 
   #ifdef SKIP_SPLIT_CELLS
    // The current grid of split/unsplit tags
    FArrayBox& split_tags = a_split_tags[dit];
   #else
    FArrayBox split_tags;
   #endif

   #if (TIME_STEPPING == RK2)
    // The current storage for flags (RK2 only)
    BaseFab<unsigned char>& flags = m_Flags[dit];
    // Local temporary storage for conserved variables
    FArrayBox& curUtmp = m_Utmp[dit];
   #else
    BaseFab<unsigned char> flags;
    FArrayBox curUtmp;
   #endif

    // The fluxes computed for this grid - used for refluxing and returning
    // other face centered quantities
    FluxBox flux;

    // Set the current box for the patch integrator
    m_patchPluto->setCurrentBox(curBox);

    Real minDtCoolGrid;
    
    grid = m_structs_grid[dit].getGrid();
 
    IBEG = grid->lbeg[IDIR]; IEND = grid->lend[IDIR];
    JBEG = grid->lbeg[JDIR]; JEND = grid->lend[JDIR];
    KBEG = grid->lbeg[KDIR]; KEND = grid->lend[KDIR];

    NX1 = grid->np_int[IDIR];
    NX2 = grid->np_int[JDIR];
    NX3 = grid->np_int[KDIR];

    NX1_TOT = grid->np_tot[IDIR];
    NX2_TOT = grid->np_tot[JDIR];
    NX3_TOT = grid->np_tot[KDIR];

    g_dt   = a_dt;
    g_time = a_time;
    g_maxRiemannIter = 0;
    PLM_CoefficientsSet (grid);  /* -- these may be needed by
                                       shock flattening algorithms */
    #if RECONSTRUCTION == PARABOLIC
     PPM_CoefficientsSet (grid);  
    #endif
    
    // reset time step coefficients 
    if (Dts.cmax == NULL) Dts.cmax = ARRAY_1D(NMAX_POINT, double);
    int id;
    Dts.invDt_hyp = 1.e-18;
    Dts.invDt_par = 1.e-18;
    Dts.dt_cool = 1.e18;
    Dts.cfl     = a_cfl;
    Where(-1, grid); /* -- store grid for subsequent calls -- */

    // Take one step
    m_patchPluto->advanceStep (curU, curUtmp, curdV, split_tags, flags, flux,
                               &Dts, curBox, grid);
 
    inv_dt = Dts.invDt_hyp + 2.0*Dts.invDt_par;
    maxWaveSpeed = Max(maxWaveSpeed, inv_dt); // Now the inverse of the timestep

    minDtCool = Min(minDtCool, Dts.dt_cool/a_cfl);

    CH_STOP(timeUpdate);

    CH_START(timeReflux);

    // Do flux register updates
    for (int idir = 0; idir < SpaceDim; idir++) {
    // Increment coarse flux register between this level and the next
    // finer level - this level is the next coarser level with respect
    // to the next finer level
      if (m_hasFiner) {
        a_finerFluxRegister.incrementCoarse(flux[idir],a_dt,dit(),
                                            UInterval, UInterval,idir);
      }

      // Increment fine flux registers between this level and the next
      // coarser level - this level is the next finer level with respect
      // to the next coarser level
       if (m_hasCoarser) {
         a_coarserFluxRegister.incrementFine(flux[idir],a_dt,dit(),
                                             UInterval, UInterval,idir);
       }
    }

    CH_STOP(timeReflux);
  }

  CH_START(timeConclude);

  {
    CH_TIME("conclude::copyU");
    // Now that we have completed the updates of all the patches, we copy the
    // contents of temporary storage, U, into the permanent storage, a_U.
    for(DataIterator dit = m_U.dataIterator(); dit.ok(); ++dit){
      a_U[dit].copy(m_U[dit]);
    }
  }

 // Find the minimum of dt's over this level
  Real dtNew;
 {
  CH_TIME("conclude::getDt");

  Real local_dtNew = 1. / maxWaveSpeed;
  #if (TIME_STEPPING == RK2) && (COOLING != NO)
   if (g_intStage == 2) local_dtNew = minDtCool;
  #else
   local_dtNew = Min(local_dtNew,minDtCool);
  #endif
  dtNew = local_dtNew;

  #ifdef CH_MPI
   #if (TIME_STEPPING == RK2) && (COOLING == NO)
    if (g_intStage == 1) {
   #endif
     int result = MPI_Allreduce(&local_dtNew, &dtNew, 1, MPI_CH_REAL,
                                    MPI_MIN, Chombo_MPI::comm);
     if(result != MPI_SUCCESS){ //bark!!!
        MayDay::Error("sorry, but I had a communcation error on new dt");
     }
   #if (TIME_STEPPING == RK2) && (COOLING == NO)
    }
   #endif
  #endif
 }

  CH_STOP(timeConclude);

  // Return the maximum stable time step
  return dtNew;
}

void LevelPluto::setGridLevel()
{

 CH_TIME("LevelPluto::setGrid");

 m_structs_grid.define(m_grids);

 #if GEOMETRY != CARTESIAN
  m_dV.define(m_grids,CHOMBO_NDV,m_numGhost*IntVect::Unit);
 #endif

 Real dlMinLoc = 1.e30;

 for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
   {
      // The current box
      Box curBox = m_grids.get(dit());
      Grid *grid = m_structs_grid[dit].getGrid();

      #if GEOMETRY != CARTESIAN 
       FArrayBox& curdV = m_dV[dit()];    
      #else
       FArrayBox  curdV; 
      #endif
       
      m_patchPluto->setGrid(curBox, grid, curdV);           
       
      for (int idir = 0; idir < SpaceDim; idir++) {
        dlMinLoc = Min(dlMinLoc,grid->dl_min[idir]);
      }
   }

#if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)

   D_EXPAND(m_dl_min = m_dx; ,
            m_dl_min = MIN(m_dl_min,m_dx*g_x2stretch); ,
            m_dl_min = MIN(m_dl_min,m_dx*g_x3stretch); )
#else

 #ifdef CH_MPI
  Real dlMin;
  int result = MPI_Allreduce(&dlMinLoc, &dlMin, 1, MPI_CH_REAL,
                             MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){ //bark!!!
   MayDay::Error("sorry, but I had a communcation error on dlMin");
  }
  m_dl_min = dlMin;
 #else
  m_dl_min = dlMinLoc;
 #endif

#endif

}

#if GEOMETRY != CARTESIAN
const LevelData<FArrayBox>& LevelPluto::getdV() const
{
  return m_dV;
}
#endif

Real LevelPluto::getDlMin()
{

 CH_TIME("LevelPluto::getDlMin");

 return m_dl_min / m_dx;
// return m_dl_min; /* If subcycling is turned off */


}

#include "NamespaceFooter.H"
