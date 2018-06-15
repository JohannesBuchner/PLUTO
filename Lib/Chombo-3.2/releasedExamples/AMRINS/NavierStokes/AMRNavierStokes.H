#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRNAVIERSTOKES_H_
#define _AMRNAVIERSTOKES_H_

#include "AMRLevel.H"
#include "LevelFluxRegister.H"
#include "CoarseAverage.H"
#include "PiecewiseLinearFillPatch.H"
#include "PatchGodunov.H"
#include "PhysIBC.H"
#include "AMRMultiGrid.H"
#include "AMRPoissonOp.H"
#include "RelaxSolver.H"
#include "LevelTGA.H"

#include "PhysBCUtil.H"
#include "CCProjector.H"
#include "VelocityAMRPoissonOp.H"

#include "UsingNamespace.H"

///
enum viscousSolverTypes
{
  ///
  backwardEuler = 0,
  ///
  CrankNicolson,
  ///
  TGA,
  ///
  NUM_SOLVER_TYPES
};

/// scalar types -- defined to faciliate lagrangian tracking
enum scalTypes
{
  ///
  xMarker = 0,
  ///
  yMarker,
#if CH_SPACEDIM >= 3
  ///
  zMarker,
#endif
  ///
  NUM_SCAL_TYPES
};

#define ADVECT_GROW 4

#ifdef CH_USE_DOUBLE
#define TIME_EPS 1.0e-10
#else
#define TIME_EPS 1.0e-5
#endif

// small parameter used in defining TGA constants for viscous solves
// #define TGA_EPS 1.0e-8
// petermc, 3 dec 2007:  same as in BaseLevelTGA.H
#define TGA_EPS 1.0e-12

///
class AMRNavierStokes: public AMRLevel
{
public:

  ///
  AMRNavierStokes();

  ///
  AMRNavierStokes(AMRLevel* a_coarser_level_ptr,
           const Box& a_prob_domain,
           int a_level, int a_ref_ratio);

  ///
  AMRNavierStokes(AMRLevel* a_coarser_level_ptr,
           const ProblemDomain& a_prob_domain,
           int a_level, int a_ref_ratio);

  ///
  virtual ~AMRNavierStokes();

  ///
  virtual AMRLevel* makeAMRLevel(AMRLevel* a_coarser_level_ptr,
                                 const Box& a_problem_domain,
                                 int a_level, int a_ref_ratio) const;

  ///
  virtual AMRLevel* makeAMRLevel(AMRLevel* a_coarser_level_ptr,
                                 const ProblemDomain& a_problem_domain,
                                 int a_level, int a_ref_ratio) const;

  ///
  virtual void define(AMRLevel* a_coarse_level_ptr,
                      const Box& a_problem_domain,
                      int a_level, int a_ref_ratio);

  ///
  virtual void define(AMRLevel* a_coarse_level_ptr,
                      const ProblemDomain& a_problem_domain,
                      int a_level, int a_ref_ratio);

  /// advance by one timestep (returns maximum allowable timestep)
  virtual Real advance();

  /// things to do after a timestep
  virtual void postTimeStep();

  /// create tags
  virtual void tagCells(IntVectSet& a_tags);

  /// create tags at initialization
  virtual void tagCellsInit(IntVectSet& a_tags);

  /// regrid
  virtual void regrid(const Vector<Box>& a_new_grids);

  /// perform any post-regridding ops -- lBase is the finest unchanged level
  virtual void postRegrid(int a_lBase);

  /// initialize grids
  virtual void initialGrid(const Vector<Box>& a_new_grids);

  /// initialize data
  virtual void initialData();

  /// things do do after initialization
  virtual void postInitialize();

#ifdef CH_USE_HDF5
  /// write checkpoint header
  virtual void writeCheckpointHeader(HDF5Handle& a_handle) const;

  /// write checkpoint data for this level
  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;

  /// read checkpoint header
  virtual void readCheckpointHeader(HDF5Handle& a_handle);

  /// read checkpoint data for this level
  virtual void readCheckpointLevel(HDF5Handle& a_handle);

  /// write plotfile header
  virtual void writePlotHeader(HDF5Handle& a_handle) const;

  /// write plotfile data for this level
  virtual void writePlotLevel(HDF5Handle& a_handle) const;
#endif

  /// compute dt
  virtual Real computeDt();

  /// compute dt with initial data
  virtual Real computeInitialDt();

  /// set CFL number
  virtual void CFL(Real a_cfl);

  /// set refinement threshold
  virtual void refinementThreshold(Real a_refine_threshold);

  /// set solver parameter
  void limitSolverCoarsening(bool a_limitSolverCoarsening);

  /// is this this finest level?
  bool finestLevel() const;

  /// is this an empty level?
  bool isEmpty() const;

  // compute vorticity
  void computeVorticity(LevelData<FArrayBox>& a_vorticity) const;

  // compute kinetic energy (only on interior cells)
  void computeKineticEnergy(LevelData<FArrayBox>& a_energy) const;

  /// access functions

  ///
  Real Dx() const;

  /// viscous coefficient
  Real Nu() const;

  /// encapsulates the finer-level-pointer casting
  AMRNavierStokes* finerNSPtr() const;

  /// encapsulates the coarser-level pointer casting
  AMRNavierStokes* crseNSPtr() const;

  ///
  LevelData<FArrayBox>& newVel();

  ///
  const LevelData<FArrayBox>& newVel() const;

  ///
  LevelData<FArrayBox>& oldVel();

  ///
  const LevelData<FArrayBox>& oldVel() const;

  /// returns velocity at time t
  void velocity(LevelData<FArrayBox>& a_vel, Real a_time) const;

  ///
  LevelData<FArrayBox>& newLambda();

  ///
  const LevelData<FArrayBox>& newLambda() const;

  ///
  LevelData<FArrayBox>& oldLambda();

  ///
  const LevelData<FArrayBox>& oldLambda() const;

  /// returns lambda at time t
  void lambda(LevelData<FArrayBox>& a_lambda, Real a_time) const;

  ///
  LevelData<FArrayBox>& newScal(const int a_comp);

  ///
  const LevelData<FArrayBox>& newScal(const int a_comp) const;

  ///
  LevelData<FArrayBox>& oldScal(const int a_comp);

  ///
  const LevelData<FArrayBox>& oldScal(const int a_comp) const;

  /// returns scalars at time t
  void scalars(LevelData<FArrayBox>& a_scalars, const Real a_time,
               const int a_comp) const;

  /// sets physical BC object
  void setPhysBC(const PhysBCUtil& a_BC);

  /// gets physical BC object
  PhysBCUtil* getPhysBCPtr() const;

protected:
  ///
  DisjointBoxLayout loadBalance(const Vector<Box>& a_grids);

  ///
  void levelSetup (const DisjointBoxLayout& level_domain);

  /// read parameters from parmParse database
  void readParameters();

  /// set whether this is the finest level or not
  void finestLevel(bool a_finest_level);

  /// move new-time state into old-time state -- also advances time by dt
  void swapOldAndNewStates();

  /// exchange old and new states, resets time to a_time
  void resetStates(const Real a_time);

  /// compute advection velocities
  void computeAdvectionVelocities(LevelData<FluxBox>& a_adv_vel);

  /// do predictor
  /**  returns approx. to (U dot del)U.
        Also updates flux registers with momentum fluxes and
        part of viscous fluxes which contains the old-time velocity.
  */
  void predictVelocities(LevelData<FArrayBox>& a_uDelU,
                         LevelData<FluxBox>& a_advVel);

  /// computes approximation to new-time velocity before projection
  /** This function computes uStar -- the approximation
      of the new-time velocity which will be projected.  The
      advection term (u-dot-del-U) is passed in in a_uStar, which
      is modified in place to contain uStar. Note that both viscous
      and inviscid cases, this returns the approximation to
      u^{n+1} + dt*gradPi */
  void computeUStar(LevelData<FArrayBox>& a_uStar);

  /// do scalar advection -- assumes all BC's already done
  void advectScalar(LevelData<FArrayBox>& a_new_scal,
                    LevelData<FArrayBox>& a_old_scal,
                    LevelData<FluxBox>&   a_adv_vel,
                    LevelFluxRegister*    a_crseFluxRegPtr,
                    LevelFluxRegister&    a_fineFluxReg,
                    PatchGodunov&         a_patchGodScalar,
                    Real                  a_dt);

  /// advect/diffuse scalar -- assumes all BC's already done?
  void advectDiffuseScalar(LevelData<FArrayBox>&       a_new_scal,
                           LevelData<FArrayBox>&       a_old_scal,
                           LevelData<FluxBox>&         a_adv_vel,
                           const Real&                 a_diffusive_coeff,
                           const LevelData<FArrayBox>* a_crse_scal_old,
                           const LevelData<FArrayBox>* a_crse_scal_new,
                           const Real                  a_old_crse_time,
                           const Real                  a_new_crse_time,
                           LevelFluxRegister*          a_crseFluxRegPtr,
                           LevelFluxRegister&          a_fineFluxReg,
                           PatchGodunov&               a_patchGodScalar,
                           BCHolder&                   a_scalPhysBC,
                           Real                        a_dt,
                           int                         a_comp);

  /// manages pressure initialization after init or regrid
  void initializeGlobalPressure();

  /// does level-based pressure initalization
  void initializeLevelPressure(Real a_currentTime, Real a_dtInit);

  /// smooths composite velocity field after regridding
  void smoothVelocityField(int lbase);

  /// compute Laplacian(velocity)
  void computeLapVel(LevelData<FArrayBox>& a_lapVel,
                     LevelData<FArrayBox>& a_vel,
                     const LevelData<FArrayBox>* a_crseVelPtr);

  /// compute Laplacian(scalars)
  void computeLapScal(LevelData<FArrayBox>& a_lapScal,
                      LevelData<FArrayBox>& a_scal,
                      const BCHolder& a_physBC,
                      const LevelData<FArrayBox>* a_crseScalPtr);

  /// fill grown velocity field using piecewise-constant interp
  void fillVelocity(LevelData<FArrayBox>& a_vel, Real a_time);

  /// fill grown velocity field using piecewise-constant interp
  void fillVelocity(LevelData<FArrayBox>& a_vel, Real a_time,
                    int vel_comp, int dest_comp, int num_comp);

  /// resize and fill temporary lambda for advection tracing computation
  /**  fills both interior and boundary cells
   */
  void fillLambda(LevelData<FArrayBox>& a_lambda, Real a_time) const;

  ///  fill interior and boundary cells (assume a_scal already sized)
  void fillScalars(LevelData<FArrayBox>& a_scal, Real a_time,
                   const int a_comp) const;

  /// number of components in plotfiles
  virtual int numPlotComps() const;

  /// stuff plotfile data
  virtual void getPlotData(LevelData<FArrayBox>& a_plot_data) const;

  /// defines operator factory used in smoothing during regridding
  void  defineRegridAMROp(AMRPoissonOpFactory& a_factory,
                          const Vector<DisjointBoxLayout>& a_grids,
                          const Vector<ProblemDomain>& a_domains,
                          const Vector<Real>& a_amrDx,
                          const Vector<int>& a_refRatios,
                          const int& a_lBase);

  /// defines multigrid solvers
  void defineViscousMGSolver(const DisjointBoxLayout& a_grids,
                             const DisjointBoxLayout* a_crseGridsPtr,
                             int                      a_refCrse);

protected:

  /// if s_set_bogus_values, assigns s_bogus_value to *m_vel_old_ptr, *m_vel_new_ptr, *m_lambda_old_ptr, *m_lambda_new_ptr, m_scal_old, m_scal_new.
  void setAllBogus();

  /// velocity vector at old time
  LevelData<FArrayBox>* m_vel_old_ptr;

  /// velocity vector at new time
  LevelData<FArrayBox>* m_vel_new_ptr;

#ifdef DEBUG
  /// velocity vector at last multilevel time
  LevelData<FArrayBox>* m_vel_save_ptr;

  Real m_saved_time;

#endif

  /// freestream preservation qty at old time
  LevelData<FArrayBox>* m_lambda_old_ptr;

  /// freestream preservation qty at new time
  LevelData<FArrayBox>* m_lambda_new_ptr;

#ifdef DEBUG
  /// velocity vector at last multilevel time
  LevelData<FArrayBox>* m_lambda_save_ptr;
#endif

  /// scalars at old time
  Vector<LevelData<FArrayBox>*> m_scal_old;

  /// scalars at new time
  Vector<LevelData<FArrayBox>*> m_scal_new;

#ifdef DEBUG
  /// scalars at last multilevel time
  Vector<LevelData<FArrayBox>*> m_scal_save;
#endif

  /// number of components of m_state (vel comps + lambda)
  static const int s_num_vel_comps = SpaceDim;

  /// number of scalars
  static int s_num_scal_comps;

  /// number of components per scalar - hardwired to 1
  static const int s_compsPerScalar = 1;

  /// maximum possible number of scalar components
  static const int s_max_scal_comps=SpaceDim;

  /// names of components
  static const char* s_vel_names[s_num_vel_comps];

  /// names of scalars
  static const char* s_scal_names[s_max_scal_comps];

  /// diffusion coefficients for scalar advection-diffusion
  static Vector<Real> s_scal_coeffs;

  /// size of physical domain
  static Vector<Real> s_domLength;

  /// have we initialized static variables from parmParse?
  static bool s_ppInit;

  /// factor to shrink initial timestep
  static Real s_init_shrink;

  /// maximum allowable timestep
  static Real s_max_dt;

  static Real s_prescribedDt;

  /// should we project initial velocity field (default is yes)
  static bool s_project_initial_vel;

  /// should we initialize pressures after grid generation & regridding?
  static bool s_initialize_pressures;

  /// do antidiffusive smoothing after regridding?
  static bool s_smooth_after_regrid;

  /// antidiffusive smoothing coefficient (multiplies nu*dtCrse)
  static Real s_regrid_smoothing_coeff;

  /// number of times to iterate when initializing pressures
  static int s_num_init_passes;

  /// should we reflux momentum?
  static bool s_reflux_momentum;

  /// if we're refluxing momentum, should we reflux the normal component?
  static bool s_reflux_normal_momentum;

  /// debugging option -- initially set arrays to s_bogus_value
  static bool s_set_bogus_values;

  /// debugging option -- value to which arrays are initially set
  static Real s_bogus_value;

  /// tag on vorticity?
  static bool s_tag_vorticity;

  /// tag on scalar (theta)?
  static bool s_tag_theta;

  /// factor for tagging -- tag if vorticity greater than factor*max_vort
  static Real s_vort_factor;

  /// amount by which to grow tagged regions before calling meshrefine
  static int s_tags_grow;

  /// do implicit refluxing for velocity?
  static bool s_implicit_reflux;

  /// apply freestream preservation correction?
  static bool s_applyFreestreamCorrection;

  /// should we reflux scalars?
  static bool s_reflux_scal;

  /// do implicit refluxing for advected/diffused scalars?
  static bool s_implicit_scal_reflux;

  /// viscous solver type: 0 = backwards euler, 1=Crank-Nicolson, 2=TGA
  static int s_viscous_solver_type;

  /// viscous solver tolerance (defaults to 1e-7)
  static Real s_viscous_solver_tol;

  /// multigrid parameter for viscous solves
  static int s_viscous_num_smooth_up;

  /// multigrid parameter for viscous solves
  static int s_viscous_num_smooth_down;

  /// viscosity
  static Real s_nu;

  /// do we specify initial grids in a gridFile?
  static bool s_specifyInitialGrids;

  /// gridfile where initial grids are spec'd, if relevant
  static string s_initialGridFile;

  /// initialize velocity field from vorticity field
  static bool s_init_vel_from_vorticity;

  /// background velocity (in the z-direction)
  static Real s_backgroundVel;

  /// -- I/O options --

  /// include divergence in plotfile?
  static bool s_write_divergence;

  /// include lambda in plotfile?
  static bool s_write_lambda;

  /// include du/dt and d(lambda)/dt in plotfile?
  static bool s_write_time_derivatives;

  /// include vorticity in plotfile?
  static bool s_write_vorticity;

  /// include scalars in plotfile?
  static bool s_write_scalars;

  /// include d(scalars)/dt in plotfile?
  static bool s_write_dScalar_dt;

  /// include strain-rates in plotfile? (2D only)
  static bool s_write_strains;

  /// include freestream preservation correcion in plotfile?
  static bool s_write_grad_eLambda;

  /// include error (if there is an exact solution) in plotfile?
  static bool s_write_error;

  /// include curl of error (if there is an exact solution) in plotfile?
  //static bool s_write_curl_error;

  /// include processor distribution visualization in plotfile?
  static bool s_write_proc_ids;

  /// compute error and report norms after composite timesteps
  static bool s_compute_scal_err;

  /// write out grids and processor assignments
  static bool s_write_grids;

  /// points where data is dumped out
  Vector<Tuple<Real, SpaceDim> > m_dataPoints;

  /// grid spacing
  Real m_dx;

  /// cfl number
  Real m_cfl;

  /// solver parameter -- passed on through when defining solvers
  bool m_limitSolverCoarsening;

  /// dt from most recent timestep
  Real m_dt_save;

  /// number of steps on this level
  int m_level_steps;

  /// averaging from fine to coarse level for velocity
  CoarseAverage m_coarse_average;

  /// averaging from fine to coarse level for freestream preservation qty
  CoarseAverage m_coarse_average_lambda;

  /// cell-centered projection object
  CCProjector m_projection;

  /// momentum flux register
  LevelFluxRegister m_flux_register;

  /// freestream preservation qty flux register
  LevelFluxRegister m_lambda_flux_reg;

  /// conserved scalar flux register
  Vector<LevelFluxRegister*> m_scal_fluxreg_ptrs;

  /// does quadratic coarse-fine interpolation for velocities
  QuadCFInterp m_velCFInterp;

  /// AMRPoissonOp for velocities
  VelocityAMRPoissonOp m_velocityAMRPoissonOp;

  /// AMRPoissonOp for advected scalars (for diffusion)
  AMRPoissonOp m_scalarsAMRPoissonOp;

  /// levelsolver for viscous solves for velocity
  Tuple< RefCountedPtr< AMRMultiGrid<LevelData<FArrayBox> > >, SpaceDim> m_velMGsolverPtrs;

  /// levelsolver for advectDiffuseScalar solves
  Vector< RefCountedPtr< AMRMultiGrid<LevelData<FArrayBox> > > > m_scalMGsolverPtrs;

  /// TGA solver for computeUstar():  LevelTGA has no default constructor
  Tuple< RefCountedPtr<LevelTGA>, SpaceDim> m_TGAsolverPtrs;

  /// TGA solver for advectDiffuseScalar():  LevelTGA has no default constructor
  Vector< RefCountedPtr<LevelTGA> > m_TGAsolverScalPtrs;

  /// Lambda advection
  PatchGodunov m_patchGodLambda;

  /// Velocity advection
  PatchGodunov m_patchGodVelocity;

  /// More scalar advection
  Vector<PatchGodunov*> m_patchGodScalars;

  /// refinement threshold for regridding
  Real m_refine_threshold;

  /// true if this is the finest extant level
  bool m_finest_level;

  /// true if this level is an empty level (no grids)
  bool m_is_empty;

  /// internal indicator whether post-regrid smoothing stuff has been done
  bool m_regrid_smoothing_done;

  /// physical BC object
  PhysBCUtil* m_physBCPtr;
};

#endif
