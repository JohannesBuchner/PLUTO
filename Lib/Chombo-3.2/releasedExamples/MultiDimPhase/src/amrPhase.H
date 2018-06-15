#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRPHASE_H_
#define _AMRPHASE_H_


//#define CH_SPACEDIM 2

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "RealVect.H"
#include "PatchAdvection.H"
#include "OldPhysIBC.H"
#include "Vector.H"
#include "ReductionCopier.H"
#include "SpreadingCopier.H"

// needed for multidim emulation -- unused in
// true multidim
#ifndef CH_MULTIDIM
#include "PoissonOp.H"
#include "AMRSolver.H"
#endif

#include "UsingNamespace.H"

#define ADVECT_GROW 4

#ifdef CH_USE_DOUBLE
#define TIME_EPS 1.0e-10
#else
#define TIME_EPS 1.0e-5
#endif

/// class to manage non-subcycled porous media AMR computation
/**
 */
class amrPHASE
{

public:
  /// Default constructor
  /** At the moment, there is only one constructor, which defines itself
      based on an inputs file.  (This may change, of course)
  */
  amrPHASE();

  /// destructor
  /** destructor
   */
  ~amrPHASE();

  /// set default values before initialization
  void setDefaults();

  /// initializes object based on inputs data passed in through parmParse
  /** initializes new object based on data in ParmParse database.  This
      may change to a more explicit initialization if it's thought to be
      a good idea.
  */
  void initialize();


  /// advance solution until either max_time or max_step are reached
  /**
   */
  void run(Real a_max_time, int a_max_step);

  /// is this object defined and initialized?
  /**
   */
  bool isDefined() const;

#ifdef CH_USE_HDF5

  /// write hdf5 plotfile to the standard location
  /** writes plotfile to the standard location; filename is
      <plot_root>.<step>.<DIM>d.hdf5
   */
  void writePlotFile();


  /// write checkpoint file out for later restarting
  /** writes checkpoint data to the standard location; filename is
      <check_root>.<step>.<DIM>d.hdf5
   */
  void writeCheckpointFile() const;


  /// read checkpoint file for restart
  /** read checkpoint data from file pointed to by a_handle
   */
  void readCheckpointFile(HDF5Handle& a_handle);

  /// set up for restart
  void restart(string& a_restart_file);

#endif

protected:
  /// compute one timestep
  void timeStep(Real a_dt);


  /// compute velocity
  void computePhaseVelocity(Vector<LevelData<FArrayBox>* >& a_CCvelocity,
                            Vector<LevelData<FArrayBox>* >& a_potential,
                            Vector<LevelData<FArrayBox>* >& a_phi,
                            Real a_time);

  /// advect scalar
  /// do scalar advection -- assumes all BC's already done
  void advectScalar(Vector<LevelData<FArrayBox>* >& a_new_scal,
                    Vector<LevelData<FArrayBox>* >& a_old_scal,
                    Vector<LevelData<FArrayBox>* >& a_CC_vel,
                    OldPhysIBC* a_scalIBC,
                    Real a_dt);

  /// fill cell-centered velocity field, including ghost cells
  void fillCCvelocity(Vector<LevelData<FArrayBox>* >& a_ccVel,
                      Vector<LevelData<FArrayBox>* >& a_CC_vel,
                      Real a_oldTime);


  /// do regridding
  void regrid();

  /// compute tags for regridding
  void tagCells(Vector<IntVectSet>& a_tags,
                Vector<IntVectSet>& a_1Dtags);

  /// compute tags for the level a_level
  void tagCellsLevel(IntVectSet& a_tags,
                     IntVectSet& a_1Dtags,
                     int a_level);

  /// compute tags at initial time
  void tagCellsInit(Vector<IntVectSet>& a_tags,
                    Vector<IntVectSet>& a_1Dtags);

  /// initialize grids at initial time
  void initGrids(int a_finest_level);

  /// set up grids from grids file
  void setupFixedGrids(const std::string& a_gridFile);

  /// initialize data on hierarchy
  void initData(Vector<LevelData<FArrayBox>*>& a_phi,
                Vector<LevelData<FArrayBox>*>& a_CCvelocity);

  /// compute timestep
  Real computeDt();

  /// compute timestep at initial time
  Real computeInitialDt();

  /// gravitational constant
  Real m_G;

  /// max number of levels
  int m_max_level;

  /// current finest level
  int m_finest_level;

  /// blocking factor
  int m_block_factor;

  /// grid efficiency
  Real m_fill_ratio;

  /// proper nesting radius
  int m_nesting_radius;

  /// max box size
  int m_max_box_size;

  /// regrid interval
  int m_regrid_interval;

  /// tagging value (undivided gradient threshold for regridding)
  Real m_tagging_val;

  /// amount to buffer tags used in regridding
  int m_tags_grow;

  /// number of components
  int m_num_comps;

  /// refinement ratios
  Vector<int>  m_refinement_ratios;

  /// cell spacing at each level
  Vector<Real> m_amrDx;

  /// problem domains at each level
  Vector<ProblemDomain> m_amrDomains;

  /// reduced problem domains at each level
  Vector<ProblemDomain> m_amrReducedDomains;

  // problem domain size
  RealVect m_domainSize;

  // location of origin
  RealVect m_origin;

  /// current grids
  Vector<DisjointBoxLayout> m_amrGrids;

  /// reduced grids
  Vector<DisjointBoxLayout> m_amrReducedGrids;

  // copiers for reducing from  phase space to physical space
  Vector<ReductionCopier*> m_reduceCopiers;

  // copiers for spreading from physical space back to phase space
  Vector<SpreadingCopier*> m_spreadCopiers;

  /// book-keeping; keeps strack of number of cells per level
  Vector<int> m_num_cells;

  /// current time
  Real m_time;

  /// most recent timestep
  Real m_dt;

  /// timestep scaling
  Real m_cfl;

  /// cfl number for initial timestep (useful if initial data needs small cfl)
  Real m_initial_cfl;

  /// maximum amount cfl number may grow in one timestep
  Real m_max_dt_grow;

  // current step
  int m_cur_step;

  /// current old-time scalar
  Vector<LevelData<FArrayBox>*> m_old_phi;

  /// new-time scalar
  Vector<LevelData<FArrayBox>*> m_new_phi;

  /// cell-centered velocity
  Vector<LevelData<FArrayBox>* > m_CCvelocity;

  /// pressure
  Vector<LevelData<FArrayBox>* > m_potential;

#ifndef CH_MULTIDIM
  // physical boundary conditions for phi
  DomainGhostBC m_potentialBC;

  /// PoissonOp upon which the solver is built
  PoissonOp m_poissOp;

  /// amr solver for potential solve
  AMRSolver m_solver;

#endif

#if 0
  // physical boundary conditions for phi
  D1::DomainGhostBC m_potentialBC;

  /// PoissonOp upon which the solver is built
  D1::PoissonOp m_poissOp;
#endif

  OldPhysIBC* m_scalIBC;

  // advection objects
  Vector<PatchAdvection*> m_patchGod;


  /// how verbose should we be?
  int m_verbosity;

  /// is this object initialized?
  bool m_is_defined;

  ///
  bool m_do_restart;

  /// if starting from a restart, timestep of restart
  int m_restart_step;

  /// initially refine entire domain (rather than generating tags)
  /** If true, the entire domain is refined to the finest possible
      resolution at the initial step, rather than generating tags, etc.
      This can be useful for random initial conditions
  */
  bool m_refine_initial_domain;

  /// which heat solver formulation to use
  /** which heat equation solver to use: 0 = backward euler, 1 = TGA
   */
  int m_solver_type;

  ///
  string m_plot_prefix;

  ///
  string m_check_prefix;

  ///
  int m_plot_interval;

  ///
  int m_check_interval;

};


#endif


