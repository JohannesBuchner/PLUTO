#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELADVECTDIFFUSE_H_
#define _AMRLEVELADVECTDIFFUSE_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "LevelFluxRegister.H"

#include "LevelAdvect.H"
#include "PhysIBC.H"
#include "ProblemDomain.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "RealVect.H"
#include "LevelFluxRegister.H"
#include "DisjointBoxLayout.H"
#include "AdvectionFunctions.H"
#include "BCFunc.H"
#include "LevelTGA.H"
#include "AMRPoissonOp.H"
#include "BiCGStabSolver.H"
#include "RelaxSolver.H"

#include "NamespaceHeader.H"

/// AMRLevel for advection-diffusion
/**
 */
class AMRLevelAdvectDiffuse : public AMRLevel
{
public:

  /// Default constructor
  AMRLevelAdvectDiffuse()
  {
    m_isDefined = false;
  }

  /// Full constructor. Arguments are same as in define()
  AMRLevelAdvectDiffuse(const AdvectPhysics&        a_gphys,
                        AdvectionVelocityFunction   a_advFunc,
                        BCHolder                    a_bcFunc,
                        const Real&                 a_cfl,
                        const Real&                 a_domainLength,
                        const Real&                 a_refineThresh,
                        const int&                  a_tagBufferSize,
                        const Real&                 a_initialDtMultiplier,
                        const bool&                 a_useLimiting,
                        const Real&                 a_nu)
  {
    define(a_gphys, a_advFunc, a_bcFunc, a_cfl, a_domainLength, a_refineThresh,
           a_tagBufferSize, a_initialDtMultiplier, a_useLimiting, a_nu);
  }

  void makeDiffusiveSource(LevelData<FArrayBox>& a_diffusiveSource);

  /// Defines this AMRLevelAdvectDiffuse
  void define(/// advection physics class
              const AdvectPhysics&        a_gphys,
              /// velocity function providing the advection velocity
              AdvectionVelocityFunction   a_advFunc,
              /// boundary condition class for diffusion solve
              BCHolder                    a_bcFunc,
              /// CFL number
              const Real&                 a_cfl,
              /// physical length of domain
              const Real&                 a_domainLength,
              /// undivided gradient size over which a cell will be tagged for refinement
              const Real&                 a_refineThresh,
              /// number of buffer cells around each tagged cell that will also be tagged
              const int&                  a_tagBufferSize,
              /// CFL number at beginning of calculation
              const Real&                 a_initialDtMultiplier,
              /// whether to use van Leer limiting
              const bool&                 a_useLimiting,
              /// diffusion coefficient
              const Real&                 a_nu);

  ///
  virtual ~AMRLevelAdvectDiffuse();

  /// Never called: historical
  virtual void define(AMRLevel*  a_coarserLevelPtr,
                      const Box& a_problemDomain,
                      int        a_level,
                      int        a_refRatio)
  {
    MayDay::Error("never called--historical");
  }

  /// Define new AMRLevelAdvectDiffuse from coarser
  /**
   */
  virtual void define(
                      AMRLevel*            a_coarserLevelPtr,
                      const ProblemDomain& a_problemDomain,
                      int                  a_level,
                      int                  a_refRatio);

  /// Advance by one timestep
  virtual Real advance();

  /// Things to do after a timestep
  virtual void postTimeStep();

  /// Create tags for regridding
  virtual void tagCells(IntVectSet& a_tags) ;

  /// Create tags at initialization
  virtual void tagCellsInit(IntVectSet& a_tags) ;

  /// Set up data on this level after regridding
  virtual void regrid(const Vector<Box>& a_newGrids);

  /// Initialize grids
  virtual void initialGrid(const Vector<Box>& a_newGrids);

  /// Initialize data
  virtual void initialData();

  /// Things to do after initialization
  virtual void postInitialize();

#ifdef CH_USE_HDF5
  /// Write checkpoint header
  virtual void writeCheckpointHeader(HDF5Handle& a_handle) const;

  /// Write checkpoint data for this level
  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;

  /// Read checkpoint header
  virtual void readCheckpointHeader(HDF5Handle& a_handle);

  /// Read checkpoint data for this level
  virtual void readCheckpointLevel(HDF5Handle& a_handle);

  /// Write plotfile header
  virtual void writePlotHeader(HDF5Handle& a_handle) const;

  /// Write plotfile data for this level
  virtual void writePlotLevel(HDF5Handle& a_handle) const;
#endif

  /// Returns the dt computed earlier for this level
  virtual Real computeDt();

  /// Compute dt using initial data
  virtual Real computeInitialDt();

  //for convergence tests
  LevelData<FArrayBox>& getStateNew()
  {
    return m_UNew;
  }

  //for convergence tests
  LevelData<FArrayBox>& getStateOld()
  {
    return m_UOld;
  }
protected:
  void setSolverCoef(Real a_alpha, Real a_beta);
  void interpolateInTime(LevelData<FArrayBox>&          a_interp,
                         const LevelData<FArrayBox>&    a_old,
                         const LevelData<FArrayBox>&    a_new,
                         Real a_time, Real a_tOld, Real a_tNew);

  void getHierarchyAndGrids(Vector<AMRLevelAdvectDiffuse*>&        a_hierarchy,
                            Vector<DisjointBoxLayout>&             a_grids,
                            Vector<int>&                           a_refRat,
                            ProblemDomain&                         a_lev0Dom,
                            Real&                                  a_lev0Dx);

  void doImplicitReflux();

  void printRefluxRHSMax(const std::string& a_preflix) ;

  Real diffusiveAdvance(LevelData<FArrayBox>& a_diffusiveSource);
  void defineSolvers();
  BCHolder                    m_bcFunc;
  // solver to do level  diffusion and flux register interaction
  static RefCountedPtr<LevelTGA>                                  s_diffuseLevTGA;
  // one amrmultigrid to rule them all and in darkness bind them
  static RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >      s_diffuseAMRMG;
  // operator factory
  static RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > s_diffuseOpFact;
  // bottom solver
  static BiCGStabSolver<LevelData<FArrayBox> >                    s_botSolver;
  // static RelaxSolver<LevelData<FArrayBox> >                       s_botSolver;

  //embarrassing but needed function
  void getCoarseDataPointers(LevelData<FArrayBox>** a_coarserDataOldPtr,
                             LevelData<FArrayBox>** a_coarserDataNewPtr,
                             LevelFluxRegister**    a_coarserFRPtr,
                             LevelFluxRegister**    a_finerFRPtr,
                             Real& a_tCoarserOld,
                             Real& a_tCoarserNew);

  void fillAdvectionVelocity();
  // Setup menagerie of data structures
  void levelSetup();

  // Get the next coarser level
  AMRLevelAdvectDiffuse* getCoarserLevel() const;

  // Get the next finer level
  AMRLevelAdvectDiffuse* getFinerLevel() const;

  // Conserved state, U, at old and new time
  LevelData<FArrayBox> m_UOld,m_UNew,m_dU;
  LevelData<FluxBox>   m_advVel;

  Vector<string> m_stateNames;
  bool m_isDefined;
  Real                 m_cfl;
  Real                 m_domainLength;
  Real                 m_refineThresh;
  int                  m_tagBufferSize;
  Real                 m_initialDtMultiplier;
  bool                 m_useLimiting;
  Real                 m_nu;
  Real                 m_dx;
  RefCountedPtr<AdvectPhysics> m_advPhys;
  AdvectionVelocityFunction    m_advFunc;

  FineInterp m_fineInterp;
  CoarseAverage m_coarseAverage;
  Real m_dtNew;
  int m_numGhost;
  LevelAdvect m_levelGodunov;
  LevelFluxRegister m_fluxRegister;

  bool m_hasCoarser;
  bool m_hasFiner;
  bool m_doImplicitReflux;
  bool m_hasDiffusion;

  DisjointBoxLayout m_grids;

private:

  // Disallowed for all the usual reasons
  void operator=(const AMRLevelAdvectDiffuse&);
  AMRLevelAdvectDiffuse(const AMRLevelAdvectDiffuse&);
};

#include "NamespaceFooter.H"

#endif
