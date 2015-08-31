#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMR_LEVEL_GODUNOV_H_
#define _AMR_LEVEL_GODUNOV_H_

#include "AMRLevel.H"
#include "CoarseAverage.H"

#include "LevelPluto.H"

#include "NamespaceHeader.H" 

/// AMR Pluto
/**
 */
class AMRLevelPluto : public AMRLevel
{
public:
  /// Constructor
  /**
   */
  AMRLevelPluto();

  /// Destructor
  /**
   */
  virtual ~AMRLevelPluto();

  /// Define the parameters the object needs
  /**
   */
  void defineParams(const Real&              a_cfl,
                    const Real&              a_domainLength,
                    const int&               a_verbosity,
                    const Real&              a_refineThresh,
                    const int&               a_tagBufferSize,
                    const Real&              a_initialDtMultiplier,
                    const PatchPluto*  const a_patchPluto);

  /// This instance should never get called - historical
  /**
   */
  virtual void define(AMRLevel*  a_coarserLevelPtr,
                      const Box& a_problemDomain,
                      int        a_level,
                      int        a_refRatio);

  /// Define new AMR level
  /**
   */
  virtual void define(AMRLevel*            a_coarserLevelPtr,
                      const ProblemDomain& a_problemDomain,
                      int                  a_level,
                      int                  a_refRatio);


  /// Advance by one timestep
  /**
   */
  virtual Real advance();

  virtual Real getDlMin();

  /// Things to do after a timestep
  /**
   */
  virtual void postTimeStep();

  /// Create tags for regridding
  /**
   */
  virtual void tagCells(IntVectSet& a_tags) ;

  /// Create tags at initialization
  /**
   */
  virtual void tagCellsInit(IntVectSet& a_tags) ;

  /// Set up data on this level after regridding
  /**
   */
  virtual void regrid(const Vector<Box>& a_newGrids);

  /// Initialize grids
  /**
   */
  virtual void initialGrid(const Vector<Box>& a_newGrids);

  /// Initialize data
  /**
   */
  virtual void initialData();

  /// Things to do after initialization
  /**
   */
  virtual void postInitialize();

#ifdef CH_USE_HDF5
  /// Write checkpoint header
  /**
   */
  virtual void writeCheckpointHeader(HDF5Handle& a_handle) const;

  /// Write checkpoint data for this level
  /**
   */
  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;

  /// Read checkpoint header
  /**
   */
  virtual void readCheckpointHeader(HDF5Handle& a_handle);

  /// Read checkpoint data for this level
  /**
   */
  virtual void readCheckpointLevel(HDF5Handle& a_handle);

  /// Write plotfile header
  /**
   */
  virtual void writePlotHeader(HDF5Handle& a_handle) const;

  /// Write plotfile data for this level
  /**
   */
  virtual void writePlotLevel(HDF5Handle& a_handle) const;

 /// Expressions for Visit
 /**
  */
  virtual void DefineExpressions(HDF5HeaderData& a_holder) const;

#endif

  /// Returns the dt computed earlier for this level
  /**
   */
  virtual Real computeDt();

  /// Compute dt using initial data
  /**
   */
  virtual Real computeInitialDt();

  /// Mark split/unsplit cells
  /**
   */
  void mark_split(const DisjointBoxLayout& finerLevelDomain);

  ///
  const LevelData<FArrayBox>& getStateNew() const;

  ///
  const LevelData<FArrayBox>& getStateOld() const;

  ///
  bool allDefined() const;

protected:
  // Create a load-balanced DisjointBoxLayout from a collection of Boxes
  DisjointBoxLayout loadBalance(const Vector<Box>& a_grids);

  // Setup menagerie of data structures
  void levelSetup();

  // Get the next coarser level
  AMRLevelPluto* getCoarserLevel() const;

  // Get the next finer level
  AMRLevelPluto* getFinerLevel() const;

  // Conserverd state, U, at old and new time
  LevelData<FArrayBox> m_UOld,m_UNew;

  // Union of FArrayBoxes marking the split/unsplit cells of this level
  LevelData<FArrayBox> m_split_tags;

  // CFL number
  Real m_cfl;

  // Grid spacing
  Real m_dx;

  // Interpolation from fine to coarse level
  FineInterpPluto m_fineInterp;

  // Averaging from coarse to fine level
  CoarseAverage m_coarseAverage;

  // New time step
  Real m_dtNew;

  // Number of converved states
  int m_numStates;

  // Names of conserved states
  Vector<string> m_ConsStateNames;

  // Names of primitive states
  Vector<string> m_PrimStateNames;

  // Number of ghost cells (in each direction)
  int m_numGhost;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;

  // Level integrator
  LevelPluto m_levelPluto;

  // Flux register
  LevelFluxRegister m_fluxRegister;

  // Patch integrator and its factory
  PatchPluto* m_patchPluto;
  PatchPluto* m_patchPlutoFactory;

  // Refinement threshold for gradient
  Real m_refineThresh;

  // Tag buffer size
  int m_tagBufferSize;

  // Flag coarser and finer levels
  bool m_hasCoarser;
  bool m_hasFiner;

  // Grid layout for this level
  DisjointBoxLayout m_grids;

  // True if all the parameters for this object are defined
  bool m_paramsDefined;

private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelPluto& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelPluto(const AMRLevelPluto& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"

#endif
