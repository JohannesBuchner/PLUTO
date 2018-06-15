#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELPOLYTROPICGAS_H_
#define _AMRLEVELPOLYTROPICGAS_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "LevelFluxRegister.H"
#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"

#include "LevelGodunov.H"
#include "PhysIBC.H"

#include "NamespaceHeader.H"

/// AMR Godunov
/**
 */
class AMRLevelPolytropicGas: public AMRLevel
{
public:
  /// Constructor
  /**
   */
  AMRLevelPolytropicGas();

  /// Destructor
  /**
   */
  virtual ~AMRLevelPolytropicGas();

  /// Define the parameters the object needs
  /**
   */
  void defineParams(/// CFL number
                    const Real&                 a_cfl,
                    /// physical dimension of the longest side of the domain
                    const Real&                 a_domainLength,
                    /// verbosity level
                    const int&                  a_verbosity,
                    /// undivided gradient size over which a cell will be tagged for refinement
                    const Real&                 a_refineThresh,
                    /// number of buffer cells around each tagged cell that will also be tagged
                    const int&                  a_tagBufferSize,
                    /// CFL number at beginning of calculation
                    const Real&                 a_initialDtMultiplier,
                    /// physics class
                    const GodunovPhysics* const a_godunovPhysics,
                    /// order of the normal predictor: 1 for PLM, 2 for PPM
                    const int&                  a_normalPredOrder,
                    /// whether to use 4th-order slopes
                    const bool&                 a_useFourthOrderSlopes,
                    /// whether to apply slope limiting to primitive variables
                    const bool&                 a_usePrimLimiting,
                    /// whether to apply slope limiting to characteristic variables
                    const bool&                 a_useCharLimiting,
                    /// whether to apply slope flattening
                    const bool&                 a_useFlattening,
                    /// whether to apply artificial viscosity
                    const bool&                 a_useArtificialViscosity,
                    /// artificial viscosity coefficient
                    const Real&                 a_artificialViscosity,
                    /// whether to apply source term
                    const bool&                 a_useSourceTerm,
                    /// scaling of source term
                    const Real&                 a_sourceTermScaling,
                    /// whether to apply 4th-order limiter
                    const bool&                 a_highOrderLimiter);

  /// This instance should never get called - historical
  /**
   */
  virtual void define(AMRLevel*  a_coarserLevelPtr,
                      const Box& a_problemDomain,
                      int        a_level,
                      int        a_refRatio);

  /// Define new AMRLevelPolytropicGas
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

  /// Things to do after a timestep
  /**
   */
  virtual void postTimeStep();

  /// Create tags for regridding
  /**
   */
  virtual void tagCells(IntVectSet& a_tags);

  /// Create tags at initialization
  /**
   */
  virtual void tagCellsInit(IntVectSet& a_tags);

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
#endif

  /// Returns the dt computed earlier for this level
  /**
   */
  virtual Real computeDt();

  /// Compute dt using initial data
  /**
   */
  virtual Real computeInitialDt();

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
  AMRLevelPolytropicGas* getCoarserLevel() const;

  // Get the next finer level
  AMRLevelPolytropicGas* getFinerLevel() const;

  // Conserved state, U, at old and new time
  LevelData<FArrayBox> m_UOld,m_UNew;

  // CFL number
  Real m_cfl;

  // Grid spacing
  Real m_dx;

  // Interpolation from coarse to fine level
  FineInterp m_fineInterp;

  // Averaging from fine to coarse level
  CoarseAverage m_coarseAverage;

  // New time step
  Real m_dtNew;

  // Number of conserved states
  int m_numStates;

  // Names of conserved states
  Vector<string> m_stateNames;

  // Number of ghost cells (in each direction)
  int m_numGhost;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;

  // Level integrator
  LevelGodunov m_levelGodunov;

  // Flux register
  LevelFluxRegister m_fluxRegister;

  // Pointer to the class defining the physics of the problem
  GodunovPhysics* m_gdnvPhysics;

  // Order of the normal predictor (1 -> PLM, 2-> PPM)
  int m_normalPredOrder;

  // Use 4th order slope computations (otherwise, use 2nd order)
  bool m_useFourthOrderSlopes;

  // Do slope limiting in the primitive or characteristic variables, respect.
  bool m_usePrimLimiting;
  bool m_useCharLimiting;

  // Do slope flattening - MUST BE USING 4th order slopes
  bool m_useFlattening;

  // Apply artificial viscosity of a set value
  bool m_useArtificialViscosity;
  Real m_artificialViscosity;

  // Supply a source term to the computation
  bool m_useSourceTerm;
  Real m_sourceTermScaling;

  // Use a high-order limiter?
  bool m_highOrderLimiter;

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
  void operator=(const AMRLevelPolytropicGas&);
  AMRLevelPolytropicGas(const AMRLevelPolytropicGas&);
};

#include "NamespaceFooter.H"

#endif
