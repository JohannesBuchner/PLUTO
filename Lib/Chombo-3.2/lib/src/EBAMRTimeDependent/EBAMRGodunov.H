#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRGODUNOV_H_
#define _EBAMRGODUNOV_H_

#include "AMRLevel.H"
#include "EBCellFAB.H"
#include "BaseIVFAB.H"
#include "LevelData.H"
#include "EBLevelRedist.H"
#include "EBCoarToCoarRedist.H"
#include "EBCoarToFineRedist.H"
#include "EBFineToCoarRedist.H"
#include "EBCoarseAverage.H"
#include "EBPWLFineInterp.H"
#include "EBFluxRegister.H"
#include "EBLevelGodunov.H"
#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "NamespaceHeader.H"

///interface for load balancing
typedef int(*LoadBalanceFunc)(Vector<int>&         a_procs,
                               const Vector<Box>&   a_boxes,
                               const ProblemDomain& a_domain,
                               bool a_verbose);


///
/**
 */
class EBAMRGodunov : public AMRLevel
{
public:

  EBAMRGodunov();

  virtual ~EBAMRGodunov();

  void sumConserved(Real& a_sumcons,
                    const int& a_ivar) const;

  virtual void define(AMRLevel*  a_coarser_level_ptr,
              const Box& a_problem_domain,
              int        a_level,
              int        a_ref_ratio);

  virtual void define(AMRLevel*  a_coarser_level_ptr,
              const ProblemDomain& a_problem_domain,
              int        a_level,
              int        a_ref_ratio);

  /// advance by one timestep
  virtual Real advance();

  /// things to do after a timestep
  virtual void postTimeStep();

  /// setup fine->coarse stuff
  /**
     stuff that needs to be setup from the finer
     level.  A bunch of objects depend upon the layouts
     from both levels and the finer level changes more
     often from regrid so this needs to be called from the finer
     level
  */
  virtual void syncWithFineLevel();

  void dumpDebug();

  /// create tags
  virtual void tagCells(IntVectSet& a_tags) ;

  /// create tags at initialization
  virtual void tagCellsInit(IntVectSet& a_tags) ;

  /// regrid
  virtual void regrid(const Vector<Box>& a_new_grids);

  /// initialize grids
  virtual void initialGrid(const Vector<Box>& a_new_grids);

  /// initialize data
  virtual void initialData();

  /// things to do after initialization
  virtual void postInitialize();

  /// object for physics integration
  virtual void patchGodunov(const EBPatchGodunovFactory* const a_ebpatchGodunovFactory);
  void dumpDebug(const string& a_debstring);
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

  void writePlotHeaderOld(HDF5Handle& a_handle) const;

  /// write plot file data for this level
  virtual void writePlotLevel(HDF5Handle& a_handle) const;

  void writePlotLevelOld(HDF5Handle& a_handle) const;

#endif

  ///
  virtual Real getDt() const;

  /// compute dt
  virtual Real computeDt();

  /// compute dt with initial data
  virtual Real computeInitialDt();

  /// CFL
  virtual void CFL(Real a_cfl);

  ///
  virtual void doSmushing(bool a_doSmushing);

  ///
  virtual   void doRZCoords(bool a_doRZCoords);

  ///
  virtual   void hasSourceTerm(bool a_hasSourceTerm);

  ///redistribution radius
  virtual   void redistRadius(int a_redistRad);

  /// longest length of the domain in real space
  virtual   void domainLength(RealVect a_domainLength);

  /// refinement threshold
  virtual   void refinementThreshold(Real a_refine_thresh);

  /// tag buffer size
  virtual   void tagBufferSize(int a_tag_buffer_size);

  virtual   void useMassRedistribution(bool a_useMassRedist);

  virtual void levelSetup();

  const EBISLayout& getEBISLayout() const;
  virtual void fillConsAndPrim(LevelData<EBCellFAB>& a_data) const;

  void tagAll(bool a_tagAll);

  // state vector at old time
  LevelData<EBCellFAB>& getStateOld();

  // state vector at new time
  LevelData<EBCellFAB>& getStateNew();

  static int s_NewPlotFile;

  static void setLoadBalance(LoadBalanceFunc a_func)
  {
    s_loadBalance      = a_func;
    s_isLoadBalanceSet = true;
  }
protected:

  static LoadBalanceFunc   s_loadBalance;
  static bool              s_isLoadBalanceSet;
  bool m_tagAll;
  bool m_useMassRedist;
  Box m_domainBox;
  Real calculateMass() const;

  EBAMRGodunov* getCoarserLevel() const;
  EBAMRGodunov* getFinerLevel() const;

  DisjointBoxLayout m_grids;
  EBISLayout m_ebisl;

  // state vector at old time
  LevelData<EBCellFAB> m_stateOld;

  // state vector at new time
  LevelData<EBCellFAB> m_stateNew;

  // number of components of m_state
  int m_nComp;

  int m_nGhost ;
  // names of components
  Vector<string> m_stateNames;
  Vector<string> m_primNames;

  // grid spacing
  RealVect m_origin;
  RealVect m_dx;
  RealVect m_aspect;

  // domain length
  RealVect m_domainLength;

  // cfl number
  Real m_cfl;

  // new time step
  Real m_dtNew;

  // averaging from coarse to fine level
  EBCoarseAverage m_ebCoarseAverage;

  // interpolation from fine to coarse level
  EBPWLFineInterp m_ebFineInterp;

  // level solver
  EBLevelGodunov m_ebLevelGodunov;

  // flux register
  EBFluxRegister m_ebFluxRegister;

  // patch integrator
  EBPatchGodunov* m_ebPatchGodunov;

  // patch integrator
  const EBPatchGodunovFactory* m_ebPatchGodunovFactory;

  // refinement threshold for gradient
  Real m_refineThresh;

  //redistribution mass
  LevelData<BaseIVFAB<Real> > m_massDiff;
  LayoutData<IntVectSet> m_sets;

  RedistStencil m_redStencil;
  EBCoarToFineRedist m_ebCoarToFineRedist;
  EBCoarToCoarRedist m_ebCoarToCoarRedist;
  EBFineToCoarRedist m_ebFineToCoarRedist;

  // tag buffer size
  int m_tagBufferSize;
  int m_redistRad;

  bool m_doRZCoords;
  bool m_hasSourceTerm;
  bool m_doSmushing;
  bool m_hasCoarser;
  bool m_hasFiner;
  bool m_isDefined;


private:
  //disallowed for all the usual reasons
  void operator=(const EBAMRGodunov& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBAMRGodunov(const EBAMRGodunov& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
