#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRNOSUBCYCLE_H_
#define _EBAMRNOSUBCYCLE_H_

#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "EBIBCFactory.H"
#include "AMRTGA.H"
#include "EBBackwardEuler.H"
#include "EBCompositeMACProjector.H"
#include "EBCompositeCCProjector.H"
#include "EBAMRPoissonOp.H"
#include "EBLevelAdvect.H"
#include "EBCoarseAverage.H"
#include "EBSimpleSolver.H"
#include "EBFastFR.H"

#include "NamespaceHeader.H"

///
/**
  Agglomeration of the many parameters of an AMR calc
*/
class AMRParameters
{
public:
  ///almost everything gets a default
  AMRParameters()
  {
    m_domainLength = 1.0;
    m_maxLevel     = 0;
    m_checkpointInterval = -1;
    m_plotInterval       = -1;
    m_maxBoxSize         = 1024;
    m_fillRatio          = 0.7;
    m_blockFactor        = 8;
    m_regridInterval     = -1;
    m_maxDtGrow          = 1.2345678e9;
    m_cfl                = 0.5;
    m_initCFL            = 0.1;
    m_refineThreshold    = 0.1;
    m_verbosity          = 1;
    m_nestingRadius      = 2;
    m_maxDt              = 1.2345678e9;
    m_gphiIterations     = 0;
    m_initIterations     = 0;
    m_copyOverOld        = false;
    m_useLimiting        = true;
    m_numFilterIterations = 1;
    m_chopXTagsHi         = false;
    m_relaxType           = 2;
    m_hang                = 1.0e-1;
    m_tolerance           = 1.0e-12;
    m_normThresh          = 1.0e-15;
    m_iterMax             = 100;
    m_mgCycle             = 1;
    m_numSmooth           = 4;
    m_numBotSmooth        = 64;
    m_doRegridSmoothing   = false;
    m_stokesFlow          = false;
    m_useStokesDt         = false;
    m_subtractOffMean     = false;
    m_tagOnScalarToo      = false;
    m_tagShrinkDomain     = 0;
    m_flowDir             = 0;
    m_initPoiseData       = false;
    m_doSteadyState       = false;
    m_orderTimeIntegration = 2;
  }

  ///
  ~AMRParameters()
  {;}

  //parameters for projections
  int m_numSmooth, m_numBotSmooth, m_iterMax, m_mgCycle, m_relaxType;
  Real m_hang, m_tolerance, m_normThresh;

  int m_numFilterIterations;
  int m_gphiIterations;

  ///
  /**
     This is always true except for truncation error test
     where we use new-old/dt to compute  du/dt.
  */
  bool m_copyOverOld;

  ///number of timesteps we iterate for gphi on startup
  int  m_initIterations;
  bool m_tagOnScalarToo; //deprecated but too much work to remove
  bool m_useLimiting;
  /// length of the domain in the 0 direction
  Real                  m_domainLength;
  Real                  m_maxDt;

  /// number of  highest AMR Level (0 - numLevels -1)
  int                   m_maxLevel;

  /// number of steps between checkpoints (-1 means no checkpoints)
  int                   m_checkpointInterval;

  /// number of steps between plotfiles (-1 means no plotfiles)
  int                   m_plotInterval;

  ///  1-D maximum size of boxes in grid hierarchy
  int                   m_maxBoxSize;

  /// minimum percentage of tagged cells per box ala B-R
  Real                  m_fillRatio;

  /// blocking factor ala B-R
  int                   m_blockFactor;

  ///time steps between regrids
  int                   m_regridInterval;

  ///maximum growth factor for time step
  Real                  m_maxDtGrow;

  ///CFL number for time steps after the first one
  Real                  m_cfl;

  ///CFL number for the first time step
  Real                  m_initCFL;

  /// factor on which to refine
  Real                  m_refineThreshold;

  /// verbosity of amr class
  int                   m_verbosity;

  ///
  int                   m_nestingRadius;

  ///
  int                   m_tagBuffer;

  ///number of cells to shrink domain by so that there is no refinement at side walls
  int                   m_tagShrinkDomain;

  ///flow direction
  int                   m_flowDir;

  /// refinement ratio between levels.  0 is the refrat between 0 and 1.
  Vector<int>           m_refRatio;

  ///
  bool m_doRegridSmoothing;

  bool m_chopXTagsHi;

  bool m_stokesFlow;
  bool m_useStokesDt;
  bool m_subtractOffMean;
  bool m_initPoiseData;
  bool m_doSteadyState;
  int  m_orderTimeIntegration;
};


/// Class to manage non-subcycled advance for EBAMRINS.
/**
*/
class EBAMRNoSubcycle
{

public:


  // weak construction is bad
  EBAMRNoSubcycle()
  {
    MayDay::Error("invalid constructor");
  }

  ///
  /**
       Note that the EBAMRNoSubcycle class defines its own grid hierarchy
       Initial and boundary conditions come through the EBIBCFactory.
       If fixedGrids != NULL, regridding gets turned off and the input
       grids are used
   */
  EBAMRNoSubcycle(const AMRParameters&      a_params,
                  const EBIBCFactory&       a_IBC,
                  const ProblemDomain&      a_coarsestDomain,
                  Real                      a_viscosity,
                  const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance());

  ///
  void
  setSolverParams(int  a_numSmooth,
                  int  a_itermax,
                  int  a_mgcycle,
                  Real a_hang,
                  Real a_tolerance);

  void
  normalVelocityPredictor(Vector<LevelData<EBFluxFAB> *>&                       a_advVel,
                          Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelLo,
                          Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelHi,
                          const Vector<LevelData<EBCellFAB> *>&                 a_velo);

  void
  computeAdvectiveDerivative(Vector<LevelData<EBCellFAB>* >                     &  a_uDotDelU,
                             Vector<LevelData<EBFluxFAB>* >                     &  a_macAdvVel,
                             Vector<LevelData<EBFluxFAB>* >                     &  a_macScalar,
                             Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelLo,
                             Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelHi,
                             Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredScalarLo,
                             Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredScalarHi,
                             bool                                                  a_nonConsOnly     = false,
                             bool                                                  a_consOnly        = false,
                             Vector<RefCountedPtr<EBLevelAdvect> >              *  a_ebLevAd = NULL);

  //fixes timestep to input
  void useFixedDt(Real a_dt);

  /// destructor.
  virtual ~EBAMRNoSubcycle();

  /// set up grids, etc for fixed-hiearchy run
  void setupForFixedHierarchyRun(const Vector<Vector<Box> >& a_vectBoxes);

  /// set up grids, etc for AMR run
  void setupForAMRRun();

  /// set up for restart
  void setupForRestart(const string& a_restartFile);

  /// advance solution.  returns time of final solution
  Real run(Real a_maxTime, int a_maxStep);

  ///for convergence tests
  Vector<DisjointBoxLayout>& getGrids()
  {
    return m_grids;
  }
  ///for convergence tests
  Vector<EBISLayout>& getEBISLayouts()
  {
    return m_ebisl;
  }

  ///for convergence tests
  Vector<LevelData<EBCellFAB>* >& getVeloNew()
  {
    return m_velo;
  }
  ///for convergence tests
  Vector<LevelData<EBCellFAB>* >& getGphiNew()
  {
    return m_gphi;
  }

  ///for convergence tests
  Vector<LevelData<EBCellFAB>* >& getVeloOld()
  {
    return m_velo;
  }
  ///for convergence tests
  Vector<LevelData<EBCellFAB>* >& getGphiOld()
  {
    return m_gphi;
  }
  ///for convergence tests
  Vector<LevelData<EBCellFAB>* >& getPresOld()
  {
    return m_pres;
  }
  ///for convergence tests
  Vector<LevelData<EBCellFAB>* >& getPresNew()
  {
    return m_pres;
  }
  ///for convergence tests
  Real getDt()
  {
    return m_dt;
  }

  Vector< LayoutData< Vector< Vector<VolIndex> > >* >& getCoveredFaceLo()
  {
    return m_coveredFaceLitLo;
  }

  Vector< LayoutData< Vector< Vector<VolIndex> > >* >& getCoveredFaceHi()
  {
    return m_coveredFaceLitHi;
  }

  Vector< LayoutData< Vector< IntVectSet > >* >& getCoveredSetsLo()
  {
    return m_coveredSetsLitLo;
  }

  Vector< LayoutData< Vector< IntVectSet > >* >& getCoveredSetsHi()
  {
    return m_coveredSetsLitHi;
  }

  Vector<Real>& getDx()
  {
    return m_dx;
  }

  Vector<ProblemDomain> & getDomain()
  {
    return m_domain;
  }

  EBCompositeMACProjector& getMACProj()
  {
    return *m_macProjector;
  }

  EBCompositeCCProjector& getCCProj()
  {
    return *m_ccProjector;
  }

protected:


  /// write final checkpoint and plot files
  void conclude();

  /// accumulate number of points updated
  void pointsUpdated();

  virtual void cacheExtraTerms(const int a_endLevel = 9999999)
  {
  }

  virtual void uncacheExtraTerms(const int a_endLevel = 9999999)
  {
  }

  virtual void interpolateExtraTerms(const int a_startLevel = -1, const int a_endLevel = 9999999)
  {
  }

  virtual void interpolateExtraTermsForRestart(const int a_startLevel = -1, const int a_endLevel = 9999999)
  {
  }

#ifdef CH_USE_HDF5
  ///
  void writeCheckpointFile();

  ///
  virtual void writeExtraDataToCheckpoint(HDF5Handle& a_handleOut, const int ilev)
  {
  }

  ///
  virtual void readExtraDataFromCheckpoint(HDF5Handle& a_handleIn, const int ilev)
  {
  }

  ///
  virtual void readCheckpointFile(const string& a_restartFile);

  ///
  virtual void writePlotFile();

  ///
  virtual void addExtraOutputData(LevelData<EBCellFAB>& a_outputData, const int a_startDstInterv, const int a_ilev)
  {
  }

  ///
  virtual void addExtraDatatoPlotFile(const string& a_filename)
  {
  }

  ///
  virtual int getNumExtraVars()
  {
    return 0;
  }

  ///
  virtual Vector<string> extraNames()
  {
    return 0;
  }

#endif

  virtual void filter( Vector<LevelData<EBCellFAB>* >&   a_veloc);
  virtual void advance();
  virtual void predictor();
  virtual void corrector();
  virtual void predictVelocity();
  virtual void correctVelocity();

  virtual void predictExtraFieldsPreVelocity()
  {
  }

  virtual void predictExtraFieldsPostVelocity()
  {
  }

  virtual void correctExtraFieldsPreVelocity()
  {
  }

  virtual void correctExtraFieldsPostVelocity()
  {
  }

  virtual void computeExtraSourceForPredictor(Vector<LevelData<EBCellFAB>*> & a_source, const int& a_dir)
  {
  }

  virtual void computeExtraSourceForCorrector(Vector<LevelData<EBCellFAB>*> & a_source)
  {
  }

  virtual void storePredictedVelocity(const Vector<LevelData<EBFluxFAB>* > &                      a_macVelocity,
                                      const Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* > & a_coveredVelocityLo,
                                      const Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* > & a_coveredVelocityHi,
                                      const int &                                                 a_velcomp)
  {
  }

  void transverseVelocityPredictor(Vector<LevelData<EBCellFAB>* >& a_uGradS,
                                   Vector<LevelData<EBCellFAB>* >& a_scalOld,
                                   bool                            a_reallyVelocity);


  /// things to do after a timestep
  virtual void postTimeStep();

  /// create tags
  void tagCells(Vector<IntVectSet>& a_tags);

  /// create tags
  virtual void tagCellsLevel(IntVectSet& a_tags, int a_level);

  /// regrid
  void regrid();

  ///
  /**
     Perform any post-regridding ops.   This includes the final smoothing and
     projecting the velocity and reinitializing the pressure.
   */
  void postRegrid();

  ///
  /**
     Perform any pre-regrdding ops.   This is vel = (I-mu*lapl)vel.
   */
  void preRegrid();

  ///
  void initialGrid(const Vector<DisjointBoxLayout>& a_new_grids);

  /// initialize data
  void initialData();

  /// initialize extra data
  virtual void initialExtraData()
  {
  }

  /// things do do after initialization
  void postInitialize();

  virtual void defineGrids(const Vector<Vector<Box> >& a_vectBoxes);
  // virtual void defineEBISLs(const int a_endLevel = 9999999);
  // virtual void defineExtraEBISLs(const int a_endLevel = 9999999){}
  virtual void defineEBISLs();
  virtual void defineExtraEBISLs()
  {
  }
  void defineNewVel(const int a_startLevel = -1);
  void defineOldVel();
  void definePressure(const int a_startLevel = -1);
  void defineProjections();
  virtual void defineExtraSolvers()
  {
  }

  ///
  virtual void defineExtraTerms(const int a_startLevel = -1)
  {
  }

  virtual void allocateDataHolders();
  /// compute dt
  virtual Real computeDt();

  /// compute dt with initial data
  Real computeInitialDt();

  // compute vorticity on level a_level
  void computeVorticity(LevelData<EBCellFAB>& a_vorticity, int a_level) ;

  /// set up the AMR hierarchy, etc for a_vectGrids.
  void setup();

  /// regrid  once new boxes have been generated
  void regrid(const Vector<Vector<Box> >& a_new_grids);

  void initialGrid(const Vector<Vector<Box> >& a_vectBoxes);

public:
  void
  extrapolateScalarCol(Vector<LevelData<EBFluxFAB>* >                     &  a_macScalar,
                       Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacLo,
                       Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacHi,
                       const RefCountedPtr<EBPhysIBCFactory>              &  a_advectBC,
                       const Vector<LevelData<EBFluxFAB>* >               &  a_advectiveVel,
                       const Vector<LevelData<EBCellFAB>*>*                  a_sourceTerm,
                       const Vector<LevelData<EBCellFAB>* >               &  a_cellScalar,
                       const Vector<LevelData<EBCellFAB>* >               &  a_cellVelocity,
                       Vector<RefCountedPtr<EBLevelAdvect> >              *  a_ebLevAd = NULL);

  //fills temporary stuff.  need to do this for some convergence tests.  do not do this unless you are sure
  //you know the context
  void allocateTemporaries();
  virtual void allocateExtraTemporaries()
  {
  }
  //delete temporary stuff.  need to do this for some convergence tests.  do not do this unless you are sure
  //you know the context
  void deleteTemporaries();
  virtual void deleteExtraTemporaries()
  {
  }
protected:
  void copyNewToOld(int a_ilev);
  void copyOldToNew(int a_ilev);

  void applyEBAMROp(Vector<LevelData<EBCellFAB>* >&       a_lap,
                    Vector<LevelData<EBCellFAB>* >&       a_phi,
                    Vector<LevelData<EBCellFAB>* >&       a_zero,
                    int                                   a_velComp,
                    Real                                  a_time);

  void
  viscousSourceForAdvect(Vector<LevelData<EBCellFAB>* >&       a_source,
                         Vector<LevelData<EBCellFAB>* >&       a_vel,
                         Vector<LevelData<EBCellFAB>* >&       a_zero,
                         int                                   a_velComp,
                         Real                                  a_time);
  void
  extrapolateToCoveredFaces(Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacLo,
                            Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacHi,
                            Vector<LevelData<EBFluxFAB>* >&                       a_macOpen,
                            Vector<LevelData<EBCellFAB>* >&                       a_cellOpen,
                            int                                                   a_idir,
                            Vector<RefCountedPtr<EBLevelAdvect> >              *  a_ebLevAd = NULL);

  void averageDown(Vector<LevelData<EBCellFAB>* >&  a_data);
  void averageDown(Vector<LevelData<EBFluxFAB>* >&  a_data);
  void setCoveredStuffToZero(LevelData<EBCellFAB>& a_vort);
  virtual void averageDownExtraTerms()
  {
  }

  //parmaters
  AMRParameters          m_params;
  EBIBC*                 m_ibc;
  const EBIndexSpace*    m_ebisPtr;

  //grid description
  Vector<ProblemDomain>     m_domain;
  Vector<Real>              m_dx;
  Vector<DisjointBoxLayout> m_grids;
  Vector<EBISLayout>        m_ebisl;
  Vector<EBLevelGrid>       m_eblg;

  EBCompositeCCProjector*        m_ccProjector;
  EBCompositeMACProjector*       m_macProjector;
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >  m_solver[SpaceDim];
  BiCGStabSolver<LevelData<EBCellFAB> >          m_bottomSolver;
  EBSimpleSolver                                 m_bottomSolverSimp;
  RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > >  m_tgaSolver[SpaceDim];
  RefCountedPtr<EBBackwardEuler>                 m_backwardSolver[SpaceDim];
  Vector<RefCountedPtr<EBQuadCFInterp> >    m_quadCFI;
  Vector<RefCountedPtr<EBCoarseAverage> >   m_aveOper;
  Vector<RefCountedPtr<EBCoarseAverage> >   m_aveSpac;

  //  finest extant level
  int m_finestLevel;
  int m_curStep;

  //  time stepping information
  Real m_time;
  Real m_dt;
  Real m_prescribedDt;
  Real m_viscosity;

  long long m_pointsUpdated;

  //switches that i would rather lose but cannot
  bool m_useFixedDt;

  bool m_isSetup;
  //if viscosity == 0 this turns off viscous solvers, etc
  bool m_viscousCalc;

  bool m_doRestart;

  bool m_advanceGphiOnly;
  bool m_steadyState;
  bool m_stopAdvance;

  //data held in memory
  Vector<LevelData<EBCellFAB>* > m_velo;
  Vector<LevelData<EBCellFAB>* > m_pres;
  Vector<LevelData<EBCellFAB>* > m_gphi;

  //making these data part of memory for decoupled calculations (i.e., steady-state velocity) that use advection in EBAMRNoSubcycle
  Vector<LevelData<EBFluxFAB>* >                       m_advVel;
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  m_coveredAdvVelLo;
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  m_coveredAdvVelHi;

  //The following is the ugly part of this stuff.
  //Covered data (advection velocities, etc) are defined with this stuff.
  //These are used to define covered data and to iterate over it.
  //The actual data is defined on the fly to avoid confusion (no really)
  //and to save on memory since all these things are transitory.
  Vector< LayoutData< Vector< Vector<VolIndex> > >* >  m_coveredFaceLitLo;
  Vector< LayoutData< Vector< Vector<VolIndex> > >* >  m_coveredFaceLitHi;
  Vector< LayoutData< Vector< IntVectSet > >* >        m_coveredSetsLitLo;
  Vector< LayoutData< Vector< IntVectSet > >* >        m_coveredSetsLitHi;
  //defines the above stuff
  void defineIrregularData();

  //temporaries that only live during advance
  Vector<LevelData<EBCellFAB>* >                       m_uDotDelU;
  Vector<LevelData<EBCellFAB>* >                       m_cellScratch;

  Vector<LevelData<EBCellFAB>* >                       m_cellScratc1;
  Vector<LevelData<EBCellFAB>* >                       m_cellScratc2;

  Vector<LevelData<EBFluxFAB>* >                       m_macGradient;
  Vector<LevelData<EBFluxFAB>* >                       m_macScratch1;
  Vector<LevelData<EBFluxFAB>* >                       m_macScratch2;

  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  m_coveredScratchLo;
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  m_coveredScratchHi;

  Vector<RefCountedPtr<EBLevelAdvect>  > m_ebLevAd;
  Vector<RefCountedPtr<EBFastFR     >  > m_fluxReg;

  void refluxUDotDelU(Vector<LevelData<EBCellFAB>* >  &  a_uDotDelU,
                      Vector<LevelData<EBFluxFAB>* >  &  a_macAdvVel,
                      Vector<LevelData<EBFluxFAB>* >  &  a_macScalar);
private:

  //we disallow copy construction and assignment of objects
  //which contain pointered data
  EBAMRNoSubcycle(const EBAMRNoSubcycle& a_input)
  {
    MayDay::Error("invalid constructor");
  }

  //we disallow copy construction and assignment of objects
  //which contain pointered data
  void operator=(const EBAMRNoSubcycle& a_input)
  {
    MayDay::Error("invalid constructor");
  }

};

#include "NamespaceFooter.H"
#endif
