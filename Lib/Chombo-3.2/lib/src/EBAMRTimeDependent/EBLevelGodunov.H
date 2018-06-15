#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtg 1-23-2002

#ifndef _EBLEVELGODUNOV_H_
#define _EBLEVELGODUNOV_H_

#include "EBCellFAB.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "PiecewiseLinearFillPatch.H"
#include "EBPWLFillPatch.H"
#include "EBFluxRegister.H"
#include "ProblemDomain.H"
#include "EBPatchGodunovFactory.H"
#include "EBLevelRedist.H"
#include "NamespaceHeader.H"

/// Level Godunov
/**
 */
class EBLevelGodunov
{
public:
  /// Default constructor
  /**
     Object requires define(..) to be called before all other functions.
  */
  EBLevelGodunov();

  /// Destructor
  /**
     destroys all objects created by define(..). Passed in
     data references of define(..) are left alone.
  */
  ~EBLevelGodunov();

  /// Actual constructor.
  /**
     For the coarsest level, an empty DisjointBoxLayout
     is passed in for coaserDisjointBoxLayout.
     Inside the routine, we cast away const-ness on the data members
     for the assignment.  The arguments passed in are maintained const.
     (coding standards).  a_nRefine is the refinement ratio between
     this level and the next coarser level.  a_numGhosts is the number
     of ghost cells in each direction.
  */
  void define(const DisjointBoxLayout&           a_thisDBL,
              const DisjointBoxLayout&           a_coarDBL,
              const EBISLayout&                  a_thisEBISL,
              const EBISLayout&                  a_coarEBISL,
              const ProblemDomain&               a_DProblem,
              const int&                         a_nRefine,
              const RealVect&                    a_dx,
              const bool&                        a_useMassRedist,
              const bool&                        a_doSmushing,
              const bool&                        a_doRZCoords,
              const bool&                        a_hasSourceTerm,
              const EBPatchGodunovFactory* const a_patchGodunov,
              const bool&                        a_hasCoarser,
              const bool&                        a_hasFiner);

  ///
  void floorConserved(LevelData<EBCellFAB>&         a_consState,
                      Real a_time, Real a_dt);

  /// Take one timestep for this grid level.
  /**
     For the coarsest level, UcoarseOld and UCoarseNew are empty
     LevelData<EBCellFAB> objects.  Also, either UCoarseOld or
     UCoarseNew might be empty to indicate that t(nf) = t(nc) the
     one grid is at the current time and interpolation is not required
     for boundary condition generation.  U must be defined on the same
     DisjointBoxLayouts as were used in define(...)
     Coarse flux register is flux register with the coarse level.
     Fine flux register is the flux register with the fine level.
     To the finer level FR, this level is the coarse level.
     To the coarser level FR, this level is the fine level.
  */
  Real step(LevelData<EBCellFAB>&        a_consState,
            LevelData<BaseIVFAB<Real> >& a_massDiff,
            EBFluxRegister&              a_fineFluxRegister,
            EBFluxRegister&              a_coarFluxRegister,
            const LevelData<EBCellFAB>&  a_consStateCoarseOld,
            const LevelData<EBCellFAB>&  a_consStateCoarseNew,
            const Real&                  a_time,
            const Real&                  a_TCold,
            const Real&                  a_TCNew,
            const Real&                  a_dt);

  void
  getDrhoDtOverRho(LevelData<EBCellFAB>& a_drhoDt,
                   const LevelData<EBCellFAB>& a_rhoNew,
                   const LevelData<EBCellFAB>& a_rhoOld,
                   const Real& a_dt);

  /// Get maximum wave speed
  Real getMaxWaveSpeed(const LevelData<EBCellFAB>& a_state);

  bool isDefined() const;


protected:
  void fillConsState(LevelData<EBCellFAB>&         a_consState,
                     const LevelData<EBCellFAB>&   a_consStateCoarseOld,
                     const LevelData<EBCellFAB>&   a_consStateCoarseNew,
                     const Real&                   a_time,
                     const Real&                   a_coarTimeOld,
                     const Real&                   a_coarTimeNew);

  void computeFlattening(Real a_time, Real a_dt,
                         LevelData<EBCellFAB>&         a_consState);

  void doRegularUpdate(EBFluxRegister&               a_fineFluxRegister,
                       EBFluxRegister&               a_coarFluxRegister,
                       Real a_time, Real a_dt,
                       LevelData<EBCellFAB>&         a_consState);

  Real doIrregularUpdate(EBFluxRegister&               a_fineFluxRegister,
                         EBFluxRegister&               a_coarFluxRegister,
                         LevelData<BaseIVFAB<Real> >&  a_massDiff,
                         Real a_time, Real a_dt,
                         LevelData<EBCellFAB>&         a_consState);


  //these are not grown by one.
  LayoutData<IntVectSet> m_irregSetsSmall;

  //these are grown by one in the directions != idir
  LayoutData<IntVectSet> m_irregSetsGrown[SpaceDim];

  LevelData<BaseIFFAB<Real> > m_fluxInterpolants[SpaceDim];
  LevelData<BaseIVFAB<Real> > m_nonConsDivergence;
  LevelData<BaseIVFAB<Real> > m_ebIrregFaceFlux;
  LayoutData<IntVectSet>      m_cfIVS;
  bool              m_hasCoarser;
  bool              m_hasFiner;
  DisjointBoxLayout m_thisGrids;
  DisjointBoxLayout m_coarGrids;
  EBISLayout        m_thisEBISL;
  EBISLayout        m_coarEBISL;
  bool              m_isDefined;
  EBPatchGodunov*   m_ebPatchGodunov;
  RealVect          m_dx;
  ProblemDomain     m_domain;
  int               m_refRatCrse;
  int               m_nCons;
  int               m_nFlux;
  int               m_nGhost;
  EBPWLFillPatch    m_patcher;
  bool              m_useMassRedist;
  EBLevelRedist     m_ebLevelRedist;
  static int        s_timestep;
  bool m_doSmushing;
  bool m_doRZCoords;
  bool m_hasSourceTerm;
  LevelData<EBCellFAB> m_flattening;
private:
  //disallowed for all the usual reasons
  void operator=(const EBLevelGodunov& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBLevelGodunov(const EBLevelGodunov& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
