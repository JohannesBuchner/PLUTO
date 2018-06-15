#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBVISCOUSTENSOROPFACTORY_H_
#define _EBVISCOUSTENSOROPFACTORY_H_

#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "Vector.H"
#include <map>
#include "RefCountedPtr.H"

#include "AMRMultiGrid.H"

#include "EBIndexSpace.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "EBViscousTensorOp.H"
#include "EBLevelDataOps.H"
#include "BaseEBBC.H"
#include "BaseDomainBC.H"
#include "CFIVS.H"
#include "EBFluxRegister.H"
#include "EBMGAverage.H"
#include "EBMGInterp.H"
#include "EBCoarsen.H"
#include "PolyGeom.H"
#include "EBAMRPoissonOp.H"
#include "EBLevelGrid.H"
#include "NamespaceHeader.H"

///
/**
   Factory class to generate EBViscousTensorOps.  This follows the AMRLevelOpFactory interface.
*/
class EBViscousTensorOpFactory: public AMRLevelOpFactory<LevelData<EBCellFAB> >
{
public:

  ///
  virtual ~EBViscousTensorOpFactory();

  //! Constructs a factory that builds EBViscousTensorOps with time-independent
  //! A and B coefficients.
  EBViscousTensorOpFactory(const Vector<EBLevelGrid>&                                  a_eblgs,
                           const Real&                                                 a_alpha,
                           const Real                                         &        a_beta,
                           const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef,
                           const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_eta,
                           const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_lambda,
                           const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_etaIrreg,
                           const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_lambdaIrreg,
                           const Real&                                                 a_dxCoarse,
                           const Vector<int>&                                          a_refRatio,
                           const RefCountedPtr<BaseDomainBCFactory>&                   a_domainBCFactory,
                           const RefCountedPtr<BaseEBBCFactory>    &                   a_ebBcFactory,
                           const IntVect&                                              a_ghostCellsPhi,
                           const IntVect&                                              a_ghostCellsRhs,
                           int a_numLevels = -1,
                           bool a_noMG = false);

  //! Constructs a factory that builds EBConductivityOps with a time-dependent
  //! A coefficient and a time-independent B coefficient. The A coefficient
  //! is interpolated between its beginning-of-step and end-of-step values.
  EBViscousTensorOpFactory(const Vector<EBLevelGrid>&                                  a_eblgs,
                           const Real&                                                 a_alpha,
                           const Real                                         &        a_beta,
                           const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef0,
                           const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef1,
                           const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_eta,
                           const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_lambda,
                           const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_etaIrreg,
                           const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_lambdaIrreg,
                           const Real&                                                 a_dxCoarse,
                           const Vector<int>&                                          a_refRatio,
                           const RefCountedPtr<BaseDomainBCFactory>&                   a_domainBCFactory,
                           const RefCountedPtr<BaseEBBCFactory>    &                   a_ebBcFactory,
                           const IntVect&                                              a_ghostCellsPhi,
                           const IntVect&                                              a_ghostCellsRhs,
                           int a_numLevels = -1,
                           bool a_noMG = false);

  ///
  virtual EBViscousTensorOp*
  MGnewOp(const ProblemDomain& a_FineindexSpace,
          int                  a_depth,
          bool                 a_homoOnly = true);

  EBViscousTensorOp* createOperator(const EBLevelGrid&             a_eblgMGLevel,
                                    const EBLevelGrid&             a_eblgCoarMG,
                                    const bool&                    a_hasMGObjects,
                                    const RealVect&                a_dxMGLevel,
                                    const RealVect&                a_dxCoar,
                                    const int&                     a_whichLevel);
  ///
  virtual void reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim);

  ///
  virtual EBViscousTensorOp*
  AMRnewOp(const ProblemDomain& a_FineindexSpace);

  ///
  virtual void AMRreclaim(EBViscousTensorOp* a_reclaim);

  ///
  /** Refinement ratio between this level and coarser level.
      Returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToFiner(const ProblemDomain& a_domain) const;


  void resetCoefficients(const Real&                                                 a_alpha,
                         const Real                                         &        a_beta,
                         const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef,
                         const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_eta,
                         const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_lambda,
                         const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_etaIrreg,
                         const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_lambdaIrreg)
  {
    m_alpha       = a_alpha      ;
    m_beta        = a_beta       ;
    m_acoef       = a_acoef      ;
    m_eta         = a_eta        ;
    m_lambda      = a_lambda     ;
    m_etaIrreg    = a_etaIrreg   ;
    m_lambdaIrreg = a_lambdaIrreg;
  }

protected:
  Vector< Vector<EBLevelGrid> >                                  m_eblgsMG;

  //! Time-independent a multigrid coefficients
  Vector< Vector< RefCountedPtr<LevelData<EBCellFAB> > > >       m_acoefMG;

  Vector< Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        > m_etaMG;
  Vector< Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        > m_lambdaMG;
  Vector< Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > > m_etaIrregMG;
  Vector< Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > > m_lambdaIrregMG;
  Vector< Vector<EBViscousTensorOp*> > m_MGops;
  std::vector< bool  >                                 m_hasMGObjects;

  Vector<EBLevelGrid>                                  m_eblgs;
  Real                                                 m_alpha;
  Real                                                 m_beta;

  //! Time-independent a coefficients.
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >        m_acoef;

  //! Beginning-of-step (time-dependent) a coefficients.
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >        m_acoef0;

  //! End-of-step (time-dependent) a coefficient.
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >        m_acoef1;

  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        m_eta;
  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        m_lambda;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_etaIrreg;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_lambdaIrreg;
  Real                                                 m_dxCoarse;
  Vector<int>                                          m_refRatio;
  RefCountedPtr<BaseDomainBCFactory>                   m_domainBCFactory;
  RefCountedPtr<BaseEBBCFactory>                       m_ebBCFactory;
  IntVect                                              m_ghostCellsPhi;
  IntVect                                              m_ghostCellsRhs;
  int                                                  m_numLevels;
  Vector<Real>                                         m_dx;
  bool m_noMG;
private:
  ///weak construction bad
  EBViscousTensorOpFactory()
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  EBViscousTensorOpFactory(const EBViscousTensorOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  void operator=(const EBViscousTensorOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
