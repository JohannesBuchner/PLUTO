#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBCONDUCTIVITYOPFACTORY_H_
#define _EBCONDUCTIVITYOPFACTORY_H_

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
#include "EBConductivityOp.H"
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


//! \class EBConductivityOpFactory
//! Factory class to generate EBConductivityOps.  This follows the
//! AMRLevelOpFactory interface.
class EBConductivityOpFactory: public AMRLevelOpFactory<LevelData<EBCellFAB> >
{
public:

  //! Constructs a factory that builds EBConductivityOps with time-independent
  //! A and B coefficients.
  EBConductivityOpFactory(const Vector<EBLevelGrid>&                                  a_eblgs,
                          const Vector<RefCountedPtr<EBQuadCFInterp> >&               a_quadCFI,
                          const Real&                                                 a_alpha,
                          const Real                                         &        a_beta,
                          const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef,
                          const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_bcoef,
                          const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_bcoefIrreg,
                          const Real&                                                 a_dxCoarse,
                          const Vector<int>&                                          a_refRatio,
                          const RefCountedPtr<BaseDomainBCFactory>&                   a_domainBCFactory,
                          const RefCountedPtr<BaseEBBCFactory>    &                   a_ebBcFactory,
                          const IntVect&                                              a_ghostCellsPhi,
                          const IntVect&                                              a_ghostCellsRhs,
                          const int&                                                  a_relaxType,
                          int a_numLevels = -1);

  //! Constructs a factory that builds EBConductivityOps with a time-dependent
  //! A coefficient and a time-independent B coefficient. The A coefficient
  //! is interpolated between its beginning-of-step and end-of-step values.
  EBConductivityOpFactory(const Vector<EBLevelGrid>&                                  a_eblgs,
                          const Vector<RefCountedPtr<EBQuadCFInterp> >&               a_quadCFI,
                          const Real&                                                 a_alpha,
                          const Real                                         &        a_beta,
                          const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef0,
                          const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef1,
                          const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_bcoef,
                          const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_bcoefIrreg,
                          const Real&                                                 a_dxCoarse,
                          const Vector<int>&                                          a_refRatio,
                          const RefCountedPtr<BaseDomainBCFactory>&                   a_domainBCFactory,
                          const RefCountedPtr<BaseEBBCFactory>    &                   a_ebBcFactory,
                          const IntVect&                                              a_ghostCellsPhi,
                          const IntVect&                                              a_ghostCellsRhs,
                          const int&                                                  a_relaxType,
                          int a_numLevels = -1);

  //! Destructor.
  virtual ~EBConductivityOpFactory();

  ///
  virtual void setData( Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_data)
  {
    m_data = a_data;
    m_dataBased = true;
  }

  ///
  virtual EBConductivityOp*
  MGnewOp(const ProblemDomain& a_FineindexSpace,
          int                  a_depth,
          bool                 a_homoOnly = true);

  EBConductivityOp* createOperator(const EBLevelGrid&             a_eblgMGLevel,
                                   const EBLevelGrid&             a_eblgCoarMG,
                                   const bool&                    a_hasMGObjects,
                                   const RealVect&                a_dxMGLevel,
                                   const RealVect&                a_dxCoar,
                                   const int&                     a_whichLevel);
  ///
  virtual void reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim);

  ///
  virtual EBConductivityOp*
  AMRnewOp(const ProblemDomain& a_FineindexSpace);

  ///
  virtual void AMRreclaim(EBConductivityOp* a_reclaim);

  ///
  /** Refinement ratio between this level and coarser level.
      Returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToFiner(const ProblemDomain& a_domain) const;


  //! Reset the coefficients for the conductivity operator (time-independent).
  void resetCoefficients(const Real&                                                 a_alpha,
                         const Real                                         &        a_beta,
                         const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef,
                         const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_bcoef,
                         const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_bcoefIrreg)
  {
    m_alpha       = a_alpha      ;
    m_beta        = a_beta       ;
    m_acoef       = a_acoef      ;
    m_bcoef       = a_bcoef      ;
    m_bcoefIrreg  = a_bcoefIrreg ;
  }

  //! Reset the coefficients for the conductivity operator (time-dependent A).
  void resetCoefficients(const Real&                                                 a_alpha,
                         const Real                                         &        a_beta,
                         const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef0,
                         const Vector<RefCountedPtr<LevelData<EBCellFAB> > >&        a_acoef1,
                         const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&        a_bcoef,
                         const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_bcoefIrreg)
  {
    m_alpha       = a_alpha      ;
    m_beta        = a_beta       ;
    m_acoef0      = a_acoef0     ;
    m_acoef1      = a_acoef1     ;
    m_bcoef       = a_bcoef      ;
    m_bcoefIrreg  = a_bcoefIrreg ;
  }

protected:

  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_data;
  bool m_dataBased;
  Vector< Vector<EBLevelGrid> >                                   m_eblgsMG;

  //! Time-independent A multigrid coefficients
  Vector< Vector< RefCountedPtr<LevelData<EBCellFAB> > > >        m_acoefMG;

  //! Beginning-of-step (time-dependent) A multigrid coefficients
  Vector< Vector< RefCountedPtr<LevelData<EBCellFAB> > > >        m_acoefMG0;

  //! End-of-step (time-dependent) A multigrid coefficients
  Vector< Vector< RefCountedPtr<LevelData<EBCellFAB> > > >        m_acoefMG1;

  Vector< Vector< RefCountedPtr<LevelData<EBFluxFAB> > > >        m_bcoefMG;
  Vector< Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > > > m_bcoefIrregMG;

  int                                                  m_relaxType;
  Vector<EBLevelGrid>                                  m_eblgs;
  Vector<RefCountedPtr<EBQuadCFInterp> >               m_quadCFI;
  Real                                                 m_alpha;
  Real                                                 m_beta;

  //! Time-independent A coefficients.
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >        m_acoef;

  //! Beginning-of-step (time-dependent) A coefficients.
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >        m_acoef0;

  //! End-of-step (time-dependent) A coefficient.
  Vector<RefCountedPtr<LevelData<EBCellFAB> > >        m_acoef1;

  Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        m_bcoef;
  Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_bcoefIrreg;
  Real                                                 m_dxCoarse;
  Vector<int>                                          m_refRatio;
  RefCountedPtr<BaseDomainBCFactory>                   m_domainBCFactory;
  RefCountedPtr<BaseEBBCFactory>                       m_ebBCFactory;
  IntVect                                              m_ghostCellsPhi;
  IntVect                                              m_ghostCellsRhs;
  int                                                  m_numLevels;
  std::vector< bool  >                                 m_hasMGObjects;
  Vector<Real>                                         m_dx;
private:
  ///weak construction bad
  EBConductivityOpFactory()
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  EBConductivityOpFactory(const EBConductivityOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  void operator=(const EBConductivityOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
