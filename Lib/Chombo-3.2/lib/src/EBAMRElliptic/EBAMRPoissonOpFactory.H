#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRPOISSONOPFACTORY_H_
#define _EBAMRPOISSONOPFACTORY_H_

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
   Factory class to generate EBAMRPoissonOps.  This follows the AMRLevelOpFactory interface.
*/
class EBAMRPoissonOpFactory: public AMRLevelOpFactory<LevelData<EBCellFAB> >
{
public:
  ///
  virtual ~EBAMRPoissonOpFactory();

  ///
  /**
     a_eblgs : layouts at each AMR level \\
     a_domainFactory : domain boundary conditions \\
     a_ebBCFactory:    eb boundary conditions \\
     a_dxCoarse:      grid spacing at coarsest level \\
     a_origin:        offset to lowest corner of the domain \\
     a_refRatio:     refinement ratios. refRatio[i] is between levels i and i+1 \\
     a_preCondIters:  number of iterations to do for pre-conditioning \\
     a_relaxType:     0 means point Jacobi, 1 is Gauss-Seidel, 2 is line solver. \\
     a_alpha:         coefficent of identity \\
     a_beta:          coefficient of laplacian.\\
     a_time:          time for boundary conditions \\
     a_ghostCellsPhi:  Number of ghost cells in phi, correction (typically one)\\
     a_ghostCellsRhs:  Number of ghost cells in RHS, residual, lphi (typically zero)\\
     Ghost cell arguments are there for caching reasons.  Once you set them, an error is thrown if
     you send in data that does not match.   Use numlevels = -1 if you want to use the
     size of the vectors for numlevels.
  */
  EBAMRPoissonOpFactory(const Vector<EBLevelGrid>&                    a_eblgs,
                        const Vector<int>&                            a_refRatio,
                        const Vector<RefCountedPtr<EBQuadCFInterp> >& a_quadCFI,
                        const RealVect&                               a_dxCoarse,
                        const RealVect&                               a_origin,
                        const int&                                    a_numPreCondIters,
                        const int&                                    a_relaxType,
                        RefCountedPtr<BaseDomainBCFactory>            a_domainBCFactory,
                        RefCountedPtr<BaseEBBCFactory>                a_ebBcFactory,
                        const Real&                                   a_alpha,
                        const Real&                                   a_beta,
                        const Real&                                   a_time,
                        const IntVect&                                a_ghostCellsPhi,
                        const IntVect&                                a_ghostCellsRhs,
                        int a_numLevels = -1);

  ///
  virtual void setData(Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_data)
  {
    m_data = a_data;
    m_dataBased = true;
  }

  ///
  virtual EBAMRPoissonOp*
  MGnewOp(const ProblemDomain& a_FineindexSpace,
          int                  a_depth,
          bool                 a_homoOnly = true);

  EBAMRPoissonOp* createOperator(const EBLevelGrid&             a_eblgMGLevel,
                                 const EBLevelGrid&             a_eblgCoarMG,
                                 const bool&                    a_hasMGObjects,
                                 const bool&                    a_layoutChanged,
                                 const RealVect&                a_dxMGLevel,
                                 const RealVect&                a_dxCoar,
                                 RefCountedPtr<EBQuadCFInterp>& a_quadCFIMGLevel,
                                 const int&                     a_whichLevel,
                                 bool a_amrop);
  ///
  virtual void reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim);

  ///
  virtual EBAMRPoissonOp*
  AMRnewOp(const ProblemDomain& a_FineindexSpace);

  ///
  virtual void AMRreclaim(EBAMRPoissonOp* a_reclaim);

  ///
  /** Refinement ratio between this level and coarser level.
      Returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToFiner(const ProblemDomain& a_domain) const;

  ///
  /**
     testRef is the size of the coarsest domain allowed in multigrid. If testRef=2,
     then the coarsest domain in multigrid will be 2x2(x2)
   **/
  static void setTestRef(int a_testRef)
  {
    s_testRef = a_testRef;
  }

  static void setMaxBoxSize(int a_maxBoxSize)
  {
    s_maxBoxSize = a_maxBoxSize;
  }

  static void setWhichReflux(int & a_whichReflux);

  static int getWhichReflux();

protected:
  static int s_testRef;
  static int s_maxBoxSize;

  bool     m_isDefined;
  int      m_numPreCondIters;
  int      m_relaxType;
  int      m_numLevels;


  Vector<EBLevelGrid>     m_eblgVec;
  Vector<RefCountedPtr<EBQuadCFInterp> > m_quadCFIVec;

  Vector< Vector<EBLevelGrid> >       m_eblgVecMG;
  std::vector< bool  >                m_hasMGObjects;
  std::vector< bool  >                m_layoutChanged;
  std::vector< std::vector<bool> >    m_layoutChangedMG;

  Vector<int>               m_refRatioVec;
  Vector<RealVect>          m_dxVec;
  RealVect                  m_origin;
  Real                      m_alpha;
  Real                      m_beta;
  Real                      m_time;
  const IntVect             m_ghostCellsPhi;
  const IntVect             m_ghostCellsRHS;

  RefCountedPtr<BaseDomainBCFactory>   m_domainBCFactory;
  RefCountedPtr<BaseEBBCFactory>       m_ebBCFactory;

  static int s_whichReflux;

private:

  Vector< RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_data;
  bool m_dataBased;

  ///weak construction bad
  EBAMRPoissonOpFactory()
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  EBAMRPoissonOpFactory(const EBAMRPoissonOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const EBAMRPoissonOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
