#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPOISSONOPFACTORY_H__
#define _EBPOISSONOPFACTORY_H__

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
#include "EBPoissonOp.H"
#include "EBLevelGrid.H"
#include "NamespaceHeader.H"


///
/**
   Factory class to generate EBPoissonOps.  This follows the AMRLevelOpFactory interface.
*/
class EBPoissonOpFactory: public MGLevelOpFactory<LevelData<EBCellFAB> >
{
public:
  ///
  virtual ~EBPoissonOpFactory();

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
     a_orderEB:       0 to not do flux interpolation at cut faces. \\
     a_alpha:         coefficent of identity \\
     a_beta:          coefficient of laplacian.\\
     a_ghostCellsPhi:  Number of ghost cells in phi, correction (typically one)\\
     a_ghostCellsRhs:  Number of ghost cells in RHS, residual, lphi (typically zero)\\
     Ghost cell arguments are there for caching reasons.  Once you set them, an error is thrown if
     you send in data that does not match.
  */
  EBPoissonOpFactory(const EBLevelGrid&                            a_eblgs,
                     const RealVect&                               a_dx,
                     const RealVect&                               a_origin,
                     const int&                                    a_orderEB,
                     const int&                                    a_numPreCondIters,
                     const int&                                    a_relaxType,
                     RefCountedPtr<BaseDomainBCFactory>            a_domainBCFactory,
                     RefCountedPtr<BaseEBBCFactory>                a_ebBcFactory,
                     const Real&                                   a_alpha,
                     const Real&                                   a_beta,
                     const IntVect&                                a_ghostCellsPhi,
                     const IntVect&                                a_ghostCellsRhs);

  ///
  virtual EBPoissonOp*
  MGnewOp(const ProblemDomain& a_FineindexSpace,
          int                  a_depth,
          bool                 a_homoOnly = true);

  EBPoissonOp* createOperator(const EBLevelGrid&             a_eblgMGLevel,
                              const EBLevelGrid&             a_eblgCoarMG,
                              const bool&                    a_hasMGObjects,
                              const RealVect&                a_dxMGLevel);


protected:
  int      m_orderEB;
  int      m_numPreCondIters;
  int      m_relaxType;


  EBLevelGrid               m_eblg;
  Vector<EBLevelGrid>       m_eblgVecMG;

  RealVect                  m_dx;
  RealVect                  m_origin;
  Real                      m_alpha;
  Real                      m_beta;
  const IntVect             m_ghostCellsPhi;
  const IntVect             m_ghostCellsRHS;

  RefCountedPtr<BaseDomainBCFactory>   m_domainBCFactory;
  RefCountedPtr<BaseEBBCFactory>       m_ebBCFactory;

private:
  ///weak construction bad
  EBPoissonOpFactory()
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  EBPoissonOpFactory(const EBPoissonOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const EBPoissonOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
