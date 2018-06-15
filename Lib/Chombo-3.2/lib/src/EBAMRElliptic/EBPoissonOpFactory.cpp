#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LoadBalance.H"
#include "EBArith.H"

#include "EBPoissonOp.H"
#include "EBAMRPoissonOp.H"
#include "EBArith.H"

#include "CH_Timer.H"
#include "EBPoissonOpFactory.H"
#include "NamespaceHeader.H"

////
EBPoissonOpFactory::~EBPoissonOpFactory()
{
}
///
EBPoissonOpFactory::
EBPoissonOpFactory(const EBLevelGrid&                            a_eblg,
                   const RealVect&                               a_dx,
                   const RealVect&                               a_origin,
                   const int&                                    a_orderEB,
                   const int&                                    a_numPreCondIters,
                   const int&                                    a_relaxType,
                   RefCountedPtr<BaseDomainBCFactory>            a_domainBCFactory,
                   RefCountedPtr<BaseEBBCFactory>                a_ebBCFactory,
                   const Real&                                   a_alpha,
                   const Real&                                   a_beta,
                   const IntVect&                                a_ghostCellsPhi,
                   const IntVect&                                a_ghostCellsRHS)
  : m_ghostCellsPhi( a_ghostCellsPhi ),
    m_ghostCellsRHS( a_ghostCellsRHS )
{
  m_domainBCFactory = a_domainBCFactory;
  m_ebBCFactory     = a_ebBCFactory;

  m_eblg            = a_eblg;
  m_orderEB         = a_orderEB;
  m_numPreCondIters = a_numPreCondIters;
  m_relaxType       = a_relaxType;
  m_origin          = a_origin;
  m_dx              = a_dx;

  m_alpha = a_alpha;
  m_beta = a_beta;

  int mgRef = 2;
  m_eblgVecMG.push_back(m_eblg);

  bool hasCoarser = true;
  while (hasCoarser)
    {
      int imgsize = m_eblgVecMG.size();
      const EBLevelGrid& eblgFine=  m_eblgVecMG[imgsize-1];
      DisjointBoxLayout dblCoarMG;
      ProblemDomain  domainCoarMG;
      int maxBoxSize = 32;
      bool dumbool;
      hasCoarser = EBAMRPoissonOp::getCoarserLayouts(dblCoarMG,
                                                     domainCoarMG,
                                                     eblgFine.getDBL(),
                                                     eblgFine.getEBISL(),
                                                     eblgFine.getDomain(),
                                                     mgRef,
                                                     eblgFine.getEBIS(),
                                                     maxBoxSize, dumbool);

      if (hasCoarser)
        {
          m_eblgVecMG.push_back(EBLevelGrid(dblCoarMG, domainCoarMG, 4, eblgFine.getEBIS()));
        }
    }

}

EBPoissonOp*
EBPoissonOpFactory::
MGnewOp(const ProblemDomain& a_domainFine,
        int                  a_depth,
        bool                 a_homoOnly)
{

  //find out if there is a real starting point here.
  if (a_domainFine != m_eblg.getDomain())
    {
      MayDay::Error("No corresponding AMRLevel to starting point of MGnewOp");
    }

  //multigrid operator.  coarn by two from depth  no
  //coarse or finer to worry about.
  EBLevelGrid eblgMGLevel;
  EBLevelGrid eblgCoarMG;
  RealVect      dxMGLevel;
  RefCountedPtr<EBQuadCFInterp> quadCFIMGLevel; //only defined if on an amr level

  bool hasCoarMGObjects = false;
  int icoar = 1;
  for (int idep = 0; idep < a_depth; idep++)
    {
      icoar *= 2;
    }
  const ProblemDomain domainFine = m_eblg.getDomain();
  ProblemDomain domainBoxMGLevel = coarsen(domainFine, icoar);
  bool foundMGLevel = false;
  int numMGLevels = m_eblgVecMG.size();
  for (int img = 0; img < numMGLevels; img++)
    {
      if (m_eblgVecMG[img].getDomain() == domainBoxMGLevel)
        {
          eblgMGLevel = m_eblgVecMG[img];
          foundMGLevel = true;

          hasCoarMGObjects = ((img+1) < (numMGLevels));
          if (hasCoarMGObjects)
            {
              eblgCoarMG = m_eblgVecMG[img+1];
            }
          break;
        }
    }
  bool coarsenable = foundMGLevel;

  dxMGLevel = m_dx;
  dxMGLevel *= Real(icoar);

  if (!coarsenable)
    {
      //not coarsenable.
      //return null
      return NULL;
    }
  //creates coarse and finer info and bcs and all that
  EBPoissonOp* op = createOperator(eblgMGLevel, eblgCoarMG, hasCoarMGObjects, dxMGLevel);

  return op;
}
//////
EBPoissonOp*
EBPoissonOpFactory::createOperator(const EBLevelGrid&             a_eblgMGLevel,
                                   const EBLevelGrid&             a_eblgCoarMG,
                                   const bool&                    a_hasMGObjects,
                                   const RealVect&                a_dxMGLevel)
{

  RefCountedPtr<BaseEBBC> ebbcMGLevel(
    m_ebBCFactory->create(    a_eblgMGLevel.getDomain(), a_eblgMGLevel.getEBISL(), a_dxMGLevel,
                              &m_ghostCellsPhi, &m_ghostCellsRHS));
  RefCountedPtr<BaseDomainBC> dombcMGLevel(
    m_domainBCFactory->create(a_eblgMGLevel.getDomain(), a_eblgMGLevel.getEBISL(), a_dxMGLevel));

  EBPoissonOp* op = new EBPoissonOp(a_eblgMGLevel, a_eblgCoarMG,
                                    dombcMGLevel, ebbcMGLevel,
                                    a_dxMGLevel, m_origin,
                                    a_hasMGObjects,
                                    m_numPreCondIters,
                                    m_relaxType, m_orderEB,
                                    m_alpha, m_beta,
                                    m_ghostCellsPhi, m_ghostCellsRHS);

  return op;
}
#include "NamespaceFooter.H"
