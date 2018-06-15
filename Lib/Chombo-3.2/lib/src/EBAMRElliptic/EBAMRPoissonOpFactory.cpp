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

#include "EBAMRPoissonOp.H"
#include "EBArith.H"

#include "CH_Timer.H"
#include "EBAMRPoissonOpFactory.H"
#include "NamespaceHeader.H"

int EBAMRPoissonOpFactory::s_testRef = 2;
int EBAMRPoissonOpFactory::s_maxBoxSize = 32;
int EBAMRPoissonOpFactory::s_whichReflux = 2;

////
EBAMRPoissonOpFactory::~EBAMRPoissonOpFactory()
{
}
///
int
EBAMRPoissonOpFactory::
refToFiner(const ProblemDomain& a_domain) const
{
  int retval = -1;
  bool found = false;
  for (int ilev = 0; ilev < m_eblgVec.size(); ilev++)
    {
      if (m_eblgVec[ilev].getDomain() == a_domain)
        {
          retval = m_refRatioVec[ilev];
          found = true;
        }
    }
  if (!found)
    {
      MayDay::Error("Domain not found in AMR hierarchy");
    }
  return retval;
}
///
EBAMRPoissonOpFactory::
EBAMRPoissonOpFactory(const Vector<EBLevelGrid>&                    a_eblgVec,
                      const Vector<int>&                            a_refRatio,
                      const Vector<RefCountedPtr<EBQuadCFInterp> >& a_quadCFI,
                      const RealVect&                               a_dxCoar,
                      const RealVect&                               a_origin,
                      const int&                                    a_numPreCondIters,
                      const int&                                    a_relaxType,
                      RefCountedPtr<BaseDomainBCFactory>            a_domainBCFactory,
                      RefCountedPtr<BaseEBBCFactory>                a_ebBCFactory,
                      const Real&                                   a_alpha,
                      const Real&                                   a_beta,
                      const Real&                                   a_time,
                      const IntVect&                                a_ghostCellsPhi,
                      const IntVect&                                a_ghostCellsRHS,
                      int a_numLevels)
  : m_ghostCellsPhi( a_ghostCellsPhi ),
    m_ghostCellsRHS( a_ghostCellsRHS )
{
  CH_assert(a_eblgVec.size() <= a_refRatio.size());
  m_dataBased = false;
  if (a_numLevels > 0)
    {
      m_numLevels = a_numLevels;
    }
  else
    {
      m_numLevels = a_eblgVec.size();
    }

  m_domainBCFactory = a_domainBCFactory;
  m_ebBCFactory     = a_ebBCFactory;

  m_eblgVec         = a_eblgVec;
  m_quadCFIVec      = a_quadCFI;
  m_refRatioVec     = a_refRatio;
  m_numPreCondIters = a_numPreCondIters;
  m_relaxType       = a_relaxType;
  m_origin          = a_origin;
  m_dxVec.resize(m_numLevels);

  m_dxVec[0] = a_dxCoar;
  for (int ilev = 1; ilev < m_numLevels; ilev++)
    {
      m_dxVec[ilev] = m_dxVec[ilev-1];
      m_dxVec[ilev] /= m_refRatioVec[ilev-1];
    }
  m_alpha = a_alpha;
  m_beta = a_beta;
  m_time = a_time;

  s_whichReflux = 2;

  m_eblgVecMG.resize(m_numLevels);
  m_hasMGObjects.resize(m_numLevels);
  m_layoutChanged.resize(m_numLevels);
  m_layoutChangedMG.resize(m_numLevels);
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      if ((ilev==0) || (m_refRatioVec[ilev] > 2))
        {
          m_hasMGObjects[ilev] = (s_testRef<s_maxBoxSize);

          int mgRef = 2;
          m_eblgVecMG[ilev].resize(0);
          m_eblgVecMG[ilev].push_back(m_eblgVec[ilev]);
          m_layoutChangedMG[ilev].resize(0);
          m_layoutChangedMG[ilev].push_back(m_layoutChanged[ilev]);

          bool hasCoarser = true;
          hasCoarser = (s_testRef<s_maxBoxSize);
          bool atAMRLevel = true;
          while (hasCoarser)
            {
              int imgsize = m_eblgVecMG[ilev].size();
              const EBLevelGrid& eblgFine=  m_eblgVecMG[ilev][imgsize-1];
              DisjointBoxLayout dblCoarMG;
              ProblemDomain  domainCoarMG;
              bool layoutChanged;
              hasCoarser = EBAMRPoissonOp::getCoarserLayouts(dblCoarMG,
                                                             domainCoarMG,
                                                             eblgFine.getDBL(),
                                                             eblgFine.getEBISL(),
                                                             eblgFine.getDomain(),
                                                             mgRef,
                                                             eblgFine.getEBIS(),
                                                             s_maxBoxSize,
                                                             layoutChanged,
                                                             s_testRef);


              if (atAMRLevel)
                {
                  m_layoutChanged[ilev] = layoutChanged;
                  atAMRLevel= false;
                }

              if (hasCoarser)
                {
                  int numGhost = 4;
                  for (int jj = 0; jj< SpaceDim; jj++)
                    numGhost = Max(numGhost, a_ghostCellsPhi[jj]);
                  m_eblgVecMG[ilev].push_back(EBLevelGrid(dblCoarMG, domainCoarMG,  numGhost, eblgFine.getEBIS()));
                  m_layoutChangedMG[ilev].push_back(layoutChanged);
                }
            }

        }
      else
        {
          m_hasMGObjects[ilev] = false;
        }
    }
}

EBAMRPoissonOp*
EBAMRPoissonOpFactory::
MGnewOp(const ProblemDomain& a_domainFine,
        int                  a_depth,
        bool                 a_homoOnly)
{

  //find out if there is a real starting point here.
  int whichlev;
  bool found = false;

  for (int ilev = 0; ilev < m_numLevels && !found; ilev++)
    {
      if (a_domainFine == m_eblgVec[ilev].getDomain())
        {
          found = true;
          whichlev = ilev ;
        }
    }
  if (!found)
    {
      MayDay::Error("No corresponding AMRLevel to starting point of MGnewOp");
    }

  //multigrid operator.  coarn by two from depth  no
  //coarse or finer to worry about.
  EBLevelGrid eblgMGLevel;
  EBLevelGrid eblgCoarMG;
  RealVect      dxMGLevel;
  RefCountedPtr<EBQuadCFInterp> quadCFIMGLevel; //only defined if on an amr level
  RealVect dxCoar = RealVect::Unit;
  dxCoar *= -1.0;
  if (whichlev > 0)
    {
      dxCoar = m_dxVec[whichlev-1];
    }
  bool hasCoarMGObjects = false;
  bool coarLayoutChanged = true;
  if (a_depth == 0)
    {
      eblgMGLevel    = m_eblgVec[whichlev];
      dxMGLevel      = m_dxVec[whichlev];
      quadCFIMGLevel = m_quadCFIVec[whichlev];

      hasCoarMGObjects = m_hasMGObjects[whichlev];
      if (hasCoarMGObjects)
        {
          eblgCoarMG = m_eblgVecMG[whichlev][1];
          coarLayoutChanged = m_layoutChangedMG[whichlev][1];
        }
    }
  else
    {
      int icoar = 1;
      for (int idep = 0; idep < a_depth; idep++)
        {
          icoar *= 2;
        }
      const ProblemDomain domainFine = m_eblgVec[whichlev].getDomain();
      ProblemDomain domainBoxMGLevel = coarsen(domainFine, icoar);
      bool foundMGLevel = false;
      int numMGLevels = m_eblgVecMG[whichlev].size();
      for (int img = 0; img < numMGLevels; img++)
        {
          if (m_eblgVecMG[whichlev][img].getDomain() == domainBoxMGLevel)
            {
              eblgMGLevel = m_eblgVecMG[whichlev][img];
              foundMGLevel = true;

              hasCoarMGObjects = ((img+1) < (numMGLevels));
              if (hasCoarMGObjects)
                {
                  eblgCoarMG = m_eblgVecMG[whichlev][img+1];
                  coarLayoutChanged = m_layoutChangedMG[whichlev][img+1];
                }
              break;
            }
        }
      bool coarsenable = foundMGLevel;

      dxMGLevel = m_dxVec[whichlev];
      dxMGLevel *= Real(icoar);

      if (!coarsenable)
        {
          //not coarsenable.
          //return null
          return NULL;
        }
    }

  //creates coarse and finer info and bcs and all that
  EBAMRPoissonOp* op = createOperator(eblgMGLevel, eblgCoarMG, hasCoarMGObjects, coarLayoutChanged,
                                      dxMGLevel, dxCoar, quadCFIMGLevel, whichlev, false);

  return op;
}
/////
EBAMRPoissonOp*
EBAMRPoissonOpFactory::
AMRnewOp(const ProblemDomain& a_domainFine)
{
  //figure out which level we are at.
  int whichlev;
  bool found = false;

  for (int ilev = 0; ilev < m_numLevels && !found; ilev++)
    {
      if (a_domainFine == m_eblgVec[ilev].getDomain())
        {
          found = true;
          whichlev = ilev ;
        }
    }
  if (!found)
    {
      MayDay::Error("No corresponding AMRLevel to starting point of AMRnewOp");
    }

  RealVect dxCoar = RealVect::Unit;
  dxCoar *= -1.0;
  if (whichlev > 0)
    {
      dxCoar = m_dxVec[whichlev-1];
    }
  //creates coarse and finer info and bcs and all that
  EBLevelGrid      eblgMGLevel = m_eblgVec[whichlev];
  RealVect           dxMGLevel =   m_dxVec[whichlev];
  RefCountedPtr<EBQuadCFInterp> quadCFIMGLevel = m_quadCFIVec[whichlev];

  bool hasCoarMGObjects = m_hasMGObjects[whichlev];
  bool coarLayoutChanged = m_layoutChanged[whichlev];
  EBLevelGrid eblgCoarMG;
  if (hasCoarMGObjects)
    {
      eblgCoarMG = m_eblgVecMG[whichlev][1];
    }

  EBAMRPoissonOp* op = createOperator(eblgMGLevel, eblgCoarMG, hasCoarMGObjects, coarLayoutChanged,
                                      dxMGLevel, dxCoar, quadCFIMGLevel, whichlev, true);

  return op;
}
//////
EBAMRPoissonOp*
EBAMRPoissonOpFactory::createOperator(const EBLevelGrid&             a_eblgMGLevel,
                                      const EBLevelGrid&             a_eblgCoarMG,
                                      const bool&                    a_hasMGObjects,
                                      const bool&                    a_layoutChanged,
                                      const RealVect&                a_dxMGLevel,
                                      const RealVect&                a_dxCoar,
                                      RefCountedPtr<EBQuadCFInterp>& a_quadCFIMGLevel,
                                      const int&                     a_whichLevel,
                                      bool a_amrop)
{

  RefCountedPtr<BaseEBBC> ebbcMGLevel(
    m_ebBCFactory->create(    a_eblgMGLevel.getDomain(), a_eblgMGLevel.getEBISL(), a_dxMGLevel,
                              &m_ghostCellsPhi, &m_ghostCellsRHS));
  RefCountedPtr<BaseDomainBC> dombcMGLevel(
    m_domainBCFactory->create(a_eblgMGLevel.getDomain(), a_eblgMGLevel.getEBISL(), a_dxMGLevel));

  //fine and coarse stuff undefined because this is a MG level
  EBLevelGrid  eblgFine,  eblgCoar;

  RefCountedPtr<BaseEBBC> ebbcCoar;
  RefCountedPtr<BaseEBBC> ebbcFine;
  RefCountedPtr<BaseDomainBC> dombcCoar;
  RefCountedPtr<BaseDomainBC> dombcFine;
  int refToFine   = 1;
  int refToCoar = 1;
  bool hasFine   = ((a_whichLevel < m_numLevels - 1)  && (a_eblgMGLevel.getDomain() == m_eblgVec[a_whichLevel].getDomain()));
  bool hasCoar   = ((a_whichLevel > 0) && (a_eblgMGLevel.getDomain() == m_eblgVec[a_whichLevel].getDomain()));
  if (hasFine)
    {
      refToFine = m_refRatioVec[a_whichLevel];
      eblgFine  = m_eblgVec[a_whichLevel+1];
    }
  if (hasCoar)
    {
      refToCoar  =  m_refRatioVec[a_whichLevel-1];
      eblgCoar   =      m_eblgVec[a_whichLevel-1];
    }
  if (a_amrop && m_dataBased)
    {
      ebbcMGLevel->setData(m_data[a_whichLevel]);
    }
  EBAMRPoissonOp* op = new EBAMRPoissonOp(eblgFine, a_eblgMGLevel, eblgCoar, a_eblgCoarMG,
                                          a_quadCFIMGLevel, dombcMGLevel, ebbcMGLevel,
                                          a_dxMGLevel, a_dxCoar, m_origin,  refToFine, refToCoar,
                                          hasFine, hasCoar, a_hasMGObjects, a_layoutChanged,
                                          m_numPreCondIters,
                                          m_relaxType,
                                          m_alpha, m_beta,
                                          m_ghostCellsPhi, m_ghostCellsRHS, s_testRef);


  return op;
}
/****/
void
EBAMRPoissonOpFactory::
reclaim(MGLevelOp<LevelData<EBCellFAB> >* a_reclaim)
{
  delete a_reclaim;
}
/****/
void
EBAMRPoissonOpFactory::
AMRreclaim(EBAMRPoissonOp* a_reclaim)
{
  delete a_reclaim;
}

void EBAMRPoissonOpFactory::setWhichReflux(int & a_whichReflux)
{
  s_whichReflux = a_whichReflux;
}

int EBAMRPoissonOpFactory::getWhichReflux()
{
  return s_whichReflux;
}

#include "NamespaceFooter.H"
