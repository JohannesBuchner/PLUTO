#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL, DTG

#include "DebugOut.H"
#include "EBISLayout.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "VoFIterator.H"
#include "EBGraphFactory.H"
#include "EBDataFactory.H"
#include "NamespaceHeader.H"
const EBIndexSpace*
EBISLayout::
getEBIS() const
{
  return m_implem->getEBIS();
}
void
EBISLayout::
setEBIS(const EBIndexSpace* const a_ebisPtr)
{
  m_implem->setEBIS(a_ebisPtr);
}
/****************/
bool EBISLayoutImplem::s_verbose = false;
/****************/
void EBISLayoutImplem::setVerbose(bool a_verbose)
{
  s_verbose = a_verbose;
}
/****************/
void
EBISLayout::define(const ProblemDomain& a_domain,
                   const DisjointBoxLayout& a_grids,
                   const int& a_nghost,
                   const LevelData<EBGraph>& a_graph,
                   const LevelData<EBData> & a_data)
{
  m_nghost = a_nghost;
  m_implem->define(a_domain, a_grids, a_nghost, a_graph, a_data);
}
/****************/
void
EBISLayoutImplem::define(const ProblemDomain& a_domain,
                         const DisjointBoxLayout& a_grids,
                         const int& a_nghost,
                         const LevelData<EBGraph>& a_graph,
                         const LevelData<EBData> & a_data)
{
  CH_TIME("EBISLayoutImplem::define");
#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
  m_domain = a_domain;
  m_nghost = a_nghost;
  m_dblInputDom = a_grids;
  m_fineLevels.resize(0);
  m_coarLevels.resize(0);
  m_maxCoarseningRatio = 2;
  m_maxRefinementRatio = 2;

  //generate the layout using the input layout and the ghost cells and domain;
  m_blGhostDom.deepCopy(a_grids);
  //they have to work with the same iterator
  m_blGhostDom.grow(a_nghost);
  m_blGhostDom &= m_domain.domainBox();
  m_blGhostDom.close();
  if (s_verbose)
    {
      pout() << "in ebislayoutimplem::define" << endl;
      pout() << "input ghost = " << a_nghost  << ", input grids = " << a_grids << endl;
    }

  EBGraphFactory graphFact(m_domain);
  //make graph ghost one bigger so we can define the data
  IntVect ivghostgraph = (a_nghost+1)*IntVect::Unit;
  IntVect ivghostdata  = (a_nghost)*IntVect::Unit;
  Interval interv(0, 0);
  LevelData<EBGraph> localGraph;
  {
    CH_TIME("LocalGraphDefine");
    localGraph.define(a_grids, 1, ivghostgraph, graphFact);
  }
  a_graph.copyTo(interv, localGraph, interv);

  if (s_verbose)
    {
      pout() << "ebislayoutimplem::define just copied graph " << endl;
      BaseIVFAB<VolData>::setVerbose(true);
    }
  EBDataFactory dataFact;
  LevelData<EBData> localData;
  {
    CH_TIME("LocalDataDefine");
    localData.define(a_grids, 1, ivghostdata, dataFact);
  }
  {
    CH_TIME("EBDataCreate");
  for (DataIterator dit= a_grids.dataIterator(); dit.ok(); ++dit)
    {
      Box localBox = grow(a_grids.get(dit()), a_nghost);
      localBox &= m_domain;
      localData[dit()].defineVoFData(localGraph[dit()], localBox);
      localData[dit()].defineFaceData(localGraph[dit()],localBox);
    }
  a_data.copyTo(interv, localData, interv);
  }

  if (s_verbose)
    {
      pout() << "EBISLayoutImplem::define just copied data " << endl;
      BaseIVFAB<VolData>::setVerbose(false);
    }
  {
    CH_TIME("EBISBoxesDefine");
  //define the layouts data with the ghosted layout.  this includes
  //the problem domain stuff
  m_ebisBoxes.define(m_blGhostDom);
  for (DataIterator dit= a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph& graphlocal = localGraph[dit()];
      Box graphregion=  graphlocal.getRegion();
      m_ebisBoxes[dit()].define(localGraph[dit()], localData[dit()]);
    }
  }
  m_defined = true;
}

bool EBISLayout::isDefined() const
{
  return m_implem->isDefined();
}

/****************/
const EBISBox&
EBISLayoutImplem::operator[](const DataIndex& a_dit) const
{
  return m_ebisBoxes[a_dit];
}
/****************/
EBISLayoutImplem::EBISLayoutImplem()
{
  m_maxCoarseningRatio = 1;
  m_maxRefinementRatio = 1;
  m_defined = false;
}
/****************/
EBISLayoutImplem::~EBISLayoutImplem()
{
  CH_TIME("EBISLayoutImplem::~EBISLayoutImplem");
}
/****************/
VolIndex
EBISLayoutImplem::coarsen(const VolIndex& a_vof,
                          const int& a_ratio,
                          const DataIndex& a_datInd) const
{
  CH_assert(a_ratio > 0);
  CH_assert(a_ratio <= m_maxCoarseningRatio);
  CH_assert(a_ratio%2 == 0);

  //for ratio of 2, just use ebisbox
  const EBISBox& ebisBoxFine = m_ebisBoxes[a_datInd];
  VolIndex coarVoF = ebisBoxFine.coarsen(a_vof);
  //for ratio > 2, chase its tail
  int icoarlev = 0;
  for (int irat = 4; irat <= a_ratio; irat *= 2)
    {
      VolIndex fineVoF = coarVoF;
      const EBISLayout& ebisl = m_coarLevels[icoarlev];
      const EBISBox& ebisBox = ebisl[a_datInd];
      coarVoF = ebisBox.coarsen(fineVoF);
      icoarlev++;
    }
  return coarVoF;
}
/****************/
FaceIndex
EBISLayoutImplem::coarsen(const FaceIndex& a_face,
                          const int& a_ratio,
                          const DataIndex& a_datInd) const
{
  CH_assert(a_ratio > 0);
  CH_assert(a_ratio <= m_maxCoarseningRatio);
  CH_assert(a_ratio%2 == 0);

  //for ratio of 2, just use ebisbox
  const EBISBox& ebisBoxFine = m_ebisBoxes[a_datInd];
  FaceIndex coarFace = ebisBoxFine.coarsen(a_face);
  //for ratio > 2, chase its tail
  int icoarlev = 0;
  for (int irat = 4; irat <= a_ratio; irat *= 2)
    {
      FaceIndex fineFace = coarFace;
      const EBISLayout& ebisl = m_coarLevels[icoarlev];
      const EBISBox& ebisBox = ebisl[a_datInd];
      coarFace = ebisBox.coarsen(fineFace);
      icoarlev++;
    }
  return coarFace;
}
/****************/
Vector<VolIndex>
EBISLayoutImplem::refine(const VolIndex& a_vof,
                         const int& a_ratio,
                         const DataIndex& a_datInd) const
{
  CH_assert(a_ratio > 0);
  CH_assert(a_ratio <= m_maxRefinementRatio);
  CH_assert(a_ratio%2 == 0);

  //for ratio of 2, just use ebisbox
  const EBISBox& ebisBoxCoar = m_ebisBoxes[a_datInd];
  Vector<VolIndex> fineVoFs = ebisBoxCoar.refine(a_vof);
  //for ratio > 2, chase its tail
  int ifinelev = 0;
  for (int irat = 4; irat <= a_ratio; irat *= 2)
    {
      Vector<VolIndex> coarVoFs = fineVoFs;

      const EBISLayout& ebisl = m_fineLevels[ifinelev];
      const EBISBox& ebisBox = ebisl[a_datInd];
      fineVoFs.resize(0);
      for (int ivof = 0; ivof < coarVoFs.size(); ivof++)
        {
          fineVoFs.append(ebisBox.refine(coarVoFs[ivof]));
        }
      ifinelev++;
    }
  return fineVoFs;

}
/****************/
Vector<FaceIndex>
EBISLayoutImplem::refine(const FaceIndex& a_face,
                         const int& a_ratio,
                         const DataIndex& a_datInd) const
{
  CH_assert(a_ratio > 0);
  CH_assert(a_ratio <= m_maxRefinementRatio);
  CH_assert(a_ratio%2 == 0);

  //for ratio of 2, just use ebisbox
  const EBISBox& ebisBoxCoar2 = m_ebisBoxes[a_datInd];
  const EBISBox& ebisBoxFine2 = m_fineLevels[0][a_datInd];
  Vector<FaceIndex> fineFaces = ebisBoxCoar2.refine(a_face,ebisBoxFine2);
  //for ratio > 2, chase its tail
  int ifinelev = 0;
  for (int irat = 4; irat <= a_ratio; irat *= 2)
    {
      Vector<FaceIndex> coarFaces = fineFaces;

      const EBISBox& ebisBoxCoar = m_fineLevels[ifinelev    ][a_datInd];
      const EBISBox& ebisBoxFine = m_fineLevels[ifinelev + 1][a_datInd];
      fineFaces.resize(0);
      for (int iface = 0; iface < coarFaces.size(); iface++)
        {
          fineFaces.append(ebisBoxCoar.refine(coarFaces[iface],ebisBoxFine));
        }
      ifinelev++;
    }
  return fineFaces;

}
/****************/
void
EBISLayoutImplem::setMaxRefinementRatio(const int& a_maxRefine, const EBIndexSpace* ebisPtr)
{
  CH_assert(a_maxRefine % 2 == 0);
  CH_assert(a_maxRefine > 0);
  if (a_maxRefine < m_maxRefinementRatio)
    {
      //MayDay::Warning("why are you turning the max ref ratio down?");
      return;
    }

  m_maxRefinementRatio = a_maxRefine;
  //figure out how many levels i will need
  int nlevels = 0;
  for (int irat = 2; irat  <= a_maxRefine; irat *= 2)
    {
      nlevels++;
    }

  m_fineLevels.resize(nlevels);
  int irat = 2;
  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      ProblemDomain fineDomain = ebrefine(m_domain, irat);
      DisjointBoxLayout fineBLDomain;
      ebrefine(fineBLDomain, m_dblInputDom, irat);
      ebisPtr->fillEBISLayout(m_fineLevels[ilev], fineBLDomain, fineDomain, m_nghost);
      irat *= 2;
    }
}
/****************/
void
EBISLayoutImplem::setMaxCoarseningRatio(const int&                a_maxCoarsen,
                                        const EBIndexSpace* const a_ebisPtr)
{
  CH_assert(a_maxCoarsen % 2 == 0);
  CH_assert(a_maxCoarsen > 0);
  if (a_maxCoarsen < m_maxCoarseningRatio)
    {
      //MayDay::Warning("why are you turning the max coarsening ratio down?");
      return;
    }

  m_maxCoarseningRatio = a_maxCoarsen;
  //figure out how many levels i will need
  int nlevels = 0;
  for (int irat = 4; irat  <= a_maxCoarsen; irat *= 2)
    {
      nlevels++;
    }

  m_coarLevels.resize(nlevels);
  int irat = 2;
  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      ProblemDomain coarDomain = ebcoarsen(m_domain, irat);
      DisjointBoxLayout coarBLDomain;
      ebcoarsen(coarBLDomain, m_dblInputDom, irat);
      a_ebisPtr->fillEBISLayout(m_coarLevels[ilev], coarBLDomain, coarDomain, m_nghost);
      irat *= 2;
    }
}
/****************/
const EBISBox&
EBISLayout::operator[](const DataIndex& a_index) const
  {
    return m_implem->operator[](a_index);
  }

/****************/
void
EBISLayout::setMaxRefinementRatio(const int& a_maxRefine,const EBIndexSpace* ebisPtr )
{
  m_implem->setMaxRefinementRatio(a_maxRefine, ebisPtr);
}
/****************/
void
EBISLayout::setMaxCoarseningRatio(const int&                a_maxCoarsen,
                                  const EBIndexSpace* const a_ebisPtr)
{
  m_implem->setMaxCoarseningRatio(a_maxCoarsen,a_ebisPtr);
}
/****************/
EBISLayout::EBISLayout()
  : m_implem( RefCountedPtr<EBISLayoutImplem>( new EBISLayoutImplem() ) )
{
}
/****************/
const BoxLayout&
EBISLayout::getGrownLayout() const
{
  return m_implem->getGrownLayout();
}
/****************/
const DisjointBoxLayout&
EBISLayout::getDisjointLayout() const
{
  return m_implem->getDisjointLayout();
}
/****************/
EBISLayout::~EBISLayout()
{
}
/****************/
VolIndex
EBISLayout::coarsen(const VolIndex& a_vof,
                    const int& a_ratio,
                    const DataIndex& a_datInd) const
{
  return m_implem->coarsen(a_vof, a_ratio, a_datInd);
}
/****************/
Vector<VolIndex>
EBISLayout::refine(const VolIndex& a_vof,
                   const int& a_ratio,
                   const DataIndex& a_datInd) const
{
  return m_implem->refine(a_vof, a_ratio, a_datInd);
}
/****************/
Vector<FaceIndex>
EBISLayout::refine(const FaceIndex& a_face,
                   const int& a_ratio,
                   const DataIndex& a_datInd) const
{
  return m_implem->refine(a_face, a_ratio, a_datInd);
}
/****************/
int
EBISLayout::getMaxCoarseningRatio() const
{
  return m_implem->getMaxCoarseningRatio();
}
/****************/
int
EBISLayout::getMaxRefinementRatio() const
{
  return m_implem->getMaxRefinementRatio();
}
/****************/
int
EBISLayoutImplem::getMaxCoarseningRatio() const
{
  return m_maxCoarseningRatio;
}
/****************/
int
EBISLayoutImplem::getMaxRefinementRatio() const
{
  return m_maxRefinementRatio;
}
/****************/
#include "NamespaceFooter.H"
