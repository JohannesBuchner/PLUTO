#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#include "EBGraphFactory.H"
#include "NamespaceHeader.H"

EBGraphFactory::EBGraphFactory(const ProblemDomain& a_domain)
{
  m_domain = a_domain;
}

EBGraphFactory::~EBGraphFactory()
{
}

EBGraph*
EBGraphFactory::create(const Box& box, int ncomps,
                       const DataIndex& a_datInd) const
{
  Box region = box & m_domain;
  EBGraph* newEBGraph =  new EBGraph(region);
  newEBGraph->setDomain(m_domain);
  return newEBGraph;
}
#include "NamespaceFooter.H"
