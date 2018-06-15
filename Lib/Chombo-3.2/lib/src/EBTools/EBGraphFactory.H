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

#ifndef _EBGRAPHFACTORY_H_
#define _EBGRAPHFACTORY_H_

#include "Box.H"
#include "BoxLayoutData.H"
#include "EBGraph.H"
#include "NamespaceHeader.H"

/// Factory class to produce EBGraph(es)
/**
   This is needed for LevelData construction.
 */
class EBGraphFactory
  : public DataFactory<EBGraph>
{
public:
  ///
  /**
      Creates an undefined EBGraph.
  */
  virtual EBGraph* create(const Box& a_box, int a_ncomps,
                          const DataIndex& a_datInd) const;

  ///
  /**
   */
  EBGraphFactory(const ProblemDomain& a_domain);

  ///
  /**
   */
  virtual ~EBGraphFactory();

private:
  EBGraphFactory()
  {
    MayDay::Error("ebisboxfactory uses strong construction");
  }
  ProblemDomain m_domain;
};
#include "NamespaceFooter.H"
#endif
