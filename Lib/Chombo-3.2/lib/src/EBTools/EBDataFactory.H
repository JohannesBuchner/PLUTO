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

#ifndef _EBDATAFACTORY_H_
#define _EBDATAFACTORY_H_

#include "Box.H"
#include "LevelData.H"
#include "EBData.H"
#include "EBGraph.H"
#include "NamespaceHeader.H"

/// Factory class to produce EBData(es)
/**
   This is needed for LevelData construction.
 */
class EBDataFactory
  : public DataFactory<EBData>
{
public:
  ///
  /**
      Creates an undefined EBData.
  */
  virtual EBData* create(const Box& a_box, int a_ncomps,
                          const DataIndex& a_datInd) const;

  ///
  /**
   */
  EBDataFactory();

  ///
  /**
   */
  virtual ~EBDataFactory();

private:
};
#include "NamespaceFooter.H"
#endif
