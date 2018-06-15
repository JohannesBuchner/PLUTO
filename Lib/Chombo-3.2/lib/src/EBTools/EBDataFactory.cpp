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

#include "EBDataFactory.H"
#include "NamespaceHeader.H"

EBDataFactory::EBDataFactory()
{
}

EBDataFactory::~EBDataFactory()
{
}

EBData*
EBDataFactory::create(const Box& box, int ncomps,
                       const DataIndex& a_datInd) const
{
  EBData* newEBData =  new EBData(box, ncomps);

  return newEBData;
}
#include "NamespaceFooter.H"
