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

#ifndef _EBFLUXFACTORY_H_
#define _EBFLUXFACTORY_H_

#include <cmath>
#include <cstdlib>
#include "Vector.H"
#include "EBISLayout.H"
#include "EBFluxFAB.H"
#include "DataIndex.H"
#include "NamespaceHeader.H"

/// Factory class to produce BaseEBFluxFABs.
/**
   Factory class to produce BaseEBFluxFABs.
   This is needed for LevelData construction.
 */
class EBFluxFactory
  : public DataFactory<EBFluxFAB>
{
public:
  /// factory function.
  /**
      Creates a new baseebcellfab object
      and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user.
  */
  virtual EBFluxFAB* create(const Box& a_box, int a_ncomps,
                            const DataIndex& a_dit) const;

  ///
  /**
     create the factory with an ispace.  calls full define function
  */
  EBFluxFactory(const EBISLayout& a_ebisl);

  ///
  virtual ~EBFluxFactory();

private:
  EBISLayout m_ebisl;

};

#include "NamespaceFooter.H"
#endif
