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

#ifndef _EBCELLFACTORY_H_
#define _EBCELLFACTORY_H_

#include <cmath>
#include <cstdlib>
#include "Vector.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "DataIndex.H"
#include "NamespaceHeader.H"

/// Factory class to produce EBCellFABs.
/**
   Factory class to produce EBCellFABs.
   This is needed for LevelData construction.
 */
class EBCellFactory
  : public DataFactory<EBCellFAB>
{
public:
  /// factory function.
  /**
      Creates a new baseebcellfab object
      and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user.
  */
  virtual EBCellFAB* create(const Box& a_box, int a_ncomps,
                            const DataIndex& a_dit) const;

  ///
  /**
     create the factory with an ispace.  calls full define function
  */
  EBCellFactory(const EBISLayout& a_ebisl);

  ///
  virtual ~EBCellFactory();

private:
  EBISLayout m_ebisl;

  //disallowed constructors and all that
  EBCellFactory(const EBCellFactory& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }
  EBCellFactory()
  {
    MayDay::Error("disallowed constructor");
  }
  void operator=(const EBCellFactory& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }

};

#include "NamespaceFooter.H"
#endif
