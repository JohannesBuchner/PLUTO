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

#ifndef _BASEEBCELLFACTORY_H_
#define _BASEEBCELLFACTORY_H_

#include <cmath>
#include <cstdlib>
#include "Vector.H"
#include "EBISLayout.H"
#include "BaseEBCellFAB.H"
#include "DataIndex.H"
#include "NamespaceHeader.H"


/// Factory class to produce BaseEBCellFABs.
/**
   Factory class to produce BaseEBCellFABs.
   This is needed for LevelData construction.
 */
template <class T>
class BaseEBCellFactory
  : public DataFactory<BaseEBCellFAB<T> >
{
public:
  /// factory function.
  /**
      Creates a new baseebcellfab object
      and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user.
  */
  virtual BaseEBCellFAB<T>* create(const Box& a_box, int a_ncomps,
                                   const DataIndex& a_dit) const;

  ///
  /**
     create the factory with an ispace.  calls full define function
  */
  BaseEBCellFactory(const EBISLayout& a_ebisl);

  ///
  virtual ~BaseEBCellFactory();

private:
  EBISLayout m_ebisl;

  //disallowed constructors and all that
  BaseEBCellFactory(const BaseEBCellFactory& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }
  BaseEBCellFactory()
  {
    MayDay::Error("disallowed constructor");
  }
  void operator=(const BaseEBCellFactory& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }

};

#include "NamespaceFooter.H"
#include "BaseEBCellFactoryI.H"

#endif
