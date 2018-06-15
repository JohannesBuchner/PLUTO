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

#ifndef _IRREGFABFACTORY_H_
#define _IRREGFABFACTORY_H_

#include <cmath>
#include <cstdlib>
#include "Vector.H"
#include "EBISLayout.H"
#include "IrregFAB.H"
#include "DataIndex.H"
#include "NamespaceHeader.H"

/// Factory class to produce IrregFABs.
/**
   Factory class to produce IrregFABs.
   This is needed for LevelData construction.
*/
class IrregFABFactory
  : public DataFactory<IrregFAB>
{
public:
  /// factory function.
  /**
     Creates a new baseIrregFAB object
     and returns a pointer to it.  Responsiblitly
     for calling operator 'delete' on this pointer is passed to the user.
  */
  virtual IrregFAB* create(const Box& a_box, int a_ncomps,
                           const DataIndex& a_dit) const;

  ///
  /**
     create the factory with an ispace.  calls full define function
  */
  IrregFABFactory(const EBISLayout& a_ebisl);

  ///
  virtual ~IrregFABFactory();

private:
  EBISLayout m_ebisl;

  //disallowed constructors and all that
  IrregFABFactory(const IrregFABFactory& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }
  IrregFABFactory()
  {
    MayDay::Error("disallowed constructor");
  }
  void operator=(const IrregFABFactory& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }

};

#include "NamespaceFooter.H"
#endif
