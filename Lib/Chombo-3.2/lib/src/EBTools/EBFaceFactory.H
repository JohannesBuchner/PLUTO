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

#ifndef _EBFACEFACTORY_H_
#define _EBFACEFACTORY_H_

#include <cmath>
#include <cstdlib>
#include "Vector.H"
#include "EBISLayout.H"
#include "EBFaceFAB.H"
#include "DataIndex.H"
#include "NamespaceHeader.H"

/// Factory class to produce BaseEBFaceFABs.
/**
   Factory class to produce BaseEBFaceFABs.
   This is needed for LevelData construction.
 */
class EBFaceFactory
  : public DataFactory<EBFaceFAB>
{
public:
  /// factory function.
  /**
      Creates a new baseebcellfab object
      and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user.
  */
  virtual EBFaceFAB* create(const Box& a_box, int a_ncomps,
                            const DataIndex& a_dit) const;

  ///
  /**
     create the factory with an ispace.  calls full define function
  */
  EBFaceFactory(const EBISLayout& a_ebisl,
                const int& a_idir);

  ///
  virtual ~EBFaceFactory();

private:
  EBISLayout m_ebisl;
  int m_idir;

  //disallowed constructors and all that
  EBFaceFactory(const EBFaceFactory& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }
  EBFaceFactory()
  {
    MayDay::Error("disallowed constructor");
  }
  void operator=(const EBFaceFactory& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }
};

#include "NamespaceFooter.H"
#endif
