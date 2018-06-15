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

#ifndef _BASEIVFACTORY_H_
#define _BASEIVFACTORY_H_

#include <cmath>
#include <cstdlib>
#include "Vector.H"
#include "EBISLayout.H"
#include "BaseIVFAB.H"
#include "DataIndex.H"
#include "NamespaceHeader.H"


/// Factory class to produce BaseIVFABs.
/**
   Factory class to produce BaseIVFABs.
   This is needed for LevelData construction.
 */
template <class T>
class BaseIVFactory
  : public DataFactory<BaseIVFAB <T> >
{
public:

  /// factory function.
  /**
      Creates a new baseivfab object
      and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer
      is passed to the user.
      Intersects a_sets[a_dit] from construction
      with input box because leveldata will expect
      that anyway when it tries to do a copy or a linearization.
  */
  virtual BaseIVFAB<T>* create(const Box& a_box, int a_ncomps,
                               const DataIndex& a_dit) const;

  ///
  /**
     fills internal data with inputs
  */
  BaseIVFactory(const EBISLayout& a_ebisl,
                const LayoutData<IntVectSet>& a_sets);

  ///
  /**
     builds BaseIVFAB objects defined at all boundary irregular cells.
  */
  BaseIVFactory(const EBISLayout& a_ebisl);

  ///
  virtual ~BaseIVFactory();

private:
  ///
  EBISLayout m_ebisl;

  ///
  LayoutData<IntVectSet> m_sets;

  //disallowed constructors and all that
  BaseIVFactory(const BaseIVFactory<T>& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }
  BaseIVFactory()
  {
    MayDay::Error("disallowed constructor");
  }
  void operator=(const BaseIVFactory<T>& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }

};

#include "NamespaceFooter.H"

//implementation
#include "BaseIVFactoryI.H"

#endif
