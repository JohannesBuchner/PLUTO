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

#ifndef _BASEIFFACTORY_H_
#define _BASEIFFACTORY_H_

#include <cmath>
#include <cstdlib>
#include "Vector.H"
#include "EBISLayout.H"
#include "BaseIFFAB.H"
#include "DataIndex.H"
#include "NamespaceHeader.H"


/// Factory class to produce BaseIVFABs.
/**
   Factory class to produce BaseIVFABs.
   This is needed for LevelData construction.
 */
template <class T>
class BaseIFFactory
  : public DataFactory<BaseIFFAB<T> >
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
  virtual BaseIFFAB<T>* create(const Box& a_box, int a_ncomps,
                               const DataIndex& a_dit) const;

  ///
  /**
     fills internal data with inputs
  */
  BaseIFFactory(const EBISLayout& a_ebisl,
                const LayoutData<IntVectSet>& a_sets,
                const int& a_facedir );

  ///
  virtual ~BaseIFFactory();

private:
  ///
  EBISLayout m_ebisl;

  ///
  LayoutData<IntVectSet> m_sets;

  ///
  int m_facedir;

  //disallowed constructors and all that
  BaseIFFactory(const BaseIFFactory<T>& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }
  BaseIFFactory()
  {
    MayDay::Error("disallowed constructor");
  }
  void operator=(const BaseIFFactory<T>& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }

};

#include "NamespaceFooter.H"

//implementation
#include "BaseIFFactoryI.H"
#endif
