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

#ifndef _VOFITERATOR_H_
#define _VOFITERATOR_H_

#include "IntVect.H"
#include "VolIndex.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "NamespaceHeader.H"
class EBGraph;

/// Iterator for all vofs within an IntVectSet and an Ebgraph
/**
   Iterator for all vofs within an IntVectSet and an Ebgraph
 */
class VoFIterator
{
public:
  ///
  VoFIterator(const IntVectSet& a_ivs,
              const EBGraph& a_ebgraph);

  ///
  VoFIterator();

  ///
  ~VoFIterator();

  ///
  void define(const IntVectSet& a_ivs,
              const EBGraph& a_ebgraph);

  ///
  void reset();

  ///
  void operator++();

  ///
  const VolIndex& operator() () const;

  ///
  bool ok() const;

  ///
  bool isDefined() const;

  const Vector<VolIndex> getVector() const;

  int size() const
  {
    return m_vols.size();
  }
private:
  bool m_isDefined;
  Vector<VolIndex> m_vols;
  int m_ivol;
};

#include "NamespaceFooter.H"
#endif
