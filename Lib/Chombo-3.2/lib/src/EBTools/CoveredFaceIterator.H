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

#ifndef _COVEREDFACEITERATOR_H_
#define _COVEREDFACEITERATOR_H_

#include "IntVect.H"
#include "VolIndex.H"
#include "IntVectSet.H"
#include "EBISBox.H"
#include "Vector.H"
#include "LoHiSide.H"
#include "NamespaceHeader.H"

/// Iterator for all vofs within an IntVectSet and an EBISBox
/**
   Iterator for all vofs within an IntVectSet and an EBISBox
 */
class CoveredFaceIterator
{
public:
  ///
  CoveredFaceIterator(const IntVectSet&         a_ivs,
                      const EBISBox&            a_ebisBox,
                      int a_dir, Side::LoHiSide a_side);

  ///
  CoveredFaceIterator();

  ///
  ~CoveredFaceIterator();

  ///
  void define(const IntVectSet& a_ivs,
              const EBISBox& a_ebisBox,
              int a_dir, Side::LoHiSide a_side);

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

private:
  CoveredFaceIterator(const CoveredFaceIterator& a_input)
  {
    MayDay::Error("Invalid Operator");
  }
  void
  operator= (const CoveredFaceIterator& a_input)
  {
    MayDay::Error("Invalid Operator");
  }
  Vector<VolIndex> m_vols;
  int m_ivol;

  bool m_isDefined;
};

#include "NamespaceFooter.H"
#endif
