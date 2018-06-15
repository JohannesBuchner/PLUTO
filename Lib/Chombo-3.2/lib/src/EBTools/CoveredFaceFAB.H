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

#ifndef _COVEREDFACEFAB_H_
#define _COVEREDFACEFAB_H_

#include "IntVect.H"
#include "VolIndex.H"
#include "IntVectSet.H"
#include "EBISBox.H"
#include "BaseIVFAB.H"
#include "Vector.H"
#include "LoHiSide.H"
#include "NamespaceHeader.H"

///holder for data over covered faces
class CoveredFaceFAB
{
public:

  ///
  int nComp() const;

  ///
  int getIndex(int a_index, Side::LoHiSide a_sd) const;

  ///
  void clear();

  ///
  CoveredFaceFAB();

  ///
  ~CoveredFaceFAB();

  ///
  CoveredFaceFAB(const IntVectSet&  a_ivs,
                 const EBISBox&     a_ebisBox,
                 int                a_nvar);

  ///
  void define(const IntVectSet&  a_ivs,
              const EBISBox&     a_ebisBox,
              int                a_nvar);

  ///
  BaseIVFAB<Real>&  operator()
    (int a_idir, Side::LoHiSide a_sd);

  ///
  const BaseIVFAB<Real>& operator()
    (int a_idir, Side::LoHiSide a_sd) const;

  ///
  const IntVectSet& getIVS(int a_idir, Side::LoHiSide a_sd) const;

  ///
  bool isDefined() const;

  ///
  void setVal(const Real& a_val);

  ///
  void
  copy(const Box& Rfrom, const Interval& Cdest,
       const Box& Rto, const CoveredFaceFAB& src,
       const Interval& Csrc);

private:

  CoveredFaceFAB(const CoveredFaceFAB& a_input)
  {
    MayDay::Error("Invalid Operator");
  }
  void
  operator= (const CoveredFaceFAB& a_input)
  {
    MayDay::Error("Invalid Operator");
  }

  EBISBox          m_ebisBox;
  IntVectSet       m_sets[2*SpaceDim];
  BaseIVFAB<Real>  m_data[2*SpaceDim];
  int              m_nComp;
  bool             m_isDefined;

};

#include "NamespaceFooter.H"
#endif
