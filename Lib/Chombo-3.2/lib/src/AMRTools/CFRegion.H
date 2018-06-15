#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CFREGION_H_
#define _CFREGION_H_

#include "REAL.H"
#include "BaseFab.H"
#include "FArrayBox.H"
#include "LayoutData.H"
#include "DisjointBoxLayout.H"
#include "CFIVS.H"
#include "NamespaceHeader.H"


/// Class that represents the edge region around a DisjointBoxLayout
/**


  */
class CFRegion
{
public:

  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  CFRegion():m_defined(false)
  {
  }

  ///
  /**
     Destructor.
  */
  virtual ~CFRegion()
  {
  }

  CFRegion(const CFRegion& a_rhs)
  {
    *this = a_rhs;
  }

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

  */
  CFRegion(const DisjointBoxLayout& a_grids,
           const ProblemDomain&     a_domain)
  {
    define(a_grids, a_domain);
  }
  ///
  /**
     Defines this object.  Existing information is overriden.
  */
  void
  define(const DisjointBoxLayout& a_grids,
         const ProblemDomain&     a_domain);


  const CFRegion& operator=(const CFRegion& a_rhs);

  const CFIVS& loCFIVS(const DataIndex& a_dit, int dir);
  const CFIVS& hiCFIVS(const DataIndex& a_dit, int dir);

  void coarsen(int refRatio);

protected:

  LayoutData<CFIVS>       m_loCFIVS[SpaceDim];
  LayoutData<CFIVS>       m_hiCFIVS[SpaceDim];
  bool m_defined;
};

#include "NamespaceFooter.H"
#endif
