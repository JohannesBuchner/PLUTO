#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INTERVAL_H_
#define _INTERVAL_H_
#include "BaseNamespaceHeader.H"

//
// Little struct for passing component ranges in code.
//

/// Structure for passing component ranges in code.
/**
   A class to encapsulate component ranges
*/
struct Interval
{
  Interval()
    :
    m_begin(0),
    m_end(-1)
  {}

  ///
  /**
      initialize with first and last component numbers
  */
  Interval(int a_firstComp, int a_lastComp)
    :
    m_begin(a_firstComp),
    m_end(a_lastComp)
  {}

  ///
  /**
   * Copy from another Interval
   */
  Interval& operator=(const Interval& p)
  {
    m_begin = p.m_begin;
    m_end = p.m_end;
    return *this;
  }

  void define(int a_firstComp,
              int a_lastComp)
  {
    m_begin = a_firstComp;
    m_end   = a_lastComp;
  }

  ///
  /**
     return first component number
  */
  int begin() const;

  ///
  /**
     return last component number
  */
  int end() const;

  ///
  /**
     return last - first + 1
   */
  int size() const
  {
    return m_end-m_begin+1;
  }

  ///
  bool contains(int a_val) const
  {
    return a_val >= m_begin && a_val <= m_end;
  }

  /// test for equality
  bool operator==(const Interval& p) const
  {
    return ( (m_begin == p.m_begin) && (m_end == p.m_end));
  }

private:
  int m_begin, m_end;
};

#ifndef WRAPPER
inline int Interval::begin() const
{
  return m_begin;
}

inline int Interval::end() const
{
  return m_end;
}
#endif

#include "BaseNamespaceFooter.H"
#endif
