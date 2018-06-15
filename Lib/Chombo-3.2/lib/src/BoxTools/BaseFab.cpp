#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstring>

#include "BaseFab.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"

Real BaseFabRealSetVal = BASEFAB_REAL_SETVAL;

template < > int BaseFab<int>::test()
{
  int retbox = testBoxAndComp();
  if (retbox != 0)
    {
      pout() << "testboxandcomp failed" << endl;
      return retbox;
    }

  Box b(IntVect::Zero, IntVect::Unit);
  BaseFab<int> blerg;
  blerg.define(b, 1);
  int val = 4;
  blerg.setVal(val);
  for (BoxIterator bit(b); bit.ok();++bit)
    {
      if (blerg(bit(), 0) != val)
        {
          pout() << "setval or index busted" << endl;
          return -1;
        }
    }

  Box bcop(IntVect::Unit, 2*IntVect::Unit);
  BaseFab<int> bfcopy(bcop, 1);
  bfcopy.setVal(2*val);
  bfcopy.copy(blerg, 0, 0, 1);
  Box binter =bcop;
  binter &= b;
  for (BoxIterator bit(binter); bit.ok();++bit)
    {
      if (bfcopy(bit(), 0) != val)
        {
          pout() << "copy busted" << endl;
          return -2;
        }
    }

  return 0;
}
template < > void BaseFab<Real>::define()
{
  CH_assert(m_nvar > 0);
  CH_assert(m_dptr == 0);
  // CH_assert(m_numpts > 0);
  // petermc, 8 Dec 2010:  Allow empty BaseFab<Real>.
  if (m_numpts == 0) return;
  CH_assert(!m_aliased);
  //CH_assert(!(The_FAB_Arena == 0)); // not a sufficient test!!!

#ifdef CH_USE_MEMORY_TRACKING
  if (s_Arena == NULL)
  {
    s_Arena = new BArena(name().c_str());
  }
#else
  if (s_Arena == NULL)
  {
    s_Arena = new BArena("");
  }
#endif

  // if (s_Arena == NULL)
  // {
  //   MayDay::Error("malloc in basefab failed");
  // }

  m_truesize = m_nvar * m_numpts;
  m_dptr     = static_cast<Real*>(s_Arena->alloc(m_truesize * sizeof(Real)));

#ifdef CH_USE_MEMORY_TRACKING
  s_Arena->bytes += m_truesize * sizeof(Real) + sizeof(BaseFab<Real>);
  CH_assert(s_Arena->bytes >= 0);
  if (s_Arena->bytes > s_Arena->peak)
  {
    s_Arena->peak = s_Arena->bytes;
  }
#endif

#ifdef CH_USE_SETVAL
  setVal(BaseFabRealSetVal);
#endif
}

template < > void BaseFab<int>::define()
{
  CH_assert(m_nvar > 0);
  CH_assert(m_dptr == 0);
  // CH_assert(m_numpts > 0);
  // petermc, 10 Apr 2012:  Allow empty BaseFab<int>.
  CH_assert(!m_aliased);
  //CH_assert(!(The_FAB_Arena == 0)); // not a sufficient test!!!

#ifdef CH_USE_MEMORY_TRACKING
  if (s_Arena == NULL)
  {
    s_Arena = new BArena(name().c_str());
  }
#else
  if (s_Arena == NULL)
  {
    s_Arena = new BArena("");
  }
#endif

  // if (s_Arena == NULL)
  // {
  //   MayDay::Error("malloc in basefab failed");
  // }

  m_truesize = m_nvar * m_numpts;
  m_dptr     = static_cast<int*>(s_Arena->alloc(m_truesize * sizeof(int)));

#ifdef CH_USE_MEMORY_TRACKING
  s_Arena->bytes += m_truesize * sizeof(int) + sizeof(BaseFab<int>);
  CH_assert(s_Arena->bytes >= 0);
  if (s_Arena->bytes > s_Arena->peak)
  {
    s_Arena->peak = s_Arena->bytes;
  }
#endif
}

template < > void BaseFab<Real>::undefine()
{
  if (m_aliased)
  {
    m_dptr = 0;
    return;
  }

  if (m_dptr == 0)
  {
    return;
  }

  s_Arena->free(m_dptr);

#ifdef CH_USE_MEMORY_TRACKING
  s_Arena->bytes -= m_truesize * sizeof(Real) + sizeof(BaseFab<Real>);
  CH_assert(s_Arena->bytes >= 0);
#endif

  m_dptr = 0;
}

template < > void BaseFab<int>::undefine()
{
  if (m_aliased)
  {
    m_dptr = 0;
    return;
  }

  if (m_dptr == 0)
  {
    return;
  }

  s_Arena->free(m_dptr);

#ifdef CH_USE_MEMORY_TRACKING
  s_Arena->bytes -= m_truesize * sizeof(int) + sizeof(BaseFab<int>);
  CH_assert(s_Arena->bytes >= 0);
#endif

  m_dptr = 0;
}

template < > void BaseFab<Real>::setVal(Real a_val)
{
  if (a_val == 0)
  {
    memset(m_dptr, 0, m_truesize*sizeof(Real));
  }
  else
  {
    Real* end = m_dptr + m_truesize;
    for (Real* v = m_dptr; v<end; v++)
    {
      *v = a_val;
    }
  }
}
#include "NamespaceFooter.H"
