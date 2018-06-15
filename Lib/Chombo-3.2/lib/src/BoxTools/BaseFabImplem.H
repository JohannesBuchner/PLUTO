#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEFABIMPLEM_H_
#define _BASEFABIMPLEM_H_

#include "BaseFabMacros.H"
#include "CH_Timer.H"
#include "SliceSpec.H"
#include "MayDay.H"
#include "parstream.H"
#include "NamespaceHeader.H"

//
// BaseFab<Real> specializations
//

template < > void BaseFab<Real>::define();

template < > void BaseFab<Real>::undefine();

template < > void BaseFab<Real>::setVal (Real val);

template < > void BaseFab<int>::define();

template < > void BaseFab<int>::undefine();

template < > int BaseFab<int>::test();

template <class T> int BaseFab<T>::testBoxAndComp()
{
  MayDay::Warning("pretty minimal test");
  Box b(IntVect::Zero, IntVect::Unit);
  BaseFab<T> blerg;
  blerg.define(b, 1);

  if (blerg.nComp() != 1)
    {
      pout() << "ncomp busted" << endl;
      return -2;
    }
  if (blerg.box() != b)
    {
      pout() << "box return busted" << endl;
      return -2;
    }
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (blerg.size()[idir] != 2)
        {
          pout() << "size busted" <<endl;
          return -3;
        }
    }

  if (blerg.smallEnd() != IntVect::Zero)
    {
      pout() << "smallend busted" <<endl;
      return -4;
    }
  if (blerg.bigEnd() != IntVect::Unit)
    {
      pout() << "bigend busted" <<endl;
      return -5;
    }
  return 0;
}///
template <class T> int BaseFab<T>::test()
{
  int retval = testBoxAndComp();
  return retval;
}

//
// Implementation.=====================================
//

template <class T> inline BaseFab<T>::BaseFab()
  :
  m_domain(Box()),
  m_nvar(0),
  m_numpts(0),
  m_truesize(0),
  m_dptr(0),
  m_aliased(false)
{
}

template <class T> inline BaseFab<T>::BaseFab(const Box& a_bx,
                                              int        a_n,
                                              T*         a_alias)
  :
  m_domain(a_bx),
  m_nvar(a_n),
  m_numpts(a_bx.numPts()),
  m_truesize(a_bx.numPts() * a_n),
  m_dptr(0),
  m_aliased(false)
{
  if (a_alias != NULL)
  {
    m_dptr    = a_alias;
    m_aliased = true;
  }
  else
  {
    define();
  }
}

template <class T> inline BaseFab<T>::BaseFab(const Interval& a_comps,
                                              BaseFab<T>&     a_original)
  :
  m_domain(a_original.m_domain),
  m_nvar(a_comps.size()),
  m_numpts(a_original.m_numpts),
  m_truesize(a_original.m_numpts * m_nvar),
  m_dptr(a_original.dataPtr(a_comps.begin())),
  m_aliased(true)
{
}

template <class T> inline BaseFab<T>::~BaseFab()
{
  undefine();
}

template <class T> void BaseFab<T>::resize(const Box& a_b,
                                           int        a_n,
                                           T*         a_alias)
{
  // m_nvar   = a_n;
  // m_domain = a_b;
  // m_numpts = m_domain.numPts();
  //
  // if (m_dptr == 0)
  // {
  //   define();
  // }
  // else if (m_nvar*m_numpts > m_truesize)
  // {
  //   undefine();
  //   define();
  // }

  undefine();

  m_nvar   = a_n;
  m_domain = a_b;
  m_numpts = m_domain.numPts();

  if (a_alias != NULL)
  {
    m_dptr    = a_alias;
    m_aliased = true;
  }
  else
  {
    m_aliased = false;
    define();
  }
}

template <class T> inline void BaseFab<T>::define(const Interval& a_comps,
                                                  BaseFab<T>&     a_original)
{
  undefine();

  m_domain   = a_original.m_domain;
  m_numpts   = a_original.m_numpts;
  m_truesize = a_original.m_numpts*a_comps.size();
  m_nvar     = a_comps.size();
  m_dptr     = a_original.dataPtr(a_comps.begin());
  m_aliased  = true;
  // resize(a_original.m_domain, a_comps.size(),
  //        a_original.dataPtr(a_comps.begin()));
}

template <class T> inline void BaseFab<T>::clear()
{
  undefine();

  m_domain = Box();
  m_nvar   = 0;
  m_numpts = 0;
}

template <class T> inline int BaseFab<T>::nComp() const
{
  return m_nvar;
}

template <class T> Arena* BaseFab<T>::s_Arena = NULL;

template <class T> inline const Box& BaseFab<T>::box() const
{
  return m_domain;
}

template <class T> inline IntVect BaseFab<T>::size() const
{
  return m_domain.size();
}

template <class T> inline const IntVect& BaseFab<T>::smallEnd() const
{
  return m_domain.smallEnd();
}

template <class T> inline const IntVect& BaseFab<T>::bigEnd() const
{
  return m_domain.bigEnd();
}

template <class T> inline T& BaseFab<T>::operator () (const IntVect& a_p,
                                                      int            a_n)
{
  CH_assert(a_n >= 0);
  CH_assert(a_n < m_nvar);
  CH_assert(!(m_dptr == 0));
  CH_assert(m_domain.contains(a_p));

  long int ind1 = m_domain.index(a_p);
  long int ind2 = a_n * m_numpts;
  long int ind  = ind1 + ind2;

  return m_dptr[ind];
}

template <class T> inline T& BaseFab<T>::operator () (const IntVect& a_p)
{
  CH_assert(!(m_dptr == 0));
  CH_assert(m_domain.contains(a_p));

  return m_dptr[m_domain.index(a_p)];
}

template <class T> inline const T& BaseFab<T>::operator () (const IntVect& a_p,
                                                            int            a_n) const
{
  CH_assert(a_n >= 0);
  CH_assert(a_n < m_nvar);
  CH_assert(!(m_dptr == 0));
  CH_assert(m_domain.contains(a_p));

  return m_dptr[m_domain.index(a_p) + a_n * m_numpts];
}

template <class T> inline const T& BaseFab<T>::operator () (const IntVect& a_p) const
{
  CH_assert(!(m_dptr == 0));
  CH_assert(m_domain.contains(a_p));

  return m_dptr[m_domain.index(a_p)];
}

template <class T> inline void BaseFab<T>::getVal(T*             a_data,
                                                  const IntVect& a_pos,
                                                  int            a_n,
                                                  int            a_numcomp) const
{
  const int  loc     = m_domain.index(a_pos);
  const long size    = m_domain.numPts();

  CH_assert(!(m_dptr == 0));
  CH_assert(a_n >= 0 && a_n + a_numcomp <= m_nvar);

  for (int k = 0; k < a_numcomp; k++)
  {
    a_data[k] = m_dptr[loc+(a_n+k)*size];
  }
}

template <class T> inline void BaseFab<T>::getVal(T*             a_data,
                                                  const IntVect& a_pos) const
{
  getVal(a_data,a_pos,0,m_nvar);
}

template <class T> inline const int* BaseFab<T>::loVect() const
{
  return m_domain.loVect();
}

template <class T> inline const int* BaseFab<T>::hiVect() const
{
  return m_domain.hiVect();
}

template <class T> inline const int* BaseFab<T>::nCompPtr() const
{
  CH_assert(!(m_dptr == 0));

  return &m_nvar;
}

template <class T> inline T* BaseFab<T>::dataPtr(int a_n)
{
  CH_assert(!(m_dptr == 0));
  CH_assert((a_n >= 0) && (a_n < m_nvar));

  return &m_dptr[a_n * m_numpts];
}

template <class T> inline const T* BaseFab<T>::dataPtr(int a_n) const
{
  CH_assert(!(m_dptr == 0));
  CH_assert((a_n >= 0) && (a_n < m_nvar));

  return &m_dptr[a_n * m_numpts];
}

template <class T> inline bool BaseFab<T>::contains(const BaseFab<T>& a_fab) const
{
  return box().contains(a_fab.box()) && m_nvar <= a_fab.m_nvar;
}

template <class T> inline bool BaseFab<T>::contains (const Box& a_bx) const
{
  return box().contains(a_bx);
}

template <class T> inline void BaseFab<T>::setVal(T          a_x,
                                                  const Box& a_bx,
                                                  int        a_nstart,
                                                  int        a_numcomp)
{
  performSetVal(a_x,a_bx,a_nstart,a_numcomp);
}

template <class T> inline void BaseFab<T>::setVal(T          a_x,
                                                  const Box& a_bx,
                                                  int        a_n)
{
  performSetVal(a_x,a_bx,a_n,1);
}

template <class T> inline void BaseFab<T>::setVal(T   a_x,
                                                  int a_n)
{
  performSetVal(a_x,m_domain,a_n,1);
}

template <class T> inline void BaseFab<T>::setVal(T a_x)
{
  performSetVal(a_x,box(),0,m_nvar);
}

template <class T>
inline BaseFab<T>& BaseFab<T>::copy(const BaseFab<T>& a_src,
                                    const Box&        a_srcbox,
                                    int               a_srccomp,
                                    const Box&        a_destbox,
                                    int               a_destcomp,
                                    int               a_numcomp)
{
  CH_assert(a_srcbox.sameSize(a_destbox));
  CH_assert(a_src.box().contains(a_srcbox));
  CH_assert(m_domain.contains(a_destbox));
  CH_assert(a_srccomp >= 0 && a_srccomp+a_numcomp <= a_src.nComp());
  CH_assert(a_destcomp >= 0 && a_destcomp+a_numcomp <= m_nvar);

  performCopy(a_src,a_srcbox,a_srccomp,a_destbox,a_destcomp,a_numcomp);

  return *this;
}

template <class T>
inline BaseFab<T>& BaseFab<T>::copy(const BaseFab<T>& a_src,
                                    int               a_srccomp,
                                    int               a_destcomp,
                                    int               a_numcomp)
{
  CH_assert(a_srccomp  >= 0 && a_srccomp  + a_numcomp <= a_src.m_nvar);
  CH_assert(a_destcomp >= 0 && a_destcomp + a_numcomp <= m_nvar);

  Box overlap(m_domain);
  overlap &= a_src.m_domain;

  if (!overlap.isEmpty())
  {
    performCopy(a_src,overlap,a_srccomp,overlap,a_destcomp,a_numcomp);
  }

  return *this;
}

template <class T>
inline BaseFab<T>& BaseFab<T>::copy(const BaseFab<T>& a_src,
                                    const Box&        a_destbox)
{
  CH_assert(m_nvar <= a_src.m_nvar);
  CH_assert(m_domain.contains(a_destbox));

  Box overlap(a_destbox);
  overlap &= a_src.m_domain;

  if (!overlap.isEmpty())
  {
    performCopy(a_src,overlap,0,overlap,0,m_nvar);
  }

  return *this;
}

template <class T>
inline BaseFab<T>& BaseFab<T>::copy(const BaseFab<T>& a_src)
{
  CH_assert(m_nvar <= a_src.m_nvar);
  CH_assert(m_domain.sameType(a_src.m_domain));

  Box overlap(m_domain);
  overlap &= a_src.m_domain;

  if (!overlap.isEmpty())
  {
    performCopy(a_src,overlap,0,overlap,0,m_nvar);
  }

  return *this;
}

template <class T> inline void BaseFab<T>::copy(const Box&        a_RegionFrom,
                                                const Interval&   a_Cdest,
                                                const Box&        a_RegionTo,
                                                const BaseFab<T>& a_src,
                                                const Interval&   a_Csrc)
{
  if ((this == &a_src) && (a_RegionFrom == a_RegionTo) && (a_Cdest == a_Csrc) )
  {
    return;
  }

  CH_assert(a_Cdest.size() == a_Csrc.size());

  copy(a_src, a_RegionFrom, a_Csrc.begin(), a_RegionTo,
       a_Cdest.begin(), a_Cdest.size());
}

template <class T> inline BaseFab<T>& BaseFab<T>::shift(const IntVect& a_v)
{
  m_domain += a_v;

  return *this;
}

template <class T> inline BaseFab<T>& BaseFab<T>::shift(int a_idir,
                                                        int a_ncells)
{
  m_domain.shift(a_idir,a_ncells);

  return *this;
}

template <class T> inline BaseFab<T> & BaseFab<T>::shiftHalf(int a_idir,
                                                             int a_numHalfs)
{
  m_domain.shiftHalf(a_idir,a_numHalfs);

  return *this;
}

template <class T> inline BaseFab<T> & BaseFab<T>::shiftHalf(const IntVect& a_v)
{
  m_domain.shiftHalf(a_v);

  return *this;
}

template <class T> inline int BaseFab<T>::size(const Box&      a_box,
                                               const Interval& a_comps) const
{
  return a_box.numPts() * sizeof(T) * a_comps.size();
}

template <class T> inline void BaseFab<T>::linearOut(void*           a_buf,
                                                     const Box&      a_R,
                                                     const Interval& a_comps) const
{
  linearOut2(a_buf, a_R, a_comps);
}
template <class T> inline void* BaseFab<T>::linearOut2(void*           a_buf,
                                                       const Box&      a_R,
                                                       const Interval& a_comps) const
{
  T* buffer = (T*)a_buf;

  ForAllThisCBNN(T,a_R,a_comps.begin(),a_comps.size())
  {
    *buffer = thisR;
    ++buffer;
  } EndFor;
  return (void*)buffer;
}

template <class T> inline void BaseFab<T>::linearIn(void*           a_buf,
                                                    const Box&      a_R,
                                                    const Interval& a_comps)
{
  //  pout() << "basefab::linearin box = " <<  a_R << "comps = (" << a_comps.begin() << "," << a_comps.end() << ")" << endl;
  linearIn2(a_buf, a_R, a_comps);
}

template <class T> inline void* BaseFab<T>::linearIn2(void*           a_buf,
                                                     const Box&      a_R,
                                                     const Interval& a_comps)
{
  T* buffer = (T*)a_buf;

  ForAllThisBNN(T,a_R,a_comps.begin(),a_comps.size())
  {
    thisR = *buffer;
    ++buffer;
  } EndFor;

  return (void*) buffer;
}

template <class T> inline void BaseFab<T>::define()
{
  CH_assert(m_nvar > 0);
  CH_assert(m_dptr == 0);
  // CH_assert(m_numpts > 0); // OK if box is empty
  CH_assert(!m_aliased);
  //CH_assert(!(The_FAB_Arena == 0));// not a sufficient test !!!

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

  if (s_Arena == NULL)
  {
    MayDay::Error("malloc in basefab failed");
  }

  m_truesize = m_nvar * m_numpts;
  if (m_truesize > 0)
    {
      m_dptr     = static_cast<T*>(s_Arena->alloc(m_truesize * sizeof(T)));

#ifdef CH_USE_MEMORY_TRACKING
      s_Arena->bytes += m_truesize * sizeof(T) + sizeof(BaseFab<T>);
      CH_assert(s_Arena->bytes > 0);
      if (s_Arena->bytes > s_Arena->peak)
        {
          s_Arena->peak = s_Arena->bytes;
        }
#endif
    }
  else
    { // m_truesize == 0: allocating no space
      m_dptr = NULL;
    }

  //
  // Now call T::T() on the raw memo'ry so we have valid Ts.
  //
  T* ptr = m_dptr;

  for (long int i = 0; i < m_truesize; i++, ptr++)
  {
    new (ptr) T;
  }
}

template <class T> inline void BaseFab<T>::undefine()
{
  //CH_assert(!(The_FAB_Arena == 0));
  //
  // Call T::~T() on the to-be-destroyed memory.
  //
  if (m_aliased)
  {
    m_dptr = 0;
    return;
  }

  if (m_dptr == 0)
  {
    return;
  }

  T* ptr = m_dptr;

  for (long int i = 0; i < m_truesize; i++, ptr++)
  {
    ptr->~T();
  }

  s_Arena->free(m_dptr);

#ifdef CH_USE_MEMORY_TRACKING
  s_Arena->bytes -= m_truesize * sizeof(T) + sizeof(BaseFab<T>);
  CH_assert(s_Arena->bytes >= 0);
#endif

  m_dptr = 0;
}

template <class T> inline std::string BaseFab<T>::name()
{
  std::string rtn = (typeid(T)).name();

  return rtn;
}

template <class T>
inline void BaseFab<T>::performCopy(const BaseFab<T>& a_src,
                                    const Box&        a_srcbox,
                                    int               a_srccomp,
                                    const Box&        a_destbox,
                                    int               a_destcomp,
                                    int               a_numcomp)
{
  CH_assert(a_src.box().contains(a_srcbox));
  CH_assert(box().contains(a_destbox));
  CH_assert(a_destbox.sameSize(a_srcbox));
  CH_assert(a_srccomp  >= 0 && a_srccomp  + a_numcomp <= a_src.nComp());
  CH_assert(a_destcomp >= 0 && a_destcomp + a_numcomp <= nComp());
  // CH_TIME("BaseFab::performCopy")
  ForAllThisBNNXCBN(T, a_destbox, a_destcomp, a_numcomp, a_src, a_srcbox, a_srccomp)
  {
    thisR = a_srcR;
  } EndForTX
}

template <class T> inline void BaseFab<T>::performSetVal(T          a_x,
                                                         const Box& a_bx,
                                                         int        a_nstart,
                                                         int        a_numcomp)
{
  CH_assert(m_domain.contains(a_bx));
  CH_assert(a_nstart >= 0 && a_nstart + a_numcomp <= m_nvar);

  if (a_bx == m_domain)
  {
    T* data = &m_dptr[a_nstart * m_numpts];

    for (long i = 0, N = a_numcomp * m_numpts; i < N; i++)
    {
      *data++ = a_x;
    }
  }
  else
    {
      ForAllThisBNN(T,a_bx,a_nstart,a_numcomp)
      {
        thisR = a_x;
      } EndFor
    }
}

template <class T>
bool BaseFab<T>::isAliased() const
{
  return m_aliased;
}

template <class T> void
BaseFab<T>::degenerate(  BaseFab<T>& a_slice,
                         const SliceSpec& a_sliceSpec ) const
{
  bool outofbounds;
  Box degenerateBox;
  this->box().degenerate( degenerateBox, a_sliceSpec, &outofbounds );
  if ( outofbounds )
    {
      MayDay::Error( "Tried to slice out-of-bounds." );
    }
  a_slice.define( degenerateBox, this->nComp() );
  a_slice.copy( *this, degenerateBox );
}

#include "NamespaceFooter.H"

#endif
