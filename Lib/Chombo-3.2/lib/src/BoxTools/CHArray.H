#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CHARRAY_H_
#define _CHARRAY_H_

/******************************************************************************/
/**
 * \file
 *
 * \brief A multi-dimensional array class for Chombo
 *
 *//***************************************************************************/
// search for 'class CHArray'

#include <iomanip>

#include "CHArray_fwd.H"

#ifdef CH_SPACEDIM
// In Chombo land
  #include "REAL.H"
  #include "Box.H"
  #include "NamespaceHeader.H"
#else
// Probably testing land
  #include <cassert>
//   #define CH_assert(cond) (void)0
  #define CH_assert(cond) assert(cond)
  #define SpaceDim 2
  #define D_DECL6(a,b,c,d,e,f) a,b
  #define D_TERM6(a,b,c,d,e,f) a b
  typedef double Real;
  class IntVect
  {
  public:
    IntVect(D_DECL6(const int i0, const int i1, const int i2,
                    const int i3, const int i4, const int i5))
      {
        D_TERM6(vec[0] = i0;, vec[1] = i1;, vec[2] = i2;,
                vec[3] = i3;, vec[4] = i4;, vec[5] = i5;)
      }
    IntVect(const IntVect &iv)
      {
        D_TERM6(vec[0] = iv.vec[0];, vec[1] = iv.vec[1];, vec[2] = iv.vec[2];,
                vec[3] = iv.vec[3];, vec[4] = iv.vec[4];, vec[5] = iv.vec[5];)
      }
    int &operator[](const int i)
    {
      return vec[i];
    }
    const int &operator[](const int i) const
    {
      return vec[i];
    }
    int product() const
      {
        int prod = 1;
        for (int i = 0; i != SpaceDim; ++i) prod *= vec[i];
        return prod;
      }
    int vec[SpaceDim];
  };
  class Box
  {
  public:
    Box(const IntVect &lb, const IntVect &ub)
      : m_lb(lb), m_ub(ub)
      { }
    int smallEnd(const int dir) const
      {
        return m_lb[dir];
      }
    int bigEnd(const int dir) const
      {
        return m_ub[dir];
      }
    IntVect size() const
      {
        return IntVect(D_DECL6(m_ub[0] - m_lb[0] + 1,
                               m_ub[1] - m_lb[1] + 1,
                               m_ub[2] - m_lb[2] + 1,
                               m_ub[3] - m_lb[3] + 1,
                               m_ub[4] - m_lb[4] + 1,
                               m_ub[5] - m_lb[5] + 1));
      }
    IntVect m_lb;
    IntVect m_ub;
  };
#endif

// How to index an IntVect.  May evaluate to, e.g., iv[0], iv[1], iv[2] or
// iv[2], iv[1], iv[0]
# define D_DEFIV D_DECL6(DimT(iv[Indexer::ixIV(0)]), \
                         DimT(iv[Indexer::ixIV(1)]), \
                         DimT(iv[Indexer::ixIV(2)]), \
                         DimT(iv[Indexer::ixIV(3)]), \
                         DimT(iv[Indexer::ixIV(4)]), \
                         DimT(iv[Indexer::ixIV(5)]))

// How to index a Box.  Same as above but determines a range.
# define D_DEFBOX D_DECL6(DimT(box.smallEnd(Indexer::ixIV(0)), \
                                 box.bigEnd(Indexer::ixIV(0))), \
                          DimT(box.smallEnd(Indexer::ixIV(1)), \
                                 box.bigEnd(Indexer::ixIV(1))), \
                          DimT(box.smallEnd(Indexer::ixIV(2)), \
                                 box.bigEnd(Indexer::ixIV(2))), \
                          DimT(box.smallEnd(Indexer::ixIV(3)), \
                                 box.bigEnd(Indexer::ixIV(3))), \
                          DimT(box.smallEnd(Indexer::ixIV(4)), \
                                 box.bigEnd(Indexer::ixIV(4))), \
                          DimT(box.smallEnd(Indexer::ixIV(5)), \
                                 box.bigEnd(Indexer::ixIV(5))))

// How to index an IntVect.  May evaluate to, e.g., iv[0], iv[1], iv[2] or
// iv[2], iv[1], iv[0]
# define D_IXIV D_DECL6(iv[Indexer::ixIV(0)], \
                        iv[Indexer::ixIV(1)], \
                        iv[Indexer::ixIV(2)], \
                        iv[Indexer::ixIV(3)], \
                        iv[Indexer::ixIV(4)], \
                        iv[Indexer::ixIV(5)])

namespace ArSp  // CHArray support
{

// An index refers to a single dimension and a size to the overall
// size of the array.  Except for 1D arrays, indexes are smaller than
// sizes by division by at least 1 stride.

typedef int      IIx_t;               ///< Type of int for an index.  Must be
                                      ///< signed to allow for lower bounds != 0
typedef unsigned USz_t;               ///< Type of unsigned int for a size
}

/* This enum had to be put in CHArray_fwd.H because it is used for default
 * template parameters.  Repeated here for informational purposes
/|/  Array configuration parameters
enum
{
  ArZeroRow,                          /|/< Zero lower bound and row-storage
  ArZeroCol,                          /|/< Zero lower bound and column-storage
  ArRangeRow,                         /|/< Range of subscripts and row-storage
  ArRangeCol                          /|/< Range of subscripts and column-
                                      /|/< storage
};
 */

/*====================================================================*/
///  Defines a range of subscripts
/*====================================================================*/

struct CHRange
{
  CHRange(const ArSp::IIx_t a_dimE)
    : dimB(0), dimE(a_dimE - 1)
    {
      CH_assert(a_dimE >= 0);
    }
  CHRange(const ArSp::IIx_t a_dimB, const ArSp::IIx_t a_dimE)
    : dimB(a_dimB), dimE(a_dimE)
    {
      CH_assert(a_dimE >= a_dimB);
    }
  ArSp::IIx_t dimB;
  ArSp::IIx_t dimE;
};

/*====================================================================*/
/**  \name Support for describing dimensions by a range of subscripts
 *//*=================================================================*/
//@{

namespace ArSp  // CHArray support
{

///  Get the size of a dimension
template <typename T>
inline USz_t sizeOfDim(const T &dim)
{
  return (USz_t)dim;
}

///
template <>
inline USz_t sizeOfDim<CHRange>(const CHRange &dim)
{
  return (USz_t)(dim.dimE - dim.dimB + 1);
}
}

///  Multiply two dimensions
inline ArSp::USz_t operator*(const CHRange &a, const CHRange &b)
{
  return (ArSp::sizeOfDim(a)*ArSp::sizeOfDim(b));
}

///
inline ArSp::USz_t operator*(const ArSp::USz_t &i, const CHRange &a)
{
  return i*ArSp::sizeOfDim(a);
}

///
inline ArSp::USz_t operator*(const CHRange &a, const ArSp::USz_t &i)
{
  return ArSp::sizeOfDim(a)*i;
}
//@}

namespace ArSp  // CHArray support
{
/*====================================================================*/
///  Test for class type
/**  Use the SFINAE principle as described in 15.2.2 "Determining
 *   Class Types" in Vandevoorde and Josuttis "C++ Templates" book to
 *   see if T is class type.  E.g.,
 *     TypeTr<int>::IsClass == false
 *     TypeTr<Box>::IsClass == true
 *//*=================================================================*/

template <typename T>
class TypeTr
{
  private:
  typedef char One;
  typedef struct
  {
    char a[2];
  } Two;
  template <typename C> static One test(int C::*);
  template <typename C> static Two test(...);
  public:
  enum
  {
    IsClass = sizeof(TypeTr<T>::template test<T>(0)) == 1
  };
};

/*====================================================================*/
/**  \name Raw address allocation policy
 *//*=================================================================*/
//@{

template <typename T, bool IsClass = TypeTr<T>::IsClass>
struct AllocRawPolicy;

///  If elements are class type, allocate using placement new
template <typename T>
struct AllocRawPolicy<T, true>
{
  static T *eval(void *const addr, const USz_t size)
    {
      return new(addr) T[size];
    }
};

///  If element are not class type, just assign the address
template <typename T>
struct AllocRawPolicy<T, false>
{
  static T *eval(void *const addr, const USz_t size)
    {
      return static_cast<T*>(addr);
    }
};
//@}

/*====================================================================*/
/**  \name Raw address release policy
 *//*=================================================================*/
//@{

template <typename T, bool IsClass = TypeTr<T>::IsClass>
struct ReleaseRawPolicy;

///  If elements are class type, invoke destructor on each element
template <typename T>
struct ReleaseRawPolicy<T, true>
{
  static void eval(T *p, USz_t size)
    {
      while (size--) p++->~T();
    }
};

///  If elements are not class type, do nothing
template <typename T>
struct ReleaseRawPolicy<T, false>
{
  static void eval(T *p, USz_t size)
  {
  }
};
//@}

/*====================================================================*/
///  Default allocator
/**  Allocators used by class array are simplified versions of the
 *   std::allocator and may contain member data.
 *   The following members are required:
 *   <ul>
 *     <li>  allocate(const USz_t size)
 *             - Allocate and construct 'size' objects
 *     <li>  deallocate(void *const addr, const USz_t size)
 *             - Deallocate and destroy 'size' objects at 'addr'
 *   </ul>
 *   The following members are optional:
 *   <ul>
 *     <li>  allocate(void *const addr, const USz_t size)
 *             - Construct 'size' objects at 'addr'.  If provided,
 *               CHArray::define(void *const addr, ...) can be used to
 *               construct on preallocated memory.  In this case,
 *               deallocate should only destroy the created objects,
 *               not delete the memory.
 *   </ul>
 *//*=================================================================*/

template <typename T>
class DefaultArrayAlloc
{
private:
  enum AllocBy
    {
      AllocUndefined,
      AllocNew,
      AllocRaw
    };
  AllocBy m_allocBy;

public:
  DefaultArrayAlloc()
    : m_allocBy(AllocUndefined)
    { }

  T* allocate(const USz_t size)
    {
      T *const p = new(std::nothrow) T[size];
      CH_assert(p != 0);
      m_allocBy = AllocNew;
      return p;
    }

  T* allocate(void *const addr, const USz_t size)
    {
      m_allocBy = AllocRaw;
      return AllocRawPolicy<T>::eval(addr, size);
    }

  void deallocate(T *p, const USz_t size)
    {
      switch (m_allocBy)
        {
        case AllocNew:
          delete[] p;
          break;
        case AllocRaw:
          ReleaseRawPolicy<T>::eval(p, size);
          break;
        case AllocUndefined:
          break;
        }
    }
};

/*====================================================================*/
/**  \name Data for the indexer classes.  This may vary in size
 *   depending on if the subscripts are zero-based or ranges.
 *//*=================================================================*/
//@{

/*--------------------------------------------------------------------*/
///  Any rank and a zero-based subscript
/*--------------------------------------------------------------------*/

template <unsigned Rank, typename DimT>
class IndexerData
{
protected:
  IndexerData()
    : m_ixStride(0)
    { }
  IndexerData(const USz_t stride, const DimT &dim)
    : m_ixStride(stride)
    { }
  void ixDefineLowerBound(const DimT &dim)
    { }
  IIx_t ixLowerBound() const
    {
      return 0;
    }
  bool ixValidBounds(const IIx_t i, const USz_t size) const
    {
      return (((USz_t)i) < size/m_ixStride);
    }
  IIx_t ixDimOffset() const
    {
      return 0;
    }
  USz_t m_ixStride;                   ///< Stride for this dimension
};

/*--------------------------------------------------------------------*/
///  Any rank and a subscript range
/*--------------------------------------------------------------------*/

template <unsigned Rank>
class IndexerData<Rank, CHRange>
{
protected:
  IndexerData()
    : m_ixStride(0), m_ixIB(0)
    { }
  IndexerData(const USz_t stride, const CHRange &dim)
    : m_ixStride(stride), m_ixIB(dim.dimB)
    { }
  void ixDefineLowerBound(const CHRange &dim)
    {
      m_ixIB = dim.dimB;
    }
  IIx_t ixLowerBound() const
    {
      return m_ixIB;
    }
  bool ixValidBounds(const IIx_t i, const USz_t size) const
    {
      return (i >= m_ixIB && i < m_ixIB + (IIx_t)(size/m_ixStride));
    }
  IIx_t ixDimOffset() const
    {
      return m_ixIB*m_ixStride;
    }
  USz_t m_ixStride;                   ///< Stride for this dimension
private:
  IIx_t m_ixIB;                       ///< Lower bound for this dimension
};

/*--------------------------------------------------------------------*/
///  Rank 1 and a zero-based subscript (no data)
/*--------------------------------------------------------------------*/

template <typename DimT>
class IndexerData<1, DimT>
{
protected:
  IndexerData()
    { }
  IndexerData(const DimT &dim)
    { }
  void ixDefineLowerBound(const DimT &dim)
    { }
  IIx_t ixLowerBound() const
    {
      return 0;
    }
  bool ixValidBounds(const IIx_t i, const USz_t size) const
    {
      return (((USz_t)i) < size);
    }
  void ixDefineTotalOffset(const IIx_t offset)
    { }
  IIx_t ixTotalOffset() const
    {
      return 0;
    }
};

/*--------------------------------------------------------------------*/
///  Rank 1 and a subscript range
/*--------------------------------------------------------------------*/

template <>
class IndexerData<1, CHRange>
{
protected:
  IndexerData()
    : m_ixIB(0)
    { }
  IndexerData(const CHRange &dim)
    : m_ixIB(dim.dimB)
    { }
  void ixDefineLowerBound(const CHRange &dim)
    {
      m_ixIB = dim.dimB;
    }
  IIx_t ixLowerBound() const
    {
      return m_ixIB;
    }
  bool ixValidBounds(const IIx_t i, const USz_t size) const
    {
      return (i >= m_ixIB && i < m_ixIB + (IIx_t)size);
    }
  void ixDefineTotalOffset(const IIx_t offset)
    {
      m_ixOffset = offset;
    }
  IIx_t ixTotalOffset() const
    {
      return m_ixOffset;
    }
private:
  IIx_t m_ixIB;                       ///< Lower bound for this dimension
  IIx_t m_ixOffset;                   ///< Total offset to an index because of
                                      ///< non-zero lower bounds for this and
                                      ///< all higher dimensions
};

/*====================================================================*/
/**  \name Provides indexing for row ordered (C) storage.
 *//*=================================================================*/
//@{

/*--------------------------------------------------------------------*/
///  Rank > 2 for zero-based subscripts and > 1 for subscript ranges
/*--------------------------------------------------------------------*/

template <unsigned Rank, typename DimT>
class RSIndexer : private IndexerData<Rank, DimT>
{
  typedef IndexerData<Rank, DimT> IxData;
public:
  RSIndexer()
    : IxData(), m_ixNext()
    { }

  RSIndexer(const IIx_t offset,
            const DimT &dim6, const DimT &dim5, const DimT &dim4,
            const DimT &dim3, const DimT &dim2, const DimT &dim1,
            const DimT &dim0)
    : IxData(dim5*dim4*dim3*dim2*dim1*dim0, dim6),
      m_ixNext(offset - IxData::ixDimOffset(),
               dim5, dim4, dim3, dim2, dim1, dim0)
    { }
  RSIndexer(const IIx_t offset,
            const DimT &dim5, const DimT &dim4, const DimT &dim3,
            const DimT &dim2, const DimT &dim1, const DimT &dim0)
    : IxData(dim4*dim3*dim2*dim1*dim0, dim5),
      m_ixNext(offset - IxData::ixDimOffset(), dim4, dim3, dim2, dim1, dim0)
    { }
  RSIndexer(const IIx_t offset,
            const DimT &dim4, const DimT &dim3, const DimT &dim2,
            const DimT &dim1, const DimT &dim0)
    : IxData(dim3*dim2*dim1*dim0, dim4),
      m_ixNext(offset - IxData::ixDimOffset(), dim3, dim2, dim1, dim0)
    { }
  RSIndexer(const IIx_t offset,
            const DimT &dim3, const DimT &dim2, const DimT &dim1,
            const DimT &dim0)
    : IxData(dim2*dim1*dim0, dim3),
      m_ixNext(offset - IxData::ixDimOffset(), dim2, dim1, dim0)
    { }
  RSIndexer(const IIx_t offset,
            const DimT &dim2, const DimT &dim1, const DimT &dim0)
    : IxData(dim1*dim0, dim2),
      m_ixNext(offset - IxData::ixDimOffset(), dim1, dim0)
    { }
  RSIndexer(const IIx_t offset,
            const DimT &dim1, const DimT &dim0)
    : IxData(sizeOfDim(dim0), dim1),
      m_ixNext(offset - IxData::ixDimOffset(), dim0)
    { }

  void ixDefine(const IIx_t offset,
                const DimT &dim6, const DimT &dim5, const DimT &dim4,
                const DimT &dim3, const DimT &dim2, const DimT &dim1,
                const DimT &dim0)
    {
      IxData::m_ixStride = dim5*dim4*dim3*dim2*dim1*dim0;
      IxData::ixDefineLowerBound(dim6);
      m_ixNext.ixDefine(offset - IxData::ixDimOffset(),
                        dim5, dim4, dim3, dim2, dim1, dim0);
    }
  void ixDefine(const IIx_t offset,
                const DimT &dim5, const DimT &dim4, const DimT &dim3,
                const DimT &dim2, const DimT &dim1, const DimT &dim0)
    {
      IxData::m_ixStride = dim4*dim3*dim2*dim1*dim0;
      IxData::ixDefineLowerBound(dim5);
      m_ixNext.ixDefine(offset - IxData::ixDimOffset(),
                        dim4, dim3, dim2, dim1, dim0);
    }
  void ixDefine(const IIx_t offset,
                const DimT &dim4, const DimT &dim3, const DimT &dim2,
                const DimT &dim1, const DimT &dim0)
    {
      IxData::m_ixStride = dim3*dim2*dim1*dim0;
      IxData::ixDefineLowerBound(dim4);
      m_ixNext.ixDefine(offset - IxData::ixDimOffset(), dim3, dim2, dim1, dim0);
    }
  void ixDefine(const IIx_t offset,
                const DimT &dim3, const DimT &dim2, const DimT &dim1,
                const DimT &dim0)
    {
      IxData::m_ixStride = dim2*dim1*dim0;
      IxData::ixDefineLowerBound(dim3);
      m_ixNext.ixDefine(offset - IxData::ixDimOffset(), dim2, dim1, dim0);
    }
  void ixDefine(const IIx_t offset,
                const DimT &dim2, const DimT &dim1, const DimT &dim0)
    {
      IxData::m_ixStride = dim1*dim0;
      IxData::ixDefineLowerBound(dim2);
      m_ixNext.ixDefine(offset - IxData::ixDimOffset(), dim1, dim0);
    }
  void ixDefine(const IIx_t offset,
                const DimT &dim1, const DimT &dim0)
    {
      IxData::m_ixStride = sizeOfDim(dim0);
      IxData::ixDefineLowerBound(dim1);
      m_ixNext.ixDefine(offset - IxData::ixDimOffset(), dim0);
    }

  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i6, const IIx_t i5, const IIx_t i4,
                  const IIx_t i3, const IIx_t i2, const IIx_t i1,
                  const IIx_t i0) const
    {
      CH_assert((IxData::ixValidBounds(i6, size)));
      return IxData::m_ixStride*i6 +
        m_ixNext.ixIndex1D(IxData::m_ixStride, i5, i4, i3, i2, i1, i0);
    }
  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i5, const IIx_t i4, const IIx_t i3,
                  const IIx_t i2, const IIx_t i1, const IIx_t i0) const
    {
      CH_assert((IxData::ixValidBounds(i5, size)));
      return IxData::m_ixStride*i5 +
        m_ixNext.ixIndex1D(IxData::m_ixStride, i4, i3, i2, i1, i0);
    }
  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i4, const IIx_t i3, const IIx_t i2,
                  const IIx_t i1, const IIx_t i0) const
    {
      CH_assert((IxData::ixValidBounds(i4, size)));
      return IxData::m_ixStride*i4 +
        m_ixNext.ixIndex1D(IxData::m_ixStride, i3, i2, i1, i0);
    }
  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i3, const IIx_t i2, const IIx_t i1,
                  const IIx_t i0) const
    {
      CH_assert((IxData::ixValidBounds(i3, size)));
      return IxData::m_ixStride*i3 +
        m_ixNext.ixIndex1D(IxData::m_ixStride, i2, i1, i0);
    }
  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i2, const IIx_t i1, const IIx_t i0) const
    {
      CH_assert((IxData::ixValidBounds(i2, size)));
      return IxData::m_ixStride*i2 +
        m_ixNext.ixIndex1D(IxData::m_ixStride, i1, i0);
    }
  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i1, const IIx_t i0) const
    {
      CH_assert((IxData::ixValidBounds(i1, size)));
      return IxData::m_ixStride*i1 +
        m_ixNext.ixIndex1D(IxData::m_ixStride, i0);
    }

  USz_t ixDimSize(const USz_t size, const unsigned dim) const
    {
      return (dim == Rank) ?
        size/IxData::m_ixStride :
        m_ixNext.ixDimSize(IxData::m_ixStride, dim);
    }
  IIx_t ixLowerBound(const unsigned dim) const
    {
      return (dim == Rank) ?
        IxData::ixLowerBound() :
        m_ixNext.ixLowerBound(dim);
    }
  IIx_t ixUpperBound(const USz_t size, const unsigned dim) const
    {
      return (dim == Rank) ?
        IxData::ixLowerBound() + size/IxData::m_ixStride - 1 :
        m_ixNext.ixUpperBound(IxData::m_ixStride, dim);
    }

  unsigned ixIV(const unsigned i) const
    {
      return SpaceDim - 1 - i;
    }

private:
  RSIndexer<Rank-1, DimT> m_ixNext;   ///< Next lower rank
};

/*--------------------------------------------------------------------*/
///  Specialization for Rank 1 (any type of subscript)
/*--------------------------------------------------------------------*/

template <typename DimT>
class RSIndexer<1, DimT> : private IndexerData<1, DimT>
{
  typedef IndexerData<1, DimT> IxData;
public:
  RSIndexer()
    : IxData()
    { }

  RSIndexer(const IIx_t offset,
            const DimT &dim0)
    : IxData(dim0)
    {
      IxData::ixDefineTotalOffset(offset - IxData::ixLowerBound());
    }

  void ixDefine(const IIx_t offset,
                const DimT &dim0)
    {
      IxData::ixDefineLowerBound(dim0);
      IxData::ixDefineTotalOffset(offset - IxData::ixLowerBound());
    }

  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i0) const
    {
      CH_assert((IxData::ixValidBounds(i0, size)));
      return i0 + IxData::ixTotalOffset();
    }

  USz_t ixDimSize(const USz_t size, const unsigned dim) const
    {
      return size;
    }
  IIx_t ixLowerBound(const unsigned dim) const
    {
      return IxData::ixLowerBound();
    }
  IIx_t ixUpperBound(const USz_t size, const unsigned dim) const
    {
      return IxData::ixLowerBound() + size - 1;
    }

  unsigned ixIV(const unsigned i) const
    {
      return SpaceDim - 1 - i;
    }
};

/*--------------------------------------------------------------------*/
///  Full specialization for Rank 2 and zero-based subscripts
/**  Saves the byte required to store the rank1 indexer (which itself
 *   holds no data if DimT=USz_t)
 *//*-----------------------------------------------------------------*/

template <>
class RSIndexer<2, USz_t> : private IndexerData<2, USz_t>
{
  typedef IndexerData<2, USz_t> IxData;
public:
  RSIndexer()
    : IxData()
    { }

  RSIndexer(const IIx_t offset,
            const USz_t &dim1, const USz_t &dim0)
    : IxData(dim0, dim1)
    { }

  void ixDefine(const IIx_t offset,
                const USz_t &dim1, const USz_t &dim0)
    {
      IxData::m_ixStride = dim0;
    }

  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i1, const IIx_t i0) const
    {
      CH_assert((IxData::ixValidBounds(i1, size)));
      CH_assert((i0 < IxData::m_ixStride));
      return IxData::m_ixStride*i1 + i0;
    }

  USz_t ixDimSize(const USz_t size, const unsigned dim) const
    {
      return (dim == 2) ? size/IxData::m_ixStride :
                          IxData::m_ixStride;
    }
  IIx_t ixLowerBound(const unsigned dim) const
    {
      return 0;
    }
  IIx_t ixUpperBound(const USz_t size, const unsigned dim) const
    {
      return ixDimSize(size, dim) - 1;
    }

  unsigned ixIV(const unsigned i) const
    {
      return SpaceDim - 1 - i;
    }
};
//@}

/*====================================================================*/
/**  \name Provides indexing for column ordered (Fortran) storage.
 *//*=================================================================*/
//@{

/*--------------------------------------------------------------------*/
///  Rank > 2 for zero-based subscripts and > 1 for subscript ranges
/*--------------------------------------------------------------------*/

template <unsigned Rank, typename DimT>
class CSIndexer : private IndexerData<Rank, DimT>
{
  typedef IndexerData<Rank, DimT> IxData;
public:
  CSIndexer()
    : IxData(), m_ixNext()
    { }

  CSIndexer(const IIx_t offset,
            const DimT &dim0, const DimT &dim1, const DimT &dim2,
            const DimT &dim3, const DimT &dim4, const DimT &dim5,
            const DimT &dim6)
    : IxData(dim0*dim1*dim2*dim3*dim4*dim5, dim6),
      m_ixNext(offset - IxData::ixDimOffset(),
               dim0, dim1, dim2, dim3, dim4, dim5)
    { }
  CSIndexer(const IIx_t offset,
            const DimT &dim0, const DimT &dim1, const DimT &dim2,
            const DimT &dim3, const DimT &dim4, const DimT &dim5)
    : IxData(dim0*dim1*dim2*dim3*dim4, dim5),
      m_ixNext(offset - IxData::ixDimOffset(), dim0, dim1, dim2, dim3, dim4)
    { }
  CSIndexer(const IIx_t offset,
            const DimT &dim0, const DimT &dim1, const DimT &dim2,
            const DimT &dim3, const DimT &dim4)
    : IxData(dim0*dim1*dim2*dim3, dim4),
      m_ixNext(offset - IxData::ixDimOffset(), dim0, dim1, dim2, dim3)
    { }
  CSIndexer(const IIx_t offset,
            const DimT &dim0, const DimT &dim1, const DimT &dim2,
            const DimT &dim3)
    : IxData(dim0*dim1*dim2, dim3),
      m_ixNext(offset - IxData::ixDimOffset(), dim0, dim1, dim2)
    { }
  CSIndexer(const IIx_t offset,
            const DimT &dim0, const DimT &dim1, const DimT &dim2)
    : IxData(dim0*dim1, dim2),
      m_ixNext(offset - IxData::ixDimOffset(), dim0, dim1)
    { }
  CSIndexer(const IIx_t offset,
            const DimT &dim0, const DimT &dim1)
    : IxData(sizeOfDim(dim0), dim1),
      m_ixNext(offset - IxData::ixDimOffset(), dim0)
    { }

  void ixDefine(const IIx_t offset,
                const DimT &dim0, const DimT &dim1, const DimT &dim2,
                const DimT &dim3, const DimT &dim4, const DimT &dim5,
                const DimT &dim6)
    {
      IxData::m_ixStride = dim0*dim1*dim2*dim3*dim4*dim5;
      IxData::ixDefineLowerBound(dim6);
      m_ixNext.ixDefine(offset - IxData::ixDimOffset(),
                        dim0, dim1, dim2, dim3, dim4, dim5);
    }
  void ixDefine(const IIx_t offset,
                const DimT &dim0, const DimT &dim1, const DimT &dim2,
                const DimT &dim3, const DimT &dim4, const DimT &dim5)
    {
      IxData::m_ixStride = dim0*dim1*dim2*dim3*dim4;
      IxData::ixDefineLowerBound(dim5);
      m_ixNext.ixDefine(offset - IxData::ixDimOffset(),
                        dim0, dim1, dim2, dim3, dim4);
    }
  void ixDefine(const IIx_t offset,
                const DimT &dim0, const DimT &dim1, const DimT &dim2,
                const DimT &dim3, const DimT &dim4)
    {
      IxData::m_ixStride = dim0*dim1*dim2*dim3;
      IxData::ixDefineLowerBound(dim4);
      m_ixNext.ixDefine(offset - IxData::ixDimOffset(), dim0, dim1, dim2, dim3);
    }
  void ixDefine(const IIx_t offset,
                const DimT &dim0, const DimT &dim1, const DimT &dim2,
                const DimT &dim3)
    {
      IxData::m_ixStride = dim0*dim1*dim2;
      IxData::ixDefineLowerBound(dim3);
      m_ixNext.ixDefine(offset - IxData::ixDimOffset(), dim0, dim1, dim2);
    }
  void ixDefine(const IIx_t offset,
                const DimT &dim0, const DimT &dim1, const DimT &dim2)
    {
      IxData::m_ixStride = dim0*dim1;
      IxData::ixDefineLowerBound(dim2);
      m_ixNext.ixDefine(offset - IxData::ixDimOffset(), dim0, dim1);
    }
  void ixDefine(const IIx_t offset,
                const DimT &dim0, const DimT &dim1)
    {
      IxData::m_ixStride = sizeOfDim(dim0);
      IxData::ixDefineLowerBound(dim1);
      m_ixNext.ixDefine(offset - IxData::ixDimOffset(), dim0);
    }

  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i0, const IIx_t i1, const IIx_t i2,
                  const IIx_t i3, const IIx_t i4, const IIx_t i5,
                  const IIx_t i6) const
    {
      CH_assert((IxData::ixValidBounds(i6, size)));
      return IxData::m_ixStride*i6 +
        m_ixNext.ixIndex1D(IxData::m_ixStride, i0, i1, i2, i3, i4, i5);
    }
  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i0, const IIx_t i1, const IIx_t i2,
                  const IIx_t i3, const IIx_t i4, const IIx_t i5) const
    {
      CH_assert((IxData::ixValidBounds(i5, size)));
      return IxData::m_ixStride*i5 +
        m_ixNext.ixIndex1D(IxData::m_ixStride, i0, i1, i2, i3, i4);
    }
  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i0, const IIx_t i1, const IIx_t i2,
                  const IIx_t i3, const IIx_t i4) const
    {
      CH_assert((IxData::ixValidBounds(i4, size)));
      return IxData::m_ixStride*i4 +
        m_ixNext.ixIndex1D(IxData::m_ixStride, i0, i1, i2, i3);
    }
  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i0, const IIx_t i1, const IIx_t i2,
                  const IIx_t i3) const
    {
      CH_assert((IxData::ixValidBounds(i3, size)));
      return IxData::m_ixStride*i3 +
        m_ixNext.ixIndex1D(IxData::m_ixStride, i0, i1, i2);
    }
  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i0, const IIx_t i1, const IIx_t i2) const
    {
      CH_assert((IxData::ixValidBounds(i2, size)));
      return IxData::m_ixStride*i2 +
        m_ixNext.ixIndex1D(IxData::m_ixStride, i0, i1);
    }
  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i0, const IIx_t i1) const
    {
      CH_assert((IxData::ixValidBounds(i1, size)));
      return IxData::m_ixStride*i1 +
        m_ixNext.ixIndex1D(IxData::m_ixStride, i0);
    }

  USz_t ixDimSize(const USz_t size, const unsigned dim) const
    {
      return (dim == Rank) ?
        size/IxData::m_ixStride :
        m_ixNext.ixDimSize(IxData::m_ixStride, dim);
    }
  IIx_t ixLowerBound(const unsigned dim) const
    {
      return (dim == Rank) ?
        IxData::ixLowerBound() :
        m_ixNext.ixLowerBound(dim);
    }
  IIx_t ixUpperBound(const USz_t size, const unsigned dim) const
    {
      return (dim == Rank) ?
        IxData::ixLowerBound() + size/IxData::m_ixStride - 1 :
        m_ixNext.ixUpperBound(IxData::m_ixStride, dim);
    }

  unsigned ixIV(const unsigned i) const
    {
      return i;
    }

private:
  CSIndexer<Rank-1, DimT> m_ixNext;   ///< Next lower rank
};

/*--------------------------------------------------------------------*/
///  Specialization for Rank 1 (any type of subscript)
/*--------------------------------------------------------------------*/

template <typename DimT>
class CSIndexer<1, DimT> : private IndexerData<1, DimT>
{
  typedef IndexerData<1, DimT> IxData;
public:
  CSIndexer()
    : IxData()
    { }

  CSIndexer(const IIx_t offset,
            const DimT &dim0)
    : IxData(dim0)
    {
      IxData::ixDefineTotalOffset(offset - IxData::ixLowerBound());
    }

  void ixDefine(const IIx_t offset, const DimT &dim0)
    {
      IxData::ixDefineLowerBound(dim0);
      IxData::ixDefineTotalOffset(offset - IxData::ixLowerBound());
    }

  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i0) const
    {
      CH_assert((IxData::ixValidBounds(i0, size)));
      return i0 + IxData::ixTotalOffset();
    }

  USz_t ixDimSize(const USz_t size, const unsigned dim) const
    {
      return size;
    }
  IIx_t ixLowerBound(const unsigned dim) const
    {
      return IxData::ixLowerBound();
    }
  IIx_t ixUpperBound(const USz_t size, const unsigned dim) const
    {
      return IxData::ixLowerBound() + size - 1;
    }

  unsigned ixIV(const unsigned i) const
    {
      return i;
    }
};

/*--------------------------------------------------------------------*/
///  Full specialization for Rank 2 and zero-based subscripts
/**  Saves the byte required to store the rank1 indexer (which itself
 *   holds no data if DimT=USz_t)
 *//*-----------------------------------------------------------------*/

template <>
class CSIndexer<2, USz_t> : private IndexerData<2, USz_t>
{
  typedef IndexerData<2, USz_t> IxData;
public:
  CSIndexer()
    : IxData()
    { }

  CSIndexer(const IIx_t offset,
            const USz_t dim0, const USz_t dim1)
    : IxData(dim0, dim1)
    { }

  void ixDefine(const IIx_t offset,
                const USz_t dim0, const USz_t dim1)
    {
      IxData::m_ixStride = dim0;
    }

  USz_t ixIndex1D(const USz_t size,
                  const IIx_t i0, const IIx_t i1) const
    {
      CH_assert((IxData::ixValidBounds(i1, size)));
      CH_assert((i0 < IxData::m_ixStride));
      return IxData::m_ixStride*i1 + i0;
    }

  IIx_t ixDimSize(const USz_t size, const unsigned dim) const
    {
      return (dim == 2) ? size/IxData::m_ixStride :
                          IxData::m_ixStride;
    }
  IIx_t ixLowerBound(const unsigned dim) const
    {
      return 0;
    }
  IIx_t ixUpperBound(const USz_t size, const unsigned dim) const
    {
      return ixDimSize(size, dim) - 1;
    }

  unsigned ixIV(const unsigned i) const
    {
      return i;
    }
};
//@}

/*====================================================================*/
/**  \name Array configuration.  Select types based on 'ArConf'
 *//*=================================================================*/
//@{

template <unsigned Rank, int ArConf>
struct ArTr;
template <unsigned Rank>
struct ArTr<Rank, ArZeroRow>
{
  typedef ArSp::USz_t DimT;
  typedef ArSp::RSIndexer<Rank, DimT> IndexerT;
};
template <unsigned Rank>
struct ArTr<Rank, ArZeroCol>
{
  typedef ArSp::USz_t DimT;
  typedef ArSp::CSIndexer<Rank, DimT> IndexerT;
};
template <unsigned Rank>
struct ArTr<Rank, ArRangeRow>
{
  typedef CHRange DimT;
  typedef ArSp::RSIndexer<Rank, DimT> IndexerT;
  typedef int BoxAllowed;
};
template <unsigned Rank>
struct ArTr<Rank, ArRangeCol>
{
  typedef CHRange DimT;
  typedef ArSp::CSIndexer<Rank, DimT> IndexerT;
  typedef int BoxAllowed;
};
//@}

}  // End of namespace ArSp

/*******************************************************************************
 */
///  Multidimensional array class
/**
 *   \tparam T          Type of element
 *   \tparam Rank       Rank of the array
 *   \tparam ArConf     Configuration of the indexing. Allowable values are:
 *                      ArZeroRow  - Row-ordered indexing with lower bound = 0
 *                      ArZeroCol  - Column-ordered indexing with lower
 *                                   bound = 0 (DEFAULT)
 *                      ArRangeRow - Row-ordered indexing with subscript ranges
 *                      ArRangeCol - Column-ordered indexing with subscript
 *                                   ranges
 *   \tparam Alloc      Allocator for the array
 *
 *   \note
 *   <ul>
 *     <li> Default template arguments specify column-ordered arrays since it
 *          is expected that this class will be mostly used as storage for
 *          data used in Fortran routines.
 *     <li> When an assertion fails, note that the indexers label the contiguous
 *          dimension as i0 or dim0 and the others from there.  This class
 *          labels all dimensions assuming row-ordered (C) storage.
 *     <li> Currently supports up to rank 7, matching the maximum permitted by
 *          the Fortran 2003 standard (the 2008 standard extends this to 15 and
 *          this class will be extended once that is commonplace).  With
 *          IntVects or Boxes used to specify some of the dimensions, 2 extra
 *          dimensions are supported, allowing a matrix in each cell.  This
 *          implies a rank of SpaceDim+2 so you can only add 1 extra dimension
 *          if you want to support SpaceDim=6.
 *     <li> Empty base-class optimizations are frequently used.  On 32-bit arch,
 *          sizeof(CHArray<int, 1, ArZeroRow, ArSp::NewArrayAlloc<int> >) = 8
 *          sizeof(CHArray<int, 2, ArZeroRow, ArSp::NewArrayAlloc<int> >) = 12
 *     <li> There is support for passing these arrays to Fortran -- see Chombo
 *          design doc
 *     <li> CHMatrix will define a real, rank 2, column-ordered matrix
 *     <li> You can also allocate an array of matrices, continguous in memory,
 *          using the ArrayOfMatrixAlloc allocator.
 *   </ul>
 *
 *   Example                                                           \verbatim
 *
 *     CHArray<int, 3> A(2, 2, 2);
 *     A = 2;
 *     A(0, 0, 0) = 4;
 *     std::cout << A << std::endl;
 *     CHArray<double, 2> B;
 *     B.define(3, 4);
 *     B(2, 3) = 5.5;
 *                                                                  \endverbatim
 *   One can also use subscript ranges but the array has to be configured to
 *   support this.  Ranges are specified with CHRange in place of a size.  Any
 *   integer types are converted to CHRange(0, int-1).  E.g.,
 *                                                                     \verbatim
 *     CHArray<int, 3, ArRangeRow> C(CHRange(-1,0), 2, 2);
 *     C = 2;
 *     C(-1, 0, 0) = 4;
 *                                                                  \endverbatim
 *   IntVects are allowed for specifying some of the dimensions so that this
 *   class can be used with boxes, e.g:                                \verbatim
 *
 *     CHArray<int, SpaceDim+1> A(box.size(), 1)
 *                                                                  \endverbatim
 *   The indices of the IntVect are rearranged so that the first index is always
 *   closest to contiguous memory no matter if the actual memory layout is row-
 *   ordered or column-ordered.  In other words, the first index of the IntVect
 *   always has the lowest stride.  As such, the allocated memory can be
 *   accessed by iterating over the box in Chombo Fortran.  Whether the
 *   remaining indices (or components) are left or right of the IntVect and the
 *   storage order determine if the components are distributed across the
 *   IntVects (like with a BaseFab) or are contiguous in each IntVect. The
 *   effect is natural; for row storage, components on the right will lead to
 *   contigous memory.  E.g.,
 *                                                                     \verbatim
 *     typedef CHArray<int, SpaceDim+1, ArZeroRow> CHArrayRS;
 *     typedef CHArray<int, SpaceDim+1, ArZeroCol> CHArrayCS;
 *     const IntVect iv(2, 2);
 *     CHArrayRS ARS_IVComp(iv, 2);
 *     CHArrayRS ARS_CompIV(2, iv);
 *     CHArrayCS ACS_IVComp(iv, 2);
 *     CHArrayCS ACS_CompIV(2, iv);
 *     int ic = 1;
 *     for (int j = 0; j != 2; ++j)
 *       for (int i = 0; i != 2; ++i)  // Column storage!
 *         {
 *           const IntVect ixv(i, j);
 *           ARS_IVComp(ixv, 0) = ic;
 *           ARS_IVComp(ixv, 1) = -ic;
 *           ARS_CompIV(0, ixv) = ic;
 *           ARS_CompIV(1, ixv) = -ic;
 *           ACS_IVComp(ixv, 0) = ic;
 *           ACS_IVComp(ixv, 1) = -ic;
 *           ACS_CompIV(0, ixv) = ic;
 *           ACS_CompIV(1, ixv) = -ic;
 *           ++ic;
 *         }
 *
 *     std::cout << "ARS_IVComp: " << ARS_IVComp << std::endl;
 *     std::cout << "ARS_CompIV: " << ARS_CompIV << std::endl;
 *     std::cout << "ACS_IVComp: " << ACS_IVComp << std::endl;
 *     std::cout << "ACS_CompIV: " << ACS_CompIV << std::endl;
 *                                                                  \endverbatim
 *   will produce
 *                                                                     \verbatim
 *     $ ARS_IVComp: 1 -1 2 -2 3 -3 4 -4
 *     $ ARS_CompIV: 1 2 3 4 -1 -2 -3 -4
 *     $ ACS_IVComp: 1 2 3 4 -1 -2 -3 -4
 *     $ ACS_CompIV: 1 -1 2 -2 3 -3 4 -4
 *                                                                  \endverbatim
 *   Note that the IntVects simply expand to dimensions and there is nothing to
 *   ensure that you access in the array in the same manner you defined it
 *   (i.e., switching the order of the IntVect and the component)!
 *
 *   IntVects are always zero-based and there is no way to specify a range of
 *   IntVects.  However, one can use Boxes instead of IntVects to define the
 *   array but only if the array is configured with subscript ranges
 *   (ArConf = ArRangeRow or ArRangeCol)
 *                                                                     \verbatim
 *     Box box(IntVect(-1,-1), IntVect(1, 1));
 *     CHArray<int, SpaceDim+1, ArRangeRow> AB(box, 2);
 *     AB = 1;
 *     AB(IntVect(-1, -1), 0) = 2;
 *     std::cout << AB << std::endl;
 *                                                                  \endverbatim
 *
 ******************************************************************************/

template <typename T,
          unsigned Rank,
          int ArConf,      // Default ArZeroCol
          typename Alloc>  // Default ArSp::DefaultArrayAlloc<T>
class CHArray : private ArSp::ArTr<Rank, ArConf>::IndexerT
{
  typedef typename ArSp::ArTr<Rank, ArConf>::DimT DimT;
  typedef typename ArSp::ArTr<Rank, ArConf>::IndexerT Indexer;
public:

  ///  Default constructor
  CHArray()
    : Indexer(), m_arrayImpl()
    { }

  ///  Construct with dimensions
  CHArray(const DimT &dim6, const DimT &dim5, const DimT &dim4,
          const DimT &dim3, const DimT &dim2, const DimT &dim1,
          const DimT &dim0)
    : Indexer(0, dim6, dim5, dim4, dim3, dim2, dim1, dim0),
      m_arrayImpl(dim6*dim5*dim4*dim3*dim2*dim1*dim0)
    { }
  CHArray(const DimT &dim5, const DimT &dim4, const DimT &dim3,
          const DimT &dim2, const DimT &dim1, const DimT &dim0)
    : Indexer(0, dim5, dim4, dim3, dim2, dim1, dim0),
      m_arrayImpl(dim5*dim4*dim3*dim2*dim1*dim0)
    { }
  CHArray(const DimT &dim4, const DimT &dim3, const DimT &dim2,
          const DimT &dim1, const DimT &dim0)
    : Indexer(0, dim4, dim3, dim2, dim1, dim0),
      m_arrayImpl(dim4*dim3*dim2*dim1*dim0)
    { }
  CHArray(const DimT &dim3, const DimT &dim2, const DimT &dim1,
          const DimT &dim0)
    : Indexer(0, dim3, dim2, dim1, dim0), m_arrayImpl(dim3*dim2*dim1*dim0)
    { }
  CHArray(const DimT &dim2, const DimT &dim1, const DimT &dim0)
    : Indexer(0, dim2, dim1, dim0), m_arrayImpl(dim2*dim1*dim0)
    { }
  CHArray(const DimT &dim1, const DimT &dim0)
    : Indexer(0, dim1, dim0), m_arrayImpl(dim1*dim0)
    { }
  CHArray(const DimT &dim0)
    : Indexer(0, dim0), m_arrayImpl(ArSp::sizeOfDim(dim0))
    { }
  // With IntVect
  CHArray(const IntVect &iv, const DimT &dimc1, const DimT &dimc0)
    : Indexer(0, D_DEFIV, dimc1, dimc0), m_arrayImpl(iv.product()*dimc1*dimc0)
    { }
  CHArray(const IntVect &iv, const DimT &dimc0)
    : Indexer(0, D_DEFIV, dimc0), m_arrayImpl(iv.product()*dimc0)
    { }
  CHArray(const DimT &dimc1, const DimT &dimc0, const IntVect &iv)
    : Indexer(0, dimc1, dimc0, D_DEFIV), m_arrayImpl(dimc1*dimc0*iv.product())
    { }
  CHArray(const DimT &dimc0, const IntVect &iv)
    : Indexer(0, dimc0, D_DEFIV), m_arrayImpl(dimc0*iv.product())
    { }
  // With Box
  CHArray(const Box &box, const DimT &dimc1, const DimT &dimc0)
    : Indexer(0, D_DEFBOX, dimc1, dimc0),
      m_arrayImpl(box.size().product()*dimc1*dimc0)
    {
      // Any error with "BoxAllowed" means you are not allowed to use boxes with
      // zero-based lower bounds.  Use ArConf=ArRangeRow or ArConf=ArRangeCol.
      typedef typename ArSp::ArTr<Rank, ArConf>::BoxAllowed Check;
    }
  CHArray(const Box &box, const DimT &dimc0)
    : Indexer(0, D_DEFBOX, dimc0), m_arrayImpl(box.size().product()*dimc0)
    {
      // Any error with "BoxAllowed" means you are not allowed to use boxes with
      // zero-based lower bounds.  Use ArConf=ArRangeRow or ArConf=ArRangeCol.
      typedef typename ArSp::ArTr<Rank, ArConf>::BoxAllowed Check;
    }
  CHArray(const DimT &dimc1, const DimT &dimc0, const Box &box)
    : Indexer(0, dimc1, dimc0, D_DEFBOX),
      m_arrayImpl(dimc1*dimc0*box.size().product())
    {
      // Any error with "BoxAllowed" means you are not allowed to use boxes with
      // zero-based lower bounds.  Use ArConf=ArRangeRow or ArConf=ArRangeCol.
      typedef typename ArSp::ArTr<Rank, ArConf>::BoxAllowed Check;
    }
  CHArray(const DimT &dimc0, const Box &box)
    : Indexer(0, dimc0, D_DEFBOX), m_arrayImpl(dimc0*box.size().product())
    {
      // Any error with "BoxAllowed" means you are not allowed to use boxes with
      // zero-based lower bounds.  Use ArConf=ArRangeRow or ArConf=ArRangeCol.
      typedef typename ArSp::ArTr<Rank, ArConf>::BoxAllowed Check;
    }

  ///  Define the dimensions
  void define(const DimT &dim6, const DimT &dim5, const DimT &dim4,
              const DimT &dim3, const DimT &dim2, const DimT &dim1,
              const DimT &dim0)
    {
      Indexer::ixDefine(0, dim6, dim5, dim4, dim3, dim2, dim1, dim0);
      m_arrayImpl.define(dim6*dim5*dim4*dim3*dim2*dim1*dim0);
    }
  void define(const DimT &dim5, const DimT &dim4, const DimT &dim3,
              const DimT &dim2, const DimT &dim1, const DimT &dim0)
    {
      Indexer::ixDefine(0, dim5, dim4, dim3, dim2, dim1, dim0);
      m_arrayImpl.define(dim5*dim4*dim3*dim2*dim1*dim0);
    }
  void define(const DimT &dim4, const DimT &dim3, const DimT &dim2,
              const DimT &dim1, const DimT &dim0)
    {
      Indexer::ixDefine(0, dim4, dim3, dim2, dim1, dim0);
      m_arrayImpl.define(dim4*dim3*dim2*dim1*dim0);
    }
  void define(const DimT &dim3, const DimT &dim2, const DimT &dim1,
              const DimT &dim0)
    {
      Indexer::ixDefine(0, dim3, dim2, dim1, dim0);
      m_arrayImpl.define(dim3*dim2*dim1*dim0);
    }
  void define(const DimT &dim2, const DimT &dim1, const DimT &dim0)
    {
      Indexer::ixDefine(0, dim2, dim1, dim0);
      m_arrayImpl.define(dim2*dim1*dim0);
    }
  void define(const DimT &dim1, const DimT &dim0)
    {
      Indexer::ixDefine(0, dim1, dim0);
      m_arrayImpl.define(dim1*dim0);
    }
  void define(const DimT &dim0)
    {
      Indexer::ixDefine(0, dim0);
      m_arrayImpl.define(ArSp::sizeOfDim(dim0));
    }
  // With IntVect
  void define(const IntVect &iv, const DimT &dimc1, const DimT &dimc0)
    {
      Indexer::ixDefine(0, D_DEFIV, dimc1, dimc0);
      m_arrayImpl.define(iv.product()*dimc1*dimc0);
    }
  void define(const IntVect &iv, const DimT &dimc0)
    {
      Indexer::ixDefine(0, D_DEFIV, dimc0);
      m_arrayImpl.define(iv.product()*dimc0);
    }
  void define(const DimT &dimc1, const DimT &dimc0, const IntVect &iv)
    {
      Indexer::ixDefine(0, dimc1, dimc0, D_DEFIV);
      m_arrayImpl.define(dimc1*dimc0*iv.product());
    }
  void define(const DimT &dimc0, const IntVect &iv)
    {
      Indexer::ixDefine(0, dimc0, D_DEFIV);
      m_arrayImpl.define(dimc0*iv.product());
    }
  // With Box
  void define(const Box &box, const DimT &dimc1, const DimT &dimc0)
    {
      // Any error with "BoxAllowed" means you are not allowed to use boxes with
      // zero-based lower bounds.  Use ArConf=ArRangeRow or ArConf=ArRangeCol.
      typedef typename ArSp::ArTr<Rank, ArConf>::BoxAllowed Check;
      Indexer::ixDefine(0, D_DEFBOX, dimc1, dimc0);
      m_arrayImpl.define(box.size().product()*dimc1*dimc0);
    }
  void define(const Box &box, const DimT &dimc0)
    {
      // Any error with "BoxAllowed" means you are not allowed to use boxes with
      // zero-based lower bounds.  Use ArConf=ArRangeRow or ArConf=ArRangeCol.
      typedef typename ArSp::ArTr<Rank, ArConf>::BoxAllowed Check;
      Indexer::ixDefine(0, D_DEFBOX, dimc0);
      m_arrayImpl.define(box.size().product()*dimc0);
    }
  void define(const DimT &dimc1, const DimT &dimc0, const Box &box)
    {
      // Any error with "BoxAllowed" means you are not allowed to use boxes with
      // zero-based lower bounds.  Use ArConf=ArRangeRow or ArConf=ArRangeCol.
      typedef typename ArSp::ArTr<Rank, ArConf>::BoxAllowed Check;
      Indexer::ixDefine(0, dimc1, dimc0, D_DEFBOX);
      m_arrayImpl.define(dimc1*dimc0*box.size().product());
    }
  void define(const DimT &dimc0, const Box &box)
    {
      // Any error with "BoxAllowed" means you are not allowed to use boxes with
      // zero-based lower bounds.  Use ArConf=ArRangeRow or ArConf=ArRangeCol.
      typedef typename ArSp::ArTr<Rank, ArConf>::BoxAllowed Check;
      Indexer::ixDefine(0, dimc0, D_DEFBOX);
      m_arrayImpl.define(dimc0*box.size().product());
    }

  ///  Define the dimensions and allocate on 'addr'
  void define(void *const addr,
              const DimT &dim6, const DimT &dim5, const DimT &dim4,
              const DimT &dim3, const DimT &dim2, const DimT &dim1,
              const DimT &dim0)
    {
      Indexer::ixDefine(0, dim6, dim5, dim4, dim3, dim2, dim1, dim0);
      m_arrayImpl.define(addr, dim6*dim5*dim4*dim3*dim2*dim1*dim0);
    }
  void define(void *const addr,
              const DimT &dim5, const DimT &dim4, const DimT &dim3,
              const DimT &dim2, const DimT &dim1, const DimT &dim0)
    {
      Indexer::ixDefine(0, dim5, dim4, dim3, dim2, dim1, dim0);
      m_arrayImpl.define(addr, dim5*dim4*dim3*dim2*dim1*dim0);
    }
  void define(void *const addr,
              const DimT &dim4, const DimT &dim3, const DimT &dim2,
              const DimT &dim1, const DimT &dim0)
    {
      Indexer::ixDefine(0, dim4, dim3, dim2, dim1, dim0);
      m_arrayImpl.define(addr, dim4*dim3*dim2*dim1*dim0);
    }
  void define(void *const addr,
              const DimT &dim3, const DimT &dim2, const DimT &dim1,
              const DimT &dim0)
    {
      Indexer::ixDefine(0, dim3, dim2, dim1, dim0);
      m_arrayImpl.define(addr, dim3*dim2*dim1*dim0);
    }
  void define(void *const addr,
              const DimT &dim2, const DimT &dim1, const DimT &dim0)
    {
      Indexer::ixDefine(0, dim2, dim1, dim0);
      m_arrayImpl.define(addr, dim2*dim1*dim0);
    }
  void define(void *const addr,
              const DimT &dim1, const DimT &dim0)
    {
      Indexer::ixDefine(0, dim1, dim0);
      m_arrayImpl.define(addr, dim1*dim0);
    }
  void define(void *const addr,
              const DimT &dim0)
    {
      Indexer::ixDefine(0, dim0);
      m_arrayImpl.define(addr, ArSp::sizeOfDim(dim0));
    }
  // With IntVect
  void define(void *const addr, const IntVect &iv,
              const DimT &dimc1, const DimT &dimc0)
    {
      Indexer::ixDefine(0, D_DEFIV, dimc1, dimc0);
      m_arrayImpl.define(addr, iv.product()*dimc1*dimc0);
    }
  void define(void *const addr, const IntVect &iv, const DimT &dimc0)
    {
      Indexer::ixDefine(0, D_DEFIV, dimc0);
      m_arrayImpl.define(addr, iv.product()*dimc0);
    }
  void define(void *const addr, const DimT &dimc1, const DimT &dimc0,
              const IntVect &iv)
    {
      Indexer::ixDefine(0, dimc1, dimc0, D_DEFIV);
      m_arrayImpl.define(addr, dimc1*dimc0*iv.product());
    }
  void define(void *const addr, const DimT &dimc0, const IntVect &iv)
    {
      Indexer::ixDefine(0, dimc0, D_DEFIV);
      m_arrayImpl.define(addr, dimc0*iv.product());
    }
  // With Box
  void define(void *const addr, const Box &box,
              const DimT &dimc1, const DimT &dimc0)
    {
      // Any error with "BoxAllowed" means you are not allowed to use boxes with
      // zero-based lower bounds.  Use ArConf=ArRangeRow or ArConf=ArRangeCol.
      typedef typename ArSp::ArTr<Rank, ArConf>::BoxAllowed Check;
      Indexer::ixDefine(0, D_DEFBOX, dimc1, dimc0);
      m_arrayImpl.define(addr, box.size().product()*dimc1*dimc0);
    }
  void define(void *const addr, const Box &box, const DimT &dimc0)
    {
      // Any error with "BoxAllowed" means you are not allowed to use boxes with
      // zero-based lower bounds.  Use ArConf=ArRangeRow or ArConf=ArRangeCol.
      typedef typename ArSp::ArTr<Rank, ArConf>::BoxAllowed Check;
      Indexer::ixDefine(0, D_DEFBOX, dimc0);
      m_arrayImpl.define(addr, box.size().product()*dimc0);
    }
  void define(void *const addr, const DimT &dimc1, const DimT &dimc0,
              const Box &box)
    {
      // Any error with "BoxAllowed" means you are not allowed to use boxes with
      // zero-based lower bounds.  Use ArConf=ArRangeRow or ArConf=ArRangeCol.
      typedef typename ArSp::ArTr<Rank, ArConf>::BoxAllowed Check;
      Indexer::ixDefine(0, dimc1, dimc0, D_DEFBOX);
      m_arrayImpl.define(addr, dimc1*dimc0*box.size().product());
    }
  void define(void *const addr, const DimT &dimc0, const Box &box)
    {
      // Any error with "BoxAllowed" means you are not allowed to use boxes with
      // zero-based lower bounds.  Use ArConf=ArRangeRow or ArConf=ArRangeCol.
      typedef typename ArSp::ArTr<Rank, ArConf>::BoxAllowed Check;
      Indexer::ixDefine(0, dimc0, D_DEFBOX);
      m_arrayImpl.define(addr, dimc0*box.size().product());
    }

  ///  Deallocate the array
  void undefine()
    {
      CH_assert(isAllocated());
      m_arrayImpl.undefine();
    }

  ///  Assign a constant to the array
  template <typename T2>
  CHArray &operator=(const T2 &val)
    {
      CH_assert(isAllocated());
      T Tval = val;
      const T *const pEnd = end();
      for (T *p = begin(); p != pEnd;) *p++ = Tval;
      return *this;
    }

  ///  Access an element
  T &operator()(const ArSp::IIx_t i6, const ArSp::IIx_t i5,
                const ArSp::IIx_t i4, const ArSp::IIx_t i3,
                const ArSp::IIx_t i2, const ArSp::IIx_t i1,
                const ArSp::IIx_t i0)
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[
        Indexer::ixIndex1D(size(), i6, i5, i4, i3, i2, i1, i0)];
    }
  T &operator()(const ArSp::IIx_t i5, const ArSp::IIx_t i4,
                const ArSp::IIx_t i3, const ArSp::IIx_t i2,
                const ArSp::IIx_t i1, const ArSp::IIx_t i0)
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[
        Indexer::ixIndex1D(size(), i5, i4, i3, i2, i1, i0)];
    }
  T &operator()(const ArSp::IIx_t i4, const ArSp::IIx_t i3,
                const ArSp::IIx_t i2, const ArSp::IIx_t i1,
                const ArSp::IIx_t i0)
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), i4, i3, i2, i1, i0)];
    }
  T &operator()(const ArSp::IIx_t i3, const ArSp::IIx_t i2,
                const ArSp::IIx_t i1, const ArSp::IIx_t i0)
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), i3, i2, i1, i0)];
    }
  T &operator()(const ArSp::IIx_t i2, const ArSp::IIx_t i1,
                const ArSp::IIx_t i0)
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), i2, i1, i0)];
    }
  T &operator()(const ArSp::IIx_t i1, const ArSp::IIx_t i0)
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), i1, i0)];
    }
  T &operator()(const ArSp::IIx_t i0)
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), i0)];
    }
  T &operator()(const IntVect& iv, const ArSp::IIx_t c1, const ArSp::IIx_t c0)
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), D_IXIV, c1, c0)];
    }
  T &operator()(const IntVect& iv, const ArSp::IIx_t c0)
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), D_IXIV, c0)];
    }
  T &operator()(const ArSp::IIx_t c1, const ArSp::IIx_t c0, const IntVect& iv)
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), c1, c0, D_IXIV)];
    }
  T &operator()(const ArSp::IIx_t c0, const IntVect& iv)
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), c0, D_IXIV)];
    }

  ///  Constant access to an element
  const T &operator()(const ArSp::IIx_t i6, const ArSp::IIx_t i5,
                      const ArSp::IIx_t i4, const ArSp::IIx_t i3,
                      const ArSp::IIx_t i2, const ArSp::IIx_t i1,
                      const ArSp::IIx_t i0) const
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[
        Indexer::ixIndex1D(size(), i6, i5, i4, i3, i2, i1, i0)];
    }
  const T &operator()(const ArSp::IIx_t i5, const ArSp::IIx_t i4,
                      const ArSp::IIx_t i3, const ArSp::IIx_t i2,
                      const ArSp::IIx_t i1, const ArSp::IIx_t i0) const
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[
        Indexer::ixIndex1D(size(), i5, i4, i3, i2, i1, i0)];
    }
  const T &operator()(const ArSp::IIx_t i4, const ArSp::IIx_t i3,
                      const ArSp::IIx_t i2, const ArSp::IIx_t i1,
                      const ArSp::IIx_t i0) const
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), i4, i3, i2, i1, i0)];
    }
  const T &operator()(const ArSp::IIx_t i3, const ArSp::IIx_t i2,
                      const ArSp::IIx_t i1, const ArSp::IIx_t i0) const
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), i3, i2, i1, i0)];
    }
  const T &operator()(const ArSp::IIx_t i2, const ArSp::IIx_t i1,
                      const ArSp::IIx_t i0) const
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), i2, i1, i0)];
    }
  const T &operator()(const ArSp::IIx_t i1, const ArSp::IIx_t i0) const
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), i1, i0)];
    }
  const T &operator()(const ArSp::IIx_t i0) const
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), i0)];
    }
  const T &operator()(const IntVect& iv,
                      const ArSp::IIx_t c1, const ArSp::IIx_t c0) const
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), D_IXIV, c1, c0)];
    }
  const T &operator()(const IntVect& iv, const ArSp::IIx_t c0) const
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), D_IXIV, c0)];
    }
  const T &operator()(const ArSp::IIx_t c1, const ArSp::IIx_t c0,
                      const IntVect& iv) const
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), c1, c0, D_IXIV)];
    }
  const T &operator()(const ArSp::IIx_t c0, const IntVect& iv) const
    {
      CH_assert(isUsable());
      return m_arrayImpl.data[Indexer::ixIndex1D(size(), c0, D_IXIV)];
    }

  ///  Get the allocator
  Alloc& getAllocator()
    {
      return *static_cast<Alloc*>(&m_arrayImpl);
    }

  ///  Access extents of the memory allocated for the array
  T* begin()
    {
      return m_arrayImpl.data;
    }
  const T* begin() const
    {
      return m_arrayImpl.data;
    }
  T* end()
    {
      return m_arrayImpl.data + size();
    }
  const T* end() const
    {
      return m_arrayImpl.data + size();
    }

  ///  Overall size of the array
  ArSp::USz_t size() const
    {
      return m_arrayImpl.size;
    }

  ///  Size of a dimension (0 is dimension with contiguous storage)
  ArSp::USz_t size(const unsigned dim) const
    {
      CH_assert(isUsable());
      CH_assert(dim < Rank);
      return Indexer::ixDimSize(size(), dim+1);
    }

  ///  Lower bound of a dimension (0 is dimension with contiguous storage)
  ArSp::IIx_t lowerBound(const unsigned dim) const
    {
      CH_assert(isUsable());
      CH_assert(dim < Rank);
      return Indexer::ixLowerBound(dim+1);
    }

  ///  Upper bound of a dimension (0 is dimension with contiguous storage)
  ArSp::IIx_t upperBound(const unsigned dim) const
    {
      CH_assert(isUsable());
      CH_assert(dim < Rank);
      return Indexer::ixUpperBound(size(), dim+1);
    }

  ///  Memory has been allocated
  bool isAllocated() const
    {
      return (m_arrayImpl.data != 0);
    }

  ///  Memory has been allocated and size is > 0
  bool isUsable() const
    {
      return isAllocated() && size() > 0;
    }

private:

  CHArray(const CHArray&);
  CHArray &operator=(const CHArray&);

  struct Array_impl : public Alloc
  {
    Array_impl()
      : size(0), data(0)
      { }
    Array_impl(const ArSp::USz_t a_size)
      : size(a_size), data(0)
      {
        data = Alloc::allocate(size);
      }
    ~Array_impl()
      {
        if (data != 0) undefine();
      }
    void define(const ArSp::USz_t a_size)
      {
        if (data != 0) undefine();
        size = a_size;
        data = Alloc::allocate(size);
      }
    void define(void *const addr, const ArSp::USz_t a_size)
      {
        if (data != 0) undefine();
        size = a_size;
        data = Alloc::allocate(addr, size);
      }
    void undefine()
      {
        Alloc::deallocate(data, size);
        data = 0;
        size = 0;
      }
    ArSp::USz_t size;                 ///< Overall size of the array
    T *data;                          ///< Data for the array
  };

  Array_impl m_arrayImpl;             ///< Data + allocation/deallocation
};

///  Vector definition
typedef CHArray<Real, 1, ArZeroCol, ArSp::DefaultArrayAlloc<Real> > CHVector;

///  Matrix defined with column-ordered storage
typedef CHArray<Real, 2, ArZeroCol, ArSp::DefaultArrayAlloc<Real> > CHMatrix;

///  Output of an array
template<typename T, unsigned Rank, int ArConf, typename Alloc>
std::ostream &operator<<(std::ostream                           &os,
                         const CHArray<T, Rank, ArConf, Alloc> &A)
{
  CH_assert(A.isAllocated());
  const T *p = A.begin();
  const T *const pEnd = A.end();
  if (p != pEnd) os << *p++;
  while (p != pEnd) os << ' ' << *p++;
  return os;
}

///  Pretty output of a matrix (should be in .cpp)
template<>
inline std::ostream &operator<<(std::ostream &os, const CHMatrix &M)
{
  CH_assert(M.isAllocated());
  const int prec = 2;
  os.setf(std::ios_base::scientific, std::ios_base::floatfield);
  os.precision(prec);
  const ArSp::USz_t iMax = M.size(0);
  const ArSp::USz_t jMax = M.size(1);
  for (ArSp::USz_t i = 0; i != iMax; ++i)
    {
      if (jMax > 0)
        {
          os << std::setw(prec+7) << M(i, 0);
          for (ArSp::USz_t j = 1; j != jMax; ++j)
            {
              os << ' ' << std::setw(prec+7) << M(i, j);
            }
          os << std::endl;
        }
    }
  os.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
  os.precision(6);
  return os;
}

namespace ArSp
{

/*====================================================================*/
/**  \name Special allocators
 *//*=================================================================*/
//@{

/*--------------------------------------------------------------------*/
///  Allocator that only permits allocation by new
/**  Use of CHArray::define(void *const addr, ...) will cause an error
 *//*-----------------------------------------------------------------*/

template <typename T>
class NewArrayAlloc
{
public:
  T* allocate(const USz_t size)
    {
      T *const p = new(std::nothrow) T[size];
      CH_assert(p != 0);
      return p;
    }

  void deallocate(T *p, const USz_t size)
    {
      delete[] p;
    }
};

/*--------------------------------------------------------------------*/
///  Allocator for an array of matrices contiguous in memory
/**  Useful for allocating an array of matrices which can be easily
 *   passed to Fortran, etc.
 *
 *   Example - 1D Array of 2 3x3 matrices                    \verbatim
 *
 *     CHArray<CHMatrix, 1, ArZeroRow, ArSp::ArrayOfMatrixAlloc> AM;
 *     AM.getAllocator().define(3, 3);
 *     AM.define(2);
 *     AM(0) = 1.1;
 *     AM(1) = 1.9;
 *     std::cout << AM << std::endl;
 *                                                        \endverbatim
 *//*-----------------------------------------------------------------*/

class ArrayOfMatrixAlloc
{
  public:
  ArrayOfMatrixAlloc()
    : defined(false)
    { }
  void define(const USz_t a_m, const USz_t a_n)
    {
      m = a_m;
      n = a_n;
      defined = true;
    }
  CHMatrix* allocate(const USz_t size)
    {
      CH_assert(defined);
      CHMatrix *const M = new(std::nothrow) CHMatrix[size];
      CH_assert(M != 0);
      const USz_t stride = m*n;
      matrixData = new(std::nothrow) Real[size*stride];
      CH_assert(matrixData != 0);
      for (USz_t i = 0; i != size; ++i)
        {
          M[i].define(matrixData + i*stride, m, n);
        }
      return M;
    }
  void deallocate(CHMatrix *p, const USz_t size)
    {
      delete[] p;
      delete[] matrixData;
    }
private:
  Real *matrixData;
  USz_t m;
  USz_t n;
  bool defined;
};
//@}
}  // End of namespace ArSp

//--Row ordered speed test

//   StopWatchLinux sw1;
//   CHArray<int, 3, ArZeroRow> AR1;
//   const unsigned nr = 200;
//   AR1.define(nr, nr, nr);
//   AR1 = 1;
//   sw1.start();
//   for (int i = 1; i != nr; ++i)
//     for (int j = 0; j != nr; ++j)
//       for (int k = 0; k != nr; ++k)
//         AR1(i, j, k) += AR1(i-1, j, k);
//   sw1.stop();
//   std::cout << AR1((int)AR1(1, 0, 0), (int)AR1(0, 1, 0), (int)AR1(0, 0, 1))
//             << std::endl;
//   std::cout << "Time: " << sw1.getTime() << std::endl;

//--Column ordered speed test

//   StopWatchLinux sw2;
//   CHArray<int, 3, ArZeroCol> AC1;
//   const unsigned nc = 200;
//   AC1.define(nc, nc, nc);
//   AC1 = 1;
//   sw2.start();
//   for (int i = 1; i != nc; ++i)
//     for (int j = 0; j != nc; ++j)
//       for (int k = 0; k != nc; ++k)
//         AC1(k, j, i) += AC1(k, j, i-1);
//   sw2.stop();
//   std::cout << AC1((int)AC1(1, 0, 0), (int)AC1(0, 1, 0), (int)AC1(0, 0, 1))
//             << std::endl;
//   std::cout << "Time: " << sw2.getTime() << std::endl;

#ifdef CH_SPACEDIM
  #include "NamespaceFooter.H"
#endif
#endif
