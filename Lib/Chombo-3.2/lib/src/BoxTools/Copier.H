#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COPIER_H_
#define _COPIER_H_

#include "DisjointBoxLayout.H"
#include "Pool.H"
#include "Vector.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

class CopyIterator;

//  These classes are public because I can't find a nice
//  way to make a class a friend of all the instantiations
//  of a template class.  These classes are not part of
//  the public interface for the Array API.
//
//  Later, if MotionItem shows up in the profiler, we
//  can start using a pool allocation scheme and placement new

class MotionItem
{
public:
  DataIndex fromIndex, toIndex;
  Box fromRegion;
  Box toRegion;
  int procID;
  // this constructor will probably eventually go away
  MotionItem(const DataIndex& a_from,
             const DataIndex& a_to,
             const Box&       a_region);
  MotionItem(const DataIndex& a_from,
             const DataIndex& a_to,
             const Box&       a_fromRegion,
             const Box&       a_toRegion);
  void reverse();

};

class CopierBuffer
{
public:

  ///null constructor, copy constructor and operator= can be compiler defined.
  CopierBuffer():m_ncomps(0), m_sendbuffer(NULL), m_sendcapacity(0),
                 m_recbuffer(NULL), m_reccapacity(0)
  {}

  ///
  virtual ~CopierBuffer();

  void clear();

  bool isDefined(int ncomps) const
  { return ncomps == m_ncomps;}

  mutable int m_ncomps;

  mutable void*  m_sendbuffer; // pointer member OK here,
                               // since LevelData<T> has no copy
  mutable size_t m_sendcapacity;
  mutable void*  m_recbuffer;  // pointer member OK here,
                               // since LevelData<T> has no copy
  mutable size_t m_reccapacity;

#ifndef DOXYGEN

  struct bufEntry
  {
    void* bufPtr; // pointer into contiguous m_buffer
    size_t size;
    const MotionItem* item;
    unsigned int procID;
    bool operator < (const bufEntry& rhs) const
      {
        if (procID == rhs.procID)
          {
            const Box& left = item->toRegion;
            const Box& right= rhs.item->toRegion;
            if (left.smallEnd() == right.smallEnd())
              {
                return left.bigEnd().lexLT(right.bigEnd());
              }
            else
              {
                return item->toRegion < rhs.item->toRegion;
              }
          }
        //else
        return procID < rhs.procID;
      }
  };

#endif
  mutable std::vector<bufEntry> m_fromMe;
  mutable std::vector<bufEntry> m_toMe;


protected:

};

/// A strange but true thing to make copying from one boxlayoutdata to another fast
/**
  class to handle the organization of data from a LevelData to a BoxLayoutData
  defined over the same index space.  The idea behind this object
  is twofold:
  A)  there is a significant amount of information that can be
  computed once and reused from one copy to the next for
  a LevelData and BoxLayoutData that have unchanging BoxLayouts
  underneath them.  In these cases, LevelData::copyTo methods
  has an optional argument.  That way, in cases where the operation
  is only performed a single time, we do not waste our time building
  optimizing data, etc.
  B)  when we interface with KeLP, this class maps quite well
  with the existing KeLP API (class MotionPlan).

  pains have been taken to not have this class be templated.  Thus
  it does not actually perform the copying itself (which would require
  knowledge of the template class type).  It provides an iterator
  interface to MotionItems that specify the required work.  LevelData
  provides the copyTo methods.
*/

class Copier
{
public:

  ///null constructor, copy constructor and operator= can be compiler defined.
  Copier():m_isDefined(false)
  {}

  Copier(const Copier& a_rhs);

  ///
  Copier(const DisjointBoxLayout& a_level,
         const BoxLayout& a_dest,
         bool a_exchange = false,
         IntVect a_shift = IntVect::Zero);

  /// this constructor contains support for periodic BC's
  Copier(const DisjointBoxLayout& a_level,
         const BoxLayout& a_dest,
         const ProblemDomain& a_domain,
         bool a_exchange  = false,
         IntVect a_shift = IntVect::Zero);

  /// Copier to copy data into the valid and invalid regions of a_dest from a_level
  Copier(const DisjointBoxLayout& a_level,
         const BoxLayout& a_dest,
         const IntVect& a_destGhost,
         bool a_exchange = false,
         IntVect a_shift = IntVect::Zero);

  /// Copier to copy data into the valid and invalid regions of a_dest from a_level with support for periodic BCs
  Copier(const DisjointBoxLayout& a_level,
         const BoxLayout& a_dest,
         const ProblemDomain& a_domain,
         const IntVect& a_destGhost,
         bool  a_exchange = false,
         IntVect a_shift = IntVect::Zero);

  ///
  virtual ~Copier();

  ///
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      bool  a_exchange = false,
                      IntVect a_shift = IntVect::Zero);

  /// contains support for periodic BCs
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      bool  a_exchange = false,
                      IntVect a_shift = IntVect::Zero);

  ///
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const IntVect& a_destGhost,
                      bool  a_exchange = false,
                      IntVect a_shift = IntVect::Zero);

  /// contains support for periodic BCs
  virtual void define(const BoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      const IntVect& a_destGhost,
                      bool  a_exchange = false,
                      IntVect a_shift = IntVect::Zero);

  ///  alternative version of define for exchange copiers that uses new optimized neighborIterator
  virtual void exchangeDefine(const DisjointBoxLayout& a_grids,
                              const IntVect& a_ghost);

  /// performs the computation necessary for moving from ghost+valid cells to all valid cells in a_dest
  /** instead of copying from the valid of a_src to the valid+ghost of a_dest,
      this routine defines a Copier that moves data from the ghost+valid region
      of a_src to  a_dest.  This is used by LevelFluxRegister.
  */
  void ghostDefine(const DisjointBoxLayout& a_src,
                   const DisjointBoxLayout& a_dest,
                   const ProblemDomain& a_domain,
                   const IntVect& a_srcGhost);

  /// performs the computation necessary for moving from ghost+valid cells to all valid+ghost cells in a_dest
  /** instead of copying from the valid of a_src to the valid+ghost of a_dest,
      this routine defines a Copier that moves data from the ghost+valid region
      of a_src to valid+ghost a_dest.  This is used by BlockBaseRegister.  Note that this does not handle 
      periodic matchings here.  That is handled by different code in the MappedMultiBlock libraries.
      Sending in a null-constructed ProblemDomain just results in a computation that ignores domain boundaries.
  */
  void ghostDefine(const DisjointBoxLayout& a_src,
                   const DisjointBoxLayout& a_dest,
                   const ProblemDomain& a_domain,
                   const IntVect& a_srcGhost,
                   const IntVect& a_destGhost);

  ///
  virtual void clear();

  ///
  /**
     The assignment operator.

  */
  Copier& operator= (const Copier& b);

  ///
  /**
     reverse the communication of this Copier.
  */
  void reverse();

  ///
  /**
    Under the right algorithmic circumstances, a user might know that they no not
    need corner cells filled in exchanges, what they really need are just sufficient
    ghost cells in each CH_SPACEDIM direction for each valid cell.  This operation
    eliminates copy operations that do not contribute to such cells.  Under reasonable
    3D configurations this can eliminate more than half the copy operations.
  */
  void trimEdges(const DisjointBoxLayout& a_exchangedLayout, const IntVect& ghost);

  ///
  /*
    OK, if you *really* know what you are doing (!)  you can coarsen a Copier in place of
    building a new Copier on a coarsened DisjointBoxLayout.  The trick is when you do have
    ghost cells and working out on paper for yourself that the coarsening is valid. One
    common configuration is an exchange Copier with *one* ghost cell.  A particular
    case where coarsen will do the *wrong* thing is an exchange Copier with more than one
    ghost cell.
  **/
  void coarsen(int a_refRatio);

  int print() const;

  bool bufferAllocated() const;
  void setBufferAllocated(bool arg) const;

  int numLocalCellsToCopy() const;
  int numFromCellsToCopy() const;
  int numToCellsToCopy() const;

  bool isDefined() const
  { return m_isDefined;}

  CopierBuffer  m_buffers;

protected:

  friend class CopyIterator;

  Vector<MotionItem*> m_localMotionPlan;
  Vector<MotionItem*> m_fromMotionPlan;
  Vector<MotionItem*> m_toMotionPlan;

  friend void dumpmemoryatexit();
  static Pool s_motionItemPool;
  mutable bool buffersAllocated;

  // keep a refcounted reference around for debugging purposes, we can
  // decide afterwards if we want to eliminate it.
  DisjointBoxLayout m_originPlan;
  BoxLayout  m_dest;

  bool m_isDefined;

  void trimMotion(const DisjointBoxLayout& a_exchangedLayout, const IntVect& a_ghost,
                  const Vector<MotionItem*>& a_oldItems, Vector<MotionItem*>& a_newItems);

  void sort();

  // sneaky end-around to problem of getting physDomains in derived classes
  const ProblemDomain& getPhysDomain(const DisjointBoxLayout& a_level) const;
};

std::ostream& operator<<(std::ostream& os, const Copier& copier);

//===========================================================================
// end of public interface for Copier.
//===========================================================================

inline MotionItem::MotionItem(const DataIndex& a_from,
                              const DataIndex& a_to,
                              const Box&       a_region)
  :fromIndex(a_from), toIndex(a_to), fromRegion(a_region),
  toRegion(a_region), procID(-1)
{
}

inline MotionItem::MotionItem(const DataIndex& a_from,
                              const DataIndex& a_to,
                              const Box&       a_fromRegion,
                              const Box&       a_toRegion)
  :fromIndex(a_from), toIndex(a_to), fromRegion(a_fromRegion),
  toRegion(a_toRegion), procID(-1)
{
}

inline
void  MotionItem::reverse()
{
  Box tmp(fromRegion);
  fromRegion=toRegion;
  toRegion=tmp;
  DataIndex tmpIndex(fromIndex);
  fromIndex = toIndex;
  toIndex = tmpIndex;
}

class CopyIterator
{
public:
  enum local_from_to
  {
    LOCAL,
    FROM,
    TO
  };

  inline CopyIterator(const Copier& a_copier, local_from_to);

  inline const MotionItem& operator()() const;

  inline const MotionItem& operator[](size_t index) const;

  inline void operator++();

  inline bool ok() const;

  inline void reset();

  inline size_t size(){ return m_motionplanPtr->size();}

private:
  const Vector<MotionItem*>* m_motionplanPtr;
  unsigned int  m_current;
};

//====== inlined functions====================================

inline CopyIterator::CopyIterator(const Copier& a_copier, local_from_to type)
  :m_current(0)
{
  switch(type)
  {
  case LOCAL:
    m_motionplanPtr = &(a_copier.m_localMotionPlan);
    break;
  case FROM:
    m_motionplanPtr = &(a_copier.m_fromMotionPlan);
    break;
  case TO:
    m_motionplanPtr = &(a_copier.m_toMotionPlan);
    break;
  default:
    MayDay::Error("illegal local_from_to option for CopyIterator");
  }
}

inline const MotionItem& CopyIterator::operator()() const
{
  CH_assert(m_current < m_motionplanPtr->size());
  return *(m_motionplanPtr->operator[](m_current));
}

inline const MotionItem& CopyIterator::operator[](size_t a_index) const
{
  CH_assert(a_index < m_motionplanPtr->size());
  return *(m_motionplanPtr->operator[](a_index));
}

inline void CopyIterator::operator++()
{
  ++m_current;
}

inline bool CopyIterator::ok() const
{
  return m_current < m_motionplanPtr->size();
}

inline void CopyIterator::reset()
{
  m_current = 0;
}

#include "NamespaceFooter.H"
#endif
