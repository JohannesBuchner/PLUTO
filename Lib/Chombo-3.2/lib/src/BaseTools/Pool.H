#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _POOL_H_
#define _POOL_H_

#include <cstdlib>
#include <list>

#include "Vector.H"
#include "BaseNamespaceHeader.H"

class Pool;
typedef std::list<Pool*> PoolList;

///  Pool is a class to optimize memory allocation.
/**
   Pool is a class to optimize memory allocation.  It is specialized to
   allocate fixed size chunks of memory specified by ptrSize in the
   constructor.  Its operation is analogous to malloc, not new.  It
   does not initialize the memory in any way.  The constructor can
   optionally specify an initial pool size, and memory alignment.  The
   pool size will grow as needed by calling ::new.  The pool size never
   shrinks.  Memory will be reclaimed at ~Pool().  The units of poolSize
   are number-of-ptrSize-chunks.  The units of alignment are bytes.

   Pool can only be used for objects of fixed size.  Typically used in
   situations where a class or struct is being constantly constructed and
   deleted.  Objects returned by value from functions, elements in a database,
   etc.  Pool's tend to be allocated statically.

   (from Copier.cpp)
    \code

    // static data member s_motionIemPool getting constructed
    Pool Copier::s_motionItemPool(sizeof(MotionItem), "Copier::MotionItem");

    // note the use of a 'placement new'
    MotionItem* item = new (s_motionItemPool.getPtr()) MotionItem(fromdi, todi, box);

    // if your object does requires it's destructor to be called.
    item->~MotionItem();

    // then return the memory chunk to your pool
    s_motionItemPool.returnPtr(item)

    \endcode

Technical note
   In the event of multi-threading Chombo, we will have to make
   pool access serialized (locks, single thread access, etc) or implement a
   fast lock-free version, or go back to new/delete and let the OS be clever
   again about memory management.

*/
class Pool
{
public:
  ///
  /**
     @param a_ptrSize Size of fixed memory allocations needed

     @param a_poolSize Size of allocations (a_ptrSize*a_poolSize) made
     of the operating system.  This controls the granularity of the
     Pool allocations.  Smaller values result in less wasted memory,
     at the cost of more calls to the operating system for heap memory.

     @param a_name optional name for this Pool.  Used in reporting memory
     by the memory tracking system in Chombo.
  */
  Pool(int         a_ptrSize,
       const char* a_name = "unnamed",
       int         a_poolSize = 100,
       int         a_alignment = sizeof(int),
       bool        a_allowUnalignedAlloc = false);

  ///
  ~Pool();

  /// request a section of memory of  ptrSize_  contiguous bytes.
  void* getPtr();

  /// return memory previous acquired with the getPtr() function.
  void returnPtr(void* a_ptr);

  /// report how much memory this Pool is currently using.
  /**
     memUsage for a Pool only grows until Pool destruction. The Pool
     object has no knowledge of what pieces of memory it has parcelled out
     to a user, so it keeps it all available. The user is responsible for
     not leaking Pool memory.
  */
  long memUsage() const;

  /**
     undocumented function, call this at own risk.  You must
     be absolutely positive that you have returned all the
     ptr's you've asked for, or you can have major seg faults.
  */
  void clear();

  // not for public consumption. used in memory tracking code.
  static PoolList* m_poolList_;
  char             m_name_[64];


protected:


private:
  ///
  Vector<char*> m_pool_;

  ///
  int m_ptrSize_;

  ///
  int m_poolSize_;

  ///
  int m_alignment_;

  bool m_allowUnalignedAlloc;         ///< T - Allows unaligned allocations
                                      ///<     (> sizof(double) if it is not
                                      ///<     known how to allocate aligned
                                      ///<     memory
                                      ///< F - Otherwise abort

  ///
  void* m_next_;

  ///
  void* getMoreMemory();

  void* getMore();

  ///
  /// Not implemented.  Compiler will not let you copy a Pool
  Pool(const Pool& a_rhs);

  /// Not implemented. Compiler will not let you copy a Pool
  const Pool& operator = (const Pool& a_rhs);

  static void clearAllPools();

  friend void dumpmemoryatexit();
};

// Stroustrup (S. 11.5.1) says a declaration in the global scope is needed
void dumpmemoryatexit();

#include "BaseNamespaceFooter.H"
#endif
