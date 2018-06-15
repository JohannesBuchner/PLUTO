#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ARENA_H_
#define _ARENA_H_

#include <cstddef>

#ifdef CH_USE_MEMORY_TRACKING
#include <list>
#include <cstddef>
#endif

#include <set>
#include <vector>
#include "BaseNamespaceHeader.H"

#ifdef CH_USE_MEMORY_TRACKING
class Arena;
typedef std::list<Arena*> ArenaList;
#endif

using std::set;

/// A Virtual Base Class for Dynamic Memory Managemen
/**
A Virtual Base Class for Dynamic Memory Management

  This is a virtual base class for objects that manage their own dynamic
  memory allocation.  Since it is a virtual base class, you have to derive
  something from it to use it.
*/
class Arena
{
public:
  /// base class constructor
  /** base class constructor.  registers Arena object with
      memory tracking system
  */
  Arena();

  /// base class destructor.
  /** base class destructor.  unregisters Arena object with
      memory tracking system
  */
  virtual ~Arena();

  /**
     Allocate a dynamic memory arena of size a_sz.
     A pointer to this memory should be returned.
  */
  virtual void* alloc(size_t a_sz) = 0;

  /**
     A pure virtual function for deleting the arena pointed to by a_pt.
  */
  virtual void free(void* a_pt) = 0;

  /**
     Given a minimum required arena size of a_sz bytes, this returns
     the next largest arena size that will hold an integral number
     of objects of the largest of the types void*, long,
     double and function pointer.
  */
  static size_t align(size_t a_sz);

  typedef void (*FP)();

#ifdef CH_USE_MEMORY_TRACKING
  /**data members used by memory tracking system.*/

  /*@{*/
  long int bytes;
  long int peak;
  static ArenaList* arenaList_;
  static int NSIZE;
  char name_[120];
  /*@}*/
#endif


protected:
    //
    // Types used by align().
    //
#ifndef DOXYGEN

    union Word
    {
        void*  p;
        double d;
        long   l;
        FP     f;
    };
#endif
};

//
// Inlines.
//

inline size_t Arena::align (size_t a_s)
{
    size_t x = a_s + sizeof(Word) - 1;
    x -= x % sizeof(Word);
    return x;
}

/// A Concrete Class for Dynamic Memory Management
/**

  This is the simplest dynamic memory management class derived from Arena.

  Makes calls to ::operator new() and ::operator delete().
*/

class BArena: public Arena
{
public:
  ///
  /**
   optional @param a_name used by memory tracker to distinguish
   between different memory Arenas
  */
  BArena(const char* a_name = "unnamed");

  BArena(const std::string& a_name);

  /**: Allocates a dynamic memory arena of size a_sz.  Returns a
       pointer to this memory.
  */
  virtual void* alloc (size_t a_sz);

  /// Deletes the arena pointed to by a_pt.
  virtual void free (void* a_pt);
};

/// A Concrete Class for Dynamic Memory Management
/**

  This is a coalescing memory manager.  It allocates (possibly) large
  chunks of heap space and apportions it out as requested.  It merges
  together neighboring chunks on each free().
*/

class CArena: public Arena
{
public:
    /**: Construct a coalescing memory manager.  `a_hunk_size' is the
               minimum size of hunks of memory to allocate from the heap.
               If a_hunk_size == 0 we use DefaultHunkSize as specified below.
    */
    CArena(size_t a_hunk_size = 0);

    /// The destructor.
    virtual ~CArena();

    /// Allocate some memory.
    virtual void* alloc(size_t a_nbytes);

    /** Free up allocated memory.  Merge neighboring free memory chunks
               into largest possible chunk.
    */
    virtual void free(void* a_vp);

#if 0
    /** Mirror the C calloc() function.  Returns zero'd memory of
               size a_nmemb*a_size.  This is here so that we can implement
               malloc(3) and family.  Users shouldn't use this function.
    */
    void* calloc(size_t a_nmemb,
                 size_t a_size);

    /** Mirror the C realloc() function.  This is here so that we can
               implement malloc(3) and family.  Users shouldn't use this
               function.
    */
    void* realloc (void*  a_ptr,
                   size_t a_size);
#endif

    /// The default memory hunk size to grab from the heap.
    enum
    {
      DefaultHunkSize = 1024*1024
    };


protected:
    //
    // The nodes in our free list and block list.
    //
#ifndef DOXYGEN
    class Node
    {
    public:
      //
      // The default constructor.
      //
      Node()
        :
        m_block(0),
        m_size(0)
      {}
      //
      // Another constructor.
      //
      Node(void* a_block, size_t a_size)
        :
        m_block(a_block),
        m_size(a_size)
      {}
      //
      // The copy constructor.
      //
      Node(const Node& a_rhs)
        :
        m_block(a_rhs.m_block),
        m_size(a_rhs.m_size)
      {}
      //
      // The copy assignment constructor.
      //
      Node& operator = (const Node& a_rhs)
      {
          m_block = a_rhs.m_block;
          m_size  = a_rhs.m_size;
          return *this;
      }
      //
      // The "less-than" operator.
      //
      bool operator < (const Node& a_rhs) const
      {
          return m_block < a_rhs.m_block;
      }
      //
      // The equality operator.
      //
      bool operator == (const Node& a_rhs) const
      {
          return m_block == a_rhs.m_block;
      }
      //
      // The block address.
      //
      void* block() const
      {
        return m_block;
      }
      //
      // Set block address.
      //
      void block (void* a_blk)
      {
        m_block = a_blk;
      }
      //
      // The size of the memory block.
      //
      size_t size() const
      {
        return m_size;
      }
      //
      // Set size.
      //
      void size(size_t a_sz)
      {
        m_size = a_sz;
      }

    private:
        //
        // The block of memory we reference.
        //
        void* m_block;
        //
        // The size of the block we represent.
        //
        size_t m_size;
    };

    //
    // The type of our freelist and blocklist.
    // We use a set sorted from lo to hi memory addresses.
    //
    typedef set < Node > NL;

    //
    // The list of blocks allocated via ::operator new().
    //
    std::vector<void*> m_alloc;

    //
    // The free list of allocated but not currently used blocks.
    // Maintained in lo to hi memory sorted order.
    //
    NL m_freelist;

    //
    // The list of busy blocks.
    // A block is either on the freelist or on the blocklist, but not on both.
    //
    NL m_busylist;

    //
    // The minimal size of hunks to request via ::operator new().
    //
    size_t m_hunk;
#endif


private:
    //
    // Disallowed.
    //
    CArena (const CArena& a_rhs);
    CArena& operator= (const CArena& a_rhs);
};

//
// The Arena used by BaseFab code.
//
extern Arena* The_FAB_Arena;

#include "BaseNamespaceFooter.H"
#endif /*CH_ARENA*/
