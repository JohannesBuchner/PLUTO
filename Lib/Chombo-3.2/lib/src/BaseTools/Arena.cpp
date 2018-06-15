#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iomanip>
#include <utility>
#include "CH_assert.H"
#include <cstdlib>
#include <cstring>
using std::pair;
#include "parstream.H"
#include "memusage.H"
//#include "memtrack.H"
#include "Arena.H"
#include "MayDay.H"
#include "BaseNamespaceHeader.H"

// DON'T include memtrack.H here, we track Arena allocation
// on the usage side, so that we don't have to keep track
// of memory size and pointer lists.  memtrack is too much
// overkill for the Arena, as it can be used cleverly and well
// without any real penalty.

#ifdef CH_USE_MEMORY_TRACKING
ArenaList* Arena::arenaList_ = NULL;
int Arena::NSIZE=120;
#endif

//
// Only really use the coalescing FAB arena if CH_COALESCE_FABS.
//
//#ifdef CH_COALESCE_FABS
//static CArena The_Static_FAB_CArena;
//Arena* The_FAB_Arena = &The_Static_FAB_CArena;
//#else
//static BArena The_Static_FAB_BArena;
//Arena* The_FAB_Arena = &The_Static_FAB_BArena;
//#endif

Arena::Arena()
{
#ifdef CH_USE_MEMORY_TRACKING
  if (arenaList_ == NULL)
  {
    arenaList_ = new ArenaList;
  }

  arenaList_->insert(arenaList_->begin(), this);

  bytes = 0;
  peak = 0;
#endif
}

Arena::~Arena()
{
#ifdef CH_USE_MEMORY_TRACKING
  arenaList_->remove(this);
#endif
}

BArena::BArena(const char* a_name)
{
#ifdef CH_USE_MEMORY_TRACKING
  strncpy(name_, a_name, NSIZE);
  name_[NSIZE-1]=0;
#endif
}

BArena::BArena(const std::string& a_name)
{
#ifdef CH_USE_MEMORY_TRACKING
  strncpy(name_, a_name.c_str(), NSIZE);
  name_[NSIZE-1]=0;
#endif
}

void* BArena::alloc(size_t a_sz)
{
  // We now assume that all allocated memory is initialized.  Our codes
  // sometimes allocates memory that is never used in computations, e.g.,
  // corner ghost cells, but is still copied and/or written to files.  Without
  // some form of initialization, these operations can fail on some machines.
  // Thus, switching back to "malloc" may cause intermittent failure on some
  // machines!
  //void* ret =  calloc(1,a_sz);
  void* ret =  malloc(a_sz);

  if (ret == NULL)
  {
    print_memory_line("Out of memory");
    //#ifdef CH_USE_MEMORY_TRACKING
    //ReportAllocatedMemory(pout());
    //#endif
    pout() << " Trying to calloc(1," << a_sz << ") in BArena::alloc()" << std::endl;
    MayDay::Error("Out of memory in BArena::alloc (BaseFab) ");
  }

  return ret;
}

void BArena::free(void* a_pt)
{
  ::free(a_pt);
}

CArena::CArena(size_t a_hunk_size)
{
  //
  // Force alignment of hunksize.
  //
  m_hunk = Arena::align(a_hunk_size == 0 ? DefaultHunkSize : a_hunk_size);

  CH_assert(m_hunk >= a_hunk_size);
  CH_assert(m_hunk%sizeof(Arena::Word) == 0);
}

CArena::~CArena()
{
  for (unsigned int i = 0; i < m_alloc.size(); i++)
  {
    ::operator delete(m_alloc[i]);
  }
}

void* CArena::alloc(size_t a_nbytes)
{
  a_nbytes = Arena::align(a_nbytes == 0 ? 1 : a_nbytes);
  //
  // Find node in freelist at lowest memory address that'll satisfy request.
  //
  NL::iterator free_it = m_freelist.begin();

  for ( ; free_it != m_freelist.end(); ++free_it)
  {
    if ((*free_it).size() >= a_nbytes)
    {
      break;
    }
  }

  void* vp = 0;

  if (free_it == m_freelist.end())
  {
    vp = ::operator new(a_nbytes < m_hunk ? m_hunk : a_nbytes);

    m_alloc.push_back(vp);

    if (a_nbytes < m_hunk)
    {
      //
      // Add leftover chunk to free list.
      //
      // Insert with a hint -- should be largest block in the set.
      //
      void* block = static_cast<char*>(vp) + a_nbytes;

      m_freelist.insert(m_freelist.end(), Node(block, m_hunk-a_nbytes));
    }
  }
  else
  {
    CH_assert((*free_it).size() >= a_nbytes);

    CH_assert(m_busylist.find(*free_it) == m_busylist.end());

    vp = (*free_it).block();

    if ((*free_it).size() > a_nbytes)
    {
      //
      // Insert remainder of free block back into freelist.
      //
      // Insert with a hint -- right after the current block being split.
      //
      Node freeblock = *free_it;

      freeblock.size(freeblock.size() - a_nbytes);

      freeblock.block(static_cast<char*>(vp) + a_nbytes);

      m_freelist.insert(free_it, freeblock);
    }

    m_freelist.erase(free_it);
  }

  m_busylist.insert(Node(vp, a_nbytes));

  CH_assert(!(vp == 0));

  return vp;
}

void CArena::free(void* a_vp)
{
  if (a_vp == 0)
  {
      //
      // Allow calls with NULL as allowed by C++ delete.
      //
      return;
  }

  //
  // `a_vp' had better be in the busy list.
  //
  NL::iterator busy_it = m_busylist.find(Node(a_vp,0));

  CH_assert(!(busy_it == m_busylist.end()));

  CH_assert(m_freelist.find(*busy_it) == m_freelist.end());

#ifndef NDEBUG
  //[NOTE: this is only used in  CH_assert() below.]
  void* freeblock = static_cast<char*>((*busy_it).block());
#endif
  //
  // Put free'd block on free list and save iterator to insert()ed position.
  //
  pair<NL::iterator,bool> pair_it = m_freelist.insert(*busy_it);

  CH_assert(pair_it.second == true);

  NL::iterator free_it = pair_it.first;

  CH_assert(free_it != m_freelist.end() && (*free_it).block() == freeblock);
  //
  // And remove from busy list.
  //
  m_busylist.erase(busy_it);
  //
  // Coalesce freeblock(s) on lo and hi side of this block.
  //
  if (!(free_it == m_freelist.begin()))
  {
    NL::iterator lo_it = free_it;
    --lo_it;
    void* addr = static_cast<char*>((*lo_it).block()) + (*lo_it).size();

    if (addr == (*free_it).block())
    {
      //
      // This cast is needed as iterators to set return const values;
      // i.e. we can't legally change an element of a set.
      // In this case I want to change the size() of a block
      // in the freelist.  Since size() is not used in the ordering
      // relations in the set, this won't effect the order;
      // i.e. it won't muck up the ordering of elements in the set.
      // I don't want to have to remove the element from the set and
      // then reinsert it with a different size() as it'll just go
      // back into the same place in the set.
      //
      Node* node = const_cast<Node*>(&(*lo_it));
      CH_assert(!(node == 0));
      node->size((*lo_it).size() + (*free_it).size());
      m_freelist.erase(free_it);
      free_it = lo_it;
    }
  }

  NL::iterator hi_it = free_it;

  void* addr = static_cast<char*>((*free_it).block()) + (*free_it).size();

  if (++hi_it != m_freelist.end() && addr == (*hi_it).block())
  {
    //
    // Ditto the above comment.
    //
    Node* node = const_cast<Node*>(&(*free_it));
    CH_assert(!(node == 0));
    node->size((*free_it).size() + (*hi_it).size());
    m_freelist.erase(hi_it);
  }
}

#if 0
void* CArena::calloc(size_t a_nmemb,
                     size_t a_size)
{
  CH_assert(!(a_size == 0));
  CH_assert(!(a_nmemb == 0));

  void* vp = CArena::alloc(a_nmemb*a_size);

  memset(vp, 0, a_nmemb*a_size);

  return vp;
}

void* CArena::realloc(void*  a_ptr,
                      size_t a_size)
{
  if (a_ptr == 0)
  {
    CH_assert(!(a_size == 0));

    return CArena::alloc(a_size);
  }
  else
  {
    if (a_size == 0)
    {
      CArena::free(a_ptr);
    }
    else
    {
      //
      // It had better be in the busy list.
      //
      NL::iterator busy_it = m_busylist.find(Node(a_ptr,0));

      CH_assert(!(busy_it == m_busylist.end()));

      CH_assert(m_freelist.find(*busy_it) == m_freelist.end());

      if (a_size > (*busy_it).size())
      {
        //
        // If the new size is larger than old allocate new chunk.
        //
        void* vp = CArena::alloc(a_size);

        memcpy(vp, a_ptr, (*busy_it).size());

        CArena::free(a_ptr);

        a_ptr = vp;
      }
    }

    return a_ptr;
  }
}
#endif

#include "BaseNamespaceFooter.H"
