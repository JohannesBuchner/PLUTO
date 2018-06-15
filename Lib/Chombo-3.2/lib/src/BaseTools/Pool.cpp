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

#include "CH_config.H"
#include "CH_System.H"
#include "Pool.H"
#include "MayDay.H"
#include "BaseNamespaceHeader.H"

PoolList* Pool::m_poolList_ = NULL;

bool PoolListInit = false;

#ifdef CH_USE_MEMORY_TRACKING
// Added so we can choose not to delete them in the clearAllPools() func
extern Pool chunkPool;
extern Pool infoPool;
#endif

// Constructor sets Pool parameters, but doesn't actually allocate any memory.
Pool::Pool(int         a_ptrSize,
           const char* a_name,
           int         a_poolSize,
           int         a_alignment,
           bool        a_allowUnalignedAlloc)
  :
  m_pool_(0),                // An array of pointers to pools
  m_ptrSize_(a_ptrSize),     // Chunk size (bytes) in each pool
  m_poolSize_(a_poolSize),   // Number of chunks in each pool
  m_alignment_(a_alignment), // allocated memory to be aligned on this boundary.
  m_allowUnalignedAlloc(a_allowUnalignedAlloc),
                             // Whether unaligned allocations > sizeof(double)
                             // are allowed
  m_next_(0)                 // Points to next available chunk or is null.
{

  if (m_poolSize_ < 1)
  {
    MayDay::Error("Pool::Pool non positive poolSize.");
  }

  // Minimum chunk size is size of pointer.
  if ((unsigned int)m_ptrSize_ < sizeof(void*))
  {
    m_ptrSize_ = sizeof(void*);
  }

  // m_alignment_ checked out only for 4 and 8.
  if ((unsigned int)m_alignment_ < sizeof(int))
  {
    m_alignment_ = sizeof(int);
  }

  // chunk size must be multiple of m_alignment_
  while ((m_ptrSize_ % m_alignment_) != 0)
  {
    m_ptrSize_++;
  }

  if (!PoolListInit)
  {
    m_poolList_ = new PoolList;
    PoolListInit = true;
  }

  m_poolList_->insert(m_poolList_->begin(), this);
  strcpy(m_name_, a_name);
}

Pool::~Pool()
{
  for (int i = 0; i < m_pool_.size(); i++)
  {
    delete [] m_pool_[i];
  }

  m_pool_.resize(0);
  m_poolList_->remove(this);
  //[NOTE: so LeakTracer wont complain. <dbs>]
  if (m_poolList_->empty())
  {
    delete m_poolList_ ;
    PoolListInit = false ;
  }
}

/// getPtr and returnPtr must be as fast as possible!
void* Pool::getPtr()
{
  if (m_next_ == 0)
  {
    m_next_ = getMoreMemory();
  }

  void* result = m_next_;      // result points to first free chunk in list
  m_next_ = *(void**)m_next_;  // point m_next_ to next free chunk in list.

  return result;
}

void Pool::returnPtr(void* ptr)
{
  *(void**)ptr = m_next_;  // Stamp *ptr with the first free chunk.
  m_next_ = ptr;           // Make ptr the first free chunk.
}

//  notes on performing getPtr and returnPtr in a lock-free synchronized manner.
//  currently this design will suffer from the 'ABA' problem of CAS1-based
//  synchronization.  I need to use a CAS2 or 'double-compare-single-swap'

// void* Pool::getPtrAtomic()
// {
//   void* newNext = 0;
//   if (m_next_ == 0)
//     {
//       newNext = getMore();
//       if ( CAS(&m_next_, 0, newNext))
//         m_pool_.push_back(po); // technically we should use a stack with atomic push
//       else
//         delete newNext; //some other thread beat me there, throw this one away.
//     }

//   void* result;
//   void* next;
//   void* check;
//   do{
//     result = m_next_;
//     next = *(void**)m_next_;
//   }while (CAS(&m_next_, result, next));

//   return result;
// }

// void Pool::returnPtrAtomic(void* ptr)
// {
//   void* next;
//   do{
//     next = m_next_;
//     *(void**)ptr = next;  // Stamp *ptr with the first free chunk.
//   }while (CAS( &m_next_, next,  ptr)); // Make ptr the first free chunk.
// }

long Pool::memUsage() const
{
  return m_ptrSize_ * m_poolSize_ * m_pool_.size();
}

void Pool::clear()
{
  for (int i = 0; i < m_pool_.size(); i++)
  {
    delete [] m_pool_[i];
  }
  m_next_ = 0;
  m_pool_.resize(0);
}

void Pool::clearAllPools()
{
  PoolList::iterator   p;

  for (p = Pool::m_poolList_->begin(); p != Pool::m_poolList_->end(); p++)
    {
#ifdef CH_USE_MEMORY_TRACKING
      // Don't clear the special memtrack Pools here.
      // This fixed some errors I was seeing in dumpmemoryatexit(). Added 3/09 (ndk)
      if (*p != &chunkPool && *p != &infoPool)
#endif
        {
          (*p)->clear();
        }
    }
}

void* Pool::getMoreMemory()
{
  void* po = getMore();
  m_pool_.push_back((char*)po);
  return po;
}

void* Pool::getMore()
{
  // Create another pool
  char* po;

  if (m_alignment_ == sizeof(int))
  {
    po = (char*)new int [m_ptrSize_ * m_poolSize_ / sizeof(int)];
  }
  else if (m_alignment_ == sizeof(double))
  {
    po = (char*)new double [m_ptrSize_ * m_poolSize_ / sizeof(double)];
  }
  else
  {
#ifndef CHDEF_SYSTEM_HAVE_ALIGNEDALLOC
    if (!m_allowUnalignedAlloc)
      {
        MayDay::Error("Pool::Pool : alignments greater than int or double are "
                      "not supported");
      }
#endif
    void* addr;
    int ier = CHSystem::memalign(&addr, m_alignment_, m_ptrSize_*m_poolSize_);
    po = static_cast<char*>(addr);
    if (ier != 0)
      {
        MayDay::Error("Pool::getMore() : System::memalign failed");
      }
  }

  // Thread pool:
  // The first few bytes in each chunk is interpreted as a pointer
  // which points to the next free chunk.  With multiple pools, and random
  // news and deletes, the pools will become threaded together.  This is ok.
  char* p = po;

  for (int i = 0; i < m_poolSize_ - 1; i++)
  {
    *(char**)p = p + m_ptrSize_;  // Chunk at i points to chunk at i+1.
    p += m_ptrSize_;
  }

  *(char**)p = 0;  // Chunk at end of list points to null
  return po;       // return first chunk.
}
#include "BaseNamespaceFooter.H"
