#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// do not include memtrack.H here, since we don't want to
// track the allocation tracker itself.
//
//  addendum.  ha ha , my tracker outfoxed me anyways.  since
//  our heavy use allocators have their own tracking mechanisms
//  (BaseFab, Pool, Vector) I end up tracking some of the trackers
//  memory anyways.  I'm aware of it all, just have to be careful
//  about things.  bvs

#ifdef CH_USE_MEMORY_TRACKING

#include <stdio.h>
#include <list>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
using std::list;
using std::cout;
using std::endl;
using std::string;
#ifndef CH_DISABLE_SIGNALS
#include <csignal>
#endif

#include "Arena.H"
#include "parstream.H"
#include "Vector.H"
#include "Pool.H"
#include "BitSet.H"
#include "MayDay.H"
#include "SPMD.H"
#include "memtrack.H"

#include "BaseNamespaceHeader.H"

static const int    BYTES_PER_MEG      = 1024*1024;
static const size_t CHUNK_FILENAME_LEN = 128;

void UnfreedMemory();

// function to register for listening to signal SIGUSR1 to provide a snapshot
// of the current memory layout
void dumpmemorymap(int a_sig)
{
  UnfreedMemory(); // act on signal

#ifndef CH_DISABLE_SIGNALS
  signal(a_sig, dumpmemorymap); // reset self for future requests
#endif
}

void dumpmemoryatexit()
{
  Pool::clearAllPools();
  //Copier::s_motionItemPool.clear();
  //IntVectSet::clearStaticMemory();

  //#ifdef PTHREAD
  //  ThreadTask::taskPool.clear();
  //#endif

  // this code was from my mad debugging of the static order
  // initialization bug in the memory tracking code that was
  // making the optimized code behave differently than the
  // debug code, because static initialization order is not
  // preserved during an optimizaed compilation.  bvs
  //
  // incr::iterator i;
  // lst::iterator  L;
  // for (i = vectorIncr_->begin(); i != vectorIncr_->end(); i++)
  // {
  //   cout << i->first << "  ";
  //   for (L = i->second->begin(); L != i->second->end(); L++)
  //   {
  //     cout << " " << *L;
  //   }
  //   cout<<"\n";
  // }

  UnfreedMemory();

  // delete  Arena::arenaList_;
  // delete  Pool::m_poolList_;
  // delete  vectorList_;
  // delete  vectorIncr_;
  // delete  Arena::s_Arena;
}

void dumpmemoryabort(int sig)
{
  UnfreedMemory();
#ifndef CH_DISABLE_SIGNALS
  signal(sig, SIG_DFL);
  raise(sig);
#endif
}

int registerMemorySignals()
{
#ifndef CH_CYGWIN
#ifndef CH_MPI
  atexit(dumpmemoryatexit);
#endif
#endif

#ifndef CH_DISABLE_SIGNALS
  // signal(SIGABRT, dumpmemoryabort);
  signal(SIGUSR2, dumpmemorymap);
#endif

  return 1;
}

VectorList* vectorList_ = NULL;
incr*       vectorIncr_ = NULL;

typedef struct
{
  void*  address;
  size_t size;
}
chunk;

typedef struct
{
  std::list<chunk*> chunks;
  char              file[CHUNK_FILENAME_LEN];
  int               line;
  long              total;
  long              peak;
  bool              mlloc;
}
ALLOC_INFO;

Pool chunkPool(sizeof(chunk),      "memtrack:chunk");
Pool infoPool (sizeof(ALLOC_INFO), "memtrack:locations");

typedef std::list<ALLOC_INFO*> AllocList;
AllocList *allocList;

unsigned long grandTotal = 0;
unsigned long grandPeak  = 0;

bool trackingOn = true;

static int registerAbortHandler = registerMemorySignals();

void ReportUnfreedMemory(ostream& a_os)
{
  AllocList::iterator  i;
  PoolList::iterator   p;
  ArenaList::iterator  a;
  VectorList::iterator v;

  long int totalSize = 0;
  long int peak = 0;

  if (trackingOn && allocList)
  {
    for (i = allocList->begin(); i != allocList->end(); i++)
    {
      peak += (*i)->peak;

      if ((*i)->total != 0)
      {
        totalSize += (*i)->total;

        // a_os<< "-----------------------------------------------------------\n";
        a_os << (*i)->file << ":" << (*i)->line << ":  " << (*i)->total << "\n";
      }
    }

    if (totalSize != 0)
    {
      a_os << "Total Unfreed from new,malloc, etc: " << totalSize << "  bytes\n";
    }
  }

  if (vectorList_ != NULL)
  {
    // a_os<<"---------------------------------\n";
    // a_os<<"Vector-based memory allocation\n";

    for (v = vectorList_->begin(); v != vectorList_->end(); v++)
    {
      long int* bytes = (v->second.first);

      if (*bytes != 0 && *bytes != 8)
      {
        a_os << "Vector " << v->first << ": "
             << *bytes << " bytes (" << (*bytes)/BYTES_PER_MEG << " Mb)\n";
        totalSize += *bytes;
      }

      peak += *(v->second.second);

      // probably wondering about my silly "8" check.  workaround for
      // now since I can't clear chunkPool and infoPool unless this
      // is an atexit call, since they need to stay fresh in case this
      // is a transient run call to memusage.  will think up a better
      // way for this at another time.  Can't kill myself for 8 bytes
      //   bvs.
    }
  }

  // a_os<<"---------------------------------\n";
  // a_os<<"Pool-based memory allocation\n";

  for (p = Pool::m_poolList_->begin(); p != Pool::m_poolList_->end(); p++)
  {
    peak += (*p)->memUsage();

    if (*p == &chunkPool || *p == &infoPool)
    {
      // a_os << (*p)->m_name_ << ": "
      //      << (*p)->memUsage() << " bytes ("
      //      << (*p)->memUsage()/BYTES_PER_MEG << " Mb) (not added to total)\n";
    }
    else if ((*p)->memUsage() != 0)
    {
      a_os << "Pool "<<(*p)->m_name_ << ": "
           << (*p)->memUsage() << " bytes ("
           << (*p)->memUsage()/BYTES_PER_MEG << " Mb)\n";
      totalSize += (*p)->memUsage();
    }
  }

  if (Arena::arenaList_ != NULL)
  {
    // a_os<<"---------------------------------\n";
    // a_os<<"Arena-based memory allocation\n";

    for (a = Arena::arenaList_->begin(); a != Arena::arenaList_->end(); a++)
    {
      peak += (*a)->peak;

      if ((*a)->bytes != 0)
      {
        a_os << "Arena " << (*a)->name_ << ": "
             << (*a)->bytes << " bytes ("
             << (*a)->bytes/BYTES_PER_MEG << " Mb)\n";
        totalSize += (*a)->bytes;
      }
    }
  }

  peak += BitSet::peak;

#if 0
  // (DFM -- 10/22/08) -- for now, comment out memtracking of
  // IntVectSet because sizeof(IntVectSet) is dimensionally dependent,
  // so we need to figure out how to keep track of it in a
  //dimensionally-independent way
  peak += IntVectSet::peakcount*sizeof(IntVectSet);
#endif

  if (BitSet::bytes != 0)
  {
    // a_os<<"---------------------------------\n";
    a_os << " BitSet allocation "
         << BitSet::bytes << " bytes ("
         << BitSet::bytes/BYTES_PER_MEG << " Mb)\n";
    totalSize += BitSet::bytes;
  }

#if 0
  // (DFM -- 10/22/08) -- for now, comment out memtracking of
  // IntVectSet because sizeof(IntVectSet) is dimensionally dependent,
  // so we need to figure out how to keep track of it in a
  //dimensionally-independent way
  if (IntVectSet::count != 0)
  {
    a_os << " IntVectSet allocation "
         << IntVectSet::count * sizeof(IntVectSet) << " bytes ("
         << IntVectSet::count * sizeof(IntVectSet)/BYTES_PER_MEG << " Mb)\n";
    totalSize += IntVectSet::count*sizeof(IntVectSet);
  }
#endif

  if (totalSize != 0)
  {
    a_os << "---------------------------------\n";
    a_os << "Total Unfreed : "
         << totalSize << " bytes ("
         << totalSize/BYTES_PER_MEG << " Mb)\n";
  }

  a_os << "peak memory usage: "
       << peak << " bytes ("
       << peak/BYTES_PER_MEG << " Mb)\n";
}

void ReportAllocatedMemory(ostream& a_os)
{
  AllocList::iterator  i;
  PoolList::iterator   p;
  ArenaList::iterator  a;
  VectorList::iterator v;

  long int totalSize = 0;
  long int peak = 0;
  char temp[1024];

  if (trackingOn && allocList)
    {
      for (i = allocList->begin(); i != allocList->end(); i++)
        {
          peak += (*i)->peak;

          if ((*i)->total != 0)
            {
              totalSize += (*i)->total;
              a_os << (*i)->file << ":" << (*i)->line << ":  " << (*i)->total << "\n";
            }
        }

      if (totalSize != 0)
        {
          a_os << "Total Unfreed from callocMT,mallocMT,reallocMT: " << totalSize << "  bytes\n";
        }
    }

  if (vectorList_ != NULL)
    {
      for (v = vectorList_->begin(); v != vectorList_->end(); v++)
        {
          long int* bytes = (v->second.first);

          if (*bytes != 0 && *bytes != 8)
            {
              string entry = "Vector";
              sprintf(temp, "%11s %-40s %12ld b  %11.4f Mb  peak=%11.4f\n", entry.data(), v->first.data(),
                      *bytes, (Real)(*bytes)/(Real)BYTES_PER_MEG,
                      (Real)(*(v->second.second))/(Real)BYTES_PER_MEG);
              a_os << temp;
              totalSize += *bytes;
            }

          peak += *(v->second.second);

        }
    }

  for (p = Pool::m_poolList_->begin(); p != Pool::m_poolList_->end(); p++)
    {
      peak += (*p)->memUsage();

      if (*p == &chunkPool || *p == &infoPool)
        {
          // what is this???
          // a_os << (*p)->m_name_ << ": "
          //      << (*p)->memUsage() << " bytes ("
          //      << (*p)->memUsage()/BYTES_PER_MEG << " Mb) (not added to total)\n";
        }
      else if ((*p)->memUsage() != 0)
        {
          string entry = "Pool";
          sprintf(temp, "%11s %-40s %12ld b  %11.4f Mb  peak=%11.4f\n", entry.data(), (*p)->m_name_,
                  (*p)->memUsage(), (Real)(*p)->memUsage()/(Real)BYTES_PER_MEG,
                  (Real)(*p)->memUsage()/(Real)BYTES_PER_MEG);
          a_os << temp;
          totalSize += (*p)->memUsage();
        }
    }

  if (Arena::arenaList_ != NULL)
    {
      for (a = Arena::arenaList_->begin(); a != Arena::arenaList_->end(); a++)
        {
          peak += (*a)->peak;

          if ((*a)->bytes != 0)
            {
              string entry = "BaseFab";
              sprintf(temp, "%11s %-40s %12ld b  %11.4f Mb  peak=%11.4f\n", entry.data(), (*a)->name_,
                      (*a)->bytes, (Real)(*a)->bytes/(Real)BYTES_PER_MEG,
                      (Real)(*a)->peak/(Real)BYTES_PER_MEG);
              a_os << temp;
              totalSize += (*a)->bytes;
            }
        }
    }

  peak += BitSet::peak;
#if 0
  // (DFM -- 10/22/08) -- for now, comment out memtracking of
  // IntVectSet because sizeof(IntVectSet) is dimensionally dependent,
  // so we need to figure out how to keep track of it in a
  //dimensionally-independent way
  peak += IntVectSet::peakcount*sizeof(IntVectSet);
#endif

  if (BitSet::bytes != 0)
    {
      string entry = "BitSet";
      string subentry = "allocation";
      sprintf(temp, "%11s %-40s %12ld b  %11.4f Mb  peak=%11.4f\n", entry.data(), subentry.data(),
              BitSet::bytes, (Real)BitSet::bytes/(Real)BYTES_PER_MEG,
              (Real)BitSet::peak/(Real)BYTES_PER_MEG);
      a_os << temp;
      totalSize += BitSet::bytes;
    }

#if 0
  // (DFM -- 10/22/08) -- for now, comment out memtracking of
  // IntVectSet because sizeof(IntVectSet) is dimensionally dependent,
  // so we need to figure out how to keep track of it in a
  //dimensionally-independent way
  if (IntVectSet::count != 0)
    {
      string entry = "IntVectSet";
      string subentry = "allocation";
      sprintf(temp, "%11s %-40s %12ld b  %11.4f Mb  peak=%11.4f\n", entry.data(), subentry.data(),
              IntVectSet::count * sizeof(IntVectSet), (Real)(IntVectSet::count * sizeof(IntVectSet))/(Real)BYTES_PER_MEG,
              (Real)(IntVectSet::peakcount*sizeof(IntVectSet))/(Real)BYTES_PER_MEG);
      a_os << temp;
      totalSize += IntVectSet::count*sizeof(IntVectSet);
    }
#endif

  //a_os << " ------------------------------------------------------------------------------------------------" << endl;
  string entry = "Total";
  string subentry = "Allocated:";
  sprintf(temp, "%11s %-40s %12ld b  %11.4f Mb  peak=%11.4f\n", entry.data(), subentry.data(),
          totalSize, (Real)totalSize/(Real)BYTES_PER_MEG,
          (Real)peak/(Real)BYTES_PER_MEG);
  a_os << temp;

}

void ReportMemoryUsage(ostream& a_os)
{
  ReportUnfreedMemory(a_os);
}

void UnfreedMemory()
{
  ReportUnfreedMemory(pout());
}

void MemoryUsage()
{
  ReportUnfreedMemory(pout());
}

void memTrackingOn()
{
  trackingOn = true;
}

void memtrackingOff()
{
  trackingOn = false;
}

void overallMemoryUsage(long long& a_currentTotal,
                        long long& a_peak)
{
  a_currentTotal = 0;
  a_peak         = 0;

  AllocList::iterator  i;
  PoolList::iterator   p;
  ArenaList::iterator  a;
  VectorList::iterator v;

  if (trackingOn && allocList)
  {
    for (i = allocList->begin(); i != allocList->end(); i++)
    {
      if ((*i)->total != 0)
      {
        a_currentTotal += (*i)->total;
        a_peak += (*i)->peak;
      }
    }
  }

  if (vectorList_ != NULL)
  {
    for (v = vectorList_->begin(); v != vectorList_->end(); v++)
    {
      long int* bytes = (v->second.first);

      a_currentTotal += *bytes;

      // little trick might break some day.  bvs
      a_peak += *(v->second.second);
    }
  }

  for (p = Pool::m_poolList_->begin(); p != Pool::m_poolList_->end(); p++)
  {
    if (*p == &chunkPool || *p == &infoPool)
    {
      // don't report on the tracker itself.
    }
    else
    {
      a_currentTotal += (*p)->memUsage();
      a_peak         += (*p)->memUsage();
    }
  }

  if (Arena::arenaList_ != NULL)
  {
    for (a = Arena::arenaList_->begin(); a != Arena::arenaList_->end(); a++)
    {
      a_currentTotal += (*a)->bytes;
      a_peak         += (*a)->peak;
    }
  }

  a_currentTotal += BitSet::bytes;
  a_peak         += BitSet::peak;

#if 0
  // (DFM -- 10/22/08) -- for now, comment out memtracking of
  // IntVectSet because sizeof(IntVectSet) is dimensionally dependent,
  // so we need to figure out how to keep track of it in a
  //dimensionally-independent way
  a_currentTotal += IntVectSet::count     * sizeof(IntVectSet);
  a_peak         += IntVectSet::peakcount * sizeof(IntVectSet);
#endif

}

void overallMemoryUsage()
{
  long long longcurrent;
  long long longpeak;

  overallMemoryUsage(longcurrent, longpeak);

  int current = (int)( (Real)(longcurrent) / (Real)(BYTES_PER_MEG) );
  int peak    = (int)( (Real)(longpeak   ) / (Real)(BYTES_PER_MEG) );

  pout()<<"Total unfreed: "<<current<<" Mb   Peak: "<<peak<<" Mb"<<std::endl;
}

void memtrackStamp(Real& a_current,
                   Real& a_peak)
{
  long long longcurrent;
  long long longpeak;

  overallMemoryUsage(longcurrent, longpeak);

  a_current = (Real)(longcurrent) / (Real)(BYTES_PER_MEG);
  a_peak    = (Real)(longpeak   ) / (Real)(BYTES_PER_MEG);

  return;
}

void AddTrack(void*       a_addr,
              size_t      a_asize,
              const char* a_fname,
              int         a_lnum,
              bool        a_mlloc)
{
  if (!allocList)
  {
    allocList = new(AllocList);
  }

  if (!trackingOn)
  {
    return;
  }

  if (a_addr == NULL)
  {
    return;
  }

  chunk* tmp = (chunk*)chunkPool.getPtr();

  tmp->address = a_addr;
  tmp->size    = a_asize;

  AllocList::iterator i;
  for (i = allocList->begin(); i != allocList->end(); i++)
  {
    ALLOC_INFO& info = *(*i);

    if (info.line == a_lnum)
    {
      if (strncmp(a_fname, info.file, strlen(a_fname)) == 0)
      {
        info.chunks.insert(info.chunks.begin(), tmp);
        info.total += a_asize;

        if (info.total > info.peak)
        {
          info.peak=info.total;
        }

        grandTotal += a_asize;

        if (grandTotal > grandPeak)
        {
          grandPeak = grandTotal;
        }

        return;
      }
    }
  }

  ALLOC_INFO* info;
  info = (ALLOC_INFO*)infoPool.getPtr();
  new (info) ALLOC_INFO;

  info->chunks.insert(info->chunks.begin(), tmp);
  info->total = a_asize;
  info->peak  = a_asize;
  strncpy(info->file, a_fname, CHUNK_FILENAME_LEN-1);
  info->line  = a_lnum;
  info->mlloc = a_mlloc;

  allocList->insert(allocList->begin(), info);
}

void RemoveTrack(void* a_addr,
                 bool  a_mlloc)
{
  if (!trackingOn)
  {
    return;
  }

  if (a_addr == NULL || !allocList)
  {
    return;
  }

  AllocList::iterator i;
  list<chunk*>::iterator c;

  for (i = allocList->begin(); i != allocList->end(); i++)
  {
    for (c = (*i)->chunks.begin(); c != (*i)->chunks.end(); c++)
    {
      if ((*c)->address == a_addr)
      {
        bool m = (*i)->mlloc;
        if (a_mlloc != m)
        {
          if (!m)
          {
            MayDay::Error("memory allocated with operator NEW returned using free");
          }
          else
          {
            MayDay::Error("memory allocated with malloc-family returned using DELETE");
          }
        }

        // At least one compiler/OS/HDF5-version combination was
        // causing segfault in the memory-tracking (MT) code.
        // Normally, one might just turn off MT at compile-time via
        // USE_MT=FALSE (EMT='').  Experimented with the following
        // code snippet and discovered that by re-arranging the order,
        // we could avoid the segfault -- at least for the case of
        // gcc3.4.0+ on linux with HDF5-1.6.2.  Not sure why this
        // happens yet.  Might be masking another problem... (ndk)

#ifdef USE_ORIGINAL_CODE_SQUIRRELY_MT_SEGFAULT
        // original code here:
        (*i)->chunks.remove(*c);
        (*i)->total -= (*c)->size;
        chunkPool.returnPtr(*c);

        if ((*i)->chunks.empty())
        {
          allocList->remove(*i);
          infoPool.returnPtr((*i));
        }
#else
        // new code here:
        (*i)->total -= (*c)->size;
        chunkPool.returnPtr(*c);
        (*i)->chunks.remove(*c);

        if ((*i)->chunks.empty())
        {
          infoPool.returnPtr((*i));
          allocList->remove(*i);
        }
#endif
        return;
      }
    }
  }

  MayDay::Error("free/DELETE called on un-malloced/NEW'd pointer, or double freed");
}

// extra piece of Brian trickier to make sure that mallocp etc.
// do not become recursive and ridiculous.  now we need access
// to the real memory functions.

// void* operator new (size_t      a_size,
//                     char const* a_file,
//                     int         a_line)
// {
//   void *ptr = (void *)malloc(a_size);
//   AddTrack(ptr, a_size, a_file, a_line, false);
//
//   return(ptr);
// }

// void operator delete (void *a_p) throw()
// {
//   RemoveTrack(a_p, false);
//   free(a_p);
// }

// void* operator new[] (size_t      a_size,
//                       char const* a_file,
//                       int         a_line)
// {
//   void *ptr = (void *)malloc(a_size);
//   AddTrack(ptr, a_size, a_file, a_line, a_false);
//
//   return(ptr);
// }

// void operator delete[] (void *a_p) throw()
// {
//   RemoveTrack(a_p, false);
//   free(a_p);
// }

#undef malloc
#undef realloc
#undef calloc
#undef free

void* mallocp(size_t      a_size,
              const char* a_file,
              int         a_line)
{
  void* ptr = (void *)malloc(a_size);
  AddTrack(ptr, a_size, a_file, a_line, true);

  return(ptr);
}

void* reallocp(void*       a_p,
               size_t      a_size,
               const char* a_file,
               int         a_line)
{
  RemoveTrack(a_p, true);

  void* ptr = (void *)realloc(a_p,a_size);
  AddTrack(ptr, a_size, a_file, a_line, true);

  return(ptr);
}

void* callocp(size_t a_nelem, size_t a_elsize, const char* a_file, int a_line)
{
  void* ptr = (void *)calloc(a_nelem, a_elsize);
  AddTrack(ptr, a_elsize * a_nelem, a_file, a_line, true);

  return(ptr);
}

void freep(void* a_p)
{
  RemoveTrack(a_p, true);
  free(a_p);
}

inline /*static*/ void
Memtrack::ReportUnfreedMemory(ostream& a_os)
{
  CH_XD::ReportUnfreedMemory(a_os);
}

inline /*static*/ void
Memtrack::UnfreedMemory()
{
  CH_XD::UnfreedMemory();
}

inline /*static*/ void
Memtrack::memTrackingOn()
{
  CH_XD::memTrackingOn();
}

inline /*static*/ void
Memtrack::memtrackingOff()
{
  CH_XD::memtrackingOff();
}

inline /*static*/ void
Memtrack::overallMemoryUsage(long long& a_currentTotal,
                             long long& a_peak)
{
  CH_XD::overallMemoryUsage(a_currentTotal, a_peak);
}

#include "BaseNamespaceFooter.H"

#else

#include "MayDay.H"
#include "BaseNamespaceHeader.H"

void dumpmemoryatexit()
{
  MayDay::Warning("Memory tracking not on; try compiling with USE_MT=TRUE");
}

void UnfreedMemory()
{
  MayDay::Warning("Memory tracking not on; try compiling with USE_MT=TRUE");
}

void MemoryUsage()
{
  MayDay::Warning("Memory tracking not on; try compiling with USE_MT=TRUE");
}

#include "BaseNamespaceFooter.H"

#endif // CH_USE_MEMORY_TRACKING
