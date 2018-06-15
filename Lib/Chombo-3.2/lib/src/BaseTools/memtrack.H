#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MEMTRACK_H_
#define _MEMTRACK_H_

#ifdef CH_USE_MEMORY_TRACKING

#include <iostream>
using std::ostream;

#include "REAL.H"

// NOTE:  These were added to avoid macro substitutions of "calloc(a,b)",
// "realloc(a,b)", and "free(a)" for member function(s) of the "Arena" class.
//
// These member function(s) are declared in "Arena.H" and used
// "BaseFabImplem.H" (which is included via "BaseFab.H").
#include "Arena.H"
//#include "BaseFab.H"

// Use these macros to include in memory tracking system
#define callocMT(a_a,a_b)   callocp (a_a, a_b, __FILE__, __LINE__)
#define reallocMT(a_a, a_b) reallocp(a_a, a_b, __FILE__, __LINE__)
#define mallocMT(a_a)       mallocp (a_a,      __FILE__, __LINE__)
#define freeMT(a_a)         freep   (a_a)

#include "BaseNamespaceHeader.H"

///
void dumpmemorymap(int a_sig);

///
void dumpmemoryatexit();

///
void dumpmemoryabort();

///
int registerMemorySignals();

///
void ReportUnfreedMemory(ostream& a_os);

///
void ReportAllocatedMemory(ostream& a_os);

///
void ReportMemoryUsage(ostream& a_os);

/// calls ReportUnfreedMemory with pout()
void UnfreedMemory();

/// calls  ReportMemoryUsage with pout()
void MemoryUsage();

///
void memTrackingOn();

///
void memtrackingOff();

void overallMemoryUsage(long long& a_currentTotal,
                        long long& a_peak);

void overallMemoryUsage();

void memtrackStamp(Real& a_current,
                   Real& a_peak);

/// Memory tracking functions
/**
   Memory tracking functions
*/

class Memtrack
{
public:
  /// Produce a formated output onto os of memory usage.
  /**
     Memory is reported one a class-by-class basis.  When
     CH_USE_MEMORY_TRACKING is turned on you also get a report from
     the atexit() function.  This can be handy for spotting memory
     leaks.  The memory tracking functionality consumes very little
     runtime and does not impede optimized performance.
  */
  static void ReportUnfreedMemory(ostream& a_os);

  /// calls ReportUnfreedMemory with pout()
  static void UnfreedMemory();

  static void memTrackingOn();

  static void memtrackingOff();

  static void overallMemoryUsage(long long& a_currentTotal,
                                 long long& a_peak);
};

// =========== end user functions.===================================

void AddTrack(void*       a_addr,
              size_t      a_asize,
              const char* a_fname,
              int         a_lnum,
              bool        a_malloc);

void RemoveTrack(void* a_addr,
                 bool  a_malloc);

// void* operator new (size_t      a_size,
//                     char const* a_file,
//                     int         a_line);

// void operator delete (void *a_p) throw();

// void* operator new[] (size_t      a_size,
//                       char const* a_file,
//                       int         a_line);

// void operator delete[] (void *a_p) throw();

// Don't call these directly -- use the callocMT/mallocMT/freeMT macros.
// Use these functions to track the source file and line number of particular
// memory alloc/frees.
void* mallocp(size_t      a_size,
              const char* a_file,
              int         a_line);

void* reallocp(void*       a_p,
               size_t      a_size,
               const char* a_file,
               int         a_line);

void* callocp(size_t      a_nelem,
              size_t      a_elsize,
              const char* a_file,
              int         a_line);

void freep(void* a_p);

#include "BaseNamespaceFooter.H"

#else // on CH_USE_MEMORY_TRACKING

// If not compiling with memory tracking, then just call the LIBC versions
#include <cstdlib>
#define callocMT(a_a,a_b)   calloc (a_a, a_b)
#define reallocMT(a_a, a_b) realloc(a_a, a_b)
#define mallocMT(a_a)       malloc (a_a)
#define freeMT(a_a)         free   (a_a)

inline void AddTrack(void*       a_addr,
              size_t      a_asize,
              const char* a_fname,
              int         a_lnum,
              bool        a_malloc)
{
  // null op
}

inline void RemoveTrack(void* a_addr,
                 bool  a_malloc)
{
  // null op
}

#endif // ifdef CH_USE_MEMORY_TRACKING

#endif // include guard
