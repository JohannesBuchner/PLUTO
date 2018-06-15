#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MEMUSAGE_H_
#define _MEMUSAGE_H_

#include "REAL.H"
#include "BaseNamespaceHeader.H"

///
/**
Print a line of memory information to pout() with a string label
 */
void print_memory_line(const char *);

///
/**
Print a line of memory information to pout()
 */
void print_memory_line(void);

///
/**
 If Linux, uses /proc/self/status to find the peak residentSetSize of process
 and the peak VM size.   Units should be MB.
 Returns 0 if not Linux.
 */
void getPeakMemoryFromOS(Real& VmPeak, Real& VmHWM);

///
/**
 Return the residentSetSize of process from either /proc/self/statm or getrusage(RUSAGE_SELF, &rus)
 Units should be MB.
 If /proc/self/statm is used (Linux only), then VM size is also returned.
 */
void getMemoryUsageFromOS(Real& residentSetSize, Real& size);

///
/**
 Maintain backward compatibility in code.
 Return the residentSetSize of process from either /proc/self/statm or getrusage(RUSAGE_SELF, &rus)
 Units should be MB.
 */
Real get_memory_usage_from_OS(void);

/**
verison that does not convert things to floating-point values
*/

void getMemoryUsageSize(unsigned int& residentSize, unsigned int& size); 
///
/**
 */
inline unsigned int getMemorySize(void)
{
  unsigned int r, s;
  getMemoryUsageSize(r, s);
  return r;
}

void gather_memory_from_procs(Real end_memory,
                              Real &avg_memory,
                              Real &min_memory,
                              Real &max_memory);

void reduce_print_avg_min_max(const char* s, Real value);
int reduce_avg_min_max(Real value, Real &avg, Real &min, Real &max);
int reduce_avg_min_max_loc(Real value, Real& avg, Real& min, Real& max, int& minloc, int& maxloc);

#include "BaseNamespaceFooter.H"

#endif
