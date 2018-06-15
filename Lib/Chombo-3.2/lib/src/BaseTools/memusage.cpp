#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "memusage.H"
#include "memtrack.H"
#include "parstream.H"
#include "MayDay.H"
#include <cstdio>
#include <cstring>
#include "SPMD.H"
#include <sys/time.h>
#include <sys/resource.h>
#ifdef CH_MPI
#include <mpi.h>
#endif

#include "BaseNamespaceHeader.H"

void print_memory_line(void)
{
  std::string s="";
  print_memory_line(s.data());
}

void print_memory_line(const char *s)
{
  char temp[1024];
#ifdef CH_USE_MEMORY_TRACKING
  Real memtrackCurrentMemory;
  Real memtrackPeakMemory;
  memtrackStamp(memtrackCurrentMemory, memtrackPeakMemory);

  Real residentSetSize=0;
  Real size=0;
  getMemoryUsageFromOS(residentSetSize, size);
#ifdef CH_Linux
  Real peakRSS, peakVirtualMem;
  getPeakMemoryFromOS(peakRSS, peakVirtualMem);
  sprintf(temp, "%26.26s|Mem Usage: OS=%8.3f (vm=%8.3f) MT_peak=%8.3f MT_current=%8.3f  OSPeakRSS=%8.3f (vm=%8.3f) (MB)\n",
          s, residentSetSize, size, memtrackPeakMemory, memtrackCurrentMemory, peakRSS, peakVirtualMem);
#else
  sprintf(temp, "%26.26s|Mem Usage: OS=%8.3f (vm=%8.3f) MT_peak=%8.3f MT_current=%8.3f (MB)\n",
          s, residentSetSize, size, memtrackPeakMemory, memtrackCurrentMemory);
#endif

#else
  sprintf(temp, "%26.26s|Mem Usage: OS=%8.3f (MB) MT is off\n",
          s, get_memory_usage_from_OS());
#endif
  pout() << temp;
}

// PAMM -- Print Average Min Max
//   reduce a value over all procs to obtain average/min/max and the ranks containing min/max
//   Print a line on rank0 containing this information (with a supplied label)
void reduce_print_avg_min_max(const char* s, Real value)
{
  Real avg,min,max;
  int minloc, maxloc;
  reduce_avg_min_max_loc(value, avg, min, max, minloc, maxloc);
  if (procID() == uniqueProc(SerialTask::compute))
  {
    char tmp[1024];
    sprintf(tmp, "PAMM %32.32s|avg=%12.1f min=%12.1f @rank%-6d max=%12.1f @rank%-6d\n",
            s, avg, min, minloc, max, maxloc);
    pout() << tmp;
  }
}

//  For linux machines, This function takes advantage of the Linux file /proc/self/status
// Some of the values in the file are described below:
//    VmPeak: 3141876 kB -- Peak virtual memory size. The largest amount of memory this process has been able to address.
//    VmSize: 3141876 kB -- Virtual memory size. Size of the process's address space less reserved regions
//    VmLck: 0 kB        -- Size of the locked pages (i.e. can't be swapped out).
//    VmHWM: 12556 kB    -- Peak resident set size ("high water mark"). Largest amount memory ever owned by the process.
//    VmRSS: 12556 kB    -- Resident set size. Size of all the page frames owned by the process.
//    VmData: 3140564 kB -- Size of all the process's addressable memory less the shared and stack sections.
//    VmStk: 88 kB       -- Size of this process's stack.
//    VmExe: 4 kB        -- Size of this process's executable code.
//    VmLib: 1204 kB     -- Shared library code size. Size of the executable section not including this process's code.
//                          i.e. the amount of memory mapped to the executable code of libraries.
//    VmPTE: 3072 kB     -- Size of the page tables this process uses.
//
//more details, take a look at the kernel code that displays this information.
//http://lxr.linux.no/linux/fs/proc/task_mmu.c#L40
//
//http://manpages.courier-mta.org/htmlman5/proc.5.html
//
void getPeakMemoryFromOS(Real& peakRSS, Real& peakVirtualMem)
{
  peakRSS=0.0;
  peakVirtualMem=0.0;
#ifdef CH_Linux
#ifndef CH_DISABLE_SIGNALS
  // if Linux(except for catamount), look at the file /proc/self/status
  FILE *fh;
  char proc_var[1024];
  char *cp;
  long long unsigned int ivmpeak=0, ivmhwm=0;

  //pout() << " Reading /proc/self/status to find memory use\n";
  fh = fopen("/proc/self/status","r");
  while (!feof(fh))
  {
    char *ret;
    ret = fgets(proc_var,80,fh);
    //pout() << proc_var;
    cp = strstr(proc_var,"VmPeak:"); // peak virtual mem
    if (cp) sscanf(cp, "VmPeak:"" %llu", &ivmpeak);

    cp = strstr(proc_var,"VmHWM:"); // resident set size
    if (cp) sscanf(cp, "VmHWM:"" %llu", &ivmhwm);
  }
  fclose(fh);
  // ivmpeak and ivmhwm are in KB at this point
  peakRSS = ivmhwm/1024.0;  //  MB
  peakVirtualMem = ivmpeak/1024.0;  //  MB
#endif
#endif
  return;
}

// For linux machines, This function takes advantage of the Linux file /proc/self/statm
//size       total program size
//           (same as VmSize in /proc/[pid]/status)
//resident   resident set size
//           (same as VmRSS in /proc/[pid]/status)
void getMemoryUsageFromOS(Real& residentSetSize, Real& size)
{
  unsigned int r, s;
  getMemoryUsageSize(r, s);
  residentSetSize = r/1024.0;
  size = s/1024.0;
}

void getMemoryUsageSize(unsigned int& residentSetSize, unsigned int& size) 
{
  size=0.0;
  residentSetSize=0.0;
#ifdef CH_Linux
#ifndef CH_DISABLE_SIGNALS
  // if Linux(except for catamount), look at the file /proc/self/statm
  FILE* f=fopen("/proc/self/statm","r");
  //fseek(f,0, SEEK_SET);       
  int isize, iresident, ishared;
  if (fscanf(f, "%d %d %d", &isize, &iresident, &ishared)==3)
    {
      //printf("%10.2f MB total size, %10.2f MB resident, %10.2f MB shared\n",
      //     (isize*4)/1024.0, (iresident*4)/1024.0, (ishared*4)/1024.0);

      // linux on intel x86 has 4KB page size...
      //return ($size * 4, $share * 4);

      size = (isize*4);
      residentSetSize = (iresident*4);
    }
  fclose(f);
#endif

#else
  static struct rusage rus;
  getrusage(RUSAGE_SELF, &rus);

  // ru_maxrss is a long int and is the "Maximum resident set size (in kilobytes)"

  //printf(" ru_maxrss= %d  %10.2fMB\n",
  // rus.ru_maxrss, rus.ru_maxrss/1024.0);
  //printf(" ru_ixrss=%d ru_idrss=%d ru_isrs=%d ru_maxrss=%d\n",
  //     rus.ru_ixrss, rus.ru_idrss, rus.ru_isrss, rus.ru_maxrss);
  // except on apple this is acutally in bytes
  residentSetSize = rus.ru_maxrss;
#ifdef __APPLE__
  residentSetSize/=1024;
#endif
#endif
}

// Maintain backward compatibility in code. Deprecate. Calls the above function.
// Return the residentSetSize of process from either /proc/self/statm or getrusage(RUSAGE_SELF, &rus)
// Units should be MB.
Real get_memory_usage_from_OS(void)
{
  Real residentSetSize=0.0;
  Real size;
  getMemoryUsageFromOS(residentSetSize, size);  // units MB
  return residentSetSize;
}

// Generic function to reduce the min,max,avg of a Real over all procs onto rank0
int reduce_avg_min_max(Real value, Real &avg, Real &min, Real &max)
{
  int eek=0;
#ifdef CH_MPI
  Real sum;
  eek = MPI_Reduce(&value, &sum, 1, MPI_CH_REAL, MPI_SUM, uniqueProc(SerialTask::compute), Chombo_MPI::comm);
  CH_assert(eek == MPI_SUCCESS);
  avg=sum/(Real)numProc();

  eek = MPI_Reduce(&value, &min, 1, MPI_CH_REAL, MPI_MIN, uniqueProc(SerialTask::compute), Chombo_MPI::comm);
  CH_assert(eek == MPI_SUCCESS);

  eek = MPI_Reduce(&value, &max, 1, MPI_CH_REAL, MPI_MAX, uniqueProc(SerialTask::compute), Chombo_MPI::comm);
  CH_assert(eek == MPI_SUCCESS);
#else
  avg=value;
  min=value;
  max=value;
#endif
  return eek;
}

// Generic function to reduce the minloc,maxloc,avg of a Real over all procs onto rank0
int reduce_avg_min_max_loc(Real value, Real& avg, Real& min, Real& max, int& minloc, int& maxloc)
{
  int eek=0;
#ifdef CH_MPI
  struct
  {
    double val;
    int rank;
  } in, out;
  in.val = value;
  in.rank = procID();

  Real sum;
  eek = MPI_Reduce(&value, &sum, 1, MPI_CH_REAL, MPI_SUM, uniqueProc(SerialTask::compute), Chombo_MPI::comm);
  CH_assert(eek == MPI_SUCCESS);
  avg=sum/(Real)numProc();

  eek = MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, uniqueProc(SerialTask::compute), Chombo_MPI::comm);
  CH_assert(eek == MPI_SUCCESS);
  min = out.val;
  minloc = out.rank;

  eek = MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, uniqueProc(SerialTask::compute), Chombo_MPI::comm);
  CH_assert(eek == MPI_SUCCESS);
  max = out.val;
  maxloc = out.rank;
#else
  avg=value;
  min=value;
  max=value;
  minloc=0;
  maxloc=0;
#endif
  return eek;
}

#ifdef CH_MPI
// Maintain backward compatibility in code.
// Can now use reduce_print_avg_min_max("label", value);
void gather_memory_from_procs(Real end_memory,
                              Real &avg_memory,
                              Real &min_memory,
                              Real &max_memory)
{
  reduce_avg_min_max(end_memory, avg_memory, min_memory, max_memory);

  // results only need to be on rank0 (and are)
  if (procID() == 0)
    {
      pout() << "Gather end memory from procs:  avg: " << avg_memory
             << "  min: " << min_memory
             << "  max: " << max_memory << " (MB)\n";
    }
}
#endif

#include "BaseNamespaceFooter.H"
