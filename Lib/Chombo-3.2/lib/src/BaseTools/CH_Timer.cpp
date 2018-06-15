#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef CH_NTIMER

#include "CH_Timer.H"

#include <iostream>
#include "memtrack.H"
#include "memusage.H"
#include <fstream>
#include <set>
#include <vector>
#include <cstdio>
#include "CH_assert.H"
#include "parstream.H"
#include <cstring>
#ifndef CH_DISABLE_SIGNALS
#include <unistd.h>
#include <csignal>
#endif

using namespace std;

#include "SPMD.H"

#include "BaseNamespaceHeader.H"

struct elem
{
  const TraceTimer* val;

  unsigned long long int time;

  elem(const TraceTimer* p, unsigned long long int t)
    :val(p),
     time(t)
  {
  }

  bool operator < (const elem& rhs) const
  {
    if (val->isPruned())
    {
      return false;
    }
    return time > rhs.time;
  }

  static void buildList(List<elem>& tlist, const TraceTimer& timer);
};

List<elem> tracerlist(false);
std::vector<TraceTimer*> TraceTimer::s_roots;
std::vector<TraceTimer*> TraceTimer::s_currentTimer;
long long int TraceTimer::s_peak = 0;
TraceTimer*  TraceTimer::s_peakTimer = NULL;
bool TraceTimer::s_memorySampling = false;
bool TraceTimer::s_tracing = false;

static int s_depth = TraceTimer::initializer();

double zeroTime = 0;
unsigned long long int zeroTicks = 0;
double secondspertick = 0;

void elem::buildList(List<elem>& tlist, const TraceTimer& timer)
{
  tlist.append(elem(&timer, timer.time()));
  const std::vector<TraceTimer*>& children = timer.children();
  for (int i=0; i<children.size(); i++)
    {
      buildList(tlist, *(children[i]));
    }
}

FILE* sampleFile;
bool sampleFileOpen = false;
int sampleFrequency = 0;

sig_atomic_t samplingOn = false;

void sampleMem(int a_sig)
{
  if (samplingOn)
    TraceTimer::sampleMemUsage2();
}

const char* currentTimer()
{
#ifdef _OPENMP
  if(onThread0())
    return TraceTimer::currentTimer();
  else
    return 0;
#else
    return TraceTimer::currentTimer();
#endif
}

unsigned int peak = 0;

void TraceTimer::sampleMemUsage2()
{
  //AD: This should be OK to run even in threaded scenario
  // whereever samplingOn is turned on it should check for thread0
  if (samplingOn) 
    samplingOn=false; // mutex to protect from double entry
  else 
    return;
  unsigned int m = getMemorySize();
  TraceTimer* current =  s_currentTimer[0];
  if (m>peak)
    {
      peak = m;
      s_peakTimer = current;
    }
  current->m_memoryMin = std::min(m,current->m_memoryMin);
  current->m_memoryMax = std::max(m,current->m_memoryMax);
  samplingOn=true;
}

int sampleCount = 0;
void TraceTimer::sampleMemUsage(const char* name)
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  if (sampleFileOpen)
    {
      unsigned int m = getMemorySize(); //from OS
    //  long long int current, peak;
#ifdef CH_USE_MEMORY_TRACKING
     // overallMemoryUsage(current, peak); // from Chombo Memory Tracking
#else
     //// current = 0;
    //  peak = 0;
#endif
    //  current/= 1024;
     // fprintf(sampleFile, "%d %lld %s\n", m, current, name);
      fprintf(sampleFile, "%d  %s\n", m, name);
      sampleCount++;
      if (sampleCount >= 5000)
      {
         fflush(sampleFile); sampleCount=0;
      }
    }
  else
    {
#ifndef CH_MPI
      sampleFile = fopen("memory.trace", "w");
      sampleFileOpen = true;
#else
      int flag_i, flag_f;
      MPI_Initialized(&flag_i);
      MPI_Finalized(&flag_f);
      if (flag_i)
        {
          char b[1024];
          int outInterv = 1;
          char* charInterv = getenv("CH_OUTPUT_INTERVAL");
          if (charInterv != NULL)
            {
              outInterv =  atoi(charInterv);
              // If zero specified, change it to numProc() which should give time.table.0 only
              if (outInterv == 0) outInterv=numProc();
            }
          int thisProc = procID();
          if ((thisProc % outInterv) != 0)
            {
              sprintf(b,"/dev/null");
            }
          else
            {
              sprintf(b, "memory.trace.%d",procID());
            }
          
          sampleFile = fopen(b, "w");
          fprintf(sampleFile, "@ s%d legend \"%d\"\n",procID(), procID());
          fprintf(sampleFile, "@ target G0.S%d\n", procID());
          sampleFileOpen = true;
        }
#endif
      sampleMemUsage(name);
    }
#ifdef _OPENMP
  }
#endif
}

void writeOnAbort(int sig)
{
  // AD: Deliberately leaving this open to any thread
   TraceTimer::report(true);
}

void writeOnExit()
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  TraceTimer::report(true);
#ifdef _OPENMP
  }
#endif
}

int TraceTimer::initializer()
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  static bool initialized = false;
  if (initialized) return -11;

#ifndef CH_NTIMER
  const char* rootName = "main";
  TraceTimer* rootTimer = new TraceTimer(rootName, NULL, 0);
  rootTimer->m_thread_id = 0;
  char mutex = 0;
  s_roots.resize(1);
  s_roots[0] = rootTimer;
  s_currentTimer.resize(1);
  s_currentTimer[0]=rootTimer;

  char* timerEnv = getenv("CH_TIMER");
  s_memorySampling = false;
  s_tracing= false;
  rootTimer->start(&mutex);
  zeroTime = TimerGetTimeStampWC();
  zeroTicks = ch_ticks();

  if (timerEnv == NULL)
    {
      rootTimer->m_pruned = true;
    }
  else if (strncmp(timerEnv, "SAMPLE=",7)==0)
    {
#ifdef CH_USE_MEMORY_TRACKING
#ifndef CH_DISABLE_SIGNALS
      sampleFrequency = atoi(timerEnv+7);
      signal(SIGALRM, sampleMem);
      ualarm( sampleFrequency, sampleFrequency );
      s_memorySampling = true;
#else 
      std::cout<<"Chombo was compiled with memory tracking, but CH_DISABLE_SIGNALS was set.  no samples created"<<std::endl;
#endif
#else
      std::cout<<"CH_USE_MEMORY_TRACKING was not used during compilation, so no memory sampling will happen"<<std::endl;
#endif
      
    } else if (strncmp(timerEnv, "TRACING",7)==0)
    {
#ifdef CH_USE_MEMORY_TRACKING
#ifndef CH_DISABLE_SIGNALS
      s_tracing = true;
#else 
      std::cout<<"Chombo was compiled with memory tracking, but CH_DISABLE_SIGNALS was set.  no samples created"<<std::endl;
#endif
#else
      std::cout<<"CH_USE_MEMORY_TRACKING was not used during compilation, so no memory sampling will happen"<<std::endl;
#endif
    }

  //  OK, I think I have the atexit vs. static objects bug under AIX worked out. we'll see.
  //#ifndef CH_AIX
  // petermc, 21 April 2006:
  // put "#ifndef CH_AIX" around this because on seaborg,
  // the presence of this line causes a segfault at the termination
  // of the program.
  //   OK, last time around the maypole.  It just seems that atexit and MPI_Finalize are
  //   not going to be cooperative.  bvs
//#ifndef CH_MPI
  if (timerEnv != NULL)
   { 
      atexit(writeOnExit);
#ifndef CH_DISABLE_SIGNALS
      signal(SIGABRT, writeOnAbort);
#endif
   }
//#endif

  //#endif
#endif // CH_NTIMER
  initialized = true;
  if (s_memorySampling)
    samplingOn = true;
#ifdef _OPENMP
  }
#endif
  return 0;
}

void normalizeMemory(unsigned int a_m, int& a_memory, char* units)
{ 
#ifdef _OPENMP
  if(onThread0()){
#endif
  int megabytes = a_m/(1024);
  if (megabytes > 15 )
    {
      strcpy(units, "M");
      a_memory = megabytes;
    }
  else
    {
      strcpy(units, "k");
      a_memory = a_m;
    }
#ifdef _OPENMP
  }
#endif
}

void TraceTimer::currentize() const
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  if (m_pruned) return;

  if (m_last_WCtime_stamp != 0)
  {
    for (int i=0; i<m_children.size(); i++)
      {
        m_children[i]->currentize();
      }
    unsigned long long int current = ch_ticks();
    (unsigned long long int&)m_accumulated_WCtime += current - m_last_WCtime_stamp;
    (unsigned long long int&)m_last_WCtime_stamp = current;

  }
#ifdef _OPENMP
  }
#endif
}

int TraceTimer::computeRank() const
{
#ifdef _OPENMP
  if(onThread0()){
#endif

  tracerlist.clear();
  elem::buildList(tracerlist, *this);
  tracerlist.sort();
  int r=0;
  ListIterator<elem> it(tracerlist);
  for (it.begin(); it.ok(); ++it)
    {
      const elem& e = *it;
      //printf("%s %e %d\n",e.val->m_name, e.time, r);
      e.val->m_rank = r;
      ++r;
    }
  return r;
#ifdef _OPENMP
  }
  return 0;
#endif
}

const TraceTimer* TraceTimer::activeChild() const
{
#ifdef _OPENMP
  if(onThread0()){
#endif

  TraceTimer* child;
  for (int i=0; i<m_children.size(); i++)
    {
      child = m_children[i];
      if (child->m_last_WCtime_stamp != 0) return child;
    }
#ifdef _OPENMP
  }
#endif
  return NULL;

}

const std::vector<TraceTimer*>& TraceTimer::children() const
{
  return m_children; 
}

void TraceTimer::report(bool a_closeAfter)
{

#ifndef CH_NTIMER
#ifdef _OPENMP
  if(onThread0()){
#endif
  if (s_memorySampling) samplingOn = false; //disable the sampling while creating a report
  if (s_tracing)
    {
      fflush(sampleFile);
    }
  char* timerEnv = getenv("CH_TIMER");
  if (timerEnv == NULL)
    {
      // pout()<<"CH_TIMER environment variable not set. Timers inactive. Not writing time.table \n";
      return;
    }

  TraceTimer& root = *(s_roots[0]); // in MThread code, loop over roots
  root.currentize();
  int numCounters = root.computeRank();

  double elapsedTime = TimerGetTimeStampWC() - zeroTime;
  unsigned long long int elapsedTicks = ch_ticks() - zeroTicks;
  secondspertick = elapsedTime/(double)elapsedTicks;

  int mpirank = 0;
#ifdef CH_MPI
  int proc = getpid();
  int finalized;
  MPI_Finalized(&finalized);
  if (finalized)
    mpirank = GetRank(proc);
  else
    mpirank = procID();
#endif

  if (mpirank >= 0)
    {
      char buf[1024];
#ifdef CH_MPI
      int outInterv = 1;
      char* charInterv = getenv("CH_OUTPUT_INTERVAL");
      if (charInterv != NULL)
        {
          outInterv =  atoi(charInterv);
          // If zero specified, change it to numProc() which should give time.table.0 only
          if (outInterv == 0) outInterv=numProc();
        }

      int thisProc = procID();
      if ((thisProc % outInterv) != 0)
        {
          sprintf(buf,"/dev/null");
        }
      else
        {
          sprintf(buf,"time.table.%d",mpirank);
        }
#else
      sprintf(buf,"time.table");
#endif
      static FILE* out = fopen(buf, "w");
      static int reportCount = 0;
      fprintf(out, "-----------\nTimer report %d (%d timers)\n--------------\n",
              reportCount, numCounters);
      reportCount++;
      if (s_memorySampling) updateMemory(root);
      ListIterator<elem> it(tracerlist);
      if (s_memorySampling) reportPeak(out);
      for (it.begin(); it.ok(); ++it)
        reportOneTree(out, *((*it).val));
      subReport(out, "FORT_", root.m_accumulated_WCtime );
      subReport(out, "MPI_", root.m_accumulated_WCtime );

      reportFullTree(out, root, root.m_accumulated_WCtime, 0); //uses recursion
      fflush(out);
      if (a_closeAfter) fclose(out);
    }

  if (s_memorySampling && !a_closeAfter) samplingOn = true; // enable sampling again.....
#endif
#ifdef _OPENMP
  }
#endif
}

void TraceTimer::reset()
{
#ifdef _OPENMP
  if(onThread0()){
#endif

  char* timerEnv = getenv("CH_TIMER");
  if (timerEnv == NULL)
    {
      // pout()<<"CH_TIMER environment variable not set. Timers inactive. Not writing time.table \n";
      return;
    }
  TraceTimer& root = *(s_roots[0]);
  root.currentize();
  reset(root);
#ifdef _OPENMP
  }
#endif
}

void TraceTimer::reset(TraceTimer& node)
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  node.m_count = 0;
  node.m_accumulated_WCtime = 0;
  for (int i=0; i<node.m_children.size(); i++)
    {
      reset(*(node.m_children[i]));
    }
#ifdef _OPENMP
  }
#endif
}

void sorterHelper(const std::vector<TraceTimer*>& children, Vector<int>& order)
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  int n = children.size();
  order.resize(n);
  for (int i=0; i<n; ++i) order[i]=i;
  bool swaps = true;
  while (swaps)
    {
      swaps = false;
      for (int i=0; i<n-1; ++i)
      {
        if (children[order[i]]->time()  < children[order[i+1]]->time())
          {
            int tmp = order[i];
            order[i] = order[i+1];
            order[i+1] = tmp;
            swaps = true;
            break;
          }
      }
    }
#ifdef _OPENMP
  }
#endif
}

void TraceTimer::subReport(FILE* out, const char* header, unsigned long long int totalTime)
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  size_t length = strlen(header);
  fprintf(out, "=======================================================\n");
  unsigned long long int subTime = 0;
  ListIterator<elem> it(tracerlist);
  for (it.begin(); it.ok(); ++it)
    {
      const char* name = (*it).val->m_name;
      if (strncmp(header, name, length) == 0)
      {
        if ((*it).val->isPruned())
        {
          //fprintf(out, "             pruned  %s  \n", name);

        }
        else
        {
          unsigned long long int t = (*it).val->time();
          int rank = (*it).val->rank();
          subTime += t;
          fprintf(out, "  %8.3f %8lld  %s [%d] \n", t*secondspertick, (*it).val->m_count, name, rank);
        }
      }
    }
  if (subTime > 0)
    fprintf(out, "  %8.3f   %4.1f%%    Total\n", subTime*secondspertick, (double)subTime/totalTime*100.0);
#ifdef _OPENMP
  }
#endif
}

void TraceTimer::updateMemory(TraceTimer& a_timer)
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  if (a_timer.m_pruned) return;

  for (int i=0; i<a_timer.m_children.size(); ++i)
    {
      TraceTimer& child = *(a_timer.m_children[i]);
      updateMemory(child);
      a_timer.m_memoryMin = std::min(child.m_memoryMin, a_timer.m_memoryMin);
      a_timer.m_memoryMax = std::max(child.m_memoryMax, a_timer.m_memoryMax);
    }
#ifdef _OPENMP
  }
#endif
}

void TraceTimer::reportFullTree(FILE* out, const TraceTimer& timer,
                                unsigned long long int totalTime, int depth)
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  if (timer.m_pruned) return;
  unsigned long long int time = timer.m_accumulated_WCtime;

  if (depth < 20)
  {
    for (int i=0; i<depth; ++i) fprintf(out,"   ");
    double percent = ((double)time)/totalTime * 100.0;
    fprintf(out, "[%d] %s %.4f %4.1f%% %lld \n", timer.m_rank, timer.m_name, time*secondspertick, percent, timer.m_count);
  }
  Vector<int> ordering;
  sorterHelper(timer.m_children, ordering);
  for (int i=0; i<timer.m_children.size(); ++i)
  {
    reportFullTree(out, *(timer.m_children[ordering[i]]), totalTime, depth+1);
  }
#ifdef _OPENMP
  }
#endif

}
void TraceTimer::reportOneTree(FILE* out, const TraceTimer& timer)
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  if (timer.m_pruned) return;
  unsigned long long int time = timer.m_accumulated_WCtime;
  unsigned long long int subTime = 0;

  fprintf(out,"---------------------------------------------------------\n");

  fprintf(out,"[%d]%s %.2f %lld", timer.m_rank, timer.m_name, time*secondspertick, timer.m_count);
  if (s_memorySampling && timer.m_memoryMax > 0)
    {
      int t_min, t_max;
      char units_min[2], units_max[2];
      normalizeMemory(timer.m_memoryMin, t_min, units_min);
      normalizeMemory(timer.m_memoryMax, t_max, units_max);
      fprintf(out," %5d%s  %5d%s \n", t_min, units_min, t_max, units_max);
    }
  else
    {
      fprintf(out, "\n");
    }
  const std::vector<TraceTimer*>& children = timer.m_children;
  Vector<int> ordering;
  sorterHelper(children, ordering);
  for (int i=0; i<children.size(); ++i)
    {
      const TraceTimer& child = *(children[ordering[i]]);
      if (!child.m_pruned)
        {
          unsigned long long int childtime = child.m_accumulated_WCtime;
          if (childtime > 0)
          {
            subTime += childtime;
            double percent = ((double)childtime) / time * 100.0;
            if (!s_memorySampling || child.m_memoryMax == 0)
              fprintf(out,"    %4.1f%% %7.2f %8lld %s [%d]\n",
                      percent, childtime*secondspertick, child.m_count, child.m_name, child.m_rank);
            else
              {
                int t_min, t_max;
                char units_min[2], units_max[2];
                normalizeMemory(child.m_memoryMin, t_min, units_min);
                normalizeMemory(child.m_memoryMax, t_max, units_max);
                fprintf(out,"    %4.1f%% %7.2f %8lld %s [%d]  %5d%s  %5d%s\n",
                        percent, childtime*secondspertick, child.m_count, 
                        child.m_name, child.m_rank, t_min, units_min, t_max, units_max);
              }
          }
        }
      else
      {
         fprintf(out,"           pruned           \n");
         i=children.size();
      }
    }
  if (time > 0 && children.size() > 0)
  {
    double totalPercent = ((double)subTime)/ time * 100.0;
    fprintf(out, "    %4.1f%%                  Total \n", totalPercent);
  }
#ifdef _OPENMP
  }
#endif
}

bool TraceTimer::find(std::list<const TraceTimer*>& trace, const TraceTimer* target)
{
#ifdef _OPENMP
  if(onThread0()){
#endif

  const TraceTimer* current = trace.back();
  if (current == target) return true;
  for (int i=0; i<current->m_children.size(); i++)
    {
      trace.push_back(current->m_children[i]);     
      bool found = find(trace, target);
      if (found) return true;
      trace.pop_back();
    }
#ifdef _OPENMP
  }
#endif
  return false;
}

void TraceTimer::reportPeak(FILE* out)
{
#ifdef _OPENMP
  if(onThread0()){
#endif

  fprintf(out,"peak memory %d K at \n", peak);
  const TraceTimer* current = s_roots[0];
  // first, have to search from root to leaves until I find the peakTimer, then I can retrace the list
  std::list<const TraceTimer*> trace;
  trace.push_back(current);
  find(trace, s_peakTimer);

  int depth = 0;
  std::list<const TraceTimer*>::const_iterator i = trace.begin();
  while (i != trace.end()) 
    {
      for (int d=0; d<depth; d++) fprintf(out," ");
      fprintf(out, "%s  \n",(*i)->m_name);
      ++depth;
      ++i;
    }
#ifdef _OPENMP
  }
#endif
}
  
/*
void TraceTimer::reportMemoryOneTree(FILE* out, const TraceTimer& timer)
{
  if (timer.m_pruned) return;
  long long int m = timer.m_memory;
  long long int p = timer.m_peak;

  if (p==0 && m==0) return;

  char units[2], punits[2];

  int  memory, peak;

  normalizeMemory(m, memory, units);
  normalizeMemory(p, peak, punits);

  fprintf(out,"---------------------------------------------------------\n");

  fprintf(out,"[%d]%s %5d%s  %5d%s %lld\n", timer.m_rank, timer.m_name, peak, punits,
          memory, units, timer.m_count);

  const std::vector<TraceTimer*>& children = timer.m_children;
  Vector<int> ordering;
  sorterHelper(children, ordering);
  for (int i=0; i<children.size(); ++i)
    {
      const TraceTimer& child = *(children[ordering[i]]);
      if (!child.m_pruned)
        {
          long long int pm = child.m_peak;
          long long int sm = child.m_memory;

          if (sm != 0 || pm != 0)
          {
            normalizeMemory(sm, memory, units);
            normalizeMemory(pm, peak, punits);
            fprintf(out,"   %5d%s   %5d%s %8lld %s [%d]\n",
                    peak, punits, memory, units, child.m_count, child.m_name, child.m_rank);
          }
        }
      else
      {
         fprintf(out,"           pruned           \n");
         i=children.size();
      }
    }
}
*/

// some compilers complain if there isn't at least 1 non-inlined
// function for every class.  These two are mostly to make those compilers happy

char AutoStart::ok = 0;

bool AutoStart::active()
{
  return true;
}

bool AutoStartLeaf::active()
{
  return true;
}

TraceTimer::~TraceTimer()
{
#ifdef _OPENMP
  if(onThread0()){
#endif

  for (int i=0; i<m_children.size(); ++i)
    {
      delete m_children[i];
    }
#ifdef _OPENMP
  }
#endif
}

void TraceTimer::reportName(std::ostream& out, const char* name)
{
#ifdef _OPENMP
  if(onThread0()){
#endif

  double elapsedTime = TimerGetTimeStampWC() - zeroTime;
  unsigned long long int elapsedTicks = ch_ticks() - zeroTicks;
  secondspertick = elapsedTime/(double)elapsedTicks;
  const TraceTimer* timer = TraceTimer::getTimer(name);
  unsigned long long int time = timer->m_accumulated_WCtime;
  out<< timer->m_name <<" "<<time*secondspertick <<" "<<timer->m_count;
#ifdef _OPENMP
  }
#endif
}
  
TraceTimer* TraceTimer::getTimer(const char* name)
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  int thread_id = 0; // this line will change in MThread-aware code.
  TraceTimer* parent = TraceTimer::s_currentTimer[thread_id];
  if (parent->m_pruned) return parent;
  std::vector<TraceTimer*>& children = parent->m_children;
  int i=0;
  for (; i<children.size(); ++i)
  {
    TraceTimer* timer =  children[i];
    if (timer->m_name == name) return timer;
  }
  TraceTimer* newTimer = new TraceTimer(name, parent, thread_id);
  children.push_back(newTimer);
  return newTimer;
#ifdef _OPENMP
  }
  return NULL;
#endif
}

#ifdef _OPENMP
TraceTimer::TraceTimer(const char* a_name, TraceTimer* parent, int thread_id)
{
  if(onThread0())
    {
      m_pruned =false;
      m_parent = parent; 
      m_name   = a_name;
      m_count  = 0;
      m_accumulated_WCtime = 0;
      m_last_WCtime_stamp = 0;
      m_thread_id = thread_id;
      //m_memory = 0;
      //m_peak = 0;
    }
}
#else
TraceTimer::TraceTimer(const char* a_name, TraceTimer* parent, int thread_id)
  :m_pruned(false), m_parent(parent), m_name(a_name), m_count(0),
   m_accumulated_WCtime(0),m_last_WCtime_stamp(0), m_thread_id(thread_id),
   m_memoryMin(0), m_memoryMax(0)
{
  m_memoryMin--;  // roll it back to the largest possible value;
}
#endif

void TraceTimer::prune()
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  int i=0;
  for (; i<m_children.size(); ++i)
  {
    TraceTimer* timer =  m_children[i];
    timer->prune();
  }
  m_pruned = true;
#ifdef _OPENMP
  }
#endif
}

void TraceTimer::start(char* mutex)
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  if (m_pruned) return;
# ifndef NDEBUG
  if (*mutex == 1)
  {
    char buf[1024];
    sprintf(buf, "double TraceTimer::start called: %s ",m_name);
    MayDay::Error(buf);
  }
# endif
  ++m_count;
  *mutex = 1;
  s_currentTimer[m_thread_id] = this;
  if (s_memorySampling) //here's to hoping a two-bit branch predictor gets this right
    {
      unsigned int m = getMemorySize();
      m_memoryMin = std::min(m_memoryMin, m);
      m_memoryMax = std::max(m_memoryMax, m);
    }
  if (s_tracing)
    {
      sampleMemUsage(m_name);
    }
  m_last_WCtime_stamp = ch_ticks();

#ifdef _OPENMP
  }
#endif
}
unsigned long long int overflowLong = (unsigned long long int)1<<50;
unsigned long long int TraceTimer::stop(char* mutex)
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  if (m_pruned) return 0;
#ifndef NDEBUG
  if (s_currentTimer[0] != this)
    {
    char buf[1024];
    sprintf(buf, "TraceTimer::stop called while not parent: %s ",m_name);
    MayDay::Error(buf);
    }
#endif
  unsigned long long int diff = ch_ticks();
  diff -= m_last_WCtime_stamp;
  if (diff > overflowLong) diff = 0;
  m_accumulated_WCtime += diff;

  if (s_memorySampling) //here's to hoping a two-bit branch predictor gets this right
    {
      unsigned int m = getMemorySize();
      m_memoryMin = std::min(m_memoryMin, m);
      m_memoryMax = std::max(m_memoryMax, m);
    }
  if (s_tracing)
    {
      sampleMemUsage("end");
    }
  m_last_WCtime_stamp=0;
  s_currentTimer[m_thread_id] = m_parent;
  *mutex=0;

  return diff;
#ifdef _OPENMP
  }
  return 0;
#endif
}

void TraceTimer::macroTest2()
{
  CH_TIME("macroTest2");
}
void TraceTimer::macroTest()
{
  CH_TIMERS("testy");
  CH_TIMER("billy", t1);
  CH_TIMER("sally", t2);
  CH_TIMER("renaldo", billy);
  CH_START(t1);
  CH_STOP(t1);
  CH_START(t2);
  CH_STOP(t2);
  CH_START(billy);
  CH_STOP(billy);
  CH_TIMER_REPORT();
  CH_TIMER_RESET();
  CH_TIMER_PRUNE(0.01);
}

void TraceTimer::PruneTimersParentChildPercent(double threshold, TraceTimer* parent)
{
#ifdef _OPENMP
  if(onThread0()){
#endif
  if (parent->isPruned()) return;
  unsigned long long int time = parent->time();
  const std::vector<TraceTimer*>& children = parent->children();

  for (int i=0; i<children.size(); ++i)
    {
      TraceTimer* child = children[i];
      if (!child->isPruned())
        {
          unsigned long long int childtime = child->time();
          if (((double)childtime)/time < threshold) child->prune();
          else PruneTimersParentChildPercent(threshold, child);
        }

    }
#ifdef _OPENMP
  }
#endif
}

void TraceTimer::PruneTimersParentChildPercent(double percent)
{
#ifdef _OPENMP
  if(onThread0()){
#endif
#ifndef CH_NTIMER
  char* timerEnv = getenv("CH_TIMER");
  if (timerEnv == NULL)
    {
      // pout()<<"CH_TIMER environment variable not set. Timers inactive. Not writing time.table \n";
      return;
    }
  TraceTimer* root = s_roots[0]; // in MThread code, loop over roots
  root->currentize();
  PruneTimersParentChildPercent(percent, root);
#endif
#ifdef _OPENMP
  }
#endif
}

#else // on CH_NTIMER
const char* currentTimer()
{
  const char* rtn ="Timers not active";
  return rtn;
}

#endif //on CH_NTIMER

#ifndef CH_NTIMER
#include "BaseNamespaceFooter.H"
#endif
