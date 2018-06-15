#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifdef MEMORY_USAGE
#include "memusage.H"
static char stuff[1024];
#endif

#include "OldTimer.H"

char dummFlag[4];
#ifndef CH_NTIMER

#include <iostream>
#include "memtrack.H"
#include "memusage.H"
#include <fstream>
#include <set>
#include <vector>
#include <cstdio>
#include "CH_assert.H"
#include <cstring>
#ifndef CH_DISABLE_SIGNALS
#include <unistd.h>
#include <csignal>
#endif

using namespace std;

#ifndef CATFISH
#include "parstream.H"
#endif

#include "SPMD.H"

#include "BaseNamespaceHeader.H"

// Must initialize the static list defined in CH_Timer.H
list<OldTimer*> *OldTimer::TimerList = NULL;

static int ID_counter=0;
static double NtotalStartStops = 0.0;

#define TIME_UNIT_FACTOR 1.0

#ifdef PAPI
static int papiEventID=PAPI_NULL;
static float CPU_MHZ = 0.1;

// These are just helper functions that call PAPI functions
//   Only used here and only for making code easier to read.
static int PAPI_initSingleEvent(const int event);
static int PAPI_initDualEvent(const int event1, const int event2);
static void PAPI_initMultiEvent(const int event1, const int event2,
                                const int event3, const int event4);
static int PAPI_checkEvent(const int event);

static int PAPI_init(void);
static void PAPI_hwinfo(const int rank, float &mhz);
static void PAPI_initMultiplex(void);
static void PAPI_cleanup();
double computeDerivedCounter(double pc0, double pc1, double pc2);
#endif

// These are helper functions to write the formatted output
//  summary.  Nothing special.
static void writeDoubleLineSeparator(FILE *out);
static void writeSingleLineSeparator(FILE *out);
static void writeTableHeader(FILE *out);
static void writeLineOfData(FILE *out,
                            const char *name,
                            const double count,
                            const bool   evenCount,
                            const double percent,
                            const double wc_avg,
                            const double wc_min,
                            const double wc_max,
                            const int number_procs);

#ifdef PAPI
static void writeLineOfData(FILE *out,
                            const char *name,
                            const double count,
                            const bool   evenCount,
                            const double percent,
                            const double wc_avg,
                            const double wc_min,
                            const double wc_max,
                            const double counter0,
                            const double counter1,
                            const int number_procs);
static void writeLineOfData(FILE *out,
                            const char *name,
                            const double count,
                            const bool   evenCount,
                            const double percent,
                            const double wc_avg,
                            const double wc_min,
                            const double wc_max,
                            const double counter0,
                            const double counter1,
                            const double counter2,
                            const double counter3,
                            const int number_procs);
#endif
static void writeTableTotals(FILE *out,
                             const int table_number,
                             const double table_count_sum,
                             const double parent_avg,
                             const double table_avg_sum,
                             const double table_min_sum,
                             const double table_max_sum,
                             const double table_pc0_sum,
                             const double table_pc1_sum,
                             const double table_pc2_sum,
                             const double table_pc3_sum,
                             const double TimerCost);

// try
//OldTimer* timers[2];

//timer[0] = new OldTimer(...);
//timer[1] = new OldTimer(...);

// OldTimer construction for the root of a tree.
OldTimer::OldTimer(const string& name, const int tableid)
  :
  m_name(name),
  m_Parent(*this)
{
#ifdef TIMER_DEBUG
  printf(" Root OldTimer: %s  ID_counter=%d tableid=%d\n", Name().c_str(), ID_counter, tableid);
#endif
  setup();
  m_tableID = tableid;

  // Make a list of pointers of all Parent OldTimers
  if (TimerList==NULL) TimerList = new list<OldTimer*>();
  TimerList->push_back(this);
}

// Non-root Managed parent/child OldTimer construction
OldTimer::OldTimer(const string& name, OldTimer& parent, const int tableid)
  :
  m_name(name),
  m_Parent(parent)
{
#ifdef TIMER_DEBUG
  printf(" OldTimer: %s  ID_counter=%d tableid=%d\n", Name().c_str(), ID_counter, tableid);
#endif
  setup();
  m_tableID = tableid;

  if (TimerList==NULL) TimerList = new list<OldTimer*>();
  TimerList->push_back(this);
}

// child-only OldTimer
OldTimer::OldTimer(const string& name, OldTimer& parent)
  :
  m_name(name),
  m_Parent(parent)
{
#ifdef TIMER_DEBUG
  printf(" child OldTimer: %s  ID_counter=%d\n", Name().c_str(), ID_counter);
#endif
  setup();

  if (TimerList==NULL) TimerList = new list<OldTimer*>();
  TimerList->push_back(this);
}

// diagnostic
OldTimer::OldTimer(const string& name, OldTimer& parent, const int tableid,
             const bool diag)
  :
  m_name(name),
  m_Parent(parent)
{
#ifdef TIMER_DEBUG
  printf(" Diagnostic OldTimer: %s  ID_counter=%d tableid=%d\n",
         Name().c_str(), ID_counter, tableid);
#endif
  setup();
  m_diagnostic = true;
  m_tableID = tableid;

  if (TimerList==NULL) TimerList = new list<OldTimer*>();
  TimerList->push_back(this);
}

// // diagnostic
// //void OldTimer::OldTimer(const string& name, OldTimer& parent,
// //               const int tableid, const bool diag):
// //  m_name(name),  m_Parent(parent)
// //{
// void OldTimer::define(const string& name, OldTimer& parent,
//                    const int tableid, const bool diag)
// {
//   m_name = name;
//   m_Parent = parent;
//   //printf(" Diagnostic OldTimer: %s  ID_counter=%d tableid=%d\n",
//   //     Name().c_str(), ID_counter, tableid);
//   setup();
//   m_diagnostic = true;
//   m_tableID = tableid;

//   if (TimerList==NULL) TimerList = new list<OldTimer*>();
//   TimerList->push_back(this);
// }

// // Counter OldTimer
// OldTimer::OldTimer(const string& name, const int table):
//   timer_name(name),  Parent(parent)
// {

//   //printf(" OldTimer: %s  ID_counter=%d\n", Name().c_str(), ID_counter);
//   setup();
//   diagnostic_table = table;

//   if (TimerList==NULL) TimerList = new list<OldTimer*>();
//   TimerList->push_back(this);
// }

// Unmanaged OldTimer construction
OldTimer::OldTimer(): m_Parent(*this)
{
  setup();
}

void OldTimer::setup()
{
  //printf(" OldTimer setup: %s  ID: %d\n", Name().c_str(), ID_counter);

  // set this OldTimer to curent id
  m_ID = ID_counter;
  // increment static integer.
  ID_counter++;

  m_count = 0;
  m_evenCountAcrossRanks = false;
  m_tableID = -1;
  m_diagnostic = false;
  //timer_on = false;

  m_accumulated_WCtime = 0.0;
  m_last_WCtime_stamp  = 0.0;

#ifdef PAPI
  m_accumulated_counter0 = 0;
  m_accumulated_counter1 = 0;

  m_values[0] = 0;
  m_values[1] = 0;

#ifdef FOUR_COUNTERS
  m_accumulated_counter2 = 0;
  m_accumulated_counter3 = 0;

  m_values[2] = 0;
  m_values[3] = 0;
#endif
#endif
}

void OldTimer::TimerInit(const int rank)
{
#ifdef PAPI

  // this only needs to be called once. by all ranks.
  // just after MPI_init, before any start/stop calls.
  // and before any PAPI calls.  not sure where to put it.
  //pout() << "before papi_init" << endl;
  PAPI_init();
  //pout() << "after papi_init" << endl;

  // I've had problems with the PAPI_init() and PAPI_hwinfo() calls on seaborg.
  //  We only need the CPU_MHZ to compute mflops.  So a hack is to hardcode
  //  the value and comment the call.
  //CPU_MHZ = 375.0;  // for seaborg
  PAPI_hwinfo(rank, CPU_MHZ);

#ifdef FOUR_COUNTERS
  PAPI_initMultiEvent(PAPI_TOT_CYC, PAPI_FP_INS,
                      PAPI_FMA_INS, PAPI_TLB_TL);

  // PAPI_TLB_TL = Total translation lookaside buffer misses

  //PAPI_FP_INS + PAPI_FMA_INS - (PM_FPU_LD_ST_ISSUES - PM_FPU_LD)

  //PAPI_initMultiEvent(PAPI_LD_INS, PAPI_L1_LDM,
  //                  PAPI_L1_STM, PAPI_L1_ICM,
  //                  &papiEventID);
  //PAPI_initMultiEvent(PAPI_TOT_CYC,
  //                  PAPI_SR_INS, PAPI_L2_LDM,
  //                PAPI_MEM_SCY,
  //                &papiEventID);

#else
  if (TIMER_COUNTER == 0)
  {
    PAPI_initDualEvent(PAPI_TOT_CYC, PAPI_FP_INS);
    //PAPI_initDualEvent(PAPI_TOT_INS, PAPI_FP_INS);
  }
  else if (TIMER_COUNTER == 1)
  {
    //PAPI_initDualEvent(PAPI_L1_TCM,  PAPI_L2_TCM);
    PAPI_initDualEvent(PAPI_L1_DCM,  PAPI_L2_DCM);
  }
  else
  {
    PAPI_initDualEvent(PAPI_TOT_INS, PAPI_BR_INS);
  }
  //PAPI_initSingleEvent(PAPI_TOT_CYC);
#endif

  PAPI_start(papiEventID);
#endif
}

// Destructor for Managed and Unmanaged OldTimers.
OldTimer::~OldTimer()
{
  //PAPI_stop(papiEventID, m_values);
  //  pout() << " OldTimer::~OldTimer() " << Name() << "\n";
}

void OldTimer::start(void)
{
  ++m_count;
  ++NtotalStartStops;
  m_last_WCtime_stamp = getTimeStampWC();

#ifdef PAPI

#ifdef NDEBUG
  PAPI_read(papiEventID, m_values);
#else
  CH_assert(PAPI_read(papiEventID, m_values) == PAPI_OK);
#endif

  m_previous_counter0 = m_values[0];
  m_previous_counter1 = m_values[1];
#ifdef FOUR_COUNTERS
  m_previous_counter2 = m_values[2];
  m_previous_counter3 = m_values[3];
#endif

#endif

  //printf(" strt OldTimer: %20s  v0=%20.10e pc0=%20.10e v1=%20.10e pc1=%20.10e\n",
  //     Name().c_str(),
  //     (double)m_values[0], (double)previous_counter0,
  //     (double)m_values[1], (double)previous_counter1);
  //if (m_values[0] < previous_counter0)
  //{
  //printf(" in strt, m_values[0] < previous_counter0\n");
  //}
#ifdef MEMORY_USAGE
  sprintf(stuff, "%30s start %10d  mem=%-10.3f\n",
          Name().c_str(), (long)Count(), get_memory_usage_from_OS());
  pout() << stuff;
#endif
}

void OldTimer::stop(void)
{
  m_accumulated_WCtime += getTimeStampWC() - m_last_WCtime_stamp;

#ifdef PAPI

#ifdef NDEBUG
  PAPI_read(papiEventID, m_values);
#else
  CH_assert(PAPI_read(papiEventID, m_values) == PAPI_OK);
#endif

  m_accumulated_counter0 += m_values[0] - m_previous_counter0;
  m_accumulated_counter1 += m_values[1] - m_previous_counter1;
#ifdef FOUR_COUNTERS
  m_accumulated_counter2 += m_values[2] - m_previous_counter2;
  m_accumulated_counter3 += m_values[3] - m_previous_counter3;
#endif

#endif //PAPI

  //printf(" stop OldTimer: %20s  papiEventID=%2d v0=%20.10e pc0=%20.10e v1=%20.10e pc1=%20.10e\n",
  //     Name().c_str(), papiEventID,
  //     (double)m_values[0], (double)m_previous_counter0,
  //     (double)m_values[1], (double)m_previous_counter1);

  //if (m_values[0] < m_previous_counter0)
  //{
  //printf(" in stop, m_values[0] < m_previous_counter0\n");
  //}
  //if (m_values[1] < m_previous_counter1)
  //{
  //  printf(" in stop, m_values[1] < m_previous_counter1\n");
  //}

  //CH_assert(m_values[0] > 0);

  //#ifndef NDEBUG
  //if (m_values[0] < m_previous_counter0 || m_values[1] < m_previous_counter1)
  //{
  //printf(" stop OldTimer: %20s  papiEventID=%2d v0=%20.10e pc0=%20.10e v1=%20.10e pc1=%20.10e\n",
  //       Name().c_str(), papiEventID,
  //       (double)m_values[0], (double)m_previous_counter0,
  //       (double)m_values[1], (double)m_previous_counter1);
  //}
  //CH_assert(m_values[0] >= m_previous_counter0);
  //CH_assert(m_values[1] >= m_previous_counter1);
  //#endif

#ifdef MEMORY_USAGE
  sprintf(stuff, "%30s stop  %10d  mem=%-10.3f\n",
          Name().c_str(), (long)Count(), get_memory_usage_from_OS());
  pout() << stuff;
#endif
}

void OldTimer::stop(Real& wc1)
{
  wc1 = getTimeStampWC() - m_last_WCtime_stamp;
  m_accumulated_WCtime += wc1;

#ifdef PAPI

#ifdef NDEBUG
  PAPI_read(papiEventID, m_values);
#else
  CH_assert(PAPI_read(papiEventID, m_values) == PAPI_OK);
#endif

  m_accumulated_counter0 += m_values[0] - m_previous_counter0;
  m_accumulated_counter1 += m_values[1] - m_previous_counter1;
#ifdef FOUR_COUNTERS
  m_accumulated_counter2 += m_values[2] - m_previous_counter2;
  m_accumulated_counter3 += m_values[3] - m_previous_counter3;
#endif

#endif //PAPI

#ifdef MEMORY_USAGE
  sprintf(stuff, "%30s stop  %10d  mem=%-10.3f\n",
          Name().c_str(), (long)Count(), get_memory_usage_from_OS());
  pout() << stuff;
#endif
}

void OldTimer::clear(void)
{
  //printf("clear() OldTimer\n");

  m_accumulated_WCtime = 0.;
  m_count = 0;

#ifdef PAPI
  m_accumulated_counter0 = 0;
  m_accumulated_counter1 = 0;
  m_previous_counter0 = 0;
  m_previous_counter1 = 0;
#ifdef FOUR_COUNTERS
  m_accumulated_counter2 = 0;
  m_accumulated_counter3 = 0;
  m_previous_counter2 = 0;
  m_previous_counter3 = 0;
#endif

#endif
}

inline double OldTimer::getTimeStampWC()
{
#ifdef CH_MPI
  return( MPI_Wtime() );
#else
  //gettimeofday(&tv, ( struct timezone * ) NULL);
  gettimeofday(&tv, &tz);
  //printf( "%d seconds, %d microseconds\n", tv.tv_sec, tv.tv_usec);
  //pout() << tvbuf->tv_sec << "\n";
  //return(  ((double)((tvbuf->tv_sec)*1000000 + tvbuf->tv_usec))*1.0e-6 );

  // this doesn't always work cuz of potential integer overflow.
  //return( (tv.tv_sec * 1000000 + tv.tv_usec) * 1.0e-6 );

  // this is what MPICH uses for MPI_Wtime on machines
  // which have gettimeofday.
  return((double)tv.tv_sec + 0.000001 * (double)tv.tv_usec);
#endif
}

double OldTimer::mflops(void)
{
  double dc=0;
#ifdef PAPI
  double pc0 = papi_counter0();
  double pc1 = papi_counter1();
#ifdef FOUR_COUNTERS
  double pc2 = papi_counter2();
  dc = computeDerivedCounter(pc0, pc1, pc2);
#else
  dc = computeDerivedCounter(pc0, pc1, 0);
#endif
#endif
  return dc;
}

// added by petermc, 18 Jan 2006
void OldTimer::writeTotalPct(const string& a_extra)
{
  float pctTime = (100. * m_accumulated_WCtime) / m_Parent.m_accumulated_WCtime;
  char outputline[1024];
  sprintf(outputline, "%16.16s%12.6f sec (%7.3f%%) %s",
          Name().c_str(), m_accumulated_WCtime, pctTime, a_extra.c_str());
  pout() << outputline << "\n";
}

// This is the Summary for the Managed OldTimers.  If running in
// parallel, i first reduce all of the accumulated times for each
// timer and get the minimum, maximum, and average.  Then i go thru
// the list of timers and make a list of parent timers.  From there i
// can make for loops that step thru the tree and print out the
// results.
void OldTimer::TimerSummary(void)
{
  // all ranks go in, but only rank0 writes the file.
  TimerSummary_(0);

  // all ranks need to clean up
#ifdef PAPI

  // not sure how important this cleanup and shutdown is...
  //for (int i=0; i<ID_counter; i++)
  //{
  //printf(" cleaning up event %d\n", i);
  //PAPI_cleanup();
  //}

  PAPI_shutdown();
#endif

  // clean up
  delete TimerList;
}

void OldTimer::TimerSummaryWithTITAfiles(void)
{
    // all ranks go in, but only rank0 writes the file.
  TimerSummary_(1);

#ifdef PAPI
  PAPI_shutdown();
#endif

  // clean up
  delete TimerList;
}

void OldTimer::TimerSummary_(const int itita)
{
  int rank=0, number_procs=1;
#ifdef CH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &number_procs);
#endif

#ifdef TIMER_DEBUG
  printf(" rank%3d in OldTimer::TimerSummary_\n", rank);
#endif

  if (rank==0)
  {
#ifdef PAPI
    printf("defined: PAPI\n");
#endif
#ifdef CH_AIX
    printf("defined: CH_AIX\n");
#endif
#ifdef FOUR_COUNTERS
    printf("defined: FOUR_COUNTERS\n");
#endif
#ifdef CH_MPI
    printf("defined: MPI\n");
#endif
#ifdef MEMORY_USAGE
    printf("defined: MEMORY_USAGE\n");
#endif
#ifdef CATFISH
    printf("defined: CATFISH\n");
#endif
#ifdef NDEBUG
    printf("defined: NDEBUG\n");
#endif
  }

  // don't write these files out all the darn time...
  if (itita == 1)
  {
    FILE *TT;
    char ttfilename[1024];
    if (number_procs < 999)
    {
      sprintf(ttfilename, "tita%03d", rank);
    }
    else
    {
      sprintf(ttfilename, "tita%04d", rank);
    }
    TT = fopen(ttfilename, "w");
    if (TT == NULL)
    {
      printf("problem opening tita\n");
      exit(0);
    }

#ifdef PAPI
    for (list<OldTimer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; ++tli)
    {
      double pc0 = (double)(*tli)->papi_counter0();
      double pc1 = (double)(*tli)->papi_counter1();
#ifdef FOUR_COUNTERS
      double pc2 = (double)(*tli)->papi_counter2();
      double pc3 = (double)(*tli)->papi_counter3();
      double dc =  computeDerivedCounter(pc0, pc1, pc2);
      fprintf(TT, " %-16s P: %-16s d%-3d n=%14.0f wc=%13.3f pc0=%17.0f pc1=%17.0f pc2=%17.0f pc3=%17.0f dc=%13.3f\n",
              (*tli)->Name().c_str(), (*tli)->m_Parent.Name().c_str(),
              (*tli)->tableID(), (double)((*tli)->Count()), (*tli)->wc_time(),
              pc0, pc1, pc2, pc3, dc);
#else
      double dc =  computeDerivedCounter(pc0, pc1, 0);
      fprintf(TT, " %-16s P: %-16s d%-3d n=%14.0f wc=%13.3f pc0=%17.0f pc1=%17.0f dc=%13.3f\n",
              (*tli)->Name().c_str(), (*tli)->m_Parent.Name().c_str(),
              (*tli)->tableID(), (double)((*tli)->Count()), (*tli)->wc_time(),
              pc0, pc1, dc);
#endif // FOUR_COUNTERS
    }
#else
    for (list<OldTimer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; ++tli)
    {
      fprintf(TT, " %-16s P: %-16s d%-3d n=%14.0f  wc=%13.3f\n",
              (*tli)->Name().c_str(), (*tli)->m_Parent.Name().c_str(),
              (*tli)->tableID(), (double)((*tli)->Count()), (*tli)->wc_time());
    }
#endif // PAPI
    fclose(TT);
  }

  double largest_time = 0.0;

  // go thru all OldTimers in static list and obtain the avg/min/max of
  // all processors onto rank0 -- obviously trivial for serial runs.
  for (list<OldTimer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; ++tli)
  {

#ifdef TIMER_DEBUG
    //printf(" OldTimer: %-16s Parent: %-16s diagnostic = %d  n=%d  wc=%10.3f\n",
    //     (*tli)->Name().c_str(), (*tli)->m_Parent.Name().c_str(),
    //     (*tli)->diagnostic_table, (*tli)->Count(), (*tli)->wc_time());
#endif

    double wc = (*tli)->wc_time();
    long long int nc = (*tli)->Count();

    if (number_procs == 1)
    {
      (*tli)->m_avgWC = wc;
      (*tli)->m_minWC = wc;
      (*tli)->m_maxWC = wc;
      (*tli)->m_totalCount = nc;
      (*tli)->m_avgCount = (double)nc;
      (*tli)->m_evenCountAcrossRanks = true;
    }
    else
    {

#ifdef CH_MPI

      // after here, rank0 will be the only guy with the right answers
      MPI_Reduce(&wc, &(*tli)->m_minWC, 1, MPI_DOUBLE,
                 MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&wc, &(*tli)->m_maxWC, 1, MPI_DOUBLE,
                 MPI_MAX, 0, MPI_COMM_WORLD);

      double temp;
      MPI_Reduce(&wc, &temp, 1, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);

      if (rank==0) (*tli)->m_avgWC = temp/number_procs;

      MPI_Reduce(&nc, &(*tli)->m_totalCount, 1, MPI_LONG_LONG,
                 MPI_SUM, 0, MPI_COMM_WORLD);

      long long int ncmin, ncmax;
      MPI_Reduce(&nc, &ncmin, 1, MPI_LONG_LONG,
                 MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&nc, &ncmax, 1, MPI_LONG_LONG,
                 MPI_MAX, 0, MPI_COMM_WORLD);
      if (ncmin == ncmax)
      {
        (*tli)->m_evenCountAcrossRanks = true;
      }
      else
      {
        (*tli)->m_evenCountAcrossRanks = false;
      }

      if (rank==0)
      {
        (*tli)->m_avgCount = (double)(*tli)->m_totalCount/number_procs;
      }
#endif

    }

    if ((*tli)->m_avgWC > largest_time) largest_time = (*tli)->m_avgWC;
  }

#ifdef PAPI
  // go thru all OldTimers in static list and obtain the sum of
  // PAPI counters onto rank0 -- obviously trivial for serial runs.
  for (list<OldTimer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; ++tli)
  {

    //printf(" OldTimer: %-16s Parent: %-16s diagnostic = %d  pc0=%14.5e pc1=%14.5e\n",
    ///     (*tli)->Name().c_str(), (*tli)->m_Parent.Name().c_str(),
    //    (*tli)->diagnostic_table,
    //     (double)(*tli)->papi_counter0(),
    //     (double)(*tli)->papi_counter1());

    double pc0 = (double)(*tli)->papi_counter0();
    double pc1 = (double)(*tli)->papi_counter1();
#ifdef FOUR_COUNTERS
    double pc2 = (double)(*tli)->papi_counter2();
    double pc3 = (double)(*tli)->papi_counter3();
#endif

    if (number_procs == 1)
    {
      (*tli)->m_totalPapiCounter0 = pc0;
      (*tli)->m_totalPapiCounter1 = pc1;
#ifdef FOUR_COUNTERS
      (*tli)->m_totalPapiCounter2 = pc2;
      (*tli)->m_totalPapiCounter3 = pc3;
#endif

    }
    else
    {

#ifdef CH_MPI
      MPI_Reduce(&pc0, &(*tli)->m_totalPapiCounter0, 1, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&pc1, &(*tli)->m_totalPapiCounter1, 1, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef FOUR_COUNTERS
      MPI_Reduce(&pc2, &(*tli)->m_totalPapiCounter2, 1, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&pc3, &(*tli)->m_totalPapiCounter3, 1, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#endif
    }

    //if (rank==0)
    //{
    // printf(" OldTimer: %-16s after reduce pc0=%14.5e pc1=%14.5e  pc0=%14.5e pc1=%14.5e\n",
    //        (*tli)->Name().c_str(),
    //        (*tli)->m_totalPapiCounter0,
    //        (*tli)->m_totalPapiCounter1,
    //        (*tli)->total_papi_counter0(),
    //        (*tli)->total_papi_counter1());
    // }
  }
#endif  // PAPI

#ifdef CH_MPI
  double combinedTotalStartStops=0;
  MPI_Reduce(&NtotalStartStops, &combinedTotalStartStops, 1, MPI_DOUBLE,
             MPI_SUM, 0, MPI_COMM_WORLD);
  double avgStartStops = combinedTotalStartStops/number_procs;
#else
  double avgStartStops = NtotalStartStops;
#endif

  // go away if you aren't rank 0
  if (rank != 0)
  {
    return;
  }

  //double largest_time = 0.0;
  //for (list<OldTimer*>::iterator pti=ParentList.begin() ; pti != ParentList.end() ; pti++)
  //{
  //double pa = (*pti)->m_avgWC*TIME_UNIT_FACTOR;
  // if (pa > largest_time) largest_time = pa;
  //}

  OldTimer TimerLoop;
  OldTimer iLooper0;
  const int Nloop = 10000;
  TimerLoop.start();
  for (int i=0; i<Nloop; i++)
  {
    iLooper0.start();    iLooper0.stop();
  }
  TimerLoop.stop();
  double TimerCost = TimerLoop.wc_time()/(double)iLooper0.Count();
  // We only do the above loop to get an estimate of the start/stop cost
  // and don't want this extra addition of start/stop calls to be counted
  // in the report, so we subtract them off here:
  NtotalStartStops -= Nloop+1 ;

  if (largest_time <= 0.0) largest_time = 0.01;

  pout() << "Total OldTimer Start/Stop Calls on rank " << rank << " = " << (int)NtotalStartStops << "\n";
  pout() << "Estimated Single OldTimer Start/Stop-combo Cost = " << TimerCost << " sec." << "\n";
  pout() << "Avg OldTimer Start/Stop Calls over all procs= " << avgStartStops << "\n";
  pout() << "Avg Estimated Total OldTimer Start/Stop-combo Cost = " << NtotalStartStops*TimerCost
         << " sec. (" << (NtotalStartStops*TimerCost)/largest_time*100.0
         << "%) " << "\n";
  // go away if you aren't rank 0
  if (rank != 0)
  {
    return;
  }

  pout() << " rank " << rank << " writing timmy.txt " << "\n";

  FILE *OUT;

  if (TIMER_COUNTER == 0)
  {
    OUT = fopen("timmy.txt", "w");
  }
  else if (TIMER_COUNTER == 1)
  {
    OUT = fopen("timmy.txt1", "w");
  }
  else
  {
    OUT = fopen("timmy.txt2", "w");
  }
  if (OUT == NULL)
  {
    printf("problem opening output file in OldTimer\n");
    exit(0);
  }

  fprintf(OUT, "Number of Processors: %d\n", number_procs );
  fprintf(OUT, "\n");
  writeParentTables(OUT, TimerCost);
  writeDiagnosticTables(OUT, TimerCost);

  fprintf(OUT, "Single OldTimer Start/Stop Cost = %15.7e sec.\n", TimerCost);
  fprintf(OUT, "Total OldTimer Start/Stop Calls = %20.10e \n", NtotalStartStops);
  fprintf(OUT, "Avg OldTimer Start/Stop per Proc= %20.10e \n", avgStartStops);
  fprintf(OUT, "Avg Est OldTimer Cost per Proc  = %15.3f sec.  (%7.2f%%)\n",
          NtotalStartStops*TimerCost/number_procs,
          (NtotalStartStops*TimerCost/number_procs)/largest_time*100.0);

  fclose(OUT);

#ifdef COUNT_PEX1
  FILE *QOUT;
  char qfilename[1024];
  sprintf(qfilename, "q%03d.dat", rank);
  QOUT = fopen(qfilename, "w");
  if (QOUT == NULL)
  {
    printf("problem opening pex out\n");
    exit(0);
  }

  for (unsigned int i=0; i<NQ; i++)
  {
    fprintf(OUT,"%20.10e  s0=%20.10e s1=%20.10e s2=%20.10e\n",
            (double)i, (double)qsubbox_len0[i],
            (double)qsubbox_len1[i], (double)qsubbox_len2[i]);
  }
  fprintf(OUT,"iter= %20.10e N=%20.10e\n",
          (double)qiter, (double)qNN);

  for (unsigned int i=0; i<qiter; i++)
  {
    fprintf(QOUT, "%12d %4d %4d %4d %4d\n",
            i, qNN[i], qsubbox_len2[i], qsubbox_len1[i], qsubbox_len0[i]);
  }
  fclose(QOUT);
#endif
}

void OldTimer::writeParentTables(FILE *out, const double TimerCost)
{
  int number_procs=1;
#ifdef CH_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &number_procs);
#endif

  // Create a list of Parent OldTimers from the list of OldTimers
  list<OldTimer*> ParentList;

#ifdef TIMER_DEBUG
  pout() << " size of ParentList = " << ParentList.size() << "\n";
#endif

  double largest_time=0.0;

  // for each timer in entire OldTimer List
  for (list<OldTimer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; tli++)
  {

    if ((*tli)->m_avgWC > largest_time) largest_time = (*tli)->m_avgWC;

#ifdef TIMER_DEBUG
    printf(" OldTimer #%3d: %-16s Parent: %-16s tableid = %d\n",
         (*tli)->m_ID,
         (*tli)->Name().c_str(),
         (*tli)->m_Parent.Name().c_str(),
         (*tli)->tableID());
#endif

    // add the Parent of this OldTimer to the Parent List
    bool add = false;
    if ((*tli)->m_Parent.tableID() >= 0) add = true;
    if ((*tli)->m_diagnostic) add = false;

    // but first, make sure it isn't already in the Parent List
    // for each parent currently in Parent List
    for (list<OldTimer*>::iterator pti=ParentList.begin() ; pti != ParentList.end() ; pti++)
    {
      if (*pti == &((*tli)->m_Parent))
      {
        add = false;
        break;
      }
    }

    // add them in ascending order based on the table ID
    if (add)
    {
#ifdef TIMER_DEBUG
      pout() << "  considering Parent:" << (*tli)->m_Parent.Name()
           << " tid=" << (*tli)->m_Parent.tableID() << "\n";
#endif
      if (ParentList.size() == 0)
      {
        ParentList.push_back( &((*tli)->m_Parent) );
#ifdef TIMER_DEBUG
        pout() << " first one" << "\n";
#endif
      }
      else
      {
        bool inserted=false;
        for (list<OldTimer*>::iterator pti=ParentList.begin() ; pti != ParentList.end() ; pti++)
        {
          if ((*tli)->m_Parent.tableID() < (*pti)->tableID())
          {
            ParentList.insert(pti, &((*tli)->m_Parent) );
            inserted = true;
#ifdef TIMER_DEBUG
            pout() << " inserted" << "\n";
#endif
            break;
          }
        }
        if (!inserted)
        {
          ParentList.push_back( &((*tli)->m_Parent) );
#ifdef TIMER_DEBUG
          pout() << " push_back" << "\n";
#endif
        }
      }
#ifdef TIMER_DEBUG
      pout() << "  Adding Parent:" << (*tli)->m_Parent.Name()
           << " tableid = " << (*tli)->tableID()
           << " wc=" << (*tli)->m_Parent.m_avgWC << "\n";
#endif
    }

#ifdef TIMER_DEBUG
    pout() << " ParentList size = " << ParentList.size() << "\n";
#endif
  }

#ifdef TIMER_DEBUG
  for (list<OldTimer*>::iterator pti=ParentList.begin() ; pti != ParentList.end() ; pti++)
  {
    printf(" Parent #%3d: %-16s Parent: %-16s tableid = %d\n",
           (*pti)->m_ID,
           (*pti)->Name().c_str(),
           (*pti)->m_Parent.Name().c_str(),
           (*pti)->tableID());
  }
#endif

  // for every Parent OldTimer -- make a new table
  for (list<OldTimer*>::iterator pti=ParentList.begin() ; pti != ParentList.end() ; pti++)
  {

    double table_avg_sum     = 0.0;
    double table_min_sum     = 0.0;
    double table_max_sum     = 0.0;
    double table_count_sum   = 0.0;
    double table_pc0_sum     = 0.0;
    double table_pc1_sum     = 0.0;
    double table_pc2_sum     = 0.0;
    double table_pc3_sum     = 0.0;
    double table_percent;

    double wc_avg = (*pti)->m_avgWC*TIME_UNIT_FACTOR;
    double wc_min = (*pti)->m_minWC*TIME_UNIT_FACTOR;
    double wc_max = (*pti)->m_maxWC*TIME_UNIT_FACTOR;
    double count_avg = (*pti)->m_avgCount;
    bool even = (*pti)->m_evenCountAcrossRanks;
    double parent_avg = wc_avg;

    if (count_avg == 0 && wc_avg == 0.0 && wc_min == 0.0)
    {
      // prolly don't want this table printed cuz it will prolly be all zeros
      continue;
    }

    fprintf(out, "\n");

    writeTableHeader(out);

    // Parent for each table
#ifdef PAPI
#ifdef FOUR_COUNTERS
    writeLineOfData(out, (*pti)->Name().c_str(), count_avg, even, 100*wc_avg/largest_time,
                    wc_avg, wc_min, wc_max,
                    (*pti)->total_papi_counter0(),
                    (*pti)->total_papi_counter1(),
                    (*pti)->total_papi_counter2(),
                    (*pti)->total_papi_counter3(),
                    number_procs);
#else
    writeLineOfData(out, (*pti)->Name().c_str(), count_avg, even, 100*wc_avg/largest_time,
                    wc_avg, wc_min, wc_max,
                    (*pti)->total_papi_counter0(),
                    (*pti)->total_papi_counter1(),
                    number_procs);
#endif

#else
    writeLineOfData(out, (*pti)->Name().c_str(), count_avg, even, 100*wc_avg/largest_time,
                    wc_avg, wc_min, wc_max, number_procs);
#endif

    writeSingleLineSeparator(out);

    // for every OldTimer (looking for children)
    for (list<OldTimer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; tli++)
    {

      if ((*tli)->m_diagnostic) continue;

      // if this OldTimer's Parent is the current Parent and
      //  if this OldTimer's Parent isn't equal to itself (Everything)
      if ( &((*tli)->m_Parent) == *pti && *tli != &((*tli)->m_Parent) )
      {

        wc_avg = (*tli)->m_avgWC*TIME_UNIT_FACTOR;
        wc_min = (*tli)->m_minWC*TIME_UNIT_FACTOR;
        wc_max = (*tli)->m_maxWC*TIME_UNIT_FACTOR;
        count_avg = (*tli)->m_avgCount;
        even = (*tli)->m_evenCountAcrossRanks;

        if (parent_avg  > 0.)
        {
          table_percent  = wc_avg/parent_avg*100.0;
        }
        else
        {
          table_percent = 0.;
        }

        table_avg_sum  += wc_avg;
        table_min_sum  += wc_min;
        table_max_sum  += wc_max;
        table_count_sum += count_avg;

        // Children
#ifdef PAPI
        table_pc0_sum += (*tli)->total_papi_counter0();
        table_pc1_sum += (*tli)->total_papi_counter1();
#ifdef FOUR_COUNTERS
        table_pc2_sum += (*tli)->total_papi_counter2();
        table_pc3_sum += (*tli)->total_papi_counter3();
        writeLineOfData(out, (*tli)->Name().c_str(), count_avg, even,
                        table_percent, wc_avg, wc_min, wc_max,
                        (*tli)->total_papi_counter0(),
                        (*tli)->total_papi_counter1(),
                        (*tli)->total_papi_counter2(),
                        (*tli)->total_papi_counter3(),
                        number_procs);
#else
        writeLineOfData(out, (*tli)->Name().c_str(), count_avg, even,
                        table_percent, wc_avg, wc_min, wc_max,
                        (*tli)->total_papi_counter0(),
                        (*tli)->total_papi_counter1(),
                        number_procs);
#endif
#else
        writeLineOfData(out, (*tli)->Name().c_str(), count_avg, even,
                        table_percent, wc_avg, wc_min, wc_max, number_procs);
#endif
      }

    } // end of table loop for this parent

    writeTableTotals(out, (*pti)->tableID(),
                     table_count_sum, parent_avg,
                     table_avg_sum, table_min_sum,  table_max_sum,
                     table_pc0_sum, table_pc1_sum,
                     table_pc2_sum, table_pc3_sum, TimerCost);

    fprintf(out,"\n");
  } // end of parent loop
}

void OldTimer::writeDiagnosticTables(FILE *out, const double TimerCost)
{
  // Create a set of Diagnostic Tables -- ie all of the OldTimers
  // with unique values of member data "diagnostic_table".
  set<int> diagTableSet;

  for (list<OldTimer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; tli++)
  {
    //printf(" OldTimer #%3d: %-16s Parent: %-16s tableid = %d diag=%d\n",
    //     (*tli)->m_ID,
    //     (*tli)->Name().c_str(),
    //     (*tli)->m_Parent.Name().c_str(),
    //     (*tli)->tableID(), (*tli)->m_diagnostic);

    if ((*tli)->m_diagnostic)
    {
      diagTableSet.insert((*tli)->tableID());
    }
  }

  //pout() << " Number of Diag Tables = " << diagTableSet.size() << "\n";
  //for (set<int>::iterator si=diagTableSet.begin(); si != diagTableSet.end(); si++)
  //{
  // pout() << " si = " << *si << "\n";
  //}

  int number_procs=1;
#ifdef CH_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &number_procs);
#endif

  // For each unique value of diagnostic_table, print a table.
  for (set<int>::iterator si=diagTableSet.begin(); si != diagTableSet.end(); si++)
  {

    fprintf(out, "\n\n Diagnostic Table %d\n", *si);

    // crude way of getting the Parent for this diagnostic table.
    // part of the crudeness is that i'm not checking that all
    // OldTimers in this diagnostic table have the same parent.
    OldTimer *DiagnosticParent = NULL;
    for (list<OldTimer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; tli++)
    {
      if ((*tli)->tableID() != *si) continue;
      DiagnosticParent = &((*tli)->m_Parent);
    }

    writeTableHeader(out);

    double table_avg_sum   = 0.0;
    double table_min_sum   = 0.0;
    double table_max_sum   = 0.0;
    double table_count_sum = 0.0;
    double table_pc0_sum   = 0.0;
    double table_pc1_sum   = 0.0;
    double table_pc2_sum   = 0.0;
    double table_pc3_sum   = 0.0;
    double table_percent;

    double wc_avg = DiagnosticParent->m_avgWC*TIME_UNIT_FACTOR;
    double wc_min = DiagnosticParent->m_minWC*TIME_UNIT_FACTOR;
    double wc_max = DiagnosticParent->m_maxWC*TIME_UNIT_FACTOR;
    double count_avg = DiagnosticParent->m_avgCount;
    bool even = DiagnosticParent->m_evenCountAcrossRanks;
    double parent_avg = wc_avg;

    // Parent for each table
#ifdef PAPI
#ifdef FOUR_COUNTERS
    writeLineOfData(out, DiagnosticParent->Name().c_str(),
                    count_avg, even, -1,
                    wc_avg, wc_min, wc_max,
                    DiagnosticParent->total_papi_counter0(),
                    DiagnosticParent->total_papi_counter1(),
                    DiagnosticParent->total_papi_counter2(),
                    DiagnosticParent->total_papi_counter3(),
                    number_procs);
#else
    writeLineOfData(out, DiagnosticParent->Name().c_str(),
                    count_avg, even, -1,
                    wc_avg, wc_min, wc_max,
                    DiagnosticParent->total_papi_counter0(),
                    DiagnosticParent->total_papi_counter1(),
                    number_procs);
#endif
#else
    writeLineOfData(out, DiagnosticParent->Name().c_str(),
                    count_avg, even, -1,
                    wc_avg, wc_min, wc_max, number_procs);
#endif

    writeSingleLineSeparator(out);

    for (list<OldTimer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; tli++)
    {

      // skip all OldTimers except ones belonging in this table
      if ((*tli)->tableID() != *si) continue;

      wc_avg = (*tli)->m_avgWC*TIME_UNIT_FACTOR;
      wc_min = (*tli)->m_minWC*TIME_UNIT_FACTOR;
      wc_max = (*tli)->m_maxWC*TIME_UNIT_FACTOR;
      count_avg = (*tli)->m_avgCount;
      even = (*tli)->m_evenCountAcrossRanks;

      if (parent_avg  > 0.)
      {
        table_percent  = wc_avg/parent_avg*100.0;
      }
      else
      {
        table_percent = 0.;
      }

      table_avg_sum  += wc_avg;
      table_min_sum  += wc_min;
      table_max_sum  += wc_max;
      table_count_sum += count_avg;

      // Children
#ifdef PAPI
      table_pc0_sum += (*tli)->total_papi_counter0();
      table_pc1_sum += (*tli)->total_papi_counter1();

#ifdef FOUR_COUNTERS
      table_pc2_sum += (*tli)->total_papi_counter2();
      table_pc3_sum += (*tli)->total_papi_counter3();
      writeLineOfData(out, (*tli)->Name().c_str(), count_avg, even,
                      table_percent, wc_avg, wc_min, wc_max,
                      (*tli)->total_papi_counter0(),
                      (*tli)->total_papi_counter1(),
                      (*tli)->total_papi_counter2(),
                      (*tli)->total_papi_counter3(),
                      number_procs);
#else
      writeLineOfData(out, (*tli)->Name().c_str(), count_avg, even,
                      table_percent, wc_avg, wc_min, wc_max,
                      (*tli)->total_papi_counter0(),
                      (*tli)->total_papi_counter1(),
                      number_procs);
#endif

#else
      writeLineOfData(out, (*tli)->Name().c_str(), count_avg, even,
                      table_percent, wc_avg, wc_min, wc_max, number_procs);
#endif

    } // end of table loop for this parent

    writeTableTotals(out, -1*(*si),
                     table_count_sum, parent_avg,
                     table_avg_sum, table_min_sum,  table_max_sum,
                     table_pc0_sum, table_pc1_sum,
                     table_pc2_sum, table_pc3_sum, TimerCost);

    fprintf(out,"\n\n");
  }
}

static void writeSingleLineSeparator(FILE *out)
{
#ifdef PAPI
  fprintf(out, "-------------------------------------------------");
  fprintf(out, "-------------------------------------------------\n");
#else
  fprintf(out, "-------------------------------------------------");
  fprintf(out, "---------------------------\n");
#endif
}

static void writeDoubleLineSeparator(FILE *out)
{
#ifdef PAPI
  fprintf(out, "=================================================");
  fprintf(out, "=================================================\n");
#else
  fprintf(out, "=================================================");
  fprintf(out, "===========================\n");
#endif
}

static void writeTableHeader(FILE *out)
{
#ifdef PAPI
  //            12345678901234567890
  fprintf(out, "         Totals      ");
  fprintf(out, "  WC per.");
  fprintf(out, "      count");

  fprintf(out, "    avg  ");
  fprintf(out, "     min  ");
  fprintf(out, "     max  ");

#ifdef FOUR_COUNTERS
  fprintf(out, "    cycles     FP INS     FMA's     TLB miss  mflops");
  //  fprintf(out, "    ");
#else
  if (TIMER_COUNTER==0)
  {
    fprintf(out, "    cycles     FP INS     mflops");
  }
  else if (TIMER_COUNTER==1)
  {
    fprintf(out, "    L1 TCM     L2 TCM     L2/L1");
  }
  else
  {
    fprintf(out, "    TOT INS    BR INS     BR/TOT");
  }
#endif

#else
  //            12345678901234567890
  fprintf(out, "         Totals      ");
  fprintf(out, "  WC per.");
  fprintf(out, "      count");

  fprintf(out, "        avg  ");
  fprintf(out, "       min  ");
  fprintf(out, "       max  ");
#endif

  fprintf(out, "\n");

  writeDoubleLineSeparator(out);
}

static void writeTableTotals(FILE *out,
                             const int table_number,
                             const double count_sum,
                             const double parent_avg,
                             const double avg_sum,
                             const double min_sum,
                             const double max_sum,
                             const double pc0_sum,
                             const double pc1_sum,
                             const double pc2_sum,
                             const double pc3_sum,
                             const double TimerCost)
{
  //   const char tcs[16] = "Timer Cost";

  //   double totalTimerCost = (double)count_sum * TimerCost;
  //   double percentTimerCost;
  //   if (parent_avg  > 0.)
  //   {
  //     percentTimerCost = totalTimerCost/parent_avg * 100.0;
  //   } else {
  //     percentTimerCost = 0.;
  //   }

  // #ifdef PAPI
  //   fprintf(out, "%16s (%6.2f%%) %9ld %8.2f [%8.2f,%8.2f]  %10.3e %10.3e %8.2f\n",
  //           tcs, percentTimerCost, (long)count_sum, totalTimerCost, 0.0, 0.0, 0.0, 0.0, 0.0);
  // #else
  //   fprintf(out, "%16s (%6.2f%%) %9ld %10.2f [%10.2f, %10.2f]\n",
  //           tcs, percentTimerCost, (long)count_sum, totalTimerCost, 0.0, 0.0);
  // #endif

  writeDoubleLineSeparator(out);

  double percent;
  if (parent_avg  > 0.)
  {
    percent = avg_sum/parent_avg * 100.0;
  }
  else
  {
    percent = 0.;
  }

  char stuff2[1024];
  if (table_number >= 0)
  {
    sprintf(stuff2, "T%-4d      table tots: ", table_number);
  }
  else
  {
    sprintf(stuff2, "DT%-4d     table tots: ", -1*table_number);
  }
#ifdef PAPI

#ifdef FOUR_COUNTERS
  double derived_counter = computeDerivedCounter(pc0_sum, pc1_sum, pc2_sum);
  fprintf(out, "%s%6.2f%%  %9ld %8.2f [%8.2f,%8.2f]  %10.3e %10.3e %10.3e %10.3e %6.2f\n",
          stuff2, percent, (long)count_sum,
          avg_sum, min_sum, max_sum,
          pc0_sum, pc1_sum,  pc2_sum, pc3_sum, derived_counter);
#else
  double derived_counter = computeDerivedCounter(pc0_sum, pc1_sum, 0);
  fprintf(out, "%s%6.2f%%  %9ld %8.2f [%8.2f,%8.2f]  %10.3e %10.3e %6.2f\n",
          stuff2, percent, (long)count_sum,
          avg_sum, min_sum, max_sum,
          pc0_sum, pc1_sum, derived_counter);
#endif
#else
  fprintf(out, "%s%6.2f%%  %9ld %10.2f [%10.2f, %10.2f]\n",
          stuff2, percent, (long)count_sum,
          avg_sum, min_sum, max_sum);
#endif
}

#ifdef PAPI
static void writeLineOfData(FILE *out,
                            const char *name,
                            const double count,
                            const bool even,
                            const double percent,
                            const double wc_avg,
                            const double wc_min,
                            const double wc_max,
                            const double counter0,
                            const double counter1,
                            const int number_procs)
{
  // PAPI counters should be coming in as sum totals of all
  //  processors.
  double derived_counter = computeDerivedCounter(counter0, counter1, 0);

  char dot;  if (even) dot = ' '; else dot = '.';

  if (percent < 0)
  {
    fprintf(out, "%21.21s           %9.0f%1c%8.2f [%8.2f,%8.2f]  %10.3e %10.3e %6.2f\n",
            name, count, dot, wc_avg, wc_min, wc_max,
            counter0, counter1,
            derived_counter);
  }
  else
  {
    fprintf(out, "%21.21s  %6.2f%%  %9.0f%1c%8.2f [%8.2f,%8.2f]  %10.3e %10.3e %6.2f\n",
            name, percent, count, dot, wc_avg, wc_min, wc_max,
            counter0, counter1,
            derived_counter);
  }
}
#endif

#if defined(PAPI) && defined(FOUR_COUNTERS)
static void writeLineOfData(FILE *out,
                            const char *name,
                            const double count,
                            const bool even,
                            const double percent,
                            const double wc_avg,
                            const double wc_min,
                            const double wc_max,
                            const double counter0,
                            const double counter1,
                            const double counter2,
                            const double counter3,
                            const int number_procs)
{
  // PAPI counters should be coming in as sum totals of all
  //  processors.
  double derived_counter = computeDerivedCounter(counter0, counter1, counter2);

  char dot;  if (even) dot = ' '; else dot = '.';

  if (percent < 0)
  {
    fprintf(out, "%21.21s           %9.0f%1c%8.2f [%8.2f,%8.2f]  %10.3e %10.3e %10.3e %10.3e %6.2f\n",
            name, count, dot, wc_avg, wc_min, wc_max,
            counter0, counter1, counter2, counter3,
            derived_counter);
  }
  else
  {
    fprintf(out, "%21.21s  %6.2f%%  %9.0f%1c%8.2f [%8.2f,%8.2f]  %10.3e %10.3e %10.3e %10.3e %6.2f\n",
            name, percent, count, dot, wc_avg, wc_min, wc_max,
            counter0, counter1, counter2, counter3,
            derived_counter);
  }
}
#endif

// regular, no PAPI
static void writeLineOfData(FILE *out,
                            const char *name,
                            const double count,
                            const bool even,
                            const double percent,
                            const double wc_avg,
                            const double wc_min,
                            const double wc_max,
                            const int number_procs)
{
  char dot;  if (even) dot = ' '; else dot = '.';

  if (percent < 0)
  {
    fprintf(out, "%21.21s           %9.0f%1c%10.2f [%10.2f, %10.2f]\n",
            name, count, dot, wc_avg, wc_min, wc_max);
  }
  else
  {
    fprintf(out, "%21.21s  %6.2f%%  %9.0f%1c%10.2f [%10.2f, %10.2f]\n",
            name, percent, count, dot, wc_avg, wc_min, wc_max);
  }
}

#ifdef PAPI

double computeDerivedCounter(double pc0, double pc1, double pc2)
{
  const double VALID_PROC_TIME = 1.0e-4;

  double proc_time = 1.0;
  double derived_counter = 0.0;

  if (TIMER_COUNTER == 0)
  {
    if (CPU_MHZ <= 0.0) CPU_MHZ = 1.0;
    proc_time = pc0/(CPU_MHZ*1.0e6);
    if (proc_time > VALID_PROC_TIME)
    {
#ifdef FOUR_COUNTERS
      // FP_INS + FP_FMA
      derived_counter = (pc1+pc2)/(proc_time*1.0e6);
#else
      derived_counter = (pc1)/(proc_time*1.0e6);
#endif
    }
    else
    {
      derived_counter = 0.0;
    }
  } else if (TIMER_COUNTER == 1)
  {
    if (pc0 > 0)
    {
      derived_counter = pc1/pc0;
    }
    else
    {
      derived_counter = 0.0;
    }
  }
  else
  {
    if (pc0 > 0)
    {
      derived_counter = pc1/pc0;
    }
    else
    {
      derived_counter = 0.0;
    }
  }

  return derived_counter;
}

static int PAPI_initSingleEvent(const int event)
{
  int retval;
  retval = PAPI_query_event(event);
  CH_assert(retval == PAPI_OK);

  retval = PAPI_create_eventset(&papiEventID);
  CH_assert(retval == PAPI_OK);

#ifdef CH_AIX
  retval = PAPI_add_event(papiEventID, event);
#else
  retval = PAPI_add_event(&papiEventID, event);
#endif
  CH_assert(retval == PAPI_OK);
  return retval; // to get rid of warning mesg
}

static int PAPI_initDualEvent(const int event1, const int event2)
{
  int retval;

  retval = PAPI_query_event(event1);
  CH_assert(retval == PAPI_OK);

  retval = PAPI_query_event(event2);
  CH_assert(retval == PAPI_OK);

  retval = PAPI_create_eventset(&papiEventID);
  CH_assert(retval == PAPI_OK);

#ifdef CH_AIX
  retval = PAPI_add_event(papiEventID, event1);
  CH_assert(retval == PAPI_OK);

  retval = PAPI_add_event(papiEventID, event2);
  CH_assert(retval == PAPI_OK);
#else
  retval = PAPI_add_event(&papiEventID, event1);
  CH_assert(retval == PAPI_OK);

  retval = PAPI_add_event(&papiEventID, event2);
  CH_assert(retval == PAPI_OK);
#endif
  return retval;  // to get rid of warning mesg
}

static void PAPI_initMultiEvent(const int event1, const int event2,
                                const int event3, const int event4)
{
  int retval;

  retval = PAPI_query_event(event1);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_query_event(event2);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_query_event(event3);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_query_event(event4);
  CH_assert(retval == PAPI_OK);

  retval = PAPI_create_eventset(&papiEventID);
  CH_assert(retval == PAPI_OK);

#ifdef CH_AIX
  retval = PAPI_add_event(papiEventID, event1);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_add_event(papiEventID, event2);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_add_event(papiEventID, event3);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_add_event(papiEventID, event4);
  CH_assert(retval == PAPI_OK);
#else
  retval = PAPI_add_event(&papiEventID, event1);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_add_event(&papiEventID, event2);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_add_event(&papiEventID, event3);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_add_event(&papiEventID, event4);
  CH_assert(retval == PAPI_OK);
#endif

  //retval = PAPI_set_multiplex(tag);
  //CH_assert(retval == PAPI_OK);
}

static int PAPI_checkEvent(const int event)
{
  long long int value = -1;
  int retval;

  if ( (retval = PAPI_query_event(event)) != PAPI_OK )
  {
    printf(" problem with PAPI_query_event retval=%d\n", retval);
    return retval;
  }
  if ( (PAPI_create_eventset(&papiEventID) != PAPI_OK) )
  {
    printf(" problem with PAPI_create_eventset\n");
    return retval;
  }

#ifdef CH_AIX
  if ( (PAPI_add_event(papiEventID, event) != PAPI_OK) )
#else
  if ( (PAPI_add_event(&papiEventID, event) != PAPI_OK) )
#endif
    {
      printf(" problem with PAPI_add_event\n");
      return retval;
    }
  if ( (PAPI_start(papiEventID) != PAPI_OK) )
  {
    printf(" problem with PAPI_start\n");
    return retval;
  }
  if ( (PAPI_stop(papiEventID, &value) != PAPI_OK) )
  {
    printf(" problem with PAPI_stop\n");
    return retval;
  }
  printf(" event looks ok.  value=%f\n", (double)value);
  return retval; // to get rid of warning mesg
}

static int PAPI_init(void)
{
  int retval;
  // only works on seaborg (newer PAPI?)
  pout() << " PAPI_library_init PAPI_VER_CURRENT=" << PAPI_VER_CURRENT << endl;
  //pout() << " PAPI_is_initialized()=" << PAPI_is_initialized() << endl;
  retval = PAPI_library_init(PAPI_VER_CURRENT);
  pout() << " PAPI_library_init retval=" << retval << endl;

  if (retval != PAPI_VER_CURRENT && retval > 0)
  {
    MayDay::Error("PAPI library version mismatch!");
  }

  CH_assert(retval == PAPI_VER_CURRENT);
  //retval = PAPI_start_counters(&Events,1))

  return retval; // to get rid of warning mesg
}

static void PAPI_hwinfo(const int rank, float &mhz)
{
  const PAPI_hw_info_t *hwinfo = NULL;
  hwinfo = PAPI_get_hardware_info();
  CH_assert( hwinfo != NULL );

  if (rank==0)
  {
    printf("Vendor string and code   : %s (%d)\n",hwinfo->vendor_string,hwinfo->vendor);
    printf("Model string and code    : %s (%d)\n",hwinfo->model_string,hwinfo->model);
    printf("CPU revision             : %f\n",hwinfo->revision);
    printf("CPU Megahertz            : %f\n",hwinfo->mhz);
    printf("CPU's in an SMP node     : %d\n",hwinfo->ncpu);
    printf("Nodes in the system      : %d\n",hwinfo->nnodes);
    printf("Total CPU's in the system: %d\n",hwinfo->totalcpus);

    int num_counters=PAPI_num_counters();
    printf("Number PAPI counters     : %d\n", num_counters);
    CH_assert (num_counters >= PAPI_OK);
  }

  mhz = hwinfo->mhz;
}

static void PAPI_initMultiplex(void)
{
  int retval = PAPI_multiplex_init();
  CH_assert(retval == PAPI_OK);
}

static void PAPI_cleanup(void)
{
#ifdef NDEBUG  // no debug
#ifdef CH_AIX
  PAPI_cleanup_eventset(papiEventID);
#else
  PAPI_cleanup_eventset(&papiEventID);
#endif
#else
#ifdef CH_AIX
  CH_assert(PAPI_cleanup_eventset(papiEventID) == PAPI_OK);
#else
  CH_assert(PAPI_cleanup_eventset(&papiEventID) == PAPI_OK);
#endif
#endif
}

#endif // on PAPI

#endif //on CH_NTIMER

#ifndef CH_NTIMER
#include "BaseNamespaceFooter.H"
#endif
