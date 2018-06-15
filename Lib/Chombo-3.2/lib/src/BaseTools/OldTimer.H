#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _OLDTIMER_H_
#define _OLDTIMER_H_

#ifndef CH_NTIMER

#ifdef CH_AIX
#define FOUR_COUNTERS
#endif

#include <cstdio>
#include "REAL.H"
#include "MayDay.H"
#include "Vector.H"
#include "ClockTicks.H"
#ifdef CH_MPI
#include "mpi.h"
#endif

extern "C"
{
  // these are for gettimeofday() wall-clock timing
#include <unistd.h>
#include <sys/time.h>

  // obviously, if you don't have PAPI installed on your system
  // (which requires kernel modification) then the counters wont
  // work, so don't compile with -DPAPI
#ifdef PAPI
#include <papi.h>
#endif
}

#ifndef TIMER_COUNTER
#define TIMER_COUNTER 0
#endif

#include <list>
#include <string>
#include "List.H"

#include <string>
#include <iostream>
#include "BaseNamespaceHeader.H"

using namespace std;

 /** A simple class for keeping track of elapsed time.

   CH_XD::OldTimer provides the ability to measure the passage of wall-clock time
   within units of code, and relate multiple measurements in a hierarchical
   manner.

   A timer can be started, stopped, reset and printed out.  Multiple timers
   can be instantiated and related to each other hierarchically: a timer instance
   may have a parent and children.  The hierarchical relationship is relevant only
   when printing out the summary report.

   In parallel, the timers on each proc operate independently.  Only the summary
   operations cause communication between procs.
 */
class OldTimer
{
  public:

    // Constructors

    /// Construct an unnamed timer that has no relation to any other instance
    OldTimer();                                                        // unmanaged
    /// Construct a named timer and add it to a table.
    OldTimer(const string& a_name,         const int a_tableID);       // root parent-only
    OldTimer(const string& a_name, OldTimer&, const int a_tableID);       // parent/child
    OldTimer(const string& a_name, OldTimer&);                            // child only
    OldTimer(const string& a_name, OldTimer&, const int a_tableID, bool); // diagnostic

    ~OldTimer();

    void setup();
    void start();
    void stop();
    void stop(Real& wc1);
    void clear();

    inline double getTimeStampWC();
    double wc_time()
    {
      return m_accumulated_WCtime;
    }
    double mflops();
    void writeTotalPct(const string& a_extra = "");

#ifdef PAPI
    long long int papi_counter0()
    {
      return m_accumulated_counter0;
    }
    long long int papi_counter1()
    {
      return m_accumulated_counter1;
    }
    double total_papi_counter0()
    {
      return m_totalPapiCounter0;
    }
    double total_papi_counter1()
    {
      return m_totalPapiCounter1;
    }
#ifdef FOUR_COUNTERS
    long long int papi_counter2()
    {
      return m_accumulated_counter2;
    }
    long long int papi_counter3()
    {
      return m_accumulated_counter3;
    }
    double total_papi_counter2()
    {
      return m_totalPapiCounter2;
    }
    double total_papi_counter3()
    {
      return m_totalPapiCounter3;
    }
#endif
#endif

    string Name()
    {
      return m_name;
    }
    long long int Count()
    {
      return m_count;
    }
    int tableID()
    {
      return m_tableID;
    }

    static list<OldTimer*> *TimerList; // must be initialized somewhere

    static void TimerSummary(const int bert, const int ernie)
    {
      TimerSummary(); // backwards compatibilty...
    }

    static void TimerInit(const int rank);
    static void TimerSummary(void);
    static void TimerSummaryWithTITAfiles(void);

  private:

    // internally called...
    static void TimerSummary_(const int itita);
    static void writeParentTables(FILE *out, const double TimerCost);
    static void writeDiagnosticTables(FILE *out, const double TimerCost);

    //bool timer_on;    // State of timer, either on(true) or off(false)
    int m_tableID;  // the table ID -- where parent goes in the summary

    double m_accumulated_WCtime, m_last_WCtime_stamp;

    // wall-clock timer data
#ifndef CH_MPI
    struct timeval tv;   //  Values from call to gettimeofday
    struct timezone tz;
#endif

    int m_ID;
    string m_name;
    bool m_diagnostic;
    OldTimer& m_Parent;

    bool m_evenCountAcrossRanks;
    double m_avgWC, m_minWC, m_maxWC, m_avgCount;

#ifdef PAPI
    long long int m_accumulated_counter0;
    long long int m_accumulated_counter1;
    long long int m_previous_counter0;
    long long int m_previous_counter1;
    double m_totalPapiCounter0;
    double m_totalPapiCounter1;
#ifdef FOUR_COUNTERS
    long long int m_values[4];
    long long int m_accumulated_counter2;
    long long int m_accumulated_counter3;
    long long int m_previous_counter2;
    long long int m_previous_counter3;
    double m_totalPapiCounter2;
    double m_totalPapiCounter3;
#else
    long long int m_values[2];
#endif
#endif

    long long int m_count;
    long long int m_totalCount;
  };

/*
#else // CH_NTIMER

// stub OldTimer class
using namespace std;

class OldTimer
{
public:

  OldTimer()
  {
  }

  OldTimer(const string&, const int)
  {
  }

  OldTimer(const string&, OldTimer&, const int)
  {
  }

  OldTimer(const string&, OldTimer&)
  {
  }

  OldTimer(const string&, OldTimer&, const int, bool)
  {
  }

  ~OldTimer()
  {
  }

  void setup()
  {
  }
  void start()
  {
  }
  void stop()
  {
  }
  void stop(Real& wc1)
  {
  }
  void clear()
  {
  }

  inline double getTimeStampWC()
  {
    return 0.0;
  }
  double wc_time()
  {
    return 0.0;
  }
  double mflops();
  void writeTotalPct(const string& a_extra = "");

#ifdef PAPI
  long long int papi_counter0()
  {
    return 0;
  }
  long long int papi_counter1()
  {
    return 0;
  }
  double total_papi_counter0()
  {
    return 0.0;
  }
  double total_papi_counter1()
  {
    return 0.0;
  }
#ifdef FOUR_COUNTERS
  long long int papi_counter2()
  {
    return 0;
  }
  long long int papi_counter3()
  {
    return 0;
  }
  double total_papi_counter2()
  {
    return 0.0;
  }
  double total_papi_counter3()
  {
    return 0.0;
  }
#endif
#endif

  string Name()
  {
    return NULL;
  }
  long long int Count()
  {
    return 0;
  }
  int tableID()
  {
    return 0;
  }

  static void TimerSummary(const int bert, const int ernie)
  {
    TimerSummary(); // backwards compatibilty...
  }

  static void TimerInit(const int rank)
  {
  }
  static void TimerSummary(void)
  {
  }
  static void TimerSummaryWithTITAfiles(void)
  {
  };

private:

  // internally called...
  static void TimerSummary_(const int itita)
  {
  }
  static void writeParentTables(FILE *out, const double TimerCost)
  {
  }
  static void writeDiagnosticTables(FILE *out, const double TimerCost)
  {
  }
};

*/

#include "BaseNamespaceFooter.H"
#endif // CH_NTIMER

#endif // CH_OLDTIMER_H
