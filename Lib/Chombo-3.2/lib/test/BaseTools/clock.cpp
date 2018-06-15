#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ClockTicks.H"
#include <unistd.h>
#include <sys/time.h>
#include <iostream>

#include "parstream.H"
#ifdef CH_MPI
#include "mpi.h"
#endif

#include "UsingBaseNamespace.H"

/// Global variables for handling output:
static const char *pgmname = "clock" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;


inline void tod(unsigned long long& seconds, int& microseconds)
{
  struct timeval tv;   //  Values from call to gettimeofday
  struct timezone tz;
  gettimeofday(&tv, &tz);

  //pout()<<tv.tv_sec<<"  "<<tv.tv_usec<<"\n";
  seconds = tv.tv_sec;
  microseconds = tv.tv_usec;
}

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  int err =0;

  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << std::endl ;

#ifdef CH_TICKS
  int m0, m1;
  unsigned long long int t0, t1, t2, seconds0, seconds1;
  tod(seconds0, m0);
  t0=ch_ticks();
  t1=0;
  for (int i=0; i<32*16384; i++)
    {
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
      ch_ticks();
    }
  t2=ch_ticks();
  tod(seconds1, m1);
  if ( verbose )
    {
      pout()<<"t0 t2:"<<t0<<" "<<t2<<"\n";
      pout()<<"seconds "<<seconds0<<" "<<seconds1<<"\n";
      pout()<<"micro   "<<m0<<" "<<m1<<"\n";
    }
  t1=t2-t0;

  if ( verbose )
    {
      int micro = (seconds1-seconds0)*1000000 + (m1-m0);
      pout()<<"ticks: "<<t1<<" microseconds: "<<micro<<std::endl;

      double resolution = ((double)(micro*1000)/t1);
      pout()<<"tick resolution="<<resolution<<" nanoseconds"<<std::endl;
      t0=t1;
      t0/=32*32*16384;
      pout()<<"ticks per ch_ticks call: "<<t0<<std::endl;
    }

#else
  pout()<<"clock.cpp (and CH_TIME* macros) need a ch_ticks function"<<std::endl;
  err = 1;
#endif

  if ( verbose )
    {
      if (err == 0)
        {
          pout() << indent << "clock tick test passed." << std::endl;
        }
      else
        {
          pout() << indent << "clock tick test  FAILED!!! "<<err<<std::endl;
        }
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return err;
}

