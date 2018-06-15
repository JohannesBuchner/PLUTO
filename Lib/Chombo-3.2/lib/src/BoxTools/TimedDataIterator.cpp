#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TimedDataIterator.H"
#include "parstream.H"
#include "ClockTicks.H"
#ifdef CH_MPI
#include "SPMD.H"
#endif

#include "NamespaceHeader.H"
inline unsigned long long  getTimeTDI()
{
  return ch_ticks();
}

TimedDataIterator::TimedDataIterator(const BoxLayout& boxlayout, const int* layoutID)
  :DataIterator(boxlayout, layoutID)
{
  m_timeEnabled = false;
  m_timeDefined = false;
}
TimedDataIterator::TimedDataIterator()
{
  m_timeEnabled = false;
  m_timeDefined = false;
}
void TimedDataIterator::clearTime()
  {
    //only defines if not defined before
    defineTime();
    for (int ibox = 0; ibox < m_time.size(); ibox++)
      {
        m_time[ibox] = 0;
      }
  }

void TimedDataIterator::defineTime()
  {
    if (!m_timeDefined)
      {
        m_timeDefined = true;
        m_time.resize(m_layout.size(), 0);
      }
  }

void TimedDataIterator::operator++()
{
  if (m_timeEnabled)
    {

      unsigned long long tdiff = getTimeTDI();
      tdiff -= m_startTime;
      const DataIndex& current = this->operator()();
      m_time[current.m_index] += tdiff;

      //      pout() << "TimedDataIterator time[" << m_current.m_index
      //             << "]="  << m_time[m_current.m_index] <<endl;
    }
  DataIterator::operator++();
}
///
bool TimedDataIterator::ok() const
{
  TimedDataIterator & nonConstDI =(TimedDataIterator&)(*this);
  if (m_timeEnabled)
    {
      nonConstDI.m_startTime = getTimeTDI();
      //      pout() << "TimedDataIterator startTime="  << m_startTime << endl;
    }
  return DataIterator::ok();
}

void TimedDataIterator::mergeTime()
{
#ifdef CH_MPI
  int count = m_time.size();
  Vector<unsigned long long> tmp(count);
  MPI_Allreduce(&(m_time[0]),&(tmp[0]), count, MPI_LONG_LONG_INT, MPI_SUM, Chombo_MPI::comm);
  m_time = tmp;
#endif
}

#include "NamespaceFooter.H"
