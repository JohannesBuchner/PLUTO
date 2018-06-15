#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TIMEDDATAITERATOR_H_
#define _TIMEDDATAITERATOR_H_

#include "DataIndex.H"
#include "BoxLayout.H"
#include "SPMD.H"
#include "LayoutIterator.H"
#include "DataIterator.H"

#include "NamespaceHeader.H"


class TimedDataIterator : public DataIterator
{
public:

  virtual ~TimedDataIterator()
  {}

  ///
  virtual bool ok() const;

  ///
  virtual void operator++();

public:
  ///sets  m_time values to zero
  void clearTime();

  ///gets current time data
  Vector<unsigned long long>& getTime()
  {
    return m_time;
  }

  Vector<Box> getBoxes()
  {
    return m_layout.boxArray();
  }

  ///gets current time data in microseconds
  const Vector<unsigned long long>& getTime() const
  {
    return m_time;
  }

  ///enables timing.  does not set to zero.  use clear time for that
  void enableTime()
  {
    //only defines if not defined before
    defineTime();
    m_timeEnabled = true;
  }

  ///turns off timing
  void disableTime()
  {
    m_timeEnabled = false;
  }

  /// After you are finished timing your local elements, call mergeTimes to fill-in off-processor times.
  void mergeTime();
private:
  int m_size;
  void defineTime();
  Vector<unsigned long long> m_time;
  bool m_timeEnabled;
  bool m_timeDefined;
  unsigned long long m_startTime;

  int m_ibox;
private:
  TimedDataIterator();

  friend class BoxLayout;
  friend class DisjointBoxLayout;

  TimedDataIterator(const BoxLayout& boxlayout, const int* layoutID);
};


#include "NamespaceFooter.H"
#endif
