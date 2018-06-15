#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LAYOUTDATAI_H_
#define _LAYOUTDATAI_H_

#include "MayDay.H"
#include "DataIterator.H"
#include "SPMD.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

// template < > void LayoutData<Real>::allocate();

// template < > void LayoutData<int>::allocate();

template<class T>
inline const T& LayoutData<T>::operator [] (const DataIndex& a_index) const
{
  CH_assert(m_boxLayout.check(a_index));

  // using a DataIndex for data on another processor -
  // if you are reading this error in a debugger, then you haven't used
  // addBox(Box, int) correctly, or you are using a LayoutIterator to
  // access a data holder.
  // CH_assert(m_boxLayout.procID(a_index) == procID());

  return *(m_vector[a_index.datInd()]);
}

template<class T>
inline const T& LayoutData<T>::operator [] (const DataIterator& a_index) const
{
  return (*this)[a_index()];
}

template<class T>
inline T& LayoutData<T>::operator [] (const DataIndex& a_index)
{
  CH_assert(m_boxLayout.check(a_index));

  // using a DataIndex for data on another processor -
  // if you are reading this error in a debugger, then you haven't used
  // addBox(Box, int) correctly, or you are using a LayoutIterator to
  // access a data holder.
  //CH_assert(m_boxLayout.procID(a_index) == procID());

  return *(m_vector[a_index.datInd()]);
}

template<class T>
inline T& LayoutData<T>::operator [] (const DataIterator& a_index)
{

  return (*this)[a_index()];
}

template<class T>
inline Box LayoutData<T>::box(const DataIndex& a_index) const
{
  return m_boxLayout.get(a_index);
}

template<class T>
inline Box LayoutData<T>::box(const DataIterator& a_index) const
{
  return m_boxLayout.get(a_index());
}

template<class T>
inline DataIterator LayoutData<T>::dataIterator() const
{
  return m_boxLayout.dataIterator();
}

template<class T>
inline TimedDataIterator LayoutData<T>::timedDataIterator() const
{
  return m_boxLayout.timedDataIterator();
}

template<class T>
inline LayoutData<T>::LayoutData()
{
  m_boxLayout.close();
  m_callDelete = true;
}

template<class T>
inline LayoutData<T>::LayoutData(const BoxLayout& a_dp)
  :
  m_boxLayout(a_dp)
{
  CH_assert(a_dp.isClosed());
  m_callDelete = true;
  allocate();
}

template<class T>
inline void LayoutData<T>::define(const BoxLayout& a_dp)
{
  CH_assert(a_dp.isClosed());
  m_boxLayout = a_dp;
  allocate();
}

template<class T>
inline LayoutData<T>::~LayoutData()
{
  CH_TIME("~LayoutData");
  if (m_callDelete == true)
  {
    DataIterator it(dataIterator());
    for (; it.ok(); ++it)
    {
      delete m_vector[it().datInd()];
    }
#if CH_USE_MEMORY_TRACKING
    m_vector.bytes -= sizeof(T)*it.size();
#endif
  }
}

template<class T>
inline void LayoutData<T>::allocate()
{
  m_callDelete = true;

  for (unsigned int i = 0; i < m_vector.size(); ++i)
  {
    delete m_vector[i];
    m_vector[i] = NULL;
  }
#if CH_USE_MEMORY_TRACKING
  this->m_vector.bytes -= sizeof(T)*m_vector.size();
#endif

  DataIterator it(dataIterator());
#if CH_USE_MEMORY_TRACKING
  m_vector.bytes += sizeof(T)*it.size();
#endif
  m_vector.resize(it.size(), NULL);

  for (; it.ok(); ++it)
  {
    unsigned int index = it().datInd();
    if (m_vector[index] == NULL)
    {
      m_vector[index] = new T;
      if (m_vector[index] == NULL)
      {
        MayDay::Error("OutOfMemory in LayoutData::allocate");
      }
    }
  }
}

#include "NamespaceFooter.H"
#endif
