#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEIVFABI_H_
#define _BASEIVFABI_H_
#include "MayDay.H"
#include "IntVectSet.H"
#include "VoFIterator.H"
#include "parstream.H"
#include "SPMD.H"
#include "DebugOut.H"
#include "EBDebugOut.H"
#include "SPMD.H"
#include "NamespaceHeader.H"

template <class T>
bool BaseIVFAB<T>::s_verbose = false;

template <class T>
void
BaseIVFAB<T>::setVerbose(bool a_verbose)
{
  s_verbose = a_verbose;
}
template <class T>
bool BaseIVFAB<T>::s_verboseDebug = false;

template <class T>
void
BaseIVFAB<T>::setVerboseDebug(bool a_verboseDebug)
{
  s_verboseDebug = a_verboseDebug;
}
/******************/
template <class T> inline
const EBGraph&
BaseIVFAB<T>::getEBGraph() const
{
  return m_ebgraph;
}
/*************/
template <class T> inline
const IntVectSet&
BaseIVFAB<T>::getIVS() const
{
  return m_ivs;
}
/******************/
template <class T> inline
BaseIVFAB<T>::BaseIVFAB()
{
  setDefaultValues();
}
/******************/
template <class T> inline
BaseIVFAB<T>::~BaseIVFAB()
{
  clear();
}
/******************/
template <class T> inline
BaseIVFAB<T>::BaseIVFAB(const IntVectSet& a_ivsin,
                        const EBGraph&    a_ebgraph,
                        const int&        a_nvarin)
{
  setDefaultValues();
  define(a_ivsin, a_ebgraph, a_nvarin);
}
/******************/
template <class T> inline
void
BaseIVFAB<T>::define(const IntVectSet& a_ivsin,
                     const EBGraph&    a_ebGraph,
                     const int&        a_nvarin)
{
  clear();
  m_isDefined = true;
  CH_assert(a_nvarin > 0);

  m_ivs = a_ivsin;
  m_nComp = a_nvarin;
  m_ebgraph = a_ebGraph;
  m_nVoFs = 0;

  if (!a_ivsin.isEmpty())
    {
      Box minbox = a_ivsin.minBox();
      m_fab.resize(minbox, 1);
      m_fab.setVal(NULL);

      //figure out how long vector has to be
      IVSIterator ivsit(m_ivs);
      for (ivsit.reset(); ivsit.ok(); ++ivsit)
        {
          m_nVoFs += m_ebgraph.numVoFs(ivsit());
        }
      //now allocate the vector set the fab to go into
      //the pool of data at the first component of the first
      //vof
      if (m_nVoFs > 0)
        {
          // Note: clear() was called above so this isn't a memory leak
          //m_data = new T[m_nVoFs*m_nComp];
          m_data.resize(m_nVoFs*m_nComp);
          T* currentLoc = &m_data[0];
          for (ivsit.reset(); ivsit.ok(); ++ivsit)
            {
              int numVoFs = m_ebgraph.numVoFs(ivsit());
              m_fab(ivsit(), 0) = currentLoc;
              currentLoc += numVoFs;
            }
        }
    }
}
/******************/
template <class T> inline
void
BaseIVFAB<T>::setVal(const T& a_value)
{
  for (int ivec = 0; ivec < m_nVoFs*m_nComp; ivec++)
    {
      m_data[ivec] = a_value;
    }
}

/******************/
template <class T> inline
void
BaseIVFAB<T>::setVal(int ivar, const T& a_value)
{
  CH_assert(isDefined());
  CH_assert(ivar >= 0);
  CH_assert(ivar < m_nComp);
  int ioffset = ivar*m_nVoFs;
  for (int ivec = 0; ivec < m_nVoFs; ivec ++)
    {
      m_data[ivec+ioffset] = a_value;
    }
}
/******************/
template <class T> inline
void
BaseIVFAB<T>::setVal(const T&   a_value,
                     const Box& a_box,
                     int        a_nstart,
                     int        a_numcomp)

{
  CH_assert(isDefined());

  if ( !m_ivs.isEmpty() )
    {
      IntVectSet ivsIntersect = m_ivs;
      ivsIntersect &= a_box;
      BaseIVFAB<T>& thisFAB = *this;
      for (VoFIterator vofit(ivsIntersect, m_ebgraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          for (int icomp = a_nstart; icomp < a_numcomp; icomp++)
            {
              thisFAB(vof, icomp) = a_value;
            }
        }
    }
}
/******************/
template <class T> inline
void
BaseIVFAB<T>::copy(const Box& a_fromBox,
                   const Interval& a_dstInterval,
                   const Box& a_toBox,
                   const BaseIVFAB<T>& a_src,
                   const Interval& a_srcInterval)
{
  CH_assert(isDefined());
  CH_assert(a_src.isDefined());
  CH_assert(a_srcInterval.size() == a_dstInterval.size());
  CH_assert(a_dstInterval.begin() >= 0);
  CH_assert(a_srcInterval.begin() >= 0);
  CH_assert(a_dstInterval.end()   < m_nComp);
  CH_assert(a_srcInterval.end()   < a_src.m_nComp);

  if ((!m_ivs.isEmpty()) && (!a_src.m_ivs.isEmpty()))
    {
      CH_assert( (a_fromBox == a_toBox) || m_ebgraph.getDomain().isPeriodic() );
      Box intBox = a_fromBox;
      IntVectSet ivsIntersect = m_ivs;
      ivsIntersect &= a_src.m_ivs; //how did this ever work without this line?
      ivsIntersect &= intBox;  //
      BaseIVFAB<T>& thisFAB = *this;
      int compSize = a_srcInterval.size();
      for (VoFIterator vofit(ivsIntersect, m_ebgraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          for (int icomp = 0; icomp < compSize; icomp++)
            {
              int isrccomp = a_srcInterval.begin() + icomp;
              int idstcomp = a_dstInterval.begin() + icomp;
              thisFAB(vof, idstcomp) = a_src(vof, isrccomp);
            }
        }
    }
}
/******************/
template <class T> inline
int BaseIVFAB<T>::size(const Box& a_region,
                       const Interval& a_comps) const
{
  CH_assert(isDefined());


  //create set of cells in fab that are also in the input region
  IntVectSet subIVS = m_ivs;
  subIVS &= a_region;

  VoFIterator vofit(subIVS, m_ebgraph);

  const Vector<VolIndex>& vofs = vofit.getVector();
  //account for vof list
  int vofsize  = linearListSize(vofs);

  //add for each data point
  int datasize = 0;
  for (int ivof = 0; ivof < vofs.size(); ivof++)
    {
      if (s_verboseDebug)
        {
          IntVect iv = vofs[ivof].gridIndex();
          pout() << "BaseIVFAB::size vof " << vofs[ivof] << " Is irregular? " << m_ebgraph.isIrregular(iv) << endl;
        }
      for (int icomp = a_comps.begin(); icomp <= a_comps.end(); icomp++)
        {
          const T& dataPt = (*this)(vofs[ivof], icomp);
          int pointsize = CH_XD::linearSize(dataPt);
          datasize += pointsize;
        }

    }

  int retval = vofsize + datasize;
  if (s_verboseDebug)
  {
    pout() << "BaseIVFAB::size " << retval << endl;
  }
  return retval;
}
/********************/
template <class T> inline
void BaseIVFAB<T>::linearOut(void* a_buf,
                             const Box& a_region,
                             const Interval& a_comps) const
{
  CH_assert(isDefined());

  if (s_verboseDebug)
    {
      pout() << "" << endl;
      pout() << "BaseIVFAB<T>::linearOut " << " for box " << a_region << endl;
    }

  //create set of cells in fab that are also in the input region
  IntVectSet subIVS = m_ivs;
  subIVS &= a_region;

  VoFIterator vofit(subIVS, m_ebgraph);
  const Vector<VolIndex>& vofs = vofit.getVector();
  //output the vofs.
  unsigned char* charbuffer = (unsigned char *) a_buf;
  linearListOut(charbuffer, vofs);
  charbuffer += linearListSize(vofs);

  //output the data
  for (int ivof = 0; ivof < vofs.size(); ivof++)
    {
      const VolIndex& vof = vofs[ivof];
      if (s_verboseDebug)
        {
          pout() << "vof " << vof << endl;
        }
      for (int icomp = a_comps.begin(); icomp <= a_comps.end(); icomp++)
        {
          const T& dataPt = (*this)(vof, icomp);
          CH_XD::linearOut(charbuffer, dataPt);
          //increment the buffer offset
          charbuffer += linearSize(dataPt);
        }
    }
}
/********************/
template <class T> inline
void BaseIVFAB<T>::linearIn(void* a_buf, const Box& a_region, const Interval& a_comps)
{
  CH_assert(isDefined());
  //input the vofs

  if (s_verboseDebug)
    {
      pout() << "" << endl;
      pout() << "BaseIVFAB<T>::linearIn " << " for box " << a_region << endl;
    }
  Vector<VolIndex> vofs;
  unsigned char* charbuffer = (unsigned char *) a_buf;
  linearListIn(vofs, charbuffer);
  charbuffer += linearListSize(vofs);

  //input the data
  for (int ivof = 0; ivof < vofs.size(); ivof++)
    {
      const VolIndex& vof = vofs[ivof];
      if (s_verboseDebug)
        {
          pout() << "vof " << vof << endl;
        }
      for (int icomp = a_comps.begin(); icomp <= a_comps.end(); icomp++)
        {
          T& dataPt = (*this)(vof, icomp);
          CH_XD::linearIn(dataPt, charbuffer) ;
          //increment the buffer offset
          charbuffer += linearSize(dataPt);
        }
    }
}
/********************/
template <class T> inline
T*
BaseIVFAB<T>::getIndex(const VolIndex& a_vof, const int& a_comp) const
{
  CH_assert(isDefined());
  CH_assert(m_ivs.contains(a_vof.gridIndex()));
  CH_assert((a_comp >= 0) && (a_comp < m_nComp));

  T* dataPtr = (T*)(m_fab(a_vof.gridIndex(), 0));
  dataPtr += a_vof.cellIndex();
  dataPtr += a_comp*m_nVoFs;
  return dataPtr;
}
/********************/
template <class T> inline
void
BaseIVFAB<T>::clear()
{
  m_nComp = 0;
  m_nVoFs = 0;
  m_ivs.makeEmpty();
  m_fab.clear();
  //if (m_data != NULL)
  //  {
  //    delete[] m_data;
  //    m_data = NULL;
  //  }
  m_isDefined = false;
}
/*************************/
template <class T> inline
bool
BaseIVFAB<T>::isDefined() const
{
  return (m_isDefined);
}
/*************************/
template <class T> inline
int
BaseIVFAB<T>::numVoFs() const
{
  return m_nVoFs;
}
/*************************/
template <class T> inline
T*
BaseIVFAB<T>::dataPtr(const int& a_comp)
{
  CH_assert(a_comp >= 0);
  CH_assert(a_comp <= m_nComp);
  static T* dummy = NULL;
  if (m_nVoFs == 0) return dummy;
  //  T* retval = m_data + a_comp*m_nVoFs;
  T* retval = &(m_data[a_comp*m_nVoFs]);
  return  retval;
}
/*************************/
template <class T> inline
const T*
BaseIVFAB<T>::dataPtr(const int& a_comp) const
{
  CH_assert(a_comp >= 0);
  CH_assert(a_comp <= m_nComp);
  static T* dummy = NULL;
  if (m_nVoFs == 0) return dummy;
  const T* retval = &(m_data[a_comp*m_nVoFs]);
  //  const T* retval = m_data + a_comp*m_nVoFs;
  return  retval;
}
/*************************/
template <class T> inline
int
BaseIVFAB<T>::nComp() const
{
  return m_nComp;
}
/*************************/
template <class T> inline
T&
BaseIVFAB<T>::operator() (const VolIndex& a_ndin,
                          const int&      a_comp)
{
  CH_assert(isDefined());
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < m_nComp);
  T* dataPtr =  getIndex(a_ndin, a_comp);
  return(*dataPtr);
}
/**************************/
template <class T> inline
const T&
BaseIVFAB<T>::operator() (const VolIndex& a_ndin,
                          const int& a_comp) const
{
  CH_assert(isDefined());
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < m_nComp);

  T* dataPtr =  getIndex(a_ndin, a_comp);
  return(*dataPtr);
}
/******************/
template <class T> inline
void
BaseIVFAB<T>::setDefaultValues()
{
  m_isDefined = false;
  m_nVoFs = 0;
  m_nComp = 0;
  //m_data = NULL;
  m_data.resize(0);
}

#include "NamespaceFooter.H"
#endif

