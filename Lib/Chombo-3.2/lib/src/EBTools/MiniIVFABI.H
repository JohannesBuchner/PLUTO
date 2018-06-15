#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MINIIVFABI_H_
#define _MINIIVFABI_H_
#include "MayDay.H"
#include "IntVectSet.H"
#include "VoFIterator.H"
#include "parstream.H"
#include "SPMD.H"
#include "DebugOut.H"
#include "EBDebugOut.H"
#include "SPMD.H"
#include "NamespaceHeader.H"


/******************/
template <class T> inline
MiniIVFAB<T>::MiniIVFAB()
{
  setDefaultValues();
}
/******************/
template <class T> inline
MiniIVFAB<T>::~MiniIVFAB()
{
  clear();
}
/******************/
template <class T> inline
MiniIVFAB<T>::MiniIVFAB(const IntVectSet& a_ivsin,
                        const EBGraph&    a_ebgraph,
                        const int&        a_nvarin)
{
  setDefaultValues();
  define(a_ivsin, a_ebgraph, a_nvarin);
}
/******************/
template <class T> inline
void
MiniIVFAB<T>::define(const IntVectSet& a_ivsin,
                     const EBGraph&    a_ebGraph,
                     const int&        a_nvarin)
{
  clear();
  this->m_isDefined = true;
  CH_assert(a_nvarin > 0);

  this->m_ivs = a_ivsin;
  this->m_nComp = a_nvarin;
  this->m_ebgraph = a_ebGraph;
  this->m_nVoFs = 0;

  if (!a_ivsin.isEmpty())
    {
      //figure out how long vector has to be
      IVSIterator ivsit(this->m_ivs);
      for (ivsit.reset(); ivsit.ok(); ++ivsit)
        {
          this->m_nVoFs += this->m_ebgraph.numVoFs(ivsit());
        }
      //now allocate the vector set the fab to go into
      //the pool of data at the first component of the first
      //vof
      if (this->m_nVoFs > 0)
        {
          //this->m_data = new T[this->m_nVoFs*this->m_nComp];
          this->m_data.resize(this->m_nVoFs*this->m_nComp);
          VoFIterator vofit(this->m_ivs, this->m_ebgraph);
          m_vofs = vofit.getVector();
        }
    }
}


/******************/
template <class T> inline
int MiniIVFAB<T>::size(const Box& a_region,
                       const Interval& a_comps) const
{
  CH_assert(this->isDefined());
  int count = 0;
  T tmp;
  //create set of cells in fab that are also in the input region
  for (int i=0; i<m_vofs.size(); i++)
    {
      if (a_region.contains(m_vofs[i].gridIndex())) count++;
    }
  if (count > 0)
    {
      return sizeof(int) + count*CH_XD::linearSize(m_vofs[0]) + count*a_comps.size()*CH_XD::linearSize(tmp);
    }
  return sizeof(int);
}
/********************/
template <class T> inline
void MiniIVFAB<T>::linearOut(void* a_buf,
                             const Box& a_region,
                             const Interval& a_comps) const
{
  CH_assert(this->isDefined());
  char* buffer = (char*)a_buf;
  buffer += sizeof(int);
  int count = 0;
  if (m_vofs.size()>0)
    {
      const T* ptr = &(this->m_data[0]);
      for (int i=0; i<m_vofs.size(); i++, ptr++)
        {
          const VolIndex& v = m_vofs[i];
          if (a_region.contains(v.gridIndex()))
            {
              count++;
              CH_XD::linearOut(buffer, v);
              buffer+= CH_XD::linearSize(v);
              for (int c=a_comps.begin(); c<=a_comps.end(); c++)
                {
                  CH_XD::linearOut(buffer, *(ptr+c*(this->m_nVoFs)));
                  buffer += CH_XD::linearSize(*ptr);
                }
            }
        }
    }
  int* b = (int*)a_buf;
  *b = count;
}
/********************/
template <class T> inline
void MiniIVFAB<T>::linearIn(void* a_buf, const Box& a_region, const Interval& a_comps)
{
  CH_assert(this->isDefined());
  int* b = (int*)a_buf;
  int count = *b;
  char* buffer = (char*)a_buf;
  buffer += sizeof(int);
  for (int i=0; i<count; i++)
    {
      VolIndex v;
      CH_XD::linearIn(v, buffer);
      buffer += linearSize(v);
      for (int c=a_comps.begin(); c<=a_comps.end(); c++)
        {
          T* ptr = getIndex(v, c);
          CH_XD::linearIn(*ptr, buffer);
          buffer+=CH_XD::linearSize(*ptr);
        }
    }
}
/********************/
template <class T> inline
T*
MiniIVFAB<T>::getIndex(const VolIndex& a_vof, const int& a_comp) const
{
  CH_assert(this->isDefined());
  CH_assert(this->m_ivs.contains(a_vof.gridIndex()));
  CH_assert((a_comp >= 0) && (a_comp < this->m_nComp));

  T* dataPtr =  (T*)&(this->m_data[0]);
  for (int i=0; i<m_vofs.size(); ++i)
    {
      if (a_vof == m_vofs[i]) break;
      dataPtr++;
    }
  dataPtr += a_comp*this->m_nVoFs;
  return dataPtr;
}
/********************/
template <class T> inline
void
MiniIVFAB<T>::clear()
{
  m_vofs.resize(0);
  BaseIVFAB<T>::clear();
}



/********************/
template <class T> inline
const Vector<VolIndex>&
MiniIVFAB<T>::getVoFs() const
{
  return m_vofs;
}


/******************/
template <class T> inline
void
MiniIVFAB<T>::setDefaultValues()
{
  BaseIVFAB<T>::setDefaultValues();
  m_vofs.resize(0);
}

#include "NamespaceFooter.H"
#endif

