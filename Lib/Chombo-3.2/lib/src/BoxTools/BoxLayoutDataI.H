#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BOXLAYOUTDATAI_H_
#define _BOXLAYOUTDATAI_H_

#include <cstdlib>
#include <algorithm>
#include <limits.h>
#include <list>
#include "CH_OpenMP.H"
#include "parstream.H"
#include "memtrack.H"
#include "Misc.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"
#include "BaseFabMacros.H"

using std::sort;

template<class T>
int BoxLayoutData<T>::s_verbosity = 0;

template <class T>
T* DefaultDataFactory<T>::create(const Box& box,
                                 int ncomps,
                                 const DataIndex& a_datInd) const
{
  return new T(box, ncomps);
}

template<class T>
inline bool BoxLayoutData<T>::isDefined() const
{
  return m_isdefined;
}

template <class T>
inline void BoxLayoutData<T>::setVector(const BoxLayoutData<T>& da,
                                        const Interval& srcComps,
                                        const Interval& destComps)
{
  if(&da != this)
    {
      DataIterator it=this->dataIterator();
      int nbox=it.size();
      //#pragma omp parallel for
      for(int box=0; box<nbox; box++)
	{
	  this->m_vector[it[box].datInd()]->copy( this->box(it[box]), destComps,
                                                  this->box(it[box]), da[it[box]], srcComps);
	}
    }
}

template<class T>
inline void BoxLayoutData<T>::define(const BoxLayoutData<T>& da, const Interval& comps,
                                     const DataFactory<T>& factory)
{
  if (this == &da)
    {
      MayDay::Error("BoxLayoutData<T>::define(const LayoutData<T>& da,.....) called with 'this'");
    }
  CH_assert(comps.size()>0);
  CH_assert(comps.end()<=m_comps);
  //AD: why are the two different 
  //  CH_assert(comps.end()<=da.m_comps);
  CH_assert(comps.begin()>=0);
  this->m_boxLayout = da.boxLayout();
  
  this->m_comps = comps.size();
  
  Interval dest(0, m_comps-1);
  allocateGhostVector(factory);
  setVector(da, comps, dest);
}

template<class T>
inline void BoxLayoutData<T>::define(const BoxLayout& boxes, int comps,
                                     const DataFactory<T>& factory)
{
  CH_assert(boxes.isClosed());
  this->m_boxLayout = boxes;
  m_comps = comps;
  m_isdefined = true;
  allocateGhostVector(factory);

}

template<class T>
inline void BoxLayoutData<T>::define(const BoxLayout& boxes)
{
  MayDay::Error("BoxLayoutData<T>::define(const BoxLayout& boxes)...needs comps");
}

template <class T>
inline   BoxLayoutData<T>::BoxLayoutData():m_comps(0) ,m_buff(0)
{
  m_isdefined = false;
#ifdef CH_MPI
  this->numSends = 0;
  this->numReceives = 0;
#endif
}
template<class T>
inline BoxLayoutData<T>::BoxLayoutData(const BoxLayout& boxes, int comps,
                                       const DataFactory<T>& factory)
  :m_comps(comps),m_buff(NULL)
{
  CH_assert(boxes.isClosed());
  this->m_boxLayout = boxes;
  m_isdefined = true;
  allocateGhostVector(factory);
#ifdef CH_MPI

  this->numSends = 0;
  this->numReceives = 0;
#endif
}

template<class T>
BoxLayoutData<T>::~BoxLayoutData()
{
  CH_TIME("~BoxLayoutData");
  completePendingSends();
}

template<class T>
inline void BoxLayoutData<T>::define(const BoxLayoutData<T>& da,
                                     const DataFactory<T>& factory)
{
  if (this != &da)
    {
      m_isdefined = da.m_isdefined;
      this->m_boxLayout = da.boxLayout();
      m_comps    = da.nComp();
      Interval srcAnddest(0, m_comps-1);
      allocateGhostVector(factory);
      setVector(da, srcAnddest, srcAnddest);
    }

}
template<class T>
inline void BoxLayoutData<T>::clear()
{
  if (this->m_callDelete == true)
    {
      for (unsigned int i=0; i<this->m_vector.size(); ++i)
        {
          delete this->m_vector[i];
          this->m_vector[i] = NULL;
        }
#if CH_USE_MEMORY_TRACKING
      this->m_vector.bytes -= sizeof(T)*this->m_vector.size();
#endif
    }
  m_isdefined = false;
}

template<class T>
inline void BoxLayoutData<T>::allocateGhostVector(const DataFactory<T>& factory, const IntVect& ghost)
{
  if (this->m_callDelete == true)
    {
      for (unsigned int i=0; i<this->m_vector.size(); ++i)
        {
          delete this->m_vector[i];
          this->m_vector[i] = NULL;
        }
#if CH_USE_MEMORY_TRACKING
      this->m_vector.bytes -= sizeof(T)*this->m_vector.size();
#endif
    }

  this->m_callDelete = factory.callDelete();

  DataIterator it(this->dataIterator()); int nbox=it.size();
#if CH_USE_MEMORY_TRACKING
  if (this->m_callDelete == true)
    {
      this->m_vector.bytes += sizeof(T)*it.size();
    }
#endif
  this->m_vector.resize(it.size(), NULL);
  //#pragma omp parallel for
  for(int i=0; i<nbox; i++)
    {
      unsigned int index = it[i].datInd();
      Box abox = this->box(it[i]);
      abox.grow(ghost);
      this->m_vector[index] = factory.create(abox, m_comps, it[i]);
      if (this->m_vector[index] == NULL)
	{
	  MayDay::Error("OutOfMemory in BoxLayoutData::allocateGhostVector");
	}
    }
}

template<class T>
inline void BoxLayoutData<T>::apply(void (*a_func)(const Box& box, int comps, T& t))
{
  DataIterator it(this->dataIterator()); int nbox=it.size();
  //there might be race conditions on 
#pragma omp parallel for
  for(int i=0; i<nbox; i++)

    {
      a_func(this->box(it[i]), m_comps, *(this->m_vector[ it[i].datInd() ]));
    }
}

//======================================================================
template <class T>
AliasDataFactory<T>::AliasDataFactory(BoxLayoutData<T>* a_original, const Interval& interval)
{
  define(a_original, interval);
}

template <class T>
void AliasDataFactory<T>::define(BoxLayoutData<T>* a_original, const Interval& interval)
{
  m_origPointer = a_original;
  m_interval    = interval;
}

template <class T>
T* AliasDataFactory<T>::create(const Box& a_box, int ncomps, const DataIndex& a_dataInd) const
{
  //CH_assert(this->box(a_dataInd) == a_box);
  CH_assert(ncomps = m_interval.size());
  T* rtn = new T(m_interval, m_origPointer->operator[](a_dataInd));
  return rtn;
}

template<class T>
void BoxLayoutData<T>::makeItSo(const Interval&   a_srcComps,
                                const BoxLayoutData<T>& a_src,
                                BoxLayoutData<T>& a_dest,
                                const Interval&   a_destComps,
                                const Copier&     a_copier,
                                const LDOperator<T>& a_op) const
{
  makeItSoBegin(a_srcComps, a_src, a_dest, a_destComps, a_copier, a_op);
  makeItSoLocalCopy(a_srcComps, a_src, a_dest, a_destComps, a_copier, a_op);
  makeItSoEnd(a_dest, a_destComps, a_op);
}

template<class T>
void BoxLayoutData<T>::makeItSoBegin(const Interval&   a_srcComps,
                                     const BoxLayoutData<T>& a_src,
                                     BoxLayoutData<T>& a_dest,
                                     const Interval&   a_destComps,
                                     const Copier&     a_copier,
                                     const LDOperator<T>& a_op) const
{
  // The following five functions are nullOps in uniprocessor mode

#ifdef CH_MPI

  allocateBuffers(a_src,  a_srcComps,
                  a_dest, a_destComps,
                  a_copier,
                  a_op);  //monkey with buffers, set up 'fromMe' and 'toMe' queues

  writeSendDataFromMeIntoBuffers(a_src, a_srcComps, a_op);

  // If there is nothing to recv/send, don't go into these functions
  // and allocate memory that will not be freed later.  (ndk)
  // The #ifdef CH_MPI is for the m_buff->m_toMe and m_buff->m_fromMe
  {
    CH_TIME("post_messages");
    this->numReceives = m_buff->m_toMe.size();

    if (this->numReceives > 0)
      {
	postReceivesToMe(); // all non-blocking
      }
  

    this->numSends = m_buff->m_fromMe.size();
    if (this->numSends > 0)
      {
        postSendsFromMe();  // all non-blocking
      }
  }    
#endif 
}
////
template<class T>
bool 
BoxLayoutData<T>::doingJustCopy(const LDOperator<T>& a_op) const
{
  const CopyLDOperator<T>* checkptr = dynamic_cast<const CopyLDOperator<T> *>(&a_op);
  bool justCopy = (checkptr != NULL);
  return justCopy;
}
////
template<class T>
void BoxLayoutData<T>::makeItSoLocalCopy(const Interval&   a_srcComps,
                                         const BoxLayoutData<T>& a_src,
                                         BoxLayoutData<T>& a_dest,
                                         const Interval&   a_destComps,
                                         const Copier&     a_copier,
                                         const LDOperator<T>& a_op) const
{
  CH_TIME("local copying");

  if(doingJustCopy(a_op))
    {
      makeItSoLocalCopyThread(a_srcComps, a_src, a_dest, a_destComps, a_copier, a_op);
    }
  else
    {
      //only thread safe operation i am sure of is copy
      makeItSoLocalCopyNoOMP(a_srcComps, a_src, a_dest, a_destComps, a_copier, a_op);
    }
}

template<class T>
void BoxLayoutData<T>::makeItSoLocalCopyThread(const Interval&   a_srcComps,
                                               const BoxLayoutData<T>& a_src,
                                               BoxLayoutData<T>& a_dest,
                                               const Interval&   a_destComps,
                                               const Copier&     a_copier,
                                               const LDOperator<T>& a_op) const
{
  CopyIterator it(a_copier, CopyIterator::LOCAL);  int items=it.size();
  //#pragma omp parallel for
  for (int n=0; n<items; n++)
    {
      const MotionItem& item = it[n];
      //const T& src = a_src[item.fromIndex];
      //T& dest = a_dest[item.toIndex];
      //int srcSize = src .size(item.fromRegion, a_srcComps);
      //int destSize= dest.size(item.toRegion, a_destComps);
      //if (srcSize != destSize) MayDay::Error("(srcSize != destSize) in makeItSo");
      a_op.op(a_dest[item.toIndex], item.fromRegion,
              a_destComps,
              item.toRegion,
              a_src[item.fromIndex],
              a_srcComps);

    }
}
template<class T>
void BoxLayoutData<T>::makeItSoLocalCopyNoOMP(const Interval&   a_srcComps,
                                              const BoxLayoutData<T>& a_src,
                                              BoxLayoutData<T>& a_dest,
                                              const Interval&   a_destComps,
                                              const Copier&     a_copier,
                                              const LDOperator<T>& a_op) const
{

  CopyIterator it(a_copier, CopyIterator::LOCAL);  int items=it.size();
  for (int n=0; n<items; n++)
    {
      const MotionItem& item = it[n];
      //const T& src = a_src[item.fromIndex];
      //T& dest = a_dest[item.toIndex];
      //int srcSize = src .size(item.fromRegion, a_srcComps);
      //int destSize= dest.size(item.toRegion, a_destComps);
      //if (srcSize != destSize) MayDay::Error("(srcSize != destSize) in makeItSo");
      a_op.op(a_dest[item.toIndex], item.fromRegion,
              a_destComps,
              item.toRegion,
              a_src[item.fromIndex],
              a_srcComps);

    }
}

template<class T>
void BoxLayoutData<T>::makeItSoEnd(BoxLayoutData<T>& a_dest,
                                   const Interval&   a_destComps,
                                   const LDOperator<T>& a_op) const
{
  // Uncomment and Move this out of unpackReceivesToMe()  (ndk)
  completePendingSends(); // wait for sends from possible previous operation

  unpackReceivesToMe(a_dest, a_destComps, a_op); // nullOp in uniprocessor mode

}

#ifndef CH_MPI
// uniprocessor version of all these nullop functions.
template<class T>
void BoxLayoutData<T>::completePendingSends() const
{
}

template<class T>
void BoxLayoutData<T>::allocateBuffers(const BoxLayoutData<T>& a_src,
                                       const Interval& a_srcComps,
                                       const BoxLayoutData<T>& a_dest,
                                       const Interval& a_destComps,
                                       const Copier&   a_copier,
                                       const LDOperator<T>& a_op
                                       ) const
{
}

template<class T>
void BoxLayoutData<T>::writeSendDataFromMeIntoBuffers(const BoxLayoutData<T>& a_src,
                                                      const Interval&     a_srcComps,
                                                      const LDOperator<T>& a_op) const
{
}

template<class T>
void BoxLayoutData<T>::postSendsFromMe() const
{
}

template<class T>
void BoxLayoutData<T>::postReceivesToMe() const
{
}

template<class T>
void BoxLayoutData<T>::unpackReceivesToMe(BoxLayoutData<T>& a_dest,
                                          const Interval&   a_destComps,
                                          const LDOperator<T>& a_op) const
{
}

template<class T>
void BoxLayoutData<T>::unpackReceivesToMe_append(LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                                 const Interval&   a_destComps,
                                                 int ncomp,
                                                 const DataFactory<T>& factory,
                                                 const LDOperator<T>& a_op) const
{
}

#else

// MPI versions of the above codes.

template<class T>
void BoxLayoutData<T>::completePendingSends() const
{
  CH_TIME("completePendingSends");
  if (this->numSends > 0)
    {
      m_sendStatus.resize(this->numSends);
      int result = MPI_Waitall(this->numSends, &(m_sendRequests[0]), &(m_sendStatus[0]));
      if (result != MPI_SUCCESS)
        {
          //hell if I know what to do about failed messaging here
        }
    }
  this->numSends = 0;
}

template<class T>
void BoxLayoutData<T>::allocateBuffers(const BoxLayoutData<T>& a_src,
                                       const Interval& a_srcComps,
                                       const BoxLayoutData<T>& a_dest,
                                       const Interval& a_destComps,
                                       const Copier&   a_copier,
                                       const LDOperator<T>& a_op) const
{
  CH_TIME("MPI_allocateBuffers");
  m_buff = &(((Copier&)a_copier).m_buffers);
  if (m_buff->isDefined(a_srcComps.size()) && T::preAllocatable()<2) return;

  m_buff->m_ncomps = a_srcComps.size();

  m_buff->m_fromMe.resize(0);
  m_buff->m_toMe.resize(0);
  size_t sendBufferSize = 0;
  size_t recBufferSize  = 0;
  // two versions of code here.  one for preAllocatable T, one not.

  T dummy;
  for (CopyIterator it(a_copier, CopyIterator::FROM); it.ok(); ++it)
    {
      const MotionItem& item = it();
      CopierBuffer::bufEntry b;
      b.item = &item;
      b.size = a_op.size(a_src[item.fromIndex], item.fromRegion, a_srcComps);
      sendBufferSize+=b.size;
      b.procID = item.procID;
      m_buff->m_fromMe.push_back(b);
    }
  sort(m_buff->m_fromMe.begin(), m_buff->m_fromMe.end());
  for (CopyIterator it(a_copier, CopyIterator::TO); it.ok(); ++it)
    {
      const MotionItem& item = it();
      CopierBuffer::bufEntry b;
      b.item = &item;
      if (T::preAllocatable() == 0)
        {
          b.size = a_op.size(dummy, item.fromRegion, a_destComps);
          recBufferSize+=b.size;
        }
      else if (T::preAllocatable() == 1)
        {
          b.size = a_op.size(a_dest[item.toIndex], item.fromRegion, a_destComps);
          recBufferSize+=b.size;
        }
      b.procID = item.procID;
      m_buff->m_toMe.push_back(b);
    }
  sort(m_buff->m_toMe.begin(), m_buff->m_toMe.end());

  if (T::preAllocatable() == 2) // dynamic allocatable, need two pass
    {
      CH_TIME("MPI_ Phase 1 of 2 Phase: preAllocatable==2");
      if (s_verbosity > 0) pout()<<"preAllocatable==2\n";

      // in the non-preallocatable case, I need to message the
      // values for the m_buff->m_toMe[*].size
      Vector<unsigned long> fdata;
      Vector<unsigned long> tdata;
      int count = 1;
      int scount = 1;
      if (m_buff->m_toMe.size() > 0)
        {
          tdata.resize(m_buff->m_toMe.size(), ULONG_MAX);
          m_receiveRequests.resize(numProc()-1);
          m_receiveStatus.resize(numProc()-1);
          MPI_Request* Rptr = &(m_receiveRequests[0]);

          int lastProc = m_buff->m_toMe[0].procID;
          int messageSize = 1;
          unsigned long * dataPtr = &(tdata[0]);
          int i = 1;

          for (;i<m_buff->m_toMe.size(); ++i)
            {
              CopierBuffer::bufEntry& b = m_buff->m_toMe[i];
              if (b.procID == lastProc)
                messageSize++;
              else
                {

                  MPI_Irecv(dataPtr, messageSize, MPI_UNSIGNED_LONG, lastProc,
                            1, Chombo_MPI::comm, Rptr);
                  Rptr++;

                  lastProc = b.procID;
                  messageSize = 1;
                  dataPtr = &(tdata[i]);
                  count++;
                }
            }

          MPI_Irecv(dataPtr, messageSize, MPI_UNSIGNED_LONG, lastProc,
                    1, Chombo_MPI::comm, Rptr );
        }

      if (m_buff->m_fromMe.size() > 0)
        {
          fdata.resize(m_buff->m_fromMe.size());
          fdata[0]=m_buff->m_fromMe[0].size;
          m_sendRequests.resize(numProc()-1);
          m_sendStatus.resize(numProc()-1);
          MPI_Request* Rptr = &(m_sendRequests[0]);

          int lastProc = m_buff->m_fromMe[0].procID;
          int messageSize = 1;
          unsigned long * dataPtr = &(fdata[0]);
          int i = 1;
          for (;i<m_buff->m_fromMe.size(); ++i)
            {
              fdata[i]    = m_buff->m_fromMe[i].size;
              CopierBuffer::bufEntry& b = m_buff->m_fromMe[i];
              if (b.procID == lastProc)
                messageSize++;
              else
                {
                  MPI_Isend(dataPtr, messageSize, MPI_UNSIGNED_LONG, lastProc,
                            1, Chombo_MPI::comm, Rptr);

                  Rptr++;
                  lastProc = b.procID;
                  messageSize = 1;
                  dataPtr = &(fdata[i]);
                  scount++;
                }
            }

          MPI_Isend(dataPtr, messageSize, MPI_UNSIGNED_LONG, lastProc,
                    1, Chombo_MPI::comm, Rptr);
        }

      if (m_buff->m_toMe.size() > 0)
        {

          int result = MPI_Waitall(count, &(m_receiveRequests[0]), &(m_receiveStatus[0]));
          if (result != MPI_SUCCESS)
            {
              MayDay::Error("First pass of two-phase communication failed");
            }

          for (int i=0; i<m_buff->m_toMe.size(); ++i)
            {
              CH_assert(tdata[i] != ULONG_MAX);
              m_buff->m_toMe[i].size = tdata[i];
              recBufferSize+= tdata[i];
            }
        }

      if (m_buff->m_fromMe.size() > 0)
        {

          int result = MPI_Waitall(scount, &(m_sendRequests[0]), &(m_sendStatus[0]));
          if (result != MPI_SUCCESS)
            {
              MayDay::Error("First pass of two-phase communication failed");
            }

        }
    }

  // allocate send and receveive buffer space.

  if (sendBufferSize > m_buff->m_sendcapacity)
    {
      freeMT((m_buff->m_sendbuffer));
      if (s_verbosity > 0) pout()<<"malloc send buffer "<<sendBufferSize<<std::endl;
      (m_buff->m_sendbuffer) = mallocMT(sendBufferSize);
      if ((m_buff->m_sendbuffer) == NULL)
        {
          MayDay::Error("Out of memory in BoxLayoutData::allocatebuffers");
        }
      m_buff->m_sendcapacity = sendBufferSize;
    }

  if (recBufferSize > m_buff->m_reccapacity)
    {
      freeMT(m_buff->m_recbuffer);
      if (s_verbosity > 0) pout()<<"malloc receive buffer "<<recBufferSize<<std::endl;
      m_buff->m_recbuffer = mallocMT(recBufferSize);
      if (m_buff->m_recbuffer == NULL)
        {
          MayDay::Error("Out of memory in BoxLayoutData::allocatebuffers");
        }
      m_buff->m_reccapacity = recBufferSize;
    }

  /*
    pout()<<"\n";
    for (int i=0; i<m_buff->m_fromMe.size(); i++)
    pout()<<m_buff->m_fromMe[i].item->region<<"{"<<m_buff->m_fromMe[i].procID<<"}"<<" ";
    pout() <<"::::";
    for (int i=0; i<m_buff->m_toMe.size(); i++)
    pout()<<m_buff->m_toMe[i].item->region<<"{"<<m_buff->m_toMe[i].procID<<"}"<<" ";
    pout() << endl;
  */

  char* nextFree = (char*)(m_buff->m_sendbuffer);
  if (m_buff->m_fromMe.size() > 0)
    {
      for (unsigned int i=0; i<m_buff->m_fromMe.size(); ++i)
        {
          m_buff->m_fromMe[i].bufPtr = nextFree;
          nextFree += m_buff->m_fromMe[i].size;
        }
    }

  nextFree = (char*)m_buff->m_recbuffer;
  if (m_buff->m_toMe.size() > 0)
    {
      for (unsigned int i=0; i<m_buff->m_toMe.size(); ++i)
        {
          m_buff->m_toMe[i].bufPtr = nextFree;
          nextFree += m_buff->m_toMe[i].size;
        }
    }

  // since fromMe and toMe are sorted based on procID, messages can now be grouped
  // together on a per-processor basis.

}
///
template<class T>
void BoxLayoutData<T>::writeSendDataFromMeIntoBuffers(const BoxLayoutData<T>& a_src,
                                                      const Interval&     a_srcComps,
                                                      const LDOperator<T>& a_op) const
{
  CH_TIME("write Data to buffers");

  if(doingJustCopy(a_op))
    {
      writeSendDataFromMeIntoBuffersThread(a_src, a_srcComps, a_op);
    }
  else
    {
      writeSendDataFromMeIntoBuffersNoOMP (a_src, a_srcComps, a_op);
    }
}
///
template<class T>
void BoxLayoutData<T>::writeSendDataFromMeIntoBuffersThread(const BoxLayoutData<T>& a_src,
                                                            const Interval&     a_srcComps,
                                                            const LDOperator<T>& a_op) const
{
#pragma omp parallel for
  for (unsigned int i=0; i<m_buff->m_fromMe.size(); ++i)
    {
      const CopierBuffer::bufEntry& entry = m_buff->m_fromMe[i];
      a_op.linearOut(a_src[entry.item->fromIndex], entry.bufPtr,
                     entry.item->fromRegion, a_srcComps);
    }
}

///
template<class T>
void BoxLayoutData<T>::writeSendDataFromMeIntoBuffersNoOMP(const BoxLayoutData<T>& a_src,
                                                            const Interval&     a_srcComps,
                                                            const LDOperator<T>& a_op) const
{
  for (unsigned int i=0; i<m_buff->m_fromMe.size(); ++i)
    {
      const CopierBuffer::bufEntry& entry = m_buff->m_fromMe[i];
      a_op.linearOut(a_src[entry.item->fromIndex], entry.bufPtr,
                     entry.item->fromRegion, a_srcComps);
    }
}
///
template<class T>
void BoxLayoutData<T>::postSendsFromMe() const
{
  CH_TIME("post_Sends");
  // now we get the magic of message coalescence
  // fromMe has already been sorted in the allocateBuffers() step.

  this->numSends = m_buff->m_fromMe.size();

  if (this->numSends > 1)
    {
      for (unsigned int i=m_buff->m_fromMe.size()-1; i>0; --i)
        {
          if (m_buff->m_fromMe[i].procID == m_buff->m_fromMe[i-1].procID)
            {
              this->numSends--;
              m_buff->m_fromMe[i-1].size = m_buff->m_fromMe[i-1].size + m_buff->m_fromMe[i].size;
              m_buff->m_fromMe[i].size = 0;
            }
        }
    }
  m_sendRequests.resize(this->numSends);
  std::list<MPI_Request> extraRequests;

  unsigned int next=0;
  long long maxSize = 0;
  for (int i=0; i<this->numSends; ++i)
    {
      const CopierBuffer::bufEntry& entry = m_buff->m_fromMe[next];
      char*  buffer = (char*)entry.bufPtr;
      size_t bsize = entry.size;
      int idtag=0;
      while (bsize > CH_MAX_MPI_MESSAGE_SIZE)
        {
          extraRequests.push_back(MPI_Request());
          {
            //CH_TIME("MPI_Isend");
            MPI_Isend(buffer, CH_MAX_MPI_MESSAGE_SIZE, MPI_BYTE, entry.procID,
                      idtag, Chombo_MPI::comm, &(extraRequests.back()));
          }
          maxSize = CH_MAX_MPI_MESSAGE_SIZE;
          bsize -= CH_MAX_MPI_MESSAGE_SIZE;
          buffer+=CH_MAX_MPI_MESSAGE_SIZE;
          idtag++;
        }
      {
        //CH_TIME("MPI_Isend");
        MPI_Isend(buffer, bsize, MPI_BYTE, entry.procID,
                  idtag, Chombo_MPI::comm, &(m_sendRequests[i]));
      }
      maxSize = Max<long long>(bsize, maxSize);
      ++next;
      while (next < m_buff->m_fromMe.size() && m_buff->m_fromMe[next].size == 0) ++next;
    }
  for (std::list<MPI_Request>::iterator it = extraRequests.begin(); it != extraRequests.end(); ++it)
    {
      m_sendRequests.push_back(*it);
    }
  this->numSends = m_sendRequests.size();

  CH_MaxMPISendSize = Max<long long>(CH_MaxMPISendSize, maxSize);

}

template<class T>
void BoxLayoutData<T>::postReceivesToMe() const
{
  CH_TIME("post_Receives");
  this->numReceives = m_buff->m_toMe.size();

  if (this->numReceives > 1)
    {
      for (unsigned int i=m_buff->m_toMe.size()-1; i>0; --i)
        {
          if (m_buff->m_toMe[i].procID == m_buff->m_toMe[i-1].procID)
            {
              this->numReceives--;
              m_buff->m_toMe[i-1].size += m_buff->m_toMe[i].size;
              m_buff->m_toMe[i].size = 0;
            }

        }
    }
  m_receiveRequests.resize(this->numReceives);
  std::list<MPI_Request> extraRequests;
  unsigned int next=0;
  long long maxSize = 0;
  for (int i=0; i<this->numReceives; ++i)
    {
      const CopierBuffer::bufEntry& entry = m_buff->m_toMe[next];
      char*  buffer = (char*)entry.bufPtr;
      size_t bsize = entry.size;
      int idtag=0;
      while (bsize > CH_MAX_MPI_MESSAGE_SIZE)
        {
          extraRequests.push_back(MPI_Request());
          {
            //CH_TIME("MPI_Irecv");
            MPI_Irecv(buffer, CH_MAX_MPI_MESSAGE_SIZE, MPI_BYTE, entry.procID,
                      idtag, Chombo_MPI::comm, &(extraRequests.back()));
          }
          maxSize = CH_MAX_MPI_MESSAGE_SIZE;
          bsize -= CH_MAX_MPI_MESSAGE_SIZE;
          buffer+=CH_MAX_MPI_MESSAGE_SIZE;
          idtag++;
        }
      {
        //CH_TIME("MPI_Irecv");
        MPI_Irecv(buffer, bsize, MPI_BYTE, entry.procID,
                  idtag, Chombo_MPI::comm, &(m_receiveRequests[i]));
      }
      ++next;
      maxSize = Max<long long>(bsize, maxSize);
      while (next < m_buff->m_toMe.size() && m_buff->m_toMe[next].size == 0) ++next;
    }
  for (std::list<MPI_Request>::iterator it = extraRequests.begin(); it != extraRequests.end(); ++it)
    {
      m_receiveRequests.push_back(*it);
    }
  this->numReceives = m_receiveRequests.size();

  CH_MaxMPIRecvSize = Max<long long>(CH_MaxMPIRecvSize, maxSize);
  //pout()<<"maxSize="<<maxSize<<" posted "<<this->numReceives<<" receives\n";

}
///
template<class T>
void BoxLayoutData<T>::unpackReceivesToMe(BoxLayoutData<T>   & a_dest,
                                          const Interval     & a_destComps,
                                          const LDOperator<T>& a_op) const
{

  if(doingJustCopy(a_op))
    {
      unpackReceivesToMeThread(a_dest, a_destComps, a_op);
    }
  else
    {
      unpackReceivesToMeNoOMP( a_dest, a_destComps, a_op);
    }
}
///
template<class T>
void BoxLayoutData<T>::unpackReceivesToMeThread(BoxLayoutData<T>   & a_dest,
                                                const Interval     & a_destComps,
                                                const LDOperator<T>& a_op) const
{
  CH_TIME("unpack_messages");
  if (this->numReceives > 0)
    {
      m_receiveStatus.resize(this->numReceives);
      int result;
      {
        CH_TIME("MPI_Waitall");
        result = MPI_Waitall(this->numReceives, &(m_receiveRequests[0]),
                             &(m_receiveStatus[0]));
      }
      if (result != MPI_SUCCESS)
        {
          //hell if I know what to do about failed messaging here
        }
#pragma omp parallel for
      for (unsigned int i=0; i<m_buff->m_toMe.size(); ++i)
        {
          const CopierBuffer::bufEntry& entry = m_buff->m_toMe[i];
          a_op.linearIn(a_dest[entry.item->toIndex], entry.bufPtr, entry.item->toRegion, a_destComps);
        }
    }
  this->numReceives = 0;
}
///
template<class T>
void BoxLayoutData<T>::unpackReceivesToMeNoOMP (BoxLayoutData<T>   & a_dest,
                                                const Interval     & a_destComps,
                                                const LDOperator<T>& a_op) const
{
  CH_TIME("unpack_messages");
  if (this->numReceives > 0)
    {
      m_receiveStatus.resize(this->numReceives);
      int result;
      {
        CH_TIME("MPI_Waitall");
        result = MPI_Waitall(this->numReceives, &(m_receiveRequests[0]),
                             &(m_receiveStatus[0]));
      }
      if (result != MPI_SUCCESS)
        {
          //hell if I know what to do about failed messaging here
        }
      for (unsigned int i=0; i<m_buff->m_toMe.size(); ++i)
        {
          const CopierBuffer::bufEntry& entry = m_buff->m_toMe[i];
          a_op.linearIn(a_dest[entry.item->toIndex], entry.bufPtr, entry.item->toRegion, a_destComps);
        }
    }
  this->numReceives = 0;
}
///
template<class T>
void BoxLayoutData<T>::unpackReceivesToMe_append(LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                                 const Interval&   a_destComps,
                                                 int ncomp,
                                                 const DataFactory<T>& factory,
                                                 const LDOperator<T>& a_op) const
{
  if(doingJustCopy(a_op))
    {
      unpackReceivesToMe_appendThread(a_dest, a_destComps, ncomp, factory,  a_op);
    }
  else
    {
      unpackReceivesToMe_appendNoOMP( a_dest, a_destComps, ncomp, factory,  a_op);
    }
}
///
template<class T>
void BoxLayoutData<T>::unpackReceivesToMe_appendThread(LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                                 const Interval&   a_destComps,
                                                 int ncomp,
                                                 const DataFactory<T>& factory,
                                                 const LDOperator<T>& a_op) const
{
  if (this->numReceives > 0)
    {
      m_receiveStatus.resize(this->numReceives);
      int result;
      {
        CH_TIME("MPI_Waitall");
        result = MPI_Waitall(this->numReceives, &(m_receiveRequests[0]),
                             &(m_receiveStatus[0]));
      }
      if (result != MPI_SUCCESS)
        {
          //hell if I know what to do about failed messaging here
        }
#pragma omp parallel for
      for (unsigned int i=0; i<m_buff->m_toMe.size(); ++i)
        {
          const CopierBuffer::bufEntry& entry = m_buff->m_toMe[i];
          const MotionItem& item = *(entry.item);
          RefCountedPtr<T> newT( factory.create(item.toRegion, ncomp, item.toIndex) );;

          a_op.linearIn(*newT, entry.bufPtr, item.toRegion, a_destComps);
          a_dest[item.toIndex].push_back(newT);
        }
    }

  this->numReceives = 0;
}

///
template<class T>
void BoxLayoutData<T>::unpackReceivesToMe_appendNoOMP(LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                                      const Interval&   a_destComps,
                                                      int ncomp,
                                                      const DataFactory<T>& factory,
                                                      const LDOperator<T>& a_op) const
{
  if (this->numReceives > 0)
    {
      m_receiveStatus.resize(this->numReceives);
      int result;
      {
        CH_TIME("MPI_Waitall");
        result = MPI_Waitall(this->numReceives, &(m_receiveRequests[0]),
                             &(m_receiveStatus[0]));
      }
      if (result != MPI_SUCCESS)
        {
          //hell if I know what to do about failed messaging here
        }
      for (unsigned int i=0; i<m_buff->m_toMe.size(); ++i)
        {
          const CopierBuffer::bufEntry& entry = m_buff->m_toMe[i];
          const MotionItem& item = *(entry.item);
          RefCountedPtr<T> newT( factory.create(item.toRegion, ncomp, item.toIndex) );;

          a_op.linearIn(*newT, entry.bufPtr, item.toRegion, a_destComps);
          a_dest[item.toIndex].push_back(newT);
        }
    }

  this->numReceives = 0;
}
#endif

template <class T>
void BoxLayoutData<T>::generalCopyTo(const BoxLayout& a_destGrids,
                                     LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                     const Interval& a_srcComps,
                                     const ProblemDomain& a_domain,
                                     const Copier& a_copier,
                                     const DataFactory<T>& factory) const
{

  CH_assert(T::preAllocatable() == 0);
  a_dest.define(a_destGrids);

  CopyLDOperator<T> copyop;

  int ncomp = a_srcComps.size();
  Interval destComps(0, ncomp-1);
  allocateBuffers(*this,  a_srcComps,
                  *this,  destComps,
                  a_copier, copyop);

  writeSendDataFromMeIntoBuffers(*this, a_srcComps, copyop);

  // If there is nothing to recv/send, don't go into these functions
  // and allocate memory that will not be freed later.  (ndk)
  // The #ifdef CH_MPI is for the m_buff->m_toMe and m_buff->m_fromMe
#ifdef CH_MPI
  this->numReceives = m_buff->m_toMe.size();
  if (this->numReceives > 0)
    {
      postReceivesToMe(); // all non-blocking
    }

  this->numSends = m_buff->m_fromMe.size();
  if (this->numSends > 0)
    {
      postSendsFromMe();  // all non-blocking
    }
#endif

  // perform local copy
  CopyIterator it(a_copier, CopyIterator::LOCAL); int items=it.size();
  //#pragma omp parallel for
  for (int i=0; i<items; i++)
    {
      const MotionItem& item = it[i];
      RefCountedPtr<T> newT( factory.create(item.toRegion, ncomp, item.toIndex) );

      copyop.op(*newT, item.fromRegion,
                destComps,
                item.toRegion,
                this->operator[](item.fromIndex),
                a_srcComps);
      a_dest[item.toIndex].push_back(newT);
    }

  // Uncomment and Move this out of unpackReceivesToMe()  (ndk)
  completePendingSends(); // wait for sends from possible previous operation

  unpackReceivesToMe_append(a_dest, destComps, ncomp, factory, copyop); // nullOp in uniprocessor mode
}

template <class T>
void BoxLayoutData<T>::generalCopyTo(const BoxLayout& a_destGrids,
                                     LayoutData<Vector<RefCountedPtr<T> > >& a_dest,
                                     const Interval& a_srcComps,
                                     const ProblemDomain& a_domain,
                                     const DataFactory<T>& factory) const
{
  Copier copier;
  copier.define(this->m_boxLayout, a_destGrids, a_domain, IntVect::Zero);

  generalCopyTo(a_destGrids, a_dest, a_srcComps, a_domain, copier, factory);
}

template <class T>
void BoxLayoutData<T>::addTo(const Interval& a_srcComps,
                             BoxLayoutData<T>& a_dest,
                             const Interval& a_destComps,
                             const ProblemDomain& a_domain) const
{
  Copier copier;
  copier.define(this->m_boxLayout, a_dest.m_boxLayout, a_domain, IntVect::Zero);
  addTo(a_srcComps, a_dest, a_destComps, a_domain, copier);
}

template <class T>
class LDaddOp : public LDOperator<T>
{
public:
  virtual void op(T& dest,
                  const Box& RegionFrom,
                  const Interval& Cdest,
                  const Box& RegionTo,
                  const T& src,
                  const Interval& Csrc) const
  {
    dest.plus(src, RegionFrom, RegionTo, Csrc.begin(), Cdest.begin(), Cdest.size());
  }
  virtual void linearIn(T& arg,  void* buf, const Box& R,
                        const Interval& comps) const
  {
    Real* buffer = (Real*)buf;

    ForAllXBNNnoindx(Real, arg, R, comps.begin(), comps.size())
      {
        argR+=*buffer;
        buffer++;
      } EndFor

          }

};

template <class T>
void BoxLayoutData<T>::addTo(const Interval& a_srcComps,
                             BoxLayoutData<T>& a_dest,
                             const Interval& a_destComps,
                             const ProblemDomain& a_domain,
                             const Copier& a_copier) const
{
  CH_TIME("addTo");
  LDaddOp<T> addOp;
  makeItSo(a_srcComps, *this, a_dest, a_destComps, a_copier, addOp);
}

#include "NamespaceFooter.H"
#endif
