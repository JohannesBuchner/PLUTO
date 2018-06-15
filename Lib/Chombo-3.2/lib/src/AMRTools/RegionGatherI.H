#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// RegionGather
// bvs, 05/26/06

#ifndef _REGIONGATHERI_H_
#define _REGIONGATHERI_H_

#ifndef _REGIONGATHER_CPP_
template < > void regionGather<Real>(const LayoutData<Real>& a_local,
                                     const RegionGather& a_copier,
                                     LayoutData<Vector<GatherObject<Real> > >& a_gatherObjects);
#endif

class GatherBuffer
{
public:
  GatherBuffer(const RegionGather::Message* a_message):message(a_message)
  {
  }

  bool operator < (const GatherBuffer& rhs) const
  {
    return rhs.message->operator<(*(message));
  }

  const RegionGather::Message* message;
  char*                     buffer;

};


template <class T>
void regionGather(const LayoutData<T>& a_local,
                  const RegionGather& a_copier,
                  LayoutData<Vector<GatherObject<T> > >& a_gatherObjects)
{

  const BoxLayout& layout_input = a_local.boxLayout();
  const BoxLayout& layout_out   = a_gatherObjects.boxLayout();

  if (!(layout_input == layout_out))
    {
      a_gatherObjects.define(layout_input);
    }

  DataIterator dit=a_local.dataIterator();

  #ifdef CH_MPI
  //  JHVH help you if you need to extend or augment this bit of code.
  //  One asumption being made here is that the gathering operation is symmetric.
  //  if my boxes on my processor are sending data to your boxes on your processor, then
  //  you are sending an equivalent set of data back to me.
  //  another subtlety here is that the messages between processors are being agglomerated.
  //  the glist.sort operation trickles down to using the
  //  RegionGather::Message::operator< function,
  //  which primary sorts on procID.
  //  the *third* subtlety here, is that every possible box-pair
  //  has a canonical global ordering. This
  //  allows two processors to know what order
  //  they are sending/receiving their agglomerated data in.
  //          bvs.

  //  Phase 1, post all sends and symmetric receives
  Vector<T> sendBuffer;
  Vector<T> recvBuffer;
  std::list<GatherBuffer> glist;
  int count = 0;
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Vector<RegionGather::Message>& m = a_copier.m_messages[dit];
      for (int i=0; i<m.size(); ++i)
        {
          glist.push_back(&(m[i]));
          count++;
        }
    }

  glist.sort();

  sendBuffer.resize(count);
  recvBuffer.resize(count);

  Vector<MPI_Request> receives;
  Vector<MPI_Request> sends;

  std::list<GatherBuffer>::iterator iter;
  iter = glist.begin();
  int messageSize = 0;
  int plast;
  int index = 0;
  {
    CH_TIME("MPI_GatherSendRecv");
    for (int i=0; i<count; )
      {
        iter->buffer = (char*)&(recvBuffer[i]);
        sendBuffer[i] = a_local[iter->message->srcIndex];
        plast = iter->message->procID;
        ++messageSize;
        ++i;
        ++iter;
        if (i==count || plast != iter->message->procID)
          {
            receives.push_back(MPI_Request());
            sends.push_back(MPI_Request());
            int numChar = messageSize*sizeof(T); // units of sizeof is "char"
            MPI_Isend(&(sendBuffer[index]), numChar, MPI_CHAR, plast,
                      0, Chombo_MPI::comm, &(sends[sends.size()-1]));
            MPI_Irecv(&(recvBuffer[index]), numChar, MPI_CHAR, plast,
                      0, Chombo_MPI::comm, &(receives[receives.size()-1]));
            messageSize = 0;
            index = i;

          }
      }
  }



  #endif

  // Phase 2: perform all local on-processor gathers

  for (dit.begin(); dit.ok(); ++dit)
    {
      const T& val = a_local[dit];
      GatherObject<T> g;
      g.m_value = val; //T object requires operator=
      const Vector<RegionGather::Message>& local = a_copier.m_local[dit];
      for (int i=0; i<local.size(); ++i)
        {
          const RegionGather::Message& message = local[i];
          g.m_offset = message.distance;
          Vector<GatherObject<T> >& gather = a_gatherObjects[message.destIndex];
          gather.push_back(g);
        }
    }

#ifdef CH_MPI

  // phase 3: receive all messages and write all off-proceesor GatherObjects

  if (sends.size() > 0)
  {
    Vector<MPI_Status> status;
    status.resize(sends.size());
    int result;
    {
      CH_TIME("MPI_GatherSendWaitall");
      result = MPI_Waitall(sends.size(), &(sends[0]),
                           &(status[0]));
    }
  }
  if (receives.size() > 0)
  {
    Vector<MPI_Status> status;
    status.resize(receives.size());
    int result;
    {
      CH_TIME("MPI_GatherRecvWaitall");
      result = MPI_Waitall(receives.size(), &(receives[0]),
                           &(status[0]));
    }
    iter = glist.begin();

    for (int i=0; i<count; ++i, ++iter)
      {
        GatherObject<T> g;
        g.m_value = *((T*)iter->buffer);
        g.m_offset = -iter->message->distance;
        Vector<GatherObject<T> >& gather = a_gatherObjects[iter->message->srcIndex];
        gather.push_back(g);
      }
  }
#endif

}

#endif
