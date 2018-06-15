#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// RegionGathers.cpp
// bvs, 05/26/06

#define _REGIONGATHER_CPP_

#include "RegionGather.H"
#include "Misc.H"
#include "parstream.H"

#include "NamespaceHeader.H"


bool  RegionGather::Message::operator < (const RegionGather::Message& rhs) const
{
  if (procID == rhs.procID)
    {
      int a1, a2, b1, b2;

      if (src < dest)
      {
        a1=src; a2=dest;
      }
      else
      {
        a1=dest; a2=src;
      }

      if (rhs.src < rhs.dest)
      {
        b1=rhs.src; b2=rhs.dest;
      }
      else
      {
        b1=rhs.dest; b2=rhs.src;
      }

      if (a1==b1)
      {
        return a2 < b2;
      }
      else
      {
        return a1 < b1;
      }
    }
  // else
  return procID < rhs.procID;
}


RegionGather::RegionGather()
{
}

void RegionGather::define(const ProblemDomain& a_domain,
                          const DisjointBoxLayout& a_layout,
                          int a_radius)
{

  if (a_domain.isPeriodic())
    {
      MayDay::Abort("Sorry, RegionGather::define not implemented for periodic");
    }

  m_messages.define(a_layout);
  m_local.define(a_layout);

  DataIterator dit = a_layout.dataIterator();
  LayoutIterator lit = a_layout.layoutIterator();

  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& b = a_layout.get(dit());
      IntVect center = (b.smallEnd()+b.bigEnd());
      center /= 2;
      int w1 = center[0] - a_radius; //moving window optimization
      int w2 = center[0] + a_radius;
      Vector<RegionGather::Message>& messages =  m_messages[dit];
      Vector<RegionGather::Message>& local    =  m_local[dit];
      int proc = a_layout.procID(dit());
      for (lit.begin(); lit.ok(); ++lit)
        {
          const Box& box = a_layout.get(lit());
          if (box.bigEnd()[0] >= w1)
            {
              IntVect distance = (box.smallEnd()+box.bigEnd());
              distance /= 2;
              distance = center - distance;
              bool connected = true;
              for (int i=0; i<CH_SPACEDIM; ++i)
              {
                if (Abs(distance[i]) > a_radius) connected = false;
              }
              if (connected)
                {
                  RegionGather::Message arc;
                  arc.distance = distance;
                  arc.src = dit().intCode();
                  arc.dest= lit().intCode();
                  arc.srcIndex  = dit();
                  arc.destIndex = DataIndex(lit());
                  arc.procID = a_layout.procID(lit());
                  if (arc.procID != proc)
                    messages.push_back(arc);
                  else
                    local.push_back(arc);
                }

              if (box.smallEnd()[0] > w2) lit.end();
            }
        }
      messages.sort();
    }
}


void RegionGather::dump() const
{
  DataIterator dit = m_messages.dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
    {
      const Vector<RegionGather::Message>& messages =  m_messages[dit];
      pout()<<m_messages.box(dit())<<" messages :"<<messages.size()<<"\n";
      for (int i=0; i<messages.size(); i++)
        {
          pout()<<messages[i].procID<<" "<<messages[i].distance<<" | ";
        }
      pout()<<std::endl;

    }
}

template void regionGather<Real>(const LayoutData<Real>& a_local,
                                  const RegionGather& a_copier,
                                  LayoutData<Vector<GatherObject<Real> > >& a_gatherObjects);


#include "NamespaceFooter.H"
