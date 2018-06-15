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

#ifndef _REGIONGATHER_H_
#define _REGIONGATHER_H_

#include "IntVect.H"
#include "REAL.H"
#include "LayoutData.H"
#include "DisjointBoxLayout.H"
#include "ProblemDomain.H"
#include <list>
#include "SPMD.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"


///  return object for the regionGather function
template <class T>
class GatherObject
{
public:
  IntVect m_offset;
  T       m_value;
};


///  Data coordination class for the regionGather function
/**
  Currently, this class does not handle periodic ProblemDomain.


  the radius is a distance cut-off for the gathering operation.  An object T will
  be gathered from a neighbouring Box in the DisjointBoxLayout only if:

   a = IntVect at center of box 1
   b = IntVect at center of box 2

   abs(a[i]-b[i]) <= radius   for all i
*/
class RegionGather
{
public:
  RegionGather();

  void define(const ProblemDomain& a_domain,
              const DisjointBoxLayout& a_layout,
              int radius);

  void dump() const;

  class Message
  {
  public:
    IntVect distance;
    int src, dest;
    DataIndex srcIndex;
    DataIndex destIndex;
    int procID;
    bool operator < (const RegionGather::Message& rhs) const;
  };


  LayoutData<Vector<RegionGather::Message> >  m_messages;
  LayoutData<Vector<RegionGather::Message> >  m_local;


};

///  Function performs a distance cut-off gather operation
/**
  the input a_local contains a value of T for every Box in the BoxLayout.

  the output a_gatherObjects contains all the T objects within the constraints
  established by the RegionGather object, and the center seperation of their
  respective Boxes.
*/
template <class T>
void regionGather(const LayoutData<T>& a_local,
                  const RegionGather& a_copier,
                  LayoutData<Vector<GatherObject<T> > >& a_gatherObjects);


#include "RegionGatherI.H"

#include "NamespaceFooter.H"
#endif
