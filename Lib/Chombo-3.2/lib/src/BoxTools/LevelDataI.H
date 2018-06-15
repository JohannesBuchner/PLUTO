#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELDATAI_H_
#define _LEVELDATAI_H_

#include <cstdlib>
#include <algorithm>
using std::sort;

#include "parstream.H"
#include "CH_Timer.H"
#include <float.h>
#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
template<class T>
LevelData<T>::LevelData()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
LevelData<T>::~LevelData()
{
  CH_TIME("~LevelData");
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
LevelData<T>::LevelData(const DisjointBoxLayout& dp, int comps, const IntVect& ghost,
                        const DataFactory<T>& a_factory)
  : m_disjointBoxLayout(dp), m_ghost(ghost)
{
#ifdef CH_MPI
  this->numSends = 0;
  this->numReceives = 0;
#endif
  this->m_boxLayout = dp;
  this->m_comps = comps;
  this->m_isdefined = true;

  if (!dp.isClosed())
    {
      MayDay::Error("non-disjoint DisjointBoxLayout: LevelData<T>::LevelData(const DisjointBoxLayout& dp, int comps)");
    }

  Interval interval(0, comps-1);
  this->allocateGhostVector(a_factory, ghost);
  this->setVector(*this, interval, interval); // Does nothing.
}
//-----------------------------------------------------------------------

// Since I need to thwart the user from invoking the
// 'define' methods that use a general BoxLayout, I cannot
// get at said functions myself now. Ha!  So, I have to recode
// them here.

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::define(const DisjointBoxLayout& dp, int comps, const IntVect& ghost,
                          const DataFactory<T> & a_factory)
{
  CH_TIME("LevelData<T>::define(dbl,comps,ghost,factory)");
  // clear exchange copier if it's already been defined
  if (this->m_isdefined)
    {
      m_exchangeCopier.clear();
    }

  this->m_isdefined = true;
  if (!dp.isClosed())
    {
      MayDay::Error("non-disjoint DisjointBoxLayout: LevelData<T>::define(const DisjointBoxLayout& dp,....)");
    }
  if (comps<=0)
    {
      MayDay::Error("LevelData::LevelData(const BoxLayout& dp, int comps)  comps<=0");
    }
  this->m_comps = comps;
  this->m_boxLayout = dp;

  m_disjointBoxLayout = dp;
  m_ghost = ghost;

  // Interval interval(0, comps-1);
  this->allocateGhostVector(a_factory, ghost);
  //  this->setVector(*this, interval, interval);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::define(const LevelData<T>& da,  const DataFactory<T> & a_factory)
{
  CH_TIME("LevelData<T>::define(LevelData<T>,factory)");
  // clear exchange copier if it's already been defined
  if (this->m_isdefined)
    {
      m_exchangeCopier.clear();
    }
  this->m_isdefined = true;
  if (this == &da) return;
  m_disjointBoxLayout = da.m_disjointBoxLayout;
  this->m_boxLayout  = da.m_disjointBoxLayout;
  this->m_comps     = da.m_comps;
  m_ghost     = da.m_ghost;

  this->allocateGhostVector(a_factory, m_ghost);
  da.copyTo(*this);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::define(const LevelData<T>& da, const Interval& comps,
                          const DataFactory<T>& a_factory)
{
  CH_TIME("LevelData<T>::define(LevelData<T>,comps,factory)");
  // clear exchange copier if it's already been defined
  if (this->m_isdefined)
    {
      m_exchangeCopier.clear();
    }
  this->m_isdefined = true;
  if (this == &da)
  {
    MayDay::Error(" LevelData<T>::define(const LevelData<T>& da, const Interval& comps) called with 'this'");
  }
  CH_assert(comps.size()>0);
  // this line doesn't make any sense!
  //CH_assert(comps.end()<=this->m_comps);
  CH_assert(comps.begin()>=0);

  m_disjointBoxLayout = da.m_disjointBoxLayout;
  this->m_boxLayout  = da.m_disjointBoxLayout;

  this->m_comps = comps.size();

  m_ghost = da.m_ghost;

  Interval dest(0, this->m_comps-1);

  this->allocateGhostVector(a_factory, m_ghost);

  this->setVector(da, comps, dest);

}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(const Interval& srcComps,
                          BoxLayoutData<T>& dest,
                          const Interval& destComps) const
{
  if ((BoxLayoutData<T>*)this == &dest) return;

  if (this->boxLayout() == dest.boxLayout())
    {
      // parallel direct copy here, no communication issues
      for (DataIterator it(this->dataIterator()); it.ok(); ++it)
        {
          dest[it()].copy(this->box(it()),
                          destComps,
                          this->box(it()),
                          this->operator[](it()),
                          srcComps);
        }
      return;
    }

  Copier copier(m_disjointBoxLayout, dest.boxLayout());
  copyTo(srcComps, dest, destComps, copier);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(BoxLayoutData<T>& dest) const
{
  CH_assert(this->nComp() == dest.nComp());
  this->copyTo(this->interval(), dest, dest.interval());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(const Interval& srcComps,
                          LevelData<T>& dest,
                          const Interval& destComps) const
{
  if (this == &dest)
  {
    MayDay::Error("src == dest in copyTo function. Perhaps you want exchange ?");
  }

  if (this->boxLayout() == dest.boxLayout()  && dest.ghostVect() == IntVect::Zero)
    {
      // parallel direct copy here, no communication issues
      for (DataIterator it(this->dataIterator()); it.ok(); ++it)
        {
          dest[it()].copy(this->box(it()),
                          destComps,
                          this->box(it()),
                          this->operator[](it()),
                          srcComps);
        }
      return;
    }

  Copier copier(m_disjointBoxLayout, dest.getBoxes(), dest.m_ghost);
  copyTo(srcComps, dest, destComps, copier);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(LevelData<T>& dest) const
{
  CH_assert(this->nComp() == dest.nComp());
  this->copyTo(this->interval(), dest, dest.interval());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(const Interval& srcComps,
                          BoxLayoutData<T>& dest,
                          const Interval& destComps,
                          const Copier& copier,
                          const LDOperator<T>& a_op) const
{
  CH_TIME("copyTo");
#ifdef CH_MPI
  {
//    CH_TIME("MPI_Barrier copyTo");
//      MPI_Barrier(Chombo_MPI::comm);
  }
#endif

  this->makeItSo(srcComps, *this, dest, destComps, copier, a_op);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(BoxLayoutData<T>& dest,
                          const Copier& copier,
                          const LDOperator<T>& a_op) const
{
  CH_assert(this->nComp() == dest.nComp());
  this->copyTo(this->interval(), dest, dest.interval(), copier, a_op);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(const Interval& srcComps,
                          LevelData<T>& dest,
                          const Interval& destComps,
                          const Copier& copier,
                          const LDOperator<T>& a_op) const
{
  CH_TIME("copyTo");
#ifdef CH_MPI
  {
//    CH_TIME("MPI_Barrier copyTo");
//      MPI_Barrier(Chombo_MPI::comm);
  }
#endif
  this->makeItSo(srcComps, *this, dest, destComps, copier, a_op);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::copyTo(LevelData<T>& dest,
                          const Copier& copier,
                          const LDOperator<T>& a_op) const
{
  CH_assert(this->nComp() == dest.nComp());
  this->copyTo(this->interval(), dest, dest.interval(), copier, a_op);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::exchange(const Interval& comps)
{
  CH_TIME("exchange+copier");
  // later on we can code this part as a direct algorithm
  // by copying and pasting the code from the Copier::define code
  // for now, just do the easy to debug approach.
    if (!m_exchangeCopier.isDefined())
    {
      m_exchangeCopier.define(m_disjointBoxLayout, m_disjointBoxLayout, m_ghost, true);
    }
  exchange(comps, m_exchangeCopier);

  // if there aren't any ghost cells, there isn't really anything
  // to do here (also, if m_ghost == Zero, m_exchangeCopier
  // wasn't defined!
  //if (m_ghost != IntVect::Zero)
  //this->makeItSo(comps, *this, *this, comps, m_exchangeCopier);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::exchange(void)
{
  exchange(this->interval());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::exchange(const Interval& comps,
                            const Copier& copier)
{
  CH_TIME("exchange");
#ifdef CH_MPI
  {
//    CH_TIME("MPI_Barrier exchange");
//      MPI_Barrier(Chombo_MPI::comm);
  }
#endif
  this->makeItSo(comps, *this, *this, comps, copier);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::exchange(const Copier& copier)
{
  exchange(this->interval(), copier);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::exchangeNoOverlap(const Copier& copier)
{
  CH_TIME("exchangeNoOverlap");
  this->makeItSoBegin(this->interval(), *this, *this, this->interval(), copier);
  this->makeItSoEnd(*this, this->interval());
  this->makeItSoLocalCopy(this->interval(), *this, *this, this->interval(), copier);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::exchangeBegin(const Copier& copier)
{
  CH_TIME("exchangeBegin");
  this->makeItSoBegin(this->interval(), *this, *this, this->interval(), copier);
  this->makeItSoLocalCopy(this->interval(), *this, *this, this->interval(), copier);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::exchangeEnd()
{
  CH_TIME("exchangeEnd");
  this->makeItSoEnd(*this, this->interval());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::define(const BoxLayout& dp, int comps,  const DataFactory<T>& a_factory)
{
  MayDay::Error("LevelData<T>::define called with BoxLayout input");
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::define(const BoxLayout& dp)
{
  MayDay::Error("LevelData<T>::define called with BoxLayout input");
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::define(const BoxLayoutData<T>& da, const DataFactory<T>& a_factory )
{
  MayDay::Error("LevelData<T>::define called with BoxLayout input");
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::define(const BoxLayoutData<T>& da, const Interval& comps,
                          const DataFactory<T>& a_factory)
{
  MayDay::Error("LevelData<T>::define called with BoxLayout input");
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::apply( void (*a_Function)(const Box& box, int comps, T& t) )
{
  for (DataIterator it(this->dataIterator()); it.ok(); ++it)
    {
      a_Function(m_disjointBoxLayout.get(it()), this->m_comps, *(this->m_vector[it().datInd()]));
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T>
void LevelData<T>::apply( const ApplyFunctor & f )
{
  for (DataIterator it(this->dataIterator()); it.ok(); ++it)
    {
      //     unsigned int index = this->m_boxLayout.lindex(it());
      f(m_disjointBoxLayout.get(it()), this->m_comps, *(this->m_vector[it().datInd()]));
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template<class T> void
LevelData<T>::degenerate( LevelData<T>& a_to,
                          const SliceSpec& a_sliceSpec ) const
{
    DisjointBoxLayout toDBL;
    m_disjointBoxLayout.degenerate( toDBL, a_sliceSpec );
    IntVect toGhost;
    for ( int i=0;i<CH_SPACEDIM;++i )
    {
        if ( i != a_sliceSpec.direction )
        {
            toGhost[i] = m_ghost[i];
        } else
        {
            toGhost[i] = 0;
        }
    }
    a_to.define( toDBL, this->nComp(), toGhost );
    copyTo( a_to );
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
#endif
