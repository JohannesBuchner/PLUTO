#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#ifndef _BASEEBCELLFABI_H_
#define _BASEEBCELLFABI_H_

#include "BoxIterator.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"
/**********************/
/**********************/
template <class T>
void
BaseEBCellFAB<T>::setCoveredCellVal(const T&    a_val,
                                    const int&  a_comp,
                                    const bool& a_doMulti)
{
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < m_nComp);

  if (m_ebisBox.isAllRegular())
    {
      return;
    }
  else if (m_ebisBox.isAllCovered())
    {
      m_regFAB.setVal(a_val, a_comp);
    }
  else
    {
      for (BoxIterator bit(m_region); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          if (m_ebisBox.isCovered(iv))
            {
              m_regFAB(iv, a_comp) = a_val;
            }
        }
      if (a_doMulti)
        {
          //also set the multivalued cells
          for (IVSIterator ivsit(getMultiCells());ivsit.ok(); ++ivsit)
            {
              const IntVect& iv = ivsit();
              m_regFAB(iv, a_comp) = a_val;
            }
        }
    }
}
/*************/
/*************/
template <class T>
BaseEBCellFAB<T>&
BaseEBCellFAB<T>::copy(const BaseEBCellFAB<T>& a_src)
{
  CH_TIMERS("BaseEBCellFAB::copy");
  CH_TIMER("sameregion_copy", t1);
  CH_TIMER("overlap_copy", t2);
  CH_assert(isDefined());
  CH_assert(a_src.isDefined());
  // CH_assert(m_ebisBox == a_src.m_ebisBox);
  CH_assert(m_nComp <= a_src.m_nComp);
  CH_assert(m_region.sameType(a_src.m_region));

  Interval comps(0,m_nComp-1);

  if (m_region == a_src.m_region)
    {
      CH_START(t1);
      m_regFAB.copy(a_src.m_regFAB);

      if (m_hasMultiCells)
        {
          const T* src  = a_src.m_irrFAB.dataPtr(0);
          T* dest = m_irrFAB.dataPtr(0);
          T* end = dest + m_irrFAB.numVoFs()*m_nComp;
          for (; dest<end ; src++, dest++)
            {
              *dest = *src;
            }
        }
      CH_STOP(t1);
    }
  else
    {
      CH_START(t2);
      Box overlap(m_region);
      overlap &= a_src.m_region;

      this->copy(overlap,comps,overlap,a_src,comps);
      CH_STOP(t2);
    }

  return *this;
}

/*************/
/*************/
template <class T>
void
BaseEBCellFAB<T>::copy(const Box& a_RegionFrom,
                       const Interval& a_dstInt,
                       const Box& a_RegionTo,
                       const BaseEBCellFAB<T>& a_src,
                       const Interval& a_srcInt)
{
  CH_TIME("BaseEBCellFAB::copy");
  CH_assert(isDefined());
  CH_assert(a_dstInt.size()==a_srcInt.size());
  CH_assert(a_dstInt.begin()>=0);
  CH_assert(a_srcInt.begin()>=0);
  CH_assert(a_srcInt.end()< a_src.m_nComp);
  CH_assert(a_dstInt.end()< m_nComp);
  ProblemDomain domain = m_ebisBox.getDomain();
  bool isPeriodic = domain.isPeriodic();
  Box intersect = a_RegionFrom & domain.domainBox();
  if ( (!isPeriodic) &&
      (!intersect.isEmpty()) &&
      (!m_ebisBox.getRegion().contains(intersect)) )
    {
      MayDay::Error("BaseEBCellFAB::copy:  I was probably defined with too small an ebisBox");
    }

  CH_assert( (a_RegionFrom == a_RegionTo) || isPeriodic);
  {
    CH_TIME("BaseEBCellFAB::regcopy");
    m_regFAB.copy(a_RegionFrom, a_dstInt, a_RegionTo, a_src.m_regFAB, a_srcInt);
  }

  if (m_hasMultiCells)
    {
      CH_TIME("BaseEBCellFAB::irrcopy");
      m_irrFAB.copy(a_RegionFrom, a_dstInt, a_RegionTo, a_src.m_irrFAB, a_srcInt);
    }
}

/*************/
/*************/
template <class T>
int
BaseEBCellFAB<T>::size(const Box& R, const Interval& comps) const
{
  int retval =  m_regFAB.size(R, comps);
  retval += m_irrFAB.size(R, comps);
  return retval;
}

/*************/
template <class T>
void
BaseEBCellFAB<T>::linearOut(void* buf, const Box& R, const Interval& comps)
  const
{
  unsigned char* buffer = (unsigned char*)buf;
  m_regFAB.linearOut(buffer, R, comps);
  buffer+= m_regFAB.size(R, comps);
  m_irrFAB.linearOut(buffer, R, comps);
}
/*************/
/*************/
template <class T>
void
BaseEBCellFAB<T>::linearIn(void* buf, const Box& R, const Interval& comps)
{
  unsigned char* buffer = (unsigned char*)buf;
  m_regFAB.linearIn(buffer, R, comps);
  buffer+= m_regFAB.size(R, comps);
  m_irrFAB.linearIn(buffer, R, comps);
}

/*************/
/*************/
template <class T>
const IntVectSet&
BaseEBCellFAB<T>::getMultiCells() const
{
  CH_assert(isDefined());
  return m_irrFAB.getIVS();
}

/*************/
/*************/
template <class T>
BaseEBCellFAB<T>::BaseEBCellFAB()
{
  setDefaultValues();
}

/*************/
/*************/
template <class T>
BaseEBCellFAB<T>::BaseEBCellFAB(const EBISBox& a_ebisBox,
                                const Box& a_region, int a_nVar)
{
  setDefaultValues();
  define(a_ebisBox, a_region, a_nVar);
}

/*************/
/*************/
template <class T>
void
BaseEBCellFAB<T>::define(const EBISBox& a_ebisBox,
                         const Box& a_region, int a_nVar)
{
  clear();
  m_isDefined = true;
  CH_assert(a_region.cellCentered());
  CH_assert(a_nVar > 0);
  CH_assert(!a_region.isEmpty());
  m_region = a_region & a_ebisBox.getRegion();
  //CH_assert(!m_region.isEmpty());
  m_nComp = a_nVar;
  m_ebisBox = a_ebisBox;

  //make reg fab the full size of the box.
  //this makes ebalias possible
  m_regFAB.resize(a_region, m_nComp);
  //add all the multi-valued
  //cells to m_multiValued
  m_multiCells =
    m_ebisBox.getMultiCells(m_region);

  m_hasMultiCells = !m_multiCells.isEmpty();

  //define irregular fab to be over multi-valued cells.
  m_irrFAB.define(m_multiCells, m_ebisBox.getEBGraph(), m_nComp);
}

/*************/
/*************/
template <class T>
BaseEBCellFAB<T>::
~BaseEBCellFAB()
{
  clear();
}

/*************/
/*************/
template <class T>
void
BaseEBCellFAB<T>::clear()
{
  m_irrFAB.clear();
  m_regFAB.clear();
  m_multiCells.makeEmpty();
  m_nComp = 0;
  m_isDefined = false;
}
/*************/
/*************/
template <class T>
void
BaseEBCellFAB<T>::setVal(const T& value)
{
  m_irrFAB.setVal(value);
  m_regFAB.setVal(value);
}

/*************/
template <class T>
void
BaseEBCellFAB<T>::setVal(int ivar,const  T& value)
{
  m_irrFAB.setVal(ivar, value);
  m_regFAB.setVal(value, m_regFAB.box(), ivar, 1);
}
/*************/
template <class T>
void
BaseEBCellFAB<T>::setVal(const T&   a_value,
                         const Box& a_box,
                         int        a_nstart,
                         int        a_numcomp)
{
  CH_assert(isDefined());
  if (!a_box.isEmpty())
    {
      m_irrFAB.setVal(a_value,a_box,a_nstart,a_numcomp);
      m_regFAB.setVal(a_value,a_box,a_nstart,a_numcomp);
    }
}

/*************/
/*************/
template <class T>
bool
BaseEBCellFAB<T>::isDefined() const
{
  return (m_isDefined);
}
/*************/
/*************/
template <class T>
int
BaseEBCellFAB<T>::nComp() const
{
  CH_assert(isDefined());
  return m_nComp;
}
/*************/
/*************/
template <class T>
const Box&
BaseEBCellFAB<T>::getRegion() const
{
  CH_assert(isDefined());
  return m_region;
}
/*************/
/*************/
template <class T>
const Box&
BaseEBCellFAB<T>::box() const
{
  return getRegion();
}

/*************/
/*************/
template <class T>
const BaseIVFAB<T>&
BaseEBCellFAB<T>::getMultiValuedFAB() const
{
  CH_assert(isDefined());
  return m_irrFAB;
}

//
/*************/
/*************/
template <class T>
BaseIVFAB<T>&
BaseEBCellFAB<T>::getMultiValuedFAB()
{
  CH_assert(isDefined());
  return m_irrFAB;
}

/*************/
/*************/
template <class T>
const BaseFab<T>&
BaseEBCellFAB<T>::getSingleValuedFAB() const
{
  CH_assert(isDefined());
  return m_regFAB;
}

//
/*************/
/*************/
template <class T>
BaseFab<T>&
BaseEBCellFAB<T>::getSingleValuedFAB()
{
  CH_assert(isDefined());
  return m_regFAB;
}
/*************/

/*************
  Arg a_isKnownMultiValued should be set to:
     -1 if arg a_ndin is known to be in a single-valued cell;
     +1 if arg a_ndin is known to be in a multi-valued cell;
      0 if you're not sure (in which case the expensive IntVectSet::contains()
        function gets called).
*/
template <class T>
const T&
BaseEBCellFAB<T>::operator()(const VolIndex& a_ndin,int  a_nVarLoc,
                             int a_isKnownMultiValued ) const
{
  CH_assert(isDefined());
  CH_assert(!m_ebisBox.isCovered(a_ndin.gridIndex()));
  //CH_assert(m_region.contains(a_ndin.gridIndex())); caught in m_regFAB
  CH_assert((a_nVarLoc >= 0)&&(a_nVarLoc < m_nComp));
  CH_assert( (a_isKnownMultiValued == -1)
    ||    (a_isKnownMultiValued ==  0)
    ||    (a_isKnownMultiValued == +1) );

  const T* returnval = NULL;

  switch( a_isKnownMultiValued )
    {
      case -1:
        returnval = &(m_regFAB(a_ndin.gridIndex(), a_nVarLoc));
        break;        // m_regFAB is a BaseFab
      case +1:
        returnval = &(m_irrFAB(a_ndin,a_nVarLoc));
        break;        // m_irrFAB is a BaseIVFAB
      case 0:
        {
          if ( m_multiCells.contains(a_ndin.gridIndex()) )
            {
              returnval = &(m_irrFAB(a_ndin,a_nVarLoc));
            } else
            {
              returnval = &(m_regFAB(a_ndin.gridIndex(), a_nVarLoc));
            }
        }
    }

  return *returnval;

}

/************
  See comment for const version of this function.
*/
template <class T>
T&
BaseEBCellFAB<T>::operator()(const VolIndex& a_ndin,int  a_nVarLoc,
                             int a_isKnownMultiValued )
{
  CH_assert(isDefined());
  CH_assert(m_region.contains(a_ndin.gridIndex()));
  CH_assert((a_nVarLoc >= 0)&&(a_nVarLoc < m_nComp));
  CH_assert(!m_ebisBox.isCovered(a_ndin.gridIndex()));
  CH_assert( (a_isKnownMultiValued == -1)
    ||    (a_isKnownMultiValued ==  0)
    ||    (a_isKnownMultiValued == +1) );

  T* returnval = NULL;

  if ( (a_isKnownMultiValued  == -1) || m_multiCells.isEmpty())
    {
        returnval = &(m_regFAB(a_ndin.gridIndex(), a_nVarLoc));
    }
  else if (a_isKnownMultiValued  == 1)
    {
        returnval = &(m_irrFAB(a_ndin,a_nVarLoc));
    }
  else if (a_isKnownMultiValued  == 0)
    {
      if ( m_multiCells.contains(a_ndin.gridIndex()) )
        {
          returnval = &(m_irrFAB(a_ndin,a_nVarLoc));
        }
      else
        {
          returnval = &(m_regFAB(a_ndin.gridIndex(), a_nVarLoc));
        }
    }
  else
    {
      MayDay::Error("bogus a_isKnownMultiValued input");
    }

  return *returnval;
}


template <class T>
void
BaseEBCellFAB<T>::fill(T* array, const VolIndex& a_ndin, const Interval& a_comps) const
{
  CH_assert(isDefined());
  CH_assert(m_region.contains(a_ndin.gridIndex()));
  CH_assert(!m_ebisBox.isCovered(a_ndin.gridIndex()));

  if (m_multiCells.contains(a_ndin.gridIndex()))
    {
      const T* val = &(m_irrFAB(a_ndin,a_comps.begin()));
      for (int i=0; i<a_comps.size(); i++)
        {
          array[i] = val[i];
        }

    }
  else
    {
      const T* val = &(m_regFAB(a_ndin.gridIndex(), a_comps.begin()));
      int numPts = m_regFAB.box().numPts();
      for (int i=0; i<a_comps.size(); i++)
      {
        array[i] = val[i*numPts];
      }

    }

}

template <class T>
void
BaseEBCellFAB<T>::assign(const T* array,
                         const VolIndex& a_ndin,
                         const Interval& a_comps)
{
  CH_assert(isDefined());
  CH_assert(m_region.contains(a_ndin.gridIndex()));
  CH_assert(!m_ebisBox.isCovered(a_ndin.gridIndex()));

  if (m_multiCells.contains(a_ndin.gridIndex()))
    {
      T* val = &(m_irrFAB(a_ndin,a_comps.begin()));
      for (int i=0; i<a_comps.size(); i++)
        {
          val[i] = array[i];
        }

    }
  else
    {
      T* val = &(m_regFAB(a_ndin.gridIndex(), a_comps.begin()));
      int numPts = m_regFAB.box().numPts();
      for (int i=0; i<a_comps.size(); i++)
      {
        val[i*numPts] = array[i];
      }

    }

}

  /*************/
/*************/
template <class T>
void
BaseEBCellFAB<T>:: setDefaultValues()
{
  m_isDefined = false;
  m_nComp = 0;
}
/*************/
/*************/
template <class T>
const EBISBox&
BaseEBCellFAB<T>::getEBISBox() const
{
  return m_ebisBox;
}

#include "NamespaceFooter.H"
#endif
