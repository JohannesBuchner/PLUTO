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

#ifndef _BASEEBFACEFABI_H_
#define _BASEEBFACEFABI_H_

#include "BoxIterator.H"
#include "NamespaceHeader.H"

/**********************/
/**********************/
template <class T>
void
BaseEBFaceFAB<T>::setCoveredFaceVal(const T&    a_val,
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
          const IntVect& ivCell = bit();
          for (SideIterator sit; sit.ok(); ++sit)
            {
              const Vector< FaceIndex >& faces = m_ebisBox.getAllFaces(ivCell,m_iDir,sit());
              if ( (faces.size()==0) || (faces.size()>1) )
                {
                  const int& hilo =  sign(sit());
                  IntVect ivFace = ivCell;
                  if (hilo == 1)
                    {
                      ivFace[m_iDir] += 1;
                    }
                  m_regFAB(ivFace,a_comp) = a_val;
                }
            }
        }
    }
}
/*************/
/*************/
template <class T>
int
BaseEBFaceFAB<T>::size(const Box& R, const Interval& comps) const
{

  Box RFace = surroundingNodes(R, m_iDir);
  int retval = m_regFAB.size(RFace, comps);
  retval    += m_irrFAB.size(R,     comps);
  return retval;
}

/*************/
template <class T>
void
BaseEBFaceFAB<T>::linearOut(void* buf, const Box& R, const Interval& comps)
  const
{
  Box RFace = surroundingNodes(R, m_iDir);
  unsigned char* buffer = (unsigned char*)buf;
  m_regFAB.linearOut(buffer, RFace, comps);
  buffer+= m_regFAB.size(RFace, comps);
  m_irrFAB.linearOut(buffer, R, comps);
}
/*************/
template <class T>
void
BaseEBFaceFAB<T>::linearIn(void* buf, const Box& R, const Interval& comps)
{
  Box RFace = surroundingNodes(R, m_iDir);
  unsigned char* buffer = (unsigned char*)buf;
  m_regFAB.linearIn(buffer, RFace, comps);
  buffer+= m_regFAB.size(RFace, comps);
  m_irrFAB.linearIn(buffer, R, comps);
}

/*************/
template <class T> inline
const EBISBox&
BaseEBFaceFAB<T>::getEBISBox() const
{
  return m_ebisBox;
}
/*************/
template <class T>
BaseEBFaceFAB<T>&
BaseEBFaceFAB<T>::copy(const BaseEBFaceFAB<T>& a_src)
{
  CH_assert(isDefined());
  CH_assert(a_src.isDefined());
  CH_assert(m_nComp <= a_src.m_nComp);
  CH_assert(m_region.sameType(a_src.m_region));

  Interval comps(0,m_nComp-1);

  Box overlap(m_region);
  overlap &= a_src.m_region;

  this->copy(overlap,comps,overlap,a_src,comps);

  return *this;
}

/*************/
template <class T>
inline void
BaseEBFaceFAB<T>::copy(const Box& a_regionFrom,
                       const Interval& a_dstInt,
                       const Box& a_regionTo,
                       const BaseEBFaceFAB<T>& a_src,
                       const Interval& a_srcInt)
{
  //CH_assert(a_regionFrom == a_regionTo);
  CH_assert(isDefined());
  CH_assert(a_dstInt.size()==a_srcInt.size());
  CH_assert(a_regionTo.cellCentered());
  CH_assert(a_src.m_iDir ==m_iDir);
  CH_assert(a_dstInt.begin()>=0);
  CH_assert(a_srcInt.begin()>=0);
  CH_assert(a_srcInt.end()< a_src.m_nComp);
  CH_assert(a_dstInt.end()< m_nComp);
  CH_assert(a_src.m_region.contains(a_regionFrom));
  CH_assert(m_region.contains(a_regionTo));
  Box RFromCell = a_regionFrom;
  Box RFromFace = surroundingNodes(RFromCell, m_iDir);
  Box RToCell = a_regionTo;
  Box RToFace = surroundingNodes(RToCell, m_iDir);

  CH_assert(m_regionFace.contains(RToFace));
  CH_assert(a_src.m_regionFace.contains(RFromFace));
  //this one has to be rface because basefab does not know
  //from cell vs edge centered
  m_regFAB.copy(RFromFace, a_dstInt, RToFace, a_src.m_regFAB, a_srcInt);
  //this has to be cell centered because that is how
  //it is defined.
  m_irrFAB.copy(RFromCell, a_dstInt, RToCell, a_src.m_irrFAB, a_srcInt);
}
/**********************/
/**********************/
template <class T> inline
BaseEBFaceFAB<T>::BaseEBFaceFAB()
{
  setDefaultValues();
}

/**********************/
/**********************/
template <class T> inline
BaseEBFaceFAB<T>::BaseEBFaceFAB(const EBISBox& a_ebisBox,
                                const Box& a_region,
                                int a_iDir, int a_nComp)
{
  setDefaultValues();
  define(a_ebisBox, a_region, a_iDir, a_nComp);
}

/**********************/
/**********************/
template <class T> inline
void
BaseEBFaceFAB<T>::define(const EBISBox&  a_ebisBox,
                         const Box& a_region,
                         int a_iDir, int a_nComp)
{
  clear();
  m_isDefined = true;
  CH_assert(a_region.cellCentered());
  CH_assert(!a_region.isEmpty());
  CH_assert(a_nComp > 0);
  CH_assert((a_iDir >= 0) && (a_iDir < SpaceDim));

  m_ebisBox = a_ebisBox;
  m_region =  a_region;
  m_region &= a_ebisBox.getDomain().domainBox();
  CH_assert(!m_region.isEmpty());
  //if this fails, you need an ebisbox with more ghost cells
  CH_assert(a_ebisBox.getRegion().contains(m_region));

  m_nComp = a_nComp;
  m_iDir = a_iDir;
  m_regionFace= surroundingNodes(m_region, a_iDir);
  Box grownRegion = m_region;
  grownRegion.grow(a_iDir, 1);
  grownRegion &= a_ebisBox.getDomain().domainBox();
  //if this fails, you need an ebisbox with more ghost cells
  CH_assert(a_ebisBox.getRegion().contains(grownRegion));

  m_regFAB.resize(m_regionFace, m_nComp);

  m_irrFAB.define(a_region, m_ebisBox.getEBGraph(), m_iDir, m_nComp);
}

/**********************/
/**********************/
template <class T> inline
BaseEBFaceFAB<T>::~BaseEBFaceFAB()
{
  clear();
}

/**********************/
/**********************/
template <class T> inline
void
BaseEBFaceFAB<T>::clear()
{
  m_irrFAB.clear();
  m_regFAB.clear();
  m_isDefined = false;
}

/**********************/
/**********************/
template <class T> inline
void
BaseEBFaceFAB<T>::setVal(const T& value)
{
  m_irrFAB.setVal(value);
  m_regFAB.setVal(value);
}

/**********************/
template <class T> inline
void
BaseEBFaceFAB<T>::setVal(int ivar, const T& value)
{
  m_irrFAB.setVal(ivar,  value);
  m_regFAB.setVal(value, m_regFAB.box(), ivar, 1);
}

/**********************/
/**********************/
template <class T> inline
bool
BaseEBFaceFAB<T>::isDefined() const
{
  return (m_isDefined);
}

/**********************/
/**********************/
template <class T> inline
int
BaseEBFaceFAB<T>::nComp() const
{
  CH_assert(isDefined());
  return m_nComp;
}

/**********************/
/**********************/
template <class T> inline
const MiniIFFAB<T>&
BaseEBFaceFAB<T>::getMultiValuedFAB() const
{
  CH_assert(isDefined());
  return m_irrFAB;
}

/**********************/
/**********************/
template <class T> inline
MiniIFFAB<T>&
BaseEBFaceFAB<T>::getMultiValuedFAB()
{
  CH_assert(isDefined());
  return m_irrFAB;
}

/**********************/
/**********************/
template <class T> inline
const BaseFab<T>&
BaseEBFaceFAB<T>::getSingleValuedFAB() const
{
  CH_assert(isDefined());
  return m_regFAB;
}

/**********************/
/**********************/
template <class T> inline
BaseFab<T>&
BaseEBFaceFAB<T>::getSingleValuedFAB()
{
  CH_assert(isDefined());
  return m_regFAB;
}

/**********************/
/**********************/
template <class T> inline
int
BaseEBFaceFAB<T>::direction() const
{
  CH_assert(isDefined());
  return m_iDir;
}

/**********************/
/**********************/
template <class T> inline
const Box&
BaseEBFaceFAB<T>::getRegion() const
{
  CH_assert(isDefined());
  return m_regionFace;
}
/**********************/
/**********************/
template <class T> inline
const Box&
BaseEBFaceFAB<T>::getCellRegion() const
{
  CH_assert(isDefined());
  return m_region;
}

/**********************/
/**********************/
template <class T> inline
T&
BaseEBFaceFAB<T>::operator()(const FaceIndex& a_facein, int  a_nCompLoc)
{
  CH_assert(isDefined());
  CH_assert((a_nCompLoc >= 0)&&(a_nCompLoc < m_nComp));
  const IntVect& ivhi = a_facein.gridIndex(Side::Hi);
  CH_assert(m_regionFace.contains(ivhi));

  const Vector<FaceIndex>& multiFaces = m_irrFAB.getFaces();
  T* returnval;

  bool isFaceHere = false;
  for (int iface = 0; iface< multiFaces.size(); iface++)
    {
      const FaceIndex& face = multiFaces[iface];
      if (face == a_facein)
        {
          isFaceHere = true;
        }
    }
  if (isFaceHere)
    {
      returnval = &m_irrFAB(a_facein, a_nCompLoc);
    }
  else
    {
      returnval = &m_regFAB(ivhi, a_nCompLoc);
    }
  return *returnval;
}

/**********************/
/**********************/
template <class T> inline
const T&
BaseEBFaceFAB<T>::operator() (const FaceIndex& a_facein, int  a_nCompLoc) const
{
  CH_assert(isDefined());
  CH_assert((a_nCompLoc >= 0)&&(a_nCompLoc < m_nComp));
  const IntVect& ivhi = a_facein.gridIndex(Side::Hi);
  CH_assert(m_regionFace.contains(ivhi));

  const Vector<FaceIndex>& multiFaces = m_irrFAB.getFaces();
  const T* returnval;

  bool isFaceHere = false;
  for (int iface = 0; iface< multiFaces.size(); iface++)
    {
      const FaceIndex& face = multiFaces[iface];
      if (face == a_facein)
        {
          isFaceHere = true;
        }
    }
  if (isFaceHere)
    {
      returnval = &m_irrFAB(a_facein, a_nCompLoc);
    }
  else
    {
      returnval = &m_regFAB(ivhi, a_nCompLoc);
    }
  return *returnval;
}

/**********************/
/**********************/
template <class T> inline
void
BaseEBFaceFAB<T>::setDefaultValues()
{
  m_isDefined = false;
  m_nComp = -1;
  m_iDir = -1;
}

#include "NamespaceFooter.H"
#endif
