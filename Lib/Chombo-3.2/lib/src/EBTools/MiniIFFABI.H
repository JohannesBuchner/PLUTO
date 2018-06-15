#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MINIIFFABI_H_
#define _MINIIFFABI_H_
#include "MayDay.H"
#include "FaceIterator.H"
#include <typeinfo>
#include "NamespaceHeader.H"

template <class T>
bool MiniIFFAB<T>::s_verbose  = false;

/******************/
template <class T> inline
MiniIFFAB<T>::MiniIFFAB()
{
  setDefaultValues();
}
/******************/
template <class T> inline
MiniIFFAB<T>::~MiniIFFAB()
{
  clear();
}
/******************/
template <class T> inline
MiniIFFAB<T>::
MiniIFFAB(const Box&     a_box,
          const EBGraph& a_ebgraph,
          const int&     a_direction,
          const int&     a_nvarin)
{
  setDefaultValues();
  define(a_box, a_ebgraph, a_direction, a_nvarin);
}
/******************/
template <class T> inline
void
MiniIFFAB<T>::
define(const Box&        a_box,
       const EBGraph&    a_ebgraph,
       const int&        a_direction,
       const int&        a_nvarin)
{
  clear();
  m_isDefined = true;
  CH_assert(a_nvarin > 0);
  CH_assert((a_direction >= 0) && (a_direction < SpaceDim));

  m_direction = a_direction;
  m_ebgraph = a_ebgraph;
  m_nComp = a_nvarin;
  m_box = a_box;

  m_faces = m_ebgraph.getMultiValuedFaces(a_direction, a_box);
  m_nFaces = m_faces.size();

  m_data = new T[m_nComp*m_nFaces];
}
/******************/
template <class T> inline
void
MiniIFFAB<T>::
setVal(const T& a_value)
{
  CH_assert(isDefined());
  for (int ivec = 0; ivec < m_nFaces*m_nComp; ivec++)
    {
      m_data[ivec] = a_value;
    }
}
/******************/
template <class T> inline
void
MiniIFFAB<T>::
setVal(int ivar, const T& a_value)
{
  CH_assert(isDefined());
  CH_assert(ivar >= 0);
  CH_assert(ivar < m_nComp);

  for (int ivec = ivar; ivec < m_nFaces*m_nComp; ivec += m_nComp)
    {
      m_data[ivec] = a_value;
    }
}
/******************/
template <class T> inline
void
MiniIFFAB<T>::
copy(const Box&          a_fromBox,
     const Interval&     a_dstInterval,
     const Box&          a_toBox,
     const MiniIFFAB<T>& a_src,
     const Interval&     a_srcInterval)
{
  CH_assert(isDefined());
  CH_assert(a_src.isDefined());
  CH_assert(a_srcInterval.size() == a_dstInterval.size());
  CH_assert(a_dstInterval.begin() >= 0);
  CH_assert(a_srcInterval.begin() >= 0);
  CH_assert(a_dstInterval.end()   < m_nComp);
  CH_assert(a_srcInterval.end()   < a_src.m_nComp);

  if ((!m_box.isEmpty()) && (!a_src.m_box.isEmpty()))
    {
      CH_assert(a_fromBox == a_toBox);
      Box intBox = a_toBox;
      Box boxIntersect = m_box;
      boxIntersect &=  a_src.m_box;
      boxIntersect &=  intBox;
      int compSize = a_srcInterval.size();

      for (int iface = 0; iface < m_faces.size(); iface++)
        {
          const FaceIndex& face = m_faces[iface];
          const IntVect& ivHi = face.gridIndex(Side::Hi);
          const IntVect& ivLo = face.gridIndex(Side::Lo);
          if (boxIntersect.contains(ivHi) || boxIntersect.contains(ivLo))
            {
              for (int icomp = 0; icomp < compSize; icomp++)
                {
                  int isrccomp = a_srcInterval.begin() + icomp;
                  int idstcomp = a_dstInterval.begin() + icomp;
                  (*this)(face, idstcomp) = a_src(face, isrccomp);
                }
            }
        }
    }
}
/*********/
template <class T> inline
T*
MiniIFFAB<T>::
getIndex(const FaceIndex& a_face, const int& a_comp) const
{
  CH_assert(m_isDefined);
  CH_assert((a_comp >= 0) && (a_comp < this->m_nComp));

  T* dataPtr =  this->m_data;
  bool found = false;
  T* retval=NULL;
  for (int iface = 0; iface < m_faces.size(); ++iface)
    {
      if (a_face == m_faces[iface])
        {
          found = true;
          retval=dataPtr;
        }
      if (found) break;
      dataPtr++;
    }
  if (!found)
    {
      MayDay::Error("MiniIFFAB:index not found");
    }
  retval += a_comp*this->m_nFaces;
  return retval;

}
/*********/
template <class T> inline
void
MiniIFFAB<T>::
clear()
{
  if (m_data != NULL)
    {
      delete[] m_data;
      m_data = NULL;
    }
  m_faces.resize(0);
  setDefaultValues();
}
/*************************/
template <class T> inline
T&
MiniIFFAB<T>::operator() (const FaceIndex& a_ndin,
                          const int& a_comp)
{
  CH_assert(isDefined());
  CH_assert(a_ndin.direction() == m_direction);

  T* dataPtr =   getIndex(a_ndin, a_comp);
  return(*dataPtr);
}
/**************************/
template <class T> inline
const T&
MiniIFFAB<T>::operator() (const FaceIndex& a_ndin,
                          const int& a_comp) const
{
  CH_assert(isDefined());
  CH_assert(a_ndin.direction() == m_direction);

  T* dataPtr =   getIndex(a_ndin, a_comp);
  return(*dataPtr);
}
/******************/
template <class T> inline
T*
MiniIFFAB<T>::dataPtr(const int& a_comp)
{
  CH_assert(a_comp >= 0);
  CH_assert(a_comp <= m_nComp);
  return m_data + a_comp*m_nFaces;
}
/******************/
template <class T> inline
const T*
MiniIFFAB<T>::dataPtr(const int& a_comp) const
{
  CH_assert(a_comp >= 0);
  CH_assert(a_comp <= m_nComp);
  return m_data + a_comp*m_nFaces;
}
/******************/
template <class T> inline
void
MiniIFFAB<T>::setDefaultValues()
{
  m_isDefined = false;
  m_data = NULL;
  m_nFaces = 0;
  m_nComp = 0;
  m_direction = -1;
}
/******************/
template <class T> inline
void
MiniIFFAB<T>::
getFaceSubset(Vector<FaceIndex>& a_subFaces, const Box& a_subBox) const
{
  a_subFaces.resize(0);
  for (int iface = 0; iface < m_faces.size(); iface++)
    {
      if ((a_subBox.contains(m_faces[iface].gridIndex(Side::Lo))) ||
         (a_subBox.contains(m_faces[iface].gridIndex(Side::Hi))))
        {
          a_subFaces.push_back(m_faces[iface]);
        }
    }
}
/*********************/
template <class T> inline
int
MiniIFFAB<T>::
size(const Box&      a_region,
     const Interval& a_comps) const
{
  CH_assert(isDefined());
  //create set of cells in fab that are also in the input region
  Vector<FaceIndex> faces;
  getFaceSubset(faces, a_region);

  //account for list of faces
  int facesize = linearListSize(faces);

  //add for each data pt
  int datasize = 0;
  for (int iface = 0; iface < faces.size(); iface++)
    {
      for (int icomp = a_comps.begin(); icomp <= a_comps.end(); icomp++)
        {
          const T& dataPt = (*this)(faces[iface], icomp);
          int pointsize = CH_XD::linearSize(dataPt);
          datasize += pointsize;
        }
    }

  int retval = facesize + datasize;

  if (s_verbose)
    {
      pout() << "###############" << std::endl;
      pout() << "MiniIFFAB size() for region " << m_ebgraph.getRegion() << std::endl;
      pout() << " a_comps  = ("   << a_comps.begin() << "," << a_comps.end()  << ")," << endl;
      pout() << "m_box = " << m_box;
      pout() << "size = " << retval << std::endl;
      pout() << "###############" << std::endl;
    }

  return retval;
}
/********************/
template <class T> inline
void
MiniIFFAB<T>::
linearOut(void* a_buf,
          const Box& a_region,
          const Interval& a_comps) const
{
  CH_assert(isDefined());
  Vector<FaceIndex> faces;
  getFaceSubset(faces, a_region);

  //output the faces.
  unsigned char* charbuffer = (unsigned char *) a_buf;
  linearListOut(charbuffer, faces);
  charbuffer += linearListSize(faces);

  //output the data
  const MiniIFFAB<T>& thisFAB = *this;
  for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      for (int icomp = a_comps.begin(); icomp <= a_comps.end(); icomp++)
        {
          //output the data into the buffer
          const T& dataPt =  thisFAB(face, icomp);
          CH_XD::linearOut(charbuffer, dataPt);
          //increment the buffer offset
          charbuffer += CH_XD::linearSize(dataPt);
        }
    }
}
/********************/
template <class T> inline
void
MiniIFFAB<T>::
linearIn(void* a_buf, const Box& a_region, const Interval& a_comps)
{
  CH_assert(isDefined());

  unsigned char* charbuffer = (unsigned char *) a_buf;
  Vector<FaceIndex> faces;
  linearListIn(faces, charbuffer);
  charbuffer += linearListSize(faces);

  for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      for (int icomp = a_comps.begin(); icomp <= a_comps.end(); icomp++)
        {
          //input the data
          T& dataPt = (*this)(face, icomp);
          CH_XD::linearIn(dataPt, charbuffer) ;
          //increment the buffer offset
          charbuffer += CH_XD::linearSize(dataPt);
        }
    }
}
/********************/
/******************/

#include "NamespaceFooter.H"
#endif
