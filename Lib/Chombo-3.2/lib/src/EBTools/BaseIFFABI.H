#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEIFFABI_H_
#define _BASEIFFABI_H_
#include "MayDay.H"
#include "FaceIterator.H"
#include <typeinfo>
#include "NamespaceHeader.H"

template <class T>
bool BaseIFFAB<T>::s_doSurroundingNodeSemantic = false;

//appalling hack to turn off surrounding node semantic
template <class T>
void
BaseIFFAB<T>::setSurroundingNodeSemantic(bool a_useSurr)
{
  s_doSurroundingNodeSemantic = a_useSurr;
}

//this if statement is to keep the above box growage from
//forcing more copying than the semantic demands.
template <class T>
bool
BaseIFFAB<T>::useThisFace(const Box& a_toBox,
                          const FaceIndex& a_face) const
{
  return ((a_toBox.contains(a_face.gridIndex(Side::Lo)))  ||
          (a_toBox.contains(a_face.gridIndex(Side::Hi))));
}

//this box growing obfuscation is because
//copy and linearization have a cell-centered semantic.
//and we want the surrounding nodes
//of the cells to be copied
//Let's see how many more times I have to rewrite this function.
//--dtg 4/2009
template <class T>
void
BaseIFFAB<T>::getBoxAndIVS(Box&        a_intBox,
                           IntVectSet& a_ivsIntersect,
                           const Box&  a_toBox) const
{
  Box grownBox = a_toBox;
  if (s_doSurroundingNodeSemantic)
    {
      grownBox.grow(m_direction, 1);
    }
  a_intBox       = grownBox;
  a_ivsIntersect = m_ivs;
  a_ivsIntersect &=  a_intBox;
}

template <class T>
bool BaseIFFAB<T>::s_verbose  = false;
/******************/
template <class T> inline   void
BaseIFFAB<T>::setVerbose(bool a_verbose)
{
  s_verbose = a_verbose;
}
/******************/
template <class T> inline
const EBGraph&
BaseIFFAB<T>::getEBGraph() const
{
  return m_ebgraph;
}
/******************/
template <class T> inline
BaseIFFAB<T>::BaseIFFAB()
{
  setDefaultValues();
}
/******************/
template <class T> inline
BaseIFFAB<T>::~BaseIFFAB()
{
  clear();
}
/******************/
template <class T> inline
BaseIFFAB<T>::BaseIFFAB(const IntVectSet& a_ivsin,
                        const EBGraph& a_ebgraph,
                        const int& a_direction,
                        const int& a_nvarin)
{
  setDefaultValues();
  define(a_ivsin, a_ebgraph, a_direction, a_nvarin);
}

template <class T> Arena* BaseIFFAB<T>::s_Arena = NULL;

/******************/
template <class T> inline
void
BaseIFFAB<T>::define(const IntVectSet& a_ivsin,
                     const EBGraph& a_ebgraph,
                     const int& a_direction,
                     const int& a_nvarin)
{
  clear();
  m_isDefined = true;
  CH_assert(a_nvarin > 0);
  CH_assert((a_direction >= 0) && (a_direction < SpaceDim));
  const ProblemDomain& domain = a_ebgraph.getDomain();
  m_direction = a_direction;
  m_ivs = a_ivsin;
  m_ebgraph = a_ebgraph;
  m_nComp = a_nvarin;
  Box minbox = a_ivsin.minBox();
  BaseFab<int> faceLocFab;
  BaseFab<bool> doneThat;
  if (!a_ivsin.isEmpty())
    {
      m_fab.resize(     surroundingNodes(minbox,a_direction), 1);
      faceLocFab.resize(surroundingNodes(minbox,a_direction), 1);
      doneThat.resize(  surroundingNodes(minbox,a_direction), 1);
      faceLocFab.setVal(-1);
      m_fab.setVal(NULL);
      doneThat.setVal(false);
    }
  //get the data size and save offsets
  //need to do this to avoid double-counting faces
  m_nFaces = 0;
  int nVoFs = 0;
  for (IVSIterator ivsit(m_ivs); ivsit.ok(); ++ivsit)
    {
      const IntVect& iv = ivsit();
      IntVect ivLo = iv-BASISV(m_direction);
      IntVect ivHi = iv+BASISV(m_direction);
      //check to see if on domain boundary
      int numVoFsHere = m_ebgraph.numVoFs(iv);

      nVoFs += numVoFsHere;
      //we hold the mapping in the high vof of the faces
      if (!doneThat(iv, 0))
        {
          doneThat(iv, 0) = true;

          bool onLoBoundary = false;
          if (!domain.isPeriodic(m_direction))
            {
              onLoBoundary = iv[m_direction] == domain.domainBox().smallEnd(m_direction);
            }
          int numFacesLoSide = numVoFsHere;
          if (!onLoBoundary)
            {
              int numVoFsLo = m_ebgraph.numVoFs(ivLo);
              numFacesLoSide = numVoFsHere*numVoFsLo;
            }

          faceLocFab(iv, 0) = m_nFaces;
          m_nFaces += numFacesLoSide;
        }
      if (!doneThat(ivHi, 0))
        {
          doneThat(ivHi, 0) = true;

          bool onHiBoundary = false;
          if (!domain.isPeriodic(m_direction))
            {
              onHiBoundary = iv[m_direction] == domain.domainBox().bigEnd(m_direction);
            }
          int numFacesHiSide = numVoFsHere;
          if (!onHiBoundary)
            {
              int numVoFsHi = m_ebgraph.numVoFs(ivHi);
              numFacesHiSide = numVoFsHere*numVoFsHi;
            }

          faceLocFab(ivHi, 0) = m_nFaces;
          m_nFaces += numFacesHiSide;
        }
    }
  // allocate the memory
  if (m_nFaces > 0)
    {
      // const IntVectSet& multi = m_ebgraph.getMultiCells( minbox);
      // pout()<< minbox << "nVoFs:"<<nVoFs<<" nFaces:"<<m_nFaces<<"\n";
#ifdef CH_USE_MEMORY_TRACKING
      // if (s_Arena == NULL)
      //   {
      //     s_Arena = new BArena( (typeid(T)).name());
      //   }
#else
      // if (s_Arena == NULL)
      //   {
      //     s_Arena = new BArena("");
      //   }
#endif

      // if (s_Arena == NULL)
      //   {
      //     MayDay::Error("malloc in basefab failed");
      //   }

      m_truesize = m_nFaces*m_nComp;

      // Note: clear() was called above so this isn't a memory leak
      //m_data = new T[m_truesize];
      m_data.resize(m_truesize);
      // m_data = static_cast<T*>(s_Arena->alloc(m_truesize * sizeof(T)));

#ifdef CH_USE_MEMORY_TRACKING
      // s_Arena->bytes += m_truesize * sizeof(T) + sizeof(BaseIFFAB<T>);
      // if (s_Arena->bytes > s_Arena->peak)
      //   {
      //     s_Arena->peak = s_Arena->bytes;
      //   }
#endif
      doneThat.setVal(false);
    }
  if (m_nFaces==0) return;
  //cache pointers to the first component of the first face into m_fab for fast indexing
  T* startLoc = &m_data[0];
  for (IVSIterator ivsit(m_ivs); ivsit.ok(); ++ivsit)
    {
      const IntVect& iv = ivsit();
      if (!doneThat(iv, 0))
        {
          doneThat(iv, 0) = true;
          int iface = faceLocFab(iv, 0);
          if (iface >= 0)
            {
              m_fab(iv, 0) = startLoc + iface;
            }
        }
      IntVect ivHi = iv+BASISV(m_direction);
      if (!doneThat(ivHi, 0))
        {
          doneThat(ivHi, 0) = true;
          int iface = faceLocFab(ivHi, 0);
          if (iface >= 0)
            {
              m_fab(ivHi, 0) = startLoc + iface;
            }
        }
    }
}
/******************/
template <class T> inline
const IntVectSet&
BaseIFFAB<T>::getIVS() const
{
  return m_ivs;
}
/******************/
template <class T> inline
int
BaseIFFAB<T>::direction() const
{
  return m_direction;
}
/******************/
template <class T> inline
void
BaseIFFAB<T>::setVal(const T& a_value)
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
BaseIFFAB<T>::setVal(int ivar, const T& a_value)
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
BaseIFFAB<T>::copy(const Box& a_fromBox,
                   const Interval& a_dstInterval,
                   const Box& a_toBox,
                   const BaseIFFAB<T>& a_src,
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

      Box intBox;
      IntVectSet ivsIntersect;
      //this does a grow and intersect because
      //copy has a cell centered semantic
      getBoxAndIVS(intBox, ivsIntersect, a_toBox);

      ivsIntersect &=  a_src.m_ivs;
      int compSize = a_srcInterval.size();

      FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;
      FaceIterator faceit(ivsIntersect, m_ebgraph, m_direction, stopCrit);
      for (faceit.reset(); faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          //this if statement is to keep the above box growage from
          //forcing more copying than the semantic demands.
          if (useThisFace(a_toBox, face))
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
BaseIFFAB<T>::getIndex(const FaceIndex& a_face, const int& a_comp) const
{
  //which face is this in the local vector (lexi on vof cell indices)

  T* dataPtr =   (T*)(m_fab(a_face.gridIndex(Side::Hi), 0));
  int faceLoc = getLocalVecIndex(a_face);
  dataPtr += faceLoc;
  dataPtr += m_nFaces*a_comp;
  return dataPtr;
}
/*********/
template <class T> inline
int
BaseIFFAB<T>::getLocalVecIndex(const FaceIndex& a_face) const
{
  int retval;
  if (!a_face.isBoundary())
    {
      //this checks that both vofs are defined
      CH_assert(a_face.cellIndex(Side::Lo) >= 0);
      CH_assert(a_face.cellIndex(Side::Hi) >= 0);
      int xlen = m_ebgraph.numVoFs(a_face.gridIndex(Side::Lo));
      int loCell = a_face.cellIndex(Side::Lo);
      int hiCell = a_face.cellIndex(Side::Hi);
      retval =  loCell + xlen*hiCell;
    }
  else
    {
      int loCell = a_face.cellIndex(Side::Lo);
      int hiCell = a_face.cellIndex(Side::Hi);
      //one should be -1, the other should be larger
      CH_assert(((loCell == -1)&&(hiCell > -1))||
                ((hiCell == -1)&&(loCell > -1)));
      //return the one that is not -1
      retval = Max(loCell, hiCell);
    }
  return retval;
}
/********************/
template <class T> inline
void
BaseIFFAB<T>::clear()
{
  m_nComp = 0;
  m_nFaces = 0;
  m_direction = -1;
  m_ivs.makeEmpty();
  m_fab.clear();
//   if (m_data != NULL)
//     {
//       delete[] m_data;
// #undef free
//       // s_Arena->free(m_data);
// #ifdef CH_USE_MEMORY_TRACKING
//       // s_Arena->bytes -= m_truesize * sizeof(T) + sizeof(BaseIFFAB<T>);
// #endif
//       m_data = NULL;
//    }
  m_isDefined = false;
}
/*************************/
template <class T> inline
bool
BaseIFFAB<T>::isDefined() const
{
  return (m_isDefined);
}
/*************************/
template <class T> inline
int
BaseIFFAB<T>::numFaces() const
{
  return m_nFaces;
}
/*************************/
template <class T> inline
int
BaseIFFAB<T>::nComp() const
{
  return m_nComp;
}
/*************************/
template <class T> inline
T&
BaseIFFAB<T>::operator() (const FaceIndex& a_ndin,
                          const int& a_comp)
{
  CH_assert(isDefined());
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < m_nComp);
  CH_assert((m_ivs.contains(a_ndin.gridIndex(Side::Lo)) ||
             m_ivs.contains(a_ndin.gridIndex(Side::Hi))));
  CH_assert(a_ndin.direction() == m_direction);

  T* dataPtr =   getIndex(a_ndin, a_comp);
  return(*dataPtr);
}
/**************************/
template <class T> inline
const T&
BaseIFFAB<T>::operator() (const FaceIndex& a_ndin,
                          const int& a_comp) const
{
  CH_assert(isDefined());
  CH_assert(a_comp >= 0);
  CH_assert(a_comp < m_nComp);
  CH_assert((m_ivs.contains(a_ndin.gridIndex(Side::Lo)) ||
             m_ivs.contains(a_ndin.gridIndex(Side::Hi))));
  CH_assert(a_ndin.direction() == m_direction);

  T* dataPtr =   getIndex(a_ndin, a_comp);
  return(*dataPtr);
}
/******************/
template <class T> inline
T*
BaseIFFAB<T>::dataPtr(const int& a_comp)
{
  CH_assert(a_comp >= 0);
  CH_assert(a_comp <= m_nComp);
  static T* dummy = NULL;
  if (m_nFaces==0) return dummy;
  //return m_data + a_comp*m_nFaces;
  return &(m_data[0]) + a_comp*m_nFaces;

}
/******************/
template <class T> inline
const T*
BaseIFFAB<T>::dataPtr(const int& a_comp) const
{
  CH_assert(a_comp >= 0);
  CH_assert(a_comp <= m_nComp);
  static T* dummy = NULL;
  if (m_nFaces==0) return dummy;
  //return m_data + a_comp*m_nFaces;
  return &(m_data[0]) + a_comp*m_nFaces;
}
/******************/
template <class T> inline
void
BaseIFFAB<T>::setDefaultValues()
{
  m_isDefined = false;
  //m_data = NULL;
  m_data.resize(0);
  m_nFaces = 0;
  m_nComp = 0;
  m_direction = -1;
}
/******************/
template <class T> inline
int BaseIFFAB<T>::size(const Box& a_region,
                       const Interval& a_comps) const
{
  CH_assert(isDefined());

  Box intBox;
  IntVectSet subIVS;
  //this does a grow and intersect because
  //linearization has a cell centered semantic
  getBoxAndIVS(intBox, subIVS, a_region);

  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;
  FaceIterator faceit(subIVS, m_ebgraph, m_direction, stopCrit);
  Vector<FaceIndex> faces = faceit.getVector();

  //account for list of faces
  int facesize = linearListSize(faces);

  //add for each data pt
  int datasize = 0;
  for (int iface = 0; iface < faces.size(); iface++)
    {
      //this makes sure that the grow done in getBoxAndIVS does
      //not add faces outside the surrounding nodes of the region
      if (useThisFace(a_region, faces[iface]))
        {
          for (int icomp = a_comps.begin(); icomp <= a_comps.end(); icomp++)
            {
              const T& dataPt = (*this)(faces[iface], icomp);
              int pointsize = CH_XD::linearSize(dataPt);
              datasize += pointsize;
            }
        }
    }

  int retval = facesize + datasize;

  if (s_verbose)
    {
      pout() << "###############" << std::endl;
      pout() << "BaseIFFAB size() for region " << m_ebgraph.getRegion() << std::endl;
      pout() << " a_comps  = ("   << a_comps.begin() << "," << a_comps.end()  << "),"
             << " a_region = "   << a_region << " subIVS = ";
      // dumpIVS(&subIVS);
      pout() << "m_ivs = ";
      // dumpIVS(&m_ivs);
      pout() << "size = " << retval << std::endl;
      pout() << "###############" << std::endl;
    }

  return retval;
 }
/********************/
template <class T> inline
void BaseIFFAB<T>::linearOut(void* a_buf,
                             const Box& a_region,
                             const Interval& a_comps) const
{
  CH_assert(isDefined());
  Box intBox;
  IntVectSet subIVS;
  //this does a grow and intersect because
  //linearization has a cell centered semantic
  getBoxAndIVS(intBox, subIVS, a_region);

  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;
  FaceIterator faceit(subIVS, m_ebgraph, m_direction, stopCrit);
  Vector<FaceIndex> faces = faceit.getVector();

  //output the faces.
  unsigned char* charbuffer = (unsigned char *) a_buf;
  linearListOut(charbuffer, faces);
  charbuffer += linearListSize(faces);

  //output the data
  const BaseIFFAB<T>& thisFAB = *this;
  for (int iface = 0; iface < faces.size(); iface++)
    {
      //this makes sure that the grow done in getBoxAndIVS does
      //not add faces outside the surrounding nodes of the region
      if (useThisFace(a_region, faces[iface]))
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
}
/********************/
template <class T> inline
void BaseIFFAB<T>::linearIn(void* a_buf, const Box& a_region, const Interval& a_comps)
{
  CH_assert(isDefined());
  //input the faces
  //input the vofs
  Box intBox;
  IntVectSet subIVS;
  //this does a grow and intersect because
  //linearization has a cell centered semantic
  getBoxAndIVS(intBox, subIVS, a_region);

  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;
  FaceIterator faceit(subIVS, m_ebgraph, m_direction, stopCrit);
  Vector<FaceIndex> faces = faceit.getVector();

  unsigned char* charbuffer = (unsigned char *) a_buf;
  linearListIn(faces, charbuffer);
  charbuffer += linearListSize(faces);

  for (int iface = 0; iface < faces.size(); iface++)
    {
      //this makes sure that the grow done in getBoxAndIVS does
      //not add faces outside the surrounding nodes of the region
      if (useThisFace(a_region, faces[iface]))
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
}
/********************/
/******************/

#include "NamespaceFooter.H"
#endif
