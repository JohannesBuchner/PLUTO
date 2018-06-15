#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CODIMBOXIMPLEM_H_
#define _CODIMBOXIMPLEM_H_

#include "CodimBox.H"

#include "NamespaceHeader.H"

/*******************************************************************************
 *
 * \file CodimBox.cpp
 *
 * \brief Non-inline definitions for CodimBox class
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*
 * Definitions of static member data (these must be in the same file
 * as the constructor definitions)
 *--------------------------------------------------------------------*/

template<typename S> int CodimBox<S>::numOrient[CodimBox<S>::numCD] =
  {
    -1
  };
template<typename S> int CodimBox<S>::totOrient[CodimBox<S>::numCD] =
  {
    -1
  };
template<typename S> unsigned CodimBox<S>::bitOrient[CodimBox<S>::numAI] =
  {
    0
  };
template<typename S> int CodimBox<S>::indexFAB[CodimBox<S>::numAI] =
  {
    -1
  };

/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
///  Default constructor
/**  The codimension can be set by the 'define' function
 *//*-----------------------------------------------------------------*/

template<typename S> CodimBox<S>::CodimBox()
{
  initialize();
}

/*--------------------------------------------------------------------*/
///  Full constructor
/**  \param[in]  a_codim
 *                      The codimension intended for the object
 *   \param[in]  a_box  The cell centered box for the FArrays
 *   \param[in]  a_nvar The number of components for the FArrays
 *//*-----------------------------------------------------------------*/

template<typename S> CodimBox<S>::CodimBox(const int a_codim, const Box& a_box, const int a_nvar)
:
m_codim(a_codim)
{
  initialize();
  define(a_box, a_nvar);
}

/*--------------------------------------------------------------------*/
///  Destructor
/*--------------------------------------------------------------------*/

template<typename S> CodimBox<S>::~CodimBox()
{
  clear();
}

/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
///  Initialize static lookup tables
/**  Construct lookup tables that give the index of the appropriate
 *   S for a given codimension and orientation.  This function
 *   is only called once to initialize the static lookup tables.  The
 *   constructors will call initialize but if any static members are
 *   used before construction, it must be called manually.  See class
 *   CodimBox for more information on the lookup tables.
 *//*-----------------------------------------------------------------*/

template<typename S> void CodimBox<S>::initialize()
{
  if (indexFAB[0] == -1)
    {
      numOrient[0] = 1;             // A codimension 0 object has 1 orientation.
      indexFAB[0] = 0;              // A codimension 0 object has 1 array with
                                    // index 0
      for (int i = 1; i != numCD; ++i) numOrient[i] = 0;
      // The following temporary is used to save the global orientation
      // sequentially for each codimension.
      int tmpCodim[numCD][numAI];
      tmpCodim[0][0] = 0;           // A codimension 0 object has global
                                    // orientation index 0
      // The bits in 'i' described the directions that a codimension object
      // is orthogonal to.  In other words, they describe the orientation of
      // the codimensional object.
      for (int i = 1; i != numAI; ++i)
        {
          // Count the number of 1 bits in 'i'.  This is the codimension
          // and is stored in j.
          unsigned k = i;
          int j = (k & 1);
          while (k >>= 1)
            if (k & 1) ++j;
          // The object with codimension 'j' and orientation 'i' is
          // associated with the FArray at index 'numOrient[j]'.  This is
          // stored for lookup in 'indexFAB'
          indexFAB[i] = numOrient[j];
          // Save the orientation
          tmpCodim[j][numOrient[j]] = i;
          ++numOrient[j];
        }
      // 'numOrient[i]' now has the number of different orientations of a
      // codimension 'i' object.
      // Now update totOrient to a the total orientation from all small
      // codimensiona and update bitOrient to sequentially store the
      // orientations for each codimension
      int k = 0;                    // Sequential index
      int cOrient = 0;              // Count of total orientations
      for (int iCodim = 0; iCodim != numCD; ++iCodim)
        {
          const int nOrient = numOrient[iCodim];
          totOrient[iCodim] = cOrient;
          cOrient += nOrient;
          for (int i = 0; i != nOrient; ++i)
            bitOrient[k++] = tmpCodim[iCodim][i];
        }
    }
}

/*--------------------------------------------------------------------*/
///  Define function allocates space
/**  \param[in]  a_codim
 *                      The codimension intended for the object
 *   \param[in]  a_box  The cell centered box for the FArrays
 *   \param[in]  a_nvar The number of components for the FArrays
 *//*-----------------------------------------------------------------*/

template<typename S> void CodimBox<S>::define(const int a_codim, const Box& a_box, const int a_nvar)
{
  CH_assert(a_codim >= 0);
  CH_assert(a_codim <= SpaceDim);
  m_codim = a_codim;
  define(a_box, a_nvar);
}

/*--------------------------------------------------------------------*/
///  Define function allocates space
/**  Codimension set by constructor
 *   \param[in]  a_box  The cell centered box for the FArrays
 *   \param[in]  a_nvar The number of components for the FArrays
 *//*-----------------------------------------------------------------*/

template<typename S> void CodimBox<S>::define(const Box& a_box, const int a_nvar)
{
  CH_assert(a_nvar > 0);
  m_box = a_box;
  m_nvar = a_nvar;
  clear();
  const int nOrient = numOrient[m_codim];
  m_S.resize(nOrient);
  for (int iOrient = 0; iOrient < nOrient; ++iOrient)
    {
      const Box cdbox = orientBox(iOrient);
      m_S[iOrient] = new S(cdbox, m_nvar);
    }
}

/*--------------------------------------------------------------------*/
///  Deallocation of all memory
/*--------------------------------------------------------------------*/

template<typename S> void CodimBox<S>::clear()
{
  if (m_S.size() != 0)
    {
      const int nOrient = numOrient[m_codim];
      for (int iOrient = 0; iOrient < nOrient; ++iOrient)
        {
          if (m_S[iOrient] != NULL)
            delete m_S[iOrient];
        }
      m_S.clear();
    }
}

/*--------------------------------------------------------------------*/
///  Get the 'i'th orthogonal direction of an orientation
/**  \param[in]  a_iOrient
 *                      The sequential orientation index from which to
 *                      get the directions
 *   \param[in]  a_i    The 'i'th direction.  0 <= a_i < m_codim
 *   \return            An integer defining the requested direction
 *                      0 <= return < SpaceDim
 *   \note
 *   <ul>
 *     <li> The directions are orderd from smallest to largest so
 *          a_i = 0 always returns the smallest direction.
 *   </ul>
 *//*-----------------------------------------------------------------*/

template<typename S> int CodimBox<S>::getDirection(const int a_iOrient, int a_i) const
{
  CH_assert(a_i < m_codim);  // Otherwise may loop forever
  unsigned bOrient = seq2bit(a_iOrient);
  int dir = 0;
  ++a_i;
  while (a_i)
    {
      if (bOrient & 1)
        --a_i;
      bOrient >>= 1;
      ++dir;
    }
  return dir - 1;
}

/*--------------------------------------------------------------------*/
///  Returns S in the directions defined by an IndexType
/**  \param[in]  a_ixType
 *                      Box index type which must define exactly
 *                      m_codim directions as 'NODE'
 *//*-----------------------------------------------------------------*/

/// Returns S in the directions defined by an IndexType
template<typename S> const S& CodimBox<S>::operator()(const IndexType &a_ixType) const
{
  // Bit representation should be the same as used internally in IndexType but
  // do a translation anyways
  int count = 0;
  unsigned bOrient = 0;
  for (int i = 0; i != SpaceDim; ++i)
    {
      if (a_ixType[i])
        {
          bOrient += (1 << i);
          ++count;
        }
    }
  CH_assert(count == m_codim);
  return getS(bOrient);
}

template<typename S> S& CodimBox<S>::operator()(const IndexType &a_ixType)
{
  // Bit representation should be the same as used internally in IndexType but
  // do a translation anyways
  int count = 0;
  unsigned bOrient = 0;
  for (int i = 0; i != SpaceDim; ++i)
    {
      if (a_ixType[i])
        {
          bOrient += (1 << i);
          ++count;
        }
    }
  CH_assert(count == m_codim);
  return getS(bOrient);
}

/*====================================================================*/
/**  \name Accessor functions (1 for each codimension)
 *   \param[in] a_dir0-a_dir5
 *                      The orthogonal directions define the
 *                      orientation of the codimensional geometry.
 *                      Directions can be specified in any order.  For
 *                      example, the edges in 3D orthogonal to
 *                      directions i and j may be specified as (0, 1)
 *                      or (1, 0).  The individual directions should
 *                      be unique but assertions are only implemented
 *                      for codimensions <= 3.
 *//*=================================================================*/
//@{

/// Returns S in the given direction (codimension 0 or SpaceDim)
/** Formally this operator is for codimension 0 objects but it is also
 *  allowed for codimension 'SpaceDim' objects since there is also only
 *  one codimension object (the vertices) in that case.  A codimension
 *  'SpaceDim' object can also be accessed by specifying all the
 *  directions as orthogonal.
 */
template<typename S> const S&CodimBox<S>::operator()() const
{
  CH_assert(SpaceDim >= 0);
  CH_assert(m_codim == 0 || m_codim == SpaceDim);
  return getS((m_codim == 0) ? 0 : numAI - 1);
}

template<typename S> S& CodimBox<S>::operator()()
{
  CH_assert(SpaceDim >= 0);
  CH_assert(m_codim == 0 || m_codim == SpaceDim);
  return getS((m_codim == 0) ? 0 : numAI - 1);
}

/// Returns S in the given direction (codimension 1)
template<typename S> const S&CodimBox<S>::operator()(const int a_dir0) const
{
  CH_assert(SpaceDim >= 1);
  CH_assert(m_codim == 1);
  CH_assert(a_dir0 < SpaceDim);
  return getS((1 << a_dir0));
}

template<typename S> S&CodimBox<S>::operator()(const int a_dir0)
{
  CH_assert(SpaceDim >= 1);
  CH_assert(m_codim == 1);
  CH_assert(a_dir0 < SpaceDim);
  return getS((1 << a_dir0));
}

/// Returns S in the given direction (codimension 2)
template<typename S> const S&CodimBox<S>::operator()(const int a_dir0,
                                            const int a_dir1) const
{
  CH_assert(SpaceDim >= 2);
  CH_assert(m_codim == 2);
  CH_assert(a_dir0 < SpaceDim);
  CH_assert(a_dir1 < SpaceDim);
  CH_assert(a_dir0 != a_dir1);
  return getS((1 << a_dir0) +
                      (1 << a_dir1));
}

template<typename S> S& CodimBox<S>::operator()(const int a_dir0,
                                            const int a_dir1)
{
  CH_assert(SpaceDim >= 2);
  CH_assert(m_codim == 2);
  CH_assert(a_dir0 < SpaceDim);
  CH_assert(a_dir1 < SpaceDim);
  CH_assert(a_dir0 != a_dir1);
  return getS((1 << a_dir0) +
                      (1 << a_dir1));
}

/// Returns S in the given direction (codimension 3)
template<typename S> const S&CodimBox<S>::operator()(const int a_dir0,
                                            const int a_dir1,
                                            const int a_dir2) const
{
  CH_assert(SpaceDim >= 3);
  CH_assert(m_codim == 3);
  CH_assert(a_dir0 < SpaceDim);
  CH_assert(a_dir1 < SpaceDim);
  CH_assert(a_dir2 < SpaceDim);
  CH_assert(a_dir0 != a_dir1);
  CH_assert(a_dir1 != a_dir2);
  CH_assert(a_dir2 != a_dir0);
  return getS((1 << a_dir0) +
                      (1 << a_dir1) +
                      (1 << a_dir2));
}

template<typename S> S& CodimBox<S>::operator()(const int a_dir0,
                                            const int a_dir1,
                                            const int a_dir2)
{
  CH_assert(SpaceDim >= 3);
  CH_assert(m_codim == 3);
  CH_assert(a_dir0 < SpaceDim);
  CH_assert(a_dir1 < SpaceDim);
  CH_assert(a_dir2 < SpaceDim);
  CH_assert(a_dir0 != a_dir1);
  CH_assert(a_dir1 != a_dir2);
  CH_assert(a_dir2 != a_dir0);
  return getS((1 << a_dir0) +
                      (1 << a_dir1) +
                      (1 << a_dir2));
}

/// Returns S in the given direction (codimension 4)
template<typename S> const S& CodimBox<S>::operator()(const int a_dir0,
                                                      const int a_dir1,
                                                      const int a_dir2,
                                                      const int a_dir3) const
{
  CH_assert(SpaceDim >= 4);
  CH_assert(m_codim == 4);
  CH_assert(a_dir0 < SpaceDim);
  CH_assert(a_dir1 < SpaceDim);
  CH_assert(a_dir2 < SpaceDim);
  CH_assert(a_dir3 < SpaceDim);
  // No more checks for unique directions
  return getS((1 << a_dir0) +
                      (1 << a_dir1) +
                      (1 << a_dir2) +
                      (1 << a_dir3));
}

template<typename S> S&CodimBox<S>::operator()(const int a_dir0,
                                            const int a_dir1,
                                            const int a_dir2,
                                            const int a_dir3)
{
  CH_assert(SpaceDim >= 4);
  CH_assert(m_codim == 4);
  CH_assert(a_dir0 < SpaceDim);
  CH_assert(a_dir1 < SpaceDim);
  CH_assert(a_dir2 < SpaceDim);
  CH_assert(a_dir3 < SpaceDim);
  // No more checks for unique directions
  return getS((1 << a_dir0) +
                      (1 << a_dir1) +
                      (1 << a_dir2) +
                      (1 << a_dir3));
}

/// Returns S in the given direction (codimension 5)
template<typename S> const S& CodimBox<S>::operator()(const int a_dir0,
                                                      const int a_dir1,
                                                      const int a_dir2,
                                                      const int a_dir3,
                                                      const int a_dir4) const
{
  CH_assert(SpaceDim >= 5);
  CH_assert(m_codim == 5);
  CH_assert(a_dir0 < SpaceDim);
  CH_assert(a_dir1 < SpaceDim);
  CH_assert(a_dir2 < SpaceDim);
  CH_assert(a_dir3 < SpaceDim);
  CH_assert(a_dir4 < SpaceDim);
  // No more checks for unique directions
  return getS((1 << a_dir0) +
                      (1 << a_dir1) +
                      (1 << a_dir2) +
                      (1 << a_dir3) +
                      (1 << a_dir4));
}

template<typename S> S& CodimBox<S>::operator()(const int a_dir0,
                                                const int a_dir1,
                                                const int a_dir2,
                                                const int a_dir3,
                                                const int a_dir4)
{
  CH_assert(SpaceDim >= 5);
  CH_assert(m_codim == 5);
  CH_assert(a_dir0 < SpaceDim);
  CH_assert(a_dir1 < SpaceDim);
  CH_assert(a_dir2 < SpaceDim);
  CH_assert(a_dir3 < SpaceDim);
  CH_assert(a_dir4 < SpaceDim);
  // No more checks for unique directions
  return getS((1 << a_dir0) +
                      (1 << a_dir1) +
                      (1 << a_dir2) +
                      (1 << a_dir3) +
                      (1 << a_dir4));
}

/// Returns S in the given direction (codimension 6)
template<typename S> const S& CodimBox<S>::operator()(const int a_dir0,
                                            const int a_dir1,
                                            const int a_dir2,
                                            const int a_dir3,
                                            const int a_dir4,
                                            const int a_dir5) const
{
  CH_assert(SpaceDim >= 6);
  CH_assert(m_codim == 6);
  CH_assert(a_dir0 < SpaceDim);
  CH_assert(a_dir1 < SpaceDim);
  CH_assert(a_dir2 < SpaceDim);
  CH_assert(a_dir3 < SpaceDim);
  CH_assert(a_dir4 < SpaceDim);
  CH_assert(a_dir5 < SpaceDim);
  // No more checks for unique directions
  return getS((1 << a_dir0) +
                      (1 << a_dir1) +
                      (1 << a_dir2) +
                      (1 << a_dir3) +
                      (1 << a_dir4) +
                      (1 << a_dir5));
}

template<typename S> S& CodimBox<S>::operator()(const int a_dir0,
                                            const int a_dir1,
                                            const int a_dir2,
                                            const int a_dir3,
                                            const int a_dir4,
                                            const int a_dir5)
{
  CH_assert(SpaceDim >= 6);
  CH_assert(m_codim == 6);
  CH_assert(a_dir0 < SpaceDim);
  CH_assert(a_dir1 < SpaceDim);
  CH_assert(a_dir2 < SpaceDim);
  CH_assert(a_dir3 < SpaceDim);
  CH_assert(a_dir4 < SpaceDim);
  CH_assert(a_dir5 < SpaceDim);
  // No more checks for unique directions
  return getS((1 << a_dir0) +
                      (1 << a_dir1) +
                      (1 << a_dir2) +
                      (1 << a_dir3) +
                      (1 << a_dir4) +
                      (1 << a_dir5));
}
//@}

/*--------------------------------------------------------------------*/
///  Set all values in all FAB
/**  \param[in]  a_x    New value
 *//*-----------------------------------------------------------------*/

template<typename S> void CodimBox<S>::setVal(const Real a_x)
{
  const int nOrient = numOrient[m_codim];
  for (int iOrient = 0; iOrient < nOrient; ++iOrient)
    {
      m_S[iOrient]->setVal(a_x);
    }
}

/*====================================================================*/
/**  \name Functions for compatability with containers like
 *         BoxLayoutData.h
 *//*=================================================================*/
//@{

/*--------------------------------------------------------------------*/
///  Copy from another CodimBox
/**  \param[in]  a_R    Centered box describing region to copy
 *   \param[in]  a_Cd   Interval for components at destination
 *   \param[in]  a_src  The source CodimBox (must be same codimension)
 *   \param[in]  a_Cs   Interval for components at source
 *   \note
 *   <ul>
 *     <li> It is assumed that both CodimBoxes have the same
 *          orientation at their intersection.
 *     <li> a_R is not assumed to be contained within m_box because
 *          higher-level copy routines generally work on cell centers.
 *   </ul>
 *//*-----------------------------------------------------------------*/

template<typename S> void CodimBox<S>::copy(const Box& a_R,
                                      const Interval& a_Cd,
                                      const CodimBox& a_src,
                                      const Interval& a_Cs)
{
  CH_assert(m_codim == a_src.m_codim);
  const int nOrient = numOrient[m_codim];
  for (int iOrient = 0; iOrient != nOrient; ++iOrient)
    {
      CH_assert(m_S[iOrient] != NULL);
      const S& srcFab = *a_src.m_S[iOrient];
      // Orient box in othogonal orientation directions
      Box cdbox = orientBox(iOrient, a_R);
      // Intersect with both source and destination
      cdbox &= srcFab.box();
      cdbox &= m_S[iOrient]->box();
      if (!cdbox.isEmpty())
        {
          m_S[iOrient]->copy(cdbox, a_Cd, cdbox, srcFab, a_Cs);
        }
    }
}

/*--------------------------------------------------------------------*/
///  Copy from another CodimBox but to a different region
/**  \param[in]  a_Rs   Centered box describing region to copy from in
 *                      source
 *   \param[in]  a_Cd   Interval for components at destination
 *   \param[in]  a_Rd   Centered box describing region to copy to in
 *                      destination
 *   \param[in]  a_src  The source CodimBox (must be same codimension)
 *   \param[in]  a_Cs   Interval for components at source
 *   \note
 *   <ul>
 *     <li> It is assumed that both CodimBoxes have the same
 *          orientation at their intersection.
 *     <li> Beware that the source and destination arguments are not
 *          grouped together.
 *     <li> The story here seems to be that copyTo in a LevelData uses
 *          centered boxes to see who to copy to.  So we need +1 ghost
 *          cells in the destination to ensure abutting boxes are
 *          matched and that we can then copy abutting codimension
 *          objects.  Hence, a_Rd may not be contained in m_box and
 *          must be intersected.  Then make assumtions about what
 *          a_Rs should be.
 *   </ul>
 *//*-----------------------------------------------------------------*/

template<typename S> void CodimBox<S>::copy(const Box& a_Rs,
                                      const Interval& a_Cd,
                                      const Box& a_Rd,
                                      const CodimBox<S>& a_src,
                                      const Interval& a_Cs)
{
  CH_assert(m_codim == a_src.m_codim);
  CH_assert(a_Rs.sameSize(a_Rd));
  const int nOrient = numOrient[m_codim];
  for (int iOrient = 0; iOrient != nOrient; ++iOrient)
    {
      CH_assert(m_S[iOrient] != NULL);
      // Orient destination box in othogonal orientation directions
      Box cdboxd = orientBox(iOrient, a_Rd);
      cdboxd &= m_S[iOrient]->box();
      // Now assume a box size for the source
      Box cdboxs = cdboxd;
      if (a_Rs != a_Rd)
        // Use 'cdboxd' but with shift from 'a_Rd' to 'a_Rs'
        {
          IntVect shiftVect(a_Rs.smallEnd());
          shiftVect -= a_Rd.smallEnd();
          cdboxs.shift(shiftVect);
        }
      if (!cdboxd.isEmpty())
        {
          m_S[iOrient]->copy(cdboxs, a_Cd, cdboxd,
                                     *a_src.m_S[iOrient], a_Cs);
        }
    }
}

/*--------------------------------------------------------------------*/
///  Returns size of serialized data
/**  Returns size, in number of bytes, of a 1D representation of data
 *   in the class.
 *   \param[in]  a_R    Centered box describing size
 *   \param[in]  a_comps
 *                      Number of components for describing size
 *//*-----------------------------------------------------------------*/

template<typename S> int CodimBox<S>::size(const Box& a_R, const Interval& a_comps) const
{
  int totalSize = 0;
  S tempFab;
  const int nOrient = numOrient[m_codim];
  for (int iOrient = 0; iOrient != nOrient; ++iOrient)
    {
      const Box cdbox = orientBox(iOrient, a_R);
      totalSize += tempFab.size(cdbox, a_comps);
    }
  return totalSize;
}

/*--------------------------------------------------------------------*/
///  Writes a serialized representation of this CodimBox
/**  Write a linear representaion of the internal data.  Assumes that
 *   sufficient memory for the buffer has already been allocated by
 *   the caller.
 *   \param[in]  a_buf  Buffer for linear representation
 *   \param[in]  a_R    Centered box from which to write
 *   \param[in]  a_comps
 *                      Components to be written
 *//*-----------------------------------------------------------------*/

template<typename S> void CodimBox<S>::linearOut(void* a_buf, const Box& a_R, const Interval& a_comps) const
{
  Real* buffer = static_cast<Real*>(a_buf);
  const int nOrient = numOrient[m_codim];
  for (int iOrient = 0; iOrient != nOrient; ++iOrient)
    {
      CH_assert(m_S[iOrient] != NULL);
      const Box cdbox = orientBox(iOrient, a_R);
      const int orientSize = m_S[iOrient]->size(cdbox, a_comps);
      m_S[iOrient]->linearOut(buffer, cdbox, a_comps);
      buffer += orientSize/sizeof(Real);
    }
}

/*--------------------------------------------------------------------*/
///  Reads a serialized representation of this CodimBox
/**  Reads in the output of linearOut.
 *   \param[in]  a_buf  Buffer for linear representation
 *   \param[in]  a_R    Centered box to read to
 *   \param[in]  a_comps
 *                      Components to be read
 *//*-----------------------------------------------------------------*/

template<typename S> void CodimBox<S>::linearIn(void* a_buf, const Box& a_R, const Interval& a_comps)
{
  Real* buffer = static_cast<Real*>(a_buf);
  const int nOrient = numOrient[m_codim];
  for (int iOrient = 0; iOrient != nOrient; ++iOrient)
    {
      CH_assert(m_S[iOrient] != NULL);
      const Box cdbox = orientBox(iOrient, a_R);
      const int orientSize = m_S[iOrient]->size(cdbox, a_comps);
      m_S[iOrient]->linearIn(buffer, cdbox, a_comps);
      buffer += orientSize/sizeof(Real);
    }
}
//@}

/*==============================================================================
 * Private member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
///  Get all the orthogonal directions of an orientation (general)
/**  General static routine for a given codimension and a bit
 *   representation of the orientation
 *   \param[in]  a_codim
 *                      The codimension
 *   \param[in]  a_bOrient
 *                      Bit representation of the orientation
 *   \param[in]  a_dir  Array with minimum size [a_codim]
 *   \param[out] a_dir  Integers defining the directions
 *                      0 <= dir[i] < SpaceDim
 *   \note
 *   <ul>
 *     <li> The directions are orderd from smallest to largest.
 *   </ul>
 *//*-----------------------------------------------------------------*/

template<typename S> void CodimBox<S>::genGetDirections(const int a_codim,
                                                  int a_bOrient,
                                                  int *const a_dir)
{
  for (int i = 0; i != a_codim; ++i)
    a_dir[i] = 0;
  int iStart = 0;
  while (a_bOrient)
    {
      if (a_bOrient & 1)
        ++iStart;
      a_bOrient >>= 1;
      for (int i = iStart; i < a_codim; ++i)
        ++a_dir[i];
    }
}

/*--------------------------------------------------------------------*/
///  Orient a centered box in the othogonal directions of the given
///  orientation (general)
/**  General static routine for a given codimension and a bit
 *   representation of the orientation
 *   \param[in]  a_bOrient
 *                      Bit representation of the orientation
 *   \param[in]  a_cbox
 *                      A centered box.
 *   \return            The orientated box
 *//*-----------------------------------------------------------------*/

template<typename S> Box CodimBox<S>::genOrientBox(int a_bOrient, Box a_cbox)
{
  int dir = 0;
  while (a_bOrient)
    {
      if (a_bOrient & 1)
        a_cbox.surroundingNodes(dir);
      a_bOrient >>= 1;
      ++dir;
    }
  return a_cbox;
}

// Only of use with a register class
// /*--------------------------------------------------------------------*/
// ///  Shift all FABs by halves in one direction.
// /**  \param[in]  a_dir  Direction to shift
//  *   \param[in]  a_numHalfs
//  *                      Number of halves to shift
//  *   \note
//  *   <ul>
//  *     <li> This really messes with the concept of the CodimBox -
//  *          hence only for private use and the shifting should be
//  *          undone as soon as possible
//  *   </ul>
//  *//*-----------------------------------------------------------------*/

// void
// CodimBox::shiftHalf(const int a_dir, const int a_numHalfs)
// {
//   const int nOrient = numOrient[m_codim];
//   for (int iOrient = 0; iOrient != nOrient; ++iOrient)
//     {
//       m_S[iOrient]->shiftHalf(a_dir, a_numHalfs);
//     }
// }

#include "NamespaceFooter.H"

#endif
