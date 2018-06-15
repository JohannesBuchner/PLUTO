#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CFIVS_H_
#define _CFIVS_H_

#include "REAL.H"
#include "DataIndex.H"
#include "IntVect.H"
#include "Box.H"
#include "IntVectSet.H"
#include "NamespaceHeader.H"

//--Forward declarations

class DisjointBoxLayout;
class ProblemDomain;
class LoHiSide;


/// Internal class to find parts of a box outside the valid region of a level
/**
 *  Typically used to find the ghost cells between this and a coarser level.
 *  These cells are represented by an IntVectSet and, if possible, packed into
 *  a box.  This class should only be used through a CFRegion which breaks down
 *  the parts of a box into regions that have a high probability of being
 *  packed.
 */

class CFIVS
{
public:
  /// Null constructor
  CFIVS();

  /// Destructor
  ~CFIVS();

  /// Explicit define based on an IntVectSet
  //  Using this means the cells to store in the IVS have been identified
  //  elsewhere!
  void define(const IntVectSet & a_IVS);

  /// General define for any box using NeighborIterators
  //  WARNING!  To use this, neighbors must be defined in the DBL.  Otherwise
  //  continue to use the legacy define functions.
  void define(const DataIndex&         a_dataIndex,
              const DisjointBoxLayout& a_grids,
              const Box&               a_ghostBox);

  /// Coarsen the stored cells
  void coarsen(int a_ref);

  /// Returns packed box for when isPacked() is true.
  const Box& packedBox() const;

  /// Returns minbox for when isPacked() is false (but can also be used if true)
  const Box& minBox() const;

  /// Get the set of IntVects
  const IntVectSet& getIVS() const;

  /// Returns "true" if this CF IntVectSet can be represented as a simple Box.
  bool isPacked() const;

  /// Returns "true" if the defined IVS is empty
  bool isEmpty() const;

  /// Returns "true" if the object has been defined
  bool isDefined() const;

protected:
  // For internal use
  void  setDefaultValues();

  /// Decrement counts (during redefine)
  void decrementCounts();

  /// Pack the IVS
  void packIVS();

public:
  static long long s_packCount, s_sparseCount;

protected:
  IntVectSet m_IVS;                   ///< IntVects that are outside the level
  Box m_packedBox;                    ///< If m_packed is true, the Box with
                                      ///< which the coarse-fine IntVectSet can
                                      ///< be represented.  If m_packed is
                                      ///< false, same as m_IVS.minBox().
  bool m_packed;                      ///< Can the coarse-fine IntVectSet be
                                      ///< represented as a Box?  Set to false
                                      ///< if the IVS is empty
  bool m_empty;                       ///< Is the IVS empty - cached here for
                                      ///< fast lookup
  bool m_defined;                     ///< Has this object been defined?


/*============================================================================*/
/** \name Legacy members
 *
 *  New uses of CFIVS should avoid these routines
 *//*=========================================================================*/

public:
//@{
  /// (Deprecated) constructor
  /** The ProblemDomain class has replaced the Box domain
      Replaces Box domain with a ProblemDomain and calls define fn.
  */
  CFIVS(const Box&               a_domain,
        const Box&               a_boxIn,
        const DisjointBoxLayout& a_fineBoxes,
        int                      a_direction,
        Side::LoHiSide           a_hiorlo);

  /// (Deprecated) full constructor
  /** CFRegion should specify the region and call general define for a box
      Calls corresponding define function.
   */
  CFIVS(const ProblemDomain&     a_domain,
        const Box&               a_boxIn,
        const DisjointBoxLayout& a_fineBoxes,
        int                      a_direction,
        Side::LoHiSide           a_hiorlo);

  /// (Deprecated) define function
  /** CFRegion should specify the region and call general define for a box
      Replaces Box domain with a ProblemDomain and calls define fn.
  */
  void define(const Box&               a_domain,
              const Box&               a_boxIn,
              const DisjointBoxLayout& a_fineBoxes,
              int                      a_direction,
              Side::LoHiSide           a_hiorlo);

  /// (Deprecated) full define function
  /** CFRegion should specify the region and call general define for a box
   */
  void define (const ProblemDomain&     a_domain,    /// Problem domain at the fine level
               const Box&               a_boxIn,     /// fine-level grid box
               const DisjointBoxLayout& a_fineBoxes, /// fine-level grids
               int                      a_direction, /// direction of the face of interest
               Side::LoHiSide           a_hiorlo);   /// is the face on the hi or lo side of the box?

  /// (Deprecated) faster define function.  Requires the use of CFStencil::buildPeriodicVector static function.
  /** CFRegion should specify the region and call general define for a box
   */
  void define(const ProblemDomain&     a_domain,
              const Box&               a_boxIn,
              const Vector<Box>&       a_periodicfineBoxes,
              int                      a_direction,
              Side::LoHiSide           a_hiorlo);

  /// (Deprecated) Updated define function that uses NeighborIterators instead of traversing the
  /// entire list of boxes.
  /** CFRegion should specify the region and call general define for a box.
   *  However, at least this one is using the recommended technology, Neighbor
   *  iterators.
   */
  void define(const ProblemDomain&     a_domain,
              const DataIndex&         a_dataIndex,
              const DisjointBoxLayout& a_grids,
              int                      a_direction,
              Side::LoHiSide           a_hiorlo);

  /// Get fine intvects which need to be interpolated.
  /** Renamed to getIVS() since CFIV has a general purpose now.
   *  This will be empty if isEmpty() returns true
  */
  const IntVectSet& getFineIVS() const;
//@}
};

// Returns packed box for when isPacked() is true.
inline const Box& CFIVS::packedBox() const
{
  CH_assert(isPacked());
  return m_packedBox;
}

// Returns minbox for when isPacked() is false (but can also be used if true)
inline const Box& CFIVS::minBox() const
{
  return m_packedBox;
}

//  Get the set of IntVects
inline const IntVectSet&
CFIVS::getIVS() const
{
  return m_IVS;
}

//  Get fine intvects which need to be interpolated (deprecated).
inline const IntVectSet&
CFIVS::getFineIVS() const
{
  return m_IVS;
}

//  Returns "true" if this CF IntVectSet can be represented as a simple Box.
/** Returns 'true' if this coarse-fine IntVectSet can be represented
    as just a Box, ie.  a user can perform a dense computation instead of
    a pointwise calculation using IVSIterator.
*/
inline bool
CFIVS::isPacked() const
{
  CH_assert(m_defined);
  return m_packed;
}

//  Returns "true" if the defined IVS is empty
inline bool
CFIVS::isEmpty() const
{
  CH_assert(m_defined);
  return m_empty;
}

//  Returns "true" if the object has been defined
inline bool
CFIVS::isDefined() const
{
  return m_defined;
}

#include "NamespaceFooter.H"
#endif
