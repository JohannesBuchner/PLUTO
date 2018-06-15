#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CODIMBOX_H_
#define _CODIMBOX_H_
#include "BaseFab.H"
#include "FArrayBox.H"
#include "BoxLayoutData.H"
#include "NamespaceHeader.H"

/*******************************************************************************
 */
///  An FArrayBox container for storage on the codimensions of a box.
/**
 *   The orientation of a box is described by specifying the orthogonal
 *   directions and internally, this is converted to a bit representation.
 *   Direction 'i' is bit 2^i (or via shifting as 1 << i).  Directions can be
 *   specified in any order.  For example, the edges in 3D orthogonal to
 *   directions i and j may be specified as (0, 1) or (1, 0).  Internally, both
 *   are converted to 0x3.  Orientations can be described by sequential indices
 *   or bit representations.  The routines seq2bit and bit2seq can be used to
 *   switch between the two using lookup tables.
 *
 *   For example, a 3 dimensional space has the following codimensions, bit
 *   representations, and sequential indices of the orientations:
 *                                                                     \verbatim
 *         Codim   Bit   Integer   Sequential   numOrient   totOrient
 *   Dir           zyx
 *             0   000         0            0           1           0
 *             1   001         1            0           3           1
 *             1   010         2            1
 *             1   100         4            2
 *             2   011         3            0           3           4
 *             2   101         5            1
 *             2   110         6            2
 *             3   111         7            0           1           7
 *                                                                  \endverbatim
 *   The 'indexFAB' and 'bitOrient' are then defined as follows:
 *                                                                     \verbatim
 *     indexFAB:  0 0 1 0 2 1 2 0
 *     bitOrient: 0 1 2 4 3 5 6 7
 *                                                                  \endverbatim
 *  \note
 *  <ul>
 *    <li> CodimBoxFactory is also defined in this file for use with
 *         BoxLayoutData and similar containers
 *    <li> A FluxBox is essentially a CodimBox with codimension 1.
 *    <li> Not much functionality exists in this class beyond creating and
 *         retrieving the FArrayBoxes.
 *    <li> Codimension 0 is the original object spanning 'SpaceDim' dimensions.
 *         Using this is no different from an FArrayBox itself.
 *    <li> Codimension 'SpaceDim' is a vertex.  This is also the same as an
 *         FArrayBox that surrounds all nodes.
 *  </ul>
 *
 ******************************************************************************/

template <class S >class CodimBox
{

public:
  enum
  {
    numCD = SpaceDim + 1,           ///< Number of codimensions
    numAI = 1 << SpaceDim           ///< Total orientations (array indices) for
                                    ///< all codimensions
  };

  /// Default constructor
  CodimBox();

  /// Full constructor
  CodimBox(const int a_codim, const Box& a_box, const int a_nvar);

  /// Destructor
  ~CodimBox();

  /// Initialize static lookup tables
  static void initialize();

  /// Full define function
  void define(const int a_codim, const Box& a_box, const int a_nvar);

  /// Define function
  void define(const Box& a_box, int a_nvar);

  /// Deallocate
  void clear();

  /// Return the codimension
  int getCodim() const;

  /// Number of components
  int nComp() const;

  /// Number of different orientations of codimensional objects
  int getNumOrient() const;

  /// Number of different orientations of codimensional objects.  Used
  /// externally for any codimension.
  static int getNumOrient(const int a_codim);

  /// Switch from a bit orientation to a sequential index of the
  /// orientation
  int bit2seq(const unsigned a_bOrient) const;

  /// Switch from a sequential index to a bit representation of the
  /// orientation
  unsigned seq2bit(const int a_iOrient) const;

  /// Switch from a sequential index to a bit representation of the orientation
  /// for any codimension
  static unsigned seq2bit(const int a_codim, const int a_iOrient);

  /// Get the 'i'th orthogonal direction of an orientation
  int getDirection(const int a_iOrient, int a_i) const;

  /// Get all the orthogonal directions of an orientation
  void getDirections(const int a_iOrient, int *const a_dir) const;

  /// Get all the orthogonal directions of an orientation.  Used externally for
  /// any codimension.
  static void getDirections(const int a_codim,
                            const int a_iOrient,
                            int *const a_dir);

  /// Returns cell-centered box which defines the CodimBox
  const Box& box() const;

  /// Returns S in the directions defined by a Box
  const S& operator()(const Box &a_box) const;
  S& operator()(const Box &a_box);

  /// Returns S in the directions defined by an IndexType
  const S& operator()(const IndexType &a_ixType) const;
  S& operator()(const IndexType &a_ixType);

  /// Returns S in the given direction (codimension 0)
  const S& operator()() const;
  S& operator()();

  /// Returns S in the given direction (codimension 1)
  const S& operator()(const int a_dir0) const;
  S& operator()(const int a_dir0);

  /// Returns S in the given direction (codimension 2)
  const S& operator()(const int a_dir0,
                              const int a_dir1) const;
  S& operator()(const int a_dir0,
                        const int a_dir1);

  /// Returns S in the given direction (codimension 3)
  const S& operator()(const int a_dir0,
                              const int a_dir1,
                              const int a_dir2) const;
  S& operator()(const int a_dir0,
                        const int a_dir1,
                        const int a_dir2);

  /// Returns S in the given direction (codimension 4)
  const S& operator()(const int a_dir0,
                              const int a_dir1,
                              const int a_dir2,
                              const int a_dir3) const;
  S& operator()(const int a_dir0,
                        const int a_dir1,
                        const int a_dir2,
                        const int a_dir3);

  /// Returns S in the given direction (codimension 5)
  const S& operator()(const int a_dir0,
                              const int a_dir1,
                              const int a_dir2,
                              const int a_dir3,
                              const int a_dir4) const;
  S& operator()(const int a_dir0,
                        const int a_dir1,
                        const int a_dir2,
                        const int a_dir3,
                        const int a_dir4);

  /// Returns S in the given direction (codimension 6)
  const S& operator()(const int a_dir0,
                              const int a_dir1,
                              const int a_dir2,
                              const int a_dir3,
                              const int a_dir4,
                              const int a_dir5) const;
  S& operator()(const int a_dir0,
                        const int a_dir1,
                        const int a_dir2,
                        const int a_dir3,
                        const int a_dir4,
                        const int a_dir5);

  /// Returns S from a sequential index
  const S& getSequential(const int a_iOrient) const;
  S&       getSequential(const int a_iOrient);

  /// Set all values in all FAB
  void setVal(const Real a_x);

  /// Orient the stored box in the othogonal directions of the given
  /// orientation
  Box orientBox(const int a_iOrient) const;

  /// Orient a centered box in the othogonal directions of the given
  /// orientation
  Box orientBox(const int a_iOrient, const Box& a_cbox) const;

  /// Orient a centered box in the othogonal directions of the given
  /// orientation.  Used externally for any codimension.
  static Box orientBox(const int a_codim,
                       const int a_iOrient,
                       const Box& a_cbox);

  /// Copy from another CodimBox
  void copy(const Box& a_R,
            const Interval& a_Cd,
            const CodimBox& a_src,
            const Interval& a_Cs);

  /// Copy from another CodimBox but to a different region
  void copy(const Box& a_Rs,        // Beware that the src and dest arguments
            const Interval& a_Cd,   // are not grouped together
            const Box& a_Rd,
            const CodimBox& a_src,
            const Interval& a_Cs);

  /// Report preallocation capability as non-static but symmetric
  static int preAllocatable();

  /// Returns size of serialized data
  int size(const Box& a_R, const Interval& a_comps) const;

  /// Writes a serialized representation of this CodimBox
  void linearOut(void* a_buf, const Box& a_R, const Interval& a_comps) const;

  /// Reads a serialized representation of this CodimBox
  void linearIn(void* a_buf, const Box& a_R, const Interval& a_comps);

protected:
  static int numOrient[numCD];        ///< Number of orientations for an object
                                      ///< of a given codimension
  static int totOrient[numCD];        ///< Total of number of orientations for
                                      ///< all smaller codimensions
  static unsigned bitOrient[CodimBox::numAI];
                                      ///< A sequential list of bit
                                      ///< representations of possible
                                      ///< orientations for each codimension
  static int indexFAB[CodimBox::numAI];
                                      ///< Returns the index of the S
                                      ///< for a given codimension and
                                      ///< orientation.  The input index
                                      ///< expresses the orientation as the sum
                                      ///< of the orthogonal directions.

private:

  Box m_box;                          ///< Cell-centered box over which this
                                      ///< CodimBox is defined
  int m_codim;                        ///< Codimension for which the Ses
                                      ///< are defined
  int m_nvar;                         ///< Number of components in the S
  Vector<S*> m_S;     ///< The Ses for each orientation

  // Copy and assignment not allowed
  CodimBox(const CodimBox&);
  CodimBox& operator=(const CodimBox&);

  /// Check that static arrays are initialized
  static bool isInitialized();

  /// Return the FArray box from a bit description of the orientation
  S&       getS(const unsigned i);
  const S& getS(const unsigned i) const;

  /// Get all the orthogonal directions of an orientation (general)
  static void genGetDirections(const int a_codim,
                               int a_bOrient,
                               int *const a_dir);

  /// Orient a centered box in the othogonal directions of the given
  /// orientation (general)
  static Box genOrientBox(int a_bOrient, Box a_cbox);

// Only of use with a register class
//   /// Shift all FABs by halves in one direction.
//   void shiftHalf(const int a_dir, const int a_numHalfs);
//  friend class LevelCodimRegister;

};

/*******************************************************************************
 *
 * Class CodimBox: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
///  Return the codimension
/*--------------------------------------------------------------------*/

template<typename S> int CodimBox<S>::getCodim() const
{
  return m_codim;
}

/*--------------------------------------------------------------------*/
///  Number of components
/*--------------------------------------------------------------------*/

template<typename S> int CodimBox<S>::nComp() const
{
  return m_nvar;
}

/*--------------------------------------------------------------------*/
///  Number of different orientations of codimensional objects
/*--------------------------------------------------------------------*/

template<typename S> int CodimBox<S>::getNumOrient() const
{
  return numOrient[m_codim];
}

/*--------------------------------------------------------------------*/
///  Number of different orientations of codimensional objects.  Used
///  externally for any codimension.
/*--------------------------------------------------------------------*/

template<typename S>int CodimBox<S>::getNumOrient(const int a_codim)
{
  CH_assert(CodimBox<S>::isInitialized());
  return numOrient[a_codim];
}

/*--------------------------------------------------------------------*/
///  Switch from a bit orientation to a sequential index of the
///  orientation
/*--------------------------------------------------------------------*/

template<typename S> int CodimBox<S>::bit2seq(const unsigned a_bOrient) const
{
  return indexFAB[a_bOrient];
}

/*--------------------------------------------------------------------*/
///  Switch from a sequential index to a bit representation of the
///  orientation
/*--------------------------------------------------------------------*/

template<typename S> unsigned CodimBox<S>::seq2bit(const int a_iOrient) const
{
  return bitOrient[totOrient[m_codim] + a_iOrient];
}

/*--------------------------------------------------------------------*/
///  Switch from a sequential index to a bit representation of the
///  orientation for any codimension
/*--------------------------------------------------------------------*/

template<typename S> unsigned CodimBox<S>::seq2bit(const int a_codim, const int a_iOrient)
{
  CH_assert(CodimBox<S>::isInitialized());
  return bitOrient[totOrient[a_codim] + a_iOrient];
}

/*--------------------------------------------------------------------*/
///  Get all the orthogonal directions of an orientation
/**  \param[in]  a_iOrient
 *                      The sequential orientation index from which to
 *                      get the directions
 *   \param[in]  a_dir  Array with minimum size [m_codim]
 *   \param[out] a_dir  Integers defining the directions
 *                      0 <= dir[i] < SpaceDim
 *   \note
 *   <ul>
 *     <li> The directions are orderd from smallest to largest.
 *   </ul>
 *//*-----------------------------------------------------------------*/

template<typename S> void CodimBox<S>::getDirections(const int a_iOrient, int *const a_dir) const
{
  genGetDirections(m_codim, seq2bit(a_iOrient), a_dir);
}

/*--------------------------------------------------------------------*/
///  Get all the orthogonal directions of an orientation.  Used
///  externally for any codimension.
/**  \param[in]  a_codim
 *                      The codimension
 *   \param[in]  a_iOrient
 *                      The sequential orientation index from which to
 *                      get the directions
 *   \param[in]  a_dir  Array with minimum size [m_codim]
 *   \param[out] a_dir  Integers defining the directions
 *                      0 <= dir[i] < SpaceDim
 *   \note
 *   <ul>
 *     <li> The directions are orderd from smallest to largest.
 *   </ul>
 *//*-----------------------------------------------------------------*/

template<typename S> void CodimBox<S>::getDirections(const int a_codim,
                                                     const int a_iOrient,
                                                     int *const a_dir)
{
  CH_assert(CodimBox<S>::isInitialized());
  genGetDirections(a_codim, seq2bit(a_codim, a_iOrient), a_dir);
}

/*--------------------------------------------------------------------*/
//  Returns cell-centered box which defines the CodimBox
/*--------------------------------------------------------------------*/
template<typename S> const Box& CodimBox<S>::box() const
{
  return m_box;
}

/*--------------------------------------------------------------------*/
/// Returns S in the directions defined by a Box
/*--------------------------------------------------------------------*/

template<typename S> S& CodimBox<S>::operator()(const Box &a_box)
{
  return this->operator()(a_box.ixType());
}

/*--------------------------------------------------------------------*/
///  Returns S from a sequential index
/**  \param[in]  a_iOrient
 *                      A sequential orientation index.
 *   \note
 *   <ul>
 *     <li> A bit description of the orientation can be obtained from
 *          \verbatim
 *            bitOrient[totOrient[a_codim] + a_iOrient]
 *          \endverbatim
 *   </ul>
 *//*-----------------------------------------------------------------*/
template<typename S> const S& CodimBox<S>::getSequential(const int a_iOrient) const
{
  CH_assert(a_iOrient < numOrient[m_codim]);
  return *m_S[a_iOrient];
}

template<typename S> S& CodimBox<S>::getSequential(const int a_iOrient)
{
  CH_assert(a_iOrient < numOrient[m_codim]);
  return *m_S[a_iOrient];
}

/*--------------------------------------------------------------------*/
/// Orient the stored box in the othogonal directions of the given
/// orientation
/**  \param[in]  a_bOrient
 *                      A bit orientation index.
 *//*-----------------------------------------------------------------*/

template<typename S> Box CodimBox<S>::orientBox(const int a_bOrient) const
{
  return orientBox(a_bOrient, m_box);
}

/*--------------------------------------------------------------------*/
///  Orient a centered box in the othogonal directions of the given
///  orientation
/**  \param[in]  a_iOrient
 *                      A sequential orientation index.  This is
 *                      converted to a bit orientation for finding the
 *                      orthogonal directions.
 *   \param[in]  a_cbox
 *                      A centered box.
 *   \return            The orientated box
 *//*-----------------------------------------------------------------*/

template<typename S> Box CodimBox<S>::orientBox(const int a_iOrient, const Box& a_cbox) const
{
  return genOrientBox(seq2bit(a_iOrient), a_cbox);
}

/*--------------------------------------------------------------------*/
///  Orient a centered box in the othogonal directions of the given
///  orientation.  Used externally for any codimension.
/**  \param[in]  a_codim
 *                      The codimension of the box
 *   \param[in]  a_iOrient
 *                      A sequential orientation index.  This is
 *                      converted to a bit orientation for finding the
 *                      orthogonal directions.
 *   \param[in]  a_cbox
 *                      A centered box.
 *   \return            The orientated box
 *//*-----------------------------------------------------------------*/

template<typename S> Box CodimBox<S>::orientBox(const int a_codim, const int a_iOrient, const Box& a_cbox)
{
  CH_assert(CodimBox<S>::isInitialized());
  return genOrientBox(seq2bit(a_codim, a_iOrient), a_cbox);
}

/*--------------------------------------------------------------------*/
///  Report preallocation capability as non-static but symmetric
/**  The size of the class is dependent upon the codimension.  See
 *   BoxLayoutData.H for more info
 *//*-----------------------------------------------------------------*/

template<typename S> int CodimBox<S>::preAllocatable()
{
  return 1;
}

/*--------------------------------------------------------------------*/
///  Check that static arrays are initialized
/**  Test assumes static arrays have at least been default
 *   initialized.
 *   \return            T - Arrays initialized
 *//*-----------------------------------------------------------------*/

template<typename S> bool CodimBox<S>::isInitialized()
{
  return (indexFAB[0] != -1);
}

/*--------------------------------------------------------------------*/
///  Return the FArray box from a bit description of the orientation
/**  The appropriate FArray box for this codimension and orientation
 *   is given by a query to 'indexFAB'
 *   \param[in]  a_bOrient
 *                      Bit description of the orientation
 *//*-----------------------------------------------------------------*/
template<typename S>  const S& CodimBox<S>::getS(const unsigned a_bOrient) const
{
  return *(m_S[bit2seq(a_bOrient)]);
}

template<typename S> S& CodimBox<S>::getS(const unsigned a_bOrient)
{
  return *(m_S[bit2seq(a_bOrient)]);
}

/*******************************************************************************
 */
///  Factory object to create CodimBox(s) for BoxLayoutData and similar
///  containers.
/*
 ******************************************************************************/
template <class S> class CodimBoxFactory : public DataFactory<CodimBox<S> >
{

public:
  /// Constructor stores codimension
  CodimBoxFactory(int a_codim)
    :
    m_codim(a_codim)
    { }

  /// Create a new CodimBox
  CodimBox<S>* create(const Box& a_box,
                      int a_ncomps,
                      const DataIndex& a_datInd) const
    {
      return new CodimBox<S>(m_codim, a_box, a_ncomps);
    }

private:
  int m_codim;                        ///< Codimension

};

#include "NamespaceFooter.H"

#include "CodimBoxImplem.H"

#endif
