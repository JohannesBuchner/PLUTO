#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BOX_H_
#define _BOX_H_

#ifdef NODEV
#undef NODEV
#endif

#include "SPACE.H"

#ifndef WRAPPER
#include <iostream>

#include "IntVect.H"
#include "Misc.H"
#include "LoHiSide.H"
#include "Vector.H"
#include "SPMD.H"
#include "NamespaceHeader.H"

struct SliceSpec;

/// Cell-Based or Node-Based Indices
/**

  The class IndexType defines an index as being cell based or node (face)
  based in each of the CH_SPACEDIM directions.  This class defines an
  enumerated type CellIndex to be either CELL or NODE; i.e. each of the
  CH_SPACEDIM dimensions must be either CELL or NODE.
*/

class IndexType
{
public:
  ///
  /**
     The cell index type: one of CELL or NODE.

  */
  enum CellIndex
  {
    CELL = 0,
    NODE = 1
  };

  ///
  /**
     The default constructor.

  */
  IndexType ();

  ///
  /**
      The copy constructor.

  */
  IndexType (const IndexType& rhs);

  ///
  /**
     Constructs an IndexType identical to an IntVect.

  */
  explicit IndexType (const IntVect& iv);

  ///
  /**
     The assignment operator.

  */
  IndexType& operator= (const IndexType& rhs);

  ///
  /**
     Another assignment operator.

  */
  IndexType& operator= (const IntVect& iv);

  ///
  /**
      Constructs an IndexType given an explicit CellIndex for each
      direction.  D_DECL6 is a macro that sets the constructor to take
      BL_SPACEDIM arguments.

  */
  IndexType (D_DECL6(CellIndex i, CellIndex j, CellIndex k,
                     CellIndex l, CellIndex m, CellIndex n));

  ///
  /**
     Sets this IndexType to be NODE based in direction dir.

  */
  void set (int dir);

  ///
  /**
     Sets this IndexType to be CELL based in direction dir.

  */
  void unset (int dir);

  ///
  /**
     True if this IndexType is NODE based in direction dir.

  */
  bool test (int dir) const;

  ///
  /**
     Set this IndexType to be NODE based in all directions.

  */
  void setall ();

  ///
  /**
     Set this IndexType to be CELL based in all directions.

  */
  void clear ();

  ///
  /**
     True if this IndexType is NODE based in any direction.

  */
  bool any () const;

  ///
  /**
     True if this IndexType is valid.

  */
  bool ok () const;

  ///
  /**
     Change this IndexType from CELL to NODE or NODE to CELL
     in direction dir.

  */
  void flip (int i);

  ///
  /**
     True if IndexTypes are identical.

  */
  bool operator== (const IndexType& t) const;

  ///
  /**
    IndexTypes comparison.
  */
  bool operator<(const IndexType& t) const;

  ///
  /**
     True if IndexTypes are not identical.

  */
  bool operator!= (const IndexType& t) const;

  ///
  /**
     True if this IndexType is CELL based in all directions.

  */
  bool cellCentered () const;

  ///
  /**
     True if this IndexType is NODE based in all directions.

  */
  bool nodeCentered () const;

  ///
  /**
     Set this IndexType to CellIndex type t in direction dir.

  */
  void setType (int       dir,
                CellIndex t);

  ///
  /**
     Returns  the CellIndex of this IndexType in direction dir.

  */
  CellIndex ixType (int dir) const;

  ///
  /**
     Return an integer representing this IndexType in direction dir.

  */
  int operator[] (int dir) const;

  ///
  /**
     Convert this IndexType to an IntVect via component-wise
     conversion of enum CellIndex to int.

  */
  IntVect ixType () const;

  ///
  /**
      This static member function returns an IndexType object all of
      whose components are of value IndexType::CELL.  It is provided as
      a convenience to our users when defining a Box all of whose faces
      should be of type IndexType::CELL.

  */
  static IndexType TheCellType ();

  ///
  /**
     This static member function returns an IndexType object all of
     whose components are of value IndexType::NODE.  It is provided as
     a convenience to our users when defining a Box all of whose faces
     should be of type IndexType::NODE.

  */
  static IndexType TheNodeType ();

  ///
  /**
     Write an IndexType to an ostream in ASCII.

  */
  friend std::ostream& operator<< (std::ostream&         os,
                                   const IndexType& itype);

  ///
  /**
     Read an IndexType from an istream.

  */
  friend std::istream& operator>> (std::istream&   is,
                                   IndexType& itype);
private:
  //
  // Returns 1<<k.
  //
  //  static int mask (int k);
  static unsigned char mask (int k);
  //
  // An integer holding the CellIndex in bits 0 - BL_SPACEDIM-1.
  //
  //unsigned int itype;
  unsigned char itype;
};

//
// Inlines.
//

inline
unsigned char
IndexType::mask (int k)
{
  return 1<<k;
}

inline
IndexType::IndexType ()
  :
  itype(0)
{
}

inline
IndexType::IndexType (const IndexType& bt)
  :
  itype(bt.itype)
{
}

inline
IndexType& IndexType::operator= (const IndexType& bt)
{
  itype = bt.itype;
  return *this;
}

inline
IndexType::IndexType (const IntVect& iv)
{
  itype = D_TERM6((iv[0]?1:0), | ((iv[1]?1:0)<<1), | ((iv[2]?1:0)<<2),
                  | (iv[3]?1:0), | ((iv[4]?1:0)<<1), | ((iv[5]?1:0)<<2));
}

inline
IndexType& IndexType::operator= (const IntVect& iv)
{
  itype = D_TERM6((iv[0]?1:0), | ((iv[1]?1:0)<<1), | ((iv[2]?1:0)<<2),
                  | (iv[3]?1:0), | ((iv[4]?1:0)<<1), | ((iv[5]?1:0)<<2));
  return *this;
}

inline
IndexType::IndexType (D_DECL6(CellIndex i, CellIndex j, CellIndex k,
                              CellIndex l, CellIndex m, CellIndex n))
{
  itype = D_TERM6(i, | (j<<1), | (k<<2),
           | (l<<3), | (m<<4), | (n<<5));
}

inline
void
IndexType::set (int dir)
{
  itype |= mask(dir);
}

inline
void
IndexType::unset (int dir)
{
  itype &= ~mask(dir);
}

inline
bool
IndexType::test (int dir) const
{
  return (itype & mask(dir)) != 0;
}

inline
void
IndexType::setall ()
{
  itype = (1 << CH_SPACEDIM) - 1;
}

inline
void
IndexType::clear ()
{
  itype = 0;
}

inline
bool
IndexType::any () const
{
  return itype != 0;
}

inline
bool
IndexType::ok () const
{
  return itype < (1 << CH_SPACEDIM);
}

inline
void
IndexType::flip (int i)
{
  itype ^= mask(i);
}

inline
bool
IndexType::operator== (const IndexType& t) const
{
  return t.itype == itype;
}

inline
bool
IndexType::operator< (const IndexType& t) const
{
  return itype < t.itype;
}

inline
bool
IndexType::operator!= (const IndexType& t) const
{
  return t.itype != itype;
}

inline
bool
IndexType::cellCentered () const
{
  return itype == 0;
}

inline
bool
IndexType::nodeCentered () const
{
  return itype == (1<<CH_SPACEDIM)-1;
}

inline
void
IndexType::setType (int       dir,
                    CellIndex t)
{
  t == CELL ? unset(dir) : set(dir);
}

inline
IndexType::CellIndex
IndexType::ixType (int dir) const
{
  return (CellIndex) ((itype & (1<<dir)) >> dir);
}

inline
int
IndexType::operator[] (int dir) const
{
  return test(dir);
}

inline
IntVect
IndexType::ixType () const
{
  return IntVect(D_DECL6(itype&1, (itype>>1)&1, (itype>>2)&1,
                         (itype>>3)&1, (itype>>4)&1, (itype>>5)&1));
}

#endif /* WRAPPER */

//
/// A Rectangular Domain on an Integer Lattice
/**

  A Box is an abstraction for defining discrete rectangular regions of
  SpaceDim-dimensioned indexing space.  Boxes have an IndexType, which
  defines IndexType::CELL or IndexType::NODE based points for each
  direction and a low and high IntVect which defines the lower and
  upper corners of the Box.  Boxes can exist in positive and negative
  indexing space.

  Box is a dimension dependent class, so SpaceDim must be
  defined as either 1, 2, or 3 when compiling.

  @ingroup chombo
*/

class Box
{
public:

  ///

  ///
  /**
     The default constructor.  The constructed Box is empty.

  */
  Box ();

  ///
  /**
     The destructor.

  */
  ~Box ()
  {}

  ///
  /**
     Construct cell-centered type Box.  It is an error if small is
     greater than big.

  */
  Box (const IntVect& small,
       const IntVect& big);

  ///
  /**
     Modify cell-centered type Box.  It is an error if small is
     greater than big.

  */
  void define(const IntVect& small, const IntVect& big);

  ///
  /**
     Construct Box with specified lengths.  It is an error
     if the lengths are negative.

  */
  Box (const IntVect& small,
       const int*     vec_len);

  ///
  /**
      Construct Box with given type.  small and big are expected to be
      consistent with given type.  It is an error if small is greater
      than big.  It is an error if typ is not convertible to an
      IndexType.

  */
  Box (const IntVect& small,
       const IntVect& big,
       const IntVect& typ);

  ///
  /**
      Modify Box with given type.  small and big are expected to be
      consistent with given type.  It is an error if small is greater
      than big.  It is an error if typ is not convertible to an
      IndexType.

  */
  void define(const IntVect& small,
              const IntVect& big,
              const IntVect& typ);

  ///
  /**
      Construct Box with given type.  small and big are expected to be
      consistent with given type.  It is an error if small is greater
      than big.

  */
  Box (const IntVect&   small,
       const IntVect&   big,
       const IndexType& t);

  ///
  /**
      Modify Box with given type.  small and big are expected to be
      consistent with given type.  It is an error if small is greater
      than big.

  */
  void define(const IntVect&   small,
              const IntVect&   big,
              const IndexType& t);

  ///
  /**
     The copy constructor.

  */
  Box (const Box& b);

  void define(const Box& b);

  Box copy() const
  {
    return *this;
  }

  /// {\bf Accessors}

  ///
  /**
      Returns the lower corner of this Box.

  */
  const IntVect& smallEnd () const;

  ///
  /**
     return smallend if side is low,
     bigend if side is high.
  */
  IntVect sideEnd(Side::LoHiSide a_side) const;

  ///
  /**
     Returns the coordinate of the low end of this Box in the given
     direction.  Directions are zero-based.  It is an error if not 0
     <= dir < SpaceDim.

  */
  int smallEnd (int dir) const;

  ///
  /**
     Returns the upper corner of this Box.

  */
  const IntVect& bigEnd () const;

  ///
  /**
     Returns the coordinate of the high end of this Box in the given
     direction.  Directions are zero-based.  It is an error if not 0
     <= dir < SpaceDim.

  */
  int bigEnd (int dir) const;

  ///
  /**
     Returns a constant pointer to the array of low end coordinates.
     Useful for calls to Fortran.  It should not be used in any other
     circumstances.

  */
  const int* loVect () const;

  ///
  /**
     Returns a constant pointer the array of high end coordinates.
     Useful for calls to Fortran.  It should not be used in any other
     circumstances.

  */
  const int* hiVect () const;

  ///
  /**
     Returns a constant pointer to the array of coordinates in the
     Box.  Useful for calls to Fortran, but otherwise too dangerous
     for use.

  */
  const int* getVect () const;

  ///
  /**
     Returns offset of point from smallend; i.e.  index(smallend) ->
     0, bigend would return volume()-1.  Is used in accessing
     FArrayBox.

  */
  long index (const IntVect& v) const;

  /// centering

  ///
  /**
     Return the indexing type of this Box.

  */
  IndexType ixType () const;

  ///
  /**
     Return the indexing type of this Box.

  */
  IntVect type () const;

  ///
  /**
     Return the indexing type of this Box in the specified direction.
     Directions are zero-based.  It is an error if not 0 <= dir <
     SpaceDim.

  */
  IndexType::CellIndex type (int dir) const;

  /// size functions

  ///
  /**
     Return an IntVect containing the size of this Box in each
     coordinate direction.

  */
  inline
  IntVect size () const;

  ///
  /**
     Return the size of this Box in the specified coordinate
     direction.  Directions are zero-based.  It is an error
     if not 0 <= dir < SpaceDim.

  */
  int size (int dir) const;

  ///
  /**
     Is the number calculated by numPts() representable in a long?

  */
  bool numPtsOK () const;

  ///
  /**
     Return the volume, in indexing space, of region enclosed by this
     Box.  This is identical to volume().  Abort()s if the number
     cannot be represented in a long.\\

  */
  long numPts () const;

  ///
  /**
     Is the number calculated by volume() representable in a long?

  */
  bool volumeOK () const;

  ///
  /**
     Return the volume, in indexing space, of region enclosed by this
     Box.  This is identical to numPts().  Abort()s if the number
     cannot be represented in a long.\\

  */
  long volume () const;

  ///
  /**
     Returns length of the longest side of this Box.  The argument dir
     is modified to give direction with longest side:
     0...SpaceDim-1. Ignores type.

  */
  int longside (int& dir) const;

  ///
  /**
     Returns length of the longest side of this Box.  Ignores type.

  */
  int longside () const;

  ///
  /**
     Returns length of shortest side of this Box.  The argument dir is
     modified to give direction with shortest side: 0...SpaceDim-1.
     Ignores type.

  */
  int shortside (int& dir) const;

  ///
  /**
     Returns length of the shortest side of this Box.  Ignores type.

  */
  int shortside () const;

  /// {\bf Comparison Functions}

  ///
  /**
     Returns true if this Box is empty.

  */
  bool isEmpty () const;

  ///
  /**
     Returns true if this Box is not empty and all
     sizes are right
  */
  bool ok() const;

  ///
  /**
     Returns true if argument is contained within this Box.  An empty
     Box does not contain and is not contained by any Box, including
     itself.

  */
  bool contains (const IntVect& p) const;

  ///
  /**
     Returns true if argument is contained within this Box.  It is an
     error if the Boxes have different types.  An empty Box does not
     contain any IntVect.

  */
  bool contains (const Box& b) const;

  ///
  /**
     Returns true if this Box and the argument have non-null
     intersections.  It is an error if the Boxes have different types.
     An empty Box does not intersect any Box, including itself.

  */
  bool intersects (const Box& b) const;

  ///
  /**
     Returns true if this Box and the argument have non-null
     intersections.  It is an error if the Boxes have different types.
     This routine does not perform the check to see if *this or b are
     empty boxes.  It is the callers responsibility to ensure that
     this never happens.  If you are unsure, the use the .intersects(..) routine.

  */
  bool intersectsNotEmpty (const Box& b) const;

  ///
  /**
     Returns true if this Box and the argument are the same size, ie
     translates of each other.  It is an error if they have different
     types.
  */
  bool sameSize (const Box& b) const;

  ///
  /**
     Returns true if this Box and the argument have same type.

  */
  bool sameType (const Box &b) const;

  ///
  /**
     Returns true if this Box and the argument are identical,
     including type.

  */
  bool operator== (const Box& b) const;

  bool eq(const Box& b) const;
  ///
  /**
     Returns true if this Box and the argument differ, including type.

  */
  bool operator!= (const Box& b) const;

  bool neq(const Box& b) const;
  ///
  /**
     Returns true if this Box is cell-centered in all indexing
     directions.

  */
  bool cellCentered () const;

  // following operators added Sept. 14, 1999.  bvs
  ///
  /**
     Returns true if this Box is lexigraphically less than rhs box.
     All comparison is based on lower box corner.  In the name of
     coding efficiency, we do not handle the case where two boxes have
     the same lower left corner. In DEBUG mode, this is checked for
     with an assert.

  */
  bool operator < (const Box& rhs) const;

  bool lt(const Box& rhs) const;

  /// {\bf Modification Functions}

  ///
  /**
     The assignment operator.

  */
  Box& operator= (const Box& b);

  ///
  /**
     Redefine the lower corner of this Box.  It is an error
     if the specified corner is greater than the exisiting upper
     corner of this Box.

  */
  Box& setSmall (const IntVect& sm);

  ///
  /**
     Redefines the lower end of this Box in the specified coordinate
     direction.  It is an error if the specified value is greater than
     the exisiting upper end of this Box.  Directions are zero-based.
     It is an error if not 0 <= dir < SpaceDim.

  */
  Box& setSmall (int dir,
                 int sm_index);

  ///
  /**
     Redefines the upper corner of this Box.  It is an error if the
     specified corner is less than the exisiting lower corner of this
     Box.

  */
  Box& setBig (const IntVect& bg);

  ///
  /**
     Redefines the upper end of this Box in the specified coordinate
     direction.  It is an error if the specified value is less than the
     exisiting lower end of this Box.  Directions are zero-based.  It
     is an error if not 0 <= dir < SpaceDim.

  */
  Box& setBig (int dir,
               int bg_index);

  ///
  /**
     Set the entire range in a given direction, starting at sm_index
     with length n_cells.  It is an error if n_cells <= 0.

  */
  Box& setRange (int dir,
                 int sm_index,
                 int n_cells = 1);

  /// centering type conversion functions

  ///
  /**
     Modifies this Box by converting from the current type into the
     argument type.  This may change the Box coordinates:\\

     type CELL -> NODE : increase coordinate by one on high end.\\

     type NODE -> CELL : reduce coordinate by one on high end.\\

     other type mappings make no change.

  */
  Box& convert (IndexType typ);

  ///
  /**
     Modifies this Box by converting from the current type into the
     argument type.  This may change the Box coordinates:\\ type CELL
     -> NODE : increase coordinate by one on high end.\\ type NODE ->
     CELL : reduce coordinate by one on high end.\\ other type
     mappings make no change.  It is an error if typ is not
     convertible to an IndexType.

  */
  Box& convert (const IntVect& typ);

  ///
  /**
     Modifies this Box by converting from the current type into the
     argument type.  This may change the Box coordinates:\\

     type CELL -> NODE : increase coordinate by one on high end.\\

     type NODE -> CELL : reduce coordinate by one on high end.\\

     Other type mappings make no change.  Directions are zero-based.
     It is an error if not 0 <= dir < SpaceDim.

  */
  Box& convert (int                  dir,
                IndexType::CellIndex typ);

  ///
  /**
     Modifies this Box by converting it to NODE type in all
     directions.  This increases all components of the upper corner by
     one.  The Empty Box is not modified by this function except for
     changing of type.

  */
  Box& surroundingNodes ();

  ///
  /**
     Modifies this Box by converting it to NODE type in given
     direction.  This increases the component of the upper corner in
     the given direction by one.  Directions are zero-based.  It is an
     error if not 0 <= dir < SpaceDim.  The Empty Box is not modified
     by this function except for changing of type.

  */
  Box& surroundingNodes (int dir);

  Box& surroundingNodes_int(int dir);
  ///
  /**
     Returns a Box with NODE based coordinates in direction dir that
     encloses the argument Box.  NOTE: equivalent to
     b.convert(dir,NODE).  NOTE: error if b.type(dir) == NODE.
     Directions are zero-based.  It is an error if not 0 <= dir <
     SpaceDim.  The Empty Box is not modified by this function except
     for changing of type.

  */
  friend  Box surroundingNodes (const Box& b,
                                int        dir);

  ///
  /**
     Returns a Box with NODE based coordinates in all directions that
     encloses argument Box.  The upper corner of the return Box will
     be one more in all components than the upper corner of the
     argument Box.  The Empty Box is not modified by this function
     except for changing of type.

  */
  friend Box surroundingNodes (const Box& b);

  ///
  /**
     Modifies this Box by converting it to CELL type in all directions.
     This decreases all components of the upper corner by one.  The
     Empty Box is not modified by this function except for changing of
     type.

  */
  Box& enclosedCells ();

  ///
  /**
     Modifies this Box by converting it to CELL type in given
     direction.  This decreases the component of the upper corner in
     the given direction by one.  Directions are zero-based.  It is an
     error if not 0 <= dir < SpaceDim.  The Empty Box is not modified
     by this function except for changing of type.

  */
  Box& enclosedCells (int dir);

  Box& enclosedCells_int (int dir);

  ///
  /**
     Returns a Box with CELL based coordinates in direction dir that
     is enclosed by argument Box.  NOTE: equivalent to
     b.convert(dir,CELL).  NOTE: error if b.type(dir) == CELL.
     Directions are zero-based.  It is an error if not 0 <= dir <
     SpaceDim.  The Empty Box is not modified by this function except
     for changing of type.

  */
  friend Box enclosedCells (const Box& b,
                            int        dir);

  ///
  /**
     Returns a Box with CELL based coordinates in all directions that
     is enclosed by argument Box.  The upper corner of the return Box
     will be one less in all components than the upper corner of the
     argument Box. The Empty Box is not modified by this function
     except for changing of type.

  */
  friend Box enclosedCells (const Box& b);

  /// shift functions

  ///
  /**
     Modifies this Box by shifting it nzones indexing positions in
     coordinate direction dir.  Directions are zero-based.  It is an
     error if not 0 <= dir < SpaceDim.  The Empty Box is not modified
     by this function.

  */
  Box& shift (int dir,
              int nzones);

  ///
  /**
     Modifies this Box by shifting.  Equivalent to
     b.shift(0,iv[0]).shift(1,iv[1]) ....  The Empty Box is not
     modified by this function.

  */
  Box& shift (const IntVect& iv);

  Box& shift_intvect (const IntVect& iv);

  ///
  /**
     This member modifies this Box by shifting by "half" indices,
     thereby converting the Box from type CELL to NODE or vice-versa.
     b.shiftHalf(0,1) shifts b to the right by 1/2 cells.
     b.shiftHalf(1,-3) shifts b in the -j direction by 3/2 cells.
     NOTE: If num_halfs is EVEN the shift is num_halfs/2 full zones
     and hence will not change the type.  This is: b.shifthalf(4) ==
     b.shift(2).  Directions are zero-based.  It is an error if not 0
     <= dir < SpaceDim.  The Empty Box is not modified by this
     function except for changing of type.

  */
  Box& shiftHalf (int dir,
                  int num_halfs);

  ///multiblock stuff.
  void convertOldToNew(const IntVect& a_permutation,
                       const IntVect& a_sign,
                       const IntVect& a_translation);

  ///multiblock stuff
  void convertNewToOld(const IntVect& a_permutation,
                       const IntVect& a_sign,
                       const IntVect& a_translation);

  ///
  /**
     Modifies this Box by shifting by half indices.  Equivalent to
     b.shiftHalf(0,iv[0]).shiftHalf(1,iv[1]) ...  The Empty Box is not
     modified by this function except for changing of type.

  */
  Box& shiftHalf (const IntVect& iv);
  Box& shiftHalf_intvect (const IntVect& iv);

  ///
  /**
     Modifies this Box by shifting it by given IntVect.  The Empty Box
     is not modified by this function.

  */
  Box& operator+= (const IntVect& v);

  ///
  /**
     Returns a Box that is this Box shifted by the given IntVect.  The
     Empty Box is not modified by this function.

  */
  Box  operator+  (const IntVect& v) const;

  ///
  /**
     Modifies this Box by shifting it by given IntVect.  The Empty Box
     is not modified by this function.

  */
  Box& operator-= (const IntVect& v);

  ///
  /**
     Returns a Box that is this Box shifted by the given IntVect.  The
     Empty Box is not modified by this function.

  */
  Box  operator-  (const IntVect& v) const;

  /// neighbor box functions

  ///
  friend Box bdryBox(const Box& b,
                     int        dir,
                     Side::LoHiSide a_sd,
                     int        len);

  ///
  /**
     Returns the edge-centered Box (in direction dir) defining the low
     side of the argument Box.  The output Box will have the given
     length in the given direction.  Directions are zero-based.  It is
     an error if not 0 <= dir < SpaceDim.  The neighbor of an Empty Box
     is an Empty Box of the appropriate type.

  */
  friend Box bdryLo (const Box& b,
                     int        dir,
                     int        len);

  ///
  /**
     Returns the edge-centered Box (in direction dir) defining the
     high side of the argument Box.  The return Box will have the
     given length in the given direction.  Directions are zero-based.
     It is an error if not 0 <= dir < SpaceDim.  The neighbor of an
     Empty Box is an Empty Box of the appropriate type.

  */
  friend Box bdryHi (const Box& b,
                     int        dir,
                     int        len);

  ///
  /**
     Returns the cell centered Box of the given length adjacent to the
     argument Box on the low end along the given coordinate direction.
     The return Box is identical to the argument Box in the other
     directions.  The return Box and the argument Box have an empty
     intersection.  \\

     NOTE: len != 0. A negative length results in cells _inside_ the original box. \\

     NOTE: Box retval = adjCellLo(b,dir,len) is equivalent to the
     following set of operations: \\

     Box retval(b); \\

     retval.convert(dir,Box::CELL); \\

     retval.setrange(dir,retval.smallEnd(dir)-len,len); \\

     Directions are zero-based.  It is an error if not 0 <= dir <
     SpaceDim.  The neighbor of an Empty Box is an Empty Box of the
     appropriate type.

  */
  friend Box adjCellLo (const Box& b,
                        int        dir,
                        int        len);

  ///
  /**
     Returns the cell centered Box of the given length adjacent to the
     argument Box on the high end along the given coordinate
     direction.  The return Box is identical to the argument Box in
     the other directions.  The return Box and the argument Box have
     an empty intersection.  \\

     NOTE: len != 0. A Negative length results in cells _inside_ the original box. \\

     NOTE: Box retval = adjCellHi(b,dir,len) is equivalent to the
     following set of operations: \\

     Box retval(b); \\

     retval.convert(dir,Box::CELL); \\

     retval.setrange(dir,retval.bigEnd(dir)+1,len);\\

     Directions are zero-based.  It is an error if not 0 <= dir <
     SpaceDim.  The neighbor of an Empty Box is an Empty Box of the
     appropriate type.

  */
  friend Box adjCellHi (const Box& b,
                        int        dir,
                        int        len);

  /// intersection functions

  ///
  /**
     Returns the Box that is intersection of this Box and the argument
     Box.  Intersection is commutative.  The Boxes MUST be of same
     type.  The intersection of the Empty Box and any Box is the Empty
     Box.

  */
  Box operator& (const Box&) const;

  ///
  /**
     Modifies this Box by intersection with the argument Box.  The
     Boxes MUST be of the same type.  The intersection of the Empty
     Box and any Box is the Empty Box.

  */
  Box& operator&= (const Box&);

  ///
  friend Box adjCellBox (const Box& b,
                         int        dir,
                         Side::LoHiSide a_side,
                         int        len);
  ///
  /**
     Modifies this Box to that of the minimum Box containing both this
     Box and the argument Box.  Both Boxes must have identical type.

  */
  Box& minBox (const Box& b);

  ///
  /**
     Returns a Box that is the minimum Box containing both the
     argument Boxes.  Both Boxes must have identical type.

  */
  friend Box minBox (const Box& b1,
                     const Box& b2);

  /// grow functions

  ///
  /**
     Modifies this Box by growing it in all directions by given
     amount.  The Empty Box is not modified by this function.  NOTE:
     n_cell negative shrinks the Box by that number of cells.  If
     shrinking produces an empty Box, the Box is transformed into the
     canonical Empty Box.

  */
  Box& grow (int i);

  ///
  /**
     Returns a Box that is the argument Box grown in all directions by
     given amount.  The Empty Box is not modified by this function.
     NOTE: n_cell negative shrinks the Box by that number of cells.
     If shrinking produces an empty Box, the Box is transformed into
     the canonical Empty Box.

  */
  friend Box grow (const Box& b,
                   int        i);

  ///
  /**
     Modifies this Box by growing it in each direction by specified
     amount.  The Empty Box is not modified by this function.  NOTE:
     components of iv may be negative, which would shrink this Box.  If
     shrinking produces an empty Box, the Box is transformed into the
     canonical Empty Box.

  */
  Box& grow (const IntVect& v);

  ///
  /**
     Returns a Box that is the argument Box grown in each
     direction by specified amount.  The Empty Box is not modified by
     this function.  NOTE: components of iv may be negative, which
     would return a shrunken Box.  If shrinking produces an empty
     Box, the Box is transformed into the canonical Empty Box.

  */
  friend  Box grow (const Box&     b,
                    const IntVect& v);

  ///
  /**
     Modifies this Box by growing it on the low and high end by n_cell
     cells in direction idir.  The Empty Box is not modified by this
     function.  NOTE: n_cell negative shrinks this Box by that number
     of cells.  If shrinking produces an empty Box, the Box is
     transformed into the canonical Empty Box.  Directions are
     zero-based.  It is an error if not 0 <= dir < SpaceDim.

  */
  Box& grow (int idir,
             int n_cell);

  ///
  /**
     Modifies this Box by growing it on the low end by n_cell cells
     in direction idir.  The Empty Box is not modified by this
     function.  NOTE: n_cell negative shrinks this Box by that number
     of cells.  If shrinking produces an empty Box, the Box is
     transformed into the canonical Empty Box.  Directions are
     zero-based.  It is an error if not 0 <= dir < SpaceDim.

  */
  Box& growLo (int idir,
               int n_cell=1);

  ///
  /**
     growLo if a_sd == Side::Lo \\
     growHi if a_sd= = Side::Hi
   */
  Box& growDir (int a_idir,
                const Side::LoHiSide& a_sd,
                int a_cell);

  ///
  /**
     Modifies this Box by growing it on the high end by n_cell cells
     in direction idir.  The Empty Box is not modified by this
     function.  NOTE: n_cell negative shrinks the Box by that number
     of cells.  If shrinking produces an empty Box, the Box is
     transformed into the canonical Empty Box.  Directions are
     zero-based.  It is an error if not 0 <= dir < SpaceDim.

  */
  Box& growHi (int idir,
               int n_cell=1);

  /// refinement

  ///
  /**
     Modifies this Box by refining it by given (positive) refinement
     ratio.  The Empty Box is not modified by this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo*ratio and hi <-
     (hi+1)*ratio - 1. \\

     NOTE: if type(dir) = NODE centered: lo <- lo*ratio and hi <-
     hi*ratio.

  */
  Box& refine (int refinement_ratio);

  ///
  /**
     Returns a Box that is the argument Box refined by given
     (positive) refinement ratio.  The Empty Box is not modified by
     this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo*ratio and hi <-
     (hi+1)*ratio - 1.\\

     NOTE: if type(dir) = NODE centered: lo <- lo*ratio and hi <-
     hi*ratio.

  */
  friend Box refine (const Box& b,
                     int        refinement_ratio);

  ///
  /**
     Modifies this Box by refining it by given (positive) refinement
     ratio.  The Empty Box is not modified by this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo*ratio and hi <-
     (hi+1)*ratio - 1.\\

     NOTE: if type(dir) = NODE centered: lo <- lo*ratio and hi <-
     hi*ratio.

  */
  Box& refine (const IntVect& refinement_ratio);

  ///
  /**
     Returns a Box that is the argument Box refined by given
     (positive) refinement ratio.  The Empty Box is not modified by
     this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo*ratio and hi <-
     (hi+1)*ratio - 1.\\

     NOTE: if type(dir) = NODE centered: lo <- lo*ratio and hi <-
     hi*ratio.

  */
  friend Box refine (const Box&     b,
                     const IntVect& refinement_ratio);

  /// coarsening

  ///
  /**
     Modifies this Box by coarsening it by given (positive) refinement
     ratio.  The Empty Box is not modified by this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo/ratio and hi <-
     hi/ratio.\\

     NOTE: if type(dir) = NODE centered: lo <- lo/ratio and hi <-
     hi/ratio + ((hi%ratio)==0 ? 0 : 1).\\

     That is, refinement of coarsened Box must contain the original
     Box.

  */

  Box& coarsen (int refinement_ratio);

  ///
  /**
     True if the if (refine(coarsen(int))) gets the same box back
  */
  bool coarsenable(int refrat) const
  {
    Box test = *this;
    test.coarsen(refrat);
    test.refine(refrat);
    return (test==*this);
  }
  ///
  /**
     Returns a Box that is the argument Box coarsened by given
     (positive) refinement ratio.  The Empty Box is not modified by
     this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo/ratio and hi <-
     hi/ratio.\\

     NOTE: if type(dir) = NODE centered: lo <- lo/ratio and hi <-
     hi/ratio + ((hi%ratio)==0 ? 0 : 1).\\

     That is, refinement of coarsened Box must contain the original
     Box.

  */
  friend Box coarsen (const Box& b,
                      int        refinement_ratio);

  ///
  /**
     Modifies this Box by coarsening by given (positive) refinement
     ratio.  The Empty Box is not modified by this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo/ratio and hi <-
     hi/ratio.\\

     NOTE: if type(dir) = NODE centered: lo <- lo/ratio and hi <-
     hi/ratio + ((hi%ratio)==0 ? 0 : 1).\\

     That is, refinement of coarsened Box must contain the original
     Box.

  */
  Box& coarsen (const IntVect& refinement_ratio);

  ///
  /**
     Returns a Box that is the argument Box coarsened by given
     (positive) refinement ratio.  The Empty Box is not modified by
     this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo/ratio and hi <-
     hi/ratio.\\

     NOTE: if type(dir) = NODE centered: lo <- lo/ratio and hi <-
     hi/ratio + ((hi%ratio)==0 ? 0 : 1).\\

     That is, refinement of coarsened Box must contain the original
     Box.

  */
  friend Box coarsen (const Box&     b,
                      const IntVect& refinement_ratio);

  // next(...) is out of favor.  use BoxIterator.
  /*
    Step through the rectangle.  It is a runtime error to give
    a point not inside rectangle.  Iteration may not be efficient.
  */
  void next (IntVect &) const;

  /*
    Scan argument IntVect over object second arg is
    increment vector.  Runtime error if IntVect is not
    contained in object Box.  Iteration may not be efficient.
  */
  void next (IntVect&   p,
             const int* shv) const;

  /// misc

  ///
  /**
     Chops the Box at the chop_pnt in the dir direction returns one
     Box, modifies the object Box.  The union of the two is the
     original Box.  The modified Box is the low end, the returned Box
     is the high end.  If type(dir) = CELL, the Boxes are disjoint with
     the chop_pnt included in the high end (new Box).  It is an ERROR
     if chop_pnt is the low end of the orig Box.  If type(dir) = NODE,
     the chop_pnt is included in both Boxes but is the only point in
     common.  It is also an error if the chop_pnt is an end node of
     the Box.  Directions are zero-based.  It is an error if not 0 <=
     dir < SpaceDim.

  */
  Box chop (int dir,
            int chop_pnt);

  /**
    Makes a_to a box that's the same as *this except that it's just one cell
    thick in the a_sliceSpec.direction-th direction, and its (one) coordinate
    in that direction is a_sliceSpec.position.

    If a_sliceSpec.position is outside the range of *this, and a_outofbounds
    isn't null, then we set *a_outofbounds to true, otherwise to false.

    It would be more natural to return a Box, than to set the value of a Box
    reference, but we do it this way to be consistent with BaseFab::degenerate()
    where we can't return an object because of the policy of disallowing copy
    constructors for heavyweight classes.
  */
  void degenerate( Box& a_to, const SliceSpec& a_sliceSpec,
                   bool* a_outofbounds=0 ) const;

  ///{\bf I/O Functions}

  ///
  /**
     Write an ASCII representation to the ostream.

  */
  friend std::ostream& operator<< (std::ostream&   os,
                                   const Box& bx);

  ///
  /**
     Read from istream.

  */
  friend std::istream& operator>> (std::istream& os,
                                   Box&     bx);

  ///
  /**
     print to pout()
  */
  void p() const;

  ///
  /**
     Gives more detail than printOn.  Useful for exiting due to an
     error.

  */
  void dumpOn (std::ostream& strm) const;

  /// {\bf Box Constants}

  ///
  /**
     This is a canonical empty Box of cell-centered type.

  */
  //static const Box Empty;

  // TheUnitBox is out of favor.
  /*
    This static member function returns a constant reference to
    an object of type Box representing the unit box in
    BL_SPACEDIM-dimensional space.
  */

  //static const Box& TheUnitBox ();

  //
  // Sets the 'len' element of the Box.  Aborts on integer overflow.
  //
  void computeBoxLen ();
  void computeBoxLenNotEmpty();

  /** Doesn't print the IndexType, and uses square outer brackets. */
  static void setTempestOutputFormat( bool );

protected:
  friend class HDF5Handle;
  //
  // A helper function for numPtsOK() and numPts().
  //
  bool numPtsOK (long& N) const;
  //
  // A helper function for volumeOK() and volume().
  //
  bool volumeOK (long& N) const;

  IntVect   smallend;
  IntVect   bigend;
  //IntVect   len;
  IndexType btype;

  static bool s_tempestOutputFormat;
};

#include "NamespaceFooter.H"

#include "BaseNamespaceHeader.H"

#include "NamespaceVar.H"

//Box specialization of linearSize
template < >
int linearSize(const CH_XDIR::Box& a_input);

//Box specialization of linearIn
template < >
void linearIn(CH_XDIR::Box& a_outputT, const void* const a_inBuf);

//Box specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const CH_XDIR::Box& a_inputT);

//Vector<Box>  specialization
template < >
int linearSize(const Vector<CH_XDIR::Box>& a_input);
template < >
void linearIn(Vector<CH_XDIR::Box>& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<CH_XDIR::Box>& a_inputT);

//Vector<Vector<Box> > specialization
template < >
int linearSize(const Vector<Vector<CH_XDIR::Box> >& a_input);
template < >
void linearIn(Vector<Vector<CH_XDIR::Box> >& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<Vector<CH_XDIR::Box> >& a_inputT);

#include "BaseNamespaceFooter.H"
#include "NamespaceHeader.H"

// global function prototypes
Box surroundingNodes (const Box& b,
                      int        dir);
Box surroundingNodes (const Box& b);
Box enclosedCells (const Box& b,
                   int        dir);
Box enclosedCells (const Box& b);
Box bdryBox(const Box& b,
            int        dir,
            Side::LoHiSide a_sd,
            int        len=1);
Box bdryLo (const Box& b,
            int        dir,
            int        len=1);
Box bdryHi (const Box& b,
            int        dir,
            int        len=1);
Box adjCellLo (const Box& b,
               int        dir,
               int        len=1);
Box adjCellHi (const Box& b,
               int        dir,
               int        len=1);
Box minBox (const Box& b1,
            const Box& b2);
Box coarsen (const Box&     b,
             const IntVect& refinement_ratio);
Box coarsen (const Box&     b,
             int refinement_ratio);
Box refine  (const Box&     b,
             const IntVect& refinement_ratio);
Box refine  (const Box&     b,
             int refinement_ratio);

//
// Inlines.
//

#ifndef WRAPPER

inline
Box::Box (const Box& b)
  : smallend(b.smallend),
    bigend(b.bigend),
    btype(b.btype)
{
//   D_EXPR6(len[0] = b.len[0],
//          len[1] = b.len[1],
//          len[2] = b.len[2]);
}

inline
Box&
Box::operator= (const Box& b)
{
  smallend = b.smallend;
  bigend = b.bigend;
  btype = b.btype;
//   D_EXPR6(len[0] = b.len[0],
//          len[1] = b.len[1],
//          len[2] = b.len[2]);
  return *this;
}

inline
IntVect
Box::sideEnd(Side::LoHiSide a_side) const
{
  IntVect retval;
  if (a_side == Side::Lo)
    retval = smallEnd();
  else
    retval = bigEnd();
  return retval;
}

inline
const IntVect&
Box::smallEnd () const
{
  return smallend;
}

inline
int
Box::smallEnd (int dir) const
{
  return smallend[dir];
}

inline
const IntVect&
Box::bigEnd () const
{
  return bigend;
}

inline
int
Box::bigEnd (int dir) const
{
  return bigend[dir];
}

inline
IndexType
Box::ixType () const
{
  return btype;
}

inline
IntVect
Box::type () const
{
  return btype.ixType();
}

inline
IndexType::CellIndex
Box::type (int dir) const
{
  return btype.ixType(dir);
}

inline
IntVect
Box::size () const
{
  return IntVect(D_DECL6(bigend[0]-smallend[0]+1,
                         bigend[1]-smallend[1]+1,
                         bigend[2]-smallend[2]+1,
                         bigend[3]-smallend[3]+1,
                         bigend[4]-smallend[4]+1,
                         bigend[5]-smallend[5]+1));
}

inline
int
Box::size (int dir) const
{
  // return len[dir];
  return bigend[dir]-smallend[dir]+1;
}

inline
const int*
Box::loVect () const
{
  return smallend.getVect();
}

inline
const int*
Box::hiVect () const
{
  return bigend.getVect();
}

inline
const int*
Box::getVect () const
{
  return smallend.getVect();
}

inline
bool
Box::numPtsOK () const
{
  long ignore;
  return numPtsOK(ignore);
}

inline
bool
Box::isEmpty () const
{
  //    return numPts() == 0;
  return D_TERM6((bigend[0] < smallend[0]),
                 || (bigend[1] < smallend[1]),
                 || (bigend[2] < smallend[2]),
                 || (bigend[3] < smallend[3]),
                 || (bigend[4] < smallend[4]),
                 || (bigend[5] < smallend[5]) );
}

inline
bool
Box::ok () const
{
  //    bool typ_ok = btype.ok();
  //     bool len_ok = D_TERM6((len[0] == (bigend[0]-smallend[0]+1)),
  //                           && (len[1] == (bigend[1]-smallend[1]+1)),
  //                           && (len[2] == (bigend[2]-smallend[2]+1)));
  //    return (bigend >= smallend) && len_ok && typ_ok;
  return (bigend >= smallend);
}

inline
bool
Box::contains (const IntVect& p) const
{

  //return ( !isEmpty() && (p >= smallend && p <= bigend) );
  return (p >= smallend && p <= bigend);
}

inline
bool
Box::sameType (const Box &b) const
{
  return btype == b.btype;
}

inline
bool
Box::contains (const Box& b) const
{
  CH_assert(sameType(b));
  return ( b.isEmpty() ||
           (!isEmpty() && b.smallend >= smallend && b.bigend <= bigend) );
}

inline
bool
Box::sameSize (const Box& b) const
{
  CH_assert(sameType(b));
  return D_TERM6(bigend[0] -smallend[0] == b.bigend[0]-b.smallend[0],
                 && bigend[1] -smallend[1] == b.bigend[1]-b.smallend[1],
                 && bigend[2] -smallend[2] == b.bigend[2]-b.smallend[2],
                 && bigend[3] -smallend[3] == b.bigend[3]-b.smallend[3],
                 && bigend[4] -smallend[4] == b.bigend[4]-b.smallend[4],
                 && bigend[5] -smallend[5] == b.bigend[5]-b.smallend[5]);
}

inline
bool
Box::operator== (const Box& b) const
{
  return smallend == b.smallend && bigend == b.bigend && b.btype == btype;
}

inline
bool
Box::operator!= (const Box& b) const
{
  return !operator==(b);
}

inline
bool
Box::cellCentered () const
{
  return !btype.any();
}

inline
bool
Box::volumeOK () const
{
  return numPtsOK();
}

inline
long
Box::index (const IntVect& v) const
{
  long result = v.vect[0]-smallend.vect[0];
#if CH_SPACEDIM == 1
  // nothing more goes here
#elif   CH_SPACEDIM==2
  int len0=(bigend[0]-smallend[0]+1);
  result += len0*(v.vect[1]-smallend.vect[1]);
#elif CH_SPACEDIM==3
  long int len0=(bigend[0]-smallend[0]+1);
  long int len1=(bigend[1]-smallend[1]+1);
  long int dif1= v.vect[1]-smallend.vect[1];
  long int dif2= v.vect[2]-smallend.vect[2];
  result += len0*(dif1
                + dif2*len1);
#elif CH_SPACEDIM==4
  int len0=(bigend[0]-smallend[0]+1);
  int len1=(bigend[1]-smallend[1]+1);
  int len2=(bigend[2]-smallend[2]+1);
  result += len0*((v.vect[1]-smallend.vect[1])
                  +len1*((v.vect[2]-smallend.vect[2])
                         +len2*(v.vect[3]-smallend[3])));
#elif CH_SPACEDIM==5
  int len0=(bigend[0]-smallend[0]+1);
  int len1=(bigend[1]-smallend[1]+1);
  int len2=(bigend[2]-smallend[2]+1);
  int len3=(bigend[3]-smallend[3]+1);
  result += len0*(v.vect[1]-smallend.vect[1]
                  +len1*(v.vect[2]-smallend.vect[2]
                         +len2*(v.vect[3]-smallend[3]
                                +len3*(v.vect[4]-smallend[4]))));
#elif CH_SPACEDIM==6
  int len0=(bigend[0]-smallend[0]+1);
  int len1=(bigend[1]-smallend[1]+1);
  int len2=(bigend[2]-smallend[2]+1);
  int len3=(bigend[3]-smallend[3]+1);
  int len4=(bigend[4]-smallend[4]+1);
  result += len0*(v.vect[1]-smallend.vect[1]
                  +len1*(v.vect[2]-smallend.vect[2]
                         +len2*(v.vect[3]-smallend[3]
                                +len3*(v.vect[4]-smallend[4]
                                       +len4*(v.vect[5]-smallend[5])))));

#else
  DIM > 6 undefined!
#endif
  return result;
}

inline
void
Box::computeBoxLen ()
{
//   if (isEmpty())
//     {
//       len = IntVect::Zero;
//     }
//   else
//     {
//       D_EXPR6(len[0] = bigend[0]-smallend[0] + 1,
//              len[1] = bigend[1]-smallend[1] + 1,
//              len[2] = bigend[2]-smallend[2] + 1);
//     }
}

inline
void
Box::computeBoxLenNotEmpty()
{
//   D_EXPR6(len[0] = bigend[0]-smallend[0] + 1,
//          len[1] = bigend[1]-smallend[1] + 1,
//          len[2] = bigend[2]-smallend[2] + 1);
}

inline
Box&
Box::setSmall (const IntVect& sm)
{
  CH_assert (sm <= bigend);

  smallend = sm;
  computeBoxLen();
  return *this;
}

inline
Box&
Box::setSmall (int dir,
               int sm_index)
{
  CH_assert (sm_index <= bigend[dir]);

  smallend.setVal(dir,sm_index);
  computeBoxLen();
  return *this;
}

inline
Box&
Box::setBig (const IntVect& bg)
{
  CH_assert (bg >= smallend);

  bigend = bg;
  computeBoxLen();
  return *this;
}

inline
Box&
Box::setBig (int dir,
             int bg_index)
{
  CH_assert (bg_index >= smallend[dir]);

  bigend.setVal(dir,bg_index);
  computeBoxLen();
  return *this;
}

inline
Box&
Box::setRange (int dir,
               int sm_index,
               int n_cells)
{
  CH_assert (n_cells > 0);

  smallend.setVal(dir,sm_index);
  bigend.setVal(dir,sm_index+n_cells-1);
  computeBoxLen();
  return *this;
}

inline
Box&
Box::shift (int dir,
            int nzones)
{
  if (!isEmpty())
    {
      smallend.shift(dir,nzones);
      bigend.shift(dir,nzones);
    }
  return *this;
}

inline
Box&
Box::shift (const IntVect& iv)
{
  if (!isEmpty())
    {
      smallend.shift(iv);
      bigend.shift(iv);
    }
  return *this;
}

inline
Box&
Box::convert (const IntVect& typ)
{
  CH_assert(typ >= IntVect::Zero && typ <= IntVect::Unit);
  if (!isEmpty())
    {
      IntVect shft(typ - btype.ixType());
      bigend += shft;
    }
  btype = IndexType(typ);
  computeBoxLen();
  return *this;
}

inline
Box&
Box::surroundingNodes (int dir)
{
  if (!(btype[dir]))
    {
      if (!isEmpty())
        {
          bigend.shift(dir,1);
        }
      //
      // Set dir'th bit to 1 = IndexType::NODE.
      //
      btype.set(dir);
      computeBoxLen();
    }
  return *this;
}

inline
Box&
Box::enclosedCells (int dir)
{
  if (btype[dir])
    {
      if (!isEmpty())
        {
          bigend.shift(dir,-1);
        }
      //
      // Set dir'th bit to 0 = IndexType::CELL.
      //
      btype.unset(dir);
      computeBoxLen();
    }
  return *this;
}

inline
Box
surroundingNodes (const Box& b,
                  int        dir)
{
  Box bx(b);
  return bx.surroundingNodes(dir);
}

inline
Box
surroundingNodes (const Box& b)
{
  Box bx(b);
  return bx.surroundingNodes();
}

inline
Box
enclosedCells (const Box& b,
               int        dir)
{
  Box bx(b);
  return bx.enclosedCells(dir);
}

inline
Box
enclosedCells (const Box& b)
{
  Box bx(b);
  return bx.enclosedCells();
}

inline
Box&
Box::operator+= (const IntVect& v)
{
  if (!isEmpty())
    {
      smallend += v;
      bigend += v;
    }
  return *this;
}

inline
Box
Box::operator+  (const IntVect& v) const
{
  if (isEmpty())
    {
      return(Box().convert(btype));
    }
  else
    {
      IntVect small(smallend);
      small += v;
      IntVect big(bigend);
      big += v;
      return Box(small,big,btype);
    }
}

inline
Box&
Box::operator-= (const IntVect& v)
{
  if (!isEmpty())
    {
      smallend -= v;
      bigend -= v;
    }
  return *this;
}

inline
bool
Box::operator < (const Box& rhs) const
{
  if (smallend == rhs.smallend) return bigend.lexLT(rhs.bigend);
  return(!isEmpty() && (rhs.isEmpty() || smallend.lexLT(rhs.smallend)));
}

inline
Box
Box::operator-  (const IntVect& v) const
{
  if (isEmpty())
    {
      return(Box().convert(btype));
    }
  else
    {
      IntVect small = smallend;
      small -= v;
      IntVect big = bigend;
      big -= v;
      return Box(small,big,btype);
    }
}

inline
Box&
Box::grow (int i)
{
  if (!isEmpty())
    {
      smallend.diagShift(-i);
      bigend.diagShift(i);
      if (!(bigend >= smallend)) *this = Box().convert(btype);
      computeBoxLen();
    }
  return *this;
}

inline
Box
grow (const Box& b,
      int        i)
{
  if (b.isEmpty())
    {
      return Box().convert(b.btype);
    }

  IntVect small = diagShift(b.smallend,-i);
  IntVect big   = diagShift(b.bigend,i);
  if (!(big >= small))
    {
      return Box().convert(b.btype);
    }
  return Box(small,big,b.btype);
}

inline
Box&
Box::grow (const IntVect& v)
{
  if (!isEmpty())
    {
      smallend -= v;
      bigend   += v;
      if (!(bigend >= smallend)) *this = Box().convert(btype);
      computeBoxLen();
    }
  return *this;
}

inline
Box
grow (const Box&     b,
      const IntVect& v)
{
  if (b.isEmpty())
    {
      return Box().convert(b.btype);
    }
  else
    {
      IntVect small = b.smallend - v;
      IntVect big   = b.bigend   + v;
      if (!(big >= small))
        {
          return Box().convert(b.btype);
        }
      else
        {
          return Box(small,big,b.btype);
        }
    }
}

inline
Box&
Box::grow (int idir,
           int n_cell)
{
  if (!isEmpty())
    {
      smallend.shift(idir, -n_cell);
      bigend.shift(idir, n_cell);
      if (!(bigend >= smallend)) *this = Box().convert(btype);
      computeBoxLen();
    }
  return *this;
}

inline
Box&
Box::growLo (int idir,
             int n_cell)
{
  if (!isEmpty())
    {
      smallend.shift(idir, -n_cell);
      if (!(bigend >= smallend)) *this = Box().convert(btype);
      computeBoxLen();
    }
  return *this;
}

inline
Box&
Box::growDir (int idir,
              const Side::LoHiSide& a_sd,
              int n_cell)
{
  if (a_sd == Side::Lo)
    {
      growLo(idir, n_cell);
    }
  else
    {
      growHi(idir, n_cell);
    }
  return *this;
}

inline
Box&
Box::growHi (int idir,
             int n_cell)
{
  if (!isEmpty())
    {
      bigend.shift(idir,n_cell);
      if (!(bigend >= smallend)) *this = Box().convert(btype);
      computeBoxLen();
    }
  return *this;
}

inline
Box
Box::operator& (const Box& rhs) const
{
  Box lhs(*this);
  return lhs &= rhs;
}

#endif /* WRAPPER */

#include "NamespaceFooter.H"
#endif /*CH_BOX_H*/
