#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PROBLEMDOMAIN_H_
#define _PROBLEMDOMAIN_H_

#ifndef WRAPPER
#include <iostream>

#include "Vector.H"
#include "IntVect.H"
#include "Box.H"
#include "Misc.H"
#endif

#include "SPACE.H"
#include "NamespaceHeader.H"

/// Class to manage box-shifting used to enforce periodic BC's
/** The ShiftIterator class contains a list of shift vectors necessary
    to enforce periodic boundary conditions (enforced by shifting the
    Box of each Fab, then copying any overlapping valid cells with
    ghost cells of the shifted Fab.  Depending on the number of periodic
    directions, the list of shift vectors over which to iterate differs.
*/
class ShiftIterator
{
public:
  /// Basic constructor
  ShiftIterator();

  /// Defining Constructor
  /**
     Builds a ShiftIterator based on the periodicity in
     \em a_isPeriodic
  */
  ShiftIterator(const bool* a_isPeriodic);

  /// Copy constructor
  ShiftIterator(const ShiftIterator& a_shiftIt);

  /// Destructor
  ~ShiftIterator();

  /// Assignment operator
  ShiftIterator& operator=(const ShiftIterator& a_src);

  /// Recompute shift vectors based on periodic directions
  void computeShifts(const bool* a_isPeriodic);

  /// Returns the current shift unit vector
  inline IntVect operator()() const;

  /// Equivalent to the () operator
  IntVect i() const
  {
    return this->operator()();
  }
  /// Increment to the next shift unit vector
  inline void operator++();

  /// Equivalent to the ++ operator
  void incr()
  {
    ++(*this);
  }

  /// Is the iterator still within its range of shift vectors?
  inline bool ok() const;

  ///which periodic image is this
  int index() const
  {
    return m_index;
  }

  const IntVect& operator[](int index) const
  {
    return  m_shift_vectors[index];
  }
  /// Reset to first shift unit vector
  inline void reset();

  /// Equivalent to reset()
  inline void begin();

  /// Skip the iterator to the end.
  /**
     ok() will return false after this method is called.
  */
  void end();

private:

  int m_index;

  Vector<IntVect> m_shift_vectors;

};

/// A class to facilitate interaction with physical boundary conditions
/**

  ProblemDomain is a class which facilitates the application of physical
  boundary conditions, both periodic and non-periodic.  This class contains
  much of the functionality of the Box class, since logically the
  computational domain is generally a Box.

  Intersection with a ProblemDomain object will result in only removing
  regions which are outside the physical domain in non-periodic directions.
  Regions outside the logical computational domain in periodic directions will
  be treated as ghost cells which can be filled with an exchange() function
  or through suitable interpolation from a coarser domain.

  Since ProblemDomain will contain a Box, it is a dimension dependent class,
  so SpaceDim must be defined as either 1, 2, or 3 when compiling.

  Note that this implementation of ProblemDomain is inherently
  cell-centered.

*/

class ProblemDomain
{
public:

  // Constructors

  /// The default constructor.  The constructed domain box is empty.
  /**
   */
  ProblemDomain ();

  /// Construct ProblemDomain with \em a_domBox as computational domain
  /**
     This constructor defaults to non-periodic domain
  */
  ProblemDomain(const Box& a_domBox);

  /// Construct ProblemDomain with a_domBox as computational domain.
  /**
     \em a_isPeriodic is a SpaceDim array of bools; if true, the
     physical boundary condition is periodic in the coordinate direction.
     False means a non-periodic BC.
  */
  ProblemDomain(const Box& a_domBox, const bool* a_isPeriodic);

  /// Construct a ProblemDomain.
  /**
     It is an error if small is greater than big.
     Defaults to non-periodic domain.
  */
  ProblemDomain (const IntVect& small,
                 const IntVect& big);

  /// Construct a ProblemDomain.
  /**
      It is an error if small is greater than big.
      \em a_isPeriodic is a SpaceDim array of bools; if true, the
      physical boundary condition is periodic in the coordinate direction.
      False means a non-periodic BC.
  */
  ProblemDomain (const IntVect& small,
                 const IntVect& big,
                 const bool* a_isPeriodic);

  /// Construct ProblemDomain with specified lengths.
  /**
     It is an error if the lengths are negative.
     Defaults to non-periodic domain.
  */
  ProblemDomain (const IntVect& small,
                 const int*     vec_len);

  /// Construct ProblemDomain with specified lengths.
  /**
     It is an error if the lengths are negative.
      \em a_isPeriodic is a SpaceDim array of bools; if true, the
      physical boundary condition is periodic in the coordinate direction.
      False means a non-periodic BC.
  */
  ProblemDomain (const IntVect& small,
                 const int*     vec_len,
                 const bool* a_isPeriodic);

  /// The copy constructor.
  /**
   */
  ProblemDomain (const ProblemDomain& a_src);

  /// Construct ProblemDomain with \em a_domBox as computational domain
  /**
     This constructor defaults to non-periodic domain
  */
  void define(const Box& a_domBox);

  /// Construct ProblemDomain with a_domBox as computational domain.
  /**
     \em a_isPeriodic is a SpaceDim array of bools; if true, the
     physical boundary condition is periodic in the coordinate direction.
     False means a non-periodic BC.
  */
  void define(const Box& a_domBox, const bool* a_isPeriodic);

  /// Construct a ProblemDomain.
  /**
     It is an error if small is greater than big.
     Defaults to non-periodic domain.
  */
  void define (const IntVect& small,
               const IntVect& big);

  /// Construct a ProblemDomain.
  /**
      It is an error if small is greater than big.
      \em a_isPeriodic is a SpaceDim array of bools; if true, the
      physical boundary condition is periodic in the coordinate direction.
      False means a non-periodic BC.
  */
  void define (const IntVect& small,
               const IntVect& big,
               const bool* a_isPeriodic);

  /// Construct ProblemDomain with specified lengths.
  /**
     It is an error if the lengths are negative.
     Defaults to non-periodic domain.
  */
  void  define (const IntVect& small,
                const int*     vec_len);

  /// Construct ProblemDomain with specified lengths.
  /**
     It is an error if the lengths are negative.
      \em a_isPeriodic is a SpaceDim array of bools; if true, the
      physical boundary condition is periodic in the coordinate direction.
      False means a non-periodic BC.
  */
  void define (const IntVect& small,
               const int*     vec_len,
               const bool* a_isPeriodic);

  /// The copy constructor.
  /**
   */
  void define (const ProblemDomain& a_src);

  //  Accessors

  /// Returns the logical computational domain
  /**
   */
  const Box& domainBox() const;

  /// Returns true if BC is periodic in direction a_dir
  /**
  */
  bool isPeriodic(int a_dir) const;

  /// Returns true is BC is periodic in _any_ direction
  /**
   */
  bool isPeriodic() const;

  /// Returns the shiftIterator for this ProblemDomain
  /**
      The shiftIterator is defined based on the periodicity of this
      ProblemDomain.
  */
  ShiftIterator shiftIterator() const;

  /// Returns true if this ProblemDomain is empty or undefined.
  /**
  */
  bool isEmpty () const;

  /// Return the size of the domainBox in the specified coordinate direction.
  /**
   */
  int size (const int& a_idir) const;

  /// Return the size of the domainBox.
  /**
   */
  IntVect size () const;

  /// Returns true if argument is contained within this ProblemDomain.
  /**
     An empty ProblemDomain does not contain and is not contained by
     any ProblemDomain, including itself.  Note also that in a periodic
     Direction, any index is valid, since the domain is infinite in that
     direction.
  */
  bool contains (const IntVect& p) const;

  /// Returns the periodic image of this IntVect inside of the ProblemDomain.
  /**
      Return true if the domain contains this IntVect, returning the image
      in 'p', otherwise returns false
  */
  bool image (IntVect& p) const;

  /// Returns true if argument is contained within this ProblemDomain.
  /**
     An empty ProblemDomain does not contain any Box.  An entirely periodic
     domain contains any Box.
  */
  bool contains (const Box& b) const;

  /// Equivalent to \em contains() function.
  bool contains_box(const Box& b) const
  {
    return contains(b);
  }

  /// Returns true if this ProblemDomain and the argument intersect.
  /**
     It is an error if a_box is not cell-centered.
     An empty ProblemDomain does not intersect any Box.
     This will do nothing in periodic directions (since a periodic domain
     is infinite in the periodic direction.  If periodic in all dimensions,
     this will always return true.
  */
  bool intersects (const Box& a_box) const;

  /// Returns true if this ProblemDomain and the argument intersect.
  /**
     It is an error if a_box is not cell-centered.
     This routine does not perform the check to see if *this or b are
     empty Boxes.  It is the callers responsibility to ensure that
     this never happens.  If you are unsure, the use the .intersects(..)
     routine.  In periodic directions, will always return true.

  */
  bool intersectsNotEmpty (const Box& a_box) const;

  /// Returns true if this argument is adjacent to a periodic boundary
  /**
   */
  bool periodicAdjacent(const Box& a_box) const;

  ///
  /**
   */
  void insertImages(std::list<Box>& a_list, const Box& a_box) const;

  /// Returns true if box1 and box2 or any of their periodic images intersect
  /**
      (useful for checking disjointness)
  */
  bool intersects(const Box& box1, const Box& box2) const;

  /// Returns true if the two domain are equivalent
  /**
   */
  bool operator== (const ProblemDomain& a_otherDomain) const;

  /// Returns true if the two domain are not equivalent
  /**
   */
  bool operator!= (const ProblemDomain& a_otherDomain) const;

  /// Modifies a_box to be the intersection of a_box and \em a_probdomain
  /**
  */
  friend void operator &=(Box& a_box, const ProblemDomain& a_probdomain);

  /// Returns a Box which is the interesection of a_box and \em a_probdomain
  /**
  */
  friend Box operator & (const Box& a_box, const ProblemDomain& a_probdomain);

  //  Modification Functions

  /// The assignment operator.
  /**
   */
  ProblemDomain& operator= (const ProblemDomain& b);

  /// Sets whether BC is periodic in direction a_dir (true == periodic)
  /**
  */
  void setPeriodic(int a_dir, bool a_isPeriodic);

  /// Grows (or shrinks) the domain Box by i in all directions
  /**
   */
  inline ProblemDomain& grow (int i);

  /// Returns a ProblemDomain with a domainBox grown by the given amount
  /**
   */
  friend inline ProblemDomain grow(const ProblemDomain& pd,
                                   int i);

  /// Grow this ProblemDomain
  /**
       Modifies this ProblemDomain by growing the domainBox in each
       direction by the specified amount
  */
  inline ProblemDomain& grow(const IntVect& v);

  /// Returns a grown version of \em pd.
  /**
     Returns a ProblemDomain that is the argument ProblemDomain
     with a DomainBox grown by the given amount.
  */
  friend inline ProblemDomain grow (const ProblemDomain& pd,
                                     const IntVect& v);

  /// Grow this ProblemDomain
  /**
      Modifies this ProblemDomain by growing it on the low and high end
      by n_cell cells in direction idir.
  */
  inline ProblemDomain& grow(int idir, int n_cell);

  /// Grow this ProblemDomain on the low side
  /**
      Modifies this ProblemDomain by growing it on the low end by n_cell
      cells in direction idir.
  */
  inline ProblemDomain& growLo(int idir, int n_cell=1);

  /// Grow this ProblemDomain on the high side
  /** Modifies this ProblemDomain by growing it on the high end by n_Cell
      cells in direction idir
  */
  inline ProblemDomain& growHi(int idir, int n_cell=1);

  /// Returns a face-centered Box at the low side of \em a_pd
  /**
     Returns the edge-centered Box (in direction \em a_dir) defining
     the low side of the argument ProblemDomain.  The output Box will
     have the given length in the given direction.  Directions are
     zero-based.  It is an error if not 0 <= dir < SpaceDim.  The neighbor
     of an Empty ProblemDomain is an Empty Box of the appropriate type.
     If \em a_dir is a periodic direction, will return an empty Box.
  */
  friend Box bdryLo (const ProblemDomain& a_pd,
                     int        a_dir,
                     int        a_len);

  friend Box bdryLo (const ProblemDomain& a_pd,
                     int        a_dir) {return bdryLo(a_pd,a_dir,1);}


  /// Returns a face-centered Box at the high side of \em a_pd
  /**
     Returns the edge-centered Box (in direction \em a_dir) defining the
     high side of the argument ProblemDomain.  The return Box will have the
     given length in the given direction.  Directions are zero-based.
     It is an error if not 0 <= dir < SpaceDim.  The neighbor of an
     Empty ProblemDomain is an Empty Box of the appropriate type.
     If \em a_dir is a periodic direction, will return an empty Box.
  */
  friend Box bdryHi (const ProblemDomain& a_pd,
                     int        a_dir,
                     int        a_len);

  friend Box bdryHi (const ProblemDomain& a_pd,
                     int        a_dir) {return bdryHi(a_pd,a_dir,1);}
  
  /// Returns the cell-centered Box adjacent to the low side of \em a_pd.
  /**
     Returns the cell centered Box of the given length adjacent to the
     argument ProblemDomain on the low end along the given coordinate direction.
     The return Box is identical to the argument ProblemDomain in the other
     directions.  The return ProblemDomain and the argument ProblemDomain
     have an empty intersection.

     \b NOTES:
     - len >= 1.
     - Box retval = adjCellLo(b,dir,len) is equivalent to the following set of
       operations:
       \code
         Box retval(b);

         retval.convert(dir,ProblemDomain::CELL);

         retval.setrange(dir,retval.smallEnd(dir)-len,len);
        \endcode

     Directions are zero-based.  It is an error if not 0 <= dir <
     SpaceDim.  The neighbor of an Empty ProblemDomain is an Empty Box of the
     appropriate type. If \em a_dir is a periodic direction, will return
     an empty Box as well.

  */
  friend Box adjCellLo (const ProblemDomain& a_pd,
                        int        a_dir,
                        int        a_len);

  friend Box adjCellLo (const ProblemDomain& a_pd,
                        int        a_dir) {return adjCellLo(a_pd,a_dir,1);}

  /// Returns the cell-centered Box adjacent to the low side of \em a_pd.
  /**
     Returns the cell centered Box of the given length adjacent to the
     argument ProblemDomain on the high end along the given coordinate
     direction.  The return Box is identical to the argument ProblemDomain in
     the other directions.  The return Box and the argument ProblemDomain have
     an empty intersection.

     \b NOTES:
     - len >= 1.

     - Box retval = adjCellHi(b,dir,len) is equivalent to the
       following set of operations:
       \code
         Box retval(b);

         retval.convert(dir,ProblemDomain::CELL);

         retval.setrange(dir,retval.bigEnd(dir)+1,len);
        \endcode

     Directions are zero-based.  It is an error if not 0 <= dir <
     SpaceDim.  The neighbor of an Empty ProblemDomain is an Empty Box of the
     appropriate type.    If \em a_dir is a periodic direction, will return
     an empty Box.

  */
  friend Box adjCellHi (const ProblemDomain& a_pd,
                        int        a_dir,
                        int        a_len);

  friend Box adjCellHi (const ProblemDomain& a_pd,
                        int        a_dir) {return adjCellHi(a_pd,a_dir,1);}


  // Intersection functions

  /// Returns the Box intersection of this ProblemDomain and \em a_b
  /**
     The intersection of the Empty ProblemDomain and any Box is the Empty
     Box. This operator does nothing in periodic directions (since
     a periodic domain is an infinite domain).

  */
  Box operator& (const Box& a_b) const;

  // refinement

  /// Refine this problem domain.
  /**
     Modifies this ProblemDomain by refining it by given (positive) refinement
     ratio.  The Empty ProblemDomain is not modified by this function.
  */
    ProblemDomain& refine (int a_refinement_ratio);

  /// Return a ProblemDomain which is a refinement of \em a_probdomain
  /**
     Returns a ProblemDomain that is the argument ProblemDomain refined by
     given (positive) refinement ratio.  The Empty ProblemDomain is not
     modified by this function.
  */
  friend ProblemDomain refine (const ProblemDomain& a_probdomain,
                               int   a_refinement_ratio);

  /// Refine this ProblemDomain
  /**
     Modifies this ProblemDomain by refining it by given (positive) refinement
     ratio.  The Empty ProblemDomain is not modified by this function.
  */
  ProblemDomain& refine (const IntVect& a_refinement_ratio);

  /// Refinement function
  /**
     Returns a ProblemDomain that is the argument ProblemDomain refined
     by given (positive) refinement ratio.  The Empty ProblemDomain is
     not modified by this function.

  */
  friend ProblemDomain refine (const ProblemDomain&     a_probdomain,
                            const IntVect& a_refinement_ratio);

  // coarsening

  /// Coarsen this ProblemDomain
  /**
     Modifies this ProblemDomain by coarsening it by given (positive)
     refinement ratio.  The Empty ProblemDomain is not modified by
     this function.
  */

  ProblemDomain& coarsen (int a_refinement_ratio);

  /// Coarsening function
  /**
     Returns a ProblemDomain that is the argument ProblemDomain coarsened
     by given (positive) refinement ratio.  The Empty ProblemDomain is not
     modified by this function.
  */
  friend ProblemDomain coarsen (const ProblemDomain& a_probdomain,
                             int        a_refinement_ratio);

  /// Coarsen this ProblemDomain
  /**
     Modifies this ProblemDomain by coarsening by given (positive) refinement
     ratio.  The Empty ProblemDomain is not modified by this function.
  */
  ProblemDomain& coarsen (const IntVect& refinement_ratio);

  /// Coarsening function
  /**
     Returns a ProblemDomain that is the argument ProblemDomain coarsened
     by given (positive) refinement ratio.  The Empty ProblemDomain is
     not modified by this function.

  */
  friend ProblemDomain coarsen (const ProblemDomain&  a_probdomain,
                                const IntVect&        a_refinement_ratio);

  ///
  void shift(const IntVect& a_shift)
  {
    m_domainBox.shift(a_shift);
  }

  bool operator<(const ProblemDomain& rhs) const
  {
    return m_domainBox < rhs.m_domainBox;
  }
  // I/O Functions

  /// Write an ASCII representation to the ostream.
  /**
  */
  friend std::ostream& operator<< (std::ostream&   os,
                                   const ProblemDomain& bx);

  /// Read from istream.
  /**
  */
  friend std::istream& operator>> (std::istream& is,
                                   ProblemDomain&     bx);

  void shiftIt(Box& a_box, int shiftIndex) const ;
  void unshiftIt(Box& a_box, int shiftIndex) const ;
  /// Gives more detail than printOn.
  /**
     Useful for exiting due to an error.
  */
  void dumpOn (std::ostream& strm) const;


protected:
  friend class HDF5Handle;

  /**
     Periodicity info
  */
  bool m_isPeriodic[SpaceDim];

  /**
     Domain index extents
  */
  Box m_domainBox;

  /**
     Shift iterator for this ProblemDomain
  */
  ShiftIterator m_shiftIt;
};

class ImageIterator
{
public:
  ImageIterator(const ProblemDomain& a_domain)
  {
    define(a_domain);
  }

  void define(const ProblemDomain& a_domain);

  void begin(const Box& a_box)
  {
    m_counter=-1;
    m_box = a_box;
    this->operator++();
  }

  void operator++();

  bool ok()
  {
    return m_shifter[m_counter] != IntVect::Zero;
  }

  const Box& box() const
  {
    return m_current;
  }

  const ProblemDomain& domain() const
  {
    return m_domain;
  }

  void checkDefine(const ProblemDomain& a_domain)
  {
    if (!(m_domain == a_domain)) define(a_domain);
  }

protected:
  ProblemDomain m_domain;
  Box m_quadrant[D_TERM6(3,*3,*3,*3,*3,*3)];
  IntVect m_shifter[D_TERM6(3,*3,*3,*3,*3,*3)];
  Box m_box;
  Box m_current;
  int m_counter;
};

//
// Inlines.
//
#ifndef WRAPPER

inline
ShiftIterator::ShiftIterator()
  : m_index(100), m_shift_vectors()
{
}

inline
ShiftIterator::ShiftIterator(const ShiftIterator& a_src)
{
  m_index = a_src.m_index;
  m_shift_vectors = a_src.m_shift_vectors;
}

inline
ShiftIterator&
ShiftIterator::operator=(const ShiftIterator& a_src)
{
  m_index = a_src.m_index;
  m_shift_vectors = a_src.m_shift_vectors;
  return *this;
}

inline
IntVect
ShiftIterator::operator()() const
{
  CH_assert(ok());
  return m_shift_vectors[m_index];
}

inline
void
ShiftIterator::operator++()
{
  m_index++;
}

inline
bool
ShiftIterator::ok() const
{
  return (m_index < m_shift_vectors.size());
}

inline
void
ShiftIterator::reset()
{
  m_index = 0;
}

inline
void
ShiftIterator::begin()
{
  m_index = 0;
}

inline
void
ShiftIterator::end()
{
  m_index = m_shift_vectors.size();
}

inline
ProblemDomain::ProblemDomain ()
{
  // default is a non-periodic domain
  for (int dir=0; dir<SpaceDim; dir++)
    {
      m_isPeriodic[dir] = false;
    }

}

inline
ProblemDomain::ProblemDomain (const ProblemDomain& b)
  : m_domainBox(b.m_domainBox), m_shiftIt(b.m_shiftIt)
{
  for (int dir=0; dir<SpaceDim; dir++)
    {
      m_isPeriodic[dir] = b.m_isPeriodic[dir];
    }
}

inline void
ProblemDomain::define(const ProblemDomain& b)

{
  m_domainBox=(b.m_domainBox);
  m_shiftIt  = (b.m_shiftIt);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      m_isPeriodic[dir] = b.m_isPeriodic[dir];
    }
}

inline
bool
ProblemDomain::operator== (const ProblemDomain& a_otherDomain) const
{
  bool result = true;

  if (m_domainBox != a_otherDomain.m_domainBox)
  {
    result = false;
  }
  else
  {
    for (int dir=0; dir<SpaceDim; dir++)
    {
      if (m_isPeriodic[dir] != a_otherDomain.m_isPeriodic[dir])
      {
        result = false;
        break;
      }
    }
  }

  return result;
}

inline
bool
ProblemDomain::operator!= (const ProblemDomain& a_otherDomain) const
{
  return !(*this == a_otherDomain);
}

inline
ProblemDomain&
ProblemDomain::operator= (const ProblemDomain& b)
{
  m_domainBox = b.m_domainBox;
  for (int dir=0; dir<SpaceDim; dir++)
    {
      m_isPeriodic[dir] = b.m_isPeriodic[dir];
    }
  m_shiftIt = b.m_shiftIt;
  return *this;
}

inline void
ProblemDomain::shiftIt(Box& a_box, int a_shiftIndex) const
{
  a_box.shift(m_shiftIt[a_shiftIndex]*m_domainBox.size());
}

inline void
ProblemDomain::unshiftIt(Box& a_box, int a_shiftIndex) const
{
  a_box.shift(- m_shiftIt[a_shiftIndex]*m_domainBox.size());
}

inline
const Box&
ProblemDomain::domainBox() const
{
  return m_domainBox;
}

inline
bool
ProblemDomain::isPeriodic(int a_dir) const
{
  return m_isPeriodic[a_dir];
}

inline
bool
ProblemDomain::isPeriodic() const
{
  return D_TERM6(m_isPeriodic[0], ||
                 m_isPeriodic[1], ||
                 m_isPeriodic[2], ||
                 m_isPeriodic[3], ||
                 m_isPeriodic[4], ||
                 m_isPeriodic[5]);
}

inline
ProblemDomain&
ProblemDomain::grow(int i)
{
  m_domainBox.grow(i);
  return *this;
}

inline
ProblemDomain
grow(const ProblemDomain& pd, int i)
{
  ProblemDomain newPd(pd);
  newPd.grow(i);
  return newPd;
}

inline
ProblemDomain&
ProblemDomain::grow(const IntVect& v)
{
  m_domainBox.grow(v);
  return *this;
}

inline
ProblemDomain
grow(const ProblemDomain& pd, const IntVect& v)
{
  ProblemDomain newPd(pd);
  newPd.grow(v);
  return newPd;
}

inline
ProblemDomain&
ProblemDomain::grow(int idir, int n_cell)
{
  m_domainBox.grow(idir, n_cell);
  return *this;
}

inline
ProblemDomain&
ProblemDomain::growLo(int idir, int n_cell)
{
  m_domainBox.growLo(idir, n_cell);
  return *this;
}

inline
ProblemDomain&
ProblemDomain::growHi(int idir, int n_cell)
{
  m_domainBox.growHi(idir, n_cell);
  return *this;
}

inline
ShiftIterator
ProblemDomain::shiftIterator() const
{
  return m_shiftIt;
}

inline
bool
ProblemDomain::isEmpty () const
{
  return (m_domainBox.isEmpty());
}

inline
int
ProblemDomain::size (const int& a_idir) const
{
  return (m_domainBox.size(a_idir));
}

inline
IntVect
ProblemDomain::size () const
{
  return (m_domainBox.size());
}

inline
bool
ProblemDomain::contains (const IntVect& p) const
{

  // boy is this ugly!
  return (   !isEmpty()
          && (D_TERM6((m_isPeriodic[0]
                       || (   p[0] >= m_domainBox.smallEnd(0)
                              && p[0] <= m_domainBox.bigEnd  (0))),
          &&         (m_isPeriodic[1]
                      || (   p[1] >= m_domainBox.smallEnd(1)
                          && p[1] <= m_domainBox.bigEnd  (1))),
          &&         (m_isPeriodic[2]
                      || (   p[2] >= m_domainBox.smallEnd(2)
                             && p[2] <= m_domainBox.bigEnd  (2))),
          &&         (m_isPeriodic[3]
                       || (   p[3] >= m_domainBox.smallEnd(3)
                              && p[3] <= m_domainBox.bigEnd  (3))),
          &&         (m_isPeriodic[4]
                      || (   p[4] >= m_domainBox.smallEnd(4)
                          && p[4] <= m_domainBox.bigEnd  (4))),
          &&         (m_isPeriodic[5]
                      || (   p[5] >= m_domainBox.smallEnd(5)
                          && p[5] <= m_domainBox.bigEnd  (5))))));
}

inline
bool
ProblemDomain::image(IntVect& p) const
{
  if (m_domainBox.contains(p)) return true;
  if (!contains(p)) return false;

  D_TERM6(
          if (m_isPeriodic[0])
  {
    if (p[0]<m_domainBox.smallEnd(0))    p[0]+= m_domainBox.size(0);
    else if (p[0]>m_domainBox.bigEnd(0)) p[0]-= m_domainBox.size(0);
  },
          if (m_isPeriodic[1])
  {
    if (p[1]<m_domainBox.smallEnd(1))    p[1]+= m_domainBox.size(1);
    else if (p[1]>m_domainBox.bigEnd(1)) p[1]-= m_domainBox.size(1);
  },
          if (m_isPeriodic[2])
  {
    if (p[2]<m_domainBox.smallEnd(2))    p[2]+= m_domainBox.size(2);
    else if (p[2]>m_domainBox.bigEnd(2)) p[2]-= m_domainBox.size(2);
  },
          if (m_isPeriodic[3])
  {
    if (p[3]<m_domainBox.smallEnd(3))    p[3]+= m_domainBox.size(3);
    else if (p[3]>m_domainBox.bigEnd(3)) p[3]-= m_domainBox.size(3);
  },
          if (m_isPeriodic[4])
  {
    if (p[4]<m_domainBox.smallEnd(4))    p[4]+= m_domainBox.size(4);
    else if (p[4]>m_domainBox.bigEnd(4)) p[4]-= m_domainBox.size(4);
  },
          if (m_isPeriodic[5])
  {
    if (p[5]<m_domainBox.smallEnd(5))    p[5]+= m_domainBox.size(5);
    else if (p[5]>m_domainBox.bigEnd(5)) p[5]-= m_domainBox.size(5);
  });

  return true;
}

inline
bool
ProblemDomain::contains (const Box& b) const
{
  // boy is this ugly!
  if (b.type() == m_domainBox.type())
  {
  return (   !isEmpty()
             && (D_TERM6((m_isPeriodic[0]
                          || (   b.smallEnd(0) >= m_domainBox.smallEnd(0)
                                 && b.bigEnd  (0) <= m_domainBox.bigEnd  (0))), &&
                         (m_isPeriodic[1]
                          || (   b.smallEnd(1) >= m_domainBox.smallEnd(1)
                                 && b.bigEnd  (1) <= m_domainBox.bigEnd  (1))), &&
                         (m_isPeriodic[2]
                          || (   b.smallEnd(2) >= m_domainBox.smallEnd(2)
                                 && b.bigEnd  (2) <= m_domainBox.bigEnd  (2))), &&
                         (m_isPeriodic[3]
                          || (   b.smallEnd(3) >= m_domainBox.smallEnd(3)
                                 && b.bigEnd  (3) <= m_domainBox.bigEnd  (3))), &&
                         (m_isPeriodic[4]
                          || (   b.smallEnd(4) >= m_domainBox.smallEnd(4)
                                 && b.bigEnd  (4) <= m_domainBox.bigEnd  (4))), &&
                         (m_isPeriodic[5]
                          || (   b.smallEnd(5) >= m_domainBox.smallEnd(5)
                                 && b.bigEnd  (5) <= m_domainBox.bigEnd  (5))))));
  }
  else
  {
    Box domainBox = m_domainBox;

    // Check b's centering and adjust intersectBox as needed
    for (int dir = 0; dir < SpaceDim; dir++)
    {
      if (b.type(dir) != domainBox.type(dir))
      {
        if (b.type(dir) == IndexType::NODE)
        {
          domainBox.surroundingNodes(dir);
        }
        else
        {
          domainBox.enclosedCells(dir);
        }
      }
    }

    return (   !isEmpty()
            && (D_TERM6((m_isPeriodic[0]
                         || (   b.smallEnd(0) >= domainBox.smallEnd(0)
                                && b.bigEnd  (0) <= domainBox.bigEnd  (0))), &&
                        (m_isPeriodic[1]
                         || (   b.smallEnd(1) >= domainBox.smallEnd(1)
                                && b.bigEnd  (1) <= domainBox.bigEnd  (1))), &&
                        (m_isPeriodic[2]
                         || (   b.smallEnd(2) >= domainBox.smallEnd(2)
                                && b.bigEnd  (2) <= domainBox.bigEnd  (2))), &&
                        (m_isPeriodic[3]
                         || (   b.smallEnd(3) >= domainBox.smallEnd(3)
                                && b.bigEnd  (3) <= domainBox.bigEnd  (3))), &&
                        (m_isPeriodic[4]
                         || (   b.smallEnd(4) >= domainBox.smallEnd(4)
                                && b.bigEnd  (4) <= domainBox.bigEnd  (4))), &&
                        (m_isPeriodic[5]
                         || (   b.smallEnd(5) >= domainBox.smallEnd(5)
                                && b.bigEnd  (5) <= domainBox.bigEnd  (5))))));
  }
}

#endif /*WRAPPER*/

#include "NamespaceFooter.H"
#endif
