#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INTVECT_H_
#define _INTVECT_H_

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "SPACE.H"
#include "Vector.H"
#include "IndexTM.H"
#include "Misc.H"

namespace BLfacade
{
  class IntVect;
}

#include "NamespaceHeader.H"

//template<typename T, int n>
//class IndexTM;

class HDF5Handle;

/// An integer Vector in SpaceDim-dimensional space
/**
  The class IntVect is an implementation of an integer vector in a
  SpaceDim-dimensional space.  It represents a point in a discrete space.
  IntVect values are accessed using the operator[] function, as for a normal
  C++ array.  In addition, the basic arithmetic operators have been overloaded
  to implement scaling and translation operations.
*/

class IntVect
{
public:

  /**
     \name Constructors and Accessors
  */
  /*@{*/

  ///
  /**
     Construct an IntVect whose components are uninitialized.
  */
  IntVect()
  {}

  explicit IntVect( const Vector<int>& vi)
  {
    D_EXPR6(vect[0]=vi[0], vect[1]=vi[1], vect[2] = vi[2],
            vect[3]=vi[3], vect[4]=vi[4], vect[5] = vi[5]);
  }

  ///
  /**
     Destructor.
  */
  ~IntVect()
  {}

  ///
  /**
     Construct an IntVect given the specific values for its
     coordinates.  D_DECL6 is a macro that sets the constructor to
     take CH_SPACEDIM arguments.
  */
  explicit IntVect (D_DECL6(int i, int j, int k,
                            int l, int m, int n));

  ///
  /**
     Construct an IntVect setting the coordinates to the corresponding
     values in the integer array <i>a</i>.
  */
  explicit IntVect (const int* a);

  ///
  /**
     The copy constructor.
  */
  IntVect (const IntVect& rhs);

  IntVect (const IndexTM<int, CH_SPACEDIM>& a_tm);

  IntVect copy() const
  {
    return *this;
  }
  ///
  /**
     The assignment operator.
  */
  IntVect& operator= (const IntVect& rhs);

  ///
  /**
     Returns a modifiable lvalue reference to the <i>i</i>'th coordinate of the
     IntVect.
  */
  inline
  int& operator[] (int i);

  ///
  /**
     Returns the <i>i</i>'th coordinate of the IntVect.
  */
  inline
  int operator[] (int i) const;

  ///
  /**
     Set <i>i</i>'th coordinate of IntVect to <i>val</i>.
  */
  void setVal (int i,
               int val);

  /*@}*/

  /**
     \name Data pointer functions
  */
  /*@{*/

  ///
  /**
     Returns a const pointer to an array of coordinates of the IntVect.
     Useful for arguments to FORTRAN calls.
  */
  const int* getVect () const;

  ///
  /**
     Only for sending to Fortran
   */
  const int*  dataPtr() const;

  ///
  /**
     Only for sending to Fortran
   */
  int*  dataPtr();

  /*@}*/

  /**
     \name Comparison Operators
  */
  /*@{*/

  ///
  /**
     Returns true if this IntVect is equivalent to argument IntVect.  All
     comparisons between analogous components must be satisfied.
  */
  bool operator== (const IntVect& p) const;

  ///
  /**
     Returns true if this IntVect is different from argument IntVect.
     All comparisons between analogous components must be satisfied.
  */
  bool operator!= (const IntVect& p) const;

  ///
  /**
     Returns true if this IntVect is less than argument IntVect.  All
     comparisons between analogous components must be satisfied.  Note
     that, since the comparison is component-wise, it is possible for
     an IntVect to be neither greater than, less than, nor equal to
     another.
  */
  bool operator< (const IntVect& p) const;

  ///
  /**
     Returns true if this IntVect is less than or equal to argument
     IntVect.  All comparisons between analogous components must be
     satisfied.  Note that, since the comparison is component-wise, it
     is possible for an IntVect to be neither greater than or equal
     to, less than or equal to, nor equal to another.
  */
  bool operator<= (const IntVect& p) const;

  ///
  /**
     Returns true if this IntVect is greater than argument IntVect.
     All comparisons between analogous components must be satisfied.
     Note that, since the comparison is component-wise, it is possible
     for an IntVect to be neither greater than, less than, nor equal
     to another.
  */
  bool operator> (const IntVect& p) const;

  ///
  /**
     Returns true if this IntVect is greater than or equal to argument
     IntVect.  All comparisons between analogous components must be
     satisfied.  Note that, since the comparison is component-wise, it
     is possible for an IntVect to be neither greater than or equal
     to, less than or equal to, nor equal to another.
  */

  bool operator>= (const IntVect& p) const;

  ///
  /**
     Returns true if this IntVect is lexically less than the argument.
     An IntVect MUST BE either lexically less than, lexically greater
     than, or equal to another IntVect.

     iv1 is lexically less than iv2 if:

     in 2-D:<br>
     (iv1[0] < iv2[0]) || ((iv1[0] == iv2[0]) && (iv1[1] < iv2[1]));

     in 3-D:<br>
     (iv1[0] < iv2[0]) || (iv1[0]==iv2[0] && ((iv1[1] < iv2[1] || ((iv1[1] == iv2[1]) && (iv1[2] < iv2[2])))));
  */
  bool lexLT (const IntVect& s) const;

  ///
  /**
     Returns true if this IntVect is lexically greater than the
     argument.  An IntVect MUST BE either lexically less than,
     lexically greater than, or equal to another IntVect.

     iv1 is lexically less than iv2 if:

     in 2-D:<br>
     (iv1[0] > iv2[0]) || ((iv1[0] == iv2[0]) && (iv1[1] > iv2[1]));

     in 3-D:<br>
     (iv1[0] > iv2[0]) || (iv1[0]==iv2[0] && ((iv1[1] > iv2[1] || ((iv1[1] == iv2[1]) && (iv1[2] > iv2[2])))));
  */
  bool lexGT (const IntVect& s) const;

  /*@}*/

  /**
     \name Unary operators
  */
  /*@{*/

  ///
  /**
     Unary plus -- for completeness.
  */
  IntVect operator+ () const;

  ///
  /**
     Unary minus -- negates all components of this IntVect.
  */
  IntVect operator- () const;

  ///
  /**
     Sum of all components of this IntVect.
  */
  int sum () const;

  ///
  /**
     Product of all components of this IntVect.
  */
  int product () const;

  /*@}*/

  /**
     \name Addition operators
  */
  /*@{*/

  ///
  /**
     Modifies this IntVect by addition of a scalar to each component.
  */
  IntVect& operator+= (int s);

  ///
  /**
     Modifies this IntVect by component-wise addition with argument.
  */
  IntVect& operator+= (const IntVect& p);

  ///
  /**
     Returns component-wise sum of this IntVect and argument.
  */
  IntVect operator+ (const IntVect& p) const;

  ///
  /**
     Return an IntVect that is this IntVect with a scalar added to
     each component.
  */
  IntVect operator+ (int s) const;

  ///
  /**
     Returns an IntVect that is an IntVect <i>p</i>
     with a scalar <i>s</i> added to each component.
  */
  friend inline IntVect operator+ (int            s,
                                   const IntVect& p);

  /*@}*/

  /**
     \name Subtraction operators
  */
  /*@{*/

  ///
  /**
     Modifies this IntVect by subtraction of a scalar from each
     component.
  */
  IntVect& operator-= (int s);

  ///
  /**
     Modifies this IntVect by component-wise subtraction by argument.
  */
  IntVect& operator-= (const IntVect& p);

  ///
  /**
     Returns an IntVect that is this IntVect with <i>p</i> subtracted
     from it component-wise.
  */
  IntVect operator- (const IntVect& p) const;

  ///
  /**
     Returns an IntVect that is this IntVect with a scalar <i>s</i> subtracted
     from each component.
  */
  IntVect operator- (int s) const;

  ///
  /**
     Returns <i>s - p</i>.
  */
  friend inline IntVect operator- (int            s,
                                   const IntVect& p);

  /*@}*/

  /**
     \name Multiplication operators
  */
  /*@{*/

  ///
  /**
     Modifies this IntVect by multiplication of each component by
     a scalar.
  */
  IntVect& operator*= (int s);

  ///
  /**
     Modifies this IntVect by component-wise multiplication by argument.
  */
  IntVect& operator*= (const IntVect& p);

  ///
  /**
     Returns component-wise product of this IntVect with argument.
  */
  IntVect operator* (const IntVect& p) const;

  ///
  /**
     Returns an IntVect that is this IntVect with each component
     multiplied by a scalar.
  */
  IntVect operator* (int s) const;

  ///
  /**
     Returns an IntVect that is an IntVect <i>p</i> with each component
     multiplied by a scalar <i>s</i>.
  */
  friend inline IntVect operator* (int            s,
                                   const IntVect& p);

  /*@}*/

  /**
     \name Division operators
  */
  /*@{*/

  ///
  /**
     Modifies this IntVect by division of each component by a scalar.
  */
  IntVect& operator/= (int s);

  ///
  /**
     Modifies this IntVect by component-wise division by IntVect
     argument.
  */
  IntVect& operator/= (const IntVect& p);

  ///
  /**
     Returns component-wise quotient of this IntVect by argument.
  */
  IntVect operator/ (const IntVect& p) const;

  ///
  /**
     Returns an IntVect that is this IntVect with each component
     divided by a scalar.
  */
  IntVect operator/ (int s) const;

  /*@}*/

  /**
     \name Other arithmetic operators
  */
  /*@{*/

  ///
  /**
     Modifies this IntVect by taking component-wise min with IntVect
     argument.
  */
  IntVect& min (const IntVect& p);

  ///
  /**
     Returns the IntVect that is the component-wise minimum of two
     argument IntVects.
  */
  friend inline IntVect min (const IntVect& p1,
                             const IntVect& p2);

  ///
  /**
     Modifies this IntVect by taking component-wise max with IntVect
     argument.
  */
  IntVect& max (const IntVect& p);

  ///
  /**
     Returns the IntVect that is the component-wise maximum of two
     argument IntVects.
  */
  friend inline IntVect max (const IntVect& p1,
                             const IntVect& p2);

  ///
  /**
     Modifies this IntVect by multiplying each component by a scalar.
  */
  IntVect& scale (int s);

  ///
  /**
     Returns an IntVect obtained by multiplying each of the components
     of the given IntVect by a scalar.
  */
  friend inline IntVect scale (const IntVect& p,
                               int            s);

  /// Returns the componentwise absolute value of the given IntVect.
  friend const IntVect absolute(const IntVect& p);

  ///
  /**
     Modifies IntVect by reflecting it in the plane defined by the
     index <i>ref_ix</i> and with normal in the direction of <i>idir</i>.
     Directions are based at zero.
  */
  IntVect& reflect (int ref_ix,
                    int idir);

  ///
  /**
     Returns an IntVect that is the reflection of the given IntVect in
     the plane which passes through <i>ref_ix</i> and normal to the
     coordinate direction idir.
  */
  friend inline IntVect reflect(const IntVect& a,
                                int            ref_ix,
                                int            idir);

  ///
  /**
     Modifies this IntVect by adding <i>s</i> to component in given coordinate
     direction.
  */
  IntVect& shift (int coord,
                  int s);

  ///
  /**
     Modifies this IntVect by component-wise addition with IntVect
     argument.
  */
  IntVect& shift (const IntVect& iv);

  ///
  /**
     Modifies this IntVect by adding a scalar <i>s</i> to each component.

  */
  IntVect& diagShift (int s);

  ///
  /**
     Returns IntVect obtained by adding a scalar to each of the
     components of the given IntVect.
  */
  friend inline IntVect diagShift (const IntVect& p,
                                   int            s);

  ///
  /**
     Modify IntVect by component-wise integer projection.
  */
  IntVect& coarsen (const IntVect& p);

  ///
  /**
     Modify IntVect by component-wise integer projection.
  */
  IntVect& coarsen (int p);

  ///
  /**
     Returns an IntVect that is the component-wise integer projection of
     <i>p</i> by <i>s</i>.
  */
  friend inline IntVect coarsen (const IntVect& p,
                                 int            s);

  ///
  /**
     Returns an IntVect which is the component-wise integer projection
     of IntVect <i>p1</i> by IntVect <i>p2</i>.
  */
  friend inline IntVect coarsen (const IntVect& p1,
                                 const IntVect& p2);

  /*@}*/

  /**
     \name I/O Functions
  */
  /*@{*/

  ///
  /**
     Print an IntVect to the ostream.
  */
  void printOn (std::ostream& os) const;

  ///
  /**
     Print an IntVect to the pout().
  */
  void p() const;

  ///
  /**
     Print an IntVect to the ostream a bit more verbosely.
  */
  void dumpOn (std::ostream& os) const;

  ///
  /**
     Print the IntVect to given output stream in ASCII.
  */
  friend std::ostream& operator<< (std::ostream&       os,
                                   const IntVect& iv);

  ///
  /**
     Read next IntVect from given input stream.
  */
  friend std::istream& operator>> (std::istream& os,
                                   IntVect& iv);

  /*@}*/

  /**
     \name Constants
  */
  /*@{*/

  ///
  /**
     Returns a basis vector in the given coordinate direction.<br>
     In 3-D:
     BASISV(0) == (1,0,0); BASISV(1) == (0,1,0); BASISV(2) == (0,0,1).<br>
     In 2-D:
     BASISV(0) == (1,0); BASISV(1) == (0,1).<br>
     Note that the coordinate directions are based at zero.
  */
  friend inline IntVect BASISV (int dir);

  /**
     This is an IntVect all of whose components are equal to zero.
  */
  static const IntVect Zero;

  /**
     This is an IntVect all of whose components are equal to one.
  */
  static const IntVect Unit;

  /**
     Initializes Zero and Unit.
  */
  static int InitStatics();

  /*@}*/

protected:
  //
  // Box is a friend of ours.
  //
  friend class Box;

  friend class HDF5Handle;
  friend class VolIndex;
  friend class FaceIndex;
  friend class BLfacade::IntVect;

  /**
     The individual components of this IntVect.
  */
  int vect[SpaceDim];

  /**
     Number of bytes of storage used by this IntVect.
  */
  static const size_t IntVectSize;
};

//
// Static initialization.  Gotta make sure there are no static object
// definitions above here (except possibly stuff in headers).  Otherwise,
// the danger is that some static object's constructor calls IntVect::Zero or
// IntVect::Unit -- the very things the following definition is supposed to
// initialize.
//
static int s_dummyForIntVectH = IntVect::InitStatics();

#ifndef WRAPPER
//
// Inlines.
//

// try uninitialized IntVect null construction for now.....

// inline
// IntVect::IntVect ()
// {
//     D_EXPR6(vect[0] = 0, vect[1] = 0, vect[2] = 0);
// }

inline
IntVect::IntVect (D_DECL6(int i, int j, int k,
                          int l, int m, int n))
{
  D_EXPR6(vect[0] = i, vect[1] = j, vect[2] = k, vect[3] = l, vect[4] = m, vect[5] = n);
}

inline
IntVect::IntVect (const int *a)
{
  D_EXPR6(vect[0] = a[0], vect[1] = a[1], vect[2] = a[2],
          vect[3] = a[3], vect[4] = a[4], vect[5] = a[5]);
}

inline
IntVect::IntVect (const IntVect &iv)
{
  //D_EXPR6(vect[0]=iv.vect[0], vect[1]=iv.vect[1], vect[2]=iv.vect[2]);
  memcpy(vect, iv.vect, IntVectSize);
}

inline
IntVect&
IntVect::operator= (const IntVect &iv)
{
  D_EXPR6(vect[0]=iv.vect[0], vect[1]=iv.vect[1], vect[2]=iv.vect[2],
          vect[3]=iv.vect[3], vect[4]=iv.vect[4], vect[5]=iv.vect[5]);
  return *this;
}

inline
int&
IntVect::operator[] (int i)
{
  CH_assert(i>=0 && i < SpaceDim);
  return vect[i];
}

inline
int
IntVect::operator[] (int i) const
{
  CH_assert(i>=0 && i < SpaceDim);
  return vect[i];
}

inline
void
IntVect::setVal (int i,
                 int val)
{
  CH_assert(i >=0 && i < SpaceDim);
  vect[i] = val;
  //  return *this;
}

inline
const int*
IntVect::dataPtr() const
{
  return vect;
}

inline
int*
IntVect::dataPtr()
{
  return vect;
}

inline
const int*
IntVect::getVect () const
{
  return vect;
}

inline
bool
IntVect::operator== (const IntVect& p) const
{
  return D_TERM6(vect[0] == p[0], && vect[1] == p[1], && vect[2] == p[2],
                 && vect[3] == p[3], && vect[4] == p[4], && vect[5] == p[5]);
}

inline
bool
IntVect::operator!= (const IntVect& p) const
{
  return D_TERM6(vect[0] != p[0], || vect[1] != p[1], || vect[2] != p[2],
                 || vect[3] != p[3], || vect[4] != p[4], || vect[5] != p[5]);
}

inline
bool
IntVect::operator< (const IntVect& p) const
{
  return D_TERM6(vect[0] < p[0], && vect[1] < p[1], && vect[2] < p[2],
                 && vect[3] < p[3], && vect[4] < p[4], && vect[5] < p[5]);
}

inline
bool
IntVect::operator<= (const IntVect& p) const
{
  return D_TERM6(vect[0] <= p[0], && vect[1] <= p[1], && vect[2] <= p[2],
                 && vect[3] <= p[3], && vect[4] <= p[4], && vect[5] <= p[5]);
}

inline
bool
IntVect::operator> (const IntVect& p) const
{
  return D_TERM6(vect[0] > p[0], && vect[1] > p[1], && vect[2] > p[2],
                 && vect[3] > p[3], && vect[4] > p[4], && vect[5] > p[5]);
}

inline
bool
IntVect::operator>= (const IntVect& p) const
{
  return D_TERM6(vect[0] >= p[0], && vect[1] >= p[1], && vect[2] >= p[2],
                 && vect[3] >= p[3], && vect[4] >= p[4], && vect[5] >= p[5]);
}

inline
bool
IntVect::lexLT (const IntVect &s) const

{
  if (vect[0] < s[0]) return true;
#if CH_SPACEDIM > 1
  if (vect[0] > s[0]) return false;
  if (vect[1] < s[1]) return true;
#endif
#if CH_SPACEDIM > 2
  if (vect[1] > s[1]) return false;
  if (vect[2] < s[2]) return true;
#endif
#if CH_SPACEDIM > 3
  if (vect[2] > s[2]) return false;
  if (vect[3] < s[3]) return true;
#endif
#if CH_SPACEDIM > 4
  if (vect[3] > s[3]) return false;
  if (vect[4] < s[4]) return true;
#endif
#if CH_SPACEDIM > 5
  if (vect[4] > s[4]) return false;
  if (vect[5] < s[5]) return true;
#endif

  return false;
}

inline
bool
IntVect::lexGT (const IntVect& s) const
{
  if (vect[0] > s[0]) return true;
#if CH_SPACEDIM > 1
  if (vect[0] < s[0]) return false;
  if (vect[1] > s[1]) return true;
#endif
#if CH_SPACEDIM > 2
  if (vect[1] < s[1]) return false;
  if (vect[2] > s[2]) return true;
#endif
#if CH_SPACEDIM > 3
  if (vect[2] < s[2]) return false;
  if (vect[3] > s[3]) return true;
#endif
#if CH_SPACEDIM > 4
  if (vect[3] < s[3]) return false;
  if (vect[4] > s[4]) return true;
#endif
#if CH_SPACEDIM > 5
  if (vect[4] < s[4]) return false;
  if (vect[5] > s[5]) return true;
#endif

  return false;
}

inline
IntVect
IntVect::operator+ () const
{
  return IntVect(*this);
}

inline
IntVect
IntVect::operator- () const
{
  return IntVect(D_DECL6(-vect[0], -vect[1], -vect[2],
                         -vect[3], -vect[4], -vect[5] ));
}

inline
int
IntVect::sum () const
{
  return D_TERM6(vect[0], + vect[1], + vect[2],
                 + vect[3], + vect[4], + vect[5]);
}

inline
int
IntVect::product () const
{
  return D_TERM6(vect[0], * vect[1], * vect[2],
                 * vect[3], * vect[4], * vect[5]);
}

inline
IntVect&
IntVect::operator+= (int s)
{
  D_EXPR6(vect[0] += s, vect[1] += s, vect[2] += s,
          vect[3] += s, vect[4] += s, vect[5] += s);
  return *this;
}

inline
IntVect&
IntVect::operator+= (const IntVect& p)
{
  D_EXPR6(vect[0] += p[0], vect[1] += p[1], vect[2] += p[2],
          vect[3] += p[3], vect[4] += p[4], vect[5] += p[5]);
  return *this;
}

inline
IntVect&
IntVect::operator*= (int s)
{
  D_EXPR6(vect[0] *= s, vect[1] *= s, vect[2] *= s,
          vect[3] *= s, vect[4] *= s, vect[5] *= s);
  return *this;
}

inline
IntVect&
IntVect::operator*= (const IntVect &p)
{
  D_EXPR6(vect[0] *= p[0], vect[1] *= p[1], vect[2] *= p[2],
          vect[3] *= p[3], vect[4] *= p[4], vect[5] *= p[5]);
  return *this;
}

inline
IntVect&
IntVect::operator/= (int s)
{
  D_EXPR6(vect[0] /= s, vect[1] /= s, vect[2] /= s,
          vect[3] /= s, vect[4] /= s, vect[5] /= s);
  return *this;
}

inline
IntVect&
IntVect::operator/= (const IntVect& p)
{
  D_EXPR6(vect[0] /= p[0], vect[1] /= p[1], vect[2] /= p[2],
          vect[3] /= p[3], vect[4] /= p[4], vect[5] /= p[5]);
  return *this;
}

inline
IntVect&
IntVect::operator-= (int s)
{
  D_EXPR6(vect[0] -= s, vect[1] -= s, vect[2] -= s,
          vect[3] -= s, vect[4] -= s, vect[5] -= s);
  return *this;
}

inline
IntVect&
IntVect::operator-= (const IntVect& p)
{
  D_EXPR6(vect[0] -= p[0], vect[1] -= p[1], vect[2] -= p[2],
          vect[3] -= p[3], vect[4] -= p[4], vect[5] -= p[5]);
  return *this;
}

inline
IntVect
IntVect::operator+ (const IntVect& p) const
{
  return IntVect(D_DECL6(vect[0] + p[0], vect[1] + p[1], vect[2] + p[2],
                         vect[3] + p[3], vect[4] + p[4], vect[5] + p[5]));
}

inline
IntVect
IntVect::operator+ (int s) const
{
  return IntVect(D_DECL6(vect[0] + s, vect[1] + s, vect[2] + s,
                         vect[3] + s, vect[4] + s, vect[5] + s));
}

inline
IntVect
IntVect::operator- (const IntVect& p) const
{
  return IntVect(D_DECL6(vect[0] - p[0], vect[1] - p[1], vect[2] - p[2],
                         vect[3] - p[3], vect[4] - p[4], vect[5] - p[5]));
}

inline
IntVect
IntVect::operator- (int s) const
{
  return IntVect(D_DECL6(vect[0] - s, vect[1] - s, vect[2] - s,
                         vect[3] - s, vect[4] - s, vect[5] - s));
}

inline
IntVect
IntVect::operator* (const IntVect& p) const
{
  return IntVect(D_DECL6(vect[0] * p[0], vect[1] * p[1], vect[2] * p[2],
                         vect[3] * p[3], vect[4] * p[4], vect[5] * p[5]));
}

inline
IntVect
IntVect::operator* (int s) const
{
  return IntVect(D_DECL6(vect[0] * s, vect[1] * s, vect[2] * s,
                         vect[3] * s, vect[4] * s, vect[5] * s));
}

inline
IntVect
IntVect::operator/ (const IntVect& p) const
{
  return IntVect(D_DECL6(vect[0] / p[0], vect[1] / p[1], vect[2] / p[2],
                         vect[3] / p[3], vect[4] / p[4], vect[5] / p[5]));
}

inline
IntVect
IntVect::operator/ (int s) const
{
  return IntVect(D_DECL6(vect[0] / s, vect[1] / s, vect[2] / s,
                         vect[3] / s, vect[4] / s, vect[5] / s));
}

inline
IntVect&
IntVect::min (const IntVect& p)
{
  D_EXPR6(vect[0] = Min(vect[0], p.vect[0]),
          vect[1] = Min(vect[1], p.vect[1]),
          vect[2] = Min(vect[2], p.vect[2]),
          vect[3] = Min(vect[3], p.vect[3]),
          vect[4] = Min(vect[4], p.vect[4]),
          vect[5] = Min(vect[5], p.vect[5]));
  return *this;
}

inline
IntVect&
IntVect::max (const IntVect& p)
{
  D_EXPR6(vect[0] = Max(vect[0], p.vect[0]),
          vect[1] = Max(vect[1], p.vect[1]),
          vect[2] = Max(vect[2], p.vect[2]),
          vect[3] = Max(vect[3], p.vect[3]),
          vect[4] = Max(vect[4], p.vect[4]),
          vect[5] = Max(vect[5], p.vect[5]));
  return *this;
}

inline
IntVect&
IntVect::scale (int s)
{
  D_EXPR6(vect[0] *= s, vect[1] *= s, vect[2] *= s,
          vect[3] *= s, vect[4] *= s, vect[5] *= s);
  return *this;
}

inline
IntVect&
IntVect::reflect (int ref_ix,
                  int idir)
{
  CH_assert(idir >= 0 && idir < SpaceDim);
  vect[idir] = -vect[idir] + 2*ref_ix;
  return *this;
}

inline
IntVect&
IntVect::shift (int coord,
                int s)
{
  CH_assert(coord >= 0 && coord < SpaceDim);
  vect[coord] += s;
  return *this;
}

inline
IntVect&
IntVect::shift (const IntVect& iv)
{
  *this += iv;
  return *this;
}

inline
IntVect&
IntVect::diagShift (int s)
{
  D_EXPR6(vect[0] += s, vect[1] += s, vect[2] += s,
          vect[3] += s, vect[4] += s, vect[5] += s);
  return *this;
}

inline
IntVect
operator+ (int            s,
           const IntVect& p)
{
  return IntVect(D_DECL6(p[0] + s, p[1] + s, p[2] + s,
                         p[3] + s, p[4] + s, p[5] + s));
}

inline
IntVect
operator- (int            s,
           const IntVect& p)
{
  return IntVect(D_DECL6(s - p[0], s - p[1], s - p[2],
                         s - p[3], s - p[4], s - p[5]));
}

inline
IntVect
operator* (int            s,
           const IntVect& p)
{
  return IntVect(D_DECL6(s * p[0], s * p[1], s * p[2],
                         s * p[3], s * p[4], s * p[5]));
}

inline
IntVect
scale (const IntVect& p,
       int            s)
{
  return IntVect(D_DECL6(s * p[0], s * p[1], s * p[2],
                         s * p[3], s * p[4], s * p[5]));
}

inline
const IntVect
absolute (const IntVect& p)
{
  return IntVect(D_DECL6(abs(p[0]), abs(p[1]), abs(p[2]),
                         abs(p[3]), abs(p[4]), abs(p[5])));
}

inline
IntVect
diagShift (const IntVect &p, int s)
{
  return IntVect(D_DECL6(p[0] + s, p[1] + s, p[2] + s,
                         p[3] + s, p[4] + s, p[5] + s));
}

inline
IntVect
min (const IntVect& p1,
     const IntVect& p2)
{
  IntVect p(p1);
  return p.min(p2);
}

inline
IntVect
max (const IntVect& p1,
     const IntVect& p2)
{
  IntVect p(p1);
  return p.max(p2);
}

inline
IntVect
BASISV (int dir)
{
  CH_assert(dir >= 0 && dir < SpaceDim);
  IntVect tmp = IntVect::Zero ;
  tmp.vect[dir] = 1;
  return tmp;
}

inline
IntVect
reflect (const IntVect& a,
         int            ref_ix,
         int            idir)
{
  CH_assert(idir >= 0 && idir < SpaceDim);
  IntVect b(a);
  b.vect[idir] = -b.vect[idir] + 2*ref_ix;
  return b;
}

inline
IntVect
coarsen (const IntVect& p,
         int            s)
{
  CH_assert(s > 0);
  return IntVect(
                 D_DECL6((p.vect[0]<0) ? -abs(p.vect[0]+1)/s-1 : p.vect[0]/s ,
                         (p.vect[1]<0) ? -abs(p.vect[1]+1)/s-1 : p.vect[1]/s ,
                         (p.vect[2]<0) ? -abs(p.vect[2]+1)/s-1 : p.vect[2]/s ,
                         (p.vect[3]<0) ? -abs(p.vect[3]+1)/s-1 : p.vect[3]/s ,
                         (p.vect[4]<0) ? -abs(p.vect[4]+1)/s-1 : p.vect[4]/s ,
                         (p.vect[5]<0) ? -abs(p.vect[5]+1)/s-1 : p.vect[5]/s ));
}

inline
IntVect
coarsen (const IntVect& p1,
         const IntVect& p2)
{
  CH_assert(p2 > IntVect::Zero);
  return IntVect(
                 D_DECL6(
                        (p1.vect[0]<0)?-abs(p1.vect[0]+1)/p2.vect[0]-1:p1.vect[0]/p2.vect[0],
                        (p1.vect[1]<0)?-abs(p1.vect[1]+1)/p2.vect[1]-1:p1.vect[1]/p2.vect[1],
                        (p1.vect[2]<0)?-abs(p1.vect[2]+1)/p2.vect[2]-1:p1.vect[2]/p2.vect[2],
                        (p1.vect[3]<0)?-abs(p1.vect[3]+1)/p2.vect[3]-1:p1.vect[3]/p2.vect[3],
                        (p1.vect[4]<0)?-abs(p1.vect[4]+1)/p2.vect[4]-1:p1.vect[4]/p2.vect[4],
                        (p1.vect[5]<0)?-abs(p1.vect[5]+1)/p2.vect[5]-1:p1.vect[5]/p2.vect[5])
                 );
}

inline
IntVect&
IntVect::coarsen (int s)
{
  CH_assert(s > 0);
  for (int i = 0; i < SpaceDim; ++i)
    vect[i] = ((vect[i]<0) ? -abs(vect[i]+1)/s-1 : vect[i]/s);
  return *this;
}

inline
IntVect&
IntVect::coarsen (const IntVect& p)
{
  CH_assert(p > IntVect::Zero);
  for (int i = 0; i <SpaceDim; ++i)
    {
      const int s = p.vect[i];
      vect[i] = ((vect[i]<0) ? -abs(vect[i]+1)/s-1 : vect[i]/s);
    }
  return *this;
}

#endif /* WRAPPER */

#include "NamespaceFooter.H"
#endif
