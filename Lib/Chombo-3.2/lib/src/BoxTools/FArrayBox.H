#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FARRAYBOX_H_
#define _FARRAYBOX_H_

#ifndef WRAPPER
#include <iostream>

// #include <Pointers.H>
#include "Box.H"
#include "BaseFab.H"
#endif

#include "REAL.H"
#include "SPACE.H"
#include "NamespaceHeader.H"

///
/**
   Fortran Array Boxes (generally called FABs) are objects constructed
   to interface with arrays in Fortran.  Useful operations can be performed
   upon FABs in C++, and they provide a convenient interface to
   Fortran when it is necessary to retreat into that language for
   doing arithmetic operations when performance matters.

   FArrayBox is derived from BaseFab<Real>.
   FArrayBox adds additional useful capabilities which make sense
   for Real types, such as I/O and L**p norms.

   The C pre-processor macro `CH_SPACEDIM' must be defined to use
   this class.  The internal precision of FArrayBox objects is
   set by defining either `CH_USE_FLOAT' or `CH_USE_DOUBLE'

   This class does NOT provide a copy constructor or assignment operator.
*/
class FArrayBox: public BaseFab<Real>
{
public:

  ///
  /**
     Constructs an invalid FArrayBox with no memory.
  */
  FArrayBox ();

  ///
  /**
     Constructs an initial FArrayBox with the data space allocated but not
     inititialized. a_ncomp is the number of components (variables) at each
     data point in the Box.
  */
  FArrayBox(const Box& a_box,
            int        a_ncomp,
            Real*      a_alias = NULL);

  ///
  /**
     Construct an aliaed FArrayBox.  See BaseFab class for details.
  */
  FArrayBox(const Interval& a_comps,
            FArrayBox&      a_original)
    :
    BaseFab<Real>(a_comps, a_original)
  {}

  ///
  /**
     Defines FArrayBox with the data space allocated but not
     inititialized. a_ncomp is the number of components (variables) at each
     data point in the Box.
  */
  virtual void define(const Box& a_box,
                      int        a_ncomp,
                      Real*      a_alias = NULL)
  {
    resize(a_box, a_ncomp, a_alias);
  }

  ///
  /**
     The (virtual) destructor.
  */
  virtual ~FArrayBox ();

  ///
  /**
     Constructs an 'aliased' BaseFab of the requested interval of the
     argument BaseFab.  This BaseFab does not allocate any memory, but
     sets its data pointer into the memory pointed to by the argument
     BaseFab.  It is the users responsiblity to ensure this aliased
     BaseFab is not used after the original BaseFab has deleted its data ptr
     (resize, define(..) called, or destruction, etc.).

     This aliased BaseFab will also generate side effects (modifying the values
     of data in one will modify the other's data).

     This aliased BaseFab will have a_comps.size() components, starting at zero.
  */
  virtual void define(const Interval& a_comps,
                      FArrayBox&      a_original)
  {
    BaseFab<Real>::define(a_comps, a_original);
  }

  ///
  /**
   * This is here only to make the Intel compiler stop warning about
   * partial override.
   */
  virtual void define(const Interval& a_comps,
              BaseFab<Real>&      a_original)
  {
    BaseFab<Real>::define(a_comps, a_original);
  }

  ///
  /**
     Constructs FArrayBox by reading it from istream.
  */
  explicit FArrayBox(std::istream& a_is);

  ///
  /**
     Returns the Lp-norm of this FAB using components
     (a_comp : a_comp+a_numcomp-1) and within the a_subbox.
     a_p < 0 -> ERROR
     a_p = 0 -> infinity norm (max norm)
     a_p = 1 -> sum of ABS(FAB)
     a_p > 1 -> Lp-norm
  */
  virtual Real norm(const Box& a_subbox,
                    int        a_p = 2,
                    int        a_comp = 0,
                    int        a_numcomp = 1) const;

  ///
  /**
     Returns the Lp-norm of this FAB using components
     (a_comp : a_comp+a_numcomp-1).
     a_p < 0  -> ERROR
     a_p = 0  -> infinity norm (max norm)
     a_p = 1  -> sum of ABS(FAB)
     a_p > 1  -> Lp-norm
  */
  virtual Real norm(int a_p = 2,
                    int a_comp = 0,
                    int a_numcomp = 1) const;

  ///
  /**
     Returns sum of pow(fab[i,c],p): i in a_subbox, a_comp <= c <
     a_comp+a_numcomp, a_p >= 2 only
  */
  virtual Real sumPow(const Box& a_subbox,
                      int        a_p = 2,
                      int        a_comp = 0,
                      int        a_numcomp = 1) const;

  /// Return the dot product of this FArrayBox with another
  /**
     Return the dot product of this FArrayBox and "a_fab2" over their common
     box and all components.
  */
  Real dotProduct(const FArrayBox& a_fab2) const;

  /// Return the dot product of this FArrayBox with another
  /**
     Return the dot product of this FArrayBox and "a_fab2" over the
     a_box box and all components.
  */
  Real dotProduct(const FArrayBox& a_fab2, const Box& a_box) const;

  ///
  /**
     Returns the minimum value of given component of this FArrayBox.
  */
  Real min(int a_comp = 0) const;

  ///
  /**
     Returns the minimum value of given component of this FArrayBox in
     given a_subbox.

  */
  Real min(const Box& a_subbox,
           int        a_comp = 0) const;

  ///
  /**
     Returns the maximum value of given component of this FArrayBox.
  */
  Real max(int a_comp = 0) const;

  ///
  /**
     Returns the maximum value of given component of this FArrayBox in
     given a_subbox.

  */
  Real max(const Box& a_subbox,
           int        a_comp = 0) const;

  ///
  /**
     Finds location of minimum value in given component of this FArrayBox.
  */
  IntVect minIndex(int a_comp = 0) const;

  ///
  /**
     Returns location of minimum value in given component  of this FArrayBox
     in given a_subbox.
  */
  IntVect minIndex(const Box& a_subbox,
                   int        a_comp = 0) const;

  ///
  /**
     Returns location of maximum value in given component  of this FArrayBox.
  */
  IntVect maxIndex(int a_comp = 0) const;

  ///
  /**
     Returns location of maximum value in given component of this FArrayBox
     in given a_subbox.
  */
  IntVect maxIndex(const Box& a_subbox,
                   int        a_comp = 0) const;

  ///
  /**
     Computes a_mask array with value of 1 in cells where this FArrayBox
     has value less than a_val, 0 otherwise.  a_mask is resized by this
     function.  The number of cells marked with 1 returned.
  */
  int maskLT(BaseFab<int>& a_mask,
             Real          a_val,
             int           a_comp = 0) const;

  ///
  /**
     Computes a_mask array with value of 1 in cells where this FArrayBox
     has value less than or equal to a_val, 0 otherwise.  a_mask is
     resized by this function.  The number of cells marked with 1
     returned.
  */
  int maskLE(BaseFab<int>& a_mask,
             Real          a_val,
             int           a_comp = 0) const;
  ///
  /**
     Computes a_mask array with value of 1 in cells where this FArrayBox
     has value equal to a_val, 0 otherwise.  a_mask is resized by this
     function.  The number of cells marked with 1 returned.

  */
  int maskEQ(BaseFab<int>& a_mask,
             Real          a_val,
             int           a_comp = 0) const;
  ///
  /**
     Computes a_mask array with value of 1 in cells where this FArrayBox
     has value greater than a_val, 0 otherwise.  a_mask is resized by this
     function.  The number of cells marked with 1 returned.
  */
  int maskGT(BaseFab<int>& a_mask,
             Real          a_val,
             int           a_comp = 0) const;

  ///
  /**
     Computes a_mask array with value of 1 in cells where this FArrayBox
     has value greater than or equal to a_val, 0 otherwise.  a_mask is
     resized by this function.  The number of cells marked with 1 returned.
  */
  int maskGE(BaseFab<int>& a_mask,
             Real          a_val,
             int           a_comp = 0) const;

  ///
  /**
     Modifies this FArrayBox by replacing each value with its absolute value.
  */
  void abs();

  ///
  /**
     Modifies this FArrayBox by replacing each value with its absolute value,
     for components (a_comp : a_comp+a_numcomp-1).
  */
  void abs(int a_comp,
           int a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by replacing eahc value with its absolute value,
     for components (a_comp : a_comp+a_numcomp-1) and within the a_subbox.
  */
  void abs (const Box& a_subbox,
            int        a_comp = 0,
            int        a_numcomp = 1);

  ///
  /**
     Returns sum of given component of FArrayBox.
  */
  Real sum(int a_comp,
           int a_numcomp = 1) const;

  ///
  /**
     Returns sum of component of this FArrayBox in given a_subbox.
  */
  Real sum(const Box& a_subbox,
           int        a_comp,
           int        a_numcomp = 1) const;

  ///
  /**
     Modifies this FArrayBox by replacing each value x with a_r/x.
  */
  FArrayBox& invert(Real a_r);

  ///
  /**
     Modifies this FArrayBox by replacing each value x with a_r/x.  For
     given range of components.
  */
  FArrayBox& invert(Real a_r,
                    int  a_comp,
                    int  a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by replacing each value x with a_r/x.  For
     given range of components and within given a_subbox.
  */
  FArrayBox& invert(Real       a_r,
                    const Box& a_subbox,
                    int        a_comp = 0,
                    int        a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by replacing each value with its additive
     inverse.  For given range of components and within given a_subbox.
  */
  FArrayBox& negate(const Box& a_subbox,
                    int        a_comp = 0,
                    int        a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by replacing each value with its additive
     inverse.  For given range of components.
  */
  FArrayBox& negate(int a_comp,
                    int a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by replacing each value with its additive
     inverse.
  */
  FArrayBox& negate();

  ///
  /**
     Modifies this FArrayBox by adding the scalar Real a_r to all values.  For
     given range of components and within given a_subbox.
  */
  FArrayBox& plus(Real       a_r,
                  const Box& a_subbox,
                  int        a_comp = 0,
                  int        a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by adding the scalar Real a_r to all values.  For
     given range of components.
  */
  FArrayBox& plus(Real a_r,
                  int  a_comp,
                  int  a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by adding the scalar Real a_r to all values.
  */
  FArrayBox& operator += (Real a_r);

  ///
  /**
     Modifies this FArrayBox by pointwise addition of the values of the
     argument FArrayBox.  You might come to grief if the domains of the
     FArrayBoxes don't match.
  */
  FArrayBox& operator += (const FArrayBox& a_x);

  ///
  /**
     Modifies this FArrayBox by adding the scalar Real a_r to all values.
  */
  FArrayBox& plus(Real a_r);

  FArrayBox& plus_real(Real a_r)
  {
    return this->plus(a_r);
  }

  ///
  /**
     Modifies this FArrayBox by pointwise addition of the values of the
     argument FArrayBox.  You might come to grief if the domains of the
     FArrayBoxes don't match.  The same as += operator.

  */
  FArrayBox& plus(const FArrayBox& a_x);

  ///
  /**
     Modifies this FArrayBox by pointwise scaled addition of the
     argument FArrayBox (a[i] <- a[i] + a_scale * a_src[i]).  Uses domain
     of the intersection of these two FArrayBoxes.
  */
  FArrayBox& plus(const FArrayBox& a_src,
                  const Real&      a_scale);

  ///
  /**
     Modifies this FArrayBox by pointwise scaled addition of the
     argument FArrayBox (a[i] <- a[i] + a_scale * a_src[i]).  Uses domain
     of the intersection of these two FArrayBoxes.
  */
  FArrayBox& plus(const FArrayBox& a_src,
                  const Real&      a_scale,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by pointwise addition of values in the argument
     FArrayBox.  Adds src's components (a_srccomp : a_srccomp+a_numcomp-1)
     to this FArrayBox's components (a_destcomp : a_destcomp+a_numcomp-1)
     where the domains of the two FArrayBoxes intersect.
  */
  FArrayBox& plus(const FArrayBox& a_src,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by pointwise addition of values in the argument
     FArrayBox.  Adds src's components (a_srccomp : a_srccomp+a_numcomp-1)
     to this FArrayBox's components (a_destcomp : a_destcomp+numcomp-1)
     where the domain of this FArrayBox intersects the a_subbox.  NOTE:
     a_subbox must be contained in this FAB.
  */
  FArrayBox& plus(const FArrayBox& a_src,
                  const Box&       a_subbox,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by pointwise addition of values in the argument
     FArrayBox.  Adds src's components (a_srccomp : a_srccomp+a_numcomp-1)
     in the Box a_srcbox to this FArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) in the Box a_destbox.  Corresponding locations
     within the two FArrayBoxes are indexed relative to a_srcbox and a_destbox,
     and will in general not be the same.  The a_srcbox and a_destbox must be
     same size.  The results are UNDEFINED if the a_src and dest FArrayBoxes
     are the same and the a_srcbox and a_destbox overlap.

  */
  FArrayBox& plus(const FArrayBox& a_src,
                  const Box&       a_srcbox,
                  const Box&       a_destbox,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  FArrayBox& plus(const FArrayBox& a_src,
                  const Box&       a_srcbox,
                  const Box&       a_destbox,
                  const Real&      a_scale,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by subtracting the scalar Real a_r to all values.
     Note: use plus(-a_r) for more general operations.
  */
  FArrayBox& operator -= (Real a_r);

  ///
  /**
     Modifies this FArrayBox by pointwise subtraction of the values of the
     argument FArrayBox.  You might come to grief if the domains of the
     FArrayBoxes don't match.
  */
  FArrayBox& operator -= (const FArrayBox& a_x);

  ///
  /**
     Modifies this FArrayBox by pointwise subtraction of the values of the
     argument FArrayBox.  You might come to grief if the domains of the
     FArrayBoxes don't match.  The same as -= operator.
  */
  FArrayBox& minus(const FArrayBox& a_x);

  /**
     Modifies this FArrayBox by pointwise subtraction of values in the
     argument FArrayBox.  Subtracts a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) from this FArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) where the domains of the two FArrayBoxes
     intersect.
  */
  FArrayBox& minus(const FArrayBox& a_src,
                   int              a_srccomp,
                   int              a_destcomp,
                   int              a_numcomp = 1);

  /**
     Modifies this FArrayBox by pointwise subtraction of values in the
     argument FArrayBox.  Subtracts a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) from this FArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) where the domain of this FArrayBox intersects
     the a_subbox.  NOTE: a_subbox must be contained in this FAB.
  */
  FArrayBox& minus(const FArrayBox& a_src,
                   const Box&       a_subbox,
                   int              a_srccomp,
                   int              a_destcomp,
                   int              a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by pointwise subtraction of values in the
     argument FArrayBox.  Subtracts a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) in the Box a_srcbox from this FArrayBox's
     components (a_destcomp : a_destcomp+a_numcomp-1) in the Box a_destbox.
     Corresponding locations within the two FArrayBoxes are indexed relative
     to a_srcbox and a_destbox, and will in general not be the same.  The
     a_srcbox and a_destbox must be same size.  The results are UNDEFINED
     if the a_src and dest FArrayBoxes are the same and the a_srcbox and
     a_destbox overlap.
  */
  FArrayBox& minus(const FArrayBox& a_src,
                   const Box&       a_srcbox,
                   const Box&       a_destbox,
                   int              a_srccomp,
                   int              a_destcomp,
                   int              a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by multiplying all values by the scalar Real a_r.
  */
  FArrayBox& operator *= (Real a_r);

  ///
  /**
     Modifies this FArrayBox by multiplying all values by the scalar Real a_r.
  */
  FArrayBox& mult(Real a_r);

  ///
  /**
     Modifies this FArrayBox by multiplying all values by the scalar
     Real a_r.  For given range of components.
  */
  FArrayBox& mult(Real a_r,
                  int  a_comp,
                  int  a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by multiplying all values by the scalar
     Real a_r.  For given range of components and within given a_subbox.
  */
  FArrayBox& mult(Real       a_r,
                  const Box& a_subbox,
                  int        a_comp = 0,
                  int        a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by pointwise multiplication of the values by the
     argument FArrayBox.  You might come to grief if the domains of the
     FArrayBoxes don't match.
  */
  FArrayBox& operator *= (const FArrayBox& a_x);

  ///
  /**
     Modifies this FArrayBox by pointwise multiplication by the values in the
     argument FArrayBox.  You might come to grief if the domains of the
     FArrayBoxes don't match.  The same as the *= operator.
  */
  FArrayBox& mult(const FArrayBox& a_x);

  ///
  /**
     Modifies this FArrayBox by pointwise multiplication by values in the
     argument FArrayBox.  Multiplies a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) by this FArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) where the domains of the two FArrayBoxes
     intersect.
  */
  FArrayBox& mult(const FArrayBox& a_src,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by pointwise multiplication by values in the
     argument FArrayBox.  Multiplies a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) by this FArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) where the domain of this FArrayBox intersects
     the a_subbox.  NOTE: a_subbox must be contained in this FAB.
  */
  FArrayBox& mult(const FArrayBox& a_src,
                  const Box&       a_subbox,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by pointwise multiplication by values in the
     argument FArrayBox.  Multiplies a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) in the Box a_srcbox by this FArrayBox's
     components (a_destcomp : a_destcomp+a_numcomp-1) in the Box a_destbox.
     Corresponding locations within the two FArrayBoxes are indexed relative
     to a_srcbox and a_destbox, and will in general not be the same.  The
     a_srcbox and a_destbox must be same size.  The results are UNDEFINED if
     the a_src and dest FArrayBoxes are the same and the a_srcbox and a_destbox
     overlap.
  */
  FArrayBox& mult(const FArrayBox& a_src,
                  const Box&       a_srcbox,
                  const Box&       a_destbox,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by dividing all values by the scalar Real a_r.
  */
  FArrayBox& operator /= (Real a_r);

  ///
  /**
     Modifies this FArrayBox by dividing all values by the scalar Real a_r.
  */
  FArrayBox& divide(Real a_r);

  ///
  /**
     Modifies this FArrayBox by dividing all values by the scalar Real a_r.
     For given range of components.
  */
  FArrayBox& divide(Real a_r,
                    int  a_comp,
                    int  a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by dividing all values by the scalar Real
     a_r.  For given range of components and within given a_subbox.
  */
  FArrayBox& divide(Real       a_r,
                    const Box& a_subbox,
                    int        a_comp = 0,
                    int        a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by pointwise division of the values by the
     argument FArrayBox.  You might come to grief if the domains of the
     FArrayBoxes don't match.
  */
  FArrayBox& operator /= (const FArrayBox& a_x);

  ///
  /**
     Modifies this FArrayBox by pointwise division by the values in the
     argument FArrayBox.  You might come to grief if the domains of the
     FArrayBoxes don't match.  The same as the /= operator.
  */
  FArrayBox& divide(const FArrayBox& a_x);

  ///
  /**
     Modifies this FArrayBox by pointwise division by values in the argument
     FArrayBox.  Divides this FArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) by a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) where the domains of the two FArrayBoxes intersect.
  */
  FArrayBox& divide(const FArrayBox& a_src,
                    int              a_srccomp,
                    int              a_destcomp,
                    int              a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by pointwise division by values in the argument
     FArrayBox.  Divides this FArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) by a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) where the domain of this FArrayBox intersects
     the a_subbox.  NOTE: a_subbox must be contained in this FAB.
  */
  FArrayBox& divide(const FArrayBox& a_src,
                    const Box&       a_subbox,
                    int              a_srccomp,
                    int              a_destcomp,
                    int              a_numcomp = 1);

  ///
  /**
     Modifies this FArrayBox by pointwise division by values in the argument
     FArrayBox.  Divides this FArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) in the Box a_destbox by a_src's components
     (a_srccomp : a_srccomp+a_numcomp-1) in the Box a_srcbox.  Corresponding
     locations within the two FArrayBoxes are indexed relative to a_srcbox and
     a_destbox, and will in general not be the same.  The a_srcbox and
     a_destbox must be same size.  The results are UNDEFINED if the a_src and
     dest FArrayBoxes are the same and the a_srcbox and a_destbox overlap.
  */
  FArrayBox& divide(const FArrayBox& a_src,
                    const Box&       a_srcbox,
                    const Box&       a_destbox,
                    int              a_srccomp,
                    int              a_destcomp,
                    int              a_numcomp = 1);

  ///
  Real get(const IntVect& a_iv,
           int            a_comp) const
  {
    return this->operator()(a_iv, a_comp);
  }

  ///
  void set(const IntVect& a_iv,
           int            a_comp,
           Real           a_val)
  {
    this->operator()(a_iv, a_comp) = a_val;
  }

  //! Computes a_A * a_X + a_B * a_Y, placing the result in this FArrayBox.
  FArrayBox& axby(const FArrayBox& a_X, const FArrayBox& a_Y,
                  Real a_A, Real a_B);

  //! Computes a_A * a_X + a_B * a_Y, placing the result in this FArrayBox.
  //! This version performs this operation only for the given component
  //! in each FArrayBox.
  FArrayBox& axby(const FArrayBox& a_X, const FArrayBox& a_Y,
                  Real a_A, Real a_B,
                  int a_destComp, int a_xComp, int a_yComp);

protected:
  virtual void performCopy(const BaseFab<Real>& a_src,
                           const Box&           a_srcbox,
                           int                  a_srccomp,
                           const Box&           a_destbox,
                           int                  a_destcomp,
                           int                  a_numcomp);


private:
  //
  // These are disallowed.
  //
  FArrayBox (const FArrayBox&);
  FArrayBox& operator = (const FArrayBox&);
};

#include "NamespaceFooter.H"
#endif
