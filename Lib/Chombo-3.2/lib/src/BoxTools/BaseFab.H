#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEFAB_H_
#define _BASEFAB_H_

#include <cstdlib>
#include "SPACE.H"
#include <string>
#include <typeinfo>
#include <cmath>

#include "Box.H"
#include "Arena.H"
#include "Interval.H"
#include "REAL.H"
#include "NamespaceHeader.H"

struct SliceSpec;

/** A Fortran Array-like Container

  BaseFab emulates the Fortran array concept.
  Useful operations can be performed upon
  BaseFab's in C++, and they provide a convenient interface to
  Fortran when it is necessary to retreat into that language.

  `BaseFab' is a template class.  Through use of the
  template, a `BaseFab' may be based upon any class.  So far at least,
  most applications have been based upon simple types like `integer's,
  `real's, or `doubleprecision's.  Most applications do not use BaseFab's
  directly, but utilize specialized classes derived from BaseFab.

  It will be easier to use a `BaseFab' if you understand the following
  concepts.  `BaseFab' objects depend on the dimensionality of space
  (indirectly through the DOMAIN `Box' member).  It is
  typical to define the macro `CH_SPACEDIM' to be 1, 2, or 3 to indicate
  the dimension of space.  See the discussion of class `Box' for more
  information.  A `BaseFab' contains a `Box' DOMAIN, which indicates the
  integer indexing space over which the array is defined.  A `BaseFab'
  also has NVAR components.  By components, we mean that for each
  point in the rectangular indexing space, there are NVAR values
  associated with that point.  A Fortran array corresponding to a
  `BaseFab' would have (CH_SPACEDIM+1) dimensions.

  By design, the array layout in a `BaseFab' mirrors that of a
  Fortran array.  The first index (x direction for example) varies
  most rapidly, the next index (y direction), if any, varies next
  fastest. The component index varies last, after all the spatial
  indices.

  It is sometimes convenient to be able to treat a sub-array within an
  existing `BaseFab' as a `BaseFab' in its own right.  This is often
  referred to as 'aliasing' the `BaseFab'.  Note that when aliasing is
  used, the BaseFab's domain will not, in general, be the same as the
  parent BaseFab's domain, nor will the number of components.
  BaseFab is a dimension dependent class, so CH_SPACEDIM must be
  defined as either 1, 2, or 3 when compiling.

  This is NOT a polymorphic class.

  It does NOT provide a copy constructor or assignment operator.

  T MUST have a default constructor and an assignment operator.
*/

extern Real BaseFabRealSetVal;

template <class T> class BaseFab
{
public:
  /// {\bf constructors, destructor and defines}

  ///for AggStencil and its minions
  virtual long offset(const IntVect& a_iv, const int& a_ivar) const
  {
    CH_assert(a_ivar >= 0);
    CH_assert(a_ivar < m_nvar);
    CH_assert(!(m_dptr == 0));
    CH_assert(m_domain.contains(a_iv));

    IntVect ivDiff = a_iv - m_domain.smallEnd();
    IntVect ivSize = m_domain.size();
    long offset = ivDiff[0];

#if CH_SPACEDIM > 1
    offset += ivDiff[1]*ivSize[0] ;
#endif
#if CH_SPACEDIM > 2
    offset += ivDiff[2]*ivSize[0]*ivSize[1];
#endif
#if CH_SPACEDIM > 3
    offset += ivDiff[3]*ivSize[0]*ivSize[1]*ivSize[2];
#endif
#if CH_SPACEDIM > 4
    offset += ivDiff[4]*ivSize[0]*ivSize[1]*ivSize[2]*ivSize[3];
#endif
#if CH_SPACEDIM > 5
    offset += ivDiff[5]*ivSize[0]*ivSize[1]*ivSize[2]*ivSize[3]*ivSize[4];
#endif
    offset += m_numpts * a_ivar;

    return offset;
  }

  /**
     Constructs an invalid `BaseFab'.  The domain is invalid, the
     number of components is zero, and no actual array memory is
     allocated.  An invalid `BaseFab' must be resize()d
     (see `BaseFab::resize') before use.
  */
  BaseFab ();

  ///
  /**
     Constructs a BaseFab with desired domain and number of components.
  */
  BaseFab(const Box& a_bx,
          int        a_n,
          T*         a_alias = NULL);

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
  BaseFab(const Interval& a_comps,
          BaseFab<T>&     a_original);

  ///
  /**
     The destructor deletes the array memory. Unless this was an aliased
     BaseFab.
  */
  virtual ~BaseFab();

  /**
     This function resizes a `BaseFab' so it covers the `Box' a_b with a_n
     components.  If a_alias is not NULL the memory it points to is used
     (without checking size); otherwise the existing data is lost and the
     new data is uninitialized.
  */
  void resize(const Box& a_b,
              int        a_n = 1,
              T*         a_alias=NULL);

  /**
     Make BaseFab with desired domain and number of components.  Existing
     data is lost.  Data is in uninialized state.
  */
  virtual void define(const Box& a_box,
                      int        a_comps,
                      T*         a_alias = NULL)
  {
    resize(a_box, a_comps, a_alias);
  }

  /**
     alias define.  no memory allocated.  this BaseFab sets its data ptr
     directly in the a_original BaseFab data.
  */
  virtual void define(const Interval& a_comps,
                      BaseFab<T>&     a_original);

  /**
     The function returns the `BaseFab' to the invalid state.  (See
     comments for constructors above.)  The memory is freed.
  */
  void clear();

  /// {\bf accessors}

  ///
  /**
     Returns the number of components.
  */
  int nComp() const;

  ///
  /**
     Returns the domain (box) where the array is defined.
  */
  const Box& box() const;

  /**
     Returns an IntVect giving the length
     of the domain in each direction.
    */
  IntVect size() const;

  /**
     Returns the lower corner of the domain.  See class `Box' for analogue.
    */
  const IntVect& smallEnd() const;

  /**
     Returns the upper corner of the domain.  See class `Box' for analogue.
    */
  const IntVect& bigEnd() const;

  ///
  /**
     Returns an Interval for the entire range on components.
  */
  Interval interval() const
  {
    return Interval(0, m_nvar-1);
  }

  /**
     Returns a modifiable lvalue reference to the Nth component value defined
     at position p in the domain.  This operator may be inefficient if the
     C++ compiler is unable to optimize the C++ code.
  */
  T& operator () (const IntVect& a_p,
                  int            a_N);

  T& operator () (const IntVect& a_p);

  /**
     Returns a conatant reference to the Nth component value defined at
     position p in the domain.  This operator may be inefficient if the C++
     compiler is unable to optimize the C++ code.
  */
  const T& operator () (const IntVect& p,
                        int            N) const;

  const T& operator () (const IntVect& p) const;

  /**
     This function puts numcomp component values, starting at component N,
     from position pos in the domain into array data, that must be allocated
     by the user.
  */
  void getVal(T*             a_data,
              const IntVect& a_pos,
              int            a_N,
              int            a_numcomp) const;

  /**
     This function puts all component values, starting at component 0,
     from position pos in the domain into array data, that must be allocated
     by the user.
  */
  void getVal(T*             a_data,
              const IntVect& a_pos) const;

  ///{\bf Fortran interface functions}

  /**
     Returns the lower corner of the domain.  Instead of returning them in
     the form of IntVects, as in smallEnd and bigEnd, it returns the values
     as a pointer to an array of constant integers.  This is useful when
     interfacing to Fortran subroutines.  It should not be used in any other
     context!!!
  */
  const int* loVect() const;

  /**
     Returns the upper corner of the domain.  Instead of returning them in
     the form of IntVects, as in smallEnd and bigEnd, it returns the values
     as a pointer to an array of constant integers.  This is useful when
     interfacing to Fortran subroutines.  It should not be used in any other
     context!!!
  */
  const int* hiVect() const;

  /**
     Returns a pointer to an integer that contains the number of components
     in the BaseFab. This is useful when interfacing to Fortran subroutines.
     It should not be used in any other context!!!
  */
  const int* nCompPtr() const;

  /**
     Returns a pointer to an object of type T that is the value of the a_nth
     component associated with the cell at the low end of the domain.  This
     is commonly used to get a pointer to data in the array which is then
     handed off to a Fortran subroutine. It should not be used in any other
     context!!!  Remember that data is stored in Fortran array order, with
     the component index coming last.  In other words, `dataPtr' returns a
     pointer to all the a_nth components.
  */
  T* dataPtr(int a_n = 0);

  /**
     Returns a constant pointer to an object of type T that is the value of
     the a_nth component associated with the cell at the low end of the domain.
     This is commonly used to get a pointer to data in the array which is
     then handed off to a Fortran subroutine.  It should not be used in any
     other context!!! Remember that data is stored in Fortran array order,
     with the component index coming last.  In other words, `dataPtr' returns
     a pointer to all the a_nth components.
  */
  const T* dataPtr(int a_n = 0) const;

  ///{\bf comparison functions}

  /**
     Returns true if the domain of a_fab is totally contained within the domain
     of this `BaseFab'.
    */
  bool contains(const BaseFab<T>& a_fab) const;

  /**
     Returns true if a_bx is totally contained within the domain of this
     `BaseFab'.
  */
  bool contains(const Box& a_bx) const;

  ///{\bf data modification functions}

  /**
     The setVal functions set subregions in the `BaseFab' to a constant value.
     This most general form specifies the subbox, the starting component
     number, and the number of components to be set.
  */
  void setVal(T          a_x,
              const Box& a_bx,
              int        a_nstart,
              int        a_numcomp);

  ///
  /**
     Modifies this BaseFab so that all values of a component, a_n, in the
     specified Box, a_bx, are set to the given value, a_x.
  */
  void setVal(T          a_x,
              const Box& a_bx,
              int        a_n);

  ///
  /**
     Modifies this BaseFab so that all values of a component, a_n, are set to
     the given value, a_x.
  */
  void setVal(T   a_x,
              int a_n);

  ///
  /**
     Modifies this BaseFab so that all values of all components are set to
     the given value, a_x.
  */
  void setVal(T a_x);

  /**
     Modifies this BaseFab by copying the contents of the argument BaseFab
     into it.  This, the most general form of copy, specifies the contents
     of any sub-box a_srcbox in `BaseFab' a_src may be copied into a (possibly
     different) a_destbox in the destination `BaseFab'.  Note that although
     the a_srcbox and the a_destbox may be disjoint, they must be the same size
     and shape.  If the sizes differ, the copy is undefined and a runtime
     error results.  This copy function is the only one of the copy
     functions to allow a copy between differing boxes. The user also
     specifies how many components are copied, starting at component
     a_srccomp in a_src and stored starting at component a_destcomp.  The
     results are UNDEFINED if the a_src and dest BaseFabs are the same and
     the a_srcbox and a_destbox overlap.

  */
  BaseFab<T>& copy(const BaseFab<T>& a_src,
                   const Box&        a_srcbox,
                   int               a_srccomp,
                   const Box&        a_destbox,
                   int               a_destcomp,
                   int               a_numcomp);

  /**
     Modifies this BaseFab by copying the contents of the argument BaseFab
     into it.  A copy within the intersecting region of the domains of the
     two BaseFabs is performed.  The user specifies how many components are
     copied, starting at component a_srccomp in a_src and stored starting at
     component a_destcomp.
  */
  BaseFab<T>& copy(const BaseFab<T>& a_src,
                   int               a_srccomp,
                   int               a_destcomp,
                   int               a_numcomp = 1);

  /**
     Modifies this BaseFab by copying the contents of the argument BaseFab
     into it.  A copy within the intersecting region of the domains of the
     two BaseFabs and the specified Box a_destbox is performed.  All
     components are copied.
  */
  BaseFab<T>& copy(const BaseFab<T>& a_src,
                   const Box&        a_destbox);

  /**
     Modifies this BaseFab by copying the contents of the argument BaseFab
     into it.  A copy within the intersecting region of the domains of the
     two BaseFabs is performed.  All components are copied.
  */
  BaseFab<T>& copy(const BaseFab<T>& a_src);

  ///
  /**
     Copy from a subsection of one box into another.  Assumes the boxes are
     both in the same index space, and that box R is completely contained
     in both the source and destination boxes.
    */
  void copy(const Box&        a_RegionFrom,
            const Interval&   a_Cdest,
            const Box&        a_RegionTo,
            const BaseFab<T>& a_src,
            const Interval&   a_Csrc);

  ///    {\bf domain modification functions}

  ///
  /**
     Modifies the domain of this BaseFab by shifting.  Equivalent to
     fab.shift(0,a_v[0]).shift(1,a_v[1])...  There is no effect upon
     the array memory.
  */
  BaseFab<T>& shift(const IntVect& a_v);

  ///
  /**
     Modifies the domain of this BaseFab by shifting it a_ncells
     indexing positions in coordinate direction a_idir.  Directions are
     zero-based.  It is an error if not 0 <= a_idir < SpaceDim.  There is
     no effect upon the array memory.
  */
  BaseFab<T>& shift(int a_idir,
                    int a_ncells);

  ///
  /**
     Modifies the domain of this BaseFab by shifting by "half" indices,
     thereby converting the Box from type CELL to NODE or vice-versa.
     fab.shiftHalf(0,1) shifts the domain to the right by 1/2 cells.
     fab.shiftHalf(1,-3) shifts the domain in the -j direction by 3/2 cells.
     NOTE: If a_numHalfs is EVEN the shift is a_numHalfs/2 full zones
     and hence will not change the type.  This is: fab.shiftHalf(1,4) ==
     fab.shift(1,2).  Directions are zero-based.  It is an error if not 0
     <= a_dir < SpaceDim.   There is no effect upon the array memory.
  */
  BaseFab<T>& shiftHalf(int a_dir,
                        int a_numHalfs);

  ///
  /**
     Modifies the domain of this BaseFab by shifting by half indices.
     Equivalent to fab.shiftHalf(0,a_v[0]).shiftHalf(1,a_v[1]) ...
     There is no effect upon the array memory.
  */
  BaseFab<T>& shiftHalf(const IntVect& a_v);

  ///{\bf linearization functions}

  ///
  /**
     Returns the size, in number of bytes, of a flat linear
     representation of the data in this object in the area defined
     by the input Box a_box and the component Interval a_comps.  The size
     does not include the size of a_box and a_comps.
  */
  virtual int size(const Box&      a_box,
                   const Interval& a_comps) const;

  ///
  /**
     Write a linear representation of the internal data.  Assumes that
     sufficient memory for the buffer has already been allocated by the
     caller.
  */
  virtual void linearOut(void*           a_buf,
                         const Box&      a_R,
                         const Interval& a_comps) const;

  /// Same as linearOut, but returns the current location in the buffer
  virtual void* linearOut2(void*           a_buf,
                           const Box&      a_R,
                           const Interval& a_comps) const;

  ///
  virtual void linearIn(void*           a_buf,
                        const Box&      a_R,
                        const Interval& a_comps);

  /// same as linearIn, but returns the current location in the buffer
  virtual void* linearIn2(void*           a_buf,
                          const Box&      a_R,
                          const Interval& a_comps);

  ///
  static int preAllocatable()
  {
    return 0; // static preAllocatable
  }

  /**
    Turns a_slice into a BaseFab that's the same as *this except that it's just
    one cell thick in the a_sliceSpec.direction-th direction, and its
    coordinate in that direction is a_sliceSpec.position.

    If a_sliceSpec.position is outside the range of *this, that's a fatal error.
  */
  void degenerate(BaseFab<T>& a_slice, const SliceSpec& a_sliceSpec) const;

  ///
  bool isAliased() const;

  ///regression test
  static int test();

  ///regression test
  static int testBoxAndComp();
protected:
  //
  // Allocates memory for the `BaseFab<T>'.
  //
  void define();

  //
  // Deallocates memory for the `BaseFab<T>'.
  //
  void undefine();

  static std::string name();

  //
  // The function called by BaseFab copy operations.
  //
  virtual void performCopy(const BaseFab<T>& a_src,
                           const Box&        a_srcbox,
                           int               a_srccomp,
                           const Box&        a_destbox,
                           int               a_destcomp,
                           int               a_numcomp);

  //
  // This function is called by the `BaseFab' setVal operations.
  //
  void performSetVal(T          a_x,
                     const Box& a_bx,
                     int        a_nstart,
                     int        a_numcomp);

  //
  // template class static data member.  not sure if this will
  // work with all compilers.  It has been in the draft standard
  // for some time
  //
  static Arena* s_Arena;

  Box  m_domain;   // My index space.
  int  m_nvar;     // Number components.
  long m_numpts;   // Cached number of points in FAB.
  long m_truesize; // m_nvar * m_numpts that was allocated on heap
                   //   (only if m_aliased == false).
  T*   m_dptr;     // The data pointer.
  bool m_aliased;  // The BaseFab is not allocated memory, but is an alias. bvs

private:
  //
  // These functions are made private to prevent use of the default
  // functions provided by the C++ compiler.
  //
  BaseFab<T>& operator= (const BaseFab<T>&);
  BaseFab (const BaseFab<T>&);
};

#include "NamespaceFooter.H"

#include "BaseFabImplem.H"

#endif
