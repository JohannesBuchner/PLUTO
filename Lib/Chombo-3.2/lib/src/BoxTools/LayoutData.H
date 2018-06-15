#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LAYOUTDATA_H_
#define _LAYOUTDATA_H_

#include "BoxLayout.H"
#include "DataIterator.H"
#include "TimedDataIterator.H"
#include "NamespaceHeader.H"

//class LayoutData

/// Data that maintains a one-to-one mapping of T to the boxes in a BoxLayout
/**
   A collection of Box-oriented objects.  The arrangement of the T objects
   is given by the underlying BoxLayout object.  LayoutData does not attempt
   to prevent users from manipulating the template class T to become out of
   synch with the boxes in the BoxLayout.  Caveat emptor.

   Non-local (off-processor access) to a LayoutData index is an error.  You
   can assure that access to the data stored in a LayoutData is local to your
   processor by using the DataIterator object for indexing.

   Data in a LayoutData CANNOT be communicated to other processors using
   the API presented in LevelData.

   class T must provide the following methods:
   <PRE>
   {
     T()
   }
   </PRE>

   This class is useful for processor-local data that needs to be indexable
   along with a BoxLayoutData.  Auxillary data objects, etc.
*/
template<class T> class LayoutData
{
public:
  ///
  LayoutData();

  ///
  /**
     constructor.  allocates a T object for every box in the BoxLayout a_dp
     using the T null constructor.  a_dp must be closed.
  */
  LayoutData(const BoxLayout& a_dp);

  ///
  virtual ~LayoutData();

  ///
  /**
     constructor.  allocates a T object for every box in the BoxLayout a_dp
     using the T null constructor.  a_dp must be closed.  previously stored
     T objects are removed
  */
  virtual void define(const BoxLayout& a_dp);

  // this approach, of using the iterator itself for indexing is my own
  // brilliance.  This way you don't need all the crappy [Const][Dependent]
  // Iterators.  It also forces all data access to be done through the
  // DataIterator, which will know how to talk to KeLP (or other parallel
  // architectures, like multithreading, etc).

  ///
  DataIterator dataIterator() const;

  TimedDataIterator timedDataIterator() const;

  /// const accessor function
  const T& operator [] (const DataIndex& a_index) const;

  /// const accessor function
  /**
     Returns member for which DataIterator is currently indexed to.
     Equivalent to layout[a_iterator()]
  */
  const T& operator [] (const DataIterator& a_iterator) const;

  /// non-const accessor function
  T& operator [] (const DataIndex& a_index);

  /// non-const accessor function
  T& operator [] (const DataIterator& a_iterator);

  ///
  Box box(const DataIndex& a_index) const;

  Box box(const DataIterator& a_iterator) const;

  // not really a protection mechanism, since a user can always
  // perform a shallow copy of the BoxLayout, and be free to manipulate
  // the data in the object.  If you wish to have an actual copy
  // of the BoxLayout, then you need to invoke the clone() method.

  ///
  const BoxLayout& boxLayout() const
  {
    return m_boxLayout;
  }


protected:

  // for BoxLayout::define(const LayoutData<Box>& a_newLayout)
  friend class BoxLayout;

  BoxLayout m_boxLayout;

  // this used to be std::vector<T>, and I could let vector handle
  // destruction for me, but vector.resize() absolutely demands that
  // I provide a copy constructor for T.  People get uncomfortable when
  // I make them provide me with copy constructors.
  Vector<T*> m_vector;

  // thinking about making this class self-documenting to a greater degree
  // and having the component names also be kept in the class and carried
  // around through various transformations.....
  //  vector<string> m_componentNames;

  bool m_callDelete;


private:
  // disallow copy constructors and assignment operators
  // to avoid very hard-to-find performance problems
  LayoutData<T>& operator= (const LayoutData<T>& a_rhs);
  LayoutData(const LayoutData& a_rhs);

  // handy usage function, tied to the implementation of this class
  // as a vector<T*>.  Assumes that m_comps and m_boxLayout have already
  // been initialized correctly.  Sets the correct size for the data
  // vector, deletes the extra T objects if m_vector is to shorten, and
  // performs either construction or define on the remaining T objects.
  void allocate();
};

//====================== inline functions ===========================
// not literally a .H file, but just an experiment in having useable,
// readable headers, while having the dependeny system work properly.
// Since LayoutData.H now depends on LayoutDataI.H, then changing either
// should cause all code that includes LayoutData.H to be recompiled.
// This way people can just read the interface in this file.
// Implementation is put in the *I.H file.

#include "NamespaceFooter.H"
#include "LayoutDataI.H"

#endif
