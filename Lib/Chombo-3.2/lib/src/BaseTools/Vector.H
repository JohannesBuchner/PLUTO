#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <cmath>
#include <cstdlib>
#include "CH_assert.H"
#include <vector>
#include <typeinfo>
#include <string>
#include <map>
#include <list>
#include <algorithm>
#include <iostream>

#include "BaseNamespaceHeader.H"

#ifdef CH_USE_MEMORY_TRACKING
//these objects are used for memory tracking system
typedef std::pair<long int*, long int*> ppair;
typedef std::map<std::string, ppair > VectorList;

typedef std::list<int> lst; //not currently used, bvs.
typedef std::map<std::string, lst*> incr; // not currently used, bvs

extern VectorList* vectorList_;
extern incr*       vectorIncr_;
#endif

using std::vector;
using std::ostream;

/// one dimensional dynamic array
/**
   Vector is a resizable one-dimensional array with constant-time
   random access and range checking.  The template type T must have a
   default constructor, a copy constructor, and an assignment
   operator.

*/
template <class T>

class Vector
{
  //friend class ParticleVector;
public:
  ///
  /**
      Default constructor.
      Creates a Vector of zero length with null data.
  */
  inline
  Vector()
  {
  }

  ///
  /**
      Destructor.
  */
  inline
  virtual ~Vector()
  {
#ifdef CH_USE_MEMORY_TRACKING
    if (size()>0) decrement(size());
#endif
  }

  ///
  /**
     Copy constructor.
  */
  inline
  Vector(const Vector<T>& invec):v(invec.v)
  {
#ifdef CH_USE_MEMORY_TRACKING
    increment(invec.size());
#endif
  }

  ///
  /**  conversion constructor
  */
  inline
  Vector(const std::vector<T>& invec):v(invec)
  {
#ifdef CH_USE_MEMORY_TRACKING
    increment(invec.size());
#endif
  }

  ///
  /**
  */
  inline
  Vector<T>& operator=(const std::vector<T>& invec)
  {
#ifdef CH_USE_MEMORY_TRACKING
    decrement(size());
#endif
    v=invec;
#ifdef CH_USE_MEMORY_TRACKING
    increment(size());
#endif
    return *this;
  }

  ///
  /**
  */
  inline
  Vector<T>& operator=(const Vector<T>& invec)
  {
#ifdef CH_USE_MEMORY_TRACKING
    decrement(size());
#endif
    v=invec.v;
#ifdef CH_USE_MEMORY_TRACKING
    increment(size());
#endif
    return *this;
  }

  /// assign a scalar to every element of the vector
  /** [NOTE: cant use operator=() for this because it would be
   *         ambiguous with the (int) constructor for number T types.]
  */
  inline
  Vector<T>& assign(const T& inval)
  {
    v.assign( size(),inval );
    return *this;
  }

  ///
  /**
      Constructs a Vector with given number of elements.\\

      {\bf Arguments:}\\
      size (not modified):  number of elements of Vector to construct.\\
      {\bf This:}\\
      -------The object is modified----------

  */

  inline
  Vector(/*long*/ unsigned int isize)
  {
    resize(isize);
  }

  ///
  /**
  */
  inline
  void clear()
  {
#ifdef CH_USE_MEMORY_TRACKING
    decrement(size());
#endif
    v.clear();
  }

  ///
  /**  size function. returns current size of Vector
  */
  inline
  size_t size() const
  {
    return v.size();
  }

  ///
  /**
      Constructs a Vector with given number of elements and constant value.\\

      {\bf Arguments:}\\
      size (not modified):  number of elements of Vector to construct.\\
      value (not modified):  value to set every element to.\\
      {\bf This:}\\
      -------The object is modified----------

  */
  inline
  Vector(/*long*/ unsigned int isize, const T& value) : v(isize, value)
  {
#ifdef CH_USE_MEMORY_TRACKING
    increment(isize);
#endif
  }

  ///
  /**
      Returns a modifiable lvalue reference to the value of the given
      element in this Vector.  It is an error if n < 0 or n >= this->size(). \\

      {\bf Arguments:}\\
      n (not modified) index of desired element.\\
      {\bf Returns:}\\
      modifiable reference to value in Vector at index n.\\
      {\bf This:}\\
      -----
      This object is modified if the returned reference is assigned a new value
      -----
  */
  inline
  T& operator[] (/*long*/ unsigned int n)
  {
    CH_assert(n < size());
    //vector<T>* svec = static_cast<vector<T>* >(this);
    //return  svec->operator[](n);
    //    return  (static_cast<vector<T>* >(this))->operator[](n);
    return v[n];
  }

  ///
  /**
      Returns a constant reference to the given element in this
      Vector.\\

      {\bf Arguments:}\\
      n (not modified) index of desired element.\\
      {\bf Returns:}\\
      constant reference to value in Vector at index n.\\
      {\bf This:}\\
      This object is not modified.
  */
  inline
  const T& operator[] (/*long*/ unsigned int n) const
  {
    CH_assert(n < size());
    //const vector<T>* svec = static_cast<const vector<T>* >(this);
    //return  svec->operator[](n);
    //   return  (static_cast<const vector<T>* >(this))->operator[](n);
    return v[n];
  }

  inline
  T &front()
  {
    return v.front();
  }

  inline
  T &back()
  {
    return v.back();
  }

  inline
  void pop_back()
  {
    v.pop_back();
  }

  inline
  const T& back() const
  {
    return v.back();
  }

  inline
  void swap(Vector<T>& other)
  {
    v.swap(other.v);
  }
  ///
  /**
  */
  inline
  void push_back(const T& in)
  {
#ifdef CH_USE_MEMORY_TRACKING
    increment(1);
#endif
    //  (static_cast<vector<T>* >(this))->push_back(in);
    v.push_back(in);
  }

  ///
  /**
     Modifies this Vector by appending the elements of the argument
     Vector.  The new Vector will have a size of this->size() +
     invec.size() (where this Vector is considered before the append is
     performed).  The first element of invec will have index
     this->size(), the second element will have index this->size()+1,
     etc.\\

     {\bf Arguments:}\\
     invec (not modified):  Vector whose elements to append to this Vector.\\
     {\bf This:}\\
     -------The object is modified----------
  */
  inline
  void  append(const Vector<T>& invec)
  {
#ifdef CH_USE_MEMORY_TRACKING
    increment(invec.size());
#endif
    for (int ivec = 0; ivec <  invec.size(); ivec++)
      {
        // (static_cast<vector<T>* >(this))->push_back(invec[ivec]);
        v.push_back(invec[ivec]);
      }
  }

  ///
  /**
  */
  inline
  void resize(/*long*/ unsigned int isize)
  {
    //    vector<T>* svec = static_cast<vector<T>* >(this);
#ifdef CH_USE_MEMORY_TRACKING
    if (isize > size()) increment(isize-size());
    else               decrement(size()-isize);
#endif
    unsigned int l= size();
    //    svec->resize(isize);
    v.resize(isize);
    for (; l<isize; ++l)
      // svec->operator[](l) = T();
      v[l] = T();
  }

  inline
  void reserve(/*long*/ size_t isize)
  {
    v.reserve(isize);
  }

  inline
  size_t capacity() const
  {
    return v.capacity();
  }
  ///
  /**
  */
  inline
  void resize(/*long*/ unsigned int isize, const T& value)
  {
    //   vector<T>* svec = static_cast<vector<T>* >(this);
#ifdef CH_USE_MEMORY_TRACKING
    if (isize > size()) increment(isize-size());
    else               decrement(size()-isize);
#endif
    //    svec->resize(isize, value);
    v.resize(isize, value);
  }

  ///
  /**
  */
  inline void sort()
  {
    std::sort(v.begin(), v.end());
  }

  /// Returns std::vector under the hood.
  inline std::vector<T>& stdVector()
  {
    return v;
  }

  inline const std::vector<T>& constStdVector() const
  {
    return v;
  }

#ifdef CH_USE_MEMORY_TRACKING
  // static template memory tracking data
  static long int bytes;
  static long int peak;
  //  static lst*  incrementList;
#endif

private:

#ifdef CH_USE_MEMORY_TRACKING
  inline
  void increment(unsigned int i)
  {
    // warning: pointless comparison of unsigned integer with zero
    // CH_assert(i>=0);  // commented out (ndk)
    static unsigned int sizzle = initFunc();
    i*=sizzle;
    bytes += i;
    //  incrementList->insert(incrementList->begin(), i);
    if (bytes > peak) peak = bytes;
  }

  inline
  void decrement(unsigned int i)
  {
    if (i==0) return;
    // CH_assert(bytes>=0);
    i*=sizeof(T);
    bytes-=i;
    //  incrementList->insert(incrementList->begin(), -i);
    //  CH_assert(bytes>=0);
  }
#endif

  unsigned int initFunc();
  std::vector<T> v;
};

#ifdef CH_USE_MEMORY_TRACKING
template <class T>
long int Vector<T>::bytes=0;

template <class T>
long int Vector<T>::peak=0;

//template <class T>
//lst*  Vector<T>::incrementList=NULL;
#endif

//class ostream;

template <class T>
ostream& operator<<(ostream& os, const Vector<T>& vec)
{
  for (int i=0; i<vec.size(); i++) os<<vec[i]<<" ";
  return os;
}

#ifdef CH_USE_MEMORY_TRACKING

template <class T>
unsigned int  Vector<T>::initFunc()
{
  if (vectorList_ == NULL)
  {
    vectorList_ = new VectorList;
    vectorIncr_ = new incr;
  }
  //  if (incrementList == NULL) incrementList = new lst;
  ppair tmp(&bytes, &peak);
  vectorList_->operator[](typeid(T).name()) = tmp;
  //  vectorIncr_->operator[](typeid(T).name()) = incrementList;
  //  bytes = 0;
  //  peak = 0;
  return sizeof(T);
}

#endif

#include "BaseNamespaceFooter.H"
#endif
