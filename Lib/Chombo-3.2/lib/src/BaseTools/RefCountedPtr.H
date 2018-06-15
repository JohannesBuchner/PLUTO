#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _REFCOUNTEDPTR_H_
#define _REFCOUNTEDPTR_H_

#include <string>
#include <typeinfo>
#include "MayDay.H"
#include "Arena.H"
#include "CH_assert.H"
#include "BaseNamespaceHeader.H"

///  A reference-counting handle class.
/**
        This is to be used as a pointer to class T.
        This will feel and smell just like a built-in pointer except:
        -#  There is no need to call delete on the pointer.
        -#  The default copy constructor and assignment implement ref-counting.
        -#  The user may call isNonUnique to determine if this pointer is
            the only pointer to the data.  This can be used to hide the
            ref-counting behavior of a class.
        -#  Checks for dereference of a null pointer.

This class is completely inlined.

typical usage:

\code

{
  Box b;
  IntVect a;
  //refCount() == 1
  RefCountedPtr<IntVectSet> set(new IntVectSet());

  // still just one IntVectSet, but refCount()==2
  RefCountedPtr<IntVectSet>  otherSet(set);

  // Like a pointer, modifying one modifies the other
  otherSet->define(b);
  (*set)|=a;
  {
    RefCountedPtr<IntVectSet> anotherSet; // null
    anotherSet = otherSet ;  //now all three have refCount()==3
  }//anotherSet out of scope, so deleted.  IntVectSet NOT deleted.

  // set.refCount() == 2
  // otherSet.refCount() == 2;
  // otherset == set;

}
// when all RefCountedPtr's for a given object are deleted, the last
// one calls 'delete' on the member pointer.
\endcode

*/
template <class T>
class RefCountedPtr
{
public:
  //
  explicit inline RefCountedPtr(T* ptr = 0);
  //
  inline RefCountedPtr(const RefCountedPtr<T>& other);
  //
  inline ~RefCountedPtr();
  /// assignement operator.  copies pointer member
  /** copies pointer member and integer pointer, decreases refcount of rhs member
      before assignment.  this refcount increased my one. */
  inline const RefCountedPtr<T>& operator =(const RefCountedPtr<T>& rhs);
  /// dereference access operator.  use like a pointer derefence access function.
  inline T* operator ->();
  inline const T* operator ->() const;
  //
  inline bool isNull() const;
  //
  /// pointer dereference.
  inline T& operator *();

  inline const T& operator *() const;

  //
  /// auto conversion to regular const pointer where required.
  inline operator const T* () const;

  /// true when refcount is one.
  inline bool isNonUnique() const;

  inline int refCount() const;

  inline void swap(RefCountedPtr<T>& b);

  /// Template Conversion
  /**
     provides a mechanism for the compiler to generate a conversion to a base class, like
     a real pointer would have
  */
  template<class T2> inline operator RefCountedPtr<T2>();

  // Will never delete ptr_.  Useful in functions that need to return a
  // pointer which, under some circumstances, we want to own, and under others
  // we don't.
  inline void neverDelete();

  template<class T2>  friend class RefCountedPtr;

  static const std::string& name()
  {
    return name_;
  }

  inline void increment();
  inline void decrement();
protected:
  T* ptr_;
  int* refCount_;
  static std::string name_;
  static BArena s_Arena;
  static int size_;
private:

};

template<class T>
std::string RefCountedPtr<T>::name_(std::string("RefCountedPtr ")+std::string(typeid(T).name()));

//template<class T>
//BArena RefCountedPtr<T>::s_Arena(name_.c_str());
template<class T>
BArena RefCountedPtr<T>::s_Arena(typeid(T).name());

template<class T>
int RefCountedPtr<T>::size_ = sizeof(T)+2*sizeof(int);

template <class T> inline
void RefCountedPtr<T>::increment()
{
#ifdef CH_USE_MEMORY_TRACKING
  s_Arena.bytes += size_;
  refCount_[1] = size_;
  if (s_Arena.bytes > s_Arena.peak)
    {
      s_Arena.peak = s_Arena.bytes;
    }
#endif
}
template <class T> inline
void RefCountedPtr<T>::decrement()
{
#ifdef CH_USE_MEMORY_TRACKING
  s_Arena.bytes -= refCount_[1];
#endif
}

template <class T> inline
RefCountedPtr<T>::RefCountedPtr(T* ptr)
  : ptr_(ptr),
    refCount_(0)
{
  if (ptr_)
    {
      refCount_ = new int[2];
      if (refCount_ == 0)
        MayDay::Error("RefCountedPtr::RefCountedPtr(T* ptr) out of memory");
      *refCount_ = 1;
      increment();
    }
}

template <class T> inline
RefCountedPtr<T>::RefCountedPtr(const RefCountedPtr<T>& other)
  : ptr_(other.ptr_),
    refCount_(other.refCount_)
{
  if (refCount_ != 0)
    {
#pragma omp atomic
      ++(*refCount_);
    }
}

template <class T> inline
RefCountedPtr<T>::~RefCountedPtr()
{
  if (refCount_ != 0)
    { 
      int count;
      //
      //#pragma omp atomic capture
      count = --(*refCount_); 
      if(count==0)
	{
	  decrement();
	  delete ptr_;
	  ptr_ = 0;
	  delete [] refCount_;
	  refCount_ = 0;	  
	}
    }
}

template <class T> inline
const RefCountedPtr<T>& RefCountedPtr<T>::operator =(const RefCountedPtr<T>& rhs)
{
  if (ptr_ != rhs.ptr_)
    {
      if (refCount_ != 0)
	{
	  int count;
          //#pragma omp atomic capture
	  count = --(*refCount_); 
	  if(count == 0)
	    {
	      decrement();
	      delete ptr_;
	      delete [] refCount_;
	    }
	}
      ptr_ = rhs.ptr_;
      refCount_ = rhs.refCount_;
      if (refCount_ != 0)
#pragma omp atomic
        ++(*refCount_);
    }
  return *this;
}

template <class T> inline
bool RefCountedPtr<T>::isNull() const
{
  return (ptr_ == 0);
}

template <class T> inline 
T* RefCountedPtr<T>::operator ->()
{
  if (ptr_ == 0)
    MayDay::Error("RefCountedPtr<T>::operator ->() on null pointer");
  return ptr_;
}

template <class T> inline
const T* RefCountedPtr<T>::operator ->() const
{
  if (ptr_ == 0)
    MayDay::Error("RefCountedPtr<T>::operator ->() on null pointer");
  return ptr_;
}

template <class T> inline
T& RefCountedPtr<T>::operator *()
{
  if (ptr_ == 0)
    MayDay::Error("RefCountedPtr<T>::operator *() on null pointer");
  return *ptr_;
}

template <class T> inline const
T& RefCountedPtr<T>::operator *() const
{
  if (ptr_ == 0)
    MayDay::Error("RefCountedPtr<T>::operator *() on null pointer");
  return *ptr_;
}

template <class T> inline
RefCountedPtr<T>::operator const T* () const
{
  return ptr_;
}

template <class T> inline
bool RefCountedPtr<T>::isNonUnique() const
{
  return refCount_ == 0 ? false : *refCount_ != 1;
}

template <class T> inline
int RefCountedPtr<T>::refCount() const
{
  return refCount_ == 0 ? 0 : *refCount_;
}

template <class T>
template <class T2>
inline RefCountedPtr<T>::operator RefCountedPtr<T2> ()
{
  T2* t2 = dynamic_cast<T2*>(ptr_); // Handles upward and downward casts... -JNJ
  CH_assert(t2 != NULL); // ... unless the cast fails, of course. :-)
  RefCountedPtr<T2> rtn;  // if you have a compiler error that points here, then T2 is probably not a T
  rtn.ptr_= t2;
#ifdef  CH_USE_MEMORY_TRACKING
  // by fiat the derived classes memory is transfered to the base class Arena
  s_Arena.bytes      -=refCount_[1];
  rtn.s_Arena.bytes  +=refCount_[1];
#endif

  ++ *refCount_;
  rtn.refCount_ = refCount_;
  return rtn;
}

template <class T> inline
void RefCountedPtr<T>::swap(RefCountedPtr<T>& b)
{
  T* s = ptr_;
  ptr_ = b.ptr_;
  b.ptr_ = s;

  int* r = refCount_;
  refCount_ = b.refCount_;
  b.refCount_ = r;

}

template <class T> inline
void RefCountedPtr<T>::neverDelete()
{
  ++ *refCount_;
}

#include "BaseNamespaceFooter.H"
#endif
