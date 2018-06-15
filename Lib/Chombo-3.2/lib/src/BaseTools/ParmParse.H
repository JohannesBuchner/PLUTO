#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PARMPARSE_H_
#define _PARMPARSE_H_

/*
  DISCLAIMER:
  This is ugly because it consists of four files from boxlib
  concatenated together and munged hopelessly.  This was done
  (1) to greatly reduce the size of the API of boxlib
  (2) to avoid the godawful task of rewriting parmparse
  (3) to avoid cluttering the global namespace
  We deeply apologize.

  Any class that begins with PP_ is a convenience
  and will not be in any way supported by anyone at ANAG.
*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <sstream>

#include "CH_assert.H"
#include "MayDay.H"
#include "Misc.H"
#include "Vector.H"
#include "BaseNamespaceHeader.H"
using std::cout;
using std::cerr;
using std::endl;
//#include <aString.H>

//#include <Pointers.H>

//#include <UseCount.H>

//
//@Man:
//@Memo: A Class Encapsulating Reference Counts for ParmParse
/*@Doc:

This class encapsulates reference counts.

This is a convenience class for ParmParse
and will not be in any way supported by anyone at ANAG.
*/
#ifndef DOXYGEN
class PP_UseCount
{
public:
  //
  ///: Construct a PP_UseCount initialized to one.
  //
  PP_UseCount ();
  //
  ///: The copy constructor -- bumps reference count.
  //
  PP_UseCount (const PP_UseCount& rhs);

  /**: The assignment operator.  Increments the reference count on
    rhs, decrements the reference count on this, and then makes
    this and rhs point to the same reference count.
  */
  PP_UseCount& operator= (const PP_UseCount& rhs);
  //
  ///: The destructor -- decrements the reference count.
  //
  ~PP_UseCount ();
  //
  ///: The PP_UseCount is unique if the reference count == 1.
  //
  bool unique () const;
  //
  ///: Returns the reference count.
  //
  int linkCount () const;

private:
  //
  // A pointer to the reference count.
  //
  unsigned int* cnt;
  //
  // Decrement the reference count and delete the reference
  // counter if there are no more references.
  //
  void decrement ();
};

//
// Inlines.
//

inline
PP_UseCount::PP_UseCount ()
  :
  cnt(new unsigned int(1))
{
}

inline
PP_UseCount::PP_UseCount (const PP_UseCount& rhs)
  :
  cnt(rhs.cnt)
{
  ++*cnt;
}

inline
bool
PP_UseCount::unique () const
{
  return *cnt == 1;
}

inline
void
PP_UseCount::decrement ()
{
  if (unique())
    {
      delete cnt;
      cnt = 0;
    }
  else
    --*cnt;
}

inline
PP_UseCount&
PP_UseCount::operator= (const PP_UseCount& rhs)
{
  ++*rhs.cnt;
  decrement();
  cnt = rhs.cnt;
  return *this;
}

inline
PP_UseCount::~PP_UseCount ()
{
  decrement();
}

inline
int
PP_UseCount::linkCount () const
{
  return *cnt;
}

///: A Smart Pointer for Intrinsic or User-Defined Types for ParmParse
/**

The template class PP_CpPtr<T> provides a simple wrapper around a pointer
to type T (T*) that builds a copy of the pointed-to object when copied
from one PP_CpPtr<T> to another.  This is in contrast to a reference-counted
pointer class that would maintain one pointed-to object with a reference
count indicating the number of references.  Hence we call this a
"copied" smart pointer class.  It is intended for use with any type
type T, including the intrinsic types.  This class does not supply
an operator->(), as such an operator on intrinsic types has only recently
become a part of the C++ language, and many compilers do not yet implement
it.

This is a convenience class  for ParmParse
and will not be in any way supported by anyone at ANAG.
*/

template <class T>
class PP_CpPtr
{
public:
  //
  ///: The default constructor.  The wrapped pointer is null.
  //
  PP_CpPtr ();
  //
  ///: Construct a PP_CpPtr<T> setting the wrapped pointer to rhs.
  //
  /*explicit*/ PP_CpPtr (T* rhs);
  //
  ///: The destructor.  Deletes the wrapped pointer.
  //
  ~PP_CpPtr ();

  /**: The copy constructor.  If the pointer wrapped by rhs is null,
    the wrapped pointer is null here as well.  Otherwise,
    the contained pointer here is set to a new'd copy of that
    wrapped by rhs, with the two pointed-to values being identical.
    This assumes that type T has a well-defined and accessible copy
    constructor.  T must also be a concrete type, not a abstract
    type.
  */
  PP_CpPtr (const PP_CpPtr<T>& rhs);

  /**: Sets the wrapped pointer to rhs.  Deletes the previously
    wrapped pointer.
  */
  PP_CpPtr<T>& operator= (T* rhs);

  /**: The copy assignment operator.  If the pointer wrapped by rhs
    is null, the wrapped pointer is null here as well.  Otherwise,
    the contained pointer here is set to a new'd copy of that
    wrapped by rhs, with the two pointed-to values being identical.
    This assumes that type T has a well-defined and accessible copy
    constructor.  T must also be a concrete type, not a abstract
    type.
  */
  PP_CpPtr<T>& operator= (const PP_CpPtr<T>& rhs);

  /**: Returns a reference to the value pointed to by the wrapped
    pointer; i.e. dereferencing this PP_CpPtr<T>, returns the
    dereferenced wrapped pointer.  It is an error if the wrapped
    pointer is null.
  */
  T& operator* () const;
  //
  ///: Returns true if the wrapped pointer null.
  //
  bool isNull () const;
  //
  ///: Sets the wrapped pointer to null and returns the previous value.
  //
  T* release ();
  //
  ///: Are the two pointers (not the values to which they point) equal?
  //
  bool operator== (const PP_CpPtr<T>& rhs) const;
  //
  ///: Are the two pointers not equal?
  //
  bool operator!= (const PP_CpPtr<T>& rhs) const;

protected:
  T* ptr;
};

/// A Smart Pointer for User-Defined Types for ParmParse
/**

The template class PP_CpClassPtr<T> is derived from PP_CpPtr<T>.  It provides a
simple wrapper around a pointer to type T (a T*) that "does the right thing"
when copied from one PP_CpPtr<T> to another.  The type T MUST be a user-defined
type, not an intrinsic type.  Given this restriction, we can supply
an operator->().

This is a convenience class  for ParmParse
and will not be in any way supported by anyone at ANAG.
*/

template<class T>
class PP_CpClassPtr :
  public PP_CpPtr<T>
{
public:
  //
  ///: The default constructor.  The wrapped pointer is null.
  //
  PP_CpClassPtr ();
  //
  ///: Construct a PP_CpPtr<T> setting the wrapped pointer to rhs.
  //
  /*explicit*/ PP_CpClassPtr (T* rhs);

  /**: The copy constructor.  If the pointer wrapped by rhs is null,
    the wrapped pointer is null here as well.  Otherwise,
    the contained pointer here is set to a new'd copy of that
    wrapped by rhs, with the two pointed-to values being identical.
    This assumes that type T has a well-defined and accessible copy
    constructor.  T must also be a concrete type, not a abstract
    type.
  */
  PP_CpClassPtr (const PP_CpClassPtr<T>& rhs);

  /**: Sets the wrapped pointer to rhs.  Deletes the previously
    wrapped pointer.
  */
  PP_CpClassPtr<T>& operator= (T* rhs);

  /**: The copy assignment operator.  If the pointer wrapped by rhs
    is null, the wrapped pointer is null here as well.  Otherwise,
    the contained pointer here is set to a new'd copy of that
    wrapped by rhs, with the two pointed-to values being identical.
    This assumes that type T has a well-defined and accessible copy
    constructor.  T must also be a concrete type, not a abstract
    type.
  */
  PP_CpClassPtr<T>& operator= (const PP_CpClassPtr<T>& rhs);
  //
  ///: Applies operator-> to the wrapped pointer.
  //
  T* operator-> () const;
};

/// A Reference Counted Smart Pointer for Intrinsic or User-Defined Types for ParmParse
/**

The template class PP_LnPtr<T> provides a reference counted wrapper around a
pointer to type T (a T*).  This "smart" pointer is intended for use with
any type type T, including the intrinsic types.  For this reason, we do
not supply an operator->(), as such an operator on intrinsic types has only
recently become a part of the C++ language and many compilers do not yet
implement it.

This is a convenience class  for ParmParse
and will not be in any way supported by anyone at ANAG.
*/

template<class T>
class PP_LnPtr
{
public:
  //
  ///: The default constructor.  The wrapped pointer is null.
  //
  PP_LnPtr ();
  //
  ///: Construct a PP_LnPtr<T> setting the wrapped pointer to rhs.
  //
  /*explicit*/ PP_LnPtr (T* rhs);

  /**: The copy assignment operator.  The contained pointer is set
    to the one wrapped by rhs.  The reference count is decremented
    on this object and the reference count is incremented for
    the newly wrapped pointer.
  */
  PP_LnPtr<T>& operator= (const PP_LnPtr<T>& rhs);

  /**: Sets the wrapped pointer to rhs.  Decrements the count
    on the previously wrapped pointer and deletes it if there
    was only one reference.
  */
  PP_LnPtr<T>& operator= (T* rhs)
  {
    //
    // This is inlined here as OCC won't have it any other way :-(
    //
    if (unique())
      delete ptr;
    ptr = rhs;
    ucnt = PP_UseCount();
    return *this;
  }

  /**: The destructor -- decrements the reference count and deletes
    the wrapped pointer if there is only one reference.
  */
  ~PP_LnPtr ();
  //
  ///: Returns true if only one reference to the wrapped pointer.
  //
  bool unique () const;
  //
  ///: Returns the number of references to the wrapped pointer.
  //
  int linkCount () const;

  /**: Returns a reference to the value pointed to by the wrapped
    pointer; i.e. dereferencing this PP_LnPtr<T>, returns the
    dereferenced wrapped pointer.  It is an error if the wrapped
    pointer is null.
  */
  T& operator* () const;
  //
  ///: Returns true if the wrapped pointer is null.
  //
  bool isNull () const;
  //
  ///: Are the two pointers (not the values to which they point) equal?
  //
  bool operator== (const PP_LnPtr<T>& rhs) const;
  //
  ///: Are the two pointers not equal?
  //
  bool operator!= (const PP_LnPtr<T>& rhs) const;

protected:
  T*       ptr;

private:
  PP_UseCount ucnt;
};

/// A Smart Reference Counted Pointer for User-Defined Types for ParmParse
/**

The template class PP_LnClassPtr<T> is derived from PP_LnPtr<T>.  It provides a
reference counted wrapper around a pointer to type T (a T*).  The type T
MUST be a user-defined type, not an intrinsic type.  Given this
restriction, we can supply an operator->().

This is a convenience class  for ParmParse
and will not be in any way supported by anyone at ANAG.
*/

template<class T>
class PP_LnClassPtr
  : public PP_LnPtr<T>
{
public:
  //
  ///: The default constructor.  The wrapped pointer is null.
  //
  PP_LnClassPtr ();
  //
  ///: Construct a PP_LnPtr<T> setting the wrapped pointer to rhs.
  //
  /*explicit*/ PP_LnClassPtr (T* rhs);

  /**: The copy assignment operator.  The contained pointer is set
    to the one wrapped by rhs.  The reference count is decremented
    on this object and the reference count is incremented for
    the newly wrapped pointer.
  */
  PP_LnClassPtr<T>& operator= (const PP_LnClassPtr<T>& rhs);

  /**: Sets the wrapped pointer to rhs.  Decrements the count
    on the previously wrapped pointer and deletes it if there
    was only one reference.
  */
  PP_LnClassPtr<T>& operator= (T* rhs);
  //
  ///: Applies operator-> to the wrapped pointer.
  //
  T* operator->() const;
};

//
// Inlines.
//

template <class T>
inline
PP_CpPtr<T>::PP_CpPtr ()
  :
  ptr(0)
{
}

template <class T>
inline
PP_CpPtr<T>::PP_CpPtr(T* rhs)
  :
  ptr(rhs)
{
}

template <class T>
inline
PP_CpPtr<T>::~PP_CpPtr()
{
  delete ptr;
}

template <class T>
inline
bool
PP_CpPtr<T>::isNull () const
{
  return ptr == 0;
}

template <class T>
inline
PP_CpPtr<T>::PP_CpPtr (const PP_CpPtr<T>& rhs)
{
  ptr = rhs.isNull() ?  0 : new T(*rhs.ptr);
}

template <class T>
inline
PP_CpPtr<T>&
PP_CpPtr<T>::operator= (const PP_CpPtr<T>& rhs)
{
  if (!(ptr == rhs.ptr))
    {
      delete ptr;
      ptr = rhs.isNull() ? 0 : new T(*rhs.ptr);
    }
  return *this;
}

template <class T>
inline
PP_CpPtr<T>&
PP_CpPtr<T>::operator= (T* rhs)
{
  delete ptr;
  ptr = rhs;
  return *this;
}

template <class T>
inline
T&
PP_CpPtr<T>::operator* () const
{
  CH_assert(ptr != 0);
  return *ptr;
}

template <class T>
inline
T*
PP_CpPtr<T>::release ()
{
  T* old = ptr;
  ptr = 0;
  return old;
}

template <class T>
inline
bool
PP_CpPtr<T>::operator== (const PP_CpPtr<T>& rhs) const
{
  return ptr == rhs.ptr;
}

template <class T>
inline
bool
PP_CpPtr<T>::operator!= (const PP_CpPtr<T>& rhs) const
{
  return ptr != rhs.ptr;
}

template <class T>
inline
PP_CpClassPtr<T>::PP_CpClassPtr ()
  :
  PP_CpPtr<T>()
{
}

template <class T>
inline
PP_CpClassPtr<T>::PP_CpClassPtr (T* rhs)
  :
  PP_CpPtr<T>(rhs)
{
}

template <class T>
inline
PP_CpClassPtr<T>::PP_CpClassPtr (const PP_CpClassPtr<T>& rhs)
  :
  PP_CpPtr<T>(rhs)
{
}

template <class T>
inline
PP_CpClassPtr<T>&
PP_CpClassPtr<T>::operator= (T* rhs)
{
  PP_CpPtr<T>::operator= (rhs);
  return *this;
}

template <class T>
inline
PP_CpClassPtr<T>&
PP_CpClassPtr<T>::operator= (const PP_CpClassPtr<T>& rhs)
{
  PP_CpPtr<T>::operator= (rhs);
  return *this;
}

template <class T>
inline
T*
PP_CpClassPtr<T>::operator-> () const
{
  return this->ptr;
}

template <class T>
inline
PP_LnPtr<T>::PP_LnPtr ()
  :
  ptr(0)
{
}

template <class T>
inline
PP_LnPtr<T>::PP_LnPtr(T* rhs)
  :
  ptr(rhs)
{
}

template <class T>
inline
bool
PP_LnPtr<T>::unique () const
{
  return ucnt.unique();
}

template <class T>
inline
PP_LnPtr<T>::~PP_LnPtr ()
{
  if (ucnt.unique())
    delete ptr;
}

template <class T>
inline
PP_LnPtr<T>&
PP_LnPtr<T>::operator= (const PP_LnPtr<T>& rhs)
{
  if (ptr != rhs.ptr)
    {
      if (unique())
        delete ptr;
      ptr = rhs.ptr;
      ucnt = rhs.ucnt;
    }
  return *this;
}

template <class T>
inline
int
PP_LnPtr<T>::linkCount () const
{
  return ucnt.linkCount();
}

template <class T>
inline
T&
PP_LnPtr<T>::operator* () const
{
  CH_assert(ptr != 0);
  return *ptr;
}

template <class T>
inline
bool
PP_LnPtr<T>::isNull () const
{
  return ptr == 0;
}

template <class T>
inline
bool
PP_LnPtr<T>::operator== (const PP_LnPtr<T>& rhs) const
{
  return ptr == rhs.ptr;
}

template <class T>
inline
bool
PP_LnPtr<T>::operator!= (const PP_LnPtr<T>& rhs) const
{
  return ptr != rhs.ptr;
}

template <class T>
inline
PP_LnClassPtr<T>::PP_LnClassPtr ()
{
}

template <class T>
inline
PP_LnClassPtr<T>::PP_LnClassPtr (T* rhs)
  :
  PP_LnPtr<T>(rhs)
{
}

template <class T>
inline
PP_LnClassPtr<T>&
PP_LnClassPtr<T>::operator= (const PP_LnClassPtr<T>& rhs)
{
  PP_LnPtr<T>::operator=(rhs);
  return *this;
}

template <class T>
inline
PP_LnClassPtr<T>&
PP_LnClassPtr<T>::operator= (T* rhs)
{
  PP_LnPtr<T>::operator=(rhs);
  return *this;
}

template <class T>
inline
T*
PP_LnClassPtr<T>::operator->() const
{
  return this->ptr;
}

//
// Helper class for class PP_String.
//

class PP_StringRep
{
  friend class PP_String;
  char* s;
  int   bufferlength;
public:
  /*explicit*/ PP_StringRep (int _len = 0);
  ~PP_StringRep ();
  //
  // Resized the buffer and copies the contents of old buffer to new one.
  //
  void resize (int n);
};

//
// PP_StringRep inlines.
//

inline
PP_StringRep::PP_StringRep (int _len)
{
  bufferlength = _len;
  s = new char [bufferlength];
}

inline
PP_StringRep::~PP_StringRep ()
{
  delete [] s;
  s = 0;
}

/// A String Class for ParmParse
/**

The class PP_String is used to store and manipulate character strings.  It
has an efficient underlying storage mechanism and some useful string
manipulation operations.

The PP_String class is implemented using a character array and reference
count. Two PP_Strings may reference the same underlying character array with
a reference count of two.  When an PP_String copy constructor or copy
operator is applied the reference count on the underlying character array
is incremented but the actual string is not copied.  That is, copying an
PP_String is an inexpensive operation.  When an PP_String is destructed, the
reference count is decremented.  The underlying character array is deleted
only when the reference count goes to zero.  Any operator that modifies an
PP_String will make its own copy of the character array before the
modification, unless it's the sole owner of the character array in the
PP_String.

This is a convenience class  for ParmParse
and will not be in any way supported by anyone at ANAG.

*/

class PP_String
{
public:
  //
  ///: Constructs an empty string.
  //
  PP_String ();

  /**: Constructs an PP_String containing the single character c.
    If c is the null character then this will be the empty string.
  */
  /*explicit*/ PP_String (char c);

  /**: Constructs an empty PP_String but allocates enough space to hold
    a character string of length len.  This may be useful when
    reading in very long lines with the `getline' function; i.e.
    it may be slightly more efficient to allocate a very large
    PP_String once and then use `getline' than to use `getline'
    on an empty string.  In general though, you won't notice any
    difference and you really shouldn't use this constructor.
  */
  /*explicit*/ PP_String (int len);

  /**: Constructs an PP_String initialized to the character string s.
    It is an error is `s' is the null string.
  */
  PP_String (const char* s);
  //
  ///: The copy constructor.
  //
  PP_String (const PP_String& rhs);
  //
  ///: The assignment operator
  //
  PP_String& operator= (const PP_String& rhs);

  /**: Catenate PP_String right onto end of this PP_String.  Return a
    reference to self to allow for operator chaining.
  */
  PP_String& operator+= (const PP_String& right);

  /**: Catenate character string right onto end of this PP_String.
    Returns a reference to self to allow for operator chaining.
    It is an error is `right' is the null string.
  */
  PP_String& operator+= (const char* right);

  /**: Catenate character c onto end of this PP_String.
    Returns a reference to self to allow for operator chaining.
    This does nothing if c is the null character.
  */
  PP_String& operator+= (char c);

  /**: Converts all characters in this PP_String to upper case.
    Returns a reference to self to allow for operator chaining.
  */
  PP_String& toUpper ();

  /**: Converts all characters in this PP_String to lower case.
    Returns a reference to self to allow for operator chaining.
  */
  PP_String& toLower ();

  /**: Read the next line from the input stream strm and store it
    as the value of the string.  The delimiter for a line of text
    is the newline character.  The newline is extracted from the
    istream, but it is NOT stored in the PP_String.  There is no
    limit to the length of string that can be extracted.
  */
  std::istream& getline (std::istream& strm);

  /**: Returns the number of characters stored in this PP_String.
    This does not include the terminating null character.
  */
  int length () const;

  //
  ///: Returns true if this is the null string.
  //
  bool isNull () const;

  /**: Returns a reference to the kth character in the string.  An
    error occurs if k < 0 or k >= length().
  */
  char& operator [] (int k);

  /**: Returns kth character in the string.  An error occurs
    if k < 0 or k >= length().
  */
  char operator[] (int k) const;

  /**: Convert an PP_String to a const char *.  This allows an PP_String
    to be used in any context where a const char* type is needed.
  */
  const char* c_str () const;

  /**: Returns the value of the string as a double.  In the case
    when the string is empty or when the initial characters
    in the string are strictly alphabetic characters, this
    returns zero.  It's equivalent to `atof(c_str())'.
  */
  double toDouble () const;

  /**: Returns the value of the string as a integer.  In the case
    when the string is empty or when the initial characters
    in the string are strictly alphabetic characters, this
    returns zero.  It's equivalent to `atoi(c_str())'.
  */
  int toInteger () const;

  /**: Returns the value of the string as a long.  In the case
    when the string is empty or when the initial characters
    in the string are strictly alphabetic characters, this
    returns zero.  It's equivalent to `atol(c_str())'.
  */
  long toLong () const;
  //
  ///: Write to an ostream in ASCII format.
  //
  friend std::ostream& operator<< (std::ostream&       os,
                                   const PP_String& str);

  /**: Read a whitespace delimited string from an istream.
    This function discards leading whitespace and then reads
    in non-whitespace character until the next whitespace
    character or end-of-file.  Note that there is no limit,
    on the length of the character that can be read in, except
    that dictated by the resources of the machine.
    Note  also that operator>> and operator<< are not completely
    symmetrical in the case where operator<< writes out a
    string that contains whitespace.  If you're trying to
    read in a string that contains whitespace, you might
    want to use getline() instead.
  */
  friend std::istream& operator>> (std::istream& is,
                                   PP_String& str);

protected:
  void copyModify ();

private:
  PP_LnClassPtr<PP_StringRep> p;
  int                   len;

  //
  // None of the following functions need to be friends.  I've made
  // them friends solely so they'll show up nicely in the HTML documentaion
  // spewed out by doc++.
  //

  //
  // Is left lexically less than right?
  //
  friend inline bool operator<  (const PP_String& left,
                                 const PP_String& right);
  //
  // Is left lexically less than or equal to right?
  //
  friend inline bool operator<= (const PP_String& left,
                                 const PP_String& right);
  //
  // Is left not equal to right?
  //
  friend inline bool operator!= (const PP_String& left,
                                 const PP_String& right);
  //
  // Is left equal to right?
  //
  friend inline bool operator== (const PP_String& left,
                                 const PP_String& right);
  //
  // Is left lexically greater than or equal to right?
  //
  friend inline bool operator>= (const PP_String& left,
                                 const PP_String& right);
  //
  // Is left lexically greater than right?
  //
  friend inline bool operator>  (const PP_String& left,
                                 const PP_String& right);
  //
  // Is left lexically less than right?
  //
  friend inline bool operator<  (const PP_String& left,
                                 const char*    right);
  //
  // Is left lexically less than or equal to right?
  //
  friend inline bool operator<= (const PP_String& left,
                                 const char*    right);
  //
  // Is left not equal to right?
  //
  friend inline bool operator!= (const PP_String& left,
                                 const char*    right);
  //
  // Is left equal to right?
  //
  friend inline bool operator== (const PP_String& left,
                                 const char*    right);
  //
  // Is left lexically greater than or equal to right?
  //
  friend inline bool operator>= (const PP_String& left,
                                 const char*    right);
  //
  // Is left lexically greater than right?
  //
  friend inline bool operator>  (const PP_String& left,
                                 const char*    right);
  //
  // Is left lexically less than right?
  //
  friend inline bool operator<  (const char*    left,
                                 const PP_String& right);
  //
  // Is left lexically less than or equal to right?
  //
  friend inline bool operator<= (const char*    left,
                                 const PP_String& right);
  //
  // Is left not equal to right?
  //
  friend inline bool operator!= (const char*    left,
                                 const PP_String& right);
  //
  // Is left equal to right?
  //
  friend inline bool operator== (const char*    left,
                                 const PP_String& right);
  //
  // Is left lexically greater than or equal to right?
  //
  friend inline bool operator>= (const char*    left,
                                 const PP_String& right);
  //
  // Is left lexically greater than right?
  //
  friend inline bool operator>  (const char*    left,
                                 const PP_String& right);
};

//
// PP_String inlines.
//

inline
bool
PP_String::isNull () const
{
  return len == 0;
}

inline
int
PP_String::length () const
{
  return len;
}

inline
double
PP_String::toDouble () const
{
  return len == 0 ? 0 : std::atof(p->s);
}

inline
int
PP_String::toInteger () const
{
  return len == 0 ? 0 : atoi(p->s);
}

inline
long
PP_String::toLong () const
{
  return len == 0 ? 0 : atol(p->s);
}

inline
const char*
PP_String::c_str () const
{
  return p->s;
}

inline
char
PP_String::operator[] (int index) const
{
  CH_assert(index >=0 && index < len);
  return p->s[index];
}

inline
PP_String
operator+ (const PP_String& left,
           const PP_String& right)
{
  PP_String result(left);
  return result += right;
}

inline
bool
operator< (const PP_String& left,
           const PP_String& right)
{
  return ::strcmp(left.c_str(), right.c_str()) < 0;
}

inline
bool
operator<= (const PP_String& left,
            const PP_String& right)
{
  return ::strcmp(left.c_str(), right.c_str()) <= 0;
}

inline
bool
operator!= (const PP_String& left,
            const PP_String& right)
{
  return ::strcmp(left.c_str(), right.c_str()) != 0;
}

inline
bool
operator== (const PP_String& left,
            const PP_String& right)
{
  return ::strcmp(left.c_str(), right.c_str()) == 0;
}

inline
bool
operator>= (const PP_String& left,
            const PP_String& right)
{
  return ::strcmp(left.c_str(), right.c_str()) >= 0;
}

inline
bool
operator>  (const PP_String& left,
            const PP_String& right)
{
  return ::strcmp(left.c_str(), right.c_str()) > 0;
}

inline
bool
operator< (const PP_String& left,
           const char*    right)
{
  return ::strcmp(left.c_str(), right) < 0;
}

inline
bool
operator<= (const PP_String& left,
            const char*    right)
{
  return ::strcmp(left.c_str(), right) <= 0;
}

inline
bool
operator!= (const PP_String& left,
            const char*    right)
{
  return ::strcmp(left.c_str(), right) != 0;
}

inline
bool
operator== (const PP_String& left,
            const char*    right)
{
  return ::strcmp(left.c_str(), right) == 0;
}

inline
bool
operator>= (const PP_String& left,
            const char*    right)
{
  return ::strcmp(left.c_str(), right) >= 0;
}

inline
bool
operator>  (const PP_String& left,
            const char*    right)
{
  return ::strcmp(left.c_str(), right) > 0;
}

inline
bool
operator< (const char*    left,
           const PP_String& right)
{
  return ::strcmp(left, right.c_str()) < 0;
}

inline
bool
operator<= (const char*    left,
            const PP_String& right)
{
  return ::strcmp(left, right.c_str()) <= 0;
}

inline
bool
operator!= (const char*    left,
            const PP_String& right)
{
  return ::strcmp(left, right.c_str()) != 0;
}

inline
bool
operator== (const char*    left,
            const PP_String& right)
{
  return ::strcmp(left, right.c_str()) == 0;
}

inline
bool
operator>= (const char*    left,
            const PP_String& right)
{
  return ::strcmp(left, right.c_str()) >= 0;
}

inline
bool
operator>  (const char*    left,
            const PP_String& right)
{
  return ::strcmp(left, right.c_str()) > 0;
}

//#include <List.H>

template <class T> class PP_ListLink;
template <class T> class PP_ListIterator;
template <class T> class PP_List;

template <class T>
class PP_ListLink
{
private:
  friend class PP_List<T>;
  friend class PP_ListIterator<T>;

  PP_ListLink (const T&     _val,
               PP_ListLink<T>* _pre,
               PP_ListLink<T>* _suc)
    :
    val(_val),
    pre(_pre),
    suc(_suc)
  {}

private:
  T            val;
  PP_ListLink<T>* pre;
  PP_ListLink<T>* suc;
};

///Iterate over a List for ParmParse
/**

The class PP_ListIterator<T> is an iterator over class PP_List<T>.

This class does NOT provide a default constructor or an assignment operator.

This is a convenience class  for ParmParse
and will not be in any way supported by anyone at ANAG.
*/

template <class T>
class PP_ListIterator
{
public:
  //
  ///: Construct a PP_ListIterator<T> to first element of aList.
  //
  /*explicit*/ PP_ListIterator (const PP_List<T>& aList);
  //
  ///: The copy constructor.
  //
  PP_ListIterator (const PP_ListIterator<T>& rhs);

  /**: Reset this PP_ListIterator<T> to point to the first element
    in the PP_List<T>.
  */
  void rewind ();

  /**: Return a constant reference to the object in the PP_List<T>
    currently pointed to by this PP_ListIterator<T>.
  */
  const T& operator() () const;

  /**: Return a constant reference to the object in the PP_List<T>
    currently pointed to by this PP_ListIterator<T>.
  */
  const T& operator* () const;

  /**: This is a conversion operator that makes the iterator look
    like a pointer.  This operator makes it easy to check if the
    iterator is pointing to an element on the PP_List<T>.  If the
    iterator has been moved off the PP_List<T> or if the PP_List<T> is
    empty, this conversion returns the NULL pointer.
  */
  operator void* ();

  /**: Returns true if PP_ListIterator<T> doesn't point to any
    element on the PP_List<T>.
  */
  bool operator! () const;

  /**: Return a constant reference to the object in the PP_List<T>
    currently pointed to by the iterator.
  */
  const T& value () const;

  /**: This is the prefix auto-increment operator.  It advances the
    PP_ListIterator<T> to point to the next element on the PP_List<T>.
    It then returns a reference to itself to allow for chaining
    with other operators.
  */
  PP_ListIterator<T>& operator++ ();

  /**: This is the prefix auto-decrement operator.  It moves the
    PP_ListIterator<T> to point to the previous element on the
    PP_List<T>.  It then returns a reference to itself to allow for
    chaining with other operators.
  */
  PP_ListIterator<T>& operator-- ();

  /**: This is the postfix auto-decrement operator.  It moves the
    PP_ListIterator<T> to point to the previous element on the
    PP_List<T>.  It then returns a PP_ListIterator<T> that points to
    the old element to allow for chaining with other operators.

  */
  PP_ListIterator<T> operator-- (int);

  /**: This is the postfix auto-increment operator.  It advances the
    PP_ListIterator<T> to point to the next element on the PP_List<T>.
    It then returns a PP_ListIterator<T> that points to the old
    element to allow for chaining with other operators.
  */
  PP_ListIterator<T> operator++ (int);

  /**: Do the two PP_ListIterator<T>s point to the same PP_List<T> and
    the same element within the PP_List<T>?
  */
  bool operator== (const PP_ListIterator<T>&) const;
  //
  ///: Are the PP_ListIterator<T>s not equal?
  //
  bool operator!= (const PP_ListIterator<T>&) const;

protected:
  //
  // Construct a PP_ListIterator<T> to a PP_List<T> and object in that PP_List<T>.
  //
  PP_ListIterator (const PP_List<T>& _list,
                   PP_ListLink<T>*   _p);
  //
  // A reference to the PP_List<T> to which we point.
  //
  const PP_List<T>& list;
  //
  // A pointer to the element in the PP_List<T> to which we point.
  //
  PP_ListLink<T>* p;

private:
  friend class PP_List<T>;
  //
  // These are disallowed.
  //
  PP_ListIterator ();
  PP_ListIterator<T>& operator= (const PP_ListIterator<T>&);
};

/// A Doubly-Linked List for ParmParse
/**

The PP_List<T> class is a template class that implements a doubly-linked list
of objects.  A PP_List<T> is a useful container class when the number of
objects in the collection is not known ahead of time.   A PP_List<T> can
contain an arbitrary number of elements; operations such as insertion,
deletion, and catenation are easily implemented and inexpensive.

The only difficulty when defining a list class is devising a mechanism to
access the elements.  In an array, an element is accessed using an
integer index.  Since the elements in a PP_List<T> are ordered by position,
we could define an integer indexing operation that walks along the
PP_List<T> links from the beginning until the numbered element is found.
Unfortunately, this would be very inefficient when accessing elements near
the  end of a long list.  Another solution is to allow user access to the
individual link objects that contain the element as well as the forward and
backward pointers.  This is not a satisfactory solution since it allows
user access to the internal representation of the class.  The solution
chosen is to define a PP_ListIterator<T> template class.

Think of a PP_ListIterator<T> as a pointer to an object in the PP_List<T>.  You
can access the element currently pointed to by the iterator, move the
iterator forward and backward through the PP_List<T>, and use it as a
mechanism to define where elements should be inserted and deleted.  If the
iterator is moved off the end of the list it behaves as a null pointer.

This is a concrete class, not a polymorphic one.

This is a convenience class  for ParmParse
and will not be in any way supported by anyone at ANAG.
*/

template <class T>
class PP_List
{
public:
  //
  ///: Construct an empty PP_List<T>.
  //
  PP_List ();
  //
  ///: The copy constructor.
  //
  PP_List (const PP_List<T>& rhs);
  //
  ///: The assignment operator.
  //
  PP_List<T>& operator= (const PP_List<T>& rhs);
  //
  ///: The destructor.
  //
  ~PP_List();
  //
  ///: Adds a copy of the value to the beginning of the PP_List<T>.
  //
  void prepend (const T& value);
  //
  ///: Adds a copy of the value to the end of the PP_List<T>.
  //
  void append (const T& value);
  //
  ///: Adds a copy of the value to the end of the PP_List<T>.
  //
  void add (const T& value);
  //
  ///: Appends a copy of all items in PP_List<T> src to this PP_List<T>.
  //
  void join (const PP_List<T>& src);

  /**: Appends a copy of all items in PP_List<T> src to this PP_List<T>.
    This differs from join() in that it unlinks the objects from
    the PP_List<T> src and glues them to the end of this PP_List<T>,
    leaving PP_List<T> src empty.  This is more efficient that
    join() if src is no longer needed.
  */
  void catenate (PP_List<T>& src);
  //
  ///: Removes all objects from the PP_List<T>.
  //
  void clear ();
  //
  ///: Returns a reference to the first element in the PP_List<T>.
  //
  T& firstElement () const;
  //
  ///: Returns a reference to the last element in the PP_List<T>.
  //
  T& lastElement () const;

  /**: Returns true if the PP_List<T> contains an object identical
    to value.  Type T must have an operator==() defined, or
    be an intrinsic type.
  */
  bool includes (const T& value) const;

  /**: Returns true if the this and rhs are memberwise equal;
    i.e. the two lists are the same size and each of the
    elements in the list compare equal. Type T must have an
    operator==() defined, or be an intrinsic type.
  */
  bool operator== (const PP_List<T>& rhs) const;
  //
  ///: Returns true if the this and rhs are not equal.
  //
  bool operator!= (const PP_List<T>& rhs) const;
  //
  ///: Returns true if the PP_List<T> is empty.
  //
  bool isEmpty () const;
  //
  ///: Returns true if the PP_List<T> is not empty.
  //
  bool isNotEmpty () const;
  //
  ///: Returns the number of objects in the PP_List<T>.
  //
  int length () const;
  //
  ///: Removes the first element in the PP_List<T>.
  //
  void removeFirst ();
  //
  ///: Removes the last element in the PP_List<T>.
  //
  void removeLast ();
  //
  ///: Returns reference to object pointed to by the PP_ListIterator<T>.
  //
  const T& operator[] (const PP_ListIterator<T>& li) const;
  //
  ///: Returns reference to object pointed to by the PP_ListIterator<T>.
  //
  T& operator[] (const PP_ListIterator<T>& li);
  //
  ///: Removes all objects in the PP_List<T> equal to value.
  //
  void remove (const T& value);

  /**: Removes all objects in the PP_List<T> equal to any of the
    values in lst.
  */
  void remove (const PP_List<T>& lst);
  //
  ///: Removes the object pointed to by the PP_ListIterator<T>.
  //
  void remove (PP_ListIterator<T>& lit);
  //
  ///: Replace the value pointed to by the PP_ListIterator<T> by val.
  //
  void replace (PP_ListIterator<T>& li,
                const T&         val);

  /**: Insert val into PP_List<T> after the object pointed to by
    PP_ListIterator<T>.
  */
  void addAfter (PP_ListIterator<T>& lit,
                 const T&         val);

  /**: Insert val into PP_List<T> before the object pointed to by
    PP_ListIterator<T>.
  */
  void addBefore (PP_ListIterator<T>& lit,
                  const T&         val);
  //
  ///: Returns a PP_ListIterator<T> to the first object in this PP_List<T>.
  //
  PP_ListIterator<T> first () const;
  //
  ///: Returns a PP_ListIterator<T> to the last object in this PP_List<T>.
  //
  PP_ListIterator<T> last () const;

protected:
  //
  // A helper function for removing nodes.
  //
  void remove (PP_ListLink<T> *ln);
  //
  // A helper function for adding nodes.
  //
  PP_ListLink<T>* addBefore (PP_ListLink<T>* ln,
                             const T&     val);
  //
  // A helper function for adding nodes.
  //
  PP_ListLink<T>* addAfter (PP_ListLink<T>* ln,
                            const T&     val);
  //
  // The head of the list.
  //
  PP_ListLink<T>* head;
  //
  // The tail of the list.
  //
  PP_ListLink<T>* tail;
  //
  // Our good and trusted friend.
  //
  friend class PP_ListIterator<T>;
};

//
// Inlines.
//

//
// The ListIterator stuff.
//

template <class T>
inline
PP_ListIterator<T>::PP_ListIterator (const PP_List<T>& _list,
                                     PP_ListLink<T>*   _p)
  :
  list(_list),
  p(_p)
{
}

template <class T>
inline
PP_ListIterator<T>::PP_ListIterator (const PP_List<T>& aList)
  :
  list(aList)
{
  p = list.head;
}

template <class T>
inline
PP_ListIterator<T>::PP_ListIterator (const PP_ListIterator<T>& li)
  :
  list(li.list),
  p(li.p)
{
}

template <class T>
inline
void
PP_ListIterator<T>::rewind ()
{
  p = list.head;
}

template <class T>
inline
const T&
PP_ListIterator<T>::operator() () const
{
  CH_assert(p != 0);
  return p->val;
}

template <class T>
inline
const T&
PP_ListIterator<T>::operator* () const
{
  CH_assert(p != 0);
  return p->val;
}

template <class T>
inline
PP_ListIterator<T>::operator void* ()
{
  return p != 0 ? this : 0;
}

template <class T>
inline
bool
PP_ListIterator<T>::operator! () const
{
  return p == 0 ? true : false;
}

template <class T>
inline
const T&
PP_ListIterator<T>::value () const
{
  CH_assert(p != 0);
  return p->val;
}

template <class T>
inline
PP_ListIterator<T>&
PP_ListIterator<T>::operator++ ()
{
  if (p)
    p = p->suc;
  return *this;
}

template <class T>
inline
PP_ListIterator<T>&
PP_ListIterator<T>::operator-- ()
{
  if (p)
    p = p->pre;
  return *this;
}

template <class T>
inline
PP_ListIterator<T>
PP_ListIterator<T>::operator++ (int)
{
  const PP_ListIterator<T> li = *this;
  ++(*this);
  return li;
}

template <class T>
inline
PP_ListIterator<T>
PP_ListIterator<T>::operator-- (int)
{
  const PP_ListIterator<T> li = *this;
  --(*this);
  return li;
}

template <class T>
inline
bool
PP_ListIterator<T>::operator== (const PP_ListIterator<T>& _li) const
{
  return (&list == &_li.list && p == _li.p) ? true : false;
}

template <class T>
inline
bool
PP_ListIterator<T>::operator!= (const PP_ListIterator<T>& _li) const
{
  return ! PP_ListIterator<T>::operator==(_li);
}

//
// List stuff.
//

template <class T>
inline
PP_List<T>::PP_List ()
  :
  head(0),
  tail(0)
{
}

template <class T>
inline
PP_List<T>::~PP_List ()
{
  clear();
}

template <class T>
inline
void
PP_List<T>::prepend (const T& value)
{
  addBefore(head, value);
}

template <class T>
inline
void
PP_List<T>::append (const T& value)
{
  addAfter(tail, value);
}

template <class T>
inline
T&
PP_List<T>::firstElement () const
{
  CH_assert(head != 0);
  return head->val;
}

template <class T>
inline
T&
PP_List<T>::lastElement () const
{
  CH_assert(tail != 0);
  return tail->val;
}

template <class T>
inline
bool
PP_List<T>::isEmpty () const
{
  return head == 0 && tail == 0;
}

template <class T>
inline
bool
PP_List<T>::isNotEmpty () const
{
  return !isEmpty();
}

template <class T>
inline
void
PP_List<T>::removeFirst ()
{
  remove(head);
}

template <class T>
inline
void
PP_List<T>::removeLast ()
{
  remove(tail);
}

template <class T>
inline
const T&
PP_List<T>::operator[] (const PP_ListIterator<T>& li) const
{
  CH_assert(li.p != 0);
  return li.p->val;
}

template <class T>
inline
T&
PP_List<T>::operator[] (const PP_ListIterator<T>& li)
{
  CH_assert(li.p != 0);
  return li.p->val;
}

template <class T>
inline
void
PP_List<T>::replace (PP_ListIterator<T>& li,
                     const T&         _val)
{
  CH_assert(li.p != 0);
  li.p->val = _val;
}

template <class T>
inline
void
PP_List<T>::addAfter (PP_ListIterator<T>& lit,
                      const T&         val)
{
  addAfter(lit.p, val);
}

template <class T>
inline
void
PP_List<T>::addBefore (PP_ListIterator<T>& lit,
                       const T&         val)
{
  addBefore(lit.p, val);
}

template <class T>
inline
PP_ListIterator<T>
PP_List<T>::first () const
{
  return PP_ListIterator<T>(*this,head);
}

template <class T>
inline
PP_ListIterator<T>
PP_List<T>::last () const
{
  return PP_ListIterator<T>(*this,tail);
}

//
// List members
//

template <class T>
inline
PP_List<T>::PP_List (const PP_List<T>& source)
  :
  head(0),
  tail(0)
{
  if (source.isEmpty())
    tail = head = 0;
  else
    for (PP_ListIterator<T> li(source); li; ++li)
      append(li());
}

//
// This isn't inlined as it's declared virtual.
//

template <class T>
inline void
PP_List<T>::add (const T& value)
{
  append(value);
}

template <class T>
inline
int
PP_List<T>::length () const
{
  int len = 0;
  for (PP_ListIterator<T> li(*this); li; ++li)
    len++;
  return len;
}

template <class T>
inline
PP_List<T>&
PP_List<T>::operator= (const PP_List<T>& source)
{
  if (!(this == &source))
    {
      clear();
      for (PP_ListIterator<T> li(source); li; ++li)
        append(li());
    }
  return *this;
}

template <class T>
inline PP_ListLink<T> *
PP_List<T>::addBefore (PP_ListLink<T>* ln,
                       const T&     val)
{
  CH_assert(ln != 0 || head == 0);

  PP_ListLink<T>* newlink;

  if (ln == head)
    {
      head = newlink = new PP_ListLink<T>(val, 0, head);

      if (tail == 0)
        tail = head;
      else
        head->suc->pre = newlink;
    }
  else
    {
      newlink = new PP_ListLink<T>(val, ln->pre, ln);

      ln->pre->suc = newlink;
      ln->pre = newlink;
    }

  return newlink;
}

template <class T>
inline
PP_ListLink<T>*
PP_List<T>::addAfter (PP_ListLink<T>* ln,
                      const T&     val)
{
  CH_assert(ln != 0 || tail == 0);

  PP_ListLink<T>* newlink;

  // trying this here to satisfy icc -O2 (ndk)
  // this worked for me -- but since we aren't actively using icc, i'll
  // leave it commented.
  //newlink = new PP_ListLink<T>(val,tail,0);
  //delete newlink;

  if (ln == tail)
    {
      tail = newlink = new PP_ListLink<T>(val,tail,0);

      if (head == 0)
        head = tail;
      else
        tail->pre->suc = newlink;
    }
  else
    {
      newlink = new PP_ListLink<T>(val, ln, ln->suc);

      ln->suc->pre = newlink;
      ln->suc = newlink;
    }

  return newlink;
}

template <class T>
inline void
PP_List<T>::join (const PP_List<T>& list2)
{
  for (PP_ListIterator<T> li2(list2); li2; ++li2)
    append(li2());
}

template <class T>
inline void
PP_List<T>::catenate (PP_List<T>& list2)
{
  if (list2.isEmpty())
    //
    // Do nothing.
    //
    ;
  else if (isEmpty())
    {
      head = list2.head;
      tail = list2.tail;
      list2.head = 0;
      list2.tail = 0;
    }
  else
    {
      tail->suc = list2.head;
      list2.head->pre = tail;
      tail = list2.tail;
      list2.head = 0;
      list2.tail = 0;
    }
}

template <class T>
inline void
PP_List<T>::clear ()
{
  PP_ListLink<T>* next = 0;

  for (PP_ListLink<T>* p = head; p != 0; p = next)
    {
      next = p->suc;
      p->suc = 0;
      delete p;
    }
  tail = head = 0;
}

template <class T>
inline bool
PP_List<T>::includes (const T& v) const
{
  bool rc = false;
  for (PP_ListIterator<T> li(*this); li && !rc; ++li)
    if (v == li())
      rc = true;
  return rc;
}

template<class T>
inline bool
PP_List<T>::operator== (const PP_List<T>& rhs) const
{
  if (length() == rhs.length())
    {
      for (PP_ListIterator<T> li(*this), ri(rhs); li; ++li, ++ri)
        if (li() != ri())
          return false;
      return true;
    }

  return false;
}

template<class T>
inline bool
PP_List<T>::operator!= (const PP_List<T>& rhs) const
{
  return !operator==(rhs);
}

template <class T>
inline void
PP_List<T>::remove (PP_ListIterator<T>& li)
{
  PP_ListLink<T> *np = li.p->suc;
  remove(li.p);
  li.p = np;
}

template <class T>
inline void
PP_List<T>::remove (const T& _v)
{
  for (PP_ListIterator<T> litr(*this); litr; ++litr)
    if (litr() == _v)
      remove(litr);
}

template <class T>
inline void
PP_List<T>::remove (const PP_List<T>& _lv)
{
  for (PP_ListIterator<T> litr(_lv); litr; ++litr)
    remove(litr());
}

template <class T>
inline void
PP_List<T>::remove (PP_ListLink<T>* ln)
{
  CH_assert(head !=0 && tail != 0);

  if (head == ln && tail == ln)
    head = tail = 0;
  else if (head == ln)
    {
      CH_assert(ln->pre == 0);
      head = ln->suc;
      head->pre = 0;
    }
  else if (tail == ln)
    {
      CH_assert(ln->suc == 0);
      tail = ln->pre;
      tail->suc  = 0;
    }
  else
    {
      CH_assert(ln->suc != 0 && ln->pre != 0);
      ln->suc->pre = ln->pre;
      ln->pre->suc = ln->suc;
    }
  delete ln;
  ln = 0;
}

template <class T> class PP_Array;

/// An Array of Objects of Type T for ParmParse
/**

This class implements an array of objects of the parameterized type
T.  In contrast with the predefined C++ array type, an `PP_Array<T>'
object knows its length, can be dynamically resized, and provides
automatic bounds checking.  The bounds checking can be turned off by
specifying the -DNDEBUG flag on the command line when compiling the
BOXLIB library.  The main reason for defining the PP_Array class is that
it is used, either by composition or inheritance, to implement many of
the other classes in the BOXLIB library.

The `PP_Array<T>' class works by storing copies of the objects it
contains.  If the objects are large, such as `FARRAYBOX's it is
probably better to use the `PArray' class which is an array class that
stores pointers to the objects (avoiding expensive copies).
The `PP_Array<T>' class destructs all objects in the array when it is
itself destructed.  If this is not the desired action, the `PArray' class
should be used.

In the PP_Array<T> class, there are two different concepts of size: the amount
of space allocated, and the amount of space actually used.  Obviously, the
allocated space must be larger than the used space.  We separate these two
concepts to allow the user to optionally avoid memory allocation costs.
Rather than construct and destroy a temporary PP_Array<T> many times, it may
be less expensive to allocate it once with enough memory to handle all uses,
and resize the PP_Array<T> for each particular use.  See the member functions
`reserve', `shrinkWrap', and `resize'.

Note that care must be taken when deriving classes from `PP_Array<T>'.  It is
a concrete class, not a polymorphic one.

This class does NOT provide an assignment operator for assigning an integer
to an PP_Array<T>.

This is a convenience class  for ParmParse
and will not be in any way supported by anyone at ANAG.
*/

template <class T>
class PP_Array
{
public:
  //
  ///: Constructs an `PP_Array<T>' with no elements
  //
  PP_Array ();

  /**: Constructs an `PP_Array<T>' of length len with the value of each
    element defined by the default constructor for `T'.
  */
  /*explicit*/ PP_Array (long len);

  /**: Constructs an `PP_Array<T>' of length len with the value of each
    elements given by initialvalue.
  */
  PP_Array (long     len,
            const T& initialvalue);

  /**: Constructs an `PP_Array<T>' of length len in which the K'th
    value is a copy of vec[K].
  */
  PP_Array (const T* vec,
            long     len);
  //
  ///: The copy constructor.
  //
  PP_Array (const PP_Array<T>& rhs);

  /**: This operator deletes the current `PP_Array<T>' and replaces it
    with a copy of rhs.
  */
  PP_Array<T>& operator= (const PP_Array<T>& rhs);
  //
  ///: Deletes the `PP_Array<T>' and all objects it contains.
  //
  ~PP_Array ();

  /**: Destructs each element in this `PP_Array<T>'.  The resulting array
    has length zero.
  */
  void clear ();
  //
  ///: Returns true if this `PP_Array<T>' is not empty.
  //
  bool ready  () const;

  /**: Reserve space for future expansion of memory.  You still must
    `resize' to address the memory.
  */
  void reserve (long _truesize);

  /**: Shrink allocated memory to be just enough for elements in
    PP_Array.  This is useful if you allocated a lot of memory to
    avoid memory allocation delays as you add elements.  Once all
    the elements are added, you can reduce memory usage to
    the minimum by calling this function.
  */
  void shrinkWrap ();

  /**: This function changes the size of this `PP_Array<T>' to newlen
    while preserving the value of as many array elements as
    possible.  If newlen is greater than the current length, the
    array is grown and the new elements have the value given by
    the default constructor for `T'.  If newlen is less than the
    current length the array is cropped with the remaining
    elements retaining their original values.
  */
  void resize (long newlen);

  /**: This function changes the size of this `PP_Array<T>' to newlen
    while preserving the value of as many array elements as
    possible.  If newlen is greater than the current length, the
    array is grown and the new elements have the value given by
    initialvalue.  If newlen is less than the current length the
    array is cropped with the remaining elements retaining their
    original values.
  */
  void resize (long     newlen,
               const T& initialvalue);
  //
  ///: Return number of elements in the array.
  //
  long length () const;

  /**: Return the maximum number of elements the array can hold
    without doing a `resize'.
  */
  long trueSize () const;

  /**: Returns a reference to the K'th element in this `PP_Array<T>'.
    The element can be modified through this reference.  The
    result can be used as an L-value.
  */
  T& operator[] (long K);
  //
  ///: Same as above, except acts on const PP_Array's.
  //
  const T& operator[] (long K) const;
  //
  ///: Different syntax for operator[] (long i).
  //
  T& get (long i);
  //
  ///: Different syntax for const operator[] (long i).
  //
  const T& get (long i) const;

  /**: Returns pointer to vector of data.  This function breaks object
    encapsulation and should only be used for interfacing to
    Fortran subroutines.
  */
  T* dataPtr ();
  //
  ///: Same as above for constant PP_Arrays.
  //
  const T* dataPtr () const;
  //
  ///: Changes the i'th element of this `PP_Array<T>' to elem.
  //
  void set (long     i,
            const T& elem);
  //
  ///: This function swaps the i'th and j'th element of the PP_Array.
  //
  void swap (long i,
             long j);
  //
  ///: Test for equality.
  //
  bool operator== (const PP_Array<T>& rhs) const;
  //
  ///: Test for inequality.
  //
  bool operator!= (const PP_Array<T>& rhs) const;

protected:
  //
  // The true size of the PP_Array.
  //
  long truesize;
  //
  // The number of elements in the PP_Array.
  //
  long nelem;
  //
  // The array itself.
  //
  T* vp;

private:
  //
  // This is disallowed.
  //
  PP_Array<T>& operator= (int);
};

//
// Inlines.
//

template <class T>
inline
PP_Array<T>::PP_Array ()
{
  nelem    = 0;
  vp       = new T[1];
  truesize = 1;
}

template <class T>
inline
PP_Array<T>::PP_Array (long len)
{
  CH_assert(len >= 0);
  nelem    = len;
  vp       = new T[len];
  truesize = nelem;
}

template <class T>
inline
void
PP_Array<T>::clear ()
{
  delete [] vp;
  vp       = 0;
  nelem    = 0;
  truesize = 0;
}

template <class T>
inline
PP_Array<T>::~PP_Array ()
{
  clear();
}

template <class T>
inline
bool
PP_Array<T>::ready () const
{
  return vp != 0 && nelem != 0;
}

template <class T>
inline
long
PP_Array<T>::length () const
{
  return nelem;
}

template <class T>
inline
long
PP_Array<T>::trueSize () const
{
  return truesize;
}

template <class T>
inline
T&
PP_Array<T>::operator[] (long i)
{
  CH_assert(vp != 0);
  CH_assert(i >= 0 && i < nelem);
  return vp[i];
}

template <class T>
inline
const T&
PP_Array<T>::operator[] (long i) const
{
  CH_assert(vp != 0);
  CH_assert(i >= 0 && i < nelem);
  return vp[i];
}

template <class T>
inline
T&
PP_Array<T>::get (long i)
{
  CH_assert(vp != 0);
  CH_assert(i >= 0 && i < nelem);
  return vp[i];
}

template <class T>
inline
const T&
PP_Array<T>::get (long i) const
{
  CH_assert(vp != 0);
  CH_assert(i >= 0 && i < nelem);
  return vp[i];
}

template <class T>
inline
void
PP_Array<T>::set (long     i,
                  const T& elem)
{
  CH_assert(vp != 0);
  CH_assert(i >= 0 && i < nelem);
  vp[i] = elem;
}

template <class T>
inline
T*
PP_Array<T>::dataPtr ()
{
  return vp;
}

template <class T>
inline
const T*
PP_Array<T>::dataPtr () const
{
  return vp;
}

template <class T>
inline
void
PP_Array<T>::swap (long i,
                   long j)
{
  CH_assert(i >= 0 && i < nelem);
  CH_assert(j >= 0 && j < nelem);
  T tmp = vp[i];
  vp[i] = vp[j];
  vp[j] = tmp;
}

template <class T>
inline
bool
PP_Array<T>::operator!= (const PP_Array<T>& rhs) const
{
  return !(operator==(rhs));
}

//
// Non-inlined stuff.
//

template <class T>
PP_Array<T>::PP_Array (long     len,
                       const T& initialValue)
{
  CH_assert(len >= 0);
  nelem = len;
  vp    = new T[len];
  truesize = nelem;
  for (long i = 0; i < nelem; ++i)
    vp[i] = initialValue;
}

template <class T>
PP_Array<T>::PP_Array (const T* vec,
                       long     len)
{
  CH_assert(len >= 0);
  nelem = len;
  vp    = new T[len];
  truesize = nelem;
  for (long i = 0; i < nelem; ++i)
    vp[i] = vec[i];
}

template <class T>
PP_Array<T>::PP_Array (const PP_Array<T>& a)
{
  nelem = a.nelem;
  vp    = new T[nelem];
  truesize = nelem;
  for (long i = 0; i < nelem; i++)
    vp[i] = a.vp[i];
}

template <class T>
PP_Array<T>&
PP_Array<T>::operator= (const PP_Array<T>& sa)
{
  if (this != &sa)
    {
      clear();
      vp       = new T[sa.nelem];
      nelem    = sa.nelem;
      truesize = nelem;
      for (long i = 0; i < nelem; i++)
        vp[i] = sa.vp[i];
    }
  return *this;
}

template <class T>
inline
void
PP_Array<T>::resize (long newlen)
{
  if (newlen == nelem)
    return;
  if (newlen <= truesize)
    {
      nelem = newlen;
      return;
    }
  T* newvp = new T[newlen];
  long len = Min(newlen,nelem);
  for (long i = 0; i < len; i++)
    newvp[i] = vp[i];
  delete [] vp;
  vp = newvp;
  nelem = newlen;
  truesize = newlen;
}

template <class T>
inline
void PP_Array<T>::resize (long     newlen,
                          const T& initialValue)
{
  if (newlen == nelem)
    return;
  if (newlen <= truesize)
    {
      for (long i = nelem; i < newlen; ++i)
        vp[i] = initialValue;
      nelem = newlen;
      return;
    }
  T* newvp = new T[newlen];
  long len = Min(newlen,nelem);
  long i;
  for (i = 0; i < len; i++)
    newvp[i] = vp[i];
  for (i = len; i < newlen; ++i)
    newvp[i] = initialValue;
  delete [] vp;
  vp = newvp;
  nelem = newlen;
  truesize = newlen;
}

template <class T>
void
PP_Array<T>::reserve (long _truesize)
{
  if (_truesize > truesize)
    {
      T* newvp = new T[_truesize];
      for (long i = 0; i < nelem; i++)
        newvp[i] = vp[i];
      delete [] vp;
      vp = newvp;
      truesize = _truesize;
    }
}

template <class T>
void
PP_Array<T>::shrinkWrap ()
{
  if (nelem != truesize)
    {
      T* newvp = new T[nelem];
      for (long i = 0; i < nelem; i++)
        newvp[i] = vp[i];
      delete [] vp;
      vp = newvp;
      truesize = nelem;
    }
}

template <class T>
bool
PP_Array<T>::operator== (const PP_Array<T>& rhs) const
{
  if (length() != rhs.length())
    return false;

  for (long i = 0; i < length(); ++i)
    if (!((*this)[i] == rhs[i]))
      return false;

  return true;
}

// -----------------------------------------------------------------
// ----------------------- COMMENTS -------------------------------
// -----------------------------------------------------------------
// The ParmParse class implements a simple database for the storage
// and retrieval of command-line and input-file arguments.  The
// entries are stored in a static table in (name,value_list) pairs.
//
// The format of the input file is a series of OPTIONS and DEFINITIONS.
//
// An OPTION is an entry of the form:  -<name> and has no associated list
// of values.  For example, the command line:
//        prog -verbose -no_opt
// has two options: "verbose" and "no_opt".
//
// A DEFINITION is of the form  <name> = <value> <value> ...
// The equal sign is important since the list of values can span multiple
// lines.
//
// Comments in an input file include all text from a '#' character to the
// end of the line.  Here is an example input file:
//
//   -no_garbage                # an OPTION
//   niter = 100                # niter is an integer
//   title = "Double Wammy"     # example of a string with spaces
//   cell_size = 0.5 0.75       # cell spacing in each dimension
//   plot.var = Density 1 10    # a list of values
//   plot.var = Energy  5 12    # another list of values
//   bigarray = 1 2 3 4 5 6 7 8 # first part of array
//              9 10 11 12      # continuation of bigarray
//   test = apple "boy blue" 10 20 30 40
//   FILE = prob_file           # insert contents of this "prob_file" here
//
// The "FILE = <filename>" definition is special.  Rather than just
// adding this entry to the database, it reads the contents of <filename>
// into the database.
//
// ParmParse stores all entries in a static table which is built the
// first time a ParmParse object is constructed (usually in main()).
// Subsequent invocations have access to this table.
// A ParmParse constructor has an optional "prefix" argument that will
// limit the searches to only those entries of the table with this prefix
// in name.  For example:
//     ParmParse pp("plot");
// will find only those entries with name given by "plot.<string>".
//
// All values in the table are stored as strings.  For example, the
// values of "cell_size" in the above input file are stored as the
// strings "0.5" and "0.75".  These strings can be returned as either
// string of numeric values by the query functions.
// Character strings with spaces must be delimited by double quotes
// in the input file but the quotes are stripped before they are entered
// into the table.  For example, 'title' in the above input file has a
// single value, the string 'Double Wammy' (without the quotes).
// Each value in the list associated with a definition can be referred to
// by its index number.  The index numbers start at 0 just like an array
// in the C programming language.  Consider the definition of "test" in
// the above input file.  The first value 'apple'is a string with index
// 0.  The second value 'boy blue' is a string with index 1.  The
// remaining four values are integers indexed 2, 3, 4, and 5.
//
// For a string value to represent an integer or float it must fit the
// following regular experssion:
//   Sign    ::= '+' | '-'
//   Digit   ::= '0' | '1' | ... | '9'
//   Integer ::= [Sign]Digit+
//   Exp     ::= ('e'|'E')Integer
//   Float   ::= ( Integer[.Digit*][Exp] | [Integer].Digit+[Exp] )
//
// Where '+' indicates one or more occurences, '*' represents zero or
// more occurences, '|' means one or the other and '[]' represents zero
// or one occurence.
//
// Note that floats and doubles have the same string representation and
// that the FORTRAN "double" exponent format is not supported.
// That is, 1.0d+3 is not a valid representation of a floating point
// number but that 1.0e+3 is acceptable.
//
// There are a host of functions allowing the user to query the database
// and retrieve values.  Here are some general rules about the names of
// the member functions:
//
// * Functions with the string "get" in their names attempt to get a
//   value or an array of values from the table.  They generate a
//   run-time error if they are not successful.
//
// * Functions with the string "query" in their names attempt to get a
//   value or an array of values from the table.  They return the value 1
//   (true) if they are successful and 0 (false) if not.
//
// * Functions with the string "arr" in their names get an Array of
//   values from the given entry in the table.  The array argument is
//   resized (if necessary) to hold all the values requested.
//
// * Functions without the string "arr" in their names get single
//   values from the given entry in the table.
//
// The following is a code sample showing how to use ParmParse:
//
// main(int argc, char **argv)
// {
//     char* in_file_name = argv[1];
//     ParmParse pp(argc-2, argv+2, 0, in_file_name);
//
//     // was the "-verbose" command line argument set?
//     int verbose = pp.contains("verbose");
//
//     // Query table for value of "niter".  If not in table
//     // then set to default value
//     if (!pp.query("niter",niter)) niter = 20;
//
//     // read array of cell sizes if in table
//     std::vector<float> dx;
//     if (nx=pp.countval("cell_size")) {
//        // get nx values starting at index 0 and store in dx.
//        // dx is automatically resized here.
//        pp.getarr("cell_size",dx,0,nx);
//     }
// }
//
// void do_graphics()
// {
//    //
//    // Will only query entries with the "plot" prefix:
//    //
//    ParmParse pp("plot");
//    //
//    // Read all variables with "plot.var" keyword.
//    //
//    std::string var_name;
//    std::vector<int> range;
//    int num = pp.countname("var");
//    //
//    // Element 0 in list is a string.
//    //
//    pp.get("var",var_name,0);
//    //
//    // Elements 1 and 2 are integers.
//    // Note that "range" will be resized to hold 2 elements.
//    //
//    pp.getarr("var",range,1,2);
//    cout << "variable = " << var_name << "lo, hi = ",
//         << range[0] << " " << range[1] << endl;
// }
// -----------------------------------------------------------------
// -----------------------  END COMMENTS ---------------------------
// -----------------------------------------------------------------

//
// Forward reference to private class.
//
class PP_entry;

#endif //for the doxygen thing

/// Parse Parameters From Command Line and Input Files
/**

The ParmParse class is used to interpret parameters passed in to a program
from the command line and an arbitrary collection of input files.  The
parameters are stored in static table that can be queried by any object
of type ParmParse.  A parameter can be either an "option" (usually
specified on the command line) or a "definition".  An option is of the form
"-<name>" and is stored in the table without the hyphen.  A definition is
of the form "<name> = <value><value>...<value>".  It is stored in the table
as a name, value-list pair.

In the following example, verbose and no_opt are stored in the table as
options.  niter is a definition with the single integer value 10; name is
a definition with the string value "big code" and dx is a definition with
the two floating point values 0.5 and 0.75.

prog -verbose -no_opt niter = 10 name = "big code" dx = 0.5 0.75

The ParmParse class has two constructors.  The first is responsible for
building the table and is usually called by the main routine of an
application.  It has arguments for the command line argc and argv parameters,
as well as an optional filename argument for reading definitions from an
input file.  The table is built by reading the input file first (if it
exists) with the command line arguments added to the end of the table.
The order of a definition in the table is significant, so command line
parameters can be used to override definitions in the input file.  A
definition of the explicit form: FILE=<filename> is not added to the table
but is a directive to include the named file at that point in the table.

The second constructor is generally used by other classes in the code.  It
permits access to the table via a large collection of query functions.
Both constructors have an optional prefix argument that narrows the search
to entries in the table with the same prefix.  For example, let PlanR be a
ParmParse object with code prefix "ope".  PlanR.get("val",v) will look for
an entry in the parameter list of the form: ope.val==<value>, and will
reject all entries not starting with the correct code prefix.

The query functions search the table for definition names that match a given
string (and prefix) and return values from the corresponding value list.
The values can be returned as ints, Array<int>s, floats, std::vector<float>s,
doubles, std::vector<double>s, std::strings, or std::vector<std::sring>s.  All values in the
table are stored as PP_String objects, but if an int, float, or double is
requested, the translation is done automatically.  In the previous example,
the value of niter could be returned as either an std::string, an int, a double,
or a float.  The values of dx can be returned as std::strings, floats, or
doubles, but the value of name can be returned only as an std::string.

Comments in an input file include all text from a `\#' character to the
end of the line.  Here is a sample input file:

-no_garbage

niter = 100

title = "Double Wammy"

cell_size = 0.5 0.75

plot.var = Density 1 10

plot.var = Energy  5 12

bigarray = 1 2 3 4 5 6 7 8

9 10 11 12

test = apple "boy blue" 10 20 30 40

FILE = prob_file
*/

class ParmParse
{
  friend class PP_entry;
public:

  ///
  /** Construct an initial ParmParse object from the argc and argv
    passed in to main().  An error will be signalled if another
    ParmParse object currently exists.  If parfile is specified,
    read the parameters in from that file first and then append
    those derived from argv to the table.  If prefix is specified,
    load this string as the code prefix for this particular
    ParmParse object.
  */
  ParmParse (int         argc,
             char**      argv,
             const char* prefix  = 0,
             const char* parfile = 0);

  /// Initialize after calling the default constructor.
  /** Same as the constructor with the same args.
   */
  void define (int         argc,
               char**      argv,
               const char* prefix  = 0,
               const char* parfile = 0);

  /** Construct an additional ParmParse object sharing the same
    internal table as any other such objects in existence.  If
    prefix is specified, load this string as the code prefix
    for this particular ParmParse object.  If there are no other
    existing objects, this doesn't do anything, and define() must
    be called before the object can be used.
  */
  /*explicit*/ ParmParse (const char* prefix = 0);

  /** The destructor.  The internal static table will only be deleted
    if there are no other ParmParse objects in existence.
  */
  ~ParmParse();

  ///
  /**
     Returns true if name is in table.
  */
  bool contains (const char* name) const;

  ///
  /**
     Returns true if name is in table.
  */
  bool contains (const std::string& name) const;

  ///
  /** Returns the number of values associated with nth occurence of
    name (prepended with the prefix) in the table.  n == -1 implies
    the last occurence.
  */
  int countval (const char* name,
                int         n = -1) const;

  ///
  /** Returns the number of times the given name (prepended with
    prefix) appears in the table.
  */
  int countname (const char* name) const;

  ///
  /** Returns the number of times the given name (prepended with
    prefix) appears in the table.
  */
  int countname (const std::string& name) const;

  ///
  /**
   * Accessor for private method that adds parameters
   * to the global table (SWS-Tech-X Research)
   * e.g.   parmParseObj.addEntries("tst_value = 2")
   */
  void addEntries (const std::string& strEntry) 
  {
    int len = (int)strEntry.size() + 1;
    bldTable(strEntry.c_str(), len, table);
  }

  ///
  /**
   * Append input parameter pairs to ParmParse
   *
   *  1. convert value to string
   *  2. construct string line syntax for ParmParse
   *  3. use special public method in ParmParse.H
   *
   * @param varName parameter string name
   * @param varVal numerical value type
   */
  template <typename TYPE>
  void setVal (const std::string& varName, TYPE varVal) 
  {
    std::string valStr =
      static_cast<std::ostringstream*>( &(std::ostringstream() << varVal) )->str();
    std::string inEntry = varName + " = " + valStr;
    addEntries(inEntry);
  }

  ///
  /**
   * Append input parameter pairs to ParmParse
   *
   *  1. convert value to string
   *  2. construct string line syntax for ParmParse
   *  3. use special public method in ParmParse.H
   *
   * @param varName parameter string name
   * @param varVal string value type
   */
  void setStr (const std::string& varName, std::string varVal) 
  {
    std::string inEntry = varName + " = " + varVal;
    addEntries(inEntry);
  }

  ///
  /**
     Write the contents of the table in ASCII to the ostream.
   */
  void dumpTable (std::ostream& os) const;

  /// access single object
  /** Get the ival'th value of last occurrence of the requested name.
    If successful, the value is converted to an int and stored
    in reference ref.  If the name does not exist or
    ival'th value does not exist, or if the printed representation
    of the value cannot be converted to an int, an error message is
    output and the program halts.   Note that ival == 0 is the first
    value in the list.
  */
  void get (const char* name,
            int&        ref,
            int         ival=0) const;

  ///
  /** Get the ival'th value of last occurrence of the
    requested name.  If successful, the value is converted
    to an int and stored in reference ref.  Returns 1 if
    successful.  Returns 0 if the name does not exist.  If
    ival'th value does not exist, or if the printed
    representation of the value cannot be converted to an
    int, an error message is output and the program halts.
    Note that ival == 0 is the first value in the list.
  */
  int query (const char* name,
             int&        ref,
             int         ival=0) const;

  /// access single object
  /** Get the ival'th value of last occurrence of the requested name.
    If successful, the value is converted to a n unsigned long and stored
    in reference ref.  If the name does not exist or
    ival'th value does not exist, or if the printed representation
    of the value cannot be converted to an int, an error message is
    output and the program halts.   Note that ival == 0 is the first
    value in the list.
  */
  void get (const char* name,
            unsigned long& ref,
            int         ival=0) const;

  ///
  /** Get the ival'th value of last occurrence of the
    requested name.  If successful, the value is converted
    to an unsigned long and stored in reference ref.  Returns 1 if
    successful.  Returns 0 if the name does not exist.  If
    ival'th value does not exist, or if the printed
    representation of the value cannot be converted to an
    int, an error message is output and the program halts.
    Note that ival == 0 is the first value in the list.
  */
  int query (const char* name,
             unsigned long& ref,
             int         ival=0) const;

  ///
  /** Get the ival'th value of last occurrence of the requested name.
    If successful, the value is converted to a float and stored
    in reference ref.  If the name does not exist or
    ival'th value does not exist, or if the printed representation
    of the value cannot be converted to a float, an error message is
    output and the program halts.   Note that ival == 0 is the first
    value in the list.
  */
  void get (const char* name,
            float&      ref,
            int         ival=0) const;

  ///
  /** Get the ival'th value of last occurrence of the
    requested name.  If successful, the value is converted
    to a float and stored in reference ref.  Returns 1 if
    successful.  Returns 0 if the name does not exist.  If
    ival'th value does not exist, or if the printed
    representation of the value cannot be converted to a float,
    an error message is output and the program halts.
    Note that ival == 0 is the first value in the list.
  */
  int query (const char* name,
             float&      ref,
             int         ival=0) const;

  ///
  /** Get the ival'th value of last occurrence of the requested name.
    If successful, the value is converted to a double and stored
    in reference ref.  If the name does not exist or
    ival'th value does not exist, or if the printed representation
    of the value cannot be converted to a double, an error message is
    output and the program halts.   Note that ival == 0 is the first
    value in the list.
  */
  void get (const char* name,
            double&     ref,
            int         ival=0) const;

  ///
  /** Get the ival'th value of last occurrence of the
    requested name.  If successful, the value is converted
    to a double and stored in reference ref.  Returns 1 if
    successful.  Returns 0 if the name does not exist.  If
    ival'th value does not exist, or if the printed
    representation of the value cannot be converted to a double,
    an error message is output and the program halts.
    Note that ival == 0 is the first value in the list.
  */
  int query (const char* name,
             double&     ref,
             int         ival=0) const;

  ///Get the ival'th value of last occurrence of the requested name.
  /** Get the ival'th value of last occurrence of the requested name.
    If successful, the value is converted to a string and stored
    in reference ref.  If the name does not exist or
    ival'th value does not exist, or if the printed representation
    of the value cannot be converted to a string, an error message is
    output and the program halts.   Note that ival == 0 is the first
    value in the list.
  */
  void get (const char*  name,
            std::string& ref,
            int          ival=0) const;

  ///Get the ival'th value of last occurrence of the requested name.
  /** Get the ival'th value of last occurrence of the
    requested name.  If successful, the value is converted
    to a string and stored in reference ref.  Returns 1 if
    successful.  Returns 0 if the name does not exist.  If
    ival'th value does not exist, or if the printed
    representation of the value cannot be converted to a string,
    an error message is output and the program halts.
    Note that ival == 0 is the first value in the list.
  */
  int query (const char*  name,
             std::string& ref,
             int          ival=0) const;

  /// Get the ival'th value of last occurrence of the requested name.
  /** Get the ival'th value of last occurrence of the requested name.
    If successful, the value is converted to a bool and stored
    in reference ref.  If the name does not exist or
    ival'th value does not exist, or if the printed representation
    of the value cannot be converted to a bool, an error message is
    output and the program halts.   Valid representations of bool
    are "true" and "false" (case is ignored), and numeric values
    1 (=true) and 0 (=false).  Note that ival == 0 is the first
    value in the list.
  */
  void get (const char* name,
            bool&       ref,
            int         ival=0) const;

  ///Get the ival'th value of last occurrence of the requested name.
  /** Get the ival'th value of last occurrence of the
    requested name.  If successful, the value is converted
    to a bool and stored in reference ref.  Returns 1 if
    successful.  Returns 0 if the name does not exist.  If
    ival'th value does not exist, or if the printed
    representation of the value cannot be converted to a bool,
    an error message is output and the program halts.
    Valid representations of bool are "true" and "false"
    (case is ignored), and numeric values 1 (=true) and 0 (=false).
    Note that ival == 0 is the first value in the list.
  */
  int query (const char* name,
             bool&       ref,
             int         ival=0) const;

  /// access an array of objects
  /** Gets a std::vector<int> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to an int and stored in the
    std::vector<int> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<int>[0], std::vector<int>[1] holds
    start_ix+1, etc.  If there are fewer than start_ix +
    num_val values associated with the last occurrence, or
    if some of the values cannot be converted to an int, an
    error message is reported and the program halts.
  */
  void getarr (const char* name,
               Vector<int>& ref,
               int         start_ix,
               int         num_val) const;

  /// access an array of objects
  /** Gets a std::vector<int> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to an int and stored in the
    std::vector<int> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<int>[0], std::vector<int>[1] holds
    start_ix+1, etc.  If there are fewer than start_ix +
    num_val values associated with the last occurrence, or
    if some of the values cannot be converted to an int, an
    error message is reported and the program halts.
  */
  void getarr (const char* name,
               std::vector<int>& ref,
               int         start_ix,
               int         num_val) const;

  /// access an array of objects
  /** Gets a std::vector<int> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to an int and stored in the
    std::vector<int> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<int>[0], std::vector<int>[1] holds
    start_ix+1, etc.  Returns 0 if the name does not
    exist.  If there are fewer than start_ix + num_val
    values associated with the last occurrence, or if some
    of the values cannot be converted to an int, an error
    message is reported and the program halts.
  */
  int queryarr (const char* name,
                Vector<int>& ref,
                int         start_ix,
                int         num_val) const;

  /// access an array of objects
  /** Gets a std::vector<int> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to an int and stored in the
    std::vector<int> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<int>[0], std::vector<int>[1] holds
    start_ix+1, etc.  Returns 0 if the name does not
    exist.  If there are fewer than start_ix + num_val
    values associated with the last occurrence, or if some
    of the values cannot be converted to an int, an error
    message is reported and the program halts.
  */
  int queryarr (const char* name,
                std::vector<int>& ref,
                int         start_ix,
                int         num_val) const;

  /// access an array of objects
  /** Gets a std::vector<unsigned long> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to an int and stored in the
    std::vector<unsigned long> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<unsigned long>[0], std::vector<unsigned long>[1] holds
    start_ix+1, etc.  If there are fewer than start_ix +
    num_val values associated with the last occurrence, or
    if some of the values cannot be converted to an unsigned long, an
    error message is reported and the program halts.
  */
  void getarr (const char* name,
               Vector<unsigned long>& ref,
               int         start_ix,
               int         num_val) const;

  /// access an array of objects
  /** Gets a std::vector<unsigned long> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to an unsigned long and stored in the
    std::vector<unsigned long> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<unsigned long>[0], std::vector<unsigned long>[1] holds
    start_ix+1, etc.  If there are fewer than start_ix +
    num_val values associated with the last occurrence, or
    if some of the values cannot be converted to an unsigned long, an
    error message is reported and the program halts.
  */
  void getarr (const char* name,
               std::vector<unsigned long>& ref,
               int         start_ix,
               int         num_val) const;

  /// access an array of objects
  /** Gets a std::vector<unsigned long> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to an unsigned long and stored in the
    std::vector<unsigned long> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<unsigned long>[0], std::vector<unsigned long>[1] holds
    start_ix+1, etc.  Returns 0 if the name does not
    exist.  If there are fewer than start_ix + num_val
    values associated with the last occurrence, or if some
    of the values cannot be converted to an unsigned long, an error
    message is reported and the program halts.
  */
  int queryarr (const char* name,
                Vector<unsigned long>& ref,
                int         start_ix,
                int         num_val) const;

  /// access an array of objects
  /** Gets a std::vector<unsigned long> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to an unsigned long and stored in the
    std::vector<unsigned long> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<unsigned long>[0], std::vector<unsigned long>[1] holds
    start_ix+1, etc.  Returns 0 if the name does not
    exist.  If there are fewer than start_ix + num_val
    values associated with the last occurrence, or if some
    of the values cannot be converted to an unsigned long, an error
    message is reported and the program halts.
  */
  int queryarr (const char* name,
                std::vector<unsigned long>& ref,
                int         start_ix,
                int         num_val) const;

  /// access an array
  /** Gets a std::vector<float> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a float and stored in the
    std::vector<float> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<float>[0], std::vector<float>[1] holds
    start_ix+1, etc.  If there are fewer than start_ix +
    num_val values associated with the last occurrence, or
    if some of the values cannot be converted to a float, an
    error message is reported and the program halts.
  */
  void getarr (const char*   name,
               Vector<float>& ref,
               int           start_ix,
               int           num_val) const;

  /// access an array
  /** Gets a std::vector<float> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a float and stored in the
    std::vector<float> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<float>[0], std::vector<float>[1] holds
    start_ix+1, etc.  If there are fewer than start_ix +
    num_val values associated with the last occurrence, or
    if some of the values cannot be converted to a float, an
    error message is reported and the program halts.
  */
  void getarr (const char*   name,
               std::vector<float>& ref,
               int           start_ix,
               int           num_val) const;

  ///
  /** Gets a std::vector<float> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a float and stored in the
    std::vector<float> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<float>[0], std::vector<float>[1] holds
    start_ix+1, etc.  Returns 0 if the name does not
    exist.  If there are fewer than start_ix + num_val
    values associated with the last occurrence, or if some
    of the values cannot be converted to a float, an error
    message is reported and the program halts.
  */
  int queryarr (const char*   name,
                Vector<float>& ref,
                int           start_ix,
                int           num_val) const;

  ///
  /** Gets a std::vector<float> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a float and stored in the
    std::vector<float> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<float>[0], std::vector<float>[1] holds
    start_ix+1, etc.  Returns 0 if the name does not
    exist.  If there are fewer than start_ix + num_val
    values associated with the last occurrence, or if some
    of the values cannot be converted to a float, an error
    message is reported and the program halts.
  */
  int queryarr (const char*   name,
                std::vector<float>& ref,
                int           start_ix,
                int           num_val) const;

  ///
  /** Gets a std::vector<double> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a double and stored in the
    std::vector<double> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<double>[0], std::vector<double>[1] holds
    start_ix+1, etc.  If there are fewer than start_ix +
    num_val values associated with the last occurrence, or
    if some of the values cannot be converted to a double, an
    error message is reported and the program halts.
  */
  void getarr (const char*    name,
               Vector<double>& ref,
               int            start_ix,
               int            num_val) const;

  ///
  /** Gets a std::vector<double> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a double and stored in the
    std::vector<double> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<double>[0], std::vector<double>[1] holds
    start_ix+1, etc.  If there are fewer than start_ix +
    num_val values associated with the last occurrence, or
    if some of the values cannot be converted to a double, an
    error message is reported and the program halts.
  */
  void getarr (const char*    name,
               std::vector<double>& ref,
               int            start_ix,
               int            num_val) const;

  ///
  /** Gets a std::vector<double> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a double and stored in the
    std::vector<double> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<double>[0], std::vector<double>[1] holds
    start_ix+1, etc.  Returns 0 if the name does not
    exist.  If there are fewer than start_ix + num_val
    values associated with the last occurrence, or if some
    of the values cannot be converted to a double, an error
    message is reported and the program halts.
  */
  int queryarr (const char*    name,
                Vector<double>& ref,
                int            start_ix,
                int            num_val) const;

  ///
  /** Gets a std::vector<double> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a double and stored in the
    std::vector<double> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<double>[0], std::vector<double>[1] holds
    start_ix+1, etc.  Returns 0 if the name does not
    exist.  If there are fewer than start_ix + num_val
    values associated with the last occurrence, or if some
    of the values cannot be converted to a double, an error
    message is reported and the program halts.
  */
  int queryarr (const char*    name,
                std::vector<double>& ref,
                int            start_ix,
                int            num_val) const;

  ///
  /** Gets a std::vector<string> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a string and stored in the
    std::vector<string> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<string>[0], std::vector<string>[1] holds
    start_ix+1, etc.  If there are fewer than start_ix +
    num_val values associated with the last occurrence, or
    if some of the values cannot be converted to a string, an
    error message is reported and the program halts.
  */
  void getarr (const char*     name,
               Vector<std::string>& ref,
               int             start_ix,
               int             num_val) const;

  ///
  /** Gets a std::vector<string> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a string and stored in the
    std::vector<string> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<string>[0], std::vector<string>[1] holds
    start_ix+1, etc.  If there are fewer than start_ix +
    num_val values associated with the last occurrence, or
    if some of the values cannot be converted to a string, an
    error message is reported and the program halts.
  */
  void getarr (const char*     name,
               std::vector<std::string>& ref,
               int             start_ix,
               int             num_val) const;

  ///
  /** Gets a std::vector<string> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a string and stored in the
    std::vector<string> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<string>[0], std::vector<string>[1] holds
    start_ix+1, etc.  Returns 0 if the name does not
    exist.  If there are fewer than start_ix + num_val
    values associated with the last occurrence, or if some
    of the values cannot be converted to a string, an error
    message is reported and the program halts.
  */
  int queryarr (const char*               name,
                Vector<std::string>& ref,
                int                       start_ix,
                int                       num_val) const;

  ///
  /** Gets a std::vector<string> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a string and stored in the
    std::vector<string> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<string>[0], std::vector<string>[1] holds
    start_ix+1, etc.  Returns 0 if the name does not
    exist.  If there are fewer than start_ix + num_val
    values associated with the last occurrence, or if some
    of the values cannot be converted to a string, an error
    message is reported and the program halts.
  */
  int queryarr (const char*               name,
                std::vector<std::string>& ref,
                int                       start_ix,
                int                       num_val) const;

///
  /** Gets a std::vector<bool> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a bool and stored in the
    std::vector<bool> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<bool>[0], std::vector<bool>[1] holds
    start_ix+1, etc.  If there are fewer than start_ix +
    num_val values associated with the last occurrence, or
    if some of the values cannot be converted to a bool, an
    error message is reported and the program halts.
  */
  void getarr (const char*   name,
               Vector<bool>& ref,
               int           start_ix,
               int           num_val) const;

  ///
  /** Gets a std::vector<bool> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a bool and stored in the
    std::vector<bool> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<bool>[0], std::vector<bool>[1] holds
    start_ix+1, etc.  If there are fewer than start_ix +
    num_val values associated with the last occurrence, or
    if some of the values cannot be converted to a bool, an
    error message is reported and the program halts.
  */
  void getarr (const char*        name,
               std::vector<bool>& ref,
               int                start_ix,
               int                num_val) const;

  ///
  /** Gets a std::vector<bool> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a bool and stored in the
    std::vector<bool> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<bool>[0], std::vector<bool>[1] holds
    start_ix+1, etc.  Returns 0 if the name does not
    exist.  If there are fewer than start_ix + num_val
    values associated with the last occurrence, or if some
    of the values cannot be converted to a bool, an error
    message is reported and the program halts.
  */
  int queryarr (const char*   name,
                Vector<bool>& ref,
                int           start_ix,
                int           num_val) const;

  ///
  /** Gets a std::vector<bool> of num_val values from last
    occurrence of given name.  If successful, the values
    are converted to a bool and stored in the
    std::vector<bool> object ref.  ref is resized (if
    necessary) to hold num_val values.  The value in the
    list indexed by start_ix is copied into
    std::vector<bool>[0], std::vector<bool>[1] holds
    start_ix+1, etc.  Returns 0 if the name does not
    exist.  If there are fewer than start_ix + num_val
    values associated with the last occurrence, or if some
    of the values cannot be converted to a bool, an error
    message is reported and the program halts.
  */
  int queryarr (const char*        name,
                std::vector<bool>& ref,
                int                start_ix,
                int                num_val) const;

  /// return prefix
  inline const char* prefix() const
  {
    return thePrefix.c_str();
  }

#ifndef DOXYGEN
  //
  // This should be protected, but cfront-derived compilers complain :-(
  //
  enum PPType
    {
      ppDefn,
      ppOption,
      ppInt,
      ppUnsignedLong,
      ppFloat,
      ppDouble,
      ppString,
      ppBool,
      ppEQ_sign,
      ppEOF
    };


protected:
  //
  // Table of entries common to all objects.
  //
  static PP_List<PP_entry*> table;
  //
  // Command line arguments.
  //
  static int    xargc;
  static char** xargv;
  //
  // Keep track of number of ParmParse objects out there.
  //
  static int num_obj;
  //
  // Parses string and builds table.
  //
  void bldTable (const char*      str,
                 int              lenstr,
                 PP_List<PP_entry*>& tab);
  //
  // Add defn to table, check for file inclusion.
  //
  void addDefn (PP_String&         def,
                PP_List<PP_String>&   val,
                PP_List<PP_entry*>& tab);
  //
  // Reads file into string then parses it with call to bldTable.
  //
  void read_file (const char*      fname,
                  PP_List<PP_entry*>& tab);
  //
  // Lexical analyser called by bldTable.
  //
  PPType getToken (const char*,
                   int&,
                   int,
                   char*);
  //
  // Reclaims memory used in table.
  //
  void rmTable ();
  //
  // Prefix used in keyword search.
  //
  PP_String thePrefix;
  //
  // Used by constructor to build table.
  //
  void ppinit (const char* parfile);
  //
  // Find n'th occurence of name in table.
  //
  const PP_entry* ppindex (int         n,
                           const char* name) const;
  //
  // Get ival_th value of k_th occurence of given name.
  // If k_th occurence does not exist or ival_th value does
  // not exist or if ival_type does not match with type, an
  // error message is reported and the program halts.
  // If successful, value is stored in ptr.
  // same as above but searches for last occurence of name.
  //
  void  getval (const char* name,
                const PPType type,
                void*        ptr,
                int          ival,
                int          k=-1) const;
  //
  // Get an array of values.
  //
  void getarr (const char*  name,
               const PPType type,
               void*        ptr,
               int          start_ix,
               int          num_val,
               int          k=-1) const;
  int queryval (const char*  name,
                const PPType type,
                void*        ptr,
                int          ival,
                int          k=-1) const;
  int queryarr (const char*  name,
                const PPType type,
                void*        ptr,
                int          start_ix,
                int          num_val,
                int          k=-1) const;

  bool isInteger (const PP_String& str, int&    val) const;
  bool isUnsigned(const PP_String& str, unsigned long& val) const;
  int  isDouble  (const PP_String& str, double& val) const;
  bool isBool    (const PP_String& str, bool&   val) const;

#endif /* DOXYGEN ENDIF */

};

#ifndef DOXYGEN
class PP_entry
{
private:
  friend class ParmParse;
  PP_entry()
  {}

  PP_entry (PP_String&          name,
            ParmParse::PPType typ,
            PP_List<PP_String>&    vals);

  ~PP_entry()
  {}

  PP_String               defname;
  ParmParse::PPType     deftype;
  PP_Array<PP_String>        val;

  void dump (std::ostream& os) const;
};
#endif //doxygen

//
// Inlines.
//

inline
int
ParmParse::countval (const char* name,
                     int         n) const
{
  //
  // First find n'th occurance of name in table.
  //
  const PP_entry* def = ppindex(n,name);
  return def == 0 ? 0 : def->val.length();
}

inline
void
ParmParse::get (const char* name,
                int&        ptr,
                int ival) const
{
  getval(name,ppInt,&ptr,ival,-1);
}

inline
int
ParmParse::query (const char* name,
                  int&        ptr,
                  int         ival) const
{
  return queryval(name,ppInt,&ptr,ival,-1);
}

inline
void
ParmParse::get (const char* name,
                unsigned long& ptr,
                int ival) const
{
  getval(name,ppUnsignedLong,&ptr,ival,-1);
}

inline
int
ParmParse::query (const char* name,
                  unsigned long& ptr,
                  int         ival) const
{
  return queryval(name,ppUnsignedLong,&ptr,ival,-1);
}

inline
void
ParmParse::get (const char* name,
                float&      ptr,
                int         ival) const
{
  getval(name,ppFloat,&ptr,ival,-1);
}

inline
int
ParmParse::query (const char* name,
                  float&      ptr,
                  int         ival) const
{
  return queryval(name,ppFloat,&ptr,ival,-1);
}

inline
void
ParmParse::get (const char* name,
                double&     ptr,
                int         ival) const
{
  getval(name,ppDouble,&ptr,ival,-1);
}

inline
int
ParmParse::query (const char* name,
                  double&     ptr,
                  int         ival) const
{
  return queryval(name,ppDouble,&ptr,ival,-1);
}

inline
void
ParmParse::get (const char*  name,
                std::string& ptr,
                int          ival) const
{
  PP_String pp_string;
  getval(name,ppString,&pp_string,ival,-1);
  ptr = pp_string.c_str();
}

inline
int
ParmParse::query (const char*  name,
                  std::string& ptr,
                  int          ival) const
{
  PP_String pp_string;
  int status = queryval(name,ppString,&pp_string,ival,-1);
  if (status != 0)
    ptr = pp_string.c_str();
  return status;
}

inline
void
ParmParse::get(const char* name,
                bool&       ptr,
                int         ival) const
{
  getval(name,ppBool,&ptr,ival,-1);
}

inline
int
ParmParse::query (const char* name,
                  bool&       ptr,
                  int         ival) const
{
  return queryval(name,ppBool,&ptr,ival,-1);
}

inline
void
ParmParse::getarr (const char* name,
                   std::vector<int>& ptr,
                   int         start_ix,
                   int         num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  int* c_array = new int[num_val];
  getarr(name,ppInt,c_array,start_ix,num_val,-1);
  for (int i = 0; i < num_val; ++i)
    {
      ptr[i] = c_array[i];
    }
  delete[] c_array;
}

inline
void
ParmParse::getarr (const char* name,
                   Vector<int>& ptr,
                   int         start_ix,
                   int         num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  int* c_array = new int[num_val];
  getarr(name,ppInt,c_array,start_ix,num_val,-1);
  for (int i = 0; i < num_val; ++i)
    {
      ptr[i] = c_array[i];
    }
  delete[] c_array;
}

inline
int
ParmParse::queryarr (const char* name,
                     std::vector<int>& ptr,
                     int         start_ix,
                     int         num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  int* c_array = new int[num_val];
  int status = queryarr(name,ppInt,c_array,start_ix,num_val,-1);
  if (status != 0 )
    {
      for (int i = 0; i < num_val; ++i)
        {
          ptr[i] = c_array[i];
        }
    }
  delete [] c_array;
  return status;
}

inline
int
ParmParse::queryarr (const char* name,
                     Vector<int>& ptr,
                     int         start_ix,
                     int         num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  int* c_array = new int[num_val];
  int status = queryarr(name,ppInt,c_array,start_ix,num_val,-1);
  if (status != 0 )
    {
      for (int i = 0; i < num_val; ++i)
        {
          ptr[i] = c_array[i];
        }
    }
  return status;
}

inline
void
ParmParse::getarr (const char*   name,
                   std::vector<float>& ptr,
                   int           start_ix,
                   int           num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  float* c_array = new float[num_val];
  getarr(name,ppFloat,c_array,start_ix,num_val,-1);
  for (int i = 0; i < num_val; ++i)
    {
      ptr[i] = c_array[i];
    }
  delete[] c_array;
}

inline
void
ParmParse::getarr (const char*   name,
                   Vector<float>& ptr,
                   int           start_ix,
                   int           num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  float* c_array = new float[num_val];
  getarr(name,ppFloat,c_array,start_ix,num_val,-1);
  for (int i = 0; i < num_val; ++i)
    {
      ptr[i] = c_array[i];
    }
  delete[] c_array;
}

inline
int
ParmParse::queryarr (const char*   name,
                     std::vector<float>& ptr,
                     int           start_ix,
                     int           num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  float* c_array = new float[num_val];
  int status = queryarr(name,ppFloat,c_array,start_ix,num_val,-1);
  if (status != 0)
    {
      for (int i = 0; i < num_val; ++i)
        {
          ptr[i] = c_array[i];
        }
    }
  delete[] c_array;
  return status;
}

inline
int
ParmParse::queryarr (const char*   name,
                     Vector<float>& ptr,
                     int           start_ix,
                     int           num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  float* c_array = new float[num_val];
  int status = queryarr(name,ppFloat,c_array,start_ix,num_val,-1);
  if (status != 0)
    {
      for (int i = 0; i < num_val; ++i)
        {
          ptr[i] = c_array[i];
        }
    }
  delete[] c_array;
  return status;
}

inline
void
ParmParse::getarr (const char*    name,
                   std::vector<double>& ptr,
                   int            start_ix,
                   int            num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  double* c_array = new double[num_val];
  int status = queryarr(name,ppDouble,c_array,start_ix,num_val,-1);
  if (status == 0)
    {
      cerr << "ParmParse::getarr(): "
           << name
           << " not found in table" << endl;
      dumpTable(cerr);
      MayDay::Abort();
    }
  for (int i = 0; i < num_val; ++i)
    {
      ptr[i] = c_array[i];
    }
  delete[] c_array;
}

inline
void
ParmParse::getarr (const char*    name,
                   Vector<double>& ptr,
                   int            start_ix,
                   int            num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  double* c_array = new double[num_val];
  int status = queryarr(name,ppDouble,c_array,start_ix,num_val,-1);
  if (status == 0)
    {
      cerr << "ParmParse::getarr(): "
           << name
           << " not found in table" << endl;
      dumpTable(cerr);
      MayDay::Abort();
    }
  for (int i = 0; i < num_val; ++i)
    {
      ptr[i] = c_array[i];
    }
  delete[] c_array;
}

inline
int
ParmParse::queryarr (const char*    name,
                     std::vector<double>& ptr,
                     int            start_ix,
                     int            num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  double* c_array = new double[num_val];
  int status = queryarr(name,ppDouble,c_array,start_ix,num_val,-1);
  if (status != 0)
    {
      for (int i = 0; i < num_val; ++i)
        {
          ptr[i] = c_array[i];
        }
    }
  delete[] c_array;
  return status;
}

inline
int
ParmParse::queryarr (const char*    name,
                     Vector<double>& ptr,
                     int            start_ix,
                     int            num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  double* c_array = new double[num_val];
  int status = queryarr(name,ppDouble,c_array,start_ix,num_val,-1);
  if (status != 0)
    {
      for (int i = 0; i < num_val; ++i)
        {
          ptr[i] = c_array[i];
        }
    }
  delete[] c_array;
  return status;
}

inline
void
ParmParse::getarr (const char*     name,
                   std::vector<std::string>& ptr,
                   int             start_ix,
                   int             num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  PP_String* c_array = new PP_String[num_val];
  getarr(name,ppString,c_array,start_ix,num_val,-1);
  for (int i = 0; i < num_val; ++i)
    {
      ptr[i] = c_array[i].c_str();
    }
  delete[] c_array;
}

inline
void
ParmParse::getarr (const char*     name,
                   Vector<std::string>& ptr,
                   int             start_ix,
                   int             num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  PP_String* c_array = new PP_String[num_val];
  getarr(name,ppString,c_array,start_ix,num_val,-1);
  for (int i = 0; i < num_val; ++i)
    {
      ptr[i] = c_array[i].c_str();
    }
  delete[] c_array;
}

inline
int
ParmParse::queryarr (const char*     name,
                     std::vector<std::string>& ptr,
                     int             start_ix,
                     int             num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  PP_String* c_array = new PP_String[num_val];
  int status = queryarr(name,ppString,c_array,start_ix,num_val,-1);
  if (status != 0)
    {
      for (int i = 0; i < num_val; ++i)
        {
          ptr[i] = c_array[i].c_str();
        }
    }
  delete[] c_array;
  return status;
}

inline
int
ParmParse::queryarr (const char*     name,
                     Vector<std::string>& ptr,
                     int             start_ix,
                     int             num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  PP_String* c_array = new PP_String[num_val];
  int status = queryarr(name,ppString,c_array,start_ix,num_val,-1);
  if (status != 0)
    {
      for (int i = 0; i < num_val; ++i)
        {
          ptr[i] = c_array[i].c_str();
        }
    }
  delete[] c_array;
  return status;
}

// ----------------------------------------
inline
void
ParmParse::getarr (const char* name,
                   std::vector<bool>& ptr,
                   int         start_ix,
                   int         num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  bool* c_array = new bool[num_val];
  getarr(name,ppBool,c_array,start_ix,num_val,-1);
  for (int i = 0; i < num_val; ++i)
    {
      ptr[i] = c_array[i];
    }
  delete[] c_array;
}

inline
int
ParmParse::queryarr (const char* name,
                     std::vector<bool>& ptr,
                     int         start_ix,
                     int         num_val) const
{
  if (ptr.size() < num_val)
    ptr.resize(num_val);
  bool* c_array = new bool[num_val];
  int status = queryarr(name,ppBool,c_array,start_ix,num_val,-1);
  if (status != 0 )
    {
      for (int i = 0; i < num_val; ++i)
        {
          ptr[i] = c_array[i];
        }
    }
  delete [] c_array;
  return status;
}

// ---------------------------------------

inline
bool
ParmParse::isInteger (const PP_String& str,
                      int&           val) const
{
  //
  // Start token scan.
  //
  char* endp = 0;
  val = (int) strtol(str.c_str(), &endp, 10);
  return *endp == 0;
}

inline
bool
ParmParse::isUnsigned (const PP_String& str,
                      unsigned long &   val) const
{
  char* endp = 0;
  val = strtoul(str.c_str(), &endp, 10);
  return *endp == 0;
}

inline
int
ParmParse::isDouble (const PP_String& str,
                     double&        val) const
{
  char* endp = 0;
  val = std::strtod(str.c_str(), &endp);
  return *endp == 0;
}

inline
bool
ParmParse::isBool (const PP_String& str,
                   bool&            val) const
{
  PP_String str_lc( str );
  str_lc.toLower();
  if ( str_lc == "true" || str_lc == "1" )
  {
    return val = true ;
  }
  else if ( str_lc == "false" || str_lc == "0" )
  {
    val = false ;
    return true ;
  }
  else
  {
    return false ;
  }
}

inline
int
ParmParse::countname (const std::string& name) const
{
  return countname(name.c_str());
}

#include "BaseNamespaceFooter.H"
#endif /*CH_PARMPARSE_H*/
