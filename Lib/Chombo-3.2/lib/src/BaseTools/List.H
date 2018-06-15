#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LIST_H_
#define _LIST_H_

#include "MayDay.H"
#include "Pool.H"

#include "BaseNamespaceHeader.H"

template <class T> class ListLink;
template <class T> class ListIterator;
template <class T> class List;

// this is so this internal class doesn't appear in doxygen documentation
#ifndef DOXYGEN

// Internal helper class for the List class
template <class T>
class ListLink
{
private:
  friend class List<T>;
  friend class ListIterator<T>;

  ListLink (const T&     _val,
            ListLink<T>* _pre,
            ListLink<T>* _suc)
    :
    val(_val),
    pre(_pre),
    suc(_suc)
  {}

  void swap(ListLink<T>* a);

  T            val;
  ListLink<T>* pre;
  ListLink<T>* suc;
};

#endif // doxygen

/// Iterator over a List
/**
   The class ListIterator<T> is an iterator over class List<T>.
   This class does NOT provide a default constructor or an assignment operator.
*/
template <class T>
class ListIterator
{
public:

  /// Construct a ListIterator<T> to first element of aList.
  inline ListIterator (const List<T>& aList);

  /// The copy constructor.
  inline ListIterator (const ListIterator<T>& rhs);

  /// Reset this ListIterator<T> to point to the first element in the List<T>.
  inline void rewind ();

  /// Reset this ListIterator<T> to point to the first element in the List<T>.
  /**
     Same as rewind(), but included to be consistent with other iterators.
  */
  inline void begin();

  /// Return a constant reference to the object in the List<T> currently pointed to by this ListIterator<T>.
  inline const T& operator() () const;

  inline T& operator() () ;

  /// Return a constant reference to the object in the List<T> currently pointed to by this ListIterator<T>.
  inline const T& operator* () const;

  /// This is a conversion operator to makes the iterator look like a pointer
  /**
     This operator makes it easy to check if the
     iterator is pointing to an element on the List<T>.  If the
     iterator has been moved off the List<T> or if the List<T> is
     empty, this conversion returns the NULL pointer.
  */
  inline operator bool () const;

  /// Return true if the iterator is not past the end of the list.
  /**
     Same as bool(), but included to be consistent with other iterators.
  */
  inline bool ok() const ;

  /// Returns true if ListIterator<T> doesn't point to any element on the List<T>.
  inline bool operator! () const;

  /// Return a constant reference to the object in the List<T> currently pointed to by the iterator.
  inline const T& value () const;

  /// Return a constant reference to the object in the List<T> currently pointed to by the iterator.
  const T& value ();

  /// The prefix auto-increment operator.
  /**
     Advances the ListIterator<T> to point to the next element on the
     List<T>. It then returns a reference to itself to allow for chaining
     with other operators.
  */
  inline ListIterator<T>& operator++ ();

  /// The prefix auto-decrement operator.
  /**
          Moves theListIterator<T> to point to the previous element on the
          List<T>.  It then returns a reference to itself to allow for
          chaining with other operators.
  */
  inline ListIterator<T>& operator-- ();

  /// The postfix auto-decrement operator.
  /**
      Moves the ListIterator<T> to point to the previous element on the
      List<T>.  It then returns a ListIterator<T> that points to
      the old element to allow for chaining with other operators.
  */
  inline ListIterator<T> operator-- (int);

  /// The postfix auto-increment operator.
  /**
     This advances the ListIterator<T> to point to the next element
     on the List<T>. It then returns a ListIterator<T> that points to
     the old element to allow for chaining with other operators.
  */
  inline ListIterator<T> operator++ (int);

  /// Equality test for two ListIterator<T>s
  /**
     Do the two ListIterator<T>s point to the same List<T> and
     the same element within the List<T>?
  */
  inline bool operator== (const ListIterator<T>&) const;

  /// Are the ListIterator<T>s not equal?
  inline bool operator!= (const ListIterator<T>&) const;

protected:

  /**
     Construct a ListIterator<T> to a List<T> and object in that List<T>.
  */
  inline ListIterator (const List<T>& _list,
                       ListLink<T>*   _p);

  /**
     A reference to the List<T> to which we point.
  */
  const List<T>& list;

  /**
     A pointer to the element in the List<T> to which we point.
  */
  ListLink<T>* p;

private:
    friend class List<T>;
    //
    // These are disallowed.
    //
    ListIterator ();
    ListIterator<T>& operator= (const ListIterator<T>&);
};

/// A Doubly-Linked List Class
/**

  The List<T> class is a template class that implements a doubly-linked list
  of objects.  A List<T> is a useful container class when the number of
  objects in the collection is not known ahead of time.   A List<T> can
  contain an arbitrary number of elements; operations such as insertion,
  deletion, and catenation are easily implemented and inexpensive.

  The only difficulty when defining a list class is devising a mechanism to
  access the elements.  In an array, an element is accessed using an
  integer index.  Since the elements in a List<T> are ordered by position,
  we could define an integer indexing operation that walks along the
  List<T> links from the beginning until the numbered element is found.
  Unfortunately, this would be very inefficient when accessing elements near
  the end of a long list.  Another solution is to allow user access to the
  individual link objects that contain the element as well as the forward and
  backward pointers.  This is not a satisfactory solution since it allows
  user access to the internal representation of the class.  The solution
  chosen is to define a ListIterator<T> template class.

  Think of a ListIterator<T> as a pointer to an object in the List<T>.  You
  can access the element currently pointed to by the iterator, move the
  iterator forward and backward through the List<T>, and use it as a
  mechanism to define where elements should be inserted and deleted.  If the
  iterator is moved off the end of the list it behaves as a null pointer.

  This is a concrete class, not a polymorphic one.
*/

template <class T>
class List
{
public:

  /// Construct an empty List<T>.
  inline List ();

  inline List (bool usePool);

  /// The copy constructor.
  List (const List<T>& rhs);

  /// The assignment operator.
  List<T>& operator= (const List<T>& rhs);

  /// The destructor.
  inline ~List();

  /// Adds a copy of the value to the beginning of the List<T>.
  inline void prepend (const T& value);

  /// Adds a copy of the value to the end of the List<T>.
  inline void append (const T& value);

  /// Adds a copy of the value to the end of the List<T>.
  void add (const T& value);

  /// Appends a copy of all items in List<T> src to this List<T>.
  void join (const List<T>& src);

  /// Appends a copy of all items in List<T> src to this List<T>.
  /**
      This differs from join() in that it unlinks the objects from
      the List<T> src and glues them to the end of this List<T>,
      leaving List<T> src empty.  This is more efficient that
      join() if src is no longer needed.
  */
  void catenate (List<T>& src);

  /// Removes all objects from the List<T>.
  void clear ();

  /// Returns a copy of this List<T> on the heap.
  /**
      It is the user's responsibility to delete this when no longer
      needed.
  */
  inline List<T>* copy () const;

  /// Returns a reference to the first element in the List<T>.
  inline T& firstElement () const;

  /// Returns a reference to the last element in the List<T>.
  inline T& lastElement () const;

  /// Returns true if the List<T> contains an object identical to \em value.
  /**
      Type T must have an operator==() defined, or be an intrinsic type.
  */
  bool includes (const T& value) const;

  /// Returns true if the this and rhs are memberwise equal.
  /**
      Lists are memberwise equal if he two lists are the same size and
      each of the elements in the list compare equal. Type T must have
      an operator==() defined, or be an intrinsic type.
  */
  bool operator== (const List<T>& rhs) const;

  /// Returns true if the this and rhs are not equal.
  bool operator!= (const List<T>& rhs) const;

  /// Returns true if the List<T> is empty.
  inline bool isEmpty () const;

  /// Returns true if the List<T> is not empty.
  inline bool isNotEmpty () const;

  /// Returns the number of objects in the List<T>.
  int length () const;

  /// Removes the first element in the List<T>.
  inline void removeFirst ();

  /// Removes the last element in the List<T>.
  inline void removeLast ();

  /// Returns reference to object pointed to by the ListIterator<T>.
  inline const T& operator[] (const ListIterator<T>& li) const;

  /// Returns reference to object pointed to by the ListIterator<T>.
  inline T& operator[] (const ListIterator<T>& li);

  /// Removes all objects in the List<T> equal to value.
  void remove (const T& value);

  /// Removes all objects in the List<T> equal to any of the values in \rm lst.
  void remove (const List<T>& lst);

  /// Removes the object pointed to by the ListIterator<T>.
  void remove (ListIterator<T>& lit);

  /// Transfer the object pointed to by lit from the List<T> lit is associated with to this one
  void transfer(ListIterator<T>& lit);

  /// Replace the value pointed to by the ListIterator<T> by val.
  inline void replace (ListIterator<T>& li,
                       const T&         val);

  /// Insert val into List<T> after the object pointed to by \em lit.
  inline void addAfter (ListIterator<T>& lit,
                        const T&         val);

  /// Insert val into List<T> before the object pointed to by \em lit
  inline void addBefore (ListIterator<T>& lit,
                         const T&         val);

  /// Returns a ListIterator<T> to the first object in this List<T>.
  inline ListIterator<T> listIterator () const;

  /// Returns a ListIterator<T> to the first object in this List<T>.
  inline ListIterator<T> first () const;

  /// Returns a ListIterator<T> to the last object in this List<T>.
  inline ListIterator<T> last () const;

  /// sort according to operator<  (note: currently implemented with BubbleSort)
  void sort();

  void checkLinks() const;
protected:

  /**
     A helper function for removing nodes.
  */
  void remove (ListLink<T> *ln);

  void removeLink(ListLink<T> *ln);

  /**
     A helper function for adding nodes.
  */
  ListLink<T>* addBefore (ListLink<T>* ln,
                          const T&     val);

  /**
     A helper function for adding nodes.
  */
  ListLink<T>* addAfter (ListLink<T>* ln,
                         const T&     val);

  /**
     The head of the list.
  */
  ListLink<T>* head;

  /**
     The tail of the list.
  */
  ListLink<T>* tail;

  /**
     Our good and trusted friend.
  */
  friend class ListIterator<T>;

  /**
     A new member that hopefully will make our List snappier.
     In particular when you have a large number of items, like in Particle code.
  */
  static Pool linkPool;

  bool m_usePool;
};

//
// Inlines.
//

//
// The ListIterator stuff.
//

template <class T>
inline
ListIterator<T>::ListIterator (const List<T>& _list,
                               ListLink<T>*   _p)
  :
  list(_list),
  p(_p)
{
}

template <class T>
inline
ListIterator<T>::ListIterator (const List<T>& aList)
  : list(aList)
{
  p = list.head;
}

template <class T>
inline
ListIterator<T>::ListIterator (const ListIterator<T>& li)
  :
  list(li.list),
  p(li.p)
{
}

template <class T>
inline
void
ListIterator<T>::rewind ()
{
  p = list.head;
}

template <class T>
inline
void
ListIterator<T>::begin ()
{
  p = list.head;
}

template <class T>
inline
const T&
ListIterator<T>::operator() () const
{
  CH_assert(p != 0);
  return p->val;
}

template <class T>
inline
T&
ListIterator<T>::operator() ()
{
  CH_assert(p != 0);
  return p->val;
}

template <class T>
inline
const T&
ListIterator<T>::operator* () const
{
  CH_assert(p != 0);
  return p->val;
}

template <class T>
inline
bool
ListIterator<T>::ok() const
{
  return p != 0 ? true : false;
}

template <class T>
inline
ListIterator<T>::operator bool () const
{
  return ok() ;
}

template <class T>
inline
bool
ListIterator<T>::operator! () const
{
  return p == 0 ? true : false;
}

template <class T>
inline
const T&
ListIterator<T>::value () const
{
  CH_assert(p != 0);
  return p->val;
}

template <class T>
inline
ListIterator<T>&
ListIterator<T>::operator++ ()
{
  if (p)
    p = p->suc;
  return *this;
}

template <class T>
inline
ListIterator<T>&
ListIterator<T>::operator-- ()
{
  if (p)
    p = p->pre;
  return *this;
}

template <class T>
inline
ListIterator<T>
ListIterator<T>::operator++ (int)
{
  const ListIterator<T> li = *this;
  ++(*this);
  return li;
}

template <class T>
inline
ListIterator<T>
ListIterator<T>::operator-- (int)
{
  const ListIterator<T> li = *this;
  --(*this);
  return li;
}

template <class T>
inline
bool
ListIterator<T>::operator== (const ListIterator<T>& _li) const
{
  return (&list == &_li.list && p == _li.p) ? true : false;
}

template <class T>
inline
bool
ListIterator<T>::operator!= (const ListIterator<T>& _li) const
{
  return ! ListIterator<T>::operator==(_li);
}

//
// List stuff.
//

template <class T>
inline
List<T>::List ()
  :
  head(0),
  tail(0),
  m_usePool(true)
{
}

template <class T>
inline
List<T>::~List ()
{
  clear();
}

template <class T>
inline
void
List<T>::prepend (const T& value)
{
  addBefore(head, value);
}

template <class T>
inline
void
List<T>::append (const T& value)
{
  addAfter(tail, value);
}

template <class T>
inline
List<T>*
List<T>::copy () const
{
  List<T>* newlist = new List<T>(*this);
  if (newlist == 0)
    MayDay::Error("Out of memory in List::copy ");
  return newlist;
}

template <class T>
inline
T&
List<T>::firstElement () const
{
  CH_assert(head != 0);
  return head->val;
}

template <class T>
inline
T&
List<T>::lastElement () const
{
  CH_assert(tail != 0);
  return tail->val;
}

template <class T>
inline
bool
List<T>::isEmpty () const
{
  return head == 0 && tail == 0;
}

template <class T>
inline
bool
List<T>::isNotEmpty () const
{
  return !isEmpty();
}

template <class T>
inline
void
List<T>::removeFirst ()
{
  remove(head);
}

template <class T>
inline
void
List<T>::removeLast ()
{
  remove(tail);
}

template <class T>
inline
const T&
List<T>::operator[] (const ListIterator<T>& li) const
{
  CH_assert(li.p != 0);
  return li.p->val;
}

template <class T>
inline
T&
List<T>::operator[] (const ListIterator<T>& li)
{
  CH_assert(li.p != 0);
  return li.p->val;
}

template <class T>
inline
void
List<T>::replace (ListIterator<T>& li,
                  const T&         _val)
{
  CH_assert(li.p != 0);
  li.p->val = _val;
}

template <class T>
inline
void
List<T>::addAfter (ListIterator<T>& lit,
                   const T&         val)
{
  addAfter(lit.p, val);
}

template <class T>
inline
void
List<T>::addBefore (ListIterator<T>& lit,
                    const T&         val)
{
  addBefore(lit.p, val);
}

template <class T>
inline
ListIterator<T>
List<T>::first () const
{
  return ListIterator<T>(*this,head);
}

template <class T>
inline
ListIterator<T>
List<T>::listIterator () const
{
  return ListIterator<T>(*this,head);  //[NOTE: same as first()]
}

template <class T>
inline
ListIterator<T>
List<T>::last () const
{
  return ListIterator<T>(*this,tail);
}

#include "BaseNamespaceFooter.H"

#include "ListImplem.H"

#endif /*_LIST_H_*/
