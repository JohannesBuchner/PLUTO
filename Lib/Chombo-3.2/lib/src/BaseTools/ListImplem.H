#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LISTIMPLEM_H_
#define _LISTIMPLEM_H_

#include "MayDay.H"
#include "BaseNamespaceHeader.H"

//
// List members
//
template <class T>
Pool List<T>::linkPool(sizeof(ListLink<T>), typeid(T).name(), 300);

template <class T>
List<T>::List (const List<T>& source)
    : head(0),
      tail(0),
      m_usePool(source.m_usePool)
{
    if (source.isEmpty())
        tail = head = 0;
    else
        for (ListIterator<T> li(source); li; ++li)
            append(li());
}
template <class T>
List<T>::List (bool usePool)
    : head(0),
      tail(0),
      m_usePool(usePool)
{
}

//
// This isn't inlined as it's declared virtual.
//

template <class T>
void
List<T>::add (const T& value)
{
    append(value);
}

template <class T>
int
List<T>::length () const
{
    int len = 0;
    for (ListIterator<T> li(*this); li; ++li)
        len++;
    return len;
}

template <class T>
List<T>&
List<T>::operator= (const List<T>& source)
{
    if (!(this == &source))
    {
        clear();
        for (ListIterator<T> li(source); li; ++li)
            append(li());
    }
    return *this;
}

template <class T>
ListLink<T> *
List<T>::addBefore (ListLink<T>* ln,
                    const T&     val)
{
    CH_assert(ln != 0 || head == 0);

    ListLink<T>* newlink;

    if (ln == head)
    {
      //head = newlink = new ListLink<T>(val, 0, head);
      if (m_usePool)
        head = newlink = new (linkPool.getPtr()) ListLink<T>(val, 0, head);
      else
        head = newlink = new ListLink<T>(val, 0, head);

      if (tail == 0)
        tail = head;
      else
        head->suc->pre = newlink;
    }
    else
    {
      //newlink = new ListLink<T>(val, ln->pre, ln);
      if (m_usePool)
        newlink = new (linkPool.getPtr()) ListLink<T>(val, ln->pre, ln);
      else
        newlink = new ListLink<T>(val, ln->pre, ln);

      ln->pre->suc = newlink;
      ln->pre = newlink;
    }

    if (newlink == 0)
      MayDay::Error("Out of memory in ListLink::addBefore");

    return newlink;
}

//  g1 <=> this <=> g2   . . .  k1 <=> a <=> k2
// to
//   g1 <=> a <=> g2   . . .  k1 <=> this <=> k2
template <class T>
void
ListLink<T>::swap(ListLink<T>* a)
{
  CH_assert(a!= NULL);
  if (a->pre == this)// if a and this are neighbours, need to use different logic
    {
      a->pre = this->pre;
      if (this->pre != NULL) pre->suc = a;
      this->pre = a;
      this->suc = a->suc;
      if (this->suc != NULL) suc->pre = this;
      a->suc = this;
    }
  else if (a->suc == this)
    {
      a->swap(this);
    }
  else
    {
      ListLink<T>* tmp = a->pre; //tmp = k1
      a->pre = this->pre;              //a.pre = g1
      if (this->pre != NULL)
        {
          pre->suc = a;          //g1.suc = a
        }
      this->pre = tmp;                 //this.pre  = k1
      tmp = a->suc;              //tmp=k2
      a->suc = this->suc;              //a.suc = g2
      if (this->suc != NULL)
        {
          suc->pre = a;          //g2.pre = a
        }
      if (this->pre!= NULL)
        {
          pre->suc = this;       //k1.suc = this
        }
      this->suc = tmp;                 //this.suc = k2
      if (suc != NULL)
        {
          suc->pre = this;       //k2.pre = this
        }
    }
}

template <class T>
ListLink<T>*
List<T>::addAfter (ListLink<T>* ln,
                   const T&     val)
{
    CH_assert(ln != 0 || tail == 0);

    ListLink<T>* newlink;

    if (ln == tail)
    {
      if (!m_usePool)
        tail = newlink = new ListLink<T>(val,tail,0);
      else
        tail  = newlink = new (linkPool.getPtr()) ListLink<T>(val,tail,0);
      if (head == 0)
        head = tail;
      else
        tail->pre->suc = newlink;
    }
    else
    {
      if (!m_usePool)
        newlink = new ListLink<T>(val, ln, ln->suc);
      else
        newlink = new (linkPool.getPtr()) ListLink<T>(val, ln, ln->suc);
      ln->suc->pre = newlink;
      ln->suc = newlink;
    }

    if (newlink == 0)
      MayDay::Error("Out of memory in ListLink::addAfter");

    return newlink;
}

template <class T>
void
List<T>::join (const List<T>& list2)
{
    for (ListIterator<T> li2(list2); li2; ++li2)
        append(li2());
}

template <class T>
void
List<T>::catenate (List<T>& list2)
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
void List<T>::checkLinks() const
{
  ListLink<T>* next = 0;

  for (ListLink<T>* p = head; p != 0; p = next)
    {
        next = p->suc;
        if (next != NULL)
          {
            if (next->pre != p)
              {
                MayDay::Abort("next->pre != this");
              }
          }
        else
        {
          if (p != tail)
            {
              MayDay::Abort("link terminated and is not tail");
            }
        }

    }
}

template <class T>
void
List<T>::clear ()
{
    ListLink<T>* next = 0;

    for (ListLink<T>* p = head; p != 0; p = next)
    {
        next = p->suc;
        p->suc = 0;
        //delete p;
        p->val.~T();
        if (m_usePool)
          linkPool.returnPtr(p);
        else
          delete p;

    }
    tail = head = 0;
}

template <class T>
bool
List<T>::includes (const T& v) const
{
    bool rc = false;
    for (ListIterator<T> li(*this); li && !rc; ++li)
        if (v == li())
            rc = true;
    return rc;
}

template<class T>
bool
List<T>::operator== (const List<T>& rhs) const
{
    if (length() == rhs.length())
    {
        for (ListIterator<T> li(*this), ri(rhs); li; ++li, ++ri)
            if (li() != ri())
                return false;
        return true;
    }

    return false;
}

template<class T>
bool
List<T>::operator!= (const List<T>& rhs) const
{
    return !operator==(rhs);
}

template <class T>
void
List<T>::transfer(ListIterator<T>& li)
{

  CH_assert(&(li.list) != this);
  List<T>& other = (List<T>&)(li.list);

  ListLink<T>* p = li.p;
  li.p = p->suc; //push iterator ahead;

  // remove p from other list.
  other.removeLink(p);

  if (head == 0)
  {
    head = tail = p;
    p->pre=0;
    p->suc=0;
  }
  else
  {
    p->suc=0;
    tail->suc = p;
    p->pre = tail;
    tail = p;
  }

}

template <class T>
void
List<T>::remove (ListIterator<T>& li)
{
    ListLink<T> *np = li.p->suc;
    remove(li.p);
    li.p = np;
}

template <class T>
void
List<T>::remove (const T& _v)
{
    for (ListIterator<T> litr(*this); litr; ++litr)
        if (litr() == _v)
            remove(litr);
}

template <class T>
void
List<T>::remove (const List<T>& _lv)
{
    for (ListIterator<T> litr(_lv); litr; ++litr)
        remove(litr());
}

template <class T>
void
List<T>::removeLink (ListLink<T>* ln)
{
    CH_assert(head !=0 && tail != 0);

    if (head == tail)
      {
        CH_assert(head == ln);
        head = tail = 0;
      }
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

}

template <class T>
void
List<T>::remove(ListLink<T>* ln)
{
  removeLink(ln);
  //delete ln;
  ln->val.~T();
  if (m_usePool)
    linkPool.returnPtr(ln);
  else
    delete ln;

  ln = 0;
}

template <class T>
void
List<T>::sort()
{
  if (head == NULL) return;

  ListLink<T>* next = 0;
  bool swapped = true;
  while (swapped)
  {
    swapped = false;
    for (ListLink<T>* p = head; p != 0;)
      {
        next = p->suc;
        if (next != NULL && next->val < p->val)
        {
          swapped = true;
          p->swap(next);
          if (head == p) head = next;
          if (tail == next) tail = p;
        }
        else
        {
          p = next;
        }
      }
  }
  //checkLinks();
}

#include "BaseNamespaceFooter.H"
#endif  /* CH_LISTIMPLEM_H */
