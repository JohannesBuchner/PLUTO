#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PRINTABLE_H_
#define _PRINTABLE_H_

#include <ostream>
#include "BaseNamespaceHeader.H"

/**
 *  This class encapsulates the concept of `printability',
 *   it automatically associates the output operator
 *    to a derived class via a virtual method
 *    to get around the problem of polymorphic operator being illegal in C++.
 *
 *  Typical usage:
 *  (1) inherit from Printable
 *     Class A : Class Printable
 *  (2) define the method print() in A
 *  (3) use '<<' operator for objects of class A
 *
 *  Variants of print can be added for different flavors.
 */
class Printable
{
public:
  virtual ~Printable()
  {
  }

  virtual void print(std::ostream& ) const = 0;
};

inline std::ostream& operator<< (std::ostream& os, const Printable& p)
{
  p.print(os);
  return os;
}

#include "BaseNamespaceFooter.H"
#endif
