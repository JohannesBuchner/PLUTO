#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CALLBACKS_H_
#define _CALLBACKS_H_

#include <vector>
#include "BaseNamespaceHeader.H"

class Callbacks
{
public:
  typedef void (*PVF)();
  void add( PVF f );
  void run() const;

private:
  std::vector<PVF> m_funcPointers;
};
#include "BaseNamespaceFooter.H"

#endif // include guard
