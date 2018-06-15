#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELFACTORY_H_
#define _AMRLEVELFACTORY_H_

#include "AMRLevel.H"
#include "NamespaceHeader.H"

/// Factory class to produce an AMRLevel.
/**
   This is needed for construction of AMRLevel-derived classes by
   the AMR class.
*/
class AMRLevelFactory
{
public:
  ///
  /**
      Factory function.  Creates a new AMRLevel.
      In your instantiation you should include defining
      all the stuff that is not part of the AMRLevel API.
  */
  virtual AMRLevel* new_amrlevel() const = 0;

  // to make gcc v4 happy
  virtual ~AMRLevelFactory()
  {
  }
};

#include "NamespaceFooter.H"
#endif
