#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#ifndef _BASEINDEX_H_
#define _BASEINDEX_H_

#include "REAL.H"
#include "IntVect.H"
#include "SPMD.H"
#include "NamespaceHeader.H"

/// index for other indicies to inherit
/**
   Aggstencil needed the concept of a unified
   index.   This is it.  Not much to it.
*/
class BaseIndex
{
public:

  BaseIndex()
  {
    ;
  }

  virtual ~BaseIndex()
  {
    ;
  }
};

#include "NamespaceFooter.H"
#endif
