#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPHYSIBCFACTORY_H_
#define _EBPHYSIBCFACTORY_H_

#include "EBPhysIBC.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

///
/**
 */
class EBPhysIBCFactory
{
public:

  ///
  virtual ~EBPhysIBCFactory();

  ///
  virtual EBPhysIBC* create() const = 0;

  ///
  EBPhysIBCFactory();
};

#include "NamespaceFooter.H"
#endif
