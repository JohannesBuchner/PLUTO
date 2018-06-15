#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _EBIBCFACTORY_H_
#define _EBIBCFACTORY_H_

#include "EBIBC.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"

#include "NamespaceHeader.H"

///
/**
 */
class EBIBCFactory
{
public:

  ///
  virtual ~EBIBCFactory();

  ///
  virtual EBIBC* create() const = 0;

  ///
  EBIBCFactory();
};

#include "NamespaceFooter.H"
#endif
