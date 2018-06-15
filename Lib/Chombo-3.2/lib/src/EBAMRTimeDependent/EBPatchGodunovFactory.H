#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPATCHGODUNOVFACTORY_H_
#define _EBPATCHGODUNOVFACTORY_H_

#include "EBPatchGodunov.H"
#include "NamespaceHeader.H"

///
/**
 */
class EBPatchGodunovFactory
{
public:
  ///
  /**
   */
  EBPatchGodunovFactory();

  ///
  /**
   */
  virtual ~EBPatchGodunovFactory();

  ///
  /**
   */
  virtual EBPatchGodunov* create() const = 0;

};
#include "NamespaceFooter.H"
#endif
