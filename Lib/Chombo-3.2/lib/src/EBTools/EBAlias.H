#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL, BVS

#ifndef _EBALIAS_H_
#define _EBALIAS_H_

#include "EBCellFAB.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "NamespaceHeader.H"

class EBAliasDataFactory : public DataFactory<FArrayBox>
{
public:
  virtual ~EBAliasDataFactory()
  {
  }

  void define(LevelData<EBCellFAB>&  aliases);

  /// factory function.  creates a new 'T' object using an aliased dataPtr for T
  /** creates a new 'T' object and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user. */
  virtual FArrayBox* create(const Box& box, int ncomps, const DataIndex& a_datInd) const;

protected:
  LayoutData<Real*> aliasPtrs;
};

void aliasEB(LevelData<FArrayBox>& a_output, LevelData<EBCellFAB>& a_input);

#include "NamespaceFooter.H"
#endif
