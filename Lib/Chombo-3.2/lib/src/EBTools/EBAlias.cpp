#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBAlias.H"
#include "NamespaceHeader.H"

void EBAliasDataFactory::define(LevelData<EBCellFAB>& aliases)
{
  aliasPtrs.define(aliases.boxLayout());
  DataIterator dit = aliases.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      aliasPtrs[dit()] = aliases[dit()].getFArrayBox().dataPtr();
    }
}

FArrayBox* EBAliasDataFactory::create(const Box& box, int ncomps, const DataIndex& a_datInd) const
{
  //const Box& b = aliasPtrs.box(a_datInd);
  //if (b != box)
  //  {
  //    MayDay::Error("Aliased data holder dimensions do not match across LevelData const.");
  //  }
  FArrayBox* rtn = new FArrayBox(box, ncomps, aliasPtrs[a_datInd]);
  return rtn;
}

void aliasEB(LevelData<FArrayBox>& a_output, LevelData<EBCellFAB>& a_input)
{
  EBAliasDataFactory factory;
  factory.define(a_input);
  a_output.define(a_input.getBoxes(), a_input.nComp(), a_input.ghostVect(),
                  factory);
}
#include "NamespaceFooter.H"
