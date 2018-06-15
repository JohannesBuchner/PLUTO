#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBSpaceTimeFunction.H"

#include "NamespaceHeader.H"

//---------------------------------------------------------------------------
EBSpaceTimeFunction::
EBSpaceTimeFunction():
  m_domain()
{
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
EBSpaceTimeFunction::
EBSpaceTimeFunction(const ProblemDomain& a_domain):
  m_domain(a_domain)
{
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
EBSpaceTimeFunction::
~EBSpaceTimeFunction()
{
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void
EBSpaceTimeFunction::
setDomain(const ProblemDomain& a_domain)
{
  m_domain = a_domain;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void 
EBSpaceTimeFunction::
evaluate(LevelData<EBCellFAB>& a_data,
         const EBISLayout& a_indexSpace,
         RealVect a_dx,
         Real a_t) const
{
  // Iterate over all grids in this level
  for(DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
  {
    // Get at the data for this grid.
    const EBISBox& ebisBox = a_indexSpace[dit()];
    EBCellFAB& fab = a_data[dit()];

    // Evaluate!
    evaluate(fab, ebisBox, a_dx, a_t);
  }
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void 
EBSpaceTimeFunction::
evaluate(LevelData<EBCellFAB>& a_data,
         const EBISLayout& a_indexSpace,
         RealVect a_dx) const
{
  // Just call the other evaluate() with time 0.
  evaluate(a_data, a_indexSpace, a_dx, 0.0);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void 
EBSpaceTimeFunction::
evaluate(EBCellFAB& a_FAB,
         const EBISBox& a_indexSpaceBox,
         RealVect a_dx,
         Real a_t) const
{
  // Set the values of all single-valued cells in the FAB.
  BaseFab<Real>& svFab = a_FAB.getSingleValuedFAB();
  Box uBox = svFab.box();
  uBox &= m_domain;
  evaluateOnSingleValuedCells(svFab, uBox, a_dx, a_t);

  // Now do the multi-valued cells.
  IntVectSet ivs = a_indexSpaceBox.getMultiCells(uBox);
  VoFIterator vofit(ivs, a_indexSpaceBox.getEBGraph());
  evaluateOnMultiValuedCells(a_FAB, vofit, a_dx, a_t);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void 
EBSpaceTimeFunction::
evaluate(EBCellFAB& a_FAB,
         const EBISBox& a_indexSpaceBox,
         RealVect a_dx) const
{
  // Just call the other evaluate() with time 0.
  evaluate(a_FAB, a_indexSpaceBox, a_dx, 0.0);
}
//---------------------------------------------------------------------------

#include "NamespaceFooter.H"
