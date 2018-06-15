#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AggStencil.H"

#include "NamespaceHeader.H"

typedef AggStencil<EBCellFAB, EBCellFAB> AggStenCellToCell;
typedef AggStencil<EBFaceFAB, EBCellFAB> AggStenFaceToCell;
typedef AggStencil<EBCellFAB, EBFaceFAB> AggStenCellToFace;

#include "NamespaceFooter.H"
