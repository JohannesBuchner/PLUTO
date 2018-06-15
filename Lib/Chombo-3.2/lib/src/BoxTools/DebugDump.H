#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DEBUGDUMP_H_
#define _DEBUGDUMP_H_

#include "DebugOut.H"
#include "NamespaceHeader.H"

void DontYouDareCallMe()
{
  dumpMaxMin(NULL);
  dumpIVS(NULL);
  dumpIVSFAB(NULL);
  dumpBL(NULL);
  dumpLDFLoc(NULL);
  dumpLDFPar(NULL);
  dumpNodeLDFLoc(NULL);
  dumpNodeLDFPar(NULL);
  dumpNodeFAB(NULL);
  dumpFAB(NULL);
  dumpDBL(NULL);
  dumpBox(NULL);
  dumpVBox(NULL);
  dumpVVBox(NULL);
  dumpLDDBL(NULL);
  dumpBFR(NULL);
  dumpBFI(NULL);
  dumpNodeLDDBL(NULL);
}

#include "NamespaceFooter.H"
#endif
