#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBDEBUGDUMP_H_
#define _EBDEBUGDUMP_H_

#include "EBDebugOut.H"
#include "NamespaceHeader.H"

void EBDontYouDareCallMe()
{
  printMaxMinEBCF(NULL);
  printMaxMinEBFF(NULL);
  dumpVVoF(NULL);
  dumpIFFAB(NULL);
  dumpIVFAB(NULL);
  dumpLDBIVF(NULL);
  dumpFace(NULL);
  dumpVFace(NULL);
  dumpVoFSten(NULL);
  dumpFaceSten(NULL);
  dumpEBLDDBL(NULL);
  dumpLDEBCF(NULL);
  getMaxEBLevel(NULL);
  getMaxEBFAB(NULL);
  dumpEBFAB(NULL);
  dumpEBFABIrreg(NULL);
  dumpEBFABIrregGeometry(NULL);
  dumpEBLevelIrreg(NULL);
  dumpEBLevel(NULL);
  dumpEBLevelAll(NULL);
  dumpEBLevelIrregThresh(NULL, 0.0);
  dumpEBLevelThresh(NULL, 0.0);
  dumpEBAMRIrregThresh(NULL, 0.0);
  dumpEBAMRIrreg(NULL);
  dumpEBAMRThresh(NULL, 0.0);
  dumpEBAMR(NULL);
  dumpEBFABIVS(NULL, NULL);
  dumpEBFaceThresh(NULL, 0.0);
  dumpEBFace(NULL);
  dumpEBFlux(NULL);
  dumpEBFluxThresh(NULL, 0.0);
  dumpEBLevelFlux(NULL);
  dumpEBLevelFluxThresh(NULL, 0.0);
  printPointEBCF(NULL);
  printPointLDCell(NULL);
  printLocalEBCF(NULL);
  printLocalLDCell(NULL);
  printPointEBFF(NULL);
  printPointLDFlux(NULL);
  printPointEBFace(NULL);
  printMaxMinLDCell(NULL);
  printMaxMinLDFlux(NULL);
}
#include "NamespaceFooter.H"
#endif
