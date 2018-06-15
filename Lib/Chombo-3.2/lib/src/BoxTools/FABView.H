#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FABVIEW_H_
#define _FABVIEW_H_

#include "AMRIO.H"
#include "NodeAMRIO.H"
#include "NamespaceHeader.H"

int InitWriteFAB()
{
#ifdef CH_USE_HDF5
  writeFAB(NULL);

  viewFAB(NULL);

  viewBFI(NULL);

  viewBFIV(NULL);

  viewBFRV(NULL);

  viewIVSFAB(NULL);

  viewIVSFABI(NULL);

  viewIVSFABIV(NULL);

  viewIVSFABRV(NULL);

  browseFAB(NULL);

  browseLevel(NULL);

  writeBFR(NULL);

  viewBFR(NULL);

  writeFABname(NULL, "label");

  writeLevel(NULL);

  viewLevel(NULL);

  viewLevelNoGhost(NULL);

  writeLevelname(NULL, "label");

  writeDBL(NULL);

  viewDBL(NULL);

  writeDBLname(NULL, "label");

  viewLevelBFI(NULL);

  viewIVS(NULL);

  viewVectorBox(NULL);

  writeNFAB(NULL);

  viewNFAB(NULL);

  writeNFABname(NULL, "label");

  writeNodeLevel(NULL);

  viewNodeLevel(NULL);

  writeNodeLevelname(NULL, "label");

#endif // CH_USE_HDF5

  return 1;
}

#include "NamespaceFooter.H"
#endif
