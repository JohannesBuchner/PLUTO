#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBFABVIEW_H_
#define _EBFABVIEW_H_

#include "EBAMRIO.H"
#include "NamespaceHeader.H"

int InitWriteEBFAB()
{
#ifdef CH_USE_HDF5
  writeEBFAB(NULL);

  viewEBFAB(NULL);

  writeEBFABname(NULL, "label");

  viewEBFace(NULL);

  writeEBLevel(NULL);

  viewEBLevel(NULL);

  writeEBLevelname(NULL, "label");

  writeEBAMR(NULL);

  viewEBAMR(NULL);

  writeEBAMRname(NULL, "label");

  writeIVFAB(NULL);

  viewIVFAB(NULL);

  browseIVFAB(NULL);

  writeIVFABname(NULL,"label");

  writeIVLevel(NULL);

  viewIVLevel(NULL);

  browseIVLevel(NULL);

  writeIVLevelname(NULL, "label");

  setWhichCellIndex(0);

  getWhichCellIndex();

  setCoveredCellValue(0.0);

  getCoveredCellValue();

  multiFaceValues(NULL,0,0,0,0);

  multiCellValues(NULL,0,0,0);


#endif
  return 1;
}

#include "NamespaceFooter.H"
#endif
