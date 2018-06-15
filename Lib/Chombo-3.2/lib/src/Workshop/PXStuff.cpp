#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <stdio.h>
#include "PXStuff.H"

#include "NamespaceHeader.H" // for Chombo

const char *PXE_ErrorCodeName[17] = { "No Error",
                                      "Memory",
                                      "Bad Input",
                                      "Nonphysical",
                                      "Read/Write",
                                      "Grid",
                                      "Search Not Found",
                                      "No Update",
                                      "Parallel",
                                      "Code Flow",
                                      "System",
                                      "Dynamic Library",
                                      "Not Converged",
                                      "Viz",
                                      "Lapack",
                                      "Hard Exit",
                                      "CGNS"};

void PXErrorReport( const char *file, const int line, const char *call, const int ierr){

/*   if (ierr == PX_NO_ERROR) */
/*     return PX_NO_ERROR; */

  printf("Error %d (%s) has occured.\n File : %s  Line : %d\n Call : %s\n", ierr,PXE_ErrorCodeName[-ierr\
], file, line, call);
  //  printf("Error %d has occured.\n   File : %s   Line : %d\n   Call : %s\n", ierr, file, line, call);
  fflush(stdout);
}

#include "NamespaceFooter.H" // for Chombo
