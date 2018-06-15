#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRBOXESANDRANKSIO_H_
#define _AMRBOXESANDRANKSIO_H_

#include "IntVect.H"
#include "Box.H"
#include "Vector.H"
#include <string>

#include "NamespaceHeader.H"

extern void readABRfile(Vector<Vector<Box> >&      a_amrBoxes,
                        Vector<Vector<int> >&      a_amrRanks,
                        Vector<int>&               a_refinementRatio,
                        Box&                       a_baseLevelBoundingBox,
                        const std::string          a_filename,
                        const bool                 a_dataPerRank=false);


extern void writeABRfile(const Vector<Vector<Box> >&   a_amrBoxes,
                         const Vector<Vector<int> >&   a_amrRanks,
                         const Vector<int>&            a_refinementRatio,
                         const int                     a_numLevels,
                         const int                     a_targetNumProcs,
                         const Box&                    a_baseLevelBoundingBox,
                         const std::string a_filename);

#include "NamespaceFooter.H"

#endif
