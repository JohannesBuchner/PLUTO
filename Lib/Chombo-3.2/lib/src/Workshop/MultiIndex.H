#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTIINDEX_H_
#define _MULTIINDEX_H_

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include <map>
using std::map;

#include "Vector.H"
#include "REAL.H"
#include "IndexTM.H"

#include "Notation.H"
#include "IFData.H"

#include "NamespaceHeader.H"

// Generate and return a vector of all multi-indices whose magnitude (the sum
// of their entries) equals "a_magnitude"
template <int dim> void generateMultiIndices(Vector<IndexTM<int,dim> > & a_indices,
                                             const int                      & a_magnitude);

#include "NamespaceFooter.H"

#include "MultiIndexImplem.H"

#endif
