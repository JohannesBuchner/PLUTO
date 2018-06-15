#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SLICESPEC_H_
#define _SLICESPEC_H_

#include "NamespaceHeader.H"

//
// Implementation note for multidim mode.  Obviously, the same SliceSpec class
// will work fine for all dimensions, so in multidim programming (e.g. in
// lib/src/MultiDim) it might seem weird to have to specify a dimensioned
// version, e.g. D1::SliceSpec.  Yet this is exactly what we have to do, the
// way multidim programming works now; we #include "SliceSpec.H.multidim", and
// that in turn declares a Chombo::D1::SliceSpec, a Chombo::D2::SliceSpec, etc
// but no plain old Chombo::SliceSpec.  One way around this is to #define
// something like "CH_DIMENSIONLESS" at the top of SliceSpec.H, add SliceSpec.H
// to the list in lib/src/BoxTools/multidim/dim-independent-headers.txt, and
// add code to NamespaceHeader.H and NamespaceFooter.H to handle it all.
// I considered that option and decided against it: on the grounds that it adds
// a complexity and, for people learning to write multidim programs, a confusing
// exception to the rule that "all regular Chombo classes are in a D<n>::
// namespace".  It's not so hard to type "D1::" or "D2::".  And in templated
// code (like that in lib/src/MultiDim/Slicing.H) we can specify the "right"
// dimension of SliceSpec by using a trait.
//

/**
  Specifies the slice we want out of an IntVect, Box, BaseFab, etc.
  For IntVect and Box slicing, only the direction matters.  For IntVect it
  refers to the element we are to remove.  For Box it refers to the dimension
  we are to remove.
*/
struct SliceSpec
{
  SliceSpec( int d, int p ) : direction(d), position(p)
  {
  }

public:
  int    direction;
  int    position;
};

#include "NamespaceFooter.H"

#endif // _SLICESPEC_H_
