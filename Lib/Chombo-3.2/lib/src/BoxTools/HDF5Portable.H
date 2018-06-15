#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _HDF5PORTABLE_H_
#define _HDF5PORTABLE_H_

#include "Vector.H"
#ifdef CH_USE_HDF5
#include <hdf5.h>
#else
  typedef long hid_t;
#endif

#include "Box.H"
#include "Interval.H"
#include "NamespaceHeader.H"

template <class T> void dataTypes(Vector<hid_t>& a_types,
                                  const T&       a_dummySpecializationArg);

template <class T> void dataSize(const T&        a_item,
                                 Vector<int>&    a_sizes,
                                 const Box&      a_box,
                                 const Interval& a_comps);

template <class T> void write(const T&        a_item,
                              Vector<void*>&  a_allocatedBuffers,
                              const Box&      a_box,
                              const Interval& a_comps);

template <class T> void read(T&              a_item,
                             Vector<void*>&  a_allocatedBuffers,
                             const Box&      a_box,
                             const Interval& a_comps);

template <class T> const char* name(const T& a_dummySpecializationArg);

// These are the default template declarations for the HDF5 I/O
// template functions needed by T objects for binary IO.  The default
// implementations use the linearOut/linearIn/size functions of the
// T template classes to provide non-portable binary IO.
//
// A user can add their own template specializations for
// portable binary IO abilities.
//
// To use this header, a user includes this file first, then
// declares their template specializations, then includes the CH_HDF5.H
// header file.
//
//   http://hdf.ncsa.uiuc.edu/HDF5/
//
// available datatypes in HDF5
//
//            Example
//                                                  Corresponding C Type
//            H5T_NATIVE_CHAR
//                                                  char
//            H5T_NATIVE_UCHAR
//                                                  unsigned char
//            H5T_NATIVE_SHORT
//                                                  short
//            H5T_NATIVE_USHORT
//                                                  unsigned short
//            H5T_NATIVE_INT
//                                                  int
//            H5T_NATIVE_UINT
//                                                  unsigned
//            H5T_NATIVE_LONG
//                                                  long
//            H5T_NATIVE_ULONG
//                                                  unsigned long
//            H5T_NATIVE_FLOAT
//                                                  float
//            H5T_NATIVE_DOUBLE
//                                                  double
//
// and some special ones defined for Chombo:
//
//            H5T_NATIVE_REAL                       Real
//
//            box_id                                Box
//
//            intvect_id                            IntVect

#include "NamespaceFooter.H"
#endif
