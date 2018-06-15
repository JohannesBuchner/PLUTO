#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _UGIO_H_
#define _UGIO_H_

#include <string>
using std::string;
using std::istream;

#include "LevelData.H"
#include "REAL.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "IntVectSet.H"
#include "NamespaceHeader.H"

#ifdef CH_USE_HDF5
/// Write single-level of data in HDF5 format
/**
   Writes a single level in HDF5 format.  Only available if the preprocessor
   macro HDF5 is defined at compilation.  This is blocking.

   \b Arguments:
   - \b a_filename: file to output to.
   - \b a_grid    : grids
   - \b a_data    : data
   - \b a_names   : names of variables.
   - \b a_domain  : domain of level.
*/
void WriteUGHDF5(const string&               a_filename,
                 const DisjointBoxLayout&    a_grids,
                 const LevelData<FArrayBox>& a_data,
                 const Box&                  a_domain);

/// Write a single level of data in HDF5 format.
/**
   Writes a single level in HDF5 format.  Only available if the preprocessor
   macro HDF5 is defined at compilation.  This is the hdf5 handle version and
   is not blocking.

   \b Arguments:
   - \b a_handle: handle to output to.
   - \b a_grid  : grids
   - \b a_data  : data
   - \b a_names : names of variables.
   - \b a_domain: domain of level.

*/
void WriteUGHDF5(HDF5Handle&                 a_handle,
                 const DisjointBoxLayout&    a_grids,
                 const LevelData<FArrayBox>& a_data,
                 const Box&                  a_domain);

/// Reads a level of data in HDF5 format.
/**
   Only available if the preprocessor macro HDF5 is defined at compilation.
   Returns 0 on success.  This is blocking.

   \b Arguments:
   - \b a_filename: file to input from.
   - \b a_grids   : grids
   - \b a_data    : data
   - \b a_domain  : domain

   \b Returns:
   - status code with values:
     -  0: success
     - -1: bogus number of levels
     - -2: bogus number of components
     - -3: error in readlevel
     - -4: file open failed
*/
int ReadUGHDF5(const string&         a_filename,
               DisjointBoxLayout&    a_grids,
               LevelData<FArrayBox>& a_data,
               Box&                  a_domain);

/// This is the handle version of ReadUGHDF5.
/**
   This is not blocking.
   Reads a  level of data in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   \b Arguments:
   - \b a_handle: hdf5 handle to input from.
   - \b a_grids : grids
   - \b a_data  : data
   - \b a_domain: domain

   \b Returns:
   - status code with values:
      -  0: success
      - -1: bogus number of levels
      - -2: bogus number of components
      - -3: error in readlevel
      - -4: file open failed
*/
int
ReadUGHDF5(HDF5Handle&           a_handle,
           DisjointBoxLayout&    a_grids,
           LevelData<FArrayBox>& a_data,
           Box&                  a_domain);
#endif

#include "NamespaceFooter.H"
#endif
