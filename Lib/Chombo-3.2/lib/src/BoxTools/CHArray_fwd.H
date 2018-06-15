#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CHARRAY_FWD_H_
#define _CHARRAY_FWD_H_

/******************************************************************************/
/**
 * \file
 *
 * \brief Forward declarations for class CHArray
 *
 * These class declarations have to be in a separate file to ensure that the
 * default template parameters are only specified once.  Include this file like
 * normal if all you want is a foward declaration of class CHArray.
 *
 *//***************************************************************************/

#ifdef CH_SPACEDIM
// In Chombo land
  #include "NamespaceHeader.H"
#endif

/// Array configuration parameters
enum
{
  ArZeroRow,                          ///< Zero lower bound and row-storage
  ArZeroCol,                          ///< Zero lower bound and column-storage
  ArRangeRow,                         ///< Range of subscripts and row-storage
  ArRangeCol                          ///< Range of subscripts and column-
                                      ///< storage
};

/// Array limits
enum
{
  CHArrayMaximumRank = 7              ///< Maximum supported rank
};

namespace ArSp
{
  template <typename T> class DefaultArrayAlloc;
}
template <typename T,
          unsigned Rank,
          int ArConf = ArZeroCol,
          typename Alloc = ArSp::DefaultArrayAlloc<T> >
class CHArray;

#ifdef CH_SPACEDIM
  #include "NamespaceFooter.H"
#endif
#endif
