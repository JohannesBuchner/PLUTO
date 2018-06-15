#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CH_SYSTEM_H_
#define _CH_SYSTEM_H_

/******************************************************************************/
/**
 * \file
 *
 * \brief System dependent functions
 *
 *//*+*************************************************************************/

#include <cstring>

#include "CH_config.H"
#include "BaseNamespaceHeader.H"

namespace CHSystem
{
  /// Allocate aligned memory
  int memalign(void **a_memptr, size_t a_alignment, size_t a_size);

/*--------------------------------------------------------------------*/
//  Find the number of bits set to 1 in a 32 bit int
/**
 *  \param[in]  a_i     Integer
 *  \return             Number of bits set to 1
 *//*-----------------------------------------------------------------*/

  inline int popcnt(unsigned i)
  {
#if CHDEF_SYSTEM_POPCNT_METHOD == 1
    int ret;
    __asm__ ("popcnt %1, %0" : "=r" (ret) : "r" (i));
    return ret;
#else
    // Known as a parallel SWAR algorithm
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return (((i + (i >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
#endif
  }
}

#include "BaseNamespaceFooter.H"
#endif
