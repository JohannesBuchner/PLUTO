#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/******************************************************************************/
/**
 * \file
 *
 * \brief System dependent functions
 *
 *//*+*************************************************************************/

#include "CH_config.H"

#ifdef CHDEF_SYSTEM_HAVE_POSIXMEMALIGN
#define _XOPEN_SOURCE 600
#endif
#include <stdlib.h>

#include "CH_System.H"
#include "BaseNamespaceHeader.H"

/*============================================================================*/
//  Allocate aligned memory
/**
 *  \param[out] a_memptr
 *                      Pointer to allocated memory
 *  \param[in]  a_alignment
 *                      Alignment in bytes.  Must be a multiple of sizeof(void*)
 *                      and a power of 2.
 *  \param[in]  a_size  Number of bytes to allocate
 *  \return             0       - Success
 *                      <posix_memalign>
 *                      EINVAL  - Invalid alignment parameter
 *                      ENOMEM  - Out of memory
 *                      <malloc>
 *                      1       - Out of memory
 *  \note
 *  <ul>
 *    <li> This function returns raw memory.  Use placement new for construction
 *         of objects if required.
 *    <li> Memory allocated with memalign should be deallocated with free()
 *  </ul>
 *//*=========================================================================*/

int CHSystem::memalign(void **a_memptr, size_t a_alignment, size_t a_size)
{
#ifdef CHDEF_SYSTEM_HAVE_POSIXMEMALIGN
    return posix_memalign(a_memptr, a_alignment, a_size);
#else
    *a_memptr = malloc(a_size);
    return (*a_memptr == 0);
#endif
}

#include "BaseNamespaceFooter.H"
