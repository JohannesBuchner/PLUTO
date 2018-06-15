#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef CH_LANG_FORT
#  ifdef CH_XD
#      undef DIM_NAMESPACE
#      undef CH_XD
#  endif
#  ifdef CH_MULTIDIM
#    if   CH_SPACEDIM == 1
#      define DIM_NAMESPACE D1
#      define CH_XD Chombo::D1
#    elif CH_SPACEDIM == 2
#      define DIM_NAMESPACE D2
#      define CH_XD Chombo::D2
#    elif CH_SPACEDIM == 3
#      define DIM_NAMESPACE D3
#      define CH_XD Chombo::D3
#    elif CH_SPACEDIM == 4
#      define DIM_NAMESPACE D4
#      define CH_XD Chombo::D4
#    elif CH_SPACEDIM == 5
#      define DIM_NAMESPACE D5
#      define CH_XD Chombo::D5
#    elif CH_SPACEDIM == 6
#      define DIM_NAMESPACE D6
#      define CH_XD Chombo::D6
#    else
Fail now: CH_SPACEDIM must be 1, 2, 3. 4, 5, or 6
#    endif
    namespace Chombo
    {
      namespace DIM_NAMESPACE
      {
    using namespace Chombo;

#  else
#    ifdef CH_NAMESPACE
#      define CH_XD Chombo
      namespace Chombo
      {
#    else
#      define CH_XD
#    endif
#  endif
#  ifdef INSIDE
    Fail now: matching NamespaceHeader.H not in this scope.  You probably have a \#include NamespaceHeader.H followed by another \#include NamespaceHeader.H, instead of followed by a \#include NamespaceFooter.H.
#  endif
#  define INSIDE
#endif
