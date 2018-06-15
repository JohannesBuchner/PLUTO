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
#     undef DIM_NAMESPACE
#       undef CH_XD
#  endif
#  ifdef CH_NAMESPACE
#     define CH_XD Chombo
      namespace Chombo
      {
#  else
#define CH_XD
#  endif
#  ifdef INSIDEBASE
    Fail now: matching BaseNamespaceHeader.H not in this scope.  You probably have a \#include BaseNamespaceHeader.H followed by another \#include BaseNamespaceHeader.H, instead of followed by a \#include BaseNamespaceFooter.H.
#  endif
#  define INSIDEBASE
#endif
