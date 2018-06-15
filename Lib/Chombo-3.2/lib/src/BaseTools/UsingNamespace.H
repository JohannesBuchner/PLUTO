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
#      define CH_XD Chombo::D1
#    elif CH_SPACEDIM == 2
#      define CH_XD Chombo::D2
#    elif CH_SPACEDIM == 3
#      define CH_XD Chombo::D3
#    elif CH_SPACEDIM == 4
#      define CH_XD Chombo::D4
#    elif CH_SPACEDIM == 5
#      define CH_XD Chombo::D5
#    elif CH_SPACEDIM == 6
#      define CH_XD Chombo::D6
#    else
      Fail now: CH_SPACEDIM must be 1, 2, 3, 4, 5 or 6.
#    endif
        // if we're in multidim, still access to Chombo::D(dim) namespace
        using namespace CH_XD;
#  else
#    ifdef CH_NAMESPACE
#      define CH_XD Chombo
#    else
#      define CH_XD
#    endif
#  endif
#  ifdef CH_NAMESPACE
        using namespace Chombo;
#  endif
#endif
