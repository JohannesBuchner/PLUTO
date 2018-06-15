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
  #ifdef CH_XDIR
      #undef DIM_NAMESPACE
      #undef CH_XDIR
  #endif
  #ifdef CH_MULTIDIM
    #if   CH_SPACEDIM == 1
      #define CH_XDIR Chombo::D1
    #elif CH_SPACEDIM == 2
      #define CH_XDIR Chombo::D2
    #elif CH_SPACEDIM == 3
      #define CH_XDIR Chombo::D3
    #elif CH_SPACEDIM == 4
      #define CH_XDIR Chombo::D4
    #elif CH_SPACEDIM == 5
      #define CH_XDIR Chombo::D5
    #elif CH_SPACEDIM == 6
      #define CH_XDIR Chombo::D6
    #else
      Fail now: CH_SPACEDIM must be 1, 2, 3, 4, 5 or 6.
    #endif
  #else
    #ifdef CH_NAMESPACE
      #define CH_XDIR Chombo
    #else
      #define CH_XDIR
    #endif
  #endif
#endif
