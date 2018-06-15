#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FORTRANNAMEMACRO_H_
#define _FORTRANNAMEMACRO_H_

#ifndef CH_MULTIDIM
      #ifdef CH_FORT_UPPERCASE
        #ifdef CH_FORT_UNDERSCORE
           #define FORTRAN_NAME( NAME ,name ) NAME ## _
        #else
           #define FORTRAN_NAME( NAME ,name ) NAME
        #endif
      #else
        #ifdef CH_FORT_UNDERSCORE
           #define FORTRAN_NAME( NAME ,name ) name ## _
        #else
           #define FORTRAN_NAME( NAME ,name ) name
        #endif
      #endif
#else
  #ifndef LAST_MULTIDIM_FortranNameMacro_SPACEDIM
    #define LAST_MULTIDIM_FortranNameMacro_SPACEDIM 0
  #endif
  #if CH_SPACEDIM != LAST_MULTIDIM_FortranNameMacro_SPACEDIM
    #undef LAST_MULTIDIM_FortranNameMacro_SPACEDIM
    #ifdef FORTRAN_NAME
      #undef FORTRAN_NAME
    #endif
    #if   CH_SPACEDIM == 0
    #elif CH_SPACEDIM == 1
      #ifdef CH_FORT_UPPERCASE
        #ifdef CH_FORT_UNDERSCORE
           #define FORTRAN_NAME( NAME ,name ) NAME ## _  ## 1d_
        #else
           #define FORTRAN_NAME( NAME ,name ) NAME  ## 1d
        #endif
      #else
        #ifdef CH_FORT_UNDERSCORE
           #define FORTRAN_NAME( NAME ,name ) name ## _  ## 1d_
        #else
           #define FORTRAN_NAME( NAME ,name ) name  ## 1d
        #endif
      #endif
      #define LAST_MULTIDIM_FortranNameMacro_SPACEDIM 1
    #elif CH_SPACEDIM == 2
      #ifdef CH_FORT_UPPERCASE
        #ifdef CH_FORT_UNDERSCORE
           #define FORTRAN_NAME( NAME ,name ) NAME ## _  ## 2d_
        #else
           #define FORTRAN_NAME( NAME ,name ) NAME  ## 2d
        #endif
      #else
        #ifdef CH_FORT_UNDERSCORE
           #define FORTRAN_NAME( NAME ,name ) name ## _  ## 2d_
        #else
           #define FORTRAN_NAME( NAME ,name ) name  ## 2d
        #endif
      #endif
      #define LAST_MULTIDIM_FortranNameMacro_SPACEDIM 2
    #elif CH_SPACEDIM == 3
      #ifdef CH_FORT_UPPERCASE
        #ifdef CH_FORT_UNDERSCORE
           #define FORTRAN_NAME( NAME ,name ) NAME ## _  ## 3d_
        #else
           #define FORTRAN_NAME( NAME ,name ) NAME  ## 3d
        #endif
      #else
        #ifdef CH_FORT_UNDERSCORE
           #define FORTRAN_NAME( NAME ,name ) name ## _  ## 3d_
        #else
           #define FORTRAN_NAME( NAME ,name ) name  ## 3d
        #endif
      #endif
      #define LAST_MULTIDIM_FortranNameMacro_SPACEDIM 3
    #elif CH_SPACEDIM == 4
      #ifdef CH_FORT_UPPERCASE
        #ifdef CH_FORT_UNDERSCORE
           #define FORTRAN_NAME( NAME ,name ) NAME ## _  ## 4d_
        #else
           #define FORTRAN_NAME( NAME ,name ) NAME  ## 4d
        #endif
      #else
        #ifdef CH_FORT_UNDERSCORE
           #define FORTRAN_NAME( NAME ,name ) name ## _  ## 4d_
        #else
           #define FORTRAN_NAME( NAME ,name ) name  ## 4d
        #endif
      #endif
      #define LAST_MULTIDIM_FortranNameMacro_SPACEDIM 4
    #elif CH_SPACEDIM == 5
      #ifdef CH_FORT_UPPERCASE
        #ifdef CH_FORT_UNDERSCORE
           #define FORTRAN_NAME( NAME ,name ) NAME ## _  ## 5d_
        #else
           #define FORTRAN_NAME( NAME ,name ) NAME  ## 5d
        #endif
      #else
        #ifdef CH_FORT_UNDERSCORE
           #define FORTRAN_NAME( NAME ,name ) name ## _  ## 5d_
        #else
           #define FORTRAN_NAME( NAME ,name ) name  ## 5d
        #endif
      #endif
      #define LAST_MULTIDIM_FortranNameMacro_SPACEDIM 5
    #elif CH_SPACEDIM == 6
      #ifdef CH_FORT_UPPERCASE
        #ifdef CH_FORT_UNDERSCORE
           #define FORTRAN_NAME( NAME ,name ) NAME ## _  ## 6d_
        #else
           #define FORTRAN_NAME( NAME ,name ) NAME  ## 6d
        #endif
      #else
        #ifdef CH_FORT_UNDERSCORE
           #define FORTRAN_NAME( NAME ,name ) name ## _  ## 6d_
        #else
           #define FORTRAN_NAME( NAME ,name ) name  ## 6d
        #endif
      #endif
      #define LAST_MULTIDIM_FortranNameMacro_SPACEDIM 6
    #else
      Fail now: CH_SPACEDIM must be 1, 2, 3, 4, 5 or 6.
    #endif
  #endif
#endif

#endif
