#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LGINTEGRATOR_H_
#define _LGINTEGRATOR_H_

#if (CH_SPACEDIM == 1)

#define WRHO   0
#define WVELX  1
#define WPRES  2

#define WNUM   3
#define WSLOPE 3

#define URHO  0
#define UMOMX 1
#define UENG  2

#define UNUM  3
#define UFLUX 3

#define MWAVE 0
#define ZWAVE 1
#define PWAVE 2

#endif

#if (CH_SPACEDIM == 2)

#define WRHO   0
#define WVELX  1
#define WVELY  2
#define WPRES  3

#define WNUM   4
#define WSLOPE 4

#define URHO  0
#define UMOMX 1
#define UMOMY 2
#define UENG  3

#define UNUM  4
#define UFLUX 4

#define MWAVE 0
#define ZWAVE 1
#define PWAVE 3

#endif

#if (CH_SPACEDIM == 3)

#define WRHO   0
#define WVELX  1
#define WVELY  2
#define WVELZ  3
#define WPRES  4

#define WNUM   5
#define WSLOPE 5

#define URHO  0
#define UMOMX 1
#define UMOMY 2
#define UMOMZ 3
#define UENG  4

#define UNUM  5
#define UFLUX 5

#define MWAVE 0
#define ZWAVE 1
#define PWAVE 4

#endif

#if (CH_SPACEDIM > 3)

#define WRHO   0
#define WPRES  (CH_SPACEDIM + 1)
#define WNUM   (CH_SPACEDIM + 2)

#define URHO  0
#define UENG  (CH_SPACEDIM + 1)
#define UNUM  (CH_SPACEDIM + 2)
#define UFLUX (CH_SPACEDIM + 2)

#endif

#endif
