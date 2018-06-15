#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SetValLevel.H"

/// -----------------------------------------------------------
void setValLevels(Vector< LevelData<FArrayBox>* >&   a_levels,
                  int                                a_min,
                  int                                a_max,
                  Real                               a_val)
{
  for (int lev = a_min; lev <= a_max; lev++)
    {
      setValLevel(*a_levels[lev], a_val);
    }
}


/// -----------------------------------------------------------
void setValLevel(LevelData<FArrayBox>&   a_level,
                 Real                    a_val)
{
  for (DataIterator dit = a_level.dataIterator(); dit.ok(); ++dit)
    {
      a_level[dit].setVal(a_val);
    }
}


/// -----------------------------------------------------------
void setValLevel(LevelData<FluxBox>&   a_level,
                 Real                  a_val)
{
  for (DataIterator dit = a_level.dataIterator(); dit.ok(); ++dit)
    {
      a_level[dit].setVal(a_val);
    }
}
