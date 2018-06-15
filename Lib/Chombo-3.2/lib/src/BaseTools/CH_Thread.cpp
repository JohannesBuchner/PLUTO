#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CH_Thread.H"
#include "NamespaceHeader.H"

bool onThread0()
{
  bool retval = true;
#ifdef _OPENMP
  int thread_num = omp_get_thread_num();
  retval = (thread_num== 0);
#endif       
  return retval;
}



#include "NamespaceFooter.H"
