#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <stdio.h>

int main(void)
{
  unsigned long long int t0, t1;
  int result;
  unsigned int ret0[2];
  unsigned int ret1[2];
  __asm__ __volatile__("rdtsc" : "=a"(ret0[0]), "=d"(ret0[1]));
  __asm__ ("xorl %ecx, %ecx \n\t"
"L1: \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "rdtsc \n\t"
        "addl   $16, %ecx \n\t"
        "cmpl   $8192, %ecx \n\t"
        "jne    L1");
   __asm__ __volatile__("rdtsc" : "=a"(ret1[0]), "=d"(ret1[1]));

   t0 = *(unsigned long long int*)ret0;
   t1 = *(unsigned long long int*)ret1;
   result = (t1-t0)/8192;
   printf("ticks per rdtsc %d \n",result);
   return result;
}
