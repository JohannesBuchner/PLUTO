#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>

#include "unidim.H"

using namespace Chombo::D2;

void Unidimmer::Func()
{
    pout() << "Unidimmer::Func(): herewith is a box: " << this->box << '\n';
}
