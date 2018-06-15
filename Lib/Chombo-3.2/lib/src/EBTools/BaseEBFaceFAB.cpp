#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifdef CH_EXPLICIT_TEMPLATES
#include "BaseEBFaceFAB.H"
#include "BoxIterator.H"
#include "BaseEBFaceFABI.H"
#include "NamespaceHeader.H"

//
// Explicit template instantiations
//
template class BaseEBFaceFAB<Real>;
template class BaseEBFaceFAB<int>;

#include "NamespaceFooter.H"
#endif // CH_EXPLICIT_TEMPLATES
