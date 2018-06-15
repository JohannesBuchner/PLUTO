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

#include "BaseIVFAB.H"
#include "BaseIVFABI.H"
#include "Stencils.H"
#include "EBData.H"
#include "NamespaceHeader.H"

//
// Explicit template instantiations
//
template class BaseIVFAB<Real>;
template class BaseIVFAB<int>;
template class BaseIVFAB<bool>;
template class BaseIVFAB<VoFStencil>;
template class BaseIVFAB<VolData>;

#include "NamespaceFooter.H"
#endif
