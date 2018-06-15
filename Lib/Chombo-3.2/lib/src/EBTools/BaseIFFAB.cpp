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

#include "BaseIFFAB.H"
#include "BaseIFFABI.H"
#include "Stencils.H"
#include "EBData.H"
#include "NamespaceHeader.H"

template class BaseIFFAB<FaceData>;
template class BaseIFFAB<Real>;
template class BaseIFFAB<int>;
template class BaseIFFAB<FaceStencil>;

#include "NamespaceFooter.H"
#endif
