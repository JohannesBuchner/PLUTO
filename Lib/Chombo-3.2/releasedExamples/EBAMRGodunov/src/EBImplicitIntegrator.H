#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBIMPLICITINTEGRATOR_H_
#define _EBIMPLICITINTEGRATOR_H_

#include "ImplicitIntegrator.H"
#include "EBImplicitIntegrationStrategy.H"

#include "NamespaceHeader.H"

typedef BasicImplicitIntegrator<EBImplicitIntegratorTraits> EBImplicitIntegrator;

#include "NamespaceFooter.H"
#endif
