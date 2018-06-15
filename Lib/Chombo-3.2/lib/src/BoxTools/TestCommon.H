#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TESTCOMMON_H_
#define _TESTCOMMON_H_

#include "IntVect.H"
#include "IntVectSet.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

///sets up circle worth of tags
void setCircleTags(IntVectSet& ivs, int circleR, int thickness, const IntVect& center = IntVect::Zero);

///
extern void buildDisjointBoxLayout(DisjointBoxLayout& plan,const IntVectSet& tags, const Box& domain);

#include "NamespaceFooter.H"

#endif
