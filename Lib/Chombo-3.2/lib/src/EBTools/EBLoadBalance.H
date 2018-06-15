#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBLOADBALANCE_H_
#define _EBLOADBALANCE_H_

#include "Box.H"
#include "ProblemDomain.H"
#include "Vector.H"
#include "REAL.H"
#include "EBIndexSpace.H"
#include "NamespaceHeader.H"

///
/**
    Load balance a Vector of Boxes.
    We  make the loads bases on a few ad hoc ideas.
    we shall make fully covered boxes = constant load = covered load
    we shall say that irregular points get irregular factor more load than
    regular points.
 */
int EBLoadBalance(Vector<int>& a_proc,
                  const Vector<Box>& a_boxes,
                  const ProblemDomain& a_domain,
                  int  a_coveredLoad = 10,
                  int  a_irregFactor = 10,
                  bool a_verbose = false,
                  const EBIndexSpace* a_ebisPtr = Chombo_EBIS::instance() );

#include "NamespaceFooter.H"
#endif
