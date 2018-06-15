#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBELLIPTICLOADBALANCE_H_
#define _EBELLIPTICLOADBALANCE_H_

#include "Box.H"
#include "ProblemDomain.H"
#include "Vector.H"
#include "REAL.H"
#include "EBIndexSpace.H"
#include "NamespaceHeader.H"

class DisjointBoxLayout;
///
/**
    Load balance a Vector of Boxes.
    We  make the loads bases on a few ad hoc ideas.
    we shall make fully covered boxes = constant load = covered load
    we shall say that irregular points get irregular factor more load than
    regular points.
 */
extern int
EBEllipticLoadBalance(Vector<int>&          a_proc,
                      const Vector<Box>&    a_boxes,
                      const ProblemDomain&  a_domain,
                      bool                  a_verbose = false,
                      const EBIndexSpace *a_ebisPtr = Chombo_EBIS::instance() );

///
/**
 */
extern void
getPoissonLoadsAndBoxes(Vector<unsigned long long>&  a_loads,
                        Vector<Box> &                a_boxes,
                        const DisjointBoxLayout&     a_dblOrig,
                        const ProblemDomain&         a_domain,
                        const EBIndexSpace *a_ebisPtr = Chombo_EBIS::instance() );

extern void
resetLoadOrder(Vector<unsigned long long>&  a_loads,
               Vector<Box>&                 a_newBoxes,
               Vector<Box>&                 a_oldBoxes);

#include "NamespaceFooter.H"
#endif
