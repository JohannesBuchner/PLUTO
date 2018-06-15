#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBELLIPTICLOHICENTER_H_
#define _EBELLIPTICLOHICENTER_H_

#include "Box.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

///
/**
   This function is used when in direction a_dir a 3 point stencil of cell-
   centered data is being used to compute something at the cell center of the
   center cell of the stencil.  The data for the stencil is valid in
   a_inputBox.  It divides a_inputBox into a cell-centered box where the full
   stencil can be used, a_centerBox, and cell-centered boxes where 2 point
   stencils can be used, a_loBox and a_hiBox, based on the current problem
   domain, a_domain, and the stencil direction, a_dir.  The union of these 2
   and 3 point stencel boxes is returned as a_entireBox.  The 2 point stencil
   boxes are one wide, at most, and if they have been defined then the
   corresponding flag, a_hasLo or a_hasHi, is set to one, otherwise these
   flags are zero.  All output boxes lie within a_resultBox and the domain.
 */
extern void ebEllipticLoHiCenter(Box&                 a_loBox,
                                 int&                 a_hasLo,
                                 Box&                 a_hiBox,
                                 int&                 a_hasHi,
                                 Box&                 a_centerBox,
                                 Box&                 a_entireBox,
                                 const Box&           a_resultBox,
                                 const Box&           a_inputBox,
                                 const ProblemDomain& a_domain,
                                 const int&           a_dir);

///
/**
   This function is used when in direction a_dir a 2 point stencil of cell-
   centered data is being used to compute something at the cell face between
   the cell centers of the stencil.  The data for the stencil is valid in
   a_inputBox.  It uses a_inputBox to compute a box (face-centered in a_dir)
   where the full stencil can be used, a_centerBox, and boxes (face-centered
   in a_dir) where a 1 point stencil can be used, a_loBox and a_hiBox based
   on the current problem domain, a_domain, and the stencil direction, a_dir.
   The union of these 1 and 2 point stencel boxes is returned as a_entireBox
   (face-centered in a_dir).  The 1 point stencil boxes are one wide, at most,
   and if they have been defined then the corresponding flag, a_hasLo or
   a_hasHi, is set to one, otherwise these flags are zero.  All output boxes
   lie within a_resultBox and the domain.
 */
extern void ebEllipticLoHiCenterFace(Box&                 a_loBox,
                                     int&                 a_hasLo,
                                     Box&                 a_hiBox,
                                     int&                 a_hasHi,
                                     Box&                 a_centerBox,
                                     Box&                 a_entireBox,
                                     const Box&           a_resultBox,
                                     const Box&           a_inputBox,
                                     const ProblemDomain& a_domain,
                                     const int&           a_dir);

extern void ebEllipticLoHiCenter(Vector<Box>&         a_loBox,
                                 Vector<Box>&         a_hiBox,
                                 Box&                 a_centerBox,
                                 Vector<Box>&         a_overLap,
                                 const Box&           a_dblBox,
                                 const ProblemDomain& a_domain);

#include "NamespaceFooter.H"
#endif
