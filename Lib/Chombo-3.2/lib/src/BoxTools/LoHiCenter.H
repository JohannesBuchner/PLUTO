#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LOHICENTER_H_
#define _LOHICENTER_H_

#include "Box.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

///
/**
   This function is used when in direction a_dir a 3-point stencil of
   cell-centered data is being used to compute something at the
   center of the central cell of the stencil.

   Inputs:
   * a_dir is the direction of the 3-point stencil;
   * a_inBox is the cell-centered box on which we are to compute a
     stencil (2-point or 3-point), expanded by 1 in a_dir direction,
     and not intersected with the domain;
   * a_domain is the problem domain.

   Output boxes are all cell-centered subboxes of a_inBox and
   are also contained in a_domain:
   * a_loBox is where a 2-point stencil must be used on the low side;
   * a_centerBox is where the full 3-point stencil can be used;
   * a_hiBox is where a 2-point stencil must be used on the high side;
   * a_entireBox is the union of a_loBox and a_centerBox and a_hiBox.

   The boxes a_loBox and a_hiBox will be at most 1 cell wide;
   here the 2-point stencil consists of a cell from this box and
   the neighboring cell from a_centerBox in direction a_dir.

   a_entireBox will be a_inBox with one layer of cells removed
   from both sides in dimension a_dir.

   Output flags:
   * a_hasLo:  1 or 0, according to whether a_loBox is defined or not;
   * a_hasHi:  1 or 0, according to whether a_hiBox is defined or not.
 */
extern void loHiCenter(Box&                 a_loBox,
                       int&                 a_hasLo,
                       Box&                 a_hiBox,
                       int&                 a_hasHi,
                       Box&                 a_centerBox,
                       Box&                 a_entireBox,
                       const Box&           a_inBox,
                       const ProblemDomain& a_domain,
                       const int&           a_dir);

///
/**
   This function is used when in direction a_dir a 2-point stencil of
   cell-centered data is being used to compute something at the
   face between cells of the stencil.

   Inputs:
   * a_dir is the direction of the 2-point stencil;
   * a_inBox is the cell-centered box (on whose faces we are to
     compute a stencil), expanded by 1 in a_dir direction, and not
     intersected with the domain.
   * a_domain is the problem domain.

   Output boxes are all face-centered subboxes of the faces of a_inBox
   in direction a_dir, and are also contained in the faces of a_domain:
   * a_loBox is where a 1-point stencil must be used on the low side;
   * a_centerBox is where the full 2-point stencil can be used;
   * a_hiBox is where a 1-point stencil must be used on the high side;
   * a_entireBox is the union of a_loBox and a_centerBox and a_hiBox.

   The boxes a_loBox and a_hiBox will be at most 1 face wide;
   here the 1-point stencil consists of only 1 of the adjacent cells
   in direction a_dir.

   Output flags:
   * a_hasLo:  1 or 0, according to whether a_loBox is defined or not;
   * a_hasHi:  1 or 0, according to whether a_hiBox is defined or not.
 */
extern void loHiCenterFace(Box&                 a_loBox,
                           int&                 a_hasLo,
                           Box&                 a_hiBox,
                           int&                 a_hasHi,
                           Box&                 a_centerBox,
                           Box&                 a_entireBox,
                           const Box&           a_inBox,
                           const ProblemDomain& a_domain,
                           const int&           a_dir);

///
/**
   This function is used when in direction a_dir a 4-point stencil of
   cell-centered data is being used to compute something at the
   face between cells of the stencil.

   Inputs:
   * a_dir is the direction of the 4-point stencil;
   * a_inBox is the cell-centered box on which we have data;
   * a_domain is the problem domain.

   Output boxes are all face-centered subboxes of the faces of a_inBox
   in direction a_dir, and are also contained in the faces of a_domain:
   * a_loBox is where a stencil must be used on the low side
     with all four points higher;
   * a_nextLoBox is where a stencil must be used on the low side
     with one point lower and three points higher;
   * a_hiBox is where a stencil must be used on the high side
     with all four points lower;
   * a_nextHiBox is where a stencil must be used on the high side
     with one point higher and three points lower;
   * a_centerBox is union of a_innerCenterBox, a_nextLoBox, and a_nextHiBox,
     where the regular 2-point stencil can be used;
   * a_innerCenterBox is where the regular 4-point stencil can be used;
   * a_entireBox is the union of a_(lo|nextLo|center|hi|nextHi)Box.

   Each of the boxes a_loBox, a_nextLoBox, a_hiBox, a_nextHiBox
   will be at most 1 face wide.

   a_loBox and a_nextLoBox will both be defined or both be undefined.
   a_hiBox and a_nextHiBox will both be defined or both be undefined.

   Output flags:
   * a_hasLo:  1 or 0, according to whether a_loBox and a_nextLoBox
     are defined or not;
   * a_hasHi:  1 or 0, according to whether a_hiBox and a_nextHiBox
     are defined or not.
 */
extern void loHiCenterFace4(Box&                 a_loBox,
                            Box&                 a_nextLoBox,
                            int&                 a_hasLo,
                            Box&                 a_hiBox,
                            Box&                 a_nextHiBox,
                            int&                 a_hasHi,
                            Box&                 a_centerBox,
                            Box&                 a_innerCenterBox,
                            Box&                 a_entireBox,
                            const Box&           a_inBox,
                            const ProblemDomain& a_domain,
                            const int&           a_dir);

///
/**
   This function is used when in direction a_dir a 5-point stencil of
   cell-centered data is being used to compute something at the
   center of the central cell of the stencil.

   Inputs:
   * a_dir is the direction of the 5-point stencil;
   * a_inBox is the cell-centered box on which we are to compute a
     stencil (2-point or 3-point or 5-point), expanded by 1 in a_dir direction,
     and not intersected with the domain;
   * a_domain is the problem domain.

   Output boxes are all cell-centered subboxes of a_inBox and
   are also contained in a_domain:
   * a_loBox is where a 2-point stencil must be used on the low side;
   * a_nextLoBox is where a 3-point stencil must be used on the low side
     with one point lower and three points higher;
   * a_hiBox is where a 2-point stencil must be used on the high side;
   * a_nextHiBox is where a 3-point stencil must be used on the high side
     with one point higher and three points lower;
   * a_centerBox is union of a_innerCenterBox, a_nextLoBox, and a_nextHiBox,
     where a 3-point stencil may be used;
   * a_innerCenterBox is where the regular 5-point stencil can be used;
     where a 5-point stencil may be used;
   * a_entireBox is the union of a_(lo|nextLo|center|hi|nextHi)Box.

   Each of the boxes a_loBox, a_nextLoBox, a_hiBox, a_nextHiBox
   will be at most 1 cell wide.

   a_loBox and a_nextLoBox will both be defined or both be undefined.
   a_hiBox and a_nextHiBox will both be defined or both be undefined.

   Output flags:
   * a_hasLo:  1 or 0, according to whether a_loBox and a_nextLoBox
     are defined or not;
   * a_hasHi:  1 or 0, according to whether a_hiBox and a_nextHiBox
     are defined or not.
 */
extern void loHiCenter5(Box&                 a_loBox,
                        Box&                 a_nextLoBox,
                        int&                 a_hasLo,
                        Box&                 a_hiBox,
                        Box&                 a_nextHiBox,
                        int&                 a_hasHi,
                        Box&                 a_centerBox,
                        Box&                 a_innerCenterBox,
                        Box&                 a_entireBox,
                        const Box&           a_inBox,
                        const ProblemDomain& a_domain,
                        const int&           a_dir);

///
/**
   This function is used when in direction a_dir a 6-point stencil of
   cell-centered data is being used to compute something at the
   face between cells of the stencil.

   Inputs:
   * a_dir is the direction of the 6-point stencil;
   * a_inBox is the cell-centered box on which we have data;
   * a_domain is the problem domain.

   Output boxes are all face-centered subboxes of the faces of a_inBox
   in direction a_dir, and are also contained in the faces of a_domain:
   * a_farLoBox is where a stencil must be used on the low side
     with all 6 points higher;
   * a_midLoBox is where a stencil must be used on the low side
     with 1 point lower and 5 points higher;
   * a_nearLoBox is where a stencil must be used on the low side
     with 2 points lower and 4 points higher;
   * a_farHiBox is where a stencil must be used on the high side
     with all 6 points lower;
   * a_midHiBox is where a stencil must be used on the high side
     with 1 point higher and 5 points lower;
   * a_nearHiBox is where a stencil must be used on the high side
     with 2 points higher and 4 points lower;
   * a_centerBox is union of a_innerCenterBox and a_near(Lo|Hi)Box,
     where the regular 2-point stencil can be used;
   * a_midCenterBox is union of a_innerCenterBox, a_near(Lo|Hi)Box, and
     a_mid(Lo|Hi)Box, where the regular 4-point stencil can be used;
   * a_innerCenterBox is where the regular 6-point stencil can be used;
   * a_entireBox is union of a_(far|mid|near)(Lo|Hi)Box and a_innerCenterBox.

   Each of the boxes a_(far|mid|near)(Lo|Hi)Box will be at most 1 face wide.

   a_farLoBox, a_midLoBox, a_nearLoBox will all be defined or all be undefined.
   a_farHiBox, a_midHiBox, a_nearHiBox will all be defined or all be undefined.

   Output flags:
   * a_hasLo:  1 or 0, set to whether a_(far|mid|near)loBox defined or not;
   * a_hasHi:  1 or 0, set to whether a_(far|mid|near)HiBox defined or not.
 */
extern void loHiCenterFace6(Box&                 a_farLoBox,
                            Box&                 a_midLoBox,
                            Box&                 a_nearLoBox,
                            int&                 a_hasLo,
                            Box&                 a_farHiBox,
                            Box&                 a_midHiBox,
                            Box&                 a_nearHiBox,
                            int&                 a_hasHi,
                            Box&                 a_centerBox,
                            Box&                 a_midCenterBox,
                            Box&                 a_innerCenterBox,
                            Box&                 a_entireBox,
                            const Box&           a_inBox,
                            const ProblemDomain& a_domain,
                            const int&           a_dir);

#include "NamespaceFooter.H"
#endif
