#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LoHiCenter.H"
#include "NamespaceHeader.H"

///
/**
   This function is used when in direction a_dir a 3-point stencil of
   cell-centered data is being used to compute something at the
   center of the central cell of the stencil.

   Inputs:
   * a_dir is the direction of the 3-point stencil;
   * a_inBox is the cell-centered box on which we have data;
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

   Output flags:
   * a_hasLo:  1 or 0, according to whether a_loBox is defined or not;
   * a_hasHi:  1 or 0, according to whether a_hiBox is defined or not.
 */
void loHiCenter(Box&                 a_loBox,
                int&                 a_hasLo,
                Box&                 a_hiBox,
                int&                 a_hasHi,
                Box&                 a_centerBox,
                Box&                 a_entireBox,
                const Box&           a_inBox,
                const ProblemDomain& a_domain,
                const int&           a_dir)
{
  // Make a copy of the input box which can be modified
  Box inBox = a_inBox;
  inBox &= a_domain;

  // The centered difference box is always one smaller in a_dir
  a_centerBox = inBox;
  a_centerBox.grow(a_dir,-1);

  // The union of all the output boxes start off equal to the center
  // difference portion on the input box (intersected with the domain)
  a_entireBox = a_centerBox;

  // See if this chops off the high side of the input box
  Box tmp = a_inBox;
  tmp.shift(a_dir,-1);
  tmp &= a_domain;
  tmp.shift(a_dir,1);

  // If so, set up the high, one-sided difference box, a_hiBox, and expand
  // the entire box to include it
  if (!a_domain.contains(tmp))
    {
      a_hasHi = 1;
      tmp.shift(a_dir,-2);
      a_hiBox = adjCellHi(tmp,a_dir);
      a_entireBox.growHi(a_dir,1);
    }
  else
    {
      a_hasHi = 0;
      a_hiBox = Box();
    }

  // See if this chops off the low side of the input box
  tmp = a_inBox;
  tmp.shift(a_dir,1);
  tmp &= a_domain;
  tmp.shift(a_dir,-1);

  // If so, set up the low, one-sided difference box, a_loBox, and expand
  // the entire box to include it
  if (!a_domain.contains(tmp))
    {
      a_hasLo = 1;
      tmp.shift(a_dir,2);
      a_loBox = adjCellLo(tmp,a_dir);
      a_entireBox.growLo(a_dir,1);
    }
  else
    {
      a_hasLo = 0;
      a_loBox = Box();
    }

  // Make some simple sanity checks
  CH_assert(a_entireBox.contains(a_centerBox));

  if (a_hasLo == 1)
    {
      CH_assert(a_entireBox.contains(a_loBox));
    }

  if (a_hasHi == 1)
    {
      CH_assert(a_entireBox.contains(a_hiBox));
    }
}

///
/**
   This function is used when in direction a_dir a 2-point stencil of
   cell-centered data is being used to compute something at the
   face between cells of the stencil.

   Inputs:
   * a_dir is the direction of the 2-point stencil;
   * a_inBox is the cell-centered box on which we have data;
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
void loHiCenterFace(Box&                 a_loBox,
                    int&                 a_hasLo,
                    Box&                 a_hiBox,
                    int&                 a_hasHi,
                    Box&                 a_centerBox,
                    Box&                 a_entireBox,
                    const Box&           a_inBox,
                    const ProblemDomain& a_domain,
                    const int&           a_dir)
{
  // Make a copy of the input box which can be modified
  Box inBox = a_inBox;
  inBox &= a_domain;

  // The centered difference box is always one smaller in a_dir
  a_centerBox = inBox;
  a_centerBox.surroundingNodes(a_dir);
  a_centerBox.grow(a_dir,-1);

  // The union of all the output boxes start off equal to the center
  // difference portion on the input box (intersected with the domain)
  a_entireBox = a_centerBox;

  // See if this chops off the high side of the input box
  Box tmp = a_inBox;
  tmp.shift(a_dir,-1);
  tmp &= a_domain;
  tmp.shift(a_dir,1);

  // If so, set up the high, one-sided difference box, a_hiBox, and expand
  // the entire box to include it
  if (!a_domain.contains(tmp))
    {
      a_hasHi = 1;
      tmp.shift(a_dir,-2);
      a_hiBox = adjCellHi(tmp,a_dir);
      a_hiBox.shiftHalf(a_dir,1);
      a_entireBox.growHi(a_dir,1);
    }
  else
    {
      a_hasHi = 0;
      a_hiBox = Box().convert(BASISV(a_dir));
    }

  // See if this chops off the low side of the input box
  tmp = a_inBox;
  tmp.shift(a_dir,1);
  tmp &= a_domain;
  tmp.shift(a_dir,-1);

  // If so, set up the low, one-sided difference box, a_loBox, and expand
  // the entire box to include it
  if (!a_domain.contains(tmp))
    {
      a_hasLo = 1;
      tmp.shift(a_dir,2);
      a_loBox = adjCellLo(tmp,a_dir);
      a_loBox.shiftHalf(a_dir,-1);
      a_entireBox.growLo(a_dir,1);
    }
  else
    {
      a_hasLo = 0;
      a_loBox = Box().convert(BASISV(a_dir));
    }

  // Make some simple sanity checks
  CH_assert(a_entireBox.contains(a_centerBox));

  if (a_hasLo == 1)
    {
      CH_assert(a_entireBox.contains(a_loBox));
    }

  if (a_hasHi == 1)
    {
      CH_assert(a_entireBox.contains(a_hiBox));
    }
}

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
   * a_centerBox is where the regular 4-point stencil can be used;
   * a_hiBox is where a stencil must be used on the high side
     with all four points lower;
   * a_nextHiBox is where a stencil must be used on the high side
     with one point higher and three points lower;
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
void loHiCenterFace4(Box&                 a_loBox,
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
                     const int&           a_dir)
{
  loHiCenterFace(a_loBox, a_hasLo, a_hiBox, a_hasHi,
                 a_centerBox, a_entireBox, a_inBox, a_domain, a_dir);
  a_innerCenterBox = a_centerBox;
  if (a_hasLo)
    {
      a_nextLoBox = a_loBox;
      a_nextLoBox.shift(a_dir, 1);
      // Remove first row of faces from a_centerBox.
      int centerMin = a_centerBox.smallEnd(a_dir);
      a_innerCenterBox.setSmall(a_dir, centerMin + 1);
    }
  if (a_hasHi)
    {
      a_nextHiBox = a_hiBox;
      a_nextHiBox.shift(a_dir, -1);
      // Remove last row of faces from a_centerBox.
      int centerMax = a_centerBox.bigEnd(a_dir);
      a_innerCenterBox.setBig(a_dir, centerMax - 1);
    }
}

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
void loHiCenter5(Box&                 a_loBox,
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
                 const int&           a_dir)
{
  loHiCenter(a_loBox, a_hasLo, a_hiBox, a_hasHi,
             a_centerBox, a_entireBox, a_inBox, a_domain, a_dir);
  a_innerCenterBox = a_centerBox;
  if (a_hasLo)
    {
      a_nextLoBox = a_loBox;
      a_nextLoBox.shift(a_dir, 1);
      // Remove first row of faces from a_centerBox.
      int centerMin = a_centerBox.smallEnd(a_dir);
      a_innerCenterBox.setSmall(a_dir, centerMin + 1);
    }
  if (a_hasHi)
    {
      a_nextHiBox = a_hiBox;
      a_nextHiBox.shift(a_dir, -1);
      // Remove last row of faces from a_centerBox.
      int centerMax = a_centerBox.bigEnd(a_dir);
      a_innerCenterBox.setBig(a_dir, centerMax - 1);
    }
}

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
void loHiCenterFace6(Box&                 a_farLoBox,
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
                     const int&           a_dir)
{
  loHiCenterFace4(a_farLoBox, a_midLoBox, a_hasLo,
                  a_farHiBox, a_midHiBox, a_hasHi,
                  a_centerBox, a_midCenterBox, 
                  a_entireBox, a_inBox, a_domain, a_dir);
  // We'll shave off the ends of a_innerCenterBox where a_hasLo or a_hasHi.
  a_innerCenterBox = a_midCenterBox;
  if (a_hasLo)
    {
      // a_nearLoBox has not been set yet: it is a_midLoBox shifted up by 1
      a_nearLoBox = a_midLoBox;
      a_nearLoBox.shift(a_dir, 1);
      // Remove first row of faces from a_innerCenterBox.
      int midCenterMin = a_midCenterBox.smallEnd(a_dir);
      a_innerCenterBox.setSmall(a_dir, midCenterMin + 1);
    }
  if (a_hasHi)
    {
      // a_nearHiBox has not been set yet: it is a_midHiBox shifted down by 1
      a_nearHiBox = a_midHiBox;
      a_nearHiBox.shift(a_dir, -1);
      // Remove last row of faces from a_innerCenterBox.
      int midCenterMax = a_midCenterBox.bigEnd(a_dir);
      a_innerCenterBox.setBig(a_dir, midCenterMax - 1);
    }
}

#include "NamespaceFooter.H"
