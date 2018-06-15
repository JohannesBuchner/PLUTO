#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBEllipticLoHiCenter.H"
#include "NamespaceHeader.H"

// This function is used when in direction a_dir a 3 point stencil of cell-
// centered data is being used to compute something at the cell center of the
// center cell of the stencil.  The data for the stencil is valid in
// a_inputBox.  It divides a_inputBox into a cell-centered box where the full
// stencil can be used, a_centerBox, and cell-centered boxes where 2 point
// stencils can be used, a_loBox and a_hiBox, based on the current problem
// domain, a_domain, and the stencil direction, a_dir.  The union of these 2
// and 3 point stencil boxes is returned as a_entireBox.  The 2 point stencil
// boxes are one wide, at most, and if they have been defined then the
// corresponding flag, a_hasLo or a_hasHi, is set to one, otherwise these
// flags are zero.  All output boxes lie within a_resultBox and the domain.

void ebEllipticLoHiCenter(Box&                 a_loBox,
                          int&                 a_hasLo,
                          Box&                 a_hiBox,
                          int&                 a_hasHi,
                          Box&                 a_centerBox,
                          Box&                 a_entireBox,
                          const Box&           a_resultBox,
                          const Box&           a_inputBox,
                          const ProblemDomain& a_domain,
                          const int&           a_dir)
{
  // Make a copy of the input box which can be modified
  Box inBox = a_inputBox;
  inBox &= a_domain;

  // The centered difference box is always one smaller in a_dir
  a_centerBox = inBox;
  a_centerBox.grow(a_dir,-1);

  // The union of all the output boxes start off equal to the center
  // difference portion on the input box (intersected with the domain)
  a_entireBox = a_centerBox;

  // See if this chops off the high side of the input box
  Box tmp = inBox;
  tmp.shift(a_dir,1);
  tmp &= a_domain;
  tmp.shift(a_dir,-1);

  // If so, set up the high, one-sided difference box, a_hiBox, and expand
  // the entire box to include it
  if (!tmp.isEmpty() && tmp != inBox)
  {
    a_hiBox = adjCellHi(tmp,a_dir);
    a_hiBox &= a_resultBox;

    // If the result is empty then there is no a_hiBox
    if (a_hiBox.isEmpty())
    {
      a_hasHi = 0;
    }
    else
    {
      a_hasHi = 1;

      // If the entire box is empty, "growHi" doesn't do anything so this is
      // necessary
      if (a_entireBox.isEmpty())
      {
        a_entireBox = a_hiBox;
      }
      else
      {
        a_entireBox.growHi(a_dir,1);
      }
    }
  }
  else
  {
    a_hasHi = 0;
  }

  // See if this chops off the low side of the input box
  tmp = inBox;
  tmp.shift(a_dir,-1);
  tmp &= a_domain;
  tmp.shift(a_dir,1);

  // If so, set up the low, one-sided difference box, a_loBox, and expand
  // the entire box to include it
  if (!tmp.isEmpty() && tmp != inBox)
  {
    a_loBox = adjCellLo(tmp,a_dir);
    a_loBox &= a_resultBox;

    // If the result is empty then there is no a_loBox
    if (a_loBox.isEmpty())
    {
      a_hasLo = 0;
    }
    else
    {
      a_hasLo = 1;

      // If the entire box is empty, "growLo" doesn't do anything so this is
      // necessary
      if (a_entireBox.isEmpty())
      {
        a_entireBox = a_loBox;
      }
      else
      {
        a_entireBox.growLo(a_dir,1);
      }
    }

  }
  else
  {
    a_hasLo = 0;
  }

  a_centerBox &= a_resultBox;
  a_entireBox &= a_resultBox;

  // Make some simple sanity checks
  CH_assert(a_entireBox.contains(a_centerBox));

  if (a_hasHi == 1)
  {
    CH_assert(a_entireBox.contains(a_hiBox));
    CH_assert(!a_hiBox.isEmpty());
  }

  if (a_hasLo == 1)
  {
    CH_assert(a_entireBox.contains(a_loBox));
    CH_assert(!a_loBox.isEmpty());
  }
}

// This function is used when in direction a_dir a 2 point stencil of cell-
// centered data is being used to compute something at the cell face between
// the cell centers of the stencil.  The data for the stencil is valid in
// a_inputBox.  It uses a_inputBox to compute a box (face-centered in a_dir)
// where the full stencil can be used, a_centerBox, and boxes (face-centered
// in a_dir) where a 1 point stencil can be used, a_loBox and a_hiBox based
// on the current problem domain, a_domain, and the stencil direction, a_dir.
// The union of these 1 and 2 point stencil boxes is returned as a_entireBox
// (face-centered in a_dir).  The 1 point stencil boxes are one wide, at most,
// and if they have been defined then the corresponding flag, a_hasLo or
// a_hasHi, is set to one, otherwise these flags are zero.  All output boxes
// lie within a_resultBox and the domain.

void ebEllipticLoHiCenterFace(Box&                 a_loBox,
                              int&                 a_hasLo,
                              Box&                 a_hiBox,
                              int&                 a_hasHi,
                              Box&                 a_centerBox,
                              Box&                 a_entireBox,
                              const Box&           a_resultBox,
                              const Box&           a_inputBox,
                              const ProblemDomain& a_domain,
                              const int&           a_dir)
{
  // Make a copy of the input box which can be modified
  Box inBox = a_inputBox;
  inBox &= a_domain;

  // The centered difference box is always one smaller in a_dir
  a_centerBox = inBox;
  a_centerBox.surroundingNodes(a_dir);
  a_centerBox.grow(a_dir,-1);

  // The union of all the output boxes start off equal to the center
  // difference portion on the input box (intersected with the domain)
  a_entireBox = a_centerBox;

  // See if this chops off the high side of the input box
  Box tmp = inBox;
  tmp.shift(a_dir,1);
  tmp &= a_domain;
  tmp.shift(a_dir,-1);

  // If so, set up the high, one-sided difference box, a_hiBox, and expand
  // the entire box to include it
  if (!tmp.isEmpty() && tmp != inBox)
  {
    a_hasHi = 1;

    a_hiBox = adjCellHi(tmp,a_dir);
    a_hiBox.shiftHalf(a_dir,1);
    a_hiBox &= a_resultBox;

    // If the result is empty then there is no a_hiBox
    if (a_hiBox.isEmpty())
    {
      a_hasHi = 0;
    }
    else
    {
      a_hasHi = 1;

      // If the entire box is empty, "growHi" doesn't do anything so this is
      // necessary
      if (a_entireBox.isEmpty())
      {
        a_entireBox = a_hiBox;
      }
      else
      {
        a_entireBox.growHi(a_dir,1);
      }
    }
  }
  else
  {
    a_hasHi = 0;
  }

  // See if this chops off the low side of the input box
  tmp = inBox;
  tmp.shift(a_dir,-1);
  tmp &= a_domain;
  tmp.shift(a_dir,1);

  // If so, set up the low, one-sided difference box, a_loBox, and expand
  // the entire box to include it
  if (!tmp.isEmpty() && tmp != inBox)
  {
    a_hasLo = 1;

    a_loBox = adjCellLo(tmp,a_dir);
    a_loBox.shiftHalf(a_dir,-1);
    a_loBox &= a_resultBox;

    // If the result is empty then there is no a_loBox
    if (a_loBox.isEmpty())
    {
      a_hasLo = 0;
    }
    else
    {
      a_hasLo = 1;

      // If the entire box is empty, "growLo" doesn't do anything so this is
      // necessary
      if (a_entireBox.isEmpty())
      {
        a_entireBox = a_loBox;
      }
      else
      {
        a_entireBox.growLo(a_dir,1);
      }
    }
  }
  else
  {
    a_hasLo = 0;
  }

  a_centerBox &= a_resultBox;
  a_entireBox &= a_resultBox;

  // Make some simple sanity checks
  CH_assert(a_entireBox.contains(a_centerBox));

  if (a_hasHi == 1)
  {
    CH_assert(a_entireBox.contains(a_hiBox));
    CH_assert(!a_hiBox.isEmpty());
  }

  if (a_hasLo == 1)
  {
    CH_assert(a_entireBox.contains(a_loBox));
    CH_assert(!a_loBox.isEmpty());
  }
}

void ebEllipticLoHiCenter(Vector<Box>&         a_loBox,
                          Vector<Box>&         a_hiBox,
                          Box&                 a_centerBox,
                          Vector<Box>&         a_overlap,
                          const Box&           a_dblBox,
                          const ProblemDomain& a_domain)
{
  CH_assert(a_loBox.size()==CH_SPACEDIM);
  CH_assert(a_hiBox.size()==CH_SPACEDIM);

  a_centerBox = a_dblBox;

  const Box& dbox = a_domain.domainBox();
  for (int i=0; i<CH_SPACEDIM; i++)
    {
      if (dbox.smallEnd(i) == a_centerBox.smallEnd(i)) //ignore periodic for now
        {
          a_centerBox.growLo(i,-1);
          a_loBox[i] = a_dblBox;
          a_loBox[i].growHi(i,-(a_dblBox.size(i)-1));
        }
      else
        {
          a_loBox[i] = Box();
        }
      if (dbox.bigEnd(i) == a_centerBox.bigEnd(i)) //ignore periodic for now
        {
          a_centerBox.growHi(i,-1);
          a_hiBox[i] = a_dblBox;
          a_hiBox[i].growLo(i,-(a_dblBox.size(i)-1));
        }
      else
        {
          a_hiBox[i] = Box();
        }
    }
  for (int i=0; i<CH_SPACEDIM; ++i)
  {
    const Box& lo = a_loBox[i];
    const Box& hi = a_hiBox[i];
    if (!lo.isEmpty())
    {
      for (int j=i+1; j<CH_SPACEDIM; ++j)
        {
          const Box& loj = a_loBox[j];
          const Box& hij = a_hiBox[j];
          if (lo.intersects(loj))
          {
            a_overlap.push_back(lo & loj);
          }
          if (lo.intersects(hij))
          {
            a_overlap.push_back(lo & hij);
          }
        }
    }
    if (!hi.isEmpty())
    {
      for (int j=i+1; j<CH_SPACEDIM; ++j)
      {
        const Box& loj = a_loBox[j];
        const Box& hij = a_hiBox[j];
        if (hi.intersects(loj))
        {
          a_overlap.push_back(hi & loj);
        }
        if (hi.intersects(hij))
        {
          a_overlap.push_back(hi & hij);
        }
      }
    }
  }
}

#include "NamespaceFooter.H"
