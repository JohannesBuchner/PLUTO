#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Mask.H"
#include "LayoutIterator.H"

// -------------------------------------------------------------
void
Mask::buildMask(BaseFab<int>& a_mask, const Box& a_dProblem,
                const BoxLayout& a_grids, const BoxLayout* a_fineGridsPtr,
                int a_nRefFine)
{
  ProblemDomain physdomain(a_dProblem);
  buildMask(a_mask, physdomain, a_grids, a_fineGridsPtr,
            a_nRefFine);
}

// -------------------------------------------------------------
void
Mask::buildMask(BaseFab<int>& a_mask, const ProblemDomain& a_dProblem,
                const BoxLayout& a_grids, const BoxLayout* a_fineGridsPtr,
                int a_nRefFine)
{
  // first set entire box to Physical BC
  a_mask.setVal(maskPhysical);

  // now set all of domain interior to coarse
  Box domainInterior(a_mask.box());
  domainInterior &= a_dProblem;
  a_mask.setVal(maskCoarse,domainInterior,0);

  // now loop over this level's boxes and set them to "copy"
  LayoutIterator lit = a_grids.layoutIterator();

  for (lit.reset(); lit.ok(); ++lit)
    {
      Box intersectBox = a_grids.get(lit());
      intersectBox &= a_mask.box();
      if (!intersectBox.isEmpty())
        {
          a_mask.setVal(maskCopy,intersectBox,0);
        }
    }

  // if finer grids exist, set them to "covered"
  if (a_fineGridsPtr != NULL)
    {
      CH_assert (a_nRefFine > 1);

      LayoutIterator litFine = a_fineGridsPtr->layoutIterator();
      for (litFine.reset(); litFine.ok(); ++litFine)
        {
          Box coarsenedBox(a_fineGridsPtr->get(litFine()));
          coarsenedBox.coarsen(a_nRefFine);
          coarsenedBox &= a_mask.box();
          if (!coarsenedBox.isEmpty())
            {
              a_mask.setVal(maskCovered,coarsenedBox,0);
            }
        }
    }
}

// -------------------------------------------------------------
void
Mask::buildMasks(LevelData<BaseFab<int> >& a_masks, const Box& a_dProblem,
                 const BoxLayout& a_grids, const BoxLayout* a_fineGridsPtr,
                 int a_nRefFine)
{
  ProblemDomain physdomain(a_dProblem);
  buildMasks(a_masks, a_dProblem, a_grids, a_fineGridsPtr, a_nRefFine);
}

// -------------------------------------------------------------
void
Mask::buildMasks(LevelData<BaseFab<int> >& a_masks,
                 const ProblemDomain& a_dProblem,
                 const BoxLayout& a_grids, const BoxLayout* a_fineGridsPtr,
                 int a_nRefFine)
{
  // simply iterate over boxes and call single-basefab function
  DataIterator dit = a_masks.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
    {
      buildMask(a_masks[dit()], a_dProblem, a_grids, a_fineGridsPtr,
                a_nRefFine);
    }
}
