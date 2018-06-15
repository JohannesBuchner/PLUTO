#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/*****************/
/*****************/

#include "timeInterp.H"
#include "NamespaceHeader.H"

void
timeInterp(LevelData<FArrayBox>& a_dest, Real a_time,
           const LevelData<FArrayBox>& a_old_phi, Real a_old_time,
           const LevelData<FArrayBox>& a_new_phi, Real a_new_time,
           const Interval& a_comps)
{
  timeInterp(a_dest, a_time, a_old_phi, a_old_time,
             a_new_phi, a_new_time, a_comps, a_comps);
}

void
timeInterp(LevelData<FArrayBox>& a_dest, Real a_time,
           const LevelData<FArrayBox>& a_old_phi, Real a_old_time,
           const LevelData<FArrayBox>& a_new_phi, Real a_new_time,
           const Interval& a_src_comps, const Interval& a_dest_comps)
{

  int ncomp = a_src_comps.size();
  CH_assert (a_old_phi.nComp() > a_src_comps.end());
  CH_assert (a_new_phi.nComp() > a_src_comps.end());
  CH_assert (a_dest.nComp() > a_dest_comps.end());

  Real alpha = (a_time - a_old_time)/(a_new_time-a_old_time);

  const DisjointBoxLayout& grids = a_dest.getBoxes();
  LevelData<FArrayBox> tempData(grids,ncomp);

  // this is potentially a little inefficient in a parallel
  // world -- might want to do interpolation on source
  // grids, then copy to new grids?  (one copyTo rather than two)
  // on the other hand, i doubt that's really a factor in the current code
  a_old_phi.copyTo(a_src_comps, a_dest, a_dest_comps);
  Interval tempComps(0, ncomp-1);
  a_new_phi.copyTo(a_src_comps, tempData, tempComps);

  // now loop over grids, scale old and new contributions, and add
  DataIterator dit = a_dest.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
    {
      FArrayBox& oldData = a_dest[dit()];
      FArrayBox& newData = tempData[dit()];

      oldData *= (1.0-alpha);
      newData *= alpha;

      oldData.plus(newData,0,a_dest_comps.begin(), ncomp);
    }

}

#include "NamespaceFooter.H"
