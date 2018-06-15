#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBNormalizeByVolumeFraction.H"
#include "EBArith.H"
#include "EBISLayout.H"
#include "NamespaceHeader.H"
#include <vector>
using namespace std;


//----------------------------------------------------------------------------
EBNormalizeByVolumeFraction::
EBNormalizeByVolumeFraction(const EBLevelGrid& a_levelGrid):
   m_levelGrid(a_levelGrid)
{
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
EBNormalizeByVolumeFraction::
~EBNormalizeByVolumeFraction()
{
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
void
EBNormalizeByVolumeFraction::
operator()(LevelData<EBCellFAB>& a_Q,
           const Interval& a_compInterval,
           const int& a_radius) const
{
  CH_TIME("EBNormalizer::operator()");
   // Endpoints of the given interval.
   int begin = a_compInterval.begin(), end = a_compInterval.end(),
       length = a_compInterval.size();

   // Loop over the EBISBoxes within our grid. The EB data structures are
   // indexed in the same manner as the non-EB data structures, so we piggy-
   // back the former on the latter.
   a_Q.exchange();
   EBISLayout ebisLayout = m_levelGrid.getEBISL();
   DisjointBoxLayout layout = m_levelGrid.getDBL();
   for (DataIterator dit = layout.dataIterator(); dit.ok(); ++dit)
   {
      const EBISBox& box = ebisLayout[dit()];
      EBCellFAB& QFAB = a_Q[dit()];

      // Go over the irregular cells in this box.
      const IntVectSet& cfivs = (*m_levelGrid.getCFIVS())[dit()];
      const IntVectSet& irregCells = box.getIrregIVS(layout[dit()]);

      // The average has to be computed from the uncorrected data from all
      // the neighbors, so we can't apply the corrections in place. For now,
      // we stash them in a map.
      map<VolIndex, vector<Real> > correctedValues;
      for (VoFIterator vit(irregCells, box.getEBGraph()); vit.ok(); ++vit)
        {
          Real kappajSum = 0.0;
          vector<Real> kappajQjSum(length, 0.0);

          // Get all of the indices of the VoFs within a monotone path
          // radius of 1.
          VolIndex vofi = vit();
          Vector<VolIndex> vofjs;
          EBArith::getAllVoFsInMonotonePath(vofjs, vofi, box, a_radius);

          // Accumulate the contributions from the neighboring cells.
          for (unsigned int j = 0; j < vofjs.size(); ++j)
            {
              VolIndex vofj = vofjs[j];

              if (!cfivs.contains(vofj.gridIndex()))
                {
                  Real kappaj = box.volFrac(vofj);
                  for (int icomp = begin; icomp <= end; ++icomp)
                    {
                      kappajQjSum[icomp] += QFAB(vofj, icomp);
                    }
                  // Add this volume fraction to the sum.
                  kappajSum += kappaj;
                }
            }

          if (kappajSum > 0.)
            {
              // Normalize the quantity and stow it.
              vector<Real> correctedValue(length);
              //         Real kappai = box.volFrac(vofi);  //unused dtg
              for (int icomp = begin; icomp <= end; ++icomp)
                {
                  // correctedValue[icomp - begin] =
                  //    QFAB(vofi, icomp) + (1.0 - kappai) * kappajQjSum[icomp] / kappajSum;
                  correctedValue[icomp - begin] = kappajQjSum[icomp] / kappajSum;
                }
              correctedValues[vofi] = correctedValue;
            }
        }

      // Apply the corrections.
      for (map<VolIndex, vector<Real> >::const_iterator
           cit = correctedValues.begin(); cit != correctedValues.end(); ++cit)
      {
         for (int icomp = begin; icomp <= end; ++icomp)
         {
            QFAB(cit->first, icomp) = cit->second[icomp-begin];
         }
      }
   }
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
void
EBNormalizeByVolumeFraction::
operator()(LevelData<EBCellFAB>& a_Q,
           const int& a_radius) const
{
  return (*this)(a_Q, Interval(0, a_Q.nComp()-1), a_radius);
}
//----------------------------------------------------------------------------


#include "NamespaceFooter.H"
