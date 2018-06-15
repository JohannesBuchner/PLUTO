#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COARSEAVERAGEEDGE_H_
#define _COARSEAVERAGEEDGE_H_

#include "REAL.H"
#include "FluxBox.H"
#include "LevelData.H"

#include "UsingNamespace.H"

/// replaces edge-centered coarse-level data w/ averaged fine-level data
/** This class replaces edge-centered data on a coarse level of refinement
    with the average of the finer-level data which overlays the edge.
    This class is similar to CoarseAverage */
class CoarseAverageEdge
{
public:
  /// Default constructor
  CoarseAverageEdge();

  /// defining constructor
  CoarseAverageEdge(const DisjointBoxLayout& a_fineGrids,
                    int a_nComp, int a_nRef);

  /// destructor
  ~CoarseAverageEdge();

  /// defines the object
  void define(const DisjointBoxLayout& a_fineGrids,
              int a_nComp, int a_nRef);

  ///
  bool isDefined() const;

  /// averages fine-level data to coarse level
  void averageToCoarse(LevelData<FluxBox>& a_coarse_data,
                       const LevelData<FluxBox>& a_fine_data);

protected:
  bool m_isDefined;

  int m_nRef;

  // work array for coarsening of fine data, same "shape" as fine data
  LevelData<FluxBox> m_coarsenedFineData;

  /// averages entire single grid data from fine->crse
  void averageGridData(FluxBox& a_coarsenedFine, const FluxBox& fine) const;
};

#endif
