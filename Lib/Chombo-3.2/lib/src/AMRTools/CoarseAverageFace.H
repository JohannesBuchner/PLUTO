#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COARSEAVERAGEFACE_H_
#define _COARSEAVERAGEFACE_H_

#include "REAL.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "NamespaceHeader.H"


/// replaces face-centered coarse-level data w/ averaged fine-level data
/** This class replaces face-centered data on a coarse level of refinement
    with the average of the finer-level data which overlays the face.
    This class is similar to CoarseAverage */
class CoarseAverageFace
{
public:
  ///
  enum averageType
  {
    arithmetic = 0,
    harmonic,
    NUM_AVERAGE_TYPES
  };

  /// Default constructor
  CoarseAverageFace();

  /// defining constructor
  CoarseAverageFace(const DisjointBoxLayout& a_fineGrids,
                    int a_nComp, int a_nRef);

  /// destructor
  ~CoarseAverageFace();

  /// defines the object
  void define(const DisjointBoxLayout& a_fineGrids,
              int a_nComp, int a_nRef);

  ///
  bool isDefined() const;

  /// averages fine-level data to internal coarse representation of fine grid
  void average(const LevelData<FluxBox>& a_fineData);

  /// harmonic averaging of fine-level data to internal coarse representation of
  /// fine grid
  void averageHarmonic(const LevelData<FluxBox>& a_fineData);

  /// sums fine-level data to internal coarse representation of fine grid
  void sum(const LevelData<FluxBox>& a_fineData,
           const int a_refFactor = 1);

  /// obtain averaged results by copying to the destination
  void copyTo(LevelData<FluxBox>& a_coarseData);

  /// averages fine-level data to coarse level
  void averageToCoarse(LevelData<FluxBox>& a_coarse_data,
                       const LevelData<FluxBox>& a_fine_data);

  /// averages fine-level data to coarse level using harmonic averaging
  void averageToCoarseHarmonic(LevelData<FluxBox>& a_coarse_data,
                               const LevelData<FluxBox>& a_fine_data);

protected:


  bool m_isDefined;
  bool m_isAveraged;
  int m_nRef;

  // work array for coarsening of fine data, same "shape" as fine data
  LevelData<FluxBox> m_coarsenedFineData;

  /// utility function to completely determine averages.
  void computeAverages(LevelData<FluxBox>& a_coarse_data,
                       const LevelData<FluxBox>& a_fine_data,
                       int a_averageType);

  /// utility for averaging fine-level data to internal coarse representation of
  /// fine grid
  void average(const LevelData<FluxBox>& a_fineData,
               const int a_averageType,
               const int a_refFactor);

  /// averages entire single grid data from fine->crse
  void averageGridData(FluxBox& a_coarsenedFine,
                       const FluxBox& a_fine,
                       int a_averageType,
                       const int a_refFactor) const;

};

#include "NamespaceFooter.H"
#endif




