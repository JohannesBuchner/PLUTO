#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COARSEAVERAGECODIM_H_
#define _COARSEAVERAGECODIM_H_

#include "CodimBox.H"
#include "LevelData.H"
#include "NamespaceHeader.H"


/*******************************************************************************
 */
/// Replaces coarse-level data with averaged fine-level data on the
/// codimensions of a box
/** This class is similar to CoarseAverage
 *
 ******************************************************************************/

class CoarseAverageCodim
{

public:

  /// Default constructor
  CoarseAverageCodim();

  /// Full constructor
  CoarseAverageCodim(const DisjointBoxLayout& a_fineGrids,
                     const int                a_codim,
                     const int                a_nComp,
                     const int                a_nRef);

  /// Destructor
  ~CoarseAverageCodim();

  /// Defines the object
  void define(const DisjointBoxLayout& a_fineGrids,
              const int                a_codim,
              const int                a_nComp,
              const int                a_nRef);

  /// Has define been called?
  bool isDefined() const;

  /// Averages fine-level data to internal coarse representation of fine grid
  void average(const LevelData<CodimBox<FArrayBox> >& a_fineData);

  /// Harmonic averaging of fine data to coarse representation of fine grid
  void averageHarmonic(const LevelData<CodimBox<FArrayBox> >& a_fineData);

  /// Summation of fine data to internal coarse representation of fine grid
  void sum(const LevelData<CodimBox<FArrayBox> >& a_fineData,
           const int                  a_refFactor = 1);

  /// Obtain averaged results by copying to the destination
  void copyTo(LevelData<CodimBox<FArrayBox> >& a_coarseData);

  /// Averages fine-level data to coarse level
  void averageToCoarse(LevelData<CodimBox<FArrayBox> >&       a_coarseData,
                       const LevelData<CodimBox<FArrayBox> >& a_fineData);

  /// Averages fine-level data to coarse level using harmonic averaging
  void averageToCoarseHarmonic(LevelData<CodimBox<FArrayBox> >&       a_coarseData,
                               const LevelData<CodimBox<FArrayBox> >& a_fineData);

protected:

  enum AverageType
  {
    arithmetic = 0,
    harmonic
  };

  LevelData<CodimBox<FArrayBox> > m_coarsenedFineData;
                                      ///< Work array for coarsening of fine
                                      ///< data, same "shape" as fine data

  bool m_isDefined;
  bool m_isAveraged;
  int m_nRef;

  /// Utility function to completely determine averages.
  void computeAverages(LevelData<CodimBox<FArrayBox> >&       a_coarseData,
                       const LevelData<CodimBox<FArrayBox> >& a_fineData,
                       const AverageType          a_averageType);

  /// Utility for averaging fine-level data to internal coarse representation of
  /// fine grid
  void average(const LevelData<CodimBox<FArrayBox> >& a_fineData,
               const AverageType          a_averageType,
               const int                  a_refFactor);

  /// Averages entire single grid data from fine->coarse
  void averageGridData(CodimBox<FArrayBox>&         a_coarsenedFine,
                       const CodimBox<FArrayBox>&   a_fine,
                       const AverageType a_averageType,
                       const int         a_refFactor) const;

};


/*******************************************************************************
 *
 * Class CoarseAverageCodim: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
/// Has define been called?
/*--------------------------------------------------------------------*/

inline bool
CoarseAverageCodim::isDefined() const
{
  return m_isDefined;
}

#include "NamespaceFooter.H"
#endif
