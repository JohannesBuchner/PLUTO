#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _KAPPASQUARENORMAL_H_
#define _KAPPASQUARENORMAL_H_

#include "EBCellFAB.H"
#include "LevelData.H"
#include "EBLevelGrid.H"

#include "NamespaceHeader.H"

/** \class KappaSquareNormal
 *  This class represents an operation that takes a volume-fraction-weighted
 *  quantity kappa*Q (where kappa is the volume fraction) and computes Q
 *  from it by normalizing it w.r.t. the volume fractions of cells within a
 *  monotone path radius. The radius may either be specified explicitly or
 *  computed from a given EBLevelGrid.
 */
class KappaSquareNormal
{
   public:

   /** Construct an KappaSquareNormal object associated with the
    *  associated with the given EBLevelGrid.
    *  \param EBLevelGrid The level grid for which data can be normalized
    *           by the resulting object.
    */
   explicit KappaSquareNormal(const EBLevelGrid& a_levelGrid);

   /** Destructor. */
   ~KappaSquareNormal();

   /** Compute the normalized quantity corresponding to the given
    *  volume-fraction-weighted quantity by volume-averaging that quantity
    *  over every monotone path radius for each cell in the level.
    *  \param a_Q A volume-fraction-weighted quantity to be normalized.
    *  \param a_compInterval An interval representing the range of components
    *                        to be normalized by this operation.
    */
   void operator()(LevelData<EBCellFAB>& a_Q,
                   const Interval& a_compInterval) const;

   /** Compute the normalized quantity corresponding to the given
    *  volume-fraction-weighted quantity by volume-averaging that quantity
    *  over every monotone path radius for each cell in the level. This
    *  version of the correction normalizes all components of \a a_Q.
    *  \param a_Q A volume-fraction-weighted quantity to be normalized.
    */
   void operator()(LevelData<EBCellFAB>& a_Q) const;

   protected:

   // Forbidden operations.
   KappaSquareNormal();
   KappaSquareNormal(const KappaSquareNormal&);
   KappaSquareNormal& operator=(const KappaSquareNormal&);

   private:

   /** The associated level grid. */
   const EBLevelGrid& m_levelGrid;
};

#include "NamespaceFooter.H"
#endif
