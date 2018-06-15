#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FINEINTERP_H_
#define _FINEINTERP_H_

#include "REAL.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "FArrayBox.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

class DisjointBoxLayout;

/// replaces fine level data with interpolation of coarse level data.
/**
   This class replaces data at a fine level of refinement with data
   interpolated from a coarser level of refinement.  Interpolation is
   piecewise bi(tri)linear, with van Leer slopes if there is room for
   the stencil, with lower-order slopes if there isn't.  See the
   design document.

 */
class FineInterp
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  FineInterp();

  ///
  /**
     Destructor.
  */
  ~FineInterp();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\
     a_fine_problem_domain (not modified): problem domain at the fine level.\\

  */
  FineInterp(const DisjointBoxLayout& a_fine_domain,
             const int& a_numcomps,
             const int& a_ref_ratio,
             const Box& a_fine_problem_domain);

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\
     a_fine_problem_domain (not modified): problem domain at the fine level.\\
  */
  FineInterp(const DisjointBoxLayout& a_fine_domain,
             const int& a_numcomps,
             const int& a_ref_ratio,
             const ProblemDomain& a_fine_problem_domain);

  ///
  /**
     Defines this object.  Existing information is overriden.

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\
     a_fine_problem_domain (not modified): problem domain at the fine level.\\

     {\bf This:}\\
     ---This object is modified.---

  */
  void
  define(const DisjointBoxLayout& a_fine_domain, // the fine level domain
         const int& a_numcomps,                   // the number of components
         const int& a_ref_ratio,                  // the refinement ratio
         const Box& a_fine_problem_domain);      // problem domain

  ///
  /**
     Defines this object.  Existing information is overriden.

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\
     a_fine_problem_domain (not modified): problem domain at the fine level.\\

     {\bf This:}\\
     ---This object is modified.---

  */
  void
  define(const DisjointBoxLayout& a_fine_domain, // the fine level domain
         const int& a_numcomps,                   // the number of components
         const int& a_ref_ratio,                  // the refinement ratio
         const ProblemDomain& a_fine_problem_domain);

  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.

     {\bf This:}\\
     This object is not modified.
  */
  bool
  isDefined() const;

  ///
  /**
     Replaces a_fine_data with data interpolated from a_coarse_data. It
     is an error to call if not this->isDefined().  The domain of
     a_fine_data should be the same as the fine domain specified in the
     most recent call to define().  It is expected that the coarse and
     fine level's domains are properly nested.  Both a_coarse_data and
     a_fine_data should have the same number of components specified in
     the most recent call to define().

     {\bf Arguments:}\\
     a_fine_data (modified): fine data. \\
     a_coarse_data (not modified): coarse data. \\
     a_averageFromDest: if true, first average data from a_fine_data down 
                        to the resolution of a_coarse_data, then interp 
                        everything back up -- necessary when the coarse 
                        grids don't cover the fine grid (i.e when flattening
                        an AMR hierarchy to a single resolution). Default is 
                        false.


     {\bf This:}\\
     Well, it's complicated.  As far as the user is concerned, this object
     is not modified.  See the design document if you care for details.

  */
  void
  interpToFine(LevelData<FArrayBox>& a_fine_data,
               const LevelData<FArrayBox>& a_coarse_data,
               bool a_averageFromDest=false);

  /// Just do piecewise-constant interpolation.
  void
  pwcinterpToFine(LevelData<FArrayBox>& a_fine_data,
                  const LevelData<FArrayBox>& a_coarse_data,
                  bool a_averageFromDest = false);

protected:
  void
  interpGridData(BaseFab<Real>& a_fine,
                 const BaseFab<Real>& a_coarse,
                 const Box& a_coarsened_fine_box,
                 int a_ref_ratio)
    const;
  void
  pwcinterpGridData(BaseFab<Real>& a_fine,
                    const BaseFab<Real>& a_coarse,
                    const Box& a_coarsened_fine_box,
                    int a_ref_ratio) const;

protected:
  bool is_defined;
  // the refinement ratio
  int m_ref_ratio;
  // work array for the coarse level data in a domain that is the
  // outline of the fine level domain on the coarse level
  LevelData<FArrayBox> m_coarsened_fine_data;
  // coarse level problem domain
  ProblemDomain m_coarse_problem_domain;
};

#include "NamespaceFooter.H"
#endif
