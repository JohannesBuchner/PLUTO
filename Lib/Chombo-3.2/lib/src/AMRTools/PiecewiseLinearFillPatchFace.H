#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PIECEWISELINEARFILLPATCHFACE_H_
#define _PIECEWISELINEARFILLPATCHFACE_H_

#include <iostream>
#include <fstream>
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "IntVectSet.H"
#include "ProblemDomain.H"

#include "NamespaceHeader.H"

/// Fills ghost cells by linear interpolation in space and time

/**
Fills some fine level ghost "faces" by piecewise linear interpolation
from the coarse level.

This class fills the first a_interp_radius layers of fine level
ghost cells that are not part of other fine level grids.  It uses
piecewise linear interpolation from the coarse level.  The slopes
are computed using van leer limiting if there is enough room for
the stencil.  It drops order and uses first-order one-sided
differences otherwise.


*/

class PiecewiseLinearFillPatchFace
{
public:
  ///
  /**
    Default constructor.  User must subsequently call define().
    */
  PiecewiseLinearFillPatchFace();

  ///
  /**
    Destructor.
    */
  ~PiecewiseLinearFillPatchFace();

  ///
  /**
    Defining constructor.

    {\bf Arguments:}\\
    a_fine_domain (not modified): domain of fine level. \\
    a_coarse_domain (not modified): domain of coarse level. \\
    a_num_comps (not modified): number of components of state vector. \\
    a_crse_problem_domain (not modified): problem domain on the coarse level.\\
    a_ref_ratio (not modified): refinement ratio. \\
    a_interp_radius (not modified): number of layers of fine ghost cells to fill by interpolation. \\
    */
  PiecewiseLinearFillPatchFace(const DisjointBoxLayout& a_fine_domain,
                               const DisjointBoxLayout& a_coarse_domain,
                               int a_num_comps,
                               const Box& a_crse_problem_domain,
                               int a_ref_ratio,
                               int a_interp_radius
                               );


  ///
  /**
    Defining constructor.

    {\bf Arguments:}\\
    a_fine_domain (not modified): domain of fine level. \\
    a_coarse_domain (not modified): domain of coarse level. \\
    a_num_comps (not modified): number of components of state vector. \\
    a_crse_problem_domain (not modified): problem domain on the coarse level.\\
    a_ref_ratio (not modified): refinement ratio. \\
    a_interp_radius (not modified): number of layers of fine ghost cells to fill by interpolation. \\
    */
  PiecewiseLinearFillPatchFace(const DisjointBoxLayout& a_fine_domain,
                               const DisjointBoxLayout& a_coarse_domain,
                               int a_num_comps,
                               const ProblemDomain& a_crse_problem_domain,
                               int a_ref_ratio,
                               int a_interp_radius
                               );

  ///
  /**
     Defines this object.  The user may call define() once and call
     fillInterp() multiple times with different valid data sets.

     {\bf Arguments:}\\
     a_fine_domain (not modified): domain of fine level. \\
     a_coarse_domain (not modified): domain of coarse level. \\
     a_num_comps (not modified): number of components of state vector. \\
     a_crse_problem_domain (not modified): problem domain on the coarse level.\\
     a_ref_ratio (not modified): refinement ratio. \\
     a_interp_radius (not modified): number of layers of fine ghost cells to fill by interpolation. \\

     {\bf This:}\\
     ---This object is modified.---

  */
  void
  define(const DisjointBoxLayout& a_fine_domain,
         const DisjointBoxLayout& a_coarse_domain,
         int a_num_comps,
         const Box& a_crse_problem_domain,
         int a_ref_ratio,
         int a_interp_radius
         );


  ///
  /**
     Defines this object.  The user may call define() once and call
     fillInterp() multiple times with different valid data sets.

     {\bf Arguments:}\\
     a_fine_domain (not modified): domain of fine level. \\
     a_coarse_domain (not modified): domain of coarse level. \\
     a_num_comps (not modified): number of components of state vector. \\
     a_crse_problem_domain (not modified): problem domain on the coarse level.\\
     a_ref_ratio (not modified): refinement ratio. \\
     a_interp_radius (not modified): number of layers of fine ghost cells to fill by interpolation. \\

     {\bf This:}\\
     ---This object is modified.---

  */
  void
  define(const DisjointBoxLayout& a_fine_domain,
         const DisjointBoxLayout& a_coarse_domain,
         int a_num_comps,
         const ProblemDomain& a_crse_problem_domain,
         int a_ref_ratio,
         int a_interp_radius
         );

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
     This is an adaptation of the Cell-centered PiecewiseLinearFillPatch
     class for face-centered data.


     Fills the first m_interp_radius layers of fine ghost cells by
     interpolation from the coarse level.  It is an error to call if not
     this->isDefined().  The range components to interpolate must be
     specified.  The corresponding components on the coarse and fine
     levels may be different.  It is required that a_fine_data's domain
     is the same as was specified in the most recent call to define().
     It is expected that the coarse and fine level's domains are
     properly nested.

     {\bf Arguments:}\\
     a_fine_data (modified): fine-level data being interpolated to.\\
     a_old_coarse_data (not modified): coarse level source data at the old time.\\
     a_new_coarse_data (not modified): coarse level source data at the new time.\\
     a_time_interp_coef (not modified): time interpolation coefficient, in the range [0:1].  0=old time, 1=new time.\\
     a_src_comp (not modifed):  starting coarse data component.\\
     a_dest_comp (not modifed):  starting fine data component.\\
     a_num_comp (not modified): number of data components to be
     interpolated.

     {\bf This:}\\
     Well, it's complicated.  As far as the user is concerned, this object
     is not modified.  See the design document if you care for details.

  */
  void
  fillInterp(LevelData<FluxBox>& a_fine_data,
             const LevelData<FluxBox>& a_old_coarse_data,
             const LevelData<FluxBox>& a_new_coarse_data,
             Real a_time_interp_coef,
             int a_src_comp,
             int a_dest_comp,
             int a_num_comp
             );

  // debugging utilities
  void
  printIntVectSets() const;


protected:
  // copy coarse data to coarsened fine work array and interpolate to
  // fine time level
  void
  timeInterp(LevelData<FluxBox>& m_coarsened_fine_data,
             const LevelData<FluxBox>& a_old_coarse_data,
             const LevelData<FluxBox>& a_new_coarse_data,
             Real a_time_interp_coef,
             int a_src_comp,
             int a_dest_comp,
             int a_num_comp
             );


  // fill the fine interpolation sites piecewise-constantly
  void
  fillConstantInterp(LevelData<FluxBox>& a_fine_data,
                     const LevelData<FluxBox>& m_coarsened_fine_data,
                     int a_src_comp,
                     int a_dest_comp,
                     int a_num_comp
                     )
    const;

  // compute slopes in specified direction
  void
  computeSlopes(FArrayBox& slope_fab,
                const FArrayBox& data_fab,
                const IntVectSet& local_centered_interp,
                const IntVectSet& local_lo_interp,
                const IntVectSet& local_hi_interp,
                int a_dir,
                int a_src_comp,
                int a_num_comp) const;

  // // compute slopes in specified direction
  // void
  // computeSlopes(int a_dir,
  //               int a_src_comp,
  //               int a_num_comp);


  // increment the fine interpolation sites with linear term for the
  // specified coordinate direction
  void
  incrementLinearInterpTangential(LevelData<FluxBox>& a_fine_data,
                     const LevelData<FluxBox>& m_coarsened_fine_data,
                                  int a_dir,
                                  int a_src_comp,
                                  int a_dest_comp,
                                  int a_num_comp)
    const;


  // now do the fine interpolation normal to the face directions
  void
  incrementLinearInterpNormal(LevelData<FluxBox>& a_fine_data,
                              int a_src_comp,
                              int a_dest_comp,
                              int a_num_comp)
    const;




protected:
  bool m_is_defined;
  // the radius of the interpolation stencil.  e.g. a stencil using
  // (i-1,j), (i,j) and (i+1,j) has a radius of 1.
  static const int s_stencil_radius;
  // refinement ratio
  int m_ref_ratio;
  // number of layers of fine ghost cells to fill by interpolation.
  int m_interp_radius;
  // // work array for coarse data in grids shaped like the fine level.
  // LevelData<FluxBox> m_coarsened_fine_data;
  // // work array for slopes
  // LevelData<FluxBox> m_slopes;
  // (FM added) coarse slope/ghost/coarsened_fine_domain
  IntVect m_coarse_slope, m_coarse_ghost;
  DisjointBoxLayout m_coarsened_fine_domain;
  // problem domain on the coarse level.
  ProblemDomain m_crse_problem_domain;
  // per-grid fine locations that you interpolate to.
  LayoutData<IntVectSet> m_fine_interp[SpaceDim];
  // per-grid coarse locations that you interpolate from, by type of
  // interpolation in the specified coordinate direction.
  // note that there are two subscripts here -- first one
  // refers to direction of slopes, second to direction of face
  LayoutData<IntVectSet> m_coarse_centered_interp[SpaceDim][SpaceDim];
  LayoutData<IntVectSet> m_coarse_lo_interp[SpaceDim][SpaceDim];
  LayoutData<IntVectSet> m_coarse_hi_interp[SpaceDim][SpaceDim];
};


#include "NamespaceFooter.H"

#endif


