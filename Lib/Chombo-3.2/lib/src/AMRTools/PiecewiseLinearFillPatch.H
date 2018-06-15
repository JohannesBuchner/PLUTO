#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PIECEWISELINEARFILLPATCH_H_
#define _PIECEWISELINEARFILLPATCH_H_

#include <iostream>
#include <fstream>
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "IntVectSet.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

/// Fills ghost cells by linear interpolation in space and time

/**
Fills some fine level ghost cells by piecewise linear interpolation
from the coarse level.

This class also presents an interface for performing piecewise
constant interpolation in space.  If this is all that is required,
define with 'a_pwconst_interp_only = true' to save memory and
expense.

This class fills the first a_interp_radius layers of fine level
ghost cells that are not part of other fine level grids.  It uses
piecewise linear interpolation from the coarse level.  The slopes
are computed using van leer limiting if there is enough room for
the stencil.  It drops order and uses first-order one-sided
differences otherwise.

Below is a picture of a fine level, with a_interp_radius = 1.  The
ghost cells of grid 1 are shown.  Cells marked X are interpolated
to from the coarse level.  Cells marked ^ are not interpolated to
because they are in the valid domain of another fine grid.  Cells
marked . are not interpolated to because they are outside the
problem domain.  All of this changes in the presence of periodic
boundary conditions -- in periodic directions, all cells are considered
to be within the domain.

<PRE>
 +=======================================+
 |                                       |
 |                                       |
 |                                       |
 |                                       |
 |         + - - - - - - - - - - - - - - | +
 |          X X X X X X X X X X X X X X X|.
 |         | +===========================+ |
 |          X|                           |.
 |         | |                           | |
 |          X|        grid 1             |.
 |         | |                           | |
 |          X|                           |.
 |         | |                           | |
 |          X|                           |.
 |   +=======+===+===============+=======+ |
 |   |      ^ ^ ^|X X X X X X X X|^ ^ ^ ^|.
 |   |     + - - | - - - - - - - | - - - | +
 |   |           |               |grid 2 |
 |   |  grid 0   |               +=======+
 |   |           |                       |
 |   |           |                       |
 |   |           |                       |
 |   +===========+                       |
 |                                       |
 |                                       |
 |              problem domain           |
 +=======================================+
</PRE>

 the picture below shows the locations of the types of slopes.  The
 coarse grid is shown, and the projection of the fine grid onto the
 coarse index space (on which slopes are computed).  Coarse cells
 marked ^ and . do not need slopes.  Van Leer slopes in both
 directions are used at coarse cells marked X.  Coarse cells marked
 V and < have one-sided differencing: cells marked V have low-sided
 differencing in the y-direction (in this example, because the
 coarse grid doesn't exist on the high side of these cells), and
 cells marked < have low-sided differencing in the x-direction
 (because the problem domain doesn't exist on the high side of
 these cells).  There can also be cells that have one-sided
 differences in both directions.

<PRE>
 +-----------------------+===============+
 |                       | coarse grid 1 |
 |                       |               |
 |                       |               |
 +=======================+ - - - - - - - | - +
 |        \ / \ / \ / \ /|\ / \ / \ / \ /|
 |       | V   V   V   V | X   X   X   < | . |
 |        / \ / \ / \ / \|/ \ / \ / \ / \|
 |       |   +---------------------------+   |
 |        \ /|           |               |
 |       | X |           |               | . |
 |        / \|           |               |
 |       |   |           |               |   |
 |        \ /|           |               |
 |       | X |           |               | . |
 |        / \|           |               |
 |   +-------+---+---------------+-------+   |
 |   |           |\ / \ /|\ / \ /|       |
 |   |   | ^   ^ | X   X | X   X | ^   ^ | . |
 |   |           |/ \ / \|/ \ / \|       |
 |   |   + - - - | - - - - - - - +-------+ - +
 |   |           |       |               |
 |   |           |       |               |
 |   |           |       |               |
 |   +-----------+       +===============+
 |                       |               |
 |                       |               |
 |      coarse grid 0    |               |
 +=======================+---------------+
</PRE>
*/

class PiecewiseLinearFillPatch
{
public:
  ///
  /**
    Default constructor.  User must subsequently call define().
    */
  PiecewiseLinearFillPatch();

  ///
  /**
    Destructor.
    */
virtual  ~PiecewiseLinearFillPatch();

  ///
  /**
    Defining constructor.

    {\bf Arguments:}\\
    a_fine_domain (not modified): domain of fine level. \\
    a_coarse_domain (not modified): domain of coarse level. \\
    a_num_comps (not modified): number of components of state vector. \\
    a_problem_domain (not modified): problem domain on the coarse level. \\
    a_ref_ratio (not modified): refinement ratio. \\
    a_interp_radius (not modified): number of layers of fine ghost cells to fill by interpolation. \\
    a_constInterpOnly (not modified): if true, only set up for piecewise-constant interpolation in space. \\
    */
  PiecewiseLinearFillPatch(const DisjointBoxLayout& a_fine_domain,
                           const DisjointBoxLayout& a_coarse_domain,
                           int a_num_comps,
                           const Box& a_crse_problem_domain,
                           int a_ref_ratio,
                           int a_interp_radius,
                           bool a_pwconst_interp_only = false
                           );

  ///
  /**
    Defining constructor.

    {\bf Arguments:}\\
    a_fine_domain (not modified): domain of fine level. \\
    a_coarse_domain (not modified): domain of coarse level. \\
    a_num_comps (not modified): number of components of state vector. \\
    a_problem_domain (not modified): problem domain on the coarse level. \\
    a_ref_ratio (not modified): refinement ratio. \\
    a_interp_radius (not modified): number of layers of fine ghost cells to fill by interpolation. \\
    a_constInterpOnly (not modified): if true, only set up for piecewise-constant interpolation in space. \\
    */
  PiecewiseLinearFillPatch(const DisjointBoxLayout& a_fine_domain,
                           const DisjointBoxLayout& a_coarse_domain,
                           int a_num_comps,
                           const ProblemDomain& a_crse_problem_domain,
                           int a_ref_ratio,
                           int a_interp_radius,
                           bool a_pwconst_interp_only = false
                           );

  ///
  /**
    Defines this object.  The user may call define() once and call
    fillInterp() multiple times with different valid data sets.

    {\bf Arguments:}\\
    a_fine_domain (not modified): domain of fine level. \\
    a_coarse_domain (not modified): domain of coarse level. \\
    a_num_comps (not modified): number of components of state vector. \\
    a_problem_domain (not modified): problem domain on the coarse level. \\
    a_ref_ratio (not modified): refinement ratio. \\
    a_interp_radius (not modified): number of layers of fine ghost cells to fill by interpolation. \\
    a_constInterpOnly (not modified): if true, only set up for piecewise-constant interpolation in space. \\

    {\bf This:}\\
    ---This object is modified.---

    */
  void
  define(const DisjointBoxLayout& a_fine_domain,
         const DisjointBoxLayout& a_coarse_domain,
         int a_num_comps,
         const Box& a_crse_problem_domain,
         int a_ref_ratio,
         int a_interp_radius,
         bool a_pwconst_interp_only = false
         );

  ///
  /**
    Defines this object.  The user may call define() once and call
    fillInterp() multiple times with different valid data sets.

    {\bf Arguments:}\\
    a_fine_domain (not modified): domain of fine level. \\
    a_coarse_domain (not modified): domain of coarse level. \\
    a_num_comps (not modified): number of components of state vector. \\
    a_problem_domain (not modified): problem domain on the coarse level. \\
    a_ref_ratio (not modified): refinement ratio. \\
    a_interp_radius (not modified): number of layers of fine ghost cells to fill by interpolation. \\
    a_constInterpOnly (not modified): if true, only set up for piecewise-constant interpolation in space. \\

    {\bf This:}\\
    ---This object is modified.---

    */
  void
  define(const DisjointBoxLayout& a_fine_domain,
         const DisjointBoxLayout& a_coarse_domain,
         int a_num_comps,
         const ProblemDomain& a_crse_problem_domain,
         int a_ref_ratio,
         int a_interp_radius,
         bool a_pwconst_interp_only = false
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
  fillInterp(LevelData<FArrayBox>& a_fine_data,
             const LevelData<FArrayBox>& a_old_coarse_data,
             const LevelData<FArrayBox>& a_new_coarse_data,
             Real a_time_interp_coef,
             int a_src_comp,
             int a_dest_comp,
             int a_num_comp
             );

  ///
  /**
    Fills the first m_interp_radius layers of fine ghost cells by
    piece-wise constant interpolation (only in space) from a coarse level
    {\bf Arguments:}\\
    a_fine_data (modified): fine-level data being interpolated to.\\
    a_coarse_data (not modified): coarse level source data.\\
    a_src_comp (not modifed):  starting coarse data component.\\
    a_dest_comp (not modifed):  starting fine data component.\\
    a_num_comp (not modified): number of data components to be
  */
  void
  fillInterpPWConstSpace(LevelData<FArrayBox>& a_fine_data,
                         const LevelData<FArrayBox>& a_coarse_data,
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
  timeInterp(const LevelData<FArrayBox>& a_old_coarse_data,
             const LevelData<FArrayBox>& a_new_coarse_data,
             Real a_time_interp_coef,
             int a_src_comp,
             int a_dest_comp,
             int a_num_comp
             );

  // fill the fine interpolation sites piecewise-constantly
virtual  void
  fillConstantInterp(LevelData<FArrayBox>& a_fine_data,
                     int a_src_comp,
                     int a_dest_comp,
                     int a_num_comp
                     )
    const;

  // compute slopes in specified direction
virtual void  computeSlopes(int a_src_comp,
                            int a_num_comp);

  void computeSimpleSlopesFab(FArrayBox & a_slopeFab,
                              const int       & a_src_comp,
                              const int       & a_num_comp,
                              const int       & a_dir,
                              const FArrayBox & a_dataFab,
                              const IntVectSet& a_local_centered_interp,
                              const IntVectSet& a_local_lo_interp,
                              const IntVectSet& a_local_hi_interp);

  void  computeMultiDimSlopes(FArrayBox      & a_slopes0,
                              FArrayBox      & a_slopes1,
                              FArrayBox      & a_slopes2,
                              const FArrayBox& a_dataFab,
                              const int      & a_src_comp,
                              const int      & a_num_comp,
                              const Box      & a_slopeBox);

  // increment the fine interpolation sites with linear term for the
  // specified coordinate direction
virtual  void
  incrementLinearInterp(LevelData<FArrayBox>& a_fine_data,
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
  // work array for coarse data in grids shaped like the fine level.
  LevelData<FArrayBox> m_coarsened_fine_data;
  // work array for slopes
  LevelData<FArrayBox> m_slopes[3];
  // problem domain on the coarse level.
  ProblemDomain m_crse_problem_domain;
  // per-grid fine locations that you interpolate to.
  LayoutData<IntVectSet> m_fine_interp;
  // per-grid coarse locations that you interpolate from, by type of
  // interpolation in the specified coordinate direction.
  LayoutData<IntVectSet> m_coarse_centered_interp[SpaceDim];
  LayoutData<IntVectSet> m_coarse_lo_interp[SpaceDim];
  LayoutData<IntVectSet> m_coarse_hi_interp[SpaceDim];
  // cached Copier
  Copier m_coarsenCopier;
};

#include "NamespaceFooter.H"
#endif
