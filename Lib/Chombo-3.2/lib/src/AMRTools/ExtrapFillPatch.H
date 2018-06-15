#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EXTRAPFILLPATCH_H_
#define _EXTRAPFILLPATCH_H_

#include "Interval.H"
#include "SPACE.H"
#include "FArrayBox.H"
#include "LayoutData.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "IntVectSet.H"
#include "NamespaceHeader.H"
class Box;

//
/// fills outer ghost cells by copying from inner ghost cells.
/**
   The picture below shows existing ghost cell data in a border around
   the valid region (including corners).  These values are
   extrapolated to the ghost cells in the outer "flaps."  The
   coordinate direction of extrapolation is a parameter (in this
   example, it is 0).

<PRE>
                  +-----+---------------------------------+-----+
                  |     | +-----------------------------+ |     |
                  |     | |                             | |     |
                  |     | |                             | |     |
                  |     | |                             | |     |
                  |     | |    InteriorBox              | |     |
                  |     | |                             |%| x   |
                  |     | |                             | |/    |
                  |  x  |%|                             | /     |
                  | /   | |                             |/|     |
                  |/    | |                             / |     |
                  /     | |                 -----------/| |     |
                 /|     | |                /            | |     |
                / |     | |               /             | |     |
               /  |     | |              /              | |     |
              /   |     | +-------------/---------------+x|     |
             /    +-----+--------------/----------------/-+-----+
            /                         /                /
           /                         /                /
          /                         /                /
         /                         /                /
  ------/------------------------ /                /
  | Ghost cells in outer flaps  |/                /
  | (locations x) are filled by |                /
  | extrapolating from existing |   ------------/------------------------
  | ghost cells (location %).   |   | Border of exisiting ghost cells.  |
  |-----------------------------|   | Assumed to be filled correctly    |
                                    | on input.                         |
                                    |-----------------------------------|
</PRE>

Limitation on grid configuration

The picture below shows the ghost cells of grid 1 for extrapolation in
direction 0.  The cells marked X were filled before extrapolation.
The cells marked > are extrapolated to from the first X to their
right.  The cells marked ? want to be filled from an X to the right,
but there isn't any.  This algorithm is not defined for such a grid
configuration, in which the outline of the level has a "step" of width
(in this case 2) less than the extrapolation radius (in this case 4).
Using this algorithm with grids that violate this condition may result
in undefined behavior.

<PRE>
            + - - - - - - - - - - - +
             > > > X X X X X X < < <
            |       +-------+       |
             > > > X|       |X < < <
            |       |       |       |
             > > > X|       |X < < <
            |       | grid 1|       |
             > > > X|       |X < < <
            |       |       |       |
             > > > X|       |X < < <
            |   +---+-------+       |
             ? ?|^ ^ ^ ^ ^ ^|X < < <
            + - - - - - - - - - - - +
                |           |
                |           |
                |  grid 0   |
                |           |
                |           |
                +-----------+
</PRE>
*/
class ExtrapFillPatch
{
public:

///
/**
   Default constructor.  User must subsequently call define().
*/
  ExtrapFillPatch();

///
/**
   Destructor.
*/
  ~ExtrapFillPatch();

///
/**
   Defining constructor.  Constructs a valid ExtrapFillPatch object.
   Equivalent to default construction followed by define().  It is a
   fatal error if any grid in a_level_domain has any dimension less
   than a_extrap_interval.end().

   There is also a restriction on grids with "jagged steps" which is
   described in "Limitation on grid configuration" in the class
   documentation above (with the ASCII art).  It is not checked
   whether the input domain complies with this restriction.  Using
   this algorithm with grids that violate this condition may result in
   undefined behavior.

   {\bf Arguments:}\\
   a_level_domain (not modified): domain of the level. \\
   a_problem_domain (not modified): problem domain at this level. \\
   a_extrap_interval (not modified): the interval of ghost cells to fill by extrapolation. \\

*/
  ExtrapFillPatch(
    const DisjointBoxLayout& a_level_domain,
    const Box& a_problem_domain,
    const Interval& a_extrap_interval
    );

///
/**
   Defining constructor.  Constructs a valid ExtrapFillPatch object.
   Equivalent to default construction followed by define().  It is a
   fatal error if any grid in a_level_domain has any dimension less
   than a_extrap_interval.end().

   There is also a restriction on grids with "jagged steps" which is
   described in "Limitation on grid configuration" in the class
   documentation above (with the ASCII art).  It is not checked
   whether the input domain complies with this restriction.  Using
   this algorithm with grids that violate this condition may result in
   undefined behavior.

   {\bf Arguments:}\\
   a_level_domain (not modified): domain of the level. \\
   a_problem_domain (not modified): problem domain at this level. \\
   a_extrap_interval (not modified): the interval of ghost cells to fill by extrapolation. \\

*/
  ExtrapFillPatch(
    const DisjointBoxLayout& a_level_domain,
    const ProblemDomain& a_problem_domain,
    const Interval& a_extrap_interval
    );

///
/**
   Defines this ExtrapFillPatch.  Existing definition is overridden.
   The user may call define() once and call fillInterp() multiple
   times with different valid data sets.

   It is a fatal error if any grid in a_level_domain has any dimension
   less than a_extrap_interval.end().

   There is also a restriction on grids with "jagged steps" which is
   described in "Limitation on grid configuration" in the class
   documentation above (with the ASCII art).  It is not checked
   whether the input domain complies with this restriction.  Using
   this algorithm with grids that violate this condition may result in
   undefined behavior.

   {\bf Arguments:}\\
   a_level_domain (not modified): domain of the level. \\
   a_problem_domain (not modified): problem domain at this level. \\
   a_extrap_interval (not modified): the interval of ghost cells to fill by extrapolation. \\

   {\bf This:}\\
   ---This object is modified.---
*/
  void
  define(
    const DisjointBoxLayout& a_level_domain,
    const Box& a_problem_domain,
    const Interval& a_extrap_interval
    );

///
/**
   Defines this ExtrapFillPatch.  Existing definition is overridden.
   The user may call define() once and call fillInterp() multiple
   times with different valid data sets.

   It is a fatal error if any grid in a_level_domain has any dimension
   less than a_extrap_interval.end().

   There is also a restriction on grids with "jagged steps" which is
   described in "Limitation on grid configuration" in the class
   documentation above (with the ASCII art).  It is not checked
   whether the input domain complies with this restriction.  Using
   this algorithm with grids that violate this condition may result in
   undefined behavior.

   {\bf Arguments:}\\
   a_level_domain (not modified): domain of the level. \\
   a_problem_domain (not modified): problem domain at this level. \\
   a_extrap_interval (not modified): the interval of ghost cells to fill by extrapolation. \\

   {\bf This:}\\
   ---This object is modified.---
*/
  void
  define(
    const DisjointBoxLayout& a_level_domain,
    const ProblemDomain& a_problem_domain,
    const Interval& a_extrap_interval
    );

///
/**
   Returns true if this object was created with the defining
   constructor or if define() has called.

   {\bf This:}\\
   This object is not modified.
*/
  bool
  isDefined() const;

///
/**
   Fills ghost cells in the outer flaps, by interpolation from data in
   the inner ghost cells.  It is an error to call if not this->isDefined().\\

   {\bf Arguments:}\\
   a_data (modified): data at this level, both the source and destination for extrapolation. \\
   a_dir (not modified): zero-based coordinate direction of extrapolation. \\
   a_dest_comp (not modified): the first component for extrapolation. \\
   a_num_comp (not modified): the number of components to
   extrapolate. \\

   {\bf This:}\\
   This object is not modified.

*/
  void
  fillExtrap(
    LevelData<FArrayBox>& a_data,
    int a_dir,
    int a_dest_comp,
    int a_num_comp
    );

  //  void
  //  printIntVectSets() const;

protected:
  bool m_is_defined;
// locations of ghost cells on the low side of the grid to extrapolate to.
  LayoutData<IntVectSet> m_lo_extrap[SpaceDim];
// locations of ghost cells on the high side of the grid to extrapolate to.
  LayoutData<IntVectSet> m_hi_extrap[SpaceDim];
// interval of ghost cells to extrapolate, where 0 are the ghost cells
// directly adjacent to the valid domain, and increasing outwards.
  Interval m_extrap_interval;
};

#include "NamespaceFooter.H"
#endif
