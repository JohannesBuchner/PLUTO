#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BRMESHREFINE_H_
#define _BRMESHREFINE_H_

#include <climits>

#include "Vector.H"
#include "Box.H"
#include "IntVectSet.H"
#include "REAL.H"
#include "MeshRefine.H"
#include "Interval.H"
#include <list>
#include "NamespaceHeader.H"

// Constants:

// Minimum acceptable ratio of tagged cells to total cells for
// the Berger-Rigoutsos algorithm in \func{makeBoxes}.
// Used as default for \var{FillRatio} optional argument.
#ifndef _BR_MIN_BOX_FILL_RATIO_
#define _BR_MIN_BOX_FILL_RATIO_ ( 0.75 )
#endif

/// Class which manages Berger-Rigoutsos grid generation
/**
   This class manages grid generation from sets of tagged cells using the
   Berger-Rigoutsos algorithm in the context of the MeshRefine class
   from which it is derived

   There are two ways grids can be defined based on tagged cells.
   one takes a single IntVectSet of tags defined on the BaseLevel
   mesh and uses that set of tags for every level to be refined;
   the other takes a Vector<IntVectSet> of tags defined on all the
   mesh levels to be refined and uses those.

   <b>  Long Description: </b>

   Create new meshes based on tagged cells on a range of levels of a mesh
   hierarchy.  Each level of tagged cells is used to generate a new mesh at
   the next finer level.  The finest level in the output mesh will be one
   level higher than the top of the range of levels given as input.  As a
   special case, use the same tags (appropriately refined) for all levels.

   \b Usage:

   Call the regrid functions after computing error estimates and tagging cells.
   To add a new mesh level, set TopLevel to the index of the finest level
   in the existing mesh and define tags on the finest level.  To keep the
   existing number of mesh levels, set TopLevel to one less than the
   index of the finest level and don't define any tags on the finest level.
   If a single IntVectSet of tags is passed (instead of a
   Vector<IntVectSet>) then the same tags (properly refined) will be used
   for all the new meshes up to level TopLevel+1.  In any case, the
   meshes at levels BaseLevel and below are not modified.  The output
   argument newmeshes will be reallocated to the necessary size before
   being used.  When this function returns, the elements of the
   newmeshes vector corresponding to the unchanged levels will be
   filled in with copies of the levels from the old mesh vector.  The
   variable tags is modified in an undefined way, so its contents
   should not be relied upon.  The variable BlockFactor specifies the
   amount by which each box will be coarsenable. Every grid box will
   have an integral multiple of BlockFactor cells in each dimension and
   also lower index values that are integral multiples.  As a side effect,
   the minimum box size will be BlockFactor.

   Expensive validations are done only when debugging is enabled
   (i.e. the DEBUG make variable is "TRUE").

   <b> Usage Notes: </b>

   All the input vectors should be defined with max index >= TopLevel.
   They should have values for indices [BaseLevel:TopLevel].
   (except for OldMeshes, which must be defined for all indices).  The
   new mesh vector newmeshes will be redefined up to index
   TopLevel+1.  RefRatios should be defined such that
   RefRatios[L] is the value to use to refine the level L mesh to
   produce the level L+1 mesh.  The tags vector is modified in an
   undefined manner.  The output variable newmeshes may not be
   completely defined if an exception occurs.
   The BlockFactor can be used to force a minimum box size.
*/
class BRMeshRefine : public MeshRefine
{
public:
  /// Default constructor -- leaves object in an unusable state
  BRMeshRefine();

  /// Full constructor -- leaves object in usable state
  BRMeshRefine(/// Level 0 domain
               const Box& a_baseDomain,
               /// Refinement ratios -- refRatio[0] is btwn levels 0 and 1
               const Vector<int>& a_refRatios,
               /// Measure of how efficiently tagged cells will be covered
               const Real a_fillRatio,
               /// Amount by which grids are guaranteed to be coarsenable
               const int a_blockFactor,
               /// Proper nesting buffer amount
               const int a_bufferSize,
               /// Maximum grid length in any direction -- 0 means no limit.
               const int a_maxSize);

  /// Full constructor -- leaves object in usable state
  BRMeshRefine(/// Level 0 domain
               const ProblemDomain& a_baseDomain,
               /// Refinement ratios -- refRatio[0] is btwn levels 0 and 1
               const Vector<int>& a_refRatios,
               /// Measure of how efficiently tagged cells will be covered
               const Real a_fillRatio,
               /// Amount by which grids are guaranteed to be coarsenable
               const int a_blockFactor,
               /// Proper nesting buffer amount
               const int a_bufferSize,
               /// Maximum grid length in any direction -- 0 means no limit.
               const int a_maxSize);

  /// Destructor
  virtual ~BRMeshRefine();

  /// Define function -- size of RefRatios will define maximum number of levels
  void define(/// Level 0 domain
              const Box& a_baseDomain,
              /// Refinement ratios -- refRatio[0] is btwn levels 0 and 1
              const Vector<int>& a_refRatios,
              /// Measure of how efficiently tagged cells will be covered
              const Real a_fillRatio,
              /// Amount by which grids are guaranteed to be coarsenable
              const int a_blockFactor,
              /// Proper nesting buffer amount
              const int a_bufferSize,
              /// Maximum grid length in any direction -- 0 means no limit
              const int a_maxSize);

  /// Define function -- size of RefRatios will define maximum number of levels
  void define(/// :evel 0 domain
              const ProblemDomain& a_baseDomain,
              /// Refinement ratios -- refRatio[0] is btwn levels 0 and 1
              const Vector<int>& a_refRatios,
              /// Measure of how efficiently tagged cells will be covered
              const Real a_fillRatio,
              /// Amount by which grids are guaranteed to be coarsenable
              const int a_blockFactor,
              /// Proper nesting buffer amount
              const int a_bufferSize,
              /// Maximum grid length in any direction -- 0 means no limit
              const int a_maxSize);

  /// Constructs a set of boxes which covers a set of tagged cells
  /** Constructs a set of boxes which covers a set of tagged cells
      by using the Berger-Rigoutsos algorithm.  Everything should
      be on the same level, and blocking factor is not applied.
      Boxes will be on the same refinement level as the tags.
      This would normally be a protected function, but it can be useful
      to call it on it's own, so it has been left public.
  */
  void
  makeBoxes(/// Putput: refined boxes at each new level
            Vector<Box>&      a_mesh,
            /// Input: set of tagged cells to cover
            const IntVectSet& a_tags,
            /// Input: proper nesting domain in which mesh boxes must live
            const IntVectSet& a_pnd,
            /// Input: physical domain
            const ProblemDomain& a_domain,
            /// Input: largest number of cells in any dimension for any box
            const int         a_maxSize,
            const int         a_totalBufferSize) const;


protected:
  /**
     Function which actually implement Berger-Rigoutsos chopping.
  */
  void
  makeBoxes(/// Output: refined boxes at each new level
            std::list<Box>&      a_mesh,
            /// Input: set of tagged cells to cover
            IntVectSet& a_tags,
            /// Input: proper nesting domain in which mesh boxes must live
            const IntVectSet& a_pnd,
            /// Input: physical domain
            const ProblemDomain& a_domain,
            /// Input: largest number of cells in any dimension for any box
            const int         a_maxSize,
            /// Input: depth of this recursion in the algorithm
            const int         a_depth,
            const int         a_totalBufferSize
            ) const;
  /**
   */
  void
  splitBox(std::list<Box> & a_boxes ,
           const std::list<Box>::iterator& a_boxindex,
           const int a_dimension ,const int a_maxboxsize ) const;

  void
  splitBox(std::list<Box> & a_boxes ,
           const std::list<Box>::iterator& a_boxindex,
           const int a_maxboxsize ) const;

  ///
  Vector<int>
  makeTrace( const IntVectSet& a_Ivs ,int a_dir ) const;

  ///
  void
  makeTraces( const IntVectSet& a_Ivs ,Vector<int>* a_traces ) const;
  ///
  int
  findSplit( const Vector<int>& a_trace ) const;

  int
  findSplit( const Vector<int>& a_trace, const int a_maxSize ) const;

  ///
  int
  findMaxInflectionPoint( const Vector<int>& a_trace ,int& a_maxVal ) const;
  int
  findMaxInflectionPoint( const Vector<int>& a_trace ,int& a_maxVal, const int a_maxSize ) const;

  ///
  void
  splitTags( const IntVectSet& a_tags,
             const int a_split_dir ,const int a_split_indx,
             IntVectSet& a_tags_lo ,IntVectSet& a_tags_hi ) const;

  void splitTagsInPlace(const int a_split_dir, const int a_split_indx,
                        IntVectSet& a_tags_inout_lo,
                        IntVectSet& a_tags_hi) const;

  void splitTagsInBestDimension(IntVectSet& a_tags_inout_lo,
                                IntVectSet& a_tags_hi,
                                const int a_maxSize) const;
  ///
  void
  breakBoxes(Vector<Box>& a_vboxin,  const int& a_maxSize,
             const int& a_idir) const;

  ///
  int
  maxloc( const int* a_V ,const int a_Size ) const;

  void makeBoxesParallel(std::list<Box>&      a_mesh,
                         IntVectSet&    a_tags,
                         const IntVectSet&    a_pnd,
                         const ProblemDomain& a_domain,
                         const int            a_maxSize,
                         const int            a_depth,
                         const int            a_totalBufferSize,
                         const int            a_minSize,
                         const Interval&      a_procInterval
                         ) const;

  void sendBoxesParallel(   const std::list<Box>& a_mesh,
                            int tag) const;

  void receiveBoxesParallel(const Interval& a_from,
                            const Interval& a_to,
                            std::list<Box>& a_mesh,
                            int tag) const;
};

/// Splits domain into vector of disjoint boxes with max size maxsize
/**
   Blocking factor is default to one.
   If you make minimum size > 1, then domain must
   be coarsenable and refineable by blockfactor
   (refine(coarsen(domain,blockfactor), minsize) == domain)
   or an error is thrown.  This would be defined in
   MeshRefine.H, except that it needs to use a BRMeshRefine object.
*/
extern void domainSplit(const ProblemDomain& a_domain, Vector<Box>& a_vbox,
                        int a_maxSize, int a_blockfactor=1);

///
/**
   Splits domain into a vector of disjoint boxes with
   maximum size maxsize.
   blocking factor is default to one.
   If you make minimum size > 1, then domain must
   be coarsenable and refineable by blockfactor
   (refine(coarsen(domain,blockfactor), minsize) == domain)
   or an error is thrown. This would be defined in
   MeshRefine.H, except that it needs to use a BRMeshRefine object.
*/
extern void domainSplit(const Box& a_domain, Vector<Box>& a_vbox,
                        int a_maxSize, int a_blockfactor=1);

///
/**
   Recursive function to enforce max size of boxes in a given direction.
   Does not call Meshrefine.
 */
extern void
breakBoxes(Vector<Box>& a_vboxin,  const int& a_maxSize, const int& a_idir);

#include "NamespaceFooter.H"
#endif
