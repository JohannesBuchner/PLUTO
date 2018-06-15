#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MESHREFINE_H_
#define _MESHREFINE_H_

#include <climits>

#include "Vector.H"
#include "Box.H"
#include "IntVectSet.H"
#include "REAL.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

// Constants:

/// Class which manages grid generation
class MeshRefine
/**
   This class manages grid generation from sets of tagged cells. It is
   designed to be a pure virtual base class from which another class
   may be derived with a specific grid-generation algorithm (for example,
   the BRMeshRefine class).

   There are two ways grids can be defined based on tagged cells.
   one takes a single IntVectSet of tags defined on the BaseLevel
   mesh and uses that set of tags for every level to be refined;
   the other takes a Vector<IntVectSet> of tags defined on all the
   mesh levels to be refined and uses those.

   <b> Long Description: </b>

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
{
public:
  /// default constructor -- leaves object in an unusable state
  MeshRefine();

  /// full constructor -- leaves object in usable state
  MeshRefine(/// level 0 domain
               const Box& a_baseDomain,
               /// refinement ratios -- refRatio[0] is btwn levels 0 and 1
               const Vector<int>& a_refRatios,
               /// measure of how efficiently tagged cells will be covered
               const Real a_fillRatio,
               /// amount by which grids are guaranteed to be coarsenable
               const int a_blockFactor,
               /// proper nesting buffer amount
               const int a_bufferSize,
               /// maximum grid length in any direction -- 0 means no limit.
               const int a_maxSize);

  /// full constructor -- leaves object in usable state
  MeshRefine(/// level 0 domain
               const ProblemDomain& a_baseDomain,
               /// refinement ratios -- refRatio[0] is btwn levels 0 and 1
               const Vector<int>& a_refRatios,
               /// measure of how efficiently tagged cells will be covered
               const Real a_fillRatio,
               /// amount by which grids are guaranteed to be coarsenable
               const int a_blockFactor,
               /// proper nesting buffer amount
               const int a_bufferSize,
               /// maximum grid length in any direction -- 0 means no limit.
               const int a_maxSize);

  /// destructor
  virtual ~MeshRefine();

  /// define function -- size of RefRatios will define maximum number of levels
  void define(/// level 0 domain
              const Box& a_baseDomain,
              /// refinement ratios -- refRatio[0] is btwn levels 0 and 1
              const Vector<int>& a_refRatios,
              /// measure of how efficiently tagged cells will be covered
              const Real a_fillRatio,
              /// amount by which grids are guaranteed to be coarsenable
              const int a_blockFactor,
              /// proper nesting buffer amount
              const int a_bufferSize,
              /// maximum grid length in any direction -- 0 means no limit
              const int a_maxSize);

  /// define function -- size of RefRatios will define maximum number of levels
  virtual void define(/// level 0 domain
                      const ProblemDomain& a_baseDomain,
                      /// refinement ratios -- refRatio[0] is btwn levels 0 and 1
                      const Vector<int>& a_refRatios,
                      /// measure of how efficiently tagged cells will be covered
                      const Real a_fillRatio,
                      /// amount by which grids are guaranteed to be coarsenable
                      const int a_blockFactor,
                      /// proper nesting buffer amount
                      const int a_bufferSize,
                      /// maximum grid length in any direction -- 0 means no limit
                      const int a_maxSize);

  ///create hierarchy of grids from a single level of tags
  /** This function creates a hierarchy of grids from a single level of
      tags on BaseLevel.  If tags exist, then all levels will have grids.
      Returns the new finest level of grids.
  */
  virtual int regrid(/// new set of grids at every level
                     Vector<Vector<Box> >& a_newmeshes,
                     /// tagged cells on baseLevel
                     const IntVectSet& a_tags,
                     /// index of base mesh level (finest unchanged level)
                     const int a_baseLevel,
                     /// top level to refine (one less than finest possible level)
                     const int a_topLevel,
                     /// existing grids (if no previous grids, set to domains)
                     const Vector<Vector<Box> >& a_oldMeshes);

  /// create hierarchy of grids from tags at all levels
  /** This function creates a hierarchy of grids from tags at all
      refinement levels.  It is possible that not all levels will
      return with grids, since there may not be tags at all levels.
      Returns the new finest level of grids.
  */
  virtual int regrid(/// new set of grids at every level
                     Vector<Vector<Box> >& a_newmeshes,
                     /// tagged cells on each existing level
                     Vector<IntVectSet>& a_tags,
                     /// index of base mesh level (finest unchanged level)
                     const int a_baseLevel,
                     /// top level to refine (one less than finest possible level)
                     const int a_topLevel,
                     /// existing grids (if no previous grids, set to domains)
                     const Vector<Vector<Box> >& a_oldMeshes);

  // Access functions

  /// returns vector of refinement ratios
  const Vector<int>& refRatios() const;

  /// returns fillRatio
  Real fillRatio() const;

  /// returns blocking factor
  int blockFactor() const;

  /// returns proper nesting buffer size
  int bufferSize() const;

  /// returns maximum box size in any dimension -- 0 means no limit
  int maxSize() const;

  /// sets vector of refinement ratios
  void refRatios(const Vector<int>& a_nRefVect);

  /// sets fillRatio
  virtual void fillRatio(const Real a_fill_ratio);

  /// sets blocking factor
  virtual void blockFactor(const int a_block_factor);

  /// sets proper nesting buffer size
  virtual void bufferSize(const int a_buffer_size);

  /// sets maximum box size in any dimension -- 0 means no limit
  virtual void  maxSize(const int a_max_size);

  /// has this object been defined properly?
  bool isDefined() const;

  /// sets proper nesting region granularity.
  void granularity(int a_granularity);

  /// constructs a set of boxes which covers a set of tagged cells
  /** constructs a set of boxes which covers a set of tagged cells
      by using the algorithm of choice.  Everything should
      be on the same level, and blocking factor is not applied.
      Boxes will be on the same refinement level as the tags.
      This would normally be a protected function, but it can be useful
      to call it on its own, so it has been left public.
  */
  virtual void
  makeBoxes(/// output: refined boxes at each new level
            Vector<Box>&      a_mesh,
            /// input: set of tagged cells to cover
            const IntVectSet& a_tags,
            /// input: proper nesting domain in which mesh boxes must live
            const IntVectSet& a_pnd,
            /// input: physical domain
            const ProblemDomain& a_domain,
            ///input: largest number of cells in any dimension for any box
            const int         a_maxSize,
            const int         a_totalBufferSize) const = 0;

  void setPNDMode(int a_mode);

protected:

  /// computes local blockFactors used internally to enforce the BlockFactor
  /** This function computes values for m_local_blockfactors array, which
      is the amount that tags on a level are coarsened in order to guarantee
      that the grids on the next finer level are coarsenable by the
      BlockFactor. */
  virtual void
  computeLocalBlockFactors();

  virtual bool properlyNested(const Box& a_box,
                      const ProblemDomain& a_domain,
                      const IntVectSet& a_pnd,
                      int   a_totalBuffer) const;

  /// Computes proper nesting domains.
  /**
     This should only be
     called by refine. it assumes that everything has already been
     coarsened by the local blocking factor
  */
  virtual void
  makePNDs(///output: proper nesting domains at each level
           Vector<IntVectSet>& a_pnds,
           Vector<int>&        a_totalBufferSize,
           /// input: index of highest AMR level not to be refined
           const int           a_baseLevel,
           /// input: index of highest AMR level in output
           const int           a_topLevel,
           /// input: (same as in \em meshRefine)
           const Vector<ProblemDomain>&  a_domains,
           /// input: boxes at mesh level \em BaseLevel
           const IntVectSet&   a_baseMesh,
           ///input: (similar to \em meshRefine; but with level-dependent coarsening factors)
           const Vector<int>&  a_bufferSize ) const;

  virtual void
  makePNDs(///output: proper nesting domains at each level
           Vector<IntVectSet>& a_pnds,
           Vector<int>&        a_totalBufferSize,
           /// input: index of highest AMR level not to be refined
           const int           a_baseLevel,
           /// input: index of highest AMR level in output
           const int           a_topLevel,
           /// input: (same as in \em meshRefine)
           const Vector<ProblemDomain>&  a_domains,
           /// input: boxes at mesh level \em BaseLevel
           const Vector<Box>&   a_baseMesh,
           ///input: (similar to \em meshRefine; but with level-dependent coarsening factors)
           const Vector<int>&  a_bufferSize ) const;

  virtual void buildSupport(const ProblemDomain& lvldomain, Vector<Box>& lvlboxes, IntVectSet& modifiedTags);

  virtual void clipBox(Box& a_box, const ProblemDomain& a_domain) const ;
  // local data members

  bool m_isDefined;

  Vector<ProblemDomain> m_vectDomains;

  Vector<IntVectSet>    m_pnds;
  int                   m_lastBase;
  int                   m_lastTop;
  int                   m_lastBuffer;

  Vector<int> m_nRefVect;

  Real m_fillRatio;

  int m_blockFactor;

  Vector<int> m_level_blockfactors;

  int m_bufferSize;

  int m_maxSize;

  int m_granularity;

  int m_PNDMode;

};

#include "NamespaceFooter.H"
#endif
