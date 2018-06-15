#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <iostream>

#include "MeshRefine.H"
#include "BoxIterator.H"
#include "MayDay.H"
#include "SPMD.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"
using std::endl;

// base Mesh refinement class
// ---------
//  class MeshRefine
//
// Short Description:
// ------------------
//  Create new meshes based on tagged cells on a range of levels of a mesh
//  hierarchy.  Each level of tagged cells is used to generate a new mesh at
//  the next finer level.  The finest level in the output mesh will be one
//  level higher than the top of the range of levels given as input.  As a
//  special case, use the same tags (appropriately refined) for all levels.
//  MeshRefine is intended as a pure virtual base class, from which should
//  be derived a class which supplied the actual box-generation algorithm
//  (for example, the BRMeshRefine class implements the Berger-Rigoutsos
//  algorithm).
//
// Usage:
// ------
//  Call the regrid functions after computing error estimates and tagging
//  cells.  To add a new mesh level, set \var{TopLevel} to the index of
//  the finest level in the existing mesh and define tags on the finest
//   level.  To keep the existing number of mesh levels, set \var{TopLevel}
//  to one less than the index of the finest level and don't define any tags
//  on the finest level. If a single IntVectSet of tags is passed (instead
//  of a Vector<IntVectSet>) then the same tags (properly refined) will be
//  used for all the new meshes up to level \var{TopLevel}+1.  In any case,
//  the meshes at levels \var{BaseLevel} and below are not modified.  The
//  output argument \var{newmeshes} will be reallocated to the necessary
//  size before being used.  When this function returns, the elements of the
//  \var{newmeshes} vector corresponding to the unchanged levels will be
//  filled in with copies of the levels from the old mesh vector.  The
//  variable \var{tags} is modified in an undefined way, so its contents
//  should not be relied upon.  The variable \var{BlockFactor} specifies a
//  minimum grid size.  Every grid box will have an integral
//  multiple of \var{RefRatio}*\var{ceiling(BlockFactor/RefRatio)} cells
//  in each dimension and also lower index values that are integral multiples.
//
//
// Arguments:
// ----------
//  \var{newmeshes}                type: Vector<Vector<Box>>
//      (output, index[0:\var{TopLevel}+1]) vector containing the new mesh hierarchy
//
//  \var{tags}                     type: Vector<IntVectSet>
//      (input/output, index [0:\var{TopLevel}]) tagged cells to be covered by the
//      new meshes.  For any level L, the new mesh at level L+1 will try to cover
//      the tagged cells at level L.  NOTE: This variable is modified in an
//      undefined manner so its values should not be used on exit
// --or--
//  \var{tags}                     type: IntVectSet
//      (input/output) tagged cells on the \var{BaseLevel} mesh.  These tags
//      will be refined and reused for the other mesh levels.
//
//  \var{BaseLevel}                type: int
//      (input, range[0:\var{TopLevel}]) index of the finest mesh level that
//      will not be modified.
//
//  \var{TopLevel}                 type: int
//      (input, range[\var{BaseLevel}:\var{OldMeshes.size()-1]) index of the
//      finest mesh level that will be refined (i.e. the finest mesh level for
//      which tags are defined); the new mesh vector \var{newmeshes} will have
//      an element at index [\var{TopLevel}+1].
//
//  \var{Domains}                  type: Vector<ProblemDomain>
//      (input, index[0:\var{TopLevel}]) boundaries of the physical domain,
//      one mesh box at each level.  Must be consistent with \var{RefRatios}
//      (i.e. Domains[levelL+1] == Domains[levelL] * RefRatios[levelL])
//
//  \var{OldMeshes}                type: Vector<Vector<Box>>
//      (input, index[0:\var{TopLevel}]) the existing mesh hierarchy vector.
//      The \var{BaseLevel} boxes should be aligned on \var{BufferSize}-sized
//      internal boundaries.
//
//  \var{RefRatios}                type: Vector<int>&
//      (input, index[0:\var{TopLevel}], range[2^i,i>0])
//      the number of cells (in each direction) in level L+1 for each cell in
//      level L.  Must be an integral power of 2.
//
//
//  \var{FillRatio}                type: Real
//      (input, range[0:1]) the minimum acceptable ratio
//      of filled cells in a box in the Berger-Rigoutsos algorithm in
//      \func{makeBoxes} (typical value = 0.75
//
//  \var{BlockFactor}              type: int
//      (input, range[1:) )
//      factor  that determines the minimum size of a box
//      and what index values a box corner may have (i.e. what alignment is
//      has). Grids are guaranteed to be coarsenable by the BlockFactor.
//
//  \var{BufferSize}               type: int
//      (input, range[1:), default=1)
//      minimum number of level L mesh cells between the boundary of a level L+1
//      mesh and the nearest level L-1 mesh cell.  This is used to guarantee
//      that stencils can be computed across coarse-fine boundaries.
//
//  \var{MaxBoxSize}               type: int
//      (input, range[2*BufferSize:), default=\infinity)
//      largest number of grid points in any dimension that will be generated
//      in an individual grid box.  Must be at least twice the BufferSize.
//
// Returns:
// --------
//  The first two arguments (\var{tags}, \var{newmeshes}) are modified.
//
//  Expensive validations are done only when debugging is enabled
//   (i.e. the DEBUG make variable is "TRUE").
//
//
// Usage Notes:
// ------------
//  All the input vectors should be defined with max index >= \var{TopLevel}.
//  They should have values for indices [\var{BaseLevel}:\var{TopLevel}].
//  (except for \var{OldMeshes}, which must be defined for all indices).  The
//  new mesh vector \var{newmeshes} will be redefined up to index
//  \var{TopLevel}+1.  \var{RefRatios} should be defined such that
//  \var{RefRatios}[L] is the value to use to refine the level L mesh to
//  produce the level L+1 mesh.  The \var{tags} vector is modified in an
//  undefined manner.  The output variable \var{newmeshes} may not be
//  completely defined if an exception occurs.  If you do not know what the
//  optional arguments do, it is probably safe to ignore them.
//  The BlockFactor can be used to force a minimum box size.
//
// Long Description:
// -----------------
//  MeshRefine replaces the meshes at levels \var{BaseLevel}+1 to
//  \var{TopLevel} and creates a new mesh at level \var{TopLevel}+1.  Levels
//  \var{BaseLevel} and below in the new mesh hierarchy are the same as in
//  the old one.  The new meshes are generated so that the tagged cells at
//  level L are refined and included in the new mesh at level L+1 (except
//  for tagged cells that aren`t properly nested in the \var{BaseLevel}
//  mesh).  If TopLevel<0, no new meshes are created and a warning code
//  is returned.
//
// How it all works:
// -----------------
//  Really quite simple.  You create the mesh on the finest level (say l+1)
//  from tags on (l).  Then you update the tags on (l) to include all the boxes
//  you just created on (l+1).  To this, you add the buffer cells to ensure
//  proper nesting.  Now you coarsen these tags to (l-1) and union with whatever
//  else is tagged by the application at (l-1).  With these tags, you create the
//  mesh on level (l).  Note that it properly nests the mesh on (l+1).
//
//  So what is a proper nesting domain (PND) all about?
//  The problem is that you cannot change the base level mesh.  In the previous
//  algorithm we expanded the tags at level (l) to properly nest level (l+1).
//  You don't want to end up in a situation where you have to expand the base
//  mesh to nest all the grids you just created.  So instead, work in the other
//  direction:  given a base-mesh, what is the maximum size of the next finer
//  mesh that is properly nested --- this is the PND.  Using the PNDs we can
//  crop the tags on the levels to enusre that proper nesting will be enforced
//  without ever having to change the base level.  This is probably how the
//  original code worked but now, for performance reasons, enforcement of the
//  PND is embedded into the box generation algorithm and the PND itself is
//  basically described by two quantities:  a refinement of the base mesh and a
//  "total" buffer describing the distance this level mesh has to be away from
//  the base PND.
//
//  What about block factors?
//  To enforce the blocking factor, the entire algorithm is performed on a
//  special coarsening of the actual mesh -- so that when this coarsening is
//  undone, the generated grids automatically enforce the blocking factor.  See
//  computeLocalBlockFactors() for how the coarsening ratio is determined.
//
// Numerical Algorithm:
// --------------------
//  The tagged cells at level L are used to generate a new (refined) mesh at
//  level L+1.  The \var{BaseLevel} mesh is unchanged.  The new meshes must
//  satisfy the proper nesting requirements, i.e. each new mesh box must be
//  contained in the proper nesting domain (PND) for the next lower level.  The
//  PNDs of all levels are derived from the \var{BaseLevel} mesh (see
//  \func{makePNDs}).  The basic idea is that the proper nesting domain of a
//  level L mesh contains the interior cells of the level L-1 PND, except for
//  faces that are on the physical boundary.  The base level PND is defined as
//  the interiors of the base level mesh boxes, except for boundary faces
//  (again) and box face cells that are surrounded by cells of other boxes.
//
//  In the following graphic, the PND for two mesh boxes with BufferSize=1 is
//  shown.
//
//        ###################################################
//        #     |   |   |   |   |   |                       #
//        #     |   | o | o | o |   |                       #
//        #     |   |   |   |   |   |                       #
//        #     +---+---+---+---+---%---+---+---+---+---+---#
//        #     |   |   |   |   |   %   |   |   |   |   |   #
//        #     |   | o | o | o |   %   |   |   |   |   |   #
//        #     |   |   |   |   |   %   |   |   |   |   |   #
//        #     +---+---+---+---+---%---+---+---+---+---+---#
//        #     |   |   |   |   |   %   |   |   |   |   |   #
//        #     |   | o | o | o | o % o | o | o | o | o | o #
//        #     |   |   |   |   |   %   |   |   |   |   |   #
//        #     +---+---+---+---+---%---+---+---+---+---+---#
//        #     |   |   |   |   |   %   |   |   |   |   |   #
//        #     |   |   |   |   |   %   | o | o | o | o | o #
//        #     |   |   |   |   |   %   |   |   |   |   |   #
//        #     +---+---+---+---+---%---+---+---+---+---+---#
//        #                         |   |   |   |   |   |   #
//        #                         |   | o | o | o | o | o #
//        #                         |   |   |   |   |   |   #
//        ###################################################
//
//  where:
//   #    are physical boundaries (non periodic)
//   | +  are cell boundaries of the two mesh boxes
//   %    are box boundaries shared between boxes
//   o    are cells in the proper nesting domain
//
//  The tagged cells outside the proper nesting domains are discarded by
//  intersecting the set of tagged cells at each level with the PND at that
//  level.  This ensures that all cells in the refined meshes are properly
//  nested.  The actual mesh generation starts at the top of the existing mesh
//  (\var{TopLevel}) and works down.  At each level a set of mesh boxes that
//  includes all the tagged cells at that level is generated (using
//  \func{makeBoxes}, which is provided in the derived class (e.g.
//  BRMeshRefine).  The goal is to minimize the number of boxes
//  and maximize the ratio of tagged to untagged cells in each box.
//  The boxes are refined to the next higher level and saved as the new
//  mesh (in \var{newmeshes}).  The unrefined boxes are coarsened and
//  projected down to the next lower level where a buffer of bufferSize cells
//  is added around all the boxes (the opposite of taking the interior
//  cells that is done in building the PND) and added to the tagged cells
//  at the coarser level.  This ensures that the next coarser level will
//  have meshes that will properly nest the meshes just generated.  This
//  is repeated until the BaseLevel+1 mesh is replaced.
//
//
// Implementation Method:
// ----------------------
//  There are 5 basic steps to perform.  Two functions (\func{makePNDs} and
//  \func{makeBoxes}) are called to handle the most complex steps and the code
//  in this function handles the rest.
//  Start
//   1) Compute the proper nesting domains (PNDs).  A \var{Vector<IntVectSet>}
//      is created and passed to \func{makePNDs}, along with the grid and
//      boundary data.  Only the part of the Vector that is active
//      [\var{BaseLevel:TopLevel}] is modified.  The other elements of the
//      vector are not touched.
//   2) Removing (clipping) the tag cells outside the PNDs is trivial:
//      intersect the two \var{IntVectSets} (PND and tag cells) at each level.
//      [Note: The PND data is not needed after this step.]
//  Loop for L from \var{TopLevel} downto \var{BaseLevel}
//   3) \func{makeBoxes} is called to find a set of boxes that cover the tag
//      cells at this level.  The output is stored in a temp \var{BoxArray}
//   4) The temp \var{BoxArray} is refined and the result is stored in the
//      output mesh vector at level L+1.
//   5) The temp \var{BoxArray} is coarsened to the next coarser level, then
//      expanded by the buffer cells in each direction (clipping at the
//      physical boundaries) and then  unioned with the tags at level L-1.
//      This step does not have to be done for the base level.  The resulting
//      tags are used to generate the mesh at the next coarser level.
//  EndLoop
//
// Implementation Notes:
// ---------------------
//  The output vector \var{newmeshes} is redefined to the necessary size
//  without necessarily reusing the existing memory space.
//
///////////////////////////////////////////////////////////////////////////////

MeshRefine::MeshRefine() : m_isDefined(false), m_granularity(1)
{
}

MeshRefine::~MeshRefine()
{
}

MeshRefine::MeshRefine(const Box& a_baseDomain,
                       const Vector<int>& a_refRatios,
                       const Real a_fillRatio,
                       const int a_blockFactor,
                       const int a_bufferSize,
                       const int a_maxBoxSize)
  :m_granularity(1)
{
  ProblemDomain crseDom(a_baseDomain);
  define(crseDom, a_refRatios, a_fillRatio, a_blockFactor,
         a_bufferSize, a_maxBoxSize);
}

MeshRefine::MeshRefine(const ProblemDomain& a_baseDomain,
                       const Vector<int>& a_refRatios,
                       const Real a_fillRatio,
                       const int a_blockFactor,
                       const int a_bufferSize,
                       const int a_maxBoxSize)
  :m_granularity(1)
{
  define(a_baseDomain, a_refRatios, a_fillRatio, a_blockFactor,
         a_bufferSize, a_maxBoxSize);
}

void
MeshRefine::define(const Box& a_baseDomain,
                   const Vector<int>& a_refRatios,
                   const Real a_fillRatio,
                   const int a_blockFactor,
                   const int a_bufferSize,
                   const int a_maxBoxSize)
{
  ProblemDomain crseDom(a_baseDomain);
  define(crseDom, a_refRatios, a_fillRatio, a_blockFactor,
         a_bufferSize, a_maxBoxSize);
}

void
MeshRefine::define(const ProblemDomain& a_baseDomain,
                   const Vector<int>& a_refRatios,
                   const Real a_fillRatio,
                   const int a_blockFactor,
                   const int a_bufferSize,
                   const int a_maxBoxSize)
{
  // Change default PND mode from 0 to 1 (ndk 8.4.2008)
  //m_PNDMode = 0; // old, tested version of using nestingRegion
  m_PNDMode = 1; // New version which does not use expensive IntVectSet nestingRegion
  int maxLevel = a_refRatios.size();
  m_nRefVect.resize(maxLevel);
  m_vectDomains.resize(maxLevel+1);
  m_pnds.resize(maxLevel);
  //m_lastBase = maxLevel+1;
  //m_lastTop  = 0;
  //m_lastBuffer = 0;
  m_level_blockfactors.resize(maxLevel);

  m_vectDomains[0] = a_baseDomain;
  for (int lev=0; lev<maxLevel; lev++)
  {
    m_nRefVect[lev] = a_refRatios[lev];
    m_vectDomains[lev+1] = refine(m_vectDomains[lev], a_refRatios[lev]);
  }

  // do some quick sanity checks
  CH_assert( a_blockFactor >= 1 );
  CH_assert( a_bufferSize >= 0 );
  CH_assert( (a_maxBoxSize >= 2*a_bufferSize) || (a_maxBoxSize == 0) );
  CH_assert( (a_blockFactor <= a_maxBoxSize) || (a_maxBoxSize == 0) );
  CH_assert( a_fillRatio > 0.0 || a_fillRatio <= 1.0 );

  m_fillRatio = a_fillRatio;
  m_blockFactor = a_blockFactor;
  m_bufferSize = a_bufferSize;
  m_maxSize = a_maxBoxSize;

  computeLocalBlockFactors();

  m_isDefined = true;

}

// this function is overloaded by MultiBlockMeshRefine
void MeshRefine::clipBox(Box& a_box, const ProblemDomain& a_domain) const
{
  a_box&= a_domain;
}

/// sets proper nesting region granularity.
void MeshRefine::granularity(int a_granularity)
{
  CH_assert(a_granularity>0);
  m_granularity = a_granularity;
}

void MeshRefine::setPNDMode(int a_mode)
{
  CH_assert(a_mode == 0 || a_mode == 1);
  m_PNDMode = a_mode;
}

//
// This regrid function takes a single IntVectSet of tags
//
int
MeshRefine::regrid(Vector<Vector<Box> >&   a_newmeshes,
                   const IntVectSet&             a_tags,
                   const int               a_baseLevel,
                   const int               a_topLevel,
                   const Vector<Vector<Box> >& a_oldMeshes)
{
  //
  // Handle special cases
  //
  CH_assert (a_topLevel >= 0 );
  CH_assert ( a_baseLevel < (a_topLevel+1) && a_baseLevel >= 0 );

  //
  // Create a vector of tags from the single set of tags that is given
  Vector<IntVectSet> tags_vector( a_topLevel+1 ,a_tags );

  // Refine the levels above the base level
  for ( int i = 1 ; i <= a_topLevel ; i++ )
    {
      tags_vector[i] = refine(tags_vector[i-1], m_nRefVect[i-1] ) ;
    }

  // run the version of refine that takes a vector of tags
  return regrid(a_newmeshes, tags_vector,
                a_baseLevel, a_topLevel, a_oldMeshes);
}

//
// This version of refine takes a Vector of tags
//
int
MeshRefine::regrid(Vector<Vector<Box> >&   a_newmeshes,
                   Vector<IntVectSet>&     a_tags,
                   const int               a_baseLevel,
                   const int               a_topLevel,
                   const Vector<Vector<Box> >& a_OldMeshes)
{
  CH_TIME("MeshRefine::regrid");
  //
  // Validate arguments and handle special cases
  //
  CH_assert( isDefined());
  CH_assert( a_topLevel >= 0 );
  CH_assert( a_baseLevel < (a_topLevel+1) && a_baseLevel >= 0 );
  CH_assert( a_OldMeshes.size() >= a_topLevel + 1 );
  CH_assert( m_vectDomains.size() >= a_topLevel + 1 );
  CH_assert( m_nRefVect.size() >= a_topLevel + 1 );
  CH_assert( a_tags.size() >= a_topLevel+1 );

  Vector<IntVectSet>& modifiedTags =a_tags;

  // set the top level to be the finest level which actually has tags
  int TopLevel = a_topLevel;
  int new_finest_level;
  int isize = a_tags.size();
  int level;
  for (level = Min(a_topLevel, isize-1); level >= a_baseLevel; --level)
    {
      if (!a_tags[level].isEmpty()) break;
    }

#ifdef CH_MPI
  int mlevel;
  MPI_Allreduce(&level, &mlevel, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
  level=mlevel;
#endif

  if (level >= a_baseLevel)
    {
      TopLevel = level;
      // pout() << "  regrid TopLevel= " << TopLevel << endl;
      // copy tags to local storage so we can mess with them.

      //modifiedTags.resize(a_tags.size());
      Vector<int>  totalBufferSize(a_tags.size(), 0);
      //for (int ilev = a_baseLevel; ilev <= TopLevel; ilev++)
      //  {
      //    modifiedTags[ilev] = a_tags[ilev];
      //  }

      // reinitialize the output array .... (noel: why +2 ?)
      a_newmeshes.resize( TopLevel+2 );

      //
      // Generate new meshes if requested.
      if ( TopLevel+1 > a_baseLevel )
        {

          Box domaint = m_vectDomains[a_baseLevel].domainBox();
          Box testdom = coarsen(domaint, m_blockFactor);
          testdom.refine(m_blockFactor);
          if (domaint != testdom)
            MayDay::Error("MeshRefine:domain and Blocking Factor incompatible");
          // [NOTE: the following validations are expensive so only do them when
          //        debugging.]

#if ! defined( NDEBUG )
          // all existing meshes should be within the boundaries
          for ( int i = a_baseLevel ; i <= TopLevel ; i++ )
          {
            if (a_OldMeshes[i].size() > 0)
            {
              Box minbox = a_OldMeshes[i][0];
              for (int ibox = 0; ibox < a_OldMeshes[i].size(); ++ibox)
                {
                  minbox.minBox(a_OldMeshes[i][ibox]);
                }
              CH_assert(m_vectDomains[i].contains(minbox));
            }
          }
#endif
          // all tagged cells must be within the existing meshes
          // make sure that each level of \var{Domains} is consistent with the
          // previous level and the refinement ratio.
          for ( int i = a_baseLevel ; i < TopLevel ; i++ )
            {
              if (m_vectDomains[i+1].domainBox()
                 != refine(m_vectDomains[i].domainBox(), m_nRefVect[i]))
                MayDay::Error("MeshRefine:domains and refratios incompatible");
            }
          //
          // coarsen old meshes, tags, problem domains, and buffers by the
          // appropriate blocking factor
          Vector<Vector<Box> > OldMeshes = a_OldMeshes;
          Vector<ProblemDomain> Domains = m_vectDomains;
          Vector<int> blocked_BufferSize(TopLevel+1);
          //for (int level = 0; level < modifiedTags.size (); ++level)
          for (int level = a_baseLevel; level <= TopLevel; ++level)
            {
              // first intersect tags with domains (note that in periodic case, this can
              // result in wrapping of tags that are outside the base domain box)
              // int nBefore = modifiedTags[level].numPts();
              // IntVectSet tmp(modifiedTags[level]);
              //tmp &= m_vectDomains[level].domainBox();
              //int nTmp = tmp.numPts();
              modifiedTags[level] &= m_vectDomains[level];
              //int nAfter = modifiedTags[level].numPts();
              //pout()<<"nBefore:"<<nBefore<<" nAfter:"<<nAfter<<" nTmp:"<<nTmp
              //            <<" diff:"<<nBefore-nAfter<<" newDiff:"<<nTmp-nAfter <<"\n";
              // coarsen by blocking factor, rounding upwards the goal of this is to
              // coarsen everything down to a level which is m_blockFactor coarser than
              // the new fine level we're going to generate.  By generating grids at
              // this fake level, we can then refine up to the new level, and the blocking
              // factor will be automatically enforced.
              modifiedTags[level].coarsen (m_level_blockfactors[level]);
              // Same block-factor coarsening for the domains
              Domains[level].coarsen(m_level_blockfactors[level]);
              // Same block-factor coarsening for the buffers
              blocked_BufferSize[level] =
                (m_bufferSize + m_level_blockfactors[level] - 1)/
                m_level_blockfactors[level];
            }
          // We only need the base mesh coarsened by the blocking factor
          Vector<Box> OldBaseMesh = a_OldMeshes[a_baseLevel];
          {
            const int crFactor = m_level_blockfactors[a_baseLevel];
            for (int i = 0; i != OldBaseMesh.size(); ++i)
              {
                OldBaseMesh[i].coarsen(crFactor);
              }
          }

          for (int i=a_baseLevel; i<=TopLevel; i++)
            {
              m_pnds[i].makeEmpty();
            }
          makePNDs(m_pnds, totalBufferSize, a_baseLevel, TopLevel, Domains,  OldBaseMesh,
                   blocked_BufferSize);

          // Clip the tagged cells to the proper nesting domains by
          // intersecting the two sets.
          // Note: if m_PNDMode = 1, this only ensures the tags are in the base
          // PND!  Remaining adherence to each level's PND is performed in the
          // box generation algorithm
          for ( int lvl = a_baseLevel ; lvl <= TopLevel ; lvl++ )
          {
            // int nBefore = modifiedTags[lvl].numPts();
            modifiedTags[lvl] &= m_pnds[lvl] ;
            //int nAfter = modifiedTags[lvl].numPts();
            //pout()<<"nBefore:"<<nBefore<<" nAfter:"<<nAfter<<" diff:"<<nBefore-nAfter<<"\n";
          }

          //
          // Generate new meshes.
          //

          // At each level, starting at the top, generate boxes that cover the tag
          // cells using makeBoxes(), refine these boxes and save them in the
          // output \var{a_newmeshes}.  Take the unrefined boxes, add a buffer zone
          // around each box, coarsen, and union with the tag cells on the next
          // coarser level.  This modifies the tags variable.  To handle
          // \var{BlockFactor}, coarsen everything before making the new
          // meshes and then refine the resulting mesh boxes.
          Vector<Box> lvlboxes ;  // new boxes on this level
          for ( int lvl = TopLevel ; lvl >= a_baseLevel ; lvl-- )
          {
            // make a new mesh at the same level as the tags

            const int dest_proc = uniqueProc(SerialTask::compute);

            Vector<IntVectSet> all_tags;
            gather(all_tags, modifiedTags[lvl], dest_proc);

            if (procID() == dest_proc)
              {
                for (int i = 0; i < all_tags.size(); ++i)
                  {
//                     modifiedTags[lvl] |= all_tags[i];
                    //**FIXME -- revert to above line when IVS is fixed.
                    //**The following works around a bug in IVS that appears if
                    //**the above line is used.  This bug is observed when there
                    //**is a coarsening of an IVS containing only IntVect::Zero
                    //**followed by an IVS |= IVS.
                    for (IVSIterator ivsit(all_tags[i]); ivsit.ok(); ++ivsit)
                      {
                        modifiedTags[lvl] |= ivsit();
                      }
                    //**FIXME -- end
                    // Regain memory used (BVS,NDK 6/30/2008)
                    all_tags[i].makeEmpty();
                  }
              }

            broadcast( modifiedTags[lvl] , dest_proc);

            // Move this union _after_ the above gather/broadcast to
            // reduce memory -- shouldn't have other effects. (BVS,NDK 6/30/2008)
            // Union the meshes from the previous level with the tags on this
            // level to guarantee that proper nesting is satisfied.  On the
            // first iteration this is a no-op because \var{lvlboxes} is empty.
            // [NOTE: for every iteration after the first, \var{lvlboxes} will
            //        already be coarsened by \var{BlockFactor} so it will be
            //        at the same refinement level as \var{tags[lvl]}, which
            //        has also been coarsened]
            // this is simple in the non-periodic case, more complicated
            // in the periodic case
            ProblemDomain lvldomain = Domains[lvl]; // domain of this level
            buildSupport(lvldomain, lvlboxes, modifiedTags[lvl]);

            // this is the maximum allowable box size at this resolution
            // which will result in satisfying the maxSize restriction when
            // everything is refined up to the new level
            const int maxBoxSizeLevel = m_maxSize/(m_level_blockfactors[lvl]*m_nRefVect[lvl]);
            makeBoxes(lvlboxes, modifiedTags[lvl], m_pnds[lvl],
                      lvldomain, maxBoxSizeLevel, totalBufferSize[lvl]);
            // After change to reduce memory, this may now be needed.
            // Previously, there were a_tags.makeEmpty() calls in BRMesh.cpp, and now,
            // if there are a few tags leftover here, they will get added onto the mix -- which is not
            // the behavior of the code before these changes (ndk) (leaving commented for now)
            modifiedTags[lvl].makeEmpty(); // <-- why is this necessary now?

            // refine the new mesh and save it
            //[NOTE: we have to undo the coarsening by \var{BlockFactor} as well
            //       as refine to the next level.]
            //[NOTE: refine() operates in-place so copy first then refine()
            //       because the unrefined mesh will be needed later.]
            a_newmeshes[lvl+1] = lvlboxes ;
            for (int ibox = 0; ibox < a_newmeshes[lvl+1].size(); ++ibox)
              {
                a_newmeshes[lvl+1][ibox].refine(m_level_blockfactors[lvl]*m_nRefVect[lvl]);
              }
            // Make the boxes ready for the next iteration.
            // Don't have to do this for the last iteration.
            if ( lvl > a_baseLevel )
              {
                // coarsen the unrefined new mesh so it matches the tags
                // at that level, then add the buffer cells, clip at the domain
                // boundaries so we can union them in the next iteration
                // allInOne_nRef does:
                //   a) refine  by BF_ref[lvl];
                //   b) coarsen by  n_ref[lvl-1];
                //   c) coarsen by BF_ref[lvl-1];
                const int allInOne_nRef =
                  (m_level_blockfactors[lvl-1]*m_nRefVect[lvl-1])/
                  m_level_blockfactors[lvl];
                for (int ibox = 0 ; ibox < lvlboxes.size() ; ++ibox)
                  {         
                    lvlboxes[ibox].grow(blocked_BufferSize[lvl]);
                    clipBox(lvlboxes[ibox], Domains[lvl]);
                    lvlboxes[ibox].coarsen(allInOne_nRef);
                  }
              }
          }
        }

      //
      // Finally, copy the old mesh levels that didn't change.
      //
      for ( int i = 0 ; i <= a_baseLevel ; i++ )
      {
        a_newmeshes[i] = a_OldMeshes[i] ;
      }

      //
      // Done generating grids
      //

      // set new finest level
      new_finest_level = TopLevel+1;
      // this is designed to catch the pathological but possible case
      // where there were tags on the TopLevel, but no grids were generated
      // (possibly if all tags were outside pnds)
      while (a_newmeshes[new_finest_level].size() == 0
             && new_finest_level > a_baseLevel)
        {
          new_finest_level -= 1;
        }

    }
  else
    {
      // if no tags on any level, just return
      new_finest_level = a_baseLevel;
      a_newmeshes.resize(a_OldMeshes.size());
      for ( int i = 0 ; i <= a_baseLevel ; i++ )
      {
        a_newmeshes[i] = a_OldMeshes[i] ;
      }
    }

  return new_finest_level;
}

void MeshRefine::buildSupport(const ProblemDomain& lvldomain, Vector<Box>& lvlboxes, IntVectSet& modifiedTags)
{
  if (lvldomain.isPeriodic())
    {
      const Box domainBox = lvldomain.domainBox();
      ShiftIterator shiftIt = lvldomain.shiftIterator();
      IntVect shiftMult(domainBox.size());
      for (int i=0; i<lvlboxes.size(); i++)
        {
          Box localBox(lvlboxes[i]);
          // will handle periodic wraparound through shifting and
          // adding shifted image to tags, which will all remain
          // within the domainBox
          localBox &= domainBox;
          modifiedTags |= localBox;
          // now do shifts to capture periodic images necessary to
          // enforce proper nesting only do this if original box was
          // not contained by domainBox
          if (localBox != lvlboxes[i])
            {
              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  IntVect shiftVect(shiftIt()*shiftMult);
                  Box localShiftedBox(lvlboxes[i]);
                  localShiftedBox.shift(shiftVect);
                  localShiftedBox &= domainBox;
                  if (!localShiftedBox.isEmpty())
                    {
                      modifiedTags |= localShiftedBox;
                    }
                } // end loop over shift directions
            } // end whether periodic checking was needed
        } // end loop over finer-level boxes to enforce nesting
    }
  else
    {
      // non periodic case is simple
      for ( int i = 0 ; i < lvlboxes.size() ; i++ )
        {
          modifiedTags |= lvlboxes[i] ;
        }
    }
}

inline bool isPower2(const int a_i)
{
  if (a_i <= 0) return false;      // Catch <= 0
  unsigned i = a_i;
  while (!(i & 1))                 // Find first bit
    {
      i >>= 1;
    }
  return (i >> 1) ? false : true;  // Can't have any other bits
}

void
MeshRefine::computeLocalBlockFactors()
{
  // computes amount that tags on a level l must be coarsened so
  // that grids on next finer level l+1 will be guaranteed to
  // satisfy the blocking factor.
  // The easiest way to do this is to coarsen everything down
  // to a level which is m_blockFactor coarser than the new
  // fine level which will be generated.  So, m_level_blockfactor[lvl]
  // is the amount by which the tags, etc at level lvl must
  // be coarsened so that they are m_blockFactor coarser than
  // the (lvl+1) grids which will be generated.  If m_BlockFactor
  // is less than the refinement ratio between levels (lvl) and (lvl+1),
  // then no coarsening needs to be done, so we default to one, in that case.
  for (int lev=0; lev<(m_level_blockfactors.size()); lev++)
    {
      // This is simply ceil(m_blockFactor / m_nRefVect[lev]).
      m_level_blockfactors[lev] = (m_blockFactor + m_nRefVect[lev] -1)/m_nRefVect[lev];
      // The coarsening ratio due to blocking needs to be a power of 2 otherwise
      // IntVectSet will complain.  Catch here to better describe the error.
      if (!isPower2(m_level_blockfactors[lev]))
        {
          pout() << "Unable to implement blocking for level " << lev+1 << ".  "
            "Blocking requires ceil(blockFactor/nRef) to be a power of 2 but "
            "for nRef[" << lev << "] = " << m_nRefVect[lev]
                 << " and blockFactor = " << m_blockFactor << ", this is "
                 << m_level_blockfactors[lev] << '.' << endl;
          MayDay::Error("aborting MeshRefine::regrid");
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
//
// Function:
// ---------
//  int makePNDs()
//
// Short Description:
// ------------------
//  Compute the proper nesting domains for a subset of the levels in a
//  multilevel adaptive mesh.
//
// Usage:
// ------
//  This is called by \func{regrid}.  It wasn't designed to be called
//  directly so if you do so you're on your own.  The output is a
//  \var{Vector<IntVectSet>}.  This function is called to compute the proper
//  nesting domains needed to do mesh refinement.  The function that calls this
//  one must have the mesh and problem domain data structures already defined
//  for all the levels and have allocated a vector<IntVectSet> for the proper
//  nesting domains.  This function computes the proper nesting domains of the
//  levels to be refined (i.e., the levels \var{a_baseLevel} and above) and
//  modifies the \var{pnds} array.  The output vector is undefined for
//  indices corresponding to mesh levels below \var{a_baseLevel}
//  (i.e. [0:\var{a_baseLevel}-1]).
//
//  NOTE: in MeshRefine::regrid, this function is called _after_ everthing
//        has been coarsened to enforce the blocking factor, so this behaves
//        differently than one might expect -- call outside of the regrid fn
//        at your own risk! (DFM 8/24/01)
//
// Arguments:
// ----------
//  \var{pnds}                     type: Vector<IntVectSet>
//      (output, index[\var{a_baseLevel}:\var{TopLevel}]) proper nesting domains
//
//  \var{a_baseLevel}                type: int
//      (input, range[0:\var{TopLevel}]) index of the mesh level to use as the
//      source for the proper nesting domains
//
//  \var{TopLevel}                 type: int
//      (input, range[\var{a_baseLevel}:\var{OldMeshes.size()-1]) index of the
//      finest mesh level to compute a PND for
//
//  \var{BaseMesh}                 type: Vector<Box>
//      (input) the boxes at level \var{a_baseLevel} in the mesh hierarchy.
//
//  \var{Domains}                  type: Vector<Box>
//      [same as for \func{MeshRefine}.]
//
//  \var{BufferSize}              type: int
//      [same as for \func{MeshRefine}.]
//
// Returns:
// --------
//  An int giving an exception number is returned.  The argument
//  \var{pnds} is modified.
//
// References:
// -----------
//  See ??? for a description of Proper Nesting Domains.
//  meshRefine -- calls this function
//
// Numerical Algorithm:
// --------------------
//  The proper nesting domain (PND) of level L+1 is defined as the interior of
//  the PND for level L (the next coarser) except at the boundaries of the
//  problem domain (ie, the level 0 box), in which case the PNDs of levels L
//  and L+1 both contain the boundary.  The PND of the base level is defined as
//  the interior of the union of all the grid boxes in the base level.  Given
//  the PND of level L (denoted by PND[L]) this function computes PND[L+1] by
//  iterating over the coordinate directions and computing the cells to remove
//  from the proper nesting domain due to the boundaries in that direction.
//  For each direction the domain is shifted (first up, then down) and only the
//  part that is shifted _out_ of the PND but remains _inside_ the problem
//  boundary is kept and then shifted back into the PND and subtracted
//  (removed) from it.  So if the shift moves part of the PND outside the
//  problem boundary then no cells will be removed.  If the cells that get
//  shifted out of the PND remain inside the boundary then they will be removed
//  from the PND when they get shifted back.  When all the directions have been
//  handled, the PND is refined to make it correct for level L+1.  The PND of
//  the base level is the union of the mesh boxes in the base level and serves
//  as the starting point.  In set notation the operation is:
//
//            D - ( ( ( (D << d) - D ) * B ) >> d )
//
//  where:
//    d is the current direction (+i,-i,+j,etc)
//    B is the boundary of the problem domain at the current level
//    - is the subtraction (removal) operator
//    * is the intersection operator
//   << is the down-shift operator (shift in negative direction)
//   >> is the up-shift operator (shift in positive direction)
//
//
// Implementation Method:
//  The PNDs are implemented as a vector of \type{IntVectSet}s.  A scratch
//  IntVectSet is used to store the PND of the current level while it is
//  being computed.  When the PND for a level is complete, it is refined
//  and used for the next level.  The base PND is computed as the union of
//  the mesh boxes in the base level mesh.  The PNDs for each level are
//  stored in the appropriate elements of the pnds vector as they are
//  computed.  Levels below a_baseLevel are not modified.  The code loops
//  through each level starting from Baselevel.  Loop through each coordinate
//  direction applying the operation defined in the "Numerical Algorithm"
//  section above.  Copy the final result into the PNDs vector.  If an
//  exception occurs, do not change PNDs on the level that took the
//   exception.
//
// Implementation Notes:
//  This assumes cell-centers are being manipulated, not vertices.  The whole
//  vector of pnds is passed even though only some are accessed
//  because it is more convenient since these variables are likely to exist
//  already.
//
///////////////////////////////////////////////////////////////////////////////

void
MeshRefine::makePNDs(Vector<IntVectSet>&          a_pnds,
                     Vector<int>&                 a_totalBufferSize,
                     const int                    a_baseLevel,
                     const int                    a_topLevel,
                     const Vector<ProblemDomain>& a_domains,
                     const IntVectSet&            a_baseMesh,
                     const Vector<int>&           a_bufferSize) const
{
  //
  // Validate inputs
  //
  CH_assert( a_baseLevel <= a_topLevel && a_baseLevel >= 0 );
  CH_assert(   a_domains.size() >= a_topLevel + 1 );
  CH_assert( m_nRefVect.size() >= a_topLevel + 1 );
  // all existing boxes on the base level must be aligned to \var{BuffSize}

  a_pnds[a_baseLevel] = a_baseMesh;
  a_totalBufferSize[a_baseLevel] = 0;

  for ( int lvl = a_baseLevel ; lvl <= a_topLevel ; lvl++)
    {
      IntVectSet& pnd = a_pnds[lvl];

      if (m_PNDMode == 0)
      {
        pnd.nestingRegion(a_bufferSize[lvl], a_domains[lvl], m_granularity);
      }
      a_totalBufferSize[lvl] += a_bufferSize[lvl];

      if ( (a_topLevel - lvl) > 0)
        {
          a_pnds[lvl+1] = pnd;

          // This does:
          //   a) refine  by BF_ref[lvl];
          //   b) refine  by  n_ref[lvl];
          //   c) coarsen by BF_ref[lvl+1];
          const int allInOne_nRef = (m_level_blockfactors[lvl]*m_nRefVect[lvl])/
            m_level_blockfactors[lvl+1];

          a_totalBufferSize[lvl+1] = a_totalBufferSize[lvl]*allInOne_nRef;
          a_pnds[lvl+1].refine(allInOne_nRef);
        }
    }
}

void
MeshRefine::makePNDs(Vector<IntVectSet>&  a_pnds,
                     Vector<int>&         a_totalBufferSize,
                     const int            a_baseLevel,
                     const int            a_topLevel,
                     const Vector<ProblemDomain>&     a_domains,
                     const Vector<Box>&   a_oldMeshes,
                     const Vector<int>&   a_bufferSize) const
{

  //pout()<<"["<<a_baseLevel<<","<<a_topLevel<<"]";

  IntVectSet mesh;

  //const ProblemDomain& domain = a_domains[a_baseLevel];
  for (int box = 0; box < a_oldMeshes.size(); ++box)
    {
      const Box& b = a_oldMeshes[box];
      mesh |= b;
    }
  makePNDs(a_pnds, a_totalBufferSize, a_baseLevel, a_topLevel, a_domains, mesh, a_bufferSize);

  return;
}

//end of makePNDs

bool MeshRefine::properlyNested(const Box& a_box,
                                const ProblemDomain& a_domain,
                                const IntVectSet& a_pnd,
                                int   a_totalBufferSize) const
{
  if (m_PNDMode == 0)
    {
      return a_pnd.contains(a_box);
    }
  else
    {
      Box growBox(a_box);
      growBox.grow(a_totalBufferSize);
      growBox &= a_domain.domainBox();
      if (!a_pnd.contains(growBox)) return false; //typical case
      if (a_domain.isPeriodic())
      {
        Box growPeriodic(a_box);
        growPeriodic.grow(a_totalBufferSize);
        growPeriodic &= a_domain; // now intersect with (possibly) periodic domain
        if (growPeriodic != growBox)
        {
          //has periodic images
          //  for (BoxIterator bit(growPeriodic); bit.ok(); ++bit)
          //        {
          //          if (!growBox.contains(bit())){
          //            IntVect image=bit();
          //            bool isImage = a_domain.image(image);
          //            CH_assert(isImage);
          //            if (!a_pnd.contains(image)) return false;
          //          }
          //        }
          static ImageIterator images(a_domain);
          images.checkDefine(a_domain);
          for (images.begin(growPeriodic); images.ok(); ++images)
          {
            if (!images.box().isEmpty())
            {
              if (!a_pnd.contains(images.box()))
              {
                return false;
              }
            }
          }
        }
      }

    }
  return true;
}

// ---------------------------------------------------------------
// Utility and access functions
// --------------------------------------------------------------

/// returns vector of refinement ratios
const Vector<int>&
MeshRefine::refRatios() const
{
  return m_nRefVect;
}

/// returns fillRatio
Real
MeshRefine::fillRatio() const
{
  return m_fillRatio;
}

/// returns blocking factor
int
MeshRefine::blockFactor() const
{
  return m_blockFactor;
}

/// returns proper nesting buffer size
int
MeshRefine::bufferSize() const
{
  return m_bufferSize;
}

/// returns maximum box size in any dimension
int
MeshRefine::maxSize() const
{
  return m_maxSize;
}

/// sets vector of refinement ratios
void
MeshRefine::refRatios(const Vector<int>& a_nRefVect)
{
  int isize = a_nRefVect.size();
  m_nRefVect.resize(isize);
  for (int i=0; i<isize; i++)
    {
      m_nRefVect[i] = a_nRefVect[i];
    }
  // will need to recompute local blockfactors in this case
  computeLocalBlockFactors();
}

/// sets fillRatio
void
MeshRefine::fillRatio(const Real a_fill_ratio)
{
  m_fillRatio = a_fill_ratio;
}

/// sets blocking factor
void
MeshRefine::blockFactor(const int a_block_factor)
{
  m_blockFactor = a_block_factor;
  // then need to recompute local blockFactors
  computeLocalBlockFactors();
}

/// sets proper nesting buffer size
void
MeshRefine::bufferSize(const int a_buffer_size)
{
  m_bufferSize = a_buffer_size;
}

/// sets maximum box size in any dimension
void
MeshRefine::maxSize(const int a_max_size)
{
  m_maxSize = a_max_size;
}

/// has this object been defined properly?
bool
MeshRefine::isDefined() const
{
  return m_isDefined;
}
#include "NamespaceFooter.H"
