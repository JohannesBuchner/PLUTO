#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _STLBOX_H_
#define _STLBOX_H_

#include "ProblemDomain.H"
#include "RealVect.H"
#include "IntVect.H"
#include "Box.H"
#include "RefCountedPtr.H"

#include "STLMesh.H"
#include "STLUtil.H"

#include "NamespaceHeader.H"

using namespace STLUtil;

/*
 * This is just a class to hold data that associates
 * a Chombo Box with an STL mesh.  It has, e.g., maps
 * between Chombo cells and the STL vertices that lie
 * within them.  The class also includes some get
 * functions and print functions for its members
 */

class STLBox
{
public:

  STLBox(RefCountedPtr<STLMesh> a_stlmesh,
         const Box&             a_region,
         const ProblemDomain&   a_domain,
         const RealVect&        a_origin,
         const RealVect&        a_dx);

  void SetMeshBox(RefCountedPtr<STLMesh> a_stlmesh,
                  const Box&             a_region,
                  const ProblemDomain&   a_domain,
                  const RealVect&        a_origin,
                  const RealVect&        a_dx);

  RefCountedPtr<STLMesh> m_msh; // pointer to mesh

  // box that we're exploring
  Box           m_region;
  ProblemDomain m_domain;
  RealVect      m_origin;
  RealVect      m_dx;

  // for each cell, a Vector of the vertices and a Vector of the triangles
  // that wholly or partially are contained by the cell
  // vertices and triangles are specified by index (see STLStruct.H for TriInCell)
  CellMap m_cellmap;

  // for each Chombo node, whether it is 'inside' or 'outside'
  // according to the normal vectors specified in the STL file and mesh.triangles.normal
  NodeMap m_nodemap;

  // for each Chombo edge that intersects the boundary, store the
  // intersection point (since we already compute it, this just makes
  // it faster to retrieve for queries)
  EdgeMap m_edgemap;

  // get maps
  void GetCellMap(CellMap** a_cellmap);
  void GetNodeMap(NodeMap** a_nodemap);
  void GetEdgeMap(EdgeMap** a_edgemap);
};

#include "NamespaceFooter.H"
#endif

