#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "STLBox.H"

#include "NamespaceHeader.H"

STLBox::STLBox(RefCountedPtr<STLMesh> a_stlmesh,
               const Box&             a_region,
               const ProblemDomain&   a_domain,
               const RealVect&        a_origin,
               const RealVect&        a_dx)
{
  SetMeshBox(a_stlmesh,a_region,a_domain,a_origin,a_dx); // just set some stuff in other method
}

void STLBox::SetMeshBox(RefCountedPtr<STLMesh> a_stlmesh,
                        const Box&             a_region,
                        const ProblemDomain&   a_domain,
                        const RealVect&        a_origin,
                        const RealVect&        a_dx)
{
  // set data, used in construction or to reset data later (e.g. after default construction)
  m_msh     = a_stlmesh;
  m_region  = a_region;
  m_domain  = a_domain;
  m_origin  = a_origin;
  m_dx      = a_dx;
}

/*
 * Functions to return data structures
 */

void STLBox::GetCellMap(CellMap** a_cellmap)
{
  *a_cellmap = &m_cellmap;
}

void STLBox::GetNodeMap(NodeMap** a_nodemap)
{
  *a_nodemap = &m_nodemap;
}

void STLBox::GetEdgeMap(EdgeMap** a_edgemap)
{
  *a_edgemap = &m_edgemap;
}

#include "NamespaceFooter.H"
