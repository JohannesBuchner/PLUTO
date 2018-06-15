#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL, DTG

#include "EBISBox.H"
#include "VoFIterator.H"
#include "BoxIterator.H"
#include "parstream.H"
#include "SPMD.H"
#include "limits.h"
#include "PolyGeom.H"
#include "NamespaceHeader.H"
/*******************************/
const Box&
EBISBox::getRegion() const
{
  return m_graph.getRegion();
}
/*******************************/
EBISBox::EBISBox()
{
}
/*******************************/
EBISBox::~EBISBox()
{
}
/*******************************/
IntVectSet
EBISBox::getMultiCells(const Box& a_subbox) const
{
  return m_graph.getMultiCells(a_subbox);
}
/*******************************/
IntVectSet
EBISBox::getIrregIVS(const Box& a_subbox) const
{
  return m_graph.getIrregCells(a_subbox);
}

///
/**
   Returns the irregular cells that have non-zero boundary area
*/

IntVectSet EBISBox::boundaryIVS(const Box& a_subbox) const
{
  IntVectSet ivs = getIrregIVS(a_subbox);
  IntVectSet rtn = ivs;
  CH_assert(rtn.isDense());
  IVSIterator it(ivs);
  for (;it.ok(); ++it)
  {
    if (m_data.bndryArea(VolIndex(it(),0)) == 0) rtn -= it();
  }
  return rtn;
}

/*******************************/
const ProblemDomain&
EBISBox::getDomain() const
{
  return m_graph.getDomain();
}
/*******************************/

/*******************************/
Vector<VolIndex>
EBISBox::getVoFs(const VolIndex& a_vof,
                 const int& a_dir,
                 const Side::LoHiSide& a_sd,
                 const int& a_steps) const
{
  CH_assert((a_dir >= 0) && (a_dir < SpaceDim));
  CH_assert(a_steps >= 0);

  Vector<VolIndex> retVoFs(1, a_vof);
  for (int irad = 1; irad <= a_steps; irad++)
    {
      Vector<VolIndex> tempVoFs(0);
      for (int ivof = 0; ivof < retVoFs.size(); ivof++)
        {
          const VolIndex& stepVoF = retVoFs[ivof];
          Vector<FaceIndex> faces = getFaces(stepVoF, a_dir, a_sd);
          for (int iface = 0; iface < faces.size(); iface++)
            {
              const FaceIndex& face = faces[iface];
              if (!face.isBoundary())
                {
                  const VolIndex& flipVoF = face.getVoF(a_sd);
                  tempVoFs.push_back(flipVoF);
                }
            }
        }
      retVoFs = tempVoFs;
    }
  return retVoFs;
}
/*******************************/

/*******************************/
int
EBISBox::numFaces(const VolIndex& a_vof,
                  const int& a_idir,
                  const Side::LoHiSide& a_sd) const
{
  Vector<FaceIndex> faces = getFaces(a_vof, a_idir, a_sd);
  int retval = faces.size();
  return retval;
}
/*******************************/
Real
EBISBox::volFrac(const VolIndex& a_vof) const
{
  Real retval;
  if (isRegular(a_vof.gridIndex()))
    {
      retval = 1.0;
    }
  else if (isCovered(a_vof.gridIndex()))
    {
      retval = 0.0;
    }
  else
    {
      retval = m_data.volFrac(a_vof);
    }

  return retval;
}
/*******************************/
Real
EBISBox::areaFracScaling(const VolIndex& a_vof) const
{
  Real alphaMax = 0;
  for (int idir=0; idir<SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          Vector<FaceIndex> faces = getFaces(a_vof, idir, sit());
          for (int iface=0; iface<faces.size(); iface++)
            {
              alphaMax = Max(areaFrac(faces[iface]), alphaMax);
            }
        }
    }
  if (alphaMax > 0)
    {
      return (1./alphaMax);
    }
  else
    {
      return 1.;
    }
}
/*******************************/
Real
EBISBox::sumArea(const VolIndex& a_vof,
                 const int& a_idir,
                 const Side::LoHiSide& a_sd) const
{
  Real retval;
  if (isRegular(a_vof.gridIndex()))
    {
      retval = 1.0;
    }
  else if (isCovered(a_vof.gridIndex()))
    {
      retval = 0.0;
    }
  else
    {
      retval = 0.0;
      Vector<FaceIndex> faces = getFaces(a_vof, a_idir, a_sd);
      for (int iface = 0; iface < faces.size(); iface++)
        {
          retval += areaFrac(faces[iface]);
        }
    }
  return retval;
}
/*******************************/
Vector<VolIndex>
EBISBox::refine(const VolIndex& a_coarVoF) const
{
  return(m_graph.refine(a_coarVoF));
}
/*******************************/
VolIndex
EBISBox::coarsen(const VolIndex& a_fineVoF) const
{
  return(m_graph.coarsen(a_fineVoF));
}
/*******************************/
bool
EBISBox::isConnected(const VolIndex& a_vof1,
                     const VolIndex& a_vof2) const
{
  return m_graph.isConnected(a_vof1, a_vof2);
}
/*******************************/
Real
EBISBox::areaFrac(const FaceIndex& a_face) const
{
  Real retval;

  Box region = m_graph.getRegion();
  const IntVect& loiv = a_face.gridIndex(Side::Lo);
  const IntVect& hiiv = a_face.gridIndex(Side::Hi);
  if (region.contains(loiv) && isRegular(loiv) )
    {
      retval = 1.0;
    }
  else if (region.contains(hiiv) && isRegular(hiiv))
    {
      retval = 1.0;
    }
  else if (region.contains(loiv) && isCovered(loiv))
    {
      retval =  0.0;
    }
  else if (region.contains(hiiv) && isCovered(hiiv))
    {
      retval =  0.0;
    }
  else
    {
      retval = m_data.areaFrac(a_face);
    }
  return retval;
}
/*******************************/
RealVect
EBISBox::normal(const VolIndex& a_vof) const
{
  //  Real bndryArea = PolyGeom::bndryArea(a_vof, *this);
  //  RealVect retval = PolyGeom::normal(a_vof, *this, bndryArea);

  RealVect retval;
  const IntVect& iv = a_vof.gridIndex();
  if (isRegular(iv))
    {
      retval = BASISREALV(0);
    }
  else if (isCovered(iv))
    {
      retval = BASISREALV(0);
    }
  else
    {
      retval = m_data.normal(a_vof);
    } //end else (vof is irregular)

  return retval;
}

RealVect
EBISBox::normal(const VolIndex& a_vof, int face) const
{
  return m_data.normal(a_vof, face);
}
/*******************************/
RealVect
EBISBox::centroid(const FaceIndex& a_face) const
{
  RealVect retval;

  const IntVect& loiv = a_face.gridIndex(Side::Lo);
  const IntVect& hiiv = a_face.gridIndex(Side::Hi);
  Box region = m_graph.getRegion();
  if (region.contains(loiv) && isRegular(loiv))
    {
      retval = RealVect::Zero;
    }
  else if (region.contains(loiv) && isCovered(loiv))
    {
      retval = RealVect::Unit;
    }
  else if (region.contains(hiiv) && isRegular(hiiv))
    {
      retval = RealVect::Zero;
    }
  else if (region.contains(hiiv) && isCovered(hiiv))
    {
      retval = RealVect::Unit;
    }
  else
    {
      retval = m_data.centroid(a_face);
    }
  return retval;

}
/*******************************/
RealVect
EBISBox::centroid(const VolIndex& a_vof) const
{
  RealVect retval;
  const IntVect& iv = a_vof.gridIndex();
  if (isRegular(iv) || (isCovered(iv)))
    {
      retval = RealVect::Zero;
    }
  else
    {
      retval = m_data.centroid(a_vof);
    }
  return retval;
}
/*******************************/
RealVect
EBISBox::bndryCentroid(const VolIndex& a_vof) const
{
  RealVect retval;
  const IntVect& iv = a_vof.gridIndex();
  if (isRegular(iv))
    {
      retval = RealVect::Unit;
      retval *= -1.0;
    }
  else if (isCovered(iv))
    {
      retval = RealVect::Unit;
      retval *= -1.0;
    }
  else
    {
      retval = m_data.bndryCentroid(a_vof);
    }
  return retval;
}

RealVect
EBISBox::bndryCentroid(const VolIndex& a_vof, int face) const
{

  return  m_data.bndryCentroid(a_vof, face);
}

/*******************************/
Real
EBISBox::bndryArea(const VolIndex& a_vof, int face) const
{
  return m_data.bndryArea(a_vof, face);
}

Real
EBISBox::bndryArea(const VolIndex& a_vof) const
{
  //  Real retval = PolyGeom::bndryArea(a_vof, *this);

  Real retval;
  const IntVect& iv = a_vof.gridIndex();
  if (isRegular(iv))
    {
      retval = 0.0;
    }
  else if (isCovered(iv))
    {
      retval = -1.0;
    }
  else
    {
      retval = m_data.bndryArea(a_vof);
    }
  return retval;
}

int
EBISBox::numFacePhase(const VolIndex& a_vof) const
{
  // I'm not going to pretent you can ask for a bndry face
  // information is the absence of a face.
  return m_data.numFacePhase(a_vof);
}

/// used by multi-fluid code
int
EBISBox::facePhase(const VolIndex& a_vof, int face) const
{
  return m_data.facePhase(a_vof, face);
}

/// used by multi-fluid code
const VolIndex&
EBISBox::faceIndex(const VolIndex& a_vof, int face) const
{
  return m_data.faceIndex(a_vof, face);
}

/*******************************/
void
EBISBox::define(const EBGraph&  a_graph,
                const EBData&   a_data)
{
  m_graph = a_graph;
  m_data  = a_data;
}
/*******************************/
void EBISBox::setDomain(const ProblemDomain& a_domain)
{
  m_graph.setDomain(a_domain);
}
/*******************************/
void EBISBox::setToAllRegular()
{
  m_graph.setToAllRegular();
}
/*******************************/
void EBISBox::setToAllCovered()
{
  m_graph.setToAllCovered();
}
/*******************************/
Vector<FaceIndex>
EBISBox::refine(const FaceIndex& a_coarFace,const EBISBox& a_fineEBISBox) const
{
  return m_graph.refine(a_coarFace, a_fineEBISBox.m_graph);
}
/*******************************/
FaceIndex
EBISBox::coarsen(const FaceIndex& a_fineFace) const
{
  return m_graph.coarsen(a_fineFace);
}
/*******************************/
EBISBox& EBISBox::operator=(const EBISBox& a_ebiin)
{
  if (&a_ebiin != this)
    {
      m_graph = a_ebiin.m_graph;
      m_data  = a_ebiin.m_data;
    }
  return *this;
}
/*******************************/
EBISBox::EBISBox(const EBISBox& a_ebiin)
{
  m_graph = a_ebiin.m_graph;
  m_data  = a_ebiin.m_data;
}
/*******************************/
const EBGraph& EBISBox::
getEBGraph() const
{
  return m_graph;
}
/*******************************/
const EBData& EBISBox::
getEBData() const
{
  return m_data;
}
/*******************************/
bool EBISBox::
operator==(const EBISBox& a_ebiin)
{
  return (m_graph == a_ebiin.m_graph) &&  (m_data == a_ebiin.m_data);
}
/*******************************/
#include "NamespaceFooter.H"
