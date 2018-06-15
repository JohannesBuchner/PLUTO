#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#include "IrregNode.H"
#include "NamespaceHeader.H"

/*******************************/
/*******************************/
IrregNode::IrregNode()
{
}
/*******************************/
/*******************************/
IrregNode::~IrregNode()
{
}


std::ostream& operator<< (std::ostream&  a_os,
                          const IrregNode& a_iv)
{
  a_os<<a_iv.m_cell<<" index:"<<a_iv.m_cellIndex<<" volFrac:"<<a_iv.m_volFrac
      <<" centroid:"<<a_iv.m_volCentroid<<" m_bndryCentroid:"<<a_iv.m_bndryCentroid<<"\n"
      <<"x-arcsLo:"<<a_iv.m_arc[0]<<" x-areaFracsLo:"<<a_iv.m_areaFrac[0]<<"\n"
      <<"y-arcsLo:"<<a_iv.m_arc[1]<<" y-areaFracsLo:"<<a_iv.m_areaFrac[1]<<"\n"
      <<"x-arcsHi:"<<a_iv.m_arc[2]<<" x-areaFracsHi:"<<a_iv.m_areaFrac[2]<<"\n"
      <<"y-arcsHi:"<<a_iv.m_arc[3]<<" y-areaFracsHi:"<<a_iv.m_areaFrac[3]<<"\n";
  return a_os;

  return a_os;
}
/*******************************/
/*******************************/
int IrregNode::
index(int a_idir, Side::LoHiSide a_sd)
{
  CH_assert(a_idir >= 0 && a_idir < SpaceDim);
  int retval;
  if (a_sd == Side::Lo)
    {
      retval = a_idir;
    }
  else
    {
      retval = a_idir + SpaceDim;
    }
  return retval;
}
/*******************************/
/*******************************/

void IrregNode::makeRegular(const IntVect& iv)
{
  m_cell = iv;
  m_volFrac = 1.0;
  m_cellIndex = 0;
  m_volCentroid = IntVect::Zero;
  m_bndryCentroid = IntVect::Zero;
  //low sides
  for (int i=0; i<SpaceDim; i++)
    {
      m_arc[i].resize(1,0);
      m_areaFrac[i].resize(1,1.0);
      RealVect faceCenter = IntVect::Zero;
      faceCenter[i] = -0.5;
      m_faceCentroid[i].resize(1,faceCenter);
    }
  //hi sides
  for (int i=0; i<SpaceDim; i++)
    {
      m_arc[i+SpaceDim].resize(1,0);
      m_areaFrac[i+SpaceDim].resize(1,1.0);
      RealVect faceCenter = IntVect::Zero;
      faceCenter[i] = 0.5;
      m_faceCentroid[i+SpaceDim].resize(1,faceCenter);
    }
}

void IrregNode::faceReserve(int location, int size)
{
  if (m_arc[location].size() < size)
    {
      m_arc[location].resize(size);
      m_areaFrac[location].resize(size);
      m_faceCentroid[location].resize(size);
    }
}
#include "NamespaceFooter.H"
