#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "STLAsciiReader.H"
#include "STLBinaryReader.H"
#include "STLExplorer.H"
#include "STLMesh.H"
#include "STLBox.H"
#include "CellEdge.H"

#include "STLIF.H"

#include "NamespaceHeader.H"

STLIF::STLIF(const char* const      a_filename,
             const STLIF::DataType& a_dataType)
{
  CH_TIME("STLIF::STLIF_file");

  m_filename = a_filename;
  m_dataType = a_dataType;

  m_explorer = NULL;

  makeExplorer();
}

STLIF::STLIF(const STLIF& a_inputIF)
{
  CH_TIME("STLIF::STLIF_copy");

  m_filename = a_inputIF.m_filename;
  m_dataType = a_inputIF.m_dataType;

  m_explorer = NULL;

  makeExplorer();
}

STLIF::~STLIF()
{
#if 0
  if (m_explorer != NULL)
  {
    delete m_explorer;
  }
#endif
}

Real STLIF::value(const RealVect& a_point) const
{
  Real retval = 0.0;

  MayDay::Error("STLIF::value should never be called");

  return retval;
}

BaseIF* STLIF::newImplicitFunction() const
{
  CH_TIME("STLIF::newImplicitFunction");

  STLIF* dataFilePtr = new STLIF(m_filename.c_str(),
                                 m_dataType);

  return static_cast<BaseIF*>(dataFilePtr);
}

STLExplorer* STLIF::getExplorer() const
{
  if (m_explorer == NULL)
  {
    MayDay::Error("STLIF::getExplorer - STLExplorer not defined yet");
  }

  return m_explorer;
}

void STLIF::makeExplorer()
{
  CH_TIME("STLIF::makeExplorer");

  RefCountedPtr<STLMesh> mesh;

  if (m_dataType == STLIF::ASCII)
  {
    STLAsciiReader reader(m_filename);
    mesh = reader.GetMesh();
  }
  else if (m_dataType == STLIF::Binary)
  {
    STLBinaryReader reader(m_filename);
    mesh = reader.GetMesh();
  }
  else
  {
    MayDay::Error("STLIF::makeExplorer - Unknown STL data type");
  }

  m_explorer = new STLExplorer(mesh);
}

#include "NamespaceFooter.H"
