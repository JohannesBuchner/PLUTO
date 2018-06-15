#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Stencils.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "NamespaceHeader.H"
Real applyVoFStencil(const VoFStencil& a_sten, const EBCellFAB& a_fab, const int& a_comp)
{
  Real retval = 0.;
  for (int isten = 0; isten < a_sten.size(); isten++)
    {
      retval += (a_sten.weight(isten))*(a_fab((a_sten.vof(isten)), a_sten.variable(isten)));
    }

  return retval;
}
Real applyFaceStencil(const FaceStencil& a_sten, const EBFaceFAB& a_fab, const int& a_comp)
{
  Real retval = 0.;
  for (int isten = 0; isten < a_sten.size(); isten++)
    {
      retval += (a_sten.weight(isten))*(a_fab((a_sten.face(isten)), a_sten.variable(isten)));
    }

  return retval;
}

/**************/
/**************/
VoFStencil::VoFStencil()
  : vofs(0),
    weights(0)
{
}
/**************/
/**************/
void
VoFStencil::clear()
{
  vofs.resize(0);
  variables.resize(0);
  weights.resize(0);
}
/**************/
/**************/
void
FaceStencil::clear()
{
  faces.resize(0);
  weights.resize(0);
  variables.resize(0);
}

/**************/
/**************/
VoFStencil::VoFStencil(const VoFStencil&  stenin)
  : vofs(stenin.vofs),
    weights(stenin.weights),
    variables(stenin.variables)
{
}
/**************/
/**************/
VoFStencil::~VoFStencil()
{
}
/**************/
/**************/
void
VoFStencil::add(const VolIndex& vof,const Real& weight, int ivar)
{
  bool alreadyhere = false;
  for (int ivof = 0; ivof < vofs.size(); ivof++)
    {
      if ((vofs[ivof] == vof) && (variables[ivof] == ivar))
        {
          alreadyhere = true;
          weights[ivof] += weight;
        }
    }
  if (!alreadyhere)
    {
      vofs.push_back(vof);
      weights.push_back(weight);
      variables.push_back(ivar);
    }
}
/**************/
/**************/
VoFStencil&
VoFStencil::operator+=(const VoFStencil& vofstenin)
{
  for (int ivof = 0; ivof < vofstenin.size(); ivof++)
    {
      add(vofstenin.vof(ivof), vofstenin.weight(ivof), vofstenin.variable(ivof));
    }
  return *this;
}
/**************/
/**************/
void
VoFStencil::operator*=(const Real& a_scaling)
{
  for (int ivof = 0; ivof < size(); ivof++)
    {
      weights[ivof] *= a_scaling;
    }
}
/**************/
/**************/
void
FaceStencil::operator*=(const Real& a_scaling)
{
  for (int iface = 0; iface < size(); iface++)
    {
      weights[iface] *= a_scaling;
    }
}
/**************/
/**************/
VoFStencil&
VoFStencil::operator=(const VoFStencil& vofstenin)
{
  clear();
  this->operator+=(vofstenin);
  return *this;
}
/**************/

/**************/
/**************/
FaceStencil::FaceStencil()
  : faces(0),
    weights(0)
{
}
/**************/
/**************/
FaceStencil::FaceStencil(const FaceStencil& facestenin)
  : faces(facestenin.faces),
    weights(facestenin.weights),
    variables(facestenin.variables)
{
}
/**************/
/**************/
// destructor
FaceStencil::~FaceStencil()
{
}
/**************/
/**************/
void
FaceStencil::add(const FaceIndex& face,const Real& weight, int ivar)
{
  bool alreadyhere = false;
  for (int iface = 0; iface < faces.size(); iface++)
    {
      if ( (faces[iface] == face) && (variables[iface] == ivar))
        {
          alreadyhere = true;
          weights[iface] += weight;
        }
    }
  if (!alreadyhere)
    {
      faces.push_back(face);
      weights.push_back(weight);
      variables.push_back(ivar);
    }
}
/**************/
/**************/
int
FaceStencil::size() const
{
  return weights.size();
}
/**************/
/**************/
const FaceIndex&
FaceStencil::face(int i) const
{
  return faces[i];
}
/**************/
/**************/
const Real&
FaceStencil::weight(int i) const
{
  return weights[i];
}
const int&
FaceStencil::variable(int i) const
{
  return variables[i];
}
int&
FaceStencil::variable(int i)
{
  return variables[i];
}
const int&
VoFStencil::variable(int i) const
{
  return variables[i];
}
int&
VoFStencil::variable(int i)
{
  return variables[i];
}
/**************/
/**************/
FaceStencil&
FaceStencil::operator+=(const FaceStencil& Facestenin)
{
  for (int iFace = 0; iFace < Facestenin.size(); iFace++)
    {
      add(Facestenin.face(iFace), Facestenin.weight(iFace));
    }
  return *this;
}
/**************/
/**************/
FaceStencil&
FaceStencil::operator=(const FaceStencil& Facestenin)
{
  clear();
  this->operator+=(Facestenin);
  return *this;
}
#include "NamespaceFooter.H"
