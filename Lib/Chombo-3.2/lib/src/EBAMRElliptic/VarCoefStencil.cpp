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
#include "VarCoefStencil.H"
#include "NamespaceHeader.H"
/**************/
VarCoefStencil::VarCoefStencil()
  : vofs(0),
    weights(0),
    variables(0),
    coefLocs(0)
{
}
/**************/
/**************/
void
VarCoefStencil::clear()
{
  vofs.resize(0);
  variables.resize(0);
  weights.resize(0);
  coefLocs.resize(0);
}
/**************/
VarCoefStencil::VarCoefStencil(const VarCoefStencil&  stenin)
  : vofs(stenin.vofs),
    weights(stenin.weights),
    variables(stenin.variables),
    coefLocs(stenin.coefLocs)
{
}
/**************/
VarCoefStencil::~VarCoefStencil()
{
}
/**************/
/**************/
void
VarCoefStencil::add(const VolIndex& vof,const FaceIndex& coefLoc, const Real& weight, int ivar)
{
  bool alreadyhere = false;
  for (int ivof = 0; (ivof < vofs.size()) && (!alreadyhere); ivof++)
    {
      if ((vofs[ivof] == vof) && (variables[ivof] == ivar) && (coefLocs[ivof] == coefLoc))
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
      coefLocs.push_back(coefLoc);
    }
}
/**************/
/**************/
VarCoefStencil&
VarCoefStencil::operator+=(const VarCoefStencil& vofstenin)
{
  for (int ivof = 0; ivof < vofstenin.size(); ivof++)
    {
      add(vofstenin.vof(ivof), vofstenin.coefLoc(ivof), vofstenin.weight(ivof),  vofstenin.variable(ivof));
    }
  return *this;
}
/**************/
void
VarCoefStencil::operator*=(const Real& a_scaling)
{
  for (int ivof = 0; ivof < size(); ivof++)
    {
      weights[ivof] *= a_scaling;
    }
}
/**************/
VarCoefStencil&
VarCoefStencil::operator=(const VarCoefStencil& vofstenin)
{
  clear();
  this->operator+=(vofstenin);
  return *this;
}
/**************/
const int&
VarCoefStencil::variable(int i) const
{
  return variables[i];
}
/**************/
int&
VarCoefStencil::variable(int i)
{
  return variables[i];
}
#include "NamespaceFooter.H"
