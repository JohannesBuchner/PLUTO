#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Stencil.H"

#include "NamespaceHeader.H"

// IndexML:operator<<
ostream&
operator<< (ostream&       os,
            const IndexML& p)
{
  os << '(' << p.m_iv << ", lev:" << p.level() << /* ", b:" << p.block() << */ ')';
  if (os.fail())
    MayDay::Error("operator<<(ostream&,IndexML&) failed");
  return os;
}
bool
IndexML::operator> (const IndexML& p) const
{
  if ( m_lev < p.m_lev ) return true;
  else if ( m_lev > p.m_lev ) return false;
 
  if ( m_blockID < p.m_blockID ) return true;
  else if ( m_blockID > p.m_blockID ) return false;

  if ( m_iv[0] > p.m_iv[0] ) return true;
  else if (m_iv[0] < p.m_iv[0] ) return false;
#if CH_SPACEDIM > 1
  if ( m_iv[1] > p.m_iv[1] ) return true;
  else if (m_iv[1] < p.m_iv[1] ) return false;
#if CH_SPACEDIM > 2
  if ( m_iv[2] > p.m_iv[2] ) return true;
  else return false;
#endif
#endif
  return false;
}

bool
IndexML::operator< (const IndexML& p) const
{
  if ( m_lev > p.m_lev ) return true;
  else if ( m_lev < p.m_lev ) return false;

  if ( m_blockID > p.m_blockID ) return true;
  else if ( m_blockID < p.m_blockID ) return false;

  if ( m_iv[0] < p.m_iv[0] ) return true;
  else if (m_iv[0] > p.m_iv[0] ) return false;
#if CH_SPACEDIM > 1
  if ( m_iv[1] < p.m_iv[1] ) return true;
  else if (m_iv[1] > p.m_iv[1] ) return false;
#if CH_SPACEDIM > 2
  if ( m_iv[2] < p.m_iv[2] ) return true;
  else return false;
#endif
#endif
  
  return false;
}

bool
IndexML::operator==(const IndexML& p) const
{
  return (p.m_iv==m_iv && p.m_lev==m_lev && p.m_blockID==m_blockID);
}

bool
IndexML::operator!= (const IndexML& p) const
{
  return (p.m_iv!=m_iv || p.m_lev!=m_lev || p.m_blockID!=m_blockID);
}

// StencilTensorValue:operator<<
ostream&
operator<< (ostream&       os,
            const StencilTensorValue& p)
{
  if (p.m_dof==1)
    {
      os << p.value();
    }
  else
    {
      for (int i=0;i<p.m_dof;i++) 
        {
          for (int j=0;j<p.m_dof-1;j++) os << p.value(i,j) << ',';
          os << p.value(i,p.m_dof-1) << std::endl;
        }
    }
  if (os.fail())
    MayDay::Error("operator<<(ostream&,StencilTensorValue&) failed");
  return os;
}

// StencilScalarValue::operator=
StencilScalarValue& StencilScalarValue::operator=(const StencilScalarValue& p)
{
  m_val = p.m_val;
  return *this;
}

// StencilScalarValue:operator<<
ostream&
operator<< (ostream&       os,
            const StencilScalarValue& p)
{
  os << p.value();
  if (os.fail())
    MayDay::Error("operator<<(ostream&,StencilScalarValue&) failed");
  return os;
}

// A specialized operator that "distributes" a_sten[jiv] with a_new: a_sten += a_new*a_sten[jiv]. 
// Would like to remove a_sten[jiv] but that messes up iterators.
void
StencilProject(IndexML a_mliv, Vector<StenScalarNode> &a_scales, StencilTensor &a_sten)
{
  //CH_TIME("PetscCompGrid::projectStencil");
  StencilTensorValue ghostNode = a_sten[a_mliv]; // node getting deleted (distributed)
  // would like to remove the root but this messes up interators
  StencilTensor::iterator root = a_sten.find(a_mliv);
  // add scaled
  for (int i=0;i<a_scales.size();i++) 
    {
      StenScalarNode &target = a_scales[i];
      CH_assert(target.first != a_mliv); // this OK in theory but never done now
      StencilTensorValue &val = a_sten[target.first]; val.define(root->second); // only place where nodes are added
      val.axpy(target.second.value(),root->second);
    }
}

#include "NamespaceFooter.H"
