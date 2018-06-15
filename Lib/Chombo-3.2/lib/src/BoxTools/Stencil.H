#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _STENCIL_H_
#define _STENCIL_H_

#include "ProblemDomain.H"
#include <map>

//
// Classes general stencil object.  Classes StencilTensor and StencilScalar are the classes defined 
// here along with helper classes for stencil nodes and for multilevel, multi-block indices (IndexML)
//
 
#include "NamespaceHeader.H"

//! \class IndexML
//! Encapsulate a multilevel index
class IndexML
{
public:
  IndexML() : m_lev(0), m_blockID(0)
  {
  }
  IndexML(const IndexML &a_ml) : m_iv(a_ml.m_iv), m_lev(a_ml.m_lev), m_blockID(0)
  {
  }
  IndexML(IntVect a_iv,short a_lev) : m_iv(a_iv), m_lev(a_lev), m_blockID(0)
  {
  }
  void define(IntVect a_iv,short a_lev, short a_blockid=0) 
  {
    m_lev = a_lev;  m_iv = a_iv; m_blockID = a_blockid;
  }
  bool operator== (const IndexML& p) const;
  bool operator!= (const IndexML& p) const;
  bool operator> (const IndexML& p) const;
  bool operator< (const IndexML& p) const;
  short level() const {return m_lev;}
  IntVect iv()const{return m_iv;}
  short block()const{return m_blockID;}
  void setIV(IntVect a_iv){m_iv = a_iv;}
  void setLevel(short a_lev){m_lev = a_lev;}
  void setBlock(short a_bid){m_blockID = a_bid;}
  // Print the IndexML to given output stream in ASCII.
  friend std::ostream& operator<< (std::ostream&       os,
                                   const IndexML& iv);
protected:
  IntVect m_iv;
  short m_lev;
  short m_blockID;
};

#define STENCIL_MAX_DOF 3
#define CHECK_DOF if (m_dof<=0) MayDay::Abort("StencilTensorValue dot defined");
//! \class StencilTensorValue
//! a 2x2 array class
class StencilTensorValue
{
public:
  StencilTensorValue(int a_dof=-1) : m_dof(a_dof) 
  {
    if (m_dof > STENCIL_MAX_DOF) MayDay::Abort("StencilTensorValue::StencilTensorValue dof > MAX");
    if (a_dof>0) for (int i=0;i<m_dof*m_dof;i++) m_val[i]=0.;
  }
  void define(int a_dof) {
    if (a_dof > STENCIL_MAX_DOF) MayDay::Abort("StencilTensorValue::define dof > MAX");
    if (a_dof <= 0) MayDay::Abort("StencilTensorValue::define dof <= 0");
    if (m_dof==-1) for (int i=0;i<a_dof*a_dof;i++) m_val[i]=0.;
    m_dof = a_dof;
  }
  void define(StencilTensorValue &a_sv) {
    define(a_sv.m_dof);
  }
  bool isDefined() const{return (m_dof > 0);}
  StencilTensorValue& operator= (const StencilTensorValue& p) { // a definer
    if (m_dof==-1) m_dof = p.m_dof; 
    if (p.m_dof != m_dof) MayDay::Abort("StencilTensorValue::operator= DOF mismatch");
    for (int i=0;i<m_dof*m_dof;i++) m_val[i] = p.m_val[i];
    return *this;
  }
  Real value(int a_idof=0, int a_jdof=0) const {CHECK_DOF; return m_val[a_idof*m_dof + a_jdof];}
  const Real *getVals() const {CHECK_DOF; return m_val;}
  void addValue(int a_idof, int a_jdof, Real a_val) { CHECK_DOF; m_val[a_idof*m_dof + a_jdof] += a_val; }
  void setValue(int a_idof, int a_jdof, Real a_val) { CHECK_DOF; m_val[a_idof*m_dof + a_jdof] = a_val; }
  void addValue(Real a_val) { CHECK_DOF; m_val[0] += a_val; }
  void setValue(Real a_val) { CHECK_DOF; m_val[0] = a_val; }
  void axpy(Real a_alpha, StencilTensorValue &a_node) {
    CHECK_DOF;
    if (a_node.m_dof != m_dof) MayDay::Abort("StencilTensorValue::axpy DOF mismatch");
    for (int i=0;i<m_dof*m_dof;i++) m_val[i] += a_alpha*a_node.m_val[i];
  }
  void scale(Real a_alpha) { 
    CHECK_DOF; 
    for (int i=0;i<m_dof*m_dof;i++) m_val[i] *= a_alpha;
  }
  // Print the StencilTensorValue to given output stream in ASCII.
  friend std::ostream& operator<< (std::ostream&       os,
                                   const StencilTensorValue& iv);
protected:
  int m_dof;
  Real m_val[STENCIL_MAX_DOF*STENCIL_MAX_DOF];
};

//! \class StencilScalarValue
//! a real array class
class StencilScalarValue
{
public:
  StencilScalarValue(int a_dof=-1) : m_val(0.)
  {}
  Real value() const {return m_val;}
  void addValue(Real a_val) { m_val += a_val; }
  void setValue(Real a_val) { m_val = a_val; }
  void axpy(Real a_alpha, StencilScalarValue &a_node) {
    MayDay::Abort("StencilScalarValue::axpy Used????");
    m_val += a_alpha*a_node.m_val;
  }
  void scale(Real a_alpha) { m_val *= a_alpha; }
  StencilScalarValue& operator= (const StencilScalarValue& p);
  // Print the StencilScalarValue to given output stream in ASCII.
  friend std::ostream& operator<< (std::ostream&       os,
                                   const StencilScalarValue& iv);
protected:
  Real m_val;
};

typedef std::map<IndexML,StencilTensorValue > StencilTensor;
typedef std::pair<IndexML,StencilScalarValue > StenScalarNode;
typedef std::map<IndexML,StencilScalarValue > StencilScalar;

// algebra

// a_sten += a_new*a_sten[jiv]. "distributes" a_sten[jiv] with a_new. Would like to remove a_sten[jiv] but that messes up iterators.
void StencilProject(IndexML jiv, Vector<StenScalarNode> &a_new, StencilTensor &a_sten);

#include "NamespaceFooter.H"

#endif
