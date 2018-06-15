#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _STENCILS_H_
#define _STENCILS_H_

#include "Vector.H"
#include "VolIndex.H"
#include "BaseIndex.H"
#include "FaceIndex.H"
#include "REAL.H"
#include "NamespaceHeader.H"
///
/**
   class for vofstencil and basestencil to inherit from
   so AggStencil can unify them.
*/
class BaseStencil
{
public:
  ///
  BaseStencil()
  {
  }

  ///
  virtual ~BaseStencil()
  {
  }

  ///
  virtual int size() const = 0;

  ///
  virtual const BaseIndex& index(int isten) const = 0;


  ///
  virtual const Real& weight(int isten) const = 0;

  ///
  virtual const int& variable(int isten) const = 0;
};

/// VoF-centered stencil
/**
   This stencil is a fundamental tool for building
   eb applications by making the ability to cache
   stencil information. This object consists of
   a vector of VoFs and a vector of weights.
*/
class VoFStencil: public BaseStencil
{
public:
  ///
  /**
     default constructor.  makes empty vectors.
  */
  VoFStencil();

  ///
  virtual ~VoFStencil();

  ///
  void clear();

  ///
  /**
     add a VoF to the Stencil, with it's associated weight
     it is required that all VoFs are in the same EBIndexSpace
     If the vof is already in the
     stencil, add the weights.
  */
  void
  add(const VolIndex& vof,const Real& weight, int ivar = 0);

  ///
  /**
     number of VoFs in the Stencil
  */
  virtual int size() const;

  ///
  /**
     access a VoF
  */
  const VolIndex&  vof(int i) const;


  ///
  virtual BaseIndex& index(int i) const
  {
    return (BaseIndex&)(vof(i));
  }

  ///
  /**
     access a weight
  */
  virtual const Real& weight(int i) const;

  ///
  /**
     access a weight
  */
  virtual Real& weight(int i);


  ///
  virtual const int& variable(int i) const;

  ///
  virtual int& variable(int i);

  ///
  /**
     add all faces and weights of inputs
     to this.  If a vof is already in the
     stencil, add the weights.
     only addition is well-defined here
     as far as arithmatic operations are
     concerned.
  */
  VoFStencil&
  operator+=(const VoFStencil& a_vofsten);

  ///
  /**
   */
  void operator*=(const Real& scaling);

  ///
  VoFStencil&
  operator=(const VoFStencil& a_vofsten);

  ///
  VoFStencil(const VoFStencil& a_vofstenin);

  void setAllVariables(int a_var)
  {
    if (vofs.size() > 0)
      {
        variables = Vector<int>(vofs.size(), a_var);
      }
  }
protected:

  /// the VoFs
  Vector<VolIndex> vofs;
  /// the weights
  Vector<Real> weights;

  Vector<int> variables;
};

///
/**
   Face-centered Stencil for embedded boundary applications.
*/
class FaceStencil: public BaseStencil
{
public:
  ///
  FaceStencil();

  ///
  virtual ~FaceStencil();

  ///
  void clear();

  ///
  /**
     add an Face and it's weight
     If the face is already in the
     stencil, add the weights.
  */
  void  add(const FaceIndex& face,const Real& weight, int variable=0);

  ///
  /**
     number of Faces in the Stencil
  */
  virtual int   size() const;

  ///
  void setAllVariables(int a_var)
  {
    if (faces.size() > 0)
      {
        variables = Vector<int>(faces.size(), a_var);
      }
  }

  ///
  /** access an Face
   */
  const FaceIndex&  face(int i) const;

  ///
  /**
     access a weight
  */
  virtual const Real& weight(int i) const;

  ///
  virtual BaseIndex& index(int i) const
  {
    return (BaseIndex&)(face(i));
  }

  ///
  virtual const int& variable(int i) const;

  ///
  virtual int& variable(int i);

  ///
  /**
     add all faces and weights of inputs
     to this.  If a face is already in the
     stencil, add the weights.
     only addition is well-defined here
     as far as arithmatic operations are
     concerned.
  */
  FaceStencil& operator+=(const FaceStencil& a_facesten);

  ///
  FaceStencil& operator=(const FaceStencil& a_facesten);

  ///
  FaceStencil(const FaceStencil& a_facesten);

  ///
  /**
   */
  void operator*=(const Real& scaling);

private:
  /// the Faces
  Vector<FaceIndex> faces;
  /// the weights
  Vector<Real> weights;
  /// the variable numbers
  Vector<int> variables;

};

class EBCellFAB;
class EBFaceFAB;
///variable argument ignored---now uses the variables that live in the stencil
extern Real applyVoFStencil(const VoFStencil& a_vofSten, const EBCellFAB& a_fab, const int& a_comp);

///variable argument ignored---now uses the variables that live in the stencil
extern Real applyFaceStencil(const FaceStencil& a_faceSten, const EBFaceFAB& a_fab, const int& a_comp);

/** inlines */

/**************/
inline int
VoFStencil::size() const
{
  return weights.size();
}
/**************/
/**************/
inline const  VolIndex&
VoFStencil::vof(int i) const
{
  return vofs[i];
}
/**************/
/**************/
inline const Real&
VoFStencil::weight(int i) const
{
  return weights[i];
}

inline Real&
VoFStencil::weight(int i)
{
  return weights[i];
}

#ifdef CH_EXPLICIT_TEMPLATES
//
// Extra no-op stuff required for successful explicit template instantiation.
//
inline int linearSize( const FaceStencil& fs )
{
  CH_assert(0);
  return -1;
}
template<> inline
void linearIn<FaceStencil>(FaceStencil& a_outputT, const void* const inBuf)
{
  CH_assert(0);
}
template<> inline
void linearOut<FaceStencil>(void* const a_outBuf, const FaceStencil& a_inputT)
{
  CH_assert(0);
}


inline int linearSize( const VoFStencil& fs )
{
  CH_assert(0);
  return -1;
}
template<> inline
void linearIn<VoFStencil>(VoFStencil& a_outputT, const void* const inBuf)
{
  CH_assert(0);
}
template<> inline
void linearOut<VoFStencil>(void* const a_outBuf, const VoFStencil& a_inputT)
{
  CH_assert(0);
}


inline int linearSize( const bool& fs )
{
  CH_assert(0);
  return -1;
}
template<> inline
void linearIn<bool>(bool& a_outputT, const void* const inBuf)
{
  CH_assert(0);
}
template<> inline
void linearOut<bool>(void* const a_outBuf, const bool& a_inputT)
{
  CH_assert(0);
}


#endif // CH_EXPLICIT_TEMPLATES


#include "NamespaceFooter.H"
#endif
