#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VARCOEFSTENCIL_H_
#define _VARCOEFSTENCIL_H_

#include "Vector.H"
#include "VolIndex.H"
#include "FaceIndex.H"
#include "REAL.H"
#include "NamespaceHeader.H"

/// VoF-centered stencil
/**
   This stencil is a fundamental tool for building
   eb applications by making the ability to cache
   stencil information. This object consists of
   a vector of VoFs and a vector of weights.
 */
class VarCoefStencil
{
public:
  ///
  /**
     default constructor.  makes empty vectors.
  */
  VarCoefStencil();

  ///
  ~VarCoefStencil();

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
  add(const VolIndex& vof,const FaceIndex& coefloc, const Real& weight, int ivar);

  ///
  /**
     number of VoFs in the Stencil
  */
  int
  size() const;

  ///
  /**
     access a VoF
  */
  const VolIndex&
  vof(int i) const;

  ///
  /**
     access a VoF
  */
  const FaceIndex&
  coefLoc(int i) const;


  ///
  /**
      access a weight
   */
  const Real&
  weight(int i) const;

  ///
  /**
      access a weight
   */
  Real& weight(int i);


  ///
  const int& variable(int i) const;

  ///
  int& variable(int i);

  ///
  /**
     add all faces and weights of inputs
     to this.  If a vof is already in the
     stencil, add the weights.
     only addition is well-defined here
     as far as arithmatic operations are
     concerned.
   */
  VarCoefStencil&
  operator+=(const VarCoefStencil& a_vofsten);

  ///
  /**
   */
  void operator*=(const Real& scaling);

  ///
  VarCoefStencil&
  operator=(const VarCoefStencil& a_vofsten);

  ///
  VarCoefStencil(const VarCoefStencil& a_vofstenin);

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

  /// Where the face-centered coefficients live
  Vector<FaceIndex> coefLocs;
};


/** inlines */

/**************/
inline int
VarCoefStencil::size() const
{
  return weights.size();
}
/**************/
inline const  FaceIndex&
VarCoefStencil::coefLoc(int i) const
{
  return coefLocs[i];
}
/**************/
inline const  VolIndex&
VarCoefStencil::vof(int i) const
{
  return vofs[i];
}
/**************/
inline const Real&
VarCoefStencil::weight(int i) const
{
  return weights[i];
}

inline Real&
VarCoefStencil::weight(int i)
{
  return weights[i];
}

#ifdef CH_EXPLICIT_TEMPLATES
//
// Extra no-op stuff required for successful explicit template instantiation.
//


inline int linearSize( const VarCoefStencil& fs )
{
    CH_assert(0);
    return -1;
}
template<> inline
void linearIn<VarCoefStencil>(VarCoefStencil& a_outputT, const void* const inBuf)
{
    CH_assert(0);
}
template<> inline
void linearOut<VarCoefStencil>(void* const a_outBuf, const VarCoefStencil& a_inputT)
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
