#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBSTENVARCOEF_H_
#define _EBSTENVARCOEF_H_

#include "Stencils.H"
#include "Vector.H"
#include "VolIndex.H"
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "EBIndexSpace.H"
#include "EBCellFAB.H"
#include "NamespaceHeader.H"

///
/**
   EBStencil variant specifically for the needs of single-level variable
   coefficent elliptic operators.   There is only one way to make it.
   Everyone has to have the same number of ghost cells.   There is
   only one way to use it.   If you want something more flexible, use
   EBStencil.
 */
class EBStenVarCoef
{
public:
  ///
  /**
     Destructor
  */
  ~EBStenVarCoef()
  {;}

  ///
  /**
    srcVofs    = list of vofs for which the stencil will be applied
    vofstencil = stencil for the divergence(F)
    validBox   = valid cells (box from the DBL)
    ghostVect  = number of ghost cells (for both phi, L(phi), rho,
                  lambda, everything)
    ebisbox    = geometric information
    varDest    = variable of lphi
  */
  EBStenVarCoef(const Vector<VolIndex>&      a_srcVofs,
                const BaseIVFAB<VoFStencil>& a_vofstencil,
                const Box&                   a_validBox,
                const EBISBox&               a_ebisbox,
                const IntVect&               a_ghostVect,
                int a_varDest);



  ///
  /**
     lphi =  alpha*alphaWeight*phi + beta*betaWeight*(div(F))
     The offsets of this object compute divF.
     lofphi      = result
     phi         = phi
     alphaWeight = variable coefficient factor of the identity term
     alpha       = constant coefficient factor of the identity term
     betaWeight  = variable coefficient factor of the div(F) term
     beta        = constant coefficient factor of the div(F) term
  */
  void apply(EBCellFAB&             a_lofphi,
             const EBCellFAB&       a_phi,
             const EBCellFAB&       a_alphaWeight,
             const Real&            a_alpha,
             const EBCellFAB&       a_betaWeight,
             const Real&            a_beta);


  ///
  /**
      phi := phi + lambda*(lphi-rhs)
      where
      lphi =  alpha*alphaWeight*phi + beta*betaWeight*(div(F))
      The offsets of this object compute divF.
      lofphi      = result
      phi         = phi
      alphaWeight = variable coefficient factor of the identity term
      alpha       = constant coefficient factor of the identity term
      betaWeight  = variable coefficient factor of the div(F) term
      beta        = constant coefficient factor of the div(F) term
      lambda      = relaxation coefficient
   */
  void relax(EBCellFAB&  a_phi,
       const EBCellFAB&  a_rhs,
       const EBCellFAB&  a_alphaWeight,
       const EBCellFAB&  a_betaWeight,
       const EBCellFAB&  a_lambda,
       Real a_alpha, Real a_beta) const;

  ///
  /**
     Cache the irregular values of the input
  */
  void cache(const EBCellFAB& a_phi, int ivar);

   ///
   /**
      take values of internal cache and put them into phi
   */
  void uncache(EBCellFAB& a_phi, int ivar) const;

  struct
  {
    int offset;
    bool multiValued;
  } typedef destTerm_t;

  struct
  {
    int offset;
    Real weight;
    //only used when debugging
    //VolIndex vof;
  } typedef stencilTerm;

  struct
  {
    Vector<stencilTerm> single;
    Vector<stencilTerm> multi;
  } typedef varcsten_t;

protected:
  void computeOffsets(const Vector<VolIndex>& a_srcVoFs, const BaseIVFAB<VoFStencil>& a_vofstencil);

  Box     m_box;
  EBISBox m_ebisBox;
  Box     m_grownBox;
  IntVect m_ghostVect;
  int     m_destVar;
  Vector<varcsten_t>   m_stencil;
  Vector<destTerm_t>   m_destTerms;
  Vector<destTerm_t>   m_sourTerms;
  mutable Vector<Real> m_cache;
  mutable Vector<Real> m_cachePhi;
  //debugging hook
  Vector<VolIndex>     m_srcVoFs;
private:
  /// weak construction is bad
  EBStenVarCoef()
  {
    MayDay::Error("invalid operator");
  }

  /// object contains pointered data.  assignment is a bad idea
  void operator=(const EBStenVarCoef& stenin)
  {
    MayDay::Error("invalid operator");
  }

  /// object contains pointered data.  copy construction is a bad idea
  EBStenVarCoef(const EBStenVarCoef& stenin)
  {
    MayDay::Error("invalid operator");
  }

};


#include "NamespaceFooter.H"
#endif
