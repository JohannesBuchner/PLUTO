#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBSTENCIL_H_
#define _EBSTENCIL_H_

#include "Stencils.H"
#include "Vector.H"
#include "VolIndex.H"
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "EBIndexSpace.H"
#include "EBCellFAB.H"
#include "NamespaceHeader.H"

/// EB stencil
/**
 */
class EBStencil
{
public:
  ///
  /**
     Destructor
  */
  ~EBStencil();

  ///
  /**
     Full constructor.  Box is unghosted.
     When setting Lphi(i,j) = L(phi(i,j))
     GhostVectPhi is the ghost vect of phi
     GhostVectLph is the ghost vect of lphi
  */
  EBStencil(const Vector<VolIndex>&      a_srcVofs,
            const BaseIVFAB<VoFStencil>& a_vofstencil,
            const Box&                   a_box,
            const EBISBox&               a_ebisbox,
            const IntVect&               a_ghostVectPhi,
            const IntVect&               a_ghostVectLph,
            int a_varDest = 0,
            bool a_doRelaxOpt = false,
            int ncomp = 1,
            IntVectSet a_setIrreg = IntVectSet(),
            bool       a_useInputSets = false);

  ///
  /**
     Second constructor for creating stencil
     where dest and src live on 2 different levels
     and dest is Lphi and src is Phi.
     no relax optimization here
  */
  EBStencil(const Vector<VolIndex>& a_srcVofs,
            const Vector<VoFStencil>& a_vofStencil,
            const Box& a_boxLph,
            const Box& a_boxPhi,
            const EBISBox& a_ebisBoxLph,
            const EBISBox& a_ebisBoxPhi,
            const IntVect& a_ghostVectLph,
            const IntVect& a_ghostVectPhi,
            int a_varDest = 0,
            int ncomp = 1,
            IntVectSet a_setIrreg = IntVectSet(),
            bool       a_useInputSets = false);



  ///
  /**
     Applies stencil to each component of phi using the stencil weights and offsets to compute L.
     If incrementOnly = true, a_lofphi is incremented without any set to zero.
     If false, a_lofphi is set to zero and set equal to L(phi).
     ivar is so you can apply a scalar ebstencil to a component of a
     larger holder
  */
  void
  apply(EBCellFAB& a_lofphi, const EBCellFAB& a_phi, bool incrementOnly = false, int ivar = 0) const;


  ///
  /**
     Applies stencil to each component of phi using the stencil weights and offsets to compute L.
     but with constants.
     a_lofphi_i = alpha*phi_i + beta*lphi_i
     If incrementOnly = true, a_lofphi is incremented without any set to zero.
     If false, a_lofphi is set to zero and set equal to a_lofphi_i
     Alpha and  beta are defined over getIrregIVS(lphBox) where lphBox = grow(a_box, a_ghostVectLph)
     where a_box are given in the constructor.
  */
  void
  apply(EBCellFAB& a_lofphi, const EBCellFAB& a_phi,
        const BaseIVFAB<Real>& a_alphaWeight,
        Real a_alpha, Real a_beta, bool incrementOnly = false) const;

  void
  apply(EBCellFAB&             a_lofphi,
        const EBCellFAB&       a_phi,
        const Real             a_lambdaFactor,
        const Real             a_alpha,
        const BaseIVFAB<Real>& a_alphaWeight,
        const Real             a_beta,
        const BaseIVFAB<Real>& a_betaWeight,
        Real                   a_one,
        bool incrementOnly = false) const;

  void
  applyInhomDomBC(EBCellFAB&             a_lofphi,
                  const EBCellFAB&       a_phi,
                  const Real             a_factor) const;

  ///
  void
  relax(EBCellFAB& a_phi,
        const EBCellFAB& a_rhs,
        const BaseIVFAB<Real>& a_alphaWeight,
        const BaseIVFAB<Real>&  a_betaWeight,
        Real a_alpha, Real a_beta, Real a_safety) const;

  ///
  void
  relaxClone(EBCellFAB& a_phi,
             const EBCellFAB& a_phiOld,
             const EBCellFAB& a_rhs,
             const BaseIVFAB<Real>& a_alphaWeight,
             const BaseIVFAB<Real>&  a_betaWeight,
             Real a_alpha, Real a_beta, Real a_safety) const;

  ///
  /**
     Cache lphi.
     ivar is so you can apply a scalar ebstencil to a component of a
     larger holder  */
  void
  cache(const EBCellFAB& a_lphi, int a_ivar = 0) const;

  ///
  /**
     Cache phi
     ivar is so you can apply a scalar ebstencil to a component of a
     larger holder  */
  void
  cachePhi(const EBCellFAB& a_lphi, int a_ivar = 0) const;


  ///
  /**
     Uncache lphi
     ivar is so you can apply a scalar ebstencil to a component of a
     larger holder  */
  void
  uncache(EBCellFAB& a_lphi, int  a_ivar = 0) const;

  ///
  /**
     Uncache phi
     ivar is so you can apply a scalar ebstencil to a component of a
     larger holder  */
  void
  uncachePhi(EBCellFAB& a_phi, int a_ivar = 0) const;


  ///
  /**
     Compute integer offsets for stencil operations
  */
  void
  computeOffsets(const Vector<VolIndex>& a_srcVoFs, const BaseIVFAB<VoFStencil>& a_vofstencil);


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
  } typedef ebstencil_t;

protected:

  Box m_box;
  EBISBox m_ebisBox;
  Box m_lphBox;
  Box m_phiBox;
  IntVect m_ghostVectPhi;
  IntVect m_ghostVectLph;
  int m_destVar;
  Vector<ebstencil_t> m_ebstencil;
  Vector<destTerm_t>  m_destTerms;
  Vector<destTerm_t>  m_sourTerms;
  Vector<int>         m_alphaBeta;
  mutable Vector<Real> m_cacheLph;
  mutable Vector<Real> m_cachePhi;
  bool m_doRelaxOpt;
  int m_nComp;
  IntVectSet m_setIrreg;
  bool       m_useInputSets;
  //only used when debugging
  //Vector<VolIndex> m_srcVoFs;
private:
  ///
  /*
     Default constructor.
  */
  EBStencil()
  {
    MayDay::Error("invalid operator");
  }

  ///
  /*
     Operator =.
  */
  void
  operator=(const EBStencil& stenin)
  {
    MayDay::Error("invalid operator");
  }

  ///
  /*
     Copy constructor.
  */
  EBStencil(const EBStencil& stenin)
  {
    MayDay::Error("invalid operator");
  }

};


#include "NamespaceFooter.H"
#endif
