#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBFLUXREGISTER_H_
#define _EBFLUXREGISTER_H_

#include "REAL.H"
#include "Vector.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "EBISLayout.H"
#include "EBISBox.H"
#include "IntVectSet.H"
#include "CFStencil.H"
#include "LoHiSide.H"
#include "LevelData.H"
#include "LayoutData.H"
#include "EBFastFR.H"
#include "NamespaceHeader.H"
class EBIndexSpace;

class EBCoarToFineRedist;
class EBCoarToCoarRedist;

///EBFluxRegister-A class to encapsulate a levels worth of flux registers.
/**
   A EBFluxRegister  handles all the data choreography
   necessary to create a levels worth of flux registers.
 */
class EBFluxRegister: public EBFastFR
{
public:

  ///
  /**
     Default constructor.  Leaves object undefined.
  */
  EBFluxRegister()
  {;}

  ///
  /**
     Full constructor.  Calls the define function which creates
     a levels worth of flux registers.
  */
  EBFluxRegister(const DisjointBoxLayout& a_dblFine,
                 const DisjointBoxLayout& a_dblCoar,
                 const EBISLayout&        a_ebislFine,
                 const EBISLayout&        a_ebislCoar,
                 const Box&               a_domainCoar,
                 const int&               a_nref,
                 const int&               a_nvar,
                 const EBIndexSpace*      ebisPtr,
                 bool a_forceNoEBCF = false);

  ///
  virtual ~EBFluxRegister();

  ///
  /**
     Full define function.  Creates a levels worth of flux registers.
  */
  void
  define(const DisjointBoxLayout& a_dblFine,
         const DisjointBoxLayout& a_dblCoar,
         const EBISLayout&        a_ebislFine,
         const EBISLayout&        a_ebislCoar,
         const ProblemDomain&     a_domainCoar,
         const int&               a_nref,
         const int&               a_nvar,
         const EBIndexSpace*      ebisPtr,
         bool a_forceNoEBCF = false);



  ///
  /**
     increments the register with data from coarseFlux, multiplied by scale.
     Increment coarse buffer with the area-weighted sum
     of the fluxes on the faces.
     buf += -scale*sum(areaFrac*flux)
     CoarseFlux must contain the coarse fluxes in the dir direction
     for the grid   m_coarseLayout[coarseDataIndex].
  */
  void
  incrementCoarseRegular(const EBFaceFAB&   a_coarseFlux,
                         const Real&        a_scale,
                         const DataIndex&   a_coarseDataIndex,
                         const Interval&    a_variables,
                         const int&         a_dir)
  {
    for (SideIterator sit; sit.ok(); ++sit)
      {
        incrementCoarRegul(a_coarseFlux,
                           a_scale,
                           a_coarseDataIndex,
                           a_variables,
                           a_dir, sit());
      }
  }


  void incrementCoarse(const EBFaceFAB&   a_coarseFlux,
                       const Real&        a_scale,
                       const DataIndex&   a_coarseDataIndex,
                       const Interval&    a_varSrc,
                       const Interval&    a_varDst,
                       const int&         a_dir)
  {
    CH_assert(a_varSrc==a_varDst);
    for (SideIterator sit; sit.ok(); ++sit)
      {
        incrementCoarseBoth(a_coarseFlux, a_scale, a_coarseDataIndex, a_varSrc, a_dir, sit());
      }
  }


  ///
  /**
     increments the register with data from coarseFlux, multiplied by scale.
     Increment coarse buffer with the area-weighted sum
     of the fluxes on the faces.
     buf += -scale*sum(areaFrac*flux)
     CoarseFlux must contain the coarse fluxes in the dir direction
     for the grid   m_coarseLayout[coarseDataIndex].
  */
  void
  incrementCoarseIrregular(const BaseIFFAB<Real>& a_coarseFlux,
                           const Real&            a_scale,
                           const DataIndex&       a_coarseDataIndex,
                           const Interval&        a_variables,
                           const int&             a_dir);

  ///
  /**
     For when you have an EBFaceFAB lying around instead of a BaseIFFAB
   */
  void incrementCoarseIrregular(const EBFaceFAB&       a_coarFlux,
                                const Real&            a_scale,
                                const DataIndex&       a_coarDatInd,
                                const Interval&        a_variables,
                                const int&             a_dir)
  {
    for (SideIterator sit; sit.ok(); ++sit)
      {
        incrementCoarIrreg(a_coarFlux,
                           a_scale,
                           a_coarDatInd,
                           a_variables,
                           a_dir,
                           sit());
      }

  }

  ///
  /**
     increments the register with data from fineFlux, multiplied by scale.
     fineFlux must contain the coarse fluxes in the dir direction
     for the grid m_fineLayout[fineDataIndex].
     This operation is local.
  */
  void
  incrementFineRegular(const EBFaceFAB&      a_fineFlux,
                       const Real&           a_scale,
                       const DataIndex&      a_fineDataIndex,
                       const Interval&       a_variables,
                       const int&            a_dir,
                       const Side::LoHiSide& a_sd)
  {
    incrementFineRegul(a_fineFlux,
                       a_scale,
                       a_fineDataIndex,
                       a_variables,
                       a_dir,
                       a_sd);
  }


  void
  incrementFine(const EBFaceFAB&      a_fineFlux,
                const Real&           a_scale,
                const DataIndex&      a_fineDataIndex,
                const Interval&       a_varSrc,
                const Interval&       a_varDst,
                const int&            a_dir)
  {
    CH_assert(a_varSrc == a_varDst);
    for (SideIterator sit; sit.ok(); ++sit)
      {
        incrementFineBoth(a_fineFlux, a_scale, a_fineDataIndex, a_varSrc, a_dir, sit());
      }
  }





  ///
  /**
     increments the register with data from fineFlux, multiplied by scale.
     fineFlux must contain the coarse fluxes in the dir direction
     for the grid m_fineLayout[fineDataIndex].
     This operation is local.
  */
  void
  incrementFineIrregular(const BaseIFFAB<Real>& a_fineFlux,
                         const Real&            a_scale,
                         const DataIndex&       a_fineDataIndex,
                         const Interval&        a_variables,
                         const int&             a_dir,
                         const Side::LoHiSide&  a_sd);

  ///
  /**
     For when you have a EBFaceFAB lying around instead of a BaseIFFAB
  */
  void
  incrementFineIrregular(const EBFaceFAB&       a_fineFlux,
                         const Real&            a_scale,
                         const DataIndex&       a_fineDataIndex,
                         const Interval&        a_variables,
                         const int&             a_dir,
                         const Side::LoHiSide&  a_sd)
  {
    incrementFineIrreg(a_fineFlux,
                       a_scale,
                       a_fineDataIndex,
                       a_variables,
                       a_dir,
                       a_sd);
  }



  ///
  /**
     Increments redistribution register with left-over mass
     from reflux divergence.
  */
  void incrementRedistRegister(EBCoarToFineRedist& a_register,
                               const Interval&     a_variables,
                               const Real&         a_scale);

  ///
  /**
     Increments redistribution register with left-over mass
     from reflux divergence.
  */
  void incrementRedistRegister(EBCoarToCoarRedist& a_register,
                               const Interval&     a_variables,
                               const Real&         a_scale);


  ///
  bool isDefined() const;

protected:


private:
  //returns true if a_src contains any cells
  bool
  copyBIFFToEBFF(EBFaceFAB&             a_dst,
                 const BaseIFFAB<Real>& a_src,
                 const Box            & a_box,
                 const EBISBox&         a_ebisBox);

  //for all the usual reasons,
  //there is no copy constructor for this class
  //and there is no operator= for this class
  //(there is no rule six).
  void operator= (const EBFluxRegister&)
  {
    MayDay::Error("invalid operator");
  }
  EBFluxRegister(const EBFluxRegister&)
  {
    MayDay::Error("invalid operator");
  }
};
#include "NamespaceFooter.H"
#endif
