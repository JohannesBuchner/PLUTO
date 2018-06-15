#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBFASTFR_H_
#define _EBFASTFR_H_

#include "REAL.H"
#include "Vector.H"
#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "EBCellFactory.H"
#include "EBLevelDataOps.H"
#include "EBISLayout.H"
#include "EBISBox.H"
#include "IntVectSet.H"
#include "CFStencil.H"
#include "LoHiSide.H"
#include "LevelData.H"
#include "EBLevelGrid.H"
#include "LayoutData.H"
#include "LevelFluxRegister.H"
#include "NamespaceHeader.H"
class EBIndexSpace;

///EBFastFR-A class to encapsulate a levels worth of flux registers.
/**
   A EBFastFR  handles all the data choreography
   necessary to create a levels worth of flux registers.
   This only does the elliptic flux register stuff for now (no interacting with
   redistribution, no RZ, no regular and regular separate).
*/
class EBFastFR
{
public:
  static bool s_verbose;
  ///
  /**
     Default constructor.  Leaves object undefined.
  */
  EBFastFR();

  ///
  /**
     Full constructor.  Calls the define function which creates
     a levels worth of flux registers.
     noebcf forces no coarse-fine intersection with EB
  */
  EBFastFR(const EBLevelGrid&       a_eblgFine,
           const EBLevelGrid&       a_eblgCoar,
           const int&               a_nref,
           const int&               a_nvar,
           bool a_forceNoEBCF = false);

  ///
  virtual ~EBFastFR();

  ///
  /**
     Full define function.  Creates a levels worth of flux registers.
  */
  virtual void
  define(const EBLevelGrid&       a_eblgFine,
         const EBLevelGrid&       a_eblgCoar,
         const int&               a_nref,
         const int&               a_nvar,
         bool  a_forceNoEBCF = false);

  ///
  /**
     Initialize values of registers  to zero.
  */
  virtual void
  setToZero();


  ///
  /**
     The irregular part of this is just
     register = input flux  (instead of increment)
     To avoid double-counting.
   */
  virtual void
  incrementCoarseBoth(const EBFaceFAB&      a_coarseFlux,
                      const Real&           a_scale,
                      const DataIndex&      a_coarseDataIndex,
                      const Interval&       a_variables,
                      const int&            a_dir,
                      const Side::LoHiSide& a_sd);


  ///
  virtual void
  incrementFineBoth(const EBFaceFAB&      a_fineFlux,
                    const Real&           a_scale,
                    const DataIndex&      a_fineDataIndex,
                    const Interval&       a_variables,
                    const int&            a_dir,
                    const Side::LoHiSide& a_sd);

  void
  compareFineSparse(const EBFaceFAB&      a_fluxOld,
                    const EBFaceFAB&      a_fluxNew,
                    const DataIndex&      a_fineDatInd,
                    const int&            a_dir,
                    const Side::LoHiSide& a_sd);

  ///to support baseiffab approach
  virtual void
  incrementFineSparse(const EBFaceFAB&      a_fineFlux,
                      const Real&           a_scale,
                      const DataIndex&      a_fineDatInd,
                      const Interval&       a_variables,
                      const int&            a_dir,
                      const Side::LoHiSide& a_sd,
                      bool a_doingFineRegular);

  ///
  virtual void reflux(LevelData<EBCellFAB>& a_uCoarse,
                      const Interval&       a_variables,
                      const Real&           a_scale,
                      bool a_multByKappaOneMinusKappa = false);

  ///
  virtual void reflux(LevelData<EBCellFAB>& a_uCoarse,
                      const Interval&       a_solutionvariables,
                      const Interval&       a_fluxvariables,
                      const Real&           a_scale,
                      bool a_multByKappaOneMinusKappa = false)
  {
    CH_assert(a_fluxvariables.size() == a_solutionvariables.size());
    EBCellFactory fact(m_eblgCoar.getEBISL());
    LevelData<EBCellFAB> diffCoar(m_eblgCoar.getDBL(), m_nComp, IntVect::Zero, fact);
    EBLevelDataOps::setToZero(diffCoar);
    reflux(diffCoar, a_fluxvariables, a_scale, a_multByKappaOneMinusKappa);
    for (DataIterator dit = m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
      {
        int isrc = a_fluxvariables.begin();
        int inco = a_fluxvariables.size();
        int idst = a_solutionvariables.begin();
        a_uCoarse[dit()].plus(diffCoar[dit()], isrc, idst, inco);
      }
  }

  ///
  /**
     Increments mass array  with left-over mass
     from reflux divergence.   this is to test this pig.
     Ordinarily this mass would go into redistribution.   this
     is used in the test suite to  check whether a constant
     flux refluxed and then unrefluxed (this function) ends up
     with a constant solution.
     Correction at each cell = (1-kappa)refluxCorrection.
  */
  void incrementDensityArray(LevelData<EBCellFAB>& a_coarMass,
                             const Interval&       a_variables,
                             const Real&           a_scale);

  ///
  bool isDefined() const;

  static int index(int a_dir, Side::LoHiSide a_side);

  void incrementFineIrreg(const EBFaceFAB&       a_fineFlux,
                          const Real&            a_newScale,
                          const DataIndex&       a_fineDatInd,
                          const Interval&        a_variables,
                          const int&             a_dir,
                          const Side::LoHiSide&  a_sd);

  ///
  /**
     The irregular part of this is just
     register = input flux  (instead of increment)
     To avoid double-counting.
     The side argument is ignored.
   */
  void incrementCoarIrreg(const EBFaceFAB&       a_coarFlux,
                          const Real&            a_scale,
                          const DataIndex&       a_coarDatInd,
                          const Interval&        a_variables,
                          const int&             a_dir,
                          const Side::LoHiSide&  a_sd);

  void incrementFineRegul(const EBFaceFAB&       a_fineFlux,
                          const Real&            a_newScale,
                          const DataIndex&       a_fineDatInd,
                          const Interval&        a_variables,
                          const int&             a_dir,
                          const Side::LoHiSide&  a_sd);

  void incrementCoarRegul(const EBFaceFAB&       a_coarFlux,
                          const Real&            a_scale,
                          const DataIndex&       a_coarDatInd,
                          const Interval&        a_variables,
                          const int&             a_dir,
                          const Side::LoHiSide&  a_sd);

  bool hasEBCF() const
  {
    return m_hasEBCF;
  }

  ///undefines object
  void clear();

  ///
  /**
     Gets the set of vofs to fill a flux for on a particular side
     and direction that correspond to an EB/CF interface.   You still
     need to fill all the other fluxes at the interface but you do not
     have to do all the fancy interpolation to face centroids
     that is needed at these key points since all the other C/F points
     are treated as regular anyway.
     Do not use this unless you really know what you are doing.   Really.
     Yes this breaks encapsulation but we do some pretty ugly things
     for performance reasons, do we not?
   */
  Vector<VoFIterator>& getVoFItCoar(const DataIndex&       a_dit,
                                    const int&             a_idir,
                                    const Side::LoHiSide&  a_sd)
  {
    int iloc = index(a_idir, a_sd);
    return (m_vofiCoar[iloc])[a_dit];
  }

  ///
  /**
     Gets the set of vofs to fill a flux for on a particular side
     and direction that correspond to an EB/CF interface.   You still
     need to fill all the other fluxes at the interface but you do not
     have to do all the fancy interpolation to face centroids
     that is needed at these key points since all the other C/F points
     are treated as regular anyway.
     Do not use this unless you really know what you are doing.   Really.
     Yes this breaks encapsulation but we do some pretty ugly things
     for performance reasons, do we not?
   */
  VoFIterator& getVoFItCoFi(const DataIndex&       a_dit,
                            const int&             a_idir,
                            const Side::LoHiSide&  a_sd)
  {
    int iloc = index(a_idir, a_sd);
    return (m_vofiCoFi[iloc])[a_dit];
  }

  ///
  const EBLevelGrid& getEBLGCoFi() const
  {
    return m_eblgCoFi;
  }

protected:

  void  irregSetToZero();
  void  defineSetsAndIterators();
  void  defineBuffers();

  void cacheOldSolution(const LevelData<EBCellFAB>& a_uCoar,
                        const Interval&             a_variables);

  void restoreOldSolution(LevelData<EBCellFAB>&       a_uCoar,
                          const Interval&             a_variables);

  void irregReflux(LevelData<EBCellFAB>& a_uCoar,
                   const Interval&       a_variables,
                   const Real&           a_scale,
                   bool a_multByKappaOneMinusKappa = false);

  void  incrementByRefluxDivergence(LevelData<EBCellFAB>& a_uCoar,
                                    LevelData<EBCellFAB>& a_fluxDiff,
                                    const Interval      & a_variables,
                                    const Real          & a_newScale,
                                    bool a_multByOneMinusKappa,
                                    bool a_multByKappaOneMinusKappa);


  LevelFluxRegister* m_levelFluxReg;
  bool m_isDefined;

  void setDefaultValues();
  EBLevelGrid       m_eblgFine;
  EBLevelGrid       m_eblgCoar;
  EBLevelGrid       m_eblgCoFi;

  int               m_refRat;
  int               m_nComp;

  //nrefdmo is the number of fine faces per coar face
  //this is intrinsically dimension-dependent
  Real m_nrefdmo;
  bool m_hasEBCF;


  bool computeHasEBCF();
  LayoutData< VoFIterator >          m_vofiCoFi[2*SpaceDim];
  LayoutData< IntVectSet  >          m_setsCoFi[2*SpaceDim];
  LayoutData< Vector<IntVectSet > >  m_setsCoar[2*SpaceDim];
  LayoutData< Vector<VoFIterator> >  m_vofiCoar[2*SpaceDim];

  Copier m_reverseCopier;

  //need to save stuff
  //so we can do this in two stages.
  LevelData<EBCellFAB> m_saveCoar;

  //fine fluxes*fine areas/nref^d-1
  LevelData<EBCellFAB> m_delUCoFi;
  LevelData<EBCellFAB> m_delUCoar;
  LevelData<EBCellFAB> m_delUDiff;
private:
  //for all the usual reasons,
  //there is no copy constructor for this class
  //and there is no operator= for this class
  //(there is no rule six).
  void operator= (const EBFastFR& out)
  {
    MayDay::Error("invalid operator");
  }
  EBFastFR(const EBFastFR& out)
  {
    MayDay::Error("invalid operator");
  }
};
#include "NamespaceFooter.H"
#endif
