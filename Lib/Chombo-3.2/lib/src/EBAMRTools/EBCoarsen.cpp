#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBCoarsen.H"
#include "EBAverageF_F.H"
#include "VoFIterator.H"
#include "BoxIterator.H"
#include "DisjointBoxLayout.H"
#include "EBCellFactory.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"
/************************************/
void
EBCoarsen::setDefaultValues()
{
  m_isDefined = false;
  m_refRat = -1;
  m_nComp = -1;
}
/************************************/
EBCoarsen::EBCoarsen()
{
  setDefaultValues();
}
/************************************/
EBCoarsen::~EBCoarsen()
{
}
/************************************/
EBCoarsen::EBCoarsen(const EBLevelGrid&             a_eblgFine,
                     const EBLevelGrid&             a_eblgCoar,
                     const int&                     a_nref,
                     const int&                     a_nvar)
{
  CH_TIME("EBCoarsen:EBCoarsen");
  setDefaultValues();

  define(a_eblgFine, a_eblgCoar,a_nref, a_nvar);
}
/************************************/
void
EBCoarsen::define(const EBLevelGrid&             a_eblgFine,
                  const EBLevelGrid&             a_eblgCoar,
                  const int&                     a_nref,
                  const int&                     a_nvar)
{
  CH_TIME("EBCoarsen:define");
  CH_assert(a_nref > 0);
  CH_assert(a_nvar > 0);
  CH_assert(a_eblgFine.getDBL().coarsenable(a_nref));
  CH_assert(a_eblgFine.getEBIS()->isDefined());
  CH_assert(a_eblgCoar.getEBIS()->isDefined());
  m_cfivsPtr = a_eblgFine.getCFIVS();
  m_isDefined  = true;
  m_refRat     = a_nref;
  m_nComp      = a_nvar;
  m_gridsCoar  = a_eblgCoar.getDBL();
  m_gridsFine  = a_eblgFine.getDBL();
  m_ebislCoar  = a_eblgCoar.getEBISL();
  m_ebislFine  = a_eblgFine.getEBISL();
  m_domainCoar = a_eblgCoar.getDomain();
  m_domainFine = a_eblgFine.getDomain();

  IntVect ghost = 4*IntVect::Unit;
  int nghost = 4;

  m_coarsenedFineGrids = DisjointBoxLayout();
  coarsen(m_coarsenedFineGrids, m_gridsFine, m_refRat);
  a_eblgFine.getEBIS()->fillEBISLayout(m_coarsenedFineEBISL,
                                       m_coarsenedFineGrids,
                                       m_domainCoar, nghost);
  m_coarsenedFineEBISL.setMaxRefinementRatio(m_refRat, a_eblgFine.getEBIS());
  EBCellFactory ebcellfact(m_coarsenedFineEBISL);
  m_coarsenedFineData.define(m_coarsenedFineGrids, m_nComp,
                             ghost, ebcellfact);

  //define the coarsening stencil for irreg vofs and
  //  for coarse vofs next to the cf interface if refRat<4
  defineStencil(*a_eblgFine.getCFIVS());
}
/************************************/
bool
EBCoarsen::isDefined() const
{
  return m_isDefined;
}
/************************************/
void
EBCoarsen::defineStencil(const LayoutData<IntVectSet>&  a_cfivs)
{
  m_coarsenStencil.define(m_coarsenedFineGrids);
  m_vofIt.define(m_coarsenedFineGrids);
  //make the coarsening stencils
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox&    ebisBoxCoar = m_coarsenedFineEBISL[dit()];
      const EBISBox&    ebisBoxFine = m_ebislFine[dit()];
      const Box&         dblBoxCoar = m_coarsenedFineGrids.get(dit());
      //      const Box&         dblBoxFine = m_gridsFine.get(dit());
      const IntVectSet&       cfivs = a_cfivs[dit()];

      //we need to do the coarsening with special care at:
      //  irregular coarse vofs
      //  if refRat < 4, then at coarse vofs (covered by finer) next to the cf interface or domain edge
      const IntVectSet     irregIVS = ebisBoxCoar.getIrregIVS(dblBoxCoar);

      BaseIVFAB<VoFStencil>& coarsenStencilBaseIVFAB = m_coarsenStencil[dit()];
      coarsenStencilBaseIVFAB.define(irregIVS,ebisBoxCoar.getEBGraph(), 1);

      //cache the coarse vof iterator
      m_vofIt[dit()].define(irregIVS,ebisBoxCoar.getEBGraph());

      //loop over the special vofs and build the coarsening stencil
      VoFIterator& vofit = m_vofIt[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& vofCoar = vofit();
          VoFStencil& coarsenStencil = coarsenStencilBaseIVFAB(vofCoar,0);

          getCoarsenVoFStencil(coarsenStencil,ebisBoxCoar,ebisBoxFine,vofCoar,dit(),cfivs);
          checkStencil(coarsenStencil,1.0);
        }
    }
}
void EBCoarsen::
checkStencil(VoFStencil& a_stencil,
             const Real& a_sum) const
{
  Real sum = 0.0;
  for (int i = 0; i < a_stencil.size(); ++i )
    {
      const Real& weight = a_stencil.weight(i);
      sum += weight;
    }

#ifdef CH_USE_FLOAT
  Real tolerance = 1.0e-6;
#else
  Real tolerance = 1.0e-9;
#endif

  CH_assert(Abs(sum-a_sum) < tolerance);
}
void EBCoarsen::
getCoarsenVoFStencil(VoFStencil&       a_stencil,
                     const EBISBox&    a_ebisBoxCoar,
                     const EBISBox&    a_ebisBoxFine,
                     const VolIndex&   a_vofCoar,
                     const DataIndex&  a_datInd,
                     const IntVectSet& a_cfivs)
{
  const IntVect ivCoar = a_vofCoar.gridIndex();

  //figure out the lower left fine iv for this coarse iv
  IntVect loFineIV = ivCoar*m_refRat;
  RealVect centerX = RealVect(loFineIV);
  centerX += ( (Real(m_refRat) - 1.0)/2.0 )*RealVect::Unit;//this is relative to the fine indexspace

  //get the fineVofs for a_vofCoar
  Vector<VolIndex> fineVoFs = m_coarsenedFineEBISL.refine(a_vofCoar, m_refRat, a_datInd);
  //go to every fine vof and do a Taylor expansion to the coarse cell center
  int numPts = 0;
  IntVect order = IntVect::Zero;//at the very least we can do piece-wise constant
  for (int ifine = 0; ifine < fineVoFs.size(); ifine++)
    {
      const VolIndex& fineVoF = fineVoFs[ifine];
      const IntVect&  ivFine = fineVoF.gridIndex();

      //figure out the distance to the coarse cell center from this fine center
      RealVect ivFineX = RealVect(ivFine);
      RealVect dist = centerX - ivFineX;

      //compute the gradient stencil
      //distances here are rigged s.t. dxfine =1
      RealVect dx = RealVect::Unit;
      VoFStencil extrapSten;

      EBArith::getExtrapolationStencil(extrapSten, dist, dx, fineVoF, a_ebisBoxFine, -1, &((*m_cfivsPtr)[a_datInd]));
      a_stencil += extrapSten;

      numPts++;
    }

  if (numPts > 0)
    {
      a_stencil *= 1.0/Real(numPts);
    }
  else
    {
      MayDay::Error("unexpected case::EBCoarsen coarsenFine");
    }
}
/***********************/
/************************************/
void
EBCoarsen::coarsenIrreg(EBCellFAB&       a_coar,
                        const EBCellFAB& a_fine,
                        const DataIndex& a_dit,
                        const Interval&  a_variables)
{
  const BaseIVFAB<VoFStencil>& stenBaseIVFAB = m_coarsenStencil[a_dit];
  VoFIterator& vofit = m_vofIt[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit)
    {
      const VolIndex&    vofCoar   = vofit();
      const VoFStencil&  stencil   = stenBaseIVFAB(vofCoar,0);

      for (int icomp=a_variables.begin();icomp ==a_variables.end();icomp++)
        {
          //coarsen irreg fine vofs to coarse
          //  compute the coarsening by using fine data
          Real phi = 0.0;
          for (int i = 0; i < stencil.size(); ++i )
            {
              const Real&      weight = stencil.weight(i);
              const VolIndex& vofFine = stencil.vof(i);
              const Real&        phiF = a_fine(vofFine,icomp);
              phi += weight * phiF;
            }

          //set the coarse value
          a_coar(vofCoar,icomp) = phi;
        }
    }
}
/************************************/
void
EBCoarsen::coarsenFine(LevelData<EBCellFAB>& a_coarData,
                       const LevelData<EBCellFAB>& a_fineData,
                       const Interval& a_variables)
{
  CH_assert(isDefined());
  //first coarsen the data onto the coarsenedFine data
  //then copy over
  for (DataIterator fineit = m_gridsFine.dataIterator();
      fineit.ok(); ++fineit)
    {
      coarsenFAB(m_coarsenedFineData[fineit()],
                 a_fineData[fineit()],
                 fineit(),
                 a_variables);
    }

  //this is the part that is blocking
  m_coarsenedFineData.copyTo(a_variables, a_coarData, a_variables);
}
/************************************/
void
EBCoarsen::coarsenFAB(EBCellFAB&       a_coar,
                      const EBCellFAB& a_fine,
                      const DataIndex& a_datInd,
                      const Interval&  a_variables)
{
  CH_assert(isDefined());
  //do all cells as if they were regular
  BaseFab<Real>& coarRegFAB =  a_coar.getSingleValuedFAB();
  const BaseFab<Real>& fineRegFAB = a_fine.getSingleValuedFAB();

  //this is how much we need to grow a coarse box to get a fine box that
  // only has the fine iv's that are neighbors of the coarse iv
  const int fac = 1 - m_refRat/2;//NOTE: fac = 0 for refRat==2...
  Box refbox(IntVect::Zero,
             (m_refRat-1)*IntVect::Unit);
  refbox.grow(fac);

  Box coarBox = m_coarsenedFineGrids.get(a_datInd);

  Box fineBox = refine(coarBox, m_refRat);
  CH_assert(coarRegFAB.box().contains(coarBox));
  CH_assert(fineRegFAB.box().contains(fineBox));

  for (int ivar = a_variables.begin();
      ivar <= a_variables.end(); ivar++)
    {
      BaseFab<Real> laplFine(fineBox, 1);
      laplFine.setVal(0.);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Box loBox, hiBox, centerBox;
          int hasLo, hasHi;
          EBArith::loHiCenter(loBox, hasLo,
                              hiBox, hasHi,
                              centerBox, m_domainFine,
                              fineBox, idir, &((*m_cfivsPtr)[a_datInd]));

          FORT_H2LAPL1DADDITIVE(CHF_FRA1(laplFine, 0),
                                CHF_FRA1(fineRegFAB, ivar),
                                CHF_CONST_INT(idir),
                                CHF_BOX(loBox),
                                CHF_CONST_INT(hasLo),
                                CHF_BOX(hiBox),
                                CHF_CONST_INT(hasHi),
                                CHF_BOX(centerBox));
        }

      FORT_EBCOARSEN(CHF_FRA1(coarRegFAB,ivar),
                     CHF_CONST_FRA1(fineRegFAB,ivar),
                     CHF_CONST_FRA1(laplFine, 0),
                     CHF_BOX(coarBox),
                     CHF_CONST_INT(m_refRat),
                     CHF_BOX(refbox));
    }

  //overwrite irregular vofs and coarse vofs next to the cfivs if refRat<4,
  //  for these vofs we coarsen based on
  //  taylor expansions from fine vofs to the coarse cell center
  coarsenIrreg(a_coar,a_fine,a_datInd,a_variables);
}
/************************************/
#include "NamespaceFooter.H"
