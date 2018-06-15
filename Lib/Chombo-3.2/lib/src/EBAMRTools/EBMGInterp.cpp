#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves fri aug 17, 2001

#include "EBMGInterp.H"
#include "EBMGInterpF_F.H"
#include "VoFIterator.H"
#include "DisjointBoxLayout.H"
#include "EBCellFactory.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "CH_Timer.H"
#include "EBLevelGrid.H"
#include "NamespaceHeader.H"

/************************************/
void
EBMGInterp::setDefaultValues()
{
  m_doLinear = false;
  m_isDefined = false;
  m_refRat = -1;
  m_nComp = -1;
}
/************************************/
EBMGInterp::EBMGInterp()
{
  setDefaultValues();
}
/************************************/
EBMGInterp::~EBMGInterp()
{
}
// /************************************/
EBMGInterp::EBMGInterp(const DisjointBoxLayout& a_dblFine,
                       const DisjointBoxLayout& a_dblCoar,
                       const EBISLayout& a_ebislFine,
                       const EBISLayout& a_ebislCoar,
                       const ProblemDomain& a_domainCoar,
                       const int& a_nref,
                       const int& a_nvar,
                       const EBIndexSpace*  ebisPtr,
                       const IntVect& a_ghost,
                       const bool& a_layoutChanged,
                       const bool& a_doLinear)
{
  setDefaultValues();

  define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar,
         a_domainCoar, a_nref, a_nvar, ebisPtr, a_ghost,
         a_layoutChanged, a_doLinear);

}
/************************************/
void
EBMGInterp::define(const DisjointBoxLayout&    a_dblFine,
                   const DisjointBoxLayout&    a_dblCoar,
                   const EBISLayout&           a_ebislFine,
                   const EBISLayout&           a_ebislCoar,
                   const ProblemDomain&        a_domainCoar,
                   const int&                  a_nref,
                   const int&                  a_nvar,
                   const EBIndexSpace*         ebisPtr,
                   const IntVect&              a_ghostCellsPhi,
                   const bool&                 a_layoutChanged,
                   const bool&                 a_doLinear)
{
  CH_TIMERS("EBMGInterp::define");
  CH_TIMER("fillEBISLayout", t1);
  m_isDefined = true;
  m_doLinear = a_doLinear;
  m_ghost = a_ghostCellsPhi;
  m_nComp = a_nvar;
  m_coarGrids = a_dblCoar;
  m_fineGrids = a_dblFine;
  m_coarEBISL = a_ebislCoar;
  m_fineEBISL = a_ebislFine;
  m_coarDomain = a_domainCoar;
  m_refRat = a_nref;
  m_fineDomain = refine(m_coarDomain, m_refRat);
  m_layoutChanged = a_layoutChanged;
  m_coarsenable   = a_dblFine.coarsenable(m_refRat);

  //only define ebislbuf and gridbuf if we are changing layouts
  if (m_layoutChanged)
    {
      ProblemDomain domebisl;
      if (m_coarsenable)
        {
          coarsen(m_buffGrids, m_fineGrids, m_refRat);
          domebisl = m_coarDomain;
        }
      else
        {
          refine(m_buffGrids,  m_coarGrids, m_refRat);
          m_copierRCtoF.define(m_buffGrids, m_fineGrids, a_ghostCellsPhi);
          m_copierFtoRC.define(m_fineGrids, m_buffGrids, a_ghostCellsPhi);
          domebisl = m_fineDomain;
        }

      CH_START(t1);
      int nghost = 4;
      ebisPtr->fillEBISLayout(m_buffEBISL,
                              m_buffGrids,
                              domebisl, nghost);
      if (m_refRat > 2)
        {
          if (m_coarsenable)
            {
              m_buffEBISL.setMaxRefinementRatio(m_refRat, ebisPtr);
            }
          else
            {
              m_buffEBISL.setMaxCoarseningRatio(m_refRat, ebisPtr);
            }
        }
      CH_STOP(t1);
    }
  defineStencils();
}

void
EBMGInterp::
defineStencils()
{
  CH_TIME("EBMGInterp::defineStencils");

  DisjointBoxLayout gridsStenCoar;
  EBISLayout        ebislStenCoar;
  DisjointBoxLayout gridsStenFine;
  EBISLayout        ebislStenFine;
  if (m_layoutChanged)
    {
      if (m_coarsenable)
        {
          gridsStenCoar = m_buffGrids;
          ebislStenCoar = m_buffEBISL;
          gridsStenFine = m_fineGrids;
          ebislStenFine = m_fineEBISL;
        }
      else
        {
          gridsStenCoar = m_coarGrids;
          ebislStenCoar = m_coarEBISL;
          gridsStenFine = m_buffGrids;
          ebislStenFine = m_buffEBISL;
        }

    }
  else
    {
      gridsStenCoar = m_coarGrids;
      ebislStenCoar = m_coarEBISL;
      gridsStenFine = m_fineGrids;
      ebislStenFine = m_fineEBISL;
    }

  {
    CH_TIME("graph walking");

    m_interpEBStencil.define(gridsStenCoar);
    if (m_doLinear)
      {
        m_linearEBStencil.define(gridsStenCoar);
      }

    for (DataIterator dit = gridsStenCoar.dataIterator(); dit.ok(); ++dit)
      {
        const EBISBox& ebisBoxCoar =  ebislStenCoar[dit()];
        const EBISBox& ebisBoxFine =  ebislStenFine[dit()];
        const Box&         boxFine =  gridsStenFine[dit()];
        const EBGraph& ebGraphFine =  ebisBoxFine.getEBGraph();

        const Box& boxCoar         = gridsStenCoar[dit()];
        IntVectSet notRegularCoar = ebisBoxCoar.getIrregIVS(boxCoar);

        IntVectSet ivsFine = refine(notRegularCoar, m_refRat);
        VoFIterator vofItIrregFine(ivsFine, ebGraphFine);
        const Vector<VolIndex>& allFineVoFs = vofItIrregFine.getVector();

        defineConstantStencil(dit(), allFineVoFs,
                              ebislStenFine, ebislStenCoar, boxFine, boxCoar);
        if (m_doLinear)
          {
            defineLinearStencil(dit(), allFineVoFs,
                                ebislStenFine, ebislStenCoar,
                                boxFine, boxCoar);
          }
      }
  }
}
void
EBMGInterp::
defineLinearStencil(const DataIndex       & a_dit,
                    const Vector<VolIndex>& a_allFineVoFs,
                    const EBISLayout      & a_ebislStenFine,
                    const EBISLayout      & a_ebislStenCoar,
                    const Box             & a_boxFine,
                    const Box             & a_boxCoar)

{

  const EBISBox& ebisBoxCoar =  a_ebislStenCoar[a_dit];
  const EBISBox& ebisBoxFine =  a_ebislStenFine[a_dit];
  Vector<VoFStencil> interpStencils(a_allFineVoFs.size());

  for (int ifine = 0; ifine < a_allFineVoFs.size(); ifine++)
    {
      VolIndex coarVoF;
      const VolIndex& fineVoF = a_allFineVoFs[ifine];
      if (m_refRat > 2)
        {
          coarVoF = a_ebislStenFine.coarsen(fineVoF, m_refRat, a_dit);
        }
      else
        {
          coarVoF = ebisBoxFine.coarsen(fineVoF);
        }

      Real dxFine = 1.;  Real dxCoar = m_refRat;
      RealVect coarLoc =  EBArith::getVofLocation(coarVoF, dxCoar*RealVect::Unit, RealVect::Zero);
      RealVect fineLoc =  EBArith::getVofLocation(fineVoF, dxFine*RealVect::Unit, RealVect::Zero);
      RealVect dist = fineLoc - coarLoc;
      VoFStencil& interpStencil = interpStencils[ifine];

      ///linear stencil holds dist*slope
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Vector<FaceIndex> facesLo = ebisBoxCoar.getFaces(coarVoF, idir, Side::Lo);
          Vector<FaceIndex> facesHi = ebisBoxCoar.getFaces(coarVoF, idir, Side::Hi);
          VoFStencil stenLo;
          for (int ivof = 0; ivof < facesLo.size(); ivof++)
            {
              stenLo.add(facesLo[ivof].getVoF(Side::Lo), -1/dxCoar);
              stenLo.add(coarVoF,                         1/dxCoar);
            }
          if (facesLo.size() > 1) stenLo *= Real(1.0/facesLo.size());
          VoFStencil stenHi;
          for (int ivof = 0; ivof < facesHi.size(); ivof++)
            {
              stenHi.add(facesHi[ivof].getVoF(Side::Hi), 1/dxCoar);
              stenHi.add(coarVoF,                       -1/dxCoar);
            }
          if (facesHi.size() > 1) stenHi *= Real(1.0/facesHi.size());
          VoFStencil slopeStencil;
          bool hasLo = (facesLo.size() > 0) && (!facesLo[0].isBoundary());
          bool hasHi = (facesHi.size() > 0) && (!facesHi[0].isBoundary());
          if (hasLo && hasHi)
            {
              slopeStencil += stenLo;
              slopeStencil += stenHi;
              slopeStencil *= 0.5;
            }
          else if (hasLo)
            {
              slopeStencil = stenLo;
            }
          else if (hasHi)
            {
              slopeStencil = stenHi;
            }

          //multiply by distance in slope direction
          slopeStencil  *= dist[idir];
          interpStencil += slopeStencil;
        }
    }

  m_linearEBStencil[a_dit] = RefCountedPtr<EBStencil>(new EBStencil(a_allFineVoFs, interpStencils,
                                                                    a_boxFine, a_boxCoar,
                                                                    ebisBoxFine, ebisBoxCoar,
                                                                    m_ghost, m_ghost));

}
void
EBMGInterp::
defineConstantStencil(const DataIndex       & a_dit,
                      const Vector<VolIndex>& a_allFineVoFs,
                      const EBISLayout      & a_ebislStenFine,
                      const EBISLayout      & a_ebislStenCoar,
                      const Box             & a_boxFine,
                      const Box             & a_boxCoar)

{

  const EBISBox& ebisBoxCoar =  a_ebislStenCoar[a_dit];
  const EBISBox& ebisBoxFine =  a_ebislStenFine[a_dit];
  Vector<VoFStencil> interpStencils(a_allFineVoFs.size());

  for (int ifine = 0; ifine < a_allFineVoFs.size(); ifine++)
    {
      VolIndex coarVoF;
      if (m_refRat > 2)
        {
          coarVoF = a_ebislStenFine.coarsen(a_allFineVoFs[ifine], m_refRat, a_dit);
        }
      else
        {
          coarVoF = ebisBoxFine.coarsen(a_allFineVoFs[ifine]);
        }

      VoFStencil& interpStencil = interpStencils[ifine];

      //weight is one because there is only one vof to the stencil
      //and we are doing piecewise constant interpolation
      //for the copy stencil the weight is one because it is a copy
      interpStencil.add(coarVoF, 1.0);
    }

  m_interpEBStencil[a_dit] = RefCountedPtr<EBStencil>(new EBStencil(a_allFineVoFs, interpStencils,
                                                                    a_boxFine, a_boxCoar,
                                                                    ebisBoxFine, ebisBoxCoar,
                                                                    m_ghost, m_ghost));

}
/************************************/
bool
EBMGInterp::isDefined() const
{
  return m_isDefined;
}

/************************************/
class EBAddOp : public LDOperator<EBCellFAB>
{
public:

  virtual void linearIn(EBCellFAB& arg,  void* buf, const Box& R,
                        const Interval& comps) const
  {
    EBCellFAB tmp;
    tmp.clone(arg);
    tmp.linearIn(buf, R, comps);
    arg.plus(tmp, R, comps.begin(), comps.begin(), comps.size());
  }

  void op(EBCellFAB& dest,
          const Box& RegionFrom,
          const Interval& Cdest,
          const Box& RegionTo,
          const EBCellFAB& src,
          const Interval& Csrc) const
  {
    dest.plus(src, RegionFrom, Csrc.begin(), Cdest.begin(), Cdest.size());
  }

};
void
EBMGInterp::pwcInterp(LevelData<EBCellFAB>&       a_fineData,
                      const LevelData<EBCellFAB>& a_coarData,
                      const Interval&             a_variables)
{
  CH_TIME("EBMGInterp::pwcInterp");
  CH_assert(a_fineData.ghostVect() == m_ghost);
  CH_assert(a_coarData.ghostVect() == m_ghost);

  if (m_layoutChanged)
    {
      if (m_coarsenable)
        {
          CH_TIME("EBMGInterp::pwcInterp::coarsenable");
          EBCellFactory ebcellfact(m_buffEBISL);

          LevelData<EBCellFAB> coarsenedFineData(m_buffGrids, m_nComp, m_ghost, ebcellfact);
          a_coarData.copyTo(a_variables, coarsenedFineData, a_variables);
          for (DataIterator dit = m_fineGrids.dataIterator(); dit.ok(); ++dit)
            {
              //does incrementonly = true
              pwcInterpFAB(a_fineData[dit()],
                           m_buffGrids[dit()],
                           coarsenedFineData[dit()],
                           dit(),
                           a_variables);
            }
        }
      else
        {
          CH_TIME("EBMGInterp::pwcInterp::uncoarsenable");
          EBCellFactory ebcellfact(m_buffEBISL);
          LevelData<EBCellFAB> refinedCoarseData(m_buffGrids, m_nComp, m_ghost, ebcellfact);
          for (DataIterator dit = m_coarGrids.dataIterator(); dit.ok(); ++dit)
            {
              refinedCoarseData[dit()].setVal(0.);
              pwcInterpFAB(refinedCoarseData[dit()],
                           m_coarGrids[dit()],
                           a_coarData[dit()],
                           dit(),
                           a_variables);
            }

          EBAddOp op;
          refinedCoarseData.copyTo(a_variables, a_fineData, a_variables, m_copierRCtoF, op);
        }
    }
  else
    {
      pwcInterpMG(a_fineData, a_coarData, a_variables);
    }
}
/*******/
void
EBMGInterp::pwlInterp(LevelData<EBCellFAB>&       a_fineData,
                      const LevelData<EBCellFAB>& a_coarData,
                      const Interval&             a_variables)
{
  CH_TIME("EBMGInterp::pwlInterp");
  CH_assert(a_fineData.ghostVect() == m_ghost);
  CH_assert(a_coarData.ghostVect() == m_ghost);
  CH_assert(m_doLinear); //otherwise stencils have not been defined

  if (m_layoutChanged)
    {
      if (m_coarsenable)
        {
          CH_TIME("EBMGInterp::pwlInterp::coarsenable");
          EBCellFactory ebcellfact(m_buffEBISL);
          LevelData<EBCellFAB> coarsenedFineData(m_buffGrids, m_nComp, m_ghost, ebcellfact);
          a_coarData.copyTo(a_variables, coarsenedFineData, a_variables);
          fillGhostCellsPWC(coarsenedFineData, m_buffEBISL, m_coarDomain);
          for (DataIterator dit = m_fineGrids.dataIterator(); dit.ok(); ++dit)
            {
              //does incrementonly = true
              pwlInterpFAB(a_fineData[dit()],
                           m_buffGrids[dit()],
                           coarsenedFineData[dit()],
                           dit(),
                           a_variables);
            }
        }
      else
        {
          CH_TIME("EBMGInterp::pwlInterp::uncoarsenable");
          EBCellFactory ebcellfact(m_buffEBISL);
          fillGhostCellsPWC((LevelData<EBCellFAB>&)a_coarData, m_coarEBISL, m_coarDomain);
          LevelData<EBCellFAB> refinedCoarseData(m_buffGrids, m_nComp, m_ghost, ebcellfact);
          for (DataIterator dit = m_coarGrids.dataIterator(); dit.ok(); ++dit)
            {
              refinedCoarseData[dit()].setVal(0.);
              pwlInterpFAB(refinedCoarseData[dit()],
                           m_coarGrids[dit()],
                           a_coarData[dit()],
                           dit(),
                           a_variables);
            }

          EBAddOp op;
          refinedCoarseData.copyTo(a_variables, a_fineData, a_variables, m_copierRCtoF, op);
        }
    }
  else
    {
      pwcInterpMG(a_fineData, a_coarData, a_variables);
    }
}
/************************************/
void
EBMGInterp::fillGhostCellsPWC(LevelData<EBCellFAB>& a_data,
                              const EBISLayout&     a_ebisl,
                              const ProblemDomain&  a_dom)
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (m_ghost[idir] < 1)
        {
          MayDay::Error("EBMGInterp:I need a ghost cell for linear interpolation");
        }
    }
  //extrapolate to every ghost cell
  const DisjointBoxLayout& dbl = a_data.disjointBoxLayout();
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = dbl.get(dit());
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Side::LoHiSide flipSide = flip(sit());
              Box ghostBox = adjCellBox(grid, idir, sit(), 1);

              for (BoxIterator boxit(ghostBox); boxit.ok(); ++boxit)
                {
                  if (a_dom.contains(boxit()))
                    {
                      Vector<VolIndex> vofs = ebgraph.getVoFs(boxit());
                      for (int ivof = 0; ivof < vofs.size(); ivof++)
                        {
                          Real extrapVal = 0;
                          Vector<FaceIndex> faces = ebgraph.getFaces(vofs[ivof], idir, flipSide);
                          for (int ivar = 0; ivar < a_data.nComp(); ivar++)
                            {
                              for (int iface = 0; iface < faces.size(); iface++)
                                {
                                  const VolIndex& flipVoF = faces[iface].getVoF(flipSide);
                                  extrapVal += a_data[dit()](flipVoF, ivar);
                                }
                              if (faces.size() > 1) extrapVal /= faces.size();
                              a_data[dit()](vofs[ivof], ivar) = extrapVal;
                            }
                        }

                    }
                  else
                    {
                      //just put in the single valued part of the data holder something sensible
                      //if not inside domain
                      int isign = sign(flipSide);
                      BaseFab<Real>& bfdata = a_data[dit()].getSingleValuedFAB();
                      IntVect otherIV = boxit() + isign*BASISV(idir);
                      for (int ivar = 0; ivar < a_data.nComp(); ivar++)
                        {
                          bfdata(boxit(), ivar) = bfdata(otherIV, ivar);
                        }
                    }
                }
            }
        }
    }
  //do exchange to overwrite ghost cells where there exists real data
  a_data.exchange();
}
/************************************/
void
EBMGInterp::pwlInterpMG(LevelData<EBCellFAB>&       a_fineData,
                        const LevelData<EBCellFAB>& a_coarData,
                        const Interval&             a_variables)
{
  CH_TIME("EBMGInterp::pwlInterpMG");
  CH_assert(a_fineData.ghostVect() == m_ghost);
  CH_assert(a_coarData.ghostVect() == m_ghost);
  CH_assert(m_doLinear); //otherwise stencils have not been defined


  for (DataIterator dit = m_coarGrids.dataIterator(); dit.ok(); ++dit)
    {
      pwlInterpFAB(a_fineData[dit()],
                   m_coarGrids[dit()],
                   a_coarData[dit()],
                   dit(),
                   a_variables);
    }
}
/************************************/
void
EBMGInterp::pwcInterpMG(LevelData<EBCellFAB>&       a_fineData,
                        const LevelData<EBCellFAB>& a_coarData,
                        const Interval&             a_variables)
{
  CH_TIME("EBMGInterp::pwcInterpMG");
  CH_assert(a_fineData.ghostVect() == m_ghost);
  CH_assert(a_coarData.ghostVect() == m_ghost);

  for (DataIterator dit = m_coarGrids.dataIterator(); dit.ok(); ++dit)
    {
      pwcInterpFAB(a_fineData[dit()],
                   m_coarGrids[dit()],
                   a_coarData[dit()],
                   dit(),
                   a_variables);
    }
}
/************************************/
void
EBMGInterp::pwcInterpFAB(EBCellFAB&       a_refCoar,
                         const Box&       a_coarBox,
                         const EBCellFAB& a_coar,
                         const DataIndex& a_datInd,
                         const Interval&  a_variables) const
{
  CH_TIMERS("EBMGInterp::interp");
  CH_TIMER("regular_interp", t1);
  CH_TIMER("irregular_interp", t2);
  CH_assert(isDefined());

  const Box& coarBox = a_coarBox;

  for (int ivar = a_variables.begin();  ivar <= a_variables.end(); ivar++)
    {
      m_interpEBStencil[a_datInd]->cache(a_refCoar, ivar);

      //do all cells as if they were regular
      Box refBox(IntVect::Zero, IntVect::Zero);
      refBox.refine(m_refRat);

      const BaseFab<Real>& coarRegFAB =    a_coar.getSingleValuedFAB();
      BaseFab<Real>& refCoarRegFAB    = a_refCoar.getSingleValuedFAB();

      CH_START(t1);

      FORT_REGPROLONG(CHF_FRA1(refCoarRegFAB,ivar),
                      CHF_CONST_FRA1(coarRegFAB,ivar),
                      CHF_BOX(coarBox),
                      CHF_BOX(refBox),
                      CHF_CONST_INT(m_refRat));

      CH_STOP(t1);

      m_interpEBStencil[a_datInd]->uncache(a_refCoar, ivar);

      CH_START(t2);
      m_interpEBStencil[a_datInd]->apply(a_refCoar, a_coar, true, ivar);
      CH_STOP(t2);
    }
}

/************************************/
void
EBMGInterp::pwlInterpFAB(EBCellFAB&       a_refCoar,
                         const Box&       a_coarBox,
                         const EBCellFAB& a_coar,
                         const DataIndex& a_datInd,
                         const Interval&  a_variables) const
{
  CH_TIMERS("EBMGInterp::pwlinterpfab");
  CH_TIMER("regular_interp", t1);
  CH_TIMER("irregular_interp", t2);
  CH_assert(isDefined());


  //first interpolate piecewise constant.
  pwcInterpFAB(a_refCoar,
               a_coarBox,
               a_coar,
               a_datInd,
               a_variables);

  //then add in slope*distance.
  const Box& coarBox = a_coarBox;

  for (int ivar = a_variables.begin();  ivar <= a_variables.end(); ivar++)
    {
      //save stuff at irreg cells because fortran result will be garbage
      //and this is an incremental process
      m_linearEBStencil[a_datInd]->cache(a_refCoar, ivar);

      //do all cells as if they were regular
      Box refBox(IntVect::Zero, IntVect::Zero);
      refBox.refine(m_refRat);

      const BaseFab<Real>& coarRegFAB =    a_coar.getSingleValuedFAB();
      BaseFab<Real>& refCoarRegFAB    = a_refCoar.getSingleValuedFAB();

      CH_START(t1);

      Real dxf = 1.0/m_coarDomain.size(0);
      Real dxc = 2.0*dxf;

      //do every cell as regular
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FORT_PROLONGADDSLOPE(CHF_FRA1(refCoarRegFAB,ivar),
                               CHF_CONST_FRA1(coarRegFAB,ivar),
                               CHF_BOX(coarBox),
                               CHF_BOX(refBox),
                               CHF_INT(idir),
                               CHF_REAL(dxf),
                               CHF_REAL(dxc),
                               CHF_CONST_INT(m_refRat));
        }
      CH_STOP(t1);

      //replace garbage fortran put in with original values (at irregular cells)
      m_linearEBStencil[a_datInd]->uncache(a_refCoar, ivar);

      CH_START(t2);
      //do what fortran should have done at irregular cells.
      m_linearEBStencil[a_datInd]->apply(a_refCoar, a_coar, true, ivar);
      CH_STOP(t2);
    }
}
/************************************/
#include "NamespaceFooter.H"
