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

#include "EBMGAverage.H"
#include "EBMGAverageF_F.H"
#include "VoFIterator.H"
#include "DisjointBoxLayout.H"
#include "EBCellFactory.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "CH_Timer.H"
#include "EBLevelDataOps.H"
#include "NamespaceHeader.H"
/************************************/
void
EBMGAverage::setDefaultValues()
{
  m_isDefined = false;
  m_refRat = -1;
  m_nComp = -1;
}
/************************************/
EBMGAverage::EBMGAverage()
{
  setDefaultValues();
}
/************************************/
EBMGAverage::~EBMGAverage()
{
}
 /************************************/
EBMGAverage::EBMGAverage(const DisjointBoxLayout& a_dblFine,
                         const DisjointBoxLayout& a_dblCoar,
                         const EBISLayout& a_ebislFine,
                         const EBISLayout& a_ebislCoar,
                         const ProblemDomain& a_domainCoar,
                         const int& a_nref,
                         const int& a_nvar,
                         const EBIndexSpace* ebisPtr,
                         const IntVect& a_ghostCellsRHS,
                         const bool& a_layoutChanged)
{
  setDefaultValues();

  define(a_dblFine, a_dblCoar, a_ebislFine, a_ebislCoar,
         a_domainCoar, a_nref, a_nvar, ebisPtr, a_ghostCellsRHS,
         a_layoutChanged);

}
/************************************/
void
EBMGAverage::define(const DisjointBoxLayout& a_dblFine,
                    const DisjointBoxLayout& a_dblCoar,
                    const EBISLayout&        a_ebislFine,
                    const EBISLayout&        a_ebislCoar,
                    const ProblemDomain&     a_domainCoar,
                    const int&               a_nref,
                    const int&               a_nvar,
                    const EBIndexSpace*      ebisPtr,
                    const IntVect&           a_ghostCellsRHS,
                    const bool&              a_layoutChanged)
{
  CH_TIMERS("EBMGAverage::define");
  CH_TIMER("fillEBISLayout", t1);
  CH_assert(a_nref > 0);
  CH_assert(a_nvar > 0);

  m_ghost = a_ghostCellsRHS;
  m_isDefined = true;
  m_refRat = a_nref;
  m_nComp = a_nvar;
  m_coarGrids = a_dblCoar;
  m_fineGrids = a_dblFine;
  m_coarEBISL = a_ebislCoar;
  m_fineEBISL = a_ebislFine;
  m_coarDomain = a_domainCoar;
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
          m_copier.define(m_buffGrids, m_coarGrids, a_ghostCellsRHS);
          domebisl = m_coarDomain;
        }
      else
        {
          refine(m_buffGrids,  m_coarGrids, m_refRat);
          m_copier.define(a_dblFine, m_buffGrids, a_ghostCellsRHS);
          domebisl = m_fineDomain;
        }

      CH_START(t1);
      int nghost = 4;
      //pout() << "mg ave dbl = " << m_buffGrids << endl;
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
EBMGAverage::
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

    LayoutData<Vector<VoFStencil> > averageStencil;
    LayoutData<VoFIterator > vofItIrregCoar;
    LayoutData<VoFIterator > vofItIrregFine;
    m_averageEBStencil.define(gridsStenCoar);
    averageStencil.define(    gridsStenCoar);
    vofItIrregFine.define(    gridsStenCoar); //not wrong.  trust me.
    vofItIrregCoar.define(    gridsStenCoar);

    for (DataIterator dit = gridsStenCoar.dataIterator(); dit.ok(); ++dit)
      {
        Box refBox(IntVect::Zero, IntVect::Zero);
        refBox.refine(m_refRat);
        int numFinePerCoar = refBox.numPts();

        const     Box&     boxFine = gridsStenFine[dit()];
        const EBISBox& ebisBoxFine = ebislStenFine[dit()];
        const EBGraph& ebGraphFine = ebisBoxFine.getEBGraph();

        const         Box& boxCoar = gridsStenCoar[dit()];
        const EBISBox& ebisBoxCoar = ebislStenCoar[dit()];
        const EBGraph& ebGraphCoar = ebisBoxCoar.getEBGraph();

        IntVectSet notRegularCoar = ebisBoxCoar.getIrregIVS(boxCoar);
        vofItIrregCoar[dit()].define(notRegularCoar, ebGraphCoar);

        IntVectSet ivsFine = refine(notRegularCoar, m_refRat);
        vofItIrregFine[dit()].define(ivsFine, ebGraphFine);

        const Vector<VolIndex>& allCoarVofs = vofItIrregCoar[dit()].getVector();

        Vector<VoFStencil>& averageStencils = averageStencil[dit()];

        averageStencils.resize(allCoarVofs.size());

        for (int icoar = 0; icoar < allCoarVofs.size(); icoar++)
          {
            Vector<VolIndex> fineVofs;
            if (m_refRat > 2)
              {
                fineVofs = ebislStenCoar.refine(allCoarVofs[icoar], m_refRat, dit());
              }
            else
              {
                fineVofs = ebislStenCoar[dit()].refine(allCoarVofs[icoar]);
              }

            VoFStencil& averageStencil = averageStencils[icoar];

            for (int ifine = 0; ifine < fineVofs.size(); ifine++)
              {
                averageStencil.add(fineVofs[ifine], 1./numFinePerCoar);
              }
          }

        m_averageEBStencil[dit()] = RefCountedPtr<EBStencil>(new EBStencil(allCoarVofs, averageStencil[dit()], boxCoar,  boxFine, ebisBoxCoar,  ebisBoxFine,  m_ghost, m_ghost));
      }
  }
}
/************************************/
bool
EBMGAverage::isDefined() const
{
  return m_isDefined;
}
/************************************/
void
EBMGAverage::average(LevelData<EBCellFAB>&       a_coarData,
                     const LevelData<EBCellFAB>& a_fineData,
                     const Interval&             a_variables)
{
  CH_TIMERS("EBMGAverage::average");
  CH_TIMER("layout_changed_coarsenable", t1);
  CH_TIMER("layout_changed_not_coarsenable", t2);
  CH_TIMER("not_layout_changed", t3);
  CH_assert(a_fineData.ghostVect() == m_ghost);
  //CH_assert(a_coarData.ghostVect() == m_ghost);

  CH_assert(isDefined());

  if (m_layoutChanged)
    {
      if (m_coarsenable)
        {
          CH_START(t1);
          EBCellFactory ebcellfact(m_buffEBISL);
          LevelData<EBCellFAB> coarsenedFineData(m_buffGrids, m_nComp, m_ghost, ebcellfact);
          EBLevelDataOps::setVal(coarsenedFineData, 0.0);

          for (DataIterator dit = m_fineGrids.dataIterator(); dit.ok(); ++dit)
            {
              averageFAB(coarsenedFineData[dit()],
                         m_buffGrids[dit()],
                         a_fineData[dit()],
                         dit(),
                         a_variables);
            }
          coarsenedFineData.copyTo(a_variables, a_coarData, a_variables, m_copier);
          CH_STOP(t1);
        }
      else
        {
          CH_START(t2);
          EBCellFactory ebcellfact(m_buffEBISL);
          LevelData<EBCellFAB> refinedCoarseData(m_buffGrids, m_nComp, m_ghost, ebcellfact);
          EBLevelDataOps::setVal(refinedCoarseData, 0.0);

          a_fineData.copyTo(a_variables, refinedCoarseData, a_variables, m_copier);

          for (DataIterator dit = m_coarGrids.dataIterator(); dit.ok(); ++dit)
            {
              averageFAB(a_coarData[dit()],
                         m_coarGrids[dit()],
                         refinedCoarseData[dit()],
                         dit(),
                         a_variables);
            }
          CH_STOP(t2);
        }
    }
  else
    {
      CH_START(t3);
      averageMG(a_coarData, a_fineData, a_variables);
      CH_STOP(t3);
    }
}
/************************************/
void
EBMGAverage::averageMG(LevelData<EBCellFAB>&       a_coarData,
                       const LevelData<EBCellFAB>& a_fineData,
                       const Interval&             a_variables)
{
  CH_assert(a_fineData.ghostVect() == m_ghost);
  CH_assert(a_coarData.ghostVect() == m_ghost);

  CH_assert(isDefined());

  for (DataIterator dit = m_coarGrids.dataIterator();
      dit.ok(); ++dit)
    {
      averageFAB(a_coarData[dit()],
                 m_coarGrids[dit()],
                 a_fineData[dit()],
                 dit(),
                 a_variables);
    }
}
/************************************/
void
EBMGAverage::averageFAB(EBCellFAB&       a_coar,
                        const Box&       a_boxCoar,
                        const EBCellFAB& a_refCoar,
                        const DataIndex& a_datInd,
                        const Interval&  a_variables) const
{
  CH_TIMERS("EBMGAverage::average");
  CH_TIMER("regular_average", t1);
  CH_TIMER("irregular_average", t2);
  CH_assert(isDefined());

  const Box& coarBox = a_boxCoar;

  //do all cells as if they were regular
  Box refBox(IntVect::Zero, IntVect::Zero);
  refBox.refine(m_refRat);
  int numFinePerCoar = refBox.numPts();

  BaseFab<Real>& coarRegFAB =             a_coar.getSingleValuedFAB();
  const BaseFab<Real>& refCoarRegFAB = a_refCoar.getSingleValuedFAB();

  //set to zero because the fortran is a bit simpleminded
  //and does stuff additively
  a_coar.setVal(0.);
  CH_START(t1);
  for (int comp = a_variables.begin();  comp <= a_variables.end(); comp++)
    {
      FORT_REGAVERAGE(CHF_FRA1(coarRegFAB,comp),
                      CHF_CONST_FRA1(refCoarRegFAB,comp),
                      CHF_BOX(coarBox),
                      CHF_BOX(refBox),
                      CHF_CONST_INT(numFinePerCoar),
                      CHF_CONST_INT(m_refRat));
    }
  CH_STOP(t1);

  //this is really volume-weighted averaging even though it does
  //not look that way.

  //so (in the traditional sense) we want to preserve
  //rhoc * volc = sum(rhof * volf)
  //this translates to
  //volfrac_C * rhoC = (1/numFinePerCoar)(sum(volFrac_F * rhoF))
  //but the data input to this routine is all kappa weigthed so
  //the volumefractions have already been multiplied
  //which means
  // rhoC = (1/numFinePerCoar)(sum(rhoF))
  //which is what this does

  CH_START(t2);
  for (int comp = a_variables.begin();  comp <= a_variables.end(); comp++)
    {
      m_averageEBStencil[a_datInd]->apply(a_coar, a_refCoar, false, comp);
    }
  CH_STOP(t2);

}
/************************************/
#include "NamespaceFooter.H"
