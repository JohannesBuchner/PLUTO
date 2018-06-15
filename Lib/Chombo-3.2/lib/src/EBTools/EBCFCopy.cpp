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

#include "CH_Timer.H"
#include "DisjointBoxLayout.H"

#include "VoFIterator.H"
#include "EBCellFactory.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "EBLevelGrid.H"
#include "EBCFCopy.H"
#include "EBCFCopyF_F.H"

#include "NamespaceHeader.H"

void EBCFCopy::setDefaultValues()
{
  m_refRat    = -1;
  m_nComp     = -1;

  m_isDefined = false;
}

EBCFCopy::EBCFCopy()
{
  setDefaultValues();
}

EBCFCopy::~EBCFCopy()
{
}

EBCFCopy::EBCFCopy(const DisjointBoxLayout & a_dblFine,
                   const DisjointBoxLayout & a_dblCoar,
                   const EBISLayout        & a_ebislFine,
                   const EBISLayout        & a_ebislCoar,
                   const ProblemDomain     & a_domainCoar,
                   const int               & a_nref,
                   const int               & a_nvar,
                   const EBIndexSpace      * a_ebisPtr,
                   const IntVect           & a_ghost,
                   const bool              & a_layoutChanged)
{
  setDefaultValues();

  define(a_dblFine,
         a_dblCoar,
         a_ebislFine,
         a_ebislCoar,
         a_domainCoar,
         a_nref,
         a_nvar,
         a_ebisPtr,
         a_ghost,
         a_layoutChanged);
}

void EBCFCopy::define(const DisjointBoxLayout & a_dblFine,
                      const DisjointBoxLayout & a_dblCoar,
                      const EBISLayout        & a_ebislFine,
                      const EBISLayout        & a_ebislCoar,
                      const ProblemDomain     & a_domainCoar,
                      const int               & a_nref,
                      const int               & a_nvar,
                      const EBIndexSpace      * a_ebisPtr,
                      const IntVect           & a_ghostCellsPhi,
                      const bool              & a_layoutChanged)
{
  CH_TIMERS("EBCFCopy::define");

  CH_TIMER("fillEBISLayout",t1);

  m_ghost = a_ghostCellsPhi;
  m_nComp = a_nvar;

  m_coarGrids = a_dblCoar;
  m_fineGrids = a_dblFine;

  m_coarEBISL = a_ebislCoar;
  m_fineEBISL = a_ebislFine;

  m_coarDomain = a_domainCoar;

  m_refRat = a_nref;
  m_fineDomain = refine(m_coarDomain,m_refRat);

  m_layoutChanged = a_layoutChanged;

  m_coarsenable   = a_dblFine.coarsenable(m_refRat);

  // Only define ebislbuf and gridbuf if we are changing layouts
  if (m_layoutChanged)
  {
    ProblemDomain domebisl;
    if (m_coarsenable)
    {
      coarsen(m_buffGrids,m_fineGrids,m_refRat);
      domebisl = m_coarDomain;
    }
    else
    {
      refine(m_buffGrids,m_coarGrids,m_refRat);
      m_copierRCtoF.define(m_buffGrids,m_fineGrids,a_ghostCellsPhi);
      m_copierFtoRC.define(m_fineGrids,m_buffGrids,a_ghostCellsPhi);
      domebisl = m_fineDomain;
    }

    CH_START(t1);

    int nghost = 4;
    a_ebisPtr->fillEBISLayout(m_buffEBISL,
                              m_buffGrids,
                              domebisl,
                              nghost);

    if (m_refRat > 2)
    {
      if (m_coarsenable)
      {
        m_buffEBISL.setMaxRefinementRatio(m_refRat,a_ebisPtr);
      }
      else
      {
        m_buffEBISL.setMaxCoarseningRatio(m_refRat,a_ebisPtr);
      }
    }

    CH_STOP(t1);
  }

  defineStencils();

  m_isDefined = true;
}

void EBCFCopy::defineStencils()
{
  CH_TIME("EBCFCopy::defineStencils");

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

    LayoutData<Vector<VoFStencil> > copyStencilLayout;
    LayoutData<VoFIterator > vofItIrregFine;

    m_copyEBStencil.  define(gridsStenCoar);
    copyStencilLayout.define(gridsStenCoar);
    vofItIrregFine.   define(gridsStenCoar);

    for (DataIterator dit = gridsStenCoar.dataIterator(); dit.ok(); ++dit)
    {
      const Box&         boxFine =  gridsStenFine[dit()];
      const EBISBox& ebisBoxFine =  ebislStenFine[dit()];
      const EBGraph& ebGraphFine =  ebisBoxFine.getEBGraph();

      const Box&         boxCoar = gridsStenCoar[dit()];
      const EBISBox& ebisBoxCoar = ebislStenCoar[dit()];
      IntVectSet  notRegularCoar = ebisBoxCoar.getIrregIVS(boxCoar);

      IntVectSet ivsFine = refine(notRegularCoar,m_refRat);
      vofItIrregFine[dit()].define(ivsFine,ebGraphFine);

      const Vector<VolIndex>& allFineVofs = vofItIrregFine[dit()].getVector();

      Vector<VoFStencil>& copyStencils = copyStencilLayout[dit()];
      copyStencils.resize(allFineVofs.size());

      for (int ifine = 0; ifine < allFineVofs.size(); ifine++)
      {
        VolIndex coarseVof;

        if (m_refRat > 2)
        {
          coarseVof = ebislStenFine.coarsen(allFineVofs[ifine],m_refRat,dit());
        }
        else
        {
          coarseVof = ebisBoxFine.coarsen(allFineVofs[ifine]);
        }

        VoFStencil& copyStencil = copyStencils[ifine];

        // weight is one because there is only one vof to the stencil
        // and we are doing piecewise constant interpolation
        // for the copy stencil the weight is one because it is a copy
        copyStencil.add(coarseVof,1.0);
      }

      m_copyEBStencil[dit()] = RefCountedPtr<EBStencil>
                                  (new EBStencil(allFineVofs,
                                                 copyStencilLayout[dit()],
                                                 boxFine,
                                                 boxCoar,
                                                 ebisBoxFine,
                                                 ebisBoxCoar,
                                                 m_ghost,
                                                 m_ghost));
    }
  }
}

bool EBCFCopy::isDefined() const
{
  return m_isDefined;
}

void EBCFCopy::copy(LevelData<EBCellFAB>       & a_fineData,
                    const LevelData<EBCellFAB> & a_coarData,
                    const Interval             & a_variables)
{
  CH_TIME("EBCFCopy::copy");

  CH_assert(a_fineData.ghostVect() == m_ghost);
  CH_assert(a_coarData.ghostVect() == m_ghost);

  if (m_layoutChanged)
  {
    if (m_coarsenable)
    {
      CH_TIME("EBCFCopy::copy::coarsenable");

      EBCellFactory ebcellfact(m_buffEBISL);
      LevelData<EBCellFAB> coarsenedFineData(m_buffGrids,m_nComp,m_ghost,ebcellfact);
      a_coarData.copyTo(a_variables,coarsenedFineData,a_variables);

      for (DataIterator dit = m_fineGrids.dataIterator(); dit.ok(); ++dit)
      {
        // does incrementonly = true
        copyFAB(a_fineData[dit()],
                m_buffGrids[dit()],
                coarsenedFineData[dit()],
                dit(),
                a_variables);
      }
    }
    else
    {
      CH_TIME("EBCFCopy::copy::uncoarsenable");

      EBCellFactory ebcellfact(m_buffEBISL);
      LevelData<EBCellFAB> refinedCoarseData(m_buffGrids,m_nComp,m_ghost,ebcellfact);

      for (DataIterator dit = m_coarGrids.dataIterator(); dit.ok(); ++dit)
      {
        copyFAB(refinedCoarseData[dit()],
                     m_coarGrids[dit()],
                     a_coarData[dit()],
                     dit(),
                     a_variables);
      }

      refinedCoarseData.copyTo(a_variables,a_fineData,a_variables,m_copierRCtoF);
    }
  }
  else
  {
    copySameLayout(a_fineData,a_coarData,a_variables);
  }
}

void EBCFCopy::copySameLayout(LevelData<EBCellFAB>       & a_fineData,
                              const LevelData<EBCellFAB> & a_coarData,
                              const Interval             & a_variables)
{
  CH_TIME("EBCFCopy::copySameLayout");

  CH_assert(a_fineData.ghostVect() == m_ghost);
  CH_assert(a_coarData.ghostVect() == m_ghost);

  for (DataIterator dit = m_coarGrids.dataIterator(); dit.ok(); ++dit)
  {
    copyFAB(a_fineData[dit()],
            m_coarGrids[dit()],
            a_coarData[dit()],
            dit(),
            a_variables);
  }
}

void EBCFCopy::copyFAB(EBCellFAB       & a_refCoar,
                       const Box       & a_coarBox,
                       const EBCellFAB & a_coar,
                       const DataIndex & a_datInd,
                       const Interval  & a_variables) const
{
  CH_TIMERS("EBCFCopy::copyFAB");
  CH_TIMER("regular_copy",t1);
  CH_TIMER("irregular_copy",t2);

  CH_assert(isDefined());

  const Box& coarBox = a_coarBox;

  for (int ivar = a_variables.begin();  ivar <= a_variables.end(); ivar++)
  {
    m_copyEBStencil[a_datInd]->cache(a_refCoar,ivar);

    // do all cells as if they were regular
    Box refBox(IntVect::Zero,IntVect::Zero);
    refBox.refine(m_refRat);

    const BaseFab<Real>& coarRegFAB = a_coar.   getSingleValuedFAB();
    BaseFab<Real>& refCoarRegFAB    = a_refCoar.getSingleValuedFAB();

    CH_START(t1);

    FORT_COPYCFFAB(CHF_FRA1(refCoarRegFAB,ivar),
                   CHF_CONST_FRA1(coarRegFAB,ivar),
                   CHF_BOX(coarBox),
                   CHF_BOX(refBox),
                   CHF_CONST_INT(m_refRat));

    CH_STOP(t1);

    m_copyEBStencil[a_datInd]->uncache(a_refCoar,ivar);

    CH_START(t2);
    m_copyEBStencil[a_datInd]->apply(a_refCoar,a_coar,false,ivar);
    CH_STOP(t2);
  }
}

#include "NamespaceFooter.H"
