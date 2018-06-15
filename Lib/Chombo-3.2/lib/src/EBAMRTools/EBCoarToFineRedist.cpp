#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBCoarToFineRedist.H"
#include "LayoutIterator.H"
#include "BaseIVFactory.H"
#include "VoFIterator.H"
#include "RedistStencil.H"
#include "EBIndexSpace.H"
#include "CH_Timer.H"
#include "EBCellFactory.H"
#include "NamespaceHeader.H"
/**********************/
EBCoarToFineRedist::
EBCoarToFineRedist()
{
  m_isDefined = false;
}
/**********************/
EBCoarToFineRedist::
~EBCoarToFineRedist()
{
}
/***********************/
void
EBCoarToFineRedist::
resetWeights(const LevelData<EBCellFAB>& a_modifierCoar,
             const int& a_ivar)
{
  CH_TIME("EBCoarToFineRedist::resetWeights");
  Interval srcInterv(a_ivar, a_ivar);
  Interval dstInterv(0,0);
  a_modifierCoar.copyTo(srcInterv, m_densityCedFine, dstInterv);

  //set the weights to mass weighting if the modifier
  //is the fine density.
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      const IntVectSet& coarSet = m_setsCedFine[dit()];
      BaseIVFAB<VoFStencil>& massStenFAB = m_stenCedFine[dit()];
      const BaseIVFAB<VoFStencil>& volStenFAB  = m_volumeStenc[dit()];
      const BaseIVFAB<VoFStencil>& stanStenFAB  = m_standardStenc[dit()];

      const EBISBox& ebisBoxCoar = m_ebislCedFine[dit()];
      const EBCellFAB& modFAB = m_densityCedFine[dit()];

      for (VoFIterator vofit(coarSet, ebisBoxCoar.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vofCoar = vofit();

          const VoFStencil& stanSten = stanStenFAB(vofCoar, 0);
          const VoFStencil& volSten  = volStenFAB(vofCoar, 0);
          VoFStencil newSten;
          for (int isten = 0; isten < volSten.size(); isten++)
            {
              const VolIndex& thatVoFCoar = volSten.vof(isten);
              Real weight = modFAB(thatVoFCoar, 0);
              //average fine weights onto coarse weight
              newSten.add(thatVoFCoar, weight);
            }
          //need normalize by the whole standard stencil
          Real sum = 0.0;
          for (int isten = 0; isten < stanSten.size(); isten++)
            {
              const VolIndex& thatVoFCoar = stanSten.vof(isten);
              Real weight = modFAB(thatVoFCoar, 0);
              Real volfrac = ebisBoxCoar.volFrac(thatVoFCoar);
              //it is weight*volfrac that is normalized
              sum += weight*volfrac;
            }
          if (Abs(sum) > 0.0)
            {
              Real scaling = 1.0/sum;
              newSten *= scaling;
            }
          massStenFAB(vofCoar, 0) = newSten;
        }
    }
}
/**********************/
void
EBCoarToFineRedist::
define(const EBLevelGrid&  a_eblgFine,
       const EBLevelGrid&  a_eblgCoar,
       const int&          a_nref,
       const int&          a_nvar,
       const int&          a_redistRad)
{
  CH_TIME("EBCoarToFineRedist::define");

  //from here we can assume the redistRad == 1
  m_isDefined  = true;
  m_nComp      = a_nvar;
  m_refRat     = a_nref;
  m_domainCoar = a_eblgCoar.getDomain().domainBox();
  m_gridsFine  = a_eblgFine.getDBL();
  m_gridsCoar  = a_eblgCoar.getDBL();
  m_ebislFine  = a_eblgFine.getEBISL();
  m_ebislCoar  = a_eblgCoar.getEBISL();
  m_redistRad  = a_redistRad;

  //created the coarsened fine layout
  m_gridsCedFine = DisjointBoxLayout();
  coarsen(m_gridsCedFine, m_gridsFine, m_refRat);
  const EBIndexSpace* const ebisPtr = a_eblgFine.getEBIS();
  CH_assert(ebisPtr->isDefined());
  int nghost = 3*m_redistRad;
  ebisPtr->fillEBISLayout(m_ebislCedFine, m_gridsCedFine,
                          m_domainCoar, nghost);
  m_ebislCedFine.setMaxRefinementRatio(m_refRat, ebisPtr);

  //define the intvectsets over which the objects live
  m_setsCoar.define(m_gridsCoar);
  m_setsCedFine.define(m_gridsCedFine);
  //the buffers for the fine level are really simple.
  //grow the coarsened fine grid by the redist radius
  //and subtract off the coarsened fine grids.
  //intersect with the irregular cells.
  //this way get irregular cells range of this grid on coarse
  //grid.
  IntVectSet coveringIVS;
  {
    CH_TIME("coarsened_fine_sets");
    coveringIVS = a_eblgFine.getCoveringIVS();
    coveringIVS.coarsen(m_refRat);
    for (DataIterator dit = m_gridsCedFine.dataIterator();
        dit.ok(); ++dit)
      {
        Box grownBox = grow(m_gridsCedFine.get(dit()), m_redistRad);
        grownBox &= m_domainCoar;
        m_setsCedFine[dit()] = m_ebislCedFine[dit()].getIrregIVS(grownBox);
        m_setsCedFine[dit()] -= coveringIVS;
      }
  }

  //the coarse buffers are more complicated.
  //Need to intersect irregular cells on this grid
  // with entire coarse/fine interface of width redist radius
  {
    CH_TIME("coarse_sets");
    IntVectSet cfInterface = coveringIVS;
    cfInterface.grow(m_redistRad);
    cfInterface &= m_domainCoar;
    cfInterface -= coveringIVS;
    for (DataIterator dit = m_gridsCoar.dataIterator();
        dit.ok(); ++dit)
      {
        Box coarBox = m_gridsCoar.get(dit());
        m_setsCoar[dit()] =  m_ebislCoar[dit()].getIrregIVS(coarBox);
        m_setsCoar[dit()] &= cfInterface;
      }
  }
  defineDataHolders();
  setToZero();
}
/************/
void
EBCoarToFineRedist::
defineDataHolders()
{
  CH_TIME("EBCoarToFineRedist::defineDataHolders");
  //make mass buffers
  BaseIVFactory<Real> factCoar(m_ebislCoar, m_setsCoar);
  IntVect noghost = IntVect::Zero;
  m_regsCoar.define(m_gridsCoar, m_nComp, noghost, factCoar);

  BaseIVFactory<Real> factCedFine(m_ebislCedFine, m_setsCedFine);
  m_regsCedFine.define(m_gridsCedFine, m_nComp, m_redistRad*IntVect::Unit, factCedFine);

  EBCellFactory ebcellfact(m_ebislCedFine);
  m_densityCedFine.define(m_gridsCedFine, 1, 2*m_redistRad*IntVect::Unit, ebcellfact);

  //define the stencils with volume weights.
  //if you want mass weights or whatever, use reset weights
  m_stenCedFine.define(m_gridsCedFine);
  m_volumeStenc.define(m_gridsCedFine);
  m_standardStenc.define(m_gridsCedFine);
  RedistStencil stenCedFine(m_gridsCedFine, m_ebislCedFine,
                            m_domainCoar, m_redistRad);

  for (DataIterator ditCedF = m_gridsCedFine.dataIterator();
      ditCedF.ok(); ++ditCedF)
    {
      BaseIVFAB<VoFStencil>& stenFAB = m_stenCedFine[ditCedF()];
      BaseIVFAB<VoFStencil>& volStenFAB = m_volumeStenc[ditCedF()];
      BaseIVFAB<VoFStencil>& stanStenFAB = m_standardStenc[ditCedF()];
      const EBISBox& ebisBox = m_ebislCedFine[ditCedF()];
      const IntVectSet& ivs  = m_setsCedFine[ditCedF()];
      const BaseIVFAB<VoFStencil>& rdStenFAB = stenCedFine[ditCedF()];
      const Box& gridCedFine = m_gridsCedFine.get(ditCedF());
      stenFAB.define(ivs, ebisBox.getEBGraph(), 1);
      volStenFAB.define(ivs, ebisBox.getEBGraph(), 1);
      stanStenFAB.define(ivs, ebisBox.getEBGraph(), 1);

      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& srcVoF = vofit();
          VoFStencil newStencil;
          const VoFStencil& stanSten = rdStenFAB(srcVoF, 0);
          for (int istan = 0; istan < stanSten.size(); istan++)
            {
              const VolIndex& dstVoF = stanSten.vof(istan);
              const Real& weight = stanSten.weight(istan);
              if (gridCedFine.contains(dstVoF.gridIndex()))
                {
                  newStencil.add(dstVoF, weight);
                }
            }
          stenFAB(srcVoF, 0)    = newStencil;
          //need to keep these around to get mass redist right
          volStenFAB(srcVoF, 0) = newStencil;
          stanStenFAB(srcVoF,0) = stanSten;

        }
    }
}
/**********************/
void
EBCoarToFineRedist::
define(const DisjointBoxLayout& a_dblFine,
       const DisjointBoxLayout& a_dblCoar,
       const EBISLayout& a_ebislFine,
       const EBISLayout& a_ebislCoar,
       const Box& a_domainCoar,
       const int& a_nref,
       const int& a_nvar,
       int a_redistRad,
       const EBIndexSpace* ebisPtr)
{
  CH_TIME("EBCoarToFineRedist::define");
  m_isDefined = true;
  m_nComp = a_nvar;
  m_refRat = a_nref;
  m_domainCoar = a_domainCoar;
  m_gridsFine = a_dblFine;
  m_gridsCoar = a_dblCoar;
  m_ebislFine = a_ebislFine;
  m_ebislCoar = a_ebislCoar;
  m_redistRad = a_redistRad;

  //created the coarsened fine layout
  m_gridsCedFine = DisjointBoxLayout();
  coarsen(m_gridsCedFine, m_gridsFine, m_refRat);

  CH_assert(ebisPtr->isDefined());
  int nghost = 3*m_redistRad;
  ebisPtr->fillEBISLayout(m_ebislCedFine, m_gridsCedFine,
                          m_domainCoar, nghost);
  m_ebislCedFine.setMaxRefinementRatio(m_refRat, ebisPtr);

  //define the intvectsets over which the objects live
  m_setsCoar.define(m_gridsCoar);
  m_setsCedFine.define(m_gridsCedFine);
  //the buffers for the fine level are really simple.
  //grow the coarsened fine grid by the redist radius
  //and subtract off the coarsened fine grids.
  //intersect with the irregular cells.
  //this way get irregular cells range of this grid on coarse
  //grid.
  {
    CH_TIME("coarsened_fine_sets");
    for (DataIterator dit = m_gridsCedFine.dataIterator();
        dit.ok(); ++dit)
      {
        Box grownBox = grow(m_gridsCedFine.get(dit()), m_redistRad);
        grownBox &= m_domainCoar;
        IntVectSet localIVS = m_ebislCedFine[dit()].getIrregIVS(grownBox);
        //subtract off all boxes on fine level so we get stuff at the
        //coarse-fine interface that will redistribute to this grid
        for (LayoutIterator lit = m_gridsCedFine.layoutIterator();
            lit.ok(); ++lit)
          {
            localIVS -= m_gridsCedFine.get(lit());
          }
        m_setsCedFine[dit()] = localIVS;
      }
  }
  //the coarse buffers are more complicated.
  //Need to intersect irregular cells on this grid
  // with entire coarse/fine interface of width redist radius
  {
    CH_TIME("coarse_sets");
    for (DataIterator dit = m_gridsCoar.dataIterator();
        dit.ok(); ++dit)
      {
        Box coarBox = m_gridsCoar.get(dit());
        //make the complement set the whole
        IntVectSet cfIntComp(coarBox);
        //subtract the CF Interface from each coarsened fine box
        //from the complement
        for (LayoutIterator lit1 = m_gridsCedFine.layoutIterator();
            lit1.ok(); ++lit1)
          {
            Box grownBox = grow(m_gridsCedFine.get(lit1()), m_redistRad);
            grownBox &= m_domainCoar;
            IntVectSet cedCFIVS(grownBox);
            //subtract off all boxes on fine level so we get stuff at the
            //coarse-fine interface that will redistribute to this grid
            for (LayoutIterator lit2 = m_gridsCedFine.layoutIterator();
                lit2.ok(); ++lit2)
              {
                cedCFIVS -= m_gridsCedFine.get(lit2());
              }

            cfIntComp -=cedCFIVS;
          }
        //coarse fine interface is whole box - complement
        IntVectSet cfInterface(coarBox);
        cfInterface -= cfIntComp;

        IntVectSet localIVS = m_ebislCoar[dit()].getIrregIVS(coarBox);
        //intersect with fat coarse-fine interface
        localIVS &= cfInterface;
        m_setsCoar[dit()] = localIVS;
      }
  }

  //define the registers and all that
  defineDataHolders();

  //initialize the buffers to zero
  setToZero();
}
/**********************/
void
EBCoarToFineRedist::
setToZero()
{
  CH_TIME("EBCoarToFineRedist::setToZero");
  for (DataIterator dit = m_gridsCedFine.dataIterator();
      dit.ok(); ++dit)
    {
      m_regsCedFine[dit()].setVal(0.0);
    }

  for (DataIterator dit = m_gridsCoar.dataIterator();
      dit.ok(); ++dit)
    {
      m_regsCoar[dit()].setVal(0.0);
    }
}
/**********************/
void
EBCoarToFineRedist::increment(const BaseIVFAB<Real>& a_coarseMass,
                              const DataIndex& a_coarseDataIndex,
                              const Interval&  a_variables)
{
  BaseIVFAB<Real>& coarBuf =  m_regsCoar[a_coarseDataIndex];
  const EBISBox&   ebisBox = m_ebislCoar[a_coarseDataIndex];
  const IntVectSet& ivsLoc =  m_setsCoar[a_coarseDataIndex];
  CH_assert(a_coarseMass.getIVS().contains(ivsLoc));
  CH_assert(coarBuf.getIVS().contains(ivsLoc));
  for (VoFIterator vofit(ivsLoc, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for (int ivar = a_variables.begin();
          ivar <= a_variables.end();  ivar++)
        {
          coarBuf(vof, ivar) += a_coarseMass(vof, ivar);
        }
    }
}
/**********************/
void
EBCoarToFineRedist::
redistribute(LevelData<EBCellFAB>& a_fineSolution,
             const Interval& a_variables)
{
  CH_TIME("EBCoarToFineRedist::redistribute");
  //copy the buffer to the fine layout
  m_regsCoar.copyTo(a_variables, m_regsCedFine, a_variables);
  //redistribute the coarsened fine registers to the fine solution
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      const BaseIVFAB<Real>& regCoar = m_regsCedFine[dit()];
      const IntVectSet& ivsCoar = m_setsCedFine[dit()];
      const EBISBox& ebisBoxCoar = m_ebislCedFine[dit()];
      const BaseIVFAB<VoFStencil>& stenFAB = m_stenCedFine[dit()];

      EBCellFAB& solFAB = a_fineSolution[dit()];

      for (VoFIterator vofitCoar(ivsCoar, ebisBoxCoar.getEBGraph());
          vofitCoar.ok(); ++vofitCoar)
        {
          const VolIndex& srcVoFCoar = vofitCoar();
          const VoFStencil& vofsten = stenFAB(srcVoFCoar, 0);
          for (int isten = 0; isten < vofsten.size(); isten++)
            {
              const Real& weight = vofsten.weight(isten);
              const VolIndex& dstVoFCoar = vofsten.vof(isten);
              Vector<VolIndex> vofsFine =
                m_ebislCedFine.refine(dstVoFCoar,m_refRat, dit());

              for (int ivar = a_variables.begin();
                  ivar <= a_variables.end();  ivar++)
                {
                  Real dmCoar = regCoar(srcVoFCoar, ivar);
                  for (int ifine = 0; ifine < vofsFine.size(); ifine++)
                    {

                      const VolIndex& dstVoFFine = vofsFine[ifine];
                      //ufine += (wcoar*dmCoar) (piecewise constant density diff)
                      Real dUFine = dmCoar*weight;
                      solFAB(dstVoFFine, ivar) += dUFine;
                    }
                }
            }
        }
    }
}
/**********************/
bool
EBCoarToFineRedist::
isDefined() const
{
  return m_isDefined;
}
/**********************/
#include "NamespaceFooter.H"
