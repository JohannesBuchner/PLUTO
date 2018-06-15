#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBCoarToCoarRedist.H"
#include "LayoutIterator.H"
#include "BaseIVFactory.H"
#include "VoFIterator.H"
#include "EBCellFactory.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"
/***********************/
/***********************/
void
EBCoarToCoarRedist::
resetWeights(const LevelData<EBCellFAB>& a_modifier,
             const int& a_ivar)
{
  CH_TIME("EBCoarToCoarRedist::resetWeights");
  CH_assert(isDefined());
  Interval srcInterv(a_ivar, a_ivar);
  Interval dstInterv(0,0);
  a_modifier.copyTo(srcInterv, m_densityCoar, dstInterv);
  //set the weights to mass weighting if the modifier
  //is the coarse density.
  for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit)
    {
      const IntVectSet& coarSet = m_setsCoar[dit()];
      BaseIVFAB<VoFStencil>& massStenFAB = m_stenCoar[dit()];
      const BaseIVFAB<VoFStencil>& volStenFAB  = m_volumeStenc[dit()];
      const BaseIVFAB<VoFStencil>& stanStenFAB  = m_standardStenc[dit()];

      const EBISBox& ebisBox = m_ebislCoar[dit()];
      const EBCellFAB& modFAB =    m_densityCoar[dit()];

      for (VoFIterator vofit(coarSet, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& thisVoF = vofit();

          VoFStencil stanSten = stanStenFAB(thisVoF, 0);
          VoFStencil oldSten =   volStenFAB(thisVoF, 0);
          VoFStencil newSten;
          for (int isten = 0; isten < oldSten.size(); isten++)
            {
              const VolIndex& thatVoF = oldSten.vof(isten);
              Real weight =modFAB(thatVoF, a_ivar);
              newSten.add(thatVoF, weight);
            }
          //need normalize by the whole standard stencil
          Real sum = 0.0;
          for (int isten = 0; isten < stanSten.size(); isten++)
            {
              const VolIndex& thatVoF = stanSten.vof(isten);
              Real weight = modFAB(thatVoF, 0);
              Real volfrac = ebisBox.volFrac(thatVoF);
              //it is weight*volfrac that is normalized
              sum += weight*volfrac;
            }
          if (Abs(sum) > 0.0)
            {
              Real scaling = 1.0/sum;
              newSten *= scaling;
            }
          massStenFAB(thisVoF, 0) = newSten;
        }
    }
}
/**********************/
EBCoarToCoarRedist::
EBCoarToCoarRedist()
{
  m_isDefined = false;
}
/**********************/
EBCoarToCoarRedist::
~EBCoarToCoarRedist()
{
}
/**********************/
void
EBCoarToCoarRedist::
define(const DisjointBoxLayout& a_dblFine,
       const DisjointBoxLayout& a_dblCoar,
       const EBISLayout& a_ebislFine,
       const EBISLayout& a_ebislCoar,
       const Box& a_domainCoar,
       const int& a_nref,
       const int& a_nvar,
       int a_redistRad)
{
  CH_TIME("EBCoarToCoarRedist::define");
  m_isDefined = true;
  m_nComp = a_nvar;
  m_refRat = a_nref;
  m_domainCoar = a_domainCoar;
  m_gridsFine = a_dblFine;
  m_gridsCoar = a_dblCoar;
  m_ebislFine = a_ebislFine;
  m_ebislCoar = a_ebislCoar;
  m_redistRad = a_redistRad;
  //define the intvectsets over which the objects live
  m_setsCoar.define(m_gridsCoar);
  m_stenCoar.define(m_gridsCoar);
  m_volumeStenc.define(m_gridsCoar);
  m_standardStenc.define(m_gridsCoar);

  //make sets
  for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit)
    {
      Box gridBox = m_gridsCoar.get(dit());
      gridBox.grow(a_redistRad);
      gridBox &= m_domainCoar;
      //find the complement of what we really want
      IntVectSet ivsComplement(gridBox);
      for (LayoutIterator litFine =
            m_gridsFine.layoutIterator(); litFine.ok(); ++litFine)
        {
          Box projBox = coarsen(m_gridsFine.get(litFine()), m_refRat);
          ivsComplement -= projBox;
        }
      //now the set we want is the gridbox - complement
      IntVectSet& coarSet = m_setsCoar[dit()];
      coarSet = IntVectSet(gridBox);
      coarSet -= ivsComplement;
      IntVectSet irregIVS = m_ebislCoar[dit()].getIrregIVS(gridBox);
      coarSet &= irregIVS;
    }
  defineDataHolders();
  //initialize the buffers to zero
  setToZero();
}
/**********************/
void
EBCoarToCoarRedist::
define(const EBLevelGrid& a_eblgFine,
       const EBLevelGrid& a_eblgCoar,
       const int&         a_nref,
       const int&         a_nvar,
       const int&         a_redistRad)
{
  CH_TIME("EBCoarToCoarRedist::define");
  m_isDefined = true;
  m_nComp     =  a_nvar;
  m_refRat    =  a_nref;
  m_redistRad =  a_redistRad;
  m_domainCoar=  a_eblgCoar.getDomain().domainBox();
  m_gridsFine =  a_eblgFine.getDBL();
  m_gridsCoar =  a_eblgCoar.getDBL();
  m_ebislFine =  a_eblgFine.getEBISL();
  m_ebislCoar =  a_eblgCoar.getEBISL();
  //define the intvectsets over which the objects live

  //make sets.
  m_setsCoar.define(m_gridsCoar);
  //need the set of points under the fine grid that
  //will redistribute to this box
  IntVectSet coveringIVS = a_eblgFine.getCoveringIVS();
  coveringIVS.coarsen(m_refRat);
  for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit)
    {
      Box gridBox = m_gridsCoar.get(dit());
      gridBox.grow(a_redistRad);
      gridBox &= m_domainCoar;
      m_setsCoar[dit()]  = m_ebislCoar[dit()].getIrregIVS(gridBox);
      m_setsCoar[dit()] &= coveringIVS;
    }
  defineDataHolders();
  //initialize the buffers to zero
  setToZero();
}
/***/
void
EBCoarToCoarRedist::
defineDataHolders()
{
  CH_TIME("EBCoarToCoarRedist::defineDataHolders");
  EBCellFactory ebcellfactcoar(m_ebislCoar);
  m_densityCoar.define(m_gridsCoar, 1, 2*m_redistRad*IntVect::Unit, ebcellfactcoar);
  m_stenCoar.define(m_gridsCoar);
  m_volumeStenc.define(m_gridsCoar);
  m_standardStenc.define(m_gridsCoar);

  RedistStencil rdStenCoar(m_gridsCoar, m_ebislCoar, m_domainCoar, m_redistRad);
  for (DataIterator dit =
        m_gridsCoar.dataIterator(); dit.ok(); ++dit)
    {
      BaseIVFAB<VoFStencil>& stenFAB = m_stenCoar[dit()];
      BaseIVFAB<VoFStencil>& volStenFAB = m_volumeStenc[dit()];
      BaseIVFAB<VoFStencil>& stanStenFAB = m_standardStenc[dit()];

      const EBISBox& ebisBox = m_ebislCoar[dit()];
      const IntVectSet& ivs  = m_setsCoar[dit()];
      const BaseIVFAB<VoFStencil>& rdStenFAB = rdStenCoar[dit()];
      const Box& gridCoar = m_gridsCoar.get(dit());
      stenFAB.define(    ivs, ebisBox.getEBGraph(), 1);
      volStenFAB.define( ivs, ebisBox.getEBGraph(), 1);
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
              if (gridCoar.contains(dstVoF.gridIndex()))
                {
                  newStencil.add(dstVoF, weight);
                }
            }
          //set both to volume weighted.  can be changed later
          stenFAB(srcVoF, 0)     = newStencil;
          volStenFAB(srcVoF, 0)  = newStencil;
          stanStenFAB(srcVoF, 0) = stanSten;

        }
    }

  BaseIVFactory<Real> factCoar(m_ebislCoar, m_setsCoar);
  IntVect ivghost = m_redistRad*IntVect::Unit;
  m_regsCoar.define(m_gridsCoar, m_nComp, ivghost, factCoar);

}
/**********************/
void
EBCoarToCoarRedist::
setToZero()
{
  CH_TIME("EBCoarToCoarRedist::setToZero");
  for (DataIterator dit = m_gridsCoar.dataIterator();
      dit.ok(); ++dit)
    {
      m_regsCoar[dit()].setVal(0.0);
    }
}
/**********************/
void
EBCoarToCoarRedist::
increment(const BaseIVFAB<Real>& a_coarMass,
          const DataIndex& a_coarDataIndex,
          const Interval&  a_variables)
{
  CH_TIME("EBCoarToCoarRedist::increment");
  BaseIVFAB<Real>& coarBuf =  m_regsCoar[a_coarDataIndex];
  const EBISBox&   ebisBox = m_ebislCoar[a_coarDataIndex];

  const Box& gridBox = m_gridsCoar.get(a_coarDataIndex);
  IntVectSet ivs = ebisBox.getIrregIVS(gridBox);
  const IntVectSet& fabIVS = a_coarMass.getIVS();
  const IntVectSet& bufIVS = m_setsCoar[a_coarDataIndex];
  CH_assert(fabIVS.contains(ivs));
  //only a subset of our points are in the set (only the ones covered by finer grids).
  ivs &= bufIVS;
  for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for (int ivar = a_variables.begin();
          ivar <= a_variables.end();  ivar++)
        {
          coarBuf(vof, ivar) += a_coarMass(vof, ivar);
        }
    }
}
/**********************/
void
EBCoarToCoarRedist::
redistribute(LevelData<EBCellFAB>& a_coarSolution,
             const Interval& a_variables)
{
  CH_TIME("EBCoarToCoarRedist::redistribute");
  m_regsCoar.exchange(a_variables);
  //at all points in the buffer, subtract off the redistributed values
  for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit)
    {
      const BaseIVFAB<Real>& regCoar = m_regsCoar[dit()];
      const IntVectSet& ivsCoar = m_setsCoar[dit()];
      const EBISBox& ebisBoxCoar = m_ebislCoar[dit()];
      const BaseIVFAB<VoFStencil>& stenFAB = m_stenCoar[dit()];
      const Box& coarBox = m_gridsCoar.get(dit());

      EBCellFAB& solFAB = a_coarSolution[dit()];

      for (VoFIterator vofit(ivsCoar,ebisBoxCoar.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& srcVoF = vofit();
          const VoFStencil& vofsten = stenFAB(srcVoF, 0);
          for (int isten = 0; isten < vofsten.size(); isten++)
            {
              const Real& weight = vofsten.weight(isten);
              const VolIndex& dstVoF = vofsten.vof(isten);
              //since we are just using the input redist stencil,
              //have to check if out of bounds
              if (coarBox.contains(dstVoF.gridIndex()))
                {
                  for (int ivar = a_variables.begin();
                      ivar <= a_variables.end();  ivar++)
                    {
                      Real dmFine = regCoar(srcVoF, ivar);
                      //ucoar+= massfine/volcoar, ie.
                      //ucoar+= (wcoar*dmCoar*volFracfine/volfraccoar)=massfine/volcoar
                      Real dUCoar = dmFine*weight;
                      //SUBTRACTING because this is re-redistribution.
                      solFAB(dstVoF, ivar) -= dUCoar;
                    }
                }
            }
        }
    }
}
/**********************/
bool
EBCoarToCoarRedist::
isDefined() const
{
  return m_isDefined;
}
/**********************/
#include "NamespaceFooter.H"
