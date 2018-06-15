#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves  11-02-2001

#include "EBFineToCoarRedist.H"
#include "LayoutIterator.H"
#include "BaseIVFactory.H"
#include "VoFIterator.H"
#include "EBIndexSpace.H"
#include "EBCellFactory.H"
#include "EBArith.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"
/***********************/
/***********************/
void
EBFineToCoarRedist::
resetWeights(const LevelData<EBCellFAB>& a_modifierCoar,
             const int& a_ivar)
{
  CH_TIME("EBFineToCoarRedist::resetWeights");
  //set the weights to mass weighting if the modifier
  //is the coarse density.
  Interval srcInterv(a_ivar, a_ivar);
  Interval dstInterv(0,0);
  a_modifierCoar.copyTo(srcInterv, m_densityCoar, dstInterv);
  for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit)
    {
      const IntVectSet& fineSet = m_setsRefCoar[dit()];
      BaseIVFAB<VoFStencil>& massStenFAB = m_stenRefCoar[dit()];
      const BaseIVFAB<VoFStencil>& volStenFAB  = m_volumeStenc[dit()];
      const BaseIVFAB<VoFStencil>& stanStenFAB  = m_standardStenc[dit()];
      const EBISBox& ebisBoxFine = m_ebislRefCoar[dit()];
      const EBCellFAB& modFAB = m_densityCoar[dit()];

      for (VoFIterator vofit(fineSet, ebisBoxFine.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vofFine = vofit();

          const VoFStencil& oldSten = volStenFAB(vofFine, 0);
          const VoFStencil& stanSten = stanStenFAB(vofFine, 0);

          VoFStencil newSten;
          for (int isten = 0; isten < oldSten.size(); isten++)
            {
              const VolIndex& thatVoFFine = oldSten.vof(isten);
              VolIndex refCoarVoF =
                m_ebislRefCoar.coarsen(thatVoFFine, m_refRat, dit());

              Real weight = modFAB(refCoarVoF, a_ivar);

              //it is weight*volfrac that is normalized
              newSten.add(thatVoFFine, weight);
            }
          //have to normalize using the whole stencil
          Real sum = 0.0;
          for (int isten = 0; isten < stanSten.size(); isten++)
            {
              const VolIndex& thatVoFFine = stanSten.vof(isten);
              VolIndex refCoarVoF =
                m_ebislRefCoar.coarsen(thatVoFFine, m_refRat, dit());

              Real weight = modFAB(refCoarVoF, a_ivar);

              Real volfrac = ebisBoxFine.volFrac(thatVoFFine);
              //it is weight*volfrac that is normalized
              sum += weight*volfrac;
            }
          if (Abs(sum) > 0.0)
            {
              Real scaling = 1.0/sum;
              newSten *= scaling;
            }
          massStenFAB(vofFine, 0) = newSten;

        }
    }
}
/**********************/
EBFineToCoarRedist::EBFineToCoarRedist()
{
  m_isDefined = false;
}
/**********************/
EBFineToCoarRedist::~EBFineToCoarRedist()
{
}
/**********************/
void
EBFineToCoarRedist::
define(const EBLevelGrid& a_eblgFine,
       const EBLevelGrid& a_eblgCoar,
       const int&         a_nref,
       const int&         a_nvar,
       const int&         a_redistRad)
{
  CH_TIME("EBFineToCoarRedist::define_with_eblevelgrids");
  m_isDefined = true;
  m_nComp     = a_nvar;
  m_refRat    = a_nref;
  m_redistRad = a_redistRad;
  m_domainCoar= a_eblgCoar.getDomain().domainBox();
  m_gridsFine = a_eblgFine.getDBL();
  m_gridsCoar = a_eblgCoar.getDBL();
  m_ebislFine = a_eblgFine.getEBISL();
  m_ebislCoar = a_eblgCoar.getEBISL();
  //created the coarsened fine layout
  m_gridsRefCoar = DisjointBoxLayout();
  refine(m_gridsRefCoar, m_gridsCoar, m_refRat);

  const EBIndexSpace* ebisPtr = a_eblgFine.getEBIS();
  CH_assert(ebisPtr->isDefined());
  int nghost = 3*m_redistRad;
  Box domainFine = refine(m_domainCoar, m_refRat);
  ebisPtr->fillEBISLayout(m_ebislRefCoar, m_gridsRefCoar,
                          domainFine, nghost);
  m_ebislRefCoar.setMaxCoarseningRatio(m_refRat,ebisPtr);

  //define the intvectsets over which the objects live
  m_setsFine.define(m_gridsFine);
  m_setsRefCoar.define(m_gridsCoar);

  //make sets
  //global set consists of redistrad on the fine side of the coarse-fine
  //interface
  //the fine set is that within one fine box.
  //the refcoar set is that set within one refined coarse box.
  IntVectSet fineShell;
  {
    CH_TIME("make_fine_sets");
    IntVectSet fineInterior = a_eblgFine.getCoveringIVS();
    EBArith::shrinkIVS(fineInterior, m_redistRad);
    fineShell = a_eblgFine.getCoveringIVS();
    fineShell -= fineInterior;
    for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
      {
        const Box& fineBox = m_gridsFine.get(dit());
        m_setsFine[dit()]  = m_ebislFine[dit()].getIrregIVS(fineBox);
        m_setsFine[dit()] &= fineShell;
      }
  }
  {
    CH_TIME("make_coar_sets");
    for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit)
      {
        Box grownBox = grow(m_gridsRefCoar.get(dit()), m_redistRad);
        grownBox &= domainFine;
        m_setsRefCoar[dit()] = m_ebislRefCoar[dit()].getIrregIVS(grownBox);
        m_setsRefCoar[dit()] &= fineShell;
      }
  }
  defineDataHolders();
  setToZero();
}
/**********************/
void
EBFineToCoarRedist::
define(const DisjointBoxLayout& a_dblFine,
       const DisjointBoxLayout& a_dblCoar,
       const EBISLayout& a_ebislFine,
       const EBISLayout& a_ebislCoar,
       const Box& a_domainCoar,
       const int& a_nref,
       const int& a_nvar,
       int a_redistRad,
       const EBIndexSpace* const a_ebisPtr)
{
  CH_TIME("EBFineToCoarRedist::stardard_define");
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
  m_gridsRefCoar = DisjointBoxLayout();
  refine(m_gridsRefCoar, m_gridsCoar, m_refRat);

  CH_assert(a_ebisPtr->isDefined());
  int nghost = 3*m_redistRad;
  Box domainFine = refine(m_domainCoar, m_refRat);
  a_ebisPtr->fillEBISLayout(m_ebislRefCoar, m_gridsRefCoar,
                            domainFine, nghost);
  m_ebislRefCoar.setMaxCoarseningRatio(m_refRat,a_ebisPtr);

  //define the intvectsets over which the objects live
  m_setsFine.define(m_gridsFine);
  m_setsRefCoar.define(m_gridsCoar);

  //make sets
  //global set consists of redistrad on the fine side of the coarse-fine
  //interface
  //the fine set is that within one fine box.
  //the refcoar set is that set within one refined coarse box.
  {
    CH_TIME("make_fine_sets");
    for (DataIterator dit =
          m_gridsFine.dataIterator(); dit.ok(); ++dit)
      {
        const Box& fineBox = m_gridsFine.get(dit());
        IntVectSet& fineSet = m_setsFine[dit()];
        //add entire fine box
        fineSet = IntVectSet(fineBox);
        IntVectSet irregIVS = m_ebislFine[dit()].getIrregIVS(fineBox);
        fineSet &= irregIVS;
        //subtract fine box shrunk by the redistribution radius
        fineSet -= grow(fineBox, -m_redistRad);
        //subtract all other the fine boxes shifted in all
        //directions
        for (LayoutIterator lit =
              m_gridsFine.layoutIterator(); lit.ok(); ++lit)
          {
            const Box& otherBox = m_gridsFine.get(lit());
            if (otherBox != fineBox)
              {
                for (int idir = 0; idir < SpaceDim; idir++)
                  {
                    for (SideIterator sit; sit.ok(); ++sit)
                      {
                        Box shiftBox = otherBox;
                        shiftBox.shift(idir, sign(sit())*m_redistRad);
                        fineSet -= shiftBox;
                      }
                  }
              }
          }
      }
  }
  {
    CH_TIME("make_coar_sets");
    for (DataIterator dit =
          m_gridsCoar.dataIterator(); dit.ok(); ++dit)
      {
        Box grownBox = grow(m_gridsRefCoar.get(dit()), m_redistRad);
        grownBox &= domainFine;
        //find the complement of what we really want
        IntVectSet ivsComplement(grownBox);
        for (LayoutIterator litFine =
              m_gridsFine.layoutIterator(); litFine.ok(); ++litFine)
          {
            const Box& fineBox = m_gridsFine.get(litFine());
            IntVectSet fineSet(fineBox);
            //add entire fine box
            fineSet = IntVectSet(fineBox);
            Box interiorBox = grow(fineBox, -m_redistRad);
            //subtract fine box shrunk by the redistribution radius
            fineSet -= interiorBox;
            //subtract all other the fine boxes shifted in all
            //directions
            for (LayoutIterator lit =
                  m_gridsFine.layoutIterator(); lit.ok(); ++lit)
              {
                const Box& otherBox = m_gridsFine.get(lit());
                if (otherBox != fineBox)
                  {
                    for (int idir = 0; idir < SpaceDim; idir++)
                      {
                        for (SideIterator sit; sit.ok(); ++sit)
                          {
                            Box shiftBox = otherBox;
                            shiftBox.shift(idir, sign(sit())*m_redistRad);
                            fineSet -= shiftBox;
                          }
                      }
                  }
              }
            //subtract the fine set from the complement
            ivsComplement -= fineSet;
          }
        //now the set we want is the grownbox - complement
        IntVectSet& refCoarSet = m_setsRefCoar[dit()];
        refCoarSet = IntVectSet(grownBox);
        refCoarSet -= ivsComplement;
        IntVectSet irregIVS = m_ebislRefCoar[dit()].getIrregIVS(grownBox);
        refCoarSet &= irregIVS;
      }
  }
  defineDataHolders();
  setToZero();
}
/***/
void
EBFineToCoarRedist::
defineDataHolders()
{
  CH_TIME("EBFineToCoarRedist::defineDataHolders");

  EBCellFactory ebcellfactcoar(m_ebislCoar);
  m_densityCoar.define(m_gridsCoar, 1, 2*m_redistRad*IntVect::Unit, ebcellfactcoar);
  //make mass buffers
  {
    CH_TIME("make_mass_buffers");
    BaseIVFactory<Real> factFine(m_ebislFine, m_setsFine);
    IntVect noghost = IntVect::Zero;
    m_regsFine.define(m_gridsFine, m_nComp, noghost, factFine);

    BaseIVFactory<Real> factRefCoar(m_ebislRefCoar, m_setsRefCoar);
    m_regsRefCoar.define(m_gridsRefCoar, m_nComp, m_redistRad*IntVect::Unit, factRefCoar);
  }
  //define the stencils with volume weights
  //use resetWeights to do anything different
  {
    CH_TIME("make_redist_stencil");
    m_stenRefCoar.define(m_gridsRefCoar);
    m_volumeStenc.define(m_gridsRefCoar);
    m_standardStenc.define(m_gridsRefCoar);
    Box domainFine = refine(m_domainCoar, m_refRat);
    RedistStencil rdstenRefCoar(m_gridsRefCoar, m_ebislRefCoar,
                                domainFine, m_redistRad);

    for (DataIterator dit =
          m_gridsRefCoar.dataIterator(); dit.ok(); ++dit)
      {
        BaseIVFAB<VoFStencil>& stenFAB = m_stenRefCoar[dit()];
        BaseIVFAB<VoFStencil>& volStenFAB = m_volumeStenc[dit()];
        BaseIVFAB<VoFStencil>& stanStenFAB = m_standardStenc[dit()];
        const EBISBox& ebisBox = m_ebislRefCoar[dit()];
        const IntVectSet& ivs  = m_setsRefCoar[dit()];
        const BaseIVFAB<VoFStencil>& rdStenFAB = rdstenRefCoar[dit()];
        const Box& gridRefCoar = m_gridsRefCoar.get(dit());
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
                if (gridRefCoar.contains(dstVoF.gridIndex()))
                  {
                    newStencil.add(dstVoF, weight);
                  }
              }
            stenFAB(srcVoF, 0)    = newStencil;
            volStenFAB(srcVoF, 0) = newStencil;
            stanStenFAB(srcVoF,0) =  stanSten;
          }
      }
  }
}
/**********************/
void
EBFineToCoarRedist::
setToZero()
{
  CH_TIME("EBFineToCoarRedist::setToZero");
  for (DataIterator dit = m_gridsRefCoar.dataIterator();
      dit.ok(); ++dit)
    {
      m_regsRefCoar[dit()].setVal(0.0);
    }

  for (DataIterator dit = m_gridsFine.dataIterator();
      dit.ok(); ++dit)
    {
      m_regsFine[dit()].setVal(0.0);
    }
}
/**********************/
void
EBFineToCoarRedist::
increment(const BaseIVFAB<Real>& a_fineMass,
          const DataIndex& a_fineDataIndex,
          const Interval&  a_variables)
{
  CH_TIME("EBFineToCoarRedist::increment");
  BaseIVFAB<Real>& fineBuf =  m_regsFine[a_fineDataIndex];
  const EBISBox&   ebisBox = m_ebislFine[a_fineDataIndex];
  const IntVectSet& ivsLoc =  m_setsFine[a_fineDataIndex];

  CH_assert(a_fineMass.getIVS().contains(ivsLoc));
  CH_assert(fineBuf.getIVS().contains(ivsLoc));

  for (VoFIterator vofit(ivsLoc, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for (int ivar = a_variables.begin();
          ivar <= a_variables.end();  ivar++)
        {
          fineBuf(vof, ivar) += a_fineMass(vof, ivar);
        }
    }
}
/**********************/
/**********************/
void
EBFineToCoarRedist::
redistribute(LevelData<EBCellFAB>& a_coarSolution,
             const Interval& a_variables)
{
  CH_TIME("EBFineToCoarRedist::redistribute");
  Real nrefD = 1.0;
  for (int idir = 0; idir < SpaceDim; idir++)
    nrefD *= m_refRat;
  //copy the buffer to the coarse layout
  m_regsFine.copyTo(a_variables, m_regsRefCoar, a_variables);
  //redistribute the refined coarse registers to the coarse solution
  int ibox = 0;
  for (DataIterator dit = m_gridsCoar.dataIterator(); dit.ok(); ++dit)
    {
      const BaseIVFAB<Real>& regRefCoar = m_regsRefCoar[dit()];
      const IntVectSet& ivsRefCoar = m_setsRefCoar[dit()];
      const EBISBox& ebisBoxRefCoar = m_ebislRefCoar[dit()];
      const EBISBox& ebisBoxCoar = m_ebislCoar[dit()];
      const BaseIVFAB<VoFStencil>& stenFAB = m_stenRefCoar[dit()];

      EBCellFAB& solFAB = a_coarSolution[dit()];

      for (VoFIterator vofit(ivsRefCoar, ebisBoxRefCoar.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& srcVoFFine = vofit();
          const VoFStencil& vofsten = stenFAB(srcVoFFine, 0);
          for (int isten = 0; isten < vofsten.size(); isten++)
            {
              const Real& weight = vofsten.weight(isten);
              const VolIndex& dstVoFFine = vofsten.vof(isten);
              VolIndex dstVoFCoar =
                m_ebislRefCoar.coarsen(dstVoFFine,m_refRat, dit());
              //by construction...
              CH_assert(m_gridsCoar.get(dit()).contains(dstVoFCoar.gridIndex()));
              Real dstVolFracFine = ebisBoxRefCoar.volFrac(dstVoFFine);
              Real dstVolFracCoar = ebisBoxCoar.volFrac(dstVoFCoar);
              Real denom = dstVolFracCoar*nrefD;
              for (int ivar = a_variables.begin();
                  ivar <= a_variables.end();  ivar++)
                {
                  Real dmFine = regRefCoar(srcVoFFine, ivar);
                  //ucoar+= massfine/volcoar, ie.
                  //ucoar+= (wcoar*dmCoar*volFracfine/volfraccoar)=massfine/volcoar
                  Real dUCoar = dmFine*weight*dstVolFracFine/denom;
                  solFAB(dstVoFCoar, ivar) += dUCoar;
                }
            }
        }
      ibox++;
    }
}
/**********************/
bool
EBFineToCoarRedist::
isDefined() const
{
  return m_isDefined;
}
/**********************/
#include "NamespaceFooter.H"
