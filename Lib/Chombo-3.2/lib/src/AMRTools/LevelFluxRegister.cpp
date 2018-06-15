#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LevelFluxRegister.H"
#include "LevelFluxRegisterF_F.H"
#include "CH_Timer.H"
#include "LayoutIterator.H"
//#include "DebugOut.H"
#include "Copier.H"
#include "parstream.H"
#include "BoxIterator.H"
#include "DebugOut.H"

#include "NamespaceHeader.H"

IntVect ivdebnoeb(D_DECL6(18, 19, 0, 0, 0, 0));
int debugdir = 1;

bool LevelFluxRegister::s_verbose = false;
void
LevelFluxRegister::poutCoarseRegisters() const
{
  for (DataIterator dit = m_coarFlux.dataIterator(); dit.ok(); ++dit)
    {
      pout() << "dumping coarse registers where non-zero" << endl;
      const FArrayBox& fab = m_coarFlux[dit()];
      const Box&       box = fab.box();
      for (BoxIterator bit(box); bit.ok(); ++bit)
        {
          if (Abs(fab(bit(), 0)) > 0.0)
            {
              pout() << bit()  << "  ";
              for (int ivar = 0; ivar < fab.nComp(); ivar++)
                {
                  pout() << fab(bit(), ivar) << "  ";
                }
              pout() << endl;
            }
        }
    }
}

void
LevelFluxRegister::poutFineRegisters() const
{
  for (DataIterator dit = m_fineFlux.dataIterator(); dit.ok(); ++dit)
    {
      pout() << "dumping fine registers where non-zero" << endl;
      const FArrayBox& fab = m_fineFlux[dit()];
      const Box&       box = fab.box();
      for (BoxIterator bit(box); bit.ok(); ++bit)
        {
          if (Abs(fab(bit(), 0)) > 0.0)
            {
              pout() << bit()  << "  ";
              for (int ivar = 0; ivar < fab.nComp(); ivar++)
                {
                  pout() << fab(bit(), ivar) << "  ";
                }
              pout() << endl;
            }
        }
    }
}
LevelFluxRegister::LevelFluxRegister(const DisjointBoxLayout& a_dblFine,
                                     const DisjointBoxLayout& a_dblCoar,
                                     const ProblemDomain&     a_dProblem,
                                     int                      a_nRefine,
                                     int                      a_nComp,
                                     bool                     a_scaleFineFluxes,
                                     int                      a_unrefinedDirection) :
  m_isDefined(FluxRegUndefined)
{
  define(a_dblFine, a_dblCoar, a_dProblem,
         a_nRefine, a_nComp, a_scaleFineFluxes,
         a_unrefinedDirection);
}

LevelFluxRegister::LevelFluxRegister(const DisjointBoxLayout& a_dblFine,
                                     const DisjointBoxLayout& a_dblCoar,
                                     const Box&               a_dProblem,
                                     int                      a_nRefine,
                                     int                      a_nComp,
                                     bool                     a_scaleFineFluxes,
                                     int                      a_unrefinedDirection) :
  m_isDefined(FluxRegUndefined)
{
  ProblemDomain physDomain(a_dProblem);
  define(a_dblFine, a_dblCoar, physDomain,
         a_nRefine, a_nComp, a_scaleFineFluxes,
         a_unrefinedDirection);
}

void
LevelFluxRegister::define(const DisjointBoxLayout& a_dbl,
                          const DisjointBoxLayout& a_dblCoarse,
                          const ProblemDomain&     a_dProblem,
                          int                      a_nRefine,
                          int                      a_nComp,
                          bool                     a_scaleFineFluxes,
                          int                      a_unrefinedDirection)
{
  CH_TIME("LevelFluxRegister::define");
  m_isDefined = FluxRegDefined;  // Basically, define was called
  m_unrefinedDirection = a_unrefinedDirection;
  m_nRefine = a_nRefine*IntVect::Unit;
  if ((m_unrefinedDirection >= 0) && (m_unrefinedDirection < SpaceDim))
    {
      m_nRefine[m_unrefinedDirection] =1;
    }

  m_scaleFineFluxes = a_scaleFineFluxes;
  DisjointBoxLayout coarsenedFine;
  coarsen(coarsenedFine, a_dbl, m_nRefine);
  m_noRealCoarseFineInterface = false;
#if 0
  // This doesn't work for multi-block calculations, which are
  // not properly nested. -JNJ

  //begin temporary optimization.  bvs
  int numPts=0;
  for (LayoutIterator lit = a_dblCoarse.layoutIterator(); lit.ok(); ++lit)
    {
      numPts += a_dblCoarse[lit].numPts();
    }
  for (LayoutIterator lit = coarsenedFine.layoutIterator(); lit.ok(); ++lit)
    {
      numPts -= coarsenedFine[lit].numPts();
    }

   if (numPts == 0)
     {
       m_coarFlux.clear();
       m_noRealCoarseFineInterface = true;
       // OK, fine region completely covers coarse region.  no registers.
       return;
     }

#endif

  //end temporary optimization.   bvs
  m_coarFlux.define( a_dblCoarse, a_nComp);
  m_isDefined |= FluxRegCoarseDefined;
  m_domain = a_dProblem;
  ProblemDomain coarsenedDomain(a_dProblem);
  coarsenedDomain.coarsen(m_nRefine);

  m_fineFlux.define( coarsenedFine, a_nComp, IntVect::Unit);
  m_isDefined |= FluxRegFineDefined;

  m_reverseCopier.ghostDefine(coarsenedFine, a_dblCoarse,
                              coarsenedDomain, IntVect::Unit);

  for (int i=0; i<CH_SPACEDIM; i++)
    {
      m_coarseLocations[i].define(a_dblCoarse);
      m_coarseLocations[i+CH_SPACEDIM].define(a_dblCoarse);
    }


  DataIterator dC   = a_dblCoarse.dataIterator();
  LayoutIterator dF = coarsenedFine.layoutIterator();

  for (dC.begin(); dC.ok(); ++dC)
    {
      const Box& cBox = a_dblCoarse.get(dC);

      for (dF.begin(); dF.ok(); ++dF)
        {
          const Box& fBox  = coarsenedFine.get(dF);

          if (fBox.bigEnd(0)+1 < cBox.smallEnd(0))
            {
              //can skip this box since they cannot intersect, due to sorting
            }
          else if (fBox.smallEnd(0)-1 > cBox.bigEnd(0))
            {
              //skip to end, since all the rest of boxes will not intersect either
              dF.end();
            }
          else
            {

              for (int i=0; i<CH_SPACEDIM; i++)
                {
                  Vector<Box>& lo = m_coarseLocations[i][dC];
                  Vector<Box>& hi = m_coarseLocations[i+CH_SPACEDIM][dC];

                  Box loBox = adjCellLo(fBox, i, 1);
                  Box hiBox = adjCellHi(fBox, i, 1);
                  if (cBox.intersectsNotEmpty(loBox)) lo.push_back(loBox & cBox);
                  if (cBox.intersectsNotEmpty(hiBox)) hi.push_back(hiBox & cBox);
                }
            }
        }
    }

  Box domainBox = coarsenedDomain.domainBox();

  if (a_dProblem.isPeriodic())
    {
      Vector<Box> periodicBoxes[2*CH_SPACEDIM];
      for (dF.begin(); dF.ok(); ++dF)
        {
          const Box& fBox  = coarsenedFine.get(dF);
          for (int i=0; i<CH_SPACEDIM; i++)
            {
              if (a_dProblem.isPeriodic(i))
                {
                  if (fBox.smallEnd(i) == domainBox.smallEnd(i))
                    periodicBoxes[i].push_back(adjCellLo(fBox, i, 1));
                  if (fBox.bigEnd(i) == domainBox.bigEnd(i))
                    periodicBoxes[i+CH_SPACEDIM].push_back(adjCellHi(fBox, i, 1));
                }
            }
        }
      for (int i=0; i<CH_SPACEDIM; i++)
        {
          Vector<Box>& loV = periodicBoxes[i];
          Vector<Box>& hiV = periodicBoxes[i+CH_SPACEDIM];
          int size = domainBox.size(i);
          for (int j=0; j<loV.size(); j++) loV[j].shift(i, size);
          for (int j=0; j<hiV.size(); j++) hiV[j].shift(i, -size);
        }
      for (dC.begin(); dC.ok(); ++dC)
        {
          const Box& cBox = a_dblCoarse.get(dC);
          for (int i=0; i<CH_SPACEDIM; i++)
            if (a_dProblem.isPeriodic(i))
              {
                Vector<Box>& loV = periodicBoxes[i];
                Vector<Box>& hiV = periodicBoxes[i+CH_SPACEDIM];

                if (cBox.smallEnd(i) == domainBox.smallEnd(i) )
                  {
                    Vector<Box>& hi = m_coarseLocations[i+CH_SPACEDIM][dC];
                    for (int j=0; j<hiV.size(); j++)
                      {
                        if (cBox.intersectsNotEmpty(hiV[j])) hi.push_back(cBox & hiV[j]);
                      }
                  }
                if (cBox.bigEnd(i) == domainBox.bigEnd(i) )
                  {
                    Vector<Box>& lo = m_coarseLocations[i][dC];
                    for (int j=0; j<loV.size(); j++)
                      {
                        if (cBox.intersectsNotEmpty(loV[j])) lo.push_back(cBox & loV[j]);
                      }
                  }

              }
        }
    }

}

void
LevelFluxRegister::define(const DisjointBoxLayout& a_dbl,
                          const DisjointBoxLayout& a_dblCoarse,
                          const ProblemDomain&     a_dProblem,
                          int                      a_nRefine,
                          int                      a_nComp)
{
  // Just call the more general version with scaling set to true.
  define(a_dbl, a_dblCoarse, a_dProblem, a_nRefine, a_nComp, true);
}


bool LevelFluxRegister::hasCF(const DataIndex& a_fineDataIndex, Side::LoHiSide) const
{
  // return m_coarFlux.isDefined();
  // Has CF if everything defined
  return isAllDefined();
}

bool LevelFluxRegister:: hasCF(const DataIndex& a_coarseIndex) const
{
  // return m_coarFlux.isDefined();
  // Has CF if everything defined
  return isAllDefined();
}


void
LevelFluxRegister::define(const DisjointBoxLayout& a_dbl,
                          const DisjointBoxLayout& a_dblCoarse,
                          const Box&               a_dProblem,
                          int                      a_nRefine,
                          int                      a_nComp,
                          bool                     a_scaleFineFluxes)
{
  ProblemDomain physDomain(a_dProblem);
  define(a_dbl, a_dblCoarse, physDomain,
         a_nRefine, a_nComp, a_scaleFineFluxes);
}

LevelFluxRegister::LevelFluxRegister()
:
  m_isDefined(FluxRegUndefined),
  m_scaleFineFluxes(true)
{
}

LevelFluxRegister::~LevelFluxRegister()
{
}

bool
LevelFluxRegister::isDefined() const
{
  return (m_isDefined & FluxRegDefined);
}

bool
LevelFluxRegister::isAllDefined() const
{
  bool isFluxRegDefined = m_isDefined & FluxRegDefined;
  bool isFluxRegFineDefined = m_isDefined & FluxRegFineDefined;
  bool isFluxRegCoarseDefined = m_isDefined & FluxRegCoarseDefined;
  return (isFluxRegDefined && isFluxRegFineDefined && isFluxRegCoarseDefined);
}

void
LevelFluxRegister::undefine()
{
  m_isDefined = FluxRegUndefined;
}

void
LevelFluxRegister::setToZero()
{
  CH_TIME("LevelFluxRegister::setToZero");
  if (m_isDefined & FluxRegCoarseDefined)
    {
      for (DataIterator d = m_coarFlux.dataIterator(); d.ok(); ++d)
        {
          m_coarFlux[d].setVal(0.0);
        }
    }
  if (m_isDefined & FluxRegFineDefined)
    {
      for (DataIterator d = m_fineFlux.dataIterator(); d.ok(); ++d)
        {
          m_fineFlux[d].setVal(0.0);
        }
    }
}

void
LevelFluxRegister::incrementCoarse(const FArrayBox& a_coarseFlux,
                                   Real a_scale,
                                   const DataIndex& a_coarseDataIndex,
                                   const Interval& a_srcInterval,
                                   const Interval& a_dstInterval,
                                   int a_dir)

{
  CH_assert(isDefined());
  incrementCoarse(a_coarseFlux, a_scale, a_coarseDataIndex, a_srcInterval,
                  a_dstInterval, a_dir, Side::Lo);
  incrementCoarse(a_coarseFlux, a_scale, a_coarseDataIndex, a_srcInterval,
                  a_dstInterval, a_dir, Side::Hi);
}

void
LevelFluxRegister::incrementCoarse(const FArrayBox& a_coarseFlux,
                                   Real a_scale,
                                   const DataIndex& a_coarseDataIndex,
                                   const Interval& a_srcInterval,
                                   const Interval& a_dstInterval,
                                   int a_dir,
                                   Side::LoHiSide a_sd)
{
  CH_assert(isDefined());
  if (!(m_isDefined & FluxRegCoarseDefined)) return;
  CH_TIME("LevelFluxRegister::incrementCoarse");

  const Vector<Box>& intersect =
     m_coarseLocations[a_dir+a_sd*CH_SPACEDIM][a_coarseDataIndex];

  FArrayBox& coarse = m_coarFlux[a_coarseDataIndex];

  // We cast away the constness in a_coarseFlux for the scope of this function. This
  // should be acceptable, since at the end of the day there is no change to it. -JNJ
  FArrayBox& coarseFlux = const_cast<FArrayBox&>(a_coarseFlux); // Muhahaha.
  coarseFlux.shiftHalf(a_dir, sign(a_sd));
  Real scale = -sign(a_sd)*a_scale;
  int s = a_srcInterval.begin();
  int d = a_dstInterval.begin();
  int size = a_srcInterval.size();

  for (int b=0; b<intersect.size(); ++b)
    {
      const Box& box = intersect[b];
      Vector<Real> regbefore(coarse.nComp());
      Vector<Real> regafter(coarse.nComp());
      if (s_verbose && (a_dir == debugdir) && box.contains(ivdebnoeb))
        {
          for (int ivar = 0; ivar < coarse.nComp(); ivar++)
            {
              regbefore[ivar] = coarse(ivdebnoeb, ivar);
            }
        }

      coarse.plus(coarseFlux, box, box, scale, s, d, size);

      if (s_verbose && (a_dir == debugdir) && box.contains(ivdebnoeb))
        {
          for (int ivar = 0; ivar < coarse.nComp(); ivar++)
            {
              regafter[ivar] = coarse(ivdebnoeb, ivar);
            }

          pout() << "levelfluxreg::incrementCoar: scale = " << scale << ", ";
          for (int ivar = 0; ivar < coarse.nComp(); ivar++)
            {
              pout() << " input flux = " << coarseFlux(ivdebnoeb, ivar) << ", ";
              pout() << " reg before = " <<             regbefore[ivar] << ", ";
              pout() << " reg after  = " <<                regafter[ivar] << ", ";
            }
          pout() << endl;

        }
    }

  coarseFlux.shiftHalf(a_dir, - sign(a_sd));
}

void
LevelFluxRegister::incrementFine(const FArrayBox& a_fineFlux,
                                 Real a_scale,
                                 const DataIndex& a_fineDataIndex,
                                 const Interval& a_srcInterval,
                                 const Interval& a_dstInterval,
                                 int a_dir)
{
  CH_assert(isDefined());
  incrementFine(a_fineFlux, a_scale, a_fineDataIndex, a_srcInterval,
                a_dstInterval, a_dir, Side::Lo);
  incrementFine(a_fineFlux, a_scale, a_fineDataIndex, a_srcInterval,
                a_dstInterval, a_dir, Side::Hi);
}

void
LevelFluxRegister::incrementFine(const FArrayBox& a_fineFlux,
                                 Real a_scale,
                                 const DataIndex& a_fineDataIndex,
                                 const Interval& a_srcInterval,
                                 const Interval& a_dstInterval,
                                 int a_dir,
                                 Side::LoHiSide a_sd)
{
  CH_assert(isDefined());
  if (!(m_isDefined & FluxRegFineDefined)) return;
  CH_assert(a_srcInterval.size() == a_dstInterval.size());

  CH_TIME("LevelFluxRegister::incrementFine");

  // We cast away the constness in a_coarseFlux for the scope of this function. This
  // should be acceptable, since at the end of the day there is no change to it. -JNJ
  FArrayBox& fineFlux = const_cast<FArrayBox&>(a_fineFlux); // Muhahaha.
  fineFlux.shiftHalf(a_dir, sign(a_sd));
  Real denom = 1.0;

  if (m_scaleFineFluxes)
  {
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (idir != a_dir)
          {
            denom *= m_nRefine[idir];
          }
      }
  }

  Real scale = sign(a_sd)*a_scale/denom;

  FArrayBox& cFine = m_fineFlux[a_fineDataIndex];
  //  FArrayBox  cFineFortran(cFine.box(), cFine.nComp());
  //  cFineFortran.copy(cFine);

  Box clipBox = m_fineFlux.box(a_fineDataIndex);
  clipBox.refine(m_nRefine);
  Box fineBox;
  if (a_sd == Side::Lo)
    {
      fineBox = adjCellLo(clipBox, a_dir, 1);
      fineBox &= fineFlux.box();
    }
  else
    {
      fineBox = adjCellHi(clipBox,a_dir,1);
      fineBox &= fineFlux.box();
    }
#if 0
  for (BoxIterator b(fineBox); b.ok(); ++b)
     {
       int s = a_srcInterval.begin();
       int d = a_dstInterval.begin();
       for (;s<=a_srcInterval.end(); ++s, ++d)
         {
           cFine(coarsen(b(), m_nRefine), d) += scale*fineFlux(b(), s);
         }
     }
#else
   // shifting to ensure fineBox is in the positive quadrant, so IntVect coarening
   // is just integer division.

  const Box& box = coarsen(fineBox, m_nRefine);
  Vector<Real> regbefore(cFine.nComp());
  Vector<Real>  regafter(cFine.nComp());
  if (s_verbose && (a_dir == debugdir) && box.contains(ivdebnoeb))
    {
      for (int ivar = 0; ivar < cFine.nComp(); ivar++)
        {
          regbefore[ivar] = cFine(ivdebnoeb, ivar);
        }
    }


   const IntVect& iv=fineBox.smallEnd();
   IntVect civ = coarsen(iv, m_nRefine);
   int srcComp = a_srcInterval.begin();
   int destComp = a_dstInterval.begin();
   int ncomp = a_srcInterval.size();
   FORT_INCREMENTFINE(CHF_CONST_FRA_SHIFT(fineFlux, iv),
                      CHF_FRA_SHIFT(cFine, civ),
                      CHF_BOX_SHIFT(fineBox, iv),
                      CHF_CONST_INTVECT(m_nRefine),
                      CHF_CONST_REAL(scale),
                      CHF_CONST_INT(srcComp),
                      CHF_CONST_INT(destComp),
                      CHF_CONST_INT(ncomp));


  if (s_verbose && (a_dir == debugdir) && box.contains(ivdebnoeb))
    {
      for (int ivar = 0; ivar < cFine.nComp(); ivar++)
        {
          regafter[ivar] = cFine(ivdebnoeb, ivar);
        }
    }

  if (s_verbose && (a_dir == debugdir) && box.contains(ivdebnoeb))
    {

      pout() << "levelfluxreg::incrementFine: scale = " << scale << endl;
      Box refbox(ivdebnoeb, ivdebnoeb);
      refbox.refine(m_nRefine);
      refbox &= fineBox;
      if (!refbox.isEmpty())
        {
          pout() << "fine fluxes = " << endl;
          for (BoxIterator  bit(refbox); bit.ok(); ++bit)
            {
              for (int ivar = 0; ivar < cFine.nComp(); ivar++)
                {
                  pout() << "iv = " << bit() << "(";
                  for (int ivar = 0; ivar < cFine.nComp(); ivar++)
                    {
                      pout() << fineFlux(bit(), ivar);
                    }
                  pout() << ")" << endl;
                }
            }
        }

      for (int ivar = 0; ivar < cFine.nComp(); ivar++)
        {
          pout() << " reg before = " <<               regbefore[ivar] << ", ";
          pout() << " reg after  = " <<                regafter[ivar] << ", ";
        }
      pout() << endl;

    }

#endif
  fineFlux.shiftHalf(a_dir, - sign(a_sd));
}

void
LevelFluxRegister::reflux(LevelData<FArrayBox>& a_uCoarse,
                          Real a_scale)
{
  Interval interval(0, a_uCoarse.nComp()-1);
  reflux(a_uCoarse, interval, interval, a_scale);
}

void LevelFluxRegister::reflux(LevelData<FArrayBox>& a_uCoarse,
                               const Interval& a_coarseVectorIntv,
                               Real a_scale)
{
  MayDay::Error("This function was never implemented as advertised and things still seem to work correctly");
  Interval interval(0, a_uCoarse.nComp()-1);
  reflux(a_uCoarse, interval, interval, a_scale);
}

class AddOp : public LDOperator<FArrayBox>
{
public:
  Real scale;
  AddOp():scale(1.0)
  {
  }

  virtual void linearIn(FArrayBox& arg,  void* buf, const Box& R,
                        const Interval& comps) const
  {
    Real* buffer = (Real*)buf;
    if (scale != 1.0)
    {
      ForAllXBNNnoindx(Real, arg, R, comps.begin(), comps.size())
      {
        argR+=*buffer * scale;
        buffer++;
      } EndFor
    }
    else
    {
      ForAllXBNNnoindx(Real, arg, R, comps.begin(), comps.size())
      {
        argR+=*buffer;
        buffer++;
      } EndFor
    }
  }

  void op(FArrayBox& dest,
          const Box& RegionFrom,
          const Interval& Cdest,
          const Box& RegionTo,
          const FArrayBox& src,
          const Interval& Csrc) const
  {
    if (scale != 1.0)
      dest.plus(src, RegionFrom, RegionTo, scale, Csrc.begin(), Cdest.begin(), Cdest.size());
    else
      dest.plus(src, RegionFrom, RegionTo, Csrc.begin(), Cdest.begin(), Cdest.size());
  }
};

void LevelFluxRegister::reflux(
                               LevelData<FArrayBox>& a_uCoarse,
                               const Interval&       a_coarse_interval,
                               const Interval&       a_flux_interval,
                               Real                  a_scale )
{
  CH_assert(isDefined() );
  if (m_noRealCoarseFineInterface)
    return;
  CH_TIME("LevelFluxRegister::reflux");
  for (DataIterator dit(a_uCoarse.dataIterator()); dit.ok(); ++dit)
    {
      FArrayBox& u = a_uCoarse[dit];
      FArrayBox& coarse = m_coarFlux[dit];
      u.plus(coarse, -a_scale, a_flux_interval.begin(),
             a_coarse_interval.begin(), a_coarse_interval.size() );
    }
  AddOp op;
  op.scale = -a_scale;
  m_fineFlux.copyTo(a_flux_interval, a_uCoarse, a_coarse_interval, m_reverseCopier, op);

}

void LevelFluxRegister::reflux(
                               LevelData<FArrayBox>& a_uCoarse,
                               Real                  a_scale,
                               const Interval&       a_coarse_interval,
                               const Interval&       a_flux_interval,
                               const LevelData<FArrayBox>& a_beta
                               )
{
  CH_assert(isAllDefined() );
  CH_TIME("LevelFluxRegister::reflux");
  CH_assert(a_beta.nComp() ==1);
  LevelData<FArrayBox> increment(a_uCoarse.disjointBoxLayout(), a_uCoarse.nComp(), a_uCoarse.ghostVect());
  for (DataIterator dit(a_uCoarse.dataIterator()); dit.ok(); ++dit)
    {
      increment[dit()].setVal(0.0);
    }

  reflux(increment, a_coarse_interval, a_flux_interval, a_scale );
  for (DataIterator dit(a_uCoarse.dataIterator()); dit.ok(); ++dit)
    {
      for (int icomp = 0; icomp < increment.nComp(); icomp++)
        {
          int srccomp = 0;
          int dstcomp = icomp;
          int ncomp = 1;
          increment[dit()].mult(a_beta[dit()], srccomp, dstcomp, ncomp);
        }
      int srccomp = 0;
      int dstcomp = 0;
      int ncomp = a_uCoarse.nComp();
      a_uCoarse[dit()].plus(increment[dit()], srccomp, dstcomp, ncomp);
    }
}
#include "NamespaceFooter.H"
