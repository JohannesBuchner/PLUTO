#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "parstream.H"
#include "ParmParse.H"
#include "PolyGeom.H"

#include "DebugOut.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "EBCellFAB.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "EBAMRIO.H"

#include "AMRLevel.H"
#include "EBAMRGodunov.H"
#include "EBCellFactory.H"
#include "BaseIVFactory.H"
#include "VoFIterator.H"
#include "EBPatchGodunovF_F.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "EBFastFR.H"
#include "NamespaceHeader.H"

int  EBAMRGodunov::s_NewPlotFile = 0;
bool EBAMRGodunov::s_isLoadBalanceSet = false;
LoadBalanceFunc EBAMRGodunov::s_loadBalance  = NULL;
IntVect ivdebamrg(D_DECL(16, 5, 0));
int debuglevel = 1;

/***************************/
void
EBAMRGodunov::
tagAll(bool a_tagAll)
{
  m_tagAll = a_tagAll;
}
/***************************/
void
EBAMRGodunov::
doSmushing(bool a_doSmushing)
{
  m_doSmushing = a_doSmushing;
}
/***************************/
void
EBAMRGodunov::
doRZCoords(bool a_doRZCoords)
{
  m_doRZCoords = a_doRZCoords;
}
/***************************/
void
EBAMRGodunov::
hasSourceTerm(bool a_hasSourceTerm)
{
  m_hasSourceTerm = a_hasSourceTerm;
}
/***************************/
void
EBAMRGodunov::
useMassRedistribution(bool a_useMassRedist)
{
  m_useMassRedist = a_useMassRedist;
}
/***************************/
EBAMRGodunov::EBAMRGodunov()
{
  m_cfl = 0.8;
  m_tagAll = false;
  m_useMassRedist = true;
  m_doRZCoords = false;
  m_hasSourceTerm = false;
  m_doSmushing = true;
  m_origin = RealVect::Zero;
  m_dx = RealVect::Unit;
  m_aspect = RealVect::Unit;
  m_domainLength = RealVect::Unit;
  m_refineThresh = 0.2;
  m_initial_dt_multiplier = 0.1;
  m_ebPatchGodunov = NULL;
  m_redistRad = 1;
  m_isDefined = false;
}
/***************************/
void EBAMRGodunov::redistRadius(int a_redistRad)
{
  m_redistRad = a_redistRad;
}
/***************************/
EBAMRGodunov::~EBAMRGodunov()
{
  if (m_ebPatchGodunov != NULL)
    delete m_ebPatchGodunov;
}
/***************************/
void EBAMRGodunov::define(AMRLevel*  a_coarser_level_ptr,
                          const Box& a_problem_domain,
                          int        a_level,
                          int        a_ref_ratio)
{
  MayDay::Error("EBAMRGodunov::define -\n\tShould never be called with a Box for a problem domain");
}
/***************************/
void EBAMRGodunov::define(AMRLevel*  a_coarser_level_ptr,
                          const ProblemDomain& a_problem_domain,
                          int        a_level,
                          int        a_ref_ratio)
{
  CH_TIME("EBAMRGodunov::define");
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::define, level=" << a_level << endl;
    }

  m_isDefined = true;
  AMRLevel::define(a_coarser_level_ptr,
                   a_problem_domain,
                   a_level,
                   a_ref_ratio);
  m_domainBox = m_problem_domain.domainBox();

  if (a_coarser_level_ptr != NULL)
    {
      EBAMRGodunov* amrg_ptr =
        dynamic_cast<EBAMRGodunov*>(a_coarser_level_ptr);
      if (amrg_ptr == NULL)
        {
          pout() << "EBAMRG::define:cannot cast  to derived class"
                 << endl;
          MayDay::Error();
        }

      m_cfl = amrg_ptr->m_cfl;
      m_domainLength = amrg_ptr->m_domainLength;
      m_refineThresh = amrg_ptr->m_refineThresh;
      m_tagBufferSize = amrg_ptr->m_tagBufferSize;
    }
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      m_dx[idir] = m_domainLength[idir]/m_domainBox.size(idir);
    }
  m_nGhost = 4;

  if (m_ebPatchGodunov != NULL)
    delete m_ebPatchGodunov;

  m_ebPatchGodunov = m_ebPatchGodunovFactory->create();
  m_ebPatchGodunov->define(m_problem_domain, m_dx);

  m_nComp = m_ebPatchGodunov->numConserved();
  m_stateNames = m_ebPatchGodunov->stateNames();
  m_primNames = m_ebPatchGodunov->primNames();
}

void
EBAMRGodunov::dumpDebug()
{
  dumpDebug("arg");
}

void
EBAMRGodunov::dumpDebug(const string& a_debstring)
{
  int ilev = m_level;
  //  if (ilev == debuglevel)
  if (0)
    {
      pout()    << setprecision(10)
                << setiosflags(ios::showpoint)
                << setiosflags(ios::scientific);
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          const Box& grid =m_grids.get(dit());
          if (grid.contains(ivdebamrg))
            {
              pout() << a_debstring << ", lev = " << ilev << ", iv=" << ivdebamrg<<  ":";
              for (int ivar = 0; ivar < m_stateNew.nComp(); ivar++)
                {
                  VolIndex vof(ivdebamrg, 0);
                  Real sold = m_stateOld[dit()](vof, ivar);
                  Real snew = m_stateNew[dit()](vof, ivar);
                  pout() << " ( " << snew << ", "  << sold << " ) ";

                }
              pout() << endl;
            }
        }
    }
}
/***************************/
Real EBAMRGodunov::advance()
{
  CH_TIME("EBAMRGodunov::advance");
  EBPatchGodunov::s_whichLev = m_level;
  dumpDebug(string("going into advance"));
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRGodunov advance for level =" << m_level << ", with dt = " << m_dt << endl;
    }
  {
    CH_TIME("copy new to old");
    m_stateNew.copyTo(m_stateNew.interval(),
                      m_stateOld,
                      m_stateOld.interval());
  }

  Real new_dt = 0.0;
  //set up arguments to step
  //undefined lfr in case we need it
  EBFluxRegister lfr;
  //undefined leveldata in case we need it
  const LevelData<EBCellFAB> ld;

  //set arguments to dummy arguments and
  //then fix if they are available
  EBFluxRegister* coarFR = &lfr;
  EBFluxRegister* fineFR = &lfr;
  const LevelData<EBCellFAB>* coarDataOld = &ld;
  const LevelData<EBCellFAB>* coarDataNew = &ld;

  Real told = 0.0;
  Real tnew = 0.0;
  if (m_hasCoarser)
    {
      EBAMRGodunov* coarPtr = getCoarserLevel();
      //recall that my flux register goes between my
      //level and the next finer level
      coarFR = &coarPtr->m_ebFluxRegister;
      coarDataOld = &coarPtr->m_stateOld;
      coarDataNew = &coarPtr->m_stateNew;
      tnew = coarPtr->m_time;
      told = tnew - coarPtr->m_dt;
      //time should never be greater than the newest coarse
      //time.  time might be very slightly smaller than
      //told because of the above subtraction.
      Real eps = 1.0e-10;
      if ((m_time > tnew) || (m_time < (told - eps)))
        {
          MayDay::Error("out of bounds time input to AMRGodunov");
        }
      //correct for said floating-point nastiness
      m_time = Max(m_time, told);
    }
  if (m_hasFiner)
    {
      //recall that my flux register goes between my
      //level and the next finer level
      fineFR = &m_ebFluxRegister;
    }

#ifndef NDEBUG
  if (!m_hasCoarser && (s_verbosity > 1))
    {
      Real summass;
      int densityIndex = m_ebPatchGodunov->densityIndex();
      sumConserved(summass,  densityIndex);

      pout() << "sum mass = " << summass << endl;
    }
#endif

  {
    CH_TIME("levelgodunov step");
    EBPatchGodunov::setCurLevel(m_level);
    new_dt = m_ebLevelGodunov.step(m_stateNew,
                                   m_massDiff,
                                   *fineFR,
                                   *coarFR,
                                   *coarDataOld,
                                   *coarDataNew,
                                   m_time,
                                   told,
                                   tnew,
                                   m_dt);
  }
  dumpDebug(string("after levelgodunov"));

  if ((s_verbosity > 2))
    {
      pout() << "for this proc max  wave speed = " << EBPatchGodunov::getMaxWaveSpeed()
             << " at cell =  " << EBPatchGodunov::getMaxWaveSpeedIV() << endl;
    }
  // flux register manipulation happens in levelgodunov
  // level redistribution happens in levelgodunov

  Interval interv(0, m_nComp-1);
  // increment redistribution register between this level
  // and next coarser level
  {
    CH_TIME("coarse-fine guano");
    if (m_hasCoarser)
      {
        for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
          {
            m_ebFineToCoarRedist.increment(m_massDiff[dit()], dit(), interv);
          }
      }

    //initialize redistirbution register between this level
    // and next finer level
    //this includes re-redistribution registers
    if (m_hasFiner)
      {

        m_ebCoarToFineRedist.setToZero();
        m_ebCoarToCoarRedist.setToZero();
        for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
          {
            BaseIVFAB<Real>& massDiffFAB = m_massDiff[dit()];
            m_ebCoarToFineRedist.increment(massDiffFAB, dit(), interv);
            m_ebCoarToCoarRedist.increment(massDiffFAB, dit(), interv);
          }
      }
  }

  //deal with time and time step
  m_time += m_dt;
  Real return_dt = m_cfl * new_dt;

  //save stable timestep to save computational effort
  m_dtNew = return_dt;

  dumpDebug(string("leaving advance   "));
  return return_dt;
}
/***************************/
void EBAMRGodunov::postTimeStep()
{
  CH_TIME("EBAMRGodunov::postTimeStep");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRGodunov postTimeStep for level " << m_level << endl;
    }
  Interval interv(0, m_nComp-1);
  if (m_hasCoarser)
    {
      if (m_doSmushing)
        {
          //redistibute to coarser level
          EBAMRGodunov* coarPtr = getCoarserLevel();
          //if use mass weighting, need to
          //fix weights of redistribution object
          //EBFineToCoarRedist::s_verbose = false;
          //if (!m_hasFiner && m_hasCoarser)
          //  {
          //    EBFineToCoarRedist::s_verbose = true;
          //    EBFineToCoarRedist::s_ivdebug = ivdebamrg;
          //  }

          if (m_useMassRedist)
            {
              int densevar = m_ebPatchGodunov->densityIndex();
              coarPtr->m_stateNew.exchange(interv);

              //if (EBFineToCoarRedist::s_verbose)
              //  {
              //    EBFineToCoarRedist::s_ivdebug = IntVect(D_DECL(36, 11, 0));
              //  }

              m_ebFineToCoarRedist.resetWeights(coarPtr->m_stateNew, densevar);
            }

          m_ebFineToCoarRedist.redistribute(coarPtr->m_stateNew, interv);

          //EBFineToCoarRedist::s_verbose = false;

          coarPtr->dumpDebug(string("after finetocoar  "));
        }

      m_ebFineToCoarRedist.setToZero();
    }
  if (m_hasFiner)
    {
      EBAMRGodunov* finePtr = getFinerLevel();
      // reflux from finer level solution
      CH_assert(Abs(m_dx[0]-m_dx[1])<1.e-9);
      Real scale = -1.0/m_dx[0];
//      int ilev = m_level;
//      if (ilev == debuglevel)
//        {
//          EBFastFR::s_verbose = true;
//        }

      m_ebFluxRegister.reflux(m_stateNew, interv, scale);

      //      EBFastFR::s_verbose = false;

      dumpDebug(string("after reflux      "));

      //do averaging from finer level
      finePtr->m_ebCoarseAverage.average(m_stateNew,
                                         finePtr->m_stateNew,
                                         interv);


      dumpDebug(string("after reflux/ave  "));


      //the flux register must modify the redistribution
      //registers


      m_ebFluxRegister.incrementRedistRegister(m_ebCoarToFineRedist,
                                               interv, scale);

//      EBFastFR::s_verbose = false;

      m_ebFluxRegister.incrementRedistRegister(m_ebCoarToCoarRedist,
                                               interv, scale);

      if (m_doSmushing)
        {

//          if (m_level == debuglevel-1) //-1 because coarse level owns it
//            {
//              EBCoarToFineRedist::s_verbose = true;
//              EBCoarToFineRedist::s_ivdebug = ivdebamrg;
//            }
          //if use mass weighting, need to
          //fix weights of redistribution object
          if (m_useMassRedist)
            {
              int densevar = m_ebPatchGodunov->densityIndex();
              m_stateNew.exchange(interv);
              dumpDebug(string("before ctofresetw "));
              m_ebCoarToFineRedist.resetWeights(m_stateNew, densevar);
              dumpDebug(string("before ctocresetw "));
              m_ebCoarToCoarRedist.resetWeights(m_stateNew, densevar);
              dumpDebug(string("after  ctocresetw "));
            }
          //redistibute to finer level
          m_ebCoarToFineRedist.redistribute(finePtr->m_stateNew, interv);

//          EBCoarToFineRedist::s_verbose = false;

          finePtr->dumpDebug(string("after coartofine  "));

          //do the re redistirubtion
          m_ebCoarToCoarRedist.redistribute(m_stateNew, interv);

          dumpDebug(string("after coartocoar  "));
        }

      // average from finer level data

 //     if (m_level == debuglevel)
 //       {
 //         EBCoarseAverage::s_verbose = true;
 //         EBCoarseAverage::s_ivdebug = ivdebamrg;
 //       }
      finePtr->m_ebCoarseAverage.average(m_stateNew,
                                         finePtr->m_stateNew,
                                         interv);

//      EBCoarseAverage::s_verbose = false;

      dumpDebug(string("after average     "));
      m_ebCoarToFineRedist.setToZero();
      m_ebCoarToCoarRedist.setToZero();
    }
}
/****************************/
void EBAMRGodunov::tagCells(IntVectSet& a_tags)
{
  CH_TIME("EBAMRGodunov::tagCells");
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::tagCells for level " << m_level << endl;
    }

  // Create tags based on undivided gradient of density
  IntVectSet localTags;
  if (!m_tagAll)
    {
      int inoz  = 0;
      ParmParse pp;
      int queryVal = pp.query("nozzle_problem", inoz);
      bool nozzleProblem = ((queryVal != 0) && (inoz == 1));
      int itagAllIrreg;
      queryVal = pp.query("tag_all_irreg", itagAllIrreg);
      bool tagAllIrregCells = ((queryVal != 0) && (itagAllIrreg == 1));

      Real xmaxTag, ymaxTag, xminTag, yminTag;
      vector<Real> problo(SpaceDim);
      bool tagTheLot = false;
      if (nozzleProblem)
        {
          pp.get("xmax_tag", xmaxTag);
          pp.get("ymax_tag", ymaxTag);
          pp.getarr("prob_lo", problo, 0, SpaceDim);
          int queryVal = pp.query("nozzle_tag_the_lot", inoz);
          tagTheLot = ((queryVal != 0) && (inoz == 1));
          if (tagTheLot)
            {
              pp.get("xmin_tag", xminTag);
              pp.get("ymin_tag", yminTag);
            }
        }
      if (tagTheLot)
        {
          CH_assert(nozzleProblem);
          for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
            {
              const Box& b = m_grids.get(dit());
              const EBISBox& ebisBox = m_ebisl[dit()];

              // Tag all cells where in the range

              IntVectSet ivsTot(b);
              for (VoFIterator vofit(ivsTot, ebisBox.getEBGraph());
                  vofit.ok(); ++vofit)
                {
                  const VolIndex& vof = vofit();
                  const IntVect& iv = vof.gridIndex();
                  Real xval=0; //compiler is being irritating
                  Real yval=0;//compiler is being irritating
                  if (nozzleProblem)
                    {
                      xval = (iv[0] + 0.5)*m_dx[0] + problo[0];
                      yval = (iv[1] + 0.5)*m_dx[1] + problo[1];
                    }

                  if ((xval < xmaxTag) && (yval < ymaxTag) && (xval > xminTag) && (yval > yminTag))
                    {
                      localTags |= iv;
                    }
                }
            }
        }
      else //actually compute gradients and all that to compute tags
        {
          // If there is a coarser level interpolate undefined ghost cells
          //only interpolate the density
          int densityIndex = m_ebPatchGodunov->densityIndex();
          Interval intervDensity(densityIndex, densityIndex);
          EBCellFactory factory(m_ebisl);
          int nCons = m_ebPatchGodunov->numConserved();
          LevelData<EBCellFAB> consTemp(m_grids, nCons, IntVect::Unit, factory);
          Interval consInterv(0, nCons-1);
          m_stateNew.copyTo(consInterv, consTemp, consInterv);
          if (m_hasCoarser)
            {
              const EBAMRGodunov* amrGodCoarserPtr = getCoarserLevel();
              int refRatCrse = amrGodCoarserPtr->refRatio();
              int nghost = 1;
              EBPWLFillPatch patcher(m_grids,
                                     amrGodCoarserPtr->m_grids,
                                     m_ebisl,
                                     amrGodCoarserPtr->m_ebisl,
                                     amrGodCoarserPtr->m_domainBox,
                                     refRatCrse, m_nComp, nghost);

              Real coarTimeOld = 0.0;
              Real coarTimeNew = 1.0;
              Real fineTime    = 0.0;
              patcher.interpolate(consTemp,
                                  amrGodCoarserPtr->m_stateOld,
                                  amrGodCoarserPtr->m_stateNew,
                                  coarTimeOld,
                                  coarTimeNew,
                                  fineTime,
                                  intervDensity);
            }
          consTemp.exchange(intervDensity);

          // Compute undivided gradient
          for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
            {
              const Box& b = m_grids.get(dit());
              const EBISBox& ebisBox = m_ebisl[dit()];
              EBCellFAB gradFab(ebisBox, b, SpaceDim);
              const EBCellFAB& stateFab = consTemp[dit()];
              BaseFab<Real>& regGradFab = gradFab.getSingleValuedFAB();
              const BaseFab<Real>& regStateFab = stateFab.getSingleValuedFAB();

              for (int idir = 0; idir < SpaceDim; ++idir)
                {
                  const Box bCenter = b & grow(m_problem_domain,-BASISV(idir));
                  const Box bLo     = b & adjCellLo(bCenter,idir);
                  const int hasLo = ! bLo.isEmpty();
                  const Box bHi     = b & adjCellHi(bCenter,idir);
                  const int hasHi = ! bHi.isEmpty();

                  FORT_GETRELATIVEGRAD(CHF_FRA1(regGradFab,idir),
                                       CHF_CONST_FRA1(regStateFab,0),
                                       CHF_CONST_INT(idir),
                                       CHF_BOX(bLo),
                                       CHF_CONST_INT(hasLo),
                                       CHF_BOX(bHi),
                                       CHF_CONST_INT(hasHi),
                                       CHF_BOX(bCenter));

                  //do one-sided diffs where necessary at irregular cells.
                  IntVectSet ivsIrreg = ebisBox.getIrregIVS(b);
                  for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph());
                      vofit.ok(); ++vofit)
                    {
                      const VolIndex& vof = vofit();
                      const IntVect&  iv = vof.gridIndex();
                      //one-sided diffs on domain bndry
                      bool onLeftDomain = iv[idir] == m_domainBox.smallEnd(idir);
                      bool onRighDomain = iv[idir] == m_domainBox.bigEnd(idir);
                      bool hasFacesLeft = (ebisBox.numFaces(vof, idir, Side::Lo) > 0) && !onLeftDomain;
                      bool hasFacesRigh = (ebisBox.numFaces(vof, idir, Side::Hi) > 0) && !onRighDomain;

                      Real valCent = stateFab(vof, densityIndex);
                      Real dpl=0;
                      Real dpr=0;
                      Real dpc=0;
                      //compute one-sided diffs where you have them
                      if (hasFacesLeft)
                        {
                          Vector<FaceIndex> facesLeft =
                            ebisBox.getFaces(vof, idir, Side::Lo);
                          Real valLeft = 0.0;
                          for (int iface = 0; iface <facesLeft.size(); iface++)
                            {
                              VolIndex vofLeft = facesLeft[iface].getVoF(Side::Lo);
                              valLeft += stateFab(vofLeft, densityIndex);
                            }
                          valLeft /= Real(facesLeft.size());
                          dpl = valCent - valLeft;
                        }
                      if (hasFacesRigh)
                        {
                          Vector<FaceIndex> facesRigh =
                            ebisBox.getFaces(vof, idir, Side::Hi);
                          Real valRigh = 0.0;
                          for (int iface = 0; iface <facesRigh.size(); iface++)
                            {
                              VolIndex vofRigh = facesRigh[iface].getVoF(Side::Hi);
                              valRigh += stateFab(vofRigh, densityIndex);
                            }
                          valRigh /= Real(facesRigh.size());
                          dpr = valRigh - valCent;
                        }
                      if (hasFacesLeft && hasFacesRigh)
                        {
                          dpc = 0.5*(dpl+dpr);
                        }
                      else if (!hasFacesLeft && !hasFacesRigh)
                        {
                          dpc = 0.0;
                        }
                      else if (hasFacesLeft && !hasFacesRigh)
                        {
                          dpc = dpl;
                        }
                      else if (hasFacesRigh && !hasFacesLeft)
                        {
                          dpc = dpr;
                        }

                      gradFab(vof, idir) = dpc/valCent;
                    }
                }

              EBCellFAB gradMagFab(ebisBox, b, 1);
              BaseFab<Real>& regGradMagFab = gradMagFab.getSingleValuedFAB();
              FORT_MAGNITUDE(CHF_FRA1(regGradMagFab,0),
                             CHF_CONST_FRA(regGradFab),
                             CHF_BOX(b));

              //pointwise op so just have to iterate over multivalued cells
              IntVectSet ivsMulti = ebisBox.getMultiCells(b);
              for (VoFIterator vofit(ivsMulti, ebisBox.getEBGraph());
                  vofit.ok(); ++vofit)
                {
                  const VolIndex& vof = vofit();
                  Real mag = 0.0;
                  for (int idir = 0; idir < SpaceDim; idir++)
                    {
                      Real graddir = gradFab(vof, idir);
                      mag += graddir*graddir;
                    }
                  mag = sqrt(mag);
                  gradMagFab(vof, 0) = mag;
                }

              // Tag where gradient exceeds threshold

              IntVectSet ivsTot(b);
              for (VoFIterator vofit(ivsTot, ebisBox.getEBGraph());
                  vofit.ok(); ++vofit)
                {
                  const VolIndex& vof = vofit();
                  const IntVect& iv = vof.gridIndex();
                  Real xval, yval;
                  if (nozzleProblem)
                    {
                      xval = (iv[0] + 0.5)*m_dx[0] + problo[0];
                      yval = (iv[1] + 0.5)*m_dx[1] + problo[1];
                    }

                  if ((!nozzleProblem) || ((xval < xmaxTag) && (yval < ymaxTag)))
                    {
                      if (gradMagFab(vof, 0) >= m_refineThresh)
                        {
                          localTags |= iv;
                        }
                    }
                }

              //refine all irregular cells
              //this is probably not ideal but is used for replicate
              //scaling tests.
              if (tagAllIrregCells)
                {
                  IntVectSet irregIVS = ebisBox.getIrregIVS(b);
                  localTags |= irregIVS;
                }
            }
        }

      localTags.grow(m_tagBufferSize);

      // Need to do this in two steps unless a IntVectSet::operator &=
      // (ProblemDomain) operator is defined
      Box localTagsBox = localTags.minBox();
      localTagsBox &= m_problem_domain;
      localTags &= localTagsBox;
    }//if !tagall
  else
    {
      localTags.makeEmpty();
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          const Box& b = m_grids.get(dit());
          const EBISBox& ebisBox = m_ebisl[dit()];
          IntVectSet ivsTot(b);
          for (VoFIterator vofit(ivsTot, ebisBox.getEBGraph());
              vofit.ok(); ++vofit)
            {
              localTags |= vofit().gridIndex();
            }
        }
    }

  a_tags = localTags;
}
/***************************/
void EBAMRGodunov::tagCellsInit(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::tagCellsInit for level " << m_level << endl;
    }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  tagCells(a_tags);
}
/***************************/
void EBAMRGodunov::regrid(const Vector<Box>& a_new_grids)
{
  CH_TIME("EBAMRGodunov::regrid");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRGodunov regrid for level " << m_level << endl;
    }
  Vector<Box>& newGrids =   (Vector<Box>&) a_new_grids;
  {
    CH_TIME("mortonordering");
    mortonOrdering(newGrids);
  }

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  // save data for later copy
  //not using m_ebisl because it gets wiped later
  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  int nGhostEBISL = 6;
  Interval interv(0,m_nComp-1);
  LevelData<EBCellFAB> stateSaved;
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  {
    CH_TIME("defines and copies and ebisls");
    EBISLayout ebislOld;
    ebisPtr->fillEBISLayout(ebislOld, m_grids, m_domainBox, nGhostEBISL);
    EBCellFactory factoryOld(ebislOld);
    stateSaved.define(m_grids, m_nComp, ivGhost, factoryOld);
    m_stateNew.copyTo(interv, stateSaved, interv);
  }
  //create grids and ebis layouts
  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  m_level_grids = a_new_grids;
  Vector<int> proc_map;

  if (s_isLoadBalanceSet)
    {
      s_loadBalance(proc_map,a_new_grids, m_domainBox, false);
    }
  else
    {
      LoadBalance(proc_map,a_new_grids);
    }

  m_grids= DisjointBoxLayout(a_new_grids,proc_map);

  ebisPtr->fillEBISLayout(m_ebisl, m_grids, m_domainBox, nGhostEBISL);

  EBCellFactory factoryNew(m_ebisl);
  // reshape state with new grids
  m_stateNew.define(m_grids,m_nComp,ivGhost, factoryNew);
  m_stateOld.define(m_grids,m_nComp,ivGhost, factoryNew);

  // set up data structures
  levelSetup();

  // interpolate to coarser level
  if (m_hasCoarser)
    {
      EBAMRGodunov* coarPtr = getCoarserLevel();
      m_ebFineInterp.interpolate(m_stateNew,
                                 coarPtr->m_stateNew,
                                 interv);
    }

  // copy from old state
  stateSaved.copyTo(interv,m_stateNew, interv);
}
/***************************/
void EBAMRGodunov::initialGrid(const Vector<Box>& a_new_grids)
{
  CH_TIME("EBAMRG::initialGrid");

  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRGodunov initialGrid for level " << m_level << endl;
    }
  Vector<Box>& newGrids =   (Vector<Box>&) a_new_grids;
  mortonOrdering(newGrids);
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  m_level_grids = a_new_grids;

  // load balance and create boxlayout
  Vector<int> proc_map;
  if (s_isLoadBalanceSet)
    {
      s_loadBalance(proc_map,a_new_grids, m_domainBox, false);
    }
  else
    {
      LoadBalance(proc_map,a_new_grids);
    }
  if (s_verbosity >= 3)
    {
      pout() << " just loadbalanced " << m_level << endl;
    }

  m_grids = DisjointBoxLayout(a_new_grids,proc_map);
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::initialgrid grids " << endl;
      dumpDBL(&m_grids);
    }

  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  int nGhostEBISL = 6;
  if (s_verbosity >= 3)
    {
      pout() << " about to fill ebislayout  in EBAMRGodunov initialGrid for level " << m_level << endl;
    }

  ebisPtr->fillEBISLayout(m_ebisl, m_grids,  m_domainBox, nGhostEBISL);

  if (s_verbosity >= 3)
    {
      pout() << " done with filling ebislayout  in EBAMRGodunov initialGrid for level " << m_level << endl;
    }

  EBCellFactory factoryNew(m_ebisl);
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_grids,m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_grids,m_nComp, ivGhost, factoryNew);

  // set up data structures
  levelSetup();
}
/***************************/
void EBAMRGodunov::initialData()
{
  CH_TIME("EBAMRG::initialData");

  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRGodunov initialData for level " << m_level << endl;
    }

  const EBPhysIBC* const ebphysIBCPtr =
    m_ebPatchGodunov->getEBPhysIBC();

  //initialize both new and old states to
  //be the same thing
  ebphysIBCPtr->initialize(m_stateNew, m_ebisl);
  ebphysIBCPtr->initialize(m_stateOld, m_ebisl);
}
/***************************/
void EBAMRGodunov::postInitialize()
{
}
void EBAMRGodunov::syncWithFineLevel()
{
  CH_TIME("EBAMRG::syncWithFineLevel");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRGodunov syncWithFineLevel for level " << m_level << endl;
    }
  //stuff that needs to be setup from the finer
  //level.  A bunch of objects depend upon the layouts
  //from both levels and the finer level changes more
  //often from regrid so this needs to be called from the finer
  //level
  CH_assert(m_hasFiner);
  if (m_hasFiner)
    {
      EBAMRGodunov* finePtr = getFinerLevel();
      int nRefFine = refRatio();
      const DisjointBoxLayout& finer_m_grids = finePtr->m_grids;
      const EBISLayout& finer_m_ebisl = finePtr->m_ebisl;
      //define fine to coarse redistribution object
      //for now set to volume weighting
      m_ebCoarToFineRedist.define(finer_m_grids, m_grids,
                                  finer_m_ebisl,  m_ebisl,
                                  m_domainBox, nRefFine , m_nComp, 1, Chombo_EBIS::instance());
      //define coarse to coarse redistribution object
      m_ebCoarToCoarRedist.define(finer_m_grids, m_grids,
                                  finer_m_ebisl,  m_ebisl,
                                  m_domainBox, nRefFine , m_nComp);
      // maintain flux registers
      m_ebFluxRegister.define(finer_m_grids,
                              m_grids,
                              finer_m_ebisl,
                              m_ebisl,
                              m_domainBox,
                              nRefFine,
                              m_nComp, Chombo_EBIS::instance());

      //set all the registers to zero
      m_ebCoarToFineRedist.setToZero();
      m_ebCoarToCoarRedist.setToZero();
      m_ebFluxRegister.setToZero();
    }

}
/***************************/
void EBAMRGodunov::patchGodunov(const EBPatchGodunovFactory* const a_ebPatchGodunovFactory)
{
  m_ebPatchGodunovFactory = a_ebPatchGodunovFactory;
}
/***************************/
Real EBAMRGodunov::computeDt()
{
  Real newDt;
  newDt = m_dtNew;

  return newDt;
}
/***************************/
Real EBAMRGodunov::computeInitialDt()
{
  Real maxwavespeed =  m_ebLevelGodunov.getMaxWaveSpeed(m_stateNew);
  CH_assert(Abs(m_dx[0]-m_dx[1])<1.e-9);
  Real newDT = m_initial_dt_multiplier * m_dx[0] /maxwavespeed;

  return newDT;
}
/***************************/
void EBAMRGodunov::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
}
/***************************/
void EBAMRGodunov::domainLength(RealVect a_domainLength)
{
  m_domainLength = a_domainLength;
}
/***************************/
void EBAMRGodunov::refinementThreshold(Real a_refineThresh)
{
  m_refineThresh = a_refineThresh;
}
/***************************/
void EBAMRGodunov::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
}

/***************************/
void
EBAMRGodunov::sumConserved(Real& a_sumcons,
                           const int& a_ivar) const
{
  CH_TIME("EBAMRG::sumConserved");
  Real sumgrid = 0;
  for (DataIterator dit= m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& state = m_stateNew[dit()];
      Box thisBox = m_grids.get(dit());
      IntVectSet uberIVS(thisBox);
      const EBISBox& ebisBox = m_ebisl[dit()];
      for (VoFIterator vofit(uberIVS, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real consVal = state(vof, a_ivar);
          Real volFrac = ebisBox.volFrac(vof);
          Real volume = volFrac;
          if (m_doRZCoords)
            {
              Real cellVol, kvol;
              EBArith::getKVolRZ(kvol, cellVol, ebisBox, m_dx[0], vof);
              volume = cellVol*kvol;
            }
          sumgrid += consVal*volume;
        }
    }

  Vector<Real> all_sum;
  gather(all_sum,sumgrid,uniqueProc(SerialTask::compute));
  Real sumallgrid = 0.;
  if (procID() == uniqueProc(SerialTask::compute))
    {
      for (int i = 0; i < all_sum.size(); ++i)
        {
          sumallgrid += all_sum[i];
        }
    }
  broadcast(sumallgrid,uniqueProc(SerialTask::compute));
  a_sumcons = sumallgrid;
}
/***************************/
EBAMRGodunov*
EBAMRGodunov::getCoarserLevel() const
{
  EBAMRGodunov* retval = NULL;
  if (m_coarser_level_ptr != NULL)
    {
      retval = dynamic_cast <EBAMRGodunov*> (m_coarser_level_ptr);

      if (retval == NULL)
        {
          pout() << "EBAMRG::getCoarserLevel: dynamic cast failed"
                 << endl;
          MayDay::Error();
        }
    }
  return retval;
}

/***************************/
EBAMRGodunov*
EBAMRGodunov::getFinerLevel() const
{
  EBAMRGodunov* retval = NULL;
  if (m_finer_level_ptr != NULL)
    {
      retval = dynamic_cast <EBAMRGodunov*> (m_finer_level_ptr);

      if (retval == NULL)
        {
          pout() << "EBAMRG::getFinerLevel: dynamic cast failed"
                 << endl;
          MayDay::Error();
        }
    }
  return retval;
}
/***************************/
void
EBAMRGodunov::levelSetup()
{
  CH_TIME("EBAMRG::levelSetup");
  if (s_verbosity >= 3)
    {
      pout() << " in EBAMRGodunov levelSetup for level " << m_level << endl;
    }

  EBAMRGodunov* coarPtr = getCoarserLevel();
  EBAMRGodunov* finePtr = getFinerLevel();

  m_hasCoarser = (coarPtr != NULL);
  m_hasFiner   = (finePtr != NULL);

  if (m_hasCoarser)
    {
      int nRefCrse = m_coarser_level_ptr->refRatio();
      const DisjointBoxLayout& coarser_m_grids = coarPtr->m_grids;
      const EBISLayout& coarser_m_ebisl = coarPtr->m_ebisl;
      const Box& domainCoar = coarPtr->m_domainBox;

      {
        CH_TIME("ave_interp_defs");
        m_ebCoarseAverage.define(m_grids,
                                 coarser_m_grids,
                                 m_ebisl,
                                 coarser_m_ebisl,
                                 domainCoar,
                                 nRefCrse,
                                 m_nComp, Chombo_EBIS::instance());
        m_ebFineInterp.define(m_grids,
                              coarser_m_grids,
                              m_ebisl,
                              coarser_m_ebisl,
                              domainCoar,
                              nRefCrse,
                              m_nComp);
      }

      // maintain levelgodunov
      m_ebLevelGodunov.define(m_grids,
                              coarser_m_grids,
                              m_ebisl,
                              coarser_m_ebisl,
                              m_problem_domain,
                              nRefCrse,
                              m_dx,
                              m_useMassRedist,
                              m_doSmushing,
                              m_doRZCoords,
                              m_hasSourceTerm,
                              m_ebPatchGodunovFactory,
                              m_hasCoarser,
                              m_hasFiner);

      //define fine to coarse redistribution object
      //for now set to volume weighting
      {
        CH_TIME("fineToCoar_defs");
        m_ebFineToCoarRedist.define(m_grids, coarser_m_grids,
                                    m_ebisl, coarser_m_ebisl,
                                    domainCoar, nRefCrse, m_nComp);
        m_ebFineToCoarRedist.setToZero();
      }

      coarPtr->syncWithFineLevel();
    }
  else
    {
      m_ebLevelGodunov.define(m_grids,
                              DisjointBoxLayout(),
                              m_ebisl,
                              EBISLayout(),
                              m_problem_domain,
                              m_ref_ratio,
                              m_dx,
                              m_useMassRedist,
                              m_doSmushing,
                              m_doRZCoords,
                              m_hasSourceTerm,
                              m_ebPatchGodunovFactory,
                              m_hasCoarser,
                              m_hasFiner);
    }

  //set up mass redistribution array
  m_sets.define(m_grids);
  for (DataIterator dit = m_grids.dataIterator();
      dit.ok(); ++dit)
    {
      Box thisBox = m_grids.get(dit());
      m_sets[dit()] = m_ebisl[dit()].getIrregIVS(thisBox);
    }
  BaseIVFactory<Real> factory(m_ebisl, m_sets);
  //the ghost cells on the mass redistribution array
  //are tied to the redistribution radius
  int nghost = 2*m_redistRad;
  IntVect ivGhost = nghost*IntVect::Unit;
  m_massDiff.define(m_grids, m_nComp, ivGhost, factory);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      m_massDiff[dit()].setVal(0.0);
    }
}
/***************************/
LevelData<EBCellFAB>&
EBAMRGodunov::getStateNew()
{
  return m_stateNew;
}
/***************************/
LevelData<EBCellFAB>&
EBAMRGodunov::getStateOld()
{
  return m_stateOld;
}
/***************************/
Real
EBAMRGodunov::getDt() const
{
  return m_dt;
}
/***************************/
const EBISLayout&
EBAMRGodunov::getEBISLayout() const
{
  return m_ebisl;
}
/***************************/
void EBAMRGodunov::fillConsAndPrim(LevelData<EBCellFAB>& a_data) const
{
  CH_TIME("EBAMRG::fillConsAndPrim");
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::fillConsAndPrim" << endl;
    }
  //fill output data including a bunch of geometric crud
  int nCons = m_ebPatchGodunov->numConserved();
  int nPrim = m_ebPatchGodunov->numPrimitives() ;
  int consAndPrim = nCons + nPrim;

  EBCellFactory ebcellfact(m_ebisl);
  a_data.define(m_grids, consAndPrim, IntVect::Zero, ebcellfact);

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisbox = m_ebisl[dit()];
      const Box& grid = m_grids.get(dit());
      const EBCellFAB& consfab = m_stateNew[dit()];
      EBCellFAB primfab(ebisbox, grid, nPrim);
      //cfivs, time and timestep fake and not used here
      Real faket = 1.0;
      int logflag = 0;
      IntVectSet emptyivs;
      m_ebPatchGodunov->setValidBox(grid, ebisbox, emptyivs, faket, faket);
      m_ebPatchGodunov->consToPrim(primfab, consfab, grid,logflag);

      EBCellFAB& outputfab = a_data[dit()];

      Interval consIntervSrc(0, nCons-1);
      Interval consIntervDst(0, nCons-1);
      Interval primIntervSrc(0, nPrim-1);
      Interval primIntervDst(nCons, consAndPrim-1);

      // copy regular data
      outputfab.copy(grid, consIntervDst,  grid, consfab, consIntervSrc);
      outputfab.copy(grid, primIntervDst,  grid, primfab, primIntervSrc);
      Real coveredVal = -10.0;
      for (int ivar = 0; ivar < consAndPrim; ivar++)
        {
          outputfab.setInvalidData(coveredVal, ivar);
        }

    }//end loop over grids
}
/***************************/
#ifdef CH_USE_HDF5
/***************************/
void EBAMRGodunov::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::writeCheckpointHeader" << endl;
    }

  //stuff in non-eb checkpoint header is already
  //set in the define function, such as
  // Setup the number of components
  // Setup the component names
  //so i will eliminate the middleman
}
/***************************/
void EBAMRGodunov::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::readCheckpointHeader" << endl;
    }

  //stuff in non-eb checkpoint header is already
  //set in the define function, such as
  // Setup the number of components
  // Setup the component names
  // So i will eliminate the middleman.
}
/***************************/
void EBAMRGodunov::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::writeCheckpointLevel" << endl;
    }

  //reasons for deviations from non-eb stuff
  // the tag buffer size is set by the factory.
  // dx is set in the define function
  // the problem domain is set in the define function

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]       = m_ref_ratio;
  header.m_real["dt"]              = m_dt;

  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::writeCheckpointLevel 2" << endl;
    }

  // Write the header for this level
  header.writeToFile(a_handle);

  // Write the data for this level
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::writeCheckpointLevel 3" << endl;
    }
  write(a_handle,m_grids);
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::writeCheckpointLevel 4" << endl;
    }
  write(a_handle,m_stateOld,"dataOld");
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::writeCheckpointLevel 5" << endl;
    }
  write(a_handle,m_stateNew,"dataNew");
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::writeCheckpointLevel 6" << endl;
    }
}
/***************************/
/***************************/
void EBAMRGodunov::readCheckpointLevel(HDF5Handle& a_handle)
{
  CH_assert(m_isDefined);
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::readCheckpointLevel" << endl;
    }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
    {
      MayDay::Error("::readCheckpointLevel: file does not contain ref_ratio");
    }
  m_ref_ratio = header.m_int["ref_ratio"];

  //reasons for deviations from non-eb stuff
  // the tag buffer size is set by the factory.
  // dx is set in the define function
  // the problem domain is set in the define function

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("readCheckpointLevel: file does not contain dt");
    }
  m_dt = header.m_real["dt"];

  // Get the grids
  Vector<Box> vboxGrids;
  const int gridStatus = read(a_handle, vboxGrids);
  if (gridStatus != 0)
    {
      MayDay::Error("readCheckpointLevel: file has no grids");
    }

  Vector<int> proc_map;
  if (s_isLoadBalanceSet)
    {
      s_loadBalance(proc_map,vboxGrids, m_domainBox, false);
    }
  else
    {
      LoadBalance(proc_map,vboxGrids);
    }
  broadcast(proc_map, uniqueProc(SerialTask::compute));

  m_grids= DisjointBoxLayout(vboxGrids,proc_map);

  //this keeps the order of the AMRLevel m_level_grids
  //consistent with m_grids
  LayoutIterator lit = m_grids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
    {
      Box b = m_grids.get(lit());
      m_level_grids.push_back(b);
    }

  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  int nGhostEBISL = 6;
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  ebisPtr->fillEBISLayout(m_ebisl, m_grids,
                          m_domainBox, nGhostEBISL);
  EBCellFactory factoryNew(m_ebisl);
  //m_nghost is set in define function
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_grids,m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_grids,m_nComp, ivGhost, factoryNew);

  //  Interval vars(0, m_nComp-1);
  //the false says to not redefine the data
  int dataStatusNew = read<EBCellFAB>(a_handle,
                                      m_stateNew,
                                      "dataNew",
                                      m_grids,
                                      Interval(),
                                      false);

  int dataStatusOld = read<EBCellFAB>(a_handle,
                                      m_stateOld,
                                      "dataOld",
                                      m_grids,
                                      Interval(),
                                      false);

  if ((dataStatusNew != 0) || (dataStatusOld != 0))
    {
      MayDay::Error("file does not contain state data");
    }
  // Set up data structures
  levelSetup();
}

void EBAMRGodunov::writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_NewPlotFile == 0)
    {
      writePlotHeaderOld(a_handle);
      return;
    }
  Vector<int> refRatios;
  const EBAMRGodunov* current = this;
  int nlevs = 0;
  while (current != NULL)
  {
    refRatios.push_back(current->refRatio());
    nlevs++;
    current = (const EBAMRGodunov*)(current-> m_finer_level_ptr);
  }

  headerEBFile(a_handle, nlevs, refRatios,
               m_problem_domain.domainBox(),
               m_origin, m_dx, m_aspect,
               m_stateNew.ghostVect());

  writeCellCenteredNames(a_handle, m_stateNames);
  a_handle.setGroup("/Expressions");
  HDF5HeaderData expressions;
  m_ebPatchGodunov->expressions(expressions);
  expressions.writeToFile(a_handle);

}
/***************************/
void EBAMRGodunov::writePlotHeaderOld(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::writePlotHeader" << endl;
    }

  HDF5HeaderData header;
  // Setup the number of components
  //have to add in a lot of geometric crap.
  // 3 norms + 6 area fracs + 1 distance + 1 volFrac
  // =  11 extra components
  //forces 3d
  int nCons = m_ebPatchGodunov->numConserved();
  int nPrim = m_ebPatchGodunov->numPrimitives() ;
  int consAndPrim = nCons + nPrim;

  int indexVolFrac = consAndPrim;
  int indexAreaFrac = indexVolFrac+1;
  int indexNormal = indexAreaFrac+ 2*SpaceDim;
  int indexDist = indexNormal+SpaceDim;
  int nCompTotal = indexDist+1;
  CH_assert(nCompTotal == consAndPrim + 3*SpaceDim+2);

  Vector<string> names(nCompTotal);

  for (int i = 0; i < nCons; i++)
    {
      names[i] = m_stateNames[i];
    }
  for (int i = 0; i < nPrim; i++)
    {
      names[nCons + i] = m_primNames[i];
    }

  string volFracName("fraction-0");
  Vector<string> normName(3);
  Vector<string> areaName(6);
  string distName("distance-0");

  normName[0] = "xnormal-0";
  normName[1] = "ynormal-0";
  normName[2] = "znormal-0";

  areaName[0] = "xAreafractionLo-0";
  areaName[1] = "xAreafractionHi-0";
  areaName[2] = "yAreafractionLo-0";
  areaName[3] = "yAreafractionHi-0";
  areaName[4] = "zAreafractionLo-0";
  areaName[5] = "zAreafractionHi-0";

  names[indexVolFrac] = volFracName;

  for (int i=0; i < 2*SpaceDim; i++)
    {
      names[indexAreaFrac+i] = areaName[i];
    }

  for (int i=0; i < SpaceDim; i++)
    {
      names[indexNormal+i] = normName[i];
    }

  names[indexDist] = distName;

  //now output this into the hdf5 handle
  header.m_int["num_components"] = nCompTotal;
  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < nCompTotal; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = names[comp];
    }

  // Write the header to the file
  header.writeToFile(a_handle);

  if (s_verbosity >= 4)
    {
      pout() << header << endl;
    }
}

void EBAMRGodunov::writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_NewPlotFile == 0)
    {
      writePlotLevelOld(a_handle);
      return;
    }
  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::writePlotLevel " << m_level<< endl;
    }
  writeCellCentered(a_handle, m_level, &m_stateNew);
}

/***************************/
void EBAMRGodunov::writePlotLevelOld(HDF5Handle& a_handle) const
{

  if (s_verbosity >= 3)
    {
      pout() << "EBAMRGodunov::writePlotLevel" << endl;
    }
  //fill output data including a bunch of geometric crud
  int nCons = m_ebPatchGodunov->numConserved();
  int nPrim = m_ebPatchGodunov->numPrimitives() ;
  int consAndPrim = nCons + nPrim;
  int indexVolFrac = consAndPrim;
  int indexAreaFrac = indexVolFrac+1;
  int indexNormal = indexAreaFrac+ 2*SpaceDim;
  int indexDist = indexNormal+SpaceDim;
  int nCompTotal = indexDist+1;
  CH_assert(nCompTotal == consAndPrim + 3*SpaceDim+2);

  Vector<Real> coveredValuesCons(nCons, -10.0);
  Vector<Real> coveredValuesPrim(nPrim, -10.0);

  Vector<Real> coveredValues;
  coveredValues.append(coveredValuesCons);
  coveredValues.append(coveredValuesPrim);

  LevelData<FArrayBox> fabData(m_grids, nCompTotal, IntVect::Zero);

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisbox = m_ebisl[dit()];
      const Box& grid = m_grids.get(dit());
      const EBCellFAB& consfab = m_stateNew[dit()];
      EBCellFAB primfab(ebisbox, grid, nPrim);
      //cfivs, time and timestep fake and not used here
      Real faket = 1.0;
      IntVectSet emptyivs;
      ParmParse pp;
      int logflag = 1;
      if (pp.contains("logflag"))
        {
          pp.get("logflag", logflag);
        }
      m_ebPatchGodunov->setValidBox(grid, ebisbox, emptyivs, faket, faket);
      m_ebPatchGodunov->consToPrim(primfab, consfab, grid, logflag);

      FArrayBox& currentFab = fabData[dit()];

      // copy regular data
      currentFab.copy(consfab.getSingleValuedFAB(),0,0,nCons);
      currentFab.copy(primfab.getSingleValuedFAB(),0,nCons,nPrim);

      // set default volume fraction
      currentFab.setVal(1.0,indexVolFrac);

      // set default area fractions
      for (int i=0; i < 2*SpaceDim; i++)
        {
          currentFab.setVal(1.0,indexAreaFrac+i);
        }

      // set default normal
      for (int i=0; i < SpaceDim; i++)
        {
          currentFab.setVal(0.0,indexNormal+i);
        }

      // set default distance of EB from corner
      currentFab.setVal(0.0,indexDist);

      // set special values
      // iterate through the current grid
      // NOTE:  this is probably an inefficient way to do this
      for (BoxIterator bit(grid); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          // set special values for covered cells
          if (ebisbox.isCovered(iv))
            {
              for (int icomp = 0; icomp < consAndPrim; icomp++)
                {
                  Real cval = coveredValues[icomp];

                  currentFab(iv,icomp) = cval;
                }
              // volume fraction is zero
              currentFab(iv,indexVolFrac) = 0.0;

              // area fractions are zero
              for (int i=0; i < 2*SpaceDim; i++)
                {
                  currentFab(iv,indexAreaFrac+i) = 0.0;
                }
            }

          // set special values for irregular cells
          if (ebisbox.isIrregular(iv))
            {
              Vector<VolIndex> vofs = ebisbox.getVoFs(iv);
              Real volFrac = ebisbox.volFrac(vofs[0]);
              RealVect normal = ebisbox.normal(vofs[0]);

              // set volume fraction
              currentFab(iv,indexVolFrac) = volFrac;

              // set area fractions--use only the first face you find
              for (int i=0; i < SpaceDim; i++)
                {
                  Vector<FaceIndex> faces;

                  faces = ebisbox.getFaces(vofs[0],i,Side::Lo);
                  if (faces.size() == 0)
                    {
                      currentFab(iv,indexAreaFrac+2*i) = 0.0;
                    }
                  else
                    {
                      currentFab(iv,indexAreaFrac+2*i) =
                        ebisbox.areaFrac(faces[0]);
                    }

                  faces = ebisbox.getFaces(vofs[0],i,Side::Hi);
                  if (faces.size() == 0)
                    {
                      currentFab(iv,indexAreaFrac+2*i+1) = 0.0;
                    }
                  else
                    {
                      currentFab(iv,indexAreaFrac+2*i+1) =
                        ebisbox.areaFrac(faces[0]);
                    }
                }

              // set normal
              for (int i=0; i < SpaceDim; i++)
                {
                  currentFab(iv,indexNormal+i) = normal[i];
                }

              // set distance unless the length of the normal is zero
              Real length = PolyGeom::dot(normal,normal);

              if (length > 0)
                {
                  Real dist = PolyGeom::computeAlpha(volFrac,normal)*m_dx[0];
                  currentFab(iv,indexDist) = -dist;
                }
            } //end if (isIrregular)
        }//end loop over cells
    }//end loop over grids

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx[0];
  header.m_real["dt"]          = m_dt;
  header.m_real["time"]        = m_time;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 4)
    {
      pout() << header << endl;
    }

  // Write the data for this level
  write(a_handle,fabData.boxLayout());
  write(a_handle,fabData,"data");
}

#endif
#include "NamespaceFooter.H"
