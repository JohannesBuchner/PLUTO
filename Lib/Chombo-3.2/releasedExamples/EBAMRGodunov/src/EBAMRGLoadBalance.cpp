#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBAMRGLoadBalance.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBISBox.H"
#include "LoadBalance.H"
#include "SPMD.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "EBLevelDataOps.H"
#include "EBPatchPolytropic.H"
#include "TimedDataIterator.H"
#include "EBExplosionIBCFactory.H"

#include "NamespaceHeader.H"

#define EBAMRG_NUM_APPLY_OPS 1
///////////////
void getEBAMRGLoadsAndBoxes(Vector<unsigned long long>&            a_loads,
                            Vector<Box> &                          a_boxes,
                            const DisjointBoxLayout&               a_dblOrig,
                            const ProblemDomain&                   a_domain)
{
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  int nghost = 4;
  EBISLayout ebisl;
  ebisPtr->fillEBISLayout(ebisl, a_dblOrig, a_domain, nghost);

  //and we need to remap
  TimedDataIterator dit = a_dblOrig.timedDataIterator();
  dit.clearTime();
  dit.enableTime();
  Real gamma = 1.4;
  Real dx = 1.0e-1;
  Real dt = 0.001;
  Real dumm = 1; //just a dummy number
  IntVectSet cfivs;
  EBExplosionIBCFactory bcfact(gamma, dumm, dumm, dumm, dumm, dumm, dumm*RealVect::Unit, false);
  EBPatchPolytropic patch;
  patch.define(a_domain, dx);
  patch.setGamma(gamma);
  patch.setEBPhysIBC(bcfact);
  patch.artificialViscosity(true);
  patch.setSlopeParameters(true, true, true);
  int ncomp =patch.numConserved();
  Real maxwave;
  for (int iapply = 0; iapply < EBAMRG_NUM_APPLY_OPS; iapply++)
    {
      for (dit.reset(); dit.ok(); ++dit)
        {
          Box grownBox = grow(a_dblOrig.get(dit()), 4);
          grownBox    &= a_domain;

          IntVectSet ivsIrreg =  ebisl[dit()].getIrregIVS(a_dblOrig.get(dit()));

          EBCellFAB    consstate(ebisl[dit()], grownBox, ncomp);
          EBCellFAB   flattening(ebisl[dit()], grownBox, 1);
          EBCellFAB       source;
          EBFluxFAB         flux(ebisl[dit()], a_dblOrig.get(dit()), ncomp);
          BaseIVFAB<Real> ebirregflux(ivsIrreg, ebisl[dit()].getEBGraph(), ncomp);
          BaseIVFAB<Real>  nonconsdiv(ivsIrreg, ebisl[dit()].getEBGraph(), ncomp);
          BaseIVFAB<Real>    massdiff(ivsIrreg, ebisl[dit()].getEBGraph(), ncomp);
          BaseIFFAB<Real> centroidflux[SpaceDim];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              centroidflux[idir].define(ivsIrreg, ebisl[dit()].getEBGraph(), idir, ncomp);
              centroidflux[idir].setVal(0.);
            }
          consstate.setVal(1.0);
          flattening.setVal(1.0);
          flux.setVal(0.);
          ebirregflux.setVal(0.);
          nonconsdiv.setVal(0.);
          massdiff.setVal(0.);
          patch.setValidBox(a_dblOrig.get(dit()), ebisl[dit()], cfivs, dumm, dt);

          patch.regularUpdate(consstate, flux, ebirregflux, nonconsdiv,
                              flattening, source,  a_dblOrig.get(dit()),
                              ivsIrreg, dit(), false);

          patch.irregularUpdate(consstate, maxwave, massdiff, centroidflux,
                                ebirregflux, nonconsdiv, a_dblOrig.get(dit()), ivsIrreg);
        }
    }
  dit.disableTime();
  dit.mergeTime();
  a_loads = dit.getTime();
  a_boxes = dit.getBoxes();

}
///////////////
void
resetLoadOrderEBAMRG(Vector<unsigned long long>&   a_loads,
                     Vector<Box>&                  a_newBoxes,
                     Vector<Box>&                  a_oldBoxes)
{
  //resets order from a_newBoxes order to a_oldBoxes order
  Vector<unsigned long long> newLoads = a_loads;
  CH_assert(a_loads.size() == a_newBoxes.size());
  CH_assert(a_loads.size() == a_oldBoxes.size());

  //now remap the loads to the original vector of boxes
  for (int ibox = 0; ibox < a_oldBoxes.size(); ibox++)
    {
      bool foundTheBox  = false;
      for (int jbox = 0; jbox < a_newBoxes.size(); jbox++)
        {
          if (a_oldBoxes[ibox] == a_newBoxes[jbox])
            {
              foundTheBox = true;
              a_loads[ibox] = newLoads[jbox];
            }
        }
      if (!foundTheBox)
        {
          MayDay::Error("Something broken in EBAMRGLoadBalance");
        }
    }
}
///////////////
int
EBAMRGLoadBalance(Vector<int>&         a_procs,
                  const Vector<Box>&   a_boxes,
                  const ProblemDomain& a_domain,
                  bool a_verbose)
{
#ifndef CH_MPI
  a_procs.resize(a_boxes.size(),0);
  int retval=0;
#else
  //first load balance the conventional way.
  Vector<Box> inBoxes = a_boxes;
  Vector<int> origProcs;
  int retval=LoadBalance(origProcs, inBoxes);
  a_procs = origProcs;

  //we shall make fully covered boxes = constant load = covered load
  //we shall say that irregular points get irregular factor more load than
  //regular points.
  //compute a load for each box.   by evaluating the poisson operator
  //use one scratch space for everything
  DisjointBoxLayout dblOrig(inBoxes, origProcs, a_domain);

  Vector<unsigned long long> loads;
  Vector<Box>  boxes;
  getEBAMRGLoadsAndBoxes(loads, boxes, dblOrig, a_domain);

  resetLoadOrderEBAMRG(loads, boxes, inBoxes);
  if (a_verbose)
    {
      pout() << "EBAMRGLoadBalance loads:" << endl;
      for (int ibox = 0; ibox < a_boxes.size(); ibox++)
        {
          pout()
            << "   box[" << ibox << "]=" <<   a_boxes[ibox]
            << ", load[" << ibox << "]=" <<     loads[ibox] << endl;
        }
    }
  //do the load balance with our EB load estimates and the original a_boxes vector
  retval = UnLongLongLoadBalance(a_procs,  loads, a_boxes);
#endif

  return retval;
}

#include "NamespaceFooter.H"
