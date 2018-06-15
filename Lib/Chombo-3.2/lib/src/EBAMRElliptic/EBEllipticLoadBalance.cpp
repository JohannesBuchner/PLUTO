#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBEllipticLoadBalance.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBISBox.H"
#include "LoadBalance.H"
#include "EBPoissonOp.H"
#include "EBPoissonOpFactory.H"
#include "SPMD.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "EBLevelDataOps.H"
#include "NeumannPoissonDomainBC.H"
#include "NeumannPoissonEBBC.H"
#include "TimedDataIterator.H"
#include "NamespaceHeader.H"

#define EBELB_NUM_APPLY_OPS 100
///////////////
void getPoissonLoadsAndBoxes(Vector<unsigned long long>&            a_loads,
                             Vector<Box> &            a_boxes,
                             const DisjointBoxLayout& a_dblOrig,
                             const ProblemDomain&     a_domain,
                             const EBIndexSpace*      a_ebisPtr )
{
  int nghost = 4;
  IntVect nghostPhi = 4*IntVect::Unit;
  EBISLayout ebisl;
  a_ebisPtr->fillEBISLayout(ebisl, a_dblOrig, a_domain, nghost);

  //and we need to remap
  TimedDataIterator dit = a_dblOrig.timedDataIterator();
  dit.clearTime();
  dit.enableTime();
  //scoping bracket
  EBCellFactory ebcellfact(ebisl);
  LevelData<EBCellFAB> phi(a_dblOrig, 1, nghostPhi, ebcellfact);
  LevelData<EBCellFAB> lph(a_dblOrig, 1, nghostPhi, ebcellfact);

  EBLevelDataOps::setToZero(lph);
  EBLevelDataOps::setToZero(phi);

  RefCountedPtr<BaseDomainBCFactory> domBC = RefCountedPtr<BaseDomainBCFactory>(new NeumannPoissonDomainBCFactory());
  RefCountedPtr<BaseEBBCFactory>      ebBC = RefCountedPtr<BaseEBBCFactory>(    new NeumannPoissonEBBCFactory());

  EBLevelGrid eblg(a_dblOrig, ebisl, a_domain);

  EBPoissonOpFactory factory(eblg, RealVect::Unit, RealVect::Zero, 2, 1, 1,domBC, ebBC, 1.0, 1.0, nghostPhi, nghostPhi);

  RefCountedPtr<EBPoissonOp>  ebpo = RefCountedPtr<EBPoissonOp>(factory.MGnewOp(a_domain, 0, false));

  //evaluate poisson operator---homogeneous bcs so i don't have to set the value
  for (int iapply = 0; iapply < EBELB_NUM_APPLY_OPS; iapply++)
    {
      ebpo->applyOp(lph, phi, true, dit, false);
    }
  Vector<Box>  boxesLocal = dit.getBoxes();
  Vector<unsigned long long> loadsLocal = dit.getTime();

  dit.disableTime();

  dit.mergeTime();
  a_loads = dit.getTime();
  a_boxes = dit.getBoxes();

}
///////////////
void
resetLoadOrder(Vector<unsigned long long>&    a_loads,
               Vector<Box>&     a_newBoxes,
               Vector<Box>&     a_oldBoxes)
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
          MayDay::Error("Something broken in EBEllipticLoadBalance");
        }
    }
}
///////////////
int
EBEllipticLoadBalance(Vector<int>&         a_procs,
                      const Vector<Box>&   a_boxes,
                      const ProblemDomain& a_domain,
                      bool a_verbose,
                      const EBIndexSpace*  a_ebis_ptr )
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
  getPoissonLoadsAndBoxes(loads, boxes, dblOrig, a_domain, a_ebis_ptr );

  resetLoadOrder(loads, boxes, inBoxes);
  if (a_verbose)
    {
      pout() << "EBEllipticLoadBalance loads:" << endl;
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
