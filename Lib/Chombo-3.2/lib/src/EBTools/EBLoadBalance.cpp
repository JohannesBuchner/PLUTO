#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBLoadBalance.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBISBox.H"
#include "LoadBalance.H"
#include "SPMD.H"
#include "NamespaceHeader.H"

int
EBLoadBalance(Vector<int>&         a_procs,
              const Vector<Box>&   a_boxes,
              const ProblemDomain& a_domain,
              int                  a_coveredLoad,
              int                  a_irregFactor,
              bool                 a_verbose,
              const EBIndexSpace*  a_ebisPtr )
{
  //first load balance the conventional way.
  Vector<Box> inBoxes = a_boxes;
  Vector<Box> outBoxes = inBoxes;
  Vector<int> origProcs;
  int retval=LoadBalance(origProcs, inBoxes);
  a_procs = origProcs;

  /**/
  if (numProc() > 1)
    {
      //we shall make fully covered boxes = constant load = covered load
      //we shall say that irregular points get irregular factor more load than
      //regular points.
      int nghost = 0;
      DisjointBoxLayout dblOrig(inBoxes, origProcs, a_domain);
      EBISLayout ebisl;
      a_ebisPtr->fillEBISLayout(ebisl, dblOrig, a_domain, nghost);

      Vector<Box>  boxesLocal;
      Vector<int>  loadsLocal;
      //compute a load for each box.  save each box in a vector too
      //because this process will screw up the box ordering
      //and we need to remap
      for (DataIterator dit = dblOrig.dataIterator(); dit.ok(); ++dit)
        {
          const Box& box = dblOrig.get(dit());
          const EBISBox& ebisBox = ebisl[dit()];
          int npts = box.numPts();
          int boxload;

          if (ebisBox.isAllRegular())
            {
              boxload = npts;
            }
          else if (ebisBox.isAllCovered())
            {
              boxload = a_coveredLoad;
            }
          else
            {
              IntVectSet ivsIrreg = ebisBox.getIrregIVS(box);
              int numIrreg = ivsIrreg.numPts();

              int irregLoad = numIrreg*a_irregFactor;
              boxload = irregLoad + npts;
            }
        boxesLocal.push_back(box);
        loadsLocal.push_back(boxload);

        }

      //gather loads and boxes into all-encompassing vectors
      Vector<Vector<Box> >   allBoxes;
      Vector<Vector<int> >   allLoads;
      gather(allBoxes, boxesLocal, uniqueProc(SerialTask::compute));
      gather(allLoads, loadsLocal, uniqueProc(SerialTask::compute));

      Vector<int> loads;
      outBoxes.resize(0);
      if (procID() == uniqueProc(SerialTask::compute))
        {
          if (allBoxes.size() != allLoads.size())
            {
              MayDay::Error("mismatch between allboxes and all loads in ebloadbalance");
            }
          for (int iproc = 0; iproc < allBoxes.size(); iproc++)
            {
              loads.append(allLoads[iproc]);
              outBoxes.append(allBoxes[iproc]);
            }
        }

      //broadcast these all-encompassing vectors to all procs
      //and call loadbalance again with our computed loads
      broadcast(outBoxes, uniqueProc(SerialTask::compute));
      broadcast(loads,    uniqueProc(SerialTask::compute));

      //for some reason, loadbalance wants longs for its loads.
      Vector<long> longLoads(loads.size());

      //now remap the loads to the original vector of boxes
      for (int ibox = 0; ibox < a_boxes.size(); ibox++)
        {
          if (a_boxes.size() != longLoads.size())
            {
              pout() << "a_boxes.size()="  << a_boxes.size()  << endl;
              pout() << "longLoads.size()="<< longLoads.size()<< endl;
              pout() << "outBoxes.size()=" << outBoxes.size() << endl;
              pout() << "allBoxes.size()=" << allBoxes.size() << endl;
              pout() << "allLoads.size()=" << allLoads.size() << endl;
              MayDay::Error("mismatch between a_boxes and longLoads  in ebloadbalance");
            }
          if (a_boxes.size() != outBoxes.size())
            {
              pout() << "a_boxes.size()="  << a_boxes.size()  << endl;
              pout() << "longLoads.size()="<< longLoads.size()<< endl;
              pout() << "outBoxes.size()=" << outBoxes.size() << endl;
              pout() << "allBoxes.size()=" << allBoxes.size() << endl;
              pout() << "allLoads.size()=" << allLoads.size() << endl;
              MayDay::Error("mismatch between a_boxes and longLoads  in ebloadbalance");
            }
          bool foundTheBox  = false;
          for (int jbox = 0; jbox < outBoxes.size(); jbox++)
            {
              if (a_boxes[ibox] == outBoxes[jbox])
                {
                  foundTheBox = true;
                  longLoads[ibox] = loads[jbox];
                }
            }
          if (!foundTheBox)
            {
              pout() << "ibox = "<< ibox << ", a_boxes[ibox] =" << a_boxes[ibox] << "not found in output vector" << endl;
              MayDay::Error("Something broken in EBLoadBalance");
            }
        }

      //do the load balance with our EB load estimates and the original a_boxes vector
      retval = LoadBalance(a_procs, longLoads, a_boxes);
      if (a_verbose)
        {
          pout() << "EBLoadBalance loads:" << endl;
          for (int ibox = 0; ibox < a_boxes.size(); ibox++)
            {
              pout()
                << "   box[" << ibox << "]=" <<   a_boxes[ibox]
                << ", load[" << ibox << "]=" << longLoads[ibox] << endl;
            }
        }
    } //end if (numProc > 1)
  /**/
  return retval;
}
#include "NamespaceFooter.H"
