#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <string>
#include "EBDebugDump.H"
#include "DebugDump.H"
#include "PolyGeom.H"
#include "DataIterator.H"
#include "SPMD.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "LayoutIterator.H"
#include "ParmParse.H"
#include "parstream.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "SlabService.H"
#include "BoxIterator.H"
#include "BaseIVFAB.H"
#include "BaseIFFAB.H"
#include "BaseEBCellFAB.H"
#include "BaseEBFaceFAB.H"
#include "VoFIterator.H"
#include "LevelData.H"
#include "EBCellFactory.H"
#include "BaseIVFactory.H"
#include "UsingNamespace.H"
#include "slab.cpp"


/***************/
/***************/
int checkIrregFabCopy(const Box& a_domain);
/***************/
/***************/
int checkRegularFabCopy(const Box& a_domain);
/***************/
/***************/
int
main(int argc, char** argv)
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //begin scoping trick
  {

    int eekflag = 0;
    //define the geometry object.
    //need the new and delete because of
    //strong construction
    Box domain, coveredDomain;
    Real dx;
    //make the gometry.  this makes the first Chombo_EBIS
    //and defines it using a geometryservice
    eekflag =  makeGeometry(domain, dx, coveredDomain);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeGeometry");
      }

    eekflag = checkIrregFabCopy(domain);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in checkIrregFabCopy");
      }

    eekflag = checkRegularFabCopy(domain);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in checkRegularFabCopy");
      }

    pout() << "fab copy test passed" << endl;
  }//end scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}

/***************/
int checkIrregFabCopy(const Box& a_domain)
{
  int eekflag = 0;
  int nvar = 1;
  Interval interv(0,0);
  int nghost= 0;
  IntVect ivgh= IntVect::Zero;
  Real tolerance = PolyGeom::getTolerance();
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  int numlevels = ebisPtr->numLevels();

  Box levelDomain = a_domain;
  for (int ilev = 0; ilev < numlevels-1; ilev++)
    {
      DisjointBoxLayout dblOne, dblTwo;
      int maxsizeOne = 4;
      int maxsizeTwo = 2;
      makeLayout(dblOne, levelDomain, maxsizeOne);
      makeLayout(dblTwo, levelDomain, maxsizeTwo);

      EBISLayout ebislOne, ebislTwo;
      makeEBISL(ebislOne, dblOne, levelDomain, nghost);
      makeEBISL(ebislTwo, dblTwo, levelDomain, nghost);

      LayoutData<IntVectSet> ldIVSOne(dblOne);
      LayoutData<IntVectSet> ldIVSTwo(dblTwo);
      for (DataIterator dit = dblOne.dataIterator(); dit.ok(); ++dit)
        ldIVSOne[dit()] = IntVectSet(dblOne.get(dit()));
      for (DataIterator dit = dblTwo.dataIterator(); dit.ok(); ++dit)
        ldIVSTwo[dit()] = IntVectSet(dblTwo.get(dit()));

      BaseIVFactory<Real> factOne(ebislOne, ldIVSOne);
      BaseIVFactory<Real> factTwo(ebislTwo, ldIVSTwo);

      LevelData<BaseIVFAB<Real> > dataOne(dblOne, nvar, ivgh, factOne);
      LevelData<BaseIVFAB<Real> > dataTwo(dblTwo, nvar, ivgh, factTwo);

      for (DataIterator dit = dblOne.dataIterator(); dit.ok(); ++dit)
        dataOne[dit()].setVal(1);
      for (DataIterator dit = dblTwo.dataIterator(); dit.ok(); ++dit)
        dataTwo[dit()].setVal(2);

      //copy dataone into datatwo.  then all datatwo should hold are ones
      dataOne.copyTo(interv, dataTwo, interv);
      Real rightVal = 1.0;
      for (DataIterator dit = dblTwo.dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBox = ebislTwo[dit()];
          const IntVectSet&  ivs = ldIVSTwo[dit()];
          const BaseIVFAB<Real>& fab = dataTwo[dit()];
          for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              Real thisVal = fab(vofit(), 0);
              if (Abs(thisVal - rightVal) > tolerance)
                {
                  eekflag = 3;
                  return eekflag;
                }
            }
        }
      levelDomain.coarsen(2);
    }
  return eekflag;
}

/***************/
/***************/
int checkRegularFabCopy(const Box& a_domain)
{
  int eekflag = 0;
  int nvar = 1;
  Interval interv(0,0);
  int nghost= 0;
  IntVect ivgh= IntVect::Zero;
  Real tolerance = PolyGeom::getTolerance();
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  int numlevels = ebisPtr->numLevels();

  Box levelDomain = a_domain;
  for (int ilev = 0; ilev < numlevels-1; ilev++)
    {
      DisjointBoxLayout dblOne, dblTwo;
      int maxsizeOne = 4;
      int maxsizeTwo = 2;
      makeLayout(dblOne, levelDomain, maxsizeOne);
      makeLayout(dblTwo, levelDomain, maxsizeTwo);

      EBISLayout ebislOne, ebislTwo;
      makeEBISL(ebislOne, dblOne, levelDomain, nghost);
      makeEBISL(ebislTwo, dblTwo, levelDomain, nghost);

      LayoutData<IntVectSet> ldIVSOne(dblOne);
      LayoutData<IntVectSet> ldIVSTwo(dblTwo);
      for (DataIterator dit = dblOne.dataIterator(); dit.ok(); ++dit)
        ldIVSOne[dit()] = IntVectSet(dblOne.get(dit()));
      for (DataIterator dit = dblTwo.dataIterator(); dit.ok(); ++dit)
        ldIVSTwo[dit()] = IntVectSet(dblTwo.get(dit()));

      EBCellFactory factOne(ebislOne);
      EBCellFactory factTwo(ebislTwo);

      LevelData<EBCellFAB> dataOne(dblOne, nvar, ivgh, factOne);
      LevelData<EBCellFAB> dataTwo(dblTwo, nvar, ivgh, factTwo);

      for (DataIterator dit = dblOne.dataIterator(); dit.ok(); ++dit)
        dataOne[dit()].setVal(1);
      for (DataIterator dit = dblTwo.dataIterator(); dit.ok(); ++dit)
        dataTwo[dit()].setVal(2);

      //copy dataone into datatwo.  then all datatwo should hold are ones
      dataOne.copyTo(interv, dataTwo, interv);
      Real rightVal = 1.0;
      for (DataIterator dit = dblTwo.dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBox = ebislTwo[dit()];
          const IntVectSet&  ivs = ldIVSTwo[dit()];
          const EBCellFAB&   fab = dataTwo[dit()];
          for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              Real thisVal = fab(vofit(), 0);
              if (Abs(thisVal - rightVal) > tolerance)
                {
                  eekflag = 3;
                  return eekflag;
                }
            }
        }
      levelDomain.coarsen(2);

    }
  return eekflag;
}
