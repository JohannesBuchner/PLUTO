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
#include "EBAlias.H"
#include "UsingNamespace.H"
#include "slab.cpp"

/***************/
/***************/
int checkLevelAlias(const Box& a_domain);

/***************/
/***************/
int checkAMRAlias(const Box& a_domain);

/***************/
/***************/
void aliasEBAMR(Vector<LevelData<FArrayBox>* >& a_output,
                Vector<LevelData<EBCellFAB>* >& a_input);
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
    makeGeometry(domain, dx, coveredDomain);

    eekflag = checkLevelAlias(domain);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in checkLevelAlias");
      }

    eekflag = checkAMRAlias(domain);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in checkAMRAlias");
      }

    pout() << "testEBAlias test passed" << endl;
  }//end scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}


/***************/
/***************/
int checkLevelAlias(const Box& a_domain)
{
  int eekflag = 0;
  int nvar = 1;
  Interval interv(0,0);
  int nghost= 0;
  IntVect ivgh= IntVect::Zero;

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



      EBCellFactory factOne(ebislOne);
      EBCellFactory factTwo(ebislTwo);

      LevelData<EBCellFAB> dataOne(dblOne, nvar, ivgh, factOne);
      LevelData<EBCellFAB> dataTwo(dblTwo, nvar, ivgh, factTwo);

      for (DataIterator dit = dblOne.dataIterator(); dit.ok(); ++dit)
        dataOne[dit].setVal(1);
      for (DataIterator dit = dblTwo.dataIterator(); dit.ok(); ++dit)
        dataTwo[dit].setVal(2);

      {
        LevelData<FArrayBox> reg1, reg2;
        aliasEB(reg1, dataOne);
        aliasEB(reg2, dataTwo);
      } // test for alias going out of scope not destroying original data.

      LevelData<FArrayBox> reg1, reg2;
      aliasEB(reg1, dataOne);
      aliasEB(reg2, dataTwo);

      for (DataIterator dit = dblOne.dataIterator(); dit.ok(); ++dit)
        {
          if (reg1[dit].dataPtr() != dataOne[dit].getFArrayBox().dataPtr())
            {
              pout()<<" error for data pointer "<<dblOne[dit]<<std::endl;
              eekflag ++;
            }

        }
      levelDomain.coarsen(2);

    }
  return eekflag;
}
/***************/
/***************/
int checkAMRAlias(const Box& a_domain)
{
  int eekflag = 0;
  int nvar = 1;
  Interval interv(0,0);
  int nghost= 0;
  IntVect ivgh= IntVect::Zero;

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  int numlevels = ebisPtr->numLevels();

  Vector<LevelData<EBCellFAB>* > dataOne(numlevels-1,NULL),dataTwo(numlevels-1,NULL);
  Vector<DisjointBoxLayout> dblOne(numlevels-1),dblTwo(numlevels-1);

  Box levelDomain = a_domain;
  for (int ilev = 0; ilev < numlevels-1; ilev++)
    {
      int maxsizeOne = 4;
      int maxsizeTwo = 2;
      makeLayout(dblOne[ilev], levelDomain, maxsizeOne);
      makeLayout(dblTwo[ilev], levelDomain, maxsizeTwo);

      EBISLayout ebislOne, ebislTwo;
      makeEBISL(ebislOne, dblOne[ilev], levelDomain, nghost);
      makeEBISL(ebislTwo, dblTwo[ilev], levelDomain, nghost);

      EBCellFactory factOne(ebislOne);
      EBCellFactory factTwo(ebislTwo);

      dataOne[ilev] = new LevelData<EBCellFAB>(dblOne[ilev], nvar, ivgh, factOne);
      dataTwo[ilev] = new LevelData<EBCellFAB>(dblTwo[ilev], nvar, ivgh, factTwo);
      levelDomain.coarsen(2);
    }

  for (int ilev = 0; ilev < numlevels-1; ilev++)
    {
      for (DataIterator dit = dblOne[ilev].dataIterator(); dit.ok(); ++dit)
        (*dataOne[ilev])[dit].setVal(1);
      for (DataIterator dit = dblTwo[ilev].dataIterator(); dit.ok(); ++dit)
        (*dataTwo[ilev])[dit].setVal(2);

      {
        LevelData<FArrayBox> reg1, reg2;
        aliasEB(reg1, *dataOne[ilev]);
        aliasEB(reg2, *dataTwo[ilev]);
      } // test for alias going out of scope not destroying original data.

      LevelData<FArrayBox> reg1, reg2;
      aliasEB(reg1, *dataOne[ilev]);
      aliasEB(reg2, *dataTwo[ilev]);

      for (DataIterator dit = dblOne[ilev].dataIterator(); dit.ok(); ++dit)
        {
          if (reg1[dit].dataPtr() != (*dataOne[ilev])[dit].getFArrayBox().dataPtr())
            {
              pout()<<" error for data pointer "<<dblOne[ilev][dit]<<std::endl;
              eekflag ++;
            }

        }
    }


  Vector<LevelData<FArrayBox>* > dataOneRegPtr(numlevels-1,NULL);
  for (int i=0; i<numlevels-1; i++)
    dataOneRegPtr[i] = new LevelData<FArrayBox>();

  aliasEBAMR(dataOneRegPtr,dataOne);

  for (int i=0; i<numlevels-1; i++)
    delete dataOneRegPtr[i] ;

  return eekflag;
}
void aliasEBAMR(Vector<LevelData<FArrayBox>* >& a_output,
                Vector<LevelData<EBCellFAB>* >& a_input)
{
  for (int ilev = 0; ilev < a_input.size(); ilev++)
    {
      //test 1
      LevelData<FArrayBox> foo;
      aliasEB(foo, *a_input[ilev]);
      //test 2
      CH_assert(a_output[ilev] != NULL);
      aliasEB(*a_output[ilev], *a_input[ilev]);
    }
}
