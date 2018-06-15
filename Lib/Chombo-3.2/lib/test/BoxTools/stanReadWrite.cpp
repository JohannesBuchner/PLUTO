#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//#if defined(NDEBUG)
//int  main()
//{
//  return -1;
//}
//#else

#include <cmath>
#include <algorithm>

#include "CH_HDF5.H"
#include "BRMeshRefine.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LoadBalance.H"
#include "DebugDump.H"
#include "UGIO.H"
#include "TestCommon.H"
#include "UsingNamespace.H"


/// Global variables for handling output:
static const char *pgmname = "stanReadWrite" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

int getVal(const IntVect& a_iv, int a_icomp)
{
  int retval = a_iv[0]*a_iv[0] + a_iv[1]*a_iv[1];
  retval *= (a_icomp+1);
  return retval;
}
///////
void
setData(LevelData<BaseFab<int> >& a_data)
{
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      for (BoxIterator bit(a_data[dit()].box()); bit.ok(); ++bit)
        {
          for (int icomp = 0; icomp < a_data.nComp(); icomp++)
            {
              a_data[dit()](bit(), icomp) = getVal(bit(), icomp);
            }
        }
    }
}
int
compareData(LevelData<BaseFab<int> >& a_data1, LevelData<BaseFab<int> >& a_data2)
{
  DataIterator dit2 = a_data2.dataIterator();
  DataIterator dit1 = a_data1.dataIterator();
  for (dit1.reset(), dit2.reset(); dit1.ok() && dit2.ok(); ++dit1, ++ dit2)
    {
      for (BoxIterator bit(a_data1[dit1()].box()); bit.ok(); ++bit)
        {
          for (int icomp = 0; icomp < a_data1.nComp(); icomp++)
            {
              if (a_data1[dit1()](bit(), icomp) != a_data2[dit2()](bit(), icomp))
                {
                  return -1;
                }
            }
        }
    }
  return 0;
}
///////
int
stanTest()
{
  //  int     domsiz = 8;
  int domsiz = 32;
  int     ncomps  = SpaceDim;
  IntVect nghost  = 3*IntVect::Unit;
  Box     domain(IntVect::Zero, (domsiz-1)*IntVect::Unit);
  int     cirrad = domsiz/4;
  IntVect center = (domsiz/2)*IntVect::Unit;
  int     thickn = 1;
  IntVectSet tags;
  setCircleTags(tags, cirrad, thickn, center);

  //DisjointBoxLayout grids(Vector<Box>(1, domain), Vector<int>(1, 0));
  DisjointBoxLayout grids;
  buildDisjointBoxLayout(grids, tags, domain);

  //create data
  LevelData<BaseFab<int> > dataOut(grids, ncomps, nghost);
  setData(dataOut);

  //for time, dx, dt all that
  Real dum = 1.;
  //sensible values for refrat and lev but i am not sure they matter
  int lev = 0;    int ref = 2;
  Interval interv(0, ncomps-1);
  //output data to file
  string filename("standfile.hdf5");
  HDF5Handle handleOut(filename.c_str(),  HDF5Handle::CREATE);

  writeLevel(handleOut, lev, dataOut, dum, dum, dum, domain, ref, nghost, interv);
  handleOut.close();


  //read data back in
  LevelData<BaseFab<int> >  dataIn(grids, ncomps, nghost);
  HDF5Handle handleIn(filename.c_str(),  HDF5Handle::OPEN_RDONLY);

  readLevel( handleIn,  lev, dataIn,  dum, dum, dum, domain, ref,  interv, false);
  handleIn.close();

  //compare the data
  int eekflag = compareData(dataOut, dataIn);
  return eekflag;
}

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int icode = stanTest();
  if (icode != 0)
    {
      pout() << indent << pgmname <<" failed"<<endl;
    }
  else
    {
      pout() << indent << pgmname <<" passed"<<endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return icode;
}
















//#endif
