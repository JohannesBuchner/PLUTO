#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BaseIFFAB.H"
#include "BaseIFFactory.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "VoFIterator.H"
#include "SlabService.H"
#include "PolyGeom.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBCellFactory.H"
#include "BaseIVFactory.H"
#include "CH_Attach.H"
#include "UsingNamespace.H"

const bool g_verbose = true;
/**********/
int makeGeometry(Box& a_domain)
{
  IntVect loiv = IntVect::Zero;
  IntVect hiiv = 31*IntVect::Unit;
  a_domain = Box(loiv, hiiv);

  IntVect loCovered = 5*IntVect::Unit;
  IntVect hiCovered = 5*IntVect::Unit;
  loCovered[0] = 0;
  hiCovered[0] = 31;
  Box coveredBox(loCovered, hiCovered);
  SlabService slab(coveredBox);

  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  RealVect origin = RealVect::Zero;
  Real dx = 1.0/32.0;
  ebisPtr->define(a_domain, origin, dx, slab);

  return 0;
}
/***************/
Real getExactValue(const VolIndex& vof)
{
  const IntVect& iv = vof.gridIndex();
  const int&     ci = vof.cellIndex();
  Real retval = 3.*Real(ci);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      retval += Real((idir+1)*iv[idir]);
    }
  return retval;
}
/***/
int
checkRegIO(const Box&               a_domain,
           const DisjointBoxLayout& a_dbl,
           const EBISLayout&        a_ebisl)
{
#ifdef CH_USE_HDF5
  EBCellFactory fact(a_ebisl);
  LevelData<EBCellFAB> outputLD(a_dbl, 1, IntVect::Zero, fact);
  LevelData<EBCellFAB>  inputLD(a_dbl, 1, IntVect::Zero, fact);
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = a_ebisl[dit()];
      EBCellFAB& outputFAB = outputLD[dit()];
      outputFAB.setCoveredCellVal(-7.77, 0);
      IntVectSet ivsGrid(a_dbl.get(dit()));
      for (VoFIterator vofit(ivsGrid, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real exactVal = getExactValue(vof);
          outputFAB(vof, 0) = exactVal;
        }
    }

  const char* const  filename= "regfabioTest.hdf5";
  HDF5Handle handleOut(filename, HDF5Handle::CREATE);
  write(handleOut, outputLD, "fileIOTestData");
  handleOut.close();

  HDF5Handle handleIn(filename, HDF5Handle::OPEN_RDONLY);

  //the false says to not redefine the data
  int dataStatusNew = read<EBCellFAB>(handleIn,
                                      inputLD,
                                      "fileIOTestData",
                                      a_dbl,
                                      Interval(),
                                      false);

  handleIn.close();
  if (dataStatusNew != 0)
    {
      pout() << "check regio: error on input" << endl;
      return dataStatusNew;

    }

  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = a_ebisl[dit()];
      const EBCellFAB& inputFAB = inputLD[dit()];
      IntVectSet ivsGrid(a_dbl.get(dit()));
      for (VoFIterator vofit(ivsGrid, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real exactVal = getExactValue(vof);
          Real inputVal = inputFAB(vof, 0);
          if (Abs(exactVal -inputVal) > PolyGeom::getTolerance())
            {
              pout() << "check regular IO error" << endl;
              pout() << "values do not match at vof = " << vof << endl;
              pout() << " exact val = " << exactVal << endl;
              pout() << " input val = " << inputVal << endl;
              return -1;
            }
        }
    }
#endif
  return 0;
}
/********/
int
checkIrregIO(const Box&               a_domain,
             const DisjointBoxLayout& a_dbl,
             const EBISLayout&        a_ebisl)
{
#ifdef CH_USE_HDF5
  LayoutData<IntVectSet> sets(a_dbl);
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = a_ebisl[dit()];
      const Box& box = a_dbl.get(dit());
      sets[dit()] = ebisBox.getIrregIVS(box);
    }
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      BaseIFFactory<Real> fact(a_ebisl, sets, idir);
      LevelData<BaseIFFAB<Real> > outputLD(a_dbl, 1, IntVect::Zero, fact);
      LevelData<BaseIFFAB<Real> >  inputLD(a_dbl, 1, IntVect::Zero, fact);
      FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;
      for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBox     =  a_ebisl[dit()];
          BaseIFFAB<Real>& outputFAB = outputLD[dit()];
          const IntVectSet& ivs      =     sets[dit()];
          for (FaceIterator faceit(ivs, ebisBox.getEBGraph(), idir, stopCrit);
              faceit.ok(); ++faceit)
            {
              const VolIndex& vof = faceit().getVoF(Side::Lo);
              Real exactVal = getExactValue(vof);
              outputFAB(faceit(), 0) = exactVal;
            }
        }

      char filename[80];
      sprintf(filename, "basiffabio%d.hdf5", idir);
      HDF5Handle handleOut(filename, HDF5Handle::CREATE);
      write(handleOut, outputLD, "fileIOTestData");
      handleOut.close();

      HDF5Handle handleIn(filename, HDF5Handle::OPEN_RDONLY);

      //the false says to not redefine the data
      int dataStatusNew = read<BaseIFFAB<Real> >(handleIn,
                                                 inputLD,
                                                 "fileIOTestData",
                                                 a_dbl,
                                                 Interval(),
                                                 false);

      if (dataStatusNew != 0)
        {
          pout() << "check face irrio: error on input" << endl;
          return dataStatusNew;

        }
      handleIn.close();

      for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBox          = a_ebisl[dit()];
          const BaseIFFAB<Real>& inputFAB = inputLD[dit()];
          const IntVectSet& ivs           =    sets[dit()];

          for (FaceIterator faceit(ivs, ebisBox.getEBGraph(), idir, stopCrit);
              faceit.ok(); ++faceit)
            {
              const VolIndex& vof = faceit().getVoF(Side::Lo);
              Real exactVal = getExactValue(vof);
              Real inputVal = inputFAB(faceit(), 0);
              if (Abs(exactVal -inputVal) > PolyGeom::getTolerance())
                {
                  pout() << "check irreg IO error" << endl;
                  pout() << "values do not match at face = " << vof << endl;
                  pout() << " exact val = " << exactVal << endl;
                  pout() << " input val = " << inputVal << endl;
                  return -2;
                }
            }
        }
    }

  BaseIVFactory<Real> fact(a_ebisl, sets);
  LevelData<BaseIVFAB<Real> > outputLD(a_dbl, 1, IntVect::Zero, fact);
  LevelData<BaseIVFAB<Real> >  inputLD(a_dbl, 1, IntVect::Zero, fact);
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox     =  a_ebisl[dit()];
      BaseIVFAB<Real>& outputFAB = outputLD[dit()];
      const IntVectSet& ivs      =     sets[dit()];
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real exactVal = getExactValue(vof);
          outputFAB(vof, 0) = exactVal;
        }
    }

  const char* const  filename= "irrfabioTest.hdf5";
  HDF5Handle handleOut(filename, HDF5Handle::CREATE);
  write(handleOut, outputLD, "fileIOTestData");
  handleOut.close();

  HDF5Handle handleIn(filename, HDF5Handle::OPEN_RDONLY);

  //the false says to not redefine the data
  int dataStatusNew = read<BaseIVFAB<Real> >(handleIn,
                                             inputLD,
                                             "fileIOTestData",
                                             a_dbl,
                                             Interval(),
                                             false);

  if (dataStatusNew != 0)
    {
      pout() << "check vof irrio: error on input" << endl;
      return dataStatusNew;

    }
  handleIn.close();

  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox          = a_ebisl[dit()];
      const BaseIVFAB<Real>& inputFAB = inputLD[dit()];
      const IntVectSet& ivs           =    sets[dit()];

      for (VoFIterator vofit(ivs, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real exactVal = getExactValue(vof);
          Real inputVal = inputFAB(vof, 0);
          if (Abs(exactVal -inputVal) > PolyGeom::getTolerance())
            {
              pout() << "check irreg IO error" << endl;
              pout() << "values do not match at vof = " << vof << endl;
              pout() << " exact val = " << exactVal << endl;
              pout() << " input val = " << inputVal << endl;
              return -3;
            }
        }
    }
#endif
  return 0;
}
/***************/
int checkIO(const Box& a_domain)
{
#ifdef CH_USE_HDF5
  Box levelDomain = a_domain;
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int numLevels = ebisPtr->numLevels();
  int eekflag = 0;
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      CH_assert(!levelDomain.isEmpty());
      Vector<Box> vbox(1, levelDomain);
      Vector<int> procAssign(1, 0);
      DisjointBoxLayout levelDBL(vbox, procAssign);
      EBISLayout levelEBISL;
      int nghost= 0;
      ebisPtr->fillEBISLayout(levelEBISL, levelDBL,
                              levelDomain, nghost);
      eekflag = checkIrregIO(levelDomain, levelDBL, levelEBISL);
      if (eekflag != 0)
        {
          pout() << "problem in irregular IO" << ilev << endl;
          return eekflag;
        }
      eekflag = checkRegIO(levelDomain, levelDBL, levelEBISL);
      if (eekflag != 0)
        {
          pout() << "problem in regular IO" << ilev << endl;
          return eekflag;
        }
     levelDomain.coarsen(2);
    }
#endif
  return 0;
}
/********/
int
main(int argc, char** argv)
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //begin forever present scoping trick
  {
    // registerDebugger();
    int eekflag = 0;
    //define the geometry object.
    Box domain;
    eekflag =  makeGeometry(domain);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeGeometry");
      }

    //check that total coarse vol = total fine vol
    //check that total coarse vol centroid = total fine vol centroid
    eekflag = checkIO(domain);
    if (eekflag != 0)
      {
        pout() << "checkIO: eek = " << eekflag << endl;
        MayDay::Error("fabIOTest: problem in checkIO");
      }


  }//end scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  pout() << "fabIO test passed" << endl;
  return 0;
}
/***************/
/***************/
/***************/
