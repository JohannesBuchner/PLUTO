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

#include "EBAMRTestCommon.H"
#include "EBTensorCFInterp.H"
#include "EBQuadCFInterp.H"

#include "UsingNamespace.H"
IntVect ivdebug(D_DECL(14,16,10));
/***************/
/**********/
int 
getError(LevelData<EBCellFAB>&    a_errorVelo,
         BaseAnalytic        &    a_exactFunc,
         const DisjointBoxLayout& a_gridsFine,
         const EBISLayout&        a_ebislFine,
         const Box&               a_domainFine,
         const RealVect&          a_dxFine,
         const DisjointBoxLayout& a_gridsCoar,
         const EBISLayout&        a_ebislCoar,
         const Box&               a_domainCoar,
         const RealVect&          a_dxCoar,
         const int & a_refrat)
{
  int eekflag = 0;
  EBCellFactory factFine(a_ebislFine);
  EBCellFactory factCoar(a_ebislCoar);
  int nvarVelo = SpaceDim;
  LevelData<EBCellFAB> veloFine(a_gridsFine, nvarVelo, IntVect::Unit, factFine);
  LevelData<EBCellFAB> veloCoar(a_gridsCoar, nvarVelo, IntVect::Unit, factCoar);
  //fill internal cells with the right answer (coarse and fine)
  for (DataIterator dit = a_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      IntVectSet ivsBox(a_gridsFine.get(dit()));
      for (VoFIterator vofit(ivsBox, a_ebislFine[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          RealVect location = EBArith::getVofLocation(vof, a_dxFine, RealVect::Zero);
          Real veloVal[SpaceDim]; 
          a_exactFunc.veloVal(veloVal, location);
          for (int ivar = 0; ivar < SpaceDim; ivar++)
            {
              veloFine[dit()](vof, ivar) = veloVal[ivar];
            }
        }
    }
  for (DataIterator dit = a_gridsCoar.dataIterator(); dit.ok(); ++dit)
    {
      IntVectSet ivsBox(a_gridsCoar.get(dit()));
      for (VoFIterator vofit(ivsBox, a_ebislCoar[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          RealVect location = EBArith::getVofLocation(vof, a_dxCoar, RealVect::Zero);
          Real veloVal[SpaceDim]; 
          a_exactFunc.veloVal(veloVal, location);
          for (int ivar = 0; ivar < SpaceDim; ivar++)
            {
              veloCoar[dit()](vof, ivar) = veloVal[ivar];
            }
        }
    }

  //now do a coarse-fine interpolation to fill ghost cells
  LayoutData<IntVectSet> cfivs;
  EBArith::defineCFIVS(cfivs, a_gridsFine, a_domainFine);
  EBQuadCFInterp interpOp(a_gridsFine,  a_gridsCoar,
                          a_ebislFine,  a_ebislCoar,
                          a_domainCoar, a_refrat, SpaceDim,
                          cfivs, Chombo_EBIS::instance(), true);


  Interval interv(0,SpaceDim-1);
  interpOp.interpolate(veloFine, veloCoar, interv);
  veloFine.exchange();
  //now compare answer everywhere (including ghos cells) with exact answer
  //will only be different in the ghost cells

  for (DataIterator dit = a_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      Box grownBox = a_gridsFine.get(dit());
      grownBox.grow(1);
      grownBox &= a_domainFine;
      IntVectSet ivsBox(grownBox);
      for (VoFIterator vofit(ivsBox, a_ebislFine[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          int ihere = 0;
          if (vof.gridIndex() == ivdebug)
            {
              ihere = 1;
            }
          RealVect location = EBArith::getVofLocation(vof, a_dxFine, RealVect::Zero);
          Real veloVal[SpaceDim]; 
          a_exactFunc.veloVal(veloVal, location);
          for (int ivar = 0; ivar < SpaceDim; ivar++)
            {
              Real veloIn = veloFine[dit()](vof, ivar);
              a_errorVelo[dit()](vof,ivar) = veloIn - veloVal[ivar];
            }
        }
    }
  return eekflag;
}
/***************/
int exactTest(const DisjointBoxLayout& a_gridsFine,
              const EBISLayout&        a_ebislFine,
              const Box&               a_domainFine,
              const RealVect&          a_dxFine,
              const DisjointBoxLayout& a_gridsCoar,
              const EBISLayout&        a_ebislCoar,
              const Box&               a_domainCoar,
              const RealVect&          a_dxCoar,
              int a_refrat, int a_halfDir)
{
  int eekflag = 0;
  int nvSolu = SpaceDim;
  EBCellFactory fact(a_ebislFine);
  IntVect ivghost = IntVect::Unit;
  LevelData<EBCellFAB> errorVelo(a_gridsFine, nvSolu, ivghost, fact);
  ConstAnalytic   constfunc;
  LinearAnalytic linearfunc;
  QuadAnalytic     quadfunc;


  pout() << "testing constant function" << endl; 
  eekflag =  getError(errorVelo,  constfunc,
                      a_gridsFine, a_ebislFine,   a_domainFine, a_dxFine,
                      a_gridsCoar, a_ebislCoar,   a_domainCoar, a_dxCoar,  a_refrat);
  if (eekflag != 0)
    {
      pout() << "exactTest: problem in geterror for constfunc" << endl;
      return eekflag;
    }
  eekflag =  EBAMRTestCommon::checkForZero(errorVelo,  a_gridsFine, a_ebislFine, a_domainFine, string("constant_vel"));
  if (eekflag != 0)
    {
      pout() << "exactTest: problem in geterror checkforzero for const" << endl;
      return -3;
    }
    
  pout() << "testing linear function" << endl; 
  eekflag =  getError(errorVelo,  linearfunc,
                      a_gridsFine, a_ebislFine,   a_domainFine, a_dxFine,
                      a_gridsCoar, a_ebislCoar,   a_domainCoar, a_dxCoar,  a_refrat);
  if (eekflag != 0)
    {
      pout() << "exactTest: problem in geterror for linearfunc" << endl;
      return eekflag;
    }

  eekflag =  EBAMRTestCommon::checkForZero(errorVelo,  a_gridsFine, a_ebislFine, a_domainFine, string("linear_vel"));
  if (eekflag != 0)
    {
      pout() << "exactTest: problem in geterror  checkforzero linearfunc" << endl;
      return -4;
    }

  pout() << "testing quadratic function" << endl; 
  eekflag =  getError(errorVelo,  quadfunc,
                      a_gridsFine, a_ebislFine,   a_domainFine, a_dxFine,
                      a_gridsCoar, a_ebislCoar,   a_domainCoar, a_dxCoar,  a_refrat);
  if (eekflag != 0)
    {
      pout() << "exactTest: problem in geterror for quad" << endl;
      return eekflag;
    }
  eekflag =  EBAMRTestCommon::checkForZero(errorVelo,  a_gridsFine, a_ebislFine, a_domainFine, string("constant_vel"));
  if (eekflag != 0)
    {
      pout() << "exactTest: problem in geterror checkforzero for quad" << endl;
      return -5;
    }
    
  return eekflag;
}
/***************/
int convergenceTest(const DisjointBoxLayout& a_gridsFine,
                    const EBISLayout&        a_ebislFine,
                    const Box&               a_domainFine,
                    const RealVect&          a_dxFine,
                    const DisjointBoxLayout& a_gridsCoar,
                    const EBISLayout&        a_ebislCoar,
                    const Box&               a_domainCoar,
                    const RealVect&          a_dxCoar,
                    int a_refrat, int a_halfDir
                    )
{
  int eekflag = 0;
  DisjointBoxLayout gridsCoFi, gridsCoCo;
  coarsen(gridsCoFi, a_gridsFine, 2);
  coarsen(gridsCoCo, a_gridsCoar, 2);
  Box domainCoFi = coarsen(a_domainFine, 2);
  Box domainCoCo = coarsen(a_domainCoar, 2);
  EBISLayout ebislCoFi, ebislCoCo;
  int nghost = 2;
  eekflag = EBAMRTestCommon::makeEBISL(ebislCoFi, gridsCoFi, domainCoFi, nghost);
  eekflag = EBAMRTestCommon::makeEBISL(ebislCoCo, gridsCoCo, domainCoCo, nghost);
  RealVect dxCoFi = 2*a_dxFine;
  RealVect dxCoCo = a_refrat*dxCoFi;
 
  int nvSolu = SpaceDim;
  EBCellFactory ebcellfactFine(a_ebislFine);
  EBCellFactory ebcellfactCoar(  ebislCoFi);
  IntVect ivghost = IntVect::Unit;

  LevelData<EBCellFAB> errorVeloFine(a_gridsFine, nvSolu, ivghost, ebcellfactFine);
  LevelData<EBCellFAB> errorVeloCoFi(  gridsCoFi, nvSolu, ivghost, ebcellfactCoar);
    
  CubeAnalytic cubefunc;
  eekflag =  getError(errorVeloFine,  cubefunc,
                      a_gridsFine, a_ebislFine,   a_domainFine, a_dxFine,
                      a_gridsCoar, a_ebislCoar,   a_domainCoar, a_dxCoar,  a_refrat);

  if (eekflag != 0) return eekflag;

  eekflag =  getError(errorVeloCoFi,  cubefunc,
                      gridsCoFi, ebislCoFi,   domainCoFi, dxCoFi,
                      gridsCoCo, ebislCoCo,   domainCoCo, dxCoCo,  a_refrat);

  if (eekflag != 0) return eekflag;

  pout() << "comparing velocity error"   << endl;
  eekflag =  EBAMRTestCommon::compareError(errorVeloFine, a_ebislFine, a_gridsFine, a_domainFine,
                                           errorVeloCoFi,   ebislCoFi,   gridsCoFi,   domainCoFi);
  if (eekflag != 0) return eekflag;


  char charVeloFine[100];
  char charVeloCoar[100];
  sprintf(charVeloFine,"bhqFineError.%dhalf.%dd.hdf5",a_halfDir,SpaceDim);
  sprintf(charVeloCoar,"bhqCoarError.%dhalf.%dd.hdf5",a_halfDir,SpaceDim);
  string fileVeloFine(charVeloFine);
  string fileVeloCoar(charVeloCoar);


  ParmParse pp;
  int output_error;
  pp.get("output_error", output_error);
  if (output_error == 1)
    {
      EBAMRTestCommon::outputError(errorVeloFine,
                                   errorVeloCoFi,
                                   a_gridsFine,
                                   gridsCoFi,
                                   a_domainFine,
                                   domainCoFi,
                                   fileVeloFine,
                                   fileVeloCoar);
    }
  return eekflag;
}
/***************/
int
main(int argc, char** argv)
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int eekflag = 0;
  {//begin forever present scoping trick

    const char* in_file = "halftensor.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    //define the geometry object.
    //need the new and delete because of
    //strong construction
    //make the gometry.  this makes the first Chombo_EBIS
    Box domainFine, domainCoar;
    RealVect dxFine, dxCoar;
    EBDebugPoint::s_ivd = ivdebug;
    //and defines it using a geometryservice
    eekflag =  EBAMRTestCommon::makeGeometry(domainFine,  dxFine);
    CH_assert(eekflag == 0);

    int refRat = 2;
    domainCoar = coarsen(domainFine, refRat);
    dxCoar = refRat*dxFine;
    //make grids
    for (int halfDir = 0; halfDir < SpaceDim; halfDir++)
      {
        pout() << " cutting grids along direction = " << halfDir << endl;
        DisjointBoxLayout gridsFine, gridsCoar;
        EBISLayout        ebislFine, ebislCoar;
        eekflag = EBAMRTestCommon::makeHalfLayouts(gridsFine, ebislFine, 
                                                   gridsCoar, ebislCoar,
                                                   domainFine, domainCoar, refRat, halfDir);
        if (eekflag != 0)
          {
            pout() << "problem in makelayout"  <<endl;
            return eekflag;
          }

        eekflag = exactTest(gridsFine,
                            ebislFine,
                            domainFine,
                            dxFine,
                            gridsCoar,
                            ebislCoar,
                            domainCoar,
                            dxCoar,
                            refRat,
                            halfDir);

        if (eekflag != 0)
          {
            pout() << "problem in exactTest"  <<endl;
            return eekflag;
          }

        eekflag = convergenceTest(gridsFine,
                                  ebislFine,
                                  domainFine,
                                  dxFine,
                                  gridsCoar,
                                  ebislCoar,
                                  domainCoar,
                                  dxCoar,
                                  refRat, halfDir);

        if (eekflag != 0)
          {
            pout() << "problem in convergence"  <<endl;
            return eekflag;
          }
      }


  }//end omnipresent scoping trick
  CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  pout() << "bigHalfQuadTest passed" << endl;
  return eekflag;
}
/***************/
