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
IntVect ivdebug(D_DECL(16,10,8));
/**********/
int 
getError(LevelData<EBCellFAB>&    a_errorVelo,
         LevelData<EBCellFAB>&    a_errorGrad,
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
  int nvarGrad = SpaceDim*SpaceDim;
  IntVect ivghost = 4*IntVect::Unit;
  LevelData<EBCellFAB> veloFine(a_gridsFine, nvarVelo, ivghost, factFine);
  LevelData<EBCellFAB> veloCoar(a_gridsCoar, nvarVelo, ivghost, factCoar);
  LevelData<EBCellFAB> gradFine(a_gridsFine, nvarGrad, ivghost, factFine);
  //fill internal cells with the right answer (coarse and fine)
  for (DataIterator dit = a_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      IntVectSet ivsBox(a_gridsFine.get(dit()));
      for (VoFIterator vofit(ivsBox, a_ebislFine[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          RealVect location = EBArith::getVofLocation(vof, a_dxFine, RealVect::Zero);
          Real veloVal[SpaceDim]; 
          Real gradVal[SpaceDim*SpaceDim]; 
          a_exactFunc.veloVal(veloVal, location);
          a_exactFunc.gradVal(gradVal, location);
          for (int ivar = 0; ivar < SpaceDim; ivar++)
            {
              veloFine[dit()](vof, ivar) = veloVal[ivar];
              for (int jdir = 0; jdir < SpaceDim; jdir++)
                {
                  int iindex = TensorCFInterp::gradIndex(ivar, jdir);
                  gradFine[dit()](vof, iindex) = gradVal[iindex];
                }
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
  EBTensorCFInterp interpOp(a_gridsFine,  a_gridsCoar,
                            a_ebislFine,  a_ebislCoar,
                            a_domainCoar, a_refrat, SpaceDim, a_dxFine[0],
                            cfivs, Chombo_EBIS::instance(), true);


  Interval interv(0,SpaceDim-1);
  interpOp.coarseFineInterp(veloFine, gradFine, veloCoar);
  veloFine.exchange();
  gradFine.exchange();
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
          Real gradVal[SpaceDim*SpaceDim]; 
          a_exactFunc.veloVal(veloVal, location);
          a_exactFunc.gradVal(gradVal, location);
          for (int ivar = 0; ivar < SpaceDim; ivar++)
            {
              Real veloIn = veloFine[dit()](vof, ivar);
              a_errorVelo[dit()](vof,ivar) = veloIn - veloVal[ivar];
              for (int jdir = 0; jdir < SpaceDim; jdir++)
                {
                  int iindex = TensorCFInterp::gradIndex(ivar, jdir);
                  a_errorGrad[dit()](vof, iindex) = 0;
                  //only tangential gradients are used
                  if (ivar != jdir)
                    {
                      Real gradIn = gradFine[dit()](vof, iindex);
                      a_errorGrad[dit()](vof, iindex) = gradIn - gradVal[iindex];
                    }
                }
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
              int a_refrat)
{
  int eekflag = 0;
  int nvSolu = SpaceDim;
  int nvGrad = SpaceDim*SpaceDim;
  EBCellFactory fact(a_ebislFine);
  IntVect ivghost = 4*IntVect::Unit;
  LevelData<EBCellFAB> errorVelo(a_gridsFine, nvSolu, ivghost, fact);
  LevelData<EBCellFAB> errorGrad(a_gridsFine, nvGrad, ivghost, fact);
  ConstAnalytic   constfunc;
  LinearAnalytic linearfunc;
  QuadAnalytic     quadfunc;
  pout() << "checking constant func" << endl;
  eekflag =  getError(errorVelo, errorGrad, constfunc,
                      a_gridsFine, a_ebislFine,   a_domainFine, a_dxFine,
                      a_gridsCoar, a_ebislCoar,   a_domainCoar, a_dxCoar,  a_refrat);
  if (eekflag != 0)
    {
      pout() << "exactTest: problem in geterror for constfunc" << endl;
      return eekflag;
    }
  eekflag = EBAMRTestCommon::checkForZero(errorVelo, a_gridsFine, a_ebislFine, a_domainFine, string("constant_vel"));
  if (eekflag != 0)
    {
      pout() << "exactTest: problem in solu checkforzero for const" << endl;
      return -31;
    }

  eekflag = EBAMRTestCommon::checkForZero(errorGrad, a_gridsFine, a_ebislFine, a_domainFine, string("constant_vel"));
  if (eekflag != 0)
    {
      pout() << "exactTest: problem in grad checkforzero for const" << endl;
      return -3;
    }

  pout() << "checking linear func" << endl;
  eekflag =  getError(errorVelo, errorGrad, linearfunc,
                      a_gridsFine, a_ebislFine,   a_domainFine, a_dxFine,
                      a_gridsCoar, a_ebislCoar,   a_domainCoar, a_dxCoar,  a_refrat);
  if (eekflag != 0)
    {
      pout() << "exactTest: problem in geterror for linearfunc" << endl;
      return eekflag;
    }

  eekflag = EBAMRTestCommon::checkForZero(errorVelo, a_gridsFine, a_ebislFine, a_domainFine, string("linear_vel"));
  if (eekflag != 0)
    {
      pout() << "exactTest: problem in solu  checkforzero linearfunc" << endl;
      return -4;
    }

  eekflag = EBAMRTestCommon::checkForZero(errorGrad, a_gridsFine, a_ebislFine, a_domainFine, string("linear_vel"));
  if (eekflag != 0)
    {
      pout() << "exactTest: problem in grad checkforzero linearfunc" << endl;
      return -41;
    }

  pout() << "checking quadratic func" << endl;
  eekflag =  getError(errorVelo, errorGrad, quadfunc,
                      a_gridsFine, a_ebislFine,   a_domainFine, a_dxFine,
                      a_gridsCoar, a_ebislCoar,   a_domainCoar, a_dxCoar,  a_refrat);
  if (eekflag != 0)
    {
      pout() << "exactTest: problem in geterror for quadfunc" << endl;
      return eekflag;
    }

  eekflag = EBAMRTestCommon::checkForZero(errorVelo, a_gridsFine, a_ebislFine, a_domainFine, string("linear_vel"));
  if (eekflag != 0)
    {
      pout() << "exactTest: problem in solu   checkforzero quadfunc" << endl;
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
  int nghost = 4;
  eekflag = EBAMRTestCommon::makeEBISL(ebislCoFi, gridsCoFi, domainCoFi, nghost);
  eekflag = EBAMRTestCommon::makeEBISL(ebislCoCo, gridsCoCo, domainCoCo, nghost);
  RealVect dxCoFi = 2*a_dxFine;
  RealVect dxCoCo = a_refrat*dxCoFi;
 
  int nvSolu = SpaceDim;
  int nvGrad = SpaceDim*SpaceDim;
  EBCellFactory ebcellfactFine(a_ebislFine);
  EBCellFactory ebcellfactCoar(  ebislCoFi);
  IntVect ivghost = 4*IntVect::Unit;

  LevelData<EBCellFAB> errorVeloFine(a_gridsFine, nvSolu, ivghost, ebcellfactFine);
  LevelData<EBCellFAB> errorVeloCoFi(  gridsCoFi, nvSolu, ivghost, ebcellfactCoar);
                                                                   
  LevelData<EBCellFAB> errorGradFine(a_gridsFine, nvGrad, ivghost, ebcellfactFine);
  LevelData<EBCellFAB> errorGradCoFi(  gridsCoFi, nvGrad, ivghost, ebcellfactCoar);
    
  CubeAnalytic cubefunc;
  eekflag =  getError(errorVeloFine, errorGradFine, cubefunc,
                      a_gridsFine, a_ebislFine,   a_domainFine, a_dxFine,
                      a_gridsCoar, a_ebislCoar,   a_domainCoar, a_dxCoar,  a_refrat);

  if (eekflag != 0) return eekflag;

  eekflag =  getError(errorVeloCoFi, errorGradCoFi, cubefunc,
                      gridsCoFi, ebislCoFi,   domainCoFi, dxCoFi,
                      gridsCoCo, ebislCoCo,   domainCoCo, dxCoCo,  a_refrat);

  if (eekflag != 0) return eekflag;

  pout() << "comparing velocity error"   << endl;
  eekflag =  EBAMRTestCommon::compareError(errorVeloFine, a_ebislFine, a_gridsFine, a_domainFine,
                                           errorVeloCoFi,   ebislCoFi,   gridsCoFi,   domainCoFi);
  if (eekflag != 0) return eekflag;

  pout() << "comparing gradient error"   << endl;
  eekflag =  EBAMRTestCommon::compareError(errorGradFine, a_ebislFine, a_gridsFine, a_domainFine,
                                           errorGradCoFi,   ebislCoFi,   gridsCoFi,   domainCoFi);
  if (eekflag != 0) return eekflag;

  char charVeloFine[100];
  char charVeloCoar[100];
  char charGradFine[100];
  char charGradCoar[100];
  sprintf(charVeloFine, "bhtFineVeloError.%dhalf.%dd.hdf5",a_halfDir, SpaceDim);
  sprintf(charVeloCoar, "bhtCoarVeloError.%dhalf.%dd.hdf5",a_halfDir, SpaceDim);
  sprintf(charGradFine, "bhtFineGradError.%dhalf.%dd.hdf5",a_halfDir, SpaceDim);
  sprintf(charGradCoar, "bhtCoarGradError.%dhalf.%dd.hdf5",a_halfDir, SpaceDim);
  string fileVeloFine(charVeloFine);
  string fileVeloCoar(charVeloCoar);
  string fileGradFine(charGradFine);
  string fileGradCoar(charGradCoar);

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

      EBAMRTestCommon::outputError(errorGradFine,
                                   errorGradCoFi,
                                   a_gridsFine,
                                   gridsCoFi,
                                   a_domainFine,
                                   domainCoFi,
                                   fileGradFine,
                                   fileGradCoar);
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
    for (int halfDir = 0; halfDir < SpaceDim; halfDir++)
      {
        pout() << " cutting grids along direction = " << halfDir << endl;
        //make grids
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
                            refRat);

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
                                  refRat, 
                                  halfDir);

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
  pout() << "halfTensorTest passed" << endl;
  return eekflag;
}
/***************/
