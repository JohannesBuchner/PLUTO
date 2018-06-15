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
#include <iostream>

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#include "BaseIVFactory.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PlaneIF.H"
#include "IntersectionIF.H"
#include "PolyGeom.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBPatchGodunovF_F.H"
#include "EBPatchAdvect.H"
#include "UsingNamespace.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

/***************/
/***************/
void
makeGeometry(Box& a_domain,
             Real& a_dx)
{
  //parse input file.  single level
  ParmParse pp;
  RealVect origin = RealVect::Zero;
  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

 CH_assert(n_cell.size() == SpaceDim);
  IntVect lo = IntVect::Zero;
  IntVect hi;
  for (int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if (n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          MayDay::Error();
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_domain = Box(lo, hi);

  Real prob_hi;
  int numOpen = n_cell[0];
  pp.get("domain_length",prob_hi);
  a_dx = prob_hi/numOpen;

  pout() << "channel geometry" << endl;
  RealVect channelNormal;
  vector<Real>  channelNormalVect(SpaceDim);
  pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      channelNormal[idir] = channelNormalVect[idir];
    }

  Real channelRadius;
  pp.get("channel_radius", channelRadius);

  Real norm = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      norm += channelNormal[idir] * channelNormal[idir];
    }
  norm = sqrt(norm);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      channelNormal[idir] /= norm;
    }

  RealVect botPoint = RealVect::Zero;
  botPoint -= channelNormal;
  botPoint *= channelRadius;

  RealVect topPoint = RealVect::Zero;
  topPoint += channelNormal;
  topPoint *= channelRadius;

  PlaneIF bot(channelNormal,botPoint,true);
  PlaneIF top(channelNormal,topPoint,false);
  IntersectionIF channel(top,bot);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(channel,0,vectDx);
  //this generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int bigboxsize=2048;
  ebisPtr->define(a_domain, origin, a_dx, workshop, bigboxsize);

  Vector<Box> vbox(1, a_domain);
  Vector<int> proc(1, 0);
  DisjointBoxLayout dbl(vbox, proc);
  EBISLayout ebisl;
  ebisPtr->fillEBISLayout(ebisl, dbl, a_domain, 0);
}
/***************/
/***************/
void
compareError(const EBCellFAB& a_errorFine,
             const EBISBox& a_ebisBoxFine,
             const EBCellFAB& a_errorCoar,
             const EBISBox& a_ebisBoxCoar)
{
  Box gridFine = a_errorFine.getRegion();
  Box gridCoar = a_errorCoar.getRegion();
  pout() << "==============================================" << endl;
  for (int comp = 0; comp < a_errorFine.nComp(); comp++)
    {
      pout() << "Comparing error in variable  " << comp << endl;
      pout() << "==============================================" << endl;
      for (int itype = 0; itype < 3; itype++)
        {
          EBNormType::NormMode normtype;
          if (itype == 0)
            {
              normtype = EBNormType::OverBoth;
              pout() << endl << "Using all uncovered cells." << endl;
            }
          else if (itype == 1)
            {
              normtype = EBNormType::OverOnlyRegular;
              pout() << endl << "Using only regular cells." << endl;
            }
          else
            {
              normtype = EBNormType::OverOnlyIrregular;
              pout() << endl << "Using only irregular cells." << endl;
            }
          for (int inorm = 0; inorm < 3; inorm++)
            {

              if (inorm == 0)
                {
                  pout() << endl << "Using max norm." << endl;
                }
              else
                {
                  pout() << endl << "Using L-" << inorm << "norm." << endl;
                }
              Real coarnorm = EBArith::norm(a_errorCoar,
                                            gridCoar, a_ebisBoxCoar,
                                            comp, inorm, normtype);
              Real finenorm = EBArith::norm(a_errorFine,
                                            gridFine, a_ebisBoxFine,
                                            comp, inorm, normtype);

              pout() << "Coarse Error Norm = " << coarnorm << endl;
              pout() << "Fine   Error Norm = " << finenorm << endl;
              if ((Abs(finenorm) > 1.0e-8) && (Abs(coarnorm) > 1.0e-8))
                {
                  Real order = log(coarnorm/finenorm)/log(2.0);
                  pout() << "Order of scheme = " << order << endl;
                }
            }
        }
      pout() << "==============================================" << endl;
    }
}
/************/
Real
exactSoln(const RealVect& a_xval)
{
  Real retval = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real x = a_xval[idir];
      retval+= x*x*x;
    }
  return retval;
}
/************/
Real
exactDeriv(const RealVect& a_xval,
           const int&      a_idir)
{
  Real x = a_xval[a_idir];
  Real retval  = 3.*x*x;
  return retval;
}
/************/
/************/
void
setToExactSolut(EBCellFAB&     a_exactSoln,
                const EBISBox& a_ebisBox,
                const Real&    a_dx,
                const Box&     a_region)
{
  a_exactSoln.setVal(0.);
  IntVectSet ivsregion(a_region);
  for (VoFIterator vofit(ivsregion, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect xval;
      IntVect iv = vof.gridIndex();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          xval[idir] = (Real(iv[idir]) + 0.5)*a_dx;
        }
      a_exactSoln(vof, 0) = exactSoln(xval);
    }
}
/************/
/************/
void
setToExactDeriv(EBCellFAB&     a_exactDeriv,
                const EBISBox& a_ebisBox,
                const Real&    a_dx,
                const Box&     a_region,
                const int&     a_dir)
{
  a_exactDeriv.setVal(0.);
  IntVectSet ivsregion(a_region);
  for (VoFIterator vofit(ivsregion, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect xval;
      IntVect iv = vof.gridIndex();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          xval[idir] = (Real(iv[idir]) + 0.5)*a_dx;
        }
      a_exactDeriv(vof, 0) = exactDeriv(xval, a_dir);
    }
}
/************/
/**************************/
void
slopeTest()
{
  int nCons = 1;
  Box domainBoxFine, domainBoxCoar;
  Real dxFine, dxCoar;
  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }

  makeGeometry(domainBoxFine, dxFine);
  domainBoxCoar = coarsen(domainBoxFine, 2);
  ProblemDomain domainFine(domainBoxFine.smallEnd(),
                           domainBoxFine.bigEnd(),
                           is_periodic);
  ProblemDomain domainCoar(domainBoxCoar.smallEnd(),
                           domainBoxCoar.bigEnd(),
                           is_periodic);
  //debug
  dxFine = 1.0;

  dxCoar = 2.*dxFine;
  Vector<Box> boxFine(1, domainBoxFine);
  Vector<Box> boxCoar(1, domainBoxCoar);
  Vector<int> proc(1, 0);
  DisjointBoxLayout dblFine(boxFine, proc);
  DisjointBoxLayout dblCoar;
  coarsen(dblCoar, dblFine, 2);

  //make ebislayouts
  //no need for ghost since we are using the domain
  //for the box

  EBISLayout ebislFine, ebislCoar;
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  ebisPtr->fillEBISLayout(ebislFine, dblFine, domainBoxFine, 0);
  ebisPtr->fillEBISLayout(ebislCoar, dblCoar, domainBoxCoar, 0);
  for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
    {
      for (DataIterator dit = dblFine.dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBoxFine = ebislFine[dit()];
          const EBISBox& ebisBoxCoar = ebislCoar[dit()];

          EBCellFAB  slopErroFine(ebisBoxFine, domainBoxFine, nCons);
          EBCellFAB  slopErroCoar(ebisBoxCoar, domainBoxCoar, nCons);
          EBCellFAB  slopNLimFine(ebisBoxFine, domainBoxFine, nCons);
          EBCellFAB  slopNLimCoar(ebisBoxCoar, domainBoxCoar, nCons);
          EBCellFAB  slopCalcFine(ebisBoxFine, domainBoxFine, nCons);
          EBCellFAB  slopCalcCoar(ebisBoxCoar, domainBoxCoar, nCons);
          EBCellFAB  slopExacFine(ebisBoxFine, domainBoxFine, nCons);
          EBCellFAB  slopExacCoar(ebisBoxCoar, domainBoxCoar, nCons);
          EBCellFAB  solnExacFine(ebisBoxFine, domainBoxFine, nCons);
          EBCellFAB  solnExacCoar(ebisBoxCoar, domainBoxCoar, nCons);

          slopErroFine.setVal(0.);
          slopErroCoar.setVal(0.);
          slopNLimFine.setVal(0.);
          slopNLimCoar.setVal(0.);
          slopCalcFine.setVal(0.);
          slopCalcCoar.setVal(0.);
          slopExacFine.setVal(0.);
          slopExacCoar.setVal(0.);
          solnExacFine.setVal(0.);
          solnExacCoar.setVal(0.);

          setToExactSolut(solnExacFine, ebisBoxFine,  dxFine, domainBoxFine);
          setToExactSolut(solnExacCoar, ebisBoxCoar,  dxCoar, domainBoxCoar);
          setToExactDeriv(slopExacFine, ebisBoxFine,  dxFine, domainBoxFine, derivDir);
          setToExactDeriv(slopExacCoar, ebisBoxCoar,  dxCoar, domainBoxCoar, derivDir);

          IntVectSet cfivs;  Real time = 0; Real dt= 1.;   //dummy values
          EBPatchAdvect ebpatchAdvectFine, ebpatchAdvectCoar;
          ebpatchAdvectFine.define(domainFine, dxFine*RealVect::Unit, true);
          ebpatchAdvectCoar.define(domainCoar, dxCoar*RealVect::Unit, true);
          ebpatchAdvectFine.setValidBox(domainBoxFine, ebisBoxFine, cfivs, time, dt);
          ebpatchAdvectCoar.setValidBox(domainBoxCoar, ebisBoxCoar, cfivs, time, dt);

          //no flattening  or fourth order here
          ebpatchAdvectFine.setSlopeParameters(false, false, true);
          ebpatchAdvectCoar.setSlopeParameters(false, false, true);
          EBCellFAB flat;
          ebpatchAdvectFine.slope(slopCalcFine, slopNLimFine, solnExacFine, flat,  derivDir, domainBoxFine, true);
          ebpatchAdvectCoar.slope(slopCalcCoar, slopNLimCoar, solnExacCoar, flat,  derivDir, domainBoxCoar, true);

          //deriv = slope/dx
          slopNLimFine /= dxFine;
          slopNLimCoar /= dxCoar;

          slopErroFine.copy(slopNLimFine);
          slopErroCoar.copy(slopNLimCoar);

          slopErroFine -= slopExacFine;
          slopErroCoar -= slopExacCoar;
          pout() << "without limiting..." << endl;
          //do the whole convergence test thing.
          compareError(slopErroFine, ebisBoxFine,
                       slopErroCoar, ebisBoxCoar);


          //deriv = slope/dx
          slopCalcFine /= dxFine;
          slopCalcCoar /= dxCoar;

          slopErroFine.copy(slopCalcFine);
          slopErroCoar.copy(slopCalcCoar);

          slopErroFine -= slopExacFine;
          slopErroCoar -= slopExacCoar;

          pout() << "with limiting..." << endl;
          //do the whole convergence test thing.
          compareError(slopErroFine, ebisBoxFine,
                       slopErroCoar, ebisBoxCoar);
        }
    }

}
/************/
/************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  {
    // setChomboMPIErrorHandler();
#endif
    const char* in_file = "diverge.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);

    slopeTest();

#ifdef CH_MPI
  }
  MPI_Finalize();
#endif
}
