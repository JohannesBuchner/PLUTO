#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif
#include <cmath>

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "LevelData.H"
#include "EBFluxFAB.H"
#include "EBFluxFactory.H"
#include "EBCoarseAverage.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "EBArith.H"
#include "AllRegularService.H"
#include "PlaneIF.H"
#include "EBAMRIO.H"
#include "AMRIO.H"

#include "UsingNamespace.H"

/******************/
/******************/
int
makeLayout(DisjointBoxLayout& dbl1,
           const Box& domain);
/***************/
/***************/
int makeGeometry(Box& domain,
                 Real& dx);
/***************/
/***************/
int makeEBISL(EBISLayout& a_ebisl,
              const DisjointBoxLayout& a_grids,
              const Box& a_domain,
              const int& nghost);

/***************/
/***************/
int getError(LevelData<EBFluxFAB>&    a_errorCoar,
             const EBISLayout&        a_ebislFine,
             const DisjointBoxLayout& a_gridsFine,
             const Box&               a_domainFine,
             const Real&              a_dxFine,
             const EBISLayout&        a_ebislCoar,
             const DisjointBoxLayout& a_gridsCoar,
             const Box&               a_domainCoar,
             const Real&              a_dxCoar,
             const int&               a_nref);
/***************/
/***************/
int compareError(const LevelData<EBFluxFAB>& a_fineError,
                 const EBISLayout& a_ebislFine,
                 const DisjointBoxLayout& a_gridsFine,
                 const Box& a_domainFine,
                 const LevelData<EBFluxFAB>& a_coarError,
                 const EBISLayout& a_ebislCoar,
                 const DisjointBoxLayout& a_gridsCoar,
                 const Box& a_domainCoar,
                 const int& a_nref);

/***************/
/***************/
Real exactFunc(const FaceIndex& a_face,
               const Real& a_dx);

/***************/
/***************/
Real maxNormOneComp(const LevelData<EBFluxFAB>& a_data,
                    const int&                  a_dir,
                    const int&                  a_comp);

/***************/
/***************/
void dumpmemoryatexit();
/***************/
/***************/
int
main(int argc, char** argv)

{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int eekflag = 0;
  //begin forever present scoping trick
  {
    const char* in_file = "ramp.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    //define the geometry object.
    //need the new and delete because of
    //strong construction
    //make the gometry.  this makes the first Chombo_EBIS
    Box domainFine, domainMedi, domainCoar;
    Real dxFine, dxMedi, dxCoar;
    //and defines it using a geometryservice
    eekflag =  makeGeometry(domainFine,  dxFine);
    CH_assert(eekflag == 0);
    int nref = 2;
    pp.get("ref_ratio",nref);
    domainMedi = coarsen(domainFine, nref);
    domainCoar = coarsen(domainMedi, nref);
    dxMedi = nref*dxFine;
    dxCoar = nref*dxMedi;
    //make grids
    DisjointBoxLayout gridsFine, gridsMedi, gridsCoar;
    eekflag = makeLayout(gridsFine, domainFine);
    CH_assert(eekflag == 0);

    coarsen(gridsMedi, gridsFine, nref);
    coarsen(gridsCoar, gridsMedi, nref);

    ///create ebislayout
    int nghost = 2;
    EBISLayout ebislFine, ebislMedi, ebislCoar;
    eekflag = makeEBISL(ebislFine, gridsFine, domainFine, nghost);
    if (eekflag != 0) return eekflag;
    eekflag = makeEBISL(ebislMedi, gridsMedi, domainMedi, nghost);
    if (eekflag != 0) return eekflag;
    eekflag = makeEBISL(ebislCoar, gridsCoar, domainCoar, nghost);
    if (eekflag != 0) return eekflag;

    int nvar = 1;
    EBFluxFactory ebfluxfactMedi(ebislMedi);
    EBFluxFactory ebfluxfactCoar(ebislCoar);
    LevelData<EBFluxFAB> errorMedi(gridsMedi, nvar,
                                   IntVect::Zero,ebfluxfactMedi);
    LevelData<EBFluxFAB> errorCoar(gridsCoar, nvar,
                                   IntVect::Zero,ebfluxfactCoar);

    eekflag = getError(errorMedi,
                       ebislFine, gridsFine, domainFine, dxFine,
                       ebislMedi, gridsMedi, domainMedi, dxMedi, nref);
    if (eekflag != 0) return eekflag;

    eekflag = getError(errorCoar,
                       ebislMedi, gridsMedi, domainMedi, dxMedi,
                       ebislCoar, gridsCoar, domainCoar, dxCoar, nref);
    if (eekflag != 0) return eekflag;

    eekflag =  compareError(errorMedi, ebislMedi, gridsMedi, domainMedi,
                            errorCoar, ebislCoar, gridsCoar, domainCoar, nref);
    if (eekflag != 0) return eekflag;

  }//end scoping trick
  CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return eekflag;
}
/***************/
/***************/
int getError(LevelData<EBFluxFAB>& a_errorCoar,
             const EBISLayout& a_ebislFine,
             const DisjointBoxLayout& a_gridsFine,
             const Box& a_domainFine,
             const Real& a_dxFine,
             const EBISLayout& a_ebislCoar,
             const DisjointBoxLayout& a_gridsCoar,
             const Box& a_domainCoar,
             const Real& a_dxCoar,
             const int&  a_nref)
{
  int eekflag = 0;
  int nvar = 1;

  EBFluxFactory ebfluxfactFine(a_ebislFine);
  EBFluxFactory ebfluxfactCoar(a_ebislCoar);
  LevelData<EBFluxFAB> phiFine(a_gridsFine, nvar,
                               IntVect::Zero,ebfluxfactFine);
  LevelData<EBFluxFAB> phiCoar(a_gridsCoar, nvar,
                               IntVect::Zero,ebfluxfactCoar);

  //fill phi fine and phiCExact
  //put phiCexact into a_errorCoar
  for (DataIterator dit = a_gridsFine.dataIterator();
      dit.ok(); ++dit)
    {
      EBFluxFAB& phiFineFluxFAB = phiFine[dit()];
      EBFluxFAB& phiCoarFluxFAB = a_errorCoar[dit()];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          EBFaceFAB& phiFineFAB = phiFineFluxFAB[idir];
          EBFaceFAB& phiCoarFAB = phiCoarFluxFAB[idir];
          {
            IntVectSet ivsBox(a_gridsFine.get(dit()));
            for (FaceIterator faceit(ivsBox, a_ebislFine[dit()].getEBGraph(),idir,
                                    FaceStop::SurroundingWithBoundary);
                faceit.ok(); ++faceit)
              {
                const FaceIndex& face = faceit();
                Real rightAns = exactFunc(face, a_dxFine);
                phiFineFAB(face, 0) = rightAns;
              }
          }

          {
            IntVectSet ivsBox(a_gridsCoar.get(dit()));
            for (FaceIterator faceit(ivsBox, a_ebislCoar[dit()].getEBGraph(),idir,
                                    FaceStop::SurroundingWithBoundary);
                faceit.ok(); ++faceit)
              {
                const FaceIndex& face = faceit();
                Real rightAns = exactFunc(face, a_dxCoar);
                phiCoarFAB(face, 0) = rightAns;
              }
          }
        }
    }

  //average phiFine onto phiC;
  EBCoarseAverage ebAveOp(a_gridsFine, a_gridsCoar,
                          a_ebislFine, a_ebislCoar,
                          a_domainCoar, a_nref, nvar, Chombo_EBIS::instance());

  Interval zeroiv(0,0);
  ebAveOp.average(phiCoar, phiFine, zeroiv);
  //error = phiC - phiCExact
  for (DataIterator dit = a_gridsCoar.dataIterator();
      dit.ok(); ++dit)
    {
      EBFluxFAB& errorFAB = a_errorCoar[dit()];
      EBFluxFAB& phiCoarFAB = phiCoar[dit()];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          errorFAB[idir] -= phiCoarFAB[idir];
        }
    }

  return eekflag;
}
/***************/
/***************/
int compareError(const LevelData<EBFluxFAB>& a_errorFine,
                 const EBISLayout& a_ebislFine,
                 const DisjointBoxLayout& a_gridsFine,
                 const Box& a_domainFine,
                 const LevelData<EBFluxFAB>& a_errorCoar,
                 const EBISLayout& a_ebislCoar,
                 const DisjointBoxLayout& a_gridsCoar,
                 const Box& a_domainCoar,
                 const int& a_nref)
{
  EBFluxFactory factFine(a_ebislFine);
  EBFluxFactory factCoar(a_ebislCoar);
  int eekflag = 0;
  int comp = 0;
  for (int inorm = 0; inorm < 1; inorm++)
    {

      if (inorm == 0)
        {
          pout() << endl << " Using max norm." << endl;
        }
      else
        {
          pout() << endl << " Using L-" << inorm << "norm." << endl;
          MayDay::Error("max norm for LevelData<EBFluxFAB>'s not implemented yet");
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {

          Real coarnorm = maxNormOneComp(a_errorCoar,
                                         idir,
                                         comp);

          Real finenorm = maxNormOneComp(a_errorFine,
                                         idir,
                                         comp);

          pout() << "Direction = " << idir << " Coarse Error Norm = " << coarnorm << endl;
          pout() << "Direction = " << idir << "   Fine Error Norm = " << finenorm << endl;
          if (Abs(finenorm) > 1.0e-8)
            {
              Real order = log(coarnorm/finenorm)/log(Real(a_nref));
              pout() << "   Order of scheme = " << order << endl;
            }
        }
    }

  return eekflag;
}

Real maxNormOneComp(const LevelData<EBFluxFAB>& a_data,
                    const int&                  a_dir,
                    const int&                  a_comp)
{
  int nComp = a_data.nComp();
  CH_assert(a_comp<nComp);

  Real max = -1.0e99;

  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;
  for (DataIterator dit = a_data.dataIterator();dit.ok();++dit)
    {
      const DisjointBoxLayout &dbl = a_data.disjointBoxLayout();
      IntVectSet ivsBox(dbl.get(dit()));
      const EBFluxFAB& dataFluxFAB = a_data[dit()];
      const EBFaceFAB& dataFaceFAB = dataFluxFAB[a_dir];
      const EBISBox&       ebisBox = dataFaceFAB.getEBISBox();
      for (FaceIterator faceit(ivsBox, ebisBox.getEBGraph(), a_dir, stopCrit);faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          const Real& val = dataFaceFAB(face, a_comp);
          if (isnan(val) || isinf(val))
            {
#ifdef CH_USE_HDF5
              char* fname = tempnam(NULL,NULL);
              const FArrayBox& reg = dataFaceFAB.getFArrayBox();
              writeFABname(&reg, fname);
              pout() << "   bad data written to an hdf5 file: " << fname << endl;
#endif
              pout() << "dir = " << a_dir << "comp = " << a_comp << "val = " << val << endl;
              MayDay::Error("face val isnan or isinf");
            }
          Real dataFace = Abs(val);
          if (dataFace>max)
            {
              max = dataFace;
            }
        }
    }

  //find the max over all processors
  //    Real pMax = parallelMax(max);
  Real pMax = max;
  return pMax;
}
/***************/
/***************/
Real exactFunc(const FaceIndex& a_face, const Real& a_dx)
{
  Real retval;
  Real probHi;
  ParmParse pp;
  pp.get("prob_hi",probHi);
  RealVect xloc;
  IntVect loIV = a_face.gridIndex(Side::Lo);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xloc[idir] = a_dx*(Real(loIV[idir]) + 0.5);
    }
  int dir = a_face.direction();
  xloc[dir] += 0.5*a_dx;
  Real x = xloc[0]/probHi;
  Real y = xloc[1]/probHi;
  retval = 1.0 + x + y + x*x + y*y;
  //  retval = x*x;
  //  retval = y*y;
  //  retval = 1.0 + x + y;
  //  retval = y;
  //  retval = x;
  //  retval = 1.0;
  //  retval = xloc[0]*xloc[0];
  return retval;
}
/***************/
/***************/
int makeEBISL(EBISLayout& a_ebisl,
              const DisjointBoxLayout& a_grids,
              const Box& a_domain,
              const int& a_nghost)
{
  const CH_XD::EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  CH_assert(ebisPtr->isDefined());
  ebisPtr->fillEBISLayout(a_ebisl, a_grids, a_domain, a_nghost);
  return 0;
}
/***************/
/***************/
int
makeLayout(DisjointBoxLayout& a_dbl,
           const Box& a_domain)
{
  ParmParse pp;
  int eekflag= 0;
  int ipieces;
  ipieces = Max(ipieces, 1);
  int maxsize;
  pp.get("maxboxsize",maxsize);
  Vector<Box> vbox(1, a_domain);
  domainSplit(a_domain, vbox,  maxsize);
  if (eekflag != 0)
    {
      pout() << "problem in domainsplit" << endl;
      return eekflag;
    }
  Vector<int>  procAssign;
  eekflag = LoadBalance(procAssign,vbox);
  if (eekflag != 0)
    {
      pout() << "problem in loadbalance" << endl;
      return eekflag;
    }
  a_dbl.define(vbox, procAssign);
  return eekflag;
}
/**********/
int makeGeometry(Box& a_domain,
                 Real& a_dx)
{
  int eekflag =  0;
  //parse input file
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
          return(-1);
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_domain.setSmall(lo);
  a_domain.setBig(hi);

  Vector<Real> prob_lo(SpaceDim, 1.0);
  Real prob_hi;
  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);
  a_dx = (prob_hi-prob_lo[0])/n_cell[0];
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      origin[idir] = prob_lo[idir];
    }
  int whichgeom;
  pp.get("which_geom",whichgeom);
  if (whichgeom == 0)
    {
      //allregular
      pout() << "all regular geometry" << endl;
      AllRegularService regserv;
      CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      ebisPtr->define(a_domain, origin, a_dx, regserv);
    }
  else if (whichgeom == 1)
    {
      pout() << "ramp geometry" << endl;
      int upDir;
      int indepVar;
      Real startPt;
      Real slope;
      pp.get("up_dir",upDir);
      pp.get("indep_var",indepVar);
      pp.get("start_pt", startPt);
      pp.get("ramp_slope", slope);

      RealVect normal = RealVect::Zero;
      normal[upDir] = 1.0;
      normal[indepVar] = -slope;

      RealVect point = RealVect::Zero;
      point[upDir] = -slope*startPt;

      bool normalInside = true;

      PlaneIF ramp(normal,point,normalInside);

      RealVect vectDx = RealVect::Unit;
      vectDx *= a_dx;

      GeometryShop workshop(ramp,0,vectDx);
      //this generates the new EBIS
      CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      ebisPtr->define(a_domain, origin, a_dx, workshop);
    }
  else
    {
      //bogus which_geom
      pout() << " bogus which_geom input = " << whichgeom;
      eekflag = 33;
    }

  return eekflag;
}
/***************/
/***************/
