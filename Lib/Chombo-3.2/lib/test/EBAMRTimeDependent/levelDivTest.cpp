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
#include "PolyGeom.H"
#include "SphereIF.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "EBFluxFactory.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "EBPatchGodunovF_F.H"

#include "DebugDump.H"
#include "EBDebugDump.H"

#include "UsingNamespace.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

/***************/
void
sphereGeometry(Box& a_coarsestDomain,
               Real& a_dx)
{
  ParmParse ppgodunov;
  //parse input file
  int max_level = 0;
  ppgodunov.get("max_level",max_level);

  int num_read_levels = Max(max_level,1);
  std::vector<int> refRatios; // (num_read_levels,1);
  // note this requires a refRatio to be defined for the
  // finest level (even though it will never be used)

  ppgodunov.getarr("ref_ratio",refRatios,0,num_read_levels+1);
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

  a_coarsestDomain = Box(lo, hi);
  Box finestDomain = a_coarsestDomain;
  for (int ilev = 0; ilev < max_level; ilev++)
    {
      finestDomain.refine(refRatios[ilev]);
    }

  Real prob_hi;
  int numOpen = n_cell[0];
  pp.get("domain_length",prob_hi);
  a_dx = prob_hi/numOpen;
  Real fineDx = a_dx;
  int ebMaxCoarsen = 2;
  for (int ilev = 0; ilev < max_level; ilev++)
    {
      fineDx /= refRatios[ilev];
      ebMaxCoarsen += refRatios[ilev]/2;
    }
  int ebMaxSize;
  ppgodunov.get("max_grid_size", ebMaxSize);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  pout() << "sphere geometry" << endl;
  vector<Real> sphere_center(SpaceDim);
  pp.getarr("sphere_center",sphere_center, 0, SpaceDim);
  RealVect sphereCenter;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      sphereCenter[idir] = sphere_center[idir];
    }
  Real sphereRadius;
  pp.get("sphere_radius", sphereRadius);
  int inside_fluid;
  pp.get("sphere_fluid_inside", inside_fluid);
  bool insideFluid = (inside_fluid==1);
  SphereIF sphereIF(sphereRadius, sphereCenter, insideFluid);
  RealVect vectFineDx(D_DECL(fineDx,fineDx,fineDx));
  GeometryShop workshop(sphereIF,0,vectFineDx);
  //this generates the new EBIS
  ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize, ebMaxCoarsen);
}
/************/
Real
exactDivergence(const RealVect& a_xval)
{
  Real retval;
  Real x = a_xval[0];
  retval = 3*x*x;

  return retval;
}
/************/
Real
exactFlux(const RealVect& a_xval, int a_faceDir)
{
  Real retval = 0.0;
  if (a_faceDir == 0)
    {
      Real x = a_xval[0];
      retval = x*x*x;
    }

  return retval;
}
/************/
void
setToExactDivF(EBCellFAB&     a_exactDivF,
               const EBISBox& a_ebisBox,
               const Box&     a_region,
               const Real&    a_dx)
{
  a_exactDivF.setVal(0.);
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
      Real solnrv = exactDivergence(xval);
      Real kappa = a_ebisBox.volFrac(vof);
      a_exactDivF(vof,0) = kappa*solnrv;
    }
}
/************/
void
setToExactFlux(EBFluxFAB&             a_flux,
               const EBISBox&         a_ebisBox,
               const Box&             a_region,
               const Real&            a_dx)
{
  IntVectSet ivsregion(a_region);
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      for (FaceIterator faceit(ivsregion, a_ebisBox.getEBGraph(), faceDir, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
        {
          RealVect xval;
          IntVect iv = faceit().gridIndex(Side::Hi);
          RealVect centroid = a_ebisBox.centroid(faceit());
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (idir == faceDir)
                {
                  xval[idir] = (Real(iv[idir]))*a_dx;
                }
              else
                {
                  xval[idir] = (Real(iv[idir]) + 0.5 + centroid[idir])*a_dx;
                }
            }
          Real fluxDir = exactFlux(xval, faceDir);
          a_flux[faceDir](faceit(), 0) = fluxDir;
        }
    } //end loop over face directions

}
/************/
void
kappaDivergence(EBCellFAB&             a_divF,
                const EBFluxFAB&       a_flux,
                const EBISBox&         a_ebisBox,
                const Box&             a_box,
                const Real&            a_dx)
{
  //set the divergence initially to zero
  //then loop through directions and increment the divergence
  //with each directions flux difference.
  a_divF.setVal(0.0);
  BaseFab<Real>&       regDivF = a_divF.getSingleValuedFAB();
  regDivF.setVal(0.);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      //update for the regular vofs in the nonconservative
      //case  works for all single valued vofs.
      /* do the regular vofs */
      /**/
      const EBFaceFAB& fluxDir = a_flux[idir];
      const BaseFab<Real>& regFluxDir = fluxDir.getSingleValuedFAB();
      int ncons = 1;
      FORT_DIVERGEF( CHF_BOX(a_box),
                     CHF_FRA(regDivF),
                     CHF_CONST_FRA(regFluxDir),
                     CHF_CONST_INT(idir),
                     CHF_CONST_INT(ncons),
                     CHF_CONST_REAL(a_dx));
      /**/
    }
  //update the irregular vofs using conservative diff
  IntVectSet ivsIrreg = a_ebisBox.getIrregIVS(a_box);
  for (VoFIterator vofit(ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      //divergence was set in regular update.  we reset it
      // to zero and recalc.
      Real update = 0.;
      for ( int idir = 0; idir < SpaceDim; idir++)
        {
          const EBFaceFAB& fluxDir = a_flux[idir];
          for (SideIterator sit; sit.ok(); ++sit)
            {
              int isign = sign(sit());
              Vector<FaceIndex> faces =
                a_ebisBox.getFaces(vof, idir, sit());
              for (int iface = 0; iface < faces.size(); iface++)
                {
                  const FaceIndex& face = faces[iface];
                  Real areaFrac = a_ebisBox.areaFrac(face);
                  Real faceFlux =fluxDir(face, 0);
                  update += isign*areaFrac*faceFlux;

                }
            }
        }
      //add EB boundary condtions in divergence
      const IntVect& iv = vof.gridIndex();
      Real bndryArea = a_ebisBox.bndryArea(vof);
      RealVect bndryCent = a_ebisBox.bndryCentroid(vof);
      RealVect normal = a_ebisBox.normal(vof);
      RealVect bndryLoc;
      RealVect exactF;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          bndryLoc[idir] = a_dx*(iv[idir] + 0.5 + bndryCent[idir]);
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          exactF[idir] = exactFlux(bndryLoc, idir);
        }
      Real bndryFlux = PolyGeom::dot(exactF, normal);

      update -= bndryFlux*bndryArea;
      update /= a_dx; //note NOT divided by volfrac

      a_divF(vof, 0) = update;
    }
}
/***************/
void
kappaDivergenceLD(LevelData<EBCellFAB>&       a_divF,
                  const LevelData<EBFluxFAB>& a_flux,
                  const EBISLayout&           a_ebisl,
                  const DisjointBoxLayout&    a_dbl,
                  const Real&                 a_dx)
{
  for (DataIterator dit= a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      kappaDivergence(a_divF[dit()],
                      a_flux[dit()],
                      a_ebisl[dit()],
                      a_dbl.get(dit()),
                      a_dx);
    }

}

/**************************/
void
setToExactFluxLD(LevelData<EBFluxFAB>&       a_flux,
                 const EBISLayout&           a_ebisl,
                 const DisjointBoxLayout&    a_dbl,
                 const Real&                 a_dx)
{
  for (DataIterator dit= a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      setToExactFlux(a_flux[dit()],
                     a_ebisl[dit()],
                     a_dbl.get(dit()),
                     a_dx);
    }
}
/**************************/
void
setToExactDivFLD(LevelData<EBCellFAB>&       a_soln,
                 const EBISLayout&           a_ebisl,
                 const DisjointBoxLayout&    a_dbl,
                 const Real&                 a_dx)
{
  for (DataIterator dit= a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      setToExactDivF(a_soln[dit()],
                     a_ebisl[dit()],
                     a_dbl.get(dit()),
                     a_dx);
    }
}
/**************************/
void
compareError(const LevelData<EBCellFAB>& a_errorFine,
             const EBISLayout& a_ebislFine,
             const DisjointBoxLayout& a_gridsFine,
             const Box& a_domainFine,
             const LevelData<EBCellFAB>& a_errorCoar,
             const EBISLayout& a_ebislCoar,
             const DisjointBoxLayout& a_gridsCoar,
             const Box& a_domainCoar)
{
 CH_assert(a_errorFine.nComp() == 1);
 CH_assert(a_errorCoar.nComp() == 1);

  pout() << "==============================================" << endl;

  EBNormType::NormMode normtype = EBNormType::OverBoth;
  pout() << endl << "Using all uncovered cells." << endl  ;

  for (int inorm = 0; inorm <= 2; inorm++)
    {

      if (inorm == 0)
        {
          pout() << endl << "Using max norm." << endl;
        }
      else
        {
          pout() << endl << "Using L-" << inorm << "norm." << endl;
        }
      int comp = 0;
      Real coarnorm = EBArith::norm(a_errorCoar,
                                    a_gridsCoar, a_ebislCoar,
                                    comp, inorm, normtype);
      Real finenorm = EBArith::norm(a_errorFine,
                                    a_gridsFine, a_ebislFine,
                                    comp, inorm, normtype);

      pout() << "Coarse Error Norm = " << coarnorm << endl;
      pout() << "Fine   Error Norm = " << finenorm << endl;

      if ((Abs(finenorm) > 1.0e-12) && (Abs(coarnorm) > 1.0e-12))
        {
          Real order = log(Abs(coarnorm/finenorm))/log(2.0);
          pout() << "Order of scheme = " << order << endl;
        }
    }
  pout() << "==============================================" << endl ;;
}
/**************************/
void
getError(LevelData<EBCellFAB>&       a_error,
         const EBISLayout&           a_ebisl,
         const DisjointBoxLayout&    a_dbl,
         const Real&                 a_dx)
{
  EBCellFactory ebcellfact(a_ebisl);
  EBFluxFactory ebfluxfact(a_ebisl);
  a_error.define(a_dbl, 1, IntVect::Zero,   ebcellfact);
  LevelData<EBCellFAB> divFCalc(a_dbl, 1, IntVect::Zero,   ebcellfact);
  LevelData<EBCellFAB> divFExac(a_dbl, 1, IntVect::Zero,   ebcellfact);
  LevelData<EBFluxFAB> faceFlux(a_dbl, 1, IntVect::Zero,   ebfluxfact);

  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      a_error[dit()].setVal(0.);
      divFCalc[dit()].setVal(0.);
      divFExac[dit()].setVal(0.);
    }

  setToExactDivFLD(divFExac,  a_ebisl, a_dbl, a_dx);
  setToExactFluxLD(faceFlux,  a_ebisl, a_dbl, a_dx);

  Interval interv(0, 0);
  faceFlux.exchange(interv);

  kappaDivergenceLD(divFCalc, faceFlux, a_ebisl, a_dbl, a_dx);

  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& errorFAB = a_error[dit()];
      EBCellFAB& exactFAB = divFExac[dit()];
      EBCellFAB& calcuFAB = divFCalc[dit()];

      errorFAB += calcuFAB;
      errorFAB -= exactFAB;
    }
}

/**************/
void
divergenceTest()
{
  int maxsize;
  ParmParse pp;
  pp.get("max_grid_size", maxsize);

  //make layouts == domain
  Box domainBoxFine, domainBoxCoar;
  Real dxFine, dxCoar;
  sphereGeometry(domainBoxFine, dxFine);
  domainBoxCoar = coarsen(domainBoxFine, 2);

  //debug
  dxFine = 1.0;
  dxCoar = 2.*dxFine;
  Vector<Box> boxFine;
  domainSplit(domainBoxFine, boxFine, maxsize);

  Vector<int> proc(boxFine.size(), 0);
  DisjointBoxLayout dblFine(boxFine, proc);
  DisjointBoxLayout dblCoar;
  coarsen(dblCoar, dblFine, 2);
  LevelData<EBCellFAB> errorFine, errorCoar;

  EBISLayout ebislFine, ebislCoar;
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  ebisPtr->fillEBISLayout(ebislFine, dblFine, domainBoxFine, 2);
  ebisPtr->fillEBISLayout(ebislCoar, dblCoar, domainBoxCoar, 2);

  getError(errorFine, ebislFine, dblFine, dxFine);
  getError(errorCoar, ebislCoar, dblCoar, dxCoar);

  for (DataIterator dit= dblFine.dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& errorFineFAB = errorFine[dit()];
      const EBCellFAB& errorCoarFAB = errorCoar[dit()];
      int i1 = errorFineFAB.getMultiCells().numPts();
      int i2 = errorCoarFAB.getMultiCells().numPts();
      pout() << i1 << i2 << endl;
    }

  compareError(errorFine, ebislFine, dblFine, domainBoxFine,
               errorCoar, ebislCoar, dblCoar, domainBoxCoar);

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
    // Check for an input file
    const char* inFile = "sphere.inputs";

    ParmParse pp(0,NULL,NULL,inFile);

    divergenceTest();

    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    ebisPtr->clear();
#ifdef CH_MPI
  }
  MPI_Finalize();
#endif
}
