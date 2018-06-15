#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
#include <iostream>
#include <fstream>
#include <stdio.h>

#include "AMRINSUtils.H"

#include "Box.H"
#include "parstream.H"
#include "ParmParse.H"
#include "Vector.H"
#include "BRMeshRefine.H"

#include "EBAMRIO.H"
#include "EBIndexSpace.H"
#include "SlabService.H"
#include "MultiSlabService.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "EBArith.H"
#include "EBEllipticLoadBalance.H"

#include "SphereIF.H"
#include "MultiSphereIF.H"
#include "SphereArrayIF.H"
#include "PlaneIF.H"
#include "IntersectionIF.H"
#include "UnionIF.H"
#include "ComplementIF.H"
#include "AllRegularService.H"
#include "PlaneIF.H"
#include "TiltedCylinderIF.H"
#include "EllipsoidIF.H"
#include "TransformIF.H"
#include "PolynomialIF.H"
#include "TylerChannelIF.H"
#include "DEMIF.H"
#include "DataFileIF.H"
#include "ArteryIF.H"
#include "HelicoilIF.H"

#include "UsingNamespace.H"

using std::ifstream;
using std::ios;

BaseIF* makeChamber(const Real& radius,
                    const Real& thick,
                    const Real& offset,
                    const Real& height);

BaseIF* makePlate(const Real& height,
                  const Real& thick,
                  const Real& radius,
                  const int&  doHoles,
                  const Real& holeRadius,
                  const Real& holeSpace);

BaseIF* makeVanes(const int&      num,
                  const Real&     thick,
                  const RealVect& normal,
                  const Real&     innerRadius,
                  const Real&     outerRadius,
                  const Real&     offset,
                  const Real&     height);

BaseIF* makeVane(const Real&     thick,
                 const RealVect& normal,
                 const Real&     innerRadius,
                 const Real&     outerRadius,
                 const Real&     offset,
                 const Real&     height,
                 const Real&     angle);

/***************/
void outputError(const Vector< LevelData<EBCellFAB>* >&   a_errorFine,
                 const Vector< LevelData<EBCellFAB>* >&   a_errorCoar,
                 const Vector< DisjointBoxLayout >&       a_gridsFine,
                 const Vector< DisjointBoxLayout >&       a_gridsCoar,
                 const ProblemDomain&                     a_level0DomainFine,
                 const ProblemDomain&                     a_level0DomainCoar,
                 const string&                            a_fileFine,
                 const string&                            a_fileCoar,
                 Vector<string>&                          a_names,
                 const AMRParameters&                     a_params)

{
#ifdef CH_USE_HDF5
  int nvar = a_names.size();
  bool replaceCovered = true;
  Vector<Real> coveredValues(nvar, 0.0);
  //values that don't matter in output file
  Real dxFine = 1.0;
  Real dxCoar = 1.0;
  Real time   = 1.0;
  Real dt     = 1.0;

  ParmParse pp;

  Vector<int> refRatio = a_params.m_refRatio;
  int numlevels = a_params.m_maxLevel + 1;
  Box domainFine = a_level0DomainFine.domainBox();
  Box domainCoar = a_level0DomainCoar.domainBox();

  writeEBHDF5(a_fileFine, a_gridsFine, a_errorFine,  a_names,
              domainFine, dxFine, dt, time, refRatio, numlevels,
              replaceCovered, coveredValues);

  writeEBHDF5(a_fileCoar, a_gridsCoar, a_errorCoar,  a_names,
              domainCoar, dxCoar, dt, time, refRatio, numlevels,
              replaceCovered, coveredValues);
#endif
}
/****/
void
getFixedLayouts(Vector<DisjointBoxLayout>&      a_gridsFine,
                Vector<DisjointBoxLayout>&      a_gridsMedi,
                Vector<DisjointBoxLayout>&      a_gridsCoar,
                Vector<EBISLayout       >&      a_ebislFine,
                Vector<EBISLayout       >&      a_ebislMedi,
                Vector<EBISLayout       >&      a_ebislCoar,
                const ProblemDomain&            a_domainFine,
                const ProblemDomain&            a_domainMedi,
                const ProblemDomain&            a_domainCoar,
                const AMRParameters&            a_params,
                int a_whichGeomForced,
                const EBIndexSpace* a_ebisPtr )
{


  Vector< Vector<Box> >  boxesFine;
  Vector< Vector<Box> >  boxesMedi;
  Vector< Vector<Box> >  boxesCoar;

  getFixedGrids(boxesFine,
                boxesMedi,
                boxesCoar,
                a_params,
                a_domainCoar,
                a_whichGeomForced);

  int nlev = boxesCoar.size();
  a_gridsFine.resize(nlev);
  a_gridsMedi.resize(nlev);
  a_gridsCoar.resize(nlev);
  a_ebislFine.resize(nlev);
  a_ebislMedi.resize(nlev);
  a_ebislCoar.resize(nlev);
  int numEBG =   4;

  ProblemDomain domLevFine = a_domainFine;
  ProblemDomain domLevMedi = a_domainMedi;
  ProblemDomain domLevCoar = a_domainCoar;

  for (int ilev = 0; ilev < nlev; ilev++)
    {
      Vector<int> procsCoar, procsFine, procsMedi;
      EBEllipticLoadBalance(procsCoar, boxesCoar[ilev],domLevCoar);
      EBEllipticLoadBalance(procsMedi, boxesMedi[ilev],domLevFine);
      EBEllipticLoadBalance(procsFine, boxesFine[ilev],domLevFine);
      a_gridsFine[ilev] = DisjointBoxLayout(boxesFine[ilev], procsFine);
      a_gridsMedi[ilev] = DisjointBoxLayout(boxesMedi[ilev], procsMedi);
      a_gridsCoar[ilev] = DisjointBoxLayout(boxesCoar[ilev], procsCoar);

      a_ebisPtr->fillEBISLayout(a_ebislFine[ilev], a_gridsFine[ilev], domLevFine.domainBox(), numEBG);
      a_ebisPtr->fillEBISLayout(a_ebislMedi[ilev], a_gridsMedi[ilev], domLevMedi.domainBox(), numEBG);
      a_ebisPtr->fillEBISLayout(a_ebislCoar[ilev], a_gridsCoar[ilev], domLevCoar.domainBox(), numEBG);

      int refRat = a_params.m_refRatio[ilev];
      domLevFine.refine(refRat);
      domLevMedi.refine(refRat);
      domLevCoar.refine(refRat);
    }
}

//kallemov 06/22/10
//adding overladed function getFixedGrids for single grid
void
getFixedGrids(Vector< Vector<Box> >&  a_coarBoxes,
              const AMRParameters&    a_params,
              const ProblemDomain&    a_coarsestDomain,
              int a_whichGeomForced,
              const EBIndexSpace* a_ebisPtr )
{
  //first generate coarsest boxes
  int whichgeom;
  if (a_whichGeomForced >= 0)
    {
      whichgeom = a_whichGeomForced;
    }
  else
    {
      ParmParse pp;
      pp.get("which_geom",whichgeom);
    }
  int numLevels = a_params.m_maxLevel + 1;
  if (whichgeom == 0)
    {
      pout() << "all regular geometry" << endl;
      pout() << "ignoring grid parameters and making simple grids" << endl;
      Vector<Box> coarBoxes(1, a_coarsestDomain.domainBox());
      a_coarBoxes.resize(numLevels);
      a_coarBoxes[0] = coarBoxes;

      Box coarBox = a_coarsestDomain.domainBox();
      for (int ilev = 1; ilev < numLevels; ilev++)
        {
          int iboxShrink = coarBox.size(0);
          iboxShrink /= 4;
          if (iboxShrink < 2)
            {
              MayDay::Error("wacky DBL generation technique failed, try making base box bigger");
            }
          coarBox.grow(-iboxShrink);
          coarBox.refine(a_params.m_refRatio[ilev-1]);
          Vector<Box> refBoxes(1, coarBox);
          a_coarBoxes[ilev] = refBoxes;
        }
    }
  else
    {
      pout() << "making grids at coarsest level by tagging irregular cells" << endl;
      pout() << "all grid parameters refer to coarsest hierarchy (max box size, etc)" << endl;
      //split up coarsest domain by max box size and
      //make a dbl at coarsest level
      Vector<Box> boxesCoarsest;
      Vector<int> procsCoarsest;

      domainSplit(a_coarsestDomain.domainBox(), boxesCoarsest,
                  a_params.m_maxBoxSize, a_params.m_blockFactor);

      a_coarBoxes.resize(numLevels);
      a_coarBoxes[0] = boxesCoarsest;
      if (numLevels > 1)
        {
          EBEllipticLoadBalance(procsCoarsest, boxesCoarsest, a_coarsestDomain);
          DisjointBoxLayout dblCoarsest(boxesCoarsest, procsCoarsest);

          //make a ebislayout at coarsest level.
          EBISLayout ebislCoarsest;
          a_ebisPtr->fillEBISLayout(ebislCoarsest,dblCoarsest,
                                  a_coarsestDomain.domainBox(), 0);

          //loop through grids and get all irregular cells into an intvectset.
          IntVectSet tagsCoarsestLocal;
          for (DataIterator dit = dblCoarsest.dataIterator(); dit.ok(); ++dit)
            {
              const EBISBox& ebisBox = ebislCoarsest[dit()];
              const Box&       box = dblCoarsest.get(dit());
              tagsCoarsestLocal |= ebisBox.getIrregIVS(box);
            }

          //generate vector of grids.
          BRMeshRefine gridder(a_coarsestDomain,     a_params.m_refRatio,
                               a_params.m_fillRatio, a_params.m_blockFactor,
                               a_params.m_tagBuffer, a_params.m_maxBoxSize);

          Vector<Vector<Box> > newMeshes(numLevels);
          Vector<Vector<Box> > oldMeshes(numLevels);
          oldMeshes[0]= boxesCoarsest;
          ProblemDomain curDomain = a_coarsestDomain;
          for (int ilev = 1; ilev < numLevels; ilev++)
            {
              domainSplit(curDomain.domainBox(), oldMeshes[ilev],
                          a_params.m_maxBoxSize, a_params.m_blockFactor);
              curDomain.refine( a_params.m_refRatio[ilev-1]);
            }
          int baseLevel = 0;
          gridder.regrid(newMeshes, tagsCoarsestLocal, baseLevel, a_params.m_maxLevel, oldMeshes);
          //   gridder.regrid(newMeshes, tagsCoarsest, baseLevel, a_params.m_maxLevel, oldMeshes);
          for (int ilev = 0; ilev < numLevels; ilev++)
            {
              a_coarBoxes[ilev] = newMeshes[ilev];
            }
        }
    }

  //kallemov
  //we need only one grid
  /*
  //now generate finer grids as refinements
  pout() << "generating finer hierarchies as refinements of coarsest hierarchy" << endl;
  a_mediBoxes = a_coarBoxes;
  a_fineBoxes = a_coarBoxes;
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      for (int ibox = 0;  ibox < a_coarBoxes[ilev].size(); ibox++)
        {
          a_mediBoxes[ilev][ibox].refine(2);
          a_fineBoxes[ilev][ibox].refine(4);
        }
    }
  */

}
//end kallemov

void
getFixedGrids(Vector< Vector<Box> >&  a_fineBoxes,
              Vector< Vector<Box> >&  a_mediBoxes,
              Vector< Vector<Box> >&  a_coarBoxes,
              const AMRParameters&    a_params,
              const ProblemDomain&    a_coarsestDomain,
              int a_whichGeomForced,
              const EBIndexSpace* a_ebisPtr )
{
  //first generate coarsest boxes
  int whichgeom;
  if (a_whichGeomForced >= 0)
    {
      whichgeom = a_whichGeomForced;
    }
  else
    {
      ParmParse pp;
      pp.get("which_geom",whichgeom);
    }
  int numLevels = a_params.m_maxLevel + 1;
  if (whichgeom == 0)
    {
      pout() << "all regular geometry" << endl;
      pout() << "ignoring grid parameters and making simple grids" << endl;
      Vector<Box> coarBoxes(1, a_coarsestDomain.domainBox());
      a_coarBoxes.resize(numLevels);
      a_coarBoxes[0] = coarBoxes;

      Box coarBox = a_coarsestDomain.domainBox();
      for (int ilev = 1; ilev < numLevels; ilev++)
        {
          int iboxShrink = coarBox.size(0);
          iboxShrink /= 4;
          if (iboxShrink < 2)
            {
              MayDay::Error("wacky DBL generation technique failed, try making base box bigger");
            }
          coarBox.grow(-iboxShrink);
          coarBox.refine(a_params.m_refRatio[ilev-1]);
          Vector<Box> refBoxes(1, coarBox);
          a_coarBoxes[ilev] = refBoxes;
        }
    }
  else
    {
      pout() << "making grids at coarsest level by tagging irregular cells" << endl;
      pout() << "all grid parameters refer to coarsest hierarchy (max box size, etc)" << endl;
      //split up coarsest domain by max box size and
      //make a dbl at coarsest level
      Vector<Box> boxesCoarsest;
      Vector<int> procsCoarsest;

      domainSplit(a_coarsestDomain.domainBox(), boxesCoarsest,
                  a_params.m_maxBoxSize, a_params.m_blockFactor);

      a_coarBoxes.resize(numLevels);
      a_coarBoxes[0] = boxesCoarsest;
      if (numLevels > 1)
        {
          EBEllipticLoadBalance(procsCoarsest, boxesCoarsest, a_coarsestDomain);
          DisjointBoxLayout dblCoarsest(boxesCoarsest, procsCoarsest);

          //make a ebislayout at coarsest level.
          EBISLayout ebislCoarsest;
          a_ebisPtr->fillEBISLayout(ebislCoarsest,dblCoarsest,
                                  a_coarsestDomain.domainBox(), 0);

          //loop through grids and get all irregular cells into an intvectset.
          IntVectSet tagsCoarsestLocal;
          for (DataIterator dit = dblCoarsest.dataIterator(); dit.ok(); ++dit)
            {
              const EBISBox& ebisBox = ebislCoarsest[dit()];
              const Box&       box = dblCoarsest.get(dit());
              tagsCoarsestLocal |= ebisBox.getIrregIVS(box);
            }

          //generate vector of grids.
          BRMeshRefine gridder(a_coarsestDomain,     a_params.m_refRatio,
                               a_params.m_fillRatio, a_params.m_blockFactor,
                               a_params.m_tagBuffer, a_params.m_maxBoxSize);

          Vector<Vector<Box> > newMeshes(numLevels);
          Vector<Vector<Box> > oldMeshes(numLevels);
          oldMeshes[0]= boxesCoarsest;
          ProblemDomain curDomain = a_coarsestDomain;
          for (int ilev = 1; ilev < numLevels; ilev++)
            {
              domainSplit(curDomain.domainBox(), oldMeshes[ilev],
                          a_params.m_maxBoxSize, a_params.m_blockFactor);
              curDomain.refine( a_params.m_refRatio[ilev-1]);
            }
          int baseLevel = 0;
          gridder.regrid(newMeshes, tagsCoarsestLocal, baseLevel, a_params.m_maxLevel, oldMeshes);
          //   gridder.regrid(newMeshes, tagsCoarsest, baseLevel, a_params.m_maxLevel, oldMeshes);
          for (int ilev = 0; ilev < numLevels; ilev++)
            {
              a_coarBoxes[ilev] = newMeshes[ilev];
            }
        }
    }

  //now generate finer grids as refinements
  pout() << "generating finer hierarchies as refinements of coarsest hierarchy" << endl;
  a_mediBoxes = a_coarBoxes;
  a_fineBoxes = a_coarBoxes;
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      for (int ibox = 0;  ibox < a_coarBoxes[ilev].size(); ibox++)
        {
          a_mediBoxes[ilev][ibox].refine(2);
          a_fineBoxes[ilev][ibox].refine(4);
        }
    }


}
/****/
void
getAMRINSParameters(AMRParameters&   a_params,
                    ProblemDomain&   a_coarsestInputDomain)
{
  // read inputs
  ParmParse ppebamr;

  ppebamr.get("domain_length"      ,a_params.m_domainLength);
  ppebamr.get("max_level"          ,a_params.m_maxLevel);
  ppebamr.get("checkpoint_interval",a_params.m_checkpointInterval);
  ppebamr.get("plot_interval"      ,a_params.m_plotInterval);
  ppebamr.get("max_grid_size"      ,a_params.m_maxBoxSize);
  ppebamr.get("fill_ratio"         ,a_params.m_fillRatio);
  ppebamr.get("block_factor"       ,a_params.m_blockFactor);
  ppebamr.get("regrid_interval"    ,a_params.m_regridInterval);
  ppebamr.get("max_dt_grow"        ,a_params.m_maxDtGrow);
  ppebamr.query("max_dt"           ,a_params.m_maxDt);
  ppebamr.get("cfl"                ,a_params.m_cfl);
  ppebamr.get("init_cfl"           ,a_params.m_initCFL);
  ppebamr.get("refine_threshold"   ,a_params.m_refineThreshold);
  ppebamr.get("verbosity"          ,a_params.m_verbosity);
  ppebamr.get("nesting_radius"     ,a_params.m_nestingRadius);
  ppebamr.get("tag_buffer"         ,a_params.m_tagBuffer);
  ppebamr.query("tag_shrinkdomain" ,a_params.m_tagShrinkDomain);
  ppebamr.query("flow_dir"         ,a_params.m_flowDir);
  ppebamr.query("order_time"       ,a_params.m_orderTimeIntegration);
  int ilimit;
  ppebamr.get("use_limiting"         ,ilimit);
  a_params.m_useLimiting = (ilimit==1);

  int num_read_levels = a_params.m_maxLevel+1;
  ppebamr.getarr("ref_ratio",a_params.m_refRatio,0,num_read_levels);

  Vector<int> n_cell(SpaceDim);
  ppebamr.getarr("n_cell",n_cell,0,SpaceDim);
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

  if (ppebamr.contains("solver_tolerance"))
    {
      ppebamr.get("solver_tolerance", a_params.m_tolerance);
    }
  if (ppebamr.contains("solver_mgcycle"))
    {
      ppebamr.get("solver_mgcycle", a_params.m_mgCycle);
    }
  if (ppebamr.contains("solver_iter_max"))
    {
      ppebamr.get("solver_iter_max", a_params.m_iterMax);
    }
  if (ppebamr.contains("solver_hang"))
    {
      ppebamr.get("solver_hang", a_params.m_hang);
    }
  if (ppebamr.contains("solver_num_smooth"))
    {
      ppebamr.get("solver_num_smooth", a_params.m_numSmooth);
    }

  if (ppebamr.contains("domain_periodicity"))
    {
      Vector<int> periodicity(SpaceDim);
      ppebamr.getarr("domain_periodicity", periodicity, 0, SpaceDim);
      bool periodicityArr[SpaceDim];
      for (int idir=0; idir<SpaceDim; idir++)
        {
          periodicityArr[idir] = periodicity[idir]==1;
        }
      a_coarsestInputDomain = ProblemDomain(lo, hi, periodicityArr);
    }
  else
    {
      a_coarsestInputDomain = ProblemDomain(lo, hi);
    }
}

/******/
int iscript(int icomp, int inorm, int ncomp)
{
  return icomp + inorm*ncomp;
}

//kallemov overload function to use output file
/**************************/
void compareError(const Vector< LevelData<EBCellFAB>* >&   a_errorFine,
                  const Vector< LevelData<EBCellFAB>* >&   a_errorCoar,
                  const Vector< DisjointBoxLayout >&       a_gridsFine,
                  const Vector< DisjointBoxLayout >&       a_gridsCoar,
                  const Vector< EBISLayout >&              a_ebislFine,
                  const Vector< EBISLayout >&              a_ebislCoar,
                  const ProblemDomain&                     a_level0DomainFine,
                  const ProblemDomain&                     a_level0DomainCoar,
                  const Vector<string>&                    a_names,
                  const string&                            a_testName,
                  const AMRParameters&                     a_params,
                  fstream*                                 a_fout)
{
  const Vector<int> refRat = a_params.m_refRatio;

  ParmParse pp;
  bool ignoreOutFlow = false;
  if (pp.contains("compare_ignore_outflow"))
    {
      pp.get("compare_ignore_outflow",ignoreOutFlow);
    }
  int ncomp = a_names.size();
  Vector< LevelData<EBCellFAB>* > errorCoar = a_errorCoar;
  Vector< LevelData<EBCellFAB>* > errorFine = a_errorFine;
  int nlev  = a_errorCoar.size();
  if (ignoreOutFlow)
    {

      Box domLevFine =  a_level0DomainFine.domainBox();
      Box domLevCoar =  a_level0DomainCoar.domainBox();
      for (int ilev = 0; ilev < nlev ; ilev++)
        {
          //now for the actual zeroing
          int ncellZero = 2;
          if (a_fout == NULL)
          {
            pout() << "compareError: ignoring last " << ncellZero << " outflow cells " << endl;
          }
          else
          {
            (*a_fout) << "compareError: ignoring last " << ncellZero << " outflow cells " << endl;
          }
          Box shrunkDomainFine = domLevFine;
          Box shrunkDomainCoar = domLevCoar;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              shrunkDomainFine.growHi(idir, -ncellZero);
              shrunkDomainCoar.growHi(idir, -ncellZero);
            }

          for (DataIterator ditFine = a_gridsFine[ilev].dataIterator(); ditFine.ok(); ++ditFine)
            {
              const Box& boxFine = a_gridsFine[ilev].get(ditFine());
              IntVectSet ivs(boxFine);
              ivs -= shrunkDomainFine;
              for (VoFIterator vofit(ivs, a_ebislFine[ilev][ditFine()].getEBGraph()); vofit.ok(); ++vofit)
                {

                  for (int ivar = 0; ivar < ncomp; ivar++)
                    {
                      (*errorFine[ilev])[ditFine()](vofit(), ivar) = 0.0;
                    }
                }
            }
          for (DataIterator ditCoar = a_gridsCoar[ilev].dataIterator(); ditCoar.ok(); ++ditCoar)
            {
              const Box& boxCoar = a_gridsCoar[ilev].get(ditCoar());
              IntVectSet ivs(boxCoar);
              ivs -= shrunkDomainCoar;
              for (VoFIterator vofit(ivs, a_ebislCoar[ilev][ditCoar()].getEBGraph()); vofit.ok(); ++vofit)
                {
                  for (int ivar = 0; ivar < ncomp; ivar++)
                    {
                      (*errorCoar[ilev])[ditCoar()](vofit(), ivar) = 0.0;
                    }
                }
            }

          refine(domLevFine, refRat[ilev]);
          refine(domLevCoar, refRat[ilev]);
        }
    }
  Vector<Real> orders;
  int testverbosity;
  pp.get("test_verbosity", testverbosity);
  EBArith::compareError(orders,
                        a_errorFine,
                        a_errorCoar,
                        a_gridsFine,
                        a_gridsCoar,
                        a_ebislFine,
                        a_ebislCoar,
                        refRat,
                        a_level0DomainCoar.domainBox(),
                        testverbosity,
                        a_fout);

}

// Needed for sphere packing geometry
typedef struct
{
  int index;
  Real max;
} MAXHEIGHT;

// Needed for sphere packing geometry
int compMax(const void* a_p1, const void* a_p2)
{
  const MAXHEIGHT* e1 = (MAXHEIGHT*)a_p1;
  const MAXHEIGHT* e2 = (MAXHEIGHT*)a_p2;

  Real diff = e1->max - e2->max;

  if (diff < 0)
  {
    return -1;
  }
  else if (diff > 0)
  {
    return 1;
  }

  return 0;
}

/************/
void
AMRINSGeometry(const AMRParameters&   a_params,
               const ProblemDomain&   a_coarsestDomain,
               int a_whichGeomForced,
               EBIndexSpace* a_ebisPtr )
{
  CH_TIME("AMRINSGeometry");
  Real coarsestDx     = a_params.m_domainLength/Real(a_coarsestDomain.size(0));

  int max_level = a_params.m_maxLevel;

  ProblemDomain finestDomain = a_coarsestDomain;
  for (int ilev = 0; ilev < max_level; ilev++)
    {
      finestDomain.refine(a_params.m_refRatio[ilev]);
    }

  Real fineDx = coarsestDx;
  int ebMaxCoarsen = -1;
  for (int ilev = 0; ilev < max_level; ilev++)
    {
      fineDx /= a_params.m_refRatio[ilev];
    }

  ParmParse pp;
  RealVect origin;
  Vector<Real> originVect(SpaceDim, 0.0);
  pp.queryarr("domain_origin", originVect, 0, SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      origin[idir] = originVect[idir];
    }

  int whichgeom;
  if (a_whichGeomForced < 0)
    {
      pp.get("which_geom",whichgeom);
    }
  else
    {
      whichgeom = a_whichGeomForced;
    }
  int ebMaxSize = a_params.m_maxBoxSize;
  if (!pp.contains("ebis_file"))
    {
      int verbosity = 0;
      if (whichgeom == 0)
        {
          //allregular
          pout() << "all regular geometry" << endl;

          AllRegularService regserv;
          a_ebisPtr->define(finestDomain, origin, fineDx, regserv, ebMaxSize);
        }
      else if (whichgeom == 1)
        {
          pout() << "ramp geometry" << endl;
          RealVect rampNormal;
          vector<Real>  rampNormalVect(SpaceDim);
          pp.getarr("ramp_normal",rampNormalVect, 0, SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              rampNormal[idir] = rampNormalVect[idir];
            }

          Real rampAlpha;
          pp.get("ramp_alpha",rampAlpha);

          RealVect rampPoint = RealVect::Zero;
          for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (rampNormal[idir] != 0.0)
            {
              rampPoint[idir] = rampAlpha / rampNormal[idir];
              break;
            }
          }

          bool inside = true;

          PlaneIF ramp(rampNormal,rampPoint,inside);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(ramp,verbosity,vectDx);
          //this generates the new EBIS
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
        }
      else if (whichgeom == 2)
        {
          pout() << "slab geometry" << endl;
          vector<int> slab_lo(SpaceDim);
          pp.getarr("slab_lo",slab_lo,0,SpaceDim);
          vector<int> slab_hi(SpaceDim);
          pp.getarr("slab_hi",slab_hi,0,SpaceDim);
          IntVect lo, hi;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              lo[idir] = slab_lo[idir];
              hi[idir] = slab_hi[idir];
            }
          Box coveredBox(lo,hi);
          SlabService slab(coveredBox);
          //this generates the new EBIS
          RealVect origin = RealVect::Zero;
          a_ebisPtr->define(finestDomain, origin, fineDx, slab, ebMaxSize);
        }
      else if (whichgeom == 3)
        {
          pout() << "multi slab geometry" << endl;
          int numSlabs;
          pp.get("num_slabs",numSlabs);
          Vector<Box> coveredBoxes(numSlabs);
          for (int islab=0; islab<numSlabs; islab++)
            {
              char slabLoString[80];
              char slabHiString[80];
              sprintf(slabLoString, "slab_lo_%d", islab);
              sprintf(slabHiString, "slab_hi_%d", islab);
              vector<int> slab_lo(SpaceDim);
              pp.getarr(slabLoString,slab_lo,0,SpaceDim);
              vector<int> slab_hi(SpaceDim);
              pp.getarr(slabHiString,slab_hi,0,SpaceDim);
              IntVect lo, hi;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  lo[idir] = slab_lo[idir];
                  hi[idir] = slab_hi[idir];
                }
              Box coveredBox(lo,hi);
              coveredBoxes.push_back(coveredBox);
            }
          MultiSlabService multiSlab(coveredBoxes);
          //this generates the new EBIS
          RealVect origin = RealVect::Zero;
          a_ebisPtr->define(finestDomain, origin, fineDx, multiSlab, ebMaxSize);
        }
      else if (whichgeom == 4)
        {

          RealVect cylinderAxis, cylinderCenterPt;
          vector<Real>  cylinderAxisVect(SpaceDim);
          vector<Real>  cylinderCentVect(SpaceDim);
          pp.getarr("cylinder_axis",     cylinderAxisVect, 0, SpaceDim);
          pp.getarr("cylinder_center_pt",cylinderCentVect, 0, SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              cylinderAxis[idir]     = cylinderAxisVect[idir];
              cylinderCenterPt[idir] = cylinderCentVect[idir];
            }
          Real sum;
          PolyGeom::unifyVector(cylinderAxis, sum);

          Real cylinderRadius;
          pp.get("cylinder_radius",cylinderRadius);

          pout() << "using a tilted cylinder implicit function" << endl;
          bool cylinderInside;
          pp.get("cylinder_inside",cylinderInside);
//           bool negativeInside = true;
          TiltedCylinderIF tunnel(cylinderRadius, cylinderAxis, cylinderCenterPt, cylinderInside);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(tunnel,verbosity,vectDx);
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
        }
      else if (whichgeom == 5)
        {
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
          bool fluidInside = false;
          pp.query("fluid_inside", fluidInside);
          SphereIF sphere(sphereRadius, sphereCenter, fluidInside);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(sphere,verbosity,vectDx);
          //this generates the new EBIS
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
        }
      else if (whichgeom == 6)
        {
          pout() << "multisphere geometry" << endl;
          int numSpheres;
          pp.get("num_spheres", numSpheres);
          Vector<Real>     radius(numSpheres);
          Vector<RealVect> center(numSpheres);
          for (int isphere = 0; isphere < numSpheres; isphere++)
            {
              char radiusString[80];
              char centerString[80];
              sprintf(radiusString, "sphere_radius_%d", isphere);
              sprintf(centerString, "sphere_center_%d", isphere);
              vector<Real> sphere_center(SpaceDim);
              Real sphereRadius;
              pp.get(radiusString, sphereRadius);
              pp.getarr(centerString,sphere_center, 0, SpaceDim);
              RealVect sphereCenter;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  sphereCenter[idir] = sphere_center[idir];
                }
              center[isphere] = sphereCenter;
              radius[isphere] = sphereRadius;
            }

          bool insideRegular = false;
          pp.query("inside",insideRegular);

          MultiSphereIF impMultiSphere(radius,center,insideRegular);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(impMultiSphere,verbosity,vectDx);

          //this generates the new EBIS
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 7)
        {
          pout() << "multiparabola geometry" << endl;
          int numParabola;
          pp.get("num_parabolas", numParabola);
          int updir;
          pp.get("parabola_updir", updir);
          Vector<Real>     amp(numParabola);
          Vector<RealVect> center(numParabola);
          for (int iparabola = 0; iparabola < numParabola; iparabola++)
            {
              char ampString[80];
              char centerString[80];
              sprintf(ampString, "parabola_amplitude_%d", iparabola);
              sprintf(centerString, "parabola_center_%d", iparabola);
              vector<Real> parabola_center(SpaceDim);
              Real parabolaAmp;
              pp.get(ampString, parabolaAmp);
              pp.getarr(centerString,parabola_center, 0, SpaceDim);
              RealVect parabolaCenter;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  parabolaCenter[idir] = parabola_center[idir];
                }
              center[iparabola] = parabolaCenter;
              amp[iparabola]    = parabolaAmp;
            }

          Vector<BaseIF*> parabolas;
          for (int iparabola = 0; iparabola < numParabola; iparabola++)
          {
            Vector<PolyTerm> poly;

            PolyTerm mono;
            Real coef;
            IntVect powers;

            if (updir != 0)
            {
              // x^2 term
              coef = amp[iparabola];
              powers = IntVect::Zero;
              powers[0] = 2;

              mono.coef   = coef;
              mono.powers = powers;

              poly.push_back(mono);

              // x term
              coef = -2.0*amp[iparabola]*center[iparabola][0];
              powers = IntVect::Zero;
              powers[0] = 1;

              mono.coef   = coef;
              mono.powers = powers;

              poly.push_back(mono);

              // y or z term
              coef = -1.0;
              powers = IntVect::Zero;
              powers[updir] = 1;

              mono.coef   = coef;
              mono.powers = powers;

              poly.push_back(mono);

              // constant
              coef = amp[iparabola]*center[iparabola][0]*center[iparabola][0] + center[iparabola][updir];
              powers = IntVect::Zero;

              mono.coef   = coef;
              mono.powers = powers;

              poly.push_back(mono);
            }
            else
            {
              // y^2 term
              coef = amp[iparabola];
              powers = IntVect::Zero;
              powers[1] = 2;

              mono.coef   = coef;
              mono.powers = powers;

              poly.push_back(mono);

              // y term
              coef = -2.0*amp[iparabola]*center[iparabola][1];
              powers = IntVect::Zero;
              powers[1] = 1;

              mono.coef   = coef;
              mono.powers = powers;

              poly.push_back(mono);

              // x term
              coef = -1.0;
              powers = IntVect::Zero;
              powers[updir] = 1;

              mono.coef   = coef;
              mono.powers = powers;

              poly.push_back(mono);

              // constant
              coef = amp[iparabola]*center[iparabola][1]*center[iparabola][1] + center[iparabola][updir];
              powers = IntVect::Zero;

              mono.coef   = coef;
              mono.powers = powers;

              poly.push_back(mono);
            }

            if (amp[iparabola] < 0)
            {
              parabolas.push_back(new PolynomialIF(poly,true));
            }
            else
            {
              parabolas.push_back(new PolynomialIF(poly,false));
            }
          }

          IntersectionIF allTogether(parabolas);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(allTogether,verbosity,vectDx);
          //this generates the new EBIS
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
        }
      else if (whichgeom == 8)
        {
          pout() << "parabolic mirror geometry" << endl;
          Real amplitude;
          RealVect center;
          vector<Real> centervec;
          int updir;
          pp.get("mirror_updir", updir);
          pp.get("mirror_amplitude", amplitude);
          pp.getarr("mirror_center",centervec, 0, SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              center[idir] = centervec[idir];
            }

          Vector<PolyTerm> poly;

          PolyTerm mono;
          Real coef;
          IntVect powers;

          if (updir != 0)
          {
            // x^2 term
            coef = amplitude;
            powers = IntVect::Zero;
            powers[0] = 2;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // x term
            coef = -2.0*amplitude*center[0];
            powers = IntVect::Zero;
            powers[0] = 1;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // y or z term
            coef = -1.0;
            powers = IntVect::Zero;
            powers[updir] = 1;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // constant
            coef = amplitude*center[0]*center[0] + center[updir];
            powers = IntVect::Zero;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);
          }
          else
          {
            // y^2 term
            coef = amplitude;
            powers = IntVect::Zero;
            powers[1] = 2;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // y term
            coef = -2.0*amplitude*center[1];
            powers = IntVect::Zero;
            powers[1] = 1;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // x term
            coef = -1.0;
            powers = IntVect::Zero;
            powers[updir] = 1;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);

            // constant
            coef = amplitude*center[1]*center[1] + center[updir];
            powers = IntVect::Zero;

            mono.coef   = coef;
            mono.powers = powers;

            poly.push_back(mono);
          }

          bool inside = (amplitude >= 0);

          PolynomialIF mirror(poly,inside);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(mirror,verbosity,vectDx);
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
        }
      else if (whichgeom == 9)
        {
          pout() << "ovoid geometry" << endl;
          Real radius, axialdist;
          RealVect center;
          vector<Real> centervec;
          pp.get("ovoid_radius", radius);
          pp.get("ovoid_axial_dist", axialdist);
          int axis;
          pp.get("ovoid_axis", axis);
          pp.getarr("ovoid_center_lo",centervec, 0, SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              center[idir] = centervec[idir];
            }

          RealVect centerLo = center;

          SphereIF sphereLo(radius,centerLo,true);

          RealVect centerHi = center;
          centerHi[axis] += axialdist;

          SphereIF sphereHi(radius,centerHi,true);

          UnionIF ends(sphereLo,sphereHi);

          RealVect cylinderAxis = RealVect::Zero;
          cylinderAxis[axis] = 1.0;

          TiltedCylinderIF tube(radius,cylinderAxis,centerLo,true);
          PlaneIF loBound(cylinderAxis,centerLo,true);
          PlaneIF hiBound(cylinderAxis,centerHi,false);

          IntersectionIF loHi(loBound,hiBound);
          IntersectionIF finiteTube(tube,loHi);

          UnionIF capsule(finiteTube,ends);

          ComplementIF ovoid(capsule,true);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(ovoid,verbosity,vectDx);
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
        }
      else if (whichgeom == 10)
        {
          pout() << "Using MEMS contraction channel geometry" << endl;

          Real a_w,        // width of contracted channel
            a_R,           // radius of curvature of "large" contraction
            a_r,           // radius of curvature of "small" contraction
            a_d;           // depth of channel
          RealVect a_c;    // center of geometry

          pp.get("w",a_w);
          pp.get("R",a_R);
          pp.get("r",a_r);
          pp.get("d",a_d);

          // ParmParse doesn't get RealVects, so work-around with Vector<Real>
          Vector<Real> center;
          pp.getarr("c",center,0,SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              a_c[idir] = center[idir];
            }

          // These are needed to create cylinder and plane objects
          RealVect cylaxis;
          RealVect planenormal;
          RealVect point;
#if (CH_SPACEDIM==2)
          // create top of triangle
          point[0] = a_c[0];
          point[1] = a_c[1];
          planenormal[0] = -1;
          planenormal[1] = -1.732050807568877;
          PlaneIF side1(planenormal,point,true);
          // create bottom of triangle
          planenormal[1] *= -1;
          PlaneIF side2(planenormal,point,true);
          // create the triangle
          IntersectionIF triangle(side1,side2);

          // create plane to chop tip of triangle
          point[0] = a_c[0] - 1.5*a_R;
          point[1] = a_c[1];
          planenormal[0] = -1;
          planenormal[1] = 0;
          PlaneIF side3(planenormal,point,true);
          // create triangle with chopped end
          IntersectionIF choppedTriangle(triangle,side3);

          // create large contracting circle
          point[0] = a_c[0] - 2.0*a_R;
          point[1] = a_c[1];
          SphereIF largeContraction(a_R,point,true);
          // create triangle with rounded end
          UnionIF roundedTriangle(choppedTriangle,largeContraction);

          // create top of contracted channel
          point[0] = a_c[0];
          point[1] = a_c[1] + 0.5*a_w;
          planenormal[0] = 0;
          planenormal[1] = -1;
          PlaneIF topChannel(planenormal,point,true);
          // create bottom of contracted channel
          point[0] = a_c[0];
          point[1] = a_c[1] - 0.5*a_w;
          planenormal[0] = 0;
          planenormal[1] = 1;
          PlaneIF bottomChannel(planenormal,point,true);
          // create contracted channel
          IntersectionIF smallChannel(topChannel,bottomChannel);

          // create the entire channel
          UnionIF entireUnroundedChannel(smallChannel,roundedTriangle);

          // now created rounded portions.
          // we are going to intersect a triangle with the intersection of
          // two not circles

          point[0] = a_c[0] - 2.0*a_R + sqrt(pow((a_R+a_r),2)-pow((0.5*a_w+a_r),2));
          point[1] = a_c[1] + 0.5*a_w + a_r;
          SphereIF roundingCircle1(a_r,point,false);

          point[0] = a_c[0] - 2.0*a_R + sqrt(pow((a_R+a_r),2)-pow((0.5*a_w+a_r),2));
          point[1] = a_c[1] - 0.5*a_w - a_r;
          SphereIF roundingCircle2(a_r,point,false);

          IntersectionIF notCircles(roundingCircle1,roundingCircle2);

          point[0] = a_c[0] - 2.0*a_R;
          point[1] = a_c[1];
          planenormal[0] = 0.5*a_w + a_r;
          planenormal[1] = -sqrt(pow((a_R+a_r),2)-pow((0.5*a_w+a_r),2));
          PlaneIF innerTri1(planenormal,point,true);

          point[0] = a_c[0] - 2.0*a_R;
          point[1] = a_c[1];
          planenormal[0] = 0.5*a_w + a_r;
          planenormal[1] = sqrt(pow((a_R+a_r),2)-pow((0.5*a_w+a_r),2));
          PlaneIF innerTri2(planenormal,point,true);

          point[0] = a_c[0] - 2.0*a_R + sqrt(pow((a_R+a_r),2)-pow((0.5*a_w+a_r),2));
          point[1] = a_c[1];
          planenormal[0] = -1.0;
          planenormal[1] = 0.0;
          PlaneIF innerTri3(planenormal,point,true);

          IntersectionIF innerTri12(innerTri1,innerTri2);
          IntersectionIF innerTri123(innerTri12,innerTri3);
          IntersectionIF roundedPart(innerTri123,notCircles);

          UnionIF entireChannel(roundedPart,entireUnroundedChannel);

#elif (CH_SPACEDIM==3)
          // create top of triangle
          point[0] = a_c[0];
          point[1] = a_c[1];
          point[2] = a_c[2];
          planenormal[0] = -1;
          planenormal[1] = -1.732050807568877;
          planenormal[2] = 0;
          PlaneIF side1(planenormal,point,true);
          // create bottom of triangle
          planenormal[1] *= -1;
          PlaneIF side2(planenormal,point,true);
          // create the triangle
          IntersectionIF triangle(side1,side2);

          // create plane to chop tip of triangle
          point[0] = a_c[0] - 1.5*a_R;
          point[1] = a_c[1];
          point[2] = a_c[2];
          planenormal[0] = -1;
          planenormal[1] = 0;
          planenormal[2] = 0;
          PlaneIF side3(planenormal,point,true);
          // create triangle with chopped end
          IntersectionIF choppedTriangle(triangle,side3);

          // create large contracting circle
          point[0] = a_c[0] - 2.0*a_R;
          point[1] = a_c[1];
          point[2] = a_c[2];
          cylaxis[0] = 0;
          cylaxis[1] = 0;
          cylaxis[2] = 1;
          TiltedCylinderIF largeContraction(a_R,cylaxis,point,true);
          // create triangle with rounded end
          UnionIF roundedTriangle(choppedTriangle,largeContraction);

          // create top of contracted channel
          point[0] = a_c[0];
          point[1] = a_c[1] + 0.5*a_w;
          point[2] = a_c[2];
          planenormal[0] = 0;
          planenormal[1] = -1;
          planenormal[2] = 0;
          PlaneIF topChannel(planenormal,point,true);
          // create bottom of contracted channel
          point[0] = a_c[0];
          point[1] = a_c[1] - 0.5*a_w;
          point[2] = a_c[2];
          planenormal[0] = 0;
          planenormal[1] = 1;
          planenormal[2] = 0;
          PlaneIF bottomChannel(planenormal,point,true);
          // create contracted channel
          IntersectionIF smallChannel(topChannel,bottomChannel);

          // create the entire channel, infinite in z
          UnionIF entireUnroundedChannel(smallChannel,roundedTriangle);

          point[0] = a_c[0] - 2.0*a_R + sqrt(pow((a_R+a_r),2)-pow((0.5*a_w+a_r),2));
          point[1] = a_c[1] + 0.5*a_w + a_r;
          SphereIF roundingCircle1(a_r,point,false);

          point[0] = a_c[0] - 2.0*a_R + sqrt(pow((a_R+a_r),2)-pow((0.5*a_w+a_r),2));
          point[1] = a_c[1] - 0.5*a_w - a_r;
          SphereIF roundingCircle2(a_r,point,false);

          IntersectionIF notCircles(roundingCircle1,roundingCircle2);

          point[0] = a_c[0] - 2.0*a_R;
          point[1] = a_c[1];
          planenormal[0] = 0.5*a_w + a_r;
          planenormal[1] = -sqrt(pow((a_R+a_r),2)-pow((0.5*a_w+a_r),2));
          PlaneIF innerTri1(planenormal,point,true);

          point[0] = a_c[0] - 2.0*a_R;
          point[1] = a_c[1];
          planenormal[0] = 0.5*a_w + a_r;
          planenormal[1] = sqrt(pow((a_R+a_r),2)-pow((0.5*a_w+a_r),2));
          PlaneIF innerTri2(planenormal,point,true);

          point[0] = a_c[0] - 2.0*a_R + sqrt(pow((a_R+a_r),2)-pow((0.5*a_w+a_r),2));
          point[1] = a_c[1];
          planenormal[0] = -1.0;
          planenormal[1] = 0.0;
          PlaneIF innerTri3(planenormal,point,true);

          IntersectionIF innerTri12(innerTri1,innerTri2);
          IntersectionIF innerTri123(innerTri12,innerTri3);
          IntersectionIF roundedPart(innerTri123,notCircles);

          UnionIF entireUnboundedChannel(roundedPart,entireUnroundedChannel);

          // create top bounding plane
          point[0] = a_c[0];
          point[1] = a_c[1];
          point[2] = a_c[2] + a_d/2.0;
          planenormal[0] = 0;
          planenormal[1] = 0;
          planenormal[2] = -1;
          PlaneIF top(planenormal,point,true);
          // create bottom bounding plane
          point[0] = a_c[0];
          point[1] = a_c[1];
          point[2] = a_c[2] - a_d/2.0;
          planenormal[0] = 0;
          planenormal[1] = 0;
          planenormal[2] = +1;
          PlaneIF bottom(planenormal,point,true);
          // create bounding plane (in z)
          IntersectionIF topbottom(top,bottom);

          // create finished channel geometry
          IntersectionIF entireChannel(entireUnboundedChannel,topbottom);
#endif

    RealVect vectDx = RealVect::Unit;
    vectDx *= fineDx;

          GeometryShop geometryshop(entireChannel,verbosity,vectDx);
          a_ebisPtr->define(finestDomain,origin,fineDx,geometryshop,ebMaxSize);
        }
      else if (whichgeom == 11)
        {
          pout() << "ellipsoid geometry" << endl;
          vector<Real> ellipsoid_center(SpaceDim);
          vector<Real> ellipsoid_radii(SpaceDim);
          vector<Real> ellipsoid_xaxis(SpaceDim);
          pp.getarr("ellipsoid_center",ellipsoid_center, 0, SpaceDim);
          pp.getarr("ellipsoid_radii", ellipsoid_radii, 0, SpaceDim);
          pp.getarr("ellipsoid_xaxis", ellipsoid_xaxis, 0, SpaceDim);
          int fluid_inside;
          pp.get("ellipsoid_fluid_inside",fluid_inside);
          bool insideCalc = (fluid_inside==1);
          RealVect ellipsoidCenter;
          RealVect ellipsoidRadii;
          RealVect ellipsoidXAxis;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              ellipsoidCenter[idir] = ellipsoid_center[idir];
              ellipsoidRadii[idir]  = ellipsoid_radii[idir];
              ellipsoidXAxis[idir]  = ellipsoid_xaxis[idir];
            }
          Real sum;
          PolyGeom::unifyVector(ellipsoidXAxis, sum);
          EllipsoidIF ellipsoid(ellipsoidRadii, ellipsoidCenter, insideCalc);
          RealVect origxaxis = BASISREALV(0);
          TransformIF rotazoid(ellipsoid);
          rotazoid.rotate(origxaxis, ellipsoidXAxis, ellipsoidCenter);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(rotazoid,verbosity,vectDx);
          //this generates the new EBIS
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
        }
      else if (whichgeom == 12)
        {
          pout() << "DEMIF geometry" << endl;
          std::string demifFile;
          pp.get("demif_file",demifFile);
          IntVect demifNCell = finestDomain.size();
          int interpType;
          pp.get("demif_interp_type", interpType);
          RealVect dxVect = RealVect::Unit;
          dxVect *= fineDx;
          Real bottomBuffer, highGround, verticalScale;
          pp.get("demif_bottom_buffer", bottomBuffer);
          pp.get("demif_high_ground", highGround);
          pp.get("demif_vertical_scale", verticalScale);


          DEMIF demifoid(demifNCell, interpType, dxVect, demifFile,
                         bottomBuffer, 1.e99, highGround, verticalScale);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(demifoid,verbosity,vectDx);
          //this generates the new EBIS
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
        }
      else if (whichgeom == 13)
        {
          pout() << "tyler channel geometry" << endl;
          Real xsize = finestDomain.size(0);
          Real ysize = finestDomain.size(1);
          Real ratio = ysize/xsize;
          Real yDomainLength = a_params.m_domainLength*ratio;
          Real x1, x2, y1, y2;
          pp.get("tyler_x1", x1);
          pp.get("tyler_x2", x2);
          pp.get("tyler_y1", y1);
          pp.get("tyler_y2", y2);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          TylerChannelIF channel(x1, x2, y1, y2, yDomainLength);
          GeometryShop workshop(channel,verbosity,vectDx);
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
        }
      else if (whichgeom == 15)
        {
          pout() << "halfsphere geometry" << endl;
          vector<Real> sphere_center(SpaceDim);
          vector<Real> cutplane_xaxis(SpaceDim);
          pp.getarr("sphere_center",sphere_center, 0, SpaceDim);
          pp.getarr("cutplane_normal", cutplane_xaxis, 0, SpaceDim);
          RealVect cutplaneXAxis;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              cutplaneXAxis[idir]  = cutplane_xaxis[idir];
            }
          RealVect sphereCenter;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              sphereCenter[idir] = sphere_center[idir];
            }
          Real sphereRadius;
          pp.get("sphere_radius", sphereRadius);
          SphereIF sphere(sphereRadius, sphereCenter, true);
          PlaneIF  plane(cutplaneXAxis, sphereCenter, false);
          IntersectionIF intersect(plane, sphere);
          ComplementIF complement(intersect);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(complement,verbosity,vectDx);
          //this generates the new EBIS
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 16)
        {
          pout() << "sphere array" << endl;
          Real sphereRadius;
          pp.get("sphere_radius", sphereRadius);

          vector<Real> first_sphere_center(SpaceDim);
          vector<Real> sphere_spacing(SpaceDim);
          pp.getarr("first_sphere_center",first_sphere_center, 0, SpaceDim);
          pp.getarr("sphere_spacing",sphere_spacing, 0, SpaceDim);
          RealVect firstCenter, sphereSpacing;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              firstCenter[idir] = first_sphere_center[idir];
              sphereSpacing[idir] = sphere_spacing[idir];
            }
          SphereArrayIF implicit(sphereRadius, firstCenter, sphereSpacing);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(implicit,verbosity,vectDx);

          //this generates the new EBIS
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 17)
        {
          pout() << "Implicit function from datafile Geometry" << endl;

          int ebMaxCoarsen = -1;
          Real levelSetValue1;
          pp.get("level",levelSetValue1);
          // Parm Parse doesn't get bools, so work-around with int
          bool inside;
          int intInsideRegular;
          pp.get("insideRegular",intInsideRegular);

          if (intInsideRegular != 0) inside = true;
          if (intInsideRegular == 0) inside = false;

          Vector<Real> prob_lo(SpaceDim, 1.0);
          pp.getarr("origin",prob_lo,0,SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              origin[idir] = prob_lo[idir];
            }

          DataFileIF* dataFile1;
          if (pp.contains("input_file"))
            {
              std::string in_file;
              pp.get("input_file",in_file);
              dataFile1 = new DataFileIF(in_file.c_str(),DataFileIF::ASCII,levelSetValue1,inside);
            }
          else
            {
              dataFile1 = new DataFileIF(DataFileIF::ASCII,levelSetValue1,inside);
            }

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop geometryShop(*dataFile1,verbosity,vectDx);
          a_ebisPtr->
            define(finestDomain, origin, fineDx, geometryShop, ebMaxSize, ebMaxCoarsen);

          delete dataFile1;
        }
      else if (whichgeom == 18)
        {
          pout() << "Low swirl burner geometry" << endl;

          //need help resolving geometry because of issue when geometry scales ~ dx
          int finestDomainRequired = 256;
          int refineFactor = 1;
          pp.query("finest_domain",finestDomainRequired);
          if (fineDx > 1./finestDomainRequired && finestDomain.size(0) < finestDomainRequired)
            {
              refineFactor = finestDomainRequired/finestDomain.size(0);
              pout() << "additional refinement of finest domain by " << refineFactor << endl;
              finestDomain.refine(refineFactor);
              fineDx /= refineFactor;
            }
//               finestDomain.refine(4);
//               fineDx /= 4;

          Box domain;
          RealVect origin;
          Vector<int> n_cell(SpaceDim);
          pp.getarr("n_cell",n_cell,0,SpaceDim);

          CH_assert(n_cell.size() == SpaceDim);

          IntVect lo = IntVect::Zero;
          IntVect hi;

          for (int ivec = 0; ivec < SpaceDim; ivec++)
            {
              if (n_cell[ivec] <= 0)
                {
                  pout() << "Bogus number of cells input = " << n_cell[ivec];
                  exit(1);
                }

              hi[ivec] = n_cell[ivec] - 1;
            }

          domain.setSmall(lo);
          domain.setBig(hi);

          Vector<Real> prob_lo(SpaceDim,1.0);

          pp.getarr("origin",prob_lo,0,SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              origin[idir] = prob_lo[idir];
            }

          Real outerRadius;
          pp.get("outer_radius",outerRadius);

          Real outerThick;
          pp.get("outer_thick",outerThick);

          Real outerHeight;
          pp.get("outer_height",outerHeight);

          Real outerOffset;
          pp.get("outer_offset",outerOffset);
          outerOffset += prob_lo[0];

          Real innerRadius;
          pp.get("inner_radius",innerRadius);

          Real innerThick;
          pp.get("inner_thick",innerThick);

          Real innerOffset = 0.0;
          innerOffset += outerOffset;

          Real innerHeight;
          pp.get("inner_height",innerHeight);

          Real plateHeight;
          pp.get("plate_height",plateHeight);
          plateHeight += outerOffset;

          Real plateThick;
          pp.get("plate_thick",plateThick);

          int doHoles;
          pp.get("do_holes",doHoles);

          Real holeRadius;
          pp.get("hole_radius",holeRadius);

          Real holeSpace;
          pp.get("hole_space",holeSpace);

          int vaneNum;
          pp.get("vane_num",vaneNum);

          Real vaneThick;
          pp.get("vane_thick",vaneThick);

          RealVect vaneNorm;

          Vector<Real> vectVaneNorm;
          pp.getarr("vane_norm",vectVaneNorm,0,SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              vaneNorm[idir] = vectVaneNorm[idir];
            }

          Real vaneOffset;
          pp.get("vane_offset",vaneOffset);

          innerOffset += vaneOffset;

          Real vaneHeight = innerHeight - 2*vaneOffset;

          vaneOffset += outerOffset;

          // Make the outer chamber
          BaseIF* outerChamber = makeChamber(outerRadius,outerThick,
                                             outerOffset,outerHeight);

          // Make the inner chamber
          BaseIF* innerChamber = makeChamber(innerRadius,innerThick,
                                             innerOffset,innerHeight);

          // Make the inner plate with holes
          BaseIF* holyPlate = makePlate(plateHeight,plateThick,innerRadius,
                                        doHoles,holeRadius,holeSpace);

          // Make the vanes
          BaseIF* vanes = makeVanes(vaneNum,vaneThick,vaneNorm,innerRadius,outerRadius,
                                    vaneOffset,vaneHeight);

          // Union all the pieces together
          Vector<BaseIF*> pieces;

          pieces.push_back(outerChamber);
          pieces.push_back(innerChamber);
          pieces.push_back(holyPlate);
          pieces.push_back(vanes);

          UnionIF swirl(pieces);
          ComplementIF outside(swirl,true);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(outside,verbosity,vectDx);

          // This generates the new EBIS
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop);

        }
      else if (whichgeom == 19)
        {
          pout() << "multi-spheres inside a tilted cylinder geometry " << endl;
          int numSpheres;
          pp.get("num_spheres", numSpheres);
          Vector<Real>     radius(numSpheres);
          Vector<RealVect> center(numSpheres);
          for (int isphere = 0; isphere < numSpheres; isphere++)
            {
              char radiusString[80];
              char centerString[80];
              sprintf(radiusString, "sphere_radius_%d", isphere);
              sprintf(centerString, "sphere_center_%d", isphere);
              vector<Real> sphere_center(SpaceDim);
              Real sphereRadius;
              pp.get(radiusString, sphereRadius);
              pp.getarr(centerString,sphere_center, 0, SpaceDim);
              RealVect sphereCenter;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  sphereCenter[idir] = sphere_center[idir];
                }
              center[isphere] = sphereCenter;
              radius[isphere] = sphereRadius;
            }

          bool insideRegular = false;
          pp.query("insideRegular",insideRegular);

          Vector<BaseIF*> spheres(numSpheres);
          for (int i = 0; i < numSpheres; i++)
          {
            spheres[i] = new SphereIF(radius[i],center[i],insideRegular);
          }

          // Take the union of the insides of the spheres
          UnionIF insideSpheres(spheres);

          // Complement to get the outside of the spheres
          ComplementIF outsideSpheres(insideSpheres,true);

          RealVect cylinderAxis, cylinderCenterPt;
          vector<Real>  cylinderAxisVect(SpaceDim);
          vector<Real>  cylinderCentVect(SpaceDim);
          pp.getarr("cylinder_axis",     cylinderAxisVect, 0, SpaceDim);
          pp.getarr("cylinder_center_pt",cylinderCentVect, 0, SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              cylinderAxis[idir]     = cylinderAxisVect[idir];
              cylinderCenterPt[idir] = cylinderCentVect[idir];
            }
          Real sum;
          PolyGeom::unifyVector(cylinderAxis, sum);

          Real cylinderRadius;
          pp.get("cylinder_radius",cylinderRadius);

          bool cylinderInside = true;

          TiltedCylinderIF  cylinder(cylinderRadius, cylinderAxis, cylinderCenterPt, cylinderInside);

          IntersectionIF finalIF(cylinder, outsideSpheres);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(finalIF,verbosity,vectDx);
          //this generates the new EBIS
          pout() << "ebisPtr->define before" << endl;
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
          pout() << "ebisPtr->define after" << endl;
        }
      else if (whichgeom == 20)
        {
          pout() << "artery geometry " << endl;

          int arteryType;
          pp.get("artery_type",arteryType);

          Real arteryRadius;
          pp.get("artery_radius",arteryRadius);

          Real arteryPerturb;
          pp.get("artery_perturb",arteryPerturb);

          Real arteryStart;
          pp.get("artery_start",arteryStart);

          Real arteryEnd;
          pp.get("artery_end",arteryEnd);

          RealVect arteryCenter;
          Vector<Real> arteryCenterVect;
          pp.getarr("artery_center",arteryCenterVect,0,SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
          {
            arteryCenter[idir] = arteryCenterVect[idir];
          }

          bool arteryInside = true;

          Vector<Real> originVect(SpaceDim, 0.0);
          pp.getarr("origin",originVect,0,SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
          {
            origin[idir] = originVect[idir];
          }

          ArteryIF arteryIF(arteryType,arteryRadius,arteryPerturb,arteryStart,arteryEnd,arteryCenter,arteryInside);

          int clotPresent;
          pp.get("clot_present",clotPresent);

          RealVect clotCenter;
          RealVect clotRadii;
          bool clotInside = false;

          Real clotAngle1;
          Real clotAngle2;

          if (clotPresent != 0)
          {
            vector<Real> clot_center(SpaceDim);
            pp.getarr("clot_center",clot_center, 0, SpaceDim);
            for (int idir = 0; idir < SpaceDim; idir++)
              {
                clotCenter[idir] = clot_center[idir];
              }

            vector<Real> clot_radii(SpaceDim);
            pp.getarr("clot_radii",clot_radii, 0, SpaceDim);
            for (int idir = 0; idir < SpaceDim; idir++)
              {
                clotRadii[idir] = clot_radii[idir];
              }

            pp.get("clot_angle1",clotAngle1);
            clotAngle1 = M_PI * clotAngle1 / 180.0;

            pp.get("clot_angle2",clotAngle2);
            clotAngle2 = M_PI * clotAngle2 / 180.0;
          }

          BaseIF* finalIF;

          if (clotPresent == 0)
          {
            finalIF = &arteryIF;
          }
          else
          {
            EllipsoidIF clotIF(clotRadii,clotCenter,clotInside);

            TransformIF rotClotIF(clotIF);

            if (SpaceDim < 3)
            {
              rotClotIF.rotate(clotAngle1,clotCenter);
            }
            else
            {
              rotClotIF.rotate(clotAngle1,clotCenter,BASISREALV(2));
              rotClotIF.rotate(clotAngle2,clotCenter,BASISREALV(0));
            }

            finalIF = new IntersectionIF(arteryIF,rotClotIF);
          }

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(*finalIF,verbosity,vectDx);

          //this generates the new EBIS
          pout() << "ebisPtr->define before" << endl;
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
          pout() << "ebisPtr->define after" << endl;
        }
      else if (whichgeom == 21)
        {
          pout() << "cylindric step " << endl;

          Real inletRadius;
          pp.get("inlet_radius",inletRadius);

          Real outletRadius;
          pp.get("outlet_radius",outletRadius);

          RealVect center;
          Vector<Real> centerVect;
          pp.getarr("center",centerVect,0,SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
          {
            center[idir] = centerVect[idir];
          }

          bool useInside = true;

          Vector<Real> originVect(SpaceDim, 0.0);
          pp.getarr("origin",originVect,0,SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
          {
            origin[idir] = originVect[idir];
          }

          RealVect axis = BASISREALV(0);

          TiltedCylinderIF  inletCyl( inletRadius, axis, center, useInside);
          TiltedCylinderIF outletCyl(outletRadius, axis, center, useInside);

          BaseIF *finalInlet;
          BaseIF *finalOutlet;

          // Truncate the larger cylinder and union the result
          if (inletRadius < outletRadius)
          {
            finalInlet = &inletCyl;

            RealVect normal = BASISREALV(0);
            PlaneIF truncPlane(normal, center, useInside);

            finalOutlet = new IntersectionIF(outletCyl, truncPlane);
          }
          else
          {
            RealVect normal = -BASISREALV(0);
            PlaneIF truncPlane(normal, center, useInside);

            finalInlet = new IntersectionIF(inletCyl, truncPlane);

            finalOutlet = &outletCyl;
          }

          UnionIF cylindricStep(*finalInlet, *finalOutlet);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop workshop(cylindricStep,verbosity,vectDx);

          //this generates the new EBIS
          pout() << "ebisPtr->define before" << endl;
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
          pout() << "ebisPtr->define after" << endl;
        }
      else if (whichgeom == 22)
        {
          pout() << "Channel with spheres" << endl;

          RealVect cylDirection;
          RealVect cylPoint;
          Real cylRadius;
          RealVect cylCenter;
          RealVect cylAxis;
          int      sphNumber;
          Real     sphMinRadius;
          Real     sphMaxRadius;
          bool     insideRegular;

          RealVect domainLengthVect;

          Vector<Real> vectorDirection;
          pp.getarr("cylDirection",vectorDirection,0,SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              cylDirection[idir] = vectorDirection[idir];
            }

          Vector<Real> vectorPoint;
          pp.getarr("cylPoint",vectorPoint,0,SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              cylPoint[idir] = vectorPoint[idir];
            }

          pp.get("cylRadius",cylRadius);

          pp.get("sphNumber",sphNumber);

          pp.get("sphMinRadius",sphMinRadius);
          pp.get("sphMaxRadius",sphMaxRadius);

          // Parm Parse doesn't get bools, so work-around with int
          int intInsideRegular;
          pp.get("insideRegular",intInsideRegular);

          if (intInsideRegular != 0) insideRegular = true;
          if (intInsideRegular == 0) insideRegular = false;

          // Stay inside until the end
          bool inside = true;

          // The vector of spheres
          Vector<BaseIF*> spheres;
          spheres.resize(sphNumber);

          // Place random spheres inside the cylinder
          int i = 0;

          Vector<int> n_cell(SpaceDim);
          pp.getarr("n_cell",n_cell,0,SpaceDim);

          Vector<Real> prob_lo(SpaceDim,1.0);

          pp.getarr("origin",prob_lo,0,SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              domainLengthVect[idir] = n_cell[idir] * fineDx;
              origin[idir] = prob_lo[idir];
            }

          int randSeed;
          pp.get("randSeed",randSeed);

          srand48(randSeed);

          while (i < sphNumber)
          {
            RealVect center;
            Real     radius;

            // Get a random center and radius for a sphere
            for (int idir = 0; idir < SpaceDim; idir++)
              {
                center[idir] = (domainLengthVect[idir] - prob_lo[idir])*drand48() + prob_lo[idir];
              }

            radius = (sphMaxRadius - sphMinRadius)*drand48() + sphMinRadius;

            // Get the distance from the center of the sphere to the axis of the
            // cylinder
            Real dist = getDistance(center,cylDirection,cylPoint);

            // If the sphere is inside cylinder, use it
            if (dist + radius <= cylRadius)
            {
              spheres[i] = new SphereIF(radius,center,inside);
              i++;
            }
          }

          // Take the union of the insides of the spheres
          UnionIF insideSpheres(spheres);

          // Complement to get the outside of the spheres
          ComplementIF outsideSpheres(insideSpheres,true);

          // Define the cylinder
          TiltedCylinderIF cylinder(cylRadius,cylDirection,cylPoint,inside);

          // Intersect the inside of the cylinder with the outside of the spheres
          IntersectionIF implicit(cylinder,outsideSpheres);

          // Complement if necessary
          ComplementIF finalIF(implicit,!insideRegular);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;
          GeometryShop workshop(finalIF,verbosity,vectDx);

          // This generates the new EBIS
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
        }
      else if (whichgeom == 23)
        {
          pout() << "Implicit function from datafile Geometry intersected with tilted cylinder" << endl;

          int ebMaxCoarsen = -1;
          Real levelSetValue1;
          pp.get("level",levelSetValue1);
          // Parm Parse doesn't get bools, so work-around with int
          bool inside;
          int intInsideRegular;
          pp.get("insideRegular",intInsideRegular);

          if (intInsideRegular != 0) inside = true;
          if (intInsideRegular == 0) inside = false;

          Vector<Real> prob_lo(SpaceDim, 1.0);
          pp.getarr("origin",prob_lo,0,SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              origin[idir] = prob_lo[idir];
            }

          DataFileIF* dataFile1;
          if (pp.contains("input_file"))
            {
              std::string in_file;
              pp.get("input_file",in_file);
              dataFile1 = new DataFileIF(in_file.c_str(),DataFileIF::ASCII,levelSetValue1,inside);
            }
          else
            {
              dataFile1 = new DataFileIF(DataFileIF::ASCII,levelSetValue1,inside);
            }

          // Define the cylinder
          RealVect cylDirection;
          RealVect cylPoint;
          Real cylRadius;
          RealVect cylCenter;
          RealVect cylAxis;
          bool cylInside;

          Vector<Real> vectorDirection;
          pp.getarr("cylDirection",vectorDirection,0,SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              cylDirection[idir] = vectorDirection[idir];
            }

          Vector<Real> vectorPoint;
          pp.getarr("cylPoint",vectorPoint,0,SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              cylPoint[idir] = vectorPoint[idir];
            }

          pp.get("cylRadius",cylRadius);

          pp.get("cylInside",cylInside);

          TiltedCylinderIF cylinder(cylRadius,cylDirection,cylPoint,cylInside);

          // Intersect the inside of the cylinder with the outside of the spheres
          IntersectionIF core(cylinder,*dataFile1);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;

          GeometryShop geometryShop(core,verbosity,vectDx);
          a_ebisPtr->
            define(finestDomain, origin, fineDx, geometryShop, ebMaxSize, ebMaxCoarsen);

          delete dataFile1;
        }
      else if (whichgeom == 24)
        {
          pout() << "Fuel element array geometry" << endl;

          int ebMaxCoarsen = -1;

          // Get common cylinder attributes
          Vector<Real> tmp(SpaceDim);
          RealVect cylAxis;
          pp.getarr("cylinder_axis", tmp, 0, SpaceDim);
          for (int idir=0; idir<SpaceDim; idir++)
            {
              cylAxis[idir] = tmp[idir];
            }
          Real cylRadius;
          pp.get("cylinder_radius", cylRadius);
          bool cylOutside;
          pp.get("cylinder_outside", cylOutside);

          // Get each cylinder's center
          Vector<RealVect> cylCenters;
          char charstr[100];
          int icyl=0;
          sprintf(charstr, "cylinder_center%d", icyl);
          while (pp.contains(charstr))
            {
              pp.getarr(charstr, tmp, 0, SpaceDim);
              RealVect ctr;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  ctr[idir] = tmp[idir];
                }
              cylCenters.push_back(ctr);
              sprintf(charstr, "cylinder_center%d", ++icyl);
            }

          // Get coil attributes
          Real coilRadius;
          pp.get("coil_radius",coilRadius);
          Real coilSeparation;
          pp.get("coil_separation",coilSeparation);
          Real coilOverlapFraction;
          pp.get("coil_overlap_fraction",coilOverlapFraction);

          // create a cylinder at the origin
          TiltedCylinderIF oneCyl(cylRadius,
                                  cylAxis,
                                  RealVect::Zero,
                                  true);

          // calculate helicoil parameters and create it
          Real helixRadius = cylRadius + coilOverlapFraction*coilRadius;
          HelicoilIF coil(helixRadius, coilSeparation, coilRadius, true);

          // crop the element to the specified height (don't crop if zero)
          Real cylHeight=0;
          pp.query("cylinder_height",cylHeight);
          BaseIF* element0;
          if (cylHeight > 0)
            {
              RealVect    cylTop =  0.5*cylHeight*cylAxis;
              RealVect cylBottom = -0.5*cylHeight*cylAxis;
              Vector<BaseIF*> planes(2, NULL);
              planes[0] = (BaseIF*) new PlaneIF(cylAxis,    cylTop, false);
              planes[1] = (BaseIF*) new PlaneIF(cylAxis, cylBottom, true);

              Vector<BaseIF*> shapes(2, NULL);
              shapes[0] = (BaseIF*) new IntersectionIF(planes);
              shapes[1] = (BaseIF*) new UnionIF(coil, oneCyl);
              element0 = (BaseIF*) new IntersectionIF(shapes);

              for (int i=0; i<2; i++)
                {
                  delete planes[i];
                  delete shapes[i];
                }
            }
          else
            {
              // take the union of the coil and cylinder
              element0 = (BaseIF*) new UnionIF(coil, oneCyl);
            }

          // make a Vector of UnionIF* to hold the fuel elements
          int nElem = cylCenters.size();
          Vector<BaseIF*> fuelElements(nElem, NULL);

          // make each fuel element
          for (int ielem=0; ielem<nElem; ielem++)
            {
              // move the element to be coaxial with the cylinder
              TransformIF* tc = new TransformIF(*element0);
              tc->translate(cylCenters[ielem]);
              fuelElements[ielem] = (BaseIF*) tc;
            }

          UnionIF array(fuelElements);
          ComplementIF outside(array,cylOutside);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;
          GeometryShop workshop(outside,verbosity,vectDx);

          // This generates the new EBIS
          a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize, ebMaxCoarsen);

          delete element0;
          for (int ielem=0; ielem<nElem; ielem++)
            {
              delete fuelElements[ielem];
            }
        }
      else if (whichgeom == 25)
        {
          pout() << "Channel with non-overlapping spheres" << endl;

          RealVect cylDirection;
          RealVect cylPoint;
          Real     cylRadius;
          RealVect cylCenter;
          RealVect cylAxis;
          int      sphNumber;
          Real     sphMinRadius;
          Real     sphMaxRadius;
          bool     insideRegular;

          RealVect domainLengthVect;

          bool useDomainBoundary = false;
          pp.query("useDomainBoundary",useDomainBoundary);

          int packingType = 1;
          pp.query("packingType",packingType);

          Vector<Real> vectorDirection;
          pp.getarr("cylDirection",vectorDirection,0,SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              cylDirection[idir] = vectorDirection[idir];
            }

          Vector<Real> vectorPoint;
          pp.getarr("cylPoint",vectorPoint,0,SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              cylPoint[idir] = vectorPoint[idir];
            }

          pp.get("cylRadius",cylRadius);

          pp.get("sphNumber",sphNumber);

          MAXHEIGHT* sortedList = new MAXHEIGHT[sphNumber];

          pp.get("sphMinRadius",sphMinRadius);
          pp.get("sphMaxRadius",sphMaxRadius);

          // Parm Parse doesn't get bools, so work-around with int
          int intInsideRegular;
          pp.get("insideRegular",intInsideRegular);

          Real minDistCyl;
          pp.get("minDistCyl",minDistCyl);

          Real minDistSph;
          pp.get("minDistSph",minDistSph);

          Real minDistInlet;
          pp.get("minDistInlet",minDistInlet);

          Real minDistOutlet;
          pp.get("minDistOutlet",minDistOutlet);

          double maxTries;
          pp.get("maxTries",maxTries);

          double maxDrops;
          pp.get("maxDrops",maxDrops);

          if (intInsideRegular != 0) insideRegular = true;
          if (intInsideRegular == 0) insideRegular = false;

          // Stay inside until the end
          bool inside = true;

          Vector<int> n_cell(SpaceDim);
          pp.getarr("n_cell",n_cell,0,SpaceDim);

          Vector<Real> prob_lo(SpaceDim,1.0);

          pp.getarr("origin",prob_lo,0,SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              domainLengthVect[idir] = n_cell[idir] * fineDx;
              origin[idir] = prob_lo[idir];
            }

          Vector<RealVect> centerArray(sphNumber);
          Vector<Real>     radiusArray(sphNumber);

          int randSeed;
          pp.get("randSeed",randSeed);

          srand48(randSeed);

          // Place random spheres inside the cylinder
          int packed = 0;
          double count = 0;

          RealVect curDirection = cylDirection;

          if (curDirection.dotProduct(BASISREALV(0)) > 0)
          {
            curDirection = -curDirection;
          }

          curDirection /= curDirection.vectorLength();

          Real a = curDirection.dotProduct(curDirection);

          while (packed < sphNumber && count < maxTries)
          {
            if (packingType == 1)
            {
              RealVect& center = centerArray[packed];
              Real&     radius = radiusArray[packed];

              // Get a random center and radius for a sphere
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  center[idir] = domainLengthVect[idir]*drand48() + origin[idir];
                }

              radius = (sphMaxRadius - sphMinRadius)*drand48() + sphMinRadius;

              bool isOkay = true;

              // Get the distance from the center of the sphere to the axis of the
              // cylinder
              Real dist = getDistance(center,cylDirection,cylPoint);

              // If the sphere is inside cylinder, use it
              if (dist + radius > cylRadius - minDistCyl)
              {
                isOkay = false;
              }

              // Check against x domain boundaries
              if (center[0] - radius - minDistInlet < origin[0])
              {
                isOkay = false;
              }

              if (center[0] + radius + minDistOutlet > origin[0] + domainLengthVect[0])
              {
                isOkay = false;
              }

              if (isOkay)
              {
                // Check if this overlaps another sphere
                for (int i = 0; i < packed; i++)
                {
                  RealVect diff = center;
                  diff -= centerArray[i];

                  Real dist = diff.vectorLength();

                  if (dist <= radius + radiusArray[i] + minDistSph)
                  {
                    isOkay = false;
                    break;
                  }
                }
              }

              if (isOkay)
              {
                packed++;
              }

              count++;
            }
            else if (packingType == 2)
#if 1
            {
              RealVect& bestCenter = centerArray[packed];
              Real&     bestRadius = radiusArray[packed];
              bool      bestOkay = false;

              for (double tries = 0; tries < maxDrops; tries++)
              {
                RealVect center;
                Real     radius;

                center[0] = domainLengthVect[0] + origin[0];

                // Get a random center and radius for a sphere
                for (int idir = 1; idir < SpaceDim; idir++)
                  {
                    center[idir] = domainLengthVect[idir]*drand48() + origin[idir];
                  }

                radius = (sphMaxRadius - sphMinRadius)*drand48() + sphMinRadius;

                bool isOkay = true;

                // Get the distance from the center of the sphere to the axis of the
                // cylinder
                Real dist = getDistance(center,cylDirection,cylPoint);

                // If the sphere is inside cylinder, use it
                if (dist + radius > cylRadius - minDistCyl)
                {
                  isOkay = false;
                }

                if (isOkay)
                {
                  Real minS = 0;
                  bool foundMin = false;
                  Real minCenterX = origin[0];

                  // Compute the minimum distance it can "fall"
                  for (int i = packed-1; i >= 0; i--)
                  {
                    if (foundMin && (sortedList[i].max < (minCenterX - radius)))
                    {
                      break;
                    }

                    RealVect& curCenter = centerArray[sortedList[i].index];
                    Real&     curRadius = radiusArray[sortedList[i].index];

                    RealVect diff = center;
                    diff -= curCenter;

                    dist = radius + curRadius + minDistSph;

                    Real b = 2*curDirection.dotProduct(diff);
                    Real c = diff.dotProduct(diff) - dist*dist;

                    Real discrim = b*b - 4*a*c;

                    if (discrim > 0)
                    {
                      Real s1 = (-b - sqrt(discrim))/(2*a);

                      if (s1 > 0 && (!foundMin || s1 < minS))
                      {
                        minS = s1;
                        foundMin = true;
                        minCenterX = center[0] + minS * curDirection[0];
                      }
                    }
                  }

                  if (!foundMin)
                  {
                    minS = (origin[0] + radius + minDistInlet - center[0]) / curDirection[0];
                    minCenterX = center[0] + minS * curDirection[0];
                  }

                  center += minS * curDirection;
                }

                // Check against high x domain boundary
                if (center[0] + radius + minDistOutlet > origin[0] + domainLengthVect[0])
                {
                  isOkay = false;
                }

                if (isOkay)
                {
                  if (bestOkay)
                  {
                    if (center[0] < bestCenter[0])
                    {
                      bestCenter = center;
                      bestRadius = radius;
                    }
                  }
                  else
                  {
                    bestCenter = center;
                    bestRadius = radius;

                    bestOkay   = true;
                  }
                }

                count++;
              }

              if (bestOkay)
              {
                sortedList[packed].index = packed;
                sortedList[packed].max   = centerArray[packed][0] + radiusArray[packed];

                packed++;

                qsort(sortedList,packed,sizeof(sortedList[0]),compMax);
              }
            }
#else
            {
              RealVect& bestCenter = centerArray[packed];
              Real&     bestRadius = radiusArray[packed];
              bool      bestOkay = false;

              for (int tries = 0; tries < maxDrops; tries++)
              {
                RealVect center;
                Real     radius;

                center[0] = domainLengthVect[0] + origin[0];

                // Get a random center and radius for a sphere
                for (int idir = 1; idir < SpaceDim; idir++)
                  {
                    center[idir] = domainLengthVect[idir]*drand48() + origin[idir];
                  }

                radius = (sphMaxRadius - sphMinRadius)*drand48() + sphMinRadius;

                bool isOkay = true;

                // Get the distance from the center of the sphere to the axis of the
                // cylinder
                Real dist = getDistance(center,cylDirection,cylPoint);

                // If the sphere is inside cylinder, use it
                if (dist + radius > cylRadius - minDistCyl)
                {
                  isOkay = false;
                }

                if (isOkay)
                {
                  RealVect curDirection = cylDirection;

                  if (curDirection.dotProduct(BASISREALV(0)) > 0)
                  {
                    curDirection = -curDirection;
                  }

                  Real minS = 0;
                  bool foundMin = false;

                  // Compute the minimum distance it can "fall"
                  for (int i = 0; i < packed; i++)
                  {
                    RealVect diff = center;
                    diff -= centerArray[i];

                    dist = radius + radiusArray[i] + minDistSph;

                    Real a = curDirection.dotProduct(curDirection);
                    Real b = 2*curDirection.dotProduct(diff);
                    Real c = diff.dotProduct(diff) - dist*dist;

                    Real discrim = b*b - 4*a*c;

                    if (discrim > 0)
                    {
                      Real s1 = (-b - sqrt(discrim))/(2*a);

                      if (s1 > 0 && (!foundMin || s1 < minS))
                      {
                        minS = s1;
                        foundMin = true;
                      }
                    }
                  }

                  if (!foundMin)
                  {
                    minS = (origin[0] + radius + minDistInlet - center[0]) / curDirection[0];
                  }

                  center += minS * curDirection / curDirection.vectorLength();
                }

                // Check against high x domain boundary
                if (center[0] + radius + minDistOutlet > origin[0] + domainLengthVect[0])
                {
                  isOkay = false;
                }

                if (isOkay)
                {
                  if (bestOkay)
                  {
                    if (center[0] < bestCenter[0])
                    {
                      bestCenter = center;
                      bestRadius = radius;
                    }
                  }
                  else
                  {
                    bestCenter = center;
                    bestRadius = radius;

                    bestOkay   = true;
                  }
                }

                count++;
              }

              if (bestOkay)
              {
                packed++;
              }
            }
#endif
            else
            {
              MayDay::Abort("Unknown packing type");
            }

            if (remainder(count,1000000) == 0)
            {
              pout() << "Tried " << count
                     << ", packed " << packed
                     << ", " << sphNumber-packed << " left"
                     << endl;
            }
          }

          if (packed == 0)
          {
            MayDay::Abort("No spheres could be packed in channel");
          }

          // The vector of spheres
          Vector<BaseIF*> spheres(packed);

          for (int i = 0; i < packed; i++)
          {
            spheres[i] = new SphereIF(radiusArray[i],centerArray[i],inside);
          }

          // Take the union of the insides of the spheres
          UnionIF insideSpheres(spheres);

          // Complement to get the outside of the spheres
          ComplementIF outsideSpheres(insideSpheres,true);

          // Define the cylinder
          TiltedCylinderIF cylinder(cylRadius,cylDirection,cylPoint,inside);

          // Intersect the inside of the cylinder with the outside of the spheres
          IntersectionIF implicit(cylinder,outsideSpheres);

          // Complement if necessary
          ComplementIF finalIF(implicit,!insideRegular);

          RealVect vectDx = RealVect::Unit;
          vectDx *= fineDx;
          if (useDomainBoundary)
            {
              GeometryShop workshop(outsideSpheres,verbosity,vectDx);

              // This generates the new EBIS
              a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
            }
          else
            {
              GeometryShop workshop(finalIF,verbosity,vectDx);

              // This generates the new EBIS
              a_ebisPtr->define(finestDomain, origin, fineDx, workshop, ebMaxSize);
            }

          delete [] sortedList;
        }
      else
        {
          //bogus which_geom
          pout() << " bogus which_geom input = "
                 << whichgeom << endl;
          MayDay::Error();
        }
    }
  else
    {
      std::string ebis_file;
      pp.get("ebis_file",ebis_file);
      pout() << " recreating  geometry from file " << ebis_file << endl;
      //define ebis anew from file input
#ifdef CH_USE_HDF5
      HDF5Handle handleIn(ebis_file, HDF5Handle::OPEN_RDONLY);
      a_ebisPtr->define(handleIn);
      handleIn.close();
#endif
    }
}

BaseIF* makeChamber(const Real& radius,
                    const Real& thick,
                    const Real& offset,
                    const Real& height)
{
  RealVect zero(D_DECL(0.0,0.0,0.0));
  RealVect xAxis(D_DECL(1.0,0.0,0.0));
  bool inside = true;

  Vector<BaseIF*> pieces;

  // Create a chamber
  TiltedCylinderIF chamberOut(radius + thick/2.0,xAxis,zero, inside);
  TiltedCylinderIF chamberIn (radius - thick/2.0,xAxis,zero,!inside);

  IntersectionIF infiniteChamber(chamberIn,chamberOut);

  pieces.push_back(&infiniteChamber);

  RealVect normal1(D_DECL(1.0,0.0,0.0));
  RealVect point1(D_DECL(offset,0.0,0.0));
  PlaneIF plane1(normal1,point1,inside);

  pieces.push_back(&plane1);

  RealVect normal2(D_DECL(-1.0,0.0,0.0));
  RealVect point2(D_DECL(offset+height,0.0,0.0));
  PlaneIF plane2(normal2,point2,inside);

  pieces.push_back(&plane2);

  IntersectionIF* chamber = new IntersectionIF(pieces);

  return chamber;
}

BaseIF* makePlate(const Real& height,
                  const Real& thick,
                  const Real& radius,
                  const int&  doHoles,
                  const Real& holeRadius,
                  const Real& holeSpace)
{
  RealVect zero(D_DECL(0.0,0.0,0.0));
  RealVect xAxis(D_DECL(1.0,0.0,0.0));
  bool inside = true;

  // Create the plate without holes
  Vector<BaseIF*> pieces;

  RealVect normal1(D_DECL(1.0,0.0,0.0));
  RealVect point1(D_DECL(height,0.0,0.0));
  PlaneIF plane1(normal1,point1,inside);

  pieces.push_back(&plane1);

  RealVect normal2(D_DECL(-1.0,0.0,0.0));
  RealVect point2(D_DECL(height+thick,0.0,0.0));
  PlaneIF plane2(normal2,point2,inside);

  pieces.push_back(&plane2);

  TiltedCylinderIF middle(radius,xAxis,zero,inside);

  pieces.push_back(&middle);

  IntersectionIF plate(pieces);

  // Make the drills
  Vector<BaseIF*> drillBits;

  // Compute how many drills are needed in each direciton - 2*num+1 -
  // conservatively
  int num = (int)((radius - holeRadius) / holeSpace + 1.0);

  if (doHoles != 0)
  {
    for (int i = -num; i <= num; i++)
    {
      for (int j = -num; j <= num; j++)
      {
        RealVect center(D_DECL(0.0,i*holeSpace,j*holeSpace));
        TiltedCylinderIF* drill = new TiltedCylinderIF(holeRadius,xAxis,center,inside);

        drillBits.push_back(drill);
      }
    }
  }

  UnionIF drills(drillBits);
  ComplementIF notDrills(drills,true);

  // Drill the plate
  IntersectionIF* holyPlate = new IntersectionIF(plate,notDrills);

  return holyPlate;
}

BaseIF* makeVanes(const int&      num,
                  const Real&     thick,
                  const RealVect& normal,
                  const Real&     innerRadius,
                  const Real&     outerRadius,
                  const Real&     offset,
                  const Real&     height)
{
  RealVect zeroVect(D_DECL(0.0,0.0,0.0));
  RealVect xAxis(D_DECL(1.0,0.0,0.0));

  Vector<BaseIF*> eachVane;

  for (int i = 0; i < num; i++)
  {
    Real angle = i*2.*M_PI/num;

    BaseIF* oneVane = makeVane(thick,normal,innerRadius,outerRadius,offset,height,angle);

    eachVane.push_back(oneVane);
  }

  UnionIF* allVanes = new UnionIF(eachVane);

  return allVanes;
}

BaseIF* makeVane(const Real&     thick,
                 const RealVect& normal,
                 const Real&     innerRadius,
                 const Real&     outerRadius,
                 const Real&     offset,
                 const Real&     height,
                 const Real&     angle)
{
  RealVect zeroVect(D_DECL(0.0,0.0,0.0));
  RealVect xAxis(D_DECL(1.0,0.0,0.0));
  bool inside = true;

  Vector<BaseIF*> vaneParts;

  Real cosTheta;
  Real sinTheta;
  cosTheta = cos(angle);
  sinTheta = sin(angle);

#if CH_SPACEDIM == 3
  // Each side of the vane (infinite)
  // rotate the normal around x-axis
  RealVect normal1(D_DECL(normal[0],cosTheta*normal[1]-sinTheta*normal[2],sinTheta*normal[1]+cosTheta*normal[2]));
  // rotate point on top of vane around x-axis
  RealVect point(D_DECL(offset+height/2.0,-thick/2.0,0.0));
  RealVect point1(D_DECL(point[0],cosTheta*point[1]-sinTheta*point[2],sinTheta*point[1]+cosTheta*point[2]));
  PlaneIF plane1(normal1,point1,inside);

  vaneParts.push_back(&plane1);

  RealVect normal2(-normal1);
  // rotate point on bottom (-point[2] of vane around x-axis
  RealVect point2(D_DECL(point[0],-cosTheta*point[1]-sinTheta*point[2],-sinTheta*point[1]+cosTheta*point[2]));
  PlaneIF plane2(normal2,point2,inside);

  vaneParts.push_back(&plane2);
#endif

  // Make sure we only get something to the right of the origin
  RealVect normal3(D_DECL(0.0,-sinTheta,cosTheta));
  RealVect point3(D_DECL(0.0,0.0,0.0));
  PlaneIF plane3(normal3,point3,inside);

  vaneParts.push_back(&plane3);

  // Cut off the top and bottom
  RealVect normal4(D_DECL(1.0,0.0,0.0));
  RealVect point4(D_DECL(offset,0.0,0.0));
  PlaneIF plane4(normal4,point4,inside);

  vaneParts.push_back(&plane4);

  RealVect normal5(D_DECL(-1.0,0.0,0.0));
  RealVect point5(D_DECL(offset+height,0.0,0.0));
  PlaneIF plane5(normal5,point5,inside);

  vaneParts.push_back(&plane5);

  // The outside of the inner cylinder
  TiltedCylinderIF inner(innerRadius,xAxis,zeroVect,!inside);

  vaneParts.push_back(&inner);

  // The inside of the outer cylinder
  TiltedCylinderIF outer(outerRadius,xAxis,zeroVect,inside);

  vaneParts.push_back(&outer);

  IntersectionIF* vane = new IntersectionIF(vaneParts);

  return vane;
}

// BaseIF* makeVanes(const int&      num,
//                   const Real&     thick,
//                   const RealVect& normal,
//                   const Real&     innerRadius,
//                   const Real&     outerRadius,
//                   const Real&     offset,
//                   const Real&     height)
// {
//   RealVect zero(D_DECL(0.0,0.0,0.0));
//   RealVect xAxis(D_DECL(1.0,0.0,0.0));

//   BaseIF* oneVane = makeVane(thick,normal,innerRadius,outerRadius,offset,height);

//   Vector<BaseIF*> eachVane;

//   for (int i = 0; i < num; i++)
//   {
//     Real angle = i * 2*M_PI / num;

//     TransformIF* currentVane = new TransformIF(*oneVane);
//     currentVane->rotate(angle,zero,xAxis);

//     eachVane.push_back(currentVane);
//   }

//   UnionIF* allVanes = new UnionIF(eachVane);

//   return allVanes;
// }

// BaseIF* makeVane(const Real&     thick,
//                  const RealVect& normal,
//                  const Real&     innerRadius,
//                  const Real&     outerRadius,
//                  const Real&     offset,
//                  const Real&     height)
// {
//   RealVect zero(D_DECL(0.0,0.0,0.0));
//   RealVect xAxis(D_DECL(1.0,0.0,0.0));
//   bool inside = true;

//   Vector<BaseIF*> vaneParts;

//   // Each side of the vane (infinite)
//   RealVect normal1(normal);
//   RealVect point1(D_DECL(offset+height/2.0,-thick/2.0,0.0));
//   PlaneIF plane1(normal1,point1,inside);

//   vaneParts.push_back(&plane1);

//   RealVect normal2(-normal);
//   RealVect point2(D_DECL(offset+height/2.0,thick/2.0,0.0));
//   PlaneIF plane2(normal2,point2,inside);

//   vaneParts.push_back(&plane2);

//   // Make sure we only get something to the right of the origin
//   RealVect normal3(D_DECL(0.0,0.0,1.0));
//   RealVect point3(D_DECL(0.0,0.0,0.0));
//   PlaneIF plane3(normal3,point3,inside);

//   vaneParts.push_back(&plane3);

//   // Cut off the top and bottom
//   RealVect normal4(D_DECL(1.0,0.0,0.0));
//   RealVect point4(D_DECL(offset,0.0,0.0));
//   PlaneIF plane4(normal4,point4,inside);

//   vaneParts.push_back(&plane4);

//   RealVect normal5(D_DECL(-1.0,0.0,0.0));
//   RealVect point5(D_DECL(offset+height,0.0,0.0));
//   PlaneIF plane5(normal5,point5,inside);

//   vaneParts.push_back(&plane5);

//   // The outside of the inner cylinder
//   TiltedCylinderIF inner(innerRadius,xAxis,zero,!inside);

//   vaneParts.push_back(&inner);

//   // The inside of the outer cylinder
//   TiltedCylinderIF outer(outerRadius,xAxis,zero,inside);

//   vaneParts.push_back(&outer);

//   IntersectionIF* vane = new IntersectionIF(vaneParts);

//   return vane;
// }

Real getDistance(const RealVect& a_point,
                 const RealVect& a_lineDirection,
                 const RealVect& a_linePoint)
{
  RealVect diff = a_linePoint;
  diff -= a_point;

  Real normDiff2 = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    normDiff2 += diff[idir] * diff[idir];
  }

  Real normDirection2 = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    normDirection2 += a_lineDirection[idir] * a_lineDirection[idir];
  }

  Real dotDiffDirection2 = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    dotDiffDirection2 += diff[idir] * a_lineDirection[idir];
  }
  dotDiffDirection2 = dotDiffDirection2 * dotDiffDirection2;

  return sqrt(normDiff2 - dotDiffDirection2 / normDirection2);
}
