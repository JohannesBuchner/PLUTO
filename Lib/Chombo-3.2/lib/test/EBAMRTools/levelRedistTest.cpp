#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves oct 16 2001

#include "BaseIVFactory.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "AllRegularService.H"
#include "PlaneIF.H"
#include "EBLevelRedist.H"
#include "RedistStencil.H"
#include "SlabService.H"

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
int testConservation(const EBISLayout& a_ebisl,
                     const DisjointBoxLayout& a_grids,
                     const Box& a_domain,
                     const int& a_redistRad);

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

    const char* in_file = "levelredist.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    //make the gometry.  this makes the first Chombo_EBIS
    Box domain;
    Real dx;
    //and defines it using a geometryservice
    eekflag =  makeGeometry(domain,  dx);
    CH_assert(eekflag == 0);

    //make grids
    DisjointBoxLayout grids;
    eekflag = makeLayout(grids, domain);
    CH_assert(eekflag == 0);
    int redistRad;
    pp.get("redist_radius", redistRad);
    ///create ebislayout
    int nghost = 3*redistRad;
    EBISLayout ebisl;
    eekflag = makeEBISL(ebisl, grids, domain, nghost);
    if (eekflag != 0)
      {
        MayDay::Abort(" levelRedistTest found a problem in makeEBISL");

      }

    eekflag = testConservation(ebisl, grids,
                               domain, redistRad);
    if (eekflag != 0)
      {
        MayDay::Abort(" levelRedistTest found a conservation problem");

      }
  }//end scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  pout() << " levelRedistTest passed" << endl;
  return 0;
}
/***************/
/***************/
Real densityFunc(const IntVect& a_iv, const Box& a_box)
{
  Real retval;
  Real problen = a_box.size(0);
  Real x = a_iv[0];
  Real y = a_iv[1];
  x /= problen;
  y /= problen;

  retval = 1.0 + x + y + x*x + y*y;
  return retval;
}
/***************/
/***************/
Real massFunc(const IntVect& a_iv, const Box& a_box)
{
  Real retval;
  Real problen = a_box.size(0);
  Real x = a_iv[0];
  Real y = a_iv[1];
  x /= problen;
  y /= problen;
  retval = x*x*x + y*y*y;

  return retval;
}
/***************/
/***************/
int testConservation(const EBISLayout& a_ebisl,
                     const DisjointBoxLayout& a_grids,
                     const Box& a_domain,
                     const int& a_redistRad)
{
  int ncomp = 1;
  IntVect ivghost = IntVect::Zero;
  EBCellFactory ebcellfact(a_ebisl);
  LayoutData<IntVectSet> sets(a_grids);
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      sets[dit()] = a_ebisl[dit()].getIrregIVS(a_grids.get(dit()));
    }
  BaseIVFactory<Real> baseivfact(a_ebisl, sets);
  LevelData<EBCellFAB> solution(a_grids, ncomp, ivghost, ebcellfact);
  LevelData<BaseIVFAB<Real> > massdiff(a_grids, ncomp, ivghost, baseivfact);

  //initialize solution and mass difference
  //and add up the total mass in the system
  Real summassdiff = 0;
  Real sumsolmass = 0;
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = a_grids.get(dit());
      const EBISBox& ebisBox = a_ebisl[dit()];
      IntVectSet ivsIrreg = ebisBox.getIrregIVS(grid);
      BaseIVFAB<Real>& massFAB = massdiff[dit()];
      for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect&  iv = vof.gridIndex();
          Real mass = massFunc(iv, a_domain);
          massFAB(vof, 0) = mass;
          summassdiff += mass;
        }
      EBCellFAB& solFAB = solution[dit()];
      IntVectSet ivsBox(grid);
      for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect&  iv = vof.gridIndex();
          Real density = densityFunc(iv, a_domain);
          Real volFrac = ebisBox.volFrac(vof);
          solFAB(vof, 0) = density;
          sumsolmass += density*volFrac;
        }
    }

  Real localSumOld = sumsolmass + summassdiff ;

  Interval interv(0,0);

  //now redistribute massdiff into solution
  EBLevelRedist distributor(a_grids, a_ebisl, a_domain, ncomp, a_redistRad);
  distributor.setToZero();
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      distributor.increment(massdiff[dit()], dit(), interv);
    }
  distributor.redistribute(solution, interv);

  //now check that the solution has all the mass
  sumsolmass = 0;
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = a_grids.get(dit());
      const EBISBox& ebisBox = a_ebisl[dit()];
      EBCellFAB& solFAB = solution[dit()];
      IntVectSet ivsBox(grid);
      for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real density = solFAB(vof, 0);
          if (density > 1.0e-12)
            {
              Real volFrac = ebisBox.volFrac(vof);
              sumsolmass += density*volFrac;
            }
        }
    }
  Real localSumNew = sumsolmass;

  //gather what each processor thinks is the sum old and new
  Vector<Real> sumVecOld, sumVecNew;
  int baseproc = 0;
  gather(sumVecOld, localSumOld, baseproc);
  gather(sumVecNew, localSumNew, baseproc);
  Real massTotOld = 0;
  Real massTotNew = 0;
  if (procID() == baseproc)
    {
      CH_assert(sumVecOld.size() == numProc());
      CH_assert(sumVecNew.size() == numProc());
      for (int ivec = 0; ivec < numProc(); ivec++)
        {
          massTotOld += sumVecOld[ivec];
          massTotNew += sumVecNew[ivec];
        }
    }
  //broadcast the sum to all processors.
  broadcast(massTotOld, baseproc);
  broadcast(massTotNew, baseproc);


  int eekflag = 0;
  if (massTotOld > 1.0e-9)
    {
      Real relDiff = Abs(massTotOld - massTotNew)/massTotOld;
      if (relDiff > 1.5e-6)
        {
          pout() << "doh! " << endl;;
          pout() << "initial solution mass + diff      = "  << massTotOld << endl;
          pout() << "total mass in solution after dist = "  << massTotNew << endl;
          pout() << "relative difference               = "  << relDiff    << endl;
          eekflag  = 3;
        }
    }
  return eekflag;
}
/***************/
/***************/
int makeEBISL(EBISLayout& a_ebisl,
              const DisjointBoxLayout& a_grids,
              const Box& a_domain,
              const int& a_nghost)
{
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
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
#if (CH_SPACEDIM==2)
  Vector<int> n_cell(SpaceDim, 64);
#else
  Vector<int> n_cell(SpaceDim, 16);
#endif

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
      EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
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
      EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      ebisPtr->define(a_domain, origin, a_dx, workshop);
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
      EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      RealVect origin = RealVect::Zero;
      ebisPtr->define(a_domain, origin, a_dx, slab);
    }
  else
    {
      //bogus which_geom
      pout() << " bogus which_geom input = "
             << whichgeom << endl;
      eekflag = 33;
    }

  return eekflag;
}
/***************/
/***************/
