#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PlaneIF.H"
#include "PolyGeom.H"

#include "UsingNamespace.H"

//global variables.  yeah, i know.
Box g_domain;
Real g_dx;
RealVect g_origin;

/******************/
/// define grids by splitting up domain
/******************/
int
makeLayout(DisjointBoxLayout& dbl1,
           const Box& domain);
/***************/
// define a ramp EBIS.
/***************/
int makeGeometry(Box& domain,
                 Real& dx,
                 RealVect& origin,
                 int& upDir,
                 int& indepVar,
                 Real& startPt,
                 Real& slope);
/***************/
//make the corresponding layout
/***************/
int makeEBISL(EBISLayout& a_ebisl,
              const DisjointBoxLayout& a_grids,
              const Box& a_domain,
              const int& nghost);
/***************/
/***************/
int checkEBISL(const EBISLayout& a_ebisl,
               const DisjointBoxLayout& a_grids,
               const Box& a_domain,
               const  Real& dx,
               const  RealVect& origin,
               const  int& upDir,
               const  int& indepVar,
               const  Real& startPt,
               const  Real& slope);
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
  //begin forever present scoping trick
  {
    const char* in_file = "ramp.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    int eekflag = 0;
    //define the geometry object.
    //need the new and delete because of
    //strong construction
    //make the gometry.  this makes the first Chombo_EBIS
    Box domain;
    Real dx;
    RealVect origin;
    int upDir;
    int indepVar;
    Real startPt;
    Real slope;
    //and defines it using a geometryservice
    eekflag =  makeGeometry(domain,  dx, origin,
                            upDir, indepVar, startPt, slope);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeGeometry");
      }

    //make grids
    DisjointBoxLayout grids;
    eekflag = makeLayout(grids, domain);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeLayouts");
      }

    ///create ebislayout
    int nghost = 2;
    EBISLayout ebisl;
    eekflag = makeEBISL(ebisl, grids, domain, nghost);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeEBISL");
      }

    //check everything i can think of on finest level
    eekflag = checkEBISL(ebisl, grids, domain, dx, origin,
                         upDir, indepVar, startPt, slope);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in checkEBISL");
      }

    pout() << "ramp test passed" << endl;
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
                 Real& a_dx,
                 RealVect& a_origin,
                 int& a_upDir,
                 int& a_indepVar,
                 Real& a_startPt,
                 Real& a_slope)
{
  int eekflag =  0;
  //parse input file
  ParmParse pp;

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
      a_origin[idir] = prob_lo[idir];
    }
  pp.get("up_dir",a_upDir);
  a_upDir = SpaceDim-1;
  pp.get("indep_var",a_indepVar);
  pp.get("start_pt", a_startPt);
  pp.get("ramp_slope", a_slope);

  RealVect normal = RealVect::Zero;
  normal[a_upDir] = 1.0;
  normal[a_indepVar] = -a_slope;

  RealVect point = RealVect::Zero;
  point[a_upDir] = -a_slope*a_startPt;

  bool normalInside = true;

  PlaneIF ramp(normal,point,normalInside);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(ramp,0,vectDx);
  //this generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain, a_origin, a_dx, workshop);

  return eekflag;
}
/***************/
/***************/
int checkEBISL(const EBISLayout& a_ebisl,
               const DisjointBoxLayout& a_grids,
               const Box& a_domain,
               const Real& a_dx,
               const RealVect& a_origin,
               const int& a_upDir,
               const int& a_indepVar,
               const Real& a_startPt,
               const Real& a_slope)
{
  //check to see that the sum of the volume fractions
  //comes to close to exactly the total volume
  int eekflag =  0;
#ifdef CH_USE_FLOAT
  Real tolerance  = 60 * PolyGeom::getTolerance();
#else
  Real tolerance  = PolyGeom::getTolerance();
#endif
  const IntVect& ivSize= a_domain.size();
  RealVect hiCorn;
  RealVect domLen;
  Real cellVolume = 1.0;
  Real totalVolume = 1.0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      hiCorn[idir] = a_origin[idir] + a_dx*Real(ivSize[idir]);
      cellVolume *= a_dx;
      domLen[idir] = hiCorn[idir] - a_origin[idir];
      totalVolume *= domLen[idir];
    }
  //find normal and alpha in the space where the length of
  //the domain in all directions is unity
  int iy = a_upDir;
  int ix = a_indepVar;
  RealVect normal = RealVect::Zero;
  normal[iy] =          domLen[iy];
  normal[ix] = -a_slope*domLen[ix];
  Real alpha = a_slope*(a_origin[ix]-a_startPt)-a_origin[iy];
  Real volExact = totalVolume*PolyGeom::computeVolume(alpha, normal);

  RealVect correctNorm = RealVect::Zero;
  Real sumSquare;
  correctNorm[a_upDir] = 1.0;
  correctNorm[a_indepVar] = -a_slope;
  PolyGeom::unifyVector(correctNorm, sumSquare);
  //
  //voltot = sum(cellVolume*volFrac) = numerical domain volume.
  Real voltot = 0;
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = a_ebisl[dit()];
      for (BoxIterator bit(a_grids.get(dit())); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          Vector<VolIndex> vofs = ebisBox.getVoFs(iv);
          if (vofs.size() > 1)
            {
              eekflag = 2;
              return eekflag;
            }
          for (int ivof = 0; ivof < vofs.size(); ivof++)
            {
              Real volFrac = ebisBox.volFrac(vofs[ivof]);
              voltot += volFrac*cellVolume;
              if (volFrac > tolerance && volFrac < 1.0-tolerance)
                {
                  if (!ebisBox.isIrregular(iv))
                    {
                      eekflag = 4;
                      return eekflag;
                    }
                  RealVect normal= ebisBox.normal(vofs[ivof]);
                  for (int idir = 0; idir < SpaceDim; idir++)
                    {
                      if (Abs(normal[idir] - correctNorm[idir]) > tolerance)
                        {
                          eekflag = 5;
                          return eekflag;
                        }
                    }
                  RealVect centroid= ebisBox.centroid(vofs[ivof]);
                  for (int idir = 0; idir < SpaceDim; idir++)
                    {
                      if (Abs(centroid[idir]) > 0.5)
                        {
                          eekflag = 6;
                          return eekflag;
                        }
                    }
                }
            }
        }
    }
#ifdef CH_MPI
  Real t=voltot;
  MPI_Allreduce ( &t, &voltot, 1,
                  MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
#endif

  if (Abs(voltot -volExact) > tolerance)
    {
      pout() << Abs(voltot - volExact)
             << " > " << tolerance << endl;
      eekflag = 3;
      return eekflag;
    }
  return eekflag;
}

/***************/
/***************/
