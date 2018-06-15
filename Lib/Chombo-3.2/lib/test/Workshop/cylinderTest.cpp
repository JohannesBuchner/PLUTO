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
#include <cstdio>
#include <iostream>

#include "ParmParse.H"
#include "EBArith.H"
#include "EBAMRIO.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "ParmParse.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "VoFIterator.H"
#include "TiltedCylinderIF.H"
#include  "RealVect.H"
#include "CH_Attach.H"
#include "UsingNamespace.H"

/***************/
RealVect getCentroidPt(const EBISBox&  a_ebisBox,
                       const VolIndex& a_vof,
                       const Real&     a_dx)
{
  RealVect centroid = a_ebisBox.bndryCentroid(a_vof);

  const IntVect&   iv = a_vof.gridIndex();
  RealVect centroidPt;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      centroidPt[idir] = Real(iv[idir]) + 0.5 + centroid[idir];
    }
  centroidPt *= a_dx;
  return centroidPt;
}

/***************/
RealVect
getTrueNorm(const EBISBox&  a_ebisBox,
            const VolIndex& a_vof,
            const Real&     a_dx)
{
  ParmParse pp;
  RealVect cylinderAxis;
  vector<Real>  cylinderAxisVect(SpaceDim);
  pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      cylinderAxis[idir] = cylinderAxisVect[idir];
    }
  Real sum;
  PolyGeom::unifyVector(cylinderAxis, sum);

  RealVect centroid = a_ebisBox.bndryCentroid(a_vof);
  //  const IntVect&   iv = a_vof.gridIndex();
  RealVect centroidPt = getCentroidPt(a_ebisBox, a_vof, a_dx);

  RealVect trueNorm, closestPt;
  RealVect cylinderCorner = RealVect::Zero;
  PolyGeom::pointToLine(closestPt, trueNorm, centroidPt,
                        cylinderCorner, cylinderAxis);
  PolyGeom::unifyVector(trueNorm, sum);
  return trueNorm;
}
/***************/
void getDebugIVS(IntVectSet& a_ivs,
                 const bool& a_thisCalcCoar)
{
  a_ivs = IntVectSet();
  if (!a_thisCalcCoar)
    {
      IntVect ivdebugloFine(D_DECL(66,10,10));
      IntVect ivdebughiFine(D_DECL(68,11,10));
      a_ivs |= ivdebugloFine;
      a_ivs |= ivdebughiFine;
    }
}

/***************/
void
getFinestDomain(Box&       a_domain,
                Real&      a_dx)
{
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
          MayDay::Error();
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_domain = Box(lo, hi);

  Real prob_hi;
  int numOpen = n_cell[0];
  pp.get("domain_length",prob_hi);
  a_dx = prob_hi/numOpen;
}
void
makeGeometry(const Box&       a_domain,
             const Real&      a_dx,
             Real&      a_cylinderRadius,
             RealVect&  a_cylinderAxis,
             int        a_cylType)
{
  //parse input file.  single level
  ParmParse pp;
  pp.get("cylinder_radius", a_cylinderRadius);
  CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int biggridsize;
  pp.get("max_grid_size", biggridsize);
  vector<Real>  cylinderAxisVect(SpaceDim);
  pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_cylinderAxis[idir] = cylinderAxisVect[idir];
    }
  Real sum;
  PolyGeom::unifyVector(a_cylinderAxis, sum);
  pout() << "using a tilted cylinder implicit function" << endl;
  RealVect corner = RealVect::Zero;
  //corner[0] = 0.0123456;
  //corner[1] = -0.00895423456;
  //corner[2] = .00234565432;
  bool negativeInside = true;
  TiltedCylinderIF tunnel(a_cylinderRadius, a_cylinderAxis, corner, negativeInside);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(tunnel,0,vectDx);
  int ebmaxcoarsen = 0;
  RealVect origin = RealVect::Zero;
  ebisPtr->define(a_domain, origin, a_dx, workshop, biggridsize, ebmaxcoarsen);
}
/************/
void
compareError(const BaseIVFAB<Real>& a_errorFine,
             const BaseIVFAB<Real>& a_errorCoar,
             const IntVectSet&      a_ivsIrregFine,
             const IntVectSet&      a_ivsIrregCoar,
             const EBISBox&         a_ebisBoxFine,
             const EBISBox&         a_ebisBoxCoar,
             const string&          a_varname)
{
  string varstring(a_varname);
  pout() << "===============================================" << endl;
  pout() <<  varstring << " at embedded boundary test " << endl;
  for (int inorm = 0; inorm <= 2; inorm++)
    {
      if (inorm == 0)
        {
          pout()  << "Using max norm." << endl;
        }
      else
        {
          pout()  << "Using L-" << inorm << " norm." << endl;
        }

      Real ebIrregNormCoar, ebIrregNormFine;

      int comp = 0;
      EBArith::irregNorm(ebIrregNormCoar,
                         a_errorCoar,
                         a_ivsIrregCoar,
                         a_ebisBoxCoar,
                         comp,  inorm);

      EBArith::irregNorm(ebIrregNormFine,
                         a_errorFine,
                         a_ivsIrregFine,
                         a_ebisBoxFine,
                         comp,  inorm);

      if (a_ivsIrregCoar.isEmpty())
        {
          pout() << "no irregular fine vofs" << endl;
        }
      else
        {
          pout() << varstring << " Error Norm Coar = " << ebIrregNormCoar << endl;
        }
      if (a_ivsIrregFine.isEmpty())
        {
          pout() << "no irregular fine vofs" << endl;
        }
      else
        {
          pout() <<  varstring << " Error Norm Fine = " << ebIrregNormFine << endl;
        }
      if ((Abs(ebIrregNormCoar) > 1.0e-12) && (Abs(ebIrregNormFine) > 1.0e-12))
        {
          Real order = log(ebIrregNormCoar/ebIrregNormFine)/log(2.0);
          pout() << "Order of " << varstring  <<" = " << order << endl << endl;
        }

    }
}
/************/
void
outputError(const BaseIVFAB<Real>& a_error,
            const IntVectSet&      a_ivsIrreg,
            const EBISBox&         a_ebisBox,
            bool                   a_isCoarse,
            int                    a_itype,
            const string&          a_varname)
{

  char filename[100];

  if (a_isCoarse)
    {
      sprintf(filename,"%s%dCoar.hdf5", a_varname.c_str(), a_itype);
    }
  else
    {
      sprintf(filename,"%s%dFine.hdf5", a_varname.c_str(), a_itype);
    }

  Box domain = a_ebisBox.getDomain().domainBox();
  EBCellFAB fabData(a_ebisBox, domain, 3);
  fabData.setVal(0.);
  for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const Real error = a_error(vof, 0);
      const Real bndryArea = a_ebisBox.bndryArea(vof);
      fabData(vof, 0) = error;
      fabData(vof, 1) = error*bndryArea;
      fabData(vof, 2) = bndryArea;
    }

  Vector<string> names(3);
  char namewithba[80];
  sprintf(namewithba, "%sXBndryArea", a_varname.c_str());

  names[0] = string(a_varname);
  names[1] = string(namewithba);
  names[2] = string("BoundaryArea");
#ifdef CH_USE_HDF5
#ifndef CH_MPI
  Real dx=1.0;
  Real dt=1.0;
  Real time=1.0;
  Vector<Real> covval(3, 0.0);
  bool repcov = false;

  writeEBHDF5(filename, domain, fabData, names, domain, dx, dt, time, repcov, covval);
#endif
#endif
}
/************/
void
getNormalDotAxis(BaseIVFAB<Real>&  a_error,
                 const IntVectSet& a_ivsIrreg,
                 const EBISBox&    a_ebisBox,
                 const RealVect&   a_cylinderAxis)
{
  //dot normal with the cylinder axis.  right answer == 0
  for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect normal = a_ebisBox.normal(vof);
      Real error = PolyGeom::dot(a_cylinderAxis, normal);
      a_error(vof, 0) = error;
    }
}
/************/
void
getNormMinTrueNormDotAxis(BaseIVFAB<Real>&  a_error,
                          const IntVectSet& a_ivsIrreg,
                          const EBISBox&    a_ebisBox,
                          const RealVect&   a_cylinderAxis,
                          const Real&       a_dx)
{
  //dot normal with truenormal.  right answer == 1 so we subtract 1
  for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect trueNorm  = getTrueNorm(a_ebisBox, vof, a_dx);

      RealVect normal = a_ebisBox.normal(vof);
      RealVect diff = normal;
      diff -= trueNorm;

      Real error = PolyGeom::dot(a_cylinderAxis, diff);
      a_error(vof, 0) = error ;
    }
}

/************/
void
getSumSqNormMinTrueNorm(BaseIVFAB<Real>&  a_error,
                        const IntVectSet& a_ivsIrreg,
                        const EBISBox&    a_ebisBox,
                        const RealVect&   a_cylinderAxis,
                        const Real&       a_dx)
{
  //dot normal with truenormal.  right answer == 1 so we subtract 1
  for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect trueNorm = getTrueNorm(a_ebisBox, vof, a_dx);

      RealVect normal = a_ebisBox.normal(vof);
      Real error = 0.0;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real diff = normal[idir] - trueNorm[idir];
          error += diff*diff;
        }
      a_error(vof, 0) = sqrt(error);
    }
}

/************/
void
getNormMinuTrueNorm(BaseIVFAB<Real>&  a_error,
                    const IntVectSet& a_ivsIrreg,
                    const EBISBox&    a_ebisBox,
                    const RealVect&   a_cylinderAxis,
                    const Real&       a_dx,
                    const int&        a_ivar)
{
  //dot normal with truenormal.  right answer == 1 so we subtract 1
  for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect trueNorm = getTrueNorm(a_ebisBox, vof, a_dx);

      RealVect normal = a_ebisBox.normal(vof);
      a_error(vof, 0) = normal[a_ivar] - trueNorm[a_ivar];
    }
}
void
getNormalDotTrueNormM1(BaseIVFAB<Real>&  a_error,
                       const IntVectSet& a_ivsIrreg,
                       const EBISBox&    a_ebisBox,
                       const RealVect&   a_cylinderAxis,
                       const Real&       a_dx)
{
  //dot normal with truenormal.  right answer == 1 so we subtract 1
  for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect trueNorm = getTrueNorm(a_ebisBox, vof, a_dx);

      RealVect normal = a_ebisBox.normal(vof);
      Real dotProd = PolyGeom::dot(trueNorm, normal);
      Real error = Abs(dotProd) - 1.0;

      a_error(vof, 0) = error ;
    }
}
/************/
void
getNormalDotCrossVec(BaseIVFAB<Real>&  a_error,
                     const IntVectSet& a_ivsIrreg,
                     const EBISBox&    a_ebisBox,
                     const RealVect&   a_cylinderAxis,
                     const Real&       a_dx)
{
  //dot normal with axis cross truenormal.  right answer == 0
  for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect trueNorm = getTrueNorm(a_ebisBox, vof, a_dx);
      RealVect crossVec = PolyGeom::cross(trueNorm, a_cylinderAxis);

      RealVect normal = a_ebisBox.normal(vof);

      Real error = PolyGeom::dot(crossVec, normal);
      a_error(vof, 0) = error;
    }
}

/************/
void
getCentroidDistError(BaseIVFAB<Real>&  a_error,
                     const IntVectSet& a_ivsIrreg,
                     const EBISBox&    a_ebisBox,
                     const RealVect&   a_cylinderAxis,
                     const Real&       a_dx,
                     const Real&       a_cylinderRadius)
{
  //dot normal with axis cross truenormal.  right answer == 0
  for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect centroidPt = getCentroidPt(a_ebisBox, vof, a_dx);

      RealVect trueNorm, closestPt;
      RealVect cylinderCorner = RealVect::Zero;
      PolyGeom::pointToLine(closestPt, trueNorm, centroidPt,
                            cylinderCorner, a_cylinderAxis);
      RealVect ebRadVec = closestPt;
      ebRadVec -= centroidPt;

      Real ebRadius = PolyGeom::dot(ebRadVec, ebRadVec);
      ebRadius = sqrt(ebRadius);
      Real error = ebRadius - a_cylinderRadius;
      a_error(vof, 0) = error;
    }
}
/************/
void
printCoveredSelect(const BaseIVFAB<Real>& a_error,
                   const EBISBox&         a_ebisBox,
                   const IntVectSet&      a_ivsDebug,
                   int ncomp,
                   const string&          a_name)
{
  pout() << setw(12)
         << setprecision(6)
         << setiosflags(ios::showpoint)
         << setiosflags(ios::scientific) ;
  IntVectSet ivs = a_ivsDebug;
  ivs &= a_error.getIVS();

  if (!ivs.isEmpty())
    {
      pout() << " printing debug " << a_name << endl;
      for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          pout() << vof << "   " ;
          for (int icomp = 0; icomp < ncomp; icomp++)
            {
              Real error = Abs(a_error(vof, icomp));
              pout() << error << "  ";
            }
          pout() << endl;
        }
    }
}
/************/
void
cylinderTest()
{
  ParmParse pp;
  //make layouts == domain
  Box domainBoxFine, domainBoxCoar;
  int verbosity, printdebug;
  pp.get("verbosity", verbosity);
  pp.get("print_debug", printdebug);
  Real dxFine, dxCoar;
  getFinestDomain(domainBoxFine, dxFine);
  dxCoar = 2.0*dxFine;
  domainBoxCoar = coarsen(domainBoxFine, 2);
  Vector<Box> boxFine(1, domainBoxFine);
  Vector<Box> boxCoar(1, domainBoxCoar);
  Vector<int> proc(1, 0);
  DisjointBoxLayout dblFine(boxFine, proc);
  DisjointBoxLayout dblCoar;
  coarsen(dblCoar, dblFine, 2);

  pout() << "==============================================" << endl;
  RealVect  cylinderAxis;
  Real  cylinderRadius;
  int itype = 0;
  makeGeometry(domainBoxFine,  dxFine, cylinderRadius, cylinderAxis,
               itype);
  EBISLayout ebislFine, ebislCoar;
  const CH_XD::EBIndexSpace* const ebisPtrFine = Chombo_EBIS::instance();
  ebisPtrFine->fillEBISLayout(ebislFine, dblFine, domainBoxFine, 0);

  const CH_XD::EBIndexSpace* const ebisPtrCoar = Chombo_EBIS::instance();
  makeGeometry(domainBoxCoar,  dxCoar, cylinderRadius, cylinderAxis,
               itype);
  ebisPtrCoar->fillEBISLayout(ebislCoar, dblCoar, domainBoxCoar, 0);

  int ifileout;
  pp.get("file_output", ifileout);
  // debug vofs== vofs to dump  errors of
  IntVectSet ivsDebugFine, ivsDebugCoar;
  getDebugIVS(ivsDebugFine, false);
  getDebugIVS(ivsDebugCoar, true);

  //do the whole convergence test thing.
  for (DataIterator dit = dblFine.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBoxFine = ebislFine[dit()];
      const EBISBox& ebisBoxCoar = ebislCoar[dit()];
      IntVectSet ivsIrregFine = ebisBoxFine.getIrregIVS(domainBoxFine);
      IntVectSet ivsIrregCoar = ebisBoxCoar.getIrregIVS(domainBoxCoar);

      BaseIVFAB<Real> errorFine(ivsIrregFine, ebisBoxFine.getEBGraph(), 1);
      BaseIVFAB<Real> errorCoar(ivsIrregCoar, ebisBoxCoar.getEBGraph(), 1);

      getNormalDotAxis(errorFine, ivsIrregFine, ebisBoxFine, cylinderAxis);
      getNormalDotAxis(errorCoar, ivsIrregCoar, ebisBoxCoar, cylinderAxis);
      string longname, shortname;

      longname = string("Normal dotted with  cylinder Axis");
      shortname = string("cylNormDotAxis");
      if (verbosity > 0)
        {
          compareError(errorFine, errorCoar,
                       ivsIrregFine, ivsIrregCoar,
                       ebisBoxFine, ebisBoxCoar, longname);
        }

      int nvar = 1;
      if (printdebug > 0)
        {
          printCoveredSelect(errorFine, ebisBoxFine, ivsDebugFine, nvar, longname);
          printCoveredSelect(errorCoar, ebisBoxCoar, ivsDebugCoar, nvar, longname);
        }

      if (ifileout == 1)
        {

          outputError(errorFine, ivsIrregFine, ebisBoxFine, false,itype, shortname);
          outputError(errorCoar, ivsIrregCoar, ebisBoxCoar, true, itype, shortname);
        }

      getNormalDotTrueNormM1(errorFine, ivsIrregFine, ebisBoxFine, cylinderAxis, dxFine);
      getNormalDotTrueNormM1(errorCoar, ivsIrregCoar, ebisBoxCoar, cylinderAxis, dxCoar);

      longname = string("norm dotted with true norm minus 1");
      shortname = string("cylNormDotTrueNormM1");
      if (printdebug > 0)
        {
          printCoveredSelect(errorFine, ebisBoxFine, ivsDebugFine, nvar,longname);
          printCoveredSelect(errorCoar, ebisBoxCoar, ivsDebugCoar, nvar,longname);
        }

      if (verbosity > 0)
        {
          compareError(errorFine, errorCoar,
                       ivsIrregFine, ivsIrregCoar,
                       ebisBoxFine, ebisBoxCoar, longname);
        }

      if (ifileout == 1)
        {
          outputError(errorFine, ivsIrregFine, ebisBoxFine, false,itype, shortname);
          outputError(errorCoar, ivsIrregCoar, ebisBoxCoar, true, itype, shortname);
        }

      getNormMinTrueNormDotAxis(errorFine, ivsIrregFine, ebisBoxFine, cylinderAxis, dxFine);
      getNormMinTrueNormDotAxis(errorCoar, ivsIrregCoar, ebisBoxCoar, cylinderAxis, dxCoar);

      longname = string("norm  minus true norm dot axis");
      shortname = string("cylNormMinuTrueNormDotAxis");
      if (printdebug > 0)
        {
          printCoveredSelect(errorFine, ebisBoxFine, ivsDebugFine, nvar, longname);
          printCoveredSelect(errorCoar, ebisBoxCoar, ivsDebugCoar, nvar, longname);
        }

      if (verbosity > 0)
        {
          compareError(errorFine, errorCoar,
                       ivsIrregFine, ivsIrregCoar,
                       ebisBoxFine, ebisBoxCoar, longname);
        }

      if (ifileout == 1)
        {
          outputError(errorFine, ivsIrregFine, ebisBoxFine, false,itype, shortname);
          outputError(errorCoar, ivsIrregCoar, ebisBoxCoar, true, itype, shortname);
        }

      getSumSqNormMinTrueNorm(errorFine, ivsIrregFine, ebisBoxFine, cylinderAxis, dxFine);
      getSumSqNormMinTrueNorm(errorCoar, ivsIrregCoar, ebisBoxCoar, cylinderAxis, dxCoar);
      longname = string("sqrt  sum sq(norm minu true norm)");
      shortname = string("cylSumSqNormMinTrueNorm");

      if (printdebug > 0)
        {
          printCoveredSelect(errorFine, ebisBoxFine, ivsDebugFine, nvar, longname);
          printCoveredSelect(errorCoar, ebisBoxCoar, ivsDebugCoar, nvar, longname);
        }

      if (verbosity > 0)
        {
          compareError(errorFine, errorCoar,
                       ivsIrregFine, ivsIrregCoar,
                       ebisBoxFine, ebisBoxCoar, longname);
        }

      if (ifileout == 1)
        {
          outputError(errorFine, ivsIrregFine, ebisBoxFine, false,itype, shortname);
          outputError(errorCoar, ivsIrregCoar, ebisBoxCoar, true, itype, shortname);
        }

#if CH_SPACEDIM==3
      getNormalDotCrossVec(errorFine, ivsIrregFine, ebisBoxFine, cylinderAxis, dxFine);
      getNormalDotCrossVec(errorCoar, ivsIrregCoar, ebisBoxCoar, cylinderAxis, dxCoar);
      longname = string("normal dotted with cross vector");
      shortname = string("cylNormDotCrossVec");

      if (printdebug > 0)
        {
          printCoveredSelect(errorFine, ebisBoxFine, ivsDebugFine, nvar, longname);
          printCoveredSelect(errorCoar, ebisBoxCoar, ivsDebugCoar, nvar, longname);
        }

      if (verbosity > 0)
        {
          compareError(errorFine, errorCoar,
                       ivsIrregFine, ivsIrregCoar,
                       ebisBoxFine, ebisBoxCoar, longname);
        }

      if (ifileout == 1)
        {
          outputError(errorFine, ivsIrregFine, ebisBoxFine, false,itype, shortname);
          outputError(errorCoar, ivsIrregCoar, ebisBoxCoar, true, itype, shortname);
        }
#endif

      getCentroidDistError(errorFine, ivsIrregFine, ebisBoxFine, cylinderAxis, dxFine, cylinderRadius);
      getCentroidDistError(errorCoar, ivsIrregCoar, ebisBoxCoar, cylinderAxis, dxCoar, cylinderRadius);
      longname = string("printing centroid Dist. from Axis - rad");
      shortname = string("cylCentroidDistError");

      if (printdebug > 0)
        {
          printCoveredSelect(errorFine, ebisBoxFine, ivsDebugFine, nvar, longname);
          printCoveredSelect(errorCoar, ebisBoxCoar, ivsDebugCoar, nvar, longname);
        }

      if (verbosity > 0)
        {
          compareError(errorFine, errorCoar,
                       ivsIrregFine, ivsIrregCoar,
                       ebisBoxFine, ebisBoxCoar, longname);
        }

      if (ifileout == 1)
        {
          outputError(errorFine, ivsIrregFine, ebisBoxFine, false,itype, shortname);
          outputError(errorCoar, ivsIrregCoar, ebisBoxCoar, true, itype, shortname);
        }

      for (int ivar = 0; ivar < SpaceDim; ivar++)
        {
          getNormMinuTrueNorm(errorFine, ivsIrregFine, ebisBoxFine, cylinderAxis, dxFine,  ivar);
          getNormMinuTrueNorm(errorCoar, ivsIrregCoar, ebisBoxCoar, cylinderAxis, dxCoar,  ivar);
          char charstr[100];
          sprintf(charstr, "printing norm minu true norm for ivar = %d", ivar);
          longname = string(charstr);
          sprintf(charstr, "normMinuTrueNormVar%d", ivar);
          shortname = string(charstr);

          if (printdebug > 0)
            {
              printCoveredSelect(errorFine, ebisBoxFine, ivsDebugFine, nvar, longname);
              printCoveredSelect(errorCoar, ebisBoxCoar, ivsDebugCoar, nvar, longname);
            }

          if (verbosity > 0)
            {
              compareError(errorFine, errorCoar,
                           ivsIrregFine, ivsIrregCoar,
                           ebisBoxFine, ebisBoxCoar, longname);
            }

          if (ifileout == 1)
            {
              outputError(errorFine, ivsIrregFine, ebisBoxFine, false,ivar, shortname);
              outputError(errorCoar, ivsIrregCoar, ebisBoxCoar, true, ivar, shortname);
            }
        }

    }
  pout() << "==============================================" << endl;
}
/************/
/************/
int
main(int a_argc, char* a_argv[])
{
  int retval = 0;
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  {
    // setChomboMPIErrorHandler();
#endif

    // Check for an input file
    ParmParse pp(0, NULL, NULL, "cylinder.inputs");

    cylinderTest();

#ifdef CH_MPI
  }
  MPI_Finalize();
#endif

  return retval;
}
