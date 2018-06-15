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

#include "Box.H"
#include "ParmParse.H"
#include "EBIndexSpace.H"
#include "Vector.H"
#include "AllRegularService.H"
#include "PlaneIF.H"
#include "SlabService.H"
#include "GeometryShop.H"
#include "BaseIF.H"
#include "SphereIF.H"
#include "UnionIF.H"
#include "IntersectionIF.H"
#include "ComplementIF.H"
#include "TiltedCylinderIF.H"
#include "WedgeIF.H"
#include "EllipsoidIF.H"
#include "SphereArrayIF.H"
#include "TransformIF.H"
#include "PolynomialIF.H"
#include "LatheIF.H"
// #include "EBMenagerieUtils.H"
#include "BRMeshRefine.H"
#include "parstream.H"
#include "PolyGeom.H"
#include "ModianoIBCFactory.H"
#include "EBAMRGodunovFactory.H"
#include "EBPatchPolytropicFactory.H"
#include "EBAMRGodunovFactory.H"
#include "DEMIF.H"
#include "EBViscousTensorOpFactory.H"
#include "DirichletPoissonEBBC.H"
#include   "NeumannPoissonEBBC.H"
#include "DirichletPoissonDomainBC.H"
#include   "NeumannPoissonDomainBC.H"
#include "DirichletViscousTensorEBBC.H"
#include   "NeumannViscousTensorEBBC.H"
#include "DirichletViscousTensorDomainBC.H"
#include   "NeumannViscousTensorDomainBC.H"
#include "DirichletConductivityDomainBC.H"
#include   "NeumannConductivityDomainBC.H"
#include "DirichletConductivityEBBC.H"
#include   "NeumannConductivityEBBC.H"
#include "GodunovGeom.H"

#include "NamespaceHeader.H"

using std::ifstream;
using std::ios;

extern "C"
{
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
  /************/
  void
  getNormalDotAxis(BaseIVFAB<Real>&  a_error,
                   const IntVectSet& a_ivsIrreg,
                   const EBISBox&    a_ebisBox,
                   const RealVect&   a_cylinderAxis,
                   const Real& a_dx)
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
  getVolFrac(BaseIVFAB<Real>&  a_error,
             const IntVectSet& a_ivsIrreg,
             const EBISBox&    a_ebisBox,
             const RealVect&   a_cylinderAxis,
             const Real& a_dx)
  {
    for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real volfrac = a_ebisBox.volFrac(vof);
        a_error(vof, 0) = volfrac;
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
        RealVect centroidPt = getCentroidPt(a_ebisBox, vof, a_dx);

        RealVect trueNorm, closestPt;
        RealVect cylinderCorner = RealVect::Zero;
        PolyGeom::pointToLine(closestPt, trueNorm, centroidPt,
                              cylinderCorner, a_cylinderAxis);

        RealVect normal = a_ebisBox.normal(vof);
        RealVect diff = normal;
        diff -= trueNorm;

        Real error = PolyGeom::dot(a_cylinderAxis, diff);
        a_error(vof, 0) = error ;
      }
  }
  /************/
  void
  getNormMinTrueNorm(BaseIVFAB<Real>&  a_error,
                     const IntVectSet& a_ivsIrreg,
                     const EBISBox&    a_ebisBox,
                     const RealVect&   a_cylinderAxis,
                     const Real&       a_dx)
  {
    //dot normal with truenormal.  right answer == 1 so we subtract 1
    for (VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        RealVect centroidPt = getCentroidPt(a_ebisBox, vof, a_dx);

        RealVect trueNorm, closestPt;
        RealVect cylinderCorner = RealVect::Zero;
        PolyGeom::pointToLine(closestPt, trueNorm, centroidPt,
                              cylinderCorner, a_cylinderAxis);

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

  /*****/
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
        RealVect centroidPt = getCentroidPt(a_ebisBox, vof, a_dx);

        RealVect trueNorm, closestPt;
        RealVect cylinderCorner = RealVect::Zero;
        PolyGeom::pointToLine(closestPt, trueNorm, centroidPt,
                              cylinderCorner, a_cylinderAxis);

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
        RealVect centroidPt = getCentroidPt(a_ebisBox, vof, a_dx);

        RealVect trueNorm, closestPt;
        RealVect cylinderCorner = RealVect::Zero;
        PolyGeom::pointToLine(closestPt, trueNorm, centroidPt,
                              cylinderCorner, a_cylinderAxis);
        RealVect crossVec = PolyGeom::cross(trueNorm, a_cylinderAxis);

        RealVect normal = a_ebisBox.normal(vof);
        Real error = PolyGeom::dot(crossVec, normal);
        a_error(vof, 0) = error;
      }
  }
}
void
coarsenBoxes(Vector< Vector<Box>      >&    a_boxesCoar,
             const Vector<Vector<Box> >&    a_boxesFine,
             int a_refToCoar)
{
  a_boxesCoar.resize(a_boxesFine.size());
  for (int ilev = 0; ilev < a_boxesFine.size(); ilev++)
    {
      a_boxesCoar[ilev].resize(a_boxesFine[ilev].size());
      for (int ibox = 0; ibox < a_boxesFine[ilev].size(); ibox++)
        {
          a_boxesCoar[ilev][ibox] = coarsen(a_boxesFine[ilev][ibox], a_refToCoar);
        }
    }
}
void
getBoxes(Vector<Vector<Box> >&   a_boxes,
         Vector<int>&            a_refRat,
         const Box&              a_domain)
{
  ParmParse pp;
  int maxlev;
  pp.get("max_level", maxlev);
  int max_grid_size;
  pp.get("max_grid_size", max_grid_size);
  if (maxlev == 1)
    {
      int amrRef = 2;
      a_refRat.resize(2, amrRef);
      a_boxes.resize(2);
      Box fineBox = refine(a_domain, amrRef);
      int ishrink = fineBox.size(0);
      //this leaves 1/4 refined.
      ishrink *= 3;
      ishrink /= 8;
      fineBox.grow(-ishrink);
      Vector<Box> fineBoxes, coarBoxes;
      // check to see if coarsest level needs to be broken up
      domainSplit(a_domain, coarBoxes, max_grid_size, max_grid_size);
      domainSplit(fineBox , fineBoxes, max_grid_size, max_grid_size);

      a_boxes[0] = coarBoxes;
      a_boxes[1] = fineBoxes;
    }
  else if (maxlev == 0)
    {
      int amrRef = 2;
      a_refRat.resize(1, amrRef);
      a_boxes.resize(1);
      Vector<Box> coarBoxes;
      // check to see if coarsest level needs to be broken up
      domainSplit(a_domain, coarBoxes, max_grid_size, max_grid_size);
      a_boxes[0] = coarBoxes;
    }
  else
    {
      MayDay::Error("getboxes only written for maxlevel = 0 or maxlevel = 1 for now");
    }
}

/*****/
int getOrderEB()
{
  CH_TIME("PoissonUtilities::getOrderEB");
  ParmParse pp;
  int orderEB;
  pp.get("order_ebbc", orderEB);
  if (orderEB == 1)
    {
      pout() << "first order EB BC" << endl;
    }
  else if (orderEB == 2)
    {
      pout() << "second order EB BC" << endl;
    }
  else
    {
      MayDay::Error("bogus EBBC order (must be 1 or 2)" );
    }
  return orderEB;
}


void
getEBAMRGFactory(RefCountedPtr<EBAMRGodunovFactory>&  a_fact,
                 const  EBPatchPolytropicFactory* const    a_patch)
{
  ParmParse pp;
  Real gamma = 1.4;
  pp.get("gamma",gamma);

  int redistRad = 0;
  pp.get("redist_radius",redistRad);

  Real domainLength;
  pp.get("domain_length",domainLength);


  //dummies
  Real refineThresh = 0.7;
  if (pp.contains("refine_thresh"))
    {
      pp.get("refine_thresh", refineThresh);
    }
  int tagBufferSize = 1;
  int iusemassredist;
  pp.get("use_mass_redist", iusemassredist);
  bool useMassRedist    = (iusemassredist ==1);
  int verbosity;
  pp.get("verbosity", verbosity);

  Real initialCFL = 0.1;
  pp.get("initial_cfl",initialCFL);
  Real cfl = 0.8;
  pp.get("cfl",cfl);


  bool doRZCoords   = false;
  bool hasSourceTerm = false;
  bool doSmushing =true;
  a_fact = RefCountedPtr<EBAMRGodunovFactory>
    (new EBAMRGodunovFactory(initialCFL,
                             cfl,
                             redistRad,
                             domainLength*RealVect::Unit,
                             refineThresh,
                             tagBufferSize,
                             verbosity,
                             useMassRedist,
                             doSmushing,
                             doRZCoords,
                             hasSourceTerm,
                             a_patch));

}

void
getEBPPFactoryXY(RefCountedPtr<EBPatchPolytropicFactory>&  a_patchGamma,
                 const  EBPhysIBCFactory*  const           a_ibc)
{
  ParmParse pp;
  Real gamma = 1.4;
  pp.get("gamma",gamma);

  int uselim;
  pp.get("use_limiting", uselim);
  bool useLimiting = (uselim==1);

  Real specHeat = 1.0;
  pp.get("specific_heat", specHeat);

  int ifourth, iflatten, iartvisc;
  pp.get("use_fourth_order_slopes", ifourth);
  pp.get("use_flattening"         , iflatten);
  pp.get("use_art_visc"           , iartvisc);
  bool useFourthOrderSlopes = (ifourth  ==1);
  bool useFlattening        = (iflatten ==1);
  bool useArtificialVisc    = (iartvisc ==1);
  bool doRZCoords = false;

  EBPatchPolytropicFactory* newfact =
    (new EBPatchPolytropicFactory(a_ibc,
                                  gamma,
                                  specHeat,
                                  useFourthOrderSlopes,
                                  useFlattening,
                                  useArtificialVisc,
                                  useLimiting,
                                  doRZCoords));

  a_patchGamma = RefCountedPtr<EBPatchPolytropicFactory>(newfact);
}
/************/
void
getModianoIBCFactory(RefCountedPtr<ModianoIBCFactory>&  a_ibc)
{
  // read inputs
  ParmParse pp;

  int testverbosity;
  pp.get("testverbosity", testverbosity);

  Real gamma = 1.4;
  pp.get("gamma",gamma);


  if (testverbosity > 1)
    pout() << "Modiano initial and boundary conditions" << endl;

  RealVect modianoAxis;

  Real waveAmp, waveWidth;
  pp.get("wave_amplitude", waveAmp);
  pp.get("wave_width", waveWidth);
  int whichGeom;
  pp.get("which_geom", whichGeom);
  if (whichGeom == 3)
    {
      if (testverbosity > 1)
        pout() << "uses cylinder axis for direction" << endl;
      vector<Real>  cylinderAxisVect(SpaceDim);
      pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
      RealVect cylinderAxis;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          cylinderAxis[idir] = cylinderAxisVect[idir];
        }
      Real sum;
      PolyGeom::unifyVector(cylinderAxis, sum);
      modianoAxis = cylinderAxis;
    }
  else
    {
      if (testverbosity > 1)
        pout() << "uses channel slope for direction" << endl;
      RealVect channelNormal;
      vector<Real>  channelNormalVect(SpaceDim);
      pp.getarr("channel_normal",channelNormalVect, 0, SpaceDim);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          channelNormal[idir] = channelNormalVect[idir];
        }
      modianoAxis = mibcGetNormalVector(channelNormal);;
    }

  vector<Real> centerpp(SpaceDim,0.5);
  RealVect center;
  pp.getarr("wave_center",centerpp,0,SpaceDim);
  for (int i = 0; i < SpaceDim; i++)
    center[i] = centerpp[i];


  int idoNegWave;
  pp.get("use_negative_wave",idoNegWave);
  bool useNegativeWave = (idoNegWave==1);
  if (useNegativeWave)
    {
      if (testverbosity > 1)
        pout() << "u - c wave " << endl;
    }
  else
    {
      if (testverbosity > 1)
        pout() << "u + c wave " << endl;
    }

  int idoFreeStream;
  pp.get("free_stream_prob",idoFreeStream);
  bool doFreeStream = (idoFreeStream==1);
  if (doFreeStream)
    {
      if (testverbosity > 1)
        pout() << "free stream modiano bump " << endl;
    }
  else
    {
      if (testverbosity > 1)
        pout() << "simple wave prob " << endl;
    }


  a_ibc = RefCountedPtr<ModianoIBCFactory> (new ModianoIBCFactory(gamma, waveAmp, waveWidth,
                                                                  center, modianoAxis, doFreeStream, useNegativeWave));
}


/************/
void
godunovGeometry(Box& a_coarsestDomain,
                RealVect& a_dx)
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

  Real domain_length;//x dir
  pp.get("domain_length",domain_length);
  for (int idir = 0;idir<SpaceDim;idir++)
    {
      a_dx[idir] = domain_length/n_cell[0];//inputs have to be set such that domain lengths and n_cell give same dx for each coordinate dir
    }
  RealVect fineDx = a_dx;
  int ebMaxCoarsen = -1;
  for (int ilev = 0; ilev < max_level; ilev++)
    {
      fineDx /= refRatios[ilev];
    }
  int whichgeom;
  pp.get("which_geom",whichgeom);
  int ebMaxSize;
  ppgodunov.get("max_grid_size", ebMaxSize);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  int verbosity = 0;
  if (!pp.contains("ebis_file"))
    {
      if (whichgeom == 0)
        {
          //allregular
          pout() << "all regular geometry" << endl;

          AllRegularService regserv;
          ebisPtr->define(finestDomain, origin, fineDx[0], regserv, ebMaxSize, ebMaxCoarsen);
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

          GeometryShop workshop(ramp,verbosity,fineDx);
          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
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
          ebisPtr->define(finestDomain, origin, fineDx[0], slab, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 4)
        {

          RealVect cylinderAxis;
          vector<Real>  cylinderAxisVect(SpaceDim);
          pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              cylinderAxis[idir] = cylinderAxisVect[idir];
            }
          Real cylinderRadius;

          pp.get("cylinder_radius",cylinderRadius);
          pout() << "using a tilted cylinder implicit function" << endl;
          RealVect corner = RealVect::Zero;
          bool negativeInside = true;
          TiltedCylinderIF tunnel(cylinderRadius, cylinderAxis, corner, negativeInside);
          GeometryShop workshop(tunnel,verbosity,fineDx);
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
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
          SphereIF sphereIF(sphereRadius, sphereCenter, false);
          GeometryShop workshop(sphereIF,verbosity,fineDx);
          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 6)
        {
          pout() << "multisphere geometry" << endl;

          int numSpheres;
          pp.get("num_spheres", numSpheres);

          Vector<Real>     radius(numSpheres);
          Vector<RealVect> center(numSpheres);
          Vector<BaseIF*>  spheres(numSpheres);

          for (int isphere = 0; isphere < numSpheres; isphere++)
            {
              char radiusString[80];
              char centerString[80];
              sprintf(radiusString, "sphere_radius_%d", isphere);
              sprintf(centerString, "sphere_center_%d", isphere);

              Real sphereRadius;
              pp.get(radiusString, sphereRadius);

              vector<Real> sphere_center(SpaceDim);
              pp.getarr(centerString,sphere_center, 0, SpaceDim);

              RealVect sphereCenter;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  sphereCenter[idir] = sphere_center[idir];
                }

              center[isphere] = sphereCenter;
              radius[isphere] = sphereRadius;

              spheres[isphere] = new SphereIF(radius[isphere],
                                              center[isphere],
                                              false);
            }

          bool insideRegular = false;
          pp.query("inside",insideRegular);

          IntersectionIF impMultisphere(spheres);
          ComplementIF sideImpMultisphere(impMultisphere,insideRegular);

          GeometryShop workshop(sideImpMultisphere,verbosity,fineDx);

          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
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

            bool inside = (amp[iparabola] < 0);

            parabolas.push_back(new PolynomialIF(poly,inside));
          }

          IntersectionIF allTogether(parabolas);

          GeometryShop workshop(allTogether,verbosity,fineDx);
          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
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

          GeometryShop workshop(mirror,verbosity,fineDx);
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
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

          GeometryShop workshop(ovoid,verbosity,fineDx);
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 10)
        {

          RealVect wedgeCenter;;
          pout() << "wedge geometry" << endl;
          vector<Real>  wedgeCenterVect(SpaceDim);
          pp.getarr("wedge_center",wedgeCenterVect, 0, SpaceDim);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              wedgeCenter[idir] = wedgeCenterVect[idir];
            }
          Real wedgeSlope;

          pp.get("wedge_slope",wedgeSlope);
          int depvar;
          int indvar;
          pp.get("wedge_indvar", indvar);
          pp.get("wedge_depvar", depvar);
          WedgeIF wedge(wedgeCenter, wedgeSlope, depvar, indvar);
          GeometryShop workshop(wedge,verbosity,fineDx);
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
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
          GeometryShop workshop(rotazoid,verbosity,fineDx);
          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 12)
        {
          pout() << "DEMIF geometry" << endl;
          std::string demifFile;
          pp.get("demif_file",demifFile);
          IntVect demifNCell = finestDomain.size();
          int interpType;
          pp.get("demif_interp_type", interpType);
          Real bottomBuffer, highGround, verticalScale;
          pp.get("demif_bottom_buffer", bottomBuffer);
          pp.get("demif_high_ground", highGround);
          pp.get("demif_vertical_scale", verticalScale);


          DEMIF demifoid(demifNCell, interpType, fineDx, demifFile,
                         bottomBuffer, 1.e99,highGround, verticalScale);
          GeometryShop workshop(demifoid,verbosity,fineDx);
          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }

      else if (whichgeom == 16)
        {
          pout() << "Sphere Array geometry" << endl;
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
          int verbosity = 0;
          GeometryShop workshop(implicit,verbosity,fineDx);

          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 20)
        {
          pout() << "Gas jet nozzle geometry" << endl;

          BaseIF* retval;

          Vector<Real> prob_lo(SpaceDim,0.0);
          Vector<Real> prob_hi(SpaceDim,1.0);

          pp.getarr("prob_lo",prob_lo,0,SpaceDim);
          pp.getarr("prob_hi",prob_hi,0,SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              origin[idir] = prob_lo[idir];
            }

          bool insideRegular;
          int intInsideRegular;
          pp.get("insideRegular",intInsideRegular);

          if (intInsideRegular != 0) insideRegular = true;
          if (intInsideRegular == 0) insideRegular = false;

          // Data for polygons making up nozzle
          Vector<Vector<RealVect> > polygons;

          // Nothing initially
          polygons.resize(0);

          // For building each polygon
          Vector<RealVect> polygon;
          RealVect point(RealVect::Zero);

          // This is how far the polygons extend - the physical domain should be small
          // than this (20 meters)
          Real far = 20000.0;

          // Add body (two pieces) which union to one

          // Piece 1 - add the points and then save the list as a polygon
          polygon.resize(0);

          point[0] =       far; point[1] =  5.64012000;
          polygon.push_back(point);

          point[0] = 1.2499200; point[1] =  5.64012000;
          polygon.push_back(point);

          point[0] = 1.2499200; point[1] = -0.95678300;
          polygon.push_back(point);

          point[0] =       far; point[1] = -0.95678300;
          polygon.push_back(point);

          polygons.push_back(polygon);

          // Piece 2 - add the points and then save the list as a polygon
          polygon.resize(0);

          point[0] =  1.2499200; point[1] =  5.64012000;
          polygon.push_back(point);

          point[0] =  1.1058800; point[1] =  5.64012000;
          polygon.push_back(point);

          point[0] =  0.3966470; point[1] =  1.62330000;
          polygon.push_back(point);

          point[0] =  0.3964730; point[1] =  0.42332500;
          polygon.push_back(point);

          point[0] =  0.7493570; point[1] =  0.00176610;
          polygon.push_back(point);

          point[0] =  1.2499200; point[1] = -0.00777410;
          polygon.push_back(point);

          polygons.push_back(polygon);

          // Add poppet - add the points and then save the list as a polygon
          polygon.resize(0);

          point[0] =  0.8528400; point[1] = -far;
          polygon.push_back(point);

          point[0] =  0.8528400; point[1] = -0.92157200;
          polygon.push_back(point);

          point[0] =  0.8521970; point[1] = -0.32181000;
          polygon.push_back(point);

          point[0] =  0.0000000; point[1] =  0.69393192;
          polygon.push_back(point);

          point[0] =  0.0000000; point[1] = -far;
          polygon.push_back(point);

          polygons.push_back(polygon);

          // Make the vector of (convex) polygons (vectors of points) into a union
          // of convex polygons, each made from the intersection of a set of half
          // planes/spaces - all represented by implicit functions.
          UnionIF* crossSection = makeCrossSection(polygons);

          if (SpaceDim == 2)
            {
              // In 2D use "as is"

              // Complement if necessary
              ComplementIF insideOut(*crossSection,!insideRegular);
              retval = insideOut.newImplicitFunction();
            }
          else
            {
              // In 3D rotate about the z-axis and complement if necessary
              LatheIF implicit(*crossSection,insideRegular);
              retval = implicit.newImplicitFunction();
            }

          GeometryShop workshop(*retval,0,fineDx);

          // This generates the new EBIS
          EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
          ebisPtr->define(finestDomain,origin,fineDx[0],workshop,ebMaxSize,ebMaxCoarsen);

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
      ebisPtr->define(handleIn);
      handleIn.close();
#endif
    }
}
UnionIF* makeCrossSection(const Vector<Vector<RealVect> >& a_polygons)
{
  // The final result
  UnionIF* retval;

  // Get the number of polygons and make this inside of the domain
  // the inside of the polygons
  int numPolys = a_polygons.size();
  bool inside = true;

  // A list of all the polygons as implicit functions
  Vector<BaseIF*> polytopes;
  polytopes.resize(0);

  // Process each polygon
  for (int p = 0; p < numPolys; p++)
  {
    // All the half planes/spaces used to make a polygon
    Vector<BaseIF*> planes;
    planes.resize(0);

    // Get the current polygon (as a vector of points)
    const Vector<RealVect>& polygon = a_polygons[p];

    // Get the number of points in the polygon
    int numPts = polygon.size();

    // Process each pair of points
    for (int n = 0; n < numPts; n++)
    {
      // The normal and point is space used to specify each half plane/space
      RealVect normal(RealVect::Zero);
      RealVect point;

      // Set the normal remembering that the last point connects to the first
      // point.
      normal[0] = -(polygon[(n+1) % numPts][1] - polygon[n][1]);
      normal[1] =  (polygon[(n+1) % numPts][0] - polygon[n][0]);

      point = polygon[n];

      // Generate the appropriate half plane/space (as an implicit function)
      PlaneIF* plane;
      plane = new PlaneIF(normal,point,inside);

      // Save the result
      planes.push_back(plane);
    }

    // Intersect all the half planes/spaces to create an implicit function
    // that represents the polygon
    IntersectionIF* polygonIF = new IntersectionIF(planes);

    polytopes.push_back(polygonIF);
  }

  // Union all the polygon implicit functions to get the implicit function
  // returned
  retval = new UnionIF(polytopes);

  return retval;
}
/***************/
void godunovFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                       const ProblemDomain&  a_domain,
                       int                   a_max_level,
                       int                   a_max_grid_size,
                       int                   a_block_factor,
                       int                   a_verbosity,
                       std::string           a_grid_file)
{
  if (procID() == uniqueProc(SerialTask::compute))
    {
      a_amrGrids.push_back(Vector<Box>(1,a_domain.domainBox()));
      // read in predefined grids
     ifstream is(a_grid_file.c_str(), ios::in);

      if (is.fail())
        {
          MayDay::Error("Cannot open grids file");
        }

      // format of file -- number of levels, then for each level starting
      // with level 1, number of grids on level, list of boxes

      int in_numLevels;
      is >> in_numLevels;

      CH_assert (in_numLevels <= a_max_level+1);

      pout() << "numLevels = " << in_numLevels << endl;

      while (is.get() != '\n');

      a_amrGrids.resize(in_numLevels);

      // check to see if coarsest level needs to be broken up
      domainSplit(a_domain,a_amrGrids[0],a_max_grid_size,a_block_factor);

      if (a_verbosity >= 3)
        {
          pout() << "level 0: ";
          for (int n = 0; n < a_amrGrids[0].size(); n++)
            {
              pout() << a_amrGrids[0][0] << endl;
            }
        }

      // now loop over levels, starting with level 1
      int ngrid;
      for (int lev = 1; lev < in_numLevels; lev++)
        {
          is >> ngrid;

          if (a_verbosity >= 3)
            {
              pout() << "level " << lev << " numGrids = " << ngrid << endl;
              pout() << "Grids: ";
            }

          while (is.get() != '\n');

          a_amrGrids[lev].resize(ngrid);

          for (int i = 0; i < ngrid; i++)
            {
              Box bx;
              is >> bx;

              while (is.get() != '\n');

              // quick check on box size
              Box bxRef(bx);

              if (bxRef.longside() > a_max_grid_size)
                {
                  pout() << "Grid " << bx << " too large" << endl;
                  MayDay::Error();
                }

              if (a_verbosity >= 3)
              {
                pout() << bx << endl;
              }

              a_amrGrids[lev][i] = bx;
            } // end loop over boxes on this level
        } // end loop over levels
    }

  broadcast(a_amrGrids,uniqueProc(SerialTask::compute));
}

#include "NamespaceFooter.H"
