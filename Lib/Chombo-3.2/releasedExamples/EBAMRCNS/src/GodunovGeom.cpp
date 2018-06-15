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
#include "EBAMRGodunovFactory.H"
#include "EBPatchPolytropicFactory.H"
#include "EBAMRGodunovFactory.H"
#include "EBAMRCNSFactory.H"
#include "EBAMRCNSParams.H"
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
#include   "EBPlanarShockSolverBC.H"
#include   "EBPlanarShockTemperatureBC.H"
#include "EBCNSVortexF_F.H"
#include "EBCNSVortexIBC.H"
#include "SchlichtingVelocityEBBC.H"
#include "SchlichtingVeloDomainBC.H"
#include "SchlichtingTempDomainBC.H"
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
    for(int idir = 0; idir < SpaceDim; idir++)
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
    for(VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
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
    for(VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
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
    for(VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
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
    for(VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        RealVect centroidPt = getCentroidPt(a_ebisBox, vof, a_dx);

        RealVect trueNorm, closestPt;
        RealVect cylinderCorner = RealVect::Zero;
        PolyGeom::pointToLine(closestPt, trueNorm, centroidPt,
                              cylinderCorner, a_cylinderAxis);

        RealVect normal = a_ebisBox.normal(vof);
        Real error = 0.0;
        for(int idir = 0; idir < SpaceDim; idir++)
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
    for(VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
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
    for(VoFIterator vofit(a_ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
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
  for(int ilev = 0; ilev < a_boxesFine.size(); ilev++)
    {
      a_boxesCoar[ilev].resize(a_boxesFine[ilev].size());
      for(int ibox = 0; ibox < a_boxesFine[ilev].size(); ibox++)
        {
          a_boxesCoar[ilev][ibox] = coarsen(a_boxesFine[ilev][ibox], a_refToCoar);
        }
    }
}
/************/
void
schlichtingGeometry(Box& a_coarsestDomain,
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
  for(int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if(n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          MayDay::Error();
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_coarsestDomain = Box(lo, hi);
  Box finestDomain = a_coarsestDomain;
  for(int ilev = 0; ilev < max_level; ilev++)
    {
      finestDomain.refine(refRatios[ilev]);
    }

  Real domain_length;//x dir
  pp.get("domain_length",domain_length);
  for(int idir = 0;idir<SpaceDim;idir++)
    {
      a_dx[idir] = domain_length/n_cell[0];
    }
  RealVect fineDx = a_dx;
  int ebMaxCoarsen = -1;
  for(int ilev = 0; ilev < max_level; ilev++)
    {
      fineDx /= refRatios[ilev];
    }
  int ebMaxSize;
  ppgodunov.get("max_grid_size", ebMaxSize);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();

  pout() << "square cylinder geometry" << endl;

  SchlichtingParams schlicht;
  ParseSchlichtingParams(schlicht);
  RealVect schlichtingAxis  = schlicht.m_axis;

  RealVect schlichtNormal[CH_SPACEDIM-1];
  PolyGeom::getTangentVectors(schlichtNormal, schlicht.m_axis);

  RealVect schlichtCorner= schlicht.m_corner;
  bool inside = true;
  
  PlaneIF ramp(schlichtNormal[0],schlichtCorner,inside);
  
  GeometryShop workshop(ramp,0,fineDx);

  ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
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

  if(maxlev == 1)
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
      int maxGrid = Min(a_domain.size(0), max_grid_size);
      domainSplit(a_domain, coarBoxes, max_grid_size, maxGrid);
      domainSplit(fineBox , fineBoxes, max_grid_size, maxGrid);
      
      a_boxes[0] = coarBoxes;
      a_boxes[1] = fineBoxes;
    }
  else if(maxlev == 0)
    {
      int amrRef = 2;
      a_refRat.resize(1, amrRef);
      a_boxes.resize(1);
      Vector<Box> coarBoxes;
      // check to see if coarsest level needs to be broken up
      int maxGrid = Min(a_domain.size(0), max_grid_size);
      domainSplit(a_domain, coarBoxes, max_grid_size, maxGrid);
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
  if(orderEB == 1)
    {
      pout() << "first order EB BC" << endl;
    }
  else if(orderEB == 2)
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
getEBAMRCNSFactory(      RefCountedPtr<EBAMRCNSFactory>&                  a_fact,
                         EBAMRCNSParams&                                  a_params,
                         const RefCountedPtr<EBPatchPolytropicFactory>&   a_patch,
                         int a_iprob,
                         Real a_reynoldsNumber
                         )
{

  ParmParse pp;

  pp.get("redist_radius",   a_params.m_redistRad);
  pp.get("domain_length",   a_params.m_domainLength);
  pp.get("refine_thresh",   a_params.m_refineThresh);
  pp.get("use_mass_redist", a_params.m_useMassRedist);
  pp.get("verbosity",       a_params.m_verbosity);
  pp.get("cfl",             a_params.m_cfl);
  pp.get("initial_cfl",     a_params.m_initialDtMultiplier);
  pp.get("do_smushing",     a_params.m_doSmushing);
  pp.get("verbosity",       a_params.m_verbosity);
  pp.get("use_air_coefficients", a_params.m_useAirCoefs);
  a_params.m_variableCoeff =  a_params.m_useAirCoefs;
  if(a_reynoldsNumber <= 0)
    {
      pp.get("mu_viscosity",    a_params.m_viscosityMu);
      pp.get("lambda_viscosity",a_params.m_viscosityLa);
    }
  else
    {
      Real radius, maxvel;
      pp.get("sphere_radius",radius);
      FORT_GETMAXVEL(CHF_REAL(maxvel));

      a_params.m_viscosityMu = (2*maxvel*radius)/a_reynoldsNumber;
      a_params.m_viscosityLa = -2.0*a_params.m_viscosityMu/3.0;
    }
  pp.get("specific_heat",   a_params.m_specHeatCv);
  pp.get("thermal_cond",    a_params.m_thermalCond);
  pp.query("backward_euler", a_params.m_backwardEuler);
  a_params.m_tagBufferSize = 1;

  fillSolverBCs(a_params, a_iprob);

  Vector<int> domainBC;
  pp.getarr("dom_bc_type", domainBC, 0, SpaceDim);
 
  // The convergence tests are set up to use IBCs for initial conditions, 
  // so we stick an adaptor in here to make sure that things still work 
  // the way they should. 
  RefCountedPtr<EBSpaceTimeFunction> a_ICs(new EBSpaceTimeFunctionIBCAdaptor());
  a_fact = RefCountedPtr<EBAMRCNSFactory> (new EBAMRCNSFactory(a_params, a_patch, a_ICs));

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
  if(pp.contains("refine_thresh"))
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
  if(useLimiting)
    {
      pout() << "limiting ON" << endl;
    }
  else
    {
      pout() << "limiting OFF" << endl;
    }
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
getVortexIBCFactory(RefCountedPtr<EBCNSVortexIBCFactory>&  a_ibc,
                    const RealVect& a_center, const Real& a_radius)
{
  // read inputs
  ParmParse pp;

  int testverbosity;
  pp.get("testverbosity", testverbosity);

  Real gamma = 1.4;
  pp.get("gamma",gamma);
  
  pout() << "vortex iniitial and boundary conditions" << endl;

  vector<Real> centerpp(SpaceDim,0.5);
  RealVect center;
  pp.getarr("sphere_center",centerpp,0,SpaceDim);
  for (int i = 0; i < SpaceDim; i++)
    center[i] = centerpp[i];

  Real rnot;
  pp.get("sphere_radius", rnot);

  Real mach;
  pp.get("mach_number", mach);
  EBCNSVortexIBCFactory* ptr =  (new EBCNSVortexIBCFactory(gamma, center, mach, rnot));
  a_ibc = RefCountedPtr<EBCNSVortexIBCFactory> (ptr);
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
  for(int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if(n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          MayDay::Error();
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_coarsestDomain = Box(lo, hi);
  Box finestDomain = a_coarsestDomain;
  for(int ilev = 0; ilev < max_level; ilev++)
    {
      finestDomain.refine(refRatios[ilev]);
    }

  Real domain_length;//x dir
  pp.get("domain_length",domain_length);
  for(int idir = 0;idir<SpaceDim;idir++)
    {
      a_dx[idir] = domain_length/n_cell[0];//inputs have to be set such that domain lengths and n_cell give same dx for each coordinate dir
    }
  RealVect fineDx = a_dx;
  int ebMaxCoarsen = -1;
  for(int ilev = 0; ilev < max_level; ilev++)
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
      if(whichgeom == 0)
        {
          //allregular
          pout() << "all regular geometry" << endl;

          AllRegularService regserv;
          ebisPtr->define(finestDomain, origin, fineDx[0], regserv, ebMaxSize, ebMaxCoarsen);
        }
      else if(whichgeom == 1)
        {
          pout() << "ramp geometry" << endl;
          RealVect rampNormal;
          vector<Real>  rampNormalVect(SpaceDim);
          pp.getarr("ramp_normal",rampNormalVect, 0, SpaceDim);
          for(int idir = 0; idir < SpaceDim; idir++)
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
      else if(whichgeom == 2)
        {
          pout() << "slab geometry" << endl;
          vector<int> slab_lo(SpaceDim);
          pp.getarr("slab_lo",slab_lo,0,SpaceDim);
          vector<int> slab_hi(SpaceDim);
          pp.getarr("slab_hi",slab_hi,0,SpaceDim);
          IntVect lo, hi;
          for(int idir = 0; idir < SpaceDim; idir++)
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
      else if(whichgeom == 4)
        {

          RealVect cylinderAxis;
          RealVect corner = RealVect::Zero;
          vector<Real>  cylinderAxisVect(SpaceDim);
          vector<Real>  axisOriginVect(SpaceDim);
          pp.getarr("cylinder_axis",cylinderAxisVect, 0, SpaceDim);
          pp.getarr("axis_origin",axisOriginVect, 0, SpaceDim);
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              cylinderAxis[idir] = cylinderAxisVect[idir];
              corner[idir] = axisOriginVect[idir];
            }
          Real cylinderRadius;

          pp.get("cylinder_radius",cylinderRadius);
          pout() << "using a tilted cylinder implicit function" << endl;
          
          bool negativeInside = true;
          TiltedCylinderIF tunnel(cylinderRadius, cylinderAxis, corner, negativeInside);
          GeometryShop workshop(tunnel,verbosity,fineDx);
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if(whichgeom == 5)
        {
          pout() << "sphere geometry" << endl;
          vector<Real> sphere_center(SpaceDim);
          pp.getarr("sphere_center",sphere_center, 0, SpaceDim);
          RealVect sphereCenter;
          for(int idir = 0; idir < SpaceDim; idir++)
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

          for(int isphere = 0; isphere < numSpheres; isphere++)
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
              for(int idir = 0; idir < SpaceDim; idir++)
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
      else if(whichgeom == 7)
        {
          pout() << "multiparabola geometry" << endl;
          int numParabola;
          pp.get("num_parabolas", numParabola);
          int updir;
          pp.get("parabola_updir", updir);
          Vector<Real>     amp(numParabola);
          Vector<RealVect> center(numParabola);
          for(int iparabola = 0; iparabola < numParabola; iparabola++)
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
              for(int idir = 0; idir < SpaceDim; idir++)
                {
                  parabolaCenter[idir] = parabola_center[idir];
                }
              center[iparabola] = parabolaCenter;
              amp[iparabola]    = parabolaAmp;
            }

          Vector<BaseIF*> parabolas;
          for(int iparabola = 0; iparabola < numParabola; iparabola++)
            {
              Vector<PolyTerm> poly;

              PolyTerm mono;
              Real coef;
              IntVect powers;

              if (updir != 0) {
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
      else if(whichgeom == 8)
        {
          pout() << "parabolic mirror geometry" << endl;
          Real amplitude;
          RealVect center;
          vector<Real> centervec;
          int updir;
          pp.get("mirror_updir", updir);
          pp.get("mirror_amplitude", amplitude);
          pp.getarr("mirror_center",centervec, 0, SpaceDim);
          for(int idir = 0; idir < SpaceDim; idir++)
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
      else if(whichgeom == 9)
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
          for(int idir = 0; idir < SpaceDim; idir++)
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
      else if(whichgeom == 10)
        {

          RealVect wedgeCenter;;
          pout() << "wedge geometry" << endl;
          vector<Real>  wedgeCenterVect(SpaceDim);
          pp.getarr("wedge_center",wedgeCenterVect, 0, SpaceDim);
          for(int idir = 0; idir < SpaceDim; idir++)
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
      else if(whichgeom == 11)
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
          for(int idir = 0; idir < SpaceDim; idir++)
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
      else if(whichgeom == 12)
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

      else if(whichgeom == 16)
        {
          pout() << "Sphere Array geometry" << endl;
          Real sphereRadius;
          pp.get("sphere_radius", sphereRadius);

          vector<Real> first_sphere_center(SpaceDim);
          vector<Real> sphere_spacing(SpaceDim);
          pp.getarr("first_sphere_center",first_sphere_center, 0, SpaceDim);
          pp.getarr("sphere_spacing",sphere_spacing, 0, SpaceDim);
          RealVect firstCenter, sphereSpacing;
          for(int idir = 0; idir < SpaceDim; idir++)
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
      else if(whichgeom == 20)
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

          bool insideRegular = false;

          // Data for polygons making up nozzle
          Vector<Vector<RealVect> > polygons;

          
          // Add body (two pieces) which union to one
          int num_polygons_in_body, num_polygons_in_poppet;
          pp.get("num_polygons_in_body", num_polygons_in_body);
          pp.get("num_polygons_in_poppet", num_polygons_in_poppet);
          int points_in_body_polygon_0;
          int points_in_body_polygon_1;
          pp.get("points_in_body_polygon_0", points_in_body_polygon_0);
          pp.get("points_in_body_polygon_1", points_in_body_polygon_1);
          int nbodypts[2]; 
          nbodypts[0] = points_in_body_polygon_0; 
          nbodypts[1] = points_in_body_polygon_1;
          

          // Piece 1 - add the points and then save the list as a polygon

          for(int ipoly = 0; ipoly < num_polygons_in_body; ipoly++)
            {
              Vector<RealVect> polygon(0);
              for(int ipt = 0; ipt < nbodypts[ipoly]; ipt++)
                {
                  Vector<Real> rvpoints(SpaceDim);
                  char pointname[100];
                  sprintf(pointname, "bodypt%d%d", ipoly,ipt);
                  pp.getarr(pointname, rvpoints, 0, SpaceDim);
                  RealVect point(RealVect::Zero);
                  for(int idir =0; idir < SpaceDim; idir++)
                    {
                      point[idir] = rvpoints[idir];
                    }
                  polygon.push_back(point);
                }
              polygons.push_back(polygon);
            }

          // Add poppet - add the points and then save the list as a polygon

          int points_in_poppet_polygon;
          pp.get("points_in_poppet_polygon", points_in_poppet_polygon);
          for(int ipoly = 0; ipoly < num_polygons_in_poppet; ipoly++)
            {
              Vector<RealVect> polygon(0);
              for(int ipt = 0; ipt < points_in_poppet_polygon; ipt++)
                {
                  Vector<Real> rvpoints(SpaceDim);
                  char pointname[100];
                  sprintf(pointname, "poppetpt%d%d", ipoly,ipt);
                  pp.getarr(pointname, rvpoints, 0, SpaceDim);
                  RealVect point(RealVect::Zero);
                  for(int idir =0; idir < SpaceDim; idir++)
                    {
                      point[idir] = rvpoints[idir];
                    }
                  polygon.push_back(point);
                }
              polygons.push_back(polygon);
            }

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
      else if (whichgeom == 21)
        {
          pout() << "multicylinder geometry" << endl;

          int numCyl;
          pp.get("num_cylinders", numCyl);

          Vector<Real>     cylradius(numCyl);
          Vector<RealVect>   cylaxis(numCyl);
          Vector<RealVect> cylorigin(numCyl);
          Vector<BaseIF*>  cylinders(numCyl);

          for(int iCyl = 0; iCyl < numCyl; iCyl++)
            {
              char radiusString[80];
              char axisString[80];
              char originString[80];
              sprintf(radiusString, "cylinder_radius_%d", iCyl);
              sprintf(  axisString, "cylinder_axis_%d"  , iCyl);
              sprintf(  originString, "cylinder_origin_%d"  , iCyl);

              Real cylRadius;
              pp.get(radiusString, cylRadius);

              vector<Real> cyl_axis(SpaceDim);
              pp.getarr(axisString,cyl_axis, 0, SpaceDim);

              vector<Real> cyl_origin(SpaceDim);
              pp.getarr(originString,cyl_origin, 0, SpaceDim);

              RealVect cylAxis;
              RealVect cylOrigin;
              for(int idir = 0; idir < SpaceDim; idir++)
                {
                  cylOrigin[idir] = cyl_origin[idir];
                  cylAxis[  idir] = cyl_axis[  idir];
                }

              cylaxis[  iCyl] = cylAxis;
              cylorigin[iCyl] = cylOrigin;
              cylradius[iCyl] = cylRadius;

              cylinders[iCyl] =  new TiltedCylinderIF(cylRadius,cylAxis,cylOrigin,true);
            }

          bool outsideRegular = false;
          pp.query("outside",outsideRegular);

          UnionIF impMultisphere(cylinders);
          ComplementIF sideImpMultisphere(impMultisphere,outsideRegular);

          GeometryShop workshop(sideImpMultisphere,verbosity,fineDx);

          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }
      else if (whichgeom == 22)
        {
          pout() << "capped multicylinder geometry" << endl;

          int numCyl;
          pp.get("num_cylinders", numCyl);

          int capping;
          // capping == 0 -> flat caps
          // capping != 0 -> spherical caps
          pp.get("capping",capping);

          Vector<Real>     cylradius(numCyl);
          Vector<Real>     cyllength(numCyl);
          Vector<RealVect>   cylaxis(numCyl);
          Vector<RealVect> cylorigin(numCyl);
          Vector<BaseIF*>  cylinders(numCyl);

          for(int iCyl = 0; iCyl < numCyl; iCyl++)
            {
              char lengthString[80];
              char radiusString[80];
              char axisString[80];
              char originString[80];
              sprintf(lengthString, "cylinder_length_%d", iCyl);
              sprintf(radiusString, "cylinder_radius_%d", iCyl);
              sprintf(  axisString, "cylinder_axis_%d"  , iCyl);
              sprintf(  originString, "cylinder_origin_%d"  , iCyl);

              Real cylRadius;
              pp.get(radiusString, cylRadius);

              Real cylLength;
              pp.get(lengthString, cylLength);

              vector<Real> cyl_axis(SpaceDim);
              pp.getarr(axisString,cyl_axis, 0, SpaceDim);

              vector<Real> cyl_origin(SpaceDim);
              pp.getarr(originString,cyl_origin, 0, SpaceDim);

              RealVect cylAxis;
              RealVect cylOrigin;
              for(int idir = 0; idir < SpaceDim; idir++)
                {
                  cylOrigin[idir] = cyl_origin[idir];
                  cylAxis[  idir] = cyl_axis[  idir];
                }

              cylaxis[  iCyl] = cylAxis;
              cylorigin[iCyl] = cylOrigin;
              cylradius[iCyl] = cylRadius;
              // Start with an infinite cylinder whose axis is the x-axis
              RealVect zero(D_DECL(0.0,0.0,0.0));
              RealVect xAxis(D_DECL(1.0,0.0,0.0));
              bool inside = true;
              TiltedCylinderIF cylinder(cylRadius,xAxis,zero,inside);
              RealVect normal1(D_DECL(1.0,0.0,0.0));
              RealVect point1(D_DECL(-cylLength/2.0,0.0,0.0));
              PlaneIF  plane1(normal1,point1,inside);

              RealVect normal2(D_DECL(-1.0,0.0,0.0));
              RealVect point2(D_DECL(cylLength/2.0,0.0,0.0));
              PlaneIF  plane2(normal2,point2,inside);

              IntersectionIF interval(plane1,plane2);
              IntersectionIF finiteCylinder(cylinder,interval);

              TransformIF* transform;

              // Add spherical caps if needed
              if (capping == 0)
                {
                  transform = new TransformIF(finiteCylinder);
                }
              else
                {
                  SphereIF sphere1(cylRadius,point1,inside);
                  SphereIF sphere2(cylRadius,point2,inside);

                  UnionIF spheres(sphere1,sphere2);
                  UnionIF cappedCylinder(finiteCylinder,spheres);

                  transform = new TransformIF(cappedCylinder);
                }
              // Rotate the x-axis to the cylinder axis
              transform->rotate(xAxis, cylAxis);

              // Translate to the center
              transform->translate(cylOrigin);

              cylinders[iCyl] =  transform;
            }

          bool outsideRegular = false;
          pp.query("outside",outsideRegular);

          UnionIF impMultisphere(cylinders);
          ComplementIF sideImpMultisphere(impMultisphere,outsideRegular);

          GeometryShop workshop(sideImpMultisphere,verbosity,fineDx);

          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);
        }

      else if (whichgeom == 23)
        {
          pout() << "sphere and cylinder geometry" << endl;
          Real cylRadius, sphRadius, cylLength;
          RealVect cylOrigin, sphCenter;
          RealVect cylAxis;
          pp.get("cylinder_radius", cylRadius);
          pp.get("sphere_radius", sphRadius);
          pp.get("cylinder_length", cylLength);

          vector<Real> cyl_axis(SpaceDim);
          pp.getarr("cylinder_axis",cyl_axis, 0, SpaceDim);
          
          vector<Real> cyl_origin(SpaceDim);
          pp.getarr("cylinder_origin",cyl_origin, 0, SpaceDim);

          vector<Real> sph_center(SpaceDim);
          pp.getarr("sphere_center",sph_center, 0, SpaceDim);
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              cylOrigin[idir] = cyl_origin[idir];
              cylAxis[  idir] = cyl_axis[  idir];
              sphCenter[idir] = sph_center[idir];
            }
          
          bool inside = true;
          SphereIF* sphereif = new SphereIF(sphRadius,sphCenter,inside);

          // Start with an infinite cylinder whose axis is the x-axis
          RealVect zero(D_DECL(0.0,0.0,0.0));
          RealVect xAxis(D_DECL(1.0,0.0,0.0));

          TiltedCylinderIF cylinder(cylRadius,xAxis,zero,inside);
          RealVect normal1(D_DECL(1.0,0.0,0.0));
          RealVect point1(D_DECL(-cylLength/2.0,0.0,0.0));
          PlaneIF  plane1(normal1,point1,inside);

          RealVect normal2(D_DECL(-1.0,0.0,0.0));
          RealVect point2(D_DECL(cylLength/2.0,0.0,0.0));
          PlaneIF  plane2(normal2,point2,inside);

          IntersectionIF interval(plane1,plane2);
          IntersectionIF finiteCylinder(cylinder,interval);

          SphereIF sphere1(cylRadius,point1,inside);
          SphereIF sphere2(cylRadius,point2,inside);

          UnionIF spheres(sphere1,sphere2);
          UnionIF cappedCylinder(finiteCylinder,spheres);

          TransformIF*    transformif = new TransformIF(cappedCylinder);
          
          // Rotate the x-axis to the cylinder axis
          transformif->rotate(xAxis, cylAxis);

          // Translate to the center
          transformif->translate(cylOrigin);

          Vector<BaseIF*> cylinders(2);
          cylinders[0] = sphereif;
          cylinders[1] = transformif;
          UnionIF impMultisphere(cylinders);
          bool outsideRegular = false;
          pp.query("outside",outsideRegular);

          ComplementIF sideImpMultisphere(impMultisphere,outsideRegular);

          GeometryShop workshop(sideImpMultisphere,verbosity,fineDx);

          //this generates the new EBIS
          ebisPtr->define(finestDomain, origin, fineDx[0], workshop, ebMaxSize, ebMaxCoarsen);

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

              if (a_verbosity >= 3) {
                pout() << bx << endl;
              }

              a_amrGrids[lev][i] = bx;
            } // end loop over boxes on this level
        } // end loop over levels
    }

  broadcast(a_amrGrids,uniqueProc(SerialTask::compute));
}
/***************/
void
fillSolverBCs(EBAMRCNSParams& a_params, const int& a_iprob)
{
  ParmParse pp;
  pp.get("do_diffusion", a_params.m_doDiffusion);
  if(a_params.m_doDiffusion)
    {
      pout() << " Homogeneous Neumann boundary conditions for our smoothing operators " << endl;
      NeumannPoissonEBBCFactory* neumbcSmooth = new NeumannPoissonEBBCFactory();
      neumbcSmooth->setValue(0.);
      a_params.m_ebBCSmooth = RefCountedPtr<BaseEBBCFactory>(neumbcSmooth);
      NeumannPoissonDomainBCFactory* neumdobcSmooth = new NeumannPoissonDomainBCFactory();
      neumdobcSmooth->setValue(0.);
      a_params.m_doBCSmooth = RefCountedPtr<BaseDomainBCFactory>(neumdobcSmooth);

      pout() << " Thermally insulated embedded boundaries "  << endl;
      NeumannConductivityEBBCFactory* neumbctemp = new NeumannConductivityEBBCFactory();
      neumbctemp->setValue(0.);
      a_params.m_ebBCTemp = RefCountedPtr<BaseEBBCFactory>(neumbctemp);
      if(a_iprob == 2)
        {
          pout() << "Schlichting problem " << endl;
          SchlichtingParams schlicht;
          ParseSchlichtingParams(schlicht);

          a_params.m_ebBCVelo = RefCountedPtr<BaseEBBCFactory>(    new SchlichtingVelocityEBBCFactory(schlicht));
          a_params.m_doBCVelo = RefCountedPtr<BaseDomainBCFactory>(new SchlichtingVeloDomainBCFactory(schlicht));
          a_params.m_doBCTemp = RefCountedPtr<BaseDomainBCFactory>(new SchlichtingTempDomainBCFactory(schlicht));
        }
      else
        {
          pout() << "No slip, no flow velocity EBBC" << endl;
          DirichletViscousTensorEBBCFactory* diribc = new DirichletViscousTensorEBBCFactory();
          diribc->setValue(0.);
          a_params.m_ebBCVelo = RefCountedPtr<BaseEBBCFactory>(diribc);

          if(a_iprob == 1)
            {
              pout() << "planar shock domain bcs" << endl;
              int inormal;

              pp.get("shock_normal",inormal);

              int ishockback;
              pp.get("shock_backward",ishockback);
              bool shockback = (ishockback == 1);
              a_params.m_doBCVelo = RefCountedPtr<BaseDomainBCFactory>(new EBPlanarShockSolverBCFactory(inormal, shockback, a_params.m_slipBoundaries));

              a_params.m_doBCTemp = RefCountedPtr<BaseDomainBCFactory>(new EBPlanarShockTemperatureBCFactory(inormal, shockback, a_params.m_slipBoundaries));
            }
          else if(a_iprob == 0)
            {
              pout() << "no slip no flow, thermally insulated domain bcs" << endl;
              DirichletViscousTensorDomainBCFactory* diribc = new DirichletViscousTensorDomainBCFactory();
              diribc->setValue(0.);
              a_params.m_doBCVelo = RefCountedPtr<BaseDomainBCFactory>(diribc);

              NeumannConductivityDomainBCFactory* neumdobc = new NeumannConductivityDomainBCFactory();
              neumdobc->setValue(0.);
              a_params.m_doBCTemp = RefCountedPtr<BaseDomainBCFactory>(neumdobc);
            }
          else
            {
              MayDay::Error("bogus problem identifier");
            }
        }    
    }
}
/***************/
void
fillAMRParams(EBAMRCNSParams& a_params, int a_iprob)
{
  // read inputs
  ParmParse ppgodunov;
  a_params.m_slipBoundaries = false;
  a_params.m_explicitHyperbolicSource = false;
  
  
  ppgodunov.query("check_max_and_min",a_params.m_checkMaxMin);
  ppgodunov.query("explicit_hyperbolic_source",a_params.m_explicitHyperbolicSource);
  if(a_params.m_explicitHyperbolicSource)
    {
      pout() << "doing explicit hyperbolic source term" << endl;
    }
  else
    {
      pout() << "doing implicit hyperbolic source term" << endl;
    }
  
  ppgodunov.query("backward_euler",a_params.m_backwardEuler);
  if(a_params.m_backwardEuler)
    {
      pout() << "advancing parabolic terms using backward euler" << endl;
    }
  else
    {
      pout() << "advancing parabolic terms using TGA" << endl;
    }

  ppgodunov.query("slip_boundaries",a_params.m_slipBoundaries);
  ppgodunov.query("use_backward_euler",a_params.m_backwardEuler);
  
  if(a_params.m_slipBoundaries)
    {
      pout() << "slip walls at domain boundaries" << endl;
    }
  else
    {
      pout() << "no slip walls at domain boundaries" << endl;
    }
  ppgodunov.get("verbosity",a_params.m_verbosity);
  CH_assert(a_params.m_verbosity >= 0);

  int nstop = 0;
  ppgodunov.get("max_step",nstop);

  ppgodunov.get("redist_radius",a_params.m_redistRad);

  ppgodunov.get("cfl",a_params.m_cfl);

  ppgodunov.get("initial_cfl", a_params.m_initialDtMultiplier);

  ppgodunov.get("domain_length",a_params.m_domainLength);


  ppgodunov.get("ebvto_do_lazy_relax", a_params.m_doLazyRelax);

  // create initial/boundary condition object

  a_params.m_doSmushing = true;
  if(ppgodunov.contains("do_smushing"))
    {
      ppgodunov.get("do_smushing", a_params.m_doSmushing);
    }
  ppgodunov.get ("refine_thresh", a_params.m_refineThresh);
  ppgodunov.get("tag_buffer_size",a_params.m_tagBufferSize);

  ppgodunov.get("use_air_coefficients", a_params.m_useAirCoefs);
  a_params.m_variableCoeff = (a_params.m_useAirCoefs);
  if(!a_params.m_useAirCoefs)
    {
      ppgodunov.get("specific_heat",        a_params.m_specHeatCv );
      ppgodunov.get("thermal_conductivity", a_params.m_thermalCond);
      ppgodunov.get("mu_viscosity",         a_params.m_viscosityMu);
      ppgodunov.get("lambda_viscosity",     a_params.m_viscosityLa);
    }
  else
    {
      ppgodunov.get("specific_heat",        a_params.m_specHeatCv );
      //the others are not used
    }
  int iusemassredist;
  ppgodunov.get("use_mass_redist", iusemassredist);
  a_params.m_useMassRedist    = (iusemassredist ==1);
  fillSolverBCs(a_params, a_iprob);
  pout () << a_params << endl;
}
/***************/
void 
setupAMR(AMR&                  a_amr, 
         const EBAMRCNSParams& a_params, 
         const EBAMRCNSFactory a_amrg_fact, 
         const Box&            a_domain,
         bool                  a_fixedDt,
         Real                  a_dt,
         Vector<Vector<Box> >  a_boxes)
{
  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }

  ProblemDomain prob_domain(a_domain.smallEnd(),
                            a_domain.bigEnd(),
                            is_periodic);

  ParmParse ppgodunov;
  int max_level = 0;
  ppgodunov.get("max_level",max_level);

  int num_read_levels = Max(max_level,1);
  std::vector<int> regrid_intervals; // (num_read_levels,1);
  ppgodunov.getarr("regrid_interval",regrid_intervals,0,num_read_levels);

  int block_factor = 4;
  ppgodunov.get("block_factor",block_factor);


  int max_grid_size = 32;
  ppgodunov.get("max_grid_size",max_grid_size);

  Real fill_ratio = 0.75;
  ppgodunov.get("fill_ratio",fill_ratio);

  int checkpoint_interval = 0;
  ppgodunov.get("checkpoint_interval",checkpoint_interval);

  int plot_interval = 0;
  ppgodunov.get("plot_interval",plot_interval);


  std::vector<int> ref_ratios; // (num_read_levels,1);
  // note this requires a ref_ratio to be defined for the
  // finest level (even though it will never be used)
  ppgodunov.getarr("ref_ratio",ref_ratios,0,num_read_levels+1);

  a_amr.define(max_level,ref_ratios,prob_domain,&a_amrg_fact);

  if (a_fixedDt)
    {
      a_amr.fixedDt(a_dt);
    }

  // set grid generation parameters
  a_amr.maxGridSize(max_grid_size);
  a_amr.blockFactor(block_factor);
  a_amr.fillRatio(fill_ratio);

  // the hyperbolic codes use a grid buffer of 1
  int gridBufferSize;
  ppgodunov.get("grid_buffer_size",gridBufferSize);
  a_amr.gridBufferSize(gridBufferSize);

  // set output parameters
  a_amr.checkpointInterval(checkpoint_interval);
  a_amr.plotInterval(plot_interval);
  a_amr.regridIntervals(regrid_intervals);
  Real max_dt_growth = 1.1;
  ppgodunov.get("max_dt_growth",max_dt_growth);

  Real dt_tolerance_factor = 1.1;
  ppgodunov.get("dt_tolerance_factor",dt_tolerance_factor);

  a_amr.maxDtGrow(max_dt_growth);
  a_amr.dtToleranceFactor(dt_tolerance_factor);

  if (ppgodunov.contains("plot_prefix"))
    {
      std::string prefix;
      ppgodunov.get("plot_prefix",prefix);
      char prefixchar[100];
      sprintf(prefixchar, "%s.nx.%d.",prefix.c_str(), a_domain.size(0));
      string prefixstr(prefixchar);
      a_amr.plotPrefix(prefixstr);
    }

  if (ppgodunov.contains("chk_prefix"))
    {
      std::string prefix;
      ppgodunov.get("chk_prefix",prefix);
      char prefixchar[100];
      sprintf(prefixchar, "%s.nx.%d.",prefix.c_str(), a_domain.size(0));
      string prefixstr(prefixchar);
      a_amr.checkpointPrefix(prefixstr);
    }

  a_amr.verbosity(a_params.m_verbosity);

  if (!ppgodunov.contains("restart_file"))
    {
      if (a_boxes.size() == 0)
        {
          // initialize from scratch for AMR run
          // initialize hierarchy of levels
          a_amr.setupForNewAMRRun();
        }
      else
        {
          a_amr.setupForFixedHierarchyRun(a_boxes,1);
        }
    }
  else
    {
      std::string restart_file;
      ppgodunov.get("restart_file",restart_file);
      pout() << " restarting from file " << restart_file << endl;

#ifdef CH_USE_HDF5
      HDF5Handle handle(restart_file,HDF5Handle::OPEN_RDONLY);
      // read from checkpoint file
      a_amr.setupForRestart(handle);
      handle.close();
#else
      MayDay::Error("amrGodunov restart only defined with hdf5");
#endif
    }
}

#include "NamespaceFooter.H"
