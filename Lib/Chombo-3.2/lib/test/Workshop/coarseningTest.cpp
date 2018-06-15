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
#include "VoFIterator.H"
#include "AllRegularService.H"
#include "GeometryShop.H"
#include "PlaneIF.H"
#include "SphereIF.H"
#include "RhodoneaIF.H"
#include "UnionIF.H"
#include "ComplementIF.H"
#include "TransformIF.H"
#include "TiltedCylinderIF.H"
#include "IntersectionIF.H"
#include "PolyGeom.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "CH_Attach.H"

#include "UsingNamespace.H"

const bool g_verbose = true;
/***************/
// define a ramp EBIS.
/***************/
int makeGeometry(Box& domain);
/***************/
int checkCoarseAssortment(const Box& a_domain);
/***************/
/***************/
void dumpmemoryatexit();
/***************/
/***************/
RealVect getVolumeCentroid(const Box& a_domain,
                           const DisjointBoxLayout& a_dbl,
                           const EBISLayout&  a_ebisl);
/***************/
/***************/
RealVect getBoundaryCentroid(const Box& a_domain,
                             const DisjointBoxLayout& a_dbl,
                             const EBISLayout&  a_ebisl);
/***************/
/***************/
Real getUncoveredVolume(const Box& a_domain,
                        const DisjointBoxLayout& a_dbl,
                        const EBISLayout&  a_ebisl);
/***************/
/***************/
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
/***************/
int
main(int argc, char** argv)
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //begin forever present scoping trick
  {
    const char* in_file;
    if (argc < 2 || argv[1][0] == '-')
      {
        in_file = "ramp.inputs";
      }
    else
      {
        in_file = argv[1];
      }
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    int eekflag = 0;
    //define the geometry object.
    Box domain;
    eekflag =  makeGeometry(domain);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeGeometry");
      }

    //check that total coarse vol = total fine vol
    //check that total coarse vol centroid = total fine vol centroid
    eekflag = checkCoarseAssortment(domain);
    if (eekflag != 0)
      {
        pout() << "checkCoarse: eek = " << eekflag << endl;
        MayDay::Error("problem in checkEBISL");
      }

    pout() << "coarsening test passed" << endl;
  }//end scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}
/**********/
int makeGeometry(Box& a_domain)
{
  Real dx;
  RealVect origin;
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
  dx = (prob_hi-prob_lo[0])/n_cell[0];
  RealVect dxVect = dx*RealVect::Unit;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      origin[idir] = prob_lo[idir];
    }
  int verbosity = 0;
  int whichgeom;
  pp.get("which_geom",whichgeom);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();

  if (whichgeom == 0)
    {
      //allregular
      pout() << "all regular geometry" << endl;
      AllRegularService regserv;
      ebisPtr->define(a_domain, origin, dx, regserv);
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

      GeometryShop workshop(ramp,verbosity,dxVect);
      //this generates the new EBIS
      ebisPtr->define(a_domain, origin, dx, workshop);
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

      bool     insideRegular = false;
      SphereIF implicit(sphereRadius,sphereCenter,insideRegular);
      GeometryShop workshop(implicit,verbosity,dxVect);
      //this generates the new EBIS
      ebisPtr->define(a_domain, origin, dx, workshop);
    }
  else if (whichgeom == 13)
    {
      pout() << "rhodonea geometry" << endl;
      vector<Real> tmp(SpaceDim);
      pp.getarr("rhodonea_center", tmp, 0, SpaceDim);
      RealVect rhodoneaCenter;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          rhodoneaCenter[idir] = tmp[idir];
        }
      Real innerRadius;
      pp.get("inner_radius", innerRadius);

      Real outerRadius;
      pp.get("outer_radius", outerRadius);

      int frequency;
      pp.get("frequency", frequency);

      bool     insideRegular = false;

      RhodoneaIF implicit(innerRadius, outerRadius, frequency,
                          rhodoneaCenter, insideRegular);

      GeometryShop workshop(implicit,verbosity,dxVect);
      //this generates the new EBIS
      ebisPtr->define(a_domain, origin, dx, workshop);
    }
  else if (whichgeom == 18)
    {
      pout() << "Low swirl burner geometry" << endl;
      //        AttachDebugger();

      Box domain;

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

      Real outerOffset = ((prob_hi - prob_lo[0]) - outerHeight) / 2.0 + prob_lo[0];

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

      GeometryShop workshop(outside,verbosity,dxVect);

      // This generates the new EBIS
      EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      ebisPtr->define(a_domain, origin, dx, workshop);

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
Real getUncoveredVolume(const Box& a_domain,
                        const DisjointBoxLayout& a_dbl,
                        const EBISLayout&  a_ebisl)
{
  CH_assert(!a_domain.isEmpty());
  Real totalVol = 0;
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = a_ebisl[dit()];
      IntVectSet domainIVS(a_domain);
      for (VoFIterator vofit(domainIVS, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real volVoF = ebisBox.volFrac(vof);
          totalVol += volVoF;
        }
    }
  Real domainSize = Real(a_domain.numPts());
  CH_assert(domainSize >= 1.0);
#ifdef CH_MPI
  Real t=totalVol;
  MPI_Allreduce ( &t, &totalVol, 1,
                  MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
#endif
  Real normedVol = totalVol/domainSize;
  return normedVol;
}
/***************/
/***************/
RealVect getVolumeCentroid(const Box& a_domain,
                           const DisjointBoxLayout& a_dbl,
                           const EBISLayout&  a_ebisl)
{
  CH_assert(!a_domain.isEmpty());
  RealVect totalVolCent = RealVect::Zero;
  Real totalVol = 0.;
  Real dx = 1.0/Real(a_domain.size(0));
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = a_ebisl[dit()];
      IntVectSet domainIVS(a_domain);
      for (VoFIterator vofit(domainIVS, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          Real volVoF = ebisBox.volFrac(vof);
          totalVol += volVoF;

          RealVect volCentVoF = ebisBox.centroid(vof);
          //make centroid relative to global coordinate system
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              volCentVoF[idir] += (Real(iv[idir])+ 0.5);
              volCentVoF[idir] *= dx;
              volCentVoF[idir] *= volVoF;
              totalVolCent[idir] += volCentVoF[idir];
            }
        }
    }
#ifdef CH_MPI

  Real t=totalVol;
  MPI_Allreduce ( &t, &totalVol, 1,
                  MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  for (int i=0; i<CH_SPACEDIM; i++)
  {
    t=totalVolCent[i];
    MPI_Allreduce ( &t, &(totalVolCent[i]), 1,
                    MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  }
#endif

  RealVect normedVolCent = totalVolCent;
  CH_assert(totalVol > PolyGeom::getTolerance());
  normedVolCent /= totalVol;
  return normedVolCent;
}
/***************/
/***************/
RealVect getBoundaryCentroid(const Box& a_domain,
                             const DisjointBoxLayout& a_dbl,
                             const EBISLayout&  a_ebisl)
{
  CH_assert(!a_domain.isEmpty());
  RealVect totalBndryCent = RealVect::Zero;
  Real totalBndryArea = 0.0;
  int  numCellsWithBoundary = 0;
  Real dx = 1.0/Real(a_domain.size(0));
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = a_ebisl[dit()];
      IntVectSet domainIVS(a_domain);
      for (VoFIterator vofit(domainIVS, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          Real bndryAreaVoF = ebisBox.bndryArea(vof);
          if (bndryAreaVoF > 0.01*PolyGeom::getTolerance())
            {
              numCellsWithBoundary++;
              totalBndryArea += bndryAreaVoF;

              RealVect bndryCentVoF = ebisBox.bndryCentroid(vof);
              //make centroid relative to global coordinate system
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  bndryCentVoF[idir] += (Real(iv[idir])+ 0.5);
                  bndryCentVoF[idir] *= dx;
                  bndryCentVoF[idir] *= bndryAreaVoF;
                  totalBndryCent[idir] += bndryCentVoF[idir];
                }
            }
        }
    }
#ifdef CH_MPI

  Real t=totalBndryArea;
  MPI_Allreduce ( &t, &totalBndryArea, 1,
                  MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  for (int i=0; i<CH_SPACEDIM; i++)
  {
    t=totalBndryCent[i];
    MPI_Allreduce ( &t, &(totalBndryCent[i]), 1,
                    MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  }
#endif

  RealVect normedBndryCent = totalBndryCent;
  if (totalBndryArea > PolyGeom::getTolerance())
    normedBndryCent /= totalBndryArea;
  return normedBndryCent;
}
/***************/
/***************/
int checkCoarseAssortment(const Box& a_domain)
{
#ifdef CH_USE_FLOAT
  Real tolerance  = 35 * PolyGeom::getTolerance();
#else
  Real tolerance  = PolyGeom::getTolerance();
#endif
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  CH_assert(ebisPtr->isDefined());
  Box levelDomain = a_domain;
  EBISLayout ebisl;
  Vector<Box> vbox(1, levelDomain);
  Vector<int>  procAssign(1, 0);
  DisjointBoxLayout dbl(vbox, procAssign);
  int nghost = 0;
  ebisPtr->fillEBISLayout(ebisl, dbl, levelDomain, nghost);

  Real fineVolume = getUncoveredVolume(levelDomain, dbl, ebisl);
  RealVect fineVolCentroid = getVolumeCentroid(levelDomain, dbl, ebisl);
  RealVect fineBndryCentroid = getBoundaryCentroid(levelDomain, dbl, ebisl);
  if (g_verbose)
    {
      pout() << "finest level:" <<  endl;
      pout() << "volume       =\t" << fineVolume << endl;
      for (int idir=0; idir < SpaceDim; idir++)
        {
          pout() << "volCentroid("<<idir <<  ") =\t"
                 << fineVolCentroid[idir] << endl;
          pout() << "bndryCentroid("<<idir <<") =\t" <<
            fineBndryCentroid[idir] << endl << endl;;
        }
    }
  int numLevels = ebisPtr->numLevels();
  for (int ilev = 1; ilev < numLevels; ilev++)
    {
      levelDomain.coarsen(2);
      CH_assert(!levelDomain.isEmpty());
      vbox[0] = levelDomain;
      DisjointBoxLayout levelDBL(vbox, procAssign);
      EBISLayout levelEBISL;
      ebisPtr->fillEBISLayout(levelEBISL, levelDBL, levelDomain, nghost);

      Real volume =
        getUncoveredVolume(levelDomain, levelDBL, levelEBISL);
      RealVect volCentroid =
        getVolumeCentroid(levelDomain, levelDBL, levelEBISL);
      RealVect bndryCentroid =
        getBoundaryCentroid(levelDomain, levelDBL, levelEBISL);
      if (g_verbose)
        {
          pout() << "ilev         =\t" << ilev  << endl;
          pout() << "volume       =\t" << volume << endl;
          for (int idir=0; idir < SpaceDim; idir++)
            {
              pout() << "volCentroid("<<idir <<  ") =\t"
                     << volCentroid[idir] << endl;
              pout() << "bndryCentroid("<<idir <<") =\t" <<
                bndryCentroid[idir] << endl << endl;;
            }
        }
      Real volDiff = Abs(volume - fineVolume);
      if (volDiff > tolerance)
        {
          pout() << "volumes don't match at level = " << ilev << ", "
                 << volDiff << " > " << tolerance << endl;
          return -1;
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real volCentDiff = Abs(volCentroid[idir] - fineVolCentroid[idir]);
          if (volCentDiff > tolerance)
            {
              pout() << "vol centroids don't match at level = "
                     << ilev << "  idir = " << idir << ", "
                     << volCentDiff << " > " << tolerance << endl;
              return -2;
            }
          Real bndryCentDiff = Abs(bndryCentroid[idir] - fineBndryCentroid[idir]);
          if (bndryCentDiff > tolerance)
            {
              pout() << "bndry centroids don't match at level = "
                     << ilev << "  idir = " << idir << ", "
                     << bndryCentDiff << " > " << tolerance << endl;
              return -3;
            }
        }
    }
  return 0;
}
/***************/
/***************/
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

  Real sinTheta = sin(angle);
#if CH_SPACEDIM == 3
  Real cosTheta = cos(angle);

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
/***************/
