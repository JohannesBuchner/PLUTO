#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include <cmath>

#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "DebugDump.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "PolyGeom.H"
#include "GeometryShop.H"

#include "EBFABView.H"

#include "SphereIF.H"
#include "TiltedCylinderIF.H"
#include "UnionIF.H"
#include "IntersectionIF.H"
#include "ComplementIF.H"

#include "EBMenagerieUtils.H"

#include "UsingNamespace.H"

/***************/
// Define an EBIS.
/***************/
BaseIF* makeGeometry(Box&      a_domain,
                     RealVect& a_origin,
                     Real&     a_dx);

/***************/
// Get the distance from a point to a line
/***************/
Real getDistance(const RealVect& a_point,
                 const RealVect& a_lineDirection,
                 const RealVect& a_linePoint);

int main (int argc, char** argv)
{
#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Begin forever present scoping trick
  {
    const char* in_file = "packedChannel.inputs";

    if (argc > 1)
    {
      in_file = argv[1];
    }

    // Parse input file
    ParmParse pp(0,NULL,NULL,in_file);

    Box domain;
    RealVect origin;
    Real dx;

    BaseIF* implicit;

    // Make geometry
    implicit = makeGeometry(domain,origin,dx);

    // Make grids
    DisjointBoxLayout grids;
    makeLayout(grids,domain);

    // Create ebislayout
    int nghost = 0;
    EBISLayout ebisl;
    makeEBISL(ebisl,grids,domain,nghost);

    // Make a LevelData
    int nComps = 1;

    IntVect ghost = IntVect::Unit;
    ghost *= nghost;

    RefCountedPtr<DataFactory<EBCellFAB> > rcpFactory(new EBCellFactory(ebisl));
    LevelData<EBCellFAB> level(grids,nComps,ghost,*rcpFactory);

    // Put some data in the data holders
    fillData(level,origin,dx,*implicit);

    // Done with this object
    delete implicit;

    // Write the data and the EB out
    const char* basename = "packedChannel";

    char name[1000];
    sprintf(name,"%s%dd.hdf5",basename,SpaceDim);

#ifdef CH_USE_HDF5
    writeEBLevelname(&level,name);
#endif
  } // End scoping trick

  // Clean up
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();

  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
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

BaseIF* makeGeometry(Box&      a_domain,
                     RealVect& a_origin,
                     Real&     a_dx)
{
  RealVect cylDirection;
  RealVect cylPoint;
  Real     cylRadius;

  int      sphNumber;
  Real     sphMinRadius;
  Real     sphMaxRadius;

  // Initialize random number generator
  srand48(time(0));

  // parse input file
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
      pout() << "Bogus number of cells input = " << n_cell[ivec];
      exit(1);
    }

    hi[ivec] = n_cell[ivec] - 1;
  }

  a_domain.setSmall(lo);
  a_domain.setBig(hi);

  Vector<Real> prob_lo(SpaceDim,1.0);
  Real prob_hi;

  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);

  RealVect domainLength;

  domainLength[0] = prob_hi - prob_lo[0];

  a_dx = domainLength[0] / n_cell[0];

  for (int idir = 1; idir < SpaceDim; idir++)
  {
    domainLength[idir] = n_cell[idir] * a_dx;
  }

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    a_origin[idir] = prob_lo[idir];
  }

  bool useDomainBoundary = false;
  pp.query("useDomainBoundary",useDomainBoundary);

  int packingType = 1;
  pp.get("packingType",packingType);

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

  Real minDistCyl;
  pp.get("minDistCyl",minDistCyl);

  Real minDistSph;
  pp.get("minDistSph",minDistSph);

  Real minDistInlet;
  pp.get("minDistInlet",minDistInlet);

  Real minDistOutlet;
  pp.get("minDistOutlet",minDistOutlet);

  Real maxTries;
  pp.get("maxTries",maxTries);

  Real maxDrops;
  pp.get("maxDrops",maxDrops);

  // Stay inside until the end
  bool inside = true;

  Vector<RealVect> centerArray(sphNumber);
  Vector<Real>     radiusArray(sphNumber);

  int randSeed;
  pp.get("randSeed",randSeed);

  srand48(randSeed);

  // Place random spheres inside the cylinder
  int packed = 0;
  Real count = 0;

  RealVect curDirection = cylDirection;

  if (curDirection.dotProduct(BASISREALV(0)) > 0)
  {
    curDirection = -curDirection;
  }

  curDirection /= curDirection.vectorLength();

  Real a = curDirection.dotProduct(curDirection);

  MAXHEIGHT* sortedList = new MAXHEIGHT[sphNumber];

  while (packed < sphNumber && count < maxTries)
  {
    if (packingType == 1)
    {
      RealVect& center = centerArray[packed];
      Real&     radius = radiusArray[packed];

      // Get a random center and radius for a sphere
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          center[idir] = domainLength[idir]*drand48() + a_origin[idir];
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
      if (center[0] - radius - minDistInlet < a_origin[0])
      {
        isOkay = false;
      }

      if (center[0] + radius + minDistOutlet > a_origin[0] + domainLength[0])
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
    {
      RealVect& bestCenter = centerArray[packed];
      Real&     bestRadius = radiusArray[packed];
      bool      bestOkay = false;

      // Get a random radius for the next sphere
      Real radius;
      radius = (sphMaxRadius - sphMinRadius)*drand48() + sphMinRadius;

      for (double tries = 0; tries < maxDrops; tries++)
      {
        RealVect center;
        center[0] = 2*domainLength[0] + a_origin[0];

        // Get a random center
        for (int idir = 1; idir < SpaceDim; idir++)
          {
            center[idir] = domainLength[idir]*drand48() + a_origin[idir];
        }

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
          Real minCenterX = a_origin[0];

          // Compute the minimum distance it can "fall"
          for (int i = packed-1; i >= 0; i--)
          {
            if (foundMin && (sortedList[i].max < (minCenterX - radius - minDistSph)))
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
            minS = (a_origin[0] + radius + minDistInlet - center[0]) / curDirection[0];
            minCenterX = center[0] + minS * curDirection[0];
          }

          center += minS * curDirection;
        }

        // Check against high x domain boundary
        if (center[0] + radius + minDistOutlet > a_origin[0] + domainLength[0])
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

  delete [] sortedList;

  pout() << "Tried " << count
         << ", packed " << packed
         << ", " << sphNumber-packed << " left"
         << endl;

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
  IntersectionIF cylinderAndSpheres(cylinder,outsideSpheres);

  BaseIF* finalIF = NULL;

  if (useDomainBoundary)
    {
      // Bounded by problem domain
      finalIF = &outsideSpheres;
    }
  else
    {
      // Bounded by a cylinder
      finalIF = &cylinderAndSpheres;
    }

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(*finalIF,0,vectDx);

  // This generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx,workshop);

  return finalIF->newImplicitFunction();
}

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
