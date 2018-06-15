#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBAMRGodunovGeomFactory.H"
#include "EBIndexSpace.H"
#include "ParmParse.H"
#include "AllRegularService.H"
#include "SphereIF.H"
#include "GeometryShop.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
void
EBAMRGodunovGeomFactory::
createRefinedDomainsFromParms(Box& a_coarsestDomain,
                              Box& a_finestDomain,
                              int& a_maxGridSize,
                              int& a_maxCoarsenings,
                              Vector<int>& a_refRatios)
{
  // Retrieve command line parameters.
  ParmParse parms;

  int maxRefinementLevel;
  parms.get("max_level", maxRefinementLevel);

  // If we find ref_ratio, read its contents into a_refRatios.
  if (maxRefinementLevel > 0)
  {
    if (parms.contains("ref_ratio"))
      parms.getarr("ref_ratio", a_refRatios,0, maxRefinementLevel+1);

    // Otherwise, assume that the refinement ratio is 2 for each level.
    else
      a_refRatios = Vector<int>(maxRefinementLevel+1, 2);
  }

  // Now figure out how many cells go in the coarsest box.
  Vector<int> nCell(SpaceDim);
  parms.getarr("n_cell",nCell,0,SpaceDim);
  IntVect lo = IntVect::Zero, hi;
  for (int ivec = 0; ivec < SpaceDim; ivec++)
  {
    if (nCell[ivec] <= 0)
    {
      pout() << "EBAMRGodunovGeomFactory: bad number of cells: " << nCell[ivec];
      MayDay::Error();
    }

    // Currently, we require the same number of cells in each direction.
    if (nCell[ivec] != nCell[0])
    {
      pout() << "EBAMRGodunovGeomFactory: n_cell should be equal in all directions";
      MayDay::Error();
    }

    hi[ivec] = nCell[ivec] - 1;
  }

  // Set up the coarse box.
  a_coarsestDomain = Box(lo, hi);

  // If we have AMR enabled, refine.
  a_finestDomain = a_coarsestDomain;
  for (int ilev = 0; ilev < maxRefinementLevel; ++ilev)
    a_finestDomain.refine(a_refRatios[ilev]);

  // Get the maximum grid size and make sure that it squares with our
  // parameters.
  parms.get("max_grid_size", a_maxGridSize);

  // If a maximum number of coarsenings is given, use it.
  if (parms.contains("max_coarsen"))
    parms.get("max_coarsen", a_maxCoarsenings);
  else
    a_maxCoarsenings = -1;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBAMRGodunovGeomFactory::
getBasicGeometry(RealVect& a_domainExtent,
                 RealVect& a_origin,
                 RealVect& a_coarsestDx,
                 RealVect& a_finestDx,
                 const Vector<int>& a_refRatios)
{
  // Retrieve command line parameters.
  ParmParse parms;

  // Get the extents of the domain.
  RealVect domainExtents = RealVect::Unit; // Unit square by default.
  if (parms.contains("domain_length")) // Square domain
  {
    parms.get("domain_length", domainExtents[0]);
    for (int idir = 1; idir < SpaceDim; ++idir)
      domainExtents[idir] = domainExtents[0];
  }
  else if (parms.contains("domain_extents"))
  {
    Vector<Real> extents(SpaceDim);
    parms.getarr("domain_extents", extents, 0, SpaceDim);
    for (int idir = 0; idir < SpaceDim; ++idir)
      domainExtents[idir] = extents[idir];
  }

  // Get the number of grid cells again.
  Vector<int> nCell(SpaceDim);
  parms.getarr("n_cell",nCell,0,SpaceDim);

  // Determine the cell spacing at the coarsest grid level.
  for (int idir = 0; idir< SpaceDim; idir++)
    a_coarsestDx[idir] = domainExtents[idir]/nCell[idir];

  // Determine the finest grid spacing.
  a_finestDx = a_coarsestDx;
  for (int ilev = 0; ilev < a_refRatios.size(); ++ilev)
    a_finestDx /= a_refRatios[ilev];

  // If the origin is given, use it. Otherwise assume zero.
  if (parms.contains("origin"))
  {
    Vector<Real> O(SpaceDim);
    parms.getarr("origin", O, 0, SpaceDim);
    for (int idir = 0; idir < SpaceDim; ++idir)
      a_origin[idir] = O[idir];
  }
  else
    a_origin = RealVect::Zero;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBAMRGodunovGeomFactory::
createRegularCellsFromParms(Box& a_coarsestDomain,
                            RealVect& a_coarsestDx,
                            Vector<int>& a_refRatios)
{
  // Retrieve grid data.
  Box finestDomain;
  int maxGridSize, maxCoarsenings;
  EBAMRGodunovGeomFactory::createRefinedDomainsFromParms(a_coarsestDomain,
                                                         finestDomain,
                                                         maxGridSize,
                                                         maxCoarsenings,
                                                         a_refRatios);

  // Retrieve basic geometric information.
  RealVect domainExtents;
  RealVect origin, finestDx;
  EBAMRGodunovGeomFactory::getBasicGeometry(domainExtents, origin,
                                            a_coarsestDx, finestDx,
                                            a_refRatios);

  // now define the regular geometry.
  AllRegularService regServ;
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(finestDomain, origin, finestDx[0], regServ, maxGridSize,
                  maxCoarsenings);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBAMRGodunovGeomFactory::
createSphericalCavityFromParms(Box& a_coarsestDomain,
                               RealVect& a_coarsestDx,
                               Vector<int>& a_refRatios)
{
  // Retrieve grid data.
  Box finestDomain;
  int maxGridSize, maxCoarsenings;
  EBAMRGodunovGeomFactory::createRefinedDomainsFromParms(a_coarsestDomain,
                                                         finestDomain,
                                                         maxGridSize,
                                                         maxCoarsenings,
                                                         a_refRatios);

  // Retrieve basic geometric information.
  RealVect domainExtents;
  RealVect origin, finestDx;
  EBAMRGodunovGeomFactory::getBasicGeometry(domainExtents, origin,
                                            a_coarsestDx, finestDx,
                                            a_refRatios);

  // now define the cavity.
  ParmParse pp;
  vector<Real> x0(SpaceDim);
  pp.getarr("sphere_center",x0, 0, SpaceDim);
  RealVect sphereCenter;
  for (int i = 0; i < SpaceDim; ++i)
    sphereCenter[i] = x0[i];
  Real sphereRadius;
  pp.get("sphere_radius", sphereRadius);
  SphereIF sphereIF(sphereRadius, sphereCenter, true);
  int verbosity = 0;
  GeometryShop workshop(sphereIF,verbosity,finestDx);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(finestDomain, origin, finestDx[0], workshop, maxGridSize,
                  maxCoarsenings);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBAMRGodunovGeomFactory::
createTorusFromParams(Box& a_coarsestDomain,
                      RealVect& a_coarsestDx,
                      Vector<int>& a_refRatios)
{
  // Retrieve grid data.
  Box finestDomain;
  int maxGridSize, maxCoarsenings;
  EBAMRGodunovGeomFactory::createRefinedDomainsFromParms(a_coarsestDomain,
                                                         finestDomain,
                                                         maxGridSize,
                                                         maxCoarsenings,
                                                         a_refRatios);

  // Retrieve basic geometric information.
  RealVect domainExtents;
  RealVect origin, finestDx;
  EBAMRGodunovGeomFactory::getBasicGeometry(domainExtents, origin,
                                            a_coarsestDx, finestDx,
                                            a_refRatios);

  // What do we have to do here???
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
