#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRGODUNOVGEOMFACTORY_H_
#define _EBAMRGODUNOVGEOMFACTORY_H_

#include <cmath>

#include "Box.H"
#include "REAL.H"
#include "RealVect.H"
#include "IntVect.H"

#include "NamespaceHeader.H"

//! \class EBAMRGodunovGeomFactory
//! This class creates computational domains for embedded boundary AMR
//! Godunov algorithms. No object of this type can be instantiated--rather,
//! one just calls its static methods to initialize geometry.
class EBAMRGodunovGeomFactory
{
  public:

  //! Create a brick of regular cells using grid information from
  //! command line parameters.
  //! \param a_domain A box that will define the coarsest computational domain.
  //! \param a_dx This will store the grid spacing at the coarsest level.
  //! \param a_refRatios This will store the refinement ratios between grid
  //!                    levels.
  static void
  createRegularCellsFromParms(Box& a_coarsestDomain,
                              RealVect& a_dx,
                              Vector<int>& a_refRatios);

  //! Create a set of regular cells within a spherical (or circular in 2D)
  //! cavity using command line parameters.
  //! \param a_domain A box that will define the coarsest computational domain.
  //! \param a_dx This will store the grid spacing at the coarsest level.
  //! \param a_refRatios This will store the refinement ratios between grid
  //!                    levels.
  static void
  createSphericalCavityFromParms(Box& a_coarsestDomain,
                                 RealVect& a_dx,
                                 Vector<int>& a_refRatios);

  //! Create a set of regular cells that is periodic in each dimension.
  //! Use command line parameters to determine the specifics.
  //! \param a_domain A box that will define the coarsest computational domain.
  //! \param a_dx This will store the grid spacing at the coarsest level.
  //! \param a_refRatios This will store the refinement ratios between grid
  //!                    levels.
  static void
  createTorusFromParams(Box& a_coarsestDomain,
                        RealVect& a_dx,
                        Vector<int>& a_refRatios);

  private:

  // --------------
  // Helper methods
  // --------------

  //! Use information in the given ParmParser to define a set of domains
  //! for a number of refinement levels.
  //! \param a_coarsestDomain This Box will be defined to the coarsest domain
  //!                         described by the given set of parameters.
  //! \param a_coarsestDx This will be set to the cell spacing at the coarsest
  //!                     grid level.
  //! \param a_finestDomain This Box will be defined to the finest domain.
  //! \param a_finestDx This will be set to the cell spacing at the finest
  //!                   grid level.
  //! \param a_maxGridSize This will store the maximum number of cells allowed
  //!                      in any direction at any refinement level.
  //! \param a_maxCoarsenings This will store the maximum number of coarsenings
  //!                         allowed to regularize the grid (or -1 if there
  //!                         is no such limit).
  //! \param a_refRatios This will store the refinement ratios for the various
  //!                    grid levels.
  static void
  createRefinedDomainsFromParms(Box& a_coarsestDomain,
                                Box& a_finestDomain,
                                int& a_maxGridSize,
                                int& a_maxCoarsenings,
                                Vector<int>& a_refRatios);

  //! Retrieve information about the basic domain geometry from the given
  //! ParmParse object.
  //! \param a_domainExtents This will be set to a RealVect containing the
  //!                        maximum spatial extents of the domain.
  //! \param a_origin This will be set to a RealVect containing the coordinates
  //!                 of the origin of the domain.
  //! \param a_coarsestDx This will be set to the cell spacing at the coarsest
  //!                     grid level.
  //! \param a_finestDx This will be set to the cell spacing at the finest
  //!                   grid level.
  //! \param a_refRatios This array defines the refinement ratios for the
  //!                    grids.
  static void
  getBasicGeometry(RealVect& a_domainExtent,
                   RealVect& a_origin,
                   RealVect& a_coarsestDx,
                   RealVect& a_finestDx,
                   const Vector<int>& a_refRatios);

  // Disallowed methods.
  EBAMRGodunovGeomFactory();
  EBAMRGodunovGeomFactory(const EBAMRGodunovGeomFactory&);
  EBAMRGodunovGeomFactory& operator=(const EBAMRGodunovGeomFactory&);
};

#include "NamespaceFooter.H"
#endif
