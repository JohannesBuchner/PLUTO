#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRINSUTILS_H_
#define _AMRINSUTILS_H_

#include <cmath>
#include <fstream>

#include "Box.H"
#include "REAL.H"
#include "PolyGeom.H"
#include "EBAMRNoSubcycle.H"

#include "UsingNamespace.H"

extern void outputError(const Vector< LevelData<EBCellFAB>* >&   a_errorFine,
                        const Vector< LevelData<EBCellFAB>* >&   a_errorCoar,
                        const Vector< DisjointBoxLayout >&       a_gridsFine,
                        const Vector< DisjointBoxLayout >&       a_gridsCoar,
                        const ProblemDomain&                     a_level0DomainFine,
                        const ProblemDomain&                     a_level0DomainCoar,
                        const string&                            a_fileFine,
                        const string&                            a_fileCoar,
                        Vector<string>&                          a_names,
                        const AMRParameters&                     a_params);
///
/**
 **/
extern void
getAMRINSParameters(AMRParameters&   a_params,
                    ProblemDomain&   a_coarsestInputDomain);

///
/**
 **/
extern void
AMRINSGeometry(const AMRParameters&   a_params,
               const ProblemDomain&   a_coarsestDomain,
               int a_whichGeomForced = -1,
               EBIndexSpace* a_ebisPtr = Chombo_EBIS::instance() );
///
/**
 **/
extern void getFixedGrids(Vector<Vector<Box> >&  a_coarBoxes,
                          const AMRParameters&   a_params,
                          const ProblemDomain&   a_coarsestDomain,
                          int a_whichGeomForced = -1,
                          const EBIndexSpace* a_ebisPtr = Chombo_EBIS::instance() );


///
/**
 **/
extern void getFixedGrids(Vector<Vector<Box> >&  a_fineBoxes,
                          Vector<Vector<Box> >&  a_mediBoxes,
                          Vector<Vector<Box> >&  a_coarBoxes,
                          const AMRParameters&   a_params,
                          const ProblemDomain&   a_coarsestDomain,
                          int a_whichGeomForced = -1,
                          const EBIndexSpace* a_ebisPtr = Chombo_EBIS::instance() ) ;

///
/**
 **/
extern void
compareError(const Vector< LevelData<EBCellFAB>* >&   a_errorFine,
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
             fstream*                                 a_fout = NULL);

extern void
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
                int a_whichGeomForced = -1,
                const EBIndexSpace* a_ebisPtr = Chombo_EBIS::instance() );

Real getDistance(const RealVect& a_point,
                 const RealVect& a_lineDirection,
                 const RealVect& a_linePoint);

#endif
