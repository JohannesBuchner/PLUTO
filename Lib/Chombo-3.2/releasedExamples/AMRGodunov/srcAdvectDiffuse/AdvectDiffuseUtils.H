#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ADVECTDIFFUSEUTILS_H_
#define _ADVECTDIFFUSEUTILS_H_

#include <cmath>

#include "AMR.H"
#include "AMRLevelAdvectDiffuseFactory.H"
#include "AdvectTestIBC.H"
#include "AdvectionFunctions.H"
#include "memtrack.H"
#include "parstream.H"

#include "NamespaceHeader.H"

///
/**
 */
extern void
makeFinestDomain(ProblemDomain& a_domain,
                 Real&          a_dx);

///
/**
**/
extern Real
dtgNorm(const Vector< LevelData<FArrayBox>* >& a_src,
        const Vector< DisjointBoxLayout >&     a_grids,
        const Vector<int>&                     a_refRatio,
        const ProblemDomain&                   a_coarseDom,
        const int& a_comp,
        const int& a_pval);

///
/**
**/
extern void
compareError(Vector<Real>&                            a_orders,
             const Vector< LevelData<FArrayBox>* >&   a_errorFine,
             const Vector< LevelData<FArrayBox>* >&   a_errorCoar,
             const Vector< DisjointBoxLayout >&       a_gridsFine,
             const Vector< DisjointBoxLayout >&       a_gridsCoar,
             const Vector<int>&                       a_refRat,
             const ProblemDomain&                     a_coarseDom,
             int a_testverbosity);


///
/**
 */
extern void
getErrorFromCoarseAndFine(Vector< LevelData<FArrayBox>* >&           a_errorCoar,
                          const Vector< LevelData<FArrayBox>* >&     a_solnCoar,
                          const Vector< DisjointBoxLayout >&         a_gridsCoar,
                          const ProblemDomain&                       a_level0DomainCoar,
                          const Vector< LevelData<FArrayBox>* >&     a_solnFine,
                          const Vector< DisjointBoxLayout >&         a_gridsFine,
                          const ProblemDomain&                       a_level0DomainFine,
                          const Vector<int>&                         a_refRat);
///
/**
**/
extern void
coarsenBoxes(Vector< Vector<Box>      >&    a_boxesCoar,
             const Vector<Vector<Box> >&    a_boxesFine,
             int a_refToCoar);
///
/*****/
extern void
getBoxes(Vector<Vector<Box> >&   a_boxes,
         Vector<int>&            a_refRat,
         const Box&              a_domain);

///
/**
 */
extern void
getProblemDomain(ProblemDomain& a_domain);

///
/**
 */
extern void
getAdvectTestIBC(RefCountedPtr<AdvectTestIBC>&  a_advPhys);

///
/**
 */
extern void
getAMRLADFactory(RefCountedPtr<AMRLevelAdvectDiffuseFactory>&  a_fact,
                 AdvectionVelocityFunction&                    a_velFunc,
                 AdvectPhysics &                               a_advPhys);

///
/**
 */
extern void
getAdvectionVelocityFunction(AdvectionVelocityFunction& a_velFunc);

///
/**
 */
extern void
defineAMR(AMR&                                          a_amr,
          RefCountedPtr<AMRLevelAdvectDiffuseFactory>&  a_fact,
          const ProblemDomain&                          a_prob_domain,
          const Vector<int>&                            a_refRat);

///
/**
 */
extern void
setupAMRForAMRRun(AMR&                                          a_amr);

/****/

#include "NamespaceFooter.H"

#endif
