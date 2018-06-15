#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _GODUNOVGEOM__
#define _GODUNOVGEOM__

#include <cmath>

#include "Box.H"
#include "REAL.H"
#include "EBAMRGodunovFactory.H"
#include "EBAMRCNSFactory.H"
#include "UnionIF.H"
#include "AMR.H"
#include "EBCNSVortexIBC.H"
#include "EBAMRCNSParams.H"
#include "EBCNSVortexIBC.H"
#include "EBAMRCNSFactory.H"

#include "NamespaceHeader.H"

extern "C"
{
  typedef void (*IrregErrorFcn)(BaseIVFAB<Real>&  a_error,
                                const IntVectSet& a_ivsIrreg,
                                const EBISBox&    a_ebisBox,
                                const RealVect&   a_cylinderAxis,
                                const Real& a_dx);
}
extern "C"
{
  extern RealVect
  getCentroidPt(const EBISBox&  a_ebisBox,
                const VolIndex& a_vof,
                const Real&     a_dx);

  extern void
  getNormalDotAxis(BaseIVFAB<Real>&  a_error,
                   const IntVectSet& a_ivsIrreg,
                   const EBISBox&    a_ebisBox,
                   const RealVect&   a_cylinderAxis,
                   const Real& a_dx);

  extern void
  getVolFrac(BaseIVFAB<Real>&  a_error,
             const IntVectSet& a_ivsIrreg,
             const EBISBox&    a_ebisBox,
             const RealVect&   a_cylinderAxis,
             const Real& a_dx);

  extern void
  getNormMinTrueNormDotAxis(BaseIVFAB<Real>&  a_error,
                            const IntVectSet& a_ivsIrreg,
                            const EBISBox&    a_ebisBox,
                            const RealVect&   a_cylinderAxis,
                            const Real&       a_dx);

  extern void
  getNormMinTrueNorm(BaseIVFAB<Real>&  a_error,
                     const IntVectSet& a_ivsIrreg,
                     const EBISBox&    a_ebisBox,
                     const RealVect&   a_cylinderAxis,
                     const Real&       a_dx);

  extern void
  getNormalDotTrueNormM1(BaseIVFAB<Real>&  a_error,
                         const IntVectSet& a_ivsIrreg,
                         const EBISBox&    a_ebisBox,
                         const RealVect&   a_cylinderAxis,
                         const Real&       a_dx);

  extern void
  getNormalDotCrossVec(BaseIVFAB<Real>&  a_error,
                       const IntVectSet& a_ivsIrreg,
                       const EBISBox&    a_ebisBox,
                       const RealVect&   a_cylinderAxis,
                       const Real&       a_dx);
}

///
/****/
extern  void
godunovFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                  const ProblemDomain&  a_coarsestDomain,
                  int                   a_max_level,
                  int                   a_max_grid_size,
                  int                   a_block_factor,
                  int                   a_verbosity,
                  std::string           a_grid_file);


extern void
getVortexIBCFactory(RefCountedPtr<EBCNSVortexIBCFactory>&  a_ibc,
                    const RealVect& a_center, const Real& a_radius);

extern void
getBoxes(Vector<Vector<Box> >&   a_boxes,
         Vector<int>&            a_refRat,
         const Box&              a_domain);

/*****/
extern void
coarsenBoxes(Vector< Vector<Box>      >&    a_boxesCoar,
             const Vector<Vector<Box> >&    a_boxesFine,
             int a_refToCoar);

///
extern void
schlichtingGeometry(Box& a_coarsestDomain,
                    RealVect& a_dx);
///
/**
**/
extern void
getEBPPFactoryXY(RefCountedPtr<EBPatchPolytropicFactory>&  a_patchGamma,
                 const  EBPhysIBCFactory*  const           a_ibc);

///
/**
 **/
extern void
getEBAMRGFactory(RefCountedPtr<EBAMRGodunovFactory>&       a_fact,
                 const  EBPatchPolytropicFactory*  const   a_patch);
///
/**
 **/
extern void
godunovGeometry(Box& a_coarsestDomain,
                RealVect& a_dx);

///
/**
 **/
extern UnionIF*
makeCrossSection(const Vector<Vector<RealVect> >& a_polygons);

///
/**
 **/
extern void
getEBAMRCNSFactory(RefCountedPtr<EBAMRCNSFactory>&                  a_fact,
                   EBAMRCNSParams&                                  a_params,
                   const RefCountedPtr<EBPatchPolytropicFactory>&   a_patch,
                   int a_iprob,
                   Real a_reynoldsNumber = -1);

///
/**
   iprob =0, no flow conditions.
   iprob =1, planar shock conditions
 */
extern void fillSolverBCs(EBAMRCNSParams& a_params, const int& a_iprob);

///
/**
   iprob =0, no flow conditions.
   iprob =1, planar shock conditions
 */
extern void fillAMRParams(EBAMRCNSParams& a_params, int a_iprob);

///
/**
   reftofinest for convergence tests.
 */
extern void  setupAMR(AMR&                  a_amr, 
                      const EBAMRCNSParams& a_params, 
                      const EBAMRCNSFactory a_amrg_fact, 
                      const Box&            a_domain,
                      bool                  a_fixedDt,
                      Real                  a_dt,
                      Vector<Vector<Box>      >  a_boxes = Vector<Vector <Box> >());

#include "NamespaceFooter.H"
#endif
