#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NODEPOISSONUTILITIES_H_
#define _NODEPOISSONUTILITIES_H_

#include <iostream>
using std::cerr;

#include "DisjointBoxLayout.H"
#include "AMRMultiGrid.H"
#include "RelaxSolver.H"
#include "ParmParse.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "NodeFArrayBox.H"
#include "PetscSolver.H"

#include "UsingNamespace.H"

///
/**
 */
class PoissonParameters
{
public:
  IntVect      nCells;
  int          maxGridSize;
  int          blockFactor;
  int          bufferSize;
  Real         fillRatio;
  Real         refineThresh;
  int          domBcType;
  int          maxLevel;
  int          numLevels;
  int          verbosity;
  Vector<int>  refRatio;
  ProblemDomain coarsestDomain;
  Real         coarsestDx;
  RealVect     domainLength;
  RealVect     probLo;
  RealVect     probHi;
  Real         alpha;
  Real         beta;
  void coarsen(int a_factor);
  void  refine(int a_factor);
};

///
/**
 */
class GlobalBCRS
{
public:
  static std::vector<bool> s_printedThatLo, s_printedThatHi;
  static std::vector<int> s_bcLo, s_bcHi;
  static RealVect s_trigvec;
  static bool s_areBCsParsed, s_valueParsed, s_trigParsed;
};

extern RealVect& getTrigRV();

extern void
outputData(const Vector<LevelData<FArrayBox>* >& vectPhi,
           const Vector<DisjointBoxLayout>& vectGrids,
           const ProblemDomain& coarsestDomain,
           const Vector<int>& vectRatio,
           Real dxCoarsest,
           int numlevels,
           string filename,
           string varname);

extern void
outputVectorData(const Vector<LevelData<FArrayBox>* >& vectPhi,
                 const Vector<DisjointBoxLayout>& vectGrids,
                 const ProblemDomain& coarsestDomain,
                 const Vector<int>& vectRatio,
                 Real dxCoarsest,
                 int numlevels,
                 string filename,
                 Vector<string>& a_vectName);

extern void
nodeOutputData(const Vector<LevelData<NodeFArrayBox>* >& vectPhi,
               const Vector<DisjointBoxLayout>& vectGrids,
               const ProblemDomain& coarsestDomain,
               const Vector<int>& vectRatio,
               Real dxCoarsest,
               int numlevels,
               string filename,
               string varname);

///get stuff from input file
extern void getPoissonParameters (PoissonParameters&  a_params);

extern int setGrids(Vector<DisjointBoxLayout>& vectGrids,
                    PoissonParameters&         a_params);

/// if numlevels is negative, use numlevels from params
extern int setRHS(Vector<LevelData<FArrayBox>* >& vectRhs,
                  PoissonParameters&              a_params,
                  int numlevels = -1);


void getDomainsAndDxes(  Vector<ProblemDomain>&     vectDomain,
                         Vector<Real>&              vectDx,
                         PoissonParameters&         a_params);

extern void
tagCells(Vector<LevelData<FArrayBox>* >& vectRHS,
         Vector<IntVectSet>& tagVect,
         Vector<Real>& vectDx,
         Vector<ProblemDomain>& vectDomain,
         const Real refine_thresh,
         const int tags_grow,
         const int baseLevel,
         int numLevels) ;

extern void TrigValueNeum(Real* pos,
                          int* dir,
                          Side::LoHiSide* side,
                          Real* a_values);

extern void TrigValueDiri(Real* pos,
                          int* dir,
                          Side::LoHiSide* side,
                          Real* a_values);

extern void ResistDiri(Real* pos,
                       int* dir,
                       Side::LoHiSide* side,
                       Real* a_values);

extern void ParseValue(Real* pos,
                       int* dir,
                       Side::LoHiSide* side,
                       Real* a_values);

extern void ParseBC(FArrayBox& a_state,
                    const Box& a_valid,
                    const ProblemDomain& a_domain,
                    Real a_dx,
                    bool a_homogeneous);

extern void NodeParseBC(NodeFArrayBox& a_state,
                        const Box& a_valid,
                        const ProblemDomain& a_domain,
                        Real a_dx,
                        bool a_homogeneous);
///
/**
 */
extern void compareError(const Vector< LevelData<FArrayBox>* >&   a_errorFine,
                         const Vector< LevelData<FArrayBox>* >&   a_errorCoar,
                         const Vector< DisjointBoxLayout >&       a_gridsFine,
                         const Vector< DisjointBoxLayout >&       a_gridsCoar,
                         const PoissonParameters&                 a_paramsFine,
                         const PoissonParameters&                 a_paramsCoar,
                         const string& a_testName);

///
/**
 */
extern void getCoarseLayoutsFromFine(Vector<DisjointBoxLayout>&       a_gridsCoar,
                                     const Vector<DisjointBoxLayout>& a_gridsFine,
                                     const PoissonParameters&         a_params);

/********/
extern void
defineSolver(AMRMultiGrid<LevelData<FArrayBox> >&         a_solver,
             const Vector<DisjointBoxLayout>&             a_grids,
             LinearSolver<LevelData<FArrayBox> >&         a_bottomSolver,
             const PoissonParameters&                     a_params);


extern void
nodeDefineSolver(AMRMultiGrid<LevelData<NodeFArrayBox> >&         a_solver,
                 const Vector<DisjointBoxLayout>&                 a_grids,
                 LinearSolver<LevelData<NodeFArrayBox> >&         a_bottomSolver,
                 const PoissonParameters&                         a_params);

///
/**
 */
extern void setTrigPhi(LevelData<FArrayBox>&     a_phi,
                       const  RealVect&          a_dx,
                       const  PoissonParameters& a_params);

extern void setTrigKappaLOfPhi(LevelData<FArrayBox>&    a_kappaLOfPhi,
                               const RealVect&          a_dx,
                               const PoissonParameters& a_params);



extern void setKLVViscous  (LevelData<FArrayBox>&    a_kappaLOfMag,
                            const RealVect&          a_dx,
                            const PoissonParameters& a_params);

extern void setVelViscous  (LevelData<FArrayBox>&    a_kappaLOfMag,
                            const RealVect&          a_dx,
                            const PoissonParameters& a_params);
extern void
defineViscousTensorCoef(Vector<RefCountedPtr<LevelData<FluxBox> > >&    a_eta,
                        Vector<RefCountedPtr<LevelData<FluxBox> > >&    a_lambda,
                        const Vector<DisjointBoxLayout>&                a_grids,
                        const PoissonParameters&                        a_params);
extern void
defineViscousTensorSolver(AMRMultiGrid<LevelData<FArrayBox> >&         a_solver,
                          const Vector<DisjointBoxLayout>&             a_grids,
                          LinearSolver<LevelData<FArrayBox> >&         a_bottomSolver,
                          const PoissonParameters&                     a_params);

extern int
definePETSCViscousTensorSolver(PetscSolverViscousTensor<LevelData<FArrayBox> >*& a_petscSolver,
                               const Vector<DisjointBoxLayout>&             a_grids,
                               const PoissonParameters&                     a_params);


#endif
