#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VCLOCALFUNCS_H_
#define _VCLOCALFUNCS_H_

#include <iostream>
using std::cerr;

#include "DisjointBoxLayout.H"
#include "AMRMultiGrid.H"
#include "VCAMRPoissonOp2.H"
#include "ParmParse.H"
#include "FArrayBox.H"
#include "functionsF_F.H"

#include "UsingNamespace.H"

///
/**
 */
class VCPoissonParameters
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
  Vector<int>  refRatio;
  ProblemDomain coarsestDomain;
  Real         coarsestDx;
  RealVect     domainLength;
  RealVect     probLo;
  RealVect     probHi;
  Real         alpha;
  Real         beta;
  int          coefficient_average_type;
  int verbosity;
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


///get stuff from input file
extern void getPoissonParameters (VCPoissonParameters&  a_params);

extern int setGrids(Vector<DisjointBoxLayout>& vectGrids,
                    VCPoissonParameters&         a_params);

extern void setRHS(LevelData<FArrayBox>&    a_rhs,
                   const RealVect&          a_dx,
                   const VCPoissonParameters& a_params);


extern
AMRLevelOpFactory<LevelData<FArrayBox> >*
defineOperatorFactory(
                      const Vector<DisjointBoxLayout>&             a_grids,
                      const Vector<ProblemDomain>&                 a_vectDom,
                      Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                      Vector<RefCountedPtr<LevelData<FluxBox> > >& a_bCoef,
                      const VCPoissonParameters&                     a_params);



void getDomainsAndDxes(  Vector<ProblemDomain>&     vectDomain,
                         Vector<Real>&              vectDx,
                         VCPoissonParameters&         a_params);


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

extern void ParseValue(Real* pos,
                       int* dir,
                       Side::LoHiSide* side,
                       Real* a_values);

extern void ParseBC(FArrayBox& a_state,
                    const Box& a_valid,
                    const ProblemDomain& a_domain,
                    Real a_dx,
                    bool a_homogeneous);


#endif
