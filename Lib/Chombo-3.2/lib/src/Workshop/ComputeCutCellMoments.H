#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COMPUTECUTCELLMOMENTS_H_
#define _COMPUTECUTCELLMOMENTS_H_

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include <map>
using std::map;

#include "Vector.H"
#include "REAL.H"
#include "IndexTM.H"

#include "Notation.H"
#include "LSquares.H"
#include "IFData.H"
#include "CutCellMoments.H"
#include "RefinementCriterion.H"

#include "NamespaceHeader.H"

template <int dim> class LSProblem;

template <int dim> class ComputeCutCellMoments
{
public:
  typedef IndexTM<int,dim>  IvDim;
  typedef IndexTM<Real,dim> RvDim;

  typedef map<IvDim,Real,LexLT <IvDim > >  PthMoment;

  typedef map<IndexTM<int,dim-1>,Real,LexLT < IndexTM<int,dim-1> > > PthMomentLesserDimension;

  typedef map<IndexTM<int,1>,Real >     OneDMoments;
  typedef map<int,IvDim>                LocPthMoment;
  typedef map<IvDim,int,LexLT <IvDim> > PthMomentLoc;

  typedef map<Iv2,CutCellMoments<dim-1>, LexLT<Iv2> > BdCutCellMoments;

  // Constructors
  ComputeCutCellMoments();
  ComputeCutCellMoments(const ComputeCutCellMoments<dim>& a_computeCutCellMoments);

  // This is used to build bd ComputeCutCellMoments
  ComputeCutCellMoments(const IFData<dim>& a_info);

  // Destructor
  ~ComputeCutCellMoments();

#if RECURSIVE_GEOMETRY_GENERATION == 0
  void computeMoments(const int                & a_order,
                      const int                & a_degreeP,
                      const bool               & a_useConstraints,
                      RefinementCriterion<dim> & a_refinementCriterion,
                      const int                & a_numberOfRefinements = 0);
#else
  void computeMoments(const int                & a_orderPmax,
                      const int                & a_degreePmax,
                      const bool               & a_useConstraints,
                      RefinementCriterion<dim> & a_refinementCriterion,
                      const int                & a_numberOfRefinements = 0);

  void computeMomentsRecursively(const int                & a_orderPmax,
                                 const int                & a_degreePmax,
                                 const bool               & a_useConstraints,
                                 RefinementCriterion<dim> & a_refinementCriterion,
                                 const int                & a_numberOfRefinements);
#endif

  Vector<Real> computeRhs(LSProblem<dim> & a_lsp,
                          const int      & a_order);

  // Solve the problem on a refined cell
#if RECURSIVE_GEOMETRY_GENERATION == 0
  Vector< CutCellMoments<dim> > refine(const int                & a_order,
                                       const int                & a_degreeP,
#else
  Vector< CutCellMoments<dim> > refine(const int                & a_orderPmax,
                                       const int                & a_degreePmax,
#endif
                                       const bool               & a_useConstraints,
                                       RefinementCriterion<dim> & a_refinementCriterion,
                                       const IndexTM<int,dim>   & a_refineInDir,
                                       const int                & a_numberOfRefinements);

  void addMomentMaps(const Vector<CutCellMoments<dim> > & a_refinedCutCellVector,
#if RECURSIVE_GEOMETRY_GENERATION == 0
                     const int                          & a_degreeP,
#else
                     const int                          & a_degreePmax,
#endif
                     const bool                         & a_useConstraints);

  void addMoments(PthMoment               & a_momentMap,
                  PthMoment               & a_refinedMomentMap,
                  const IndexTM<Real,dim> & a_refinedCenterDelta);

  void addBdMoments(CutCellMoments<dim>     & a_coarseCutCell,
                    const IFData<dim+1>     & a_IFData,
#if RECURSIVE_GEOMETRY_GENERATION == 0
                    const int               & a_degreeP,
#else
                    const int               & a_degreePmax,
#endif
                    const bool              & a_useConstraints,
                    const IndexTM<Real,dim> & a_refinedCenterDelta,
                    const IndexTM<int,dim>  & a_localHilo);

#if RECURSIVE_GEOMETRY_GENERATION == 0
  void computeResiduals(const int  & a_order,
                        const int  & a_degreeP,
#else
  void computeResiduals(const int  & a_orderPmax,
                        const int  & a_degreePmax,
#endif
                        const bool & a_useConstraints);

  void computeResiduals(const Vector<CutCellMoments<dim> > & a_refinedCCMoms,
#if RECURSIVE_GEOMETRY_GENERATION == 0
                        const int & a_degreeP);
#else
                        const int & a_degreePmax);
#endif

  // Output
  void print(ostream& out) const;

  void dump() const;

  // Operators
  void operator=(const ComputeCutCellMoments<dim>& a_computeCutCellMoments);

  Real factorial(const IvDim & a_multiIndex) const;

  // All the moments (and other data)
  CutCellMoments<dim> m_cutCellMoments;

  // Flag whether the boundary moments have been computed already
  bool m_boundaryMomentsComputed;
};

// One dimensional ComputeCutCellMoments
template <> class ComputeCutCellMoments<1>
{
public:
  typedef map<IndexTM<int,1>,Real> OneDMoments;

  // Constructors
  ComputeCutCellMoments();
  ComputeCutCellMoments(const ComputeCutCellMoments<1> & a_computeCutCellMoments);
  ComputeCutCellMoments(const IFData<1>& a_info);

  // Destructor
  ~ComputeCutCellMoments();

#if RECURSIVE_GEOMETRY_GENERATION == 0
  void computeMoments(const int              & a_order,
                      const int              & a_degree,
#else
  void computeMoments(const int              & a_orderPmax,
                      const int              & a_degreePmax,
#endif
                      const bool             & a_useConstraints,
                      RefinementCriterion<1> & a_refinementCriterion);

  void simpleComputeMoments(const Real & a_loPt,
                            const Real & a_hiPt,
#if RECURSIVE_GEOMETRY_GENERATION == 0
                            const int  & a_degree);
#else
                            const int  & a_degreePmax);
#endif

  void computeMomentsUsingBinomial(const Real & a_loPt,
                                   const Real & a_hiPt,
                                   const int  & a_loSign,
                                   const int  & a_hiSign,
#if RECURSIVE_GEOMETRY_GENERATION == 0
                                   const int  & a_degree);
#else
                                   const int  & a_degreePmax);
#endif

  // Output
  void print(ostream& out) const;

  void dump() const;

  // Operators
  void operator=(const ComputeCutCellMoments<1>& a_computeCutCellMoments);

  // Member data
  CutCellMoments<1> m_cutCellMoments;
};

#include "NamespaceFooter.H"

#include "ComputeCutCellMomentsImplem.H"

#endif
