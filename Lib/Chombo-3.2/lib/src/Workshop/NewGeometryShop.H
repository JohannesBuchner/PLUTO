#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEWGEOMETRYSHOP_H_
#define _NEWGEOMETRYSHOP_H_

#include "REAL.H"
#include "RealVect.H"
#include "Box.H"
#include "IntVect.H"
#include "EBISBox.H"
#include "GeometryService.H"
#include "GeometryShop.H"
#include "LSquares.H"
#include "BaseIF.H"
#include "AMRIO.H"

#include "CutCellMoments.H"

#include "NamespaceHeader.H"
///
/**
   This is the base class for the workshop algorithm.
   It forms the interface between the shallow water workshop classes
   and the geometryservice class. If globalDim < 3, then this code will break
 */
class NewGeometryShop: public GeometryService
{

public:
  ///
  /**
     Define the workshop using the local geometry description
  */
  NewGeometryShop(const BaseIF        & a_baseIF,
                  const RealVect      & a_origin,
                  const RealVect      & a_vectDx,
                  const ProblemDomain & a_domain,
#if RECURSIVE_GEOMETRY_GENERATION == 0
                  const int           & a_order          = GEOMETRY_ORDER,
                  const int           & a_degreeP        = GEOMETRY_DEGREEP,
#else
                  const int           & a_orderPmax      = GEOMETRY_ORDERPMAX,
                  const int           & a_degreePmax     = GEOMETRY_DEGREEPMAX,
#endif
                  const bool          & a_useConstraints = GEOMETRY_CONSTRAINTS);

  ///
  ~NewGeometryShop();

  /**
      Return true if every cell in region is regular at the
      refinement described by dx.
  */
  bool isRegular(const Box           & a_region,
                 const ProblemDomain & a_domain,
                 const RealVect      & a_origin,
                 const Real          & a_dx) const;

  ///
  /**
      Return true if every cell in region is covered at the
      refinement described by dx.
  */
  bool isCovered(const Box           & a_region,
                 const ProblemDomain & a_domain,
                 const RealVect      & a_origin,
                 const Real          & a_dx) const;


  virtual bool isIrregular(const Box           & a_region,
                           const ProblemDomain & a_domain,
                           const RealVect      & a_origin,
                           const Real          & a_dx) const;

  virtual bool canGenerateMultiCells() const
  {
    return false;
  }

  ///
  /**
     Define the internals of the input ebisRegion.
  */
  virtual void fillGraph(BaseFab<int>        & a_regIrregCovered,
                         Vector<IrregNode>   & a_nodes,
                         const Box           & a_validRegion,
                         const Box           & a_ghostRegion,
                         const ProblemDomain & a_domain,
                         const RealVect      & a_origin,
                         const Real          & a_dx) const;

  /**
  */
  void computeVoFInternals(Real                & a_volFrac,
                           Vector<int>           a_loArc[SpaceDim],
                           Vector<int>           a_hiArc[SpaceDim],
                           Vector<Real>          a_loAreaFrac[SpaceDim],
                           Vector<Real>          a_hiAreaFrac[SpaceDim],
                           Real                & a_bndryArea,
                           RealVect            & a_normal,
                           RealVect            & a_volCentroid,
                           RealVect            & a_bndryCentroid,
                           Vector<RealVect>      a_loFaceCentroid[SpaceDim],
                           Vector<RealVect>      a_hiFaceCentroid[SpaceDim],
                           const IntVectSet    & a_ivsIrreg,
                           const VolIndex      & a_vof,
                           const ProblemDomain & a_domain,
                           const RealVect      & a_origin,
                           const Real          & a_dx,
                           const RealVect      & a_vectDx,
                           const IntVect       & a_iv) const;

  //computeVofInternals caches a_iv a m_currIv
  void settCurrIv(const IntVect& a_iv)const;

  //fillResiduals
#if RECURSIVE_GEOMETRY_GENERATION == 0
  void fillResiduals(int & a_degreeP)const;
#else
  void fillResiduals(int & a_degreePmax)const;
#endif

  //fillVolFrac
  Real fillVolFrac()const;

  //fill arc records the connectivitry between vofs
  void fillArc(Vector<int>        a_arc[SpaceDim],
               const int        & a_hilo,
               const IntVectSet & a_ivsIrreg)const;

  void fillLoArc(Vector<int>        a_loArc[SpaceDim],
                 const IntVectSet & a_ivsIrreg)const;

  void fillHiArc(Vector<int>        a_hiArc[SpaceDim],
                 const IntVectSet & a_ivsIrreg)const;

  void fillAreaFrac(Vector<Real>  a_AreaFrac[SpaceDim],
                    const int   & a_hilo)const;
  //face area:lo
  void fillLoAreaFrac(Vector<Real> a_loAreaFrac[SpaceDim])const;

  //face area:hi
  void fillHiAreaFrac(Vector<Real> a_hiAreaFrac[SpaceDim])const;

  //fill boundary area
  Real fillBndryArea()const;

  //fillNormal
  RealVect fillNormal()const;

  //fillvolCentroid:
  RealVect fillvolCentroid()const;

  //fillbndryCentroid
  RealVect fillBndryCentroid()const;

  //faceCentroids
  void fillFaceCentroid(Vector<RealVect>a_loFaceCentroid[SpaceDim],
                        const int       & a_hilo)const;

  void fillLoFaceCentroid(Vector<RealVect>a_LOFaceCentroid[SpaceDim])const;

  void fillHiFaceCentroid(Vector<RealVect>a_hiFaceCentroid[SpaceDim])const;



  void clipComputedVal(Real            &    a_volFrac,
                       Vector<Real>         a_loAreaFrac[SpaceDim],
                       Vector<Real>         a_hiAreaFrac[SpaceDim],
                       Real            &    a_bndryArea,
                       RealVect        &    a_volCentroid,
                       RealVect        &    a_bndryCentroid,
                       Vector<RealVect>     a_loFaceCentroid[SpaceDim],
                       Vector<RealVect>     a_hiFaceCentroid[SpaceDim],
                       const IntVect   &    a_iv)const;

  //takes IndexTM<Real,SpaceDim> to RealVect
  RealVect convertIndexTM2RealVect(const IndexTM<Real,SpaceDim>& a_indexTm)const;

 //takes  RealVect to IndexTM<Real,SpaceDim>
  IndexTM<Real,SpaceDim>convertRealVect2IndexTM(const RealVect& a_realVect)const;


  //takes IndexTM<int,SpaceDim> to IntVect
  IntVect convertIndexTM2IntVect(const IndexTM<int,SpaceDim>& a_indexTm)const;

 //takes  IntVect to IndexTM<int,SpaceDim>
  IndexTM<int,SpaceDim>convertIntVect2IndexTM(const IntVect& a_intVect)const;

  //converts a RealVect in physical coordinates to a RealVect in coordinates relative to a cell center
  RealVect convert2RelativeCoord(const RealVect& a_rVect)const;

  //converts a IndexTM<Real,SpaceDim> in physical coordinates to a
  //RealVect in coordinates relative to a cell center
  RealVect convert2RelativeCoord(const IndexTM<Real,SpaceDim>& a_rVect)const;
  int m_phase;

  //output residuals
#if RECURSIVE_GEOMETRY_GENERATION == 0
  void outputResidual(int & type,int & a_degreeP)const;
#else
  void outputResidual(int & type,int & a_degreePmax) const;
#endif
  void outputGradNormal()const;

private:
  int                 m_numCellsClipped;
  Real                m_threshold;

  //origin
  RealVect            m_origin;
  //dx as both a RealVect and IndexTM
  RealVect            m_vectDx;
  RvgDim              m_dxVect;

#if RECURSIVE_GEOMETRY_GENERATION == 0
  int                 m_order;
  int                 m_degreeP;
#else
  int                 m_orderPmax;
  int                 m_degreePmax;
#endif
  bool                m_useConstraints;

  //reciprocal of product of components used for scaling face areas
  Real                m_volScaleFactor;

  //dx^(SpaceDim-1)
  Real m_bndryAreaScaleFactor;

  //FAB member data to store residuals
  mutable FArrayBox m_residuals;
  mutable FArrayBox m_gradNormal;

  //used by computeVofInternals
  mutable IntVect     m_currIv;

  ProblemDomain m_domain;

  const BaseIF *m_baseIF;

#if USING_TOP_FACE_MOMENTS
  CutCellMoments<GLOBALDIM-1> m_cutCellMoments;
#else
  CutCellMoments<GLOBALDIM>   m_cutCellMoments;
#endif

  static bool s_verbose;

  /**
      Return true if every cell in region is regular at the
      refinement described by dx.
  */
  bool isRegularEveryPoint(const Box&           a_region,
                           const ProblemDomain& a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const;

  ///
  /**
      Return true if every cell in region is covered at the
      refinement described by dx.
  */
  bool isCoveredEveryPoint(const Box&           a_region,
                           const ProblemDomain& a_domain,
                           const RealVect&      a_origin,
                           const Real&          a_dx) const;


  virtual bool isIrregularEveryPoint(const Box&           a_region,
                                     const ProblemDomain& a_domain,
                                     const RealVect&      a_origin,
                                     const Real&          a_dx,
                                     const Real&          a_originVal) const ;

  int getNumCellsClipped();

  //stuff disallowed for all the usual reasons.
  NewGeometryShop()
  {
    MayDay::Abort("GeometryShop uses strong construction only");
  }
  NewGeometryShop(const NewGeometryShop& a_workshopin)
  {
    MayDay::Abort("GeometryShop disallows copy contruction");
  }
  void operator=(const NewGeometryShop& a_workshopin)
  {
    MayDay::Abort("GeometryShop disallows the assignment operator");
  }

};
#include "NamespaceFooter.H"
#endif
