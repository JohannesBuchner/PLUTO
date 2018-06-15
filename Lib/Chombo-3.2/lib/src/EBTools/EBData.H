#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#ifndef _EBDATA_H_
#define _EBDATA_H_

#include "REAL.H"
#include "RealVect.H"

#include "EBGraph.H"
#include "IrregNode.H"
#include "BaseIVFAB.H"
#include "BaseIFFAB.H"

#include "NamespaceHeader.H"

class EBIndexSpace;

template <class T>
int linearSize(const T& inputT);

template <class T>
void linearIn(T& a_outputT, const void* const inBuf);

template <class T>
void linearOut(void* const a_outBuf, const T& inputT);

class BoundaryData
{
public:
  BoundaryData();

  BoundaryData(const Real& bndryArea, const RealVect& normal,
               const RealVect& bndryCentroid, int phase,const VolIndex& index)
    :
    m_bndryArea(bndryArea),
    m_normal(normal),
    m_bndryCentroid(bndryCentroid),
    m_bndryPhase(phase),
    m_volIndex(index)
  {
  }

  Real     m_bndryArea;
  RealVect m_normal;
  RealVect m_bndryCentroid;
  int      m_bndryPhase;
  VolIndex m_volIndex;
};

///
/**
 */
class VolData
{
public:
  VolData();
  Real                 m_volFrac;
  RealVect             m_volCentroid;
  BoundaryData         m_averageFace;
  Vector<BoundaryData> m_phaseFaces;
};

/// VolData specializations of linear* functions used by LevelData
//[NOTE: see also SPMD.H. <dbs>]
template < >
int linearSize(const VolData& vdata);
template < >
void linearIn(VolData& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const VolData& a_inputT);



std::ostream& operator<< (std::ostream&  a_os,
                          const VolData& a_iv);
///
/**
 */
class FaceData
{
public:
  Real     m_areaFrac;
  RealVect m_faceCentroid;
};

template < >
int linearSize(const FaceData& vdata);
template < >
void linearIn(FaceData& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const FaceData& a_inputT);

///
/**
   This class contains all the geometric information
   for an ebisbox.
 */
class EBDataImplem
{

public:

  ///
  EBDataImplem();

  ///
  ~EBDataImplem();

  ///
  /**
     Noop.
  */
  void define(const Box& box, int comps);

  ///
  /**
     Noop
  */
  EBDataImplem(const Box& a_box, int a_comps);

  ///
  /**
     Copy the information from a_source to the over the intersection
     of the box a_region, the box of the current EBDataImplem and
     the box of a_source.  The Interval arguments are ignored.
     This function is required by LevelData.
  */
  void copy(const Box&      a_regionFrom,
            const Interval& a_Cd,
            const Box&      a_regionto,
            const EBDataImplem&   a_source,
            const Interval& a_Cs);

  /// each data holder is defined over the irregular cells of the graph
  void
  defineVoFData(const EBGraph& a_graph, const Box& a_region);

  /// each data holder is defined over the irregular cells of the graph
  void
  defineFaceData(const EBGraph& a_graph, const Box& a_region);

  ///define the whole thing
  void
  define(const EBGraph&           a_graph,
         const Vector<IrregNode>& a_irregData,
         const Box&               a_validBox);

  ///
  void
  coarsenVoFs(const EBDataImplem&   a_fineEBDataImplem,
              const EBGraph&        a_fineGraph,
              const EBGraph&        a_coarGraph,
               const Box&           a_validRegion);

  ///
  void
  coarsenFaces(const EBDataImplem& a_fineEBDataImplem,
               const EBGraph&      a_fineGraph,
               const EBGraph&      a_coarGraph,
               const Box&          a_validRegion);
  ///
  const Real& volFrac(const VolIndex& a_vof) const;

  ///
  const Real& areaFrac(const FaceIndex& a_face1) const;

  ///
  const RealVect& centroid(const FaceIndex& facein) const;

  ///
  const RealVect& centroid(const VolIndex& a_vof) const;

  ///
  const RealVect& bndryCentroid(const VolIndex& a_vof) const;
  const RealVect& bndryCentroid(const VolIndex& a_vof, int face) const;

  ///
  const Real& bndryArea(const VolIndex& a_vof) const;
  const Real& bndryArea(const VolIndex& a_vof, int face) const;

  ///
  const RealVect& normal(const VolIndex& a_vof) const;
  const RealVect& normal(const VolIndex& a_vof, int face) const;

  //  New set of boundary inquiry functions to support multi-fluid applications
  /// used by multi-fluid applications

  /// used by multi-fluid code
  int facePhase(const VolIndex& a_vof, int face) const ;

  /// used by multi-fluid code
  const VolIndex& faceIndex(const VolIndex& a_vof, int face) const ;

  /// used by multi-fluid code
  void setFacePhase(const VolIndex& a_vof, int face, int phase);
  ///
  /// used by multi-fluid code
  void setFaceIndex(const VolIndex& a_vof, int face, const VolIndex& index);

  ///
  int numFacePhase(const VolIndex& a_vof) const ;

  void clearMultiBoundaries();
  void setBoundaryPhase(int phase);
  /**
     This stuff required by LevelData in parallel:
  */
  int size(const Box& R, const Interval& comps) const;

  ///
  void linearOut(void* buf, const Box& R, const Interval& comps) const;

  ///
  void linearIn(void* buf, const Box& R, const Interval& comps);

  ///
  static int preAllocatable()
  {
    return 2; // dyanmic allocatable.
  }

  BaseIVFAB<VolData>& getVolData()
  {
    return m_volData;
  }
  const BaseIVFAB<VolData>& getVolData() const
  {
    return m_volData;
  }

  ///
  void addFullIrregularVoFs(const IntVectSet& a_vofsToChange,
                            const EBGraph&    a_newGhostGraph,
                            const BaseIVFAB<VolData>&     a_newGhostData,
                            const EBGraph&    a_oldGhostGraph);

  ///multifluid angels dancing on the heads of pins.
  void addEmptyIrregularVoFs(const IntVectSet& a_vofsToChange,
                             const EBGraph&    a_newGraph);
  ///
  static void setVerbose(bool a_verbose);
  ///
  static void setVerboseDebug(bool a_verboseDebug);
private:

  static   bool s_verbose;
  static   bool s_verboseDebug;
  ///
  BaseIVFAB<VolData>   m_volData;

  ///
  BaseIFFAB<FaceData>  m_faceData[SpaceDim];

  ///
  bool m_isFaceDataDefined;

  ///
  bool m_isVoFDataDefined;

  void operator=(const EBDataImplem& ebiin)
  {;}

  EBDataImplem(const EBDataImplem& ebiin)
  {;}

  void
  coarsenFaceCentroid(RealVect&                a_centroidCoar,
                      const Vector<RealVect>&  a_centroidsFine,
                      const Vector<Real>&      a_areaFracFine,
                      const Vector<FaceIndex>& a_facesFine,
                      const FaceIndex&         a_faceCoar);
  void
  coarsenAreaFrac(Real& a_areaFracCoar,
                  const Vector<Real>& a_areaFracFine);

  void
  coarsenVolFracAndCentroid(Real&                   a_volFracCoar,
                            RealVect&               a_volCentroidCoar,
                            const Vector<Real>&     a_volFracFine,
                            const Vector<RealVect>& a_volCentroidFine,
                            const Vector<VolIndex>& a_fineVoFs,
                            const VolIndex&         a_coarVoF);

  void
  coarsenBoundaryAreaAndNormal(Real&                    a_bndryAreaCoar,
                               RealVect&                a_normalCoar,
                               const Vector<Real>&      a_bndryAreaFine,
                               const Vector<RealVect>&  a_normalFine);

  RealVect
  fineToCoarseTransform(const RealVect& a_finePoint,
                        const IntVect&  a_coarCell,
                        const IntVect&  a_fineCell);

  void
  coarsenBndryCentroid(RealVect&               a_bndryCentroidCoar,
                       const Vector<RealVect>& a_bndryCentroidFine,
                       const Vector<Real>&     a_bndryAreaFine,
                       const Vector<VolIndex>& a_fineVoFs,
                       const VolIndex&         a_coarVoF);

  void fetch(std::list<const VolData*>& fineVols, const Vector<VolIndex>& vofsFine) const;

  friend class EBIndexSpace;
  friend class EBISLevel;
};

///
/**
   Ref-counted version of EBDataImplem.
 */
class EBData
{
public:

  ///
  EBData();

  ///
  ~EBData();

  ///
  /**
     Noop.
  */
  void define(const Box& box, int comps);

  ///
  /**
     Noop
  */
  EBData(const Box& a_box, int a_comps);

  ///
  /**
     Copy the information from a_source to the over the intersection
     of the box a_region, the box of the current EBData and
     the box of a_source.  The Interval arguments are ignored.
     This function is required by LevelData.
  */
  void copy(const Box&      a_regionFrom,
            const Interval& a_Cd,
            const Box&      a_regionto,
            const EBData&   a_source,
            const Interval& a_Cs);

  /// each data holder is defined over the irregular cells of the graph
  void
  defineVoFData(const EBGraph& a_graph, const Box& a_region);

  /// each data holder is defined over the irregular cells of the graph
  void
  defineFaceData(const EBGraph& a_graph, const Box& a_region);

  ///define the whole thing
  void
  define(const EBGraph&           a_graph,
         const Vector<IrregNode>& a_irregData,
         const Box&               a_validBox);

  ///
  void
  coarsenVoFs(const EBData&   a_fineEBData,
              const EBGraph&  a_fineGraph,
              const EBGraph&  a_coarGraph,
              const Box&      a_validRegion);

  void
  coarsenFaces(const EBData&  a_fineEBData,
               const EBGraph& a_fineGraph,
               const EBGraph& a_coarGraph,
               const Box&     a_validRegion);
  ///
  const Real& volFrac(const VolIndex& a_vof) const;

  ///
  const Real& areaFrac(const FaceIndex& a_face1) const;

  ///
  const RealVect& centroid(const FaceIndex& facein) const;

  ///
  const RealVect& centroid(const VolIndex& a_vof) const;

  ///
  const RealVect& bndryCentroid(const VolIndex& a_vof) const;
  const RealVect& bndryCentroid(const VolIndex& a_vof, int face) const;

  ///
  const Real& bndryArea(const VolIndex& a_vof) const;
  const Real& bndryArea(const VolIndex& a_vof, int face) const;

  ///
  const RealVect& normal(const VolIndex& a_vof) const;
  const RealVect& normal(const VolIndex& a_vof, int face) const;

  /// used by multi-fluid code
  int facePhase(const VolIndex& a_vof, int face) const ;

  /// used by multi-fluid code
  const VolIndex& faceIndex(const VolIndex& a_vof, int face) const ;

  /// used by multi-fluid code
  void setFacePhase(const VolIndex& a_vof, int face, int phase);

  /// used by multi-fluid code
  void setFaceIndex(const VolIndex& a_vof, int face, int index);
  ///
  int numFacePhase(const VolIndex& a_vof) const ;

  void clearMultiBoundaries();
  void setBoundaryPhase(int phase);
  ///
  EBData(const EBData& a_ebiin);

  ///
  /**
     This is a pointer comparison.
  */
  bool operator==(const EBData& a_ebiin);

  ///
  EBData& operator=(const EBData& a_ebiin);

  ///
  /**
     This stuff required by LevelData in parallel:
  */
  int size(const Box& R, const Interval& comps) const;

  ///
  void linearOut(void* buf, const Box& R, const Interval& comps) const;

  ///
  void linearIn(void* buf, const Box& R, const Interval& comps);

  ///
  void addFullIrregularVoFs(const IntVectSet& a_vofsToChange,
                            const EBGraph&    a_newGhostGraph,
                            const BaseIVFAB<VolData>& a_grownData,
                            const EBGraph&    a_oldGhostGraph);

  ///multifluid angels dancing on the heads of pins.
  void addEmptyIrregularVoFs(const IntVectSet& a_vofsToChange,
                             const EBGraph&    a_newGraph);

  ///
  void
  computeNormalsAndBoundaryAreas(const EBGraph& a_graph,
                                 const Box&     a_validRegion);

  BaseIVFAB<VolData>& getVolData()
  {
    return m_implem->getVolData();
  }

  const BaseIVFAB<VolData>& getVolData() const
  {
    return m_implem->getVolData();
  }

  ///
  static int preAllocatable()
  {
    return 2; // dyanmic allocatable.
  }
private:

  ///
  RefCountedPtr<EBDataImplem> m_implem;

  friend class EBIndexSpace;
  friend class EBISLevel;
};

/*******************************/
inline const Real& EBData::volFrac(const VolIndex& a_vof) const
{
  return m_implem->volFrac(a_vof);
}
/*******************************/
inline const Real& EBData::bndryArea(const VolIndex& a_vof, int face) const
{
  return m_implem->bndryArea(a_vof, face);
}
inline const Real& EBData::bndryArea(const VolIndex& a_vof) const
{
  return m_implem->bndryArea(a_vof);
}
/*******************************/
inline const RealVect& EBData::normal(const VolIndex& a_vof, int face) const
{
  return m_implem->normal(a_vof, face);
}
inline const RealVect& EBData::normal(const VolIndex& a_vof) const
{
  return m_implem->normal(a_vof);
}
/*******************************/
inline const RealVect& EBData::centroid(const VolIndex& a_vof) const
{
  return m_implem->centroid(a_vof);
}
/*******************************/
inline const RealVect& EBData::bndryCentroid(const VolIndex& a_vof, int face) const
{
  return m_implem->bndryCentroid(a_vof, face);
}
inline const RealVect& EBData::bndryCentroid(const VolIndex& a_vof) const
{
  return m_implem->bndryCentroid(a_vof);
}
/*******************************/
inline const RealVect& EBData::centroid(const FaceIndex& a_face) const
{
  return m_implem->centroid(a_face);
}
/*******************************/
inline const Real& EBData::areaFrac(const FaceIndex& a_face) const
{
  return m_implem->areaFrac(a_face);
}
/*******************************/

#include "NamespaceFooter.H"
#endif
