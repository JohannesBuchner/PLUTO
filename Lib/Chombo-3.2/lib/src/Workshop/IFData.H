#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _IFDATA_H_
#define _IFDATA_H_

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
#include "IFSlicer.H"
#include "CoordinateSystem.H"

using std::endl;
using std::cerr;
#include "NamespaceHeader.H"

template <int dim> class IFData
{
public:
  // data types to record vertex signs and edge intersections
  typedef IndexTM<int,dim>                      EdgeIndex;
  typedef map<EdgeIndex,Real,LexLT<EdgeIndex> > EdgeIntersections;

  typedef IndexTM<int,dim>                      Vertex;
  typedef map<Vertex,int,LexLT<Vertex> >        CornerSigns;

  typedef IndexTM<int,dim>                      IvDim;
  typedef IndexTM<Real,dim>                     RvDim;

  typedef map<IvDim,RvDim,LexLT<IvDim> >        NormalDerivatives;

  // member data
  CornerSigns           m_cornerSigns;
  EdgeIntersections     m_intersections;
  IFSlicer<dim>*        m_function;

  CoordinateSystem<dim> m_globalCoord;
  CoordinateSystem<dim> m_cellCenterCoord;
  CoordinateSystem<dim> m_parentCoord;
  CoordinateSystem<dim> m_localCoord;

  int                   m_maxOrder;
  NormalDerivatives     m_normalDerivatives;
  bool                  m_badNormal;

  bool                  m_allVerticesIn;
  bool                  m_allVerticesOut;
  bool                  m_allVerticesOn;

  // empty constructor
  IFData();

  // copy constructor
  IFData(const IFData<dim>& a_IFData);

  // makes an initial IFData from a BaseIF
  IFData(const BaseIF & a_function,
         const RvDim  & a_dx,
         const RvDim  & a_origin,
         const int    & a_maxOrder);

  // makes an lower dimensional IFData from a IFSlicer; used for refinement
  IFData(IFSlicer<dim> * a_function,
         const RvDim   & a_dx,
         const RvDim   & a_origin,
         const int     & a_maxOrder);

  // makes a reduced IFData
  IFData(const IFData<dim+1> & a_hInfo,
#if RECURSIVE_GEOMETRY_GENERATION != 0
         const int           & a_maxOrder,
#endif
         const int           & a_idir,
         const int           & a_hilo);

  // destructor
  ~IFData();

  void setNormalDerivatives();
  void makeCornerSigns();
  void findIntersectionPts();
  void defineLocalCoords();

  // used by the constructor
  bool isConnected(int          & a_edgeDir,
                   const Vertex & a_vertex1,
                   const Vertex & a_vertex2);

  void  makeEdgeKey(const int    & a_edgeDir,
                    const Vertex & a_vertex1,
                    const Vertex & a_vertex2);

  //this function calls the midPtRootFinder
  Real rootFinder(const EdgeIndex & a_thisEdge);

  Real BrentRootFinder(const RvDim & a_loPt,
                       const RvDim & a_hiPt,
                       const int   & a_edgeDir) const;
//                       const Real  & a_smallestRoot,
//                       const Real  & a_biggestRoot) const;

  void checkIntersection(bool       & a_hiOn,
                         bool       & a_loOn,
                         const Real & a_pt) const;

  void remakeCornerSigns();

  // output functions
  void print(ostream& out) const;

  // operators
  void operator=(const IFData & a_ifData);

  // this function sets the changing direction given the normal direction
  void setChangingDirection();

  bool searchNormalDir(const int i) const;
};

template<> class IFData<1>
{
public:
  // true is in the fluid
  typedef Real            EdgeIntersection;

  typedef int             Vertex;
  typedef map<Vertex,int> CornerSigns;

  typedef int             IvDim;
  typedef Real            RvDim;

  // member data
  CornerSigns         m_cornerSigns;
  EdgeIntersection    m_intersection;

  CoordinateSystem<1> m_globalCoord;
  CoordinateSystem<1> m_cellCenterCoord;
  CoordinateSystem<1> m_parentCoord;

  bool                m_allVerticesIn;
  bool                m_allVerticesOut;
  bool                m_allVerticesOn;
  bool                m_badNormal;

  IFData();

  IFData(const IFData<1> & a_ifData);

  IFData(const IFData<2> & a_hIFData,
#if RECURSIVE_GEOMETRY_GENERATION != 0
         const int       & a_maxOrder,
#endif
         const int       & a_idir,
         const int       & a_hilo);

  ~IFData();

  void print(ostream& a_out) const;

  void operator=(const IFData<1> & a_ifData);
};

#include "NamespaceFooter.H"

#include "IFDataImplem.H"

#endif
