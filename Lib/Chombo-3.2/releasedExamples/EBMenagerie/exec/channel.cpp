#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <stdio.h>

#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "DebugDump.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "PolyGeom.H"
#include "GeometryShop.H"

#include "EBFABView.H"

#include "BaseIF.H"

#include "EBMenagerieUtils.H"

#include "UsingNamespace.H"

/***************/
// Define an EBIS.
/***************/
BaseIF* makeGeometry(Box&      a_domain,
                     RealVect& a_origin,
                     Real&     a_dx);

int main (int argc, char** argv)
{
#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Begin forever present scoping trick
  {
    const char* in_file = "channel.inputs";

    if (argc > 1)
    {
      in_file = argv[1];
    }

    // Parse input file
    ParmParse pp(0,NULL,NULL,in_file);

    Box domain;
    RealVect origin;
    Real dx;

    BaseIF* implicit;

    // Make geometry
    implicit = makeGeometry(domain,origin,dx);

    // Make grids
    DisjointBoxLayout grids;
    makeLayout(grids,domain);

    // Create ebislayout
    int nghost = 0;
    EBISLayout ebisl;
    makeEBISL(ebisl,grids,domain,nghost);

    // Make a LevelData
    int nComps = 1;

    IntVect ghost = IntVect::Unit;
    ghost *= nghost;

    RefCountedPtr<DataFactory<EBCellFAB> > rcpFactory(new EBCellFactory(ebisl));
    LevelData<EBCellFAB> level(grids,nComps,ghost,*rcpFactory);

    // Put some data in the data holders
    fillData(level,origin,dx,*implicit);

    // Done with this object
    delete implicit;

    // Write the data and the EB out
    const char* basename = "channel";

    char name[1000];
    sprintf(name,"%s%dd.hdf5",basename,SpaceDim);
#ifdef CH_USE_HDF5
    writeEBLevelname(&level,name);
#endif
  } // End scoping trick

  // Clean up
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();

  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}

class ChannelIF: public BaseIF
{
public:
  // A contracting and then expanding channel
  ChannelIF(const Real& a_length1,
            const Real& a_length2,
            const Real& a_transitionLength,
            const Real& a_largeWidth,
            const Real& a_smallWidth,
            const Real& a_depth);

  ChannelIF(const ChannelIF& a_inputIF);

  virtual ~ChannelIF();

  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  // Distance from the inlet to the first inflection point
  Real m_length1;

  // Distance from the first inflection point to the second inflection point
  Real m_length2;

  // The length of the transition from the large width to the small width
  Real m_transitionLength;

  // Width at the inlet and outlet
  Real m_largeWidth;

  // Width at the throat
  Real m_smallWidth;

  // In 3D, the depth of the channel (centered around 0.0)
  Real m_depth;

private:
  ChannelIF()
  {
    MayDay::Error("ChannelIF uses strong construction");
  }

  void operator=(const ChannelIF& a_inputIF)
  {
    MayDay::Error("ChannelIF doesn't allow assignment");
  }
};

ChannelIF::ChannelIF(const Real& a_length1,
                     const Real& a_length2,
                     const Real& a_transitionLength,
                     const Real& a_largeWidth,
                     const Real& a_smallWidth,
                     const Real& a_depth)
{
  // Remember the parameters
  m_length1          = a_length1;
  m_length2          = a_length2;
  m_transitionLength = a_transitionLength;
  m_largeWidth       = a_largeWidth;
  m_smallWidth       = a_smallWidth;
  m_depth            = a_depth;
}

ChannelIF::ChannelIF(const ChannelIF& a_inputIF)
{
  // Remember the parameters
  m_length1          = a_inputIF.m_length1;
  m_length2          = a_inputIF.m_length2;
  m_transitionLength = a_inputIF.m_transitionLength;
  m_largeWidth       = a_inputIF.m_largeWidth;
  m_smallWidth       = a_inputIF.m_smallWidth;
  m_depth            = a_inputIF.m_depth;
}

ChannelIF::~ChannelIF()
{
}

Real ChannelIF::value(const RealVect& a_point) const
{
  Real retval;

  Real x = a_point[0];

  // Symmetric in y
  Real y = Abs(a_point[1]);

  Real z;
  if (SpaceDim > 2)
  {
    // Symmetric in z
    z = Abs(a_point[2]);
  }
  else
  {
    z = 0.0;
  }

  if (x < m_length1 - m_transitionLength/2.0)
  {
    retval = y - m_largeWidth/2.0;
  }
  else if (x < m_length1 + m_transitionLength/2.0)
  {
    Real localx = x - m_length1;
    Real scaling = (m_largeWidth/2.0 - m_smallWidth/2.0) / 2.0;
    Real offset  = (m_largeWidth/2.0 + m_smallWidth/2.0) / 2.0;

    retval = y - (-scaling * sin(M_PI * (localx / m_transitionLength)) + offset);
  }
  else if (x < m_length1 + m_length2 - m_transitionLength/2.0)
  {
    retval = y - m_smallWidth/2.0;
  }
  else if (x < m_length1 + m_length2 + m_transitionLength/2.0)
  {
    Real localx = x - (m_length1 + m_length2);
    Real scaling = (m_largeWidth/2.0 - m_smallWidth/2.0) / 2.0;
    Real offset  = (m_largeWidth/2.0 + m_smallWidth/2.0) / 2.0;

    retval = y - (scaling * sin(M_PI * (localx / m_transitionLength)) + offset);
  }
  else
  {
    retval = y - m_largeWidth/2.0;
  }

  if (SpaceDim > 2)
  {
    retval = Max(retval,z - m_depth/(Real)2.0);
  }

  return retval;
}

BaseIF* ChannelIF::newImplicitFunction() const
{
  ChannelIF* ChannelPtr = new ChannelIF(*this);

  return static_cast<BaseIF*>(ChannelPtr);
}

BaseIF* makeGeometry(Box&      a_domain,
                     RealVect& a_origin,
                     Real&     a_dx)
{
  RealVect point;
  RealVect normal;

  // parse input file
  ParmParse pp;

  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

  CH_assert(n_cell.size() == SpaceDim);

  IntVect lo = IntVect::Zero;
  IntVect hi;

  for (int ivec = 0; ivec < SpaceDim; ivec++)
  {
    if (n_cell[ivec] <= 0)
    {
      pout() << "Bogus number of cells input = " << n_cell[ivec];
      exit(1);
    }

    hi[ivec] = n_cell[ivec] - 1;
  }

  a_domain.setSmall(lo);
  a_domain.setBig(hi);

  Vector<Real> prob_lo(SpaceDim,1.0);
  Real prob_hi;

  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);

  a_dx = (prob_hi-prob_lo[0])/n_cell[0];

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    a_origin[idir] = prob_lo[idir];
  }

  Real length1;
  pp.get("length1",length1);

  Real length2;
  pp.get("length2",length2);

  Real transitionLength;
  pp.get("transition_length",transitionLength);

  Real largeWidth;
  pp.get("large_width",largeWidth);

  Real smallWidth;
  pp.get("small_width",smallWidth);

  Real depth;
  pp.get("depth",depth);

  ChannelIF implicit(length1,
                     length2,
                     transitionLength,
                     largeWidth,
                     smallWidth,
                     depth);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(implicit,0,vectDx);

  // This generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx,workshop);

  return implicit.newImplicitFunction();
}
