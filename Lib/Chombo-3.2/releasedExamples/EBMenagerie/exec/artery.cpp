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
    const char* in_file = "artery.inputs";

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
    const char* basename = "artery";

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

// A specialized geometry class to create a few examples of something which
// looks like a segment of a artery.  For this example the domain is expected
// to be [0,1]^SpaceDim.
//
// If "type" is 1, the geometry is a surface of revolution (revolved about
// the axis defined by (x,0.5,0.5)).  Thus, the centerline is a line in the x
// direction and radius varies along that centerline.
//
// If "type" is 2, the geometry has a varying centerline and radius.
//
// In both cases, the artery starts at the plane x = 0 and ends at the plane
// x = 1.
class ArteryIF: public BaseIF
{
public:
  // "a_type" it the artery type (see above)
  // "a_radius" is the initial radius (at the plane x = 0)
  // "a_amplitude" is the amplitude of the perturbation (zero -> cylinder)
  // "a_minX" is the minimum value of x in the domain
  // "a_maxX" is the maximum value of x in the domain
  // "a_inside" being true means the domain is inside the artery, otherwise it
  //            is outside
  ArteryIF(const int&  a_type,
           const Real& a_radius,
           const Real& a_amplitude,
           const Real& a_minX,
           const Real& a_maxX,
           const bool& a_inside);

  ArteryIF(const ArteryIF& a_inputIF);

  virtual ~ArteryIF();

  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  // The artery type (see above)
  int  m_type;

  // The initial radius (at the plane x = 0)
  Real m_radius;

  // The amplitude of the perturbation (zero -> cylinder)
  Real m_amplitude;

  // The length of the unperturbed inlet tube
  Real m_inlet;

  // The minimum value of x in the domain
  Real m_minX;

  // The maximum value of x in the domain
  Real m_maxX;

  // True means the domain is inside the artery, otherwise it is outside
  bool m_inside;

private:
  ArteryIF()
  {
    MayDay::Error("ArteryIF uses strong construction");
  }

  void operator=(const ArteryIF& a_inputIF)
  {
    MayDay::Error("ArteryIF doesn't allow assignment");
  }
};

ArteryIF::ArteryIF(const int&  a_type,
                   const Real& a_radius,
                   const Real& a_amplitude,
                   const Real& a_minX,
                   const Real& a_maxX,
                   const bool& a_inside)
{
  if (a_type != 1 && a_type != 2)
  {
    MayDay::Error("ArteryIF::ArteryIF: Artery type must be 1 or 2");
  }

  // Remember the parameters
  m_type      = a_type;
  m_radius    = a_radius;
  m_amplitude = a_amplitude;
  m_minX      = a_minX;
  m_maxX      = a_maxX;
  m_inside    = a_inside;
}

ArteryIF::ArteryIF(const ArteryIF& a_inputIF)
{
  // Remember the parameters
  m_type      = a_inputIF.m_type;
  m_radius    = a_inputIF.m_radius;
  m_amplitude = a_inputIF.m_amplitude;
  m_minX      = a_inputIF.m_minX;
  m_maxX      = a_inputIF.m_maxX;
  m_inside    = a_inputIF.m_inside;
}

ArteryIF::~ArteryIF()
{
}

Real ArteryIF::value(const RealVect& a_point) const
{
  Real retval;

  Real x = a_point[0];

  if ((x - m_minX) < ((m_maxX - m_minX) - 1.0)/2.0)
  {
    x = 0.0;
  }
  else if ((x - m_minX) < ((m_maxX - m_minX) + 1.0)/2.0)
  {
    x = (x - m_minX) - ((m_maxX - m_minX) - 1.0)/2.0;
  }
  else
  {
    x = 1.0;
  }

  Real func = 0.5 * (cos(2*M_PI*    x)
                   - cos(2*M_PI*3.5*x));

  Real perturb = m_amplitude * func;

  if (m_type == 1)
  {
    Real radius = 0.0;

    for (int idir = 1; idir < SpaceDim; idir++)
    {
      radius += (a_point[idir]-0.5)*(a_point[idir]-0.5);
    }

    radius = sqrt(radius);

    retval = radius - (m_radius + perturb);
  }

  if (m_type == 2)
  {
    Real radius = 0.0;

    if (SpaceDim == 2)
    {
      radius += ((a_point[1]-0.5) - perturb*sin(2*M_PI*x))
              * ((a_point[1]-0.5) - perturb*sin(2*M_PI*x));
    }

    if (SpaceDim == 3)
    {
      radius += ((a_point[1]-0.5) - perturb*cos(2*M_PI*x))
              * ((a_point[1]-0.5) - perturb*cos(2*M_PI*x));

      radius += ((a_point[2]-0.5) - perturb*sin(2*M_PI*x))
              * ((a_point[2]-0.5) - perturb*sin(2*M_PI*x));
    }

    radius = sqrt(radius);

    retval = radius - (m_radius + perturb);
  }

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

BaseIF* ArteryIF::newImplicitFunction() const
{
  ArteryIF* arteryPtr = new ArteryIF(m_type,
                                     m_radius,
                                     m_amplitude,
                                     m_minX,
                                     m_maxX,
                                     m_inside);

  return static_cast<BaseIF*>(arteryPtr);
}

BaseIF* makeGeometry(Box&      a_domain,
                     RealVect& a_origin,
                     Real&     a_dx)
{
  RealVect point;
  RealVect normal;
  bool     insideRegular;

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

  int arteryType;
  pp.get("artery_type",arteryType);

  Real radius;
  pp.get("radius",radius);

  Real amplitude;
  pp.get("amplitude",amplitude);

  // Parm Parse doesn't get bools, so work-around with int
  int intInsideRegular;
  pp.get("insideRegular",intInsideRegular);

  if (intInsideRegular != 0) insideRegular = true;
  if (intInsideRegular == 0) insideRegular = false;

  ArteryIF implicit(arteryType,radius,amplitude,prob_lo[0],prob_hi,insideRegular);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(implicit,0,vectDx);

  // This generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx,workshop);

  return implicit.newImplicitFunction();
}
