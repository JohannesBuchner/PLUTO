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
    const char* in_file = "spiral.inputs";

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
    const char* basename = "spiral";

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

class SpiralIF: public BaseIF
{
public:
  // "a_number" the numberr of  spirals
  // "a_inside" choose between the two spirals
  SpiralIF(const int&  a_number,
           const Real& a_variation,
           const bool& a_inside);

  SpiralIF(const SpiralIF& a_inputIF);

  virtual ~SpiralIF();

  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  // The number of spirals
  int  m_number;

  // Amplitude of variation
  Real m_variation;

  // True means the domain is inside the spiral, otherwise it is outside
  bool m_inside;

private:
  SpiralIF()
  {
    MayDay::Error("SpiralIF uses strong construction");
  }

  void operator=(const SpiralIF& a_inputIF)
  {
    MayDay::Error("SpiralIF doesn't allow assignment");
  }
};

SpiralIF::SpiralIF(const int&  a_number,
                   const Real& a_variation,
                   const bool& a_inside)
{
  // Remember the parameters
  m_number    = a_number;
  m_variation = a_variation;
  m_inside    = a_inside;
}

SpiralIF::SpiralIF(const SpiralIF& a_inputIF)
{
  // Remember the parameters
  m_number    = a_inputIF.m_number;
  m_variation = a_inputIF.m_variation;
  m_inside    = a_inputIF.m_inside;
}

SpiralIF::~SpiralIF()
{
}

Real SpiralIF::value(const RealVect& a_point) const
{
  Real retval;

  Real x = a_point[0];
  Real y = a_point[1];

  // Get the polar coordinates radius
  Real r = sqrt(x*x + y*y);

  if (r == 0.0)
  {
    retval = 0.0;
  }
  else
  {
    // Get the polar coordinates angle and make between 0 and 2*pi
    Real theta = atan2(y,x);
    if (theta < 0)
    {
      theta += 2*M_PI;
    }

    retval = sin(-m_number * 2*M_PI*r + theta);

    Real frac = theta/(2*M_PI);
    Real steepness = 32.0;
    if (m_number*r - frac > m_number - 1.25)
    {
      if (m_number*r - frac <= m_number - 0.25 && frac <= 1.0/steepness)
      {
        retval = (1.0 - steepness*frac)*retval + steepness*frac*1.0;
      }
      else
      {
        retval = 1.0;
      }
    }
  }

  retval *= m_variation / 2.0;

  if (SpaceDim > 2)
  {
    retval -= a_point[2];
  }

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

BaseIF* SpiralIF::newImplicitFunction() const
{
  SpiralIF* spiralPtr = new SpiralIF(*this);

  return static_cast<BaseIF*>(spiralPtr);
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

  int number;
  pp.get("number",number);

  Real variation;
  pp.get("variation",variation);

  // Parm Parse doesn't get bools, so work-around with int
  int intInsideRegular;
  pp.get("insideRegular",intInsideRegular);

  if (intInsideRegular != 0) insideRegular = true;
  if (intInsideRegular == 0) insideRegular = false;

  SpiralIF implicit(number,variation,insideRegular);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(implicit,0,vectDx);

  // This generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx,workshop);

  return implicit.newImplicitFunction();
}
