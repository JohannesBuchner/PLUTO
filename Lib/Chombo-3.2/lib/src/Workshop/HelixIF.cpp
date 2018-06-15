#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "HelixIF.H"

#include "NamespaceHeader.H"

HelixIF::HelixIF(const BaseIF& a_impFunc,
                 const Real&   a_rate,
                 const bool&   a_inside,
                 const bool&   a_vertical)
{
  // Make a copy of the implicit function and rate
  m_impFunc = a_impFunc.newImplicitFunction();
  m_rate    = a_rate;

  // Save inside flag
  m_inside = a_inside;

  // Treat the impliction function are vertical (in the x-z plane) and not
  // horizontal (in the x-y plane)
  m_vertical = a_vertical;

  // This is the valid range of z for sweeping the vertical IF
  if (m_vertical)
  {
    m_minz = -M_PI / m_rate;
    m_maxz =  M_PI / m_rate;
  }
}

HelixIF::HelixIF(const HelixIF& a_inputIF)
{
  // Make a copy of the implicit function and rate
  m_impFunc = a_inputIF.m_impFunc->newImplicitFunction();
  m_rate    = a_inputIF.m_rate;

  // Save inside flag
  m_inside = a_inputIF.m_inside;

  // Treat the impliction function are vertical (in the x-z plane) and not
  // horizontal (in the x-y plane)
  m_vertical = a_inputIF.m_vertical;

  // This is the valid range of z for sweeping the vertical IF
  m_minz = a_inputIF.m_minz;
  m_maxz = a_inputIF.m_maxz;
}

HelixIF::~HelixIF()
{
  delete m_impFunc;
}

Real HelixIF::value(const RealVect& a_point) const
{
  Real retval;

#if CH_SPACEDIM == 2
  // In 2D, return the given implicit function
  retval =  m_impFunc->value(a_point);
#elif CH_SPACEDIM == 3
  // Get the point coordinates
  Real x = a_point[0];
  Real y = a_point[1];
  Real z = a_point[2];

  if (!m_vertical)
  {
    // Compute the rotation angle
    Real theta = m_rate * z;

    // Rotate the x,y point back to z = 0
    Real x2,y2;
    x2 =  cos(-theta)*x + sin(-theta)*y;
    y2 = -sin(-theta)*x + cos(-theta)*y;

    // Construct this point in x-y plane (z = 0)
    RealVect point(x2,y2,0.0);

    // Return the given implicit function at that point
    retval = m_impFunc->value(point);
  }
  else
  {
    Real r = sqrt(x*x + y*y);
    Real theta = atan2(y,x);

    Real unrotateZ = z - theta/m_rate;
    Real numZ = (unrotateZ - m_minz) / (m_maxz - m_minz);
    Real pullbackZ = (numZ - floor(numZ)) * (m_maxz - m_minz) + m_minz;

    RealVect point(r,pullbackZ,0.0);

    retval = m_impFunc->value(point);
  }
#else
  MayDay::Abort("need higher dim in HelixIF\n");
#endif

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

BaseIF* HelixIF::newImplicitFunction() const
{
  HelixIF* helixPtr;

  helixPtr = new HelixIF(*m_impFunc,m_rate,m_inside,m_vertical);

  return static_cast<BaseIF*>(helixPtr);
}

#include "NamespaceFooter.H"
