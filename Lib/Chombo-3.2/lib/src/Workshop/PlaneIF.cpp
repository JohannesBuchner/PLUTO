#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PlaneIF.H"

#include "NamespaceHeader.H"

PlaneIF::PlaneIF(const RealVect& a_normal,
                 const RealVect& a_point,
                 const bool&     a_inside)
{
  // Remember the parameters
  m_normal = a_normal;
  m_point = a_point;
  m_inside = a_inside;
}

PlaneIF::PlaneIF(const PlaneIF& a_inputIF)
{
  // Remember the parameters
  m_normal = a_inputIF.m_normal;
  m_point = a_inputIF.m_point;
  m_inside = a_inputIF.m_inside;
}

PlaneIF::~PlaneIF()
{
}

void PlaneIF::GetParams(RealVect& a_normal,
                        RealVect& a_point,
                        bool&     a_inside) const
{
  // Copy parameter information over
  a_normal = m_normal;
  a_point = m_point;
  a_inside = m_inside;
}

void PlaneIF::SetParams(const RealVect& a_normal,
                        const RealVect& a_point,
                        const bool&     a_inside)
{
  // Set parameter information
  m_normal = a_normal;
  m_point = a_point;
  m_inside = a_inside;
}

Real PlaneIF::value(const RealVect& a_point) const
{
  Real retval;

  // Vector from a_point to a point on the plane
  RealVect direction(m_point);
  direction -= a_point;

  // Compute (and return) the dot product of the vector above and the normal
  retval = direction.dotProduct(m_normal);


  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

GeometryService::InOut PlaneIF::InsideOutside(const RealVect& lo, const RealVect& hi) const
{

  RealVect len = hi-lo;
  Real firstValue = value(lo);
  GeometryService::InOut rtn;
  if ( firstValue < 0 )
    {
      rtn = GeometryService::Regular;
    }
  else
    {
      rtn = GeometryService::Covered;
    }
  Box unit(IntVect::Zero, IntVect::Unit);
  for (BoxIterator b(unit); b.ok(); ++b)
    {
      RealVect corner = lo + RealVect(b())*len;
      Real functionValue = value(corner);

      if (functionValue * firstValue <= 0.0 )
      {
        return GeometryService::Irregular;
      }
    }
  return rtn;
}


BaseIF* PlaneIF::newImplicitFunction() const
{
  PlaneIF* planePtr = new PlaneIF(m_normal,
                                  m_point,
                                  m_inside);

  return static_cast<BaseIF*>(planePtr);
}

#include "NamespaceFooter.H"
