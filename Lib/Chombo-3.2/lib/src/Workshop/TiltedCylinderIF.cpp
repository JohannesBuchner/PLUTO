#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PolyGeom.H"
#include "TiltedCylinderIF.H"

#include "NamespaceHeader.H"

TiltedCylinderIF::TiltedCylinderIF(const Real&     a_radius,
                                   const RealVect& a_direction,
                                   const RealVect& a_point,
                                   const bool&     a_inside)
{
  SetParams(a_radius,a_direction,a_point,a_inside);
}

TiltedCylinderIF::TiltedCylinderIF(const TiltedCylinderIF& a_inputIF)
{
  // Remember the parameters
  m_radius = a_inputIF.m_radius;
  m_direction = a_inputIF.m_direction;
  m_tiltedAxis = a_inputIF.m_tiltedAxis;
  m_coordDir = a_inputIF.m_coordDir;
  m_point = a_inputIF.m_point;
  m_inside = a_inputIF.m_inside;

  // Precompute the radius squared
  m_radius2 = m_radius * m_radius;

  // Precompute the length squared of the direction vector
  m_length2 = PolyGeom::dot(m_direction,m_direction);
}

TiltedCylinderIF::~TiltedCylinderIF()
{
}

GeometryService::InOut TiltedCylinderIF::InsideOutside(const RealVect& lo,
                                                       const RealVect& hi) const
{

  // Fast Sphere-Box intersection from "Graphics Gems" pp 335-338, 1990

  Real dmin = 0; //distance from cylinder axis to nearest  point in box (squared)
  Real dmax = 0; //distance from cylinder axis to farthest point in box (squared)
  Real ai, bi, a, b;

  Tuple<int,CH_SPACEDIM-1> tanDirs = PolyGeom::computeTanDirs(m_coordDir);
  for (int i=0; i<CH_SPACEDIM-1; ++i)
    {
      ai = m_point[tanDirs[i]] - lo[tanDirs[i]];
      bi = m_point[tanDirs[i]] - hi[tanDirs[i]];
      a = ai * ai;
      b = bi * bi;
      dmax = dmax + Max(a,b);
      if (m_point[tanDirs[i]] < lo[tanDirs[i]] || m_point[tanDirs[i]] > hi[tanDirs[i]]) dmin = dmin + Min(a,b);
    }

  if (m_inside)
    {
      if (dmin >= m_radius2) return GeometryService::Covered;
      if (dmax <  m_radius2) return GeometryService::Regular;
    }
  else
    {
      if (dmin >  m_radius2) return GeometryService::Regular;
      if (dmax <= m_radius2) return GeometryService::Covered;
    }
  return GeometryService::Irregular;
}


void TiltedCylinderIF::GetParams(Real&     a_radius,
                                 RealVect& a_direction,
                                 RealVect& a_point,
                                 bool&     a_inside) const
{
  // Copy parameter information over
  a_radius = m_radius;
  a_direction = m_direction;
  a_point = m_point;
  a_inside = m_inside;
}

void TiltedCylinderIF::SetParams(const Real&     a_radius,
                                 const RealVect& a_direction,
                                 const RealVect& a_point,
                                 const bool&     a_inside)
{
  // Set parameter information
  m_radius = a_radius;
  m_direction = a_direction;
  m_point = a_point;
  m_inside = a_inside;

  // Precompute the radius squared
  m_radius2 = m_radius * m_radius;

  // Precompute the length squared of the direction vector
  m_length2 = PolyGeom::dot(m_direction,m_direction);

  m_tiltedAxis = true;

  for (int i=0; i<CH_SPACEDIM; ++i)
    {
      if (m_direction[i] != 0.0)
        {
          bool idDirection = true;
          Tuple<int,CH_SPACEDIM-1> tanDirs = PolyGeom::computeTanDirs(i);
          for (int itan = 0; itan < CH_SPACEDIM-1; itan++)
            {
              if (m_direction[tanDirs[itan]] != 0.0)
                {
                  idDirection = false;
                }
            }

          if (idDirection)
            {
              m_tiltedAxis = false;
              m_coordDir = i;
            }
        }
    }
}

Real TiltedCylinderIF::value(const RealVect& a_point) const
{
  Real retval;

  // Get a_point relative to a point on the cylinder's axis, m_point
  RealVect delta;

  delta = a_point;
  delta -= m_point;

  // Get the dot product of the relative vector (above) and the cylinder's
  // axis vector and normalize by the square of the length of the direction
  // vector
  Real dot = PolyGeom::dot(m_direction,delta);
  dot /= m_length2;

  // Find the vector from a_point to the point on the cylinder's axis closest
  // to a_point
  RealVect closest;

  closest = m_direction;
  closest *= dot;
  closest -= delta;

  // Find the length squared of this vector and compare it to the radius
  // squared
  Real length2 = PolyGeom::dot(closest,closest);

  // Return the difference between the sqaures (zero on the cylinder)
  retval = length2 - m_radius2;

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

BaseIF* TiltedCylinderIF::newImplicitFunction() const
{
  TiltedCylinderIF* tiltedCylinderPtr = new TiltedCylinderIF(m_radius,
                                                             m_direction,
                                                             m_point,
                                                             m_inside);

  return static_cast<BaseIF*>(tiltedCylinderPtr);
}

#include "NamespaceFooter.H"
