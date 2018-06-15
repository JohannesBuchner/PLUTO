#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SPACE.H"

#include "MitochondriaIF.H"

#include "UsingNamespace.H"

Real generalDimpleValue(const RealVect & a_point,
                        const bool     & a_cylinder)
{
  Real x;
  Real y;
  Real z;

  x = D_TERM(a_point[0], +        0.0, +        0.0);
  y = D_TERM(       0.0, + a_point[1], +        0.0);
  z = D_TERM(       0.0, +        0.0, + a_point[2]);

  Real const1 = 0.540901699437495;
  Real r1     = 0.031250886819824;
  Real r2     = 0.1;

  Real part1 = Max(
                   Abs(x - const1) - r2,
                   D_TERM((Real)0.0, + y*y, + z*z) - r1*r1
                  );

  Real part2 = Min(
                   part1,
                   D_TERM((x-(Real)0.3)*(x-(Real)0.3), + y*y, + z*z) - r2*r2
                  );

  Real const2 = 0.3809017;
  Real const3 = 0.4409017;

  Real part3 = Max(
                   -(x - const2),
                    (x - const3)
                  );

  Real const4 =  0.058778;
  Real const5 =  1.376381666666666666666666666;
  Real const6 =  0.0825829;
  Real const7 =  0.0275276;
  Real const8 = 16.66666666666666666666666666;

  Real part4 = sqrt(D_TERM(0.0, + y*y, + z*z)) -
               (
                const4 -
                const5 * pow(         (x - const2),1) +
                const6 * pow(const8 * (x - const2),2) -
                const7 * pow(const8 * (x - const2),3)
               );

  Real part5 = Max(part3,part4);

  Real part6 = Min(part2,part5);

  Real const9 = 0.640901699437495;

  Real part7 = -(x - const9);

  Real const10 =  1.010079120233248;
  Real const11 = 10.038444719705328;

  Real const12 = 0.710643486023442;
  Real const13 = 0.703552297823303;

  Real const14 = const1 + 0.11;
  Real rcut = 0.16;

  if (a_cylinder)
  {
    const14 += 0.0061 * (y/rcut)*(y/rcut);
  }

  Real part8 = const10 *
               (
                (const13 * (x - const14)) +
                (const12 * (sqrt(D_TERM(0.0, + y*y, + z*z)) - r1))
               ) *
               (
                1.0 - (
                       const11 * pow(
                                     const13 * (x - const14) +
                                     const12 * (sqrt(D_TERM(0.0, + y*y, + z*z)) - r1),
                                     2
                                    )
                      )
               ) -
               (
                (const12 * (x - const14)) -
                (const13 * (sqrt(D_TERM(0.0, + y*y, + z*z)) - r1))
               );

  Real part9 = Max(part7,part8);

  Real part10 = Min(part6,part9);

  Real r4 = 0.80;
  Real xoff = 0.20;

  Real part11 = D_TERM((x-xoff)*(x-xoff), + y*y, + z*z) - r4*r4;

  Real part12 = Max(part10,part11);

  Real cylinder = D_TERM(0.0, + y*y, + z*z) - rcut*rcut;

  Real cut = Max(part12,cylinder);

  Real final = cut;

  return final;
}

SphereDimpleIF::SphereDimpleIF()
{
}

SphereDimpleIF::~SphereDimpleIF()
{
}

Real SphereDimpleIF::value(const RealVect & a_point) const
{
  bool cylinderDimple = false;

  return generalDimpleValue(a_point,cylinderDimple);
}

BaseIF* SphereDimpleIF::newImplicitFunction() const
{
  return new SphereDimpleIF();
}


CylinderDimpleIF::CylinderDimpleIF()
{
}

CylinderDimpleIF::~CylinderDimpleIF()
{
}

Real CylinderDimpleIF::value(const RealVect & a_point) const
{
  bool cylinderDimple = true;

  return generalDimpleValue(a_point,cylinderDimple);
}

BaseIF* CylinderDimpleIF::newImplicitFunction() const
{
  return new CylinderDimpleIF();
}


MitochondriaIF::MitochondriaIF()
{
  m_dimple = new SphereDimpleIF();
}

MitochondriaIF::~MitochondriaIF()
{
  delete m_dimple;
}

Real MitochondriaIF::value(const RealVect & a_point) const
{
  Real x;
  Real y;
  Real z;

  x = D_TERM(a_point[0], +        0.0, +        0.0);
  y = D_TERM(       0.0, + a_point[1], +        0.0);
  z = D_TERM(       0.0, +        0.0, + a_point[2]);

  Real r3 = 0.9;

  Real sphere = D_TERM(x*x, + y*y, + z*z) - r3*r3;
  // Real sphere = D_TERM(x*x, + 0.0, + z*z) - r3*r3;

  Real cut1 = Max(sphere,-m_dimple->value(a_point));
#if 0
  Real paste1 = Min(cut1,m_dimple2->value(a_point));

  Real cut2 = Max(paste1,-m_dimple3->value(a_point));

  Real final = cut2;
#else
  Real final = cut1;
#endif

  return final;
}

BaseIF* MitochondriaIF::newImplicitFunction() const
{
  return new MitochondriaIF();
}
