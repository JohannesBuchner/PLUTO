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

#include "EBFaceFAB.H"
#include "EBArithF_F.H"
#include "NamespaceHeader.H"
/**********************/
/**********************/
EBFaceFAB::EBFaceFAB():BaseEBFaceFAB<Real>()
{
}

/**********************/
/**********************/
EBFaceFAB::EBFaceFAB(const EBISBox& a_ebisBox,
                     const Box& a_region,
                     int a_iDir, int a_nComp)
  :BaseEBFaceFAB<Real>(a_ebisBox, a_region, a_iDir, a_nComp)
{

#ifdef CH_USE_SETVAL
  setVal(BaseFabRealSetVal);
#else
  setVal(0.);
#endif
}
void
EBFaceFAB::define(const EBISBox&  a_ebisBox,
            const Box& a_region,
            int a_iDir, int a_nComp)
{
  BaseEBFaceFAB<Real>::define(a_ebisBox, a_region, a_iDir, a_nComp);

#ifdef CH_USE_SETVAL
  setVal(BaseFabRealSetVal);
#else
  setVal(0.);
#endif
}

/**********************/
/**********************/
EBFaceFAB::~EBFaceFAB()
{
}

/**********************/
/**********************/
const FArrayBox&
EBFaceFAB::getFArrayBox() const
{
  CH_assert(isDefined());
  return (const FArrayBox&)m_regFAB;
}
/**********************/
/**********************/
FArrayBox&
EBFaceFAB::getFArrayBox()
{
  CH_assert(isDefined());
  return (FArrayBox&)m_regFAB;
}

/**********************/
/**********************/
EBFaceFAB&
EBFaceFAB::operator+=(const EBFaceFAB& a_src)
{
  CH_assert(a_src.m_nComp == m_nComp);
  plus(a_src, 0, 0, m_nComp);

  return *this;
}

EBFaceFAB&
EBFaceFAB::plus(const EBFaceFAB& a_src,
                int a_srccomp,
                int a_destcomp,
                int a_numcomp)
{
  CH_assert(isDefined());
  CH_assert(a_src.isDefined());
  //  CH_assert(m_ebisBox == a_src.m_ebisBox);
  CH_assert(a_srccomp + a_numcomp <= a_src.m_nComp);
  CH_assert(a_destcomp + a_numcomp <= m_nComp);

  Box locRegionFace = a_src.m_regionFace & m_regionFace;
  Box locRegion = a_src.m_region & m_region;
  if (!locRegion.isEmpty())
    {
      FORT_ADDTWOFAB(CHF_FRA(m_regFAB),
                     CHF_CONST_FRA(a_src.m_regFAB),
                     CHF_BOX(locRegionFace),
                     CHF_INT(a_srccomp),
                     CHF_INT(a_destcomp),
                     CHF_INT(a_numcomp));

      const Vector<FaceIndex>& faces = m_irrFAB.getFaces();
      for (int iface = 0; iface < faces.size(); iface++)
        {
          const FaceIndex& face = faces[iface];
          if (locRegionFace.contains(face.gridIndex(Side::Hi)))
            {
              for (int icomp = 0; icomp < a_numcomp; ++icomp)
                {
                  m_irrFAB(face, a_destcomp+icomp) +=
                    a_src.m_irrFAB(face, a_srccomp+icomp);
                }
            }
        }
    }
  return *this;
}

/**********************/
/**********************/
EBFaceFAB&
EBFaceFAB::operator-=(const EBFaceFAB& a_src)
{
  CH_assert(a_src.m_nComp == m_nComp);

  minus(a_src, 0, 0, m_nComp);

  return *this;
}

EBFaceFAB&
EBFaceFAB::minus(const EBFaceFAB& a_src,
                 int a_srccomp,
                 int a_destcomp,
                 int a_numcomp)
{
  CH_assert(isDefined());
  CH_assert(a_src.isDefined());
  // (DFM 7/28/05) This assertion fails for multifluid cases
  // where this and src have been defined using different
  // factories; what we really need is a better implementation
  // of the operator== for EBISBox
  //CH_assert(m_ebisBox == a_src.m_ebisBox);
  CH_assert(a_srccomp + a_numcomp <= a_src.m_nComp);
  CH_assert(a_destcomp + a_numcomp <= m_nComp);

  Box locRegionFace = a_src.m_regionFace & m_regionFace;
  Box locRegion = a_src.m_region & m_region;
  if (!locRegion.isEmpty())
    {
      FORT_SUBTRACTTWOFAB(CHF_FRA(m_regFAB),
                          CHF_CONST_FRA(a_src.m_regFAB),
                          CHF_BOX(locRegionFace),
                          CHF_INT(a_srccomp),
                          CHF_INT(a_destcomp),
                          CHF_INT(a_numcomp));

      const Vector<FaceIndex>& faces = m_irrFAB.getFaces();
      for (int iface = 0; iface < faces.size(); iface++)
        {
          const FaceIndex& face = faces[iface];
          if (locRegionFace.contains(face.gridIndex(Side::Hi)))
            {
              for (int icomp = 0; icomp < a_numcomp; ++icomp)
                {
                  m_irrFAB(face, a_destcomp+icomp) -=
                    a_src.m_irrFAB(face, a_srccomp+icomp);
                }
            }
        }
    }
  return *this;
}

/**********************/
/**********************/
EBFaceFAB&
EBFaceFAB::operator*=(const EBFaceFAB& a_src)
{
  CH_assert(a_src.m_nComp == m_nComp);

  mult(a_src, 0,0,m_nComp);

  return *this;
}



EBFaceFAB&
EBFaceFAB::mult(const EBFaceFAB& a_src,
                int a_srccomp,
                int a_destcomp,
                int a_numcomp)
{
  CH_assert(isDefined());
  CH_assert(a_src.isDefined());
  CH_assert(m_ebisBox == a_src.m_ebisBox);
  CH_assert(a_srccomp + a_numcomp <= a_src.m_nComp);
  CH_assert(a_destcomp + a_numcomp <= m_nComp);

  Box locRegionFace = a_src.m_regionFace & m_regionFace;
  Box locRegion = a_src.m_region & m_region;
  if (!locRegion.isEmpty())
    {
      FORT_MULTIPLYTWOFAB(CHF_FRA(m_regFAB),
                          CHF_CONST_FRA(a_src.m_regFAB),
                          CHF_BOX(locRegionFace),
                          CHF_INT(a_srccomp),
                          CHF_INT(a_destcomp),
                          CHF_INT(a_numcomp));

      const Vector<FaceIndex>& faces = m_irrFAB.getFaces();
      for (int iface = 0; iface < faces.size(); iface++)
        {
          const FaceIndex& face = faces[iface];
          if (locRegionFace.contains(face.gridIndex(Side::Hi)))
            {
              for (int icomp = 0; icomp < a_numcomp; ++icomp)
                {
                  m_irrFAB(face, a_destcomp+icomp) *=
                    a_src.m_irrFAB(face, a_srccomp+icomp);
                }
            }
        }
    }
  return *this;
}

/**********************/
/**********************/
EBFaceFAB&
EBFaceFAB::operator/=(const EBFaceFAB& a_src)
{
  CH_assert(a_src.m_nComp == m_nComp);

  divide(a_src, 0, 0, m_nComp);

  return *this;
}


EBFaceFAB&
EBFaceFAB::divide(const EBFaceFAB& a_src,
                  int a_srccomp,
                  int a_destcomp,
                  int a_numcomp)
{
  CH_assert(isDefined());
  CH_assert(a_src.isDefined());
  CH_assert(m_ebisBox == a_src.m_ebisBox);
  CH_assert(a_srccomp + a_numcomp <= a_src.m_nComp);
  CH_assert(a_destcomp + a_numcomp <= m_nComp);

  Box locRegionFace = a_src.m_regionFace & m_regionFace;
  Box locRegion = a_src.m_region & m_region;
  if (!locRegion.isEmpty())
    {
      FORT_DIVIDETWOFAB(CHF_FRA(m_regFAB),
                        CHF_CONST_FRA(a_src.m_regFAB),
                        CHF_BOX(locRegionFace),
                        CHF_INT(a_srccomp),
                        CHF_INT(a_destcomp),
                        CHF_INT(a_numcomp));

      const Vector<FaceIndex>& faces = m_irrFAB.getFaces();
      for (int iface = 0; iface < faces.size(); iface++)
        {
          const FaceIndex& face = faces[iface];
          if (locRegionFace.contains(face.gridIndex(Side::Hi)))
            {
              for (int icomp = 0; icomp < a_numcomp; ++icomp)
                {
                  m_irrFAB(face, a_destcomp+icomp) /=
                    a_src.m_irrFAB(face, a_srccomp+icomp);
                }
            }
        }
    }
  return *this;
}

/**********************/
/**********************/
EBFaceFAB&
EBFaceFAB::operator+=(const Real& a_src)
{
  CH_assert(isDefined());

  FORT_ADDFABR(CHF_FRA(m_regFAB),
               CHF_CONST_REAL(a_src),
               CHF_BOX(m_regionFace));

  const Vector<FaceIndex>& faces = m_irrFAB.getFaces();
  for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      for (int icomp = 0; icomp < m_nComp; ++icomp)
        {
          m_irrFAB(face, icomp) += a_src;
        }
    }
  return *this;
}

/**********************/
/**********************/
EBFaceFAB&
EBFaceFAB::operator-=(const Real& a_src)
{
  CH_assert(isDefined());

  FORT_SUBTRACTFABR(CHF_FRA(m_regFAB),
                    CHF_CONST_REAL(a_src),
                    CHF_BOX(m_regionFace));

  const Vector<FaceIndex>& faces = m_irrFAB.getFaces();
  for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      for (int icomp = 0; icomp < m_nComp; ++icomp)
        {
          m_irrFAB(face, icomp) -= a_src;
        }
    }
  return *this;
}
/**********************/
/**********************/
EBFaceFAB&
EBFaceFAB::operator*=(const Real& a_src)
{
  CH_assert(isDefined());

  FORT_MULTIPLYFABR(CHF_FRA(m_regFAB),
                    CHF_CONST_REAL(a_src),
                    CHF_BOX(m_regionFace));

  const Vector<FaceIndex>& faces = m_irrFAB.getFaces();
  for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      for (int icomp = 0; icomp < m_nComp; ++icomp)
        {
          m_irrFAB(face, icomp) *= a_src;
        }
    }
  return *this;
}
/**********************/
/**********************/
EBFaceFAB&
EBFaceFAB::operator/=(const Real& a_src)
{
  CH_assert(isDefined());

  FORT_DIVIDEFABR(CHF_FRA(m_regFAB),
                  CHF_CONST_REAL(a_src),
                  CHF_BOX(m_regionFace));

  const Vector<FaceIndex>& faces = m_irrFAB.getFaces();
  for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      for (int icomp = 0; icomp < m_nComp; ++icomp)
        {
          m_irrFAB(face, icomp) /= a_src;
        }
    }
  return *this;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
EBFaceFAB::max(int a_comp) const
{
  CH_assert(isDefined());
  Real val = -DBL_MAX;

  // Find the min on regular cells.
  FORT_MAXFAB(CHF_REAL(val), CHF_FRA(m_regFAB), CHF_BOX(m_regionFace),
              CHF_CONST_INT(a_comp));

  // Find the max on irregular faces.
  const Vector<FaceIndex>& faces = m_irrFAB.getFaces();
  for (int iface = 0; iface < faces.size(); iface++)
  {
    const FaceIndex& face = faces[iface];
    for (int icomp = 0; icomp < m_nComp; ++icomp)
    {
      val = Max(val, m_irrFAB(face, icomp));
    }
  }
  return val;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
EBFaceFAB::min(int a_comp) const
{
  CH_assert(isDefined());
  Real val = DBL_MAX;

  // Find the min on regular cells.
  FORT_MINFAB(CHF_REAL(val), CHF_FRA(m_regFAB), CHF_BOX(m_regionFace),
              CHF_CONST_INT(a_comp));

  // Find the max on irregular faces.
  const Vector<FaceIndex>& faces = m_irrFAB.getFaces();
  for (int iface = 0; iface < faces.size(); iface++)
  {
    const FaceIndex& face = faces[iface];
    for (int icomp = 0; icomp < m_nComp; ++icomp)
    {
      val = Min(val, m_irrFAB(face, icomp));
    }
  }
  return val;
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
