#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _REALTENSOR_H_
#define _REALTENSOR_H_

#include <algorithm>

#include "SPACE.H"
#include "REAL.H"
#include "MayDay.H"

#include "NamespaceHeader.H"

//! \class RealTensor
//! This class stores tensors that represent outer products of RealVects.
class RealTensor
{
  public:

  RealTensor()
  {
    std::fill(m_T, m_T + SpaceDim*SpaceDim, 0.0);
  }

#if CH_SPACEDIM == 1
  explicit RealTensor(Real a_t11)
  {
    m_T[0] = a_t11;
  }
#elif CH_SPACEDIM == 2
  RealTensor(Real a_t11, Real a_t12,
             Real a_t21, Real a_t22)
  {
    m_T[0] = a_t11; m_T[1] = a_t12;
    m_T[2] = a_t21; m_T[3] = a_t22;
  }
#elif CH_SPACEDIM == 3
  RealTensor(Real a_t11, Real a_t12, Real a_t13,
             Real a_t21, Real a_t22, Real a_t23,
             Real a_t31, Real a_t32, Real a_t33)
  {
    m_T[0] = a_t11; m_T[1] = a_t12; m_T[2] = a_t13;
    m_T[3] = a_t21; m_T[4] = a_t22; m_T[5] = a_t23;
    m_T[6] = a_t31; m_T[7] = a_t32; m_T[8] = a_t33;
  }
#elif CH_SPACEDIM == 4
  RealTensor(Real a_t11, Real a_t12, Real a_t13, Real a_t14,
             Real a_t21, Real a_t22, Real a_t23, Real a_t24,
             Real a_t31, Real a_t32, Real a_t33, Real a_t34,
             Real a_t41, Real a_t42, Real a_t43, Real a_t44)
  {
    m_T[ 0] = a_t11; m_T[ 1] = a_t12; m_T[ 2] = a_t13; m_T[ 3] = a_t14;
    m_T[ 4] = a_t21; m_T[ 5] = a_t22; m_T[ 6] = a_t23; m_T[ 7] = a_t24;
    m_T[ 8] = a_t31; m_T[ 9] = a_t32; m_T[10] = a_t33; m_T[11] = a_t34;
    m_T[12] = a_t31; m_T[13] = a_t32; m_T[14] = a_t33; m_T[15] = a_t44;
  }
#elif CH_SPACEDIM == 5
  RealTensor(Real a_t11, Real a_t12, Real a_t13, Real a_t14, Real a_t15,
             Real a_t21, Real a_t22, Real a_t23, Real a_t24, Real a_t25,
             Real a_t31, Real a_t32, Real a_t33, Real a_t34, Real a_t35,
             Real a_t41, Real a_t42, Real a_t43, Real a_t44, Real a_t45,
             Real a_t51, Real a_t52, Real a_t53, Real a_t54, Real a_t55)
  {
    m_T[ 0] = a_t11; m_T[ 1] = a_t12; m_T[ 2] = a_t13; m_T[ 3] = a_t14; m_T[ 4] = a_t15;
    m_T[ 5] = a_t21; m_T[ 6] = a_t22; m_T[ 7] = a_t23; m_T[ 8] = a_t24; m_T[ 9] = a_t25;
    m_T[10] = a_t31; m_T[11] = a_t32; m_T[12] = a_t33; m_T[13] = a_t34; m_T[14] = a_t35;
    m_T[15] = a_t41; m_T[16] = a_t42; m_T[17] = a_t43; m_T[18] = a_t44; m_T[19] = a_t45;
    m_T[20] = a_t51; m_T[21] = a_t52; m_T[22] = a_t53; m_T[23] = a_t54; m_T[24] = a_t55;
  }
#elif CH_SPACEDIM == 6
  RealTensor(Real a_t11, Real a_t12, Real a_t13, Real a_t14, Real a_t15, Real a_t16,
             Real a_t21, Real a_t22, Real a_t23, Real a_t24, Real a_t25, Real a_t26,
             Real a_t31, Real a_t32, Real a_t33, Real a_t34, Real a_t35, Real a_t36,
             Real a_t41, Real a_t42, Real a_t43, Real a_t44, Real a_t45, Real a_t46,
             Real a_t51, Real a_t52, Real a_t53, Real a_t54, Real a_t55, Real a_t56,
             Real a_t61, Real a_t62, Real a_t63, Real a_t64, Real a_t65, Real a_t66)
  {
    m_T[ 0] = a_t11; m_T[ 1] = a_t12; m_T[ 2] = a_t13; m_T[ 3] = a_t14; m_T[ 4] = a_t15; m_T[ 5] = a_t16;
    m_T[ 6] = a_t21; m_T[ 7] = a_t22; m_T[ 8] = a_t23; m_T[ 9] = a_t24; m_T[10] = a_t25; m_T[11] = a_t26;
    m_T[12] = a_t31; m_T[13] = a_t32; m_T[14] = a_t33; m_T[15] = a_t34; m_T[16] = a_t35; m_T[17] = a_t36;
    m_T[18] = a_t41; m_T[19] = a_t42; m_T[20] = a_t43; m_T[21] = a_t44; m_T[22] = a_t45; m_T[23] = a_t46;
    m_T[24] = a_t51; m_T[25] = a_t52; m_T[26] = a_t53; m_T[27] = a_t54; m_T[28] = a_t55; m_T[29] = a_t56;
    m_T[30] = a_t61; m_T[31] = a_t62; m_T[32] = a_t63; m_T[33] = a_t64; m_T[34] = a_t65; m_T[35] = a_t66;
  }
#else
#error "RealTensor is only implemented for SpaceDim <= 6."
#endif

  // Copy constructor.
  RealTensor(const RealTensor& a_rhs)
  {
    std::copy(a_rhs.m_T, a_rhs.m_T + SpaceDim*SpaceDim, m_T);
  }

  //! Destructor.
  ~RealTensor()
  {
  }

  //! Assignment operator.
  RealTensor& operator=(const RealTensor& a_rhs)
  {
    if (this != &a_rhs)
    {
      std::copy(a_rhs.m_T, a_rhs.m_T + SpaceDim*SpaceDim, m_T);
    }

    return *this;
  }

  //! Returns the (i,j) component of this tensor.
  Real operator()(int a_i, int a_j) const
  {
    return m_T[SpaceDim*a_j + a_i];
  }

  Real& operator()(int a_i, int a_j)
  {
    return m_T[SpaceDim*a_j + a_i];
  }

  //! Computes the determinant of this tensor using cofactors.
  Real det() const
  {
    const Real* tensor = m_T;

    return recurDet(tensor,SpaceDim);
  }

  //! Computes the transpose of this tensor.
  RealTensor transpose() const
  {
    RealTensor t = *this;
#if CH_SPACEDIM >= 2
    std::swap(t(0, 1), t(1, 0));
#endif
#if CH_SPACEDIM >= 3
    std::swap(t(0, 2), t(2, 0));
    std::swap(t(1, 2), t(2, 1));
#endif
#if CH_SPACEDIM >= 4
    std::swap(t(0, 3), t(3, 0));
    std::swap(t(1, 3), t(3, 1));
    std::swap(t(2, 3), t(3, 2));
#endif
#if CH_SPACEDIM >= 5
    std::swap(t(0, 4), t(4, 0));
    std::swap(t(1, 4), t(4, 1));
    std::swap(t(2, 4), t(4, 2));
    std::swap(t(3, 4), t(4, 3));
#endif
#if CH_SPACEDIM == 6
    std::swap(t(0, 5), t(5, 0));
    std::swap(t(1, 5), t(5, 1));
    std::swap(t(2, 5), t(5, 2));
    std::swap(t(3, 5), t(5, 3));
    std::swap(t(4, 5), t(5, 4));
#endif
    return t;
  }

  private:

  // Compute the determinant of this tensor using cofactors.
  Real recurDet(const Real* a_tensor,
                const int&  a_size) const
  {
    Real retval = 0.0;

    if (a_size == 1)
    {
      retval += *a_tensor;
    }
    else
    {
      int subSize = a_size - 1;
      Real* subTensor = new Real[subSize*subSize];

      Real sign = 1.0;

      int row1 = 0;
      for (int col1 = 0; col1 < a_size; col1++)
      {
        int subRow = 0;
        for (int row2 = 0; row2 < a_size; row2++)
        {
          if (row2 != row1)
          {
            int subCol = 0;
            for (int col2 = 0; col2 < a_size; col2++)
            {
              if (col2 != col1)
              {
                subTensor[subSize*subRow + subCol] = a_tensor[a_size*row2 + col2];

                subCol++;
              }
            }

            subRow++;
          }
        }

        retval += sign * a_tensor[a_size*row1 + col1] * recurDet(subTensor,subSize);

        sign *= -1.0;
      }

      delete [] subTensor;
    }

    return retval;
  }

  Real m_T[CH_SPACEDIM*CH_SPACEDIM];
};

#include "NamespaceFooter.H"

#endif
