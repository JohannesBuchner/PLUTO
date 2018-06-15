#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include <cstring>
#include "SPACE.H"
using std::cin;
using std::cout;
using std::cerr;
using std::setw;
using std::setprecision;
using std::ios;
using std::pow;
using std::sqrt;

#include "Misc.H"
#include "FArrayBox.H"
#include "MayDay.H"
#include "NamespaceHeader.H"

FArrayBox::FArrayBox()
  :
  BaseFab<Real>()
{
}

FArrayBox::FArrayBox(const Box& a_box,
                     int        a_n,
                     Real*      a_alias)
  :
  BaseFab<Real>(a_box,a_n,a_alias)
{
  // Note: All work is done in BaseFab<Real> constructor
}

FArrayBox::~FArrayBox()
{
}

Real FArrayBox::norm(const Box& a_subbox,
                     int        a_p,
                     int        a_comp,
                     int        a_numcomp) const
{
  CH_assert(a_comp >= 0 && a_comp + a_numcomp <= nComp());
  CH_assert(a_numcomp > 0);
  CH_assert(a_p >= 0);
  CH_assert(m_domain.contains(a_subbox));

  Real* tmp = 0;
  int tmplen = 0;
  Real nrm = 0;

  if (a_p == 0)
    {
      // here begins a normed fab function piece
      ForAllThisCPencil(Real,a_subbox,a_comp,a_numcomp)
        {
          const Real* row = &thisR;
          if (tmp == 0)
            {
              tmp = new Real[thisLen];
              tmplen = thisLen;
              for (int i = 0; i < thisLen; i++)
                {
                  tmp[i] = Abs(Real(row[i]));
                }
            }
          else
            {
              for (int i = 0; i < thisLen; i++)
                {
                  tmp[i] = Max(tmp[i],Real(Abs(row[i])));
                }
            }
        } EndForPencil

      nrm = tmp[0];
      for (int i = 1; i < tmplen; i++)
        {
          nrm = Max(nrm, tmp[i]);
        }
      // here it ends
    }
  else if (a_p == 1)
    {
      // here begins a normed fab function piece
      ForAllThisCPencil(Real,a_subbox,a_comp,a_numcomp)
        {
          const Real* row = &thisR;
          if (tmp == 0)
            {
              tmp = new Real[thisLen];
              tmplen = thisLen;
              for (int i = 0; i < thisLen; i++)
                {
                  tmp[i] = Abs(Real(row[i]));
                }
            }
          else
            {
              for (int i = 0; i < thisLen; i++)
                {
                  tmp[i] += Abs(Real(row[i]));
                }
            }
        } EndForPencil

      nrm = tmp[0];
      for (int i = 1; i < tmplen; i++)
        {
          nrm += tmp[i];
        }
      // here it ends
    }
  else if (a_p == 2)
    {
      nrm = sqrt(sumPow(a_subbox, 2, a_comp, a_numcomp));
    }
  else
    {
      // so standard norms weren't good enough for you?
      Real invpwr = 1.0/a_p;
      nrm = pow(sumPow(a_subbox, a_p, a_comp, a_numcomp),invpwr);
    }

  delete [] tmp;

  return nrm;
}

Real FArrayBox::norm(int a_p,
                     int a_comp,
                     int a_numcomp) const
{
  return norm(m_domain,a_p,a_comp,a_numcomp);
}

// utility function used in norms and such only works for a_p >= 2
Real FArrayBox::sumPow(const Box& a_subbox,
                       int        a_p,
                       int        a_comp,
                       int        a_numcomp) const
{
  CH_assert(a_p >= 2);
  CH_assert(a_numcomp > 0);
  CH_assert(m_domain.contains(a_subbox));

  Real sum = 0;
  Real* tmp = NULL;
  int tmplen = 0;

  if (a_p == 2)
  {
    ForAllThisCPencil(Real,a_subbox,a_comp,a_numcomp)
    {
      const Real* row = &thisR;
      if (tmp == 0)
      {
        tmp = new Real[thisLen];
        tmplen = thisLen;
        for (int i = 0; i < thisLen; i++)
        {
          tmp[i] = row[i]*row[i];
        }
      }
      else
      {
        for (int i = 0; i < thisLen; i++)
        {
          tmp[i] += row[i]*row[i];
        }
      }
    } EndForPencil

    sum = tmp[0];
    for (int i = 1; i < tmplen; i++)
    {
      sum += tmp[i];
    }
  }
  else
  {
    // so standard norms weren't good enough for you?
    Real pwr = Real(a_p);

    ForAllThisCPencil(Real,a_subbox,a_comp,a_numcomp)
    {
      const Real* row = &thisR;
      if (tmp == 0)
      {
        tmp = new Real[thisLen];
        tmplen = thisLen;
        for (int i = 0; i < thisLen; i++)
        {
          tmp[i] = pow(row[i],pwr);
        }
      }
      else
      {
        for (int i = 0; i < thisLen; i++)
        {
          tmp[i] += pow(row[i],pwr);
        }
      }
    } EndForPencil

    sum = tmp[0];
    for (int i = 1; i < tmplen; i++)
    {
      sum += tmp[i];
    }
  }

  delete [] tmp;

  return sum;
}

// Take the dot product of "this" and "a_fab2" over their common box and
// all components.
Real FArrayBox::dotProduct(const FArrayBox& a_fab2) const
{
  const FArrayBox& fab1 = *this;
  Box commonBox = fab1.box() & a_fab2.box();
  return dotProduct(a_fab2, commonBox);
}

Real FArrayBox::dotProduct(const FArrayBox& a_fab2, const Box& a_box) const
{
  Real dot = 0.0;
  const FArrayBox& fab1 = *this;

  int startcomp = 0;
  int endcomp = fab1.nComp()-1;
  int numcomp = endcomp+1;

  CH_assert(fab1.nComp() == a_fab2.nComp());

  ForAllThisCBNNXC(Real, a_box, startcomp, numcomp, a_fab2, startcomp)
    {
      dot += thisR * a_fab2R;
    } EndForTX

  return dot;
}

//here is where the orderedfab functions start

Real FArrayBox::min(int a_comp) const
{
  Real *_min_row = 0;
  int _X_len = 0;

  ForAllThisCPencil(Real,m_domain,a_comp,1)
  {
    const Real* _row = &thisR;
    if (_min_row == 0)
    {
      _min_row = new Real[thisLen];
      _X_len = thisLen;
      for (int i = 0; i < thisLen; i++)
      {
        _min_row[i] = _row[i];
      }
    }
    else
    {
      for (int i = 0; i < thisLen; i++)
      {
        _min_row[i] = Min(_row[i],_min_row[i]);
      }
    }
  } EndForPencil;

  Real _min = _min_row[0];
  for (int i = 1; i < _X_len; i++)
  {
    _min = Min(_min,_min_row[i]);
  }

  delete [] _min_row;

  return _min;
}

Real FArrayBox::min(const Box& a_subbox,
                    int        a_comp) const
{
  CH_assert(m_domain.contains(a_subbox));

  Real *_min_row = 0;
  int _X_len = 0;

  ForAllThisCPencil(Real,a_subbox,a_comp,1)
  {
    const Real* _row = &thisR;
    if (_min_row == 0)
    {
      _min_row = new Real[thisLen];
      _X_len = thisLen;
      for (int i = 0; i < thisLen; i++)
      {
        _min_row[i] = _row[i];
      }
    }
    else
    {
      for (int i = 0; i < thisLen; i++)
      {
        _min_row[i] = Min(_row[i],_min_row[i]);
      }
    }
  } EndForPencil;

  Real _min = _min_row[0];
  for (int i = 1; i < _X_len; i++)
  {
    _min = Min(_min,_min_row[i]);
  }

  delete [] _min_row;

  return _min;
}

Real FArrayBox::max(int a_comp) const
{
  Real *_max_row = 0;
  int _X_len = 0;

  ForAllThisCPencil(Real,m_domain,a_comp,1)
  {
    const Real* _row = &thisR;
    if (_max_row== 0)
    {
      _max_row = new Real[thisLen];
      _X_len = thisLen;
      for (int i = 0; i < thisLen; i++)
      {
        _max_row[i] = _row[i];
      }
    }
    else
    {
      for (int i = 0; i < thisLen; i++)
      {
        _max_row[i] = Max(_row[i],_max_row[i]);
      }
    }
  } EndForPencil;

  Real _max = _max_row[0];
  for (int i = 1; i < _X_len; i++)
  {
    _max = Max(_max,_max_row[i]);
  }

  delete [] _max_row;

  return _max;
}

Real FArrayBox::max(const Box& a_subbox,
                    int        a_comp) const
{
  CH_assert(m_domain.contains(a_subbox));

  Real *_max_row = 0;
  int _X_len = 0;

  ForAllThisCPencil(Real,a_subbox,a_comp,1)
  {
    const Real* _row = &thisR;
    if (_max_row == 0)
    {
      _max_row = new Real[thisLen];
      _X_len = thisLen;
      for (int i = 0; i < thisLen; i++)
      {
        _max_row[i] = _row[i];
      }
    }
    else
    {
      for (int i = 0; i < thisLen; i++)
      {
        _max_row[i] = Max(_row[i],_max_row[i]);
      }
    }
  } EndForPencil;

  Real _max = _max_row[0];
  for (int i = 1; i < _X_len; i++)
  {
    _max = Max(_max,_max_row[i]);
  }

  delete [] _max_row;

  return _max;
}

IntVect FArrayBox::minIndex(int a_comp) const
{
  IntVect _min_loc(m_domain.smallEnd());
  Real _min_val = (*this).operator()(_min_loc,a_comp);

  ForAllThisCBNN(Real,m_domain,a_comp,1)
  {
    if (thisR < _min_val)
    {
      _min_val = thisR;
      D_EXPR6(_min_loc[0] = iR,
              _min_loc[1] = jR,
              _min_loc[2] = kR,
              _min_loc[3] = _iv[3],
              _min_loc[4] = _iv[4],
              _min_loc[5] = _iv[5]);
    }
  } EndFor

  return _min_loc;
}

IntVect FArrayBox::minIndex(const Box& a_subbox,
                            int        a_comp) const
{
  CH_assert(m_domain.contains(a_subbox));

  IntVect _min_loc(a_subbox.smallEnd());
  Real _min_val = (*this).operator()(_min_loc,a_comp);

  ForAllThisCBNN(Real,a_subbox,a_comp,1)
  {
    if (thisR < _min_val)
    {
      _min_val = thisR;
      D_EXPR6(_min_loc[0] = iR,
              _min_loc[1] = jR,
              _min_loc[2] = kR,
              _min_loc[3] = _iv[3],
              _min_loc[4] = _iv[4],
              _min_loc[5] = _iv[5]);
    }
  } EndFor

  return _min_loc;
}

IntVect FArrayBox::maxIndex(int a_comp) const
{
  IntVect _max_loc(m_domain.smallEnd());
  Real _max_val = (*this).operator()(_max_loc,a_comp);

  ForAllThisCBNN(Real,m_domain,a_comp,1)
  {
    if (thisR > _max_val)
    {
      _max_val = thisR;
      D_EXPR6(_max_loc[0] = iR,
              _max_loc[1] = jR,
              _max_loc[2] = kR,
              _max_loc[3] = _iv[3],
              _max_loc[4] = _iv[4],
              _max_loc[5] = _iv[5]);
    }
  } EndFor

  return _max_loc;
}

IntVect FArrayBox::maxIndex(const Box& a_subbox,
                            int        a_comp) const
{
  CH_assert(m_domain.contains(a_subbox));

  IntVect _max_loc(a_subbox.smallEnd());
  Real _max_val = (*this).operator()(_max_loc,a_comp);

  ForAllThisCBNN(Real,a_subbox,a_comp,1)
  {
    if (thisR > _max_val)
    {
      _max_val = thisR;
      D_EXPR6(_max_loc[0] = iR,
              _max_loc[1] = jR,
              _max_loc[2] = kR,
              _max_loc[3] = _iv[3],
              _max_loc[4] = _iv[4],
              _max_loc[5] = _iv[5]);
    }
  } EndFor

  return _max_loc;
}

int FArrayBox::maskLT(BaseFab<int>& a_mask,
                      Real          a_val,
                      int           a_comp) const
{
  a_mask.resize(m_domain,1);
  a_mask.setVal(0);

  int *mptr = a_mask.dataPtr();
  int cnt = 0;

  // because this is a BaseFab<Real>, and mask is a
  // BaseFab<int>, we can't use the BaseFabMacros
  // involving two BaseFabs, since T is different.
  // so, access the mask value by incrementing the
  // mask pointer with each pass through the loop
  ForAllThisCBNN(Real,m_domain,a_comp,1)
  {
    if (thisR < a_val)
    {
      *mptr = 1;
      cnt++;
    }
    mptr++;
  } EndFor

  return cnt;
}

int FArrayBox::maskLE(BaseFab<int>& a_mask,
                      Real          a_val,
                      int           a_comp) const
{
  a_mask.resize(m_domain,1);
  a_mask.setVal(0);

  int *mptr = a_mask.dataPtr();
  int cnt = 0;

  // because this is a BaseFab<Real>, and mask is a
  // BaseFab<int>, we can't use the BaseFabMacros
  // involving two BaseFabs, since T is different.
  // so, access the mask value by incrementing the
  // mask pointer with each pass through the loop
  ForAllThisCBNN(Real,m_domain,a_comp,1)
  {
    if (thisR <= a_val)
    {
      *mptr = 1;
      cnt++;
    }
    mptr++;
  } EndFor

  return cnt;
}

int FArrayBox::maskEQ(BaseFab<int>& a_mask,
                      Real          a_val,
                      int           a_comp) const
{
  a_mask.resize(m_domain,1);
  a_mask.setVal(0);

  int *mptr = a_mask.dataPtr();
  int cnt = 0;

  // because this is a BaseFab<Real>, and mask is a
  // BaseFab<int>, we can't use the BaseFabMacros
  // involving two BaseFabs, since T is different.
  // so, access the mask value by incrementing the
  // mask pointer with each pass through the loop
  ForAllThisCBNN(Real,m_domain,a_comp,1)
  {
    if (thisR == a_val)
    {
      *mptr = 1;
      cnt++;
    }
    mptr++;
  } EndFor

  return cnt;
}

int FArrayBox::maskGT(BaseFab<int>& a_mask,
                      Real          a_val,
                      int           a_comp) const
{
  a_mask.resize(m_domain,1);
  a_mask.setVal(0);

  int *mptr = a_mask.dataPtr();
  int cnt = 0;

  // because this is a BaseFab<Real>, and mask is a
  // BaseFab<int>, we can't use the BaseFabMacros
  // involving two BaseFabs, since T is different.
  // so, access the mask value by incrementing the
  // mask pointer with each pass through the loop
  ForAllThisCBNN(Real,m_domain,a_comp,1)
  {

    if (thisR > a_val)
    {
      *mptr = 1;
      cnt++;
    }
    mptr++;
  } EndFor

  return cnt;
}

int FArrayBox::maskGE(BaseFab<int>& a_mask,
                      Real          a_val,
                      int           a_comp) const
{
  a_mask.resize(m_domain,1);
  a_mask.setVal(0);

  int *mptr = a_mask.dataPtr();
  int cnt = 0;

  // because this is a BaseFab<Real>, and mask is a
  // BaseFab<int>, we can't use the BaseFabMacros
  // involving two BaseFabs, since T is different.
  // so, access the mask value by incrementing the
  // mask pointer with each pass through the loop
  ForAllThisCBNN(Real,m_domain,a_comp,1)
  {
    if (thisR >= a_val)
    {
      *mptr = 1;
      cnt++;
    }
    mptr++;
  } EndFor

  return cnt;
}

//here is where the orderedfab functions end

//here is where the normedfab functions begin

void FArrayBox::abs()
{
  ForAllThis(Real)
  {
    thisR = Abs(thisR);
  } EndFor
}

void FArrayBox::abs(int a_comp,
                    int a_numcomp)
{
  ForAllThisNN(Real,a_comp,a_numcomp)
  {
    thisR = Abs(thisR);
  } EndFor
}

void FArrayBox::abs(const Box& a_subbox,
                    int        a_comp,
                    int        a_numcomp)
{
  CH_assert(m_domain.contains(a_subbox));

  ForAllThisBNN(Real,a_subbox,a_comp,a_numcomp)
  {
    thisR = Abs(thisR);
  } EndFor
}

//here is where the normedfab functions end

//here is where the arithfab functions begin

Real FArrayBox::sum(int a_comp,
                    int a_numcomp) const
{
  Real *_sum_row = 0;
  int _sum_len = 0;

  ForAllThisCPencil(Real,m_domain,a_comp,a_numcomp)
  {
    const Real* _row = &thisR;
    if (_sum_row == 0)
    {
      _sum_row = new Real[thisLen];
      _sum_len = thisLen;
      for (int i = 0; i < thisLen; i++)
      {
        _sum_row[i] = _row[i];
      }
    }
    else
    {
      for (int i = 0; i < thisLen; i++)
      {
        _sum_row[i] += _row[i];
      }
    }
  } EndForPencil;

  Real _sum = _sum_row[0];
  for (int i = 1; i < _sum_len; i++)
  {
    _sum += _sum_row[i];
  }

  delete [] _sum_row;

  return _sum;
}

Real FArrayBox::sum(const Box& a_subbox,
                    int        a_comp,
                    int        a_numcomp) const
{
  CH_assert(m_domain.contains(a_subbox));

  Real *_sum_row = 0;
  int _sum_len = 0;

  ForAllThisCPencil(Real,a_subbox,a_comp,a_numcomp)
  {
    const Real* _row = &thisR;
    if (_sum_row == 0)
    {
      _sum_row = new Real[thisLen];
      _sum_len = thisLen;
      for (int i = 0; i < thisLen; i++)
      {
        _sum_row[i] = _row[i];
      }
    }
    else
      {
        for (int i = 0; i < thisLen; i++)
        {
          _sum_row[i] += _row[i];
        }
      }
  } EndForPencil;

  Real _sum = _sum_row[0];
  for (int i = 1; i < _sum_len; i++)
  {
    _sum += _sum_row[i];
  }

  delete [] _sum_row;

  return _sum;
}

FArrayBox& FArrayBox::invert(Real a_r)
{
  ForAllThis(Real)
  {
    thisR = a_r / thisR;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::invert(Real a_r,
                             int  a_comp,
                             int  a_numcomp)
{
  ForAllThisNN(Real,a_comp,a_numcomp)
  {
    thisR = a_r / thisR;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::invert(Real       a_r,
                             const Box& a_subbox,
                             int        a_comp,
                             int        a_numcomp)
{
  ForAllThisBNN(Real,a_subbox,a_comp,a_numcomp)
  {
    thisR = a_r / thisR;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::negate(const Box& a_subbox,
                             int        a_comp,
                             int        a_numcomp)
{
  ForAllThisBNN(Real,a_subbox,a_comp,a_numcomp)
  {
    thisR = - thisR;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::negate(int a_comp,
                             int a_numcomp)
{
  ForAllThisNN(Real,a_comp,a_numcomp)
  {
    thisR = - thisR;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::negate()
{
  ForAllThis(Real)
  {
    thisR = - thisR;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::plus(Real       a_r,
                           const Box& a_subbox,
                           int        a_comp,
                           int        a_numcomp)
{
  ForAllThisBNN(Real,a_subbox,a_comp,a_numcomp)
  {
    thisR += a_r;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::plus(Real a_r,
                           int  a_comp,
                           int  a_numcomp)
{
  ForAllThisNN(Real,a_comp,a_numcomp)
  {
    thisR += a_r;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::operator += (Real a_r)
{
  ForAllThis(Real)
  {
    thisR += a_r;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::operator += (const FArrayBox& a_x)
{
  ForAllThisXC(Real,a_x)
  {
    thisR += a_xR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::plus(Real a_r)
{
  return operator += (a_r);
}

FArrayBox& FArrayBox::plus(const FArrayBox& a_x)
{
  return operator += (a_x);
}

// added bvs Tue May 18, PDT 1999
FArrayBox& FArrayBox::plus(const FArrayBox& a_src,
                           const Real&      a_scale)
{
  ForAllThisBNNXC(Real, a_src.box() & m_domain, 0, m_nvar, a_src, 0)
  {
    thisR += a_srcR * a_scale;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::plus(const FArrayBox& a_src,
                           const Real&      a_scale,
                           int              a_srccomp,
                           int              a_destcomp,
                           int              a_numcomp)
{
  ForAllThisBNNXC(Real, a_src.box() & m_domain, a_destcomp, a_numcomp, a_src, a_srccomp)
  {
    thisR += a_srcR * a_scale;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::plus(const FArrayBox& a_src,
                           int              a_srccomp,
                           int              a_destcomp,
                           int              a_numcomp)
{
  ForAllThisBNNXC(Real,m_domain,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR += a_srcR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::plus(const FArrayBox& a_src,
                           const Box&       a_subbox,
                           int              a_srccomp,
                           int              a_destcomp,
                           int              a_numcomp)
{
  CH_assert(m_domain.contains(a_subbox));

  ForAllThisBNNXC(Real,a_subbox,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR += a_srcR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::plus(const FArrayBox& a_src,
                           const Box&       a_srcbox,
                           const Box&       a_destbox,
                           int              a_srccomp,
                           int              a_destcomp,
                           int              a_numcomp)
{
  ForAllThisBNNXCBN(Real,a_destbox,a_destcomp,a_numcomp,a_src,a_srcbox,a_srccomp)
  {
    thisR += a_srcR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::plus(const FArrayBox& a_src,
                           const Box&       a_srcbox,
                           const Box&       a_destbox,
                           const Real&      a_scale,
                           int              a_srccomp,
                           int              a_destcomp,
                           int              a_numcomp)
{
  ForAllThisBNNXCBN(Real,a_destbox,a_destcomp,a_numcomp,a_src,a_srcbox,a_srccomp)
  {
    thisR += a_srcR * a_scale;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::operator -= (Real a_r)
{
  return operator += (-a_r);
}

FArrayBox& FArrayBox::operator -= (const FArrayBox& a_x)
{
  ForAllThisXC(Real,a_x)
  {
    thisR -= a_xR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::minus(const FArrayBox& a_x)
{
  return operator -= (a_x);
}

FArrayBox& FArrayBox::minus(const FArrayBox& a_src,
                           int               a_srccomp,
                           int               a_destcomp,
                           int               a_numcomp)
{
  ForAllThisBNNXC(Real,m_domain,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR -= a_srcR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::minus(const FArrayBox& a_src,
                            const Box&       a_subbox,
                            int              a_srccomp,
                            int              a_destcomp,
                            int              a_numcomp)
{
  CH_assert(m_domain.contains(a_subbox));

  ForAllThisBNNXC(Real,a_subbox,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR -= a_srcR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::minus(const FArrayBox& a_src,
                            const Box&       a_srcbox,
                            const Box&       a_destbox,
                            int              a_srccomp,
                            int              a_destcomp,
                            int              a_numcomp)
{
  ForAllThisBNNXCBN(Real,a_destbox,a_destcomp,a_numcomp,a_src,a_srcbox,a_srccomp)
  {
    thisR -= a_srcR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::operator *= (Real a_r)
{
  ForAllThis(Real)
  {
    thisR *= a_r;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::mult(Real a_r)
{
  return operator *= (a_r);
}

FArrayBox& FArrayBox::mult(Real a_r,
                           int  a_comp,
                           int  a_numcomp)
{
  ForAllThisNN(Real,a_comp,a_numcomp)
  {
    thisR *= a_r;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::mult(Real       a_r,
                           const Box& a_subbox,
                           int        a_comp,
                           int        a_numcomp)
{
  ForAllThisBNN(Real,a_subbox,a_comp,a_numcomp)
  {
    thisR *= a_r;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::operator *= (const FArrayBox &a_x)
{
  ForAllThisXC(Real,a_x)
  {
    thisR *= a_xR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::mult(const FArrayBox& a_x)
{
  return operator *= (a_x);
}

FArrayBox& FArrayBox::mult(const FArrayBox& a_src,
                           int              a_srccomp,
                           int              a_destcomp,
                           int              a_numcomp)
{
  ForAllThisBNNXC(Real,m_domain,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR *= a_srcR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::mult(const FArrayBox& a_src,
                           const Box&       a_subbox,
                           int              a_srccomp,
                           int              a_destcomp,
                           int              a_numcomp)
{
  CH_assert(m_domain.contains(a_subbox));

  ForAllThisBNNXC(Real,a_subbox,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR *= a_srcR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::mult(const FArrayBox& a_src,
                           const Box&       a_srcbox,
                           const Box&       a_destbox,
                           int              a_srccomp,
                           int              a_destcomp,
                           int              a_numcomp)
{
  ForAllThisBNNXCBN(Real,a_destbox,a_destcomp,a_numcomp,a_src,a_srcbox,a_srccomp)
  {
    thisR *= a_srcR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::operator /= (Real a_r)
{
  ForAllThis(Real)
  {
    thisR /= a_r;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::divide(Real a_r)
{
  return operator /= (a_r);
}

FArrayBox& FArrayBox::divide(Real a_r,
                             int  a_comp,
                             int  a_numcomp)
{
  ForAllThisNN(Real,a_comp,a_numcomp)
  {
    thisR /= a_r;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::divide(Real       a_r,
                             const Box& a_subbox,
                             int        a_comp,
                             int        a_numcomp)
{
  ForAllThisBNN(Real,a_subbox,a_comp,a_numcomp)
  {
    thisR /= a_r;
  } EndFor

  return *this;
}

FArrayBox& FArrayBox::operator /= (const FArrayBox &a_x)
{
  ForAllThisXC(Real,a_x)
  {
    thisR /= a_xR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::divide(const FArrayBox& a_x)
{
  return operator /= (a_x);
}

FArrayBox& FArrayBox::divide(const FArrayBox& a_src,
                             int              a_srccomp,
                             int              a_destcomp,
                             int              a_numcomp)
{
  ForAllThisBNNXC(Real,m_domain,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR /= a_srcR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::divide(const FArrayBox& a_src,
                             const Box&       a_subbox,
                             int              a_srccomp,
                             int              a_destcomp,
                             int              a_numcomp)
{
  CH_assert(m_domain.contains(a_subbox));
  ForAllThisBNNXC(Real,a_subbox,a_destcomp,a_numcomp,a_src,a_srccomp)
  {
    thisR /= a_srcR;
  } EndForTX

  return *this;
}

FArrayBox& FArrayBox::divide(const FArrayBox& a_src,
                             const Box&       a_srcbox,
                             const Box&       a_destbox,
                             int              a_srccomp,
                             int              a_destcomp,
                             int              a_numcomp)
{
  ForAllThisBNNXCBN(Real,a_destbox,a_destcomp,a_numcomp,a_src,a_srcbox,a_srccomp)
  {
    thisR /= a_srcR;
  } EndForTX

  return *this;
}

//-----------------------------------------------------------------------
FArrayBox&
FArrayBox::
axby(const FArrayBox& a_X, const FArrayBox& a_Y,
     Real a_A, Real a_B)
{
  ForAllThisBNNXCBNYCBN(Real, m_domain, 0, nComp(), a_X, a_X.m_domain, \
                        0, a_Y, a_Y.m_domain, 0)
  {
    thisR = a_A * a_XR + a_B * a_YR;
  } EndForTX

  return *this;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
FArrayBox&
FArrayBox::
axby(const FArrayBox& a_X, const FArrayBox& a_Y,
     Real a_A, Real a_B,
     int a_destComp, int a_xComp, int a_yComp)
{
  ForAllThisBNNXCBNYCBN(Real, m_domain, a_destComp, 1, a_X, a_X.m_domain, \
                        a_xComp, a_Y, a_Y.m_domain, a_yComp)
  {
    thisR = a_A * a_XR + a_B * a_YR;
  } EndForTX

  return *this;
}
//-----------------------------------------------------------------------

void FArrayBox::performCopy(const BaseFab<Real>& a_src,
                            const Box&           a_srcbox,
                            int                  a_srccomp,
                            const Box&           a_destbox,
                            int                  a_destcomp,
                            int                  a_numcomp)
{
  BaseFab<Real>::performCopy(a_src, a_srcbox, a_srccomp, a_destbox, a_destcomp, a_numcomp);

  // FArrayBox& dest = *this;
  // FORT_PERFORMCOPY(CHF_FRA(dest),
  //                  CHF_CONST_FRA(a_src),
  //                  CHF_BOX(a_destbox),
  //                  CHF_BOX(a_srcbox),
  //                  CHF_CONST_INT(a_srccomp),
  //                  CHF_CONST_INT(a_destcomp),
  //                  CHF_CONST_INT(a_numcomp));
}
#include "NamespaceFooter.H"
