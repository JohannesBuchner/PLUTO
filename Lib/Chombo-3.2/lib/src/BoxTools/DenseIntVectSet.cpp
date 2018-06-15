#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DenseIntVectSet.H"
#include "BoxIterator.H"
#include "MayDay.H"
#include "ProblemDomain.H"
#include "SPACE.H"
#include "SPMD.H"
#include "CH_Timer.H"

#include "NamespaceHeader.H"

DenseIntVectSet::DenseIntVectSet(const Box& a_domain, bool init)
  :
  m_domain(a_domain),
  m_bits(a_domain.numPts(), init),
  m_minBox(a_domain)
{
}

bool DenseIntVectSet::operator[](const IntVect& index) const
{
  if (!m_domain.contains(index)) return false;

  return m_bits[m_domain.index(index)];
}

DenseIntVectSet& DenseIntVectSet::operator-=(const Box& b)
{
  if (isEmpty()) return *this;
  if (!m_domain.intersects(b)) return *this;
  BoxIterator it(b & m_domain);
  for (it.begin(); it.ok(); ++it) this->operator-=(it());
  return *this;
}

DenseIntVectSet& DenseIntVectSet::operator|=(const DenseIntVectSet& d)
{
  if (m_domain.contains(d.m_domain))
        {
          BoxIterator bit(d.m_domain);
          int i=0;
          for (bit.begin(); bit.ok(); ++bit, ++i)
          {
                if (d.m_bits[i]) m_bits.setTrue(m_domain.index(bit()));
          }
        }
  else if (d.m_domain.contains(m_domain))
        {
          DenseIntVectSet newSet = d;
          BoxIterator bit(m_domain);
          int i=0;
          for (bit.begin(); bit.ok(); ++bit, ++i)
          {
                if (m_bits[i]) newSet.m_bits.setTrue(newSet.m_domain.index(bit()));
          }
          *this = newSet;
        }
  else
        {
          Box newDomain = minBox(m_domain, d.m_domain);
          DenseIntVectSet newSet(newDomain, false);
          newSet |= *this;
          newSet |= d;
          *this = newSet;
        }
  return *this;
}

DenseIntVectSet& DenseIntVectSet::operator|=(const Box& b)
{
  //if (!m_domain.contains(b))
  //{
  //      MayDay::Error("Box union with DenseIntVectSet outside m_domain");
  //}
  CH_assert(m_domain.contains(b));

  BoxIterator bit(b);
  for (bit.begin(); bit.ok(); ++bit)
        {
          m_bits.setTrue(m_domain.index(bit()));
        }
  return *this;
}

DenseIntVectSet& DenseIntVectSet::operator|=(const IntVect& intvect)
{
  //if (!m_domain.contains(intvect))
  //{
  //      MayDay::Error("union with DenseIntVectSet outside m_domain");
  //}
  CH_assert(m_domain.contains(intvect));

  m_bits.setTrue(m_domain.index(intvect));
  return *this;
}

DenseIntVectSet& DenseIntVectSet::operator&=(const ProblemDomain& domain)
{
  if (domain.domainBox().contains(m_domain))
  {
    return *this;
  }

  IntVect chopSml = domain.domainBox().smallEnd();
  IntVect chopBig = domain.domainBox().bigEnd();
  if (domain.isPeriodic())
    {
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          if (domain.isPeriodic(dir))
            {
              chopSml[dir] = m_domain.smallEnd(dir);
              chopBig[dir] = m_domain.bigEnd(dir);
            }
        }
    }
  Box chopper(chopSml, chopBig);
  *this &= chopper;

  return *this;
}

DenseIntVectSet& DenseIntVectSet::operator&=(const Box& b)
{
  if (b.contains(m_domain) || b == m_domain || isEmpty())
  {
    return *this;
  }

  if (!m_domain.intersects(b))
  {
    *this = DenseIntVectSet();
    return *this;
  }

  if (m_domain.numPts() < b.numPts())
  {
    int i;
    BoxIterator bit(m_domain);
    for (i = 0, bit.begin(); bit.ok(); ++i, ++bit)
    {
      if (!b.contains(bit()))
      {
        m_bits.setFalse(i);
      }
    }
  }
  else
  {
    Box btmp(b);
    btmp &= m_domain;

    DenseIntVectSet tmp = DenseIntVectSet(btmp, true);

    int i;
    BoxIterator bit(btmp);
    for (i = 0, bit.begin(); bit.ok(); ++i, ++bit)
    {
      if (!this->operator[](bit()))
      {
        tmp.m_bits.setFalse(i);
      }
    }

    makeEmptyBits();
    this->operator|=(tmp);
  }

  return *this;
}

Vector<Box> DenseIntVectSet::createBoxes() const
{
  Vector<Box> boxes;
  DenseIntVectSetIterator it(*this);
  for (it.begin(); it.ok(); ++it)
  {
        boxes.push_back(Box(it(), it()));
  }
  return boxes;
}

void DenseIntVectSet::recalcMinBox() const
{
  if (isEmpty())
  {
    (Box&)m_minBox = Box();
    return;
  }
  DenseIntVectSetIterator it(*this);
  it.begin();
  (Box&)m_minBox = Box(it(), it());
  ++it;
  for (; it.ok(); ++it)
  {
    Box m(it(), it());
    ((Box&)m_minBox).minBox(m);
  }
}

void DenseIntVectSet::makeEmptyBits()
{
  m_bits.setAllFalse();
}

DenseIntVectSet DenseIntVectSet::chop(int dir, int chop_pnt)
{
  if (m_minBox.smallEnd(dir) >= chop_pnt)
    {
      DenseIntVectSet rtn(*this);
      *this = DenseIntVectSet();
      return rtn;
    }
  if (m_minBox.bigEnd(dir) < chop_pnt)
    {
      return DenseIntVectSet();
    }
  Box chop = m_domain;
  Box chopThis = chop.chop(dir, chop_pnt);
  DenseIntVectSet left(chop);
  DenseIntVectSet right(chopThis);
  left.intersect(*this);
  right.intersect(*this);
  *this = right;
  return left;

}

bool DenseIntVectSet::contains(const Box& box) const
{
  if (!m_minBox.contains(box)) return false;
  for (BoxIterator bit(box); bit.ok(); ++bit)
    if (!this->operator[](bit())) return false;
  return true;
}

DenseIntVectSet& DenseIntVectSet::operator&=(const DenseIntVectSet& ivs)
{
  if (isEmpty()) return *this;
  if (ivs.isEmpty())
    {
      *this = DenseIntVectSet();
      return *this;
    }
  if (!m_domain.intersectsNotEmpty(ivs.m_domain))
    {
      *this = DenseIntVectSet();
      return *this;
    }
  return intersect(ivs);
}

DenseIntVectSet& DenseIntVectSet::intersect(const DenseIntVectSet& ivs)
{
  BoxIterator bit(m_domain);
  int i=0;
  for (bit.begin(); bit.ok();++i, ++bit)
    {
      if (m_bits[i])
        {
          if (ivs[bit()])
          {
            // do nothing
          }
          else
          {
            m_bits.setFalse(i);
          }
        }
    }
  return *this;
}

void DenseIntVectSet::coarsen(int iref)
{
  if (iref == 1) return;
  CH_assert(iref >= 1);
  // int refinements = iref/2;
  CH_assert((iref/2)*2 == iref); // check iref for power of 2

  Box newDomain(m_domain);
  newDomain.coarsen(iref);
  DenseIntVectSet newSet(newDomain, false);
  BoxIterator bit(m_domain);
  int count=0;
  for (bit.begin(); bit.ok(); ++bit, ++count)
    {
      if (m_bits[count])
        {
          IntVect iv(bit());
          iv.coarsen(iref);
          long index = newDomain.index(iv);
          newSet.m_bits.setTrue(index);
        }
    }

  *this = newSet;
}

void DenseIntVectSet::refine(int iref)
{
  if (iref == 1) return;
  if (isEmpty()) return;
  CH_assert(iref >= 1);
  //int refinements = iref/2;
  CH_assert((iref/2)*2 == iref); // check iref for power of 2
  Box newDomain(m_domain);
  newDomain.refine(iref);
  DenseIntVectSet newSet(newDomain, false);
  IntVect iv;
  BoxIterator bit(newDomain);
  int count=0;
  for (bit.begin(); bit.ok(); ++bit, ++count)
    {
      iv = bit();
      iv.coarsen(iref);
      if (this->operator[](iv))
        {
          newSet.m_bits.setTrue(count);
        }
    }
  *this = newSet;
}

void DenseIntVectSet::nestingRegion(int radius, const Box& domain)
{
  CH_assert(radius >= 0);
  if (radius == 0) return;

  DenseIntVectSet tmp(*this);

  {
    IntVect lo = m_domain.smallEnd();
    IntVect hi = m_domain.bigEnd();
    for (int i=0; i<SpaceDim; ++i)
      {
        if (lo[i] != domain.smallEnd()[i]) lo[i] += radius;
        if (hi[i] != domain.bigEnd()[i])   hi[i] -= radius;
      }
    Box shrink(lo, hi);
    *this &= shrink;
  }

  Box clobberBox(IntVect::Zero, IntVect::Zero);
  IntVect center(IntVect::Zero);
  clobberBox.grow(radius);
  BoxIterator bit(m_domain);
  int i=0;
  for (bit.begin(); bit.ok();++i, ++bit)
    {
      if (!(tmp.m_bits[i])) // was there a zero at this bit ?
        {
          if (domain.contains(bit()))
            {
              clobberBox.shift(bit()-center);
              *this -= clobberBox;
              center=bit();
            }
        }
    }
}

void DenseIntVectSet::nestingRegion(int radius, const ProblemDomain& a_domain)
{
  CH_assert(radius >= 0);
  if (radius == 0) return;

  DenseIntVectSet tmp;
  Box region = a_domain.domainBox();
  if (!a_domain.isPeriodic())
    {
      tmp = *this;
    }
  else
    {
      D_TERM6(if (a_domain.isPeriodic(0)) region.grow(0, radius);,
              if (a_domain.isPeriodic(1)) region.grow(1, radius);,
              if (a_domain.isPeriodic(2)) region.grow(2, radius);,
              if (a_domain.isPeriodic(3)) region.grow(3, radius);,
              if (a_domain.isPeriodic(4)) region.grow(4, radius);,
              if (a_domain.isPeriodic(5)) region.grow(5, radius);)
        tmp = DenseIntVectSet(region, false);
      tmp |= *this;
      for (int i=0; i<SpaceDim; ++i)
        {
          if (a_domain.isPeriodic(i))
            {
              int size = a_domain.domainBox().size(i);
              Box hiRegion = adjCellHi(a_domain.domainBox(), i, radius);
              for (BoxIterator bit(hiRegion); bit.ok(); ++bit)
                {
                  IntVect image = bit();
                  image[i] -= size;
                  if (this->operator[](image))
                    {
                      tmp |= bit();
                    }
                  image[i] += (size - radius);
                  if (this->operator[](image))
                    {
                      image[i] -= size;
                      tmp |= image;
                    }
                }
            }
        }
    }

  {
    IntVect lo = m_domain.smallEnd();
    IntVect hi = m_domain.bigEnd();
    for (int i=0; i<SpaceDim; ++i)
      {
        if (lo[i] != region.smallEnd()[i]) lo[i] += radius;
        if (hi[i] != region.bigEnd()[i])   hi[i] -= radius;

      }
    Box shrink(lo, hi);
    *this &= shrink;
  }

  Box clobberBox(IntVect::Zero, IntVect::Zero);
  IntVect center(IntVect::Zero);
  clobberBox.grow(radius);
  BoxIterator bit(region);
  int i=0;
  for (bit.begin(); bit.ok();++i, ++bit)
    {
      if (!(tmp.m_bits[i])) // was there a zero at this bit ?
        {
                  clobberBox.shift(bit()-center);
                  *this -= clobberBox;
                  center=bit();
        }
    }
}

bool DenseIntVectSet::isHighEmpty() const
{
  return isHighEmpty(0, SpaceDim-1);
}

bool DenseIntVectSet::isHighEmpty(const int a_dir) const
{
  return isHighEmpty(a_dir, a_dir);
}

bool DenseIntVectSet::isHighEmpty(const int a_dirS, const int a_dirE) const
{
  CH_assert(a_dirS >= 0 && a_dirE < SpaceDim && a_dirS <= a_dirE);
  // Stride for each direction index
  IntVect stride;
  D_EXPR6(stride[0] = 1,
          stride[1] = stride[0]*m_domain.size(0),
          stride[2] = stride[1]*m_domain.size(1),
          stride[3] = stride[2]*m_domain.size(2),
          stride[4] = stride[3]*m_domain.size(3),
          stride[5] = stride[4]*m_domain.size(4));

  // Loop over each hyperface on the high side of 'dir'
  for (int dir = a_dirS; dir <= a_dirE; ++dir)
    {
      // Start locations
      IntVect is(D_DECL6(0, 0, 0, 0, 0, 0));
      // In direction 'dir', only iterate once
      is[dir] = m_domain.size(dir) - 1;

      // The index is really
      //   index = il0*stride[0]  // = index[0]
      //         + il1*stride[1]  // = index[1]
      //         + il2*stride[2]  // = ...
      //         + il3*stride[3]
      //         + il4*stride[4]
      //         + il5*stride[5]);
      // But instead of multiplying throughout every time, just increment the
      // appropriate index (shown in comments on the right) with the
      // corresponding stride for each iteration.  When a loop is initialized,
      // it takes the index from the nearest outer loop plus the contribution
      // from its start value
      IntVect index;

      // We want to loop over directions SpaceDim->0 but space macros go the
      // other way.  Use reverse indices for the loop macros.
      D_TERM6(const int ir0 = SpaceDim - 1;,
              const int ir1 = SpaceDim - 2;,
              const int ir2 = SpaceDim - 3;,
              const int ir3 = SpaceDim - 4;,
              const int ir4 = SpaceDim - 5;,
              const int ir5 = SpaceDim - 6;)
      D_TERM6(int il0;, int il1;, int il2;, int il3;, int il4;, int il5;)
      // So within the loop, must pair the loop index 'il#' with 'ir#'.  Note
      // il# is used to control the loop and 'index' is being calculated.
      D_TERM6(for (il0 = is[ir0], index[ir0] =          0 + is[ir0]*stride[ir0]; il0 < m_domain.size(ir0); ++il0, index[ir0] += stride[ir0]),
              for (il1 = is[ir1], index[ir1] = index[ir0] + is[ir1]*stride[ir1]; il1 < m_domain.size(ir1); ++il1, index[ir1] += stride[ir1]),
              for (il2 = is[ir2], index[ir2] = index[ir1] + is[ir2]*stride[ir2]; il2 < m_domain.size(ir2); ++il2, index[ir2] += stride[ir2]),
              for (il3 = is[ir3], index[ir3] = index[ir2] + is[ir3]*stride[ir3]; il3 < m_domain.size(ir3); ++il3, index[ir3] += stride[ir3]),
              for (il4 = is[ir4], index[ir4] = index[ir3] + is[ir4]*stride[ir4]; il4 < m_domain.size(ir4); ++il4, index[ir4] += stride[ir4]),
              for (il5 = is[ir5], index[ir5] = index[ir4] + is[ir5]*stride[ir5]; il5 < m_domain.size(ir5); ++il5, index[ir5] += stride[ir5]))
      {
        // Note that index is still stored normally so 0 has the total value.
        if (m_bits[index[0]]) return false;
      }
    }
  return true;
}

void DenseIntVectSet::grow(int igrow)
{
  if (igrow >= 1)
    {
      IntVect range(IntVect::Unit);
      range*=igrow;
      grow(range);
    }
  else if (igrow < 0)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          grow(idir, igrow);
        }
    }
  else
    {
      //nothing to do if igrow == 0
      return;
    }

}

void DenseIntVectSet::grow(int idir, int igrow)
{
  CH_assert(idir >= 0);
  CH_assert(idir < SpaceDim);
  if (igrow >= 1)
    {
      IntVect range(IntVect::Zero);
      range[idir] = igrow;
      grow(range);
    }
  else if (igrow < 0)
    {
      IntVect shiftvec= igrow*BASISV(idir);
      DenseIntVectSet ivsShiftPlus = *this;
      DenseIntVectSet ivsShiftMinu = *this;
      ivsShiftPlus.shift(shiftvec);
      ivsShiftMinu.shift(-shiftvec);
      *this &= ivsShiftPlus;
      *this &= ivsShiftMinu;
    }
  else
    {
      //nothing to do if igrow == 0
      return;
    }

}

void DenseIntVectSet::grow(const IntVect& range)
{
  Box newDomain(m_domain);
  newDomain.grow(range);
  DenseIntVectSet newSet(newDomain, false);

  int index = 0;
#if CH_SPACEDIM == 1
  for (int i=0; i<m_domain.size(0); ++i, ++index)
    {
      if (m_bits[index])
        {
          int indexNew = i;
          for (int ii = 0; ii<2*range[0] + 1; ++ii, ++indexNew)
            newSet.m_bits.setTrue(indexNew);
        }
    }
#endif

#if CH_SPACEDIM == 2

  for (int j = 0; j<m_domain.size(1); ++j)
    {
      for (int i=0; i<m_domain.size(0); ++i, ++index)
        {
          if (m_bits[index])
            {
              for (int jj=0; jj<2*range[1] +1 ; ++jj)
                {
                  int indexNew = i + newDomain.size(0)*(j+jj);
                  for (int ii = 0; ii<2*range[0] + 1; ++ii, ++indexNew)
                    newSet.m_bits.setTrue(indexNew);
                }
            }
        }
    }
#endif

#if CH_SPACEDIM == 3
  for (int k=0; k<m_domain.size(2); ++k)
    for (int j = 0; j<m_domain.size(1); ++j)
      for (int i = 0; i<m_domain.size(0); ++i, ++index)
        {
          if (m_bits[index])
            {
              for (int kk=0; kk<2*range[2]+1; ++kk)
                for (int jj=0; jj<2*range[1] +1 ; ++jj)
                  {
                    int indexNew = i+
                      newDomain.size(0)*(j+jj) +
                      newDomain.size(1)*newDomain.size(0)*(k+kk);
                    for (int ii = 0; ii<2*range[0] + 1; ++ii, ++indexNew)
                      newSet.m_bits.setTrue(indexNew);
                  }
            }
        }

#endif

#if CH_SPACEDIM == 4
  for (int u=0; u<m_domain.size(3); u++)
    for (int k=0; k<m_domain.size(2); ++k)
      for (int j = 0; j<m_domain.size(1); ++j)
        for (int i = 0; i<m_domain.size(0); ++i, ++index)
          {
            if (m_bits[index])
              {
                for (int uu=0; uu<2*range[3]+1; ++uu)
                  for (int kk=0; kk<2*range[2]+1; ++kk)
                    for (int jj=0; jj<2*range[1] +1 ; ++jj)
                      {
                        int indexNew = i+
                          newDomain.size(0)*(j+jj) +
                          newDomain.size(0)*newDomain.size(1)*(k+kk) +
                          newDomain.size(0)*newDomain.size(1)*newDomain.size(2)*(u+uu);
                        for (int ii = 0; ii<2*range[0] + 1; ++ii, ++indexNew)
                          newSet.m_bits.setTrue(indexNew);
                      }
              }
          }

#endif

#if CH_SPACEDIM == 5
  for (int v=0; v<m_domain.size(4); v++)
    for (int u=0; u<m_domain.size(3); u++)
      for (int k=0; k<m_domain.size(2); ++k)
        for (int j = 0; j<m_domain.size(1); ++j)
          for (int i = 0; i<m_domain.size(0); ++i, ++index)
            {
              if (m_bits[index])
                {
                  for (int vv=0; vv<2*range[4]+1; ++vv)
                    for (int uu=0; uu<2*range[3]+1; ++uu)
                      for (int kk=0; kk<2*range[2]+1; ++kk)
                        for (int jj=0; jj<2*range[1] +1 ; ++jj)
                          {
                            int indexNew = i+
                              newDomain.size(0)*(j+jj) +
                              newDomain.size(0)*newDomain.size(1)*(k+kk) +
                              newDomain.size(0)*newDomain.size(1)*newDomain.size(2)*(u+uu) +
                              newDomain.size(0)*newDomain.size(1)*newDomain.size(2)*newDomain.size(3)*(v+vv);
                            for (int ii = 0; ii<2*range[0] + 1; ++ii, ++indexNew)
                              newSet.m_bits.setTrue(indexNew);
                          }
                }
            }

#endif

  // this is starting to look pretty ugly..
#if CH_SPACEDIM == 6
  for (int w=0; w<m_domain.size(5); w++)
    for (int v=0; v<m_domain.size(4); v++)
      for (int u=0; u<m_domain.size(3); u++)
        for (int k=0; k<m_domain.size(2); ++k)
          for (int j = 0; j<m_domain.size(1); ++j)
            for (int i = 0; i<m_domain.size(0); ++i, ++index)
              {
                if (m_bits[index])
                  {
                    for (int ww=0; ww<2*range[5]+1; ++ww)
                      for (int vv=0; vv<2*range[4]+1; ++vv)
                        for (int uu=0; uu<2*range[3]+1; ++uu)
                          for (int kk=0; kk<2*range[2]+1; ++kk)
                            for (int jj=0; jj<2*range[1] +1 ; ++jj)
                              {
                                int indexNew = i+
                                  newDomain.size(0)*(j+jj) +
                                  newDomain.size(0)*newDomain.size(1)*(k+kk) +
                                  newDomain.size(0)*newDomain.size(1)*newDomain.size(2)*(u+uu) +
                                  newDomain.size(0)*newDomain.size(1)*newDomain.size(2)*newDomain.size(3)*(v+vv) +
                                  newDomain.size(0)*newDomain.size(1)*newDomain.size(2)*newDomain.size(3)*newDomain.size(4)*(v+vv);
                                for (int ii = 0; ii<2*range[0] + 1; ++ii, ++indexNew)
                                  newSet.m_bits.setTrue(indexNew);
                              }
                  }
              }

#endif

#if (CH_SPACEDIM > 6)
  MayDay::Error("DenseIntVectSet::grow undefined for DIM>3");
#endif

  *this = newSet;
}

void DenseIntVectSet::growHi()
{
  if (!isHighEmpty())
    {
      Box newDomain(m_domain.smallEnd(), m_domain.bigEnd() + 1);
      DenseIntVectSet newSet(newDomain, false);
      newSet |= *this;
      *this = newSet;
    }

  // Stride for each direction index
  IntVect stride;
  D_EXPR6(stride[0] = 1,
          stride[1] = stride[0]*m_domain.size(0),
          stride[2] = stride[1]*m_domain.size(1),
          stride[3] = stride[2]*m_domain.size(2),
          stride[4] = stride[3]*m_domain.size(3),
          stride[5] = stride[4]*m_domain.size(4));

  // We can skip down through a layer of 0 bits (in the SpaceDim-1 direction)
  // before starting (this is m_domain.size - stride for last index)
  int start = m_domain.size().product() -
    m_domain.size().product()/m_domain.size(SpaceDim-1);
  BitSetIterator itb(m_bits);
  itb.setpos(start-1);
  // 'i' is used to control the loop
  for (int i = start; i--; --itb)
    {
      if (itb())
        {
          // Set a box going from i to i+1 in each direction.  The macro nests
          // the loops in column-order but this is really scattered anyways.
          // See DenseIntVectSet::isHighEmpty(int, int) for more information on
          // how 'index' is updated.
          IntVect index;
          D_TERM6(int il0;, int il1;, int il2;, int il3;, int il4;, int il5;)
          D_TERM6(for (il0 = 0, index[0] =        0; il0 <= 1; ++il0, index[0] += stride[0]),
                  for (il1 = 0, index[1] = index[0]; il1 <= 1; ++il1, index[1] += stride[1]),
                  for (il2 = 0, index[2] = index[1]; il2 <= 1; ++il2, index[2] += stride[2]),
                  for (il3 = 0, index[3] = index[2]; il3 <= 1; ++il3, index[3] += stride[3]),
                  for (il4 = 0, index[4] = index[3]; il4 <= 1; ++il4, index[4] += stride[4]),
                  for (il5 = 0, index[5] = index[4]; il5 <= 1; ++il5, index[5] += stride[5]))
            {
              m_bits.setTrue(i + index[SpaceDim-1]);
            }
        }
    }
}

void DenseIntVectSet::growHi(const int a_dir)
{
  if (!isHighEmpty(a_dir))
    {
      Box newDomain(m_domain);
      newDomain.growHi(a_dir);
      DenseIntVectSet newSet(newDomain, false);
      newSet |= *this;
      *this = newSet;
    }
  const int offset = D_TERM6(1,
                             *std::max(1, (1 <= a_dir)*m_domain.size(0)),
                             *std::max(1, (2 <= a_dir)*m_domain.size(1)),
                             *std::max(1, (3 <= a_dir)*m_domain.size(2)),
                             *std::max(1, (4 <= a_dir)*m_domain.size(3)),
                             *std::max(1, (5 <= a_dir)*m_domain.size(4)));
  // We can skip down by 'offset' before starting which is a layer of 0 bits in
  // direction 'a_dir'
  int start = m_domain.size().product() - offset;
  BitSetIterator itb(m_bits);
  itb.setpos(start-1);
  // 'i' is used to control the loop
  for (int i = start; i--; --itb)
    {
      if (itb())
        {
          // Set the bit at i+1 in direction 'a_dir'
          m_bits.setTrue(i + offset);
        }
    }
}

DenseIntVectSet& DenseIntVectSet::operator-=(const DenseIntVectSet& ivs)
{
  BoxIterator bit(m_domain);
  int i=0;
  for (bit.begin(); bit.ok();++i, ++bit)
    {
      if (m_bits[i])
        {
          if (ivs[bit()])
            m_bits.setFalse(i);
        }
    }
  return *this;
}

bool DenseIntVectSet::isEmpty() const
{
  return m_bits.isEmpty();
}
bool DenseIntVectSet::isFull() const
{
  return m_bits.isFull();
}

int DenseIntVectSet::numPts() const
{
  int n = m_domain.numPts();
  int c=0;
  for (int i=0; i<n; ++i)
    if (m_bits[i]) ++c;
  return c;
}

bool DenseIntVectSet::operator==(const DenseIntVectSet& a_lhs) const
{
  if (numPts() != a_lhs.numPts()) return false;
  DenseIntVectSetIterator it(*this);
  for (; it.ok(); ++it)
    {
      if (!a_lhs[it()]) return false;
    }
  return true;
}

bool DenseIntVectSet::operator<(const DenseIntVectSet& a_ivs) const
{
  // Primary criterion: Box::operator<() applied to m_domain.
  // Secondary criterion: BitSet::operator<() applied to m_bits.

  if ( m_domain < a_ivs.m_domain )
  {
    return true;
  } else
  if ( a_ivs.m_domain < m_domain )
  {
    return false;
  }

  if ( m_bits < a_ivs.m_bits )
  {
    return true;
  }

  return false;
}

int DenseIntVectSet::linearSize() const
{
  if (isEmpty())
  {
    return CH_XD::linearSize<int>(0);
  }
  return m_bits.linearSize() + 2 * CH_XD::linearSize<Box>(m_minBox);
}

void DenseIntVectSet::linearIn(const void* const inBuf)
{
  static DenseIntVectSet emptySet;
  *this = emptySet;
  int* b = (int*)inBuf;
  if (*b == 0) return;

  char* buf = (char*)inBuf;
  m_bits.linearIn(buf);
  buf+=m_bits.linearSize();
  CH_XD::linearIn<Box>(m_domain, buf);
  buf+=CH_XD::linearSize<Box>(m_domain);
  CH_XD::linearIn<Box>(m_minBox, buf);
}

void DenseIntVectSet::linearOut(void* const a_outBuf) const
{
  if (isEmpty())
    {
      int* b = (int*)a_outBuf;
      *b = 0;
      return;
    }
  char* buf = (char*)a_outBuf;
  m_bits.linearOut(buf);
  buf+=m_bits.linearSize();
  CH_XD::linearOut<Box>(buf, m_domain);
  buf+=CH_XD::linearSize<Box>(m_domain);
  CH_XD::linearOut<Box>(buf, m_minBox);
}

void DenseIntVectSet::compact() const
{
  DenseIntVectSet* nthis = (DenseIntVectSet*)this;
  if (isEmpty())
  {
    *nthis = DenseIntVectSetIterator::emptyDenseIntVectSet;
    return;
  } else if (isFull())
  {
    return;
  }
  DenseIntVectSetIterator it(*this);
  Box nDomain(it(), it());
  IntVect& lo = (IntVect&)(nDomain.smallEnd());
  IntVect& hi = (IntVect&)(nDomain.bigEnd());
  for (++it; it.ok(); ++it)
    {
      lo.min(it());
      hi.max(it());
    }
  nDomain.computeBoxLenNotEmpty();
  DenseIntVectSet nset(nDomain);
  nset.intersect(*this);
  *nthis = nset;
}

DenseIntVectSet  DenseIntVectSetIterator::emptyDenseIntVectSet;

void
DenseIntVectSetIterator::thisIntVect(const int a_linearPos)
{
  if (a_linearPos - m_prevLinearPos == 1)
    {
      D_TERM6(
        ++(m_current[0]);,

        const IntVect& big = m_ivsPtr->m_domain.bigEnd();
        const IntVect& small = m_ivsPtr->m_domain.smallEnd();
        if (m_current[0] > big[0])
          {
            m_current[0] = small[0];
            ++m_current[1];
          },

        if (m_current[1] > big[1])
          {
            m_current[1] = small[1];
            ++m_current[2];
          },

        if (m_current[2] > big[2])
          {
            m_current[2] = small[2];
            ++m_current[3];
          },

        if (m_current[3] > big[3])
          {
            m_current[3] = small[3];
            ++m_current[4];
          },

        if (m_current[4] > big[4])
          {
            m_current[4] = small[4];
            ++m_current[5];
          }
      )
    }
  else
    {
      const IntVect& small = m_ivsPtr->m_domain.smallEnd();
      D_TERM6(
        const int d0 = SpaceDim - 1;
        const int p0 = a_linearPos/m_stride[d0];
        m_current[d0] = small[d0] + p0;,

        int used = p0*m_stride[d0];
        const int d1 = SpaceDim - 2;
        const int p1 = (a_linearPos - used)/m_stride[d1];
        m_current[d1] = small[d1] + p1;,

        used += p1*m_stride[d1];
        const int d2 = SpaceDim - 3;
        const int p2 = (a_linearPos - used)/m_stride[d2];
        m_current[d2] = small[d2] + p2;,

        used += p2*m_stride[d2];
        const int d3 = SpaceDim - 4;
        const int p3 = (a_linearPos - used)/m_stride[d3];
        m_current[d3] = small[d3] + p3;,

        used += p3*m_stride[d3];
        const int d4 = SpaceDim - 5;
        const int p4 = (a_linearPos - used)/m_stride[d4];
        m_current[d4] = small[d4] + p4;,

        used += p4*m_stride[d4];
        const int d5 = SpaceDim - 6;
        const int p5 = (a_linearPos - used)/m_stride[d5];
        m_current[d5] = small[d5] + p5;
      )
    }
  m_prevLinearPos = a_linearPos;
}

#include "NamespaceFooter.H"
