#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "IntVectSet.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

long int IntVectSet::count = 0;
long int IntVectSet::peakcount = 0;
int      IntVectSet::s_maxDense = 6400000;

IntVectSet::~IntVectSet()
{
  count--;
}

void IntVectSet::define()
{
  m_ivs.clear();
  m_dense = DenseIntVectSet();
  m_isdense = true;
}

void IntVectSet::define(const DenseIntVectSet& a_dense)
{
  m_ivs.clear();
  m_dense = a_dense;
  m_isdense = true;
}

void IntVectSet::define(const TreeIntVectSet& a_tree)
{
  m_ivs = a_tree;
  m_dense = DenseIntVectSet();
  m_isdense = false;
}

IntVectSet::IntVectSet(const DenseIntVectSet& a_dense)
{
  count++;
  define(a_dense);
}

IntVectSet::IntVectSet(const TreeIntVectSet& a_tree)
{
  count++;
  define(a_tree);
}

IntVectSet::IntVectSet(const IntVect& iv_in)
{
  count++;
  if (count > peakcount) peakcount = count;
  define(iv_in);
}

void IntVectSet::define(const IntVect& iv_in)
{
  define(Box(iv_in, iv_in));
}

IntVectSet::IntVectSet(const Box& b)
{
  count++;
  if (count > peakcount) peakcount = count;
  define(b);
}

void
IntVectSet::setMaxDense(const int& a_maxDense)
{
  s_maxDense = a_maxDense;
}
void IntVectSet::define(const Box& b)
{
  if (b.numPts() < s_maxDense)
    {
      m_ivs.clear();
      m_isdense = true;
      m_dense = DenseIntVectSet(b);
    }
  else
    {
      m_isdense = false;
      m_ivs.define(b);
    }
}

void IntVectSet::define_boxCorners(const Box& b)
{
  const IntVect ends[2] =
  {
    b.smallEnd(),
    b.bigEnd()
  };
  int nCorners = 1;
  for (int d=0; d<SpaceDim; d++)
    nCorners *= 2;
  Vector<IntVect> corners(nCorners, IntVect::Zero);
  // each component of a corner IntVect has two choices: lo or hi.
  // For all test IntVects, split them into two groups for dimension d.
  //  different d leads to different ways of splitting,
  //  which is realized by alternating remainders of (i/stride)
  int stride = nCorners/2;
  for (int d=0; d<SpaceDim; d++)
    {
      for (int i=0; i<corners.size(); i++)
        {
          const int endId = (i/stride)%2; // only 2 ends.
          corners[i][d] = ends[endId][d];
        }
      stride /= 2;
    }
  // add the results into this set
  for (int i=0; i<corners.size(); i++)
    this->operator|=(corners[i]);

  /// This is too slow due to the implementation of IntvectSet.
  // define(b);
  // Box testBox(b);
  // testBox.grow(-IntVect::Unit);
  // for (int d=0; d<SpaceDim; d++)
  //   {
  //     testBox.grow(+BASISV(d));
  //     this->operator-=(testBox);
  //     testBox.grow(-BASISV(d));
  //   }
}

IntVectSet operator|(const IntVectSet& ivs1, const IntVectSet& ivs2)
{
  IntVectSet rtn = ivs1;
  rtn |= ivs2;
  return rtn;
}

IntVectSet operator|(const IntVectSet& ivs, const IntVect& iv)
{
  IntVectSet rtn = ivs;
  rtn |= iv;
  return rtn;
}

IntVectSet operator|(const IntVect& iv, const IntVectSet& ivs)
{
  IntVectSet rtn = ivs;
  rtn |= iv;
  return rtn;
}

IntVectSet operator|(const IntVectSet& ivs, const Box& b)
{
  IntVectSet rtn = ivs;
  rtn |= b;
  return rtn;
}

IntVectSet operator|(const Box& b, const IntVectSet& ivs)
{
  IntVectSet rtn = ivs;
  rtn |= b;
  return rtn;
}

IntVectSet& IntVectSet::operator|=(const IntVect& iv)
{

  if (m_isdense)
  {
        if (!m_dense.box().contains(iv))
        {
          convert();
          m_ivs |= iv;
        }
        else
        {
          m_dense |= iv;
        }
  }
  else
  {
        m_ivs |= iv;
  }
  return *this;
}

IntVectSet& IntVectSet::operator|=(const Box& b)
{
  if (m_isdense)
  {
        if (!m_dense.box().contains(b))
        {
          convert();
          m_ivs |= b;
        }
        else
        {
          m_dense |= b;
        }
  }
  else
  {
        m_ivs |= b;
  }
  return *this;

}

IntVectSet& IntVectSet::operator|=(const IntVectSet& ivs)
{
  if (m_isdense)
    {
      if (ivs.m_isdense)
        {
          if (CH_XD::minBox(m_dense.box(), ivs.m_dense.box()).numPts() <=
              std::max(m_dense.box().numPts(), ivs.m_dense.box().numPts()))
            {
              m_dense |= ivs.m_dense;
              return *this;
            }
          ivs.convert();
        }
      convert();
    }
  else if (ivs.m_isdense)
    {
      ivs.convert();
    }
  m_ivs |= ivs.m_ivs;
  return *this;
}

IntVectSet IntVectSet::operator-(const IntVectSet& ivs) const
{
  IntVectSet rtn = *this;
  rtn -= ivs;
  return rtn;
}

IntVectSet IntVectSet::operator-(const IntVect& iv) const
{
  IntVectSet rtn = *this;
  rtn -= iv;
  return rtn;
}

IntVectSet IntVectSet::operator-(const Box& b) const
{
  IntVectSet rtn = *this;
  rtn -= b;
  return rtn;
}

IntVectSet& IntVectSet::operator-=(const IntVectSet& ivs)
{
  if (!(minBox().intersects(ivs.minBox()))) return *this;
  if (m_isdense)
    {
      if (ivs.m_isdense) m_dense-=ivs.m_dense;
      else
        {
          for (TreeIntVectSetIterator it(ivs.m_ivs); it.ok(); ++it) m_dense -= it();
        }
    }
  else
    {
      if (ivs.m_isdense)
        {
          for (DenseIntVectSetIterator it(ivs.m_dense);it.ok(); ++it) m_ivs -= it();
        }
      else
        m_ivs -= ivs.m_ivs;
    }

  return *this;
}

bool IntVectSet::operator==(const IntVectSet& a_lhs) const
{
  if (m_isdense)
  {
    if (!a_lhs.m_isdense) return false;
    return m_dense == a_lhs.m_dense;
  }
  if (a_lhs.m_isdense) return false;
  return m_ivs == a_lhs.m_ivs;
}

bool IntVectSet::operator<(const IntVectSet& a_ivs) const
{
  if ( m_isdense )
  {
    if ( !a_ivs.m_isdense )
    {
      return true;
    } else
    {
      return m_dense < a_ivs.m_dense;
    }
  } else
  {
    if ( !a_ivs.m_isdense )
    {
      return m_ivs < a_ivs.m_ivs;
    } else
    {
      return false;
    }
  }
}

int IntVectSet::linearSize() const
{
  if (m_isdense) return m_dense.linearSize() + sizeof(int);
  return m_ivs.linearSize() + sizeof(int);
}

void IntVectSet::linearIn(const void* const inBuf)
{
  int* b = (int*)inBuf;
  char* buf = (char*)inBuf;
  buf+=sizeof(int);
  if (*b == 0)
  {
    m_isdense = true;
    m_dense.linearIn(buf);
  }
  else
  {
    m_isdense = false;
    m_ivs.linearIn(buf);
  }
}
void IntVectSet::linearOut(void* const a_outBuf) const
{
  int* b = (int*)a_outBuf;
  char* buf = (char*)a_outBuf;
  buf+=sizeof(int);
  if (m_isdense)
  {
    *b=0;
    m_dense.linearOut(buf);
  }
  else
  {
    *b=1;
    m_ivs.linearOut(buf);
  }
}

IntVectSet& IntVectSet::operator-=(const IntVect& iv)
{
  if (m_isdense) m_dense -= iv;
  else          m_ivs   -= iv;
  return *this;
}

IntVectSet& IntVectSet::operator-=(const Box& b)
{
  if (m_isdense) m_dense -= b;
  else          m_ivs   -= b;
  return *this;
}

IntVectSet operator&(const IntVectSet& ivs1, const IntVectSet& ivs2)
{
  IntVectSet rtn = ivs1;
  rtn &= ivs2;
  return rtn;
}

IntVectSet operator&(const IntVectSet& ivs, const Box& b)
{
  IntVectSet rtn = ivs;
  rtn &= b;
  return rtn;
}

IntVectSet operator&(const Box& b, const IntVectSet& ivs)
{
  IntVectSet rtn = ivs;
  rtn &= b;
  return rtn;
}

IntVectSet& IntVectSet::operator&=(const Box& b)
{
  if (m_isdense) m_dense &= b;
  else          m_ivs   &= b;
  return *this;
}

IntVectSet& IntVectSet::operator&=(const ProblemDomain& d)
{
  if (m_isdense) m_dense &= d;
  else          m_ivs   &= d;
  return *this;
}

IntVectSet& IntVectSet::operator&=(const IntVectSet& ivs)
{
  if (!(minBox().intersects(ivs.minBox())))
    {
      m_ivs.clear();
      m_isdense = true;
      m_dense = DenseIntVectSet();
      return *this;
    }
  if (m_isdense)
    {
      if (ivs.m_isdense) m_dense&=ivs.m_dense;
      else
        {
          convert();
          m_ivs &= ivs.m_ivs;
        }
    }
  else
    {
      if (ivs.m_isdense) ivs.convert();
        m_ivs &= ivs.m_ivs;
    }

  return *this;
}

IntVectSet grow(const IntVectSet& ivs, int igrow)
{
  IntVectSet rtn(ivs);
  rtn.grow(igrow);
  return rtn;
}

void IntVectSet::grow(int igrow)
{
  if (m_isdense) m_dense.grow(igrow);
  else          m_ivs.grow(igrow);
  //  return *this;
}

void IntVectSet::nestingRegion(int radius, const Box& domain, int granularity)
{
  if (m_isdense) m_dense.nestingRegion(radius, domain);
  else          m_ivs.nestingRegion(radius, domain, granularity);
}

void IntVectSet::nestingRegion(int radius, const ProblemDomain& domain, int granularity)
{
  if (m_isdense) m_dense.nestingRegion(radius, domain);
  else          m_ivs.nestingRegion(radius, domain, granularity);
}

IntVectSet& IntVectSet::grow(int idir, int igrow)
{
  CH_assert(idir >= 0);
  CH_assert(idir < SpaceDim);
  if (m_isdense) m_dense.grow(idir, igrow);
  else          m_ivs.grow(idir, igrow);
  return *this;
}

void IntVectSet::growHi()
{
  if (m_isdense) m_dense.growHi();
  else          m_ivs.growHi();
}

void IntVectSet::growHi(const int a_dir)
{
  if (m_isdense) m_dense.growHi(a_dir);
  else          m_ivs.growHi(a_dir);
}

IntVectSet refine(const IntVectSet& ivs, int iref)
{
  IntVectSet rtn(ivs);
  rtn.refine(iref);
  return rtn;
}

IntVectSet& IntVectSet::refine(int iref)
{
  if (m_isdense) m_dense.refine(iref);
  else          m_ivs.refine(iref);
  return *this;
}

IntVectSet coarsen(const IntVectSet& ivs, int iref)
{
  IntVectSet rtn(ivs);
  rtn.coarsen(iref);
  return rtn;
}

IntVectSet& IntVectSet::coarsen(int iref)
{
  if (m_isdense) m_dense.coarsen(iref);
  else          m_ivs.coarsen(iref);
  return *this;
}

void IntVectSet::shift(const IntVect& iv)
{
  if (m_isdense) m_dense.shift(iv);
  else          m_ivs.shift(iv);
}

void IntVectSet::makeEmpty()
{
  if (m_isdense) m_dense = DenseIntVectSet();
  else          m_ivs.clear();
}

void IntVectSet::makeEmptyBits()
{
  if (m_isdense) m_dense.makeEmptyBits();
  else          m_ivs.clear();
}

void IntVectSet::compact() const
{
  if (m_isdense) m_dense.compact();
  else          m_ivs.compact();
}

IntVectSet IntVectSet::chop(int dir, int chop_pnt)
{
  if (m_isdense)
  {
        DenseIntVectSet r = m_dense.chop(dir, chop_pnt);
        IntVectSet rtn(r);
        return rtn;
  }
  else
  {
        TreeIntVectSet t = m_ivs.chop(dir, chop_pnt);
        IntVectSet rtn(t);
        return rtn;
  }

}

void IntVectSet::chop(int dir, int chop_pnt, IntVectSet& a_hi)
{
  if (m_isdense)
    a_hi = chop(dir, chop_pnt);
  else
    {
      m_ivs.chop(dir, chop_pnt, a_hi.m_ivs);
      a_hi.m_isdense = false;
    }

}
const Box& IntVectSet::minBox() const
{
  if (m_isdense) return m_dense.mBox();
  m_ivs.recalcMinBox();
  return m_ivs.minBox();
}

void IntVectSet::recalcMinBox() const
{
   if (m_isdense)
     m_dense.recalcMinBox();
   else
     m_ivs.recalcMinBox();
}

bool IntVectSet::isEmpty() const
{
  if (m_isdense) return m_dense.isEmpty();
  return m_ivs.isEmpty();
}

int IntVectSet::numPts() const
{
  if (m_isdense) return m_dense.numPts();
  return m_ivs.numPts();
}

bool IntVectSet::contains(const IntVect& iv) const
{
  if (m_isdense) return m_dense[iv];
  return m_ivs.contains(iv);
}

bool IntVectSet::contains(const IntVectSet& ivs) const
{
  if (&ivs == this) return true;
  for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit)
    {
      if (!contains(ivsit()))
        return false;
    }
  //return true for empty input
  return true;
}

bool IntVectSet::contains(const Box& box) const
{
  if (m_isdense) return m_dense.contains(box);
  return m_ivs.contains(box);
}

Vector<Box> IntVectSet::boxes() const
{
  if (m_isdense) return m_dense.createBoxes();
  return m_ivs.createBoxes();
}

void IntVectSet::printBoxes(std::ostream& a_ostrm) const
{
  Vector<Box> b = boxes();
  for (int i=0; i<b.size(); ++i) a_ostrm<<b[i]<<std::endl;
}

std::ostream& operator<<(std::ostream& os, const IntVectSet& ivs)
{
  ivs.printBoxes(os);
  return os;
}

void IntVectSet::convert() const
{
  if (!m_isdense) return; //already converted
  if (m_dense.isEmpty())
  {
    // do nothing
  }
  else if (m_dense.isFull())
    {
      ((TreeIntVectSet&)m_ivs).define(m_dense.box());
    }
  else
    {
      DenseIntVectSetIterator it(m_dense);
      for (it.begin(); it.ok(); ++it)
        {
          ((TreeIntVectSet&)m_ivs) |= it();
        }
      m_ivs.compact();
    }
  (bool&)m_isdense = false;
}

// if you are in this function, you better really
// have returned all the TreeNodes to the pool before
// you call this function!!

void IntVectSet::clearStaticMemory()
{
  TreeIntVectSet::treeNodePool->clear();
  TreeIntVectSet::index.clear();
  TreeIntVectSet::parents.clear();
  TreeIntVectSet::boxes.clear();
  TreeIntVectSet::bufferOffset.clear();
}

//====================================================================
IVSIterator::IVSIterator(const IntVectSet& ivs)
{
  define(ivs);
}

void IVSIterator::define(const IntVectSet& ivs)
{
  if (ivs.m_isdense)
    {
      m_isdense = true;
      m_dense.define(ivs.m_dense);
    }
  else
    {
      m_isdense = false;
      m_tree.define(ivs.m_ivs);
    }
}
#include "NamespaceFooter.H"
