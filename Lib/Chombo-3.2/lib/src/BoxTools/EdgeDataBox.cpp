#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EdgeDataBox.H"

#include "NamespaceHeader.H"

// first do simple access functions
// ---------------------------------------------------------
int
EdgeDataBox::nComp() const
{
  return m_nvar;
}

// ---------------------------------------------------------
const Box&
EdgeDataBox::box() const
{
  return m_bx;
}

// ---------------------------------------------------------
FArrayBox&
EdgeDataBox::getData(const int dir)
{
  CH_assert(m_nvar >0);
  CH_assert(dir < SpaceDim);
  CH_assert(m_data[dir] != NULL);

  return *m_data[dir];
}

// ---------------------------------------------------------
const FArrayBox&
EdgeDataBox::getData(const int dir) const
{
  CH_assert(m_nvar >0);
  CH_assert(dir < SpaceDim);
  CH_assert(m_data[dir] != NULL);

  return *m_data[dir];
}

// ---------------------------------------------------------
FArrayBox&
EdgeDataBox::operator[] (const int dir)
{
  CH_assert(m_nvar >0);
  CH_assert(dir < SpaceDim);
  CH_assert(m_data[dir] != NULL);

  return *m_data[dir];
}

// ---------------------------------------------------------
const FArrayBox&
EdgeDataBox::operator[] (const int dir)  const
{
  CH_assert(m_nvar >0);
  CH_assert(dir < SpaceDim);
  CH_assert(m_data[dir] != NULL);

  return *m_data[dir];
}

// ---------------------------------------------------------
// constructors and destructors
// ---------------------------------------------------------
EdgeDataBox::EdgeDataBox() : m_data(SpaceDim,NULL)
{
  m_nvar = -1;
}

// ---------------------------------------------------------
EdgeDataBox::EdgeDataBox(const Box& a_bx, int a_nComp)
  : m_data(SpaceDim,NULL)
{
  define(a_bx, a_nComp);
}

// ---------------------------------------------------------
EdgeDataBox::~EdgeDataBox()
{

  clear();
}

// ---------------------------------------------------------
void
EdgeDataBox::clear()
{
  // first delete storage
  for (int dir=0; dir < SpaceDim; dir++)
  {
    if (m_data[dir] != NULL)
    {
      delete m_data[dir];
      m_data[dir] = NULL;
    }
  }

  // now reset all other variables
  m_nvar = -1;
  // set the box to the empty box...
  m_bx = Box();
}

// ---------------------------------------------------------
// define function
void
EdgeDataBox::define(const Box& a_bx, int a_nComp)
{
  CH_assert(a_nComp > 0);

  m_bx = a_bx;
  m_nvar = a_nComp;

  m_data.resize(SpaceDim);

  for (int dir=0; dir<SpaceDim; dir++)
  {
    Box edgeBox(surroundingNodes(m_bx));
    edgeBox.enclosedCells(dir);
    FArrayBox* newFabPtr = new FArrayBox(edgeBox,m_nvar);
    m_data[dir] = newFabPtr;
  }

}

// ---------------------------------------------------------
// should resize fluxes in space (could be faster than re-allocating
// storage)
void
EdgeDataBox::resize(const Box& a_bx, int a_nComp)
{
  // if this object has not already been defined, call define fn.
  if (m_nvar < 0)
  {
    define(a_bx, a_nComp);
  }
  else
  {
    CH_assert(a_nComp > 0);
    m_bx = a_bx;
    m_nvar = a_nComp;

    for (int dir=0; dir<SpaceDim; dir++)
    {
      Box edgeBox(surroundingNodes(m_bx));
      edgeBox.enclosedCells(dir);

      if (m_data[dir] != NULL)
      {
        m_data[dir]->resize(edgeBox,m_nvar);
      }
      else
      {
        FArrayBox* newFabPtr = new FArrayBox(edgeBox, m_nvar);
        m_data[dir] = newFabPtr;
      }
    }

  }

}

// ---------------------------------------------------------
void
EdgeDataBox::setVal(const Real val)
{
  CH_assert(m_nvar > 0);

  for (int dir = 0; dir < SpaceDim; dir++)
  {
    CH_assert(m_data[dir] != NULL);
    m_data[dir]->setVal(val);
  }
}

// ---------------------------------------------------------
void
EdgeDataBox::setVal(const Real val, int dir)
{
  CH_assert(dir < SpaceDim);

  CH_assert(m_data[dir] != NULL);
  m_data[dir]->setVal(val);
}

// ---------------------------------------------------------
void
EdgeDataBox::setVal(const Real val, const int dir, const int startComp,
                    const int nComp)
{
  CH_assert(startComp >-1);
  CH_assert(startComp + nComp < m_nvar);
  CH_assert(dir < SpaceDim);
  CH_assert(m_data[dir] != NULL);

  for (int comp=startComp; comp < startComp+nComp; comp++)
  {
    m_data[dir]->setVal(val,comp);
  }

}

// ---------------------------------------------------------
void
EdgeDataBox::setVal(const Real val, const Box& bx)
{
  CH_assert(m_bx.contains(bx));

  for (int dir = 0; dir<SpaceDim; dir++)
  {
    CH_assert(m_data[dir] != NULL);
    // move cell-centered box to appropriate edge
    Box edgeBox(surroundingNodes(bx));
    edgeBox.enclosedCells(dir);
    m_data[dir]->setVal(val,edgeBox,0,m_nvar);
  }

}

// ---------------------------------------------------------
void
EdgeDataBox::setVal(const Real val, const Box& bx, const int dir,
                const int startComp, const int nComp)
{
  CH_assert(m_bx.contains(bx));
  CH_assert(m_data[dir] != NULL);
  CH_assert(startComp > -1);
  CH_assert(startComp +nComp <=m_nvar);

  Box edgeBox(surroundingNodes(bx));
  edgeBox.enclosedCells(dir);
  m_data[dir]->setVal(val, edgeBox, startComp, nComp);
}

// ---------------------------------------------------------
void
EdgeDataBox::copy(const EdgeDataBox& src)
{
  CH_assert(src.box() == m_bx);
  CH_assert(src.nComp() == m_nvar);

  for (int dir=0; dir < SpaceDim; dir++)
  {
    CH_assert(m_data[dir] != NULL);
    m_data[dir]->copy(src[dir]);
  }

}

// ---------------------------------------------------------
void
EdgeDataBox::copy(const EdgeDataBox& src, const int srcComp,
                  const int destComp, const int numComp)
{
  // to ensure that neither comp is negative
  CH_assert(srcComp*destComp > -1);
  CH_assert(srcComp+numComp <= src.nComp());
  CH_assert(destComp+numComp <= m_nvar);

  for (int dir=0; dir<SpaceDim; dir++)
  {
    CH_assert(m_data[dir] != 0);
    const FArrayBox& srcFab = src[dir];
    m_data[dir]->copy(srcFab,srcComp, destComp, numComp);
  }
}

// ---------------------------------------------------------
void
EdgeDataBox::copy(const EdgeDataBox& src, const int dir,
                  const int srcComp, const int destComp,
                  const int numComp)
{
  // to ensure that neither comp is negative
  CH_assert(srcComp*destComp > -1);
  CH_assert(srcComp+numComp <= src.nComp());
  CH_assert(destComp+numComp <= m_nvar);
  CH_assert(dir < SpaceDim);
  CH_assert(dir > -1);
  CH_assert(m_data[dir] != NULL);

  const FArrayBox& srcFab = src[dir];
  m_data[dir]->copy(srcFab,srcComp, destComp, numComp);
}

// ---------------------------------------------------------
void
EdgeDataBox::copy(const Box& RegionFrom,
                  const Interval& Cdest,
                  const Box& RegionTo,
                  const EdgeDataBox& src,
                  const Interval& Csrc)
{

  for (int dir=0; dir<SpaceDim; dir++)
  {
    CH_assert(m_data[dir] != NULL);
    Box RedgeFrom(RegionFrom);
    RedgeFrom.surroundingNodes();
    RedgeFrom.enclosedCells(dir);
    Box RedgeTo(RegionTo);
    RedgeTo.surroundingNodes();
    RedgeTo.enclosedCells(dir);
    const FArrayBox& srcFab = src[dir];
    m_data[dir]->copy(RedgeFrom,Cdest,RedgeTo,srcFab,Csrc);
  }

}

// ---------------------------------------------------------
int
EdgeDataBox::size(const Box& bx, const Interval& comps) const
{
  int totalSize = 0;

  FArrayBox tempFab;
  for (int dir =0; dir<SpaceDim; dir++)
  {
    // why did i put this line in???
    //    CH_assert(m_fluxes[dir] != NULL);
    Box edgeBox(surroundingNodes(bx));
    edgeBox.enclosedCells(dir);
    const int dirSize = tempFab.size(edgeBox,comps);
    totalSize += dirSize;
  }

  return totalSize;
}

// ---------------------------------------------------------
void
EdgeDataBox::linearOut(void* buf, const Box& R, const Interval& comps) const
{
  Real* buffer = (Real*) buf;
  for (int dir=0; dir<SpaceDim; dir++)
  {
    CH_assert(m_data[dir] != NULL);
    Box dirBox(surroundingNodes(R));
    dirBox.enclosedCells(dir);
    int dirSize = m_data[dir]->size(dirBox, comps);
    m_data[dir]->linearOut(buffer, dirBox, comps);
    buffer += dirSize/sizeof(Real);
  }
}

// ---------------------------------------------------------
void
EdgeDataBox::linearIn(void* buf, const Box& R, const Interval& comps)
{
  Real* buffer = (Real*) buf;
  for (int dir=0; dir<SpaceDim; dir++)
  {
    CH_assert(m_data[dir] != NULL);
    Box dirBox(surroundingNodes(R));
    dirBox.enclosedCells(dir);
    int dirSize = m_data[dir]->size(dirBox, comps);
    m_data[dir]->linearIn(buffer,dirBox,comps);
    buffer += dirSize/sizeof(Real);
  }
}

#include "NamespaceFooter.H"
