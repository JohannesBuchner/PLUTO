#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxIterator.H"
#include "RefCountedPtr.H"

#include "channelBC.H"
#include "channelBCF_F.H"
#include "probF_F.H"

class ConstValueFunction: public BCValueFunction
{
  public:
    Real m_value;
    int m_nComp;

    ConstValueFunction()
      :
      m_nComp(0)
    {
    }

    ConstValueFunction(Real a_value,int a_nComp)
      :
      m_value(a_value),
      m_nComp(a_nComp)
    {
    }

    virtual void operator()(Real*           a_pos,
                            int*            a_dir,
                            Side::LoHiSide* a_side,
                            Real*           a_value)
    {
      if (m_nComp > 0)
      {
        for (int comp = 0; comp < m_nComp; comp++)
        {
          a_value[comp] = m_value;
        }
      }
      else
      {
        MayDay::Error("undefined ConstValueFunction object");
      }
    }
};

channelBC::channelBC()
{
  // initialize to bogus values
  for (int dir=0; dir<SpaceDim; dir++)
    {
      m_loBC[dir] = bogusBC;
      m_hiBC[dir] = bogusBC;
    }

  m_maxInflowVel = 1.0e8;

  // now initialize BCs
  setBCs();
}

channelBC::~channelBC()
{
}

PhysBCUtil*
channelBC::newPhysBCUtil() const
{
  channelBC* newBCPtr = new channelBC;
  newBCPtr->setMaxInflowVel(m_maxInflowVel);
  return static_cast<PhysBCUtil*>(newBCPtr);
}

void
channelBC::setMaxInflowVel(Real a_maxInflowVel)
{
  m_maxInflowVel = a_maxInflowVel;
}

Real
channelBC::maxInflowVel() const
{
  return m_maxInflowVel;
}

void
channelBC::setBCs()
{
  ParmParse ppBC("physBC");

  vector<int> tempBC(SpaceDim);

  ppBC.getarr("lo", tempBC, 0, SpaceDim);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      m_loBC[dir] = tempBC[dir];
    }

  ppBC.getarr("hi", tempBC, 0, SpaceDim);

  for (int dir=0; dir<SpaceDim; dir++)
    {
      m_hiBC[dir] = tempBC[dir];
    }

  ppBC.get("maxInflowVel", m_maxInflowVel);
}

void
channelBC::computeBoundaryDt(Real& a_dt, Real a_cfl, Real a_Dx) const
{
  Real maxInflowDt = a_cfl*a_Dx/m_maxInflowVel;

  if (maxInflowDt < a_dt)
  {
    a_dt = maxInflowDt;
  }
}
