#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CHANNELBC_H_
#define _CHANNELBC_H_

#include "PhysBCUtil.H"

#include "UsingNamespace.H"

/// specific BC's for flow in a channel
class channelBC : public PhysBCUtil
{

public:
  /// constructor
  channelBC();

  /// destructor
  virtual ~channelBC();

  /// "virtual" constructor
  virtual PhysBCUtil* newPhysBCUtil() const;

  ///
  void setMaxInflowVel(Real a_maxInflowVel);

  ///
  Real maxInflowVel() const;

  /// interact with ParmParse to set physical BC types
  virtual void setBCs();

  /// compute effect of boundary conditions on dt
  virtual void computeBoundaryDt(Real& a_dt, Real a_cfl, Real a_dx) const;

protected:

  ///
  //  virtual BCHolder basicScalarFuncBC(bool a_isHomogeneous,
  //                                     int a_scalType) const;

  //  virtual BCHolder basicCCVelBC(bool a_isHomogeneous,
  //                                bool a_isViscous,
  //                                int a_velComponent) const;

  /// max inflow velocity
  Real m_maxInflowVel;

private:

  /// copy constructor
  channelBC(const channelBC& rhs)
  {
  }

};

#endif
