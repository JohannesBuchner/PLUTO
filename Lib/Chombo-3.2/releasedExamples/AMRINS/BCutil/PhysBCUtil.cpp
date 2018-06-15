#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/*****************/
/*****************/

#include "BoxIterator.H"
#include "PhysBCUtil.H"
#include "RefCountedPtr.H"
#include "viscousBCF_F.H"

#include "AdvectIBC.H"
#include "AdvectScalarIBC.H"
#include "VelIBC.H"

#include "MayDay.H"

// ---------------------------------------------------------------
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
    : m_value(a_value),
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

// ---------------------------------------------------------------
extern "C"
{
  static void ExtraBC(FArrayBox&     a_state,
                      const Box&     a_valid,
                      int            a_dir,
                      Side::LoHiSide a_side,
                      int            a_order)
  {
    int isign = sign(a_side);
    Box to_region = adjCellBox(a_valid, a_dir, a_side, 1);

    for (BoxIterator bit(to_region); bit.ok(); ++bit)
    {
      const IntVect& ivTo = bit();
      IntVect ivFrom1 = ivTo -   isign*BASISV(a_dir);
      IntVect ivFrom2 = ivTo - 2*isign*BASISV(a_dir);
      IntVect ivFrom3 = ivTo - 3*isign*BASISV(a_dir);

      for (int icomp = 0; icomp < a_state.nComp(); icomp++)
      {
        Real ghostVal;
        if (a_order == 1)
        {
          Real val1 = a_state(ivFrom1, icomp);
          Real val2 = a_state(ivFrom2, icomp);

          ghostVal = 2.0*val1 - val2;
        }
        else if (a_order == 2)
        {
          Real val1 = a_state(ivFrom1, icomp);
          Real val2 = a_state(ivFrom2, icomp);
          Real val3 = a_state(ivFrom3, icomp);

          ghostVal = 3.0*(val1 - val2) + val3;
        }
        else
        {
          MayDay::Error("ExtrapBC - order argument not 1 or 2");
        }

        a_state(ivTo, icomp) = ghostVal;
      }
    }
  }

  static void streamBC(FArrayBox&           a_state,
                       const Box&           a_valid,
                       const ProblemDomain& a_domain,
                       Real                 a_dx,
                       bool                 a_homogeneous)
  {
    const Box& domainBox = a_domain.domainBox();
    PhysBCUtil bcInfo; // sets BCs from ParmParse table
    RefCountedPtr<ConstValueFunction>
      zeroFunc(new ConstValueFunction(0.0, a_state.nComp()));

    for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (!a_domain.isPeriodic(idir))
          {
            SideIterator sit;
            for (sit.reset(); sit.ok(); ++sit)
              {
                Side::LoHiSide side = sit();
                if (a_valid.sideEnd(side)[idir] ==
                    domainBox.sideEnd(side)[idir])
                  {
                    int bctype = bcInfo.getBC(idir, side);
                    switch (bctype)
                      {
                      case PhysBCUtil::SolidWall :
                      case PhysBCUtil::noShear :
                        {
                          // Doesn't work with higher order
                          int order = 1;
                          DiriBC(a_state, a_valid, a_dx,
                                 a_homogeneous,
                                 BCValueHolder(zeroFunc),
                                 idir, side, order);
                          break;
                        }
                      case PhysBCUtil::Inflow :
                        {
                          MayDay::Error("streamBC - Don't know how to set BC inflow");
                          break;
                        }
                      case PhysBCUtil::Outflow :
                        {
                          int order = 1;
                          ExtraBC(a_state, a_valid,
                                  idir, side, order);
                          break;
                        }
                      default :
                        {
                          MayDay::Error("streamBC - unknown BC type");
                        }
                      } // end switch
                  } // if ends match
              } // end loop over sides
          } // if not periodic in this direction
      } // end loop over directions
  }
}


// ---------------------------------------------------------------
class BasicScalarBCFunction: public BCFunction
{
public:

  int m_scalarType;

  BasicScalarBCFunction()
    :
    m_scalarType(-1)
  {
  }

  BasicScalarBCFunction(int a_scalarType)
    :
    m_scalarType(a_scalarType)
  {
  }

  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_scalarType >= 0)
      {
        const Box& domainBox = a_domain.domainBox();
        PhysBCUtil bcInfo; // sets BCs from ParmParse table

        for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (!a_domain.isPeriodic(idir))
              {
                SideIterator sit;
                for (sit.reset(); sit.ok(); ++sit)
                  {
                    Side::LoHiSide side = sit();
                    if (a_valid.sideEnd(side)[idir] ==
                        domainBox.sideEnd(side)[idir])
                      {
                        int bctype = bcInfo.getBC(idir, side);
                        switch(bctype)
                          {
                          case PhysBCUtil::SolidWall :
                          case PhysBCUtil::Outflow :
                          case PhysBCUtil::noShear :
                            {
                              int order = 1;
                              ExtraBC(a_state, a_valid, idir, side, order);
                              break;
                            }
                          case PhysBCUtil::Inflow :
                            {
                              MayDay::Error("BasicScalarBCFunction - BC inflow should not be used here");
                              break;
                            }
                          default :
                            {
                              MayDay::Error("BasicScalarBCFunction - unknown BC type");
                            }
                          } // end switch
                      } // end condition of matching coordinates
                  } // end loop over both sides
              } // end if not periodic in this direction
          } // end loop over all directions
      } // if (m_scalarType >= 0)
    else // m_scalarType < 0
      {
        MayDay::Error("undefined BasicScalarBCFunction object");
      } // if (m_scalarType >= 0)
  }
};


// ---------------------------------------------------------------
class FreestreamCorrBCFunction: public BCFunction
{
public:

  bool m_isDefined;

  FreestreamCorrBCFunction()
    :
    m_isDefined(false)
  {
  }

  FreestreamCorrBCFunction(bool a_isDefined)
    :
    m_isDefined(a_isDefined)
  {
  }

  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_isDefined)
      {
        const Box& domainBox = a_domain.domainBox();
        PhysBCUtil bcInfo; // sets BCs from ParmParse table
        RefCountedPtr<ConstValueFunction>
          zeroFunc(new ConstValueFunction(0.0, a_state.nComp()));

        for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (!a_domain.isPeriodic(idir))
              {
                SideIterator sit;
                for (sit.reset(); sit.ok(); ++sit)
                  {
                    Side::LoHiSide side = sit();
                    if (a_valid.sideEnd(side)[idir] ==
                        domainBox.sideEnd(side)[idir])
                      {
                        int bctype = bcInfo.getBC(idir, side);
                        switch(bctype)
                          {
                          case PhysBCUtil::Inflow :
                          case PhysBCUtil::Outflow :
                          case PhysBCUtil::SolidWall :
                          case PhysBCUtil::noShear :
                            {
                              NeumBC(a_state, a_valid, a_dx,
                                     a_homogeneous,
                                     BCValueHolder(zeroFunc),
                                     idir, side);
                              break;
                            }
                          default:
                            {
                              MayDay::Error("FreestreamCorrBCFunction - unknown BC type");
                            }
                          } // end switch
                      } // end condition of matching coordinates
                  } // end loop over both sides
              } // end if not periodic in this direction
          } // end loop over all directions
      }
    else
      {
        MayDay::Error("undefined FreestreamCorrBCFunction object");
      }
  }
};

// ---------------------------------------------------------------
class BasicLambdaBCFunction: public BCFunction
{
public:

  bool m_isDefined;

  BasicLambdaBCFunction()
    :
    m_isDefined(false)
  {
  }

  BasicLambdaBCFunction(bool a_isDefined)
    :
    m_isDefined(a_isDefined)
  {
  }

  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_isDefined)
      {
        const Box& domainBox = a_domain.domainBox();
        PhysBCUtil bcInfo; // sets BCs from ParmParse table
        RefCountedPtr<ConstValueFunction>
          oneFunc(new ConstValueFunction(1.0, a_state.nComp()));

        for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (!a_domain.isPeriodic(idir))
              {
                SideIterator sit;
                for (sit.reset(); sit.ok(); ++sit)
                  {
                    Side::LoHiSide side = sit();
                    if (a_valid.sideEnd(side)[idir] ==
                        domainBox.sideEnd(side)[idir])
                      {
                        int bctype = bcInfo.getBC(idir, side);
                        switch(bctype)
                          {
                          case PhysBCUtil::SolidWall :
                          case PhysBCUtil::Outflow :
                          case PhysBCUtil::noShear :
                            {
                              int order = 1;
                              ExtraBC(a_state, a_valid,
                                      idir, side, order);
                              break;
                            }
                          case PhysBCUtil::Inflow :
                            {
                              // Doesn't work with higher order
                              int order = 1;
                              DiriBC(a_state, a_valid, a_dx,
                                     a_homogeneous,
                                     BCValueHolder(oneFunc),
                                     idir, side, order);
                              break;
                            }
                          default:
                            {
                              MayDay::Error("BasicLambdaBCFunction - unknown BC type");
                            }
                          } // end switch
                      } // end condition of matching coordinates
                  } // end loop over both sides
              } // end if not periodic in this direction
          } // end loop over all directions
      }
    else
      {
        MayDay::Error("undefined BasicLambdaBCFunction object");
      }
  }
};

// ---------------------------------------------------------------
class BasicPressureBCFunction: public BCFunction
{
public:

  bool m_isDefined;

  bool m_isHomogeneous;

  BasicPressureBCFunction()
    :
    m_isDefined(false)
  {
  }

  BasicPressureBCFunction(bool a_isDefined,
                          bool a_isHomogeneous)
    :
    m_isDefined(a_isDefined),
    m_isHomogeneous(a_isHomogeneous)
  {
  }

  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_isDefined)
      {
        const Box& domainBox = a_domain.domainBox();
        PhysBCUtil bcInfo; // sets BCs from ParmParse table
        RefCountedPtr<ConstValueFunction>
          zeroFunc(new ConstValueFunction(0.0, a_state.nComp()));

        for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (!a_domain.isPeriodic(idir))
              {
                SideIterator sit;
                for (sit.reset(); sit.ok(); ++sit)
                  {
                    Side::LoHiSide side = sit();
                    if (a_valid.sideEnd(side)[idir] ==
                        domainBox.sideEnd(side)[idir])
                      {
                        int bctype = bcInfo.getBC(idir, side);
                        switch(bctype)
                          {
                          case PhysBCUtil::SolidWall :
                          case PhysBCUtil::Inflow :
                          case PhysBCUtil::noShear :
                            {
                              NeumBC(a_state, a_valid, a_dx,
                                     a_homogeneous,
                                     BCValueHolder(zeroFunc),
                                     idir, side);
                              break;
                            }
                          case PhysBCUtil::Outflow :
                            {
                              // Doesn't work with higher order
                              int order = 1;
                              DiriBC(a_state, a_valid, a_dx,
                                     a_homogeneous,
                                     BCValueHolder(zeroFunc),
                                     idir, side, order);
                              break;
                            }
                          case PhysBCUtil::Symmetry :
                            {
                              MayDay::Error("BasicPressureBCFunction - ReflectBC not implemented");
                              break;
                            }
                          default :
                            {
                              MayDay::Error("BasicPressureBCFunction - unknown BC type");
                            }
                          } // end switch
                      } // end condition of matching coordinates
                  } // end loop over both sides
              } // end if not periodic in this direction
          } // end loop over all directions
      }
    else
      {
        MayDay::Error("undefined BasicPressureBCFunction object");
      }
  }
};

// ---------------------------------------------------------------
class ExtrapBCFunction: public BCFunction
{
public:

  bool m_isDefined;

  int m_order;

  ExtrapBCFunction()
    :
    m_isDefined(false)
  {
  }

  ExtrapBCFunction(bool a_isDefined,
                   int a_order)
    :
    m_isDefined(a_isDefined),
    m_order(a_order)
  {
  }

  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_isDefined)
      {
        const Box& domainBox = a_domain.domainBox();
        PhysBCUtil bcInfo; // sets BCs from ParmParse table
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (!a_domain.isPeriodic(idir))
              {
                SideIterator sit;
                for (sit.reset(); sit.ok(); ++sit)
                  {
                    Side::LoHiSide side = sit();
                    if (a_valid.sideEnd(side)[idir] ==
                        domainBox.sideEnd(side)[idir])
                      {
                        ExtraBC(a_state, a_valid,
                                idir, side, m_order);
                      } // end condition of matching coordinates
                  } // end loop over both sides
              } // end if not periodic in this direction
          } // end loop over all directions
      }
    else
      {
        MayDay::Error("undefined ExtrapBCFunction object");
      }
  }
};

// ---------------------------------------------------------------
class BasicGradPressureBCFunction: public BCFunction
{
public:

  bool m_isDefined;

  BasicGradPressureBCFunction()
    :
    m_isDefined(false)
  {
  }

  BasicGradPressureBCFunction(bool a_isDefined)
    :
    m_isDefined(a_isDefined)
  {
  }

  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_isDefined)
      {
        const Box& domainBox = a_domain.domainBox();
        PhysBCUtil bcInfo; // sets BCs from ParmParse table

        for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (!a_domain.isPeriodic(idir))
              {
                SideIterator sit;
                for (sit.reset(); sit.ok(); ++sit)
                  {
                    Side::LoHiSide side = sit();
                    if (a_valid.sideEnd(side)[idir] ==
                        domainBox.sideEnd(side)[idir])
                      {
                        int bctype = bcInfo.getBC(idir, side);
                        switch(bctype)
                          {
                          case PhysBCUtil::SolidWall :
                          case PhysBCUtil::Inflow :
                          case PhysBCUtil::Outflow :
                          case PhysBCUtil::noShear :
                            {
                              int order = 2; // equivalent of HOExtrapBC
                              ExtraBC(a_state, a_valid,
                                      idir, side, order);
                              break;
                            }
                          case PhysBCUtil::Symmetry :
                            {
                              MayDay::Error("BasicGradPressureBCFunction - ReflectOddBC not implemented");
                              break;
                            }
                          default:
                            {
                              MayDay::Error("BasicGradPressureBCFunction - unknown BC type");
                            }
                          } // end switch
                      } // end condition of matching coordinates
                  } // end loop over both sides
              } // end if not periodic in this direction
          } // end loop over all directions
      }
    else
      {
        MayDay::Error("undefined BasicGradPressureBCFunction object");
      }
  }
};

// ---------------------------------------------------------------
class BasicCCVelBCFunction: public BCFunction
{
public:

  Real m_bcVal;

  bool m_isHomogeneous;

  bool m_isViscous;

  int m_comp;

  Interval m_interval;

  BasicCCVelBCFunction()
    :
    m_comp(-1)
  {
  }

  BasicCCVelBCFunction(Real a_bcVal,
                       bool a_isHomogeneous,
                       bool a_isViscous,
                       int  a_comp,
                       const Interval& a_interval)
    : m_bcVal(a_bcVal),
      m_isHomogeneous(a_isHomogeneous),
      m_isViscous(a_isViscous),
      m_comp(a_comp),
      m_interval(a_interval)
  {
  }

  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_comp >= 0)
      {
        const Box& domainBox = a_domain.domainBox();
        PhysBCUtil bcInfo; // sets BCs from ParmParse table

        FArrayBox aliasStateFab(m_interval, a_state);

        RefCountedPtr<ConstValueFunction>
          zeroFunc(new ConstValueFunction(0.0, aliasStateFab.nComp()));
        RefCountedPtr<ConstValueFunction>
          bcValueFunc(new ConstValueFunction(m_bcVal, aliasStateFab.nComp()));

        // loop over directions
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (!a_domain.isPeriodic(idir))
              {
                SideIterator sit;
                for (sit.reset(); sit.ok(); ++sit)
                  {
                    Side::LoHiSide side = sit();
                    if (a_valid.sideEnd(side)[idir] ==
                        domainBox.sideEnd(side)[idir])
                      {
                        int bctype = bcInfo.getBC(idir, side);
                        switch (bctype)
                          {
                          case PhysBCUtil::SolidWall :
                            {
                              // normal velocity BCs:
                              // always no-flow
                              if (idir == m_comp)
                                {
                                  int order = 1;
                                  DiriBC(aliasStateFab, a_valid, a_dx,
                                         a_homogeneous,
                                         BCValueHolder(zeroFunc),
                                         idir, side, order);
                                }
                              else
                                {
                                  // tangential BCs:
                                  // no-slip if viscous, extrap if inviscid
                                  if (m_isViscous)
                                    {
                                      int order = 1;
                                      DiriBC(aliasStateFab, a_valid, a_dx,
                                             a_homogeneous,
                                             BCValueHolder(zeroFunc),
                                             idir, side, order);
                                    }
                                  else // inviscid
                                    {
                                      int order = 1;
                                      ExtraBC(aliasStateFab, a_valid,
                                              idir, side, order);
                                    }
                                } // end if tangential
                              break;
                            }
                          case PhysBCUtil::Inflow :
                            {
                              // this will most likely get overwritten in a derived class
                              if (!m_isHomogeneous && idir == m_comp)
                                {
                                  int order = 1;
                                  DiriBC(aliasStateFab, a_valid, a_dx,
                                         a_homogeneous,
                                         BCValueHolder(bcValueFunc),
                                         idir, side, order);
                                }
                              else
                                {
                                  int order = 1;
                                  DiriBC(aliasStateFab, a_valid, a_dx,
                                         a_homogeneous,
                                         BCValueHolder(zeroFunc),
                                         idir, side, order);
                                }
                              break;
                            }
                          case PhysBCUtil::Outflow :
                            {
                              NeumBC(aliasStateFab, a_valid, a_dx,
                                     a_homogeneous,
                                     BCValueHolder(zeroFunc),
                                     idir, side);
                              break;
                            }
                          case PhysBCUtil::noShear :
                            {
                              // normal velocity BC's still no-flow
                              if (idir == m_comp)
                                {
                                  int order = 1;
                                  DiriBC(aliasStateFab, a_valid, a_dx,
                                         a_homogeneous,
                                         BCValueHolder(zeroFunc),
                                         idir, side, order);
                                }
                              else
                                {
                                  NeumBC(aliasStateFab, a_valid, a_dx,
                                         a_homogeneous,
                                         BCValueHolder(zeroFunc),
                                         idir, side);
                                } // end if tangential
                              break;
                            }
                          default :
                            {
                              MayDay::Error("BasicCCVelBCFunction - unknown BC type");
                            }
                          } // end switch
                      } // if ends match
                  } // end iteration over sides
              } // if not periodic in this direction
          } // end iteration over directions
      } // if m_comp >= 0
    else // m_comp < 0
      {
        MayDay::Error("undefined BasicCCVelBCFunction object");
      }
  }
};

// ---------------------------------------------------------------
class BasicECVelBCFunction: public BCFunction
{
public:

  Real m_bcVal;

  bool m_isHomogeneous;

  bool m_isViscous;

  int m_comp;

  Interval m_interval;

  BasicECVelBCFunction()
    :
    m_comp(-1)
  {
  }

  BasicECVelBCFunction(Real a_bcVal,
                       bool a_isHomogeneous,
                       bool a_isViscous,
                       int  a_comp,
                       const Interval& a_interval)
    : m_bcVal(a_bcVal),
      m_isHomogeneous(a_isHomogeneous),
      m_isViscous(a_isViscous),
      m_comp(a_comp),
      m_interval(a_interval)
  {
  }

  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    // a_state is FACE-centered, now in m_comp direction;
    // a_valid is CELL-centered
    if (m_comp >= 0)
      {
        const Box& domainBox = a_domain.domainBox();
        PhysBCUtil bcInfo; // sets BCs from ParmParse table

        RefCountedPtr<ConstValueFunction>
          zeroFunc(new ConstValueFunction(0.0, a_state.nComp()));
        RefCountedPtr<ConstValueFunction>
          bcValueFunc(new ConstValueFunction(m_bcVal, a_state.nComp()));

        // loop over directions
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (!a_domain.isPeriodic(idir))
              {
                SideIterator sit;
                for (sit.reset(); sit.ok(); ++sit)
                  {
                    Side::LoHiSide side = sit();
                    if (a_valid.sideEnd(side)[idir] ==
                        domainBox.sideEnd(side)[idir])
                      {
                        int bctype = bcInfo.getBC(idir, side);
                        switch (bctype)
                          {
                          case PhysBCUtil::SolidWall :
                            {
                              // normal velocity BCs:
                              // always no-flow
                              if (idir == m_comp)
                                {
                                  DiriEdgeBC(a_state, a_valid, a_dx,
                                             a_homogeneous,
                                             BCValueHolder(zeroFunc),
                                             idir, side);
                                }
                              else
                                {
                                  // tangential BCs:
                                  // no-slip if viscous, extrap if inviscid
                                  if (m_isViscous)
                                    {
                                      int order = 1;
                                      // need to fake this a bit
                                      // want to use Cell-centered
                                      // BC function to set tangential BC's
                                      // on face-centered data. Do this by
                                      // shifting valid-region to
                                      // face-centering.  DiriBC function
                                      // only really cares that valid box
                                      // and state box have the same centering
                                      // (DFM -- 9/22/08)
                                      Box validFace(a_valid);
                                      validFace.surroundingNodes(m_comp);
                                      DiriBC(a_state, validFace, a_dx,
                                             a_homogeneous,
                                             BCValueHolder(zeroFunc),
                                             idir, side, order);
                                    }
                                  else // inviscid
                                    {
                                      int order = 1;
                                      ExtraBC(a_state, a_valid,
                                              idir, side, order);
                                    }
                                } // end if tangential
                              break;
                            }
                          case PhysBCUtil::Inflow :
                            {
                              // this will most likely get overwritten in a derived class
                              if (!m_isHomogeneous && idir == m_comp)
                                {
                                  DiriEdgeBC(a_state, a_valid, a_dx,
                                             a_homogeneous,
                                             BCValueHolder(bcValueFunc),
                                             idir, side);
                                }
                              else
                                {
                                  int order = 1;
                                  DiriBC(a_state, a_valid, a_dx,
                                         a_homogeneous,
                                         BCValueHolder(zeroFunc),
                                         idir, side, order);
                                }
                              break;
                            }
                          case PhysBCUtil::Outflow :
                            {
                              // this is set to whatever it is set to by MAC
                              // NoOp
                              break;
                            }
                          case PhysBCUtil::noShear :
                            {
                              // normal velocity BC's still no-flow
                              if (idir == m_comp)
                                {
                                  DiriEdgeBC(a_state, a_valid, a_dx,
                                             a_homogeneous,
                                             BCValueHolder(zeroFunc),
                                             idir, side);
                                }
                              else
                                {
                                  // tangential BC is do-nothing BC???
                                  // NoOp
                                } // end if tangential
                              break;
                            }
                          default :
                            {
                              MayDay::Error("BasicECVelBCFunction - unknown BC type");
                            }
                          } // end switch
                      } // if ends match
                  } // end iteration over sides
              } // if not periodic in this direction
          } // end iteration over directions
      } // if m_comp >= 0
    else // m_comp < 0
      {
        MayDay::Error("undefined BasicECVelBCFunction object");
      }
  }

  virtual void DiriEdgeBC(FArrayBox&           a_state,
                          const Box&           a_valid,
                          Real                 a_dx,
                          bool                 a_homogeneous,
                          BCValueHolder        a_value,
                          int                  a_idir,
                          Side::LoHiSide       a_side)
  {
    // a_state is FACE-centered on face idir;
    // a_valid is CELL-centered.
    Box face(a_valid);
    face.surroundingNodes(a_idir);
    int coord = face.sideEnd(a_side)[a_idir];
    face.setRange(a_idir, coord);

    if (a_homogeneous)
      {
        a_state.setVal(0.);
      }
    else
      {
        RealVect facePos;
        int ncomp = a_state.nComp();
        Real* value = new Real[ncomp];
        RealVect junk;
        a_value(junk.dataPtr(), &a_idir, &a_side, value);
        for (int comp = 0; comp < ncomp; comp++)
          {
            a_state.setVal(value[comp], face, comp);
          }
      }
  }
};

// ---------------------------------------------------------------
class ViscousBCFunction: public BCFunction
{
public:

  bool m_isDefined;

  ViscousBCFunction()
    :
    m_isDefined(false)
  {
  }

  ViscousBCFunction(bool a_isDefined)
    :
    m_isDefined(a_isDefined)
  {
  }

  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
  {
    if (m_isDefined)
      {
        int nComp = a_state.nComp();
        const Box& grownBox = a_state.box();
        PhysBCUtil bcInfo; // sets BCs from ParmParse table

        for (int idir = 0; idir < SpaceDim; idir++)
          {
            // Huh, these will be undefined if conditions below don't hold.
            Box loBox, hiBox;

            if (grownBox.smallEnd(idir) < a_valid.smallEnd(idir))
              {
                loBox = grownBox;
                loBox.setBig(idir, a_valid.smallEnd(idir)-1);
              }
            if (grownBox.bigEnd(idir) > a_valid.bigEnd(idir))
              {
                hiBox = grownBox;
                hiBox.setSmall(idir, a_valid.bigEnd(idir)+1);
              }

            FORT_VISCOUSBC(CHF_FRA(a_state),
                           CHF_BOX(a_valid),
                           CHF_BOX(loBox),
                           CHF_BOX(hiBox),
                           CHF_INT(idir),
                           CHF_INT(nComp));
          } // end loop over all directions
      }
    else
      {
        MayDay::Error("undefined ViscousBCFunction object");
      }
  }
};

// ---------------------------------------------------------------
PhysBCUtil::PhysBCUtil()
{
  // initialize to bogus values
  for (int idir=0; idir<SpaceDim; idir++)
    {
      m_loBC[idir] = bogusBC;
      m_hiBC[idir] = bogusBC;
    }

  // now initialize BC's
  setBCs();
}

// ---------------------------------------------------------------
PhysBCUtil::~PhysBCUtil()
{
}

// ---------------------------------------------------------------
int
PhysBCUtil::loBC(int a_dir) const
{
  return m_loBC[a_dir];
}

// ---------------------------------------------------------------
int
PhysBCUtil::hiBC(int a_dir) const
{
  return m_hiBC[a_dir];
}

// ---------------------------------------------------------------
int
PhysBCUtil::getBC(int a_dir, const Side::LoHiSide& a_side) const
{
  int bcVal;
  if (a_side == Side::Lo)
    {
      bcVal = loBC(a_dir);
    }
  else if (a_side == Side::Hi)
    {
      bcVal = hiBC(a_dir);
    }
  else
    {
      MayDay::Error("PhysBCUtil::getBC -- invalid side");
    }
  return bcVal;
}

// ---------------------------------------------------------------
PhysBCUtil&
PhysBCUtil::operator= (const PhysBCUtil& rhs)
{
  m_loBC = rhs.m_loBC;
  m_hiBC = rhs.m_hiBC;
  return *this;
}

// ---------------------------------------------------------------
PhysBCUtil::PhysBCUtil(const PhysBCUtil& rhs)
{
  m_loBC = rhs.m_loBC;
  m_hiBC = rhs.m_hiBC;
}

// ---------------------------------------------------------------
void
PhysBCUtil::setBCs()
{
  ParmParse ppBC("physBC");

  vector<int> tempBC(SpaceDim);

  ppBC.getarr("lo", tempBC, 0, SpaceDim);
  for (int idir=0; idir<SpaceDim; idir++)
    {
      m_loBC[idir] = tempBC[idir];
    }

  ppBC.getarr("hi", tempBC, 0, SpaceDim);
  for (int idir=0; idir<SpaceDim; idir++)
    {
      m_hiBC[idir] = tempBC[idir];
    }
}

// ---------------------------------------------------------------
PhysBCUtil*
PhysBCUtil::newPhysBCUtil() const
{
  PhysBCUtil* newBCPtr = new PhysBCUtil;
  return newBCPtr;
}

// ----------------------------------------------------------
// basic utility functions
void
PhysBCUtil::Dx(const Real a_dx)
{
  m_dx = a_dx;
}

// ---------------------------------------------------------------
Real
PhysBCUtil::Dx() const
{
  return m_dx;
}

// ---------------------------------------------------------------
void
PhysBCUtil::Time(const Real a_time)
{
  m_time = a_time;
}

// ---------------------------------------------------------------
Real
PhysBCUtil::Time() const
{
  return m_time;
}

// -----------------------------------------------------------
// Now begin basic BC functions.

// ---------------------------------------------------------------
Tuple<BCHolder, SpaceDim>
PhysBCUtil::viscousVelFuncBC() const
{
  Tuple<BCHolder, SpaceDim> bcVec;
  bool isHomogeneous = false;
  bool isViscous = true;
  for (int idir=0; idir<SpaceDim; idir++)
    {
      Interval intvl(idir, idir);
      bcVec[idir] = basicCCVelFuncBC(isHomogeneous, isViscous, idir,
                                     intvl);
    }
  return bcVec;
}

// ---------------------------------------------------------------
Tuple<BCHolder, SpaceDim>
PhysBCUtil::uDelUFuncBC(bool a_isViscous) const
{
  Tuple<BCHolder, SpaceDim> bcVec;
  bool isHomogeneous = true;
  for (int idir=0; idir<SpaceDim; idir++)
    {
      Interval intvl(idir, idir);
      bcVec[idir] = basicCCVelFuncBC(isHomogeneous, a_isViscous, idir,
                                     intvl);
    }
  return bcVec;
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::tracingSolveFuncBC(int a_idir) const
{
  bool isHomogeneous = false;
  bool isViscous = false;
  Interval intvl(0, 0);
  return basicCCVelFuncBC(isHomogeneous, isViscous, a_idir,
                          intvl);
}

// ---------------------------------------------------------------
Tuple<BCHolder, SpaceDim>
PhysBCUtil::tracingVelFuncBC() const
{
  Tuple<BCHolder, SpaceDim> bcVec;
  bool isHomogeneous = false;
  bool isViscous = false;
  for (int idir=0; idir<SpaceDim; idir++)
    {
      Interval intvl(idir, idir);
      bcVec[idir] = basicCCVelFuncBC(isHomogeneous, isViscous, idir,
                                     intvl);
    }
  return bcVec;
}

// ---------------------------------------------------------------
Tuple<BCHolder, SpaceDim>
PhysBCUtil::advectionVelFuncBC(bool a_isViscous) const
{
  Tuple<BCHolder, SpaceDim> bcVec;
  bool isHomogeneous = false;
  for (int idir=0; idir<SpaceDim; idir++)
    {
      // data will have one component in each space dimension
      Interval intvl(0, 0);
      bcVec[idir] = basicECVelFuncBC(isHomogeneous, a_isViscous, idir,
                                     intvl);
    }
  return bcVec;
}

// ---------------------------------------------------------------
Tuple<BCHolder, SpaceDim>
PhysBCUtil::uStarFuncBC(bool a_isViscous) const
{
  Tuple<BCHolder, SpaceDim> bcVec;
  bool isHomogeneous = false;
  for (int idir=0; idir<SpaceDim; idir++)
    {
      Interval intvl(idir, idir);
      bcVec[idir] = basicCCVelFuncBC(isHomogeneous, a_isViscous, idir,
                                     intvl);
    }
  return bcVec;
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::viscousRefluxBC(int a_dir) const
{
  bool isViscous = true;
  bool isHomogeneous = true;
  Interval intvl(0, 0); // 13 nov 2007:  not sure about this
  return basicCCVelFuncBC(isHomogeneous, isViscous, a_dir,
                          intvl);
}

// ---------------------------------------------------------------
/// this is a BC object used in the PatchGodunov stuff
PhysIBC*
PhysBCUtil::advectionVelIBC() const
{
  VelIBC* newIBCPtr = new VelIBC;
  // default to solid wall BC's everywhere
  // (i.e. set normal velocities to 0)
  // if something different (like inflow) is desired, do this in the
  // derived BC server class.
  for (int idir=0; idir<SpaceDim; idir++)
    {
      SideIterator sit;
      for (sit.reset(); sit.ok(); ++sit)
        {
          Side::LoHiSide side = sit();
          newIBCPtr->setNormalWallVel(0.0, idir, side);
        }
    }

  return newIBCPtr;
}

// -------------------------------------------------------
/// this is a BC object used in the PatchGodunov stuff
PhysIBC*
PhysBCUtil::scalarTraceExtrapIBC(int a_scalarType) const
{
  AdvectScalarIBC* newIBCPtr = new AdvectScalarIBC;
  return newIBCPtr;
}

// -------------------------------------------------------
BCHolder PhysBCUtil::scalarRefluxSolveBC(int a_scalarType) const
{
  bool isHomogeneous = true;
  return basicScalarFuncBC(isHomogeneous, a_scalarType);
}

// -------------------------------------------------------
BCHolder PhysBCUtil::scalarTraceFuncBC(int a_scalarType) const
{
  bool isHomogeneous = false;
  return basicScalarFuncBC(isHomogeneous, a_scalarType);
}

// -------------------------------------------------------
BCHolder PhysBCUtil::lambdaFuncBC() const
{
  return basicLambdaFuncBC();
}

// -------------------------------------------------------
/// this is a BC object used in the PatchGodunov stuff
PhysIBC*
PhysBCUtil::scalarTraceIBC(int a_scalarType) const
{
  AdvectIBC* newIBCPtr = new AdvectIBC;
  return newIBCPtr;
}

// ---------------------------------------------------------------
/// this is a BC object used in the PatchGodunov stuff
PhysIBC*
PhysBCUtil::lambdaTraceIBC() const
{
  AdvectIBC* newIBCPtr = new AdvectIBC;
  // for lambda, all boundary values are 1
  for (int idir=0; idir<SpaceDim; idir++)
    {
      SideIterator sit;
      for (sit.reset(); sit.ok(); ++sit)
        {
          Side::LoHiSide side = sit();
          newIBCPtr->setBoundaryValue(1.0, idir, side);
        }
    }

  return newIBCPtr;
}

// ---------------------------------------------------------------
BCFunc PhysBCUtil::streamBC() const
{
  return ::streamBC;
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::basicCCVelFuncBC(bool a_isHomogeneous,
                                      bool a_isViscous,
                                      int  a_comp,
                                      const Interval& a_interval) const
{
  Real bcVal = 1.0;
  RefCountedPtr<BasicCCVelBCFunction>
    basicCCVelBCFunction(new BasicCCVelBCFunction(bcVal,
                                                  a_isHomogeneous,
                                                  a_isViscous,
                                                  a_comp,
                                                  a_interval));

  return BCHolder(basicCCVelBCFunction);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::basicECVelFuncBC(bool a_isHomogeneous,
                                      bool a_isViscous,
                                      int  a_comp,
                                      const Interval& a_interval) const
{
  Real bcVal = 1.0;
  RefCountedPtr<BasicECVelBCFunction>
    basicECVelBCFunction(new BasicECVelBCFunction(bcVal,
                                                  a_isHomogeneous,
                                                  a_isViscous,
                                                  a_comp,
                                                  a_interval));

  return BCHolder(basicECVelBCFunction);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::basicScalarFuncBC(bool a_isHomogeneous,
                                       int  a_scalarType) const
{
  RefCountedPtr<BasicScalarBCFunction>
    basicScalarBCFunction(new BasicScalarBCFunction(a_scalarType));

  return BCHolder(basicScalarBCFunction);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::basicLambdaFuncBC() const
{
  RefCountedPtr<BasicLambdaBCFunction>
    basicLambdaBCFunction(new BasicLambdaBCFunction(true));

  return BCHolder(basicLambdaBCFunction);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::extrapFuncBC(int a_order) const
{
  RefCountedPtr<ExtrapBCFunction>
    extrapBCFunction(new ExtrapBCFunction(true, a_order));

  return BCHolder(extrapBCFunction);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::FreestreamCorrFuncBC() const
{
  RefCountedPtr<FreestreamCorrBCFunction>
    freestreamCorrBCFunction(new FreestreamCorrBCFunction(true));

  return BCHolder(freestreamCorrBCFunction);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::BasicPressureFuncBC(bool a_isHomogeneous) const
{
  RefCountedPtr<BasicPressureBCFunction>
    basicPressureBCFunction(new BasicPressureBCFunction(true,
                                                        a_isHomogeneous));

  return BCHolder(basicPressureBCFunction);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::BasicGradPressureFuncBC() const
{
  RefCountedPtr<BasicGradPressureBCFunction>
    basicGradPressureBCFunction(new BasicGradPressureBCFunction(true));

  return BCHolder(basicGradPressureBCFunction);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::viscousSolveFuncBC(int a_idir) const
{
  bool isHomogeneous = false;
  bool isViscous = true;
  Interval intvl(0, 0);
  return basicCCVelFuncBC(isHomogeneous, isViscous, a_idir,
                          intvl);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::viscousFuncBC() const
{
  RefCountedPtr<ViscousBCFunction>
    viscousBCFunction(new ViscousBCFunction(true));

  return BCHolder(viscousBCFunction);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::LevelPressureFuncBC() const
{
  return BasicPressureFuncBC(false);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::SyncProjFuncBC() const
{
  return BasicPressureFuncBC(true);
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::gradMacPressureFuncBC() const
{
  return BasicGradPressureFuncBC();
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::gradPiFuncBC() const
{
  return BasicGradPressureFuncBC();
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::gradESyncFuncBC() const
{
  return BasicGradPressureFuncBC();
}

// ---------------------------------------------------------------
BCHolder PhysBCUtil::gradELambdaFuncBC() const
{
  return BasicGradPressureFuncBC();
}

// ---------------------------------------------------------------
void
PhysBCUtil::computeBoundaryDt(Real& a_dt, Real a_cfl, Real a_dx) const
{
  // in the default case, nothing is done here.
}
