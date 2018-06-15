#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BCFunc.H"
#include "RealVect.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"
//#include "BCFuncF_F.H"

void doNothingBC(FArrayBox&           a_state,
                 const Box&           a_valid,
                 const ProblemDomain& a_domain,
                 Real                 a_dx,
                 bool                 a_homogeneous)
{
}

static inline Real linearInterp(const Real& a_inhomogVal,
                                const Real& a_nearVal)
{

  // ignore value at x=1
  Real retval = 2*a_inhomogVal-a_nearVal;
  return retval;
}

static inline Real quadraticInterp(const Real& a_inhomogVal,
                                   const Real& a_nearVal,
                                   const Real& a_farVal)
{
  Real retval = (8.0/3.0)*a_inhomogVal + (1.0/3.0)*(a_farVal) -2*(a_nearVal);
  return retval;
}

ConstBCFunction::ConstBCFunction(const IntVect&  a_loSideType,
                                 const RealVect& a_loSideValue,
                                 const IntVect&  a_hiSideType,
                                 const RealVect& a_hiSideValue)
{
  m_loSideType  = a_loSideType;
  m_loSideValue = a_loSideValue;

  m_hiSideType  = a_hiSideType;
  m_hiSideValue = a_hiSideValue;
}

ConstBCFunction::~ConstBCFunction()
{
}

void ConstBCFunction::operator()(FArrayBox&           a_state,
                                 const Box&           a_valid,
                                 const ProblemDomain& a_domain,
                                 Real                 a_dx,
                                 bool                 a_homogeneous)
{
  const Box& domainBox = a_domain.domainBox();

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    if (!a_domain.isPeriodic(idir))
    {
      for (SideIterator sit; sit.ok(); ++sit)
      {
        Side::LoHiSide side = sit();

        if (a_valid.sideEnd(side)[idir] == domainBox.sideEnd(side)[idir])
        {
          int  bcType;
          Real bcValue;

          if (side == Side::Lo)
          {
            bcType  = m_loSideType [idir];
            bcValue = m_loSideValue[idir];
          }
          else
          {
            bcType  = m_hiSideType [idir];
            bcValue = m_hiSideValue[idir];
          }

          if (bcType == 0)
          {
            // Neumann BC
            int isign = sign(side);

            Box toRegion = adjCellBox(a_valid, idir, side, 1);
            toRegion &= a_state.box();

            Box fromRegion = toRegion;
            fromRegion.shift(idir, -isign);

            a_state.copy(a_state, fromRegion, 0, toRegion, 0, a_state.nComp());

            if (!a_homogeneous)
            {
              for (BoxIterator bit(toRegion); bit.ok(); ++bit)
              {
                const IntVect& ivTo = bit();
                // IntVect ivClose = ivTo - isign*BASISV(idir);

                for (int icomp = 0; icomp < a_state.nComp(); icomp++)
                {
                  a_state(ivTo, icomp) += Real(isign)*a_dx*bcValue;
                }
              }
            }
          }
          else if (bcType == 1)
          {
            // Dirichlet BC
            int isign = sign(side);

            Box toRegion = adjCellBox(a_valid, idir, side, 1);
            toRegion &= a_state.box();

            for (BoxIterator bit(toRegion); bit.ok(); ++bit)
            {
              const IntVect& ivTo = bit();

              IntVect ivClose = ivTo -   isign*BASISV(idir);
              // IntVect ivFar   = ivTo - 2*isign*BASISV(idir);

              Real inhomogVal = 0.0;

              if (!a_homogeneous)
              {
                inhomogVal = bcValue;
              }

              for (int icomp = 0; icomp < a_state.nComp(); icomp++)
              {
                Real nearVal = a_state(ivClose, icomp);
                // Real farVal  = a_state(ivFar,   icomp);

                Real ghostVal = linearInterp(inhomogVal, nearVal);

                a_state(ivTo, icomp) = ghostVal;
              }
            }
          }
          else
          {
            MayDay::Abort("ConstBCFunction::operator() - unknown BC type");
          }
        } // if ends match
      } // end loop over sides
    } // if not periodic in this direction
  } // end loop over directions
}

RefCountedPtr<BCFunction> ConstDiriNeumBC(const IntVect&  a_loSideType,
                                          const RealVect& a_loSideValue,
                                          const IntVect&  a_hiSideType,
                                          const RealVect& a_hiSideValue)
{
  RefCountedPtr<BCFunction> retval(new ConstBCFunction(a_loSideType,
                                                         a_loSideValue,
                                                         a_hiSideType,
                                                         a_hiSideValue));

  return retval;
}

void NeumBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            BCValueHolder   a_value)
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          NeumBC(a_state, a_valid, a_dx, a_homogeneous, a_value, idir,sit());
        }
    }
}

static void getDomainFacePosition(RealVect&             a_retval,
                                  const IntVect&        a_validIV,
                                  const Real&           a_dx,
                                  const int&            a_dir,
                                  const Side::LoHiSide& a_side)
{
  Real* dataPtr = a_retval.dataPtr();

  D_TERM( dataPtr[0] = a_dx*(a_validIV[0] + 0.5);,\
          dataPtr[1] = a_dx*(a_validIV[1] + 0.5);,\
          dataPtr[2] = a_dx*(a_validIV[2] + 0.5);)

  int isign = sign(a_side);
  dataPtr[a_dir] += 0.5*Real(isign)*a_dx;
}

void NeumBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            const BCValueHolder&   a_valueA,
            int             a_dir,
            Side::LoHiSide  a_side)
{
  Interval stateInterval = a_state.interval();
  NeumBC(a_state, a_valid, a_dx, a_homogeneous,
         a_valueA, a_dir,a_side, stateInterval);
}

void NeumBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            const BCValueHolder&   a_valueA,
            int             a_dir,
            Side::LoHiSide  a_side,
            Interval&        a_interval)
{
  BCValueHolder& a_value = (BCValueHolder&)a_valueA;
  int isign = sign(a_side);
  RealVect facePos;

  Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
  toRegion &= a_state.box();

  Box fromRegion = toRegion;
  fromRegion.shift(a_dir, -isign);

  a_state.copy(a_state, fromRegion, 0, toRegion, 0, a_state.nComp());

  if (!a_homogeneous)
    {
      Real* value = new Real[a_state.nComp()];

      for (BoxIterator bit(toRegion); bit.ok(); ++bit)
        {
          const IntVect& ivTo = bit();
          IntVect ivClose = ivTo - isign*BASISV(a_dir);

          getDomainFacePosition(facePos, ivClose, a_dx, a_dir, a_side);

          a_value(facePos.dataPtr(), &a_dir, &a_side, value);
          for (int icomp = a_interval.begin(); icomp <= a_interval.end(); icomp++)
            {
              a_state(ivTo, icomp) += Real(isign)*a_dx*value[icomp];
            }
        }

      delete[] value;
    }
}

void DiriBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            BCValueHolder   a_value,
            int             a_order)
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {

          DiriBC(a_state, a_valid, a_dx, a_homogeneous, a_value, idir, sit(), a_order);
        }
    }
}

void DiriBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            BCValueHolder   a_value,
            int             a_dir,
            Side::LoHiSide  a_side,
            int             a_order)
{
  Interval stateInterval = a_state.interval();
  DiriBC(a_state,a_valid, a_dx, a_homogeneous,
         a_value,a_dir, a_side, stateInterval, a_order);
}

void DiriBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            BCValueHolder   a_value,
            int             a_dir,
            Side::LoHiSide  a_side,
            Interval&        a_interval,
            int             a_order)
{
  int isign = sign(a_side);
  RealVect facePos;

  Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
  toRegion &= a_state.box();

  Real* value = new Real[a_state.nComp()];

  for (BoxIterator bit(toRegion); bit.ok(); ++bit)
    {
      const IntVect& ivTo = bit();

      IntVect ivClose = ivTo -   isign*BASISV(a_dir);
      IntVect ivFar   = ivTo - 2*isign*BASISV(a_dir);

      if (!a_homogeneous)
        {
          getDomainFacePosition(facePos, ivClose, a_dx, a_dir, a_side);
          a_value(facePos.dataPtr(), &a_dir, &a_side, value);
        }

      for (int icomp = a_interval.begin(); icomp <= a_interval.end(); icomp++)
        {
          Real nearVal = a_state(ivClose, icomp);
          Real farVal  = a_state(ivFar,   icomp);

          Real inhomogVal = 0.0;

          if (!a_homogeneous)
            {
              inhomogVal = value[icomp];
            }

          Real ghostVal=0;

          if (a_order == 1)
            {
              ghostVal = linearInterp(inhomogVal, nearVal);
            }
          else if (a_order == 2)
            {
              ghostVal = quadraticInterp(inhomogVal, nearVal, farVal);
            }
          else
            {
              MayDay::Error("bogus order argument");
            }

          a_state(ivTo, icomp) = ghostVal;
        }
    }

  delete[] value;
}

void NoSlipVectorBC(FArrayBox&     a_state,
                    const Box&     a_valid,
                    Real           a_dx,
                    int            a_dir,
                    Side::LoHiSide a_side,
                    int            a_order)
{
  int isign = sign(a_side);
  RealVect facePos;

  Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);

  for (BoxIterator bit(toRegion); bit.ok(); ++bit)
    {
      const IntVect& ivTo = bit();
      IntVect ivClose = ivTo -   isign*BASISV(a_dir);
      IntVect ivFar   = ivTo - 2*isign*BASISV(a_dir);

      Real inhomogVal = 0.0;

      for (int icomp = 0; icomp < a_state.nComp(); icomp++)
        {
          Real nearVal = a_state(ivClose, icomp);
          Real farVal  = a_state(ivFar,   icomp);

          Real ghostVal=0;

          if (a_order == 1)
            {
              ghostVal = linearInterp(inhomogVal, nearVal);
            }
          else if (a_order == 2)
            {
              ghostVal = quadraticInterp(inhomogVal, nearVal, farVal);
            }

          a_state(ivTo, icomp) = ghostVal;
        }
    }
}

void ReflectiveVectorBC(FArrayBox&     a_state,
                        const Box&     a_valid,
                        Real           a_dx,
                        int            a_dir,
                        Side::LoHiSide a_side,
                        int            a_order)
{
  int isign = sign(a_side);
  RealVect facePos;

  Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);

  for (BoxIterator bit(toRegion); bit.ok(); ++bit)
    {
      const IntVect& ivTo = bit();
      IntVect ivClose = ivTo -   isign*BASISV(a_dir);
      IntVect ivFar   = ivTo - 2*isign*BASISV(a_dir);

      Real inhomogVal = 0.0;

      for (int icomp = 0; icomp < a_state.nComp(); icomp++)
        {
          Real ghostVal=0;
          Real nearVal = a_state(ivClose, icomp);

          // zero dirichlet for normal dir
          if (icomp == a_dir)
            {
              Real farVal = a_state(ivFar,   icomp);
              if (a_order == 1)
                {
                  ghostVal = linearInterp(inhomogVal, nearVal);
                }
              else if (a_order == 2)
                {
                  ghostVal = quadraticInterp(inhomogVal, nearVal, farVal);
                }
            }
          else
            {
              // zero neumann for tang dir
              ghostVal = nearVal;
            }
          a_state(ivTo, icomp) = ghostVal;
        }
    }
}

///
/**
   Extrapolation boundary conditions for a side, specified component interval
   For use in AMRPoissonOp.
 */
void ExtrapolateBC(FArrayBox&      a_state,
                   const Box&      a_valid,
                   Real            a_dx,
                   int             a_dir,
                   Side::LoHiSide  a_side,
                   Interval&       a_interval,
                   int             a_order)
{
  // for now, only 1st-order (linear) extrapolation implemented
  CH_assert(a_order == 1);

  int isign = sign(a_side);

  Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
  toRegion &= a_state.box();

  for (BoxIterator bit(toRegion); bit.ok(); ++bit)
    {
      const IntVect& ivTo = bit();

      IntVect ivClose = ivTo -   isign*BASISV(a_dir);
      IntVect ivFar   = ivTo - 2*isign*BASISV(a_dir);

      for (int icomp = a_interval.begin(); icomp <= a_interval.end(); icomp++)
        {
          Real nearVal = a_state(ivClose, icomp);
          Real farVal  = a_state(ivFar,   icomp);

          Real ghostVal=0;

          if (a_order == 1)
            {
              // for the linear extrap case, can use the "linearInterp" function
              ghostVal = linearInterp(nearVal, farVal);
            }
          else
            {
              MayDay::Error("bogus order argument");
            }

          a_state(ivTo, icomp) = ghostVal;
        }
    }

}

///
/**
   Extrapolation boundary conditions for a side, all components.
   For use in AMRPoissonOp.
 */
void ExtrapolateBC(FArrayBox&      a_state,
                   const Box&      a_valid,
                   Real            a_dx,
                   int             a_dir,
                   Side::LoHiSide  a_side,
                   int             a_order)
{
  // for now, only 1st-order (linear) extrapolation implemented
  CH_assert(a_order == 1);

  Interval stateInterval = a_state.interval();
  ExtrapolateBC(a_state,a_valid, a_dx, a_dir, a_side,
                stateInterval, a_order);

}

///
/**
   Extrapolation boundary conditions for one side.
   For use in AMRPoissonOp.
 */
void ExtrapolateBC(FArrayBox&      a_state,
                   const Box&      a_valid,
                   Real            a_dx,
                   int             a_order)
{
  // for now, only 1st-order (linear) extrapolation implemented
  CH_assert(a_order == 1);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          ExtrapolateBC(a_state, a_valid, a_dx, idir, sit(), a_order);
        }
    }
}

#include "NamespaceFooter.H"
