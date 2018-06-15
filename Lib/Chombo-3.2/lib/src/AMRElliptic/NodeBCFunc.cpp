#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NodeBCFunc.H"
#include "RealVect.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"

void NodeNeumBC(NodeFArrayBox&  a_state,
                const Box&      a_valid,
                Real            a_dx,
                bool            a_homogeneous,
                BCValueHolder   a_value)
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          NodeNeumBC(a_state, a_valid, a_dx, a_homogeneous, a_value, idir,sit());
        }
    }
}

static void getDomainFacePositionNode(RealVect&             a_retval,
                                      const IntVect&        a_nodeIV,
                                      const Real&           a_dx,
                                      const int&            a_dir,
                                      const Side::LoHiSide& a_side)
{
  Real* dataPtr = a_retval.dataPtr();

  D_TERM( dataPtr[0] = a_dx*(a_nodeIV[0]);,\
          dataPtr[1] = a_dx*(a_nodeIV[1]);,\
          dataPtr[2] = a_dx*(a_nodeIV[2]);)

}
void NodeNeumBC(NodeFArrayBox& a_state,
                const Box&     a_valid,
                Real           a_dx,
                bool           a_homogeneous,
                BCValueHolder  a_value,
                int            a_dir,
                Side::LoHiSide a_side)
{

  int isign = sign(a_side);
  Box toRegionCell = adjCellBox(a_valid, a_dir, a_side, 1);
  Box toRegionNode = surroundingNodes(toRegionCell);
  toRegionNode.growDir(a_dir, a_side, -1);

  Real* value  = new Real[a_state.nComp()];
  for (BoxIterator bit(toRegionNode); bit.ok(); ++bit)
    {
      const IntVect& ivTo = bit();
      IntVect ivClose = ivTo -   isign*BASISV(a_dir);
      if (!a_homogeneous)
        {
          RealVect facePos;
          getDomainFacePositionNode(facePos, ivTo, a_dx, a_dir, a_side);
          a_value(facePos.dataPtr(), &a_dir, &a_side, value);
        }
      for (int icomp = 0; icomp < a_state.nComp(); icomp++)
        {
          Real inhomogValue = 0;
          if (!a_homogeneous)
            {
              inhomogValue = value[icomp];
            }
          // on high side
          //dphi/dx = (valTo - valClose)/dx = value
          // valTo = value*dx + valClose
          //on low side
          //dphi/dx = -(valTo - valClose)/dx = value
          // valTo = -value*dx + valClose
          Real valClose = a_state(ivClose, icomp);
          Real valTo = isign*(inhomogValue*a_dx) + valClose;
          a_state(ivTo, icomp) = valTo;
        }
    }
  delete[] value;
}

void NodeDiriBC(NodeFArrayBox&  a_state,
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

          NodeDiriBC(a_state, a_valid, a_dx, a_homogeneous, a_value, idir, sit(), a_order);
        }
    }
}
void NodeDiriBC(NodeFArrayBox& a_state,
                const Box&     a_valid,
                Real           a_dx,
                bool           a_homogeneous,
                BCValueHolder  a_value,
                int            a_dir,
                Side::LoHiSide a_side,
                int            a_order)
{
  Box toRegionCell = adjCellBox(a_valid, a_dir, a_side, 1);
  Box toRegionNode = surroundingNodes(toRegionCell);
  toRegionNode.growDir(a_dir, a_side, -1);

  Real* value  = new Real[a_state.nComp()];
  for (BoxIterator bit(toRegionNode); bit.ok(); ++bit)
    {
      const IntVect& ivTo = bit();
      if (!a_homogeneous)
        {
          RealVect facePos;
          getDomainFacePositionNode(facePos, ivTo, a_dx, a_dir, a_side);
          a_value(facePos.dataPtr(), &a_dir, &a_side, value);
        }
      for (int icomp = 0; icomp < a_state.nComp(); icomp++)
        {
          const IntVect& ivTo = bit();
          Real inhomogValue = 0;
          if (!a_homogeneous)
            {
              inhomogValue = value[icomp];
            }
          // The one upside of node-based solvers.
          // Boundary conditions are what they are.
          a_state(ivTo, icomp) = inhomogValue;
        }
    }
  delete[] value;
}

#include "NamespaceFooter.H"
