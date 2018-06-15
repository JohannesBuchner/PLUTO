#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NoFlowAdvectBC.H"
#include "EBISLayout.H"
#include "EBLoHiCenter.H"
#include "VoFIterator.H"

#include "NamespaceHeader.H"

/****************************/
NoFlowAdvectBC::
NoFlowAdvectBC(int a_velComp)
  :EBPhysIBC()
{
  m_velComp = a_velComp;
  m_isDefined = false;
}
/****************************/
void
NoFlowAdvectBC::
fluxBC(EBFluxFAB&            a_primGdnv,
       const EBCellFAB&      a_primCenter,
       const EBCellFAB&      a_primExtrap,
       const Side::LoHiSide& a_side,
       const Real&           a_time,
       const EBISBox&        a_ebisBox,
       const DataIndex&      a_dit,
       const Box&            a_box,
       const Box&            a_faceBox,
       const int&            a_dir)
{
  CH_assert(m_isDefined);

  Box FBox = a_faceBox;
  Box cellBox = FBox;
  CH_assert(a_primGdnv[a_dir].nComp()==1);

  // Determine which side and thus shifting directions
  int isign = sign(a_side);
  cellBox.shiftHalf(a_dir,isign);

  // Is there a domain boundary next to this grid
  if (!m_domain.contains(cellBox))
    {
      cellBox &= m_domain;
      // Find the strip of cells next to the domain boundary
      Box bndryBox = adjCellBox(cellBox, a_dir, a_side, 1);

      // Shift things to all line up correctly
      bndryBox.shift(a_dir,-isign);

      IntVectSet ivs(bndryBox);
      for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();

          Vector<FaceIndex> bndryFaces = a_ebisBox.getFaces(vof, a_dir, a_side);
          for (int iface= 0; iface < bndryFaces.size(); iface++)
            {
              const FaceIndex& face = bndryFaces[iface];
              //set all fluxes to zero then fix momentum flux
              //solid wall
              if (a_dir == m_velComp)
                {
                  a_primGdnv[a_dir](face, 0) = 0.0;
                }
              else
                {
                  a_primGdnv[a_dir](face, 0) = a_primExtrap(vof, 0);
                }
            }
        }
    }
}
/****************************/
void
NoFlowAdvectBC::define(const ProblemDomain&  a_domain,
                       const RealVect&       a_dx)
{
  m_domain = a_domain;
  m_dx = a_dx;
  m_isDefined = true;
}
/****************************/
void
NoFlowAdvectBC::
initialize(LevelData<EBCellFAB>& a_conState,
           const EBISLayout& a_ebisl) const
{
  MayDay::Error("should not be called");
}
/****************************/
NoFlowAdvectBC::
~NoFlowAdvectBC()
{
}
/******************/
NoFlowAdvectBCFactory::
NoFlowAdvectBCFactory(int a_velComp)
  :EBPhysIBCFactory()
{
  m_velComp = a_velComp;
}
/******************/
NoFlowAdvectBCFactory::
~NoFlowAdvectBCFactory()
{
}
/******************/
EBPhysIBC*
NoFlowAdvectBCFactory::
create() const
{
  NoFlowAdvectBC* retval = new NoFlowAdvectBC(m_velComp);

  return static_cast<EBPhysIBC*>(retval);
}
/******************/

#include "NamespaceFooter.H"
