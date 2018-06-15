#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ExtrapAdvectBC.H"
#include "EBISLayout.H"
#include "EBLoHiCenter.H"
#include "VoFIterator.H"

#include "NamespaceHeader.H"

/****************************/
ExtrapAdvectBC::
ExtrapAdvectBC()
  :EBPhysIBC()
{
  m_isDefined = false;
}
/****************************/
void
ExtrapAdvectBC::
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

  Box FBox = a_primGdnv[a_dir].getRegion();
  Box cellBox = a_faceBox;
  int numFlux = a_primGdnv[a_dir].nComp();

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
              for (int ivar = 0; ivar < numFlux; ivar++)
                {
                  a_primGdnv[a_dir](face, ivar) = a_primExtrap(vof,ivar);
                }
            }
        }
    }
}
/****************************/
void
ExtrapAdvectBC::define(const ProblemDomain&  a_domain,
                       const RealVect& a_dx)
{
  m_domain = a_domain;
  m_dx = a_dx;
  m_isDefined = true;
}
/****************************/
void
ExtrapAdvectBC::
initialize(LevelData<EBCellFAB>& a_conState,
           const EBISLayout& a_ebisl) const
{
  MayDay::Error("should not be called");
}
/****************************/
ExtrapAdvectBC::
~ExtrapAdvectBC()
{
}
/******************/
ExtrapAdvectBCFactory::
ExtrapAdvectBCFactory()
  :EBPhysIBCFactory()
{
}
/******************/
ExtrapAdvectBCFactory::
~ExtrapAdvectBCFactory()
{
}
/******************/
EBPhysIBC*
ExtrapAdvectBCFactory::
create() const
{
  ExtrapAdvectBC* retval = new ExtrapAdvectBC();

  return static_cast<EBPhysIBC*>(retval);
}
/******************/

#include "NamespaceFooter.H"
