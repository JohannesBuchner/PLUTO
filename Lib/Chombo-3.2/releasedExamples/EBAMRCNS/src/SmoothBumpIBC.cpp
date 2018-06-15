#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <cmath>

#include "SmoothBumpIBC.H"
#include "SmoothBumpF_F.H"
#include "EBISLayout.H"
#include "EBLoHiCenter.H"
#include "VoFIterator.H"
#include "EBLGIntegrator.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

/****************************/
/****************************/
SmoothBumpIBC::
SmoothBumpIBC(const Real&     a_gamma,
               const RealVect& a_center,
               const Real&     a_mach,
               const Real&     a_rnot)
  :EBPhysIBC()
{
  Real smallness;
  ParmParse pp;
  pp.get("smallness", smallness);
  /**/
  FORT_SETSMOOTHBUMP(CHF_CONST_REAL(a_gamma),
                     CHF_CONST_REALVECT(a_center),
                     CHF_CONST_REAL(a_mach),
                     CHF_CONST_REAL(a_rnot),
                     CHF_CONST_REAL(smallness));
  /**/

  m_gamma    = a_gamma;
  m_mach       = a_mach;
  m_rnot       = a_rnot;
  m_center   = a_center;
  m_isFortranCommonSet = true;
  m_isDefined = false;
}
/****************************/
/****************************/
void
SmoothBumpIBC::
setBndrySlopes(EBCellFAB&       a_deltaPrim,
               const EBCellFAB& a_primState,
               const EBISBox&   a_ebisBox,
               const Box&       a_box,
               const int&       a_dir)
{
  //already set to one-sided
}
/****************************/
/****************************/
void
SmoothBumpIBC::
fluxBC(EBFluxFAB&            a_flux,
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
  //none of this matters because we should not have any domain BCs.
 CH_assert(m_isDefined);
 CH_assert(m_isFortranCommonSet);
 CH_assert(!m_domain.isPeriodic(a_dir));

  Box FBox = a_flux[a_dir].getSingleValuedFAB().box();
  Box cellBox = FBox;
  int numFlux = a_flux[a_dir].nComp();

  // Determine which side and thus shifting directions
  int isign = sign(a_side);
  cellBox.shiftHalf(a_dir,isign);

  // Is there a domain boundary next to this grid
  if (!m_domain.contains(cellBox))
    {
      cellBox &= m_domain;
      // Find the strip of cells next to the domain boundary
      Box boundaryBox = bdryBox(cellBox, a_dir, a_side, 1);

      // Shift things to all line up correctly
      boundaryBox.shiftHalf(a_dir,-isign);
      BaseFab<Real>& regFlux = a_flux[a_dir].getSingleValuedFAB();
      const BaseFab<Real>& regPrimExtrap = a_primExtrap.getSingleValuedFAB();
      regFlux.shiftHalf(a_dir,-isign);

      // Set the boundary fluxes
      /**/
      FORT_SOLIDBC(CHF_FRA(regFlux),
                   CHF_CONST_FRA(regPrimExtrap),
                   CHF_CONST_INT(isign),
                   CHF_CONST_REAL(m_dx),
                   CHF_CONST_INT(a_dir),
                   CHF_BOX(boundaryBox));

      // Shift returned fluxes to be face centered
      regFlux.shiftHalf(a_dir,isign);
      int inormMomVar = CMOMX + a_dir;
      int inormVelVar = QVELX + a_dir;
      //now for the multivalued cells.  Since it is pointwise,
      //the regular calc is correct for all single-valued cells.
      IntVectSet ivs = a_ebisBox.getMultiCells(boundaryBox);
      for(VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          //set all fluxes except normal momentum  to zero.
          //set normal momemtum flux to sign(normal)*pressure
          Vector<FaceIndex> bndryFaces = a_ebisBox.getFaces(vof, a_dir, a_side);
          for(int iface= 0; iface < bndryFaces.size(); iface++)
            {
              const FaceIndex& face = bndryFaces[iface];
              //set all fluxes to zero then fix momentum flux
              for(int ivar = 0; ivar < numFlux; ivar++)
                {
                  a_flux[a_dir](face, ivar) = 0;
                }
              Real press = a_primExtrap(vof, QPRES);
              Real dense = a_primExtrap(vof, QRHO);
              Real unorm = a_primExtrap(vof, inormVelVar);
              Real speed = sqrt(m_gamma*press/dense);

              //blessed cartesian coords.
              a_flux[a_dir](face, inormMomVar) =
                press + isign*dense*unorm*speed;
            }
        }
    }
}
/****************************/
/****************************/
void
SmoothBumpIBC::define(const ProblemDomain&  a_domain,
                      const RealVect&       a_dx)
{
  m_domain = a_domain;
  m_dx = a_dx[0];
  m_isDefined = true;
}
/****************************/
/****************************/
void
SmoothBumpIBC::
initialize(LevelData<EBCellFAB>& a_conState,
           const EBISLayout& a_ebisl) const
{
 CH_assert(m_isDefined);
 CH_assert(m_isFortranCommonSet);

  // Iterator of all grids in this level
  for(DataIterator dit = a_conState.dataIterator(); dit.ok(); ++dit)
  {
    const EBISBox& ebisBox = a_ebisl[dit()];
    // Storage for current grid
    EBCellFAB& conFAB = a_conState[dit()];

    BaseFab<Real>& regConFAB = conFAB.getSingleValuedFAB();
    // Box of current grid
    Box uBox = regConFAB.box();
    uBox &= m_domain;
    // Set up initial condition in this grid
    /**/
    FORT_SMOOTHBUMPINIT(CHF_CONST_FRA(regConFAB),
                        CHF_CONST_REAL(m_dx),
                        CHF_BOX(uBox));
    /**/

    //now for the multivalued cells.  Since it is pointwise,
    //the regular calc is correct for all single-valued cells.
    IntVectSet ivs = ebisBox.getMultiCells(uBox);
    for(VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        const IntVect& iv = vof.gridIndex();
        RealVect momentum;
        Real energy, density;
        /**/
        FORT_POINTSMOOTHBUMPINIT(CHF_REAL(density),
                                 CHF_REALVECT(momentum),
                                 CHF_REAL(energy),
                                 CHF_CONST_INTVECT(iv),
                                 CHF_CONST_REAL(m_dx));
        /**/
        conFAB(vof, CRHO) = density;
        conFAB(vof, CENG) = energy;
        for(int idir = 0; idir < SpaceDim; idir++)
          {
            conFAB(vof, CMOMX+idir) = momentum[idir];
          }
      }//end loop over multivalued cells
  } //end loop over boxes
}
/****************************/
/****************************/
SmoothBumpIBC::
~SmoothBumpIBC()
{
}
/****************************/
/****************************/

#include "NamespaceFooter.H"
