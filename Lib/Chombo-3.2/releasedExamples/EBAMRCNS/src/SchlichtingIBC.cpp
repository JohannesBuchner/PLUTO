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

#include "SchlichtingIBC.H"
#include "EBISLayout.H"
#include "EBLoHiCenter.H"
#include "ParmParse.H"
#include "VoFIterator.H"
#include "EBLGIntegrator.H"
#include "EBLevelDataOps.H"
#include "EBPatchPolytropicF_F.H"

#include "NamespaceHeader.H"

/****************************/
/****************************/
SchlichtingIBC::
SchlichtingIBC(const SchlichtingParams&  a_params)
{
  m_params = a_params;
  FORT_SETGAMMAANDSMALL(CHF_CONST_REAL(a_params.m_gamma));

  m_isDefined = false;
}
/****************************/
/****************************/
void
SchlichtingIBC::
setBndrySlopes(EBCellFAB&       a_deltaPrim,
               const EBCellFAB& a_primState,
               const EBISBox&   a_ebisBox,
               const Box&       a_box,
               const int&       a_dir)
{
}
/****************************/
/****************************/
void
SchlichtingIBC::
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

  Box FBox = a_flux[a_dir].getSingleValuedFAB().box();
  Box cellBox = FBox;

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
      IntVectSet ivs(boundaryBox);
      VoFIterator vofit(ivs, a_ebisBox.getEBGraph());
      // Set the boundary fluxes 
      for(vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Vector<Real> qgdnv(QNUM);
          Vector<Real> fluxv(FNUM);
          if(a_side == Side::Lo)  //inflow side
            {
              qgdnv[QRHO] = m_params.m_initDense;
              for(int idir = 0; idir < SpaceDim; idir++)
                {
                  qgdnv[QVELX+idir] = m_params.m_velMag*m_params.m_axis[idir];
                }
              qgdnv[QPRES] = m_params.m_initPress;
            }
          else //outflow side --extrapolation
            {
              for(int ivar = 0; ivar < QNUM; ivar++)
                {
                  qgdnv[ivar] = a_primExtrap(vof, ivar);
                }
            }
          Vector<FaceIndex> bndryFaces = a_ebisBox.getFaces(vof, a_dir, a_side);
          for(int iface= 0; iface < bndryFaces.size(); iface++)
            {
              /**/
              FORT_POINTGETFLUX(CHF_VR(fluxv),
                                CHF_VR(qgdnv),
                                CHF_CONST_INT(a_dir));
              /**/
              const FaceIndex& face = bndryFaces[iface];
              for(int ivar = 0; ivar < FNUM; ivar++)
                {
                  a_flux[a_dir](face, ivar) = fluxv[ivar];
                }
            }
        }
    }
}
/****************************/
/****************************/
void
SchlichtingIBC::define(const ProblemDomain&  a_domain,
                       const RealVect& a_dx)
{
  m_domain = a_domain;
  m_dx = a_dx[0];
  m_isDefined = true;
}
/****************************/
/****************************/
void
SchlichtingIBC::
initialize(LevelData<EBCellFAB>& a_conState,
           const EBISLayout& a_ebisl) const
{
  CH_assert(m_isDefined);
  Real initialR = m_params.m_initDense;
  Real initialV = m_params.m_velMag;


  RealVect initialMom = m_params.m_initVel;
  initialMom *= initialR;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      EBLevelDataOps::setVal(a_conState, initialMom[idir], CMOMX+idir);
    }

  Real initKE = 0.5*initialR*initialV*initialV;
  Real initIE = initialR*m_params.m_specificHeat*m_params.m_initTemp;

  Real initialE = initKE + initIE;
  EBLevelDataOps::setVal(a_conState, initialR, CRHO);
  EBLevelDataOps::setVal(a_conState, initialE, CENG);
}
/****************************/
/****************************/
SchlichtingIBC::
~SchlichtingIBC()
{
}
/****************************/
/****************************/

#include "NamespaceFooter.H"
