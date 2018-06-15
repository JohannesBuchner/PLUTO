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

#include "EBNoFlowIBC.H"
#include "EBISLayout.H"
#include "EBLoHiCenter.H"
#include "VoFIterator.H"
#include "EBLGIntegrator.H"
#include "ParmParse.H"

// These are used here because we have cut-n-pasted from the 
// EBExplosionIBC class. :-(
#include "EBExplosionF_F.H"
#include "EBSolidF_F.H"

#include "NamespaceHeader.H"

//------------------------------------------------------------------------------
EBNoFlowIBC::
EBNoFlowIBC(Real a_gamma)
  :EBPhysIBC(),
  m_dx(),
  m_domain(),
  m_isDefined(false)
{
  // Because the model relies on Fortran code to set the pressure, we must 
  // initialize the polytropic index gamma within the Fortran common block.
  // Currently, this is done via the IBC classes. Since I don't want to write 
  // a ChF just for this boundary condition, I choose to call the initializer 
  // for the EBExplosionIBC below. We should probably revisit how this 
  // initializaton is done in a boundary-condition-independent way. -JNJ
  Real zero = 0.0;
  RealVect vecZero;
  FORT_SETEXPLOSION(CHF_CONST_REAL(a_gamma),
                    CHF_CONST_REALVECT(vecZero),
                    CHF_CONST_REAL(zero),
                    CHF_CONST_REAL(zero),
                    CHF_CONST_REAL(zero),
                    CHF_CONST_REAL(zero),
                    CHF_CONST_REAL(zero),
                    CHF_CONST_REAL(zero));
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void
EBNoFlowIBC::
setBndrySlopes(EBCellFAB&       a_deltaPrim,
               const EBCellFAB& a_primState,
               const EBISBox&   a_ebisBox,
               const Box&       a_box,
               const int&       a_dir)
{
  // This should never be used with periodic domains.
  CH_assert(!m_domain.isPeriodic(a_dir));

  Box loBox,hiBox,centerBox,domain;
  int hasLo,hasHi;
  Box slopeBox = a_deltaPrim.getRegion()&m_domain;
  int numSlop = a_deltaPrim.nComp();
  // Generate the domain boundary boxes, loBox and hiBox, if there are
  // domain boundaries there
  eblohicenter(loBox,hasLo,hiBox,hasHi,centerBox,domain,
               slopeBox,m_domain,a_dir);

  // Set the boundary slopes if necessary
  if ((hasLo != 0) || (hasHi != 0))
  {
    BaseFab<Real>& regDeltaPrim = a_deltaPrim.getSingleValuedFAB();
    const BaseFab<Real>& regPrimState = a_primState.getSingleValuedFAB();
    /**/
    FORT_SLOPEBCS(CHF_FRA(regDeltaPrim),
        CHF_CONST_FRA(regPrimState),
        CHF_CONST_REAL(m_dx[0]),
        CHF_CONST_INT(a_dir),
        CHF_BOX(loBox),
        CHF_CONST_INT(hasLo),
        CHF_BOX(hiBox),
        CHF_CONST_INT(hasHi));
    /**/
  }
  for(SideIterator sit; sit.ok(); ++sit)
  {
    bool doThisSide;
    Box thisBox;
    if(sit() == Side::Lo)
    {
      doThisSide = (hasLo != 0);
      thisBox = loBox;
    }
    else
    {
      doThisSide = (hasHi != 0);
      thisBox = hiBox;
    }
    if(doThisSide)
    {
      // the cells for which the regular calculation
      // are incorrect are the grown set of the multi-valued
      // cells intersected with the appropriate bndry box
      // and added to the irregular cells on the boundary.
      Box boxGrown = thisBox;
      boxGrown.grow(a_dir, 1);
      boxGrown &= m_domain;
      IntVectSet ivs = a_ebisBox.getMultiCells(boxGrown);
      ivs &= loBox;
      IntVectSet ivsIrreg = a_ebisBox.getIrregIVS(thisBox);
      ivs |= ivsIrreg;
      for(VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        //all slopes at boundary get set to zero
        //except normal velocity.  just following the
        //fortran here
        int inormVelVar = a_dir + QVELX;
        for(int ivar = 0; ivar < numSlop; ivar++)
        {
          if(ivar != inormVelVar)
          {
            a_deltaPrim(vof, ivar) = 0.0;
          }
        }

        //for normal velocity
        //do strangely limited slope
        //just lifted from the fortran.
        Side::LoHiSide otherSide = flip(sit());
        Vector<FaceIndex> faces =
          a_ebisBox.getFaces(vof, a_dir, otherSide);

        Real thisVel = a_primState(vof, inormVelVar);
        Real otherVel = 0.;
        for(int iface = 0; iface < faces.size(); iface++)
        {
          VolIndex otherVoF = faces[iface].getVoF(otherSide);
          otherVel += a_primState(otherVoF,inormVelVar);
        }
        if(faces.size() > 0)
        {
          otherVel /= Real(faces.size());
          Real slope;
          if(sit() == Side::Lo)
          {
            slope = otherVel - thisVel;
          }
          else
          {
            slope = thisVel  - otherVel;
          }
          //trick to make sure state will not change sign
          if(slope*thisVel < 0)
          {
            a_deltaPrim(vof, inormVelVar) = 0.0;
          }
          else
          { //slope and vel same sign
            Real rsign = 1.0;
            if(thisVel < 0.0)
              rsign = -1.0;
            //told you this was odd.
            Real dvmin = Min(Abs(slope), Abs((Real)2.*thisVel));
            a_deltaPrim(vof, inormVelVar) = rsign*dvmin;
          }
        }
        else //no faces on high side of low boundary vof
        {
          a_deltaPrim(vof, inormVelVar) = 0.0;
        }
      } //end loop over irregular vofs at boundary
    }
  } //end loop over boundary sides
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void
EBNoFlowIBC::
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
  // No periodic domains allowed.
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
                 CHF_CONST_REAL(m_dx[0]),
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
        Real speed = 0.0; //sqrt(m_gamma*press/dense);


          //blessed cartesian coords.
          a_flux[a_dir](face, inormMomVar) =
            press + isign*dense*unorm*speed;
      }
    }
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void
EBNoFlowIBC::define(const ProblemDomain&  a_domain,
                    const RealVect&       a_dx)
{
  m_domain = a_domain;
  m_dx = a_dx;
  CH_assert(Abs(a_dx[0]-a_dx[1]) < 1.e-9);
  m_isDefined = true;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void
EBNoFlowIBC::
initialize(LevelData<EBCellFAB>& a_conState,
           const EBISLayout& a_ebisl) const
{
  // Shouldn't be called!
  CH_assert(false);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
EBNoFlowIBC::
~EBNoFlowIBC()
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
EBNoFlowIBCFactory::
EBNoFlowIBCFactory(Real a_gamma):
  m_gamma(a_gamma)
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
EBNoFlowIBCFactory::
~EBNoFlowIBCFactory()
{
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
EBPhysIBC*
EBNoFlowIBCFactory::
create() const
{
  return new EBNoFlowIBC(m_gamma);
}
//------------------------------------------------------------------------------

#include "NamespaceFooter.H"
