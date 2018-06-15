#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AdvectPhysics.H"
#include "AdvectPhysicsF_F.H"

#include "LoHiSide.H"

#include "NamespaceHeader.H"

/***/
AdvectPhysics::
AdvectPhysics():GodunovPhysics()
{
  m_isVelSet = false;
}
/***/
AdvectPhysics::
~AdvectPhysics()
{
}

// Compute the maximum wave speed
Real
AdvectPhysics::
getMaxWaveSpeed(const FArrayBox& a_U,
                const Box&       a_box)
{
  CH_assert(isDefined());

  // note that for this advection problem, the max wavespeed
  // is based entirely on the cell-centered velocity
  Real maxSpeed = 0.0;
  Real speed    = 0.0;
  for (int dir = 0; dir < SpaceDim; dir++)
    {
      FArrayBox& vel = *m_cellVelPtr;
      speed = vel.norm(0,dir,1);
      if (maxSpeed < speed) maxSpeed = speed;
    }

  return maxSpeed;
}

// Factory method - this object is its own factory.  It returns a pointer
// to new AdvectPhysics object with the same definition as this object.
GodunovPhysics*
AdvectPhysics::
new_godunovPhysics() const
{
  // Make the new object
  AdvectPhysics* newPtr = new AdvectPhysics();

  // Pass the current initial and boundary condition (IBC) object to this new
  // patch integrator so it can create its own IBC object and define it as it
  // needs to (i.e. new domain and grid spacing if necessary)

  newPtr->setPhysIBC(getPhysIBC());
  // Define the same domain and dx
  newPtr->define(m_domain,m_dx);

  GodunovPhysics* retval = static_cast<GodunovPhysics*>(newPtr);

  // Return the new object
  return retval;
}
/****/
void
AdvectPhysics::
getFlux(FArrayBox&       a_flux,
        const FArrayBox& a_WHalf,
        const int&       a_dir,
        const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_flux.nComp() == a_WHalf.nComp());

  FArrayBox& advVelDir = (*m_advVelPtr)[a_dir];
  // flux is WHalf*advVel
  a_flux.copy(a_WHalf, a_box);
  // componentwise multiplication, since all components of flux
  // get multiplied by 0th componenent of advVel
  for (int comp=0; comp<a_flux.nComp(); comp++)
    {
      a_flux.mult(advVelDir, a_box, 0, comp, 1);
    }
}
/****/
void
AdvectPhysics::
charValues(FArrayBox&       a_lambda,
           const FArrayBox& a_W,
           const int&       a_dir,
           const Box&       a_box)
{

  CH_assert(isDefined());
  CH_assert(m_isVelSet);
  CH_assert(a_lambda.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  CH_assert(m_cellVelPtr != NULL);
  FArrayBox& cellVel = *m_cellVelPtr;

  int nComp = a_lambda.nComp();
  int copyComp = 1;

  for (int comp = 0; comp < nComp; comp++)
  {
    a_lambda.copy(cellVel,a_box,a_dir,a_box,comp,copyComp);
  }
}
/***/
// Compute a Riemann problem and generate fluxes at the faces
void
AdvectPhysics::
postNormalPred(FArrayBox&       a_dWMinus,
               FArrayBox&       a_dWPlus,
               const FArrayBox& a_W,
               const Real&      a_dt,
               const Real&      a_dx,
               const int&       a_dir,
               const Box&       a_box)
{
  CH_assert(isDefined());
  int numPrim = numPrimitives();

  CH_assert(a_dWPlus.nComp() == numPrim);
  CH_assert(a_dWMinus.nComp() == numPrim);
  CH_assert(a_dWPlus.box().contains(a_box));
  CH_assert(m_advVelPtr != NULL);
  FluxBox& advVel = *m_advVelPtr;
  FArrayBox& advVelDir = advVel[a_dir];

  advVelDir.shiftHalf(a_dir,1);
  FORT_POSTNORMALSOURCE(
                CHF_FRA(a_dWPlus),
                CHF_FRA(a_dWMinus),
                CHF_CONST_FRA(a_W),
                CHF_CONST_FRA1(advVelDir,0),
                CHF_CONST_REAL(a_dt),
                CHF_CONST_REAL(a_dx),
                CHF_CONST_INT(a_dir),
                CHF_BOX(a_box));
  advVelDir.shiftHalf(a_dir,-1);
}
void
AdvectPhysics::
riemann(FArrayBox&       a_WGdnv,
        const FArrayBox& a_WLeft,
        const FArrayBox& a_WRight,
        const FArrayBox& a_W,
        const Real&      a_time,
        const int&       a_dir,
        const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_WGdnv.box().contains(a_box));

  CH_assert(m_advVelPtr != NULL);
  FluxBox& advVel = *m_advVelPtr;
  // CH_assert(advVel.box().contains(a_box));

  // Get the numbers of relevant variables
  int numPrim = numPrimitives();

  CH_assert(a_WGdnv .nComp() == numPrim);
  CH_assert(a_WLeft .nComp() == numPrim);
  CH_assert(a_WRight.nComp() == numPrim);

  // Cast away "const" inputs so their boxes can be shifted left or right
  // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&)a_WLeft;
  FArrayBox& shiftWRight = (FArrayBox&)a_WRight;

  // Solution to the Riemann problem

  // Shift the left and right primitive variable boxes 1/2 cell so they are
  // face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);

  CH_assert(shiftWLeft .box().contains(a_box));
  CH_assert(shiftWRight.box().contains(a_box));

  FArrayBox& advVelDir = advVel[a_dir];
  CH_assert(advVelDir.box().contains(a_box));

  // Riemann solver computes WGdnv all edges that are not on the physical
  // boundary.
  FORT_RIEMANNF(CHF_FRA(a_WGdnv),
                CHF_CONST_FRA(shiftWLeft),
                CHF_CONST_FRA(shiftWRight),
                CHF_CONST_FRA1(advVelDir,0),
                CHF_CONST_INT(a_dir),
                CHF_BOX(a_box));

  // Call boundary Riemann solver (note: periodic BC's are handled there).
  m_bc->primBC(a_WGdnv,shiftWLeft ,a_W,a_dir,Side::Hi,a_time);
  m_bc->primBC(a_WGdnv,shiftWRight,a_W,a_dir,Side::Lo,a_time);

  // Shift the left and right primitive variable boxes back to their original
  // position
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}

/******/
void
AdvectPhysics::
quasilinearUpdate(FArrayBox&       a_AdWdx,
                  const FArrayBox& a_WHalf,
                  const FArrayBox& a_W,
                  const Real&      a_scale,
                  const int&       a_dir,
                  const Box&       a_box)
{
  CH_assert(isDefined());

  CH_assert(a_AdWdx.box().contains(a_box));
  CH_assert(a_W    .box().contains(a_box));

  // Get the numbers of relevant variables
  //int numPrim = numPrimitives();

  const FArrayBox& cellVel = *m_cellVelPtr;

  FORT_QUASILINEARUPDATE(CHF_FRA(a_AdWdx),
                         CHF_CONST_FRA(a_WHalf),
                         CHF_CONST_FRA1(cellVel,a_dir),
                         CHF_CONST_REAL(a_scale),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(a_box));
}
/****/
// Compute the primitive variables from the conserved variables
void
AdvectPhysics::
consToPrim(FArrayBox&       a_W,
           const FArrayBox& a_U,
           const Box&       a_box)
{
  CH_assert(isDefined());
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  // since the primitive and conserved variables are the same in
  // this set of equations, even I can do this one...
  a_W.copy(a_U, a_box);
}

#include "NamespaceFooter.H"
