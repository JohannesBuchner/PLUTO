#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include <string>
using std::string;

#include "LoHiSide.H"

#include "GodunovUtilities.H"
#include "LoHiCenter.H"
#include "GodunovPhysics.H"
#include "GodunovUtilitiesF_F.H"
#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Flag everything as not defined or set
GodunovUtilities::GodunovUtilities()
{
  m_highOrderLimiter = false;
  m_isDefined = false;
}

//////////////////////////////////////////////////////////////////////////////
GodunovUtilities::~GodunovUtilities()
{
}

//////////////////////////////////////////////////////////////////////////////
// Define this object and the boundary condition object
void GodunovUtilities::define(const ProblemDomain& a_domain,
                              const Real&          a_dx)
{
  // Store the domain and grid spacing
  m_domain    = a_domain;
  m_dx        = a_dx;
  m_isDefined = true;
}

//////////////////////////////////////////////////////////////////////////////
// Compute the flattening coefficients from the primitive variables
void GodunovUtilities::computeFlattening(FArrayBox&       a_flattening,
                                         const FArrayBox& a_W,
                                         const Interval&  a_velInt,
                                         const int&       a_presInd,
                                         const Real&      a_smallPres,
                                         const int&       a_bulkModulusInd,
                                         const Box&       a_box)
{
  CH_assert(m_isDefined);
  CH_assert(a_W.box().contains(a_box));

  // The current direction
  int idir;

  // The directional flattening coefficients
  FArrayBox zetaDir(a_box,SpaceDim);

  // The divergence of the velocity
  FArrayBox dVel(a_box,SpaceDim);

  // The interval of the primitive variables corresponding to the velocity
  Interval velInterval= a_velInt;
  int v0index = velInterval.begin();

  // Get the directional flattening coefficients in each direction
  for (idir = 0; idir < SpaceDim; idir++)
  {
    // A box one larger (in direction "idir") than the final result box
    Box box1 = a_box;
    box1.grow(idir,1);

    // A box two larger (in direction "idir") than the final result box
    Box box2 = a_box;
    box2.grow(idir,2);

    // A box three larger (in direction "idir") than the final result box
    Box box3 = a_box;
    box3.grow(idir,3);

    // The primitive variables need to be defined over the largest box.
    // This assert may fail when a_W has been defined on a box that's grown
    // by 3 but then intersected with the domain.
    // CH_assert(a_W.box().contains(box3));

    // *** Compute delta1p, the first differences in "pressure"

    // Compute where centered differences can be used and where one-sided
    // differences need to be used.  The data used for the differences is
    // defined on cell-centered "box3".
    Box loBox, hiBox, centerBox, entireBox;
    int hasLo, hasHi;

    loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
               box3, m_domain, idir);

    // Set delta1p to cell-centered difference of cell-centered pressure
    // in direction idir.  Result is on entireBox.
    // For i in centerBox, delta1p[i] = (pressure[i+1] - pressure[i-1])/2.
    // For i in loBox, delta1p[i] = pressure[i+1] - pressure[i].
    // For i in hiBox, delta1p[i] = pressure[i] - pressure[i-1].

    FArrayBox delta1p(entireBox,1);
    FORT_GETGRADF(CHF_FRA1(delta1p,0),
                  CHF_CONST_FRA1(a_W,a_presInd),
                  CHF_CONST_INT(idir),
                  CHF_BOX(loBox),
                  CHF_CONST_INT(hasLo),
                  CHF_BOX(hiBox),
                  CHF_CONST_INT(hasHi),
                  CHF_BOX(centerBox));

    // *** Compute delta2p, the second differences in "pressure"

    // Compute where centered differences can be used and where one-sided
    // differences need to be used.  The data used for the differences is
    // defined on cell-centered "box2"
    loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
               box2, m_domain, idir);

    // Set delta2p to twice the cell-centered average of cell-centered delta1p
    // in direction idir.  Result is on entireBox.
    // For i in centerBox, delta2p[i] = delta1p[i-1] + delta1p[i+1].
    // For i in loBox, delta2p[i] = delta1p[i] + delta1p[i+1].
    // For i in hiBox, delta2p[i] = delta1p[i-1] + delta1p[i].

    FArrayBox delta2p(entireBox,1);
    FORT_GETDPTWOF(CHF_FRA1(delta2p,0),
                   CHF_CONST_FRA1(delta1p,0),
                   CHF_CONST_INT(idir),
                   CHF_BOX(loBox),
                   CHF_CONST_INT(hasLo),
                   CHF_BOX(hiBox),
                   CHF_CONST_INT(hasHi),
                   CHF_BOX(centerBox));

    // *** Compute bulkMin, a 3-way minimum of the "bulk modulus"

    // Same loBox, centerBox, hiBox, and entireBox as delta2p.

    // Set bulkMin to minimum of cell-centered bulk in neighborhood.
    // Result is on entireBox.
    // For i in centerBox, bulkMin[i] = min{bulk[i-1], bulk[i], bulk[i+1]}.
    // For i in loBox, bulkMin[i] = min{bulk[i], bulk[i+1]}.
    // For i in hiBox, bulkMin[i] = min{bulk[i-1], bulk[i]}.

    FArrayBox bulkMin(entireBox,1);
    int bulkIndex = a_bulkModulusInd;
    FORT_MIN3PTSF(CHF_FRA1(bulkMin,0),
                  CHF_CONST_FRA1(a_W,bulkIndex),
                  CHF_CONST_INT(idir),
                  CHF_BOX(loBox),
                  CHF_CONST_INT(hasLo),
                  CHF_BOX(hiBox),
                  CHF_CONST_INT(hasHi),
                  CHF_BOX(centerBox));

    // *** Use the first and second differences normalized by the
    // 3-way minimum, bulkMin, to generate zetaTwiddleDir,
    // flattening coefficients in direction idir.

    // Same loBox, centerBox, hiBox, and entireBox as delta2p and bulkMin.

    // if abs(delta1p[i] / bulkMin[i]) < 0.33,
    //    zetaTwiddleDir[i] = 1.;
    // else
    //    ratio[i] = abs(delta1p[i]) / max{abs(delta2p[i]), smallp};
    //    if ratio[i] <= 0.75, zetaTwiddleDir[i] = 1.;
    //    elseif ratio[i] >= 0.85, zetaTwiddleDir[i] = 0.;
    //    else zetaTwiddleDir[i] = 1. - (ratio[i] - 0.75) / (0.85 - 0.75).

    FArrayBox zetaTwiddleDir(entireBox,1);
    FORT_GETFLATF(CHF_FRA1(zetaTwiddleDir,0),
                  CHF_CONST_FRA1(delta1p,0),
                  CHF_CONST_FRA1(delta2p,0),
                  CHF_CONST_REAL(a_smallPres),
                  CHF_CONST_FRA1(bulkMin,0),
                  CHF_BOX(entireBox));

    // *** Compute component idir of zetaDir, a 3-way minimum
    // of the directional flattening coefficients, zetaTwiddleDir.

    // Compute where centered differences can be used and where one sided
    // differences need to be used.  The data used for the differences is
    // defined on cell-centered "box1"
    loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
               box1, m_domain, idir);

    // Set zetaDir to minimum of cell-centered zetTwiddleDir in neighborhood.
    // Result is on entireBox.
    // For i in centerBox, zetaDir[i] = min{zTD[i-1], zTD[i], zTD[i+1]}.
    // For i in loBox, zetaDir[i] = min{zTD[i], zTD[i+1]}.
    // For i in hiBox, zetaDir[i] = min{zTD[i-1], zTD[i]}.
    FORT_MIN3PTSF(CHF_FRA1(zetaDir,idir),
                  CHF_CONST_FRA1(zetaTwiddleDir,0),
                  CHF_CONST_INT(idir),
                  CHF_BOX(loBox),
                  CHF_CONST_INT(hasLo),
                  CHF_BOX(hiBox),
                  CHF_CONST_INT(hasHi),
                  CHF_BOX(centerBox));

    // *** Compute component idir of dVel, the divergence of the velocity

    // Set dVel to cell-centered difference of cell-centered pressure
    // in direction idir.  Result is on entireBox.
    // For i in centerBox, dVel[i] = (velocity[i+1] - velocity[i-1])/2.
    // For i in loBox, dVel[i] = velocity[i+1] - velocity[i].
    // For i in hiBox, dVel[i] = velocity[i] - velocity[i-1].
    FORT_GETGRADF(CHF_FRA1(dVel,idir),
                  CHF_CONST_FRA1(a_W,v0index+idir),
                  CHF_CONST_INT(idir),
                  CHF_BOX(loBox),
                  CHF_CONST_INT(hasLo),
                  CHF_BOX(hiBox),
                  CHF_CONST_INT(hasHi),
                  CHF_BOX(centerBox));
  }

  // *** Where divergence of velocity is negative, set a_flattening
  // to minimum of directional flattening coefficients zetaDir at that point

  // sum{dVel[i]} is divergence of dVel[i] (multiplied by mesh spacing).
  //
  // if sum{dVel[i]} >= 0,
  //    a_flattening[i] = 1.;
  // else
  //    a_flattening[i] = min{zetaDir[i]}.

  // At each point, set the flattening coefficient to the minimum of all
  // the directional flattening coefficients if the divergence of the velocity
  // is negative, otherwise set it to 1 (no flattening).
  FORT_MINFLATF(CHF_FRA1(a_flattening,0),
                CHF_CONST_FRA(zetaDir),
                CHF_CONST_FRA(dVel),
                CHF_BOX(a_box));
}

//////////////////////////////////////////////////////////////////////////////
void GodunovUtilities::applyFlattening(      FArrayBox& a_dW,
                                       const FArrayBox& a_flat,
                                       const Box&       a_box)
{
  CH_assert(a_dW.box().contains(a_box));
  CH_assert(a_flat.box().contains(a_box));

  int numSlopes = a_dW.nComp();

  for (int islope = 0;islope < numSlopes;islope++)
  {
    a_dW.mult(a_flat,a_box,0,islope);
  }
}

//////////////////////////////////////////////////////////////////////////////
void GodunovUtilities::vanLeerSlopes(FArrayBox&       a_dW,
                                     const FArrayBox& a_W,
                                     const int&       a_numSlopes,
                                     const bool&      a_useLimiting,
                                     const int&       a_dir,
                                     const Box&       a_box)
{
  // A box one larger (in direction "a_dir") than the final result box
  // 2 Sep 2008:  For vanLeerSlopesExtPreserving, expand by 2?  17 Sep 2008
  Box box1 = a_box;
  int ghostbox1 = 1; // FIX, 19 sep 2008 (m_highOrderLimiter) ? 2 : 1;
  // int ghostbox1 = (m_highOrderLimiter) ? 2 : 1;
  box1.grow(a_dir, ghostbox1);

  // Compute where centered differences can be used and where one-sided
  // differences need to be used.
  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
             box1,m_domain,a_dir);

  // Compute 2nd order slopes - including one-sided differences
  FArrayBox dWMinus(entireBox,a_numSlopes);
  FArrayBox dWPlus (entireBox,a_numSlopes);

  // We calculate a_dW and dWMinus and dWPlus on centerBox,
  // and we need a_W on centerBox grown by 1 in direction a_dir.
  slopes(a_dW, dWMinus, dWPlus, a_W, a_numSlopes, a_dir,
         loBox, hiBox, centerBox, entireBox, hasLo, hasHi);

  // Apply the slope limiter if requested
  if (a_useLimiting)
  {
    // Apply slopeLimiter only on centerBox; elsewhere, a_dW is unchanged.

    // 2 Sep 2008:  replace slopeLimiter with slopeLimiterExtPreserving

    if (m_highOrderLimiter)
      {
        // Modifies a_dW on centerBox,
        // and needs dWMinus on centerBox including shift down by 1 in direction a_dir,
        // and needs dWPlus on centerBox including shift up by 1 in direction a_dir.
        slopeLimiterExtPreserving(a_dW, dWMinus, dWPlus, a_numSlopes, centerBox, a_dir);
      }
    else
      {
        slopeLimiter(a_dW, dWMinus, dWPlus, a_numSlopes, centerBox);
      }
  }
}

//////////////////////////////////////////////////////////////////////////////
void GodunovUtilities::fourthOrderSlopes(FArrayBox&       a_dW,
                                         const FArrayBox& a_W,
                                         const FArrayBox& a_dWvL,
                                         const int&       a_numSlopes,
                                         const int&       a_dir,
                                         const Box&       a_box)
{
  // Number of slopes to compute
  int numSlope = a_numSlopes;

  CH_assert(a_dW.nComp() == numSlope);
  CH_assert(a_W.nComp() >= numSlope);

  // A box one larger (in direction "a_dir") than the final result box
  Box box1 = a_box;
  box1.grow(a_dir,1);

  // Compute where centered differences can be used and where one sided
  // differences need to be used.
  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
             box1,m_domain,a_dir);

  CH_assert(a_dW.box().contains(entireBox));

  FORT_FOURTHSLOPEDIFFSF(CHF_FRA(a_dW),
                         CHF_CONST_FRA(a_W),
                         CHF_CONST_FRA(a_dWvL),
                         CHF_CONST_INT(numSlope),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(loBox),
                         CHF_CONST_INT(hasLo),
                         CHF_BOX(hiBox),
                         CHF_CONST_INT(hasHi),
                         CHF_BOX(centerBox));
}

//////////////////////////////////////////////////////////////////////////////
void GodunovUtilities::oneSidedDifferences(FArrayBox&       a_dWMinus,
                                           FArrayBox&       a_dWPlus,
                                           const FArrayBox& a_W,
                                           const int&       a_dir,
                                           const Box&       a_box)
{
  int numSlopes = a_dWMinus.nComp();

  // A box one larger (in direction "a_dir") than the final result box
  Box box1 = a_box;
  box1.grow(a_dir,1);

  // Compute where centered differences can be used and where one sided
  // differences need to be used.
  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
             box1,m_domain,a_dir);

  // Compute 2nd order slopes - including one sided differences
  FArrayBox deltaWC(entireBox, numSlopes);

  slopes(deltaWC, a_dWMinus, a_dWPlus, a_W, numSlopes, a_dir,
         loBox, hiBox, centerBox, entireBox, hasLo, hasHi);
}

//////////////////////////////////////////////////////////////////////////////
void GodunovUtilities::slopes(FArrayBox&       a_dWCent,
                              FArrayBox&       a_dWMinus,
                              FArrayBox&       a_dWPlus,
                              const FArrayBox& a_W,
                              const int&       a_numSlopes,
                              const int&       a_dir,
                              const Box&       a_loBox,
                              const Box&       a_hiBox,
                              const Box&       a_centerBox,
                              const Box&       a_entireBox,
                              const int&       a_hasLo,
                              const int&       a_hasHi)

{
  CH_assert(a_dWCent .nComp() == a_numSlopes);
  CH_assert(a_dWMinus.nComp() == a_numSlopes);
  CH_assert(a_dWPlus .nComp() == a_numSlopes);
  CH_assert(a_W.nComp() >= a_numSlopes);

  CH_assert(a_dWCent .box().contains(a_entireBox));
  CH_assert(a_dWMinus.box().contains(a_entireBox));
  CH_assert(a_dWPlus .box().contains(a_entireBox));
  CH_assert(a_W.box().contains( a_entireBox));

  CH_assert((a_dir >= 0 )&& (a_dir < SpaceDim));

  FORT_SECONDSLOPEDIFFSF(CHF_FRA(a_dWCent),
                         CHF_FRA(a_dWMinus),
                         CHF_FRA(a_dWPlus),
                         CHF_CONST_FRA(a_W),
                         CHF_CONST_INT(a_numSlopes),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(a_loBox),
                         CHF_CONST_INT(a_hasLo),
                         CHF_BOX(a_hiBox),
                         CHF_CONST_INT(a_hasHi),
                         CHF_BOX(a_centerBox));
}

//////////////////////////////////////////////////////////////////////////////
// Apply a van Leer limiter directly to the slopes.
void GodunovUtilities::slopeLimiter(FArrayBox&       a_dW,
                                    const FArrayBox& a_dWLeft,
                                    const FArrayBox& a_dWRigh,
                                    const int&       a_numSlopes,
                                    const Box&       a_box)
{
  CH_assert(m_isDefined);
  CH_assert(a_dW.nComp()     == a_numSlopes);
  CH_assert(a_dWLeft.nComp() == a_numSlopes);
  CH_assert(a_dWRigh.nComp() == a_numSlopes);
  CH_assert(a_dW.box().contains(a_box));
  CH_assert(a_dWLeft.box().contains(a_box));
  CH_assert(a_dWRigh.box().contains(a_box));

  FORT_VANLEERLIMITERF(CHF_FRA(a_dW),
                       CHF_CONST_FRA(a_dWLeft),
                       CHF_CONST_FRA(a_dWRigh),
                       CHF_CONST_INT(a_numSlopes),
                       CHF_BOX(a_box));
}

//////////////////////////////////////////////////////////////////////////////
// Apply an extremum-preserving van Leer limiter directly to the slopes.
void GodunovUtilities::slopeLimiterExtPreserving(FArrayBox&       a_dW,
                                                 const FArrayBox& a_dWLeft,
                                                 const FArrayBox& a_dWRigh,
                                                 const int&       a_numSlopes,
                                                 const Box&       a_box,
                                                 const int&       a_dir)
{
  CH_assert(m_isDefined);
  CH_assert(a_dW.nComp()     == a_numSlopes);
  CH_assert(a_dWLeft.nComp() == a_numSlopes);
  CH_assert(a_dWRigh.nComp() == a_numSlopes);
  CH_assert(a_dW.box().contains(a_box));
  CH_assert(a_dWLeft.box().contains(a_box));
  CH_assert(a_dWRigh.box().contains(a_box));

  // Compute where centered differences can be used and where one-sided
  // differences need to be used.
  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;

  loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
             a_box, m_domain, a_dir);

  if (hasLo)
    {
      FORT_VANLEERLIMITERF(CHF_FRA(a_dW),
                           CHF_CONST_FRA(a_dWLeft),
                           CHF_CONST_FRA(a_dWRigh),
                           CHF_CONST_INT(a_numSlopes),
                           CHF_BOX(loBox));
    }
  if (hasHi)
    {
      FORT_VANLEERLIMITERF(CHF_FRA(a_dW),
                           CHF_CONST_FRA(a_dWLeft),
                           CHF_CONST_FRA(a_dWRigh),
                           CHF_CONST_INT(a_numSlopes),
                           CHF_BOX(hiBox));
    }
  if (!centerBox.isEmpty())
    {
      // Modifies a_dW on centerBox,
      // and needs a_dWLeft on centerBox including shift down by 1 in direction a_dir,
      // and needs a_dWRigh on centerBox including shift up by 1 in direction a_dir.
      FORT_EXTPRESERVINGVANLEERLIMITERF(CHF_FRA(a_dW),
                                        CHF_CONST_FRA(a_dWLeft),
                                        CHF_CONST_FRA(a_dWRigh),
                                        CHF_CONST_INT(a_numSlopes),
                                        CHF_CONST_INT(a_dir),
                                        CHF_BOX(centerBox));
    }
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

//////////////////////////////////////////////////////////////////////////////
// Piecewise linear normal predictor. Computes increments in the
// characteristic amplitudes.
void GodunovUtilities::PLMNormalPred(FArrayBox&       a_dWCharLo,
                                     FArrayBox&       a_dWCharHi,
                                     const FArrayBox& a_dWChar,
                                     const FArrayBox& a_Lambda,
                                     const Real&      a_dtbydx,
                                     const Box&       a_box)
{
  int numPrim = a_dWChar.nComp();
  CH_assert(a_dWCharLo.nComp() == numPrim);
  CH_assert(a_dWCharHi.nComp() == numPrim);
  CH_assert(a_Lambda.nComp() == numPrim);
  CH_assert(a_dWCharLo.box().contains(a_box));
  CH_assert(a_dWCharHi.box().contains(a_box));
  CH_assert(a_dWChar.box().contains(a_box));
  CH_assert(a_Lambda.box().contains(a_box));


  FORT_PLMNORMALPREDF(CHF_FRA(a_dWCharLo),
                      CHF_FRA(a_dWCharHi),
                      CHF_CONST_FRA(a_dWChar),
                      CHF_CONST_FRA(a_Lambda),
                      CHF_CONST_REAL(a_dtbydx),
                      CHF_CONST_INT(numPrim),
                      CHF_BOX(a_box));
}

//////////////////////////////////////////////////////////////////////////////
void GodunovUtilities::PPMFaceValues(FArrayBox&            a_WFace,
                                     const FArrayBox&      a_W,
                                     const int&            a_numSlopes,
                                     const bool&           a_useLimiting,
                                     const int&            a_dir,
                                     const Box&            a_box,
                                     const Real&           a_time,
                                     const GodunovPhysics* a_physPtr)
{
  // a_box is the face-centered box on which a_WFace is computed.
  CH_assert(a_WFace.box().contains(a_box));

  if (m_highOrderLimiter)
    {
      // petermc, 7 Jan 2010: changed from a_box to grown a_box.
      Box box1cells = grow(a_box, BASISV(a_dir));
      box1cells.enclosedCells();

      Box loFaces, nextLoFaces;
      Box hiFaces, nextHiFaces;
      Box centerFaces, innerCenterFaces, entireFaces;
      int hasLoFaces, hasHiFaces;
      loHiCenterFace4(loFaces, nextLoFaces, hasLoFaces,
                      hiFaces, nextHiFaces, hasHiFaces,
                      centerFaces, innerCenterFaces, entireFaces,
                      box1cells, m_domain, a_dir);

      // For i-e/2 in innerCenterFaces, set a_WFace[i-e/2] from a_W[i-2e:i+e].
      // If hasLoFaces:
      // For i-e/2 in loFaces, set a_WFace[i-e/2] from a_W[i:i+3e].
      //           in nextLoFaces,                from a_W[i-e:i+2e].
      // If hasHiFaces:
      // For i-e/2 in hiFaces, set a_WFace[i-e/2] from a_W[i-4e:i-e].
      //           in nextHiFaces,                from a_W[i-3e:i].
      FORT_FOURTHINTERPFACES(CHF_FRA(a_WFace),
                             CHF_CONST_FRA(a_W),
                             CHF_CONST_INT(a_numSlopes),
                             CHF_CONST_INT(a_dir),
                             CHF_BOX(loFaces),
                             CHF_BOX(nextLoFaces),
                             CHF_CONST_INT(hasLoFaces),
                             CHF_BOX(hiFaces),
                             CHF_BOX(nextHiFaces),
                             CHF_CONST_INT(hasHiFaces),
                             CHF_BOX(innerCenterFaces));

      // WAS if (a_useLimiting) call face limiter;
      // this removed by petermc, 7 Oct 2010

      // dummy statement in order to get around gdb bug
      int dummy_unused = 0; dummy_unused = 0;
    }
  else // !m_highOrderLimiter :  this is the old method
    {
      // A box one larger (in direction "a_dir") than the final result box
      // petermc, 14 Aug 2009:  first grow(), then enclosedCells(),
      // rather than reverse order, so you don't end up with empty box1cells.
      Box box1cells(a_box);
      int ghostbox1 = 1;
      box1cells.grow(a_dir, ghostbox1);
      box1cells.enclosedCells();

      FArrayBox dW(box1cells, a_numSlopes);
      vanLeerSlopes(dW, a_W, a_numSlopes, a_useLimiting, a_dir, box1cells);

      if (a_physPtr != NULL)
        {
          Real a_time = 0.0;
          a_physPtr->getPhysIBC()->setBdrySlopes(dW, a_W, a_dir, a_time);
        }

      Box loBox,hiBox,centerBox,entireBox;
      int hasLo,hasHi;

      loHiCenterFace(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                     box1cells, m_domain, a_dir);

      // a_Wface[i-e/2] = (a_W[i-e] + dW[i-e]/3)/2 + (a_W[i] - dW[i]/3)/2
      FORT_PPMFACEVALUESF(CHF_FRA(a_WFace),
                          CHF_CONST_FRA(a_W),
                          CHF_CONST_FRA(dW),
                          CHF_CONST_INT(a_numSlopes),
                          CHF_CONST_INT(a_dir),
                          CHF_BOX(loBox),
                          CHF_CONST_INT(hasLo),
                          CHF_BOX(hiBox),
                          CHF_CONST_INT(hasHi),
                          CHF_BOX(centerBox));
    }
}

//////////////////////////////////////////////////////////////////////////////
void GodunovUtilities::PPMLimiter(FArrayBox& a_dWMinus,
                                  FArrayBox& a_dWPlus,
                                  const FArrayBox& a_W,
                                  const int& a_numSlopes,
                                  const int& a_dir,
                                  const Box& a_box)
{
  // Called by PatchGodunov::PPMNormalPred,
  // which is called by PatchGodunov::computeWHalf.
  // a_dWMinus[i] = WFace[i - e/2] - a_W[i]
  // a_dWPlus[i] = WFace[i + e/2] - a_W[i]
  // where e is unit vector in dimension a_dir.

  if (m_highOrderLimiter)
    {
      // this option added by petermc, 5 Sep 2008
      // Will need to recalculate some D^2's.

      // We calculate d2Wfcf on a_box,
      // and we need a_dWMinus on a_box,
      // and we need a_dWPlus on a_box.
      // d2Wfcf[i] = 6 * (a_dWMinus[i] + a_dWPlus[i])
      //           = 6 * (thisFaceWDir[i-e/2] - a_cellW[i] +
      //                  thisFaceWDir[i+e/2] - a_cellW[i])
      FArrayBox d2Wfcf(a_box, a_numSlopes);
      d2Wfcf.copy(a_dWMinus);
      d2Wfcf.plus(a_dWPlus, 0, 0, a_numSlopes);
      d2Wfcf *= 6.;

      // petermc, 21 Sep 2010:
      // In order to get a_dWMinus and a_dWPlus on a_box,
      // we need d2W on a_box grown by 3 in a_dir directions.
      Box box3 = a_box;
      box3.grow(a_dir, 3);

      Box loBox, hiBox, centerBox, entireBox;
      int hasLo, hasHi;
      loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                 box3, m_domain, a_dir);

      FArrayBox d2W(entireBox, a_numSlopes);
      // On centerBox, use 3-point stencil for d2W;
      // on loBox and hiBox, copy result from neighbor.
      // In the case of centerBox == entireBox == box1,
      // where box1 is a_box grown by 1 in dimension a_dir:
      // We calculate d2W on a_box grown by 2 (was 1) in dimension a_dir,
      // and we need a_W on a_box grown by 3 (was 2) in dimension a_dir.

      // petermc, 21 Sep 2010, changed layer of d2W from 2 to 1,
      // and of a_W from 3 to 2.
      FORT_GETSECONDDIFF(CHF_FRA(d2W),
                         CHF_CONST_FRA(a_W),
                         CHF_CONST_INT(a_numSlopes),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(loBox),
                         CHF_CONST_INT(hasLo),
                         CHF_BOX(hiBox),
                         CHF_CONST_INT(hasHi),
                         CHF_BOX(centerBox));

      Box box1 = a_box;
      box1.grow(a_dir, 1);
      Box nextLoBox, nextHiBox, innerCenterBox;
      loHiCenter5(loBox, nextLoBox, hasLo,
                  hiBox, nextHiBox, hasHi,
                  centerBox, innerCenterBox, entireBox,
                  box1, m_domain, a_dir);
      CH_assert(entireBox == a_box);

      Real limitC = 1.25;
      Real eps = 1.0e-12;
      Real c3 = 0.1;
      FORT_CHECKCUBICLIMITERF(// <W>_(i-1/2) - <W>_i, on a_box
                              CHF_FRA(a_dWMinus),
                              // <W>_(i+1/2) - <W>_i, on a_box
                              CHF_FRA(a_dWPlus),
                              CHF_CONST_FRA(a_W),
                              CHF_CONST_FRA(d2W),
                              CHF_CONST_FRA(d2Wfcf),
                              CHF_CONST_INT(a_numSlopes),
                              CHF_CONST_INT(a_dir),
                              CHF_BOX(loBox),
                              CHF_BOX(nextLoBox),
                              CHF_CONST_INT(hasLo),
                              CHF_BOX(hiBox),
                              CHF_BOX(nextHiBox),
                              CHF_CONST_INT(hasHi),
                              CHF_BOX(innerCenterBox),
                              CHF_CONST_REAL(limitC),
                              CHF_CONST_REAL(c3),
                              CHF_CONST_REAL(eps));
      // dummy statement in order to get around gdb bug
      int dummy_unused = 0; dummy_unused = 0;
    }
  else
    {
      FORT_PPMLIMITERF(// <W>_(i-1/2) - <W>_i, on a_box
                       CHF_FRA(a_dWMinus),
                       // <W>_(i+1/2) - <W>_i, on a_box
                       CHF_FRA(a_dWPlus),
                       CHF_CONST_INT(a_numSlopes),
                       CHF_BOX(a_box));
    }
}

//////////////////////////////////////////////////////////////////////////////
void GodunovUtilities::PPMNormalPred(FArrayBox&       a_dWMinus,
                                     FArrayBox&       a_dWPlus,
                                     const FArrayBox& a_Lambda,
                                     const Real&      a_dtbydx,
                                     const int&       a_numSlopes,
                                     const Box&       a_box)
{
  FORT_PPMNORMALPREDF(CHF_FRA(a_dWMinus),
                      CHF_FRA(a_dWPlus),
                      CHF_CONST_FRA(a_Lambda),
                      CHF_CONST_REAL(a_dtbydx),
                      CHF_CONST_INT(a_numSlopes),
                      CHF_BOX(a_box));
}

//////////////////////////////////////////////////////////////////////////////
// Compute a face centered divergence of the velocity
void GodunovUtilities::divVel(FArrayBox&       a_divVel,
                              const FArrayBox& a_W,
                              const Interval&  a_velInt,
                              const int&       a_dir,
                              const Box&       a_box)
{
  // Quick description of this function:
  // Let v[d] denote d-th component of velocity (from a_velInt within a_W).
  // Then this function returns, on face-centered a_box,
  // a_divVel[i + e[d]/2] := v[d][i+e[d]] - v[d][i] +
  //                 sum_{d' ~= d} ( Dv[d'][i+e[d]] + Dv[d'][i] ) / 2
  // where Dv[d'][i] = (v[d'][i+e[d']] - v[d'][i-e[d']])/2  in center
  //                or  v[d'][i+e[d']] - v[d'][i]  on low  d'-side
  //                or  v[d'][i] - v[d'][i-e[d']]  on high d'-side
  CH_assert((a_dir >= 0 )&& (a_dir < SpaceDim));
  CH_assert(a_divVel.box().contains(a_box));

  Box dveltanBox = a_box;
  dveltanBox.enclosedCells(a_dir);
  dveltanBox.grow(1);
  dveltanBox &= m_domain;

  // First, we need to calculate the directional derivatives of
  // the tangential components of velocity at the cell-centers.
  // want to make sure numTanComps is > 0 even in 1D, since
  // we pass it in to the divergence (although it isn't used there)
  int numTanComps = std::max(1, SpaceDim-1);
  FArrayBox dveltan(dveltanBox,numTanComps);

  // Get the interval of the primitive variables corresponding to the velocity
  int v0index = a_velInt.begin();

  // Go through the tangential directions
  for (int i = 0, dir = 0; dir < SpaceDim; ++dir)
  {
    if (dir != a_dir)
    {
      // This velocity component is tangential.  Build the box in which
      // d(v[dir])/d(x[dir]) is to be computed.
      Box primBox = a_box;
      primBox.enclosedCells(a_dir);
      primBox.grow(dir,1).grow(a_dir,1);

      // Compute where centered differences can be used and where one sided
      // differences need to be used.  The data used for the differences is
      // defined on "dveltanBox"
      Box gradBox,hiBox,loBox,centerBox;
      int hasLo,hasHi;

      loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,gradBox,
                 primBox,m_domain,dir);

      // Compute d(v[dir])/d(x[dir]).
      FORT_GETGRADF(CHF_FRA1(dveltan,i),
                    CHF_CONST_FRA1(a_W,v0index+dir),
                    CHF_CONST_INT(dir),
                    CHF_BOX(loBox),
                    CHF_CONST_INT(hasLo),
                    CHF_BOX(hiBox),
                    CHF_CONST_INT(hasHi),
                    CHF_BOX(centerBox));
      i++;
    }
  }

  // Now, we can calculate the divergence of the normal velocity
  // at the center normal-direction faces. To do this, we determine
  // the faces at which we have sufficient data to compute centered
  // estimates of h*(div(u)). At the remaining faces. i.e. those
  // corresponding to the physical boundaries, we use zeroth-order
  // extrapolation.

  Box divBox = a_box;
  divBox.enclosedCells(a_dir);
  divBox.grow(a_dir,1);

  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  loHiCenterFace(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                 divBox,m_domain,a_dir);

  // All of the boxes computed above are shifted so as to be cell-centered,
  // with the index of the cell center being identified with the low face.
  // We then shift a_divVel to be compatible with that convention on input,
  // and undo the shift on output.  Basically, everything is made compatible
  // a_W (which is cell-centered).

  loBox.shiftHalf(a_dir,1);
  centerBox.shiftHalf(a_dir,1);
  hiBox.shiftHalf(a_dir,1);

  a_divVel.shiftHalf(a_dir,1);

  FORT_DIVUEDGEF(CHF_FRA1(a_divVel,0),
                 CHF_CONST_FRA1(a_W,v0index+a_dir),
                 CHF_CONST_FRA(dveltan),
                 CHF_CONST_INT(a_dir),
                 CHF_BOX(loBox),
                 CHF_CONST_INT(hasLo),
                 CHF_BOX(hiBox),
                 CHF_CONST_INT(hasHi),
                 CHF_BOX(centerBox));

  a_divVel.shiftHalf(a_dir,-1);
}

//////////////////////////////////////////////////////////////////////////////
// Increment fluxes with artificial viscosity.
void GodunovUtilities::artificialViscosity(FArrayBox&       a_F,
                                           const FArrayBox& a_U,
                                           const FArrayBox& a_divVel,
                                           const Real&      a_coeff,
                                           const int&       a_dir,
                                           const Box&       a_box)
{
  CH_assert((a_dir >= 0 )&& (a_dir < SpaceDim));
  CH_assert(a_F.box().contains(a_box));
  CH_assert(a_divVel.box().contains(a_box));
  CH_assert(a_coeff >= 0.0);

  // Note: a_box is face centered in the a_dir direction and is set up for
  // updating a_F

  // Remove any boundary faces in direction a_dir from a_box
  Box noBoundaryBox = a_box;
  noBoundaryBox.enclosedCells(a_dir);
  noBoundaryBox.grow(a_dir,1);
  noBoundaryBox &= m_domain;
  noBoundaryBox.grow(a_dir,-1);
  noBoundaryBox.surroundingNodes(a_dir);

  // WHAT GETS CALCULATED:
  // with d = a_dir,
  // a_F[i+e[d]/2] += a_coeff*min{a_divVel[i+e[d]/2], 0}*(a_U[i+e[d]]-a_U[i]) .
  //
  // That is, if a_divVel[i] >= 0 then no change.
  // Otherwise,
  // a_F[i+e[d]/2] += a_coeff * a_divVel[i+e[d]/2] * (a_U[i+e[d]]-a_U[i]) .
  FORT_ARTVISCF(CHF_FRA(a_F),
                CHF_CONST_FRA(a_U),
                CHF_CONST_FRA1(a_divVel,0),
                CHF_CONST_REAL(a_coeff),
                CHF_CONST_INT(a_dir),
                CHF_BOX(noBoundaryBox));
}

//////////////////////////////////////////////////////////////////////////////
void GodunovUtilities::highOrderLimiter(bool a_highOrderLimiter)
{
  m_highOrderLimiter = a_highOrderLimiter;
}

//////////////////////////////////////////////////////////////////////////////

void GodunovUtilities::deconvolve(FArrayBox&        a_avgFab,
                                  const FArrayBox&  a_cellFab,
                                  int               a_sign)
{
  const Box& bx = a_cellFab.box();
  int numSlopes = a_cellFab.nComp();
  FArrayBox d2Fab(bx, numSlopes);
  Real scale = (a_sign * 1.) /24.;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box loBox, hiBox, centerBox, entireBox;
      int hasLo, hasHi;
      // Generate the domain boundary boxes, loBox and hiBox, if there are
      // domain boundaries there
      loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                 bx, m_domain, idir);

      // d2Fab[i] = a_cellFab[i-e] - 2*a_cellFab[i] + a_cellFab[i+e]
      FORT_GETSECONDDIFF(CHF_FRA(d2Fab),
                         CHF_CONST_FRA(a_cellFab),
                         CHF_CONST_INT(numSlopes),
                         CHF_CONST_INT(idir),
                         CHF_BOX(loBox),
                         CHF_CONST_INT(hasLo),
                         CHF_BOX(hiBox),
                         CHF_CONST_INT(hasHi),
                         CHF_BOX(centerBox));
      // We have started with a_avgFab set to a_cellFab.
      a_avgFab.plus(d2Fab, scale);
    }
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

//////////////////////////////////////////////////////////////////////////////

void GodunovUtilities::deconvolve(FArrayBox&        a_avgFab,
                                  const FArrayBox&  a_cellFab,
                                  const Box&        a_box,
                                  int               a_sign)
{
  int numSlopes = a_cellFab.nComp();
  FArrayBox d2Fab(a_box, numSlopes);
  Real scale = (a_sign * 1.) /24.;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box bx1(a_box);
      bx1.grow(idir, 1);
      Box loBox, hiBox, centerBox, entireBox;
      int hasLo, hasHi;
      // Generate the domain boundary boxes, loBox and hiBox, if there are
      // domain boundaries there
      loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                 bx1, m_domain, idir);

      // d2Fab[i] = a_cellFab[i-e] - 2*a_cellFab[i] + a_cellFab[i+e]
      FORT_GETSECONDDIFF(CHF_FRA(d2Fab),
                         CHF_CONST_FRA(a_cellFab),
                         CHF_CONST_INT(numSlopes),
                         CHF_CONST_INT(idir),
                         CHF_BOX(loBox),
                         CHF_CONST_INT(hasLo),
                         CHF_BOX(hiBox),
                         CHF_CONST_INT(hasHi),
                         CHF_BOX(centerBox));
      // We have started with a_avgFab set to a_cellFab.
      a_avgFab.plus(d2Fab, scale);
    }
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

//////////////////////////////////////////////////////////////////////////////

void GodunovUtilities::deconvolveFace(FluxBox&        a_avgFlux,
                                      const FluxBox&  a_cellFlux,
                                      int             a_sign)
{
  int numSlopes = a_cellFlux.nComp();
  Real scale = (a_sign * 1.) /24.;
  const Box& bxCells = a_cellFlux.box();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox& avgFabDir = a_avgFlux[idir];
      const FArrayBox& cellFabDir = a_cellFlux[idir];
      for (int tanDir = 0; tanDir < SpaceDim; tanDir++)
        {
          if (tanDir != idir)
            {
              Box loBox, hiBox, centerBox, entireBox;
              int hasLo, hasHi;
              // Generate the domain boundary boxes, loBox and hiBox,
              // if there are domain boundaries there.
              // Don't use loHiCenterFace, because that's for 2-point stencils.
              loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                         bxCells, m_domain, tanDir);
              centerBox.surroundingNodes(idir);
              if (hasLo) loBox.surroundingNodes(idir); // lowest layer
              if (hasHi) hiBox.surroundingNodes(idir); // highest layer

              Box bxFaces = bxCells; // from loHiCenterFace
              bxFaces.surroundingNodes(idir);
              FArrayBox d2Fab(bxFaces, numSlopes);
              // for each face i in direction idir:
              // d2Fab[i] = a_cellFab[i-e] - 2*a_cellFab[i] + a_cellFab[i+e]
              FORT_GETSECONDDIFF(CHF_FRA(d2Fab),
                                 CHF_CONST_FRA(cellFabDir),
                                 CHF_CONST_INT(numSlopes),
                                 CHF_CONST_INT(tanDir),
                                 CHF_BOX(loBox),
                                 CHF_CONST_INT(hasLo),
                                 CHF_BOX(hiBox),
                                 CHF_CONST_INT(hasHi),
                                 CHF_BOX(centerBox));
              // We have started with avgFabDir set to cellFabDir.
              avgFabDir.plus(d2Fab, scale);
            }
        }
    }
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

//////////////////////////////////////////////////////////////////////////////

void GodunovUtilities::deconvolveFace(FluxBox&        a_avgFlux,
                                      const FluxBox&  a_cellFlux,
                                      const Box&      a_box,
                                      int             a_sign)
{
  int numSlopes = a_cellFlux.nComp();
  Real scale = (a_sign * 1.) /24.;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox& avgFabDir = a_avgFlux[idir];
      const FArrayBox& cellFabDir = a_cellFlux[idir];
      for (int tanDir = 0; tanDir < SpaceDim; tanDir++)
        {
          if (tanDir != idir)
            {
              Box bx1(a_box);
              bx1.grow(tanDir, 1);
              Box loBox, hiBox, centerBox, entireBox;
              int hasLo, hasHi;
              // Generate the domain boundary boxes, loBox and hiBox,
              // if there are domain boundaries there.
              // Don't use loHiCenterFace, because that's for 2-point stencils.
              loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                         bx1, m_domain, tanDir);
              // We'll set d2Fab on idir-faces (NOT tanDir-faces).
              centerBox.surroundingNodes(idir);
              if (hasLo) loBox.surroundingNodes(idir); // lowest layer
              if (hasHi) hiBox.surroundingNodes(idir); // highest layer

              Box bxFaces = a_box; // from loHiCenterFace
              bxFaces.surroundingNodes(idir);
              FArrayBox d2Fab(bxFaces, numSlopes);
              // for each face i in direction idir:
              // d2Fab[i] = a_cellFab[i-e] - 2*a_cellFab[i] + a_cellFab[i+e]
              FORT_GETSECONDDIFF(CHF_FRA(d2Fab),
                                 CHF_CONST_FRA(cellFabDir),
                                 CHF_CONST_INT(numSlopes),
                                 CHF_CONST_INT(tanDir),
                                 CHF_BOX(loBox),
                                 CHF_CONST_INT(hasLo),
                                 CHF_BOX(hiBox),
                                 CHF_CONST_INT(hasHi),
                                 CHF_BOX(centerBox));
              // We have started with avgFabDir set to cellFabDir.
              avgFabDir.plus(d2Fab, scale);
            }
        }
    }
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

//////////////////////////////////////////////////////////////////////////////
// Compute a face-centered nonlinear function of the divergence suitable for use
// as an artificial viscosity coefficient for a fourth-order method.
// F_visc = -dx*K*dU/dx
// K = max(-dx*div(u)*min(1,h^2*|div(u)/(c^2)/a_M0sq|),0)

void GodunovUtilities::divVelHO(FArrayBox&       a_divVel,
                                const FArrayBox& a_W,
                                const int&       a_dir,
                                const Box&       a_box,
                                GodunovPhysics* a_physPtr)
{
  CH_assert((a_dir >= 0 )&& (a_dir < SpaceDim));
  CH_assert(a_divVel.box().contains(a_box));
  CH_assert(a_physPtr->fourthOrderArtificialViscosityIsDefined() == true);
  // Compute face-centered divergence.
  Interval velInt = a_physPtr->velocityInterval();
  // a_W on valid cells + 1 layer of ghosts;
  // a_box = a_dir-faces of valid cells
  divVel(a_divVel, a_W, velInt, a_dir, a_box);

  // box1cells = valid cells + 1 layer of ghosts;
  // boxcsq = box1cells intersected with domain.
  Box box1cells = a_box;
  box1cells.enclosedCells();
  int ghostbox1 = 1;
  box1cells.grow(ghostbox1);
  Box boxcsq = box1cells & m_domain;

  // Compute cell-centered (bulk modulus)/(density).
  int bulkIndex = a_physPtr->bulkModulusIndex();
  int densityIndex = a_physPtr->densityIndex();
  Real M0sq = a_physPtr->getFourthOrderArtificialViscosityParameter();
  FArrayBox csq1(boxcsq, 1);
  FArrayBox csq2(boxcsq, 1);
  // csq1 = a_W[bulkIndex] / a_W[densityIndex] on boxcsq,
  // which is valid cells + 1 ghost layer, intersected with domain.
  csq1.setVal(0.);
  csq2.setVal(0.);
  csq1.copy(a_W, boxcsq, bulkIndex, boxcsq, 0, 1);
  csq1.divide(a_W, boxcsq, boxcsq, densityIndex, 0, 1);
  Box hiBox,loBox,centerBox,entireBox;
  int hasLo,hasHi;
  FArrayBox* csqin_ptr = &csq1;
  FArrayBox* csqout_ptr = &csq2;

  // Compute min of csq in transverse direction(s).
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (dir != a_dir)
        {
          // WHAT GETS CALCULATED:
          // Pass on transverse direction dir
          // takes *csqin_ptr and stores result in *csqout_ptr.
          // csqout[i] =
          // min{csqin[i-e[dir]], csqin[i], csqin[i+e[dir]]}  in center;
          // min{csqin[i-e[dir]], csqin[i]}                   at high end;
          //                  min{csqin[i], csqin[i+e[dir]]}  at low end.
          FArrayBox& csqin = *csqin_ptr;
          FArrayBox& csqout = *csqout_ptr;
          loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                     box1cells, m_domain,dir);

          FORT_MIN3PTSF(CHF_FRA1(csqout,0),
                        CHF_CONST_FRA1(csqin,0),
                        CHF_CONST_INT(dir),
                        CHF_BOX(loBox),
                        CHF_CONST_INT(hasLo),
                        CHF_BOX(hiBox),
                        CHF_CONST_INT(hasHi),
                        CHF_BOX(centerBox));
          box1cells.grow(dir,-1);
          FArrayBox* csqtmp = csqin_ptr;
          csqin_ptr = csqout_ptr;
          csqout_ptr = csqtmp;
        }
    }
  // divBox = valid cells + 1 ghost layer in direction a_dir
  Box divBox = a_box;
  divBox.enclosedCells(a_dir);
  divBox.grow(a_dir, 1);
  loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
             divBox, m_domain, a_dir);

  // All of the boxes computed above are shifted so as to be cell-centered,
  // with the index of the cell center being identified with the low face.
  // We then shift a_divVel to be compatible with that convention on input,
  // and undo the shift on output.  Basically, everything is made compatible
  // with a_W (which is cell-centered).

  // WHAT GETS CALCULATED:
  // Here we multiply face-centered a_divVel by min(a_divVel**2/csqmin/M0sq, 1)
  // where csqmin is min of csq on the two neighboring cells
  // (or csq on the one neighboring cell if on a boundary).
  // That is, if a_divVel**2 >= csqmin*M0sq, then no change in a_divVel.
  // Otherwise, a_divVel *= a_divVel**2/csqmin/M0sq.

  // shift a_divVel from a_dir-faces of valid cells
  // to valid cells + 1 ghost layer in Hi direction a_dir
  a_divVel.shiftHalf(a_dir, 1);
  // if (!hasLo && !hasHi)
  //   centerBox = valid cells - 1 layer in Lo direction a_dir
  // centerBox.growDir(a_dir, Side::Lo, -1);
  // petermc, 10 Jul 2009, changed to this
  centerBox.growDir(a_dir, Side::Hi, +1);
  // csq on valid cells + 1 ghost layer, intersected with domain.
  // csq contains pressure / density, after passes along transverse
  // directions taking minimum with neighbors.
  FArrayBox& csq = *csqin_ptr;
  FORT_HODIVCOEF(CHF_FRA1(a_divVel,0), // set this on centerBox, loBox, hiBox
                 CHF_CONST_FRA1(csq,0),
                 CHF_CONST_INT(a_dir),
                 CHF_CONST_REAL(M0sq),
                 CHF_BOX(loBox),
                 CHF_CONST_INT(hasLo),
                 CHF_BOX(hiBox),
                 CHF_CONST_INT(hasHi),
                 CHF_BOX(centerBox));
  a_divVel.shiftHalf(a_dir,-1);

  // if (!hasLo && !hasHi)
  //   a_divVel has been set on a_dir-faces BETWEEN valid cells
}


#include "NamespaceFooter.H"
