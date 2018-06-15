#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBPatchAdvect.H"
#include "EBPatchAdvectF_F.H"
#include "EBLoHiCenter.H"
#include "EBArith.H"
#include "PolyGeom.H"
#include "FaceIterator.H"
#include "VoFIterator.H"
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>
#include "CH_Timer.H"
#include "EBPatchGodunovF_F.H"
#include "NamespaceHeader.H"
/******/
void
EBPatchAdvect::
setValidBox(const Box&        a_validBox,
            const EBISBox&    a_ebisBox,
            const IntVectSet& a_coarseFineIVS,
            const Real&       a_cumulativeTime,
            const Real&       a_timeStep)
{
  CH_TIME("EBPatchAdvect::setValidBox");
  EBPatchGodunov::setValidBox(a_validBox, a_ebisBox, a_coarseFineIVS, a_cumulativeTime, a_timeStep);

}

void
EBPatchAdvect::
slope(EBCellFAB&       a_slopePrim,
      EBCellFAB&       a_slopeNLim,
      const EBCellFAB& a_primState,
      const EBCellFAB& a_flattening,
      const int&       a_dir,
      const Box&       a_box,
      bool a_doAgg)
{
  CH_TIME("EBPatchAdvect::slope");
  if (a_doAgg)
    {
      EBPatchGodunov::slope(a_slopePrim, a_slopeNLim, a_primState, a_flattening, a_dir, a_box, a_doAgg);
    }
  else
    {
      CH_assert(a_slopePrim.nComp() == numSlopes());
      CH_assert(a_primState.nComp() >= numSlopes());
      EBCellFAB deltaWL, deltaWR;

      {
        doSecondOrderSlopes(a_slopePrim,
                            deltaWL,
                            deltaWR,
                            a_slopeNLim,
                            a_primState,
                            a_dir,
                            a_box);
      }
    }
}
void
EBPatchAdvect::
doSecondOrderSlopes(EBCellFAB&       a_delta2W,
                    EBCellFAB&       a_deltaWL,
                    EBCellFAB&       a_deltaWR,
                    EBCellFAB&       a_deltaWC,
                    const EBCellFAB& a_primState,
                    const int&       a_dir,
                    const Box&       a_box)
{
  CH_TIME("EBPatchAdvect::doSecondOrderSlopes");
  CH_assert(m_isSlopeSet);
  int numSlope = numSlopes();
  Box box1 = a_box;
  box1.grow(a_dir, 1);
  box1 &= m_domain;

  Box box2 = a_box;
  box2.grow(a_dir, 2);
  box2 &= m_domain;

  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;

  if (usesFourthOrderSlopes())
    {
      eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                   box2, m_domain, a_dir);
    }
  else
    {
      eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                   box1, m_domain, a_dir);
    }

  {
    CH_TIME("defines_and_setvals");
    a_deltaWL.define(m_ebisBox, entireBox, numSlope);
    a_deltaWR.define(m_ebisBox, entireBox, numSlope);
    a_delta2W.setVal(0.);
    a_deltaWL.setVal(0.);
    a_deltaWR.setVal(0.);
    a_deltaWC.setVal(0.);
  }


  {
    CH_TIME("fortran");
    BaseFab<Real>& regDelta2W = a_delta2W.getSingleValuedFAB();
    BaseFab<Real>& regDeltaWC = a_deltaWC.getSingleValuedFAB();
    BaseFab<Real>& regDeltaWL = a_deltaWL.getSingleValuedFAB();
    BaseFab<Real>& regDeltaWR = a_deltaWR.getSingleValuedFAB();
    const BaseFab<Real>& regPrimState = a_primState.getSingleValuedFAB();

    FORT_SECONDSLOPEDIFFS(CHF_FRA(regDeltaWC),
                          CHF_FRA(regDeltaWL),
                          CHF_FRA(regDeltaWR),
                          CHF_CONST_FRA(regPrimState),
                          CHF_CONST_INT(numSlope),
                          CHF_CONST_INT(a_dir),
                          CHF_BOX(loBox),
                          CHF_CONST_INT(hasLo),
                          CHF_BOX(hiBox),
                          CHF_CONST_INT(hasHi),
                          CHF_BOX(centerBox));

    //apply van leer limiter to regular cells.
    //limiting for irregular cells is more complicated
    //now that we are doing higher order limited
    //one-sided diffs
    regDelta2W.copy(regDeltaWC);
    if (m_useLimiting)
      {
        FORT_VLLIMITER(CHF_FRA(regDelta2W),
                       CHF_CONST_FRA(regDeltaWL),
                       CHF_CONST_FRA(regDeltaWR),
                       CHF_BOX(centerBox));
      }
  }

  {
    CH_TIME("irreg_slopes");
    //at irregular cells, if one side does not exist,
    //set all to one-sided diffs.  try to do higher order
    //limited diffs if possible.
    for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        if (entireBox.contains(vof.gridIndex()))
          {
            int ivar = 0;
            bool hasFacesLeft, hasFacesRigh;
            Real dql, dqr, dqc;
            bool verbose =false;
            pointGetSlopes(dql, dqr, dqc,
                           hasFacesLeft,
                           hasFacesRigh,
                           vof, a_primState, a_dir, ivar, verbose);
            Real dqSec = 0;
            if (!m_useLimiting)
              {
                dqSec = dqc;
              }
            else
              {
                if (hasFacesLeft && hasFacesRigh)
                  {
                    Real dqlim = dqc;

                    FORT_POINTVLLIMITER(CHF_REAL(dqlim),
                                        CHF_CONST_REAL(dql),
                                        CHF_CONST_REAL(dqr));
                    dqSec = dqlim;
                  }
                else if (!hasFacesLeft && !hasFacesRigh)
                  {
                    dqSec = 0.0;
                  }
                else if (hasFacesLeft && !hasFacesRigh)
                  {
                    if (dqc*dql > TOL)
                      {
                        Real rsign = 1.0;
                        if (dqc < -TOL)
                          {
                            rsign = -1.0;
                          }
                        dqSec = rsign*Min(Abs(dqc), Abs(dql));
                      }
                    else
                      {
                        dqSec = 0.0;
                      }

                  }
                else if (hasFacesRigh && !hasFacesLeft)
                  {
                    if (dqc*dqr > TOL)
                      {
                        Real rsign = 1.0;
                        if (dqc < -TOL)
                          {
                            rsign = -1.0;
                          }
                        dqSec = rsign*Min(Abs(dqc), Abs(dqr));
                      }
                    else
                      {
                        dqSec = 0.0;
                      }
                  }
                else
                  {
                    MayDay::Error("doh! missed a case");
                  }
              }//end if using limiting

            if (usesFlattening())
              {

                int pindex = pressureIndex();
                Real press = Max(a_primState(vof, pindex), Real(1.0e-10));
                Real dqp, dqlp, dqrp;
                verbose = false;
                pointGetSlopes(dqlp, dqrp, dqp,
                               hasFacesLeft,
                               hasFacesRigh,
                               vof, a_primState, a_dir, pindex,verbose);
                //try to detect pressure discontinutity
                dqp = Max(Abs(dqp), Abs(dqlp));
                dqp = Max(Abs(dqp), Abs(dqrp));
                if (Abs(dqp/press) > 0.1)
                  {
                    dqSec = 0.0;
                  }
              }
            a_delta2W(vof, ivar) = dqSec;
            a_deltaWL(vof, ivar) = dql;
            a_deltaWR(vof, ivar) = dqr;
            a_deltaWC(vof, ivar) = dqc;
          } //if vof is in this box
      }//end loop over vofs
  }
  {
    CH_TIME("bndry_slopes");
    //want to be able to call this in test codes
    if (m_isBCSet)
      {
        m_bc->setBndrySlopes(a_delta2W, a_primState, m_ebisBox, entireBox, a_dir);
      }
  }
}


/******************/
EBPatchAdvectFactory::
EBPatchAdvectFactory(RefCountedPtr<EBPhysIBCFactory>& a_bcFactoryPtr,
                     const bool&                      a_useLimiting,
                     const Real&                      a_maxVal,
                     const Real&                      a_minVal,
                     const bool&                      a_setMaxMin)
  :EBPatchGodunovFactory()
{
  m_bcFactoryPtr = a_bcFactoryPtr;
  m_useLimiting  = a_useLimiting;
  m_setMaxMin    = a_setMaxMin;
  m_maxVal       = a_maxVal;
  m_minVal       = a_minVal;
}
/******************/
EBPatchAdvectFactory::
~EBPatchAdvectFactory()
{
}
/******************/
EBPatchGodunov*
EBPatchAdvectFactory::
create() const
{
  EBPatchAdvect* retval= new EBPatchAdvect();
  retval->setEBPhysIBC(*m_bcFactoryPtr);
  retval->useLimiting(m_useLimiting);
  if (m_setMaxMin)
    {
      retval->setMaxMin(m_maxVal, m_minVal);
    }
  return static_cast<EBPatchGodunov*>(retval);
}
/******************/

/******/
void
EBPatchAdvect::
extrapToCoveredFaces(BaseIVFAB<Real>&        a_extendedPrim,
                     const EBFaceFAB&        a_primFace,
                     const EBCellFAB&        a_primState,
                     const Vector<VolIndex>& a_coveredFaces,
                     const int&              a_faceDir,
                     const Side::LoHiSide&   a_sd,
                     const Box&              a_box)
{
  CH_TIME("EBPatchAdvect::extrapToCoveredFaces");
  for (int ivof = 0; ivof < a_coveredFaces.size(); ivof++)
    {
      const VolIndex& vof = a_coveredFaces[ivof];

      RealVect normal = m_ebisBox.normal(vof);
      if (a_box.contains(vof.gridIndex()))
        {
          const int numPrim = numPrimitives();
          Vector<Real> extPrim(numPrim, 0.0);
#if CH_SPACEDIM==2
          pointExtrapToCovered2D(extPrim,
                                 a_primFace,
                                 a_primState,
                                 a_faceDir,
                                 vof,
                                 normal,
                                 a_sd,
                                 numPrim);
#elif CH_SPACEDIM==3
          pointExtrapToCovered3D(extPrim,
                                 a_primFace,
                                 a_primState,
                                 a_faceDir,
                                 vof,
                                 normal,
                                 a_sd,
                                 numPrim);
#else
          bogus_ch_spacedim();
#endif
          for (int ivar = 0; ivar < numPrim; ivar++)
            {
              a_extendedPrim(vof, ivar) = extPrim[ivar];
            }
        }
    }
  const IntVectSet& ivs  = a_extendedPrim.getIVS();
  floorPrimitives(a_extendedPrim, ivs);
}

/********/
void
EBPatchAdvect::
pointExtrapToCovered2D(Vector<Real>&           a_extrapVal,
                       const EBFaceFAB&        a_primFace,
                       const EBCellFAB&        a_primState,
                       const int&              a_faceDir,
                       const VolIndex&         a_vof,
                       const RealVect&         a_normal,
                       const Side::LoHiSide&   a_sd,
                       const int&              a_numPrim)
{
  CH_assert(SpaceDim==2);
  int tangenDir = 1 - a_faceDir;

  int signNorm = 1;
  int signTang = 1;
  if (a_normal[a_faceDir] < 0.0) signNorm = -1;
  if (Abs(a_normal[tangenDir]) == 0.0)
    {//since signTang is arbitrary for this case,
      //check if we are next to the high domain edge
      IntVect ivHi = a_vof.gridIndex();;
      ivHi[tangenDir] += 1;
      if (!m_domain.contains(ivHi))
        {
          signTang = -1;
        }
    }
  else if (a_normal[tangenDir] < 0.0)
    {
      signTang = -1;
    }

  //iv[0][0] is the vof. iv[1][0] is vofside.  iv[0][1] is vofup iv[1][1] is vofcorner
  IntVect   ivSten[2][2];
  bool      hasVoF[2][2];
  VolIndex     vof[2][2];
  FaceIndex   face[2][2];
  Real         val[2][2];
  ivSten[0][0] = a_vof.gridIndex();
  ivSten[1][0] = ivSten[0][0] + signNorm*BASISV(a_faceDir);
  ivSten[0][1] = ivSten[0][0] - signNorm*BASISV(a_faceDir) + signTang*BASISV(tangenDir) ;;
  ivSten[1][1] = ivSten[0][0]                              + signTang*BASISV(tangenDir) ;

  int radius = 1;
  Vector<VolIndex> vofsStencil;
  EBArith::getAllVoFsInMonotonePath(vofsStencil, a_vof, m_ebisBox, radius);

  for (int ivar = 0; ivar < a_numPrim; ivar++)
    {
      for (int ix = 0; ix < 2; ix++)
        {
          for (int iy = 0; iy < 2; iy++)
            {
              Vector<FaceIndex> facesTmp;
              hasVoF[ix][iy] = EBArith::isVoFHere(vof[ix][iy], vofsStencil, ivSten[ix][iy]);
              if (hasVoF[ix][iy])
                {
                  facesTmp = m_ebisBox.getFaces(vof[ix][iy], a_faceDir, flip(a_sd));
                  hasVoF[ix][iy] = (hasVoF[ix][iy] && (facesTmp.size() == 1));
                }
              if (hasVoF[ix][iy])
                {
                  face[ix][iy] = facesTmp[0];
                  val [ix][iy] = a_primFace(face[ix][iy], ivar);
                }
              else
                {
                  val[ix][iy] = 0.0;
                }
            }
        }

      int dirBigNorm, dirLitNorm;
      VolIndex vofBigNorm, vofLitNorm;
      bool hasVoFBigNorm,  hasVoFLitNorm;
      int  signBigNorm, signLitNorm;
      Real wBigNorm, wLitNorm;
      Real dWBigNorm[2],dWLitNorm[2];
      if (Abs(a_normal[tangenDir]/a_normal[a_faceDir]) > m_dx[tangenDir]/m_dx[a_faceDir])
        {
          dirBigNorm = tangenDir;
          dirLitNorm = a_faceDir;

          signBigNorm   = signTang;
          signLitNorm   = signNorm;
          hasVoFBigNorm = hasVoF[0][1];
          if (hasVoF[0][1])
            {
              wBigNorm      =    val[0][1];
              vofBigNorm    =    vof[0][1];
            }
          else
            {
              wBigNorm      =    val[0][0];
              vofBigNorm    =    vof[0][0];
            }
        }
      else
        {
          dirBigNorm = a_faceDir;
          dirLitNorm = tangenDir;

          signBigNorm   = signNorm;
          signLitNorm   = signTang;
          hasVoFBigNorm = hasVoF[1][0];

          if (hasVoF[1][0])
            {
              //              Real deltaW = a_slopesSeco[a_faceDir](vof[1][0]  , ivar);
              Real deltaW;
              coveredExtrapSlopes(deltaW, vof[1][0],a_primState, a_faceDir, ivar);
              wBigNorm      =    val[1][0] - signBigNorm*deltaW;
              vofBigNorm    =    vof[1][0];
            }
          else
            {
              wBigNorm    =    val[0][0];
              vofBigNorm  =    vof[0][0];
            }

        }
      hasVoFLitNorm = hasVoF[1][1];
      if (hasVoF[1][1])
        {
          wLitNorm      =    val[1][1];
          vofLitNorm    =    vof[1][1];
        }
      else
        {
          wLitNorm      =    val[0][0];
          vofLitNorm    =    vof[0][0];
        }

      if (hasVoFLitNorm && hasVoFBigNorm)
        {
          coveredExtrapSlopes(dWBigNorm[0],vofBigNorm,a_primState,0,ivar);
          coveredExtrapSlopes(dWBigNorm[1],vofBigNorm,a_primState,1,ivar);
          coveredExtrapSlopes(dWLitNorm[0],vofLitNorm,a_primState,0,ivar);
          coveredExtrapSlopes(dWLitNorm[1],vofLitNorm,a_primState,1,ivar);

          const Real xLit = Abs( (m_dx[dirBigNorm]*a_normal[dirLitNorm])/ (m_dx[dirLitNorm]*a_normal[dirBigNorm]) );
          const Real xBig = 1.0 - xLit;

          const Real dWBigExtrap = xLit*dWLitNorm[dirBigNorm] + xBig*dWBigNorm[dirBigNorm];
          const Real dWLitExtrap = xLit*dWLitNorm[dirLitNorm] + xBig*dWBigNorm[dirLitNorm];

          const Real wc = xLit*wLitNorm + xBig*wBigNorm;
          a_extrapVal[ivar] = wc -signBigNorm*dWBigExtrap - xLit*signLitNorm*dWLitExtrap;
          // a_extrapVal[ivar] = wc;
        }
      else if (hasVoF[1][0])
        {//change methods, still second order
          Real deltaW;
          coveredExtrapSlopes(deltaW,vof[1][0],a_primState,a_faceDir,ivar);

          a_extrapVal[ivar] = val[1][0] - 2.0*signNorm*deltaW;
          // a_extrapVal[ivar] = val[1][0];
        }
      else if (hasVoFLitNorm)
        {//drop order
          //pout() << "EBPatchAdvect::pointExtrapToCovered2D:lit: vof " << a_vof << " drops order " << endl;
          a_extrapVal[ivar] = wLitNorm;
        }
      else if (hasVoFBigNorm)
        {
          //pout() << "EBPatchAdvect::pointExtrapToCovered2D:big: vof " << a_vof << " drops order " << endl;
          a_extrapVal[ivar] = wBigNorm;
        }
      else
        {
          //pout() << "EBPatchAdvect::pointExtrapToCovered2D:none: vof " << a_vof << " drops order " << endl;
          a_extrapVal[ivar] = a_primState(a_vof,ivar);
        }
    }
}

void
EBPatchAdvect::
pointExtrapToCovered3D(Vector<Real>&           a_extrapVal,
                       const EBFaceFAB&        a_primFace,
                       const EBCellFAB&        a_primState,
                       const int&              a_faceDir,
                       const VolIndex&         a_vof,
                       const RealVect&         a_normal,
                       const Side::LoHiSide&   a_sd,
                       const int&              a_numPrim)
{
  RealVect normal = a_normal;
  //HACK put in true normal
  //RealVect truenorm = getTrueNorm(a_vof,a_faceDir,a_sd);
  //normal = truenorm;
  //END HACK

  CH_assert(SpaceDim==3);
  Tuple<int,CH_SPACEDIM-1> tangenDir = PolyGeom::computeTanDirs(a_faceDir);
  Real anormNorm = Abs(normal[a_faceDir]);
  Real anormTan[2];
  int signNorm = 1;
  if (normal[a_faceDir] < 0.0) signNorm = -1;
  int signTang[2];
  for (int itan = 0; itan < 2; itan++)
    {
      int tanDir = tangenDir[itan];
      anormTan[itan] = Abs(normal[tanDir]);
      signTang[itan] =  1;
      if (anormTan[itan] == 0.0)
        {//since signTang is arbitrary for this case,
          //check if we are next to the high domain edge
          IntVect ivHi = a_vof.gridIndex();;
          ivHi[tanDir] += 1;
          if (!m_domain.contains(ivHi))
            {
              signTang[itan] = -1;
            }
        }
      else if (normal[tanDir] < 0.0)
        {
          signTang[itan] = -1;
        }
    }

  const IntVect& ivVoF= a_vof.gridIndex();

  //whice one of the tangential directions has the largest normal
  int d1, d2;
  if (anormTan[0]/m_dx[tangenDir[0]] > anormTan[1]/m_dx[tangenDir[1]])
    {
      d1 = 0;
      d2 = 1;
    }
  else
    {
      d1 = 1;
      d2 = 0;
    }

  // figure out in which plane we are extrapolating
  bool faceDirOut = ((anormNorm/m_dx[a_faceDir] > anormTan[0]/m_dx[tangenDir[0]]) &&
                     (anormNorm/m_dx[a_faceDir] > anormTan[1]/m_dx[tangenDir[1]]));

  IntVect ivSten[2][2];
  bool hasVoF[2][2];
  VolIndex    vofSten[2][2];
  FaceIndex  faceSten[2][2];

  if (faceDirOut)
    {// face direction has largest normal.

      ivSten[0][0] = ivVoF + signNorm*BASISV(a_faceDir);
      ivSten[1][0] = ivVoF + signTang[0]*BASISV(tangenDir[0]);
      ivSten[0][1] = ivVoF + signTang[1]*BASISV(tangenDir[1]);
      ivSten[1][1] = ivVoF + signTang[0]*BASISV(tangenDir[0]) + signTang[1]*BASISV(tangenDir[1]) ;
      d1 = 0;
      d2 = 1;

    }
  else
    { //tandir[d1] is the biggest

      ivSten[0][0] = ivVoF + signTang[d1]*BASISV(tangenDir[d1])                                      - signNorm*BASISV(a_faceDir);
      ivSten[1][0] = ivVoF + signTang[d1]*BASISV(tangenDir[d1])                                      ;
      ivSten[0][1] = ivVoF + signTang[d1]*BASISV(tangenDir[d1]) + signTang[d2]*BASISV(tangenDir[d2]) - signNorm*BASISV(a_faceDir);
      ivSten[1][1] = ivVoF + signTang[d1]*BASISV(tangenDir[d1]) + signTang[d2]*BASISV(tangenDir[d2]) ;

    }

  int radius = 1;
  Vector<VolIndex> vofsStencil;
  EBArith::getAllVoFsInMonotonePath(vofsStencil, a_vof, m_ebisBox, radius);

  //get the face info for a 1D extrap stencil (if needed)
  VolIndex vofSten1D;
  FaceIndex faceSten1D;
  const IntVect ivSten1D = ivVoF + signNorm*BASISV(a_faceDir);
  bool has1DSten = EBArith::isVoFHere(vofSten1D, vofsStencil, ivSten1D);
  Vector<FaceIndex> facesTmp1D;
  if (has1DSten)
    {
      facesTmp1D = m_ebisBox.getFaces(vofSten1D, a_faceDir, flip(a_sd));
      if (facesTmp1D.size() == 1)
        {
          faceSten1D = facesTmp1D[0];
        }
      else
        {
          has1DSten = false;
        }
    }

  for (int ix = 0; ix < 2; ix++)
    {
      for (int iy = 0; iy < 2; iy++)
        {
          hasVoF[ix][iy] =  EBArith::isVoFHere(vofSten[ix][iy], vofsStencil, ivSten[ix][iy]);
          Vector<FaceIndex> facesTmp;
          if (hasVoF[ix][iy])
            {
             facesTmp = m_ebisBox.getFaces(vofSten[ix][iy], a_faceDir, flip(a_sd));
             hasVoF[ix][iy] = (hasVoF[ix][iy] && (facesTmp.size() == 1));
            }
          if (hasVoF[ix][iy])
            {
              faceSten[ix][iy] = facesTmp[0];
            }
        }
    }
  bool hasAllVoFs = hasVoF[0][0] && hasVoF[1][0] && hasVoF[0][1] && hasVoF[1][1];

  for (int ivar = 0; ivar < a_numPrim; ivar++)
    {
      if (hasAllVoFs)
        {
          Real WVal[2][2];
          for (int ix = 0; ix < 2; ix++)
            {
              for (int iy = 0; iy < 2; iy++)
                {
                  WVal[ix][iy] = a_primFace(faceSten[ix][iy], ivar);
                }
            }
          Real delWNorm[2][2];
          Real delWTan0[2][2];
          Real delWTan1[2][2];
          for (int ix = 0; ix < 2; ix++)
            {
              for (int iy = 0; iy < 2 ; iy++)
                {
                  coveredExtrapSlopes(delWNorm[ix][iy], vofSten[ix][iy], a_primState, a_faceDir   , ivar);
                  coveredExtrapSlopes(delWTan0[ix][iy], vofSten[ix][iy], a_primState, tangenDir[d1], ivar);
                  coveredExtrapSlopes(delWTan1[ix][iy], vofSten[ix][iy], a_primState, tangenDir[d2], ivar);
                }
            }
          if (faceDirOut)
            {
              // case 00 is special because we have to move it backwards in face direction to get it to plane
              Real deltaW;
              coveredExtrapSlopes(deltaW, vofSten[0][0],a_primState, a_faceDir, ivar);
              //Real slipSlope = a_slopesSeco[a_faceDir](vofSten[0][0], ivar);
              Real slipSlope = deltaW;
              WVal[0][0] -= signNorm*slipSlope;
              //debug to make compatible with 2d
              //WVal[0][1] = WVal[0][0];
              //end debug

              Real x0Val = (m_dx[a_faceDir]*anormTan[0])/(m_dx[tangenDir[0]]*anormNorm);
              Real x1Val = (m_dx[a_faceDir]*anormTan[1])/(m_dx[tangenDir[1]]*anormNorm);
              Real    funcVal = bilinearFunc(    WVal,x0Val,x1Val);
              Real deltaWNorm = bilinearFunc(delWNorm,x0Val,x1Val);
              Real deltaWTan0 = bilinearFunc(delWTan0,x0Val,x1Val);
              Real deltaWTan1 = bilinearFunc(delWTan1,x0Val,x1Val);
              Real dWdn =  -signNorm*deltaWNorm - signTang[d2]*x1Val*deltaWTan1 - signTang[d1]*x0Val*deltaWTan0;
              a_extrapVal[ivar] = funcVal + dWdn;

            }
          else
            {
              //tanDir[d1] the biggest normal
              Real x0Val = (m_dx[tangenDir[d1]]*anormNorm   )/(m_dx[    a_faceDir]*anormTan[d1]);
              Real x1Val = (m_dx[tangenDir[d1]]*anormTan[d2])/(m_dx[tangenDir[d2]]*anormTan[d1]);
              Real    funcVal = bilinearFunc(WVal,    x0Val,x1Val);
              Real deltaWNorm = bilinearFunc(delWNorm,x0Val,x1Val);
              Real deltaWTan0 = bilinearFunc(delWTan0,x0Val,x1Val);
              Real deltaWTan1 = bilinearFunc(delWTan1,x0Val,x1Val);

              Real dWdn =  -signNorm*x0Val*deltaWNorm - signTang[d2]*x1Val*deltaWTan1 - signTang[d1]*deltaWTan0;
              a_extrapVal[ivar] = funcVal + dWdn;

            }
        }
      else if (has1DSten)
        {//try 1D extrap
          const Real WVal1D = a_primFace(faceSten1D, ivar);
          Real deltaW;
          coveredExtrapSlopes(deltaW,vofSten1D,a_primState,a_faceDir,ivar);
          a_extrapVal[ivar] = WVal1D - 2.0*signNorm*deltaW;
        }
      else
        // At least one of the vofs in the stencil exists - we drop order
        // by using a weighted sum of the the wVal's.
        {
          //pout() << "EBPatchAdvect::pointExtrapToCovered3D: vof " << a_vof << " drops order " << endl;
          a_extrapVal[ivar] = 0.;
          Real volTot = 0.;
          bool hasAVoF = false;
          for (int ix = 0; ix < 2; ix++)
            {
              for (int iy = 0; iy < 2; iy++)
                {
                  if (hasVoF[ix][iy])
                    {
                      hasAVoF = true;
                      Real wVal = a_primFace(faceSten[ix][iy], ivar);
                      Real volFrac = m_ebisBox.volFrac(vofSten[ix][iy]);
                      a_extrapVal[ivar] += wVal*volFrac;
                      volTot += volFrac;
                    }
                }
            }
          if (hasAVoF)
            {
              CH_assert(volTot > 0.0);
              a_extrapVal[ivar] = a_extrapVal[ivar]/volTot;
            }
          else
            {
              a_extrapVal[ivar] = a_primState(a_vof, ivar);
            }
        }
    }
}
void
EBPatchAdvect::
interpolateFluxToCentroids(BaseIFFAB<Real>      a_centroidFlux[SpaceDim],
                           const EBFluxFAB&     a_fluxInterpolant,
                           const IntVectSet&    a_irregIVS)
{
  CH_TIME("EBPatchAdvect::interpolateFluxToCentroids");
  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;

  //could just use numFluxes but this allows unit testing
  //with a single variable interpolant
  int nflux = a_centroidFlux[0].nComp();
  //now loop through the irregular faces
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      BaseIFFAB<Real>& fluxDir = a_centroidFlux[faceDir];
      const EBFaceFAB& interpol = a_fluxInterpolant[faceDir];
      const BaseIFFAB<FaceStencil>& stencils = m_interpStencils[faceDir];

      for (FaceIterator faceit(a_irregIVS, m_ebisBox.getEBGraph(), faceDir, stopCrit);
          faceit.ok(); ++faceit)
        {
          const FaceIndex& centFace = faceit();
          //          Real     areaFrac = m_ebisBox.areaFrac(centFace);
          //          RealVect centroid = m_ebisBox.centroid(centFace);
          const FaceStencil& sten = stencils(centFace, 0);
          for (int ivar = 0; ivar < nflux; ivar++)
            {
              //              Real centVal = interpol(centFace, ivar);
              Real newFlux = 0.0;
              for (int isten = 0; isten < sten.size(); isten++)
                {
                  const FaceIndex& stenFace= sten.face(isten);
                  Real weight = sten.weight(isten);
                  Real interpFlux = interpol(stenFace, ivar);
                  newFlux += weight*interpFlux;
                }
              fluxDir(centFace, ivar) = newFlux;
            }

          for (int ivar = 0; ivar < numFluxes(); ivar++)
            {
              //debug!  turn off interpolation
              //fluxDir(centFace, ivar) = interpol(centFace, ivar);
            }
        }
    }
}

/*****************************/
void EBPatchAdvect::
normalPred(EBCellFAB&       a_rhoLo,
           EBCellFAB&       a_rhoHi,
           const EBCellFAB& a_rho,
           const EBCellFAB& a_dRho,
           const Real&      a_dtbydx,
           const int&       a_dir,
           const Box&       a_box)
{
  CH_TIME("EBPatchAdvect::normalPred");
  CH_assert(m_isVelSet);
  const EBCellFAB& velcc = *m_normalVelPtr;


  int ivar = 0;

  /**/
  BaseFab<Real>& regRhoLo    = a_rhoLo.getSingleValuedFAB();
  BaseFab<Real>& regRhoHi    = a_rhoHi.getSingleValuedFAB();
  const BaseFab<Real>& regVel =  velcc.getSingleValuedFAB();
  const BaseFab<Real>& regRho =  a_rho.getSingleValuedFAB();
  const BaseFab<Real>& regDRho= a_dRho.getSingleValuedFAB();

  FORT_PREDADVECT(CHF_BOX(a_box),
                  CHF_CONST_FRA1(regRho,ivar),
                  CHF_CONST_FRA1(regDRho, ivar),
                  CHF_CONST_FRA(regVel),
                  CHF_FRA1(regRhoLo, ivar),
                  CHF_FRA1(regRhoHi, ivar),
                  CHF_CONST_INT(a_dir),
                  CHF_CONST_REAL(a_dtbydx));
  /**/

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      Real dense, denlo, denhi, denslope;
      RealVect veloc;

      dense    = a_rho(vof, ivar);
      denslope = a_dRho(vof, ivar);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          veloc[idir] = velcc(vof, idir);
        }

      FORT_POINTPREDADVECT(CHF_REAL(dense),
                           CHF_REAL(denlo),
                           CHF_REAL(denhi),
                           CHF_REAL(denslope),
                           CHF_REALVECT(veloc),
                           CHF_CONST_INT(a_dir),
                           CHF_CONST_REAL(a_dtbydx));

      a_rhoLo(vof, ivar) = denlo;
      a_rhoHi(vof, ivar) = denhi;
    }
}
/*****************************/
void EBPatchAdvect::
transversePred(EBCellFAB&       a_rhoLo,
               EBCellFAB&       a_rhoHi,
               const EBCellFAB& a_rho,
               const EBCellFAB& a_dRho,
               const Real&      a_dtbydx,
               const int&       a_dir,
               const Box&       a_box)
{
  CH_TIME("EBPatchAdvect::normalPred");
  CH_assert(m_isVelSet);
  const EBCellFAB& velcc = *m_normalVelPtr;

  int ivar = 0;

  /**/
  BaseFab<Real>& regRhoLo    = a_rhoLo.getSingleValuedFAB();
  BaseFab<Real>& regRhoHi    = a_rhoHi.getSingleValuedFAB();
  const BaseFab<Real>& regVel =  velcc.getSingleValuedFAB();
  const BaseFab<Real>& regRho =  a_rho.getSingleValuedFAB();
  const BaseFab<Real>& regDRho= a_dRho.getSingleValuedFAB();

  FORT_PREDADVECTTRANS(CHF_BOX(a_box),
                       CHF_CONST_FRA1(regRho,ivar),
                       CHF_CONST_FRA1(regDRho, ivar),
                       CHF_CONST_FRA(regVel),
                       CHF_FRA1(regRhoLo, ivar),
                       CHF_FRA1(regRhoHi, ivar),
                       CHF_CONST_INT(a_dir),
                       CHF_CONST_REAL(a_dtbydx));
  /**/

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      Real dense, denlo, denhi, denslope;
      RealVect veloc;

      dense    = a_rho(vof, ivar);
      denslope = a_dRho(vof, ivar);
      denlo    = a_rhoLo(vof, ivar);
      denhi    = a_rhoHi(vof, ivar);
     for (int idir = 0; idir < SpaceDim; idir++)
        {
          veloc[idir] = velcc(vof, idir);
        }

      FORT_POINTPREDADVECTTRANS(CHF_REAL(dense),
                                CHF_REAL(denlo),
                                CHF_REAL(denhi),
                                CHF_REAL(denslope),
                                CHF_REALVECT(veloc),
                                CHF_CONST_INT(a_dir),
                                CHF_CONST_REAL(a_dtbydx));

      a_rhoLo(vof, ivar) = denlo;
      a_rhoHi(vof, ivar) = denhi;
    }
}
/*******************/
void
EBPatchAdvect::
extrapolatePrim(EBFluxFAB&                       a_flux,
                EBCellFAB&                       a_primState,
                EBCellFAB                        a_slopePrim[SpaceDim],
                EBCellFAB                        a_slopeNLim[SpaceDim],
                Vector<BaseIVFAB<Real> * >&      a_coveredFluxMinu,
                Vector<BaseIVFAB<Real> * >&      a_coveredFluxPlus,
                const Vector<IntVectSet >&       a_coveredSetsMinu,
                const Vector<IntVectSet >&       a_coveredSetsPlus,
                const Vector<Vector<VolIndex> >& a_coveredFaceMinu,
                const Vector<Vector<VolIndex> >& a_coveredFacePlus,
                const EBCellFAB&                 a_flattening,
                const EBCellFAB&                 a_consState,
                const EBCellFAB&                 a_source,
                const Box&                       a_box,
                const DataIndex&                 a_dit,
                bool                             a_verbose)
{
  CH_TIME("EBPatchAdvect::extrapolatePrim");
  CH_assert(a_coveredFluxMinu.size() == SpaceDim);
  CH_assert(a_coveredFluxPlus.size() == SpaceDim);
  CH_assert(a_coveredFaceMinu.size() == SpaceDim);
  CH_assert(a_coveredFacePlus.size() == SpaceDim);
  CH_assert(a_coveredSetsMinu.size() == SpaceDim);
  CH_assert(a_coveredSetsPlus.size() == SpaceDim);

#if CH_SPACEDIM==2

  extrapolatePrim2D(m_primMinu, m_primPlus,
                    m_primState, a_slopePrim, a_slopeNLim,
                    a_flattening, a_consState, a_source,
                    a_box, a_dit, a_verbose);


#elif CH_SPACEDIM==3

  extrapolatePrim3D(m_primMinu,   m_primPlus,
                    m_primState,  a_slopePrim, a_slopeNLim,
                    a_flattening, a_consState, a_source,
                    a_box, a_dit,a_verbose);

#else
  MayDay::Error("bogus SpaceDim");
#endif

  Box faceBox[SpaceDim];
  Box bndryFaceBox[SpaceDim];
  //this keeps the fluxes from being calculated
  //on boundaries
  for (int dir1 = 0; dir1 < SpaceDim; ++dir1)
    {
      bndryFaceBox[dir1] = a_box;
      bndryFaceBox[dir1] &= m_domain;
      bndryFaceBox[dir1].surroundingNodes(dir1);

      faceBox[dir1] = a_box;
      faceBox[dir1].grow(dir1,1);
      faceBox[dir1] &= m_domain;
      faceBox[dir1].grow(dir1,-1);
      faceBox[dir1].surroundingNodes(dir1);
    }
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_extendStatePlusG4[idir].setVal(0.);
      m_extendStateMinuG4[idir].setVal(0.);
    }
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      //extrapolate to covered faces using updated extrapolated state
      EBPatchGodunov::extrapToCoveredFaces(m_extendStateMinuG4[faceDir],
                                           m_primMinu[faceDir],
                                           m_primPlus[faceDir],
                                           m_primState,
                                           a_coveredFaceMinu[faceDir],
                                           faceDir, Side::Lo, a_box);

      EBPatchGodunov::extrapToCoveredFaces(m_extendStatePlusG4[faceDir],
                                           m_primMinu[faceDir],
                                           m_primPlus[faceDir],
                                           m_primState,
                                           a_coveredFacePlus[faceDir],
                                           faceDir, Side::Hi, a_box);


      riemann(a_flux[faceDir], m_primPlus[faceDir], m_primMinu[faceDir],
              faceDir, faceBox[faceDir]);

      // Use the user supplied PhysBC object to obtain boundary fluxes
      m_bc->fluxBC(a_flux, a_primState, m_primMinu[faceDir],
                   Side::Lo,  m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[faceDir], faceDir);
      m_bc->fluxBC(a_flux, a_primState, m_primPlus[faceDir],
                   Side::Hi,  m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[faceDir], faceDir);
      //solve riemann problem between extended state and
      //value in vof to get covered flux.
      riemann(*a_coveredFluxMinu[faceDir],
              m_extendStateMinuG4[faceDir], m_primMinu[faceDir],
              a_coveredFaceMinu[faceDir], faceDir, Side::Lo, a_box);
      riemann(*a_coveredFluxPlus[faceDir],
               m_extendStatePlusG4[faceDir], m_primPlus[faceDir],
              a_coveredFacePlus[faceDir], faceDir, Side::Hi, a_box);
    }
}

/*******************/
void
EBPatchAdvect::
extrapolateBCG(EBFluxFAB&                       a_flux,
               // EBCellFAB&                       a_primState,
               EBCellFAB                        a_slopePrim[SpaceDim],
               EBCellFAB                        a_slopeNLim[SpaceDim],
               const EBCellFAB&                 a_flattening,
               const EBCellFAB&                 a_consState,
               const EBCellFAB&                 a_source,
               const Box&                       a_box,
               const DataIndex&                 a_dit,
               bool                             a_verbose)
{
  CH_TIME("EBPatchAdvect::extrapolateBCG");

  //Dan does this in EBPatchGodunov::extrapolatePrim which is called by EBLevelAdvect
  //we never go there for BCG so just use the virtual function
  //we don't need a_primState in interface any more but need a non-const in extrapolateBCG
  //use m_primState which is just sitting there as scratch space
  int logflag = 0;
  Box primBox = a_consState.box();
  // consToPrim(a_primState, a_consState, primBox, logflag);
  consToPrim(m_primState, a_consState, primBox, logflag);

  //Taylor extrapolation
  extrapolateBCG(m_primMinu,
                 m_primPlus,
                 // a_primState,
                 a_slopePrim,
                 a_slopeNLim,
                 a_flattening,
                 a_consState,
                 a_source,
                 a_box,
                 a_dit,
                 a_verbose);

  Box faceBox[SpaceDim];
  Box bndryFaceBox[SpaceDim];
  //this keeps the fluxes from being calculated on boundaries
  for (int dir1 = 0; dir1 < SpaceDim; ++dir1)
    {
      bndryFaceBox[dir1] = a_box;
      bndryFaceBox[dir1] &= m_domain;
      bndryFaceBox[dir1].surroundingNodes(dir1);

      faceBox[dir1] = a_box;
      faceBox[dir1].grow(dir1,1);
      faceBox[dir1] &= m_domain;
      faceBox[dir1].grow(dir1,-1);
      faceBox[dir1].surroundingNodes(dir1);
    }

  //solve the Riemann problem for plus and minus states at a face
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      riemann(a_flux[faceDir], m_primPlus[faceDir], m_primMinu[faceDir],
              faceDir, faceBox[faceDir]);

      // Use the user supplied PhysBC object to obtain boundary fluxes
      // m_bc->fluxBC(a_flux, a_primState, m_primMinu[faceDir],
      m_bc->fluxBC(a_flux, a_consState, m_primMinu[faceDir],
                   Side::Lo,  m_time, m_ebisBox, a_dit, a_box,
                   bndryFaceBox[faceDir], faceDir);
      // m_bc->fluxBC(a_flux, a_primState, m_primPlus[faceDir],
      m_bc->fluxBC(a_flux, a_consState, m_primPlus[faceDir],
                   Side::Hi,  m_time, m_ebisBox, a_dit, a_box,
                   bndryFaceBox[faceDir], faceDir);
    }
}

/******/
void
EBPatchAdvect::
extrapolateBCG(EBCellFAB           a_primMinu[SpaceDim],
               EBCellFAB           a_primPlus[SpaceDim],
               // EBCellFAB&          a_primState,
               EBCellFAB           a_slopesPrim[SpaceDim],
               EBCellFAB           a_slopesSeco[SpaceDim],
               const EBCellFAB&    a_flattening,
               const EBCellFAB&    a_consState,
               const EBCellFAB&    a_source,
               const Box&          a_box,
               const DataIndex&    a_dit,
               bool                a_verbose)
{
  CH_TIME("EBPatchAdvect::extrapolateBCG");

  Box slopeBoxG1 = grow(a_box, 1);
  Box slopeBoxG2 = grow(a_box, 2);
  Box slopeBoxG3 = grow(a_box, 3);
  slopeBoxG1 &= m_domain;
  slopeBoxG2 &= m_domain;
  slopeBoxG3 &= m_domain;
  Box cellBoxG2[SpaceDim];
  Box cellBoxG1[SpaceDim];
  Box faceBox[SpaceDim];
  Box bndryFaceBox[SpaceDim];
  int numSlop = numSlopes();
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      cellBoxG2[idir] = a_box;
      cellBoxG2[idir].grow(2);
      cellBoxG2[idir].grow(idir, -1);
      cellBoxG2[idir] &= m_domain;

      cellBoxG1[idir] = a_box;
      cellBoxG1[idir].grow(1);
      cellBoxG1[idir].grow(idir, -1);
      cellBoxG1[idir] &= m_domain;

      bndryFaceBox[idir] = a_box;
      bndryFaceBox[idir].grow(1);
      bndryFaceBox[idir] &= m_domain;
      bndryFaceBox[idir].surroundingNodes(idir);


      faceBox[idir] = a_box;
      faceBox[idir].grow(2);
      faceBox[idir] &= m_domain;
      faceBox[idir].grow(idir, -1);
      faceBox[idir].surroundingNodes(idir);
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      // a_slopesPrim[idir].define(m_ebisBox, slopeBoxG2, numSlop);
      // a_slopesSeco[idir].define(m_ebisBox, slopeBoxG2, numSlop);
      a_slopesPrim[idir].define(m_ebisBox, slopeBoxG1, numSlop);
      a_slopesSeco[idir].define(m_ebisBox, slopeBoxG1, numSlop);
      // a_slopesPrim[idir].define(m_ebisBox, a_box, numSlop);
      // a_slopesSeco[idir].define(m_ebisBox, a_box, numSlop);
      a_primPlus[idir].setVal(0.);
      a_primMinu[idir].setVal(0.);
    }

  //check flattening coefficients have correct box
  if (usesFlattening())
    {
      if (!a_flattening.getRegion().contains(slopeBoxG2))
        {
          MayDay::Error("flattening defined over too small a region");
        }
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      CH_assert(a_slopesPrim[idir].nComp() == numSlopes());
      // CH_assert(a_primState.nComp() >= numSlopes());
      CH_assert(a_consState.nComp() >= numSlopes());

      EBCellFAB normSlopes;
      EBCellFAB tanSlopes[SpaceDim-1];

      Box normSlopeBox;
      Box tanSlopeBox[SpaceDim-1];
      normSlopeBox = a_box;
      normSlopeBox.grow(idir, 1);
      normSlopeBox &= m_domain;

      // normSlopes.define(m_ebisBox, normSlopeBox, numSlop);
      normSlopes.define(m_ebisBox, normSlopeBox, numSlop);

      EBCellFAB deltaWL, deltaWR;
      // doSecondOrderSlopes(a_slopesPrim[idir],//second-order limited slope
      doSecondOrderSlopes(normSlopes,//second-order limited slope
                          deltaWL,//one-sided diff left
                          deltaWR,//one-sided diff right
                          a_slopesSeco[idir],//centered difference, non-limited
                          // a_primState,
                          a_consState,
                          idir,
                          normSlopeBox);
                          // slopeBoxG1);
                          // slopeBoxG2);
                          // a_box);//send the unghosted box which gets grown and intersected and then chopped with eblohicenter

      normalPred(a_primMinu[idir],
                 a_primPlus[idir],
                 // a_primState,
                 a_consState,
                 // a_slopesPrim[idir],
                 normSlopes,
                 m_dt/m_dx[idir],
                 idir,
                 // slopeBoxG2);
                 // a_box);
                 normSlopeBox);

      if (a_source.isDefined())
        {
          incrementWithSource(a_primMinu[idir],
                              a_source,
                              0.5*m_dt,
                              // slopeBoxG2);
                              // slopeBoxG1);
                              // a_box);
                              normSlopeBox);
          incrementWithSource(a_primPlus[idir],
                              a_source,
                              0.5*m_dt,
                              // slopeBoxG2);
                              // slopeBoxG1);
                              // a_box);
                              normSlopeBox);
        }

      Tuple<int,SpaceDim-1> tanDirs = PolyGeom::computeTanDirs(idir);
      for (int itan = 0; itan < SpaceDim-1; itan++)
        {
          int tanDir = tanDirs[itan];
          // tanSlopeBox[itan] = a_box;
          tanSlopeBox[itan] = normSlopeBox;
          tanSlopeBox[itan].grow(tanDir, 1);
          tanSlopeBox[itan] &= m_domain;

          tanSlopes[itan].define(m_ebisBox, tanSlopeBox[itan], numSlop);

          //Minion stability fix
          EBCellFAB& statePlusSource = m_primState;
          // EBCellFAB statePlusSource;
          // statePlusSource.copy(a_consState);
          statePlusSource.plus(a_source, 0.5*m_dt);

          // upwindSlope(a_slopesPrim[tanDir],
          upwindSlope(tanSlopes[itan],
                      statePlusSource,
                      tanDir,
                      // slopeBoxG3);
                      tanSlopeBox[itan]);
                      // a_box);

          transversePred(a_primMinu[idir],
                         a_primPlus[idir],
                         // a_primState,
                         a_consState,
                         // a_slopesPrim[tanDir],
                         tanSlopes[itan],
                         m_dt/m_dx[tanDir],
                         tanDir,
                         // slopeBoxG2);
                         // a_box);
                         normSlopeBox);
        }

      floorPrimitives(a_primMinu[idir], normSlopeBox);
      floorPrimitives(a_primPlus[idir], normSlopeBox);

    }
}

void
EBPatchAdvect::
upwindSlope(EBCellFAB&       a_slopeUpWi,
            const EBCellFAB& a_primState,
            const int&       a_dir,
            const Box&       a_box)
{
  CH_TIME("EBPatchAdvect::upwindSlope");
  int numSlope = numSlopes();
  CH_assert(a_slopeUpWi.nComp() == numSlopes());
  CH_assert(a_primState.nComp() >= numSlopes());

  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;
  eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
               a_box, m_domain, a_dir);

  // a_slopeUpWi.define(m_ebisBox, entireBox, numSlope);

  BaseFab<Real>& regSlopeUpWi       = a_slopeUpWi.getSingleValuedFAB();
  const BaseFab<Real>& regPrimState = a_primState.getSingleValuedFAB();
  const BaseFab<Real>& regNormalVel = m_normalVelPtr->getSingleValuedFAB();


  FORT_UPWINDDIFFS(CHF_FRA(regSlopeUpWi),
                   CHF_CONST_FRA(regPrimState),
                   CHF_CONST_FRA(regNormalVel),
                   CHF_CONST_INT(numSlope),
                   CHF_CONST_INT(a_dir),
                   CHF_BOX(loBox),
                   CHF_CONST_INT(hasLo),
                   CHF_BOX(hiBox),
                   CHF_CONST_INT(hasHi),
                   CHF_BOX(centerBox));

  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          for (int ivar = 0; ivar < numSlopes(); ivar++)
            {
              bool hasFacesLeft, hasFacesRigh;
              Real dql, dqr, dqc;
              bool verbose =false;
              pointGetSlopesUpwind(dql, dqr, dqc,
                                   hasFacesLeft,
                                   hasFacesRigh,
                                   vof, a_primState, a_dir, ivar, verbose);


              Real velPt = (*m_normalVelPtr)(vof, a_dir);
              if (velPt > TOL)
                {
                  dqc = dql;
                }
              else if (velPt < -TOL)
                {
                  dqc = dqr;
                }
              else
                {
                  dqc = 0.5*(dql+dqr);
                }

              a_slopeUpWi(vof, ivar) = dqc;
            } //end loop over variables
        }
    }//end loop over vofs
  {
    CH_TIME("bndry_slopes");
    //want to be able to call this in test codes
    if (m_isBCSet)
      {
        m_bc->setBndrySlopes(a_slopeUpWi, a_primState, m_ebisBox, entireBox, a_dir);
      }
  }
}
/*************Need this because don't want to limit transverse slopes in BCG****************/
void EBPatchAdvect::
pointGetSlopesUpwind(Real&            a_dql,
                     Real&            a_dqr,
                     Real&            a_dqc,
                     bool&            a_hasFacesLeft,
                     bool&            a_hasFacesRigh,
                     const VolIndex&  a_vof,
                     const EBCellFAB& a_primState,
                     const int&       a_dir,
                     const int&       a_ivar,
                     const bool&      a_verbose)
{
  CH_TIME("EBPG::pointGetSlopes");
  const IntVect&  iv = a_vof.gridIndex();
  //one-sided diffs on domain bndry
  const Box& domainBox = m_domain.domainBox();
  // const EBCellFAB& velcc = *m_normalVelPtr;
  bool onLeftDomain = (iv[a_dir] == domainBox.smallEnd(a_dir));
  bool onRighDomain = (iv[a_dir] == domainBox.bigEnd(a_dir)  );
  VolIndex vofLeft;
  VolIndex vofRigh;
  a_hasFacesLeft = (m_ebisBox.numFaces(a_vof, a_dir, Side::Lo)==1) && !onLeftDomain;
  a_hasFacesRigh = (m_ebisBox.numFaces(a_vof, a_dir, Side::Hi)==1) && !onRighDomain;

  if (a_hasFacesLeft)
    {
      Vector<FaceIndex> facesLeft =
        m_ebisBox.getFaces(a_vof, a_dir, Side::Lo);
      vofLeft = facesLeft[0].getVoF(Side::Lo);
    }

  if (a_hasFacesRigh)
    {
      Vector<FaceIndex> facesRigh =
        m_ebisBox.getFaces(a_vof, a_dir, Side::Hi);
      vofRigh = facesRigh[0].getVoF(Side::Hi);
    }

  Real valCent = a_primState(a_vof, a_ivar);
  Real valLeft = 0.0;
  Real valRigh = 0.0;
  if (a_hasFacesLeft)
    {
      valLeft = a_primState(vofLeft, a_ivar);
    }
  if (a_hasFacesRigh)
    {
      valRigh = a_primState(vofRigh, a_ivar);
    }

  //at irregular cells, if one side does not exist,
  //set all to one-sided diffs.
  //if neither exists, set all slopes to zero
  if (a_hasFacesLeft)
    {
      a_dql = valCent - valLeft;
    }
  if (a_hasFacesRigh)
    {
      a_dqr = valRigh - valCent;
    }
  if (a_hasFacesLeft && a_hasFacesRigh)
    {
      a_dqc = 0.5*(a_dql+a_dqr);
    }
  else if (!a_hasFacesLeft && !a_hasFacesRigh)
    {
      a_dql = 0.0;
      a_dqr = 0.0;
      a_dqc = 0.0;
    }
  else if (a_hasFacesLeft && !a_hasFacesRigh)
    {
      a_dqr = a_dql;
      //no higher-order one sided diffs
      a_dqc = a_dql;
    }
  else if (a_hasFacesRigh && !a_hasFacesLeft)
    {
      a_dql = a_dqr;
      //no higher-order one sided diffs
      a_dqc = a_dqr;
    }
  else
    {
      MayDay::Error("EBPatchGodunov::pointGetSlopes -- missed a case");
    }
  if (a_verbose)
    {
      pout() << "  a_vof="    << a_vof
             << ", vofLeft= " << vofLeft
             << ", vofRigh= " << vofRigh  << endl;
      pout() << "  hasFacesLeft=" << a_hasFacesLeft
             << ", hasFacesRigh=" << a_hasFacesRigh;
      pout() << ", valLeft=" << valLeft
             << ", valRigh=" << valRigh
             << ", valCent=" << valCent << endl;
      pout() << "  a_dql=" << a_dql
             << ", a_dqr=" << a_dqr
             << ", a_dqc=" << a_dqc << endl;
    }
  if (a_dqc != a_dqc)
    {
      MayDay::Error("EBPatchAdvect::pointGetSlopesUpwind -- a_dqc != a_dqc");
    }
}
/****************/
void
EBPatchAdvect::
advectiveDerivative(EBCellFAB&                       a_uDotDelRho,
                    const EBFluxFAB&                 a_faceRho,
                    const EBFluxFAB&                 a_faceVel,
                    const Vector<BaseIVFAB<Real>*> & a_coveredRhoLo,
                    const Vector<BaseIVFAB<Real>*> & a_coveredRhoHi,
                    const Vector<BaseIVFAB<Real>*> & a_coveredVelLo,
                    const Vector<BaseIVFAB<Real>*> & a_coveredVelHi,
                    const Vector<Vector<VolIndex> >& a_coveredFaceLo,
                    const Vector<Vector<VolIndex> >& a_coveredFaceHi,
                    const Box&                       a_box)
{
  CH_TIME("EBPatchAdvect::advectiveDerivative");
  int ncomp = a_faceRho.nComp();
  CH_assert(m_isDefined);
  CH_assert(m_isBoxSet);

  CH_assert(a_coveredRhoLo.size() == SpaceDim);
  CH_assert(a_coveredVelLo.size() == SpaceDim);
  CH_assert(a_coveredRhoHi.size() == SpaceDim);
  CH_assert(a_coveredVelHi.size() == SpaceDim);
  CH_assert(a_coveredFaceHi.size() == SpaceDim);
  CH_assert(a_coveredFaceLo.size() == SpaceDim);

  //set udotdelrho to zero.  the fortran is additive
  a_uDotDelRho.setVal(0.);
  BaseFab<Real>& regUDotDelRho = a_uDotDelRho.getSingleValuedFAB();
  //compute udotdelu everywhere as regular.
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      const EBFaceFAB& faceRho = a_faceRho[faceDir];
      const EBFaceFAB& faceVel = a_faceVel[faceDir];

      const BaseFab<Real>& regFaceRho = faceRho.getSingleValuedFAB();
      const BaseFab<Real>& regFaceVel = faceVel.getSingleValuedFAB();

      //this does udotdelvel += 0.5*(uhigh+ulow)*(velhigh-vellow)/dx (non-conservative)
      //this does udotdelrho +=   (velhigh*rhohigh-vellow*rholow)/dx (    conservative)
      //where high == ivec + half*e^facedir
      FORT_ADVECTIVEF(CHF_FRA(regUDotDelRho),
                      CHF_CONST_FRA(regFaceRho),
                      CHF_CONST_FRA1(regFaceVel, 0),
                      CHF_CONST_INT(faceDir),
                      CHF_CONST_INT(ncomp),
                      CHF_CONST_REAL(m_dx[faceDir]),
                      CHF_BOX(a_box),
                      CHF_INT(s_doingVel));
    }

  //update the irregular vofs
  for (int ivof = 0; ivof < m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      //the set can be bigger than the box for performance reasons.
      if (a_box.contains(vof.gridIndex()))
        {
          // if (s_verbose && vof.gridIndex()[0]==12 && vof.gridIndex()[1]==8){pout() << "nonconservative advection " << vof << endl;}
          for (int ivar = 0; ivar < ncomp; ivar++)
            {
              //udelrho was set in regular update.  we reset it
              // to zero and recalc.
              Real uDelRhoPt = 0.0;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, idir, Side::Lo);
                  Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, idir, Side::Hi);
                  Real rhoLo = 0.0;
                  Real velLo = 0.0;
                  if (facesLo.size() > 0)
                    {
                      for (int iface = 0; iface < facesLo.size(); iface++)
                        {
                          rhoLo += a_faceRho[idir](facesLo[iface], ivar);
                          velLo += a_faceVel[idir](facesLo[iface], 0);
                        }
                      velLo /=  facesLo.size();
                      rhoLo /=  facesLo.size();
                    }
                  else
                    {
                      const BaseIVFAB<Real>& coveredRhoLo = *a_coveredRhoLo[idir];
                      const BaseIVFAB<Real>& coveredVelLo = *a_coveredVelLo[idir];
                      rhoLo = coveredRhoLo(vof, ivar);
                      velLo = coveredVelLo(vof, 0);
                    }
                  Real rhoHi = 0.0;
                  Real velHi = 0.0;
                  if (facesHi.size() > 0)
                    {
                      for (int iface = 0; iface < facesHi.size(); iface++)
                        {
                          rhoHi += a_faceRho[idir](facesHi[iface], ivar);
                          velHi += a_faceVel[idir](facesHi[iface], 0);
                        }
                      velHi /=  facesHi.size();
                      rhoHi /=  facesHi.size();
                    }
                  else
                    {
                      const BaseIVFAB<Real>& coveredRhoHi = *a_coveredRhoHi[idir];
                      const BaseIVFAB<Real>& coveredVelHi = *a_coveredVelHi[idir];
                      rhoHi = coveredRhoHi(vof, ivar);
                      velHi = coveredVelHi(vof, 0);
                    }

                  Real velAve = 0.5*(velHi + velLo);
                  Real rhoDiff = rhoHi - rhoLo;
                  uDelRhoPt += velAve*rhoDiff/m_dx[idir];
                } //end loop over directions
              a_uDotDelRho(vof, ivar) = uDelRhoPt;
              // if (s_verbose && vof.gridIndex()[0]==12 && vof.gridIndex()[1]==8){pout() << "   uDotDelRho " << uDelRhoPt << endl;}
            }//end loop over variables
        } //if vof is in this box
    }//end loop over vofs
}
/****************/
void
EBPatchAdvect::
advectiveDerivative(EBCellFAB&                       a_uDotDelRho,
                    const EBFluxFAB&                 a_faceRho,
                    const EBFluxFAB&                 a_faceVel,
                    const Box&                       a_box)
{
  CH_TIME("EBPatchAdvect::advectiveDerivative");
  int ncomp = a_faceRho.nComp();
  CH_assert(m_isDefined);
  CH_assert(m_isBoxSet);

  //set udotdelrho to zero.  the fortran is additive
  a_uDotDelRho.setVal(0.);
  BaseFab<Real>& regUDotDelRho = a_uDotDelRho.getSingleValuedFAB();
  //compute udotdelu everywhere as regular.
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      const EBFaceFAB& faceRho = a_faceRho[faceDir];
      const EBFaceFAB& faceVel = a_faceVel[faceDir];

      const BaseFab<Real>& regFaceRho = faceRho.getSingleValuedFAB();
      const BaseFab<Real>& regFaceVel = faceVel.getSingleValuedFAB();

      //this does udotdelvel += 0.5*(uhigh+ulow)*(velhigh-vellow)/dx (non-conservative)
      //this does udotdelrho +=   (velhigh*rhohigh-vellow*rholow)/dx (    conservative)
      //where high == ivec + half*e^facedir
      FORT_ADVECTIVEF(CHF_FRA(regUDotDelRho),
                      CHF_CONST_FRA(regFaceRho),
                      CHF_CONST_FRA1(regFaceVel, 0),
                      CHF_CONST_INT(faceDir),
                      CHF_CONST_INT(ncomp),
                      CHF_CONST_REAL(m_dx[faceDir]),
                      CHF_BOX(a_box),
                      CHF_INT(s_doingVel));
    }

}
/****************/
void
EBPatchAdvect::
setVelocities(const EBCellFAB& a_normalVel,
              const EBFluxFAB& a_advectionVel)

{
  m_normalVelPtr    = &a_normalVel;
  m_advectionVelPtr = &a_advectionVel;
  m_isVelSet = true;
}
/******/
void
EBPatchAdvect::
updatePrim(EBCellFAB&              a_primMinu,
           EBCellFAB&              a_primPlus,
           const EBFaceFAB&        a_primFace,
           const BaseIVFAB<Real>&  a_coveredPrimMinu,
           const BaseIVFAB<Real>&  a_coveredPrimPlus,
           const Vector<VolIndex>& a_coveredFaceMinu,
           const Vector<VolIndex>& a_coveredFacePlus,
           const int&              a_faceDir,
           const Box&              a_box,
           const Real&             a_scale)
{
  CH_TIME("EBPatchAdvect::updatePrim");
  //uses flux as prims
  //update the regular vofs.
  //first we have to copy the original state
  //so that the irregular stuff can be updated later.

  int numPrim = numPrimitives();
  //save state so that the regular can overwrite irregular cells
  Vector<Vector<Real> > cacheMinu;
  Vector<Vector<Real> > cachePlus;
  {
    CH_TIME("EBPatchAdvect::cache");
    cacheEBCF(cacheMinu, a_primMinu);
    cacheEBCF(cachePlus, a_primPlus);
  }

  const EBCellFAB& normalVel = *m_normalVelPtr;
  BaseFab<Real>& regPrimMinu = a_primMinu.getSingleValuedFAB();
  BaseFab<Real>& regPrimPlus = a_primPlus.getSingleValuedFAB();
  const BaseFab<Real>& regPrimFace      =    a_primFace.getSingleValuedFAB();

  {
    CH_TIME("EBPatchAdvect::fortran");
    const BaseFab<Real>& regNormVel = normalVel.getSingleValuedFAB();
    FORT_ADVECTUPDATE(CHF_BOX(a_box),
                      CHF_FRA(regPrimMinu),
                      CHF_FRA(regPrimPlus),
                      CHF_CONST_FRA(regPrimFace),
                      CHF_CONST_FRA(regNormVel),
                      CHF_CONST_INT(a_faceDir),
                      CHF_CONST_INT(numPrim),
                      CHF_CONST_REAL(a_scale),
                      CHF_BOX(a_box));
  }

  //put irregular cells back to their original state
  {
    CH_TIME("EBPatchAdvect::uncache");
    uncacheEBCF(a_primMinu ,cacheMinu);
    uncacheEBCF(a_primPlus ,cachePlus);
  }

  {
    CH_TIME("EBPatchAdvect::irregularLoop");
    for (int ivof = 0; ivof < m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        //the set can be bigger than the box for performance reasons.
        if (a_box.contains(vof.gridIndex()))
          {
            Real vel = normalVel(vof, a_faceDir);

            for (int ivar = 0; ivar < numPrim; ivar++)
              {

                Real primLo = 0;
                {
                  Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, a_faceDir, Side::Lo);
                  if (facesLo.size() > 0)
                    {
                      for (int iface = 0; iface < facesLo.size(); iface++)
                        {
                          primLo += a_primFace(facesLo[iface], ivar);
                        }
                      primLo /= facesLo.size();
                    }
                  else
                    {
                      primLo = a_coveredPrimMinu(vof, ivar);
                    }
                }
                Real primHi = 0;
                {
                  Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, a_faceDir, Side::Hi);
                  if (facesHi.size() > 0)
                    {
                      for (int iface = 0; iface < facesHi.size(); iface++)
                        {
                          primHi += a_primFace(facesHi[iface], ivar);
                        }
                      primHi /= facesHi.size();
                    }
                  else
                    {
                      primHi = a_coveredPrimPlus(vof, ivar);
                    }
                }
                Real primDiff = primHi - primLo;

                //state changed in regular update.
                //but we cached and uncached
                //scale holds 0.5*dt/dx in 2d.
                a_primMinu(vof, ivar) -=  a_scale*vel*primDiff;
                a_primPlus(vof, ivar) -=  a_scale*vel*primDiff;

                //debug set to inputs
                //a_primMinu(vof, ivar) = origPrimMinu;
                //a_primPlus(vof, ivar) = origPrimPlus;
                //end debug
              }
          }
      }
  }

  {
    CH_TIME("EBPatchAdvect::floors");
    floorPrimitives(a_primMinu, a_box);
    floorPrimitives(a_primPlus, a_box);
  }
}
/*****/
bool
EBPatchAdvect::
usesFlattening() const
{
  return false;
}
/******/
bool
EBPatchAdvect::
usesArtificialViscosity() const
{
  return false;
}
/******/
bool
EBPatchAdvect::
usesFourthOrderSlopes() const
{
  return false;
}
/******/
void
EBPatchAdvect::
consToPrim(EBCellFAB&       a_primState,
           const EBCellFAB& a_consState,
           const Box&       a_box,
           int              a_logflag,
           bool             a_verbose)
{
  CH_TIME("EBPatchAdvect::consToPrim");
  CH_assert(numPrimitives() == numConserved());
  Interval interv(0, numConserved()-1);
  a_primState.copy(a_box, interv, a_box, a_consState, interv);
}
/******/
void
EBPatchAdvect::
consToPrim(BaseIVFAB<Real>&  a_primState,
           const EBCellFAB&  a_consState,
           const IntVectSet& a_ivs)
{
  CH_TIME("EBPatchAdvect::consToPrimIrr");
  CH_assert(numPrimitives() == numConserved());
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      for (int ivar = 0; ivar < numPrimitives(); ivar++)
        {
          a_primState(vofit(), ivar) = a_consState(vofit(), ivar);
        }
    }
}
/******/
void
EBPatchAdvect::
consToPrim(BaseIVFAB<Real>&        a_primState,
           const BaseIVFAB<Real>&  a_consState,
           const IntVectSet&       a_ivs)
{
  CH_TIME("EBPatchAdvect::consToPrimIrrIrr");
  CH_assert(numPrimitives() == numConserved());
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      for (int ivar = 0; ivar < numPrimitives(); ivar++)
        {
          a_primState(vofit(), ivar) = a_consState(vofit(), ivar);
        }
    }
}
/******/
void
EBPatchAdvect::
primToCons(EBCellFAB&       a_consState,
           const EBCellFAB& a_primState,
           const Box&       a_box)
{
  CH_TIME("EBPatchAdvect::primToConsRegReg");
  CH_assert(numPrimitives() == numConserved());
  Interval interv(0, numConserved()-1);
  a_consState.copy(a_box, interv, a_box, a_primState, interv);
}
/******/
void
EBPatchAdvect::
primToCons(BaseIVFAB<Real>&       a_consState,
           const BaseIVFAB<Real>& a_primState,
           const IntVectSet&      a_ivs)
{
  CH_TIME("EBPatchAdvect::primToConsRegIrr");
  CH_assert(numPrimitives() == numConserved());
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      for (int ivar = 0; ivar < numPrimitives(); ivar++)
        {
          a_consState(vofit(), ivar) = a_primState(vofit(), ivar);
        }
    }
}
///returns primitives not fluxes
/***/
void
EBPatchAdvect::
riemann(EBFaceFAB&       a_primGdnv,
        const EBCellFAB& a_primLeft,
        const EBCellFAB& a_primRigh,
        const int&       a_faceDir,
        const Box&       a_box)
{
  CH_TIME("EBPatchAdvect::riemann");
  //these two need to be set
  CH_assert((s_doingVel == 1) || (s_doingVel == 0));
  //  CH_assert((s_curComp >= 0) && (s_curComp < SpaceDim);
  CH_assert(s_curComp >= 0);
  CH_assert(m_isVelSet);
  CH_assert(!m_domain.isPeriodic(a_faceDir));

  // const EBCellFAB& normalVel = *m_normalVelPtr;
  const EBFluxFAB& advectVel = *m_advectionVelPtr;
  const EBFaceFAB& advectVelFaceFAB = advectVel[a_faceDir];
  int nPrim = numPrimitives();

  Box cellBox = enclosedCells(a_box);

  //explicitly cast the left and right states to modifiable references
  //because we need to shift them to faces and then shift them back.
  BaseFab<Real>& regPrimRigh = (BaseFab<Real>&)a_primRigh.getSingleValuedFAB();
  BaseFab<Real>& regPrimLeft = (BaseFab<Real>&)a_primLeft.getSingleValuedFAB();
  BaseFab<Real>& regPrimGdnv = a_primGdnv.getSingleValuedFAB();
  // const BaseFab<Real>& regNormalVel = normalVel.getSingleValuedFAB();
  const BaseFab<Real>& regAdvectVel = advectVelFaceFAB.getSingleValuedFAB();

  //find the regular part of qgdnv.
  FORT_ADVECTRIEMANN(CHF_BOX(a_box),
                     CHF_FRA(regPrimGdnv),
                     CHF_CONST_FRA(regPrimLeft),
                     CHF_CONST_FRA(regPrimRigh),
                     // CHF_CONST_FRA(regNormalVel),
                     CHF_CONST_FRA1(regAdvectVel, 0),//only 1 (normal) component at face
                     CHF_CONST_INT(a_faceDir),
                     CHF_CONST_INT(nPrim),
                     CHF_CONST_INT(s_curComp),
                     CHF_CONST_INT(s_doingVel));

  //the box sent into this is face-centered.
  //we need to use the cell-centered one it surrounds.
  Vector<FaceIndex> multiFaces =  m_ebisBox.getEBGraph().getMultiValuedFaces(a_faceDir, cellBox);
  for (int iface = 0; iface< multiFaces.size(); iface++)
    {
      const FaceIndex& face = multiFaces[iface];
      if (!face.isBoundary())
        {
          VolIndex vofl = face.getVoF(Side::Lo);
          VolIndex vofr = face.getVoF(Side::Hi);

          Real velFace = advectVelFaceFAB(face, 0);
          for (int ivar = 0; ivar < nPrim; ivar++)
            {
              if (velFace > TOL)
                {
                  a_primGdnv(face, ivar) = a_primLeft(vofl, ivar);
                }
              else if (velFace < -TOL)
                {
                  a_primGdnv(face, ivar) = a_primRigh(vofr, ivar);
                }
              else
                {
                  if ((s_doingVel == 1) && (s_curComp == a_faceDir))
                    {
                      a_primGdnv(face, ivar) = 0.0;
                    }
                  else
                    {
                      a_primGdnv(face, ivar) = 0.5*(a_primRigh(vofr, ivar) + a_primLeft(vofl, ivar));
                    }
                }
            }

        }
    }
}
/***/
///returns primitives not fluxes
/*****************************/
void
EBPatchAdvect::
riemann(BaseIVFAB<Real>&        a_coveredPrim,
        const BaseIVFAB<Real>&  a_exteState,
        const EBCellFAB&        a_primState,
        const Vector<VolIndex>& a_vofset,
        const int&              a_faceDir,
        const Side::LoHiSide&   a_sd,
        const Box&       a_box)
{
  CH_TIME("EBPatchAdvect::riemannIrr");
  //this holds velocity at time n
  const EBCellFAB& normalVel = *m_normalVelPtr;
  int nPrim = numPrimitives();
  for (int ivof = 0; ivof < a_vofset.size(); ivof++)
    {
      const VolIndex& vof = a_vofset[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          Real  vel = normalVel(vof, a_faceDir);

          for (int ivar = 0; ivar < nPrim; ivar++)
            {
              Real leftState, righState;
              if (a_sd == Side::Lo)
                {
                  leftState = a_exteState(vof, ivar);
                  righState = a_primState(vof, ivar);
                }
              else
                {
                  leftState = a_primState(vof, ivar);
                  righState = a_exteState(vof, ivar);
                }

              if (vel > TOL)
                {
                  a_coveredPrim(vof, ivar) = leftState;
                }
              else if (vel < -TOL)
                {
                  a_coveredPrim(vof, ivar) = righState;
                }
              else
                {
                  a_coveredPrim(vof, ivar) = 0.5*(righState+leftState);
                }
            }
        }
    }
}
/*****************************/
EBPatchAdvect::
EBPatchAdvect():EBPatchGodunov()
{
  m_isVelSet        = false;
  m_isMaxMinSet     = false;
  m_advectionVelPtr = NULL;
  m_normalVelPtr    = NULL;
  m_useLimiting            = true;
  m_useFourthOrderSlopes   = false;
  m_useFlattening          = false;
}
/*****************************/
void
EBPatchAdvect::
useLimiting(bool a_limiting)
{
  //sets flattening and fourth order slopes to false
  setSlopeParameters(false, false, a_limiting);
}
/*****************************/
EBPatchAdvect::
~EBPatchAdvect()
{
}
/*****************************/
Vector<string>
EBPatchAdvect::
stateNames()
{
  Vector<string> retval(numConserved());
  for (int ivar = 0; ivar < numConserved(); ivar++)
    {
      char strname[100];
      sprintf(strname, "consname%d",ivar);
      retval[ivar] = string(strname);
    }
  return retval;
}
/*****************************/
Vector<string>
EBPatchAdvect::
primNames()
{
  Vector<string> retval(numConserved());
  for (int ivar = 0; ivar < numConserved(); ivar++)
    {
      char strname[100];
      sprintf(strname, "primname%d",ivar);
      retval[ivar] = string(strname);
    }
  return retval;
}
/*****************************/
void
EBPatchAdvect::
computeEBIrregFlux(BaseIVFAB<Real>&  a_ebIrregFlux,
                   const EBCellFAB&  a_primState,
                   const EBCellFAB   a_slopePrim[SpaceDim],
                   const IntVectSet& a_irregIVS,
                   const EBCellFAB&  a_source)
{
  CH_TIME("EBPatchAdvect::computeEBIrregFlux");
  // loop over the cells with irregular faces and extrapolate
  // half-time variables to the faces using crude Taylor series
  // approximation.
  int numPrim = numPrimitives();
  const EBCellFAB& ccVel = *m_normalVelPtr;
  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if (m_validBox.contains(vof.gridIndex()))
        {
          IntVect iv = vof.gridIndex();
          RealVect bndryCentroid = m_ebisBox.bndryCentroid(vof);
          bndryCentroid *= m_dx;
          for (int comp=0; comp<numPrim; comp++)
            {
              Real traceVal = 0.0;
              for (int dir=0; dir<SpaceDim; dir++)
                {
                  Real gradMult = bndryCentroid[dir] - 0.5*m_dt*ccVel(vof, dir);
                  traceVal -= gradMult*(a_slopePrim[dir])(vof, comp)/m_dx[dir];
                }
              traceVal += 0.5*m_dt*a_source(vof, comp);
              traceVal += a_primState(vof, comp);
              a_ebIrregFlux(vof, comp) = traceVal;
            } // end loop over components
        }
    } // end loop over vofs
}
/*****************************/
void
EBPatchAdvect::
averageVelToCC(EBCellFAB&                        a_normalVel,
               const EBFluxFAB&                  a_advectionVel,
               const Vector<BaseIVFAB<Real> * >& a_coveredVeloLo,
               const Vector<BaseIVFAB<Real> * >& a_coveredVeloHi,
               const Vector<Vector<VolIndex> >&  a_coveredFaceLo,
               const Vector<Vector<VolIndex> >&  a_coveredFaceHi,
               const Box&                        a_box) const
{
  CH_TIME("EBPatchAdvect::averageVelToCC");
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      const EBFaceFAB& faceNormVel = a_advectionVel[faceDir];

      //treat every cell as regular
      BaseFab<Real>&       regCellVel = a_normalVel.getSingleValuedFAB();
      const BaseFab<Real>& regFaceVel = faceNormVel.getSingleValuedFAB();
      FORT_AVEFACETOCELL(CHF_FRA(regCellVel),
                         CHF_CONST_FRA1(regFaceVel,0),
                         CHF_CONST_INT(faceDir),
                         CHF_BOX(a_box));
    }
  //correct way
  {
    CH_TIME("EBPatchAdvect::irregularLoop");
    for (int ivof = 0; ivof < m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        //the set can be bigger than the box for performance reasons.
        if (a_box.contains(vof.gridIndex()))
          {
            for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
              {
                Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, faceDir, Side::Lo);
                Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, faceDir, Side::Hi);
                Real velLo= 0.0;
                if (facesLo.size() > 0)
                  {
                    for (int iface = 0; iface < facesLo.size(); iface++)
                      {
                        velLo += a_advectionVel[faceDir](facesLo[iface], 0);
                      }
                    velLo /=  facesLo.size();
                  }
                else
                  {
                    const BaseIVFAB<Real>& coveredVelLo = *a_coveredVeloLo[faceDir];
                    velLo = coveredVelLo(vof, 0);
                  }
                Real velHi= 0.0;
                if (facesHi.size() > 0)
                  {
                    for (int iface = 0; iface < facesHi.size(); iface++)
                      {
                        velHi += a_advectionVel[faceDir](facesHi[iface], 0);
                      }
                    velHi /=  facesHi.size();
                  }
                else
                  {
                    const BaseIVFAB<Real>& coveredVelHi = *a_coveredVeloHi[faceDir];
                    velHi = coveredVelHi(vof, 0);
                  }

                Real velAve = 0.5*(velHi + velLo);
                a_normalVel(vof, faceDir) = velAve;
              }//end loop over directions
          }
      }//end loop over irregular vofs
  }
}
#include "NamespaceFooter.H"
