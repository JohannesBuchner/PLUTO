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

#include "EBDebugOut.H"
#include "DebugOut.H"
#include "EBPatchPolytropic.H"
#include "EBPatchPolytropicF_F.H"
#include "EBPatchGodunovF_F.H"
#include "EBLGIntegrator.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "EBLoHiCenter.H"
#include "EBPatchGodunov.H"
#include "Stencils.H"
#include "EBArith.H"
#include "EBAMRIO.H"
#include "PolyGeom.H"
#include "TensorCFInterp.H" //just for gradIndex
#include "ParmParse.H"
#include <cstdio>

#include "NamespaceHeader.H"

void
EBPatchPolytropic::
setValidBox(const Box& a_validBox,
            const EBISBox& a_ebisBox,
            const IntVectSet& a_coarseFineIVS,
            const Real& a_time,
            const Real& a_dt)
{
  EBPatchGodunov::setValidBox(a_validBox, a_ebisBox, a_coarseFineIVS, a_time, a_dt);
  Box slopeBoxG2 = grow(a_validBox, 2);
  slopeBoxG2 &= m_domain;
  IntVectSet ivs = m_ebisBox.getIrregIVS(slopeBoxG2);
  m_vofsStencil.define(ivs, m_ebisBox.getEBGraph(), 1);
  for(VoFIterator vofit(ivs, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      int radius = 1;
      EBArith::getAllVoFsInMonotonePath(m_vofsStencil(vofit(), 0), vofit(), m_ebisBox, radius);
    }
}
bool
EBPatchPolytropic::
checkCoveredFaceStencil(const VolIndex& a_vof, 
                        const Vector<VolIndex>& a_vofsStencil,
                        const RealVect& a_normal, 
                        const int&      a_faceDir,
                        const Side::LoHiSide& a_sd)
{
  bool hasAllVoFs = false;
  if(SpaceDim==2)
    {
      hasAllVoFs = checkCoveredFaceStencil2d(a_vof, a_vofsStencil, a_normal, a_faceDir, a_sd);
    }
  else
    {
      hasAllVoFs = checkCoveredFaceStencil3d(a_vof, a_vofsStencil, a_normal, a_faceDir, a_sd);
    }
  return hasAllVoFs;
}
bool
EBPatchPolytropic::
checkCoveredFaceStencil2d(const VolIndex& a_vof, 
                          const Vector<VolIndex>& a_vofsStencil,
                          const RealVect& a_normal, 
                          const int&      a_faceDir,
                          const Side::LoHiSide& a_sd)
{
  CH_assert(SpaceDim==2);
  bool hasAllVoFs;
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
  ivSten[0][0] = a_vof.gridIndex();
  ivSten[1][0] = ivSten[0][0] + signNorm*BASISV(a_faceDir);
  ivSten[0][1] = ivSten[0][0] - signNorm*BASISV(a_faceDir) + signTang*BASISV(tangenDir) ;;
  ivSten[1][1] = ivSten[0][0]                              + signTang*BASISV(tangenDir) ;

  for (int ix = 0; ix < 2; ix++)
    {
      for (int iy = 0; iy < 2; iy++)
        {
          hasVoF[ix][iy] = EBArith::isVoFHere(vof[ix][iy], a_vofsStencil, ivSten[ix][iy]);
        }
    }

  int dirBigNorm, dirLitNorm;
  bool hasVoFBigNorm,  hasVoFLitNorm;
  int  signBigNorm, signLitNorm;
  Real eps = 1.0e-12;
  if ((Abs(a_normal[a_faceDir]) < eps) || (Abs(a_normal[tangenDir]/a_normal[a_faceDir]) > m_dx[tangenDir]/m_dx[a_faceDir]))
    {
      dirBigNorm = tangenDir;
      dirLitNorm = a_faceDir;

      signBigNorm   = signTang;
      signLitNorm   = signNorm;
      hasVoFBigNorm = hasVoF[0][1];
    }
  else
    {
      dirBigNorm = a_faceDir;
      dirLitNorm = tangenDir;

      signBigNorm   = signNorm;
      signLitNorm   = signTang;
      hasVoFBigNorm = hasVoF[1][0];

    }
  hasVoFLitNorm = hasVoF[1][1];
  hasAllVoFs = (hasVoFLitNorm && hasVoFBigNorm);

  return hasAllVoFs;
}
/*********/
bool
EBPatchPolytropic::
checkCoveredFaceStencil3d(const VolIndex& a_vof, 
                          const Vector<VolIndex>& a_vofsStencil,
                          const RealVect& a_normal, 
                          const int&      a_faceDir,
                          const Side::LoHiSide& a_sd)
{
  CH_assert(SpaceDim==3);
  bool hasAllVoFs;
  Tuple<int,CH_SPACEDIM-1> tangenDir = PolyGeom::computeTanDirs(a_faceDir);
  Real anormNorm = Abs(a_normal[a_faceDir]);
  Real anormTan[2];
  int signNorm = 1;
  if (a_normal[a_faceDir] < 0.0) signNorm = -1;
  int signTang[2];
  for (int itan = 0; itan < 2; itan++)
    {
      int tanDir = tangenDir[itan];
      anormTan[itan] = Abs(a_normal[tanDir]);
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
      else if (a_normal[tanDir] < 0.0)
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
  VolIndex vofSten[2][2];

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

  for (int ix = 0; ix < 2; ix++)
    {
      for (int iy = 0; iy < 2; iy++)
        {
          hasVoF[ix][iy] =  EBArith::isVoFHere(vofSten[ix][iy], a_vofsStencil, ivSten[ix][iy]);
        }
    }
  hasAllVoFs = hasVoF[0][0] && hasVoF[1][0] && hasVoF[0][1] && hasVoF[1][1];

  return hasAllVoFs;
}
void
EBPatchPolytropic::
pointExtrapToCovered3D(Vector<Real>&           a_extrapVal,
                       const Vector<VolIndex>& a_vofsStencil,
                       const EBCellFAB&        a_primMinu,
                       const EBCellFAB&        a_primPlus,
                       const EBCellFAB&        a_primState,
                       const int&              a_faceDir,
                       const VolIndex&         a_vof,
                       const RealVect&         a_normal,
                       const Side::LoHiSide&   a_sd,
                       const int&              a_numPrim)
{
  RealVect normal = a_normal;
  CH_assert(SpaceDim==3);
  bool dropExtrap = false;
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
  VolIndex vofSten[2][2];

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

  for (int ix = 0; ix < 2; ix++)
    {
      for (int iy = 0; iy < 2; iy++)
        {
          hasVoF[ix][iy] =  EBArith::isVoFHere(vofSten[ix][iy], a_vofsStencil, ivSten[ix][iy]);
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
                  if (a_sd == Side::Hi)
                    {
                      WVal[ix][iy] = a_primMinu(vofSten[ix][iy], ivar);
                    }
                  else
                    {
                      WVal[ix][iy] = a_primPlus(vofSten[ix][iy], ivar);
                    }
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

              if (usesFlattening())
                {
                  int pindex = pressureIndex();
                  int dindex = densityIndex();
                  if ((ivar == dindex) && (a_extrapVal[ivar] < 0.5*funcVal))
                    {
                      dropExtrap = true;
                    }

                  if ((ivar == pindex) || (ivar == dindex))
                    {
                      //if extrapolated to a negative density or pressure then turn
                      //off the slopes of the extrapolation
                      if (a_extrapVal[ivar] < 1.0e-8)
                        {
                          a_extrapVal[ivar] = funcVal;
                        }
                    }
                  if (dropExtrap)
                    {
                      a_extrapVal[ivar] = funcVal;
                    }
                }
            }
          else
            {
              //tanDir[d1] the biggestnormal
              Real x0Val = (m_dx[tangenDir[d1]]*anormNorm   )/(m_dx[    a_faceDir]*anormTan[d1]);
              Real x1Val = (m_dx[tangenDir[d1]]*anormTan[d2])/(m_dx[tangenDir[d2]]*anormTan[d1]);
              Real    funcVal = bilinearFunc(WVal,    x0Val,x1Val);
              Real deltaWNorm = bilinearFunc(delWNorm,x0Val,x1Val);
              Real deltaWTan0 = bilinearFunc(delWTan0,x0Val,x1Val);
              Real deltaWTan1 = bilinearFunc(delWTan1,x0Val,x1Val);

              Real dWdn =  -signNorm*x0Val*deltaWNorm - signTang[d2]*x1Val*deltaWTan1 - signTang[d1]*deltaWTan0;
              a_extrapVal[ivar] = funcVal + dWdn;

              if (usesFlattening())
                {
                  int pindex = pressureIndex();
                  int dindex = densityIndex();
                  if ((ivar == dindex) && (a_extrapVal[ivar] < 0.5*funcVal))
                    {
                      dropExtrap = true;
                    }
                  if ((ivar == pindex) || (ivar == dindex))
                    {
                      //if extrapolated to a negative density or pressure then turn
                      //off the slopes of the extrapolation
                      if (a_extrapVal[ivar] < 0.0)
                        {
                          a_extrapVal[ivar] = funcVal;
                        }
                    }
                  if (dropExtrap)
                    {
                      a_extrapVal[ivar] = funcVal;
                    }
                }
            }
        }
      else
        {
          MayDay::Error("this should not get here because I precomputed hasAllVoFs");
        }
    } //end loop over variables
}
void
EBPatchPolytropic::
pemberExtrap(Vector<Real>&           a_extrapVal,
             const Vector<VolIndex>& a_vofsStencil,
             const EBCellFAB&        a_primMinu,
             const EBCellFAB&        a_primPlus,
             const EBCellFAB&        a_primState,
             const int&              a_faceDir,
             const VolIndex&         a_vof,
             const RealVect&         a_normal,
             const Side::LoHiSide&   a_sd,
             const int&              a_numPrim)
{
  for (int ivar = 0; ivar < a_numPrim; ivar++)
    {
      Real volsum = 0;
      Real valsum = 0;
      for(int ivof = 0; ivof < a_vofsStencil.size(); ivof++)
        {
          Real volfrac = m_ebisBox.volFrac(a_vofsStencil[ivof]);
          volsum += volfrac;
          if(a_sd == Side::Hi)
            {
              valsum += volfrac*a_primMinu(a_vofsStencil[ivof], ivar);
            }
          else
            {
              valsum += volfrac*a_primPlus(a_vofsStencil[ivof], ivar);
            }
        }
      if(volsum > 1.0e-6)
        {
          valsum /= volsum;
        }
      a_extrapVal[ivar] = valsum;
    }
}
/********/
void
EBPatchPolytropic::
pointExtrapToCovered2D(Vector<Real>&           a_extrapVal,
                       const Vector<VolIndex>& a_vofsStencil,
                       const EBCellFAB&        a_primMinu,
                       const EBCellFAB&        a_primPlus,
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
  Real         val[2][2];
  ivSten[0][0] = a_vof.gridIndex();
  ivSten[1][0] = ivSten[0][0] + signNorm*BASISV(a_faceDir);
  ivSten[0][1] = ivSten[0][0] - signNorm*BASISV(a_faceDir) + signTang*BASISV(tangenDir) ;;
  ivSten[1][1] = ivSten[0][0]                              + signTang*BASISV(tangenDir) ;

  for (int ivar = 0; ivar < a_numPrim; ivar++)
    {
      for (int ix = 0; ix < 2; ix++)
        {
          for (int iy = 0; iy < 2; iy++)
            {
              hasVoF[ix][iy] = EBArith::isVoFHere(vof[ix][iy], a_vofsStencil, ivSten[ix][iy]);
              if (hasVoF[ix][iy])
                {
                  if (a_sd == Side::Hi)
                    {
                      val[ix][iy] = a_primMinu(vof[ix][iy], ivar);
                    }
                  else
                    {
                      val[ix][iy] = a_primPlus(vof[ix][iy], ivar);
                    }
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
      Real eps = 1.0e-12;
      if ((Abs(a_normal[a_faceDir]) < eps) || (Abs(a_normal[tangenDir]/a_normal[a_faceDir]) > m_dx[tangenDir]/m_dx[a_faceDir]))
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
          a_extrapVal[ivar] = wc - xBig*signBigNorm*dWBigExtrap - xLit*signLitNorm*dWLitExtrap;
          //debug turn off extrapolation bit
          //a_extrapVal[ivar] = wc;
          //end debug

          if (usesFlattening())
            {
              int pindex = pressureIndex();
              int dindex = densityIndex();

              if ((ivar == pindex) || (ivar == dindex))
                {
                  //if extrapolated to a negative density or pressure then turn
                  //off the slopes of the extrapolation
                  if (a_extrapVal[ivar] < 0.0)
                    {
                      a_extrapVal[ivar] = xLit*wLitNorm + xBig*wBigNorm;
                    }
                }
            }
        }
      else
        {
          MayDay::Error("should not get here since I precalculated these cases");
        }
    }
}
/******/
void
EBPatchPolytropic::
extrapToCoveredFaces(BaseIVFAB<Real>&        a_extendedPrim,
                     const EBCellFAB&        a_primMinu,
                     const EBCellFAB&        a_primPlus,
                     const EBCellFAB&        a_primState,
                     const Vector<VolIndex>& a_coveredFaces,
                     const int&              a_faceDir,
                     const Side::LoHiSide&   a_sd,
                     const Box&              a_box)
{
  CH_TIME("EBPatchPolytropic::extrapToCoveredFaces");
  for (int ivof = 0; ivof < a_coveredFaces.size(); ivof++)
    {
      const VolIndex& vof = a_coveredFaces[ivof];

      RealVect normal = m_ebisBox.normal(vof);
      if (a_box.contains(vof.gridIndex()))
        {
          const int numPrim = numPrimitives();
          Vector<Real> extPrim(numPrim, 0.0);

          Box localBox(vof.gridIndex(), vof.gridIndex());
          localBox.grow(1);
          bool nearBoundary = (!m_domain.contains(localBox));
          bool hasAllVoFs = checkCoveredFaceStencil(vof, m_vofsStencil(vof, 0), normal, a_faceDir, a_sd);
          if((!hasAllVoFs) || nearBoundary)
            {
              pemberExtrap(extPrim,
                           m_vofsStencil(vof,0),
                           a_primMinu,
                           a_primPlus,
                           a_primState,
                           a_faceDir,
                           vof,
                           normal,
                           a_sd,
                           numPrim);
            }
          else
            {
              if (SpaceDim== 2)
                {
                  pointExtrapToCovered2D(extPrim,
                                         m_vofsStencil(vof,0),
                                         a_primMinu,
                                         a_primPlus,
                                         a_primState,
                                         a_faceDir,
                                         vof,
                                         normal,
                                         a_sd,
                                         numPrim);
                }
              else if (SpaceDim==3)
                {
                  pointExtrapToCovered3D(extPrim,
                                         m_vofsStencil(vof,0),
                                         a_primMinu,
                                         a_primPlus,
                                         a_primState,
                                         a_faceDir,
                                         vof,
                                         normal,
                                         a_sd,
                                         numPrim);
                }
              else
                {
                  MayDay::Error("Bogus SpaceDim");
                }
            }
          for (int ivar = 0; ivar < numPrim; ivar++)
            {
              a_extendedPrim(vof, ivar) = extPrim[ivar];
            }
        }
    }

  const IntVectSet& ivs  = a_extendedPrim.getIVS();
  floorPrimitives(a_extendedPrim, ivs);
}
//-----------------------------------------------------------------------
void fillValues(Vector<Real>&            a_prim,
                const Vector<FaceIndex>& a_faces,
                const EBFaceFAB&         a_facePrim,
                const BaseIVFAB<Real>&   a_coveredPrim,
                const VolIndex           a_vof)
{
  CH_assert(a_prim.size() == a_facePrim.nComp());
  CH_assert(a_prim.size() == a_coveredPrim.nComp());

  for(int ivar = 0; ivar < a_prim.size(); ivar++)
    {
      if(a_faces.size() == 0)
        {
          a_prim[ivar] = a_coveredPrim(a_vof, ivar);
        }
      else
        {
          Real primVal = 0;
          for(int iface = 0; iface < a_faces.size(); iface++)
            {
              primVal += a_facePrim(a_faces[iface], ivar);
            }
          if(a_faces.size() > 1) primVal /= Real(a_faces.size());
          a_prim[ivar] = primVal;
        }
    }
}
//-----------------------------------------------------------------------
void
EBPatchPolytropic::
getRhoUDotDele(EBCellFAB&          a_rhoUDotDele,
               EBFluxFAB&          a_facePrim,
               BaseIVFAB<Real>     a_coveredPrimLo[SpaceDim],
               BaseIVFAB<Real>     a_coveredPrimHi[SpaceDim],
               const Box&          a_box,
               const IntVectSet&   a_ivsSmall)
{
  CH_TIME("EBPatchPolytropic::getUDotDele");

  //this is necessary since we are doing increments from here
  a_rhoUDotDele.setVal(0.); 


  for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      FORT_INCRRHOUDOTDELE(CHF_FRA1(a_rhoUDotDele.getSingleValuedFAB(), 0),
                           CHF_CONST_FRA(a_facePrim[faceDir].getSingleValuedFAB()),
                           CHF_CONST_REAL(m_specHeat),
                           CHF_CONST_REAL(m_dx[faceDir]),
                           CHF_CONST_INT(faceDir),
                           CHF_BOX(a_box));
    }

  //update the irregular vofs
  for(int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if(a_box.contains(vof.gridIndex()))
        {
          Real irregRHSVelo = 0;
          for(int diffDir = 0; diffDir < SpaceDim; diffDir++)
            {
              Vector<Real> primLo(QNUM);
              Vector<Real> primHi(QNUM);
              Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, diffDir, Side::Lo);
              Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, diffDir, Side::Hi);
              fillValues(primLo, facesLo, a_facePrim[diffDir], a_coveredPrimLo[diffDir], vof);
              fillValues(primHi, facesHi, a_facePrim[diffDir], a_coveredPrimHi[diffDir], vof);
              Real eLo, eHi;
              FORT_POINTGETTEMPFROMPRIM(CHF_REAL(eLo),CHF_CONST_VR(primLo), CHF_REAL(m_specHeat));
              FORT_POINTGETTEMPFROMPRIM(CHF_REAL(eHi),CHF_CONST_VR(primHi), CHF_REAL(m_specHeat));
              eLo *= m_specHeat;
              eHi *= m_specHeat;

              int dirc = QVELX + diffDir;
              Real rhoave   =   0.5*(primHi[QRHO ] +  primLo[QRHO ]);
              Real ufaceave =   0.5*(primHi[dirc ] +  primLo[dirc ]);

              irregRHSVelo  += rhoave*ufaceave*(eHi - eLo)/m_dx[diffDir];
            }
          a_rhoUDotDele(vof, 0) = irregRHSVelo;
        }
    }//end loop over irreg vofs
}
////////
void
EBPatchPolytropic::
getGradU(EBCellFAB&          a_gradU,
         EBFluxFAB&          a_facePrim,
         BaseIVFAB<Real>     a_coveredPrimLo[SpaceDim],
         BaseIVFAB<Real>     a_coveredPrimHi[SpaceDim],
         const Box&          a_box,
         const IntVectSet&   a_ivsSmall)
{
      //always difference in the facedir direction
  for(int velDir = 0; velDir < SpaceDim; velDir++)
    {
      for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          int gradIndex = TensorCFInterp::gradIndex(velDir, faceDir);

          FORT_EBPPGRADVEL(CHF_FRA(a_gradU.getSingleValuedFAB()),
                           CHF_CONST_FRA(a_facePrim[faceDir].getSingleValuedFAB()),
                           CHF_CONST_REAL(m_dx[faceDir]),
                           CHF_CONST_INT(faceDir),
                           CHF_CONST_INT(velDir),
                           CHF_CONST_INT(gradIndex),
                           CHF_BOX(a_box));
          
          //update the irregular vofs
          for(int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
            {
              const VolIndex& vof = m_irregVoFs[ivof];
              if(a_box.contains(vof.gridIndex()))
                {
                  Vector<Real> primLo(QNUM);
                  Vector<Real> primHi(QNUM);
                  Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, faceDir, Side::Lo);
                  Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, faceDir, Side::Hi);
                  fillValues(primLo, facesLo, a_facePrim[faceDir], a_coveredPrimLo[faceDir], vof);
                  fillValues(primHi, facesHi, a_facePrim[faceDir], a_coveredPrimHi[faceDir], vof);
                      
                  a_gradU(vof, gradIndex) = (primHi[QVELX+velDir] - primLo[QVELX+velDir])/m_dx[faceDir];
                }
            }//end loop over irreg vofs
        }
    }

}

////////
void
EBPatchPolytropic::
getVeloHalf(EBCellFAB&          a_veloHalf,
            EBFluxFAB&          a_facePrim,
            BaseIVFAB<Real>     a_coveredPrimLo[SpaceDim],
            BaseIVFAB<Real>     a_coveredPrimHi[SpaceDim],
            const Box&          a_box,
            const IntVectSet&   a_ivsSmall)
{
  //always average the normal (facedir) velocities
  for(int velDir = 0; velDir < SpaceDim; velDir++)
    {

      FORT_EBPPVELOHALF(CHF_FRA(a_veloHalf.getSingleValuedFAB()),
                        CHF_CONST_FRA(a_facePrim[velDir].getSingleValuedFAB()),
                        CHF_CONST_INT(velDir),
                        CHF_BOX(a_box));
          
      //update the irregular vofs
      for(int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
        {
          const VolIndex& vof = m_irregVoFs[ivof];
          if(a_box.contains(vof.gridIndex()))
            {
              Vector<Real> primLo(QNUM);
              Vector<Real> primHi(QNUM);
              Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, velDir, Side::Lo);
              Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, velDir, Side::Hi);
              fillValues(primLo, facesLo, a_facePrim[velDir], a_coveredPrimLo[velDir], vof);
              fillValues(primHi, facesHi, a_facePrim[velDir], a_coveredPrimHi[velDir], vof);
                      
              a_veloHalf(vof, velDir) = 0.5*(primHi[QVELX+velDir] + primLo[QVELX+velDir]);
            }
        }//end loop over irreg vofs
    }

}
////////
void
EBPatchPolytropic::
getPDivU(EBCellFAB&          a_pDivU,
         EBFluxFAB&          a_facePrim,
         BaseIVFAB<Real>     a_coveredPrimLo[SpaceDim],
         BaseIVFAB<Real>     a_coveredPrimHi[SpaceDim],
         const Box&          a_box,
         const IntVectSet&   a_ivsSmall)
{
  //have to do this because we are incrementing
  a_pDivU.setVal(0.);
  for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      FORT_INCRPDIVU(CHF_FRA1(a_pDivU.getSingleValuedFAB(), 0),
                     CHF_CONST_FRA(a_facePrim[faceDir].getSingleValuedFAB()),
                     CHF_CONST_REAL(m_dx[faceDir]),
                     CHF_CONST_INT(faceDir),
                     CHF_BOX(a_box));
    }
          
  //update the irregular vofs
  for(int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if(a_box.contains(vof.gridIndex()))
        {
          Real irregVal = 0;
          for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
            {
              int vind = QVELX + faceDir;

              Vector<Real> primLo(QNUM);
              Vector<Real> primHi(QNUM);
              Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, faceDir, Side::Lo);
              Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, faceDir, Side::Hi);
              fillValues(primLo, facesLo, a_facePrim[faceDir], a_coveredPrimLo[faceDir], vof);
              fillValues(primHi, facesHi, a_facePrim[faceDir], a_coveredPrimHi[faceDir], vof);
                  
              Real pave  = 0.5*(primLo[QPRES] + primHi[QPRES]);
              Real udiff =     (primLo[vind ] + primHi[vind ])/m_dx[faceDir];
              irregVal = irregVal + pave*udiff;
            }
          a_pDivU(vof, 0) = irregVal;
        }
    }//end loop over irreg vofs
}


//-----------------------------------------------------------------------
void
EBPatchPolytropic::
primitivesAndDivergences(EBCellFAB&          a_nonConsDivF,
                         EBCellFAB&          a_consState,
                         EBFluxFAB&          a_facePrim,
                         BaseIVFAB<Real>     a_coveredPrimMinu[SpaceDim],
                         BaseIVFAB<Real>     a_coveredPrimPlus[SpaceDim],
                         Vector<VolIndex>    a_coveredFaceMinu[SpaceDim],
                         Vector<VolIndex>    a_coveredFacePlus[SpaceDim],
                         EBFluxFAB&          a_flux,
                         BaseIVFAB<Real>&    a_ebIrregFlux,
                         BaseIVFAB<Real>&    a_nonConservativeDivergence,
                         const EBCellFAB&    a_flattening,
                         const EBCellFAB&    a_source,
                         const Box&          a_box,
                         const IntVectSet&   a_ivsSmall,
                         const DataIndex&    a_dit,
                         bool                a_verbose)
{
  CH_TIME("EBPatchPolytropic::regularUpdate");
  CH_assert(isDefined());
  int numCons = numConserved();

  CH_assert(a_flux.nComp() >= numConserved());
  CH_assert(a_consState.nComp() >= numConserved());
  setCoveredConsVals(a_consState);

  EBCellFAB slopePrim[SpaceDim];
  EBCellFAB slopeNLim[SpaceDim];

  computeFluxes(a_flux,
                m_coveredFluxMinuG4, m_coveredFluxPlusG4,
                a_coveredFaceMinu,   a_coveredFacePlus,
                m_primState,    slopePrim, slopeNLim,
                a_flattening, a_consState, a_source,
                a_box, a_dit,a_verbose);

  //now I happen to know that m_primgdnv holds faceprim
  //and m_extendedStateMinuG4 holds  coveredPrimMinu and so on
  a_facePrim.define(m_ebisBox, a_box, numPrimitives());
  Interval inter(0, numPrimitives()-1);
  a_facePrim.copy(a_box, inter, a_box, m_primGdnv, inter);

  for(int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      IntVectSet ivsMinu = m_coveredSetsMinuG4[faceDir] & a_box;
      IntVectSet ivsPlus = m_coveredSetsPlusG4[faceDir] & a_box;
      a_coveredPrimMinu[faceDir].define(ivsMinu, m_ebisBox.getEBGraph(), numPrimitives());
      a_coveredPrimPlus[faceDir].define(ivsPlus, m_ebisBox.getEBGraph(), numPrimitives());
      a_coveredPrimMinu[faceDir].copy(a_box, inter, a_box, m_extendStateMinuG4[faceDir], inter);
      a_coveredPrimPlus[faceDir].copy(a_box, inter, a_box, m_extendStatePlusG4[faceDir], inter);
    }

  //now that we have the fluxes, modify them with
  //artificial viscosity if that is called for.
  //the artificial viscosity correction to the
  //embedded boundary flux is computed inside
  //computeebirregflux
  if(usesArtificialViscosity())
    {
      EBFluxFAB openDivU(m_ebisBox, a_box, 1);

      getFaceDivergence(openDivU,
                        m_primState, slopeNLim,
                        a_box, a_ivsSmall);

      applyArtificialViscosity(a_flux,
                               m_coveredFluxMinuG4,
                               m_coveredFluxPlusG4,
                               m_coveredFaceMinuG4,
                               m_coveredFacePlusG4,
                               a_consState,
                               openDivU,
                               a_box,
                               a_ivsSmall);
    }

  //this is the stable, non-conservative estimate of the solution update
  nonconservativeDivergence(a_nonConsDivF, a_flux,
                            m_coveredFluxMinuG4,
                            m_coveredFluxPlusG4,
                            a_coveredFaceMinu,
                            a_coveredFacePlus,
                            a_box);

  //copy nonconservative div f into sparse output thingy
  IntVectSet ivsIrreg = a_ivsSmall;
  for(VoFIterator vofit(ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for(int ivar = 0; ivar < numCons; ivar++)
        {
          a_nonConservativeDivergence(vof, ivar) = a_nonConsDivF(vof, ivar);
        }
    }

  //compute irregular boundary flux.  this includes
  //an artificial viscosity modification if appropriate
  computeEBIrregFlux( a_ebIrregFlux, m_primState,
                      slopeNLim, ivsIrreg, a_source);
}
/*****************************/
void
EBPatchPolytropic::
doRZCoords(bool a_doRZCoords)
{
  m_doRZCoords = a_doRZCoords;
}
/*****************************/
void
EBPatchPolytropic::
setSource(EBCellFAB&       a_source,
          const EBCellFAB& a_consState,
          const Box&       a_box)
{
}

/*****************************/
void
EBPatchPolytropic::getCoveredValuesCons(Vector<Real>& a_covValues)
{
  a_covValues.resize(numConserved());
  a_covValues[CRHO] = 1.0;
  a_covValues[CENG] = 1.0;
  a_covValues[CMOMX] = 0.0;
  a_covValues[CMOMY] = 0.0;
#if CH_SPACEDIM==3
  a_covValues[CMOMZ] = 0.0;
#endif
}
/*****************************/
void
EBPatchPolytropic::getCoveredValuesPrim(Vector<Real>& a_covValues)
{
  a_covValues.resize(numPrimitives());
  a_covValues[QRHO] = 1.0;
  a_covValues[QVELX] = 0.0;
  a_covValues[QVELY] = 0.0;
  a_covValues[QPRES] = 1.0;
  a_covValues[QENTR] = 0.0;
  a_covValues[QINTERN] = 0.0;
  a_covValues[QCVTEMP] = 0.0;
  a_covValues[QC] = 1.0;
#if CH_SPACEDIM==3
  a_covValues[QVELZ] = 0.0;
#endif
}
/*****************************/
void
EBPatchPolytropic::
setGamma(const Real& a_gamma)
{
  m_gamma = a_gamma;
  m_isGammaSet = true;
}
/*****************************/
Real
EBPatchPolytropic::
getGamma() const
{
  CH_assert(m_isGammaSet);
  return m_gamma;
}
/******/
void
EBPatchPolytropic::
setSpecHeat(const Real& a_specHeat)
{
  m_specHeat = a_specHeat;
  m_isSpecHeatSet = true;
}
/*****************************/
Real
EBPatchPolytropic::
getSpecHeat() const
{
  CH_assert(m_isSpecHeatSet);
  return m_specHeat;
}
/******/
bool
EBPatchPolytropic::
usesFlattening() const
{
  CH_assert(m_isSlopeSet);

  return m_useFlattening;
}
/******/
bool
EBPatchPolytropic::
usesArtificialViscosity() const
{
  return m_useArtificialVisc;
}
/******/
bool
EBPatchPolytropic::
usesFourthOrderSlopes() const
{
  CH_assert(m_isSlopeSet);
  return m_useFourthOrderSlopes;
}
/******/
Vector<string>
EBPatchPolytropic::
stateNames()
{
  Vector<string> retval;

  retval.push_back("mass-density");
  retval.push_back("x-momentum");
  if (SpaceDim >= 2)
    {
      retval.push_back("y-momentum");
    }

  if (SpaceDim >= 3)
    {
      retval.push_back("z-momentum");
    }
  retval.push_back("energy-density");
  return retval;
}
/******/
Vector<string>
EBPatchPolytropic::
primNames()
{
  Vector<string> retval;

  ParmParse pp;
  int logflag = 0;
  if(pp.contains("logflag"))
    {
      pp.get("logflag", logflag);
    }
  if(logflag == 1)
    {
      retval.push_back("log10density");
    }
  else
    {
      retval.push_back("density");
    }
  retval.push_back("x-velocity");
  retval.push_back("y-velocity");


#if CH_SPACEDIM==3
  retval.push_back("z-velocity");
#endif

  if(logflag == 1)
    {
      retval.push_back("log10pressure");
      retval.push_back("log10entropy");
    }
  else
    {
      retval.push_back("pressure");
      retval.push_back("entropy");
    }
  retval.push_back("internal_energy");
  retval.push_back("cv_temperature");
  retval.push_back("soundspeed");


  return retval;
}

/******/
Vector<string>
EBPatchPolytropic::
primNamesNoLog()
{
  Vector<string> retval;

  retval.push_back("density");
  retval.push_back("x-velocity");
  retval.push_back("y-velocity");

#if CH_SPACEDIM==3
  retval.push_back("z-velocity");
#endif

  retval.push_back("pressure");
  retval.push_back("entropy");
  retval.push_back("internal_energy");
  retval.push_back("cv_temperature");
  retval.push_back("soundspeed");


  return retval;
}
/******/
int
EBPatchPolytropic::
numPrimitives() const
{
  return QNUM;
}
/******/
int
EBPatchPolytropic::
numFluxes() const
{
  return FNUM;
}
/******/
int
EBPatchPolytropic::
numConserved() const
{
  return CNUM;
}
/******/
int
EBPatchPolytropic::
numSlopes() const
{
  return QSLOPE;
}
/******/
Interval
EBPatchPolytropic::
velocityInterval() const
{
#if CH_SPACEDIM==2
  Interval retval(QVELX, QVELY);
#elif CH_SPACEDIM==3
  Interval retval(QVELX, QVELZ);
#else
  bogus_spacedim();
#endif
  return retval;
}
/******/
int
EBPatchPolytropic::
pressureIndex() const
{
  return QPRES;
}
/******/
int
EBPatchPolytropic::
densityIndex() const
{
  return QRHO;
}
/******/
int
EBPatchPolytropic::
bulkModulusIndex() const
{
  //phil said this was OK 1-2-2002
  //  MayDay::Warning("returning pressure for bulk modulus");
  return QPRES;
}
/******/
Real
EBPatchPolytropic::
artificialViscosityCoefficient() const
{
  ParmParse pp;
  Real retval;
  pp.get("artificial_viscosity", retval);
  return retval;
}
/******/
Real
EBPatchPolytropic::
getMaxWaveSpeed(const EBCellFAB& a_consState,
                const Box& a_box)
{
  CH_TIME("EBPatchPolytropic::getMaxWaveSpeeed");
  CH_assert(m_isDefined && m_isBoxSet && m_isGammaSet);
  Real speed = 0.0;
  const EBISBox& ebisBox = a_consState.getEBISBox();
  IntVectSet ivs(a_box);
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  //cannot just send to fortran because covered cell values
  //would get in and multiply-valued cells would use bogus
  //values that could also get in.  The could be avoided by
  //either masks or some clever manipulation of the underlying
  //basefab
  for(VoFIterator vofit(ivs, ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      //      const  IntVect& iv = vof.gridIndex();
      Real dense = Max(smallr, a_consState(vof, CRHO));
      Real eng   = Max(small,  a_consState(vof, CENG));
      Real vmax = 0.0;
      Real kinetic = 0.0;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          Real momen = a_consState(vof, CMOMX + idir);
          Real vel = Abs(momen/dense);
          kinetic += 0.5*vel*vel;
          vmax = Max(vmax, vel);
        }
      Real intEng = eng/dense - kinetic;
      intEng = Max(intEng, (Real)0.001*small);
      Real sqrtarg = m_gamma*(m_gamma -1.0)*intEng;
      CH_assert(sqrtarg >= 0.0);
      Real cspeed = sqrt(sqrtarg);
      if(Abs(cspeed + vmax) > speed)
        {
          speed = cspeed+ vmax;
          if(speed > EBPatchGodunov::getMaxWaveSpeed())
            {
              EBPatchGodunov::setMaxWaveSpeed(speed);
              EBPatchGodunov::setMaxWaveSpeedIV(vof.gridIndex());
            }
        }
    }
  return speed;

}
/******/
void
EBPatchPolytropic::
consToPrim(EBCellFAB&       a_primState,
           const EBCellFAB& a_consState,
           const Box&       a_box,
           int              a_logflag,
           bool             a_verbose)
{
  CH_TIME("EBPatchPolytropic::consToPrim");
  CH_assert(m_isDefined && m_isBoxSet);
  CH_assert(m_isGammaSet);
  CH_assert(a_consState.getRegion().contains(a_box));
  //have to set the covered cell vals so need the cast
  //does not change real data
  //  EBCellFAB& conData = (EBCellFAB&) a_consState;
  //  setCoveredConsVals(conData);
  //debug
  if(a_verbose)
    {
      pout()  << "constoprim " << endl;
    }
  //end debug
  const BaseFab<Real>&   regCons = a_consState.getSingleValuedFAB();
  BaseFab<Real>&   regPrim = a_primState.getSingleValuedFAB();

  int iverbose = 0;
  if(a_verbose) iverbose = 1;
  FORT_CONS2PRM(CHF_BOX(a_box),
                CHF_CONST_FRA(regCons),
                CHF_FRA(regPrim),
                CHF_CONST_INT(a_logflag),
                CHF_CONST_INT(iverbose));


  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for(VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Vector<Real> conserved(CNUM);
      Vector<Real> primitive(QNUM);
      for(int ivar = 0; ivar < CNUM; ivar++)
        {
          conserved[ivar] = a_consState(vof, ivar);
        }

      FORT_POINTCONS2PRM(CHF_VR(conserved),
                         CHF_VR(primitive),
                         CHF_CONST_INT(a_logflag));

      for(int ivar = 0; ivar < QNUM; ivar++)
        {
          a_primState(vof, ivar)  = primitive[ivar];
        }
      //  setCoveredPrimVals(a_primState);
    }
}
/******/
void
EBPatchPolytropic::
consToPrim(BaseIVFAB<Real>&  a_primState,
           const EBCellFAB&  a_consState,
           const IntVectSet& a_ivs)
{
  CH_TIME("EBPatchPolytropic::consToPrimIrrReg");
  CH_assert(isDefined());
  for(VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Vector<Real> conserved(CNUM);
      Vector<Real> primitive(QNUM);
      for(int ivar = 0; ivar < CNUM; ivar++)
        {
          conserved[ivar] = a_consState(vof, ivar);
        }

      int logflag = 0;
      FORT_POINTCONS2PRM(CHF_VR(conserved),
                         CHF_VR(primitive),
                         CHF_CONST_INT(logflag));

      for(int ivar = 0; ivar < QNUM; ivar++)
        {
          a_primState(vof, ivar)  = primitive[ivar];
        }
    }
}
#ifdef CH_USE_HDF5

void EBPatchPolytropic::expressions(HDF5HeaderData& a_expressions)
{
  char gammaStr[1024], specHeatStr[1024];
  snprintf(gammaStr, 1024, "%g", m_gamma);
  snprintf(specHeatStr, 1024, "%g", m_specHeat);
  a_expressions.m_string["scalar gamma"] = static_cast<const char*>(gammaStr);
  a_expressions.m_string["scalar specHeat"] = static_cast<const char*>(specHeatStr);
  a_expressions.m_string["vector velocity"] = "momentum/<mass-density>";
  a_expressions.m_string["scalar kinetic_energy"] = "dot(velocity,velocity)/2";
  a_expressions.m_string["scalar pressure"] = "(gamma-1)*(<energy-density>-kinetic_energy*<mass-density>";
  a_expressions.m_string["scalar soundspeed"] = "sqrt(gamma*(pressure/<mass-density>))";
  a_expressions.m_string["scalar log10entropy"] = "log10(pressure) - gamma*log10(<mass-density>)";
  a_expressions.m_string["scalar temperature"] = "(<energy-density>-kinetic_energy)/(<mass-density>*specHeat)";
}

#endif

/******/
void
EBPatchPolytropic::
consToPrim(BaseIVFAB<Real>&  a_primState,
           const BaseIVFAB<Real>&  a_consState,
           const IntVectSet& a_ivs)
{
  CH_TIME("EBPatchPolytropic::consToPrimIrr2");
  CH_assert(isDefined());
  for(VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      Vector<Real> conserved(CNUM);
      Vector<Real> primitive(QNUM);
      for(int ivar = 0; ivar < CNUM; ivar++)
        {
          conserved[ivar] = a_consState(vof, ivar);
        }
      int logflag = 0;
      FORT_POINTCONS2PRM(CHF_VR(conserved),
                         CHF_VR(primitive),
                         CHF_CONST_INT(logflag));

      for(int ivar = 0; ivar < QNUM; ivar++)
        {
          a_primState(vof, ivar)  = primitive[ivar];
        }
    }
}
/******/
void
EBPatchPolytropic::
primToCons(EBCellFAB&       a_consState,
           const EBCellFAB& a_primState,
           const Box&       a_box)
{
  CH_TIME("EBPatchPolytropic::primToCons");
  CH_assert(isDefined());
  const Box& regBox = a_box;
  CH_assert(a_primState.getRegion().contains(regBox));
  //set covered vals. does not change real data.
  //  EBCellFAB& primData = (EBCellFAB&) a_primState;
  //  setCoveredPrimVals(primData);

  BaseFab<Real>&   regCons = a_consState.getSingleValuedFAB();
  const BaseFab<Real>&  regPrim = a_primState.getSingleValuedFAB();


  FORT_PRM2CONS(CHF_BOX(regBox),
                CHF_FRA(regCons),
                CHF_CONST_FRA(regPrim));


  IntVectSet ivsMulti = m_ebisBox.getMultiCells(regBox);
  for(VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Vector<Real> conserved(CNUM);
      Vector<Real> primitive(QNUM);
      for(int ivar = 0; ivar < QNUM; ivar++)
        {
          primitive[ivar] = a_primState(vof, ivar);
        }

      FORT_POINTPRM2CONS(CHF_VR(conserved),
                         CHF_VR(primitive));

      for(int ivar = 0; ivar < CNUM; ivar++)
        {
          a_consState(vof, ivar)  = conserved[ivar];
        }
    }
  //  setCoveredConsVals(a_consState);
}
/******/
void
EBPatchPolytropic::
primToCons(BaseIVFAB<Real>&       a_consState,
           const BaseIVFAB<Real>& a_primState,
           const IntVectSet&      a_ivs)
{
  CH_TIME("EBPatchPolytropic::primToConsIrr");
  CH_assert(isDefined());
  for(VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Vector<Real> conserved(CNUM);
      Vector<Real> primitive(QNUM);
      for(int ivar = 0; ivar < QNUM; ivar++)
        {
          primitive[ivar] = a_primState(vof, ivar);
        }

      FORT_POINTPRM2CONS(CHF_VR(conserved),
                         CHF_VR(primitive));

      for(int ivar = 0; ivar < CNUM; ivar++)
        {
          a_consState(vof, ivar)  = conserved[ivar];
        }
    }
}
/*****************************/
void EBPatchPolytropic::
normalPred(EBCellFAB&       a_primLo,
           EBCellFAB&       a_primHi,
           const EBCellFAB& a_primState,
           const EBCellFAB& a_slopePrim,
           const Real&      a_dtbydx,
           const int&       a_dir,
           const Box&       a_box)
{
  CH_TIME("EBPatchPolytropic::normalPred");
  int nslope =  numSlopes();
  CH_assert(isDefined());
  CH_assert(a_primLo.nComp()    >= nslope);
  CH_assert(a_primHi.nComp()    >= nslope);
  CH_assert(a_primState.nComp() >= nslope);
  CH_assert(a_slopePrim.nComp() == nslope);

  Real dtbydx = a_dtbydx;
  const BaseFab<Real>& regState = a_primState.getSingleValuedFAB();
  const BaseFab<Real>& regSlope = a_slopePrim.getSingleValuedFAB();
  BaseFab<Real>& regPrimLo = a_primLo.getSingleValuedFAB();
  BaseFab<Real>& regPrimHi = a_primHi.getSingleValuedFAB();
  int useflat = 0;
  if(usesFlattening())
    useflat = 1;


  FORT_PRED(CHF_BOX(a_box),
            CHF_CONST_FRA(regState),
            CHF_CONST_FRA(regSlope),
            CHF_FRA(regPrimLo),
            CHF_FRA(regPrimHi),
            CHF_CONST_INT(a_dir),
            CHF_CONST_REAL(dtbydx),
            CHF_CONST_INT(useflat));

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for(VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      Vector<Real> primlo(QNUM), primhi(QNUM), primit(QNUM), pslope(QNUM);
      for(int ivar = 0; ivar < QNUM; ivar++)
        {
          primit[ivar] = a_primState(vof, ivar);
          pslope[ivar] = a_slopePrim(vof, ivar);
        }


      FORT_POINTPRED(CHF_VR(primit),
                     CHF_VR(pslope),
                     CHF_VR(primlo),
                     CHF_VR(primhi),
                     CHF_CONST_INT(a_dir),
                     CHF_CONST_REAL(dtbydx),
                     CHF_CONST_INT(useflat));


      for(int ivar = 0; ivar < QNUM; ivar++)
        {
          a_primLo(vof, ivar) = primlo[ivar];
          a_primHi(vof, ivar) = primhi[ivar];
        }
    }
}
/******/
void
EBPatchPolytropic::
riemann(EBFaceFAB&       a_flux,
        const EBCellFAB& a_primLeft,
        const EBCellFAB& a_primRigh,
        const int&       a_dir,
        const Box&       a_box)
{
  CH_TIME("EBPatchPolytropic::riemann");
  //int nPrim = numPrimitives();
  Box cellBox = enclosedCells(a_box);
  EBFaceFAB& primGdnv = m_primGdnv[a_dir];

  //explicitly cast the left and right states to modifiable references
  //because we need to shift them to faces and then shift them back.
  BaseFab<Real>& regPrimRigh = (BaseFab<Real>&)a_primRigh.getSingleValuedFAB();
  BaseFab<Real>& regPrimLeft = (BaseFab<Real>&)a_primLeft.getSingleValuedFAB();
  BaseFab<Real>& regPrimGdnv = primGdnv.getSingleValuedFAB();
  //shift the cell centered stuff to edges
  regPrimRigh.shiftHalf(a_dir, -1);
  regPrimLeft.shiftHalf(a_dir,  1);

  CH_assert(regPrimRigh.box().contains(a_box));
  CH_assert(regPrimLeft.box().contains(a_box));
  CH_assert(regPrimGdnv.box().contains(a_box));

  //find the regular part of qgdnv.

  FORT_RIEMANN(CHF_BOX(a_box),
               CHF_CONST_FRA(regPrimLeft),
               CHF_CONST_FRA(regPrimRigh),
               CHF_FRA(regPrimGdnv),
               CHF_CONST_INT(a_dir));


  //shift the cell centered stuff back to cells
  regPrimRigh.shiftHalf(a_dir,  1);
  regPrimLeft.shiftHalf(a_dir, -1);

  //the box sent into this is face-centered.
  //we need to use the cell-centered one it surrounds.
  //this can be more than multivalued cells if there
  //are neighbors that are multivalued
  
  Box grownBox = cellBox;
  grownBox.grow(a_dir, 1);
  grownBox &= m_domain;
  IntVectSet ivsMulti = m_ebisBox.getIrregIVS(grownBox);

  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingNoBoundary;
  for(FaceIterator faceit(ivsMulti, m_ebisBox.getEBGraph(), a_dir, stopCrit);
      faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      if(a_box.contains(face.gridIndex(Side::Hi)))
        {
          VolIndex vofl = face.getVoF(Side::Lo);
          VolIndex vofr = face.getVoF(Side::Hi);
          Vector<Real> ql(QNUM), qr(QNUM), qgod(QNUM);
          for(int ivar = 0; ivar < QNUM; ivar++)
            {
              ql[ivar] = a_primLeft(vofl, ivar);
              qr[ivar] = a_primRigh(vofr, ivar);
            }


          FORT_POINTRIEMANN(CHF_VR(ql), CHF_VR(qr), CHF_VR(qgod),
                            CHF_CONST_INT(a_dir));


          for(int ivar = 0; ivar < QNUM; ivar++)
            {
              primGdnv(face, ivar) = qgod[ivar];
            }
        }
    }
  //from the godunov state, get the flux
  //fix the prim gdnv boundary face.
  for(SideIterator sit; sit.ok(); ++sit)
    {
      Box boundBox; 
      if(sit() == Side::Lo)
        {
          boundBox = adjCellLo(m_domain, a_dir, 1);
        }
      else
        {
          boundBox = adjCellHi(m_domain, a_dir, 1);
        }
      boundBox.shift(a_dir, -sign(sit()));
      IntVectSet ivsBound(boundBox);
      ivsBound &= m_validBox;
      if(!ivsBound.isEmpty())
        {
          FaceStop::WhichFaces stopCrit = FaceStop::AllBoundaryOnly;
          for(FaceIterator faceit(ivsBound, m_ebisBox.getEBGraph(), a_dir, stopCrit);
              faceit.ok(); ++faceit)
            {
              for(int ivar = 0; ivar < QNUM; ivar++)
                {
                  Real primExtrap;
                  if(sit() == Side::Lo)
                    {
                      primExtrap = a_primRigh(faceit().getVoF(Side::Hi), ivar);
                    }
                  else
                    {
                      primExtrap = a_primLeft(faceit().getVoF(Side::Lo), ivar);
                    }
                  primGdnv(faceit(), ivar) = primExtrap;
                }
            }
        }
    }
  getFlux(a_flux, primGdnv, a_dir, a_box);
}
/*****************************/
void
EBPatchPolytropic::
riemann(BaseIVFAB<Real>&        a_coveredFlux,
        const BaseIVFAB<Real>&  a_exteState,
        const EBCellFAB&        a_primState,
        const Vector<VolIndex>& a_vofset,
        const int&              a_dir,
        const Side::LoHiSide&   a_sd,
        const Box&       a_box)
{
  CH_TIME("EBPatchPolytropic::riemannIrr");
  for(int ivof = 0; ivof < a_vofset.size(); ivof++)
    {
      const VolIndex& vof = a_vofset[ivof];
      if(a_box.contains(vof.gridIndex()))
        {
          Vector<Real> ql(QNUM), qr(QNUM), qgod(QNUM), flux(FNUM);
          int iloc = vof.gridIndex()[0];
          Real radius;

          if(a_sd == Side::Hi)
            {
              radius = (iloc+1)*m_dx[0];
              for(int ivar = 0; ivar < QNUM; ivar++)
                {
                  ql[ivar] = a_primState(vof, ivar);
                  qr[ivar] = a_exteState(vof, ivar);
                }
            }
          else
            {
              radius = iloc*m_dx[0];
              for(int ivar = 0; ivar < QNUM; ivar++)
                {
                  ql[ivar] = a_exteState(vof, ivar);
                  qr[ivar] = a_primState(vof, ivar);
                }
            }


          FORT_POINTRIEMANN(CHF_VR(ql), CHF_VR(qr), CHF_VR(qgod),
                            CHF_CONST_INT(a_dir));

          //only takes effect if idir == 0
          FORT_POINTGETFLUX(CHF_VR(flux), CHF_VR(qgod),
                            CHF_CONST_INT(a_dir));


          for(int ivar = 0; ivar < FNUM; ivar++)
            {
              a_coveredFlux(vof, ivar) = flux[ivar];
            }
        }
    }
}
/******/
void
EBPatchPolytropic::
getFlux(EBFaceFAB&       a_flux,
        const EBFaceFAB& a_prim,
        const int&       a_dir,
        const Box&       a_box)
{
  CH_TIME("EBPatchPolytropic::getFlux");
  CH_assert(a_flux.direction() == a_dir);
  CH_assert(a_flux.getRegion().contains(a_box));
  CH_assert(a_prim.getRegion().contains(a_box));
  CH_assert(a_prim.direction() == a_dir);

  const BaseFab<Real>& regPrim = a_prim.getSingleValuedFAB();
  BaseFab<Real>& regFlux = a_flux.getSingleValuedFAB();

  FORT_GETFLUX(CHF_BOX(a_box),
               CHF_CONST_FRA(regPrim),
               CHF_CONST_INT(a_dir),
               CHF_FRA(regFlux),
               CHF_CONST_REAL(m_dx[0]));

  //box sent in is  face-centered
  //this grow dance fixes the problem of stuff being next to multivalued
  //cells
  Box cellBox = enclosedCells(a_box);
  cellBox.grow(a_dir, 1);
  cellBox &= m_domain;

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(cellBox);
  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingNoBoundary;
  for(FaceIterator faceit(ivsMulti, m_ebisBox.getEBGraph(), a_dir, stopCrit);
      faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      if(a_box.contains(face.gridIndex(Side::Hi)))
        {
          Vector<Real> qgdnv(QNUM), fluxvec(FNUM);
          for(int ivar = 0; ivar < QNUM; ivar++)
            {
              qgdnv[ivar] = a_prim(face, ivar);
            }

          FORT_POINTGETFLUX(CHF_VR(fluxvec),
                            CHF_VR(qgdnv),
                            CHF_CONST_INT(a_dir));

          for(int ivar = 0; ivar < FNUM; ivar++)
            {
              a_flux(face, ivar) = fluxvec[ivar];
            }
        }
    }
}
/******/
void
EBPatchPolytropic::
updateCons(EBCellFAB&            a_consState,
           const EBFaceFAB&        a_flux,
           const BaseIVFAB<Real>&  a_coveredFluxMinu,
           const BaseIVFAB<Real>&  a_coveredFluxPlus,
           const Vector<VolIndex>& a_coveredFaceMinu,
           const Vector<VolIndex>& a_coveredFacePlus,
           const int&              a_dir,
           const Box&              a_box,
           const Real&             a_scale)
{
  CH_TIME("EBPatchPolytropic::updateCons");
  EBPatchGodunov::updateCons( a_consState,
                              a_flux,
                              a_coveredFluxMinu,
                              a_coveredFluxPlus,
                              a_coveredFaceMinu,
                              a_coveredFacePlus,
                              a_dir,
                              a_box,
                              a_scale);
}

/******/
void
EBPatchPolytropic::floorPrimitives(EBCellFAB&  a_primState,
                                   const Box&  a_box)
{
  CH_TIME("EBPatchPolytropic::floorPrimitives");
  BaseFab<Real>&       primReg = a_primState.getSingleValuedFAB();

  FORT_FLOORPRIM( CHF_BOX(a_box),
                  CHF_FRA(primReg));
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  for(int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if(a_box.contains(vof.gridIndex()))
        {
          a_primState(vof, QRHO)  = Max(a_primState(vof, QRHO) ,  smallr);
          a_primState(vof, QPRES) = Max(a_primState(vof, QPRES),  smallp);
        }
    }
}
/******/
void
EBPatchPolytropic::floorConserved(EBCellFAB&  a_consState,
                                  const Box&  a_box)
{
  CH_TIME("EBPatchPolytropic::floorConserved");
  BaseFab<Real>&       consReg = a_consState.getSingleValuedFAB();

  FORT_FLOORCONS( CHF_BOX(a_box),
                  CHF_FRA(consReg));

  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  IntVectSet ivsIrreg = m_ebisBox.getIrregIVS(a_box);
  for(VoFIterator vofit(ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      a_consState(vof, CRHO) = Max(a_consState(vof, CRHO), smallr);
      a_consState(vof, CENG) = Max(a_consState(vof, CENG), small );
      //      Real rmom =       a_consState(vof, CMOMX) ;
      //      Real zmom =       a_consState(vof, CMOMY) ;
      //      Real ratio = Abs(rmom/zmom);
      //      if((ratio > 1.0) && (Abs(zmom) > 1.0e-3))
      //        {
      //          cout << vof << " has wacked ratio =" << ratio << endl;
      //        }
    }
}

/******/
void
EBPatchPolytropic::floorTemperature(EBCellFAB&  a_temperat,
                                    const Box&  a_box)
{
  CH_TIME("EBPatchPolytropic::floorConserved");
  FORT_FLOORTEMPERATURE( CHF_BOX(a_box),
                         CHF_FRA1(a_temperat.getSingleValuedFAB(), 0));

  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  IntVectSet ivsIrreg = m_ebisBox.getIrregIVS(a_box);
  for(VoFIterator vofit(ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      a_temperat(vof, 0) = Max(a_temperat(vof, 0), small );
    }
}
/******/
void
EBPatchPolytropic::floorConserved(BaseIVFAB<Real>& a_consState,
                                  const IntVectSet& a_ivsIrreg)
{
  CH_TIME("EBPatchPolytropic::floorConservedIrr");
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  for(VoFIterator vofit(a_ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      a_consState(vof, CRHO) = Max(a_consState(vof, CRHO), smallr);
      a_consState(vof, CENG) = Max(a_consState(vof, CENG), small );
    }

}
/******/
void
EBPatchPolytropic::floorPrimitives(BaseIVFAB<Real>& a_primState,
                                   const IntVectSet& a_ivsIrreg)
{
  CH_TIME("EBPatchPolytropic::floorPrimitivesIrr");
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  for(VoFIterator vofit(a_ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      //floors
      a_primState(vof, QRHO)  = Max(a_primState(vof, QRHO) ,  smallr);
      a_primState(vof, QPRES) = Max(a_primState(vof, QPRES),  smallp);
    }
}
/*****************************/
EBPatchPolytropic::
EBPatchPolytropic():EBPatchGodunov()
{
  m_isGammaSet    = false;
  m_isSpecHeatSet = false;
  m_isArtViscSet  = false;
  m_doRZCoords    = false;
}
/******/
bool
EBPatchPolytropic::
isDefined() const
{
  bool retval = EBPatchGodunov::isDefined() &&  m_isGammaSet;
  return retval;
}

/*****************************/
EBPatchPolytropic::
~EBPatchPolytropic()
{
}
/*****************************/
/*****************************/
void
EBPatchPolytropic::
computeEBIrregFlux(BaseIVFAB<Real>&  a_ebIrregFlux,
                   const EBCellFAB&  a_primState,
                   const EBCellFAB   a_slopePrim[SpaceDim],
                   const IntVectSet& a_irregIVS,
                   const EBCellFAB&  a_source)
{
  CH_TIME("EBPatchPolytropic::computeEBIrregFlux");
  EBCellFAB* primPtr = (EBCellFAB*)(&a_primState);
  EBCellFAB primT;
  if(a_source.isDefined() && (s_conservativeSource))
    {
      CH_assert(a_source.box().contains(m_validBox));
      primT.define(m_ebisBox, a_source.box(), QNUM);
      primT.copy(a_primState);
      EBCellFAB cons(m_ebisBox, a_source.box(), CNUM);
      primToCons(cons, primT, m_validBox);
      cons.plus(a_source, 0.5*m_dt);
      consToPrim(primT, cons, m_validBox, 0);

      primPtr = &primT;
    }
  EBCellFAB& prim = *primPtr;
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  int numCons = numConserved();
  int numPrim = numPrimitives();
  CH_assert(isDefined());
  CH_assert(a_ebIrregFlux.nComp() == numCons);
  CH_assert(prim.nComp() == numPrim);
  CH_assert(a_ebIrregFlux.getIVS().contains(a_irregIVS));

  for(VoFIterator vofit(a_irregIVS, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect centroid = m_ebisBox.bndryCentroid(vof);

      CH_assert(prim(vof, QRHO) > 0.0);

      Real     dense  = prim(vof, QRHO);
      Real     press  = prim(vof, QPRES);
      Real     sound  = sqrt(m_gamma*press/dense);
      RealVect veloc;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          veloc[idir] = prim(vof, QVELX+idir);
        }
      Real     deltaDense = 0.0;
      Real     deltaPress = 0.0;
      RealVect deltaVeloc = RealVect::Zero;
      Real dtOver2H = 0.5*m_dt/m_dx[0];
      Real    rhoc2 = dense*sound*sound;
      for(int extrapDir = 0; extrapDir < SpaceDim; extrapDir++)
        {
          Tuple<int, CH_SPACEDIM-1> tanDirs = PolyGeom::computeTanDirs(extrapDir);
          Real centDir = centroid[extrapDir];
          const EBCellFAB& slopePrimDir =  a_slopePrim[extrapDir];

          Real     slopeDense = slopePrimDir(vof, QRHO);
          Real     slopePress = slopePrimDir(vof, QPRES);
          RealVect slopeVeloc;
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              slopeVeloc[idir]= slopePrimDir(vof, QVELX + idir);
            }
          Real slopeUNorm = slopeVeloc[extrapDir];
          Real      velNo =      veloc[extrapDir];

          deltaDense += centDir*slopeDense - dtOver2H*(velNo*slopeDense + dense*slopeUNorm);
          deltaPress += centDir*slopePress - dtOver2H*(velNo*slopePress + rhoc2*slopeUNorm);

          deltaVeloc[extrapDir] += centDir*slopeUNorm - dtOver2H*(velNo*slopeUNorm + slopePress/dense);


          for(int itan = 0; itan < SpaceDim-1; itan++)
            {
              int tanDir = tanDirs[itan];
              Real slopeUTang = slopeVeloc[tanDir];

              deltaVeloc[tanDir] += centDir*slopeUTang - dtOver2H*(velNo*slopeUTang);
            }
        }

      Real xstate[QNUM];
      xstate[QRHO]  = dense + deltaDense;
      xstate[QPRES] = press + deltaPress;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          xstate[QVELX + idir] = veloc[idir] + deltaVeloc[idir];
        }

      if(a_source.isDefined() && (!s_conservativeSource))
        {
          for(int ivar = 0; ivar < QNUM; ivar++)
            {
              Real sourceval = a_source(vof, ivar);
              xstate[ivar] += 0.5*m_dt*sourceval;
            }
        }


      //enforce positivity
      Real pext = xstate[QPRES];
      Real rext = xstate[QRHO];
      if((pext < 0.0) || (rext < 0.0))
        {
          for(int ivar = 0; ivar < numPrim; ivar++)
            {
              xstate[ivar] = prim(vof, ivar);
            }
        }

      //now compute the boundary condition assuming
      //solid walls at the embedded boundary
      RealVect velocity;
      Real     density, pressure;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          velocity[idir] = xstate[QVELX + idir];
        }

      density  = xstate[QRHO];
      pressure = xstate[QPRES];
      RealVect normal = m_ebisBox.normal(vof);
      Real normalvel = PolyGeom::dot(normal, velocity);
      //HACK putin true normal for computing normal vel
      //RealVect truenorm = getTrueNorm(vof);
      //Real normalvel1 = PolyGeom::dot(truenorm, velocity);
      //normalvel = normalvel1;
      //END HACK
      if (density <= 0.0) 
        {
          pout() << "Density is non-positive at " << vofit() << endl;
          // FIXME: Doesn't work in parallel.
          abort();
        }
      CH_assert(density > 0.0);
      pressure = Max(pressure, smallp);
      Real soundspeed = sqrt(m_gamma*pressure/density);
      Real pstar = pressure - density*normalvel*soundspeed;
      if(usesFlattening())
        {
          //limit in case pressure gets too big or too small
          Real ratio = Abs(pstar/pressure);
          if((ratio > 1.1) || (ratio < 0.9))
            {
              pstar = pressure;
            }
        }

      //modification for artificial viscosity
#if 1
      if(usesArtificialViscosity())
        {
          Real coeff = artificialViscosityCoefficient();
          //approximation to the velocity divergence
          Real divu = 0.0;
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              const EBCellFAB& slopeDir = a_slopePrim[idir];
              int ivel = QVELX + idir;
              divu += slopeDir(vof, ivel);
            }
          Real dp = -2.0*density*coeff*normalvel*abs(divu);
          pstar = pstar + dp;
        }
#endif
      //debug
      pstar = Max(pstar, smallp);
      //end debug

      //fluxes as per modiano colella
      a_ebIrregFlux(vof, CRHO) = 0.0;
      a_ebIrregFlux(vof, CENG) = 0.0;
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          int imom = CMOMX + idir;
          //negative because our normal points into the fluid
          //these are not separated out for RZ since the pressure flux
          //is the only flux at the EB
          a_ebIrregFlux(vof, imom) = -pstar*normal[idir];
        }
    }
}
/********/
/*****************************/
void
EBPatchPolytropic::
setCoveredConsVals(EBCellFAB& a_consState)
{
  CH_TIME("EBPatchPolytropic::setCoveredConsVals");
  a_consState.setInvalidData(1.0, CRHO);
  a_consState.setInvalidData(0.0, CMOMX);
  a_consState.setInvalidData(0.0, CMOMY);
#if CH_SPACEDIM==3
  a_consState.setInvalidData(0.0, CMOMZ);
#endif
  a_consState.setInvalidData(1.0, CENG);
}
/**************/
//this is the stable, non-conservative estimate of the solution update
void
EBPatchPolytropic::
nonconservativeDivergence(EBCellFAB&             a_divF,
                          const EBFluxFAB&       a_flux,
                          const BaseIVFAB<Real>  a_coveredFluxMinu[SpaceDim],
                          const BaseIVFAB<Real>  a_coveredFluxPlus[SpaceDim],
                          const Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
                          const Vector<VolIndex> a_coveredFacePlus[SpaceDim],
                          const Box&             a_box)
{
  CH_TIME("EBPatchPolytropic::nonconservativeDiv");
  EBPatchGodunov::nonconservativeDivergence(a_divF, a_flux,
                                            a_coveredFluxMinu,
                                            a_coveredFluxPlus,
                                            a_coveredFaceMinu,
                                            a_coveredFacePlus,
                                            a_box);
}
/*****************************/
void
EBPatchPolytropic::
consUndividedDivergence(BaseIVFAB<Real>&       a_divF,
                        const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                        const BaseIVFAB<Real>& a_ebIrregFlux,
                        const IntVectSet&      a_ivs)
{
  CH_TIME("EBPatchPolytropic::consDiv");
  EBPatchGodunov::consUndividedDivergence(a_divF, a_centroidFlux, a_ebIrregFlux, a_ivs);
}
/*****************************/

#include "NamespaceFooter.H"
