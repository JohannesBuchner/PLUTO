#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBPatchGodunov.H"
#include "EBPatchGodunovF_F.H"
#include "EBArith.H"
#include "PolyGeom.H"
#include "EBDebugOut.H"
#include "DebugOut.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "EBLoHiCenter.H"
#include "EBPatchGodunov.H"
#include "Stencils.H"
#include "EBArith.H"
#include "EBAMRIO.H"
#include "PolyGeom.H"
#include "ParmParse.H"
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>
#include "CH_Timer.H"
#include "NamespaceHeader.H"

bool EBPatchGodunov::s_conservativeSource  = false;
bool EBPatchGodunov::s_verbose  = false;
int  EBPatchGodunov::s_curLevel = -1;
int  EBPatchGodunov::s_curComp  = -1;
int  EBPatchGodunov::s_doingVel  = -1;
int  EBPatchGodunov::s_doingAdvVel  = -1;
Real EBPatchGodunov::s_maxWaveSpeed  = -1.0;
IntVect EBPatchGodunov::s_maxWaveSpeedIV   = IntVect(D_DECL(-1,-1,-1));
IntVect EBPatchGodunov::s_debugIV          = IntVect(D_DECL(16, 3, 0));
int     EBPatchGodunov::s_whichLev=-1;

/*****************************/
void
EBPatchGodunov::
setValidBox(const Box&        a_validBox,
            const EBISBox&    a_ebisBox,
            const IntVectSet& a_coarseFineIVS,
            const Real&       a_cumulativeTime,
            const Real&       a_timestep)
{
  m_isBoxSet = true;
  m_dt       = a_timestep;
  m_time     = a_cumulativeTime;
  m_validBox = a_validBox;
  m_ebisBox  = a_ebisBox;

  //define the interpolation stencils
  IntVectSet ivsIrreg = m_ebisBox.getIrregIVS(m_validBox);
  FaceStop::WhichFaces facestop = FaceStop::SurroundingWithBoundary;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      BaseIFFAB<FaceStencil>& stenDir = m_interpStencils[idir];
      stenDir.define(ivsIrreg, m_ebisBox.getEBGraph(), idir, 1);
      for (FaceIterator faceit(ivsIrreg, m_ebisBox.getEBGraph(), idir, facestop);
          faceit.ok(); ++faceit)
        {
          const  FaceIndex& face = faceit();
          FaceStencil sten = EBArith::getInterpStencil(face, a_coarseFineIVS,
                                                       m_ebisBox, m_domain.domainBox());
          stenDir(face, 0) = sten;
        }
    }
  //I know that a lot of these objects are slightly larger than they
  //need to be.   Be warned, however, that it is
  //really important that all these things are the same size.
  //This allows me to make all kinds of wacky stenciling assumptions down
  //the road.  Mess with these definitions at thy peril.
  m_validBoxG4  = grow(m_validBox, 4);
  m_validBoxG4 &= m_domain;
  int numFlux = numFluxes();
  int numPrim = numPrimitives();
  m_primState.define(   m_ebisBox, m_validBoxG4, numPrim);
  m_primMinuTemp.define(m_ebisBox, m_validBoxG4, numPrim);
  m_primPlusTemp.define(m_ebisBox, m_validBoxG4, numPrim);
  m_primGdnv.define(m_ebisBox, m_validBoxG4, numPrim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_primPlus[idir].define(m_ebisBox, m_validBoxG4,       numPrim);
      m_primMinu[idir].define(m_ebisBox, m_validBoxG4,       numPrim);
      m_fluxOne [idir].define(m_ebisBox, m_validBoxG4, idir, numFlux);

      IntVectSet irregIVSPlus, irregIVSMinu;
      computeCoveredFaces(m_coveredFacePlusG4[idir],
                          m_coveredSetsPlusG4[idir],
                          irregIVSPlus,
                          idir, Side::Hi, m_validBoxG4);
      computeCoveredFaces(m_coveredFaceMinuG4[idir],
                          m_coveredSetsMinuG4[idir],
                          irregIVSMinu,
                          idir, Side::Lo, m_validBoxG4);

      m_extendStatePlusG4  [idir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), numPrim);
      m_extendStateMinuG4  [idir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), numPrim);
      m_extendStateNormPlus[idir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), numPrim);
      m_extendStateNormMinu[idir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), numPrim);
      m_coveredFluxPlusG4  [idir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), numFlux);
      m_coveredFluxMinuG4  [idir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), numFlux);
      m_coveredFluxNormPlus[idir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), numFlux);
      m_coveredFluxNormMinu[idir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), numFlux);

      m_extendStatePlusG4[idir].setVal(0.);
      m_extendStateMinuG4[idir].setVal(0.);
      m_coveredFluxPlusG4[idir].setVal(0.);
      m_coveredFluxMinuG4[idir].setVal(0.);

      if (SpaceDim==3)
        {
          for (int jdir = 0; jdir < SpaceDim; jdir++)
            {
              m_fluxTwo[idir][jdir].define(m_ebisBox, m_validBoxG4, idir, numFlux);

              m_extendStatePlus3D[idir][jdir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), numPrim);
              m_extendStateMinu3D[idir][jdir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), numPrim);
              m_coveredFluxPlus3D[idir][jdir].define(m_coveredSetsPlusG4[idir], m_ebisBox.getEBGraph(), numFlux);
              m_coveredFluxMinu3D[idir][jdir].define(m_coveredSetsMinuG4[idir], m_ebisBox.getEBGraph(), numFlux);
              m_extendStatePlus3D[idir][jdir].setVal(0.);
              m_extendStateMinu3D[idir][jdir].setVal(0.);
              m_coveredFluxPlus3D[idir][jdir].setVal(0.);
              m_coveredFluxMinu3D[idir][jdir].setVal(0.);
            }
        }
    }
  IntVectSet        ivsIrregG4 =  m_ebisBox.getIrregIVS(m_validBoxG4);
  VoFIterator vofit(ivsIrregG4,   m_ebisBox.getEBGraph());
  m_irregVoFs = vofit.getVector();
  m_updateStencil.resize(m_irregVoFs.size());
  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      fillUpdateStencil(m_updateStencil[ivof], m_irregVoFs[ivof]);
    }

  if (m_useAgg)
    setSlopeStuff();

}
/******************/
void
EBPatchGodunov::
cacheEBCF( Vector<Vector<Real> >& a_cache, const  EBCellFAB& a_input)
{
  CH_assert(a_input.getRegion() == m_validBoxG4);

  int nvar = a_input.nComp();
  a_cache.resize(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    {
      a_cache[ivar].resize(m_irregVoFs.size());
    }
  for (int ivar = 0; ivar < nvar; ivar++)
    {
      const Real* inputSV = a_input.getSingleValuedFAB().dataPtr(ivar);
      const Real* inputMV =  a_input.getMultiValuedFAB().dataPtr(ivar);
      for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
        {
          pointerOffset_t& vofOff = m_updateStencil[ivof].m_vofOffset;
          const Real* minuPtr;
          if (vofOff.m_multiValued)
            {
              minuPtr = inputMV + vofOff.m_offset;
            }
          else
            {
              minuPtr = inputSV + vofOff.m_offset;
            }

          a_cache[ivar][ivof] = *minuPtr;
        }
    }
}
/******************/
void
EBPatchGodunov::
uncacheEBCF(EBCellFAB& a_output, const Vector<Vector<Real> >& a_cache)
{
  CH_assert(a_output.getRegion() == m_validBoxG4);
  CH_assert(a_output.nComp()     == a_cache.size());

  for (int ivar = 0; ivar < a_cache.size(); ivar++)
    {
      Real* outputSV = a_output.getSingleValuedFAB().dataPtr(ivar);
      Real* outputMV = a_output.getMultiValuedFAB().dataPtr(ivar);

      for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
        {
          //const VolIndex& vof = m_irregVoFs[ivof];
          pointerOffset_t& vofOff = m_updateStencil[ivof].m_vofOffset;
          Real* minuPtr;
          if (vofOff.m_multiValued)
            {
              minuPtr = outputMV + vofOff.m_offset;
            }
          else
            {
              minuPtr = outputSV + vofOff.m_offset;
            }
          *minuPtr = a_cache[ivar][ivof];
        }
    }
}
/******************/
void
EBPatchGodunov::
fillUpdateStencil(EBPatchGodunov::updateStencil_t& a_stencil, const VolIndex& a_vof)
{
  //i am going to use this assumption in many places
  int numVoFs = m_ebisBox.numVoFs(a_vof.gridIndex());
  if (numVoFs > 1)
    {
      a_stencil.m_vofOffset.m_multiValued = true;
      const BaseIVFAB<Real>& baseivfabPhi= m_primState.getMultiValuedFAB();
      a_stencil.m_vofOffset.m_offset = baseivfabPhi.getIndex(a_vof, 0) - baseivfabPhi.dataPtr(0);
    }
  else
    {
      a_stencil.m_vofOffset.m_multiValued = false;
      IntVect ncells  = m_validBoxG4.size();
      IntVect iv = a_vof.gridIndex() - m_validBoxG4.smallEnd();

      a_stencil.m_vofOffset.m_offset = iv[0] + iv[1]*ncells[0] ;
      if (SpaceDim==3)
        {
          a_stencil.m_vofOffset.m_offset +=    iv[2]*ncells[0]*ncells[1];
        }
    }
}
/****/
void
EBPatchGodunov::
setMaxWaveSpeed(Real a_maxWaveSpeed)
{
  s_maxWaveSpeed = a_maxWaveSpeed;
}
/**/
Real
EBPatchGodunov::
getMaxWaveSpeed()
{
  return s_maxWaveSpeed;
}
/****/
void
EBPatchGodunov::
setMaxWaveSpeedIV(const IntVect& a_maxWaveSpeedIV)
{
  s_maxWaveSpeedIV = a_maxWaveSpeedIV;
}
/**/
IntVect
EBPatchGodunov::
getMaxWaveSpeedIV()
{
  return s_maxWaveSpeedIV;
}
/****/
void
EBPatchGodunov::
setVerbose(bool a_verbose)
{
  s_verbose = a_verbose;
}
/*****/
int
EBPatchGodunov::
getCurComp()
{
  return s_curComp;
}
/*****/
void
EBPatchGodunov::
setCurComp(int a_curComp)
{
  s_curComp = a_curComp;
}
/*****/
/*****/
bool
EBPatchGodunov::
getVerbose()
{
  return s_verbose;
}
/*****/
int
EBPatchGodunov::
getDoingVel()
{
  return s_doingVel;
}
/*****/
int
EBPatchGodunov::
getDoingAdvVel()
{
  return s_doingAdvVel;
}
/*****/
void
EBPatchGodunov::
setDoingVel(int a_doingVel)
{
  s_doingVel = a_doingVel;
}
/*****/
/*****/
void
EBPatchGodunov::
setDoingAdvVel(int a_doingAdvVel)
{
  s_doingAdvVel = a_doingAdvVel;
}
/*****/

/*****/
int
EBPatchGodunov::
getCurLevel()
{
  return s_curLevel;
}
/*****/
void
EBPatchGodunov::
setCurLevel(int a_curLevel)
{
  s_curLevel = a_curLevel;
}
/*****/

EBPatchGodunov::EBPatchGodunov()
{
  m_isDefined     = false;
  m_isBCSet       = false;
  m_isBoxSet      = false;
  m_useAgg = false;
  m_bc = NULL;
}

EBPatchGodunov::~EBPatchGodunov()
{
  if (m_bc != NULL)
    {
      delete m_bc;
    }
}
/******/
void
EBPatchGodunov::
incrementWithSource(EBCellFAB&       a_state,
                    const EBCellFAB& a_source,
                    const Real&      a_scale,
                    const Box&       a_box)
{
  //if conservative source is set to true,  state is conservative state
  //if conservative source is set to false, state is primitive    state
  CH_TIME("EBPatchGodunov::incrmentWithSource");
  CH_assert(m_isBoxSet);
  CH_assert(a_source.nComp() == a_state.nComp());
  CH_assert(a_source.getRegion().contains(a_box));
  CH_assert(a_state.getRegion().contains(a_box));

  BaseFab<Real>&       primReg = a_state.getSingleValuedFAB();
  const BaseFab<Real>& sourReg =    a_source.getSingleValuedFAB();

  FORT_INCSOURCE(CHF_FRA(primReg),
                 CHF_CONST_FRA(sourReg),
                 CHF_CONST_REAL(a_scale),
                 CHF_BOX(a_box));

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for (int ivar = 0; ivar < a_source.nComp(); ivar++)
        {
          Real increm = a_scale*a_source(vof, ivar);
          a_state(vof, ivar) +=  increm;
        }
    }
}
/******/
void
EBPatchGodunov::
extrapToCoveredFaces(BaseIVFAB<Real>&        a_extendedPrim,
                     const EBCellFAB&        a_primMinu,
                     const EBCellFAB&        a_primPlus,
                     const EBCellFAB&        a_primState,
                     const Vector<VolIndex>& a_coveredFaces,
                     const int&              a_faceDir,
                     const Side::LoHiSide&   a_sd,
                     const Box&              a_box)
{
  CH_TIME("EBPatchGodunov::extrapToCoveredFaces");
  for (int ivof = 0; ivof < a_coveredFaces.size(); ivof++)
    {
      const VolIndex& vof = a_coveredFaces[ivof];

      RealVect normal = m_ebisBox.normal(vof);
      if (a_box.contains(vof.gridIndex()))
        {
          const int numPrim = numPrimitives();
          Vector<Real> extPrim(numPrim, 0.0);
          if (SpaceDim== 2)
            {
              pointExtrapToCovered2D(extPrim,
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
EBPatchGodunov::
pointExtrapToCovered2D(Vector<Real>&           a_extrapVal,
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

  int radius = 1;
  Vector<VolIndex> vofsStencil;
  EBArith::getAllVoFsInMonotonePath(vofsStencil, a_vof, m_ebisBox, radius);

  for (int ivar = 0; ivar < a_numPrim; ivar++)
    {
      for (int ix = 0; ix < 2; ix++)
        {
          for (int iy = 0; iy < 2; iy++)
            {
              hasVoF[ix][iy] = EBArith::isVoFHere(vof[ix][iy], vofsStencil, ivSten[ix][iy]);
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
      else if (hasVoF[1][0])
        {//change methods, still second order
          Real deltaW;
          coveredExtrapSlopes(deltaW,vof[1][0],a_primState,a_faceDir,ivar);

          a_extrapVal[ivar] = val[1][0] - 2.0*signNorm*deltaW;
        }
      else if (hasVoFLitNorm)
        {//drop order
          a_extrapVal[ivar] = wLitNorm;
        }
      else if (hasVoFBigNorm)
        {
          a_extrapVal[ivar] = wBigNorm;
        }
      else
        {
          a_extrapVal[ivar] = a_primState(a_vof,ivar);
        }
    }
}

/********/
Real
EBPatchGodunov::pointLimiter(const Real& a_deltaW1,
                             const Real& a_deltaW2)
{
  Real retval;
  if (a_deltaW1*a_deltaW2 > 0.0)
    {
      retval = Min(Abs(a_deltaW1), Abs(a_deltaW2));
      if (a_deltaW1 < 0.0)
        {
          retval *= -1.0;
        }
    }
  else
    {
      retval = 0.0;
    }
  return retval;
}
/******/
Real
EBPatchGodunov::
bilinearFunc(const Real  a_WVal[2][2],
             const Real& a_xd1,
             const Real& a_xd2)
{
  Real D =                a_WVal[0][0];
  Real A =(a_WVal[1][0] - D);
  Real B =(a_WVal[0][1] - D);
  Real C =(a_WVal[1][1] - D - A - B);
  Real retval = A*a_xd1 + B*a_xd2 + C*a_xd1*a_xd2 + D;
  return retval;
}
Real
EBPatchGodunov::
maxFunction(const Real  a_WVal[2][2],
            const Real& a_xd1,
            const Real& a_xd2)
{
  Real retval;
  if (a_xd1 > .5)
    {
      if (a_xd2 > .5)
        {
          retval = a_WVal[1][1];
        }
      else
        {
          retval = a_WVal[1][0];
        }
    }
  else
    {
      if (a_xd2 > .5)
        {
          retval = a_WVal[0][1];
        }
      else
        {
          retval = a_WVal[0][0];
        }
    }
  return retval;
}
void
EBPatchGodunov::
pointExtrapToCovered3D(Vector<Real>&           a_extrapVal,
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

  int radius = 1;
  Vector<VolIndex> vofsStencil;
  EBArith::getAllVoFsInMonotonePath(vofsStencil, a_vof, m_ebisBox, radius);

  for (int ix = 0; ix < 2; ix++)
    {
      for (int iy = 0; iy < 2; iy++)
        {
          hasVoF[ix][iy] =  EBArith::isVoFHere(vofSten[ix][iy], vofsStencil, ivSten[ix][iy]);
        }
    }
  bool hasAllVoFs = hasVoF[0][0] && hasVoF[1][0] && hasVoF[0][1] && hasVoF[1][1];

  //get the vof info for a 1D extrap stencil (if needed)
  VolIndex vofSten1D;
  const IntVect ivSten1D = ivVoF + signNorm*BASISV(a_faceDir);
  const bool has1DVoF =  EBArith::isVoFHere(vofSten1D, vofsStencil, ivSten1D);

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
      else if (has1DVoF)
        {//try 1D extrap
          Real WVal1D;
          if (a_sd == Side::Hi)
            {
              WVal1D = a_primMinu(vofSten1D, ivar);
            }
          else
            {
              WVal1D = a_primPlus(vofSten1D, ivar);
            }
          Real deltaW;
          coveredExtrapSlopes(deltaW,vofSten1D,a_primState,a_faceDir,ivar);
          a_extrapVal[ivar] = WVal1D - 2.0*signNorm*deltaW;
        }
      else
        // At least one of the vofs in the stencil exists - we drop order
        // by using a weighted sum of the the wVal's.
        {
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
                      Real wVal;
                      if (a_sd == Side::Hi)
                        {
                          wVal = a_primMinu(vofSten[ix][iy], ivar);
                        }
                      else
                        {
                          wVal = a_primPlus(vofSten[ix][iy], ivar);
                        }
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
      // None of the vofs in the stencil exists. We use the value at the
      // cell center adjacent to the face.
    }
  //    if (!hasAllVoFs)
  //          cout << "Dropping order at covered face = " << a_vof <<
  //          "direction = " << a_faceDir <<  "\n";
}
/******/
/*****/
void
EBPatchGodunov::
updatePrim(EBCellFAB&              a_primMinu,
           EBCellFAB&              a_primPlus,
           const EBFaceFAB&        a_flux,
           const BaseIVFAB<Real>&  a_coveredFluxMinu,
           const BaseIVFAB<Real>&  a_coveredFluxPlus,
           const Vector<VolIndex>& a_coveredFaceMinu,
           const Vector<VolIndex>& a_coveredFacePlus,
           const int&              a_dir,
           const Box&              a_box,
           const Real&             a_scale)
{
  CH_TIME("EBPatchGodunov::updatePrim");
  CH_assert(isDefined());
  CH_assert(a_primPlus.nComp() == numPrimitives());
  CH_assert(a_primMinu.nComp() == numPrimitives());
  CH_assert(a_flux.nComp() == numFluxes());
  CH_assert((a_dir >= 0) && (a_dir < SpaceDim));
  CH_assert(a_flux.direction() == a_dir);
  CH_assert(a_flux.getRegion().contains(surroundingNodes(a_box, a_dir)));
  CH_assert(a_primPlus.getRegion().contains(a_box));
  CH_assert(a_primMinu.getRegion().contains(a_box));
  CH_assert(a_primPlus.getRegion() == a_primMinu.getRegion());

  const Box& primBox = a_primPlus.getRegion();
  int nCons =  numConserved();
  EBCellFAB consTemp(m_ebisBox, m_validBoxG4, nCons);
  primToCons(consTemp, a_primMinu, primBox);

  updateCons(consTemp, a_flux,
             a_coveredFluxMinu, a_coveredFluxPlus,
             a_coveredFaceMinu, a_coveredFacePlus,
             a_dir, a_box, a_scale);
  int logflag = 0;
  consToPrim(a_primMinu, consTemp, m_validBoxG4, logflag);

  primToCons(consTemp, a_primPlus, primBox);
  updateCons(consTemp, a_flux,
             a_coveredFluxMinu, a_coveredFluxPlus,
             a_coveredFaceMinu, a_coveredFacePlus,
             a_dir, a_box, a_scale);
  consToPrim(a_primPlus, consTemp, primBox, logflag);

}
/*****/
void
EBPatchGodunov::
updateCons(EBCellFAB&              a_consState,
           const EBFaceFAB&        a_flux,
           const BaseIVFAB<Real>&  a_coveredFluxMinu,
           const BaseIVFAB<Real>&  a_coveredFluxPlus,
           const Vector<VolIndex>& a_coveredFaceMinu,
           const Vector<VolIndex>& a_coveredFacePlus,
           const int&              a_dir,
           const Box&              a_box,
           const Real&             a_scale)
{
  CH_TIME("EBPatchGodunov::updateCons");
  CH_assert(isDefined());
  CH_assert(a_consState.nComp() <= a_flux.nComp());
  CH_assert(a_consState.nComp() == numConserved() );
  CH_assert((a_dir >= 0) && (a_dir < SpaceDim));
  CH_assert(a_flux.direction() == a_dir);
  CH_assert(a_flux.getRegion().contains(surroundingNodes(a_box, a_dir)));
  CH_assert(a_consState.getRegion().contains(a_box));

  Vector<Vector<Real> > cache;
  {
    CH_TIME("update cons cache");
    cacheEBCF(cache, a_consState);
  }
  {
    CH_TIME("EBPatchGodunov::updateConsRegular");
    int ncons = numConserved();
    BaseFab<Real>&       consReg = a_consState.getSingleValuedFAB();
    const BaseFab<Real>& fluxReg = a_flux.getSingleValuedFAB();
    FORT_UPDATE( CHF_BOX(a_box),
                 CHF_FRA(consReg),
                 CHF_CONST_FRA(fluxReg),
                 CHF_CONST_INT(a_dir),
                 CHF_CONST_INT(ncons),
                 CHF_CONST_REAL(a_scale));
  }

  {
    CH_TIME("update cons uncache");
    uncacheEBCF(a_consState, cache);
  }
  //update the irregular vofs
  {
    CH_TIME("EBPatchGodunov::updateConsIrregular");
    for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        if (a_box.contains(vof.gridIndex()))
          {
            for (SideIterator sit; sit.ok(); ++sit)
              {
                int isign = sign(sit());
                Vector<FaceIndex> faces = m_ebisBox.getFaces(vof, a_dir, sit());
                for (int ivar = 0; ivar < a_consState.nComp(); ivar++)
                  {
                    //aggregate the fluxes over the faces for a_dir and sit()
                    Real aggFlux = 0.0;
                    for (int iface = 0; iface < faces.size(); iface++)
                      {
                        const FaceIndex& face = faces[iface];
                        Real faceFlux = a_flux(face, ivar);
                        aggFlux += faceFlux;
                      }

                    if (faces.size() > 1)
                      {
                        aggFlux /= Real(faces.size());
                      }

                    //dx is already divided out (part of scale)
                    aggFlux *= -isign*a_scale;

                     //a_consState modified in regular update, but restored in uncache
                    a_consState(vof, ivar) += aggFlux;
                  }
              }
          }
      }
  }
  {
    CH_TIME("EBPatchGodunov::updateCons::Covered");
    //add in covered face fluxes
    for (SideIterator sit; sit.ok(); ++sit)
      {
        int isign = sign(sit());
        const BaseIVFAB<Real>*  coveredFluxPtr;
        const Vector<VolIndex>* coveredFacePtr;
        if (sit() == Side::Lo)
          {
            coveredFluxPtr = &(a_coveredFluxMinu);
            coveredFacePtr = &(a_coveredFaceMinu);
          }
        else
          {
            coveredFluxPtr = &(a_coveredFluxPlus);
            coveredFacePtr = &(a_coveredFacePlus);
          }
        const BaseIVFAB<Real>&  coveredFlux = *coveredFluxPtr;
        const Vector<VolIndex>& coveredFace = *coveredFacePtr;
        for (int ivof = 0; ivof < coveredFace.size(); ivof++)
          {
            const VolIndex& vof = coveredFace[ivof];
            //the input set can be bigger than the box over which we are integrating.
            //this was done for peformance reasonas so the sets to not get
            //calculated over and over.
            if (a_box.contains(vof.gridIndex()))
              {
                for (int ivar = 0; ivar < a_consState.nComp(); ivar++)
                  {

                    Real covFlux = coveredFlux(vof, ivar);
                    Real update =   isign*covFlux;

                    //dx is already divided out (part of scale)
                    update *= -a_scale;
                    Real state = a_consState(vof, ivar);
                    a_consState(vof, ivar) = state + update;

                  }
              }
          }
      }
  }

  {
    CH_TIME("EBPatchGodunov::updateCons::floors");
    floorConserved(a_consState, a_box);
  }
}
/******/
bool
EBPatchGodunov::
isDefined() const
{
  return
    m_isDefined &&
    m_isBCSet &&
    m_isBoxSet &&
    m_isSlopeSet &&
    m_isArtViscSet;
}
/*****************************/
void
EBPatchGodunov::
setSlopeParameters(bool a_useFourthOrderSlopes,
                   bool a_useFlattening,
                   bool a_useLimiting)
{
  // Slope flattening is only allowed with 4th order slopes
  CH_assert(a_useFourthOrderSlopes || !a_useFlattening);

  // Store the slope computation parameters
  m_useFourthOrderSlopes = a_useFourthOrderSlopes;
  m_useFlattening        = a_useFlattening;
  m_useLimiting          = a_useLimiting;

  m_isSlopeSet = true;
}
/*****************************/
const EBPhysIBC*
EBPatchGodunov::
getEBPhysIBC() const
{
  CH_assert(m_isBCSet);
  return m_bc;
}
/*****************************/
void
EBPatchGodunov::
setEBPhysIBC(const EBPhysIBCFactory& a_bcFact)
{
  // Delete old boundary condition object - if any
  if (m_bc != NULL)
    {
      delete m_bc;
    }
  m_bc = a_bcFact.create();
  if (m_isDefined)
    m_bc->define(m_domain, m_dx);
  m_isBCSet = true;
}
/*****************************/
void
EBPatchGodunov::
define(const ProblemDomain&  a_domain,
       const RealVect& a_dx,
       bool a_useAgg)
{
  // Store the domain and grid spacing
  m_isDefined = true;
  m_domain = a_domain;
  // Set the domain and grid spacing in the boundary condition object
  m_dx = a_dx;
  m_useAgg = a_useAgg;
  //figure out dxScale (used for unscaling the boundary area)
  {
    Real maxDx=0.0,dv=1.0;
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        dv *= m_dx[idir];
        if ( m_dx[idir] > maxDx )
          {
            maxDx = m_dx[idir];
          }
      }
    m_dxScale = pow(maxDx,SpaceDim-1)/dv;
  }

  if (m_isBCSet)
    m_bc->define(a_domain, a_dx);
}
/******/
void
EBPatchGodunov::
applyLimiter(EBCellFAB&       a_slopePrim,
             const EBCellFAB& a_slopeLeft,
             const EBCellFAB& a_slopeRigh,
             const int&       a_dir,
             const Box&       a_box)
{
  CH_TIME("EBPatchGodunov::applyLimiter");
  CH_assert(isDefined());
  CH_assert(a_slopePrim.getRegion().contains(a_box));
  CH_assert(a_slopeLeft.getRegion().contains(a_box));
  CH_assert(a_slopeRigh.getRegion().contains(a_box));
  CH_assert((a_dir >= 0) && (a_dir < SpaceDim));
  CH_assert(a_slopePrim.nComp() == numSlopes());
  CH_assert(a_slopeLeft.nComp() == numSlopes());
  CH_assert(a_slopeRigh.nComp() == numSlopes());

  const BaseFab<Real>& regSlopeLeft = a_slopeLeft.getSingleValuedFAB();
  const BaseFab<Real>& regSlopeRigh = a_slopeRigh.getSingleValuedFAB();
  BaseFab<Real>& regSlopePrim       = a_slopePrim.getSingleValuedFAB();

  {
    CH_TIME("regular limiting");
    FORT_VLLIMITER(CHF_FRA(regSlopePrim),
                   CHF_CONST_FRA(regSlopeLeft),
                   CHF_CONST_FRA(regSlopeRigh),
                   CHF_BOX(a_box));
  }

  {
    CH_TIME("irregular limiting");
    IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
    for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
        vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        for (int ivar = 0; ivar < numSlopes(); ivar++)
          {
            //the output array arrives with the centered diff value
            const Real& dql = a_slopeLeft(vof, ivar);
            const Real& dqr = a_slopeRigh(vof, ivar);
            Real& dqlim = a_slopePrim(vof, ivar);

            FORT_POINTVLLIMITER(CHF_REAL(dqlim),
                                CHF_CONST_REAL(dql),
                                CHF_CONST_REAL(dqr));
          }
      }
  }
}
/*****************************/
void
EBPatchGodunov::
slope(EBCellFAB&       a_slopePrim,
      EBCellFAB&       a_slopeNLim,
      const EBCellFAB& a_primState,
      const EBCellFAB& a_flattening,
      const int&       a_dir,
      const Box&       a_box,
      bool a_doAggregated)
{
  CH_TIME("EBPatchGodunov::slope");
  int numSlope = numSlopes();
  CH_assert(a_slopePrim.nComp() == numSlopes());
  CH_assert(a_primState.nComp() >= numSlopes());
  EBCellFAB delta2W, deltaWL, deltaWR, deltaWC;

  {
    CH_TIME("second_order_slopes");
    doSecondOrderSlopes(delta2W,
                        deltaWL,
                        deltaWR,
                        deltaWC,
                        a_primState,
                        a_dir,
                        a_box,
                        a_doAggregated);
  }

  a_slopeNLim.copy(deltaWC);

  if (usesFourthOrderSlopes())
    {
      EBCellFAB delta4W;
      {
        CH_TIME("forth_order_slopes");
        doFourthOrderSlopes(delta4W,
                            deltaWC,
                            delta2W,
                            deltaWL,
                            deltaWR,
                            a_primState,
                            a_dir,
                            a_box);
      }

      //apply flattening coefficient
      //would be nice to just use *= operator here
      //but they have different numbers of variables
      if (usesFlattening())
        {
          CH_TIME("apply_flatttening");
          Box entireBox = delta4W.getRegion();
          const BaseFab<Real>& regFlattening = a_flattening.getSingleValuedFAB();
          BaseFab<Real>& regDelta4W = delta4W.getSingleValuedFAB();

          FORT_APPLYFLAT(CHF_FRA(regDelta4W),
                         CHF_CONST_FRA1(regFlattening, 0),
                         CHF_CONST_INT(numSlope),
                         CHF_BOX(entireBox));

          // all single valued cells multiplied in fortran
          IntVectSet multiIVS  = m_ebisBox.getMultiCells(entireBox);
          for (VoFIterator vofit(multiIVS, m_ebisBox.getEBGraph());
              vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              for (int ivar = 0; ivar < numSlope; ivar++)
                {
                  delta4W(vof, ivar) *= a_flattening(vof, 0);
                }
            }
        }
        {
          CH_TIME("slope_copy1");
          a_slopePrim.copy(delta4W);
        }
    }//end if (usefourthorder slopes())
  else
    {
      {
        CH_TIME("slope_copy2");
        a_slopePrim.copy(delta2W);
      }
    }
  //debug
  //a_slopePrim.setVal(0.);
  //end debug
}
/*****************************/
void
EBPatchGodunov::
doSecondOrderSlopes(EBCellFAB&       a_delta2W,
                    EBCellFAB&       a_deltaWL,
                    EBCellFAB&       a_deltaWR,
                    EBCellFAB&       a_deltaWC,
                    const EBCellFAB& a_primState,
                    const int&       a_dir,
                    const Box&       a_box,
                    bool a_doAggregated)
{
  CH_TIME("EBPG::doSecondOrderSlopes");

  int numSlope;
  Box box1, box2;
  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;
  {
    //    CH_TIME("doSecondOrderSlopes_prelims");
    CH_assert(m_isSlopeSet);
    numSlope = numSlopes();
    box1 = a_box;
    box1.grow(a_dir, 1);
    box1 &= m_domain;

    box2 = a_box;
    box2.grow(a_dir, 2);
    box2 &= m_domain;

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
  }
  {
    CH_TIME("EBPG::doSecondOrderSlopes_fab_defs");
    a_delta2W.define(m_ebisBox, entireBox, numSlope);
    a_deltaWL.define(m_ebisBox, entireBox, numSlope);
    a_deltaWR.define(m_ebisBox, entireBox, numSlope);
    a_deltaWC.define(m_ebisBox, entireBox, numSlope);
  }
  //  { not necessary
  //  CH_TIME("EBPG::doSecondOrderSlopes_zeros");
  //  a_delta2W.setVal(0.);
  //  a_deltaWL.setVal(0.);
  //  a_deltaWR.setVal(0.);
  //  a_deltaWC.setVal(0.);
  //  }

    BaseFab<Real>& regDelta2W = a_delta2W.getSingleValuedFAB();
    BaseFab<Real>& regDeltaWC = a_deltaWC.getSingleValuedFAB();
    BaseFab<Real>& regDeltaWL = a_deltaWL.getSingleValuedFAB();
    BaseFab<Real>& regDeltaWR = a_deltaWR.getSingleValuedFAB();
    const BaseFab<Real>& regPrimState = a_primState.getSingleValuedFAB();

    {
      CH_TIME("EBPG::doSecondOrderSlopes_regular_calc");
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
    }
  //apply van leer limiter to regular cells.
  //limiting for irregular cells is more complicated
  //now that we are doing higher order limited
  //one-sided diffs
    {
      CH_TIME("copy_wc_to_2w");
      regDelta2W.copy(regDeltaWC);
    }
  if (m_useLimiting)
    {
      CH_TIME("limiting");
      FORT_VLLIMITER(CHF_FRA(regDelta2W),
                     CHF_CONST_FRA(regDeltaWL),
                     CHF_CONST_FRA(regDeltaWR),
                     CHF_BOX(centerBox));
    }

  if (a_doAggregated && m_useAgg)
    {
      //the slow version has some funky flattening logic that I do
      //not have the time to program into the fast version right now as I
      //am trying to make a deadline for a code that does not use flattening.
      if (usesFlattening())
        {
          MayDay::Warning("agg irreg slopes does not contain flattening limiter");
        }
      aggIrregSecondOrderSlopes(a_delta2W,
                                a_deltaWL,
                                a_deltaWR,
                                a_deltaWC,
                                a_primState,
                                a_dir,
                                entireBox);
    }
  else
    {
      irregSecondOrderSlopes(a_delta2W,
                             a_deltaWL,
                             a_deltaWR,
                             a_deltaWC,
                             a_primState,
                             a_dir,
                             a_box);
    }

  //want to be able to call this in test codes
  if (m_isBCSet)
    {
      //      CH_TIME("boundary_slopes");
      m_bc->setBndrySlopes(a_delta2W, a_primState, m_ebisBox, entireBox, a_dir);
    }
}
void
EBPatchGodunov::
irregSecondOrderSlopes(EBCellFAB&       a_delta2W,
                       EBCellFAB&       a_deltaWL,
                       EBCellFAB&       a_deltaWR,
                       EBCellFAB&       a_deltaWC,
                       const EBCellFAB& a_primState,
                       const int&       a_dir,
                       const Box&       a_entireBox)
{
  //at irregular cells, if one side does not exist,
  //set all to one-sided diffs.  try to do higher order
  //limited diffs if possible.
  {
    CH_TIME("irreg_slopes");
    for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        if (a_entireBox.contains(vof.gridIndex()))
          {
            for (int ivar = 0; ivar < numSlopes(); ivar++)
              {
                bool hasFacesLeft, hasFacesRigh;
                Real dql, dqr, dqc;
                bool verbose =false;
                pointGetSlopes(dql, dqr,dqc,
                               hasFacesLeft,
                               hasFacesRigh,
                               vof, a_primState, a_dir, ivar, verbose);
                Real dqSec=0;
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
                        if (dqc*dql > 0.0)
                          {
                            Real rsign = 1.0;
                            if (dqc < 0.0)
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
                        if (dqc*dqr > 0.0)
                          {
                            Real rsign = 1.0;
                            if (dqc < 0.0)
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
                        MayDay::Error("EBPatchGodunov::doSecondOrderSlopes -- missed a case");
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
              } //end loop over variables
          }
      }//end loop over vofs
  }
}
/*****************************/
void
EBPatchGodunov::
getArgBox(Box a_argBox[SpaceDim])
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (SpaceDim==2)
        {
          Box slopeBoxG2 = grow(m_validBox, 2);
          slopeBoxG2 &= m_domain;
          a_argBox[idir] = slopeBoxG2;
        }
      else
        {
          Box modBoxOpen = m_validBox;
          modBoxOpen.grow(3);
          modBoxOpen &= m_domain;
          a_argBox[idir] = modBoxOpen;
        }
    }
}
/*****************************/
void
EBPatchGodunov::
setSlopeStuff()
{
  //first get the elusive and irritating m_entirebox

  //this is the box sent to slope routines
  Box argBox[SpaceDim];
  getArgBox(argBox);

  //this is the raindance the slope routine uses to generate entirebox
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box box1, box2;
      Box loBox, hiBox, centerBox;
      int hasLo, hasHi;

      box1 = argBox[idir];
      box1.grow(idir, 1);
      box1 &= m_domain;

      box2 = argBox[idir];
      box2.grow(idir, 2);
      box2 &= m_domain;

      if (usesFourthOrderSlopes())
        {
          eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, m_entireBox[idir],
                       box2, m_domain, idir);
        }
      else
        {
          eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, m_entireBox[idir],
                       box1, m_domain, idir);
        }
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      IntVectSet        ivsIrregG4 =  m_ebisBox.getIrregIVS(m_validBoxG4);
      ivsIrregG4 &= m_entireBox[idir] ;
      VoFIterator vofit(ivsIrregG4,   m_ebisBox.getEBGraph());
      const Vector<VolIndex>& vofs = vofit.getVector();
      m_slopVec[idir].resize(vofs.size());
      Vector<VoFStencil> slowStenLo(vofs.size());
      Vector<VoFStencil> slowStenHi(vofs.size());
      //for getting offsets
      EBCellFAB dumprim(m_ebisBox,      m_validBoxG4, numPrimitives());
      EBCellFAB dumslop(m_ebisBox, m_entireBox[idir], numSlopes());
      for (int ivof = 0; ivof < vofs.size(); ivof++)
        {
          const VolIndex& vof = vofs[ivof];
          const IntVect&  iv = vof.gridIndex();
          bool onLeftDomain = (iv[idir] == m_domain.domainBox().smallEnd(idir));
          bool onRighDomain = (iv[idir] == m_domain.domainBox().bigEnd(  idir));

          Vector<FaceIndex> facesLo = (m_ebisBox.getFaces(vof, idir, Side::Lo));
          Vector<FaceIndex> facesHi = (m_ebisBox.getFaces(vof, idir, Side::Hi));

          //the faces.size() == 1 thing should probably be done as >= 0
          //and we figure out what we want to do in the case of multiple faces.
          //as it is, I do not want to change the algorithm in the midst of trying
          //to do optimization.   Since pointgetslopes does it this way, this is the
          //way we shall do it.
          bool hasFacesLo = (facesLo.size()==1) && !onLeftDomain;
          bool hasFacesHi = (facesHi.size()==1) && !onRighDomain;

          m_slopVec[idir][ivof].hasLo = hasFacesLo;
          m_slopVec[idir][ivof].hasHi = hasFacesHi;
          m_slopVec[idir][ivof].slop_access.offset = dumslop.offset(  vof, 0);
          m_slopVec[idir][ivof].slop_access.dataID = dumslop.dataType(vof);
          if (hasFacesLo)
            {
              slowStenLo[ivof].add(facesLo[0].getVoF(Side::Lo), -1.0);
              slowStenLo[ivof].add(                        vof,  1.0);
            }
          if (hasFacesHi)
            {
              slowStenHi[ivof].add(facesHi[0].getVoF(Side::Hi),  1.0);
              slowStenHi[ivof].add(                        vof, -1.0);
            }
        }

      Vector<RefCountedPtr<BaseIndex   > > baseindice(vofs.size());
      Vector<RefCountedPtr<BaseStencil > > basestenlo(vofs.size());
      Vector<RefCountedPtr<BaseStencil > > basestenhi(vofs.size());
      for (int ivof= 0; ivof < vofs.size(); ivof++)
        {
          baseindice[ivof] =   RefCountedPtr<BaseIndex>(new   VolIndex(      vofs[ivof]));
          basestenlo[ivof] = RefCountedPtr<BaseStencil>(new VoFStencil(slowStenLo[ivof]));
          basestenhi[ivof] = RefCountedPtr<BaseStencil>(new VoFStencil(slowStenHi[ivof]));
        }

      m_slopStenLo[idir] = RefCountedPtr< AggStencil <EBCellFAB, EBCellFAB > >
        (new AggStencil <EBCellFAB, EBCellFAB >(baseindice, basestenlo, dumprim, dumslop));
      m_slopStenHi[idir] = RefCountedPtr< AggStencil <EBCellFAB, EBCellFAB > >
        (new AggStencil <EBCellFAB, EBCellFAB >(baseindice, basestenhi, dumprim, dumslop));

    }



}

/******************/
void
EBPatchGodunov::
aggIrregSecondOrderSlopes(EBCellFAB&       a_delta2W,
                          EBCellFAB&       a_deltaWL,
                          EBCellFAB&       a_deltaWH,
                          EBCellFAB&       a_deltaWC,
                          const EBCellFAB& a_primState,
                          const int&       a_dir,
                          const Box&       a_entireBox)
{
  CH_TIME("agg_irreg_second_order_slopes");
  CH_assert(m_useAgg);
  CH_assert(a_entireBox == m_entireBox[a_dir]);
  CH_assert(a_delta2W.getRegion() == m_entireBox[a_dir]);
  CH_assert(a_deltaWL.getRegion() == m_entireBox[a_dir]);
  CH_assert(a_deltaWH.getRegion() == m_entireBox[a_dir]);
  CH_assert(a_deltaWC.getRegion() == m_entireBox[a_dir]);
  CH_assert(a_primState.getRegion() == m_validBoxG4);
  //the slow version has some funky flattening logic that I do
  //not have the time to program into the fast version right now as I
  //am trying to make a deadline for a code that does not use flattening.

  //get low and high slopes (or left and right if you prefer)
  m_slopStenLo[a_dir]->apply(a_deltaWL, a_primState, 0, 0,  numSlopes(), false);
  m_slopStenHi[a_dir]->apply(a_deltaWH, a_primState, 0, 0,  numSlopes(), false);

  //at irregular cells, if one side does not exist,
  //set all to one-sided diffs.  try to do higher order
  //limited diffs if possible.

    for (int ivof = 0; ivof< m_slopVec[a_dir].size(); ivof++)
      {
        //slopes have all the same box so they will have the same offset and type
        const size_t& offsetS = m_slopVec[a_dir][ivof].slop_access.offset;
        const int   & dataIDS = m_slopVec[a_dir][ivof].slop_access.dataID;

        const bool& hasFacesLo = m_slopVec[a_dir][ivof].hasLo;
        const bool& hasFacesHi = m_slopVec[a_dir][ivof].hasHi;

        for (int ivar = 0; ivar < numSlopes(); ivar++)
          {
            Real*       ptrLo =   a_deltaWL.dataPtr(dataIDS,ivar);
            Real*       ptrHi =   a_deltaWH.dataPtr(dataIDS,ivar);
            Real*       ptrCe =   a_deltaWC.dataPtr(dataIDS,ivar);
            Real*       ptrSe =   a_delta2W.dataPtr(dataIDS,ivar);

            Real& dql = *(ptrLo + offsetS);
            Real& dqr = *(ptrHi + offsetS);
            Real& dqc = *(ptrCe + offsetS);
            Real& dqS = *(ptrSe + offsetS);
            //this logic from pointgetslopes to get dqc
            if (hasFacesLo && hasFacesHi)
              {
                dqc = 0.5*(dql+dqr);
                if (m_useLimiting)
                  {
                    Real dqmin = Min(Abs(dql)*2.0,Abs(dqr)*2.0);
                    if (dql*dqr < 0.)
                      {
                        dqmin = 0.;
                      }
                    if (dqc < 0.)
                      {
                        dqc = -Min(dqmin,Abs(dqc));
                      }
                    else
                      {
                        dqc = Min(dqmin,Abs(dqc));
                      }
                  }
              }
            else if (!hasFacesLo && !hasFacesHi)
              {
                dql = 0.0;
                dqr = 0.0;
                dqc = 0.0;
              }
            else if (hasFacesLo && !hasFacesHi)
              {
                dqr = dql;
                //no higher-order one sided diffs
                dqc = dql;
              }
            else if (hasFacesHi && !hasFacesLo)
              {
                dql = dqr;
                //no higher-order one sided diffs
                dqc = dqr;
              }
            else
              {
                MayDay::Error("EBPatchGodunov::aggirregpointGetSlopes -- missed a case");
              }

            //now for the limiting bits from irregSlopes
            //looks like a lot of this is redundant
            if (!m_useLimiting)
              {
                dqS = dqc;
              }
            else
              {
                if (hasFacesLo && hasFacesHi)
                  {
                    Real dqlim = dqc;

                    FORTNT_POINTVLLIMITER(CHF_REAL(dqlim),
                                          CHF_CONST_REAL(dql),
                                          CHF_CONST_REAL(dqr));
                    dqS = dqlim;
                  }
                else if (!hasFacesLo && !hasFacesHi)
                  {
                    dqS = 0.0;
                  }
                else if (hasFacesLo && !hasFacesHi)
                  {
                    if (dqc*dql > 0.0)
                      {
                        Real rsign = 1.0;
                        if (dqc < 0.0)
                          {
                            rsign = -1.0;
                          }
                        dqS = rsign*Min(Abs(dqc), Abs(dql));
                      }
                    else
                      {
                        dqS = 0.0;
                      }

                  }
                else if (hasFacesHi && !hasFacesLo)
                  {
                    if (dqc*dqr > 0.0)
                      {
                        Real rsign = 1.0;
                        if (dqc < 0.0)
                          {
                            rsign = -1.0;
                          }
                        dqS = rsign*Min(Abs(dqc), Abs(dqr));
                      }
                    else
                      {
                        dqS = 0.0;
                      }
                  }
                else
                  {
                    MayDay::Error("EBPatchGodunov::aggirreg2-- missed a case");
                  }
              }//end if using limiting

          } //end loop over variables
      }//end loop over vofs
}


/*****************************/
void EBPatchGodunov::
coveredExtrapSlopes(Real&            a_dq,
                    const VolIndex&  a_vof,
                    const EBCellFAB& a_primState,
                    const int&       a_dir,
                    const int&       a_ivar)
{
  CH_TIME("EBPatchGodunov::coveredExtrapSlopes");
  bool hasFacesLeft, hasFacesRigh;
  Real dql, dqr;
  bool verbose = false;
  pointGetSlopes(dql, dqr,a_dq,
                 hasFacesLeft,
                 hasFacesRigh,
                 a_vof, a_primState, a_dir, a_ivar, verbose);
  if (usesFlattening())
    {

      int pindex = pressureIndex();
      int rindex = densityIndex();
      Real press = Max(a_primState(a_vof, pindex), Real(1.0e-10));
      Real dense = Max(a_primState(a_vof, rindex), Real(1.0e-10));
      Real dqp, dqd;

      verbose = false;
      pointGetSlopes(dql, dqr, dqp,
                     hasFacesLeft,
                     hasFacesRigh,
                     a_vof, a_primState, a_dir, pindex,verbose);

      //try to detect pressure discontinutity
      dqp = Max(Abs(dqp), Abs(dql));
      dqp = Max(Abs(dqp), Abs(dqr));
      pointGetSlopes(dql, dqr, dqd,
                     hasFacesLeft,
                     hasFacesRigh,
                     a_vof, a_primState, a_dir, rindex,verbose);

      //try to detect denisty discontinutity
      dqd = Max(Abs(dqd), Abs(dql));
      dqd = Max(Abs(dqd), Abs(dqr));
      if ((Abs(dqp/press) > 0.1) || (Abs(dqd/dense) > 0.1))
        {
          a_dq = 0.0;
        }
    }
  //debug! set slopes to zero!
  // a_dq = 0;
  //end debug
}
/*****************************/
void EBPatchGodunov::
pointGetSlopes(Real&            a_dql,
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
   const IntVect&  iv = a_vof.gridIndex();
  //one-sided diffs on domain bndry
  const Box& domainBox = m_domain.domainBox();
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
  //limiting sprinked in here to make higher order
  //one-sided diffs possible
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
      if (m_useLimiting)
        {
          Real dqmin = Min(Abs(a_dql)*2.0,Abs(a_dqr)*2.0);
          if (a_dql*a_dqr < 0.)
            {
              dqmin = 0.;
            }
          if (a_dqc < 0.)
            {
              a_dqc = -Min(dqmin,Abs(a_dqc));
            }
          else
            {
              a_dqc = Min(dqmin,Abs(a_dqc));
            }
        }
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
      MayDay::Error("EBPatchGodunov::pointGetSlopes -- a_dqc != a_dqc");
    }
  //atttempts to detect discontinuties in constrained vars
  if (usesFlattening())
    {
      if ((Abs(a_dqc) > Abs(0.1*valCent)) && ((a_ivar == densityIndex()) || (a_ivar == pressureIndex())))
        {
          a_dql = 0.0;
          a_dqr = 0.0;
          a_dqc = 0.0;
        }
    }
  //
}
/*****************************/
void EBPatchGodunov::
doFourthOrderSlopes(EBCellFAB&       a_delta4W,
                    EBCellFAB&       a_deltaWC,
                    const EBCellFAB& a_delta2W,
                    const EBCellFAB& a_deltaWL,
                    const EBCellFAB& a_deltaWR,
                    const EBCellFAB& a_primState,
                    const int&       a_dir,
                    const Box&       a_box)
{
  CH_TIME("dofourthorderslopes");
  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;
  int numSlope = numSlopes();
  const Box& domainBox = m_domain.domainBox();
  Box box1 = a_box;
  box1.grow(a_dir, 1);
  box1 &= m_domain;
  eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
               box1, m_domain, a_dir);

  a_delta4W.define(m_ebisBox, entireBox, numSlope);

  BaseFab<Real>& regDelta4W = a_delta4W.getSingleValuedFAB();

  const BaseFab<Real>& regDelta2W = a_delta2W.getSingleValuedFAB();
  //const BaseFab<Real>& regDeltaWC = a_deltaWC.getSingleValuedFAB();
  //const BaseFab<Real>& regDeltaWL = a_deltaWL.getSingleValuedFAB();
  //const BaseFab<Real>& regDeltaWR = a_deltaWR.getSingleValuedFAB();
  const BaseFab<Real>& regPrimState = a_primState.getSingleValuedFAB();
  {
    CH_TIME("reg4thOrderSlopes");
    FORT_FORTHSLOPEDIFFS(CHF_FRA(regDelta4W),
                         CHF_CONST_FRA(regPrimState),
                         CHF_CONST_FRA(regDelta2W),
                         CHF_CONST_INT(numSlope),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(loBox),
                         CHF_CONST_INT(hasLo),
                         CHF_BOX(hiBox),
                         CHF_CONST_INT(hasHi),
                         CHF_BOX(centerBox));
  }

  {
    CH_TIME("irreg4thOrderSlopes");
    for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        if (a_box.contains(vof.gridIndex()))
          {
            const IntVect&  iv = vof.gridIndex();
            //one-sided diffs on domain bndry
            bool onLeftDomain = iv[a_dir] == domainBox.smallEnd(a_dir);
            bool onRighDomain = iv[a_dir] == domainBox.bigEnd(a_dir);
            bool hasFacesLeft = (m_ebisBox.numFaces(vof, a_dir, Side::Lo) > 0) && !onLeftDomain;
            bool hasFacesRigh = (m_ebisBox.numFaces(vof, a_dir, Side::Hi) > 0) && !onRighDomain;

            for (int ivar = 0; ivar < numSlopes(); ivar++)
              {
                Real dq4;
                if (hasFacesRigh && hasFacesLeft)
                  {
                    Vector<FaceIndex> facesLeft =
                      m_ebisBox.getFaces(vof, a_dir, Side::Lo);
                    Real valLeft = 0.0;
                    Real sloLeft = 0.0;
                    for (int iface = 0; iface < facesLeft.size(); iface++)
                      {
                        VolIndex vofLeft = facesLeft[iface].getVoF(Side::Lo);
                        valLeft += a_primState(vofLeft, ivar);
                        sloLeft +=      a_delta2W(vofLeft, ivar);
                      }
                    valLeft /= Real(facesLeft.size());
                    sloLeft /= Real(facesLeft.size());

                    Vector<FaceIndex> facesRigh =
                      m_ebisBox.getFaces(vof, a_dir, Side::Hi);
                    Real valRigh = 0.0;
                    Real sloRigh = 0.0;
                    for (int iface = 0; iface <facesRigh.size(); iface++)
                      {
                        VolIndex vofRigh = facesRigh[iface].getVoF(Side::Hi);
                        valRigh += a_primState(vofRigh, ivar);
                        sloRigh +=     a_delta2W(vofRigh, ivar);
                      }
                    valRigh /= Real(facesRigh.size());
                    sloRigh /= Real(facesRigh.size());

                    Real dwr = valRigh - 0.25*sloRigh;
                    Real dwl = valLeft + 0.25*sloLeft;
                    dq4 = (2.0/3.0)*(dwr - dwl);

                  }
                else if (hasFacesRigh || hasFacesLeft)
                  {
                    dq4 = a_delta2W(vof, ivar);
                  }
                else
                  {
                    //no faces on either side.  set the slope to zero
                    dq4 = 0.0;
                  }
                a_delta4W(vof, ivar) = dq4;
              } //end loop over variables
          }
      }//end loop over vofs
  }

  if (m_useLimiting)
    {
      CH_TIME("4th order limiting");
      applyLimiter(a_delta4W, a_deltaWL, a_deltaWR,  a_dir, centerBox);
    }

}
/*****************************/
void
EBPatchGodunov::
regularUpdate(EBCellFAB&          a_consState,
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
  int numCons = numConserved();

  CH_TIME("EBPatchGodunov::regularUpdate");
  CH_assert(isDefined());

  CH_assert(a_flux.nComp() >= numConserved());
  CH_assert(a_consState.nComp() >= numConserved());
  setCoveredConsVals(a_consState);

  EBCellFAB slopePrim[SpaceDim];
  EBCellFAB slopeNLim[SpaceDim];
  Vector<VolIndex> coveredFacePlus[SpaceDim];
  Vector<VolIndex> coveredFaceMinu[SpaceDim];

  computeFluxes(a_flux,
                m_coveredFluxMinuG4, m_coveredFluxPlusG4,
                coveredFaceMinu,     coveredFacePlus,
                m_primState,    slopePrim, slopeNLim,
                a_flattening, a_consState, a_source,
                a_box, a_dit,a_verbose);

  //now that we have the fluxes, modify them with
  //artificial viscosity if that is called for.
  //the artificial viscosity correction to the
  //embedded boundary flux is computed inside
  //computeebirregflux
  if (usesArtificialViscosity())
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

  //compute nonconservative divergences of fluxes
  EBCellFAB nonConsDivF(m_ebisBox, a_box, numConserved());

  //this is the stable, non-conservative estimate of the solution update
  nonconservativeDivergence(nonConsDivF, a_flux,
                            m_coveredFluxMinuG4,
                            m_coveredFluxPlusG4,
                            coveredFaceMinu,
                            coveredFacePlus,
                            a_box);

  //copy nonconservative div f into sparse output thingy
  IntVectSet ivsIrreg = a_ivsSmall;
  for (VoFIterator vofit(ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for (int ivar = 0; ivar < numCons; ivar++)
        {
          a_nonConservativeDivergence(vof, ivar) = nonConsDivF(vof, ivar);
        }
    }

  //compute irregular boundary flux.  this includes
  //an artificial viscosity modification if appropriate
  computeEBIrregFlux( a_ebIrregFlux, m_primState,
                      slopeNLim, ivsIrreg, a_source);

  BaseFab<Real>& regConsState   =
    a_consState.getSingleValuedFAB();
  const BaseFab<Real>& regNonConsDivF =
    nonConsDivF.getSingleValuedFAB();
  //increment solution with dt*source everywhere
  if (a_source.isDefined())
    {
      if (s_conservativeSource)
        {
          incrementWithSource(a_consState, a_source,m_dt, a_box);
        }
      else
        {
          int logflag = 0;
          EBCellFAB  primState(m_ebisBox, a_box, numPrimitives());
          consToPrim(primState, a_consState, a_box, logflag);
          incrementWithSource(primState, a_source,m_dt, a_box);
          primToCons(a_consState, primState, a_box);
        }

    }
  //this updates all single valued cells with the non-conservative divergence.
  FORT_REGUPDATE( CHF_BOX(a_box),
                  CHF_CONST_FRA(regConsState),
                  CHF_CONST_FRA(regNonConsDivF),
                  CHF_CONST_INT(numCons),
                  CHF_CONST_REAL(m_dt));

  //this does the multi-valued cells.
  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for (int ivar = 0; ivar < numCons; ivar++)
        {
          Real update = -m_dt*nonConsDivF(vof, ivar);
          a_consState(vof, ivar) += update;
        }
    }

  //irregular cells get updated in irregularUpdate
  //floors also happen there.
  // with the linear combination of
  //divergences thing, taking into account that nonConsDivF has
  //already been included once.
}

/*****************************/
void
EBPatchGodunov::
regularDivergences(EBCellFAB&          a_nonConsDivF,
                   EBCellFAB&          a_consState,
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
  CH_TIME("EBPatchGodunov::regularUpdate");
  CH_assert(isDefined());
  int numCons = numConserved();

  CH_assert(a_flux.nComp() >= numConserved());
  CH_assert(a_consState.nComp() >= numConserved());
  setCoveredConsVals(a_consState);

  EBCellFAB slopePrim[SpaceDim];
  EBCellFAB slopeNLim[SpaceDim];
  Vector<VolIndex> coveredFacePlus[SpaceDim];
  Vector<VolIndex> coveredFaceMinu[SpaceDim];

  computeFluxes(a_flux,
                m_coveredFluxMinuG4, m_coveredFluxPlusG4,
                  coveredFaceMinu,     coveredFacePlus,
                m_primState,    slopePrim, slopeNLim,
                a_flattening, a_consState, a_source,
                a_box, a_dit,a_verbose);

  //now that we have the fluxes, modify them with
  //artificial viscosity if that is called for.
  //the artificial viscosity correction to the
  //embedded boundary flux is computed inside
  //computeebirregflux
  if (usesArtificialViscosity())
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
                            coveredFaceMinu,
                            coveredFacePlus,
                            a_box);

  //copy nonconservative div f into sparse output thingy
  IntVectSet ivsIrreg = a_ivsSmall;
  for (VoFIterator vofit(ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for (int ivar = 0; ivar < numCons; ivar++)
        {
          a_nonConservativeDivergence(vof, ivar) = a_nonConsDivF(vof, ivar);
        }
    }

  //compute irregular boundary flux.  this includes
  //an artificial viscosity modification if appropriate
  computeEBIrregFlux( a_ebIrregFlux, m_primState,
                      slopeNLim, ivsIrreg, a_source);
}
void
EBPatchGodunov::
nonconservativeDivergence(EBCellFAB&             a_divF,
                          const EBFluxFAB&       a_flux,
                          const BaseIVFAB<Real>  a_coveredFluxMinu[SpaceDim],
                          const BaseIVFAB<Real>  a_coveredFluxPlus[SpaceDim],
                          const Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
                          const Vector<VolIndex> a_coveredFacePlus[SpaceDim],
                          const Box&             a_box)
{
  CH_TIME("EBPatchGodunov::nonconservativeDivergence");
  CH_assert(m_isDefined);
  CH_assert(m_isBoxSet);
  CH_assert(a_divF.nComp() >= numConserved());
  CH_assert(a_flux.nComp() >= numConserved());
  CH_assert(a_flux.getRegion().contains(a_box));
  CH_assert(a_divF.getRegion().contains(a_box));

  //set the divergence initially to zero
  //then loop through directions and increment the divergence
  //with each directions flux difference.
  int ncons = numConserved();
  a_divF.setVal(0.0);
  BaseFab<Real>&       regDivF = a_divF.getSingleValuedFAB();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      //update for the regular vofs in the nonconservative
      //case  works for all single valued vofs.
      const EBFaceFAB& fluxDir = a_flux[idir];

      const BaseFab<Real>& regFluxDir = fluxDir.getSingleValuedFAB();

      /* do the regular vofs */

      FORT_DIVERGEF( CHF_BOX(a_box),
                     CHF_FRA(regDivF),
                     CHF_CONST_FRA(regFluxDir),
                     CHF_CONST_INT(idir),
                     CHF_CONST_INT(ncons),
                     CHF_CONST_REAL(m_dx[idir]));

    }
  //update the irregular vofsn
    for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
      {
        const VolIndex& vof = m_irregVoFs[ivof];
        if (a_box.contains(vof.gridIndex()))
          {
            //divergence was set in regular update.  we reset it
            // to zero and recalc.
            for (int ivar = 0; ivar < ncons; ivar++)
              {
                Real irregDiv = 0.0;
                for (int idir = 0; idir < SpaceDim; idir++)
                  {
                    const EBFaceFAB& fluxDir = a_flux[idir];
                    for (SideIterator sit; sit.ok(); ++sit)
                      {
                        int isign = sign(sit());
                        Vector<FaceIndex> faces =
                          m_ebisBox.getFaces(vof, idir, sit());
                        Real update = 0.;
                        for (int iface = 0; iface < faces.size(); iface++)
                          {
                            const FaceIndex& face = faces[iface];
                            Real flux = fluxDir(face, ivar);
                            update += isign*flux;

                          }
                        if (faces.size() > 1)
                          update /= Real(faces.size());
                        irregDiv += update/m_dx[idir];
                      } //end loop over sides
                  }//end loop over directions
                a_divF(vof, ivar) = irregDiv;
              }//end loop over variables
          }
    }//end loop over irreg vofs

  //now correct for the covered fluxes
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int isign = sign(sit());
          const BaseIVFAB<Real>*  coveredFluxPtr;
          const Vector<VolIndex>* coveredFacePtr;
          if (sit() == Side::Lo)
            {
              coveredFluxPtr = &a_coveredFluxMinu[idir];
              coveredFacePtr = &a_coveredFaceMinu[idir];
            }
          else
            {
              coveredFluxPtr = &a_coveredFluxPlus[idir];
              coveredFacePtr = &a_coveredFacePlus[idir];
            }
          const BaseIVFAB<Real>&  coveredFlux = *coveredFluxPtr;
          const Vector<VolIndex>& coveredFace = *coveredFacePtr;
          for (int ivof = 0; ivof < coveredFace.size(); ivof++)
            {
              const VolIndex& vof = coveredFace[ivof];
              //the sets can be bigger than the box for performance reasons.
              if (a_box.contains(vof.gridIndex()))
                {
                  //face on this side is covered.  use covered flux.
                  for (int ivar = 0; ivar < ncons; ivar++)
                    {
                      Real flux = coveredFlux(vof, ivar);
                      Real update = isign*flux/m_dx[idir];
                      a_divF(vof, ivar) += update;
                    }
                }
            }
        }
    }
}
/*****************************/
void
EBPatchGodunov::
consUndividedDivergence(BaseIVFAB<Real>&       a_divF,
                        const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                        const BaseIVFAB<Real>& a_ebIrregFlux,
                        const IntVectSet&      a_ivs)
{
  CH_TIME("EBPatchGodunov::consUndividedDivergence");
  CH_assert(m_isDefined);
  CH_assert(m_isBoxSet);
  int ncons = a_divF.nComp();
  a_divF.setVal(0.0);
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Real bndryArea = m_ebisBox.bndryArea(vof);
      for (int ivar = 0; ivar < ncons ; ivar++)
        {
          Real bndryFlux = a_ebIrregFlux(vof,ivar);
          Real update = 0.;
          for ( int idir = 0; idir < SpaceDim; idir++)
            {
              const BaseIFFAB<Real>& fluxDir = a_centroidFlux[idir];
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  int isign = sign(sit());
                  Vector<FaceIndex> faces =
                    m_ebisBox.getFaces(vof, idir, sit());
                  for (int iface = 0; iface < faces.size(); iface++)
                    {
                      const FaceIndex& face = faces[iface];
                      Real areaFrac = m_ebisBox.areaFrac(face);
                      Real faceFlux = fluxDir(face, ivar);
                      update += isign*areaFrac*faceFlux/m_dx[idir];
                    }
                }
            }

          //add EB boundary conditions in divergence
          update += bndryFlux*bndryArea*m_dxScale;
          //note NOT divided by volfrac
          a_divF(vof, ivar) = update;
        } //end loop over variables
    } //end loop over vofs
}
void
EBPatchGodunov::
interpolateFluxToCentroids(BaseIFFAB<Real>              a_centroidFlux[SpaceDim],
                           const BaseIFFAB<Real>* const a_fluxInterpolant[SpaceDim],
                           const IntVectSet&            a_irregIVS)
{
  CH_TIME("EBPatchGodunov::interpolateFluxToCentroids");
  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;

  //could just use numFluxes but this allows unit testing
  //with a single variable interpolant
  int nflux = a_centroidFlux[0].nComp();
  //now loop through the irregular faces
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      BaseIFFAB<Real>& fluxDir = a_centroidFlux[faceDir];
      const BaseIFFAB<Real>& interpol = *(a_fluxInterpolant[faceDir]);
      const BaseIFFAB<FaceStencil>& stencils = m_interpStencils[faceDir];

      for (FaceIterator faceit(a_irregIVS, m_ebisBox.getEBGraph(), faceDir, stopCrit);
          faceit.ok(); ++faceit)
        {
          const FaceIndex& centFace = faceit();
          const FaceStencil& sten = stencils(centFace, 0);
          for (int ivar = 0; ivar < nflux; ivar++)
            {
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

          //debug!  turn off interpolation
          // for (int ivar = 0; ivar < numFluxes(); ivar++)
          //   {
          //     fluxDir(centFace, ivar) = interpol(centFace, ivar);
          //   }
        }
    }
}
/*****************************/
void
EBPatchGodunov::
finalUpdate(EBCellFAB&              a_consState,
            BaseIVFAB<Real>&        a_massDiff,
            const BaseIVFAB<Real>&  a_nonConsDivF,
            const BaseIVFAB<Real>&  a_conservDivF,
            const IntVectSet&       a_ivs)
{
  CH_TIME("EBPatchGodunov::finalUpdate");
  CH_assert(a_nonConsDivF.nComp() == numConserved());
  CH_assert(a_conservDivF.nComp() == numConserved());
  CH_assert(a_consState.nComp()   == numConserved());
  CH_assert(a_massDiff.nComp()    == numConserved());

  //initialize mass difference to zero
  a_massDiff.setVal(0.0);
  int ncons = numConserved();
  //update the irreg vofs
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      //      const IntVect& iv = vof.gridIndex();
      Real volFrac = m_ebisBox.volFrac(vof);

      for (int ivar = 0; ivar < ncons ; ivar++)
        {
          //volume fraction not divided out when divergence
          //calculated.  so we do not multiply here.

          //conservative divergence is already
          //multiplied by volfrac
          Real kapConsDiv = a_conservDivF(vof, ivar);
          Real ncDivergeF = a_nonConsDivF(vof, ivar);

          //dt*DF^NC has already been added in regular update
          //subtract back off to get back to the old state
          Real updateIncrementNC =  -m_dt*ncDivergeF;
          Real oldState = a_consState(vof, ivar) - updateIncrementNC;

          //un+1 = un -dt(kappa*divfcons + (1-kapp)divfnc)
          Real newState = oldState
            - m_dt*kapConsDiv
            - m_dt*(1.0-volFrac)*ncDivergeF;
          //dm =-(1-kappa)*kappa*dt*divfcons- kappa*dt*(1-kappa)*divfnc
          //dm =-(1-kappa)*dt*(kappa*divfcons- kappa*divfnc)
          Real deltaM =
            -(1.0-volFrac)*m_dt*(kapConsDiv - volFrac*ncDivergeF);

          a_consState(vof, ivar) = newState;
          a_massDiff(vof, ivar) = deltaM;
        }
    }
}

void
EBPatchGodunov::
computeFluxes(EBFluxFAB&       a_flux,
              BaseIVFAB<Real>  a_coveredFluxMinu[SpaceDim],
              BaseIVFAB<Real>  a_coveredFluxPlus[SpaceDim],
              Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
              Vector<VolIndex> a_coveredFacePlus[SpaceDim],
              EBCellFAB&       a_primState,
              EBCellFAB        a_slopePrim[SpaceDim],
              EBCellFAB        a_slopeNLim[SpaceDim],
              const EBCellFAB& a_flattening,
              const EBCellFAB& a_consState,
              const EBCellFAB& a_source,
              const Box&       a_box,
              const DataIndex& a_dit,
              bool             a_verbose)
{

  if (SpaceDim==2)
    {
      extrapolatePrim2D(m_primMinu, m_primPlus,
                        a_primState,    a_slopePrim, a_slopeNLim,
                        a_flattening, a_consState, a_source,
                        a_box, a_dit,a_verbose);
    }
  else if (SpaceDim==3)
    {

      extrapolatePrim3D(m_primMinu,   m_primPlus,
                        a_primState,  a_slopePrim, a_slopeNLim,
                        a_flattening, a_consState, a_source,
                        a_box, a_dit,a_verbose);
    }
  else
    {
      MayDay::Error("bogus SpaceDim");
    }

  IntVectSet        coveredSetsPlus[SpaceDim];
  IntVectSet        coveredSetsMinu[SpaceDim];
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
      IntVectSet irregIVSPlus, irregIVSMinu;
      computeCoveredFaces(a_coveredFacePlus[idir],
                          coveredSetsPlus[idir],
                          irregIVSPlus,
                          idir, Side::Hi, a_box);
      computeCoveredFaces(a_coveredFaceMinu[idir],
                          coveredSetsMinu[idir],
                          irregIVSMinu,
                          idir, Side::Lo, a_box);
    }
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      //extrapolate to covered faces using updated extrapolated state
      extrapToCoveredFaces(m_extendStateMinuG4[faceDir],
                           m_primMinu[faceDir],
                           m_primPlus[faceDir],
                           a_primState,
                           a_coveredFaceMinu[faceDir],
                           faceDir, Side::Lo, a_box);

      extrapToCoveredFaces(m_extendStatePlusG4[faceDir],
                           m_primMinu[faceDir],
                           m_primPlus[faceDir],
                           a_primState,
                           a_coveredFacePlus[faceDir],
                           faceDir, Side::Hi, a_box);

      //equation 1.19 top
      riemann(a_flux[faceDir], m_primPlus[faceDir], m_primMinu[faceDir],
              faceDir, faceBox[faceDir]);

      // Use the user supplied PhysBC object to obtain boundary fluxes
      m_bc->fluxBC(a_flux, a_primState, m_primMinu[faceDir],
                   Side::Lo,  m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[faceDir], faceDir);
      m_bc->fluxBC(a_flux, a_primState, m_primPlus[faceDir],
                   Side::Hi,  m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[faceDir], faceDir);

      //solve riemann problem between extended state and
      //value in vof to get covered flux.
      //Equation 1.19 bottom.
      riemann(a_coveredFluxMinu[faceDir],
              m_extendStateMinuG4[faceDir], m_primMinu[faceDir],
              a_coveredFaceMinu[faceDir], faceDir, Side::Lo, a_box);
      riemann(a_coveredFluxPlus[faceDir],
              m_extendStatePlusG4[faceDir], m_primPlus[faceDir],
              a_coveredFacePlus[faceDir], faceDir, Side::Hi, a_box);
    }
}
/******/
/******/
void
EBPatchGodunov::
extrapolatePrim2D(EBCellFAB           a_primMinu[SpaceDim],
                  EBCellFAB           a_primPlus[SpaceDim],
                  EBCellFAB&          a_primState,
                  EBCellFAB           a_slopesPrim[SpaceDim],
                  EBCellFAB           a_slopesSeco[SpaceDim],
                  const EBCellFAB&    a_flattening,
                  const EBCellFAB&    a_consState,
                  const EBCellFAB&    a_source,
                  const Box&          a_box,
                  const DataIndex&    a_dit,
                  bool                a_verbose)
{
  CH_TIME("EBPatchGodunov::extrapolatePrim2D");

  //now define the plethora of data holders that I need


  //now do the actual computation.
  //1. transform to primitive state
  int logflag = 0;
  Box primBox = a_consState.box();
  consToPrim(a_primState, a_consState, primBox, logflag);

  //2.  compute slopes delta^d w
  //compute normal derivative stuff

  doNormalDerivativeExtr2D(a_primMinu,
                           a_primPlus,
                           m_fluxOne,
                           m_coveredFluxNormMinu,
                           m_coveredFluxNormPlus,
                           m_coveredFaceMinuG4,
                           m_coveredFacePlusG4,
                           a_slopesPrim,
                           a_slopesSeco,
                           a_flattening,
                           a_primState,
                           a_source,
                           a_dit,
                           a_box);


  /**/
  // Do the final corrections to the fluxes
  //this only happens on the non-ghosted box
  finalExtrap2D(a_primMinu,
                a_primPlus,
                m_coveredFluxNormMinu,
                m_coveredFluxNormPlus,
                m_coveredFaceMinuG4,
                m_coveredFacePlusG4,
                m_fluxOne,
                a_primState,
                a_slopesPrim,
                a_slopesSeco,
                a_box
                );

  /**/



}

void
EBPatchGodunov::
extrapolatePrim3D(EBCellFAB           a_primMinu[SpaceDim],
                  EBCellFAB           a_primPlus[SpaceDim],
                  EBCellFAB&          a_primState,
                  EBCellFAB           a_slopesPrim[SpaceDim],
                  EBCellFAB           a_slopesSeco[SpaceDim],
                  const EBCellFAB&    a_flattening,
                  const EBCellFAB&    a_consState,
                  const EBCellFAB&    a_source,
                  const Box&          a_box,
                  const DataIndex&    a_dit,
                  bool                a_verbose)
{
  CH_TIME("EBPatchGodunov::extrapolatePrim3D");
  //define the plethora of data holders that I need

  //now do the actual computation.
  //1. transform to primitive state
  int numPrim = numPrimitives();
  Box primBox = a_consState.getRegion();
  CH_assert(m_domain.contains(primBox));
  a_primState.define(m_ebisBox, primBox, numPrim);
  int logflag = 0;
  consToPrim(a_primState, a_consState, primBox, logflag);

  //2.  compute slopes delta^d w
  //compute normal derivative stuff
  doNormalDerivativeExtr3D(a_primMinu,
                           a_primPlus,
                           m_fluxOne,
                           m_coveredFluxNormMinu,
                           m_coveredFluxNormPlus,
                           m_coveredFaceMinuG4,
                           m_coveredFacePlusG4,
                           a_slopesPrim,
                           a_slopesSeco,
                           a_flattening,
                           a_primState,
                           a_source,
                           a_dit,
                           a_box);

  //6. In 3D compute corrections to u corresponding to one set
  //   of transverse derivatives appropriate to obtain (1, 1, 1)
  //   coupling.
  // In 3D, compute some additional intermediate fluxes
  // NOTE:  The diagonal entries of this array of fluxes are not
  // used and will not be defined.

  do111coupling(m_fluxTwo,
                m_coveredFluxMinu3D,
                m_coveredFluxPlus3D,
                a_primMinu,
                a_primPlus,
                m_coveredFluxNormMinu,
                m_coveredFluxNormPlus,
                m_coveredFaceMinuG4,
                m_coveredFacePlusG4,
                m_fluxOne,
                a_primState,
                a_slopesPrim,
                a_slopesSeco,
                a_dit,
                a_box);

  finalExtrap3D(a_primMinu,
                a_primPlus,
                m_coveredFluxMinu3D,
                m_coveredFluxPlus3D,
                m_fluxTwo,
                a_primState,
                a_slopesPrim,
                a_slopesSeco,
                a_box);
}
/******************/
void EBPatchGodunov::
finalExtrap2D(EBCellFAB              a_primMinu[SpaceDim],
              EBCellFAB              a_primPlus[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxNormMinu[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxNormPlus[SpaceDim],
              const Vector<VolIndex> a_coveredFaceNormMinu[SpaceDim],
              const Vector<VolIndex> a_coveredFaceNormPlus[SpaceDim],
              const EBFaceFAB        a_fluxOne[SpaceDim],
              const EBCellFAB&       a_primState,
              const EBCellFAB        a_slopesPrim[SpaceDim],
              const EBCellFAB        a_slopesSeco[SpaceDim],
              const Box&             a_box)

{
  CH_TIME("EBPatchGodunov::finalExtrap2D");
  Box slopeBoxG1 = grow(a_box, 1);
  slopeBoxG1 &= m_domain;

  //7. compute fnal corrections to W due to the
  //final transverse derivatives.
  // Do the final corrections to the fluxes
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      // Correct the flux using fluxes in the remaining direction(s)
      for (int diffDir = 0; diffDir < SpaceDim; diffDir++)
        {
          // A different direction has been found
          if (diffDir != faceDir)
            {

              // In 2D,
              //the current primitive state is updated by a flux in
              // the other direction
              //a_coveredFlux holds the covered flux
              //due to the normal derivative
              //equation 1.18 line 3.
              updatePrim(a_primMinu[faceDir],
                         a_primPlus[faceDir], a_fluxOne[diffDir],
                         a_coveredFluxNormMinu[diffDir], a_coveredFluxNormPlus[diffDir],
                         a_coveredFaceNormMinu[diffDir], a_coveredFaceNormPlus[diffDir],
                         diffDir,  slopeBoxG1, (0.5)*m_dt/m_dx[diffDir]);

            } //end if dir2 != faceDir
        } //loop over dir2
    } //loop over facedir
}
/***********/
void EBPatchGodunov::
finalExtrap3D(EBCellFAB              a_primMinu[SpaceDim],
              EBCellFAB              a_primPlus[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxMinu3D[SpaceDim][SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxPlus3D[SpaceDim][SpaceDim],
              const EBFaceFAB        a_fluxTwo[SpaceDim][SpaceDim],
              const EBCellFAB&       a_primState,
              const EBCellFAB        a_slopesPrim[SpaceDim],
              const EBCellFAB        a_slopesSeco[SpaceDim],
              const Box&             a_box)

{
  CH_TIME("EBPatchGodunov::finalExtrap3D");
  Box slopeBoxG1 = grow(a_box, 1);
  slopeBoxG1 &= m_domain;

  //7. compute fnal corrections to W due to the
  //final transverse derivatives.
  // Do the final corrections to the fluxes
  for (int dir1 = 0; dir1 < SpaceDim; dir1++)
    {
      // Correct the flux using fluxes in the remaining direction(s)
      for (int dir2 = 0; dir2 < SpaceDim; dir2++)
        {
          // A different direction has been found
          if (dir2 != dir1)
            {
              // In 3D,  find a direction different from the two above
              // here we need to use the coveredFlux1D for the covered flux

              int dir3 = 3 - dir2  - dir1;
              // A different direction has been found
              //const EBFluxFAB& fluxTwo = a_fluxTwo[dir3];
              // Update the conservative state
              //using both corrected fluxes in
              // the other two directions.  Equation 1.18 line 4
              updatePrim(a_primMinu[dir1], a_primPlus[dir1],
                         a_fluxTwo[dir2][dir3],
                         a_coveredFluxMinu3D[dir2][dir3],
                         a_coveredFluxPlus3D[dir2][dir3],
                         m_coveredFaceMinuG4[dir2],
                         m_coveredFacePlusG4[dir2],
                         dir2,  slopeBoxG1,  (0.5)*m_dt/m_dx[dir2]);

            }
        } //loop over dir2

    } //loop over dir1

}
/*******************/
void EBPatchGodunov::
do111coupling(EBFaceFAB              a_fluxTwo[SpaceDim][SpaceDim],
              BaseIVFAB<Real>        a_coveredFluxMinu3D[SpaceDim][SpaceDim],
              BaseIVFAB<Real>        a_coveredFluxPlus3D[SpaceDim][SpaceDim],
              const EBCellFAB        a_primMinu[SpaceDim],
              const EBCellFAB        a_primPlus[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxNormMinu[SpaceDim],
              const BaseIVFAB<Real>  a_coveredFluxNormPlus[SpaceDim],
              const Vector<VolIndex> a_coveredFaceNormMinu[SpaceDim],
              const Vector<VolIndex> a_coveredFaceNormPlus[SpaceDim],
              const EBFaceFAB        a_fluxOne[SpaceDim],
              const EBCellFAB&       a_primState,
              const EBCellFAB        a_slopesPrim[SpaceDim],
              const EBCellFAB        a_slopesSeco[SpaceDim],
              const DataIndex&       a_dit,
              const Box&             a_box)
{
  CH_TIME("EBPatchGodunov::do111coupling");
  Box slopeBoxG1 = grow(a_box, 1);
  Box slopeBoxG2 = grow(a_box, 2);
  slopeBoxG1 &= m_domain;
  slopeBoxG2 &= m_domain;
  Box faceBox[SpaceDim][SpaceDim];
  Box bndryFaceBox[SpaceDim][SpaceDim];
  Box modBoxCov[SpaceDim];
  Box modBoxOpen[SpaceDim][SpaceDim];
  for (int faceDir = 0; faceDir < SpaceDim; ++faceDir)
    {
      for (int diffDir = 0; diffDir < SpaceDim; ++diffDir)
        {

          modBoxOpen[faceDir][diffDir] = grow(a_box, 2);
          modBoxOpen[faceDir][diffDir].grow(diffDir, -1);
          modBoxOpen[faceDir][diffDir] &= m_domain;
        }
    }

  for (int modDir = 0; modDir < SpaceDim; ++modDir)
    {
      modBoxCov[modDir] = grow(a_box, 1);
      modBoxCov[modDir] &= m_domain;
      for (int faceDir = 0; faceDir < SpaceDim; ++faceDir)
        {
          bndryFaceBox[faceDir][modDir] = modBoxCov[modDir];
          bndryFaceBox[faceDir][modDir].surroundingNodes(faceDir);

          faceBox[faceDir][modDir] = modBoxCov[modDir];
          faceBox[faceDir][modDir].grow(faceDir, 1);
          faceBox[faceDir][modDir] &= m_domain;
          faceBox[faceDir][modDir].grow(faceDir, -1);
          faceBox[faceDir][modDir].surroundingNodes(faceDir);
        }
    }

  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      for (int diffDir = 0; diffDir < SpaceDim; diffDir++)
        {
          if (diffDir != faceDir)
            {
              a_coveredFluxPlus3D[faceDir][diffDir].setVal(0.);
              a_coveredFluxMinu3D[faceDir][diffDir].setVal(0.);
              a_fluxTwo[faceDir][diffDir].setVal(0.);
              m_extendStateMinu3D[faceDir][diffDir].setVal(0.);
              m_extendStatePlus3D[faceDir][diffDir].setVal(0.);
            }
        }
    }
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      for (int diffDir = 0; diffDir < SpaceDim; diffDir++)
        {
          if (diffDir != faceDir)
            {
              int modDir = 3 - faceDir - diffDir;

              // Copy data for in place modification
              m_primMinuTemp.copy(a_primMinu[faceDir]);
              m_primPlusTemp.copy(a_primPlus[faceDir]);

              // Update the current, extrapolated primitive
              //variable using a flux
              // in a different direction.  uses covered fluxes calculated
              // above.  Equation 1.12
              updatePrim(m_primMinuTemp, m_primPlusTemp, a_fluxOne[diffDir],
                         a_coveredFluxNormMinu[diffDir], a_coveredFluxNormPlus[diffDir],
                         a_coveredFaceNormMinu[diffDir], a_coveredFaceNormPlus[diffDir],
                         diffDir,  modBoxOpen[faceDir][diffDir],  (1.0/3.0)*m_dt/m_dx[diffDir]);

              extrapToCoveredFaces(m_extendStateMinu3D[faceDir][diffDir],
                                   m_primMinuTemp,
                                   m_primPlusTemp,
                                   a_primState,
                                   m_coveredFaceMinuG4[faceDir],
                                   faceDir, Side::Lo, modBoxCov[modDir]);

              extrapToCoveredFaces(m_extendStatePlus3D[faceDir][diffDir],
                                   m_primMinuTemp,
                                   m_primPlusTemp,
                                   a_primState,
                                   m_coveredFacePlusG4[faceDir],
                                   faceDir, Side::Hi, modBoxCov[modDir]);

              // Solve the Riemann problem and get fluxes.  Eqution 1.15
              riemann(a_fluxTwo[faceDir][diffDir], m_primPlusTemp, m_primMinuTemp,
                      faceDir, faceBox[faceDir][modDir]);

              //some wackiness to work around the fact that fluxfab did
              //not exist when a lot of this was written.
              Vector<EBFaceFAB*> bcfluxes(SpaceDim);
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  bcfluxes[idir] = &(a_fluxTwo[idir][diffDir]);
                }
              EBFluxFAB fluxAlias;
              fluxAlias.alias(bcfluxes);

              m_bc->fluxBC(fluxAlias,  a_primState, m_primMinuTemp,
                           Side::Lo,  m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[faceDir][diffDir], faceDir);
              m_bc->fluxBC(fluxAlias,  a_primState, m_primPlusTemp,
                           Side::Hi,  m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[faceDir][diffDir], faceDir);

              //solve riemann problem between extended state and
              //value in vof to get covered flux.  Equation 1.14
              riemann(a_coveredFluxMinu3D[faceDir][diffDir],
                      m_extendStateMinu3D[faceDir][diffDir], m_primMinuTemp,
                      m_coveredFaceMinuG4[faceDir], faceDir, Side::Lo, modBoxOpen[faceDir][diffDir]);
              riemann(a_coveredFluxPlus3D[faceDir][diffDir],
                      m_extendStatePlus3D[faceDir][diffDir], m_primPlusTemp,
                      m_coveredFacePlusG4[faceDir], faceDir, Side::Hi, modBoxOpen[faceDir][diffDir]);

            }
        }
    }
}
/*******************/
void EBPatchGodunov::
doNormalDerivativeExtr3D(EBCellFAB              a_primMinu[SpaceDim],
                         EBCellFAB              a_primPlus[SpaceDim],
                         EBFaceFAB              a_fluxOne[SpaceDim],
                         BaseIVFAB<Real>        a_coveredFluxNormMinu[SpaceDim],
                         BaseIVFAB<Real>        a_coveredFluxNormPlus[SpaceDim],
                         Vector<VolIndex>       a_coveredFaceNormMinu[SpaceDim],
                         Vector<VolIndex>       a_coveredFaceNormPlus[SpaceDim],
                         EBCellFAB              a_slopesPrim[SpaceDim],
                         EBCellFAB              a_slopesSeco[SpaceDim],
                         const EBCellFAB&       a_flattening,
                         const EBCellFAB&       a_primState,
                         const EBCellFAB&       a_source,
                         const DataIndex&       a_dit,
                         const Box&             a_box)

{
  CH_TIME("EBPatchGodunov::doNormalDerivative3D");
  Box modBoxCov[SpaceDim];
  Box modBoxOpen[SpaceDim];
  Box faceBox[SpaceDim];
  Box bndryFaceBox[SpaceDim];
  int numSlop = numSlopes();

  //set up the data structures
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      modBoxOpen[idir] = a_box;
      modBoxOpen[idir].grow(3);
      modBoxOpen[idir].grow(idir, -1);
      modBoxOpen[idir] &= m_domain;

      modBoxCov[idir] = a_box;
      modBoxCov[idir].grow(2);
      modBoxCov[idir].grow(idir, -1);
      modBoxCov[idir] &= m_domain;

      bndryFaceBox[idir] = modBoxCov[idir];
      bndryFaceBox[idir].surroundingNodes(idir);

      faceBox[idir] = modBoxCov[idir];
      faceBox[idir].grow(idir, 1);
      faceBox[idir] &= m_domain;
      faceBox[idir].grow(idir, -1);
      faceBox[idir].surroundingNodes(idir);

      a_slopesPrim[idir].define(m_ebisBox, modBoxOpen[idir], numSlop);
      a_slopesSeco[idir].define(m_ebisBox, modBoxOpen[idir], numSlop);


      m_extendStateNormPlus[idir].setVal(0.);
      m_extendStateNormMinu[idir].setVal(0.);
      a_coveredFluxNormPlus[idir].setVal(0.);
      a_coveredFluxNormMinu[idir].setVal(0.);
      a_primPlus[idir].setVal(0.);
      a_primMinu[idir].setVal(0.);
      a_fluxOne[idir].setVal(0.);
    }
  //all one big loop since no cross terms here
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      //check flattening coefficients have correct box and compute slopes
      if (usesFlattening())
        CH_assert(a_flattening.getRegion().contains(modBoxOpen[idir]));

      //compute slopes
      slope(a_slopesPrim[idir], a_slopesSeco[idir], a_primState,
            a_flattening, idir, modBoxOpen[idir]);

      //3. Compute the effect of the normal derivative terms and the
      //source term on the extrapolation in space and time from cell
      //centers to faces. Equation 1.7 (top).
      normalPred(a_primMinu[idir], a_primPlus[idir],
                 a_primState, a_slopesPrim[idir],
                 m_dt/m_dx[idir], idir, modBoxOpen[idir]);

      // If the source term is valid add it to the primitive quantities
      //if source term is not defined, assume it is zero.  lots of
      //problems will have no sources.  not the most elegant method
      //but it beats passing a null pointer.
      if (a_source.isDefined())
        {
          if (!s_conservativeSource)
            {
              incrementWithSource(a_primMinu[idir], a_source,
                                  0.5*m_dt, modBoxOpen[idir]);
              incrementWithSource(a_primPlus[idir], a_source,
                                  0.5*m_dt, modBoxOpen[idir]);
            }
          else
            {
              int logflag = 0;
              const Box& modBox = modBoxOpen[idir];
              EBCellFAB  consTemp(m_ebisBox, modBox, numConserved());

              primToCons(consTemp,  a_primMinu[idir],  modBox);
              incrementWithSource(consTemp,  a_source, 0.5*m_dt,  modBox); //
              consToPrim(a_primMinu[idir],  consTemp,    modBox, logflag);

              primToCons(consTemp,  a_primPlus[idir],  modBox);
              incrementWithSource(consTemp,  a_source,  0.5*m_dt, modBox); //
              consToPrim(a_primPlus[idir],  consTemp,   modBox, logflag);
            }

        }

      //4. compute estimates of the flux suitable for computing
      // 1D flux derivatives using a riemann solver for the interior R
      //and for the boundary RB.  Equation 1.8
      riemann(a_fluxOne[idir], a_primPlus[idir], a_primMinu[idir],
              idir, faceBox[idir]);

      //some wackiness to work around the fact that fluxfab did
      //not exist when a lot of this was written.
      Vector<EBFaceFAB*> bcfluxes(SpaceDim);
      for (int kdir = 0; kdir < SpaceDim; kdir++)
        {
          bcfluxes[kdir] = &(a_fluxOne[kdir]);
        }
      EBFluxFAB fluxAlias;
      fluxAlias.alias(bcfluxes);

      // Use the user supplied PhysBC object to obtain boundary fluxes
      m_bc->fluxBC(fluxAlias, a_primState, a_primMinu[idir],
                   Side::Lo, m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[idir], idir);
      m_bc->fluxBC(fluxAlias, a_primState, a_primPlus[idir],
                   Side::Hi, m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[idir], idir);

      //extrapolate to covered faces
      extrapToCoveredFaces(m_extendStateNormMinu[idir],
                           a_primMinu[idir],
                           a_primPlus[idir],
                           a_primState,
                           a_coveredFaceNormMinu[idir],
                           idir, Side::Lo, modBoxCov[idir]);

      extrapToCoveredFaces(m_extendStateNormPlus[idir],
                           a_primMinu[idir],
                           a_primPlus[idir],
                           a_primState,
                           a_coveredFaceNormPlus[idir],
                           idir, Side::Hi, modBoxCov[idir]);

      //solve riemann problem between extended state and
      //value in vof to get covered flux  Equation 1.9
      riemann(a_coveredFluxNormMinu[idir],
              m_extendStateNormMinu[idir], a_primMinu[idir],
              a_coveredFaceNormMinu[idir], idir, Side::Lo, modBoxOpen[idir]);
      riemann(a_coveredFluxNormPlus[idir],
              m_extendStateNormPlus[idir], a_primPlus[idir],
              a_coveredFaceNormPlus[idir], idir, Side::Hi, modBoxOpen[idir]);
    }
}

/*******************/
void EBPatchGodunov::
doNormalDerivativeExtr2D(EBCellFAB              a_primMinu[SpaceDim],
                         EBCellFAB              a_primPlus[SpaceDim],
                         EBFaceFAB              a_fluxOne[SpaceDim],
                         BaseIVFAB<Real>        a_coveredFluxNormMinu[SpaceDim],
                         BaseIVFAB<Real>        a_coveredFluxNormPlus[SpaceDim],
                         Vector<VolIndex>       a_coveredFaceNormMinu[SpaceDim],
                         Vector<VolIndex>       a_coveredFaceNormPlus[SpaceDim],
                         EBCellFAB              a_slopesPrim[SpaceDim],
                         EBCellFAB              a_slopesSeco[SpaceDim],
                         const EBCellFAB&       a_flattening,
                         const EBCellFAB&       a_primState,
                         const EBCellFAB&       a_source,
                         const DataIndex&       a_dit,
                         const Box&             a_box)

{
  CH_TIME("EBPatchGodunov::doNormalDerivative2D");
  Box slopeBoxG1 = grow(a_box, 1);
  Box slopeBoxG2 = grow(a_box, 2);
  slopeBoxG1 &= m_domain;
  slopeBoxG2 &= m_domain;
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
      a_slopesPrim[idir].define(m_ebisBox, slopeBoxG2, numSlop);
      a_slopesSeco[idir].define(m_ebisBox, slopeBoxG2, numSlop);

      m_extendStateNormPlus[idir].setVal(0.);
      m_extendStateNormMinu[idir].setVal(0.);
      a_coveredFluxNormPlus[idir].setVal(0.);
      a_coveredFluxNormMinu[idir].setVal(0.);
      a_primPlus[idir].setVal(0.);
      a_primMinu[idir].setVal(0.);
      a_fluxOne[idir].setVal(0.);
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
      slope(a_slopesPrim[idir], a_slopesSeco[idir], a_primState,
            a_flattening, idir, slopeBoxG2, m_useAgg);
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {

      //3. Compute the effect of the normal derivative terms and the
      //source term on the extrapolation in space and time from cell
      //centers to faces. Equation 1.7 (top).
      normalPred(a_primMinu[idir], a_primPlus[idir],
                 a_primState, a_slopesPrim[idir],
                 m_dt/m_dx[idir], idir, slopeBoxG2);

      // If the source term is valid add it to the primitive quantities
      //if source term is not defined, assume it is zero.  lots of
      //problems will have no sources.  not the most elegant method
      //but it beats passing a null pointer.
      if (a_source.isDefined())
        {
          if (!s_conservativeSource)
            {
              incrementWithSource(a_primMinu[idir], a_source,
                                  0.5*m_dt, slopeBoxG2);
              incrementWithSource(a_primPlus[idir], a_source,
                                  0.5*m_dt, slopeBoxG2);
            }
          else
            {
              int logflag = 0;
              const Box& modBox = slopeBoxG2;
              EBCellFAB  consTemp(m_ebisBox, modBox, numConserved());

              primToCons(consTemp,  a_primMinu[idir],  modBox);
              incrementWithSource(consTemp,  a_source, 0.5*m_dt,  modBox); //
              consToPrim(a_primMinu[idir],  consTemp,    modBox, logflag);

              primToCons(consTemp,  a_primPlus[idir],  modBox);
              incrementWithSource(consTemp,  a_source,  0.5*m_dt, modBox); //
              consToPrim(a_primPlus[idir],  consTemp,   modBox, logflag);
            }
        }

      //4. compute estimates of the flux suitable for computing
      // 1D flux derivatives using a riemann solver for the interior R
      //and for the boundary RB.  Equation 1.8
      riemann(a_fluxOne[idir], a_primPlus[idir], a_primMinu[idir],
              idir, faceBox[idir]);

      //some wackiness to work around the fact that fluxfab did
      //not exist when a lot of this was written.
      Vector<EBFaceFAB*> bcfluxes(SpaceDim);
      for (int kdir = 0; kdir < SpaceDim; kdir++)
        {
          bcfluxes[kdir] = &(a_fluxOne[kdir]);
        }
      EBFluxFAB fluxAlias;
      fluxAlias.alias(bcfluxes);

      // Use the user supplied PhysBC object to obtain boundary fluxes
      m_bc->fluxBC(fluxAlias, a_primState, a_primMinu[idir],
                   Side::Lo, m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[idir], idir);
      m_bc->fluxBC(fluxAlias, a_primState, a_primPlus[idir],
                   Side::Hi, m_time, m_ebisBox, a_dit,a_box, bndryFaceBox[idir], idir);

      //extrapolate to covered faces
      extrapToCoveredFaces(m_extendStateNormMinu[idir],
                           a_primMinu[idir],
                           a_primPlus[idir],
                           a_primState,
                           a_coveredFaceNormMinu[idir],
                           idir, Side::Lo, slopeBoxG1);

      extrapToCoveredFaces(m_extendStateNormPlus[idir],
                           a_primMinu[idir],
                           a_primPlus[idir],
                           a_primState,
                           a_coveredFaceNormPlus[idir],
                           idir, Side::Hi, slopeBoxG1);

      //solve riemann problem between extended state and
      //value in vof to get covered flux  Equation 1.9
      riemann(a_coveredFluxNormMinu[idir],
              m_extendStateNormMinu[idir], a_primMinu[idir],
              a_coveredFaceNormMinu[idir], idir, Side::Lo, slopeBoxG2);
      riemann(a_coveredFluxNormPlus[idir],
              m_extendStateNormPlus[idir], a_primPlus[idir],
              a_coveredFaceNormPlus[idir], idir, Side::Hi, slopeBoxG2);
    }
}

/*****************************/
/*****************************/
void
EBPatchGodunov::
computeCoveredFaces(Vector<VolIndex>&     a_coveredFace,
                    IntVectSet&           a_coveredSets,
                    IntVectSet&           a_irregIVS,
                    const int&            a_idir,
                    const Side::LoHiSide& a_sd,
                    const Box&            a_region)
{
  CH_TIME("EBPatchGodunov::computeCoveredFaces");

  EBArith::computeCoveredFaces(a_coveredFace, a_coveredSets,
                               a_irregIVS, a_idir, a_sd,
                               m_ebisBox, a_region);
}
/*****************************/
/*****************************/
void EBPatchGodunov::
getFaceDivergence(EBFluxFAB&        a_openDivU,
                  const EBCellFAB&  a_primState,
                  const EBCellFAB   a_slopePrim[SpaceDim],
                  const Box&        a_box,
                  const IntVectSet& a_ivsIrreg)
{
  CH_TIME("EBPatchGodunov::getFaceDivergence");
  //first set open divergence to zero because
  //we will compute the divergence additively
  CH_assert(a_openDivU.nComp() == 1);
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      a_openDivU[faceDir].setVal(0.);
    }

  //compute divergence on all cells as though they
  //were regular.  Because clarity is not the only virtue
  //in life, we will do this in fortran
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      //this box calc stuff poached directly from
      //the non-eb version.  The fortran used here
      //is slightly different since i am using the
      //existing slopes instead of recalculating the
      //gradients.

      // Now, we can calculate the divergence of the normal velocity
      // at the center normal-direction edges. To do this, we determine
      // which edges at which we have sufficient data to compute centered
      // estimates of h*(div(u)). At the remaining edges. i.e. those
      // corresponding to the physical boundaries, we use zeroth-order
      // extrapolation.
      EBFaceFAB& divUDir = a_openDivU[faceDir];

      Box divBox = divUDir.getRegion();
      divBox.enclosedCells(faceDir);
      divBox.grow(faceDir,1);
      divBox &= m_domain;

      Box loBox,hiBox,centerBox,entireBox;
      int hasLo,hasHi;

      eblohicenterFace(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                     divBox,m_domain,faceDir);

      // All of the boxes computed above are shifted so as to be cell-centered,
      // with the index of the cell center being identified with the low edge.
      // We then shift a_divVel to be compatible with that convention on input,
      // and undo the shift on output.  Nutty though this may seem,
      // it still beats the crap out of rewriting it.

      loBox.shiftHalf(faceDir,1);
      centerBox.shiftHalf(faceDir,1);
      hiBox.shiftHalf(faceDir,1);

      BaseFab<Real>& regOpenDivU = divUDir.getSingleValuedFAB();
      const BaseFab<Real>& regPrimState = a_primState.getSingleValuedFAB();
      Interval velInt = velocityInterval();
      int inormVel = velInt.begin() + faceDir;

      //this does nothing
      regOpenDivU.shiftHalf(faceDir,1);

      //put in normal direction velocity comp
      FORT_DIVUONED(CHF_FRA1(regOpenDivU,0),
                    CHF_CONST_FRA1(regPrimState,inormVel),
                    CHF_CONST_INT(faceDir),
                    CHF_BOX(centerBox));

      //add in slopes
      for (int tranDir = 0; tranDir < SpaceDim; tranDir++)
        {
          if (tranDir != faceDir)
            {
              const BaseFab<Real>& regSlopeDir =  a_slopePrim[tranDir].getSingleValuedFAB();
              Interval velInt = velocityInterval();
              int itranVel = velInt.begin() + tranDir;
              FORT_DIVUTRAN(CHF_FRA1(regOpenDivU,0),
                            CHF_CONST_FRA1(regSlopeDir,itranVel),
                            CHF_CONST_INT(faceDir),
                            CHF_BOX(centerBox));
            }
        }

      //enforce bcs
      FORT_DIVUEDGE(CHF_FRA1(regOpenDivU,0),
                    CHF_CONST_INT(faceDir),
                    CHF_BOX(loBox),
                    CHF_CONST_INT(hasLo),
                    CHF_BOX(hiBox),
                    CHF_CONST_INT(hasHi));

      //this does nothing.
      regOpenDivU.shiftHalf(faceDir,-1);

    }

  //now for the irregular cells
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      EBFaceFAB& openDiv = a_openDivU[faceDir];
      FaceStop::WhichFaces stopCrit;
      if (m_domain.isPeriodic(faceDir))
        {
          stopCrit = FaceStop::SurroundingWithBoundary;
        }
      else
        {
          stopCrit = FaceStop::SurroundingNoBoundary;
        }
      FaceIterator faceit(a_ivsIrreg, m_ebisBox.getEBGraph(), faceDir,
                          stopCrit);

      //reset the divergence to zero
      for (faceit.reset(); faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          openDiv(face, 0) = 0.0;
        }
      //loop through divergence directions and add
      //in diffs
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          //use the normal velocity diff if the direction
          //is the same direction as the face normal.
          //use the ave of velocity slopes otherwise

          Interval velInt = velocityInterval();
          int velInd = velInt.begin() + idir;
          if (idir == faceDir) //use velocity diff
            {
              for (faceit.reset(); faceit.ok(); ++faceit)
                {
                  const FaceIndex& face = faceit();
                  Real normalDiv = 0.0;
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      VolIndex vof = face.getVoF(sit());
                      Real velDir = a_primState(vof, velInd);
                      normalDiv += sign(sit())*velDir;
                    }
                  openDiv(face, 0) += normalDiv;
                }
            }
          else //use ave of slopes
            {
              const EBCellFAB& slopeDir = a_slopePrim[idir];
              for (faceit.reset(); faceit.ok(); ++faceit)
                {
                  const FaceIndex& face = faceit();
                  Real tranDiv = 0.0;
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      VolIndex vof = face.getVoF(sit());
                      Real velSlop = slopeDir(vof, velInd);
                      //the 0.5 is there because we are doing an average
                      tranDiv += 0.5*velSlop;
                    }
                  openDiv(face, 0) += tranDiv;
                }
            }

        }
    }
}
/*****************************/
void EBPatchGodunov::
applyArtificialViscosity(EBFluxFAB&             a_openFlux,
                         BaseIVFAB<Real>        a_coveredFluxMinu[SpaceDim],
                         BaseIVFAB<Real>        a_coveredFluxPlus[SpaceDim],
                         const Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
                         const Vector<VolIndex> a_coveredFacePlus[SpaceDim],
                         const EBCellFAB&       a_consState,
                         const EBFluxFAB&       a_divVel,
                         const Box&             a_box,
                         const IntVectSet&      a_ivsIrreg)
{
  CH_TIME("EBPatchGodunov::applyArtificialViscosity");
  // Get the artificial viscosity coefficient
  Real coeff = artificialViscosityCoefficient();

  Box faceBox[SpaceDim];
  for (int dir1 = 0; dir1 < SpaceDim; ++dir1)
    {
      faceBox[dir1] = a_box;
      faceBox[dir1].grow(dir1,1);
      faceBox[dir1] &= m_domain;
      faceBox[dir1].grow(dir1,-1);
      faceBox[dir1].surroundingNodes(dir1);
    }
  //copy the flux into a temp array so to
  //not screw things up at irregular cells
  //since this is an update in place
  int numCons = numConserved();
  CH_assert(a_openFlux.getRegion() == a_box);
  EBFluxFAB fluxSave(m_ebisBox, a_box, numCons);

  Interval interv(0, numCons-1);
  fluxSave.copy(a_box, interv, a_box, a_openFlux, interv);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      const Box& faceBoxDir = faceBox[idir];
      BaseFab<Real>& regFluxDir = a_openFlux[idir].getSingleValuedFAB();
      const BaseFab<Real>& regConsState = a_consState.getSingleValuedFAB();
      const BaseFab<Real>& regDivVel = a_divVel[idir].getSingleValuedFAB();
      FORT_ARTVISC(CHF_FRA(regFluxDir),
                   CHF_CONST_FRA(regConsState),
                   CHF_CONST_FRA1(regDivVel,0),
                   CHF_CONST_REAL(coeff),
                   CHF_CONST_INT(idir),
                   CHF_BOX(faceBoxDir),
                   CHF_CONST_INT(numCons),
                   CHF_CONST_REAL(m_dx[idir]));
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      EBFaceFAB& fluxDir = a_openFlux[idir];
      const EBFaceFAB& fluxSaveDir = fluxSave[idir];
      const EBFaceFAB& divergeUDir = a_divVel[idir];
      //save the change in the flux so we can
      //do 0th order extrapolation to the covered faces
      EBFaceFAB fluxChange(m_ebisBox, a_box, idir, numCons);
      fluxChange.setVal(0.);
      FaceStop::WhichFaces stopCrit;
      if (m_domain.isPeriodic(idir))
        {
          stopCrit = FaceStop::SurroundingWithBoundary;
        }
      else
        {
          stopCrit = FaceStop::SurroundingNoBoundary;
        }
      for (FaceIterator  faceit(a_ivsIrreg, m_ebisBox.getEBGraph(), idir,
                               stopCrit);
          faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          VolIndex voflo = face.getVoF(Side::Lo);
          VolIndex vofhi = face.getVoF(Side::Hi);
          for (int ivar = 0; ivar < numCons; ivar++)
            {
              Real dv = divergeUDir(face, 0);
              Real uhi= a_consState(vofhi, ivar);
              Real ulo= a_consState(voflo, ivar);

              Real fsave = fluxSaveDir(face, ivar);
              Real deltaF = -coeff*std::max(-dv, (Real)0.)*(uhi-ulo);

              fluxDir(face, ivar) = fsave + deltaF;
              fluxChange(face, ivar) = deltaF;

            }
        }
      //now do the covered face fluxes
      for (SideIterator sit; sit.ok(); ++sit)
        {
          BaseIVFAB<Real>*  coveredFluxPtr;
          const Vector<VolIndex>* coveredFacePtr;
          if (sit() == Side::Lo)
            {
              coveredFacePtr = &a_coveredFaceMinu[idir];
              coveredFluxPtr = &a_coveredFluxMinu[idir];
            }
          else
            {
              coveredFacePtr = &a_coveredFacePlus[idir];
              coveredFluxPtr = &a_coveredFluxPlus[idir];
            }
          BaseIVFAB<Real>&  coveredFlux = *coveredFluxPtr;
          const Vector<VolIndex>& coveredFace = *coveredFacePtr;
          for (int ivof = 0; ivof < coveredFace.size(); ivof++)
            {
              const VolIndex& vof = coveredFace[ivof];
              if (a_box.contains(vof.gridIndex()))
                {
                  Vector<FaceIndex> flipFaces =
                    m_ebisBox.getFaces(vof, idir, flip(sit()));
                  for (int ivar = 0; ivar < numCons; ivar++)
                    {
                      Real areaTot = 0.;
                      Real diffTot = 0.;
                      for (int iface = 0; iface < flipFaces.size(); iface++)
                        {
                          const FaceIndex& face = flipFaces[iface];
                          areaTot += m_ebisBox.areaFrac(face);
                          diffTot += areaTot*fluxChange(face, ivar);
                        }
                      if (areaTot > 0.0)
                        {
                          diffTot /= areaTot;
                        }
                      coveredFlux(vof, ivar) += diffTot;
                    }
                }
            }
        }
    }
}
/*****************************/
void
EBPatchGodunov::
irregularUpdate(EBCellFAB&             a_consState,
                Real&                  a_maxWaveSpeed,
                BaseIVFAB<Real>&       a_massDiff,
                const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                const BaseIVFAB<Real>& a_ebIrregFlux,
                const BaseIVFAB<Real>& a_nonConservativeDivergence,
                const Box&             a_box,
                const IntVectSet&      a_ivs)
{
  CH_TIME("EBPatchGodunov::irregularUpdate");
  CH_assert(isDefined());

  //compute conservative divergences of fluxes
  BaseIVFAB<Real> conservativeDivergence(a_ivs, m_ebisBox.getEBGraph(), numConserved());

  consUndividedDivergence(conservativeDivergence,
                          a_centroidFlux, a_ebIrregFlux, a_ivs);

  // Update conserved variables to the next time step using the final
  // flux in each direction.  As directions progress, increment a_massDiff
  finalUpdate(a_consState, a_massDiff, a_nonConservativeDivergence,
              conservativeDivergence,  a_ivs);

  floorConserved(a_consState, a_box);
  // Get and return the maximum wave speed on this patch/grid
  a_maxWaveSpeed = getMaxWaveSpeed(a_consState, a_box);
}
/*****************************/

/*****************************/
void
EBPatchGodunov::
hybridDivergence(EBCellFAB&             a_hybridDiv,
                 EBCellFAB&             a_consState,
                 BaseIVFAB<Real>&       a_massDiff,
                 const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                 const BaseIVFAB<Real>& a_ebIrregFlux,
                 const BaseIVFAB<Real>& a_nonConsDivF,
                 const Box&             a_box,
                 const IntVectSet&      a_ivs)
{
  CH_TIME("EBPatchGodunov::irregularUpdate");
  CH_assert(isDefined());

  //compute conservative divergences of fluxes
  BaseIVFAB<Real>         conservatDivF(a_ivs, m_ebisBox.getEBGraph(), numConserved());

  consUndividedDivergence(conservatDivF, a_centroidFlux, a_ebIrregFlux, a_ivs);


  //initialize mass difference to zero
  a_massDiff.setVal(0.0);
  int ncons = numConserved();
  //update the irreg vofs
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      //      const IntVect& iv = vof.gridIndex();
      Real volFrac = m_ebisBox.volFrac(vof);

      for (int ivar = 0; ivar < ncons ; ivar++)
        {
          //volume fraction not divided out when divergence
          //calculated.  so we do not multiply here.

          //conservative divergence is already
          //multiplied by volfrac
          Real kapConsDiv = conservatDivF(vof, ivar);
          Real ncDivergeF = a_nonConsDivF(vof, ivar);;

          //dm =-(1-kappa)*kappa*dt*divfcons- kappa*dt*(1-kappa)*divfnc
          //dm =-(1-kappa)*dt*(kappa*divfcons- kappa*divfnc)
          Real deltaM =
            -(1.0-volFrac)*m_dt*(kapConsDiv - volFrac*ncDivergeF);

          a_hybridDiv(vof, ivar) = kapConsDiv + (1.0-volFrac)*ncDivergeF;
          a_massDiff(vof, ivar) = deltaM;
        }
    }
}
/*****************************/
void
EBPatchGodunov::
artificialViscosity(bool a_useArtificialVisc)
{
  m_useArtificialVisc = a_useArtificialVisc;
  m_isArtViscSet = true;
}
/*****************************/

void EBPatchGodunov::
computeFlattening(EBCellFAB&       a_flattening,
                  const EBCellFAB& a_primState,
                  const Box&       a_box)
{
  CH_TIME("EBPatchGodunov::computeFlattening");
  CH_assert(isDefined());
  CH_assert(isDefined());
  CH_assert(a_primState.getRegion().contains(a_box));
  CH_assert(usesFourthOrderSlopes());
  CH_assert(a_flattening.nComp() == 1);
  CH_assert(a_flattening.getRegion().contains(a_box));

  BaseFab<Real>& regFlattening = a_flattening.getSingleValuedFAB();
  const BaseFab<Real>& regPrimState  = a_primState.getSingleValuedFAB();

  const Box& domainBox = m_domain.domainBox();
  EBCellFAB zetaDir(m_ebisBox, a_box, SpaceDim);
  EBCellFAB delta1U(m_ebisBox, a_box, SpaceDim);
  BaseFab<Real>& regZetaDir = zetaDir.getSingleValuedFAB();
  BaseFab<Real>& regDelta1U = delta1U.getSingleValuedFAB();

  Interval velInterval= velocityInterval();
  int v0index = velInterval.begin();

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box box1 = a_box;
      box1.grow(idir, 1);
      box1 &= m_domain;

      Box box2 = a_box;
      box2.grow(idir, 2);
      box2 &= m_domain;

      Box box3 = a_box;
      box3.grow(idir, 3);
      box3 &= m_domain;

      CH_assert(a_primState.getRegion().contains(box3));

      Box loBox, hiBox, centerBox, entireBox;
      int hasLo, hasHi;

      eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                 box3, m_domain, idir);

      EBCellFAB delta1P(m_ebisBox, entireBox, 1);
      BaseFab<Real>& regDelta1P = delta1P.getSingleValuedFAB();

      int pressIndex = pressureIndex();
      //compute delta1P

      FORT_GETGRAD(CHF_FRA1(regDelta1P, 0),
                   CHF_CONST_FRA1(regPrimState, pressIndex),
                   CHF_CONST_INT(idir),
                   CHF_BOX(loBox),
                   CHF_CONST_INT(hasLo),
                   CHF_BOX(hiBox),
                   CHF_CONST_INT(hasHi),
                   CHF_BOX(centerBox));

      //update the irregular vofsn
      for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
        {
          const VolIndex& vof = m_irregVoFs[ivof];
          if (entireBox.contains(vof.gridIndex()))
            {
              const IntVect&  iv = vof.gridIndex();
              //one-sided diffs on domain bndry
              bool onLeftDomain = iv[idir] == domainBox.smallEnd(idir);
              bool onRighDomain = iv[idir] == domainBox.bigEnd(idir);
              bool hasFacesLeft = (m_ebisBox.numFaces(vof, idir, Side::Lo) > 0) && !onLeftDomain;
              bool hasFacesRigh = (m_ebisBox.numFaces(vof, idir, Side::Hi) > 0) && !onRighDomain;

              Real valCent = a_primState(vof, pressIndex);
              Real dpl = 0;
              Real dpr = 0;
              Real dpc = 0;
              //compute one-sided diffs where you have them
              if (hasFacesLeft)
                {
                  Vector<FaceIndex> facesLeft =
                    m_ebisBox.getFaces(vof, idir, Side::Lo);
                  Real valLeft = 0.0;
                  for (int iface = 0; iface <facesLeft.size(); iface++)
                    {
                      VolIndex vofLeft = facesLeft[iface].getVoF(Side::Lo);
                      valLeft += a_primState(vofLeft, pressIndex);
                    }
                  valLeft /= Real(facesLeft.size());
                  dpl = valCent - valLeft;
                }
              if (hasFacesRigh)
                {
                  Vector<FaceIndex> facesRigh =
                    m_ebisBox.getFaces(vof, idir, Side::Hi);
                  Real valRigh = 0.0;
                  for (int iface = 0; iface <facesRigh.size(); iface++)
                    {
                      VolIndex vofRigh = facesRigh[iface].getVoF(Side::Hi);
                      valRigh += a_primState(vofRigh, pressIndex);
                    }
                  valRigh /= Real(facesRigh.size());
                  dpr = valRigh - valCent;
                }
              if (hasFacesLeft && hasFacesRigh)
                {
                  dpc = 0.5*(dpl+dpr);
                }
              else if (!hasFacesLeft && !hasFacesRigh)
                {
                  dpc = 0.0;
                }
              else if (hasFacesLeft && !hasFacesRigh)
                {
                  dpc = dpl;
                }
              else if (hasFacesRigh && !hasFacesLeft)
                {
                  dpc = dpr;
                }

              delta1P(vof, 0) = dpc;
            }
        }

      eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                   box2, m_domain, idir);

      EBCellFAB delta2P(m_ebisBox, entireBox, 1);
      EBCellFAB bulkMin(m_ebisBox, entireBox, 1);
      EBCellFAB zetaTwi(m_ebisBox, entireBox, 1);
      BaseFab<Real>& regDelta2P = delta2P.getSingleValuedFAB();
      BaseFab<Real>& regBulkMin = bulkMin.getSingleValuedFAB();
      BaseFab<Real>& regZetaTwi = zetaTwi.getSingleValuedFAB();

      //compute delta2P

      FORT_GETDPTWO(CHF_FRA1(regDelta2P, 0),
                    CHF_CONST_FRA1(regDelta1P, 0),
                    CHF_CONST_INT(idir),
                    CHF_BOX(loBox),
                    CHF_CONST_INT(hasLo),
                    CHF_BOX(hiBox),
                    CHF_CONST_INT(hasHi),
                    CHF_BOX(centerBox));

      //compute min3(press, d) = bulkMin
      int bulkIndex = bulkModulusIndex();

      FORT_MIN3PTS(CHF_FRA1(regBulkMin, 0),
                   CHF_CONST_FRA1(regPrimState, bulkIndex),
                   CHF_CONST_INT(idir),
                   CHF_BOX(loBox),
                   CHF_CONST_INT(hasLo),
                   CHF_BOX(hiBox),
                   CHF_CONST_INT(hasHi),
                   CHF_BOX(centerBox));

      //compute preliminary flattening coeff

      FORT_GETFLAT(CHF_FRA1(regZetaTwi, 0),
                   CHF_CONST_FRA1(regDelta1P, 0),
                   CHF_CONST_FRA1(regDelta2P, 0),
                   CHF_CONST_FRA1(regBulkMin, 0),
                   CHF_BOX(entireBox));

      //update the irregular vofsn
      for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
        {
          const VolIndex& vof = m_irregVoFs[ivof];
          if (entireBox.contains(vof.gridIndex()))
            {
              const IntVect&  iv = vof.gridIndex();
              //one-sided diffs on domain bndry
              bool onLeftDomain = iv[idir] == domainBox.smallEnd(idir);
              bool onRighDomain = iv[idir] == domainBox.bigEnd(idir);
              bool hasFacesLeft = (m_ebisBox.numFaces(vof, idir, Side::Lo) > 0) && !onLeftDomain;
              bool hasFacesRigh = (m_ebisBox.numFaces(vof, idir, Side::Hi) > 0) && !onRighDomain;

              Real dp2c;
              Real bulkMinL, bulkMinR, bulkMinC;

              Real dp1Cent = delta1P(vof, 0);
              Real bulkCent = a_primState(vof, bulkIndex);
              Real dp1Left=0;
              Real dp1Righ=0;
              Real bulkLeft=0;
              Real bulkRigh=0;
              //compute one-sided diffs where you have them
              if (hasFacesLeft)
                {
                  Vector<FaceIndex> facesLeft =
                    m_ebisBox.getFaces(vof, idir, Side::Lo);
                  dp1Left = 0.0;
                  bulkLeft = 0.0;
                  for (int iface = 0; iface <facesLeft.size(); iface++)
                    {
                      VolIndex vofLeft = facesLeft[iface].getVoF(Side::Lo);
                      dp1Left += delta1P(vofLeft, 0);
                      bulkLeft +=a_primState(vofLeft, bulkIndex);
                    }
                  dp1Left /= Real(facesLeft.size());
                  bulkLeft /= Real(facesLeft.size());

                  bulkMinL = Min(bulkCent, bulkLeft);
                }
              if (hasFacesRigh)
                {
                  Vector<FaceIndex> facesRigh =
                    m_ebisBox.getFaces(vof, idir, Side::Hi);
                  dp1Righ = 0.0;
                  bulkRigh = 0.0;
                  for (int iface = 0; iface <facesRigh.size(); iface++)
                    {
                      VolIndex vofRigh = facesRigh[iface].getVoF(Side::Hi);
                      dp1Righ += delta1P(vofRigh, 0);
                      bulkRigh += a_primState(vofRigh, bulkIndex);
                    }
                  dp1Righ /= Real(facesRigh.size());
                  bulkRigh /= Real(facesRigh.size());
                  bulkMinR = Min(bulkCent, bulkRigh);
                }

              if (hasFacesLeft && hasFacesRigh)
                {
                  dp2c = (dp1Left+dp1Righ);
                  bulkMinC = Min(bulkMinR, bulkMinL);
                }
              else if (!hasFacesLeft && !hasFacesRigh)
                {
                  dp2c = dp1Cent;
                  bulkMinC = bulkCent;
                }
              else if (hasFacesLeft && !hasFacesRigh)
                {
                  dp2c = dp1Cent + dp1Left;
                  bulkMinC = bulkMinL;
                }
              else if (hasFacesRigh && !hasFacesLeft)
                {
                  dp2c = dp1Cent + dp1Righ;
                  bulkMinC = bulkMinR;
                }

              delta2P(vof, 0) = dp2c;
              bulkMin(vof, 0) = bulkMinC;

              //this is that whole zeta(dp1, dp2, p0) func
              Real r0 = 0.75;  Real r1= 0.85; Real d = 0.33;
              Real d1pVoF = Abs(delta1P(vof, 0)) ;
              //
              //     bad idea among many
              Real smallp = 1.0e-7;
              Real d2pVoF = Max(Abs(delta2P(vof, 0)), smallp);

              Real strength = Abs(d1pVoF/bulkMinC);
              Real zetaFuncDP;
              if (strength >= d)
                {
                  Real ratio =  d1pVoF/d2pVoF;
                  if ( ratio <= r0)
                    {
                      zetaFuncDP = 1.0;
                    }
                  else if (ratio >= r1)
                    {
                      zetaFuncDP = 0.0;
                    }
                  else
                    {
                      zetaFuncDP= 1.0 - (ratio - r0)/(r1-r0);;
                    }
                }
              else //strength < d
                {
                  zetaFuncDP = 1.0;
                }
              zetaTwi(vof, 0) = zetaFuncDP;
            }
        }

      eblohicenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                   box1, m_domain, idir);

      //compute zeta = min3(zetatwid) and delta1U

      FORT_MIN3PTS(CHF_FRA1(regZetaDir, idir),
                   CHF_CONST_FRA1(regZetaTwi, 0),
                   CHF_CONST_INT(idir),
                   CHF_BOX(loBox),
                   CHF_CONST_INT(hasLo),
                   CHF_BOX(hiBox),
                   CHF_CONST_INT(hasHi),
                   CHF_BOX(centerBox));

      int velIndex = v0index + idir;

      FORT_GETGRAD(CHF_FRA1(regDelta1U, idir),
                   CHF_CONST_FRA1(regPrimState, velIndex),
                   CHF_CONST_INT(idir),
                   CHF_BOX(loBox),
                   CHF_CONST_INT(hasLo),
                   CHF_BOX(hiBox),
                   CHF_CONST_INT(hasHi),
                   CHF_BOX(centerBox));

      for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
        {
          const VolIndex& vof = m_irregVoFs[ivof];
          if (entireBox.contains(vof.gridIndex()))
            {
              const IntVect&  iv = vof.gridIndex();
              //one-sided diffs on domain bndry
              bool onLeftDomain = iv[idir] == domainBox.smallEnd(idir);
              bool onRighDomain = iv[idir] == domainBox.bigEnd(idir);
              bool hasFacesLeft = (m_ebisBox.numFaces(vof, idir, Side::Lo) > 0) && !onLeftDomain;
              bool hasFacesRigh = (m_ebisBox.numFaces(vof, idir, Side::Hi) > 0) && !onRighDomain;

              Real velLeft=0;
              Real velRigh=0;
              Real du1c=0;
              Real zetaMinL, zetaMinR, zetaMinC;
              Real zetaLeft, zetaRigh;
              Real velCent = a_primState(vof, velIndex);
              Real zetaCent = zetaTwi(vof, 0);
              //compute one-sided diffs where you have them
              if (hasFacesLeft)
                {
                  Vector<FaceIndex> facesLeft =
                    m_ebisBox.getFaces(vof, idir, Side::Lo);
                  velLeft  = 0.0;
                  zetaLeft = 0.0;
                  for (int iface = 0; iface <facesLeft.size(); iface++)
                    {
                      VolIndex vofLeft = facesLeft[iface].getVoF(Side::Lo);
                      velLeft += a_primState(vofLeft, velIndex);
                      zetaLeft +=    zetaTwi(vofLeft, 0);
                    }
                  velLeft /= Real(facesLeft.size());
                  zetaLeft /= Real(facesLeft.size());

                  zetaMinL = Min(zetaCent, zetaLeft);
                }
              if (hasFacesRigh)
                {
                  Vector<FaceIndex> facesRigh =
                    m_ebisBox.getFaces(vof, idir, Side::Hi);
                  velRigh = 0.0;
                  zetaRigh = 0.0;
                  for (int iface = 0; iface <facesRigh.size(); iface++)
                    {
                      VolIndex vofRigh = facesRigh[iface].getVoF(Side::Hi);
                      velRigh += a_primState(vofRigh, velIndex);
                      zetaRigh += zetaTwi(vofRigh, 0);
                    }
                  velRigh /= Real(facesRigh.size());
                  zetaRigh /= Real(facesRigh.size());
                  zetaMinR = Min(zetaCent, zetaRigh);
                }

              if (hasFacesLeft && hasFacesRigh)
                {
                  du1c = 0.5*(velRigh - velLeft);
                  zetaMinC = Min(zetaMinR, zetaMinL);
                }
              else if (!hasFacesLeft && !hasFacesRigh)
                {
                  du1c = 0.0;
                  zetaMinC = zetaCent;
                }
              else if (hasFacesLeft && !hasFacesRigh)
                {
                  du1c = velCent-velLeft;
                  zetaMinC = zetaMinL;
                }
              else if (hasFacesRigh && !hasFacesLeft)
                {
                  du1c = velRigh-velCent;
                  zetaMinC = zetaMinR;
                }

              delta1U(vof, idir) = du1c;
              zetaDir(vof, idir) = zetaMinC;
            }

        }
    }

  // take the minimum of all directions if the sum of velocity diffs
  // are negative.  unity (no flattening) otherwise.

  FORT_MINFLAT(CHF_FRA1(regFlattening, 0),
               CHF_CONST_FRA(regZetaDir),
               CHF_CONST_FRA(regDelta1U),
               CHF_BOX(a_box));

  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          Real velDiffSum = 0.0;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              velDiffSum += delta1U(vof, idir);
            }
          Real finalFlat;
          if (velDiffSum < 0)
            {
              finalFlat = zetaDir(vof, 0);
              for (int idir = 1; idir < SpaceDim; idir++)
                {
                  finalFlat = Min(finalFlat, zetaDir(vof, idir));
                }
              CH_assert(finalFlat >= 0.0);
            }
          else
            {
              finalFlat = 1.0;
            }
          a_flattening(vof, 0) = finalFlat;
        }
    }
}
#include "NamespaceFooter.H"
