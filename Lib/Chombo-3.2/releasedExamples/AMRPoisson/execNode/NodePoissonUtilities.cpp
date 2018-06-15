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

#include "CONSTANTS.H"
#include "parstream.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "LevelData.H"
#include "RealVect.H"
#include "NodePoissonUtilities.H"
#include "functionsF_F.H"
#include "PoissProbF_F.H"
#include "AMRIO.H"
#include "NodeAMRIO.H"
#include "BoxIterator.H"
#include "BiCGStabSolver.H"
#include "computeNorm.H"
#include "BCFunc.H"
#include "NodeBCFunc.H"
#include "AMRNodeOp.H"


std::vector<bool> GlobalBCRS::s_printedThatLo = std::vector<bool>(SpaceDim, false);
std::vector<bool> GlobalBCRS::s_printedThatHi = std::vector<bool>(SpaceDim, false);
std::vector<int>  GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int>  GlobalBCRS::s_bcHi = std::vector<int>();
RealVect          GlobalBCRS::s_trigvec = RealVect::Zero;
bool              GlobalBCRS::s_areBCsParsed= false;
bool              GlobalBCRS::s_valueParsed= false;
bool              GlobalBCRS::s_trigParsed= false;

void
locGetAcoef(LevelData<FArrayBox>& a_acoef,
            const DisjointBoxLayout& a_grids,
            const Real& a_dx)
{
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = a_grids.get(dit());
      for (BoxIterator bit(grid); bit.ok(); ++bit)
        {
          const IntVect& pt = bit();
          RealVect loc(pt);
          loc += 0.5*RealVect::Unit;
          loc *= a_dx;
          Real val = 0;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              Real x = loc[idir];
              val += sin(x)*sin(x);
            }
        }
    }
}


void
BISICLESGetAcoef(LevelData<FArrayBox>& a_acoef,
                 const DisjointBoxLayout& a_grids,
                 const Real& a_dx,
                 const RealVect& a_domLen)
{
  // parameters hardwired for now...
  RealVect omega = RealVect::Unit;
  RealVect twoPiOmega = 2.0*Pi*omega/a_domLen;
  Real magOffset = 0.25;
  Real eps = 0.001;
  Real aVal = 1000.0;

  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = a_grids.get(dit());
      FArrayBox& thisA = a_acoef[dit];
      for (BoxIterator bit(grid); bit.ok(); ++bit)
        {
          const IntVect& pt = bit();
          RealVect loc(pt);
          loc += 0.5*RealVect::Unit;
          loc *= a_dx;
          Real offset = magOffset * sin(twoPiOmega[0]*loc[0]);
          Real val = aVal*(1.0 + eps + sin(twoPiOmega[1]*loc[1] + offset));
          thisA(pt,0) = val;
        }

    }
}

void setBISICLESEta(FArrayBox&      a_eta,
                    const RealVect&          a_dx,
                    const PoissonParameters& a_params,
                    int a_faceDir)
{
  // this is a bit of a kluge -- based on an FFT of the parameters from
  // a BISICLES test run
  Box faceBox = a_eta.box();

  RealVect domainSize = a_params.domainLength;
  int numModes = 64;
  int modeSize = 128;
  Vector<Real> fftA(modeSize,0.0);
  Vector<Real> fftB(modeSize,0.0);

  fftA[0]  =  1.10059014e+13        ;
  fftA[1]  =  -7.30038415e+10       ;
  fftA[2]  =  -3.74986386e+12       ;
  fftA[3]  =  4.18457256e+10        ;
  fftA[4]  =  1.89276504e+12        ;
  fftA[5]  =  5.05355032e+09        ;
  fftA[6]  =  -1.11304968e+12       ;
  fftA[7]  =  -5.14254244e+10       ;
  fftA[8]  =  6.83708396e+11        ;
  fftA[9]  =  9.16484155e+10        ;
  fftA[10] =  -4.31380737e+11       ;
  fftA[11] =  -1.21827704e+11       ;
  fftA[12] =  2.84005884e+11        ;
  fftA[13] =  1.40244626e+11        ;
  fftA[14] =  -2.02288219e+11       ;
  fftA[15] =  -1.46977463e+11       ;
  fftA[16] =  1.61329964e+11        ;
  fftA[17] =  1.43497321e+11        ;
  fftA[18] =  -1.44107633e+11       ;
  fftA[19] =  -1.32233653e+11       ;
  fftA[20] =  1.38756314e+11        ;
  fftA[21] =  1.16115103e+11        ;
  fftA[22] =  -1.37251092e+11       ;
  fftA[23] =  -9.81218723e+10       ;
  fftA[24] =  1.34612917e+11        ;
  fftA[25] =  8.08917333e+10        ;
  fftA[26] =  -1.28292069e+11       ;
  fftA[27] =  -6.64211523e+10       ;
  fftA[28] =  1.17600133e+11        ;
  fftA[29] =  5.58888127e+10        ;
  fftA[30] =  -1.03154523e+11       ;
  fftA[31] =  -4.96140264e+10       ;
  fftA[32] =  8.63480100e+10        ;
  fftA[33] =  4.71416791e+10        ;
  fftA[34] =  -6.88736391e+10       ;
  fftA[35] =  -4.74318691e+10       ;
  fftA[36] =  5.23373356e+10        ;
  fftA[37] =  4.91146881e+10        ;
  fftA[38] =  -3.79847994e+10       ;
  fftA[39] =  -5.07675639e+10       ;
  fftA[40] =  2.65523152e+10        ;
  fftA[41] =  5.11703682e+10        ;
  fftA[42] =  -1.82373366e+10       ;
  fftA[43] =  -4.94998512e+10       ;
  fftA[44] =  1.27695239e+10        ;
  fftA[45] =  4.54386817e+10        ;
  fftA[46] =  -9.55077407e+09       ;
  fftA[47] =  -3.91861688e+10       ;
  fftA[48] =  7.82921588e+09        ;
  fftA[49] =  3.13776153e+10        ;
  fftA[50] =  -6.86881112e+09       ;
  fftA[51] =  -2.29307581e+10       ;
  fftA[52] =  6.08626264e+09        ;
  fftA[53] =  1.48506983e+10        ;
  fftA[54] =  -5.13160689e+09       ;
  fftA[55] =  -8.03099509e+09       ;
  fftA[56] =  3.90745203e+09        ;
  fftA[57] =  3.08351894e+09        ;
  fftA[58] =  -2.52938614e+09       ;
  fftA[59] =  -2.28274977e+08       ;
  fftA[60] =  1.24470188e+09        ;
  fftA[61] =  -7.39570342e+08       ;
  fftA[62] =  -3.30870630e+08       ;
  fftA[63] =  4.01811318e+08        ;
  //          fftA[64] =  1.00000000e+04        ;

  fftB[0]  =  0.00000000e+00     ;
  fftB[1]  =  2.97384561e+12     ;
  fftB[2]  =  -1.84219001e+11    ;
  fftB[3]  =  -5.67288982e+11    ;
  fftB[4]  =  1.86421083e+11     ;
  fftB[5]  =  -4.09732194e+10    ;
  fftB[6]  =  -1.65105324e+11    ;
  fftB[7]  =  2.96372001e+11     ;
  fftB[8]  =  1.35998058e+11     ;
  fftB[9]  =  -4.08128910e+11    ;
  fftB[10] =  -1.08055243e+11    ;
  fftB[11] =  4.40228782e+11     ;
  fftB[12] =  8.61522516e+10     ;
  fftB[13] =  -4.24526524e+11    ;
  fftB[14] =  -7.23798748e+10    ;
  fftB[15] =  3.81024746e+11     ;
  fftB[16] =  6.68250609e+10     ;
  fftB[17] =  -3.23724788e+11    ;
  fftB[18] =  -6.81578418e+10    ;
  fftB[19] =  2.62703118e+11     ;
  fftB[20] =  7.41667778e+10     ;
  fftB[21] =  -2.04971857e+11    ;
  fftB[22] =  -8.22651642e+10    ;
  fftB[23] =  1.54951337e+11     ;
  fftB[24] =  8.99454472e+10     ;
  fftB[25] =  -1.14857461e+11    ;
  fftB[26] =  -9.51479083e+10    ;
  fftB[27] =  8.51108240e+10     ;
  fftB[28] =  9.65119197e+10     ;
  fftB[29] =  -6.47906775e+10    ;
  fftB[30] =  -9.34938124e+10    ;
  fftB[31] =  5.21112317e+10     ;
  fftB[32] =  8.63480200e+10     ;
  fftB[33] =  -4.48826203e+10    ;
  fftB[34] =  -7.59903328e+10    ;
  fftB[35] =  4.09149992e+10     ;
  fftB[36] =  6.37732454e+10     ;
  fftB[37] =  -3.83294918e+10    ;
  fftB[38] =  -5.12165650e+10    ;
  fftB[39] =  3.57545418e+10     ;
  fftB[40] =  3.97383395e+10     ;
  fftB[41] =  -3.24032839e+10    ;
  fftB[42] =  -3.04271861e+10    ;
  fftB[43] =  2.80413205e+10     ;
  fftB[44] =  2.38900529e+10     ;
  fftB[45] =  -2.28719031e+10    ;
  fftB[46] =  -2.01934721e+10    ;
  fftB[47] =  1.73700441e+10     ;
  fftB[48] =  1.89013409e+10     ;
  fftB[49] =  -1.21036732e+10    ;
  fftB[50] =  -1.91970914e+10    ;
  fftB[51] =  7.57530204e+09     ;
  fftB[52] =  2.00636823e+10     ;
  fftB[53] =  -4.10974751e+09    ;
  fftB[54] =  -2.04865464e+10    ;
  fftB[55] =  1.80340591e+09     ;
  fftB[56] =  1.96440702e+10     ;
  fftB[57] =  -5.35061174e+08    ;
  fftB[58] =  -1.70517642e+10    ;
  fftB[59] =  2.81347269e+07     ;
  fftB[60] =  1.26376039e+10     ;
  fftB[61] =  5.45388811e+07     ;
  fftB[62] =  -6.73520941e+09    ;
  fftB[63] =  -9.86948912e+06    ;

  fftA[0] = fftA[0]/128;
  for (int k=1; k<fftA.size(); k++)
    {
      fftA[k] = 2.0*fftA[k]/128;
      fftB[k] = -2.0*fftB[k]/128;
    }

  // parameters hardwired for now...
  RealVect omega = RealVect::Unit;
  //RealVect twoPiOmega = 2.0*Pi*omega/domainSize;
  RealVect twoPiOmega = omega;
  Real magOffset = 0.1;
  //Real eps = 0.001;
  //Real aVal = 1000.0;

  RealVect offsetVect = RealVect::Unit;
  offsetVect[a_faceDir] = 0.0;
  BoxIterator bit(faceBox);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect loc(RealVect::Unit);
      loc *= iv;
      loc += 0.5*offsetVect;
      loc *= a_dx;
      loc /= domainSize;
      // this bit of silliness is due to the fact that
      // visit reversed the order of the coefficients when
      // I wrote them to a file.
      loc[1] = 1.0-loc[1];
      loc *= 2.0*Pi;

      Real offset = magOffset*sin(twoPiOmega[0]*loc[0]);
      //offset = 0.0;
      loc[1] -= offset;

      // this was Mark's original eta
      //Real etaVal = 1.3e11 + 1.2e11*sin(2*Pi*pow(loc[0],8))*sin(0.5*loc[1]);
      Real etaVal = fftA[0];
      Real maxEta = 2.3863e11;

      for (int k=1; k<numModes; k++)
        {
          etaVal += fftA[k]*cos(k*loc[1]) + fftB[k]*sin(k*loc[1]);
        }
      a_eta(iv,0) = min(maxEta, etaVal);
    }

}



/***************/
RealVect& getTrigRV()
{
  Real pi = 4.*atan(1.0);
  if (!GlobalBCRS::s_trigParsed)
    {
      GlobalBCRS::s_trigParsed = true;
      ParmParse pp;
      std::vector<Real> trigvec(SpaceDim);
      pp.getarr("trig",trigvec,0,SpaceDim);

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          GlobalBCRS::s_trigvec[idir] = pi*trigvec[idir];
        }
    }
  return GlobalBCRS::s_trigvec;
}
void TrigValueNeum(Real* pos,
                   int* dir,
                   Side::LoHiSide* side,
                   Real* a_values)
{
  RealVect& trig = getTrigRV();
  RealVect xval;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xval[idir] = pos[idir];
    }
  RealVect gradPhi;
  FORT_GETGRADPHIPOINT(CHF_REALVECT(gradPhi),
                       CHF_CONST_REALVECT(trig),
                       CHF_CONST_REALVECT(xval));

  a_values[0] = gradPhi[*dir];
}

void ResistDiri(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  RealVect& trig = getTrigRV();
  RealVect xval;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xval[idir] = pos[idir];
    }
  ParmParse pp;
  int whichMag;
  pp.get("which_mag", whichMag);
  for (int icomp = 0; icomp < 3; icomp++)
    {
      Real value;
      FORT_GETMAGPOINTRESIST(CHF_REAL(value),
                             CHF_CONST_REALVECT(trig),
                             CHF_CONST_REALVECT(xval),
                             CHF_INT(icomp),
                             CHF_INT(whichMag));
      a_values[icomp] = value;
    }
}

void ViscousDiri(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  RealVect& trig = getTrigRV();
  RealVect xval;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xval[idir] = pos[idir];
    }
  ParmParse pp;
  int whichMag;
  pp.get("which_mag", whichMag);
  for (int icomp = 0; icomp < SpaceDim; icomp++)
    {
      Real value;
      FORT_GETMAGPOINTRESIST(CHF_REAL(value),
                             CHF_CONST_REALVECT(trig),
                             CHF_CONST_REALVECT(xval),
                             CHF_INT(icomp),
                             CHF_INT(whichMag));
      a_values[icomp] = value;
    }
}


void TrigValueDiri(Real* pos,
                   int* dir,
                   Side::LoHiSide* side,
                   Real* a_values)
{
  RealVect& trig = getTrigRV();
  RealVect xval;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xval[idir] = pos[idir];
    }
  Real value;
  FORT_GETPHIPOINT(CHF_REAL(value),
                   CHF_CONST_REALVECT(trig),
                   CHF_CONST_REALVECT(xval));
  a_values[0] = value;
}

void
nodeOutputData(const Vector<LevelData<NodeFArrayBox>* >& vectPhi,
               const Vector<DisjointBoxLayout>& vectGrids,
               const ProblemDomain&     domain,
               const Vector<int>& vectRatio,
               Real dxCoarsest,
               int numlevels,
               string filename,
               string varname)
{
  // for now, make a placeholders for these -- if more than
  // one comp, probably want to change vectNames accordingly
  Vector<string> vectNames(vectPhi[0]->nComp(), varname);
  Real dt = 1.0;
  Real time = 0.0;

#ifdef CH_USE_HDF5
  WriteAMRHierarchyHDF5(filename,
                        vectGrids,
                        vectPhi,
                        vectNames,
                        domain.domainBox(),
                        dxCoarsest,
                        dt, time,
                        vectRatio,
                        numlevels);
#endif

}

void
outputData(const Vector<LevelData<FArrayBox>* >& vectPhi,
           const Vector<DisjointBoxLayout>& vectGrids,
           const ProblemDomain&     domain,
           const Vector<int>& vectRatio,
           Real dxCoarsest,
           int numlevels,
           string filename,
           string varname)
{
#ifdef CH_USE_HDF5
  Vector<string> vectName(vectPhi[0]->nComp(), varname);
  Real time = 1;  //placeholder
  WriteAMRHierarchyHDF5(filename,
                        vectGrids,
                        vectPhi,
                        vectName,
                        domain.domainBox(),
                        dxCoarsest, time, time,
                        vectRatio,
                        numlevels);
#endif
}


void
outputVectorData(const Vector<LevelData<FArrayBox>* >& vectPhi,
                 const Vector<DisjointBoxLayout>& vectGrids,
                 const ProblemDomain&     domain,
                 const Vector<int>& vectRatio,
                 Real dxCoarsest,
                 int numlevels,
                 string filename,
                 Vector<string>& a_vectName)
{
#ifdef CH_USE_HDF5
  Real time = 1;  //placeholder
  WriteAMRHierarchyHDF5(filename,
                        vectGrids,
                        vectPhi,
                        a_vectName,
                        domain.domainBox(),
                        dxCoarsest, time, time,
                        vectRatio,
                        numlevels);
#endif
}

int setRHS(Vector<LevelData<FArrayBox>* >& vectRhs,
           PoissonParameters&              a_params,
           int  a_numLevels)

{
  int numlevels = a_numLevels;
  if (numlevels < 0)
    {
      numlevels = a_params.numLevels;
    }
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  Real rhono, rno;
  int iprob;
  ParmParse pp;
  pp.get("iprob", iprob);

  rhono = 0.75;
  rno = 0.5;

  Real dxlev = a_params.coarsestDx;
  ProblemDomain domlev = a_params.coarsestDomain;
  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      if (iprob == -1)
        {
          //trig prob for convergence tests
          RealVect vectDx = dxlev*RealVect::Unit;
          setTrigKappaLOfPhi(*vectRhs[ilev], vectDx, a_params);
        }
      else
        {
          LevelData<FArrayBox>& rhsLD = *vectRhs[ilev];
          DataIterator dit =  rhsLD.dataIterator();
          for (dit.reset(); dit.ok(); ++dit)
            {
              FArrayBox& rhsFab = rhsLD[dit()];
              rhsFab.setVal(0.);
              bool useEBGrids;
              pp.get("use_eb_grids", useEBGrids);

              if (useEBGrids)
                {
                  rhsFab.setVal(1.0);
                }
              else
                {
                  FORTNT_GETRHSPOIS(CHF_FRA(rhsFab),
                                  CHF_BOX(rhsFab.box()),
                                  CHF_BOX(domlev.domainBox()),
                                  CHF_CONST_REAL(dxlev),
                                  CHF_CONST_REAL(rhono),
                                  CHF_CONST_REAL(rno),
                                  CHF_CONST_INT(iprob));
                }
            }
        }
      dxlev /= a_params.refRatio[ilev];
      domlev.refine(a_params.refRatio[ilev]);

    }
  return 0;
}

/*
  tag cells for refinement based on magnitude(RHS)
*/
void
tagCells(Vector<LevelData<FArrayBox>* >& vectRHS,
         Vector<IntVectSet>& tagVect,
         Vector<Real>& vectDx,
         Vector<ProblemDomain>& vectDomain,
         const Real refine_thresh,
         const int tags_grow,
         const int baseLevel,
         int numLevels)
{
  for (int lev=baseLevel; lev!= numLevels; lev++)
    {
      IntVectSet local_tags;
      LevelData<FArrayBox> & levelRhs = *vectRHS[lev];
      DisjointBoxLayout level_domain = levelRhs.getBoxes();
      DataIterator dit = levelRhs.dataIterator();

      Real maxRHS = 0;

      maxRHS = norm(levelRhs, levelRhs.interval(), 0);

      Real tagVal = maxRHS * refine_thresh;

      // now loop through grids and tag cells where RHS > tagVal
      for (dit.reset(); dit.ok(); ++dit)
        {
          const Box thisBox = level_domain.get(dit());
          const FArrayBox& thisRhs = levelRhs[dit()];
          BoxIterator bit(thisBox);
          for (bit.begin(); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if (abs(thisRhs(iv)) >= tagVal)
                local_tags |= iv;
            }
        } // end loop over grids on this level

      local_tags.grow(tags_grow);
      const Box& domainBox = vectDomain[lev].domainBox();
      local_tags &= domainBox;

      tagVect[lev] = local_tags;

    } // end loop over levels
}

/*
  Set grid hierarchy from input file
*/
void getDomainsAndDxes(  Vector<ProblemDomain>&     vectDomain,
                         Vector<Real>&              vectDx,
                         PoissonParameters&         a_params)
{

  vectDomain.resize(a_params.numLevels);
  vectDx.resize(    a_params.numLevels);
  vectDx[0] = a_params.coarsestDx;
  for (int ilev = 1; ilev < a_params.numLevels; ilev++)
    {
      vectDx[ilev] = vectDx[ilev-1]/a_params.refRatio[ilev-1];
    }


  vectDomain[0] = a_params.coarsestDomain;
  for (int ilev = 1;ilev < a_params.numLevels; ilev++)
    {
      vectDomain[ilev] = refine(vectDomain[ilev-1],a_params.refRatio[ilev-1]);
    }
}

int setGrids(Vector<DisjointBoxLayout>& vectGrids,
             PoissonParameters&         a_params)
{
  Vector<ProblemDomain>     vectDomain;
  Vector<Real>              vectDx;
  getDomainsAndDxes(vectDomain, vectDx, a_params);

  int numlevels = a_params.numLevels;

  ParmParse pp;
  bool useEBGrids;
  // grid generation parameters

  vectGrids.resize(numlevels);
  bool readInGrids = false;
  pp.query("use_eb_grids", useEBGrids);
  if (pp.contains("read_in_grids"))
    {
      pp.get("read_in_grids", readInGrids);
    }


  if (readInGrids)
    {

      ProblemDomain levDomain = a_params.coarsestDomain;
      for (int ilev = 0; ilev < a_params.numLevels; ilev++)
        {
          Vector<Box>   boxes;
          char boxCountVar[100];
          int boxCount;
          sprintf(boxCountVar, "level_%d_box_count", ilev);
          pp.get(boxCountVar, boxCount);
          boxes.resize(boxCount);
          for (int ibox = 0; ibox < boxCount; ibox++)
            {
              char boxLoVar[100];
              char boxHiVar[100];
              sprintf(boxLoVar, "level_%d_box_%d_lo", ilev, ibox);
              sprintf(boxHiVar, "level_%d_box_%d_hi", ilev, ibox);
              Vector<int> boxLo, boxHi;
              pp.getarr(boxLoVar, boxLo, 0, SpaceDim);
              pp.getarr(boxHiVar, boxHi, 0, SpaceDim);
              IntVect ivLo(D_DECL(boxLo[0], boxLo[1], boxLo[2]));
              IntVect ivHi(D_DECL(boxHi[0], boxHi[1], boxHi[2]));
              boxes[ibox] = Box(ivLo, ivHi);
              if (!levDomain.contains(boxes[ibox]))
                {
                  MayDay::Error("box outside of domain");
                }
            }
          //check to see if level 0 domain is covered
          if (ilev == 0)
            {
              IntVectSet ivDom(levDomain.domainBox());
              for (int ibox = 0; ibox < boxes.size(); ibox++)
                {
                  ivDom -= boxes[ibox];
                }
              if (!ivDom.isEmpty())
                {
                  MayDay::Error("level 0 boxes must cover the domain");
                }
            }
          Vector<int>  proc(a_params.numLevels);
          LoadBalance(proc,boxes);
          vectGrids[ilev] = DisjointBoxLayout(boxes, proc, levDomain);
          levDomain.refine(a_params.refRatio[ilev]);
        }

    }
  else if (useEBGrids)
    {
      pout() << "all regular geometry" << endl;
      pout() << "ignoring grid parameters and making simple  grids" << endl;
      Vector<int> proc(1, 0);
      Box coarBox = vectDomain[0].domainBox();
      Vector<Box> coarBoxes(1, coarBox);
      vectGrids[0] = DisjointBoxLayout(coarBoxes, proc, vectDomain[0]);

      for (int ilev = 1; ilev < numlevels; ilev++)
        {
          int iboxShrink = coarBox.size(0);
          iboxShrink /= 4;
          if (iboxShrink < 2)
            {
              MayDay::Error("wacky DBL generation technique failed, try making base box bigger");
            }
          coarBox.grow(-iboxShrink);
          coarBox.refine(a_params.refRatio[ilev-1]);
          Vector<Box> refBoxes(1, coarBox);
          vectGrids[ilev] = DisjointBoxLayout(refBoxes, proc,
                                              vectDomain[ilev]);
        }
    }
  else
    {
      pout() << "tagging on gradient of RHS" << endl;
      int maxLevel = numlevels-1;
      Vector<Vector<Box> > newBoxes(numlevels);
      Vector<Vector<Box> > oldBoxes(numlevels);

      // determine grids dynamically, based on grad(RHS)
      // will need temp storage for RHS
      Vector<LevelData<FArrayBox>* > vectRHS(maxLevel+1,NULL);
      int ncomps = 1;

      // define base level first
      Vector< Vector<int> > procAssign(maxLevel+1);
      domainSplit(vectDomain[0], oldBoxes[0], a_params.maxGridSize, a_params.blockFactor);
      procAssign[0].resize(oldBoxes[0].size());
      LoadBalance(procAssign[0],oldBoxes[0]);

      vectGrids[0].define(oldBoxes[0],procAssign[0],vectDomain[0]);

      vectRHS[0] = new LevelData<FArrayBox>(vectGrids[0], ncomps,
                                            IntVect::Zero);

      int topLevel = 0;

      bool moreLevels = (maxLevel > 0);

      int nesting_radius = 2;
      // create grid generation object
      BRMeshRefine meshrefine(vectDomain[0], a_params.refRatio,
                              a_params.fillRatio,
                              a_params.blockFactor, nesting_radius,
                              a_params.maxGridSize);

      while (moreLevels)
        {
          // default is moreLevels = false
          // (only repeat loop in the case where a new level
          // is generated which is still less than maxLevel)
          moreLevels = false;

          int baseLevel = 0;
          int oldTopLevel = topLevel;

          // now initialize RHS for this existing hierarchy
          setRHS(vectRHS, a_params, topLevel+1);

          Vector<IntVectSet> tagVect(topLevel+1);
          int tags_grow = 1;
          tagCells(vectRHS, tagVect, vectDx, vectDomain,
                   a_params.refineThresh,
                   tags_grow, baseLevel, topLevel+1);

          int new_finest = meshrefine.regrid(newBoxes, tagVect,
                                             baseLevel,
                                             topLevel, oldBoxes);

          if (new_finest > topLevel)
            {
              topLevel++;
            }

          oldBoxes = newBoxes;

          //  no need to do this for the base level (already done)
          for (int lev=1; lev<= topLevel; lev++)
            {
              // do load balancing
              procAssign[lev].resize(newBoxes[lev].size());
              LoadBalance(procAssign[lev], newBoxes[lev]);
              const DisjointBoxLayout newDBL(newBoxes[lev], procAssign[lev],
                                             vectDomain[lev]);
              vectGrids[lev] = newDBL;
              delete vectRHS[lev];
              vectRHS[lev] = new LevelData<FArrayBox>(vectGrids[lev], ncomps,
                                                      IntVect::Zero);
            } // end loop over levels for initialization

          // figure out whether we need another pass through grid generation
          if ((topLevel<maxLevel) && (topLevel > oldTopLevel))
            moreLevels = true;

        } // end while moreLevels loop
      // clean up temp storage
      for (int ilev=0; ilev <vectRHS.size(); ilev++)
        {
          if (vectRHS[ilev] != NULL)
            {
              delete vectRHS[ilev];
              vectRHS[ilev] = NULL;
            }
        }
    }
  return 0;
}




void ParseValue(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  ParmParse pp;
  Real bcVal;
  pp.get("bc_value",bcVal);
  a_values[0]=bcVal;
}


void ParseBC(FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             bool a_homogeneous)
{
  if (!a_domain.domainBox().contains(a_state.box()))
    {

      if (!GlobalBCRS::s_areBCsParsed)
        {
          ParmParse pp;
          pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
          pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
          GlobalBCRS::s_areBCsParsed = true;
        }

      Box valid = a_valid;
      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))
            {
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if (!a_domain.domainBox().contains(ghostBoxLo))
                {
                  if (GlobalBCRS::s_bcLo[i] == 1)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "const neum bcs lo for direction " << i << endl;
                        }
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Lo);
                    }
                  else if (GlobalBCRS::s_bcLo[i] == 2)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "trig neum bcs lo for direction " << i << endl;
                        }
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             TrigValueNeum,
                             i,
                             Side::Lo);
                    }
                  else if (GlobalBCRS::s_bcLo[i] == 0)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "const diri bcs lo for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Lo);

                    }
                  else if (GlobalBCRS::s_bcLo[i] == 3)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "trig diri bcs lo for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             TrigValueDiri,
                             i,
                             Side::Lo);

                    }
                  else if (GlobalBCRS::s_bcLo[i] == 4)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "periodic bcs lo for direction " << i << endl;
                        }

                    }
                  else if (GlobalBCRS::s_bcLo[i] == 5)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() == 1)
                            {
                              MayDay::Error("using vector bc function for scalar");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "reflective slip bcs lo for direction " << i << endl;
                        }
                      ReflectiveVectorBC(a_state,
                                         valid,
                                         a_dx,
                                         i,
                                         Side::Lo);
                    }
                  else if (GlobalBCRS::s_bcLo[i] == 6)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() == 1)
                            {
                              MayDay::Error("using vector bc function for scalar");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "no slip bcs lo for direction " << i << endl;
                        }
                      NoSlipVectorBC(a_state,
                                     valid,
                                     a_dx,
                                     i,
                                     Side::Lo, 1);
                    }
                  else if (GlobalBCRS::s_bcLo[i] == 7)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() != 3)
                            {
                              MayDay::Error("bc function hardwired for ncomp == 3");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "Resistive diri bcs lo for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ResistDiri,
                             i,
                             Side::Lo);

                    }
                  else if (GlobalBCRS::s_bcLo[i] == 8)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() != SpaceDim)
                            {
                              MayDay::Error("bc function hardwired for ncomp == SpaceDim");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "Viscous diri bcs lo for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ViscousDiri,
                             i,
                             Side::Lo);

                    }
                  else
                    {
                      MayDay::Error("bogus bc flag lo");
                    }
                }

              if (!a_domain.domainBox().contains(ghostBoxHi))
                {
                  if (GlobalBCRS::s_bcHi[i] == 1)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "const neum bcs hi for direction " << i << endl;
                        }
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 2)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "trig neum bcs hi for direction " << i << endl;
                        }
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             TrigValueNeum,
                             i,
                             Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 0)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "const diri bcs hi for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 3)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "trig diri bcs hi for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             TrigValueDiri,
                             i,
                             Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 4)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "periodic bcs hi for direction " << i << endl;
                        }
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 5)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() == 1)
                            {
                              MayDay::Error("using vector bc function for scalar");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "reflective slip bcs hi for direction " << i << endl;
                        }
                      ReflectiveVectorBC(a_state,
                                         valid,
                                         a_dx,
                                         i,
                                         Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 6)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() == 1)
                            {
                              MayDay::Error("using vector bc function for scalar");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "no slip bcs hi for direction " << i << endl;
                        }
                      NoSlipVectorBC(a_state,
                                     valid,
                                     a_dx,
                                     i,
                                     Side::Hi, 1);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 7)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() != 3)
                            {
                              MayDay::Error("this bc hardwired for  ncomp=3");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "resisitive mhd diri bcs hi for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ResistDiri,
                             i,
                             Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 8)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() != SpaceDim)
                            {
                              MayDay::Error("this bc hardwired for spacedim ncomp");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "resisitive mhd diri bcs hi for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ViscousDiri,
                             i,
                             Side::Hi);
                    }
                  else
                    {
                      MayDay::Error("bogus bc flag hi");
                    }
                }
            } // end if is not periodic in ith direction
        }
    }
}

void NodeParseBC(NodeFArrayBox& a_state,
                 const Box& a_valid,
                 const ProblemDomain& a_domain,
                 Real a_dx,
                 bool a_homogeneous)
{
  if (!a_domain.domainBox().contains(a_state.box()))
    {

      if (!GlobalBCRS::s_areBCsParsed)
        {
          ParmParse pp;
          pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
          pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
          GlobalBCRS::s_areBCsParsed = true;
        }

      Box valid = a_valid;
      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          if (!a_domain.isPeriodic(i))
            {
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if (!a_domain.domainBox().contains(ghostBoxLo))
                {
                  if (GlobalBCRS::s_bcLo[i] == 1)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "const neum bcs lo for direction "
                                 << i << endl;
                        }
                      NodeNeumBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 ParseValue,
                                 i,
                                 Side::Lo);
                    }
                  else if (GlobalBCRS::s_bcLo[i] == 2)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "trig neum bcs lo for direction "
                                 << i << endl;
                        }
                      NodeNeumBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 TrigValueNeum,
                                 i,
                                 Side::Lo);
                    }
                  else if (GlobalBCRS::s_bcLo[i] == 0)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "const diri bcs lo for direction "
                                 << i << endl;
                        }
                      NodeDiriBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 ParseValue,
                                 i,
                                 Side::Lo);

                    }
                  else if (GlobalBCRS::s_bcLo[i] == 3)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "trig diri bcs lo for direction "
                                 << i << endl;
                        }
                      NodeDiriBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 TrigValueDiri,
                                 i,
                                 Side::Lo);

                    }
                  else
                    {
                      MayDay::Error("bogus bc flag lo");
                    }
                }

              if (!a_domain.domainBox().contains(ghostBoxHi))
                {
                  if (GlobalBCRS::s_bcHi[i] == 1)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "const neum bcs hi for direction "
                                 << i << endl;
                        }
                      NodeNeumBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 ParseValue,
                                 i,
                                 Side::Hi);
                    }
                  if (GlobalBCRS::s_bcHi[i] == 2)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "trig neum bcs hi for direction "
                                 << i << endl;
                        }
                      NodeNeumBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 TrigValueNeum,
                                 i,
                                 Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 0)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "const diri bcs hi for direction "
                                 << i << endl;
                        }
                      NodeDiriBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 ParseValue,
                                 i,
                                 Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 3)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "trig diri bcs hi for direction "
                                 << i << endl; //
                        }
                      NodeDiriBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 TrigValueDiri,
                                 i,
                                 Side::Hi);
                    }
                  else
                    {
                      MayDay::Error("bogus bc flag hi");
                    }
                }
            } // end if periodic in ith direction
        }

    }

}


void
nodeDefineSolver(AMRMultiGrid<LevelData<NodeFArrayBox> >&         a_solver,
                 const Vector<DisjointBoxLayout>&                 a_grids,
                 LinearSolver<LevelData<NodeFArrayBox> >&         a_bottomSolver,
                 const PoissonParameters&                         a_params)
{
  ParmParse pp2;
  AMRNodeOpFactory opFactory;
  opFactory.define(a_params.coarsestDomain,
                   a_grids,
                   a_params.refRatio,
                   a_params.coarsestDx,
                   &NodeParseBC, a_params.alpha, a_params.beta);

  ProblemDomain coarsestDomain(a_params.coarsestDomain);
  a_solver.define(coarsestDomain, opFactory,  &a_bottomSolver, a_params.numLevels);

  int numSmooth, numMG, maxIter;
  Real eps, hang;
  pp2.get("num_smooth", numSmooth);
  pp2.get("num_mg",     numMG);
  pp2.get("max_iterations", maxIter);
  pp2.get("tolerance", eps);
  pp2.get("hang",      hang);
  Real normThresh = 1.0e-30;
  a_solver.setSolverParameters(numSmooth, numSmooth, numSmooth,
                               numMG, maxIter, eps, hang, normThresh);
  a_solver.m_verbosity = 3;

}


void setKLVViscous(LevelData<FArrayBox>&    a_klv,
                   const RealVect&          a_dx,
                   const PoissonParameters& a_params)
{
  int whichVel, whichLambda, whichEta,  whichBeta;
  Real eps;
  ParmParse pp;
  pp.get("which_vel", whichVel);
  pp.get("which_eta", whichEta);
  pp.get("which_lambda", whichLambda);
  pp.get("which_beta", whichBeta);
  pp.get("eta_eps",   eps);

  //which_lambda == 2 forces lambda = -factor*eta
  Real factor = 1;
  if (whichLambda == 2)
    {
      pp.get("lambda_factor", factor);
    }
  for (DataIterator dit = a_klv.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& curFAB = a_klv[dit()];
      Box curBox = curFAB.box();

      const RealVect&     trig = getTrigRV();

      for (int comp = 0; comp < a_klv.nComp(); comp++)
        {
          FORT_GETKLVVISCOUS(CHF_FRA1(curFAB,comp),
                             CHF_CONST_REALVECT(trig),
                             CHF_CONST_REALVECT(a_dx),
                             CHF_CONST_REALVECT(a_params.probLo),
                             CHF_CONST_REAL(a_params.alpha),
                             CHF_CONST_REAL(a_params.beta),
                             CHF_BOX(curBox),
                             CHF_INT(comp),
                             CHF_REAL(eps),
                             CHF_INT(whichVel),
                             CHF_INT(whichEta),
                             CHF_INT(whichLambda),
                             CHF_REAL(factor)
                             );
        }
    }

}

/********/
void setEtaResistive(LevelData<FluxBox>&      a_eta,
                     const RealVect&          a_dx,
                     const PoissonParameters& a_params,
                     int a_whichEta)
{
  CH_assert(a_eta.nComp() == 1);

  int whichEta;
  ParmParse pp;
  if (a_whichEta < 0)
    {
      pp.get("which_eta", whichEta);
    }
  else
    {
      whichEta = a_whichEta;
    }
  Real eps;
  pp.get("eta_eps", eps);

  int ibox = 0;
  for (DataIterator dit = a_eta.dataIterator(); dit.ok(); ++dit)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FArrayBox& curFAB = a_eta[dit()][idir];

          if (whichEta == 0)
            {
              curFAB.setVal(eps);
            }
          else if (whichEta == 2)
            {
              // BISICLES test case
              setBISICLESEta(curFAB, a_dx, a_params, idir);

              // eta(x,y) = 1.3e11 + 1.2e11*sin(2*pi*(x/L)8) * sin(pi*(y/L));

            }
          else
            {

              Box curBox = curFAB.box();
              const RealVect&     trig = getTrigRV();

              FORT_GETETARESIST(CHF_FRA1(curFAB, 0),
                                CHF_CONST_REALVECT(trig),
                                CHF_CONST_REALVECT(a_dx),
                                CHF_CONST_REALVECT(a_params.probLo),
                                CHF_BOX(curBox),
                                CHF_INT(idir),
                                CHF_REAL(eps),
                                CHF_INT(whichEta));
              ibox++;
            }
        }
    }
}

/********/
void setTrigPhi(LevelData<FArrayBox>&    a_phi,
                const RealVect&          a_dx,
                const PoissonParameters& a_params)
{
  CH_assert(a_phi.nComp() == 1);

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& curPhiFAB = a_phi[dit()];
      Box curPhiBox = curPhiFAB.box();

      const RealVect&     trig = getTrigRV();

      for (int comp = 0; comp < a_phi.nComp(); comp++)
        {
          FORT_GETPHI(CHF_FRA1(curPhiFAB,comp),
                      CHF_CONST_REALVECT(trig),
                      CHF_CONST_REALVECT(a_dx),
                      CHF_CONST_REALVECT(a_params.probLo),
                      CHF_CONST_REALVECT(a_params.probHi),
                      CHF_BOX(curPhiBox));
        }
    }
}

/********/
void setTrigKappaLOfPhi(LevelData<FArrayBox>&    a_kappaLOfPhi,
                        const RealVect&          a_dx,
                        const PoissonParameters& a_params)
{
  CH_assert(a_kappaLOfPhi.nComp() == 1);

  const RealVect&  trig = getTrigRV();

  for (DataIterator dit = a_kappaLOfPhi.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& curKappaLOfPhiFAB = a_kappaLOfPhi[dit()];
      Box curKappaLOfPhiBox = curKappaLOfPhiFAB.box();
      for (int comp = 0; comp < a_kappaLOfPhi.nComp(); comp++)
        {
          FORT_GETLOFPHI(CHF_FRA1(curKappaLOfPhiFAB,comp),
                         CHF_CONST_REALVECT(trig),
                         CHF_CONST_REALVECT(a_dx),
                         CHF_CONST_REALVECT(a_params.probLo),
                         CHF_CONST_REALVECT(a_params.probHi),
                         CHF_CONST_REAL(a_params.alpha),
                         CHF_CONST_REAL(a_params.beta),
                         CHF_BOX(curKappaLOfPhiBox));
        }

    }
}
/******/
int iscript(int icomp, int inorm, int ncomp)
{
  return icomp + inorm*ncomp;
}
/******/
void compareError(const Vector< LevelData<FArrayBox>* >&   a_errorFine,
                  const Vector< LevelData<FArrayBox>* >&   a_errorCoar,
                  const Vector< DisjointBoxLayout >&       a_gridsFine,
                  const Vector< DisjointBoxLayout >&       a_gridsCoar,
                  const PoissonParameters&                 a_paramsFine,
                  const PoissonParameters&                 a_paramsCoar,
                  const string& a_testName)
{
  const Vector<int> refRat = a_paramsCoar.refRatio;
  const int ncomp = a_errorFine[0]->nComp();
  const int nnorm = 3;
  Real* normsCoar = new Real[ncomp*nnorm];
  Real* normsFine = new Real[ncomp*nnorm];
  Real* orders    = new Real[ncomp*nnorm];
  for (int icomp = 0; icomp < ncomp; icomp++)
    {
      orders[icomp] = 0.0;
      for (int inorm = 0; inorm < nnorm; inorm++)
        {
          normsCoar[iscript(icomp, inorm, ncomp)] = 0;
          normsFine[iscript(icomp, inorm, ncomp)] = 0;
        }
    }
  ParmParse pp;
  //pout() << "==============================================" << endl;
  for (int comp = 0; comp < ncomp; comp++)
    {
      //pout() << "Comparing error in variable  " << comp << endl;
      //pout() << "==============================================" << endl;
      for (int inorm = 0; inorm <= 2; inorm++)
        {

          if (inorm == 0)
            {
              //pout() << endl << "Using max norm." << endl;
            }
          else
            {
              //pout() << endl << "Using L-" << inorm << "norm." << endl;
            }
          Real dxCoar = a_paramsCoar.coarsestDx;
          Real dxFine = a_paramsFine.coarsestDx;
          Interval comps(comp,comp);
          int lbase = 0;
          Real coarnorm = computeNorm(a_errorCoar, refRat, dxCoar, comps, inorm, lbase);

          Real finenorm = computeNorm(a_errorFine, refRat, dxFine, comps, inorm, lbase);

          //pout() << "Coarse Error Norm = " << coarnorm << endl;
          //pout() << "Fine   Error Norm = " << finenorm << endl;

          normsCoar[iscript(comp,inorm,ncomp)] = coarnorm;
          normsFine[iscript(comp,inorm,ncomp)] = finenorm;

          if ((Abs(finenorm) > 1.0e-10) && (Abs(coarnorm) > 1.0e-10))
            {
              Real order = log(Abs(coarnorm/finenorm))/log(2.0);
              ////pout() << "Order of scheme = " << order << endl;
              orders[iscript(comp,inorm,ncomp)] = order;
            }
        }
      //pout() << "==============================================" << endl ;;
    }


  //output in latex format to be safe
//  int nfine = a_paramsFine.coarsestDomain.size(0);
//  pout() << setw(12)
//         << setprecision(6)
//         << setiosflags(ios::showpoint)
//         << setiosflags(ios::scientific) ;

  for (int inorm = 0; inorm <= 2; inorm++)
    {
//      pout() << "\\begin{table}[p]" << endl;
//      pout() << "\\begin{center}" << endl;
//      pout() << "\\begin{tabular}{|c|c|c|c|} \\hline" << endl;
//      pout() << "Variable & Coarse Error & Fine Error & Order\\\\" << endl;;
//      pout() << "\\hline \\hline " << endl;
      for (int icomp = 0; icomp < ncomp; icomp++)
        {
          int iindex = iscript(icomp,inorm,ncomp);
          pout() << "var" << icomp << " &    \t "
                 << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << normsCoar[iindex]  << " & "
                 << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << normsFine[iindex] << " & "
                 << setw(12)
                 << setprecision(2)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << orders[iindex];
          pout() << endl;
          //          pout() << " \\\\ " << endl <<   "\\hline"  <<  endl;
        }
//      pout() << "\\end{tabular}" << endl;
//      pout() << "\\end{center}" << endl;
//      pout() << "\\caption{";
//      pout() << a_testName ;
//      pout() << " convergence rates using L-" << inorm << " norm. " << endl;
//      pout() << "$h_f = \\frac{1}{" << nfine << "}$ and $h_c = 2 h_f$, $D = " << SpaceDim << "$ }"  << endl;
//      pout() << "\\end{table}" << endl;
//      pout() << endl << endl;
    }
  pout() << "data for paper: ";
  int nx = a_paramsCoar.coarsestDomain.domainBox().size(0);
  pout() << nx << " ";
  for (int inorm = 0; inorm <= 2; inorm++)
    {
      int icomp = 0;
      int iindex = iscript(icomp,inorm,ncomp);
      pout()     << setw(12)
                 << setprecision(6)
                 << normsCoar[iindex]  << "  ";
    }
  pout () << endl;
  delete[] normsCoar;
  delete[] normsFine;
  delete[] orders   ;
}
/********/
void getCoarseLayoutsFromFine(Vector<DisjointBoxLayout>&       a_gridsCoar,
                              const Vector<DisjointBoxLayout>& a_gridsFine,
                              const PoissonParameters&         a_paramsCoar)
{
  int nlevels = a_paramsCoar.numLevels;
  a_gridsCoar.resize(nlevels);
  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      CH_assert(a_gridsFine[ilev].coarsenable(2));
      coarsen(a_gridsCoar[ilev], a_gridsFine[ilev], 2);
    }

}
/********/
void PoissonParameters::coarsen(int a_factor)
{
  coarsestDx *= a_factor;
  coarsestDomain.coarsen(a_factor);
}
/********/
void PoissonParameters::refine(int a_factor)
{
  coarsestDx /= a_factor;
  coarsestDomain.refine(a_factor);
}
/********/
void getPoissonParameters(PoissonParameters&  a_params)
{
  ParmParse pp;

  std::vector<int> nCellsArray(SpaceDim);
  pp.getarr("n_cells",nCellsArray,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.nCells[idir] = nCellsArray[idir];
    }

  Vector<int> is_periodic(SpaceDim, false);
  pp.queryarr("periodic", is_periodic, 0, SpaceDim);

  pp.get("refine_threshold",a_params.refineThresh);
  pp.get("block_factor",a_params.blockFactor);
  pp.get("fill_ratio",a_params.fillRatio);
  pp.get("buffer_size",a_params.bufferSize);
  pp.get("alpha",a_params.alpha);
  pp.get("beta", a_params.beta);

  pp.get("verbose", a_params.verbosity);

  pp.get("max_level", a_params.maxLevel);
  a_params.numLevels = a_params.maxLevel + 1;
  pp.getarr("ref_ratio",a_params.refRatio,0,a_params.numLevels);

  IntVect lo = IntVect::Zero;
  IntVect hi = a_params.nCells;
  hi -= IntVect::Unit;

  Box crseDomBox(lo,hi);
  ProblemDomain crseDom(crseDomBox);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      crseDom.setPeriodic(dir, is_periodic[dir]);
    }
  a_params.coarsestDomain = crseDom;

  std::vector<Real> dLArray(SpaceDim);
  pp.getarr("domain_length",dLArray,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.domainLength[idir] = dLArray[idir];
    }

  pp.get("max_grid_size",a_params.maxGridSize);

  //derived stuff
  a_params.coarsestDx = a_params.domainLength[0]/a_params.nCells[0];

  a_params.probLo = RealVect::Zero;
  a_params.probHi = RealVect::Zero;
  a_params.probHi += a_params.domainLength;


}

