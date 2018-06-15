#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <fstream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Vector.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "BCFunc.H"
#include "ViscousTensorOp.H"
#include "functionsF_F.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "memusage.H"

#include "UsingNamespace.H"

int s_verbosity = 1;

enum probTypes {zeroRHS = 0,
                unityRHS,
                sinusoidal,
                gaussians,
                numProbTypes};

//int s_probtype = zeroRHS;
//int s_probtype = sinusoidal;
int s_probtype = gaussians;

//  -----------------------------------------
// boundary condition stuff
//  -----------------------------------------
///
/**
 */
class GlobalBCRS
{
public:
  static std::vector<bool> s_printedThatLo, s_printedThatHi;
  static std::vector<int> s_bcLo, s_bcHi;
  static RealVect s_trigvec;
  static bool s_areBCsParsed, s_valueParsed, s_trigParsed;
};

std::vector<bool> GlobalBCRS::s_printedThatLo = std::vector<bool>(SpaceDim, false);
std::vector<bool> GlobalBCRS::s_printedThatHi = std::vector<bool>(SpaceDim, false);
std::vector<int>  GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int>  GlobalBCRS::s_bcHi = std::vector<int>();
RealVect          GlobalBCRS::s_trigvec = RealVect::Zero;
bool              GlobalBCRS::s_areBCsParsed= false;
bool              GlobalBCRS::s_valueParsed= false;
bool              GlobalBCRS::s_trigParsed= false;

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

                  if (GlobalBCRS::s_bcLo[i] == 5)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() != SpaceDim)
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
                          if (a_state.nComp() != SpaceDim)
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
                  if (GlobalBCRS::s_bcHi[i] == 5)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() != SpaceDim)
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
                          if (a_state.nComp() != SpaceDim)
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

// ----------------------------------------
// end BC stuff
// ----------------------------------------

void
setVTcoeffs(Vector<RefCountedPtr<LevelData<FluxBox> > >&  a_vectEta,
            Vector<RefCountedPtr<LevelData<FluxBox> > >&  a_vectLambda,
            Vector<RefCountedPtr<LevelData<FArrayBox> > >&  a_vectAcoef,
            const Vector<DisjointBoxLayout>& a_amrGrids, 
            const Vector<Real> & a_amrDx, 
            int a_finestLevel)
  
{
  int numLevels = a_finestLevel+1;
  a_vectEta.resize(numLevels);
  a_vectLambda.resize(numLevels);
  a_vectAcoef.resize(numLevels);

  for (int lev=0; lev<= a_finestLevel; lev++)
    {
      const DisjointBoxLayout& levelGrids = a_amrGrids[lev];
      a_vectEta[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(levelGrids,1,IntVect::Zero));
      a_vectLambda[lev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(levelGrids,1,IntVect::Zero));
      a_vectAcoef[lev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(levelGrids,1,IntVect::Zero));

      // now set values
      DataIterator dit = levelGrids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& thisA = (*a_vectAcoef[lev])[dit];
          FluxBox& thisEta = (*a_vectEta[lev])[dit];
          FluxBox& thisLambda = (*a_vectLambda[lev])[dit];

          // first do cell-centered A
          
          RealVect cellOffset = (0.5*a_amrDx[lev])*RealVect::Unit;

          BoxIterator cellBit(thisA.box());
          for (cellBit.begin(); cellBit.ok(); ++cellBit)
            {
              IntVect iv = cellBit();
              RealVect loc(iv);
              loc *= a_amrDx[lev];
              loc += cellOffset;

              // loopinv over cells and defining loc is kind of overkill, 
              // since we're only setting A to 1, 
              // but this will make it easier to extend the example
              thisA(iv,0) = 1.0;
            }

          // now do face-centered coefficients
          for (int dir=0; dir<SpaceDim; dir++)
            {
              RealVect faceOffset(cellOffset);
              faceOffset[dir] = 0.0;

              FArrayBox& thisEtaDir = thisEta[dir];
              FArrayBox& thisLambdaDir = thisLambda[dir];

              BoxIterator faceBit(thisEtaDir.box());
              for (faceBit.begin(); faceBit.ok(); ++faceBit)
                {
                  IntVect iv = faceBit();
                  RealVect loc(iv);
                  loc *= a_amrDx[lev];
                  loc += faceOffset;


                  thisEtaDir(iv,0) = 1.0;
                  thisLambdaDir(iv,0) = 1.0;
                } // end loop over faces
            } // end loop over face directions
        } // end loop over patches on this level
    } // end loop over levels

}


void setRHS(Vector<LevelData<FArrayBox>* > a_rhs,
            Vector<ProblemDomain>& a_amrDomains,
            Vector<int>& a_refRatios,
            Vector<Real>& a_amrDx,
            int a_finestLevel)
{

  for (int lev=0; lev<=a_finestLevel; lev++)
    {
      LevelData<FArrayBox>& levelRhs = *(a_rhs[lev]);
      const DisjointBoxLayout& levelGrids = levelRhs.getBoxes();

      // rhs is cell-centered...
      RealVect ccOffset = 0.5*a_amrDx[lev]*RealVect::Unit;

      DataIterator levelDit = levelGrids.dataIterator();
      for (levelDit.begin(); levelDit.ok(); ++levelDit)
        {
          FArrayBox& thisRhs = levelRhs[levelDit];

          if (s_probtype == zeroRHS)
            {
              thisRhs.setVal(0.0);
            }
          else if (s_probtype == unityRHS)
            {
              thisRhs.setVal(1.0);
            }
          else if (s_probtype == sinusoidal)
            {

              BoxIterator bit(thisRhs.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  RealVect loc(iv);
                  loc *= a_amrDx[lev];
                  loc += ccOffset;
                  loc *= 2.0*Pi;


                  thisRhs(iv,0) = D_TERM(sin(loc[0]),
                                         *cos(loc[1]),
                                         *sin(loc[2]));
                  if (SpaceDim > 1)
                    {
                      // can't embed macros 
                      thisRhs(iv,1) = D_TERM(cos(loc[0]),
                                             *sin(loc[1]),
                                             *cos(loc[2]));
                      if (SpaceDim > 2)
                        {
                          thisRhs(iv,2) = D_TERM(sin(loc[0]),
                                                 *sin(loc[1]),
                                                 *sin(loc[2]));
                        }
                    }
                        
                 }

            }
          else if (s_probtype == gaussians)
            {
              int numGaussians = 3;
              Vector<RealVect> center(numGaussians,RealVect::Zero);
              Vector<Real> scale(numGaussians, 1.0);
              Vector<Real> strength(numGaussians, 1.0);

              for (int n=0; n<numGaussians; n++)
                {
                  if (n==0)
                    {
                      strength[0] = 1.0;
                      scale[0] = 1.0e-2;
                      center[0] = 0.25*RealVect::Unit;
                    }
                  else if (n == 1)
                    {
                      strength[1] = 3.0;
                      scale[1] = 1.0e-2;
                      center[1] = RealVect(D_DECL(0.5,0.75, 0.75));
                    }
                  else if (n == 2)
                    {
                      strength[2] = 2.0;
                      scale[2] = 1.0e-2;
                      center[2] = RealVect(D_DECL(0.75,0.5, 0.5));
                    }
                  else
                    {
                      MayDay::Error("too many Gaussian sources attempted");
                    }
                }

              thisRhs.setVal(0.0);

              BoxIterator bit(thisRhs.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  RealVect loc(iv);
                  loc *= a_amrDx[lev];
                  loc += ccOffset;

                  for (int n=0; n<numGaussians; n++)
                    {
                      RealVect dist = loc - center[n];
                      Real radSqr = D_TERM(dist[0]*dist[0],
                                           +dist[1]*dist[1],
                                            +dist[2]*dist[2]);

                       Real val = strength[n]*exp(-radSqr/scale[n]);
                       D_TERM(
                              thisRhs(iv,0) += val;,
                              thisRhs(iv,1) += val;,
                              thisRhs(iv,2) += val;)
                     }
                 }
             }
           else
             {
               MayDay::Error("undefined problem type");
             }
         } // end loop over grids on this level
     } // end loop over levels
 }




void
setupGrids(Vector<DisjointBoxLayout>& a_amrGrids,
           Vector<ProblemDomain>& a_amrDomains,
           Vector<int>& a_refRatios,
           Vector<Real>& a_amrDx,
           int& a_finestLevel)
{
  a_finestLevel = 0;
  ParmParse ppGrids("grids");

  // get grid generation parameters
  int maxLevel, maxBoxSize, blockFactor;
  Real fillRatio;

  ppGrids.get("max_level", maxLevel);

  ppGrids.get("max_box_size",maxBoxSize);

  ppGrids.get("block_factor", blockFactor);

  ppGrids.get("fillRatio", fillRatio);

  // note that there only need to be numLevels-1 refinement ratios
  a_refRatios.resize(maxLevel);
  ppGrids.getarr("ref_ratio", a_refRatios, 0, maxLevel);

  Vector<int>  is_periodic_int;
  bool is_periodic[SpaceDim];
  ppGrids.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      is_periodic[dir] = (is_periodic_int[dir] == 1);
    }

  IntVect numCells;
  Vector<int> incells(SpaceDim);
  ppGrids.getarr("num_cells", incells, 0, SpaceDim);
  numCells = IntVect(D_DECL6(incells[0],incells[1],incells[2],
                             incells[3],incells[4],incells[5]) );

  RealVect domainSize = RealVect::Unit;
  if (ppGrids.contains("domain_size"))
    {
      Vector<Real> insize(SpaceDim);
      ppGrids.getarr("domain_size", insize, 0, SpaceDim);
      domainSize = RealVect(D_DECL6(insize[0],insize[1],insize[2],
                              insize[3],insize[4],insize[5]) );
    }

  // resize dataholders
  int maxNumLevels = maxLevel +1;
  a_amrGrids.resize(maxNumLevels);
  a_amrDomains.resize(maxNumLevels);
  a_amrDx.resize(maxNumLevels,-1);
  a_finestLevel = 0;

  // assumes dx=dy=dz
  a_amrDx[0] = domainSize[0]/numCells[0];

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = numCells - IntVect::Unit;

  ProblemDomain baseDomain(domLo, domHi, is_periodic);
  a_amrDomains[0] = baseDomain;

  // set up refined domains, etc
  for (int lev=1; lev<= maxLevel; lev++)
    {
      a_amrDomains[lev] = a_amrDomains[lev-1];
      a_amrDomains[lev].refine(a_refRatios[lev-1]);
      a_amrDx[lev] = a_amrDx[lev-1]/a_refRatios[lev-1];
    }

  Vector<Vector<Box> > vectBoxes(maxLevel+1);

  // local scope. for base-level grid generation
  {
    // generate base level grids

    domainSplit(baseDomain, vectBoxes[0], maxBoxSize, blockFactor);

    Vector<int> procAssign(vectBoxes[0].size(), 0);

    LoadBalance(procAssign, vectBoxes[0]);

    DisjointBoxLayout baseGrids(vectBoxes[0], procAssign, baseDomain);

    a_amrGrids[0] = baseGrids;
  }


  if (maxLevel > 0)
    {
      bool read_grids = false;
      ppGrids.query("read_in_grids", read_grids);
      if (read_grids)
        {
          for (int ilev = 1; ilev <= maxLevel; ilev++)
            {
              const ProblemDomain& levDomain = a_amrDomains[ilev];

              Vector<Box>   boxes;
              char boxCountVar[100];
              int boxCount;
              sprintf(boxCountVar, "level_%d_box_count", ilev);
              ppGrids.get(boxCountVar, boxCount);
              boxes.resize(boxCount);
              for (int ibox = 0; ibox < boxCount; ibox++)
                {
                  char boxLoVar[100];
                  char boxHiVar[100];
                  sprintf(boxLoVar, "level_%d_box_%d_lo", ilev, ibox);
                  sprintf(boxHiVar, "level_%d_box_%d_hi", ilev, ibox);
                  Vector<int> boxLo, boxHi;
                  ppGrids.getarr(boxLoVar, boxLo, 0, SpaceDim);
                  ppGrids.getarr(boxHiVar, boxHi, 0, SpaceDim);
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
              Vector<int>  proc(boxes.size());
              LoadBalance(proc,boxes);
              a_amrGrids[ilev] = DisjointBoxLayout(boxes, proc, levDomain);
              a_finestLevel++;
            }

        }
      else
        {
          // tag on grad(rhs)
          int bufferSize = 1;
          BRMeshRefine meshGen(a_amrDomains[0],
                               a_refRatios,
                               fillRatio,
                               blockFactor,
                               bufferSize,
                               maxBoxSize);

          // to be used by MeshRefine...
          Vector<Vector<Box> > oldMeshes(maxLevel+1);
          oldMeshes[0] = vectBoxes[0];
          for (int lev=1; lev<oldMeshes.size(); lev++)
            {
              oldMeshes[lev].push_back(a_amrDomains[lev].domainBox());
            }

          Real refineThresh;
          ppGrids.get("refine_threshold", refineThresh);

          Real threshSqr = refineThresh*refineThresh;

          bool moreLevels = true;
          while (moreLevels)
            {
              // tag based on grad(rhs)
              // first need to allocate RHS
              Vector<LevelData<FArrayBox>* > tempRHS(a_finestLevel+1, NULL);
              for (int lev=0; lev<= a_finestLevel; lev++)
                {
                  // note that we add a ghost cell to simplify gradients
                  tempRHS[lev] = new LevelData<FArrayBox>(a_amrGrids[lev],
                                                          SpaceDim,
                                                          IntVect::Unit);
                }

              setRHS(tempRHS, a_amrDomains, a_refRatios, a_amrDx, 
                     a_finestLevel);

              Vector<IntVectSet> tags(a_finestLevel+1);

              for (int lev=0; lev<a_finestLevel+1; lev++)
                {
                  const DisjointBoxLayout& levelGrids = a_amrGrids[lev];
                  const LevelData<FArrayBox>& levelRHS = *tempRHS[lev];
                  IntVectSet& levelTags = tags[lev];

                  // compute mag(gradient)
                  DataIterator dit = levelGrids.dataIterator();
                  for (dit.begin(); dit.ok(); ++dit)
                    {
                      const FArrayBox& rhsFab = levelRHS[dit];
                      // local storage foer gradient
                      FArrayBox gradFab(levelGrids[dit],SpaceDim);
                      gradFab.setVal(0.0);
                      Real thisGrad;

                      BoxIterator bit(levelGrids[dit]);
                      for (bit.begin(); bit.ok(); ++bit)
                        {
                          IntVect iv=bit();
                          for (int dir=0; dir<SpaceDim; dir++)
                            {
                              // use mag(undivided gradient)
                              IntVect hi = iv + BASISV(dir);
                              IntVect lo = iv - BASISV(dir);
                              for (int comp=0; comp<gradFab.nComp(); comp++)
                                {
                                  thisGrad = rhsFab(hi,comp) - rhsFab(lo,comp);
                                  gradFab(iv,comp) += (thisGrad*thisGrad);
                                }
                            } // end loop over directions
                        } // end loop over cells

                      //gradFab now has mag(grad*dx)^2

                      // tag where mag(gradient) > tolerance^2
                      for (bit.begin(); bit.ok(); ++bit)
                        {
                          IntVect iv = bit();
                          for (int comp=0; comp<gradFab.nComp(); comp++)
                            {
                              if (gradFab(iv,comp) > threshSqr)
                                {
                                  levelTags |= iv;
                                }
                            } // end loop over grad comps
                        } // end loop over cells
                    } // end loop over grids on this level

                } // end loop over levels


              // call meshRefine.
              for (int lev=1; lev<=a_finestLevel; lev++)
                {
                  oldMeshes[lev] = vectBoxes[lev];
                }

              int topLevel = a_finestLevel;
              int newFinestLevel =  meshGen.regrid(vectBoxes,
                                                   tags,
                                                   0,
                                                   topLevel,
                                                   oldMeshes);


              // define new grids if necessary and test to see if we're done
              if (newFinestLevel > a_finestLevel)
                {
                  a_finestLevel = newFinestLevel;

                  // setup new grid hierarchy
                  for (int lev=1; lev<=a_finestLevel; lev++)
                    {
                      Vector<int> procAssign(vectBoxes[lev].size(),0);
                      LoadBalance(procAssign, vectBoxes[lev]);
                      DisjointBoxLayout levelGrids(vectBoxes[lev],
                                                   procAssign,
                                                   a_amrDomains[lev]);
                      a_amrGrids[lev] = levelGrids;
                    }
                }
              else
                {
                  moreLevels = false;
                }

              if (a_finestLevel == maxLevel)
                {
                  moreLevels = false;
                }

              // clean up before starting again
              for (int lev=0; lev<tempRHS.size(); lev++)
                {
                  delete tempRHS[lev];
                }

            } // end while (moreLevels)

        }

      // fill in remaining levels with empty DisjointBoxLayouts
      for (int lev= a_finestLevel+1; lev<=maxLevel; lev++)
        {
          a_amrGrids[lev] = DisjointBoxLayout();
        }

    }


}


void
setupSolver(AMRMultiGrid<LevelData<FArrayBox> >& a_amrSolver,
            LinearSolver<LevelData<FArrayBox> >& a_bottomSolver,
            const Vector<DisjointBoxLayout>& a_amrGrids,
            const Vector<ProblemDomain>& a_amrDomains,
            const Vector<int>& a_refRatios,
            const Vector<Real>& a_amrDx,
            int a_finestLevel)
{

  ParmParse ppSolver("solver");

  int numLevels = a_finestLevel+1;

  // solving poisson problem here
  Real alpha =0.0;
  Real beta = 1.0;

  Vector<RefCountedPtr<LevelData<FluxBox> > > vectEta;
  Vector<RefCountedPtr<LevelData<FluxBox> > > vectLambda;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > vectACoef;

  setVTcoeffs(vectEta, vectLambda, vectACoef,
              a_amrGrids, a_amrDx, a_finestLevel);

  ViscousTensorOpFactory opFactory(a_amrGrids,
                                   vectEta,
                                   vectLambda,
                                   vectACoef,
                                   alpha, 
                                   beta,
                                   a_refRatios,
                                   a_amrDomains[0],
                                   a_amrDx[0],
                                   &ParseBC);

  AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) opFactory;

  a_amrSolver.define(a_amrDomains[0], castFact,
                     &a_bottomSolver, numLevels);

  // multigrid solver parameters
  int numSmooth, numMG, maxIter;
  Real eps, hang;
  ppSolver.get("num_smooth", numSmooth);
  ppSolver.get("num_mg",     numMG);
  ppSolver.get("max_iterations", maxIter);
  ppSolver.get("tolerance", eps);
  ppSolver.get("hang",      hang);

  // optional parameters
  ppSolver.query("num_pre", a_amrSolver.m_pre);
  ppSolver.query("num_post", a_amrSolver.m_post);
  ppSolver.query("num_bottom", a_amrSolver.m_bottom);

  Real normThresh = 1.0e-30;
  a_amrSolver.setSolverParameters(numSmooth, numSmooth, numSmooth,
                               numMG, maxIter, eps, hang, normThresh);
  a_amrSolver.m_verbosity = s_verbosity-1;

}

 int runSolver()
 {
   int status = 0;
   ParmParse ppMain("main");

   ppMain.query("verbosity", s_verbosity);

   // set up grids&
   Vector<DisjointBoxLayout> amrGrids;
   Vector<ProblemDomain> amrDomains;
   Vector<int> refRatios;
   Vector<Real> amrDx;
   int finestLevel;

   setupGrids(amrGrids, amrDomains, refRatios, amrDx, finestLevel);

   // initialize solver
   AMRMultiGrid<LevelData<FArrayBox> > amrSolver;
   BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
   bottomSolver.m_verbosity = s_verbosity-2;
   setupSolver(amrSolver, bottomSolver, amrGrids, amrDomains,
               refRatios, amrDx, finestLevel);


   // allocate solution and RHS, initialize RHS
   int numLevels = amrGrids.size();
   Vector<LevelData<FArrayBox>* > phi(numLevels, NULL);
   Vector<LevelData<FArrayBox>* > rhs(numLevels, NULL);
   // this is for convenience
   //Vector<LevelData<FArrayBox>* > resid(numLevels, NULL);

   for (int lev=0; lev<=finestLevel; lev++)
     {
       const DisjointBoxLayout& levelGrids = amrGrids[lev];
       phi[lev] = new LevelData<FArrayBox>(levelGrids, SpaceDim, IntVect::Unit);
       rhs[lev] = new LevelData<FArrayBox>(levelGrids, SpaceDim, IntVect::Zero);
       //resid[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
     }

   setRHS(rhs, amrDomains, refRatios, amrDx,finestLevel);


   // do solve
   int iterations = 1;
   ppMain.get("iterations", iterations);
   for (int iiter = 0; iiter < iterations; iiter++)
     {
       bool zeroInitialGuess = true;
       amrSolver.solve(phi, rhs, finestLevel, 0, zeroInitialGuess);

     }

   // write results to file

   bool writePlots = true;
   ppMain.query("writePlotFiles", writePlots);

#ifdef CH_USE_HDF5
   if (writePlots)
     {
       int numLevels = finestLevel +1;
       Vector<LevelData<FArrayBox>* > plotData(numLevels, NULL);

       for (int lev=0; lev<numLevels; lev++)
         {
           plotData[lev] = new LevelData<FArrayBox>(amrGrids[lev],
                                                    2*SpaceDim, IntVect::Zero);

           Interval phiInterval(0,SpaceDim-1);
           phi[lev]->copyTo(phiInterval, *plotData[lev], phiInterval);
           Interval rhsInterval(SpaceDim,2*SpaceDim-1);
           rhs[lev]->copyTo(phiInterval, *plotData[lev], rhsInterval);
         }

       string fname = "viscousTensorOut.";

       char suffix[30];
       sprintf(suffix, "%dd.hdf5",SpaceDim);
       fname += suffix;

       Vector<string> varNames(2*SpaceDim);
       D_TERM(
              varNames[0] = "phi0";,
              varNames[1] = "phi1";,
              varNames[2] = "phi2";
              )
       D_TERM(
              varNames[SpaceDim] = "rhs0";,
              varNames[SpaceDim+1] = "rhs1";,
              varNames[SpaceDim+2] = "rhs2";
              )

       Real bogusVal = 1.0;

       WriteAMRHierarchyHDF5(fname,
                             amrGrids,
                             plotData,
                             varNames,
                             amrDomains[0].domainBox(),
                             amrDx[0],
                             bogusVal,
                             bogusVal,
                             refRatios,
                             numLevels);

       // clean up
       for (int lev=0; lev<plotData.size(); lev++)
         {
           delete plotData[lev];
         }
     } // end if writing plots
#endif // end if HDF5

   // clean up
   for (int lev=0; lev<phi.size(); lev++)
     {
       delete phi[lev];
       delete rhs[lev];
       //delete resid[lev];
     }

  return status;
}



/*****/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int status = 0;

  // scoping...
  {
    if (argc < 2)
      {
        cerr<< " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }
    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);

    int solverStatus = runSolver();
    status += solverStatus;


  }//end scoping trick

#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif

  return(status);
}
