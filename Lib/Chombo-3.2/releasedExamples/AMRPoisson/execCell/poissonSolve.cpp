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
#include "AMRPoissonOp.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "memusage.H"
#include "CH_Attach.H"

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

void ParseValue(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  //  ParmParse pp;
  //Real bcVal;
  //pp.get("bc_value",bcVal);
  a_values[0]=0.;
}

void ParseBC(FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             bool a_homogeneous)
{

  if (!a_domain.domainBox().contains(a_state.box()))
    {

      // if (!GlobalBCRS::s_areBCsParsed)
      //   {
      //     ParmParse pp;
      //     pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
      //     pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
      //     GlobalBCRS::s_areBCsParsed = true;
      //   }

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
                  // if (GlobalBCRS::s_bcLo[i] == 1)
                  //   {
                  //     if (!GlobalBCRS::s_printedThatLo[i])
                  //       {
                  //         if (a_state.nComp() != 1)
                  //           {
                  //             MayDay::Error("using scalar bc function for vector");
                  //           }
                  //         GlobalBCRS::s_printedThatLo[i] = true;
                  //         if (s_verbosity>5)pout() << "const neum bcs lo for direction " << i << endl;
                  //       }
                  //     NeumBC(a_state,
                  //            valid,
                  //            a_dx,
                  //            a_homogeneous,
                  //            ParseValue,
                  //            i,
                  //            Side::Lo);
                  //   }
                  // else if (GlobalBCRS::s_bcLo[i] == 0)
                  //   {
                  //     if (!GlobalBCRS::s_printedThatLo[i])
                  //       {
                  //         if (a_state.nComp() != 1)
                  //           {
                  //             MayDay::Error("using scalar bc function for vector");
                  //           }
                  //         GlobalBCRS::s_printedThatLo[i] = true;
                  //         if (s_verbosity>5)pout() << "const diri bcs lo for direction " << i << endl;
                  //       }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             true,
                             ParseValue,
                             i,
                             Side::Lo,
                             1);
                  //   }
                  // else
                  //   {
                  //     MayDay::Error("bogus bc flag lo");
                  //   }
                }

              if (!a_domain.domainBox().contains(ghostBoxHi))
                {
                  // if (GlobalBCRS::s_bcHi[i] == 1)
                  //   {
                  //     if (!GlobalBCRS::s_printedThatHi[i])
                  //       {
                  //         if (a_state.nComp() != 1)
                  //           {
                  //             MayDay::Error("using scalar bc function for vector");
                  //           }
                  //         GlobalBCRS::s_printedThatHi[i] = true;
                  //         if (s_verbosity>5)pout() << "const neum bcs hi for direction " << i << endl;
                  //       }
                  //     NeumBC(a_state,
                  //            valid,
                  //            a_dx,
                  //            a_homogeneous,
                  //            ParseValue,
                  //            i,
                  //            Side::Hi);
                  //   }
                  // else if (GlobalBCRS::s_bcHi[i] == 0)
                  //   {
                  //     if (!GlobalBCRS::s_printedThatHi[i])
                  //       {
                  //         if (a_state.nComp() != 1)
                  //           {
                  //             MayDay::Error("using scalar bc function for vector");
                  //           }
                  //         GlobalBCRS::s_printedThatHi[i] = true;
                  //         if (s_verbosity>5)pout() << "const diri bcs hi for direction " << i << endl;
                  //       }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             true,
                             ParseValue,
                             i,
                             Side::Hi,
                             1);
                  //   }
                  // else
                  //   {
                  //     MayDay::Error("bogus bc flag hi");
                  //   }
                }
            } // end if is not periodic in ith direction
        }
    }
}
 
//#ifdef FAS_HACKS
#if 0
////////////////////////////////////////////////////////////////////////
// DampBC helper
////////////////////////////////////////////////////////////////////////
static void DampDiriBC( FArrayBox&      a_state,
                        const Box&      a_valid,
                        const ProblemDomain& a_domain,
                        int             a_ratio,
                        int             a_dir,
                        Side::LoHiSide  a_side,
                        Interval&       a_interval
                        )
{
  Real dc = 1.0, df = dc/(Real)a_ratio;
  for (int kk=1 ; kk <= 2 ; kk++, df += 1., dc += 1. )
    {
      Real fact = df/dc;
      
      Box region = adjCellBox( a_valid, a_dir, a_side, 1 );
      region.shift( a_dir, -kk*sign(a_side) );

      for (BoxIterator bit(region); bit.ok(); ++bit)
        {
          const IntVect& ivTo = bit();
          for (int icomp = a_interval.begin(); icomp <= a_interval.end(); icomp++)
            {
              Real ghostVal = a_state(ivTo, icomp);
              a_state(ivTo, icomp) = fact*ghostVal;
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////
// DampBC -- method to damp (residual) at Dirchlet BCs
////////////////////////////////////////////////////////////////////////
void DampBC( FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             int a_ratio )
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

      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))
            {
              Box valid = a_valid;
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if ( !a_domain.domainBox().contains(ghostBoxLo) && GlobalBCRS::s_bcLo[i]==0 )
                {
                  Interval stateInterval = a_state.interval();
                  DampDiriBC(a_state,
                             valid,
                             a_domain, a_ratio,
                             i,
                             Side::Lo,
                             stateInterval
                             );
                }
              if (!a_domain.domainBox().contains(ghostBoxHi) && GlobalBCRS::s_bcHi[i]==0 )
                {
                  Interval stateInterval = a_state.interval();
                  DampDiriBC(a_state,
                             valid,
                             a_domain, a_ratio,
                             i,
                             Side::Hi,
                             stateInterval
                             );
                }
            } // end if is not periodic in ith direction
        }
    }
}

////////////////////////////////////////////////////////////////////////
// convergeGS_BC helper
////////////////////////////////////////////////////////////////////////
extern void convDiriBC_RBGS( FArrayBox&      a_state,
                             const FArrayBox& a_rhs,
                             const Box&      a_valid,
                             const ProblemDomain& a_domain,
                             Real a_dx,
                             int             a_whichpass,
                             int             a_dir,
                             Side::LoHiSide  a_side);

////////////////////////////////////////////////////////////////////////
// convergeGS_BC -- method to converge G-S on Dirchlet BCs
////////////////////////////////////////////////////////////////////////
void convergeGS_BC( FArrayBox& a_state,
                    const FArrayBox& a_rhs,
                    const Box& a_valid,
                    const ProblemDomain& a_domain,
                    Real a_dx,
                    int a_whichpass )
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

      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))
            {
              Box ghostBoxLo = adjCellBox(a_valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(a_valid, i, Side::Hi, 1);
              if ( !a_domain.domainBox().contains(ghostBoxLo) && GlobalBCRS::s_bcLo[i]==0 )
                {
                  Interval stateInterval = a_state.interval();
                  // iterate 
                  for (int kk=0;kk<10;kk++)
                    {
                      // apply BC
                      DiriBC(a_state,
                             a_valid,
                             a_dx,
                             true,
                             ParseValue,
                             i,
                             Side::Lo,
                             1);
                      // apply G-S to boundary
                      convDiriBC_RBGS(a_state, 
                                      a_rhs,
                                      a_valid,
                                      a_domain, 
                                      a_dx,
                                      a_whichpass,
                                      i,
                                      Side::Lo);
                    }
                }
              if (!a_domain.domainBox().contains(ghostBoxHi) && GlobalBCRS::s_bcHi[i]==0 )
                {
                  Interval stateInterval = a_state.interval();
                  // iterate 
                  for (int kk=0;kk<10;kk++)
                    {
                      // apply BC
                      DiriBC(a_state,
                             a_valid,
                             a_dx,
                             true,
                             ParseValue,
                             i,
                             Side::Hi,
                             1);
                      // apply G-S to boundary
                      convDiriBC_RBGS(a_state, 
                                      a_rhs,
                                      a_valid,
                                      a_domain, 
                                      a_dx,
                                      a_whichpass,
                                      i,
                                      Side::Hi);
                    }
                }
            } // end if is not periodic in ith direction
        }
    }
}
#endif

// ----------------------------------------
// end BC stuff
// ----------------------------------------

void setRHS(Vector<LevelData<FArrayBox>* > a_rhs,
            Vector<ProblemDomain>& a_amrDomains,
            Vector<int>& a_refRatios,
            Vector<Real>& a_amrDx,
            int a_finestLevel)
{
  CH_TIME("setRHS");

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
                       thisRhs(iv,0) += val;
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
  CH_TIME("setupGrids");

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
    CH_TIME("BaseGridCreation");
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
                                                          1, IntVect::Unit);
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
                      FArrayBox gradFab(levelGrids[dit],1);
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
                              thisGrad = rhsFab(hi,0) - rhsFab(lo,0);
                              gradFab(iv,0) += (thisGrad*thisGrad);
                            } // end loop over directions
                        } // end loop over cells

                      //gradFab now has mag(grad*dx)^2

                      // tag where mag(gradient) > tolerance^2
                      for (bit.begin(); bit.ok(); ++bit)
                        {
                          IntVect iv = bit();
                          if (gradFab(iv,0) > threshSqr)
                            {
                              levelTags |= iv;
                            }
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
                  if (s_verbosity>2) pout() << "setupGrids: "<< a_finestLevel <<") size " << a_amrGrids[a_finestLevel].size() << endl;
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
setupSolver(AMRMultiGrid<LevelData<FArrayBox> > *a_amrSolver,
            LinearSolver<LevelData<FArrayBox> >& a_bottomSolver,
            const Vector<DisjointBoxLayout>& a_amrGrids,
            const Vector<ProblemDomain>& a_amrDomains,
            const Vector<int>& a_refRatios,
            const Vector<Real>& a_amrDx,
            int a_finestLevel)
{
  CH_TIME("setupSolver");

  ParmParse ppSolver("solver");

  int numLevels = a_finestLevel+1;

  AMRPoissonOpFactory opFactory;

  // solving poisson problem here
  Real alpha =0.0;
  Real beta = 1.0;

  opFactory.define(a_amrDomains[0],
                   a_amrGrids,
                   a_refRatios,
                   a_amrDx[0],
                   &ParseBC, alpha, beta);

  AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) opFactory;

  a_amrSolver->define(a_amrDomains[0], castFact,
                     &a_bottomSolver, numLevels);

  // multigrid solver parameters
  int numSmooth, numMG, maxIter;
  Real eps, hang;
  ppSolver.get("num_smooth", numSmooth);
  ppSolver.get("num_mg",     numMG);
  ppSolver.get("max_iterations", maxIter);
  ppSolver.get("tolerance", eps);
  ppSolver.get("hang",      hang);

  Real normThresh = 1.0e-30;
  a_amrSolver->setSolverParameters(numSmooth, numSmooth, numSmooth,
                               numMG, maxIter, eps, hang, normThresh);
  a_amrSolver->m_verbosity = s_verbosity-1;

  // optional parameters
  ppSolver.query("num_pre", a_amrSolver->m_pre);
  ppSolver.query("num_post", a_amrSolver->m_post);
  ppSolver.query("num_bottom", a_amrSolver->m_bottom);
}

 int runSolver()
 {
   CH_TIME("runSolver");

   int status = 0, mg_type = 0;
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
   AMRMultiGrid<LevelData<FArrayBox> > *amrSolver;
   if ( mg_type==0 ) 
     {
       amrSolver = new AMRMultiGrid<LevelData<FArrayBox> >();
     }
   else 
     {
       MayDay::Error("FAS not supported");
       // int type = (int)VCYCLE;
       // ParmParse ppSolver("solver");
       
       // AMRFASMultiGrid<LevelData<FArrayBox> > *psol;
       // psol = new AMRFASMultiGrid<LevelData<FArrayBox> >();
       // ppSolver.query("cycle_type", type);
       // psol->setCycleType( (FASMG_type)type );
       // bool avoid_norms = false;
       // ppSolver.query("avoid_norms", avoid_norms);
       // psol->setAvoidNorms(avoid_norms);
       // int numv=1;
       // ppSolver.query("num_v_cycles", numv);
       // psol->setNumVcycles(numv);
       // amrSolver = psol;
     }
   BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
   bottomSolver.m_verbosity = s_verbosity-2;
   setupSolver(amrSolver, bottomSolver, amrGrids, amrDomains,
               refRatios, amrDx, finestLevel);


   // allocate solution and RHS, initialize RHS
   int numLevels = amrGrids.size();
   Vector<LevelData<FArrayBox>* > phi(numLevels, NULL);
   Vector<LevelData<FArrayBox>* > rhs(numLevels, NULL);
   // this is for convenience
   Vector<LevelData<FArrayBox>* > resid(numLevels, NULL);

   for (int lev=0; lev<=finestLevel; lev++)
     {
       const DisjointBoxLayout& levelGrids = amrGrids[lev];
       phi[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);
       rhs[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
       resid[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Zero);
     }

   setRHS(rhs, amrDomains, refRatios, amrDx, finestLevel );

   // do solve
   int iterations = 1;
   ppMain.get("iterations", iterations);

   for (int iiter = 0; iiter < iterations; iiter++)
     {
       bool zeroInitialGuess = true;
       pout() << "about to go into solve" << endl;
       amrSolver->solve(phi, rhs, finestLevel, 0, zeroInitialGuess);
       pout() << "done solve" << endl;
     }

   // write results to file

   bool writePlots = true;
   ppMain.query("writePlotFiles", writePlots);

#ifdef CH_USE_HDF5

   if (writePlots)
     {
       int numLevels = finestLevel +1;
       Vector<LevelData<FArrayBox>* > plotData(numLevels, NULL);
       
       pout() << "Write Plots. norm=" << amrSolver->computeAMRResidual(resid,phi,rhs,finestLevel,0) << endl;
       
       for (int lev=0; lev<numLevels; lev++)
         {
           plotData[lev] = new LevelData<FArrayBox>(amrGrids[lev],
                                                    3, IntVect::Zero);

           Interval phiInterval(0,0);
           phi[lev]->copyTo(phiInterval, *plotData[lev], phiInterval);
           Interval rhsInterval(1,1);
           rhs[lev]->copyTo(phiInterval, *plotData[lev], rhsInterval);
           Interval resInterval(2,2);
           resid[lev]->copyTo(phiInterval, *plotData[lev], resInterval);
         }

       string fname = "poissonOut.";

       char suffix[30];
       sprintf(suffix, "%dd.hdf5",SpaceDim);
       fname += suffix;

       Vector<string> varNames(3);
       varNames[0] = "phi";
       varNames[1] = "rhs";
       varNames[2] = "res";

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
       delete resid[lev];
     }

   delete amrSolver;

   return status;
}


/*****/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int status = 0;
  //  AttachDebugger();
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
  }
  //end scoping trick
  
#ifdef CH_MPI
  dumpmemoryatexit();
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif

  return(status);
}
