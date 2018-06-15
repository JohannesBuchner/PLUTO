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
#include "parstream.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "VClocalFuncs.H"
#include "functionsF_F.H"
#include "AMRIO.H"
#include "BoxIterator.H"
#include "BiCGStabSolver.H"
#include "computeNorm.H"
#include "BCFunc.H"
#include "CoarseAverage.H"

#include "UsingNamespace.H"

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
          Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
          Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
          if (!a_domain.domainBox().contains(ghostBoxLo))
            {
              if (GlobalBCRS::s_bcLo[i] == 1)
                {
                  if (!GlobalBCRS::s_printedThatLo[i])
                    {
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
              else
                {
                  MayDay::Error("bogus bc flag hi");
                }
            }
        }
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

/********/
void getPoissonParameters(VCPoissonParameters&  a_params)
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

  // set to a bogus default value, so we only break from solver
  // default if it's set to something real
  a_params.coefficient_average_type = -1;
  if (pp.contains("coefficient_average_type"))
    {
      std::string tempString;
      pp.get("coefficient_average_type", tempString);
      if (tempString == "arithmetic")
        {
          a_params.coefficient_average_type = CoarseAverage::arithmetic;
        }
      else if (tempString == "harmonic")
        {
          a_params.coefficient_average_type = CoarseAverage::harmonic;
        }
      else
        {
          MayDay::Error("bad coefficient_average_type in input");
        }
    } // end if an average_type is present in inputs

  pp.get("max_level", a_params.maxLevel);
  a_params.numLevels = a_params.maxLevel + 1;
  pp.getarr("ref_ratio",a_params.refRatio,0,a_params.numLevels);
  a_params.verbosity = 3;
  pp.query("verbosity", a_params.verbosity);

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


int setGrids(Vector<DisjointBoxLayout>& vectGrids,
             VCPoissonParameters&         a_params)
{
  Vector<ProblemDomain>     vectDomain;
  Vector<Real>              vectDx;
  getDomainsAndDxes(vectDomain, vectDx, a_params);

  int numlevels = a_params.numLevels;

  ParmParse pp;
  // grid generation parameters

  vectGrids.resize(numlevels);
  bool readInGrids = false;
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

          for (int level=0; level<=topLevel; level++)
            {
              RealVect dxLevel = vectDx[level]*RealVect::Unit;
              setRHS(*vectRHS[level], dxLevel, a_params);
            }

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

/********/
void setRHS(LevelData<FArrayBox>&    a_rhs,
            const RealVect&          a_dx,
            const VCPoissonParameters& a_params)
{
  CH_assert(a_rhs.nComp() == 1);

  int comp = 0;
  const RealVect&  trig = getTrigRV();

  for (DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& thisRHS = a_rhs[dit()];
      Box thisBox = thisRHS.box();
      FORT_GETLOFPHI(CHF_FRA1(thisRHS,comp),
                     CHF_CONST_REALVECT(trig),
                     CHF_CONST_REALVECT(a_dx),
                     CHF_CONST_REALVECT(a_params.probLo),
                     CHF_CONST_REALVECT(a_params.probHi),
                     CHF_CONST_REAL(a_params.alpha),
                     CHF_CONST_REAL(a_params.beta),
                     CHF_BOX(thisBox));

    }
}


extern
AMRLevelOpFactory<LevelData<FArrayBox> >*
defineOperatorFactory(
                      const Vector<DisjointBoxLayout>&             a_grids,
                      const Vector<ProblemDomain>&                 a_vectDomain,
                      Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                      Vector<RefCountedPtr<LevelData<FluxBox> > >& a_bCoef,
                      const VCPoissonParameters&                     a_params)
{
  ParmParse pp2;

  VCAMRPoissonOp2Factory* opFactory = new VCAMRPoissonOp2Factory;

  opFactory->define(a_params.coarsestDomain,
                    a_grids,
                    a_params.refRatio,
                    a_params.coarsestDx,
                    &ParseBC,
                    a_params.alpha,
                    a_aCoef,
                    a_params.beta,
                    a_bCoef);

  if (a_params.coefficient_average_type >= 0)
    {
      opFactory->m_coefficient_average_type
        = a_params.coefficient_average_type;
    }

  return (AMRLevelOpFactory<LevelData<FArrayBox> >*) opFactory;
}

/*
  Set grid hierarchy from input file
 */
void getDomainsAndDxes(  Vector<ProblemDomain>&     vectDomain,
                         Vector<Real>&              vectDx,
                         VCPoissonParameters&         a_params)
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
