#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AdvectDiffuseUtils.H"
#include "ParmParse.H"
#include "CoarseAverage.H"
#include "computeNorm.H"

#include "NamespaceHeader.H"

void ADParseValue(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  ParmParse pp;
  Real bcVal;
  pp.get("bc_value",bcVal);
  a_values[0]=bcVal;
}

void ADParseBC(FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             bool a_homogeneous)
{
  if (!a_domain.domainBox().contains(a_state.box()))
    {

      Box valid = a_valid;
      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))
            {
              ParmParse pp;
              std::vector<int>  bcLo = std::vector<int>();
              std::vector<int>  bcHi = std::vector<int>();
              pp.getarr("bc_lo", bcLo, 0, SpaceDim);
              pp.getarr("bc_hi", bcHi, 0, SpaceDim);
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if (!a_domain.domainBox().contains(ghostBoxLo))
                {
                  if (bcLo[i] == 1)
                    {
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ADParseValue,
                             i,
                             Side::Lo);
                    }
                  else if (bcLo[i] == 0)
                    {
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ADParseValue,
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
                  if (bcHi[i] == 1)
                    {
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ADParseValue,
                             i,
                             Side::Hi);
                    }
                  else if (bcHi[i] == 0)
                    {
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ADParseValue,
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
int
orderScript(int icomp, int inorm, int ncomp)
{
  return icomp + inorm*ncomp;
}

void
getErrorFromCoarseAndFine(Vector< LevelData<FArrayBox>* >&           a_errorCoar,
                          const Vector< LevelData<FArrayBox>* >&     a_solnCoar,
                          const Vector< DisjointBoxLayout >&         a_gridsCoar,
                          const ProblemDomain&                       a_level0DomainCoar,
                          const Vector< LevelData<FArrayBox>* >&     a_solnFine,
                          const Vector< DisjointBoxLayout >&         a_gridsFine,
                          const ProblemDomain&                       a_level0DomainFine,
                          const Vector<int>&                         a_refRat)
{

  int nlevels = a_solnFine.size();
  a_errorCoar.resize(nlevels);
  int nref = 2;  //nothing to do with param refinement ratio. this is the refinement between the two solutions
  int nvar = a_solnFine[0]->nComp();
  Interval interv(0, nvar-1);

  ProblemDomain domLevCoar = a_level0DomainCoar;
  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      int nvar = a_solnCoar[ilev]->nComp();

      a_errorCoar[ilev] = new LevelData<FArrayBox>(a_gridsCoar[ilev], nvar,  IntVect::Zero);

      CoarseAverage averageOp(a_gridsFine[ilev], nvar, nref);
      //here make error = Ave(fine)
      averageOp.averageToCoarse(*a_errorCoar[ilev], *a_solnFine[ilev]);
      //now subtract off coarse so error= Ave(Fine) - coar
      for (DataIterator dit = a_gridsCoar[ilev].dataIterator(); dit.ok(); ++dit)
        {
          FArrayBox&       errorFAB = (*a_errorCoar[ilev])[dit()];
          const FArrayBox& solnFAB  = (*a_solnCoar[ilev])[dit()];

          errorFAB -= solnFAB;
        }
      domLevCoar.refine(a_refRat[ilev]);
    }
}
Real
dtgNorm(const Vector< LevelData<FArrayBox>* >& a_src,
        const Vector< DisjointBoxLayout >&     a_grids,
        const Vector<int>&                     a_refRatio,
        const ProblemDomain&                   a_domain,
        const int& a_comp,
        const int& a_pval)
{
  Real dx = 1;
  Interval interv (a_comp, a_comp);
  int lbase = 0;
  Real volWeightedSum = computeNorm(a_src, a_refRatio, dx, interv, a_pval, lbase);

  //now unweight it from the volume.
  Real norm = volWeightedSum;
  if (a_pval != 0)
    {
      Real volume = 1;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          volume *= a_domain.size(idir);
        }
      norm /= volume;
    }

  return norm;
}
/******/
void
compareError(Vector<Real>&                            a_orders,
             const Vector< LevelData<FArrayBox>* >&   a_errorFine,
             const Vector< LevelData<FArrayBox>* >&   a_errorCoar,
             const Vector< DisjointBoxLayout >&       a_gridsFine,
             const Vector< DisjointBoxLayout >&       a_gridsCoar,
             const Vector<int>&                       a_refRat,
             const ProblemDomain &                    a_coarseDom,
             int a_testverbosity)
{
  ProblemDomain fineDom = a_coarseDom;
  fineDom.refine(2);
  const Vector<int>& refRat = a_refRat;
  const int ncomp = a_errorFine[0]->nComp();
  const int nnorm = 3;
  Real* normsCoar = new Real[ncomp*nnorm];
  Real* normsFine = new Real[ncomp*nnorm];
  a_orders.resize(ncomp*nnorm, 0.0);
  Real* orders    = &(a_orders[0]);
  for (int icomp = 0; icomp < ncomp; icomp++)
    {
      orders[icomp] = 0.0;
      for (int inorm = 0; inorm < nnorm; inorm++)
        {
          normsCoar[orderScript(icomp, inorm, ncomp)] = 0;
          normsFine[orderScript(icomp, inorm, ncomp)] = 0;
        }
    }
  int testverbosity = a_testverbosity;
  for (int comp = 0; comp < ncomp; comp++)
    {
      for (int inorm = 0; inorm <= 2; inorm++)
        {

          Real coarnorm = dtgNorm(a_errorCoar,
                                  a_gridsCoar,
                                  refRat, a_coarseDom,
                                  comp, inorm);
          Real finenorm = dtgNorm(a_errorFine,
                                  a_gridsFine,
                                  refRat, fineDom,
                                  comp, inorm);
          normsCoar[orderScript(comp,inorm,ncomp)] = coarnorm;
          normsFine[orderScript(comp,inorm,ncomp)] = finenorm;
          Real order = 0;
          if ((Abs(finenorm) > 1.0e-10) && (Abs(coarnorm) > 1.0e-10))
            {
              order = log(Abs(coarnorm/finenorm))/log(2.0);
            }
          orders[orderScript(comp,inorm,ncomp)] = order;
        }
    }

  //output in latex format
  int nfine = a_coarseDom.size(0);
  nfine *= 2;
  if (testverbosity > 0)
    {
      pout() << setw(12)
             << setprecision(6)
             << setiosflags(ios::showpoint)
             << setiosflags(ios::scientific) ;
      for (int inorm = 0; inorm <= 2; inorm++)
        {
          pout() << "\\begin{table}[p]" << endl;
          pout() << "\\begin{center}" << endl;
          pout() << "\\begin{tabular}{|c|c|c|c|} \\hline" << endl;
          pout() << "Variable & Coarse Error & Fine Error & Order\\\\" << endl;;
          pout() << "\\hline \\hline " << endl;
          for (int icomp = 0; icomp < ncomp; icomp++)
            {
              int iindex = orderScript(icomp,inorm,ncomp);
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
                     << setprecision(6)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << orders[iindex];
              pout() << " \\\\ " << endl <<   "\\hline"  <<  endl;
            }
          pout() << "\\end{tabular}" << endl;
          pout() << "\\end{center}" << endl;
          pout() << "\\caption{Solution error convergence rates using L-" << inorm << " norm. " << endl;
          pout() << "$h_f = \\frac{1}{" << nfine << "}$ and $h_c = 2 h_f$, $D = " << SpaceDim << "$ }"  << endl;
          pout() << "\\end{table}" << endl;
          pout() << endl << endl;
        }
    }
  //latex output

  delete normsFine;
  delete normsCoar;
}
void
makeFinestDomain(ProblemDomain& a_domain,
                 Real&          a_dx)
{
  ParmParse pp;
  Real domainLength;
  pp.get("domain_length",domainLength);
  getProblemDomain(a_domain);
  int ncells = a_domain.size(0);
  a_dx = domainLength/ncells;
}
///
/**
**/
void
coarsenBoxes(Vector< Vector<Box>      >&    a_boxesCoar,
             const Vector<Vector<Box> >&    a_boxesFine,
             int a_refToCoar)
{
  a_boxesCoar.resize(a_boxesFine.size());
  for (int ilev = 0; ilev < a_boxesFine.size(); ilev++)
    {
      a_boxesCoar[ilev].resize(a_boxesFine[ilev].size());
      for (int ibox = 0; ibox < a_boxesFine[ilev].size(); ibox++)
        {
          a_boxesCoar[ilev][ibox] = coarsen(a_boxesFine[ilev][ibox], a_refToCoar);
        }
    }
}
/*****/
void
getBoxes(Vector<Vector<Box> >&   a_boxes,
         Vector<int>&            a_refRat,
         const Box&              a_domain)
{
  ParmParse pp;
  int maxLevel;
  pp.get("max_level", maxLevel);
  if (maxLevel == 1)
    {
      int amrRef = 2;
      a_refRat.resize(2, amrRef);
      a_boxes.resize(2);
      Box fineBox = refine(a_domain, amrRef);
      int ishrink = fineBox.size(0);
      //this leaves 1/4 refined.
      ishrink *= 3;
      ishrink /= 8;
      fineBox.grow(-ishrink);
      a_boxes[0] = Vector<Box>(1, a_domain);
      a_boxes[1] = Vector<Box>(1, fineBox);
    }
  else if (maxLevel == 0)
    {
      a_refRat.resize(1, 2);
      a_boxes.resize(1);
      a_boxes[0] = Vector<Box>(1, a_domain);
    }
  else
    {
      MayDay::Error("can only deal with two levels now");
    }
}
void
getProblemDomain(ProblemDomain& a_domain)
{
  ParmParse pp;
  Vector<int> ncell(SpaceDim);
  pp.getarr("n_cell", ncell, 0, SpaceDim);
  IntVect hiEnd(D_DECL(ncell[0]-1,ncell[1]-1,ncell[2]-1));
  Box level0Domain(IntVect::Zero, hiEnd);

  Vector<int> v_is_periodic(SpaceDim);
  pp.getarr("periodic_bc", v_is_periodic, 0, SpaceDim);


  bool is_periodic[SpaceDim];
  for (int idir = 0; idir < SpaceDim; idir++) is_periodic[idir] = (v_is_periodic[idir]==1);

  ProblemDomain prob_domain(level0Domain.smallEnd(),
                            level0Domain.bigEnd(),
                            &(is_periodic[0]));

  a_domain = prob_domain;
}

void
getAdvectTestIBC(RefCountedPtr<AdvectTestIBC>& a_ibc)
{
  ParmParse pp;
  Vector<Real> center;
  RealVect centRV;
  Real radius;
  pp.get(   "blob_radius", radius);
  pp.getarr("blob_center", center, 0, SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    centRV[idir] =center[idir];

  a_ibc = RefCountedPtr<AdvectTestIBC>(new AdvectTestIBC(centRV, radius));
}

void
getAdvectionVelocityFunction(AdvectionVelocityFunction& a_velFunc)
{
  ParmParse pp;
  int ivel;
  pp.get("advection_vel_type", ivel);
  if (ivel == 0)
    {
      pout() << "constant advection velocity" << endl;
      a_velFunc = constantAdvection;
    }
  else if (ivel == 1)
    {
      pout() << "rotating advection velocity" << endl;
      a_velFunc = rotatingAdvection;
    }
  else
    {
      MayDay::Error("bogus velocity flag");
    }
}

void
getAMRLADFactory(RefCountedPtr<AMRLevelAdvectDiffuseFactory>&  a_fact,
                 AdvectionVelocityFunction&                    a_velFunc,
                 AdvectPhysics &                               a_advPhys)
{
  ParmParse pp;
  Real cfl = 0.8;
  pp.get("cfl",cfl);

  Real initialCFL = 0.1;
  pp.get("initial_cfl",initialCFL);

  Real domainLength;
  pp.get("domain_length",domainLength);

  Real refineThresh = 0.3;
  pp.get ("refine_thresh",refineThresh);

  int tagBufferSize = 3;
  pp.get("tag_buffer_size",tagBufferSize);

  bool useLimiting = true;
  pp.get("use_limiting", useLimiting);

  Real nu;
  pp.get("diffusion_coef", nu);

  a_fact = RefCountedPtr<AMRLevelAdvectDiffuseFactory>
    (new AMRLevelAdvectDiffuseFactory(a_advPhys, a_velFunc,
                                      ADParseBC, cfl, domainLength,
                                      refineThresh, tagBufferSize,
                                      initialCFL, useLimiting, nu));

}
void
defineAMR(AMR&                                          a_amr,
          RefCountedPtr<AMRLevelAdvectDiffuseFactory>&  a_fact,
          const ProblemDomain&                          a_prob_domain,
          const Vector<int>&                            a_refRat)
{
  ParmParse pp;
  int max_level = 0;
  pp.get("max_level",max_level);

  int num_read_levels = Max(max_level,1);
  std::vector<int> regrid_intervals; // (num_read_levels,1);
  pp.getarr("regrid_interval",regrid_intervals,0,num_read_levels);

  int block_factor = 1;
  pp.get("block_factor",block_factor);

  int max_grid_size = 32;
  pp.get("max_grid_size",max_grid_size);

  Real fill_ratio = 0.75;
  pp.get("fill_ratio",fill_ratio);

  int checkpoint_interval = 0;
  pp.get("checkpoint_interval",checkpoint_interval);

  int plot_interval = 0;
  pp.get("plot_interval",plot_interval);

  Real max_dt_growth = 1.1;
  pp.get("max_dt_growth",max_dt_growth);

  Real dt_tolerance_factor = 1.1;
  pp.get("dt_tolerance_factor",dt_tolerance_factor);
  AMR amr;
  a_amr.define(max_level, a_refRat,
               a_prob_domain,&(*a_fact));

  // set grid generation parameters
  a_amr.maxGridSize(max_grid_size);
  a_amr.blockFactor(block_factor);
  a_amr.fillRatio(fill_ratio);

  // the hyperbolic codes use a grid buffer of 1
  int gridBufferSize;
  pp.get("grid_buffer_size",gridBufferSize);
  a_amr.gridBufferSize(gridBufferSize);

  // set output parameters
  a_amr.checkpointInterval(checkpoint_interval);
  a_amr.plotInterval(plot_interval);
  a_amr.regridIntervals(regrid_intervals);
  a_amr.maxDtGrow(max_dt_growth);
  a_amr.dtToleranceFactor(dt_tolerance_factor);
  if (pp.contains("use_subcycling"))
    {
      bool useSubcycling;
      pp.get("use_subcycling", useSubcycling);
      if (!useSubcycling)
        {
          pout() << "SUBCYCLING IN TIME TURNED OFF!!!"  << endl;
        }
      a_amr.useSubcyclingInTime(useSubcycling);
    }
  if (pp.contains("plot_prefix"))
    {
      std::string prefix;
      pp.get("plot_prefix",prefix);
      a_amr.plotPrefix(prefix);
    }

  if (pp.contains("chk_prefix"))
    {
      std::string prefix;
      pp.get("chk_prefix",prefix);
      a_amr.checkpointPrefix(prefix);
    }

  int verbosity;
  pp.get("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  a_amr.verbosity(verbosity);
}
void
setupAMRForAMRRun(AMR& a_amr)
{
  ParmParse pp;

  if (!pp.contains("restart_file"))
    {
      // initialize from scratch for AMR run
      // initialize hierarchy of levels
      a_amr.setupForNewAMRRun();
    }
  else
    {
      std::string restart_file;
      pp.get("restart_file",restart_file);
      pout() << " restarting from file " << restart_file << endl;

#ifdef CH_USE_HDF5
      HDF5Handle handle(restart_file,HDF5Handle::OPEN_RDONLY);
      // read from checkpoint file
      a_amr.setupForRestart(handle);
      handle.close();
#else
      MayDay::Error("amrGodunov restart only defined with hdf5");
#endif
    }

}

#include "NamespaceFooter.H"
