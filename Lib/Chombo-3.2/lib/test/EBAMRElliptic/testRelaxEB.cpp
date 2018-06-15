#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SphereIF.H"
#include "PlaneIF.H"


#include "EBISLayout.H"
#include "EBAMRIO.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
//#include "DebugDump.H"
#include "CH_Attach.H"

#include "BiCGStabSolver.H"
#include "DirichletPoissonDomainBC.H"
#include "DirichletPoissonEBBC.H"
#include "RelaxSolver.H"
#include "EBCellFactory.H"
#include "EBAMRPoissonOp.H"
#include "EBAMRPoissonOpFactory.H"
#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testRelaxEB" ;
static const char* indent = "";
static const char* indent2 = "" ;


class DirichletBC: public BaseBCValue
{
public:
  Real value (const RealVect &a_point,
              const RealVect &a_normal,
              const Real&     a_time,
              const int&      a_comp) const;
};

class RHS: public BaseBCValue
{
public:
  Real value (const RealVect &a_point,
              const RealVect &a_normal,
              const Real&     a_time,
              const int&      a_comp) const;
};


Real DirichletBC::value (const RealVect &a_point,
                         const RealVect &a_normal,
                         const Real&     a_time,
                         const int&      a_comp) const
{
  //  return 0.0;
  return 3*a_point[0]*a_point[0] + 2*a_point[1]+1 ;
  //return a_point[0]*2.0/3.0 + a_point[1]*1.5;
  //return 2*a_point[0];
}


Real RHS::value (const RealVect &a_point,
                 const RealVect &a_normal,
                 const Real&     a_time,
                 const int&      a_comp) const
{
  return 6;
  //return 0;

}


int readGeometryInfo(Box& a_domain,
                     RealVect& a_dx,
                     RealVect& a_origin,
                     RealVect& a_center,
                     Real& a_radius);


void makeHierarchy(Vector<DisjointBoxLayout>& dbl,
                   const ProblemDomain& baseDomain,
                   const IntVectSet&baseTags);



void swapvector(Vector<LevelData<EBCellFAB>* >& left,
                Vector<LevelData<EBCellFAB>* >& right);

void setValue(LevelData<EBCellFAB>& phase, const BaseBCValue& value,
              const Box& domain,
              const RealVect& dx, const RealVect& origin,
              bool  useKappa = false);

/***************/
// define a sphere EBIS.
/***************/
int makeGeometry(
                 const Box& domain,
                 const RealVect& dx,
                 const RealVect& origin,
                 const RealVect& center,
                 const Real& radius);
/***************/

/***************/
/***************/
void dumpmemoryatexit();
/***************/
/***************/

char iter_str[80];

int
main(int argc,char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // registerDebugger();

  // begin forever present scoping trick
  {
    pout()<<std::endl;

    Vector<std::string> names0(1, "phi");
    Vector<int> refRatio(3,2);
    Vector<Real> coveredVal(1,3.0);

    const char* in_file = "sphere.inputs";
    // read in an input file or use default file
    // if (argc > 1)
    // {
    //   in_file = argv[1];
    // }

    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    RealVect center;
    Real radius;
    RealVect origin;
    RealVect dx;
    Box domain;
    ProblemDomain pDomain(domain);

    int eekflag = 0;



    LevelData<EBCellFAB> fine, med, coarse;

    LevelData<EBCellFAB> fineRHS, medRHS, coarseRHS;
    LevelData<EBCellFAB> fineResidual, mediumResidual, coarseResidual;


    Vector<LevelData<EBCellFAB>* > ebvector(3,NULL);
    Vector<LevelData<EBCellFAB>* > vresidual(3,NULL);
    Vector<LevelData<EBCellFAB>* > rhsvector(3,NULL);
    ebvector[0]=&coarse;
    ebvector[1]=&med;
    ebvector[2]=&fine;
    vresidual[0]=&coarseResidual;
    vresidual[1]=&mediumResidual;
    vresidual[2]=&fineResidual;
    rhsvector[0] = &coarseRHS;
    rhsvector[1] = &medRHS;
    rhsvector[2] = &fineRHS;




    readGeometryInfo(domain,
                     dx,
                     origin,
                     center,
                     radius);

    Box domainFine(domain), domainMedi, domainCoar;
    ProblemDomain pFine(domain);
    RealVect dxFine(dx), dxMedi, dxCoar;

    CH_assert(eekflag == 0);

    domainMedi = coarsen(domainFine, 2);
    domainCoar = coarsen(domainMedi, 2);
    dxMedi = 2.0*dxFine;
    dxCoar = 2.0*dxMedi;
    Vector<RealVect> xVec(3, IntVect::Unit);
    xVec[0]*= dxCoar;
    xVec[1]*= dxMedi;
    xVec[2]*= dxFine;

    Vector<DisjointBoxLayout> grids(3);
    ProblemDomain baseDomain(domainCoar);
    ProblemDomain pMed(domainMedi);
    Vector<ProblemDomain> pd(3);
    pd[0] = baseDomain;
    pd[1] = pMed;
    pd[2] = ProblemDomain(domainFine);


    RefCountedPtr<BaseBCValue>value(new DirichletBC());
    DirichletPoissonDomainBC*  domainBC = new DirichletPoissonDomainBC();
    domainBC->setFunction(value);
    RefCountedPtr<BaseDomainBC> bc(domainBC);



    //make data holders
    Vector<int> comps(2,1);

    int steps= 5;
    int step = 0;


    while (step < steps)
      {


        eekflag = makeGeometry(
                               domain,
                               dx,
                               origin,
                               center,
                               radius);


        //make grids
        //IntVectSet tags = mfIndexSpace->interfaceRegion(2);
        IntVectSet   tags(domainCoar);
        tags.grow(1);
        makeHierarchy(grids, baseDomain, tags);


        const CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();

        Vector<EBISLayout> layouts(3);
        EBISLayout& fineLayout = layouts[2];
        ebisPtr->fillEBISLayout(fineLayout, grids[2], domainFine, 2);
        EBCellFactory fineFactory(fineLayout);
        ebvector[2]->define(grids[2], 1, IntVect::Unit, fineFactory);
        rhsvector[2]->define(grids[2], 1, IntVect::Zero, fineFactory);

        EBISLayout& medLayout = layouts[1];
        ebisPtr->fillEBISLayout(medLayout, grids[1], domainMedi, 2);
        EBCellFactory medFactory(medLayout);
        ebvector[1]->define(grids[1], 1, IntVect::Unit, medFactory);
        rhsvector[1]->define(grids[1], 1, IntVect::Zero, medFactory);

        EBISLayout& coarseLayout = layouts[0];
        ebisPtr->fillEBISLayout(coarseLayout, grids[0], domainCoar, 2);
        EBCellFactory coarseFactory(coarseLayout);
        ebvector[0]->define(grids[0], 1, IntVect::Unit, coarseFactory);
        rhsvector[0]->define(grids[0], 1, IntVect::Zero, coarseFactory);



        for (int lev=0; lev<3; lev++)
          {
            setValue(*rhsvector[lev], RHS(), pd[lev].domainBox(), xVec[lev], origin, true);
          }



        Vector<int> refRatio(3,2);

        int max_iter = 40;
        pp.get("max_iter", max_iter);
        Real eps = 1.e-6;
        pp.get("eps", eps);
        int relaxType;
        pp.get("relaxType",relaxType);

        DirichletPoissonDomainBCFactory* domDirBC = new DirichletPoissonDomainBCFactory();
        domDirBC->setFunction(value);
        RefCountedPtr<BaseDomainBCFactory> domBC( domDirBC );
        DirichletPoissonEBBCFactory* ebDirBC = new DirichletPoissonEBBCFactory();
        ebDirBC->setFunction(value);
        RefCountedPtr<BaseEBBCFactory> ebBC( ebDirBC );

        Vector<EBLevelGrid> eblgs(3);
        Vector<RefCountedPtr<EBQuadCFInterp> > quadCFI(3, RefCountedPtr<EBQuadCFInterp>());
        for (int i=0; i<3; i++)
          {
            eblgs[i] = EBLevelGrid(grids[i], layouts[i], pd[i]);
            if (i > 0)
              {
                quadCFI[i] = RefCountedPtr<EBQuadCFInterp>(
                             new EBQuadCFInterp(grids[i],
                                                grids[i-1],
                                                layouts[i],
                                                layouts[i-1],
                                                pd[i-1],
                                                refRatio[i-1],
                                                1, *eblgs[i].getCFIVS()));
              }
          }

        EBAMRPoissonOpFactory opFact(eblgs, refRatio, quadCFI, xVec[0], RealVect::Zero,
                                     4, relaxType, domBC, ebBC, 0.0, 1.0, 0.0,
                                     IntVect::Unit, IntVect::Zero);
        for (int i=0; i<3; i++)
          {

            LevelData<EBCellFAB> &phi=*ebvector[i], &rhs=*rhsvector[i], &residual=*vresidual[i];
            LevelData<EBCellFAB> correction;

            DisjointBoxLayout dblMGCoar;
            EBISLayout ebislMGCoar;

            EBAMRPoissonOp* opPtr = opFact.AMRnewOp(pd[i]);
            EBAMRPoissonOp& op = *opPtr;


            RelaxSolver<LevelData<EBCellFAB> > solver;
            solver.define(&op, false);
            solver.m_imax = max_iter;
            solver.m_eps  = eps;

            op.create(residual, rhs);
            op.create(correction, phi);
            op.setToZero(residual);
            op.setToZero(phi);

            op.residual(residual, phi, rhs);
            Real r2norm = op.norm(residual, 2);
            Real r0norm = op.norm(residual, 0);


            pout()<<indent<<"Residual L2 norm "<<r2norm<<"Residual max norm = "
                  <<r0norm<<std::endl;

            solver.solve(phi, rhs);

            op.residual(residual, phi, rhs);
            r2norm = op.norm(residual, 2);
            r0norm = op.norm(residual, 0);


            pout()<<indent2<<"Residual L2 norm "<<r2norm<<" Residual max norm = "
                  <<r0norm<<std::endl;

            delete opPtr;
          }


#ifdef CH_USE_HDF5
        sprintf(iter_str, "residual.%03d.%dd.hdf5",step, SpaceDim);
        Vector<std::string> names(1); names[0]="residual";

        writeEBHDF5(iter_str, grids, vresidual, names, domainCoar,
                    dxCoar[0], 1, step, refRatio, 3, true, coveredVal);

        sprintf(iter_str, "phi.%03d.%dd.hdf5",step, SpaceDim);
        names[0]="phi";
        writeEBHDF5(iter_str, grids, ebvector ,names, domainCoar,
                    dxCoar[0], 1, step, refRatio, 3, true, coveredVal);
#endif
        step++;

        center[0]-= dx[0]/3.0;
        center[1]-= dx[1]/2.0;
        radius += dx[0]/6.0;

        Chombo_EBIS::instance()->clear();
        pout()<<step<<std::endl;
      }

    pout() <<"\n "<<indent2<<pgmname<<" test passed " << endl;



  } // end scoping trick





#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}


void swapvector(Vector<LevelData<EBCellFAB>* >& left,
                Vector<LevelData<EBCellFAB>* >& right)

{

  LevelData<EBCellFAB>* swap;
  for (int i=0; i<left.size(); ++i)
    {
      swap = left[i];
      left[i]=right[i];
      right[i]=swap;
    }
}

int readGeometryInfo(Box& a_domain,
                     RealVect& a_dx,
                     RealVect& a_origin,
                     RealVect& a_center,
                     Real& a_radius)
{

  // parse input file
  ParmParse pp;

  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

  CH_assert(n_cell.size() == SpaceDim);

  IntVect lo = IntVect::Zero;
  IntVect hi;

  for (int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if (n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          return(-1);
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_domain.setSmall(lo);
  a_domain.setBig(hi);

  Vector<Real> prob_lo(SpaceDim, 1.0);
  Vector<Real> prob_hi(SpaceDim, 1.0);

  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.getarr("prob_hi",prob_hi,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_dx[idir]     = (prob_hi[idir]-prob_lo[idir])/n_cell[0];
      a_origin[idir] =  prob_lo[idir];
    }

  pp.get("radius",a_radius);

  // ParmParse doesn't get RealVects, so work-around with Vector<Real>
  Vector<Real> vectorCenter;
  pp.getarr("center",vectorCenter,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_center[idir] = vectorCenter[idir];
    }

  return 0;
} // end read from file

int makeGeometry(
                 const Box& a_domain,
                 const RealVect& a_dx,
                 const RealVect& a_origin,
                 const RealVect& a_center,
                 const Real& a_radius)
{

  int eekflag = 0;

  RealVect normal = BASISV(0);

  bool     insideRegular = false;
  SphereIF implicit(a_radius,a_center,insideRegular);

  int verbosity = 0;
  GeometryShop workshop(implicit,verbosity,a_dx);

  CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx[0], workshop);

  return eekflag;
}

void setValue(LevelData<EBCellFAB>& phase, const BaseBCValue& value,
              const Box& domain,
              const RealVect& dx, const RealVect& origin,
              bool useKappa)
{
  RealVect loc;
  IntVect  originIV = domain.smallEnd();
  DataIterator dit = phase.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      EBCellFAB& efab = phase[dit];

      FArrayBox& fab  = efab.getFArrayBox();
      ForAllX(Real, fab)
        {
          D_EXPR(loc[0]=iR+nR-nR+originIV[0], loc[1]=jR+originIV[1], loc[2]=kR+originIV[2]);
          loc+=0.5;
          loc*=dx;
          loc+=origin;
          fabR = value.value(loc, RealVect::Zero, 1.e99,0);
        } EndFor ;
      if (useKappa)
      {
        const EBISBox& ebox = efab.getEBISBox();
        IntVectSet ivs = ebox.getIrregIVS(efab.getRegion());
        IVSIterator it(ivs);
        for (;it.ok(); ++it)
        {
          const IntVect iv = it();
          VolIndex vi(iv,0);
          if (ebox.bndryArea(vi) != 0)
            {
              loc = iv;
              loc+=0.5;
              loc+= ebox.centroid(vi);
              loc*=dx;
              loc+=origin;
              efab(vi, 0) = value.value(loc, RealVect::Zero, 1.e99,0)*ebox.volFrac(vi);
            }
        }
      }
    }
}

void makeHierarchy(Vector<DisjointBoxLayout>& dbl,
                   const ProblemDomain& baseDomain,
                   const IntVectSet&baseTags)
{

  dbl.resize(3);
  Vector<int> refRatios(2,2);
  Real fillRatio = 0.85;
  int blockingFactor= 4;
  int bufferSize = 2;
  int maxSize = 32;
  BRMeshRefine regridder(baseDomain, refRatios, fillRatio, blockingFactor,
                         bufferSize, maxSize);

  Vector<Vector<Box> > oldGrids(3,1), newGrids(3);
  oldGrids[0][0]=baseDomain.domainBox();
  oldGrids[1][0]=coarsen(oldGrids[0][0], 2);
  oldGrids[2][0]=coarsen(oldGrids[1][0], 2);

  regridder.regrid(newGrids, baseTags, 0, 1, oldGrids);

  Vector<int> procs;
  for (int i=0; i<3; i++)
  {
    LoadBalance(procs, newGrids[i]);
    dbl[i] = DisjointBoxLayout(newGrids[i], procs);
  }

}
