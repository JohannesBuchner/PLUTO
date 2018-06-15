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

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "VoFIterator.H"
#include "AllRegularService.H"
#include "GeometryShop.H"
#include "PlaneIF.H"
#include "SphereIF.H"
#include "RhodoneaIF.H"
#include "UnionIF.H"
#include "ComplementIF.H"
#include "TransformIF.H"
#include "TiltedCylinderIF.H"
#include "IntersectionIF.H"
#include "DirichletViscousTensorEBBC.H"
#include "BaseBCFuncEval.H"
#include "EBCellFactory.H"
#include "EBArith.H"
#include "PolyGeom.H"
#include "EBDebugDump.H"
#include "DebugDump.H"
#include "EBAMRIO.H"
#include "EBFABView.H"
#include "UsingNamespace.H"
//each vel comp  = constant
int s_nghost = 4;
class ConstBCFunc : public BaseBCFuncEval
{
public:
  //just set the darn thing
  RealVect m_value;

  virtual Real value(const RealVect&       a_point,
                     const int&            a_comp) const
  {
    return m_value[a_comp];
  }

  virtual Real derivative(const RealVect&       a_point,
                          const int&            a_comp,
                          const int&            a_derivDir
                          ) const
  {
    return  0.0;
  }
};

//each vel comp  = constant + grad*dist
class LinearBCFunc : public BaseBCFuncEval
{
public:
  //just set the darn thing
  RealVect m_cons;
  RealVect m_origin;
  LinearBCFunc()
  {
  }
  //indexed by grad[component][direction]
  Real m_grad[SpaceDim][SpaceDim];

  virtual Real value(const RealVect&       a_point,
                     const int&            a_comp) const
  {
    Real retval = m_cons[a_comp];
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        retval += m_grad[a_comp][idir]*(a_point[idir]-m_origin[idir]);
      }
    return retval;
  }

  virtual Real derivative(const RealVect&       a_point,
                          const int&            a_comp,
                          const int&            a_derivDir
                          ) const
  {
    return  m_grad[a_comp][a_derivDir];
  }
};

//each vel comp  = constant + grad*dist
class DistBCFunc : public BaseBCFuncEval
{
public:
  //just set the darn thing
  RealVect m_cons;
  RealVect m_lineNorm;
  RealVect m_lineOrigin;
  DistBCFunc()
  {
    m_lineNorm   = RealVect::Unit;
    m_lineOrigin = RealVect::Zero;
  }

  virtual Real value(const RealVect&       a_point,
                     const int&            a_comp) const
  {
    Real retval = m_cons[a_comp];
    Real distanceToPlane =
      PolyGeom::distancePointToPlane(a_point, m_lineNorm, m_lineOrigin);

    retval += distanceToPlane;
    return retval;
  }

  virtual Real derivative(const RealVect&       a_point,
                          const int&            a_comp,
                          const int&            a_derivDir
                          ) const
  {
    return  0;
  }
};
bool nearDomain(const VolIndex& a_vof, const Box& a_domain)
{
  bool retval = false;
  const IntVect& loside = a_domain.smallEnd();
  const IntVect& hiside = a_domain.bigEnd();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if ((a_vof.gridIndex()[idir] == loside[idir]) || (a_vof.gridIndex()[idir] == hiside[idir]))
        {
          retval = true;
        }
    }
  return retval;
}
/****/
int makeLayout(DisjointBoxLayout& a_dbl,
               EBISLayout&        a_ebisl,
               const Box& a_domain)
{
  a_dbl = DisjointBoxLayout(Vector<Box>(1, a_domain), Vector<int>(1, 0));
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->fillEBISLayout(a_ebisl, a_dbl, a_domain, s_nghost);
  return 0;
}
/****/
int makeGeometry(Box& a_domain, RealVect& a_point, RealVect& a_normal)
{
  pout() << "ramp geometry" << endl;
  Real ncell = 16;
  a_domain = Box(IntVect::Zero, (ncell-1)*IntVect::Unit);
  RealVect normal = RealVect::Zero;
  normal[0] = -1;
  normal[1] = 2;
  RealVect point = RealVect::Zero;
  point[0]=0.123456789;
  bool normalInside = true;
  a_normal = normal;
  a_point =  point;
  PlaneIF ramp(normal,point,normalInside);

  Real dx = 1.0/ncell;
  GeometryShop workshop(ramp, 0,dx*RealVect::Unit);
  //this generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain, RealVect::Zero, dx, workshop);
  return 0;
}
/***/
void
setPhi(LevelData<EBCellFAB>&    a_phi,
       const BaseBCFuncEval&    a_func,
       const DisjointBoxLayout& a_dbl,
       const EBISLayout&        a_ebisl,
       const Box&               a_domain,
       const Real&              a_dx,
       const RealVect&          a_origin)
{
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      IntVectSet ivs(a_dbl[dit()]);
      for (VoFIterator vofit(ivs, a_ebisl[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          RealVect point = EBArith::getVofLocation(vofit(), a_dx*RealVect::Unit, a_origin);
          for (int icomp = 0; icomp < SpaceDim; icomp++)
            {
              Real output = a_func.value(point, icomp);
              a_phi[dit()](vofit(), icomp) = output;
            }
        }
    }
  a_phi.exchange();
}
/***/
int
constantFunctionCheck(Box& a_domain, RealVect& a_origin, RealVect& a_normal)
{
  DisjointBoxLayout dbl;
  EBISLayout ebisl;
  makeLayout(dbl, ebisl, a_domain);
  Real dx = 1.0/a_domain.size(0);
  IntVect ghostiv = s_nghost*IntVect::Unit;
  DirichletViscousTensorEBBC bc(ProblemDomain(a_domain), ebisl, dx*RealVect::Unit, &ghostiv, &ghostiv);
  EBCellFactory fact(ebisl);
  LevelData<EBCellFAB> phi(dbl, SpaceDim, ghostiv, fact);
  RealVect consRV = RealVect::Unit;
  for (int idir = 0; idir<SpaceDim; idir++)
    consRV[idir] *= (7.0+idir);
  ConstBCFunc* constfunc = new ConstBCFunc();
  RefCountedPtr<BaseBCFuncEval> func = RefCountedPtr<BaseBCFuncEval>(constfunc);
  constfunc->m_value = consRV;

  Real rightAns = 0;
  Real tol = 1.0e-8;
#ifdef CH_USE_FLOAT
    tol = 1.0e-2;
#endif
  setPhi(phi, *func, dbl, ebisl, a_domain, dx, a_origin);

  bc.setFunction(func);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      IntVectSet ivs = ebisl[dit()].getIrregIVS(dbl[dit()]);
      for (VoFIterator vofit(ivs, ebisl[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          //only check real irregular cells
          Real bndryArea = ebisl[dit()].bndryArea(vofit());
          if (bndryArea > 0.1)
            {
              Real grad[SpaceDim][SpaceDim];
              bc.getGradient(grad, vofit(), phi[dit()], ebisl[dit()], dit(), dx, false);
              for (int idir = 0; idir < SpaceDim; idir ++)
                {
                  for (int icomp = 0; icomp < SpaceDim; icomp++)
                    {
                      if (Abs(grad[icomp][idir] - rightAns) > tol)
                        {
                          pout() << "failed constant function test" << endl;
                          return -1;
                        }
                    }
                }
            }
        }
    }

  Real value = 7;
  consRV = value*RealVect::Unit;
  constfunc->m_value = consRV;
  setPhi(phi, *func, dbl, ebisl, a_domain, dx, a_origin);

  bc.setValue(value);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      IntVectSet ivs = ebisl[dit()].getIrregIVS(dbl[dit()]);
      for (VoFIterator vofit(ivs, ebisl[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          //only check real irregular cells
          Real bndryArea = ebisl[dit()].bndryArea(vofit());
          if (bndryArea > 0.1)
            {
              Real grad[SpaceDim][SpaceDim];
              bc.getGradient(grad, vofit(), phi[dit()], ebisl[dit()], dit(), dx, false);
              for (int idir = 0; idir < SpaceDim; idir ++)
                {
                  for (int icomp = 0; icomp < SpaceDim; icomp++)
                    {
                      if (Abs(grad[icomp][idir] - rightAns) > tol)
                        {
                          pout() << "failed constant value test" << endl;
                          return -2;
                        }
                    }
                }
            }
        }
    }

  return 0;
}

/***/
int
linearFunctionCheck(Box& a_domain, RealVect& a_origin, RealVect& a_normal)
{
  //first try the same function as the last one.
  {
    DistBCFunc* linearfunc = new DistBCFunc;
    RefCountedPtr<BaseBCFuncEval> func(linearfunc);
    RealVect consRV = RealVect::Unit;
    linearfunc->m_cons = consRV;
    linearfunc->m_lineNorm = a_normal;
    linearfunc->m_lineOrigin = a_origin;

    Real rightAns[SpaceDim][SpaceDim];
    RealVect unitNormal = a_normal;
    Real sum;
    PolyGeom::unifyVector(unitNormal, sum);
    for (int idir = 0; idir<SpaceDim; idir++)
      {
        for (int icomp = 0; icomp < SpaceDim; icomp++)
          {
            rightAns[icomp][idir] = unitNormal[idir];
          }
      }

    DisjointBoxLayout dbl;
    EBISLayout ebisl;
    makeLayout(dbl, ebisl, a_domain);
    Real dx = 1.0/a_domain.size(0);
    IntVect ghostiv = s_nghost*IntVect::Unit;
    DirichletViscousTensorEBBC bc(a_domain, ebisl, dx*RealVect::Unit, &ghostiv, &ghostiv);
    EBCellFactory fact(ebisl);
    LevelData<EBCellFAB> phi(dbl, SpaceDim, ghostiv, fact);

    Real tol = 1.0e-8;
#ifdef CH_USE_FLOAT
    tol = 1.0e-2;
#endif
    setPhi(phi, *func, dbl, ebisl, a_domain, dx, RealVect::Zero);

    //should be a constant along line
    bc.setFunction(func);
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
      {
        IntVectSet ivs = ebisl[dit()].getIrregIVS(dbl[dit()]);
        for (VoFIterator vofit(ivs, ebisl[dit()].getEBGraph()); vofit.ok(); ++vofit)
          {
            const VolIndex& vof = vofit();
            //only check real irregular cells
            Real bndryArea = ebisl[dit()].bndryArea(vofit());
            if (bndryArea > 0.1)
              {
                Real grad[SpaceDim][SpaceDim];
                bc.getGradient(grad, vof, phi[dit()], ebisl[dit()], dit(), dx, false);
                for (int idir = 0; idir < SpaceDim; idir ++)
                  {
                    for (int icomp = 0; icomp < SpaceDim; icomp++)
                      {
                        if (Abs(grad[icomp][idir] - rightAns[icomp][idir]) > tol)
                          {
                            pout() << "failed dist function calculated gradient test" << endl;
                            return -3;
                          }
                      }
                  }


              }
          }
      }

  }
  pout() << "passed distance function bit" << endl;
  {

    //
    LinearBCFunc* linearfunc = new LinearBCFunc;
    RefCountedPtr<BaseBCFuncEval> func(linearfunc);
    RealVect consRV = RealVect::Unit;
    for (int idir = 0; idir<SpaceDim; idir++)
      consRV[idir] *= (7.0+idir);
    //    linearfunc->m_cons = consRV;
    //    consRV = RealVect::Zero;
    //    linearfunc->m_cons = RealVect::Zero;

    Real rightAns[SpaceDim][SpaceDim];
    for (int idir = 0; idir<SpaceDim; idir++)
      {
        for (int icomp = 0; icomp < SpaceDim; icomp++)
          {
            rightAns[icomp][idir] = idir;
            linearfunc->m_grad[icomp][idir] = rightAns[icomp][idir];
          }
      }

    DisjointBoxLayout dbl;
    EBISLayout ebisl;
    makeLayout(dbl, ebisl, a_domain);
    Real dx = 1.0/a_domain.size(0);
    IntVect ghostiv = s_nghost*IntVect::Unit;
    DirichletViscousTensorEBBC bc(a_domain, ebisl, dx*RealVect::Unit, &ghostiv, &ghostiv);
    EBCellFactory fact(ebisl);
    LevelData<EBCellFAB> phi(dbl, SpaceDim, ghostiv, fact);

    Real tol = 1.0e-8;
#ifdef CH_USE_FLOAT
    tol = 1.0e-2;
#endif
    setPhi(phi, *func, dbl, ebisl, a_domain, dx, a_origin);

    bc.setFunction(func);
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
      {
        IntVectSet ivs = ebisl[dit()].getIrregIVS(dbl[dit()]);
        for (VoFIterator vofit(ivs, ebisl[dit()].getEBGraph()); vofit.ok(); ++vofit)
          {
            //only check real irregular cells
            Real bndryArea = ebisl[dit()].bndryArea(vofit());
            if (bndryArea > 0.1)
              {
                Real grad[SpaceDim][SpaceDim];
                bc.getGradient(grad, vofit(), phi[dit()], ebisl[dit()], dit(), dx, false);
                for (int idir = 0; idir < SpaceDim; idir ++)
                  {
                    for (int icomp = 0; icomp < SpaceDim; icomp++)
                      {
                        if (Abs(grad[icomp][idir] - rightAns[icomp][idir]) > tol)
                          {
                            pout() << "failed linear function calculated gradient test " << endl;
                            return -4;
                          }
                      }
                  }
              }
          }
      }
  }
  pout() << "passed linear function bit" << endl;
  return 0;
}

/***/
int
linearValueCheck(Box& a_domain, RealVect& a_origin, RealVect& a_normal)
{
  DistBCFunc* linearfunc = new DistBCFunc;
  RefCountedPtr<BaseBCFuncEval> func(linearfunc);
  RealVect consRV = RealVect::Unit;
  linearfunc->m_cons = consRV;
  linearfunc->m_lineNorm = a_normal;
  linearfunc->m_lineOrigin = a_origin;

  Real rightAns[SpaceDim][SpaceDim];
  RealVect unitNormal = a_normal;
  Real sum;
  PolyGeom::unifyVector(unitNormal, sum);
  for (int idir = 0; idir<SpaceDim; idir++)
    {
      for (int icomp = 0; icomp < SpaceDim; icomp++)
        {
          rightAns[icomp][idir] = unitNormal[idir];
        }
    }

  DisjointBoxLayout dbl;
  EBISLayout ebisl;
  makeLayout(dbl, ebisl, a_domain);
  Real dx = 1.0/a_domain.size(0);
  IntVect ghostiv = s_nghost*IntVect::Unit;
  DirichletViscousTensorEBBC bc(a_domain, ebisl, dx*RealVect::Unit, &ghostiv, &ghostiv);
  EBCellFactory fact(ebisl);
  LevelData<EBCellFAB> phi(dbl, SpaceDim, ghostiv, fact);

  Real tol = 1.0e-8;
#ifdef CH_USE_FLOAT
  tol = 1.0e-2;
#endif
  setPhi(phi, *func, dbl, ebisl, a_domain, dx, RealVect::Zero);

  //should be a constant along line
  bc.setValue(consRV[0]);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      IntVectSet ivs = ebisl[dit()].getIrregIVS(dbl[dit()]);
      for (VoFIterator vofit(ivs, ebisl[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          //only check real irregular cells
          Real bndryArea = ebisl[dit()].bndryArea(vofit());
          bool vofNearDomain = nearDomain(vof, a_domain);
          if (!vofNearDomain && (bndryArea > 0.1))
            {
              Real grad[SpaceDim][SpaceDim];
              bc.getGradient(grad, vof, phi[dit()], ebisl[dit()], dit(), dx, false);
              for (int idir = 0; idir < SpaceDim; idir ++)
                {
                  for (int icomp = 0; icomp < SpaceDim; icomp++)
                    {
                      if (Abs(grad[icomp][idir] - rightAns[icomp][idir]) > tol)
                        {
                          pout() << "failed linear value calulated gradient test" << endl;
                          return -7;
                        }
                    }
                }

                Real gradSten[SpaceDim][SpaceDim];
                bc.getGradientStenValue(gradSten, vofit(), phi[dit()], ebisl[dit()], dit(), dx, false);
                for (int idir = 0; idir < SpaceDim; idir ++)
                  {
                    for (int icomp = 0; icomp < SpaceDim; icomp++)
                      {
                        if (Abs(gradSten[icomp][idir] - rightAns[icomp][idir]) > tol)
                          {
                            pout() << "failed linear value stencil gradient test" << endl;
                            return -5;
                          }
                      }
                  }

            }
        }
    }

  return 0;
}
/***/
int
main(int argc,char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // begin forever present scoping trick
  {
    Box domain;
    int eekflag = 0;
    RealVect point, normal;
    eekflag =  makeGeometry(domain, point, normal);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeGeometry");
      }
    eekflag = constantFunctionCheck(domain, point, normal);

    if (eekflag == 0)
      {
        pout() << "testDirVTEBBC: constantFunctionCheck passed" << endl;
      }
    else
      {
        pout() << "constantFunctionCheck: eek = " << eekflag << endl;
        MayDay::Error("problem in constantFunctionCheck");
      }

    eekflag = linearValueCheck(domain, point, normal);
    if (eekflag == 0)
      {
        pout() << "testDirVTEBBC: linearValueCheck passed" << endl;
      }
    else
      {
        pout() << "linearValueCheck: eek = " << eekflag << endl;
        MayDay::Error("problem in linearValueCheck");
      }

    /**
    eekflag = linearFunctionCheck(domain, point, normal);
    if (eekflag == 0)
      {
        pout() << "testDirVTEBBC: linearFunctionCheck passed" << endl;
      }
    else
      {
        pout() << "linearFunctionCheck: eek = " << eekflag << endl;
        MayDay::Error("problem in linearFunctionCheck");
      }
    **/
    pout() << "DirichletViscousTensorEBBC gradient test passed" << endl;
  } // end scoping trick

  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}



