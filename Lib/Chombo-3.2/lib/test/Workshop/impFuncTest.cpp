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

#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "DebugDump.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "PolyGeom.H"
#include "GeometryShop.H"

#include "EBFABView.H"

#include "SphereIF.H"
#include "CH_Attach.H"

// #include "TorusIF.H"
// #include "TiltedCylinderIF.H"
// #include "DataFileIF.H"
// #include "UnionIF.H"
// #include "IntersectionIF.H"
// #include "ComplementIF.H"
#include "UsingNamespace.H"

/***************/
// define a sphere EBIS.
/***************/
int makeGeometry(Box&      a_domain,
                 Real&     a_dx,
                 RealVect& a_origin,
                 RealVect& a_center,
                 Real&     a_radius,
                 bool&     a_insideRegular);

/******************/
/// define grids by splitting up domain
/******************/
int
makeLayout(DisjointBoxLayout& a_dbl1,
           const Box&         a_domain);

/***************/
// make the corresponding layout
/***************/
int makeEBISL(EBISLayout&              a_ebisl,
              const DisjointBoxLayout& a_grids,
              const Box&               a_domain,
              const int&               a_nghost);

/***************/
// check the corresponding layout
/***************/
int checkEBISL(const EBISLayout&        a_ebisl,
               const DisjointBoxLayout& a_grids,
               const Box&               a_domain,
               const Real&              a_dx,
               const RealVect&          a_origin,
               const Real&              a_radius,
               const RealVect&          a_center,
               const bool&              a_insideRegular);

int main (int argc, char** argv)
{
#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // begin forever present scoping trick
  {
    const char* in_file = "impFuncTest.inputs";

    // parse input file
    ParmParse pp(0,NULL,NULL,in_file);

    RealVect center;
    Real radius;

    bool insideRegular;

    RealVect origin;
    Real dx;

    Box domain;

    int eekflag = 0;

    // make geometry
    eekflag = makeGeometry(domain,
                           dx,
                           origin,
                           center,
                           radius,
                           insideRegular);

    if (eekflag != 0)
    {
      pout() << "non zero eek detected = " << eekflag << endl;
      MayDay::Error("problem in makeGeometry");
    }

    // make grids
    DisjointBoxLayout grids;
    eekflag = makeLayout(grids,domain);

    if (eekflag != 0)
    {
      pout() << "non zero eek detected = " << eekflag << endl;
      MayDay::Error("problem in makeLayouts");
    }

    // create ebislayout
    int nghost = 2;
    EBISLayout ebisl;
    eekflag = makeEBISL(ebisl,grids,domain,nghost);

    // make a LevelData
    int nComps = 2;

    IntVect ghost = IntVect::Unit;
    ghost *= nghost;

    RefCountedPtr<DataFactory<EBCellFAB> > rcpFactory(new EBCellFactory(ebisl));
    LevelData<EBCellFAB> level(grids,nComps,ghost,*rcpFactory);

    // writeEBLevelname(&level,"level.hdf5");

    if (eekflag != 0)
    {
      pout() << "non zero eek detected = " << eekflag << endl;
      MayDay::Error("problem in makeEBISL");
    }

    // check volume and surface area of approximate sphere
    eekflag = checkEBISL(ebisl,
                         grids,
                         domain,
                         dx,
                         origin,
                         radius,
                         center,
                         insideRegular);

    if (eekflag != 0)
    {
      pout() << "non zero eek detected = " << eekflag << endl;
      MayDay::Error("problem in checkEBISL");
    }

    pout() << "implicit function test passed" << endl;
  } // end scoping trick

  CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}

int makeGeometry(Box&      a_domain,
                 Real&     a_dx,
                 RealVect& a_origin,
                 RealVect& a_center,
                 Real&     a_radius,
                 bool&     a_insideRegular)
{
  int eekflag =  0;

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

  Vector<Real> prob_lo(SpaceDim,1.0);
  Real prob_hi;

  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);

  a_dx = (prob_hi-prob_lo[0])/n_cell[0];

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    a_origin[idir] = prob_lo[idir];
  }

  pp.get("radius",a_radius);

  // ParmParse doesn't get RealVects, so work-around with Vector<Real>
  Vector<Real> vectorCenter;
  pp.getarr("center",vectorCenter,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    a_center[idir] = vectorCenter[idir];
  }

  // Parm Parse doesn't get bools, so work-around with int
  int intInsideRegular;
  pp.get("intInsideRegular",intInsideRegular);

  if (intInsideRegular ==  1) a_insideRegular = true;
  if (intInsideRegular == -1) a_insideRegular = false;

  SphereIF implicit(a_radius,a_center,a_insideRegular);

  // TorusIF implicit(a_radius,0.5*a_radius,a_center,a_insideRegular);
  //
  // TiltedCylinderIF implicit(0.5*a_radius,a_center,a_origin,a_insideRegular);
  //
  // RealVect spacing(D_DECL(1,2,3));
  // RealVect origin(D_DECL(1,2,3));
  // IntVect num(D_DECL(1,2,3));
  // DataFileIF implicit(DataFileIF::ASCII,0.0,true);
  //
  // RealVect c1;
  // c1 = a_center;
  // c1[0] -= a_radius/3.0;
  // SphereIF sphere1(a_radius/1.5,c1,a_insideRegular);
  // ComplementIF invSph1(sphere1);
  //
  // RealVect c2;
  // c2 = a_center;
  // c2[0] += a_radius/3.0;
  // SphereIF sphere2(a_radius/1.5,c2,a_insideRegular);
  // ComplementIF invSph2(sphere2);
  //
  // SphereIF sphere3(a_radius/2.5,a_center,a_insideRegular);
  //
  // IntersectionIF inter(sphere1,invSph2);
  // UnionIF implicit(inter,sphere3);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(implicit,0,vectDx);

  // this generates the new EBIS
  CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx,workshop);

  return eekflag;
}

int makeLayout(DisjointBoxLayout& a_dbl,
               const Box&         a_domain)
{
  ParmParse pp;
  int eekflag = 0;

  int ipieces;
  ipieces = Max(ipieces,1);

  int maxsize;
  Vector<Box> vbox(1,a_domain);

  pp.get("maxboxsize",maxsize);

  domainSplit(a_domain,vbox,maxsize);
  if (eekflag != 0)
  {
    pout() << "problem in domainsplit" << endl;
    return eekflag;
  }

  Vector<int> procAssign;
  eekflag = LoadBalance(procAssign,vbox);

  if (eekflag != 0)
  {
    pout() << "problem in loadbalance" << endl;
    return eekflag;
  }

  a_dbl.define(vbox,procAssign);

  return eekflag;
}

int makeEBISL(EBISLayout&              a_ebisl,
              const DisjointBoxLayout& a_grids,
              const Box&               a_domain,
              const int&               a_nghost)
{
  const CH_XD::EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
 CH_assert(ebisPtr->isDefined());

  ebisPtr->fillEBISLayout(a_ebisl,a_grids,a_domain,a_nghost);
  return 0;
}

int checkEBISL(const EBISLayout&        a_ebisl,
               const DisjointBoxLayout& a_grids,
               const Box&               a_domain,
               const Real&              a_dx,
               const RealVect&          a_origin,
               const Real&              a_radius,
               const RealVect&          a_center,
               const bool&              a_insideRegular)
{
  //check to see that the sum of the volume fractions
  //comes to close to exactly the total volume
  int eekflag =  0;
  //First: calculate volume of domain

  const IntVect& ivSize = a_domain.size();
  RealVect hiCorn;
  RealVect domLen;
  Real cellVolume = 1.0;
  Real totalVolume = 1.0;

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    hiCorn[idir] = a_origin[idir] + a_dx*Real(ivSize[idir]);
    cellVolume *= a_dx;
    domLen[idir] = hiCorn[idir] - a_origin[idir];
    totalVolume *= domLen[idir];
  }

  // Calculate the exact volume of the regular region
  Real volExact;
  Real volError;
  Real bndryExact;
  Real bndryError;

  if (SpaceDim == 2)
  {
    volExact = acos(-1.0) * a_radius*a_radius;
    volError = 0.0001851;

    bndryExact = 2*acos(-1.0) * a_radius;
    bndryError = 0.000047;
  }

  if (SpaceDim == 3)
  {
    volExact = 4.0/3.0 * acos(-1.0) * a_radius*a_radius*a_radius;
    volError = 0.000585;

    bndryExact = 4.0 * acos(-1.0) * a_radius*a_radius;
    bndryError = 0.000287;
  }

  if (a_insideRegular == false)
  {
    volExact = totalVolume - volExact;
  }

  // Now calculate the volume of approximate sphere
  Real tolerance = 0.001;
  Real volTotal = 0.0;
  Real bndryTotal = 0.0;


  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
  {
    const EBISBox& ebisBox = a_ebisl[dit()];
    for (BoxIterator bit(a_grids.get(dit())); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      Vector<VolIndex> vofs = ebisBox.getVoFs(iv);

      if (vofs.size() > 1)
      {
        eekflag = 2;
        return eekflag;
      }

      for (int ivof = 0; ivof < vofs.size(); ivof++)
      {
        Real volFrac = ebisBox.volFrac(vofs[ivof]);
        if (ebisBox.isIrregular(iv))
          {
            pout() << "iv = " << iv << ", volFrac = " << volFrac << endl;
          }
        volTotal += volFrac * cellVolume;

        Real bndryArea = ebisBox.bndryArea(vofs[ivof]);
        bndryTotal += bndryArea * cellVolume / a_dx;

        if (volFrac > tolerance && volFrac < 1.0 - tolerance)
        {
          if (!ebisBox.isIrregular(iv))
          {
            eekflag = 4;
            return eekflag;
          }
        }
      }
    }
  }
#ifdef CH_MPI
  Real t=volTotal;
  MPI_Allreduce ( &t, &volTotal, 1,
                  MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  t=bndryTotal;
  MPI_Allreduce ( &t, &bndryTotal, 1,
                   MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
#endif

  //check how close is the answer
  pout() << endl;

  Real error;
  error = volTotal - volExact;

  if (Abs(error/volExact) > volError)
  {
    eekflag = 5;
  }

  pout() << "volTotal  = " << volTotal             << endl;
  pout() << "volExact  = " << volExact             << endl;
  pout() << "error     = " << error                << endl;
  pout() << "relError  = " << Abs(error/volExact) << endl;
  pout() << "tolerance = " << volError             << endl;
  pout() << endl;

  error = bndryTotal - bndryExact;

  if (Abs(error/bndryExact) > bndryError)
  {
    eekflag = 6;
  }

  pout() << "bndryTotal = " << bndryTotal             << endl;
  pout() << "bndryExact = " << bndryExact             << endl;
  pout() << "error      = " << error                  << endl;
  pout() << "relError   = " << Abs(error/bndryExact) << endl;
  pout() << "tolerance  = " << bndryError             << endl;
  pout() << endl;

  return eekflag;
}
