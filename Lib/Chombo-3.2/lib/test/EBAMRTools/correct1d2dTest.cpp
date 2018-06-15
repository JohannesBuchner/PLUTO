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
#include "RealVect.H"
#include "EBFluxFAB.H"
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
#include "EBArith.H"
#include "PolyGeom.H"
#include "Correct1D2D.H"
#include "EBDebugDump.H"
#include "EBAMRIO.H"
#include "AMRIO.H"
#include "DebugDump.H"
#include "EBCellFactory.H"
#include "EBLevelDataOps.H"
#include "EBFABView.H"
#include "UsingNamespace.H"


int makeGeometry(Box& a_domain)
{
  Real dx;
  RealVect origin;
  int eekflag =  0;
  //parse input file
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
  Real prob_hi;
  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);
  dx = (prob_hi-prob_lo[0])/n_cell[0];
  RealVect dxVect = dx*RealVect::Unit;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      origin[idir] = prob_lo[idir];
    }
  int verbosity = 0;
  int whichgeom;
  pp.get("which_geom",whichgeom);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();

  if (whichgeom == 0)
    {
      //allregular
      pout() << "all regular geometry" << endl;
      AllRegularService regserv;
      ebisPtr->define(a_domain, origin, dx, regserv);
    }
  else if (whichgeom == 1)
    {
      pout() << "ramp geometry" << endl;
      int upDir;
      int indepVar;
      Real startPt;
      Real slope;
      pp.get("up_dir",upDir);
      pp.get("indep_var",indepVar);
      pp.get("start_pt", startPt);
      pp.get("ramp_slope", slope);

      RealVect normal = RealVect::Zero;
      normal[upDir] = 1.0;
      normal[indepVar] = -slope;

      RealVect point = RealVect::Zero;
      point[upDir] = -slope*startPt;

      bool normalInside = true;

      PlaneIF ramp(normal,point,normalInside);

      GeometryShop workshop(ramp,verbosity,dxVect);
      //this generates the new EBIS
      ebisPtr->define(a_domain, origin, dx, workshop);
    }
  else if (whichgeom == 5)
    {
      pout() << "sphere geometry" << endl;
      vector<Real> sphere_center(SpaceDim);
      pp.getarr("sphere_center",sphere_center, 0, SpaceDim);
      RealVect sphereCenter;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          sphereCenter[idir] = sphere_center[idir];
        }
      Real sphereRadius;
      pp.get("sphere_radius", sphereRadius);

      bool     insideRegular = false;
      SphereIF implicit(sphereRadius,sphereCenter,insideRegular);
      GeometryShop workshop(implicit,verbosity,dxVect);
      //this generates the new EBIS
      ebisPtr->define(a_domain, origin, dx, workshop);
    }
  else if (whichgeom == 13)
    {
      pout() << "rhodonea geometry" << endl;
      vector<Real> tmp(SpaceDim);
      pp.getarr("rhodonea_center", tmp, 0, SpaceDim);
      RealVect rhodoneaCenter;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          rhodoneaCenter[idir] = tmp[idir];
        }
      Real innerRadius;
      pp.get("inner_radius", innerRadius);

      Real outerRadius;
      pp.get("outer_radius", outerRadius);

      int frequency;
      pp.get("frequency", frequency);

      bool     insideRegular = false;

      RhodoneaIF implicit(innerRadius, outerRadius, frequency,
                          rhodoneaCenter, insideRegular);

      GeometryShop workshop(implicit,verbosity,dxVect);
      //this generates the new EBIS
      ebisPtr->define(a_domain, origin, dx, workshop);
    }
  else
    {
      //bogus which_geom
      pout() << " bogus which_geom input = " << whichgeom;
      eekflag = 33;
    }

  return eekflag;
}

Real
divergence(const EBFluxFAB&  a_func,
           const RealVect&   a_bndryFlux,
           const EBISBox&    a_ebisBox,
           const VolIndex&   a_vof,
           const Real&       a_dx)
{
  Real retval = 0;

  Real bndryArea = a_ebisBox.bndryArea(a_vof);
  RealVect normal = a_ebisBox.normal(a_vof);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real bndryFlux    = a_bndryFlux[idir]; //normal already dealt with
      Real bndryContrib = -bndryFlux*bndryArea*normal[idir];
      Real openContrib = 0;
      if (bndryArea > 1.0e-3)
        {
          openContrib = 0;
        }
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int isign = sign(sit());
          Real rsign = isign;
          Vector<FaceIndex> faces = a_ebisBox.getFaces(a_vof, idir, sit());
          for (int iface = 0; iface < faces.size(); ++iface)
            {
              Real areaFrac = a_ebisBox.areaFrac(faces[iface]);
              Real faceFlux = a_func[idir](faces[iface], 0);
              openContrib += rsign*faceFlux*areaFrac;
            }
        }
      retval += openContrib + bndryContrib;

    }
  retval /= a_dx;
  return retval;
}
/****/
void divergence(EBCellFAB&                a_divF,
                const EBFluxFAB&          a_flux,
                const EBISBox&            a_ebisBox,
                const Box&                a_box,
                const RealVect&           a_fluxVal,
                const Real&               a_dx)
{
  IntVectSet ivs(a_box);
  for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      Real divval = divergence(a_flux, a_fluxVal, a_ebisBox, vofit(), a_dx);
      a_divF(vofit(), 0) = divval;
    }

}

/****/
void setBorkedFlux(EBFluxFAB&          a_flux,
                   const EBISBox&      a_ebisBox,
                   const Box&          a_box,
                   const RealVect&     a_fluxVal,
                   const BaseFab<int>& a_map)
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real bogval = 0;
      //the bogus value will be zero so it will not
      //do to have the correct value be zero
      if (Abs(a_fluxVal[idir]) < 1.0e-3) bogval = 1.0;

      a_flux[idir].setVal(a_fluxVal[idir]);
      //if i am in a 1d box and the box next to me
      //is 2d, set the border flux to
      for (SideIterator sit; sit.ok(); ++sit)
        {
          Box borderBox;
          if (sit() == Side::Lo)
            {
              borderBox = adjCellLo(a_box, idir,  -1);
            }
          else
            {
              borderBox = adjCellHi(a_box, idir,  -1);
            }
          for (BoxIterator bit(borderBox); bit.ok(); ++bit)
            {
              const IntVect& thisIV = bit();
              IntVect thatIV = thisIV + sign(sit())*BASISV(idir);
              if ((a_map(thisIV, 0) == 1) && (a_map(thatIV, 0) == 2))
                {
                  Vector<FaceIndex> faces = a_ebisBox.getAllFaces(thisIV, idir, sit());
                  for (int iface = 0; iface < faces.size(); iface++)
                    {
                      a_flux[idir](faces[iface],0) = bogval;
                    }
                }
            }
        }
    }

}
/****/
int checkCorrection(const Box& a_domain)
{
  int eekflag = 0;
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  CH_assert(ebisPtr->isDefined());
  //create a dbl of split up into 2^D boxes
  Vector<Box> vbox;
  Vector<int> proc;
  domainSplit(a_domain, vbox, a_domain.size(0)/2);
  LoadBalance(proc, vbox);
  DisjointBoxLayout dbl(vbox, proc);
  LayoutData<bool>  map1d(dbl);
  int ibox = 0;

  //make every other map 1d
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      map1d[dit()] = (ibox%2 == 0);
      ibox++;
    }

  LevelData< BaseFab<int> > intmap(dbl, 1 , IntVect::Unit);
  Correct1D2D::makeIntMap(intmap, map1d, dbl);

  //this defines an ebisl with 2 ghost cells
  EBLevelGrid eblg(dbl, ProblemDomain(a_domain), 2, Chombo_EBIS::instance()) ;

  RealVect fluxVal = RealVect::Unit;
  Real dx = 1.0/a_domain.size(0);

  EBCellFactory fact(eblg.getEBISL());
  LevelData<EBCellFAB> divU(dbl, 1, IntVect::Zero, fact);
  Correct1D2D correctinator(eblg, map1d, 1);
  correctinator.setToZero(); //this is important
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      //make a flux that is borked on the 1d side of box boundaries
      EBFluxFAB fluxFunc(     eblg.getEBISL()[dit()], dbl[dit()], 1);
      setBorkedFlux(fluxFunc, eblg.getEBISL()[dit()], dbl[dit()], fluxVal, intmap[dit()]);
      //increment 1d buffers in corrector
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          correctinator.increment1D(fluxFunc[idir], 1./dx, dit());
          correctinator.increment2D(fluxFunc[idir], 1./dx, dit());
        }

      //now take the divergence of the flux.   it should be zero except
      //at box boundaries on the 1d side
      divergence(divU[dit()], fluxFunc, eblg.getEBISL()[dit()], dbl[dit()], fluxVal, dx);
    }
  //find max value of divergence before corrections

  Real max, min;
  EBLevelDataOps::getMaxMin(max, min, divU, 0);
  pout() << "before correction max divU = " << max << ", min divU = " << min << endl;

  //correct solution due to mismatch at 1d-2d boundaries
  correctinator.correctSolution(divU);

  //recalc max min
  EBLevelDataOps::getMaxMin(max, min, divU, 0);
  pout() << "after  correction max divU = " << max << ", min divU = " << min << endl;
  if (Abs(max-min) > 1.0e-3)
    {
      pout() << "corrector did not seem to correct" << endl;
      return -7;
    }
  return eekflag;
}
int
main(int argc,char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // begin forever present scoping trick
  {
    const char* in_file;
    if (argc < 2 || argv[1][0] == '-')
      {
        in_file = "correct12.inputs";
      }
    else
      {
        in_file = argv[1];
      }
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    Box domain;
    int eekflag = 0;
    eekflag =  makeGeometry(domain);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeGeometry");
      }

    eekflag = checkCorrection(domain);
    if (eekflag != 0)
      {
        pout() << "checkCorrection: eek = " << eekflag << endl;
        MayDay::Error("problem in checkCorrection");
      }
    pout() << "correct1d2d test passed" << endl;
  } // end scoping trick

  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}



