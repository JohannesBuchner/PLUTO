#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBGradDivFilter.H"
#include "EBGradDivFilterF_F.H"
#include "EBCellFactory.H"
#include "EBFluxFactory.H"
#include "FaceIterator.H"
#include "VoFIterator.H"
#include "PolyGeom.H"
#include "EBArith.H"
#include "EBLevelDataOps.H"
#include "EBLevelCCProjector.H"

#include "NamespaceHeader.H"

int
EBGradDivFilter::
getGradComp(int a_velDir, int a_derivDir)
{
  int gradcomp = SpaceDim*a_velDir + a_derivDir;
  return gradcomp;
}

EBGradDivFilter::
~EBGradDivFilter()
{
  if (m_tensorCFI != NULL)
    {
      delete m_tensorCFI;
    }
}
/***********************/
EBGradDivFilter::
EBGradDivFilter(const DisjointBoxLayout&                  a_gridsFine,
                const DisjointBoxLayout*                  a_gridsCoarPtr,
                const EBISLayout&                         a_ebislFine,
                const EBISLayout*                         a_ebislCoarPtr,
                const ProblemDomain&                      a_domainFine,
                const RealVect&                           a_dxFine,
                const int&                                a_refRat,
                const Real&                               a_lambdaScale,
                const EBIndexSpace* const                 a_ebisPtr)
{
  CH_TIME("EBGradDivFilter::EBGradDivFilter");
  m_gridsFine   =  a_gridsFine    ;
  m_gridsCoarPtr=  a_gridsCoarPtr ;
  m_ebislFine   =  a_ebislFine    ;
  m_eblgFine.define(m_gridsFine, m_ebislFine, m_domainFine);

  m_ebislCoarPtr=  a_ebislCoarPtr ;
  m_domainFine  =  a_domainFine   ;
  m_refRat      =  a_refRat       ;
  m_dxFine      =  a_dxFine;
  m_hasCoarser = (m_gridsCoarPtr != NULL);
  EBCellFactory ebcellfact(m_ebislFine);
  EBFluxFactory ebfluxfact(m_ebislFine);
  m_gradVel.define(m_gridsFine, SpaceDim*SpaceDim, IntVect::Unit, ebcellfact);
  m_lambda.define( m_gridsFine,                 1, IntVect::Zero, ebcellfact);
  m_faceDivCell.define(m_gridsFine,             1,4*IntVect::Unit, ebfluxfact);
  m_dropOrder.define( m_gridsFine);
  m_johanStencil.define(m_gridsFine);
  m_distanceAlongLine.define(m_gridsFine);

  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      BaseIVFAB<VoFStencil>& curJohanStencil      = m_johanStencil[dit()];
      BaseIVFAB<char>&       curDropOrder         = m_dropOrder[dit()];
      BaseIVFAB<Real>&       curDistanceAlongLine = m_distanceAlongLine[dit()];

      EBGraph ebgraph =  m_ebislFine[dit()].getEBGraph();
      IntVectSet ivsIrreg = m_ebislFine[dit()].getIrregIVS(m_gridsFine.get(dit()));
      curJohanStencil.define(     ivsIrreg,ebgraph, 2);
      curDistanceAlongLine.define(ivsIrreg,ebgraph, 2);
      curDropOrder.define(        ivsIrreg,ebgraph, 1);

      for (VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
        {
          Vector<VoFStencil> pointStencils;
          Vector<Real>       distanceAlongLine;
          bool dropOrder = false;
          EBArith::johanStencil(dropOrder,
                                pointStencils,
                                distanceAlongLine,
                                vofit(), m_ebislFine[dit()], m_dxFine,
                                (*m_eblgFine.getCFIVS())[dit()]);

          curDropOrder(vofit(), 0) = dropOrder;
          if (!dropOrder)
            {
              for (int ipoint = 0; ipoint < 2; ipoint++)
                {
                  curJohanStencil(vofit(), ipoint)      = pointStencils[ipoint];
                  curDistanceAlongLine(vofit(), ipoint) = distanceAlongLine[ipoint];
                }
            }
        }
    }
  m_tensorCFI  = NULL;
  if (m_hasCoarser)
    {
      Real dx = a_dxFine[0];
      m_domainCoar  =  coarsen(m_domainFine, m_refRat);
      m_tensorCFI= new EBTensorCFInterp(m_gridsFine, *m_gridsCoarPtr,
                                        m_ebislFine, *m_ebislCoarPtr,
                                        m_domainCoar, m_refRat, SpaceDim, dx,
                                        *m_eblgFine.getCFIVS(), a_ebisPtr);

    }

  m_lambdaScale = a_lambdaScale;
  fillLambda();
}
/***********************/
void
EBGradDivFilter::
filter(LevelData<EBCellFAB>&       a_velFine,
       const LevelData<EBFluxFAB>& a_fluxVelFine,
       const LevelData<EBCellFAB>* a_velCoar,
       bool a_lowOrderOneSidedGrad,
       bool a_noExtrapToCovered)
{
  CH_TIME("EBGradDivFilter::filter");
  CH_assert(a_velFine.nComp() == SpaceDim);
  Interval velInterv(0, SpaceDim-1);
  //do coarse-fine interpolation if necessary
  if (m_hasCoarser)
    {
      CH_assert(a_velCoar != NULL);
      CH_assert(a_velCoar->nComp() == SpaceDim);
      m_tensorCFI->coarseFineInterp(a_velFine, m_gradVel, *a_velCoar);
    }

  //fill in ghost cells
  a_velFine.exchange(velInterv);
  //compute gradient of velocity for parts NOT in ghost cells
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      gradVel(m_gradVel[dit()],
              a_velFine[dit()],
              m_gridsFine.get(dit()),
              m_ebislFine[dit()],
              a_lowOrderOneSidedGrad);
    }
  m_gradVel.exchange();

  int ibox = 0;
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      EBFluxFAB& faceDiv = m_faceDivCell[dit()];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          faceDivergence(faceDiv[idir],
                         m_gradVel[dit()],
                         a_velFine[dit()],
                         a_fluxVelFine[dit()],
                         m_gridsFine.get(dit()),
                         m_ebislFine[dit()],
                         idir);
        }
      ibox++;
    }
  m_faceDivCell.exchange();

  //apply filter grid by grid
  //the lambda gets multiplied in on the fly
  //(the true tells grad div to do this multiplication)
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB lambdaGradDivVel(m_ebislFine[dit()], m_gridsFine.get(dit()), SpaceDim);
      gradDiv(lambdaGradDivVel,
              m_gradVel[dit()],
              a_velFine[dit()],
              a_fluxVelFine[dit()],
              m_faceDivCell[dit()],
              m_gridsFine.get(dit()),
              m_ebislFine[dit()],
              dit(),
              true,
              a_noExtrapToCovered);
      a_velFine[dit()] += lambdaGradDivVel;
    }
}
/***********************/
void
EBGradDivFilter::
gradDiv(LevelData<EBCellFAB>&       a_gradDivVel,
        LevelData<EBCellFAB>&       a_velFine,
        const LevelData<EBFluxFAB>& a_fluxVelFine,
        const LevelData<EBCellFAB>* a_velCoar,
        bool a_lowOrderOneSidedGrad,
        bool a_noExtrapToCovered)
{
  CH_TIME("EBGradDivFilter::gradDiv");
  CH_assert(a_velFine.nComp() == SpaceDim);
  Interval velInterv(0, SpaceDim-1);
  //do coarse-fine interpolation if necessary
  if (m_hasCoarser)
    {
      CH_assert(a_velCoar != NULL);
      CH_assert(a_velCoar->nComp() == SpaceDim);
      m_tensorCFI->coarseFineInterp(a_velFine, m_gradVel, *a_velCoar);
    }
  //fill in ghost cells
  a_velFine.exchange(velInterv);
  //compute gradient of velocity for parts NOT in ghost cells
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      m_gradVel[dit()].setVal(0.);
      gradVel(m_gradVel[dit()],
              a_velFine[dit()],
              m_gridsFine.get(dit()),
              m_ebislFine[dit()],
              a_lowOrderOneSidedGrad);
    }
  m_gradVel.exchange();

  int ibox = 0;
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      EBFluxFAB& fluxVel = m_faceDivCell[dit()];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          faceDivergence(fluxVel[idir],
                         m_gradVel[dit()],
                         a_velFine[dit()],
                         a_fluxVelFine[dit()],
                         m_gridsFine.get(dit()),
                         m_ebislFine[dit()],
                         idir);
        }
      ibox++;
    }
  m_faceDivCell.exchange();

  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      //the false tells gradDiv to not premultiply by lambda
      gradDiv(a_gradDivVel[dit()],
              m_gradVel[dit()],
              a_velFine[dit()],
              a_fluxVelFine[dit()],
              m_faceDivCell[dit()],
              m_gridsFine.get(dit()),
              m_ebislFine[dit()],
              dit(),
              false,
              a_noExtrapToCovered);
    }
}

/***********************/
void
EBGradDivFilter::
gradDiv(EBCellFAB&             a_gradDivVel,
        const EBCellFAB&       a_gradVel,
        const EBCellFAB&       a_velFine,
        const EBFluxFAB&       a_fluxVelFine,
        const EBFluxFAB&       a_divUFaceCent,
        const Box&             a_gridFine,
        const EBISBox&         a_ebisBoxFine,
        const DataIndex&       a_dataIndex,
        bool a_multiplyByLambda,
        bool a_noExtrapToCovered)
{
  cellGradient(a_gradDivVel, a_gradVel, a_velFine, a_fluxVelFine,
               a_divUFaceCent, a_gridFine, a_ebisBoxFine,
               a_dataIndex, a_multiplyByLambda, a_noExtrapToCovered);
}
/***********************/
Real
EBGradDivFilter::
getDomainDivergence(const EBCellFAB&       a_gradVel,
                    const EBCellFAB&       a_vel,
                    const EBFluxFAB&       a_fluxVel,
                    const Box&             a_grid,
                    const EBISBox&         a_ebisBox,
                    const int&             a_faceDir,
                    const FaceIndex&       a_face,
                    const Side::LoHiSide&  a_side)
{
  CH_TIME("EBGradDivFilter::getDomainDivergence");
  Real divVel;
  CH_assert(a_fluxVel.nComp() == 1);
  Real velB = a_fluxVel[a_faceDir](a_face, 0);
  //compute normal derivative.
  //need second-order derivative at the boundary so
  //given u_b, u_i, u_i+1,
  //ux_b = (9u_i - u_i+1 - 8u_b)/(6 dx)
  VolIndex  vofNear, vofStart;
  Real valNear=0;
  Real valStart =0;
  bool hasVoFNear;
  vofStart = a_face.getVoF(flip(a_side));
  valStart = a_vel(vofStart, a_faceDir);
  int iflipsign = sign(flip(a_side));
  Real rsign = Real(iflipsign);

  Vector<FaceIndex> facesNear = a_ebisBox.getFaces(vofStart, a_faceDir, flip(a_side));
  hasVoFNear = (facesNear.size() == 1);
  if (hasVoFNear)
    {
      vofNear = facesNear[0].getVoF(flip(a_side));
      valNear = a_vel(vofNear, a_faceDir);
    }

  Real hsign = rsign*m_dxFine[a_faceDir];
  if (hasVoFNear)
    {
      divVel  = (2.25*(valStart-velB)  - 0.25*(valNear-velB))/(0.75*hsign);
    }
  else
    {
      //if there is no farther away value, drop order of derivative
      divVel = (valStart - velB)/(0.5*hsign);
    }

  //now add in tangential derivatives of tangential velocities
  //using extrapolated gradients
  for (int tanDir = 0; tanDir< SpaceDim; tanDir++)
    {
      if (tanDir != a_faceDir)
        {
          Real gradVel = 0.0;
          int gradcomp = getGradComp(tanDir, tanDir);
          if (hasVoFNear)
            {
              Real gradNear = a_gradVel(vofNear,   gradcomp);
              Real gradStart = a_gradVel(vofStart, gradcomp);
              gradVel = 1.5*gradStart - 0.5*gradNear;
            }
          else
            {
              Real gradStart = a_gradVel(vofStart, gradcomp);
              gradVel = gradStart;
            }
          divVel += gradVel;
        }
    }

  return divVel;

}
/***********************/
void
EBGradDivFilter::
faceDivergence(EBFaceFAB&             a_divVel,
               const EBCellFAB&       a_gradVel,
               const EBCellFAB&       a_vel,
               const EBFluxFAB&       a_fluxVel,
               const Box&             a_grid,
               const EBISBox&         a_ebisBox,
               const int&             a_faceDir)
{
  CH_TIME("EBGradDivFilter::faceDivergence");
  CH_assert(a_divVel.nComp() == 1);
  CH_assert(a_vel.nComp() == SpaceDim);

  a_divVel.setVal(0.0);
  BaseFab<Real>& regDivVel       =  a_divVel.getSingleValuedFAB();
  const BaseFab<Real>& regVel    =     a_vel.getSingleValuedFAB();
  const BaseFab<Real>& regGradVel= a_gradVel.getSingleValuedFAB();

  Box interiorFaceBox = a_grid;
  interiorFaceBox.grow(a_faceDir, 1);
  interiorFaceBox &= m_domainFine;
  interiorFaceBox.grow(a_faceDir, -1);
  interiorFaceBox.surroundingNodes(a_faceDir);

  for (int divDir = 0; divDir < SpaceDim; divDir++)
    {
      FORT_EBGDFFACEDIVINCR(CHF_FRA1(regDivVel, 0),
                            CHF_FRA(regVel),
                            CHF_FRA(regGradVel),
                            CHF_BOX(interiorFaceBox),
                            CHF_REAL(m_dxFine[divDir]),
                            CHF_INT(a_faceDir),
                            CHF_INT(divDir));
    }

  IntVectSet irregIVS = a_ebisBox.getIrregIVS(a_grid);
  FaceStop::WhichFaces stopCrit;
  if (m_domainFine.isPeriodic(a_faceDir))
    {
      stopCrit = FaceStop::SurroundingWithBoundary;
    }
  else
    {
      stopCrit = FaceStop::SurroundingNoBoundary;
    }
  for (FaceIterator faceit(irregIVS, a_ebisBox.getEBGraph(), a_faceDir,
                          stopCrit);
      faceit.ok(); ++faceit)
    {
      Real divVal = 0.0;
      for (int divDir=0; divDir<SpaceDim; divDir++)
        {
          if (divDir == a_faceDir)
            {
              divVal += (a_vel(faceit().getVoF(Side::Hi), a_faceDir) -
                         a_vel(faceit().getVoF(Side::Lo), a_faceDir))/m_dxFine[divDir];
            }
          else
            {
              // take average of cell-centered tangential derivatives
              // and increment div
              //so this is partial (vel_divdir)/partial (x_divdir)
              int velcomp = divDir;
              int gradcomp = getGradComp(velcomp, divDir);

              divVal += 0.5*(a_gradVel(faceit().getVoF(Side::Hi), gradcomp) +
                             a_gradVel(faceit().getVoF(Side::Lo), gradcomp));
            }
        }

      a_divVel(faceit(), 0) = divVal;
    } // end loop over interior irregular faces

  if (!m_domainFine.isPeriodic(a_faceDir))
    {
      //set the boundary face divergence to an extrapolation of the neighboring face divergences
      for (SideIterator sit; sit.ok(); ++sit)
        {
          Box bndryBox = adjCellBox(a_grid, a_faceDir, sit(), 1);
          int ishift = -sign(sit());
          bndryBox.shift(a_faceDir, ishift);
          IntVectSet bndryIVS(bndryBox);
          for (FaceIterator faceit(bndryIVS, a_ebisBox.getEBGraph(), a_faceDir,
                                  FaceStop::AllBoundaryOnly);
              faceit.ok(); ++faceit)
            {
              Real faceDiv = getDomainDivergence(a_gradVel,
                                                 a_vel,
                                                 a_fluxVel,
                                                 a_grid,
                                                 a_ebisBox,
                                                 a_faceDir,
                                                 faceit(),
                                                 sit());
              a_divVel(faceit(), 0) = faceDiv;
            }
        }
    }
}
/**/
/***********************/
void
EBGradDivFilter::
gradVel(EBCellFAB&             a_gradVel,
        const  EBCellFAB&      a_vel,
        const Box&             a_grid,
        const EBISBox&         a_ebisBox,
        bool a_lowOrderOneSide)
{
  CH_TIME("EBGradDivFilter::gradVel");
  CH_assert(a_gradVel.nComp() == SpaceDim*SpaceDim);
  CH_assert(a_vel.nComp() == SpaceDim);

  int iLowOrderOneSide = 0;
  if (a_lowOrderOneSide)
    {
      iLowOrderOneSide = 1;
    }

  for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
    {
      Box loBox, hiBox, centerBox;
      int hasLo, hasHi;
      EBArith::loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, m_domainFine, a_grid, derivDir);
      for (int velDir = 0; velDir < SpaceDim; velDir++)
        {
          BaseFab<Real>&   regGradVel  = a_gradVel.getSingleValuedFAB();
          const BaseFab<Real>& regVel  =     a_vel.getSingleValuedFAB();
          int gradcomp = getGradComp(velDir, derivDir);

          FORT_EBGDFGRADVEL(CHF_FRA1(regGradVel, gradcomp),
                            CHF_CONST_FRA1(regVel, velDir),
                            CHF_BOX(loBox),
                            CHF_INT(hasLo),
                            CHF_BOX(hiBox),
                            CHF_INT(hasHi),
                            CHF_BOX(centerBox),
                            CHF_CONST_REAL(m_dxFine[derivDir]),
                            CHF_CONST_INT(derivDir),
                            CHF_CONST_INT(iLowOrderOneSide));

          IntVectSet ivsIrreg = a_ebisBox.getIrregIVS(a_grid);
          if (!a_lowOrderOneSide)
            {
              ivsIrreg.grow(1);
              ivsIrreg &= a_grid;
            }
          for (VoFIterator vofit(ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              bool hasHi, hasLo;
              VolIndex vofHi, vofLo;
              Vector<FaceIndex> facesHi = a_ebisBox.getFaces(vof, derivDir, Side::Hi);
              Vector<FaceIndex> facesLo = a_ebisBox.getFaces(vof, derivDir, Side::Lo);

              hasHi = (facesHi.size() == 1) && (!facesHi[0].isBoundary());
              hasLo = (facesLo.size() == 1) && (!facesLo[0].isBoundary());
              if (hasLo)
                {
                  vofLo = facesLo[0].getVoF(Side::Lo);
                }
              if (hasHi)
                {
                  vofHi = facesHi[0].getVoF(Side::Hi);
                }

              Real gradVal;
              if (hasHi && hasLo)
                {
                  gradVal = (a_vel(vofHi, velDir) - a_vel(vofLo, velDir))/(2.*m_dxFine[derivDir]);
                }
              else if (hasHi || hasLo)
                {
                  //do one-sided diff
                  CH_assert(!(hasHi && hasLo));
                  Side::LoHiSide side;
                  VolIndex vofNear;
                  if (hasLo)
                    {
                      side = Side::Lo;
                      vofNear = vofLo;
                    }
                  else
                    {
                      side = Side::Hi;
                      vofNear = vofHi;
                    }
                  int isign = sign(side);

                  Vector<FaceIndex> facesFar = a_ebisBox.getFaces(vofNear, derivDir, side);
                  bool hasVoFFar = (facesFar.size()==1) && (!facesFar[0].isBoundary());
                  Real valStart = a_vel(vof,     velDir);
                  Real valNear  = a_vel(vofNear, velDir);
                  if (hasVoFFar && !a_lowOrderOneSide)
                    {
                      VolIndex vofFar = facesFar[0].getVoF(side);
                      Real valFar = a_vel(vofFar, velDir);
                      gradVal = (4*(valNear - valStart)- (valFar-valStart))/(2.*isign*m_dxFine[derivDir]);
                    }
                  else
                    {
                      //do not have another point for stencil or specified low order one-sided.
                      //drop order of derivative
                      gradVal = (valNear - valStart)/(isign*m_dxFine[derivDir]);
                    }
                }
              else
                {
                  //only have one point, can only set gradient to zero
                  gradVal = 0.0;
                }

              a_gradVel(vof, gradcomp) = gradVal;
            }
        }
    }
}
/***********************/
void
EBGradDivFilter::
getAreaFracs(Vector<FaceIndex> a_facesLo[SpaceDim],
             Vector<FaceIndex> a_facesHi[SpaceDim],
             bool              a_hasFacesLo[SpaceDim],
             bool              a_hasFacesHi[SpaceDim],
             RealVect&         a_areaFracLo,
             RealVect&         a_areaFracHi,
             const VolIndex&   a_vof,
             const EBISBox&    a_ebisBox)
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real areaSumLo = 0;
      Real areaSumHi = 0;
      a_facesLo[idir] = a_ebisBox.getFaces(a_vof, idir, Side::Lo);
      a_facesHi[idir] = a_ebisBox.getFaces(a_vof, idir, Side::Hi);
      a_hasFacesLo[idir] = ( ( a_facesLo[idir].size() > 0) &&
                             (!a_facesLo[idir][0].isBoundary()));
      a_hasFacesHi[idir] = ( ( a_facesHi[idir].size() > 0) &&
                             (!a_facesHi[idir][0].isBoundary()));
      //want boundary faces to look covered
      if (a_hasFacesLo[idir])
        {
          for (int iface = 0; iface < a_facesLo[idir].size(); iface++)
            {
              areaSumLo += a_ebisBox.areaFrac(a_facesLo[idir][iface]);
            }
        }
      //want boundary faces to look covered
      if (a_hasFacesHi[idir])
        {
          for (int iface = 0; iface < a_facesHi[idir].size(); iface++)
            {
              areaSumHi += a_ebisBox.areaFrac(a_facesHi[idir][iface]);
            }
        }

      a_areaFracLo[idir] = areaSumLo;
      a_areaFracHi[idir] = areaSumHi;
    }
}
/***********************/
void
EBGradDivFilter::
fillLambda()
{
  CH_TIME("EBGradDivFilter::fillLambda");
  Real lambdaVal = 1.0e10;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real lambdaDir  = m_lambdaScale*m_dxFine[idir]*m_dxFine[idir];
      lambdaVal = Min(lambdaVal, lambdaDir);
    }
  for (DataIterator dit = m_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      m_lambda[dit()].setVal(lambdaVal);
    }
}
/***********************/
void
EBGradDivFilter::
cellGradient(EBCellFAB&             a_gradDiv,
             const EBCellFAB&       a_gradVel,
             const  EBCellFAB&      a_vel,
             const EBFluxFAB&       a_fluxVel,
             const EBFluxFAB&       a_div,
             const Box&             a_grid,
             const EBISBox&         a_ebisBox,
             const DataIndex&       a_dataIndex,
             bool a_multiplyByLambda,
             bool a_noExtrapToCovered)
{
  CH_TIME("EBGradDivFilter::cellGradient");
  int icomp = 0;
  EBCellFAB& lambdaFAB = m_lambda[a_dataIndex];
  IntVectSet irregIVS = a_ebisBox.getIrregIVS(a_grid);

  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      const EBFaceFAB&  faceDiv       =  a_div[faceDir];
      const BaseFab<Real>& regFaceDiv =   faceDiv.getSingleValuedFAB();
      const BaseFab<Real>& regLambda =   lambdaFAB.getSingleValuedFAB();
      BaseFab<Real>& regGradDiv       = a_gradDiv.getSingleValuedFAB();

      int imultByLambda = 0;
      if (a_multiplyByLambda)
        {
          imultByLambda = 1;
        }
      FORT_EBGDFCELLGRAD(CHF_FRA1(regGradDiv, faceDir),
                         CHF_CONST_FRA1(regFaceDiv, icomp),
                         CHF_CONST_FRA1(regLambda, icomp),
                         CHF_BOX(a_grid),
                         CHF_CONST_REAL(m_dxFine[faceDir]),
                         CHF_CONST_INT(faceDir),
                         CHF_CONST_INT(imultByLambda));

      //nonconservative gradient
      for (VoFIterator vofit(irregIVS, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real gradVal = 0.0;

          Vector<FaceIndex> facesHi = a_ebisBox.getFaces(vof, faceDir, Side::Hi);
          Vector<FaceIndex> facesLo = a_ebisBox.getFaces(vof, faceDir, Side::Lo);

          bool zeroGrad = a_noExtrapToCovered && ((facesHi.size() == 0) || (facesLo.size() == 0));
          if (zeroGrad)
            {
              gradVal = 0;
            }
          else
            {
              Real valHi= 0;
              if (facesHi.size() > 0)
                {
                  for (int iface = 0; iface < facesHi.size(); iface++)
                    {
                      valHi += faceDiv(facesHi[iface], icomp);
                    }
                  valHi /= facesHi.size();
                }
              else
                {
                  valHi = ccpGetCoveredExtrapValue(vof, faceDir, Side::Hi,
                                                   faceDiv, a_ebisBox, a_grid,
                                                   m_domainFine, m_dxFine, icomp);
                }

              Real valLo= 0;
              if (facesLo.size() > 0)
                {
                  for (int iface = 0; iface < facesLo.size(); iface++)
                    {
                      valLo += faceDiv(facesLo[iface], 0);
                    }
                  valLo /= facesLo.size();
                }
              else
                {
                  valLo = ccpGetCoveredExtrapValue(vof, faceDir, Side::Lo ,
                                                   faceDiv, a_ebisBox, a_grid,
                                                   m_domainFine, m_dxFine, icomp);

                }
              gradVal = (valHi-valLo)/(m_dxFine[faceDir]);

              //multiply the lambda in since we have the area fractions lying about
              if (a_multiplyByLambda)
                {
                  Real lambda = lambdaFAB(vof, 0);
                  gradVal *= lambda;
                }
            }

          a_gradDiv(vof, faceDir) = gradVal;
        }
    }
}
/***********************/

#include "NamespaceFooter.H"
