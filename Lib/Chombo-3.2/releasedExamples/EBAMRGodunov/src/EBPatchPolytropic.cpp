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

#include "EBDebugOut.H"
#include "DebugOut.H"
#include "EBPatchPolytropic.H"
#include "EBPatchPolytropicF_F.H"
#include "EBPatchGodunovF_F.H"
#include "EBLGIntegrator.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "EBLoHiCenter.H"
#include "EBPatchGodunov.H"
#include "Stencils.H"
#include "EBArith.H"
#include "EBAMRIO.H"
#include "PolyGeom.H"
#include "TensorCFInterp.H" //just for gradIndex
#include "ParmParse.H"
#include <cstdio>

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
void fillValues(Vector<Real>&            a_prim,
                const Vector<FaceIndex>& a_faces,
                const EBFaceFAB&         a_facePrim,
                const BaseIVFAB<Real>&   a_coveredPrim,
                const VolIndex           a_vof)
{
  CH_assert(a_prim.size() == a_facePrim.nComp());
  CH_assert(a_prim.size() == a_coveredPrim.nComp());

  for (int ivar = 0; ivar < a_prim.size(); ivar++)
    {
      if (a_faces.size() == 0)
        {
          a_prim[ivar] = a_coveredPrim(a_vof, ivar);
        }
      else
        {
          Real primVal = 0;
          for (int iface = 0; iface < a_faces.size(); iface++)
            {
              primVal += a_facePrim(a_faces[iface], ivar);
            }
          if (a_faces.size() > 1) primVal /= Real(a_faces.size());
          a_prim[ivar] = primVal;
        }
    }
}
//-----------------------------------------------------------------------
void
EBPatchPolytropic::
getRhoUDotDele(EBCellFAB&          a_rhoUDotDele,
               EBFluxFAB&          a_facePrim,
               BaseIVFAB<Real>     a_coveredPrimLo[SpaceDim],
               BaseIVFAB<Real>     a_coveredPrimHi[SpaceDim],
               const Box&          a_box,
               const IntVectSet&   a_ivsSmall)
{
  CH_TIME("EBPatchPolytropic::getUDotDele");

  //this is necessary since we are doing increments from here
  a_rhoUDotDele.setVal(0.);


  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      FORT_INCRRHOUDOTDELE(CHF_FRA1(a_rhoUDotDele.getSingleValuedFAB(), 0),
                           CHF_CONST_FRA(a_facePrim[faceDir].getSingleValuedFAB()),
                           CHF_CONST_REAL(m_specHeat),
                           CHF_CONST_REAL(m_dx[faceDir]),
                           CHF_CONST_INT(faceDir),
                           CHF_BOX(a_box));
    }

  //update the irregular vofs
  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          Real irregRHSVelo = 0;
          for (int diffDir = 0; diffDir < SpaceDim; diffDir++)
            {
              Vector<Real> primLo(QNUM);
              Vector<Real> primHi(QNUM);
              Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, diffDir, Side::Lo);
              Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, diffDir, Side::Hi);
              fillValues(primLo, facesLo, a_facePrim[diffDir], a_coveredPrimLo[diffDir], vof);
              fillValues(primHi, facesHi, a_facePrim[diffDir], a_coveredPrimHi[diffDir], vof);
              Real eLo, eHi;
              FORT_POINTGETTEMPFROMPRIM(CHF_REAL(eLo),CHF_CONST_VR(primLo), CHF_REAL(m_specHeat));
              FORT_POINTGETTEMPFROMPRIM(CHF_REAL(eHi),CHF_CONST_VR(primHi), CHF_REAL(m_specHeat));
              eLo *= m_specHeat;
              eHi *= m_specHeat;

              int dirc = QVELX + diffDir;
              Real rhoave   =   0.5*(primHi[QRHO ] +  primLo[QRHO ]);
              Real ufaceave =   0.5*(primHi[dirc ] +  primLo[dirc ]);

              irregRHSVelo  += rhoave*ufaceave*(eHi - eLo)/m_dx[diffDir];
            }
          a_rhoUDotDele(vof, 0) = irregRHSVelo;
        }
    }//end loop over irreg vofs
}
////////
void
EBPatchPolytropic::
getGradU(EBCellFAB&          a_gradU,
         EBFluxFAB&          a_facePrim,
         BaseIVFAB<Real>     a_coveredPrimLo[SpaceDim],
         BaseIVFAB<Real>     a_coveredPrimHi[SpaceDim],
         const Box&          a_box,
         const IntVectSet&   a_ivsSmall)
{
      //always difference in the facedir direction
  for (int velDir = 0; velDir < SpaceDim; velDir++)
    {
      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          int gradIndex = TensorCFInterp::gradIndex(velDir, faceDir);

          FORT_EBPPGRADVEL(CHF_FRA(a_gradU.getSingleValuedFAB()),
                           CHF_CONST_FRA(a_facePrim[faceDir].getSingleValuedFAB()),
                           CHF_CONST_REAL(m_dx[faceDir]),
                           CHF_CONST_INT(faceDir),
                           CHF_CONST_INT(velDir),
                           CHF_CONST_INT(gradIndex),
                           CHF_BOX(a_box));

          //update the irregular vofs
          for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
            {
              const VolIndex& vof = m_irregVoFs[ivof];
              if (a_box.contains(vof.gridIndex()))
                {
                  Vector<Real> primLo(QNUM);
                  Vector<Real> primHi(QNUM);
                  Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, faceDir, Side::Lo);
                  Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, faceDir, Side::Hi);
                  fillValues(primLo, facesLo, a_facePrim[faceDir], a_coveredPrimLo[faceDir], vof);
                  fillValues(primHi, facesHi, a_facePrim[faceDir], a_coveredPrimHi[faceDir], vof);

                  a_gradU(vof, gradIndex) = (primHi[QVELX+velDir] - primLo[QVELX+velDir])/m_dx[faceDir];
                }
            }//end loop over irreg vofs
        }
    }

}

////////
void
EBPatchPolytropic::
getVeloHalf(EBCellFAB&          a_veloHalf,
            EBFluxFAB&          a_facePrim,
            BaseIVFAB<Real>     a_coveredPrimLo[SpaceDim],
            BaseIVFAB<Real>     a_coveredPrimHi[SpaceDim],
            const Box&          a_box,
            const IntVectSet&   a_ivsSmall)
{
  //always average the normal (facedir) velocities
  for (int velDir = 0; velDir < SpaceDim; velDir++)
    {

      FORT_EBPPVELOHALF(CHF_FRA(a_veloHalf.getSingleValuedFAB()),
                        CHF_CONST_FRA(a_facePrim[velDir].getSingleValuedFAB()),
                        CHF_CONST_INT(velDir),
                        CHF_BOX(a_box));

      //update the irregular vofs
      for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
        {
          const VolIndex& vof = m_irregVoFs[ivof];
          if (a_box.contains(vof.gridIndex()))
            {
              Vector<Real> primLo(QNUM);
              Vector<Real> primHi(QNUM);
              Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, velDir, Side::Lo);
              Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, velDir, Side::Hi);
              fillValues(primLo, facesLo, a_facePrim[velDir], a_coveredPrimLo[velDir], vof);
              fillValues(primHi, facesHi, a_facePrim[velDir], a_coveredPrimHi[velDir], vof);

              a_veloHalf(vof, velDir) = 0.5*(primHi[QVELX+velDir] + primLo[QVELX+velDir]);
            }
        }//end loop over irreg vofs
    }

}
////////
void
EBPatchPolytropic::
getPDivU(EBCellFAB&          a_pDivU,
         EBFluxFAB&          a_facePrim,
         BaseIVFAB<Real>     a_coveredPrimLo[SpaceDim],
         BaseIVFAB<Real>     a_coveredPrimHi[SpaceDim],
         const Box&          a_box,
         const IntVectSet&   a_ivsSmall)
{
  //have to do this because we are incrementing
  a_pDivU.setVal(0.);
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      FORT_INCRPDIVU(CHF_FRA1(a_pDivU.getSingleValuedFAB(), 0),
                     CHF_CONST_FRA(a_facePrim[faceDir].getSingleValuedFAB()),
                     CHF_CONST_REAL(m_dx[faceDir]),
                     CHF_CONST_INT(faceDir),
                     CHF_BOX(a_box));
    }

  //update the irregular vofs
  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          Real irregVal = 0;
          for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
            {
              int vind = QVELX + faceDir;

              Vector<Real> primLo(QNUM);
              Vector<Real> primHi(QNUM);
              Vector<FaceIndex> facesLo = m_ebisBox.getFaces(vof, faceDir, Side::Lo);
              Vector<FaceIndex> facesHi = m_ebisBox.getFaces(vof, faceDir, Side::Hi);
              fillValues(primLo, facesLo, a_facePrim[faceDir], a_coveredPrimLo[faceDir], vof);
              fillValues(primHi, facesHi, a_facePrim[faceDir], a_coveredPrimHi[faceDir], vof);

              Real pave  = 0.5*(primLo[QPRES] + primHi[QPRES]);
              Real udiff =     (primLo[vind ] + primHi[vind ])/m_dx[faceDir];
              irregVal = irregVal + pave*udiff;
            }
          a_pDivU(vof, 0) = irregVal;
        }
    }//end loop over irreg vofs
}


//-----------------------------------------------------------------------
void
EBPatchPolytropic::
primitivesAndDivergences(EBCellFAB&          a_nonConsDivF,
                         EBCellFAB&          a_consState,
                         EBFluxFAB&          a_facePrim,
                         BaseIVFAB<Real>     a_coveredPrimMinu[SpaceDim],
                         BaseIVFAB<Real>     a_coveredPrimPlus[SpaceDim],
                         Vector<VolIndex>    a_coveredFaceMinu[SpaceDim],
                         Vector<VolIndex>    a_coveredFacePlus[SpaceDim],
                         EBFluxFAB&          a_flux,
                         BaseIVFAB<Real>&    a_ebIrregFlux,
                         BaseIVFAB<Real>&    a_nonConservativeDivergence,
                         const EBCellFAB&    a_flattening,
                         const EBCellFAB&    a_source,
                         const Box&          a_box,
                         const IntVectSet&   a_ivsSmall,
                         const DataIndex&    a_dit,
                         bool                a_verbose)
{
  CH_TIME("EBPatchGodunov::regularUpdate");
  CH_assert(isDefined());
  int numCons = numConserved();

  CH_assert(a_flux.nComp() >= numConserved());
  CH_assert(a_consState.nComp() >= numConserved());
  setCoveredConsVals(a_consState);

  EBCellFAB slopePrim[SpaceDim];
  EBCellFAB slopeNLim[SpaceDim];

  computeFluxes(a_flux,
                m_coveredFluxMinuG4, m_coveredFluxPlusG4,
                a_coveredFaceMinu,   a_coveredFacePlus,
                m_primState,    slopePrim, slopeNLim,
                a_flattening, a_consState, a_source,
                a_box, a_dit,a_verbose);

  //now I happen to know that m_primgdnv holds faceprim
  //and m_extendedStateMinuG4 holds  coveredPrimMinu and so on
  a_facePrim.define(m_ebisBox, a_box, numPrimitives());
  Interval inter(0, numPrimitives()-1);
  a_facePrim.copy(a_box, inter, a_box, m_primGdnv, inter);

  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      IntVectSet ivsMinu = m_coveredSetsMinuG4[faceDir] & a_box;
      IntVectSet ivsPlus = m_coveredSetsPlusG4[faceDir] & a_box;
      a_coveredPrimMinu[faceDir].define(ivsMinu, m_ebisBox.getEBGraph(), numPrimitives());
      a_coveredPrimPlus[faceDir].define(ivsPlus, m_ebisBox.getEBGraph(), numPrimitives());
      a_coveredPrimMinu[faceDir].copy(a_box, inter, a_box, m_extendStateMinuG4[faceDir], inter);
      a_coveredPrimPlus[faceDir].copy(a_box, inter, a_box, m_extendStatePlusG4[faceDir], inter);
    }

  //now that we have the fluxes, modify them with
  //artificial viscosity if that is called for.
  //the artificial viscosity correction to the
  //embedded boundary flux is computed inside
  //computeebirregflux
  if (usesArtificialViscosity())
    {
      EBFluxFAB openDivU(m_ebisBox, a_box, 1);

      getFaceDivergence(openDivU,
                        m_primState, slopeNLim,
                        a_box, a_ivsSmall);

      applyArtificialViscosity(a_flux,
                               m_coveredFluxMinuG4,
                               m_coveredFluxPlusG4,
                               m_coveredFaceMinuG4,
                               m_coveredFacePlusG4,
                               a_consState,
                               openDivU,
                               a_box,
                               a_ivsSmall);
    }

  //this is the stable, non-conservative estimate of the solution update
  nonconservativeDivergence(a_nonConsDivF, a_flux,
                            m_coveredFluxMinuG4,
                            m_coveredFluxPlusG4,
                            a_coveredFaceMinu,
                            a_coveredFacePlus,
                            a_box);

  //copy nonconservative div f into sparse output thingy
  IntVectSet ivsIrreg = a_ivsSmall;
  for (VoFIterator vofit(ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      for (int ivar = 0; ivar < numCons; ivar++)
        {
          a_nonConservativeDivergence(vof, ivar) = a_nonConsDivF(vof, ivar);
        }
    }

  //compute irregular boundary flux.  this includes
  //an artificial viscosity modification if appropriate
  computeEBIrregFlux( a_ebIrregFlux, m_primState,
                      slopeNLim, ivsIrreg, a_source);
}
/*****************************/
void
EBPatchPolytropic::
doRZCoords(bool a_doRZCoords)
{
  m_doRZCoords = a_doRZCoords;
}
/*****************************/
void
EBPatchPolytropic::
setSource(EBCellFAB&       a_source,
          const EBCellFAB& a_consState,
          const Box&       a_box)
{
  if (m_doRZCoords)
    {
      setRZSource(a_source,
                  a_consState,
                  a_box);
    }
}
void
EBPatchPolytropic::
setRZSource(EBCellFAB&       a_source,
            const EBCellFAB& a_consState,
            const Box&       a_box)
{
  a_source.setVal(0.);
  const BaseFab<Real>& regConsState = a_consState.getSingleValuedFAB();
  BaseFab<Real>& regSource          =    a_source.getSingleValuedFAB();
  FORT_SETSOURCERZ(CHF_FRA(regSource),
                   CHF_CONST_FRA(regConsState),
                   CHF_CONST_REAL(m_dx[0]),
                   CHF_BOX(a_box));

  //pointwise operation so just have to do the multi valued cells
  //as an irregular iteration
  IntVectSet multiIVS = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(multiIVS, m_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();
      Real dense = a_consState(vof, CRHO);
      Real energy = a_consState(vof, CENG);
      RealVect momen;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          momen[idir] = a_consState(vof, CMOMX + idir);
        }
      Real radius = (Real(iv[0]) + 0.5)*m_dx[0];
      Real densitySource, pressureSource;
      FORT_POINTSETSOURCERZ(CHF_REAL(densitySource),
                            CHF_REAL(pressureSource),
                            CHF_CONST_REAL(dense),
                            CHF_CONST_REALVECT(momen),
                            CHF_CONST_REAL(energy),
                            CHF_CONST_REAL(radius));
      //only the density and energy equations have non-zero source
      //so set the others to zero
      for (int ivar = 0; ivar < numConserved(); ivar++)
        {
          a_source(vof, ivar) =  0.0;
        }
      a_source(vof, QRHO) = densitySource;
      a_source(vof, QPRES) = pressureSource;
    }
}
/*****************************/
void
EBPatchPolytropic::
assembleFluxReg(EBFaceFAB&       a_fluxRegFlux,
                const EBFaceFAB& a_godunovFlux,
                const int&       a_idir,
                const Box&       a_cellBox)
{
  CH_assert(a_idir==a_fluxRegFlux.direction());
  CH_assert(a_idir==a_godunovFlux.direction());
  CH_assert(a_fluxRegFlux.getCellRegion().contains(a_cellBox));
  CH_assert(a_godunovFlux.getCellRegion().contains(a_cellBox));
  CH_assert(a_godunovFlux.nComp() == numFluxes());
  CH_assert(a_fluxRegFlux.nComp() == numConserved());
  CH_assert(SpaceDim==2);

  int rdir = 0;
  BaseFab<Real>& regFluxRegFlux       = a_fluxRegFlux.getSingleValuedFAB();
  const BaseFab<Real>& regGodunovFlux = a_godunovFlux.getSingleValuedFAB();
  const Box& faceRegion = regGodunovFlux.box();
  CH_assert(regFluxRegFlux.box().contains(faceRegion));

  FORT_FLUXASSEMBLE(CHF_FRA(regFluxRegFlux),
                    CHF_CONST_FRA(regGodunovFlux),
                    CHF_CONST_REAL(m_dx[0]),
                    CHF_CONST_INT(a_idir),
                    CHF_BOX(faceRegion));

  IntVectSet ivsirreg = m_ebisBox.getIrregIVS(a_cellBox);
  FaceStop::WhichFaces facestop = FaceStop::SurroundingWithBoundary;
  for (FaceIterator faceit(ivsirreg, m_ebisBox.getEBGraph(), a_idir, facestop);
      faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      Real rind = face.gridIndex(Side::Hi)[rdir];
      Real faceRad;
      if (a_idir == rdir)
        {
          faceRad = m_dx[0]*rind;
        }
      else
        {
          Real rcent = m_ebisBox.centroid(face)[rdir];
          faceRad = m_dx[0]*(rind + 0.5 + rcent);
        }

      for (int ivar=0; ivar < numConserved(); ivar++)
        {
          Real flux = faceRad*a_godunovFlux(face, ivar);
          //add in pressure term
          if (((ivar == CMOMX) && (a_idir==0)) ||
             ((ivar == CMOMY) && (a_idir==1)))
            {
              flux += faceRad*a_godunovFlux(face, CPRES);
            }
          a_fluxRegFlux(face, ivar) = flux;
        }
    }
}
/*****************************/
void
EBPatchPolytropic::
assembleFluxIrr(BaseIFFAB<Real>&       a_fluxRegFlux,
                const BaseIFFAB<Real>& a_godunovFlux,
                const int&             a_idir,
                const Box&             a_cellBox,
                const IntVectSet&      a_set)
{
  int rdir = 0;
  IntVectSet ivsirreg = m_ebisBox.getIrregIVS(a_cellBox);
  FaceStop::WhichFaces facestop = FaceStop::SurroundingWithBoundary;
  for (FaceIterator faceit(ivsirreg, m_ebisBox.getEBGraph(), a_idir, facestop);
      faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      Real rind = face.gridIndex(Side::Hi)[rdir];
      Real faceRad;
      if (a_idir == rdir)
        {
          faceRad = m_dx[0]*rind;
        }
      else
        {
          Real rcent = m_ebisBox.centroid(face)[rdir];
          faceRad = m_dx[0]*(rind + 0.5 + rcent);
        }
      for (int ivar=0; ivar < numConserved(); ivar++)
        {
          Real flux = faceRad*a_godunovFlux(face, ivar);
          //add in pressure term
          if (((ivar == CMOMX) && (a_idir==0)) ||
             ((ivar == CMOMY) && (a_idir==1)))
            {
              flux += faceRad*a_godunovFlux(face, CPRES);
            }
          a_fluxRegFlux(face, ivar) = flux;
        }
    }
}
/*****************************/
void
EBPatchPolytropic::
nonconservativeDivergenRZ(EBCellFAB&             a_divF,
                          const EBFluxFAB&       a_flux,
                          const BaseIVFAB<Real>  a_coveredFluxMinu[SpaceDim],
                          const BaseIVFAB<Real>  a_coveredFluxPlus[SpaceDim],
                          const Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
                          const Vector<VolIndex> a_coveredFacePlus[SpaceDim],
                          const Box&             a_box)
{
  CH_assert(SpaceDim==2);
  CH_assert(isDefined());
  CH_assert(a_divF.nComp() >= numConserved());
  CH_assert(a_flux.nComp() >= numConserved());
  CH_assert(a_flux.getRegion().contains(a_box));
  CH_assert(a_divF.getRegion().contains(a_box));

  //set the divergence initially to zero
  //then loop through directions and increment the divergence
  //with each directions flux difference.
  int ncons = numConserved();
  a_divF.setVal(0.0);
  BaseFab<Real>&       regDivF = a_divF.getSingleValuedFAB();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      //update for the regular vofs in the nonconservative
      //case  works for all single valued vofs.
      const EBFaceFAB& fluxDir = a_flux[idir];

      const BaseFab<Real>& regFluxDir = fluxDir.getSingleValuedFAB();

      /* do the regular vofs */
      /**/
      FORT_DIVERGERZ( CHF_BOX(a_box),
                      CHF_FRA(regDivF),
                      CHF_CONST_FRA(regFluxDir),
                      CHF_CONST_INT(idir),
                      CHF_CONST_INT(ncons),
                      CHF_CONST_REAL(m_dx[0]));
      /**/

    } // end loop over directions

  //update the irregular vofs
  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          //divergence was set in regular update.  we reset it
          // to zero and recalc.
          for (int ivar = 0; ivar < ncons; ivar++)
            {
              //r diffs
              Real rupdate = 0.0;
              {
                int rdir = 0;
                const EBFaceFAB& rflux = a_flux[rdir];
                for (SideIterator sit; sit.ok(); ++sit)
                  {
                    int isign = sign(sit());
                    Vector<FaceIndex> faces = m_ebisBox.getFaces(vof, rdir, sit());
                    int rindex   =  vof.gridIndex()[rdir];
                    Real cellRad = m_dx[0]*(Real(rindex) + 0.5);
                    Real faceRad = cellRad + 0.5*isign*m_dx[0];
                    Real faceupdate = 0.0;
                    for (int iface = 0; iface < faces.size(); iface++)
                      {
                        const FaceIndex& face = faces[iface];
                        Real flux = rflux(face, ivar);
                        faceupdate += isign*flux*faceRad/(cellRad*m_dx[0]);
                      }
                    if (faces.size() > 1)
                      {
                        faceupdate /= Real(faces.size());
                      }
                    rupdate += faceupdate;
                    //add in pressure terms
                    if (ivar == CMOMX)
                      {
                        Real pupdate = 0.0;
                        for (int iface = 0; iface < faces.size(); iface++)
                          {
                            const FaceIndex& face = faces[iface];
                            Real flux = rflux(face, CPRES);
                            pupdate += isign*flux/m_dx[0];
                          }
                        if (faces.size() > 1)
                          {
                            pupdate /= Real(faces.size());
                          }
                        rupdate += pupdate;
                      }
                  }
              }

              //z diffs
              Real zupdate = 0.0;
              {
                int zdir = 1;
                const EBFaceFAB& zflux = a_flux[zdir];
                for (SideIterator sit; sit.ok(); ++sit)
                  {
                    int isign = sign(sit());
                    Vector<FaceIndex> faces = m_ebisBox.getFaces(vof, zdir, sit());
                    Real faceupdate = 0.0;
                    for (int iface = 0; iface < faces.size(); iface++)
                      {
                        const FaceIndex& face = faces[iface];
                        Real flux = zflux(face, ivar);
                        faceupdate += isign*flux/m_dx[0];
                      }
                    if (faces.size() > 1)
                      {
                        faceupdate /= Real(faces.size());
                      }
                    zupdate += faceupdate;
                    //add in pressure terms
                    if (ivar == CMOMY)
                      {
                        Real pupdate = 0.0;
                        for (int iface = 0; iface < faces.size(); iface++)
                          {
                            const FaceIndex& face = faces[iface];
                            Real flux = zflux(face, CPRES);
                            pupdate += isign*flux/m_dx[0];
                          }
                        if (faces.size() > 1)
                          {
                            pupdate /= Real(faces.size());
                          }
                        zupdate += pupdate;
                      }
                  } //end loop over sides
              }
              a_divF(vof, ivar) = rupdate + zupdate;
            }//end loop over variables
        }
    }//end loop over irreg vofs

  //now correct for the covered fluxes
  for (SideIterator sit; sit.ok(); ++sit)
    {
      int isign = sign(sit());
      const BaseIVFAB<Real>*  rcoveredFluxPtr;
      const Vector<VolIndex>* rcoveredFacePtr;
      const BaseIVFAB<Real>*  zcoveredFluxPtr;
      const Vector<VolIndex>* zcoveredFacePtr;
      int rdir = 0;
      int zdir = 1;
      if (sit() == Side::Lo)
        {
          rcoveredFluxPtr = &a_coveredFluxMinu[rdir];
          rcoveredFacePtr = &a_coveredFaceMinu[rdir];
          zcoveredFluxPtr = &a_coveredFluxMinu[zdir];
          zcoveredFacePtr = &a_coveredFaceMinu[zdir];
        }
      else
        {
          rcoveredFluxPtr = &a_coveredFluxPlus[rdir];
          rcoveredFacePtr = &a_coveredFacePlus[rdir];
          zcoveredFluxPtr = &a_coveredFluxPlus[zdir];
          zcoveredFacePtr = &a_coveredFacePlus[zdir];
        }
      //r update
      {
        const BaseIVFAB<Real>&  rflux = *rcoveredFluxPtr;
        const Vector<VolIndex>& rface = *rcoveredFacePtr;
        for (int ivof = 0; ivof < rface.size(); ivof++)
          {
            const VolIndex& vof = rface[ivof];
            int rindex   =  vof.gridIndex()[rdir];
            Real cellRad = m_dx[0]*(Real(rindex) + 0.5);
            Real faceRad = cellRad + 0.5*isign*m_dx[0];

            if (a_box.contains(vof.gridIndex()))
              {
                //face on this side is covered.  use covered flux.
                for (int ivar = 0; ivar < ncons; ivar++)
                  {
                    Real flux =  rflux(vof, ivar);
                    Real update = isign*flux*faceRad/(cellRad*m_dx[0]);
                    //add in pressure terms
                    if (ivar == CMOMX)
                      {
                        flux = rflux(vof, CPRES);
                        Real pupdate  = isign*flux/m_dx[0];
                        update += pupdate;
                      }
                    a_divF(vof, ivar) += update;
                  }
              }
          }
      }

      //z update
      {
        const BaseIVFAB<Real>&  zflux = *zcoveredFluxPtr;
        const Vector<VolIndex>& zface = *zcoveredFacePtr;
        for (int ivof = 0; ivof < zface.size(); ivof++)
          {
            const VolIndex& vof = zface[ivof];
            if (a_box.contains(vof.gridIndex()))
              {
                //face on this side is covered.  use covered flux.
                for (int ivar = 0; ivar < ncons; ivar++)
                  {
                    Real flux =  zflux(vof, ivar);
                    Real update = isign*flux/(m_dx[0]);
                    //add in pressure terms
                    if (ivar == CMOMY)
                      {
                        flux = zflux(vof, CPRES);
                        Real pupdate = isign*flux/m_dx[0];
                        update += pupdate;
                      }
                    a_divF(vof, ivar) += update;
                  }
              }
          }
      }
    }
}
/*****************************/
void
EBPatchPolytropic::
consUndividedDivergenRZ(BaseIVFAB<Real>&       a_divF,
                        const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                        const BaseIVFAB<Real>& a_ebIrregFlux,
                        const IntVectSet&      a_ivs)
{
  CH_assert(SpaceDim==2);
  CH_assert(isDefined());
  CH_assert(a_divF.nComp() >= numConserved());

  //this will break horribly in anisotropic-land
  CH_assert(Abs(m_dx[0]-m_dx[1]) < 1.0e-10);

  int ncons = numConserved();
  a_divF.setVal(0.0);
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Real cellVol, kvol;
      EBArith::getKVolRZNoDx(kvol, cellVol, m_ebisBox, vof);

      CH_assert(kvol >= 0.0);
      CH_assert(kvol <= 1.0);
      //      int rindex   =  vof.gridIndex()[0];
      //     Real cellRad = m_dx[0]*(Real(rindex) + 0.5);
      Real volFrac = m_ebisBox.volFrac(vof);

      Real bndryArea = m_ebisBox.bndryArea(vof);

      Real vofFactor = volFrac/(kvol * cellVol);

      int rdir = 0;
      for (int ivar = 0; ivar < ncons ; ivar++)
        {
          //do r diffs
          Real rupdate = 0.;
          {
            const BaseIFFAB<Real>& rflux = a_centroidFlux[rdir];
            Real faceupdate = 0.0;
            for (SideIterator sit; sit.ok(); ++sit)
              {
                int isign = sign(sit());
                Vector<FaceIndex> faces = m_ebisBox.getFaces(vof, rdir, sit());
                for (int iface = 0; iface < faces.size(); iface++)
                  {
                    const FaceIndex& face = faces[iface];
                    Real areaFrac = m_ebisBox.areaFrac(face);
                    Real faceFlux = rflux(face, ivar);
                    int faceInd =  face.gridIndex(Side::Hi)[rdir];
                    Real faceRad    = Real(faceInd);
                    Real faceFactor = faceRad*areaFrac;

                    faceupdate += isign*vofFactor*faceFactor*faceFlux;

                  }
              }
            rupdate += faceupdate/m_dx[0];
            //lets not forget the pressure
            if (ivar == CMOMX)
              {
                Real pupdate = 0.0;
                for (SideIterator sit; sit.ok(); ++sit)
                  {
                    int isign = sign(sit());
                    Vector<FaceIndex> faces = m_ebisBox.getFaces(vof, rdir, sit());
                    for (int iface = 0; iface < faces.size(); iface++)
                      {
                        const FaceIndex& face = faces[iface];
                        Real pflux = rflux(face, CPRES);
                        Real areaFrac = m_ebisBox.areaFrac(face);
                        pupdate += isign*areaFrac*pflux;
                      }
                  }
                rupdate += pupdate/m_dx[0];
              }
          }
          //zdiffs look like cartesian
          Real zupdate = 0.;
          {
            int zdir = 1;
            const BaseIFFAB<Real>& zflux = a_centroidFlux[zdir];
            Real faceupdate = 0.0;
            for (SideIterator sit; sit.ok(); ++sit)
              {
                int isign = sign(sit());
                Vector<FaceIndex> faces = m_ebisBox.getFaces(vof, zdir, sit());
                for (int iface = 0; iface < faces.size(); iface++)
                  {
                    const FaceIndex& face = faces[iface];
                    Real areaFrac = m_ebisBox.areaFrac(face);
                    Real faceFlux = zflux(face, ivar);
                    int faceInd =  face.gridIndex(Side::Hi)[rdir];
                    RealVect faceCent = m_ebisBox.centroid(face);
                    Real faceRad    = (Real(faceInd) + faceCent[0]);
                    Real faceFactor = faceRad*areaFrac;

                    faceupdate += isign*vofFactor*faceFactor*faceFlux;
                  }
              }
            zupdate += faceupdate/m_dx[0];
            if (ivar == CMOMY)
              {

                Real pupdate = 0.0;
                for (SideIterator sit; sit.ok(); ++sit)
                  {
                    int isign = sign(sit());
                    Vector<FaceIndex> faces = m_ebisBox.getFaces(vof, zdir, sit());
                    for (int iface = 0; iface < faces.size(); iface++)
                      {
                        const FaceIndex& face = faces[iface];
                        Real areaFrac = m_ebisBox.areaFrac(face);
                        Real pflux = zflux(face, CPRES);
                        pupdate += isign*areaFrac*pflux;
                      }
                  }
                zupdate += pupdate/m_dx[0];
              }
          }

          //add in the stinking boundary flux.
          //only pressure terms here so always diff as if cartesian.
          //the rmomentum term in the ebirregflux fab is left in the rmom
          //slot (so we dont have to futz with the CPRES crap)
          //since the pressure term is all there is.
          Real pflux = a_ebIrregFlux(vof, ivar);
          Real bupdate =  pflux*bndryArea/m_dx[0];

          a_divF(vof, ivar) = rupdate + zupdate + bupdate;
        } //end loop over variables
    } //end loop over vofs
}
/*****************************/

/*****************************/
void
EBPatchPolytropic::getCoveredValuesCons(Vector<Real>& a_covValues)
{
  a_covValues.resize(numConserved());
  a_covValues[CRHO] = 1.0;
  a_covValues[CENG] = 1.0;
  a_covValues[CMOMX] = 0.0;
  a_covValues[CMOMY] = 0.0;
#if CH_SPACEDIM==3
  a_covValues[CMOMZ] = 0.0;
#endif
}
/*****************************/
void
EBPatchPolytropic::getCoveredValuesPrim(Vector<Real>& a_covValues)
{
  a_covValues.resize(numPrimitives());
  a_covValues[QRHO] = 1.0;
  a_covValues[QVELX] = 0.0;
  a_covValues[QVELY] = 0.0;
  a_covValues[QPRES] = 1.0;
  a_covValues[QENTR] = 0.0;
  a_covValues[QINTERN] = 0.0;
  a_covValues[QCVTEMP] = 0.0;
  a_covValues[QC] = 1.0;
#if CH_SPACEDIM==3
  a_covValues[QVELZ] = 0.0;
#endif
}
/*****************************/
void
EBPatchPolytropic::
setGamma(const Real& a_gamma)
{
  m_gamma = a_gamma;
  m_isGammaSet = true;
}
/*****************************/
Real
EBPatchPolytropic::
getGamma() const
{
  CH_assert(m_isGammaSet);
  return m_gamma;
}
/******/
void
EBPatchPolytropic::
setSpecHeat(const Real& a_specHeat)
{
  m_specHeat = a_specHeat;
  m_isSpecHeatSet = true;
}
/*****************************/
Real
EBPatchPolytropic::
getSpecHeat() const
{
  CH_assert(m_isSpecHeatSet);
  return m_specHeat;
}
/******/
bool
EBPatchPolytropic::
usesFlattening() const
{
  CH_assert(m_isSlopeSet);

  return m_useFlattening;
}
/******/
bool
EBPatchPolytropic::
usesArtificialViscosity() const
{
  return m_useArtificialVisc;
}
/******/
bool
EBPatchPolytropic::
usesFourthOrderSlopes() const
{
  CH_assert(m_isSlopeSet);
  return m_useFourthOrderSlopes;
}
/******/
Vector<string>
EBPatchPolytropic::
stateNames()
{
  Vector<string> retval;

  retval.push_back("mass-density");
  if (m_doRZCoords)
    {
      retval.push_back("r-momentum");
      if (SpaceDim >= 2)
        {
          retval.push_back("z-momentum");
        }
      if (SpaceDim >= 3)
        {
          MayDay::Error();
        }
    }
  else
    {
      retval.push_back("x-momentum");
      if (SpaceDim >= 2)
        {
          retval.push_back("y-momentum");
        }

      if (SpaceDim >= 3)
        {
          retval.push_back("z-momentum");
        }
    }
  retval.push_back("energy-density");
  return retval;
}
/******/
Vector<string>
EBPatchPolytropic::
primNames()
{
  Vector<string> retval;

  ParmParse pp;
  int logflag = 0;
  if (pp.contains("logflag"))
    {
      pp.get("logflag", logflag);
    }
  if (logflag == 1)
    {
      retval.push_back("log10density");
    }
  else
    {
      retval.push_back("density");
    }
  if (m_doRZCoords)
    {
      retval.push_back("r-velocity");
      retval.push_back("z-velocity");
    }
  else
    {
      retval.push_back("x-velocity");
      retval.push_back("y-velocity");
    }

#if CH_SPACEDIM==3
  retval.push_back("z-velocity");
#endif

  if (logflag == 1)
    {
      retval.push_back("log10pressure");
      retval.push_back("log10entropy");
    }
  else
    {
      retval.push_back("pressure");
      retval.push_back("entropy");
    }
  retval.push_back("internal_energy");
  retval.push_back("cv_temperature");
  retval.push_back("soundspeed");

#ifdef MODIANO_PROBLEM
  retval.push_back("modiano-velocity-axial");
  retval.push_back("modiano-velocity-tangent0");

#if CH_SPACEDIM==3
  retval.push_back("modiano-velocity-tangent1");
#endif
#endif

  return retval;
}

/******/
Vector<string>
EBPatchPolytropic::
primNamesNoLog()
{
  Vector<string> retval;

  retval.push_back("density");
  if (m_doRZCoords)
    {
      retval.push_back("r-velocity");
      retval.push_back("z-velocity");
    }
  else
    {
      retval.push_back("x-velocity");
      retval.push_back("y-velocity");
    }

#if CH_SPACEDIM==3
  retval.push_back("z-velocity");
#endif

  retval.push_back("pressure");
  retval.push_back("entropy");
  retval.push_back("internal_energy");
  retval.push_back("cv_temperature");
  retval.push_back("soundspeed");

#ifdef MODIANO_PROBLEM
  retval.push_back("modiano-velocity-axial");
  retval.push_back("modiano-velocity-tangent0");

#if CH_SPACEDIM==3
  retval.push_back("modiano-velocity-tangent1");
#endif
#endif

  return retval;
}
/******/
int
EBPatchPolytropic::
numPrimitives() const
{
  return QNUM;
}
/******/
int
EBPatchPolytropic::
numFluxes() const
{
  return FNUM;
}
/******/
int
EBPatchPolytropic::
numConserved() const
{
  return CNUM;
}
/******/
int
EBPatchPolytropic::
numSlopes() const
{
  return QSLOPE;
}
/******/
Interval
EBPatchPolytropic::
velocityInterval() const
{
#if CH_SPACEDIM==2
  Interval retval(QVELX, QVELY);
#elif CH_SPACEDIM==3
  Interval retval(QVELX, QVELZ);
#else
  bogus_spacedim();
#endif
  return retval;
}
/******/
int
EBPatchPolytropic::
pressureIndex() const
{
  return QPRES;
}
/******/
int
EBPatchPolytropic::
densityIndex() const
{
  return QRHO;
}
/******/
int
EBPatchPolytropic::
bulkModulusIndex() const
{
  //phil said this was OK 1-2-2002
  //  MayDay::Warning("returning pressure for bulk modulus");
  return QPRES;
}
/******/
Real
EBPatchPolytropic::
artificialViscosityCoefficient() const
{
  ParmParse pp;
  Real retval;
  pp.get("artificial_viscosity", retval);
  return retval;
}
/******/
Real
EBPatchPolytropic::
getMaxWaveSpeed(const EBCellFAB& a_consState,
                const Box& a_box)
{
  CH_TIME("EBPatchPolytropic::getMaxWaveSpeeed");
  CH_assert(m_isDefined && m_isBoxSet && m_isGammaSet);
  Real speed = 0.0;
  const EBISBox& ebisBox = a_consState.getEBISBox();
  IntVectSet ivs(a_box);
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  //cannot just send to fortran because covered cell values
  //would get in and multiply-valued cells would use bogus
  //values that could also get in.  The could be avoided by
  //either masks or some clever manipulation of the underlying
  //basefab
  for (VoFIterator vofit(ivs, ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      //      const  IntVect& iv = vof.gridIndex();
      Real dense = Max(smallr, a_consState(vof, CRHO));
      Real eng   = Max(small,  a_consState(vof, CENG));
      Real vmax = 0.0;
      Real kinetic = 0.0;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real momen = a_consState(vof, CMOMX + idir);
          Real vel = Abs(momen/dense);
          kinetic += 0.5*vel*vel;
          vmax = Max(vmax, vel);
        }
      Real intEng = eng/dense - kinetic;
      intEng = Max(intEng, (Real)0.001*small);
      Real sqrtarg = m_gamma*(m_gamma -1.0)*intEng;
      CH_assert(sqrtarg >= 0.0);
      Real cspeed = sqrt(sqrtarg);
      if (Abs(cspeed + vmax) > speed)
        {
          speed = cspeed+ vmax;
          if (speed > EBPatchGodunov::getMaxWaveSpeed())
            {
              EBPatchGodunov::setMaxWaveSpeed(speed);
              EBPatchGodunov::setMaxWaveSpeedIV(vof.gridIndex());
            }
        }
    }
  return speed;

}
/******/
void
EBPatchPolytropic::
consToPrim(EBCellFAB&       a_primState,
           const EBCellFAB& a_consState,
           const Box&       a_box,
           int              a_logflag,
           bool             a_verbose)
{
  CH_TIME("EBPatchPolytropic::consToPrim");
  CH_assert(m_isDefined && m_isBoxSet);
  CH_assert(m_isGammaSet);
  CH_assert(a_consState.getRegion().contains(a_box));
  //have to set the covered cell vals so need the cast
  //does not change real data
  //  EBCellFAB& conData = (EBCellFAB&) a_consState;
  //  setCoveredConsVals(conData);
  //debug
  if (a_verbose)
    {
      pout()  << "constoprim " << endl;
    }
  //end debug
  const BaseFab<Real>&   regCons = a_consState.getSingleValuedFAB();
  BaseFab<Real>&   regPrim = a_primState.getSingleValuedFAB();

  int iverbose = 0;
  if (a_verbose) iverbose = 1;
  FORT_CONS2PRM(CHF_BOX(a_box),
                CHF_CONST_FRA(regCons),
                CHF_FRA(regPrim),
                CHF_CONST_INT(a_logflag),
                CHF_CONST_INT(iverbose));


  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Vector<Real> conserved(CNUM);
      Vector<Real> primitive(QNUM);
      for (int ivar = 0; ivar < CNUM; ivar++)
        {
          conserved[ivar] = a_consState(vof, ivar);
        }

      FORT_POINTCONS2PRM(CHF_VR(conserved),
                         CHF_VR(primitive),
                         CHF_CONST_INT(a_logflag));

      for (int ivar = 0; ivar < QNUM; ivar++)
        {
          a_primState(vof, ivar)  = primitive[ivar];
        }
      //  setCoveredPrimVals(a_primState);
    }
}
/******/
void
EBPatchPolytropic::
consToPrim(BaseIVFAB<Real>&  a_primState,
           const EBCellFAB&  a_consState,
           const IntVectSet& a_ivs)
{
  CH_TIME("EBPatchPolytropic::consToPrimIrrReg");
  CH_assert(isDefined());
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Vector<Real> conserved(CNUM);
      Vector<Real> primitive(QNUM);
      for (int ivar = 0; ivar < CNUM; ivar++)
        {
          conserved[ivar] = a_consState(vof, ivar);
        }

      int logflag = 0;
      FORT_POINTCONS2PRM(CHF_VR(conserved),
                         CHF_VR(primitive),
                         CHF_CONST_INT(logflag));

      for (int ivar = 0; ivar < QNUM; ivar++)
        {
          a_primState(vof, ivar)  = primitive[ivar];
        }
    }
}
#ifdef CH_USE_HDF5

void EBPatchPolytropic::expressions(HDF5HeaderData& a_expressions)
{
  char gammaStr[1024], specHeatStr[1024];
  snprintf(gammaStr, 1024, "%g", m_gamma);
  snprintf(specHeatStr, 1024, "%g", m_specHeat);
  a_expressions.m_string["scalar gamma"] = static_cast<const char*>(gammaStr);
  a_expressions.m_string["scalar specHeat"] = static_cast<const char*>(specHeatStr);
  a_expressions.m_string["vector velocity"] = "momentum/<mass-density>";
  a_expressions.m_string["scalar kinetic_energy"] = "dot(velocity,velocity)/2";
  a_expressions.m_string["scalar pressure"] = "(gamma-1)*(<energy-density>-kinetic_energy*<mass-density>";
  a_expressions.m_string["scalar soundspeed"] = "sqrt(gamma*(pressure/<mass-density>))";
  a_expressions.m_string["scalar log10entropy"] = "log10(pressure) - gamma*log10(<mass-density>)";
  a_expressions.m_string["scalar temperature"] = "(<energy-density>-kinetic_energy)/(<mass-density>*specHeat)";
}

#endif

/******/
void
EBPatchPolytropic::
consToPrim(BaseIVFAB<Real>&  a_primState,
           const BaseIVFAB<Real>&  a_consState,
           const IntVectSet& a_ivs)
{
  CH_TIME("EBPatchPolytropic::consToPrimIrr2");
  CH_assert(isDefined());
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      Vector<Real> conserved(CNUM);
      Vector<Real> primitive(QNUM);
      for (int ivar = 0; ivar < CNUM; ivar++)
        {
          conserved[ivar] = a_consState(vof, ivar);
        }
      int logflag = 0;
      FORT_POINTCONS2PRM(CHF_VR(conserved),
                         CHF_VR(primitive),
                         CHF_CONST_INT(logflag));

      for (int ivar = 0; ivar < QNUM; ivar++)
        {
          a_primState(vof, ivar)  = primitive[ivar];
        }
    }
}
/******/
void
EBPatchPolytropic::
primToCons(EBCellFAB&       a_consState,
           const EBCellFAB& a_primState,
           const Box&       a_box)
{
  CH_TIME("EBPatchPolytropic::primToCons");
  CH_assert(isDefined());
  const Box& regBox = a_box;
  CH_assert(a_primState.getRegion().contains(regBox));
  //set covered vals. does not change real data.
  //  EBCellFAB& primData = (EBCellFAB&) a_primState;
  //  setCoveredPrimVals(primData);

  BaseFab<Real>&   regCons = a_consState.getSingleValuedFAB();
  const BaseFab<Real>&  regPrim = a_primState.getSingleValuedFAB();


  FORT_PRM2CONS(CHF_BOX(regBox),
                CHF_FRA(regCons),
                CHF_CONST_FRA(regPrim));


  IntVectSet ivsMulti = m_ebisBox.getMultiCells(regBox);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Vector<Real> conserved(CNUM);
      Vector<Real> primitive(QNUM);
      for (int ivar = 0; ivar < QNUM; ivar++)
        {
          primitive[ivar] = a_primState(vof, ivar);
        }

      FORT_POINTPRM2CONS(CHF_VR(conserved),
                         CHF_VR(primitive));

      for (int ivar = 0; ivar < CNUM; ivar++)
        {
          a_consState(vof, ivar)  = conserved[ivar];
        }
    }
  //  setCoveredConsVals(a_consState);
}
/******/
void
EBPatchPolytropic::
primToCons(BaseIVFAB<Real>&       a_consState,
           const BaseIVFAB<Real>& a_primState,
           const IntVectSet&      a_ivs)
{
  CH_TIME("EBPatchPolytropic::primToConsIrr");
  CH_assert(isDefined());
  for (VoFIterator vofit(a_ivs, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Vector<Real> conserved(CNUM);
      Vector<Real> primitive(QNUM);
      for (int ivar = 0; ivar < QNUM; ivar++)
        {
          primitive[ivar] = a_primState(vof, ivar);
        }

      FORT_POINTPRM2CONS(CHF_VR(conserved),
                         CHF_VR(primitive));

      for (int ivar = 0; ivar < CNUM; ivar++)
        {
          a_consState(vof, ivar)  = conserved[ivar];
        }
    }
}
/*****************************/
void EBPatchPolytropic::
normalPred(EBCellFAB&       a_primLo,
           EBCellFAB&       a_primHi,
           const EBCellFAB& a_primState,
           const EBCellFAB& a_slopePrim,
           const Real&      a_dtbydx,
           const int&       a_dir,
           const Box&       a_box)
{
  CH_TIME("EBPatchPolytropic::normalPred");
  int nslope =  numSlopes();
  CH_assert(isDefined());
  CH_assert(a_primLo.nComp()    >= nslope);
  CH_assert(a_primHi.nComp()    >= nslope);
  CH_assert(a_primState.nComp() >= nslope);
  CH_assert(a_slopePrim.nComp() == nslope);

  Real dtbydx = a_dtbydx;
  const BaseFab<Real>& regState = a_primState.getSingleValuedFAB();
  const BaseFab<Real>& regSlope = a_slopePrim.getSingleValuedFAB();
  BaseFab<Real>& regPrimLo = a_primLo.getSingleValuedFAB();
  BaseFab<Real>& regPrimHi = a_primHi.getSingleValuedFAB();
  int useflat = 0;
  if (usesFlattening())
    useflat = 1;


  FORT_PRED(CHF_BOX(a_box),
            CHF_CONST_FRA(regState),
            CHF_CONST_FRA(regSlope),
            CHF_FRA(regPrimLo),
            CHF_FRA(regPrimHi),
            CHF_CONST_INT(a_dir),
            CHF_CONST_REAL(dtbydx),
            CHF_CONST_INT(useflat));

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(ivsMulti, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      Vector<Real> primlo(QNUM), primhi(QNUM), primit(QNUM), pslope(QNUM);
      for (int ivar = 0; ivar < QNUM; ivar++)
        {
          primit[ivar] = a_primState(vof, ivar);
          pslope[ivar] = a_slopePrim(vof, ivar);
        }


      FORT_POINTPRED(CHF_VR(primit),
                     CHF_VR(pslope),
                     CHF_VR(primlo),
                     CHF_VR(primhi),
                     CHF_CONST_INT(a_dir),
                     CHF_CONST_REAL(dtbydx),
                     CHF_CONST_INT(useflat));


      for (int ivar = 0; ivar < QNUM; ivar++)
        {
          a_primLo(vof, ivar) = primlo[ivar];
          a_primHi(vof, ivar) = primhi[ivar];
        }
    }
}
/******/
void
EBPatchPolytropic::
riemann(EBFaceFAB&       a_flux,
        const EBCellFAB& a_primLeft,
        const EBCellFAB& a_primRigh,
        const int&       a_dir,
        const Box&       a_box)
{
  CH_TIME("EBPatchPolytropic::riemann");
  //int nPrim = numPrimitives();
  Box cellBox = enclosedCells(a_box);
  EBFaceFAB& primGdnv = m_primGdnv[a_dir];

  //explicitly cast the left and right states to modifiable references
  //because we need to shift them to faces and then shift them back.
  BaseFab<Real>& regPrimRigh = (BaseFab<Real>&)a_primRigh.getSingleValuedFAB();
  BaseFab<Real>& regPrimLeft = (BaseFab<Real>&)a_primLeft.getSingleValuedFAB();
  BaseFab<Real>& regPrimGdnv = primGdnv.getSingleValuedFAB();
  //shift the cell centered stuff to edges
  regPrimRigh.shiftHalf(a_dir, -1);
  regPrimLeft.shiftHalf(a_dir,  1);

  CH_assert(regPrimRigh.box().contains(a_box));
  CH_assert(regPrimLeft.box().contains(a_box));
  CH_assert(regPrimGdnv.box().contains(a_box));

  //find the regular part of qgdnv.

  FORT_RIEMANN(CHF_BOX(a_box),
               CHF_CONST_FRA(regPrimLeft),
               CHF_CONST_FRA(regPrimRigh),
               CHF_FRA(regPrimGdnv),
               CHF_CONST_INT(a_dir));


  //shift the cell centered stuff back to cells
  regPrimRigh.shiftHalf(a_dir,  1);
  regPrimLeft.shiftHalf(a_dir, -1);

  //the box sent into this is face-centered.
  //we need to use the cell-centered one it surrounds.
  //this can be more than multivalued cells if there
  //are neighbors that are multivalued

  Box grownBox = cellBox;
  grownBox.grow(a_dir, 1);
  grownBox &= m_domain;
  IntVectSet ivsMulti = m_ebisBox.getIrregIVS(grownBox);

  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingNoBoundary;
  for (FaceIterator faceit(ivsMulti, m_ebisBox.getEBGraph(), a_dir, stopCrit);
      faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      if (a_box.contains(face.gridIndex(Side::Hi)))
        {
          VolIndex vofl = face.getVoF(Side::Lo);
          VolIndex vofr = face.getVoF(Side::Hi);
          Vector<Real> ql(QNUM), qr(QNUM), qgod(QNUM);
          for (int ivar = 0; ivar < QNUM; ivar++)
            {
              ql[ivar] = a_primLeft(vofl, ivar);
              qr[ivar] = a_primRigh(vofr, ivar);
            }


          FORT_POINTRIEMANN(CHF_VR(ql), CHF_VR(qr), CHF_VR(qgod),
                            CHF_CONST_INT(a_dir));


          for (int ivar = 0; ivar < QNUM; ivar++)
            {
              primGdnv(face, ivar) = qgod[ivar];
            }
        }
    }
  //from the godunov state, get the flux
  //fix the prim gdnv boundary face.
  for (SideIterator sit; sit.ok(); ++sit)
    {
      Box boundBox;
      if (sit() == Side::Lo)
        {
          boundBox = adjCellLo(m_domain, a_dir, 1);
        }
      else
        {
          boundBox = adjCellHi(m_domain, a_dir, 1);
        }
      boundBox.shift(a_dir, -sign(sit()));
      IntVectSet ivsBound(boundBox);
      ivsBound &= m_validBox;
      if (!ivsBound.isEmpty())
        {
          FaceStop::WhichFaces stopCrit = FaceStop::AllBoundaryOnly;
          for (FaceIterator faceit(ivsBound, m_ebisBox.getEBGraph(), a_dir, stopCrit);
              faceit.ok(); ++faceit)
            {
              for (int ivar = 0; ivar < QNUM; ivar++)
                {
                  Real primExtrap;
                  if (sit() == Side::Lo)
                    {
                      primExtrap = a_primRigh(faceit().getVoF(Side::Hi), ivar);
                    }
                  else
                    {
                      primExtrap = a_primLeft(faceit().getVoF(Side::Lo), ivar);
                    }
                  primGdnv(faceit(), ivar) = primExtrap;
                }
            }
        }
    }
  getFlux(a_flux, primGdnv, a_dir, a_box);
}
/*****************************/
void
EBPatchPolytropic::
riemann(BaseIVFAB<Real>&        a_coveredFlux,
        const BaseIVFAB<Real>&  a_exteState,
        const EBCellFAB&        a_primState,
        const Vector<VolIndex>& a_vofset,
        const int&              a_dir,
        const Side::LoHiSide&   a_sd,
        const Box&       a_box)
{
  CH_TIME("EBPatchPolytropic::riemannIrr");
  for (int ivof = 0; ivof < a_vofset.size(); ivof++)
    {
      const VolIndex& vof = a_vofset[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          Vector<Real> ql(QNUM), qr(QNUM), qgod(QNUM), flux(FNUM);
          int iloc = vof.gridIndex()[0];
          Real radius;

          if (a_sd == Side::Hi)
            {
              radius = (iloc+1)*m_dx[0];
              for (int ivar = 0; ivar < QNUM; ivar++)
                {
                  ql[ivar] = a_primState(vof, ivar);
                  qr[ivar] = a_exteState(vof, ivar);
                }
            }
          else
            {
              radius = iloc*m_dx[0];
              for (int ivar = 0; ivar < QNUM; ivar++)
                {
                  ql[ivar] = a_exteState(vof, ivar);
                  qr[ivar] = a_primState(vof, ivar);
                }
            }


          FORT_POINTRIEMANN(CHF_VR(ql), CHF_VR(qr), CHF_VR(qgod),
                            CHF_CONST_INT(a_dir));

          //only takes effect if idir == 0
          if (!m_doRZCoords)
            {
              FORT_POINTGETFLUX(CHF_VR(flux), CHF_VR(qgod),
                                CHF_CONST_INT(a_dir));
            }
          else
            {
              FORT_POINTGETFLUXRZ(CHF_VR(flux), CHF_VR(qgod),
                                  CHF_CONST_INT(a_dir),
                                  CHF_CONST_REAL(radius));
            }

          for (int ivar = 0; ivar < FNUM; ivar++)
            {
              a_coveredFlux(vof, ivar) = flux[ivar];
            }
        }
    }
}
/******/
void
EBPatchPolytropic::
getFlux(EBFaceFAB&       a_flux,
        const EBFaceFAB& a_prim,
        const int&       a_dir,
        const Box&       a_box)
{
  CH_TIME("EBPatchPolytropic::getFlux");
  CH_assert(a_flux.direction() == a_dir);
  CH_assert(a_flux.getRegion().contains(a_box));
  CH_assert(a_prim.getRegion().contains(a_box));
  CH_assert(a_prim.direction() == a_dir);

  const BaseFab<Real>& regPrim = a_prim.getSingleValuedFAB();
  BaseFab<Real>& regFlux = a_flux.getSingleValuedFAB();

  if (!m_doRZCoords)
    {
      FORT_GETFLUX(CHF_BOX(a_box),
                   CHF_CONST_FRA(regPrim),
                   CHF_CONST_INT(a_dir),
                   CHF_FRA(regFlux),
                   CHF_CONST_REAL(m_dx[0]));
    }
  else
    {
      FORT_GETFLUXRZ(CHF_BOX(a_box),
                     CHF_CONST_FRA(regPrim),
                     CHF_CONST_INT(a_dir),
                     CHF_FRA(regFlux),
                     CHF_CONST_REAL(m_dx[0]));
    }

  //box sent in is  face-centered
  //this grow dance fixes the problem of stuff being next to multivalued
  //cells
  Box cellBox = enclosedCells(a_box);
  cellBox.grow(a_dir, 1);
  cellBox &= m_domain;

  IntVectSet ivsMulti = m_ebisBox.getMultiCells(cellBox);
  FaceStop::WhichFaces stopCrit = FaceStop::SurroundingNoBoundary;
  for (FaceIterator faceit(ivsMulti, m_ebisBox.getEBGraph(), a_dir, stopCrit);
      faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      if (a_box.contains(face.gridIndex(Side::Hi)))
        {
          Vector<Real> qgdnv(QNUM), fluxvec(FNUM);
          for (int ivar = 0; ivar < QNUM; ivar++)
            {
              qgdnv[ivar] = a_prim(face, ivar);
            }
          int iloc = face.gridIndex(Side::Hi)[0];
          if (!m_doRZCoords)
            {
              FORT_POINTGETFLUX(CHF_VR(fluxvec),
                                CHF_VR(qgdnv),
                                CHF_CONST_INT(a_dir));
            }
          else
            {
              Real radius = iloc*m_dx[0];
              FORT_POINTGETFLUXRZ(CHF_VR(fluxvec),
                                  CHF_VR(qgdnv),
                                  CHF_CONST_INT(a_dir),
                                  CHF_CONST_REAL(radius));
            }


          for (int ivar = 0; ivar < FNUM; ivar++)
            {
              a_flux(face, ivar) = fluxvec[ivar];
            }
        }
    }
}
/******/
void
EBPatchPolytropic::
updateCons(EBCellFAB&            a_consState,
           const EBFaceFAB&        a_flux,
           const BaseIVFAB<Real>&  a_coveredFluxMinu,
           const BaseIVFAB<Real>&  a_coveredFluxPlus,
           const Vector<VolIndex>& a_coveredFaceMinu,
           const Vector<VolIndex>& a_coveredFacePlus,
           const int&              a_dir,
           const Box&              a_box,
           const Real&             a_scale)
{
  CH_TIME("EBPatchPolytropic::updateCons");
  if (m_doRZCoords)
    {
      updateConsRZ( a_consState,
                    a_flux,
                    a_coveredFluxMinu,
                    a_coveredFluxPlus,
                    a_coveredFaceMinu,
                    a_coveredFacePlus,
                    a_dir,
                    a_box,
                    a_scale);
    }
  else
    {
      EBPatchGodunov::updateCons( a_consState,
                                  a_flux,
                                  a_coveredFluxMinu,
                                  a_coveredFluxPlus,
                                  a_coveredFaceMinu,
                                  a_coveredFacePlus,
                                  a_dir,
                                  a_box,
                                  a_scale);
    }

}
void
EBPatchPolytropic::
updateConsRZ(EBCellFAB&              a_consState,
             const EBFaceFAB&        a_flux,
             const BaseIVFAB<Real>&  a_coveredFluxMinu,
             const BaseIVFAB<Real>&  a_coveredFluxPlus,
             const Vector<VolIndex>& a_coveredFaceMinu,
             const Vector<VolIndex>& a_coveredFacePlus,
             const int&              a_dir,
             const Box&              a_box,
             const Real&             a_scale)
{
  CH_assert(m_doRZCoords);
  CH_assert(isDefined());
  CH_assert(a_consState.nComp() <= a_flux.nComp());
  CH_assert(a_consState.nComp() == numConserved() );
  CH_assert((a_dir >= 0) && (a_dir < SpaceDim));
  CH_assert(a_flux.direction() == a_dir);
  CH_assert(a_flux.getRegion().contains(surroundingNodes(a_box, a_dir)));
  CH_assert(a_consState.getRegion().contains(a_box));

  //update the regular vofs.
  //first we have to copy the original state
  //so that the irregular stuff can be updated later.
  EBCellFAB inputConsState(m_ebisBox,
                           a_consState.getRegion(),
                           a_consState.nComp());
  //  setCoveredConsVals(a_consState);
  Interval wholeInterv(0, a_consState.nComp()-1);
  inputConsState.copy(a_consState.getRegion(),  wholeInterv,
                      a_consState.getRegion(),  a_consState, wholeInterv);

  int ncons = numConserved();
  BaseFab<Real>&       consReg = a_consState.getSingleValuedFAB();
  const BaseFab<Real>& fluxReg = a_flux.getSingleValuedFAB();


  FORT_UPDATERZ( CHF_BOX(a_box),
                 CHF_FRA(consReg),
                 CHF_CONST_FRA(fluxReg),
                 CHF_CONST_INT(a_dir),
                 CHF_CONST_INT(ncons),
                 CHF_CONST_REAL(a_scale));

  //update the irregular vofs
  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          for (int ivar = 0; ivar < a_consState.nComp(); ivar++)
            {
              Real update = 0.0;
              //state changed in regular update.
              //need to overwrite this.
              Real originalState= inputConsState(vof, ivar);
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  int isign = sign(sit());
                  Vector<FaceIndex> faces =
                    m_ebisBox.getFaces(vof, a_dir, sit());
                  //if there are faces, use the flux at them
                  Real aggFlux = 0.0;
                  for (int iface = 0; iface < faces.size(); iface++)
                    {
                      const FaceIndex& face = faces[iface];
                      Real faceFlux = a_flux(face, ivar);
                      if (m_doRZCoords)
                        {
                          //add in pressure term if we have taken it out
                          if (((ivar==CMOMX) && (a_dir==0)) ||
                             ((ivar==CMOMY) && (a_dir==1)))
                            {
                              faceFlux += a_flux(face,CPRES);
                            }
                        }
                      aggFlux += faceFlux;
                    }
                  if (faces.size() > 1)
                    aggFlux /= Real(faces.size());

                  update += isign*aggFlux;
                }
              //dx is already divided out (part of scale)
              update *= -a_scale;
              a_consState(vof, ivar) = originalState + update;
            }
        }
    }
  //add in covered face fluxes
  for (SideIterator sit; sit.ok(); ++sit)
    {
      int isign = sign(sit());
      const BaseIVFAB<Real>*  coveredFluxPtr;
      const Vector<VolIndex>* coveredFacePtr;
      if (sit() == Side::Lo)
        {
          coveredFluxPtr = &(a_coveredFluxMinu);
          coveredFacePtr = &(a_coveredFaceMinu);
        }
      else
        {
          coveredFluxPtr = &(a_coveredFluxPlus);
          coveredFacePtr = &(a_coveredFacePlus);
        }
      const BaseIVFAB<Real>&  coveredFlux = *coveredFluxPtr;
      const Vector<VolIndex>& coveredFace = *coveredFacePtr;
      for (int ivof = 0; ivof < coveredFace.size(); ivof++)
        {
          const VolIndex& vof = coveredFace[ivof];
          //          const IntVect& iv = vof.gridIndex();
          if (a_box.contains(vof.gridIndex()))
            {
              for (int ivar = 0; ivar < a_consState.nComp(); ivar++)
                {
                  Real covFlux = coveredFlux(vof, ivar);
                  if (m_doRZCoords)
                    {
                      //add in pressure term if we have taken it out
                      if (((ivar==CMOMX) && (a_dir==0)) ||
                         ((ivar==CMOMY) && (a_dir==1)))
                        {
                          covFlux += coveredFlux(vof, CPRES);
                        }
                    }
                  Real update =   isign*covFlux;
                  //dx is already divided out (part of scale)
                  update *= -a_scale;
                  Real state = a_consState(vof, ivar);
                  a_consState(vof, ivar) = state + update;
                }
            }
        }
    }

  floorConserved(a_consState, a_box);
}
/******/
void
EBPatchPolytropic::floorPrimitives(EBCellFAB&  a_primState,
                                   const Box&  a_box)
{
  CH_TIME("EBPatchPolytropic::floorPrimitives");
  BaseFab<Real>&       primReg = a_primState.getSingleValuedFAB();

  FORT_FLOORPRIM( CHF_BOX(a_box),
                  CHF_FRA(primReg));
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  for (int ivof = 0; ivof<m_irregVoFs.size(); ivof++)
    {
      const VolIndex& vof = m_irregVoFs[ivof];
      if (a_box.contains(vof.gridIndex()))
        {
          a_primState(vof, QRHO)  = Max(a_primState(vof, QRHO) ,  smallr);
          a_primState(vof, QPRES) = Max(a_primState(vof, QPRES),  smallp);
        }
    }
}
/******/
void
EBPatchPolytropic::floorConserved(EBCellFAB&  a_consState,
                                  const Box&  a_box)
{
  CH_TIME("EBPatchPolytropic::floorConserved");
  BaseFab<Real>&       consReg = a_consState.getSingleValuedFAB();

  FORT_FLOORCONS( CHF_BOX(a_box),
                  CHF_FRA(consReg));

  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  IntVectSet ivsIrreg = m_ebisBox.getIrregIVS(a_box);
  for (VoFIterator vofit(ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      a_consState(vof, CRHO) = Max(a_consState(vof, CRHO), smallr);
      a_consState(vof, CENG) = Max(a_consState(vof, CENG), small );
      //      Real rmom =       a_consState(vof, CMOMX) ;
      //      Real zmom =       a_consState(vof, CMOMY) ;
      //      Real ratio = Abs(rmom/zmom);
      //      if ((ratio > 1.0) && (Abs(zmom) > 1.0e-3))
      //        {
      //          cout << vof << " has wacked ratio =" << ratio << endl;
      //        }
    }
}
/******/
void
EBPatchPolytropic::floorConserved(BaseIVFAB<Real>& a_consState,
                                  const IntVectSet& a_ivsIrreg)
{
  CH_TIME("EBPatchPolytropic::floorConservedIrr");
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  for (VoFIterator vofit(a_ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      a_consState(vof, CRHO) = Max(a_consState(vof, CRHO), smallr);
      a_consState(vof, CENG) = Max(a_consState(vof, CENG), small );
    }

}
/******/
void
EBPatchPolytropic::floorPrimitives(BaseIVFAB<Real>& a_primState,
                                   const IntVectSet& a_ivsIrreg)
{
  CH_TIME("EBPatchPolytropic::floorPrimitivesIrr");
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  for (VoFIterator vofit(a_ivsIrreg, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      //floors
      a_primState(vof, QRHO)  = Max(a_primState(vof, QRHO) ,  smallr);
      a_primState(vof, QPRES) = Max(a_primState(vof, QPRES),  smallp);
    }
}
/*****************************/
EBPatchPolytropic::
EBPatchPolytropic():EBPatchGodunov()
{
  m_isGammaSet    = false;
  m_isSpecHeatSet = false;
  m_isArtViscSet  = false;
  m_doRZCoords    = false;
}
/******/
bool
EBPatchPolytropic::
isDefined() const
{
  bool retval = EBPatchGodunov::isDefined() &&  m_isGammaSet;
  return retval;
}

/*****************************/
EBPatchPolytropic::
~EBPatchPolytropic()
{
}
/*****************************/
/*****************************/
void
EBPatchPolytropic::
computeEBIrregFlux(BaseIVFAB<Real>&  a_ebIrregFlux,
                   const EBCellFAB&  a_primState,
                   const EBCellFAB   a_slopePrim[SpaceDim],
                   const IntVectSet& a_irregIVS,
                   const EBCellFAB&  a_source)
{
  CH_TIME("EBPatchPolytropic::computeEBIrregFlux");
  EBCellFAB* primPtr = (EBCellFAB*)(&a_primState);
  EBCellFAB primT;
  if (a_source.isDefined() && (s_conservativeSource))
    {
      CH_assert(a_source.box().contains(m_validBox));
      primT.define(m_ebisBox, a_source.box(), QNUM);
      primT.copy(a_primState);
      EBCellFAB cons(m_ebisBox, a_source.box(), CNUM);
      primToCons(cons, primT, m_validBox);
      cons.plus(a_source, 0.5*m_dt);
      consToPrim(primT, cons, m_validBox, 0);

      primPtr = &primT;
    }
  EBCellFAB& prim = *primPtr;
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  int numCons = numConserved();
  int numPrim = numPrimitives();
  CH_assert(isDefined());
  CH_assert(a_ebIrregFlux.nComp() == numCons);
  CH_assert(prim.nComp() == numPrim);
  CH_assert(a_ebIrregFlux.getIVS().contains(a_irregIVS));

  for (VoFIterator vofit(a_irregIVS, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect centroid = m_ebisBox.bndryCentroid(vof);

      CH_assert(prim(vof, QRHO) > 0.0);

      Real     dense  = prim(vof, QRHO);
      Real     press  = prim(vof, QPRES);
      Real     sound  = sqrt(m_gamma*press/dense);
      RealVect veloc;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          veloc[idir] = prim(vof, QVELX+idir);
        }
      Real     deltaDense = 0.0;
      Real     deltaPress = 0.0;
      RealVect deltaVeloc = RealVect::Zero;
      Real dtOver2H = 0.5*m_dt/m_dx[0];
      Real    rhoc2 = dense*sound*sound;
      for (int extrapDir = 0; extrapDir < SpaceDim; extrapDir++)
        {
          Tuple<int, CH_SPACEDIM-1> tanDirs = PolyGeom::computeTanDirs(extrapDir);
          Real centDir = centroid[extrapDir];
          const EBCellFAB& slopePrimDir =  a_slopePrim[extrapDir];

          Real     slopeDense = slopePrimDir(vof, QRHO);
          Real     slopePress = slopePrimDir(vof, QPRES);
          RealVect slopeVeloc;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              slopeVeloc[idir]= slopePrimDir(vof, QVELX + idir);
            }
          Real slopeUNorm = slopeVeloc[extrapDir];
          Real      velNo =      veloc[extrapDir];

          deltaDense += centDir*slopeDense - dtOver2H*(velNo*slopeDense + dense*slopeUNorm);
          deltaPress += centDir*slopePress - dtOver2H*(velNo*slopePress + rhoc2*slopeUNorm);

          deltaVeloc[extrapDir] += centDir*slopeUNorm - dtOver2H*(velNo*slopeUNorm + slopePress/dense);


          for (int itan = 0; itan < SpaceDim-1; itan++)
            {
              int tanDir = tanDirs[itan];
              Real slopeUTang = slopeVeloc[tanDir];

              deltaVeloc[tanDir] += centDir*slopeUTang - dtOver2H*(velNo*slopeUTang);
            }
        }

      Real xstate[QNUM];
      xstate[QRHO]  = dense + deltaDense;
      xstate[QPRES] = press + deltaPress;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          xstate[QVELX + idir] = veloc[idir] + deltaVeloc[idir];
        }

      if (a_source.isDefined() && (!s_conservativeSource))
        {
          for (int ivar = 0; ivar < QNUM; ivar++)
            {
              Real sourceval = a_source(vof, ivar);
              xstate[ivar] += 0.5*m_dt*sourceval;
            }
        }


      //enforce positivity
      Real pext = xstate[QPRES];
      Real rext = xstate[QRHO];
      if ((pext < 0.0) || (rext < 0.0))
        {
          for (int ivar = 0; ivar < numPrim; ivar++)
            {
              xstate[ivar] = prim(vof, ivar);
            }
        }

      //now compute the boundary condition assuming
      //solid walls at the embedded boundary
      RealVect velocity;
      Real     density, pressure;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          velocity[idir] = xstate[QVELX + idir];
        }

      density  = xstate[QRHO];
      pressure = xstate[QPRES];
      RealVect normal = m_ebisBox.normal(vof);
      Real normalvel = PolyGeom::dot(normal, velocity);
      //HACK putin true normal for computing normal vel
      //RealVect truenorm = getTrueNorm(vof);
      //Real normalvel1 = PolyGeom::dot(truenorm, velocity);
      //normalvel = normalvel1;
      //END HACK
      if (density <= 0.0)
        {
          pout() << "Density is non-positive at " << vofit() << endl;
          // FIXME: Doesn't work in parallel.
          abort();
        }
      CH_assert(density > 0.0);
      pressure = Max(pressure, smallp);
      Real soundspeed = sqrt(m_gamma*pressure/density);
      Real pstar = pressure - density*normalvel*soundspeed;
      if (usesFlattening())
        {
          //limit in case pressure gets too big or too small
          Real ratio = Abs(pstar/pressure);
          if ((ratio > 1.1) || (ratio < 0.9))
            {
              pstar = pressure;
            }
        }

      //modification for artificial viscosity
#if 1
      if (usesArtificialViscosity())
        {
          Real coeff = artificialViscosityCoefficient();
          //approximation to the velocity divergence
          Real divu = 0.0;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              const EBCellFAB& slopeDir = a_slopePrim[idir];
              int ivel = QVELX + idir;
              divu += slopeDir(vof, ivel);
            }
          Real dp = -2.0*density*coeff*normalvel*abs(divu);
          pstar = pstar + dp;
        }
#endif
      //debug
      pstar = Max(pstar, smallp);
      //end debug

      //fluxes as per modiano colella
      a_ebIrregFlux(vof, CRHO) = 0.0;
      a_ebIrregFlux(vof, CENG) = 0.0;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          int imom = CMOMX + idir;
          //negative because our normal points into the fluid
          //these are not separated out for RZ since the pressure flux
          //is the only flux at the EB
          a_ebIrregFlux(vof, imom) = -pstar*normal[idir];
        }
    }
}
/********/
/********/
void
EBPatchPolytropic::
computeBoundaryPress(BaseIVFAB<Real>&  a_pstar,
                     BaseIVFAB<Real>&  a_vnorm,
                     BaseIVFAB<Real>&  a_normDotModiano,
                     BaseIVFAB<Real>&  a_normError,
                     BaseIVFAB<Real>&  a_pstarNoVNorm,
                     BaseIVFAB<Real>&  a_vDotTrueNorm,
                     const RealVect&   a_modianoAxis,
                     const RealVect&   a_modianoCorner,
                     const EBCellFAB&  a_primState,
                     const EBCellFAB&  a_nonConsDivF,
                     const EBCellFAB   a_slopePrim[SpaceDim],
                     const IntVectSet& a_irregIVS,
                     const EBCellFAB&  a_source)
{
  CH_TIME("EBPatchPolytropic::computeBoundaryPress");
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  int numCons = numConserved();
  int numPrim = numPrimitives();
  CH_assert(isDefined());
  CH_assert(a_nonConsDivF.nComp() == numCons);
  CH_assert(a_primState.nComp() == numPrim);
  CH_assert(a_pstar.getIVS().contains(a_irregIVS));

  for (VoFIterator vofit(a_irregIVS, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      RealVect centroid = m_ebisBox.bndryCentroid(vof);

      CH_assert(a_primState(vof, QRHO) > 0.0);

      Real     dense  = a_primState(vof, QRHO);
      Real     press  = a_primState(vof, QPRES);
      Real     sound  = sqrt(m_gamma*press/dense);
      RealVect veloc;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          veloc[idir] = a_primState(vof, QVELX+idir);
        }
      Real     deltaDense = 0.0;
      Real     deltaPress = 0.0;
      RealVect deltaVeloc = RealVect::Zero;
      Real dtOver2H = 0.5*m_dt/m_dx[0];
      Real    rhoc2 = dense*sound*sound;
      for (int extrapDir = 0; extrapDir < SpaceDim; extrapDir++)
        {
          Tuple<int, CH_SPACEDIM-1> tanDirs = PolyGeom::computeTanDirs(extrapDir);
          Real centDir = centroid[extrapDir];
          const EBCellFAB& slopePrimDir =  a_slopePrim[extrapDir];

          Real     slopeDense = slopePrimDir(vof, QRHO);
          Real     slopePress = slopePrimDir(vof, QPRES);
          RealVect slopeVeloc;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              slopeVeloc[idir]= slopePrimDir(vof, QVELX + idir);
            }
          Real slopeUNorm = slopeVeloc[extrapDir];
          Real      velNo =      veloc[extrapDir];

          deltaDense += centDir*slopeDense - dtOver2H*(velNo*slopeDense + dense*slopeUNorm);
          deltaPress += centDir*slopePress - dtOver2H*(velNo*slopePress + rhoc2*slopeUNorm);

          deltaVeloc[extrapDir] += centDir*slopeUNorm - dtOver2H*(velNo*slopeUNorm + slopePress/dense);


          for (int itan = 0; itan < SpaceDim-1; itan++)
            {
              int tanDir = tanDirs[itan];
              Real slopeUTang = slopeVeloc[tanDir];

              deltaVeloc[tanDir] += centDir*slopeUTang - dtOver2H*(velNo*slopeUTang);
            }
        }
      Real xstate[QNUM];
      xstate[QRHO]  = dense + deltaDense;
      xstate[QPRES] = press + deltaPress;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          xstate[QVELX + idir] = veloc[idir] + deltaVeloc[idir];
        }

      if (a_source.isDefined())
        {
          for (int ivar = 0; ivar < QNUM; ivar++)
            {
              Real sourceval = a_source(vof, ivar);
              xstate[ivar] += 0.5*m_dt*sourceval;
            }
        }


      //enforce positivity
      Real pext = xstate[QPRES];
      Real rext = xstate[QRHO];
      if ((pext < 0.0) || (rext < 0.0))
        {
          for (int ivar = 0; ivar < numPrim; ivar++)
            {
              xstate[ivar] = a_primState(vof, ivar);
            }
        }

      //now compute the boundary condition assuming
      //solid walls at the embedded boundary
      RealVect velocity;
      Real     density, pressure;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          velocity[idir] = xstate[QVELX + idir];
        }

      density  = xstate[QRHO];
      pressure = xstate[QPRES];
      RealVect normal = m_ebisBox.normal(vof);
      Real normalvel = PolyGeom::dot(normal, velocity);

      CH_assert(density > 0.0);
      pressure = Max(pressure, smallp);
      Real soundspeed = sqrt(m_gamma*pressure/density);
      Real pstar = pressure - density*normalvel*soundspeed;
      //modification for artificial viscosity
      if (usesArtificialViscosity())
        {
          Real coeff = artificialViscosityCoefficient();
          //approximation to the velocity divergence
          Real divu = 0.0;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              const EBCellFAB& slopeDir = a_slopePrim[idir];
              int ivel = QVELX + idir;
              divu += slopeDir(vof, ivel);
            }
          Real dp = -2.0*coeff*normalvel*Max(-divu, (Real)0.0);
          pstar = pstar + dp;
        }
      a_pstar(vof, 0) = pstar;
      a_vnorm(vof, 0) = normalvel;
      Real normdotmod = PolyGeom::dot(normal, a_modianoAxis);

      a_normDotModiano(vof, 0) = normdotmod;

      a_pstarNoVNorm(vof, 0)  = pressure;

      const IntVect& iv = vof.gridIndex();
      RealVect centroidPt;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          centroidPt[idir] = Real(iv[idir]) + 0.5 + centroid[idir];
        }
      centroidPt *= m_dx[0];
      RealVect trueNorm, closestPt;
      PolyGeom::pointToLine(closestPt, trueNorm, centroidPt,
                            a_modianoCorner, a_modianoAxis);

      Real vdottruenorm = PolyGeom::dot(velocity, trueNorm);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_normError(vof, idir) = normal[idir] - trueNorm[idir];
        }
      a_vDotTrueNorm(vof, 0) = vdottruenorm;

    }
}
/*****************************/
/********/
void
EBPatchPolytropic::
computeBoundaryWDivU(BaseIVFAB<Real>&  a_pstar,
                     BaseIVFAB<Real>&  a_vnorm,
                     BaseIVFAB<Real>&  a_normDotModiano,
                     const RealVect&   a_modianoAxis,
                     const EBCellFAB&  a_primState,
                     const EBCellFAB&  a_nonConsDivF,
                     const EBCellFAB   a_slopePrim[SpaceDim],
                     const IntVectSet& a_irregIVS,
                     const EBCellFAB&  a_source)
{
  CH_TIME("EBPatchPolytropic::computeBoundaryWDivU");
  Real small, smallp, smallu, smallr;
  FORT_GETSMALL(CHF_REAL(small),
                CHF_REAL(smallp),
                CHF_REAL(smallu),
                CHF_REAL(smallr));
  int numCons = numConserved();
  int numPrim = numPrimitives();
  CH_assert(isDefined());
  CH_assert(a_nonConsDivF.nComp() == numCons);
  CH_assert(a_primState.nComp() == numPrim);
  CH_assert(a_pstar.getIVS().contains(a_irregIVS));

  for (VoFIterator vofit(a_irregIVS, m_ebisBox.getEBGraph());
      vofit.ok(); ++vofit)
    {
      //extrapolate in space
      const VolIndex& vof = vofit();
      RealVect centroid = m_ebisBox.bndryCentroid(vof);

      CH_assert(a_primState(vof, QRHO) > 0.0);

      Vector<Real> xstatePrim(QNUM);
      for (int ivar = 0; ivar < numPrim; ivar++)
        {
          xstatePrim[ivar] = a_primState(vof,ivar);
        }
      for (int ivar = 0; ivar < numPrim; ivar++)
        {
          for (int extrapDir = 0; extrapDir < SpaceDim; extrapDir++)
            {
              Real centDir = centroid[extrapDir];
              const EBCellFAB& slopePrimDir =  a_slopePrim[extrapDir];

              xstatePrim[ivar] += centDir*slopePrimDir(vof, ivar);
            }
        }
      if (a_source.isDefined())
        {
          for (int ivar = 0; ivar < QNUM; ivar++)
            {
              Real sourceval = a_source(vof, ivar);
              xstatePrim[ivar] += 0.5*m_dt*sourceval;
            }
        }


      //convert to conserved
      Vector<Real> xstateCons(CNUM);
      FORT_POINTPRM2CONS(CHF_VR(xstateCons),
                         CHF_VR(xstatePrim));


      //extrapolate in time
      for (int ivar = 0; ivar <  CNUM; ivar++)
        {
          Real divf = a_nonConsDivF(vof, ivar);
          xstateCons[ivar] -= 0.5*m_dt*divf;
        }

      //back into primitives
      int logflag = 0;
      FORT_POINTCONS2PRM(CHF_VR(xstateCons),
                         CHF_VR(xstatePrim),
                         CHF_CONST_INT(logflag));

      //enforce positivity
      Real pext = xstatePrim[QPRES];
      Real rext = xstatePrim[QRHO];
      if ((pext < 0.0) || (rext < 0.0))
        {
          for (int ivar = 0; ivar < numPrim; ivar++)
            {
              xstatePrim[ivar] = a_primState(vof, ivar);
            }
        }


      //now compute the boundary condition assuming
      //solid walls at the embedded boundary
      RealVect velocity;
      Real     density, pressure;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          velocity[idir] = xstatePrim[QVELX + idir];
        }

      density  = xstatePrim[QRHO];
      pressure = xstatePrim[QPRES];
      RealVect normal = m_ebisBox.normal(vof);
      Real normalvel = PolyGeom::dot(normal, velocity);
      Real soundspeed = xstatePrim[QC];
      CH_assert(density > 0.0);
      pressure = Max(pressure, smallp);
      Real pstar = pressure - density*normalvel*soundspeed;
      //modification for artificial viscosity
      if (usesArtificialViscosity())
        {
          Real coeff = artificialViscosityCoefficient();
          //approximation to the velocity divergence
          Real divu = 0.0;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              const EBCellFAB& slopeDir = a_slopePrim[idir];
              int ivel = QVELX + idir;
              divu += slopeDir(vof, ivel);
            }
          Real dp = -2.0*coeff*normalvel*Max(-divu, (Real)0.0);
          pstar = pstar + dp;
        }
      a_pstar(vof, 0) = pstar;
      a_vnorm(vof, 0) = normalvel;
      Real normdotmod = PolyGeom::dot(normal, a_modianoAxis);

      a_normDotModiano(vof, 0) = normdotmod;
    }
}
/*****************************/
/*****************************/
void
EBPatchPolytropic::
setCoveredConsVals(EBCellFAB& a_consState)
{
  CH_TIME("EBPatchPolytropic::setCoveredConsVals");
  a_consState.setInvalidData(1.0, CRHO);
  a_consState.setInvalidData(0.0, CMOMX);
  a_consState.setInvalidData(0.0, CMOMY);
#if CH_SPACEDIM==3
  a_consState.setInvalidData(0.0, CMOMZ);
#endif
  a_consState.setInvalidData(1.0, CENG);
}
/**************/
//this is the stable, non-conservative estimate of the solution update
void
EBPatchPolytropic::
nonconservativeDivergence(EBCellFAB&             a_divF,
                          const EBFluxFAB&       a_flux,
                          const BaseIVFAB<Real>  a_coveredFluxMinu[SpaceDim],
                          const BaseIVFAB<Real>  a_coveredFluxPlus[SpaceDim],
                          const Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
                          const Vector<VolIndex> a_coveredFacePlus[SpaceDim],
                          const Box&             a_box)
{
  CH_TIME("EBPatchPolytropic::nonconservativeDiv");
  if (m_doRZCoords)
    {
      nonconservativeDivergenRZ(a_divF, a_flux,
                                a_coveredFluxMinu,
                                a_coveredFluxPlus,
                                a_coveredFaceMinu,
                                a_coveredFacePlus,
                                a_box);
    }
  else
    {
      EBPatchGodunov::nonconservativeDivergence(a_divF, a_flux,
                                                a_coveredFluxMinu,
                                                a_coveredFluxPlus,
                                                a_coveredFaceMinu,
                                                a_coveredFacePlus,
                                                a_box);
    }
}
/*****************************/
void
EBPatchPolytropic::
consUndividedDivergence(BaseIVFAB<Real>&       a_divF,
                        const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                        const BaseIVFAB<Real>& a_ebIrregFlux,
                        const IntVectSet&      a_ivs)
{
  CH_TIME("EBPatchPolytropic::consDiv");
  if (m_doRZCoords)
    {
      consUndividedDivergenRZ(a_divF, a_centroidFlux, a_ebIrregFlux, a_ivs);
    }
  else
    {
      EBPatchGodunov::consUndividedDivergence(a_divF, a_centroidFlux, a_ebIrregFlux, a_ivs);
    }
}
/*****************************/

#include "NamespaceFooter.H"
