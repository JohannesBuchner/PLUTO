#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ModianoIBC.H"
#include "PolyGeom.H"
#include "ModianoF_F.H"
#include "PolyGeom.H"
#include "EBLGIntegrator.H"
#include "VoFIterator.H"
#include "EBPatchPolytropicF_F.H"
#include "EBSolidF_F.H"
#include "EBLoHiCenter.H"
#include "FaceIterator.H"

#include "NamespaceHeader.H"

///
/**
   IBC for doing the modiano problem.
   The lower left corner of the domain corresponds to the bottom corner
   of the channel, so it sort of looks like this:
*/
ModianoIBC::~ModianoIBC()
{
}
ModianoIBC::ModianoIBC(const Real&     a_gamma,
                       const Real&     a_waveAmp,
                       const Real&     a_waveWidth,
                       const RealVect& a_center,
                       const RealVect& a_waveDir,
                       bool            a_freeStreamOnly,
                       bool            a_useNegativeWave)
  :EBPhysIBC()
{
  m_gamma       = a_gamma       ;
  m_waveAmp     = a_waveAmp     ;
  m_waveWidth   = a_waveWidth   ;
  m_waveDir     = a_waveDir     ;

  m_iNegWave    = 0;
  if (a_useNegativeWave)
    m_iNegWave  = 1;

  Real sumSquare;
  PolyGeom::unifyVector(m_waveDir, sumSquare);
  m_center      = a_center      ;
  m_idoFreeStreamOnly = 0;
  if (a_freeStreamOnly)
    m_idoFreeStreamOnly = 1;

  /**/
  FORT_SETMODIANO(CHF_CONST_REAL(m_gamma),
                  CHF_CONST_REAL(m_waveAmp),
                  CHF_CONST_REAL(m_waveWidth),
                  CHF_CONST_REALVECT(m_center),
                  CHF_CONST_REALVECT(m_waveDir),
                  CHF_CONST_INT(m_idoFreeStreamOnly),
                  CHF_CONST_INT(m_iNegWave));
  /**/
  m_freeStreamOnly = a_freeStreamOnly;
  m_isFortranCommonSet = true;
  m_isDefined = false;
}

void
ModianoIBC::define(const ProblemDomain&  a_domain,
                   const RealVect&           a_dx)
{
  m_isDefined = true;;
  m_domain = a_domain;
  m_dx = a_dx[0];
 CH_assert(Abs(a_dx[0]-a_dx[1]) < 1.e-9);
}
/***************/
/***************/
RealVect
mibcGetNormalVector(const RealVect& channelNormal)
{
  RealVect channelSlope;
  //make slope perp to normal
  RealVect channelSortNormal = channelNormal;
  IntVect ivsort, ivsigns;
  PolyGeom::posifyVector(channelSortNormal, ivsigns);
  PolyGeom::sortVector(channelSortNormal, ivsort);
  //divide sorted normal by the largest one
  channelSortNormal /= channelSortNormal[SpaceDim-1];
#if CH_SPACEDIM==2
  if (Abs(channelSortNormal[0]) < PolyGeom::getTolerance())
    {
      //really a one-d channel.  make the slope perp
      //to the normal
      channelSlope = BASISREALV(ivsort[0]);
    }
  else
    {
      //we have a genuine 2D channel.
      channelSlope = RealVect::Unit;
      channelSlope[ivsort[0]] =-1.0/channelSortNormal[0];
      channelSlope[ivsort[0]] *= ivsigns[ivsort[0]];
      channelSlope[ivsort[1]] *= ivsigns[ivsort[1]];
    }
#elif CH_SPACEDIM==3
  if (Abs(channelSortNormal[0]) < PolyGeom::getTolerance() &&
     Abs(channelSortNormal[1]) < PolyGeom::getTolerance())
    {
      //really a one-d channel.  make the slope perp
      //to the normal
      channelSlope = BASISREALV(ivsort[0]);
    }
  else if (Abs(channelSortNormal[0]) < PolyGeom::getTolerance())
    {
      //really a 2-d channel.  make the slope perp
      //to the normal
      //further confuse the issue by dividing out the biggest one
      //this way there is only one non-unity number in the slope
      //channelSortNormal /= channelSortNormal[1];
      channelSlope = RealVect::Unit;
      channelSlope[ivsort[0]] = 0.0;
      channelSlope[ivsort[1]] = -1.0/channelSortNormal[1];
      channelSlope[ivsort[1]] *= ivsigns[ivsort[1]];
      channelSlope[ivsort[2]] *= ivsigns[ivsort[2]];
    }
  else
    {
      //genuine 3d channel.
      channelSlope = RealVect::Unit;
      channelSlope[ivsort[0]] = -0.5/channelSortNormal[0];
      channelSlope[ivsort[0]] *= ivsigns[ivsort[0]];
      channelSlope[ivsort[1]] = -0.5/channelSortNormal[1];
      channelSlope[ivsort[1]] *= ivsigns[ivsort[1]];
    }
#else
  bogus_ch_spacedim_macro();
#endif

  Real dotprod = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      dotprod += channelSlope[idir]*channelNormal[idir];
    }
  if (Abs(dotprod > PolyGeom::getTolerance()))
    {
      MayDay::Error("getnormalvector failed");
    }
  return channelSlope;
}
/***************/
/***************/

void
ModianoIBC::fluxBC(EBFluxFAB&            a_flux,
                   const EBCellFAB&      a_primCenter,
                   const EBCellFAB&      a_primExtrap,
                   const Side::LoHiSide& a_side,
                   const Real&           a_time,
                   const EBISBox&        a_ebisBox,
                   const DataIndex&      a_dit,
                   const Box&            a_box,
                   const Box&            a_faceBox,
                   const int&            a_dir)
{
 CH_assert(m_isDefined);
 CH_assert(m_isFortranCommonSet);
 CH_assert(!m_domain.isPeriodic(a_dir));

  Box FBox = a_flux[a_dir].getSingleValuedFAB().box();
  Box cellBox = FBox;

  // Determine which side and thus shift directions
  int isign = sign(a_side);
  cellBox.shiftHalf(a_dir,isign);

  // Is there a domain boundary next to this grid
  if (!m_domain.contains(cellBox))
    {
      cellBox &= m_domain;
      // Find the strip of cells next to the domain boundary
      Box boundaryBox = bdryBox(cellBox, a_dir, a_side, 1);

      // Shift things to all line up correctly
      boundaryBox.shiftHalf(a_dir,-isign);
      BaseFab<Real>& regFlux = a_flux[a_dir].getSingleValuedFAB();
      const BaseFab<Real>& regPrimExtrap = a_primExtrap.getSingleValuedFAB();
      regFlux.shiftHalf(a_dir,-isign);

      // Set the boundary fluxes
      /**/
      FORT_EXTRAPBC(CHF_FRA(regFlux),
                    CHF_CONST_FRA(regPrimExtrap),
                    CHF_CONST_REAL(m_dx),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(boundaryBox));
      /**/

      // Shift returned fluxes to be face centered
      regFlux.shiftHalf(a_dir,isign);
      //now for the multivalued cells.  Since it is pointwise,
      //the regular calc is correct for all single-valued cells.
      IntVectSet ivs = a_ebisBox.getMultiCells(boundaryBox);
      for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Vector<Real> qgdnv(QNUM);
          Vector<Real> fluxv(FNUM);
          for (int ivar = 0; ivar < QNUM; ivar++)
            {
              qgdnv[ivar] = a_primExtrap(vof, ivar);
            }
          Vector<FaceIndex> bndryFaces =
            a_ebisBox.getFaces(vof, a_dir, a_side);
          for (int iface= 0; iface < bndryFaces.size(); iface++)
            {
              /**/
              FORT_POINTGETFLUX(CHF_VR(fluxv),
                                CHF_VR(qgdnv),
                                CHF_CONST_INT(a_dir));
              /**/
              const FaceIndex& face = bndryFaces[iface];
              for (int ivar = 0; ivar < FNUM; ivar++)
                {
                  a_flux[a_dir](face, ivar) = fluxv[ivar];
                }
            }
        }
    }
}

void
ModianoIBC::initialize(LevelData<EBCellFAB>& a_conState,
  const EBISLayout& a_ebisl) const
{
CH_assert(m_isDefined);
CH_assert(m_isFortranCommonSet);

  // Iterator of all grids in this level
  for (DataIterator dit = a_conState.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = a_ebisl[dit()];
      // Storage for current grid
      EBCellFAB& conFAB = a_conState[dit()];

      initialize(conFAB, ebisBox);
    }
}

void
ModianoIBC::initialize(EBCellFAB& a_conFAB,
                       const EBISBox& a_ebisBox) const
{
 CH_assert(m_isDefined);
 CH_assert(m_isFortranCommonSet);

  BaseFab<Real>& regConFAB = a_conFAB.getSingleValuedFAB();
  // Box of current grid
  Box uBox = regConFAB.box();
  uBox &= m_domain;
  // Set up initial condition in this grid
  /**/
  FORT_MODIANOINIT(CHF_CONST_FRA(regConFAB),
                   CHF_CONST_REAL(m_dx),
                   CHF_BOX(uBox));
  /**/

  //now for the multivalued cells.  Since it is pointwise,
  //the regular calc is correct for all single-valued cells.
  IntVectSet ivs = a_ebisBox.getMultiCells(uBox);
  for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();
      RealVect momentum;
      Real energy, density;
      RealVect xval;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          xval[idir] = m_dx*(Real(iv[idir]) + 0.5);
        }
      /**/
      Real time = 0.0;
      FORT_POINTMODIANOEXACT(CHF_REAL(density),
                             CHF_REALVECT(momentum),
                             CHF_REAL(energy),
                             CHF_CONST_REALVECT(xval),
                             CHF_CONST_REAL(time));
      /**/
      a_conFAB(vof, CRHO) = density;
      a_conFAB(vof, CENG) = energy;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_conFAB(vof, CMOMX+idir) = momentum[idir];
        }
    }//end loop over multivalued cells
}

void
ModianoIBC::computeExactFluxes(EBFluxFAB&       a_exact,
                               BaseIVFAB<Real>  a_coveredExactMinu[SpaceDim],
                               BaseIVFAB<Real>  a_coveredExactPlus[SpaceDim],
                               Vector<VolIndex> a_coveredFacesMinu[SpaceDim],
                               Vector<VolIndex> a_coveredFacesPlus[SpaceDim],
                               IntVectSet       a_coveredSetsMinu[SpaceDim],
                               IntVectSet       a_coveredSetsPlus[SpaceDim],
                               BaseIVFAB<Real>& a_boundaryPress,
                               IntVectSet&      a_irregIVS,
                               const EBISBox&   a_ebisBox,
                               const Box&       a_region,
                               const Box&       a_domainBox,
                               const Real&      a_time)
{
 CH_assert(m_isDefined);
 CH_assert(m_isFortranCommonSet);
  IntVectSet ivsWhole(a_region);
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      //first set the ebfluxfab
      EBFaceFAB& fluxDir = a_exact[faceDir];
      for (FaceIterator faceit(ivsWhole, a_ebisBox.getEBGraph(), faceDir,
                              FaceStop::SurroundingWithBoundary);
          faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          const IntVect& iv = face.gridIndex(Side::Hi);
          RealVect momentum;
          Real energy, density;
          RealVect xval;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (idir != faceDir)
                {
                  xval[idir] = m_dx*(Real(iv[idir]) + 0.5);
                }
              else
                {
                  xval[idir] = m_dx*(Real(iv[idir]));
                }
            }
          /**/
          FORT_POINTMODIANOFLUX(CHF_REAL(density),
                                CHF_REALVECT(momentum),
                                CHF_REAL(energy),
                                CHF_REALVECT(xval),
                                CHF_CONST_INT(faceDir),
                                CHF_CONST_REAL(a_time));
          /**/
          fluxDir(face, CRHO) = density;
          fluxDir(face, CENG) = energy;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fluxDir(face, CMOMX+idir) = momentum[idir];
            }
        }
    }

  int cnum = a_exact.nComp();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_coveredExactPlus[idir].define(a_coveredSetsPlus[idir], a_ebisBox.getEBGraph(), cnum);
      a_coveredExactMinu[idir].define(a_coveredSetsMinu[idir], a_ebisBox.getEBGraph(), cnum);
    }

  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      //now set the covered faces on each side.
      for (int ivof = 0; ivof < a_coveredFacesPlus[faceDir].size(); ivof++)
        {
          const VolIndex& vof = a_coveredFacesPlus[faceDir][ivof];
          IntVect iv = vof.gridIndex();
          RealVect xval;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (idir != faceDir)
                {
                  xval[idir] = m_dx*(Real(iv[idir]) + 0.5);
                }
              else
                {
                  //because this is the face on the plus side of the vof,
                  //need to add one to the grid index
                  xval[idir] = m_dx*(Real(iv[idir]) + 1.0);
                }
            }
          RealVect momentum;
          Real energy, density;
          /**/
          FORT_POINTMODIANOFLUX(CHF_REAL(density),
                                CHF_REALVECT(momentum),
                                CHF_REAL(energy),
                                CHF_REALVECT(xval),
                                CHF_CONST_INT(faceDir),
                                CHF_CONST_REAL(a_time));
          /**/
          a_coveredExactPlus[faceDir](vof, CRHO) = density;
          a_coveredExactPlus[faceDir](vof, CENG) = energy;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              a_coveredExactPlus[faceDir](vof, CMOMX+idir) = momentum[idir];
            }
        }

      //now set the covered faces on each side.
      for (int ivof = 0; ivof < a_coveredFacesMinu[faceDir].size(); ivof++)
        {
          const VolIndex& vof = a_coveredFacesMinu[faceDir][ivof];
          //because this is the face on the minus side of the vof,
          //no need to add one to the grid index
          IntVect iv = vof.gridIndex();
          RealVect momentum;
          Real energy, density;
          RealVect xval;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (idir != faceDir)
                {
                  xval[idir] = m_dx*(Real(iv[idir]) + 0.5);
                }
              else
                {
                  xval[idir] = m_dx*(Real(iv[idir]));
                }
            }
          /**/
          FORT_POINTMODIANOFLUX(CHF_REAL(density),
                                CHF_REALVECT(momentum),
                                CHF_REAL(energy),
                                CHF_REALVECT(xval),
                                CHF_CONST_INT(faceDir),
                                CHF_CONST_REAL(a_time));
          /**/
          a_coveredExactMinu[faceDir](vof, CRHO) = density;
          a_coveredExactMinu[faceDir](vof, CENG) = energy;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              a_coveredExactMinu[faceDir](vof, CMOMX+idir) = momentum[idir];
            }
        }
    }

  //now set the irregular boundary flux
  for (VoFIterator vofit(a_irregIVS, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      //irreg fluxes live at centroids
      IntVect iv = vof.gridIndex();
      RealVect momentum, velocity;
      Real energy, density, pressure;
      RealVect centroid = a_ebisBox.bndryCentroid(vof);
      RealVect xval;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          xval[idir] = m_dx*(Real(iv[idir]) + 0.5 + centroid[idir]);
        }
      FORT_POINTMODIANOEXACT(CHF_REAL(density),
                             CHF_REALVECT(momentum),
                             CHF_REAL(energy),
                             CHF_CONST_REALVECT(xval),
                             CHF_CONST_REAL(a_time));

      Vector<Real> conserved(CNUM), primitive(QNUM);
      conserved[CRHO] = density;
      conserved[CENG] = energy;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          conserved[CMOMX + idir] = momentum[idir];
        }
      int logflag = 0;
      FORT_POINTCONS2PRM(CHF_VR(conserved),
                         CHF_VR(primitive),
                         CHF_CONST_INT(logflag));
      pressure = primitive[QPRES];
      a_boundaryPress(vof, 0) = pressure;
    }
}
void
ModianoIBC::setToExactConsAndPrim(EBCellFAB&     a_consState,
                                  const EBISBox& a_ebisBox,
                                  const Real&    a_time) const
{
  BaseFab<Real>& regConFAB = a_consState.getSingleValuedFAB();
  // Box of current grid
  Box uBox = regConFAB.box();
  uBox &= m_domain;
  // Set up initial condition in this grid
  /**/
  FORT_MODIANOEXACTCONSANDPRIM(CHF_CONST_FRA(regConFAB),
                               CHF_CONST_REAL(m_dx),
                               CHF_CONST_REAL(a_time),
                               CHF_BOX(uBox));
  /**/

  //now for the multivalued cells.  Since it is pointwise,
  //the regular calc is correct for all single-valued cells.
  IntVectSet ivs = a_ebisBox.getMultiCells(uBox);
  for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();
      RealVect momentum, velocity, modianovel;
      Real energy, density, entropy, soundspeed, pressure;
      RealVect xval;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          xval[idir] = m_dx*(Real(iv[idir]) + 0.5);
        }
      /**/
      FORT_POINTMODIANOEXACTCONSANDPRIM(CHF_REAL(density),
                                        CHF_REALVECT(momentum),
                                        CHF_REAL(energy),
                                        CHF_REALVECT(velocity),
                                        CHF_REAL(pressure),
                                        CHF_REAL(entropy),
                                        CHF_REAL(soundspeed),
                                        CHF_REALVECT(modianovel),
                                        CHF_CONST_REALVECT(xval),
                                        CHF_CONST_REAL(a_time));
      /**/
      a_consState(vof, CRHO) = density;
      a_consState(vof, CENG) = energy;
      a_consState(vof, CNUM+QRHO) = density;
      a_consState(vof, CNUM+QPRES) = pressure;
      a_consState(vof, CNUM+QENTR+1)= entropy;
      a_consState(vof, CNUM+QC)= soundspeed;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_consState(vof, CMOMX+idir) = momentum[idir];
          a_consState(vof, CNUM+QVELX+idir) = velocity[idir];
#ifdef MODIANO_PROBLEM
          a_consState(vof, CNUM+QMVAX+idir) = modianovel[idir];
#endif
        }
    }//end loop over multivalued cells
}
void
ModianoIBC::setToExact(EBCellFAB&     a_consState,
                       const EBISBox& a_ebisBox,
                       const Real&    a_time) const
{
  BaseFab<Real>& regConFAB = a_consState.getSingleValuedFAB();
  // Box of current grid
  Box uBox = regConFAB.box();
  uBox &= m_domain;
  // Set up initial condition in this grid
  /**/
  FORT_MODIANOEXACT(CHF_CONST_FRA(regConFAB),
                    CHF_CONST_REAL(m_dx),
                    CHF_CONST_REAL(a_time),
                    CHF_BOX(uBox));
  /**/

  //now for the multivalued cells.  Since it is pointwise,
  //the regular calc is correct for all single-valued cells.
  IntVectSet ivs = a_ebisBox.getMultiCells(uBox);
  for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();
      RealVect momentum;
      Real energy, density;
      RealVect xval;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          xval[idir] = m_dx*(Real(iv[idir]) + 0.5);
        }
      /**/
      FORT_POINTMODIANOEXACT(CHF_REAL(density),
                             CHF_REALVECT(momentum),
                             CHF_REAL(energy),
                             CHF_CONST_REALVECT(xval),
                             CHF_CONST_REAL(a_time));
      /**/
      a_consState(vof, CRHO) = density;
      a_consState(vof, CENG) = energy;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_consState(vof, CMOMX+idir) = momentum[idir];
        }
    }//end loop over multivalued cells
}
void
ModianoIBC::setToExactConsAndPrim(LevelData<EBCellFAB>& a_consState,
                                  const EBISLayout& a_ebisl,
                                  const Real&       a_time) const
{
 CH_assert(m_isDefined);
 CH_assert(m_isFortranCommonSet);

  // Iterator of all grids in this level
  for (DataIterator dit = a_consState.dataIterator(); dit.ok(); ++dit)
    {
      setToExactConsAndPrim(a_consState[dit()], a_ebisl[dit()], a_time);
    } //end loop over boxes
}
void
ModianoIBC::setToExact(LevelData<EBCellFAB>& a_consState,
                       const EBISLayout& a_ebisl,
                       const Real&       a_time) const
{
 CH_assert(m_isDefined);
 CH_assert(m_isFortranCommonSet);

  // Iterator of all grids in this level
  for (DataIterator dit = a_consState.dataIterator(); dit.ok(); ++dit)
    {
      setToExact(a_consState[dit()], a_ebisl[dit()], a_time);
    } //end loop over boxes
}
void
ModianoIBC::setBndrySlopes(EBCellFAB&       a_deltaPrim,
                           const EBCellFAB& a_primState,
                           const EBISBox&   a_ebisBox,
                           const Box&       a_box,
                           const int&       a_dir)
{
  // don't want to deal with the periodic case
 CH_assert(!m_domain.isPeriodic(a_dir));
 CH_assert(m_isDefined);
 CH_assert(m_isFortranCommonSet);

  Box loBox,hiBox,centerBox,domain;
  int hasLo,hasHi;
  Box slopeBox = a_deltaPrim.getRegion()&m_domain;
  int numSlop = a_deltaPrim.nComp();
  // Generate the domain boundary boxes, loBox and hiBox, if there are
  // domain boundarys there
  eblohicenter(loBox,hasLo,hiBox,hasHi,centerBox,domain,
             slopeBox,m_domain,a_dir);

  // Set the boundary slopes if necessary
  if ((hasLo != 0) || (hasHi != 0))
    {
      BaseFab<Real>& regDeltaPrim = a_deltaPrim.getSingleValuedFAB();
      const BaseFab<Real>& regPrimState = a_primState.getSingleValuedFAB();
      /**/
      FORT_SLOPEBCS(CHF_FRA(regDeltaPrim),
                    CHF_CONST_FRA(regPrimState),
                    CHF_CONST_REAL(m_dx),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(loBox),
                    CHF_CONST_INT(hasLo),
                    CHF_BOX(hiBox),
                    CHF_CONST_INT(hasHi));
      /**/
    }
  for (SideIterator sit; sit.ok(); ++sit)
    {
      bool doThisSide;
      Box thisBox;
      if (sit() == Side::Lo)
        {
          doThisSide = (hasLo != 0);
          thisBox = loBox;
        }
      else
        {
          doThisSide = (hasHi != 0);
          thisBox = hiBox;
        }
      if (doThisSide)
        {

          // the cells for which the regular calculation
          //are incorrect are the grown set of the multivalued
          //cells intersected with the appropriate bndry box
          //and added to the irregular cells on the boundary.
          Box boxGrown = thisBox;
          boxGrown.grow(a_dir, 1);
          boxGrown &= m_domain;
          IntVectSet ivs = a_ebisBox.getMultiCells(boxGrown);
          ivs &= loBox;
          IntVectSet ivsIrreg = a_ebisBox.getIrregIVS(thisBox);
          ivs |= ivsIrreg;
          for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              //all slopes at boundary get set to zero
              //except normal velocity.  just following the
              //fortran here
              int inormVelVar = a_dir + QVELX;
              for (int ivar = 0; ivar < numSlop; ivar++)
                {
                  if (ivar != inormVelVar)
                    {
                      a_deltaPrim(vof, ivar) = 0.0;
                    }
                }

              //for normal velocity
              //do strangely limited slope
              //just lifted from the fortran.
              Side::LoHiSide otherSide = flip(sit());
              Vector<FaceIndex> faces =
                a_ebisBox.getFaces(vof, a_dir, otherSide);

              Real thisVel = a_primState(vof, inormVelVar);
              Real otherVel = 0.;
              for (int iface = 0; iface < faces.size(); iface++)
                {
                  VolIndex otherVoF = faces[iface].getVoF(otherSide);
                  otherVel += a_primState(otherVoF,inormVelVar);
                }
              if (faces.size() > 0)
                {
                  otherVel /= Real(faces.size());
                  Real slope;
                  if (sit() == Side::Lo)
                    {
                      slope = otherVel - thisVel;
                    }
                  else
                    {
                      slope = thisVel  - otherVel;
                    }
                  //trick to make sure state will not change sign
                  if (slope*thisVel < 0)
                    {
                      a_deltaPrim(vof, inormVelVar) = 0.0;
                    }
                  else
                    { //slope and vel same sign
                      Real rsign = 1.0;
                      if (thisVel < 0.0)
                        rsign = -1.0;
                      //told you this was odd.
                      Real dvmin = Min(Abs(slope), Abs((Real)2.*thisVel));
                      a_deltaPrim(vof, inormVelVar) = rsign*dvmin;
                    }
                }
              else //no faces on high side of low boundary vof
                {
                  a_deltaPrim(vof, inormVelVar) = 0.0;
                }
            } //end loop over irregular vofs at boundary
        }
    } //end loop over boundary sides
}

#include "NamespaceFooter.H"
