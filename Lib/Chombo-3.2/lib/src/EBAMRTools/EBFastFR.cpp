#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBArith.H"
#include "EBFastFR.H"
#include "LayoutIterator.H"
#include "BaseIVFactory.H"
#include "EBCoarToCoarRedist.H"
#include "EBCoarToFineRedist.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "EBIndexSpace.H"
#include "parstream.H"
#include "EBAlias.H"
#include "EBCellFactory.H"
#include "EBLevelDataOps.H"
#include "CH_Timer.H"
#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>
#include "EBCFData.H"
#include "NamespaceHeader.H"
IntVect   ivdebebffr(D_DECL(5,1,0));
VolIndex vofdebebffr(ivdebebffr, 0);
bool EBFastFR::s_verbose = false;

void
dumpEBFFR(const LevelData<EBCellFAB>& a_data, const string& a_ident)
{
  pout() << a_ident << ":" << endl;
  bool found = false;
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {

      if (a_data[dit()].box().contains(ivdebebffr))
        {
          found = true;
          pout() << ivdebebffr << " found in data grid box " << a_data.disjointBoxLayout()[dit()] << " = " ;
          for (int ivar = 0; ivar < a_data.nComp(); ivar++)
            {

              pout() << a_data[dit()](vofdebebffr, ivar);
              if (ivar < a_data.nComp()-1) pout() << ", ";
            }
          pout() << endl;
        }
    }

  if (!found)
    pout() << ivdebebffr << " not found" << endl;
}
void
EBFastFR::setDefaultValues()
{
  m_levelFluxReg = NULL;
  m_isDefined = false;
  m_nComp = -1;
  m_refRat = -1;
}

class EBAddOpEBFFR : public LDOperator<EBCellFAB>
{
public:
  EBAddOpEBFFR()
  {
  }

  virtual void linearIn(EBCellFAB& arg,  void* buf, const Box& R,
                        const Interval& comps) const
  {
    EBCellFAB incr;
    incr.clone(arg);
    incr.setVal(0.);
    incr.linearIn(buf, R, comps);

    EBCellFAB&       dst = arg;
    const EBCellFAB& src = incr;

    int isrc = comps.begin();
    int idst = comps.begin();
    int inco = comps.size();
    dst.plus(src, isrc, idst, inco);
  }

  void op(EBCellFAB& dst,
          const Box& RegionFrom,
          const Interval& Cdst,
          const Box& RegionTo,
          const EBCellFAB& src,
          const Interval& Csrc) const
  {
    CH_assert(Csrc.size() == Cdst.size());
    int isrc = Csrc.begin();
    int idst = Cdst.begin();
    int inco = Csrc.size();
    dst.plus(src, isrc, idst, inco);
  }
};
/*******************/
int
EBFastFR::
index(int dir, Side::LoHiSide side)
{
  CH_assert(dir >= 0);
  CH_assert(dir < SpaceDim);
  CH_assert((side == Side::Lo) || (side == Side::Hi));

  int ioffset;
  if (side == Side::Lo)
    ioffset = 0;
  else
    ioffset = 1;

  return ioffset*SpaceDim+dir;
}
/*******************/
EBFastFR::
EBFastFR(const EBLevelGrid& a_eblgFine,
         const EBLevelGrid& a_eblgCoar,
         const int&         a_refRat,
         const int&         a_nvar,
         bool a_forceNoEBCF)
{
  setDefaultValues();
  define(a_eblgFine, a_eblgCoar, a_refRat, a_nvar, a_forceNoEBCF);
}
/*******************/
bool
EBFastFR::computeHasEBCF()
{
  CH_TIME("EBFastFR::computeHasEBCF");
  const ProblemDomain&            domai =   m_eblgFine.getDomain();
  const EBISLayout&               ebisl =   m_eblgFine.getEBISL();
  const DisjointBoxLayout&        grids =   m_eblgFine.getDBL();
  const LayoutData<IntVectSet>&   cfivs = *(m_eblgFine.getCFIVS());
  bool doEBCrossing = false;

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
        {
          Box grid  = grids.get(dit());
          IntVect ivgrow = IntVect::Zero;
          IntVectSet ebcfivsLo, ebcfivsHi;
          EBCFData::getEBCFIVSGrid(ebcfivsLo, grid, idir, Side::Lo, ivgrow, domai, cfivs[dit()], ebisl[dit()]);
          EBCFData::getEBCFIVSGrid(ebcfivsHi, grid, idir, Side::Hi, ivgrow, domai, cfivs[dit()], ebisl[dit()]);

          //keep track if this  crossing ever happens
          if ((!ebcfivsLo.isEmpty()) ||
             (!ebcfivsHi.isEmpty()))
            {
              doEBCrossing = true;
            }
        }
    }

  //in the case of parallel, need to check if ANY of the procs
  //have ebcrossing
#ifdef CH_MPI
  int gatherint = 0;
  if (doEBCrossing) gatherint = 1;
  //  pout() << "before gather = " << gatherint << endl;

  int idoebcf;
  MPI_Allreduce(&gatherint, &idoebcf, 1, MPI_INT,
                MPI_MAX, Chombo_MPI::comm);
  doEBCrossing = (idoebcf==1);


#endif

  return doEBCrossing;
}
/*******************/
void
EBFastFR::
define(const EBLevelGrid& a_eblgFine,
       const EBLevelGrid& a_eblgCoar,
       const int&         a_refRat,
       const int&         a_nvar,
       bool a_forceNoEBCF)
{
  CH_TIME("EBFastFR::define");
  clear();
  m_refRat   = a_refRat;
  m_nComp    = a_nvar;
  m_eblgFine = a_eblgFine;
  m_eblgCoar = a_eblgCoar;
  if (!m_eblgFine.coarsenable(a_refRat))
    {
      MayDay::Error("EBFastFR::define: dbl not coarsenable by refrat");
    }

  coarsen(m_eblgCoFi, a_eblgFine, a_refRat);
  m_eblgCoFi.getEBISL().setMaxRefinementRatio(a_refRat, m_eblgCoFi.getEBIS());
  m_reverseCopier.ghostDefine(m_eblgCoFi.getDBL(), m_eblgCoar.getDBL(),m_eblgCoar.getDomain(),IntVect::Unit);

#if (CH_SPACEDIM == 2)
  m_nrefdmo = m_refRat;
#elif (CH_SPACEDIM == 3)
  m_nrefdmo = m_refRat*m_refRat;
#else
  bogus spacedim;
#endif

  m_levelFluxReg = new LevelFluxRegister();

  //  pout() << "before regular define" << endl;
  m_levelFluxReg->define(a_eblgFine.getDBL(),
                         a_eblgCoar.getDBL(),
                         a_eblgFine.getDomain(),
                         a_refRat, a_nvar);
  if (a_forceNoEBCF)
    {
      m_hasEBCF = false;
    }
  else
    {
      m_hasEBCF = computeHasEBCF();
    }

  //if no EBCF--nothing happens here but calls to level flux register
  if (m_hasEBCF)
    {
      defineSetsAndIterators();
      defineBuffers();
    }
  //set all registers to zero to start.
  setToZero();
  m_isDefined = true;
  //  pout() << "leaving define" << endl;
}
/*********/
void
getCFIVSSubset(Vector<IntVectSet>&      a_cfivsCoFi,
               const Box&               a_subreCoar,
               const DisjointBoxLayout& a_layouCoFi,
               const int&               a_idir,
               const Side::LoHiSide&    a_sd)
{

  a_cfivsCoFi.resize(0);
  for (LayoutIterator lito = a_layouCoFi.layoutIterator(); lito.ok(); ++lito)
    {
      const Box& gridCoFiO = a_layouCoFi[lito()];
      //the - 1 is due to the adjacent cell box thing
      if (gridCoFiO.bigEnd(0) < a_subreCoar.smallEnd(0) - 1 )
        {
          //can skip rest cuz we haven't gotten to something interesting
          continue;
        }

      Box cfboxCoFi = adjCellBox(gridCoFiO, a_idir, a_sd, 1);

      //only want the cells within the input subregion
      cfboxCoFi &= a_subreCoar;

      IntVectSet ivsCoFi(cfboxCoFi);
      for (LayoutIterator liti = a_layouCoFi.layoutIterator(); liti.ok(); ++liti)
        {
          const Box& gridCoFiI = a_layouCoFi[liti()];
          if (gridCoFiI.bigEnd(0) < cfboxCoFi.smallEnd(0))
            {
              //can skip rest cuz we haven't gotten to something interesting
              continue;
            }

          //remove from the set boxes that overlap on the own grid
          //(fine-fine interface points)
          ivsCoFi -= gridCoFiI;

          if ((gridCoFiI.smallEnd(0) > cfboxCoFi.bigEnd(0)) || ivsCoFi.isEmpty())
            {
              //can break out of loop, since we know that the smallEnd
              // of all the remaining boxes are lexigraphically beyond this ghosted box.
              break;
            }
        }
      if (!ivsCoFi.isEmpty())
        {
          a_cfivsCoFi.push_back(ivsCoFi);
        }

      //the + 1 is due to the adjacent cell box thing
      if (gridCoFiO.smallEnd(0) > a_subreCoar.bigEnd(0) + 1 )
        {
          //can break out of loop, since we know that the smallEnd
          // of all the remaining boxes are lexigraphically beyond this ghosted box.
          break;
        }

    }
}
/*********/
void
EBFastFR::
defineSetsAndIterators()
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {

      for (SideIterator sit; sit.ok(); ++sit)
        {
          int iindex = index(idir, sit());
          m_setsCoFi[iindex].define(m_eblgCoFi.getDBL());
          m_vofiCoFi[iindex].define(m_eblgCoFi.getDBL());

          for (DataIterator dit = m_eblgCoFi.getDBL().dataIterator(); dit.ok(); ++dit)
            {
              const Box& grid = m_eblgCoFi.getDBL().get(dit());
              Box boxCoFi;
              if (sit() == Side::Lo)
                {
                  boxCoFi = adjCellLo(grid, idir, 1);
                }
              else
                {
                  boxCoFi = adjCellHi(grid, idir, 1);
                }
              const EBGraph& ebgraphCoFi =    m_eblgCoFi.getEBISL()[dit()].getEBGraph();
              const IntVectSet& ivsCF    =  (*m_eblgCoFi.getCFIVS())[dit()];
              IntVectSet ivsEBCF         = ebgraphCoFi.getIrregCells(boxCoFi);

              ivsEBCF &= ivsCF;
              ivsEBCF &= boxCoFi;

              (m_setsCoFi[iindex])[dit()]      = ivsEBCF;
              (m_vofiCoFi[iindex])[dit()].define(ivsEBCF, ebgraphCoFi);
            }

          m_setsCoar[iindex].define(m_eblgCoar.getDBL());
          m_vofiCoar[iindex].define(m_eblgCoar.getDBL());

          for (DataIterator dit = m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
            {
              const Box& boxCoar = m_eblgCoar.getDBL()[dit()];
              IntVectSet irregIVS = m_eblgCoar.getEBISL()[dit()].getIrregIVS(boxCoar);
              Vector<IntVectSet>  cfIVSVec;
              int iindex = index(idir, sit());
              getCFIVSSubset(cfIVSVec, boxCoar, m_eblgCoFi.getDBL(), idir, sit());
              const EBGraph& ebgraphCoar =  m_eblgCoar.getEBISL()[dit()].getEBGraph();
              (m_setsCoar[iindex])[dit()].resize(cfIVSVec.size());
              (m_vofiCoar[iindex])[dit()].resize(cfIVSVec.size());
              for (int ivec = 0; ivec < cfIVSVec.size(); ivec++)
                {
                  IntVectSet ivsEBCF = irregIVS;
                  ivsEBCF &= cfIVSVec[ivec];

                  //need to grow by one to get on both sides of
                  //coarse-fine interface
                  (m_setsCoar[iindex])[dit()][ivec]      = ivsEBCF;
                  (m_vofiCoar[iindex])[dit()][ivec].define(ivsEBCF, ebgraphCoar);
                }
            }
        }
    }
}
/************/
void
EBFastFR::
defineBuffers()
{
  EBCellFactory ebcellfactCoar(m_eblgCoar.getEBISL());
  EBCellFactory ebcellfactCoFi(m_eblgCoFi.getEBISL());
  m_saveCoar.define(m_eblgCoar.getDBL(), m_nComp, IntVect::Unit, ebcellfactCoar);
  m_delUCoar.define(m_eblgCoar.getDBL(), m_nComp, IntVect::Unit, ebcellfactCoar);
  m_delUDiff.define(m_eblgCoar.getDBL(), m_nComp, IntVect::Unit, ebcellfactCoar);
  m_delUCoFi.define(m_eblgCoFi.getDBL(), m_nComp, IntVect::Unit, ebcellfactCoFi);
}

/*********/
EBFastFR::
EBFastFR()
{
  setDefaultValues();
}
/*********/
void
EBFastFR::
clear()
{
  if (m_isDefined)
    {
      m_saveCoar.clear();
      m_delUCoar.clear();
      m_delUDiff.clear();
      m_delUCoFi.clear();
      if (m_levelFluxReg != NULL)
        {
          delete m_levelFluxReg;
          m_levelFluxReg = NULL;
        }
    }
  m_isDefined = false;
}
/*******************/
EBFastFR::~EBFastFR()
{
  clear();
}
/*******************/
bool
EBFastFR::
isDefined() const
{
  return m_isDefined;
}
/*******************/
void
EBFastFR::
setToZero()
{
  CH_TIME("EBFastFR::setToZero");
  m_levelFluxReg->setToZero();

  if (m_hasEBCF)
    {
      irregSetToZero();
    }
}
void
EBFastFR::
irregSetToZero()
{
  CH_TIME("EBFastFR::irregSetToZero");
  EBLevelDataOps::setVal(m_saveCoar, 0.0);
  EBLevelDataOps::setVal(m_delUCoar, 0.0);
  EBLevelDataOps::setVal(m_delUDiff, 0.0);
  EBLevelDataOps::setVal(m_delUCoFi, 0.0);
}
/*******************/
void
EBFastFR::
incrementCoarseBoth(const EBFaceFAB&      a_coarFlux,
                    const Real&           a_scale,
                    const DataIndex&      a_coarDatInd,
                    const Interval&       a_variables,
                    const int&            a_dir,
                    const Side::LoHiSide& a_sd)
{
  incrementCoarRegul(a_coarFlux,
                     a_scale,
                     a_coarDatInd,
                     a_variables,
                     a_dir,
                     a_sd);

  incrementCoarIrreg(a_coarFlux,
                     a_scale,
                     a_coarDatInd,
                     a_variables,
                     a_dir,
                     a_sd);

}
/*************/
void
EBFastFR::
incrementCoarRegul(const EBFaceFAB& a_coarFlux,
                   const Real&      a_scale,
                   const DataIndex& a_coarDatInd,
                   const Interval&  a_variables,
                   const int&       a_dir,
                   const Side::LoHiSide& a_sd)
{
  CH_assert(m_isDefined);
  CH_TIME("EBFastFR::incrementCoarseBoth");
  FArrayBox& coarFluxFAB = (FArrayBox&)(a_coarFlux.getFArrayBox());

  //increment  as if there were no EB
  m_levelFluxReg->incrementCoarse(coarFluxFAB, a_scale, a_coarDatInd,
                                  a_variables, a_variables, a_dir, a_sd);
}
/*******************/
void
EBFastFR::
incrementCoarIrreg(const EBFaceFAB&       a_coarFlux,
                   const Real&            a_scale,
                   const DataIndex&       a_coarDatInd,
                   const Interval&        a_variables,
                   const int&             a_dir,
                   const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFastFR::incrementCoarseIrreg");
  CH_assert(m_isDefined);
  if (m_hasEBCF)
    {

      const EBISBox& ebisBox     = m_eblgCoar.getEBISL()[a_coarDatInd];
      int iindex = index(a_dir, a_sd);

      Vector<VoFIterator>& vofits = (m_vofiCoar[iindex])[a_coarDatInd];

      for (int ivofits = 0; ivofits < vofits.size(); ivofits++)
        {
          VoFIterator& vofit = vofits[ivofits];

          //increment coar buffer with the area-weighted sum
          //of the fluxes on the faces.
          //buf += sum(areaFrac*flux)
          for (vofit.reset(); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              //we are increment the coarse side so we have to get the faces on the
              //flip of sd
              Vector<FaceIndex> faces = ebisBox.getFaces(vof, a_dir, flip(a_sd));
              for (int iface = 0; iface < faces.size(); iface++)
                {
                  const FaceIndex& face = faces[iface];
                  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
                    {
                      Real area = ebisBox.areaFrac(face);
                      Real flux = a_coarFlux(face, ivar);
                      Real areaSum = area*flux;
                      //consistent with non-eb version
                      //              Real oldFlux = m_delUCoar[a_coarDatInd][a_dir](face, ivar);
                      areaSum *= -a_scale;
                      m_delUCoar[a_coarDatInd](vof, ivar) +=  areaSum;
                    } //end loop over variables
                }//end is this face in the holders
            }//loop over vofs
        }//Be bloody, bold and resolute
    } //for no man of woman born shall harm macbeth
}
/*******************/
void
EBFastFR::
incrementFineBoth(const EBFaceFAB&      a_fineFlux,
                  const Real&           a_scale,
                  const DataIndex&      a_fineDatInd,
                  const Interval&       a_variables,
                  const int&            a_dir,
                  const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFastFR::incrementFineBoth");
  incrementFineRegul(a_fineFlux,
                     a_scale,
                     a_fineDatInd,
                     a_variables,
                     a_dir,
                     a_sd);


  incrementFineIrreg(a_fineFlux,
                     a_scale,
                     a_fineDatInd,
                     a_variables,
                     a_dir,
                     a_sd);

}
/*******************/
void
EBFastFR::
incrementFineRegul(const EBFaceFAB&      a_fineFlux,
                   const Real&           a_scale,
                   const DataIndex&      a_fineDatInd,
                   const Interval&       a_variables,
                   const int&            a_dir,
                   const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFastFR::incrementFineRegul");
  //increment  as if there were no EB
  FArrayBox& fineFluxFAB = (FArrayBox&)(a_fineFlux.getFArrayBox());
  m_levelFluxReg->incrementFine(fineFluxFAB, a_scale, a_fineDatInd,
                                a_variables, a_variables, a_dir, a_sd);

  incrementFineSparse(a_fineFlux, a_scale, a_fineDatInd, a_variables,
                      a_dir, a_sd, true);
}
/*******************/
void
EBFastFR::
incrementFineIrreg(const EBFaceFAB&      a_fineFlux,
                   const Real&           a_scale,
                   const DataIndex&      a_fineDatInd,
                   const Interval&       a_variables,
                   const int&            a_dir,
                   const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFastFR::incrementFineIrreg");
  incrementFineSparse(a_fineFlux, a_scale, a_fineDatInd, a_variables,
                      a_dir, a_sd, false);

}// mmm deep loops

void
EBFastFR::
compareFineSparse(const EBFaceFAB&      a_fluxOld,
                  const EBFaceFAB&      a_fluxNew,
                  const DataIndex&      a_fineDatInd,
                  const int&            a_dir,
                  const Side::LoHiSide& a_sd)
{
  CH_TIME("EBFastFR::incrementFineSparse");
  if (m_hasEBCF)
    {
      int iindex = index(a_dir, a_sd);

      VoFIterator& vofit       = (m_vofiCoFi[iindex])[a_fineDatInd];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          //remember the registers live at the coar level
          const VolIndex& coarVoF = vofit();
          const EBISBox& ebisBoxCoFi = m_eblgCoFi.getEBISL()[a_fineDatInd];

          Vector<FaceIndex> facesCoar = ebisBoxCoFi.getFaces(coarVoF, a_dir, flip(a_sd));
          for (int ifacec = 0; ifacec < facesCoar.size(); ifacec++)
            {
              Vector<FaceIndex> facesFine =
                m_eblgCoFi.getEBISL().refine(facesCoar[ifacec], m_refRat, a_fineDatInd);
              for (int ifacef = 0; ifacef < facesFine.size(); ifacef++)
                {
                  VolIndex vofFine = facesFine[ifacef].getVoF(flip(a_sd));
                  for (int ivar = 0; ivar < a_fluxOld.nComp(); ivar++)
                    {
                      Real fluxOld  = a_fluxOld(facesFine[ifacef], ivar);
                      Real fluxNew  = a_fluxNew(facesFine[ifacef], ivar);
                      Real eps = 1.0e-9;
                      if (Abs(fluxOld-fluxNew) > eps)
                        {
                          pout() << "hurm at fine face" << facesFine[ifacef] << endl;
                        }
                    } //this is the way the world ends
                }//this is the way the world ends
            }//this is the way the world ends
        }//not with a bang
    } //but with a
}//whimper
/*******************/
void
EBFastFR::
incrementFineSparse(const EBFaceFAB&      a_fineFlux,
                    const Real&           a_scale,
                    const DataIndex&      a_fineDatInd,
                    const Interval&       a_variables,
                    const int&            a_dir,
                    const Side::LoHiSide& a_sd,
                    bool a_doingFineRegular)
{
  CH_TIME("EBFastFR::incrementFineSparse");
  if (m_hasEBCF)
    {
      //sign because of which side of the
      //vof we are on is opposite of which side
      // of the fine grid we are on
      //int isign = -sign(a_sd);
      //Real rsign = isign;
      Real newScale = a_scale/m_nrefdmo;
      int iindex = index(a_dir, a_sd);

      const EBISBox& ebisBoxFine= m_eblgFine.getEBISL() [a_fineDatInd];
      VoFIterator& vofit       = (m_vofiCoFi[iindex])[a_fineDatInd];
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          //remember the registers live at the coar level
          const VolIndex& coarVoF = vofit();
          const EBISBox& ebisBoxCoFi = m_eblgCoFi.getEBISL()[a_fineDatInd];

          EBCellFAB& delUCoFi = m_delUCoFi[a_fineDatInd];
          //we are on the coarse side of the face so we have to look back on it---ergo flip
          Vector<FaceIndex> facesCoar = ebisBoxCoFi.getFaces(coarVoF, a_dir, flip(a_sd));
          for (int ifacec = 0; ifacec < facesCoar.size(); ifacec++)
            {
              Vector<FaceIndex> facesFine =
                m_eblgCoFi.getEBISL().refine(facesCoar[ifacec], m_refRat, a_fineDatInd);
              for (int ifacef = 0; ifacef < facesFine.size(); ifacef++)
                {
                  Real area = ebisBoxFine.areaFrac(facesFine[ifacef]);
                  VolIndex vofFine = facesFine[ifacef].getVoF(flip(a_sd));
                  bool vofIsRegular = ebisBoxFine.isRegular(vofFine.gridIndex());
                  if (( a_doingFineRegular &&  vofIsRegular) || (!a_doingFineRegular && !vofIsRegular))
                    {
                      for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
                        {
                          Real flux          = a_fineFlux(facesFine[ifacef], ivar);
                          delUCoFi(coarVoF, ivar) +=  newScale*area*flux;
                        } //end loop over variables
                    }
                } //end loop over fine faces
            }// end loop over coarse faces
        }//end loop over EBCF vofs
    } //end if we have EBCF in this box
}// mmm deep loops
/*******************/
void
EBFastFR::
reflux(LevelData<EBCellFAB>& a_uCoar,
       const Interval&       a_variables,
       const Real&           a_scale,
       bool a_multByKappaOneMinusKappa)
{
  CH_TIME("EBFastFR::reflux");
  LevelData<FArrayBox> uCoarLDF;
  //save initial values of ucoar because the non-eb flux  reg will
  //change it in its ignorance
  if (m_hasEBCF)
    {
      cacheOldSolution(a_uCoar, a_variables);
    }

  //reflux as if there were no EB
  aliasEB(uCoarLDF, a_uCoar);
  m_levelFluxReg->reflux(uCoarLDF, a_variables, a_variables, a_scale);

  //correct at irregular cells
  if (m_hasEBCF)
    {
      restoreOldSolution(a_uCoar, a_variables);

      irregReflux(a_uCoar, a_variables, a_scale, a_multByKappaOneMinusKappa);
    }
}
/*******************/
void
EBFastFR::
irregReflux(LevelData<EBCellFAB>& a_uCoar,
            const Interval&       a_variables,
            const Real&           a_scale,
            bool a_multByKappaOneMinusKappa)
{
  CH_assert(m_isDefined);
  CH_TIME("EBFastFR::irregReflux");
  if (m_hasEBCF)
    {
      //coming into this routine, coar holds -coarFlux*area
      //and fine holds area*fineFlux
      //make diff = area(fineFlux-coarFlux)  (with scaling stuff)

      //dumpEBFFR(m_delUCoar, string("ducoar holds"));
      //for (int idir = 0; idir < SpaceDim; idir++)
      // {
      //   for (SideIterator sit; sit.ok(); ++sit)
      //     {
      //       int iside = sign(sit());
      //       pout() << "idir = " << idir << ", side = " << iside;
      //       dumpEBFFR(m_delUCoFi[index(idir, sit())], string(", duCoFi holds"));
      //     }
      // }


      EBLevelDataOps::setVal(m_delUDiff, 0.0);
      m_delUCoar.copyTo(a_variables, m_delUDiff, a_variables);

      //      dumpEBFFR(m_delUDiff, string("dudiff holds after first   copy"));


      EBAddOpEBFFR op;
      m_delUCoFi.copyTo(a_variables, m_delUDiff, a_variables, m_reverseCopier, op);

      //dumpEBFFR(m_delUDiff, string("dudiff holds after reverse copy"));

      //dumpEBFFR(a_uCoar, string("uCoar holds before reflux divergence"));

      //add refluxDivergence to solution u -= a_scale*(area*(fineFlux-coarFlux))
      incrementByRefluxDivergence(a_uCoar, m_delUDiff, a_variables, a_scale, false,
                                  a_multByKappaOneMinusKappa);

      //dumpEBFFR(a_uCoar, string("uCoar holds after reflux divergence"));

    }
}
/*******************/
void
EBFastFR::
incrementByRefluxDivergence(LevelData<EBCellFAB>& a_uCoar,
                            LevelData<EBCellFAB>& a_delUDiff,
                            const Interval      & a_variables,
                            const Real          & a_newScale,
                            bool a_multByOneMinusKappa,
                            bool a_multByKappaOneMinusKappa)
{
  if (m_hasEBCF)
    {
      //cannot both be true
      CH_assert(!a_multByOneMinusKappa  || !a_multByKappaOneMinusKappa);
      for (DataIterator dit= m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  int iindex = index(idir, sit());
                  const EBISBox& ebisBox       = m_eblgCoar.getEBISL()[dit()];
                  Vector<VoFIterator> vofitVec = m_vofiCoar[iindex][dit()];
                  for (int ivec = 0; ivec < vofitVec.size(); ivec++)
                    {
                      VoFIterator& vofit = vofitVec[ivec];
                      for (vofit.reset(); vofit.ok(); ++vofit)
                        {
                          const VolIndex& vof = vofit();
                          Real scale = sign(sit())*a_newScale;
                          Real volFrac = ebisBox.volFrac(vof);
                          if (a_multByOneMinusKappa)
                            {
                              scale *= (1.0-volFrac);
                            }
                          if (a_multByKappaOneMinusKappa)
                            {
                              scale *= volFrac*(1.0-volFrac);
                            }
                          for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++)
                            {
                              Real udelDif = a_delUDiff[dit()](vof, ivar);
                              Real soluVal = a_uCoar[dit()](vof, ivar);
                              Real soluNew = soluVal - scale*udelDif;

                              a_uCoar[dit()](vof, ivar) = soluNew;
                            }
                        }
                    }
                }
            }
        }
    }

}
/*******************/
void
EBFastFR::
restoreOldSolution(LevelData<EBCellFAB>&       a_uCoar,
                   const Interval&             a_variables)
{

  CH_TIME("EBFastFR::saveOldSolution");
  for (DataIterator dit = m_eblgCoar.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              int iindex = index(idir, sit());
              for (int ivec = 0; ivec < (m_vofiCoar[iindex])[dit()].size(); ivec++)
                {
                  VoFIterator& vofit = (m_vofiCoar[iindex])[dit()][ivec];
                  for (vofit.reset(); vofit.ok(); ++vofit)
                    {
                      for (int icomp = a_variables.begin(); icomp <= a_variables.end(); icomp++)
                        {
                          a_uCoar[dit()](vofit(), icomp) =  m_saveCoar[dit()](vofit(), icomp) ;
                        }
                    }
                }
            }
        }
    }
}
/*******************/
void
EBFastFR::
cacheOldSolution(const LevelData<EBCellFAB>& a_uCoar,
                 const Interval&             a_variables)
{
  CH_TIME("EBFastFR::saveOldSolution");
  a_uCoar.copyTo(a_variables, m_saveCoar, a_variables);
}
/*******************/
void
EBFastFR::
incrementDensityArray(LevelData<EBCellFAB>& a_coarDense,
                      const Interval&       a_variables,
                      const Real&           a_scale)
{

  //only do stuff on irregular cells because on non-irregular cells 1-kappa=0
  //does not need to be fast because this is a debugging tool
  if (m_hasEBCF)
    {
      //coming into this routine, coarFlux holds -coarFlux*area
      //and fineFlux holds area*fineFlux
      //make fluxDiff = area(fineFlux-coarFlux)  (with scaling stuff)
      EBLevelDataOps::setToZero(m_delUDiff);

      m_delUCoar.copyTo(a_variables, m_delUDiff, a_variables);
      EBAddOpEBFFR op;

      m_delUCoFi.copyTo(a_variables, m_delUDiff, a_variables, m_reverseCopier, op);

      //add refluxDivergence to solution u = (1-kappa)*a_scale*(area*(fineFlux-coarFlux))
      incrementByRefluxDivergence(a_coarDense, m_delUDiff, a_variables, a_scale, true, false);
    }
}
/*************************/
#include "NamespaceFooter.H"
