#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeQCFI.cpp
// petermc, Mon, Apr 23, 2001
// petermc, 1 Nov 2001, redid this so that it uses NodeQCFI2
#include "NodeQCFI.H"
#include "NamespaceHeader.H"

// ---------------------------------------------------------
void
NodeQCFI::setDefaultValues()
{
  m_isDefined = false;
  m_verbose = false;
}

// ---------------------------------------------------------
NodeQCFI::NodeQCFI()
{
  setDefaultValues();
}

// ---------------------------------------------------------
NodeQCFI::NodeQCFI(const DisjointBoxLayout& a_grids,
                   Real a_dx,
                   const Box& a_domain,
                   const LayoutData<NodeCFIVS>* const a_loCFIVS,
                   const LayoutData<NodeCFIVS>* const a_hiCFIVS,
                   int a_refToCoarse,
                   NodeBCFunc               a_bc,
                   int a_interpolationDegree,
                   int a_ncomp,
                   bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  setDefaultValues();
  define(a_grids, a_dx, probdomain,
         a_loCFIVS, a_hiCFIVS,
         a_refToCoarse, a_bc, a_interpolationDegree, a_ncomp, a_verbose);
}

// ---------------------------------------------------------
NodeQCFI::NodeQCFI(const DisjointBoxLayout& a_grids,
                   Real a_dx,
                   const ProblemDomain& a_domain,
                   const LayoutData<NodeCFIVS>* const a_loCFIVS,
                   const LayoutData<NodeCFIVS>* const a_hiCFIVS,
                   int a_refToCoarse,
                   NodeBCFunc               a_bc,
                   int a_interpolationDegree,
                   int a_ncomp,
                   bool a_verbose)
{
  setDefaultValues();
  define(a_grids, a_dx, a_domain,
         a_loCFIVS, a_hiCFIVS,
         a_refToCoarse, a_bc, a_interpolationDegree, a_ncomp, a_verbose);
}

// ---------------------------------------------------------
NodeQCFI::~NodeQCFI()
{
  clearMemory();
}

// ---------------------------------------------------------
bool
NodeQCFI::isDefined() const
{
  return(m_isDefined);
}

// ---------------------------------------------------------
void
NodeQCFI::define(const DisjointBoxLayout& a_grids,
                 Real a_dx,
                 const Box& a_domain,
                 const LayoutData<NodeCFIVS>* const a_loCFIVS,
                 const LayoutData<NodeCFIVS>* const a_hiCFIVS,
                 int a_refToCoarse,
                 NodeBCFunc               a_bc,
                 int a_interpolationDegree,
                 int a_ncomp,
                 bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  define(a_grids, a_dx, probdomain,
         a_loCFIVS, a_hiCFIVS,
         a_refToCoarse, a_bc, a_interpolationDegree, a_ncomp, a_verbose);
}

// ---------------------------------------------------------
void
NodeQCFI::define(const DisjointBoxLayout& a_grids,
                 Real a_dx,
                 const ProblemDomain& a_domain,
                 const LayoutData<NodeCFIVS>* const a_loCFIVS,
                 const LayoutData<NodeCFIVS>* const a_hiCFIVS,
                 int a_refToCoarse,
                 NodeBCFunc               a_bc,
                 int a_interpolationDegree,
                 int a_ncomp,
                 bool a_verbose)
{
  CH_assert(a_refToCoarse >= 2);
  m_bc = a_bc;
  m_grids = a_grids;
  m_dx = a_dx;
  m_verbose = a_verbose;

  // m_coarsenings == log2(a_refToCoarse);
  m_coarsenings = 0;
  for (int interRatio = a_refToCoarse; interRatio >= 2; interRatio /= 2)
    m_coarsenings++;

  m_qcfi2.resize(m_coarsenings);
  m_inter.resize(m_coarsenings-1);
  m_loCFIVScoarser.resize(m_coarsenings-1);
  m_hiCFIVScoarser.resize(m_coarsenings-1);
  // for i=1:m_c-2, m_qcfi2[i] interpolates m_inter[i] from m_inter[i-1].

  if (m_coarsenings == 1)
    {
      m_qcfi2[0] = new NodeQuadCFInterp2(a_grids, a_domain,
                                         a_loCFIVS, a_hiCFIVS,
                                         false, a_interpolationDegree, a_ncomp);
    }
  else
    {
      // if m_coarsenings == 2:
      // m_qcfi2[1] = new NodeQuadCFInterp2(a_grids, true);
      // m_qcfi2[0] = new NodeQuadCFInterp2(a_grids.coarsen(2), false);
      // m_inter[0] = new LevelData(a_grids.coarsen(2));

      ProblemDomain domainLevel(a_domain);
      Real dxLevel = m_dx;
      DisjointBoxLayout gridsInterFine(a_grids);
      DisjointBoxLayout gridsInterCoarse;
      // m_qcfi2[m_coarsenings-1] interpolates from a 2-coarsening of
      // the fine grids onto the fine grids (a_grids).
      // Interpolation is from the interface only.
      m_qcfi2[m_coarsenings-1] =
        new NodeQuadCFInterp2(a_grids, domainLevel,
                              a_loCFIVS, a_hiCFIVS,
                              true, a_interpolationDegree, a_ncomp);
      for (int interlev = m_coarsenings - 2; interlev >= 0; interlev--)
        {
          // m_coarsenings >= 2, so this will be executed at least once.
          coarsen(gridsInterCoarse, gridsInterFine, 2);
          domainLevel.coarsen(2);
          dxLevel *= 2.;
          bool interfaceOnly = (interlev > 0);
          // m_qcfi2[interlev] interpolates from some coarsening of
          // the fine grids onto gridsInterCoarse.
          // Except for the first time (interlev == 0),
          // interpolation is from the interface only.

          // Define objects containing the nodes of the coarse/fine interfaces.
          m_loCFIVScoarser[interlev] = new LayoutData<NodeCFIVS>[SpaceDim];
          m_hiCFIVScoarser[interlev] = new LayoutData<NodeCFIVS>[SpaceDim];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              LayoutData<NodeCFIVS>& loCFIVS =
                m_loCFIVScoarser[interlev][idir];
              LayoutData<NodeCFIVS>& hiCFIVS =
                m_hiCFIVScoarser[interlev][idir];

              loCFIVS.define(gridsInterCoarse);
              hiCFIVS.define(gridsInterCoarse);
              for (DataIterator dit = gridsInterCoarse.dataIterator(); dit.ok(); ++dit)
                {
                  const Box& bx = gridsInterCoarse.get(dit());
                  loCFIVS[dit()].define(domainLevel, bx, gridsInterCoarse,
                                        idir, Side::Lo);
                  hiCFIVS[dit()].define(domainLevel, bx, gridsInterCoarse,
                                        idir, Side::Hi);
                }
            }

          m_qcfi2[interlev] =
            new NodeQuadCFInterp2(gridsInterCoarse, domainLevel,
                                  m_loCFIVScoarser[interlev],
                                  m_hiCFIVScoarser[interlev],
                                  interfaceOnly, a_interpolationDegree, a_ncomp);

          m_inter[interlev] =
            new LevelData<NodeFArrayBox>(gridsInterCoarse, a_ncomp, IntVect::Zero);
          gridsInterFine = gridsInterCoarse;
        }
      m_domainPenultimate = domainLevel;
      m_dxPenultimate = dxLevel;
    }

  m_ncomp = a_ncomp;

  m_isDefined = true;
}

// ---------------------------------------------------------
void
NodeQCFI::clearMemory()
{
  if (isDefined())
    {
      for (int interlev = 0; interlev < m_coarsenings; interlev++)
        {
          delete m_qcfi2[interlev];
        }
      for (int interlev = 0; interlev < m_coarsenings-1; interlev++)
        {
          delete m_inter[interlev];
          delete [] m_loCFIVScoarser[interlev];
          delete [] m_hiCFIVScoarser[interlev];
        }
    }
  m_refToCoarse = -1;
  m_isDefined = false;
}

// ---------------------------------------------------------
void
NodeQCFI::setVerbose(bool a_verbose)
{
  m_verbose = a_verbose;
}

// ---------------------------------------------------------
void
NodeQCFI::coarseFineInterp(LevelData<NodeFArrayBox>& a_phiFine,
                           const LevelData<NodeFArrayBox>& a_phiCoarse,
                           bool a_inhomogeneous)
{
  CH_assert(isDefined());
  CH_assert(a_phiFine.nComp() == m_ncomp);
  CH_assert(a_phiCoarse.nComp() == m_ncomp);

  if (m_coarsenings == 1)
    {
      m_qcfi2[0]->coarseFineInterp(a_phiFine, a_phiCoarse);
    }
  else // m_coarsenings >= 2
    {
      // if m_coarsenings == 2:
      // m_qcfi2[1] = new NodeQuadCFInterp2(a_grids, true);
      // m_qcfi2[0] = new NodeQuadCFInterp2(a_grids.coarsen(2), false);
      // m_inter[0] = new LevelData(a_grids.coarsen(2));

      // m_qcfi2[0] interpolates m_inter[0] from a_phiCoarse on all of it.
      // m_qcfi2[1] interpolates a_phiFine from m_inter[0] on interface only.

      m_qcfi2[0]->coarseFineInterp(*m_inter[0], a_phiCoarse);
      m_inter[0]->exchange(m_inter[0]->interval());

      // Set domain and dx for coarsest level refined by 2.
      ProblemDomain domain(m_domainPenultimate);
      Real dx = m_dxPenultimate;

      bool homogeneous = !a_inhomogeneous;
      const DisjointBoxLayout& dbl0 = m_inter[0]->disjointBoxLayout();
      for (DataIterator dit = dbl0.dataIterator(); dit.ok(); ++dit)
        {
          m_bc((*m_inter[0])[dit()],  dbl0.get(dit()), domain,  dx, homogeneous);
        }

      for (int interlev = 1; interlev < m_coarsenings - 1; interlev++)
        {
          m_qcfi2[interlev]->coarseFineInterp(*m_inter[interlev],
                                              *m_inter[interlev-1]);

          m_inter[interlev]->exchange(m_inter[interlev]->interval());

          domain.refine(2);
          dx *= 0.5;
          const DisjointBoxLayout& dbl = m_inter[interlev]->disjointBoxLayout();
          for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
            {
              m_bc((*m_inter[interlev])[dit()],  dbl.get(dit()), domain,  dx, homogeneous);
            }
        }

      m_qcfi2[m_coarsenings-1]->coarseFineInterp(a_phiFine,
                                                 *m_inter[m_coarsenings-2]);
      // Don't need to set boundary conditions for a_phiFine
      // because the calling function will do that.
    }
}

#include "NamespaceFooter.H"
