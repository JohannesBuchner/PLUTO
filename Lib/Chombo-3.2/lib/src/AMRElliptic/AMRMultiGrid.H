#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// BVS, June 18, 2003

// We can assume that template class T has null construction.

#ifndef _AMRMULTIGRID_H_
#define _AMRMULTIGRID_H_

#include "MultiGrid.H"
#include "REAL.H"
#include "Box.H"
#include "NoOpSolver.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "Copier.H"
#include "SPMD.H"
#include "Misc.H"
#include "FArrayBox.H"
#include "LevelDataOps.H"
#include "CornerCopier.H"

#include "NamespaceHeader.H"

///
/**
   Operator class for AMR Multigrid
 */
template <typename T>
class AMRLevelOp : public MGLevelOp<T>
{
public:

  ///
  virtual void dumpAMR(Vector<T*>& a_data, string name)
  {
  }

  virtual void dumpLevel(T& a_data, string name)
  {
  }

  //! Constructor.
  AMRLevelOp()
    :MGLevelOp<T>()
  {
  }

  virtual void dumpStuff(Vector<T*> data, string filename)
  {
  }

  //! Destructor.
  virtual ~AMRLevelOp()
  {
  }

  virtual Real AMRNorm(const T&   a_coarResid,
                       const T&   a_fineResid,
                       const int& a_refRat,
                       const int& a_ord)
  {
    return this->norm(a_coarResid, 0);
  }

  //debugging hook
  virtual void outputLevel(T& a_rhs, string& a_name)
  {
  }

  //debugging hook
  virtual void outputAMR(Vector<T*>& a_rhs, string& a_name)
  {
  }
  ///
  /**
     return the refinement ratio to next coarser level.
     return 1 when there are no coarser AMRLevelOp objects.
   */
  virtual int refToCoarser() = 0;

  ///
  /**
      a_residual = a_rhs - L(a_phiFine, a_phi, a_phiCoarse)
  */
  virtual void AMRResidual(T& a_residual, const T& a_phiFine, const T& a_phi,
                           const T& a_phiCoarse, const T& a_rhs,
                           bool a_homogeneousDomBC,
                           AMRLevelOp<T>*  a_finerOp) = 0;

  ///
  /**
      a_residual = a_rhs - L^nf(a_phi, a_phiCoarse)
      assume no finer AMR level
  */
  virtual void AMRResidualNF(T& a_residual, const T& a_phi, const T& a_phiCoarse,
                             const T& a_rhs, bool a_homogeneousBC) = 0;

  ///
  /**
      a_residual = a_rhs - L(a_phiFine, a_phi)
      assume no coarser AMR level
  */
  virtual void AMRResidualNC(T& a_residual, const T& a_phiFine, const T& a_phi,
                             const T& a_rhs, bool a_homogeneousBC,
                             AMRLevelOp<T>* a_finerOp) = 0;

  ///
  /**
     Apply the AMR operator, including coarse-fine matching
  */
  virtual void AMROperator(T& a_LofPhi,
                           const T& a_phiFine, const T& a_phi,
                           const T& a_phiCoarse,
                           bool a_homogeneousDomBC,
                           AMRLevelOp<T>*  a_finerOp) = 0;

  ///
  /**
      Apply the AMR operator, including coarse-fine matching.
      assume no finer AMR level
  */
  virtual void AMROperatorNF(T& a_LofPhi,
                             const T& a_phi,
                             const T& a_phiCoarse,
                             bool a_homogeneousBC) = 0;

  ///
  /**
      Apply the AMR operator, including coarse-fine matching
      assume no coarser AMR level
  */
  virtual void AMROperatorNC(T& a_LofPhi,
                             const T& a_phiFine,
                             const T& a_phi,
                             bool a_homogeneousBC,
                             AMRLevelOp<T>* a_finerOp) = 0;

  ///
  /** a_resCoarse = I[h-2h]( a_residual - L(a_correction, a_coarseCorrection)) */
  virtual void AMRRestrict(T& a_resCoarse, const T& a_residual, const T& a_correction,
                           const T& a_coarseCorrection, bool a_skip_res) = 0;

  ///
  /** a_correction += I[2h->h](a_coarseCorrection) */
  virtual void AMRProlong(T& a_correction, const T& a_coarseCorrection) = 0;

  ///
  /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
  virtual void AMRUpdateResidual(T& a_residual, const T& a_correction,
                                 const T& a_coarseCorrection) = 0;

  ///
  /**
   */
  virtual void createCoarsened(T&       a_lhs,
                               const T& a_rhs,
                               const int&     a_refRat) = 0;

  //===================================================================
  // optional optimizations for an AMRLevelOp.  These are not pure virtual
  // functions, since we can build the equivalent algorithmic components from
  // pure virtual functions.  The AMRMultiGrid algorithm actually calls *these*
  // functions, which a smart operator can perform in faster ways.
  //===================================================================

  virtual void buildCopier(Copier& a_copier, const T& a_lhs, const T& a_rhs)
  {
  }

  virtual void assignCopier(T& a_lhs, const T& a_rhs, const Copier& a_copier)
  {
    this->assign(a_lhs, a_rhs);
  }

  virtual void zeroCovered(T& a_lhs, T& a_rhs, const Copier& a_copier)
  {
    this->setToZero(a_rhs);
    this->assignCopier(a_lhs, a_rhs, a_copier);
  }
  virtual Real localMaxNorm(const T& a_phi)
  {
    return this->norm(a_phi, 0);
  }

  /** optimization of AMRProlong that sends in the existing temporary and copier */
  virtual void AMRProlongS(T& a_correction, const T& a_coarseCorrection,
                           T& a_temp, const Copier& a_copier)
  {
    AMRProlong(a_correction, a_coarseCorrection);
  }

  /** optimization of AMRProlong that sends in the existing temporary and copier -- high order */
  virtual void AMRProlongS_2( T& a_correction, const T& a_coarseCorrection,
                              T& a_temp, const Copier& a_copier,
                              const Copier& a_cornerCopier,
                              const AMRLevelOp<LevelData<FArrayBox> >*  a_crsOp )
  {
    AMRProlong( a_correction, a_coarseCorrection );
  }

  virtual void AMRRestrictS(T& a_resCoarse, const T& a_residual, const T& a_correction,
                            const T& a_coarseCorrection, T& scratch, bool a_skip_res = false )
  {
    AMRRestrict(a_resCoarse, a_residual, a_correction, a_coarseCorrection, a_skip_res);
  }

  virtual unsigned int orderOfAccuracy(void) const
  {
    return 2;
  }

  /// This routine is for operators with orderOfAccuracy()>2.
  virtual void enforceCFConsistency(T& a_coarseCorrection, const T& a_correction)
  {
  }
};

///
/**
   Factory to create AMRLevelOps
 */
template <class T>
class AMRLevelOpFactory : public MGLevelOpFactory<T>
{
public:
  virtual ~AMRLevelOpFactory()
  {
  }

  ///
  /**
     return a new operator.  this is done with a new call.
     caller is responsible for deletion
   */
  virtual AMRLevelOp<T>* AMRnewOp(const ProblemDomain& a_indexSpace)=0;

  ///
  /**
     return refinement ratio  to next finer level.
   */
  virtual int refToFiner(const ProblemDomain& a_indexSpace) const =0;

};

//! \class AMRMultiGridInspector
//! This base class allows one to construct methods for inspecting the
//! multigrid algorithm at each of its steps.
template <class T>
class AMRMultiGridInspector
{
  public:

  //! Base class constructor. This must be called by all subclasses.
  AMRMultiGridInspector()
  {
  }

  //! Destructor.
  virtual ~AMRMultiGridInspector()
  {
  }

  //! Override this method to keep track of a multigrid residual computed
  //! during a multigrid iteration at the given levels.
  //! \param a_residuals An array containing the residuals computed by the multigrid
  //!                    algorithm at each level in the range [a_minLevel, a_maxLevel].
  //! \param a_minLevel The lowest AMR level at which residuals were computed.
  //! \param a_maxLevel The highest AMR level at which residuals were computed.
  //! \param a_iter The multigrid iteration number.
  virtual void recordResiduals(const Vector<T*>& a_residuals,
                               int a_minLevel,
                               int a_maxLevel,
                               int a_iter) = 0;

  //! Override this method to keep track of a multigrid correction computed
  //! during a V cycle at the given level.
  //! \param a_corrections An array containing the corrections computed during a V cycle
  //!                      at each level in the range [a_minLevel, a_maxLevel].
  //! \param a_minLevel The lowest AMR level at which corrections were computed.
  //! \param a_maxLevel The highest AMR level at which corrections were computed.
  //! \param a_iter The multigrid iteration number.
  virtual void recordCorrections(const Vector<T*>& a_corrections,
                                 int a_minLevel,
                                 int a_maxLevel,
                                 int a_iter) = 0;

  private:

  AMRMultiGridInspector(const AMRMultiGridInspector&);
  AMRMultiGridInspector& operator=(const AMRMultiGridInspector&);
};

///
/**
   Class to solve elliptic equations using the Martin and Cartwright algorithm.
 */
template <class T>
class AMRMultiGrid
{
public:

  AMRMultiGrid();

  virtual ~AMRMultiGrid();

  void outputAMR(Vector<T*>& a_data, string& a_name, int a_lmax, int a_lbase);

  ///
  /**
     Define the solver.
     a_coarseDomain is the index space on the coarsest AMR level.
     a_factory is the operator factory through which all special information is conveyed.
     a_bottomSolver is the solver to be used at the termination of multigrid coarsening.
       It is the client's responsibility to free up the dynamically-allocated memory.
     a_numLevels is the number of AMR levels.
   */
  virtual void define(const ProblemDomain& a_coarseDomain,
                      AMRLevelOpFactory<T>& a_factory,
                      LinearSolver<T>* a_bottomSolver,
                      int a_numLevels);

  //! Add an inspector to the list of inspectors maintained by this AMRMultiGrid
  //! instance. It will be given the opportunity to record intermediate data.
  void addInspector(RefCountedPtr<AMRMultiGridInspector<T> >& a_inspector)
  {
    CH_assert(!a_inspector.isNull());
    m_inspectors.push_back(a_inspector);
  }

  ///
  /**
     Solve L(phi) = rho from l_base to l_max.  To solve over all levels,
     l_base = 0 and l_max = max_level = numLevels-1.
   */
  virtual void solve( Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                      int l_max, int l_base, bool a_zeroPhi=true,
                      bool forceHomogeneous = false );

  ///
  /** same as "solve" except user has taken the reponsibility of having previously
      called "init" so solver can allocate temporary holders.
  */
  virtual void solveNoInit(Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                           int l_max, int l_base, bool a_zeroPhi=true,
                           bool forceHomogeneous = false);

  ///use if you want final residual
  virtual void solveNoInitResid(Vector<T*>& a_phi, Vector<T*>& a_finalResid,
                                const Vector<T*>& a_rhs,
                                int l_max, int l_base, bool a_zeroPhi=true,
                                bool forceHomogeneous = false);

  void relaxOnlyHomogeneous(Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                            int l_max, int l_base);

  virtual void AMRVCycle(Vector<T*>& a_correction,
                         Vector<T*>& a_residual,
                         int l, int l_max, int l_base);

  void setMGCycle(int a_numMG);

  virtual void init(const Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                    int l_max, int l_base);
  //init messes with multigrid depth.  this puts it back
  void revert(const Vector<T*>& a_phi, const Vector<T*>& a_rhs,
            int l_max, int l_base);

  Real m_eps, m_hang, m_normThresh;
  bool   m_solverParamsSet;
  int m_imin, m_iterMax, m_iterMin, m_verbosity, m_exitStatus;
  int m_pre, m_post, m_bottom, m_numMG;
  /// max no. of coarsenings -- -1 (default) means coarsen as far as possible
  /** If using a value besides the default, need to set it _before_
      define function is called */
  int m_maxDepth;
  AMRLevelOp<T>& levelOp(int level);
  // default m_convergenceMetric = 0.:  initial residual will be set to
  // result of computeAMRResidual.
  // if m_convergenceMetric > 0., then initial residual will be set to
  // m_convergenceMetric.
  Real m_convergenceMetric;

  // used to give an additional cushion in the EPS used for bottom solves
  Real m_bottomSolverEpsCushion;

  ///
  /**
     resid = L(phi) - rhs
   */
  Real computeAMRResidual(Vector<T*>&       a_resid,
                          Vector<T*>&       a_phi,
                          const Vector<T*>& a_rhs,
                          int               l_max,
                          int               l_base,
                          bool              a_homogeneousBC=false,
                          bool              a_computeNorm=true);

  /** just return the normed value of computeAMRResidual.  used for benchmarking */
  Real computeAMRResidual(Vector<T*>&      a_phi,
                          const Vector<T*>& a_rhs,
                          int l_max,
                          int l_min);

  ///
  /**
     lph = L(phi)
   */
  void computeAMROperator(Vector<T*>&       a_lph,
                          Vector<T*>&       a_phi,
                          int               l_max,
                          int               l_base,
                          bool              a_homogeneousBC=false);

  ///
  /**
     For changing coefficients.  Use at thy own peril.
  */
  Vector< MGLevelOp<T>* > getAllOperators();
  Vector< MGLevelOp<T>* > getOperatorsOp();
  Vector< Vector< MGLevelOp<T>* > > getOperatorsMG();
  Vector< AMRLevelOp<T> * >& getAMROperators()
  {
    return m_op;
  }

  ///
  /**
     Set parameters of the solve.
     a_pre is the number of smoothings before averaging.
     a_post is the number of smoothings after averaging.
     a_bottom is the number of smoothings at the bottom level.
     a_numMG = 1 for vcycle, =2 for wcycle (use 1).
     a_itermax is the max number of v cycles.
     a_hang is the minimum amount of change per vcycle.
     a_eps is the solution tolerance.
     a_normThresh is how close to zero eps*resid is allowed to get.
   */
  void setSolverParameters(const int&   a_pre,
                           const int&   a_post,
                           const int&   a_bottom,
                           const int&   a_numMG,
                           const int&   a_iterMax,
                           const Real&  a_eps,
                           const Real&  a_hang,
                           const Real&  a_normThresh);

  /// set up bottom solver for internal MG solver
  /**
     This function is normally called by the solve(...) function.
     However, it must be called if solve will not be called (in particular,
     if only the V-cycle is being used)
  */
  void setBottomSolver(int l_max, int l_base);

  void setBottomSolverEpsCushion(Real a_bottomSolverEpsCushion);

protected:

  void relax(T& phi, T& R, int depth, int nRelax = 2);

  void computeAMRResidualLevel(Vector<T*>&       a_resid,
                               Vector<T*>&       a_phi,
                               const Vector<T*>& a_rhs,
                               int l_max, int l_base, int ilev,
                               bool a_homogeneousBC);

  Vector<AMRLevelOp<T>*>          m_op;
  Vector<MultiGrid<T> *>          m_mg;
  Vector<T*>  m_correction;
  Vector<T*>  m_residual;
  Vector<T*>  m_resC;
  Vector<Copier> m_resCopier;
  Vector<Copier> m_reverseCopier;

  NoOpSolver<T>    m_nosolve;

  LinearSolver<T>* m_bottomSolver;

  Vector<char> m_hasInitBeenCalled;

  void clear();

private:

  // A list of inspectors maintained by this instance.
  Vector<RefCountedPtr<AMRMultiGridInspector<T> > > m_inspectors;

  // Forbidden copiers.
  AMRMultiGrid(const AMRMultiGrid<T>&);
  AMRMultiGrid& operator=(const AMRMultiGrid<T>&);
};

//*******************************************************
// AMRMultigrid Implementation
//*******************************************************

//===================================================================

template <class T>
void
AMRMultiGrid<T>::outputAMR(Vector<T*>& a_data, string& a_name, int a_lmax, int a_lbase)
{
  Vector<T*> outputData;
  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      outputData.push_back(a_data[ilev]);
    }
  m_op[a_lbase]->outputAMR(outputData, a_name);
}

template <class T>
void
AMRMultiGrid<T>::setSolverParameters(const int&   a_pre,
                                     const int&   a_post,
                                     const int&   a_bottom,
                                     const int&   a_numMG,
                                     const int&   a_iterMax,
                                     const Real&  a_eps,
                                     const Real&  a_hang,
                                     const Real&  a_normThresh)
{
  m_solverParamsSet = true;
  m_pre        =    a_pre;
  m_post       =    a_post;
  m_bottom     =    a_bottom;
  m_eps        =    a_eps;
  m_hang       =    a_hang;
  m_normThresh =    a_normThresh;
  m_iterMax    =    a_iterMax;
  for (int img = 0; img < m_mg.size(); img++)
    {
      m_mg[img]->m_pre    = a_pre;
      m_mg[img]->m_post   = a_post;
      m_mg[img]->m_bottom   = a_bottom;
    }
  setMGCycle(a_numMG);
  m_bottomSolverEpsCushion = 1.0;
}
template <class T>
Vector< MGLevelOp<T> * >
AMRMultiGrid<T>::getAllOperators()
{
  Vector< MGLevelOp<T>* > retval;
  for (int iop = 0;  iop < m_op.size(); iop++)
    {
      MGLevelOp<T>* operPtr = (MGLevelOp<T>*) m_op[iop];
      retval.push_back(operPtr);
    }

  for (int img = 0; img < m_mg.size(); img++)
    {
      Vector< MGLevelOp<T>* > mgOps = m_mg[img]->getAllOperators();
      retval.append(mgOps);
    }
  return retval;
}

template <class T>
Vector< MGLevelOp<T> * >
AMRMultiGrid<T>::getOperatorsOp()
{
  Vector< MGLevelOp<T>* > retval;
  for (int iop = 0;  iop < m_op.size(); iop++)
    {
      MGLevelOp<T>* operPtr = (MGLevelOp<T>*) m_op[iop];
      retval.push_back(operPtr);
    }
  return retval;
}

template <class T>
Vector<Vector< MGLevelOp<T> * > >
AMRMultiGrid<T>::getOperatorsMG()
{
  Vector< Vector< MGLevelOp<T>* > > retval(m_mg.size());

  for (int img = 0; img < m_mg.size(); img++)
    {
      retval[img] = m_mg[img]->getAllOperators();
    }
  return retval;
}

template <class T>
AMRMultiGrid<T>::AMRMultiGrid()
 :m_eps(1E-6),
  m_hang(1E-15),
  m_normThresh(1E-30),
  m_imin(5),
  m_iterMax(20),
  m_iterMin(-1),
  m_verbosity(3),
  m_pre(2),
  m_post(2),
  m_bottom(2),
  m_numMG(1),
  m_maxDepth(-1),
  m_convergenceMetric(0.),
  m_bottomSolverEpsCushion(1.0),
  m_bottomSolver(NULL),
  m_inspectors()
{
  m_solverParamsSet = false;

  //m_maxDepth = 4;
  //pout() << "       AMRMultiGrid<T>::AMRMultiGrid() = " << m_maxDepth << endl;
}
template <class T>
void AMRMultiGrid<T>::setMGCycle(int a_numMG)
{
  for (int ilev = 0; ilev < m_op.size(); ilev++)
    {
      m_mg[ilev]->m_numMG = a_numMG;
      m_mg[ilev]->m_cycle = a_numMG;
    }
  m_numMG = a_numMG;
}
template <class T>
void AMRMultiGrid<T>::relax(T& a_correction, T& a_residual, int depth, int a_numSmooth)
{
  CH_TIME("AMRMultiGrid::relax");

  if (m_op[depth]->refToCoarser() > 2)
    {
      // intermediate multigrid levels exist between this level and the
      // next coarser AMR level, so do a mini v-cyle to smooth on those

      int intermediateDepth = 0;
      int r = m_op[depth]->refToCoarser();
      while (r>2)
        {
          r/=2;
          intermediateDepth++;
        }
      int tmp = m_mg[depth]->m_depth;
      m_mg[depth]->m_depth = intermediateDepth;
      //ok, there is an intermediate multigrid level that is
      // not an AMR level.  We use regular MultiGrid smoothing algorithm for this.
      // note that MultiGrid::cycle handles relaxation on this level
      // as well, so there's no need to call the LinearOp's relax function.
      m_mg[depth]->cycle(0, a_correction, a_residual);
      m_mg[depth]->m_depth = tmp;
    }
  else
    {
      // no intermediate multigrid levels between AMR levels, so
      // just call relax on this level and be done with it.
      m_op[depth]->relax(a_correction, a_residual, a_numSmooth);  //numSmoothDown
    }

}
/************/
template <class T>
AMRMultiGrid<T>::~AMRMultiGrid()
{
  CH_TIME("~AMRMultiGrid");
  clear();
}

template <class T>
AMRLevelOp<T>&  AMRMultiGrid<T>::levelOp(int level)
{
  return *(m_op[level]);
}

/************/
template <class T>
Real AMRMultiGrid<T>::computeAMRResidual(Vector<T*>&       a_resid,
                                         Vector<T*>&       a_phi,
                                         const Vector<T*>& a_rhs,
                                         int               l_max,
                                         int               l_base,
                                         bool              a_homogeneousBC,
                                         bool              a_computeNorm)
{
  CH_TIME("AMRMultiGrid::computeAMRResidual");

  Real rnorm = 0;
  Real localNorm = 0;
  for (int ilev = l_base; ilev <= l_max; ilev++)
    {
      //always used at top level where bcs are inhomogeneous
      computeAMRResidualLevel(a_resid,
                              a_phi,
                              a_rhs,
                              l_max, l_base, ilev, a_homogeneousBC);
      if (a_computeNorm)
        {
          if (ilev == l_max)
            {
              localNorm =m_op[ilev]->localMaxNorm(*a_resid[ilev]);
            }
          else
            {
              m_op[ilev]->zeroCovered(*a_resid[ilev], *m_resC[ilev+1], m_resCopier[ilev+1]);
              localNorm = m_op[ilev]->localMaxNorm(*a_resid[ilev]);
            }
          const int normType = 0;
          if (normType==2)
            localNorm = m_op[ilev]->norm(*a_resid[ilev],normType);
          rnorm = Max(localNorm, rnorm);
        }
    }
#ifdef CH_MPI
  if (a_computeNorm)
    {
      CH_TIME("MPI_Allreduce");
      Real recv;
      int result = MPI_Allreduce(&rnorm, &recv, 1, MPI_CH_REAL,
                                 MPI_MAX, Chombo_MPI::comm);
      if (result != MPI_SUCCESS)
      {
        //bark!!!
        MayDay::Error("sorry, but I had a communcation error on norm");
      }
      rnorm = recv;
    }
#endif

  return rnorm; // if a_computeNorm is false, then this just returns zero.
}

template <class T>
Real AMRMultiGrid<T>::computeAMRResidual(Vector<T*>&       a_phi,
                                         const Vector<T*>&       a_rhs,
                                         int l_max,
                                         int l_min)
{
  return computeAMRResidual(m_residual,
                           a_phi,
                           a_rhs,
                           l_max,
                           l_min,
                           false,
                           true);
}

/************/
template <class T>
void AMRMultiGrid<T>::computeAMROperator(Vector<T*>&       a_lph,
                                         Vector<T*>&       a_phi,
                                         int               l_max,
                                         int               l_base,
                                         bool              a_homogeneousBC)
{
  for (int ilev = l_base; ilev <= l_max; ilev++)
    {
      m_op[ilev]->create(*m_residual[ilev], *a_lph[ilev]);
      m_op[ilev]->setToZero(*m_residual[ilev]);
    }
  computeAMRResidual(a_lph, a_phi, m_residual, l_max, l_base, a_homogeneousBC, false);
  // fixing the bug of not negating the result --- Qinghai Zhang 12/10/2009
  // In loving memory of the two weeks I lost over this bug!!!
  for (int ilev = l_base; ilev <= l_max; ilev++)
    {
      m_op[ilev]->scale(*a_lph[ilev], -1.);
    }
}

/************/
template <class T>
void AMRMultiGrid<T>::computeAMRResidualLevel(Vector<T*>&       a_resid,
                                              Vector<T*>&       a_phi,
                                              const Vector<T*>& a_rhs,
                                              int l_max, int l_base, int ilev,
                                              bool a_homogeneousBC)
{
  CH_TIME("AMRMultiGrid<T>::computeAMRResidualLevel");
  CH_assert(m_hasInitBeenCalled[ilev]=='t');

  //m_op[ilev]->setToZero(*(a_resid[l_max]));
  if (l_max != l_base)
    {
      if (ilev == l_max)
        {
          m_op[l_max]->AMRResidualNF(*(a_resid[l_max]), *(a_phi[l_max]),
                                     *(a_phi[l_max-1]), *(a_rhs[l_max]),
                                     a_homogeneousBC);
        }
      else if (ilev == l_base)
        {
          if (l_base == 0)
            {
              m_op[l_base]->AMRResidualNC(*(a_resid[l_base]), *(a_phi[l_base+1]),
                                          *(a_phi[l_base]),  *(a_rhs[l_base]),
                                           a_homogeneousBC, m_op[l_base+1]);
            }
          else
            {
              m_op[l_base]->AMRResidual(*a_resid[l_base], *a_phi[l_base+1], *a_phi[l_base],
                                        *a_phi[l_base-1], *a_rhs[l_base],
                                         a_homogeneousBC, m_op[l_base+1]);
            }
        }
      else
        {
          m_op[ilev]->AMRResidual(*a_resid[ilev], *a_phi[ilev+1], *a_phi[ilev],
                                  *a_phi[ilev-1], *a_rhs[ilev],
                                   a_homogeneousBC, m_op[ilev+1]);
        }
    }
  else
    {
      CH_assert(ilev == l_base);
      if (l_base == 0)
        {
          m_op[l_max]->residual(*a_resid[l_max], *a_phi[l_max], *a_rhs[l_max],a_homogeneousBC);
        }
      else
        {
          m_op[l_max]->AMRResidualNF(*(a_resid[l_max]), *(a_phi[l_max]),
                                     *(a_phi[l_max-1]), *(a_rhs[l_max]),
                                      a_homogeneousBC);
        }
    }

}

template<class T>
void AMRMultiGrid<T>::solve(Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                            int l_max, int l_base, bool a_zeroPhi,
                            bool a_forceHomogeneous)
{
  CH_TIME("AMRMultiGrid::solve");
  init(a_phi, a_rhs, l_max, l_base);
  solveNoInit(a_phi, a_rhs, l_max, l_base, a_zeroPhi, a_forceHomogeneous);
  //puts multigrid depth back where it started so the solver can be reused
  revert(a_phi, a_rhs, l_max, l_base);
}

template<class T>
void AMRMultiGrid<T>::solveNoInit(Vector<T*>& a_phi,
                                  const Vector<T*>& a_rhs,
                                  int l_max, int l_base, bool a_zeroPhi,
                                  bool a_forceHomogeneous)
{
  CH_TIMERS("AMRMultiGrid::solveNoInit");
  Vector<T*> uberResidual(a_rhs.size());
  int lowlim = l_base;
  if (l_base > 0)  // we need an ubercorrection one level lower than l_base
     lowlim--;
  // for (int ilev = l_base; ilev <= l_max; ilev++)
  for (int ilev = lowlim; ilev <= l_max; ilev++)
    {
      uberResidual[ilev] = new T();
      if (ilev >= l_base)
        {
          m_op[ilev]->create(*uberResidual[ilev], *a_rhs[ilev]);
        }
    }
  solveNoInitResid(a_phi, uberResidual, a_rhs, l_max, l_base, a_zeroPhi, a_forceHomogeneous);
  for (int i = lowlim; i <= l_max; i++)
    {
      m_op[i]->clear(*uberResidual[i]);
      delete uberResidual[i];
    }
}

template<class T>
void AMRMultiGrid<T>::solveNoInitResid(Vector<T*>& a_phi,  Vector<T*>& uberResidual,
                                       const Vector<T*>& a_rhs,
                                       int l_max, int l_base, bool a_zeroPhi,
                                       bool a_forceHomogeneous)
{
  CH_TIMERS("AMRMultiGrid::solveNo-InitResid");
  CH_TIMER("AMRMultiGrid::AMRVcycle", vtimer);
  setBottomSolver(l_max, l_base);

  CH_assert(l_base <= l_max);
  CH_assert(a_rhs.size() == a_phi.size());

  //these correspond to the residual and correction
  //that live in AMRSolver
  Vector<T*> uberCorrection(a_rhs.size());

  int lowlim = l_base;
  if (l_base > 0)  // we need an ubercorrection one level lower than l_base
     lowlim--;
  // for (int ilev = l_base; ilev <= l_max; ilev++)
  bool outputIntermediates = false;

  for (int ilev = lowlim; ilev <= l_max; ilev++)
    { CH_TIME("uberCorrection");
      uberCorrection[ilev] = new T();
      m_op[ilev]->create(*uberCorrection[ilev], *a_phi[ilev]);
      if (ilev >= l_base)
        {
          m_op[ilev]->create(*uberResidual[ilev], *a_rhs[ilev]);
          m_op[ilev]->setToZero(*(uberResidual[ilev]));
        }
      m_op[ilev]->setToZero(*(uberCorrection[ilev]));
    }
  //  m_op[0]->dumpStuff(uberResidual, string("initialRes.hdf5"));
  if (a_zeroPhi)
    for (int ilev = l_base; ilev <=l_max; ++ilev)
      {
        m_op[ilev]->setToZero(*(a_phi[ilev]));
      }
  //compute initial residual and initialize internal residual to it

  Real initial_rnorm = 0;
  {
    CH_TIME("Initial AMR Residual");
    initial_rnorm=computeAMRResidual(uberResidual, a_phi, a_rhs, l_max, l_base, a_forceHomogeneous, true);
  }

  if (m_convergenceMetric != 0.)
    {
      initial_rnorm = m_convergenceMetric;
    }

  Real rnorm = initial_rnorm;
  Real norm_last = 2*initial_rnorm;

  /// set bottom solver convergence norm and solver tolerance
  m_bottomSolver->setConvergenceMetrics(initial_rnorm, m_bottomSolverEpsCushion*m_eps);

  int iter=0;
  if (m_verbosity >= 2) ////
    {
      pout() << "    AMRMultiGrid:: iteration = " << iter << ", residual norm = " << rnorm << std::endl;
    }

  bool goNorm = rnorm > m_normThresh;                        //iterate if norm is not small enough
  bool goRedu = rnorm > m_eps*initial_rnorm;                 //iterate if initial norm is not reduced enough
  bool goIter = iter < m_iterMax;                            //iterate if iter < max iteration count
  bool goHang = iter < m_imin || rnorm <(1-m_hang)*norm_last;//iterate if we didn't hang
  bool goMin = iter < m_iterMin ; // iterate if iter < min
  while (goMin || (goIter && goRedu && goHang && goNorm))
    {

      // Report the residuals to the inspectors.
      for (int i = 0; i < m_inspectors.size(); ++i)
        m_inspectors[i]->recordResiduals(uberResidual, l_base, l_max, iter);

      if (outputIntermediates)
        {
          char strresname[100];
          sprintf(strresname, "amrmg.res.iter.%03d", iter);
          string nameres(strresname);
          outputAMR(uberResidual, nameres, l_max, l_base);
        }

      norm_last = rnorm;

      //this generates a correction from the current residual
      CH_START(vtimer);
      AMRVCycle(uberCorrection, uberResidual, l_max, l_max, l_base);
      CH_STOP(vtimer);

      char charname[100];
      sprintf(charname, "resid_iter%03d.%dd.hdf5", iter, SpaceDim);
      string sname(charname);
      //      m_op[0]->dumpAMR(uberResidual, sname);
      // Report the corrections to the inspectors.
      for (int i = 0; i < m_inspectors.size(); ++i)
        m_inspectors[i]->recordCorrections(uberCorrection, l_base, l_max, iter);

      //increment phi by correction and reset correction to zero
      for (int ilev = l_base; ilev <= l_max; ilev++)
        {
          if (outputIntermediates)
            {
              char strcorname[100];
              sprintf(strcorname, "amrmg.phi.iter.%03d", iter);
              string namecor(strcorname);
              outputAMR(a_phi, namecor, l_max, l_base);
            }

          m_op[ilev]->incr(*(a_phi[ilev]), *(uberCorrection[ilev]), 1.0);

          if (outputIntermediates)
            {
              char strcorname[100];
              sprintf(strcorname, "amrmg.cor.iter.%03d", iter);
              string namecor(strcorname);
              outputAMR(uberCorrection, namecor, l_max, l_base);
            }

          m_op[ilev]->setToZero(*(uberCorrection[ilev]));
          // clean const out
        }

      // For solvers with accuracy higher than 2nd order
      //  consistency between levels has to be explicitly enforced.
      // Qinghai Zhang
      if (m_op[0]->orderOfAccuracy()>2)
        {
          for (int ilev=l_max; ilev>l_base; ilev--)
            {
              m_op[ilev]->enforceCFConsistency(*a_phi[ilev-1], *a_phi[ilev]);
            }
        }
      //------end enforcing consistency.

      //recompute residual
      rnorm = computeAMRResidual(uberResidual, a_phi, a_rhs, l_max, l_base, a_forceHomogeneous, true);
      iter++;
      if (m_verbosity >= 3) ////
        {
          pout() << "    AMRMultiGrid:: iteration = " << iter << ", residual norm = " << rnorm;
          if (rnorm > 0.0)
          {
            pout() << ", rate = " << norm_last/rnorm;
          }
          pout() << std::endl;
        }

      goNorm = rnorm > m_normThresh;                        //keep iterating if norm is not small enough
      goRedu = rnorm > m_eps*initial_rnorm;                 //keep iterating if initial norm is not reduced enough
      goIter = iter < m_iterMax;                            //keep iterating if iter < max iteration count
      goHang = iter < m_imin || rnorm <(1-m_hang)*norm_last;//keep iterating if we didn't hang
      goMin = iter < m_iterMin ; // keep iterating if iter < min
    }
  // if ((rnorm > 10.*initial_rnorm) && (rnorm > 10.*m_eps))
  //   {
  //     pout() << "solver seems to have blown up" << endl;
  //     MayDay::Error("kaboom");
  //   }
  m_exitStatus = int(!goRedu) + int(!goIter)*2 + int(!goHang)*4 + int(!goNorm)*8;
  if (m_verbosity >= 2)
    {
      pout() << "    AMRMultiGrid:: iteration = " << iter << ", residual norm = " << rnorm << std::endl;
    }
  if (m_verbosity > 1)
    {
      if (!goIter && goRedu && goNorm) // goRedu=T, goIter=F, goHang=?, goNorm=T
        { // m_exitStatus == 0 + 2 + 0|4 + 0 = 2|6
          pout() << "    AMRMultiGrid:: WARNING: Exit because max iteration count exceeded" << std::endl;
        }
      if (!goHang && goRedu && goNorm) // goRedu=T, goIter=?, goHang=F, goNorm=T
        { // m_exitStatus == 0 + 0|2 + 4 + 0 = 4|6
          pout() << "    AMRMultiGrid:: WARNING: Exit because of solver hang" << std::endl;
        }
      if (m_verbosity > 4)
        {
          pout() << "    AMRMultiGrid:: exitStatus = " << m_exitStatus << std::endl;
        }
    }
  for (int i = lowlim; i <= l_max; i++)
    {
      m_op[i]->clear(*uberCorrection[i]);
      delete uberCorrection[i];
    }
}

template<class T>
void AMRMultiGrid<T>::relaxOnlyHomogeneous(Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                                           int l_max, int l_base)
{
  CH_TIME("AMRMultiGrid::relaxOnly");
  CH_assert(l_max == l_base);
  CH_assert(l_max == 0);
  init(a_phi, a_rhs, l_max, l_base);

  CH_assert(a_rhs.size() == a_phi.size());

  Vector<T*> uberResidual(a_rhs.size());

  for (int ilev = 0; ilev < m_op.size(); ilev++)
    {
      uberResidual[ilev] = new T();
      m_op[ilev]->create(*uberResidual[ilev], *a_rhs[ilev]);
    }

  //compute initial residual and initialize internal residual to it
  m_op[0]->residual(*uberResidual[0], *a_phi[0], *a_rhs[0], true);
  Real initial_rnorm = m_op[0]->norm(*uberResidual[0], 0);
  Real rnorm = initial_rnorm;
  Real norm_last = 2*initial_rnorm;

  int iter=0;
  if (m_verbosity >= 3)
    {
      pout() << "    AMRMultiGrid::relaxOnly iteration = " << iter << ", residual norm = " << rnorm  << std::endl;
    }
  bool goNorm = rnorm > m_normThresh;                        //iterate if norm is not small enough
  bool goRedu = rnorm > m_eps*initial_rnorm;                 //iterate if initial norm is not reduced enough
  bool goIter = iter < m_iterMax;                            //iterate if iter < max iteration count
  bool goHang = iter < m_imin || rnorm <(1-m_hang)*norm_last;//iterate if we didn't hang
  while (goIter && goRedu && goHang && goNorm)
    {
      norm_last = rnorm;
      m_op[0]->relax(*a_phi[0], *a_rhs[0], 1);

      iter++;
      //recompute residual
      m_op[0]->residual(*uberResidual[0], *a_phi[0], *a_rhs[0], true);
      rnorm = m_op[0]->norm(*uberResidual[0], 2);
      if (m_verbosity >= 4)
        {
          pout() << "    AMRMultiGrid::relaxOnly iteration = " << iter << ", residual norm = " << rnorm
                 << ", rate = " << norm_last/rnorm << std::endl;
        }
      goNorm = rnorm > m_normThresh;                        //keep iterating if norm is not small enough
      goRedu = rnorm > m_eps*initial_rnorm;                 //keep iterating if initial norm is not reduced enough
      goIter = iter < m_iterMax;                            //keep iterating if iter < max iteration count
      goHang = iter < m_imin || rnorm <(1-m_hang)*norm_last;//keep iterating if we didn't hang
    }
  m_exitStatus = 0;
  for (int i = 0; i < m_op.size(); i++)
    {
      m_op[i]->clear(*uberResidual[i]);
      delete uberResidual[i];
    }
}

template <class T>
void AMRMultiGrid<T>::clear()
{
  for (int i = 0; i < m_op.size(); i++)
    {
      m_op[i]->clear(*m_correction[i]);
      m_op[i]->clear(*m_residual[i]);
      m_op[i]->clear(*m_resC[i]);

      delete m_correction[i];
      delete m_residual[i];
      delete m_resC[i];
      delete m_op[i];
      delete m_mg[i];

      m_correction[i] = NULL;
      m_residual[i] = NULL;
      m_resC[i] = NULL;
      m_op[i] = NULL;
      m_mg[i] = NULL;
    }
  m_hasInitBeenCalled.resize(0);
}

template <class T>
void AMRMultiGrid<T>::init(const Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                           int l_max, int l_base)
{
  //if (m_hasInitBeenCalled[l_max]=='t' && m_hasInitBeenCalled[l_base]=='t') return;
  CH_TIME("AMRMultiGrid::init");
//   CH_assert(a_phi.size()>=m_op.size());
//   CH_assert(a_rhs.size()>=m_op.size());
  for (int i = l_base; i <= l_max; i++)
    {
      m_hasInitBeenCalled[i] = 't';
      AMRLevelOp<T>& op = *(m_op[i]);

      op.create(*m_correction[i], *a_phi[i]);
      op.create(*m_residual[i],   *a_rhs[i]);
      int r = op.refToCoarser();
      if (i!= l_base)
        {
          int intermediateDepth = 0;

          while (r>2)
            {
              r/=2;
              intermediateDepth++;
            }
          m_mg[i]->m_depth = intermediateDepth;
        }
      m_mg[i]->init(*a_phi[i],     *a_rhs[i]);

      if (i != l_base)
      {
        r = op.refToCoarser();
        op.createCoarsened(*m_resC[i], *a_rhs[i], r);
        op.buildCopier(m_resCopier[i],  *a_rhs[i-1], *m_resC[i]);
        m_reverseCopier[i] = m_resCopier[i];
        m_reverseCopier[i].reverse();
      }
    }
}

template <class T>
void AMRMultiGrid<T>::revert(const Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                             int l_max, int l_base)
{
  CH_TIME("AMRMultiGrid::revert");
  for (int i = l_base; i <= l_max; i++)
    {
      if (i!= l_base)
        {
          m_mg[i]->m_depth = m_mg[i]->m_defaultDepth;
        }
    }
}

template <class T>
void AMRMultiGrid<T>::setBottomSolver(int l_max, int l_base)
{
  //  for (int ilev = 0; ilev <= m_mg.size(); ilev++)
  for (int ilev = l_base; ilev <= l_max; ilev++)
    {
      // only do this if we haven't already set these
      if (!m_solverParamsSet)
        {
          m_mg[ilev]->m_pre = m_pre;
          m_mg[ilev]->m_post = m_post;
          m_mg[ilev]->m_bottom = m_bottom;
        }
      m_mg[ilev]->m_bottomSolver = &m_nosolve;
    }

  m_mg[l_base]->setBottomSolver(m_bottomSolver);
}

template <class T>
void AMRMultiGrid<T>::setBottomSolverEpsCushion(Real a_bottomSolverEpsCushion)
{
  CH_assert((a_bottomSolverEpsCushion >  0.) &&
            (a_bottomSolverEpsCushion <= 1.));
  m_bottomSolverEpsCushion = a_bottomSolverEpsCushion;
}

template <class T>
void AMRMultiGrid<T>::define(const ProblemDomain& a_coarseDomain,
                             AMRLevelOpFactory<T>& a_factory,
                             LinearSolver<T>* a_bottomSolver,
                             int a_maxAMRLevels)
{
  CH_TIME("AMRMultiGrid::define");
  this->clear();
  m_op.resize( a_maxAMRLevels, NULL);
  m_mg.resize( a_maxAMRLevels, NULL);
  m_hasInitBeenCalled.resize(a_maxAMRLevels, 'f');

  m_correction.resize( a_maxAMRLevels, NULL);
  m_residual.  resize( a_maxAMRLevels, NULL);
  m_resC.      resize( a_maxAMRLevels, NULL);
  m_resCopier. resize( a_maxAMRLevels);
  m_reverseCopier.resize(  a_maxAMRLevels);
  m_bottomSolver = a_bottomSolver;

  ProblemDomain current = a_coarseDomain;
  for (int i = 0; i < a_maxAMRLevels; i++)
    {
      m_correction[i] = new T();
      m_residual[i]   = new T();
      m_resC[i]       = new T();
      m_mg[i]= new MultiGrid<T>();
      m_op[i] = a_factory.AMRnewOp(current);
      m_mg[i]->define(a_factory, &m_nosolve, current, m_maxDepth, m_op[i]);

      // Only do this if it will be used (avoiding a reference to invalid
      // and/or unavailable refinement ratios)
      if (i < a_maxAMRLevels-1)
      {
        current.refine(a_factory.refToFiner(current));
      }
    }
}

template<class T>
void AMRMultiGrid<T>::AMRVCycle(Vector<T*>& a_uberCorrection,
                                Vector<T*>& a_uberResidual,
                                int ilev, int l_max, int l_base)
{
  if (ilev == l_max)
    {
      for (int level = l_base; level <= l_max; level++)
        {
          m_op[level]->assignLocal(*m_residual[level], *a_uberResidual[level]);
          m_op[level]->setToZero(*m_correction[level]);
        }
    }

  if (l_max == l_base)
    {
      CH_assert(ilev == l_base);
      m_mg[l_base]->oneCycle(*(a_uberCorrection[ilev]), *(a_uberResidual[ilev]));
    }
  else if (ilev == l_base)
    {
      m_mg[l_base]->oneCycle(*(m_correction[ilev]), *(m_residual[ilev]));
      m_op[ilev]->incr(*(a_uberCorrection[ilev]), *(m_correction[ilev]), 1.0);
    }
  else
    {
      //============= Downsweep ========================

      this->relax(*(m_correction[ilev]), *(m_residual[ilev]), ilev, m_pre);
      m_op[ilev]->incr(*(a_uberCorrection[ilev]), *(m_correction[ilev]), 1.0);

      // Set next coarser level correction to zero
      m_op[ilev-1]->setToZero(*(m_correction[ilev-1]));

      // Recompute residual on next coarser level
      //  for the valid region NOT covered by this level.
      computeAMRResidualLevel(m_residual,
                              a_uberCorrection,
                              a_uberResidual,
                              l_max, l_base, ilev-1,
                              true);

      // Compute the restriction of the residual to the coarser level resC.
      m_op[ilev]->AMRRestrictS(*(m_resC[ilev]),
                               *(m_residual[ilev]),
                               *(m_correction[ilev]),
                               *(m_correction[ilev-1]),
                               *(a_uberCorrection[ilev]));

      // Overwrite residual on the valid region of the next coarser level
      //  with coarsened residual from this level
      m_op[ilev-1]->assignCopier(*m_residual[ilev-1], *(m_resC[ilev]), m_resCopier[ilev]);

      //============finish Compute residual for the next coarser level======

      for (int img = 0; img < m_numMG; img++)
        {
          AMRVCycle(a_uberCorrection, a_uberResidual, ilev-1, l_max, l_base);
        }

      //================= Upsweep ======================
      //increment the correction with coarser version
      //m_op[ilev]->AMRProlong(*(m_correction[ilev]), *(m_correction[ilev-1]));
      m_op[ilev]->AMRProlongS(*(m_correction[ilev]), *(m_correction[ilev-1]),
                              *m_resC[ilev], m_reverseCopier[ilev]);
      //recompute residual
      m_op[ilev]->AMRUpdateResidual(*(m_residual[ilev]), *(m_correction[ilev]), *(m_correction[ilev-1]));

      //compute correction to the correction
      T& dCorr = *(a_uberCorrection[ilev]); // user uberCorrection as holder for correction to correction
      m_op[ilev]->setToZero(dCorr);
      this->relax(dCorr, *(m_residual[ilev]), ilev, m_post);

      //correct the correction with the correction to the correction
      m_op[ilev]->incr(*(m_correction[ilev]), dCorr, 1.0);

      m_op[ilev]->assignLocal(*(a_uberCorrection[ilev]), *(m_correction[ilev]));
    }
}

enum OLD_FASMG_type {FULL=0,VCYCLE=1,FCYCLE=2};

///
/**
   Class to solve elliptic equations using the FAS multigrid
 */
template <class T>
class AMRFASMultiGrid : public AMRMultiGrid<T>
{
public:
  AMRFASMultiGrid();
  virtual ~AMRFASMultiGrid();

  virtual void define(const ProblemDomain& a_coarseDomain,
                      AMRLevelOpFactory<T>& a_factory,
                      LinearSolver<T>* a_bottomSolver,
                      int a_numLevels);

  virtual void solveNoInit(Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                           int l_max, int l_base, bool a_zeroPhi=true,
                           bool forceHomogeneous = false);

  void setCycleType( OLD_FASMG_type a_type )
  {
    m_type = a_type;
  }
  void setAvoidNorms( bool b = true )
  {
    m_avoid_norms = b;
  }
  void setNumVcycles( int n )
  {
    m_numVcycles = n;
  }

private:
  virtual void FMG(Vector<T*>& a_phi,
                   const Vector<T*>& a_rhs,
                   int l, int l_max, int l_base);
  
  virtual void VCycle(Vector<T*>& a_phi,
                      const Vector<T*>& a_rhs,
                      int l, int l_max, int l_base);

  virtual void init(const Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                    int l_max, int l_base);

  void clear_private();

  // data
  OLD_FASMG_type m_type;
  bool m_avoid_norms;  // flag to avoid norms and residuals (for convergence checking)
  int m_numVcycles;
  Vector<Copier> m_HOCopier;
  Vector<Copier> m_HOCornerCopier;
  ProblemDomain m_coarseDomain; // need to cache this
};

//*******************************************************
// AMRFASMultigrid Implementation
//*******************************************************

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::AMRFASMultiGrid
////////////////////////////////////////////////////////////////////////
template <class T>
AMRFASMultiGrid<T>::AMRFASMultiGrid() 
  : AMRMultiGrid<T>()
{
  m_numVcycles = 1;
  m_type = VCYCLE;
  m_avoid_norms = false;
}

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::~AMRFASMultiGrid
////////////////////////////////////////////////////////////////////////
template <class T>
AMRFASMultiGrid<T>::~AMRFASMultiGrid()
{
  CH_TIME("~AMRFASMultiGrid");
  clear_private();
}

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::define
////////////////////////////////////////////////////////////////////////
template <class T>
void AMRFASMultiGrid<T>::define( const ProblemDomain& a_coarseDomain,
                                 AMRLevelOpFactory<T>& a_factory,
                                 LinearSolver<T>* a_bottomSolver,
                                 int a_maxAMRLevels )
{
  CH_TIME("AMRFASMultiGrid::define");
  
  AMRMultiGrid<T>::define( a_coarseDomain,a_factory,a_bottomSolver,a_maxAMRLevels);
  
  m_HOCopier.resize( a_maxAMRLevels );
  m_HOCornerCopier.resize( a_maxAMRLevels );
  m_coarseDomain = a_coarseDomain;
}

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::init
////////////////////////////////////////////////////////////////////////
template <class T>
void AMRFASMultiGrid<T>::init( const Vector<T*>& a_phi, const Vector<T*>& a_rhs,
                               int l_max, int l_base )
{
  CH_TIME("AMRFASMultiGrid::init");

  AMRMultiGrid<T>::init( a_phi, a_rhs, l_max,l_base );

  // set new copiers for high order interp
  ProblemDomain dom = m_coarseDomain;
  for (int i = l_base+1; i <= l_max; i++)
    {
      AMRLevelOp<T>& op = *(this->m_op[i]);
      int r = op.refToCoarser();
      m_HOCopier[i].define( a_phi[i-1]->disjointBoxLayout(), this->m_resC[i]->disjointBoxLayout(), 
                            a_phi[i-1]->ghostVect() );
      dom = dom.refine(r);
      //m_HOCornerCopier[i].define( this->m_resC[i]->disjointBoxLayout(), this->m_resC[i]->disjointBoxLayout(), dom, this->m_resC[i]->ghostVect(), true ); -- smart corner copier -- needed for AMR
      //m_HOCornerCopier[i].define( this->m_resC[i]->disjointBoxLayout(), this->m_resC[i]->disjointBoxLayout(), this->m_resC[i]->ghostVect() ); -- dumb copier - not needed w/o AMR
    }
}

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::clear_private
////////////////////////////////////////////////////////////////////////
template <class T>
void AMRFASMultiGrid<T>::clear_private()
{
  CH_TIME("AMRFASMultiGrid::clear_private");

  for (int i = 0; i < this->m_op.size(); i++)
    {
    }
}

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::solveNoInit
////////////////////////////////////////////////////////////////////////
template<class T>
void AMRFASMultiGrid<T>::solveNoInit( Vector<T*>& a_phi,
                                      const Vector<T*>& a_rhs,
                                      int l_max, int l_base, bool a_zeroPhi,
                                      bool a_forceHomogeneous)
{
  CH_TIMERS("AMRFASMultiGrid::solveNoInit");
  CH_TIMER("AMRFASMultiGrid::cycle", vtimer);
  CH_assert(!a_forceHomogeneous);
  bool outputIntermediates = false;

  this->setBottomSolver( l_max, l_base );

  CH_assert(l_base <= l_max);
  CH_assert(a_rhs.size() == a_phi.size());

  if (a_zeroPhi)
    for (int ilev = l_base; ilev <=l_max; ++ilev)
      {
        this->m_op[ilev]->setToZero(*(a_phi[ilev]));
      }

  Real initial_rnorm = 0;
  if ( !m_avoid_norms )
  {
    CH_TIME("Initial AMR Residual");
    initial_rnorm = this->computeAMRResidual( a_phi, a_rhs, l_max, l_base );
  }

  if (this->m_convergenceMetric != 0.)
    {
      initial_rnorm = this->m_convergenceMetric;
    }

  // set bottom solver convergence norm and solver tolerance
  this->m_bottomSolver->setConvergenceMetrics(initial_rnorm, this->m_bottomSolverEpsCushion*this->m_eps);

  Real rnorm = initial_rnorm;
  Real norm_last = 2*initial_rnorm;

  int iter=0;
  if (this->m_verbosity >= 2 && !m_avoid_norms )
    {
      pout() << "    AMRFASMultiGrid:: iteration = " << iter << ", residual norm = " << rnorm << std::endl;
    }

  bool goNorm = rnorm > this->m_normThresh || m_avoid_norms;       //iterate if norm is not small enough
  bool goRedu = rnorm > this->m_eps*initial_rnorm || m_avoid_norms;//iterate if initial norm is not reduced enough
  bool goIter = iter < this->m_iterMax;                            //iterate if iter < max iteration count
  bool goHang = iter < this->m_imin || rnorm <(1-this->m_hang)*norm_last; //iterate if we didn't hang
  bool goMin = iter < this->m_iterMin ; // iterate if iter < min
  while (goMin || (goIter && goRedu && goHang && goNorm))
    {
      norm_last = rnorm;

      //this generates a correction from the current residual
      CH_START(vtimer);
      if ( m_type == FULL )
        {
          FMG( a_phi, a_rhs, l_max, l_max, l_base);
        }
      else if ( m_type == VCYCLE )
        {
          this->m_op[l_max]->assignLocal( *(this->m_residual[l_max]), *(a_rhs[l_max]) );
          VCycle( a_phi, a_rhs, l_max, l_max, l_base);
        }
      else 
        {
          MayDay::Error("unknown FAS type");
        }
      CH_STOP(vtimer);
      //increment phi by correction and reset correction to zero
      for (int ilev = l_base; ilev <= l_max; ilev++)
        {
          if (outputIntermediates)
            {
              char strcorname[100];
              sprintf(strcorname, "amrmg.phi.iter.%03d", iter);
              string namecor(strcorname);
              this->outputAMR(a_phi, namecor, l_max, l_base);
            }
        }

      // For solvers with accuracy higher than 2nd order
      //  consistency between levels has to be explicitly enforced.
      if (this->m_op[0]->orderOfAccuracy()>2)
        {
          for (int ilev=l_max; ilev>l_base; ilev--)
            {
              this->m_op[ilev]->enforceCFConsistency(*a_phi[ilev-1], *a_phi[ilev]);
            }
        }
      
      // recompute residual
      iter++;
      if ( this->m_verbosity >= 2 && !m_avoid_norms )
        {
          rnorm = this->computeAMRResidual( a_phi, a_rhs, l_max, l_base );
          
          pout() << "    AMRFASMultiGrid:: iteration = " << iter << ", residual norm = " << rnorm;
          if (rnorm > 0.0)
            {
              pout() << ", rate = " << norm_last/rnorm;
            }
          pout() << std::endl;

          goNorm = rnorm > this->m_normThresh;                        //keep iterating if norm is not small enough
          goRedu = rnorm > this->m_eps*initial_rnorm;                 //keep iterating if initial norm is not reduced enough
          goHang = iter < this->m_imin || rnorm <(1-this->m_hang)*norm_last;//keep iterating if we didn't hang
        }
      goIter = iter < this->m_iterMax;                            //keep iterating if iter < max iteration count
      goMin = iter < this->m_iterMin ; // keep iterating if iter < min
    }

  this->m_exitStatus = int(!goRedu) + int(!goIter)*2 + int(!goHang)*4 + int(!goNorm)*8;
  if (this->m_verbosity >= 2 && !m_avoid_norms)
    {
      pout() << "    AMRFASMultiGrid:: iteration = " << iter << ", residual norm = " << rnorm << std::endl;
    }
  if (this->m_verbosity > 1)
    {
      if (!goIter && goRedu && goNorm) // goRedu=T, goIter=F, goHang=?, goNorm=T
        { // m_exitStatus == 0 + 2 + 0|4 + 0 = 2|6
          pout() << "    AMRMultiGrid:: WARNING: Exit because max iteration count exceeded" << std::endl;
        }
      if (!goHang && goRedu && goNorm) // goRedu=T, goIter=?, goHang=F, goNorm=T
        { // m_exitStatus == 0 + 0|2 + 4 + 0 = 4|6
          pout() << "    AMRMultiGrid:: WARNING: Exit because of solver hang" << std::endl;
        }
      if (this->m_verbosity > 4)
        {
          pout() << "    AMRMultiGrid:: exitStatus = " << this->m_exitStatus << std::endl;
        }
    }
}

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::FMG
////////////////////////////////////////////////////////////////////////
template<class T>
void AMRFASMultiGrid<T>::FMG( Vector<T*>& a_phi,
                              const Vector<T*>& a_rhs,
                              int i, int l_max, int l_base)
{
  // coarse grid
  this->m_op[l_base]->assignLocal( *(this->m_residual[l_base]), *(a_rhs[l_base]) );
  this->m_mg[l_base]->oneCycle( *(a_phi[l_base]), *(this->m_residual[l_base]) );

  for ( int ilev = l_base+1 ; ilev <= l_max ; ilev++ )
    {
      // FMG interpolation -- high order
      this->m_op[ilev]->AMRProlongS_2( *(a_phi[ilev]), *(a_phi[ilev-1]),
                                       *(this->m_resC[ilev]), this->m_HOCopier[ilev],
                                       this->m_HOCornerCopier[ilev], this->m_op[ilev-1] );
      // V-cycle
      this->m_op[ilev]->assignLocal( *(this->m_residual[ilev]), *(a_rhs[ilev]) );
      for (int i=0;i<m_numVcycles;i++)
        {
          VCycle( a_phi, a_rhs, ilev, l_max, l_base );
        }
    }
}

////////////////////////////////////////////////////////////////////////
// AMRFASMultiGrid<T>::VCycle
//   'm_residual' is really the full RHS ala FAS.
//   'm_correction' is just used as a temp
////////////////////////////////////////////////////////////////////////
//#define FAS_PRINT
template<class T>
void AMRFASMultiGrid<T>::VCycle( Vector<T*>& a_phi,
                                 const Vector<T*>& a_rhs_dummy,
                                 int ilev, int l_max, int l_base )
{
  if (ilev == l_base)
    {
      CH_TIME("AMRFASMultiGrid::VCycle:coarse_grid_solver");
      // exact solver

#ifdef FAS_PRINT
      this->m_op[ilev]->write(a_phi[ilev],"zh_coarsest_phi_0.hdf5");
      this->m_op[ilev]->write(this->m_residual[ilev],"zh_coarsest_rhs.hdf5");
#endif

      this->m_mg[l_base]->oneCycle( *(a_phi[l_base]), *(this->m_residual[l_base]) );

#ifdef FAS_PRINT
      this->m_op[ilev]->write(a_phi[ilev],"zh_coarsest_phi_1.hdf5");
#endif
    }
  else
    {
      //============= Downsweep ========================

#ifdef FAS_PRINT
this->m_op[ilev]->write(this->m_residual[ilev],"za_fine_res_0.hdf5");
this->m_op[ilev]->write(a_phi[ilev],"zb_fine_phi_0.hdf5");
#endif

      this->m_op[ilev]->relax( *(a_phi[ilev]), *(this->m_residual[ilev]), this->m_pre );

#ifdef FAS_PRINT
this->m_op[ilev]->write(this->m_residual[ilev],"za_fine_res_1.hdf5");
this->m_op[ilev]->write(a_phi[ilev],"zb_fine_phi_1.hdf5");
this->m_op[ilev]->residual( *(this->m_correction[ilev]),  *a_phi[ilev], *(this->m_residual[ilev]), false);
this->m_op[ilev]->write(this->m_correction[ilev],"zz_fine_res_0.hdf5");
#endif

      // Compute the restriction of the solution to the coarser level phiC -- R(u_f) 
      this->m_op[ilev]->AMRRestrictS(*(this->m_resC[ilev]),      // output
                                     *(a_phi[ilev]),             // input
                                     *(a_phi[ilev]),             // dummy
                                     *(a_phi[ilev-1]),           // dummy
                                     *(this->m_correction[ilev]),// scratch
                                     true );                     // skip residuals

#ifdef FAS_PRINT
this->m_op[ilev]->write(a_phi[ilev],"zb_fine_phi_2.hdf5");
this->m_op[ilev-1]->write(this->m_resC[ilev],"zb_coarse_phi_1.hdf5");
#endif
      // Overwrite R(u_f) on the valid region of the next coarser level a_phi[ilev-1]
      //this->m_op[ilev-1]->assignCopier( *(a_phi[ilev-1]), *(m_phiC[ilev]), m_phiCopier[ilev] );
      this->m_resC[ilev]->copyTo(this->m_resC[ilev]->interval(), *(a_phi[ilev-1]), a_phi[ilev-1]->interval(), this->m_resCopier[ilev]);

      // Recompute residual on next coarser level
      //  for the valid region NOT covered by this level.
      // Compute the restriction of the _residual_ to the coarser level resC -- R(f - L_f(u_f))
      this->m_op[ilev]->AMRRestrictS(*(this->m_resC[ilev]),      // output
                                     *(this->m_residual[ilev]),  // const input
                                     *(a_phi[ilev]),             // const but C-F interp modified
                                     *(a_phi[ilev-1]),           // coarse phi, for C-F interp.
                                     *(this->m_correction[ilev])); // scratch

#ifdef FAS_PRINT
this->m_op[ilev]->write(this->m_resC[ilev],"zd_resC.hdf5");
this->m_op[ilev]->write(a_phi[ilev-1],"ze_coarse_phi_2.hdf5");
this->m_op[ilev]->write(this->m_residual[ilev],"za_fine_res_2.hdf5");
#endif

      // Overwrite residual on the valid region of the next coarser level
      //  with coarsened residual from this level
      //this->m_op[ilev-1]->assignCopier( *(this->m_residual[ilev-1]), *(this->m_resC[ilev]), this->m_resCopier[ilev]);
      this->m_resC[ilev]->copyTo(this->m_resC[ilev]->interval(), *(this->m_residual[ilev-1]), this->m_residual[ilev-1]->interval(), this->m_resCopier[ilev] );

#ifdef FAS_PRINT
this->m_op[ilev-1]->write(this->m_residual[ilev-1],"zf_coarse_res_1.hdf5");
#endif

      // put tau correction in m_correction[ilev-1] -- R( f - L_f(u_f) ) + L_c( R(u_f) )
      this->m_op[ilev-1]->AMROperatorNC( *(this->m_correction[ilev-1]), *(a_phi[ilev]), *(a_phi[ilev-1]), true, this->m_op[ilev]);

#ifdef FAS_PRINT
this->m_op[ilev-1]->write(this->m_correction[ilev-1],"zf_coarse_Lu.hdf5");
this->m_op[ilev]->write(a_phi[ilev],"zb_fine_phi_3.hdf5");
#endif

      this->m_op[ilev-1]->axby( *(this->m_residual[ilev-1]), *(this->m_residual[ilev-1]), *(this->m_correction[ilev-1]), 1.0, 1.0);

#ifdef FAS_PRINT
this->m_op[ilev-1]->write(this->m_residual[ilev-1],"zf_coarse_res_3.hdf5");
#endif

      {
        // store correction in R_u_f -- R(u_f)
        T R_u_f;
        this->m_op[ilev-1]->create( R_u_f, *a_phi[ilev-1]);     
        this->m_op[ilev-1]->assignLocal( R_u_f, *(a_phi[ilev-1]) ); 
        //============finish Compute residual for the next coarser level======
        for (int img = 0; img < this->m_numMG; img++)
          {
            VCycle( a_phi, a_rhs_dummy, ilev-1, l_max, l_base );
          }

#ifdef FAS_PRINT
this->m_op[ilev-1]->write(a_phi[ilev-1],"zi_coarse_phi_2.hdf5");
#endif

        // subtract off initial solution to get an increment (a_phi[ilev-1] is an increment)
        this->m_op[ilev-1]->axby( *a_phi[ilev-1], *a_phi[ilev-1], R_u_f, 1.0, -1.0 );   
        this->m_op[ilev-1]->clear( R_u_f );
      }

#ifdef FAS_PRINT
this->m_op[ilev-1]->residual( *(this->m_correction[ilev-1]),  *a_phi[ilev-1], *(this->m_residual[ilev-1]), true);
this->m_op[ilev-1]->write(this->m_correction[ilev-1],"zz_coarse_res.hdf5");
#endif

      //================= Upsweep ======================
      // increment the correction with coarser version
      this->m_op[ilev]->AMRProlongS_2( *(a_phi[ilev]), *(a_phi[ilev-1]),
                                       *(this->m_resC[ilev]), this->m_HOCopier[ilev],
                                       this->m_HOCornerCopier[ilev], this->m_op[ilev-1] );

#ifdef FAS_PRINT 
this->m_op[ilev]->write(a_phi[ilev],"zl_fine_phi_1.hdf5");
this->m_op[ilev]->residual( *(this->m_correction[ilev]),  *a_phi[ilev], *(this->m_residual[ilev]), false);
this->m_op[ilev]->write(this->m_correction[ilev],"zz_fine_res_1.hdf5");
#endif

      this->m_op[ilev]->relax( *(a_phi[ilev]), *(this->m_residual[ilev]), this->m_post );

#ifdef FAS_PRINT
this->m_op[ilev]->write(a_phi[ilev],"zl_fine_phi_2.hdf5");
this->m_op[ilev]->write(this->m_residual[ilev],"za_fine_res_3.hdf5");
this->m_op[ilev]->write(a_rhs_dummy[ilev],"za_rhs.hdf5");
this->m_op[ilev]->residual( *(this->m_correction[ilev]),  *a_phi[ilev], *(this->m_residual[ilev]), false);
this->m_op[ilev]->write(this->m_correction[ilev],"zz_fine_res_2.hdf5");
#endif
this->m_op[ilev]->write(a_phi[ilev],"exact.hdf5"); // debug
    }
}

#include "NamespaceFooter.H"

#endif
