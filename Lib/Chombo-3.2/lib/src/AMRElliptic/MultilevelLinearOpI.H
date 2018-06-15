#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTILEVELLINEAROPI_H_
#define _MULTILEVELLINEAROPI_H_

#include "AMRIO.H"

#include "NamespaceHeader.H"

// set a couple of defaults here (do before define so that
// they can be overridden before calling define)
template<class T>
MultilevelLinearOp<T>::MultilevelLinearOp()
{
  // set some default values
  m_use_multigrid_preconditioner = true;
  m_num_mg_iterations = 1;
  m_num_mg_smooth = 4;
  // default is to coarsen all the way down
  m_preCondSolverDepth = -1;
  m_precondBottomSolverPtr = NULL;
}

/// define function
template<class T>
void
MultilevelLinearOp<T>::define(const Vector<DisjointBoxLayout>& a_vectGrids,
                              const Vector<int>& a_refRatios,
                              const Vector<ProblemDomain>& a_domains,
                              const Vector<RealVect>& a_vectDx,
                              RefCountedPtr<AMRLevelOpFactory<LevelData<T> > >& a_opFactory,
                              int a_lBase)
{
  CH_TIME("MultilevelLinearOp::define");

  int numLevels = a_vectGrids.size();

  // couple of sanity checks...
  CH_assert (a_lBase >= 0);
  CH_assert (a_lBase < numLevels);
  CH_assert ( a_domains.size() == numLevels);
  // since you technically only need numLevels -1 refinement ratios...
  CH_assert (a_refRatios.size() >= numLevels -1);

  m_lBase = a_lBase;

  m_vectGrids = a_vectGrids;
  m_refRatios = a_refRatios;
  m_domains = a_domains;
  m_vectDx = a_vectDx;
  m_vectOperators.resize(numLevels);

  // need to at least define the levels to lBase -1 (for coarser BC's)
  int coarsestLevel = Max(0, m_lBase-1);
  for (int level=coarsestLevel; level<numLevels; level++)
    {
      RefCountedPtr<AMRLevelOp<LevelData<T> > > op = RefCountedPtr<AMRLevelOp<LevelData<T> > >(a_opFactory->AMRnewOp(m_domains[level]));
      m_vectOperators[level] = op;
    }

  // finally, define AMRMultigrid preconditioner if required
  if (m_use_multigrid_preconditioner)
    {
      // preconditioner requires a bottom smoother
      BiCGStabSolver<LevelData<T> >* newPtr = new BiCGStabSolver<LevelData<T> >;

      m_precondBottomSolverPtr = newPtr;
      int bottomSolverVerbosity = 1;

      newPtr->m_verbosity = bottomSolverVerbosity;

      if (m_preCondSolverDepth >= 0)
        {
          m_preCondSolver.m_maxDepth = m_preCondSolverDepth;
        }
      m_preCondSolver.define(m_domains[0],
                             *a_opFactory,
                             m_precondBottomSolverPtr,
                             numLevels);

      // most of these have no real meaning since we're only doina a
      // single V-cycle, rather than a full solve
      Real solverEps = 1.0e-7;
      int maxIterations = 1;
      Real hang = 1.0e-7;
      Real normThresh = 1.0e-30;

      int num_mg = 1;
      m_preCondSolver.setSolverParameters(m_num_mg_smooth,
                                          m_num_mg_smooth,
                                          m_num_mg_smooth,
                                          num_mg,
                                          maxIterations,
                                          solverEps,
                                          hang,
                                          normThresh);

      // set preconditioner solver to be _really_ quiet...
      m_preCondSolver.m_verbosity = 1;

      // AMRMultiGrid::init has not yet been called
      // (will call later, during actual solve)
      m_isPrecondSolverInitialized = false;

    }
}

  ///
template<class T>
MultilevelLinearOp<T>::~MultilevelLinearOp()
{
  CH_TIME("MultilevelLinearOp::~MultilevelLinearOp");

  if (m_precondBottomSolverPtr != NULL)
    {
      delete m_precondBottomSolverPtr;
      m_precondBottomSolverPtr = NULL;
    }
}

///
/**
   Say you are  solving L(phi) = rhs.   Make a_lhs = L(a_phi) - a_rhs.
   If a_homogeneous is true, evaluate the operator using homogeneous
   boundary conditions.
*/
template<class T>
void
MultilevelLinearOp<T>::residual(Vector<LevelData<T>* >& a_lhs,
                                const Vector<LevelData<T>* >& a_phi,
                                const Vector<LevelData<T>* >& a_rhs,
                                bool a_homogeneous )
{
  CH_TIME("MultilevelLinearOp::residual");

  // do this using the operators in AMRLevelOp
  LevelData<T> dummyT;

  int numLevels = a_lhs.size();
  CH_assert (a_phi.size() == numLevels);
  CH_assert (a_rhs.size() == numLevels);
  CH_assert (m_vectOperators.size() >= numLevels);

  int maxLevel = numLevels -1;
  for (int level = m_lBase; level < numLevels; level++)
    {
      const LevelData<T>& thisPhi = *(a_phi[level]);
      LevelData<T>& thisResid = *(a_lhs[level]);
      const LevelData<T>& thisRHS = *(a_rhs[level]);
      // theres gotta be a cleaner way to do this...
      if (level >0)
        {
          const LevelData<T>& phiCrse = *a_phi[level-1];

          if (level < maxLevel)
            {
              const LevelData<T>& phiFine = *a_phi[level+1];
              m_vectOperators[level]->AMRResidual(thisResid,
                                                  phiFine,
                                                  thisPhi,
                                                  phiCrse,
                                                  thisRHS,
                                                  a_homogeneous,
                                                  (m_vectOperators[level+1].operator->()));
            }
          else
            {
              // no finer level exists
              m_vectOperators[level]->AMRResidualNF(thisResid,
                                                    thisPhi,
                                                    phiCrse,
                                                    thisRHS,
                                                    a_homogeneous);
            } // end if no finer level

        } // end if a coarser level exists
      else
        {
          // we're on level 0
          if (level < maxLevel)
            {
              const LevelData<T>& phiFine = *a_phi[level+1];
              m_vectOperators[level]->AMRResidualNC(thisResid,
                                                    phiFine,
                                                    thisPhi,
                                                    thisRHS,
                                                    a_homogeneous,
                                                    (m_vectOperators[level+1].operator->()));
            }
          else
            {
              // no finer level exists
              m_vectOperators[level]->residual(thisResid,
                                               thisPhi,
                                               thisRHS,
                                               a_homogeneous);
            } // end if no finer level
        } // end if on level 0
    } // end loop over levels
}

///
/**
   Given the current state of the residual the correction, apply
   your preconditioner to a_cor.
*/
template<class T>
void
MultilevelLinearOp<T>::preCond(Vector<LevelData<T>* >& a_cor,
                               const Vector<LevelData<T>* >& a_residual)

{
  CH_TIME("MultilevelLinearOp::preCond");

  // sanity checks
  int numLevels = a_cor.size();
  CH_assert (a_residual.size() == numLevels);
  CH_assert (m_vectOperators.size() <= numLevels);

  // use AMR multigrid for a preconditioner
  if (m_use_multigrid_preconditioner)
    {
      CH_TIME("MultilevelLinearOp::preCond::Multigrid");

      Vector<LevelData<T>* > localCorr;
      Vector<LevelData<T>* > localResid;
      create(localCorr, a_cor);
      create(localResid, a_residual);

      // this is kinda silly, but I need to get rid of RefCountedPtr
      // part of things
      Vector<LevelData<T>* > tempCorr(numLevels, NULL);
      Vector<LevelData<T>* > tempResid(numLevels, NULL);

      {
        CH_TIME("MultilevelLinearOp::preCond::Multigrid::NewStorage");

        // in place. Otherwise, need to allocate new storage for them
        for (int lev=0; lev<numLevels; lev++)
          {
            tempCorr[lev] = localCorr[lev];
            tempResid[lev] = const_cast<LevelData<T>*>(localResid[lev]);
          }
      }

      int finestLevel = numLevels -1;
      {
        CH_TIME("MultilevelLinearOp::preCond::Multigrid::Initialize");

        if (!m_isPrecondSolverInitialized)
          {
            // because we're not calling solve, need to call init directly
            m_preCondSolver.init(tempCorr, tempResid, finestLevel,
                                 m_lBase);
            //  also need to set bottom solver
            m_preCondSolver.setBottomSolver(finestLevel, m_lBase);
            m_isPrecondSolverInitialized = true;
          }
      }

      for (int iter = 0; iter<m_num_mg_iterations; iter++)
        {
          CH_TIME("MultilevelLinearOp::preCond::Multigrid::Iterate");

          // use AMRVCycle instead of solve because we're already
          // in residual-correction form, so we want to use
          // homogeneous form of boundary conditions.

          // however, AMRVCycle requires that the initial correction be zero
          // to do this properly, recompute residual before calling V-Cycle

          bool homogeneous = true;
          residual(localResid, a_cor, a_residual, homogeneous);
          setToZero(localCorr);

          m_preCondSolver.AMRVCycle(tempCorr, tempResid,
                                    finestLevel, finestLevel,
                                    m_lBase);

          // now increment a_cor with local correction
          incr(a_cor, localCorr, 1.0);
        }

      {
        CH_TIME("MultilevelLinearOp::preCond::Multigrid::Clear");

        clear(localCorr);
        clear(localResid);
      }
    }
  else
    {
      CH_TIME("MultilevelLinearOp::preCond::LevelBased");

      // just use whatever the AMRLevelOp provides for levels finer than lBase
      for (int level = m_lBase; level < numLevels; level++)
        {
          m_vectOperators[level]->preCond(*a_cor[level], *a_residual[level]);
        }
    }
}

///
/**
   In the context of solving L(phi) = rhs, set a_lhs = L(a_phi).
   If a_homogeneous is true,
   evaluate the operator using homogeneous boundary conditions.
*/
template<class T>
void
MultilevelLinearOp<T>::applyOp(Vector<LevelData<T>* >& a_lhs,
                               const Vector<LevelData<T>* >& a_phi,
                               bool a_homogeneous)
{
  CH_TIME("MultilevelLinearOp::applyOp");

  int numLevels = a_phi.size();
  int maxLevel = numLevels-1;
  LevelData<T> dummyLDF;
  for (int level = m_lBase; level<numLevels; level++)
    {
      if (level == 0)
        {
          if (level < maxLevel)
            {
              m_vectOperators[level]->AMROperatorNC(*a_lhs[level],
                                                    *a_phi[level+1],
                                                    *a_phi[level],
                                                    a_homogeneous,
                                                    m_vectOperators[level+1].operator->());
            }
          else
            {
              // this is the single-level case...
              m_vectOperators[level]->applyOp(*a_lhs[level],
                                              *a_phi[level],
                                              a_homogeneous);
            }
        } // end if level is 0
      else
        {
          // coarser level exists
          if (level < maxLevel)
            {
              m_vectOperators[level]->AMROperator(*a_lhs[level],
                                                  *a_phi[level+1],
                                                  *a_phi[level],
                                                  *a_phi[level-1],
                                                  a_homogeneous,
                                                  m_vectOperators[level+1].operator->());
            }
          else
            {
              m_vectOperators[level]->AMROperatorNF(*a_lhs[level],
                                                    *a_phi[level],
                                                    *a_phi[level-1],
                                                    a_homogeneous);
            }
        } // end if coarser level exists
    } // end loop over levels

}

///
/**
   Create data holder a_lhs that mirrors a_rhs.   You do not need
   to copy the data of a_rhs, just  make a holder the same size.
*/
template<class T>
void
MultilevelLinearOp<T>::create(Vector<LevelData<T>* >& a_lhs,
                              const Vector<LevelData<T>* >& a_rhs)
{
  CH_TIME("MultilevelLinearOp::create");

  a_lhs.resize(a_rhs.size());
  // need to create an lBase-1 just in case it's needed for C/F BCs
  int coarsestLevel = Max(0, m_lBase-1);
  for (int level =coarsestLevel; level < a_rhs.size(); level++)
    {
      a_lhs[level] = new LevelData<T>;
      m_vectOperators[level]->create(*a_lhs[level],
                                     *a_rhs[level]);
    }
}

///
/**
   Clean up data holder before it goes out of scope.  This is necessary
   because create calls new.
*/
template<class T>
void
MultilevelLinearOp<T>::clear(Vector<LevelData<T>* >& a_lhs)
{
  CH_TIME("MultilevelLinearOp::clear");

  for (int level =0; level < a_lhs.size(); level++)
    {
      if (a_lhs[level] != NULL)
        {
          delete a_lhs[level];
          a_lhs[level] = NULL;
        }
    }
}

///
/**
   Set a_lhs  equal to a_rhs.
*/
template<class T>
void
MultilevelLinearOp<T>::assign(Vector<LevelData<T>* >& a_lhs,
                              const Vector<LevelData<T>* >& a_rhs)
{
  CH_TIME("MultilevelLinearOp::assign");

  CH_assert (a_lhs.size() == a_rhs.size());
  // only do this for lBase and finer levels
  for (int level = m_lBase; level < a_lhs.size(); level++)
    {
      if (!(a_lhs[level] ==NULL) && !(a_rhs[level] == NULL))
        {
          m_vectOperators[level]->assign(*a_lhs[level], *a_rhs[level]);
        }
    }
}

///
/**
   Compute and return the dot product of a_1 and a_2.   In most
   contexts, this means return the sum over all data points of a_1*a_2.
*/
template<class T>
Real
MultilevelLinearOp<T>::dotProduct(const Vector<LevelData<T>* >& a_1,
                                  const Vector<LevelData<T>* >& a_2)
{
  CH_TIME("MultilevelLinearOp::dotProduct");

  // want to do this in an AMR way, so need to generate a temporary
  Vector<LevelData<T>* > temp1, temp2;
  create(temp1, a_1);
  create (temp2, a_2);

  // first set to zero, then call assign, since assign only sets
  // valid regions (that way ghost cells are set to zero)

  setToZero(temp1);
  setToZero(temp2);

  assign(temp1, a_1);
  assign(temp2, a_2);

  // now set covered regions to zero
  for (int level =m_lBase; level<temp1.size()-1; level++)
    {
      LevelData<T>& temp1Level = *temp1[level];
      LevelData<T>& temp2Level = *temp2[level];

      CH_assert(temp1[level]->getBoxes() == temp2[level]->getBoxes());
      CH_assert(temp1[level+1] != NULL);

      int nRefFine = m_refRatios[level];
      const DisjointBoxLayout& finerGrids = temp1[level+1]->getBoxes();
      const DisjointBoxLayout& levelGrids = temp1[level]->getBoxes();
      DataIterator levelDit = levelGrids.dataIterator();
      LayoutIterator finerLit = finerGrids.layoutIterator();

      for (levelDit.begin(); levelDit.ok(); ++levelDit)
        {
          const Box& thisBox = levelGrids[levelDit];
          for (finerLit.begin(); finerLit.ok(); ++finerLit)
            {
              Box testBox = finerGrids[finerLit];
              testBox.coarsen(nRefFine);
              testBox &= thisBox;
              if (!testBox.isEmpty())
                {
                  temp1Level[levelDit].setVal(0.0, testBox,
                                              0, temp1Level.nComp());

                  temp2Level[levelDit].setVal(0.0, testBox,
                                              0, temp2Level.nComp());
                }
            } // end loop over finer boxes
        } // end loop over boxes on this level
    } // end loop over levels for setting covered regions to zero;

  // now loop over levels and call AMRLevelOp dotProduct
  Real prod = 0.0;
  for (int level =m_lBase; level<temp1.size(); level++)
    {
      Real levelProd = m_vectOperators[level]->dotProduct(*temp1[level],
                                                          *temp2[level]);

      // incorporate scaling for AMR dot products
      RealVect dxLev = m_vectDx[level];
      Real scale = D_TERM(dxLev[0], *dxLev[1], *dxLev[2]);
      levelProd *= scale;
      prod += levelProd;
    }

  clear(temp1);
  clear (temp2);

  return prod;
}

///
/**
   Increment by scaled amount (a_lhs += a_scale*a_x).
*/
template<class T>
void
MultilevelLinearOp<T>::incr  (Vector<LevelData<T>* >& a_lhs,
                              const Vector<LevelData<T>* >& a_x,
                              Real a_scale)
{
  CH_TIME("MultilevelLinearOp::incr");

  // int this case, only operate on lBase and finer
  for (int level=m_lBase; level<a_lhs.size(); level++)
    {
      m_vectOperators[level]->incr(*a_lhs[level], *a_x[level], a_scale);
    }
}

///
/**
   Set input to a scaled sum (a_lhs = a_a*a_x + a_b*a_y).
*/
template<class T>
void
MultilevelLinearOp<T>::axby(Vector<LevelData<T>* >& a_lhs,
                            const Vector<LevelData<T>* >& a_x,
                            const Vector<LevelData<T>* >& a_y,
                            Real a_a,
                            Real a_b)
{
  CH_TIME("MultilevelLinearOp::axby");

  // only do this for lBase and finer
  for (int level=m_lBase; level<a_lhs.size(); ++level)
    {
      m_vectOperators[level]->axby(*a_lhs[level],
                                   *a_x[level],
                                   *a_y[level],
                                   a_a, a_b);
    }
}

///
/**
   Multiply the input by a given scale (a_lhs *= a_scale).
*/
template<class T>
void
MultilevelLinearOp<T>::scale(Vector<LevelData<T>* >& a_lhs,
                             const Real& a_scale)
{
  CH_TIME("MultilevelLinearOp::scale");

  // only do this for lBase and finer
  for (int level =m_lBase; level<a_lhs.size(); level++)
    {
      m_vectOperators[level]->scale(*a_lhs[level], a_scale);
    }

}

///
/**
   Return the norm of  a_rhs.
   a_ord == 0  max norm, a_ord == 1 sum(abs(a_rhs)), else, L(a_ord) norm.
*/
template<class T>
Real
MultilevelLinearOp<T>::norm(const Vector<LevelData<T>* >& a_rhs,
                            int a_ord)
{
  CH_TIME("MultilevelLinearOp::norm");

  // now loop over levels and call AMRLevelOp AMRNorm
  int maxLevel = a_rhs.size()-1;
  Real thisNorm = 0.0;
  for (int level = m_lBase; level<a_rhs.size(); level++)
    {
      Real normLevel;

      if (level < maxLevel)
        {
          // finer level exists
          int refRatio = m_refRatios[level];
          normLevel = m_vectOperators[level]->AMRNorm(*a_rhs[level],
                                                      *a_rhs[level+1],
                                                      refRatio,
                                                      a_ord);
        }
      else
        {
          // no finer level exists
          LevelData<T> tempLDF;
          int refRatio = -1;
          normLevel = m_vectOperators[level]->AMRNorm(*a_rhs[level],
                                                      tempLDF,
                                                      refRatio,
                                                      a_ord);
        }

      // now need to do a bit of rescaling to make things work out right
      if (a_ord > 2)
        {
          MayDay::Error("MultilevelLinearOp::norm -- undefined for order > 2");
        }

      if (a_ord == 2)
        {
          normLevel = normLevel*normLevel;
        }

      if (a_ord != 0)
        {
          thisNorm += normLevel;
        }
      else if (normLevel > thisNorm)
        {
          thisNorm = normLevel;
        }
    } // end

  if (a_ord == 2)
    {
      thisNorm = sqrt(thisNorm);
    }

  return thisNorm;
}

///
/**
   Set a_lhs to zero.
*/
template<class T>
void
MultilevelLinearOp<T>::setToZero(Vector<LevelData<T>* >& a_lhs)
{
  CH_TIME("MultilevelLinearOp::setToZero");

  // do this for lBase -1 if necessary, to account for cases where correction
  // lBase-1 needs to be zero
  int coarsestLevel = Max(0, m_lBase-1);
  for (int level = coarsestLevel; level< a_lhs.size(); level++)
    {
      m_vectOperators[level]->setToZero(*a_lhs[level]);
    }
}

template<class T>
void
MultilevelLinearOp<T>::write(const Vector<LevelData<T>* >* a_data,
                             const char*                   a_filename)
{
#ifdef CH_USE_HDF5
  writeVectorLevelName(a_data, &m_refRatios, a_filename);
#else
  MayDay::Warning("MultilevelLinearOp<T>::write unimplemented since CH_USE_HDF5 undefined");
#endif
}

#include "NamespaceFooter.H"

#endif
