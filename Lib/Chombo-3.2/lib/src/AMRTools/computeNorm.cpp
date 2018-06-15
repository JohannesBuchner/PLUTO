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
using std::pow;
using std::sqrt;
#ifdef CH_MPI
#include <mpi.h>
#endif

#include "LayoutIterator.H"
#include "BoxLayoutData.H"

#include "computeNorm.H"
#include "NamespaceHeader.H"

Real computeNorm(const Vector<LevelData<FArrayBox>* >& a_phi,
                 const Vector<int>&                    a_nRefFine,
                 const Real                            a_dxCrse,
                 const Interval                        a_comps,
                 const int                             a_p,
                 const int                             a_lBase)
{
  Real dxLevel = a_dxCrse;
  int numLevels = a_phi.size();
  Real norm;
  Real normLevel;

  // it is often the case that while a_phi has many possible
  // levels, only a subset of them are defined -- check that
  // just to be sure
  if (a_phi[numLevels-1] == NULL)
  {
    int lev = numLevels-1;
    while (a_phi[lev] == NULL)
    {
      lev--;
    }
    numLevels = lev+1;
  }

  norm = 0.0;

  // loop over levels
  for (int lev = a_lBase; lev < numLevels; lev++)
  {
    int refRatio = -1;
    //in case there are extra levels which are not defined
    if (a_phi[lev] != NULL)
    {
      CH_assert(a_phi[lev]->isDefined());

      LevelData<FArrayBox>& thisPhi = *(a_phi[lev]);
      const DisjointBoxLayout* finerGridsPtr = NULL;

      if (lev < numLevels-1)
      {
        finerGridsPtr = &(a_phi[lev+1]->getBoxes());
        refRatio = a_nRefFine[lev];
      }

      normLevel = computeNorm(thisPhi, finerGridsPtr, refRatio,
                              dxLevel, a_comps, a_p);

      if (a_p != 0)
      {
        Real p = a_p;
        normLevel = pow(normLevel,p);

        norm += normLevel;
      }
      else if (normLevel > norm)
      {
        norm = normLevel;
      }
    }

    // update  dxLevel
    dxLevel = dxLevel/refRatio;

  }

  // shouldn't need to do broadcast/gather thing

  // now take the inverse power of the norm
  if (a_p != 0)
  {
    Real invp = 1.0/a_p;
    norm = pow(norm,invp);
  }

  return norm;
}

Real computeNorm(const LevelData<FArrayBox>& a_phi,
                 const DisjointBoxLayout*    a_finerGridsPtr,
                 const int                   a_nRefFine,
                 const Real                  a_dx,
                 const Interval              a_comps,
                 const int                   a_p)
{
  Real normLevel;
  normLevel = 0.0;

  const DisjointBoxLayout& levelGrids = a_phi.getBoxes();

  LevelData<FArrayBox> temp(levelGrids, a_comps.size());

  Interval tempComps(0,a_comps.size()-1);

  DataIterator dit = a_phi.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    temp[dit()].copy(temp[dit()].box(), tempComps,
                     temp[dit()].box(), a_phi[dit()], a_comps);

    if (a_finerGridsPtr != NULL)
    {
      LayoutIterator litFine = a_finerGridsPtr->layoutIterator();

      // now loop over fine boxes and set covered regions to 0
      for (litFine.reset(); litFine.ok(); ++litFine)
      {
        Box coveredBox(a_finerGridsPtr->get(litFine()));
        coveredBox.coarsen(a_nRefFine);
        coveredBox &= temp[dit()].box();

        if (!coveredBox.isEmpty())
        {
          temp[dit()].setVal(0.0, coveredBox, 0, tempComps.size());
        }
      }  // end loop over fine-grid boxes
    } // end if there is a finer level
  } // end loop over this level's grids

  // compute unscaled norm
  normLevel = norm(temp, tempComps, a_p);

  if (a_p != 0)
  {
    Real exponent = SpaceDim;
    exponent /= a_p;
    Real scale = pow(a_dx,exponent);
    normLevel *= scale;
  }

  return normLevel;
}

Real computeNorm(const Vector<LevelData<FArrayBox>* >& a_phi,
                 const Vector<IntVect>&                a_nRefFine,
                 const RealVect&                       a_dxCrse,
                 const Interval                        a_comps,
                 const int                             a_p,
                 const int                             a_lBase)
{
  RealVect dxLevel = a_dxCrse;
  int numLevels = a_phi.size();
  Real norm;
  Real normLevel;

  // it is often the case that while a_phi has many possible
  // levels, only a subset of them are defined -- check that
  // just to be sure
  if (a_phi[numLevels-1] == NULL)
  {
    int lev = numLevels-1;
    while (a_phi[lev] == NULL)
    {
      lev--;
    }
    numLevels = lev+1;
  }

  norm = 0.0;

  // loop over levels
  for (int lev = a_lBase; lev < numLevels; lev++)
  {
    IntVect refRatio = -IntVect::Unit;
    //in case there are extra levels which are not defined
    if (a_phi[lev] != NULL)
    {
      CH_assert(a_phi[lev]->isDefined());

      LevelData<FArrayBox>& thisPhi = *(a_phi[lev]);
      const DisjointBoxLayout* finerGridsPtr = NULL;

      if (lev < numLevels-1)
      {
        finerGridsPtr = &(a_phi[lev+1]->getBoxes());
        refRatio = a_nRefFine[lev];
      }

      normLevel = computeNorm(thisPhi, finerGridsPtr, refRatio,
                              dxLevel, a_comps, a_p);

      if (a_p != 0)
      {
        Real p = a_p;
        normLevel = pow(normLevel,p);

        norm += normLevel;
      }
      else if (normLevel > norm)
      {
        norm = normLevel;
      }
    }

    // update  dxLevel
    dxLevel = dxLevel/RealVect(refRatio);

  }

  // shouldn't need to do broadcast/gather thing

  // now take the inverse power of the norm
  if (a_p != 0)
  {
    Real invp = 1.0/a_p;
    norm = pow(norm,invp);
  }

  return norm;
}

Real computeNorm(const LevelData<FArrayBox>& a_phi,
                 const DisjointBoxLayout*    a_finerGridsPtr,
                 const IntVect&              a_nRefFine,
                 const RealVect&             a_dx,
                 const Interval              a_comps,
                 const int                   a_p)
{
  Real normLevel;
  normLevel = 0.0;

  const DisjointBoxLayout& levelGrids = a_phi.getBoxes();

  LevelData<FArrayBox> temp(levelGrids, a_comps.size());

  Interval tempComps(0,a_comps.size()-1);

  DataIterator dit = a_phi.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    temp[dit()].copy(temp[dit()].box(), tempComps,
                     temp[dit()].box(), a_phi[dit()], a_comps);

    if (a_finerGridsPtr != NULL)
    {
      LayoutIterator litFine = a_finerGridsPtr->layoutIterator();

      // now loop over fine boxes and set covered regions to 0
      for (litFine.reset(); litFine.ok(); ++litFine)
      {
        Box coveredBox(a_finerGridsPtr->get(litFine()));
        coveredBox.coarsen(a_nRefFine);
        coveredBox &= temp[dit()].box();

        if (!coveredBox.isEmpty())
        {
          temp[dit()].setVal(0.0, coveredBox, 0, tempComps.size());
        }
      }  // end loop over fine-grid boxes
    } // end if there is a finer level
  } // end loop over this level's grids

  // compute unscaled norm
  normLevel = norm(temp, tempComps, a_p);

  if (a_p != 0)
  {
    Real dV = a_dx.product();
    Real exponent = 1. / Real(a_p);
    Real scale = pow(dV, exponent);

    normLevel *= scale;
  }

  return normLevel;
}

Real computeMax(const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<int>&                    a_nRefFine,
                const Interval                        a_comps,
                const int                             a_lBase)
{
  int numLevels = a_phi.size();

  // it is often the case that while a_phi has many possible
  // levels, only a subset of them are defined -- check that
  // just to be sure
  if (a_phi[numLevels-1] == NULL)
  {
    int lev= numLevels-1;

    while (a_phi[lev] == NULL)
    {
      lev--;
    }

    numLevels = lev+1;
  }

  Real max;
  Real maxOnLevel;

  max = -HUGE_VAL;

  // loop over levels
  for (int lev = a_lBase; lev < numLevels; lev++)
  {
    //in case there are extra levels which are not defined
    if (a_phi[lev] != NULL)
    {
      CH_assert(a_phi[lev]->isDefined());

      LevelData<FArrayBox>& thisPhi = *(a_phi[lev]);
      const DisjointBoxLayout* finerGridsPtr = NULL;

      if (lev < numLevels-1)
      {
        finerGridsPtr = &(a_phi[lev+1]->getBoxes());

        maxOnLevel = computeMax(thisPhi, finerGridsPtr, a_nRefFine[lev],
                                a_comps);
      }
      else
      {
        int bogusRefRatio = 100000;
        maxOnLevel = computeMax(thisPhi, finerGridsPtr, bogusRefRatio,
                                a_comps);
      }

      if (maxOnLevel > max)
      {
        max = maxOnLevel;
      }
    }
  }

  // shouldn't need to do broadcast/gather thing
  return max;
}

Real computeMax(const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGridsPtr,
                const int                   a_nRefFine,
                const Interval              a_comps)
{

  Real levelMax;
  Real thisMax;

  levelMax = -99999999999.9;

  const DisjointBoxLayout& levelGrids = a_phi.getBoxes();

  LevelData<FArrayBox> temp(levelGrids, a_comps.size());

  Interval tempComps(0,a_comps.size()-1);

  DataIterator dit = a_phi.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    temp[dit()].copy(temp[dit()].box(), tempComps,
                     temp[dit()].box(), a_phi[dit()], a_comps);

    if (a_finerGridsPtr != NULL)
    {
      LayoutIterator litFine = a_finerGridsPtr->layoutIterator();

      // now loop over fine boxes and set covered regions to 0
      for (litFine.reset(); litFine.ok(); ++litFine)
      {
        Box coveredBox(a_finerGridsPtr->get(litFine()));
        coveredBox.coarsen(a_nRefFine);
        coveredBox &= temp[dit()].box();

        if (!coveredBox.isEmpty())
        {
          temp[dit()].setVal(levelMax, coveredBox, 0, tempComps.size());
        }
      }  // end loop over fine-grid boxes
    } // end if there is a finer level

    // while we're looping over the grids, get max as well
    // need to loop over comps
    for (int comp = tempComps.begin(); comp <= tempComps.end(); ++comp)
    {
      thisMax = temp[dit()].max(comp);
      if (thisMax > levelMax)
      {
        levelMax = thisMax;
      }
    }
  } // end loop over this level's grids

  // need to do MPI gather here
#ifdef CH_MPI
  Real recv;
  int result = MPI_Allreduce(&levelMax, &recv, 1, MPI_CH_REAL,
                             MPI_MAX, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
  {
    MayDay::Error("Sorry, but I had a communication error in computeNorm");
  }

  levelMax = recv;
#endif

  return levelMax;
}

Real computeMin(const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<int>&                    a_nRefFine,
                const Interval                        a_comps,
                const int                             a_lBase)
{
  int numLevels = a_phi.size();

  // it is often the case that while a_phi has many possible
  // levels, only a subset of them are defined -- check that
  // just to be sure
  if (a_phi[numLevels-1] == NULL)
  {
    int lev = numLevels-1;

    while (a_phi[lev] == NULL)
    {
      lev--;
    }

    numLevels = lev+1;
  }

  Real min;
  Real minOnLevel;

  min = HUGE_VAL;

  // loop over levels
  for (int lev = a_lBase; lev < numLevels; lev++)
  {
    //in case there are extra levels which are not defined
    if (a_phi[lev] != NULL)
    {
      CH_assert(a_phi[lev]->isDefined());

      LevelData<FArrayBox>& thisPhi = *(a_phi[lev]);
      const DisjointBoxLayout* finerGridsPtr = NULL;

      if (lev < numLevels-1)
      {
        finerGridsPtr = &(a_phi[lev+1]->getBoxes());

        minOnLevel = computeMin(thisPhi, finerGridsPtr, a_nRefFine[lev],
                                a_comps);
      }
      else
      {
        int bogusRefRatio = 1000000;
        minOnLevel = computeMin(thisPhi, finerGridsPtr, bogusRefRatio,
                                a_comps);
      }

      if (minOnLevel < min)
      {
        min = minOnLevel;
      }
    }
  }

  // shouldn't need to do broadcast/gather thing
  return min;
}

Real computeMin(const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGridsPtr,
                const int                   a_nRefFine,
                const Interval              a_comps)
{

  Real levelMin;
  Real thisMin;

  levelMin = 99999999999.9;

  const DisjointBoxLayout& levelGrids = a_phi.getBoxes();

  LevelData<FArrayBox> temp(levelGrids, a_comps.size());

  Interval tempComps(0,a_comps.size()-1);

  DataIterator dit = a_phi.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    temp[dit()].copy(temp[dit()].box(), tempComps,
                     temp[dit()].box(), a_phi[dit()], a_comps);

    if (a_finerGridsPtr != NULL)
    {
      LayoutIterator litFine = a_finerGridsPtr->layoutIterator();

      // now loop over fine boxes and set covered regions to 0
      for (litFine.reset(); litFine.ok(); ++litFine)
      {
        Box coveredBox(a_finerGridsPtr->get(litFine()));
        coveredBox.coarsen(a_nRefFine);
        coveredBox &= temp[dit()].box();

        if (!coveredBox.isEmpty())
        {
          temp[dit()].setVal(levelMin, coveredBox, 0, tempComps.size());
        }
      }  // end loop over fine-grid boxes
    } // end if there is a finer level

    // while we're looping over the grids, get min as well
    // need to loop over comps
    for (int comp = tempComps.begin(); comp <= tempComps.end(); ++comp)
    {
      thisMin = temp[dit()].min(comp);
      if (thisMin < levelMin)
      {
        levelMin = thisMin;
      }
    }
  } // end loop over this level's grids

#ifdef CH_MPI
  Real recv;
  int result = MPI_Allreduce(&levelMin, &recv, 1, MPI_CH_REAL,
                             MPI_MIN, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in computeNorm");
    }

  levelMin = recv;
#endif

  return levelMin;
}

Real norm(const BoxLayoutData<FluxBox>& a_layout,
          const Interval&               a_interval,
          const int&                    a_p)
{
  int startComp = a_interval.begin();
  int numComp = a_interval.size();
  const BoxLayout& validRegion = a_layout.boxLayout();
  DataIterator dit = validRegion.dataIterator();
  Real norm = 0;
  if (a_p == 0)
    {
      for (dit.reset(); dit.ok(); ++dit)
        {
          const Box& validBox = validRegion[dit];
          const FluxBox& thisFlux = a_layout[dit];
          for (int dir = 0; dir < SpaceDim; dir++)
            {
              Box faceBox = surroundingNodes(validBox, dir);
              norm = Max(norm, thisFlux[dir].norm(faceBox, 0,
                                                  startComp, numComp));
            }
        }
#ifdef CH_MPI
      Real recv;
      int result = MPI_Allreduce(&norm, &recv, 1, MPI_CH_REAL,
                                 MPI_MAX, Chombo_MPI::comm);
      if (result != MPI_SUCCESS)
        {
          MayDay::Error("sorry, but I had a communication error in norm");
        }
      norm = recv;
#endif
    }
  else if (a_p == 1) // abs sum norm
    {
      for (dit.reset(); dit.ok(); ++dit)
        {
          const Box& validBox = validRegion[dit];
          const FluxBox& thisFlux = a_layout[dit];
          for (int dir = 0; dir < SpaceDim; ++dir)
            {
              Box faceBox = surroundingNodes(validBox, dir);
              norm += thisFlux[dir].norm(faceBox, 1,
                                         startComp, numComp);
            }
        }
#ifdef CH_MPI
      Real recv;
      int result = MPI_Allreduce(&norm, &recv, 1, MPI_CH_REAL,
                                 MPI_SUM, Chombo_MPI::comm);
      if (result != MPI_SUCCESS)
        {
          MayDay::Error("sorry, but I had a communication error in norm");
        }
      norm = recv;
#endif
    }
  else
    {
      for (dit.reset(); dit.ok(); ++dit)
        {
          const Box& validBox = validRegion[dit];
          const FluxBox& thisFlux = a_layout[dit];
          for (int dir = 0; dir < SpaceDim; dir++)
            {
              Box faceBox = surroundingNodes(validBox, dir);
              norm += thisFlux[dir].sumPow(faceBox, a_p,
                                           startComp, numComp);
            }
        }
#ifdef CH_MPI
      Real recv;
      int result = MPI_Allreduce(&norm, &recv, 1, MPI_CH_REAL,
                                 MPI_SUM, Chombo_MPI::comm);
      if (result != MPI_SUCCESS)
        {
          MayDay::Error("sorry, but I had a communication error on norm");
        }
      norm = recv;
#endif
      Real invpwr = 1.0/a_p;
      norm = pow(norm, invpwr);
    }
  return norm;
}
#include "NamespaceFooter.H"
