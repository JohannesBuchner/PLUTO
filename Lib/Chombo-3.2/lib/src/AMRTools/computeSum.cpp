#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LayoutIterator.H"

#include "computeSum.H"
#include "NamespaceHeader.H"

Real computeSum(const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<int>&                    a_nRefFine,
                const Real&                           a_dxCrse,
                const Interval&                       a_comps,
                const int&                            a_lBase)
{
  Real volume;

  return computeSum(volume,a_phi,a_nRefFine,a_dxCrse,a_comps,a_lBase);
}

Real computeSum(Real&                                 a_volume,
                const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<int>&                    a_nRefFine,
                const Real&                           a_dxCrse,
                const Interval&                       a_comps,
                const int&                            a_lBase)
{
  Real dxLevel = a_dxCrse;
  int numLevels = a_phi.size();
  Real sum, sumLevel;
  Real vol, volLevel;

  sum = 0.0;
  vol = 0.0;

  // loop over levels
  for (int lev = a_lBase; lev < numLevels; lev++)
  {
    int refRatio = -1;
    //in case there are extra levels which are not defined
    if (a_phi[lev] != NULL)
    {
      CH_assert(a_phi[lev]->isDefined());

      LevelData<FArrayBox>& thisPhi = *(a_phi[lev]);
      const DisjointBoxLayout* finerGridsPtr;

      if ((lev < numLevels-1) && a_phi[lev+1] !=NULL)
      {
        finerGridsPtr = &(a_phi[lev+1]->getBoxes());
        refRatio = a_nRefFine[lev];
      }
      else
      {
        finerGridsPtr = NULL;
      }

      sumLevel = computeSum(volLevel,thisPhi,finerGridsPtr,refRatio,
                            dxLevel,a_comps,false);
      sum += sumLevel;
      vol += volLevel;

      dxLevel = dxLevel/refRatio;
    }
  }

  // do broadcast/gather thing here
#ifdef CH_MPI
  Real recv;
  int result = MPI_Allreduce(&sum, &recv, 1, MPI_CH_REAL,
                             MPI_SUM, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
  {
    MayDay::Error("Sorry, but I had a communication error in computeSum");
  }

  sum = recv;

  result = MPI_Allreduce(&vol, &recv, 1, MPI_CH_REAL,
                         MPI_SUM, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
  {
    MayDay::Error("Sorry, but I had a communication error in computeSum");
  }

  vol = recv;
#endif

  a_volume = vol;

  return sum;
}

Real computeSum(const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGridsPtr,
                const int&                  a_nRefFine,
                const Real&                 a_dx,
                const Interval&             a_comps,
                const bool                  a_global)
{
  Real volume;

  return computeSum(volume,a_phi,a_finerGridsPtr,a_nRefFine,a_dx,a_comps,
                    a_global);
}

Real computeSum(Real&                       a_volume,
                const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGridsPtr,
                const int&                  a_nRefFine,
                const Real&                 a_dx,
                const Interval&             a_comps,
                const bool                  a_global)
{
  Real sum, sumLevel;
  Real vol, volLevel;

  sum = 0.0;
  vol = 0.0;

  const DisjointBoxLayout& levelGrids = a_phi.getBoxes();

  LevelData<FArrayBox> temp(levelGrids,a_comps.size());
  LevelData<FArrayBox> volTemp(levelGrids,1);

  Interval tempComps(0,a_comps.size()-1);

  DataIterator dit = a_phi.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    Box copyBox = temp[dit()].box();
    temp[dit()].copy(copyBox,tempComps,copyBox,a_phi[dit()],a_comps);
    volTemp[dit()].setVal(1.0);

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
          temp[dit()].setVal(0.0,coveredBox,0,tempComps.size());
          volTemp[dit()].setVal(0.0,coveredBox,0,1);
        }
      }  // end loop over fine-grid boxes
    } // end if there is a finer level

    sumLevel = temp[dit()].sum(0,tempComps.size());
    sum += sumLevel;

    volLevel = volTemp[dit()].sum(0,1);
    vol += volLevel;
  } // end loop over this level's grids

  for (int i = 0; i < SpaceDim; i++)
  {
    sum *= a_dx;
    vol *= a_dx;
  }

  // do broadcast/gather thing here
#ifdef CH_MPI
  if (a_global)
  {
    Real recv;
    int result = MPI_Allreduce(&sum, &recv, 1, MPI_CH_REAL,
                               MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in computeSum");
    }

    sum = recv;

    result = MPI_Allreduce(&vol, &recv, 1, MPI_CH_REAL,
                           MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in computeSum");
    }

    vol = recv;
  }
#endif

  a_volume = vol;

  return sum;
}

// ANISOTROPIC

Real computeSum(const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<IntVect>&                a_nRefFine,
                const RealVect&                       a_dxCrse,
                const Interval&                       a_comps,
                const int&                            a_lBase)
{
  Real volume;

  return computeSum(volume,a_phi,a_nRefFine,a_dxCrse,a_comps,a_lBase);
}

Real computeSum(Real&                                 a_volume,
                const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<IntVect>&                a_nRefFine,
                const RealVect&                       a_dxCrse,
                const Interval&                       a_comps,
                const int&                            a_lBase)
{
  RealVect dxLevel = a_dxCrse;
  int numLevels = a_phi.size();
  Real sum, sumLevel;
  Real vol, volLevel;

  sum = 0.0;
  vol = 0.0;

  // loop over levels
  for (int lev = a_lBase; lev < numLevels; lev++)
  {
    IntVect refRatio = -IntVect::Unit;
    //in case there are extra levels which are not defined
    if (a_phi[lev] != NULL)
    {
      CH_assert(a_phi[lev]->isDefined());

      LevelData<FArrayBox>& thisPhi = *(a_phi[lev]);
      const DisjointBoxLayout* finerGridsPtr;

      if ((lev < numLevels-1) && a_phi[lev+1] !=NULL)
      {
        finerGridsPtr = &(a_phi[lev+1]->getBoxes());
        refRatio = a_nRefFine[lev];
      }
      else
      {
        finerGridsPtr = NULL;
      }

      sumLevel = computeSum(volLevel,thisPhi,finerGridsPtr,refRatio,
                            dxLevel,a_comps,false);
      sum += sumLevel;
      vol += volLevel;

      dxLevel = dxLevel/RealVect(refRatio);
    }
  }

  // do broadcast/gather thing here
#ifdef CH_MPI
  Real recv;
  int result = MPI_Allreduce(&sum, &recv, 1, MPI_CH_REAL,
                             MPI_SUM, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
  {
    MayDay::Error("Sorry, but I had a communication error in computeSum");
  }

  sum = recv;

  result = MPI_Allreduce(&vol, &recv, 1, MPI_CH_REAL,
                         MPI_SUM, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
  {
    MayDay::Error("Sorry, but I had a communication error in computeSum");
  }

  vol = recv;
#endif

  a_volume = vol;

  return sum;
}

Real computeSum(const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGridsPtr,
                const IntVect&              a_nRefFine,
                const RealVect&             a_dx,
                const Interval&             a_comps,
                const bool                  a_global)
{
  Real volume;

  return computeSum(volume,a_phi,a_finerGridsPtr,a_nRefFine,a_dx,a_comps,
                    a_global);
}

Real computeSum(Real&                       a_volume,
                const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGridsPtr,
                const IntVect&              a_nRefFine,
                const RealVect&             a_dx,
                const Interval&             a_comps,
                const bool                  a_global)
{
  Real sum, sumLevel;
  Real vol, volLevel;

  sum = 0.0;
  vol = 0.0;

  const DisjointBoxLayout& levelGrids = a_phi.getBoxes();

  LevelData<FArrayBox> temp(levelGrids,a_comps.size());
  LevelData<FArrayBox> volTemp(levelGrids,1);

  Interval tempComps(0,a_comps.size()-1);

  DataIterator dit = a_phi.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    Box copyBox = temp[dit()].box();
    temp[dit()].copy(copyBox,tempComps,copyBox,a_phi[dit()],a_comps);
    volTemp[dit()].setVal(1.0);

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
          temp[dit()].setVal(0.0,coveredBox,0,tempComps.size());
          volTemp[dit()].setVal(0.0,coveredBox,0,1);
        }
      }  // end loop over fine-grid boxes
    } // end if there is a finer level

    sumLevel = temp[dit()].sum(0,tempComps.size());
    sum += sumLevel;

    volLevel = volTemp[dit()].sum(0,1);
    vol += volLevel;
  } // end loop over this level's grids

  for (int i = 0; i < SpaceDim; i++)
  {
    sum *= a_dx[i];
    vol *= a_dx[i];
  }

  // do broadcast/gather thing here
#ifdef CH_MPI
  if (a_global)
  {
    Real recv;
    int result = MPI_Allreduce(&sum, &recv, 1, MPI_CH_REAL,
                               MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in computeSum");
    }

    sum = recv;

    result = MPI_Allreduce(&vol, &recv, 1, MPI_CH_REAL,
                           MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in computeSum");
    }

    vol = recv;
  }
#endif

  a_volume = vol;

  return sum;
}

#include "NamespaceFooter.H"
