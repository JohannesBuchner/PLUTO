#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VELOCITYAMRPOISSONOP_H_
#define _VELOCITYAMRPOISSONOP_H_

#include "AMRPoissonOp.H"
#include "Tuple.H"
#include "VelBCHolder.H"

#include "UsingNamespace.H"

class VelocityAMRPoissonOp: public AMRPoissonOp
{

public:

  ///
  VelocityAMRPoissonOp();

  ///
  ~VelocityAMRPoissonOp();

  ///
  /** The raison d'etre for this class -- want to be able to
      apply operator with different BCs for each component
      in a_phi. */
  void applyOpI(LevelData<FArrayBox>&        a_lhs,
                const LevelData<FArrayBox>&  a_phi,
                bool                         a_homogeneous,
                VelBCHolder&                 a_velBC);

  ///
  void AMROperatorNF(LevelData<FArrayBox>&        a_LofPhi,
                     const LevelData<FArrayBox>&  a_phi,
                     const LevelData<FArrayBox>&  a_phiCoarse,
                     bool                         a_homogeneousPhysBC,
                     VelBCHolder&                 a_velBC);

  ///
  void applyOp(LevelData<FArrayBox>&        a_lhs,
               const LevelData<FArrayBox>&  a_phi,
               const LevelData<FArrayBox>*  a_phiCoarsePtr,
               bool                         a_homogeneous,
               VelBCHolder&                 a_velBC);

protected:

private:

};

#endif
