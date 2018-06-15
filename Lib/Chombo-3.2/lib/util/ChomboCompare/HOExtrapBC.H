#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Dan Martin, Fri, Jan 14, 2000

#ifndef _HOEXTRAPBC_H_
#define _HOEXTRAPBC_H_

#include "GhostBC.H"

#include "UsingNamespace.H"

/// class to do extrapolation BC's
class HOExtrapBC : public BoxGhostBC
{
public:
  ///
  friend class DomainGhostBC;

  ///
  HOExtrapBC();

  ///
  ~HOExtrapBC();

  ///
  HOExtrapBC(int dir, Side::LoHiSide sd,
             const Interval& a_comps);

  ///
  HOExtrapBC(int dir, Side::LoHiSide sd);


protected:

  ///
  virtual void applyHomogeneousBCs(FArrayBox& a_state,
                                   const ProblemDomain& a_domain,
                                   Real a_dx) const;

  ///
  virtual void applyInhomogeneousBCs(FArrayBox& a_state,
                                     const ProblemDomain& a_domain,
                                     Real a_dx) const;

  ///
  virtual void applyHOExtrapBCs(FArrayBox& a_state,
                              const ProblemDomain& a_domain,
                              Real a_dx) const;

  /// this doesn't do anything since we don't use the factors
  virtual void fillBCValues(FArrayBox& a_neumfac,
                            FArrayBox& a_dircfac,
                            FArrayBox& a_inhmval,
                            Real a_dx,
                            const ProblemDomain& a_domain) const;

  /// deprecated interface
  virtual void applyHomogeneousBCs(FArrayBox& a_state,
                                   const Box& a_domain,
                                   Real a_dx) const;

  /// deprecated interface
  virtual void applyInhomogeneousBCs(FArrayBox& a_state,
                                     const Box& a_domain,
                                     Real a_dx) const;

  /// deprecated interface
  virtual void applyHOExtrapBCs(FArrayBox& a_state,
                              const Box& a_domain,
                              Real a_dx) const;

  /// this doesn't do anything since we don't use the factors
  virtual void fillBCValues(FArrayBox& a_neumfac,
                            FArrayBox& a_dircfac,
                            FArrayBox& a_inhmval,
                            Real a_dx,
                            const Box& a_domain) const;

  ///
  BoxGhostBC* new_boxghostbc() const;

private:
  HOExtrapBC(const HOExtrapBC& extrapin)
  {
  }

  void operator=(const HOExtrapBC& extrapin)
  {
  }
};

#endif
