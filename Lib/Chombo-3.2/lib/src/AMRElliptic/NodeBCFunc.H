#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NODEBCFUNC_H_
#define _NODEBCFUNC_H_

#include "BCFunc.H"
#include "NodeFArrayBox.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

///
/**
 */
typedef void(*NodeBCFunc)(NodeFArrayBox&       a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous);

///
/**
   neum bc for a particular side
*/
void NodeNeumBC(NodeFArrayBox& a_state,
                const Box&     a_valid,
                Real           a_dx,
                bool           a_homogeneous,
                BCValueHolder  a_value,
                int            a_dir,
                Side::LoHiSide a_side);

///
/**
   neum bc for all sides
*/
void NodeNeumBC(NodeFArrayBox& a_state,
                const Box&     a_valid,
                Real           a_dx,
                bool           a_homogeneous,
                BCValueHolder  a_value);

///
/**
 */
void NodeDiriBC(NodeFArrayBox& a_state,
                const Box&     a_valid,
                Real           a_dx,
                bool           a_homogeneous,
                BCValueHolder  a_value,
                int            a_dir,
                Side::LoHiSide a_side,
                int order = 2);

///
/**
   diri bc for all  sides
*/
void NodeDiriBC(NodeFArrayBox& a_state,
                const Box&     a_valid,
                Real           a_dx,
                bool           a_homogeneous,
                BCValueHolder  a_value,
                int order = 2);

#include "NamespaceFooter.H"
#endif
