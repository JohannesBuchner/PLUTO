#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
using namespace std;

#include "MayDay.H"
#include "CHOMBO_VERSION.H"
#ifdef CH_MPI
#include "mpi.h"
#endif

#ifdef CH_MULTIDIM
#undef CH_MULTIDIM
#endif

#include "FortranNameMacro.H"
#include "BaseNamespaceHeader.H"

extern "C"
{
  void
  FORTRAN_NAME(MAYDAYERROR,maydayerror) (void)
  {
    MayDay::Error( "A ChomboFortran routine called MAYDAYERROR().  Rerun with the debugger to find where." ) ;
  }

  void
  FORTRAN_NAME(MAYDAYABORT,maydayabort) (void)
  {
    MayDay::Abort( "A ChomboFortran routine called MAYDAYABORT().  Rerun with the debugger to find where." ) ;
  }

  void
  FORTRAN_NAME(MAYDAYEXIT,maydayexit) (int* exit_code)
  {
    MayDay::Error( "A ChomboFortran routine called MAYDAYEXIT().  Rerun with the debugger to find where." ,*exit_code ) ;
  }

  void
  FORTRAN_NAME(MAYDAY_ERROR,mayday_error) (void)
  {
    MayDay::Error( "A ChomboFortran routine called MAYDAY_ERROR().  Rerun with the debugger to find where." ) ;
  }

  void
  FORTRAN_NAME(MAYDAY_ABORT,mayday_abort) (void)
  {
    MayDay::Abort( "A ChomboFortran routine called MAYDAY_ABORT().  Rerun with the debugger to find where." ) ;
  }

  void
  FORTRAN_NAME(MAYDAY_EXIT,mayday_exit) (int* exit_code)
  {
    MayDay::Error( "A ChomboFortran routine called MAYDAY_EXIT().  Rerun with the debugger to find where." ,*exit_code ) ;
  }

}
#include "BaseNamespaceFooter.H"
