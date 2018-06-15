#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PARSTREAM_H_
#define _PARSTREAM_H_

#include <iostream>
#include <string>
#include "BaseNamespaceHeader.H"

/// Use this in place of std::cout for program output.
/**
 * Replaces std::cout in most of the Chombo code.  In serial this just
 * returns std::cout.  In parallel, this creates a separate file for each
 * proc called <basename>.n where n is the procID and <basename> defaults to
 * "pout" but can be set by calling setPoutBaseName().  Output is then
 * directed to these files.  This keeps the output from different processors
 * from getting all jumbled up.
     If    you want fewer files, you can do
     setenv CH_OUTPUT_INTERVAL nproc
     and it will only output every nproc processors pout.n files
     (where n%nproc == 0).
*/
std::ostream& pout();

/// Changes the base part of the filename for pout() files.
/**
 * When in parallel, changes the base name of the pout() files.  If pout()
 * has already been called, it closes the current output file and opens a
 * new one (unless the name is the same, in which case it does nothing).
 * In serial, ignores the argument and does nothing.
*/
void setPoutBaseName(const std::string & a_Name);

/// Accesses the filename for the local pout() file.
/**
 * Returns the name used for the local pout() file.  In parallel this is
 * "<pout_basename>.<procID>", where <pout_basename> defaults to "pout" and
 * can be modified by calling setPoutBaseName(), and <procID> is the local
 * proc number.  In serial, this always returns the string "cout".  It is an
 * error (exit code 111) to call this in parallel before MPI_Initialize().
*/
const std::string & poutFileName();

#include "BaseNamespaceFooter.H"
#endif // _PARSTREAM_H_
