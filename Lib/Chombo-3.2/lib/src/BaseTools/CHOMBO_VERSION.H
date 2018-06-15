#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CHOMBO_VERSION_H_
#define _CHOMBO_VERSION_H_
#include "BaseNamespaceHeader.H"

#define CHOMBO_VERSION 3.2

/** signal that this is the svn trunk version of Chombo, rather than 
    the release version. This should be changed in the release version 
    of this file */
#define CHOMBO_TRUNK FALSE

/** which patch of the release is this?
 */
#define CHOMBO_PATCH_NUMBER 2

/**
   \mainpage Chombo: Block-Structured Adaptive Mesh Refinement Library

   <a href="http://seesar.lbl.gov/anag/chombo">Chombo Web Page</a>

   <table WIDTH=640 HEIGHT=26 BORDER=0 CELLPADDING=0 CELLSPACING=0>
   <tr>
   <td WIDTH=106 BGCOLOR=#303784>&nbsp;</td>
   <td WIDTH=534 BGCOLOR=#FFEF89>
   <a NAME="What is Chombo?">What is Chombo?</a></b></td>
   </tr>
   </table><br CLEAR=ALL>

   The Chombo package provides a set of tools for implementing finite
   difference methods for the solution of partial differential equations
   on block-structured adaptively refined rectangular grids.  Both elliptic
   and time-dependent modules are included.  Support for parallel platforms
   and standardized self-describing file formats are included.

   <p>

   Chombo provides a distributed infrastructure for parallel calculations
   over block-structured, adaptively refined grids.  Chombo's design is
   uniquely flexible and accessible.  Any collaborator will be able to
   develop parallel applications to solve the partial differential equations
   in which she is interested with far shorter development times than would
   be possible without the infrastructure.  Very careful design and
   documentation allows said collaborator to enter the software at many
   levels.  She  will be able to use Chombo to investigate deep technical
   issues of adaptive mesh refinement algorithms or to simply adapt the
   example applications to solve different scientific problems.
*/

/**
 *  @defgroup chombo Chombo
 *
 */

/** @defgroup eb EBChombo
 *
 */

/** @defgroup mf MFChombo
 *
 */

/**
   @ingroup macros
   The major version number of Chombo -- an integer.

   The Chombo version number has the form:

               CHOMBO_VERSION_MAJOR.CHOMBO_VERSION_MINOR
*/
#define CHOMBO_VERSION_MAJOR 3

/**
   @ingroup macros
   The minor version number of Chombo -- an integer.
*/
#define CHOMBO_VERSION_MINOR 2

#include "BaseNamespaceFooter.H"
#endif
