#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SPACE_H_
#define _SPACE_H_

#ifdef HAVE_CONFIG_H
#include <autoconfig.H>
#endif

#ifdef CH_LANG_CC
# include "MayDay.H"
#endif

# include "CH_assert.H"

#include "NamespaceHeader.H"

#if !defined(CH_LANG_FORT) || defined(CH_LANG_CC)

/**
  The SpaceDim C++ integer constant specifies the dimension of the
  Chombo library to be built.  It is initialized to the value of
  the macro CH_SPACEDIM.  Hence it will have one of the values 1, 2, or 3.

  C++ code that depends on the dimensionality of the problem that
  you're trying to solve should be cast in terms of this constant,
  instead of the CH_SPACEDIM macro, if possible, as constants are
  visible to debuggers while macros are not.
*/
   const int SpaceDim = CH_SPACEDIM;
#endif

#ifdef DOXYGEN
/** \defgroup macros Chombo Macros
 */

/**
@ingroup macros
compile-time cpp macro that controls the dimensionality of
    the Chombo package.  */
#define CH_SPACEDIM
#undef  CH_SPACEDIM
#define CH_SPACEDIM 3
#endif

#if ! defined(CH_SPACEDIM)
#error CH_SPACEDIM must be defined
#endif

#if (CH_SPACEDIM < 1 || CH_SPACEDIM > 6)
#error CH_SPACEDIM must be either 1, 2, 3, 4, 5, or 6
#endif

#ifdef DOXYGEN
#undef CH_SPACEDIM
#endif

#include "NamespaceFooter.H"

#endif
#ifdef CH_LANG_CC
// That matches the include guard.  The rest of this file is all macros
// whose definitions depend on CH_SPACEDIM.  In CH_MULTIDIM mode, we want
// to rescan these macros after every time CH_SPACEDIM has changed.
// Of course, if CH_MULTIDIM is not defined, we want to scan just once.
#endif

#ifndef LAST_SPACE_H_SPACEDIM
#define LAST_SPACE_H_SPACEDIM 0
#endif

#if CH_SPACEDIM != LAST_SPACE_H_SPACEDIM
#undef  LAST_SPACE_H_SPACEDIM

#ifdef D_EXPR6
#undef D_EXPR
#undef D_DECL
#undef D_TERM
#undef D_INVTERM
#undef D_SELECT
#undef D_EXPR6
#undef D_DECL6
#undef D_TERM6
#undef D_INVTERM6
#undef D_SELECT6
#endif

#if CH_SPACEDIM==1
#  define D_EXPR(a,b,c) ((void)((a),0))
#  define D_DECL(a,b,c) a
#  define D_TERM(a,b,c) a
#  define D_INVTERM(a,b,c) a
#  define D_SELECT(a,b,c) a
#  define D_EXPR6(a,b,c,d,e,f) ((void)((a),0))
#  define D_DECL6(a,b,c,d,e,f) a
#  define D_TERM6(a,b,c,d,e,f) a
#  define D_INVTERM6(a,b,c,d,e,f) a
#  define D_SELECT6(a,b,c,d,e,f) a
#  define LAST_SPACE_H_SPACEDIM 1

#elif CH_SPACEDIM==2
#  define D_EXPR(a,b,c) ((void)((a),(b),0))
#  define D_DECL(a,b,c) a,b
#  define D_TERM(a,b,c) a b
#  define D_INVTERM(a,b,c) b a
#  define D_SELECT(a,b,c) b
#  define D_EXPR6(a,b,c,d,e,f) ((void)((a),(b),0))
#  define D_DECL6(a,b,c,d,e,f) a,b
#  define D_TERM6(a,b,c,d,e,f) a b
#  define D_INVTERM6(a,b,c,d,e,f) b a
#  define D_SELECT6(a,b,c,d,e,f) b
#  define LAST_SPACE_H_SPACEDIM 2

#elif CH_SPACEDIM==3
#ifdef DOXYGEN
/**
@ingroup macros

  The D_EXPR(a,b,c) macro expands to a comma expression that will
  evaluate SpaceDim of its arguments, counting from the left; i.e.
  CH_SPACEDIM==1 implies only `a' is evaluated, CH_SPACEDIM==2
  implies both `a' and `b', and CH_SPACEDIM==3 implies all three
  arguments are evaluated.  This macro can be used to enhance
  the portability of code that works with structures whose size
  depends on CH_SPACEDIM.

  For example:

     D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s);

  More explicitly, this macro evaluates to one of

     ((void)((a),0))

     ((void)((a),(b),0))

     ((void)((a),(b),(c),0))

  corresponding to CH_SPACEDIM values of 1, 2, or 3.
*/
#endif
#  define D_EXPR(a,b,c) ((void)((a),(b),(c),0))
#  define D_EXPR6(a,b,c,d,e,f) ((void)((a),(b),(c),0))

#ifdef DOXYGEN
/**
@ingroup macros

  The D_DECL(a,b,c) macro expands to a comma-separated list of
  1, 2, or all 3 of the arguments of the call, depending on the
  value of CH_SPACEDIM.  This can be used to write portable
  function calls that depend on 1, 2, or 3 arguments, corresponding
  to the value of CH_SPACEDIM.

  For example:

     return IntVect(D_DECL(p[0] + s, p[1] + s, p[2] + s));

  More explicitly, this macro evaluates to one of

     a

     a,b

     a,b,c

  corresponding to CH_SPACEDIM values of 1, 2, or 3.
*/
#endif
#  define D_DECL(a,b,c) a,b,c
#  define D_DECL6(a,b,c,d,e,f) a,b,c

#ifdef DOXYGEN
/**
@ingroup macros
   define D_TERM(a,b,c) a b c

  More explicitly, this macro evaluates to one of

     a

     a b

     a b c

    corresponding to CH_SPACEDIM values of 1, 2, or 3.

    handy for situations like:

    Real norm = sqrt(D_TERM(u[0]*u[0],+u[1]*u[1],+u[2]*u[2]));
*/
#endif
#  define D_TERM(a,b,c) a b c
#  define D_TERM6(a,b,c,d,e,f) a b c

#ifdef DOXYGEN
/**
@ingroup macros
   define D_INVTERM(a,b,c) c b a

  More explicitly, this macro evaluates to one of

     a

     b a

     c b a

    corresponding to CH_SPACEDIM values of 1, 2, or 3.

    handy for situations like declaring nested loops where you want (a)
    to be the innermost or last for loop statement:

    {
      D_INVTERM(
        for (int i0 = 0; i0 != n0; ++i0)
          {
            data[(i2*n1 + i1)*n0 + i0] = val;
          },
        for (int i1 = 0; i1 != n1; ++i1) {,
        for (int i2 = 0; i2 != n2; ++i2) {)
    D_TERM(},},})
*/
#endif
#  define D_INVTERM(a,b,c) c b a
#  define D_INVTERM6(a,b,c,d,e,f) c b a

#ifdef DOXYGEN
/**
@ingroup macros
   define D_SELECT(a,b,c) a b c

  More explicitly, this macro evaluates to one of

     a

     b

     c

    corresponding to CH_SPACEDIM values of 1, 2, or 3.
*/
#endif
#  define D_SELECT(a,b,c) c
#  define D_SELECT6(a,b,c,d,e,f) c

#  define LAST_SPACE_H_SPACEDIM 3

#ifdef DOXYGEN
/** higher-dimensional (up to 6D) extensions of the original macros
   (note that the original macros are no longer defined once dim>3)
 */
#endif
#elif CH_SPACEDIM==4
#  define D_EXPR6(a,b,c,d,e,f) ((void)((a),(b),(c),(d),0))
#  define D_DECL6(a,b,c,d,e,f) a,b,c,d
#  define D_TERM6(a,b,c,d,e,f) a b c d
#  define D_INVTERM6(a,b,c,d,e,f) d c b a
#  define D_SELECT6(a,b,c,d,e,f) d
#  define LAST_SPACE_H_SPACEDIM 4

#elif CH_SPACEDIM==5
#  define D_EXPR6(a,b,c,d,e,f) ((void)((a),(b),(c),(d),(e),0))
#  define D_DECL6(a,b,c,d,e,f) a,b,c,d,e
#  define D_TERM6(a,b,c,d,e,f) a b c d e
#  define D_INVTERM6(a,b,c,d,e,f) e d c b a
#  define D_SELECT6(a,b,c,d,e,f) e
#  define LAST_SPACE_H_SPACEDIM 5

#elif CH_SPACEDIM==6
#  define D_EXPR6(a,b,c,d,e,f) ((void)((a),(b),(c),(d),(e),(f),0))
#  define D_DECL6(a,b,c,d,e,f) a,b,c,d,e,f
#  define D_TERM6(a,b,c,d,e,f) a b c d e f
#  define D_INVTERM6(a,b,c,d,e,f) f e d c b a
#  define D_SELECT6(a,b,c,d,e,f) f
#  define LAST_SPACE_H_SPACEDIM 6

#else /* matches elif CH_SPACEDIM==3 */
// note that D_TERM, D_DECL, D_EXPR, D_SELECT are not defined for DIM>3
#  define LAST_SPACE_H_SPACEDIM CH_SPACEDIM
#endif
#endif /* matches CH_SPACEDIM != LAST_SPACE_H_SPACEDIM */
